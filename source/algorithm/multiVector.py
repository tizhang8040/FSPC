from .algorithm import Algorithm
from .. import tools
import numpy as np

# %% Interface Quasi-Newton with Inverse Least Squares

class IQN_MVJ(Algorithm):
    def __init__(self,input,param,com):
        Algorithm.__init__(self,input,param)

        if com.rank == 1: 

            self.omega = param['omega']
            self.Vprev = list()
            self.Wprev = list()

# %% Coupling at Each Time Step

    def couplingAlgo(self,com):

        verified = False
        self.iteration = 0
        self.converg.epsilon = np.inf
        timeFrame = self.step.timeFrame()

        if com.rank == 1:

            self.V = list()
            self.W = list()

        while True:

            # Solid to fluid mechanical transfer

            self.clock['Communication'].start()
            self.transferDispSF(com)
            self.clock['Communication'].end()

            # Fluid solver call for FSI subiteration
            
            if com.rank == 0:

                self.clock['Solver Run'].start()
                verified = self.log.exec(self.solver.run,*timeFrame)
                self.clock['Solver Run'].end()

            verified = tools.scatterFS(verified,com)
            if not verified: return False

            # Fluid to solid mechanical transfer

            self.clock['Communication'].start()
            self.transferLoadFS(com)
            self.clock['Communication'].end()

            # Solid solver call for FSI subiteration
            
            if com.rank == 1:

                self.clock['Solver Run'].start()
                verified = self.log.exec(self.solver.run,*timeFrame)
                self.clock['Solver Run'].end()

            verified = tools.scatterSF(verified,com)
            if not verified: return False

            # Compute the mechanical residual
            
            if com.rank == 1:
            
                self.residualDispS()
                self.converg.update(self.residual)
                self.logGen.printRes()

                # Use the relaxation for solid displacement
            
                self.clock['Relax IQN-MVJ'].start()
                self.relaxation()
                self.clock['Relax IQN-MVJ'].end()
            
            # Check the converence of the FSI

            if com.rank == 1: verified = self.converg.isVerified()
            verified = tools.scatterSF(verified,com)
            self.iteration += 1

            # End of the coupling iteration

            if verified == True:
                if com.rank == 1: self.Vprev = self.V.copy()
                if com.rank == 1: self.Wprev = self.W.copy()
                return True

            if self.iteration > self.iterMax:
                if com.rank == 1: self.Vprev = list()
                if com.rank == 1: self.Wprev = list()
                return False

# %% IQN Relaxation of Solid Displacement

    def relaxation(self):

            disp = self.solver.getDisplacement()
            R = np.concatenate(self.residual.T)

            # Performs either BGS or IQN iteration

            if self.iteration > 0:

                self.V.insert(0,np.concatenate((self.residual-self.prevResidual).T))
                self.W.insert(0,np.concatenate((disp-self.prevDisp).T))
                V = np.transpose(self.V)
                W = np.transpose(self.W)

                if not self.Vprev:

                    C = np.linalg.lstsq(V,-R,rcond=-1)[0]
                    correction = np.split(np.dot(W,C)+R,self.dim)
                    self.interp.disp += np.transpose(correction)

                else:

                    Vprev = np.transpose(self.Vprev)
                    Wprev = np.transpose(self.Wprev)

                    # Make the actual multi-vector Jacobian computation

                    C1 = np.linalg.lstsq(V,R,rcond=-1)[0]
                    C2 = np.linalg.lstsq(Vprev,R,rcond=-1)[0]
                    C3 = np.linalg.lstsq(Vprev,V.dot(C1),rcond=-1)[0]
                    correction = np.split(Wprev.dot(C3-C2)-W.dot(C1)+R,self.dim)
                    self.interp.disp += np.transpose(correction)

            elif self.Vprev:

                    V = np.transpose(self.Vprev)
                    W = np.transpose(self.Wprev)
                    C = np.linalg.lstsq(V,-R,rcond=-1)[0]
                    correction = np.split(np.dot(W,C)+R,self.dim)
                    self.interp.disp += np.transpose(correction)

            else: self.interp.disp += self.omega*self.residual

            # Updates the residuals and displacement

            self.prevDisp = disp.copy()
            self.prevResidual = self.residual.copy()