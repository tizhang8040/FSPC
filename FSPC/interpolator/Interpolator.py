from mpi4py.MPI import COMM_WORLD as CW
from ..general import Toolbox as tb
from scipy import sparse as sp
import numpy as np

# %% Parent Interpolator Class

class Interpolator(object):
    def __init__(self):

        # Share the position vectors between solvers

        if CW.rank == 0:

            self.recvPos = CW.recv(source=1,tag=1)
            CW.send(tb.solver.getPosition(),1,tag=2)
            self.recvNode = self.recvPos.shape[0]

        if CW.rank == 1:

            CW.send(tb.solver.getPosition(),0,tag=1)
            self.recvPos = CW.recv(source=0,tag=2)
            self.recvNode = self.recvPos.shape[0]

        # Initialize the interpolation matrix

        self.nbrNode = tb.solver.nbrNode
        self.H = sp.dok_matrix((self.nbrNode,self.recvNode))

# %% Facets from the Target Interface Mesh

    def getFaceList(self):

        if CW.rank == 0:
            
            CW.send(tb.solver.getFace(),1,tag=7)
            self.recvFace = CW.recv(source=1,tag=8)

        if CW.rank == 1:

            self.recvFace = CW.recv(source=0,tag=7)
            CW.send(tb.solver.getFace(),0,tag=8)

# %% Interpolate RecvData and Return the Result

    @tb.compute_time
    def interpData(self,recvData):
        return self.H.dot(recvData)
    
    @tb.only_solid
    def initialize(self):

        if tb.convMech: self.pos = tb.solver.getPosition()
        if tb.convTher: self.temp = tb.solver.getTemperature()

# %% Apply Actual Loading on Solid

    @tb.conv_mecha
    def applyLoadFS(self):

        if CW.rank == 0: CW.send(tb.solver.getLoading(),1,tag=3)
        if CW.rank == 1:

            load = CW.recv(source=0,tag=3)
            tb.solver.applyLoading(self.interpData(load))

# Apply predicted displacement on fluid

    @tb.conv_mecha
    def applyDispSF(self):

        if CW.rank == 1: CW.send(self.pos,0,tag=4)
        if CW.rank == 0:

            pos = CW.recv(source=1,tag=4)
            tb.solver.applyPosition(self.interpData(pos))

# Apply actual heat flux on solid

    @tb.conv_therm
    def applyHeatFS(self):

        if CW.rank == 0: CW.send(tb.solver.getHeatFlux(),1,tag=5)
        if CW.rank == 1:
            
            heat = CW.recv(source=0,tag=5)
            tb.solver.applyHeatFlux(self.interpData(heat))

# Apply predicted temperature on fluid

    @tb.conv_therm
    def applyTempSF(self):

        if CW.rank == 1: CW.send(self.temp,0,tag=6)
        if CW.rank == 0:
            
            temp = CW.recv(source=1,tag=6)
            tb.solver.applyTemperature(self.interpData(temp))

# %% Predict the Solution for Next Time Step

    @tb.conv_mecha
    def predPosition(self,verified):
        
        if verified:

            self.prevPos = np.copy(self.pos)
            self.velocityP = tb.solver.getVelocity()

        else: self.pos = np.copy(self.prevPos)
        self.pos += tb.step.dt*self.velocityP

    # Predictor for the temparature coupling

    @tb.conv_therm
    def predTemperature(self,verified):

        if verified:

            self.prevTemp = np.copy(self.temp)
            self.velocityT = tb.solver.getTempVeloc()

        else: self.temp = np.copy(self.prevTemp)
        self.temp += tb.step.dt*self.velocityT
            
