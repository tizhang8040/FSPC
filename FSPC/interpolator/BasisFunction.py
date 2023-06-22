from .Interpolator import Interpolator
from ..general import Toolbox as tb
import numpy as np

# %% Mesh Interpolation with Radial Basis Functions

class RBF(Interpolator):
    def __init__(self,func):
        Interpolator.__init__(self)

        # Compute the FS mesh interpolation matrix

        self.function = func
        position = tb.solver.getPosition()
        self.computeMapping(position)

# %% Mapping Matrix from RecvPos to Position

    @tb.compute_time
    def computeMapping(self,pos):

        size = self.recvNode+tb.solver.dim+1
        B = np.ones((self.nbrNode,size))
        A = np.zeros((size,size))

        # Fill the matrices A,B with nodal positions

        B[:,self.recvNode+1:] = pos
        A[:self.recvNode,self.recvNode+1:] = self.recvPos
        A[:self.recvNode,self.recvNode] = 1
        A += np.transpose(A)

        # Fill the matrices A,B using the basis function

        for i,position in enumerate(self.recvPos):
            
            rad = np.linalg.norm(position-self.recvPos,axis=1)
            A[i,:self.recvNode] = self.function(rad)

            rad = np.linalg.norm(pos-position,axis=1)
            B[:,i] = self.function(rad)

        # Compute the interpolation H matrix

        try: H = np.linalg.lstsq(A.T,B.T,-1)[0].T
        except: H = np.linalg.solve(A.T,B.T).T
        self.H = H[:self.nbrNode,:self.recvNode]
