from .Interpolator import Interpolator
from ..general import Toolbox as tb
from scipy import sparse as sp
import numpy as np

# %% Mesh Interpolation K-Nearest Neighbours

class KNN(Interpolator):
    def __init__(self,K):
        self.K = int(abs(K))

    # Compute the FS mesh interpolation matrix

    def initialize(self):

        Interpolator.__init__(self)
        position = tb.solver.getPosition()
        self.computeMapping(position)
        self.H = self.H.tocsr()

    # Interpolate recvData and return the result

    @tb.compute_time
    def interpData(self,recvData):
        return self.H.dot(recvData)

# %% Mapping Matrix from RecvPos to Position

    @tb.compute_time
    def computeMapping(self,position):

        size = tb.solver.nbrNod,len(self.recvPos)
        self.H = sp.dok_matrix(size)

        if self.K == 1: self.search(position)
        else: self.interpolate(position)

# %% Find the K Nearest Neighbours
 
    def search(self,position):

        for i,pos in enumerate(position):

            dist = np.linalg.norm(pos-self.recvPos,axis=1)
            self.H[i,np.argmin(dist)] = 1

    def interpolate(self,position):

        for i,pos in enumerate(position):

            dist = np.linalg.norm(pos-self.recvPos,axis=1)
            index = np.argsort(dist)[range(self.K)]
            weight = np.zeros(self.K)
            dist = dist[index]

            for j in range(self.K):

                val = [R for k,R in enumerate(dist) if k != j]
                weight[j] = np.prod(val)

            self.H[i,index] = weight/np.sum(weight)
