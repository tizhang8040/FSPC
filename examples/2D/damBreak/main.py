import os.path as path
import numpy as np
import FSPC

# %% Input Parameters for FSPC

pathF = path.dirname(__file__)+'/inputF.lua'
pathS = path.dirname(__file__)+'/inputS.py'
RBF = lambda r: np.square(r)*np.ma.log(r)

# %% Initialize the Manager Module

FSPC.setConvMech(1e-8)
FSPC.setStep(1e-3,1e-2)
FSPC.setSolver(pathF,pathS)
FSPC.setInterp(FSPC.interpolator.RBF,RBF)

# Configure the algorithm

algorithm = FSPC.algorithm.MVJ()
algorithm.maxIter = 25
algorithm.endTime = 1
algorithm.omega = 0.5

# Start the FSPC simulation

algorithm.simulate()
FSPC.general.printClock()