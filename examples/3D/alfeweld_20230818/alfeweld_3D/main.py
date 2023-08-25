import os.path as path
import numpy as np
import FSPC

# %% Input Parameters for FSPC

pathF = path.dirname(__file__)+'/inputF.lua'
pathS = path.dirname(__file__)+'/inputS.py'
RBF = lambda r: np.square(r)*np.ma.log(r)

# %% Initialize the Simulation

FSPC.setConvTher(1e-6)
FSPC.setStep(0.02,0.1)
FSPC.setSolver(pathF,pathS)
FSPC.setInterp(FSPC.interpolator.RBF,RBF)
#FSPC.setInterp(FSPC.interpolator.ETM,9)

# Configure the algorithm

algorithm = FSPC.algorithm.ILS()
algorithm.maxIter = 25
algorithm.endTime = 0.1
algorithm.omega = 0.5

# Start the FSPC simulation

algorithm.simulate()
FSPC.general.printClock()