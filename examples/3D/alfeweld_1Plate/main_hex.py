import os.path as path
import FSPC
#import numpy as np

# %% Input Parameters for FSPC

pathF = path.dirname(__file__)+'/inputF.lua'
pathS = path.dirname(__file__)+'/inputS_hex.py'
#RBF = lambda r: np.square(r)*np.ma.log(r)

# %% Initialize the Simulation

FSPC.setConvTher(1e-6)
FSPC.setStep(0.01,0.1)
FSPC.setSolver(pathF,pathS)
FSPC.setInterp(FSPC.interpolator.ETM,9)
#FSPC.setInterp(FSPC.interpolator.RBF,RBF)

# Configure the algorithm

algorithm = FSPC.algorithm.MVJ()
#algorithm = FSPC.algorithm.ILS()
algorithm.maxIter = 25
algorithm.endTime = 26
algorithm.omega = 0.5

# Start the FSPC simulation

algorithm.simulate()
FSPC.general.printClock()
