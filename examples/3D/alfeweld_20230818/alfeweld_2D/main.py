import os.path as path
import numpy as np
import FSPC

# %% Paths to the input files

pathF = path.dirname(__file__)+'/inputF.lua'
pathS = path.dirname(__file__)+'/inputS.py'

# Configure the algorithm
RBF = lambda r: np.square(r)*np.ma.log(r)

FSPC.setConvTher(1e-6)
FSPC.setStep(0.02,0.1)
FSPC.setSolver(pathF,pathS)
FSPC.setInterp(FSPC.interpolator.RBF,RBF)

algorithm = FSPC.algorithm.ILS()

algorithm.endTime = 15 # 60 for the full simulation
algorithm.omega = 0.5
algorithm.maxIter = 20

# Start the FSPC simulation

algorithm.simulate()
FSPC.general.printClock()