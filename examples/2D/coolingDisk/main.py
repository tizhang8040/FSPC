import os.path as path
import FSPC

# %% Input Parameters for FSPC

pathF = path.dirname(__file__)+'/inputF.lua'
pathS = path.dirname(__file__)+'/inputS.py'

# %% Initialize the Manager Module

FSPC.setConvMech(1e-8)
FSPC.setConvTher(1e-6)
FSPC.setStep(1e-2,0.01)
FSPC.setSolver(pathF,pathS)
FSPC.setInterp(FSPC.interpolator.KNN,2)

# Configure the algorithm

algorithm = FSPC.algorithm.ILS()
algorithm = FSPC.algorithm.MVJ()
algorithm.maxIter = 25
algorithm.endTime = 8
algorithm.omega = 0.5

# Start the FSPC simulation

algorithm.simulate()
FSPC.general.printClock()