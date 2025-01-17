import os.path as path
import FSPC

# %% Paths to the input files

pathF = path.dirname(__file__)+'/inputF.lua'
pathS = path.dirname(__file__)+'/inputS.py'

# %% Fluid Structure Coupling

process = FSPC.Process()
solver = process.getSolver(pathF,pathS)

# Configure the algorithm

algorithm = FSPC.MVJ(solver)
algorithm.interp = FSPC.KNN(solver,2)
algorithm.convergT = FSPC.Convergence(1e-6)
algorithm.step = FSPC.TimeStep(0.02,0.1)

algorithm.endTime = 60
algorithm.omega = 0.5
algorithm.maxIter = 20

# Start the FSPC simulation

algorithm.simulate()
FSPC.printClock()