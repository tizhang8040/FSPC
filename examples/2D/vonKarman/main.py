import os.path as path
import numpy as np
import FSPC

# %% Paths to the input files

pathF = path.dirname(__file__)+'/inputF.lua'
pathS = path.dirname(__file__)+'/inputS.py'

# %% Fluid Structure Coupling

process = FSPC.Process()
solver = process.getSolver(pathF,pathS)
com = process.com

RBF = lambda r: np.square(r/1e-3)*np.ma.log(r/1e-3)

# Configure the algorithm

algorithm = FSPC.IQN_MVJ(solver,com)
algorithm.interp = FSPC.RBF(solver,RBF,com)
algorithm.converg = FSPC.Convergence(1e-6)
algorithm.step = FSPC.TimeStep(1e-3)

algorithm.endTime = 10
algorithm.omega = 0.5
algorithm.iterMax = 25
algorithm.dtWrite = 1e-3

# Start the FSPC simulation

algorithm.simulate(com)
FSPC.printClock(com)