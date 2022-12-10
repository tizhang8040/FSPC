from matplotlib import pyplot as plt
from Externals import tools
from os import listdir
import numpy as np
import os

# %% Reads the Data

folder = os.getcwd()+'/validation/MarcoLucio/'
data = [np.loadtxt(folder+file).T for file in listdir(folder)]

# %% Print the Mass

workspace = os.getcwd()+'/workspace'
os.chdir(workspace)

# Reads the results

mass = np.loadtxt('mass.txt',delimiter=',').T
mass[1] = 100*mass[1]/mass[1,0]

# Save Results

plt.figure(1)
plt.gca().set_title('Total Mass')
np.savetxt('mass.dat',np.transpose(mass),fmt=['%.6f','%.6f'])
plt.plot(*mass)
plt.grid()
plt.show()

# %% Print the Output

workspace = os.getcwd()+'/metafor'
os.chdir(workspace)

# Reads the results

tag = tools.getIndex([0.292,0.1,0.08])
time,disp = tools.readNode(tag)
results = [time,disp[:,0]]

# Moves to main folder

workspace = os.path.split(workspace)[0]
os.chdir(workspace)

# Save Results

plt.figure(2)
plt.gca().set_title('Output Data')
np.savetxt('output.dat',np.transpose(results),fmt=['%.6f','%.6f'])
for curve in data: plt.plot(*curve)
plt.plot(*results,'k--')
plt.grid()
plt.show()