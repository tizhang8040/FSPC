import os,gmsh
from gmsh import model as sh
gmsh.initialize()

# %% Parameters

L = 1
HF = 0.2
HS = 0.02
d = 0.01

# %% Points List

p = list()

p.append(sh.occ.addPoint(L,HS,0,d))
p.append(sh.occ.addPoint(0,HS,0,d))
p.append(sh.occ.addPoint(L,HS+HF,0,d))
p.append(sh.occ.addPoint(0,HS+HF,0,d))

# %% Lines List

l = list()

l.append(sh.occ.addLine(p[0],p[1]))
l.append(sh.occ.addLine(p[0],p[2]))
l.append(sh.occ.addLine(p[2],p[3]))
l.append(sh.occ.addLine(p[3],p[1]))

# %% Fluid Surface

k = sh.occ.addCurveLoop(l)
s = sh.occ.addPlaneSurface([k])
sh.occ.synchronize()

# %% Physical Boundary

sh.addPhysicalGroup(2,[s],name='Fluid')
sh.addPhysicalGroup(1,[l[0]],name='FSInterface')
sh.addPhysicalGroup(1,[l[2]],name='FreeSurface')
sh.addPhysicalGroup(1,[l[1],l[3]],name='Wall')

# %% Save the Mesh

sh.mesh.generate(2)
gmsh.option.setNumber('Mesh.Binary',1)
gmsh.write(os.path.dirname(__file__)+'/geometryF.msh')
gmsh.fltk.run()
gmsh.finalize()