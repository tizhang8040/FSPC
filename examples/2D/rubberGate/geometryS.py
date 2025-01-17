import os,gmsh
from gmsh import model as sh
gmsh.initialize()

# %% Parameters

H = 0.079
S = 0.005
L1 = 0.1

eps = 1e-5
N = 60
M = 5

# %% Point List

p = list()

p.append(sh.occ.addPoint(L1,H,0))
p.append(sh.occ.addPoint(L1+S,H,0))
p.append(sh.occ.addPoint(L1,eps,0))
p.append(sh.occ.addPoint(L1+S,eps,0))

# %% Line List

l = list()

l.append(sh.occ.addLine(p[3],p[1]))
l.append(sh.occ.addLine(p[1],p[0]))
l.append(sh.occ.addLine(p[0],p[2]))
l.append(sh.occ.addLine(p[2],p[3]))

# %% Solid Mesh

k = sh.occ.addCurveLoop(l)
s = sh.occ.addPlaneSurface([k])
sh.occ.synchronize()

sh.mesh.setTransfiniteCurve(l[0],N)
sh.mesh.setTransfiniteCurve(l[2],N)
sh.mesh.setTransfiniteCurve(l[1],M)
sh.mesh.setTransfiniteCurve(l[3],M)

sh.mesh.setTransfiniteSurface(s)
sh.mesh.setRecombine(2,s)

# %% Physycal Boundary

sh.addPhysicalGroup(2,[s],name='Solid')
sh.addPhysicalGroup(1,[l[2]],name='FSInterface')
sh.addPhysicalGroup(1,[l[1]],name='SolidBase')

# %% Save the Mesh

sh.mesh.generate(2)
gmsh.option.setNumber('Mesh.Binary',1)
gmsh.write(os.path.dirname(__file__)+'/geometryS.msh')
gmsh.fltk.run()
gmsh.finalize()