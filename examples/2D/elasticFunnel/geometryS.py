import os,gmsh
from gmsh import model as sh
gmsh.initialize()

# %% Parameters

R = 2.25
H = 3.75
s = 0.2

d = 0.05
N = 160
M = 80
P = 5

# %% Points List

p = list()

p.append(sh.occ.addPoint(-R,H,0))
p.append(sh.occ.addPoint(-(R+s),H,0))
p.append(sh.occ.addPoint(-R,0,0))
p.append(sh.occ.addPoint(-(R+s),0,0))
p.append(sh.occ.addPoint(R,H,0))
p.append(sh.occ.addPoint(R+s,H,0))
p.append(sh.occ.addPoint(R,0,0))
p.append(sh.occ.addPoint(R+s,0,0))
p.append(sh.occ.addPoint(0,0,0))

# %% Lines List

l = list()

l.append(sh.occ.addLine(p[0],p[1]))
l.append(sh.occ.addLine(p[1],p[3]))
l.append(sh.occ.addLine(p[2],p[3]))
l.append(sh.occ.addLine(p[2],p[0]))

l.append(sh.occ.addLine(p[6],p[7]))
l.append(sh.occ.addLine(p[7],p[5]))
l.append(sh.occ.addLine(p[5],p[4]))
l.append(sh.occ.addLine(p[4],p[6]))

l.append(sh.occ.addCircleArc(p[7],p[8],p[3]))
l.append(sh.occ.addCircleArc(p[6],p[8],p[2]))

# %% Solid Surface

k = list()
s = list()

k.append(sh.occ.addCurveLoop(l[0:4]))
k.append(sh.occ.addCurveLoop(l[4:8]))
k.append(sh.occ.addCurveLoop([l[2],l[8],l[4],l[9]]))
sh.occ.synchronize()

sh.mesh.setTransfiniteCurve(l[8],N)
sh.mesh.setTransfiniteCurve(l[9],N)

sh.mesh.setTransfiniteCurve(l[1],M)
sh.mesh.setTransfiniteCurve(l[3],M)
sh.mesh.setTransfiniteCurve(l[7],M)
sh.mesh.setTransfiniteCurve(l[5],M)

sh.mesh.setTransfiniteCurve(l[0],P)
sh.mesh.setTransfiniteCurve(l[2],P)
sh.mesh.setTransfiniteCurve(l[4],P)
sh.mesh.setTransfiniteCurve(l[6],P)

for a in k: s.append(sh.occ.addPlaneSurface([a]))
sh.occ.synchronize()

for a in s: sh.mesh.setTransfiniteSurface(a)
for a in s: sh.mesh.setRecombine(2,a)

# %% Physical Boundary

sh.addPhysicalGroup(2,s,name='Solid')
sh.addPhysicalGroup(1,[l[3],l[7],l[9]],name='FSInterface')
sh.addPhysicalGroup(1,[l[0],l[6]],name='SolidBase')

# %% Save the Mesh

sh.mesh.generate(2)
gmsh.option.setNumber('Mesh.Binary',1)
gmsh.write(os.path.dirname(__file__)+'/geometryS.msh')
gmsh.fltk.run()
gmsh.finalize()