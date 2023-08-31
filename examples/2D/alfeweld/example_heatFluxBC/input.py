import toolbox.gmsh as gmsh
import numpy as np
import wrap as w
import os

# %% Main Function

metafor = None
def getMetafor(input):

    global metafor
    if metafor: return metafor

    w.StrVectorBase.useTBB()
    w.StrMatrixBase.useTBB()
    w.ContactInteraction.useTBB()
    
    # Group and interaction sets

    metafor = w.Metafor()
    domain = metafor.getDomain()
    tsm = metafor.getTimeStepManager()
    materset = domain.getMaterialSet()
    tim = metafor.getThermalIterationManager()
    solvermanager = metafor.getSolverManager()
    interactionset = domain.getInteractionSet()
    mim = metafor.getMechanicalIterationManager()
    initcondset = metafor.getInitialConditionSet()

    # Dimension and DSS solver

    domain.getGeometry().setDimPlaneStrain(1)
    solvermanager.setSolver(w.DSSolver())
    
    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__),'geometry.msh')
    importer = gmsh.GmshImport(mshFile,domain)
    groups = importer.groups
    importer.binary = True
    importer.execute()

    # Defines the ball domain

    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
    interactionset.add(app)

    # Solid material parameters

    materset.define(1,w.TmElastHypoMaterial)
    materset(1).put(w.ELASTIC_MODULUS,1e7)
    materset(1).put(w.THERM_EXPANSION,0)
    materset(1).put(w.HEAT_CAPACITY,100)
    materset(1).put(w.MASS_DENSITY,950)
    materset(1).put(w.POISSON_RATIO,0)
    materset(1).put(w.CONDUCTIVITY,5)
    materset(1).put(w.DISSIP_TE,0)
    materset(1).put(w.DISSIP_TQ,0)

    # Finite element properties

    prp1 = w.ElementProperties(w.TmTriangleVolume2DElement)
    prp1.put(w.CAUCHYMECHVOLINTMETH,w.VES_CMVIM_SRIPR)
    prp1.put(w.STIFFMETHOD,w.STIFF_ANALYTIC)
    prp1.put(w.MATERIAL,1)
    app.addProperty(prp1)

    # Elements for surface heat flux

    # prp2 = w.ElementProperties(w.NodHeatFlux2DElement)
    # heat = w.NodInteraction(2)
    # heat.push(groups['FSInterface'])
    # heat.addProperty(prp2)
    # interactionset.add(heat)

    # for i in range(groups['FSInterface'].getNumberOfMeshPoints()):

    #     node = groups['FSInterface'].getMeshPoint(i)
    #     x = node.getValue(w.Field1D(w.TX, w.AB))
    #     y = node.getValue(w.Field1D(w.TY, w.AB))
    #     norm = np.array([0.5,0.5])-np.array([x,y])
    #     flux =  2e4*norm/np.linalg.norm(norm)
    #     heat.setNodVector(node,*flux)

    # Elements for classical heat flux

    fun = lambda x : np.exp(-np.square(x-0.5)/np.square(0.1)/2)
    F = w.PythonOneParameterFunction(fun)

    prp2 = w.ElementProperties(w.TmHeatFlux2DElement)
    prp2.put(w.HEATEL_VALUE,1e3)
    prp2.depend(w.HEATEL_VALUE,F,w.Field1D(w.TX,w.AB))

    heat = w.HeatInteraction(3)
    heat.push(groups['FSInterface'])
    heat.addProperty(prp2)
    interactionset.add(heat)

    # Boundary conditions

    initcondset.define(groups['Solid'],w.Field1D(w.TO,w.AB),180)
    initcondset.define(groups['FSInterface'],w.Field1D(w.TO,w.AB),180)

    # Mechanical and thermal time integration

    ti_M = w.QuasiStaticTimeIntegration(metafor)
    ti_T = w.TrapezoidalThermalTimeIntegration(metafor)

    ti = w.StaggeredTmTimeIntegration(metafor)
    ti.setMechanicalTimeIntegration(ti_M)
    ti.setThermalTimeIntegration(ti_T) 
    metafor.setTimeIntegration(ti)

    # Mechanical and thermal iterations

    mim.setMaxNbOfIterations(25)
    mim.setResidualTolerance(1e-6)

    tim.setMaxNbOfIterations(25)
    tim.setResidualTolerance(1e-6)

    # Time step iterations
    
    tscm = w.NbOfStaggeredTmNRIterationsTimeStepComputationMethod(metafor)
    tsm.setTimeStepComputationMethod(tscm)
    tscm.setTimeStepDivisionFactor(2)
    tscm.setNbOptiIte(25)

    # Time step manager

    Tend = 50
    dTmax = 1e-2

    tsm.setInitialTime(0,dTmax)
    tsm.setNextTime(Tend,100,dTmax)
    return metafor