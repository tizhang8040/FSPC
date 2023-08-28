import toolbox.gmsh as gmsh
import wrap as w
import numpy as np
import math
import os


# ===================
# ===== tianyu ======
# ===================

# ===================
# ===== tianyu ======
# ===================

def createRigidCone(geometry, idGeo, X0, Y0, Z0, rExt, rInt, hCyld, rf, angCone):
    #rf: radius of the rounded fillet
    #rh: radius of the rounded head
    sinA = math.sin(math.radians(45))  # = cosA
    sinB = math.sin(math.radians(angCone))
    cosB = math.cos(math.radians(angCone))
    tanB = math.tan(math.radians(angCone))
    sinC = math.sin(math.radians(angCone/2))
    cosC = math.cos(math.radians(angCone/2))
    
    rh = rInt/sinB
    
    dx7 = 0
    dz7 = 0
    
    dx6 = rh*sinC
    dz6 = rh*(1-cosC)
    
    dx5 = rInt
    dz5 = rh*(1-cosB)
    
    dx4 = rExt-rf*(1-sinB)
    dz4 = dx4*tanB - (rInt*tanB-dz5)
    
    dx2 = rExt
    dz2 = dz4 + rf*cosB
        
    dx3 = rExt - rf*(1-sinA)
    dz3 = dz2 - rf*sinA
    
    dx1 = rExt
    dz1 = dz2 + hCyld
    # sets
    pointset = geometry.getPointSet()
    curveset = geometry.getCurveSet()
    wireset = geometry.getWireSet()
    surfaceset = geometry.getSurfaceSet()
    # -- Point --
    #generatrix
    pointset.define(idGeo+1, X0-dx1, Y0, Z0+dz1)
    pointset.define(idGeo+2, X0-dx2, Y0, Z0+dz2)
    pointset.define(idGeo+3, X0-dx3, Y0, Z0+dz3)
    pointset.define(idGeo+4, X0-dx4, Y0, Z0+dz4)
    pointset.define(idGeo+5, X0-dx5, Y0, Z0+dz5)
    pointset.define(idGeo+6, X0-dx6, Y0, Z0+dz6)
    pointset.define(idGeo+7, X0-dx7, Y0, Z0+dz7)
    #axis
    pointset.define(idGeo+51, X0, Y0, Z0)
    pointset.define(idGeo+52, X0, Y0, Z0+dz1)

    # -- Line - Wire --
    #curveset.add(Line(idGeo+51, pointset(idGeo+51), pointset(idGeo+52)))
    curveset.add(w.Line(idGeo+51, pointset(idGeo+52), pointset(idGeo+51)))
    curveset.add(w.Line(idGeo+1, pointset(idGeo+1), pointset(idGeo+2)))
    curveset.add(w.Arc (idGeo+2, pointset(idGeo+2), pointset(idGeo+3), pointset(idGeo+4)))
    curveset.add(w.Line(idGeo+3, pointset(idGeo+4), pointset(idGeo+5)))
    curveset.add(w.Arc (idGeo+4, pointset(idGeo+5), pointset(idGeo+6), pointset(idGeo+7)))
    wireset.add(w.Wire(idGeo+1, [curveset(idGeo+1),curveset(idGeo+2),curveset(idGeo+3),curveset(idGeo+4)]))
    # -- Surface --
    surfRC = surfaceset.add(w.RevolutionSurface(idGeo+1, curveset(idGeo+51), wireset(idGeo+1)))
    
    return surfRC


def defineTmFrictionlessContactMaterial(materset, materIndex, normalPeno, prof, hNom, hardVickers, wExp):
    materset.define(materIndex, w.TmFrictionlessContactMaterial)
    materset(materIndex).put(w.PEN_NORMALE,    normalPeno)
    materset(materIndex).put(w.PROF_CONT,      prof)
    materset(materIndex).put(w.CTM_H_NOMINAL,  hNom)
    materset(materIndex).put(w.CTM_HARDNESS,   hardVickers)
    materset(materIndex).put(w.CTM_EXPONENT_E, wExp)


def defineTmRDContactElePropInter(eleType, materIndex, stiffMethod, AICMethod, T_Tool, fctT_Tool, interset, interIndex, gTool, gSlave, setSN):
    prpElem = w.ElementProperties(eleType)
    prpElem.put(w.MATERIAL,      materIndex)
    prpElem.put(w.STIFFMETHOD,   stiffMethod)
    prpElem.put(w.AREAINCONTACT, AICMethod)
    prpElem.put(w.TOOLTEMP,      T_Tool)
    prpElem.depend(w.TOOLTEMP, fctT_Tool, w.Field1D(w.TM,w.RE))
    
    rdci = w.RdContactInteraction(interIndex)
    rdci.setTool(gTool)
    rdci.push(gSlave)
    rdci.setSmoothNormals(setSN)
    rdci.addProperty(prpElem)
    interset.add(rdci)
    
    return rdci

# ===================
# ===================
# ===================

# %% Main Function

metafor = None
metaforStandalone = False

def getMetafor(parm):

    global metafor
    if metafor: return metafor

    w.StrVectorBase.useTBB()
    w.StrMatrixBase.useTBB()
    w.ContactInteraction.useTBB()
    
    # ======= Tianyu =======
    # Parameters
    
    #ids ------------------------
    idGeo_RT = 1000             #   
    idMatCI = 1000              #   
    #geo ------------------------
    L = 110 		                #   solid domain length
    HF = 3.0  		            #   fluid domain height
    HS = 1.5  		            #   solid domain height
    W = 80                      #   solid domain width
    
    rExt_RT = 10.          #   rigid tool radius
    rf_RT = rExt_RT/10          #   radius of the rounded fillet
    hCyld_RT = 3.            #
    rInt_RT = rExt_RT/10        #
    angCone_RT = 0.1            #   degree
    
    X0_RT = 0.- L/2 + 20.          #   origin x of the rigid tool center
    Y0_RT = 0.         #   origin y of the rigid tool center with a small inter-penetration
    Z0_RT = HF+HS-1.e-3                 #   origin z of the rigid tool center
    #damped ---------------------
    #dampStiff = 1.e-5
    #dampMass = 0.0
    #contact --------------------
    peno_CI = 1.5e3             #
    facProf_CI = 0.95           #
    hc0_CI = 250.               #
    HV_CI = 500.                #
    w_CI = 1.25                 #
    eleType_CI = w.TmContact3DElement   #   
    T_Tool = 1150.                      #   Temperature of the tool 
    setSN_CI = True #False              #   
    #time step manager ------------------   
    time0 = 0.0                         #   initial time
    time1 = time0+5.0                   #   vertical displacement time (tool)
    time2 = time1+15.0                  #   horizontal displacement time (tool)
    dt0 = 5.e-5
    dtMax1 = 5.e-1
    dtMax2 = 5.e-1
    nbArchi1 = 10
    nbArchi2 = 50
    # ======================

    # Group and interaction sets

    metafor = w.Metafor()
    domain = metafor.getDomain()
    # ======= tianyu =======
    geometry = domain.getGeometry()
    # ======================
    tsm = metafor.getTimeStepManager()
    materset = domain.getMaterialSet()
    # ====== tianyu ========
    loadingset = domain.getLoadingSet()
    # ======================
    tim = metafor.getThermalIterationManager()
    solvermanager = metafor.getSolverManager()
    interactionset = domain.getInteractionSet()
    mim = metafor.getMechanicalIterationManager()
    initcondset = metafor.getInitialConditionSet()

    # Dimension and DSS solver

    domain.getGeometry().setDim3D()
    solvermanager.setSolver(w.DSSolver())
    
    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__),'geometryS_hex.msh')
    importer = gmsh.GmshImport(mshFile,domain)
    groups = importer.groups
    importer.execute()

    # ===== tianyu =======
    # Defines the rigid tool
    surf_RT = createRigidCone(geometry, idGeo_RT, X0_RT, Y0_RT, Z0_RT, rExt_RT, rInt_RT, hCyld_RT, rf_RT, angCone_RT)
    # ====================

    # Defines the solid domain

    app = w.FieldApplicator(1)
    app.push(groups['FeTop'])
    interactionset.add(app)

    # Solid material parameters

    materset.define(1,w.TmElastHypoMaterial)
    materset(1).put(w.ELASTIC_MODULUS,210.e3)
    materset(1).put(w.THERM_EXPANSION,0)
    materset(1).put(w.HEAT_CAPACITY,500.0e6)
    materset(1).put(w.MASS_DENSITY,7.85e-9)
    materset(1).put(w.CONDUCTIVITY,30.0)
    materset(1).put(w.POISSON_RATIO,0.28)
    materset(1).put(w.DISSIP_TE,0)
    materset(1).put(w.DISSIP_TQ,0)

    # ===== tianyu =======
    # Contact material parameters
    defineTmFrictionlessContactMaterial(materset, idMatCI, peno_CI, rf_RT*facProf_CI, hc0_CI, HV_CI, w_CI)
    # ====================

    # Finite element properties

    prp1 = w.ElementProperties(w.TmVolume3DElement)
    prp1.put(w.CAUCHYMECHVOLINTMETH,w.VES_CMVIM_SRIPR)
    prp1.put(w.STIFFMETHOD,w.STIFF_ANALYTIC)
    # ===== tianyu =======
    #prp1.put(w.DAMPSTIFF, dampStiff)
    #prp1.put(w.DAMPMASS, dampMass)
    # ====================
    prp1.put(w.MATERIAL,1)
    app.addProperty(prp1)

    # ===== tianyu =======
    # Contact finite element properties
    
    fctT_Tool = w.PieceWiseLinearFunction()
    fctT_Tool.setData(time0, 1.0)
    fctT_Tool.setData(time1, 1.0)
    fctT_Tool.setData(time2, 1.0)
    ci = defineTmRDContactElePropInter(eleType_CI, idMatCI, w.STIFF_ANALYTIC, w.AIC_ONCE, T_Tool, fctT_Tool, interactionset, idMatCI, surf_RT, groups['FeTop_Up'], setSN_CI)
    # ====================

    # Elements for surface heat flux

    prp2 = w.ElementProperties(w.NodHeatFlux3DElement)
    heat = w.NodInteraction(2)
    heat.push(groups['FSInterface'])
    heat.addProperty(prp2)
    interactionset.add(heat)

    # ===== tianyu =======
    # Top sliding RT boundary conditions

    AmplDispX_RT = 2*abs(X0_RT)  # horizontal displacement
    fct_dispX = w.PieceWiseLinearFunction()
    fct_dispX.setData(time0, 0.0)
    fct_dispX.setData(time1, 0.0)
    fct_dispX.setData(time2, 1.0)

    AmplDispZ_RT = -0.03            # vertical displacement
    fct_dispZ = w.PieceWiseLinearFunction()
    fct_dispZ.setData(time0, 0.0)
    fct_dispZ.setData(time1, 1.0)
    fct_dispZ.setData(time2, 1.0)

    loadingset.define(surf_RT, w.Field1D(w.TX,w.RE), AmplDispX_RT, fct_dispX)
    loadingset.define(surf_RT, w.Field1D(w.TY,w.RE))
    loadingset.define(surf_RT, w.Field1D(w.TZ,w.RE), AmplDispZ_RT, fct_dispZ)
    # ====================

    # Initial and boundary conditions

    
    initcondset.define(groups['FeTop'],w.Field1D(w.TO,w.AB),300)
    #initcondset.define(groups['Bottom'],w.Field1D(w.TO,w.AB),300.0)
    #initcondset.define(groups['FSInterface'],w.Field1D(w.TO,w.AB),300)

    #loadingset.define(groups['FeTop_Down'],w.Field1D(w.TO,w.RE))
    loadingset.define(groups['FeTop_Down'],w.Field1D(w.TX,w.RE))
    loadingset.define(groups['FeTop_Down'],w.Field1D(w.TY,w.RE))
    loadingset.define(groups['FeTop_Down'],w.Field1D(w.TZ,w.RE))
    
    loadingset.define(groups['FSInterface'],w.Field1D(w.TX,w.RE))
    loadingset.define(groups['FSInterface'],w.Field1D(w.TY,w.RE))
    loadingset.define(groups['FSInterface'],w.Field1D(w.TZ,w.RE))
    
    #if metaforStandalone:
    #    loadingset.define(groups['FSInterface'],w.Field1D(w.TO,w.RE))


    # Mechanical and thermal time integration

    ti_M = w.AlphaGeneralizedTimeIntegration(metafor)
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

    # Parameters for FSPC
    if metaforStandalone:
        # time step manager
        tsm.setInitialTime(time0, dt0)
        tsm.setNextTime(time1, nbArchi1, dtMax1)
        tsm.setNextTime(time2, nbArchi2, dtMax2)
        
        # archiving
        valuesmanager = metafor.getValuesManager()
        #1-time step
        valuesmanager.add(1000, w.MiscValueExtractor(metafor,w.EXT_T), 'time')
        valuesmanager.add(1001, w.MiscValueExtractor(metafor,w.EXT_DT), 'dt')
        
        dataCurve1 = w.VectorDataCurve(1001, valuesmanager.getDataVector(1000), valuesmanager.getDataVector(1001),'dt')
        dataCurveSet1 = w.DataCurveSet()
        dataCurveSet1.add(dataCurve1)
        winc1 = w.VizWin()
        winc1.add(dataCurveSet1)
        metafor.addObserver(winc1)

    if not metaforStandalone:
        parm['interacT'] = heat
        parm['FSInterface'] = groups['FSInterface']
        parm['exporter'] = gmsh.GmshExport('metafor/output.msh',metafor)
        parm['exporter'].addDataBaseField([w.TO])
    
    return metafor
