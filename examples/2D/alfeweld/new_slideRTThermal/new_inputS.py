import toolbox.gmsh as gmsh
import numpy as np
import wrap as w
import os,math



def create2DRigidCylinder(geometry, idGeo, X0, Y0, Z0, rExt, h, rf):
    #rf: radius of the rounded fillet
    piDemi = math.asin(1)
    sinA = math.sin(piDemi/2.0)
    # -- Point --
    #left
    pointset = geometry.getPointSet()
    pointset.define(idGeo+1, X0-rExt,             Y0+h+rf,        Z0)
    pointset.define(idGeo+2, X0-rExt,             Y0+rf,          Z0)
    pointset.define(idGeo+3, X0-rExt+rf*(1-sinA), Y0+rf*(1-sinA), Z0)
    pointset.define(idGeo+4, X0-rExt+rf,          Y0,             Z0)
    #right
    pointset.define(idGeo+5, X0+rExt-rf,          Y0,             Z0)
    pointset.define(idGeo+6, X0+rExt-rf*(1-sinA), Y0+rf*(1-sinA), Z0)
    pointset.define(idGeo+7, X0+rExt,             Y0+rf,          Z0)
    pointset.define(idGeo+8, X0+rExt,             Y0+h+rf,        Z0)

    # -- Line - Wire --
    curveset = geometry.getCurveSet()
    curveset.add(w.Line(idGeo+1, pointset(idGeo+1), pointset(idGeo+2)))
    curveset.add(w.Arc (idGeo+2, pointset(idGeo+2), pointset(idGeo+3), pointset(idGeo+4)))
    curveset.add(w.Line(idGeo+3, pointset(idGeo+4), pointset(idGeo+5)))
    curveset.add(w.Arc (idGeo+4, pointset(idGeo+5), pointset(idGeo+6), pointset(idGeo+7)))
    curveset.add(w.Line(idGeo+5, pointset(idGeo+7), pointset(idGeo+8)))
    wireset = geometry.getWireSet()
    surfRC = wireset.add(w.Wire(idGeo+1, [curveset(idGeo+1),curveset(idGeo+2),curveset(idGeo+3),curveset(idGeo+4),curveset(idGeo+5)]))
    
    return surfRC



def defineTmCoulombContactMaterial(materset, materIndex, normalPeno, tangPeno, prof, statCoefFrot, dynCoefFrot, hNom, hardVickers, wExp):
    materset.define(materIndex, w.TmCoulombContactMaterial)
    materset(materIndex).put(w.PEN_NORMALE,    normalPeno)
    materset(materIndex).put(w.PEN_TANGENT,    tangPeno)
    materset(materIndex).put(w.PROF_CONT,      prof)
    materset(materIndex).put(w.COEF_FROT_STA,  statCoefFrot)
    materset(materIndex).put(w.COEF_FROT_DYN,  dynCoefFrot)
    materset(materIndex).put(w.CTM_H_NOMINAL,  hNom)
    materset(materIndex).put(w.CTM_HARDNESS,   hardVickers)
    materset(materIndex).put(w.CTM_EXPONENT_E, wExp)



def defineTmRDContElePropInter(eleType, materIndex, stiffMethod, AICMethod, T_Tool, fctT_Tool, interset, interIndex, gTool, gSlave, setSN):
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



# %% Main Function

metaforStandalone = True

metafor = None
def getMetafor(input):

    global metafor
    if metafor: return metafor

    w.StrVectorBase.useTBB()
    w.StrMatrixBase.useTBB()
    w.ContactInteraction.useTBB()
    
    # Parameters
    
    #ids
    idGeo_RT = 1000
    idMatCI = 1000
    #geo
    geoRT = 'cylinder'  #cylinder/cone
    L = 1  #solid domain length
    HF = 0.2  #fluid domain height
    HS = 0.1  #solid domain height
    N = 101  #number of nodes in x direction
    clx = L/(N-1)  #mesh size in x direction
    rExt_RT = 0.05  #rigid tool radius
    h_RT = 0.025  #rigid tool height
    facRf_RT = 1.5  #factor of the rounded fillet radius
    rf_RT = facRf_RT*clx  #radius of the rounded fillet
    X0_RT = rExt_RT+0.05  #origin x of the rigid tool center
    Y0_RT = HF+HS-1.e-6  #origin y of the rigid tool center with a small inter-penetration
    Z0_RT = 0  #origin z of the rigid tool center
    #contact    
    mu_CI = 0.
    peno_CI = 1.e4
    peta_CI = peno_CI*mu_CI
    facProf_CI = 0.95
    hc0_CI = 250.
    HV_CI = 500.
    w_CI = 1.25
    eleType_CI = w.TmContact2DElement
    T_Tool = 1200.
    setSN_CI = False
    #time step manager
    time0 = 0.0
    time1 = time0+3.0
    time2 = time1+17.0
    dt0 = 5.e-5
    dtMax1 = 5.e-1
    dtMax2 = 5.e-1
    nbArchi1 = 5
    nbArchi2 = 30
    
    # Group and interaction sets

    metafor = w.Metafor()
    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    tsm = metafor.getTimeStepManager()
    materset = domain.getMaterialSet()
    loadingset = domain.getLoadingSet()
    tim = metafor.getThermalIterationManager()
    solvermanager = metafor.getSolverManager()
    interactionset = domain.getInteractionSet()
    mim = metafor.getMechanicalIterationManager()
    initcondset = metafor.getInitialConditionSet()

    # Dimension and DSS solver

    domain.getGeometry().setDimPlaneStrain(1)
    solvermanager.setSolver(w.DSSolver())
    
    # Imports the mesh

    mshFile = os.path.join(os.path.dirname(__file__),'geometryS.msh')
    importer = gmsh.GmshImport(mshFile,domain)
    groups = importer.groups
    importer.binary = True
    importer.execute()
    
    # Defines the solid domain

    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
    interactionset.add(app)
    
    # Defines the rigid tool
    
    if geoRT == 'cylinder':
        surf_RT = create2DRigidCylinder(geometry, idGeo_RT, X0_RT, Y0_RT, Z0_RT, rExt_RT, h_RT, rf_RT)

    # check geometry & mesh
    if 0:
        win = VizWin()
        win.add(geometry.getPointSet())
        win.add(geometry.getCurveSet())
        #win.add(geometry.getMesh().getPointSet())
        win.add(geometry.getMesh().getCurveSet())
        #win.add(geometry.getGroupSet())
        win.open()
        input()
    
    # Solid material parameters

    materset.define(1,w.TmElastHypoMaterial)
    materset(1).put(w.ELASTIC_MODULUS,210.e3)
    materset(1).put(w.THERM_EXPANSION,0)
    materset(1).put(w.HEAT_CAPACITY,500.e6)
    materset(1).put(w.MASS_DENSITY,7.85e-9)
    materset(1).put(w.POISSON_RATIO,0)
    materset(1).put(w.CONDUCTIVITY,30)
    materset(1).put(w.DISSIP_TE,0)
    materset(1).put(w.DISSIP_TQ,0)
    
    # Contact material parameters
    
    defineTmCoulombContactMaterial(materset, idMatCI, peno_CI, peta_CI, rf_RT*facProf_CI, mu_CI, mu_CI, hc0_CI, HV_CI, w_CI)
    
    # Finite element properties

    prp1 = w.ElementProperties(w.TmVolume2DElement)
    prp1.put(w.CAUCHYMECHVOLINTMETH,w.VES_CMVIM_SRIPR)
    prp1.put(w.STIFFMETHOD,w.STIFF_ANALYTIC)
    prp1.put(w.MATERIAL,1)
    app.addProperty(prp1)
    
    # Contact finite element properties
    
    fctT_Tool = w.PieceWiseLinearFunction()
    fctT_Tool.setData(time0, 1.0)
    fctT_Tool.setData(time1, 1.0)
    fctT_Tool.setData(time2, 1.0)
    ci = defineTmRDContElePropInter(eleType_CI, idMatCI, w.STIFF_ANALYTIC, w.AIC_ONCE, T_Tool, fctT_Tool, interactionset, idMatCI, surf_RT, groups['Solid'], setSN_CI)
    
    # Elements for surface heat flux

    prp2 = w.ElementProperties(w.NodHeatFlux2DElement)
    heat = w.NodInteraction(2)
    heat.push(groups['FSInterface'])
    heat.addProperty(prp2)
    interactionset.add(heat)

    # Top sliding RT boundary conditions

    AmplDispX_RT = 1-0.2
    fct_dispX = w.PieceWiseLinearFunction()
    fct_dispX.setData(time0, 0.0)
    fct_dispX.setData(time1, 0.0)
    fct_dispX.setData(time2, 1.0)
    
    AmplDispY_RT = -HS*2.5e-2
    fct_dispY = w.PieceWiseLinearFunction()
    fct_dispY.setData(time0, 0.0)
    fct_dispY.setData(time1, 1.0)
    fct_dispY.setData(time2, 1.0)

    loadingset.define(surf_RT, w.Field1D(w.TX,w.RE), AmplDispX_RT, fct_dispX)
    loadingset.define(surf_RT, w.Field1D(w.TY,w.RE), AmplDispY_RT, fct_dispY)

    # Other boundary conditions

    initcondset.define(groups['Solid'],w.Field1D(w.TO,w.AB),300)
    loadingset.define(groups['FSInterface'],w.Field1D(w.TX,w.RE))
    loadingset.define(groups['FSInterface'],w.Field1D(w.TY,w.RE))

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
    
    # Parameters for MetaforStandalone
    if metaforStandalone:
        # Time step manager
        tsm.setInitialTime(time0, dt0)
        tsm.setNextTime(time1, nbArchi1, dtMax1)
        tsm.setNextTime(time2, nbArchi2, dtMax2)

    # Parameters for FSPC
    if not metaforStandalone:
        input['interacT'] = heat
        input['FSInterface'] = groups['FSInterface']
        input['exporter'] = gmsh.GmshExport('metafor/output.msh',metafor)
        input['exporter'].addDataBaseField([w.TO])
        input['exporter'].binary = True

    return metafor
