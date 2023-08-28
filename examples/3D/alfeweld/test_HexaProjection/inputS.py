# -*- coding: utf-8 -*-

# - 3D complet model
# - rigid tool: cone without rotation
# - material type: TmElastHypoMaterial/TmEvpIsoHHypoMaterial
# - hardening type for stainless steel: KrupkowskyIsotropicHardening
# - hardening type for AA1050: SaturatedIsotropicHardening
# - contact type: positively unilateral
# - contact law:  Frictionless & TmCoulomb
# - method used to compute the stiffness matrix: STIFF_ANALYTIC
# - method used to integrate the Cauchy stresses: VES_CMVIM_SRIPR
# - method to take into account the contact surface: AIC_ONCE
# - method to take into account the sticking surface: SA_ONCE
# - setSinglePass: True for Defo-Defo contact interactions
#(double pass: artefact in contact force computation if different mesh size between master and slave)
#(double pass: to be verified but seems heat flux will be counted twice)
# - N-R residual computation: mechanical-Method4; thermal-Method3
# - time integration scheme": staggered; mechanical-(Damped)AlphaGeneralized; thermal-Trapezoidal
# - solver: DSS

# importation modules
from wrap import *
import toolbox.gmsh as gmsh
import os, math


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
    curveset.add(Line(idGeo+51, pointset(idGeo+52), pointset(idGeo+51)))
    curveset.add(Line(idGeo+1, pointset(idGeo+1), pointset(idGeo+2)))
    curveset.add(Arc (idGeo+2, pointset(idGeo+2), pointset(idGeo+3), pointset(idGeo+4)))
    curveset.add(Line(idGeo+3, pointset(idGeo+4), pointset(idGeo+5)))
    curveset.add(Arc (idGeo+4, pointset(idGeo+5), pointset(idGeo+6), pointset(idGeo+7)))
    wireset.add(Wire(idGeo+1, [curveset(idGeo+1),curveset(idGeo+2),curveset(idGeo+3),curveset(idGeo+4)]))
    # -- Surface --
    surfRC = surfaceset.add(RevolutionSurface(idGeo+1, curveset(idGeo+51), wireset(idGeo+1)))
    
    return surfRC

def cstSpringMater(materset, materIndex, kMech):
    matSetI = materset.define(materIndex,ConstantSpringMaterial)
    matSetI.put(SPRING_FK, kMech)

def tmElastHypoMater(materset, materIndex, E, Nu, Alpha, K, C, Rho, dissipTe, dissipTq):
    matSetI = materset.define(materIndex,TmElastHypoMaterial)
    matSetI.put(ELASTIC_MODULUS, E)
    matSetI.put(POISSON_RATIO,   Nu)
    matSetI.put(THERM_EXPANSION, Alpha)
    matSetI.put(CONDUCTIVITY,    K)
    matSetI.put(MASS_DENSITY,    Rho)
    matSetI.put(HEAT_CAPACITY,   C)
    matSetI.put(DISSIP_TE,       dissipTe)
    matSetI.put(DISSIP_TQ,       dissipTq)

def defineTmFrictionlessContactMaterial(materset, materIndex, normalPeno, prof, hNom, hardVickers, wExp):
    materset.define(materIndex, TmFrictionlessContactMaterial)
    materset(materIndex).put(PEN_NORMALE,    normalPeno)
    materset(materIndex).put(PROF_CONT,      prof)
    materset(materIndex).put(CTM_H_NOMINAL,  hNom)
    materset(materIndex).put(CTM_HARDNESS,   hardVickers)
    materset(materIndex).put(CTM_EXPONENT_E, wExp)

def defineTmCoulombContactMaterial(materset, materIndex, normalPeno, tangPeno, prof, statCoefFrot, dynCoefFrot, hNom, hardVickers, wExp):
    materset.define(materIndex, TmCoulombContactMaterial)
    materset(materIndex).put(PEN_NORMALE,    normalPeno)
    materset(materIndex).put(PEN_TANGENT,    tangPeno)
    materset(materIndex).put(PROF_CONT,      prof)
    materset(materIndex).put(COEF_FROT_STA,  statCoefFrot)
    materset(materIndex).put(COEF_FROT_DYN,  dynCoefFrot)
    materset(materIndex).put(CTM_H_NOMINAL,  hNom)
    materset(materIndex).put(CTM_HARDNESS,   hardVickers)
    materset(materIndex).put(CTM_EXPONENT_E, wExp)

def defineElePropInter(eleType, materIndex, stiffMethod, mechIntegMethod, interset, interIndex, gObject):
    prpElem = ElementProperties(eleType)
    prpElem.put(MATERIAL,             materIndex)
    prpElem.put(STIFFMETHOD,          stiffMethod)
    prpElem.put(CAUCHYMECHVOLINTMETH, mechIntegMethod)
    
    fa = FieldApplicator(interIndex)
    fa.push(gObject)
    fa.addProperty(prpElem)
    interset.add(fa)

def defineTmStickElePropInter(eleType, materIndex, stiffMethod, SAMethod, kTher, interset, interIndex, gTool, gSlave):
    prpElem = ElementProperties(eleType)
    prpElem.put(MATERIAL,     materIndex)
    prpElem.put(STIFFMETHOD,  stiffMethod)
    prpElem.put(STICKINGAREA, SAMethod)
    prpElem.put(KTHER,        kTher)
    
    si = StickingInteraction(interIndex)
    si.setTool(gTool)
    si.push(gSlave)
    si.addProperty(prpElem)
    interset.add(si)
    
    return si

def defineTmRDContactElePropInter(eleType, materIndex, stiffMethod, AICMethod, T_Tool, fctT_Tool, interset, interIndex, gTool, gSlave, setSN):
    prpElem = ElementProperties(eleType)
    prpElem.put(MATERIAL,      materIndex)
    prpElem.put(STIFFMETHOD,   stiffMethod)
    prpElem.put(AREAINCONTACT, AICMethod)
    prpElem.put(TOOLTEMP,      T_Tool)
    prpElem.depend(TOOLTEMP, fctT_Tool, Field1D(TM,RE))
    
    rdci = RdContactInteraction(interIndex)
    rdci.setTool(gTool)
    rdci.push(gSlave)
    rdci.setSmoothNormals(setSN)
    rdci.addProperty(prpElem)
    interset.add(rdci)
    
    return rdci

def defineDDContactElePropInter(eleType, materIndex, stiffMethod, AICMethod, interset, interIndex, gTool, gSlave, setSN, setSP):
    prpElem = ElementProperties(eleType)
    prpElem.put(MATERIAL,      materIndex)
    prpElem.put(STIFFMETHOD,   stiffMethod)
    prpElem.put(AREAINCONTACT, AICMethod)
    
    ddci = DdContactInteraction(interIndex)
    ddci.setTool(gTool)
    ddci.push(gSlave)
    ddci.setSmoothNormals(setSN)
    if setSP:
        ddci.setSinglePass()    
    ddci.addProperty(prpElem)
    interset.add(ddci)
    
    return ddci



# interface with Metafor
_metafor = None
# for parallel computation
StrVectorBase.useTBB()
StrMatrixBase.useTBB()
ContactInteraction.useTBB()
#ValuesManager.useTBB()  # lead to "Segmentation fault" if used with "IFGeoPointValueExtractor"


def getMetafor(_parameters={}):
    global _metafor
    if _metafor == None:
        _metafor = buildDomain(_parameters)
        
    return _metafor


# construction of parameters dictionary
def getParameters(_parameters={}):

    parameters = {}
    
    # material
    parameters['alpha_FeTop'] = 0.0	#10.0e-6  #thermal expansion coefficient #K-1
    parameters['alpha_Al']    = 0.0	#23.5e-6  #thermal expansion coefficient #K-1
    
    # contact-RTFeTop
    parameters['peno_RTFeTop'] = 1.5e3
    parameters['hc0_RTFeTop'] = 250.0	#nominal thermal conductance under conduction #mW/(mm^2 K) #150.#45
    parameters['w_RTFeTop'] = 1.25	#0.95#1.25  #exponent coeficient #1.0/1.35
    parameters['HV_RTFeTop'] = 500.0	#Vickers's material hardness #MPa #900.#1500.
    parameters['setSN_RTFeTop'] = True
    # sticking-FeTopAl
    parameters['stickKMech_FeTopAl'] = 1.e5
    parameters['stickKTher_FeTopAl'] = 1.e5
    
    # time step manager & N-R & time integration scheme
    # simulation
    parameters['time0'] = 0.  # initial time (sec)
    parameters['time1'] = parameters['time0']+32.0  # next time (sec) / total time 17.0 (sec)
    parameters['dt0'] = 0.01  # initial time step (sec)
    parameters['dtMax1'] = 0.5  # maximal time step (sec)
    parameters['nbArchi1'] = 20  # number of archives
        
    # BCs
    parameters['uz_RT'] = 0.10
    
    
    parameters.update(_parameters)

    return parameters


def buildDomain(_parameters={}):
    # get parameters
    parameters  = getParameters(_parameters)
    print("parameters = ", parameters)
    
    # get sets
    metafor = Metafor()
    domain = metafor.getDomain()
    geometry = domain.getGeometry()
    materset = domain.getMaterialSet()
    hardset = domain.getMaterialLawSet()
    interset = domain.getInteractionSet()    
    initcondset = metafor.getInitialConditionSet()
    loadingset = domain.getLoadingSet()    
    tsm = metafor.getTimeStepManager()
    mim = metafor.getMechanicalIterationManager()
    tim = metafor.getThermalIterationManager()
    sm  = metafor.getSolverManager()
        
    # set dimension and solver
    geometry.setDim3D()    
    try:
        if metaSol == 'DSS':
            solver = DSSolver()
            sm.setSolver(solver)
            print("OK : DSSolver found.")
        elif metaSol == 'MUMPS':
            solver = MUMPSolver()
            sm.setSolver(solver)
            print("OK : MUMPSolver found.")
        #elif metaSol == 'ISS':
        #    solver = ISSolver()
        #    solver.setRestart(d['issRestart'])
        #    solver.setItMax(d['issItMax'])
        #    solver.setTol(d['issTol'])
        #    solver.useILUT(d['issILUT'])
        #    sm.setSolver(solver)
        #    print("OK : ISSolver found.")
    except NameError:
        print("Asked solver not found. Use of default skyline solver instead.")
        pass

    # fixed parameters
    # ids
    idGeo_FeTop = 0
    idGeo_Al = 100
    idGeo_RT = 1000
    # material & interaction
    idMatFA_FeTop = 1
    idMatFA_Al = 2
    idMatCI_RTFeTop = 11
    idMatSI_FeTopAl = 12
    # group
    idGrpGeo_FeTop = 0		#1-6: surfaces; 7: volume
    idGrpGeo_Al = 100
    idGrpCI = 1000
    idGrpBC = 2000
    idGrpArchi = 3000

    # geometry & caracteristic lengths & number of elements
    #geometrical extension to avoid rebuild contact system each time step
    gX0_RT = 20.0		# [mm]
    gZ0_RT = 1.e-3		# pre-inter-penetration RT-FeTop
    # rigid tool
    rExt_RT = 10.0		# [mm]
    rf_RT = rExt_RT/10		# radius of the rounded fillet
    hCyld_RT = 3.0
    rInt_RT = rExt_RT/10
    angCone_RT = 0.5  		# [degree]
    # Fe-top
    Lx_FeTop = 120.		# [mm]
    Ly_FeTop = 80.		# [mm]
    Lz_FeTop = 1.5		# [mm]
    # Al-middle
    Lx_Al = Lx_FeTop
    Ly_Al = Ly_FeTop
    Lz_Al = 3.0
    # origin
    X0 = 0.0
    Y0 = 0.0
    # Al-middle
    Z0_Al = 0.0
    # Fe-top
    Z0_FeTop = Z0_Al+Lz_Al
    # rigid tool
    X0_RT = X0-Lx_FeTop/2+rExt_RT+gX0_RT
    Y0_RT = Y0
    Z0_RT = Z0_FeTop+Lz_FeTop-gZ0_RT

    # material & element & interaction
    # Fe-top
    eleType_FeTop = TmVolume3DElement
    matType_FeTop = TmElastHypoMaterial
    rho_FeTop = 7.85e-9					# [T/mm3]
    E_FeTop = 210.e3					# [MPa]
    nu_FeTop = 0.28
    alpha_FeTop = parameters['alpha_FeTop']		#thermal expansion coefficient [K-1]
    k_FeTop = 30.0					#thermal conductivity [mW/(mm.K)]
    cp_FeTop = 500.0e6					#specific heat capacity [mJ/(T.K)]
    betaE_FeTop = 0.0  #coefficient d'énergie élastique dissipée sous forme de chaleur
    betaP_FeTop = 0.0  #coefficient d'énergie plastique dissipée sous forme de chaleur (Taylor–Quinney coefficient)
    # Al-middle
    eleType_Al = TmVolume3DElement
    matType_Al = TmElastHypoMaterial
    rho_Al = 2.71e-9					# [T/mm3]
    E_Al = 69.e3					# [MPa]
    nu_Al = 0.33
    alpha_Al = parameters['alpha_Al']			#thermal expansion coefficient [K-1]
    k_Al = 229.0					#thermal conductivity [mW/(mm.K)]
    cp_Al = 900.0e6					#specific heat capacity [mJ/(T.K)]
    betaE_Al = 0.0  #coefficient d'énergie élastique dissipée sous forme de chaleur
    betaP_Al = 0.0  #coefficient d'énergie plastique dissipée sous forme de chaleur (Taylor–Quinney coefficient)
    # contact-RTFeTop
    eleType_RTFeTop = TmContact3DElement
    muSta_RTFeTop = 0.0
    muDyn_RTFeTop = muSta_RTFeTop			#for numerical stability
    peno_RTFeTop = parameters['peno_RTFeTop']
    peta_RTFeTop = peno_RTFeTop*muSta_RTFeTop		#for the same accuracy
    facProf_RTFeTop = 0.95				#thickness factor to compute the depth from which contact is detected
    hc0_RTFeTop = parameters['hc0_RTFeTop']
    HV_RTFeTop = parameters['HV_RTFeTop']
    w_RTFeTop = parameters['w_RTFeTop']
    setSN_RTFeTop = parameters['setSN_RTFeTop']
    # sticking-FeTopAl
    stickKMech_FeTopAl = parameters['stickKMech_FeTopAl']
    stickEleType_FeTopAl = TmSticking3DElement
    stickKTher_FeTopAl = parameters['stickKTher_FeTopAl']
    
    # N-R & time integration scheme & others
    # N-R iteration manager
    # mechanical
    mechResMethod = Method4ResidualComputation()
    #mechMaxIter = 10
    mechResTol = 1.0E-4
    # thermal
    therResMethod = Method3ResidualComputation()
    #therMaxIter = 10
    therResTol = 1.0E-6
    # time integration scheme
    staggeredTimeInteg = 1  # ?(attention matrice de raideur doit être numérique !)?
    mechTimeInteg = AlphaGeneralizedTimeIntegration	# AlphaGeneralized/QuasiStatic
    therTimeInteg = TrapezoidalThermalTimeIntegration	# Mpg/Trapezoidal
    #mechAlphaM = -0.97#-1.0
    #mechAlphaF = (mechAlphaM+1)/3 #+0.02#+0.01#+0.0
    #mechTheta = 1.0
    #therTheta = 1.0
    # others
    metaSol = 'DSS'
    
    
    # BCs
    AmplDispZ_RT = -parameters['uz_RT']

    T_Tool = 273.15+1200.0				# [K]
    T0_plates = 300.0#273.15+20.0				# [K]


###----------------------------------------------------------------------------------------------------------###
    # geometry
    # rigid tool
    surf_RT = createRigidCone(geometry, idGeo_RT, X0_RT, Y0_RT, Z0_RT, rExt_RT, rInt_RT, hCyld_RT, rf_RT, angCone_RT)
    # Imports the mesh
    mshFile = os.path.join(os.path.dirname(__file__),'geometryS.msh')
    importer = gmsh.GmshImport(mshFile,domain)
    groups = importer.groups
    importer.execute()

    # constitutive law/materials properties & element & interaction (material property - FE linking)
    # Fe-top
    tmElastHypoMater(materset, idMatFA_FeTop, E_FeTop, nu_FeTop, alpha_FeTop, k_FeTop, cp_FeTop, rho_FeTop, betaE_FeTop, betaP_FeTop)
    defineElePropInter(eleType_FeTop, idMatFA_FeTop, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_FeTop, groups['FeTop'])
    # Al-middle
    tmElastHypoMater(materset, idMatFA_Al, E_Al, nu_Al, alpha_Al, k_Al, cp_Al, rho_Al, betaE_Al, betaP_Al)
    defineElePropInter(eleType_Al, idMatFA_Al, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_Al, groups['Al'])
    # contact-RT/FeTop
    fctT_Tool = PieceWiseLinearFunction()
    fctT_Tool.setData(parameters['time0'], 1.0)
    fctT_Tool.setData(parameters['time1'], 1.0)
    if muSta_RTFeTop == 0.0:
        defineTmFrictionlessContactMaterial(materset, idMatCI_RTFeTop, peno_RTFeTop, rf_RT*facProf_RTFeTop, hc0_RTFeTop, HV_RTFeTop, w_RTFeTop)
    else:
        defineTmCoulombContactMaterial(materset, idMatCI_RTFeTop, peno_RTFeTop, peta_RTFeTop, rf_RT*facProf_RTFeTop, muSta_RTFeTop, muDyn_RTFeTop, hc0_RTFeTop, HV_RTFeTop, w_RTFeTop)
    ci_RTFeTop = defineTmRDContactElePropInter(eleType_RTFeTop, idMatCI_RTFeTop, STIFF_ANALYTIC, AIC_ONCE, T_Tool, fctT_Tool, interset, idMatCI_RTFeTop, surf_RT, groups['FeTop_Up'], setSN_RTFeTop)  
    # sticking-FeTop/Al
    cstSpringMater(materset, idMatSI_FeTopAl, stickKMech_FeTopAl)
    defineTmStickElePropInter(stickEleType_FeTopAl, idMatSI_FeTopAl, STIFF_ANALYTIC, SA_ONCE, stickKTher_FeTopAl, interset, idMatSI_FeTopAl, groups['Al_Up'], groups['FeTop_Down'])


    # BCs
    # initial BCs
    initialconditionset = metafor.getInitialConditionSet()
    initialconditionset.define(groups['FeTop'], Field1D(TO,AB), 0.0)
    initialconditionset.define(groups['Al'],    Field1D(TO,AB), 0.0)
    initialconditionset.define(groups['FeTop'], Field1D(TO,RE), T0_plates)
    initialconditionset.define(groups['Al'],    Field1D(TO,RE), T0_plates)
    # loading BCs
    loadingset=domain.getLoadingSet()
    # fixed end-displacement
    # rigid tool
    fct_dispZ = PieceWiseLinearFunction()
    fct_dispZ.setData(parameters['time0'], 0.0)
    fct_dispZ.setData(parameters['time1'], 1.0)

    loadingset.define(surf_RT, Field1D(TX,RE))
    loadingset.define(surf_RT, Field1D(TY,RE))
    loadingset.define(surf_RT, Field1D(TZ,RE), AmplDispZ_RT, fct_dispZ)
    # bottom
    loadingset.define(groups['Al_Down'], Field1D(TX,RE))
    loadingset.define(groups['Al_Down'], Field1D(TY,RE))
    loadingset.define(groups['Al_Down'], Field1D(TZ,RE))
    # contact interfaces
    loadingset.define(groups['Al_Up'],Field1D(TX,RE))
    loadingset.define(groups['Al_Up'],Field1D(TY,RE))
    loadingset.define(groups['Al_Up'],Field1D(TZ,RE))
    loadingset.define(groups['FeTop_Down'],Field1D(TX,RE))
    loadingset.define(groups['FeTop_Down'],Field1D(TY,RE))
    loadingset.define(groups['FeTop_Down'],Field1D(TZ,RE))


    #N-R & time integration scheme
    # N-R iteration manager
    # mechanical
    mim.setResidualComputationMethod(mechResMethod)
    #mim.setMaxNbOfIterations(mechMaxIter)
    mim.setResidualTolerance(mechResTol)
    #mim.setVerbose()
    # thermal
    tim.setResidualComputationMethod(therResMethod)
    #tim.setMaxNbOfIterations(therMaxIter)
    tim.setResidualTolerance(therResTol)
    #tim.setVerbose()    
    # time integration scheme
    # mechanical
    tiMech = mechTimeInteg(metafor)
    #tiMech.setTheta(mechTheta)
    #tiMech.setAlphaM(mechAlphaM)
    #tiMech.setAlphaF(mechAlphaF)
    # thermal
    tiTher = therTimeInteg(metafor)
    #tiTher.setTheta(therTheta)
    # staggered
    if not staggeredTimeInteg:
        ti = CoupledTmTimeIntegration(metafor)
    elif staggeredTimeInteg:
        ti = StaggeredTmTimeIntegration(metafor)
        ti.setIsAdiabatic(False)
        ti.setWithStressReevaluation(False)
    ti.setMechanicalTimeIntegration(tiMech)
    ti.setThermalTimeIntegration(tiTher)
    metafor.setTimeIntegration(ti)

    # time step iterations
    #tscm = NbOfStaggeredTmNRIterationsTimeStepComputationMethod(metafor)
    #tsm.setTimeStepComputationMethod(tscm)
    #tscm.setTimeStepDivisionFactor(2)
    #tscm.setNbOptiIte(25)

    # time step manager
    tsm.setInitialTime(parameters['time0'], parameters['dt0'])
    tsm.setNextTime(parameters['time1'], parameters['nbArchi1'], parameters['dtMax1'])



    return metafor
