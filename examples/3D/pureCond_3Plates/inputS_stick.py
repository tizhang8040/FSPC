# -*- coding: utf-8 -*-

# - 3D complet model
# - interaction mode for interface heat transfert: sticking
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

def tmFrictionlessContactMaterial(materset, materIndex, normalPeno, prof, hNom, hardVickers, wExp):
    materset.define(materIndex, TmFrictionlessContactMaterial)
    materset(materIndex).put(PEN_NORMALE,    normalPeno)
    materset(materIndex).put(PROF_CONT,      prof)
    materset(materIndex).put(CTM_H_NOMINAL,  hNom)
    materset(materIndex).put(CTM_HARDNESS,   hardVickers)
    materset(materIndex).put(CTM_EXPONENT_E, wExp)

def tmCoulombContactMaterial(materset, materIndex, normalPeno, tangPeno, prof, statCoefFrot, dynCoefFrot, hNom, hardVickers, wExp):
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

def defineTmStickElePropInter(eleType, materIndex, outTol, stiffMethod, SAMethod, kTher, interset, interIndex, gTool, gSlave):
    prpElem = ElementProperties(eleType)
    prpElem.put(MATERIAL,     materIndex)
    prpElem.put(OUTSIDETOL,   outTol)
    prpElem.put(STIFFMETHOD,  stiffMethod)
    prpElem.put(STICKINGAREA, SAMethod)
    prpElem.put(KTHER,        kTher)
    
    si = StickingInteraction(interIndex)
    si.setTool(gTool)
    si.push(gSlave)
    si.addProperty(prpElem)
    interset.add(si)
    
    return si

def defineNodHeatFluxElePropInter(eleType, interset, interIndex, gObject):
    prpElem = ElementProperties(eleType)
    
    ni = NodInteraction(interIndex)
    ni.push(gObject)
    ni.addProperty(prpElem)
    interset.add(ni)
    
    return ni

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

def defineConvectionElePropInter(eleType, stiffMethod, hCvt, TInf, interset, interIndex, gObjectList):
    prpElem = ElementProperties(eleType)
    prpElem.put(STIFFMETHOD, stiffMethod)
    prpElem.put(CONV_COEF,   hCvt)
    prpElem.put(TEMP_FLUIDE, TInf)
    
    li = LoadingInteraction(interIndex)
    for ig in range(len(gObjectList)):
        li.push(gObjectList[ig])
    li.addProperty(prpElem)
    interset.add(li)

def defineRadiationElePropInter(eleType, stiffMethod, cstBoltzmann, eRad, TInf, interset, interIndex, gObjectList):
    prpElem = ElementProperties(eleType)
    prpElem.put(STIFFMETHOD,    stiffMethod)
    prpElem.put(BOLTZMANN_CST,  cstBoltzmann)
    prpElem.put(RAY_EMISSIVITY, eRad)
    prpElem.put(RAY_TEMP_AMB,   TInf)
    
    li = LoadingInteraction(interIndex)
    for ig in range(len(gObjectList)):
        li.push(gObjectList[ig])
    li.addProperty(prpElem)
    interset.add(li)






# interface with Metafor
_metafor = None
# for parallel computation
StrVectorBase.useTBB()
StrMatrixBase.useTBB()
ContactInteraction.useTBB()
#ValuesManager.useTBB()  # lead to "Segmentation fault" if used with "IFGeoPointValueExtractor"


def getMetafor(_p={}):
    global _metafor
    if _metafor == None:
        _metafor = buildDomain(_p)
        
    return _metafor


# construction of parameters dictionary
def getParameters(_p={}):

    p = {}

    # material
    p['materType'] = "linear"	#linear/nonlinear
    p['alpha_FeTop'] = 0.0	#10.0e-6  #thermal expansion coefficient #K-1
    p['alpha_Al']    = 0.0	#23.5e-6  #thermal expansion coefficient #K-1
    p['alpha_FeBot'] = 0.0	#16.5e-6  #thermal expansion coefficient #K-1

    # contact-RTFeTop
    p['peno_RTFeTop'] = 1.5e3
    p['hc0_RTFeTop'] = 250.0	#nominal thermal conductance under conduction #mW/(mm^2 K) #150.#45
    p['w_RTFeTop'] = 1.25	#0.95#1.25  #exponent coeficient #1.0/1.35
    p['HV_RTFeTop'] = 500.0	#Vickers's material hardness #MPa #900.#1500.
    p['setSN_RTFeTop'] = True
    # sticking-FeTopAl
    p['stickKMech_FeTopAl'] = 1.e5	# to be determined
    p['stickKTher_FeTopAl'] = 1.e5	# to be determined
    # sticking-AlFeBot
    p['stickKMech_AlFeBot'] = 1.e5	# to be determined
    p['stickKTher_AlFeBot'] = 1.e5	# to be determined
    # convection
    p['cvtTop'] = True
    p['cvtLat'] = True
    p['cvtBot'] = True
    # radiation
    p['radTop'] = True
    p['radLat'] = True

    # time step manager & N-R & time integration scheme
    # simulation
    p['time0'] = 0.  			# initial time [s]
    p['time1'] = p['time0']+5.0  	# next time [s]
    p['time2'] = p['time1']+15.0   	# next time [s]
    p['time3'] = p['time2']+1.0   	# next time [s]

    # BCs
    p['uz_RT'] = 0.05
    p['FeTopUpEdgesBC'] = False
    p['endSurfBC'] = 'none'  		# back/front/back and front/none

    # run the simulation with Metafor/FSPC
    p['metaforStandalone'] = False

    p.update(_p)

    if p['metaforStandalone']:
        p['dt0'] = 0.01  		# initial time step [s]
        p['dtMax1'] = 0.01  		# maximal time step [s]
        p['dtMax2'] = 0.01  		# maximal time step [s]
        p['dtMax3'] = 0.01  		# maximal time step [s]
        # archiving
        p['nbArchi1'] = 10  		# number of archives
        p['nbArchi2'] = 30   		# number of archives
        p['nbArchi3'] = 2   		# number of archives


    p.update(_p)

    return p


def buildDomain(_p={}):
    # get parameters
    p = getParameters(_p)    
    print("parameters = ", p)
    
    # get sets
    metafor = Metafor()
    initcondset = metafor.getInitialConditionSet()
    tsm         = metafor.getTimeStepManager()
    mim         = metafor.getMechanicalIterationManager()
    tim         = metafor.getThermalIterationManager()
    sm          = metafor.getSolverManager()
    domain      = metafor.getDomain()
    geometry    = domain.getGeometry()
    materset    = domain.getMaterialSet()
    hardset     = domain.getMaterialLawSet()
    interset    = domain.getInteractionSet()
    loadingset  = domain.getLoadingSet()


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
    idGeo_RT = 1000
    # material & interaction
    idMatFA_FeTop = 1
    idMatFA_Al = 2
    idMatFA_FeBot = 3
    idMatCI_RTFeTop = 11
    idMatSI_FeTopAl = 12
    idMatSI_AlFeBot = 13    
    idNI_FSI = 14
    idCvt_FeTopAirUp = 21
    idCvt_FeTopAirLat = 22
    idCvt_AlAirLat = 23
    idCvt_FeBotAirLat = 24
    idCvt_FeBotSolDown = 25
    idRad_FeTopAirUp = 31
    idRad_FeTopAirLat = 32
    idRad_AlAirLat = 33
    idRad_FeBotAirLat = 34

    # geometry & caracteristic lengths & number of elements
    gZ0_RT = 1.e-3		# pre-inter-penetration RT-FeTop
    # rigid tool
    rExt_RT = 12.0		# [mm]
    rf_RT = rExt_RT/10		# radius of the rounded fillet
    hCyld_RT = 3.0
    rInt_RT = rExt_RT/10
    angCone_RT = 0.1  		# [degree]
    # Fe-top
    Lx_FeTop = 135.		# [mm]
    Ly_FeTop = 60.		# [mm]
    Lz_FeTop = 1.5		# [mm]
    # Al-middle
    Lx_Al = Lx_FeTop
    Ly_Al = Ly_FeTop
    Lz_Al = 3.0
    gx_AlF = 12.5		# gap in x direction between left edges of AlS-AlF
    # Fe-bottom
    Lx_FeBot = Lx_Al
    Ly_FeBot = Ly_Al
    Lz_FeBot = 5.0
    # origin
    X0 = 0.0
    Y0 = 0.0
    # Fe-bottom
    Z0_FeBot = 0.0
    # Al-middle
    Z0_Al = Z0_FeBot+Lz_FeBot
    # Fe-top
    Z0_FeTop = Z0_Al+Lz_Al
    # rigid tool
    X0_RT = X0-Lx_FeTop/2+rExt_RT+gx_AlF
    Y0_RT = Y0
    Z0_RT = Z0_FeTop+Lz_FeTop-gZ0_RT

    # material & element & interaction
    # Fe-top
    eleType_FeTop = TmVolume3DElement
    if p['materType'] == "linear":
        matType_FeTop = TmElastHypoMaterial
    elif p['materType'] == "nonlinear":
        matType_FeTop = TmEvpIsoHHypoMaterial
        K_FeTop = 1300.0				#Krupkowsky Isotropic Hardening parameter [MPa]
        epsPl0_FeTop = 0.05				#Krupkowsky Isotropic Hardening parameter
        n_FeTop = 0.48					#Krupkowsky Isotropic Hardening parameter
    rho_FeTop = 7.85e-9					#[T/mm3]
    E_FeTop = 210.e3					#[MPa]
    nu_FeTop = 0.28
    alpha_FeTop = p['alpha_FeTop']			#thermal expansion coefficient [K-1]
    k_FeTop = 30.0					#thermal conductivity [mW/(mm.K)]
    cp_FeTop = 500.0e6					#specific heat capacity [mJ/(T.K)]
    betaE_FeTop = 0.0  #coefficient d'énergie élastique dissipée sous forme de chaleur
    betaP_FeTop = 0.0  #coefficient d'énergie plastique dissipée sous forme de chaleur (Taylor–Quinney coefficient)
    # Al-middle
    eleType_Al = TmVolume3DElement
    if p['materType'] == "linear":
        matType_Al = TmElastHypoMaterial
    elif p['materType'] == "nonlinear":
        matType_Al = TmEvpIsoHHypoMaterial
        sigY0_Al = 45.0					#Saturated Isotropic Hardening parameter [MPa]
        Q_Al = 46.866					#Saturated Isotropic Hardening parameter [MPa]
        ksi_Al = 12.0					#Saturated Isotropic Hardening parameter
    rho_Al = 2.71e-9					#[T/mm3]
    E_Al = 69.e3					#[MPa]
    nu_Al = 0.33
    alpha_Al = p['alpha_Al']				#thermal expansion coefficient [K-1]
    k_Al = 229.0					#thermal conductivity [mW/(mm.K)]
    cp_Al = 900.0e6					#specific heat capacity [mJ/(T.K)]
    betaE_Al = 0.0  #coefficient d'énergie élastique dissipée sous forme de chaleur
    betaP_Al = 0.0  #coefficient d'énergie plastique dissipée sous forme de chaleur (Taylor–Quinney coefficient)
    # Fe-bottom
    eleType_FeBot = TmVolume3DElement
    if p['materType'] == "linear":
        matType_FeBot = TmElastHypoMaterial
    elif p['materType'] == "nonlinear":
        matType_FeBot = TmEvpIsoHHypoMaterial
        K_FeBot = 1350.0				#Krupkowsky Isotropic Hardening parameter [MPa]
        epsPl0_FeBot = 0.04				#Krupkowsky Isotropic Hardening parameter
        n_FeBot = 0.47					#Krupkowsky Isotropic Hardening parameter
    rho_FeBot = 8.00e-9					#[T/mm3]
    E_FeBot = 193.e3					#[MPa]
    nu_FeBot = 0.27
    alpha_FeBot = p['alpha_FeBot']			#thermal expansion coefficient [K-1]
    k_FeBot = 16.3					#thermal conductivity [mW/(mm.K)]
    cp_FeBot = 500.0e6					#specific heat capacity [mJ/(T.K)]
    betaE_FeBot = 0.0  #coefficient d'énergie élastique dissipée sous forme de chaleur
    betaP_FeBot = 0.0  #coefficient d'énergie plastique dissipée sous forme de chaleur (Taylor–Quinney coefficient)
    # contact-RTFeTop
    eleType_RTFeTop = TmContact3DElement
    muSta_RTFeTop = 0.0
    muDyn_RTFeTop = muSta_RTFeTop			#for numerical stability
    peno_RTFeTop = p['peno_RTFeTop']
    peta_RTFeTop = peno_RTFeTop*muSta_RTFeTop		#for the same accuracy
    facProf_RTFeTop = 0.95				#thickness factor to compute the depth from which contact is detected
    hc0_RTFeTop = p['hc0_RTFeTop']
    HV_RTFeTop = p['HV_RTFeTop']
    w_RTFeTop = p['w_RTFeTop']
    setSN_RTFeTop = p['setSN_RTFeTop']
    # sticking-FeTopAl
    stickEleType_FeTopAl = TmSticking3DElement
    stickKMech_FeTopAl = p['stickKMech_FeTopAl']
    stickKTher_FeTopAl = p['stickKTher_FeTopAl']
    stickOutTol_FeTopAl = 0.3
    # sticking-AlFeBot
    stickEleType_AlFeBot = TmSticking3DElement
    stickKMech_AlFeBot = p['stickKMech_AlFeBot']   
    stickKTher_AlFeBot = p['stickKTher_AlFeBot']
    stickOutTol_AlFeBot = 0.3
    # convections
    eleType_cvt = TmConvection3DElement
    hCvt_FeTopAirUp = 35.e-3				#thermal convection coefficient between solid and ambient air [mW/(mm2.K)]
    hCvt_FeTopAirLat = 7.5e-3				#thermal convection coefficient between solid and ambient air [mW/(mm2.K)]
    hCvt_AlAirLat = 15.e-3				#thermal convection coefficient between solid and ambient air [mW/(mm2.K)]
    hCvt_FeBotAirLat = 5.e-3				#thermal convection coefficient between solid and ambient air [mW/(mm2.K)]
    hCvt_FeBotSolDown = 2.5e-3				#thermal convection coefficient between solid and solid [mW/(mm2.K)]
    # radiaton
    eleType_rad = TmRayonnement3DElement
    cstBoltzmann = 5.67e-11				#Boltzmann radiation constant [mW/(mm2.K4)]
    eRad_FeTopAirUp = 0.50				#thermal radiation emissivity between solid and ambient air
    eRad_FeTopAirLat = 0.30				#thermal radiation emissivity between solid and ambient air
    eRad_AlAirLat = 0.05				#thermal radiation emissivity between solid and ambient air
    eRad_FeBotAirLat = 0.20				#thermal radiation emissivity between solid and ambient air
    
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
    staggeredTimeInteg = 1				#?(attention matrice de raideur doit être numérique !)?
    mechTimeInteg = AlphaGeneralizedTimeIntegration	#(Damped)AlphaGeneralized/QuasiStatic
    therTimeInteg = TrapezoidalThermalTimeIntegration	#Mpg/Trapezoidal
    #mechAlphaM = -0.97#-1.0
    #mechAlphaF = (mechAlphaM+1)/3 #+0.02#+0.01#+0.0
    #mechTheta = 1.0
    #therTheta = 1.0
    # others
    metaSol = 'DSS'
    showFig = True
    
    
    # BCs
    AmplDispX_RT = abs(X0_RT)*2
    AmplDispZ_RT = -p['uz_RT']
    T_Tool = 1480.#273.15+1200.0			#[K]
    T0_plates = 300.#273.15+20.0			#[K]
    T_air = 300.#273.15+20.0  				#convection-ambient air temperature [K]
    T_sol = 300.#273.15+20.0				#convection-backing solid temperature [K]


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
    if p['materType'] == "linear":
        # Fe-top
        tmElastHypoMater(materset, idMatFA_FeTop, E_FeTop, nu_FeTop, alpha_FeTop, k_FeTop, cp_FeTop, rho_FeTop, betaE_FeTop, betaP_FeTop)
        defineElePropInter(eleType_FeTop, idMatFA_FeTop, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_FeTop, groups['FeTop'])
        # Al-middle
        tmElastHypoMater(materset, idMatFA_Al, E_Al, nu_Al, alpha_Al, k_Al, cp_Al, rho_Al, betaE_Al, betaP_Al)
        defineElePropInter(eleType_Al, idMatFA_Al, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_Al, groups['AlS'])
        # Fe-bottom
        tmElastHypoMater(materset, idMatFA_FeBot, E_FeBot, nu_FeBot, alpha_FeBot, k_FeBot, cp_FeBot, rho_FeBot, betaE_FeBot, betaP_FeBot)
        defineElePropInter(eleType_FeBot, idMatFA_FeBot, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_FeBot, groups['FeBot'])    
    elif p['materType'] == "nonlinear":
        # Fe-top
        tmMaterYield(materset, idMatFA_FeTop, idMatFA_FeTop, E_FeTop, nu_FeTop, alpha_FeTop, k_FeTop, cp_FeTop, rho_FeTop, betaE_FeTop, betaP_FeTop)
        materLawIHKrupkowsky(hardset, idMatFA_FeTop, K_FeTop, epsPl0_FeTop, n_FeTop)
        defineElePropInter(eleType_FeTop, idMatFA_FeTop, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_FeTop, groups['FeTop'])
        # Al-middle
        tmMaterYield(materset, idMatFA_Al, idMatFA_Al, E_Al, nu_Al, alpha_Al, k_Al, cp_Al, rho_Al, betaE_Al, betaP_Al)
        materLawIHSat(hardset, idMatFA_Al, sigY0_Al, Q_Al, ksi_Al)
        defineElePropInter(eleType_Al, idMatFA_Al, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_Al, groups['AlS'])
        # Fe-bottom
        tmMaterYield(materset, idMatFA_FeBot, idMatFA_FeBot, E_FeBot, nu_FeBot, alpha_FeBot, k_FeBot, cp_FeBot, rho_FeBot, betaE_FeBot, betaP_FeBot)
        materLawIHKrupkowsky(hardset, idMatFA_FeBot, K_FeBot, epsPl0_FeBot, n_FeBot)
        defineElePropInter(eleType_FeBot, idMatFA_FeBot, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_FeBot, groups['FeBot'])
    
    # contact-RT/FeTop
    fctT_Tool = PieceWiseLinearFunction()
    fctT_Tool.setData(p['time0'], 1.0)
    fctT_Tool.setData(p['time1'], 1.0)
    fctT_Tool.setData(p['time2'], 1.0)
    fctT_Tool.setData(p['time3'], 1.0)
    tmFrictionlessContactMaterial(materset, idMatCI_RTFeTop, peno_RTFeTop, rf_RT*facProf_RTFeTop, hc0_RTFeTop, HV_RTFeTop, w_RTFeTop)
    ci_RTFeTop = defineTmRDContactElePropInter(eleType_RTFeTop, idMatCI_RTFeTop, STIFF_ANALYTIC, AIC_ONCE, T_Tool, fctT_Tool, interset, idMatCI_RTFeTop, surf_RT, groups['FeTop_UpCentral'], setSN_RTFeTop)
    # sticking-FeTop/Al
    cstSpringMater(materset, idMatSI_FeTopAl, stickKMech_FeTopAl)
    si_FeTopAl = defineTmStickElePropInter(stickEleType_FeTopAl, idMatSI_FeTopAl, stickOutTol_FeTopAl, STIFF_ANALYTIC, SA_ONCE, stickKTher_FeTopAl, interset, idMatSI_FeTopAl, groups['AlS_Up'], groups['FeTop_Down'])
    # sticking-Al/FeBot
    cstSpringMater(materset, idMatSI_AlFeBot, stickKMech_AlFeBot)
    si_AlFeBot = defineTmStickElePropInter(stickEleType_AlFeBot, idMatSI_AlFeBot, stickOutTol_AlFeBot, STIFF_ANALYTIC, SA_ONCE, stickKTher_AlFeBot, interset, idMatSI_AlFeBot, groups['FeBot_Up'], groups['AlS_Down'])    
    
    # convection
    if p['cvtTop']:
        defineConvectionElePropInter(eleType_cvt, STIFF_ANALYTIC, hCvt_FeTopAirUp,   T_air, interset, idCvt_FeTopAirUp,   [groups['FeTop_Up']])
    if p['cvtLat']:
        defineConvectionElePropInter(eleType_cvt, STIFF_ANALYTIC, hCvt_FeTopAirLat,  T_air, interset, idCvt_FeTopAirLat,  [groups['FeTop_Front'],groups['FeTop_Back'],groups['FeTop_Sides']])
        defineConvectionElePropInter(eleType_cvt, STIFF_ANALYTIC, hCvt_AlAirLat,     T_air, interset, idCvt_AlAirLat,     [groups['AlS_Front'],groups['AlS_Back'],groups['AlS_Sides']])
        defineConvectionElePropInter(eleType_cvt, STIFF_ANALYTIC, hCvt_FeBotAirLat,  T_air, interset, idCvt_FeBotAirLat,  [groups['FeBot_Front'],groups['FeBot_Back'],groups['FeBot_Sides']])
    if p['cvtBot']:
        defineConvectionElePropInter(eleType_cvt, STIFF_ANALYTIC, hCvt_FeBotSolDown, T_sol, interset, idCvt_FeBotSolDown, [groups['FeBot_Down']])
    # radiation
    if p['radTop']:
        defineRadiationElePropInter(eleType_rad, STIFF_ANALYTIC, cstBoltzmann, eRad_FeTopAirUp,  T_air, interset, idRad_FeTopAirUp,  [groups['FeTop_Up']])
    if p['radLat']:
        defineRadiationElePropInter(eleType_rad, STIFF_ANALYTIC, cstBoltzmann, eRad_FeTopAirLat, T_air, interset, idRad_FeTopAirLat, [groups['FeTop_Front'],groups['FeTop_Back'],groups['FeTop_Sides']])
        defineRadiationElePropInter(eleType_rad, STIFF_ANALYTIC, cstBoltzmann, eRad_AlAirLat,    T_air, interset, idRad_AlAirLat,    [groups['AlS_Front'],groups['AlS_Back'],groups['AlS_Sides']])
        defineRadiationElePropInter(eleType_rad, STIFF_ANALYTIC, cstBoltzmann, eRad_FeBotAirLat, T_air, interset, idRad_FeBotAirLat, [groups['FeBot_Front'],groups['FeBot_Back'],groups['FeBot_Sides']])
    
    # FSInterface heat flux transfer
    if not p['metaforStandalone']:
        ni_FSI = defineNodHeatFluxElePropInter(NodHeatFlux3DElement, interset, idNI_FSI, groups['FSInterface'])


    # BCs
    # initial BCs
    initcondset.define(groups['FeTop'], Field1D(TO,AB), 0.0)
    initcondset.define(groups['AlS'],   Field1D(TO,AB), 0.0)
    initcondset.define(groups['FeBot'], Field1D(TO,AB), 0.0)
    initcondset.define(groups['FeTop'], Field1D(TO,RE), T0_plates)
    initcondset.define(groups['AlS'],   Field1D(TO,RE), T0_plates)
    initcondset.define(groups['FeBot'], Field1D(TO,RE), T0_plates)
    # loading BCs
    # imposed rigid tool displacement
    #Z
    fct_dispZ = PieceWiseLinearFunction()
    fct_dispZ.setData(p['time0'], 0.0)
    fct_dispZ.setData(p['time1'], 100.0/100)
    fct_dispZ.setData(p['time2'], 80.0/100)
    fct_dispZ.setData(p['time3'], 80.0/100)
    #X
    fct_dispX = PieceWiseLinearFunction()
    fct_dispX.setData(p['time0'], 0.0)
    fct_dispX.setData(p['time1'], 0.0)
    fct_dispX.setData(p['time2'], 1.0)
    fct_dispX.setData(p['time3'], 1.0)
    #
    loadingset.define(surf_RT, Field1D(TX,RE), AmplDispX_RT, fct_dispX)
    loadingset.define(surf_RT, Field1D(TY,RE))
    loadingset.define(surf_RT, Field1D(TZ,RE), AmplDispZ_RT, fct_dispZ)
    # fixed interfaces
    loadingset.define(groups['FeBot_Up'],Field1D(TX,RE))
    loadingset.define(groups['FeBot_Up'],Field1D(TY,RE))
    loadingset.define(groups['FeBot_Up'],Field1D(TZ,RE))
    loadingset.define(groups['AlS_Down'],Field1D(TX,RE))
    loadingset.define(groups['AlS_Down'],Field1D(TY,RE))
    loadingset.define(groups['AlS_Down'],Field1D(TZ,RE))
    loadingset.define(groups['AlS_Up'],Field1D(TX,RE))
    loadingset.define(groups['AlS_Up'],Field1D(TY,RE))
    loadingset.define(groups['AlS_Up'],Field1D(TZ,RE))
    loadingset.define(groups['FeTop_Down'],Field1D(TX,RE))
    loadingset.define(groups['FeTop_Down'],Field1D(TY,RE))
    loadingset.define(groups['FeTop_Down'],Field1D(TZ,RE))
    # fixed FSInterfaces
    loadingset.define(groups['FSInterface'],Field1D(TX,RE))
    loadingset.define(groups['FSInterface'],Field1D(TY,RE))
    loadingset.define(groups['FSInterface'],Field1D(TZ,RE))
    # fixed displacement
    # bottom
    loadingset.define(groups['FeBot_Down'], Field1D(TX,RE))
    loadingset.define(groups['FeBot_Down'], Field1D(TY,RE))
    loadingset.define(groups['FeBot_Down'], Field1D(TZ,RE))
    # up lines
    if p['FeTopUpEdgesBC']:
        loadingset.define(groups['FeTop_UpEdges'],  Field1D(TZ,RE))
    # end surfaces
    if p['endSurfBC'] == 'back':
        loadingset.define(groups['FeTop_Back'], Field1D(TX,RE))
        loadingset.define(groups['AlS_Back'],   Field1D(TX,RE))
        loadingset.define(groups['FeBot_Back'], Field1D(TX,RE))
    elif p['endSurfBC'] == 'front':
        loadingset.define(groups['FeTop_Front'], Field1D(TX,RE))
        loadingset.define(groups['AlS_Front'],   Field1D(TX,RE))
        loadingset.define(groups['FeBot_Front'], Field1D(TX,RE))
    elif p['endSurfBC'] == 'back and front':
        loadingset.define(groups['FeTop_Back'], Field1D(TX,RE))
        loadingset.define(groups['AlS_Back'],   Field1D(TX,RE))
        loadingset.define(groups['FeBot_Back'], Field1D(TX,RE))
        loadingset.define(groups['FeTop_Front'], Field1D(TX,RE))
        loadingset.define(groups['AlS_Front'],   Field1D(TX,RE))
        loadingset.define(groups['FeBot_Front'], Field1D(TX,RE))
    elif p['endSurfBC'] == 'none':
        pass


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
    tscm = NbOfStaggeredTmNRIterationsTimeStepComputationMethod(metafor)
    tsm.setTimeStepComputationMethod(tscm)
    tscm.setTimeStepDivisionFactor(2)
    tscm.setNbOptiIte(25)
    
    # parameters for metaforStandalone
    if p['metaforStandalone']:
        # time step manager
        tsm.setInitialTime(p['time0'], p['dt0'])
        tsm.setNextTime(p['time1'], p['nbArchi1'], p['dtMax1'])
        tsm.setNextTime(p['time2'], p['nbArchi2'], p['dtMax2'])
        tsm.setNextTime(p['time3'], p['nbArchi3'], p['dtMax3'])
    # parameters for FSPC
    if not p['metaforStandalone']:
        _p['interacT'] = ni_FSI
        _p['FSInterface'] = groups['FSInterface']
        _p['exporter'] = gmsh.GmshExport('metafor/output.msh',metafor)
        _p['exporter'].addDataBaseField([TO])



    return metafor
