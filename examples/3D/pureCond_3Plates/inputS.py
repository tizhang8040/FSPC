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
# - setSinglePass: True for Defo-Defo contact interactions
#(double pass: artefact in contact force computation if different mesh size between master and slave)
#(double pass: to be verified but seems heat flux will be counted twice)
# - N-R residual computation: mechanical-Method4; thermal-Method3
# - time integration scheme": staggered; mechanical-(Damped)AlphaGeneralized; thermal-Trapezoidal
# - solver: DSS

# importation modules
from wrap import *
import toolbox.gmsh as gmsh
import numpy as np
import os, math
#import sys.path


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
    dy7 = 0
    
    dx6 = rh*sinC
    dy6 = rh*(1-cosC)
    
    dx5 = rInt
    dy5 = rh*(1-cosB)
    
    dx4 = rExt-rf*(1-sinB)
    dy4 = dx4*tanB - (rInt*tanB-dy5)
    
    dx2 = rExt
    dy2 = dy4 + rf*cosB
        
    dx3 = rExt - rf*(1-sinA)
    dy3 = dy2 - rf*sinA
    
    dx1 = rExt
    dy1 = dy2 + hCyld
    # sets
    pointset = geometry.getPointSet()
    curveset = geometry.getCurveSet()
    wireset = geometry.getWireSet()
    surfaceset = geometry.getSurfaceSet()
    # -- Point --
    #generatrix
    pointset.define(idGeo+1, X0-dx1, Y0+dy1, Z0)
    pointset.define(idGeo+2, X0-dx2, Y0+dy2, Z0)
    pointset.define(idGeo+3, X0-dx3, Y0+dy3, Z0)
    pointset.define(idGeo+4, X0-dx4, Y0+dy4, Z0)
    pointset.define(idGeo+5, X0-dx5, Y0+dy5, Z0)
    pointset.define(idGeo+6, X0-dx6, Y0+dy6, Z0)
    pointset.define(idGeo+7, X0-dx7, Y0+dy7, Z0)
    #axis
    pointset.define(idGeo+51, X0, Y0, Z0)
    pointset.define(idGeo+52, X0, Y0+dy1, Z0)

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
# run the simulation with Metafor/FSPC
metaforStandalone = False


def getMetafor(_parameters={}):
    global _metafor
    if _metafor == None:
        _metafor = buildDomain(_parameters)
        
    return _metafor


# construction of parameters dictionary
def getParameters(_parameters={}):

    parameters = {}
    
    # interaction mode for interface heat transfert
    parameters['heatMode'] = "contact"	#contact/sticking
    
    # material
    parameters['materType'] = "linear"	#linear/nonlinear
    parameters['alpha_FeTop'] = 0.0	#10.0e-6  #thermal expansion coefficient #K-1
    parameters['alpha_Al']    = 0.0	#23.5e-6  #thermal expansion coefficient #K-1
    parameters['alpha_FeBot'] = 0.0	#16.5e-6  #thermal expansion coefficient #K-1
    
    # contact-RTFeTop
    parameters['peno_RTFeTop'] = 1.5e3
    parameters['hc0_RTFeTop'] = 250.0	#nominal thermal conductance under conduction #mW/(mm^2 K) #150.#45
    parameters['w_RTFeTop'] = 1.25	#0.95#1.25  #exponent coeficient #1.0/1.35
    parameters['HV_RTFeTop'] = 500.0	#Vickers's material hardness #MPa #900.#1500.
    parameters['setSN_RTFeTop'] = True
    if parameters['heatMode'] == "contact":
        # contact-FeTopAl
        parameters['peno_FeTopAl'] = 1.5e3
        parameters['hc0_FeTopAl'] = 350.0	#nominal thermal conductance under conduction #mW/(mm^2 K) #150.#45
        parameters['w_FeTopAl'] = 0.0		#0.95#1.25  #exponent coeficient #1.0/1.35
        parameters['HV_FeTopAl'] = 0.0		#50.0  #Vickers's material hardness #MPa #900./1500./230.
        parameters['setSN_FeTopAl'] = True
        # contact-AlFeBot
        parameters['peno_AlFeBot'] = 1.5e3
        parameters['hc0_AlFeBot'] = 100.0	#nominal thermal conductance under conduction #mW/(mm^2 K) #150.#45
        parameters['w_AlFeBot'] = 0.0		#0.95#1.25  #exponent coeficient #1.0/1.35
        parameters['HV_AlFeBot'] = 0.0		#300.0  #Vickers's material hardness #MPa #900.#1500.#300.
        parameters['setSN_AlFeBot'] = True
    
    # time step manager & N-R & time integration scheme
    # simulation
    parameters['time0'] = 0.  # initial time (sec)
    parameters['time1'] = parameters['time0']+32.0  # next time (sec) / total time 17.0 (sec)
    parameters['time2'] = parameters['time1']+8.0   # next time (sec) / total time 19.0 (sec)
    parameters['time3'] = parameters['time2']+12.0  # next time (sec) / total time 37.0 (sec)
    parameters['time4'] = parameters['time3']+78.0  # next time (sec) / total time 65.0 (sec)
    #parameters['time5'] = parameters['time4']+45.0   # next time (sec) / total time 67.0 (sec)
    parameters['dt0'] = 5.e-5  # initial time step (sec)
    parameters['dtMax1'] = 5.e-1  # maximal time step (sec)
    parameters['dtMax2'] = 5.e-1  # maximal time step (sec)
    parameters['dtMax3'] = 5.e-1  # maximal time step (sec)
    parameters['dtMax4'] = 5.e-1  # maximal time step (sec)
    #parameters['dtMax5'] = 5.e-1  # maximal time step (sec)
    # archiving
    parameters['nbArchi1'] = 20  # number of archives
    parameters['nbArchi2'] = 5   # number of archives
    parameters['nbArchi3'] = 5  # number of archives
    parameters['nbArchi4'] = 30  # number of archives
    #parameters['nbArchi5'] = 30   # number of archives
        
    # BCs
    parameters['endSurfBC'] = 'none'  # back/front/back and front/none
    parameters['latSurfBC'] = 'none'  # lat tot/lat par/none
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
    idGeo_FeBot = 200
    idGeo_RT = 1000
    # material & interaction
    idMatFA_FeTop = 1
    idMatFA_Al = 2
    idMatFA_FeBot = 3
    idMatCI_RTFeTop = 11
    idMatCI_FeTopAl = 12
    idMatCI_AlFeBot = 13
    # group
    idGrpGeo_FeTop = 0		#1-6: surfaces; 7: volume
    idGrpGeo_Al = 100
    idGrpGeo_FeBot = 200
    idGrpCI = 1000
    idGrpBC = 2000
    idGrpArchi = 3000

    # geometry & caracteristic lengths & number of elements
    #geometrical extension to avoid rebuild contact system each time step
    if parameters['heatMode'] == "contact":
        gExten = 0.;
    elif parameters['heatMode'] == "sticking":
        gExten = 0.;
    gX0_RT = 20.0		# [mm]
    gZ0_RT = 1.e-3		# pre-inter-penetration RT-FeTop
    # rigid tool
    rExt_RT = 10.0		# [mm]
    rf_RT = rExt_RT/10		# radius of the rounded fillet
    hCyld_RT = 3.0
    rInt_RT = rExt_RT/10
    angCone_RT = 1.0  		# [degree]
    # Fe-top
    Lx_FeTop = 110.		# [mm]
    Ly_FeTop = 80.		# [mm]
    Lz_FeTop = 1.5		# [mm]
    facLy1_FeTop = 1.0		# factor for the Ly1 length
    Ly1_FeTop = facLy1_FeTop*rExt_RT
    # Al-middle
    Lx_Al = Lx_FeTop+gExten
    Ly_Al = Ly_FeTop+gExten
    Lz_Al = 3.0
    # Fe-bottom
    Lx_FeBot = Lx_Al+gExten
    Ly_FeBot = Ly_Al+gExten
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
    X0_RT = X0-Lx_FeTop/2+rExt_RT+gX0_RT
    Y0_RT = Y0
    Z0_RT = Z0_FeTop+Lz_FeTop-gZ0_RT

    # material & element & interaction
    # Fe-top
    eleType_FeTop = TmVolume3DElement
    if parameters['materType'] == "linear":
        matType_FeTop = TmElastHypoMaterial
    elif parameters['materType'] == "nonlinear":
        matType_FeTop = TmEvpIsoHHypoMaterial
        K_FeTop = 1300.0				#Krupkowsky Isotropic Hardening parameter [MPa]
        epsPl0_FeTop = 0.05				#Krupkowsky Isotropic Hardening parameter
        n_FeTop = 0.48					#Krupkowsky Isotropic Hardening parameter
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
    if parameters['materType'] == "linear":
        matType_Al = TmElastHypoMaterial
    elif parameters['materType'] == "nonlinear":
        matType_Al = TmEvpIsoHHypoMaterial
        sigY0_Al = 45.0					#Saturated Isotropic Hardening parameter [MPa]
        Q_Al = 46.866					#Saturated Isotropic Hardening parameter [MPa]
        ksi_Al = 12.0					#Saturated Isotropic Hardening parameter
    rho_Al = 2.71e-9					# [T/mm3]
    E_Al = 69.e3					# [MPa]
    nu_Al = 0.33
    alpha_Al = parameters['alpha_Al']			#thermal expansion coefficient [K-1]
    k_Al = 229.0					#thermal conductivity [mW/(mm.K)]
    cp_Al = 900.0e6					#specific heat capacity [mJ/(T.K)]
    betaE_Al = 0.0  #coefficient d'énergie élastique dissipée sous forme de chaleur
    betaP_Al = 0.0  #coefficient d'énergie plastique dissipée sous forme de chaleur (Taylor–Quinney coefficient)
    # Fe-bottom
    eleType_FeBot = TmVolume3DElement
    if parameters['materType'] == "linear":
        matType_FeBot = TmElastHypoMaterial
    elif parameters['materType'] == "nonlinear":
        matType_FeBot = TmEvpIsoHHypoMaterial
        K_FeBot = 1350.0				#Krupkowsky Isotropic Hardening parameter [MPa]
        epsPl0_FeBot = 0.04				#Krupkowsky Isotropic Hardening parameter
        n_FeBot = 0.47					#Krupkowsky Isotropic Hardening parameter
    rho_FeBot = 8.00e-9					# [T/mm3]
    E_FeBot = 193.e3					# [MPa]
    nu_FeBot = 0.27
    alpha_FeBot = parameters['alpha_FeBot']		#thermal expansion coefficient [K-1]
    k_FeBot = 16.3					#thermal conductivity [mW/(mm.K)]
    cp_FeBot = 500.0e6					#specific heat capacity [mJ/(T.K)]
    betaE_FeBot = 0.0  #coefficient d'énergie élastique dissipée sous forme de chaleur
    betaP_FeBot = 0.0  #coefficient d'énergie plastique dissipée sous forme de chaleur (Taylor–Quinney coefficient)
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
    if parameters['heatMode'] == "contact":
        # contact-FeTopAl
        eleType_FeTopAl = TmContact3DElement
        muSta_FeTopAl = 0.45
        muDyn_FeTopAl = muSta_FeTopAl			#for numerical stability
        peno_FeTopAl = parameters['peno_FeTopAl']
        peta_FeTopAl = peno_FeTopAl*muSta_FeTopAl	#for the same accuracy
        facProf_FeTopAl = 0.65				#thickness factor to compute the depth from which contact is detected
        hc0_FeTopAl = parameters['hc0_FeTopAl']
        HV_FeTopAl = parameters['HV_FeTopAl']
        w_FeTopAl = parameters['w_FeTopAl']
        setSN_FeTopAl = parameters['setSN_FeTopAl']
        # contact-AlFeBot
        eleType_AlFeBot = TmContact3DElement
        muSta_AlFeBot = 0.45
        muDyn_AlFeBot = muSta_AlFeBot			#for numerical stability
        peno_AlFeBot = parameters['peno_AlFeBot']
        peta_AlFeBot = peno_AlFeBot*muSta_AlFeBot	#for the same accuracy
        facProf_AlFeBot = 0.65  			#thickness factor to compute the depth from which contact is detected
        hc0_AlFeBot = parameters['hc0_AlFeBot']
        HV_AlFeBot = parameters['HV_AlFeBot']
        w_AlFeBot = parameters['w_AlFeBot']
        setSN_AlFeBot = parameters['setSN_AlFeBot']
    
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
    showFig = True
    
    
    # BCs
    AmplDispX_RT = abs(X0_RT)*2
    AmplDispZ_RT = -parameters['uz_RT']
    '''
    if parameters['RT'] == 'cone':
        AmplDispZ_RT = -(Z0_RT-Z0_FeTop-Lz_FeTop)-hCone_RT-parameters['uz_RT']
    elif parameters['RT'] == 'cylinder':
        AmplDispZ_RT = -(Z0_RT-Z0_FeTop-Lz_FeTop)-parameters['uz_RT']
    '''
    T_Tool = 273.15+1200.0				# [K]
    T0_plates = 273.15+20.0				# [K]


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
    if parameters['materType'] == "linear":
        # Fe-top
        tmElastHypoMater(materset, idMatFA_FeTop, E_FeTop, nu_FeTop, alpha_FeTop, k_FeTop, cp_FeTop, rho_FeTop, betaE_FeTop, betaP_FeTop)
        defineElePropInter(eleType_FeTop, idMatFA_FeTop, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_FeTop, groups['FeTop'])
        # Al-middle
        tmElastHypoMater(materset, idMatFA_Al, E_Al, nu_Al, alpha_Al, k_Al, cp_Al, rho_Al, betaE_Al, betaP_Al)
        defineElePropInter(eleType_Al, idMatFA_Al, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_Al, groups['AlS'])
        # Fe-bottom
        tmElastHypoMater(materset, idMatFA_FeBot, E_FeBot, nu_FeBot, alpha_FeBot, k_FeBot, cp_FeBot, rho_FeBot, betaE_FeBot, betaP_FeBot)
        defineElePropInter(eleType_FeBot, idMatFA_FeBot, STIFF_ANALYTIC, VES_CMVIM_SRIPR, interset, idMatFA_FeBot, groups['FeBot'])    
    elif parameters['materType'] == "nonlinear":
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
    fctT_Tool.setData(parameters['time0'], 1.0)
    fctT_Tool.setData(parameters['time1'], 1.0)
    fctT_Tool.setData(parameters['time2'], 1.0)
    fctT_Tool.setData(parameters['time3'], 1.0)
    fctT_Tool.setData(parameters['time4'], 1.0)
    #fctT_Tool.setData(parameters['time5'], 1.0)
    if muSta_RTFeTop == 0.0:
        defineTmFrictionlessContactMaterial(materset, idMatCI_RTFeTop, peno_RTFeTop, rf_RT*facProf_RTFeTop, hc0_RTFeTop, HV_RTFeTop, w_RTFeTop)
    else:
        defineTmCoulombContactMaterial(materset, idMatCI_RTFeTop, peno_RTFeTop, peta_RTFeTop, rf_RT*facProf_RTFeTop, muSta_RTFeTop, muDyn_RTFeTop, hc0_RTFeTop, HV_RTFeTop, w_RTFeTop)
    ci_RTFeTop = defineTmRDContactElePropInter(eleType_RTFeTop, idMatCI_RTFeTop, STIFF_ANALYTIC, AIC_ONCE, T_Tool, fctT_Tool, interset, idMatCI_RTFeTop, surf_RT, groups['FeTop_Up'], setSN_RTFeTop)
    if parameters['heatMode'] == "contact":    
        # contact-FeTop/Al
        defineTmCoulombContactMaterial(materset, idMatCI_FeTopAl, peno_FeTopAl, peta_FeTopAl, Lz_Al*facProf_FeTopAl, muSta_FeTopAl, muDyn_FeTopAl, hc0_FeTopAl, HV_FeTopAl, w_FeTopAl)
        ci_FeTopAl = defineDDContactElePropInter(eleType_FeTopAl, idMatCI_FeTopAl, STIFF_ANALYTIC, AIC_ONCE, interset, idMatCI_FeTopAl, groups['AlS_Up'], groups['FeTop_Down'], setSN_FeTopAl, True)
        # contact-Al/FeBot
        defineTmCoulombContactMaterial(materset, idMatCI_AlFeBot, peno_AlFeBot, peta_AlFeBot, Lz_FeBot*facProf_AlFeBot, muSta_AlFeBot, muDyn_AlFeBot, hc0_AlFeBot, HV_AlFeBot, w_AlFeBot)
        ci_AlFeBot = defineDDContactElePropInter(eleType_AlFeBot, idMatCI_AlFeBot, STIFF_ANALYTIC, AIC_ONCE, interset, idMatCI_AlFeBot, groups['FeBot_Up'], groups['AlS_Down'], setSN_AlFeBot, True)
        
    # Elements for surface heat flux
    prp200 = ElementProperties(NodHeatFlux3DElement)
    heat = NodInteraction(200)
    heat.push(groups['FSInterface'])
    heat.addProperty(prp200)
    interset.add(heat)


    # BCs
    # initial BCs
    initialconditionset = metafor.getInitialConditionSet()
    initialconditionset.define(groups['FeTop'], Field1D(TO,AB), 0.0)
    initialconditionset.define(groups['AlS'],    Field1D(TO,AB), 0.0)
    initialconditionset.define(groups['FeBot'], Field1D(TO,AB), 0.0)
    initialconditionset.define(groups['FeTop'], Field1D(TO,RE), T0_plates)
    initialconditionset.define(groups['AlS'],    Field1D(TO,RE), T0_plates)
    initialconditionset.define(groups['FeBot'], Field1D(TO,RE), T0_plates)
    # loading BCs
    loadingset=domain.getLoadingSet()
    # fixed end-displacement
    # rigid tool
    fct_dispZ = PieceWiseLinearFunction()
    fct_dispZ.setData(parameters['time0'], 0.0)
    fct_dispZ.setData(parameters['time1'], 55.0/100)
    fct_dispZ.setData(parameters['time2'], 100.0/100)
    fct_dispZ.setData(parameters['time3'], 60.0/100)
    fct_dispZ.setData(parameters['time4'], 60.0/100)
    #fct_dispZ.setData(parameters['time5'], 60.0/100)
    
    fct_dispX = PieceWiseLinearFunction()
    fct_dispX.setData(parameters['time0'], 0.0)
    fct_dispX.setData(parameters['time1'], 0.0)
    fct_dispX.setData(parameters['time2'], 0.0)
    fct_dispX.setData(parameters['time3'], 0.0)
    fct_dispX.setData(parameters['time4'], 100.0/100)
    #fct_dispX.setData(parameters['time5'], 100/100)

    loadingset.define(surf_RT, Field1D(TX,RE), AmplDispX_RT, fct_dispX)
    loadingset.define(surf_RT, Field1D(TY,RE))
    loadingset.define(surf_RT, Field1D(TZ,RE), AmplDispZ_RT, fct_dispZ)
    # bottom
    loadingset.define(groups['FeBot_Down'], Field1D(TX,RE))
    loadingset.define(groups['FeBot_Down'], Field1D(TY,RE))
    loadingset.define(groups['FeBot_Down'], Field1D(TZ,RE))
    # up lines
    #loadingset.define(grpBC_FeTopUp,  Field1D(TZ,RE))
    # end surfaces
    if parameters['endSurfBC'] == 'back':
        loadingset.define(grpSurf6_FeTop, Field1D(TX,RE))
        loadingset.define(grpSurf6_Al,    Field1D(TX,RE))
        loadingset.define(grpSurf6_FeBot, Field1D(TX,RE))
    elif parameters['endSurfBC'] == 'front':
        loadingset.define(grpSurf4_FeTop, Field1D(TX,RE))
        loadingset.define(grpSurf4_Al,    Field1D(TX,RE))
        loadingset.define(grpSurf4_FeBot, Field1D(TX,RE))
    elif parameters['endSurfBC'] == 'back and front':
        loadingset.define(grpSurf6_FeTop, Field1D(TX,RE))
        loadingset.define(grpSurf6_Al,    Field1D(TX,RE))
        loadingset.define(grpSurf6_FeBot, Field1D(TX,RE))
        loadingset.define(grpSurf4_FeTop, Field1D(TX,RE))
        loadingset.define(grpSurf4_Al,    Field1D(TX,RE))
        loadingset.define(grpSurf4_FeBot, Field1D(TX,RE))
    elif parameters['endSurfBC'] == 'none':
        pass
    # lateral surfaces
    if parameters['latSurfBC'] == 'lat tot':
        loadingset.define(grpSurf5_FeTop, Field1D(TX,RE))
        loadingset.define(grpSurf5_Al,    Field1D(TX,RE))
        loadingset.define(grpSurf5_FeBot, Field1D(TX,RE))
    elif parameters['latSurfBC'] == 'lat par':
        loadingset.define(grpBC_FeTopExt, Field1D(TX,RE))
        loadingset.define(grpBC_AlExt,    Field1D(TX,RE))
        loadingset.define(grpBC_FeBotExt, Field1D(TX,RE))
    elif parameters['latSurfBC'] == 'none':
        pass
    # FSInterfaces
    loadingset.define(groups['FSInterface'],Field1D(TX,RE))
    loadingset.define(groups['FSInterface'],Field1D(TY,RE))
    loadingset.define(groups['FSInterface'],Field1D(TZ,RE))
    # contact interfaces
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
    
    # parameters for metaforStandalone
    if metaforStandalone:
        # time step manager
        tsm.setInitialTime(parameters['time0'], parameters['dt0'])
        tsm.setNextTime(parameters['time1'], parameters['nbArchi1'], parameters['dtMax1'])
        tsm.setNextTime(parameters['time2'], parameters['nbArchi2'], parameters['dtMax2'])
        tsm.setNextTime(parameters['time3'], parameters['nbArchi3'], parameters['dtMax3'])
        tsm.setNextTime(parameters['time4'], parameters['nbArchi4'], parameters['dtMax4'])
    # parameters for FSPC
    if not metaforStandalone:
        _parameters['interacT'] = heat
        _parameters['FSInterface'] = groups['FSInterface']
        _parameters['exporter'] = gmsh.GmshExport('metafor/output.msh',metafor)
        _parameters['exporter'].addDataBaseField([TO])



    return metafor
