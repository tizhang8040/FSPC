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

metaforStandalone = False

metafor = None
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
    L = 40 		        #   solid domain length
    HF = 3.0  		        #   fluid domain height
    HS = 1.5  		        #   solid domain height
    W = 30                      #   solid domain width
    
    rExt_RT = 12.275            #   rigid tool radius
    rf_RT = rExt_RT/10          #   radius of the rounded fillet
    hCyld_RT = 3.               #
    rInt_RT = rExt_RT/10        #
    angCone_RT = 1.0            #   degree
    
    X0_RT = rExt_RT+3           #   origin x of the rigid tool center
    Y0_RT = HF+HS-1.e-3         #   origin y of the rigid tool center with a small inter-penetration
    Z0_RT = W/2                 #   origin z of the rigid tool center
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
    time3 = time2+1.0                   #   horizontal displacement time (tool)
    dt0 = 5.e-5
    dtMax1 = 5.e-1
    dtMax2 = 5.e-1
    dtMax3 = 5.e-1
    nbArchi1 = 10
    nbArchi2 = 50
    nbArchi3 = 2
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

    #mshFile = os.path.join(os.path.dirname(__file__),'geometryS_hex.msh')
    mshFile = os.path.join(os.path.dirname(__file__),'geometryS_traHex.msh')
    importer = gmsh.GmshImport(mshFile,domain)
    groups = importer.groups
    importer.execute()

    # ===== tianyu =======
    # Defines the rigid tool
    surf_RT = createRigidCone(geometry, idGeo_RT, X0_RT, Y0_RT, Z0_RT, rExt_RT, rInt_RT, hCyld_RT, rf_RT, angCone_RT)
    # ====================

    # Defines the solid domain

    app = w.FieldApplicator(1)
    app.push(groups['Solid'])
    interactionset.add(app)

    # Solid material parameters

    materset.define(1,w.TmElastHypoMaterial)
    materset(1).put(w.ELASTIC_MODULUS,2e5)
    materset(1).put(w.THERM_EXPANSION,0)
    materset(1).put(w.HEAT_CAPACITY,500e6)
    materset(1).put(w.MASS_DENSITY,7.85e-9)
    materset(1).put(w.CONDUCTIVITY,30)
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
    prp1.put(w.MATERIAL,1)
    app.addProperty(prp1)

    # ===== tianyu =======
    # Contact finite element properties
    
    fctT_Tool = w.PieceWiseLinearFunction()
    fctT_Tool.setData(time0, 1.0)
    fctT_Tool.setData(time1, 1.0)
    fctT_Tool.setData(time2, 1.0)
    fctT_Tool.setData(time3, 1.0)
    ci = defineTmRDContactElePropInter(eleType_CI, idMatCI, w.STIFF_ANALYTIC, w.AIC_ONCE, T_Tool, fctT_Tool, interactionset, idMatCI, surf_RT, groups['Solid'], setSN_CI)
    # ====================

    # Elements for surface heat flux

    #prp2 = w.ElementProperties(w.NodHeatFlux3DElement)
    prp2 = w.ElementProperties(w.NodHeatFlux3DElement)
    heat = w.NodInteraction(2)
    heat.push(groups['FSInterface'])
    heat.addProperty(prp2)
    interactionset.add(heat)

#    # Top heat flux boundary conditions
#    prp3 = w.ElementProperties(w.NodTriangleHeatFlux3DElement)
#    inter = w.NodInteraction(3)
#    inter.push(groups['Top'])
#    inter.addProperty(prp3)
#    interactionset.add(inter)
#
#    for i in range(groups['Top'].getNumberOfMeshPoints()):
#        node = groups['Top'].getMeshPoint(i)
#        x = node.getValue(w.Field1D(w.TX, w.AB))
#        y = node.getValue(w.Field1D(w.TY, w.AB))
#        norm = np.array([2.0,2.0])-np.array([x,y])
#        flux =  10000000*norm/np.linalg.norm(norm)
#        inter.setNodVector(node,*flux)


#    fun = lambda x : np.exp(-np.square(x-2)/np.square(0.1)/2)/100000
#    F = w.PythonOneParameterFunction(fun)
#
#    prp3 = w.ElementProperties(w.NodTriangleHeatFlux3DElement)
#    prp3.put(w.HEATEL_VALUE,1e3)
#    prp3.depend(w.HEATEL_VALUE,F,w.Field1D(w.TX,w.AB))
#
#    inter = w.HeatInteraction(3)
#    inter.push(groups['Top'])
#    inter.addProperty(prp3)
#    interactionset.add(inter)

    # ===== tianyu =======
    # Top sliding RT boundary conditions

    AmplDispX_RT = 9  # horizontal displacement
    fct_dispX = w.PieceWiseLinearFunction()
    fct_dispX.setData(time0, 0.0)
    fct_dispX.setData(time1, 0.0)
    fct_dispX.setData(time2, 1.0)
    fct_dispX.setData(time3, 1.0)

    AmplDispY_RT = -HS*2.5e-2  # vertical displacement
    fct_dispY = w.PieceWiseLinearFunction()
    fct_dispY.setData(time0, 0.0)
    fct_dispY.setData(time1, 1.0)
    fct_dispY.setData(time2, 1.0)
    fct_dispY.setData(time3, 1.0)

    loadingset.define(surf_RT, w.Field1D(w.TX,w.RE), AmplDispX_RT, fct_dispX)
    loadingset.define(surf_RT, w.Field1D(w.TY,w.RE), AmplDispY_RT, fct_dispY)
    loadingset.define(surf_RT, w.Field1D(w.TZ,w.RE))
    # ====================

    # Initial and boundary conditions

    initcondset.define(groups['Solid'],w.Field1D(w.TO,w.AB),0.0)
    initcondset.define(groups['Solid'],w.Field1D(w.TO,w.RE),300)
    initcondset.define(groups['FSInterface'],w.Field1D(w.TO,w.AB),0.0)
    initcondset.define(groups['FSInterface'],w.Field1D(w.TO,w.RE),300)


    # ==== Tianyu =====
    loadingset.define(groups['FSInterface'],w.Field1D(w.TX,w.RE))
    loadingset.define(groups['FSInterface'],w.Field1D(w.TY,w.RE))
    loadingset.define(groups['FSInterface'],w.Field1D(w.TZ,w.RE))
    # =================

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
        tsm.setNextTime(time3, nbArchi3, dtMax3)

    # Parameters for FSPC
    if not metaforStandalone:
        parm['interacT'] = heat
        parm['FSInterface'] = groups['FSInterface']
        parm['exporter'] = gmsh.GmshExport('metafor/output.msh',metafor)
        parm['exporter'].addDataBaseField([w.TO])
    
    return metafor
