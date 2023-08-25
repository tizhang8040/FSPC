-- Problem Parameters

Problem = {}
Problem.verboseOutput = true
Problem.autoRemeshing = false
Problem.simulationTime = math.huge
Problem.id = 'Conduction'

-- Mesh Parameters

Problem.Mesh = {}
Problem.Mesh.alpha = 5.5
Problem.Mesh.omega = 0.8
Problem.Mesh.gamma = 0.3
Problem.Mesh.hchar = 1.2
Problem.Mesh.gammaFS = 0.2
Problem.Mesh.addOnFS = true
Problem.Mesh.minAspectRatio = 1e-3
Problem.Mesh.keepFluidElements = true
Problem.Mesh.deleteFlyingNodes = true
Problem.Mesh.deleteBoundElements = false
Problem.Mesh.laplacianSmoothingBoundaries = false
Problem.Mesh.boundingBox = {-1000,-1000,-1000,1000,1000,1000}
Problem.Mesh.exclusionZones = {}

Problem.Mesh.localHcharGroups = {'FSInterface'}

Problem.Mesh.remeshAlgo = 'GMSH'
Problem.Mesh.mshFile = 'geometryF.msh'
Problem.Mesh.exclusionGroups = {}
Problem.Mesh.ignoreGroups = {}

-- Extractor Parameters

Problem.Extractors = {}

Problem.Extractors[0] = {}
Problem.Extractors[0].kind = 'GMSH'
Problem.Extractors[0].writeAs = 'NodesElements'
Problem.Extractors[0].outputFile = 'pfem/output.msh'
Problem.Extractors[0].whatToWrite = {'T'}
Problem.Extractors[0].timeBetweenWriting = math.huge

Problem.Extractors[1] = {}
Problem.Extractors[1].kind = 'Global'
Problem.Extractors[1].whatToWrite = 'mass'
Problem.Extractors[1].outputFile = 'mass.txt'
Problem.Extractors[1].timeBetweenWriting = math.huge

-- Material Parameters

Problem.Material = {}
Problem.Material.mu = 1.3e-9
Problem.Material.gamma = 0
Problem.Material.rho = 2.71e-9
Problem.Material.epsRad = 0
Problem.Material.sigmaRad = 5.670374419e-11
Problem.Material.R = 8.31446261815324
Problem.Material.alphaLin = 0.0002
Problem.Material.DgammaDT = 0
Problem.Material.Tinf = 300
Problem.Material.cp = 900e6
Problem.Material.DmuDT = 0
Problem.Material.DcpDT = 0
Problem.Material.DkDT = 0
Problem.Material.k = 200
Problem.Material.h = 1

Problem.Material.Tr = 933.5                     -- K^(-1)
Problem.Material.Tl = 933.5+15                  -- K^(-1) 482
Problem.Material.Ts = 933.5-15                  -- K^(-1)
Problem.Material.Lm = 400e9                     -- 400 kJ/kg = 400 * 10^9 mm^2/s^2 

-- Boussinesq PC
Problem.Material.C = 1e-3                       -- 
Problem.Material.eps = 1e-3                     -- 1e-3
Problem.Material.epsRad = 0                     -- 
Problem.Material.sigmaRad = 5.670374419e-11     -- 5.67e-8 W/(m^2 * K^4) = 5.67e-11 ton/(s^3 * K^4)

Problem.Material.mu_s = 1.3e-9
Problem.Material.Tr_s = 930                     
Problem.Material.k_s = 200
Problem.Material.Dk_sDT = 0
Problem.Material.cp_s = 0
Problem.Material.Dcp_sDT = 0
Problem.Material.Dmu_sDT = 0
Problem.Material.gamma_s = 0
Problem.Material.Dgamma_sDT = 0

-- Solver Parameters

Problem.Solver = {}
Problem.Solver.id = 'PSPG'

Problem.Solver.adaptDT = true
Problem.Solver.maxDT = math.huge
Problem.Solver.initialDT = math.huge
Problem.Solver.coeffDTDecrease = math.huge
Problem.Solver.coeffDTincrease = math.huge
Problem.Solver.solveHeatFirst = true

-- Momentum Continuity Equation

Problem.Solver.MomContEq = {}
Problem.Solver.MomContEq.nlAlgo   = 'Picard'
Problem.Solver.MomContEq.residual = 'U_P'
Problem.Solver.MomContEq.minRes   = 1e-6
Problem.Solver.MomContEq.maxIter  = 10
Problem.Solver.MomContEq.bodyForce = {0, -9810, 0}

Problem.Solver.MomContEq.gammaFS = 0.5
Problem.Solver.MomContEq.cgTolerance = 1e-16
Problem.Solver.MomContEq.sparseSolverLib = 'MKL'

-- Heat Equation

Problem.Solver.HeatEq = {}
Problem.Solver.HeatEq.residual = 'Ax_f'
Problem.Solver.HeatEq.nlAlgo = 'Picard'
Problem.Solver.HeatEq.sparseSolver = 'CG'

Problem.Solver.HeatEq.maxIter = 25
Problem.Solver.HeatEq.minRes = 1e-6
Problem.Solver.HeatEq.cgTolerance = 1e-16

-- Heat Momentum Continuity BC

Problem.IC = {}
Problem.Solver.HeatEq.BC = {}
Problem.Solver.MomContEq.BC = {}
Problem.Solver.HeatEq.BC['FSInterfaceTExt'] = true

function Problem.IC.initStates(x,y,z)
	return {300}
end

function Problem.Solver.HeatEq.BC.WallQ(x,y,z,t)
    return 0,0,0
end

function Problem.Solver.HeatEq.BC.FreeSurfaceQ(x,y,z,t)
    return 0,0,0
end

function Problem.Solver.MomContEq.BC.WallV(x, y, z, t) 
	return 0, 0, 0
end

function Problem.Solver.MomContEq.BC.FSInterfaceV(x, y, z, t) 
	return 0, 0, 0
end
