-- Problem Parameters

Problem = {}
Problem.autoRemeshing = false
Problem.verboseOutput = false
Problem.simulationTime = math.huge
Problem.id = 'WCompNewtonNoT'

-- FSPC Parameters

Problem.interface = 'FSInterface'
Problem.maxFactor = 100

-- Mesh Parameters

Problem.Mesh = {}
Problem.Mesh.alpha = 1.2
Problem.Mesh.omega = 0.5
Problem.Mesh.gamma = 0.6
Problem.Mesh.hchar = 0.01
Problem.Mesh.gammaFS = 0.2
Problem.Mesh.addOnFS = false
Problem.Mesh.minAspectRatio = 1e-2
Problem.Mesh.keepFluidElements = true
Problem.Mesh.deleteFlyingNodes = false
Problem.Mesh.deleteBoundElements = true
Problem.Mesh.laplacianSmoothingBoundaries = false
Problem.Mesh.boundingBox = {-0.8,-1,0.8,1}
Problem.Mesh.exclusionZones = {}

Problem.Mesh.remeshAlgo = 'GMSH'
Problem.Mesh.mshFile = 'geometryF.msh'
Problem.Mesh.exclusionGroups = {'Poly1','Poly2','Poly3','Poly4'}
Problem.Mesh.ignoreGroups = {}

-- Extractor Parameters

Problem.Extractors = {}

Problem.Extractors[0] = {}
Problem.Extractors[0].kind = 'GMSH'
Problem.Extractors[0].writeAs = 'NodesElements'
Problem.Extractors[0].outputFile = 'pfem/output.msh'
Problem.Extractors[0].whatToWrite = {'p','velocity'}
Problem.Extractors[0].timeBetweenWriting = math.huge

Problem.Extractors[1] = {}
Problem.Extractors[1].kind = 'Global'
Problem.Extractors[1].whatToWrite = 'mass'
Problem.Extractors[1].outputFile = 'mass.txt'
Problem.Extractors[1].timeBetweenWriting = math.huge

-- Material Parameters

Problem.Material = {}
Problem.Material.K0p = 1
Problem.Material.K0 = 100
Problem.Material.mu = 1e-4
Problem.Material.gamma = 0
Problem.Material.rhoStar = 1e-3

-- Solver Parameters

Problem.Solver = {}
Problem.Solver.id = 'CDS_dpdt'
Problem.Solver.securityCoeff = 0.2

Problem.Solver.adaptDT = true
Problem.Solver.maxDT = math.huge
Problem.Solver.initialDT = math.huge
Problem.Solver.maxRemeshDT = math.huge

-- Momentum Continuity Equation

Problem.Solver.MomEq = {}
Problem.Solver.ContEq = {}
Problem.Solver.ContEq.pExt = 0
Problem.Solver.MomEq.bodyForce = {0,0}
Problem.Solver.ContEq.stabilization = 'Meduri'

-- Momentum Continuity BC

Problem.IC = {}
Problem.Solver.MomEq.BC = {}
Problem.Solver.ContEq.BC = {}
Problem.Solver.MomEq.BC['FSInterfaceVExt'] = true

function Problem.IC:initStates(pos)
	return {0,0,0,Problem.Material.rhoStar,0,0}
end

function Problem.Solver.MomEq.BC:Poly1V(pos,t)
	return 0,0
end

function Problem.Solver.MomEq.BC:Poly2V(pos,t)
	return 0,0
end

function Problem.Solver.MomEq.BC:Poly3V(pos,t)
	return 0,0
end

function Problem.Solver.MomEq.BC:Poly4V(pos,t)
	return 0,0
end

function Problem.Solver.MomEq.BC:BorderV(pos,t)
	return 0,0
end

function Problem.Solver.MomEq.BC:OutletP(pos,t)
	return 0
end

function Problem.Solver.MomEq.BC:InletV(pos,t)

	local tmax = 1
	local acc = 10

	if (t<tmax) then
		return acc,0
	else
		return 0,0
	end
end
