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
Problem.Mesh.alpha = 1e3
Problem.Mesh.omega = 0.5
Problem.Mesh.gamma = 0.6
Problem.Mesh.hchar = 0.02
Problem.Mesh.gammaFS = 0.2
Problem.Mesh.addOnFS = true
Problem.Mesh.minAspectRatio = 1e-2
Problem.Mesh.keepFluidElements = true
Problem.Mesh.deleteFlyingNodes = false
Problem.Mesh.deleteBoundElements = true
Problem.Mesh.laplacianSmoothingBoundaries = false
Problem.Mesh.boundingBox = {-1,-2.625,1,2.375}
Problem.Mesh.exclusionZones = {}

Problem.Mesh.remeshAlgo = 'GMSH'
Problem.Mesh.mshFile = 'geometryF.msh'
Problem.Mesh.exclusionGroups = {'Polytope'}
Problem.Mesh.localHcharGroups = {'FSInterface'}
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
Problem.Material.mu = 1e-5
Problem.Material.K0p = 1
Problem.Material.gamma = 0
Problem.Material.K0 = 100
Problem.Material.rhoStar = 1e-6

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

function Problem.Solver.MomEq.BC:ReservoirV(pos,t)
	return 0,0
end

function Problem.Solver.MomEq.BC:PolytopeV(pos,t)
	return 0,0
end

function Problem.Solver.MomEq.BC:OutletP(pos,t)
	return 0
end

function Problem.Solver.MomEq.BC:InletVEuler(pos,t)

	local amax = -4e5
	local tmax = 2.5e-4
	local r = math.abs(pos[1])

	if (t<tmax) then
		return 0,amax*(1-r^2)
	else
		return 0,0
	end
end

function Problem.Mesh:computeHcharFromDistance(pos,t,dist)

	local hchar = Problem.Mesh.hchar
	return hchar+dist*0.1
end