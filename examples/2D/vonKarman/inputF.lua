-- Problem Parameters

Problem = {}
Problem.verboseOutput = true
Problem.autoRemeshing = false
Problem.simulationTime = math.huge
Problem.id = 'IncompNewtonNoT'

-- Mesh Parameters

Problem.Mesh = {}
Problem.Mesh.alpha = 1e3
Problem.Mesh.omega = 0.7
Problem.Mesh.gamma = 0.9
Problem.Mesh.hchar = 4e-3
Problem.Mesh.gammaFS = 0.3
Problem.Mesh.addOnFS = true
Problem.Mesh.minAspectRatio = 1e-3
Problem.Mesh.keepFluidElements = true
Problem.Mesh.deleteFlyingNodes = false
Problem.Mesh.deleteBoundElements = false
Problem.Mesh.laplacianSmoothingBoundaries = false
Problem.Mesh.boundingBox = {0,0,1.5,0.41}
Problem.Mesh.exclusionZones = {}

Problem.Mesh.remeshAlgo = 'GMSH'
Problem.Mesh.mshFile = 'geometryF.msh'
Problem.Mesh.exclusionGroups = {'Polytope'}
--Problem.Mesh.localHcharGroups = {'Polytope','Wall','Inlet'}
Problem.Mesh.localHcharGroups = {}
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
Problem.Material.mu = 1
Problem.Material.gamma = 0
Problem.Material.rho = 1000

-- Solver Parameters

Problem.Solver = {}
Problem.Solver.id = 'FracStep'

Problem.Solver.adaptDT = true
Problem.Solver.maxDT = math.huge
Problem.Solver.initialDT = math.huge
Problem.Solver.coeffDTDecrease = math.huge
Problem.Solver.coeffDTincrease = math.huge

-- Momentum Continuity Equation

Problem.Solver.MomContEq = {}
Problem.Solver.MomContEq.residual = 'U_P'
Problem.Solver.MomContEq.nlAlgo = 'Picard'
Problem.Solver.MomContEq.sparseSolverPstep = 'CG'

Problem.Solver.MomContEq.pExt = 0
Problem.Solver.MomContEq.maxIter = 25
Problem.Solver.MomContEq.gammaFS = 0.5
Problem.Solver.MomContEq.minRes = 1e-8
Problem.Solver.MomContEq.cgTolerance = 1e-16
Problem.Solver.MomContEq.bodyForce = {0,0}

-- Momentum Continuity BC

Problem.IC = {}
Problem.Solver.MomContEq.BC = {}
Problem.Solver.MomContEq.BC['FSInterfaceVExt'] = true

function Problem.IC.initStates(x,y,z)
    return {0,0,0}
end

function Problem.Solver.MomContEq.BC.WallV(x,y,z,t)
    return 0,0
end

function Problem.Solver.MomContEq.BC.PolytopeV(x,y,z,t)
    return 0,0
end

function Problem.Solver.MomContEq.BC.OutletP(x,y,z,t)
    return 0
end

function Problem.Solver.MomContEq.BC.InletVEuler(x,y,z,t)

    local vbar = 1
    local tmax = 2
    local H = 0.41/2
    local vx = 1.5*vbar*y*(2*H-y)/(H*H)

	if (t<tmax) then
        return vx*(1-math.cos(math.pi*t/2))/2,0
    else
        return vx,0
    end
end

function Problem.Mesh.computeHcharFromDistance(x,y,z,t,dist)

	local hchar = Problem.Mesh.hchar
	return hchar+dist*0.1
end
