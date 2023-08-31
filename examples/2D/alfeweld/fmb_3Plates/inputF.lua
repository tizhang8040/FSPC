-- Problem Parameters

Problem = {}
Problem.verboseOutput = true
Problem.autoRemeshing = false
Problem.simulationTime = math.huge
Problem.id = 'BoussinesqPC'

-- Mesh Parameters

Problem.Mesh = {}
Problem.Mesh.alpha = 4.7
Problem.Mesh.omega = 0.7
Problem.Mesh.gamma = 0.4
Problem.Mesh.hchar = 0.5
Problem.Mesh.gammaFS = 0.2
Problem.Mesh.addOnFS = true
Problem.Mesh.minAspectRatio = 1e-3
Problem.Mesh.keepFluidElements = true
Problem.Mesh.deleteFlyingNodes = true
Problem.Mesh.deleteBoundElements = false
Problem.Mesh.laplacianSmoothingBoundaries = false
Problem.Mesh.boundingBox = {0,-10,1000,3.001}
Problem.Mesh.exclusionZones = {}

Problem.Mesh.localHcharGroups = {'FSInterface'}

Problem.Mesh.remeshAlgo = 'CGALConstrainedChew'
Problem.Mesh.mshFile = 'geometryF_refine.msh'
--Problem.Mesh.mshFile = 'geometryF_coarse.msh'
Problem.Mesh.exclusionGroups = {}
Problem.Mesh.ignoreGroups = {}

-- Extractor Parameters

Problem.Extractors = {}

Problem.Extractors[0] = {}
Problem.Extractors[0].kind = 'GMSH'
Problem.Extractors[0].writeAs = 'NodesElements'
Problem.Extractors[0].outputFile = 'pfem/output.msh'
Problem.Extractors[0].whatToWrite = {'T','fl','velocity','ke'}
Problem.Extractors[0].timeBetweenWriting = math.huge

Problem.Extractors[1] = {}
Problem.Extractors[1].kind = 'Global'
Problem.Extractors[1].whatToWrite = 'mass'
Problem.Extractors[1].outputFile = 'mass.txt'
Problem.Extractors[1].timeBetweenWriting = math.huge

-- Material Parameters

Problem.Material = {}
--Problem.Material.id = 'Alfeweld'
Problem.Material.mu = 1.3e-9                    -- Pa*s = kg/(m*s) = 10^(-6) ton/(mm*s)
Problem.Material.gamma = 0                      -- 
Problem.Material.rho = 2.71e-9                  -- ton/mm^(3)

Problem.Material.alphaLin = 0.0002              -- K^(-1)
Problem.Material.DgammaDT = 0                   -- 
Problem.Material.Tinf = 300                     -- K^(-1)
Problem.Material.k = 200                        -- ton * mm * s^(-3) * K^(-1)
Problem.Material.cp = 900e6                     -- mm^(2) * s^(-2) * K^(-1)
Problem.Material.h = 0                          -- W/(K m^2) = kg / (K s^3) =  10^(-3) ton / (K s^3)

-- phase change
Problem.Material.Tr = 933.5                     -- K^(-1)
Problem.Material.Tl = 933.5 + 5                 -- K^(-1) 482
Problem.Material.Ts = 933.5 - 5                 -- K^(-1)
Problem.Material.Lm = 400e9                     -- 400 kJ/kg = 400 * 10^9 mm^2/s^2 

-- Boussinesq PC
Problem.Material.C = 1e-3                       -- 
Problem.Material.eps = 1e-3                     -- 
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
Problem.Material.DmuDT = 0
Problem.Material.DcpDT = 0
Problem.Material.DkDT = 0


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
Problem.Solver.MomContEq.residual = 'Ax_f'
Problem.Solver.MomContEq.minRes   = 1e-6
Problem.Solver.MomContEq.maxIter  = 10
Problem.Solver.MomContEq.bodyForce = {0, -9810}

-- Heat Equation

Problem.Solver.HeatEq = {}
Problem.Solver.HeatEq.residual = 'Ax_f'
Problem.Solver.HeatEq.nlAlgo = 'Picard'
Problem.Solver.HeatEq.sparseSolverLib = 'Eigen'

Problem.Solver.HeatEq.maxIter = 10
Problem.Solver.HeatEq.minRes = 1e-6
Problem.Solver.HeatEq.cgTolerance = 1e-8

-- Heat Momentum Continuity BC

Problem.IC = {}
Problem.Solver.HeatEq.BC = {}
Problem.Solver.MomContEq.BC = {}
Problem.Solver.HeatEq.BC['FSInterfaceTExt'] = true

function Problem.IC.initStates(x,y,z)
	return {0, 0, 0, 300}
end


-- function Problem.Mesh.computeHcharFromDistance(x,y,z,t,T)	
--     
--     -- element sizes
--     local hmaxS = 0.5   -- size at room temperature         (at T0 -- T1)
--     local hminS = 0.3   -- size at higher temperature       (at T2)
--     local hmin  = 0.1   -- size at melting temperature      (at Trn--Tr--Trp)
--     local hmaxF = 0.3   -- maximum size at high temperature (at Tmax)
--     
--     -- reference temperatures
--     local T1   = 850
--     local T2   = 900
--     local Tr   = 933.5
--     local Trn  = Tr  - 5  
--     local Trp  = Tr  + 5
--     local Tmax = Trp + 30
-- 
--     -- element sizes by segments
--     local h0 = hmaxS
--     local h1 = (h0    - hminS)/(T1 - T2 ) * (T - T2 ) + hminS
--     local h2 = (hminS - hmin )/(T2 - Trn) * (T - Trn) + hmin
--     local h3 = hmin
--     local h4 = (hmin - hmaxF)/(Trp - Tmax) * (T - Tmax) + hmaxF
--     local h5 = hmaxF
-- 
--     if (T < T1) then
--         return h0
--     elseif (T < T2) then 
--         return h1
--     elseif (T < Trn) then 
--         return h2
--     elseif (T < Trp) then 
--         return h3
--     elseif (T < Tmax) then 
--         return h4
--     else
--         return h5            
--     end
-- end


function Problem.Solver.MomContEq.BC.FSInterfaceV(x, y, z, t) 
	return 0, 0
end
