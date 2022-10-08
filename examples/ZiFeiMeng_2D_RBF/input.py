# %% Input Parameters

def getParam(path):

    # Metafor and PFEM solvers
    
    param = dict()
    param['inputS'] = 'input_meta'
    param['inputF'] = path+'/input_pfem.lua'
    
    # Algorithm parameters

    param['RBF'] = 'IMQ'
    param['radius'] = 0.01
    param['interp'] = 'RBF'
    param['algo'] = 'IQN_MVJ'
    param['omega'] = 0.5
    param['maxIt'] = 25
    param['tol'] = 1e-8

    # Time Parameters

    param['dt'] = 1e-3
    param['dtWrite'] = 1e-3
    param['tTot'] = 0.4

    return param