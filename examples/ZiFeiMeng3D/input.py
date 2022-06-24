# %% Input Parameters

def getParam(path):

    # Metafor and PFEM solvers
    
    param = dict()
    param['inputS'] = 'input_meta'
    param['inputF'] = path+'/input_pfem.lua'
    
    # Algorithm parameters

    param['algo'] = 'IQN_ILS'
    param['retainStep'] = 0
    param['keepStep'] = 0
    param['omega'] = 0.5
    param['maxIt'] = 25
    param['tol'] = 1e-6

    # Time Parameters

    param['dt'] = 1e-3
    param['dtWrite'] = 1e-3
    param['tTot'] = 0.4

    return param