# %% Input Parameters

def getParam(path):

    # Metafor and PFEM solvers
    
    param = dict()
    param['inputS'] = 'input_meta'
    param['inputF'] = path+'/input_pfem.lua'

    # FSI parameters

    param['interp'] = 'MM_CNS'
    param['algo'] = 'IQN_MVJ'
    param['omega'] = 0.5
    param['maxIt'] = 25
    param['tol'] = 1e-6

    # Time Parameters

    param['dt'] = 1e-3
    param['dtWrite'] = 0.01
    param['tTot'] = 4

    return param