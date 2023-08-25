# -*- coding: utf-8 -*-
#
# test of
from wrap import *

def getMetafor(p={}):
    d={}
    # interaction mode for interface heat transfert
    p['heatMode'] = "contact"	#contact/sticking
    
    d.update(p)
    import inputS as m
    return m.getMetafor(d)
