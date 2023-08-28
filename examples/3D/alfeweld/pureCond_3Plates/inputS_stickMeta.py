# -*- coding: utf-8 -*-
#
# test of
from wrap import *

def getMetafor(p={}):
    d={}
    
    p['metaforStandalone'] = True
    
    d.update(p)
    import inputS_stick as m
    return m.getMetafor(d)
