#!/bin/env python

import os, sys
import traceback
import math

from fabber import self_test, FabberException

TEST_DATA = {
    "dce_tofts" : [
        {
            "convergence" : "trialmode",
            "max-trials" : 20,
            "max-iterations" : "30", 
            "aif" : "orton", 
            "delt" : 0.1, 
            "fa" : 12, 
            "tr" : 0.0045, 
            "r1" : 4.5,
            "delay" : 0.5,
            "sig0" : 10000,
        }, # Model options
        {
            "ktrans" : [0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8],
            "ve" : [0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8],
            #"vp" : 0.5,
        }, 
        {
            "nt" : 25,
            "patchsize" : 2,
            "noise" : 10,
        } # Other options
    ]
}

save = "--save" in sys.argv

try:
    for model, test_data in TEST_DATA.items():
        rundata, params, kwargs = test_data
        log = self_test(model, rundata, params, save_input=save, save_output=save, 
                        invert=True, **kwargs)
        #print(log)
except FabberException, e:
    print e.log
    traceback.print_exc()
except:
    traceback.print_exc()

