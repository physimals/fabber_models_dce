#!/bin/env python

import os, sys
import traceback
import math

FSLDIR = os.environ["FSLDIR"]
sys.path.insert(0, FSLDIR + "/lib/python")
from fabber import self_test, FabberException

LOG = False

def proc(v):
    if LOG: return math.log(v)
    else: return v

TEST_DATA = {
    "dce_tofts" : [
        {"max-iterations" : "30", "aif" : "orton", "delt" : 0.1, "fa" : 12, "tr" : 0.0045, "r1" : 4.5}, # Model options
        {"Ktrans" : [proc(v) for v in [0.1, 0.2, 0.3, 0.4, 0.5]],
         "Ve" : [proc(v) for v in [0.1, 0.2, 0.3, 0.4, 0.5]]}, # Parameter values - at most 3 can vary
        {"nt" : 25} # Other options
    ]
}

save = "--save" in sys.argv

try:
    for model, test_data in TEST_DATA.items():
        rundata, params, kwargs = test_data
        if LOG: rundata["use-log-params"] = ""
        log = self_test(model, rundata, params, save_input=save, save_output=save, invert=True, noise=0.01, **kwargs)
        #print(log)
except FabberException, e:
    #print e.log
    traceback.print_exc()
except:
    traceback.print_exc()

