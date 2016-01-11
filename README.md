This folder includes popular Compartmental models for analysis of Dynamic contrast enhanced MRI for in the FABBER framework
Currently four models are included:

    1. The One compartment model often refered to as the standard Tofts model 
        model-name= "dce"
    2. The Extended Tofts model
        model-name= "dce_ETM"
    3. The Compartmental Tissue Uptake model
        model-name= "dce_CTU"
    4. The 2 Compartment Exchange model.
        model-name= "dce_2CXM"
    
These models all need an estimate of an input function which must be supplied in the when calling FABBER in the command line
