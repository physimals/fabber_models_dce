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
    5. The Adiabatic Approximation to the Tissue Homogeniety model.
        model-name= "dce_AATH"
    6. The Linear One Compartment model.
        model-name= "dce_LLS"
    7. The Linear Compartmental Tissue Uptake model.
        model-name= "dce_CTU_LLS"

The linear model only work when it is supplied with Concentration curves rather then raw data.
    
These models all need an estimate of an input function which must be supplied in the when calling FABBER in the command line

example of command line execution

--output=/home/fsl/Desktop/Data_out/Data
--data=/home/fsl/Desktop/Data/DCE_signal_intensity.nii
--mask=/home/fsl/Desktop/Data/mask.nii
--aifconc
--aif=/home/fsl/Desktop/Data/parker.dat
--method=vb
--model=dce_ETM
--inferdelay
--delt=0.035
--noise=white
--data-order=singlefile
--save-model-fit
--Acq_tech=SRTF
--Tsat=0.025
--FA=10
--TR=0.0029
--r1=3.6
--PSP_byname1=T10
--PSP_byname1_type=I
--PSP_byname1_image=/home/fsl/Desktop/Data/T1_map.nii
--PSP_byname2=sig0
--PSP_byname2_type=I
--PSP_byname2_image=/home/fsl/Desktop/Data/DCE_signal_baseline.nii
--mcsteps=2