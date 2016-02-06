This folder includes popular Compartmental models for analysis of Dynamic contrast enhanced MRI for in the FABBER framework
Currently four models are included:

    1. The One compartment model often refered to as the standard Tofts model 
        model-name= "dce"
    2. The Patlak model
        model-name= "dce_Patlak"
    3. The Extended Tofts model
        model-name= "dce_ETM"
    4. The Compartmental Tissue Uptake model
        model-name= "dce_CTU"
    5. The 2 Compartment Exchange model.
        model-name= "dce_2CXM"
    6. The Adiabatic Approximation to the Tissue Homogeniety model.
        model-name= "dce_AATH"
    7. The Linear One Compartment model.
        model-name= "dce_LLS"
    8. The Linear Extended Tofts model.
        model-name= "dce_ETM_LLS"
    9. The Linear Compartmental Tissue Uptake model.
        model-name= "dce_CTU_LLS"
    10. The Linear 2 Compartmental Exchange model.
        model-name= "dce_2CXM_LLS"

The linear model only work when it is supplied with Concentration curves rather then raw data. And It does not infer any delay between the arrival of the AIF and the tissue curves. A word of caution is that the 2CXM requires both high temporal sampling and good CNR!!
    
These models all need an estimate of an input function which must be supplied in the when calling FABBER in the command line

example of command line execution

    --output=/home/fsl/Desktop/Data_out/Data
    --data=/home/fsl/Desktop/Data/DCE_signal_intensity.nii
    --mask=/home/fsl/Desktop/Data/mask.nii
    --aifconc
    --aif=/home/fsl/Desktop/Data/parker.dat
    --method=vb       //[nlls, vb, spatialvb] 
    --model=dce       // See above
    --convmtx=expConv //[expConv, simple,  voltera]
    --inferdelay
    --delt=0.035      // in Minutes
    --noise=white
    --data-order=singlefile
    --save-model-fit
    --Acq_tech=SRTF
    --Tsat=0.025      // in Seconds
    --FA=10           // in degrees    
    --TR=0.0029       // in Seconds
    --r1=3.6          // s^(-1) mM^(-1)
    --PSP_byname1=T10
    --PSP_byname1_type=I
    --PSP_byname1_image=/home/fsl/Desktop/Data/T1_map.nii
    --PSP_byname2=sig0
    --PSP_byname2_type=I
    --PSP_byname2_image=/home/fsl/Desktop/Data/DCE_signal_baseline.nii
    --mcsteps=2