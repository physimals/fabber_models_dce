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

Mandatory:

    --output=/home/fsl/Desktop/Data_out/Data
    --data=/home/fsl/Desktop/Data/DCE_signal_intensity.nii
    --mask=/home/fsl/Desktop/Data/mask.nii
    --aif=/home/fsl/Desktop/Data/parker.dat  //in mM for MR data
    --aifconc  // this just makes sure you consider the above mM concentration of the AIF 
    --method=vb       //[nlls, vb, spatialvb]   NB! spatialvb is recommended!
    --model=dce       // See above
    --delt=0.035      // in Minutes
    --noise=white
    --data-order=singlefile
    
Optional:    

    --convmtx=expConv //[expConv, simple,  voltera]
    --inferdelay
    --save-model-fit
    --Acq_tech=SRTF   //[none, SPGR, SRTF] 
    --Tsat=0.025      // in Seconds (only for SRTF) 
    --FA=10           // in degrees (For both SRTF and SPGR)   
    --TR=0.0029       // in Seconds (For both SRTF and SPGR)
    --r1=3.6          // s^(-1) mM^(-1) (For both SRTF and SPGR)
    --PSP_byname1=T10 // Spatial prior variable
    --PSP_byname1_type=I // Spatial prior data type
    --PSP_byname1_image=/home/fsl/Desktop/Data/T1_map.nii // Spatial prior data location
    --PSP_byname2=sig0  // Spatial prior variable
    --PSP_byname2_type=I // Spatial prior data type
    --PSP_byname2_image=/home/fsl/Desktop/Data/DCE_signal_baseline.nii  // Spatial prior data location
    --mcsteps=2     //number of motion correction step