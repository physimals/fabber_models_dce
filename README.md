# FABBER_MODELS_DCE

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

    --convmtx=expConv //[expConv, simple,  voltera] (expConv works only for models: [dce, dce_ETM, dce_CTU, dce_2CXM])
    --inferdelay
    --save-model-fit
    --Acq_tech=SRTF   //[none, SPGR, SRTF, CT]  //CT is only implemented for [dce, dce_ETM, dce_CTU, dce_2CXM, dce_AATH]
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
    
# Building the DCE models

The current version of the DCE models are designed to build with the latest 
version of Fabber. Although it should be possible to build them against FSL 5.0 
you may need to edit the source slightly. So your first step should be to build 
[fabber_core](https://ibme-gitcvs.eng.ox.ac.uk/fabber/fabber_core) from Git. 
Build instructions for fabber_core are in the Wiki on its Git page. 

You will make the process straightforward if you create a directory (e.g. 
'fabber') and download the Git repositories for fabber_core and fabber_models_dce
into this directory. Then the model build should be able to find the Fabber 
libraries without you having to install them or tell the build scripts where 
they are.

### Building in a standard FSL environment

The process should be:

    source $FSLDIR/etc/fslconf/fsl-devel.sh
    make

This will produce an executable `fabber_dce`.

### Building with non-standard FSL

In this case we can use the cross-platform `cmake` system. Convenience scripts
are provided, so the build should be:

    scripts/build.sh relwithdebinfo
    
This creates a release build (fully optimized) but including debugging symbols
so crashes can be traced. The output is in the `build_relwithdebinfo` directory
and will include the executable `fabber_dce` as well as a shared library
(e.g. `libfabber_models_dce.so` on Linux)

#### It failed with something about `recompile with -fPIC`

Short answer: if you were using the `cmake` build method, try the `standard
FSL environment` method.

You can also go into the build directory and just build the executable:

    cd build_relwithdebinfo
    make fabber_dce
    
Long answer: In order to build the shared library, the FSL libraries have to 
contain 'Position-independent code'. If they don't you can't build the library 
and will need to use the fabber_dce executable instead. The only way around 
this is to recompile the FSL libraries with the -fPIC flag which you can
do using the CMake-enabled code [here](https://ibme-gitcvs.eng.ox.ac.uk/fsl/fsl).
It's not worth doing this unless you really need the shared library.

## Successful build?

You can verify that the DCE models are available whichever method you choose:

    fabber_dce --listmodels
    
    dce
    dce_2CXM
    dce_2CXM_LLS
    dce_AATH
    dce_CTU
    dce_CTU_LLS
    dce_ETM
    dce_ETM_LLS
    dce_LLS
    dce_Patlak
    linear
    poly

or if you used cmake and built the shared library, you can also do:

    fabber --loadmodels=libfabber_models_dce.so --listmodels
    
    dce
    dce_2CXM
    dce_2CXM_LLS
    dce_AATH
    dce_CTU
    dce_CTU_LLS
    dce_ETM
    dce_ETM_LLS
    dce_LLS
    dce_Patlak
    linear
    poly
    