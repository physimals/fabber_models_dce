/*  fwdmodel_dce_Patlak.cc - Implements a convolution based model for DCE analysis

    Jesper Kallehauge, IBME

    Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_Patlak.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "fabbercore/easylog.h"
#include "miscmaths/miscprob.h"

FactoryRegistration<FwdModelFactory, DCE_Patlak_FwdModel>
  DCE_Patlak_FwdModel::registration("dce_Patlak");

string DCE_Patlak_FwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_dce_Patlak.cc,v 1.11 2016/01/06 15:20:47 Kallehauge Exp $";
}

void DCE_Patlak_FwdModel::HardcodedInitialDists(MVNDist& prior,
    MVNDist& posterior) const
{
    Tracer_Plus tr("DCE_Patlak_FwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());

     SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Fp or Ktrans whatever you belive
     prior.means(Ktrans_index()) = 0.01;
     precisions(Ktrans_index(),Ktrans_index()) = 1e-12;

     prior.means(Vp_index()) = 0.01;
     precisions(Vp_index(),Vp_index()) = 1e-12;

     if (Acq_tech != "none") {
     precisions(sig0_index(),sig0_index())=1e-12;
     precisions(T10_index(),T10_index())=10;
     }


     if (inferdelay) {
       // delay parameter
       prior.means(delta_index()) = 0;
       precisions(delta_index(),delta_index()) = 0.04; //[0.1]; //<1>;
     }

   
    // Set precsions on priors
    prior.SetPrecisions(precisions);
    
    // Set initial posterior
    posterior = prior;

    // For parameters with uniformative prior chosoe more sensible inital posterior
    // Tissue perfusion
    posterior.means(Ktrans_index()) = 0.1;
    precisions(Ktrans_index(),Ktrans_index()) = 0.1;

    posterior.means(Vp_index()) = 0.1;
    precisions(Vp_index(),Vp_index()) = 0.1;

    posterior.SetPrecisions(precisions);
    
}    
    
    

void DCE_Patlak_FwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
  Tracer_Plus tr("DCE_Patlak_FwdModel::Evaluate");

    // ensure that values are reasonable
    // negative check
   ColumnVector paramcpy = params;
    for (int i=1;i<=NumParams();i++) {
      if (params(i)<0) { paramcpy(i) = 0; }
      }
  
   // parameters that are inferred - extract and give sensible names
   float Ktrans;
   float Vp; //mean of the transit time distribution
   float delta;
   float sig0; //'inital' value of the signal
   float T10;
   float FA_radians;

   // extract values from params
   Ktrans = paramcpy(Ktrans_index());
   Vp = paramcpy(Vp_index());
   if (Vp<1e-8) Vp=1e-8;
   if (Vp>1) Vp=1;
   if (Ktrans<1e-8) Ktrans=1e-8;

   if (inferdelay) {
   delta = params(delta_index()); // NOTE: delta is allowed to be negative
   }
   else {
     delta = 0;
   }
   if (Acq_tech != "none") {
   sig0 = paramcpy(sig0_index());
   T10 = paramcpy(T10_index());
   FA_radians=FA*3.1415926/180;
   }

     

   ColumnVector artsighere; // the arterial signal to use for the analysis
   if (artsig.Nrows()>0) {
     artsighere = artsig; //use the artsig that was loaded by the model
   }
   else {
     //use an artsig from supplementary data
     if (suppdata.Nrows()>0) {
       artsighere = suppdata;
     }
     else {
       cout << "No valid AIF found" << endl;
       throw;
     }
   }
   // use length of the aif to determine the number of time points
   int ntpts = artsighere.Nrows();



   // sensible limits on delta (beyond which it gets silly trying to estimate it)
   if (delta > ntpts/2*delt) {delta = ntpts/2*delt;}
   if (delta < -ntpts/2*delt) {delta = -ntpts/2*delt;}   



  //upsampled timeseries
  int upsample;
  int nhtpts;
  float hdelt;
  ColumnVector htsamp;


      // Create vector of sampled times
      ColumnVector tsamp(ntpts);
      for (int i=1; i<=ntpts; i++) {
	tsamp(i) = (i-1)*delt;
      }

      upsample=1;
      nhtpts=(ntpts-1)*upsample+1;
      htsamp.ReSize(nhtpts); 
      htsamp(1) = tsamp(1);
      hdelt = delt/upsample;
      for (int i=2; i<=nhtpts-1; i++) {
	htsamp(i) = htsamp(i-1)+hdelt;
      }
      htsamp(nhtpts)=tsamp(ntpts);

   // calculate the arterial input function (from upsampled artsig)
   ColumnVector aif_low(ntpts);
   aif_low=artsighere;
   
   // upsample the signal
   ColumnVector aif; 
   aif.ReSize(nhtpts);
   aif(1) = aif_low(1);
   int j=0; int ii=0;
   for (int i=2; i<=nhtpts-1; i++) {
     j = floor((i-1)/upsample)+1;
     ii= i - upsample*(j-1)-1;
     aif(i) = aif_low(j) + ii/upsample*(aif_low(j+1) - aif_low(j));
   }
   aif(nhtpts) = aif_low(ntpts);     

   // create the AIF matrix - empty for the time being
   LowerTriangularMatrix A(nhtpts); A=0.0;

   // deal with delay parameter - this shifts the aif
   ColumnVector aifnew(aif);
   aifnew = aifshift(aif,delta,hdelt);

     // populate AIF matrix
     createconvmtx(A,aifnew);

   ColumnVector Ktrans_vec(aifnew);
   for (int i=1; i<=ntpts; i++){
       Ktrans_vec(i)=Ktrans;
   }

   ColumnVector C;

   C = Vp*aifnew+hdelt*A*Ktrans_vec;

   ColumnVector C_low(ntpts);
   for (int i=1; i<=ntpts; i++) {
     C_low(i) = C((i-1)*upsample+1);
     //C_low(i) = interp1(htsamp,C,tsamp(i));
//     if (inferart && !artoption) { //add in arterial contribution
//       C_low(i) += C_art((i-1)*upsample+1);
//     }
     } 

    ColumnVector S_low(ntpts);
   if (Acq_tech == "SPGR") {
            for (int i=1; i<=ntpts; i++){
            S_low(i)=sig0*(1-exp(-TR*(1/T10+r1*C_low(i))))/(1-cos(FA_radians)*exp(-TR*(1/T10+r1*C_low(i))));//SPGR
            }
       }
   if (Acq_tech == "SRTF") {
           S_low=sig0*(1-exp(-Tsat*(1/T10+r1*C_low)))/(1-exp(-Tsat/T10));
           // for (int i=1; i<=ntpts; i++){
           // S_low(i)=sig0*(1-exp(-Tsat*(1/T10+r1*C_low(i))))*(1-exp(-TR*(1/T10+r1*C_low(i))))/(1-cos(FA_radians)*exp(-TR*(1/T10+r1*C_low(i))));//SRTF
           //      }
            }
    if (Acq_tech == "none") {
   S_low=C_low;
        }


    result.ReSize(ntpts);
    //cout<<C.t()<<endl;

    result=S_low;
//   ColumnVector sig_art(ntpts);
//   result.ReSize(ntpts);
//   for (int i=1; i<=ntpts; i++) {
     

//     if (inferart && artoption) {
//       sig_art(i) = C_art((i-1)*upsample+1);
//       sig_art(i) = exp(-sig_art(i)*te);

//       /*
//       float cbv = gmu*cbf;
//       float sumbv = artmag+cbv;
//       if (sumbv<1e-12) sumbv=1e-12; //catch cases where both volumes are zero
//       float ratio = artmag/sumbv;
//       result(i) = sig0*(1 + ratio*(sig_art(i)-1) + (1-ratio)*(exp(-C_low(i)*te)-1) ); //assume relative scaling is based on the relative proportions of blood volume
//       */
//       result(i) = sig0*(1 + (sig_art(i)-1) + (exp(-C_low(i)*te)-1) );
//     }
//     else {
//       result(i) = sig0*exp(-C_low(i)*te);
//     }
//   }

   for (int i=1; i<=ntpts; i++) {
     if (isnan(result(i)) || isinf(result(i))) {
       exit(0);
       LOG << "Warning NaN of inf in result" << endl;
       LOG << "result: " << result.t() << endl;
       LOG << "params: " << params.t() << endl;

       result=0.0; 
       break;
	 }
   }


   // downsample back to normal time points
   //cout << estsig.t() << endl;
   //result.ReSize(ntpts);
   //result=estsig;
   /*for (int i=1; i<=ntpts; i++) {
     result(i) = interp1(htsamp,estsig,tsamp(i));
     }
   if ((result-estsig).SumAbsoluteValue()>0.1){
     cout << result.t() << endl;
     cout << estsig.t() << endl;
     }*/

   //cout << result.t()<< endl;
}

FwdModel* DCE_Patlak_FwdModel::NewInstance()
{
  return new DCE_Patlak_FwdModel();
}

void DCE_Patlak_FwdModel::Initialize(ArgsType& args)
{
  Tracer_Plus tr("DCE_Patlak_FwdModel::DCE_Patlak_FwdModel");
    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    
    if (scanParams == "cmdline")
    {
      // specify command line parameters here
      delt = convertTo<double>(args.Read("delt"));


      // specify options of the model
      inferdelay = args.ReadBool("inferdelay");

      convmtx = args.ReadWithDefault("convmtx","simple");
      
      // Read in the arterial signal (this will override an image supplied as supplementary data)
      //ColumnVector artsig;
      string artfile = args.Read("aif");
      if (artfile != "none") {
	artsig = read_ascii_matrix( artfile );

    //cout<<data<<"  \n";
    //exit(0);
      }

       Acq_tech = args.ReadWithDefault("Acq_tech","none");
       cout<<Acq_tech<<endl;
       if (Acq_tech != "none") {
           if (Acq_tech == "SPGR") {
                    FA = convertTo<double>(args.Read("FA"));
                    TR = convertTo<double>(args.Read("TR"));
                    r1 = convertTo<double>(args.Read("r1"));
                  }
                    if (Acq_tech == "SRTF") {
                             FA = convertTo<double>(args.Read("FA"));
                             TR = convertTo<double>(args.Read("TR"));
                             r1 = convertTo<double>(args.Read("r1"));
                             Tsat    = convertTo<double>(args.Read("Tsat"));
                           }

       }

      aifconc = args.ReadBool("aifconc"); // indicates that the AIF is a CTC not signal curve
        cout<<aifconc<<"  \n";
      doard=false;
     // if (inferart) doard=true;

      imageprior = args.ReadBool("imageprior"); //temp way to indicate we have some image priors (very fixed meaning!)             

      // add information about the parameters to the log
      /* do logging here*/
   
    }

    else
        throw invalid_argument("Only --scan-params=cmdline is accepted at the moment");    
    
 
}

vector<string> DCE_Patlak_FwdModel::GetUsage() const
{
  vector<string> usage;

  usage.push_back( "\nThis model is a one compartment model\n");
  usage.push_back( "It returns  2 parameters :\n");
  usage.push_back( " Ktrans: the Forward transport constant\n");
  usage.push_back( " Vp: the plasma volume fraction\n");

  return usage;
}

void DCE_Patlak_FwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    
}

void DCE_Patlak_FwdModel::NameParams(vector<string>& names) const
{
  names.clear();
  
  names.push_back("Ktrans");
  names.push_back("Vp");
  if (inferdelay)
  names.push_back("delay");
   if (Acq_tech != "none") {
    names.push_back("T10");
    names.push_back("sig0");
   }



  //names.push_back("sig0");
  
}



ColumnVector DCE_Patlak_FwdModel::aifshift( const ColumnVector& aif, const float delta, const float hdelt ) const
{
  // Shift a vector in time by interpolation (linear)
  // NB Makes assumptions where extrapolation is called for.
   int nshift = floor(delta/hdelt); // number of time points of shift associated with delta
   float minorshift = delta - nshift*hdelt; // shift within the sampled time points (this is always a 'forward' shift)
      
   ColumnVector aifnew(aif);
   int index;
   int nhtpts = aif.Nrows();
   for (int i=1; i<=nhtpts; i++) {
     index = i-nshift;
     if (index==1) { aifnew(i) = aif(1)*minorshift/hdelt; } //linear interpolation with zero as 'previous' time point
     else if (index < 1) { aifnew(i) = 0; } // extrapolation before the first time point - assume aif is zero
     else if (index>nhtpts) { aifnew(i) = aif(nhtpts); } // extrapolation beyond the final time point - assume aif takes the value of the final time point
     else {
       //linear interpolation
       aifnew(i) = aif(index) + (aif(index-1)-aif(index))*minorshift/hdelt;
     }
   }
   return aifnew;
}

void DCE_Patlak_FwdModel::createconvmtx( LowerTriangularMatrix& A, const ColumnVector aifnew ) const
{
  // create the convolution matrix
  int nhtpts = aifnew.Nrows();

   if (convmtx=="simple")
     {
       // Simple convolution matrix
       for (int i=1; i<=nhtpts; i++) {//
     for (int j=1; j <= i; j++) {
	   A(i,j) = aifnew(i-j+1); //note we are using the local aifnew here! (i.e. it has been suitably time shifted)

     }
       }
    }

   //cout << A << endl;

   else if (convmtx=="voltera")
     {
       ColumnVector aifextend(nhtpts+2);
       ColumnVector zero(1);
       zero=0;
       aifextend = zero & aifnew & zero;
       int x, y, z;
       //voltera convolution matrix (as defined by Sourbron 2007) - assume zeros outside aif range
       for (int i=1; i<=nhtpts; i++) {
	 for (int j=1; j <= i; j++) {
	   //cout << i << "  " << j << endl;
	   x = i+1;y=j+1; z = i-j+1;
	   if (j==1) { A(i,j) =(2*aifextend(x) + aifextend(x-1))/6; }
	   else if (j==i) { A(i,j) = (2*aifextend(2) + aifextend(3))/6; }
	   else { 
	     A(i,j) =  (4*aifextend(z) + aifextend(z-1) + aifextend(z+1))/6; 
	     //cout << x << "  " << y << "  " << z << "  " << ( 4*aifextend(z) + aifextend(z-1) + aifextend(z+1) )/6 << "  " << 1/6*(4*aifextend(z) + aifextend(z-1) + aifextend(z+1)) << endl;
	     // cout << aifextend(z) << "  " << aifextend(z-1) << "  " << aifextend(z+1) << endl;
	   }
	   //cout << i << "  " << j << "  " << aifextend(z) << "  " << A(i,j) << endl<<endl;
	 }
       }
     }
}