/*  fwdmodel_dce_LLS.cc - Implements a convolution based model for DCE analysis

 Jesper Kallehauge, IBME

 Copyright (C) 2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_LLS.h"

#include <iostream>
#include <complex>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "fabber_core/easylog.h"
#include "miscmaths/miscprob.h"

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCE_LLS_FwdModel> DCE_LLS_FwdModel::registration("dce_LLS");

std::string DCE_LLS_FwdModel::GetDescription() const
{
	return "The Linear One Compartment model";
}

void DCE_LLS_FwdModel::HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const
{
	assert(prior.means.Nrows() == NumParams());

	SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

	// Set priors
	// Fp or Ktrans whatever you belive
	prior.means(fp_idx) = 0.01;
	precisions(fp_idx, fp_idx) = 1e-12;

	prior.means(vp_idx) = 0.01;
	precisions(vp_idx, vp_idx) = 1e-12;

	if (Acq_tech != "none")
	{
		precisions(sig0_idx, sig0_idx) = 1e-2;
		precisions(t10_idx, t10_idx) = 10;
	}

	if (inferdelay)
	{
		// delay parameter
		prior.means(delta_idx) = 0;
		precisions(delta_idx, delta_idx) = 0.04; //[0.1]; //<1>;
	}

	// Set precsions on priors
	prior.SetPrecisions(precisions);

	// Set initial posterior
	posterior = prior;

	// For parameters with uniformative prior chosoe more sensible inital posterior
	// Tissue perfusion
	posterior.means(fp_idx) = 0.1;
	precisions(fp_idx, fp_idx) = 0.1;

	posterior.means(vp_idx) = 0.1;
	precisions(vp_idx, vp_idx) = 0.1;

	posterior.SetPrecisions(precisions);

}

void DCE_LLS_FwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
	// ensure that values are reasonable
	// negative check
	ColumnVector paramcpy = params;
	for (int i = 1; i <= NumParams(); i++)
	{
		if (params(i) < 0)
		{
			paramcpy(i) = 0;
		}
	}

	// parameters that are inferred - extract and give sensible names
	float Fp;
	float Vp; //mean of the transit time distribution
	float delta;
	float sig0; //'inital' value of the signal
	float T10;
	float FA_radians;
	int baseline;

	ColumnVector Cvec(data);
	ColumnVector C(data);
	ColumnVector denom;
	ColumnVector numer;
	ColumnVector R1_nl;

	// extract values from params
	Fp = paramcpy(fp_idx);
	Vp = paramcpy(vp_idx);
	if (Vp < 1e-8)
		Vp = 1e-8;
	if (Vp > 1)
		Vp = 1;
	if (Fp < 1e-8)
		Fp = 1e-8;

	if (inferdelay)
	{
		delta = params(delta_idx); // NOTE: delta is allowed to be negative
	}
	else
	{
		delta = 0;
	}

	ColumnVector artsighere; // the arterial signal to use for the analysis
	if (artsig.Nrows() > 0)
	{
		artsighere = artsig; //use the artsig that was loaded by the model
	}
	else
	{
		//use an artsig from supplementary data
		if (suppdata.Nrows() > 0)
		{
			artsighere = suppdata;
		}
		else
		{
			cout << "No valid AIF found" << endl;
			throw;
		}
	}
	// use length of the aif to determine the number of time points
	int ntpts = artsighere.Nrows();

	// sensible limits on delta (beyond which it gets silly trying to estimate it)
	if (delta > ntpts / 2 * delt)
	{
		delta = ntpts / 2 * delt;
	}
	if (delta < -ntpts / 2 * delt)
	{
		delta = -ntpts / 2 * delt;
	}

	int aifpeak;
	float deltpeak;
	artsighere.MaximumAbsoluteValue1(aifpeak);
	deltpeak = (10.0 / 60.0) / delt; //10 seconds earlier than tmax is the end of the baseline
	baseline = aifpeak - ceil(deltpeak);

	if (Acq_tech != "none")
	{
		sig0 = paramcpy(sig0_idx);
		if (sig0 > 1)
		{
			sig0 = paramcpy(sig0_idx);
		}
		else
		{
			float mean_base;
			mean_base = 0.0;
			for (int i = 1; i <= baseline; i++)
			{
				mean_base = mean_base + data(i);
			}
			sig0 = mean_base / baseline;
		}
		T10 = paramcpy(t10_idx);
		FA_radians = FA * 3.1415926 / 180;
	}
	if (Acq_tech == "SPGR")
	{
		for (int i = 1; i <= ntpts; i++)
		{
			denom(i) = ((data(i) - sig0) / sig0) * (exp(-TR / T10) - 1.0) + exp(-TR / T10) * (1.0 - cos(FA_radians));
			numer(i) = 1.0 + cos(FA_radians) * (((data(i) - sig0) / sig0) * (exp(-TR / T10) - 1.0) - 1.0);
			R1_nl(i) = (-1.0 / TR) * (log(numer(i) / denom(i)));
			Cvec(i) = (1.0 / r1) * (R1_nl(i) - 1.0 / T10);
		}
	}
	if (Acq_tech == "SRTF")
	{
		for (int i = 1; i <= ntpts; i++)
		{
			float ctemp;
			ctemp = (1.0 - exp(Tsat / T10)) * (((data(i)) / sig0) + (exp(Tsat / T10) / (1 - exp(Tsat / T10))));
			//cout<<sig0<<endl;
			if (ctemp > 0)
			{
				//cout<<ctemp<<endl;
				Cvec(i) = ((-1.0
						* log(
								(1.0 - exp(Tsat / T10))
										* (((data(i)) / sig0) + (exp(Tsat / T10) / (1 - exp(Tsat / T10))))))
						/ (r1 * Tsat));
			}
			else
			{
				//cout<<"Cvec(i)"<<endl;
				Cvec(i) = 0.0;
			}

		}
	}
	if (Acq_tech == "none")
	{
		Cvec = data;
	}
	//cout<<Cvec<<endl;
	// exit(0);

	//upsampled timeseries
	int upsample;
	int nhtpts;
	float hdelt;
	ColumnVector htsamp;

	// Create vector of sampled times
	ColumnVector tsamp(ntpts);
	for (int i = 1; i <= ntpts; i++)
	{
		tsamp(i) = (i - 1) * delt;
	}

	upsample = 1;
	nhtpts = (ntpts - 1) * upsample + 1;
	htsamp.ReSize(nhtpts);
	htsamp(1) = tsamp(1);
	hdelt = delt / upsample;
	for (int i = 2; i <= nhtpts - 1; i++)
	{
		htsamp(i) = htsamp(i - 1) + hdelt;
	}
	htsamp(nhtpts) = tsamp(ntpts);

	// calculate the arterial input function (from upsampled artsig)
	ColumnVector aif_low(ntpts);
	aif_low = artsighere;

	// upsample the signal
	ColumnVector aif;
	aif.ReSize(nhtpts);
	aif(1) = aif_low(1);
	int j = 0;
	int ii = 0;
	for (int i = 2; i <= nhtpts - 1; i++)
	{
		j = floor((i - 1) / upsample) + 1;
		ii = i - upsample * (j - 1) - 1;
		aif(i) = aif_low(j) + ii / upsample * (aif_low(j + 1) - aif_low(j));
	}
	aif(nhtpts) = aif_low(ntpts);

	// create the AIF matrix - empty for the time being
	Matrix A(nhtpts, 2);
	A = 0.0;

	// deal with delay parameter - this shifts the aif
	ColumnVector aifnew(aif);
	aifnew = aifshift(aif, delta, hdelt);

	//----fill the matrix for linear inversion-----

	float temp;
	float temp2;
	for (int i = baseline; i <= ntpts; i++)
	{
		temp = 0.0;
		temp2 = 0.0;
		for (int j = baseline; j <= i; j++)
		{
			temp = temp + aifnew(j);
			temp2 = temp2 + Cvec(j);
		}
		A(i, 2) = -temp2;
		A(i, 1) = temp;
	}

	ColumnVector B(2);
	B(1) = Fp * delt;
	B(2) = B(1) / Vp;
	C = A * B;

	//exit(0);

	// cout<<A(1,2)<<endl;
	// C = hdelt*A*B;
	//convert to the DCE signal

	//cout  << "C: " << C.t() << endl;
	// exit(0);
	//cout << "sig0: " << sig0 << " r2: " << r2 << " te: " << te << endl;

	//cout<< htsamp.t() << endl;

	ColumnVector C_low(ntpts);
	for (int i = 1; i <= ntpts; i++)
	{
		C_low(i) = C((i - 1) * upsample + 1);
		//C_low(i) = interp1(htsamp,C,tsamp(i));
//     if (inferart && !artoption) { //add in arterial contribution
//       C_low(i) += C_art((i-1)*upsample+1);
//     }
	}

	ColumnVector S_low(ntpts);
	if (Acq_tech == "SPGR")
	{
		for (int i = 1; i <= ntpts; i++)
		{
			S_low(i) = sig0 * (1 - exp(-TR * (1 / T10 + r1 * C_low(i))))
					/ (1 - cos(FA_radians) * exp(-TR * (1 / T10 + r1 * C_low(i))));   //SPGR
		}
	}
	if (Acq_tech == "SRTF")
	{
		S_low = sig0 * (1 - exp(-Tsat * (1 / T10 + r1 * C_low))) / (1 - exp(-Tsat / T10));
		// for (int i=1; i<=ntpts; i++){
		// S_low(i)=sig0*(1-exp(-Tsat*(1/T10+r1*C_low(i))))*(1-exp(-TR*(1/T10+r1*C_low(i))))/(1-cos(FA_radians)*exp(-TR*(1/T10+r1*C_low(i))));//SRTF
		//      }
	}
	if (Acq_tech == "none")
	{
		S_low = C_low;
	}

	result.ReSize(ntpts);
	//cout<<C.t()<<endl;

	result = S_low;
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

	for (int i = 1; i <= ntpts; i++)
	{
		if (isnan(result(i)) || isinf(result(i)))
		{
			//exit(0);
			LOG << "Warning NaN of inf in result" << endl;
			LOG << "result: " << result.t() << endl;
			LOG << "params: " << params.t() << endl;

			result = 0.0;
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

FwdModel* DCE_LLS_FwdModel::NewInstance()
{
	return new DCE_LLS_FwdModel();
}

vector<string> DCE_LLS_FwdModel::GetUsage() const
{
	vector<string> usage;

	usage.push_back("\nThis model is a one compartment model\n");
	usage.push_back("It returns  2 parameters :\n");
	usage.push_back(" Fp: the Plasma flow constant\n");
	usage.push_back(" Vp: the plasma volume fraction\n");

	return usage;
}

void DCE_LLS_FwdModel::NameParams(vector<string>& names) const
{
	names.clear();

	names.push_back("Fp");
	names.push_back("Vp");
	if (inferdelay)
		names.push_back("delay");
	if (Acq_tech != "none")
	{
		names.push_back("T10");
		names.push_back("sig0");
	}

	//names.push_back("sig0");

}

