/*  fwdmodel_dce_ETM.cc - Implements a convolution based model for DCE analysis

 Jesper Kallehauge, IBME

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_dce_ETM.h"

#include "newimage/newimageall.h"
#include <iostream>
#include "armawrap/newmat.h"
#include <stdexcept>
using namespace NEWIMAGE;
#include "fabber_core/easylog.h"
#include "miscmaths/miscprob.h"

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, DCE_ETM_FwdModel> DCE_ETM_FwdModel::registration("dce_ETM");

std::string DCE_ETM_FwdModel::GetDescription() const
{
    return "The Extended Tofts model";
}

void DCE_ETM_FwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());

    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // Set priors
    // Ktrans or Ktrans whatever you belive
    prior.means(ktrans_idx) = 0.01;
    precisions(ktrans_idx, ktrans_idx) = 1e-12;

    prior.means(ve_idx) = 0.01;
    precisions(ve_idx, ve_idx) = 1e-12;

    prior.means(vp_idx) = 0.01;
    precisions(vp_idx, vp_idx) = 1e-12;

    if (Acq_tech != "none")
    {
        precisions(sig0_idx, sig0_idx) = 1e-12;
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

    // For parameters with uninformative prior chosoe more sensible inital posterior
    // Tissue perfusion
    posterior.means(ktrans_idx) = 0.1;
    precisions(ktrans_idx, ktrans_idx) = 0.1;

    posterior.means(vp_idx) = 0.1;
    precisions(vp_idx, vp_idx) = 0.1;

    posterior.means(ve_idx) = 0.1;
    precisions(ve_idx, ve_idx) = 0.1;

    posterior.SetPrecisions(precisions);
}

void DCE_ETM_FwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
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
    float Ktrans;
    float Vp; //mean of the transit time distribution
    float Ve;
    float delta;
    float sig0; //'inital' value of the signal
    float T10;
    float FA_radians;

    // extract values from params
    Ktrans = paramcpy(ktrans_idx);
    Vp = paramcpy(vp_idx);
    Ve = paramcpy(ve_idx);
    //cout<<"Ktrans = "<< Ktrans<<endl;
    //cout<<"Vp = "<< Vp<<endl;
    //cbf = exp(params(cbf_index()));
    if (Vp < 1e-8)
        Vp = 1e-8;
    if (Vp > 1)
        Vp = 1;
    if (Ve < 1e-8)
        Ve = 1e-8;
    if (Ve > 1)
        Ve = 1;

    if (Ktrans < 1e-8)
        Ktrans = 1e-8;

    if (inferdelay)
    {
        delta = params(delta_idx); // NOTE: delta is allowed to be negative
    }
    else
    {
        delta = 0;
    }
    if (Acq_tech != "none")
    {
        if (Acq_tech == "CT")
        {
            sig0 = paramcpy(sig0_idx);
        }
        else
        {
            sig0 = paramcpy(sig0_idx);
            T10 = paramcpy(t10_idx);
            FA_radians = FA * 3.1415926 / 180;
        }
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
    LowerTriangularMatrix A(nhtpts);
    A = 0.0;

    // deal with delay parameter - this shifts the aif
    ColumnVector aifnew(aif);
    aifnew = aifshift(aif, delta, hdelt);

    // populate AIF matrix
    createconvmtx(A, aifnew);

    // --- Residue Function ----
    ColumnVector residue(nhtpts);
    residue = 0.0;
    residue = exp(-(Ktrans / Ve) * htsamp); // the 1 compartment model.

    // do the multiplication
    ColumnVector C;

    if (convmtx == "expConv")
    {
        float T;
        T = Ve / Ktrans;
        C = Ktrans * T * expConv(aifnew, T, htsamp) + Vp * aifnew;
    }
    else
    {
        C = Ktrans * hdelt * A * residue + Vp * aifnew;
    }

    //convert to the DCE signal
    ColumnVector C_low(ntpts);
    for (int i = 1; i <= ntpts; i++)
    {
        C_low(i) = C((i - 1) * upsample + 1);
    }

    ColumnVector S_low(ntpts);
    if (Acq_tech == "SPGR")
    {
        for (int i = 1; i <= ntpts; i++)
        {
            S_low(i) = sig0 * (1 - exp(-TR * (1 / T10 + r1 * C_low(i))))
                / (1 - cos(FA_radians) * exp(-TR * (1 / T10 + r1 * C_low(i)))); //SPGR
        }
    }
    if (Acq_tech == "SRTF")
    {
        S_low = sig0 * (1 - exp(-Tsat * (1 / T10 + r1 * C_low))) / (1 - exp(-Tsat / T10));
    }
    if (Acq_tech == "CT")
    {
        S_low = C_low + sig0;
    }
    if (Acq_tech == "none")
    {
        S_low = C_low;
    }

    result.ReSize(ntpts);
    result = S_low;

    for (int i = 1; i <= ntpts; i++)
    {
        if (isnan(result(i)) || isinf(result(i)))
        {
            LOG << "Warning NaN of inf in result" << endl;
            LOG << "result: " << result.t() << endl;
            LOG << "params: " << params.t() << endl;

            result = 0.0;
            break;
        }
    }
}

FwdModel *DCE_ETM_FwdModel::NewInstance()
{
    return new DCE_ETM_FwdModel();
}

vector<string> DCE_ETM_FwdModel::GetUsage() const
{
    vector<string> usage;

    usage.push_back("\nThis is the extended Tofts model (JMRI 1999;10:223-232)\n");
    usage.push_back("It returns  3 parameters :\n");
    usage.push_back(" Ktrans: the forward transfer constant\n");
    usage.push_back(" Ve: the extravascular extracellular volume fraction\n");
    usage.push_back(" Vp: the plasma volume fraction\n");

    return usage;
}

void DCE_ETM_FwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    names.push_back("Ktrans");
    names.push_back("Ve");
    names.push_back("Vp");
    if (inferdelay)
        names.push_back("delay");
    if (Acq_tech != "none")
    {
        if (Acq_tech == "CT")
        {
            names.push_back("sig0");
        }
        else
        {
            names.push_back("sig0");
            names.push_back("T10");
        }
    }
}
