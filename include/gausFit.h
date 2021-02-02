#ifndef gausFit_h
#define gausFit_h

#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include <iostream>

TF1* cleverGaus(TH1* h, const char* title="h", Float_t c = 1.5, bool quietMode=true)
{
    if ( h->GetEntries() == 0 )
    {
        TF1 *fit0 = new TF1(title,"gaus",-1,1);
        fit0->SetParameters(0,0,0);
        return fit0;
    }

    Int_t peakBin  = h->GetMaximumBin();
    Double_t peak =  h->GetBinCenter(peakBin);
    Double_t sigma = h->GetRMS();

    TF1 *fit1 = new TF1(title,"gaus",peak-c*sigma,peak+c*sigma);
    fit1->SetParameter(1, peak);
    fit1->SetParameter(2, sigma);
    //fit1->SetParameter(1, 1.0);
    //fit1->SetParameter(1, 0.0005);
    if (quietMode) h->Fit(fit1,"LL M O Q R");
    else    h->Fit(fit1,"LL M O Q R");
    return fit1;
}

#endif
