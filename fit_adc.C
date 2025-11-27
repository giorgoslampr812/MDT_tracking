#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>

// ---- Skew-normal function ----
// params: [0]=Norm, [1]=Mean, [2]=Sigma, [3]=Alpha(skew)
double skewnormal(double *x, double *par) {
    double norm  = par[0];
    double mean  = par[1];
    double sigma = par[2];
    double alpha = par[3];

    double t = (x[0] - mean) / sigma;
    double gaus = TMath::Gaus(x[0], mean, sigma, true);
    double erf  = TMath::Erf(alpha * t / TMath::Sqrt2());

    return norm * gaus * (1 + erf);
}

void fit_adc() {

    // ---- Load histogram ----
    TFile *f = new TFile("adc_time_hist.root", "READ");
    if (!f || f->IsZombie()) { 
        std::cerr << "❌ Cannot open adc_time_hist.root\n"; 
        return; 
    }

    TH1F *h = (TH1F*)f->Get("h_adc");
    if (!h) { std::cerr << "❌ Histogram h_adc not found\n"; return; }

    // ---- Initial statistics ----
    double n    = h->GetEntries();
    double rms  = h->GetRMS();
    double mean = h->GetMean();

    double ymin = mean - 2*rms;
    if (ymin < 35) ymin = 35;
    double ymax = mean + 1.5*rms;

    std::cout << "Histogram mean=" << mean << " rms=" << rms
              << "   Fit Range: " << ymin << " to " << ymax << std::endl;

    // ---- Define skew-normal ADC fit ----
    TF1 *fadc = new TF1("adc_fun", skewnormal, ymin, ymax, 4);

    // If MDT — apply explicit fit window
    bool IsMDT = true;
    if (IsMDT) {
        fadc->SetRange(100, 300);
    }

    // Set initial params
    fadc->SetParameters(n, mean, rms, 1.5);

    // Parameter limits
    fadc->SetParLimits(1, 100, 300);       // mean window
    fadc->SetParLimits(2, 0.5*rms, 1.5*rms);
    fadc->SetParLimits(3, 1.0, 4.0);       // skew range

    fadc->SetLineColor(kRed);

    // ---- Perform fit ----
    h->Fit("adc_fun", "R");

    // ---- Check skew parameter and refit if needed ----
    
    if (fadc->GetParameter(3) > 2.5) {
        std::cout << "⚠️ Refitting with fixed skew=1.6 ..." << std::endl;
        fadc->SetParameter(3, 1.6);
        fadc->SetParLimits(3, 1.6, 1.6);
        h->Fit("adc_fun", "R");
    }
    
    // ---- Show results ----
    double fitMean  = fadc->GetParameter(1);
    double fitSigma = fadc->GetParameter(2);
    double skew     = fadc->GetParameter(3);

    std::cout << "✅ Fit results:" << std::endl;
    std::cout << "Mean = " << fitMean 
              << "   Sigma = " << fitSigma
              << "   Skew = " << skew << std::endl;

    // ---- Plot ----
    TCanvas *c = new TCanvas("c_adc", "ADC Fit", 800, 600);
    h->Draw();
    fadc->Draw("same");
    c->SaveAs("adc_fit.png");

    std::cout << "📁 Output saved: adc_fit.png" << std::endl;
}

