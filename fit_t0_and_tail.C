#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

// t0 Fermi function: f(t) = A / (1 + exp((B-t)/C)) + k
double mt_t0_fermi(double *x, double *par)
{
    double t = x[0];
    double t0   = par[0];
    double slope = par[1];
    double back  = par[2];
    double ampl  = par[3];
    return ampl / (1.0 + exp((t0 - t)/slope)) + back;
}

// Tail function: (A + D*t) / (1 + exp((B - t)/C)) + k
double fermi_tail(double *x, double *par)
{
    double t = x[0];
    double A = par[0];
    double B = par[1];
    double C = par[2];
    double D = par[3];
    double k = par[4];
    return (A + D*t) / (1.0 + exp((B - t)/C)) + k;
}

void fit_t0_and_tail()
{
    // Open ROOT file and histogram
    TFile *f = TFile::Open("drift_time_hist.root");
    if(!f || f->IsZombie()) { std::cout << "❌ File not found\n"; return; }

    TH1F *h = (TH1F*)f->Get("h_drift");
    if(!h) { std::cout << "❌ Histogram h_drift not found\n"; return; }

    // --- 1) T0 FIT -----------------------------------------------------

    double maxval = h->GetMaximum();
    double t0     = h->GetBinCenter(h->FindFirstBinAbove(0.45*maxval));
    double slope  = 2.5;
    double ampl   = maxval/1.1;
    double xmin   = t0 - 300.;
    double xmax   = t0 + 35.;

    TF1 *ft0 = new TF1("ft0", mt_t0_fermi, xmin, xmax, 4);
    ft0->SetParameters(t0, slope, 0., ampl);
    ft0->SetParLimits(2,0.,ampl);
    ft0->SetLineColor(kRed);

    h->Fit("ft0","R"); // quiet + range

    double B0 = ft0->GetParameter(0);  // B = t0 location
    double DBrise = ft0->GetParError(0);
    // if background k < 0, refit with fixed k=0
    if(ft0->GetParameter(2) < 0) {
        ft0->FixParameter(2,0.);
        h->Fit("ft0","R");
    }

    // --- 2) TAIL FIT ---------------------------------------------------

    double xmin2 = 900;
    double xmax2 = 1300;

    TF1 *ftail = new TF1("ftail", fermi_tail, xmin2, xmax2, 5);
    ftail->SetParNames("A","B","C","D","k");
    ftail->SetParameters(480, 1130, -6.3, -0.55, 1.4);
    ftail->SetLineColor(kBlue);
    ftail->SetNpx(1000);

    h->Fit("ftail","R");

    double Btail = ftail->GetParameter(1);
    double DBtail =  ftail->GetParError(1);
    double tmax  = Btail - B0;
    double DB= sqrt(DBtail*DBtail + DBrise*DBrise);

    // --- OUTPUT --------------------------------------------------------

    printf("\n=========== Fit Results ===========\n");
    printf("t0 fit B0       = %.3f ns\n", B0);
    printf("Tail fit Btail  = %.3f ns\n", Btail);
    printf("-----------------------------------\n");
    printf("tmax = Btail - B0 = %.3f ns\n", tmax);
    printf("DB=%.3f ns\n",DB);
    printf("===================================\n");

    // Plot & save
    TCanvas *c = new TCanvas("c","T0 and Tail Fit",1200,800);
    h->Draw();
    ft0->Draw("same");
    ftail->Draw("same");
    c->SaveAs("drift_t0_tail_fit.png");
}

