#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>

void drift_hist() {
    // Histogram settings
    int nbins = (1800) / 3;
    TH1F *h_drift = new TH1F("h_drift", "drift Time;drift Time [ns];Counts", nbins, 0, 1800);

    // Open CSV
    std::ifstream file("hits_0.csv");
    if (!file.is_open()) {
        std::cerr << "ERROR: File not found: perpendicular_hits_all_pairs.csv\n";
        return;
    }

    std::string line;
    bool first = true;

    // Loop through rows
    while (std::getline(file, line)) {
        if (first) { first = false; continue; } // skip header

        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> cols;

        while (std::getline(ss, item, ',')) cols.push_back(item);

        if (cols.size() < 1) continue;

        double drift_time = atof(cols[15].c_str()); // adjust if drift_time column index is different
        h_drift->Fill(drift_time);
    }

    file.close();

    // Save ROOT file
    TFile *outfile = new TFile("drift_time_hist.root", "RECREATE");
    h_drift->Write();
    outfile->Close();

    // Draw histogram
    TCanvas *c = new TCanvas("c", "Drift Time", 800, 600);
    h_drift->Draw();
    c->SaveAs("drift_time_hist.png");

   // std::cout << "✅ Histogram saved to adc_time_hist.root and adc_time_hist.png\n";
}

