// test_psimodel_station_scan_range.cc
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <algorithm>

#include "cpol_fit.h"

using namespace std;

int main(int argc, char** argv) {
    // Expect:
    // argv[1]  = start station (1..5)
    // argv[2]  = end station   (1..5)
    // argv[3]..argv[18] = 16 fit parameters
    if (argc < 19) {
        cerr << "Usage:\n"
             << "  " << argv[0]
             << " start_station end_station"
             << " p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16\n";
        return 1;
    }

    // ------------------------------------------------------------------------
    // Read station range
    // ------------------------------------------------------------------------
    int start_station_human = stoi(argv[1]);
    int end_station_human   = stoi(argv[2]);

    if (start_station_human < 1 || start_station_human > 5 ||
        end_station_human   < 1 || end_station_human   > 5) {
        cerr << "Error: station numbers must be between 1 and 5.\n";
        return 1;
    }

    if (start_station_human > end_station_human) {
        std::swap(start_station_human, end_station_human);
    }

    // Internal station convention: 0=A1, 1=A2, 2=A3, 3=A4, 4=A5
    int start_station = start_station_human - 1;
    int end_station   = end_station_human - 1;

    // ------------------------------------------------------------------------
    // Read fit parameters
    // ------------------------------------------------------------------------
    const int NPAR = 16;
    Double_t fit_params[NPAR];
    for (int i = 0; i < NPAR; i++) {
        fit_params[i] = stod(argv[i + 3]);
    }

    cout << "Station range: A" << start_station_human
         << " to A" << end_station_human << "\n";

    cout << "Input fit parameters:\n";
    for (int i = 0; i < NPAR; i++) {
        cout << "  par[" << i << "] = " << fit_params[i] << "\n";
    }

    // ------------------------------------------------------------------------
    // Build depth/height array
    // Use pulser height from -1600 m to -600 m in steps of 10 m
    // If psiModel expects positive depth, convert with depth = -height
    // ------------------------------------------------------------------------
    const double height_min  = -1600.0;
    const double height_max  = -600.0;
    const double height_step = 10.0;

    const int nDepths = static_cast<int>((height_max - height_min) / height_step) + 1;

    Double_t* pulserDepth = new Double_t[nDepths];
    vector<double> pulserHeight(nDepths);

    for (int i = 0; i < nDepths; i++) {
        pulserHeight[i] = height_min + i * height_step;  // -1600 ... -600
        pulserDepth[i]  = -pulserHeight[i];              // 1600 ... 600
    }

    // ------------------------------------------------------------------------
    // Colors and labels
    // ------------------------------------------------------------------------
    int colors[5] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+2};
    string labels[5] = {"A1", "A2", "A3", "A4", "A5"};

    vector<TGraph*> graphs;

    double global_ymin = 1e99;
    double global_ymax = -1e99;

    // ------------------------------------------------------------------------
    // Run psiModel for each requested station
    // ------------------------------------------------------------------------
    for (int station = start_station; station <= end_station; station++) {
        Double_t* psi_model = new Double_t[nDepths];
        for (int i = 0; i < nDepths; i++) psi_model[i] = 0.0;

        int result = psiModel(argc, argv, fit_params, pulserDepth, psi_model, station, nDepths);

        if (result != 0) {
            cerr << "psiModel returned non-zero status for station A"
                 << (station + 1) << ": " << result << endl;
            delete[] psi_model;
            continue;
        }

        // Build graph versus pulser height
        TGraph* g = new TGraph(nDepths, pulserHeight.data(), psi_model);
        g->SetLineColor(colors[station]);
        g->SetMarkerColor(colors[station]);
        g->SetLineWidth(3);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.6);
        g->SetName(labels[station].c_str());

        for (int i = 0; i < nDepths; i++) {
            if (psi_model[i] < global_ymin) global_ymin = psi_model[i];
            if (psi_model[i] > global_ymax) global_ymax = psi_model[i];
        }

        graphs.push_back(g);

        cout << "\nStation A" << (station + 1) << ":\n";
        cout << "Height (m)    Depth (m)    Psi_model\n";
        cout << "------------------------------------\n";
        for (int i = 0; i < nDepths; i++) {
            cout << fixed << setprecision(2)
                 << setw(10) << pulserHeight[i] << "    "
                 << setw(10) << pulserDepth[i]  << "    "
                 << setw(10) << psi_model[i]    << "\n";
        }

        delete[] psi_model;
    }

    if (graphs.empty()) {
        cerr << "No valid graphs were produced.\n";
        delete[] pulserDepth;
        return 2;
    }

    // ------------------------------------------------------------------------
    // Plot all stations together
    // ------------------------------------------------------------------------
    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "psiModel scan by station", 1000, 700);
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    c1->SetGrid();

    double yrange = global_ymax - global_ymin;
    if (yrange <= 0) yrange = 1.0;

    TH1F* frame = new TH1F("frame", "", 100, height_min, height_max);
    frame->SetMinimum(global_ymin - 0.1 * yrange);
    frame->SetMaximum(global_ymax + 0.1 * yrange);
    frame->GetXaxis()->SetTitle("Pulser height (m)");
    frame->GetYaxis()->SetTitle("#psi_{model} (deg)");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetTitleOffset(1.2);
    frame->Draw();

    TLegend* leg = new TLegend(0.72, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);

    for (auto* g : graphs) {
        g->Draw("LP SAME");
        leg->AddEntry(g, g->GetName(), "lp");
    }

    leg->Draw();

    c1->SaveAs("psiModel_station_scan_range.pdf");
    c1->SaveAs("psiModel_station_scan_range.png");

    // ------------------------------------------------------------------------
    // Cleanup
    // ------------------------------------------------------------------------
    delete leg;
    delete frame;
    for (auto* g : graphs) delete g;
    delete c1;
    delete[] pulserDepth;

    return 0;
}