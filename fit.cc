#include <TROOT.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMinuit.h>
#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TStyle.h>

#include "cpol_fit.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// ============================================================
// User configuration
// ============================================================
static const double FIXED_DELTA_DEG = 10.0;   // par[9], fixed
static const int NPOLY_PARAMS = 9;            // vary only par[0]..par[8]
static const int NMODEL_PARAMS = 10;          // par[0]..par[9], with par[9]=delta fixed
static const int NSTARTS = 1;                // number of parallel Minuit fits
static const int FIT_MODE = 0;                // 0 = use measurement variances
static const int MAX_MINUIT_ITER = 2000;

// Output files
static const string OUT_SUMMARY_TXT = "fit_summary.txt";
static const string OUT_ROOT_FILE   = "fit_results.root";

// A2 and A4 only
static const vector<string> DATA_FILES = {
    "/data/user/alansalgo/forAlex/spiceData_5meterResolution/A2_spiceReco.root",
    "/data/user/alansalgo/forAlex/spiceData_5meterResolution/A4_spiceReco.root"
};

// psiModel station convention:
// 0=A1, 1=A2, 2=A3, 3=A4, 4=A5
static const vector<int> STATION_INDICES = {1, 3};
static const vector<string> STATION_LABELS = {"A2", "A4"};

// ============================================================
// Data container
// ============================================================
struct StationData {
    string label;
    string filename;
    int station_index = -1;

    vector<Double_t> depth;
    vector<Double_t> psi;
    vector<Double_t> err;

    Long64_t nEntries = 0;
    double avg_err = 1.0;
};

struct FitResult {
    bool valid = false;
    double chi2 = numeric_limits<double>::infinity();
    int ndof = 0;
    array<double, NMODEL_PARAMS> pars{};
    array<double, NPOLY_PARAMS> errs{};
    int minuit_status = -1;
    int seed_id = -1;
};

// ============================================================
// Utilities
// ============================================================
double wrap360(double x) {
    double y = fmod(x, 360.0);
    if (y < 0) y += 360.0;
    return y;
}

// Your convention:
// angle measured from vertical, so:
// 10 = 170 = 190 = 350
double canonicalVerticalAngle(double psi_deg) {
    double a = wrap360(psi_deg);

    if (a > 180.0) a = 360.0 - a;
    if (a > 90.0) a = 180.0 - a;

    return fabs(a);
}

double wrappedResidualVertical(double model_deg, double data_deg) {
    return canonicalVerticalAngle(model_deg) - canonicalVerticalAngle(data_deg);
}

// For plotting as pulser "height" instead of depth:
// depth 600 -> height -600
// depth 1600 -> height -1600
double depthToHeight(double depth_m) {
    return -depth_m;
}

// ============================================================
// Read Justin's data
// ============================================================
bool extractPsi(const string& filename,
                vector<Double_t>& pulserDepth,
                vector<Double_t>& psi_median,
                vector<Double_t>& psi_errors)
{
    TFile* file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << filename << endl;
        return false;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("polReco"));
    if (!tree) {
        cerr << "Error getting tree 'polReco' from file: " << filename << endl;
        file->Close();
        return false;
    }

    Long64_t nEntries = tree->GetEntries();
    if (nEntries <= 0) {
        cerr << "No entries found in tree: " << filename << endl;
        file->Close();
        return false;
    }

    Double_t pulserDepthValue = 0.0;
    Double_t psi_medianValue = 0.0;
    Double_t psi_upperLimitValue = 0.0;

    tree->SetBranchAddress("pulserDepth", &pulserDepthValue);
    tree->SetBranchAddress("psi_median", &psi_medianValue);
    tree->SetBranchAddress("psi_upperLimit", &psi_upperLimitValue);

    if (filename.find("A3_spiceReco.root") != string::npos ||
        filename.find("A5_spiceReco.root") != string::npos) {
        nEntries -= 1;
    }

    pulserDepth.resize(nEntries);
    psi_median.resize(nEntries);
    psi_errors.resize(nEntries);

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        pulserDepth[i] = pulserDepthValue;
        psi_median[i] = psi_medianValue;
        psi_errors[i] = fabs(psi_medianValue - psi_upperLimitValue);
    }

    file->Close();
    return true;
}

// ============================================================
// Build station data
// ============================================================
vector<StationData> loadStations() {
    vector<StationData> stations;

    for (size_t i = 0; i < DATA_FILES.size(); ++i) {
        StationData s;
        s.label = STATION_LABELS[i];
        s.filename = DATA_FILES[i];
        s.station_index = STATION_INDICES[i];

        if (!extractPsi(s.filename, s.depth, s.psi, s.err)) {
            cerr << "Failed to load station " << s.label << endl;
            continue;
        }

        s.nEntries = static_cast<Long64_t>(s.depth.size());

        double sum_err = 0.0;
        int count_err = 0;
        for (double e : s.err) {
            if (e > 0) {
                sum_err += e;
                count_err++;
            }
        }
        s.avg_err = (count_err > 0) ? (sum_err / count_err) : 5.0;

        cout << "Loaded " << s.label
             << " with " << s.nEntries
             << " entries, avg_err = " << s.avg_err << endl;

        stations.push_back(std::move(s));
    }

    return stations;
}

// ============================================================
// Evaluate total chi2
// ============================================================
double evaluateTotalChi2(const vector<StationData>& stations,
                         const double* pars10,
                         int argc,
                         char** argv,
                         int& totalPoints)
{
    totalPoints = 0;
    double totalChi2 = 0.0;

    for (const auto& s : stations) {
        vector<Double_t> model(s.nEntries, 0.0);

        Double_t* depth_ptr = const_cast<Double_t*>(s.depth.data());
        Double_t* model_ptr = model.data();

        int rc = psiModel(argc, argv,
                          pars10,
                          depth_ptr,
                          model_ptr,
                          s.station_index,
                          s.nEntries);

        if (rc != 0) {
            return numeric_limits<double>::infinity();
        }

        for (Long64_t i = 0; i < s.nEntries; ++i) {
            double sigma = s.err[i] > 0 ? s.err[i] : s.avg_err;
            double resid = wrappedResidualVertical(model[i], s.psi[i]);

            if (FIT_MODE == 0) {
                totalChi2 += (resid * resid) / (sigma * sigma);
            } else {
                double denom = fabs(model[i]);
                if (denom <= 0) denom = sigma * sigma;
                totalChi2 += (resid * resid) / denom;
            }
        }

        totalPoints += static_cast<int>(s.nEntries);
    }

    return totalChi2;
}

// ============================================================
// Evaluate model for one station
// ============================================================
bool evaluateStationModel(const StationData& s,
                          const double* pars10,
                          int argc,
                          char** argv,
                          vector<double>& model_out)
{
    model_out.assign(s.nEntries, 0.0);

    Double_t* depth_ptr = const_cast<Double_t*>(s.depth.data());
    Double_t* model_ptr = model_out.data();

    int rc = psiModel(argc, argv,
                      pars10,
                      depth_ptr,
                      model_ptr,
                      s.station_index,
                      s.nEntries);

    return (rc == 0);
}

// ============================================================
// Minuit context
// ============================================================
struct MinuitContext {
    const vector<StationData>* stations = nullptr;
    int argc = 0;
    char** argv = nullptr;
    double fixed_delta = FIXED_DELTA_DEG;
};

static thread_local MinuitContext* gMinuitContext = nullptr;

// ============================================================
// FCN for Minuit
// ============================================================
void minuitFCN(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t iflag)
{
    (void)npar;
    (void)grad;
    (void)iflag;

    if (!gMinuitContext || !gMinuitContext->stations) {
        fval = 1e30;
        return;
    }

    double fullPars[NMODEL_PARAMS] = {0};

    for (int i = 0; i < NPOLY_PARAMS; ++i) {
        fullPars[i] = par[i];
    }
    fullPars[9] = gMinuitContext->fixed_delta;

    int totalPoints = 0;
    fval = evaluateTotalChi2(*(gMinuitContext->stations),
                             fullPars,
                             gMinuitContext->argc,
                             gMinuitContext->argv,
                             totalPoints);

    if (!std::isfinite(fval)) fval = 1e30;
}

// ============================================================
// One Minuit fit from one start
// ============================================================
FitResult runOneFit(const vector<StationData>& stations,
                    int argc,
                    char** argv,
                    int seed_id)
{
    FitResult out;
    out.seed_id = seed_id;

    MinuitContext ctx;
    ctx.stations = &stations;
    ctx.argc = argc;
    ctx.argv = argv;
    ctx.fixed_delta = FIXED_DELTA_DEG;

    gMinuitContext = &ctx;

    TMinuit minuit(NPOLY_PARAMS);
    minuit.SetFCN(minuitFCN);
    minuit.SetPrintLevel(-1);
    minuit.SetErrorDef(1.0);
    minuit.SetMaxIterations(MAX_MINUIT_ITER);

    const double lo = -90.0;
    const double hi =  90.0;

    srand(12345 + 97 * seed_id);

    for (int i = 0; i < NPOLY_PARAMS; ++i) {
        double init = lo + (hi - lo) * (rand() / (double)RAND_MAX);
        double step = 1.0;
        string pname = "p" + to_string(i);
        minuit.DefineParameter(i, pname.c_str(), init, step, lo, hi);
    }

    int ierr = 0;
    minuit.mnexcm("MIGRAD", nullptr, 0, ierr);

    Double_t fmin, fedm, errdef;
    Int_t npari, nparx, istat;
    minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);

    out.minuit_status = istat;

    if (!std::isfinite(fmin)) {
        out.valid = false;
        return out;
    }

    for (int i = 0; i < NPOLY_PARAMS; ++i) {
        Double_t val = 0.0, err = 0.0;
        minuit.GetParameter(i, val, err);
        out.pars[i] = val;
        out.errs[i] = err;
    }

    out.pars[9] = FIXED_DELTA_DEG;

    int totalPoints = 0;
    out.chi2 = evaluateTotalChi2(stations, out.pars.data(), argc, argv, totalPoints);
    out.ndof = totalPoints - NPOLY_PARAMS;
    out.valid = std::isfinite(out.chi2) && out.ndof > 0;

    return out;
}

// ============================================================
// Print best result
// ============================================================
void printBestResult(const FitResult& best)
{
    cout << "\n================ BEST FIT ================\n";
    cout << "seed_id   : " << best.seed_id << "\n";
    cout << "valid     : " << best.valid << "\n";
    cout << "chi2      : " << best.chi2 << "\n";
    cout << "ndof      : " << best.ndof << "\n";
    cout << "chi2/ndof : " << (best.ndof > 0 ? best.chi2 / best.ndof : -1) << "\n";
    cout << "status    : " << best.minuit_status << "\n\n";

    for (int i = 0; i < NPOLY_PARAMS; ++i) {
        cout << "par[" << i << "] = "
             << setw(12) << best.pars[i]
             << " +/- " << best.errs[i] << "\n";
    }
    cout << "par[9]  = " << best.pars[9] << "   (fixed delta)\n";
    cout << "==========================================\n";
}

// ============================================================
// Save text summary
// ============================================================
void saveTextSummary(const string& outname,
                     const vector<StationData>& stations,
                     const vector<FitResult>& results,
                     const FitResult& best,
                     double runtime_sec)
{
    ofstream fout(outname);
    if (!fout) {
        cerr << "Could not open " << outname << " for writing.\n";
        return;
    }

    fout << fixed << setprecision(6);

    fout << "==================================================\n";
    fout << "Fit summary\n";
    fout << "==================================================\n";
    fout << "Fixed delta (par[9]) = " << FIXED_DELTA_DEG << " deg\n";
    fout << "Number of starts     = " << NSTARTS << "\n";
    fout << "FIT_MODE             = " << FIT_MODE << "\n";
    fout << "Runtime (s)          = " << runtime_sec << "\n\n";

    fout << "Stations used:\n";
    for (const auto& s : stations) {
        fout << "  " << s.label
             << "  file=" << s.filename
             << "  nEntries=" << s.nEntries
             << "  avg_err=" << s.avg_err << "\n";
    }
    fout << "\n";

    fout << "================ BEST FIT ================\n";
    fout << "seed_id   = " << best.seed_id << "\n";
    fout << "valid     = " << best.valid << "\n";
    fout << "chi2      = " << best.chi2 << "\n";
    fout << "ndof      = " << best.ndof << "\n";
    fout << "chi2/ndof = " << (best.ndof > 0 ? best.chi2 / best.ndof : -1.0) << "\n";
    fout << "status    = " << best.minuit_status << "\n";
    for (int i = 0; i < NPOLY_PARAMS; ++i) {
        fout << "par[" << i << "] = " << best.pars[i]
             << " +/- " << best.errs[i] << "\n";
    }
    fout << "par[9] = " << best.pars[9] << " (fixed delta)\n\n";

    fout << "================ ALL STARTS ================\n";
    fout << "seed_id valid chi2 ndof chi2_ndof status";
    for (int i = 0; i < NMODEL_PARAMS; ++i) fout << " par" << i;
    fout << "\n";

    for (const auto& r : results) {
        fout << r.seed_id << " "
             << r.valid << " "
             << r.chi2 << " "
             << r.ndof << " "
             << (r.ndof > 0 ? r.chi2 / r.ndof : -1.0) << " "
             << r.minuit_status;

        for (int i = 0; i < NMODEL_PARAMS; ++i) {
            fout << " " << r.pars[i];
        }
        fout << "\n";
    }

    fout.close();
    cout << "Saved summary text to " << outname << endl;
}

// ============================================================
// Save ROOT outputs + overlay plots
// ============================================================
void saveROOTAndPlots(const string& rootname,
                      const vector<StationData>& stations,
                      const FitResult& best,
                      int argc,
                      char** argv)
{
    TFile* fout = new TFile(rootname.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        cerr << "Could not create ROOT file " << rootname << endl;
        return;
    }

    gStyle->SetOptStat(0);

    for (const auto& s : stations) {
        vector<double> model;
        if (!evaluateStationModel(s, best.pars.data(), argc, argv, model)) {
            cerr << "Could not evaluate best-fit model for " << s.label << endl;
            continue;
        }

        vector<double> x_height(s.nEntries);
        vector<double> y_data(s.nEntries);
        vector<double> y_model(s.nEntries);
        vector<double> y_err(s.nEntries);
        vector<double> y_resid(s.nEntries);

        for (Long64_t i = 0; i < s.nEntries; ++i) {
            x_height[i] = depthToHeight(s.depth[i]);
            y_data[i]   = canonicalVerticalAngle(s.psi[i]);
            y_model[i]  = canonicalVerticalAngle(model[i]);
            y_err[i]    = (s.err[i] > 0 ? s.err[i] : s.avg_err);
            y_resid[i]  = wrappedResidualVertical(model[i], s.psi[i]);
        }

        // Reverse so -1600 is on the left and -600 on the right
        std::reverse(x_height.begin(), x_height.end());
        std::reverse(y_data.begin(), y_data.end());
        std::reverse(y_model.begin(), y_model.end());
        std::reverse(y_err.begin(), y_err.end());
        std::reverse(y_resid.begin(), y_resid.end());

        fout->mkdir(s.label.c_str());
        fout->cd(s.label.c_str());

        TGraphErrors* gData = new TGraphErrors(
            (int)s.nEntries,
            x_height.data(),
            y_data.data(),
            nullptr,
            y_err.data()
        );
        gData->SetName((s.label + "_data").c_str());
        gData->SetTitle((s.label + " data").c_str());
        gData->SetMarkerStyle(20);
        gData->SetMarkerSize(1.0);
        gData->GetXaxis()->SetTitle("Pulser height (m)");
        gData->GetYaxis()->SetTitle("#psi (deg)");

        TGraph* gModel = new TGraph(
            (int)s.nEntries,
            x_height.data(),
            y_model.data()
        );
        gModel->SetName((s.label + "_model_bestfit").c_str());
        gModel->SetTitle((s.label + " best-fit model").c_str());
        gModel->SetLineWidth(3);
        gModel->SetLineColor(kRed);
        gModel->GetXaxis()->SetTitle("Pulser height (m)");
        gModel->GetYaxis()->SetTitle("#psi (deg)");

        TGraph* gResid = new TGraph(
            (int)s.nEntries,
            x_height.data(),
            y_resid.data()
        );
        gResid->SetName((s.label + "_residuals").c_str());
        gResid->SetTitle((s.label + " residuals").c_str());
        gResid->SetLineWidth(2);
        gResid->SetMarkerStyle(20);
        gResid->GetXaxis()->SetTitle("Pulser height (m)");
        gResid->GetYaxis()->SetTitle("model - data (deg)");

        gData->Write();
        gModel->Write();
        gResid->Write();

        // Overlay plot
        TCanvas* c1 = new TCanvas((s.label + "_c1").c_str(), (s.label + " overlay").c_str(), 900, 700);
        c1->SetLeftMargin(0.12);
        c1->SetBottomMargin(0.12);
        c1->SetGrid();

        gData->SetTitle("");
        gData->GetXaxis()->SetTitle("Pulser height (m)");
        gData->GetYaxis()->SetTitle("#psi (deg)");
        gData->Draw("AP");

        gModel->Draw("L SAME");

        TLegend* leg = new TLegend(0.60, 0.75, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->AddEntry(gData, (s.label + " data").c_str(), "lep");
        leg->AddEntry(gModel, (s.label + " best fit").c_str(), "l");
        leg->Draw();

        c1->Write((s.label + "_overlay_canvas").c_str());
        c1->SaveAs((s.label + "_bestfit_overlay.pdf").c_str());
        c1->SaveAs((s.label + "_bestfit_overlay.png").c_str());

        // Residual plot
        TCanvas* c2 = new TCanvas((s.label + "_c2").c_str(), (s.label + " residuals").c_str(), 900, 700);
        c2->SetLeftMargin(0.12);
        c2->SetBottomMargin(0.12);
        c2->SetGrid();
        gResid->Draw("ALP");
        c2->Write((s.label + "_residual_canvas").c_str());

        delete leg;
        delete c1;
        delete c2;

        fout->cd();
    }

    fout->Close();
    delete fout;

    cout << "Saved ROOT outputs to " << rootname << endl;
}

// ============================================================
// Main
// ============================================================
int main(int argc, char** argv)
{
    // if (argc < 10) {
    //     cerr << "Usage:\n"
    //          << "  " << argv[0]
    //          << " gamma0 gamma1 gamma2 alpha0 alpha1 alpha2 beta0 beta1 beta2\n\n"
    //          << "These 9 arguments are kept so psiModel gets argv in the format it expects.\n"
    //          << "This program runs many random Minuit starts in parallel.\n";
    //     return 1;
    // }

    auto t0 = chrono::high_resolution_clock::now();

    cout << "Loading A2 and A4 only\n";
    cout << "Fixed delta = " << FIXED_DELTA_DEG << " deg\n";
    cout << "Number of random starts = " << NSTARTS << "\n";

    vector<StationData> stations = loadStations();
    if (stations.empty()) {
        cerr << "No station data loaded.\n";
        return 2;
    }

    vector<FitResult> results(NSTARTS);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < NSTARTS; ++i) {
        results[i] = runOneFit(stations, argc, argv, i);
    }

    FitResult best;
    for (const auto& r : results) {
        if (!r.valid) continue;
        if (r.chi2 < best.chi2) best = r;
    }

    if (!best.valid) {
        cerr << "No valid Minuit fit found.\n";
        return 3;
    }

    printBestResult(best);

    auto t1 = chrono::high_resolution_clock::now();
    double dt = chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() / 1000.0;

    saveTextSummary(OUT_SUMMARY_TXT, stations, results, best, dt);
    saveROOTAndPlots(OUT_ROOT_FILE, stations, best, argc, argv);

    cout << "Total runtime: " << dt << " s\n";
    return 0;
}