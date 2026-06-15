/*

This code describes the effects of birefringence at the South Pole

*/
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TFile.h"
#include "TText.h"
#include "TGaxis.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "TRotation.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <string>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TSystem.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>

// Local / project includes
#include "IceRayTracing.h"
//#include "/data/user/alansalgo/ARA/IceRayTracing/IceRayTracing.h"
#include "TLegend.h"
#include "birefringence.hh"
// Refactoring into .hh files
#include "constants.hh"
#include "config.hh"
#include "geometry.hh"
#include "station_data.hh"
#include "ice_profile.hh"
#include "antenna_measurements.hh"
//#include "plotting.hh"


using namespace std;

// ----------------------------------------------------------------------------
// Constants and global variables
// ----------------------------------------------------------------------------

// Principal indices used for a simple "toy" dielectric tensor if not depth-dependent
// (often interpreted as principal axes of the indicatrix)
vector<double> nvec{1.77,1.7805,1.7815};

// ROOT color palettes used for plotting
int icolors[6]={kGreen+2,kRed+1,kOrange+1,kViolet+1,kBlue,kBlack};
int icolors_dave[6]={kBlack,kRed,kGreen,kBlue,kYellow,kBlack};

// -----------------------------------------
// Free helper functions (geometry / plotting style) live in cpol_helpers.
// The core per-station physics loop lives in psi_model.
// -----------------------------------------
#include "cpol_helpers.hh"
#include "psi_model.hh"

// MACHTAY define some colors for printing text
#define RESET   "\x1b[0m"
#define BOLD    "\x1b[1m"
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define CYAN    "\x1b[36m"

// ----------------------------------------------------------------------------
// main
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
    cout << fixed << setprecision(8);
  // Read the terminal inputs
  Config cfg(argc, argv);   // Get the run details
  Geometry geom;            // Set up the geometry of stations and pulsers

	using clock = std::chrono::system_clock;
	using sec = std::chrono::duration<double>;
	const auto start = clock::now();
	// MACHTAY: to check speed, let's turn off print statements:
	//cout.setstate(ios_base::failbit);

	// MACHTAY: skip middle ray-path steps (only compute TX and RX endpoints).
	// Set true to skip intermediate istep iterations -- saves time but zeroes
	// out sumphase/model_data.vtimediff. Set false for full stepping.
	const bool SKIP_MIDDLE_STEPS = true;

	// MACHTAY: skip GetFull*RayPath entirely -- build rhat_thisstep directly
	// from launch/receive angles. Fastest mode; only polarization rotations
	// (epsilon/Psi/Omega) are valid. sumphase, atten, beam factors are zeroed.
	const bool SKIP_RAYTRACE_PATH = true;

    // Get the station geometry using these nvecs
    geom.compute_station_geometry(cfg.BIAXIAL, nvec, cfg.phi, cfg.theta, cfg.gamma);
    // Initialize lots of vectors for data
    StationData model_data;
    // Read in the index of refraction data and smooth it
    Ice_Profile ice(cfg, nvec);
    // Read in the actual data measured in A1-A5
    load_antenna_measurements(cfg, geom, model_data);

    // ------------------------------------------------------------------------
    // Build the pulser-depth scan used in the model.
    // For all stations, use a common grid:
    //   start at -400 m
    //   step downward by 10 m
    //   take 200 depth values
    //
    // model_data.vdepth[station] stores the pulser depths,
    // model_data.videpth[station] stores the corresponding integer-like indices.
    // ------------------------------------------------------------------------
    double arianna_minpulserdepth=-400.;
    double arianna_pulserstep=10.;
    int NARIANNA_PULSER=200;
    for (int istations=0;istations<geom.NSTATIONS;istations++) {
        for (int i=0;i<NARIANNA_PULSER;i++) {
            model_data.vdepth[istations].push_back(arianna_minpulserdepth-(double)arianna_pulserstep*(double)i);
            model_data.videpth[istations].push_back((double)i);
        }
    }

    // Build a reversed-depth version, useful for plots where increasing pulser depth
    // (positive number) is more natural than negative height.
    for (int istations=geom.minstation;istations<=geom.maxstation;istations++) {
        for (int j=0;j<model_data.vdepth[istations].size();j++) {
            model_data.vreversedepth[istations].push_back(-1.*model_data.vdepth[istations][model_data.vdepth[istations].size()-1-j]);
        }
    }

    // ------------------------------------------------------------------------
    // ROOT graphs:
    //   model_data.gsnrmax[i] = measured SNR vs pulser depth for station i
    //   model_data.g_idepth[i] = mapping from physical depth -> index in the depth scan
    //
    // Note: station index 5 (ARIANNA) has no model_data.gsnrmax built here.
    // ------------------------------------------------------------------------
    for (int i=geom.minstation;i<=geom.maxstation;i++) {
        if (i!=5)
        model_data.gsnrmax[i]=new TGraph(model_data.vsnrmax[i].size(),&model_data.vdepth_data[i][0],&model_data.vsnrmax[i][0]);
        model_data.g_idepth[i]=new TGraph(model_data.vdepth[i].size(),&model_data.vdepth[i][0],&model_data.videpth[i][0]);


    }

    // ------------------------------------------------------------------------
    // For each station, find:
    //   igreatestdepth   = index of deepest pulser depth (most negative)
    //   imostshallowdepth= index of shallowest pulser depth (least negative)
    // ------------------------------------------------------------------------
    for (int istations=geom.minstation;istations<=geom.maxstation;istations++) {
        geom.igreatestdepth[istations]=0;
        geom.imostshallowdepth[istations]=0;
        for (int idepth=0;idepth<model_data.vdepth[istations].size();idepth++) {
            if (model_data.vdepth[istations][idepth]<model_data.vdepth[istations][geom.igreatestdepth[istations]])
            geom.igreatestdepth[istations]=idepth;
            if (model_data.vdepth[istations][idepth]>model_data.vdepth[istations][geom.imostshallowdepth[istations]])
            geom.imostshallowdepth[istations]=idepth;
        }

    }

    // ------------------------------------------------------------------------
    // Distance-domain "big picture" scan parameters (per-station pseudo-depths).
    // ------------------------------------------------------------------------
    const int NDISTANCES_BIGPIC=5000;
    const double STEP=1.;

    int jmin[geom.NSTATIONS];

    // ------------------------------------------------------------------------
    // Main per-station physics loop.
    //
    // For each station:
    //   - allocate per-depth arrays
    //   - build a distance-domain "big picture" scan
    //   - ray trace from each pulser depth to the station
    //   - follow the ray through the ice
    //   - accumulate attenuation and birefringent phase
    //   - compute TX/RX polarization projections and beam factors
    // ------------------------------------------------------------------------    
    // MACHTAY: time the idepth loop
    PsiTiming timing;
    double n_e1 = 0.0;
    double n_e2 = 0.0;
    TVector3 p_e1;
    TVector3 p_e2;
    for (int i=geom.minstation;i<=geom.maxstation;i++) {
        process_station(cfg, geom, ice, model_data, i,
                        p_e1, p_e2, n_e1, n_e2,
                        SKIP_MIDDLE_STEPS, SKIP_RAYTRACE_PATH, timing);
    }

	    string filename;

	    // MACHTAY: End function here (skip plotting)
    gROOT->SetBatch(true);  // prevent hang in headless environment
    cout.clear();
    const sec duration = clock::now() - start;
    std::cout << "Uzair time: " << timing.uzair_time << "s" << std::endl;
    std::cout << "Overhead time: " << timing.overhead << "ms" << std::endl;
    std::cout << "Eval time: " << timing.eval_time << "s" << std::endl;
    std::cout << "Raytracing time: " << timing.ray_time << "s" << std::endl;
    std::cout << "idepth loop duration: " << timing.loop_time << "s" << std::endl;
    std::cout << "Total duration: " << duration.count() << "s" << std::endl;

    // MACHTAY: Reproduce Fig. 6 -- epsilon_1^T and epsilon_1^R vs pulser height, all stations
    //plot_epsilons(geom, model_data);

    plot_epsilons_tx_rx(geom, model_data);
    plot_epsilons_differences(geom, model_data);

    return 0;
}
