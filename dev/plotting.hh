#pragma once

#include "TF1.h"
#include "TGraph.h"

struct PlotData {
    TF1* f1[6] = {};
    TF1* f1_nointerference[6] = {};
    TF1* f1_distances[6] = {};
    TF1* f1_nointerference_distances[6] = {};

    TGraph* g1_distances[6] = {};
    TGraph* g_atten_beam_distances[6] = {};
    TGraph* g_atten_beam_crosspol_distances[6] = {};
    TGraph* g_atten_beam_crosspol_nointerferencefunc_distances[6] = {};
    TGraph* g_atten_beam_crosspol_func_distances[6] = {};
    TGraph* g_nsolutions_distances[6] = {};
    TGraph* g_sumphase_distances[6] = {};
    const int NSPECIALDEPTHS=2;

    double specialdepths[NSPECIALDEPTHS]={-850.,-875.};

    // g_spectra[station][depth_index] will hold the frequency spectrum
    vector< vector<TGraph*> > g_spectra;

    g_spectra.resize(geom.NSTATIONS);
    string sfunc;

    // Envelope with interference term
    // [0] is an overall scale
    // [1], [2] are the two contributions
    // [3] is the interference product term
    // [4] is the phase
    sfunc="[0]*sqrt(([1]+[2])*([1]+[2])-2*[3]*sin([4])*sin([4]))";

    string sfunc_nointerference;
    // Same envelope, but dropping interference
    sfunc_nointerference="[0]*sqrt(([1]+[2])*([1]+[2]))";


    string sfunc_distances;
    sfunc_distances="[0]*sqrt(([1]+[2])*([1]+[2])-2*[3]*sin([4])*sin([4]))";


    string sfunc_nointerference_distances;
    sfunc_nointerference_distances="[0]*sqrt(([1]+[2])*([1]+[2]))";

    // ------------------------------------------------------------------------
    // Graph containers for many derived quantities vs depth
    // ------------------------------------------------------------------------
    TGraph *g_parameter0[geom.NSTATIONS];
    TGraph *g_atten[geom.NSTATIONS];
    TGraph *g_atten_power[geom.NSTATIONS];

    TGraph *gfunc_noadjust[geom.NSTATIONS];
    TGraph *g_atten_beam[geom.NSTATIONS];
    TGraph *g_atten_beam_power[geom.NSTATIONS];
    TGraph *g_atten_beam_crosspol[geom.NSTATIONS];

    TGraph *g_atten_beam_crosspol_nointerferencefunc[geom.NSTATIONS];
    TGraph *g_atten_beam_crosspol_func[geom.NSTATIONS];
    TGraph *g_sumphase[geom.NSTATIONS];
    TGraph *g_notflipped[geom.NSTATIONS];
    TGraph *g_deltan[geom.NSTATIONS];
    TGraph *g_notflipped_alongpath[geom.NSTATIONS];
    TGraph *g_theta1_alongpath[geom.NSTATIONS];
    TGraph *g_theta2_alongpath[geom.NSTATIONS];
    TGraph *g_thetape1_alongpath[geom.NSTATIONS];
    TGraph *g_thetape2_alongpath[geom.NSTATIONS];
    TGraph *g_thetape1_phipe1_alongpath[geom.NSTATIONS];
    TGraph *g_thetape2_phipe2_alongpath[geom.NSTATIONS];
    TGraph *g_phipe1_alongpath[geom.NSTATIONS];
    TGraph *g_phipe2_alongpath[geom.NSTATIONS];
    TGraph *g_deltan_pulserdepth[geom.NSTATIONS];
    TGraph *g_depth_istep[geom.NSTATIONS];
    TGraph *g_path[geom.NSTATIONS];
    TGraph *g_sumlength[geom.NSTATIONS];
    TGraph *g_attenlengths[geom.NSTATIONS];

    // ------------------------------------------------------------------------
    // "Big picture" scans parameterized by total distance instead of pulser depth.
    // These are used later to make distance-domain versions of the model curves.
    // ------------------------------------------------------------------------
    vector< vector<double> > vdistances_bigpic;
    vector< vector<double> > vpseudodepths_bigpic;
    vector< vector<double> > vmag_atten_beam_bigpic;
    vector< vector<double> > vmag_atten_beam_crosspol_bigpic;
    vector< vector<double> > vmag_atten_beam_crosspol_nointerferencefunc_bigpic;
    vector< vector<double> > vmag_atten_beam_crosspol_func_bigpic;
    vector< vector<double> > vtimediff_bigpic;

    vector< vector<double> > vspectrum_bigpic;

    vmag_atten_beam_bigpic.resize(geom.NSTATIONS);
    vmag_atten_beam_crosspol_bigpic.resize(geom.NSTATIONS);
    vmag_atten_beam_crosspol_nointerferencefunc_bigpic.resize(geom.NSTATIONS);
    vmag_atten_beam_crosspol_func_bigpic.resize(geom.NSTATIONS);
    vtimediff_bigpic.resize(geom.NSTATIONS);
    vspectrum_bigpic.resize(geom.NSTATIONS);

    const int NDISTANCES_BIGPIC=5000;
    const double STEP=1.;

    // ------------------------------------------------------------------------
    // Main depth-dependent model outputs
    // ------------------------------------------------------------------------
    vector< vector<double> > vmag_atten_func;
    vector< vector<double> > vmag_parameter0;

    vector< vector<double> > vmag_atten;
    vector< vector<double> > vmag_atten_beam;
    vector< vector<double> > vmag_atten_beam_crosspol;
    vector< vector<double> > vmag_atten_beam_crosspol_nointerferencefunc;
    vector< vector<double> > vmag_atten_beam_crosspol_func;

    vector< vector<double> > vmag_func_noadjust;
    vector< vector<double> > vtimediff;
    vector< vector<double> > vsumlength;
    vector< vector<double> > vdeltan;
    vector< vector<double> > vnotflipped_alongpath;
    vector< vector<double> > vtheta1_alongpath;
    vector< vector<double> > vtheta2_alongpath;
    vector< vector<double> > vthetape1_alongpath;
    vector< vector<double> > vthetape2_alongpath;
    vector< vector<double> > vphipe1_alongpath;
    vector< vector<double> > vphipe2_alongpath;
    vector< vector<double> > vdepth_step;
    vector< vector<double> > vnotflipped;

    vector< vector<double> > vlengths;
    vector< vector<double> > vattenlengths;

    // Frequency-dependent spectrum / attenuation storage:
    //   vspectra[station][depth_index][freq_index]
    //   vattens [station][depth_index][freq_index]
    vector< vector< vector<double> > > vspectra;
    vector< vector< vector<double> > > vattens;

    // Resize outer station dimension
    model_data.vistep.resize(geom.NSTATIONS);
    vmag_parameter0.resize(geom.NSTATIONS);
    vmag_atten.resize(geom.NSTATIONS);
    vmag_atten_beam.resize(geom.NSTATIONS);
    vmag_atten_beam_crosspol.resize(geom.NSTATIONS);
    vmag_atten_beam_crosspol_nointerferencefunc.resize(geom.NSTATIONS);
    vmag_atten_beam_crosspol_func.resize(geom.NSTATIONS);
    vmag_func_noadjust.resize(geom.NSTATIONS);
    vtimediff.resize(geom.NSTATIONS);
    vnotflipped.resize(geom.NSTATIONS);
    vsumlength.resize(geom.NSTATIONS);
    vdeltan.resize(geom.NSTATIONS);
    vnotflipped_alongpath.resize(geom.NSTATIONS);
    vtheta1_alongpath.resize(geom.NSTATIONS);
    vtheta2_alongpath.resize(geom.NSTATIONS);
    vthetape1_alongpath.resize(geom.NSTATIONS);
    vthetape2_alongpath.resize(geom.NSTATIONS);
    vphipe1_alongpath.resize(geom.NSTATIONS);
    vphipe2_alongpath.resize(geom.NSTATIONS);
    vdepth_step.resize(geom.NSTATIONS);
    vlengths.resize(geom.NSTATIONS);
    vattenlengths.resize(geom.NSTATIONS);
    vspectra.resize(geom.NSTATIONS);
    vattens.resize(geom.NSTATIONS);

    vdistances_bigpic.resize(geom.NSTATIONS);
    vpseudodepths_bigpic.resize(geom.NSTATIONS);

    int jmin[geom.NSTATIONS];


};
