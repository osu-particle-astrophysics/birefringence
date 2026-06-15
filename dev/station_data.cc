#include <vector>
#include "TVector3.h"
#include "TGraph.h"

using namespace std;
// MACHTAY define some colors for printing text
#define RESET   "\x1b[0m"
#define BOLD    "\x1b[1m"
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define CYAN    "\x1b[36m"


//StationData::StationData(const Config& cfg, const vector<double>& nvec) {
//  // Read in the file
//  //
//}
//
void convert_to_TGraphs(const Config& cfg, const int station, StationData& model_data){
    // --------------------------------------------------------------------
		// Convert all accumulated per-depth vectors into ROOT TGraphs for plotting
		// and for later writing to output files.
		// --------------------------------------------------------------------

		model_data.gtxdepth_beam1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_beam1[station][0]);
		model_data.gtxdepth_beam2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_beam2[station][0]);

		model_data.grxdepth_beam1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_beam1[station][0]);
		model_data.grxdepth_beam2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_beam2[station][0]);

		model_data.grxdepth_atten[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_atten[station][0]);
		model_data.grxdepth_atten_beam[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_atten_beam[station][0]);
		model_data.grxdepth_atten_power[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_atten_power[station][0]);
		model_data.grxdepth_atten_beam_power[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_atten_beam_power[station][0]);

		model_data.gtxdepth_theta1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_theta1[station][0]);
		model_data.gtxdepth_theta2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_theta2[station][0]);

		model_data.grxdepth_theta1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_theta1[station][0]);
		model_data.grxdepth_theta2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_theta2[station][0]);

		model_data.grxdepthE_theta1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepthE_theta1[station][0]);
		model_data.grxdepthE_theta2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepthE_theta2[station][0]);

		model_data.gtxdepth_theta1_Sclock[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_theta1_Sclock[station][0]);
		model_data.gtxdepth_theta2_Sclock[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_theta2_Sclock[station][0]);

		model_data.grxdepth_theta1_Sclock[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_theta1_Sclock[station][0]);
		model_data.grxdepth_theta2_Sclock[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_theta2_Sclock[station][0]);

		model_data.gtxdepth_dispersion1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_dispersion1[station][0]);
		model_data.gtxdepth_dispersion2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_dispersion2[station][0]);

		// Note: these two graphs use model_data.vtxdepth_theta1/2 rather than model_data.vtxdepthE_theta1/2
		model_data.gtxdepthE_theta1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_theta1[station][0]);
		model_data.gtxdepthE_theta2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepth_theta2[station][0]);

		model_data.gtxdepthE_theta1_Sclock[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepthE_theta1_Sclock[station][0]);
		model_data.gtxdepthE_theta2_Sclock[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtxdepthE_theta2_Sclock[station][0]);

		model_data.grxdepthE_theta1_Sclock[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepthE_theta1_Sclock[station][0]);
		model_data.grxdepthE_theta2_Sclock[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepthE_theta2_Sclock[station][0]);

		model_data.gdotShats_tx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vdotShats_tx[station][0]);
		model_data.gdotEhats_tx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vdotEhats_tx[station][0]);
		model_data.gdotDhats_tx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vdotDhats_tx[station][0]);

		model_data.g_parameter0[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vmag_parameter0[station][0]);
		model_data.g_atten[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vmag_atten[station][0]);
		model_data.g_atten_beam[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vmag_atten_beam[station][0]);
		model_data.g_atten_beam_crosspol[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vmag_atten_beam_crosspol[station][0]);
		model_data.gfunc_noadjust[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vmag_func_noadjust[station][0]);
		model_data.g_atten_beam_crosspol_nointerferencefunc[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vmag_atten_beam_crosspol_nointerferencefunc[station][0]);
		model_data.g_atten_beam_crosspol_func[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vmag_atten_beam_crosspol_func[station][0]);
		model_data.g_sumphase[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vtimediff[station][0]);
		model_data.g_notflipped[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vnotflipped[station][0]);
		model_data.g_sumlength[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vsumlength[station][0]);

		// One frequency spectrum graph per pulser depth
		for (int idepth=0;idepth<model_data.vdepth[station].size();idepth++) {
		    model_data.g_spectra[station][idepth]=new TGraph(cfg.vfreqs.size(),&cfg.vfreqs[0],&model_data.vspectra[station][idepth][0]);
		}

		model_data.g_atten_power[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_atten_power[station][0]);
		model_data.g_atten_beam_power[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vrxdepth_atten_beam_power[station][0]);

		// Voltage / field / interference component graphs
		model_data.gV1_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV1_r1[station][0]);
		model_data.gV1_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV1_r2[station][0]);
		model_data.gV1squared_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV1squared_r1[station][0]);
		model_data.gV1squared_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV1squared_r2[station][0]);
		model_data.gpower_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vpower_r1[station][0]);
		model_data.gpower_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vpower_r2[station][0]);
		model_data.gvoltage_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vvoltage_r1[station][0]);
		model_data.gvoltage_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vvoltage_r2[station][0]);
		model_data.gfield_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vfield_r1[station][0]);
		model_data.gfield_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vfield_r2[station][0]);
		model_data.gV2_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV2_r1[station][0]);
		model_data.gV2_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV2_r2[station][0]);
		model_data.gV1V2_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV1V2_r1[station][0]);
		model_data.gV1V2_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV1V2_r2[station][0]);
		model_data.goppositeV1V2_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.voppositeV1V2_r1[station][0]);
		model_data.goppositeV1V2_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.voppositeV1V2_r2[station][0]);
		model_data.gV2squared_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV2squared_r1[station][0]);
		model_data.gV2squared_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vV2squared_r2[station][0]);
		model_data.genvelope_minus_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.venvelope_minus_r1[station][0]);
		model_data.genvelope_minus_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.venvelope_minus_r2[station][0]);
		model_data.genvelope_plus_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.venvelope_plus_r1[station][0]);
		model_data.genvelope_plus_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.venvelope_plus_r2[station][0]);
		model_data.gvenvelope_minus_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vvenvelope_minus_r1[station][0]);
		model_data.gvenvelope_minus_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vvenvelope_minus_r2[station][0]);
		model_data.gvenvelope_plus_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vvenvelope_plus_r1[station][0]);
		model_data.gvenvelope_plus_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vvenvelope_plus_r2[station][0]);
		model_data.gEenvelope_minus_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vEenvelope_minus_r1[station][0]);
		model_data.gEenvelope_minus_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vEenvelope_minus_r2[station][0]);
		model_data.gEenvelope_plus_r1[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vEenvelope_plus_r1[station][0]);
		model_data.gEenvelope_plus_r2[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vEenvelope_plus_r2[station][0]);
		model_data.gepsilon1_tx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vepsilon1_tx[station][0]);
		model_data.gepsilon2_tx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][1],&model_data.vepsilon2_tx[station][1]);
//    model_data.epsilon_differences[station]=new TGraph(model_data.vdepth[station.size(), &model_data.vdepth[station][0], &(model_data.vepsilon1_tx[station][0] - model_data.vepsilon1_rx[station][0]));

		// Polarization angle graphs
		model_data.gpolarization_Omega_rx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vpolarization_Omega_rx[station][0]);
		model_data.gpolarization_Psi_rx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vpolarization_Psi_rx[station][0]);

		model_data.gEpolarization_Omega_rx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vEpolarization_Omega_rx[station][0]);
		model_data.gEpolarization_Psi_rx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vEpolarization_Psi_rx[station][0]);

		// Build reversed-depth versions for plots that want positive pulser depth
		for (int j=0;j<model_data.vdepth[station].size();j++) {
		    model_data.vpolarization_reversedepth_Omega_rx[station].push_back(model_data.vpolarization_Omega_rx[station][model_data.vpolarization_Omega_rx[station].size()-1-j]);
		    model_data.vpolarization_reversedepth_Psi_rx[station].push_back(model_data.vpolarization_Psi_rx[station][model_data.vpolarization_Psi_rx[station].size()-1-j]);
		    model_data.vEpolarization_reversedepth_Omega_rx[station].push_back(model_data.vEpolarization_Omega_rx[station][model_data.vEpolarization_Omega_rx[station].size()-1-j]);
		    model_data.vEpolarization_reversedepth_Psi_rx[station].push_back(model_data.vEpolarization_Psi_rx[station][model_data.vEpolarization_Psi_rx[station].size()-1-j]);

		}

		model_data.gpolarization_reversedepth_Omega_rx[station]=new TGraph(model_data.vreversedepth[station].size(),&model_data.vreversedepth[station][0],&model_data.vpolarization_reversedepth_Omega_rx[station][0]);
		model_data.gpolarization_reversedepth_Psi_rx[station]=new TGraph(model_data.vreversedepth[station].size(),&model_data.vreversedepth[station][0],&model_data.vpolarization_reversedepth_Psi_rx[station][0]);

		model_data.gEpolarization_reversedepth_Omega_rx[station]=new TGraph(model_data.vreversedepth[station].size(),&model_data.vreversedepth[station][0],&model_data.vEpolarization_reversedepth_Omega_rx[station][0]);
		model_data.gEpolarization_reversedepth_Psi_rx[station]=new TGraph(model_data.vreversedepth[station].size(),&model_data.vreversedepth[station][0],&model_data.vEpolarization_reversedepth_Psi_rx[station][0]);

		model_data.gdiffepsilon_tx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vdiffepsilon_tx[station][0]);

		model_data.gepsilon1_rx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vepsilon1_rx[station][0]);
		model_data.gepsilon2_rx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][1],&model_data.vepsilon2_rx[station][1]);
		model_data.gdiffepsilon_rx[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vdiffepsilon_rx[station][0]);

		// Along-path diagnostic graphs for the special pulser depth
		model_data.g_deltan[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vdeltan[station][0]);
		model_data.g_notflipped_alongpath[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vnotflipped_alongpath[station][0]);
		model_data.g_theta1_alongpath[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vtheta1_alongpath[station][0]);
		model_data.g_theta2_alongpath[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vtheta2_alongpath[station][0]);
		model_data.g_thetape1_alongpath[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vthetape1_alongpath[station][0]);
		model_data.g_thetape2_alongpath[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vthetape2_alongpath[station][0]);
		model_data.g_thetape1_phipe1_alongpath[station]=new TGraph(model_data.vphipe1_alongpath[station].size(),&model_data.vphipe1_alongpath[station][0],&model_data.vthetape1_alongpath[station][0]);
		model_data.g_thetape2_phipe2_alongpath[station]=new TGraph(model_data.vphipe2_alongpath[station].size(),&model_data.vphipe2_alongpath[station][0],&model_data.vthetape2_alongpath[station][0]);
		model_data.g_phipe1_alongpath[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vphipe1_alongpath[station][0]);
		model_data.g_phipe2_alongpath[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vphipe2_alongpath[station][0]);
		model_data.g_deltan_pulserdepth[station]=new TGraph(model_data.vdepth_step[station].size(),&model_data.vdepth_step[station][0],&model_data.vdeltan[station][0]);

		model_data.g_depth_istep[station]=new TGraph(model_data.vdepth_step[station].size(),&model_data.vdepth_step[station][0],&model_data.vistep[station][0]);

		model_data.g_attenlengths[station]=new TGraph(model_data.vlengths[station].size(),&model_data.vlengths[station][0],&model_data.vattenlengths[station][0]);
		model_data.g_receive_launch[station]=new TGraph(model_data.vreceiveangle[station].size(),&model_data.vreceiveangle[station][0],&model_data.vlaunchangle[station][0]);
		model_data.g_receive[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vreceiveangle[station][0]);
		model_data.g_launch[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.vlaunchangle[station][0]);
		model_data.g_output6[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.voutput6[station][0]);
		model_data.g_output7[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.voutput7[station][0]);
		model_data.g_output8[station]=new TGraph(model_data.vdepth[station].size(),&model_data.vdepth[station][0],&model_data.voutput8[station][0]);

}


void plot_epsilons_tx_rx(const Geometry& geom, StationData& model_data) {
    // Side-by-side plot of transmitted (left) and received (right) epsilon_1
    string sdir_e = "output_plots/";
    gSystem->mkdir(sdir_e.c_str(), true);

    TCanvas *c = new TCanvas("c_tx_rx", "c_tx_rx", 1600, 800);
    c->Divide(2, 1);

    // --- Left panel: epsilon_1^T (TX) ---
    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    TH2D *hL = new TH2D("hL", "", 100, 600., 1600., 100, 0., 91.0);
    titles(hL, "", "Pulser height (m)", "#epsilon_{ 1}^{ T} (degrees)");
    hL->GetXaxis()->SetNdivisions(504);
    hL->GetYaxis()->SetNdivisions(504);
    hL->Draw();

    for (int ist = geom.minstation; ist <= geom.maxstation; ist++) {
        TGraph *g = model_data.gepsilon1_tx[ist];
        int n = g->GetN();
        TGraph *gplot = new TGraph(n);
        for (int i = 0; i < n; i++) {
            double x, y;
            g->GetPoint(i, x, y);
            gplot->SetPoint(i, -1.0 * x, y);
        }
        
        gplot->SetLineColor(icolors[ist]);
        gplot->SetLineWidth(2);
        gplot->SetLineStyle(kSolid);
        gplot->Draw("lsame");
    }

    // --- Right panel: epsilon_1^R (RX) with legend ---
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    TH2D *hR = new TH2D("hR", "", 100, 600., 1600., 100, -5., 95.);
    titles(hR, "", "Pulser height (m)", "#epsilon_{ 1}^{ R} (degrees)");
    hR->GetXaxis()->SetNdivisions(504);
    hR->GetYaxis()->SetNdivisions(504);
    hR->Draw();

    auto leg = new TLegend(0.55, 0.55, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    for (int ist = geom.minstation; ist <= geom.maxstation; ist++) {
        TGraph *g = model_data.gepsilon1_rx[ist];
        int n = g->GetN();
        TGraph *gplot = new TGraph(n);
        for (int i = 0; i < n; i++) {
            double x, y;
            g->GetPoint(i, x, y);
            gplot->SetPoint(i, -1.0 * x, y);
        }
        gplot->SetLineColor(icolors[ist]);
        gplot->SetLineWidth(2);
        gplot->SetLineStyle(kSolid);
        gplot->Draw("lsame");
        gplot->SetName(geom.snames[ist].c_str());
        leg->AddEntry(gplot, geom.snames[ist].c_str(), "l");
    }
    leg->Draw("same");

    c->Print((sdir_e + "fig6_epsilons_tx_rx.pdf").c_str());
    std::cout << "Printed TX/RX epsilons to: " << sdir_e + "fig6_epsilons_tx_rx.pdf" << std::endl;
    delete c;
}


void plot_epsilons_differences(const Geometry& geom, StationData& model_data) {
    // Plot of epsilon_1^T - epsilon_1^R (the difference) vs pulser height
    string sdir_e = "output_plots/";
    gSystem->mkdir(sdir_e.c_str(), true);

    TCanvas *c = new TCanvas("c_diff", "c_diff", 800, 800);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);

    TH2D *h = new TH2D("hDiff", "", 1000, 600., 1600., 100, -5., 95.);
    titles(h, "", "Pulser height (m)", "#epsilon_{ 1}^{ T} - #epsilon_{ 1}^{ R} (degrees)");
    h->GetXaxis()->SetNdivisions(504);
    h->GetYaxis()->SetNdivisions(504);
    h->Draw();

    auto leg = new TLegend(0.55, 0.55, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);

    for (int ist = geom.minstation; ist <= geom.maxstation; ist++) {
        TGraph *g_tx = model_data.gepsilon1_tx[ist];
        TGraph *g_rx = model_data.gepsilon1_rx[ist];
        int n = g_tx->GetN();
        model_data.epsilon_difference[ist] = new TGraph(n);

        for (int i = 0; i < n; i++) {
            double x_tx, y_tx, x_rx, y_rx;
            g_tx->GetPoint(i, x_tx, y_tx);
            g_rx->GetPoint(i, x_rx, y_rx);
            model_data.epsilon_difference[ist]->SetPoint(i, -1.0 * x_tx, y_tx - y_rx);
            //model_data.epsilon_difference[ist]->SetPoint(i, 1.0 * x_tx, y_tx - y_rx);
            if (x_tx == -600 && ist == 3 || x_tx == -1600 && ist == 3){
              cout << GREEN << "Depth: " << x_tx << endl;
              cout << GREEN << "Index (i): " << i << endl;
              cout << GREEN << "y_tx: " << y_tx << RESET << endl;
              cout << GREEN << "y_rx: " << y_rx << RESET << endl;
              cout << GREEN << "y_tx - y_rx: " << y_tx - y_rx << RESET << endl;
            }
        }

        model_data.epsilon_difference[ist]->SetLineColor(icolors[ist]);
        model_data.epsilon_difference[ist]->SetLineWidth(2);
        model_data.epsilon_difference[ist]->SetLineStyle(kSolid);
        model_data.epsilon_difference[ist]->Draw("lsame");
        model_data.epsilon_difference[ist]->SetName(geom.snames[ist].c_str());
        leg->AddEntry(model_data.epsilon_difference[ist], geom.snames[ist].c_str(), "l");
    }
    leg->Draw("same");

    c->Print((sdir_e + "fig6_epsilons_difference.pdf").c_str());
    std::cout << "Printed epsilon differences to: " << sdir_e + "fig6_epsilons_difference.pdf" << std::endl;
    delete c;
}


