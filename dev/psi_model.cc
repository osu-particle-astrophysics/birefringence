// psi_model.cc
//
// Definition of process_station(): the core per-station birefringence model,
// moved verbatim out of main() in cpol_main.cc.
//
// NOTE: this file is part of the unity build (M.cpol_main concatenates the
// .cc files into tmp.cc and compiles them as a single translation unit, with
// cpol_main.cc first). The ray-tracing routines used below (IceRayTracing,
// GetDirectRayPar, GetIceAttenuationLength, GetFull*RayPath, ...) are DEFINED
// in IceRayTracing.h, which has no include guard, so it is included exactly
// once -- by cpol_main.cc, earlier in the same translation unit. We therefore
// deliberately do NOT re-include IceRayTracing.h here.
#include "psi_model.hh"
#include "constants.hh"
#include "birefringence.hh"
#include "cpol_helpers.hh"

#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>

#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"

using namespace std;

// Text colors for diagnostic printouts (consistent with cpol_main.cc).
#ifndef RESET
#define RESET   "\x1b[0m"
#define BOLD    "\x1b[1m"
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define CYAN    "\x1b[36m"
#endif

void process_station(const Config& cfg, Geometry& geom, Ice_Profile& ice,
                     StationData& model_data, int i,
                     TVector3& p_e1, TVector3& p_e2, double& n_e1, double& n_e2,
                     bool SKIP_MIDDLE_STEPS, bool SKIP_RAYTRACE_PATH,
                     PsiTiming& timing) {

    using clock = std::chrono::system_clock;
    using sec   = std::chrono::duration<double>;

    // Distance-domain "big picture" scan parameters (were main() locals).
    const int NDISTANCES_BIGPIC=5000;
    const double STEP=1.;
    int jmin[geom.NSTATIONS];

    // Bind the timing accumulators to the names used by the moved loop body.
    double& uzair_time = timing.uzair_time;
    double& overhead   = timing.overhead;
    double& loop_time  = timing.loop_time;
    double& eval_time  = timing.eval_time;
    double& ray_time   = timing.ray_time;

    // ===================== begin moved per-station body =====================
	const auto overhead_before = clock::now(); 
        model_data.vmag_atten_beam[i].resize(model_data.vdepth[i].size());
        model_data.vmag_atten_beam_crosspol[i].resize(model_data.vdepth[i].size());
        model_data.vmag_atten_beam_crosspol_nointerferencefunc[i].resize(model_data.vdepth[i].size());
        model_data.vmag_atten_beam_crosspol_func[i].resize(model_data.vdepth[i].size());
        model_data.vspectra[i].resize(model_data.vdepth[i].size());
        model_data.vattens[i].resize(model_data.vdepth[i].size());
        model_data.g_spectra[i].resize(model_data.vdepth[i].size());
        model_data.vmag_atten[i].resize(model_data.vdepth[i].size());
        model_data.vmag_parameter0[i].resize(model_data.vdepth[i].size());
        model_data.vmag_func_noadjust[i].resize(model_data.vdepth[i].size());

        // Measured data vectors only exist for A1-A5, not ARIANNA
        if (i!=5) {
            model_data.vtotal_distances[i].resize(model_data.vdepth_data[i].size());
            model_data.vtotal_distances_err[i].resize(model_data.vdepth_data[i].size());
        }
    
        // Smallest possible source-receiver separation is the horizontal distance
        double mindistance=geom.horizontal_distances[i];
        jmin[i]=(int)(mindistance/STEP);

        // Build a distance scan and convert each distance into an equivalent
        // pseudo-depth for plotting/model-comparison purposes.
        for (int j=jmin[i]+1;j<NDISTANCES_BIGPIC;j++) {
            double thisdistance=STEP*(double)j;
            model_data.vdistances_bigpic[i].push_back(thisdistance);
            model_data.vpseudodepths_bigpic[i].push_back(geom.station_depths[i]-1.*sqrt(thisdistance*thisdistance-mindistance*mindistance));

        }

        // --------------------------------------------------------------------
        // Loop over all pulser depths for this station
        // --------------------------------------------------------------------
        const auto before = clock::now();
        for (int idepth=0;idepth<model_data.vdepth[i].size();idepth++) {

            // Reduced 2D geometry for ray tracing:
            // x = horizontal separation, z = depth
            double posstation[2];
            posstation[0]=sqrt(pow(geom.station_coords[i][0]-geom.pulser_coords[0],2)+pow(geom.station_coords[i][1]-geom.pulser_coords[1],2));
            posstation[1]=geom.station_depths[i];
            double pospulser[2];
            pospulser[0]=0.;
            pospulser[1]=model_data.vdepth[i][idepth];

            // Full 3D pulser position, useful for later geometric interpretation
            TVector3 pospulser3D;
            pospulser3D.SetX(geom.pulser_coords[0]);
            pospulser3D.SetY(geom.pulser_coords[1]);
            pospulser3D.SetZ(model_data.vdepth[i][idepth]);

            // 3D vector from pulser to station
            geom.pulsertostation[i][0]=(geom.station_coords[i][0]-geom.pulser_coords[0]);
            geom.pulsertostation[i][1]=(geom.station_coords[i][1]-geom.pulser_coords[1]);
            geom.pulsertostation[i][2]=geom.station_depths[i]-model_data.vdepth[i][idepth];

            // Legacy / diagnostic "special" geometry, based on station 0 and the first pulser depth
            double atten_special;
            double pospulser_special[2];
            pospulser_special[0]=0.;
            pospulser_special[1]=model_data.vdepth[0][0];
            double posstation_special[2];
            posstation_special[0]=sqrt(pow(geom.station_coords[0][0]-geom.pulser_coords[0],2)+pow(geom.station_coords[0][1]-geom.pulser_coords[1],2));
            posstation_special[1]=geom.station_depths[0];

            // Ray-tracing endpoints
            double x0=0;
            double z0=pospulser[1];
            double x1=posstation[0];
            double z1=posstation[1];

	    // MACHTAY look here:
	    const auto ray_before = clock::now(); 
            // Solve ray-tracing problem for this pulser depth and station
            double *getresults=IceRayTracing(x0,z0,x1,z1);

            double lvalue;
            vector<double> res;     // path x-coordinates along the ray
            vector<double> zs;      // path z-coordinates along the ray
            double launch_angle;
            double receive_angle;
            double *paramsd;
            double *paramsra;
            double *paramsre;

            // Store outputs from the ray tracer
            model_data.voutput6[i].push_back(getresults[6]);
            model_data.voutput7[i].push_back(getresults[7]);
            model_data.voutput8[i].push_back(getresults[8]);

	    const sec overhead_duration = clock::now() - overhead_before;
            overhead += static_cast<double>( overhead_duration.count() );
            // ------------------------------------------------------------
            // Determine which ray type exists and build the full ray path
            // ------------------------------------------------------------
	    // MACHTAY let's measure the time to do raytracing
//	    const auto ray_before = clock::now(); 
            if (getresults[6]!=-1000) {
                // Direct ray solution
                paramsd=GetDirectRayPar(z0,x1,z1);
                launch_angle=paramsd[1]/DEGRAD;
                receive_angle=paramsd[0]/DEGRAD;
                if (!SKIP_RAYTRACE_PATH)
                    GetFullDirectRayPath(z0, x1, z1, paramsd[3], res, zs);
            }
            else if (getresults[8]!=1000) {

                // Refracted ray solution
                paramsre=GetReflectedRayPar(z0, x1 ,z1);
                double LangR=paramsre[1];
                double RangR=paramsre[0];
                paramsra=GetRefractedRayPar(z0, x1 ,z1,LangR,RangR);
                launch_angle=paramsra[1]/DEGRAD;

                receive_angle=paramsra[0]/DEGRAD;
                if (!SKIP_RAYTRACE_PATH)
                    GetFullRefractedRayPath(z0, x1, z1, paramsra[7], paramsra[3], res, zs);

            }
            else if (getresults[7]!=1000) {

                // Reflected ray solution
                paramsre=GetReflectedRayPar(z0, x1 ,z1);
                double LangR=paramsre[1];
                double RangR=paramsre[0];
                launch_angle=paramsre[1]/DEGRAD;
                receive_angle=paramsre[0]/DEGRAD;
                if (!SKIP_RAYTRACE_PATH)
                    GetFullReflectedRayPath(z0, x1, z1, LangR, res, zs);
            }
	    const sec ray_duration = clock::now() - ray_before;
            ray_time += static_cast<double>( ray_duration.count() );


            // Running amplitude attenuation factor along the ray
            double atten=1.;

            // Frequency-dependent attenuation starts at unity for all frequencies
            model_data.vattens[i][idepth].clear();
            for (int ifreq=0;ifreq<NFREQ;ifreq++) {
                model_data.vattens[i][idepth].push_back(1.);
            }

            // Running path-length and birefringent phase accumulators
            double sumlength=0.;
            double sumphase=0.;

            // Proceed only if a usable ray type exists
            if (getresults[6]!=-1000 || getresults[8]!=-1000 || getresults[7]!=0) {

                // Horizontal unit vector toward the station.
                // Used to embed the 2D ray-tracing solution back into 3D.
                TVector3 yhat(geom.station_coords[i][0]-geom.pulser_coords[0],
                geom.station_coords[i][1]-geom.pulser_coords[1],
                0.);
                if (yhat.Mag()<HOWSMALLISTOOSMALL)
                cout << "yhat mag is " << yhat.Mag() << "\n";
                yhat.SetMag(1.);

                double angle_yhat=atan2(yhat[1],yhat[0]);


                // Declare variables needed by both fast and normal paths
			vector<double> nvec_thisstep(3,0.);
			TVector3 rhat_thisstep;
			double deltan_alongpath=0.;
			double deltantimeslength_alongpath=0.;
			TVector3 p_e1_previous=p_e1;
			TVector3 p_e2_previous=p_e2;
			double notflipped_previous=1.;
			double notflipped_atend=1.;
			double theta_e1_start=0.;
			double notflipped=1.;

			if (!SKIP_RAYTRACE_PATH) {
			    // Principal axis at start of path (needs zs[0] from full ray path)
			    nvec_thisstep[0]=ice.gn1->Eval(zs[0]);
			    nvec_thisstep[1]=ice.gn2->Eval(zs[0]);
			    nvec_thisstep[2]=ice.gn3->Eval(zs[0]);
			    rhat_thisstep[0]=-1.*(res[UZAIRSTEP]-res[0])*yhat[0];
			    rhat_thisstep[1]=-1.*(res[UZAIRSTEP]-res[0])*yhat[1];
			    rhat_thisstep[2]=-1.*(zs[UZAIRSTEP]-zs[0]);
			    if (rhat_thisstep.Mag()<1.E-8)
			        cout << "before calling getDeltaN at place 1, rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
			    deltan_alongpath=getDeltaN(cfg.BIAXIAL,nvec_thisstep,rhat_thisstep,geom.angle_iceflow,n_e1,n_e2,p_e1,p_e2, cfg.phi, cfg.theta, cfg.gamma);
			    if (p_e2.Mag()<HOWSMALLISTOOSMALL)
			        cout << "1, p_e2 is " << p_e2.Mag();
			    p_e1_previous=p_e1;
			    p_e2_previous=p_e2;
			}

			// --------------------------------------------------------
			// MACHTAY: SKIP_RAYTRACE_PATH fast path
			// Build rhat_thisstep for RX and TX directly from launch/receive
			// angles, bypassing the istep loop entirely.
			// --------------------------------------------------------
			if (SKIP_RAYTRACE_PATH) {
			    TVector3 plusz(0.,0.,1.);
			    TVector3 vrotate = plusz.Cross(yhat);
			    vrotate.SetMag(1.);

			    // --- RX endpoint (receive side) ---
			    TVector3 rhat_rx = plusz;
			    rhat_rx.Rotate(receive_angle, vrotate);
			    rhat_rx.SetMag(1.);

			    vector<double> nvec_rx(3);
			    nvec_rx[0]=ice.gn1->Eval(geom.station_depths[i]);
			    nvec_rx[1]=ice.gn2->Eval(geom.station_depths[i]);
			    nvec_rx[2]=ice.gn3->Eval(geom.station_depths[i]);
			    getDeltaN(cfg.BIAXIAL,nvec_rx,rhat_rx,geom.angle_iceflow,n_e1,n_e2,p_e1,p_e2, cfg.phi, cfg.theta, cfg.gamma);
          //if (i == 3 && idepth == 0 || i == 3 && idepth == 199){
          if (i == 3 && idepth == 20 || i == 3 && idepth == 120){
          cout << "station: " << i + 1 << " Depth: " << idepth * 5 + 600 << " DELTA N: " << getDeltaN(cfg.BIAXIAL,nvec_rx,rhat_rx,geom.angle_iceflow,n_e1,n_e2,p_e1,p_e2, cfg.phi, cfg.theta, cfg.gamma) << endl;
          cout << RED << "launch angle: " << GetDirectRayPar(z0,x1,z1)[1] << endl;
          cout << RED << "receive angle: " << GetDirectRayPar(z0,x1,z1)[0] << RESET << endl;
    cout << CYAN << "nvec_rx: " << nvec_rx[0] << " " << nvec_rx[1] << " " << nvec_rx[2] << RESET << endl;

          }
			    TVector3 epsilon_rx; epsilon_rx[0]=nvec_rx[0]*nvec_rx[0]; epsilon_rx[1]=nvec_rx[1]*nvec_rx[1]; epsilon_rx[2]=nvec_rx[2]*nvec_rx[2];
			    TVector3 E_e1_rx=rotateD(epsilon_rx,geom.angle_iceflow,p_e1);
			    TVector3 E_e2_rx=rotateD(epsilon_rx,geom.angle_iceflow,p_e2);

			    {
			        double theta_e1,theta_e2,thetaE_e1,thetaE_e2,theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock;
			        TVector3 Shat_e1,Shat_e2;
			        double E_e1_tc,E_e2_tc,E_e1_pc,E_e2_pc;
			        getManyAnglesontheClock(cfg.BIAXIAL,cfg.CROSSPOLANGLE_RX,
			            rhat_rx,p_e1,p_e2,E_e1_rx,E_e2_rx,
			            theta_e1,theta_e2,thetaE_e1,thetaE_e2,
			            theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock,
			            Shat_e1,Shat_e2,E_e1_tc,E_e2_tc,E_e1_pc,E_e2_pc);
              //cout << "CROSSPOLANGLE_RX: " << cfg.CROSSPOLANGLE_RX << endl;
			        // push zeros/defaults for unused quantities
			        model_data.vrxdepth_theta1[i].push_back(theta_e1*DEGRAD);
			        model_data.vrxdepth_theta2[i].push_back(theta_e2*DEGRAD);
			        model_data.vrxdepthE_theta1[i].push_back(theta_e1*DEGRAD);
			        model_data.vrxdepthE_theta2[i].push_back(theta_e2*DEGRAD);
			        model_data.vrxdepth_theta1_Sclock[i].push_back(theta_e1_Sclock*DEGRAD);
			        model_data.vrxdepth_theta2_Sclock[i].push_back(theta_e2_Sclock*DEGRAD);
			        model_data.vrxdepthE_theta1_Sclock[i].push_back(thetaE_e1_Sclock*DEGRAD);
			        model_data.vrxdepthE_theta2_Sclock[i].push_back(thetaE_e2_Sclock*DEGRAD);
			        double eps1_rx=0.,eps2_rx=0.;
			        thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,eps1_rx,eps2_rx);
              if (i == 3 && idepth == 20){
              cout << "idepth: " << idepth << endl;
              cout << "thetaE1_e1_Sclock: " << thetaE_e1_Sclock * DEGRAD << endl;
              cout << RED << "eps1_rx: " << eps1_rx * DEGRAD << endl;
              cout << RED << "eps2_rx: " << eps2_rx * DEGRAD << RESET << endl;
    cout << CYAN << "nvec_rx: " << nvec_rx[0] << " " << nvec_rx[1] << " " << nvec_rx[2] << RESET << endl;

              }
              if (i == 3 && idepth == 120){
              cout << "idepth: " << idepth << endl;
              cout << "thetaE1_e1_Sclock: " << thetaE_e1_Sclock * DEGRAD << endl;
              cout << RED << "eps1_rx: " << eps1_rx * DEGRAD << endl;
              cout << RED << "eps2_rx: " << eps2_rx * DEGRAD << RESET << endl;
    cout << CYAN << "nvec_rx: " << nvec_rx[0] << " " << nvec_rx[1] << " " << nvec_rx[2] << RESET << endl;

              }
			        model_data.vepsilon1_rx[i].push_back(eps1_rx*DEGRAD);
			        model_data.vepsilon2_rx[i].push_back(eps2_rx*DEGRAD);
			        model_data.vdiffepsilon_rx[i].push_back((eps2_rx-eps1_rx)*DEGRAD);
			        double beam_rx_e1=sin(Shat_e1.Theta());
			        double beam_rx_e2=sin(Shat_e2.Theta());
			        model_data.vrxdepth_beam1[i].push_back(beam_rx_e1);
			        model_data.vrxdepth_beam2[i].push_back(beam_rx_e2);
			        TVector3 plusz2(0.,0.,1.);
			        model_data.Pr2[i]=Shat_e1.Cross(plusz2); model_data.Pr2[i].SetMag(1.);
			        model_data.Pr1[i]=model_data.Pr2[i].Cross(Shat_e1); model_data.Pr1[i].SetMag(1.);
			    }

			    // --- TX endpoint (launch side) ---
			    TVector3 rhat_tx = plusz;
			    rhat_tx.Rotate(launch_angle, vrotate);
			    rhat_tx.SetMag(1.);

			    vector<double> nvec_tx(3);
			    nvec_tx[0]=ice.gn1->Eval(model_data.vdepth[i][idepth]);
			    nvec_tx[1]=ice.gn2->Eval(model_data.vdepth[i][idepth]);
			    nvec_tx[2]=ice.gn3->Eval(model_data.vdepth[i][idepth]);
			    getDeltaN(cfg.BIAXIAL,nvec_tx,rhat_tx,geom.angle_iceflow,n_e1,n_e2,p_e1,p_e2,cfg.phi,cfg.theta,cfg.gamma);
			    TVector3 epsilon_tx; epsilon_tx[0]=nvec_tx[0]*nvec_tx[0]; epsilon_tx[1]=nvec_tx[1]*nvec_tx[1]; epsilon_tx[2]=nvec_tx[2]*nvec_tx[2];
			    TVector3 E_e1_tx=rotateD(epsilon_tx,geom.angle_iceflow,p_e1);
			    TVector3 E_e2_tx=rotateD(epsilon_tx,geom.angle_iceflow,p_e2);

			    {
			        double theta_e1,theta_e2,thetaE_e1,thetaE_e2,theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock;
			        TVector3 Shat_e1,Shat_e2;
			        double E_e1_tc,E_e2_tc,E_e1_pc,E_e2_pc;
			        getManyAnglesontheClock(cfg.BIAXIAL,cfg.CROSSPOLANGLE_TX,
			            rhat_tx,p_e1,p_e2,E_e1_tx,E_e2_tx,
			            theta_e1,theta_e2,thetaE_e1,thetaE_e2,
			            theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock,
			            Shat_e1,Shat_e2,E_e1_tc,E_e2_tc,E_e1_pc,E_e2_pc);
			        if (i==5 && cfg.BIAXIAL==1) {
			            switchThem(thetaE_e1_Sclock,thetaE_e2_Sclock);
			            switchThem(theta_e1_Sclock,theta_e2_Sclock);
			            switchThem(theta_e1,theta_e2);
			            switchThem(thetaE_e1,thetaE_e2);
			        }
			        model_data.vangle_Shat_e1_khat[i].push_back(acos(Shat_e1.Dot(rhat_tx)/Shat_e1.Mag()/rhat_tx.Mag())*DEGRAD);
			        model_data.vangle_Shat_e2_khat[i].push_back(acos(Shat_e2.Dot(rhat_tx)/Shat_e2.Mag()/rhat_tx.Mag())*DEGRAD);
			        model_data.vtxdepth_beam1[i].push_back(sin(Shat_e1.Theta()));
			        model_data.vtxdepth_beam2[i].push_back(sin(Shat_e2.Theta()));
			        model_data.vtxdepth_theta1[i].push_back(theta_e1*DEGRAD);
			        model_data.vtxdepth_theta2[i].push_back(theta_e2*DEGRAD);
			        model_data.vtxdepth_theta1_Sclock[i].push_back(theta_e1_Sclock*DEGRAD);
			        model_data.vtxdepth_theta2_Sclock[i].push_back(theta_e2_Sclock*DEGRAD);
			        model_data.vtxdepthE_theta1[i].push_back(thetaE_e1*DEGRAD);
			        model_data.vtxdepthE_theta2[i].push_back(thetaE_e2*DEGRAD);
			        model_data.vtxdepthE_theta1_Sclock[i].push_back(thetaE_e1_Sclock*DEGRAD);
			        model_data.vtxdepthE_theta2_Sclock[i].push_back(thetaE_e2_Sclock*DEGRAD);
			        model_data.vtxdepth_dispersion1[i].push_back(acos(E_e1_tx.Dot(p_e1)/E_e1_tx.Mag()/p_e1.Mag())*DEGRAD);
			        model_data.vtxdepth_dispersion2[i].push_back(acos(E_e2_tx.Dot(p_e2)/E_e2_tx.Mag()/p_e2.Mag())*DEGRAD);
			        double eps1_tx=0.,eps2_tx=0.;
			        thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,eps1_tx,eps2_tx);
              if (i == 3 && idepth == 120) {
    cout << "p_e1: " << p_e1[0] << " " << p_e1[1] << " " << p_e1[2] << endl;
    cout << "p_e2: " << p_e2[0] << " " << p_e2[1] << " " << p_e2[2] << endl;
    cout << "E_e1_tx: " << E_e1_tx[0] << " " << E_e1_tx[1] << " " << E_e1_tx[2] << endl;
    cout << "E_e2_tx: " << E_e2_tx[0] << " " << E_e2_tx[1] << " " << E_e2_tx[2] << endl;
    cout << "rhat_tx: " << rhat_tx[0] << " " << rhat_tx[1] << " " << rhat_tx[2] << endl;
    cout << fixed << setprecision(8);
    cout << CYAN << "nvec_tx: " << nvec_tx[0] << " " << nvec_tx[1] << " " << nvec_tx[2] << endl;
    cout << CYAN << "nvec_rx: " << nvec_rx[0] << " " << nvec_rx[1] << " " << nvec_rx[2] << RESET << endl;
}
              if (i == 3 && idepth == 20){
              cout << "idepth: " << idepth << endl;
              cout << "thetaE1_e1_Sclock: " << thetaE_e1_Sclock * DEGRAD << endl;
              cout << RED << "eps1_tx: " << eps1_tx * DEGRAD << endl;
              cout << RED << "eps2_tx: " << eps2_tx * DEGRAD << RESET << endl;
              }
              if (i == 3 && idepth == 120){
              cout << "idepth: " << idepth << endl;
              cout << "thetaE1_e1_Sclock: " << thetaE_e1_Sclock * DEGRAD << endl;
              cout << RED << "eps1_tx: " << eps1_tx * DEGRAD << endl;
              cout << RED << "eps2_tx: " << eps2_tx * DEGRAD << RESET << endl;
              }
			        if (eps1_tx>PI/2.) eps1_tx-=PI;
			        model_data.vepsilon1_tx[i].push_back(eps1_tx*DEGRAD);
			        model_data.vepsilon2_tx[i].push_back(eps2_tx*DEGRAD);
			        model_data.vdiffepsilon_tx[i].push_back((eps2_tx-eps1_tx)*DEGRAD);
			        model_data.vdotShats_tx[i].push_back(Shat_e1.Dot(Shat_e2)/Shat_e1.Mag()/Shat_e2.Mag());
			        model_data.vdotEhats_tx[i].push_back(E_e1_tx.Dot(E_e2_tx)/E_e1_tx.Mag()/E_e2_tx.Mag());
			        model_data.vdotDhats_tx[i].push_back(p_e1.Dot(p_e2)/p_e1.Mag()/p_e2.Mag());
			        // atten/beam-weighted quantities: push zeros (not computed in fast mode)
			        model_data.vrxdepth_atten[i].push_back(0.);
			        model_data.vrxdepth_atten_beam[i].push_back(0.);
			        model_data.vrxdepth_atten_power[i].push_back(0.);
			        model_data.vrxdepth_atten_beam_power[i].push_back(0.);
			    }

			    // push zeros for all accumulated-path quantities
			    model_data.vtimediff[i].push_back(0.);
			    model_data.vnotflipped[i].push_back(1.);
			    model_data.vsumlength[i].push_back(0.);
			    // push zeros for voltage/power quantities (depend on atten which is 0)
			    for (auto* v : {&model_data.vV1_r1[i],&model_data.vV2_r1[i],&model_data.vV1_r2[i],&model_data.vV2_r2[i],
			                    &model_data.vE1_r1[i],&model_data.vE2_r1[i],&model_data.vE1_r2[i],&model_data.vE2_r2[i],
			                    &model_data.vV1_r1_lpda[i],&model_data.vV2_r1_lpda[i]}) v->push_back(0.);
			    for (auto* v : {&model_data.vV1squared_r1[i],&model_data.vV2squared_r1[i],&model_data.vV1V2_r1[i],
			                    &model_data.vV1squared_r2[i],&model_data.vV2squared_r2[i],&model_data.vV1V2_r2[i],
			                    &model_data.voppositeV1V2_r1[i],&model_data.voppositeV1V2_r2[i]}) v->push_back(0.);
			    for (auto* v : {&model_data.venvelope_minus_r1[i],&model_data.venvelope_plus_r1[i],
			                    &model_data.venvelope_minus_r2[i],&model_data.venvelope_plus_r2[i],
			                    &model_data.vSenvelope_minus_r1[i],&model_data.vSenvelope_plus_r1[i],
			                    &model_data.vSenvelope_minus_r2[i],&model_data.vSenvelope_plus_r2[i],
			                    &model_data.venvelope_minus_r1_lpda[i],&model_data.venvelope_plus_r1_lpda[i]}) v->push_back(0.);
			    for (auto* v : {&model_data.vvenvelope_minus_r1[i],&model_data.vvenvelope_plus_r1[i],
			                    &model_data.vvenvelope_minus_r2[i],&model_data.vvenvelope_plus_r2[i],
			                    &model_data.vEenvelope_minus_r1[i],&model_data.vEenvelope_plus_r1[i],
			                    &model_data.vEenvelope_minus_r2[i],&model_data.vEenvelope_plus_r2[i]}) v->push_back(0.);
			    model_data.vpower_r1[i].push_back(0.); model_data.vpower_r2[i].push_back(0.);
			    model_data.vpoynting_r1[i].push_back(0.); model_data.vpoynting_r2[i].push_back(0.);
			    model_data.vpower_r1_lpda[i].push_back(0.);
			    model_data.vvoltage_r1[i].push_back(0.); model_data.vvoltage_r2[i].push_back(0.);
			    model_data.vfield_r1[i].push_back(0.); model_data.vfield_r2[i].push_back(0.);
			    model_data.vvoltage_r1_lpda[i].push_back(0.);
			    model_data.vpolarization_Psi_rx[i].push_back(0.);
			    model_data.vpolarization_Omega_rx[i].push_back(0.);
			    model_data.vEpolarization_Psi_rx[i].push_back(0.);
			    model_data.vEpolarization_Omega_rx[i].push_back(0.);
			    for (int ifreq=0;ifreq<NFREQ;ifreq++) model_data.vspectra[i][idepth].push_back(0.);
			    // rhat_launch/receive from angles (same as normal path bottom)
			    model_data.rhat_launch[i] = plusz; model_data.rhat_launch[i].Rotate(launch_angle,vrotate); model_data.rhat_launch[i].SetMag(1.);
			    model_data.rhat_receive[i] = plusz; model_data.rhat_receive[i].Rotate(receive_angle,vrotate); model_data.rhat_receive[i].SetMag(1.);
			    model_data.vreceiveangle[i].push_back(model_data.rhat_receive[i].Theta()*DEGRAD);
			    model_data.vlaunchangle[i].push_back(model_data.rhat_launch[i].Theta()*DEGRAD);
			} else {
			// --------------------------------------------------------
			// Normal path: step along the ray
			// --------------------------------------------------------

			// --------------------------------------------------------
			// Step along the traced ray path in chunks of UZAIRSTEP
			// --------------------------------------------------------
			// MACHTAY time this loop

			const auto uzair_before = clock::now(); 
			for (int istep=UZAIRSTEP;istep<res.size();istep+=UZAIRSTEP) {

			    // MACHTAY: optionally skip middle steps -- only process RX (first) and TX (last)
			    bool is_rx_step = (istep == UZAIRSTEP);
			    bool is_tx_step = (abs((double)(istep-(int)res.size())) <= UZAIRSTEP);
			    if (SKIP_MIDDLE_STEPS && !is_rx_step && !is_tx_step) continue;

			    nvec_thisstep.resize(3);

			    // Refractive indices at this path point
			const auto eval_before = clock::now(); 
			    nvec_thisstep[0]=ice.gn1->Eval(zs[istep]);
			    nvec_thisstep[1]=ice.gn2->Eval(zs[istep]);
			    nvec_thisstep[2]=ice.gn3->Eval(zs[istep]);
			const sec eval_duration = clock::now() - eval_before;
			eval_time += static_cast<double>( eval_duration.count() );


			    if (istep>0) {

				// Local propagation direction reconstructed from the discrete path
				rhat_thisstep[0]=-1.*(res[istep]-res[istep-UZAIRSTEP])*yhat[0];
				rhat_thisstep[1]=-1.*(res[istep]-res[istep-UZAIRSTEP])*yhat[1];
				rhat_thisstep[2]=-1.*(zs[istep]-zs[istep-UZAIRSTEP]);

				// Debugging / illustration printout for A1 near -1000 m
				if (i==0 && idepth==model_data.g_idepth[i]->Eval(-1000.)) {
				    if (istep==50) {
					cout << "receive angle is " << rhat_thisstep.Theta()*DEGRAD << "\n";
					cout << "recieve vector is " << res[istep-50]-res[istep] << "\t" << zs[istep-50]-zs[istep] << "\n";
					TVector3 v3dtemp(-1.*(zs[istep-50.]-zs[istep])/500.,0.,(res[istep-50.]-res[istep])/500.);
					if (v3dtemp.Mag()<HOWSMALLISTOOSMALL)
					cout << "v3dtemp is " << v3dtemp.Mag() << "\n";
					v3dtemp.SetMag(0.075);
					cout << "polarization vector is " << v3dtemp[0] << "\t" << v3dtemp[2] << "\n";
				    }
				    if (abs((double)(istep-(int)res.size()))<=50) {

					cout << RED << "launch angle is " << rhat_thisstep.Theta()*DEGRAD << "\n";
					cout << RED << "launch vector is " << res[istep-50]-res[istep] << "\t" << zs[istep-50]-zs[istep] << RESET << "\n";
					TVector3 v3dtemp(-1.*(zs[istep-50.]-zs[istep])/500.,0.,(res[istep-50.]-res[istep])/500.);
					if (v3dtemp.Mag()<HOWSMALLISTOOSMALL)
					cout << "v3dtemp is " << v3dtemp.Mag() << "\n";
					v3dtemp.SetMag(0.075);
					cout << "polarization vector is " << v3dtemp[0] << "\t" << v3dtemp[2] << "\n";
				    }
				    if (istep%50==0)
				    cout << "\\draw[very thick] (" << res[istep-50]/1000. << ",0.," << zs[istep-50]/1000. << ") -- (" << res[istep]/1000. << ",0.," << zs[istep]/1000. << ");\n";


				}

				// Physical segment length of this path step
				double length=rhat_thisstep.Mag();

				if (rhat_thisstep.Mag()<HOWSMALLISTOOSMALL)
				cout << "rhat_thisstep mag is " << rhat_thisstep.Mag() << "\n";

				// Convert to unit direction for angular calculations
				rhat_thisstep.SetMag(1.);

				// Attenuation length at this depth for the chosen frequency
				double atten_length=GetIceAttenuationLength(zs[istep], freq/1.E9);

				// Amplitude attenuation accumulated along the path
				atten*=exp(-1.*length/atten_length);

				if (rhat_thisstep.Mag()<1.E-8){
				    cout << "before calling getDeltaN at place 2, rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
				}
				// Update local birefringence splitting and eigenvectors/eigenvalues
				deltan_alongpath=getDeltaN(cfg.BIAXIAL,nvec_thisstep,rhat_thisstep,geom.angle_iceflow,n_e1,n_e2,p_e1,p_e2,cfg.phi,cfg.theta,cfg.gamma);

				if (p_e2.Mag()<HOWSMALLISTOOSMALL)
				cout << "2, p_e2 is " << p_e2.Mag() << "\n";

				// Save previous eigenvectors in case continuity tracking is needed
				p_e1_previous=p_e1;
				p_e2_previous=p_e2;

				// Local dielectric tensor diagonal entries in the principal basis
				TVector3 epsilon_thisstep;

				epsilon_thisstep[0]=nvec_thisstep[0]*nvec_thisstep[0];
				epsilon_thisstep[1]=nvec_thisstep[1]*nvec_thisstep[1];
				epsilon_thisstep[2]=nvec_thisstep[2]*nvec_thisstep[2];


				// Convert D eigenvectors to electric-field directions
				TVector3 E_e1=rotateD(epsilon_thisstep,geom.angle_iceflow,p_e1);
				TVector3 E_e2=rotateD(epsilon_thisstep,geom.angle_iceflow,p_e2);

				// ----------------------------------------------------
				// At the final step (TX side in this path convention),
				// compute transmitter-side polarization / beam quantities
				// ----------------------------------------------------
				if (abs((double)(istep-(int)res.size()))<=UZAIRSTEP) {


				    if (i==5 && idepth==model_data.g_idepth[i]->Eval(-400.))
				    cout << "Station 5 launch angle at -400:" << rhat_thisstep.Theta()*DEGRAD << "\n";

				    double theta_e1,theta_e2;
				    double thetaE_e1,thetaE_e2;
				    double theta_e1_Sclock,theta_e2_Sclock;
				    double thetaE_e1_Sclock,thetaE_e2_Sclock;

				    TVector3 Shat_e1,Shat_e2;

				    double E_e1_thetacomponent,E_e2_thetacomponent;
				    double E_e1_phicomponent,E_e2_phicomponent;

				    // Decompose the eigenmodes into the chosen TX cross-pol frame
				    getManyAnglesontheClock(cfg.BIAXIAL,cfg.CROSSPOLANGLE_TX,
					rhat_thisstep,
					p_e1,p_e2,E_e1,E_e2,
					theta_e1,theta_e2,thetaE_e1,thetaE_e2,
					theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock,
					Shat_e1,Shat_e2,
					E_e1_thetacomponent,E_e2_thetacomponent,
					E_e1_phicomponent,E_e2_phicomponent);

				    // Store angle between each Poynting vector and propagation direction
				    model_data.vangle_Shat_e1_khat[i].push_back(acos(Shat_e1.Dot(rhat_thisstep)/Shat_e1.Mag()/rhat_thisstep.Mag())*DEGRAD);
				    model_data.vangle_Shat_e2_khat[i].push_back(acos(Shat_e2.Dot(rhat_thisstep)/Shat_e2.Mag()/rhat_thisstep.Mag())*DEGRAD);

				    // Special ordering swap for ARIANNA in the biaxial case
				    if (i==5 && cfg.BIAXIAL==1) {
					switchThem(thetaE_e1_Sclock,thetaE_e2_Sclock);
					switchThem(theta_e1_Sclock,theta_e2_Sclock);
					switchThem(theta_e1,theta_e2);
					switchThem(thetaE_e1,thetaE_e2);
				    }

				    // Verbose diagnostics for selected reference depths
				    if (i==0 && idepth==model_data.g_idepth[i]->Eval(-1000.) ||
				    i==5 && idepth==model_data.g_idepth[i]->Eval(-400.)) {
/*
					cout << "At Tx:\n";
					cout << "A" << i+1 << ", depth is " << model_data.vdepth[i][idepth] << "\n";

					cout << "angle of yhat is " << angle_yhat << "\n";

					cout << "p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
					cout << "mag of p_e1 is " << p_e1.Mag() << "\n";
					cout << "p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";
					cout << "mag of p_e2 is " << p_e2.Mag() << "\n";
					cout << "dot product is " << p_e1[0]*p_e2[0]+p_e1[1]*p_e2[1]+p_e1[2]*p_e2[2] << "\n";
					cout << "epsilon is " << epsilon_thisstep[0] << "\t" << epsilon_thisstep[1] << "\t" << epsilon_thisstep[2] << "\n";
					cout << "E_e1 is " << E_e1[0] << "\t" << E_e1[1] << "\t" << E_e1[2] << "\n";
					cout << "mag of E_e1 is " << E_e1.Mag() << "\n";
					cout << "E_e2 is " << E_e2[0] << "\t" << E_e2[1] << "\t" << E_e2[2] << "\n";
					cout << "mag of E_e2 is " << E_e2.Mag() << "\n";
					cout << "theta component of e1: " << E_e1_thetacomponent << "\t theta component of e2: " << E_e2_thetacomponent << "\n";
					cout << "phi component of e1: " << E_e1_phicomponent << "\t phi component of e2: " << E_e2_phicomponent << "\n";
					cout << "dot product is " << E_e1[0]*E_e2[0]+E_e1[1]*E_e2[1]+E_e1[2]*E_e2[2] << "\n";
					cout << "theta_e1, theta_e2 are " << theta_e1 << "\t" << theta_e2 << "\n";
					cout << "diff is " << (theta_e1-theta_e2)*DEGRAD << "\n";
					cout << "thetaE_e1, thetaE_e2 are " << thetaE_e1 << "\t" << thetaE_e2 << "\n";
					cout << "diff is " << (thetaE_e1-thetaE_e2)*DEGRAD << "\n";
					cout << "thetas on the Sclock are " << thetaE_e1_Sclock*DEGRAD << "\t" << thetaE_e2_Sclock*DEGRAD << "\n";
					cout << "diff is " << (thetaE_e1_Sclock-thetaE_e2_Sclock)*DEGRAD << "\n";
					cout << "depth is " << model_data.vdepth[i][idepth] << "\n";
					cout << "theta of rhat_thisstep is " << rhat_thisstep.Theta()*DEGRAD << "\n";
					cout << "rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
					cout << "Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
					cout << "Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";
*/
					TVector3 E_e1_temp=E_e1;
					TVector3 E_e2_temp=E_e2;
					TVector3 Shat_e1_temp=Shat_e1;
					TVector3 Shat_e2_temp=Shat_e2;
					TVector3 rhat_thisstep_temp=rhat_thisstep;

					// Total E field from adding the two eigenmode contributions
					TVector3 E_total_temp=E_e1_temp+E_e2_temp;

					TVector3 zaxis(0.,0.,1.);

					// Rotate so the line of sight lies in the plane of the page
					E_e1_temp.Rotate(-1.*angle_yhat,zaxis);
					E_e2_temp.Rotate(-1.*angle_yhat,zaxis);
					E_total_temp.Rotate(-1.*angle_yhat,zaxis);

					Shat_e1_temp.Rotate(-1.*angle_yhat,zaxis);
					Shat_e2_temp.Rotate(-1.*angle_yhat,zaxis);
					rhat_thisstep_temp.Rotate(-1.*angle_yhat,zaxis);

					// Rescale for prettier printing / drawing
					double scalefactor=0.2/E_total_temp.Mag();
					E_total_temp=scalefactor*E_total_temp;
					E_e1_temp=scalefactor*E_e1_temp;
					E_e2_temp=scalefactor*E_e2_temp;

					Shat_e1_temp=scalefactor*Shat_e1_temp;
					Shat_e2_temp=scalefactor*Shat_e2_temp;
					rhat_thisstep_temp=scalefactor*rhat_thisstep_temp;

          /*
					cout << "These are rotated so the line of sight from pulser to station is in the plane of the page.\n";
					cout << "E_e1_temp is " << E_e1_temp[0] << "\t" << E_e1_temp[1] << "\t" << E_e1_temp[2] << "\n";
					cout << "mag of E_e1_temp is " << E_e1_temp.Mag() << "\n";
					cout << "E_e2_temp is " << E_e2_temp[0] << "\t" << E_e2_temp[1] << "\t" << E_e2_temp[2] << "\n";
					cout << "mag of E_e2_temp is " << E_e2_temp.Mag() << "\n";

					cout << "E_total_temp is " << E_total_temp[0] << "\t" << E_total_temp[1] << "\t" << E_total_temp[2] << "\n";
					cout << "mag of E_total_temp is " << E_total_temp.Mag() << "\t" << 1/sqrt(E_total_temp.Mag()) << "\n";
					cout << "polarization angle is " << DEGRAD*atan(E_total_temp[1]/sqrt(E_total_temp[0]*E_total_temp[0]+E_total_temp[2]*E_total_temp[2])) << "\n";

					cout << "Shat_e1_temp is " << Shat_e1_temp[0] << "\t" << Shat_e1_temp[1] << "\t" << Shat_e1_temp[2] << "\n";
					cout << "Shat_e2_temp is " << Shat_e2_temp[0] << "\t" << Shat_e2_temp[1] << "\t" << Shat_e2_temp[2] << "\n";
					cout << "rhat_thisstep_temp is " << rhat_thisstep_temp[0] << "\t" << rhat_thisstep_temp[1] << "\t" << rhat_thisstep_temp[2] << "\n";

					TVector3 D_e1_temp=cos(theta_e1)*p_e1;
					TVector3 D_e2_temp=cos(theta_e2)*p_e2;

					D_e1_temp.Rotate(-1.*angle_yhat,zaxis);
					D_e2_temp.Rotate(-1.*angle_yhat,zaxis);

					cout << "mag of p_e1 is " << p_e1.Mag() << "\n";
					cout << "mag of p_e2 is " << p_e2.Mag() << "\n";

					cout << "D_e1_temp is " << D_e1_temp[0] << "\t" << D_e1_temp[1] << "\t" << D_e1_temp[2] << "\n";
					cout << "D_e2_temp is " << D_e2_temp[0] << "\t" << D_e2_temp[1] << "\t" << D_e2_temp[2] << "\n";
					cout << "mag of D_e1_temp is " << D_e1_temp.Mag() << "\n";
					cout << "mag of D_e2_temp is " << D_e2_temp.Mag() << "\n";
          */
				    }

				    // Store orthogonality / alignment diagnostics
				    model_data.vdotShats_tx[i].push_back(Shat_e1.Dot(Shat_e2)/Shat_e1.Mag()/Shat_e2.Mag());
				    model_data.vdotEhats_tx[i].push_back(E_e1.Dot(E_e2)/E_e1.Mag()/E_e2.Mag());
				    model_data.vdotDhats_tx[i].push_back(p_e1.Dot(p_e2)/p_e1.Mag()/p_e2.Mag());

				    if (i==0 && idepth==model_data.g_idepth[i]->Eval(-1000.) ||
				    i==5 && idepth==model_data.g_idepth[i]->Eval(-400.)) {
					//cout << "Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
					//cout << "Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";
					//cout << "theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";
				    }

				    // Beam factors at the transmitter from the Poynting-vector polar angle
				    double beam_tx_e1=sin(Shat_e1.Theta());
				    double beam_tx_e2=sin(Shat_e2.Theta());

				    model_data.vtxdepth_beam1[i].push_back(beam_tx_e1);
				    model_data.vtxdepth_beam2[i].push_back(beam_tx_e2);

				    // Store TX angular diagnostics
				    model_data.vtxdepth_theta1[i].push_back(theta_e1*DEGRAD);
				    model_data.vtxdepth_theta2[i].push_back(theta_e2*DEGRAD);

				    model_data.vtxdepth_theta1_Sclock[i].push_back(theta_e1_Sclock*DEGRAD);
				    model_data.vtxdepth_theta2_Sclock[i].push_back(theta_e2_Sclock*DEGRAD);

				    model_data.vtxdepthE_theta1[i].push_back(thetaE_e1*DEGRAD);
				    model_data.vtxdepthE_theta2[i].push_back(thetaE_e2*DEGRAD);


				    model_data.vtxdepthE_theta1_Sclock[i].push_back(thetaE_e1_Sclock*DEGRAD);
				    model_data.vtxdepthE_theta2_Sclock[i].push_back(thetaE_e2_Sclock*DEGRAD);

				    // Angle between E and D for each eigenmode
				    model_data.vtxdepth_dispersion1[i].push_back(acos(E_e1.Dot(p_e1)/E_e1.Mag()/p_e1.Mag())*DEGRAD);
				    model_data.vtxdepth_dispersion2[i].push_back(acos(E_e2.Dot(p_e2)/E_e2.Mag()/p_e2.Mag())*DEGRAD);

				    // Convert S-clock angles into epsilon parameters
				    double epsilon1_tx=0.;
				    double epsilon2_tx=0.;
				    thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,
				    epsilon1_tx,epsilon2_tx);

				    if (epsilon2_tx>PI/2.){
					epsilon2_tx-=PI;
				    }

				    model_data.vepsilon1_tx[i].push_back(epsilon1_tx*DEGRAD);
				    model_data.vepsilon2_tx[i].push_back(epsilon2_tx*DEGRAD);

				    model_data.vdiffepsilon_tx[i].push_back((epsilon2_tx-epsilon1_tx)*DEGRAD);

				    // Store attenuation-only and attenuation×beam estimates at the receiver
				    model_data.vrxdepth_atten[i].push_back(VOLTAGENORM*atten);
				    model_data.vrxdepth_atten_beam[i].push_back(VOLTAGENORM*model_data.vrxdepth_beam1[i][idepth]*model_data.vtxdepth_beam1[i][idepth]*atten);
				    model_data.vrxdepth_atten_power[i].push_back(VOLTAGENORM*VOLTAGENORM*atten*atten);
				    model_data.vrxdepth_atten_beam_power[i].push_back(VOLTAGENORM*VOLTAGENORM*model_data.vrxdepth_beam1[i][idepth]*model_data.vrxdepth_beam1[i][idepth]*model_data.vtxdepth_beam1[i][idepth]*model_data.vtxdepth_beam1[i][idepth]*atten*atten);
				    
				    notflipped_atend=notflipped;

				}

				// ----------------------------------------------------
				// At the first step (RX side in this convention),
				// compute receiver-side polarization / beam quantities
				// ----------------------------------------------------
				if (istep==UZAIRSTEP) {

				    double theta_e1,theta_e2;
				    double thetaE_e1,thetaE_e2;
				    double theta_e1_Sclock,theta_e2_Sclock;
				    double thetaE_e1_Sclock,thetaE_e2_Sclock;

				    TVector3 Shat_e1,Shat_e2;

				    double E_e1_thetacomponent,E_e2_thetacomponent;
				    double E_e1_phicomponent,E_e2_phicomponent;

				    // Decompose eigenmodes into the chosen RX cross-pol basis
				    getManyAnglesontheClock(cfg.BIAXIAL,cfg.CROSSPOLANGLE_RX,
					rhat_thisstep,
					p_e1,p_e2,E_e1,E_e2,
					theta_e1,theta_e2,thetaE_e1,thetaE_e2,
					theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock,
					Shat_e1,Shat_e2,
					E_e1_thetacomponent,E_e2_thetacomponent,
					E_e1_phicomponent,E_e2_phicomponent);

				    // Reference angle used to determine whether an eigenvector "flips"
				    theta_e1_start=theta_e1;

				    // Debug diagnostics at selected reference depths
				    if (i==0 && idepth==model_data.g_idepth[i]->Eval(-1000.) ||
				    i==5 && idepth==model_data.g_idepth[i]->Eval(-400.)) {
              /*
					cout << "At Rx:\n";
					cout << "A" << i+1 << ", depth is " << model_data.vdepth[i][idepth] << "\n";
					cout << "thetas on the Sclock are " << thetaE_e1_Sclock*DEGRAD << "\t" << thetaE_e2_Sclock*DEGRAD << "\n";
					cout << "depth is " << model_data.vdepth[i][idepth] << "\n";

					cout << "theta of rhat_thisstep is " << rhat_thisstep.Theta()*DEGRAD << "\n";
					cout << "rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
					cout << "This points from pulser to station: " << geom.station_coords[i][0]-geom.pulser_coords[0] << "\t" << geom.station_coords[i][1]-geom.pulser_coords[1] << "\t" << geom.station_depths[i]-model_data.vdepth[i][idepth] << "\n";

					cout << "Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
					cout << "Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";
					cout << "theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";
          */
					TVector3 E_e1_temp=E_e1;

					TVector3 E_e2_temp=E_e2;

					TVector3 Shat_e1_temp=Shat_e1;
					TVector3 Shat_e2_temp=Shat_e2;
					TVector3 rhat_thisstep_temp=rhat_thisstep;

					TVector3 E_total_temp=E_e1_temp+E_e2_temp;

					TVector3 zaxis(0.,0.,1.);

					E_e1_temp.Rotate(-1.*angle_yhat,zaxis);
					E_e2_temp.Rotate(-1.*angle_yhat,zaxis);
					E_total_temp.Rotate(-1.*angle_yhat,zaxis);

					Shat_e1_temp.Rotate(-1.*angle_yhat,zaxis);
					Shat_e2_temp.Rotate(-1.*angle_yhat,zaxis);
					rhat_thisstep_temp.Rotate(-1.*angle_yhat,zaxis);

					double scalefactor=0.2/E_total_temp.Mag();
					E_total_temp=scalefactor*E_total_temp;
					E_e1_temp=scalefactor*E_e1_temp;
					E_e2_temp=scalefactor*E_e2_temp;

					Shat_e1_temp=scalefactor*Shat_e1_temp;
					Shat_e2_temp=scalefactor*Shat_e2_temp;
					rhat_thisstep_temp=scalefactor*rhat_thisstep_temp;

					cout << "E_e1_temp is " << E_e1_temp[0] << "\t" << E_e1_temp[1] << "\t" << E_e1_temp[2] << "\n";
					cout << "mag of E_e1_temp is " << E_e1_temp.Mag() << "\n";
					cout << "E_e2_temp is " << E_e2_temp[0] << "\t" << E_e2_temp[1] << "\t" << E_e2_temp[2] << "\n";
					cout << "mag of E_e2_temp is " << E_e2_temp.Mag() << "\n";
					cout << "E_total_temp is " << E_total_temp[0] << "\t" << E_total_temp[1] << "\t" << E_total_temp[2] << "\n";
					cout << "theta component of e1: " << E_e1_thetacomponent << "\t theta component of e2: " << E_e2_thetacomponent << "\n";
					cout << "phi component of e1: " << E_e1_phicomponent << "\t phi component of e2: " << E_e2_phicomponent << "\n";

					cout << "Shat_e1_temp is " << Shat_e1_temp[0] << "\t" << Shat_e1_temp[1] << "\t" << Shat_e1_temp[2] << "\n";
					cout << "Shat_e2_temp is " << Shat_e2_temp[0] << "\t" << Shat_e2_temp[1] << "\t" << Shat_e2_temp[2] << "\n";
					cout << "rhat_thisstep_temp is " << rhat_thisstep_temp[0] << "\t" << rhat_thisstep_temp[1] << "\t" << rhat_thisstep_temp[2] << "\n";

					cout << "mag of E_total_temp is " << E_total_temp.Mag() << "\t" << 1/sqrt(E_total_temp.Mag()) << "\n";
					cout << "polarization angle is " << DEGRAD*atan(E_total_temp[1]/sqrt(E_total_temp[0]*E_total_temp[0]+E_total_temp[2]*E_total_temp[2])) << "\n";

					TVector3 D_e1_temp=cos(theta_e1)*p_e1;
					TVector3 D_e2_temp=cos(theta_e2)*p_e2;

					cout << "mag of p_e1 is " << p_e1.Mag() << "\n";
					cout << "mag of p_e2 is " << p_e2.Mag() << "\n";

					cout << "D_e1_temp is " << D_e1_temp[0] << "\t" << D_e1_temp[1] << "\t" << D_e1_temp[2] << "\n";
					cout << "D_e2_temp is " << D_e2_temp[0] << "\t" << D_e2_temp[1] << "\t" << D_e2_temp[2] << "\n";
					cout << "mag of D_e1_temp is " << D_e1_temp.Mag() << "\n";
					cout << "mag of D_e2_temp is " << D_e2_temp.Mag() << "\n";
				    }

				    // Store RX-side angular diagnostics
				    model_data.vrxdepth_theta1[i].push_back(theta_e1*DEGRAD);
				    model_data.vrxdepth_theta2[i].push_back(theta_e2*DEGRAD);

				    // Note: these two lines store theta_e1/theta_e2 rather than thetaE_e1/thetaE_e2.
				    // That may be intentional or may deserve a later check.
				    model_data.vrxdepthE_theta1[i].push_back(theta_e1*DEGRAD);
				    model_data.vrxdepthE_theta2[i].push_back(theta_e2*DEGRAD);

				    model_data.vrxdepth_theta1_Sclock[i].push_back(theta_e1_Sclock*DEGRAD);
				    model_data.vrxdepth_theta2_Sclock[i].push_back(theta_e2_Sclock*DEGRAD);


				    model_data.vrxdepthE_theta1_Sclock[i].push_back(thetaE_e1_Sclock*DEGRAD);
				    model_data.vrxdepthE_theta2_Sclock[i].push_back(thetaE_e2_Sclock*DEGRAD);

				    // Convert RX S-clock angles into epsilon parameters
				    double epsilon1_rx=0.;
				    double epsilon2_rx=0.;

				    thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,
					epsilon1_rx,epsilon2_rx);

				    model_data.vepsilon1_rx[i].push_back(epsilon1_rx*DEGRAD);
				    model_data.vepsilon2_rx[i].push_back(epsilon2_rx*DEGRAD);
				    model_data.vdiffepsilon_rx[i].push_back((epsilon2_rx-epsilon1_rx)*DEGRAD);

				    // Beam factors at the receiver
				    double beam_rx_e1=sin(Shat_e1.Theta());
				    double beam_rx_e2=sin(Shat_e2.Theta());

				    model_data.vrxdepth_beam1[i].push_back(beam_rx_e1);
				    model_data.vrxdepth_beam2[i].push_back(beam_rx_e2);

				    // Build an orthonormal receiver polarization basis:
				    //   Pr2 perpendicular to Shat_e1 and +z
				    //   Pr1 perpendicular to both Pr2 and Shat_e1
				    TVector3 plusz(0.,0.,1.);
				    model_data.Pr2[i]=Shat_e1.Cross(plusz);
				    if (model_data.Pr2[i].Mag()<HOWSMALLISTOOSMALL){
					cout << "model_data.Pr2[i] is " << model_data.Pr2[i].Mag() << "\n";
				    }
				    model_data.Pr2[i].SetMag(1.);
				    model_data.Pr1[i]=model_data.Pr2[i].Cross(Shat_e1);
				    if (model_data.Pr1[i].Mag()<HOWSMALLISTOOSMALL){
					cout << "model_data.Pr1[i] is " << model_data.Pr1[i].Mag() << "\n";
				    }
				    model_data.Pr1[i].SetMag(1.);
				}

				// ----------------------------------------------------
				// Update the frequency-dependent attenuation spectrum
				// along this path segment
				// ----------------------------------------------------
				for (int ifreq=0;ifreq<NFREQ;ifreq++) {


				    double this_atten_length=GetIceAttenuationLength(zs[istep], cfg.vfreqs[ifreq]/1.E9);

				    // Power-like attenuation: two exponential factors
				    model_data.vattens[i][idepth][ifreq] = model_data.vattens[i][idepth][ifreq]*exp(-1.*length/this_atten_length)*exp(-1.*length/this_atten_length);
				}

				// Recompute angular decomposition in the "neutral" clock convention
				double theta_e1,theta_e2;
				double thetaE_e1,thetaE_e2;
				double theta_e1_Sclock,theta_e2_Sclock;
				double thetaE_e1_Sclock,thetaE_e2_Sclock;

				TVector3 Shat_e1,Shat_e2;

				double E_e1_thetacomponent,E_e2_thetacomponent;
				double E_e1_phicomponent,E_e2_phicomponent;
				getManyAnglesontheClock(cfg.BIAXIAL,0.,
				    rhat_thisstep,
				    p_e1,p_e2,E_e1,E_e2,
				    theta_e1,theta_e2,thetaE_e1,thetaE_e2,
				    theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock,
				    Shat_e1,Shat_e2,
				    E_e1_thetacomponent,E_e2_thetacomponent,
				    E_e1_phicomponent,E_e2_phicomponent);

				// Determine whether the eigenvector basis has effectively flipped
				// relative to the starting RX-side convention
				notflipped=Flipped(theta_e1, theta_e1_start);

				// Accumulate signed birefringence phase integral
				deltantimeslength_alongpath+=deltan_alongpath*length*notflipped;
				notflipped_previous=notflipped;

				// Accumulate geometric path length
				sumlength+=length;

				// Store detailed along-path diagnostics only for a special pulser depth
				if (idepth==(int)(model_data.g_idepth[i]->Eval(geom.depth_special))) {

				    model_data.vtheta1_alongpath[i].push_back(theta_e1*DEGRAD);
				    model_data.vtheta2_alongpath[i].push_back(theta_e2*DEGRAD);
				    model_data.vthetape1_alongpath[i].push_back(p_e1.Theta()*DEGRAD);
				    model_data.vthetape2_alongpath[i].push_back(p_e2.Theta()*DEGRAD);
				    model_data.vphipe1_alongpath[i].push_back(p_e1.Phi()*DEGRAD);
				    model_data.vphipe2_alongpath[i].push_back(p_e2.Phi()*DEGRAD);
				    model_data.vnotflipped_alongpath[i].push_back(notflipped);
				    model_data.vdeltan[i].push_back(deltan_alongpath);
				    model_data.vdepth_step[i].push_back(zs[istep]);
				    model_data.vistep[i].push_back((double)istep);
				    model_data.vlengths[i].push_back(sumlength);
				    model_data.vattenlengths[i].push_back(atten_length);
				}
			    }
			    else {                        
				// Fallback for the first step if needed:
				// use the previously stored launch-direction unit vector.
				rhat_thisstep=model_data.rhat_launch[i];
			    }
			}
			const sec uzair_duration = clock::now() - uzair_before;
			uzair_time += static_cast<double>( uzair_duration.count() );

			} // end else (normal ray path stepping)

			// ----------------------------------------------------------------
			// Convert accumulated birefringence phase integral into a phase.
			//
			// deltantimeslength_alongpath has units of (delta n) * length.
			// Multiplying by (pi / c) * f converts this into the phase used
			// in the two-mode interference terms below.
			// ----------------------------------------------------------------
			sumphase=deltantimeslength_alongpath*PI/TMath::C()*freq;

			// Receiver-side S-clock angles for the two eigenmodes
			double theta1_Sclock_atrx,theta2_Sclock_atrx;


			// ----------------------------------------------------------------
			// Build the contributions of the two birefringent eigenmodes
			// to two receiver polarization channels (r1 and r2).
			//
			// Interpretation:
			//   - vV*_r1 / vV*_r2 are voltage-like projected amplitudes
			//   - vE*_r1 / vE*_r2 are field-like projected amplitudes
			//
			// Each contribution includes:
			//   - attenuation along the path
			//   - TX projection into a given eigenmode
			//   - RX projection into channel r1 or r2
			//   - beam-pattern factors at TX and RX (for voltages)
			// ----------------------------------------------------------------

			// Mode 1 contribution into receiver channel r1
			theta1_Sclock_atrx=model_data.vrxdepthE_theta1_Sclock[i][idepth];
			theta2_Sclock_atrx=model_data.vrxdepthE_theta2_Sclock[i][idepth];


			model_data.vV1_r1[i].push_back(model_data.vrxdepth_atten[i][idepth]*cos(model_data.vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(theta1_Sclock_atrx/DEGRAD)*model_data.vtxdepth_beam1[i][idepth]*model_data.vrxdepth_beam1[i][idepth]);
			model_data.vE1_r1[i].push_back(model_data.vrxdepth_atten[i][idepth]*cos(model_data.vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(theta1_Sclock_atrx/DEGRAD)*model_data.vtxdepth_beam1[i][idepth]);

			// Mode 2 contribution into receiver channel r1
			model_data.vV2_r1[i].push_back(model_data.vrxdepth_atten[i][idepth]*cos(model_data.vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD)*model_data.vtxdepth_beam2[i][idepth]*model_data.vrxdepth_beam2[i][idepth]);
			model_data.vE2_r1[i].push_back(model_data.vrxdepth_atten[i][idepth]*cos(model_data.vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD)*model_data.vtxdepth_beam2[i][idepth]);

			// LPDA-like projection: remove the component parallel to the beam axis
			model_data.vV1_r1_lpda[i].push_back(model_data.vV1_r1[i][idepth]*sqrt(1.-model_data.vrxdepth_beam1[i][idepth]*model_data.vrxdepth_beam1[i][idepth]));
			model_data.vV2_r1_lpda[i].push_back(model_data.vV2_r1[i][idepth]*sqrt(1.-model_data.vrxdepth_beam2[i][idepth]*model_data.vrxdepth_beam2[i][idepth]));

			// Store quadratic combinations for interference calculations
			model_data.vV1squared_r1[i].push_back(model_data.vV1_r1[i][idepth]*model_data.vV1_r1[i][idepth]);
			model_data.vV2squared_r1[i].push_back(model_data.vV2_r1[i][idepth]*model_data.vV2_r1[i][idepth]);
			model_data.vV1V2_r1[i].push_back(model_data.vV1_r1[i][idepth]*model_data.vV2_r1[i][idepth]);

			// Mode 1 contribution into receiver channel r2
			model_data.vV1_r2[i].push_back(model_data.vrxdepth_atten[i][idepth]*cos(model_data.vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(theta1_Sclock_atrx/DEGRAD)*model_data.vtxdepth_beam1[i][idepth]*model_data.vrxdepth_beam1[i][idepth]);
			model_data.vE1_r2[i].push_back(model_data.vrxdepth_atten[i][idepth]*cos(model_data.vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(theta1_Sclock_atrx/DEGRAD)*model_data.vtxdepth_beam1[i][idepth]);

			// Mode 2 contribution into receiver channel r2
			model_data.vV2_r2[i].push_back(model_data.vrxdepth_atten[i][idepth]*cos(model_data.vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(theta2_Sclock_atrx/DEGRAD)*model_data.vtxdepth_beam2[i][idepth]*model_data.vrxdepth_beam2[i][idepth]);
			model_data.vE2_r2[i].push_back(model_data.vrxdepth_atten[i][idepth]*cos(model_data.vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(theta2_Sclock_atrx/DEGRAD)*model_data.vtxdepth_beam2[i][idepth]);

			// More quadratic combinations for the second receiver channel
			model_data.vV1squared_r2[i].push_back(model_data.vV1_r2[i][idepth]*model_data.vV1_r2[i][idepth]);
			model_data.vV2squared_r2[i].push_back(model_data.vV2_r2[i][idepth]*model_data.vV2_r2[i][idepth]);
			model_data.vV1V2_r2[i].push_back(model_data.vV1_r2[i][idepth]*model_data.vV2_r2[i][idepth]);

			// Debug print for a reference pulser depth
//			if (idepth==model_data.g_idepth[i]->Eval(-1000.)){
//			    cout << "station, depth, V1squared_r2, V2squared_r2, V1V2_r2 are " << i << "\t" << model_data.vV1squared_r2[i][idepth] << "\t" << model_data.vV2squared_r2[i][idepth] << "\t" << model_data.vV1V2_r2[i][idepth] << "\n";
//			}

			// Negative cross terms, convenient for expressions written as
			// A+B-2*sqrt(AB)*sin^2(phi) or equivalent forms
			model_data.voppositeV1V2_r2[i].push_back(-1.*model_data.vV1_r2[i][idepth]*model_data.vV2_r2[i][idepth]);
			model_data.voppositeV1V2_r1[i].push_back(-1.*model_data.vV1_r1[i][idepth]*model_data.vV2_r1[i][idepth]);


			// ----------------------------------------------------------------
			// Build envelope quantities:
			//   minus envelope = destructive combination
			//   plus envelope  = constructive combination
			//
			// Versions are stored both for voltage-like quantities and
			// field/Poynting-like quantities.
			// ----------------------------------------------------------------
			model_data.venvelope_minus_r1[i].push_back((model_data.vV1_r1[i][idepth]-model_data.vV2_r1[i][idepth])*(model_data.vV1_r1[i][idepth]-model_data.vV2_r1[i][idepth]));
			model_data.vSenvelope_minus_r1[i].push_back((model_data.vE1_r1[i][idepth]-model_data.vE2_r1[i][idepth])*(model_data.vE1_r1[i][idepth]-model_data.vE2_r1[i][idepth]));

			model_data.venvelope_minus_r1_lpda[i].push_back((model_data.vV1_r1_lpda[i][idepth]-model_data.vV2_r1_lpda[i][idepth])*(model_data.vV1_r1_lpda[i][idepth]-model_data.vV2_r1_lpda[i][idepth]));

			model_data.venvelope_minus_r2[i].push_back((model_data.vV1_r2[i][idepth]-model_data.vV2_r2[i][idepth])*(model_data.vV1_r2[i][idepth]-model_data.vV2_r2[i][idepth]));
			model_data.vSenvelope_minus_r2[i].push_back((model_data.vE1_r2[i][idepth]-model_data.vE2_r2[i][idepth])*(model_data.vE1_r2[i][idepth]-model_data.vE2_r2[i][idepth]));
			model_data.venvelope_plus_r1[i].push_back((model_data.vV2_r1[i][idepth]+model_data.vV1_r1[i][idepth])*(model_data.vV1_r1[i][idepth]+model_data.vV2_r1[i][idepth]));
			model_data.vSenvelope_plus_r1[i].push_back((model_data.vE2_r1[i][idepth]+model_data.vE1_r1[i][idepth])*(model_data.vE1_r1[i][idepth]+model_data.vE2_r1[i][idepth]));
			model_data.venvelope_plus_r1_lpda[i].push_back((model_data.vV2_r1_lpda[i][idepth]+model_data.vV1_r1_lpda[i][idepth])*(model_data.vV1_r1_lpda[i][idepth]+model_data.vV2_r1_lpda[i][idepth]));

			model_data.venvelope_plus_r2[i].push_back((model_data.vV2_r2[i][idepth]+model_data.vV1_r2[i][idepth])*(model_data.vV1_r2[i][idepth]+model_data.vV2_r2[i][idepth]));
			model_data.vSenvelope_plus_r2[i].push_back((model_data.vE2_r2[i][idepth]+model_data.vE1_r2[i][idepth])*(model_data.vE1_r2[i][idepth]+model_data.vE2_r2[i][idepth]));

			model_data.vvenvelope_minus_r1[i].push_back(sqrt((model_data.vV1_r1[i][idepth]-model_data.vV2_r1[i][idepth])*(model_data.vV1_r1[i][idepth]-model_data.vV2_r1[i][idepth])));
			model_data.vvenvelope_minus_r2[i].push_back(sqrt((model_data.vV1_r2[i][idepth]-model_data.vV2_r2[i][idepth])*(model_data.vV1_r2[i][idepth]-model_data.vV2_r2[i][idepth])));
			model_data.vvenvelope_plus_r1[i].push_back(sqrt((model_data.vV2_r1[i][idepth]+model_data.vV1_r1[i][idepth])*(model_data.vV1_r1[i][idepth]+model_data.vV2_r1[i][idepth])));
			model_data.vvenvelope_plus_r2[i].push_back(sqrt((model_data.vV2_r2[i][idepth]+model_data.vV1_r2[i][idepth])*(model_data.vV1_r2[i][idepth]+model_data.vV2_r2[i][idepth])));

			model_data.vEenvelope_minus_r1[i].push_back(sqrt((model_data.vE1_r1[i][idepth]-model_data.vE2_r1[i][idepth])*(model_data.vE1_r1[i][idepth]-model_data.vE2_r1[i][idepth])));
			model_data.vEenvelope_minus_r2[i].push_back(sqrt((model_data.vE1_r2[i][idepth]-model_data.vE2_r2[i][idepth])*(model_data.vE1_r2[i][idepth]-model_data.vE2_r2[i][idepth])));
			model_data.vEenvelope_plus_r1[i].push_back(sqrt((model_data.vE2_r1[i][idepth]+model_data.vE1_r1[i][idepth])*(model_data.vE1_r1[i][idepth]+model_data.vE2_r1[i][idepth])));
			model_data.vEenvelope_plus_r2[i].push_back(sqrt((model_data.vE2_r2[i][idepth]+model_data.vE1_r2[i][idepth])*(model_data.vE1_r2[i][idepth]+model_data.vE2_r2[i][idepth])));

			// ----------------------------------------------------------------
			// Build a frequency spectrum for this station and pulser depth.
			//
			// The phase scales linearly with frequency, so a spectrum is obtained
			// by re-evaluating the same interference expression at each frequency.
			// ----------------------------------------------------------------
			for (int ifreq=0;ifreq<NFREQ;ifreq++) {
			    double thisfreq=cfg.vfreqs[ifreq];
			    double thissumphase=sumphase*thisfreq/freq;

			    model_data.vspectra[i][idepth].push_back(model_data.vattens[i][idepth][ifreq]/model_data.vrxdepth_atten[i][idepth]*(model_data.venvelope_plus_r1[i][idepth]-4*model_data.vV1_r1[i][idepth]*model_data.vV2_r1[i][idepth]*sin(thissumphase)*sin(thissumphase)));
//			    if (i==5 && idepth==model_data.g_idepth[i]->Eval(-1000.))
//			    cout << "i, model_data.vattens[i][idepth], model_data.vrxdepth_atten[i][idepth], model_data.vspectra are " << i << "\t" << model_data.vattens[i][idepth][ifreq] << "\t" << model_data.vrxdepth_atten[i][idepth] << "\t" << model_data.vspectra[i][idepth][ifreq] << "\n";
			}

			// Time delay corresponding to the accumulated phase difference
			//   sumphase = pi * f * dt
			// so dt = sumphase / (pi*f)
			model_data.vtimediff[i].push_back(sumphase/(PI)*1./freq*1.E9);

			// Whether the eigenvector tracking ended in the same orientation sign
			model_data.vnotflipped[i].push_back(notflipped_atend);

			// ----------------------------------------------------------------
			// Final interference-modified power / Poynting-like quantities
			// for the two receiver channels.
			// ----------------------------------------------------------------
			model_data.vpower_r1[i].push_back(model_data.venvelope_plus_r1[i][idepth]-4*model_data.vV1_r1[i][idepth]*model_data.vV2_r1[i][idepth]*sin(sumphase)*sin(sumphase));
			model_data.vpoynting_r1[i].push_back(model_data.vSenvelope_plus_r1[i][idepth]-4*model_data.vE1_r1[i][idepth]*model_data.vE2_r1[i][idepth]*sin(sumphase)*sin(sumphase));

			model_data.vpower_r1_lpda[i].push_back(model_data.venvelope_plus_r1_lpda[i][idepth]-4*model_data.vV1_r1_lpda[i][idepth]*model_data.vV2_r1_lpda[i][idepth]*sin(sumphase)*sin(sumphase));
			model_data.vpower_r2[i].push_back(model_data.venvelope_plus_r2[i][idepth]-4*model_data.vV1_r2[i][idepth]*model_data.vV2_r2[i][idepth]*sin(sumphase)*sin(sumphase));
			model_data.vpoynting_r2[i].push_back(model_data.vSenvelope_plus_r2[i][idepth]-4*model_data.vE1_r2[i][idepth]*model_data.vE2_r2[i][idepth]*sin(sumphase)*sin(sumphase));

			// Convert power-like quantities into amplitude-like quantities
			model_data.vvoltage_r1[i].push_back(sqrt(model_data.vpower_r1[i][idepth]));
			model_data.vfield_r1[i].push_back(sqrt(model_data.vpoynting_r1[i][idepth]));
			model_data.vvoltage_r1_lpda[i].push_back(sqrt(model_data.vpower_r1_lpda[i][idepth]));
			model_data.vvoltage_r2[i].push_back(sqrt(model_data.vpower_r2[i][idepth]));
			model_data.vfield_r2[i].push_back(sqrt(model_data.vpoynting_r2[i][idepth]));

			// Detailed printout at selected reference depths
			if (i==0 && idepth==model_data.g_idepth[i]->Eval(-1000.) ||
			i==5 && idepth==model_data.g_idepth[i]->Eval(-400.)) {

			    cout << "A " << i+1 << ", depth is " << model_data.vdepth[i][idepth] << "\n";

          /*
			    cout << "before scalefactors:\n";
			    cout << "powers are " << model_data.vpower_r1[i][idepth] << "\t" << model_data.vpower_r2[i][idepth] << "\n";
			    cout << "voltages are " << model_data.vvoltage_r1[i][idepth] << "\t" << model_data.vvoltage_r2[i][idepth] << "\n";
			    cout << "lpda voltages are " << model_data.vvoltage_r1_lpda[i][idepth] << "\t" << model_data.vvoltage_r2[i][idepth] << "\n";

			    // Scale to a convenient display size
			    double scalefactor=0.2/sqrt(model_data.vvoltage_r1[i][idepth]*model_data.vvoltage_r1[i][idepth]+model_data.vvoltage_r2[i][idepth]*model_data.vvoltage_r2[i][idepth]);

			    // Common attenuation × beam prefactor
			    double prefactor=model_data.vrxdepth_atten[i][idepth]*model_data.vtxdepth_beam1[i][idepth]*model_data.vrxdepth_beam1[i][idepth];
			    cout << "voltage r1 is " << scalefactor*model_data.vvoltage_r1[i][idepth] << "\n";
			    cout << "x, z components of voltage r1 are " << model_data.vtxdepth_beam1[i][idepth]*scalefactor*model_data.vvoltage_r1[i][idepth] << "\t" << sqrt(1.-model_data.vtxdepth_beam1[i][idepth]*model_data.vtxdepth_beam1[i][idepth])*scalefactor*model_data.vvoltage_r1[i][idepth] << "\n";

			    cout << "voltage r2 is " << scalefactor*model_data.vvoltage_r2[i][idepth] << "\n";
			    cout << "polarization angle is " << DEGRAD*atan2(model_data.vvoltage_r2[i][idepth],model_data.vvoltage_r1[i][idepth]) << "\n";
			    cout << "thetas_Sclock_atrx are " << theta1_Sclock_atrx << "\t" << theta2_Sclock_atrx << "\n";
			    cout << "diff is " << (theta1_Sclock_atrx-theta2_Sclock_atrx) << "\n";
			    cout << "thetas_Sclock_attx are " << model_data.vtxdepthE_theta1_Sclock[i][idepth] << "\t" << model_data.vtxdepthE_theta2_Sclock[i][idepth] << "\n";
			    cout << "diff is " << (model_data.vtxdepthE_theta1_Sclock[i][idepth]-model_data.vtxdepthE_theta2_Sclock[i][idepth]) << "\n";

			    cout << "model_data.vV1_r1, model_data.vV2_r1 are " << model_data.vV1_r1[i][idepth]/prefactor << "\t" << model_data.vV2_r1[i][idepth]/prefactor << "\n";
			    cout << "model_data.vV1_r2, model_data.vV2_r2 are " << model_data.vV1_r2[i][idepth]/prefactor << "\t" << model_data.vV2_r2[i][idepth]/prefactor << "\n";

			    cout << "ray 1 at tx is " << scalefactor*model_data.vrxdepth_atten[i][idepth]*model_data.vtxdepth_beam1[i][idepth]*cos(model_data.vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD) << "\n";
			    cout << "ray 2 at tx is " << scalefactor*model_data.vrxdepth_atten[i][idepth]*model_data.vtxdepth_beam1[i][idepth]*cos(model_data.vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD) << "\n";
			    cout << "fraction of wavelength is " << sumphase/(2.*PI) << "\n";
			    cout << "model_data.vpower_r1+model_data.vpower_r2 \t" <<model_data.vpower_r1[i][idepth]+model_data.vpower_r2[i][idepth] << "\n";
			    cout << "prefactor is " << pow(model_data.vrxdepth_atten[i][idepth]*model_data.vtxdepth_beam1[i][idepth]*model_data.vrxdepth_beam1[i][idepth],2) << "\n";
			    cout << "three factors are " << model_data.vrxdepth_atten[i][idepth] << "\t" << model_data.vtxdepth_beam1[i][idepth] << "\t" << model_data.vrxdepth_beam1[i][idepth] << "\n";
          */


			}

			// ----------------------------------------------------------------
			// Derived polarization observables at the receiver.
			//
			// Psi   ~ arctangent of the channel ratio
			// Omega ~ complementary representation using arccos of normalized r1
			//
			// Same quantities are computed both for voltage-like and field-like
			// amplitudes.
			// ----------------------------------------------------------------
			model_data.vpolarization_Psi_rx[i].push_back(atan2(model_data.vvoltage_r2[i][idepth],model_data.vvoltage_r1[i][idepth])*DEGRAD);
			model_data.vpolarization_Omega_rx[i].push_back(acos(model_data.vvoltage_r1[i][idepth]/(sqrt(model_data.vvoltage_r1[i][idepth]*model_data.vvoltage_r1[i][idepth]+model_data.vvoltage_r2[i][idepth]*model_data.vvoltage_r2[i][idepth])))*DEGRAD);
			model_data.vEpolarization_Psi_rx[i].push_back(atan2(model_data.vfield_r2[i][idepth],model_data.vfield_r1[i][idepth])*DEGRAD);
			model_data.vEpolarization_Omega_rx[i].push_back(acos(model_data.vfield_r1[i][idepth]/(sqrt(model_data.vfield_r1[i][idepth]*model_data.vfield_r1[i][idepth]+model_data.vfield_r2[i][idepth]*model_data.vfield_r2[i][idepth])))*DEGRAD);

			// Store total geometric path length
			model_data.vsumlength[i].push_back(sumlength);

			// Refractive indices at the pulser depth itself
			vector<double> nvec_thisdepth;
			nvec_thisdepth.resize(3);
			nvec_thisdepth[0]=ice.gn1->Eval(model_data.vdepth[i][idepth]);
			nvec_thisdepth[1]=ice.gn2->Eval(model_data.vdepth[i][idepth]);
			nvec_thisdepth[2]=ice.gn3->Eval(model_data.vdepth[i][idepth]);

			// ----------------------------------------------------------------
			// Build launch and receive direction unit vectors either from
			// ray tracing or from simple straight-line geometry.
			// ----------------------------------------------------------------
			if (geom.DORAYTRACING) {

			    TVector3 plusz(0.,0.,1.);

			    // Rotation axis that brings +z into the vertical plane of the ray
			    TVector3 vrotate=plusz.Cross(yhat);

			    if (vrotate.Mag()<HOWSMALLISTOOSMALL){
				cout << "vrotate mag is " << vrotate.Mag() << "\n";
			    }

			    // Launch direction
			    vrotate.SetMag(1.);
			    model_data.rhat_launch[i]=plusz;
			    model_data.rhat_launch[i].Rotate(launch_angle,vrotate);

			    // Receive direction
			    model_data.rhat_receive[i]=plusz;
			    model_data.rhat_receive[i].Rotate(receive_angle,vrotate);

			}
			else {

			    // Straight-line approximation if ray tracing is disabled
			    model_data.rhat_launch[i].SetX(geom.station_coords[i][0]-geom.pulser_coords[0]);
			    model_data.rhat_launch[i].SetY(geom.station_coords[i][1]-geom.pulser_coords[1]);
			    model_data.rhat_launch[i].SetZ(geom.station_depths[i]-model_data.vdepth[i][idepth]);

			    model_data.rhat_receive[i].SetX(model_data.rhat_launch[i][0]);
			    model_data.rhat_receive[i].SetY(model_data.rhat_launch[i][1]);
			    model_data.rhat_receive[i].SetZ(model_data.rhat_launch[i][2]);
			}

			// Normalize and sanity check
			if (model_data.rhat_launch[i].Mag()<HOWSMALLISTOOSMALL){
			    cout << "rhat[i] mag is " << model_data.rhat_launch[i].Mag() << "\n";
			}
			model_data.rhat_launch[i].SetMag(1.);

			if (model_data.rhat_receive[i].Mag()<HOWSMALLISTOOSMALL){
			    cout << "model_data.rhat_receive[i] mag is " << model_data.rhat_receive[i].Mag() << "\n";
			}
			model_data.rhat_receive[i].SetMag(1.);

		    }
		    else {

			// ----------------------------------------------------------------
			// If no ray solution was found, fall back to simple straight-line
			// source-to-station geometry.
			// ----------------------------------------------------------------
			model_data.rhat_launch[i].SetX(geom.station_coords[i][0]-geom.pulser_coords[0]);
			model_data.rhat_launch[i].SetY(geom.station_coords[i][1]-geom.pulser_coords[1]);
			model_data.rhat_launch[i].SetZ(geom.station_depths[i]-model_data.vdepth[i][idepth]);

			model_data.rhat_receive[i].SetX(model_data.rhat_launch[i][0]);
			model_data.rhat_receive[i].SetY(model_data.rhat_launch[i][1]);
			model_data.rhat_receive[i].SetZ(model_data.rhat_launch[i][2]);
		    }

		    // Store launch and receive polar angles (degrees)
		    model_data.vreceiveangle[i].push_back(model_data.rhat_receive[i].Theta()*DEGRAD);
		    model_data.vlaunchangle[i].push_back(model_data.rhat_launch[i].Theta()*DEGRAD);

		    // Final normalization / safety checks
		    if (model_data.rhat_launch[i].Mag()<HOWSMALLISTOOSMALL){
			cout << "rhat[i] mag is " << model_data.rhat_launch[i].Mag() << "\n";
		    }

		    model_data.rhat_launch[i].SetMag(1.);

		    if (model_data.rhat_receive[i].Mag()<HOWSMALLISTOOSMALL){
			cout << "model_data.rhat_receive[i] mag is " << model_data.rhat_receive[i].Mag() << "\n";
		    }

		    model_data.rhat_receive[i].SetMag(1.);
		}

		const sec loop_duration = clock::now() - before;
		loop_time += static_cast<double>( loop_duration.count() );
		// --------------------------------------------------------------------
		// Convert all accumulated per-depth vectors into ROOT TGraphs for plotting
		// and for later writing to output files.
		// --------------------------------------------------------------------
    convert_to_TGraphs(cfg, i, model_data);
    // ====================== end moved per-station body ======================
}
