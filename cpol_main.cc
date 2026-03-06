/*

This code describes the effects of birefringence at the South Pole

*/

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

#include <iostream>
#include <fstream>
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
#include "/data/user/alansalgo/ARA/IceRayTracing/IceRayTracing.h"
#include "TLegend.h"
#include "birefringence.hh"

using namespace std;

// ----------------------------------------------------------------------------
// Constants and global variables
// ----------------------------------------------------------------------------
const double PI=3.1415926;
const double CLIGHT=3.E8;
const double DEGRAD=180./PI;
const double NICE=1.78;

// Small numeric tolerances / discretization knobs
const double DELTAN=0.0005;
const double HOWSMALLISTOOSMALL=1.e-8;
const int UZAIRSTEP=5;

// Principal indices used for a simple "toy" dielectric tensor if not depth-dependent
// (often interpreted as principal axes of the indicatrix)
vector<double> nvec{1.77,1.7805,1.7815};

// Attenuation length (used for power/beam propagation)
const double L_ATTEN=2200.;

// Feet-to-meters conversion factor
// MFT = (1/100)*2.54*12 = 0.3048, i.e. 1 ft -> 0.3048 m
const double MFT=1./100.*2.54/1.*12.;

// ROOT color palettes used for plotting
int icolors[6]={kGreen+2,kRed+1,kOrange+1,kViolet+1,kBlue,kBlack};
int icolors_dave[6]={kBlack,kRed,kGreen,kBlue,kYellow,kBlack};

// -----------------------------------------
// Forward declarations
// -----------------------------------------    
void titles(TH2 *inH, TString title, TString xtitle, TString ytitle);
void titles(TGraph *inGr, TString title, TString xtitle, TString ytitle);

TCanvas *makePretty2Panel();
double Flipped(double theta_e1,double theta_e1_start);
void switchThem(double &thetaE_e1_Sclock,double &thetaE_e2_Sclock);


double getDeltaN(double alpha);
TVector3 getNewkandE(vector<double> nvec,double angle_iceflow,double n,TVector3 D,TVector3 kguess,
TVector3 &E);

TVector3 directionNextStep(vector<double>nvec, TVector3 E);

TVector3 rotateE(TVector3 epsilon, double angle_iceflow, TVector3 E);


// ----------------------------------------------------------------------------
// main
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {

    // ------------------------------------------------------------------------
    // Build dielectric tensor epsilon and its inverse from principal indices.
    // Here epsilon is represented as a TVector3 holding diagonal elements
    // (epsilon_x, epsilon_y, epsilon_z) in the principal basis.
    // ------------------------------------------------------------------------
    TVector3 inverse_epsilon;
    inverse_epsilon[0]=1/(nvec[0]*nvec[0]);
    inverse_epsilon[1]=1/(nvec[1]*nvec[1]);
    inverse_epsilon[2]=1/(nvec[2]*nvec[2]);

    TVector3 epsilon;
    epsilon[0]=nvec[0]*nvec[0];
    epsilon[1]=nvec[1]*nvec[1];
    epsilon[2]=nvec[2]*nvec[2];

    // ------------------------------------------------------------------------
    // Geometry setup
    // ------------------------------------------------------------------------
    const double DEPTH=-600.;
    const int ispecial=0;

    const double depth_special=-1100.;

    // Depth grid used for "special" evaluations (e.g., diagnostic depths)
    const double stepspecial=20.;
    const int NSPECIAL=50;
    double whichspecial[NSPECIAL];
    double startspecial=-1600.;
    for (int ispecial=0;ispecial<NSPECIAL;ispecial++) {
        whichspecial[ispecial]=startspecial+stepspecial*(double)ispecial;
    }

    const int WHICHPOL=0;
    const int DORAYTRACING=1;

    // ------------------------------------------------------------------------
    // Station definitions: A1..A5 plus ARIANNA
    // ------------------------------------------------------------------------
    const int NSTATIONS=6;
    const int minstation=0;
    const int maxstation=5;
    int igreatestdepth[NSTATIONS];
    int imostshallowdepth[NSTATIONS];

    // Horizontal pulser->station distances (in meters after MFT conversion)
    double horizontal_distances[NSTATIONS]={1257.,2353.,3146.,3199.,5179.8892,653.804525};
    string snames[NSTATIONS]={"A1","A2","A3","A4","A5","ARIANNA"};

    // ARIANNA azimuth/phi angle (radians)
    double PHI_ARIANNA=(90.+36.+46./60.+23./3600.+1.4)/DEGRAD;


    const int NSTEPS=100.;
    double min_altitude=-1000.;
    double max_altitude=-600.;

    double N_ICE=1.78;
    double DELTA_N=0.427;


    static double VOLTAGENORM=150.*11./7.;

    // ------------------------------------------------------------------------
    // Command-line options (getopt):
    //   -b : BIAXIAL flag
    //   -f : frequency in MHz
    //   -t : TX cross-pol angle in degrees
    //   -u : RX cross-pol angle in degrees
    //   -c : CONSTANTINDICATRIX flag
    //
    // BIAXIAL meaning:
    //   1  => fully biaxial: n1(z), n2(z), n3(z) all independent from files
    //   0  => uniaxial: n2(z) forced equal to n1(z), n3(z) read from file
    //  -1  => nearly isotropic: n2(z)=n1(z), n3(z)=n1(z)+1e-5 (break degeneracy)
    //
    // CONSTANTINDICATRIX:
    //   1 => enforce depth-independent indicatrix (take first depth sample and repeat)
    //   0 => keep depth dependence from files
    // ------------------------------------------------------------------------
    char clswitch;

    double freq=160.E6;             // default frequency [Hz]
    int CROSSPOLANGLE_TX_INT=7;
    int CROSSPOLANGLE_RX_INT=-7;
    int BIAXIAL=1;
    int CONSTANTINDICATRIX=0;

    if (argc>1) {
        while ((clswitch = getopt(argc, argv, "b:f:t:u:c:")) != EOF) {
            switch(clswitch) {
            case 'b':
                cout << optarg << "\n";
                BIAXIAL=atoi(optarg);
                cout << "biaxial " << BIAXIAL << endl;
                break;
            case 'f':
                cout << optarg << "\n";
                freq=(double)atof(optarg)*1.E6;
                cout << "freq " << freq << endl;
                break;
            case 't':
                CROSSPOLANGLE_TX_INT=atoi(optarg);
                cout << "CROSSPOLANGLE_TX_INT " << CROSSPOLANGLE_TX_INT << endl;
                break;
            case 'u':
                CROSSPOLANGLE_RX_INT=atoi(optarg);
                cout << "CROSSPOLANGLE_RX_INT " << CROSSPOLANGLE_RX_INT << endl;
                break;
            case 'c':
                CONSTANTINDICATRIX=atoi(optarg);
                cout << "CONSTANTINDICATRIX " << CONSTANTINDICATRIX << endl;
                break;
            }
        }
    }

    // Convert cross-pol angles to radians
    double CROSSPOLANGLE_TX=(double)CROSSPOLANGLE_TX_INT/DEGRAD;
    double CROSSPOLANGLE_RX=(double)CROSSPOLANGLE_RX_INT/DEGRAD;

    // ------------------------
    // Frequency 
    // ------------------------
    double freqmin=0.;
    double freqmax=1.E9;
    const int NFREQ=100;
    vector<double> vfreqs;


    for (int i=0;i<NFREQ;i++) {
        vfreqs.push_back(freqmin+(freqmax-freqmin)/(double)NFREQ*(double)i);
    }

    // ------------------------------------------------------------------------
    // Station coordinates and pulser coordinates (originally in feet; converted to meters)
    // ------------------------------------------------------------------------
    double pulser_coords[2]={42358.94,48974.2};

    double station_coords[NSTATIONS][2]={


        {38754., 51051.},
        {35481., 45369.},
        {32200., 51053.},
        {35478., 56737.},
        {32356., 39746.},


        {41153.,50381.75}
    };

    cout << "station_coords of ARIANNA is " << station_coords[5][0] << "\t" << station_coords[5][1] << "\n";
    for (int i=0;i<NSTATIONS;i++) {
        // Print horizontal distance in meters after conversion
        cout << "distance to station " << i+1 << " is " << sqrt((station_coords[i][0]-pulser_coords[0])*(station_coords[i][0]-pulser_coords[0]) + (station_coords[i][1]-pulser_coords[1])*(station_coords[i][1]-pulser_coords[1]))*MFT << "\n";
    }

    // Convert pulser and station coordinates to meters
    for (int j=0;j<2;j++) {
        pulser_coords[j]=pulser_coords[j]*MFT;
    }
    for (int i=0;i<6;i++) {
        for (int j=0;j<2;j++) {
            station_coords[i][j]=station_coords[i][j]*MFT;
        }
    }

    // ------------------------------------------------------------------------
    // Ice-flow direction and "ordinary" axis direction in horizontal plane.
    // angle_iceflow: azimuth defining principal axis orientation relative to x-y.
    // ------------------------------------------------------------------------
    const double angle_iceflow=(36.+ (46./60.) + (23./3600.) + 90.)/DEGRAD;
    TVector3 ordinary(cos(angle_iceflow),sin(angle_iceflow),0.);

    // Default TX/RX orientations (both along +z)
    TVector3 tx_orientation(0.,0.,1.);
    TVector3 rx1_orientation(0.,0.,1.);

    // ------------------------------------------------------------------------
    // Arrays that store station-wide results / diagnostics
    // ------------------------------------------------------------------------
    double beam[2]={0.};
    double deltan_exp[6];
    double deltan_obs[6];
    double station[6]={1,2,3,4,5,6};
    double station_depths[6]={-80.,-180.,-180.,-180.,-180.,-1.}; // receiver depths (m)
    double zeroes[6]={0.};
    double deltan_obs_err[6]={0.};

    double alpha[6];    // station azimuth relative to ice-flow axis
    double alpha_deg[6];
    TVector3 pulsertostationhat_specialdepth[6];    
    TVector3 pulsertostation[6];

    // ------------------------------------------------------------------------
    // Compute geometry-dependent alpha for each station and expected delta-n
    // ------------------------------------------------------------------------
    for (int i=minstation;i<=maxstation;i++) {

        // Vector from pulser to station, using station depth and pulser depth
        pulsertostationhat_specialdepth[i][0]=(station_coords[i][0]-pulser_coords[0]);
        pulsertostationhat_specialdepth[i][1]=(station_coords[i][1]-pulser_coords[1]);
        pulsertostationhat_specialdepth[i][2]=station_depths[i]-(DEPTH);

        if (pulsertostationhat_specialdepth[i].Mag()<HOWSMALLISTOOSMALL)
            cout << "pulsertostationhat_specialdepth mag is " << pulsertostationhat_specialdepth[i].Mag() << "\n";

        // Normalize to unit vector
        pulsertostationhat_specialdepth[i].SetMag(1.);

        // Station azimuth in x-y plane
        double y=station_coords[i][1]-pulser_coords[1];
        double x=station_coords[i][0]-pulser_coords[0];
        double angle_thisstation=atan2(y,x);

        // alpha = (station azimuth) - (ice-flow azimuth), wrapped to [-pi, +pi]
        double thisalpha=angle_thisstation-angle_iceflow;

        if (thisalpha>PI)
        thisalpha-=2.*PI;
        if (thisalpha<-1.*PI)
        thisalpha+=2.*PI;
        alpha[i]=thisalpha;
        cout << "A" << i+1 << ": angle_thisstation, angle_iceflow, alpha are " << angle_thisstation*DEGRAD << "\t" << angle_iceflow*DEGRAD << "\t" << alpha[i]*DEGRAD << ".\n";

        alpha_deg[i]=alpha[i]*DEGRAD;

        // Outputs from getDeltaN: eigenstates p_e1, p_e2 and indices n_e1, n_e2
        TVector3 p_e1;
        TVector3 p_e2;
        double n_e1;
        double n_e2;

        // Expected birefringence delta n magnitude for that propagation direction
        deltan_exp[i]=getDeltaN(BIAXIAL,nvec,pulsertostationhat_specialdepth[i],angle_iceflow,
        n_e1,n_e2,p_e1,p_e2);
    }

    // ------------------------------------------------------------------------
    // Special handling for ARIANNA: place it ~650 m away along ice-flow direction
    // and compute alpha relative to "ordinary" axis.
    // ------------------------------------------------------------------------
    double x_arianna=650.*cos(angle_iceflow);
    double y_arianna=650.*sin(angle_iceflow);
    double z_arianna=0.;
    pulsertostationhat_specialdepth[5][0]=x_arianna-pulser_coords[0];
    pulsertostationhat_specialdepth[5][1]=y_arianna-pulser_coords[1];
    pulsertostationhat_specialdepth[5][2]=z_arianna-(DEPTH);

    if (pulsertostationhat_specialdepth[5].Mag()<HOWSMALLISTOOSMALL)
    cout << "pulsertostationhat_specialdepth[5] mag is " << pulsertostationhat_specialdepth[5].Mag() << "\n";

    pulsertostationhat_specialdepth[5].SetMag(1.);

    double thisalpha=acos(pulsertostationhat_specialdepth[5].Dot(ordinary));
    if (thisalpha>PI)
    thisalpha-=2.*PI;
    if (thisalpha<-1.*PI)
    thisalpha+=2.*PI;
    alpha[5]=thisalpha;

    TVector3 p_e1;
    TVector3 p_e2;
    double n_e1;
    double n_e2;


    // ------------------------------------------------------------------------
    // Containers for reading data and for storing outputs.
    // Each outer vector index corresponds to a station.
    // ------------------------------------------------------------------------
    string line;
    std::vector< std::vector<double> > videpth;
    std::vector< std::vector<double> > vdepth;

    std::vector< std::vector<double> > vdepth_data;
    std::vector< std::vector<double> > vdepth_data_err;

    std::vector< std::vector<double> > vreversedepth;


    std::vector< std::vector<double> > vreceiveangle;
    std::vector< std::vector<double> > vlaunchangle;
    std::vector< std::vector<double> > voutput6;
    std::vector< std::vector<double> > voutput7;
    std::vector< std::vector<double> > voutput8;


    std::vector< std::vector<double> > vtotal_distances;
    std::vector< std::vector<double> > vtotal_distances_err;

    std::vector< std::vector<double> > vsnrmax;
    std::vector< std::vector<double> > vsnrmax_err;

    // Various angle diagnostics between propagation direction / eigenvectors
    std::vector< std::vector<double> > vangle_khat_0_khat_1_2;
    std::vector< std::vector<double> > vangle_khat_0_khat_2_2;
    std::vector< std::vector<double> > vangle_khat_1_2_khat_2_2;
    std::vector< std::vector<double> > vangle_Shat_e1_khat;
    std::vector< std::vector<double> > vangle_Shat_e2_khat;

    // Ray-tracing path containers (two rays per station?)
    TVector3 raypath[2][6];
    TVector3 raypath_n[2][6];
    std::vector< std::vector< double > > vistep;
    std::vector< std::vector< std::vector<double> > > vraypos;
    std::vector< std::vector<double> > vraypath_ne1;
    std::vector< std::vector<double> > vraypath_ne2;

    // Beam / attenuation vs depth containers
    std::vector< std::vector<double> > vrxdepth_beam1;
    std::vector< std::vector<double> > vrxdepth_beam2;
    std::vector< std::vector<double> > vtxdepth_beam1;
    std::vector< std::vector<double> > vtxdepth_beam2;

    std::vector< std::vector<double> > vrxdepth_atten;
    std::vector< std::vector<double> > vrxdepth_notflipped;
    std::vector< std::vector<double> > vrxdepth_atten_beam;
    std::vector< std::vector<double> > vrxdepth_atten_power;
    std::vector< std::vector<double> > vrxdepth_atten_beam_power;

    // Angular evolution of E-field components (theta1/theta2) along the path
    std::vector< std::vector<double> > vtxdepth_theta1;
    std::vector< std::vector<double> > vtxdepth_theta2;
    std::vector< std::vector<double> > vrxdepth_theta1;
    std::vector< std::vector<double> > vrxdepth_theta2;

    // Same angles, but expressed in a "S-clock" coordinate convention
    std::vector< std::vector<double> > vtxdepth_theta1_Sclock;
    std::vector< std::vector<double> > vtxdepth_theta2_Sclock;
    std::vector< std::vector<double> > vrxdepth_theta1_Sclock;
    std::vector< std::vector<double> > vrxdepth_theta2_Sclock;

    // E-field polarization angles (and S-clock versions)
    std::vector< std::vector<double> > vtxdepthE_theta1;
    std::vector< std::vector<double> > vtxdepthE_theta2;
    std::vector< std::vector<double> > vrxdepthE_theta1;
    std::vector< std::vector<double> > vrxdepthE_theta2;
    std::vector< std::vector<double> > vtxdepthE_theta1_Sclock;
    std::vector< std::vector<double> > vtxdepthE_theta2_Sclock;
    std::vector< std::vector<double> > vrxdepthE_theta1_Sclock;
    std::vector< std::vector<double> > vrxdepthE_theta2_Sclock;

    // Dispersion contributions along the path
    std::vector< std::vector<double> > vtxdepth_dispersion1;
    std::vector< std::vector<double> > vtxdepth_dispersion2;

    // Dot products of S/E/D unit vectors with something (TX frame)
    std::vector< std::vector<double> > vdotShats_tx;
    std::vector< std::vector<double> > vdotEhats_tx;
    std::vector< std::vector<double> > vdotDhats_tx;

    // ------------------------------------------------------------------------
    // A very large set of per-station waveform / voltage / power bookkeeping
    // containers (r1 and r2 likely refer to two rays / eigenmodes)
    // ------------------------------------------------------------------------
    std::vector<std::vector<double>> vV1_r1;
    std::vector<std::vector<double>> vV1_r1_lpda;
    std::vector<std::vector<double>> vV2_r1;
    std::vector<std::vector<double>> vE1_r1;
    std::vector<std::vector<double>> vE2_r1;
    std::vector<std::vector<double>> vE1_r2;
    std::vector<std::vector<double>> vE2_r2;
    std::vector<std::vector<double>> vV2_r1_lpda;
    std::vector<std::vector<double>> vV1squared_r1;
    std::vector<std::vector<double>> vV2squared_r1;
    std::vector<std::vector<double>> vV1V2_r1;
    std::vector<std::vector<double>> vV1V2_r1_lpda;
    std::vector<std::vector<double>> vV1V2_r2;
    std::vector<std::vector<double>> voppositeV1V2_r2;
    std::vector<std::vector<double>> voppositeV1V2_r1;
    std::vector<std::vector<double>> venvelope_minus_r1;
    std::vector<std::vector<double>> venvelope_minus_r1_lpda;
    std::vector<std::vector<double>> venvelope_plus_r1;
    std::vector<std::vector<double>> venvelope_plus_r1_lpda;
    std::vector<std::vector<double>> vvenvelope_minus_r1;
    std::vector<std::vector<double>> vvenvelope_plus_r1;
    std::vector<std::vector<double>> vEenvelope_minus_r1;
    std::vector<std::vector<double>> vEenvelope_plus_r1;
    std::vector<std::vector<double>> vSenvelope_minus_r1;
    std::vector<std::vector<double>> vSenvelope_plus_r1;
    std::vector<std::vector<double>> vpower_r1;
    std::vector<std::vector<double>> vpoynting_r1;
    std::vector<std::vector<double>> vpower_r1_lpda;

    std::vector<std::vector<double>> vV1_r2;
    std::vector<std::vector<double>> vV2_r2;
    std::vector<std::vector<double>> vV1squared_r2;
    std::vector<std::vector<double>> vV2squared_r2;
    std::vector<std::vector<double>> venvelope_minus_r2;
    std::vector<std::vector<double>> venvelope_plus_r2;
    std::vector<std::vector<double>> vvenvelope_minus_r2;
    std::vector<std::vector<double>> vvenvelope_plus_r2;
    std::vector<std::vector<double>> vEenvelope_minus_r2;
    std::vector<std::vector<double>> vEenvelope_plus_r2;
    std::vector<std::vector<double>> vSenvelope_minus_r2;
    std::vector<std::vector<double>> vSenvelope_plus_r2;
    std::vector<std::vector<double>> vpower_r2;
    std::vector<std::vector<double>> vpoynting_r2;
    std::vector<std::vector<double>> vvoltage_r1;
    std::vector<std::vector<double>> vfield_r1;
    std::vector<std::vector<double>> vvoltage_r1_lpda;
    std::vector<std::vector<double>> vvoltage_r2;
    std::vector<std::vector<double>> vfield_r2;

    // Polarization vectors along the ray (per station)
    std::vector<std::vector<TVector3>> pol_r1;
    std::vector<std::vector<TVector3>> pol_r2;

    // Dielectric tensor eigenvalues / polarization angles at TX/RX
    std::vector<std::vector<double>> vepsilon1_tx;
    std::vector<std::vector<double>> vepsilon2_tx;
    std::vector<std::vector<double>> vdiffepsilon_tx;
    std::vector<std::vector<double>> vpolarization_Psi_rx;
    std::vector<std::vector<double>> vpolarization_Omega_rx;
    std::vector<std::vector<double>> vEpolarization_Psi_rx;
    std::vector<std::vector<double>> vEpolarization_Omega_rx;
    std::vector<std::vector<double>> vEpolarization_reversedepth_Psi_rx;
    std::vector<std::vector<double>> vEpolarization_reversedepth_Omega_rx;
    std::vector<std::vector<double>> vpolarization_reversedepth_Psi_rx;
    std::vector<std::vector<double>> vpolarization_reversedepth_Omega_rx;


    std::vector<std::vector<double>> vepsilon1_rx;
    std::vector<std::vector<double>> vepsilon2_rx;
    std::vector<std::vector<double>> vdiffepsilon_rx;

    // ------------------------------------------------------------------------
    // Resize all per-station containers to NSTATIONS (=6).
    // This sets up empty vectors per station, ready to push_back later.
    // ------------------------------------------------------------------------
    vdotShats_tx.resize(6);
    vdotEhats_tx.resize(6);
    vdotDhats_tx.resize(6);

    vepsilon1_tx.resize(6);
    vepsilon2_tx.resize(6);
    vdiffepsilon_tx.resize(6);

    vepsilon1_rx.resize(6);
    vepsilon2_rx.resize(6);
    vdiffepsilon_rx.resize(6);

    vV1_r1.resize(6);
    vV2_r1.resize(6);
    vE1_r1.resize(6);
    vE2_r1.resize(6);
    vV1_r1_lpda.resize(6);
    vV2_r1_lpda.resize(6);
    vV1_r2.resize(6);
    vV2_r2.resize(6);
    vE1_r2.resize(6);
    vE2_r2.resize(6);

    vV1V2_r1.resize(6);
    vV1V2_r1_lpda.resize(6);
    vV1V2_r2.resize(6);
    voppositeV1V2_r2.resize(6);
    voppositeV1V2_r1.resize(6);


    vV1squared_r1.resize(6);
    vV2squared_r1.resize(6);
    vV1squared_r2.resize(6);
    vV2squared_r2.resize(6);

    vpower_r1.resize(6);
    vpoynting_r1.resize(6);
    vpower_r1_lpda.resize(6);
    vpower_r2.resize(6);
    vpoynting_r2.resize(6);
    vvoltage_r1.resize(6);
    vfield_r1.resize(6);
    vvoltage_r1_lpda.resize(6);
    vvoltage_r2.resize(6);
    vfield_r2.resize(6);

    vpolarization_Psi_rx.resize(6);
    vpolarization_Omega_rx.resize(6);
    vEpolarization_Psi_rx.resize(6);
    vEpolarization_Omega_rx.resize(6);
    vEpolarization_reversedepth_Psi_rx.resize(6);
    vEpolarization_reversedepth_Omega_rx.resize(6);
    vpolarization_reversedepth_Psi_rx.resize(6);
    vpolarization_reversedepth_Omega_rx.resize(6);

    venvelope_minus_r1.resize(6);
    vSenvelope_minus_r1.resize(6);
    venvelope_minus_r1_lpda.resize(6);
    venvelope_plus_r1.resize(6);
    venvelope_plus_r1_lpda.resize(6);
    venvelope_minus_r2.resize(6);
    vSenvelope_minus_r2.resize(6);
    venvelope_plus_r2.resize(6);
    vSenvelope_plus_r2.resize(6);
    vvenvelope_minus_r1.resize(6);
    vvenvelope_plus_r1.resize(6);
    vEenvelope_minus_r1.resize(6);
    vEenvelope_plus_r1.resize(6);
    vSenvelope_plus_r1.resize(6);
    vvenvelope_minus_r2.resize(6);
    vvenvelope_plus_r2.resize(6);
    vEenvelope_minus_r2.resize(6);
    vEenvelope_plus_r2.resize(6);

    // Ray path / stepping containers
    vistep.resize(6);
    vraypos.resize(6);
    vraypath_ne1.resize(6);
    vraypath_ne2.resize(6);
    vrxdepth_beam1.resize(6);
    vrxdepth_beam2.resize(6);
    vtxdepth_beam1.resize(6);
    vtxdepth_beam2.resize(6);

    vrxdepth_atten.resize(6);
    vrxdepth_atten_beam.resize(6);
    vrxdepth_atten_power.resize(6);
    vrxdepth_atten_beam_power.resize(6);
    vtxdepth_theta1.resize(6);
    vtxdepth_theta2.resize(6);
    vrxdepth_theta1.resize(6);
    vrxdepthE_theta2.resize(6);
    vrxdepthE_theta1.resize(6);
    vrxdepth_theta2.resize(6);
    vtxdepth_theta1_Sclock.resize(6);
    vtxdepth_theta2_Sclock.resize(6);
    vrxdepth_theta1_Sclock.resize(6);
    vrxdepth_theta2_Sclock.resize(6);
    vtxdepth_dispersion1.resize(6);
    vtxdepth_dispersion2.resize(6);
    vtxdepthE_theta1.resize(6);
    vtxdepthE_theta2.resize(6);
    vtxdepthE_theta1_Sclock.resize(6);
    vtxdepthE_theta2_Sclock.resize(6);
    vrxdepthE_theta1_Sclock.resize(6);
    vrxdepthE_theta2_Sclock.resize(6);

    // vraypos[station] is a vector of 3 vectors: x(z), y(z), z(z)
    for (int istations=0;istations<6;istations++) {
        vraypos[istations].resize(3);
    }

    // ------------------------------------------------------------------------
    // Input data files (Dave logs + n1/n2/n3 refractive index profiles)
    // ------------------------------------------------------------------------
    string sfile;
    if (WHICHPOL==0)
    sfile="dave_data/day359_pol0.log";
    if (WHICHPOL==1)
    sfile="dave_data/day359_pol1.log";


    ifstream myfile(sfile.c_str());
    ifstream davea5file("dave_data/a5_amyformat.txt");

    // Depth-dependent principal indices (text files)
    string sn1file="data/n1.txt";
    string sn2file="data/n2.txt";
    string sn3file="data/n3.txt";

    ifstream n1file(sn1file.c_str());
    ifstream n2file(sn2file.c_str());
    ifstream n3file(sn3file.c_str());

    // Depth arrays and refractive index arrays for each principal axis
    vector<double> vdepths_n1;
    vector<double> vdepths_n2;
    vector<double> vdepths_n3;

    vector<double> n1vec;
    vector<double> n2vec;
    vector<double> n3vec;

    // Derived quantity vs depth (e.g., V(n1,n2,n3) from birefringence.hh)
    std::vector<double>  vV;

    // ------------------------------------------------------------------------
    // Read in depth->n profiles from files.
    // Each file is assumed to have 1 header token, then 81 rows of: depth n
    // Convention: store z = -depth so deeper points have more negative z.
    // ------------------------------------------------------------------------
    int NDEPTHS_NS=81;

    string stemp;
    double thisdepth,thisn;
    double firstn1, firstn2, firstn3;

    // ---- n1(z) ----
    n1file >> stemp; // discard header token
    for (int i=0;i<NDEPTHS_NS;i++) {
        n1file >> thisdepth >> thisn;
        vdepths_n1.push_back(-1.*thisdepth);   // flip sign
        n1vec.push_back(thisn);
    }

    // ---- n2(z) ----
    n2file >> stemp;
    for (int i=0;i<NDEPTHS_NS;i++) {
        n2file >> thisdepth >> thisn;
        vdepths_n2.push_back(-1.*thisdepth);

        // If fully biaxial, use n2 file.
        // If uniaxial or isotropic-like modes, force n2=n1.
        if (BIAXIAL==1)
            n2vec.push_back(thisn);
        else if (BIAXIAL==0 || BIAXIAL==-1)
            n2vec.push_back(n1vec[i]);
    }

    // ---- n3(z) ----
    n3file >> stemp;
    for (int i=0;i<NDEPTHS_NS;i++) {
        n3file >> thisdepth >> thisn;
        vdepths_n3.push_back(-1.*thisdepth);

        // If uniaxial or biaxial, use n3 file.
        // If BIAXIAL==-1 (nearly isotropic), set n3=n1+1e-5 to break degeneracy.
        if (BIAXIAL==0 || BIAXIAL==1)
            n3vec.push_back(thisn);
        else if (BIAXIAL==-1)
            n3vec.push_back(n1vec[i]+1.E-5);
    }

    // ------------------------------------------------------------------------
    // Optional: "constant indicatrix" -> set n1,n2,n3 to their first depth values
    // at all depths. This removes depth dependence but retains birefringence.
    // ------------------------------------------------------------------------
    if (CONSTANTINDICATRIX==1) {
        firstn1=n1vec[0];
        firstn2=n2vec[0];
        firstn3=n3vec[0];

        int thissize=(int)n1vec.size();

        n1vec.clear();
        n2vec.clear();
        n3vec.clear();

        for (int i=0;i<thissize;i++) {
            n1vec.push_back(firstn1);
            n2vec.push_back(firstn2);
            n3vec.push_back(firstn3);
        }
    }

    // Debug printout of sizes and raw indices
    cout << "sizes are " << n1vec.size() << "\t" << n2vec.size() << "\t" << n3vec.size() << "\n";
    cout << "n's are \n";
    for (int i=0;i<NDEPTHS_NS;i++) {
        cout << "n1, n2, n3 are " << n1vec[i] << "\t" << n2vec[i] << "\t" << n3vec[i] << "\n";
    }

    // ------------------------------------------------------------------------
    // Smooth n1/n2/n3 vs depth using a simple moving average (boxcar) filter.
    // NSMOOTH=5 -> replace each interior point with average over i-2..i+2.
    // Edges are left unchanged for the first/last 2 points.
    // ------------------------------------------------------------------------
    vector<double> tmp;
    tmp.resize(n1vec.size());
    int NSMOOTH=5;
    int min=(int)(((double)NSMOOTH)/2.);    // = 2 for NSMOOTH=5

    // Copy edges unchanged
    for (int i=0;i<min;i++) {
        tmp[i]=n1vec[i];
    }
    for (int i=n1vec.size()-(NSMOOTH-min);i<n1vec.size();i++) {
        tmp[i]=n1vec[i];
    }

    // Smooth interior points
    for (int i=min;i<n1vec.size()-(NSMOOTH-min);i++) {
        double tmpdouble=0.;
        for (int j=i-min;j<i+(NSMOOTH-min);j++) {
            tmpdouble+=n1vec[j];
        }
        tmpdouble=tmpdouble/(double)NSMOOTH;
        tmp[i]=tmpdouble;
    }
    n1vec=tmp;

    // Repeat smoothing for n2
    tmp.clear();
    tmp.resize(n2vec.size());

    min=(int)(((double)NSMOOTH)/2.);
    for (int i=0;i<min;i++) {
        tmp[i]=n2vec[i];
    }
    for (int i=n2vec.size()-(NSMOOTH-min);i<n2vec.size();i++) {
        tmp[i]=n2vec[i];
    }
    for (int i=min;i<n2vec.size()-(NSMOOTH-min);i++) {
        double tmpdouble=0.;
        for (int j=i-min;j<i+(NSMOOTH-min);j++) {
            tmpdouble+=n2vec[j];
        }
        tmpdouble=tmpdouble/(double)NSMOOTH;
        tmp[i]=tmpdouble;
    }
    n2vec=tmp;

    // Repeat smoothing for n3
    tmp.clear();
    tmp.resize(n3vec.size());

    min=(int)(((double)NSMOOTH)/2.);
    for (int i=0;i<min;i++) {
        tmp[i]=n3vec[i];
    }
    for (int i=n3vec.size()-(NSMOOTH-min);i<n3vec.size();i++) {
        tmp[i]=n3vec[i];
    }
    for (int i=min;i<n3vec.size()-(NSMOOTH-min);i++) {
        double tmpdouble=0.;
        for (int j=i-min;j<i+(NSMOOTH-min);j++) {
            tmpdouble+=n3vec[j];
        }
        tmpdouble=tmpdouble/(double)NSMOOTH;
        tmp[i]=tmpdouble;
    }
    n3vec=tmp;

    // ------------------------------------------------------------------------
    // Build ROOT graphs of n1(z), n2(z), n3(z) for interpolation via Eval(z)
    // ------------------------------------------------------------------------
    TGraph *gn1=new TGraph(n1vec.size(),&vdepths_n1[0],&n1vec[0]);
    TGraph *gn2=new TGraph(n2vec.size(),&vdepths_n2[0],&n2vec[0]);
    TGraph *gn3=new TGraph(n3vec.size(),&vdepths_n3[0],&n3vec[0]);
    TGraph *g_V;

    // Compute derived quantity V(z) from the birefringence model (getV in birefringence.hh)
    vector<double> nvec_tmp;
    nvec_tmp.resize(3);

    for (int i=0;i<n1vec.size();i++) {
        // Evaluate indices at each depth (using the smoothed/interpolated curves)
        nvec_tmp[0]=gn1->Eval(vdepths_n1[i]);
        nvec_tmp[1]=gn2->Eval(vdepths_n2[i]);
        nvec_tmp[2]=gn3->Eval(vdepths_n3[i]);
        vV.push_back(getV(nvec_tmp));
    }

    // Graph of V(z)
    g_V=new TGraph(vdepths_n1.size(),&vdepths_n1[0],&vV[0]);

    // ------------------------------------------------------------------------
    // Declare lots of ROOT TGraphs that will be filled later 
    // ------------------------------------------------------------------------
    TGraph *graypath_z_x[6];
    TGraph *graypath_z_y[6];
    TGraph *graypath_y_x[6];

    TGraph *graypath_n[6];

    TGraph *grxdepth_atten[6];
    TGraph *grxdepth_atten_beam[6];
    TGraph *grxdepth_atten_power[6];
    TGraph *grxdepth_atten_beam_power[6];
    TGraph *grxdepth_beam1[6];
    TGraph *grxdepth_beam2[6];
    TGraph *gtxdepth_beam1[6];
    TGraph *gtxdepth_beam2[6];
    TGraph *gtxdepth_theta1[6];
    TGraph *gtxdepth_theta2[6];
    TGraph *grxdepth_theta1[6];
    TGraph *grxdepth_theta2[6];
    TGraph *grxdepthE_theta1[6];
    TGraph *grxdepthE_theta2[6];
    TGraph *gtxdepth_theta1_Sclock[6];
    TGraph *gtxdepth_theta2_Sclock[6];
    TGraph *grxdepth_theta1_Sclock[6];
    TGraph *grxdepth_theta2_Sclock[6];
    TGraph *gtxdepth_dispersion1[6];
    TGraph *gtxdepth_dispersion2[6];
    TGraph *gtxdepthE_theta1[6];
    TGraph *gtxdepthE_theta2[6];
    TGraph *gtxdepthE_theta1_Sclock[6];
    TGraph *gtxdepthE_theta2_Sclock[6];
    TGraph *grxdepthE_theta1_Sclock[6];
    TGraph *grxdepthE_theta2_Sclock[6];

    TGraph *gdotShats_tx[6];
    TGraph *gdotEhats_tx[6];
    TGraph *gdotDhats_tx[6];

    TGraph *gsnrmax[6];
    TGraph *g_idepth[6];

    // Waveform-related graphs (per station)
    TGraph *gV1_r1[6];
    TGraph *gV2_r1[6];
    TGraph *gV1squared_r1[6];
    TGraph *gV2squared_r1[6];
    TGraph *gV1V2_r1[6];
    TGraph *gV1V2_r2[6];
    TGraph *goppositeV1V2_r1[6];
    TGraph *goppositeV1V2_r2[6];
    TGraph *gpower_r1[6];
    TGraph *gpower_r2[6];
    TGraph *gvoltage_r1[6];
    TGraph *gvoltage_r2[6];
    TGraph *gfield_r1[6];
    TGraph *gfield_r2[6];
    TGraph *genvelope_minus_r1[6];
    TGraph *genvelope_minus_r2[6];
    TGraph *genvelope_plus_r1[6];
    TGraph *genvelope_plus_r2[6];
    TGraph *gvenvelope_minus_r1[6];
    TGraph *gvenvelope_minus_r2[6];
    TGraph *gvenvelope_plus_r1[6];
    TGraph *gvenvelope_plus_r2[6];
    TGraph *gEenvelope_minus_r1[6];
    TGraph *gEenvelope_minus_r2[6];
    TGraph *gEenvelope_plus_r1[6];
    TGraph *gEenvelope_plus_r2[6];

    TGraph *gV1_r2[6];
    TGraph *gV2_r2[6];
    TGraph *gV1squared_r2[6];
    TGraph *gV2squared_r2[6];

    // Dielectric eigenvalues at TX/RX (per station)
    TGraph *gepsilon1_tx[6];
    TGraph *gepsilon2_tx[6];
    TGraph *gdiffepsilon_tx[6];

    TGraph *gepsilon1_rx[6];
    TGraph *gepsilon2_rx[6];
    TGraph *gdiffepsilon_rx[6];

    // Polarization angles (per station)
    TGraph *gpolarization_Omega_rx[6];
    TGraph *gpolarization_Psi_rx[6];
    TGraph *gEpolarization_Omega_rx[6];
    TGraph *gEpolarization_Psi_rx[6];

    TGraph *gpolarization_reversedepth_Omega_rx[6];
    TGraph *gpolarization_reversedepth_Psi_rx[6];

    TGraph *gEpolarization_reversedepth_Omega_rx[6];
    TGraph *gEpolarization_reversedepth_Psi_rx[6];

    // Misc outputs from ray tracer / solver
    TGraph *g_receive_launch[6];
    TGraph *g_receive[6];
    TGraph *g_launch[6];
    TGraph *g_output6[6];
    TGraph *g_output7[6];
    TGraph *g_output8[6];

    // ------------------------------------------------------------------------
    // Read Dave's data: parse station/day/pol/depth/snr
    // Store depth, snr, distance, and estimate an error from last-3-point RMS.
    // ------------------------------------------------------------------------
    int NSHOTS=640;

    // Resize containers for data (outer size = number of stations)
    vdepth_data.resize(5);
    vdepth.resize(6);
    vreceiveangle.resize(6);
    vlaunchangle.resize(6);
    voutput6.resize(6);
    voutput7.resize(6);
    voutput8.resize(6);
    videpth.resize(6);
    vreversedepth.resize(6);
    vangle_Shat_e1_khat.resize(6);
    vangle_Shat_e2_khat.resize(6);
    vangle_khat_0_khat_1_2.resize(6);
    vangle_khat_0_khat_2_2.resize(6);
    vangle_khat_1_2_khat_2_2.resize(6);
    vsnrmax.resize(6);
    vtotal_distances.resize(5);
    vdepth_data_err.resize(5);
    vsnrmax_err.resize(6);
    vtotal_distances_err.resize(5);

    double running_rms=0.;
    double running_mean=0.;


    if (myfile.is_open())
    {
        for (int i=0;i<NSHOTS;i++) {

            int this_station, this_day, this_pol;
            double this_depth, this_snrmax;

            // Format: station day pol depth snr
            myfile >> this_station >> this_day >> this_pol >> this_depth >> this_snrmax;

            // Rescale SNR to a common reference distance (A1 at -1000 m) for comparisons
            this_snrmax=this_snrmax*sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1])/sqrt((-1000.-station_depths[0])*(-1000.-station_depths[0])+horizontal_distances[0]*horizontal_distances[0]);

            // Keep positive SNR, and exclude A5 pol0 special case
            if (this_snrmax>0. && !(this_station==5 && this_pol==0)) {
                
                // Store depth and SNR for that station
                vdepth_data[this_station-1].push_back(this_depth);
                vsnrmax[this_station-1].push_back(this_snrmax);

                // Store total geometric distance (straight-line) from pulser to receiver
                vtotal_distances[this_station-1].push_back(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1]));

                // Estimate an uncertainty using RMS of last 3 SNR points (when available)
                if (i>2) {
                    running_mean=0.;
                    running_rms=0.;
                    for (int j=0;j<3;j++) {
                        running_mean+=vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1];
                    }
                    running_mean=running_mean/3.;

                    for (int j=0;j<3;j++) {
                        running_rms+=(vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1]-running_mean)*(vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1]-running_mean);
                    }
                    running_rms=sqrt(running_rms/2.);
                }

                // Store errors (depth/distance errors set to 0 here)
                vsnrmax_err[this_station-1].push_back(running_rms);
                vdepth_data_err[this_station-1].push_back(0.);
                vtotal_distances_err[this_station-1].push_back(0.);
            }
        }

        myfile.close();
    }

    // ------------------------------------------------------------------------
    // Special handling for A5: for WHICHPOL==0, read a separate file in Amy format.
    // ------------------------------------------------------------------------
    NSHOTS=33;

    if (WHICHPOL==0) {
        if (davea5file.is_open())
        {
            cout << "i'm reading dave's a5 file.\n";
            for (int i=0;i<NSHOTS;i++) {

                int this_station, this_day, this_pol;
                double this_depth, this_snrmax;
                string stemp;
                int this_channel;
                this_station=5;

                // Format: <token> <channel> <depth> <snr>
                davea5file >> stemp >> this_channel >> this_depth >> this_snrmax;

                // Keep only positive-SNR entries
                if (this_snrmax>0.) {
                    vdepth_data[this_station-1].push_back(this_depth);

                    // Store measured SNR and corresponding straight-line pulser->station distance
                    vsnrmax[this_station-1].push_back(this_snrmax);
                    vtotal_distances[this_station-1].push_back(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1]));

                    // Same running RMS estimate as above:
                    // use the last 3 stored SNR values to estimate a local scatter / uncertainty
                    if (i>2) {
                        running_mean=0.;
                        running_rms=0.;
                        for (int j=0;j<3;j++) {
                            running_mean+=vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1];
                        }
                        running_mean=running_mean/3.;

                        for (int j=0;j<3;j++) {
                            running_rms+=(vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1]-running_mean)*(vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1]-running_mean);
                        }
                        // Divide by 2 since this is an RMS estimate over 3 points (N-1 denominator)
                        running_rms=sqrt(running_rms/2.);

                    }
                    // Store per-point uncertainties
                    vsnrmax_err[this_station-1].push_back(running_rms);

                    // Depth and distance uncertainties are taken as zero here
                    vdepth_data_err[this_station-1].push_back(0.);
                    vtotal_distances_err[this_station-1].push_back(0.);
                }
            }

            davea5file.close();

        }
    }

    // ------------------------------------------------------------------------
    // Build the pulser-depth scan used in the model.
    // For all stations, use a common grid:
    //   start at -400 m
    //   step downward by 10 m
    //   take 200 depth values
    //
    // vdepth[station] stores the pulser depths,
    // videpth[station] stores the corresponding integer-like indices.
    // ------------------------------------------------------------------------
    double arianna_minpulserdepth=-400.;
    double arianna_pulserstep=10.;
    int NARIANNA_PULSER=200;
    for (int istations=0;istations<NSTATIONS;istations++) {
        for (int i=0;i<NARIANNA_PULSER;i++) {
            vdepth[istations].push_back(arianna_minpulserdepth-(double)arianna_pulserstep*(double)i);
            videpth[istations].push_back((double)i);
        }
    }

    // Build a reversed-depth version, useful for plots where increasing pulser depth
    // (positive number) is more natural than negative height.
    for (int istations=minstation;istations<=maxstation;istations++) {
        for (int j=0;j<vdepth[istations].size();j++) {
            vreversedepth[istations].push_back(-1.*vdepth[istations][vdepth[istations].size()-1-j]);
        }
    }

    // ------------------------------------------------------------------------
    // ROOT graphs:
    //   gsnrmax[i] = measured SNR vs pulser depth for station i
    //   g_idepth[i] = mapping from physical depth -> index in the depth scan
    //
    // Note: station index 5 (ARIANNA) has no gsnrmax built here.
    // ------------------------------------------------------------------------
    for (int i=minstation;i<=maxstation;i++) {
        if (i!=5)
        gsnrmax[i]=new TGraph(vsnrmax[i].size(),&vdepth_data[i][0],&vsnrmax[i][0]);
        g_idepth[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&videpth[i][0]);


    }

    // ------------------------------------------------------------------------
    // For each station, find:
    //   igreatestdepth   = index of deepest pulser depth (most negative)
    //   imostshallowdepth= index of shallowest pulser depth (least negative)
    // ------------------------------------------------------------------------
    for (int istations=minstation;istations<=maxstation;istations++) {
        igreatestdepth[istations]=0;
        imostshallowdepth[istations]=0;
        for (int idepth=0;idepth<vdepth[istations].size();idepth++) {
            if (vdepth[istations][idepth]<vdepth[istations][igreatestdepth[istations]])
            igreatestdepth[istations]=idepth;
            if (vdepth[istations][idepth]>vdepth[istations][imostshallowdepth[istations]])
            imostshallowdepth[istations]=idepth;
        }

    }

    // ------------------------------------------------------------------------
    // Geometric / polarization bookkeeping vectors used later during ray tracing
    // and field decomposition.
    //
    // A brief interpretation of some names:
    //   rhat         = launch-direction unit vector
    //   rhat_receive = receive-direction unit vector
    //   p_o, p_e     = ordinary / extraordinary polarization directions
    //   Pt, Pr1, Pr2 = transmitter / receiver basis vectors
    // ------------------------------------------------------------------------
    vector <TVector3> rhat;
    vector <TVector3> rhat_receive;
    vector <TVector3> p_o;
    vector <TVector3> p_e;
    vector <TVector3> p_e_rx;
    vector <TVector3> p_o_rx;
    vector <TVector3> acrossp_o;
    vector <TVector3> acrossp_e;
    vector <TVector3> extraordinary;
    vector <TVector3> Pt;
    vector <TVector3> Pt_extraordinary;
    vector <TVector3> Pt_ordinary;
    vector <TVector3> Pr1;

    vector <TVector3> Pr2;

    // These term1/term2/term3 arrays appear to be placeholders for decomposed
    // contributions to receiver projections or signal terms.
    vector <vector <double>> term1_Pr1;
    vector <vector <double>> term2_Pr1;
    vector <vector <double>> term3_Pr1;
    vector <vector <double>> term1_Pr2;
    vector <vector <double>> term2_Pr2;
    vector <vector <double>> term3_Pr2;

    // Same kind of decomposition, but for distance-based rather than depth-based scans
    vector <vector <double>> term1_Pr1_distances;
    vector <vector <double>> term2_Pr1_distances;
    vector <vector <double>> term3_Pr1_distances;
    vector <vector <double>> term1_Pr2_distances;
    vector <vector <double>> term2_Pr2_distances;
    vector <vector <double>> term3_Pr2_distances;

    TVector3 temp;

    // One entry per station
    rhat.resize(6);
    rhat_receive.resize(6);
    extraordinary.resize(6);
    p_o.resize(6);
    p_e.resize(6);
    p_o_rx.resize(6);
    p_e_rx.resize(6);
    acrossp_o.resize(6);
    acrossp_e.resize(6);
    Pt.resize(6);
    Pt_extraordinary.resize(6);
    Pt_ordinary.resize(6);
    Pr1.resize(6);

    Pr2.resize(6);
    term1_Pr1.resize(6);
    term2_Pr1.resize(6);
    term3_Pr1.resize(6);
    term1_Pr2.resize(6);
    term2_Pr2.resize(6);
    term3_Pr2.resize(6);

    term1_Pr1_distances.resize(6);
    term2_Pr1_distances.resize(6);
    term3_Pr1_distances.resize(6);
    term1_Pr2_distances.resize(6);
    term2_Pr2_distances.resize(6);
    term3_Pr2_distances.resize(6);


    // ------------------------------------------------------------------------
    // Fit / model functions and graphs used later for comparing the predicted
    // interference/attenuation envelopes to the measured data.
    // ------------------------------------------------------------------------
    TF1 *f1[NSTATIONS];

    TF1 *f1_nointerference[NSTATIONS];
    TF1 *f1_distances[NSTATIONS];
    TF1 *f1_nointerference_distances[NSTATIONS];
    TGraph *g1_distances[NSTATIONS];
    TGraph *g_atten_beam_distances[NSTATIONS];
    TGraph *g_atten_beam_crosspol_distances[NSTATIONS];
    TGraph *g_atten_beam_crosspol_nointerferencefunc_distances[NSTATIONS];
    TGraph *g_atten_beam_crosspol_func_distances[NSTATIONS];

    TGraph *g_nsolutions_distances[NSTATIONS];
    TGraph *g_sumphase_distances[NSTATIONS];

    // Diagnostic / output spectra are later written out for selected depths
    const int NSPECIALDEPTHS=2;

    double specialdepths[NSPECIALDEPTHS]={-850.,-875.};

    // g_spectra[station][depth_index] will hold the frequency spectrum
    vector< vector<TGraph*> > g_spectra;

    g_spectra.resize(NSTATIONS);
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
    TGraph *g_parameter0[NSTATIONS];
    TGraph *g_atten[NSTATIONS];
    TGraph *g_atten_power[NSTATIONS];

    TGraph *gfunc_noadjust[NSTATIONS];
    TGraph *g_atten_beam[NSTATIONS];
    TGraph *g_atten_beam_power[NSTATIONS];
    TGraph *g_atten_beam_crosspol[NSTATIONS];

    TGraph *g_atten_beam_crosspol_nointerferencefunc[NSTATIONS];
    TGraph *g_atten_beam_crosspol_func[NSTATIONS];
    TGraph *g_sumphase[NSTATIONS];
    TGraph *g_notflipped[NSTATIONS];
    TGraph *g_deltan[NSTATIONS];
    TGraph *g_notflipped_alongpath[NSTATIONS];
    TGraph *g_theta1_alongpath[NSTATIONS];
    TGraph *g_theta2_alongpath[NSTATIONS];
    TGraph *g_thetape1_alongpath[NSTATIONS];
    TGraph *g_thetape2_alongpath[NSTATIONS];
    TGraph *g_thetape1_phipe1_alongpath[NSTATIONS];
    TGraph *g_thetape2_phipe2_alongpath[NSTATIONS];
    TGraph *g_phipe1_alongpath[NSTATIONS];
    TGraph *g_phipe2_alongpath[NSTATIONS];
    TGraph *g_deltan_pulserdepth[NSTATIONS];
    TGraph *g_depth_istep[NSTATIONS];
    TGraph *g_path[NSTATIONS];
    TGraph *g_sumlength[NSTATIONS];
    TGraph *g_attenlengths[NSTATIONS];

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

    vmag_atten_beam_bigpic.resize(NSTATIONS);
    vmag_atten_beam_crosspol_bigpic.resize(NSTATIONS);
    vmag_atten_beam_crosspol_nointerferencefunc_bigpic.resize(NSTATIONS);
    vmag_atten_beam_crosspol_func_bigpic.resize(NSTATIONS);
    vtimediff_bigpic.resize(NSTATIONS);
    vspectrum_bigpic.resize(NSTATIONS);

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
    vistep.resize(NSTATIONS);
    vmag_parameter0.resize(NSTATIONS);
    vmag_atten.resize(NSTATIONS);
    vmag_atten_beam.resize(NSTATIONS);
    vmag_atten_beam_crosspol.resize(NSTATIONS);
    vmag_atten_beam_crosspol_nointerferencefunc.resize(NSTATIONS);
    vmag_atten_beam_crosspol_func.resize(NSTATIONS);
    vmag_func_noadjust.resize(NSTATIONS);
    vtimediff.resize(NSTATIONS);
    vnotflipped.resize(NSTATIONS);
    vsumlength.resize(NSTATIONS);
    vdeltan.resize(NSTATIONS);
    vnotflipped_alongpath.resize(NSTATIONS);
    vtheta1_alongpath.resize(NSTATIONS);
    vtheta2_alongpath.resize(NSTATIONS);
    vthetape1_alongpath.resize(NSTATIONS);
    vthetape2_alongpath.resize(NSTATIONS);
    vphipe1_alongpath.resize(NSTATIONS);
    vphipe2_alongpath.resize(NSTATIONS);
    vdepth_step.resize(NSTATIONS);
    vlengths.resize(NSTATIONS);
    vattenlengths.resize(NSTATIONS);
    vspectra.resize(NSTATIONS);
    vattens.resize(NSTATIONS);

    vdistances_bigpic.resize(NSTATIONS);
    vpseudodepths_bigpic.resize(NSTATIONS);

    int jmin[NSTATIONS];

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
    for (int i=minstation;i<=maxstation;i++) {
        vmag_atten_beam[i].resize(vdepth[i].size());
        vmag_atten_beam_crosspol[i].resize(vdepth[i].size());
        vmag_atten_beam_crosspol_nointerferencefunc[i].resize(vdepth[i].size());
        vmag_atten_beam_crosspol_func[i].resize(vdepth[i].size());
        vspectra[i].resize(vdepth[i].size());
        vattens[i].resize(vdepth[i].size());
        g_spectra[i].resize(vdepth[i].size());
        vmag_atten[i].resize(vdepth[i].size());
        vmag_parameter0[i].resize(vdepth[i].size());
        vmag_func_noadjust[i].resize(vdepth[i].size());

        // Measured data vectors only exist for A1-A5, not ARIANNA
        if (i!=5) {
            vtotal_distances[i].resize(vdepth_data[i].size());
            vtotal_distances_err[i].resize(vdepth_data[i].size());
        }
    
        // Smallest possible source-receiver separation is the horizontal distance
        double mindistance=horizontal_distances[i];
        jmin[i]=(int)(mindistance/STEP);

        // Build a distance scan and convert each distance into an equivalent
        // pseudo-depth for plotting/model-comparison purposes.
        for (int j=jmin[i]+1;j<NDISTANCES_BIGPIC;j++) {
            double thisdistance=STEP*(double)j;
            vdistances_bigpic[i].push_back(thisdistance);
            vpseudodepths_bigpic[i].push_back(station_depths[i]-1.*sqrt(thisdistance*thisdistance-mindistance*mindistance));

        }

        // TF1 wrappers used later to draw or fit the interference model
        f1[i]=new TF1("f1",sfunc.c_str(),-2000.,-100.);
        f1_nointerference[i]=new TF1("f1_nointerference",sfunc_nointerference.c_str(),-2000.,-100.);

        // --------------------------------------------------------------------
        // Loop over all pulser depths for this station
        // --------------------------------------------------------------------
        for (int idepth=0;idepth<vdepth[i].size();idepth++) {

            // Reduced 2D geometry for ray tracing:
            // x = horizontal separation, z = depth
            double posstation[2];
            posstation[0]=sqrt(pow(station_coords[i][0]-pulser_coords[0],2)+pow(station_coords[i][1]-pulser_coords[1],2));
            posstation[1]=station_depths[i];
            double pospulser[2];
            pospulser[0]=0.;
            pospulser[1]=vdepth[i][idepth];

            // Full 3D pulser position, useful for later geometric interpretation
            TVector3 pospulser3D;
            pospulser3D.SetX(pulser_coords[0]);
            pospulser3D.SetY(pulser_coords[1]);
            pospulser3D.SetZ(vdepth[i][idepth]);

            // 3D vector from pulser to station
            pulsertostation[i][0]=(station_coords[i][0]-pulser_coords[0]);
            pulsertostation[i][1]=(station_coords[i][1]-pulser_coords[1]);
            pulsertostation[i][2]=station_depths[i]-vdepth[i][idepth];

            // Legacy / diagnostic "special" geometry, based on station 0 and the first pulser depth
            double atten_special;
            double pospulser_special[2];
            pospulser_special[0]=0.;
            pospulser_special[1]=vdepth[0][0];
            double posstation_special[2];
            posstation_special[0]=sqrt(pow(station_coords[0][0]-pulser_coords[0],2)+pow(station_coords[0][1]-pulser_coords[1],2));
            posstation_special[1]=station_depths[0];

            // Ray-tracing endpoints
            double x0=0;
            double z0=pospulser[1];
            double x1=posstation[0];
            double z1=posstation[1];

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
            voutput6[i].push_back(getresults[6]);
            voutput7[i].push_back(getresults[7]);
            voutput8[i].push_back(getresults[8]);

            // ------------------------------------------------------------
            // Determine which ray type exists and build the full ray path
            // ------------------------------------------------------------
            if (getresults[6]!=-1000) {
                // Direct ray solution
                paramsd=GetDirectRayPar(z0,x1,z1);
                launch_angle=paramsd[1]/DEGRAD;
                receive_angle=paramsd[0]/DEGRAD;
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
                GetFullRefractedRayPath(z0, x1, z1, paramsra[7], paramsra[3], res, zs);

            }
            else if (getresults[7]!=1000) {

                // Reflected ray solution
                paramsre=GetReflectedRayPar(z0, x1 ,z1);
                double LangR=paramsre[1];
                double RangR=paramsre[0];
                launch_angle=paramsre[1]/DEGRAD;
                receive_angle=paramsre[0]/DEGRAD;
                GetFullReflectedRayPath(z0, x1, z1, LangR, res, zs);
            }

            // Running amplitude attenuation factor along the ray
            double atten=1.;

            // Frequency-dependent attenuation starts at unity for all frequencies
            vattens[i][idepth].clear();
            for (int ifreq=0;ifreq<NFREQ;ifreq++) {
                vattens[i][idepth].push_back(1.);
            }

            // Running path-length and birefringent phase accumulators
            double sumlength=0.;
            double sumphase=0.;

            // Proceed only if a usable ray type exists
            if (getresults[6]!=-1000 || getresults[8]!=-1000 || getresults[7]!=0) {

                // Horizontal unit vector toward the station.
                // Used to embed the 2D ray-tracing solution back into 3D.
                TVector3 yhat(station_coords[i][0]-pulser_coords[0],
                station_coords[i][1]-pulser_coords[1],
                0.);
                if (yhat.Mag()<HOWSMALLISTOOSMALL)
                cout << "yhat mag is " << yhat.Mag() << "\n";
                yhat.SetMag(1.);

                double angle_yhat=atan2(yhat[1],yhat[0]);


                vector<double> nvec_thisstep;
                nvec_thisstep.resize(3);

                // Principal axis at the start of the path
                nvec_thisstep[0]=gn1->Eval(zs[0]);
                nvec_thisstep[1]=gn2->Eval(zs[0]);
                nvec_thisstep[2]=gn3->Eval(zs[0]);

                TVector3 rhat_thisstep;

                // First estimate of local ray direction from the first path segment
                rhat_thisstep[0]=-1.*(res[UZAIRSTEP]-res[0])*yhat[0];
                rhat_thisstep[1]=-1.*(res[UZAIRSTEP]-res[0])*yhat[1];
                rhat_thisstep[2]=-1.*(zs[UZAIRSTEP]-zs[0]);

                if (rhat_thisstep.Mag()<1.E-8){
                    cout << "before calling getDeltaN at place 1, rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
                }

                // Solve local birefringence eigenvalues and eigenvectors for this propagation direction
                double deltan_alongpath=getDeltaN(BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1,p_e2);


                if (p_e2.Mag()<HOWSMALLISTOOSMALL)
                cout << "1, p_e2 is " << p_e2.Mag();

                // Integral of (+/- delta n)*ds along the ray
                double deltantimeslength_alongpath=0.;

                // Keep track of previous eigenvectors / flipping state
                TVector3 p_e1_previous=p_e1;
                TVector3 p_e2_previous=p_e2;
                double notflipped_previous=1.;
                double notflipped_atend=1.;
                double theta_e1_start=0.;
                double notflipped=1.;

                // --------------------------------------------------------
                // Step along the traced ray path in chunks of UZAIRSTEP
                // --------------------------------------------------------
                for (int istep=UZAIRSTEP;istep<res.size();istep+=UZAIRSTEP) {

                    nvec_thisstep.resize(3);

                    // Refractive indices at this path point
                    nvec_thisstep[0]=gn1->Eval(zs[istep]);
                    nvec_thisstep[1]=gn2->Eval(zs[istep]);
                    nvec_thisstep[2]=gn3->Eval(zs[istep]);

                    if (istep>0) {

                        // Local propagation direction reconstructed from the discrete path
                        rhat_thisstep[0]=-1.*(res[istep]-res[istep-UZAIRSTEP])*yhat[0];
                        rhat_thisstep[1]=-1.*(res[istep]-res[istep-UZAIRSTEP])*yhat[1];
                        rhat_thisstep[2]=-1.*(zs[istep]-zs[istep-UZAIRSTEP]);

                        // Debugging / illustration printout for A1 near -1000 m
                        if (i==0 && idepth==g_idepth[i]->Eval(-1000.)) {
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

                                cout << "launch angle is " << rhat_thisstep.Theta()*DEGRAD << "\n";
                                cout << "launch vector is " << res[istep-50]-res[istep] << "\t" << zs[istep-50]-zs[istep] << "\n";
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
                        deltan_alongpath=getDeltaN(BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1,p_e2);

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
                        TVector3 E_e1=rotateD(epsilon_thisstep,angle_iceflow,p_e1);
                        TVector3 E_e2=rotateD(epsilon_thisstep,angle_iceflow,p_e2);

                        // ----------------------------------------------------
                        // At the final step (TX side in this path convention),
                        // compute transmitter-side polarization / beam quantities
                        // ----------------------------------------------------
                        if (abs((double)(istep-(int)res.size()))<=UZAIRSTEP) {


                            if (i==5 && idepth==g_idepth[i]->Eval(-400.))
                            cout << "launch angle is " << rhat_thisstep.Theta()*DEGRAD << "\n";

                            double theta_e1,theta_e2;
                            double thetaE_e1,thetaE_e2;
                            double theta_e1_Sclock,theta_e2_Sclock;
                            double thetaE_e1_Sclock,thetaE_e2_Sclock;

                            TVector3 Shat_e1,Shat_e2;

                            double E_e1_thetacomponent,E_e2_thetacomponent;
                            double E_e1_phicomponent,E_e2_phicomponent;

                            // Decompose the eigenmodes into the chosen TX cross-pol frame
                            getManyAnglesontheClock(BIAXIAL,CROSSPOLANGLE_TX,
                                rhat_thisstep,
                                p_e1,p_e2,E_e1,E_e2,
                                theta_e1,theta_e2,thetaE_e1,thetaE_e2,
                                theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock,
                                Shat_e1,Shat_e2,
                                E_e1_thetacomponent,E_e2_thetacomponent,
                                E_e1_phicomponent,E_e2_phicomponent);

                            // Store angle between each Poynting vector and propagation direction
                            vangle_Shat_e1_khat[i].push_back(acos(Shat_e1.Dot(rhat_thisstep)/Shat_e1.Mag()/rhat_thisstep.Mag())*DEGRAD);
                            vangle_Shat_e2_khat[i].push_back(acos(Shat_e2.Dot(rhat_thisstep)/Shat_e2.Mag()/rhat_thisstep.Mag())*DEGRAD);

                            // Special ordering swap for ARIANNA in the biaxial case
                            if (i==5 && BIAXIAL==1) {
                                switchThem(thetaE_e1_Sclock,thetaE_e2_Sclock);
                                switchThem(theta_e1_Sclock,theta_e2_Sclock);
                                switchThem(theta_e1,theta_e2);
                                switchThem(thetaE_e1,thetaE_e2);
                            }

                            // Verbose diagnostics for selected reference depths
                            if (i==0 && idepth==g_idepth[i]->Eval(-1000.) ||
                            i==5 && idepth==g_idepth[i]->Eval(-400.)) {

                                cout << "At Tx:\n";
                                cout << "A" << i+1 << ", depth is " << vdepth[i][idepth] << "\n";

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
                                cout << "depth is " << vdepth[i][idepth] << "\n";
                                cout << "theta of rhat_thisstep is " << rhat_thisstep.Theta()*DEGRAD << "\n";
                                cout << "rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
                                cout << "Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
                                cout << "Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";

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
                            }

                            // Store orthogonality / alignment diagnostics
                            vdotShats_tx[i].push_back(Shat_e1.Dot(Shat_e2)/Shat_e1.Mag()/Shat_e2.Mag());
                            vdotEhats_tx[i].push_back(E_e1.Dot(E_e2)/E_e1.Mag()/E_e2.Mag());
                            vdotDhats_tx[i].push_back(p_e1.Dot(p_e2)/p_e1.Mag()/p_e2.Mag());

                            if (i==0 && idepth==g_idepth[i]->Eval(-1000.) ||
                            i==5 && idepth==g_idepth[i]->Eval(-400.)) {
                                cout << "Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
                                cout << "Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";
                                cout << "theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";
                            }

                            // Beam factors at the transmitter from the Poynting-vector polar angle
                            double beam_tx_e1=sin(Shat_e1.Theta());
                            double beam_tx_e2=sin(Shat_e2.Theta());

                            vtxdepth_beam1[i].push_back(beam_tx_e1);
                            vtxdepth_beam2[i].push_back(beam_tx_e2);

                            // Store TX angular diagnostics
                            vtxdepth_theta1[i].push_back(theta_e1*DEGRAD);
                            vtxdepth_theta2[i].push_back(theta_e2*DEGRAD);

                            vtxdepth_theta1_Sclock[i].push_back(theta_e1_Sclock*DEGRAD);
                            vtxdepth_theta2_Sclock[i].push_back(theta_e2_Sclock*DEGRAD);

                            vtxdepthE_theta1[i].push_back(thetaE_e1*DEGRAD);
                            vtxdepthE_theta2[i].push_back(thetaE_e2*DEGRAD);


                            vtxdepthE_theta1_Sclock[i].push_back(thetaE_e1_Sclock*DEGRAD);
                            vtxdepthE_theta2_Sclock[i].push_back(thetaE_e2_Sclock*DEGRAD);

                            // Angle between E and D for each eigenmode
                            vtxdepth_dispersion1[i].push_back(acos(E_e1.Dot(p_e1)/E_e1.Mag()/p_e1.Mag())*DEGRAD);
                            vtxdepth_dispersion2[i].push_back(acos(E_e2.Dot(p_e2)/E_e2.Mag()/p_e2.Mag())*DEGRAD);

                            // Convert S-clock angles into epsilon parameters
                            double epsilon1_tx=0.;
                            double epsilon2_tx=0.;
                            thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,
                            epsilon1_tx,epsilon2_tx);

                            if (epsilon2_tx>PI/2.){
                                epsilon2_tx-=PI;
                            }

                            vepsilon1_tx[i].push_back(epsilon1_tx*DEGRAD);
                            vepsilon2_tx[i].push_back(epsilon2_tx*DEGRAD);

                            vdiffepsilon_tx[i].push_back((epsilon2_tx-epsilon1_tx)*DEGRAD);

                            // Store attenuation-only and attenuation×beam estimates at the receiver
                            vrxdepth_atten[i].push_back(VOLTAGENORM*atten);
                            vrxdepth_atten_beam[i].push_back(VOLTAGENORM*vrxdepth_beam1[i][idepth]*vtxdepth_beam1[i][idepth]*atten);
                            vrxdepth_atten_power[i].push_back(VOLTAGENORM*VOLTAGENORM*atten*atten);
                            vrxdepth_atten_beam_power[i].push_back(VOLTAGENORM*VOLTAGENORM*vrxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]*vtxdepth_beam1[i][idepth]*vtxdepth_beam1[i][idepth]*atten*atten);
                            
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
                            getManyAnglesontheClock(BIAXIAL,CROSSPOLANGLE_RX,
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
                            if (i==0 && idepth==g_idepth[i]->Eval(-1000.) ||
                            i==5 && idepth==g_idepth[i]->Eval(-400.)) {
                                cout << "At Rx:\n";
                                cout << "A" << i+1 << ", depth is " << vdepth[i][idepth] << "\n";
                                cout << "thetas on the Sclock are " << thetaE_e1_Sclock*DEGRAD << "\t" << thetaE_e2_Sclock*DEGRAD << "\n";
                                cout << "depth is " << vdepth[i][idepth] << "\n";

                                cout << "theta of rhat_thisstep is " << rhat_thisstep.Theta()*DEGRAD << "\n";
                                cout << "rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
                                cout << "This points from pulser to station: " << station_coords[i][0]-pulser_coords[0] << "\t" << station_coords[i][1]-pulser_coords[1] << "\t" << station_depths[i]-vdepth[i][idepth] << "\n";

                                cout << "Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
                                cout << "Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";
                                cout << "theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";
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
                            vrxdepth_theta1[i].push_back(theta_e1*DEGRAD);
                            vrxdepth_theta2[i].push_back(theta_e2*DEGRAD);

                            // Note: these two lines store theta_e1/theta_e2 rather than thetaE_e1/thetaE_e2.
                            // That may be intentional or may deserve a later check.
                            vrxdepthE_theta1[i].push_back(theta_e1*DEGRAD);
                            vrxdepthE_theta2[i].push_back(theta_e2*DEGRAD);

                            vrxdepth_theta1_Sclock[i].push_back(theta_e1_Sclock*DEGRAD);
                            vrxdepth_theta2_Sclock[i].push_back(theta_e2_Sclock*DEGRAD);


                            vrxdepthE_theta1_Sclock[i].push_back(thetaE_e1_Sclock*DEGRAD);
                            vrxdepthE_theta2_Sclock[i].push_back(thetaE_e2_Sclock*DEGRAD);

                            // Convert RX S-clock angles into epsilon parameters
                            double epsilon1_rx=0.;
                            double epsilon2_rx=0.;

                            thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,
                                epsilon1_rx,epsilon2_rx);

                            vepsilon1_rx[i].push_back(epsilon1_rx*DEGRAD);
                            vepsilon2_rx[i].push_back(epsilon2_rx*DEGRAD);
                            vdiffepsilon_rx[i].push_back((epsilon2_rx-epsilon1_rx)*DEGRAD);

                            // Beam factors at the receiver
                            double beam_rx_e1=sin(Shat_e1.Theta());
                            double beam_rx_e2=sin(Shat_e2.Theta());

                            vrxdepth_beam1[i].push_back(beam_rx_e1);
                            vrxdepth_beam2[i].push_back(beam_rx_e2);

                            // Build an orthonormal receiver polarization basis:
                            //   Pr2 perpendicular to Shat_e1 and +z
                            //   Pr1 perpendicular to both Pr2 and Shat_e1
                            TVector3 plusz(0.,0.,1.);
                            Pr2[i]=Shat_e1.Cross(plusz);
                            if (Pr2[i].Mag()<HOWSMALLISTOOSMALL){
                                cout << "Pr2[i] is " << Pr2[i].Mag() << "\n";
                            }
                            Pr2[i].SetMag(1.);
                            Pr1[i]=Pr2[i].Cross(Shat_e1);
                            if (Pr1[i].Mag()<HOWSMALLISTOOSMALL){
                                cout << "Pr1[i] is " << Pr1[i].Mag() << "\n";
                            }
                            Pr1[i].SetMag(1.);
                        }

                        // ----------------------------------------------------
                        // Update the frequency-dependent attenuation spectrum
                        // along this path segment
                        // ----------------------------------------------------
                        for (int ifreq=0;ifreq<NFREQ;ifreq++) {


                            double this_atten_length=GetIceAttenuationLength(zs[istep], vfreqs[ifreq]/1.E9);

                            // Power-like attenuation: two exponential factors
                            vattens[i][idepth][ifreq] = vattens[i][idepth][ifreq]*exp(-1.*length/this_atten_length)*exp(-1.*length/this_atten_length);
                        }

                        // Recompute angular decomposition in the "neutral" clock convention
                        double theta_e1,theta_e2;
                        double thetaE_e1,thetaE_e2;
                        double theta_e1_Sclock,theta_e2_Sclock;
                        double thetaE_e1_Sclock,thetaE_e2_Sclock;

                        TVector3 Shat_e1,Shat_e2;

                        double E_e1_thetacomponent,E_e2_thetacomponent;
                        double E_e1_phicomponent,E_e2_phicomponent;
                        getManyAnglesontheClock(BIAXIAL,0.,
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
                        if (idepth==(int)(g_idepth[i]->Eval(depth_special))) {

                            vtheta1_alongpath[i].push_back(theta_e1*DEGRAD);
                            vtheta2_alongpath[i].push_back(theta_e2*DEGRAD);
                            vthetape1_alongpath[i].push_back(p_e1.Theta()*DEGRAD);
                            vthetape2_alongpath[i].push_back(p_e2.Theta()*DEGRAD);
                            vphipe1_alongpath[i].push_back(p_e1.Phi()*DEGRAD);
                            vphipe2_alongpath[i].push_back(p_e2.Phi()*DEGRAD);
                            vnotflipped_alongpath[i].push_back(notflipped);
                            vdeltan[i].push_back(deltan_alongpath);
                            vdepth_step[i].push_back(zs[istep]);
                            vistep[i].push_back((double)istep);
                            vlengths[i].push_back(sumlength);
                            vattenlengths[i].push_back(atten_length);
                        }
                    }
                    else {                        
                        // Fallback for the first step if needed:
                        // use the previously stored launch-direction unit vector.
                        rhat_thisstep=rhat[i];
                    }
                }

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
                theta1_Sclock_atrx=vrxdepthE_theta1_Sclock[i][idepth];
                theta2_Sclock_atrx=vrxdepthE_theta2_Sclock[i][idepth];


                vV1_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
                vE1_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]);

                // Mode 2 contribution into receiver channel r1
                vV2_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);
                vE2_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]);

                // LPDA-like projection: remove the component parallel to the beam axis
                vV1_r1_lpda[i].push_back(vV1_r1[i][idepth]*sqrt(1.-vrxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]));
                vV2_r1_lpda[i].push_back(vV2_r1[i][idepth]*sqrt(1.-vrxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]));

                // Store quadratic combinations for interference calculations
                vV1squared_r1[i].push_back(vV1_r1[i][idepth]*vV1_r1[i][idepth]);
                vV2squared_r1[i].push_back(vV2_r1[i][idepth]*vV2_r1[i][idepth]);
                vV1V2_r1[i].push_back(vV1_r1[i][idepth]*vV2_r1[i][idepth]);

                // Mode 1 contribution into receiver channel r2
                vV1_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
                vE1_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]);

                // Mode 2 contribution into receiver channel r2
                vV2_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);
                vE2_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]);

                // More quadratic combinations for the second receiver channel
                vV1squared_r2[i].push_back(vV1_r2[i][idepth]*vV1_r2[i][idepth]);
                vV2squared_r2[i].push_back(vV2_r2[i][idepth]*vV2_r2[i][idepth]);
                vV1V2_r2[i].push_back(vV1_r2[i][idepth]*vV2_r2[i][idepth]);

                // Debug print for a reference pulser depth
                if (idepth==g_idepth[i]->Eval(-1000.)){
                    cout << "station, depth, V1squared_r2, V2squared_r2, V1V2_r2 are " << i << "\t" << vV1squared_r2[i][idepth] << "\t" << vV2squared_r2[i][idepth] << "\t" << vV1V2_r2[i][idepth] << "\n";
                }

                // Negative cross terms, convenient for expressions written as
                // A+B-2*sqrt(AB)*sin^2(phi) or equivalent forms
                voppositeV1V2_r2[i].push_back(-1.*vV1_r2[i][idepth]*vV2_r2[i][idepth]);
                voppositeV1V2_r1[i].push_back(-1.*vV1_r1[i][idepth]*vV2_r1[i][idepth]);


                // ----------------------------------------------------------------
                // Build envelope quantities:
                //   minus envelope = destructive combination
                //   plus envelope  = constructive combination
                //
                // Versions are stored both for voltage-like quantities and
                // field/Poynting-like quantities.
                // ----------------------------------------------------------------
                venvelope_minus_r1[i].push_back((vV1_r1[i][idepth]-vV2_r1[i][idepth])*(vV1_r1[i][idepth]-vV2_r1[i][idepth]));
                vSenvelope_minus_r1[i].push_back((vE1_r1[i][idepth]-vE2_r1[i][idepth])*(vE1_r1[i][idepth]-vE2_r1[i][idepth]));

                venvelope_minus_r1_lpda[i].push_back((vV1_r1_lpda[i][idepth]-vV2_r1_lpda[i][idepth])*(vV1_r1_lpda[i][idepth]-vV2_r1_lpda[i][idepth]));

                venvelope_minus_r2[i].push_back((vV1_r2[i][idepth]-vV2_r2[i][idepth])*(vV1_r2[i][idepth]-vV2_r2[i][idepth]));
                vSenvelope_minus_r2[i].push_back((vE1_r2[i][idepth]-vE2_r2[i][idepth])*(vE1_r2[i][idepth]-vE2_r2[i][idepth]));
                venvelope_plus_r1[i].push_back((vV2_r1[i][idepth]+vV1_r1[i][idepth])*(vV1_r1[i][idepth]+vV2_r1[i][idepth]));
                vSenvelope_plus_r1[i].push_back((vE2_r1[i][idepth]+vE1_r1[i][idepth])*(vE1_r1[i][idepth]+vE2_r1[i][idepth]));
                venvelope_plus_r1_lpda[i].push_back((vV2_r1_lpda[i][idepth]+vV1_r1_lpda[i][idepth])*(vV1_r1_lpda[i][idepth]+vV2_r1_lpda[i][idepth]));

                venvelope_plus_r2[i].push_back((vV2_r2[i][idepth]+vV1_r2[i][idepth])*(vV1_r2[i][idepth]+vV2_r2[i][idepth]));
                vSenvelope_plus_r2[i].push_back((vE2_r2[i][idepth]+vE1_r2[i][idepth])*(vE1_r2[i][idepth]+vE2_r2[i][idepth]));

                vvenvelope_minus_r1[i].push_back(sqrt((vV1_r1[i][idepth]-vV2_r1[i][idepth])*(vV1_r1[i][idepth]-vV2_r1[i][idepth])));
                vvenvelope_minus_r2[i].push_back(sqrt((vV1_r2[i][idepth]-vV2_r2[i][idepth])*(vV1_r2[i][idepth]-vV2_r2[i][idepth])));
                vvenvelope_plus_r1[i].push_back(sqrt((vV2_r1[i][idepth]+vV1_r1[i][idepth])*(vV1_r1[i][idepth]+vV2_r1[i][idepth])));
                vvenvelope_plus_r2[i].push_back(sqrt((vV2_r2[i][idepth]+vV1_r2[i][idepth])*(vV1_r2[i][idepth]+vV2_r2[i][idepth])));

                vEenvelope_minus_r1[i].push_back(sqrt((vE1_r1[i][idepth]-vE2_r1[i][idepth])*(vE1_r1[i][idepth]-vE2_r1[i][idepth])));
                vEenvelope_minus_r2[i].push_back(sqrt((vE1_r2[i][idepth]-vE2_r2[i][idepth])*(vE1_r2[i][idepth]-vE2_r2[i][idepth])));
                vEenvelope_plus_r1[i].push_back(sqrt((vE2_r1[i][idepth]+vE1_r1[i][idepth])*(vE1_r1[i][idepth]+vE2_r1[i][idepth])));
                vEenvelope_plus_r2[i].push_back(sqrt((vE2_r2[i][idepth]+vE1_r2[i][idepth])*(vE1_r2[i][idepth]+vE2_r2[i][idepth])));

                // ----------------------------------------------------------------
                // Build a frequency spectrum for this station and pulser depth.
                //
                // The phase scales linearly with frequency, so a spectrum is obtained
                // by re-evaluating the same interference expression at each frequency.
                // ----------------------------------------------------------------
                for (int ifreq=0;ifreq<NFREQ;ifreq++) {
                    double thisfreq=vfreqs[ifreq];
                    double thissumphase=sumphase*thisfreq/freq;

                    vspectra[i][idepth].push_back(vattens[i][idepth][ifreq]/vrxdepth_atten[i][idepth]*(venvelope_plus_r1[i][idepth]-4*vV1_r1[i][idepth]*vV2_r1[i][idepth]*sin(thissumphase)*sin(thissumphase)));
                    if (i==5 && idepth==g_idepth[i]->Eval(-1000.))
                    cout << "i, vattens[i][idepth], vrxdepth_atten[i][idepth], vspectra are " << i << "\t" << vattens[i][idepth][ifreq] << "\t" << vrxdepth_atten[i][idepth] << "\t" << vspectra[i][idepth][ifreq] << "\n";
                }

                // Time delay corresponding to the accumulated phase difference
                //   sumphase = pi * f * dt
                // so dt = sumphase / (pi*f)
                vtimediff[i].push_back(sumphase/(PI)*1./freq*1.E9);

                // Whether the eigenvector tracking ended in the same orientation sign
                vnotflipped[i].push_back(notflipped_atend);

                // ----------------------------------------------------------------
                // Final interference-modified power / Poynting-like quantities
                // for the two receiver channels.
                // ----------------------------------------------------------------
                vpower_r1[i].push_back(venvelope_plus_r1[i][idepth]-4*vV1_r1[i][idepth]*vV2_r1[i][idepth]*sin(sumphase)*sin(sumphase));
                vpoynting_r1[i].push_back(vSenvelope_plus_r1[i][idepth]-4*vE1_r1[i][idepth]*vE2_r1[i][idepth]*sin(sumphase)*sin(sumphase));

                vpower_r1_lpda[i].push_back(venvelope_plus_r1_lpda[i][idepth]-4*vV1_r1_lpda[i][idepth]*vV2_r1_lpda[i][idepth]*sin(sumphase)*sin(sumphase));
                vpower_r2[i].push_back(venvelope_plus_r2[i][idepth]-4*vV1_r2[i][idepth]*vV2_r2[i][idepth]*sin(sumphase)*sin(sumphase));
                vpoynting_r2[i].push_back(vSenvelope_plus_r2[i][idepth]-4*vE1_r2[i][idepth]*vE2_r2[i][idepth]*sin(sumphase)*sin(sumphase));

                // Convert power-like quantities into amplitude-like quantities
                vvoltage_r1[i].push_back(sqrt(vpower_r1[i][idepth]));
                vfield_r1[i].push_back(sqrt(vpoynting_r1[i][idepth]));
                vvoltage_r1_lpda[i].push_back(sqrt(vpower_r1_lpda[i][idepth]));
                vvoltage_r2[i].push_back(sqrt(vpower_r2[i][idepth]));
                vfield_r2[i].push_back(sqrt(vpoynting_r2[i][idepth]));

                // Detailed printout at selected reference depths
                if (i==0 && idepth==g_idepth[i]->Eval(-1000.) ||
                i==5 && idepth==g_idepth[i]->Eval(-400.)) {

                    cout << "A " << i+1 << ", depth is " << vdepth[i][idepth] << "\n";

                    cout << "before scalefactors:\n";
                    cout << "powers are " << vpower_r1[i][idepth] << "\t" << vpower_r2[i][idepth] << "\n";
                    cout << "voltages are " << vvoltage_r1[i][idepth] << "\t" << vvoltage_r2[i][idepth] << "\n";
                    cout << "lpda voltages are " << vvoltage_r1_lpda[i][idepth] << "\t" << vvoltage_r2[i][idepth] << "\n";

                    // Scale to a convenient display size
                    double scalefactor=0.2/sqrt(vvoltage_r1[i][idepth]*vvoltage_r1[i][idepth]+vvoltage_r2[i][idepth]*vvoltage_r2[i][idepth]);

                    // Common attenuation × beam prefactor
                    double prefactor=vrxdepth_atten[i][idepth]*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth];
                    cout << "voltage r1 is " << scalefactor*vvoltage_r1[i][idepth] << "\n";
                    cout << "x, z components of voltage r1 are " << vtxdepth_beam1[i][idepth]*scalefactor*vvoltage_r1[i][idepth] << "\t" << sqrt(1.-vtxdepth_beam1[i][idepth]*vtxdepth_beam1[i][idepth])*scalefactor*vvoltage_r1[i][idepth] << "\n";

                    cout << "voltage r2 is " << scalefactor*vvoltage_r2[i][idepth] << "\n";
                    cout << "polarization angle is " << DEGRAD*atan2(vvoltage_r2[i][idepth],vvoltage_r1[i][idepth]) << "\n";
                    cout << "thetas_Sclock_atrx are " << theta1_Sclock_atrx << "\t" << theta2_Sclock_atrx << "\n";
                    cout << "diff is " << (theta1_Sclock_atrx-theta2_Sclock_atrx) << "\n";
                    cout << "thetas_Sclock_attx are " << vtxdepthE_theta1_Sclock[i][idepth] << "\t" << vtxdepthE_theta2_Sclock[i][idepth] << "\n";
                    cout << "diff is " << (vtxdepthE_theta1_Sclock[i][idepth]-vtxdepthE_theta2_Sclock[i][idepth]) << "\n";

                    cout << "vV1_r1, vV2_r1 are " << vV1_r1[i][idepth]/prefactor << "\t" << vV2_r1[i][idepth]/prefactor << "\n";
                    cout << "vV1_r2, vV2_r2 are " << vV1_r2[i][idepth]/prefactor << "\t" << vV2_r2[i][idepth]/prefactor << "\n";

                    cout << "ray 1 at tx is " << scalefactor*vrxdepth_atten[i][idepth]*vtxdepth_beam1[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD) << "\n";
                    cout << "ray 2 at tx is " << scalefactor*vrxdepth_atten[i][idepth]*vtxdepth_beam1[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD) << "\n";
                    cout << "fraction of wavelength is " << sumphase/(2.*PI) << "\n";
                    cout << "vpower_r1+vpower_r2 \t" <<vpower_r1[i][idepth]+vpower_r2[i][idepth] << "\n";
                    cout << "prefactor is " << pow(vrxdepth_atten[i][idepth]*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth],2) << "\n";
                    cout << "three factors are " << vrxdepth_atten[i][idepth] << "\t" << vtxdepth_beam1[i][idepth] << "\t" << vrxdepth_beam1[i][idepth] << "\n";


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
                vpolarization_Psi_rx[i].push_back(atan2(vvoltage_r2[i][idepth],vvoltage_r1[i][idepth])*DEGRAD);
                vpolarization_Omega_rx[i].push_back(acos(vvoltage_r1[i][idepth]/(sqrt(vvoltage_r1[i][idepth]*vvoltage_r1[i][idepth]+vvoltage_r2[i][idepth]*vvoltage_r2[i][idepth])))*DEGRAD);
                vEpolarization_Psi_rx[i].push_back(atan2(vfield_r2[i][idepth],vfield_r1[i][idepth])*DEGRAD);
                vEpolarization_Omega_rx[i].push_back(acos(vfield_r1[i][idepth]/(sqrt(vfield_r1[i][idepth]*vfield_r1[i][idepth]+vfield_r2[i][idepth]*vfield_r2[i][idepth])))*DEGRAD);

                // Store total geometric path length
                vsumlength[i].push_back(sumlength);

                // Refractive indices at the pulser depth itself
                vector<double> nvec_thisdepth;
                nvec_thisdepth.resize(3);
                nvec_thisdepth[0]=gn1->Eval(vdepth[i][idepth]);
                nvec_thisdepth[1]=gn2->Eval(vdepth[i][idepth]);
                nvec_thisdepth[2]=gn3->Eval(vdepth[i][idepth]);

                // ----------------------------------------------------------------
                // Build launch and receive direction unit vectors either from
                // ray tracing or from simple straight-line geometry.
                // ----------------------------------------------------------------
                if (DORAYTRACING) {

                    TVector3 plusz(0.,0.,1.);

                    // Rotation axis that brings +z into the vertical plane of the ray
                    TVector3 vrotate=plusz.Cross(yhat);

                    if (vrotate.Mag()<HOWSMALLISTOOSMALL){
                        cout << "vrotate mag is " << vrotate.Mag() << "\n";
                    }

                    // Launch direction
                    vrotate.SetMag(1.);
                    rhat[i]=plusz;
                    rhat[i].Rotate(launch_angle,vrotate);

                    // Receive direction
                    rhat_receive[i]=plusz;
                    rhat_receive[i].Rotate(receive_angle,vrotate);

                }
                else {

                    // Straight-line approximation if ray tracing is disabled
                    rhat[i].SetX(station_coords[i][0]-pulser_coords[0]);
                    rhat[i].SetY(station_coords[i][1]-pulser_coords[1]);
                    rhat[i].SetZ(station_depths[i]-vdepth[i][idepth]);

                    rhat_receive[i].SetX(rhat[i][0]);
                    rhat_receive[i].SetY(rhat[i][1]);
                    rhat_receive[i].SetZ(rhat[i][2]);
                }

                // Normalize and sanity check
                if (rhat[i].Mag()<HOWSMALLISTOOSMALL){
                    cout << "rhat[i] mag is " << rhat[i].Mag() << "\n";
                }
                rhat[i].SetMag(1.);

                if (rhat_receive[i].Mag()<HOWSMALLISTOOSMALL){
                    cout << "rhat_receive[i] mag is " << rhat_receive[i].Mag() << "\n";
                }
                rhat_receive[i].SetMag(1.);

            }
            else {

                // ----------------------------------------------------------------
                // If no ray solution was found, fall back to simple straight-line
                // source-to-station geometry.
                // ----------------------------------------------------------------
                rhat[i].SetX(station_coords[i][0]-pulser_coords[0]);
                rhat[i].SetY(station_coords[i][1]-pulser_coords[1]);
                rhat[i].SetZ(station_depths[i]-vdepth[i][idepth]);

                rhat_receive[i].SetX(rhat[i][0]);
                rhat_receive[i].SetY(rhat[i][1]);
                rhat_receive[i].SetZ(rhat[i][2]);
            }

            // Store launch and receive polar angles (degrees)
            vreceiveangle[i].push_back(rhat_receive[i].Theta()*DEGRAD);
            vlaunchangle[i].push_back(rhat[i].Theta()*DEGRAD);

            // Final normalization / safety checks
            if (rhat[i].Mag()<HOWSMALLISTOOSMALL){
                cout << "rhat[i] mag is " << rhat[i].Mag() << "\n";
            }

            rhat[i].SetMag(1.);

            if (rhat_receive[i].Mag()<HOWSMALLISTOOSMALL){
                cout << "rhat_receive[i] mag is " << rhat_receive[i].Mag() << "\n";
            }

            rhat_receive[i].SetMag(1.);
        }

        // --------------------------------------------------------------------
        // Convert all accumulated per-depth vectors into ROOT TGraphs for plotting
        // and for later writing to output files.
        // --------------------------------------------------------------------
        gtxdepth_beam1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_beam1[i][0]);
        gtxdepth_beam2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_beam2[i][0]);

        grxdepth_beam1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_beam1[i][0]);
        grxdepth_beam2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_beam2[i][0]);

        grxdepth_atten[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_atten[i][0]);
        grxdepth_atten_beam[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_atten_beam[i][0]);
        grxdepth_atten_power[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_atten_power[i][0]);
        grxdepth_atten_beam_power[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_atten_beam_power[i][0]);

        gtxdepth_theta1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_theta1[i][0]);
        gtxdepth_theta2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_theta2[i][0]);

        grxdepth_theta1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_theta1[i][0]);
        grxdepth_theta2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_theta2[i][0]);

        grxdepthE_theta1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepthE_theta1[i][0]);
        grxdepthE_theta2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepthE_theta2[i][0]);

        gtxdepth_theta1_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_theta1_Sclock[i][0]);
        gtxdepth_theta2_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_theta2_Sclock[i][0]);

        grxdepth_theta1_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_theta1_Sclock[i][0]);
        grxdepth_theta2_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_theta2_Sclock[i][0]);

        gtxdepth_dispersion1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_dispersion1[i][0]);
        gtxdepth_dispersion2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_dispersion2[i][0]);

        // Note: these two graphs use vtxdepth_theta1/2 rather than vtxdepthE_theta1/2
        gtxdepthE_theta1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_theta1[i][0]);
        gtxdepthE_theta2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_theta2[i][0]);

        gtxdepthE_theta1_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepthE_theta1_Sclock[i][0]);
        gtxdepthE_theta2_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepthE_theta2_Sclock[i][0]);

        grxdepthE_theta1_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepthE_theta1_Sclock[i][0]);
        grxdepthE_theta2_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepthE_theta2_Sclock[i][0]);

        gdotShats_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vdotShats_tx[i][0]);
        gdotEhats_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vdotEhats_tx[i][0]);
        gdotDhats_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vdotDhats_tx[i][0]);

        g_parameter0[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vmag_parameter0[i][0]);
        g_atten[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vmag_atten[i][0]);
        g_atten_beam[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vmag_atten_beam[i][0]);
        g_atten_beam_crosspol[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vmag_atten_beam_crosspol[i][0]);
        gfunc_noadjust[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vmag_func_noadjust[i][0]);
        g_atten_beam_crosspol_nointerferencefunc[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vmag_atten_beam_crosspol_nointerferencefunc[i][0]);
        g_atten_beam_crosspol_func[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vmag_atten_beam_crosspol_func[i][0]);
        g_sumphase[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtimediff[i][0]);
        g_notflipped[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vnotflipped[i][0]);
        g_sumlength[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vsumlength[i][0]);

        // One frequency spectrum graph per pulser depth
        for (int idepth=0;idepth<vdepth[i].size();idepth++) {
            g_spectra[i][idepth]=new TGraph(vfreqs.size(),&vfreqs[0],&vspectra[i][idepth][0]);
        }

        g_atten_power[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_atten_power[i][0]);
        g_atten_beam_power[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_atten_beam_power[i][0]);

        // Voltage / field / interference component graphs
        gV1_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV1_r1[i][0]);
        gV1_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV1_r2[i][0]);
        gV1squared_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV1squared_r1[i][0]);
        gV1squared_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV1squared_r2[i][0]);
        gpower_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vpower_r1[i][0]);
        gpower_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vpower_r2[i][0]);
        gvoltage_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vvoltage_r1[i][0]);
        gvoltage_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vvoltage_r2[i][0]);
        gfield_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vfield_r1[i][0]);
        gfield_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vfield_r2[i][0]);
        gV2_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV2_r1[i][0]);
        gV2_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV2_r2[i][0]);
        gV1V2_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV1V2_r1[i][0]);
        gV1V2_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV1V2_r2[i][0]);
        goppositeV1V2_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&voppositeV1V2_r1[i][0]);
        goppositeV1V2_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&voppositeV1V2_r2[i][0]);
        gV2squared_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV2squared_r1[i][0]);
        gV2squared_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vV2squared_r2[i][0]);
        genvelope_minus_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&venvelope_minus_r1[i][0]);
        genvelope_minus_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&venvelope_minus_r2[i][0]);
        genvelope_plus_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&venvelope_plus_r1[i][0]);
        genvelope_plus_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&venvelope_plus_r2[i][0]);
        gvenvelope_minus_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vvenvelope_minus_r1[i][0]);
        gvenvelope_minus_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vvenvelope_minus_r2[i][0]);
        gvenvelope_plus_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vvenvelope_plus_r1[i][0]);
        gvenvelope_plus_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vvenvelope_plus_r2[i][0]);
        gEenvelope_minus_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vEenvelope_minus_r1[i][0]);
        gEenvelope_minus_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vEenvelope_minus_r2[i][0]);
        gEenvelope_plus_r1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vEenvelope_plus_r1[i][0]);
        gEenvelope_plus_r2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vEenvelope_plus_r2[i][0]);
        gepsilon1_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vepsilon1_tx[i][0]);
        gepsilon2_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vepsilon2_tx[i][0]);

        // Polarization angle graphs
        gpolarization_Omega_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vpolarization_Omega_rx[i][0]);
        gpolarization_Psi_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vpolarization_Psi_rx[i][0]);

        gEpolarization_Omega_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vEpolarization_Omega_rx[i][0]);
        gEpolarization_Psi_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vEpolarization_Psi_rx[i][0]);

        // Build reversed-depth versions for plots that want positive pulser depth
        for (int j=0;j<vdepth[i].size();j++) {
            vpolarization_reversedepth_Omega_rx[i].push_back(vpolarization_Omega_rx[i][vpolarization_Omega_rx[i].size()-1-j]);
            vpolarization_reversedepth_Psi_rx[i].push_back(vpolarization_Psi_rx[i][vpolarization_Psi_rx[i].size()-1-j]);
            vEpolarization_reversedepth_Omega_rx[i].push_back(vEpolarization_Omega_rx[i][vEpolarization_Omega_rx[i].size()-1-j]);
            vEpolarization_reversedepth_Psi_rx[i].push_back(vEpolarization_Psi_rx[i][vEpolarization_Psi_rx[i].size()-1-j]);

        }

        gpolarization_reversedepth_Omega_rx[i]=new TGraph(vreversedepth[i].size(),&vreversedepth[i][0],&vpolarization_reversedepth_Omega_rx[i][0]);
        gpolarization_reversedepth_Psi_rx[i]=new TGraph(vreversedepth[i].size(),&vreversedepth[i][0],&vpolarization_reversedepth_Psi_rx[i][0]);

        gEpolarization_reversedepth_Omega_rx[i]=new TGraph(vreversedepth[i].size(),&vreversedepth[i][0],&vEpolarization_reversedepth_Omega_rx[i][0]);
        gEpolarization_reversedepth_Psi_rx[i]=new TGraph(vreversedepth[i].size(),&vreversedepth[i][0],&vEpolarization_reversedepth_Psi_rx[i][0]);

        gdiffepsilon_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vdiffepsilon_tx[i][0]);

        gepsilon1_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vepsilon1_rx[i][0]);
        gepsilon2_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vepsilon2_rx[i][0]);
        gdiffepsilon_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vdiffepsilon_rx[i][0]);

        // Along-path diagnostic graphs for the special pulser depth
        g_deltan[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vdeltan[i][0]);
        g_notflipped_alongpath[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vnotflipped_alongpath[i][0]);
        g_theta1_alongpath[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vtheta1_alongpath[i][0]);
        g_theta2_alongpath[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vtheta2_alongpath[i][0]);
        g_thetape1_alongpath[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vthetape1_alongpath[i][0]);
        g_thetape2_alongpath[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vthetape2_alongpath[i][0]);
        g_thetape1_phipe1_alongpath[i]=new TGraph(vphipe1_alongpath[i].size(),&vphipe1_alongpath[i][0],&vthetape1_alongpath[i][0]);
        g_thetape2_phipe2_alongpath[i]=new TGraph(vphipe2_alongpath[i].size(),&vphipe2_alongpath[i][0],&vthetape2_alongpath[i][0]);
        g_phipe1_alongpath[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vphipe1_alongpath[i][0]);
        g_phipe2_alongpath[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vphipe2_alongpath[i][0]);
        g_deltan_pulserdepth[i]=new TGraph(vdepth_step[i].size(),&vdepth_step[i][0],&vdeltan[i][0]);

        g_depth_istep[i]=new TGraph(vdepth_step[i].size(),&vdepth_step[i][0],&vistep[i][0]);

        g_attenlengths[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vattenlengths[i][0]);
        g_receive_launch[i]=new TGraph(vreceiveangle[i].size(),&vreceiveangle[i][0],&vlaunchangle[i][0]);
        g_receive[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vreceiveangle[i][0]);
        g_launch[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vlaunchangle[i][0]);
        g_output6[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&voutput6[i][0]);
        g_output7[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&voutput7[i][0]);
        g_output8[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&voutput8[i][0]);

    }


    string filename;


    TGraph *g1[NSTATIONS];

    // ------------------------------------------------------------------------
    // First main canvas: measured SNR-vs-depth data for each station
    // with room to overlay model curves.
    // ------------------------------------------------------------------------
    TCanvas *c1=new TCanvas("c1","c1",800,800);
    c1->Divide(2,3);
    TH2D *h2[NSTATIONS];

    // Plot limits for SNR/voltage-like quantities by station and polarization
    double vmax[2][NSTATIONS]={{100.,40.,40.,40.,40.,100.},
        {100.,40.,40.,40.,40.,100.}};

    // Plot limits for power-like quantities
    double pmax[2][NSTATIONS]={{1.E5,1.E5,1.E5,1.E5,1.E5,1.E5},
        {1.E5,1.E5,1.E5,1.E5,1.E5,1.E5}};


    for (int i=minstation;i<=maxstation;i++) {

        h2[i]=new TH2D("h2","h2",200,-1500.,-850,100,0.,vmax[WHICHPOL][i]);
        h2[i]->SetXTitle("SpiceCore pulser height (m)");
        h2[i]->SetYTitle("V_{SNR}^{max}");
    }


    for (int i=minstation;i<=maxstation;i++) {

        c1->cd(i+1);

        h2[i]->Draw();

        // Draw measured SNR points where available
        if (i!=5) {
            g1[i]=new TGraphErrors(vdepth_data[i].size(),&vdepth_data[i][0],&vsnrmax[i][0],&vdepth_data_err[i][0],&vsnrmax_err[i][0]);
            g1[i]->SetMarkerColor(icolors[i]);
            g1[i]->SetMarkerSize(1.);
            g1[i]->SetMarkerStyle(20);
            g1[i]->Draw("pesame");
        }

        // Placeholder for fit-parameter uncertainty extraction
        deltan_obs_err[i]=f1[i]->GetParError(1);

        // Style model components
        g_parameter0[i]->SetLineWidth(3);
        g_parameter0[i]->SetLineColor(kGreen);

        g_atten[i]->SetLineWidth(3);
        g_atten[i]->SetLineColor(kOrange);

        g_atten_beam[i]->SetLineWidth(3);
        g_atten_beam[i]->SetLineStyle(kDashed);
        g_atten_beam[i]->SetLineColor(kOrange);

        g_atten_beam_crosspol[i]->SetLineWidth(3);
        g_atten_beam_crosspol[i]->SetLineStyle(kDashed);
        g_atten_beam_crosspol[i]->SetLineColor(kOrange);

        gfunc_noadjust[i]->SetLineWidth(3);
        gfunc_noadjust[i]->SetLineColor(kBlack);

        g_atten_beam_crosspol_nointerferencefunc[i]->SetLineWidth(3);
        g_atten_beam_crosspol_nointerferencefunc[i]->SetLineStyle(kDashed);
        g_atten_beam_crosspol_nointerferencefunc[i]->SetLineColor(kBlack);

        g_atten_beam_crosspol_func[i]->SetLineWidth(3);
        g_atten_beam_crosspol_func[i]->SetLineColor(kBlack);

        // Same interference model, but parameterized in total distance instead of depth
        f1_distances[i]=new TF1("f1",sfunc_distances.c_str(),0.,5000.);
        f1_nointerference_distances[i]=new TF1("f1_nointerference",sfunc_nointerference_distances.c_str(),0.,5000.);

        f1_distances[i]->SetLineColor(kBlack);


        // --------------------------------------------------------------------
        // Distance-domain "big picture" scan for this station.
        // This mirrors the depth-based computation above, but over total
        // source-receiver distance rather than a prescribed pulser depth list.
        // --------------------------------------------------------------------
        #pragma omp parallel for schedule(dynamic)
        for (int j=0;j<NDISTANCES_BIGPIC-(jmin[i]+1);j++) {


            double posstation[2];
            posstation[0]=sqrt(pow(station_coords[i][0]-pulser_coords[0],2)+pow(station_coords[i][1]-pulser_coords[1],2));
            posstation[1]=station_depths[i];

            double pospulser[2];
            pospulser[0]=0.;
            pospulser[1]=vpseudodepths_bigpic[i][j];


            double x0=0;
            double z0=pospulser[1];
            double x1=posstation[0];
            double z1=posstation[1];

            double *getresults=IceRayTracing(x0,z0,x1,z1);

            double lvalue;
            vector<double> res;
            vector<double> zs;
            double launch_angle;
            double receive_angle;
            double *paramsd;
            double *paramsra;

            // Solve direct or refracted branch if available
            if (getresults[6]!=-1000) {
                paramsd=GetDirectRayPar(z0,x1,z1);
                launch_angle=paramsd[1]/DEGRAD;
                receive_angle=paramsd[0]/DEGRAD;
                GetFullDirectRayPath(z0, x1, z1, paramsd[3], res, zs);
            }
            else if (getresults[8]!=1000) {
                paramsd=GetDirectRayPar(z0,x1,z1);
                GetReflectedRayPar(z0, x1 ,z1);
                double LangR=paramsd[1];
                double RangR=paramsd[0];
                paramsra=GetRefractedRayPar(z0, x1 ,z1,LangR,RangR);
                launch_angle=paramsra[1]/DEGRAD;

                receive_angle=paramsra[0]/DEGRAD;
                GetFullRefractedRayPath(z0, x1, z1, paramsra[7], paramsra[3], res, zs);

            }


            double atten=1.;
            double sumphase=0;

            if (getresults[6]!=-1000 || getresults[8]!=-1000) {

                TVector3 yhat(station_coords[i][0]-pulser_coords[0],
                    station_coords[i][1]-pulser_coords[1],
                    0.);

                if (yhat.Mag()<HOWSMALLISTOOSMALL){
                    cout << "yhat mag is " << yhat.Mag() << "\n";
                }

                yhat.SetMag(1.);

                for (int istep=UZAIRSTEP;istep<res.size();istep+=UZAIRSTEP) {

                    vector<double> nvec_thisstep;
                    TVector3 rhat_thisstep;

                    nvec_thisstep.resize(3);

                    nvec_thisstep[0]=gn1->Eval(zs[istep]);
                    nvec_thisstep[1]=gn2->Eval(zs[istep]);
                    nvec_thisstep[2]=gn3->Eval(zs[istep]);

                    if (istep>0) {
                        rhat_thisstep[0]=-1.*(res[istep]-res[istep-UZAIRSTEP])*yhat[0];
                        rhat_thisstep[1]=-1.*(res[istep]-res[istep-UZAIRSTEP])*yhat[1];
                        rhat_thisstep[2]=-1.*(zs[istep]-zs[istep-UZAIRSTEP]);

                        double length=rhat_thisstep.Mag();

                        if (rhat_thisstep.Mag()<HOWSMALLISTOOSMALL)
                        cout << "rhat_thissetp mag is " << rhat_thisstep.Mag() << "\n";

                        rhat_thisstep.SetMag(1.);

                        double atten_length=GetIceAttenuationLength(zs[istep], freq/1.E9);

                        atten*=exp(-1.*length/atten_length);

                        // This computes local delta n but is not yet accumulated into sumphase
                        double deltan_alongpath=getDeltaN(BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1,p_e2);

                    }
                    else {
                        rhat_thisstep=rhat[i];

                    }

                }

                // Store time/phase proxy for this distance-domain point
                vtimediff_bigpic[i].push_back(sumphase);


                if (DORAYTRACING) {

                    TVector3 plusz(0.,0.,1.);

                    TVector3 vrotate=plusz.Cross(yhat);
                    if (vrotate.Mag()<HOWSMALLISTOOSMALL)
                    cout << "vrotate mag is " << vrotate.Mag() << "\n";

                    vrotate.SetMag(1.);
                    rhat[i]=plusz;
                    rhat[i].Rotate(launch_angle,vrotate);


                    rhat_receive[i]=plusz;
                    rhat_receive[i].Rotate(receive_angle,vrotate);
                }
                else {

                    rhat[i].SetX(station_coords[i][0]-pulser_coords[0]);
                    rhat[i].SetY(station_coords[i][1]-pulser_coords[1]);
                    rhat[i].SetZ(station_depths[i]-vpseudodepths_bigpic[i][j]);

                    rhat_receive[i].SetX(rhat[i][0]);
                    rhat_receive[i].SetY(rhat[i][1]);
                    rhat_receive[i].SetZ(rhat[i][2]);
                }

            }
            else {

                // Straight-line fallback if no ray solution exists
                rhat[i].SetX(station_coords[i][0]-pulser_coords[0]);
                rhat[i].SetY(station_coords[i][1]-pulser_coords[1]);
                rhat[i].SetZ(station_depths[i]-vpseudodepths_bigpic[i][j]);

                rhat_receive[i].SetX(rhat[i][0]);
                rhat_receive[i].SetY(rhat[i][1]);
                rhat_receive[i].SetZ(rhat[i][2]);
            }
            if (rhat[i].Mag()<HOWSMALLISTOOSMALL){
                cout << "rhat[i] mag is " << rhat[i].Mag() << "\n";
            }
            
            rhat[i].SetMag(1.);
            
            if (rhat_receive[i].Mag()<HOWSMALLISTOOSMALL){
                cout << "rhat_receive[i] mag is " << rhat_receive[i].Mag() << "\n";
            }

            rhat_receive[i].SetMag(1.);

        }

        // Distance-based measured-data graph
        if (i!=5){
            g1_distances[i]=new TGraphErrors(vtotal_distances[i].size(),&vtotal_distances[i][0],&vsnrmax[i][0],&vtotal_distances_err[i][0],&vsnrmax_err[i][0]);
        }

        // Distance-domain model curves
        g_atten_beam_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_beam_bigpic[i][0]);
        g_atten_beam_crosspol_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_beam_crosspol_bigpic[i][0]);
        g_atten_beam_crosspol_nointerferencefunc_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_beam_crosspol_nointerferencefunc_bigpic[i][0]);
        g_atten_beam_crosspol_func_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_beam_crosspol_func_bigpic[i][0]);

        g_sumphase_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vtimediff_bigpic[i][0]);
    }

    // ------------------------------------------------------------------------
    // Output directory naming convention:
    // includes whether the indicatrix is constant, whether the ice is isotropic/
    // uniaxial/biaxial, the frequency, and the TX/RX cross-pol angles.
    // ------------------------------------------------------------------------
    const int NAXIAL=3;
    const int NCONSTANT=2;
    string sbiaxial[NAXIAL]={"isotropic_","uniaxial_",""};
    string sconstant[NCONSTANT]={"","constant_"};

    string sdir="cpol_main_plots_" + sconstant[CONSTANTINDICATRIX] + sbiaxial[BIAXIAL+1] + to_string((int)(freq/1.E6)) + "MHz_angletx" + to_string(CROSSPOLANGLE_TX_INT) + "_anglerx" + to_string(CROSSPOLANGLE_RX_INT) + "/";
    gSystem->mkdir(sdir.c_str(), /*recursive=*/true);

    // Save the first data canvas
    string sname=sdir+"c1.pdf";
    c1->Print(sname.c_str());
    cout << "just printed c1.\n";

    // Overlay of one particular model component across stations
    TCanvas *c1b=new TCanvas("c1b","c1b",800,800);
    h2[0]->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {
        g_atten_beam_crosspol_func[istations]->SetLineColor(icolors_dave[istations]);
        g_atten_beam_crosspol_func[istations]->Draw("same");
    }

    sname=sdir+"c1b.pdf";
    c1b->Print(sname.c_str());
    cout << "printed c1b.\n";

    // ------------------------------------------------------------------------
    // Distance-domain comparison plots: measured SNR vs total path length
    // together with selected model curves.
    // ------------------------------------------------------------------------
    TCanvas *c3=new TCanvas("c3","c3",800,800);
    c3->Divide(2,3);
    TH2D *h3[6];


    for (int i=minstation;i<=maxstation;i++) {

        c3->cd(i+1);
        h3[i]=new TH2D("h2","h2",3500,1000.,4500,1000,0.,vmax[WHICHPOL][i]);

        h3[i]->Draw();

        if (i!=5) {
            g1_distances[i]->SetMarkerColor(icolors[i]);
            g1_distances[i]->SetMarkerSize(0.5);
            g1_distances[i]->SetMarkerStyle(20);
            g1_distances[i]->Draw("pesame");
        }

        g_atten_beam_distances[i]->SetLineColor(kOrange);
        g_atten_beam_distances[i]->SetLineStyle(kDashed);
        g_atten_beam_distances[i]->SetLineWidth(2);
        g_atten_beam_distances[i]->Draw("lsame");

        g_atten_beam_crosspol_nointerferencefunc_distances[i]->SetLineColor(kBlack);
        g_atten_beam_crosspol_nointerferencefunc_distances[i]->SetLineWidth(2);
        g_atten_beam_crosspol_nointerferencefunc_distances[i]->SetLineStyle(kDashed);
        g_atten_beam_crosspol_nointerferencefunc_distances[i]->Draw("lsame");

        g_atten_beam_crosspol_func_distances[i]->SetLineColor(kBlack);
        g_atten_beam_crosspol_func_distances[i]->SetLineWidth(2);
        g_atten_beam_crosspol_func_distances[i]->SetLineStyle(kSolid);
        g_atten_beam_crosspol_func_distances[i]->Draw("lsame");

        cout << "f1_distances evaluated at 0 is " << f1_distances[i]->Eval(0.) << "\n";
    }

    sname=sdir+"c3.pdf";
    c3->Print(sname.c_str());


    // ------------------------------------------------------------------------
    // Compare observed and expected birefringence delta-n as a function of
    // station angle alpha relative to the ice-flow direction.
    // ------------------------------------------------------------------------
    TGraphErrors *g_obs=new TGraphErrors(6,alpha_deg,deltan_obs,zeroes,deltan_obs_err);
    TGraph *g_exp=new TGraph(6,alpha_deg,deltan_exp);

    // Analytic expectation for delta n vs alpha in the toy model
    TF1 *fgetDeltaN=new TF1("f1","2*([0]-([0]-[1]/2.)/sqrt(1-[1]/[0]*cos(x*3.14159/180.)*cos(x*3.14159/180.)*(1-[1]/(4.*[0]))))",0.,90.);
    fgetDeltaN->SetParameter(0,NICE);
    fgetDeltaN->SetParameter(1,DELTAN);

    cout << "function evaluated is " << fgetDeltaN->Eval(10.) << "\t" << fgetDeltaN->Eval(90.) << "\n";

    TH2D *h2_compare=new TH2D("h2_compare","h2_compare",100,0.,180.,100,1.E-4,1.0);

    h2_compare->SetXTitle("#alpha (degrees)");
    h2_compare->SetYTitle("#Delta n");
    TCanvas *c2=new TCanvas("c2","c2",800,800);
    c2->SetLogy();
    gStyle->SetOptStat(0);
    h2_compare->Draw();
    fgetDeltaN->SetLineColor(kBlack);

    g_exp->SetMarkerStyle(21);
    g_exp->SetMarkerColor(kBlack);
    g_exp->SetMarkerSize(1.0);
    g_exp->Draw("psame");
    g_obs->SetMarkerStyle(21);
    g_obs->SetMarkerColor(kRed);
    g_obs->SetMarkerSize(1.5);
    g_obs->SetLineColor(kBlack);
    g_obs->SetLineWidth(3);
    g_obs->Draw("pesame");
    sname=sdir+"compare.pdf";
    c2->Print(sname.c_str());


    // ------------------------------------------------------------------------
    // Plot the angle between each Poynting vector S1/S2 and the propagation
    // direction k as a function of pulser depth.
    // ------------------------------------------------------------------------
    TGraph *g_angleS1_k[6];
    TGraph *g_angleS2_k[6];

    TCanvas *c7=new TCanvas("c7","c7",800,800);

    TH2D *h7=new TH2D("","",100,-1500.,-500.,100,0.00011,0.2);
    h7->SetXTitle("Pulser height (m)");
    h7->SetYTitle("Angle between S_{1,2} and k (^{ o })");
    h7->GetYaxis()->SetTitleOffset(1.55);
    h7->GetXaxis()->SetTitleOffset(1.2);
    h7->GetXaxis()->SetNdivisions(504);

    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLogy();
    h7->Draw();

    auto legend7 = new TLegend(0.67,0.17,0.93,0.4);
    legend7->SetBorderSize();
    legend7->SetTextSize(0.04);

    auto legend7b = new TLegend(0.27,0.17,0.43,0.3);
    legend7b->SetBorderSize();
    legend7b->SetTextSize(0.04);

    for (int i=0;i<NSTATIONS;i++) {

        g_angleS1_k[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vangle_Shat_e1_khat[i][0]);

        g_angleS1_k[i]->SetMarkerColor(icolors[i]);
        g_angleS1_k[i]->SetLineColor(icolors[i]);
        g_angleS1_k[i]->SetLineStyle(kSolid);
        g_angleS1_k[i]->SetLineWidth(3);

        g_angleS1_k[i]->Draw("lsame");

        g_angleS2_k[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vangle_Shat_e2_khat[i][0]);

        g_angleS2_k[i]->SetMarkerColor(icolors[i]);
        g_angleS2_k[i]->SetLineColor(icolors[i]);
        g_angleS2_k[i]->SetLineStyle(kDashed);
        g_angleS2_k[i]->SetLineWidth(3);

        g_angleS2_k[i]->Draw("lsame");

        legend7->AddEntry(g_angleS1_k[i],snames[i].c_str(),"l");
    }

    legend7b->AddEntry(g_angleS1_k[0],"S_{1}","l");
    legend7b->AddEntry(g_angleS2_k[0],"S_{2}","l");

    legend7->Draw("same");
    legend7b->Draw("same");

    sname=sdir+"angle_Shat_khat.pdf";
    c7->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Time-delay / phase-difference vs pulser depth
    // ------------------------------------------------------------------------
    TCanvas *c9=new TCanvas("c9","c9",800,800);
    TH2D *h9=new TH2D("","",100,-1600.,-600.,100,-20.,90.);
    auto legend26a = new TLegend(0.55,0.55,0.88,0.88);
    legend26a->SetBorderSize(0);
    legend26a->SetTextSize(0.05);

    h9->GetXaxis()->SetTitle("Pulser height (m)");
    h9->GetYaxis()->SetTitle("Time difference (ns)");
    h9->GetXaxis()->SetTitleOffset(1.1);
    h9->GetYaxis()->SetTitleOffset(1.2);
    h9->GetXaxis()->SetTitleSize(0.04);
    h9->GetYaxis()->SetTitleSize(0.04);

    h9->GetXaxis()->SetNdivisions(504);
    h9->Draw();


    for (int istations=minstation;istations<=maxstation;istations++) {

        g_sumphase[istations]->SetMarkerColor(icolors[istations]);
        g_sumphase[istations]->SetLineColor(icolors[istations]);
        g_sumphase[istations]->SetLineStyle(kSolid);
        g_sumphase[istations]->SetLineWidth(3);


        g_sumphase[istations]->Draw("lsame");

        legend26a->AddEntry(g_sumphase[istations],snames[istations].c_str(),"l");


    }
    legend26a->Draw("same");
    sname=sdir+"sumphase.pdf";
    c9->Print(sname.c_str());

    // Plot of whether the eigenvector tracking indicates a sign flip
    TCanvas *c9a=new TCanvas("c9a","c9a",800,800);
    c9a->Divide(2,3);
    TH2D *h9a=new TH2D("h9a","h9a",100,-1800.,0.,100,-1.1,1.1);

    h9a->GetXaxis()->SetTitle("Pulser height (m)");
    h9a->GetYaxis()->SetTitle("Not flipped");
    h9a->GetXaxis()->SetTitleOffset(1.1);
    h9a->GetYaxis()->SetTitleOffset(1.45);


    for (int istations=minstation;istations<=maxstation;istations++) {
        c9a->cd(istations+1);
        h9a->Draw();
        g_notflipped[istations]->SetMarkerColor(icolors[istations]);
        g_notflipped[istations]->SetLineColor(icolors[istations]);
        g_notflipped[istations]->SetLineStyle(kSolid);
        g_notflipped[istations]->SetLineWidth(2);

        g_notflipped[istations]->Draw("lsame");
    }
    sname=sdir+"notflipped.pdf";
    c9a->Print(sname.c_str());

    // Total geometric path length for one selected station (currently only minstation)
    TCanvas *c9b=new TCanvas("c9b","c9b",800,800);
    TH2D *h9b=new TH2D("h9b","h9b",100,0.,5000.,100,0.,5000.);
    h9b->Draw();

    for (int istations=minstation;istations<=minstation;istations++) {

        g_sumlength[istations]->SetMarkerColor(icolors[istations]);
        g_sumlength[istations]->SetLineColor(icolors[istations]);
        g_sumlength[istations]->SetLineStyle(kSolid);
        g_sumlength[istations]->SetLineWidth(2);

        g_sumlength[istations]->Draw("lsame");

    }
    sname=sdir+"sumlength.pdf";
    c9b->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Frequency spectra at one selected pulser depth (depth_special)
    // ------------------------------------------------------------------------
    TCanvas *c10=new TCanvas("c10","c10",800,800);
    c10->Divide(2,3);
    TH2D *h10=new TH2D("h10","h10",100,freqmin,freqmax,100,1.,pmax[0][0]);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c10->cd(istations+1);

        g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->SetMarkerColor(icolors[istations]);
        g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->SetLineColor(icolors[istations]);
        g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->SetLineStyle(kSolid);
        g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->SetLineWidth(2);

        g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->Draw("al");
    }
    sname=sdir+"spectra.pdf";
    c10->Print(sname.c_str());


    // ------------------------------------------------------------------------
    // Along-path delta n for the special pulser depth
    // ------------------------------------------------------------------------
    TCanvas *c11=new TCanvas("c11","c11",800,800);

    TH2D *h11=new TH2D("h11","h11",100,-1200.,0.,100,1.E-5,1.E-2);
    string stitle ="Depth along path for pulser depth of " + to_string(depth_special) + "\n";
    h11->SetXTitle("Depth along path (m)");
    h11->SetYTitle("#Delta n");

    h11->Draw();
    gPad->SetLogy();

    for (int istations=minstation;istations<=maxstation;istations++) {

        g_deltan_pulserdepth[istations]->SetMarkerColor(icolors[istations]);
        g_deltan_pulserdepth[istations]->SetLineColor(icolors[istations]);
        g_deltan_pulserdepth[istations]->SetLineStyle(kSolid);
        g_deltan_pulserdepth[istations]->SetLineWidth(2);

        g_deltan_pulserdepth[istations]->Draw("lsame");
    }
    sname=sdir+"deltan.pdf";
    c11->Print(sname.c_str());

    // Attenuation length along the special path
    TCanvas *c11b=new TCanvas("c11b","c11b",800,800);

    TH2D *h11b=new TH2D("h11b","h11b",100,0.,4000.,100,1.,1.E5);

    gPad->SetLogy();

    for (int istations=minstation;istations<=maxstation;istations++) {

        g_attenlengths[istations]->SetMarkerColor(icolors[istations]);
        g_attenlengths[istations]->SetLineColor(icolors[istations]);
        g_attenlengths[istations]->SetLineStyle(kSolid);
        g_attenlengths[istations]->SetLineWidth(2);

        g_attenlengths[istations]->Draw("al");
    }
    sname=sdir+"attenlengths.pdf";
    c11b->Print(sname.c_str());

    // Flipping diagnostic along the path
    TCanvas *c11c=new TCanvas("c11c","c11c",800,800);

    TH2D *h11c=new TH2D("h11c","h11c",100,0.,2000.,100,-1.1,1.1);
    h11c->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {

        g_notflipped_alongpath[istations]->SetLineColor(icolors[istations]);
        g_notflipped_alongpath[istations]->SetLineStyle(kSolid);
        g_notflipped_alongpath[istations]->SetLineWidth(2);
        g_notflipped_alongpath[istations]->Draw("lsame");

    }
    sname=sdir+"notflipped_alongpath.pdf";
    c11c->Print(sname.c_str());

    // Theta angles of the two eigenmodes along the path
    TCanvas *c11d=new TCanvas("c11d","c11d",800,800);

    TH2D *h11d=new TH2D("h11d","h11d",100,0.,2000.,100,-180.,180.);
    h11d->Draw();


    for (int istations=minstation;istations<=maxstation;istations++) {

        g_theta1_alongpath[istations]->SetLineColor(icolors[istations]);
        g_theta1_alongpath[istations]->SetLineStyle(kSolid);
        g_theta1_alongpath[istations]->SetLineWidth(2);
        g_theta1_alongpath[istations]->Draw("lsame");
        g_theta2_alongpath[istations]->SetLineColor(icolors[istations]);
        g_theta2_alongpath[istations]->SetLineStyle(kDashed);
        g_theta2_alongpath[istations]->SetLineWidth(2);
        g_theta2_alongpath[istations]->Draw("lsame");

    }
    sname=sdir+"thetas_alongpath.pdf";
    c11d->Print(sname.c_str());

    // Theta of p_e1 and p_e2 along the path
    TCanvas *c11e=new TCanvas("c11e","c11e",800,800);

    TH2D *h11e=new TH2D("h11e","h11e",100,0.,2000.,100,0.,180.);
    h11e->Draw();
    h11e->SetXTitle("Phi (degrees)");
    h11e->SetYTitle("Theta (degrees)");

    for (int istations=minstation;istations<=maxstation;istations++) {

        g_thetape1_alongpath[istations]->SetLineColor(icolors[istations]);
        g_thetape1_alongpath[istations]->SetLineStyle(kSolid);
        g_thetape1_alongpath[istations]->SetLineWidth(2);
        g_thetape1_alongpath[istations]->Draw("lsame");
        g_thetape2_alongpath[istations]->SetLineColor(icolors[istations]);
        g_thetape2_alongpath[istations]->SetLineStyle(kDashed);
        g_thetape2_alongpath[istations]->SetLineWidth(2);
        g_thetape2_alongpath[istations]->Draw("lsame");

    }
    sname=sdir+"thetas_pe_alongpath.pdf";
    c11e->Print(sname.c_str());

    // Phi-vs-theta trajectory of the two eigenvectors along the path
    TCanvas *c11f=new TCanvas("c11f","c11f",800,800);

    TH2D *h11f=new TH2D("h11f","h11f",100,-180.,180.,100,0.,180.);
    h11f->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {

        g_thetape1_phipe1_alongpath[istations]->SetLineColor(icolors[istations]);
        g_thetape1_phipe1_alongpath[istations]->SetLineStyle(kSolid);
        g_thetape1_phipe1_alongpath[istations]->SetLineWidth(2);
        g_thetape1_phipe1_alongpath[istations]->Draw("lsame");
        g_thetape2_phipe2_alongpath[istations]->SetLineColor(icolors[istations]);
        g_thetape2_phipe2_alongpath[istations]->SetLineStyle(kDashed);
        g_thetape2_phipe2_alongpath[istations]->SetLineWidth(2);
        g_thetape2_phipe2_alongpath[istations]->Draw("lsame");

    }
    sname=sdir+"thetas_phis_pe_alongpath.pdf";
    c11f->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot the depth->index mapping used throughout the code.
    // This is useful because many later selections use g_idepth->Eval(depth)
    // to convert a physical pulser depth into the corresponding vector index.
    // ------------------------------------------------------------------------
    TCanvas *c12=new TCanvas("c12","c12",800,800);

    TH2D *h12=new TH2D("h12","h12",100,0.,2000.,100,0.,200.);

    gPad->SetLogy();

    for (int istations=minstation;istations<=maxstation;istations++) {

        g_idepth[istations]->SetMarkerColor(icolors[istations]);
        g_idepth[istations]->SetLineColor(icolors[istations]);
        g_idepth[istations]->SetLineStyle(kSolid);
        g_idepth[istations]->SetLineWidth(2);

        g_idepth[istations]->Draw("al");
    }
    sname=sdir+"index.pdf";
    c12->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot launch and receive angles:
    //   panel 1: receive angle vs pulser depth
    //   panel 2: launch angle vs pulser depth
    //   panel 3: launch angle vs receive angle
    // ------------------------------------------------------------------------
    TCanvas *c13=new TCanvas("c13","c13",800,800);

    TH2D *h13a=new TH2D("h13a","h13a",100,-2000.,0.,100,0.,180.);
    TH2D *h13b=new TH2D("h13b","h13b",100,-2000.,0.,100,0.,180.);
    TH2D *h13c=new TH2D("h13c","h13c",100,0.,180.,100,0.,180.);

    c13->Divide(2,2);
    gPad->SetLogy();


    for (int istations=minstation;istations<=maxstation;istations++) {

        g_receive_launch[istations]->SetMarkerColor(icolors[istations]);
        g_receive_launch[istations]->SetLineColor(icolors[istations]);
        g_receive_launch[istations]->SetLineStyle(kSolid);
        g_receive_launch[istations]->SetLineWidth(2);

        g_launch[istations]->SetMarkerColor(icolors[istations]);
        g_launch[istations]->SetLineColor(icolors[istations]);
        g_launch[istations]->SetLineStyle(kSolid);
        g_launch[istations]->SetLineWidth(2);

        g_receive[istations]->SetMarkerColor(icolors[istations]);
        g_receive[istations]->SetLineColor(icolors[istations]);
        g_receive[istations]->SetLineStyle(kSolid);
        g_receive[istations]->SetLineWidth(2);
    }

    c13->cd(1);
    h13a->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {
        h13a->SetXTitle("Receive angles (degrees)");
        g_receive[istations]->Draw("lsame");
    }

    c13->cd(2);
    h13b->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {
        h13b->SetXTitle("Launch angles (degrees)");
        g_launch[istations]->Draw("lsame");
    }

    c13->cd(3);
    h13c->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {
        h13c->SetXTitle("Receive angle (degrees)");
        h13c->SetYTitle("Launch angle (degrees)");

        g_receive_launch[istations]->Draw("lsame");
    }

    sname=sdir+"receive_launch.pdf";
    c13->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot the derived quantity V(z) computed from n1(z), n2(z), n3(z).
    // This comes from getV(...) and summarizes some property of the indicatrix
    // or birefringence model as a function of depth.
    // ------------------------------------------------------------------------
    TCanvas *c14=new TCanvas("c14","c14",800,800);

    TH2D *h14=new TH2D("h14","h14",100,-1800.,0.,100,0.,90.);
    h14->Draw();
    h14->SetTitle("");
    h14->GetXaxis()->SetTitleOffset(1.2);
    h14->SetXTitle("Height (meters)");
    h14->SetYTitle("V_{z} (degrees)");

    g_V->SetLineColor(kBlack);
    g_V->SetLineStyle(kSolid);
    g_V->SetLineWidth(2);

    g_V->Draw("lsame");

    sname=sdir+"V.pdf";
    c14->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot the three principal refractive indices vs depth.
    // gn1, gn2, gn3 were built from the smoothed depth-dependent profiles.
    // ------------------------------------------------------------------------
    TCanvas *c15=new TCanvas("c15","c15",800,800);
    c15->SetLeftMargin(0.15);
    c15->SetBottomMargin(0.15);

    gn1->SetLineColor(kGray+1);
    gn1->SetLineStyle(kSolid);
    gn1->SetLineWidth(3);

    gn2->SetLineColor(kGray+2);
    gn2->SetLineStyle(kSolid);
    gn2->SetLineWidth(3);

    gn3->SetLineColor(kGray+3);
    gn3->SetLineStyle(kSolid);
    gn3->SetLineWidth(3);



    // ------------------------------------------------------------------------
    // Plot raw ray-tracing diagnostic outputs stored in voutput6/7/8.
    // These come directly from IceRayTracing(...) and help diagnose which
    // solution branches are present or how the solver is behaving.
    // ------------------------------------------------------------------------
    TH2D *h16=new TH2D("h16","h16",100,-1800.,0.,100,1.77,1.79);

    h16->SetTitle("");
    h16->SetXTitle("Height (meters)");
    h16->GetYaxis()->SetTitleOffset(1.4);
    h16->GetYaxis()->SetMaxDigits(2);
    h16->GetXaxis()->SetTitleOffset(1.0);
    h16->GetXaxis()->SetTitleSize(0.05);
    h16->GetYaxis()->SetTitleSize(0.05);
    h16->GetXaxis()->SetNdivisions(504);
    h16->SetYTitle("Principal axis");
    h16->Draw();

    gn1->Draw("lsame");
    gn2->Draw("lsame");
    gn3->Draw("lsame");

    auto legend2 = new TLegend(0.6,0.2,0.8,0.4);

    legend2->SetTextSize(0.05);
    legend2->SetBorderSize(0);
    gn1->SetName("gn1");
    gn2->SetName("gn2");
    gn3->SetName("gn3");
    legend2->AddEntry("gn1","n_{ #alpha}","l");
    legend2->AddEntry("gn2","n_{ #beta}","l");
    legend2->AddEntry("gn3","n_{ #gamma}","l");
    legend2->Draw("same");

    sname=sdir+"n123_zoomed.pdf";
    c15->Print(sname.c_str());


    TCanvas *c17=new TCanvas("c17","c17",800,800);
    TH2D *h17=new TH2D("h17","h17",100,-1800.,0.,100,-1100.,200.);
    c17->Divide(1,3);
    c17->cd(1);
    h17->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {
        g_output6[istations]->SetLineColor(icolors[istations]);
        g_output6[istations]->SetLineWidth(2);
        g_output6[istations]->Draw("lsame");
    }

    c17->cd(2);
    h17->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {
        g_output7[istations]->SetLineColor(icolors[istations]);
        g_output7[istations]->SetLineWidth(2);
        g_output7[istations]->Draw("lsame");
    }

    c17->cd(3);
    h17->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {
        g_output8[istations]->SetLineColor(icolors[istations]);
        g_output8[istations]->SetLineWidth(2);
        g_output8[istations]->Draw("lsame");
    }

    sname=sdir+"outputs.pdf";
    c17->Print(sname.c_str());


    // ------------------------------------------------------------------------
    // Plot TX-side mode angles in the "clock" convention:
    //   top row    : theta1 and thetaE1
    //   bottom row : theta2 and thetaE2
    //
    // Solid and dashed lines compare D-like and E-like angle conventions.
    // ------------------------------------------------------------------------
    TCanvas *c18=new TCanvas("c18","c18",800,800);
    TH2D *h18a=new TH2D("h18","h18",100,-1800.,0.,100,0.,180.);
    TH2D *h18b=new TH2D("h18","h18",100,-1800.,0.,100,-90.,90.);
    c18->Divide(3,4);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c18->cd(istations+1);
        h18a->Draw();
        gtxdepth_theta1[istations]->SetLineColor(icolors[istations]);
        gtxdepth_theta1[istations]->SetLineWidth(2);
        gtxdepth_theta1[istations]->Draw("lsame");
        gtxdepthE_theta1[istations]->SetLineColor(icolors[istations]);
        gtxdepthE_theta1[istations]->SetLineWidth(2);
        gtxdepthE_theta1[istations]->SetLineStyle(kDashed);
        gtxdepthE_theta1[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c18->cd(NSTATIONS+istations+1);
        h18b->Draw();
        gtxdepth_theta2[istations]->SetLineColor(icolors[istations]);
        gtxdepth_theta2[istations]->SetLineWidth(2);
        gtxdepth_theta2[istations]->Draw("lsame");
        gtxdepthE_theta2[istations]->SetLineColor(icolors[istations]);
        gtxdepthE_theta2[istations]->SetLineWidth(2);
        gtxdepthE_theta2[istations]->SetLineStyle(kDashed);
        gtxdepthE_theta2[istations]->Draw("lsame");
    }
    sname=sdir+"anglesontheclock.pdf";
    c18->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot the TX-side dispersion angles:
    // angle between E and D for each eigenmode as a function of pulser depth.
    // ------------------------------------------------------------------------
    TCanvas *c19=new TCanvas("c19","c19",800,800);
    TH2D *h19=new TH2D("h19","h19",100,-1800.,0.,100,-1.,1.);
    c19->Divide(3,4);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c19->cd(istations+1);
        h19->Draw();
        gtxdepth_dispersion1[istations]->SetLineColor(icolors[istations]);
        gtxdepth_dispersion1[istations]->SetLineWidth(2);
        gtxdepth_dispersion1[istations]->Draw("lsame");
    }
    for (int istations=minstation;istations<=maxstation;istations++) {
        c19->cd(NSTATIONS+istations+1);
        h19->Draw();
        gtxdepth_dispersion2[istations]->SetLineColor(icolors[istations]);
        gtxdepth_dispersion2[istations]->SetLineWidth(2);
        gtxdepth_dispersion2[istations]->Draw("lsame");
    }
    sname=sdir+"dispersionangles.pdf";
    c19->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Same as c18, but using the S-clock angular convention.
    // ------------------------------------------------------------------------
    TCanvas *c20=new TCanvas("c20","c20",800,800);
    TH2D *h20a=new TH2D("h20","h20",100,-1800.,0.,100,0.,180.);
    TH2D *h20b=new TH2D("h20","h20",100,-1800.,0.,100,-90.,90.);
    c20->Divide(3,4);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c20->cd(istations+1);
        h20a->Draw();
        gtxdepth_theta1_Sclock[istations]->SetLineColor(icolors[istations]);
        gtxdepth_theta1_Sclock[istations]->SetLineWidth(2);
        gtxdepth_theta1_Sclock[istations]->Draw("lsame");
        gtxdepthE_theta1_Sclock[istations]->SetLineColor(icolors[istations]);
        gtxdepthE_theta1_Sclock[istations]->SetLineWidth(2);
        gtxdepthE_theta1_Sclock[istations]->SetLineStyle(kDashed);
        gtxdepthE_theta1_Sclock[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c20->cd(NSTATIONS+istations+1);

        gtxdepth_theta2_Sclock[istations]->SetLineColor(icolors[istations]);
        gtxdepth_theta2_Sclock[istations]->SetLineWidth(2);

        gtxdepthE_theta2_Sclock[istations]->SetLineColor(icolors[istations]);
        gtxdepthE_theta2_Sclock[istations]->SetLineWidth(2);
        gtxdepthE_theta2_Sclock[istations]->SetLineStyle(kDashed);

        h20b->Draw();
        gtxdepth_theta2_Sclock[istations]->Draw("lsame");
        gtxdepthE_theta2_Sclock[istations]->Draw("lsame");
    }
    sname=sdir+"anglesontheSclock.pdf";
    c20->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot the individual mode-1 voltage contributions into receiver channels
    // r1 and r2.
    // ------------------------------------------------------------------------
    TCanvas *c21=new TCanvas("c21","c21",800,800);
    TH2D *h21a=new TH2D("h21","h21",100,-1800.,0.,100,-1.1,1.1);
    TH2D *h21b=new TH2D("h21","h21",100,-1800.,0.,100,-1.1,1.1);
    c21->Divide(3,4);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c21->cd(istations+1);
        h21a->Draw();
        gV1_r1[istations]->SetLineColor(icolors[istations]);
        gV1_r1[istations]->SetLineWidth(2);
        gV1_r1[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c21->cd(NSTATIONS+istations+1);
        h21b->Draw();
        gV1_r2[istations]->SetLineColor(icolors[istations]);
        gV1_r2[istations]->SetLineWidth(2);
        gV1_r2[istations]->Draw("lsame");
    }

    sname=sdir+"V1.pdf";
    c21->Print(sname.c_str());


    // ------------------------------------------------------------------------
    // Same as c21, but for the mode-2 voltage contributions.
    // ------------------------------------------------------------------------
    TCanvas *c22=new TCanvas("c22","c22",800,800);
    TH2D *h22a=new TH2D("h22","h22",100,-1800.,0.,100,0.,100.);
    TH2D *h22b=new TH2D("h22","h22",100,-1800.,0.,100,-1.1,1.1);
    c22->Divide(3,4);
    for (int istations=minstation;istations<=maxstation;istations++) {
        c22->cd(istations+1);
        h22a->Draw();
        gV2_r1[istations]->SetLineColor(icolors[istations]);
        gV2_r1[istations]->SetLineWidth(2);
        gV2_r1[istations]->Draw("lsame");

    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c22->cd(NSTATIONS+istations+1);

        h22b->Draw();
        gV2_r2[istations]->SetLineColor(icolors[istations]);
        gV2_r2[istations]->SetLineWidth(2);
        gV2_r2[istations]->Draw("lsame");
    }

    sname=sdir+"V2.pdf";
    c22->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot constructive and destructive envelopes for the two receiver channels.
    // ------------------------------------------------------------------------
    TCanvas *c23=new TCanvas("c23","c23",800,800);
    TH2D *h23a=new TH2D("h23","h23",100,-1800.,0.,100,0.,1.1);
    TH2D *h23b=new TH2D("h23","h23",100,-1800.,0.,100,0.,1.1);
    c23->Divide(3,4);
    for (int istations=minstation;istations<=maxstation;istations++) {
        c23->cd(istations+1);
        h23a->Draw();

        genvelope_plus_r1[istations]->SetLineColor(icolors[istations]);
        genvelope_plus_r1[istations]->SetLineWidth(2);
        genvelope_plus_r1[istations]->Draw("lsame");
        genvelope_minus_r1[istations]->SetLineColor(icolors[istations]);
        genvelope_minus_r1[istations]->SetLineWidth(2);
        genvelope_minus_r1[istations]->Draw("lsame");

    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c23->cd(NSTATIONS+istations+1);

        h23b->Draw();
        genvelope_plus_r2[istations]->SetLineColor(icolors[istations]);
        genvelope_plus_r2[istations]->SetLineWidth(2);
        genvelope_plus_r2[istations]->Draw("lsame");
        genvelope_minus_r2[istations]->SetLineColor(icolors[istations]);
        genvelope_minus_r2[istations]->SetLineWidth(2);
        genvelope_minus_r2[istations]->Draw("lsame");
    }
    sname=sdir+"envelopes.pdf";
    c23->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot TX-side epsilon angles for all stations, side-by-side:
    //   left  panel: epsilon_1^T
    //   right panel: epsilon_1^R
    //
    // These epsilon quantities are derived from the S-clock angles of the
    // eigenmodes at the transmitter and receiver.
    // ------------------------------------------------------------------------
    TCanvas *c24=new TCanvas("c24","c24",1600,800);
    c24->Divide(2,1);

    c24->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    TH2D *h24a=new TH2D("","",100,-1600.,-600.,100,-10.,70.);
    titles(h24a, "", "Pulser height (m)", "#epsilon_{ 1}^{ T}");

    h24a->GetXaxis()->SetNdivisions(504);
    h24a->GetYaxis()->SetNdivisions(504);

    h24a->Draw();
    auto legend24 = new TLegend(0.18,0.2,0.36,0.54);
    legend24->SetBorderSize(0);
    legend24->SetTextSize(0.04);

    for (int istations=minstation;istations<=maxstation;istations++) {

        gepsilon1_tx[istations]->SetLineColor(icolors[istations]);
        gepsilon1_tx[istations]->SetLineWidth(2);
        gepsilon1_tx[istations]->Draw("lsame");
        gepsilon2_tx[istations]->SetLineColor(icolors[istations]);
        gepsilon2_tx[istations]->SetLineWidth(2);
        gepsilon2_tx[istations]->SetLineStyle(kDashed);

        gepsilon1_tx[istations]->SetName(snames[istations].c_str());
        legend24->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");
    }

    c24->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    TH2D *h26a=new TH2D("","",100,-1600.,-600.,100,-10.,70.);

    auto legend26 = new TLegend(0.55,0.5,0.88,0.8);

    titles(h26a, "", "Pulser height (m)", "#epsilon_{ 1}^{ R}");

    h26a->GetXaxis()->SetNdivisions(504);
    h26a->GetYaxis()->SetNdivisions(504);

    h26a->Draw();
    legend26->SetBorderSize(0);
    legend26->SetTextSize(0.05);
    for (int istations=minstation;istations<=maxstation;istations++) {

        gepsilon1_rx[istations]->SetLineColor(icolors[istations]);
        gepsilon1_rx[istations]->SetLineWidth(2);
        gepsilon1_rx[istations]->Draw("lsame");

        gepsilon2_rx[istations]->SetLineColor(icolors[istations]);
        gepsilon2_rx[istations]->SetLineWidth(2);
        gepsilon2_rx[istations]->SetLineStyle(kDashed);

        gepsilon1_rx[istations]->SetName(snames[istations].c_str());
        legend26->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");

    }
    legend26->Draw("same");
    sname=sdir+"epsilons_sidebyside.pdf";
    c24->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot epsilon differences:
    //   left  panel: epsilon_1^T - epsilon_2^T
    //   right panel: epsilon_1^R - epsilon_2^R
    //
    // These directly show the angular separation between the two eigenmodes
    // in the chosen epsilon representation.
    // ------------------------------------------------------------------------
    TCanvas *c24b=new TCanvas("c24b","c24b",1600,800);
    c24b->Divide(2,1);

    c24b->cd(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.17);
    gPad->SetRightMargin(0.01);
    TH2D *h24b=new TH2D("","",100,-1600.,-600.,100,-0.15,0.15);
    titles(h24b, "", "Pulser height (m)", "#epsilon_{ 1}^{ T} - #epsilon_{ 2}^{ T} (^{o})");

    h24b->GetXaxis()->SetNdivisions(504);
    h24b->GetYaxis()->SetNdivisions(504);

    h24b->GetXaxis()->SetTitleOffset(1.2);
    h24b->GetYaxis()->SetTitleOffset(1.7);

    h24b->Draw();
    auto legend24b = new TLegend(0.18,0.2,0.36,0.64);
    auto legend26b = new TLegend(0.67,0.52,0.93,0.89);
    auto legend26c = new TLegend(0.67,0.12,0.93,0.22);
    legend24b->SetBorderSize(0);
    legend24b->SetTextSize(0.04);
    legend26b->SetBorderSize(0);
    legend26b->SetTextSize(0.05);
    legend26c->SetBorderSize(0);
    legend26c->SetTextSize(0.05);

    for (int istations=minstation;istations<=maxstation;istations++) {

        gdiffepsilon_tx[istations]->SetLineColor(icolors[istations]);
        gdiffepsilon_tx[istations]->SetLineWidth(2);
        gdiffepsilon_tx[istations]->Draw("lsame");

        legend24b->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");
        if (istations<5){
            legend26b->AddEntry(gdiffepsilon_rx[istations],snames[istations].c_str(),"l");
        }
        else if (istations==5){
            legend26c->AddEntry(gdiffepsilon_rx[istations],snames[istations].c_str(),"l");
        }
    }

    c24b->cd(2);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.17);
    gPad->SetRightMargin(0.01);
    TH2D *h26b=new TH2D("","",100,-1600.,-600.,100,-0.05,0.05);

    titles(h26b, "", "Pulser height (m)", "#epsilon_{ 1}^{ R} - #epsilon_{ 2}^{ R} (^{o})");

    h26b->GetXaxis()->SetNdivisions(504);
    h26b->GetYaxis()->SetNdivisions(504);

    h26b->GetYaxis()->SetTitleOffset(1.70);
    h26b->GetXaxis()->SetTitleOffset(1.2);

    h26b->Draw();

    for (int istations=minstation;istations<=maxstation;istations++) {

        gdiffepsilon_rx[istations]->SetLineColor(icolors[istations]);
        gdiffepsilon_rx[istations]->SetLineWidth(2);
        gdiffepsilon_rx[istations]->Draw("lsame");
    }
    legend26b->Draw("same");
    legend26c->Draw("same");
    sname=sdir+"diffepsilons_sidebyside.pdf";
    c24b->Print(sname.c_str());



    // ------------------------------------------------------------------------
    // Plot dot products between:
    //   panel 1: S-hat vectors of the two modes
    //   panel 2: E-hat vectors of the two modes
    //   panel 3: D-hat / p vectors of the two modes
    //
    // These are orthogonality / alignment diagnostics versus pulser depth.
    // ------------------------------------------------------------------------
    TCanvas *c27=new TCanvas("c27","c27",800,800);
    TH2D *h27a=new TH2D("h27a","h27a",100,-1800.,0.,100,-1.,1.);
    TH2D *h27b=new TH2D("h27b","h27b",100,-1800.,0.,100,-1.,1.);
    TH2D *h27c=new TH2D("h27c","h27c",100,-1800.,0.,100,-1.,1.);

    c27->Divide(1,3);
    c27->cd(1);
    h27a->Draw();
    for (int istations=minstation;istations<=maxstation;istations++) {
        gdotShats_tx[istations]->SetLineColor(icolors[istations]);
        gdotShats_tx[istations]->SetLineWidth(2);
        gdotShats_tx[istations]->Draw("lsame");
    }
    c27->cd(2);
    h27b->Draw();
    for (int istations=minstation;istations<=maxstation;istations++) {
        gdotEhats_tx[istations]->SetLineColor(icolors[istations]);
        gdotEhats_tx[istations]->SetLineWidth(2);
        gdotEhats_tx[istations]->Draw("lsame");
    }
    c27->cd(3);
    h27c->Draw();
    for (int istations=minstation;istations<=maxstation;istations++) {
        gdotDhats_tx[istations]->SetLineColor(icolors[istations]);
        gdotDhats_tx[istations]->SetLineWidth(2);
        gdotDhats_tx[istations]->Draw("lsame");
    }
    sname=sdir+"dotproducts.pdf";
    c27->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot path depth vs internal path-step index for the special ray path.
    // This is a diagnostic of the discretized ray-tracing path.
    // ------------------------------------------------------------------------
    TCanvas *c28=new TCanvas("c28","c28",800,800);
    TH2D *h28a=new TH2D("h28a","h28a",100,-1800.,0.,100,0.,100000.);
    h28a->Draw();
    for (int istations=minstation;istations<=maxstation;istations++) {
        g_depth_istep[istations]->Draw("lsame");
    }
    sname=sdir+"depth_idepth.pdf";
    c28->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot TX beam factors for the two modes.
    // Top panels    : beam1 at TX
    // Bottom panels : beam2 at TX
    // ------------------------------------------------------------------------
    TCanvas *c29=new TCanvas("c29","c29",800,800);
    TH2D *h29a=new TH2D("h29","h29",100,-1800.,0.,100,0.,1.);
    TH2D *h29b=new TH2D("h29","h29",100,-1800.,0.,100,0.,1.);
    c29->Divide(3,4);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c29->cd(istations+1);
        h29a->Draw();
        gtxdepth_beam1[istations]->SetLineColor(icolors[istations]);
        gtxdepth_beam1[istations]->SetLineWidth(2);
        gtxdepth_beam1[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c29->cd(NSTATIONS+istations+1);
        h29b->Draw();
        gtxdepth_beam2[istations]->SetLineColor(icolors[istations]);
        gtxdepth_beam2[istations]->SetLineWidth(2);
        gtxdepth_beam2[istations]->Draw("lsame");
    }
    sname=sdir+"beams_tx.pdf";
    c29->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot RX beam factors for the two modes.
    // Top panels    : beam1 at RX
    // Bottom panels : beam2 at RX
    // ------------------------------------------------------------------------
    TCanvas *c30=new TCanvas("c30","c30",800,800);
    TH2D *h30a=new TH2D("h30","h30",100,-1800.,0.,100,0.,1.);
    TH2D *h30b=new TH2D("h30","h30",100,-1800.,0.,100,0.,1.);
    c30->Divide(3,4);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c30->cd(istations+1);
        h30a->Draw();
        grxdepth_beam1[istations]->SetLineColor(icolors[istations]);
        grxdepth_beam1[istations]->SetLineWidth(2);
        grxdepth_beam1[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c30->cd(NSTATIONS+istations+1);
        h30b->Draw();
        grxdepth_beam2[istations]->SetLineColor(icolors[istations]);
        grxdepth_beam2[istations]->SetLineWidth(2);
        grxdepth_beam2[istations]->Draw("lsame");
    }
    sname=sdir+"beams_rx.pdf";
    c30->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot attenuation-only and attenuation×beam factors vs pulser depth.
    // ------------------------------------------------------------------------
    TCanvas *c31=new TCanvas("c31","c31",800,800);
    TH2D *h31=new TH2D("h31","h31",100,-1800.,0.,100,1.E-5,1.);
    c31->SetLogy();
    h31->Draw();
    for (int istations=minstation;istations<=maxstation;istations++) {
        grxdepth_atten[istations]->SetLineColor(kGray);
        grxdepth_atten[istations]->SetLineWidth(2);
        grxdepth_atten[istations]->Draw("lsame");

        grxdepth_atten_beam[istations]->SetLineColor(kGray+1);
        grxdepth_atten_beam[istations]->SetLineWidth(2);
        grxdepth_atten_beam[istations]->SetLineStyle(kDashed);
        grxdepth_atten_beam[istations]->Draw("lsame");


    }

    sname=sdir+"dispersionangles.pdf";
    c31->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot RX-side mode angles in the ordinary "clock" convention.
    // ------------------------------------------------------------------------
    TCanvas *c32=new TCanvas("c32","c32",800,800);
    TH2D *h32a=new TH2D("h32","h32",100,-1800.,0.,100,85.,95.);
    TH2D *h32b=new TH2D("h32","h32",100,-1800.,0.,100,-5.,5.);
    c32->Divide(3,4);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c32->cd(istations+1);
        h32a->Draw();
        grxdepth_theta1[istations]->SetLineColor(icolors[istations]);
        grxdepth_theta1[istations]->SetLineWidth(2);
        grxdepth_theta1[istations]->Draw("lsame");
        grxdepthE_theta1[istations]->SetLineColor(icolors[istations]);
        grxdepthE_theta1[istations]->SetLineWidth(2);
        grxdepthE_theta1[istations]->SetLineStyle(kDashed);
        grxdepthE_theta1[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c32->cd(NSTATIONS+istations+1);
        h32b->Draw();
        grxdepth_theta2[istations]->SetLineColor(icolors[istations]);
        grxdepth_theta2[istations]->SetLineWidth(2);
        grxdepth_theta2[istations]->Draw("lsame");
        grxdepthE_theta2[istations]->SetLineColor(icolors[istations]);
        grxdepthE_theta2[istations]->SetLineWidth(2);
        grxdepthE_theta2[istations]->SetLineStyle(kDashed);
        grxdepthE_theta2[istations]->Draw("lsame");
    }
    sname=sdir+"anglesontheclock_rx.pdf";
    c32->Print(sname.c_str());


    // ------------------------------------------------------------------------
    // Same as c32, but using the S-clock convention at the receiver.
    // ------------------------------------------------------------------------
    TCanvas *c33=new TCanvas("c33","c33",800,800);
    TH2D *h33a=new TH2D("h33","h33",100,-1800.,0.,100,0.,180.);
    TH2D *h33b=new TH2D("h33","h33",100,-1800.,0.,100,-90.,90.);
    c33->Divide(3,4);

    for (int istations=minstation;istations<=maxstation;istations++) {
        c33->cd(istations+1);
        h33a->Draw();
        grxdepth_theta1_Sclock[istations]->SetLineColor(icolors[istations]);
        grxdepth_theta1_Sclock[istations]->SetLineWidth(2);
        grxdepth_theta1_Sclock[istations]->Draw("lsame");
        grxdepthE_theta1_Sclock[istations]->SetLineColor(icolors[istations]);
        grxdepthE_theta1_Sclock[istations]->SetLineWidth(2);
        grxdepthE_theta1_Sclock[istations]->SetLineStyle(kDashed);
        grxdepthE_theta1_Sclock[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c33->cd(NSTATIONS+istations+1);

        grxdepth_theta2_Sclock[istations]->SetLineColor(icolors[istations]);
        grxdepth_theta2_Sclock[istations]->SetLineWidth(2);

        grxdepthE_theta2_Sclock[istations]->SetLineColor(icolors[istations]);
        grxdepthE_theta2_Sclock[istations]->SetLineWidth(2);
        grxdepthE_theta2_Sclock[istations]->SetLineStyle(kDashed);

        h33b->Draw();
        grxdepth_theta2_Sclock[istations]->Draw("lsame");
        grxdepthE_theta2_Sclock[istations]->Draw("lsame");
    }
    sname=sdir+"anglesontheSclock_rx.pdf";
    c33->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot total power in the two receiver channels.
    // Top panels    : r^2 P_theta-like quantity (r1)
    // Bottom panels : r^2 P_phi-like quantity   (r2)
    // ------------------------------------------------------------------------
    TCanvas *c34=new TCanvas("c34","c34",800,800);
    TH2D *h34a[NSTATIONS];
    TH2D *h34b[NSTATIONS];

    for (int i=minstation;i<=maxstation;i++) {
        h34a[i]=new TH2D("h34a","h34a",200,-1800.,0.,100,0.,pmax[0][i]);
        h34a[i]->SetXTitle("SpiceCore pulser height (m)");
        h34a[i]->SetYTitle("r^{2} P_{#theta} (arb. units)");
        h34b[i]=new TH2D("h34b","h34b",200,-1800.,0.,100,0.,pmax[1][i]);
        h34b[i]->SetXTitle("SpiceCore pulser height (m)");
        h34b[i]->SetYTitle("r^{2} P_{#phi} (arb. units)");
    }

    c34->Divide(3,4);
    for (int istations=minstation;istations<=maxstation;istations++) {
        c34->cd(istations+1);
        h34a[istations]->Draw();
        gpower_r1[istations]->SetLineColor(icolors[istations]);
        gpower_r1[istations]->SetLineWidth(2);
        gpower_r1[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c34->cd(NSTATIONS+istations+1);
        h34b[istations]->Draw();
        gpower_r2[istations]->SetLineColor(icolors[istations]);
        gpower_r2[istations]->SetLineWidth(2);
        gpower_r2[istations]->Draw("lsame");
    }
    sname=sdir+"powers.pdf";
    c34->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot square-rooted constructive/destructive voltage envelopes.
    // These are amplitude-like versions of the envelope quantities.
    // ------------------------------------------------------------------------
    TCanvas *c35=new TCanvas("c35","c35",800,800);
    TH2D *h35a[NSTATIONS];
    TH2D *h35b[NSTATIONS];

    for (int i=minstation;i<=maxstation;i++) {
        h35a[i]=new TH2D("h35a","h35a",200,-1500.,-850,100,0.,vmax[0][i]);
        h35a[i]->SetXTitle("SpiceCore pulser height (m)");
        h35a[i]->SetYTitle("r^{2} P_{#theta} (arb. units)");
        h35b[i]=new TH2D("h35b","h35b",200,-1500.,-850,100,0.,vmax[1][i]);
        h35b[i]->SetXTitle("SpiceCore pulser height (m)");
        h35b[i]->SetYTitle("r^{2} P_{#phi} (arb. units)");
    }

    c35->Divide(3,4);
    for (int istations=minstation;istations<=maxstation;istations++) {
        c35->cd(istations+1);
        h35a[istations]->Draw();

        gvenvelope_plus_r1[istations]->SetLineColor(icolors[istations]);
        gvenvelope_plus_r1[istations]->SetLineWidth(2);
        gvenvelope_plus_r1[istations]->Draw("lsame");
        gvenvelope_minus_r1[istations]->SetLineColor(icolors[istations]);
        gvenvelope_minus_r1[istations]->SetLineWidth(2);
        gvenvelope_minus_r1[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c35->cd(NSTATIONS+istations+1);

        h35b[istations]->Draw();
        gvenvelope_plus_r2[istations]->SetLineColor(icolors[istations]);
        gvenvelope_plus_r2[istations]->SetLineWidth(2);
        gvenvelope_plus_r2[istations]->Draw("lsame");
        gvenvelope_minus_r2[istations]->SetLineColor(icolors[istations]);
        gvenvelope_minus_r2[istations]->SetLineWidth(2);
        gvenvelope_minus_r2[istations]->Draw("lsame");
    }
    sname=sdir+"venvelopes.pdf";
    c35->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Compare attenuation-only, attenuation×beam, envelope, and final voltage.
    // Top panels    : r1/theta-like channel
    // Bottom panels : r2/phi-like channel
    // ------------------------------------------------------------------------
    TCanvas *c36=new TCanvas("c36","c36",800,800);
    TH2D *h36a[NSTATIONS];
    TH2D *h36b[NSTATIONS];

    for (int i=minstation;i<=maxstation;i++) {
        h36a[i]=new TH2D("","",200,-1600.,-600.,100,0.,vmax[0][i]);
        h36a[i]->SetXTitle("Pulser height (m)");
        h36a[i]->SetYTitle("Voltage @ 300 MHz (arb. units)");
        h36a[i]->GetXaxis()->SetTitleOffset(0.8);
        h36a[i]->GetYaxis()->SetTitleOffset(0.8);
        h36a[i]->GetXaxis()->SetTitleSize(0.05);
        h36a[i]->GetYaxis()->SetTitleSize(0.05);
        h36a[i]->GetXaxis()->SetNdivisions(504);

        h36b[i]=new TH2D("h36b","h36b",200,-1600.,-600.,100,0.,vmax[1][i]);
        h36b[i]->SetXTitle("SpiceCore pulser height (m)");
        h36b[i]->SetYTitle("sqrt(Power)");
    }

    c36->Divide(3,4);
    for (int istations=minstation;istations<=maxstation;istations++) {
        c36->cd(istations+1);
        h36a[istations]->Draw();
        grxdepth_atten[istations]->Draw("lsame");
        grxdepth_atten_beam[istations]->Draw("lsame");
        gvenvelope_plus_r1[istations]->SetLineColor(icolors[istations]);
        gvenvelope_plus_r1[istations]->SetLineStyle(kDashed);
        gvenvelope_plus_r1[istations]->Draw("lsame");
        gvoltage_r1[istations]->SetLineColor(icolors[istations]);
        gvoltage_r1[istations]->SetLineWidth(2);
        gvoltage_r1[istations]->Draw("lsame");
    }

    for (int istations=minstation;istations<=maxstation;istations++) {
        c36->cd(NSTATIONS+istations+1);
        h36b[istations]->Draw();
        grxdepth_atten[istations]->Draw("lsame");
        grxdepth_atten_beam[istations]->Draw("lsame");
        gvenvelope_plus_r2[istations]->SetLineColor(icolors[istations]);
        gvenvelope_plus_r2[istations]->SetLineStyle(kDashed);
        gvenvelope_plus_r2[istations]->Draw("lsame");
        gvoltage_r2[istations]->SetLineColor(icolors[istations]);
        gvoltage_r2[istations]->SetLineWidth(2);
        gvoltage_r2[istations]->Draw("lsame");
    }
    sname=sdir+"voltages.pdf";
    c36->Print(sname.c_str());

    // Special zoomed version for A5 (station index 4)
    TCanvas *c36_a5=new TCanvas("c36_a5","c36_a5",800,800);
    h36a[4]->Draw();
    h36a[4]->GetYaxis()->SetRangeUser(0.,25.);
    grxdepth_atten[4]->Draw("lsame");
    grxdepth_atten_beam[4]->Draw("lsame");
    gvenvelope_plus_r1[4]->Draw("lsame");
    gvoltage_r1[4]->Draw("lsame");

    sname=sdir+"voltages_a5.pdf";
    c36_a5->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Side-by-side voltage plots:
    //   left  panel: theta-pol / r1
    //   right panel: phi-pol   / r2
    //
    // Each panel overlays constructive/destructive envelopes and the final
    // interference-modified voltage at 300 MHz.
    // ------------------------------------------------------------------------
    TCanvas *c37=new TCanvas("c37","c37",1600,800);

    // Right panel: phi-polarized channel
    c37->Divide(2,1);
    c37->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetLogy();

    string stext="{#phi}-pol";
    TText *texttemp_hpol=new TText(-1000.,160.,stext.c_str());
    texttemp_hpol->SetTextFont(42);
    texttemp_hpol->SetTextSize(0.08);
    
    TH2D *h37=new TH2D("","",100,-1600.,-600.,100,0.11,200.);
    titles(h37, "", "Pulser height (m)", "r V_{#phi} (arb. units)");
    
    auto legend37 = new TLegend(0.22,0.52,0.36,0.9);
    legend37->SetBorderSize();
    legend37->SetTextSize(0.05);

    h37->GetYaxis()->SetTitleOffset(1.2);

    auto legend38b = new TLegend(0.17,0.76,0.62,0.89);
    legend38b->SetBorderSize();
    legend38b->SetTextSize(0.05);

    string snametemp="Voltage @ 300 MHz";
    h37->Draw();

    h37->GetXaxis()->SetNdivisions(504);

    for (int istations=minstation;istations<=maxstation;istations++) {

        gvenvelope_plus_r2[istations]->SetLineStyle(kSolid);
        gvenvelope_minus_r2[istations]->SetLineStyle(kSolid);
        gvoltage_r2[istations]->SetLineStyle(kDashed);

        // Skip ARIANNA in these overlays
        if (istations!=5) {
            gvenvelope_minus_r2[istations]->Draw("lsame");
            gvenvelope_plus_r2[istations]->Draw("lsame");
            gvoltage_r2[istations]->Draw("lsame");
        }

        gvoltage_r2[istations]->SetName(snames[istations].c_str());
        legend37->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");

        if (istations==0) {
            legend38b->AddEntry(gvenvelope_plus_r2[istations],"Voltage envelope","l");
            legend38b->AddEntry(gvoltage_r2[istations],snametemp.c_str(),"l");
        }
    }

    legend38b->Draw("same");

    // Left panel: theta-polarized channel
    c37->cd(1);

    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetLogy();

    TH2D *h38=new TH2D("","",100,-1600.,-600.,100,1.,100.);
    titles(h38, "", "Pulser height (m)", "r V_{#theta} (arb. units)");

    auto legend38 = new TLegend(0.62,0.2,0.88,0.5);
    legend38->SetBorderSize();
    legend38->SetTextSize(0.05);

    h38->GetYaxis()->SetTitleOffset(1.35);
    h38->GetXaxis()->SetNdivisions(504);

    h38->Draw();
    stext="{#theta}-pol";
    TText* texttemp_vpol=new TText(-1000.,160.,stext.c_str());
    texttemp_vpol->SetTextFont(42);

    // Optional axis helper for alternate scaling
    TF1 *fy2=new TF1("fy2","x",0.,200.);
    TGaxis *A2 = new TGaxis(-600.,0.,-600.,200.,"fy2",510,"L+");
    A2->SetTitle("r V_{#theta} @ 300 MHz (arb. units)");
    A2->SetTextFont(42);
    A2->SetLabelFont(42);
    A2->SetTitleSize(0.05);

    for (int istations=minstation;istations<=maxstation;istations++) {

        gvenvelope_plus_r1[istations]->SetLineStyle(kSolid);
        gvenvelope_plus_r1[istations]->SetName(snames[istations].c_str());

        gvenvelope_minus_r1[istations]->SetLineStyle(kSolid);
        gvoltage_r1[istations]->SetLineStyle(kDashed);

        if (istations==0) {
            gvoltage_r1[istations]->SetName(snametemp.c_str());
        }

        if (istations!=5) {
            gvoltage_r1[istations]->Draw("lsame");
            gvenvelope_minus_r1[istations]->Draw("lsame");
            gvenvelope_plus_r1[istations]->Draw("lsame");
        }

        legend38->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");
    }

    legend38->Draw("same");

    sname=sdir+"HPolVPolvoltages_sidebyside.pdf";
    c37->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Same as c37, but for electric-field amplitudes instead of voltages.
    // ------------------------------------------------------------------------
    TCanvas *c37b=new TCanvas("c37b","c37b",1600,800);

    c37b->Divide(2,1);

    // Right panel: phi-like electric field
    c37b->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetLogy();

    TH2D *h37b=new TH2D("","",100,-1600.,-600.,100,0.11,300.);
    titles(h37b, "", "Pulser height (m)", "r |E_{#phi}| (arb. units)");
    auto legend37b = new TLegend(0.22,0.52,0.36,0.9);
    legend37b->SetBorderSize();
    legend37b->SetTextSize(0.05);

    h37b->GetYaxis()->SetTitleOffset(1.2);

    auto legend38c = new TLegend(0.17,0.76,0.62,0.89);
    legend38c->SetBorderSize();
    legend38c->SetTextSize(0.05);

    string snametempb="E field @ 300 MHz";
    h37b->Draw();

    h37b->GetXaxis()->SetNdivisions(504);

    for (int istations=minstation;istations<=maxstation;istations++) {

        gEenvelope_plus_r2[istations]->SetLineStyle(kSolid);
        gEenvelope_plus_r2[istations]->SetLineColor(icolors[istations]);
        gEenvelope_plus_r2[istations]->SetLineWidth(2);

        gEenvelope_plus_r2[istations]->Draw("lsame");

        gEenvelope_minus_r2[istations]->SetLineStyle(kSolid);
        gEenvelope_minus_r2[istations]->SetLineColor(icolors[istations]);
        gEenvelope_minus_r2[istations]->SetLineWidth(2);

        gEenvelope_minus_r2[istations]->Draw("lsame");

        gfield_r2[istations]->SetLineStyle(kDashed);
        gfield_r2[istations]->SetLineColor(icolors[istations]);
        gfield_r2[istations]->SetLineWidth(2);
        gfield_r2[istations]->Draw("lsame");
        gfield_r2[istations]->SetName(snames[istations].c_str());
        legend37b->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");

        if (istations==0) {
            legend38c->AddEntry(gEenvelope_plus_r2[istations],"E field envelope","l");
            legend38c->AddEntry(gfield_r2[istations],snametempb.c_str(),"l");
        }

    }

    legend38c->Draw("same");

    // Left panel: theta-like electric field
    c37b->cd(1);

    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    gPad->SetLogy();

    TH2D *h38b=new TH2D("","",100,-1600.,-600.,100,1.,100.);
    titles(h38b, "", "Pulser height (m)", "r |E_{#theta}| (arb. units)");

    auto legend38d = new TLegend(0.62,0.2,0.88,0.5);
    legend38d->SetBorderSize();
    legend38d->SetTextSize(0.05);

    h38b->GetYaxis()->SetTitleOffset(1.35);
    h38b->GetXaxis()->SetNdivisions(504);

    h38b->Draw();
    stext="{#theta}-pol";

    TF1 *fy2b=new TF1("fy2b","x",0.,200.);
    TGaxis *A2b = new TGaxis(-600.,0.,-600.,200.,"fy2b",510,"L+");
    A2b->SetTitle("r |E_{#theta}| @ 300 MHz (arb. units)");
    A2b->SetTextFont(42);
    A2b->SetLabelFont(42);
    A2b->SetTitleSize(0.05);

    for (int istations=minstation;istations<=maxstation;istations++) {

        gEenvelope_plus_r1[istations]->SetLineColor(icolors[istations]);
        gEenvelope_plus_r1[istations]->SetLineWidth(2);
        gEenvelope_plus_r1[istations]->SetLineStyle(kSolid);
        gEenvelope_plus_r1[istations]->SetName(snames[istations].c_str());
        gEenvelope_plus_r1[istations]->Draw("lsame");

        gEenvelope_minus_r1[istations]->SetLineColor(icolors[istations]);
        gEenvelope_minus_r1[istations]->SetLineWidth(2);
        gEenvelope_minus_r1[istations]->SetLineStyle(kSolid);
        gEenvelope_minus_r1[istations]->Draw("lsame");

        gfield_r1[istations]->SetLineStyle(kDashed);
        gfield_r1[istations]->SetLineColor(icolors[istations]);
        gfield_r1[istations]->SetLineWidth(2);

        if (istations==0) {
            gfield_r1[istations]->SetName(snametemp.c_str());
        }
        gfield_r1[istations]->Draw("lsame");

        legend38d->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");
    }

    legend38d->Draw("same");

    sname=sdir+"HPolVPolfields_sidebyside.pdf";
    c37b->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Plot receiver polarization angles Omega and Psi for all stations.
    // Top panels    : Omega
    // Bottom panels : Psi
    // ------------------------------------------------------------------------
    TCanvas *c40=new TCanvas("c40","c40",800,800);

    TH2D *h40a=new TH2D("","",100,-1800.,0.,100,0.,90.);
    TH2D *h40b=new TH2D("","",100,-1800.,0.,100,-90.,90.);

    c40->Divide(3,4);

    for (int istations=0;istations<NSTATIONS;istations++) {
        c40->cd(istations+1);
        h40a->Draw();

        gpolarization_Omega_rx[istations]->SetLineColor(icolors[istations]);
        gpolarization_Omega_rx[istations]->SetLineStyle(kSolid);
        gpolarization_Omega_rx[istations]->Draw("lsame");

    }
    for (int istations=0;istations<NSTATIONS;istations++) {
        c40->cd(NSTATIONS+istations+1);
        h40b->Draw();
        gpolarization_Psi_rx[istations]->SetLineColor(icolors[istations]);
        gpolarization_Psi_rx[istations]->SetLineStyle(kSolid);
        gpolarization_Psi_rx[istations]->Draw("lsame");

    }

    sname=sdir+"polarizations.pdf";
    c40->Print(sname.c_str());


    // ------------------------------------------------------------------------
    // Special ARIANNA-only plot of Psi versus positive pulser depth.
    // Uses the reversed-depth graph so the x-axis runs from ~800 to 1700 m.
    // ------------------------------------------------------------------------
    TCanvas *c41=new TCanvas("c41","c41",800,800);

    TH2D *h41=new TH2D("","",100,800.,1700.,100,0.,30.);
    h41->GetXaxis()->SetTitle("Pulser depth (m)");
    h41->GetYaxis()->SetTitle("Polarization angle #Psi (^{o})");
    h41->GetXaxis()->SetTitleOffset(1.1);
    h41->GetYaxis()->SetTitleOffset(1.45);
    h41->Draw();

    gpolarization_reversedepth_Psi_rx[5]->SetLineColor(kBlack);
    gpolarization_reversedepth_Psi_rx[5]->SetLineStyle(kSolid);
    gpolarization_reversedepth_Psi_rx[5]->SetLineWidth(3);
    gpolarization_reversedepth_Psi_rx[5]->Draw("lsame");

    sname=sdir+"polarization_arianna.pdf";
    c41->Print(sname.c_str());

    // ------------------------------------------------------------------------
    // Prepare output ROOT file that stores selected summary graphs for later use.
    // The filename encodes the TX and RX cross-polarization settings.
    // ------------------------------------------------------------------------
    char name[100];

    string sfilenames=sdir+"myoutputs_tx" + to_string(CROSSPOLANGLE_TX_INT) + "_rx" + to_string(CROSSPOLANGLE_RX_INT) + ".root";

    TFile *fout=new TFile(sfilenames.c_str(),"RECREATE");

    // ------------------------------------------------------------------------
    // Write a selected subset of output graphs to the ROOT file.
    //
    // The stored objects include:
    //   - receiver polarization-vs-depth summaries
    //   - voltage and field amplitudes
    //   - attenuation-only and attenuation×beam factors
    //   - constructive envelopes
    //   - selected spectra at "special" pulser depths
    //   - interference building blocks (V1^2, V2^2, V1V2, etc.)
    //
    // Each object name is tagged by station index so the file can be reused
    // later without ambiguity.
    // ------------------------------------------------------------------------
    for (int istations=0;istations<NSTATIONS;istations++) {

        // Receiver polarization angle Psi vs reversed depth
        sprintf(name,"gpolarization_reversedepth_Psi_rx_%d",istations);
        gpolarization_reversedepth_Psi_rx[istations]->Write(name);
        
        // Same, but for field-based polarization instead of voltage-based        
        sprintf(name,"gEpolarization_reversedepth_Psi_rx_%d",istations);
        gEpolarization_reversedepth_Psi_rx[istations]->Write(name);

        // Final voltage amplitudes in receiver channels r1 and r2
        sprintf(name,"gvoltage_r1_%d",istations);
        gvoltage_r1[istations]->Write(name);

        sprintf(name,"gvoltage_r2_%d",istations);
        gvoltage_r2[istations]->Write(name);
        
        // Attenuation-only and attenuation×beam transfer factors
        sprintf(name,"grxdepth_atten_%d",istations);
        grxdepth_atten[istations]->Write(name);

        sprintf(name,"grxdepth_atten_beam_%d",istations);
        grxdepth_atten_beam[istations]->Write(name);

        // Constructive voltage envelope in channel r1
        sprintf(name,"gvenvelope_plus_r1_%d",istations);
        gvenvelope_plus_r1[istations]->Write(name);

        // Store spectra for a fixed list of "special" pulser depths
        for (int ispecial=0;ispecial<NSPECIAL;ispecial++) {
            sprintf(name,"g_spectra_%d_%d",istations,ispecial);
            g_spectra[istations][(int)g_idepth[istations]->Eval(whichspecial[ispecial])]->Write(name);
        }

        // Constructive voltage envelope in channel r2
        sprintf(name,"gvenvelope_plus_r2_%d",istations);
        gvenvelope_plus_r2[istations]->Write(name);

        // Final voltage in r2 is written again here (same object as above)
        sprintf(name,"gvoltage_r2_%d",istations);
        gvoltage_r2[istations]->Write(name);

        // Final field amplitude in r2
        sprintf(name,"gfield_r2_%d",istations);
        gfield_r2[istations]->Write(name);

        // Final field amplitude in r2
        sprintf(name,"g_atten_power_%d",istations);
        g_atten_power[istations]->Write(name);

        sprintf(name,"g_atten_beam_power_%d",istations);
        g_atten_beam_power[istations]->Write(name);

        // Final field amplitude in r2
        sprintf(name,"gV1squared_r1_%d",istations);
        gV1squared_r1[istations]->Write(name);

        sprintf(name,"gV2squared_r1_%d",istations);
        gV2squared_r1[istations]->Write(name);

        sprintf(name,"gV1V2_r1_%d",istations);
        gV1V2_r1[istations]->Write(name);

        // Final field amplitude in r2
        sprintf(name,"gV1V2_r2_%d",istations);
        gV1V2_r2[istations]->Write(name);

        sprintf(name,"gV1squared_r2_%d",istations);
        gV1squared_r2[istations]->Write(name);

        sprintf(name,"gV2squared_r2_%d",istations);
        gV2squared_r2[istations]->Write(name);

        // Final field amplitude in r2
        sprintf(name,"goppositeV1V2_r2_%d",istations);
        goppositeV1V2_r2[istations]->Write(name);

        sprintf(name,"goppositeV1V2_r1_%d",istations);
        goppositeV1V2_r1[istations]->Write(name);

        // Debug printout of channel-r2 interference terms at a reference depth
        int idepth_temp=g_idepth[istations]->Eval(-1000.);
        cout << "station, depth, V1squared_r2, V2squared_r2, V1V2_r2 are " << istations << "\t" << vV1squared_r2[istations][idepth_temp] << "\t" << vV2squared_r2[istations][idepth_temp] << "\t" << vV1V2_r2[istations][idepth_temp] << "\n";

    }

    // Close the output ROOT file
    fout->Close();

    return (0);
}



// --------------------------------------------------------------------------
// directionNextStep
//
// Given a refractive-index triplet nvec and a polarization-like direction p_e1,
// build a new direction vector by weighting the Cartesian components by n_i.
// The output is normalized before being returned.
//
// Note:
//   p_e2 is passed in but not used in the current implementation.
// --------------------------------------------------------------------------
TVector3 directionNextStep(vector<double>nvec, TVector3 p_e1, TVector3 p_e2) {

    p_e1.SetMag(1.);

    TVector3 temp1;
    temp1[0]=nvec[0]*p_e1[0];
    temp1[1]=nvec[1]*p_e1[1];
    temp1[2]=nvec[2]*p_e1[2];
    temp1.SetMag(1.);

    return temp1;
}


// --------------------------------------------------------------------------
// rotateE
//
// Rotate a vector E into the frame where x is aligned with the ice-flow axis,
// apply the diagonal dielectric tensor epsilon there, then rotate back to the
// original laboratory frame.
//
// In other words, this computes:
//   E_lab -> E_principal -> epsilon * E_principal -> back to lab
// --------------------------------------------------------------------------
TVector3 rotateE(TVector3 epsilon, double angle_iceflow, TVector3 E) {

    double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
        {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
        {0.,0.,1.}};

    TVector3 tempvec;
    for (int i=0;i<3;i++) {
        double sum=0.;
        for (int j=0;j<3;j++) {
            sum+=rotate_toxalongiceflow[i][j]*E[j];
        }
        tempvec[i]=sum;
    }
    E=tempvec;

    for (int i=0;i<3;i++) {
        tempvec[i]=epsilon[i]*E[i];
    }
    E=tempvec;

    double rotate_backtonormal[3][3];

    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) {
            rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
        }
    }

    for (int i=0;i<3;i++) {
        double sum=0.;

        for (int j=0;j<3;j++) {
            sum+=rotate_backtonormal[i][j]*E[j];
        }

        tempvec[i]=sum;
    }

    E=tempvec;

    return E;
}


// --------------------------------------------------------------------------
// getNewkandE
//
// Starting from:
//   - principal refractive indices nvec
//   - ice-flow direction angle
//   - mode refractive index n
//   - displacement-field-like vector D
//   - a guessed propagation direction kguess
//
// this function:
//
//   1. rotates D into the principal-axis frame
//   2. computes E = epsilon^{-1} D
//   3. constructs k from the relation k = -D + n^2 E
//   4. rotates both k and E back into the lab frame
//   5. enforces the sign of k to be aligned with kguess
//
// E is returned by reference, and k is returned as the function value.
// --------------------------------------------------------------------------
TVector3 getNewkandE(vector<double> nvec,double angle_iceflow,double n,TVector3 D,TVector3 kguess,
TVector3 &E) {

    double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
        {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
        {0.,0.,1.}};

    TVector3 tempvec;
    for (int i=0;i<3;i++) {
        double sum=0.;
        for (int j=0;j<3;j++) {
            sum+=rotate_toxalongiceflow[i][j]*D[j];
        }
        tempvec[i]=sum;
    }
    D=tempvec;

    vector<double> epsilon;
    epsilon.resize(3);

    vector<double> epsilon_inverse;
    epsilon_inverse.resize(3);

    epsilon[0]=nvec[0]*nvec[0];
    epsilon[1]=nvec[1]*nvec[1];
    epsilon[2]=nvec[2]*nvec[2];

    epsilon_inverse[0]=1/epsilon[0];
    epsilon_inverse[1]=1/epsilon[1];
    epsilon_inverse[2]=1/epsilon[2];

    E[0]=epsilon_inverse[0]*D[0];
    E[1]=epsilon_inverse[1]*D[1];
    E[2]=epsilon_inverse[2]*D[2];

    TVector3 kvec=-1.*D+n*n*E;

    if (kvec.Mag()<HOWSMALLISTOOSMALL)
    cout << "kvec mag is " << kvec.Mag() << "\n";

    kvec.SetMag(1.);

    double rotate_backtonormal[3][3];

    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) {
            rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
        }
    }

    TVector3 tempvec2;
    for (int i=0;i<3;i++) {
        double sum1=0.;
        double sum2=0.;

        for (int j=0;j<3;j++) {
            sum1+=rotate_backtonormal[i][j]*kvec[j];
            sum2+=rotate_backtonormal[i][j]*E[j];
        }

        tempvec[i]=sum1;
        tempvec2[i]=sum2;
    }

    kvec=tempvec;
    E=tempvec2;

    if (kvec.Dot(kguess)<0.) {
        kvec=-1.*kvec;
    }

    return kvec;
}


// --------------------------------------------------------------------------
// getDeltaN(alpha)
//
// Analytic toy-model expression for delta n as a function of angle alpha
// relative to the ice-flow axis, using NICE and DELTAN as global constants.
//
// This is the simplified scalar model used for comparison plots, not the more
// general path-dependent solver used elsewhere in the main calculation.
// --------------------------------------------------------------------------
double getDeltaN(double alpha) {

    double n_L_prime=(NICE-DELTAN/2.)/sqrt(1-DELTAN/NICE*cos(alpha)*cos(alpha)*(1-DELTAN/(4*NICE)));
    return 2*(NICE-n_L_prime);
}


// --------------------------------------------------------------------------
// switchThem
//
// Swap two angles in place.
// Used in special cases (e.g. ARIANNA / biaxial mode ordering conventions)
// when the two eigenmodes need to be relabeled.
// --------------------------------------------------------------------------
void switchThem(double &thetaE_e1_Sclock,double &thetaE_e2_Sclock) {

    double tempangle=thetaE_e1_Sclock;
    thetaE_e1_Sclock=thetaE_e2_Sclock;
    thetaE_e2_Sclock=tempangle;

}

// --------------------------------------------------------------------------
// Flipped
//
// Compare the current angle theta_e1 to its starting value theta_e1_start.
// If the change is close to pi/2 (after wrapping to [-pi,pi]), interpret this
// as a sign flip / branch flip and return -1. Otherwise return +1.
//
// This is used to keep track of sign changes in the eigenvector transport along
// the ray path, so the accumulated birefringent phase can be corrected by a
// factor of ±1.
// --------------------------------------------------------------------------
double Flipped(double theta_e1,double theta_e1_start) {

    double diff=theta_e1-theta_e1_start;
    if (diff>PI)
    diff=diff-2.*PI;
    if (diff<-1.*PI)
    diff=diff+2.*PI;

    double notflipped=1.;

    if (abs(diff)>PI/2.*0.9)
    notflipped=-1.;

    return notflipped;
}

// --------------------------------------------------------------------------
// makePretty2Panel
//
// Utility function to create a two-panel canvas with customized pad geometry,
// margins, and grid settings.
// --------------------------------------------------------------------------
TCanvas * makePretty2Panel(){
    TCanvas *ccc = new TCanvas("plotEvent","plotEvent", 800, 400);
    ccc->Divide(2, 0);

    ccc->GetPad(1)->SetPad(.005, .005, .4975, .995);
    ccc->GetPad(1)->Divide(0, 2);

    ccc->GetPad(2)->SetPad(.5025, .005, .995, .995);
    ccc->SetWindowSize(1200, 700);
    ccc->cd(1)->cd(1)->SetGrid();


    gPad->SetLeftMargin(.15);
    gPad->SetBottomMargin(.12);
    ccc->cd(1)->cd(2);


    gPad->SetBottomMargin(.12);
    gPad->SetRightMargin(.19);
    gPad->SetLeftMargin(.15);
    return ccc;
}

// --------------------------------------------------------------------------
// titles(TGraph*)
//
// Helper to apply a consistent style to TGraph axis labels and titles.
// --------------------------------------------------------------------------
void titles(TGraph *inGr, TString title, TString xtitle, TString ytitle){
    auto sizeT=.055;
    inGr->SetTitle(title);
    inGr->GetXaxis()->SetTitle(xtitle);
    inGr->GetYaxis()->SetTitle(ytitle);
    inGr->GetXaxis()->SetTitleSize(sizeT);
    inGr->GetYaxis()->SetTitleSize(sizeT);
    inGr->GetXaxis()->SetLabelSize(sizeT);
    inGr->GetYaxis()->SetLabelSize(sizeT);

    inGr->GetYaxis()->SetLabelOffset(.01);
    inGr->GetYaxis()->SetTitleOffset(1.2);
}


// --------------------------------------------------------------------------
// titles(TH2*)
//
// Same helper as above, but for TH2 histogram frames used as empty axes.
// --------------------------------------------------------------------------
void titles(TH2 *inH, TString title, TString xtitle, TString ytitle){
    auto sizeT=.055;
    inH->SetTitle(title);
    inH->GetXaxis()->SetTitle(xtitle);
    inH->GetYaxis()->SetTitle(ytitle);
    inH->GetXaxis()->SetTitleSize(sizeT);
    inH->GetYaxis()->SetTitleSize(sizeT);
    inH->GetXaxis()->SetLabelSize(sizeT);
    inH->GetYaxis()->SetLabelSize(sizeT);

    inH->GetYaxis()->SetLabelOffset(.01);
    inH->GetYaxis()->SetTitleOffset(1.2);
}
