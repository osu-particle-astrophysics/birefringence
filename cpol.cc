/////////////////////////////////////////////////////////////////////
// C++ example program to demonstrate usage of the MSTW 2008 PDFs. //
// Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>.               //
/////////////////////////////////////////////////////////////////////

// things to try
// global fit
// try changing axes to see if it could be interference between other two
// ray bending

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
//#include "y.hh" 
#include <iostream>
#include <fstream>
#include <string>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TVector2.h"
//#include "units.h"
//#include "raytracing_tools.h"
#include "IceRayTracing/IceRayTracing.C"
#include "TLegend.h"
#include "birefringence.hh"

using namespace std;

const double PI=3.1415926;
const double CLIGHT=3.E8;
const double DEGRAD=180./PI;
const double NICE=1.78;
//const double DELTAN=35.E-9/20.E-6;
//const double DELTAN=35.E-9/20.E-6;
const double DELTAN=0.0005; // for ice flow
const double HOWSMALLISTOOSMALL=1.e-8;
const int UZAIRSTEP=5;




//const double CROSSPOLANGLE_TX=0.;
//const double CROSSPOLANGLE_RX=0.;

vector<double> nvec{1.77,1.7805,1.7815};

//vector<double> nvec{1.77,1.77,1.7815};
//const double DELTAN=0.003; // for vertical compression
//const double L_ATTEN=690.; // meters, from Eugene's thesis.
//const double L_ATTEN=1660.; // from ARA prototype paper.  1660+255/-120 m @300 MHz
const double L_ATTEN=2200.; // from ARA prototype paper.  1660+255/-120 m @300 MHz
const double MFT=1./100.*2.54/1.*12.; // (something in ft.)*1m/100cm*2.54cm/1in*12in

int icolors[6]={kGreen+2,kRed+1,kOrange+1,kViolet+1,kBlue,kBlack};
int icolors_dave[6]={kBlack,kRed,kGreen,kBlue,kYellow,kBlack};

void titles(TH2 *inH, TString title, TString xtitle, TString ytitle);
void titles(TGraph *inGr, TString title, TString xtitle, TString ytitle);
//double Flipped(TVector3 p_e1_previous, TVector3 p_e2_previous, TVector3 &p_e1, TVector3 &p_e2);
TCanvas *makePretty2Panel();
double Flipped(double theta_e1,double theta_e1_start);
void switchThem(double &thetaE_e1_Sclock,double &thetaE_e2_Sclock);
/*
void getManyAnglesontheClock(int BIAXIAL,double crosspolangle_tx,
			     TVector3 rhat_thisstep,
			     TVector3 p_e1,TVector3 p_e2,TVector3 E_e1,TVector3 E_e2,
			     double &theta_e1,double &theta_e2,double &thetaE_e1,double &thetaE_e2,
			     double &theta_e1_Sclock,double &theta_e2_Sclock,
			     double &thetaE_e1_Sclock,double &thetaE_e2_Sclock,
			     TVector3 &Shat_e1,TVector3 &Shat_e2,
			     double &E_e1_thetacomponent, double &E_e2_thetacomponent,
			     double &E_e1_phicomponent, double &E_e2_phicomponent);
*/
/*
void getAnglesontheClock(TVector3 rhat_thisstep1, TVector3 rhat_thisstep2,
			 TVector3 p_e1, TVector3 p_e2,
			 double &theta_e1,double &theta_e2,
			 double &p_e1_thetacomponent,double &p_e2_thetacomponent,
			 double &p_e1_phicomponent,double &p_e2_phicomponent);
/*
double getDeltaN(int BIAXIAL,vector<double> nvec,TVector3 rhat,double angle_iceflow, 
		 double &n_e1, double &n_e2,TVector3 &p_e1,TVector3 &p_e2);
*/

double getDeltaN(double alpha);
TVector3 getNewkandE(vector<double> nvec,double angle_iceflow,double n,TVector3 D,TVector3 kguess,
		     TVector3 &E);

TVector3 directionNextStep(vector<double>nvec, TVector3 E);

TVector3 rotateE(TVector3 epsilon, double angle_iceflow, TVector3 E);
//TVector3 rotateD(TVector3 epsilon, double angle_iceflow, TVector3 D);

//double getV(vector<double> nvec);

//void thetastoEpsilons(double thetaE_e1_Sclock,double thetaE_e2_Sclock,
//			double &epsilon1,double &epsilon2);


int main(int argc, char** argv) {  


  TVector3 inverse_epsilon;
  inverse_epsilon[0]=1/(nvec[0]*nvec[0]);
  inverse_epsilon[1]=1/(nvec[1]*nvec[1]);
  inverse_epsilon[2]=1/(nvec[2]*nvec[2]);

  TVector3 epsilon;
  epsilon[0]=nvec[0]*nvec[0];
  epsilon[1]=nvec[1]*nvec[1];
  epsilon[2]=nvec[2]*nvec[2];


  //  cout << "V is " << getV(nvec) << "\n";

  const double DEPTH=-600.; // take an average pulser depth for now
  const int ispecial=0;
  

  const double depth_special=-1100.;
  //const int depth_special=-850.;
  //const int depth_special=-1100.;
  
  const double stepspecial=20.;
  const int NSPECIAL=50;
  double whichspecial[NSPECIAL];
  double startspecial=-1600.;
  for (int ispecial=0;ispecial<NSPECIAL;ispecial++) {
    whichspecial[ispecial]=startspecial+stepspecial*(double)ispecial;
  }

  const int WHICHPOL=0;
  const int DORAYTRACING=1;
  
  const int NSTATIONS=6; // 5 ARA stations + ARIANNA SP
  const int minstation=0; // which to start with
  const int maxstation=5; // which to end with
  int igreatestdepth[NSTATIONS];
  int imostshallowdepth[NSTATIONS];
  double horizontal_distances[NSTATIONS]={1257.,2353.,3146.,3199.,5179.8892,653.804525};
  string snames[NSTATIONS]={"A1","A2","A3","A4","A5","ARIANNA"};
  // from geoff's github, for ARIANNA station 51:
  // r = 653.804524 # meters
  // seems to be distance to arianna from spicecore.

  //# Took azimuth measuremtn of 205 degrees with repsect to magnetic north (magnetic north is 32.7 degrees west of grid north; but we measure WRT grid east so 90 + 205 + 32.7)
    //actualTiltDirection = 327.7
  // this is called Direction_from_Spice
  //  phi = 312.448284

  // paper quotes 1.4deg off ice flow direction
  //  double PHI_ARIANNA=(90.+36.+46./60.+23./3600.)/DEGRAD;
  double PHI_ARIANNA=(90.+36.+46./60.+23./3600.+1.4)/DEGRAD;
  // 126.77 + 1.4

  const int NSTEPS=100.; // steps in depth of pulser
  double min_altitude=-1000.;  // minimum altitude (meaning it's a negative number)
  double max_altitude=-600.;
  
  double N_ICE=1.78; // index of ref of ice
  double DELTA_N=0.427; // parameters of ice model
  //double Z_0=0.71*utl::m; //meters
  //double Z_0=71.*utl::m; //meters
  
  static double VOLTAGENORM=150.*11./7.;

  //  cout << "index of refraction at 200m depth is " << Getnz(-200.) << "\n";
  //cout << "index of refraction at 1000m depth is " << Getnz(-1000.) << "\n"; 

  //  cout << "A2 distance at the bottom is " << sqrt((min_altitude-(-200.))*(min_altitude-(-200.))+horizontal_distances[1]*horizontal_distances[1]) << "\n";
  //cout << "A2 distance at the top is " << sqrt((max_altitude-(-200.))*(max_altitude-(-200.))+horizontal_distances[1]*horizontal_distances[1]) << "\n";


  char clswitch; // command line switch

  double freq=160.E6; // what frequency are we using 
  int CROSSPOLANGLE_TX_INT=7;
  int CROSSPOLANGLE_RX_INT=-7;
  int BIAXIAL=1; // if this is zero then it's uniaxial.  if it is -1 then it is isotropic.  the -1 option doesn't work yet. 
  int CONSTANTINDICATRIX=0; // try out making the indicatrix constant with depth and see which effects are still there.   default is zero.

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
      } // end switch
    } // end while
  } // end if arg>1

  double CROSSPOLANGLE_TX=(double)CROSSPOLANGLE_TX_INT/DEGRAD;
  double CROSSPOLANGLE_RX=(double)CROSSPOLANGLE_RX_INT/DEGRAD;




  double freqmin=0.;
  double freqmax=1.E9;
  const int NFREQ=100;
  vector<double> vfreqs;


  for (int i=0;i<NFREQ;i++) {
    vfreqs.push_back(freqmin+(freqmax-freqmin)/(double)NFREQ*(double)i);
  }

  //  double pulser_coords[2]={42500, 48815.78947368421}; // my guess off the map
  double pulser_coords[2]={42358.94,48974.2}; // dzb said this came from Leah Street.

  double station_coords[NSTATIONS][2]={ // this time the numbering is 5 ARA stations + ARIANNA station 51
    // these i had digitized
    //{38881.57894736842, 50921.05263157895}, // in feet
    //{35460.52631578947, 45328.94736842105},
    //{32302.631578947367, 51052.631578947374},
    //{35460.52631578947, 56710.52631578947},

    // these are reading off the surveyor's figure
    {38754., 51051.}, // in feet
    {35481., 45369.},
    {32200., 51053.},
    {35478., 56737.},
    {32356., 39746.}, // from Kaeli
    //{pulser_coords[0]+horizontal_distances[5]/MFT*cos(PHI_ARIANNA) ,
    //pulser_coords[1]+horizontal_distances[5]/MFT*sin(PHI_ARIANNA) }
    // these are from geoff's github.
    {41153.,50381.75}
  };

  cout << "station_coords of ARIANNA is " << station_coords[5][0] << "\t" << station_coords[5][1] << "\n";
  for (int i=0;i<NSTATIONS;i++) {
    cout << "distance to station " << i+1 << " is " << sqrt((station_coords[i][0]-pulser_coords[0])*(station_coords[i][0]-pulser_coords[0]) + (station_coords[i][1]-pulser_coords[1])*(station_coords[i][1]-pulser_coords[1]))*MFT << "\n"; 
    }
  for (int j=0;j<2;j++) {
    pulser_coords[j]=pulser_coords[j]*MFT;
  }  
  for (int i=0;i<6;i++) {
    for (int j=0;j<2;j++) {
      station_coords[i][j]=station_coords[i][j]*MFT;
    }
  }


  // ice flow is N36deg46'23''
  const double angle_iceflow=(36.+ (46./60.) + (23./3600.) + 90.)/DEGRAD; // in degrees, counter-clockwise from traditional phi=0.
  TVector3 ordinary(cos(angle_iceflow),sin(angle_iceflow),0.); // ordinary axis along the direction of ice flow
  //TVector3 ordinary(0.,0.,1.); // ordinary axis along the direction of "up"
  
  TVector3 tx_orientation(0.,0.,1.);
  //TVector3 tx_orientation(1.,0.,0.);
  TVector3 rx1_orientation(0.,0.,1.);
  
  double beam[2]={0.};
  double deltan_exp[6];
  double deltan_obs[6];
  double station[6]={1,2,3,4,5,6};
  double station_depths[6]={-80.,-180.,-180.,-180.,-180.,-1.}; // always forget which is 100m.
  double zeroes[6]={0.};
  double deltan_obs_err[6]={0.};
  
  double alpha[6];
  double alpha_deg[6];
  TVector3 pulsertostationhat_specialdepth[6];
  TVector3 pulsertostation[6];


  //  for (int i=0;i<NSTATIONS;i++) {
  for (int i=minstation;i<=maxstation;i++) {

    pulsertostationhat_specialdepth[i][0]=(station_coords[i][0]-pulser_coords[0]); // this will be a negative number
    pulsertostationhat_specialdepth[i][1]=(station_coords[i][1]-pulser_coords[1]); 
    pulsertostationhat_specialdepth[i][2]=station_depths[i]-(DEPTH);
    if (pulsertostationhat_specialdepth[i].Mag()<HOWSMALLISTOOSMALL)
      cout << "pulsertostationhat_specialdepth mag is " << pulsertostationhat_specialdepth[i].Mag() << "\n";

    pulsertostationhat_specialdepth[i].SetMag(1.);

    double y=station_coords[i][1]-pulser_coords[1];
    double x=station_coords[i][0]-pulser_coords[0];
    double angle_thisstation=atan2(y,x);
    //cout << "i, angle_iceflow, angle is " << i << "\t" << angle_iceflow*DEGRAD << "\t" << angle_thisstation*DEGRAD << "\n";
    double thisalpha=angle_thisstation-angle_iceflow; // in degrees
    //double thisalpha=acos(pulsertostationhat_specialdepth[i].Dot(ordinary)); // in degrees
    if (thisalpha>PI)
      thisalpha-=2.*PI;
    if (thisalpha<-1.*PI)
      thisalpha+=2.*PI;
    alpha[i]=thisalpha;
    cout << "A" << i+1 << ": angle_thisstation, angle_iceflow, alpha are " << angle_thisstation*DEGRAD << "\t" << angle_iceflow*DEGRAD << "\t" << alpha[i]*DEGRAD << ".\n";
    
    alpha_deg[i]=alpha[i]*DEGRAD;

    TVector3 p_e1;
    TVector3 p_e2;
    double n_e1;
    double n_e2;
    //    deltan_exp[i]=getDeltaN(alpha[i]);
    deltan_exp[i]=getDeltaN(BIAXIAL,nvec,pulsertostationhat_specialdepth[i],angle_iceflow,
			    n_e1,n_e2,p_e1,p_e2);
    //cout << "alpha, deltan_exp are " << alpha[i]*DEGRAD << "\t" << deltan_exp[i] << "\n";    
    //cout << "for alpha=PI/2, deltan is " << getDeltaN(PI/2.) << "\n";

    //    cout << "At a 1.5km pulser depth, I would expect oscillations in A" << i << " spectra every "  << 2*CLIGHT/(PI*deltan_exp[i]*sqrt(1500*1500+horizontal_distances[i]*horizontal_distances[i]))/1.E6 << " MHz.\n";

  } // end loop over 5 stations
  double x_arianna=650.*cos(angle_iceflow);
  double y_arianna=650.*sin(angle_iceflow);
  double z_arianna=0.;
  pulsertostationhat_specialdepth[5][0]=x_arianna-pulser_coords[0];
  pulsertostationhat_specialdepth[5][1]=y_arianna-pulser_coords[1];
  pulsertostationhat_specialdepth[5][2]=z_arianna-(DEPTH);

if (pulsertostationhat_specialdepth[5].Mag()<HOWSMALLISTOOSMALL)
      cout << "pulsertostationhat_specialdepth[5] mag is " << pulsertostationhat_specialdepth[5].Mag() << "\n";

  pulsertostationhat_specialdepth[5].SetMag(1.);

    double thisalpha=acos(pulsertostationhat_specialdepth[5].Dot(ordinary)); // in degrees
    if (thisalpha>PI)
      thisalpha-=2.*PI;
    if (thisalpha<-1.*PI)
      thisalpha+=2.*PI;
    alpha[5]=thisalpha;

    TVector3 p_e1;
    TVector3 p_e2;
    double n_e1;
    double n_e2;
  //  alpha[5]=1.4/DEGRAD; // from ARIANNA paper

    //    cout << "delta n from ARIANNA is " << getDeltaN(nvec,pulsertostationhat_specialdepth[5],angle_iceflow,n_e1,n_e2,p_e1,p_e2) << "\n";
    //cout << "At 1.5km pulser depth, I would expect oscillations in ARIANNA spectrum every "  << 2*CLIGHT/(PI*getDeltaN(nvec,pulsertostationhat_specialdepth[5],angle_iceflow,n_e1,n_e2,p_e1,p_e2)*sqrt(1500*1500+horizontal_distances[5]*horizontal_distances[5]))/1.E6 << " MHz.\n";


  string line;
  std::vector< std::vector<double> > videpth; // these are negative numbers
  std::vector< std::vector<double> > vdepth; // these are negative numbers

  std::vector< std::vector<double> > vdepth_data; // these are negative numbers
  std::vector< std::vector<double> > vdepth_data_err;

  std::vector< std::vector<double> > vreversedepth; // these are negative numbers
  



  std::vector< std::vector<double> > vreceiveangle;
  std::vector< std::vector<double> > vlaunchangle;
  std::vector< std::vector<double> > voutput6;
  std::vector< std::vector<double> > voutput7;
  std::vector< std::vector<double> > voutput8;


  std::vector< std::vector<double> > vtotal_distances; // these are negative numbers
  std::vector< std::vector<double> > vtotal_distances_err; // these are negative numbers

  std::vector< std::vector<double> > vsnrmax;
  std::vector< std::vector<double> > vsnrmax_err;

  std::vector< std::vector<double> > vangle_khat_0_khat_1_2;
  std::vector< std::vector<double> > vangle_khat_0_khat_2_2;

  std::vector< std::vector<double> > vangle_khat_1_2_khat_2_2;

  std::vector< std::vector<double> > vangle_Shat_e1_khat;
  std::vector< std::vector<double> > vangle_Shat_e2_khat;



  TVector3 raypath[2][6]; // for each ray, each station
  TVector3 raypath_n[2][6]; // for each ray, each station
  std::vector< std::vector< double > > vistep; // first element is station, second element is x, y, z, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector< std::vector<double> > > vraypos; // first element is station, second element is x, y, z, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vraypath_ne1; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vraypath_ne2; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.


  std::vector< std::vector<double> > vrxdepth_beam1; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vrxdepth_beam2; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vtxdepth_beam1; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vtxdepth_beam2; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vrxdepth_atten; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vrxdepth_notflipped; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vrxdepth_atten_beam; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vrxdepth_atten_power; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vrxdepth_atten_beam_power; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.


  std::vector< std::vector<double> > vtxdepth_theta1; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vtxdepth_theta2; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vrxdepth_theta1; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vrxdepth_theta2; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vtxdepth_theta1_Sclock; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vtxdepth_theta2_Sclock; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vrxdepth_theta1_Sclock; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vrxdepth_theta2_Sclock; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vtxdepthE_theta1; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vtxdepthE_theta2; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vrxdepthE_theta1; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vrxdepthE_theta2; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vtxdepthE_theta1_Sclock; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vtxdepthE_theta2_Sclock; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vrxdepthE_theta1_Sclock; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vrxdepthE_theta2_Sclock; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.

  std::vector< std::vector<double> > vtxdepth_dispersion1; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  std::vector< std::vector<double> > vtxdepth_dispersion2; // first element is station, third element is steps along the ray.  for now, this is just for ray 1.
  
  std::vector< std::vector<double> > vdotShats_tx;
  std::vector< std::vector<double> > vdotEhats_tx;
  std::vector< std::vector<double> > vdotDhats_tx;



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

  std::vector<std::vector<TVector3>> pol_r1;
  std::vector<std::vector<TVector3>> pol_r2;


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
  for (int istations=0;istations<6;istations++) {
    vraypos[istations].resize(3);
  }

  string sfile;
  if (WHICHPOL==0) 
    sfile="dave_data/day359_pol0.log";
  if (WHICHPOL==1) 
    sfile="dave_data/day359_pol1.log";

  //  ifstream myfile("dave_data/day359_pol0.log");
  ifstream myfile(sfile.c_str());
  //ifstream myfile("dave_data/day360_pol0.log");
//  ifstream myfile("dave_data/alldays_pol0.log");
  ifstream davea5file("dave_data/a5_amyformat.txt");

  string sn1file="birefringence_paper/n1.txt";
  string sn2file="birefringence_paper/n2.txt";
  string sn3file="birefringence_paper/n3.txt";

  ifstream n1file(sn1file.c_str());
  ifstream n2file(sn2file.c_str());
  ifstream n3file(sn3file.c_str());

  vector<double> vdepths_n1;
  vector<double> vdepths_n2;
  vector<double> vdepths_n3;

  vector<double> n1vec;
  vector<double> n2vec;
  vector<double> n3vec;
  std::vector<double>  vV; // the special angle V

  birefringenceFileRead(n1file, n2file, n3file, stemp, thisn, thisdepth, firstn1, firstn2, firstn3, 81)
 
  // smooth them.
   smoothVecs(n1vec,n2vec,n3vec,tmp)

  //What do I do with this?//
  TGraph *gn1=new TGraph(n1vec.size(),&vdepths_n1[0],&n1vec[0]);
  TGraph *gn2=new TGraph(n2vec.size(),&vdepths_n2[0],&n2vec[0]);
  TGraph *gn3=new TGraph(n3vec.size(),&vdepths_n3[0],&n3vec[0]);
  TGraph *g_V;  

 
  

  //v angle func
  Vangle(nvec_tmp,vdepths_n1, vdepths_n2,vdepths_n3,3)

  //find angle v, make function setVVector?, make graph

    g_V=new TGraph(vdepths_n1.size(),&vdepths_n1[0],&vV[0]);

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

    TGraph *gepsilon1_tx[6];
    TGraph *gepsilon2_tx[6];
    TGraph *gdiffepsilon_tx[6];

    TGraph *gepsilon1_rx[6];
    TGraph *gepsilon2_rx[6];
    TGraph *gdiffepsilon_rx[6];

    TGraph *gpolarization_Omega_rx[6];
    TGraph *gpolarization_Psi_rx[6];
    TGraph *gEpolarization_Omega_rx[6];
    TGraph *gEpolarization_Psi_rx[6];

    TGraph *gpolarization_reversedepth_Omega_rx[6];
    TGraph *gpolarization_reversedepth_Psi_rx[6];

    TGraph *gEpolarization_reversedepth_Omega_rx[6];
    TGraph *gEpolarization_reversedepth_Psi_rx[6];

    TGraph *g_receive_launch[6];
    TGraph *g_receive[6];
    TGraph *g_launch[6];
    TGraph *g_output6[6];
    TGraph *g_output7[6];
    TGraph *g_output8[6];

    voutput7.resize(6);
    voutput8.resize(6);
    videpth.resize(6);
    vreversedepth.resize(6);
    vangle_Shat_e1_khat.resize(6);
    vangle_Shat_e2_khat.resize(6);
    vangle_khat_0_khat_1_2.resize(6);
    vangle_khat_0_khat_2_2.resize(6);

  //  cout << "n1, n2, and n3 at 1000m depth is " << gn1->Eval(-1000.) << "\t" << gn2->Eval(-1000.) << "\t" << gn3->Eval(-1000.) << "\n";

  int NSHOTS=640;
 // const int NSHOTS=1536;{
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
  }
  //also make into a function!
  //while (!myfile.eof())
  readDepthFiles(myfile,Dave5afile,this_station,this_day,this_pol,this_depth,this_snrmax)
  void readDepthFiles(ifstream myfile,ifstream Dave5afile,int this_station, int this_day,int this_pol, double this_depth, double this_snrmax){ //take myfile and dave's data and read in values
    if (myfile.is_open())
      {
	for (int i=0;i<NSHOTS;i++) {
	  //cout << "i'm here. \n";
	  // this_depth is a negative number
	  myfile >> this_station >> this_day >> this_pol >> this_depth >> this_snrmax;
	  this_snrmax=this_snrmax*sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1])/sqrt((-1000.-station_depths[0])*(-1000.-station_depths[0])+horizontal_distances[0]*horizontal_distances[0]); // correct for 1/r
	//if (this_station==1)
	//	if (i>640)
	//	cout << "station, this_pol, depth, snrmax are " << this_station << "\t" << this_pol << "\t" << this_depth << "\t" << this_snrmax << "\n";
	//	this_snrmax=this_snrmax*exp(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1])/L_ATTEN);
	//cout << "this_depth is " << this_depth << "\n";
	
	  if (this_snrmax>0. && !(this_station==5 && this_pol==0)) {
	    vdepth_data[this_station-1].push_back(this_depth);
	    //videpth[this_station-1].push_back((double)(vdepth_data[this_station-1].size()-1));
	    vsnrmax[this_station-1].push_back(this_snrmax);
	    vtotal_distances[this_station-1].push_back(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1]));

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
	    vsnrmax_err[this_station-1].push_back(running_rms);
	    vdepth_data_err[this_station-1].push_back(0.);
	    vtotal_distances_err[this_station-1].push_back(0.);
	  } // if this_snrmax>0.
	}
      
	myfile.close();
      
      }
    NSHOTS=33;

    if (WHICHPOL==0) {
      if (davea5file.is_open())
	{
	  cout << "i'm reading dave's a5 file.\n";
	  for (int i=0;i<NSHOTS;i++) {
	    //cout << "i'm here. \n";
	    int this_station, this_day, this_pol;
	    double this_depth, this_snrmax;
	    string stemp;
	    int this_channel;
	    this_station=5; // station 5
	  
	  
	    davea5file >> stemp >> this_channel >> this_depth >> this_snrmax;
	  //if (this_station==1)
	  //if (i>640)
	  //	  cout << "station, depth, snrmax are " << this_station << "\t" << this_depth << "\t" << this_snrmax << "\n";
	  //this_snrmax=this_snrmax*exp(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1])/L_ATTEN);
	  //cout << "this_depth is " << this_depth << "\n";
	  
	    if (this_snrmax>0.) {
	      vdepth_data[this_station-1].push_back(this_depth);
	    //videpth[this_station-1].push_back((double)(vdepth_data[this_station-1].size()-1));
	      vsnrmax[this_station-1].push_back(this_snrmax);
	      vtotal_distances[this_station-1].push_back(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1]));
	    
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
	      vsnrmax_err[this_station-1].push_back(running_rms);
	      vdepth_data_err[this_station-1].push_back(0.);
	      vtotal_distances_err[this_station-1].push_back(0.);
	    } // if this_snrmax>0.
	  }
	
	  davea5file.close();
	
	}
    }
  }///end this function: should it be one or two???

  //  cout << "size of vdepth is " << vdepth_data.size() << "\n";
  //cout << "size of vdepth[4] is " << vdepth_data[4].size() << "\n";
  double minpulserdepth=-400.; //remove the arianna from the pulser variables
  double pulserstep=10.;
  int N_PULSER=200;
  for (int istations=0;istations<NSTATIONS;istations++) {
    for (int i=0;i<N_PULSER;i++) {
      vdepth[istations].push_back(minpulserdepth-(double)pulserstep*(double)i);
      videpth[istations].push_back((double)i);   
    }
  }
  for (int istations=minstation;istations<=maxstation;istations++) {
    for (int j=0;j<vdepth[istations].size();j++) {
      vreversedepth[istations].push_back(-1.*vdepth[istations][vdepth[istations].size()-1-j]);
    }
  }


  //  for (int i=0;i<NSTATIONS;i++) {
  for (int i=minstation;i<=maxstation;i++) {
    if (i!=5)
      gsnrmax[i]=new TGraph(vsnrmax[i].size(),&vdepth_data[i][0],&vsnrmax[i][0]);
    g_idepth[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&videpth[i][0]);
    //cout << "i, ispecial, size are " << i << "\t" << g_idepth[i]->Eval(depth_special) << "\t" << vdepth[i].size() << "\n";

  }
  // else cout << "Unable to open file.";
 

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


  vector <TVector3> rhat;
  vector <TVector3> rhat_receive;
  vector <TVector3> p_o; // this is at transmitter
  vector <TVector3> p_e; // at transmitter
  vector <TVector3> p_e_rx; // at rx
  vector <TVector3> p_o_rx; // at rx
  vector <TVector3> acrossp_o;
  vector <TVector3> acrossp_e;
  vector <TVector3> extraordinary;
  vector <TVector3> Pt;
  vector <TVector3> Pt_extraordinary;
  vector <TVector3> Pt_ordinary;
  vector <TVector3> Pr1;

  vector <TVector3> Pr2;
  vector <vector <double>> term1_Pr1;
  vector <vector <double>> term2_Pr1;
  vector <vector <double>> term3_Pr1;
  vector <vector <double>> term1_Pr2;
  vector <vector <double>> term2_Pr2;
  vector <vector <double>> term3_Pr2;

  vector <vector <double>> term1_Pr1_distances;
  vector <vector <double>> term2_Pr1_distances;
  vector <vector <double>> term3_Pr1_distances;
  vector <vector <double>> term1_Pr2_distances;
  vector <vector <double>> term2_Pr2_distances;
  vector <vector <double>> term3_Pr2_distances;

  TVector3 temp;

  //  vector <vector<double>> r;
  //vector <vector<double>> p_o;
  //vector <vector<double>> p_e;
  //vector <vector<double>> extraordinary;
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


  const int NSPECIALDEPTHS=2;
  //double specialdepths[NSPECIALDEPTHS]={-600.,-700.,-800.,-900.,-1000.,-1100.,-1200.,-1300.,-1400.,-1500.};

  //  double specialdepths[NSPECIALDEPTHS]={-850.};
  double specialdepths[NSPECIALDEPTHS]={-850.,-875.};
  //  TGraph *g_power[NSTATIONS][NSPECIALDEPTHS];
  vector< vector<TGraph*> > g_spectra;
  //vector< vector<TGraph*> > g_power;
  g_spectra.resize(NSTATIONS);
    string sfunc;
    //      sfunc="[0]*sqrt([1]+[2]-2*[3]*sin([4])*sin([4]))";
    sfunc="[0]*sqrt(([1]+[2])*([1]+[2])-2*[3]*sin([4])*sin([4]))";

    string sfunc_nointerference;
    sfunc_nointerference="[0]*sqrt(([1]+[2])*([1]+[2]))";
    //sfunc_nointerference="[0]*sqrt([1]*[1])+0.*[2]";

    string sfunc_distances;
    //    sfunc_distances="[0]*sqrt(([1]+[2]-2.*[3]*sin([4])*sin([4])))";
    sfunc_distances="[0]*sqrt(([1]+[2])*([1]+[2])-2*[3]*sin([4])*sin([4]))";
    //sfunc_distances="[0]*sqrt([1]+[2])";

    string sfunc_nointerference_distances;
    //    sfunc_nointerference_distances="[0]*sqrt([1]+[2])";
    sfunc_nointerference_distances="[0]*sqrt(([1]+[2])*([1]+[2]))";
    


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
    

    vector< vector<double> > vdistances_bigpic;
    vector< vector<double> > vpseudodepths_bigpic;
    vector< vector<double> > vmag_atten_beam_bigpic;
    vector< vector<double> > vmag_atten_beam_crosspol_bigpic;
    vector< vector<double> > vmag_atten_beam_crosspol_nointerferencefunc_bigpic;
    vector< vector<double> > vmag_atten_beam_crosspol_func_bigpic;
    vector< vector<double> > vtimediff_bigpic;
    //    vector< vector<double> > vnsolutions_bigpic;
    vector< vector<double> > vspectrum_bigpic;

    vmag_atten_beam_bigpic.resize(NSTATIONS);
    vmag_atten_beam_crosspol_bigpic.resize(NSTATIONS);
    vmag_atten_beam_crosspol_nointerferencefunc_bigpic.resize(NSTATIONS);
    vmag_atten_beam_crosspol_func_bigpic.resize(NSTATIONS);
    vtimediff_bigpic.resize(NSTATIONS);
    vspectrum_bigpic.resize(NSTATIONS);

    const int NDISTANCES_BIGPIC=5000;
    const double STEP=1.;


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
    vector< vector<double> > vnotflipped; // 1 if not flipped, -1 if flipped


    vector< vector<double> > vlengths;
    vector< vector<double> > vattenlengths;

    vector< vector< vector<double> > > vspectra;
    vector< vector< vector<double> > > vattens;
    
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
    //vnsolutions_bigpic.resize(NSTATIONS);
    
    
    int jmin[NSTATIONS];
    //    for (int i=0;i<NSTATIONS;i++) {
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
   
    if (i!=5) {
      vtotal_distances[i].resize(vdepth_data[i].size());
      vtotal_distances_err[i].resize(vdepth_data[i].size());
    }
    double mindistance=horizontal_distances[i];
    jmin[i]=(int)(mindistance/STEP);
    for (int j=jmin[i]+1;j<NDISTANCES_BIGPIC;j++) {
      double thisdistance=STEP*(double)j;
      vdistances_bigpic[i].push_back(thisdistance);      
      vpseudodepths_bigpic[i].push_back(station_depths[i]-1.*sqrt(thisdistance*thisdistance-mindistance*mindistance));
      //cout << "thisdistance, mindistance, vpseudodepths_bigpic are " << thisdistance << "\t" << mindistance << "\t" << vpseudodepths_bigpic[i][vpseudodepths_bigpic[i].size()-1] << "\n";
    }



    
    f1[i]=new TF1("f1",sfunc.c_str(),-2000.,-100.);
  
    f1_nointerference[i]=new TF1("f1_nointerference",sfunc_nointerference.c_str(),-2000.,-100.);

    //cout << "length for station " << i << " is " << vdepth[i].size() << "\n";
  
    for (int idepth=0;idepth<vdepth[i].size();idepth++) {

      //      cout << "idepth is " << idepth << "\n";


      double posstation[2];
      posstation[0]=sqrt(pow(station_coords[i][0]-pulser_coords[0],2)+pow(station_coords[i][1]-pulser_coords[1],2));
      posstation[1]=station_depths[i];
      double pospulser[2];
      pospulser[0]=0.;
      pospulser[1]=vdepth[i][idepth];

      TVector3 pospulser3D;
      pospulser3D.SetX(pulser_coords[0]);
      pospulser3D.SetY(pulser_coords[1]);
      pospulser3D.SetZ(vdepth[i][idepth]);


      pulsertostation[i][0]=(station_coords[i][0]-pulser_coords[0]); // this will be a negative number
      pulsertostation[i][1]=(station_coords[i][1]-pulser_coords[1]); 
      pulsertostation[i][2]=station_depths[i]-vdepth[i][idepth];


      double atten_special;
      double pospulser_special[2];
      pospulser_special[0]=0.;
      pospulser_special[1]=vdepth[0][0];
      double posstation_special[2];
      posstation_special[0]=sqrt(pow(station_coords[0][0]-pulser_coords[0],2)+pow(station_coords[0][1]-pulser_coords[1],2));
      posstation_special[1]=station_depths[0];
  
  
      // cout << "posstation is " << posstation[0] << "\t" << posstation[1] << "\n";
      // cout << "pospulser is " << pospulser[0] << "\t" << pospulser[1] << "\n";

      // cout << "about to find solutions.\n";

      double x0=0;
      double z0=pospulser[1];
      double x1=posstation[0];
      double z1=posstation[1];

      //      cout << "z0, x1, z1 are " << z0 << "\t" << x1 << "\t" << z1 << "\n";

      


      double *getresults=IceRayTracing(x0,z0,x1,z1);

      double lvalue;
      vector<double> res; 	//res are y coordinates
      vector<double> zs;	//zs are z coordinates
      double launch_angle;
      double receive_angle;
      double *paramsd;
      double *paramsra;
      double *paramsre;

      voutput6[i].push_back(getresults[6]);
      voutput7[i].push_back(getresults[7]);
      voutput8[i].push_back(getresults[8]);

     if (getresults[6]!=-1000) { 
	paramsd=GetDirectRayPar(z0,x1,z1);
	launch_angle=paramsd[1]/DEGRAD;
	receive_angle=paramsd[0]/DEGRAD;
	GetFullDirectRayPath(z0, x1, z1, paramsd[3], res, zs);
      }
      else if (getresults[8]!=1000) {
	//	paramsd=GetDirectRayPar(z0,x1,z1);
	paramsre=GetReflectedRayPar(z0, x1 ,z1);
	double LangR=paramsre[1];
	double RangR=paramsre[0];
	paramsra=GetRefractedRayPar(z0, x1 ,z1,LangR,RangR);
	launch_angle=paramsra[1]/DEGRAD;
	//	receive_angle=(180.-paramsra[0])/DEGRAD;
	receive_angle=paramsra[0]/DEGRAD;
	GetFullRefractedRayPath(z0, x1, z1, paramsra[7], paramsra[3], res, zs);
	// need to get fullrefractedraypath
      }
      else if (getresults[7]!=1000) {
	//      	paramsd=GetDirectRayPar(z0,x1,z1);
	paramsre=GetReflectedRayPar(z0, x1 ,z1);
	double LangR=paramsre[1];
	double RangR=paramsre[0];
	launch_angle=paramsre[1]/DEGRAD;
	receive_angle=paramsre[0]/DEGRAD;
	GetFullReflectedRayPath(z0, x1, z1, LangR, res, zs);
      }

      //      vector< vector<double> > solutions = find_solutions(pospulser,posstation,
      //						  N_ICE,DELTA_N,Z_0);

    //vector< vector<double> > solutions_special = find_solutions(pospulser_special,posstation_special,
      //						  N_ICE,DELTA_N,Z_0);



       
    //  cout << "size of solutions is " << solutions.size() << "\n";

      //      cout << "i=" << i << " about to get attenuation.\n";


      double atten=1.;
      //      if (idepth==imostshallowdepth[i]) {
      //if (idepth==igreatestdepth[i]) {
      //      if (idepth==(int)(g_idepth[i]->Eval(depth_special))) {
      //      if (idepth==vdepth[i].size()-1) {
	//if (idepth==ispecial) {
	//cout << "this should be the special one.  depth is " << vdepth[i][idepth] << "\n";
	vattens[i][idepth].clear();
	for (int ifreq=0;ifreq<NFREQ;ifreq++) {
	  vattens[i][idepth].push_back(1.);
	}

	//} // if it's the special depth
      // inside loop over idepth

      double sumlength=0.;
      double sumphase=0.;
      //if (solutions.size()>0) {
      if (getresults[6]!=-1000 || getresults[8]!=-1000 || getresults[7]!=0) {
	//      cout << "attenuation length is " << get_attenuation_length(zs[0],freq/1.E9,1) << "\n";

	//	atten=get_attenuation_along_path(pospulser, posstation, solutions[0][1],
	//				 freq/1000., N_ICE, DELTA_N, Z_0, 1);
	//cout << "atten is " << atten << "\n";
	//double path_length=get_path_length(pospulser, posstation, solutions[0][1], N_ICE, DELTA_N, Z_0);
	//cout << "path_length is " << path_length << "\n";
	


	//	get_path(N_ICE, DELTA_N, Z_0, pospulser, posstation, solutions[0][1], res, zs, 200);
	


	  // the "y" direction is the horizontal direction from the pulser to the station.
	TVector3 yhat(station_coords[i][0]-pulser_coords[0],
		      station_coords[i][1]-pulser_coords[1],
		      0.); // yhat points from pulser to station
	if (yhat.Mag()<HOWSMALLISTOOSMALL)
	  cout << "yhat mag is " << yhat.Mag() << "\n";
	yhat.SetMag(1.);
	//cout << "A" << i+1 << ": res size is " << res.size() << "\n";
	double angle_yhat=atan2(yhat[1],yhat[0]);


	vector<double> nvec_thisstep;
	nvec_thisstep.resize(3);
	
	nvec_thisstep[0]=gn1->Eval(zs[0]);
	nvec_thisstep[1]=gn2->Eval(zs[0]);
	nvec_thisstep[2]=gn3->Eval(zs[0]);

	TVector3 rhat_thisstep;

	rhat_thisstep[0]=-1.*(res[UZAIRSTEP]-res[0])*yhat[0];
	rhat_thisstep[1]=-1.*(res[UZAIRSTEP]-res[0])*yhat[1];
	rhat_thisstep[2]=-1.*(zs[UZAIRSTEP]-zs[0]);

	if (rhat_thisstep.Mag()<1.E-8)
	  cout << "before calling getDeltaN at place 1, rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
	double deltan_alongpath=getDeltaN(BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1,p_e2);

	
	


	if (p_e2.Mag()<HOWSMALLISTOOSMALL) 
	  cout << "1, p_e2 is " << p_e2.Mag();
	
	
	//	 cout << "p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
	//cout << "p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";
      
    
	double deltantimeslength_alongpath=0.;

	TVector3 p_e1_previous=p_e1;
	TVector3 p_e2_previous=p_e2;
	double notflipped_previous=1.;
	double notflipped_atend=1.;
	double theta_e1_start=0.;
	double notflipped=1.;


	for (int istep=UZAIRSTEP;istep<res.size();istep+=UZAIRSTEP) {
	  
	  nvec_thisstep.resize(3);


	  nvec_thisstep[0]=gn1->Eval(zs[istep]);
	  nvec_thisstep[1]=gn2->Eval(zs[istep]);
	  nvec_thisstep[2]=gn3->Eval(zs[istep]);

	  if (istep>0) {

	    rhat_thisstep[0]=-1.*(res[istep]-res[istep-UZAIRSTEP])*yhat[0];
	    rhat_thisstep[1]=-1.*(res[istep]-res[istep-UZAIRSTEP])*yhat[1];
	    rhat_thisstep[2]=-1.*(zs[istep]-zs[istep-UZAIRSTEP]);

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
	      //	      cout << "zs, res this are " << zs[istep] << "\t" << res[istep] << "\n";
	      //cout << "zs, res previous are " << zs[istep-1] << "\t" << res[istep-1] << "\n";
	    }

	    double length=rhat_thisstep.Mag();
    
	    if (rhat_thisstep.Mag()<HOWSMALLISTOOSMALL)
	      cout << "rhat_thisstep mag is " << rhat_thisstep.Mag() << "\n";

	    rhat_thisstep.SetMag(1.);

	    //	    double atten_length=get_attenuation_length(zs[istep],freq/1.E9,1);
	    double atten_length=GetIceAttenuationLength(zs[istep], freq/1.E9);

	    //	    cout << "zs, length, atten_length, atten are " << zs[istep] << "\t" << length << "\t" << atten_length << "\t" << atten << "\n";

	    atten*=exp(-1.*length/atten_length);
	    
	    if (rhat_thisstep.Mag()<1.E-8)
	      cout << "before calling getDeltaN at place 2, rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
	    deltan_alongpath=getDeltaN(BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1,p_e2);
	    //cout << "deltan_alongpath 1 is " << deltan_alongpath << "\n";
	    if (p_e2.Mag()<HOWSMALLISTOOSMALL) 
	      cout << "2, p_e2 is " << p_e2.Mag() << "\n";


	    // if (notflipped*notflipped_previous<0.) {
	    //   p_e1=-1.*p_e1;
	    //   p_e2=-1.*p_e2;
	    // }

	    // if (i==5 || i==0) {

	    //   cout << "A" << i+1 << " depth " << vdepth[i][idepth] << "\n";
	    //   cout << "notflipped is " << notflipped << "\n";
	    //   cout << "p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
	    //   cout << "p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";
	    //   cout << "angle between p_e1_previous and p_e1 is " << acos(p_e1_previous.Dot(p_e1))*DEGRAD << "\n";
	    //   cout << "angle between p_e2_previous and p_e2 is " << acos(p_e2_previous.Dot(p_e2))*DEGRAD << "\n";
	    //   cout << "angle between p_e1_previous and p_e2 is " << acos(p_e1_previous.Dot(p_e2))*DEGRAD << "\n";
	    //   cout << "angle between p_e2_previous and p_e1 is " << acos(p_e2_previous.Dot(p_e1))*DEGRAD << "\n";
	      
	    // }
	    
	    p_e1_previous=p_e1;
	    p_e2_previous=p_e2;

	    TVector3 epsilon_thisstep;
	  

	    epsilon_thisstep[0]=nvec_thisstep[0]*nvec_thisstep[0];
	    epsilon_thisstep[1]=nvec_thisstep[1]*nvec_thisstep[1];
	    epsilon_thisstep[2]=nvec_thisstep[2]*nvec_thisstep[2];
	    
	    // electric fields associated with ray1 and ray2.
	    // these seem to be unnormalized.
	    TVector3 E_e1=rotateD(epsilon_thisstep,angle_iceflow,p_e1);
	    TVector3 E_e2=rotateD(epsilon_thisstep,angle_iceflow,p_e2);

	    if (abs((double)(istep-(int)res.size()))<=UZAIRSTEP) {
	      // this is the transmitter.
	      // lower step numbers are more shallow
	    if (i==5 && idepth==g_idepth[i]->Eval(-400.)) 	      
		cout << "launch angle is " << rhat_thisstep.Theta()*DEGRAD << "\n";

	      double theta_e1,theta_e2;
	      double thetaE_e1,thetaE_e2;
	      double theta_e1_Sclock,theta_e2_Sclock;
	      double thetaE_e1_Sclock,thetaE_e2_Sclock;
	      
	      // here i want to plot where p_e1 and p_2 are on the clock, with 12 o'clock being the in the plane of rhat at launch and the z axis.
	      
	      TVector3 Shat_e1,Shat_e2;
	      
	      double E_e1_thetacomponent,E_e2_thetacomponent;
	      double E_e1_phicomponent,E_e2_phicomponent;

	      getManyAnglesontheClock(BIAXIAL,CROSSPOLANGLE_TX,
				      rhat_thisstep,
				      p_e1,p_e2,E_e1,E_e2,
				      theta_e1,theta_e2,thetaE_e1,thetaE_e2,
				      theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock,
				      Shat_e1,Shat_e2,
				      E_e1_thetacomponent,E_e2_thetacomponent,
				      E_e1_phicomponent,E_e2_phicomponent);

	      vangle_Shat_e1_khat[i].push_back(acos(Shat_e1.Dot(rhat_thisstep)/Shat_e1.Mag()/rhat_thisstep.Mag())*DEGRAD);
	      vangle_Shat_e2_khat[i].push_back(acos(Shat_e2.Dot(rhat_thisstep)/Shat_e2.Mag()/rhat_thisstep.Mag())*DEGRAD);

	      //	      cout << "theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";
	      //cout << "diff is " << theta_e1_Sclock-theta_e2_Sclock << "\n";

	  
	      //    E_e1.SetMag(fabs(E_e1_component));
	      //E_e2.SetMag(fabs(E_e2_component));
	      //if (E_e1_component<0.)
	      //E_e1=-1.*E_e1;
	      //if (E_e2_component<0.)
	      //E_e2=-1.*E_e2;

	      if (i==5 && BIAXIAL==1) { // for arianna, in the biaxial crystal the 2 rays trade places mid-flight
		switchThem(thetaE_e1_Sclock,thetaE_e2_Sclock);
		switchThem(theta_e1_Sclock,theta_e2_Sclock);
		switchThem(theta_e1,theta_e2);
		switchThem(thetaE_e1,thetaE_e2);
		//	thetaE_e2_Sclock=thetaE_e2_Sclock-PI;
	      }

	      if (i==0 && idepth==g_idepth[i]->Eval(-1000.) ||
		  i==5 && idepth==g_idepth[i]->Eval(-400.)) {

		cout << "At Tx:\n";
		cout << "A" << i+1 << ", depth is " << vdepth[i][idepth] << "\n";
		// yhat is the intersection of the horizontal plane and the plane of the ray (remember ray is always tangent to k-hat)
		cout << "angle of yhat is " << angle_yhat << "\n";
		//cout << "alpha is " << alpha[i] << "\n";

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
		//E_e1_temp.SetMag(1.);
		
		TVector3 E_e2_temp=E_e2;
		//E_e2_temp.SetMag(1.);
		
		TVector3 Shat_e1_temp=Shat_e1;
		TVector3 Shat_e2_temp=Shat_e2;
		TVector3 rhat_thisstep_temp=rhat_thisstep;

		TVector3 E_total_temp=E_e1_temp+E_e2_temp;
		 //		 E_total_temp.SetMag(1.);

		// rotate them so the line of sight from pulser to station
		// is in the plane of the page.

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


		 //cout << "components are " << E_e1_component << "\t" << E_e2_component << "\n";
		cout << "These are rotated so the line of sight from pulser to station is in the plane of the page.\n";
		cout << "E_e1_temp is " << E_e1_temp[0] << "\t" << E_e1_temp[1] << "\t" << E_e1_temp[2] << "\n";
		cout << "mag of E_e1_temp is " << E_e1_temp.Mag() << "\n";
		cout << "E_e2_temp is " << E_e2_temp[0] << "\t" << E_e2_temp[1] << "\t" << E_e2_temp[2] << "\n";
		cout << "mag of E_e2_temp is " << E_e2_temp.Mag() << "\n";
		// if S-hat were in the same direction as k-hat, then E_total_temp should not have a y-component.  but they are slightly different
		// so E_total_temp does have a slight y component.
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

	      vdotShats_tx[i].push_back(Shat_e1.Dot(Shat_e2)/Shat_e1.Mag()/Shat_e2.Mag());
	      vdotEhats_tx[i].push_back(E_e1.Dot(E_e2)/E_e1.Mag()/E_e2.Mag());
	      vdotDhats_tx[i].push_back(p_e1.Dot(p_e2)/p_e1.Mag()/p_e2.Mag());

	      if (i==0 && idepth==g_idepth[i]->Eval(-1000.) ||
		  i==5 && idepth==g_idepth[i]->Eval(-400.)) {
		cout << "Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
		cout << "Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";
		cout << "theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";
	      }

	      double beam_tx_e1=sin(Shat_e1.Theta());
	      double beam_tx_e2=sin(Shat_e2.Theta());

	      vtxdepth_beam1[i].push_back(beam_tx_e1);
	      vtxdepth_beam2[i].push_back(beam_tx_e2);

	      vtxdepth_theta1[i].push_back(theta_e1*DEGRAD);
	      vtxdepth_theta2[i].push_back(theta_e2*DEGRAD);
	      
	      vtxdepth_theta1_Sclock[i].push_back(theta_e1_Sclock*DEGRAD);
	      vtxdepth_theta2_Sclock[i].push_back(theta_e2_Sclock*DEGRAD);
	      
	      vtxdepthE_theta1[i].push_back(thetaE_e1*DEGRAD);
	      vtxdepthE_theta2[i].push_back(thetaE_e2*DEGRAD);

	      // if (BIAXIAL==-1) {
	      // 	vtxdepthE_theta1_Sclock[i].push_back(vrxdepthE_theta1_Sclock[i][idepth]);
	      // 	vtxdepthE_theta2_Sclock[i].push_back(vrxdepthE_theta2_Sclock[i][idepth]);

	      // } 
	      //	      else {
	      vtxdepthE_theta1_Sclock[i].push_back(thetaE_e1_Sclock*DEGRAD);
	      vtxdepthE_theta2_Sclock[i].push_back(thetaE_e2_Sclock*DEGRAD);
		//}
	      vtxdepth_dispersion1[i].push_back(acos(E_e1.Dot(p_e1)/E_e1.Mag()/p_e1.Mag())*DEGRAD);
	      vtxdepth_dispersion2[i].push_back(acos(E_e2.Dot(p_e2)/E_e2.Mag()/p_e2.Mag())*DEGRAD);
	      
	      //if (i==5) {
	      //cout << "depth is " << vdepth[i][idepth] << "\n";
	      //cout << "before switching, thetas are " << thetaE_e1_Sclock << "\t" << thetaE_e2_Sclock << "\n";
	      //}

	      if (i==5 && BIAXIAL==1) {
		//double tempangle=thetaE_e1_Sclock;
		//thetaE_e1_Sclock=thetaE_e2_Sclock;
		//thetaE_e2_Sclock=tempangle;
		// //		thetaE_e1_Sclock-=PI/2.;

		//		thetaE_e2_Sclock-=PI;

		//	if (thetaE_e1_Sclock<0.)
		//thetaE_e1_Sclock+=PI;
		//if (thetaE_e2_Sclock<0.)
		//thetaE_e2_Sclock+=PI;
		

	      }
	      else {
		//		thetaE_e1_Sclock-=PI/2.;
	      }
	      //	      if (i==5) {
	      //cout << "after switching, theta_e1's at the tx are " << theta_e1*DEGRAD << "\t" << theta_e2*DEGRAD << "\n";
	      //cout << "after switching, thetaE_e1_Sclock's at the tx are " << thetaE_e1_Sclock*DEGRAD << "\t" << thetaE_e2_Sclock*DEGRAD << "\n";
	      //}

	      // instead of doing the switcheroo above, try just
	      // having a conditional on the next two lines.
	      // this will be used to get our epsilon
	      double epsilon1_tx=0.;
	      double epsilon2_tx=0.;
	      thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,
			       epsilon1_tx,epsilon2_tx);
	
	      if (epsilon2_tx>PI/2.)
		epsilon2_tx-=PI;

	      //	      double epsilon1_tx=thetaE_e1_Sclock-PI/2.;
	      //double epsilon2_tx=thetaE_e2_Sclock;
	   
	      //	      if (i==5)
	      //cout << "i, epsilons are " << i << "\t" << epsilon1_tx*DEGRAD << "\t" << epsilon2_tx*DEGRAD << "\n";
	   

	      //	      if (i==5) {
	      //epsilon1_tx=epsilon1_tx+PI;
	      //epsilon2_tx=epsilon2_tx+PI;		
		//cout << "i, epsilons at the tx are " << i << "\t" << epsilon1_tx << "\t" << epsilon2_tx << "radians. \n";
	      //}

	      // if (i==0 || i==5) {
	      // cout << "i, depth are " << i << "\t" << vdepth[i][idepth] << "\n";
	      // cout << "thetas at tx are " << thetaE_e1_Sclock << "\t" << thetaE_e2_Sclock << "\n";
	      // cout << "epsilons at tx are " << epsilon1_tx << "\t" << epsilon2_tx << "\n";
	      // }

	      vepsilon1_tx[i].push_back(epsilon1_tx*DEGRAD);
	      vepsilon2_tx[i].push_back(epsilon2_tx*DEGRAD);




	      // if (i==5) {
	      // 	cout << "at the tx, epsilon1, epsilon2 (in degrees), sizes, idepth are " << epsilon1_tx*DEGRAD << "\t" << epsilon2_tx*DEGRAD << "\t" << vepsilon1_tx[i].size() << "\t" << vepsilon2_tx[i].size() << "\t" << idepth << "\n";
	      // }

	      vdiffepsilon_tx[i].push_back((epsilon2_tx-epsilon1_tx)*DEGRAD);
	 	 
	      //	      cout << "attenuation is " << atten << "\n";
	      vrxdepth_atten[i].push_back(VOLTAGENORM*atten);
	      vrxdepth_atten_beam[i].push_back(VOLTAGENORM*vrxdepth_beam1[i][idepth]*vtxdepth_beam1[i][idepth]*atten);
	      vrxdepth_atten_power[i].push_back(VOLTAGENORM*VOLTAGENORM*atten*atten);
	      vrxdepth_atten_beam_power[i].push_back(VOLTAGENORM*VOLTAGENORM*vrxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]*vtxdepth_beam1[i][idepth]*vtxdepth_beam1[i][idepth]*atten*atten);
	      notflipped_atend=notflipped;
	    } // if step at the transmitter
	    if (istep==UZAIRSTEP) {
	      // this is the receiver.
	      if (i==5 && idepth==g_idepth[i]->Eval(-400.)) { 
		
	      }
	      double theta_e1,theta_e2;
	      double thetaE_e1,thetaE_e2;
	      double theta_e1_Sclock,theta_e2_Sclock;
	      double thetaE_e1_Sclock,thetaE_e2_Sclock;
	      
	      // here i want to plot where p_e1 and p_2 are on the clock, with 12 o'clock being the in the plane of rhat at launch and the z axis.
	      
	      TVector3 Shat_e1,Shat_e2;


	      double E_e1_thetacomponent,E_e2_thetacomponent;
	      double E_e1_phicomponent,E_e2_phicomponent;

	      getManyAnglesontheClock(BIAXIAL,CROSSPOLANGLE_RX,
				      rhat_thisstep,
				      p_e1,p_e2,E_e1,E_e2,
				      theta_e1,theta_e2,thetaE_e1,thetaE_e2,
				      theta_e1_Sclock,theta_e2_Sclock,thetaE_e1_Sclock,thetaE_e2_Sclock,
				      Shat_e1,Shat_e2,
				      E_e1_thetacomponent,E_e2_thetacomponent,
				      E_e1_phicomponent,E_e2_phicomponent);
	      
	      //	      double E_e1_component=sqrt(E_e1_thetacomponent*E_e1_thetacomponent+E_e1_phicomponent*E_e1_phicomponent);
	      //double E_e2_component=sqrt(E_e2_thetacomponent*E_e2_thetacomponent+E_e2_phicomponent*E_e2_phicomponent);
	      //E_e1.SetMag(fabs(E_e1_component));
	      //E_e2.SetMag(fabs(E_e2_component));
	      //if (E_e1_component<0.)
	      //E_e1=-1.*E_e1;
	      //if (E_e2_component<0.)
	      //E_e2=-1.*E_e2;

	      // uzair's ray tracing starts at the receiver, so this is the theta of ray 1 at the start.
	      theta_e1_start=theta_e1;
	      

	      if (i==0 && idepth==g_idepth[i]->Eval(-1000.) ||
		  i==5 && idepth==g_idepth[i]->Eval(-400.)) {
		cout << "At Rx:\n";
		cout << "A" << i+1 << ", depth is " << vdepth[i][idepth] << "\n";
		cout << "thetas on the Sclock are " << thetaE_e1_Sclock*DEGRAD << "\t" << thetaE_e2_Sclock*DEGRAD << "\n";
		cout << "depth is " << vdepth[i][idepth] << "\n";
		// }
		cout << "theta of rhat_thisstep is " << rhat_thisstep.Theta()*DEGRAD << "\n";
		cout << "rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
		cout << "This points from pulser to station: " << station_coords[i][0]-pulser_coords[0] << "\t" << station_coords[i][1]-pulser_coords[1] << "\t" << station_depths[i]-vdepth[i][idepth] << "\n";

		cout << "Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
		cout << "Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";
		cout << "theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";
		TVector3 E_e1_temp=E_e1;
		//E_e1_temp.SetMag(1.);
		
		TVector3 E_e2_temp=E_e2;
		//E_e2_temp.SetMag(1.);

		TVector3 Shat_e1_temp=Shat_e1;
		TVector3 Shat_e2_temp=Shat_e2;
		TVector3 rhat_thisstep_temp=rhat_thisstep;
		
		TVector3 E_total_temp=E_e1_temp+E_e2_temp;
		//		 E_total_temp.SetMag(1.);
		
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


		//		cout << "components are " << E_e1_component << "\t" << E_e2_component << "\n";
		
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


	      if (i==5 && BIAXIAL==1) {
		//cout << "depth, theta_e1_start is " << vdepth[i][idepth] << "\t" << theta_e1_start << "\n";

	      //cout << "A" << i+1 << ": thetas on the Sclock are " << thetaE_e1_Sclock*DEGRAD << "\t" << thetaE_e2_Sclock*DEGRAD << "\n";
		//double tempangle=thetaE_e1_Sclock;
		//thetaE_e1_Sclock=thetaE_e2_Sclock;
		//thetaE_e2_Sclock=tempangle-PI;


	      }


	      // if (i==5) {
	      //  	cout << "depth is " << vdepth[i][idepth] << "\n";
	      // 	cout << "at the rx, after switching, theta_e1's are " << theta_e1*DEGRAD << "\t" << theta_e2*DEGRAD << "\n";
	      //  	cout << "at the rx, after switching, thetaE_e1_Sclock's are " << thetaE_e1_Sclock*DEGRAD << "\t" << thetaE_e2_Sclock*DEGRAD << "\n";
		
	      // }


	      vrxdepth_theta1[i].push_back(theta_e1*DEGRAD);
	      vrxdepth_theta2[i].push_back(theta_e2*DEGRAD);

	      vrxdepthE_theta1[i].push_back(theta_e1*DEGRAD);
	      vrxdepthE_theta2[i].push_back(theta_e2*DEGRAD);
	    
	      vrxdepth_theta1_Sclock[i].push_back(theta_e1_Sclock*DEGRAD);
	      vrxdepth_theta2_Sclock[i].push_back(theta_e2_Sclock*DEGRAD);



	      vrxdepthE_theta1_Sclock[i].push_back(thetaE_e1_Sclock*DEGRAD);
	      vrxdepthE_theta2_Sclock[i].push_back(thetaE_e2_Sclock*DEGRAD);


	      //	      if (i==5) {
	      //thetaE_e2_Sclock+=PI/2.;
	      //}
	      //else {
	      //thetaE_e1_Sclock-=PI/2.;
	      //}
	      double epsilon1_rx=0.;
	      double epsilon2_rx=0.;
				   
	      thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,
			       epsilon1_rx,epsilon2_rx);
	      //	      double epsilon1_rx=thetaE_e1_Sclock-PI/2.;
	      //double epsilon2_rx=thetaE_e2_Sclock;

	      // if (i==0 || i==5) {
	      // cout << "i, depth is " << i << "\t" << vdepth[i][idepth] << "\n"; 
	      // cout << "thetaE_e1_Sclock's at the rx are " << thetaE_e1_Sclock*DEGRAD << "\t" << thetaE_e2_Sclock*DEGRAD << "\n";
	      // cout << "epsilons at the rx are " << epsilon1_rx*DEGRAD << "\t" << epsilon2_rx*DEGRAD << " degrees. \n";
	      // }


	      //	      if (i==1 || i==5) {
	      //cout << "epsilons at rx are " << epsilon1_rx << "\t" << epsilon2_rx << "\n";
	      //	      }

	      vepsilon1_rx[i].push_back(epsilon1_rx*DEGRAD);
	      vepsilon2_rx[i].push_back(epsilon2_rx*DEGRAD);
	      vdiffepsilon_rx[i].push_back((epsilon2_rx-epsilon1_rx)*DEGRAD);

	      double beam_rx_e1=sin(Shat_e1.Theta());
	      double beam_rx_e2=sin(Shat_e2.Theta());

	      vrxdepth_beam1[i].push_back(beam_rx_e1);
	      vrxdepth_beam2[i].push_back(beam_rx_e2);

	      //cout << "A" << i+1 << ": entered this if statement.\n";
	      TVector3 plusz(0.,0.,1.);
	      Pr2[i]=Shat_e1.Cross(plusz);
	      if (Pr2[i].Mag()<HOWSMALLISTOOSMALL)
		  cout << "Pr2[i] is " << Pr2[i].Mag() << "\n";
	      Pr2[i].SetMag(1.);
	      Pr1[i]=Pr2[i].Cross(Shat_e1);
	      if (Pr1[i].Mag()<HOWSMALLISTOOSMALL)
		  cout << "Pr1[i] is " << Pr1[i].Mag() << "\n";
	      Pr1[i].SetMag(1.);

	    } // if step at the receiver


	  
	
	
	    //	    if (idepth==igreatestdepth[i])
	    //   if (idepth==imostshallowdepth[i])
	    //if (idepth==(int)(g_idepth[i]->Eval(depth_special)))
      //	    if (idepth==vdepth[i].size()-1)
	    //if (idepth==ispecial)
	      for (int ifreq=0;ifreq<NFREQ;ifreq++) {
		//		double this_atten_length=get_attenuation_length(zs[istep],vfreqs[ifreq]/1.E9,1);

		double this_atten_length=GetIceAttenuationLength(zs[istep], vfreqs[ifreq]/1.E9);

		vattens[i][idepth][ifreq] = vattens[i][idepth][ifreq]*exp(-1.*length/this_atten_length)*exp(-1.*length/this_atten_length);
	      }
	    // cout << "about to call the offending getDeltaN.\n";
	    // cout << "rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
	 
	  
	
	    //	    cout << "coarse: deltan_alongpath, one period is " << deltan_alongpath << "\t" << 1./(deltan_alongpath*PI/TMath::C()*freq) << "\n";
	    
	    //	    cout << "A" << i+1 << ": istep, res.size(), abs are " << istep << "\t" << res.size() << "\t" << abs(double(istep-(int)res.size())) << "\n";


	    
	  
	    
	    //if (idepth==igreatestdepth[i]) {
	      //if (idepth==imostshallowdepth[i]) {
	    



	      double theta_e1,theta_e2;
	      double thetaE_e1,thetaE_e2;
	      double theta_e1_Sclock,theta_e2_Sclock;
	      double thetaE_e1_Sclock,thetaE_e2_Sclock;

	      // here i want to plot where p_e1 and p_2 are on the clock, with 12 o'clock being the in the plane of rhat at launch and the z axis.
	      //cout<< "position 1.45\n";
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
	      //cout<< "position 1.445\n";
	      //	      double E_e1_component=sqrt(E_e1_thetacomponent*E_e1_thetacomponent+E_e1_phicomponent*E_e1_phicomponent);
	      //double E_e2_component=sqrt(E_e2_thetacomponent*E_e2_thetacomponent+E_e2_phicomponent*E_e2_phicomponent);		
		     
	      //E_e1.SetMag(fabs(E_e1_component));
	      //E_e2.SetMag(fabs(E_e2_component));
	      //if (E_e1_component<0.)
	      //	E_e1=-1.*E_e1;
	      //if (E_e2_component<0.)
	      //	E_e2=-1.*E_e2;

	      //	      if (idepth==g_idepth[i]->Eval(-1000.)) {
	      //cout << "depth is " << vdepth[i][idepth] << "\n";
	      //}

	      //	      if (i==5) {
		//cout << "istep is " << istep << "\n";
		//cout << "step depth is " << zs[istep] << "\n";
		notflipped=Flipped(theta_e1, theta_e1_start);
		//}
       
		deltantimeslength_alongpath+=deltan_alongpath*length*notflipped;
		notflipped_previous=notflipped;

		//double term=deltan_alongpath*PI*length/TMath::C()*freq*notflipped;
		//if (i==5)
		//term=-1.*term;
		
		
	      //sumphase+=term;
		
		
		sumlength+=length;

		//	    cout<< "position 1.4\n";

	      if (idepth==(int)(g_idepth[i]->Eval(depth_special))) {
		//these are filled for every step along the path.
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
	   
	      //	      if (i==5 || i==0) {
		//		cout << "A" << i+1 << ", istep " << istep << ", depth " << zs[istep] << " lengths " << vlengths[i][vlengths.size()-1] << ": notflipped is " << notflipped << "\n";
	      //if (Flipped(p_e1_previous, p_e2_previous, p_e1, p_e2)<0. || (zs[istep]<-100. && zs[istep]>-200.)) {
	      //cout << "p_e1_previous is " << p_e1_previous[0] << "\t" << p_e1_previous[1] << "\t" << p_e1_previous[2] << "\n";
	      //cout << "p_e2_previous is " << p_e2_previous[0] << "\t" << p_e2_previous[1] << "\t" << p_e2_previous[2] << "\n";
	      //cout << "p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
		//   cout << "p_e1 theta is " << p_e1.Theta()*DEGRAD << "\n";
	      
		//cout << "p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";
		//	      cout << "p_e2 theta is " << p_e2.Theta()*DEGRAD << "\n";
		//}

	      //	    }
		//	      } // if station 1 or arianna



	    } // if it's the right depth

	  } // if number of steps >0
	  else {
	    rhat_thisstep=rhat[i];
	  }
	} // end loop over steps along path

	sumphase=deltantimeslength_alongpath*PI/TMath::C()*freq;

	  //vV1_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vepsilon1_tx[i][idepth]/DEGRAD)*cos(vepsilon1_rx[i][idepth]/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
	  //vV2_r1[i].push_back(vrxdepth_atten[i][idepth]*sin(vepsilon2_tx[i][idepth]/DEGRAD)*sin(vepsilon2_rx[i][idepth]/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);


	double theta1_Sclock_atrx,theta2_Sclock_atrx;
	//if (notflipped_atend>0.) {
	//if (i!=5) {
	  theta1_Sclock_atrx=vrxdepthE_theta1_Sclock[i][idepth];
	  theta2_Sclock_atrx=vrxdepthE_theta2_Sclock[i][idepth];
      
	  //}
	  //else {
	  //theta2_Sclock_atrx=vrxdepthE_theta1_Sclock[i][idepth];
	  //theta1_Sclock_atrx=vrxdepthE_theta2_Sclock[i][idepth];
	  //}


	//	vV1_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(vrxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
	//vV2_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(vrxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);

	// if (i==5 || i==0) {
	//   cout << "A" << i+1 << ", depth is " << vdepth[i][idepth] << ": time to fill the V's.\n";
	//   cout << "rx thetas are " << theta1_Sclock_atrx << "\t" << theta2_Sclock_atrx << "\n";
	//   cout << "tx thetas are " << vtxdepthE_theta1_Sclock[i][idepth] << "\t" << vtxdepthE_theta2_Sclock[i][idepth] << "\n";
	// }

	  // this is the real one, for r1.
	//**I think this is the real one.
	  vV1_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
	  vE1_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]);
	  //** this is the alternate
	  //vV1_r1[i].push_back(vrxdepth_atten[i][idepth]*sin(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);

	//**I think this is the real one.
	  vV2_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);
	  vE2_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]);
	  //** this is the alternate
	  //vV2_r1[i].push_back(vrxdepth_atten[i][idepth]*sin(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);

	  //	  cout << "vV2_r1 is " << vV2_r1[i][idepth] << "\n";


	  //vV1_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(theta1_Sclock_atrx/DEGRAD));
	  //vV2_r1[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD));

	  // this is an attempt to take into account that arianna
	  // only uses the LPDAs, which are only measuring polarizations
	  // in the horizontal plane.
	  // so for each component we take the cos(theta) where theta is
	  // zenith angle of shat
	  // since beam factors are sin(theta) we take cos=sqrt(1-sin*sin).
	  // I think this only needs to be applied to r1.
	  //cout << "push back the first one.\n";
	  vV1_r1_lpda[i].push_back(vV1_r1[i][idepth]*sqrt(1.-vrxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]));
	  vV2_r1_lpda[i].push_back(vV2_r1[i][idepth]*sqrt(1.-vrxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]));
	  //	  cout << "done pushing back the first one.\n";
	  



	vV1squared_r1[i].push_back(vV1_r1[i][idepth]*vV1_r1[i][idepth]);
	vV2squared_r1[i].push_back(vV2_r1[i][idepth]*vV2_r1[i][idepth]);
	vV1V2_r1[i].push_back(vV1_r1[i][idepth]*vV2_r1[i][idepth]);

	//	vV1_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vepsilon1_tx[i][idepth]/DEGRAD)*sin(vepsilon1_rx[i][idepth]/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
	//vV2_r2[i].push_back(vrxdepth_atten[i][idepth]*sin(vepsilon2_tx[i][idepth]/DEGRAD)*cos(vepsilon2_rx[i][idepth]/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);
   
	//	vV1_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(vrxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
	//vV2_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(vrxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);

	// this is the real one, for r2.
	//vV1_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
	//vV2_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);

	//	cout << "about to do the r2s.\n";

	//**I think this is the real one.
	vV1_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);
	vE1_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*sin(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]);
	//** this is the alternate
	//vV1_r2[i].push_back(vrxdepth_atten[i][idepth]*sin(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD)*cos(theta1_Sclock_atrx/DEGRAD)*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth]);

	//	vV2_r2[i].push_back(vrxdepth_atten[i][idepth]*sin(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);
	//**I think this is the real one.
	vV2_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);
	vE2_r2[i].push_back(vrxdepth_atten[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*sin(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]);
	//** this is the alternate
	//vV2_r2[i].push_back(vrxdepth_atten[i][idepth]*sin(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD)*cos(theta2_Sclock_atrx/DEGRAD)*vtxdepth_beam2[i][idepth]*vrxdepth_beam2[i][idepth]);


	vV1squared_r2[i].push_back(vV1_r2[i][idepth]*vV1_r2[i][idepth]);
	vV2squared_r2[i].push_back(vV2_r2[i][idepth]*vV2_r2[i][idepth]);
	vV1V2_r2[i].push_back(vV1_r2[i][idepth]*vV2_r2[i][idepth]);

	if (idepth==g_idepth[i]->Eval(-1000.))
	  cout << "station, depth, V1squared_r2, V2squared_r2, V1V2_r2 are " << i << "\t" << vV1squared_r2[i][idepth] << "\t" << vV2squared_r2[i][idepth] << "\t" << vV1V2_r2[i][idepth] << "\n"; 

	voppositeV1V2_r2[i].push_back(-1.*vV1_r2[i][idepth]*vV2_r2[i][idepth]);
	voppositeV1V2_r1[i].push_back(-1.*vV1_r1[i][idepth]*vV2_r1[i][idepth]);
      
	//	if (voppositeV1V2_r2[i][idepth]<0.) {
	//cout << "vopposite is negative.  i, idepth are " << i << "\t" << idepth << "\n";
	//}

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




	  for (int ifreq=0;ifreq<NFREQ;ifreq++) {
	    double thisfreq=vfreqs[ifreq];
	    double thissumphase=sumphase*thisfreq/freq;
	    
	    vspectra[i][idepth].push_back(vattens[i][idepth][ifreq]/vrxdepth_atten[i][idepth]*(venvelope_plus_r1[i][idepth]-4*vV1_r1[i][idepth]*vV2_r1[i][idepth]*sin(thissumphase)*sin(thissumphase)));
	    if (i==5 && idepth==g_idepth[i]->Eval(-1000.)) 
	      cout << "i, vattens[i][idepth], vrxdepth_atten[i][idepth], vspectra are " << i << "\t" << vattens[i][idepth][ifreq] << "\t" << vrxdepth_atten[i][idepth] << "\t" << vspectra[i][idepth][ifreq] << "\n";
	  }
      
	  // these are filled one per each pulser depth
	vtimediff[i].push_back(sumphase/(PI)*1./freq*1.E9); // this is in ns
     
	vnotflipped[i].push_back(notflipped_atend);
      
	vpower_r1[i].push_back(venvelope_plus_r1[i][idepth]-4*vV1_r1[i][idepth]*vV2_r1[i][idepth]*sin(sumphase)*sin(sumphase));
	vpoynting_r1[i].push_back(vSenvelope_plus_r1[i][idepth]-4*vE1_r1[i][idepth]*vE2_r1[i][idepth]*sin(sumphase)*sin(sumphase));
	
	vpower_r1_lpda[i].push_back(venvelope_plus_r1_lpda[i][idepth]-4*vV1_r1_lpda[i][idepth]*vV2_r1_lpda[i][idepth]*sin(sumphase)*sin(sumphase));
	vpower_r2[i].push_back(venvelope_plus_r2[i][idepth]-4*vV1_r2[i][idepth]*vV2_r2[i][idepth]*sin(sumphase)*sin(sumphase));
	vpoynting_r2[i].push_back(vSenvelope_plus_r2[i][idepth]-4*vE1_r2[i][idepth]*vE2_r2[i][idepth]*sin(sumphase)*sin(sumphase));


	//	if (vvenvelope_plus_r2[i][idepth]>0.1 || vvenvelope_minus_r2[i][idepth]>0.1) {
	//cout << "greater than 0.1.\n";
	//cout << "i, idepth, vvenvelops_plus_r1, vvenvelope_minus_r1, vV1_r1, vV2_r1 are " << vvenvelope_plus_r1[i][idepth] << "\t" << vvenvelope_minus_r1[i][idepth] << "\t" << vV1_r1[i][idepth] << "\t" << vV2_r1[i][idepth] << "\n";
	//cout << "i, idepth, vvenvelops_plus_r2, vvenvelope_minus_r2, vV1_r2, vV2_r2 are " << vvenvelope_plus_r2[i][idepth] << "\t" << vvenvelope_minus_r2[i][idepth] << "\t" << vV1_r2[i][idepth] << "\t" << vV2_r2[i][idepth] << "\n";
	//cout << "vpower_r1, vpower_r2 are " << vpower_r1[i][idepth] << "\t" << vpower_r2[i][idepth] << "\n";
	//cout << "ratio is " << sqrt(vpower_r2[i][idepth]/vpower_r1[i][idepth]) << "\n";
	//cout << "sumphase is " << sumphase << "\n";
	//cout << "vtxdepthE_theta1_Sclock, theta1_Sclock_atrx are " << vtxdepthE_theta1_Sclock[i][idepth] << "\t" << theta1_Sclock_atrx << "\n";
	//cout << "vtxdepthE_theta2_Sclock, theta2_Sclock_atrx are " << vtxdepthE_theta2_Sclock[i][idepth] << "\t" << theta2_Sclock_atrx << "\n";
	  //}


	vvoltage_r1[i].push_back(sqrt(vpower_r1[i][idepth]));
	vfield_r1[i].push_back(sqrt(vpoynting_r1[i][idepth]));
	vvoltage_r1_lpda[i].push_back(sqrt(vpower_r1_lpda[i][idepth]));
	vvoltage_r2[i].push_back(sqrt(vpower_r2[i][idepth]));
	vfield_r2[i].push_back(sqrt(vpoynting_r2[i][idepth]));

	if (i==0 && idepth==g_idepth[i]->Eval(-1000.) ||
	    i==5 && idepth==g_idepth[i]->Eval(-400.)) {
	  
	  cout << "A " << i+1 << ", depth is " << vdepth[i][idepth] << "\n";


	  cout << "before scalefactors:\n";
	  cout << "powers are " << vpower_r1[i][idepth] << "\t" << vpower_r2[i][idepth] << "\n";
	  cout << "voltages are " << vvoltage_r1[i][idepth] << "\t" << vvoltage_r2[i][idepth] << "\n";
	  cout << "lpda voltages are " << vvoltage_r1_lpda[i][idepth] << "\t" << vvoltage_r2[i][idepth] << "\n";



	  double scalefactor=0.2/sqrt(vvoltage_r1[i][idepth]*vvoltage_r1[i][idepth]+vvoltage_r2[i][idepth]*vvoltage_r2[i][idepth]);
	  double prefactor=vrxdepth_atten[i][idepth]*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth];
	  cout << "voltage r1 is " << scalefactor*vvoltage_r1[i][idepth] << "\n";
	  cout << "x, z components of voltage r1 are " << vtxdepth_beam1[i][idepth]*scalefactor*vvoltage_r1[i][idepth] << "\t" << sqrt(1.-vtxdepth_beam1[i][idepth]*vtxdepth_beam1[i][idepth])*scalefactor*vvoltage_r1[i][idepth] << "\n";

	  cout << "voltage r2 is " << scalefactor*vvoltage_r2[i][idepth] << "\n";
	  cout << "polarization angle is " << DEGRAD*atan2(vvoltage_r2[i][idepth],vvoltage_r1[i][idepth]) << "\n";
	  cout << "thetas_Sclock_atrx are " << theta1_Sclock_atrx << "\t" << theta2_Sclock_atrx << "\n";
	  cout << "diff is " << (theta1_Sclock_atrx-theta2_Sclock_atrx) << "\n";
	  cout << "thetas_Sclock_attx are " << vtxdepthE_theta1_Sclock[i][idepth] << "\t" << vtxdepthE_theta2_Sclock[i][idepth] << "\n";
	  cout << "diff is " << (vtxdepthE_theta1_Sclock[i][idepth]-vtxdepthE_theta2_Sclock[i][idepth]) << "\n";

	  // these are for drawing
	  cout << "vV1_r1, vV2_r1 are " << vV1_r1[i][idepth]/prefactor << "\t" << vV2_r1[i][idepth]/prefactor << "\n";
	  cout << "vV1_r2, vV2_r2 are " << vV1_r2[i][idepth]/prefactor << "\t" << vV2_r2[i][idepth]/prefactor << "\n";
	  
	  cout << "ray 1 at tx is " << scalefactor*vrxdepth_atten[i][idepth]*vtxdepth_beam1[i][idepth]*cos(vtxdepthE_theta1_Sclock[i][idepth]/DEGRAD) << "\n";
	  cout << "ray 2 at tx is " << scalefactor*vrxdepth_atten[i][idepth]*vtxdepth_beam1[i][idepth]*cos(vtxdepthE_theta2_Sclock[i][idepth]/DEGRAD) << "\n";
	  cout << "fraction of wavelength is " << sumphase/(2.*PI) << "\n";
	  cout << "vpower_r1+vpower_r2 \t" <<vpower_r1[i][idepth]+vpower_r2[i][idepth] << "\n";
	  cout << "prefactor is " << pow(vrxdepth_atten[i][idepth]*vtxdepth_beam1[i][idepth]*vrxdepth_beam1[i][idepth],2) << "\n";
	  cout << "three factors are " << vrxdepth_atten[i][idepth] << "\t" << vtxdepth_beam1[i][idepth] << "\t" << vrxdepth_beam1[i][idepth] << "\n";
	  
	  
	}
	if (i==5) {

	  // cout << "i, depth are " << i << "\t" << vdepth[i][idepth] << "\n";
	  // cout << "angles at rx are " << theta1_Sclock_atrx << "\t" << theta2_Sclock_atrx << "\n";
	  // cout << "angles at tx are " << vtxdepthE_theta1_Sclock[i][idepth] << "\t" << vtxdepthE_theta2_Sclock[i][idepth] << "\n"; 
	  // cout << "vV2_r1, vV1_r1 are " << vV2_r1[i][idepth] << "\t" << vV1_r1[i][idepth] << "\n";
	  // cout << "vV2_r2, vV1_r2 are " << vV2_r2[i][idepth] << "\t" << vV1_r2[i][idepth] << "\n";
	  // cout << "venvelope_plus_r1, venvelope_plus_r2 are " << venvelope_plus_r1[i][idepth] << "\t" << venvelope_plus_r2[i][idepth] << "\n";
	  // cout << "cross terms, sumphase are " << 4*vV1_r1[i][idepth]*vV2_r1[i][idepth] << "\t" <<  4*vV1_r2[i][idepth]*vV2_r2[i][idepth] << "\t" << sumphase << "\n";
	  // cout << "vvoltage_r2, vvoltage_r1 are " << vvoltage_r2[i][idepth] << "\t" << vvoltage_r1[i][idepth] << "\n";

	  // cout << "polarization angle is " << atan2(vvoltage_r2[i][idepth],vvoltage_r1[i][idepth])*DEGRAD << "\n";


	}

	//	if (i!=5) {
	vpolarization_Psi_rx[i].push_back(atan2(vvoltage_r2[i][idepth],vvoltage_r1[i][idepth])*DEGRAD);
	vpolarization_Omega_rx[i].push_back(acos(vvoltage_r1[i][idepth]/(sqrt(vvoltage_r1[i][idepth]*vvoltage_r1[i][idepth]+vvoltage_r2[i][idepth]*vvoltage_r2[i][idepth])))*DEGRAD);
	vEpolarization_Psi_rx[i].push_back(atan2(vfield_r2[i][idepth],vfield_r1[i][idepth])*DEGRAD);
	vEpolarization_Omega_rx[i].push_back(acos(vfield_r1[i][idepth]/(sqrt(vfield_r1[i][idepth]*vfield_r1[i][idepth]+vfield_r2[i][idepth]*vfield_r2[i][idepth])))*DEGRAD);


	  //	}
	  //else {
	  //vpolarization_Psi_rx[i].push_back(atan2(vvoltage_r2[i][idepth],vvoltage_r1_lpda[i][idepth])*DEGRAD);
	  //vpolarization_Omega_rx[i].push_back(acos(vvoltage_r1_lpda[i][idepth]/(sqrt(vvoltage_r1_lpda[i][idepth]*vvoltage_r1_lpda[i][idepth]+vvoltage_r2[i][idepth]*vvoltage_r2[i][idepth])))*DEGRAD);


	  //	}
	
	//	cout << "A" << i+1 << ": r1 envelope plus, minus are " << venvelope_plus_r1[i][idepth] << "\t" << venvelope_minus_r1[i][idepth] << "\n";

	//	cout << "A" << i+1 << ": ended loop over steps along the path.\n";
	// end loop over solutions
	// end loop over idepth
	
	vsumlength[i].push_back(sumlength);


	//	double *getangles=GetDirectRayPar(z0, x1, z1);

  //	double launch_angle=get_angle(posstation,pospulser, solutions[0][1], N_ICE, DELTA_N, Z_0);
  //	double receive_angle=get_receive_angle(pospulser,posstation, solutions[0][1], N_ICE, DELTA_N, Z_0);
      

      

	//receive_angle=PI-receive_angle;

	vector<double> nvec_thisdepth;
	nvec_thisdepth.resize(3);
	nvec_thisdepth[0]=gn1->Eval(vdepth[i][idepth]);
	nvec_thisdepth[1]=gn2->Eval(vdepth[i][idepth]);
	nvec_thisdepth[2]=gn3->Eval(vdepth[i][idepth]);




	// cout << "at " << vdepth[i][idepth] << " depth, n's are " << nvec_thisdepth[0] << "\t" << nvec_thisdepth[1] << "\t" << nvec_thisdepth[2] << "\n";
	// cout << "posstation is " << posstation[0] << "\t" << posstation[1] << "\t" << posstation[2] << "\n";
	// cout << "pospulser is " << pospulser[0] << "\t" << pospulser[1] << "\t" << pospulser[2] << "\n";

	// cout << "A" << i << ", depth " << vdepth[i][idepth] << ": launch_angle, receive_angle, and V are " << launch_angle*DEGRAD << "\t" << receive_angle*DEGRAD << "\t" << getV(nvec_thisdepth) << "\n";

	// cout << "depth, sumphase are " << vdepth[i][idepth] << "\t" << sumphase << "\n";
	

	if (DORAYTRACING) {


	  TVector3 plusz(0.,0.,1.);
	  //	  TVector3 minusz(0.,0.,-1.);
	  TVector3 vrotate=plusz.Cross(yhat);
	  //	  cout << "vrotate is " << vrotate[0] << "\t" << vrotate[1] << "\t" << vrotate[2] << "\n";
	  if (vrotate.Mag()<HOWSMALLISTOOSMALL)
	    cout << "vrotate mag is " << vrotate.Mag() << "\n";

	  vrotate.SetMag(1.);
	  rhat[i]=plusz;
	  rhat[i].Rotate(launch_angle,vrotate);
	  //cout << "rhat after rotation is " << rhat[i][0] << "\t" << rhat[i][1] << "\t" << rhat[i][2] << "\n";
	  

	  rhat_receive[i]=plusz;
	  rhat_receive[i].Rotate(receive_angle,vrotate);
	  //	  cout << "rhat_receive after rotation is " << rhat_receive[i][0] << "\t" << rhat_receive[i][1] << "\t" << rhat_receive[i][2] << "\n";

	} // end do ray tracing
	else {
	  // this is if no ray tracing.

	  rhat[i].SetX(station_coords[i][0]-pulser_coords[0]); // this will be a negative number
	  rhat[i].SetY(station_coords[i][1]-pulser_coords[1]);
	  rhat[i].SetZ(station_depths[i]-vdepth[i][idepth]);
	  //rhat[i].SetZ(station_depths[i]-vdepth[i][0]); // this is a difference between two negative numbers.
	  
	//}
	  
	  rhat_receive[i].SetX(rhat[i][0]);
	  rhat_receive[i].SetY(rhat[i][1]);
	  rhat_receive[i].SetZ(rhat[i][2]);
	  //rhat_receive[i].SetZ(0.);

	}


	if (rhat[i].Mag()<HOWSMALLISTOOSMALL)
	  cout << "rhat[i] mag is " << rhat[i].Mag() << "\n";
	rhat[i].SetMag(1.);

	if (rhat_receive[i].Mag()<HOWSMALLISTOOSMALL)
	  cout << "rhat_receive[i] mag is " << rhat_receive[i].Mag() << "\n";
	rhat_receive[i].SetMag(1.);

	//	cout << "rhat after rotate is " << rhat[i][0] << "\t" << rhat[i][1] << "\t" << rhat[i][2] << "\n";
	//cout << "rhat_receive after rotate is " << rhat_receive[i][0] << "\t" << rhat_receive[i][1] << "\t" << rhat_receive[i][2] << "\n";

	//	cout << "path_length, sumphase are " << path_length << "\t" << sumphase << "\n";
      } // uzair's ray tracer gives solutions.
      else {


      //if (idepth==64) {

      //vdepth[i][idepth]=-1000.;
	
      //      }
      
  //  r[i].clear();
  //  extraordinary[i].clear();
      //for (int j=0;j<3;j++) {
  	rhat[i].SetX(station_coords[i][0]-pulser_coords[0]); // this will be a negative number
  	rhat[i].SetY(station_coords[i][1]-pulser_coords[1]);
  	rhat[i].SetZ(station_depths[i]-vdepth[i][idepth]);
	//rhat[i].SetZ(station_depths[i]-vdepth[i][0]); // this is a difference between two negative numbers.

	//}

      rhat_receive[i].SetX(rhat[i][0]);
      rhat_receive[i].SetY(rhat[i][1]);
      rhat_receive[i].SetZ(rhat[i][2]);
      //rhat_receive[i].SetZ(0.);



      }



      vreceiveangle[i].push_back(rhat_receive[i].Theta()*DEGRAD);
      vlaunchangle[i].push_back(rhat[i].Theta()*DEGRAD);
    
      // still in a loop over depths.

      if (rhat[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "rhat[i] mag is " << rhat[i].Mag() << "\n";

      rhat[i].SetMag(1.);

      if (rhat_receive[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "rhat_receive[i] mag is " << rhat_receive[i].Mag() << "\n";

      rhat_receive[i].SetMag(1.);


      /*
      // cout << "nominal rhat is " << rhat[i][0] << "\t" << rhat[i][1] << "\t" << rhat[i][2] << "\n";      
      // cout << "nominal rhat_receive is " << rhat_receive[i][0] << "\t" << rhat_receive[i][1] << "\t" << rhat_receive[i][2] << "\n";      

      double averagedepth=(station_depths[i]+vdepth[i][idepth])/2.; // the -1 in front makes this a positive number.
      //cout << "averagedepth, vdepthsmin, vdepthsmax are " << averagedepth << "\t" << vdepths_n1[1] << "\t" << vdepths_n1[vdepths_n1.size()-2] << "\n";

      if (averagedepth>vdepths_n1[1] && averagedepth<vdepths_n1[vdepths_n1.size()-2]) {
	nvec[0]=gn1->Eval(averagedepth);
	nvec[1]=gn2->Eval(averagedepth);
	nvec[2]=gn3->Eval(averagedepth);
      }
      else if (averagedepth<vdepths_n1[1]) {
	nvec[0]=gn1->Eval(vdepths_n1[0]);
	nvec[1]=gn2->Eval(vdepths_n2[0]);
	nvec[2]=gn3->Eval(vdepths_n3[0]);
      }
      else if (averagedepth>vdepths_n1[vdepths_n1.size()-2]) {
	nvec[0]=gn1->Eval(vdepths_n1[vdepths_n1.size()-1]);
	nvec[1]=gn2->Eval(vdepths_n2[vdepths_n1.size()-1]);
	nvec[2]=gn3->Eval(vdepths_n3[vdepths_n1.size()-1]);
      }

      if (rhat[i].Mag()<1.E-8)
	cout << "before calling getDeltaN at place 3, rhat is " << rhat[i][0] << "\t" << rhat[i][1] << "\t" << rhat[i][2] << "\n";
      double deltan=getDeltaN(nvec,rhat[i],angle_iceflow,n_e1,n_e2,p_e1,p_e2);
    
      deltan_exp[i]=deltan;

      //      cout << "A" << i+1 << ": for biaxial case, nvec, delta n is " << nvec[0] << "\t" << nvec[1] << "\t" << nvec[2] << "\t" << deltan << "\n";
      
      // first find the directions of the e and o components at the rx for a signal traveling in a straight line

      nvec[0]=gn1->Eval(station_depths[i]);
      nvec[1]=gn2->Eval(station_depths[i]);
      nvec[2]=gn3->Eval(station_depths[i]);


      TVector3 khat1_atrx=rhat_receive[i]; // This is just a first guess for khat at the rx. rhat_receive is Shat which is not the same as khat
      TVector3 khat2_atrx=rhat_receive[i]; // This is just a first guess for khat at the rx. rhat_receive is Shat which is not the same as khat

      if (khat1_atrx.Mag()<1.E-8)
	cout << "before calling getDeltaN at place 4, khat1_atrx is " << khat1_atrx[0] << "\t" << khat1_atrx[1] << "\t" << khat1_atrx[2] << "\n";
      double deltan_attherx=getDeltaN(nvec,khat1_atrx,angle_iceflow,n_e1,n_e2,p_e1,p_e2); // this is just a guess, it isn't right since rhat_receive isn't k. p_e1, p_e2 are a guess for the direction of the D's.

      TVector3 Ehat_khat1_e1_atrx;
      TVector3 Ehat_khat2_e1_atrx;
      TVector3 Ehat_khat1_e2_atrx;
      TVector3 Ehat_khat2_e2_atrx;

      khat1_atrx=getNewkandE(nvec,angle_iceflow,n_e1,p_e1,khat1_atrx,
				     Ehat_khat1_e1_atrx); // Find E_e1 given p_e1.   khat_atrx as the input is just a guess.
      khat2_atrx=getNewkandE(nvec,angle_iceflow,n_e1,p_e1,khat2_atrx,
				     Ehat_khat2_e1_atrx); // Find E_e1 given p_e1.   khat_atrx as the input is just a guess.

      khat1_atrx=getNewkandE(nvec,angle_iceflow,n_e2,p_e2,khat1_atrx,
				     Ehat_khat1_e2_atrx); // Find E_e1 given p_e1.   khat_atrx as the input is just a guess.
      khat2_atrx=getNewkandE(nvec,angle_iceflow,n_e2,p_e2,khat2_atrx,
				     Ehat_khat2_e2_atrx); // Find E_e1 given p_e1.   khat_atrx as the input is just a guess.

      //find what Shat would be
      if (Ehat_khat1_e1_atrx.Mag()<HOWSMALLISTOOSMALL)
	cout << "Ehat_khat1_e1_atrx mag is " << Ehat_khat1_e1_atrx.Mag() << "\n";
      Ehat_khat1_e1_atrx.SetMag(1.);

      if (Ehat_khat2_e1_atrx.Mag()<HOWSMALLISTOOSMALL)
	cout << "Ehat_khat2_e1_atrx mag is " << Ehat_khat2_e1_atrx.Mag() << "\n";
      Ehat_khat2_e1_atrx.SetMag(1.);

      if (Ehat_khat1_e2_atrx.Mag()<HOWSMALLISTOOSMALL)
	cout << "Ehat_khat1_e2_atrx mag is " << Ehat_khat1_e2_atrx.Mag() << "\n";
      Ehat_khat1_e2_atrx.SetMag(1.);

      if (Ehat_khat2_e2_atrx.Mag()<HOWSMALLISTOOSMALL)
	cout << "Ehat_khat2_e2_atrx mag is " << Ehat_khat2_e2_atrx.Mag() << "\n";
      Ehat_khat2_e2_atrx.SetMag(1.);

      TVector3 Hhat1_atrx=p_e1.Cross(Ehat_khat1_e1_atrx);

      if (Hhat1_atrx.Mag()<HOWSMALLISTOOSMALL)
	cout << "Hhat1_atrx mag is " << Hhat1_atrx.Mag() << "\n";
      Hhat1_atrx.SetMag(1.);


      TVector3 Hhat2_atrx=p_e2.Cross(Ehat_khat2_e1_atrx);

      if (Hhat2_atrx.Mag()<HOWSMALLISTOOSMALL)
	cout << "Hhat2_atrx mag is " << Hhat2_atrx.Mag() << "\n";
      Hhat2_atrx.SetMag(1.);



      TVector3 Shat1_atrx=Ehat_khat1_e1_atrx.Cross(Hhat1_atrx);
      if (Shat1_atrx.Mag()<HOWSMALLISTOOSMALL)
	cout << "Shat1_atrx mag is " << Shat1_atrx.Mag() << "\n";
      Shat1_atrx.SetMag(1.);

      TVector3 rotatethem=Shat1_atrx.Cross(rhat_receive[i]); 
      double angle_rotatethem=acos(Shat1_atrx.Dot(rhat_receive[i])/rhat_receive[i].Mag()/Shat1_atrx.Mag());
      //      cout << "angle_rotatethem is " << angle_rotatethem << "\n";
      Shat1_atrx.Rotate(angle_rotatethem,rotatethem); // Step 7.  Rotate E_e1 into E_e1^prime
      Ehat_khat1_e1_atrx.Rotate(angle_rotatethem,rotatethem); // Step 7.  Rotate E_e1 into E_e1^prime


      TVector3 Shat2_atrx=Ehat_khat2_e1_atrx.Cross(Hhat2_atrx);
      if (Shat2_atrx.Mag()<HOWSMALLISTOOSMALL)
	cout << "Shat2_atrx mag is " << Shat2_atrx.Mag() << "\n";
      Shat2_atrx.SetMag(1.);

      rotatethem=Shat2_atrx.Cross(rhat_receive[i]); 
      angle_rotatethem=acos(Shat2_atrx.Dot(rhat_receive[i])/rhat_receive[i].Mag()/Shat2_atrx.Mag());
      Shat2_atrx.Rotate(angle_rotatethem,rotatethem); // Step 7.  Rotate E_e1 into E_e1^prime
      Ehat_khat2_e1_atrx.Rotate(angle_rotatethem,rotatethem); // Step 7.  Rotate E_e1 into E_e1^prime

      
      TVector3 p_e1_khat1,p_e1_khat2;
      TVector3 p_e2_khat1,p_e2_khat2;
      double n_e1_khat1,n_e2_khat1,n_e1_khat2,n_e2_khat2;

      if (khat1_atrx.Mag()<1.E-8)
	cout << "before calling getDeltaN at place 5, khat1_atrx is " << khat1_atrx[0] << "\t" << khat1_atrx[1] << "\t" << khat1_atrx[2] << "\n";
      deltan_attherx=getDeltaN(nvec,khat1_atrx,angle_iceflow,n_e1_khat1,n_e2_khat1,p_e1_khat1,p_e2_khat1); // this is just a guess, it isn't right since rhat_receive isn't k
      if (khat2_atrx.Mag()<1.E-8)
	cout << "before calling getDeltaN at place 6, khat2_atrx is " << khat2_atrx[0] << "\t" << khat2_atrx[1] << "\t" << khat2_atrx[2] << "\n";
      deltan_attherx=getDeltaN(nvec,khat2_atrx,angle_iceflow,n_e1_khat2,n_e2_khat2,p_e1_khat2,p_e2_khat2); // this is just a guess, it isn't right since rhat_receive isn't k

      TVector3 vtemp=n_e1_khat1*p_e1_khat1+n_e1_khat2*p_e1_khat2;
      if (vtemp.Mag()<HOWSMALLISTOOSMALL)
	cout << "vtemp mag is " << vtemp.Mag() << "\n";

      vtemp.SetMag(1.);
      p_e1=vtemp;

      vtemp=n_e2_khat1*p_e2_khat1+n_e2_khat2*p_e2_khat2;
      if (vtemp.Mag()<HOWSMALLISTOOSMALL)
	cout << "vtemp mag is " << vtemp.Mag() << "\n";
      vtemp.SetMag(1.);
      p_e2=vtemp;

      if (WHICHPOL==1) {
	p_o_rx[i]=p_e2;
	p_e_rx[i]=p_e1;
      }
      else if (WHICHPOL==0) {
	p_o_rx[i]=p_e1;
	p_e_rx[i]=p_e2;
      }

      nvec[0]=gn1->Eval(vdepth[i][idepth]);
      nvec[1]=gn2->Eval(vdepth[i][idepth]);
      nvec[2]=gn3->Eval(vdepth[i][idepth]);


      if (rhat[i].Mag()<1.E-8)
	cout << "before calling getDeltaN at place 7, rhat is " << rhat[i][0] << "\t" << rhat[i][1] << "\t" << rhat[i][2] << "\n";
      double deltan_atthetx=getDeltaN(nvec,rhat[i],angle_iceflow,n_e1,n_e2,p_e1,p_e2);

      if (WHICHPOL==1) {
	p_o[i]=p_e2;
	p_e[i]=p_e1;
      }
      else if (WHICHPOL==0) {
	p_o[i]=p_e1;
	p_e[i]=p_e2;
      }

      //      p_o[i] = ordinary.Cross(rhat[i]);
      //p_o[i].SetMag(1.);

      //      p_e[i] = rhat[i].Cross(p_o[i]);
      //p_e[i] = (rhat[i].Dot(ordinary))*ordinary;
      //p_e[i].SetMag(1.);
      //p_o[i] = rhat[i]-p_e[i];
      //p_o[i].SetMag(1.);

      temp = tx_orientation.Cross(rhat[i]);
      Pt[i]=rhat[i].Cross(temp);

      if (Pt[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "Pt[i] mag is " << Pt[i].Mag() << "\n";

      Pt[i].SetMag(1.);



      //begin iterative stuff
      double n_e1,n_e2;
      TVector3 khat_0=rhat[i]; // Step 1.  This is rhat_launch.


      TVector3 E=Pt[i];
      TVector3 D=rotateE(epsilon, angle_iceflow, E);  // D=epsilon*E



      TVector3 khat_ray1_0=rotateE(epsilon, angle_iceflow, khat_0); // Step 2.  First guesses khat_ray1_0 = khat_ray2_0 = khat_0 are rhat_launch rotated by epsilon. 

      TVector3 khat_ray2_0=khat_ray1_0;

      if (khat_ray1_0.Mag()<1.E-8)
	cout << "before calling getDeltaN at place 8, rhat is " << khat_ray1_0[0] << "\t" << khat_ray1_0[1] << "\t" << khat_ray1_0[2] << "\n";

      getDeltaN(nvec,khat_ray1_0,angle_iceflow,
		n_e1,n_e2,p_e1,p_e2); // Step 3.  p_e1 and p_e2 are the directions of the eigenvectors.     

      TVector3 D_norm=D;
      if (D_norm.Mag()<HOWSMALLISTOOSMALL)
	cout << "D_norm mag is " << D_norm.Mag() << "\n";
      D_norm.SetMag(1.);  // direction of D

      double component_d1=p_e1.Dot(D_norm); // component of D direction that is parallel with eigenvector 1 of ray 1
      //cout << "component_d1 is " << component_d1 << "\n";

      double component_d2=p_e2.Dot(D_norm); // component of D direction that is parallel to eigenvector 2 of ray 1
      //cout << "component_d2 are " << component_d2 << "\n";      

      TVector3 D_e1=n_e1*n_e1*component_d1*p_e1; // Finish Step 3.  This is the displacement vector D for eigenvalue 1 of ray 1

      TVector3 E_e1,E_e2;

      // cout << "first iteration:\n";
      // cout << "khat_0 (perp to Pt) is " << khat_0[0] << "\t" << khat_0[1] << "\t" << khat_0[2] << "\n";
      // cout << "e1:\n";
      // cout << "finding khat_ray1_1 and E_e1.\n";
      TVector3 khat_ray1_1=getNewkandE(nvec,angle_iceflow,n_e1,D_e1,khat_ray1_0,
				     E_e1); // Step 4.  Find E_e1 given khat_ray1_1 and D_e1 of ray 1

      // cout << "khat_ray1_1 is " << khat_ray1_1[0] << "\t" << khat_ray1_1[1] << "\t" << khat_ray1_1[2] << "\n";      

      // cout << "E_e1 is " << E_e1[0] << "\t" << E_e1[1] << "\t" << E_e1[2] << "\n";

      TVector3 D_e2=rotateE(epsilon,angle_iceflow,E)-D_e1; // Step 5.  Find D_e2^prime, which is the difference between the "true" D and D for eigenvalue 1 of ray 1.  this is not in general the other eigenvalue of the first ray.

      // this was all just for the first eigenvalue.  see if that is the one we should be using.


      //      cout << "finding khat_ray2_1 and E_e2.\n";
      // find khat and E_e2 corresponding to the eigenvalue of interest of the second ray
      TVector3 khat_ray2_1=getNewkandE(nvec,angle_iceflow,n_e2,D_e2,khat_ray2_0,
				     E_e2); // Step 6.  Find E_e2^prime

      //cout << "khat_ray2_1 is " << khat_ray2_1[0] << "\t" << khat_ray2_1[1] << "\t" << khat_ray2_1[2] << "\n";      


      //     cout << "angle between two rays is " << acos(khat_ray1_1.Dot(khat_ray2_1)/khat_ray1_1.Mag()/khat_ray2_1.Mag())*DEGRAD << "\n";
      // this seems out of place.
     TVector3 khat_ray1_0_unit=khat_ray1_0;
     if (khat_ray1_0_unit.Mag()<HOWSMALLISTOOSMALL)
       cout << "khat_ray1_0_unit mag is " << khat_ray1_0_unit.Mag() << "\n";

     khat_ray1_0_unit.SetMag(1.);

     //     cout << "recall that khat_ray1_0_unit is " << khat_ray1_0_unit[0] << "\t" << khat_ray1_0_unit[1] << "\t" << khat_ray1_0_unit[2] << "\n";



      // cout << "E_e2 is " << E_e2[0] << "\t" << E_e2[1] << "\t" << E_e2[2] << "\n";

      // cout << "D_e1 is " << D_e1[0] << "\t" << D_e1[1] << "\t" << D_e1[2] << "\n";
      // cout << "D_e2 is " << D_e2[0] << "\t" << D_e2[1] << "\t" << D_e2[2] << "\n";

      // cout << "angle between D_e1 and D_e2 is " << DEGRAD*acos(D_e1.Dot(D_e2)) << "\n";

      // cout << "angle between E_e1 and E_e2 is " << DEGRAD*acos(E_e1.Dot(E_e2)) << "\n";


      TVector3 EminusE_e1=E-E_e1;
      TVector3 EminusE_e2=E-E_e2; // Finish Step 6.  Find E_e2^prime

      // cout << "E-E_e1 is " << EminusE_e1[0] << "\t" << EminusE_e1[1] << "\t" << EminusE_e1[2] << "\n";

      // cout << "E-E_e2 is " << EminusE_e2[0] << "\t" << EminusE_e2[1] << "\t" << EminusE_e2[2] << "\n";

      TVector3 E_sum=E_e1+E_e2;


      // cout << "ray 1: satisfies wave equation if this is zero: \t" << (-1.*n_e1*n_e1*E_e1+(khat_ray1_1.Dot(E_e1))*n_e1*n_e1*khat_ray1_1 + D_e1).Mag() << "\n";

      // cout << "ray 2: satisfies wave equation if this is zero: \t" << (-1.*n_e2*n_e2*E_e2+(khat_ray2_1.Dot(E_e2))*n_e2*n_e2*khat_ray2_1 + D_e2).Mag() << "\n";

      // cout << "but I think we knew this wouldn't work because D is not equal to D_sum?  like D_e2=D-D_e1 wasn't quite right?\n";

      // rotate E_e1 into EminusE_e2.
      
      TVector3 rotatedE=E_e1.Cross(EminusE_e2); // make sure E_e1 aligns with E-E_e2. 
      double angle_rotateE=acos(E_e1.Dot(EminusE_e2)/E_e1.Mag()/EminusE_e2.Mag());
      E_e1.Rotate(angle_rotateE,rotatedE); // Step 7.  Rotate E_e1 into E_e1^prime
      // now E_e1+E_e2 will not longer be in the direction of Pt[i]
      D_e1.Rotate(angle_rotateE,rotatedE); // Finish Step 7.  Rotate D_e1 by the same matrix.

      // cout << "after rotating E_e1 into EminuE_e2, E_e1 is " << E_e1[0] << "\t" << E_e1[1] << "\t" << E_e1[2] << "\n";
      
      // cout << "after rotating D_e1 by same matrix, D_e1 is " << D_e1[0] << "\t" << D_e1[1] << "\t" << D_e1[2] << "\n";


      // cout << "from rotating E, D is " << D[0] << "\t" << D[1] << "\t" << D[2] << "\n";
      // cout << "mag of D is " << D.Mag() << "\n";

      // cout << "angle between E and D is " << acos(E.Dot(D)/E.Mag()/D.Mag())*DEGRAD << "\n";

      TVector3 D_sum=D_e1+D_e2; // Step 8.  Find D^prime(0)

      //      cout << "from summing, D is " << D_sum[0] << "\t" << D_sum[1] << "\t" << D_sum[2] << "\n";

      // now rotate D_sum into D

      //      cout << "mag of D_sum is " << D_sum.Mag() << "\n";


      TVector3 rotateD=D_sum.Cross(D);
      double angle_rotateD=acos(D_sum.Dot(D)/D_sum.Mag()/D.Mag());

      D_sum.Rotate(angle_rotateD,rotateD); // Finish Step 8.  Rotate D^prime(0) to align with D
      E_e1.Rotate(angle_rotateD,rotateD); // Finish Step 8.  Rotate D^prime(0) to align with D
      E_e2.Rotate(angle_rotateD,rotateD); // Finish Step 8.  Rotate D^prime(0) to align with D

      //     cout << "after rotation, D_sum is " << D_sum[0] << "\t" << D_sum[1] << "\t" << D_sum[2] << "\n";

      TVector3 tempvec=khat_ray1_1;
      tempvec.Rotate(angle_rotateE,rotatedE);
      tempvec.Rotate(angle_rotateD,rotateD);
      TVector3 khat_ray1_2=tempvec; // Step 9.  Our new guess for wave vector of first ray

      tempvec=khat_ray2_1;
      tempvec.Rotate(angle_rotateD,rotateD);
      TVector3 khat_ray2_2=tempvec; // Step 9.  Our new guess for wave vector of second ray

      // should rotate D_e1 and D_e2 by the same angle_rotateD about vector rotateD, find E_e1 and E_e2 given D_e1, D_e2, khat_ray1_2, and khat_ray2_2, and see how close E_e1+E_e2 is to Pt[i].

      D_e1.Rotate(angle_rotateD,rotateD);
      D_e2.Rotate(angle_rotateD,rotateD);

      TVector3 E_e1_test;
      TVector3 E_e2_test;

      TVector3 khat_ray1_2_test=getNewkandE(nvec,angle_iceflow,n_e1,D_e1,khat_ray1_2,
					  E_e1_test); // Step 6.  Find E_e2^prime
      TVector3 khat_ray2_2_test=getNewkandE(nvec,angle_iceflow,n_e2,D_e2,khat_ray2_2,
					  E_e2_test); // Step 6.  Find E_e2^prime

      TVector3 E_test=E_e1_test+E_e2_test;
      TVector3 E_sum_new=E_e1+E_e2;
      TVector3 D_sum_test=D_e1+D_e2;

      double field1=E_e1.Mag();
      double field2=E_e2.Mag();

      // cout << "angle between D and D_sum_test is " << acos(D.Dot(D_sum_test)/D.Mag()/D_sum_test.Mag())*DEGRAD << "\n";
      // cout << "angle between E and E_test is " << acos(E.Dot(E_test)/E.Mag()/E_test.Mag())*DEGRAD << "\n";
      // cout << "angle between E and E_sum_new is " << acos(E.Dot(E_sum_new)/E.Mag()/E_sum_new.Mag())*DEGRAD << "\n";
      // cout << "angle between E and E_sum was " << acos(E.Dot(E_sum)/E.Mag()/E_sum.Mag())*DEGRAD << "\n";
      // cout << "khat_ray1_2 is " << khat_ray1_2[0] << "\t" << khat_ray1_2[1] << "\t" << khat_ray1_2[2] << "\n";
      // cout << "khat_ray2_2 is " << khat_ray2_2[0] << "\t" << khat_ray2_2[1] << "\t" << khat_ray2_2[2] << "\n";
      // cout << "angle between two rays is " << acos(khat_ray1_2.Dot(khat_ray2_2)/khat_ray1_2.Mag()/khat_ray2_2.Mag())*DEGRAD << "\n";
      
      // cout << "angle between khat_0 and khat_ray1_2 are " << acos(khat_0.Dot(khat_ray1_2)/khat_0.Mag()/khat_ray1_2.Mag())*DEGRAD << "\n";
      

      vangle_khat_0_khat_1_2[i].push_back(acos(khat_0.Dot(khat_ray1_2)/khat_0.Mag()/khat_ray1_2.Mag())*DEGRAD);
      vangle_khat_0_khat_2_2[i].push_back(acos(khat_0.Dot(khat_ray2_2)/khat_0.Mag()/khat_ray2_2.Mag())*DEGRAD);
      vangle_khat_1_2_khat_2_2[i].push_back(acos(khat_ray1_2.Dot(khat_ray2_2)/khat_ray1_2.Mag()/khat_ray2_2.Mag())*DEGRAD);
    

      // cout << "E is " << E[0] << "\t" << E[1] << "\t" << E[2] << "\n";

      // cout << "mag of E is " << E.Mag() << "\n";

      // cout << "E_sum is " << E_sum[0] << "\t" << E_sum[1] << "\t" << E_sum[2] << "\n";
      // cout << "mag of E_sum is " << E_sum.Mag() << "\n";

      // cout << "angle between E and E_sum is " << acos(E.Dot(E_sum)/E.Mag()/E_sum.Mag())*DEGRAD << "\n";

      TVector3 Hhat_e1=1./E_e1.Mag()/D_e1.Mag()*D_e1.Cross(E_e1);
      Hhat_e1.SetMag(1.);
      TVector3 Shat_e1=E_e1.Cross(Hhat_e1);
      Shat_e1.SetMag(1.);

      //      cout << "Shat_e1 before is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
      if (Shat_e1.Dot(khat_ray1_2)<0.)
	Shat_e1 = -1.*Shat_e1;
      //cout << "Shat_e1 after is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";

      vangle_Shat_e1_khat_1_2[i].push_back(acos(Shat_e1.Dot(khat_ray1_2)/Shat_e1.Mag()/khat_ray1_2.Mag())*DEGRAD);


      TVector3 Hhat_e2=1./E_e2.Mag()/D_e2.Mag()*D_e2.Cross(E_e2);
      Hhat_e2.SetMag(1.);
      TVector3 Shat_e2=E_e2.Cross(Hhat_e2);
      Shat_e2.SetMag(1.);

      //      cout << "Shat_e2 before is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";
      if (Shat_e2.Dot(khat_ray2_2)<0.)
	Shat_e2 = -1.*Shat_e2;
      //cout << "Shat_e2 after is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";

      vangle_Shat_e2_khat_2_2[i].push_back(acos(Shat_e2.Dot(khat_ray2_2)/Shat_e2.Mag()/khat_ray2_2.Mag())*DEGRAD);

      //      if (idepth==igreatestdepth[i]){
	//if (idepth==imostshallowdepth[i]){
      if (idepth==(int)(g_idepth[i]->Eval(depth_special))){
      //      if (idepth==vdepth[i].size()-1) {
  //      if (idepth==ispecial) {

	//	raypath[0][i].SetX(station_coords[5][0]);
	//raypath[0][i].SetY(station_coords[5][1]);
	raypath[0][i].SetX(0.); // indices are ray, station
	raypath[0][i].SetY(0.);
	raypath[0][i].SetZ(vdepth[i][idepth]);
	
	double radial_distance=0.;
	double stepray=200.; // 1m for now
	
	vector<double> nvec_thisstep;
	vector<double> nvec_nextstep;
	
	nvec_thisstep.resize(3);
	nvec_nextstep.resize(3);
	
	TVector3 E_e1_thisstep=E_e1;
	TVector3 D_e1_thisstep=D_e1;
	TVector3 Hhat_e1_thisstep=Hhat_e1;
	
	//	cout << "A" << i+1 << ": radial_distance, horizontal_distances[i] are " << radial_distance << "\t" << horizontal_distances[i] << "\n";
	while (radial_distance<horizontal_distances[i]) {

	  
	  vraypos[i][0].push_back(raypath[0][i][0]); // 0th ray, i station, x
	  vraypos[i][1].push_back(raypath[0][i][1]);
	  vraypos[i][2].push_back(raypath[0][i][2]);
	  
	  nvec_thisstep[0]=gn1->Eval(raypath[0][i][2]); // find indicatrix for each depth
	  nvec_thisstep[1]=gn2->Eval(raypath[0][i][2]);
	  nvec_thisstep[2]=gn3->Eval(raypath[0][i][2]);
	  



	  double n_e1_thisstep,n_e2_thisstep;
	  TVector3 p_e1_temp,p_e2_temp;
	  if (khat_ray1_2.Mag()<1.E-8)
	    cout << "before calling getDeltaN at place 9, khat_ray1_2 is " << khat_ray1_2[0] << "\t" << khat_ray1_2[1] << "\t" << khat_ray1_2[2] << "\n";

	  double deltan_thisstep=getDeltaN(nvec_thisstep,khat_ray1_2,angle_iceflow,
					   n_e1_thisstep,n_e2_thisstep,p_e1_temp,p_e2_temp); // using khat for ray 1 that we got from the iterative method, what is delta n, and eigenvectors for D.
	  // this khat should be adjusted for every Shat
	  // and actually it's khat that observes snell's law not Shat...


	  vraypath_ne1[i].push_back(n_e1_thisstep);
	  vraypath_ne2[i].push_back(n_e2_thisstep);

	  TVector3 zhat(0.,0.,1.);
	  double theta_i=acos(khat_ray1_2.Dot(zhat)); // incident angle for snell's law
	  
	  raypath[0][i]=raypath[0][i]+stepray*khat_ray1_2;
	  
	  nvec_nextstep[0]=gn1->Eval(raypath[0][i][2]);
	  nvec_nextstep[1]=gn2->Eval(raypath[0][i][2]);
	  nvec_nextstep[2]=gn3->Eval(raypath[0][i][2]);
	  
	  TVector3 epsilon_nextstep;
	  

	  epsilon_nextstep[0]=nvec_nextstep[0]*nvec_nextstep[0];
	  epsilon_nextstep[1]=nvec_nextstep[1]*nvec_nextstep[1];
	  epsilon_nextstep[2]=nvec_nextstep[2]*nvec_nextstep[2];

	  double n_e1_nextstep,n_e2_nextstep;
	  double deltan_nextstep=getDeltaN(nvec_nextstep,khat_ray1_2,angle_iceflow,
					   n_e1_nextstep,n_e2_nextstep,p_e1_temp,p_e2_temp);
	  
	  // Snell's law
	  double theta_r=asin(n_e1_thisstep*sin(theta_i)/n_e1_nextstep);
	  
	  TVector3 rotate_Snell=khat_ray1_2.Cross(zhat);
	  if (rotate_Snell.Mag()<HOWSMALLISTOOSMALL)
	    cout << "rotate_Snell mag is " << rotate_Snell.Mag() << "\n";

	  rotate_Snell.SetMag(1.);
	  TVector3 khat_e1_nextstep=khat_ray1_2; // problem is that the D we'd get from epsilon*E won't in general be an eigenvalue of khat_ray1_2.  so that's why we need the iterative stuff.
	  khat_e1_nextstep.Rotate(theta_r-theta_i,rotate_Snell);
	  //cout << "ne1_thisstep, ne1_nextstep, angle diff are " << n_e1_thisstep << "\t" << n_e1_nextstep << "\t" << DEGRAD*(theta_r-theta_i) << "\n";
	  
	  
	  // E_e1_thisstep.Rotate(theta_r-theta_i,rotate_Snell);
	  // Hhat_e1_thisstep.Rotate(theta_r-theta_i,rotate_Snell);
	  

	  // D_e1_thisstep=rotateE(epsilon_nextstep, angle_iceflow, E_e1_thisstep);  // D=epsilon*E
	  // TVector3 E_e1_nextstep;
	  // TVector3 khat_e1_nextstep=getNewkandE(nvec_nextstep,angle_iceflow,n_e1_nextstep,D_e1_thisstep,khat_ray1_2,
	  // 					E_e1_nextstep);



	  // Shat_e1_nextstep=1./E_e1_nextstep.Mag()*E_e1_nextstep.Cross(Hhat_e1_thisstep);
	  khat_ray1_2=khat_e1_nextstep;
	  radial_distance=sqrt(raypath[0][i][0]*raypath[0][i][0]+raypath[0][i][1]*raypath[0][i][1]);
//	  radial_distance+=stepray;
	  // radial_distance=(raypath[0][i]).Mag();
	} // end while radial_distance=5000
      } // emd if special depth
      

      // now rotate D_sum into D



      /*     

      getDeltaN(nvec,khat_ray2_1,angle_iceflow,
		n_e1,n_e2,p_e1,p_e2); // These are the eigenvectors of D_e1,0 and D_e2,0



      TVector3 D_e2=n_e2*n_e2*component_e2*p_e2; 









      cout << "e2:\n";
      //TVector3 khat_ray2_1=getNewkandE(nvec,angle_iceflow,n_e2,D_e2,
      //			     E_e2);
      
      TVector3 E=E_e1+E_e2;



      cout << "after first iteration, E is " << E[0] << "\t" << E[1] << "\t" << E[2] << "\n";
      cout << "remember Pt was " << Pt[i][0] << "\t" << Pt[i][1] << "\t" << Pt[i][2] << "\n";
      cout << "after first iteration, D is " << D[0] << "\t" << D[1] << "\t" << D[2] << "\n";
      cout << "remember D was (we took it to be Pt)" << Pt[i][0] << "\t" << Pt[i][1] << "\t" << Pt[i][2] << "\n"; 
      cout << "E_e1 is " << E_e1[0] << "\t" << E_e1[1] << "\t" << E_e1[2] << "\n";





      TVector3 testvec=-1.*D_e1+n_e1*n_e1*E_e1;
      cout << "mag of -D+n*n*E is " << testvec.Mag() << "\n";
      testvec.SetMag(1.);
      cout << "direction of -D+n*n*E is " << testvec[0] << "\t" << testvec[1] << "\t" << testvec[2] << "\n";

      cout << "E_e2 is " << E_e2[0] << "\t" << E_e2[1] << "\t" << E_e2[2] << "\n";
      cout << "D_e2 is " << D_e2[0] << "\t" << D_e2[1] << "\t" << D_e2[2] << "\n";




      E.SetMag(1.);

      TVector3 RotateE;
      RotateE=E.Cross(Pt[i]);
      RotateE.SetMag(1.);
      cout << "RotateE is " << RotateE[0] << "\t" << RotateE[1] << "\t" << RotateE[2] << "\n";
      double angle_rotate=acos(E.Dot(Pt[i]));
      cout << "angle_rotate is " << angle_rotate << "\n";
      E.Rotate(angle_rotate,RotateE);

      cout << "after rotation, E is " << E[0] << "\t" << E[1] << "\t" << E[2] << "\n";
      khat_ray1_1.Rotate(angle_rotate,RotateE);
      khat_ray2_1.Rotate(angle_rotate,RotateE);
      D.Rotate(angle_rotate,RotateE);

     cout << "after rotation, khat_ray1_1 is " << khat_ray1_1[0] << "\t" << khat_ray1_1[1] << "\t" << khat_ray1_1[2] << "\n";
     cout << "after rotation, khat_ray2_1 is " << khat_ray2_1[0] << "\t" << khat_ray2_1[1] << "\t" << khat_ray2_1[2] << "\n";
     cout << "recall that khat_0 is " << khat_0[0] << "\t" << khat_0[1] << "\t" << khat_0[2] << "\n";

     cout << "after rotation, D is " << D[0] << "\t" << D[1] << "\t" << D[2] << "\n";
      // rotate to Pt
      // use same rotation vector to rotate khats and p_e's
      
      
      getDeltaN(nvec,khat_ray1_1,angle_iceflow,
		n_e1,n_e2,p_e1,p_e2); // These are the eigenvectors of D_e1,0 and D_e2,0
      
      D.SetMag(1.);
      component_e1=p_e1.Dot(D);
      component_e2=p_e2.Dot(D);

      cout << "second iteration:\n";


      D_e1=n_e1*n_e1*component_e1*p_e1;

      cout << "e1:\n";
      TVector3 khat_ray1_2=getNewkandE(nvec,angle_iceflow,n_e1,D_e1,
				     E_e1);


      cout << "remember that after the first iteration khat_ray2_1 was " << khat_ray1_1[0] << "\t" << khat_ray1_1[1] << "\t" << khat_ray1_1[2] << "\n";

      D_e2=n_e2*n_e2*component_e2*p_e2;
      
      cout << "e2:\n";
      TVector3 khat_ray2_2=getNewkandE(nvec,angle_iceflow,n_e2,D_e2,
				     E_e2);

     cout << "khat_ray2_2 is " << khat_ray2_2[0] << "\t" << khat_ray2_2[1] << "\t" << khat_ray2_2[2] << "\n";
      cout << "remember that after the first iteration khat_ray2_1 was " << khat_ray2_1[0] << "\t" << khat_ray2_1[1] << "\t" << khat_ray2_1[2] << "\n";

      E=E_e1+E_e1;

      cout << "after second iteration, E is " << E[0] << "\t" << E[1] << "\t" << E[2] << "\n";
      cout << "remember Pt was " << Pt[i][0] << "\t" << Pt[i][1] << "\t" << Pt[i][2] << "\n";


      getDeltaN(nvec,khat_ray1_2,angle_iceflow,
		n_e1,n_e2,p_e1,p_e2); // These are the eigenvectors of D_e1,0 and D_e2,0
      
      cout << "second iteration:\n";

      D_e1=n_e1*n_e1*component_e1*p_e1;

      cout << "e1:\n";
      TVector3 khat_e1_3=getNewkandE(nvec,angle_iceflow,n_e1,D_e1,
				     E_e1);

      D_e2=n_e2*n_e2*component_e2*p_e2;      
      cout << "e2:\n";
      TVector3 khat_e2_3=getNewkandE(nvec,angle_iceflow,n_e2,D_e2,
				     E_e2);

      E=E_e1+E_e1;




      Pt[i]=E;
      rhat[i]=khat_ray1_2;
      */ // end iterative stuff
      /*

      temp = rx1_orientation.Cross(rhat_receive[i]);
      Pr1[i]=rhat_receive[i].Cross(temp);

      //      cout << "about to setmag for Pt and Pr1.\n";
      
      if (Pr1[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "Pr1[i] mag is " << Pr1[i].Mag() << "\n";
      Pr1[i].SetMag(1.);


      //double magnitude_of_pe=Pt[i].Dot(ordinary);
      //double magnitude_of_po=sqrt(1.0-magnitude_of_pe*magnitude_of_pe);
      //p_e[i]=magnitude_of_pe*ordinary;
      //p_o[i]=Pt[i]-p_e[i];

      //     cout << "dot product is " << rhat[i].Dot(ordinary) << "\n";
      //p_e[i]=Pt[i]-(rhat[i].Dot(ordinary))*ordinary; // component of Pt in ordinary direction
      //p_e[i]=(rhat[i].Dot(ordinary))*ordinary; // component of Pt in ordinary direction
      //p_o[i]=Pt[i]-p_e[i];

      //      temp=ordinary.Cross(Pt[i]);
      //TVector3 temp2=temp.Cross(ordinary);
      //temp2.SetMag(1.);
      //temp.SetMag(1.);
      //p_o[i]=(Pt[i].Dot(temp2))*temp2+(Pt[i].Dot(temp))*temp;

      //      p_o[i]=Pt[i]-p_e[i];

      //      cout << "about to setmag for p_o and p_e.\n";

      //      p_e[i]=Pt[i]-(rhat[i].Dot(ordinary))*ordinary; // component of Pt in ordinary direction                                                                                              
      //p_o[i]=Pt[i]-p_e[i];

      //      p_e[i]=Pt[i].Dot(ordinary)*ordinary;
      //p_o[i]=Pt[i]-p_e[i];

      //      p_o[i].SetMag(1.);
      //p_e[i].SetMag(1.);      

      if (WHICHPOL==1) {
	p_o[i]=p_e2;
	p_e[i]=p_e1;
      }
      else if (WHICHPOL==0) {
	p_o[i]=p_e1;
	p_e[i]=p_e2;
      }


      //      Pt_extraordinary[i]=Pt[i].Dot(p_e[i])*p_e[i];
      //Pt_ordinary[i]=Pt[i].Dot(p_o[i])*p_o[i];
      
      //      acrossp_o[i]=rhat_receive[i].Cross(p_o[i]);
      //acrossp_e[i]=rhat_receive[i].Cross(p_e[i]);

      //      cout << "about to setmag for acrossp_o and acrossp_e.\n";
      //acrossp_o[i].SetMag(1.);
      //acrossp_e[i].SetMag(1.);

      //      if (idepth==vdepth[i].size()-1) {


      Pr2[i] = rhat_receive[i].Cross(Pr1[i]);
      //cout << "last.\n";
      if (Pr2[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "Pr2[i] mag is " << Pr2[i].Mag() << "\n";
      Pr2[i].SetMag(1.);


      beam[0]=sin(acos(Pt[i].Dot(tx_orientation)))
	*sin(acos(rx1_orientation.Dot(rhat[i])));
      
      beam[1]=sin(acos(Pt[i].Dot(tx_orientation)))
	*sin(acos(rx1_orientation.Dot(rhat_receive[i])));



      //      if (vdepth[i][idepth]>-1001. && vdepth[i][idepth]<-998.) {
      //cout << "vdepth, idepth are " << vdepth[i][idepth] << "\t" << idepth << "\n";
	//}

      //      term1_Pr1[i].push_back(pow(Pt[i].Dot(p_e[i])*Pr1[i].Dot(p_e_rx[i]),2));
      //term2_Pr1[i].push_back(pow(Pt[i].Dot(p_o[i])*Pr1[i].Dot(p_o_rx[i]),2));

      term1_Pr1[i].push_back(Pr1[i].Dot(Ehat_khat1_e1_atrx)*field1);
      term2_Pr1[i].push_back(Pr1[i].Dot(Ehat_khat2_e1_atrx)*field2);

      term1_Pr2[i].push_back(Pr2[i].Dot(Ehat_khat1_e1_atrx)*field1);
      term2_Pr2[i].push_back(Pr2[i].Dot(Ehat_khat2_e1_atrx)*field2);

      //      term3_Pr1[i].push_back(Pt[i].Dot(p_e[i])*Pr1[i].Dot(p_e_rx[i])*Pt[i].Dot(p_o[i])*Pr1[i].Dot(p_o_rx[i]));
      //term3_Pr2[i].push_back(Pt[i].Dot(p_e[i])*Pr2[i].Dot(p_e_rx[i])*Pt[i].Dot(p_o[i])*Pr2[i].Dot(p_o_rx[i]));

      term3_Pr1[i].push_back(Pr1[i].Dot(Ehat_khat1_e1_atrx)*Pr1[i].Dot(Ehat_khat2_e1_atrx)*field1*field2);
      term3_Pr2[i].push_back(Pr2[i].Dot(Ehat_khat1_e1_atrx)*Pr2[i].Dot(Ehat_khat2_e1_atrx)*field1*field2);
    

      f1[i]->FixParameter(0,150.); // amplitude
      //f1[i]->FixParameter(1,term1_Pr1[i][term1_Pr1[i].size()-1]);
      //f1[i]->FixParameter(2,term2_Pr1[i][term2_Pr1[i].size()-1]);
      //f1[i]->FixParameter(3,term3_Pr1[i][term3_Pr1[i].size()-1]);
      //f1[i]->FixParameter(1,term1_Pr1[i][0]);
      //f1[i]->FixParameter(2,term2_Pr1[i][0]);
      //f1[i]->FixParameter(3,term3_Pr1[i][0]);
      if (WHICHPOL==1) {
	//f1[i]->FixParameter(0,300.); // amplitude
	f1[i]->FixParameter(1,term1_Pr2[i][term1_Pr2[i].size()-1]);
	f1[i]->FixParameter(2,term2_Pr2[i][term2_Pr2[i].size()-1]);
	f1[i]->FixParameter(3,term3_Pr2[i][term3_Pr2[i].size()-1]);
	//f1[i]->FixParameter(1,term1_Pr2[i][term1_Pr2[i].size()-1]);
	//f1[i]->FixParameter(2,term2_Pr2[i][term2_Pr2[i].size()-1]);
	//f1[i]->FixParameter(3,term3_Pr2[i][term3_Pr2[i].size()-1]);
	//f1[i]->FixParameter(3,0.5);
	//	cout << "fixing parameters for f1.  terms are " << term1_Pr2[i][term1_Pr2[i].size()-1] << "\t" << term2_Pr2[i][term1_Pr2[i].size()-1] << "\t" << term3_Pr2[i][term1_Pr2[i].size()-1] << "\n";
      }
      else if (WHICHPOL==0) {
	f1[i]->FixParameter(1,term1_Pr1[i][term1_Pr1[i].size()-1]);
	f1[i]->FixParameter(2,term2_Pr1[i][term2_Pr1[i].size()-1]);
	f1[i]->FixParameter(3,term3_Pr1[i][term3_Pr1[i].size()-1]);
	//f1[i]->FixParameter(1,term1_Pr1[i][term1_Pr1[i].size()-1]);
	//f1[i]->FixParameter(2,term2_Pr1[i][term2_Pr1[i].size()-1]);
	//f1[i]->FixParameter(3,term3_Pr1[i][term3_Pr1[i].size()-1]);
	//f1[i]->FixParameter(3,0.5);
	
      }
      //      f1[i]->FixParameter(4,freq);
      f1[i]->SetParameter(4,sumphase);

      //      f1[i]->SetParameter(5,deltan_exp[i]);
      //f1[i]->FixParameter(6,horizontal_distances[i]);
      //f1[i]->FixParameter(7,station_depths[i]);
      
      for (int jparam=0;jparam<3;jparam++) {
	f1_nointerference[i]->SetParameter(jparam,f1[i]->GetParameter(jparam));
      }

      //     cout << "vdepth, idepth are " << vdepth[i][idepth] << "\t" << idepth << "\n";
	//}
      // if (idepth==64) {
      // 	cout << "idepth=64 and parameters are " << idepth << "\t";
      // 	//	for (int iparam=0;iparam<8;iparam++) {
      // 	  for (int iparam=0;iparam<5;iparam++) {
      // 	  cout << f1[i]->GetParameter(iparam) << "\t";
      // 	}
      // 	cout << "\n";
      // cout << "idepth=64, rhat is " << rhat[i][0] << "\t" << rhat[i][1] << "\t" << rhat[i][2] << "\n";
      // cout << "idepth=64, p_e is " << p_e[i][0] << "\t" << p_e[i][1] << "\t" << p_e[i][2] << "\n";
      // cout << "idepth=64, p_o is " << p_o[i][0] << "\t" << p_o[i][1] << "\t" << p_o[i][2] << "\n";
      // cout << "idepth=64, p_e_rxatt is " << p_e_rx[i][0] << "\t" << p_e_rx[i][1] << "\t" << p_e_rx[i][2] << "\n";
      // cout << "idepth=64, p_o_rx is " << p_o_rx[i][0] << "\t" << p_o_rx[i][1] << "\t" << p_o_rx[i][2] << "\n";
      // cout << "idepth=64, Pt is " << Pt[i][0] << "\t" << Pt[i][1] << "\t" << Pt[i][2] << "\n";
      // cout << "idepth=64, function evaluated at 1000 depth is " << f1[i]->Eval(1000.) << "\n";
      // cout << "idepth=64, sumphase at 1000 depth is " << sumphase << "\n";

      // }
    
  


	  // get_attenuation_along_path(pospulser_special, posstation_special, solutions_special[0][1],
	  //			     freq, N_ICE, DELTA_N, Z_0, 1); // using model=1, should ask Uzair what that is.
      //cout << "got attenuation.\n";

      //      vmag_atten[i][idepth]=f1[i]->GetParameter(0)*exp(-1.*vtotal_distances[i][idepth]/L_ATTEN)/exp(-1.*vtotal_distances[0][0]/L_ATTEN);  

      // if (i==5) {

      // cout << "i, idepth are " << i << "\t" << idepth << "\t" << vmag_atten.size() << "\t" << vmag_atten[i].size() << "\t" << vmag_atten_beam[i].size() << "\n";
      // //      cout << "i, idepth, sizes are " << i << "\t" << idepth << "\t" << vmag_atten.size() << "\t" << vmag_atten[i].size() << "\t" << vmag_atten_beam[i].size() << "\t" << vmag_atten_beam_crosspol[i].size() << "\t" << f1[5]->GetParameter(0) << "\t" << beam[WHICHPOL] << "\n";
      // cout << "i, atten, beam are " << i << "\t" << atten << "\t" << beam[WHICHPOL] << "\t" << f1[i]->GetParameter(0) << "\t" << f1_nointerference[i]->Eval(vdepth[i][idepth]) << "\n";
      // cout << "parameters of f1, f1_nointerference are \n";
      // cout << f1[i]->GetParameter(0) << "\t" << f1_nointerference[i]->GetParameter(0) << "\n";
      // cout << f1[i]->GetParameter(1) << "\t" << f1_nointerference[i]->GetParameter(1) << "\n";
      // cout << f1[i]->GetParameter(2) << "\t" << f1_nointerference[i]->GetParameter(2) << "\n";

      // }
    
      vmag_parameter0[i][idepth]=f1[i]->GetParameter(0);
      vmag_atten[i][idepth]=f1[i]->GetParameter(0)*atten;
      //cout << "made it here 1.\n";
      vmag_atten_beam[i][idepth]=f1[i]->GetParameter(0)*atten*beam[WHICHPOL];
      //cout << "made it here 2.\n";
      vmag_atten_beam_crosspol[i][idepth]=f1[i]->GetParameter(0)*atten*beam[WHICHPOL]; // gonna have to read in both polarizatoins no matter which polarization we are interested in.
      //vmag_func[i][idepth]=f1[i]->Eval(vdepth[i][idepth]);
      //cout << "made it here 3.\n";

      //      cout << "f1[i]->GetParameter(0), vtotal_distances, vtotal_distances are " << f1[i]->GetParameter(0) << "\t" << vtotal_distances[i][idepth] << exp(-1.*vtotal_distances[0][0]/L_ATTEN) << "\n";

      //      vmag_atten_func[i][idepth]=f1[i]->Eval(vdepth[i][idepth])*exp(-1.*vtotal_distances[i][idepth]/L_ATTEN)/exp(-1.*vtotal_distances[0][0]/L_ATTEN);
      vmag_atten_beam_crosspol_nointerferencefunc[i][idepth]=f1_nointerference[i]->Eval(vdepth[i][idepth])*atten*beam[WHICHPOL];
      vmag_atten_beam_crosspol_func[i][idepth]=f1[i]->Eval(vdepth[i][idepth])*atten*beam[WHICHPOL];

      //      for (int ispecialdepths=0;ispecialdepths<NSPECIALDEPTHS;ispecialdepths++) {



      


      //      if (idepth==vdepth[i].size()-1) {
      if (idepth==0) {
	// //cout << "i is " << i << "\n";
	// cout << "ordinary is " << ordinary.X() << "\t" << ordinary.Y() << "\t" << ordinary.Z() << "\n";
	// cout << "A" << i+1 << ": r is " << rhat[i].X() << "\t" << rhat[i].Y() << "\t" << rhat[i].Z() << "\n";
	// cout << "A" << i+1 << ": p_o is " << p_o[i].X() << "\t" << p_o[i].Y() << "\t" << p_o[i].Z() << "\n";
	// cout << "A" << i+1 << ": p_e is " << p_e[i].X() << "\t" << p_e[i].Y() << "\t" << p_e[i].Z() << "\n";
	// //cout << "A" << i+1 << ": acrossp_o is " << acrossp_o[i].X() << "\t" << acrossp_o[i].Y() << "\t" << acrossp_o[i].Z() << "\n";
	// //cout << "A" << i+1 << ": acrossp_e is " << acrossp_e[i].X() << "\t" << acrossp_e[i].Y() << "\t" << acrossp_e[i].Z() << "\n";
	// cout << "A" << i+1 << ": Pt is " << Pt[i].X() << "\t" << Pt[i].Y() << "\t" << Pt[i].Z() << "\n";
	// cout << "A" << i+1 << ": Pt*p_e and Pt*p_o are " << Pt[i].Dot(p_e[i]) << "\t" << Pt[i].Dot(p_o[i]) << "\n";
	// cout << "A" << i+1 << ": Pt component of extraordinary that is HPol is " << Pt_extraordinary[i].Dot(Pr2[i]) << "\n";
	// cout << "A" << i+1 << ": Pt component of ordinary that is HPol is " << Pt_ordinary[i].Dot(Pr2[i]) << "\n";
	// cout << "A" << i+1 << ": Pt component of extraordinary that is VPol is " << Pt_extraordinary[i].Dot(Pr1[i]) << "\n";
	// cout << "A" << i+1 << ": Pt component of ordinary that is VPol is " << Pt_ordinary[i].Dot(Pr1[i]) << "\n";
      
	// cout << "A" << i+1 << ": Sum of projections of Pt is " << pow(Pt[i].Dot(p_e[i]),2) +pow(Pt[i].Dot(p_o[i]),2) << "\n";
	// //	cout << "A" << i+1 << ": Pt*acrossp_e and Pt*acrossp_o are " << Pt[i].Dot(acrossp_e[i]) << "\t" << Pt[i].Dot(acrossp_o[i]) << "\n";
	// //cout << "A" << i+1 << ": Sum of projections of Pt is " << pow(Pt[i].Dot(acrossp_e[i]),2) +pow(Pt[i].Dot(acrossp_o[i]),2) << "\n";
	
	// cout << "A" << i+1 << ": Pr1 is " << Pr1[i].X() << "\t" << Pr1[i].Y() << "\t" << Pr1[i].Z() << "\n";
	// cout << "A" << i+1 << ": Pr1*p_e and Pr1*p_o are " << Pr1[i].Dot(p_e[i]) << "\t" << Pr1[i].Dot(p_o[i]) << "\n";
	// cout << "A" << i+1 << ": Sum of projections of Pr1 is " << pow(Pr1[i].Dot(p_e[i]),2) +pow(Pr1[i].Dot(p_o[i]),2) << "\n";
	// //	cout << "A" << i+1 << ": Pr1*acrossp_e and Pr1*acrossp_o are " << Pr1[i].Dot(acrossp_e[i]) << "\t" << Pr1[i].Dot(acrossp_o[i]) << "\n";
	// //	cout << "A" << i+1 << ": Sum of projections of Pr1 is " << pow(Pr1[i].Dot(acrossp_e[i]),2) +pow(Pr1[i].Dot(acrossp_o[i]),2) << "\n";

	// cout << "A" << i+1 << ": Pr2 is " << Pr2[i].X() << "\t" << Pr2[i].Y() << "\t" << Pr2[i].Z() << "\n";
	// cout << "A" << i+1 << ": Pr2*p_e and Pr2*p_o are " << Pr2[i].Dot(p_e[i]) << "\t" << Pr2[i].Dot(p_o[i]) << "\n";
	// cout << "A" << i+1 << ": Sum of projections of Pr2 is " << pow(Pr2[i].Dot(p_e[i]),2) +pow(Pr2[i].Dot(p_o[i]),2) << "\n";
	// //	cout << "A" << i+1 << ": Pr2*acrossp_e and Pr2*acrossp_o are " << Pr2[i].Dot(acrossp_e[i]) << "\t" << Pr2[i].Dot(acrossp_o[i]) << "\n";
	// //cout << "A" << i+1 << ": Sum of projections of Pr2 is " << pow(Pr2[i].Dot(acrossp_e[i]),2) +pow(Pr2[i].Dot(acrossp_o[i]),2) << "\n";

	// cout << "A" << i+1 << ": For station " << i+1 << ", magnitude of oscillation in V is " << Pt[i].Dot(p_e[i])*Pr1[i].Dot(p_e[i])*Pt[i].Dot(p_o[i])*Pr1[i].Dot(p_o[i]) << "\n";
	// cout << "A" << i+1 << ": For station " << i+1 << ", magnitude of oscillation in H is " << Pt[i].Dot(p_e[i])*Pr2[i].Dot(p_e[i])*Pt[i].Dot(p_o[i])*Pr1[i].Dot(p_o[i]) << "\n";
	// cout << "A" << i+1 << ": Pr1 terms are " << term1_Pr1[i][idepth] << "\t" << term2_Pr1[i][idepth] << "\t" << term3_Pr1[i][idepth] << "\n";
	// cout << "A" << i+1 << ": Pr2 terms are " << term1_Pr2[i][idepth] << "\t" << term2_Pr2[i][idepth] << "\t" << term3_Pr2[i][idepth] << "\n";
      }
      */

    } // end loop over depths
      

    //    graypath_z_x[i]=new TGraph (vraypos[i][0].size(),&vraypos[i][0][0],&vraypos[i][2][0]);
    //graypath_z_x[i]->SetLineColor(icolors[i]);
    //graypath_z_y[i]=new TGraph (vraypos[i][1].size(),&vraypos[i][1][0],&vraypos[i][2][0]);
    //graypath_z_y[i]->SetLineColor(icolors[i]);
    //graypath_y_x[i]=new TGraph (vraypos[i][0].size(),&vraypos[i][0][0],&vraypos[i][1][0]);
    //graypath_y_x[i]->SetLineColor(icolors[i]);


    //    cout << "about to make graphs.\n";
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

    gtxdepthE_theta1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_theta1[i][0]);
    gtxdepthE_theta2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepth_theta2[i][0]);

    gtxdepthE_theta1_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepthE_theta1_Sclock[i][0]);
    gtxdepthE_theta2_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vtxdepthE_theta2_Sclock[i][0]);

    grxdepthE_theta1_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepthE_theta1_Sclock[i][0]);
    grxdepthE_theta2_Sclock[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepthE_theta2_Sclock[i][0]);
    
    gdotShats_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vdotShats_tx[i][0]);
    gdotEhats_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vdotEhats_tx[i][0]);
    gdotDhats_tx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vdotDhats_tx[i][0]);



    //    cout << "sizes are " << vdepth[i].size() << "\t" << vmag_atten[i].size() << "\n";
    //cout << "vdepths are \n";
    //for (int itest=0;itest<vdepth[i].size();itest++) {
    //cout << vdepth[i][itest] << "\n";
    //}
    //cout << "i is " << i << ", starting.\n";
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
    for (int idepth=0;idepth<vdepth[i].size();idepth++) {
      g_spectra[i][idepth]=new TGraph(vfreqs.size(),&vfreqs[0],&vspectra[i][idepth][0]);
    }
    g_atten_power[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_atten_power[i][0]);
    g_atten_beam_power[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vrxdepth_atten_beam_power[i][0]);

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

    // if (i==5) {
    //   cout << "index at -1000m is " << g_idepth[i]->Eval(-1000.) << "\n";
    //   cout << "eval gepsilon1_tx at -1000m is " << gepsilon1_tx[i]->Eval(-1000.) << "\n";
    //   cout << "eval gepsilon2_tx at -1000m is " << gepsilon2_tx[i]->Eval(-1000.) << "\n";
    // }
    
    gpolarization_Omega_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vpolarization_Omega_rx[i][0]);
    gpolarization_Psi_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vpolarization_Psi_rx[i][0]);

    gEpolarization_Omega_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vEpolarization_Omega_rx[i][0]);
    gEpolarization_Psi_rx[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vEpolarization_Psi_rx[i][0]);

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
    //cout << "i is " << i << ", got half way through.\n";
    //g_path[i]=new TGraph(vres[i].size(),&vres[i].size(),&vzs[i].size());
    g_depth_istep[i]=new TGraph(vdepth_step[i].size(),&vdepth_step[i][0],&vistep[i][0]);

    g_attenlengths[i]=new TGraph(vlengths[i].size(),&vlengths[i][0],&vattenlengths[i][0]);
    g_receive_launch[i]=new TGraph(vreceiveangle[i].size(),&vreceiveangle[i][0],&vlaunchangle[i][0]);
    g_receive[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vreceiveangle[i][0]);
    g_launch[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vlaunchangle[i][0]);
    g_output6[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&voutput6[i][0]);
    g_output7[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&voutput7[i][0]);
    g_output8[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&voutput8[i][0]);
    //cout << "i is " << i << ", got all the way through.\n";
  } // loop over stations


  string filename;






  TGraph *g1[NSTATIONS];

  TCanvas *c1=new TCanvas("c1","c1",800,800);
  c1->Divide(2,3);
  TH2D *h2[NSTATIONS];
    

  double vmax[2][NSTATIONS]={{100.,40.,40.,40.,40.,100.},
  //double vmax[2][NSTATIONS]={{400.,400.,400.,400.,100.,25.},
  //double vmax[2][NSTATIONS]={{10.,10.,10.,10.,1.,1.},
			       //		       {100.,100.,100.,100.,100.,100.}};
		       {100.,40.,40.,40.,40.,100.}};
    //double vmax[6]={500.,1500.,1000.,1500.,4000.,10.};
   
  double pmax[2][NSTATIONS]={{1.E5,1.E5,1.E5,1.E5,1.E5,1.E5},
			     {1.E5,1.E5,1.E5,1.E5,1.E5,1.E5}};


  //  for (int i=0;i<NSTATIONS;i++) {
  for (int i=minstation;i<=maxstation;i++) {

    h2[i]=new TH2D("h2","h2",200,-1500.,-850,100,0.,vmax[WHICHPOL][i]);
    h2[i]->SetXTitle("SpiceCore pulser height (m)");
    h2[i]->SetYTitle("V_{SNR}^{max}");
  }
    


  //  for (int i=0;i<5;i++) {
  //  for (int i=0;i<NSTATIONS;i++) {
  for (int i=minstation;i<=maxstation;i++) {
    //if (i==0 || i==3) {
  //    if (i==0 || i==1) {
    c1->cd(i+1);
   
    h2[i]->Draw();
    //f1[i]=new TF1("f1","[3]*abs(cos(3.14159*[0]*[1]/300000000.*sqrt([2]*[2]+x*x)))",-1500.,-600.);
    // this was circular polarization
    //    f1[i]=new TF1("f1","[3]*abs(sin(3.14159*[0]*[1]/300000000.*sqrt([2]*[2]+(x-[5])*(x-[5]))))",-1500.,-600.);
    //f1[i]->SetLineColor(icolors[i]);
    //    f1[i]->SetLineColor(kBlack);
    //f1[i]->FixParameter(0,freq);    
    //f1[i]->SetParameter(1,deltan_exp[i]);
    //f1[i]->SetParameter(1,deltan_exp[0]);
    //f1[i]->FixParameter(2,horizontal_distances[i]);
    //f1[i]->SetParameter(3,100.);
    //f1[i]->FixParameter(4,L_ATTEN);
    //f1[i]->FixParameter(5,station_depths[i]);

    //    cout << "A" << i+1 << ": Pr2 terms are " << f1[i]->GetParameter(1) << "\t" << f1[i]->GetParameter(2) << "\t" << f1[i]->GetParameter(3) << "\n";
    //cout << "A" << i+1 << ": length of vectors are " << vdepth[i].size() << "\t" << vsnrmax[i].size() << "\n";
    //   cout << "A" << i+1 << ": factor in front are " << TMath::Pi()*f1[i]->GetParameter(4)*f1[i]->GetParameter(5)/TMath::C() << "\n";
    //cout << "A" << i+1 << ": pi / factor in front are " << TMath::Pi()/(TMath::Pi()*f1[i]->GetParameter(4)*f1[i]->GetParameter(5)/TMath::C()) << "\n";
    //    cout << "i, factor in front are " << TMath::Pi()*f1[i]->GetParameter(4)*deltan_exp[i]/TMath::C() << "\n";
    //cout << "i, pi / factor in front are " << TMath::Pi()/(TMath::Pi()*f1[i]->GetParameter(4)*2.79/TMath::C()) << "\n";
    //g1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vsnrmax[i][0]);
    if (i!=5) {
      g1[i]=new TGraphErrors(vdepth_data[i].size(),&vdepth_data[i][0],&vsnrmax[i][0],&vdepth_data_err[i][0],&vsnrmax_err[i][0]);
      g1[i]->SetMarkerColor(icolors[i]);
      g1[i]->SetMarkerSize(1.);
      g1[i]->SetMarkerStyle(20);
      g1[i]->Draw("pesame");
    }
    //g1[i]->Fit(f1[i]);
    //    deltan_obs[i]=f1[i]->GetParameter(5);
    //cout << "A" << i+1 << ": alpha is " << alpha_deg[i] << "\n";
    //cout << "A" << i+1 << ": expected deltan is " << deltan_exp[i] << "\n";
    //cout << "A" << i+1 << ": observed deltan is " << deltan_obs[i] << "\n";
    deltan_obs_err[i]=f1[i]->GetParError(1);
    //cout << "A" << i+1 << ": observed deltan error is " << deltan_obs_err[i] << "\n";
 

    //    f1[i]->Draw("same");


    g_parameter0[i]->SetLineWidth(3);
    g_parameter0[i]->SetLineColor(kGreen);
    //g_parameter0[i]->Draw("same");

    g_atten[i]->SetLineWidth(3);
    g_atten[i]->SetLineColor(kOrange);
    //g_atten[i]->Draw("same");

    g_atten_beam[i]->SetLineWidth(3);
    g_atten_beam[i]->SetLineStyle(kDashed);
    g_atten_beam[i]->SetLineColor(kOrange);
    //g_atten_beam[i]->Draw("same");

    g_atten_beam_crosspol[i]->SetLineWidth(3);
    g_atten_beam_crosspol[i]->SetLineStyle(kDashed);
    g_atten_beam_crosspol[i]->SetLineColor(kOrange);
    //g_atten_beam_crosspol[i]->Draw("same");

    gfunc_noadjust[i]->SetLineWidth(3);
    gfunc_noadjust[i]->SetLineColor(kBlack);

    g_atten_beam_crosspol_nointerferencefunc[i]->SetLineWidth(3);
    g_atten_beam_crosspol_nointerferencefunc[i]->SetLineStyle(kDashed);
    g_atten_beam_crosspol_nointerferencefunc[i]->SetLineColor(kBlack);
    //g_atten_beam_crosspol_nointerferencefunc[i]->Draw("same");

    g_atten_beam_crosspol_func[i]->SetLineWidth(3);
    g_atten_beam_crosspol_func[i]->SetLineColor(kBlack);
    //g_atten_beam_crosspol_func[i]->Draw("same");


    //gfunc_noadjust[i]->SetLineWidth(3);
    //gfunc_noadjust[i]->SetLineColor(kOrange);
    
    //gfunc_noadjust[i]->Draw("same");
    //    gfunc_adjust_latten[i]->Draw("same");

    //    f1[i]->SetParameter(1,deltan_exp[i]);
    //f1[i]->SetLineColor(kBlack);
    //f1[i]->SetLineStyle(2);
    //f1[i]->Draw("same");

    //    f1_distances[i]=new TF1("f1","[0]*sqrt([1]+[2]-2*[3]*cos(TMath::Pi()*[4]*[5]/TMath::C()*x))",0.,5000.);
    //    f1_distances[i]=new TF1("f1","[0]*sqrt([1]+[2]-2*[3]*cos(TMath::Pi()*[4]*[5]/TMath::C()*x))",0.,5000.);
    //f1_distances[i]=new TF1("f1","[0]*([1]+[2]-2.*[3]*cos(TMath::Pi()*[4]*[5]/TMath::C()*x))",0.,5000.);
    //f1_distances[i]=new TF1("f1","[0]*(2.*[3]*sqrt(abs(cos(TMath::Pi()*[4]*[5]/TMath::C()*x))))",0.,5000.);
    f1_distances[i]=new TF1("f1",sfunc_distances.c_str(),0.,5000.);
    f1_nointerference_distances[i]=new TF1("f1_nointerference",sfunc_nointerference_distances.c_str(),0.,5000.);

    //    f1_distances[i]=new TF1("f1","[3]*abs(sin(3.14159*[0]*[1]/300000000.*x))",0.,5000.);
    f1_distances[i]->SetLineColor(kBlack);
  


    for (int j=0;j<NDISTANCES_BIGPIC-(jmin[i]+1);j++) {

      //      if (j==117)
      //vpseudodepths_bigpic[i][j]=-1000.;

      double posstation[2];
      posstation[0]=sqrt(pow(station_coords[i][0]-pulser_coords[0],2)+pow(station_coords[i][1]-pulser_coords[1],2));
      posstation[1]=station_depths[i];

      double pospulser[2];
      pospulser[0]=0.;
      pospulser[1]=vpseudodepths_bigpic[i][j];

      //      vector< vector<double> > solutions = find_solutions(pospulser,posstation,
      //						  N_ICE,DELTA_N,Z_0);


      //      vnsolutions_bigpic[i].push_back((double)solutions.size());

      double x0=0;
      double z0=pospulser[1];
      double x1=posstation[0];
      double z1=posstation[1];

      //      double *params_bigpic=GetDirectRayPar(z0,x1,z1);

      //    double *getresults_bigpic=IceRayTracing(x0,z0,x1,z1);

     double *getresults=IceRayTracing(x0,z0,x1,z1);

      double lvalue;
      vector<double> res; 	//res are y coordinates
      vector<double> zs;	//zs are z coordinates
      double launch_angle;
      double receive_angle;
      double *paramsd;
      double *paramsra;

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
	//	receive_angle=(180.-paramsra[0])/DEGRAD;
	receive_angle=paramsra[0]/DEGRAD;
	GetFullRefractedRayPath(z0, x1, z1, paramsra[7], paramsra[3], res, zs);
	// need to get fullrefractedraypath
      }


      double atten=1.;
      double sumphase=0;
      //      if (solutions.size()>0) {
      if (getresults[6]!=-1000 || getresults[8]!=-1000) {
	//cout << "attenuation length is " << get_attenuation_length(-200.,0.3,1) << "\n";
	//atten=get_attenuation_along_path(pospulser, posstation, solutions[0][1],
	//freq/1000., N_ICE, DELTA_N, Z_0, 1);
      //cout << "atten is " << atten << "\n";
	//double path_length=get_path_length(pospulser, posstation, solutions[0][1], N_ICE, DELTA_N, Z_0);
	//cout << "path_length is " << path_length << "\n";
	


	//	get_path(N_ICE, DELTA_N, Z_0, pospulser, posstation, solutions[0][1], res, zs, 200);
	  // the "y" direction is the horizontal direction from the pulser to the station.
	
	//	GetFullDirectRayPath(z0, x1, z1, params_bigpic[3], res, zs);


	TVector3 yhat(station_coords[i][0]-pulser_coords[0],
		      station_coords[i][1]-pulser_coords[1],
		      0.);
	if (yhat.Mag()<HOWSMALLISTOOSMALL)
	  cout << "yhat mag is " << yhat.Mag() << "\n";
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
	    
	    //	    double atten_length=get_attenuation_length(zs[istep],freq/1.E9,1);
	    
	    double atten_length=GetIceAttenuationLength(zs[istep], freq/1.E9);

	    //	    cout << "zs, length, atten_length, atten are " << zs[istep] << "\t" << length << "\t" << atten_length << "\t" << atten << "\n";
	    
	    atten*=exp(-1.*length/atten_length);
	   

	    double deltan_alongpath=getDeltaN(BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1,p_e2);
	   

	    //	    cout << "deltan_alongpath 2 is " << deltan_alongpath << "\n";
	    // need to add the stuff about notflipped
	    //   cout << "fine:  deltan_alongpath, one period is " << deltan_alongpath << "\t" << 1./(deltan_alongpath*PI/TMath::C()*freq) << "\n";
	    

	    

	  }
	  else {
	    rhat_thisstep=rhat[i];
	    
	  }


	} // end stepping over the ray
	// inside loop over solutions
	// inside loop over distances_bigpic
	vtimediff_bigpic[i].push_back(sumphase);
    
	//	cout << "pseudodepth, sumphase are " << vpseudodepths_bigpic[i][j] << "\t" << sumphase << "\n";	
	//	double launch_angle=get_angle(posstation,pospulser, solutions[0][1], N_ICE, DELTA_N, Z_0);
	//double receive_angle=get_receive_angle(pospulser,posstation, solutions[0][1], N_ICE, DELTA_N, Z_0);

	//	double *getangles_bigpic=GetDirectRayPar(z0, x1, z1);

  //	double launch_angle=get_angle(posstation,pospulser, solutions[0][1], N_ICE, DELTA_N, Z_0);
  //	double receive_angle=get_receive_angle(pospulser,posstation, solutions[0][1], N_ICE, DELTA_N, Z_0);
      




	//	cout << "A" << i << ": posstation is " << posstation[0] << "\t" << posstation[1] << "\n";
	//cout << "A" << i << ": pospulser is " << pospulser[0] << "\t" << pospulser[1] << "\n";
	//cout << "A" << i << ": launch_angle, receive_angle are " << launch_angle << "\t" << receive_angle << "\n";
      
	if (DORAYTRACING) {

	  TVector3 plusz(0.,0.,1.);
	  //	  TVector3 minusz(0.,0.,-1.);
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

 	rhat[i].SetX(station_coords[i][0]-pulser_coords[0]); // this will be a negative number
  	rhat[i].SetY(station_coords[i][1]-pulser_coords[1]);
  	rhat[i].SetZ(station_depths[i]-vpseudodepths_bigpic[i][j]);
	//rhat[i].SetZ(station_depths[i]-vdepth[i][0]); // this is a difference between two negative numbers.

	//}

	rhat_receive[i].SetX(rhat[i][0]);
	rhat_receive[i].SetY(rhat[i][1]);
	rhat_receive[i].SetZ(rhat[i][2]);
      //rhat_receive[i].SetZ(0.);



	}

      }
      else {




      //      for (int j=0;j<3;j++) {
      rhat[i].SetX(station_coords[i][0]-pulser_coords[0]); // this will be a negative number
      rhat[i].SetY(station_coords[i][1]-pulser_coords[1]);
      rhat[i].SetZ(station_depths[i]-vpseudodepths_bigpic[i][j]);
      //cout << "station_depths[i], vpseudodepths_bigpic[i][j] are " << station_depths[i] << "\t" << vpseudodepths_bigpic[i][j] << "\n";
      rhat_receive[i].SetX(rhat[i][0]);
      rhat_receive[i].SetY(rhat[i][1]);
      rhat_receive[i].SetZ(rhat[i][2]);
      //rhat_receive[i].SetZ(0.);
	//}
      }
      if (rhat[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "rhat[i] mag is " << rhat[i].Mag() << "\n";
      rhat[i].SetMag(1.);
      if (rhat_receive[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "rhat_receive[i] mag is " << rhat_receive[i].Mag() << "\n";

      rhat_receive[i].SetMag(1.);
      /*
      temp = tx_orientation.Cross(rhat[i]);
      Pt[i]=rhat[i].Cross(temp);

      temp = rx1_orientation.Cross(rhat_receive[i]);
      Pr1[i]=rhat_receive[i].Cross(temp);

      //      cout << "about to setmag for Pt and Pr1.\n";
      if (Pt[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "Pt[i] mag is " << Pt[i].Mag() << "\n";

      Pt[i].SetMag(1.);
      
      if (Pr1[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "Pr1[i] mag is " << Pr1[i].Mag() << "\n";
      Pr1[i].SetMag(1.);

      Pr2[i] = rhat_receive[i].Cross(Pr1[i]);
      //cout << "last.\n";
      if (Pr2[i].Mag()<HOWSMALLISTOOSMALL)
	cout << "Pr2[i] mag is " << Pr2[i].Mag() << "\n";
      Pr2[i].SetMag(1.);


	beam[0]=sin(acos(Pt[i].Dot(tx_orientation)))
	  *sin(acos(rx1_orientation.Dot(rhat[i])));

	beam[1]=sin(acos(Pt[i].Dot(tx_orientation)))
	  *sin(acos(rx1_orientation.Dot(rhat_receive[i])));



      nvec[0]=gn1->Eval(station_depths[i]);
      nvec[1]=gn2->Eval(station_depths[i]);
      nvec[2]=gn3->Eval(station_depths[i]);

      double deltan_attherx=getDeltaN(nvec,rhat[i],angle_iceflow,n_e1,n_e2,p_e1,p_e2);

      if (WHICHPOL==1) {
	p_o_rx[i]=p_e2;
	p_e_rx[i]=p_e1;
      }
      else if (WHICHPOL==0) {
	p_o_rx[i]=p_e1;
	p_e_rx[i]=p_e2;
      }

      nvec[0]=gn1->Eval(vpseudodepths_bigpic[i][j]);
      nvec[1]=gn2->Eval(vpseudodepths_bigpic[i][j]);
      nvec[2]=gn3->Eval(vpseudodepths_bigpic[i][j]);


      //      cout << "vpseudodepths, nvecs are " << vpseudodepths_bigpic[i][j] << "\t" << nvec[0] << "\t" << nvec[1] << "\t" << nvec[2] << "\n";

      double deltan_atthetx=getDeltaN(nvec,rhat[i],angle_iceflow,n_e1,n_e2,p_e1,p_e2);



      if (WHICHPOL==1) {
	p_o[i]=p_e2;
	p_e[i]=p_e1;
      }
      else if (WHICHPOL==0) {
	p_o[i]=p_e1;
	p_e[i]=p_e2;
      }


      term1_Pr1_distances[i].push_back(pow(Pt[i].Dot(p_e[i])*Pr1[i].Dot(p_e_rx[i]),2));
      term2_Pr1_distances[i].push_back(pow(Pt[i].Dot(p_o[i])*Pr1[i].Dot(p_o_rx[i]),2));

      term1_Pr2_distances[i].push_back(pow(Pt[i].Dot(p_e[i])*Pr2[i].Dot(p_e_rx[i]),2));
      term2_Pr2_distances[i].push_back(pow(Pt[i].Dot(p_o[i])*Pr2[i].Dot(p_o_rx[i]),2));

      term3_Pr1_distances[i].push_back(Pt[i].Dot(p_e[i])*Pr1[i].Dot(p_e_rx[i])*Pt[i].Dot(p_o[i])*Pr1[i].Dot(p_o_rx[i]));
      term3_Pr2_distances[i].push_back(Pt[i].Dot(p_e[i])*Pr2[i].Dot(p_e_rx[i])*Pt[i].Dot(p_o[i])*Pr2[i].Dot(p_o_rx[i]));



      double averagedepth=(station_depths[i]+vpseudodepths_bigpic[i][j])/2.;
      nvec[0]=gn1->Eval(averagedepth);
      nvec[1]=gn2->Eval(averagedepth);
      nvec[2]=gn3->Eval(averagedepth);

      double deltan=getDeltaN(nvec,rhat[i],angle_iceflow,n_e1,n_e2,p_e1,p_e2);

      deltan_exp[i]=deltan;

      //      cout << "For _distances plot, A" << i+1 << ": expected deltan is " << deltan_exp[i] << "\n";
      //      if (vpseudodepths_bigpic[i][j]>-1010. && vpseudodepths_bigpic[i][j]<-990.) {
      //cout << "Depth, j are " << vpseudodepths_bigpic[i][j] << "\t" << j << "\n";
	//}

    if (WHICHPOL==0) {
      //      f1_distances[i]->SetParameter(0,1000.);
      f1_distances[i]->FixParameter(0,f1[i]->GetParameter(0));
      //      f1_distances[i]->FixParameter(1,f1[i]->GetParameter(1));
      //f1_distances[i]->FixParameter(2,f1[i]->GetParameter(2));
      //f1_distances[i]->FixParameter(3,f1[i]->GetParameter(3));
      f1_distances[i]->FixParameter(1,term1_Pr1_distances[i][term1_Pr1_distances[i].size()-1]);
      f1_distances[i]->FixParameter(2,term2_Pr1_distances[i][term2_Pr1_distances[i].size()-1]);
      f1_distances[i]->FixParameter(3,term3_Pr1_distances[i][term3_Pr1_distances[i].size()-1]);
      //f1_distances[i]->FixParameter(3,0.2);
      //cout << "fixing parameters for f1.  terms are " << f1_distances[i]->GetParameter(1) << "\t" << f1_distances[i]->GetParameter(2) << "\t" << f1_distances[i]->GetParameter(3) << "\n";
      //      f1_distances[i]->FixParameter(4,freq);
      f1_distances[i]->FixParameter(4,sumphase);
      //      f1_distances[i]->SetParameter(5,deltan_exp[i]);
      //f1_distances[i]->FixParameter(5,deltan_exp[i]);
    }
    else if (WHICHPOL==1) {
      //f1_distances[i]->SetParameter(0,150.);
      //      f1_distances[i]->FixParameter(0,150.);
      //      f1_distances[i]->SetParameter(0,f1[i]->GetParameter(0));
      //      f1_distances[i]->FixParameter(1,f1[i]->GetParameter(1));
      //f1_distances[i]->FixParameter(2,f1[i]->GetParameter(2));
      //f1_distances[i]->FixParameter(3,f1[i]->GetParameter(3));
      f1_distances[i]->FixParameter(1,term1_Pr2_distances[i][term1_Pr2_distances[i].size()-1]);
      f1_distances[i]->FixParameter(2,term2_Pr2_distances[i][term2_Pr2_distances[i].size()-1]);
      f1_distances[i]->FixParameter(3,term3_Pr2_distances[i][term3_Pr2_distances[i].size()-1]);
      //f1_distances[i]->FixParameter(3,0.2);
      //      f1_distances[i]->FixParameter(4,freq);
      f1_distances[i]->FixParameter(4,sumphase);
      //      f1_distances[i]->SetParameter(5,deltan_exp[i]);
      //f1_distances[i]->FixParameter(5,deltan_exp[i]);
    }
      for (int jparam=0;jparam<3;jparam++) {
	f1_nointerference_distances[i]->SetParameter(jparam,f1_distances[i]->GetParameter(jparam));
      }
    
    if (j==117) {
      cout << "j=117 and parameters are " << j << "\t";
      //      for (int iparam=0;iparam<6;iparam++) {
      for (int iparam=0;iparam<5;iparam++) {
	cout << f1_distances[i]->GetParameter(iparam) << "\t";
      }
      cout << "\n";
      cout << "j=117 rhat is " << rhat[i][0] << "\t" << rhat[i][1] << "\t" << rhat[i][2] << "\n";
      cout << "j=117 p_e is " << p_e[i][0] << "\t" << p_e[i][1] << "\t" << p_e[i][2] << "\n";
      cout << "j=117 p_o is " << p_o[i][0] << "\t" << p_o[i][1] << "\t" << p_o[i][2] << "\n";
      cout << "j=117 p_e_rx is " << p_e_rx[i][0] << "\t" << p_e_rx[i][1] << "\t" << p_e_rx[i][2] << "\n";
      cout << "j=117 p_o_rx is " << p_o_rx[i][0] << "\t" << p_o_rx[i][1] << "\t" << p_o_rx[i][2] << "\n";
      cout << "j=117 Pt is " << Pt[i][0] << "\t" << Pt[i][1] << "\t" << Pt[i][2] << "\n";
      cout << "j=117, function evaluated at 1000 depth is " << f1_distances[i]->Eval(sqrt(pow(-1000.-station_depths[i],2) + pow(horizontal_distances[i],2) )) << "\n";
      cout << "j=117, sumphase is " << sumphase << "\n";
    }
    //    f1_distances[i]->FixParameter(0,freq);    
    //f1_distances[i]->SetParameter(1,f1[i]->GetParameter(1));
    //f1[i]->SetParameter(1,deltan_exp[0]);
    //f1_distances[i]->FixParameter(2,horizontal_distances[i]);
    //f1_distances[i]->SetParameter(3,f1[i]->GetParameter(3));
    //f1_distances[i]->FixParameter(4,L_ATTEN);
    //f1_distances[i]->FixParameter(5,station_depths[i]);

    //    cout << "f1_distances Eval is " << i << "\t" << j << "\t" << f1_distances[i]->Eval(vdistances_bigpic[i][j]) << "\n";
    vmag_atten_beam_bigpic[i].push_back(f1_distances[i]->GetParameter(0)*atten*beam[WHICHPOL]);
    vmag_atten_beam_crosspol_bigpic[i].push_back(f1_distances[i]->GetParameter(0)*atten*beam[WHICHPOL]);
    vmag_atten_beam_crosspol_nointerferencefunc_bigpic[i].push_back(f1_distances[i]->Eval(vdistances_bigpic[i][j])*atten*beam[WHICHPOL]);
    vmag_atten_beam_crosspol_func_bigpic[i].push_back(f1_distances[i]->Eval(vdistances_bigpic[i][j])*atten*beam[WHICHPOL]);

    //vmag_atten_func_bigpic[i].push_back(f1_distances[i]->Eval(vdistances_bigpic[i][j]));
    */ 
    } // end loop over NDISTANCES_BIGPIC
    if (i!=5)
      g1_distances[i]=new TGraphErrors(vtotal_distances[i].size(),&vtotal_distances[i][0],&vsnrmax[i][0],&vtotal_distances_err[i][0],&vsnrmax_err[i][0]);
    //  g_atten_func_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_func_bigpic[i][0]);
    g_atten_beam_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_beam_bigpic[i][0]);
    g_atten_beam_crosspol_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_beam_crosspol_bigpic[i][0]);
    g_atten_beam_crosspol_nointerferencefunc_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_beam_crosspol_nointerferencefunc_bigpic[i][0]);
    g_atten_beam_crosspol_func_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vmag_atten_beam_crosspol_func_bigpic[i][0]);
 

    //    g_nsolutions_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vnsolutions_bigpic[i][0]);

    g_sumphase_distances[i]=new TGraph(vdistances_bigpic[i].size(),&vdistances_bigpic[i][0],&vtimediff_bigpic[i][0]);
    }

    

     // end loop 0 to 4
  //}
  
    
  //  
  const int NAXIAL=3;
  const int NCONSTANT=2;
  string sbiaxial[NAXIAL]={"isotropic_","uniaxial_",""};
  string sconstant[NCONSTANT]={"","constant_"};
 
  string sdir="plots_" + sconstant[CONSTANTINDICATRIX] + sbiaxial[BIAXIAL+1] + to_string((int)(freq/1.E6)) + "MHz_angletx" + to_string(CROSSPOLANGLE_TX_INT) + "_anglerx" + to_string(CROSSPOLANGLE_RX_INT) + "/";
  string sname=sdir+"c1.pdf";
  c1->Print(sname.c_str());
  cout << "just printed c1.\n";
  TCanvas *c1b=new TCanvas("c1b","c1b",800,800);
  h2[0]->Draw();
  //  for (int istations=0;istations<NSTATIONS;istations++) {
  for (int istations=minstation;istations<=maxstation;istations++) {
    g_atten_beam_crosspol_func[istations]->SetLineColor(icolors_dave[istations]);
    g_atten_beam_crosspol_func[istations]->Draw("same");
  }
    
  sname=sdir+"c1b.pdf";
  c1b->Print(sname.c_str());
  cout << "printed c1b.\n";


  TCanvas *c3=new TCanvas("c3","c3",800,800);
  c3->Divide(2,3);
  TH2D *h3[6];    

  //  for (int i=0;i<NSTATIONS;i++) {
  for (int i=minstation;i<=maxstation;i++) {
    //if (i==0 || i==3) {
      //  if (i==0 || i==1) {
    c3->cd(i+1);
    h3[i]=new TH2D("h2","h2",3500,1000.,4500,1000,0.,vmax[WHICHPOL][i]);
    //h3[i]=new TH2D("h2","h2",500,1250.,1750,1000,0.,vmax[WHICHPOL][i]/2.);


    h3[i]->Draw();
    //    g1_distances[i]->Fit(f1_distances[i]);
    
    //    f1_distances[i]->Draw("same");
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
    //g_atten_func_distances[i]->Draw("al");

    cout << "f1_distances evaluated at 0 is " << f1_distances[i]->Eval(0.) << "\n";

  }
  
  sname=sdir+"c3.pdf";
  c3->Print(sname.c_str());
    
    
  //  TGraphErrors *g_obs=new TGraphErrors(6,station,deltan_obs,zeroes,deltan_obs_err);
  TGraphErrors *g_obs=new TGraphErrors(6,alpha_deg,deltan_obs,zeroes,deltan_obs_err);
  //TGraph *g_exp=new TGraph(6,station,deltan_exp);
  TGraph *g_exp=new TGraph(6,alpha_deg,deltan_exp);


  TF1 *fgetDeltaN=new TF1("f1","2*([0]-([0]-[1]/2.)/sqrt(1-[1]/[0]*cos(x*3.14159/180.)*cos(x*3.14159/180.)*(1-[1]/(4.*[0]))))",0.,90.);
  fgetDeltaN->SetParameter(0,NICE);
  fgetDeltaN->SetParameter(1,DELTAN);


    cout << "function evaluated is " << fgetDeltaN->Eval(10.) << "\t" << fgetDeltaN->Eval(90.) << "\n";

  //  TH2D *h2_compare=new TH2D("h2_compare","h2_compare",100,0.,5.5,100,1.E-9,0.01);
    TH2D *h2_compare=new TH2D("h2_compare","h2_compare",100,0.,180.,100,1.E-4,1.0);
  //TH2D *h2_compare=new TH2D("","",100,0.,90.,100,1.E-5,0.01);
  h2_compare->SetXTitle("#alpha (degrees)");
  h2_compare->SetYTitle("#Delta n");
  TCanvas *c2=new TCanvas("c2","c2",800,800);
  c2->SetLogy();
  gStyle->SetOptStat(0);
  h2_compare->Draw();
  fgetDeltaN->SetLineColor(kBlack);
  //fgetDeltaN->Draw("same");

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


  // TCanvas *c4=new TCanvas("c4","c4",800,800);
  // c4->Divide(2,3);
  // TH2D *h4[6];
  //   for (int i=0;i<NSTATIONS;i++) {
  //     c4->cd(i+1);

  //     h4[i]=new TH2D("h2","h2",3500,1000.,4500,1000,-0.5,2.5);
  //   //h3[i]=new TH2D("h2","h2",500,1250.,1750,1000,0.,vmax[WHICHPOL][i]/2.);

  //     h4[i]->Draw();

  //     g_nsolutions_distances[i]->SetMarkerStyle(21);
  //     g_nsolutions_distances[i]->SetMarkerColor(kBlack);
  //     g_nsolutions_distances[i]->Draw("psame");
  //   }

  //   c4->Print("nsolutions.pdf");

//     TCanvas *c5=new TCanvas("c5","c5",800,800);
  
//     TH2D *h5=new TH2D("","",100,-1500.,-500.,100,0.01,50.);
//     h5->SetXTitle("Pulser height (m)");
//     h5->SetYTitle("Angle between usual k (no birefring.) and k_{1}, k_{2} (degrees)");
//     h5->GetYaxis()->SetTitleOffset(1.2);
//     gPad->SetLogy();
//     h5->Draw();
   
//     TGraph *g_angle1[6];
//     TGraph *g_angle2[6];
//     TGraph *g_angle1_2[6];
  TGraph *g_angleS1_k[6];
  TGraph *g_angleS2_k[6];


//     //    for (int i=0;i<NSTATIONS;i++) {
//     for (int i=minstation;i<=maxstation;i++) {
//       g_angle1[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vangle_khat_0_khat_1_2[i][0]);
//       g_angle2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vangle_khat_0_khat_2_2[i][0]);
     
//       g_angle1[i]->SetMarkerColor(icolors[i]);
//       g_angle1[i]->SetLineColor(icolors[i]);
//       g_angle1[i]->SetLineWidth(2);

//       g_angle2[i]->SetMarkerColor(icolors[i]);
//       g_angle2[i]->SetLineColor(icolors[i]);
//       g_angle2[i]->SetLineStyle(kDashed);
//       g_angle2[i]->SetLineWidth(2);



//       g_angle1[i]->Draw("lsame");    
//       g_angle2[i]->Draw("lsame");    
//     }


//     sname=sdir+"angle_khat_0_khat_1_2.pdf";
//     c5->Print(sname.c_str());
    
//     TCanvas *c6=new TCanvas("c6","c6",800,800);
  
//     TH2D *h6=new TH2D("","",100,-1500.,-500.,100,0.01,50.);
//     h6->SetXTitle("Pulser height (m)");
//     h6->SetYTitle("Angle between k_{1} and k_{2} (degrees)");
//     h6->GetYaxis()->SetTitleOffset(1.2);
//     gPad->SetLogy();
//     h6->Draw();

//     //    for (int i=0;i<NSTATIONS;i++) {
//     for (int i=minstation;i<=maxstation;i++) {

//       g_angle1_2[i]=new TGraph(vdepth[i].size(),&vdepth[i][0],&vangle_khat_1_2_khat_2_2[i][0]);

//       g_angle1_2[i]->SetMarkerColor(icolors[i]);
//       g_angle1_2[i]->SetLineColor(icolors[i]);
//       g_angle1_2[i]->SetLineStyle(kSolid);
//       g_angle1_2[i]->SetLineWidth(2);

//       g_angle1_2[i]->Draw("lsame");    

//     }
//     sname=sdir+"angle_khat_1_1_khat_1_2.pdf";
//     c6->Print(sname.c_str());

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
//     for (int i=minstation;i<=maxstation;i++) {

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


// TCanvas *c8=new TCanvas("c8","c8",800,800);
//     c8->Divide(2,2);
//     TH2D *h8_z_x=new TH2D("h8_z_x","h8_z_x",100,-5000.,0.,100,-1800.,200.);
//     c8->cd(1);
//     h8_z_x->Draw();
//     //    for (int istations=0;istations<NSTATIONS;istations++) {
//     for (int istations=minstation;istations<=maxstation;istations++) {
//       //if (istations==0)
//       //graypath_z_x[istations]->Draw("al");
//       //else
// 	graypath_z_x[istations]->Draw("lsame");
	
//     }

//     c8->cd(2);
//     TH2D *h8_z_y=new TH2D("h8_z_y","h8_z_y",100,-5000.,5000.,100,-1800.,200.);
//     //h6->Draw();
//     h8_z_y->SetXTitle("y (m)");
//     h8_z_y->SetYTitle("z (m)");

//     h8_z_y->Draw();
//     //    for (int istations=0;istations<NSTATIONS;istations++) {
//     for (int istations=minstation;istations<=minstation;istations++) {
//     // if (istations==0)
//     //	graypath_z_y[istations]->Draw("al");
//     //else
//     graypath_z_y[istations]->Draw("lsame");
//     }
    
//     c8->cd(3);
//     //h6->Draw();
//     TH2D *h8_y_x=new TH2D("h8_y_x","h8_y_x",100,-5000.,5000.,100,-5000.,5000.);
//     h8_y_x->SetXTitle("x (m)");
//     h8_y_x->SetYTitle("y (m)");
//     h8_y_x->Draw();
//     //    for (int istations=0;istations<NSTATIONS;istations++) {
//     for (int istations=minstation;istations<=minstation;istations++) {
//       //if (istations==0)
//       //graypath_y_x[istations]->Draw("al");
//       //else
// 	graypath_y_x[istations]->Draw("lsame");	
//     }
//     sname=sdir+"raypath.pdf";
//     c8->Print(sname.c_str());


    TCanvas *c9=new TCanvas("c9","c9",800,800);    
    TH2D *h9=new TH2D("","",100,-1600.,-600.,100,-20.,90.);
    auto legend26a = new TLegend(0.55,0.55,0.88,0.88);
    legend26a->SetBorderSize(0);
    legend26a->SetTextSize(0.05);

    //    titles(h9, "", "Pulser height (m)", "Time difference (ns)");


    h9->GetXaxis()->SetTitle("Pulser height (m)");
    h9->GetYaxis()->SetTitle("Time difference (ns)");
    h9->GetXaxis()->SetTitleOffset(1.1);
    h9->GetYaxis()->SetTitleOffset(1.2);
    h9->GetXaxis()->SetTitleSize(0.04);
    h9->GetYaxis()->SetTitleSize(0.04);

    h9->GetXaxis()->SetNdivisions(504);
    h9->Draw();
    //    gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {

      g_sumphase[istations]->SetMarkerColor(icolors[istations]);
      g_sumphase[istations]->SetLineColor(icolors[istations]);
      g_sumphase[istations]->SetLineStyle(kSolid);
      g_sumphase[istations]->SetLineWidth(3);

      //      if (istations==0)
      //g_sumphase[istations]->Draw("al");
	//else
      g_sumphase[istations]->Draw("lsame");
	
      legend26a->AddEntry(g_sumphase[istations],snames[istations].c_str(),"l");	
	//}
      


    }
    legend26a->Draw("same");
    sname=sdir+"sumphase.pdf";
    c9->Print(sname.c_str());

    TCanvas *c9a=new TCanvas("c9a","c9a",800,800);   
    c9a->Divide(2,3);
    TH2D *h9a=new TH2D("h9a","h9a",100,-1800.,0.,100,-1.1,1.1);
    
    h9a->GetXaxis()->SetTitle("Pulser height (m)");
    h9a->GetYaxis()->SetTitle("Not flipped");
    h9a->GetXaxis()->SetTitleOffset(1.1);
    h9a->GetYaxis()->SetTitleOffset(1.45);
    //    gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      c9a->cd(istations+1);
      h9a->Draw();  
      g_notflipped[istations]->SetMarkerColor(icolors[istations]);
      g_notflipped[istations]->SetLineColor(icolors[istations]);
      g_notflipped[istations]->SetLineStyle(kSolid);
      g_notflipped[istations]->SetLineWidth(2);

      //      if (istations==0)
      //g_sumphase[istations]->Draw("al");
	//else
      g_notflipped[istations]->Draw("lsame");
	
    }
    sname=sdir+"notflipped.pdf";
    c9a->Print(sname.c_str());

    TCanvas *c9b=new TCanvas("c9b","c9b",800,800);    
    TH2D *h9b=new TH2D("h9b","h9b",100,0.,5000.,100,0.,5000.);
    h9b->Draw();
    //    gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=minstation;istations++) {

      g_sumlength[istations]->SetMarkerColor(icolors[istations]);
      g_sumlength[istations]->SetLineColor(icolors[istations]);
      g_sumlength[istations]->SetLineStyle(kSolid);
      g_sumlength[istations]->SetLineWidth(2);

      //      if (istations==0)
      //g_sumphase[istations]->Draw("al");
	//else
      g_sumlength[istations]->Draw("lsame");
	
    }
    sname=sdir+"sumlength.pdf";
    c9b->Print(sname.c_str());

    TCanvas *c10=new TCanvas("c10","c10",800,800);
    c10->Divide(2,3);
    TH2D *h10=new TH2D("h10","h10",100,freqmin,freqmax,100,1.,pmax[0][0]);
    //  h10->Draw();
    //    gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    //    for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {
    //    for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      c10->cd(istations+1);
      //      for (int ispecialdepths=0;ispecialdepths<NSPECIALDEPTHS;ispecialdepths++) {
      //for (int ispecialdepths=0;ispecialdepths<g_power[istations].size();ispecialdepths++) {
      
      g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->SetMarkerColor(icolors[istations]);
      g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->SetLineColor(icolors[istations]);
      g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->SetLineStyle(kSolid);
      g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->SetLineWidth(2);
	
	// //      if (istations==0)
	// g_power[istations][ispecialdepths]->Draw("lsame");
      //}
	//else
      g_spectra[istations][(int)g_idepth[istations]->Eval(depth_special)]->Draw("al");
      //}
	
    }
    sname=sdir+"spectra.pdf";
    c10->Print(sname.c_str());
   TCanvas *c11=new TCanvas("c11","c11",800,800);
    
   //    TH2D *h11=new TH2D("h11","h11",100,0.,4000.,100,1.E-5,1.E-2);
    TH2D *h11=new TH2D("h11","h11",100,-1200.,0.,100,1.E-5,1.E-2);
    string stitle ="Depth along path for pulser depth of " + to_string(depth_special) + "\n";
    h11->SetXTitle("Depth along path (m)");
    h11->SetYTitle("#Delta n");

    h11->Draw();
    gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      //for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {

      // g_deltan[istations]->SetMarkerColor(icolors[istations]);
      // g_deltan[istations]->SetLineColor(icolors[istations]);
      // g_deltan[istations]->SetLineStyle(kSolid);
      // g_deltan[istations]->SetLineWidth(2);

      // //      if (istations==0)
      // g_deltan[istations]->Draw("lsame");

      g_deltan_pulserdepth[istations]->SetMarkerColor(icolors[istations]);
      g_deltan_pulserdepth[istations]->SetLineColor(icolors[istations]);
      g_deltan_pulserdepth[istations]->SetLineStyle(kSolid);
      g_deltan_pulserdepth[istations]->SetLineWidth(2);

      //      if (istations==0)
      g_deltan_pulserdepth[istations]->Draw("lsame");
	//else
      //g_deltan[istations]->Draw("lsame");
	
    }
    sname=sdir+"deltan.pdf";
    c11->Print(sname.c_str());

    TCanvas *c11b=new TCanvas("c11b","c11b",800,800);

   TH2D *h11b=new TH2D("h11b","h11b",100,0.,4000.,100,1.,1.E5);
    //h11->Draw();
    gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      //for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {

      g_attenlengths[istations]->SetMarkerColor(icolors[istations]);
      g_attenlengths[istations]->SetLineColor(icolors[istations]);
      g_attenlengths[istations]->SetLineStyle(kSolid);
      g_attenlengths[istations]->SetLineWidth(2);

      //      if (istations==0)
      g_attenlengths[istations]->Draw("al");
	//else
      //g_deltan[istations]->Draw("lsame");
	
    }
    sname=sdir+"attenlengths.pdf";
    c11b->Print(sname.c_str());

    TCanvas *c11c=new TCanvas("c11c","c11c",800,800);

    TH2D *h11c=new TH2D("h11c","h11c",100,0.,2000.,100,-1.1,1.1);
    h11c->Draw();
    //gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      //for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {

      g_notflipped_alongpath[istations]->SetLineColor(icolors[istations]);
      g_notflipped_alongpath[istations]->SetLineStyle(kSolid);
      g_notflipped_alongpath[istations]->SetLineWidth(2);
      g_notflipped_alongpath[istations]->Draw("lsame");
	
    }
    sname=sdir+"notflipped_alongpath.pdf";
    c11c->Print(sname.c_str());

   TCanvas *c11d=new TCanvas("c11d","c11d",800,800);

    TH2D *h11d=new TH2D("h11d","h11d",100,0.,2000.,100,-180.,180.);
    h11d->Draw();
    //gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      //for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {

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

  TCanvas *c11e=new TCanvas("c11e","c11e",800,800);

    TH2D *h11e=new TH2D("h11e","h11e",100,0.,2000.,100,0.,180.);
    h11e->Draw();
    h11e->SetXTitle("Phi (degrees)");
    h11e->SetYTitle("Theta (degrees)");
    
    //gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      //for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {

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

TCanvas *c11f=new TCanvas("c11f","c11f",800,800);

    TH2D *h11f=new TH2D("h11f","h11f",100,-180.,180.,100,0.,180.);
    h11f->Draw();
    //gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      //for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {

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

    TCanvas *c12=new TCanvas("c12","c12",800,800);
    
    TH2D *h12=new TH2D("h12","h12",100,0.,2000.,100,0.,200.);
    //h11->Draw();
    gPad->SetLogy();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      //for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {
      
      g_idepth[istations]->SetMarkerColor(icolors[istations]);
      g_idepth[istations]->SetLineColor(icolors[istations]);
      g_idepth[istations]->SetLineStyle(kSolid);
      g_idepth[istations]->SetLineWidth(2);
      
      //      if (istations==0)
      g_idepth[istations]->Draw("al");
      //else
      //g_deltan[istations]->Draw("lsame");
      
    }
    sname=sdir+"index.pdf";
    c12->Print(sname.c_str());

   TCanvas *c13=new TCanvas("c13","c13",800,800);
    
    TH2D *h13a=new TH2D("h13a","h13a",100,-2000.,0.,100,0.,180.);
    TH2D *h13b=new TH2D("h13b","h13b",100,-2000.,0.,100,0.,180.);
    TH2D *h13c=new TH2D("h13c","h13c",100,0.,180.,100,0.,180.);
    //h11->Draw();
    c13->Divide(2,2);
    gPad->SetLogy();


    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      //for (int istations=NSTATIONS-1;istations<NSTATIONS;istations++) {
      
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
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      h13a->SetXTitle("Receive angles (degrees)");
      g_receive[istations]->Draw("lsame");
    }
    c13->cd(2);
    h13b->Draw();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      h13b->SetXTitle("Launch angles (degrees)");
      g_launch[istations]->Draw("lsame");
    }
    c13->cd(3);
    h13c->Draw();
    //    for (int istations=0;istations<NSTATIONS;istations++) {
    for (int istations=minstation;istations<=maxstation;istations++) {
      h13c->SetXTitle("Receive angle (degrees)");
      h13c->SetYTitle("Launch angle (degrees)");

      g_receive_launch[istations]->Draw("lsame");
    }  

    sname=sdir+"receive_launch.pdf";
    c13->Print(sname.c_str());


    TCanvas *c14=new TCanvas("c14","c14",800,800);
    
    TH2D *h14=new TH2D("h14","h14",100,-1800.,0.,100,0.,90.);
    h14->Draw();
    h14->SetTitle("");
    h14->GetXaxis()->SetTitleOffset(1.2);
    h14->SetXTitle("Height (meters)");
    h14->SetYTitle("V_{z} (degrees)");

    //gPad->SetLogy();
    
     
    g_V->SetLineColor(kBlack);
    g_V->SetLineStyle(kSolid);
    g_V->SetLineWidth(2);
    
    g_V->Draw("lsame");
      
    
    sname=sdir+"V.pdf";
    c14->Print(sname.c_str());


   TCanvas *c15=new TCanvas("c15","c15",800,800);
   c15->SetLeftMargin(0.15);
   c15->SetBottomMargin(0.15);
   //c15->Divide(2,1,0,0);
   // c15->cd(2);

   // //    gPad->SetLeftMargin(0.14);
   // //gPad->SetBottomMargin(0.14);
   // //gPad->SetTopMargin(0.14);
   // //gPad->SetRightMargin(0.14);
     // TH2D *h15=new TH2D("h15","h15",100,-1800.,0.,100,1.30,1.79);
    
     // h15->SetTitle("");
     // h15->SetXTitle("Height (meters)");
     // h15->GetYaxis()->SetTitleOffset(1.45);
     // h15->GetXaxis()->SetTitleOffset(1.1);
     // h15->GetXaxis()->SetTitleSize(0.04);
     // h15->GetYaxis()->SetTitleSize(0.04);
     // h15->GetXaxis()->SetNdivisions(504);
     // h15->SetYTitle("Principal axis");
     // //h15->GetYaxis()->SetTickLength(0.);
     // h15->Draw();
   //  //gPad->SetLogy();
   
      
   //gn1->SetMarkerColor(kBlue);
      gn1->SetLineColor(kGray+1);
      gn1->SetLineStyle(kSolid);
      gn1->SetLineWidth(3);      
      //      gn1->Draw("lsame");

      // gn2->SetMarkerColor(kBlue);
      gn2->SetLineColor(kGray+2);
      gn2->SetLineStyle(kSolid);
      gn2->SetLineWidth(3);      
      //gn2->Draw("lsame");

      //gn3->SetMarkerColor(kBlue);
      gn3->SetLineColor(kGray+3);
      gn3->SetLineStyle(kSolid);
      gn3->SetLineWidth(3);      
      //gn3->Draw("lsame");
      
      // // 1.78-0.43*EXP(0.0132*A3)
      // TF1 *fuzair=new TF1("fuzair","[0]+[1]*exp(x*[2])",-1800.,0.);
      // fuzair->SetParameter(0,1.78);
      // fuzair->SetParameter(1,-0.43);
      // fuzair->SetParameter(2,0.0132);
      // fuzair->SetLineColor(kBlack);
      // fuzair->SetLineStyle(kDashed);
      // fuzair->SetLineWidth(3);

      // fuzair->Draw("lsame");

      // TF1 *fuzair_amyadjusted=new TF1("fuzair_amyadjusted","[0]+[1]*exp(x*[2])",-1800.,0.);
      // fuzair_amyadjusted->SetParameter(0,1.78);
      // fuzair_amyadjusted->SetParameter(1,-0.43);
      // fuzair_amyadjusted->SetParameter(2,0.03624);
      // fuzair_amyadjusted->SetLineColor(kBlack);
      // fuzair_amyadjusted->SetLineStyle(kDotted);
      // fuzair_amyadjusted->SetLineWidth(3);

      // fuzair_amyadjusted->Draw("lsame");

      // // auto legend = new TLegend(0.1,0.7,0.48,0.9);
      // auto legend = new TLegend(0.2,0.15,0.4,0.42);
      
      // //  C=0.0132;
      // // ***********  I'm changing this so that n is always between n1,n2 and n3
      // //C=0.03624; 
      // //      legend->AddEntry("fuzair","1.78-0.45 e^{0.0132 x}","l");
      // //legend->AddEntry("fuzair_amyadjusted","1.78-0.45 e^{0.03624 x}","l");
      

      // legend->AddEntry("gn1","n1","l");
      // legend->AddEntry("gn2","n2","l");
      // legend->AddEntry("gn3","n3","l");
      // legend->AddEntry("fuzair","n (Latif)","l");
      // legend->AddEntry("fuzair_amyadjusted","n (Latif, adjusted)","l");
      // //legend->SetTextSize(0.03);
      // legend->SetTextSize(0.04);
      // legend->SetBorderSize(0);
      // legend->Draw("same");

      // TF1 *fy=new TF1("fy","x",1.60,1.79);
      // TGaxis *A1 = new TGaxis(0.,1.57,0.,1.79,"fy",510,"L+");
      // A1->SetTitle("Index of refraction");
      // //A1->SetTitleOffset(0.8);
      // A1->SetTextFont(42);
      // A1->SetLabelFont(42);
      // A1->SetTitleSize(0.05);
      // A1->Draw("same");


      //      c15->cd(1);
      //    gPad->SetLeftMargin(0.14);
      //gPad->SetBottomMargin(0.14);
      //gPad->SetTopMargin(0.14);
      //gPad->SetRightMargin(0.14);

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

      // TCanvas *c16b=new TCanvas("c16b","c16b",800,400);
      // c16b->Divide(2,1);
      // c16b->cd(1);
      // h15->Draw();

      // c16b->cd(2);
      // h16->Draw();


      // c16b->Print("n123_sidebyside.pdf");






      TCanvas *c17=new TCanvas("c17","c17",800,800);
      TH2D *h17=new TH2D("h17","h17",100,-1800.,0.,100,-1100.,200.);
      c17->Divide(1,3);
      c17->cd(1);
      h17->Draw();
      //      for (int istations=0;istations<NSTATIONS;istations++) {
      for (int istations=minstation;istations<=maxstation;istations++) {
	g_output6[istations]->SetLineColor(icolors[istations]);
	g_output6[istations]->SetLineWidth(2);
	g_output6[istations]->Draw("lsame");
      }

      c17->cd(2);
      h17->Draw();
      //      for (int istations=0;istations<NSTATIONS;istations++) {
      for (int istations=minstation;istations<=maxstation;istations++) {
	g_output7[istations]->SetLineColor(icolors[istations]);
	g_output7[istations]->SetLineWidth(2);
	g_output7[istations]->Draw("lsame");
      }

      c17->cd(3);
      h17->Draw();
      //      for (int istations=0;istations<NSTATIONS;istations++) {
      for (int istations=minstation;istations<=maxstation;istations++) {
	g_output8[istations]->SetLineColor(icolors[istations]);
	g_output8[istations]->SetLineWidth(2);
	g_output8[istations]->Draw("lsame");
      }

      sname=sdir+"outputs.pdf";
      c17->Print(sname.c_str());

      TCanvas *c18=new TCanvas("c18","c18",800,800);
      TH2D *h18a=new TH2D("h18","h18",100,-1800.,0.,100,0.,180.);
      TH2D *h18b=new TH2D("h18","h18",100,-1800.,0.,100,-90.,90.);
      c18->Divide(3,4);


      //      for (int istations=0;istations<NSTATIONS;istations++) {
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

      //h18->Draw();
      //      for (int istations=0;istations<NSTATIONS;istations++) {
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

      TCanvas *c20=new TCanvas("c20","c20",800,800);
      TH2D *h20a=new TH2D("h20","h20",100,-1800.,0.,100,0.,180.);
      TH2D *h20b=new TH2D("h20","h20",100,-1800.,0.,100,-90.,90.);
      c20->Divide(3,4);


      //      for (int istations=0;istations<NSTATIONS;istations++) {
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

      //h20->Draw();
      //      for (int istations=0;istations<NSTATIONS;istations++) {
      for (int istations=minstation;istations<=maxstation;istations++) {
	c20->cd(NSTATIONS+istations+1);

	gtxdepth_theta2_Sclock[istations]->SetLineColor(icolors[istations]);
	gtxdepth_theta2_Sclock[istations]->SetLineWidth(2);
	
	gtxdepthE_theta2_Sclock[istations]->SetLineColor(icolors[istations]);
	gtxdepthE_theta2_Sclock[istations]->SetLineWidth(2);
	gtxdepthE_theta2_Sclock[istations]->SetLineStyle(kDashed);
	
	//	if (!(istations==0 || istations==5)) {
	h20b->Draw();
	gtxdepth_theta2_Sclock[istations]->Draw("lsame");
	gtxdepthE_theta2_Sclock[istations]->Draw("lsame");
	  //}
	  //else {
	  //gtxdepth_theta2_Sclock[istations]->Draw("al");
	  //gtxdepthE_theta2_Sclock[istations]->Draw("lsame");

	  //	}
      }
      sname=sdir+"anglesontheSclock.pdf";
      c20->Print(sname.c_str());
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
 

 


      //      TCanvas *c24=makePretty2Panel();
      TCanvas *c24=new TCanvas("c24","c24",1600,800);
      c24->Divide(2,1);
      //c24->SetLeftMargin(0.15);
      //c24->SetBottomMargin(0.15);
      c24->cd(1);
gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.01);
      TH2D *h24a=new TH2D("","",100,-1600.,-600.,100,-10.,70.);
      titles(h24a, "", "Pulser height (m)", "#epsilon_{ 1}^{ T}");

      //      h24a->GetXaxis()->SetTitle("Pulser height (m)");
      h24a->GetXaxis()->SetNdivisions(504);
      h24a->GetYaxis()->SetNdivisions(504);
      // h24a->GetYaxis()->SetTitle("#epsilon_{1,T}, #epsilon_{2,T}");
      //h24a->GetXaxis()->SetTitleOffset(0.6);
      //h24a->GetYaxis()->SetTitleOffset(0.6);
      //h24a->GetXaxis()->SetTitleSize(0.08);
      //h24a->GetYaxis()->SetTitleSize(0.07);

      h24a->Draw();
      auto legend24 = new TLegend(0.18,0.2,0.36,0.54);
      legend24->SetBorderSize(0);
      legend24->SetTextSize(0.04);
      for (int istations=minstation;istations<=maxstation;istations++) {
	//if (istations!=5) {
	  gepsilon1_tx[istations]->SetLineColor(icolors[istations]);
	  gepsilon1_tx[istations]->SetLineWidth(2);
	  gepsilon1_tx[istations]->Draw("lsame");
	  gepsilon2_tx[istations]->SetLineColor(icolors[istations]);
	  gepsilon2_tx[istations]->SetLineWidth(2);
	  gepsilon2_tx[istations]->SetLineStyle(kDashed);
	  //gepsilon2_tx[istations]->Draw("al");
	  gepsilon1_tx[istations]->SetName(snames[istations].c_str());	
	  legend24->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");	
	  //}
      }
      //      legend24->Draw("same");
      
      c24->cd(2);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.01);
      TH2D *h26a=new TH2D("","",100,-1600.,-600.,100,-10.,70.);
      
      auto legend26 = new TLegend(0.55,0.5,0.88,0.8);

      titles(h26a, "", "Pulser height (m)", "#epsilon_{ 1}^{ R}");

      //      h26a->GetXaxis()->SetTitle("Pulser height (m)");
      h26a->GetXaxis()->SetNdivisions(504);
      h26a->GetYaxis()->SetNdivisions(504);
      // h26a->GetYaxis()->SetTitle("#epsilon_{1,R}, #epsilon_{2,R}");
      // h26a->GetXaxis()->SetTitleOffset(0.6);
      // h26a->GetYaxis()->SetTitleOffset(0.6);
      // h26a->GetXaxis()->SetTitleSize(0.08);
      // h26a->GetYaxis()->SetTitleSize(0.08);
      h26a->Draw();
      legend26->SetBorderSize(0);
      legend26->SetTextSize(0.05);
      for (int istations=minstation;istations<=maxstation;istations++) {
	//if (istations!=5) {
	gepsilon1_rx[istations]->SetLineColor(icolors[istations]);
	gepsilon1_rx[istations]->SetLineWidth(2);
	gepsilon1_rx[istations]->Draw("lsame");
	gepsilon2_rx[istations]->SetLineColor(icolors[istations]);
	gepsilon2_rx[istations]->SetLineWidth(2);
	gepsilon2_rx[istations]->SetLineStyle(kDashed);
	//gepsilon2_rx[istations]->Draw("lsame");	
	gepsilon1_rx[istations]->SetName(snames[istations].c_str());	
	legend26->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");	
	//}
      }
      legend26->Draw("same");
      sname=sdir+"epsilons_sidebyside.pdf";
      c24->Print(sname.c_str());

      TCanvas *c24b=new TCanvas("c24b","c24b",1600,800);
      c24b->Divide(2,1);
      //c24->SetLeftMargin(0.15);
      //c24->SetBottomMargin(0.15);
      c24b->cd(1);
      gPad->SetLeftMargin(0.18);
      gPad->SetBottomMargin(0.17);
      gPad->SetRightMargin(0.01);
      TH2D *h24b=new TH2D("","",100,-1600.,-600.,100,-0.15,0.15);
      titles(h24b, "", "Pulser height (m)", "#epsilon_{ 1}^{ T} - #epsilon_{ 2}^{ T} (^{o})");

      //      h24a->GetXaxis()->SetTitle("Pulser height (m)");
      h24b->GetXaxis()->SetNdivisions(504);
      h24b->GetYaxis()->SetNdivisions(504);
      
      // h24a->GetYaxis()->SetTitle("#epsilon_{1,T}, #epsilon_{2,T}");
      h24b->GetXaxis()->SetTitleOffset(1.2);
      h24b->GetYaxis()->SetTitleOffset(1.7);
      //h24a->GetXaxis()->SetTitleSize(0.08);
      //h24a->GetYaxis()->SetTitleSize(0.07);

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
	//if (istations==5) {
	  gdiffepsilon_tx[istations]->SetLineColor(icolors[istations]);
	  gdiffepsilon_tx[istations]->SetLineWidth(2);
	  gdiffepsilon_tx[istations]->Draw("lsame");
	  //	  gdiffepsilon_tx[istations]->Draw("al");
	  legend24b->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");	
	  if (istations<5)
	    legend26b->AddEntry(gdiffepsilon_rx[istations],snames[istations].c_str(),"l");
	  else if (istations==5)
	    legend26c->AddEntry(gdiffepsilon_rx[istations],snames[istations].c_str(),"l");	
	  //}
      }
      //      legend24->Draw("same");
      
      c24b->cd(2);
      gPad->SetLeftMargin(0.18);
      gPad->SetBottomMargin(0.17);
      gPad->SetRightMargin(0.01);
      TH2D *h26b=new TH2D("","",100,-1600.,-600.,100,-0.05,0.05);
      
      titles(h26b, "", "Pulser height (m)", "#epsilon_{ 1}^{ R} - #epsilon_{ 2}^{ R} (^{o})");

      //      h26a->GetXaxis()->SetTitle("Pulser height (m)");

      h26b->GetXaxis()->SetNdivisions(504);
      h26b->GetYaxis()->SetNdivisions(504);
      // h26a->GetYaxis()->SetTitle("#epsilon_{1,R}, #epsilon_{2,R}");
      // h26a->GetXaxis()->SetTitleOffset(0.6);
      h26b->GetYaxis()->SetTitleOffset(1.70);
      h26b->GetXaxis()->SetTitleOffset(1.2);
      // h26a->GetYaxis()->SetTitleSize(0.08);
      h26b->Draw();

      for (int istations=minstation;istations<=maxstation;istations++) {
	//if (istations!=5) {
	gdiffepsilon_rx[istations]->SetLineColor(icolors[istations]);
	gdiffepsilon_rx[istations]->SetLineWidth(2);
	gdiffepsilon_rx[istations]->Draw("lsame");
	

	//}
      }
      legend26b->Draw("same");      
      legend26c->Draw("same");      
      sname=sdir+"diffepsilons_sidebyside.pdf";
      c24b->Print(sname.c_str());

















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

      TCanvas *c28=new TCanvas("c28","c28",800,800);
      TH2D *h28a=new TH2D("h28a","h28a",100,-1800.,0.,100,0.,100000.);
      h28a->Draw();
      for (int istations=minstation;istations<=maxstation;istations++) {
	g_depth_istep[istations]->Draw("lsame");
      }
      sname=sdir+"depth_idepth.pdf";
      c28->Print(sname.c_str());

      TCanvas *c29=new TCanvas("c29","c29",800,800);
      TH2D *h29a=new TH2D("h29","h29",100,-1800.,0.,100,0.,1.);
      TH2D *h29b=new TH2D("h29","h29",100,-1800.,0.,100,0.,1.);
      c29->Divide(3,4);

      //      for (int istations=0;istations<NSTATIONS;istations++) {
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

      TCanvas *c30=new TCanvas("c30","c30",800,800);
      TH2D *h30a=new TH2D("h30","h30",100,-1800.,0.,100,0.,1.);
      TH2D *h30b=new TH2D("h30","h30",100,-1800.,0.,100,0.,1.);
      c30->Divide(3,4);

      //      for (int istations=0;istations<NSTATIONS;istations++) {
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

     TCanvas *c32=new TCanvas("c32","c32",800,800);
      TH2D *h32a=new TH2D("h32","h32",100,-1800.,0.,100,85.,95.);
      TH2D *h32b=new TH2D("h32","h32",100,-1800.,0.,100,-5.,5.);
      c32->Divide(3,4);


      //      for (int istations=0;istations<NSTATIONS;istations++) {
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

      //h32->Draw();
      //      for (int istations=0;istations<NSTATIONS;istations++) {
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

     TCanvas *c33=new TCanvas("c33","c33",800,800);
      TH2D *h33a=new TH2D("h33","h33",100,-1800.,0.,100,0.,180.);
      TH2D *h33b=new TH2D("h33","h33",100,-1800.,0.,100,-90.,90.);
      c33->Divide(3,4);


      //      for (int istations=0;istations<NSTATIONS;istations++) {
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

      //h33->Draw();
      //      for (int istations=0;istations<NSTATIONS;istations++) {
      for (int istations=minstation;istations<=maxstation;istations++) {
	c33->cd(NSTATIONS+istations+1);

	grxdepth_theta2_Sclock[istations]->SetLineColor(icolors[istations]);
	grxdepth_theta2_Sclock[istations]->SetLineWidth(2);
	
	grxdepthE_theta2_Sclock[istations]->SetLineColor(icolors[istations]);
	grxdepthE_theta2_Sclock[istations]->SetLineWidth(2);
	grxdepthE_theta2_Sclock[istations]->SetLineStyle(kDashed);
	
	//	if (!(istations==0 || istations==5)) {
	h33b->Draw();
	grxdepth_theta2_Sclock[istations]->Draw("lsame");
	grxdepthE_theta2_Sclock[istations]->Draw("lsame");
	  //}
	  //else {
	  //grxdepth_theta2_Sclock[istations]->Draw("al");
	  //grxdepthE_theta2_Sclock[istations]->Draw("lsame");

	  //	}
      }
      sname=sdir+"anglesontheSclock_rx.pdf";
      c33->Print(sname.c_str());

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

      TCanvas *c36_a5=new TCanvas("c36_a5","c36_a5",800,800);
      h36a[4]->Draw();
      h36a[4]->GetYaxis()->SetRangeUser(0.,25.);
      grxdepth_atten[4]->Draw("lsame");
      grxdepth_atten_beam[4]->Draw("lsame");
      gvenvelope_plus_r1[4]->Draw("lsame");
      gvoltage_r1[4]->Draw("lsame");

      sname=sdir+"voltages_a5.pdf";
      c36_a5->Print(sname.c_str());

      //      TCanvas *c37=makePretty2Panel();
      TCanvas *c37=new TCanvas("c37","c37",1600,800);
      //c37->SetLeftMargin(0.15);
      //c37->SetBottomMargin(0.15);
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

      //h37->GetXaxis()->SetTitleOffset(0.8);
      h37->GetYaxis()->SetTitleOffset(1.2);
      //h37->GetXaxis()->SetTitleSize(0.05);
      //h37->GetYaxis()->SetTitleSize(0.05);
      //h37->GetXaxis()->SetNdivisions(504);

      auto legend38b = new TLegend(0.17,0.76,0.62,0.89);
      legend38b->SetBorderSize();
      legend38b->SetTextSize(0.05);

      string snametemp="Voltage @ 300 MHz";
      h37->Draw();
      //      texttemp_hpol->Draw("same");
      // h37->GetXaxis()->SetTitle("Pulser height (m)");
      // h37->GetYaxis()->SetTitle("r V_{#phi} (arb. units)");
      // h37->GetXaxis()->SetTitleOffset(0.8);
      // h37->GetYaxis()->SetTitleOffset(0.85);
      // h37->GetXaxis()->SetTitleSize(0.09);
      // h37->GetYaxis()->SetTitleSize(0.09);
      h37->GetXaxis()->SetNdivisions(504);
      
      for (int istations=minstation;istations<=maxstation;istations++) {
	//if (istations!=5) {
	gvenvelope_plus_r2[istations]->SetLineStyle(kSolid);
	//if (istations!=5)

	gvenvelope_minus_r2[istations]->SetLineStyle(kSolid);
	//if (istations==5)

	gvoltage_r2[istations]->SetLineStyle(kDashed);

	if (istations!=5) {
	  gvenvelope_minus_r2[istations]->Draw("lsame");
	  gvenvelope_plus_r2[istations]->Draw("lsame");	
	  gvoltage_r2[istations]->Draw("lsame");
	}
	gvoltage_r2[istations]->SetName(snames[istations].c_str());
	legend37->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");
	//	}


	if (istations==0) {
	  legend38b->AddEntry(gvenvelope_plus_r2[istations],"Voltage envelope","l");
	  legend38b->AddEntry(gvoltage_r2[istations],snametemp.c_str(),"l");
	}

      }
      //      legend37->Draw("same");
      legend38b->Draw("same");


      
      c37->cd(1);


 
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetLogy();
      //      TH2D *h38=new TH2D("","",100,-1600.,-600.,100,1.1,100.);
      TH2D *h38=new TH2D("","",100,-1600.,-600.,100,1.,100.);
      titles(h38, "", "Pulser height (m)", "r V_{#theta} (arb. units)");

      auto legend38 = new TLegend(0.62,0.2,0.88,0.5);
      legend38->SetBorderSize();
      legend38->SetTextSize(0.05);
      
      // h38->GetXaxis()->SetTitle("Pulser height (m)");
      // h38->GetYaxis()->SetTitle("r V (arb. units)");
      // h38->GetXaxis()->SetTitleOffset(0.8);
      h38->GetYaxis()->SetTitleOffset(1.35);
      // h38->GetXaxis()->SetTitleSize(0.09);
      //h38->GetYaxis()->SetTitleSize(1.1);
      //h38->GetXaxis()->SetTitleOffset(0.8);
      //h38->GetYaxis()->SetTitleOffset(0.8);
      //h38->GetXaxis()->SetTitleSize(0.05);
      //h38->GetYaxis()->SetTitleSize(0.05);
      //h38->GetXaxis()->SetNdivisions(504);
      h38->GetXaxis()->SetNdivisions(504);
      //gPad->SetTicks(2);
      //h38->GetYaxis()->SetTickLength(0.);
      h38->Draw();
      stext="{#theta}-pol";
      TText* texttemp_vpol=new TText(-1000.,160.,stext.c_str());
      texttemp_vpol->SetTextFont(42);
      //texttemp_vpol->Draw("same");
      TF1 *fy2=new TF1("fy2","x",0.,200.);
      TGaxis *A2 = new TGaxis(-600.,0.,-600.,200.,"fy2",510,"L+");
      A2->SetTitle("r V_{#theta} @ 300 MHz (arb. units)");
      //A2->SetTitleOffset(0.8);
      A2->SetTextFont(42);
      A2->SetLabelFont(42);
      A2->SetTitleSize(0.05);
      //A2->Draw("same");

      for (int istations=minstation;istations<=maxstation;istations++) {
	//if (istations!=5) {
	gvenvelope_plus_r1[istations]->SetLineStyle(kSolid);
	//	if (istations==5)


	gvenvelope_plus_r1[istations]->SetName(snames[istations].c_str());


	gvenvelope_minus_r1[istations]->SetLineStyle(kSolid);

	//if (istations!=5)
	

	gvoltage_r1[istations]->SetLineStyle(kDashed);

	if (istations==0) {
	  gvoltage_r1[istations]->SetName(snametemp.c_str());
	  
	}
	if (istations!=5) {
	  gvoltage_r1[istations]->Draw("lsame");
	  gvenvelope_minus_r1[istations]->Draw("lsame");
	  gvenvelope_plus_r1[istations]->Draw("lsame");
	}
	//	}

	legend38->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l"); 
     }

      legend38->Draw("same");

      sname=sdir+"HPolVPolvoltages_sidebyside.pdf";

      c37->Print(sname.c_str());

      //      TCanvas *c37=makePretty2Panel();
      TCanvas *c37b=new TCanvas("c37b","c37b",1600,800);
      //c37->SetLeftMargin(0.15);
      //c37->SetBottomMargin(0.15);
      c37b->Divide(2,1);
      c37b->cd(2);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetLogy();
      //string stext="{#phi}-pol";
      //TText *texttemp_hpol=new TText(-1000.,160.,stext.c_str());
      //texttemp_hpol->SetTextFont(42);
      //texttemp_hpol->SetTextSize(0.08);
      TH2D *h37b=new TH2D("","",100,-1600.,-600.,100,0.11,300.);
      titles(h37b, "", "Pulser height (m)", "r |E_{#phi}| (arb. units)");
      auto legend37b = new TLegend(0.22,0.52,0.36,0.9);
      legend37b->SetBorderSize();
      legend37b->SetTextSize(0.05);

      //h37->GetXaxis()->SetTitleOffset(0.8);
      h37b->GetYaxis()->SetTitleOffset(1.2);
      //h37->GetXaxis()->SetTitleSize(0.05);
      //h37->GetYaxis()->SetTitleSize(0.05);
      //h37->GetXaxis()->SetNdivisions(504);

      auto legend38c = new TLegend(0.17,0.76,0.62,0.89);
      legend38c->SetBorderSize();
      legend38c->SetTextSize(0.05);

      string snametempb="E field @ 300 MHz";
      h37b->Draw();
      //      texttemp_hpol->Draw("same");
      // h37->GetXaxis()->SetTitle("Pulser height (m)");
      // h37->GetYaxis()->SetTitle("r V_{#phi} (arb. units)");
      // h37->GetXaxis()->SetTitleOffset(0.8);
      // h37->GetYaxis()->SetTitleOffset(0.85);
      // h37->GetXaxis()->SetTitleSize(0.09);
      // h37->GetYaxis()->SetTitleSize(0.09);
      h37b->GetXaxis()->SetNdivisions(504);
      
      for (int istations=minstation;istations<=maxstation;istations++) {
	//if (istations!=5) {
	gEenvelope_plus_r2[istations]->SetLineStyle(kSolid);
	gEenvelope_plus_r2[istations]->SetLineColor(icolors[istations]);
	gEenvelope_plus_r2[istations]->SetLineWidth(2);
	//if (istations!=5)
	gEenvelope_plus_r2[istations]->Draw("lsame");	

	gEenvelope_minus_r2[istations]->SetLineStyle(kSolid);
	gEenvelope_minus_r2[istations]->SetLineColor(icolors[istations]);
	gEenvelope_minus_r2[istations]->SetLineWidth(2);
	//if (istations==5)
	gEenvelope_minus_r2[istations]->Draw("lsame");

	gfield_r2[istations]->SetLineStyle(kDashed);
	gfield_r2[istations]->SetLineColor(icolors[istations]);
	gfield_r2[istations]->SetLineWidth(2);
	gfield_r2[istations]->Draw("lsame");
	gfield_r2[istations]->SetName(snames[istations].c_str());
	legend37b->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l");
	//	}


	if (istations==0) {
	  legend38c->AddEntry(gEenvelope_plus_r2[istations],"E field envelope","l");
	  legend38c->AddEntry(gfield_r2[istations],snametempb.c_str(),"l");
	}

      }
      //      legend37->Draw("same");
      legend38c->Draw("same");


      
      c37b->cd(1);


 
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetLogy();
      //      TH2D *h38=new TH2D("","",100,-1600.,-600.,100,1.1,100.);
      TH2D *h38b=new TH2D("","",100,-1600.,-600.,100,1.,100.);
      titles(h38b, "", "Pulser height (m)", "r |E_{#theta}| (arb. units)");

      auto legend38d = new TLegend(0.62,0.2,0.88,0.5);
      legend38d->SetBorderSize();
      legend38d->SetTextSize(0.05);
      
      // h38->GetXaxis()->SetTitle("Pulser height (m)");
      // h38->GetYaxis()->SetTitle("r V (arb. units)");
      // h38->GetXaxis()->SetTitleOffset(0.8);
      h38b->GetYaxis()->SetTitleOffset(1.35);
      // h38->GetXaxis()->SetTitleSize(0.09);
      //h38->GetYaxis()->SetTitleSize(1.1);
      //h38->GetXaxis()->SetTitleOffset(0.8);
      //h38->GetYaxis()->SetTitleOffset(0.8);
      //h38->GetXaxis()->SetTitleSize(0.05);
      //h38->GetYaxis()->SetTitleSize(0.05);
      //h38->GetXaxis()->SetNdivisions(504);
      h38b->GetXaxis()->SetNdivisions(504);
      //gPad->SetTicks(2);
      //h38->GetYaxis()->SetTickLength(0.);
      h38b->Draw();
      stext="{#theta}-pol";
      //TText* texttemp_vpol=new TText(-1000.,160.,stext.c_str());
      //texttemp_vpol->SetTextFont(42);
      //texttemp_vpol->Draw("same");
      TF1 *fy2b=new TF1("fy2b","x",0.,200.);
      TGaxis *A2b = new TGaxis(-600.,0.,-600.,200.,"fy2b",510,"L+");
      A2b->SetTitle("r |E_{#theta}| @ 300 MHz (arb. units)");
      //A2->SetTitleOffset(0.8);
      A2b->SetTextFont(42);
      A2b->SetLabelFont(42);
      A2b->SetTitleSize(0.05);
      //A2->Draw("same");

      for (int istations=minstation;istations<=maxstation;istations++) {
	//if (istations!=5) {
	gEenvelope_plus_r1[istations]->SetLineColor(icolors[istations]);
	gEenvelope_plus_r1[istations]->SetLineWidth(2);
	gEenvelope_plus_r1[istations]->SetLineStyle(kSolid);

	gEenvelope_plus_r1[istations]->SetName(snames[istations].c_str());
	//	if (istations==5)
	gEenvelope_plus_r1[istations]->Draw("lsame");



	gEenvelope_minus_r1[istations]->SetLineColor(icolors[istations]);
	gEenvelope_minus_r1[istations]->SetLineWidth(2);
	gEenvelope_minus_r1[istations]->SetLineStyle(kSolid);

	//if (istations!=5)
	gEenvelope_minus_r1[istations]->Draw("lsame");

	gfield_r1[istations]->SetLineStyle(kDashed);
	gfield_r1[istations]->SetLineColor(icolors[istations]);
	gfield_r1[istations]->SetLineWidth(2);

	if (istations==0) {
	  gfield_r1[istations]->SetName(snametemp.c_str());
	  
	}
	gfield_r1[istations]->Draw("lsame");
	
	//	}

	legend38d->AddEntry(snames[istations].c_str(),snames[istations].c_str(),"l"); 
     }

      legend38d->Draw("same");

      sname=sdir+"HPolVPolfields_sidebyside.pdf";

      c37b->Print(sname.c_str());


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

     char name[100];
     //     sprintf(name,"myoutputs_%dMHz_tx%d_rx%d.root",(int)(freq/1.E6),CROSSPOLANGLE_TX_INT,CROSSPOLANGLE_RX_INT);
     
     string sfilenames=sdir+"myoutputs_tx" + to_string(CROSSPOLANGLE_TX_INT) + "_rx" + to_string(CROSSPOLANGLE_RX_INT) + ".root";

     //     sprintf(name,"myoutputs_tx%d_rx%d.root",CROSSPOLANGLE_TX_INT,CROSSPOLANGLE_RX_INT);
     TFile *fout=new TFile(sfilenames.c_str(),"RECREATE");

     for (int istations=0;istations<NSTATIONS;istations++) {
       sprintf(name,"gpolarization_reversedepth_Psi_rx_%d",istations);
       gpolarization_reversedepth_Psi_rx[istations]->Write(name);
       sprintf(name,"gEpolarization_reversedepth_Psi_rx_%d",istations);
       gEpolarization_reversedepth_Psi_rx[istations]->Write(name);
       sprintf(name,"gvoltage_r1_%d",istations);
       gvoltage_r1[istations]->Write(name);
       sprintf(name,"gvoltage_r2_%d",istations);
       gvoltage_r2[istations]->Write(name);
       sprintf(name,"grxdepth_atten_%d",istations);
       grxdepth_atten[istations]->Write(name);
       sprintf(name,"grxdepth_atten_beam_%d",istations);
       grxdepth_atten_beam[istations]->Write(name);
       sprintf(name,"gvenvelope_plus_r1_%d",istations);
       gvenvelope_plus_r1[istations]->Write(name);
       //       if (istations==5) {
       for (int ispecial=0;ispecial<NSPECIAL;ispecial++) {
	 sprintf(name,"g_spectra_%d_%d",istations,ispecial);	 
	 
	 g_spectra[istations][(int)g_idepth[istations]->Eval(whichspecial[ispecial])]->Write(name);
       }
	 //}
       sprintf(name,"gvenvelope_plus_r2_%d",istations);
       gvenvelope_plus_r2[istations]->Write(name);

       sprintf(name,"gvoltage_r2_%d",istations);
       gvoltage_r2[istations]->Write(name);

       sprintf(name,"gfield_r2_%d",istations);
       gfield_r2[istations]->Write(name);

       sprintf(name,"g_atten_power_%d",istations);
       g_atten_power[istations]->Write(name);

       sprintf(name,"g_atten_beam_power_%d",istations);
       g_atten_beam_power[istations]->Write(name);

       sprintf(name,"gV1squared_r1_%d",istations);
       gV1squared_r1[istations]->Write(name);

       sprintf(name,"gV2squared_r1_%d",istations);
       gV2squared_r1[istations]->Write(name);

       sprintf(name,"gV1V2_r1_%d",istations);
       gV1V2_r1[istations]->Write(name);

       sprintf(name,"gV1V2_r2_%d",istations);
       gV1V2_r2[istations]->Write(name);

       sprintf(name,"gV1squared_r2_%d",istations);
       gV1squared_r2[istations]->Write(name);

       sprintf(name,"gV2squared_r2_%d",istations);
       gV2squared_r2[istations]->Write(name);

       sprintf(name,"goppositeV1V2_r2_%d",istations);
       goppositeV1V2_r2[istations]->Write(name);

       sprintf(name,"goppositeV1V2_r1_%d",istations);
       goppositeV1V2_r1[istations]->Write(name);
       
       int idepth_temp=g_idepth[istations]->Eval(-1000.);
       cout << "station, depth, V1squared_r2, V2squared_r2, V1V2_r2 are " << istations << "\t" << vV1squared_r2[istations][idepth_temp] << "\t" << vV2squared_r2[istations][idepth_temp] << "\t" << vV1V2_r2[istations][idepth_temp] << "\n"; 

     }
     fout->Close();

     



  return (0);

}

TVector3 directionNextStep(vector<double>nvec, TVector3 p_e1, TVector3 p_e2) {

  p_e1.SetMag(1.);

  TVector3 temp1;
  temp1[0]=nvec[0]*p_e1[0];
  temp1[1]=nvec[1]*p_e1[1];
  temp1[2]=nvec[2]*p_e1[2];
  temp1.SetMag(1.);

  

  
  return temp1;
}
/*
TVector3 rotateD(TVector3 epsilon, double angle_iceflow, TVector3 D) {

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

  TVector3 inverseepsilon(1./epsilon[0],1./epsilon[1],1./epsilon[2]);

  for (int i=0;i<3;i++) {
    tempvec[i]=inverseepsilon[i]*D[i];
  }
  D=tempvec;

  double rotate_backtonormal[3][3];
  
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
    }    
  }

  for (int i=0;i<3;i++) {
    double sum=0.;
    
    for (int j=0;j<3;j++) {
      sum+=rotate_backtonormal[i][j]*D[j];
    }

    tempvec[i]=sum;
   }

  D=tempvec;

  return D; // this is actually returning an electric field


}
*/
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

  //  double epsilon_perp_c= (nvec[0]*nvec[0]+nvec[1]*nvec[1]+nvec[2]*nvec[2])/3.;

  //  cout << "input D is " << D[0] << "\t" << D[1] << "\t" << D[2] << "\n";

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

  // TVector3 D;

  // D[0]=epsilon[0]*Pt[0];
  // D[1]=epsilon[1]*Pt[1];
  // D[2]=epsilon[2]*Pt[2];

  E[0]=epsilon_inverse[0]*D[0];
  E[1]=epsilon_inverse[1]*D[1];
  E[2]=epsilon_inverse[2]*D[2];

  
  // cout << "input D is " << D[0] << "\t" << D[1] << "\t" << D[2] << "\n";
  // cout << "epsilon is " << epsilon[0] << "\t" << epsilon[1] << "\t" << epsilon[2] << "\n";
  // cout << "E is " << E[0] << "\t" << E[1] << "\t" << E[2] << "\n";
  // cout << "n*n*E is " << n*n*E[0] << "\t" << n*n*E[1] << "\t" << n*n*E[2] << "\n";


  TVector3 kvec=-1.*D+n*n*E;
  //cout << "kvec dot E is " << kvec.Dot(E) << "\n";


  //  if (kvec.Dot(E)<0.)
  //kvec=-1.*kvec;


  if (kvec.Mag()<HOWSMALLISTOOSMALL)
    cout << "kvec mag is " << kvec.Mag() << "\n";

  kvec.SetMag(1.);


  double rotate_backtonormal[3][3];

  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
    }

  }




  //  cout << "kvec is " << kvec[0] << "\t" << kvec[1] << "\t" << kvec[2] << "\n";

  
  //cout << "kvec dot E is " << kvec.Dot(Pt) << "\n";

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
    //E=-1.*E;
  }

  return kvec;



}
/*
double getDeltaN(int BIAXIAL,vector<double> nvec,TVector3 rhat,double angle_iceflow, double &n_e1, double &n_e2,TVector3 &p_e1,TVector3 &p_e2) {

  //  for testing
  //  rhat[0]=0.001*cos(angle_iceflow);
  //rhat[1]=0.001*sin(angle_iceflow);
  //rhat[2]=sqrt(1.-0.001*0.001);

  //  rhat[0]=0.9999999*cos(angle_iceflow);
  //rhat[1]=0.9999999*sin(angle_iceflow);
  //rhat[2]=sqrt(1.-0.9999999*0.9999999);

  //  rhat[0]=0.999999;
  //rhat[1]=0.;
  //rhat[2]=sqrt(1.-0.999999*0.999999);

  int FLIPPED=0;

  if (rhat[2]<0.) {
    FLIPPED=1;
    rhat[2]=-1.*rhat[2];
  }




  TVector3 myy;
  myy[0]=0.;
  myy[1]=-1.;
  myy[2]=0.;

  TVector3 myz;
  myz[0]=0.;
  myz[1]=0.;
  myz[2]=1.;
 
  //  rhat[0]=1/sqrt(2);
  //rhat[1]=0.;
  //rhat[2]=1/sqrt(2.);


  double phi_rhat=atan2(rhat[1],rhat[0]);
  double phi_wrticeflow=angle_iceflow-phi_rhat;
  if (phi_wrticeflow<-PI)
    phi_wrticeflow+=2.*PI;
  if (phi_wrticeflow>PI)
    phi_wrticeflow-=2.*PI;


  double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
				       {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
				       {0.,0.,1.}};

  TVector3 rhat_iceflowalongx;
  TVector3 myy_iceflowalongx;
  TVector3 myz_iceflowalongx;
  TVector3 x_iceflowalongx(1.,0.,0.);
  TVector3 y_iceflowalongx(0.,1.,0.);
  TVector3 z_iceflowalongx(0.,0.,1.);

  //  cout << "rhat is " << rhat[0] << "\t" << rhat[1] << "\t" << rhat[2] << "\n";
  for (int i=0;i<3;i++) {
    double sum=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_toxalongiceflow[i][j]*rhat[j];
      }
    rhat_iceflowalongx[i]=sum;
  }



  TVector3 nominal_pe1=myz.Cross(rhat_iceflowalongx);
  


  if (nominal_pe1.Mag()<HOWSMALLISTOOSMALL) {
    //cout << "myz is " << myz[0] << "\t" << myz[1] << "\t" << myz[2] << "\n";
    //cout << "rhat_iceflowalongx is " << rhat_iceflowalongx[0] << "\t" << rhat_iceflowalongx[1] << "\t" << rhat_iceflowalongx[2] << "\n";
    //cout << "rhat mag is " << rhat.Mag() << "\n";

    cout << "myz is " << myz[0] << "\t" << myz[1] << "\t" << myz[2] << "\n";
    cout << "rhat is " << rhat[0] << "\t" << rhat[1] << "\t" << rhat[2] << "\n";
    cout << "rhat_iceflowalongx is " << rhat_iceflowalongx[0] << "\t" << rhat_iceflowalongx[1] << "\t" << rhat_iceflowalongx[2] << "\n";
    cout << "cross of them is " << nominal_pe1[0] << "\t" << nominal_pe1[1] << "\t" << nominal_pe1[2] << "\n";
    cout << "nominal_pe1 mag is " << nominal_pe1.Mag() << "\n";
  }
  nominal_pe1.SetMag(1.);
  TVector3 nominal_pe2=rhat_iceflowalongx.Cross(nominal_pe1);

  double theta_nominal=atan2(nominal_pe1[1],nominal_pe1[0]);


  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_toxalongiceflow[i][j]*myy[j];
      sum2+=rotate_toxalongiceflow[i][j]*myz[j];
      }
    myy_iceflowalongx[i]=sum;
    myz_iceflowalongx[i]=sum2;
  }




  double a=nvec[0];
  double b=nvec[1];
  double c=nvec[2];
  

  
  //double Ax=sqrt(1.-rhat[2]*rhat[2])*cos(phi_wrticeflow);
  //double Ay=sqrt(1.-rhat[2]*rhat[2])*sin(phi_wrticeflow);
  //double Az=rhat[2];

  double Ax=rhat_iceflowalongx[0];
  double Ay=rhat_iceflowalongx[1];
  double Az=rhat_iceflowalongx[2];

  //  cout << "At beginning of getDeltaN, rhat is " << rhat[0] << "\t" << rhat[1] << "\t" << rhat[2] << "\n";
  //cout << "a, b, c are " << a << "\t" << b << "\t" << c << "\n";
  //cout << "phi_wrticeflow is " << phi_wrticeflow << "\n";
  //cout << "Ax, Ay, Az are " << Ax << "\t" << Ay << "\t" << Az << "\n";

  double A=1/(a*a)+(Ax*Ax)/(Az*Az*c*c);
  double B=(2.*Ax*Ay)/(Az*Az*c*c);
  double C=1/(b*b)+(Ay*Ay)/(Az*Az*c*c);
  double F=-1.;

  double M=A;
  double P=B/2.;
  double Q=B/2.;
  double R=C;

  double theta_initial=atan2( B , A - C )/2.; // not sure this is rotated in the right direction - check this.

  //  cout << "theta_initial is " << theta_initial << "\n";


  

  double lambda2=(1.*(M+R)+sqrt((M-R)*(M-R)+4*P*Q))/2.;
  double lambda1=(1.*(M+R)-sqrt((M-R)*(M-R)+4*P*Q))/2.;
  
  // these are only the n's for the scenario where the plane arrives straight from above, a test scenario
  double ne2=sqrt(-1./lambda2);
  double ne1=sqrt(-1./lambda1);



  double Psi=PI/2.-atan2(abs(Az),sqrt(Ax*Ax+Ay*Ay));
  double omega=-1.*(PI - atan2(Ay,Ax));
  double epsilon=0.;
  
  //  cout << "Ax, Ay, atan2, omega are " << Ax << "\t" << Ay << "\t" << atan2(Ay,Ax) << "\t" << omega << "\n";
  //cout << "rhat[2], acos, Psi, omega are " << rhat[2] << "\t" << acos(rhat[2]) << "\t" << Psi << "\t" << omega << "\n";

  //    double rotate[3][3]={{cos(Psi)*cos(omega),sin(omega),-sin(Psi)*cos(omega)},
  //{-1.*cos(Psi)*sin(omega),cos(omega),sin(Psi)*sin(omega)},
  //		       {sin(Psi),0.,cos(Psi)}};

  //  double rotate[3][3]={{cos(Psi)*cos(omega),sin(omega),sin(Psi)*cos(omega)},
  //		       {-1.*cos(Psi)*sin(omega),cos(omega),-1.*sin(Psi)*sin(omega)},
  //		       {-1.*sin(Psi),0.,cos(Psi)}};
  
  //  double rotate[3][3]={{1.*cos(omega),-1.*sin(omega),0.},
  //		       {1.*sin(omega),cos(omega),0.},
  //		       {0.,0.,1.}};

  double rotate[3][3]={{cos(Psi)*cos(omega),cos(Psi)*sin(omega),sin(Psi)},
  		       {-1.*sin(omega),cos(omega),0.},
  		       {-1.*sin(Psi)*cos(omega),-1.*sin(Psi)*sin(omega),cos(Psi)}};

  // cout << "rotate is \n";
  // for (int i=0;i<3;i++) {
  //   for (int j=0;j<3;j++) {

  //     cout << rotate[i][j] << "\t";

  //   }
  //   cout << "\n";
  // }

  
  double myy_rotate[3];
  double myz_rotate[3];

  TVector3 nominal_pe1_rotate;
  TVector3 nominal_pe2_rotate;

  TVector3 x_iceflowalongx_rotate;
  TVector3 y_iceflowalongx_rotate;

  TVector3 rhat_rotate;
  TVector3 tmpvec;
  TVector3 tmpvec2;
  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum1=0.;
    double sum2=0.;
    double sum3=0.;
    double sum4=0.;
    double sum5=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate[i][j]*rhat_iceflowalongx[j];
      sum1+=rotate[i][j]*nominal_pe1[j];
      sum2+=rotate[i][j]*nominal_pe2[j];
      sum3+=rotate[i][j]*x_iceflowalongx[j];
      sum4+=rotate[i][j]*y_iceflowalongx[j];
      //sum5+=rotate[i][j]*findthefreakingaxis[j];
    }
    rhat_rotate[i]=sum;
    nominal_pe1_rotate[i]=sum1;
    nominal_pe2_rotate[i]=sum2;
    x_iceflowalongx_rotate[i]=sum3;
    y_iceflowalongx_rotate[i]=sum4;
    //tmpvec[i]=sum5;
  }
  //  findthefreakingaxis=tmpvec;

  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate[i][j]*myy_iceflowalongx[j];
      sum2+=rotate[i][j]*myz_iceflowalongx[j];
      }
    myy_rotate[i]=sum;
    myz_rotate[i]=sum2;
  }

  //  cout << "after rotating, findthefreakingaxis is " << findthefreakingaxis[0] << "\t" << findthefreakingaxis[1] << "\t" << findthefreakingaxis[2] << "\n";









  // cout << "rhat is ";

  // for (int i=0;i<3;i++) {
  //   cout << rhat[i] << "\t";
  // }
  // cout << "\n";




  // cout << "rhat_iceflowalongx is ";
  // for (int i=0;i<3;i++) {
  //   cout << rhat_iceflowalongx[i] << "\t";
  // }
  // cout << "\n";

  // cout << "nominal_pe1 is ";
  // for (int i=0;i<3;i++) {
  //   cout << nominal_pe1[i] << "\t";
  // }
  // cout << "\n";

  // cout << "nominal_pe2 is ";
  // for (int i=0;i<3;i++) {
  //   cout << nominal_pe2[i] << "\t";
  // }
  // cout << "\n";

  // cout << "theta_nominal is " << theta_nominal << "\n";


  // cout << "myy_iceflowalongx is ";
  // for (int i=0;i<3;i++) {
  //   cout << myy_iceflowalongx[i] << "\t";
  // }
  // cout << "\n";

  // cout << "myz_iceflowalongx is ";
  // for (int i=0;i<3;i++) {
  //   cout << myz_iceflowalongx[i] << "\t";
  // }
  // cout << "\n";

  // cout << "rhat_rotate is ";
  // for (int i=0;i<3;i++) {
  //   cout << rhat_rotate[i] << "\t";
  // }
  // cout << "\n";

  // cout << "myy_rotate is ";
  // for (int i=0;i<3;i++) {
  //   cout << myy_rotate[i] << "\t";
  // }
  // cout << "\n";


  // cout << "myz_rotate is ";
  // for (int i=0;i<3;i++) {
  //   cout << myz_rotate[i] << "\t";
  // }
  // cout << "\n";

  // cout << "nominal_pe1_rotate is ";
  // for (int i=0;i<3;i++) {
  //   cout << nominal_pe1_rotate[i] << "\t";
  // }
  // cout << "\n";


  // cout << "nominal_pe2_rotate is ";
  // for (int i=0;i<3;i++) {
  //   cout << nominal_pe2_rotate[i] << "\t";
  // }
  // cout << "\n";

  // cout << "x_iceflowalongx_rotate is ";
  // for (int i=0;i<3;i++) {
  //   cout << x_iceflowalongx_rotate[i] << "\t";
  // }
  // cout << "\n";

  // cout << "y_iceflowalongx_rotate is ";
  // for (int i=0;i<3;i++) {
  //   cout << y_iceflowalongx_rotate[i] << "\t";
  // }
  // cout << "\n";



  //  double Psi_T=-1.*asin(rotate[0][2]);
  //double omega_T=-1.*atan2(rotate[1][0],rotate[1][1]);
  //double Psi_T=-1.*Psi;
  //double omega_T=1.*omega;
  //double omega_T=0.;

  //  cout << "rotate[0][2] is " << rotate[0][2] << "\n";

  //  double rotate_T[3][3]={{cos(Psi_T)*cos(omega_T),sin(omega_T),-sin(Psi_T)*cos(omega_T)},
  //		       {-1.*cos(Psi_T)*sin(omega_T),cos(omega_T),sin(Psi_T)*sin(omega_T)},
  //		       {sin(Psi_T),0.,cos(Psi_T)}};
  
  //  double rotate_T[3][3]={{cos(Psi_T)*cos(omega_T),sin(omega_T),sin(Psi_T)*cos(omega_T)},
  //		       {-1.*cos(Psi_T)*sin(omega_T),cos(omega_T),-1.*sin(Psi_T)*sin(omega_T)},
  //		       {-1.*sin(Psi_T),0.,cos(Psi_T)}};

  double rotate_T[3][3]={{0.}};
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_T[i][j]=rotate[j][i];
    }
  }

  //  double epsilon_T=-1.*atan2(rotate_T[2][1],rotate_T[2][2]);
  //double Psi_T= asin(rotate_T[2][0]);
  //double omega_T=-1.*atan2(rotate_T[1][0],rotate_T[0][0]);

  double epsilon_T=-1.*atan2(rotate_T[1][2],rotate_T[2][2]);
  double Psi_T= asin(rotate_T[0][2]);
  double omega_T=-1.*atan2(rotate_T[0][1],rotate_T[0][0]);
  
  //epsilon_T=-1.*epsilon_T;
  //Psi_T=-1.*Psi_T;
  //omega_T=-1.*omega_T;


  double rotate_T_usingformula[3][3]={{cos(Psi_T)*cos(omega_T),cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T),sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T)},
				      {-1.*cos(Psi_T)*sin(omega_T),cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T),sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T)},
				      {sin(Psi_T),-1.*sin(epsilon_T)*cos(Psi_T),cos(epsilon_T)*cos(Psi_T)}};






//   cout << "Psi_T, omega_T are " << Psi_T << "\t" << omega_T << "\t" << epsilon_T << "\n";
//   cout << "rotate_T is \n";
//   for (int i=0;i<3;i++) {
//     for (int j=0;j<3;j++) {

//       cout << rotate_T[i][j] << "\t";

//     }
//     cout << "\n";
//   }

// cout << "rotate_T_usingformula is \n";
//   for (int i=0;i<3;i++) {
//     for (int j=0;j<3;j++) {

//       cout << rotate_T_usingformula[i][j] << "\t";

//     }
//     cout << "\n";
//   }

  TVector3 rhat_rotateback;
  TVector3 myy_rotateback;
  TVector3 myz_rotateback;

  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    double sum3=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_T[i][j]*rhat_rotate[j];
      sum2+=rotate_T[i][j]*myy_rotate[j];
      sum3+=rotate_T[i][j]*myz_rotate[j];
    }
    rhat_rotateback[i]=sum;
    myy_rotateback[i]=sum2;
    myz_rotateback[i]=sum3;
  }
  // cout << "rhat_rotateback is ";
  // for (int i=0;i<3;i++) {
  //   cout << rhat_rotateback[i] << "\t";
  // }
  // cout << "\n";
  // cout << "myy_rotateback is ";
  // for (int i=0;i<3;i++) {
  //   cout << myy_rotateback[i] << "\t";
  // }
  // cout << "\n";
  // cout << "myz_rotateback is ";
  // for (int i=0;i<3;i++) {
  //   cout << myz_rotateback[i] << "\t";
  // }
  // cout << "\n";

  double Anew=1/(a*a)*(cos(Psi_T)*cos(Psi_T)*cos(omega_T)*cos(omega_T)) +
    1/(b*b)*pow(cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T),2) +
    1/(c*c)*pow(sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T),2);

  double Bnew=1/(a*a)*(-2.*cos(Psi_T)*cos(Psi_T)*cos(omega_T)*sin(omega_T)) +
    1/(b*b)*2.*(cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T))*(cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T)) +
    1/(c*c)*2.*(sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T))*(sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T));

  double Cnew=1/(a*a)*(cos(Psi_T)*cos(Psi_T)*sin(omega_T)*sin(omega_T)) +
    1/(b*b)*pow(cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T),2) +
    1/(c*c)*pow(sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T),2);

  double Dnew=0.;
  double Enew=0.;
  

  //  cout << "Anew, Bnew, Cnew are " << Anew << "\t" << Bnew << "\t" << Cnew << "\n";

  double Fnew=-1.;

  double Mnew=Anew;
  double Pnew=Bnew/2.;
  double Qnew=Bnew/2.;
  double Rnew=Cnew;
  

  double lambda2_new=(1.*(Mnew+Rnew)+sqrt((Mnew-Rnew)*(Mnew-Rnew)+4*Pnew*Qnew))/2.;
  double lambda1_new=(1.*(Mnew+Rnew)-sqrt((Mnew-Rnew)*(Mnew-Rnew)+4*Pnew*Qnew))/2.;
  

  double ne1new=sqrt(1./lambda2_new);
  double ne2new=sqrt(1./lambda1_new);

  //  cout << "lambda1, lambda2 are " << lambda1_new << "\t" << lambda2_new << "\n";
  //cout << "ne1new, ne2new are " << ne1new << "\t" << ne2new << "\n";

  p_e1[0]=1.;
  p_e1[1]=0.;
  p_e1[2]=0.;

  TVector3 rhat_rotatetheta=rhat_rotate;

  //  p_e2[0]=0.;
  //p_e2[1]=1.;
  //p_e2[2]=0.;
  
  //  double theta=atan2( Bnew , Anew - Cnew )/2.-theta_initial; // not sure this is rotated in the right direction - check this.
  double theta2=atan2( Bnew , Anew - Cnew )/2.; // not sure this is rotated in the right direction - check this.
  //double theta=0.; // not sure this is rotated in the right direction - check this.


  TVector3 findthefreakingaxis(1.,0.,0.);
  TVector3 findthefreakingaxis_perp=findthefreakingaxis;

  findthefreakingaxis.RotateZ(theta2);
  findthefreakingaxis_perp.RotateZ(theta2+PI/2.);

  // cout << "theta2 is " << theta2 << "\n";
  // cout << "after rotating by theta2, findthefreakingaxis is " << findthefreakingaxis[0] << "\t" << findthefreakingaxis[1] << "\t" << findthefreakingaxis[2] << "\n";
  // cout << "after rotating by theta2, findthefreakingaxis_perp is " << findthefreakingaxis_perp[0] << "\t" << findthefreakingaxis_perp[1] << "\t" << findthefreakingaxis_perp[2] << "\n";
  // cout << "rhat_rotate is " << rhat_rotate[0] << "\t" << rhat_rotate[1] << "\t" << rhat_rotate[2] << "\n";

  TVector3 rhat_unrotate=rhat_rotate;

  TVector3 tmpvec3;

  for (int i=0;i<3;i++) {
    double sum1=0.;
    double sum2=0.;
    double sum3=0.;
    for (int j=0;j<3;j++) {
      sum1+=rotate_T[i][j]*findthefreakingaxis[j];
      sum2+=rotate_T[i][j]*findthefreakingaxis_perp[j];
      sum3+=rotate_T[i][j]*rhat_rotate[j];
    }
    tmpvec[i]=sum1;
    tmpvec2[i]=sum2;
    tmpvec3[i]=sum3;
  }
  findthefreakingaxis=tmpvec;
  findthefreakingaxis_perp=tmpvec2;
  rhat_unrotate=tmpvec3;
  // this should be the +y axis in the unrotated coordinate system (iceflowalongx) rotated by theta_initial.  find the difference between the angle this makes with the +y axis and theta_initial and rotate it by that difference (I the uniaxial case this difference might always be zero.)
  // cout << "after unrotating, findthefreakingaxis is " << findthefreakingaxis[0] << "\t" << findthefreakingaxis[1] << "\t" << findthefreakingaxis[2] << "\n";
  // cout << "after unrotating, findthefreakingaxis_perp is " << findthefreakingaxis_perp[0] << "\t" << findthefreakingaxis_perp[1] << "\t" << findthefreakingaxis_perp[2] << "\n";
  // cout << "after unrotating, rhat_unrotate is " << rhat_unrotate[0] << "\t" << rhat_unrotate[1] << "\t" << rhat_unrotate[2] << "\n";

  TVector3 findthefreakingaxis_projecttoXY(findthefreakingaxis[0],findthefreakingaxis[1],0.);
  TVector3 findthefreakingaxis_perp_projecttoXY(findthefreakingaxis_perp[0],findthefreakingaxis_perp[1],0.);
  TVector3 yaxis(0.,1.,0.);

  double anglebetweenthem=findthefreakingaxis_projecttoXY.Angle(yaxis);



  double diffangle=theta_initial-anglebetweenthem;
  // cout << "theta_inital is " << theta_initial << "\n";
  // cout << "angle between findthefreakingaxis_projectXY and yaxis is " << anglebetweenthem << "\n";
  // cout << "diffangle is " << diffangle << "\n";

  diffangle=0.;
  findthefreakingaxis_projecttoXY.RotateZ(diffangle);
  findthefreakingaxis_perp_projecttoXY.RotateZ(diffangle);
  rhat_unrotate.RotateZ(diffangle);

  findthefreakingaxis[0]=findthefreakingaxis_projecttoXY[0];
  findthefreakingaxis[1]=findthefreakingaxis_projecttoXY[1];

  findthefreakingaxis_perp[0]=findthefreakingaxis_perp_projecttoXY[0];
  findthefreakingaxis_perp[1]=findthefreakingaxis_perp_projecttoXY[1];


  // cout << "after possible correction, findthefreakingaxis is " << findthefreakingaxis[0] << "\t" << findthefreakingaxis[1] << "\t" << findthefreakingaxis[2] << "\n";
  // cout << "after possible correction, findthefreakingaxis_perp is " << findthefreakingaxis_perp[0] << "\t" << findthefreakingaxis_perp[1] << "\t" << findthefreakingaxis_perp[2] << "\n";
  // cout << "after possible correction, rhat_rotate is " << rhat_unrotate[0] << "\t" << rhat_unrotate[1] << "\t" << rhat_unrotate[2] << "\n";
  
  // cout << "looks good here.\n";

  // cout << "findthefreakingaxis dot rhat_unrotate is " << findthefreakingaxis.Dot(rhat_unrotate) << "\n";
  // cout << "rhat_unrotate dot findthefreakingaxis_perp is " << rhat_unrotate.Dot(findthefreakingaxis_perp) << "\n";
  // cout << "findthefreakingaxis_perp dot findthefreakingaxis is " << findthefreakingaxis_perp.Dot(findthefreakingaxis) << "\n";

  // double rotate_theta[3][3]={ { cos(theta2) , -1.*sin(theta2) , 0.     },
  // 			      { 1.*sin(theta2) , cos(theta2), 0.     },
  // 			      { 0. , 0. , 1.     }  };

  // TVector3 p_e1_tmp;
  // TVector3 p_e2_tmp;

  // for (int i=0;i<3;i++) {
  //   double sum1=0.;
  //   double sum2=0.;

  //   for (int j=0;j<3;j++) {
  //     sum1+=rotate_theta[i][j]*p_e1[j];
  //     sum2+=rotate_theta[i][j]*rhat_rotate[j];
      
  //   }
  //   p_e1_tmp[i]=sum1;
  //   rhat_rotatetheta[i]=sum2;
  // }

  // p_e1=p_e1_tmp;
  // p_e2=rhat_rotate.Cross(p_e1);

  // //  cout << "after theta rotation, p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
  // //cout << "after theta rotation, p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";
  // // cout << "after theta rotation, rhat_rotatetheta is " << rhat_rotatetheta[0] << "\t" << rhat_rotatetheta[1] << "\t" << rhat_rotatetheta[2] << "\n";





  // cout << "phi_wrticeflow is " << phi_wrticeflow << "\n";


  // TVector3 rhat_rotateawayfromiceflow=rhat_rotatetheta;
  TVector3 rhat_rotateawayfromiceflow;


  double rotate_backtonormal[3][3];

  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
    }

  }

  //  TVector3 tmpvec;
  
  for (int i=0;i<3;i++) {
    double sum1=0.;
    double sum2=0.;
    double sum3=0.;
    double sum4=0.;

    for (int j=0;j<3;j++) {
      sum1+=rotate_T[i][j]*p_e1[j];
      sum2+=rotate_T[i][j]*p_e2[j];
      sum3+=rotate_T[i][j]*rhat_rotatetheta[j];
      //sum4+=rotate_T[i][j]*findthefreakingaxis[j];
    }
    p_e1_tmp[i]=sum1;
    p_e2_tmp[i]=sum2;
    rhat_unrotate[i]=sum3;
    //tmpvec[i]=sum4;
  }

  p_e1=p_e1_tmp;
  p_e2=p_e2_tmp;
  

  //findthefreakingaxis=tmpvec;

  //  cout << "after rotate away from iceflow, rhat_rotateawayfromiceflow is " << rhat_rotateawayfromiceflow[0] << "\t" << rhat_rotateawayfromiceflow[1] << "\t" << rhat_rotateawayfromiceflow[2] << "\n";
  //cout << "after unrotating, rhat_unrotate is " << rhat_unrotate[0] << "\t" << rhat_unrotate[1] << "\t" << rhat_unrotate[2] << "\n";
  //cout << "after unrotating, p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
  //cout << "after unrotating, p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";
  //cout << "after unrotating, findthefreakingaxis is " << findthefreakingaxis[0] << "\t" << findthefreakingaxis[1] << "\t" << findthefreakingaxis[2] << "\n";
  //cout << "after unrotating, findthefreakingaxis_perp is " << findthefreakingaxis_perp[0] << "\t" << findthefreakingaxis_perp[1] << "\t" << findthefreakingaxis_perp[2] << "\n";
  
  p_e1=findthefreakingaxis;
  p_e2=findthefreakingaxis_perp;




  //  cout << "p_e1 dot rhat_unrotate is " << p_e1.Dot(rhat_unrotate) << "\n";
  //cout << "rhat_unrotate dot p_e2 is " << rhat_unrotate.Dot(p_e2) << "\n";
  //cout << "p_e2 dot p_e1 is " << p_e2.Dot(p_e1) << "\n";

/*
  for (int i=0;i<3;i++) {
    double sum1=0.;
    double sum2=0.;
    for (int j=0;j<3;j++) {
      sum1+=rotate_T[i][j]*p_e1[j];
      sum2+=rotate_T[i][j]*p_e2[j];
      //sum1+=rotate[i][j]*p_e1[j];
      //sum2+=rotate[i][j]*p_e2[j];
    }
    p_e1_tmp[i]=sum1;
    p_e2_tmp[i]=sum2;
  }
  p_e1=p_e1_tmp;
  p_e2=p_e2_tmp;

// p_e1 is in 1st or 4th quadrant in coordinate system where
// ice is along the x axis
  if (!(p_e1.Phi()>-1*PI/2. && p_e1.Phi()<PI/2.))
    p_e1=-1.*p_e1;
  // p_e1 crossed with p_e2 should be in the vertical z direction
  if ((p_e1.Cross(p_e2)).Dot(myz)<0.)
    p_e2=-1.*p_e2;

 

  			    
  for (int i=0;i<3;i++) {
    double sum1=0.;
    double sum2=0.;
    double sum3=0.;
    double sum4=0.;
    double sum5=0.;
    for (int j=0;j<3;j++) {
      //      sum1+=rotate_backtonormal[i][j]*p_e1[j];
      //sum2+=rotate_backtonormal[i][j]*p_e2[j];
      sum3+=rotate_backtonormal[i][j]*rhat_unrotate[j];
      sum4+=rotate_backtonormal[i][j]*findthefreakingaxis[j];
      sum5+=rotate_backtonormal[i][j]*findthefreakingaxis_perp[j];
    }
    //    p_e1_tmp[i]=sum1;
    //p_e2_tmp[i]=sum2;
    rhat_rotateawayfromiceflow[i]=sum3;
    tmpvec[i]=sum4;
    tmpvec2[i]=sum5;
   }


  findthefreakingaxis=tmpvec;
  findthefreakingaxis_perp=tmpvec2;

  //  p_e1=p_e1_tmp;
  //p_e2=p_e2_tmp;

  //  cout << "after rotate away from iceflow, p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
  //cout << "after rotate away from iceflow, p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";

  // cout << "doesn't look good here.\n";
  // cout << "after rotate away from iceflow, findthefreakingaxis is " << findthefreakingaxis[0] << "\t" << findthefreakingaxis[1] << "\t" << findthefreakingaxis[2] << "\n";
  // cout << "after rotate away from iceflow, findthefreakingaxis_perp is " << findthefreakingaxis_perp[0] << "\t" << findthefreakingaxis_perp[1] << "\t" << findthefreakingaxis_perp[2] << "\n";
  // cout << "after rotate away from iceflow, rhat_rotateawayfromiceflow is " << rhat_rotateawayfromiceflow[0] << "\t" << rhat_rotateawayfromiceflow[1] << "\t" << rhat_rotateawayfromiceflow[2] << "\n";
  // cout << "remember input rhat was " << rhat[0] << "\t" << rhat[1] << "\t" << rhat[2] << "\n";

  p_e1=findthefreakingaxis;
  p_e2=findthefreakingaxis_perp;


  if (FLIPPED==1) {
    p_e1[2]=-1.*p_e1[2];
    p_e2[2]=-1.*p_e2[2];
    rhat_rotateawayfromiceflow[2]=-1.*rhat_rotateawayfromiceflow[2];
  }

  // cout << "p_e1 dot rhat is " << p_e1.Dot(rhat_rotateawayfromiceflow) << "\n";
  // cout << "rhat dot p_e2 is " << rhat_rotateawayfromiceflow.Dot(p_e2) << "\n";
  // cout << "p_e2 dot p_e1 is " << p_e2.Dot(p_e1) << "\n";

  n_e1=ne1new;
  n_e2=ne2new;

  // ray 1 is the one with the shortest index of refraction.
  if (n_e2<n_e1) {

    double n_e1_temp=n_e1;
    n_e1=n_e2;
    n_e2=n_e1_temp;
    
    TVector3 p_e1_temp=p_e1;
    p_e1=p_e2;
    p_e2=p_e1_temp;

  }

  // for an isotropic medium it can tend to pick the orientation of the axes
  // somewhat randomly.
  // this is why when an isotropic medium is chosen, I pick the n1
  // principal axis to be very slightly longer than the other two.
  // here is make sure that the 2st eigenvector is the one at 12 o'clock
  if (BIAXIAL==-1) {
    TVector3 temp1=myz.Cross(rhat);
    TVector3 twelveoclock=rhat.Cross(temp1);
    twelveoclock.SetMag(1.);
    double mindotproduct=1000.;
    TVector3 threeoclock=-1.*temp1;
    threeoclock.SetMag(1.);
    mindotproduct=fabs(p_e2.Dot(twelveoclock));
    if (fabs(p_e1.Dot(twelveoclock)>mindotproduct)) {
      
      double n_e1_temp=n_e1;
      n_e1=n_e2;
      n_e2=n_e1_temp;
      
      TVector3 p_e1_temp=p_e1;
      p_e1=p_e2;
      p_e2=p_e1_temp;

    }
    p_e1=threeoclock;
    p_e2=twelveoclock;
      

  }

  double deltan=0.;
  if (BIAXIAL==-1)
    deltan=0.;
  else
    deltan=ne2new-ne1new;

  return deltan;

}
*/

double getDeltaN(double alpha) {

  //  const double NICE=1.78;
  //const double DELTAN=35.E-9/20.E-6;

  //  cout << "inside function, DELTAN is " << DELTAN << "\n";
  //cout << "sqrt is " << sqrt(1-sin(alpha)*sin(alpha)) << "\n";
  //  double n_L_prime=(NICE-DELTAN/2.)/sqrt(1-DELTAN/NICE*sin(alpha)*sin(alpha)*(1-DELTAN/(4*NICE)));
  double n_L_prime=(NICE-DELTAN/2.)/sqrt(1-DELTAN/NICE*cos(alpha)*cos(alpha)*(1-DELTAN/(4*NICE)));
  return 2*(NICE-n_L_prime);

}
/*
double getV(vector<double> nvec) {

  double alpha=nvec[0];
  double beta=nvec[1];
  double gamma=nvec[2];

  return 90.-acos((gamma*gamma)*(beta*beta-alpha*alpha)/((beta*beta)*(gamma*gamma-alpha*alpha)  ))*DEGRAD;

}
*/
/*
void thetastoEpsilons(double thetaE_e1_Sclock,double thetaE_e2_Sclock,
			double &epsilon1,double &epsilon2) {

  epsilon1=thetaE_e1_Sclock-PI/2.;
  epsilon2=thetaE_e2_Sclock;

}
*/
/*
void getAnglesontheClock(TVector3 rhat_thisstep1, TVector3 rhat_thisstep2,
			 TVector3 p_e1, TVector3 p_e2,
			 double &theta_e1,double &theta_e2,
			 double &p_e1_thetacomponent,double &p_e2_thetacomponent,
			 double &p_e1_phicomponent,double &p_e2_phicomponent) {

  double theta_I_1=rhat_thisstep1.Theta();
  double theta_I_2=rhat_thisstep2.Theta();

  TVector3 zaxis(0.,0.,1.);
  TVector3 yaxis(0.,1.,0.);
  // rhat_thisstep is the direction of khat.
  // rotate khat to be pointing in the +z direction, 
  // and rotate p_e1 and p_e2 along with them.
  double phi1=-1.*rhat_thisstep1.Phi();
  rhat_thisstep1.Rotate(phi1,zaxis);
  p_e1.Rotate(phi1,zaxis);

  double phi2=-1.*rhat_thisstep2.Phi();
  rhat_thisstep2.Rotate(phi2,zaxis);
  p_e2.Rotate(phi2,zaxis);

  double theta1=-1.*rhat_thisstep1.Theta();
  rhat_thisstep1.Rotate(theta1,yaxis);
  p_e1.Rotate(theta1,yaxis);

  double theta2=-1.*rhat_thisstep2.Theta();
  rhat_thisstep2.Rotate(theta2,yaxis);
  p_e2.Rotate(theta2,yaxis);

  // now find the angles on the clock that p_e1 and p_e2 sit at.
  theta_e1=p_e1.Phi();
  theta_e2=p_e2.Phi();

  double epsilon1=0.;
  double epsilon2=0.;
  thetastoEpsilons(theta_e1,theta_e2,
		   epsilon1,epsilon2);


  double scalefactor=0.2/(sin(theta_I_1)+sin(theta_I_2))*2.;

  p_e1_thetacomponent=scalefactor*sin(theta_I_1)*cos(epsilon1)*cos(epsilon1);
  p_e2_thetacomponent=scalefactor*sin(theta_I_2)*sin(epsilon2)*sin(epsilon2);

  p_e1_phicomponent=scalefactor*sin(theta_I_1)*sin(epsilon1)*cos(epsilon1);
  p_e2_phicomponent=scalefactor*sin(theta_I_2)*-1.*cos(epsilon2)*sin(epsilon2);

}
*/
void switchThem(double &thetaE_e1_Sclock,double &thetaE_e2_Sclock) {

  double tempangle=thetaE_e1_Sclock;
  thetaE_e1_Sclock=thetaE_e2_Sclock;
  thetaE_e2_Sclock=tempangle;
  
}
/*
void getManyAnglesontheClock(int BIAXIAL,double crosspolangle_tx,
			     TVector3 rhat_thisstep,
			     TVector3 p_e1,TVector3 p_e2,
			     TVector3 E_e1,TVector3 E_e2,
			     double &theta_e1,double &theta_e2,
			     double &thetaE_e1,double &thetaE_e2,
			     double &theta_e1_Sclock,double &theta_e2_Sclock,
			     double &thetaE_e1_Sclock,double &thetaE_e2_Sclock,
			     TVector3 &Shat_e1,TVector3 &Shat_e2,
			     double &E_e1_thetacomponent_Sclock,double &E_e2_thetacomponent_Sclock,
			     double &E_e1_phicomponent_Sclock,double &E_e2_phicomponent_Sclock) {

  // these aren't currently used for anything but they are the theta and phi components of D-hat for each ray, on the clock where k is in the page
  double p_e1_thetacomponent=0.;
  double p_e2_thetacomponent=0.;
  double p_e1_phicomponent=0.;
  double p_e2_phicomponent=0.;
  getAnglesontheClock(rhat_thisstep,rhat_thisstep, // we are using the same k vector for both rays
		      p_e1,p_e2,
		      theta_e1,theta_e2,
		      p_e1_thetacomponent,p_e2_thetacomponent,
		      p_e1_phicomponent,p_e2_phicomponent);


  if (theta_e1<0.)
    theta_e1+=PI;
  if (theta_e2<-1.*PI/2.)
    theta_e2+=PI;
  if (theta_e2>PI/2.)
    theta_e2-=PI;
	    

  theta_e1+=crosspolangle_tx;
  theta_e2+=crosspolangle_tx;

  // these aren't currently used for anything but they are the theta and phi components of the E field for each ray, on the clock where k is in the page
  double E_e1_thetacomponent=0.;
  double E_e2_thetacomponent=0.;
  double E_e1_phicomponent=0.;
  double E_e2_phicomponent=0.;

  getAnglesontheClock(rhat_thisstep,rhat_thisstep,
		      E_e1,E_e2,
		      thetaE_e1,thetaE_e2,
		      E_e1_thetacomponent,E_e2_thetacomponent,
		      E_e1_phicomponent,E_e2_phicomponent);
  
  if (thetaE_e1<0.)
    thetaE_e1+=PI;
  if (thetaE_e2<-1.*PI/2.)
    thetaE_e2+=PI;
  if (thetaE_e2>PI/2.)
    thetaE_e2-=PI;

  thetaE_e1+=crosspolangle_tx;
  thetaE_e2+=crosspolangle_tx;

  // 09/05/21 for uniaxial, problem is that p and E are parallel I think.
  TVector3 Hhat_e1;
  TVector3 Hhat_e2;
  
   // for an isotropic medium, I think Hhat_e1 keeps getting flipped
  // back and forth 
  // this part makes sure it stays oriented the say way relative to
  // each eigenvector and the direction of the ray.
  if (BIAXIAL==-1) {
    
//    TVector3 perpvec1=p_e1.Cross(rhat_thisstep);
  //  TVector3 perpvec2=p_e2.Cross(rhat_thisstep);
   // if (Hhat_e1.Dot(perpvec1)<0.)
    //  Hhat_e1=-1.*Hhat_e1;
    //if (Hhat_e2.Dot(perpvec2)<0.)
    //  Hhat_e2=-1.*Hhat_e2;
    
    //    Shat_e1=rhat_thisstep;
    //Shat_e2=rhat_thisstep;
    Hhat_e1=rhat_thisstep.Cross(p_e1);
    Hhat_e2=rhat_thisstep.Cross(p_e2);
    
    if (Hhat_e1.Mag()<HOWSMALLISTOOSMALL)
      cout << "Hhat_e1 is " << Hhat_e1.Mag() << "\n";

    Hhat_e1.SetMag(1);

    if (Hhat_e2.Mag()<HOWSMALLISTOOSMALL) {
      cout << "Shat_e2 is " << Shat_e2.Mag() << "\n";
      cout << "p_e2 is " << p_e2.Mag() << "\n";
      cout << "Hhat_e2 is " << Hhat_e2.Mag() << "\n";
    }
    Hhat_e2.SetMag(1);

    Shat_e1=E_e1.Cross(Hhat_e1);
    Shat_e2=E_e2.Cross(Hhat_e2);

    Shat_e1.SetMag(1.);
    Shat_e2.SetMag(1.);
  }
  else if (BIAXIAL==0) {

    Hhat_e1=rhat_thisstep.Cross(p_e1);
    Hhat_e2=rhat_thisstep.Cross(E_e2);

    Hhat_e1.SetMag(1.);
    Hhat_e2.SetMag(1.);

    Shat_e1=E_e1.Cross(Hhat_e1);
    Shat_e2=E_e2.Cross(Hhat_e2);

    Shat_e1.SetMag(1.);
    Shat_e2.SetMag(1.);


  }
  else {

    Hhat_e1=rhat_thisstep.Cross(p_e1);
    Hhat_e2=rhat_thisstep.Cross(p_e2);

    //    Hhat_e1=p_e1.Cross(E_e1);

    if (Hhat_e1.Mag()<HOWSMALLISTOOSMALL) {
      cout << "E_e1 is " << E_e1[0] << "\t" << E_e1[1] << "\t" << E_e1[2] << "\n";
      cout << "p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
      cout << "Hhat_e1 is " << Hhat_e1.Mag() << "\n";
    }

    //    Hhat_e2=p_e2.Cross(E_e2);
  
    Hhat_e1.SetMag(1.);
    Hhat_e2.SetMag(1.);

    Shat_e1=E_e1.Cross(Hhat_e1);
    Shat_e2=E_e2.Cross(Hhat_e2);

    Shat_e1.SetMag(1.);
    Shat_e2.SetMag(1.);

  }
  
  
  //cout << "before, E_e1 is " << E_e1[0] << "\t" << E_e1[1] << "\t" << E_e1[2] << "\n";
  //cout << "before, E_e2 is " << E_e2[0] << "\t" << E_e2[1] << "\t" << E_e2[2] << "\n";

  //cout << "before, Hhat_e1 is " << Hhat_e1[0] << "\t" << Hhat_e1[1] << "\t" << Hhat_e1[2] << "\n";
  //cout << "before, Hhat_e2 is " << Hhat_e2[0] << "\t" << Hhat_e2[1] << "\t" << Hhat_e2[2] << "\n";
  

  


  // here i want to plot where p_e1 and p_2 are on the clock, with 12 o'clock being the in the plane of rhat at launch and the z axis.
  // these aren't currently used for anything but they are the theta and phi components of the D-hat eigenvector for each ray, on the clock where S is in the page
  double p_e1_thetacomponent_Sclock=0.;
  double p_e2_thetacomponent_Sclock=0.;
  double p_e1_phicomponent_Sclock=0.;
  double p_e2_phicomponent_Sclock=0.;
  getAnglesontheClock(Shat_e1,Shat_e2,
		      p_e1,p_e2,
		      theta_e1_Sclock,theta_e2_Sclock,
		      p_e1_thetacomponent_Sclock,p_e2_thetacomponent_Sclock,
		      p_e1_phicomponent_Sclock,p_e2_phicomponent_Sclock);
  
  
  //cout << "before, p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
  //cout << "before, p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";
  //cout << "before, Shat_e1 is " << Shat_e1[0] << "\t" << Shat_e1[1] << "\t" << Shat_e1[2] << "\n";
  //cout << "before, Shat_e2 is " << Shat_e2[0] << "\t" << Shat_e2[1] << "\t" << Shat_e2[2] << "\n";

//  cout << "before, theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";
  
  
  if (theta_e1_Sclock<0.)
    theta_e1_Sclock+=PI;
  if (theta_e2_Sclock<-1.*PI/2.)
    theta_e2_Sclock+=PI;
  if (theta_e2_Sclock>PI/2.)
    theta_e2_Sclock-=PI;

  //  cout << "after, theta_e1_Sclock, theta_e2_Sclock are " << theta_e1_Sclock << "\t" << theta_e2_Sclock << "\n";

  theta_e1_Sclock+=crosspolangle_tx;
  theta_e2_Sclock+=crosspolangle_tx;

  getAnglesontheClock(Shat_e1,Shat_e2,
		      E_e1,E_e2,
		      thetaE_e1_Sclock,thetaE_e2_Sclock,
		      E_e1_thetacomponent_Sclock,E_e2_thetacomponent_Sclock,
		      E_e1_phicomponent_Sclock,E_e2_phicomponent_Sclock);

  if (thetaE_e1_Sclock<0.)
    thetaE_e1_Sclock+=PI;
  if (thetaE_e2_Sclock<-1.*PI/2.)
    thetaE_e2_Sclock+=PI;
  if (thetaE_e2_Sclock>PI/2.)
    thetaE_e2_Sclock-=PI;
  
  thetaE_e1_Sclock+=crosspolangle_tx;
  thetaE_e2_Sclock+=crosspolangle_tx;


  //  E_e1_thetacomponent_Sclock = E_e1.Mag()*cos(thetaE_e1_Sclock); // theta component of E_e1
  //E_e2_thetacomponent_Sclock = -1.*E_e2.Mag()*sin(thetaE_e1_Sclock); // theta component of E_e2

  //E_e1_phicomponent_Sclock = E_e1.Mag()*sin(thetaE_e1_Sclock); // phi component of E_e1
  //E_e2_phicomponent_Sclock = E_e2.Mag()*cos(thetaE_e1_Sclock); // phi component of E_e22

}
*/

double Flipped(double theta_e1,double theta_e1_start) {


  double diff=theta_e1-theta_e1_start;
  if (diff>PI)
    diff=diff-2.*PI;
  if (diff<-1.*PI)
    diff=diff+2.*PI;



  double notflipped=1.;

  if (abs(diff)>PI/2.*0.9)
    notflipped=-1.;

  //  cout << "theta_e1, theta_e1_start, diff are " << theta_e1 << "\t" << theta_e1_start << "\t" << diff << "\t" << notflipped << "\n";

  return notflipped;

}
// double Flipped(TVector3 p_e1_previous, TVector3 p_e2_previous, TVector3 &p_e1, TVector3 &p_e2) {
  
//   double notflipped=1.;
  
//   double angle1=acos(p_e1_previous.Dot(p_e1)/p_e1.Mag()/p_e1_previous.Mag());
//   double angle2b=acos(p_e2_previous.Dot(p_e2)/p_e2.Mag()/p_e2_previous.Mag());
  
//   if (fabs(angle1)*DEGRAD>170.)
//     p_e1=-1.*p_e1;

//   if (fabs(angle2b)*DEGRAD>170.)
//     p_e2=-1.*p_e2;

//   angle1=acos(p_e1_previous.Dot(p_e1)/p_e1.Mag()/p_e1_previous.Mag());
//   double angle2=acos(p_e2_previous.Dot(p_e1)/p_e1.Mag()/p_e2_previous.Mag());

//   //  if (angle2<angle1 && fabs(angle1)<170./DEGRAD)
//   if (fabs(angle2)<fabs(angle1)) {
//     notflipped=-1.;
//     //TVector3 temp=p_e1;
//     //p_e1=p_e2;
//     //p_e2=temp;
//   }

//   // angle1=acos(p_e1_previous.Dot(p_e1)/p_e1.Mag()/p_e1_previous.Mag());
//   // angle2b=acos(p_e2_previous.Dot(p_e2)/p_e2.Mag()/p_e2_previous.Mag());
  
//   // if (fabs(angle1)*DEGRAD>175.)
//   //   p_e1=-1.*p_e1;

//   // if (fabs(angle2b)*DEGRAD>175.)
//   //   p_e2=-1.*p_e2;
//   return notflipped;


// }
TCanvas * makePretty2Panel(){
   TCanvas *ccc = new TCanvas("plotEvent","plotEvent", 800, 400);
   ccc->Divide(2, 0);
   //   ccc->GetPad(1)->SetPad(.005, .005, .3475, .995);
   ccc->GetPad(1)->SetPad(.005, .005, .4975, .995);
   ccc->GetPad(1)->Divide(0, 2);
   //   ccc->GetPad(2)->SetPad(.3525, .005, .995, .995);
   ccc->GetPad(2)->SetPad(.5025, .005, .995, .995);
   ccc->SetWindowSize(1200, 700);
   ccc->cd(1)->cd(1)->SetGrid();
//stuff for first plot here
//these are good margins for something like a graph with no colorbar
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.12);
  ccc->cd(1)->cd(2);
//stuff for second plot here  
  //these are good margins for something with a colorbar
  gPad->SetBottomMargin(.12);
  gPad->SetRightMargin(.19);
  gPad->SetLeftMargin(.15);
  return ccc;
}
void titles(TGraph *inGr, TString title, TString xtitle, TString ytitle){
  auto sizeT=.055;
  inGr->SetTitle(title);
  inGr->GetXaxis()->SetTitle(xtitle);
  inGr->GetYaxis()->SetTitle(ytitle);
  inGr->GetXaxis()->SetTitleSize(sizeT);
  inGr->GetYaxis()->SetTitleSize(sizeT);
  inGr->GetXaxis()->SetLabelSize(sizeT);
  inGr->GetYaxis()->SetLabelSize(sizeT);
  //inGr->GetXaxis()->SetTitleOffset(1.2);                                                                                                            
  inGr->GetYaxis()->SetLabelOffset(.01);
  inGr->GetYaxis()->SetTitleOffset(1.2);
}
void titles(TH2 *inH, TString title, TString xtitle, TString ytitle){
  auto sizeT=.055;
  inH->SetTitle(title);
  inH->GetXaxis()->SetTitle(xtitle);
  inH->GetYaxis()->SetTitle(ytitle);
  inH->GetXaxis()->SetTitleSize(sizeT);
  inH->GetYaxis()->SetTitleSize(sizeT);
  inH->GetXaxis()->SetLabelSize(sizeT);
  inH->GetYaxis()->SetLabelSize(sizeT);
  //inGr->GetXaxis()->SetTitleOffset(1.2);                                                                                                            
  inH->GetYaxis()->SetLabelOffset(.01);
  inH->GetYaxis()->SetTitleOffset(1.2);
}
