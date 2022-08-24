/////////////////////////////////////////////////////////////////                                                                                                  
// C++ example program to demonstrate usage of birefringence code. //                                                                                              
// Comments to Amy Connolly <connolly(at)physics.osu.edu>.               //                                                                                               
/////////////////////////////////////////////////////////////////////                                                                                              


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
#include <vector>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TVector2.h"
//#include "units.h"                                                                                                                                               
//#include "raytracing_tools.h"                                                                                                                                    
//#include "IceRayTracing-master/IceRayTracing.h"
#include "TLegend.h"


using namespace std;
#include "birefringence.hh"

const double PI=3.1415926;
const double CLIGHT=3.E8;
const double DEGRAD=180./PI;
const double HOWSMALLISTOOSMALL=1.e-8;

int main(int argc, char** argv) {

  // so far this just calculates epsilon angles for a k vector and an indicatrix that you give it, but more to come.


  vector<double> nvec_thisstep;
  nvec_thisstep.resize(3);

  // this is just for Fig. 5
  // set n_alpha, n_beta, and n_gamma here.
  nvec_thisstep[0]=1.0;
  nvec_thisstep[1]=1.5;
  nvec_thisstep[2]=1.3;

  // this is the epsilon tensor.
  TVector3 epsilon_thisstep;
  epsilon_thisstep[0]=nvec_thisstep[0]*nvec_thisstep[0];
  epsilon_thisstep[1]=nvec_thisstep[1]*nvec_thisstep[1];
  epsilon_thisstep[2]=nvec_thisstep[2]*nvec_thisstep[2];

  // 1=biaxial, 0=uniaxial, -1=isotropic
  int BIAXIAL=1;

  // here let's just be in a frame where the iceflow is along the x axis
  double angle_iceflow=0.;

  // the azimuth and elevation angles of the incident k vector
  double az=25./DEGRAD;
  double el=25./DEGRAD;

  // the k vector
  TVector3 rhat_thisstep(sin(-1.*az)*cos(el),cos(-1.*az)*cos(el),sin(el));

  // these will be the indices of refraction seen by the two rays
  double n_e1,n_e2;
  // the eigenvectors of D for the two rays.
  TVector3 p_e1;
  TVector3 p_e2;

  double deltan_alongpath=getDeltaN(BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,
				    n_e1,n_e2,p_e1,p_e2); // outputs

  cout << "n_e1, n_e2 are " << n_e1 << "\t" << n_e2 << "\n";
  cout << "delta n is " << deltan_alongpath << "\n";
  cout << "p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
  cout << "p_e2 is " << p_e2[0] << "\t" << p_e2[1] << "\t" << p_e2[2] << "\n";

  // from D's, find E's.
  TVector3 E_e1=rotateD(epsilon_thisstep,angle_iceflow,p_e1);
  TVector3 E_e2=rotateD(epsilon_thisstep,angle_iceflow,p_e2);

  // choice of cross-pol angle delta
  double crosspolangle_tx=0.;
  
  // these will all be outputs to getManyAnglesontheClock
  double theta_e1,theta_e2;
  double thetaE_e1,thetaE_e2;
  double theta_e1_Sclock,theta_e2_Sclock;
  double thetaE_e1_Sclock,thetaE_e2_Sclock;
  TVector3 Shat_e1,Shat_e2;
  double E_e1_thetacomponent,E_e2_thetacomponent;
  double E_e1_phicomponent,E_e2_phicomponent;
  
  getManyAnglesontheClock(BIAXIAL,crosspolangle_tx, // inputs
			  rhat_thisstep,
			  p_e1,p_e2,E_e1,E_e2,
			  theta_e1,theta_e2,thetaE_e1,thetaE_e2, // outputs
			  theta_e1_Sclock,theta_e2_Sclock,
			  thetaE_e1_Sclock,thetaE_e2_Sclock,
			  Shat_e1,Shat_e2,
			  E_e1_thetacomponent,E_e2_thetacomponent,
			  E_e1_phicomponent, E_e2_phicomponent);
  

  cout << "thetas are " << theta_e1*DEGRAD << "\t" << theta_e2*DEGRAD << "\n";
  
  // now find the epsilon angles (relative to 12 o'clock and 3 o'clock)
  double epsilon1,epsilon2;
  thetastoEpsilons(thetaE_e1_Sclock,thetaE_e2_Sclock,                                                            
		   epsilon1,epsilon2);   
  
  cout << "epsilons are " << epsilon1*DEGRAD << "\t" << epsilon2*DEGRAD << "\n";

}
