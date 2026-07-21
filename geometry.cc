#include <iostream>
#include <cmath>
#include "constants.hh"
#include "geometry.hh"

using namespace std;

Geometry::Geometry() {
  for (int ispecial=0;ispecial<NSPECIAL;ispecial++) {
      whichspecial[ispecial]=startspecial+stepspecial*(double)ispecial;
  }
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
}

void Geometry::compute_station_geometry(int BIAXIAL, const std::vector<double>& nvec){
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

  if (thisalpha>PI) thisalpha-=2.*PI;
  if (thisalpha<-1.*PI) thisalpha+=2.*PI;
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
  deltan_exp[5]=getDeltaN(BIAXIAL,nvec,pulsertostationhat_specialdepth[5],angle_iceflow,n_e1,n_e2,p_e1,p_e2);
  alpha_deg[5]=alpha[5]*DEGRAD;

}
