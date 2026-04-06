#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "TVector3.h"
#include "constants.hh"

// ------------------------------------------------------------------------
// Geometry setup
// ------------------------------------------------------------------------
struct Geometry {
  const double DEPTH=-600.;
  const double depth_special=-1100.;

  // Depth grid used for "special" evaluations (e.g., diagnostic depths)
  const double stepspecial=20.;
  static constexpr int NSPECIAL=50;
  double whichspecial[NSPECIAL];
  double startspecial=-1600.;

  const int WHICHPOL=0;
  const int DORAYTRACING=1;

  // ------------------------------------------------------------------------
  // Station definitions: A1..A5 plus ARIANNA
  // ------------------------------------------------------------------------
  static constexpr int NSTATIONS=6;
  const int minstation=0;
  const int maxstation=5;


  int igreatestdepth[NSTATIONS];
  int imostshallowdepth[NSTATIONS];

  // Horizontal pulser->station distances (in meters after MFT conversion)
  double horizontal_distances[NSTATIONS]={1257.,2353.,3146.,3199.,5179.8892,653.804525};
  std::string snames[NSTATIONS]={"A1","A2","A3","A4","A5","ARIANNA"};

  // ARIANNA azimuth/phi angle (radians)
  double PHI_ARIANNA=(90.+36.+46./60.+23./3600.+1.4)/DEGRAD;


  const int NSTEPS=100;
  double min_altitude=-1000.;
  double max_altitude=-600.;

  // ------------------------------------------------------------------------
  // Station coordinates and pulser coordinates (originally in feet; converted to meters)
  // ------------------------------------------------------------------------
  double pulser_coords[2]={42358.94,48974.2};

  double station_coords[NSTATIONS][2]={
      {38754., 51051.}, // A1
      {35481., 45369.}, // A2
      {32200., 51053.}, // A3
      {35478., 56737.}, // A4
      {32356., 39746.}, // A5
      {41153.,50381.75} // ARIANNA
  };
  // ------------------------------------------------------------------------
  // Ice-flow direction and "ordinary" axis direction in horizontal plane.
  // angle_iceflow: azimuth defining principal axis orientation relative to x-y.
  // ------------------------------------------------------------------------
  const double angle_iceflow=(36.+ (46./60.) + (23./3600.) + 90.)/DEGRAD;
  TVector3 ordinary = TVector3 (std::cos(angle_iceflow),std::sin(angle_iceflow),0.);

  // Default TX/RX orientations (both along +z)
  TVector3 tx_orientation = TVector3 (0.,0.,1.);
  TVector3 rx1_orientation = TVector3 (0.,0.,1.);

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

  Geometry();

  void compute_station_geometry(int BIAXIAL, const std::vector<double>& vec);
};
