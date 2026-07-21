// cpol_helpers.hh
//
// Small, self-contained free functions that were previously defined at the
// bottom of cpol_main.cc. These are generic geometry / plotting-style helpers
// that do not depend on any of the program's run-time state objects
// (Config, Geometry, StationData, ...), so they live in their own module.
#pragma once

#include <vector>

#include "TVector3.h"
#include "TGraph.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"

// Build a normalized direction by weighting p_e1's Cartesian components by nvec.
// (p_e2 is currently unused but kept for signature compatibility.)
TVector3 directionNextStep(std::vector<double> nvec, TVector3 p_e1, TVector3 p_e2);

// Rotate E into the ice-flow-aligned frame, apply the diagonal dielectric
// tensor epsilon, then rotate back to the lab frame.
TVector3 rotateE(TVector3 epsilon, double angle_iceflow, TVector3 E);

// From principal indices, ice-flow angle, mode index n and a D vector, compute
// the propagation direction k (returned) and the electric field E (by ref).
TVector3 getNewkandE(std::vector<double> nvec, double angle_iceflow, double n,
                     TVector3 D, TVector3 kguess, TVector3 &E);

// Analytic toy-model delta n as a function of angle alpha to the ice-flow axis.
double getDeltaN(double alpha);

// Swap two angles in place (used for ARIANNA / biaxial mode relabeling).
void switchThem(double &thetaE_e1_Sclock, double &thetaE_e2_Sclock);

// Detect an eigenvector branch flip (~pi/2 change); returns -1 if flipped, else +1.
double Flipped(double theta_e1, double theta_e1_start);

// Create a styled two-panel canvas.
TCanvas *makePretty2Panel();

// Apply consistent axis label / title styling.
void titles(TGraph *inGr, TString title, TString xtitle, TString ytitle);
void titles(TH2 *inH, TString title, TString xtitle, TString ytitle);
