// psi_model.hh
//
// Core per-station birefringence ("psi model") computation extracted from
// main() in cpol_main.cc. For one station, process_station() ray-traces from
// every pulser depth to the station, transports the two birefringent
// eigenmodes along the path, and fills the per-depth observables in StationData
// (epsilons, beam factors, voltages, envelopes, polarization angles, ...).
#pragma once

#include "config.hh"
#include "geometry.hh"
#include "ice_profile.hh"
#include "station_data.hh"

#include "TVector3.h"

// Lightweight accumulator for the wall-clock timing diagnostics that main()
// prints after the station loop. Carried across stations by reference.
struct PsiTiming {
    double uzair_time = 0.0;
    double overhead   = 0.0;
    double loop_time  = 0.0;
    double eval_time  = 0.0;
    double ray_time   = 0.0;
};

// Run the full per-depth model for a single station i.
//
// p_e1, p_e2, n_e1, n_e2 carry the current eigenmode state and are updated in
// place (they persist across stations, exactly as in the original main loop).
// SKIP_MIDDLE_STEPS / SKIP_RAYTRACE_PATH select the fast code paths.
void process_station(const Config& cfg, Geometry& geom, Ice_Profile& ice,
                     StationData& model_data, int i,
                     TVector3& p_e1, TVector3& p_e2, double& n_e1, double& n_e2,
                     bool SKIP_MIDDLE_STEPS, bool SKIP_RAYTRACE_PATH,
                     PsiTiming& timing);
