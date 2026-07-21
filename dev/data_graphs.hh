#pragma once

#include <vector>
#include "TVector3.h"

struct Graphs {
    // --------------------------------------------------------------------
    // Containers for reading data and for storing outputs.
    // Each outer vector index corresponds to a station.
    // --------------------------------------------------------------------
    std::vector<std::vector<double>> videpth;
    std::vector<std::vector<double>> vdepth;

    std::vector<std::vector<double>> vdepth_data;
    std::vector<std::vector<double>> vdepth_data_err;

    std::vector<std::vector<double>> vreversedepth;

    std::vector<std::vector<double>> vreceiveangle;
    std::vector<std::vector<double>> vlaunchangle;
    std::vector<std::vector<double>> voutput6;
    std::vector<std::vector<double>> voutput7;
    std::vector<std::vector<double>> voutput8;

    std::vector<std::vector<double>> vtotal_distances;
    std::vector<std::vector<double>> vtotal_distances_err;

    std::vector<std::vector<double>> vsnrmax;
    std::vector<std::vector<double>> vsnrmax_err;

    // Various angle diagnostics between propagation direction / eigenvectors
    std::vector<std::vector<double>> vangle_khat_0_khat_1_2;
    std::vector<std::vector<double>> vangle_khat_0_khat_2_2;
    std::vector<std::vector<double>> vangle_khat_1_2_khat_2_2;
    std::vector<std::vector<double>> vangle_Shat_e1_khat;
    std::vector<std::vector<double>> vangle_Shat_e2_khat;

    // Ray-tracing path containers
    TVector3 raypath[2][6];
    TVector3 raypath_n[2][6];
    std::vector<std::vector<double>> vistep;
    std::vector<std::vector<std::vector<double>>> vraypos;
    std::vector<std::vector<double>> vraypath_ne1;
    std::vector<std::vector<double>> vraypath_ne2;

    // Beam / attenuation vs depth containers
    std::vector<std::vector<double>> vrxdepth_beam1;
    std::vector<std::vector<double>> vrxdepth_beam2;
    std::vector<std::vector<double>> vtxdepth_beam1;
    std::vector<std::vector<double>> vtxdepth_beam2;

    std::vector<std::vector<double>> vrxdepth_atten;
    std::vector<std::vector<double>> vrxdepth_notflipped;
    std::vector<std::vector<double>> vrxdepth_atten_beam;
    std::vector<std::vector<double>> vrxdepth_atten_power;
    std::vector<std::vector<double>> vrxdepth_atten_beam_power;

    // Angular evolution of E-field components (theta1/theta2) along the path
    std::vector<std::vector<double>> vtxdepth_theta1;
    std::vector<std::vector<double>> vtxdepth_theta2;
    std::vector<std::vector<double>> vrxdepth_theta1;
    std::vector<std::vector<double>> vrxdepth_theta2;

    // Same angles, but expressed in a "S-clock" coordinate convention
    std::vector<std::vector<double>> vtxdepth_theta1_Sclock;
    std::vector<std::vector<double>> vtxdepth_theta2_Sclock;
    std::vector<std::vector<double>> vrxdepth_theta1_Sclock;
    std::vector<std::vector<double>> vrxdepth_theta2_Sclock;

    // E-field polarization angles (and S-clock versions)
    std::vector<std::vector<double>> vtxdepthE_theta1;
    std::vector<std::vector<double>> vtxdepthE_theta2;
    std::vector<std::vector<double>> vrxdepthE_theta1;
    std::vector<std::vector<double>> vrxdepthE_theta2;
    std::vector<std::vector<double>> vtxdepthE_theta1_Sclock;
    std::vector<std::vector<double>> vtxdepthE_theta2_Sclock;
    std::vector<std::vector<double>> vrxdepthE_theta1_Sclock;
    std::vector<std::vector<double>> vrxdepthE_theta2_Sclock;

    // Dispersion contributions along the path
    std::vector<std::vector<double>> vtxdepth_dispersion1;
    std::vector<std::vector<double>> vtxdepth_dispersion2;

    // Dot products of S/E/D unit vectors with something (TX frame)
    std::vector<std::vector<double>> vdotShats_tx;
    std::vector<std::vector<double>> vdotEhats_tx;
    std::vector<std::vector<double>> vdotDhats_tx;

    // --------------------------------------------------------------------
    // Per-station waveform / voltage / power bookkeeping
    // --------------------------------------------------------------------
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

    void resize_station_vectors(std::size_t nstations = 6) {
        vdotShats_tx.resize(nstations);
        vdotEhats_tx.resize(nstations);
        vdotDhats_tx.resize(nstations);

        vepsilon1_tx.resize(nstations);
        vepsilon2_tx.resize(nstations);
        vdiffepsilon_tx.resize(nstations);

        vepsilon1_rx.resize(nstations);
        vepsilon2_rx.resize(nstations);
        vdiffepsilon_rx.resize(nstations);

        vV1_r1.resize(nstations);
        vV2_r1.resize(nstations);
        vE1_r1.resize(nstations);
        vE2_r1.resize(nstations);
        vV1_r1_lpda.resize(nstations);
        vV2_r1_lpda.resize(nstations);
        vV1_r2.resize(nstations);
        vV2_r2.resize(nstations);
        vE1_r2.resize(nstations);
        vE2_r2.resize(nstations);

        vV1V2_r1.resize(nstations);
        vV1V2_r1_lpda.resize(nstations);
        vV1V2_r2.resize(nstations);
        voppositeV1V2_r2.resize(nstations);
        voppositeV1V2_r1.resize(nstations);

        vV1squared_r1.resize(nstations);
        vV2squared_r1.resize(nstations);
        vV1squared_r2.resize(nstations);
        vV2squared_r2.resize(nstations);

        vpower_r1.resize(nstations);
        vpoynting_r1.resize(nstations);
        vpower_r1_lpda.resize(nstations);
        vpower_r2.resize(nstations);
        vpoynting_r2.resize(nstations);
        vvoltage_r1.resize(nstations);
        vfield_r1.resize(nstations);
        vvoltage_r1_lpda.resize(nstations);
        vvoltage_r2.resize(nstations);
        vfield_r2.resize(nstations);

        vpolarization_Psi_rx.resize(nstations);
        vpolarization_Omega_rx.resize(nstations);
        vEpolarization_Psi_rx.resize(nstations);
        vEpolarization_Omega_rx.resize(nstations);
        vEpolarization_reversedepth_Psi_rx.resize(nstations);
        vEpolarization_reversedepth_Omega_rx.resize(nstations);
        vpolarization_reversedepth_Psi_rx.resize(nstations);
        vpolarization_reversedepth_Omega_rx.resize(nstations);

        venvelope_minus_r1.resize(nstations);
        vSenvelope_minus_r1.resize(nstations);
        venvelope_minus_r1_lpda.resize(nstations);
        venvelope_plus_r1.resize(nstations);
        venvelope_plus_r1_lpda.resize(nstations);
        venvelope_minus_r2.resize(nstations);
        vSenvelope_minus_r2.resize(nstations);
        venvelope_plus_r2.resize(nstations);
        vSenvelope_plus_r2.resize(nstations);
        vvenvelope_minus_r1.resize(nstations);
        vvenvelope_plus_r1.resize(nstations);
        vEenvelope_minus_r1.resize(nstations);
        vEenvelope_plus_r1.resize(nstations);
        vSenvelope_plus_r1.resize(nstations);
        vvenvelope_minus_r2.resize(nstations);
        vvenvelope_plus_r2.resize(nstations);
        vEenvelope_minus_r2.resize(nstations);
        vEenvelope_plus_r2.resize(nstations);

        vistep.resize(nstations);
        vraypos.resize(nstations);
        vraypath_ne1.resize(nstations);
        vraypath_ne2.resize(nstations);
        vrxdepth_beam1.resize(nstations);
        vrxdepth_beam2.resize(nstations);
        vtxdepth_beam1.resize(nstations);
        vtxdepth_beam2.resize(nstations);

        vrxdepth_atten.resize(nstations);
        vrxdepth_notflipped.resize(nstations);
        vrxdepth_atten_beam.resize(nstations);
        vrxdepth_atten_power.resize(nstations);
        vrxdepth_atten_beam_power.resize(nstations);
        vtxdepth_theta1.resize(nstations);
        vtxdepth_theta2.resize(nstations);
        vrxdepth_theta1.resize(nstations);
        vrxdepthE_theta2.resize(nstations);
        vrxdepthE_theta1.resize(nstations);
        vrxdepth_theta2.resize(nstations);
        vtxdepth_theta1_Sclock.resize(nstations);
        vtxdepth_theta2_Sclock.resize(nstations);
        vrxdepth_theta1_Sclock.resize(nstations);
        vrxdepth_theta2_Sclock.resize(nstations);
        vtxdepth_dispersion1.resize(nstations);
        vtxdepth_dispersion2.resize(nstations);
        vtxdepthE_theta1.resize(nstations);
        vtxdepthE_theta2.resize(nstations);
        vtxdepthE_theta1_Sclock.resize(nstations);
        vtxdepthE_theta2_Sclock.resize(nstations);
        vrxdepthE_theta1_Sclock.resize(nstations);
        vrxdepthE_theta2_Sclock.resize(nstations);

        // Basic scan/data containers
        videpth.resize(nstations);
        vdepth.resize(nstations);
        vdepth_data.resize(nstations);
        vdepth_data_err.resize(nstations);
        vreversedepth.resize(nstations);

        vreceiveangle.resize(nstations);
        vlaunchangle.resize(nstations);
        voutput6.resize(nstations);
        voutput7.resize(nstations);
        voutput8.resize(nstations);

        vtotal_distances.resize(nstations);
        vtotal_distances_err.resize(nstations);

        vsnrmax.resize(nstations);
        vsnrmax_err.resize(nstations);

        vangle_khat_0_khat_1_2.resize(nstations);
        vangle_khat_0_khat_2_2.resize(nstations);
        vangle_khat_1_2_khat_2_2.resize(nstations);
        vangle_Shat_e1_khat.resize(nstations);
        vangle_Shat_e2_khat.resize(nstations);

        pol_r1.resize(nstations);
        pol_r2.resize(nstations);

        for (std::size_t istation = 0; istation < nstations; istation++) {
            vraypos[istation].resize(3);
        }
    }
};
