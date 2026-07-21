// config.hh
#pragma once

#include <vector>

struct Config {
    // Data members with defaults
    double freq = 160.E6;
    int CROSSPOLANGLE_TX_INT = 0.0;
    int CROSSPOLANGLE_RX_INT = 0.0;
    int BIAXIAL = 1;
    int CONSTANTINDICATRIX = 0;

    // Derived quantities
    double CROSSPOLANGLE_TX = 0.0;
    double CROSSPOLANGLE_RX = 0.0;

    double freqmin = 0.;
    double freqmax = 1.E9;
    static const int NFREQ = 100;
    std::vector<double> vfreqs;

    // Initialize rotation quantities
    double phi = 0.0;
    double theta = 0.0;
    double gamma = 0.0;

    // Constructor takes argc/argv and does all the parsing
    Config(int argc, char** argv);
};

