// Constants needed for the psi model
#pragma once

// Constants/conversions
const double N_ICE=1.78;
const double NICE=1.78;
const double DELTA_N=0.427;
const double VOLTAGENORM=150.0 * 11.0 / 7.0;
const double PI=3.1415926;
const double CLIGHT=3.E8;
const double DEGRAD=180./PI;
const double MFT=1./100.*2.54/1.*12.;
const double L_ATTEN=2200.;

// Small numeric tolerances / discretization knobs
const double DELTAN=0.0005;
const double HOWSMALLISTOOSMALL=1.e-8;
const int UZAIRSTEP=5;
double freq=160.E6;             // default frequency [Hz]

// Birefringence parameters
int CROSSPOLANGLE_TX_INT=0;
int CROSSPOLANGLE_RX_INT=0;
int BIAXIAL=1;
int CONSTANTINDICATRIX=0;

// Frequencies to step through 
double freqmin=0.;
double freqmax=1.E9;
const int NFREQ=100;

