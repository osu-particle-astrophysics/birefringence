#include "config.hh"
#include "constants.hh"
#include <getopt.h>
#include <iostream>

using namespace std;

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


Config::Config(int argc, char** argv) {
  // Parse command line
  
  char clswitch;
  if (argc > 1) {
      while ((clswitch = getopt(argc, argv, "b:f:t:u:c:")) != EOF) {
          switch(clswitch) {
          case 'b':
              BIAXIAL = atoi(optarg);
              break;
          case 'f':
              freq = (double)atof(optarg) * 1.E6;
              break;
          case 't':
              CROSSPOLANGLE_TX_INT = atoi(optarg);
              break;
          case 'u':
              CROSSPOLANGLE_RX_INT = atoi(optarg);
              break;
          case 'c':
              CONSTANTINDICATRIX = atoi(optarg);
              break;
          }
      }
  }
  // convert crosspol angles to radians
  CROSSPOLANGLE_TX = (double)CROSSPOLANGLE_TX_INT / DEGRAD;
  CROSSPOLANGLE_RX = (double)CROSSPOLANGLE_RX_INT / DEGRAD;

  vfreqs.reserve(NFREQ);
  for (int i = 0; i < NFREQ; i++) {
      vfreqs.push_back(freqmin + (freqmax - freqmin) / (double)NFREQ * (double)i);
  }
}
