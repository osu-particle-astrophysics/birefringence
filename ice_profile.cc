#include <iostream>
#include <cmath>
#include "constants.hh"
#include "ice_profile.hh"
#include "config.hh"

using namespace std;

Ice_Profile::Ice_Profile(const Config& cfg, const vector<double>& nvec) {
  // Read in the files
  ifstream n1file(sn1file.c_str());
  ifstream n2file(sn2file.c_str());
  ifstream n3file(sn3file.c_str());

  std::string stemp1, stemp2, stemp3;
  double thisdepth1, thisdepth2, thisdepth3;
  double thisn1, thisn2, thisn3;


  // ---- n1(z) ----
  n1file >> stemp1; // discard header token
  n2file >> stemp2; // discard header token
  n3file >> stemp3; // discard header token
  // read in all the index data in one loop

  for (int i=0;i<NDEPTHS_NS;i++) {
    // read the data at this line
    n1file >> thisdepth1 >> thisn1;
    n2file >> thisdepth2 >> thisn2;
    n3file >> thisdepth3 >> thisn3;
    // set the depth for each direction
    vdepths_n1.push_back(-1.*thisdepth1);   // flip sign
    vdepths_n2.push_back(-1.*thisdepth2);   // flip sign
    vdepths_n3.push_back(-1.*thisdepth3);   // flip sign
    // set n for each direction
    n1vec.push_back(thisn1);
    if (cfg.BIAXIAL==1) { // Everything different
        n2vec.push_back(thisn2);
        n3vec.push_back(thisn3);
    }
    else if (cfg.BIAXIAL==0) { // Two identical
      n2vec.push_back(thisn1);
      n3vec.push_back(thisn3);
    }
    else if (cfg.BIAXIAL==-1){ // All identical
      n2vec.push_back(thisn1);
      n3vec.push_back(thisn1 + 1.E-5); // Need slight difference to break degeneracy
    }
  }
 
  if (cfg.CONSTANTINDICATRIX==1) {
    n1vec.assign(NDEPTHS_NS, nvec[0]);
    n2vec.assign(NDEPTHS_NS, nvec[1]);
    n3vec.assign(NDEPTHS_NS, nvec[2]);
  }

  // Debug printout of sizes and raw indices
  cout << "sizes are " << n1vec.size() << "\t" << n2vec.size() << "\t" << n3vec.size() << "\n";
  cout << "n's are \n";
  for (int i=0;i<NDEPTHS_NS;i++) {
      cout << "n1, n2, n3 are " << n1vec[i] << "\t" << n2vec[i] << "\t" << n3vec[i] << "\n";
  }

  // Smooth the indices
  smooth_indices(n1vec);
  smooth_indices(n2vec);
  smooth_indices(n3vec);
}

void Ice_Profile::smooth_indices(std::vector<double>& n_vec) {
  // ------------------------------------------------------------------------
  // Smooth n1/n2/n3 vs depth using a simple moving average (boxcar) filter.
  // NSMOOTH=5 -> replace each interior point with average over i-2..i+2.
  // Edges are left unchanged for the first/last 2 points.
  // ------------------------------------------------------------------------
  vector<double> tmp;
  tmp.resize(n_vec.size());
  int NSMOOTH=5;
  int min=(int)(((double)NSMOOTH)/2.);    // = 2 for NSMOOTH=5

  // Copy edges unchanged
  for (int i=0;i<min;i++) {
      tmp[i]=n_vec[i];
  }
  for (int i=n_vec.size()-(NSMOOTH-min);i<n_vec.size();i++) {
      tmp[i]=n_vec[i];
  }

  // Smooth interior points
  for (int i=min;i<n_vec.size()-(NSMOOTH-min);i++) {
      double tmpdouble=0.;
      for (int j=i-min;j<i+(NSMOOTH-min);j++) {
          tmpdouble+=n_vec[j];
      }
      tmpdouble=tmpdouble/(double)NSMOOTH;
      tmp[i]=tmpdouble;
  }
  n_vec=tmp;

}
