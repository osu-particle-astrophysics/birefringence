#pragma once

#include <vector>
#include <string>
#include <fstream>
#include "config.hh"
#include "birefringence.hh"

struct Ice_Profile{

  static constexpr int NDEPTHS_NS=81;

  // Depth-dependent principal indices (text files)
  std::string sn1file="data/n1.txt";
  std::string sn2file="data/n2.txt";
  std::string sn3file="data/n3.txt";

  // Depth arrays and refractive index arrays for each principal axis
  std::vector<double> vdepths_n1, vdepths_n2, vdepths_n3;
  std::vector<double> n1vec, n2vec, n3vec;

  // Derived quantity vs depth (e.g., V(n1,n2,n3) from birefringence.hh)
  std::vector<double>  vV;

  Ice_Profile(const Config& cfg, const std::vector<double>& nvec);

  void smooth_indices(std::vector<double>& n_vec); 
  // Declare some TGraphs for storing data
  TGraph* gn1 = nullptr;
  TGraph* gn2 = nullptr;
  TGraph* gn3 = nullptr;
  TGraph* g_V = nullptr;

};

