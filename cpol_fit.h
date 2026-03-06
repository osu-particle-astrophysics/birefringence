#ifndef CPOL_H
#define CPOL_H

#include <iostream>

int psiModel(int argc, char** argv, const Double_t* par_fit, Double_t* &fit_depths, Double_t* &psi_median_model, int Station_Fit, Long64_t num_entries_fit);

#endif // CPOL_H