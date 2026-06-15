// cpol_helpers.cc
//
// Definitions of the small free helper functions declared in cpol_helpers.hh.
// Moved verbatim out of cpol_main.cc to keep main() focused on the model logic.
#include "cpol_helpers.hh"
#include "constants.hh"

#include <vector>
#include <cmath>
#include <iostream>

#include "TVirtualPad.h"

using namespace std;

// --------------------------------------------------------------------------
// directionNextStep
//
// Given a refractive-index triplet nvec and a polarization-like direction p_e1,
// build a new direction vector by weighting the Cartesian components by n_i.
// The output is normalized before being returned.
//
// Note:
//   p_e2 is passed in but not used in the current implementation.
// --------------------------------------------------------------------------
TVector3 directionNextStep(vector<double>nvec, TVector3 p_e1, TVector3 p_e2) {

    p_e1.SetMag(1.);

    TVector3 temp1;
    temp1[0]=nvec[0]*p_e1[0];
    temp1[1]=nvec[1]*p_e1[1];
    temp1[2]=nvec[2]*p_e1[2];
    temp1.SetMag(1.);

    return temp1;
}


// --------------------------------------------------------------------------
// rotateE
//
// Rotate a vector E into the frame where x is aligned with the ice-flow axis,
// apply the diagonal dielectric tensor epsilon there, then rotate back to the
// original laboratory frame.
//
// In other words, this computes:
//   E_lab -> E_principal -> epsilon * E_principal -> back to lab
// --------------------------------------------------------------------------
TVector3 rotateE(TVector3 epsilon, double angle_iceflow, TVector3 E) {

    double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
        {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
        {0.,0.,1.}};

    TVector3 tempvec;
    for (int i=0;i<3;i++) {
        double sum=0.;
        for (int j=0;j<3;j++) {
            sum+=rotate_toxalongiceflow[i][j]*E[j];
        }
        tempvec[i]=sum;
    }
    E=tempvec;

    for (int i=0;i<3;i++) {
        tempvec[i]=epsilon[i]*E[i];
    }
    E=tempvec;

    double rotate_backtonormal[3][3];

    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) {
            rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
        }
    }

    for (int i=0;i<3;i++) {
        double sum=0.;

        for (int j=0;j<3;j++) {
            sum+=rotate_backtonormal[i][j]*E[j];
        }

        tempvec[i]=sum;
    }

    E=tempvec;

    return E;
}


// --------------------------------------------------------------------------
// getNewkandE
//
// Starting from:
//   - principal refractive indices nvec
//   - ice-flow direction angle
//   - mode refractive index n
//   - displacement-field-like vector D
//   - a guessed propagation direction kguess
//
// this function:
//
//   1. rotates D into the principal-axis frame
//   2. computes E = epsilon^{-1} D
//   3. constructs k from the relation k = -D + n^2 E
//   4. rotates both k and E back into the lab frame
//   5. enforces the sign of k to be aligned with kguess
//
// E is returned by reference, and k is returned as the function value.
// --------------------------------------------------------------------------
TVector3 getNewkandE(vector<double> nvec,double angle_iceflow,double n,TVector3 D,TVector3 kguess,
TVector3 &E) {

    double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
        {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
        {0.,0.,1.}};

    TVector3 tempvec;
    for (int i=0;i<3;i++) {
        double sum=0.;
        for (int j=0;j<3;j++) {
            sum+=rotate_toxalongiceflow[i][j]*D[j];
        }
        tempvec[i]=sum;
    }
    D=tempvec;

    vector<double> epsilon;
    epsilon.resize(3);

    vector<double> epsilon_inverse;
    epsilon_inverse.resize(3);

    epsilon[0]=nvec[0]*nvec[0];
    epsilon[1]=nvec[1]*nvec[1];
    epsilon[2]=nvec[2]*nvec[2];

    epsilon_inverse[0]=1/epsilon[0];
    epsilon_inverse[1]=1/epsilon[1];
    epsilon_inverse[2]=1/epsilon[2];

    E[0]=epsilon_inverse[0]*D[0];
    E[1]=epsilon_inverse[1]*D[1];
    E[2]=epsilon_inverse[2]*D[2];

    TVector3 kvec=-1.*D+n*n*E;

    if (kvec.Mag()<HOWSMALLISTOOSMALL)
    cout << "kvec mag is " << kvec.Mag() << "\n";

    kvec.SetMag(1.);

    double rotate_backtonormal[3][3];

    for (int i=0;i<3;i++) {
        for (int j=0;j<3;j++) {
            rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
        }
    }

    TVector3 tempvec2;
    for (int i=0;i<3;i++) {
        double sum1=0.;
        double sum2=0.;

        for (int j=0;j<3;j++) {
            sum1+=rotate_backtonormal[i][j]*kvec[j];
            sum2+=rotate_backtonormal[i][j]*E[j];
        }

        tempvec[i]=sum1;
        tempvec2[i]=sum2;
    }

    kvec=tempvec;
    E=tempvec2;

    if (kvec.Dot(kguess)<0.) {
        kvec=-1.*kvec;
    }

    return kvec;
}


// --------------------------------------------------------------------------
// getDeltaN(alpha)
//
// Analytic toy-model expression for delta n as a function of angle alpha
// relative to the ice-flow axis, using NICE and DELTAN as global constants.
//
// This is the simplified scalar model used for comparison plots, not the more
// general path-dependent solver used elsewhere in the main calculation.
// --------------------------------------------------------------------------
double getDeltaN(double alpha) {

    double n_L_prime=(NICE-DELTAN/2.)/sqrt(1-DELTAN/NICE*cos(alpha)*cos(alpha)*(1-DELTAN/(4*NICE)));
    return 2*(NICE-n_L_prime);
}


// --------------------------------------------------------------------------
// switchThem
//
// Swap two angles in place.
// Used in special cases (e.g. ARIANNA / biaxial mode ordering conventions)
// when the two eigenmodes need to be relabeled.
// --------------------------------------------------------------------------
void switchThem(double &thetaE_e1_Sclock,double &thetaE_e2_Sclock) {

    double tempangle=thetaE_e1_Sclock;
    thetaE_e1_Sclock=thetaE_e2_Sclock;
    thetaE_e2_Sclock=tempangle;

}

// --------------------------------------------------------------------------
// Flipped
//
// Compare the current angle theta_e1 to its starting value theta_e1_start.
// If the change is close to pi/2 (after wrapping to [-pi,pi]), interpret this
// as a sign flip / branch flip and return -1. Otherwise return +1.
//
// This is used to keep track of sign changes in the eigenvector transport along
// the ray path, so the accumulated birefringent phase can be corrected by a
// factor of ±1.
// --------------------------------------------------------------------------
double Flipped(double theta_e1,double theta_e1_start) {

    double diff=theta_e1-theta_e1_start;
    if (diff>PI)
    diff=diff-2.*PI;
    if (diff<-1.*PI)
    diff=diff+2.*PI;

    double notflipped=1.;

    if (abs(diff)>PI/2.*0.9)
    notflipped=-1.;

    return notflipped;
}

// --------------------------------------------------------------------------
// makePretty2Panel
//
// Utility function to create a two-panel canvas with customized pad geometry,
// margins, and grid settings.
// --------------------------------------------------------------------------
TCanvas * makePretty2Panel(){
    TCanvas *ccc = new TCanvas("plotEvent","plotEvent", 800, 400);
    ccc->Divide(2, 0);

    ccc->GetPad(1)->SetPad(.005, .005, .4975, .995);
    ccc->GetPad(1)->Divide(0, 2);

    ccc->GetPad(2)->SetPad(.5025, .005, .995, .995);
    ccc->SetWindowSize(1200, 700);
    ccc->cd(1)->cd(1)->SetGrid();


    gPad->SetLeftMargin(.15);
    gPad->SetBottomMargin(.12);
    ccc->cd(1)->cd(2);


    gPad->SetBottomMargin(.12);
    gPad->SetRightMargin(.19);
    gPad->SetLeftMargin(.15);
    return ccc;
}

// --------------------------------------------------------------------------
// titles(TGraph*)
//
// Helper to apply a consistent style to TGraph axis labels and titles.
// --------------------------------------------------------------------------
void titles(TGraph *inGr, TString title, TString xtitle, TString ytitle){
    auto sizeT=.055;
    inGr->SetTitle(title);
    inGr->GetXaxis()->SetTitle(xtitle);
    inGr->GetYaxis()->SetTitle(ytitle);
    inGr->GetXaxis()->SetTitleSize(sizeT);
    inGr->GetYaxis()->SetTitleSize(sizeT);
    inGr->GetXaxis()->SetLabelSize(sizeT);
    inGr->GetYaxis()->SetLabelSize(sizeT);

    inGr->GetYaxis()->SetLabelOffset(.01);
    inGr->GetYaxis()->SetTitleOffset(1.2);
}


// --------------------------------------------------------------------------
// titles(TH2*)
//
// Same helper as above, but for TH2 histogram frames used as empty axes.
// --------------------------------------------------------------------------
void titles(TH2 *inH, TString title, TString xtitle, TString ytitle){
    auto sizeT=.055;
    inH->SetTitle(title);
    inH->GetXaxis()->SetTitle(xtitle);
    inH->GetYaxis()->SetTitle(ytitle);
    inH->GetXaxis()->SetTitleSize(sizeT);
    inH->GetYaxis()->SetTitleSize(sizeT);
    inH->GetXaxis()->SetLabelSize(sizeT);
    inH->GetYaxis()->SetLabelSize(sizeT);

    inH->GetYaxis()->SetLabelOffset(.01);
    inH->GetYaxis()->SetTitleOffset(1.2);
}
