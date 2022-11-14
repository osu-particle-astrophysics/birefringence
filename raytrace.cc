/////////////////////////////////////////////////////////////////                                                                                                  
// C++ example program to demonstrate usage of birefringence code. //                                                                                              
// Comments to Amy Connolly <connolly(at)physics.osu.edu>.               //                                                                                               
/////////////////////////////////////////////////////////////////////                                                                                              
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TFile.h"
#include "TText.h"
#include "TGaxis.h"
#include "TRandom3.h"
#include "TCanvas.h"
//#include "y.hh"                                                                                                                                                  
#include <iostream>
#include <fstream>
#include <string>
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TVector2.h"
//#include "units.h"                                                                                                                                               
//#include "raytracing_tools.h"                                                                                                                                    
#include "IceRayTracing/IceRayTracing.C"
#include "TLegend.h"
#include "birefringence.hh"

using namespace std;

const double PI=3.1415926;
const double CLIGHT=3.E8;
const double DEGRAD=180./PI;
const double HOWSMALLISTOOSMALL=1.e-8;

double returnf(Double_t *val, Double_t *par);
double returnf1(Double_t *val, Double_t *par);

void findNextGuess(TVector3 raydir,double x_guess,double y_guess,TGraph2D *g2_upper,
		   double &x_nextguess,double &y_nextguess);

void finddel(double x_guess, double y_guess, TF2 *f,
	     double del[2]);

//void finddel(TVector3 raydir,double x_guess,double y_guess,vector<double> nvec,TF2 *fA, TF2 *fB, TF2 *fC,
//	       double del[2]);

TF2* findfA(vector<double> nvec,double XMIN,double XMAX,double YMIN,double YMAX);
TF2* findfB(vector<double> nvec,double XMIN,double XMAX,double YMIN,double YMAX);
TF2* findfC(vector<double> nvec,double XMIN,double XMAX,double YMIN,double YMAX);

//double findLaplacian(TVector3 raydir,double x_guess,double y_guess,TGraph2D *g2_upper);

void getabc(TF2 *fA, TF2 *fB, TF2 *fC, double x, double y,
	    double &a, double &b, double &c);

double evalUpper(TF2 *fA, TF2 *fB, TF2 *fC, double x, double y);
double evalLower(TF2 *fA, TF2 *fB, TF2 *fC, double x, double y);

double evalUpper(Double_t *val, Double_t *par);
double evalLower(Double_t *val, Double_t *par);

//TVector3 findvperp(TVector3 raydir,double x_guess,double y_guess,vector<double> nvec,TF2 *fA,TF2 *fB,TF2 *fC);
TVector3 findvperp(TVector3 raydir,double x_guess,double y_guess,vector<double> nvec);


int main(int argc, char** argv) {


  // this is the indicatrix at this depth
  vector<double> nvec_thisstep;
  nvec_thisstep.resize(3);

  // this is just for Fig. 5
  // set n_alpha, n_beta, and n_gamma here.
  // these are to represent shallow ice (200m)
  nvec_thisstep[0]=1.7785;
  nvec_thisstep[1]=1.779;
  nvec_thisstep[2]=1.7815;

  // these are to represent deep ice (-1700m)
  /* 
  nvec_thisstep[0]=1.777;
  nvec_thisstep[1]=1.78125;
  nvec_thisstep[2]=1.78175;
  */

  /*
  nvec_thisstep[0]=1.1;
  nvec_thisstep[1]=1.3;
  nvec_thisstep[2]=1.5;
  */

  // this is the velocity if c=1.
  double v[3]={0.};
  for (int i=0;i<3;i++) {
    v[i]=1./nvec_thisstep[i];
  }
  
  // This is V_z, which gives the direction of the optic axis
  // this angle is wrt the z axis
  double Vz=getV(nvec_thisstep);

  //  cout << "Vz is " << Vz << "\n";

  // Eq. 2 of Latorre et al. is of the form Az^4+Bz^2+C.
  // where A, B, and C are all functions of x and y.
  // the functions below represent A, B, and C.

  
  double MINX=-1.0;
  double MAXX=1.0;
  double MINY=-1.0;
  double MAXY=1.0;


  TF2 *fA = findfA(nvec_thisstep,MINX,MAXX,MINY,MAXY);
  TF2 *fB = findfB(nvec_thisstep,MINX,MAXX,MINY,MAXY);
  TF2 *fC = findfC(nvec_thisstep,MINX,MAXX,MINY,MAXY);
  

  TGraph2D *g2_upper=new TGraph2D();
  TGraph2D *g2_lower=new TGraph2D();

  TGraph2D *g2_sqrt=new TGraph2D();

  TGraph2D *g2_f=new TGraph2D();

  TGraph2D *g2_a=new TGraph2D();


  const int NX=100;
  const int NY=100;

  double stepx=(MAXX-MINX)/(double)NX;
  double stepy=(MAXY-MINY)/(double)NY;

  int index=0;

  TVector3 raydir[2]; // follow two solutions?

  raydir[0].SetXYZ(0.2,0.2,0.2); // this is just chosen at random for now
  raydir[0].SetMag(1.); // make it a unit vector

  // this loop just fills some graphs
  for (int i=0;i<NX;i++) {
    for (int j=0;j<NY;j++) {

      double x=MINX+(double)i*stepx;
      double y=MINY+(double)j*stepy;
      
      double a,b,c;

      getabc(fA,fB,fC,x,y,
	     a,b,c);
      
      if (b*b>4*a*c) {
	//cout << "x, y, a, b, c are " << x << "\t" << y << "\t" << a << "\t" << b << "\t" << c << "\n";
	if ((-1.*b-sqrt(b*b-4*a*c))/(2*a)>0.) {
	  g2_upper->SetPoint(index,x,y,sqrt((-1.*b+sqrt(b*b-4*a*c))/(2*a)));
	  g2_lower->SetPoint(index,x,y,sqrt((-1.*b-sqrt(b*b-4*a*c))/(2*a)));
	  g2_sqrt->SetPoint(index,x,y,sqrt(b*b-4*a*c));
	  
	  index++;
	  //cout << "upper, lower are " << ((-1.*b+sqrt(b*b-4*a*c))/(2*a)) << "\t" << ((-1.*b-sqrt(b*b-4*a*c))/(2*a)) << "\n";
	}
      }

    }
  }
  // this loop also just fills some graphs
  index=0;

  // this function is going to be our figure of merit that we are trying to minimize
  // it is the magnitude of the cross product between the vector perp to the surface and raydir[0].
  TF2 *func=new TF2("func",returnf,MINX,MAXX,MINY,MAXY,6,2);

  func->SetParameter(3,raydir[0][0]);
  func->SetParameter(4,raydir[0][1]);
  func->SetParameter(5,raydir[0][2]);
  func->SetParameter(0,nvec_thisstep[0]);
  func->SetParameter(1,nvec_thisstep[1]);
  func->SetParameter(2,nvec_thisstep[2]);


  for (int i=0;i<NX;i++) {
    for (int j=0;j<NY;j++) {

      double x=MINX+(double)i*stepx;
      double y=MINY+(double)j*stepy;

      double a,b,c;

      getabc(fA,fB,fC,x,y,
	     a,b,c);
      /*
      cout << "b^2-4ac is " << b*b-4*a*c << "\n";
      cout << "(-b+sqrt(b^2-4ac))/(2a) is " << (-1.*b-sqrt(b*b+4*a*c))/(2*a) << "\n";
      cout << "(-b-sqrt(b^2-4ac))/(2a) is " << (-1.*b-sqrt(b*b-4*a*c))/(2*a) << "\n";
      */
    
      if (b*b>4*a*c) {
	if ((-1.*b-sqrt(b*b-4*a*c))/(2*a)>0.) {
	  
	  g2_f->SetPoint(index,x,y,func->Eval(x,y));
	  cout << "x, y, returnf are " << x << "\t" << y << "\t" << func->Eval(x,y) << "\n";
	  g2_a->SetPoint(index,x,y,a);
	  index++;
	}
      }

    }
  }


  // this is the x where the two surfaces meet (b^2-4ac=0)
  double x_min=v[2]*sqrt((v[1]*v[1]-v[0]*v[0])/(v[2]*v[2]-v[0]*v[0])); // this is from my math, where b^2-4ac=0
  cout << "I think x should be " << x_min << "\n";  
  cout << "sqrt at x_min,0 is " << g2_sqrt->Interpolate(x_min,0.) << "\n"; // this should be close to zero then because g2_sqrt is plotting sqrt(b^2-4ac)
  // this is 90.-atan(z/x) at x_min which should be V_z
  cout << "angle where they intersect is " << 90.-atan2((g2_upper->Interpolate(x_min,0.)),x_min)*DEGRAD << "\n"; // this should be close to V.

  

  // starting some ray tracing
  // we're going to start with an initial guess raydir[0] for the starting ray direction

  // we are going to find the place on the surface where raydir[0] is perp to the surface.
  // this is just an initial guess for the x and y where this condition is satisfied.
  // we are going to use the method of steepest descent
  double x_guess=0.0;
  double y_guess=0.0;

  TVector3 vperp = findvperp(raydir[0],x_guess,y_guess,nvec_thisstep); // this is the vector perpendicular to the surface

  cout << "for the guess, f is " << func->Eval(x_guess,y_guess) << "\n";

  // how many tries.
  // maybe this should instead stop after some condition is met.
  int NTRY=1000;
  for (int itry=0;itry<NTRY;itry++) { // come up with some condition for when to find a next guess

    // these next few lines get the next guess using the method of gradient descent 
    double del[2]={0.};

    finddel(x_guess,y_guess,func,
	    del); // find the del of the function

    double alpha=0.002;

    x_guess=x_guess-alpha*del[0]; // x_j+1 = x_j - alpha*del_x(f)
    y_guess=y_guess-alpha*del[1]; // y_j+1 = y_j - alpha*del_y(f)
      
    cout << "del is " << del[0] << "\t" << del[1] << "\n";
    cout << "x_nextguess is " << x_guess << "\n";
    cout << "y_nextguess is " << y_guess << "\n";

    cout << "for the next guess, f is " << func->Eval(x_guess,y_guess) << "\n"; // this is the figure of merit

    // this is just for printint out and checking the answer makes sense.
    vperp = findvperp(raydir[0],x_guess,y_guess,nvec_thisstep); // this is the vector perpendicular to the surface
    cout << "vperp for this guess is " << vperp[0] << "\t" << vperp[1] << "\t" << vperp[2] << "\n";
    cout << "theta, phi for this guess is " << vperp.Theta()*180./PI << "\t" << vperp.Phi()*180./PI << "\n";


  }

  cout << "vperp is " << vperp[0] << "\t" << vperp[1] << "\t" << vperp[2] << "\n";

  // make this next part into a function called findpointB
  // find point B.
  TVector3 zvec(0.,0.,1.);
  TVector3 intothepage=zvec.Cross(vperp);
  intothepage.SetMag(1.); // this points into the page

  cout << "intothepage is " << intothepage[0] << "\t" << intothepage[1] << "\t" << intothepage[2] << "\n";

  TVector3 towardsB=raydir[0].Cross(intothepage);
  towardsB.SetMag(1.); // this points from the tangent point to point B

  cout << "towardsB is " << towardsB[0] << "\t" << towardsB[1] << "\t" << towardsB[2] << "\n";

  TVector3 pointofinterest(x_guess,y_guess,evalLower(fA, fB, fC, x_guess, y_guess));

  // this is where the plane wave hits the surface
  cout << "pointofinterest is " << pointofinterest[0] << "\t" << pointofinterest[1] << "\t" << pointofinterest[2] << "\n";
    
  double xpointB = pointofinterest[0]-pointofinterest[2]*towardsB[0]/towardsB[2];
  double ypointB = pointofinterest[1]-pointofinterest[2]*towardsB[1]/towardsB[2];

  TVector3 pointB(xpointB,ypointB,0.);

  cout << "pointB is " << pointB[0] << "\t" << pointB[1] << "\t" << pointB[2] << "\n";

  cout << "phi of pointB is " << pointB.Phi()*DEGRAD << "\n";

  // for the next step, grab new indicatrix, then
  // find the new wave surfaces,
  // then find the new transmitted wave: line connecting point of interest to point B should be perp to plane of incidence, tangent to wave surface
  // loop over x, y
 

  // loop over theta_T (transmitted theta), and angle phi (angle ray makes with plane of incidence in plane defined by theta_T).
  // use method of gradient descent to find where line from center of the wave surface to the point on the wave surface makes a 90deg. angle with the line from the point on the wave surface to point B.

  int NTHETAT=100;
  double THETATMIN=0.;
  double THETATMAX=90.;
  double STEPTHETAT=(THETATMAX-THETATMIN)/(double)NTHETAT;

  int NPHIT=100;
  double PHITMIN=0.;
  double PHITMAX=0.;
  double STEPPHIT=(PHITMAX-PHITMIN)/(double)NPHIT;

  for (int ithetat=0;ithetat<NTHETAT;ithetat++) {
      double thetat=THETATMIN+(double)ithetat*STEPTHETAT;
    for (int iphit=0;iphit<NPHIT;iphit++) {
      double phit=PHITMIN+(double)iphit*STEPPHIT;


    }
  }


  TCanvas *c2=new TCanvas("c2","c2",800,800);


  c2->Divide(2,2);


  TH2D *h2=new TH2D("h2","h2",100,MINX,MAXX,100,MINY,MAXY);
  h2->SetMaximum(0.1);

  gStyle->SetPalette(1);

  c2->cd(1);
  gPad->SetTheta(75.);
  gPad->SetPhi(45.);
  g2_upper->Draw("surf1");
  

  c2->cd(2);
  gPad->SetTheta(75.);
  gPad->SetPhi(45.);
  g2_lower->Draw("surf1");

  c2->cd(3);
  gPad->SetTheta(75.);
  gPad->SetPhi(45.);

  //h2->Draw("surf1");
  g2_upper->SetMarkerStyle(20);
  g2_upper->SetMarkerColor(kRed);
  g2_upper->SetMarkerSize(1.0);
  //g2_f->Draw("surf1");
  g2_f->Draw("colz");
  //g2_upper->Draw("colz");
//  g2_upper->Draw("colz");

  c2->cd(4);
  gPad->SetTheta(75.);
  gPad->SetPhi(65.);
  //h2->Draw("surf1");
  g2_lower->SetMarkerStyle(20);
  g2_lower->SetMarkerColor(kBlue);
  g2_lower->SetMarkerSize(1.0);
  //g2_lower->Draw("colz");
  //  g2_lower->Draw("cont1");
  //g2_lower->Draw("surf1");

  // g2_sqrt->Draw("colz");
  //g2_sqrt->Draw("surf1");
  //g2_a->Draw("surf1");
  c2->Print("c2.png");


}
double evalUpper(TF2 *fA, TF2 *fB, TF2 *fC, double x, double y) {
  
  double a,b,c;
  getabc(fA, fB, fC, x, y,
	 a, b, c);

  //return (-1.*b+sqrt(b*b-4*a*c))/(2*a);
  if (b*b-4*a*c>0.)
    return (-1.*b+sqrt(b*b-4*a*c))/(2*a);
  else
    return -1.;

}
double evalUpper(Double_t *par, Double_t *var) {
  
  double x_guess=var[0];
  double y_guess=var[1];

  vector<double> nvec{par[0],par[1],par[2]};

  double MINX=par[3];
  double MAXX=par[4];
  double MINY=par[5];
  double MAXY=par[6];
  
  TF2* fA=findfA(nvec,MINX,MAXX,MINY,MAXY);
  TF2* fB=findfB(nvec,MINX,MAXX,MINY,MAXY);
  TF2* fC=findfC(nvec,MINX,MAXX,MINY,MAXY);

  double a,b,c;
  getabc(fA, fB, fC, x, y,
	 a, b, c);

  //return (-1.*b+sqrt(b*b-4*a*c))/(2*a);
  if (b*b-4*a*c>0.)
    return (-1.*b+sqrt(b*b-4*a*c))/(2*a);
  else
    return -1.;

}
double evalLower(TF2 *fA, TF2 *fB, TF2 *fC, double x, double y) {
  
  double a,b,c;
  getabc(fA, fB, fC, x, y,
	 a, b, c);

  if (b*b-4*a*c>0.)
    return (-1.*b-sqrt(b*b-4*a*c))/(2*a);
  else
    return -1.;

}

void getabc(TF2 *fA, TF2 *fB, TF2 *fC, double x, double y,
	    double &a, double &b, double &c) {
  a=fA->Eval(x,y);
  b=fB->Eval(x,y);
  c=fC->Eval(x,y);
}
void findNextGuess(TVector3 raydir,double x_guess,double y_guess,TGraph2D *g2_upper,
		   double &x_nextguess,double &y_nextguess) {

}
/*
double findLaplacian(TVector3 raydir,double x_guess,double y_guess,TGraph2D *g2_upper) {

  return 0.;
}
*/
void finddel(double x_guess, double y_guess, TF2 *f,
	     double del[2]) {

  double xstep=0.02;
  double ystep=0.02;
  
  double xnext=x_guess+xstep/2.;
  double ynext=y_guess+ystep/2.;
  
  double xprevious=x_guess-xstep/2.;
  double yprevious=y_guess-ystep/2.;
  
  cout << "xnext, xprevious, next f, previous f, diffx are " << xnext << "\t" << xprevious << "\t" << f->Eval(xnext,y_guess) << "\t" << f->Eval(xprevious,y_guess) << "\t" << xnext-xprevious << "\n";

  del[0]=(f->Eval(xnext,y_guess)-f->Eval(xprevious,y_guess))/(xnext-xprevious);
  del[1]=(f->Eval(x_guess,ynext)-f->Eval(x_guess,yprevious))/(ynext-yprevious);

}
/*
void finddel(TVector3 raydir, double x_guess, double y_guess, vector<double> nvec,TF2 *fA, TF2 *fB, TF2 *fC,
	       double del[2]) {

  double xstep=0.02;
  double ystep=0.02;
  
  double xnext=x_guess+xstep/2.;
  double ynext=y_guess+ystep/2.;
  
  double xprevious=x_guess-xstep/2.;
  double yprevious=y_guess-ystep/2.;
  
  cout << "xnext, xprevious, next f, previous f, diffx are " << xnext << "\t" << xprevious << "\t" << returnf(raydir,xnext,y_guess,nvec,fA,fB,fC) << "\t" << returnf(raydir,xprevious,y_guess,nvec,fA,fB,fC) << "\t" << xnext-xprevious << "\n";

  del[0]=(returnf(raydir,xnext,y_guess,nvec)-returnf(raydir,xprevious,y_guess,nvec,fA,fB,fC))/(xnext-xprevious);
  del[1]=(returnf(raydir,x_guess,ynext,nvec)-returnf(raydir,x_guess,yprevious,nvec,fA,fB,fC))/(ynext-yprevious);

}
*/
TF2* findfA(vector<double> nvec,double MINX,double MAXX,double MINY,double MAXY) {

  double v[3]={0.};

  for (int i=0;i<3;i++) {
    v[i]=1./nvec[i];
  }

  TF2 *fA = new TF2("fA","[0]*[0]",MINX,MAXX,MINY,MAXY);
  fA->SetParameter(0,v[2]);

  return fA;

}
TF2* findfB(vector<double> nvec,double MINX,double MAXX,double MINY,double MAXY) {

  double v[3]={0.};
  
  for (int i=0;i<3;i++) {
    v[i]=1./nvec[i];
  }
  
  TF2 *fB = new TF2("fB","[1]*[1]*x*x+[2]*[2]*y*y+[0]*[0]*(x*x+y*y)-[0]*[0]*([1]*[1]+[2]*[2])",MINX,MAXX,MINY,MAXY);
  fB->SetParameter(0,v[2]);
  fB->SetParameter(1,v[0]);
  fB->SetParameter(2,v[1]);

  return fB;

}
TF2* findfC(vector<double> nvec,double MINX,double MAXX,double MINY,double MAXY) {

  double v[3]={0.};

  for (int i=0;i<3;i++) {
    v[i]=1./nvec[i];
  }

  TF2 *fC = new TF2("fC","(x*x+y*y)*([1]*[1]*x*x+[2]*[2]*y*y)-[1]*[1]*([2]*[2]+[0]*[0])*x*x-[2]*[2]*([0]*[0]+[1]*[1])*y*y+[1]*[1]*[2]*[2]*[0]*[0]",MINX,MAXX,MINY,MAXY);
  fC->SetParameter(0,v[2]);
  fC->SetParameter(1,v[0]);
  fC->SetParameter(2,v[1]);

  return fC;



}

TVector3 findvperp(TVector3 raydir,double x_guess,double y_guess,vector<double> nvec) {

  double MINX=x_guess-1.;
  double MAXX=x_guess+1.;
  double MINY=y_guess-1.;
  double MAXY=y_guess+1.;

  TF2 *fA=findfA(nvec,MINX,MAXX,MINY,MAXY);
  TF2 *fB=findfB(nvec,MINX,MAXX,MINY,MAXY);
  TF2 *fC=findfC(nvec,MINX,MAXX,MINY,MAXY);

  double z_guess=evalLower(fA,fB,fC,x_guess,y_guess);
  
  double vx=nvec[0];
  double vy=nvec[1];
  double vz=nvec[2];

  double dgdx=4*vx*vx*x_guess*x_guess*x_guess + 2*vy*vy*y_guess*y_guess*x_guess + 2*vz*vz*z_guess*z_guess*x_guess 
    + 2*vx*vx*y_guess*y_guess*x_guess + 2*vx*vx*z_guess*z_guess*x_guess - 2*vx*vx*(vy*vy + vz*vz)*x_guess;

  double dgdy=2*vy*vy*x_guess*x_guess*y_guess + 2*vx*vx*x_guess*x_guess*y_guess + 4*vy*vy*y_guess*y_guess*y_guess
    + 2*vz*vz*z_guess*z_guess*y_guess + 2*vy*vy*z_guess*z_guess*y_guess - 2*vy*vy*(vx*vx+vz*vz)*y_guess;

  double dgdz=2*vz*vz*x_guess*x_guess*z_guess + 2*vz*vz*y_guess*y_guess*z_guess + 2*vx*vx*x_guess*x_guess*z_guess 
    + 2*vy*vy*y_guess*y_guess*z_guess + 4*vz*vz*z_guess*z_guess*z_guess - 2*vz*vz*(vx*vx+vy*vy)*z_guess;

  TVector3 vperp(dgdx,dgdy,dgdz);


  vperp.SetMag(1.);

  return vperp;

}

double returnf(Double_t *val, Double_t *par) {
  // now find x and y such that the tangent plane to the g2_upper curve perp to raydir
  
  Double_t x_guess = val[0];
  Double_t y_guess = val[1];

  vector<double> nvec; // the three principal axes
  nvec.resize(3);

  nvec[0]=par[0];
  nvec[1]=par[1];
  nvec[2]=par[2];

  TVector3 raydir; // direction of the incident ray
  raydir.SetXYZ(par[3],par[4],par[5]);


  TVector3 vperp=findvperp(raydir,x_guess,y_guess,nvec); // direction perpendicular to contour
  
  double dot=vperp.Dot(raydir)/vperp.Mag()/raydir.Mag(); // cos of angle between vperp and raydir
  
  //  cout << "dot product is " << dot << "\n";
  
  double f=sqrt(1-dot*dot); // sin of the angle between vperp and raydir.  want to minimize this quantity.
  
  /*
  cout << "x_guess, y_guess, dot, f are " << x_guess << "\t" << y_guess << "\t" << dot << "\t" << f << "\n";
  cout << "g2_upper is " << evalUpper(fA,fB,fC,x_guess,y_guess) << "\n";
  cout << "vperp is " << vperp[0] << "\t" << vperp[1] << "\t" << vperp[2] << "\n";
  cout << "raydir is " << raydir[0] << "\t" << raydir[1] << "\t" << raydir[2] << "\n";
  */
  
  return f;
}
double returnf2(Double_t *val, Double_t *par) {
  // now find x and y such that the tangent plane to the g2_upper curve perp to raydir
  
  Double_t x_guess = val[0];
  Double_t y_guess = val[1];


  TF2 *func=new TF2("func",returnf,MINX,MAXX,MINY,MAXY,6,2);

  func->SetParameter(3,par[3]); // these are the direction of the ray
  func->SetParameter(4,par[4]);
  func->SetParameter(5,par[5]);

  func->SetParameter(0,par[0]); // these are the three principal axes
  func->SetParameter(1,par[1]);
  func->SetParameter(2,par[2]);

  // do we really want to evaluate this here?
  double f=func->Eval(x_guess,y_guess); // this is the quantity to be minimized
  f=f*f;

  TVector3 pointB;
  pointB.SetXYZ(par[6],par[7],par[8]);

  // need to turn evalUpper into a function
  // these is the higher of the two closed surfaces
  TF2 *f_upper=new TF2("f_upper",evalUpper,MINX,MAXX,MINY,MAXY,7,2);
  f_upper->SetParameter(0,par[3]); // these are the direction of the ray?  or should they be the principle axes
  f_upper->SetParameter(1,par[4]);
  f_upper->SetParameter(2,par[5]);

  TVector3 vperp=findvperp(raydir,x_guess,y_guess,nvec);
  
  double dot=vperp.Dot(raydir)/vperp.Mag()/raydir.Mag();
  
  //  cout << "dot product is " << dot << "\n";
  
  double f=sqrt(1-dot*dot);  // this looks like the thing we want to minimize
  
  /*
  cout << "x_guess, y_guess, dot, f are " << x_guess << "\t" << y_guess << "\t" << dot << "\t" << f << "\n";
  cout << "g2_upper is " << evalUpper(fA,fB,fC,x_guess,y_guess) << "\n";
  cout << "vperp is " << vperp[0] << "\t" << vperp[1] << "\t" << vperp[2] << "\n";
  cout << "raydir is " << raydir[0] << "\t" << raydir[1] << "\t" << raydir[2] << "\n";
  */
  
  return f;
}
