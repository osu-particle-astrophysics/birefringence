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
#include "IceRayTracing-master/IceRayTracing.h"
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
double getV(vector<double> nvec) {

  double alpha=nvec[0];
  double beta=nvec[1];
  double gamma=nvec[2];

  //  return 90.-acos(sqrt((gamma*gamma)*(beta*beta-alpha*alpha)/((beta*beta)*(gamma*gamma-alpha*alpha)  )))*DEGRAD;
  return 90.-acos(((gamma*gamma)*(beta*beta-alpha*alpha)/((beta*beta)*(gamma*gamma-alpha*alpha)  )))*DEGRAD;

}
void thetastoEpsilons(double thetaE_e1_Sclock,double thetaE_e2_Sclock,
		      double &epsilon1,double &epsilon2) {

  epsilon1=thetaE_e1_Sclock-PI/2.;
  epsilon2=thetaE_e2_Sclock;

}
void getAnglesontheClock(TVector3 rhat_thisstep1, TVector3 rhat_thisstep2,
                         TVector3 p_e1, TVector3 p_e2,
                         double &theta_e1,double &theta_e2,
                         double &p_e1_thetacomponent,double &p_e2_thetacomponent,
                         double &p_e1_phicomponent,double &p_e2_phicomponent) {

  double theta_I_1=rhat_thisstep1.Theta();
  double theta_I_2=rhat_thisstep2.Theta();

  TVector3 zaxis(0.,0.,1.);
  TVector3 yaxis(0.,1.,0.);
  // rhat_thisstep is the direction of khat.                                                                                        
  // rotate khat to be pointing in the +z direction,                                                                                
  // and rotate p_e1 and p_e2 along with them.                                                                                      
  double phi1=-1.*rhat_thisstep1.Phi();
  rhat_thisstep1.Rotate(phi1,zaxis);
  p_e1.Rotate(phi1,zaxis);

  double phi2=-1.*rhat_thisstep2.Phi();
  rhat_thisstep2.Rotate(phi2,zaxis);
  p_e2.Rotate(phi2,zaxis);

  double theta1=-1.*rhat_thisstep1.Theta();
  rhat_thisstep1.Rotate(theta1,yaxis);
  p_e1.Rotate(theta1,yaxis);

  double theta2=-1.*rhat_thisstep2.Theta();
  rhat_thisstep2.Rotate(theta2,yaxis);
  p_e2.Rotate(theta2,yaxis);

  // now find the angles on the clock that p_e1 and p_e2 sit at.                                                                    
  theta_e1=p_e1.Phi();
  theta_e2=p_e2.Phi();

  double epsilon1=0.;
  double epsilon2=0.;
  thetastoEpsilons(theta_e1,theta_e2,
                   epsilon1,epsilon2);


  double scalefactor=0.2/(sin(theta_I_1)+sin(theta_I_2))*2.;

  p_e1_thetacomponent=scalefactor*sin(theta_I_1)*cos(epsilon1)*cos(epsilon1);
  p_e2_thetacomponent=scalefactor*sin(theta_I_2)*sin(epsilon2)*sin(epsilon2);

  p_e1_phicomponent=scalefactor*sin(theta_I_1)*sin(epsilon1)*cos(epsilon1);
  p_e2_phicomponent=scalefactor*sin(theta_I_2)*-1.*cos(epsilon2)*sin(epsilon2);

}

void getManyAnglesontheClock(int BIAXIAL,double crosspolangle_tx,
                             TVector3 rhat_thisstep,
                             TVector3 p_e1,TVector3 p_e2,
                             TVector3 E_e1,TVector3 E_e2,
                             double &theta_e1,double &theta_e2,
                             double &thetaE_e1,double &thetaE_e2,
                             double &theta_e1_Sclock,double &theta_e2_Sclock,
                             double &thetaE_e1_Sclock,double &thetaE_e2_Sclock,
                             TVector3 &Shat_e1,TVector3 &Shat_e2,
                             double &E_e1_thetacomponent_Sclock,double &E_e2_thetacomponent_Sclock,
                             double &E_e1_phicomponent_Sclock,double &E_e2_phicomponent_Sclock) {

  // these aren't currently used for anything but they are the theta and phi components of D-hat for each ray, on the clock where k is in the page                                                                                                                     
   double p_e1_thetacomponent=0.;
   double p_e2_thetacomponent=0.;
   double p_e1_phicomponent=0.;
   double p_e2_phicomponent=0.;
   getAnglesontheClock(rhat_thisstep,rhat_thisstep, // we are using the same k vector for both rays                                  
		       p_e1,p_e2,
		       theta_e1,theta_e2,
		       p_e1_thetacomponent,p_e2_thetacomponent,
		       p_e1_phicomponent,p_e2_phicomponent);
   
   if (theta_e1<0.)
     theta_e1+=PI;
   if (theta_e2<-1.*PI/2.)
     theta_e2+=PI;
   if (theta_e2>PI/2.)
     theta_e2-=PI;
   
   
   theta_e1+=crosspolangle_tx;
   theta_e2+=crosspolangle_tx;
   
   // these aren't currently used for anything but they are the theta and phi components of the E field for each ray, on the clock where k is in the page                                                                                                               
   double E_e1_thetacomponent=0.;
   double E_e2_thetacomponent=0.;
   double E_e1_phicomponent=0.;
   double E_e2_phicomponent=0.;
   
   getAnglesontheClock(rhat_thisstep,rhat_thisstep,
		       E_e1,E_e2,
		       thetaE_e1,thetaE_e2,
		       E_e1_thetacomponent,E_e2_thetacomponent,
		       E_e1_phicomponent,E_e2_phicomponent);
   
   if (thetaE_e1<0.)
     thetaE_e1+=PI;
   if (thetaE_e2<-1.*PI/2.)
     thetaE_e2+=PI;
   if (thetaE_e2>PI/2.)
     thetaE_e2-=PI;
   
   thetaE_e1+=crosspolangle_tx;
   thetaE_e2+=crosspolangle_tx;
   
   // 09/05/21 for uniaxial, problem is that p and E are parallel I think.                                                           
   TVector3 Hhat_e1;
   TVector3 Hhat_e2;
   
   // for an isotropic medium, I think Hhat_e1 keeps getting flipped                                                                
   // back and forth                                                                                                                 
   // this part makes sure it stays oriented the say way relative to                                                                 
   // each eigenvector and the direction of the ray.                                                                                 
   if (BIAXIAL==-1) {

     Hhat_e1=rhat_thisstep.Cross(p_e1);
     Hhat_e2=rhat_thisstep.Cross(p_e2);
     
     if (Hhat_e1.Mag()<HOWSMALLISTOOSMALL)
       cout << "Hhat_e1 is " << Hhat_e1.Mag() << "\n";
     
     Hhat_e1.SetMag(1);
     
     if (Hhat_e2.Mag()<HOWSMALLISTOOSMALL) {
       cout << "Shat_e2 is " << Shat_e2.Mag() << "\n";
       cout << "p_e2 is " << p_e2.Mag() << "\n";
       cout << "Hhat_e2 is " << Hhat_e2.Mag() << "\n";
     }
     Hhat_e2.SetMag(1);
     
     Shat_e1=E_e1.Cross(Hhat_e1);
     Shat_e2=E_e2.Cross(Hhat_e2);
     
     Shat_e1.SetMag(1.);
     Shat_e2.SetMag(1.);
   }
   else if (BIAXIAL==0) {
     
     Hhat_e1=rhat_thisstep.Cross(p_e1);
     Hhat_e2=rhat_thisstep.Cross(E_e2);
     
     Hhat_e1.SetMag(1.);
     Hhat_e2.SetMag(1.);
     Shat_e1=E_e1.Cross(Hhat_e1);
     Shat_e2=E_e2.Cross(Hhat_e2);
     
     Shat_e1.SetMag(1.);
     Shat_e2.SetMag(1.);
     
     
   }
   else {
     
     Hhat_e1=rhat_thisstep.Cross(p_e1);
     Hhat_e2=rhat_thisstep.Cross(p_e2);
     
     
     if (Hhat_e1.Mag()<HOWSMALLISTOOSMALL) {
       cout << "E_e1 is " << E_e1[0] << "\t" << E_e1[1] << "\t" << E_e1[2] << "\n";
       cout << "p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
       cout << "Hhat_e1 is " << Hhat_e1.Mag() << "\n";
     }
     
     
     Hhat_e1.SetMag(1.);
     Hhat_e2.SetMag(1.);
     Shat_e1=E_e1.Cross(Hhat_e1);
     Shat_e2=E_e2.Cross(Hhat_e2);
     
     Shat_e1.SetMag(1.);
     Shat_e2.SetMag(1.);
     
   }
   

   // here i want to plot where p_e1 and p_2 are on the clock, with 12 o'clock being the in the plane of rhat at launch and the z axis.                                                                                                                                 
   // these aren't currently used for anything but they are the theta and phi components of the D-hat eigenvector for each ray, on the clock where S is in the page                                                                                                     
   double p_e1_thetacomponent_Sclock=0.;
   double p_e2_thetacomponent_Sclock=0.;
   double p_e1_phicomponent_Sclock=0.;
   double p_e2_phicomponent_Sclock=0.;
   
   getAnglesontheClock(Shat_e1,Shat_e2,
		       p_e1,p_e2,
		       theta_e1_Sclock,theta_e2_Sclock,
		       p_e1_thetacomponent_Sclock,p_e2_thetacomponent_Sclock,
		       p_e1_phicomponent_Sclock,p_e2_phicomponent_Sclock);
   

   if (theta_e1_Sclock<0.)
     theta_e1_Sclock+=PI;
   if (theta_e2_Sclock<-1.*PI/2.)
     theta_e2_Sclock+=PI;
   if (theta_e2_Sclock>PI/2.)
     theta_e2_Sclock-=PI;
   

   
   theta_e1_Sclock+=crosspolangle_tx;
   theta_e2_Sclock+=crosspolangle_tx;
   
   getAnglesontheClock(Shat_e1,Shat_e2,
		       E_e1,E_e2,
		       thetaE_e1_Sclock,thetaE_e2_Sclock,
		       E_e1_thetacomponent_Sclock,E_e2_thetacomponent_Sclock,
		       E_e1_phicomponent_Sclock,E_e2_phicomponent_Sclock);
   
   if (thetaE_e1_Sclock<0.)
     thetaE_e1_Sclock+=PI;
   if (thetaE_e2_Sclock<-1.*PI/2.)
     thetaE_e2_Sclock+=PI;
   if (thetaE_e2_Sclock>PI/2.)
     thetaE_e2_Sclock-=PI;
   
   thetaE_e1_Sclock+=crosspolangle_tx;
   thetaE_e2_Sclock+=crosspolangle_tx;
   
   
}




TVector3 rotateD(TVector3 epsilon, double angle_iceflow, TVector3 D) {
  
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
  
  TVector3 inverseepsilon(1./epsilon[0],1./epsilon[1],1./epsilon[2]);

  for (int i=0;i<3;i++) {
    tempvec[i]=inverseepsilon[i]*D[i];
  }
  D=tempvec;
  
  double rotate_backtonormal[3][3];
  
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
    }
  }
  for (int i=0;i<3;i++) {
    double sum=0.;
    
    for (int j=0;j<3;j++) {
      sum+=rotate_backtonormal[i][j]*D[j];
    }
    
    tempvec[i]=sum;
  }
  
  D=tempvec;
  
    return D; // this is actually returning an electric field                                                                         
    
    
}

double getDeltaN(int BIAXIAL,vector<double> nvec,TVector3 rhat,double angle_iceflow, double &n_e1, double &n_e2,TVector3 &p_e1,TVector3 &p_e2) {                   
                                                                                                                                                                   
  int FLIPPED=0;                                                                                                                                                   
  
  if (rhat[2]<0.) {                                                                                                                                                
    FLIPPED=1;                                                                                                                                                     
    rhat[2]=-1.*rhat[2];                                                                                                                                           
  }                 
  
  TVector3 myy;
  myy[0]=0.;
  myy[1]=-1.;
  myy[2]=0.;
  
  TVector3 myz;
  myz[0]=0.;
  myz[1]=0.;
  myz[2]=1.;


  double phi_rhat=atan2(rhat[1],rhat[0]);
  double phi_wrticeflow=angle_iceflow-phi_rhat;
  if (phi_wrticeflow<-PI)
    phi_wrticeflow+=2.*PI;
  if (phi_wrticeflow>PI)
    phi_wrticeflow-=2.*PI;

  double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
                                       {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
                                       {0.,0.,1.}};

  TVector3 rhat_iceflowalongx;
  TVector3 myy_iceflowalongx;
  TVector3 myz_iceflowalongx;
  TVector3 x_iceflowalongx(1.,0.,0.);
  TVector3 y_iceflowalongx(0.,1.,0.);
  TVector3 z_iceflowalongx(0.,0.,1.);


  for (int i=0;i<3;i++) {
    double sum=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_toxalongiceflow[i][j]*rhat[j];
    }
    rhat_iceflowalongx[i]=sum;
  }



  TVector3 nominal_pe1=myz.Cross(rhat_iceflowalongx);



  if (nominal_pe1.Mag()<HOWSMALLISTOOSMALL) {
    cout << "myz is " << myz[0] << "\t" << myz[1] << "\t" << myz[2] << "\n";
    cout << "rhat is " << rhat[0] << "\t" << rhat[1] << "\t" << rhat[2] << "\n";
    cout << "rhat_iceflowalongx is " << rhat_iceflowalongx[0] << "\t" << rhat_iceflowalongx[1] << "\t" << rhat_iceflowalongx[2] << "\n";
    cout << "cross of them is " << nominal_pe1[0] << "\t" << nominal_pe1[1] << "\t" << nominal_pe1[2] << "\n";
    cout << "nominal_pe1 mag is " << nominal_pe1.Mag() << "\n";
  }
  nominal_pe1.SetMag(1.);
  TVector3 nominal_pe2=rhat_iceflowalongx.Cross(nominal_pe1);



  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_toxalongiceflow[i][j]*myy[j];
      sum2+=rotate_toxalongiceflow[i][j]*myz[j];
    }
    myy_iceflowalongx[i]=sum;
    myz_iceflowalongx[i]=sum2;
  }


  double a=nvec[0];
  double b=nvec[1];
  double c=nvec[2];


  double Ax=rhat_iceflowalongx[0];
  double Ay=rhat_iceflowalongx[1];
  double Az=rhat_iceflowalongx[2];


  double A=1/(a*a)+(Ax*Ax)/(Az*Az*c*c);
  double B=(2.*Ax*Ay)/(Az*Az*c*c);
  double C=1/(b*b)+(Ay*Ay)/(Az*Az*c*c);
  //double F=-1.;

  /*
  double M=A;
  double P=B/2.;
  double Q=B/2.;
  double R=C;
  */

  double theta_initial=atan2( B , A - C )/2.; // not sure this is rotated in the right direction - check this.                                                     

  //  double lambda2=(1.*(M+R)+sqrt((M-R)*(M-R)+4*P*Q))/2.;
  //double lambda1=(1.*(M+R)-sqrt((M-R)*(M-R)+4*P*Q))/2.;

  // these are only the n's for the scenario where the plane arrives straight from above, a test scenario                                                          
  //double ne2=sqrt(-1./lambda2);
  //double ne1=sqrt(-1./lambda1);


  double Psi=PI/2.-atan2(abs(Az),sqrt(Ax*Ax+Ay*Ay));
  double omega=-1.*(PI - atan2(Ay,Ax));
  //double epsilon=0.;


  double rotate[3][3]={{cos(Psi)*cos(omega),cos(Psi)*sin(omega),sin(Psi)},
                       {-1.*sin(omega),cos(omega),0.},
                       {-1.*sin(Psi)*cos(omega),-1.*sin(Psi)*sin(omega),cos(Psi)}};


  double myy_rotate[3];
  double myz_rotate[3];

  TVector3 nominal_pe1_rotate;
  TVector3 nominal_pe2_rotate;

  TVector3 x_iceflowalongx_rotate;
  TVector3 y_iceflowalongx_rotate;

  TVector3 rhat_rotate;
  TVector3 tmpvec;
  TVector3 tmpvec2;
  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum1=0.;
    double sum2=0.;
    double sum3=0.;
    double sum4=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate[i][j]*rhat_iceflowalongx[j];
      sum1+=rotate[i][j]*nominal_pe1[j];
      sum2+=rotate[i][j]*nominal_pe2[j];
      sum3+=rotate[i][j]*x_iceflowalongx[j];
      sum4+=rotate[i][j]*y_iceflowalongx[j];

    }
    rhat_rotate[i]=sum;
    nominal_pe1_rotate[i]=sum1;
    nominal_pe2_rotate[i]=sum2;
    x_iceflowalongx_rotate[i]=sum3;
    y_iceflowalongx_rotate[i]=sum4;

  }

  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate[i][j]*myy_iceflowalongx[j];
      sum2+=rotate[i][j]*myz_iceflowalongx[j];
    }
    myy_rotate[i]=sum;
    myz_rotate[i]=sum2;
  }


  double rotate_T[3][3]={{0.}};
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_T[i][j]=rotate[j][i];
    }
  }

  double epsilon_T=-1.*atan2(rotate_T[1][2],rotate_T[2][2]);
  double Psi_T= asin(rotate_T[0][2]);
  double omega_T=-1.*atan2(rotate_T[0][1],rotate_T[0][0]);

  TVector3 rhat_rotateback;
  TVector3 myy_rotateback;
  TVector3 myz_rotateback;

  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    double sum3=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_T[i][j]*rhat_rotate[j];
      sum2+=rotate_T[i][j]*myy_rotate[j];
      sum3+=rotate_T[i][j]*myz_rotate[j];
    }
    rhat_rotateback[i]=sum;
    myy_rotateback[i]=sum2;
    myz_rotateback[i]=sum3;
  }

  double Anew=1/(a*a)*(cos(Psi_T)*cos(Psi_T)*cos(omega_T)*cos(omega_T)) +
    1/(b*b)*pow(cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T),2) +
    1/(c*c)*pow(sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T),2);

  double Bnew=1/(a*a)*(-2.*cos(Psi_T)*cos(Psi_T)*cos(omega_T)*sin(omega_T)) +
    1/(b*b)*2.*(cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T))*(cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T)) +
    1/(c*c)*2.*(sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T))*(sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T));

  double Cnew=1/(a*a)*(cos(Psi_T)*cos(Psi_T)*sin(omega_T)*sin(omega_T)) +
    1/(b*b)*pow(cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T),2) +
    1/(c*c)*pow(sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T),2);

  /*
  double Dnew=0.;
  double Enew=0.;
  double Fnew=-1.;
  */

  double Mnew=Anew;
  double Pnew=Bnew/2.;
  double Qnew=Bnew/2.;
  double Rnew=Cnew;


  double lambda2_new=(1.*(Mnew+Rnew)+sqrt((Mnew-Rnew)*(Mnew-Rnew)+4*Pnew*Qnew))/2.;
  double lambda1_new=(1.*(Mnew+Rnew)-sqrt((Mnew-Rnew)*(Mnew-Rnew)+4*Pnew*Qnew))/2.;


  double ne1new=sqrt(1./lambda2_new);
  double ne2new=sqrt(1./lambda1_new);


  p_e1[0]=1.;
  p_e1[1]=0.;
  p_e1[2]=0.;

  TVector3 rhat_rotatetheta=rhat_rotate;

  double theta2=atan2( Bnew , Anew - Cnew )/2.; 


  TVector3 findthefreakingaxis(1.,0.,0.);
  TVector3 findthefreakingaxis_perp=findthefreakingaxis;

  findthefreakingaxis.RotateZ(theta2);
  findthefreakingaxis_perp.RotateZ(theta2+PI/2.);


  TVector3 rhat_unrotate=rhat_rotate;

  TVector3 tmpvec3;

  for (int i=0;i<3;i++) {
    double sum1=0.;
    double sum2=0.;
    double sum3=0.;
    for (int j=0;j<3;j++) {
      sum1+=rotate_T[i][j]*findthefreakingaxis[j];
      sum2+=rotate_T[i][j]*findthefreakingaxis_perp[j];
      sum3+=rotate_T[i][j]*rhat_rotate[j];
    }
    tmpvec[i]=sum1;
    tmpvec2[i]=sum2;
    tmpvec3[i]=sum3;
  }
  findthefreakingaxis=tmpvec;
  findthefreakingaxis_perp=tmpvec2;
  rhat_unrotate=tmpvec3;

  TVector3 findthefreakingaxis_projecttoXY(findthefreakingaxis[0],findthefreakingaxis[1],0.);
  TVector3 findthefreakingaxis_perp_projecttoXY(findthefreakingaxis_perp[0],findthefreakingaxis_perp[1],0.);
  TVector3 yaxis(0.,1.,0.);

  double anglebetweenthem=findthefreakingaxis_projecttoXY.Angle(yaxis);



  double diffangle=theta_initial-anglebetweenthem;

  diffangle=0.;
  findthefreakingaxis_projecttoXY.RotateZ(diffangle);
  findthefreakingaxis_perp_projecttoXY.RotateZ(diffangle);
  rhat_unrotate.RotateZ(diffangle);

  findthefreakingaxis[0]=findthefreakingaxis_projecttoXY[0];
  findthefreakingaxis[1]=findthefreakingaxis_projecttoXY[1];

  findthefreakingaxis_perp[0]=findthefreakingaxis_perp_projecttoXY[0];
  findthefreakingaxis_perp[1]=findthefreakingaxis_perp_projecttoXY[1];

  TVector3 rhat_rotateawayfromiceflow;


  double rotate_backtonormal[3][3];

  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
    }

  }

  p_e1=findthefreakingaxis;
  p_e2=findthefreakingaxis_perp;

  // p_e1 is in 1st or 4th quadrant in coordinate system where                                                                                                       
  // ice is along the x axis                                                                                                                                         
  if (!(p_e1.Phi()>-1*PI/2. && p_e1.Phi()<PI/2.))
    p_e1=-1.*p_e1;
  // p_e1 crossed with p_e2 should be in the vertical z direction                                                                                                  
  if ((p_e1.Cross(p_e2)).Dot(myz)<0.)
    p_e2=-1.*p_e2;




  for (int i=0;i<3;i++) {
    double sum3=0.;
    double sum4=0.;
    double sum5=0.;
    for (int j=0;j<3;j++) {
      sum3+=rotate_backtonormal[i][j]*rhat_unrotate[j];
      sum4+=rotate_backtonormal[i][j]*findthefreakingaxis[j];
      sum5+=rotate_backtonormal[i][j]*findthefreakingaxis_perp[j];
    }
    rhat_rotateawayfromiceflow[i]=sum3;
    tmpvec[i]=sum4;
    tmpvec2[i]=sum5;
  }


  findthefreakingaxis=tmpvec;
  findthefreakingaxis_perp=tmpvec2;

  p_e1=findthefreakingaxis;
  p_e2=findthefreakingaxis_perp;


  if (FLIPPED==1) {
    p_e1[2]=-1.*p_e1[2];
    p_e2[2]=-1.*p_e2[2];
    rhat_rotateawayfromiceflow[2]=-1.*rhat_rotateawayfromiceflow[2];
  }

  n_e1=ne1new;
  n_e2=ne2new;

  // ray 1 is the one with the shortest index of refraction.                                                                                                       
  if (n_e2<n_e1) {

    double n_e1_temp=n_e1;
    n_e1=n_e2;
    n_e2=n_e1_temp;

    TVector3 p_e1_temp=p_e1;
    p_e1=p_e2;
    p_e2=p_e1_temp;

  }

  // for an isotropic medium it can tend to pick the orientation of the axes                                                                                       
  // somewhat randomly.                                                                                                                                            
  // this is why when an isotropic medium is chosen, I pick the n1                                                                                                 
  // principal axis to be very slightly longer than the other two.                                                                                                 
  // here is make sure that the 2st eigenvector is the one at 12 o'clock                                                                                           
  if (BIAXIAL==-1) {
    TVector3 temp1=myz.Cross(rhat);
    TVector3 twelveoclock=rhat.Cross(temp1);
    twelveoclock.SetMag(1.);
    double mindotproduct=1000.;
    TVector3 threeoclock=-1.*temp1;
    threeoclock.SetMag(1.);
    mindotproduct=fabs(p_e2.Dot(twelveoclock));
    if (fabs(p_e1.Dot(twelveoclock)>mindotproduct)) {

      double n_e1_temp=n_e1;
      n_e1=n_e2;
      n_e2=n_e1_temp;
      
      TVector3 p_e1_temp=p_e1;
      p_e1=p_e2;
      p_e2=p_e1_temp;
      
    }
    p_e1=threeoclock;
    p_e2=twelveoclock;
    

  }
  
  double deltan=0.;
  if (BIAXIAL==-1)
    deltan=0.;
  else
    deltan=ne2new-ne1new;

  return deltan;
  
}
