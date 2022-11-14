// Find a way to make this inherit from the IceModel class IF being used in AraSim!
// Search for things like <#ifdef> and makefile primers

#include "TRandom3.h"
#include "Constants.h"
#include "Primaries.h"

#include "IceModel.h"
#include "EarthModel.h"
#include "Vector.h"
#include "Ray.h"
#include "Settings.h"
//#include "icemodel.hh"
////#include "earthmodel.hh"
////#include "vector.hh"
//#include "ray.hh"
//
#include "Tools.h"
//
#include <iostream>
#include <fstream>
#include <cstdlib>


// Note: I think we should also add BirefringenceModel.o to the makefile
//  in the same manner as things like IceModel.o
// Make a function to read in the n1, n2, n3 values 
// Would need a .h file that inherits from IceModel (and maybe EarthModel ?)
// Make n1, n2, n3 are private members of the Birefringence class
// The setup of the class looks like this:
//
//
BirefringenceModel::BirefringenceModel() : IceModel(){ // Need arguments for these!

	// Need something like setupBirefringenceModel(model)

}

BirefringenceModel::~BirefringenceModel (){

}

void BirefringenceModel::setUpBirefringenceModel() {

}
// Include functions like getDeltaN from birefringence.cc
// rhat is the usual k vector, AraSim currently has k as the ray direction
//   but we'll need to allow for also having the S vector
//   k is used for propagating the through the ice, S is the used for the detector Beam pattern of the detector
// Both k and S are used in birefringence.cc, but AraSim doesn't currently has S
//
//
// Make a new repository for shared Birefringence files
// Branch off of AraSim to make edits to allow for birefringence
// Add relevant functions to this class
// Make AraSim capable of using both k and S vectors for each ray
//
// Here are the first few functions: 

TVector3 getkfromS(TVector3 S,vector<double> nvec) {

  S.SetMag(1.);

  TVector3 k;

  TVector3 epsilon(nvec[0]*nvec[0],nvec[1]*nvec[1],nvec[2]*nvec[2]);
  TVector3 epsilonS(epsilon[0]*S[0],epsilon[1]*S[1],epsilon[2]*S[2]);
  double magofepsilonS=epsilonS.Mag();

  TVector3 inverseepsilon(1./epsilon[0],1./epsilon[1],1./epsilon[2]);

  k.SetX(inverseepsilon[0]*S[0]);
  k.SetY(inverseepsilon[1]*S[1]);
  k.SetZ(inverseepsilon[2]*S[2]);

  k=magofepsilonS*k;

  return k;

}
TVector3 getSfromk(TVector3 k,vector<double> nvec) {

  k.SetMag(1.);

  TVector3 S;

  TVector3 epsilon(nvec[0]*nvec[0],nvec[1]*nvec[1],nvec[2]*nvec[2]);
  TVector3 epsilonk(epsilon[0]*k[0],epsilon[1]*k[1],epsilon[2]*k[2]);
  
  S=epsilonk;

  S.SetMag(1.);

  return S;

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
