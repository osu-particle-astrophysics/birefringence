
double getDeltaN(int BIAXIAL,vector<double> nvec,TVector3 rhat,double angle_iceflow,
                 double &n_e1, double &n_e2,TVector3 &p_e1,TVector3 &p_e2); // finds indices of refraction of two rays and unit vectors in the direction of the eigenvectors of D
TVector3 rotateD(TVector3 epsilon, double angle_iceflow, TVector3 D); // given D, outputs E, given the epsilon tensor.
void getManyAnglesontheClock(int BIAXIAL,double crosspolangle_tx,
                             TVector3 rhat_thisstep,
                             TVector3 p_e1,TVector3 p_e2,
			     TVector3 E_e1,TVector3 E_e2,
                             double &theta_e1,double &theta_e2, // angles on the clock of the D eigenvectors from the perspective of k 
			     double &thetaE_e1,double &thetaE_e2, // angles on the clock of the E eigenvectors from the perspective of k 
                             double &theta_e1_Sclock,double &theta_e2_Sclock, // angles on the clock of the D eigenvectors from the perspective of S 
                             double &thetaE_e1_Sclock,double &thetaE_e2_Sclock, // angles on the clock of the E eigenvectors from the perspective of S (from these we get the epsilon angles) 
                             TVector3 &Shat_e1,TVector3 &Shat_e2,
                             double &E_e1_thetacomponent, double &E_e2_thetacomponent,
                             double &E_e1_phicomponent, double &E_e2_phicomponent); // this calculates angles on the clock from the perspective of a k vector and an S vector.
void getAnglesontheClock(TVector3 rhat_thisstep1, TVector3 rhat_thisstep2,
                         TVector3 p_e1, TVector3 p_e2,
                         double &theta_e1,double &theta_e2,
                         double &p_e1_thetacomponent,double &p_e2_thetacomponent,
                         double &p_e1_phicomponent,double &p_e2_phicomponent); // this is used many times by getManyAnglesontheClock
void thetastoEpsilons(double thetaE_e1_Sclock,double thetaE_e2_Sclock,
		      double &epsilon1,double &epsilon2); // converts angles on the clock to epsilon angles (relative to 12 o'clock and 3 o'clock)
double getV(vector<double> nvec);

//add my function declarations here!

void birefringenceFileRead(ifstream n1file,ifstream n2file,ifstream n3file, string stemp, double thisn, double thisdepth, double firstn1, double firstn2, double firstn3, int NDEPTHS_NS) //reads in files and gets depth vectors

double smoothVecs(vector<double>  n1vec,vector<double> n2vec,vector<double> n3vec,vector<double>tmp)//smooths vectors 1,2,3

double VAngle(vector<double> nvec_tmp,vector<double> vdepths_n1,vector<double> vdepths_n2,vector<double> vdepths_n3,int n)//finds the angles between the three vectors

void readDepthFiles(ifstream myfile,ifstream Dave5afile,int this_station, int this_day,int this_pol, double this_depth, double this_snrmax){ //take myfile and dave's data and read in values

