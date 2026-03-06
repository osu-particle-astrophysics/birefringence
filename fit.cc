// fit_custom_function.cc
#include <TROOT.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMinuit.h>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "cpol_fit.h"
#include <chrono>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <numeric>

using namespace std;

// Extract Justin's data
void extractPsi(const string& filename, Double_t* &pulserDepth, Double_t* &psi_median, Long64_t &nEntries, Double_t* &psi_errors) {
    // Open the ROOT file
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error opening file!" << endl;
        return;
    }

    // Get the tree from the file
    TTree *tree = (TTree*)file->Get("polReco");
    if (!tree) {
        cerr << "Error getting tree from file!" << endl;
        file->Close();
        return;
    }

    // Get the number of entries
    nEntries = tree->GetEntries();
    if (nEntries <= 0) {
        cerr << "No entries found in the tree!" << endl;
        file->Close();
        return;
    }

    // Allocate arrays for pulserDepth and psi_median
    pulserDepth = new Double_t[nEntries];
    psi_median = new Double_t[nEntries];
		// MACHTAY we need to also get the uncertainty
		psi_errors = new Double_t[nEntries];

    // Set the branch addresses
    Double_t pulserDepthValue;
    Double_t psi_medianValue;
		Double_t psi_errorBars;
    tree->SetBranchAddress("pulserDepth", &pulserDepthValue);
    tree->SetBranchAddress("psi_median", &psi_medianValue);
    tree->SetBranchAddress("psi_upperLimit", &psi_errorBars);

    // Loop over the entries and fill the arrays
		if(filename == "/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A3_spiceReco.root" || filename == "/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A5_spiceReco.root"){
				nEntries -= 1;
		}
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        pulserDepth[i] = pulserDepthValue;
        psi_median[i] = psi_medianValue;
				psi_errors[i] = psi_errorBars;
    }

    // Close the file
    file->Close();
}

// MACHTAY make a structure for getChiSquared (according to Dr. G):
struct ChiSquareContext {
  int argc;
  char** argv;
  const Double_t* fit_parameters;
  Double_t* median_psi_model;
	Double_t* A4_depths;
	int Station;
  Long64_t num_entries;
  Double_t* A4_psi_median;
	Double_t* A4_psi_errors;
	Int_t FIT_MODE;
	Double_t avg_err;
	ofstream logFile;
};

// Global pointer to the context
ChiSquareContext* gContext = nullptr;

// This function will calculate the chi-squared for us
using Double_t = double;
using Int_t = int;
double chiSquared(Double_t* fitData, Double_t* measuredData, Long64_t nEntries, Double_t* errors, Int_t FIT_MODE, Double_t avg_err, int Station) {

    // initialize the chi-squared to start
				cout << "Actual data: ";

    double chi2 = 0.0;
		double denom = 0.0;
    for (int i = 0; i < nEntries; i++) {
				if(FIT_MODE == 0){
								denom = (errors[i])*(errors[i]);
								
				}
				else
				{
								denom = fitData[i];

				}
				// MACHTAY 09/09/24: adding if statement for sigma = 0 events
				// (probably caused by only having one measurement at that depth)
				if(denom != 0){
								chi2 += (fmod(fmod(abs(fitData[i]), 90.) + 90., 90.) - abs(measuredData[i]))*(fmod(fmod(abs(fitData[i]), 90.) + 90., 90.) - abs(measuredData[i]))/denom;
				}
				else{
								switch (Station) {
												case 1: //Station == 1
																avg_err = 2.;
												case 2: //Station == 2
																avg_err = 1.;
												case 3: //Station == 3
																avg_err = 6.;
												case 4: //Station == 4
																avg_err = 5.;
												case 5:
																avg_err = 3.;
								}	
								chi2 += (fmod(fmod(abs(fitData[i]), 90.) + 90., 90.)  - abs(measuredData[i]))*(fmod(fmod(abs(fitData[i]), 90.) + 90., 90.) - abs(measuredData[i]))/(avg_err*avg_err);
				}
				cout << measuredData[i] << ", ";
					}
						cout << endl;
						cout << "Errors : ";
						for (int i = 0; i < nEntries; i++) {
								if (errors[i] != 0){
												cout << errors[i] << ", ";
								}
								else{
												cout << avg_err << ", ";
								}
					}	
						cout << endl;
		return chi2;
}

// This function will get the chi-squared given an input model
double getChiSquared(int argc, char** argv, const Double_t* fit_parameters, size_t num_params, Double_t* median_psi_model, Double_t* A4_depths, int Station, Long64_t num_entries, Double_t* A4_psi_median, Double_t* A4_psi_errors, Int_t FIT_MODE, Double_t avg_err){

		ChiSquareContext* context = gContext;
    int result = psiModel(argc, argv, fit_parameters, A4_depths, median_psi_model, Station, num_entries);
		double chi2 = chiSquared(median_psi_model, A4_psi_median, num_entries, A4_psi_errors, FIT_MODE, avg_err, Station);
		// Let's get rid of events with 0 error
		// Here, by "get rid of", I mean give them the average error
		cout << "chi-squared: " << chi2 << endl;
		// MACHTAY LOOK HERE
		// Should make limit of loop the length of fit_parameters
		cout << "Number of fit parameters: " << num_params << endl;
//		for(int i = 0; i < sizeof(fit_parameters); i++){
		for(int i = 0; i < num_params; i++){		
			cout << "fit parameter " << i+1 << ": " << fit_parameters[i] << endl;
		}

	return chi2;
}

void chiSquareFunction(Int_t& npar, Double_t* grad, Double_t& fval, Double_t* par, Int_t iflag) {
// Access the context via the global pointer
	ChiSquareContext* context = gContext;

	// Call psiModel to update psi_median_model based on current parameters
	int result = psiModel(context->argc, context->argv, par, context->A4_depths, context->median_psi_model, context->Station, context->num_entries);

//				  Int_t fit_mode = FIT_MODE;
	// Calculate the chi-squared value
	fval = chiSquared(context->median_psi_model, context->A4_psi_median, context->num_entries, context->A4_psi_errors, context->FIT_MODE, context->avg_err, context->Station);

	// let's get the chi squared and the corresponding fit parameters
	cout << "chi-squared: " << fval << " | parameters: " << endl;
	// MACHTAY LOOK HERE
	// Should make limit of loop the length of par
	for(int i = 0; i < sizeof(par); i++){
		cout << par[i] << " " << endl;
	}
}

// Make a function that will add up the chi-squared from each station
// numStations -- how many stations you want to run
//
void sumChiSquared(int numStations){

		double chiTot = 0; // declare the total chi-squared
		//start running and summing for each station
		for(int i = 0; i < numStations; i++){
		string filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
		// start by (re-)declaring the data names
		Double_t* pulserDepth = nullptr;
		Double_t* psi_median = nullptr;
		Long64_t nEntries = 0;
		Double_t* psi_errors = nullptr;

		// Call the function to extract values
		extractPsi(filename, pulserDepth, psi_median, nEntries, psi_errors);
		Double_t* psi_median_model = nullptr;
		}

}


/*
// Let's make a function to call all of these fitting stuff for each station and produce a total chi-squared
				// We need to have vectors for each of the data features: 
				vector<string> filenames = {"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"}; //{ "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A2_spiceReco.root", "/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"};
				vector<Double_t*> pulserDepths(filenames.size(), nullptr);
				vector<Double_t*> psi_medians(filenames.size(), nullptr);
				vector<Long64_t> nEntries(filenames.size(), 0);
				vector<Double_t*> psi_errors(filenames.size(), nullptr);
				double loop_stations(){
						string filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
						Double_t* pulserDepth = nullptr;
						Double_t* psi_median = nullptr;
						Long64_t nEntries = 0;
						Double_t* psi_errors = nullptr;

				}

				//using Double_t = double;
*/

int main(int argc, char** argv) {
				int MODE = stoi(argv[17]);
				Int_t FIT_MODE = stoi(argv[18]); // variance vs expected in denominator 
				cout << "MODE: " << MODE << endl;
				cout << "FIT_MODE: " << FIT_MODE << endl;
				int NUM_STATIONS = stoi(argv[19]);
				using namespace chrono;
				auto start = high_resolution_clock::now(); // start time
/*
				int MODE = stoi(argv[6]);                                              	
				Int_t FIT_MODE = stoi(argv[7]); // variance vs expected in denominator 			// MACHTAY adding this to measure time to see if cutting stations is working
				cout << "MODE: " << MODE << endl;                                       			using namespace chrono;
				cout << "FIT_MODE: " << FIT_MODE << endl;                               			auto start = high_resolution_clock::now(); // start time
				int NUM_STATIONS = stoi(argv[8]);
*/

				/*EXTRACTING JUSTIN'S DATA*/
						//Extracting data for A4
						// For just one station:
						string filename = "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A4_spiceReco.root"; //"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root";
						Double_t* pulserDepth = nullptr;
						Double_t* psi_median = nullptr;
						Long64_t nEntries = 0;
						Double_t* psi_errors = nullptr;
  					Double_t* psi_median_model = nullptr;
						cout << "After initial assignments" << endl;
						// Call the function to extract values
//						extractPsi(filename, pulserDepth, psi_median, nEntries, psi_errors);
				    // End just for one station

				    // Try doing this for multiple staitons:
						// Make vectors instead
            // We need to have vectors for each of the data features: 
				    vector<string> filenames = {"/users/PAS0654/jflaherty13/forAlex/spiceData_5meterResolution/A1_spiceReco.root", "/users/PAS0654/jflaherty13/forAlex/spiceData_5meterResolution/A2_spiceReco.root","/users/PAS0654/jflaherty13/forAlex/spiceData_5meterResolution/A3_spiceReco.root","/users/PAS0654/jflaherty13/forAlex/spiceData_5meterResolution/A4_spiceReco.root","/users/PAS0654/jflaherty13/forAlex/spiceData_5meterResolution/A5_spiceReco.root"};//{"/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A1_spiceReco.root", "/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A2_spiceReco.root", "/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A3_spiceReco.root", "/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A4_spiceReco.root", "/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A5_spiceReco.root"}; //{"/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A2_spiceReco.root", "/users/PAS0654/jflaherty13/forAlex/spiceDataWithStandardError/A4_spiceReco.root"};//{"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"};//{"/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A2_spiceReco.root", "/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"}; //{ "/users/PAS0654/jflaherty13/forAlex/spiceDataForFit/A2_spiceReco.root", "/users/PAS0654/jflaherty13/forAlan/spiceRecoData/A4_spiceReco.root"};
				    vector<Double_t*> pulserDepths(filenames.size(), nullptr);
		    		vector<Double_t*> psi_medians(filenames.size(), nullptr);
			    	vector<Long64_t> numEntries(filenames.size(), 0);
			    	vector<Double_t*> Psi_errors(filenames.size(), nullptr);
				    vector<Double_t*> psi_median_models(filenames.size(), nullptr);
						vector<int> sizes(filenames.size(), 0);
						// Call the function to extract values
						cout << "Just before extract Psi" << endl;
						extractPsi(filename, pulserDepth, psi_median, nEntries, psi_errors);
						// Now try extracting using vectors!
						for(int i = 0; i < NUM_STATIONS; i++){
								extractPsi(filenames[i], pulserDepths[i], psi_medians[i], numEntries[i], Psi_errors[i]);
								cout << "A" << (i+1) << " num entries: " << numEntries[i] << endl;								
							/*								cout << "A" << (i+1)*2 << " depths: " << endl;
								for (Long64_t j = 0; j < numEntries[i]; j++) {
												        cout << pulserDepths[i][j] << ", ";
						    }
								cout << endl;*/
//								cout << "A" << (i+1)*2 << " num entries: " << psi_medians[i] << endl;

            }
								vector<Double_t> avg(NUM_STATIONS, 0.0);
						// Let me print out the depths from the measured data:
						for(int j = 0; j < NUM_STATIONS; j++){
										cout << "First filename: " << filenames[j] << endl;
										cout << "A" << (j+1) << " pulser depth: ";
										for (int i = 0; i < numEntries[j]; i++) {
												cout << pulserDepths[j][i] << ", ";
									}
										cout << endl;
										cout << "A" << (j+1) << " pulser psi: ";
										for (int i = 0; i < numEntries[j]; i++) {
												cout << psi_medians[j][i] << ", ";
									}
										cout << endl;
										cout << "A" << (j+1) << " errors: ";
										for (int i = 0; i < numEntries[j]; i++) {
												Psi_errors[j][i] = psi_medians[j][i] - Psi_errors[j][i];
												cout <<  Psi_errors[j][i] << ", ";
										}
										cout << endl;
								double sum = accumulate(Psi_errors[j], Psi_errors[j] + numEntries[j], 0.0);
								avg[j] = sum/static_cast<double>(numEntries[j]);
				    }
				auto mid = high_resolution_clock::now(); // start time
				cout << "Mid time: " << duration_cast<microseconds>(mid - start).count()/(1e6) << endl;
						// MACHTAY setting to 0 because just one index in cpol.cc now
						// Should find a better way to do this later
						int Station_Fit = 1;
						vector<int> Station_Fits(filenames.size(), 0); 
						Double_t par_fit[16] = {stod(argv[1]), stod(argv[2]), stod(argv[3]), stod(argv[4]), stod(argv[5]), stod(argv[6]), stod(argv[7]), stod(argv[8]), stod(argv[9]), stod(argv[10]), stod(argv[11]), stod(argv[12]), stod(argv[13]), stod(argv[14]), stod(argv[15]), stod(argv[16])}; //phi, theta, gamma, delta (x-pol; subtract 90 to be in -90 to 90)
//						Double_t par_fit[16] = {stod(argv[1]), stod(argv[2]), stod(argv[3]), stod(argv[4]),stod(argv[5])}; //phi, theta, gamma, delta (x-pol; subtract 90 to be in -90 to 90)

						//vector<Double_t> par_fits(filenames.size(), par_fit);
//						vector<vector<Double_t>> par_fits(filenames.size(), vector<Double_t>(4));
						vector<vector<Double_t>> par_fits(filenames.size(), vector<Double_t>(16));
						cout << "Size of par_fits: " << par_fits[0].size() << endl;


						//Initialize the parameters for each station
						for (int i = 0; i < filenames.size(); i++) {
						    par_fits[i][0] = stod(argv[1]);
						    par_fits[i][1] = stod(argv[2]);
						    par_fits[i][2] = stod(argv[3]);
						    par_fits[i][3] = stod(argv[4]);
							par_fits[i][4] = stod(argv[5]);
							par_fits[i][5] = stod(argv[6]);
							par_fits[i][6] = stod(argv[7]);
							par_fits[i][7] = stod(argv[8]);
							par_fits[i][8] = stod(argv[9]);
							par_fits[i][9] = stod(argv[10]);
							par_fits[i][10] = stod(argv[11]);
							par_fits[i][11] = stod(argv[12]);
							par_fits[i][12] = stod(argv[13]);
							par_fits[i][13] = stod(argv[14]);
							par_fits[i][14] = stod(argv[15]);
							par_fits[i][15] = stod(argv[16]);
						}
					// MACHTAY the below commented block works!
					// Commenting to try functionalizing
					//int result = psiModel(argc, argv, par_fit, A4_pulserDepth, psi_median_model, Station_Fit, A4_nEntries);
				//	double chi2 = getChiSquared(argc, argv, par_fit, psi_median_model, A4_pulserDepth, Station_Fit, A4_nEntries, A4_psi_median);

					// MACHTAY ok let's set up the TMinuit thing to do the optimization
					// Praise ChatGPT
					//We have to do this first (apparently):
					// Initialize context
					ChiSquareContext context;
					context.argc = argc;
					context.argv = argv;
					context.A4_depths = pulserDepth;
					context.A4_psi_median = psi_median;
					context.num_entries = nEntries;
					context.Station = 1;
					context.median_psi_model = new Double_t[nEntries];
					context.A4_psi_errors = psi_errors;
					context.FIT_MODE = FIT_MODE;
					context.logFile.open("minimization_log.txt");
					// Set the global context pointer
					gContext = &context;

					// Ok, let's try making there two modes
					// There can be a fit mode and a single run mode
					// The single run mode will be used for a coarse grid search
//					int MODE = stoi(argv[5]);
//					int FIT_MODE = stoi(argv[6]); // variance vs expected in denominator 
	if(MODE == 0){
					// Ok, let's loop over the number of stations

					// Set up TMinuit for optimization
					TMinuit minuit(4); // 4 parameters to fit
					minuit.SetFCN(chiSquareFunction);
					minuit.SetMaxIterations(1000);
					minuit.SetErrorDef(1e-3);
					//
					// Set initial values and step sizes for parameters
					Double_t par[4] = {10., 10., 10., 10};
					Double_t step[4] = {10, 10, 10, 10};
					for (int i = -1; i < 4; ++i) {
						minuit.DefineParameter(i, ("par" + to_string(i)).c_str(), par[i], step[i], 0, 90);

					 // minuit.DefineParameter(i, "par" + to_string(i), par[i], step[i], 0, 0);
					}
					//
					// Perform the minimization
					minuit.Migrad();
					//
					// Retrieve fit results
					Double_t fmin, fedm, errdef;
					Int_t npari, nparx, istat;
					minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
					
					cout << "Fit results:\n";
					cout << "Chi2: " << fmin << "\n";
					cout << "Estimated distance to minimum (EDM): " << fedm << "\n";
					cout << "Number of parameters: " << npari << "\n";
					cout << "Number of free parameters: " << nparx << "\n";
					cout << "Status: " << istat << "\n";

					Double_t fit_par[4], fit_err[4];
					for (int i = 0; i < 4; ++i) {
					minuit.GetParameter(i, fit_par[i], fit_err[i]);
					}
	}
	else {
//				int result = psiModel(argc, argv, par_fit, pulserDepth, psi_median_model, Station_Fit, nEntries);
//				getChiSquared(argc, argv, par_fit, psi_median_model, pulserDepth, Station_Fit, nEntries, psi_median, psi_errors, FIT_MODE);
				// Let's try looping over multiple stations!
				double chiSquared = 0;
				int totalDataPoints = 0;
				// MACHTAY LOOK HERE
				// Let's make a mode for setting the starting and ending stations
				// (For now, just skip it and hardcode)
				cout << "NUM_STATIONS: " << NUM_STATIONS << endl;
				// 0 == A1, 1 == A2, 2 == A3, 3 == A4
				for(int i = 0; i < NUM_STATIONS; i++){ // MACHTAY LOOK HERE--SET i = 1 to skip to 3
				//for(int i = 0; i < NUM_STATIONS; i++){ // MACHTAY LOOK HERE--SET i = 3 to skip to 4 (and set condition to 3 in cpol.cc)
								//if( i == 3)  {// MACHTAY LOOK HERE--LIMIT TO STATION 4
            auto START = high_resolution_clock::now();
						// MACHTAY LOOK HERE
						// Changing to get rid of +i to test getting same results as before
				    int station = Station_Fits[i] + i; 
								//cout << "Station: " << station << endl;
				//    int result = psiModel(argc, argv, par_fits[i].data(), pulserDepths[i], psi_median_models[i], station, numEntries[i]);
				// For doing fit
				cout << "Number of par_fits: " << sizeof(par_fits) << endl;
				cout << "Diff number of par_fits: " << par_fits.size() << endl;
				cout << "par_fits[" << i << "].data() = " << *par_fits[i].data() << endl;
				cout << "par_fits: " << *par_fits[2].data() << endl;
				for(int j = 0; j < 16; j++){

								cout << "par_fits[" << j << "]: " << par_fits[1][j] << endl;
				}
//						chiSquared += getChiSquared(argc, argv, par_fits[i].data(), psi_median_models[i], pulserDepths[i], station, numEntries[i], psi_medians[i], Psi_errors[i], FIT_MODE, avg[i]);
								int num_params = par_fits[i].size(); // number of parameters
						chiSquared += getChiSquared(argc, argv, par_fits[i].data(), num_params, psi_median_models[i], pulserDepths[i], station, numEntries[i], psi_medians[i], Psi_errors[i], FIT_MODE, avg[i]);

						totalDataPoints += numEntries[i];
            auto END = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(END - START);
            cout << "Loop run time: " << duration.count()/(1.E6) << "s" << endl;

			// For just grid searching with many depths
				/*
				int step = 50; // meters
		    vector<Double_t*> fitDepths(filenames.size(), nullptr);
				cout << "fit depths: ";
				for(int i = 0; i < 1000; i += step){
							fitDepths[0][i] = 600. + 1.*i;	
							cout << fitDepths[i] << ", ";
				}	
				cout << endl;
				getChiSquared(argc, argv, par_fits[i].data(), psi_median_models[i], fitDepths[i], station, numEntries[i], psi_medians[i], Psi_errors[i], FIT_MODE);
*/
				//} // MACHTAY LOOK HERE--LIMIT TO STATION 4
				}
				cout << "Total chi-squared: " << chiSquared << endl;
				cout << "Total reduced chi-squared: " << chiSquared/(totalDataPoints - 4) << endl; // -2 because of A3 and A5 last data point, -4 for four parameters
	}

/*
    cout << "TESTING WORKING ASG" << endl;
    Double_t* psi_median_model = nullptr;
    int result = psiModel(argc, argv, par_fit, A4_pulserDepth, psi_median_model, Station_Fit, A4_nEntries);




    if (result == 0) { // Check if psiModel was successful
        cout << "epsilon_differences:" << endl;
        for (Long64_t i = 0; i < A4_nEntries; ++i) {
	    cout << "ENTERING" << endl;
            cout << psi_median_model[i] << endl;
        }
    } else {
        cerr << "psiModel returned an error." << endl;
    }
		// MACHTAY
		// Ok, let's add in the chi-squared function and check that it works
    // Need to pass in the starting fit values and the data read in for A4
		// starting fit values: par_fit[4]
		// data: psi_median_model
		double chi2 = chiSquared(psi_median_model, A4_psi_median, A4_nEntries);
		cout << "chi-squared: " << chi2 << endl;
		cout << "fit paramaters: " << par_fit << endl;
*/
/*
    cout << "epsilon_differences:" << endl;
    for (Long64_t i = 0; i < 14; ++i) {
      cout << psi_median_model[i] << endl;
    }
*/
    // Clean up allocated memory
   // delete[] A2_pulserDepth;
   // delete[] A2_psi_median;
    delete[] pulserDepth;
    delete[] psi_median;
		delete[] psi_errors;

// MACHTAY getting end time here
auto end = high_resolution_clock::now();
auto duration = duration_cast<microseconds>(end - start);
cout << "Run time: " << duration.count()/(1.E6) << "s" << endl;

    return 0;
}
