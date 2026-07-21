// A file to hold the configuration parameters for cpol_main.cc
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
//
//
struct Config {
	char clswitch;

	double freq=160.E6;             // default frequency [Hz]
	int CROSSPOLANGLE_TX_INT=0;
	int CROSSPOLANGLE_RX_INT=0;
	int BIAXIAL=1;
	int CONSTANTINDICATRIX=0;

	if (argc>1) {
	while ((clswitch = getopt(argc, argv, "b:f:t:u:c:")) != EOF) {
	    switch(clswitch) {
	    case 'b':
		cout << optarg << "\n";
		BIAXIAL=atoi(optarg);
		cout << "biaxial " << BIAXIAL << endl;
		break;
	    case 'f':
		cout << optarg << "\n";
		freq=(double)atof(optarg)*1.E6;
		cout << "freq " << freq << endl;
		break;
	    case 't':
		CROSSPOLANGLE_TX_INT=atoi(optarg);
		cout << "CROSSPOLANGLE_TX_INT " << CROSSPOLANGLE_TX_INT << endl;
		break;
	    case 'u':
		CROSSPOLANGLE_RX_INT=atoi(optarg);
		cout << "CROSSPOLANGLE_RX_INT " << CROSSPOLANGLE_RX_INT << endl;
		break;
	    case 'c':
		CONSTANTINDICATRIX=atoi(optarg);
		cout << "CONSTANTINDICATRIX " << CONSTANTINDICATRIX << endl;
		break;
	    }
	}
	}

	// Convert cross-pol angles to radians
	double CROSSPOLANGLE_TX=(double)CROSSPOLANGLE_TX_INT/DEGRAD;
	double CROSSPOLANGLE_RX=(double)CROSSPOLANGLE_RX_INT/DEGRAD;

	// ------------------------
	// Frequency 
	// ------------------------
	double freqmin=0.;
	double freqmax=1.E9;
	const int NFREQ=100;
	vector<double> vfreqs;


	for (int i=0;i<NFREQ;i++) {
	vfreqs.push_back(freqmin+(freqmax-freqmin)/(double)NFREQ*(double)i);
	}

}
