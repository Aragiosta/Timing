//Digitizer data from the LAB

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};

TH1D* getHistoFromTree(const char *name_file, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h_spectrum->Fill(indata.qlong);
	}
	// return
	return h_spectrum;
}

TH1D* getHistoForChannelFromTree(const char *name_file, short chan, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h_spectrum->Fill(indata.qlong);
	}
	// return
	return h_spectrum;
}


TH1D* getHistoWithFilter(const char *name_file, int numBins, double minX, double maxX, double lowThr = 0, double highThr = 999999) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);
	TH1D *h_spectrum = new TH1D("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		if (indata.qlong>lowThr && indata.qlong<highThr) {
			h_spectrum->Fill(indata.qlong);
		}
	}
	// return
	return h_spectrum;
}


TGraph *getSignal(const char *name_file, int nSamples=250, short chan=0, int nrEv=1){
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TGraph *graph = new TGraph();
	
	//Setting the desired event
	inbranch->GetEntry(nrEv);

	//Looping over the samples
	for (int i=0; i<nSamples; ++i){
		graph->SetPoint(i, i, indata.samples[i]);
	}
	return graph;
}

/////////////////////////////////////////////////////////////////////////////////////
//			ADDED PARTS AS OF 29/10,	~GIOELE


TBranch * getBranchFromTree(const char *name_file, const char *name_tree, const char *name_branch) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get(name_tree);
	TBranch *inbranch = intree->GetBranch(name_branch);

	return inbranch;
}

double min_val=9999999;
double max_val=-1;
//if in doubt on what to use, leave default params
TH1D* createSpectrum(TBranch* inbranch, int numBins = 1e3, double minX = 1, double maxX = -1, int nSamples = 250) {
	// variables
	slimport_data_t indata;
	//TFile *infile = new TFile(name_file);
	//TTree *intree = (TTree*)infile->Get("acq_tree_0");
	//TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);

    TH1D * spectrum = new TH1D("quantized_spectrum","Total spectrum",numBins,minX,maxX);
	//TH1D * spectrum = new TH1D();

    // cycle over all events, and then calculate integral
    for (int event = 0; event < inbranch -> GetEntries(); event++){
        inbranch -> GetEntry(event);

        //integration via a trapezoid, uniformely-discrete grid approximation
        double integral = 0;

        // cycle over quantized samples of waveform and add area of infinitesimal trapezoid
        for(int i = 0; i < nSamples - 1; i++){
            //UShort_t f1,f2;
            //UShort_t dx = 1;    //grid spacing, CHECK IF RIGHT VALUE (1 channel)
			
			double f1,f2;
			double dx = 1;

            // access quantized values
            f1 = indata.samples[i];
            f2 = indata.samples[i+1];

            // calculate integral via trapezoid approx.
            integral += (f1 + f2) * dx / 2.0;
        }

        // check thresholds and populate histogram
        //if (integral >= lowThr && integral <= highThr){
            spectrum -> Fill(integral);
            //if(integral > max_val) max_val = integral;
            //if(integral < min_val) min_val = integral;
        //}
    }

	// return
	return spectrum;
}

