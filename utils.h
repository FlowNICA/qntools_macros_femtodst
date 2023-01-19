TChain* makeChain(string& fileName, const char* treename) {
  cout << "Adding files to chain:" << endl;
  TChain *chain = new TChain(treename);
  if (fileName.rfind(".root") < fileName.size())
    chain->Add(fileName.data());
  else {
    TFileCollection fc("fc", "", fileName.c_str());
    chain->AddFileInfoList((TCollection*)fc.GetList());
  }
  chain->ls();
  return chain;
}

TTree* makeTree4RDF(std::string fileName) {
	StFemtoDstReader* femtoReader = new StFemtoDstReader(fileName.c_str());
  femtoReader->Init();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  femtoReader->SetStatus("*",0);
  femtoReader->SetStatus("Event",1);
  femtoReader->SetStatus("Track",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !femtoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
	Long64_t events2read = femtoReader->chain()->GetEntries();

	float ev_cent, ev_vtxX, ev_vtxY, ev_vtxZ, ev_vpdZ, ev_runId;
	std::vector<int> tr_Nhits, tr_NhitsFit, tr_NhitsPoss, tr_isTofTrack;
	std::vector<double> tr_px, tr_py, tr_pz;
	std::vector<float> tr_dca, tr_chi2, tr_dedx, tr_tofm2, tr_ch;
	std::vector<float> tr_nSigEl, tr_nSigPi, tr_nSigKa, tr_nSigPr;

	TTree *tree = new TTree("tStarData", "Simplified tree with STAR data needed for QnTools");
	tree->Branch("ev_cent", &ev_cent);
	tree->Branch("ev_vtxX", &ev_vtxX);
	tree->Branch("ev_vtxY", &ev_vtxY);
	tree->Branch("ev_vtxZ", &ev_vtxZ);
	tree->Branch("ev_vpdZ", &ev_vpdZ);
	//tree->Branch("ev_runId", &ev_runId);
	tree->Branch("tr_Nhits", &tr_Nhits);
	tree->Branch("tr_NhitsFit", &tr_NhitsFit);
	tree->Branch("tr_NhitsPoss", &tr_NhitsPoss);
	tree->Branch("tr_isTofTrack", &tr_isTofTrack);
	tree->Branch("tr_px", &tr_px);
	tree->Branch("tr_py", &tr_py);
	tree->Branch("tr_pz", &tr_pz);
	tree->Branch("tr_dca", &tr_dca);
	tree->Branch("tr_chi2", &tr_chi2);
	tree->Branch("tr_dedx", &tr_dedx);
	tree->Branch("tr_tofm2", &tr_tofm2);
	tree->Branch("tr_ch", &tr_ch);
	tree->Branch("tr_nSigEl", &tr_nSigEl);
	tree->Branch("tr_nSigPi", &tr_nSigPi);
	tree->Branch("tr_nSigKa", &tr_nSigKa);
	tree->Branch("tr_nSigPr", &tr_nSigPr);

  std::cout << "Number of events to read: " << events2read << std::endl;
	// Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

		if (iEvent % 1000 == 0)
		{
			std::cout << "Preparing event #[" << (iEvent+1)
					      << "/" << events2read << "] for RDataFrame" << std::endl;
		}

    Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      break;
    }

    // Retrieve femtoDst
    StFemtoDst *dst = femtoReader->femtoDst();

    // Retrieve event information
    StFemtoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    // Return primary vertex position
    TVector3 pVtx = event->primaryVertex();

		ev_cent = (15.-(float)event->cent16())*5.+2.5;
		ev_vtxX = pVtx.X();
		ev_vtxY = pVtx.Y();
		ev_vtxZ = pVtx.Z();
		ev_vpdZ = event->vpdVz();
		
		// Track analysis
    Int_t nTracks = dst->numberOfTracks();
		
		tr_Nhits.clear();
		tr_NhitsFit.clear();
		tr_NhitsPoss.clear();
		tr_isTofTrack.clear();
		tr_px.clear();
		tr_py.clear();
		tr_pz.clear();
		tr_dca.clear();
		tr_chi2.clear();
		tr_dedx.clear();
		tr_tofm2.clear();
		tr_ch.clear();
		tr_nSigEl.clear();
		tr_nSigPi.clear();
		tr_nSigKa.clear();
		tr_nSigPr.clear();

    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);

      if (!femtoTrack) continue;

      // Must be a primary track
      if ( !femtoTrack->isPrimary() ) continue;

      TVector3 mom = femtoTrack->gMom();
      TVector3 dca = femtoTrack->gDCA(pVtx);

			tr_px.push_back(mom.X());
			tr_py.push_back(mom.Y());
			tr_pz.push_back(mom.Z());
			tr_dca.push_back(dca.Mag());
			tr_Nhits.push_back(femtoTrack->nHits());
			tr_NhitsFit.push_back(femtoTrack->nHitsFit());
			tr_NhitsPoss.push_back(femtoTrack->nHitsPoss());
			tr_chi2.push_back(femtoTrack->chi2());
			tr_isTofTrack.push_back((int)femtoTrack->isTofTrack());
			tr_dedx.push_back(femtoTrack->dEdx());
			tr_tofm2.push_back(femtoTrack->massSqr());
			tr_ch.push_back(femtoTrack->charge());
			tr_nSigEl.push_back(femtoTrack->nSigmaElectron());
			tr_nSigPi.push_back(femtoTrack->nSigmaPion());
			tr_nSigKa.push_back(femtoTrack->nSigmaKaon());
			tr_nSigPr.push_back(femtoTrack->nSigmaProton());

		} // track loop
		tree->Fill();
	} // event loop
	femtoReader->Finish();
	std::cout << "Chain of femtoDst is prepared to be read by RDataFrame" << std::endl;
	return tree;
}

inline std::function<bool(double)> range(double lo, double hi) 
  {return [=](double x) { return lo <= x && x <= hi; };}

inline std::function<bool(double)> rangeStrict(double lo, double hi) 
  {return [=](double x) { return lo < x && x < hi; };}

inline std::function<bool(double)> equal(int val) 
  {return [=](double _val) { return TMath::Nint(_val) == val; };}

inline bool is (double x) 
  {return fabs(x - 1) < 1e-5;}

inline bool isNot (double x) 
  {return fabs(x) < 1e-5;}

inline bool positive (double x) 
  {return x > 0;}

inline bool negative (double x) 
  {return x < 0;}
