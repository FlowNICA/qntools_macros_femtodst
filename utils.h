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

//--------------------------------------------------------------------------------------------------------------------------//
//this function simply connects the gain values read in to the BBC azimuthal distribution
//since tiles 7 and 9 (+ 13 and 15) share a gain value it is ambiguous how to assign the geometry here
//I prefer assigning the angle between the tiles thus "greying out" the adcs. 
//Others have assigned all of the adc to one (exclusive) or the the other. 
Float_t BBC_GetPhi(const Int_t eastWest, const Int_t tileId)
{
  //float GetPhiInBBC(int eastWest, int bbcN) { //tileId=0 to 23
  const float Pi = TMath::Pi();
  const float Phi_div = Pi / 6;
  float bbc_phi = Phi_div;
  switch(tileId) {
    case 0: bbc_phi = 3.*Phi_div;
  break;
    case 1: bbc_phi = Phi_div;
  break;
    case 2: bbc_phi = -1.*Phi_div;
  break;
    case 3: bbc_phi = -3.*Phi_div;
  break;
    case 4: bbc_phi = -5.*Phi_div;
  break;
    case 5: bbc_phi = 5.*Phi_div;
  break;
    //case 6: bbc_phi= (mRndm.Rndm() > 0.5) ? 2.*Phi_div:4.*Phi_div;	//tiles 7 and 9 are gained together we randomly assign the gain to one XOR the other
    case 6: bbc_phi = 3.*Phi_div;
  break;
    case 7: bbc_phi = 3.*Phi_div;
  break;
    case 8: bbc_phi = Phi_div;
  break;
    case 9: bbc_phi = 0.;
  break;
    case 10: bbc_phi = -1.*Phi_div;
  break;
    //case 11: bbc_phi = (mRndm.Rndm() > 0.5) ? -2.*Phi_div:-4.*Phi_div;	//tiles 13 and 15 are gained together
    case 11: bbc_phi = -3.*Phi_div;
  break;
    case 12: bbc_phi = -3.*Phi_div;
  break;
    case 13: bbc_phi = -5.*Phi_div;
  break;
    case 14: bbc_phi = Pi;
  break;
    case 15: bbc_phi = 5.*Phi_div;
  break;
  }

  //if we're looking at the east BBC we need to flip around x in the STAR coordinates, 
  //a line parallel to the beam would go through tile 1 on the W BBC and tile 3 on the 
  if(0 == eastWest){
    if (bbc_phi > -0.001){ //this is not a >= since we are talking about finite adcs -- not to important
      bbc_phi = Pi - bbc_phi;
    }
    else {
      bbc_phi= -Pi - bbc_phi;
    }
  }

  if(bbc_phi < 0.0) bbc_phi += 2.*Pi;
  if(bbc_phi > 2.*Pi) bbc_phi -= 2.*Pi;

  return bbc_phi;
}

Double_t GetZDCPosition(Int_t eastwest, Int_t verthori, Int_t strip)
// Get position of each slat;strip starts from 0
{

  std::vector<Double_t> zdcsmd_x = {0.5,2,3.5,5,6.5,8,9.5};
  std::vector<Double_t> zdcsmd_y = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  Double_t mZDCSMDCenterex = 4.72466;
  Double_t mZDCSMDCenterey = 5.53629;
  Double_t mZDCSMDCenterwx = 4.39604;
  Double_t mZDCSMDCenterwy = 5.19968;

  if(eastwest==0 && verthori==0) return zdcsmd_x.at(strip)-mZDCSMDCenterex;
  if(eastwest==1 && verthori==0) return mZDCSMDCenterwx-zdcsmd_x.at(strip);
  if(eastwest==0 && verthori==1) return zdcsmd_y.at(strip)/sqrt(2.)-mZDCSMDCenterey;
  if(eastwest==1 && verthori==1) return zdcsmd_y.at(strip)/sqrt(2.)-mZDCSMDCenterwy;

  return -999.;
}

Double_t GetZDCPhi(Int_t eastwest, Int_t verthori, Int_t strip)
{
  double position = GetZDCPosition(eastwest, verthori, strip);
  if (position == -999.) return -999.;
  
  double cos, sin, phi;
  cos = -999.;
  sin = -999.;
  if(eastwest==0 && verthori==0)
  {
    cos = position;
    sin = 0.;
  }
  if(eastwest==1 && verthori==0)
  {
    cos = position;
    sin = 0.;
  }
  if(eastwest==0 && verthori==1)
  {
    cos = 0.;
    sin = position;
  }
  if(eastwest==1 && verthori==1)
  {
    cos = 0.;
    sin = position;
  }

  if (cos == -999.) return -999.;
  if (sin == -999.) return -999.;

  phi = atan2(sin,cos);

  return phi;
}

pair<double, double> GetXY(double _val_m2, double _val_sig, 
  double _fit_mean_m2_pi, double _fit_mean_m2_ka,
  double _fit_mean_sig_pi, double _fit_mean_sig_ka,
  double _fit_sig_m2_pi,
  double _fit_sig_sig_pi)
{
  double f = _fit_sig_sig_pi / _fit_sig_m2_pi;
  double alpha = -1.*atan2( (_fit_mean_m2_ka - _fit_mean_m2_pi), ((_fit_mean_sig_ka - _fit_mean_sig_pi)/f) );

  double x0 = (_val_sig - _fit_mean_sig_pi)/f;
  double y0 = _val_m2 - _fit_mean_m2_pi;

  double x = cos(alpha)*x0 - sin(alpha)*y0;
  double y = sin(alpha)*x0 + cos(alpha)*y0;
  return {x,y};
}

TTree* makeTree4RDF(std::string fileName="",
                    TGraphErrors *const& gr_mean_m2_pi=nullptr, TGraphErrors *const& gr_mean_ns_pi=nullptr, TGraphErrors *const& gr_sigm_m2_pi=nullptr, TGraphErrors *const& gr_sigm_ns_pi=nullptr,
                    TGraphErrors *const& gr_mean_m2_ka=nullptr, TGraphErrors *const& gr_mean_ns_ka=nullptr, TGraphErrors *const& gr_sigm_m2_ka=nullptr, TGraphErrors *const& gr_sigm_ns_ka=nullptr,
                    TGraphErrors *const& gr_mean_m2_pr=nullptr, TGraphErrors *const& gr_mean_ns_pr=nullptr, TGraphErrors *const& gr_sigm_m2_pr=nullptr, TGraphErrors *const& gr_sigm_ns_pr=nullptr,
                    TGraphErrors *const& gr_mean_x_pi=nullptr,  TGraphErrors *const& gr_mean_y_pi=nullptr,  TGraphErrors *const& gr_sigm_x_pi=nullptr,  TGraphErrors *const& gr_sigm_y_pi=nullptr,
                    TGraphErrors *const& gr_mean_x_ka=nullptr,  TGraphErrors *const& gr_mean_y_ka=nullptr,  TGraphErrors *const& gr_sigm_x_ka=nullptr,  TGraphErrors *const& gr_sigm_y_ka=nullptr,
                    TGraphErrors *const& gr_mean_x_pr=nullptr,  TGraphErrors *const& gr_mean_y_pr=nullptr,  TGraphErrors *const& gr_sigm_x_pr=nullptr,  TGraphErrors *const& gr_sigm_y_pr=nullptr)
{
  bool isPidXY = true;
  if (!gr_mean_m2_pi) { isPidXY = false; std::cerr << "Failed to open gr_mean_m2_pi!" << std::endl; }
  if (!gr_mean_ns_pi) { isPidXY = false; std::cerr << "Failed to open gr_mean_ns_pi!" << std::endl; }
  if (!gr_sigm_m2_pi) { isPidXY = false; std::cerr << "Failed to open gr_sigm_m2_pi!" << std::endl; }
  if (!gr_sigm_ns_pi) { isPidXY = false; std::cerr << "Failed to open gr_sigm_ns_pi!" << std::endl; }
  if (!gr_mean_m2_ka) { isPidXY = false; std::cerr << "Failed to open gr_mean_m2_ka!" << std::endl; }
  if (!gr_mean_ns_ka) { isPidXY = false; std::cerr << "Failed to open gr_mean_ns_ka!" << std::endl; }
  if (!gr_sigm_m2_ka) { isPidXY = false; std::cerr << "Failed to open gr_sigm_m2_ka!" << std::endl; }
  if (!gr_sigm_ns_ka) { isPidXY = false; std::cerr << "Failed to open gr_sigm_ns_ka!" << std::endl; }
  if (!gr_mean_m2_pr) { isPidXY = false; std::cerr << "Failed to open gr_mean_m2_pr!" << std::endl; }
  if (!gr_mean_ns_pr) { isPidXY = false; std::cerr << "Failed to open gr_mean_ns_pr!" << std::endl; }
  if (!gr_sigm_m2_pr) { isPidXY = false; std::cerr << "Failed to open gr_sigm_m2_pr!" << std::endl; }
  if (!gr_sigm_ns_pr) { isPidXY = false; std::cerr << "Failed to open gr_sigm_ns_pr!" << std::endl; }
  if (!gr_mean_x_pi)  { isPidXY = false; std::cerr << "Failed to open gr_mean_x_pi!" << std::endl; }
  if (!gr_mean_y_pi)  { isPidXY = false; std::cerr << "Failed to open gr_mean_y_pi!" << std::endl; }
  if (!gr_sigm_x_pi)  { isPidXY = false; std::cerr << "Failed to open gr_sigm_x_pi!" << std::endl; }
  if (!gr_sigm_y_pi)  { isPidXY = false; std::cerr << "Failed to open gr_sigm_y_pi!" << std::endl; }
  if (!gr_mean_x_ka)  { isPidXY = false; std::cerr << "Failed to open gr_mean_x_ka!" << std::endl; }
  if (!gr_mean_y_ka)  { isPidXY = false; std::cerr << "Failed to open gr_mean_y_ka!" << std::endl; }
  if (!gr_sigm_x_ka)  { isPidXY = false; std::cerr << "Failed to open gr_sigm_x_ka!" << std::endl; }
  if (!gr_sigm_y_ka)  { isPidXY = false; std::cerr << "Failed to open gr_sigm_y_ka!" << std::endl; }
  if (!gr_mean_x_pr)  { isPidXY = false; std::cerr << "Failed to open gr_mean_x_pr!" << std::endl; }
  if (!gr_mean_y_pr)  { isPidXY = false; std::cerr << "Failed to open gr_mean_y_pr!" << std::endl; }
  if (!gr_sigm_x_pr)  { isPidXY = false; std::cerr << "Failed to open gr_sigm_x_pr!" << std::endl; }
  if (!gr_sigm_y_pr)  { isPidXY = false; std::cerr << "Failed to open gr_sigm_y_pr!" << std::endl; }

  const vector<double> pid_Msqr = {0.13957*0.13957, 0.493677*0.493677, 0.938272*0.938272}; //pi, K, p m^2 in (GeV/c^2)^2

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

  int ev_tofMatched;
	float ev_cent, ev_vtxX, ev_vtxY, ev_vtxZ, ev_vpdZ, ev_runId;
	std::vector<int> tr_Nhits, tr_NhitsFit, tr_NhitsPoss, tr_isTofTrack;
	std::vector<double> tr_px, tr_py, tr_pz;
	std::vector<float> tr_dca, tr_chi2, tr_dedx, tr_tofm2, tr_ch;
	std::vector<float> tr_nSigEl, tr_nSigPi, tr_nSigKa, tr_nSigPr;
	std::vector<int> bbc_id, bbc_side; // id, side (west, east)
	std::vector<float> bbc_en, bbc_phi;
	std::vector<int> zdc_id, zdc_side, zdc_type; // id, side (west, east), type (vertical, horizontal)
	std::vector<float> zdc_en, zdc_phi;
  std::vector<float> tr_pidX, tr_pidY; // PID related information
  std::vector<float> tr_pidMeanXPi, tr_pidSigmXPi, tr_pidMeanYPi, tr_pidSigmYPi; // PID related information for pions
  std::vector<float> tr_pidMeanXKa, tr_pidSigmXKa, tr_pidMeanYKa, tr_pidSigmYKa; // PID related information for kaons
  std::vector<float> tr_pidMeanXPr, tr_pidSigmXPr, tr_pidMeanYPr, tr_pidSigmYPr; // PID related information for protons

  const int fNBbcModules = 16;
  const int fNZdcSmdStripsHorizontal = 8;
  const int fNZdcSmdStripsVertical   = 7;

	TTree *tree = new TTree("tStarData", "Simplified tree with STAR data needed for QnTools");
	tree->Branch("ev_cent", &ev_cent);
	tree->Branch("ev_vtxX", &ev_vtxX);
	tree->Branch("ev_vtxY", &ev_vtxY);
	tree->Branch("ev_vtxZ", &ev_vtxZ);
	tree->Branch("ev_vpdZ", &ev_vpdZ);
	tree->Branch("ev_tofMatched", &ev_tofMatched);
	// tree->Branch("ev_runId", &ev_runId);
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
	tree->Branch("tr_pidX", &tr_pidX);
	tree->Branch("tr_pidY", &tr_pidY);
	tree->Branch("tr_pidMeanXPi", &tr_pidMeanXPi);
	tree->Branch("tr_pidSigmXPi", &tr_pidSigmXPi);
  tree->Branch("tr_pidMeanYPi", &tr_pidMeanYPi);
	tree->Branch("tr_pidSigmYPi", &tr_pidSigmYPi);
  tree->Branch("tr_pidMeanXKa", &tr_pidMeanXKa);
	tree->Branch("tr_pidSigmXKa", &tr_pidSigmXKa);
  tree->Branch("tr_pidMeanYKa", &tr_pidMeanYKa);
	tree->Branch("tr_pidSigmYKa", &tr_pidSigmYKa);
  tree->Branch("tr_pidMeanXPr", &tr_pidMeanXPr);
	tree->Branch("tr_pidSigmXPr", &tr_pidSigmXPr);
  tree->Branch("tr_pidMeanYPr", &tr_pidMeanYPr);
	tree->Branch("tr_pidSigmYPr", &tr_pidSigmYPr);
	tree->Branch("bbc_id", &bbc_id);
	tree->Branch("bbc_side", &bbc_side);
	tree->Branch("bbc_en", &bbc_en);
	tree->Branch("bbc_phi", &bbc_phi);
  tree->Branch("zdc_id", &zdc_id);
	tree->Branch("zdc_side", &zdc_side);
	tree->Branch("zdc_type", &zdc_type);
	tree->Branch("zdc_en", &zdc_en);
  tree->Branch("zdc_phi", &zdc_phi);
  
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
    tr_pidX.clear();
    tr_pidY.clear();
    tr_pidMeanXPi.clear();
    tr_pidSigmXPi.clear();
    tr_pidMeanYPi.clear();
    tr_pidSigmYPi.clear();
    tr_pidMeanXKa.clear();
    tr_pidSigmXKa.clear();
    tr_pidMeanYKa.clear();
    tr_pidSigmYKa.clear();
    tr_pidMeanXPr.clear();
    tr_pidSigmXPr.clear();
    tr_pidMeanYPr.clear();
    tr_pidSigmYPr.clear();

    ev_tofMatched = 0;

    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);

      if (!femtoTrack) continue;

      // Must be a primary track
      if ( !femtoTrack->isPrimary() ) continue;
      ev_tofMatched++;

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
			tr_isTofTrack.push_back((int) (femtoTrack->isTofTrack()) ? 1 : 0);
			tr_dedx.push_back(femtoTrack->dEdx());
			tr_tofm2.push_back(femtoTrack->massSqr());
			tr_ch.push_back(femtoTrack->charge());
			tr_nSigEl.push_back(femtoTrack->nSigmaElectron());
			tr_nSigPi.push_back(femtoTrack->nSigmaPion());
			tr_nSigKa.push_back(femtoTrack->nSigmaKaon());
			tr_nSigPr.push_back(femtoTrack->nSigmaProton());
      auto pt = (float) mom.Pt();

      //PID related info
      if (isPidXY && pt > 0.3 && pt < 3.2) {
        double mean_ns_pi = gr_mean_ns_pi->Eval(pt,nullptr,"S");
        double mean_ns_ka = gr_mean_ns_ka->Eval(pt,nullptr,"S");
        double sigm_m2_pi = gr_sigm_m2_pi->Eval(pt,nullptr,"S");
        double sigm_ns_pi = gr_sigm_ns_pi->Eval(pt,nullptr,"S");
        double mean_x_pi  = gr_mean_x_pi->Eval(pt,nullptr,"S");
        double mean_y_pi  = gr_mean_y_pi->Eval(pt,nullptr,"S");
        double sigm_x_pi  = gr_sigm_x_pi->Eval(pt,nullptr,"S");
        double sigm_y_pi  = gr_sigm_y_pi->Eval(pt,nullptr,"S");
        double mean_x_ka  = gr_mean_x_ka->Eval(pt,nullptr,"S");
        double mean_y_ka  = gr_mean_y_ka->Eval(pt,nullptr,"S");
        double sigm_x_ka  = gr_sigm_x_ka->Eval(pt,nullptr,"S");
        double sigm_y_ka  = gr_sigm_y_ka->Eval(pt,nullptr,"S");
        double mean_x_pr  = gr_mean_x_pr->Eval(pt,nullptr,"S");
        double mean_y_pr  = gr_mean_y_pr->Eval(pt,nullptr,"S");
        double sigm_x_pr  = gr_sigm_x_pr->Eval(pt,nullptr,"S");
        double sigm_y_pr  = gr_sigm_y_pr->Eval(pt,nullptr,"S");

        auto xy_coord = GetXY(tr_tofm2.at(iTrk), tr_nSigPi.at(iTrk),
                              pid_Msqr.at(0), pid_Msqr.at(1),
                              mean_ns_pi, mean_ns_ka, sigm_m2_pi, sigm_ns_pi);
        tr_pidX.push_back(xy_coord.first);
        tr_pidY.push_back(xy_coord.second);
        tr_pidMeanXPi.push_back(mean_x_pi);
        tr_pidMeanYPi.push_back(mean_y_pi);
        tr_pidSigmXPi.push_back(sigm_x_pi);
        tr_pidSigmYPi.push_back(sigm_y_pi);
        tr_pidMeanXKa.push_back(mean_x_ka);
        tr_pidMeanYKa.push_back(mean_y_ka);
        tr_pidSigmXKa.push_back(sigm_x_ka);
        tr_pidSigmYKa.push_back(sigm_y_ka);
        tr_pidMeanXPr.push_back(mean_x_pr);
        tr_pidMeanYPr.push_back(mean_y_pr);
        tr_pidSigmXPr.push_back(sigm_x_pr);
        tr_pidSigmYPr.push_back(sigm_y_pr);
      } else {
        tr_pidX.push_back(-999.);
        tr_pidY.push_back(-999.);
        tr_pidMeanXPi.push_back(-999.);
        tr_pidMeanYPi.push_back(-999.);
        tr_pidSigmXPi.push_back(-999.);
        tr_pidSigmYPi.push_back(-999.);
        tr_pidMeanXKa.push_back(-999.);
        tr_pidMeanYKa.push_back(-999.);
        tr_pidSigmXKa.push_back(-999.);
        tr_pidSigmYKa.push_back(-999.);
        tr_pidMeanXPr.push_back(-999.);
        tr_pidMeanYPr.push_back(-999.);
        tr_pidSigmXPr.push_back(-999.);
        tr_pidSigmYPr.push_back(-999.);
      }

		} // track loop

    // BBC loop
    bbc_id.clear();
    bbc_side.clear();
    bbc_en.clear();
    bbc_phi.clear();
    bbc_id.resize(fNBbcModules*2);
    bbc_side.resize(fNBbcModules*2);
    bbc_en.resize(fNBbcModules*2);
    bbc_phi.resize(fNBbcModules*2);
    for (int iBbcMod=0; iBbcMod<fNBbcModules; iBbcMod++)
    {

    } // Bbc loop

    // ZDC loop
    zdc_id.clear();
    zdc_side.clear();
    zdc_type.clear();
    zdc_en.clear();
    zdc_phi.clear();
    zdc_id.resize(2*(fNZdcSmdStripsHorizontal+fNZdcSmdStripsVertical));
    zdc_side.resize(2*(fNZdcSmdStripsHorizontal+fNZdcSmdStripsVertical));
    zdc_type.resize(2*(fNZdcSmdStripsHorizontal+fNZdcSmdStripsVertical));
    zdc_en.resize(2*(fNZdcSmdStripsHorizontal+fNZdcSmdStripsVertical));
    zdc_phi.resize(2*(fNZdcSmdStripsHorizontal+fNZdcSmdStripsVertical));
    // Horizontal
    for (int iZdcStrip=0; iZdcStrip<fNZdcSmdStripsHorizontal; iZdcStrip++)
    {

    } // Bbc loop - horizontal
    // Vertical
    for (int iZdcStrip=0; iZdcStrip<fNZdcSmdStripsVertical; iZdcStrip++)
    {

    } // Bbc loop - vertical

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
