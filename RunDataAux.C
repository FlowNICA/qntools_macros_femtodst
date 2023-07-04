#include "utils.h"

bool isGoodEvent(StFemtoEvent *const& event);
bool isGoodTrack(StFemtoTrack *const& track, StFemtoEvent *const& event);

void RunDataAux(std::string iFileName, std::string oFileName)
{
  TStopwatch timer;
  timer.Start();

  const int Ncentralities = 16;
  const int NptBins = 15;
  const double ptBins[NptBins+1] = {0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2};

  StFemtoDstReader* femtoReader = new StFemtoDstReader(iFileName.c_str());
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
    return;
  }
	Long64_t events2read = femtoReader->chain()->GetEntries();

	float ev_cent, ev_vtxX, ev_vtxY, ev_vtxZ, ev_vpdZ, ev_runId;

  TFile *fo = new TFile(oFileName.c_str(), "recreate");
  TH1D *h_pt[Ncentralities];
  for (int icent = 0; icent < Ncentralities; icent++) {
    h_pt[icent] = new TH1D(Form("h_pt_cent%i",icent), Form("h_pt for %3.0f-%3.0f%s", (double)icent*5., (double)icent*5+5., "%"), 5000, 0., 5.);
  }

  TH2D *h_NsigPiMsqr[NptBins];
  TH2D *h_NsigKaMsqr[NptBins];
  TH2D *h_NsigPrMsqr[NptBins];

  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_NsigPiMsqr[ipt] = new TH2D(Form("h_NsigPiMsqr_pt%i", ipt), Form("h_NsigPiMsqr for %1.1f < p_{T} < %1.1f GeV/c;m^{2}, (Gev/c^{2})^{2}; n#sigma_{#pi}", ptBins[ipt], ptBins[ipt+1]), 2000, -0.5, 1.5, 6000, -30., 30.);
    h_NsigKaMsqr[ipt] = new TH2D(Form("h_NsigKaMsqr_pt%i", ipt), Form("h_NsigKaMsqr for %1.1f < p_{T} < %1.1f GeV/c;m^{2}, (Gev/c^{2})^{2}; n#sigma_{K}",   ptBins[ipt], ptBins[ipt+1]), 2000, -0.5, 1.5, 6000, -30., 30.);
    h_NsigPrMsqr[ipt] = new TH2D(Form("h_NsigPrMsqr_pt%i", ipt), Form("h_NsigPrMsqr for %1.1f < p_{T} < %1.1f GeV/c;m^{2}, (Gev/c^{2})^{2}; n#sigma_{p}",   ptBins[ipt], ptBins[ipt+1]), 2000, -0.5, 1.5, 6000, -30., 30.);
  }

  std::cout << "Number of events to read: " << events2read << std::endl;
	// Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

		if (iEvent % 1000 == 0)
		{
			std::cout << "Preparing event #[" << (iEvent+1)
					      << "/" << events2read << "]" << std::endl;
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

    if (!isGoodEvent(event)) continue;

		ev_cent = (15.-(float)event->cent16())*5.+2.5;
    int centBin = (int)(15-(int)event->cent16());
    if (centBin < 0) continue;
		
		// Track analysis
    Int_t nTracks = dst->numberOfTracks();

    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);

      if (!isGoodTrack(femtoTrack, event)) continue;

      TVector3 mom = femtoTrack->gMom();

      h_pt[centBin]->Fill(mom.Pt());

      int ipt_bin = -1;
      for (int ipt = 0; ipt < NptBins; ipt++) {
        if (mom.Pt() > ptBins[ipt] && mom.Pt() <= ptBins[ipt+1]) {
          ipt_bin = ipt;
          // break;
        }
      }
      if (ipt_bin < 0) continue;
      h_NsigPiMsqr[ipt_bin]->Fill(femtoTrack->massSqr(), femtoTrack->nSigmaPion());
      h_NsigKaMsqr[ipt_bin]->Fill(femtoTrack->massSqr(), femtoTrack->nSigmaKaon());
      h_NsigPrMsqr[ipt_bin]->Fill(femtoTrack->massSqr(), femtoTrack->nSigmaProton());
    }
  }

  fo->cd();
  for (int icent = 0; icent < Ncentralities; icent++) {
    h_pt[icent]->Write();
  }
  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_NsigPiMsqr[ipt]->Write();
    h_NsigKaMsqr[ipt]->Write();
    h_NsigPrMsqr[ipt]->Write();
  }
  fo->Close();

	std::cout << "Neccessary data collected!" << std::endl;
  femtoReader->Finish();
  timer.Stop();
  timer.Print();
}


bool isGoodEvent(StFemtoEvent *const& event)
{
  TVector3 pVtx = event->primaryVertex();

  if (abs(pVtx.Z()) >= 30) return false;
  if (sqrt(pVtx.X()*pVtx.X()+pVtx.Y()*pVtx.Y()) >= 2.) return false;
  if (abs(pVtx.Z()-event->vpdVz()) >= 3.) return false;

  return true;
}

bool isGoodTrack(StFemtoTrack *const& track, StFemtoEvent *const& event)
{
  TVector3 mom = track->gMom();
  TVector3 pVtx = event->primaryVertex();
  TVector3 dca = track->gDCA(pVtx);

  if ( !track->isPrimary() ) return false;
  if (abs(mom.Eta()) >= 1.) return false;
  if (mom.Mag() >= 10.) return false;
  if (track->nHitsFit() <= 15.) return false;
  if (track->nHitsPoss() <= 0.) return false;
  if ((double)track->nHitsFit()/(double)track->nHitsPoss() <= 0.51) return false;
  if (dca.Mag() >= 1.) return false;

  return true;
}
