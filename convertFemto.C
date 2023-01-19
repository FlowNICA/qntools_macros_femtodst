#include "utils.h"

void convertFemto(std::string iFileName, std::string oFileName)
{
	TStopwatch timer;
	timer.Start();
	TTree *tree = (TTree*) makeTree4RDF(iFileName);
	TFile *fo = new TFile(oFileName.c_str(),"recreate");
	fo->cd();
	tree->Write();
	fo->Close();
	timer.Stop();
	timer.Print();
}
