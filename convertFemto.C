#include "utils.h"

void convertFemto(std::string iFileName, std::string oFileName)
{
	TTree *tree = (TTree*) makeTree4RDF(iFileName);
	TFile *fo = new TFile(oFileName.c_str(),"recreate");
	fo->cd();
	tree->Write();
	fo->Close();
}
