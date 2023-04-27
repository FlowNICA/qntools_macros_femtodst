#include "utils.h"
#include <fstream>

void convertFemto(std::string iFileName, std::string oFileName)
{
	TStopwatch timer;
	timer.Start();
	TTree *tree;
	TFile *fo;
	/*
	tree = (TTree*) makeTree4RDF(iFileName);
	fo = new TFile(oFileName.c_str(),"recreate");
	fo->cd();
	tree->Write();
	fo->Close();
	*/
	
	std::string line;
	std::ifstream infile(iFileName);
	std::fstream outfile;
	outfile.open(Form("%s.list",oFileName.c_str()), ios::out);
	int icount = 0;
	while (std::getline(infile, line))
	{
		std::cout << "Processing file: " << line.c_str() << std::endl;
		tree = (TTree*) makeTree4RDF(line);
		fo = new TFile(Form("%s_%i.root",oFileName.c_str(), icount),"recreate");
		fo->cd();
		tree->Write();
		fo->Close();
		outfile << Form("%s_%i.root",oFileName.c_str(), icount) << std::endl;
		icount++;
	}
	outfile.close();
	infile.close();
	

	timer.Stop();
	timer.Print();
}
