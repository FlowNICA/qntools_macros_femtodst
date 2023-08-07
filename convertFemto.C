#include "utils.h"
#include <fstream>

void convertFemto(std::string iFileName, std::string oFileName, std::string pidFileName)
{
	TStopwatch timer;
	timer.Start();
  TFile *fo;
	TTree *tree;
	/*
	tree = (TTree*) makeTree4RDF(iFileName);
	fo = new TFile(oFileName.c_str(),"recreate");
	fo->cd();
	tree->Write();
	fo->Close();
	*/

  TFile *fiPid = new TFile(pidFileName.c_str(),"read");
  TGraphErrors *gr_mean_m2[3], *gr_mean_ns[3], *gr_sigm_m2[3], *gr_sigm_ns[3];
  TGraphErrors *gr_mean_x[3],  *gr_mean_y[3],  *gr_sigm_x[3],  *gr_sigm_y[3];

  for (int i=0; i<3; i++) {
    gr_mean_m2[i] = (TGraphErrors*)fiPid->Get(Form("nsig_m2_1D/gr_orig_signal2_mean_m2_pid%i_smoothed", i));
    gr_mean_ns[i] = (TGraphErrors*)fiPid->Get(Form("nsig_m2_1D/gr_orig_signal2_mean_ns_pid%i_smoothed", i));
    gr_sigm_m2[i] = (TGraphErrors*)fiPid->Get(Form("nsig_m2_1D/gr_orig_signal2_sigm_m2_pid%i_smoothed", i));
    gr_sigm_ns[i] = (TGraphErrors*)fiPid->Get(Form("nsig_m2_1D/gr_orig_signal2_sigm_ns_pid%i_smoothed", i));
    gr_mean_x[i]  = (TGraphErrors*)fiPid->Get(Form("x_y_1D/gr_orig_signal3_mean_x_pid%i_smoothed", i));
    gr_mean_y[i]  = (TGraphErrors*)fiPid->Get(Form("x_y_1D/gr_orig_signal2_mean_y_pid%i_smoothed", i));
    gr_sigm_x[i]  = (TGraphErrors*)fiPid->Get(Form("x_y_1D/gr_orig_signal3_sigm_x_pid%i_smoothed", i));
    gr_sigm_y[i]  = (TGraphErrors*)fiPid->Get(Form("x_y_1D/gr_orig_signal2_sigm_y_pid%i_smoothed", i));
  }
	
	std::string line;
	std::ifstream infile(iFileName);
	std::fstream outfile;
	outfile.open(Form("%s.list",oFileName.c_str()), ios::out);
	int icount = 0;
	while (std::getline(infile, line))
	{
		std::cout << "Processing file " << icount+1 << ": " << line.c_str() << std::endl;
    TFile *fo = new TFile(Form("%s_%i.root",oFileName.c_str(), icount),"recreate");
    fo->cd();
		TTree *tree = (TTree*) makeTree4RDF(line,
                                        gr_mean_m2[0], gr_mean_ns[0], gr_sigm_m2[0], gr_sigm_ns[0],
                                        gr_mean_m2[1], gr_mean_ns[1], gr_sigm_m2[1], gr_sigm_ns[1],
                                        gr_mean_m2[2], gr_mean_ns[2], gr_sigm_m2[2], gr_sigm_ns[2],
                                        gr_mean_x[0],  gr_mean_y[0],  gr_sigm_x[0],  gr_sigm_y[0],
                                        gr_mean_x[1],  gr_mean_y[1],  gr_sigm_x[1],  gr_sigm_y[1],
                                        gr_mean_x[2],  gr_mean_y[2],  gr_sigm_x[2],  gr_sigm_y[2]);
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
