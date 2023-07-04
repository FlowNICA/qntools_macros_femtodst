#include <iostream>
#include <vector>

#include <TStopwatch.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

void RunPidFits(string iFileName, string oFileName)
{
  TStopwatch timer;
  timer.Start();

  const int Ncentralities = 16;
  const int NptBins = 15;
  const double ptBins[NptBins+1] = {0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2};

  const vector<double> pid_mass = {0.13957, 0.493677, 0.938272}; //pi, K, p masses in GeV/c2
  const vector<double> pid_Msqr = {0.13957*0.13957, 0.493677*0.493677, 0.938272*0.938272}; //pi, K, p m^2 in (GeV/c^2)^2
  const int Npid = 3; // pi, K, p

  const int Niter = 1; // for more precise fits

  const vector<vector<double>> NsigMeanInitial = {{0., 10., 20.}, {0., 4., 12.}, {0., 2., 8.}, {0., 1., 5.},                     // 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.0
                                                  {0., 0., 3.}, {0., -1., 1.}, {0., -1., 0.}, {0., -1.2, -0.5}, {0., -1.5, -1.}, // 1.0-1.2, 1.2-1.4, 1.4-1.6, 1.6-1.8, 1.8-2.0
                                                  {0., -1.5, -1.5}, {0.,  -1.5, -1.5}, {0.,  -1.5, -2.}, {0.,  -1.5, -2.}, {0.,  -1.5, -2.}, // 2.0-2.2, 2.2-2.4, 2.4-2.6, 2.6-2.8, 2.8-3.0
                                                  {0., -2., -2.5}};                                                              // 3.0-3.2

  // Open aux file
  cout << "Open input file: " << iFileName << endl;
  TFile *fi = new TFile(iFileName.c_str(), "read");
  TFile *fo = new TFile(oFileName.c_str(), "recreate");

  TH2D *h_NsigPiMsqr[NptBins];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_NsigPiMsqr[ipt] = (TH2D*) fi->Get(Form("h_NsigPiMsqr_pt%i", ipt));
    if (!h_NsigPiMsqr[ipt]) {
      cerr << "Warning: " << Form("h_NsigPiMsqr_pt%i", ipt) << " was not found in the input file! Abort!" << endl;
      return;
    }
  }

  // Do a 1D projections and 1D fits
  cout << "Fitting 1D projections of NsigmaPi and M^2..." << endl;
  TH1D *h_NsigPi[NptBins][Npid]; // slice NsigPi for each particle
  TH1D *h_Msqr[NptBins];

  TF1 *f1_gaus1_Msqr[NptBins][Npid];
  TF1 *f1_gaus2_Msqr[NptBins][Npid];
  TF1 *f1_gaus3_Msqr[NptBins];
  TF1 *f1_gaus1_Nsig[NptBins][Npid];
  TF1 *f1_gaus2_Nsig[NptBins][Npid];
  TF1 * f_bg_Msqr[NptBins];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    for (int ipid=0; ipid<Npid; ipid++) {
      f1_gaus1_Msqr[ipt][ipid] = new TF1(Form("f1_gaus1_Msqr_pt%i_pid%i", ipt, ipid), "gaus", pid_Msqr.at(ipid)*0.9, pid_Msqr.at(ipid)*1.1);
      f1_gaus1_Msqr[ipt][ipid]->SetParameter(1, pid_Msqr.at(ipid));
    }
  }

  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_Msqr[ipt] = (TH1D*) h_NsigPiMsqr[ipt]->ProjectionX(Form("h_Msqr_pt%i", ipt), 0, -1);
    for (int ipid=0; ipid<Npid; ipid++) {
      for (int iter=0; iter<Niter; iter++) {
        h_Msqr[ipt]->Fit(f1_gaus1_Msqr[ipt][ipid], "MR0Q");
      }
      h_NsigPi[ipt][ipid] = (TH1D*) h_NsigPiMsqr[ipt]->ProjectionY(Form("h_NsigPi_pt%i_pid%i", ipt, ipid), 
        h_NsigPiMsqr[ipt]->GetXaxis()->FindBin(f1_gaus1_Msqr[ipt][ipid]->GetParameter(1)-0.5*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2)),
        h_NsigPiMsqr[ipt]->GetXaxis()->FindBin(f1_gaus1_Msqr[ipt][ipid]->GetParameter(1)+0.5*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2)));
      f1_gaus1_Nsig[ipt][ipid] = new TF1(Form("f1_gaus1_Nsig_pt%i_pid%i", ipt, ipid), "gaus", 
        NsigMeanInitial.at(ipt).at(ipid) - 1., 
        NsigMeanInitial.at(ipt).at(ipid) + 1.);
      for (int iter=0; iter<Niter; iter++) {
        h_NsigPi[ipt][ipid]->Fit(f1_gaus1_Nsig[ipt][ipid], "MR0Q");
      }
    }
  }

  // Do a 2x1D fits
  cout << "Fitting 2x1D projections of NsigmaPi and M^2..." << endl;
  for (int ipt = 0; ipt < NptBins; ipt++) {
    f1_gaus3_Msqr[ipt] = new TF1(Form("f1_gaus1_Msqr_pt%i", ipt), "gaus(0)+gaus(3) + gaus(6)+gaus(9) + gaus(12)+gaus(15)", -0.2, 1.5);
    for (int ipid=0; ipid<Npid; ipid++) {
      f1_gaus2_Msqr[ipt][ipid] = new TF1(Form("f1_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "gaus(0)+gaus(3)", pid_Msqr.at(ipid)*0.7, pid_Msqr.at(ipid)*1.3);

      f1_gaus2_Msqr[ipt][ipid]->SetParName(0, "SglConstant");
      f1_gaus2_Msqr[ipt][ipid]->SetParName(1, "SglMean");
      f1_gaus2_Msqr[ipt][ipid]->SetParName(2, "SglSigma");
      f1_gaus2_Msqr[ipt][ipid]->SetParName(3, "BgdConstant");
      f1_gaus2_Msqr[ipt][ipid]->SetParName(4, "BgdMean");
      f1_gaus2_Msqr[ipt][ipid]->SetParName(5, "BgdSigma");

      f1_gaus2_Msqr[ipt][ipid]->SetParameter(0, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0));
      f1_gaus2_Msqr[ipt][ipid]->SetParameter(1, f1_gaus1_Msqr[ipt][ipid]->GetParameter(1));
      f1_gaus2_Msqr[ipt][ipid]->SetParameter(2, f1_gaus1_Msqr[ipt][ipid]->GetParameter(2));
      f1_gaus2_Msqr[ipt][ipid]->SetParameter(3, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0));
      f1_gaus2_Msqr[ipt][ipid]->SetParameter(4, f1_gaus1_Msqr[ipt][ipid]->GetParameter(1));
      f1_gaus2_Msqr[ipt][ipid]->SetParameter(5, 3.*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2));

      f1_gaus2_Msqr[ipt][ipid]->SetParLimits(0, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0)*1.5);
      if (ipid == 0) f1_gaus2_Msqr[ipt][ipid]->SetParLimits(1, pid_Msqr.at(ipid)*0., pid_Msqr.at(ipid)*5.);
      else f1_gaus2_Msqr[ipt][ipid]->SetParLimits(1, pid_Msqr.at(ipid)*0.5, pid_Msqr.at(ipid)*1.5);
      f1_gaus2_Msqr[ipt][ipid]->SetParLimits(2, f1_gaus1_Msqr[ipt][ipid]->GetParameter(2)*0.5, f1_gaus1_Msqr[ipt][ipid]->GetParameter(2)*1.5);
      f1_gaus2_Msqr[ipt][ipid]->SetParLimits(3, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0)*1.5);
      if (ipid == 0) f1_gaus2_Msqr[ipt][ipid]->SetParLimits(4, pid_Msqr.at(ipid)*0., pid_Msqr.at(ipid)*5.);
      else f1_gaus2_Msqr[ipt][ipid]->SetParLimits(4, pid_Msqr.at(ipid)*0.5, pid_Msqr.at(ipid)*1.5);
      f1_gaus2_Msqr[ipt][ipid]->SetParLimits(5, 1.5*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2), 1000.*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2));

      for (int iter=0; iter<Niter; iter++) {
        h_Msqr[ipt]->Fit(f1_gaus2_Msqr[ipt][ipid], "MR0Q");
      }

      f1_gaus2_Nsig[ipt][ipid] = new TF1(Form("f1_gaus2_Nsig_pt%i_pid%i", ipt, ipid), "gaus(0)+gaus(3)", 
        f1_gaus1_Nsig[ipt][ipid]->GetParameter(1)-3.*f1_gaus1_Nsig[ipt][ipid]->GetParameter(2), 
        f1_gaus1_Nsig[ipt][ipid]->GetParameter(1)+3.*f1_gaus1_Nsig[ipt][ipid]->GetParameter(2));
      
      f1_gaus2_Nsig[ipt][ipid]->SetParName(0, "SglConstant");
      f1_gaus2_Nsig[ipt][ipid]->SetParName(1, "SglMean");
      f1_gaus2_Nsig[ipt][ipid]->SetParName(2, "SglSigma");
      f1_gaus2_Nsig[ipt][ipid]->SetParName(3, "BgdConstant");
      f1_gaus2_Nsig[ipt][ipid]->SetParName(4, "BgdMean");
      f1_gaus2_Nsig[ipt][ipid]->SetParName(5, "BgdSigma");

      f1_gaus2_Nsig[ipt][ipid]->SetParameter(0, f1_gaus1_Nsig[ipt][ipid]->GetParameter(0));
      if (ipid == 0) f1_gaus2_Nsig[ipt][ipid]->SetParameter(1, 0.);
      else f1_gaus2_Nsig[ipt][ipid]->SetParameter(1, f1_gaus1_Nsig[ipt][ipid]->GetParameter(1));
      f1_gaus2_Nsig[ipt][ipid]->SetParameter(2, f1_gaus1_Nsig[ipt][ipid]->GetParameter(2));
      f1_gaus2_Nsig[ipt][ipid]->SetParameter(3, f1_gaus1_Nsig[ipt][ipid]->GetParameter(0));
      if (ipid == 0) f1_gaus2_Nsig[ipt][ipid]->SetParameter(4, 0.);
      else f1_gaus2_Nsig[ipt][ipid]->SetParameter(4, f1_gaus1_Nsig[ipt][ipid]->GetParameter(1));
      f1_gaus2_Nsig[ipt][ipid]->SetParameter(5, 3.*f1_gaus1_Nsig[ipt][ipid]->GetParameter(2));

      // f1_gaus2_Nsig[ipt][ipid]->SetParLimits(1, f1_gaus1_Nsig[ipt][ipid]->GetParameter(1)*0.5, f1_gaus1_Nsig[ipt][ipid]->GetParameter(1)*2.);
      // f1_gaus2_Nsig[ipt][ipid]->SetParLimits(2, f1_gaus1_Nsig[ipt][ipid]->GetParameter(2)*0.5, f1_gaus1_Nsig[ipt][ipid]->GetParameter(2)*2.);
      // f1_gaus2_Nsig[ipt][ipid]->SetParLimits(3, f1_gaus1_Nsig[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_Nsig[ipt][ipid]->GetParameter(0)*2.);
      // f1_gaus2_Nsig[ipt][ipid]->SetParLimits(4, f1_gaus1_Nsig[ipt][ipid]->GetParameter(1)*0.5, f1_gaus1_Nsig[ipt][ipid]->GetParameter(1)*2.);
      // f1_gaus2_Nsig[ipt][ipid]->SetParLimits(5, 1.5*f1_gaus1_Nsig[ipt][ipid]->GetParameter(2), 10.*f1_gaus1_Nsig[ipt][ipid]->GetParameter(2));

      for (int iter=0; iter<Niter; iter++) {
        h_NsigPi[ipt][ipid]->Fit(f1_gaus2_Nsig[ipt][ipid], "MR0Q");
      }
    }
  }

  // // Do a 2x1D 3gaus fits - check
  // cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  // cout << "Fitting 2x1D 3gaus M^2... - check" << endl;
  // for (int ipt = 0; ipt < NptBins; ipt++) {
  //   // cout << "\tPt bin " << ipt <<endl;
  //   f1_gaus3_Msqr[ipt] = new TF1(Form("f1_gaus3_Msqr_pt%i", ipt), "gaus(0)+gaus(3) + gaus(6)+gaus(9) + gaus(12)+gaus(15)", 
  //     // -0.2, 1.5
  //     f1_gaus2_Msqr[ipt][0]->GetParameter(1)-3.*f1_gaus2_Msqr[ipt][0]->GetParameter(5),
  //     f1_gaus2_Msqr[ipt][2]->GetParameter(1)+3.*f1_gaus2_Msqr[ipt][2]->GetParameter(5)
  //   );

  //   f1_gaus3_Msqr[ipt]->SetParName(0, "Peak1SglConstant");
  //   f1_gaus3_Msqr[ipt]->SetParName(1, "Peak1SglMean");
  //   f1_gaus3_Msqr[ipt]->SetParName(2, "Peak1SglSigma");
  //   f1_gaus3_Msqr[ipt]->SetParName(3, "Peak1BgdConstant");
  //   f1_gaus3_Msqr[ipt]->SetParName(4, "Peak1BgdMean");
  //   f1_gaus3_Msqr[ipt]->SetParName(5, "Peak1BgdSigma");
  //   f1_gaus3_Msqr[ipt]->SetParName(6, "Peak2SglConstant");
  //   f1_gaus3_Msqr[ipt]->SetParName(7, "Peak2SglMean");
  //   f1_gaus3_Msqr[ipt]->SetParName(8, "Peak2SglSigma");
  //   f1_gaus3_Msqr[ipt]->SetParName(9, "Peak2BgdConstant");
  //   f1_gaus3_Msqr[ipt]->SetParName(10, "Peak2BgdMean");
  //   f1_gaus3_Msqr[ipt]->SetParName(11, "Peak2BgdSigma");
  //   f1_gaus3_Msqr[ipt]->SetParName(12, "Peak3SglConstant");
  //   f1_gaus3_Msqr[ipt]->SetParName(13, "Peak3SglMean");
  //   f1_gaus3_Msqr[ipt]->SetParName(14, "Peak3SglSigma");
  //   f1_gaus3_Msqr[ipt]->SetParName(15, "Peak3BgdConstant");
  //   f1_gaus3_Msqr[ipt]->SetParName(16, "Peak3BgdMean");
  //   f1_gaus3_Msqr[ipt]->SetParName(17, "Peak3BgdSigma");

  //   f1_gaus3_Msqr[ipt]->SetParameter(0, f1_gaus2_Msqr[ipt][0]->GetParameter(0));
  //   f1_gaus3_Msqr[ipt]->SetParameter(1, f1_gaus2_Msqr[ipt][0]->GetParameter(1));
  //   f1_gaus3_Msqr[ipt]->SetParameter(2, f1_gaus2_Msqr[ipt][0]->GetParameter(2));
  //   f1_gaus3_Msqr[ipt]->SetParameter(3, f1_gaus2_Msqr[ipt][0]->GetParameter(3));
  //   f1_gaus3_Msqr[ipt]->SetParameter(4, f1_gaus2_Msqr[ipt][0]->GetParameter(4));
  //   f1_gaus3_Msqr[ipt]->SetParameter(5, f1_gaus2_Msqr[ipt][0]->GetParameter(5));
  //   f1_gaus3_Msqr[ipt]->SetParameter(6, f1_gaus2_Msqr[ipt][1]->GetParameter(0));
  //   f1_gaus3_Msqr[ipt]->SetParameter(7, f1_gaus2_Msqr[ipt][1]->GetParameter(1));
  //   f1_gaus3_Msqr[ipt]->SetParameter(8, f1_gaus2_Msqr[ipt][1]->GetParameter(2));
  //   f1_gaus3_Msqr[ipt]->SetParameter(9, f1_gaus2_Msqr[ipt][1]->GetParameter(3));
  //   f1_gaus3_Msqr[ipt]->SetParameter(10, f1_gaus2_Msqr[ipt][1]->GetParameter(4));
  //   f1_gaus3_Msqr[ipt]->SetParameter(11, f1_gaus2_Msqr[ipt][1]->GetParameter(5));
  //   f1_gaus3_Msqr[ipt]->SetParameter(12, f1_gaus2_Msqr[ipt][2]->GetParameter(0));
  //   f1_gaus3_Msqr[ipt]->SetParameter(13, f1_gaus2_Msqr[ipt][2]->GetParameter(1));
  //   f1_gaus3_Msqr[ipt]->SetParameter(14, f1_gaus2_Msqr[ipt][2]->GetParameter(2));
  //   f1_gaus3_Msqr[ipt]->SetParameter(15, f1_gaus2_Msqr[ipt][2]->GetParameter(3));
  //   f1_gaus3_Msqr[ipt]->SetParameter(16, f1_gaus2_Msqr[ipt][2]->GetParameter(4));
  //   f1_gaus3_Msqr[ipt]->SetParameter(17, f1_gaus2_Msqr[ipt][2]->GetParameter(5));

  //   f1_gaus3_Msqr[ipt]->SetParLimits(0,  f1_gaus2_Msqr[ipt][0]->GetParameter(0)*0.5, f1_gaus2_Msqr[ipt][0]->GetParameter(0)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(1,  f1_gaus2_Msqr[ipt][0]->GetParameter(1)*0.5, f1_gaus2_Msqr[ipt][0]->GetParameter(1)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(2,  f1_gaus2_Msqr[ipt][0]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][0]->GetParameter(2)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(3,  f1_gaus2_Msqr[ipt][0]->GetParameter(3)*0.5, f1_gaus2_Msqr[ipt][0]->GetParameter(3)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(4,  f1_gaus2_Msqr[ipt][0]->GetParameter(4)*0.5, f1_gaus2_Msqr[ipt][0]->GetParameter(3)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(5,  f1_gaus2_Msqr[ipt][0]->GetParameter(5)*0.5, f1_gaus2_Msqr[ipt][0]->GetParameter(5)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(6,  f1_gaus2_Msqr[ipt][1]->GetParameter(0)*0.5, f1_gaus2_Msqr[ipt][1]->GetParameter(0)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(7,  f1_gaus2_Msqr[ipt][1]->GetParameter(1)*0.5, f1_gaus2_Msqr[ipt][1]->GetParameter(1)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(8,  f1_gaus2_Msqr[ipt][1]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][1]->GetParameter(2)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(9,  f1_gaus2_Msqr[ipt][1]->GetParameter(3)*0.5, f1_gaus2_Msqr[ipt][1]->GetParameter(3)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(10, f1_gaus2_Msqr[ipt][1]->GetParameter(4)*0.5, f1_gaus2_Msqr[ipt][1]->GetParameter(4)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(11, f1_gaus2_Msqr[ipt][1]->GetParameter(5)*0.5, f1_gaus2_Msqr[ipt][1]->GetParameter(5)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(12, f1_gaus2_Msqr[ipt][2]->GetParameter(0)*0.5, f1_gaus2_Msqr[ipt][2]->GetParameter(0)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(13, f1_gaus2_Msqr[ipt][2]->GetParameter(1)*0.5, f1_gaus2_Msqr[ipt][2]->GetParameter(1)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(14, f1_gaus2_Msqr[ipt][2]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][2]->GetParameter(2)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(15, f1_gaus2_Msqr[ipt][2]->GetParameter(3)*0.5, f1_gaus2_Msqr[ipt][2]->GetParameter(3)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(16, f1_gaus2_Msqr[ipt][2]->GetParameter(4)*0.5, f1_gaus2_Msqr[ipt][2]->GetParameter(4)*1.5);
  //   f1_gaus3_Msqr[ipt]->SetParLimits(17, f1_gaus2_Msqr[ipt][2]->GetParameter(5)*0.5, f1_gaus2_Msqr[ipt][2]->GetParameter(5)*1.5);

  //   for (int iter=0; iter<Niter; iter++) {
  //     h_Msqr[ipt]->Fit(f1_gaus3_Msqr[ipt], "MR0Q");
  //   }
  //   // timer.Print();
  //   // timer.Continue();
  // }

  // // Do 3 1x2D gaus fits - check
  // cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  // cout << "Fitting 3 2D of NsigmaPi and M^2... - check" << endl;
  // TF2 *f2_gaus1_Msqr[NptBins][Npid];
  // for (int ipt = 0; ipt < NptBins; ipt++) {
  //   for (int ipid=0; ipid<Npid; ipid++) {
  //     f2_gaus1_Msqr[ipt][ipid] = new TF2(Form("f2_gaus1_Msqr_pt%i_pid%i", ipt, ipid), "xygaus", 
  //       f1_gaus2_Msqr[ipt][ipid]->GetParameter(1)-0.5*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
  //       f1_gaus2_Msqr[ipt][ipid]->GetParameter(1)+0.5*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
  //       f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)-0.5*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2),
  //       f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)+0.5*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));

  //     f2_gaus1_Msqr[ipt][ipid]->SetParName(0, "Constant");
  //     f2_gaus1_Msqr[ipt][ipid]->SetParName(1, "MeanX");
  //     f2_gaus1_Msqr[ipt][ipid]->SetParName(2, "SigmaX");
  //     f2_gaus1_Msqr[ipt][ipid]->SetParName(3, "MeanY");
  //     f2_gaus1_Msqr[ipt][ipid]->SetParName(4, "SigmaY");
      
  //     f2_gaus1_Msqr[ipt][ipid]->SetParameter(0, f1_gaus2_Msqr[ipt][ipid]->GetParameter(0));
  //     f2_gaus1_Msqr[ipt][ipid]->SetParameter(1, f1_gaus2_Msqr[ipt][ipid]->GetParameter(1));
  //     // f2_gaus1_Msqr[ipt][ipid]->FixParameter(1, pid_Msqr.at(ipid));
  //     f2_gaus1_Msqr[ipt][ipid]->SetParameter(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2));
  //     f2_gaus1_Msqr[ipt][ipid]->SetParameter(3, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1));
  //     f2_gaus1_Msqr[ipt][ipid]->SetParameter(4, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));

  //     f2_gaus1_Msqr[ipt][ipid]->SetParLimits(0, 0., h_Msqr[ipt]->GetEntries());
  //     f2_gaus1_Msqr[ipt][ipid]->SetParLimits(1, f1_gaus2_Msqr[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(1)*1.5);
  //     f2_gaus1_Msqr[ipt][ipid]->SetParLimits(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*1.5);
  //     f2_gaus1_Msqr[ipt][ipid]->SetParLimits(3, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*1.5);
  //     f2_gaus1_Msqr[ipt][ipid]->SetParLimits(4, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2)*1.5);
  //   }
  // }
  // for (int ipt = 0; ipt < NptBins; ipt++) {
  //   cout << "\tPt bin " << ipt <<endl;
  //   for (int ipid=0; ipid<Npid; ipid++) {
  //     // cout << "\t\tipid bin " << ipid << endl;
  //     for (int iter=0; iter<Niter; iter++) {
  //       // if (iter != Niter-1) 
  //         h_NsigPiMsqr[ipt]->Fit(f2_gaus1_Msqr[ipt][ipid],"MR0Q");
  //       // else h_NsigPiMsqr[ipt]->Fit(f2_gaus1_Msqr[ipt][ipid],"MR");
  //     }
  //   }
  //   timer.Print();
  //   timer.Continue();
  // }

  // Do 3 2x2D gaus fits
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fitting 3 2x2D of NsigmaPi and M^2..." << endl;
  TF2 *f2_gaus2_Msqr[NptBins][Npid];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    for (int ipid=0; ipid<Npid; ipid++) {
      if (ipid == 0 && ipt < 8)
        f2_gaus2_Msqr[ipt][ipid] = new TF2(Form("f2_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          pid_Msqr.at(ipid)-4.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          pid_Msqr.at(ipid)+4.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)-4.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)+4.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));
      if (ipid == 0 && ipt >= 8)
        f2_gaus2_Msqr[ipt][ipid] = new TF2(Form("f2_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          pid_Msqr.at(ipid)-3.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          pid_Msqr.at(ipid)+3.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)-4.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)+4.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));
      else if (ipid == 1 && ipt < 8)
        f2_gaus2_Msqr[ipt][ipid] = new TF2(Form("f2_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          pid_Msqr.at(ipid)-4.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          pid_Msqr.at(ipid)+4.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));
      else if (ipid == 1 && ipt >= 8)
        f2_gaus2_Msqr[ipt][ipid] = new TF2(Form("f2_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          pid_Msqr.at(ipid)-2.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          pid_Msqr.at(ipid)+2.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));
      else if (ipid == 2 && ipt < 8)
        f2_gaus2_Msqr[ipt][ipid] = new TF2(Form("f2_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          pid_Msqr.at(ipid)-4.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          pid_Msqr.at(ipid)+4.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));
      else
        f2_gaus2_Msqr[ipt][ipid] = new TF2(Form("f2_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          pid_Msqr.at(ipid)-2.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          pid_Msqr.at(ipid)+2.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));

      f2_gaus2_Msqr[ipt][ipid]->SetParName(0, "SglConstant");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(1, "SglMeanX");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(2, "SglSigmaX");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(3, "SglMeanY");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(4, "SglSigmaY");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(5, "BgdConstant");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(6, "BgdMeanX");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(7, "BgdSigmaX");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(8, "BgdMeanY");
      f2_gaus2_Msqr[ipt][ipid]->SetParName(9, "BgdSigmaY");
      
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(0, f1_gaus2_Msqr[ipt][ipid]->GetParameter(0));
      // if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParameter(1, pid_Msqr.at(ipid));
      // else f2_gaus2_Msqr[ipt][ipid]->SetParameter(1, pid_Msqr.at(ipid));
      f2_gaus2_Msqr[ipt][ipid]->FixParameter(1, pid_Msqr.at(ipid));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(3, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(4, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(3));
      // if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParameter(6, pid_Msqr.at(ipid));
      // else f2_gaus2_Msqr[ipt][ipid]->SetParameter(6, pid_Msqr.at(ipid));
      f2_gaus2_Msqr[ipt][ipid]->FixParameter(6, pid_Msqr.at(ipid));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(7, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(8, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(9, f1_gaus2_Nsig[ipt][ipid]->GetParameter(5));

      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(0, 0., h_Msqr[ipt]->GetEntries());
      // if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(1, pid_Msqr.at(ipid)*0.5, pid_Msqr.at(ipid)*1.5);
      // else f2_gaus2_Msqr[ipt][ipid]->SetParLimits(1, pid_Msqr.at(ipid)*0.9, pid_Msqr.at(ipid)*1.1);
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*5.);
      else f2_gaus2_Msqr[ipt][ipid]->SetParLimits(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(3, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(4, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(5, 0., h_Msqr[ipt]->GetEntries());
      // if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(6, pid_Msqr.at(ipid)*0.5, pid_Msqr.at(ipid)*1.5);
      // else f2_gaus2_Msqr[ipt][ipid]->SetParLimits(6, pid_Msqr.at(ipid)*0.9, pid_Msqr.at(ipid)*1.1);
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(7, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*5.);
      else f2_gaus2_Msqr[ipt][ipid]->SetParLimits(7, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*0.9, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(8, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(9, f1_gaus2_Nsig[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(5)*1.5);
    }
  }
  for (int ipt = 0; ipt < NptBins; ipt++) {
    cout << "\tPt bin " << ipt <<endl;
    for (int ipid=0; ipid<Npid; ipid++) {
      // cout << "\t\tipid bin " << ipid << endl;
      for (int iter=0; iter<Niter; iter++) {
        // if (iter != Niter-1) 
          h_NsigPiMsqr[ipt]->Fit(f2_gaus2_Msqr[ipt][ipid],"MR0Q");
        // else h_NsigPiMsqr[ipt]->Fit(f2_gaus2_Msqr[ipt][ipid],"MR");
      }
    }
    timer.Print();
    timer.Continue();
  }

  // Do 1 2x2D 3gaus fits
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fitting 2x2D 3gaus of NsigmaPi and M^2..." << endl;
  TF2 *f2_gaus3_Msqr[NptBins];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    f2_gaus3_Msqr[ipt] = new TF2(Form("f2_gaus3_Msqr_pt%i", ipt), "xygaus(0)+xygaus(5) + xygaus(10)+xygaus(15) + xygaus(20)+xygaus(25)", -0.2, 1.5, -5.,30.);

    f2_gaus3_Msqr[ipt]->SetParName(0, "Peak1SglConstant");
    f2_gaus3_Msqr[ipt]->SetParName(1, "Peak1SglMeanX");
    f2_gaus3_Msqr[ipt]->SetParName(2, "Peak1SglSigmaX");
    f2_gaus3_Msqr[ipt]->SetParName(3, "Peak1SglMeanY");
    f2_gaus3_Msqr[ipt]->SetParName(4, "Peak1SglSigmaY");
    f2_gaus3_Msqr[ipt]->SetParName(5, "Peak1BgdConstant");
    f2_gaus3_Msqr[ipt]->SetParName(6, "Peak1BgdMeanX");
    f2_gaus3_Msqr[ipt]->SetParName(7, "Peak1BgdSigmaX");
    f2_gaus3_Msqr[ipt]->SetParName(8, "Peak1BgdMeanY");
    f2_gaus3_Msqr[ipt]->SetParName(9, "Peak1BgdSigmaY");
    f2_gaus3_Msqr[ipt]->SetParName(10, "Peak2SglConstant");
    f2_gaus3_Msqr[ipt]->SetParName(11, "Peak2SglMeanX");
    f2_gaus3_Msqr[ipt]->SetParName(12, "Peak2SglSigmaX");
    f2_gaus3_Msqr[ipt]->SetParName(13, "Peak2SglMeanY");
    f2_gaus3_Msqr[ipt]->SetParName(14, "Peak2SglSigmaY");
    f2_gaus3_Msqr[ipt]->SetParName(15, "Peak2BgdConstant");
    f2_gaus3_Msqr[ipt]->SetParName(16, "Peak2BgdMeanX");
    f2_gaus3_Msqr[ipt]->SetParName(17, "Peak2BgdSigmaX");
    f2_gaus3_Msqr[ipt]->SetParName(18, "Peak2BgdMeanY");
    f2_gaus3_Msqr[ipt]->SetParName(19, "Peak2BgdSigmaY");
    f2_gaus3_Msqr[ipt]->SetParName(20, "Peak3SglConstant");
    f2_gaus3_Msqr[ipt]->SetParName(21, "Peak3SglMeanX");
    f2_gaus3_Msqr[ipt]->SetParName(22, "Peak3SglSigmaX");
    f2_gaus3_Msqr[ipt]->SetParName(23, "Peak3SglMeanY");
    f2_gaus3_Msqr[ipt]->SetParName(24, "Peak3SglSigmaY");
    f2_gaus3_Msqr[ipt]->SetParName(25, "Peak3BgdConstant");
    f2_gaus3_Msqr[ipt]->SetParName(26, "Peak3BgdMeanX");
    f2_gaus3_Msqr[ipt]->SetParName(27, "Peak3BgdSigmaX");
    f2_gaus3_Msqr[ipt]->SetParName(28, "Peak3BgdMeanY");
    f2_gaus3_Msqr[ipt]->SetParName(29, "Peak3BgdSigmaY");
    
    f2_gaus3_Msqr[ipt]->SetParameter(0, f2_gaus2_Msqr[ipt][0]->GetParameter(0));
    // f2_gaus3_Msqr[ipt]->FixParameter(1, pid_Msqr.at(0)); 
    f2_gaus3_Msqr[ipt]->SetParameter(1, f2_gaus2_Msqr[ipt][0]->GetParameter(1));
    f2_gaus3_Msqr[ipt]->SetParameter(2, f2_gaus2_Msqr[ipt][0]->GetParameter(2));
    // f2_gaus3_Msqr[ipt]->FixParameter(3, 0.); 
    f2_gaus3_Msqr[ipt]->SetParameter(3, f2_gaus2_Msqr[ipt][0]->GetParameter(3));
    f2_gaus3_Msqr[ipt]->SetParameter(4, f2_gaus2_Msqr[ipt][0]->GetParameter(4));
    f2_gaus3_Msqr[ipt]->SetParameter(5, f2_gaus2_Msqr[ipt][0]->GetParameter(5));
    // f2_gaus3_Msqr[ipt]->FixParameter(6, pid_Msqr.at(0)); 
    f2_gaus3_Msqr[ipt]->SetParameter(6, f2_gaus2_Msqr[ipt][0]->GetParameter(6));
    f2_gaus3_Msqr[ipt]->SetParameter(7, f2_gaus2_Msqr[ipt][0]->GetParameter(7));
    // f2_gaus3_Msqr[ipt]->FixParameter(8, 0.); 
    f2_gaus3_Msqr[ipt]->SetParameter(8, f2_gaus2_Msqr[ipt][0]->GetParameter(8));
    f2_gaus3_Msqr[ipt]->SetParameter(9, f2_gaus2_Msqr[ipt][0]->GetParameter(9));
    f2_gaus3_Msqr[ipt]->SetParameter(10, f2_gaus2_Msqr[ipt][1]->GetParameter(0));
    // f2_gaus3_Msqr[ipt]->FixParameter(11, pid_Msqr.at(1)); 
    f2_gaus3_Msqr[ipt]->SetParameter(11, f2_gaus2_Msqr[ipt][1]->GetParameter(1));
    f2_gaus3_Msqr[ipt]->SetParameter(12, f2_gaus2_Msqr[ipt][1]->GetParameter(2));
    f2_gaus3_Msqr[ipt]->SetParameter(13, f2_gaus2_Msqr[ipt][1]->GetParameter(3));
    f2_gaus3_Msqr[ipt]->SetParameter(14, f2_gaus2_Msqr[ipt][1]->GetParameter(4));
    f2_gaus3_Msqr[ipt]->SetParameter(15, f2_gaus2_Msqr[ipt][1]->GetParameter(5));
    // f2_gaus3_Msqr[ipt]->FixParameter(16, pid_Msqr.at(1)); 
    f2_gaus3_Msqr[ipt]->SetParameter(16, f2_gaus2_Msqr[ipt][1]->GetParameter(6));
    f2_gaus3_Msqr[ipt]->SetParameter(17, f2_gaus2_Msqr[ipt][1]->GetParameter(7));
    f2_gaus3_Msqr[ipt]->SetParameter(18, f2_gaus2_Msqr[ipt][1]->GetParameter(8));
    f2_gaus3_Msqr[ipt]->SetParameter(19, f2_gaus2_Msqr[ipt][1]->GetParameter(9));
    f2_gaus3_Msqr[ipt]->SetParameter(20, f2_gaus2_Msqr[ipt][2]->GetParameter(0));
    // f2_gaus3_Msqr[ipt]->FixParameter(21, pid_Msqr.at(2)); 
    f2_gaus3_Msqr[ipt]->SetParameter(21, f2_gaus2_Msqr[ipt][2]->GetParameter(1));
    f2_gaus3_Msqr[ipt]->SetParameter(22, f2_gaus2_Msqr[ipt][2]->GetParameter(2));
    f2_gaus3_Msqr[ipt]->SetParameter(23, f2_gaus2_Msqr[ipt][2]->GetParameter(3));
    f2_gaus3_Msqr[ipt]->SetParameter(24, f2_gaus2_Msqr[ipt][2]->GetParameter(4));
    f2_gaus3_Msqr[ipt]->SetParameter(25, f2_gaus2_Msqr[ipt][2]->GetParameter(5));
    // f2_gaus3_Msqr[ipt]->FixParameter(26, pid_Msqr.at(2)); 
    f2_gaus3_Msqr[ipt]->SetParameter(26, f2_gaus2_Msqr[ipt][2]->GetParameter(6));
    f2_gaus3_Msqr[ipt]->SetParameter(27, f2_gaus2_Msqr[ipt][2]->GetParameter(7));
    f2_gaus3_Msqr[ipt]->SetParameter(28, f2_gaus2_Msqr[ipt][2]->GetParameter(8));
    f2_gaus3_Msqr[ipt]->SetParameter(29, f2_gaus2_Msqr[ipt][2]->GetParameter(9));

    f2_gaus3_Msqr[ipt]->SetParLimits(0, f2_gaus2_Msqr[ipt][0]->GetParameter(0)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(0)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(1, f2_gaus2_Msqr[ipt][0]->GetParameter(1)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(1)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(2, f2_gaus2_Msqr[ipt][0]->GetParameter(2)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(2)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(3, f2_gaus2_Msqr[ipt][0]->GetParameter(3)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(3)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(4, f2_gaus2_Msqr[ipt][0]->GetParameter(4)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(4)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(5, f2_gaus2_Msqr[ipt][0]->GetParameter(5)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(5)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(6, f2_gaus2_Msqr[ipt][0]->GetParameter(6)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(6)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(7, f2_gaus2_Msqr[ipt][0]->GetParameter(7)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(7)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(8, f2_gaus2_Msqr[ipt][0]->GetParameter(8)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(8)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(9, f2_gaus2_Msqr[ipt][0]->GetParameter(9)*0.5, f2_gaus2_Msqr[ipt][0]->GetParameter(9)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(10, f2_gaus2_Msqr[ipt][1]->GetParameter(0)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(0)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(11, f2_gaus2_Msqr[ipt][1]->GetParameter(1)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(1)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(12, f2_gaus2_Msqr[ipt][1]->GetParameter(2)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(2)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(13, f2_gaus2_Msqr[ipt][1]->GetParameter(3)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(3)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(14, f2_gaus2_Msqr[ipt][1]->GetParameter(4)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(4)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(15, f2_gaus2_Msqr[ipt][1]->GetParameter(5)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(5)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(16, f2_gaus2_Msqr[ipt][1]->GetParameter(6)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(6)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(17, f2_gaus2_Msqr[ipt][1]->GetParameter(7)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(7)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(18, f2_gaus2_Msqr[ipt][1]->GetParameter(8)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(8)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(19, f2_gaus2_Msqr[ipt][1]->GetParameter(9)*0.5, f2_gaus2_Msqr[ipt][1]->GetParameter(9)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(20, f2_gaus2_Msqr[ipt][2]->GetParameter(0)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(0)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(21, f2_gaus2_Msqr[ipt][2]->GetParameter(1)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(1)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(22, f2_gaus2_Msqr[ipt][2]->GetParameter(2)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(2)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(23, f2_gaus2_Msqr[ipt][2]->GetParameter(3)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(3)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(24, f2_gaus2_Msqr[ipt][2]->GetParameter(4)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(4)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(25, f2_gaus2_Msqr[ipt][2]->GetParameter(5)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(5)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(26, f2_gaus2_Msqr[ipt][2]->GetParameter(6)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(6)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(27, f2_gaus2_Msqr[ipt][2]->GetParameter(7)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(7)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(28, f2_gaus2_Msqr[ipt][2]->GetParameter(8)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(8)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(29, f2_gaus2_Msqr[ipt][2]->GetParameter(9)*0.5, f2_gaus2_Msqr[ipt][2]->GetParameter(9)*1.5);
  }
  // for (int ipt = 0; ipt < NptBins; ipt++) {
  //   cout << "\tPt bin " << ipt <<endl;
  //   for (int iter=0; iter<Niter; iter++) {
  //     // if (iter != Niter-1) 
  //       h_NsigPiMsqr[ipt]->Fit(f2_gaus3_Msqr[ipt],"MR0Q");
  //     // else h_NsigPiMsqr[ipt]->Fit(f2_gaus3_Msqr[ipt],"MR");
  //   }
  //   timer.Print();
  //   timer.Continue();
  // }

  // Writing output
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Writing output..." << endl;
  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_Msqr[ipt]->Write();
    // f1_gaus3_Msqr[ipt]->Write();
    for (int ipid=0; ipid<Npid; ipid++) {
      h_NsigPi[ipt][ipid]->Write();
      f1_gaus1_Msqr[ipt][ipid]->Write();
      f1_gaus1_Nsig[ipt][ipid]->Write();
      f1_gaus2_Msqr[ipt][ipid]->Write();
      f1_gaus2_Nsig[ipt][ipid]->Write();
    }
  }
  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_NsigPiMsqr[ipt]->Write();
    f2_gaus3_Msqr[ipt]->Write();
    for (int ipid=0; ipid<Npid; ipid++) {
      // f2_gaus1_Msqr[ipt][ipid]->Write();
      f2_gaus2_Msqr[ipt][ipid]->Write();
    }
  }

  timer.Stop();
  timer.Print();
}