#include <iostream>
#include <vector>

#include <TStopwatch.h>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <ROOT/TProcessExecutor.hxx>
#include <ROOT/TSeq.hxx>

using std::string;
using std::vector;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;

pair<double, double> GetXY(double _val_m2, double _val_sig, 
  double _fit_mean_m2_pi, double _fit_mean_m2_ka,
  double _fit_mean_sig_pi, double _fit_mean_sig_ka,
  double _fit_sig_m2_pi,
  double _fit_sig_sig_pi);

void DoLaplaceSmooth(vector<double> &vect, int istep, vector<double> weight = {});
TGraphErrors *DoLaplaceSmooth(TGraphErrors *const& gr, int niter);

void RunPidFits(string iFileName, string oFileName, bool is_multithread = false)
{
  TStopwatch timer;
  timer.Start();

  bool fast_check = false; //true or false;

  const int Ncentralities = 16;
  const int NptBins = 15;
  const double ptBins[NptBins+1] = {0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2};

  const vector<double> pid_mass = {0.13957, 0.493677, 0.938272}; //pi, K, p masses in GeV/c2
  const vector<double> pid_Msqr = {0.13957*0.13957, 0.493677*0.493677, 0.938272*0.938272}; //pi, K, p m^2 in (GeV/c^2)^2
  const int Npid = 3; // pi, K, p

  const int Niter1D = 5; // for more precise fits
  const int Niter2D = 1; // for more precise fits

  const vector<vector<double>> NsigMeanInitial = {{0., 10., 20.}, {0., 4., 12.}, {0., 2., 8.}, {0., 1., 5.},                     // 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.0
                                                  {0., 0., 3.}, {0., -1., 1.}, {0., -1., 0.}, {0., -1.2, -0.5}, {0., -1.5, -1.}, // 1.0-1.2, 1.2-1.4, 1.4-1.6, 1.6-1.8, 1.8-2.0
                                                  {0., -1.5, -1.5}, {0.,  -1.5, -1.5}, {0.,  -1.5, -2.}, {0.,  -1.5, -2.}, {0.,  -1.5, -2.}, // 2.0-2.2, 2.2-2.4, 2.4-2.6, 2.6-2.8, 2.8-3.0
                                                  {0., -2., -2.5}};                                                              // 3.0-3.2

  string fitOptions;
  if (is_multithread) fitOptions = "MRNQ MULTITHREAD";
  else fitOptions = "MRNQ";

  // Open aux file
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
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
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
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
      for (int iter=0; iter<Niter1D; iter++) {
        h_Msqr[ipt]->Fit(f1_gaus1_Msqr[ipt][ipid], "MRNQ");
      }
      h_NsigPi[ipt][ipid] = (TH1D*) h_NsigPiMsqr[ipt]->ProjectionY(Form("h_NsigPi_pt%i_pid%i", ipt, ipid), 
        h_NsigPiMsqr[ipt]->GetXaxis()->FindBin(f1_gaus1_Msqr[ipt][ipid]->GetParameter(1)-0.25*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2)),
        h_NsigPiMsqr[ipt]->GetXaxis()->FindBin(f1_gaus1_Msqr[ipt][ipid]->GetParameter(1)+0.25*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2)));
      f1_gaus1_Nsig[ipt][ipid] = new TF1(Form("f1_gaus1_Nsig_pt%i_pid%i", ipt, ipid), "gaus", 
        NsigMeanInitial.at(ipt).at(ipid) - 0.5, 
        NsigMeanInitial.at(ipt).at(ipid) + 0.5);
      for (int iter=0; iter<Niter1D; iter++) {
        h_NsigPi[ipt][ipid]->Fit(f1_gaus1_Nsig[ipt][ipid], "MRNQ");
      }
    }
  }

  // Do a 2x1D fits
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fitting 2x1D projections of NsigmaPi and M^2..." << endl;
  for (int ipt = 0; ipt < NptBins; ipt++) {
    f1_gaus3_Msqr[ipt] = new TF1(Form("f1_gaus1_Msqr_pt%i", ipt), "gaus(0)+gaus(3) + gaus(6)+gaus(9) + gaus(12)+gaus(15)", -0.2, 1.5);
    for (int ipid=0; ipid<Npid; ipid++) {
      if (ipid == 0) f1_gaus2_Msqr[ipt][ipid] = new TF1(Form("f1_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "gaus(0)+gaus(3)", pid_Msqr.at(ipid)*0.1, pid_Msqr.at(ipid)*2.);
      else f1_gaus2_Msqr[ipt][ipid] = new TF1(Form("f1_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "gaus(0)+gaus(3)", pid_Msqr.at(ipid)*0.7, pid_Msqr.at(ipid)*1.3);

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
      f1_gaus2_Msqr[ipt][ipid]->SetParameter(5, 10.*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2));

      f1_gaus2_Msqr[ipt][ipid]->SetParLimits(0, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0)*1.5);
      if (ipid == 0) f1_gaus2_Msqr[ipt][ipid]->SetParLimits(1, pid_Msqr.at(ipid)*0., pid_Msqr.at(ipid)*5.);
      else f1_gaus2_Msqr[ipt][ipid]->SetParLimits(1, pid_Msqr.at(ipid)*0.5, pid_Msqr.at(ipid)*1.5);
      f1_gaus2_Msqr[ipt][ipid]->SetParLimits(2, f1_gaus1_Msqr[ipt][ipid]->GetParameter(2)*0.5, f1_gaus1_Msqr[ipt][ipid]->GetParameter(2)*1.5);
      f1_gaus2_Msqr[ipt][ipid]->SetParLimits(3, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_Msqr[ipt][ipid]->GetParameter(0)*1.5);
      if (ipid == 0) f1_gaus2_Msqr[ipt][ipid]->SetParLimits(4, pid_Msqr.at(ipid)*0., pid_Msqr.at(ipid)*5.);
      else f1_gaus2_Msqr[ipt][ipid]->SetParLimits(4, pid_Msqr.at(ipid)*0.5, pid_Msqr.at(ipid)*1.5);
      f1_gaus2_Msqr[ipt][ipid]->SetParLimits(5, 1.5*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2), 100.*f1_gaus1_Msqr[ipt][ipid]->GetParameter(2));

      for (int iter=0; iter<Niter1D; iter++) {
        h_Msqr[ipt]->Fit(f1_gaus2_Msqr[ipt][ipid], "MRNQ");
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
      f1_gaus2_Nsig[ipt][ipid]->SetParameter(5, 10.*f1_gaus1_Nsig[ipt][ipid]->GetParameter(2));

      for (int iter=0; iter<Niter1D; iter++) {
        h_NsigPi[ipt][ipid]->Fit(f1_gaus2_Nsig[ipt][ipid], "MRNQ");
      }
    }
  }

  // Do 3 2x2D gaus fits
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fitting 3 2x2D of NsigmaPi and M^2..." << endl;
  TF2 *f2_gaus2_Msqr[NptBins][Npid];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    for (int ipid=0; ipid<Npid; ipid++) {
      if (ipid == 0 && ipt < 5)
        f2_gaus2_Msqr[ipt][ipid] = new TF2(Form("f2_gaus2_Msqr_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          pid_Msqr.at(ipid)-4.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          pid_Msqr.at(ipid)+4.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)-8.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2),
          f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)+8.*f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));
      if (ipid == 0 && ipt >= 5)
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
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParameter(1, pid_Msqr.at(ipid));
      else f2_gaus2_Msqr[ipt][ipid]->SetParameter(1, pid_Msqr.at(ipid));
      // f2_gaus2_Msqr[ipt][ipid]->FixParameter(1, pid_Msqr.at(ipid));
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParameter(2, 5.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(2));
      else f2_gaus2_Msqr[ipt][ipid]->SetParameter(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(3, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(4, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(3));
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParameter(6, pid_Msqr.at(ipid));
      else f2_gaus2_Msqr[ipt][ipid]->SetParameter(6, pid_Msqr.at(ipid));
      // f2_gaus2_Msqr[ipt][ipid]->FixParameter(6, pid_Msqr.at(ipid));
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParameter(7, 10.*f1_gaus2_Msqr[ipt][ipid]->GetParameter(5));
      else f2_gaus2_Msqr[ipt][ipid]->SetParameter(7, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(8, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1));
      f2_gaus2_Msqr[ipt][ipid]->SetParameter(9, f1_gaus2_Nsig[ipt][ipid]->GetParameter(5));

      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(0, 0., h_Msqr[ipt]->GetEntries());
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(1, pid_Msqr.at(ipid)*0.9, pid_Msqr.at(ipid)*1.1);
      else f2_gaus2_Msqr[ipt][ipid]->SetParLimits(1, pid_Msqr.at(ipid)*0.9, pid_Msqr.at(ipid)*1.1);
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*20.);
      else if (ipid == 1) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*5.);
      else f2_gaus2_Msqr[ipt][ipid]->SetParLimits(2, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(2)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(3, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(4, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(2)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(5, 0., h_Msqr[ipt]->GetEntries());
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(6, pid_Msqr.at(ipid)*0.5, pid_Msqr.at(ipid)*1.5);
      else f2_gaus2_Msqr[ipt][ipid]->SetParLimits(6, pid_Msqr.at(ipid)*0.5, pid_Msqr.at(ipid)*1.5);
      if (ipid == 0) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(7, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*20.);
      else if (ipid == 1) f2_gaus2_Msqr[ipt][ipid]->SetParLimits(7, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*5.);
      else f2_gaus2_Msqr[ipt][ipid]->SetParLimits(7, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_Msqr[ipt][ipid]->GetParameter(5)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(8, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(1)*1.5);
      f2_gaus2_Msqr[ipt][ipid]->SetParLimits(9, f1_gaus2_Nsig[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_Nsig[ipt][ipid]->GetParameter(5)*1.5);
    }
  }
  for (int ipt = 0; ipt < NptBins; ipt++) {
    cout << "\tPt bin " << ipt <<endl;
    for (int ipid=0; ipid<Npid; ipid++) {
      for (int iter=0; iter<Niter2D; iter++) {
        h_NsigPiMsqr[ipt]->Fit(f2_gaus2_Msqr[ipt][ipid],fitOptions.c_str());
      }
    }
    timer.Print();
    timer.Continue();
  }

  // Do polinomial fits for parameter smoothening for 3gaus
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Do polinomial fits for parameter smoothening for 3gaus..." << endl;
  TGraphErrors *gr_signal1_mean_m2[Npid];
  TGraphErrors *gr_signal1_mean_ns[Npid];
  TGraphErrors *gr_signal1_sigm_m2[Npid];
  TGraphErrors *gr_signal1_sigm_ns[Npid];
  TGraphErrors *gr_orig_signal1_mean_m2[Npid];
  TGraphErrors *gr_orig_signal1_mean_ns[Npid];
  TGraphErrors *gr_orig_signal1_sigm_m2[Npid];
  TGraphErrors *gr_orig_signal1_sigm_ns[Npid];

  TF1 *f_signal1_mean_m2[Npid];
  TF1 *f_signal1_mean_ns[Npid];
  TF1 *f_signal1_sigm_m2[Npid];
  TF1 *f_signal1_sigm_ns[Npid];

  f_signal1_mean_m2[0] = new TF1(Form("f_signal1_mean_m2_pid%i", 0), "pol2", 0.3, 3.1);
  f_signal1_mean_m2[1] = new TF1(Form("f_signal1_mean_m2_pid%i", 1), "pol2", 0.3, 3.1);
  f_signal1_mean_m2[2] = new TF1(Form("f_signal1_mean_m2_pid%i", 2), "pol2", 0.3, 3.1);

  f_signal1_mean_ns[0] = new TF1(Form("f_signal1_mean_ns_pid%i", 0), "pol2", 0.3, 3.1);
  f_signal1_mean_ns[1] = new TF1(Form("f_signal1_mean_ns_pid%i", 1), "pol2", 0.3, 3.1);
  f_signal1_mean_ns[2] = new TF1(Form("f_signal1_mean_ns_pid%i", 2), "pol2", 0.3, 3.1);

  f_signal1_sigm_m2[0] = new TF1(Form("f_signal1_sigm_m2_pid%i", 0), "pol3", 0.3, 3.1);
  f_signal1_sigm_m2[1] = new TF1(Form("f_signal1_sigm_m2_pid%i", 1), "pol2", 0.3, 3.1);
  f_signal1_sigm_m2[2] = new TF1(Form("f_signal1_sigm_m2_pid%i", 2), "pol5", 0.3, 3.1);

  f_signal1_sigm_ns[0] = new TF1(Form("f_signal1_sigm_ns_pid%i", 0), "pol2", 0.3, 3.1);
  f_signal1_sigm_ns[1] = new TF1(Form("f_signal1_sigm_ns_pid%i", 1), "pol2", 0.3, 3.1);
  f_signal1_sigm_ns[2] = new TF1(Form("f_signal1_sigm_ns_pid%i", 2), "pol2", 0.3, 3.1);

  vector<double> vx_tmp, vy_tmp, vex_tmp, vey_tmp;
  for (int ipid=0; ipid<Npid; ipid++) {
    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus2_Msqr[ipt][ipid]->GetParameter(1));
      vey_tmp.push_back(f2_gaus2_Msqr[ipt][ipid]->GetParError(1));
    }
    gr_orig_signal1_mean_m2[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal1_mean_m2[ipid]->SetName(Form("gr_orig_signal1_mean_m2_pid%i", ipid));
    gr_orig_signal1_mean_m2[ipid]->SetTitle(Form("gr_orig_signal1_mean_m2_pid%i", ipid));
    gr_orig_signal1_mean_m2[ipid]->SetMarkerStyle(kFullSquare);
    gr_orig_signal1_mean_m2[ipid]->Fit(f_signal1_mean_m2[ipid], "RNQ");

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus2_Msqr[ipt][ipid]->GetParameter(3));
      vey_tmp.push_back(f2_gaus2_Msqr[ipt][ipid]->GetParError(3));
    }
    gr_orig_signal1_mean_ns[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal1_mean_ns[ipid]->SetName(Form("gr_orig_signal1_mean_ns_pid%i", ipid));
    gr_orig_signal1_mean_ns[ipid]->SetTitle(Form("gr_orig_signal1_mean_ns_pid%i", ipid));
    gr_orig_signal1_mean_ns[ipid]->SetMarkerStyle(kFullSquare);
    gr_orig_signal1_mean_ns[ipid]->Fit(f_signal1_mean_ns[ipid], "RNQ");

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus2_Msqr[ipt][ipid]->GetParameter(2));
      vey_tmp.push_back(f2_gaus2_Msqr[ipt][ipid]->GetParError(2));
    }
    gr_orig_signal1_sigm_m2[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal1_sigm_m2[ipid]->SetName(Form("gr_orig_signal1_sigm_m2_pid%i", ipid));
    gr_orig_signal1_sigm_m2[ipid]->SetTitle(Form("gr_orig_signal1_sigm_m2_pid%i", ipid));
    gr_orig_signal1_sigm_m2[ipid]->SetMarkerStyle(kFullSquare);
    gr_orig_signal1_sigm_m2[ipid]->Fit(f_signal1_sigm_m2[ipid], "RNQ");

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus2_Msqr[ipt][ipid]->GetParameter(4));
      vey_tmp.push_back(f2_gaus2_Msqr[ipt][ipid]->GetParError(4));
    }
    gr_orig_signal1_sigm_ns[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal1_sigm_ns[ipid]->SetName(Form("gr_orig_signal1_sigm_ns_pid%i", ipid));
    gr_orig_signal1_sigm_ns[ipid]->SetTitle(Form("gr_orig_signal1_sigm_ns_pid%i", ipid));
    gr_orig_signal1_sigm_ns[ipid]->SetMarkerStyle(kFullSquare);
    gr_orig_signal1_sigm_ns[ipid]->Fit(f_signal1_sigm_ns[ipid], "RNQ");
  }

  // Smoothen the pT-dependencies
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_signal1_mean_m2[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal1_mean_m2[ipid], -1);
    gr_signal1_mean_ns[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal1_mean_ns[ipid], -1);
    gr_signal1_sigm_m2[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal1_sigm_m2[ipid], -1);
    gr_signal1_sigm_ns[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal1_sigm_ns[ipid], -1);
  }

  // Do 1 2x2D 3gaus fits
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fitting 2x2D 3gaus of NsigmaPi and M^2..." << endl;
  TF2 *f2_gaus3_Msqr[NptBins];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    f2_gaus3_Msqr[ipt] = new TF2(Form("f2_gaus3_Msqr_pt%i", ipt), "xygaus(0)+xygaus(5) + xygaus(10)+xygaus(15) + xygaus(20)+xygaus(25)", 0., 1.2, -5.,25.);

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
    f2_gaus3_Msqr[ipt]->FixParameter(1, pid_Msqr.at(0)); 
    // f2_gaus3_Msqr[ipt]->SetParameter(1, f_signal1_mean_m2[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(2, gr_signal1_sigm_m2[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->FixParameter(3, 0.); 
    // f2_gaus3_Msqr[ipt]->SetParameter(3, gr_signal1_mean_ns[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(4, gr_signal1_sigm_ns[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(5, f2_gaus2_Msqr[ipt][0]->GetParameter(5));
    f2_gaus3_Msqr[ipt]->FixParameter(6, pid_Msqr.at(0)); 
    // f2_gaus3_Msqr[ipt]->SetParameter(6, gr_signal1_mean_m2[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(7, f2_gaus2_Msqr[ipt][0]->GetParameter(7));
    f2_gaus3_Msqr[ipt]->FixParameter(8, 0.); 
    // f2_gaus3_Msqr[ipt]->SetParameter(8, f2_gaus2_Msqr[ipt][0]->GetParameter(8));
    f2_gaus3_Msqr[ipt]->SetParameter(9, f2_gaus2_Msqr[ipt][0]->GetParameter(9));
    f2_gaus3_Msqr[ipt]->SetParameter(10, f2_gaus2_Msqr[ipt][1]->GetParameter(0));
    f2_gaus3_Msqr[ipt]->FixParameter(11, pid_Msqr.at(1)); 
    // f2_gaus3_Msqr[ipt]->SetParameter(11, gr_signal1_mean_m2[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(12, gr_signal1_sigm_m2[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(13, gr_signal1_mean_ns[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(14, gr_signal1_sigm_ns[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(15, f2_gaus2_Msqr[ipt][1]->GetParameter(5));
    f2_gaus3_Msqr[ipt]->FixParameter(16, pid_Msqr.at(1)); 
    // f2_gaus3_Msqr[ipt]->SetParameter(16, gr_signal1_mean_m2[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(17, f2_gaus2_Msqr[ipt][1]->GetParameter(7));
    f2_gaus3_Msqr[ipt]->SetParameter(18, gr_signal1_mean_ns[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(19, f2_gaus2_Msqr[ipt][1]->GetParameter(9));
    f2_gaus3_Msqr[ipt]->SetParameter(20, f2_gaus2_Msqr[ipt][2]->GetParameter(0));
    f2_gaus3_Msqr[ipt]->FixParameter(21, pid_Msqr.at(2)); 
    // f2_gaus3_Msqr[ipt]->SetParameter(21, gr_signal1_mean_m2[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(22, gr_signal1_sigm_m2[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(23, f2_gaus2_Msqr[ipt][2]->GetParameter(3));
    f2_gaus3_Msqr[ipt]->SetParameter(24, f2_gaus2_Msqr[ipt][2]->GetParameter(4));
    f2_gaus3_Msqr[ipt]->SetParameter(25, f2_gaus2_Msqr[ipt][2]->GetParameter(5));
    f2_gaus3_Msqr[ipt]->FixParameter(26, pid_Msqr.at(2)); 
    // f2_gaus3_Msqr[ipt]->SetParameter(26, gr_signal1_mean_m2[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_Msqr[ipt]->SetParameter(27, f2_gaus2_Msqr[ipt][2]->GetParameter(7));
    f2_gaus3_Msqr[ipt]->SetParameter(28, f2_gaus2_Msqr[ipt][2]->GetParameter(8));
    f2_gaus3_Msqr[ipt]->SetParameter(29, f2_gaus2_Msqr[ipt][2]->GetParameter(9));

    f2_gaus3_Msqr[ipt]->SetParLimits(0, f2_gaus3_Msqr[ipt]->GetParameter(0)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(0)*1.5);
    // f2_gaus3_Msqr[ipt]->SetParLimits(1, f2_gaus3_Msqr[ipt]->GetParameter(1)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(1)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(2, f2_gaus3_Msqr[ipt]->GetParameter(2)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(2)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(3, f2_gaus3_Msqr[ipt]->GetParameter(3)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(3)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(4, f2_gaus3_Msqr[ipt]->GetParameter(4)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(4)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(5, f2_gaus3_Msqr[ipt]->GetParameter(5)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(5)*1.5);
    // f2_gaus3_Msqr[ipt]->SetParLimits(6, f2_gaus3_Msqr[ipt]->GetParameter(6)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(6)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(7, f2_gaus3_Msqr[ipt]->GetParameter(7)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(7)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(8, f2_gaus3_Msqr[ipt]->GetParameter(8)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(8)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(9, f2_gaus3_Msqr[ipt]->GetParameter(9)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(9)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(10, f2_gaus3_Msqr[ipt]->GetParameter(10)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(10)*1.5);
    // f2_gaus3_Msqr[ipt]->SetParLimits(11, f2_gaus3_Msqr[ipt]->GetParameter(11)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(11)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(12, f2_gaus3_Msqr[ipt]->GetParameter(12)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(12)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(13, f2_gaus3_Msqr[ipt]->GetParameter(13)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(13)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(14, f2_gaus3_Msqr[ipt]->GetParameter(14)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(14)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(15, f2_gaus3_Msqr[ipt]->GetParameter(15)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(15)*1.5);
    // f2_gaus3_Msqr[ipt]->SetParLimits(16, f2_gaus3_Msqr[ipt]->GetParameter(16)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(16)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(17, f2_gaus3_Msqr[ipt]->GetParameter(17)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(17)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(18, f2_gaus3_Msqr[ipt]->GetParameter(18)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(18)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(19, f2_gaus3_Msqr[ipt]->GetParameter(19)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(19)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(20, f2_gaus3_Msqr[ipt]->GetParameter(20)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(20)*1.5);
    // f2_gaus3_Msqr[ipt]->SetParLimits(21, f2_gaus3_Msqr[ipt]->GetParameter(21)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(21)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(22, f2_gaus3_Msqr[ipt]->GetParameter(22)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(22)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(23, f2_gaus3_Msqr[ipt]->GetParameter(23)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(23)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(24, f2_gaus3_Msqr[ipt]->GetParameter(24)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(24)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(25, f2_gaus3_Msqr[ipt]->GetParameter(25)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(25)*1.5);
    // f2_gaus3_Msqr[ipt]->SetParLimits(26, f2_gaus3_Msqr[ipt]->GetParameter(26)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(26)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(27, f2_gaus3_Msqr[ipt]->GetParameter(27)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(27)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(28, f2_gaus3_Msqr[ipt]->GetParameter(28)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(28)*1.5);
    f2_gaus3_Msqr[ipt]->SetParLimits(29, f2_gaus3_Msqr[ipt]->GetParameter(29)*0.5, f2_gaus3_Msqr[ipt]->GetParameter(29)*1.5);
  }
  for (int ipt = 0; ipt < NptBins; ipt++) {
    cout << "\tPt bin " << ipt <<endl;
    for (int iter=0; iter<Niter2D; iter++) {
        if (!fast_check) h_NsigPiMsqr[ipt]->Fit(f2_gaus3_Msqr[ipt],Form("%s B",fitOptions.c_str()));
    }
    timer.Print();
    timer.Continue();
  }

  // Do polinomial fits for x,y transformation
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Do polinomial fits for x,y transformation..." << endl;
  TGraphErrors *gr_signal2_mean_m2[Npid];
  TGraphErrors *gr_signal2_mean_ns[Npid];
  TGraphErrors *gr_signal2_sigm_m2[Npid];
  TGraphErrors *gr_signal2_sigm_ns[Npid];
  TGraphErrors *gr_orig_signal2_mean_m2[Npid];
  TGraphErrors *gr_orig_signal2_mean_ns[Npid];
  TGraphErrors *gr_orig_signal2_sigm_m2[Npid];
  TGraphErrors *gr_orig_signal2_sigm_ns[Npid];

  TF1 *f_signal2_mean_m2[Npid];
  TF1 *f_signal2_mean_ns[Npid];
  TF1 *f_signal2_sigm_m2[Npid];
  TF1 *f_signal2_sigm_ns[Npid];

  f_signal2_mean_m2[0] = new TF1(Form("f_signal2_mean_m2_pid%i", 0), "pol2", 0.3, 3.1);
  f_signal2_mean_m2[1] = new TF1(Form("f_signal2_mean_m2_pid%i", 1), "pol2", 0.3, 3.1);
  f_signal2_mean_m2[2] = new TF1(Form("f_signal2_mean_m2_pid%i", 2), "pol2", 0.3, 3.1);

  f_signal2_mean_ns[0] = new TF1(Form("f_signal2_mean_ns_pid%i", 0), "pol2", 0.3, 3.1);
  f_signal2_mean_ns[1] = new TF1(Form("f_signal2_mean_ns_pid%i", 1), "pol2", 0.3, 3.1);
  f_signal2_mean_ns[2] = new TF1(Form("f_signal2_mean_ns_pid%i", 2), "pol2", 0.3, 3.1);

  f_signal2_sigm_m2[0] = new TF1(Form("f_signal2_sigm_m2_pid%i", 0), "pol2", 0.3, 3.1);
  f_signal2_sigm_m2[1] = new TF1(Form("f_signal2_sigm_m2_pid%i", 1), "pol2", 0.3, 3.1);
  f_signal2_sigm_m2[2] = new TF1(Form("f_signal2_sigm_m2_pid%i", 2), "pol2", 0.3, 3.1);

  f_signal2_sigm_ns[0] = new TF1(Form("f_signal2_sigm_ns_pid%i", 0), "pol2", 0.3, 3.1);
  f_signal2_sigm_ns[1] = new TF1(Form("f_signal2_sigm_ns_pid%i", 1), "pol2", 0.3, 3.1);
  f_signal2_sigm_ns[2] = new TF1(Form("f_signal2_sigm_ns_pid%i", 2), "pol2", 0.3, 3.1);

  for (int ipid=0; ipid<Npid; ipid++) {
    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus3_Msqr[ipt]->GetParameter(1+10*ipid));
      vey_tmp.push_back(f2_gaus3_Msqr[ipt]->GetParError(1+10*ipid));
    }
    gr_orig_signal2_mean_m2[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal2_mean_m2[ipid]->SetName(Form("gr_orig_signal2_mean_m2_pid%i", ipid));
    gr_orig_signal2_mean_m2[ipid]->SetTitle(Form("gr_orig_signal2_mean_m2_pid%i", ipid));
    gr_orig_signal2_mean_m2[ipid]->SetMarkerStyle(kFullSquare);
    gr_orig_signal2_mean_m2[ipid]->Fit(f_signal2_mean_m2[ipid], "RNQ");

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus3_Msqr[ipt]->GetParameter(3+10*ipid));
      vey_tmp.push_back(f2_gaus3_Msqr[ipt]->GetParError(3+10*ipid));
    }
    gr_orig_signal2_mean_ns[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal2_mean_ns[ipid]->SetName(Form("gr_orig_signal2_mean_ns_pid%i", ipid));
    gr_orig_signal2_mean_ns[ipid]->SetTitle(Form("gr_orig_signal2_mean_ns_pid%i", ipid));
    gr_orig_signal2_mean_ns[ipid]->SetMarkerStyle(kFullSquare);
    gr_orig_signal2_mean_ns[ipid]->Fit(f_signal2_mean_ns[ipid], "RNQ");

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus3_Msqr[ipt]->GetParameter(2+10*ipid));
      vey_tmp.push_back(f2_gaus3_Msqr[ipt]->GetParError(2+10*ipid));
    }
    gr_orig_signal2_sigm_m2[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal2_sigm_m2[ipid]->SetName(Form("gr_orig_signal2_sigm_m2_pid%i", ipid));
    gr_orig_signal2_sigm_m2[ipid]->SetTitle(Form("gr_orig_signal2_sigm_m2_pid%i", ipid));
    gr_orig_signal2_sigm_m2[ipid]->SetMarkerStyle(kFullSquare);
    gr_orig_signal2_sigm_m2[ipid]->Fit(f_signal2_sigm_m2[ipid], "RNQ");

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus3_Msqr[ipt]->GetParameter(4+10*ipid));
      vey_tmp.push_back(f2_gaus3_Msqr[ipt]->GetParError(4+10*ipid));
    }
    gr_orig_signal2_sigm_ns[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal2_sigm_ns[ipid]->SetName(Form("gr_orig_signal2_sigm_ns_pid%i", ipid));
    gr_orig_signal2_sigm_ns[ipid]->SetTitle(Form("gr_orig_signal2_sigm_ns_pid%i", ipid));
    gr_orig_signal2_sigm_ns[ipid]->SetMarkerStyle(kFullSquare);
    gr_orig_signal2_sigm_ns[ipid]->Fit(f_signal2_sigm_ns[ipid], "RNQ");
  }

  // Smoothen the pT-dependencies
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_signal2_mean_m2[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal2_mean_m2[ipid], -1);
    gr_signal2_mean_ns[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal2_mean_ns[ipid], -1);
    gr_signal2_sigm_m2[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal2_sigm_m2[ipid], -1);
    gr_signal2_sigm_ns[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal2_sigm_ns[ipid], -1);
  }
  // Resampling and transforming m2,nsigma -> x,y
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Resampling and doing (m2,nsig)->(x,y) transformation..." << endl;
  TH2D *h_pidXY[NptBins];
  TH2D *h_pikaXY[NptBins];
  double x_old, y_old;
  std::pair<double,double> new_coord;
  for (int ipt = 0; ipt < NptBins; ipt++) {
    cout << "\tPt bin " << ipt <<endl;

    h_pidXY[ipt] = new TH2D(Form("h_pidXY_pt%i", ipt), Form("h_pidXY for %1.1f < p_{T} < %1.1f GeV/c;X;Y", ptBins[ipt], ptBins[ipt+1]), 6000, -1.5, 1.5, 6000, -1.5, 1.5);
    h_pikaXY[ipt] = new TH2D(Form("h_pikaXY_pt%i", ipt), Form("h_pikaXY for %1.1f < p_{T} < %1.1f GeV/c (#pi, K only);X;Y", ptBins[ipt], ptBins[ipt+1]), 6000, -1.5, 1.5, 6000, -1.5, 1.5);

    long nentries = (long)h_NsigPiMsqr[ipt]->GetEntries()/(long)(exp((double)(NptBins-ipt+1)/3.));
    for (long ientry=0; ientry<nentries; ientry++) {
      h_NsigPiMsqr[ipt]->GetRandom2(x_old,y_old);
      new_coord = GetXY(x_old, y_old, 
        pid_Msqr.at(0), pid_Msqr.at(1),
        // gr_signal2_mean_m2[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.), gr_signal2_mean_m2[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.),
        gr_signal2_mean_ns[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.), gr_signal2_mean_ns[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.),
        gr_signal2_sigm_m2[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.), gr_signal2_sigm_ns[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
      
      h_pidXY[ipt]->Fill(new_coord.first, new_coord.second);
      if (x_old < 0.65) h_pikaXY[ipt]->Fill(new_coord.first, new_coord.second);
    }

    timer.Print();
    timer.Continue();
  }

  // Fit (x,y): 3 1D gaus
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fit (x,y): 3 1D gaus..." << endl;
  TH1D *h_x[NptBins];
  TH1D *h_y[NptBins][Npid]; // slice Y for each particle

  TF1 *f1_gaus1_x[NptBins][Npid];
  TF1 *f1_gaus2_x[NptBins][Npid];
  TF1 *f1_gaus1_y[NptBins][Npid];
  TF1 *f1_gaus2_y[NptBins][Npid];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    for (int ipid=0; ipid<Npid; ipid++) {
      if (ipid == 0) f1_gaus1_x[ipt][ipid] = new TF1(Form("f1_gaus1_x_pt%i_pid%i", ipt, ipid), "gaus", -0.1, 0.1);
      else if (ipid == 1) f1_gaus1_x[ipt][ipid] = new TF1(Form("f1_gaus1_x_pt%i_pid%i", ipt, ipid), "gaus", 0.15, 0.35);
      else f1_gaus1_x[ipt][ipid] = new TF1(Form("f1_gaus1_x_pt%i_pid%i", ipt, ipid), "gaus", 0.7, 1.3);
      f1_gaus1_x[ipt][ipid]->SetParameter(1, 0.);
    }
  }

  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_x[ipt] = (TH1D*) h_pidXY[ipt]->ProjectionX(Form("h_x_pt%i", ipt), 0, -1);
    for (int ipid=0; ipid<Npid; ipid++) {
      for (int iter=0; iter<Niter1D; iter++) {
        h_x[ipt]->Fit(f1_gaus1_x[ipt][ipid], "MRNQ");
      }
      h_y[ipt][ipid] = (TH1D*) h_pidXY[ipt]->ProjectionY(Form("h_y_pt%i_pid%i", ipt, ipid), 
        h_pidXY[ipt]->GetXaxis()->FindBin(f1_gaus1_x[ipt][ipid]->GetParameter(1)-0.25*f1_gaus1_x[ipt][ipid]->GetParameter(2)),
        h_pidXY[ipt]->GetXaxis()->FindBin(f1_gaus1_x[ipt][ipid]->GetParameter(1)+0.25*f1_gaus1_x[ipt][ipid]->GetParameter(2)));
      f1_gaus1_y[ipt][ipid] = new TF1(Form("f1_gaus1_y_pt%i_pid%i", ipt, ipid), "gaus", 
        h_y[ipt][ipid]->GetMean() - 0.1, 
        h_y[ipt][ipid]->GetMean() + 0.1);
      for (int iter=0; iter<Niter1D; iter++) {
        h_y[ipt][ipid]->Fit(f1_gaus1_y[ipt][ipid], "MRNQ");
      }
    }
  }

  // Fit (x,y): 3 2x1D gaus
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fit (x,y): 3 2x1D gaus..." << endl;
  for (int ipt = 0; ipt < NptBins; ipt++) {
    for (int ipid=0; ipid<Npid; ipid++) {
      if (ipid == 0) f1_gaus2_x[ipt][ipid] = new TF1(Form("f1_gaus2_x_pt%i_pid%i", ipt, ipid), "gaus(0)+gaus(3)", -0.3, 0.3);
      else if (ipid == 1) f1_gaus2_x[ipt][ipid] = new TF1(Form("f1_gaus2_x_pt%i_pid%i", ipt, ipid), "gaus(0)+gaus(3)", 0.1, 0.4);
      else f1_gaus2_x[ipt][ipid] = new TF1(Form("f1_gaus2_x_pt%i_pid%i", ipt, ipid), "gaus(0)+gaus(3)", 0.5, 1.6);

      f1_gaus2_x[ipt][ipid]->SetParName(0, "SglConstant");
      f1_gaus2_x[ipt][ipid]->SetParName(1, "SglMean");
      f1_gaus2_x[ipt][ipid]->SetParName(2, "SglSigma");
      f1_gaus2_x[ipt][ipid]->SetParName(3, "BgdConstant");
      f1_gaus2_x[ipt][ipid]->SetParName(4, "BgdMean");
      f1_gaus2_x[ipt][ipid]->SetParName(5, "BgdSigma");

      f1_gaus2_x[ipt][ipid]->SetParameter(0, f1_gaus1_x[ipt][ipid]->GetParameter(0));
      if (ipid == 0) f1_gaus2_x[ipt][ipid]->FixParameter(1, 0.);
      else f1_gaus2_x[ipt][ipid]->SetParameter(1, f1_gaus1_x[ipt][ipid]->GetParameter(1));
      // f1_gaus2_x[ipt][ipid]->SetParameter(1, f1_gaus1_x[ipt][ipid]->GetParameter(1));
      f1_gaus2_x[ipt][ipid]->SetParameter(2, f1_gaus1_x[ipt][ipid]->GetParameter(2));
      f1_gaus2_x[ipt][ipid]->SetParameter(3, f1_gaus1_x[ipt][ipid]->GetParameter(0));
      if (ipid == 0) f1_gaus2_x[ipt][ipid]->FixParameter(4, 0.);
      else f1_gaus2_x[ipt][ipid]->SetParameter(4, f1_gaus1_x[ipt][ipid]->GetParameter(1));
      f1_gaus2_x[ipt][ipid]->SetParameter(5, 10.*f1_gaus1_x[ipt][ipid]->GetParameter(2));

      f1_gaus2_x[ipt][ipid]->SetParLimits(0, f1_gaus1_x[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_x[ipt][ipid]->GetParameter(0)*1.5);
      if (ipid == 0) f1_gaus2_x[ipt][ipid]->SetParLimits(1, -0.2, 0.2);
      else if (ipid == 1) f1_gaus2_x[ipt][ipid]->SetParLimits(1, f1_gaus1_x[ipt][ipid]->GetParameter(1)*0.95, f1_gaus1_x[ipt][ipid]->GetParameter(1)*1.05);
      else f1_gaus2_x[ipt][ipid]->SetParLimits(1, f1_gaus1_x[ipt][ipid]->GetParameter(1)*0.5, f1_gaus1_x[ipt][ipid]->GetParameter(1)*1.5);
      f1_gaus2_x[ipt][ipid]->SetParLimits(2, f1_gaus1_x[ipt][ipid]->GetParameter(2)*0.5, f1_gaus1_x[ipt][ipid]->GetParameter(2)*1.5);
      f1_gaus2_x[ipt][ipid]->SetParLimits(3, f1_gaus1_x[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_x[ipt][ipid]->GetParameter(0)*1.5);
      if (ipid == 0) f1_gaus2_x[ipt][ipid]->SetParLimits(4, -0.2, 0.2);
      else if (ipid == 1) f1_gaus2_x[ipt][ipid]->SetParLimits(4, f1_gaus1_x[ipt][ipid]->GetParameter(1)*0.9, f1_gaus1_x[ipt][ipid]->GetParameter(1)*1.1);
      else f1_gaus2_x[ipt][ipid]->SetParLimits(4, f1_gaus1_x[ipt][ipid]->GetParameter(1)*0.5, f1_gaus1_x[ipt][ipid]->GetParameter(1)*1.5);
      f1_gaus2_x[ipt][ipid]->SetParLimits(5, 1.5*f1_gaus1_x[ipt][ipid]->GetParameter(2), 100.*f1_gaus1_x[ipt][ipid]->GetParameter(2));

      for (int iter=0; iter<Niter1D; iter++) {
        h_x[ipt]->Fit(f1_gaus2_x[ipt][ipid], "MRNQ");
      }

      f1_gaus2_y[ipt][ipid] = new TF1(Form("f1_gaus2_y_pt%i_pid%i", ipt, ipid), "gaus(0)+gaus(3)", 
        h_y[ipt][ipid]->GetMean() - 0.3, 
        h_y[ipt][ipid]->GetMean() + 0.3);
      
      f1_gaus2_y[ipt][ipid]->SetParName(0, "SglConstant");
      f1_gaus2_y[ipt][ipid]->SetParName(1, "SglMean");
      f1_gaus2_y[ipt][ipid]->SetParName(2, "SglSigma");
      f1_gaus2_y[ipt][ipid]->SetParName(3, "BgdConstant");
      f1_gaus2_y[ipt][ipid]->SetParName(4, "BgdMean");
      f1_gaus2_y[ipt][ipid]->SetParName(5, "BgdSigma");

      f1_gaus2_y[ipt][ipid]->SetParameter(0, f1_gaus1_y[ipt][ipid]->GetParameter(0));
      if (ipid == 0) f1_gaus2_y[ipt][ipid]->FixParameter(1, 0.);
      else if (ipid == 1) f1_gaus2_y[ipt][ipid]->FixParameter(1, 0.);
      else f1_gaus2_y[ipt][ipid]->SetParameter(1, f1_gaus1_y[ipt][ipid]->GetParameter(1));
      f1_gaus2_y[ipt][ipid]->SetParameter(2, f1_gaus1_y[ipt][ipid]->GetParameter(2));
      f1_gaus2_y[ipt][ipid]->SetParameter(3, f1_gaus1_y[ipt][ipid]->GetParameter(0));
      if (ipid == 0) f1_gaus2_y[ipt][ipid]->FixParameter(4, 0.);
      else if (ipid == 1) f1_gaus2_y[ipt][ipid]->FixParameter(4, 0.);
      else f1_gaus2_y[ipt][ipid]->SetParameter(4, f1_gaus1_y[ipt][ipid]->GetParameter(1));
      f1_gaus2_y[ipt][ipid]->SetParameter(5, 50.*f1_gaus1_y[ipt][ipid]->GetParameter(2));

      f1_gaus2_y[ipt][ipid]->SetParLimits(0, f1_gaus1_y[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_y[ipt][ipid]->GetParameter(0)*1.5);
      if (ipid != 2) f1_gaus2_y[ipt][ipid]->SetParLimits(1, -0.005, 0.005);
      else f1_gaus2_y[ipt][ipid]->SetParLimits(1, f1_gaus1_y[ipt][ipid]->GetParameter(1)*0.5, f1_gaus1_y[ipt][ipid]->GetParameter(1)*1.5);
      f1_gaus2_y[ipt][ipid]->SetParLimits(2, f1_gaus1_y[ipt][ipid]->GetParameter(2)*0.5, f1_gaus1_y[ipt][ipid]->GetParameter(2)*1.5);
      f1_gaus2_y[ipt][ipid]->SetParLimits(3, f1_gaus1_y[ipt][ipid]->GetParameter(0)*0.5, f1_gaus1_y[ipt][ipid]->GetParameter(0)*1.5);
      if (ipid != 2) f1_gaus2_y[ipt][ipid]->SetParLimits(4, -0.005, 0.005);
      else f1_gaus2_y[ipt][ipid]->SetParLimits(4, f1_gaus1_y[ipt][ipid]->GetParameter(1)*0.5, f1_gaus1_y[ipt][ipid]->GetParameter(1)*1.5);
      f1_gaus2_y[ipt][ipid]->SetParLimits(4, f1_gaus1_y[ipt][ipid]->GetParameter(2)*1.5, f1_gaus1_y[ipt][ipid]->GetParameter(2)*100.);

      for (int iter=0; iter<Niter1D; iter++) {
        h_y[ipt][ipid]->Fit(f1_gaus2_y[ipt][ipid], "MRNQ");
      }
    }
  }

  // Fit (x,y): 3 2x2D gaus
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fit (x,y): 3 2x2D gaus..." << endl;
  TF2 *f2_gaus2_xy[NptBins][Npid];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    for (int ipid=0; ipid<Npid; ipid++) {
      if (ipid == 0 && ipt < 5)
        f2_gaus2_xy[ipt][ipid] = new TF2(Form("f2_gaus2_xy_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          -0.2, //f1_gaus2_x[ipt][ipid]->GetParameter(1)-4.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          0.2, // f1_gaus2_x[ipt][ipid]->GetParameter(1)+4.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          -0.2, //f1_gaus2_y[ipt][ipid]->GetParameter(1)-8.*f1_gaus2_y[ipt][ipid]->GetParameter(2),
          0.2); //f1_gaus2_y[ipt][ipid]->GetParameter(1)+8.*f1_gaus2_y[ipt][ipid]->GetParameter(2));
      if (ipid == 0 && ipt >= 5)
        f2_gaus2_xy[ipt][ipid] = new TF2(Form("f2_gaus2_xy_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          -0.2, // f1_gaus2_x[ipt][ipid]->GetParameter(1)-3.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          0.2, // f1_gaus2_x[ipt][ipid]->GetParameter(1)+3.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          -0.2, // f1_gaus2_y[ipt][ipid]->GetParameter(1)-4.*f1_gaus2_y[ipt][ipid]->GetParameter(2),
          0.2); // f1_gaus2_y[ipt][ipid]->GetParameter(1)+4.*f1_gaus2_y[ipt][ipid]->GetParameter(2));
      else if (ipid == 1 && ipt < 8)
        f2_gaus2_xy[ipt][ipid] = new TF2(Form("f2_gaus2_xy_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          f1_gaus2_x[ipt][ipid]->GetParameter(1)-4.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          f1_gaus2_x[ipt][ipid]->GetParameter(1)+4.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          -0.2, //f1_gaus2_y[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_y[ipt][ipid]->GetParameter(2),
          0.2); //f1_gaus2_y[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_y[ipt][ipid]->GetParameter(2));
      else if (ipid == 1 && ipt >= 8)
        f2_gaus2_xy[ipt][ipid] = new TF2(Form("f2_gaus2_xy_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          f1_gaus2_x[ipt][ipid]->GetParameter(1)-1.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          f1_gaus2_x[ipt][ipid]->GetParameter(1)+1.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          -0.2, // f1_gaus2_y[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_y[ipt][ipid]->GetParameter(2),
          0.2); // f1_gaus2_y[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_y[ipt][ipid]->GetParameter(2));
      else if (ipid == 2 && ipt < 8)
        f2_gaus2_xy[ipt][ipid] = new TF2(Form("f2_gaus2_xy_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          f1_gaus2_x[ipt][ipid]->GetParameter(1)-4.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          f1_gaus2_x[ipt][ipid]->GetParameter(1)+4.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          f1_gaus2_y[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_y[ipt][ipid]->GetParameter(2),
          f1_gaus2_y[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_y[ipt][ipid]->GetParameter(2));
      else
        f2_gaus2_xy[ipt][ipid] = new TF2(Form("f2_gaus2_xy_pt%i_pid%i", ipt, ipid), "xygaus(0)+xygaus(5)", 
          f1_gaus2_x[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          f1_gaus2_x[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_x[ipt][ipid]->GetParameter(2),
          f1_gaus2_y[ipt][ipid]->GetParameter(1)-2.*f1_gaus2_y[ipt][ipid]->GetParameter(2),
          f1_gaus2_y[ipt][ipid]->GetParameter(1)+2.*f1_gaus2_y[ipt][ipid]->GetParameter(2));

      f2_gaus2_xy[ipt][ipid]->SetParName(0, "SglConstant");
      f2_gaus2_xy[ipt][ipid]->SetParName(1, "SglMeanX");
      f2_gaus2_xy[ipt][ipid]->SetParName(2, "SglSigmaX");
      f2_gaus2_xy[ipt][ipid]->SetParName(3, "SglMeanY");
      f2_gaus2_xy[ipt][ipid]->SetParName(4, "SglSigmaY");
      f2_gaus2_xy[ipt][ipid]->SetParName(5, "BgdConstant");
      f2_gaus2_xy[ipt][ipid]->SetParName(6, "BgdMeanX");
      f2_gaus2_xy[ipt][ipid]->SetParName(7, "BgdSigmaX");
      f2_gaus2_xy[ipt][ipid]->SetParName(8, "BgdMeanY");
      f2_gaus2_xy[ipt][ipid]->SetParName(9, "BgdSigmaY");
      
      f2_gaus2_xy[ipt][ipid]->SetParameter(0, f1_gaus2_x[ipt][ipid]->GetParameter(0));
      if (ipid == 0) f2_gaus2_xy[ipt][ipid]->SetParameter(1, 0.);
      else f2_gaus2_xy[ipt][ipid]->SetParameter(1, f1_gaus2_x[ipt][ipid]->GetParameter(1));
      // f2_gaus2_xy[ipt][ipid]->FixParameter(1, f1_gaus2_x[ipt][ipid]->GetParameter(1));
      if (ipid == 0) f2_gaus2_xy[ipt][ipid]->SetParameter(2, 5.*f1_gaus2_x[ipt][ipid]->GetParameter(2));
      else f2_gaus2_xy[ipt][ipid]->SetParameter(2, f1_gaus2_x[ipt][ipid]->GetParameter(2));
      if (ipid != 2) f2_gaus2_xy[ipt][ipid]->SetParameter(3, 0.);
      else f2_gaus2_xy[ipt][ipid]->SetParameter(3, f1_gaus2_y[ipt][ipid]->GetParameter(1));
      f2_gaus2_xy[ipt][ipid]->SetParameter(4, f1_gaus2_y[ipt][ipid]->GetParameter(2));
      f2_gaus2_xy[ipt][ipid]->SetParameter(5, f1_gaus2_x[ipt][ipid]->GetParameter(3));
      if (ipid == 0) f2_gaus2_xy[ipt][ipid]->SetParameter(6, 0.);
      else f2_gaus2_xy[ipt][ipid]->SetParameter(6, f1_gaus2_x[ipt][ipid]->GetParameter(1));
      // f2_gaus2_xy[ipt][ipid]->FixParameter(6, f1_gaus2_x[ipt][ipid]->GetParameter(1));
      if (ipid == 0) f2_gaus2_xy[ipt][ipid]->SetParameter(7, 10.*f1_gaus2_x[ipt][ipid]->GetParameter(5));
      else f2_gaus2_xy[ipt][ipid]->SetParameter(7, f1_gaus2_x[ipt][ipid]->GetParameter(5));
      if (ipid != 2) f2_gaus2_xy[ipt][ipid]->SetParameter(8, 0.);
      else f2_gaus2_xy[ipt][ipid]->SetParameter(8, f1_gaus2_y[ipt][ipid]->GetParameter(1));
      f2_gaus2_xy[ipt][ipid]->SetParameter(9, f1_gaus2_y[ipt][ipid]->GetParameter(5));

      f2_gaus2_xy[ipt][ipid]->SetParLimits(0, 0., h_x[ipt]->GetEntries());
      if (ipid == 0) f2_gaus2_xy[ipt][ipid]->SetParLimits(1, -0.1, 0.1);
      else if (ipid == 1) f2_gaus2_xy[ipt][ipid]->SetParLimits(1, f1_gaus2_x[ipt][ipid]->GetParameter(1)*0.9, f1_gaus2_x[ipt][ipid]->GetParameter(1)*1.1);
      else f2_gaus2_xy[ipt][ipid]->SetParLimits(1, f1_gaus2_x[ipt][ipid]->GetParameter(1)*0.9, f1_gaus2_x[ipt][ipid]->GetParameter(1)*1.1);
      if (ipid == 0) f2_gaus2_xy[ipt][ipid]->SetParLimits(2, f1_gaus2_x[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_x[ipt][ipid]->GetParameter(2)*20.);
      else if (ipid == 1) f2_gaus2_xy[ipt][ipid]->SetParLimits(2, f1_gaus2_x[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_x[ipt][ipid]->GetParameter(2)*5.);
      else f2_gaus2_xy[ipt][ipid]->SetParLimits(2, f1_gaus2_x[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_x[ipt][ipid]->GetParameter(2)*1.5);
      if (ipid != 2) f2_gaus2_xy[ipt][ipid]->SetParLimits(3, -0.1, 0.1);
      else f2_gaus2_xy[ipt][ipid]->SetParLimits(3, f1_gaus2_y[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_y[ipt][ipid]->GetParameter(1)*1.5);
      f2_gaus2_xy[ipt][ipid]->SetParLimits(4, f1_gaus2_y[ipt][ipid]->GetParameter(2)*0.5, f1_gaus2_y[ipt][ipid]->GetParameter(2)*1.5);
      f2_gaus2_xy[ipt][ipid]->SetParLimits(5, 0., h_x[ipt]->GetEntries());
      if (ipid == 0) f2_gaus2_xy[ipt][ipid]->SetParLimits(6, f1_gaus2_x[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_x[ipt][ipid]->GetParameter(1)*1.5);
      else if (ipid == 1) f2_gaus2_xy[ipt][ipid]->SetParLimits(6, f1_gaus2_x[ipt][ipid]->GetParameter(1)*0.9, f1_gaus2_x[ipt][ipid]->GetParameter(1)*1.1);
      else f2_gaus2_xy[ipt][ipid]->SetParLimits(6, -0.1, 0.1);
      if (ipid == 0) f2_gaus2_xy[ipt][ipid]->SetParLimits(7, f1_gaus2_x[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_x[ipt][ipid]->GetParameter(5)*20.);
      else if (ipid == 1) f2_gaus2_xy[ipt][ipid]->SetParLimits(7, f1_gaus2_x[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_x[ipt][ipid]->GetParameter(5)*5.);
      else f2_gaus2_xy[ipt][ipid]->SetParLimits(7, f1_gaus2_x[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_x[ipt][ipid]->GetParameter(5)*1.5);
      if (ipid != 2) f2_gaus2_xy[ipt][ipid]->SetParLimits(8, -0.1, 0.1);
      else f2_gaus2_xy[ipt][ipid]->SetParLimits(8, f1_gaus2_y[ipt][ipid]->GetParameter(1)*0.5, f1_gaus2_y[ipt][ipid]->GetParameter(1)*1.5);
      f2_gaus2_xy[ipt][ipid]->SetParLimits(9, f1_gaus2_y[ipt][ipid]->GetParameter(5)*0.5, f1_gaus2_y[ipt][ipid]->GetParameter(5)*1.5);
    }
  }
  for (int ipt = 0; ipt < NptBins; ipt++) {
    cout << "\tPt bin " << ipt << endl;
    for (int ipid=0; ipid<Npid; ipid++) {
      for (int iter=0; iter<Niter2D; iter++) {
        h_pidXY[ipt]->Fit(f2_gaus2_xy[ipt][ipid], fitOptions.c_str());
      }
    }
    timer.Print();
    timer.Continue();
  }

  // Do polinomial fits for parameter smoothening for 3gaus - x,y
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fit (x,y): Do polinomial fits for parameter smoothening for 3gaus..." << endl;
  TGraphErrors *gr_signal1_mean_x[Npid];
  TGraphErrors *gr_signal1_mean_y[Npid];
  TGraphErrors *gr_signal1_sigm_x[Npid];
  TGraphErrors *gr_signal1_sigm_y[Npid];
  TGraphErrors *gr_orig_signal1_mean_x[Npid];
  TGraphErrors *gr_orig_signal1_mean_y[Npid];
  TGraphErrors *gr_orig_signal1_sigm_x[Npid];
  TGraphErrors *gr_orig_signal1_sigm_y[Npid];

  for (int ipid=0; ipid<Npid; ipid++) {
    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus2_xy[ipt][ipid]->GetParameter(1));
      vey_tmp.push_back(f2_gaus2_xy[ipt][ipid]->GetParError(1));
    }
    gr_orig_signal1_mean_x[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal1_mean_x[ipid]->SetName(Form("gr_orig_signal1_mean_x_pid%i", ipid));
    gr_orig_signal1_mean_x[ipid]->SetTitle(Form("gr_orig_signal1_mean_x_pid%i", ipid));
    gr_orig_signal1_mean_x[ipid]->SetMarkerStyle(kFullSquare);

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus2_xy[ipt][ipid]->GetParameter(3));
      vey_tmp.push_back(f2_gaus2_xy[ipt][ipid]->GetParError(3));
    }
    gr_orig_signal1_mean_y[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal1_mean_y[ipid]->SetName(Form("gr_orig_signal1_mean_y_pid%i", ipid));
    gr_orig_signal1_mean_y[ipid]->SetTitle(Form("gr_orig_signal1_mean_y_pid%i", ipid));
    gr_orig_signal1_mean_y[ipid]->SetMarkerStyle(kFullSquare);

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus2_xy[ipt][ipid]->GetParameter(2));
      vey_tmp.push_back(f2_gaus2_xy[ipt][ipid]->GetParError(2));
    }
    gr_orig_signal1_sigm_x[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal1_sigm_x[ipid]->SetName(Form("gr_orig_signal1_sigm_x_pid%i", ipid));
    gr_orig_signal1_sigm_x[ipid]->SetTitle(Form("gr_orig_signal1_sigm_x_pid%i", ipid));
    gr_orig_signal1_sigm_x[ipid]->SetMarkerStyle(kFullSquare);

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus2_xy[ipt][ipid]->GetParameter(4));
      vey_tmp.push_back(f2_gaus2_xy[ipt][ipid]->GetParError(4));
    }
    gr_orig_signal1_sigm_y[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal1_sigm_y[ipid]->SetName(Form("gr_orig_signal1_sigm_y_pid%i", ipid));
    gr_orig_signal1_sigm_y[ipid]->SetTitle(Form("gr_orig_signal1_sigm_y_pid%i", ipid));
    gr_orig_signal1_sigm_y[ipid]->SetMarkerStyle(kFullSquare);
  }

  // Smoothen the pT-dependencies
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_signal1_mean_x[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal1_mean_x[ipid], -1);
    gr_signal1_mean_y[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal1_mean_y[ipid], -1);
    gr_signal1_sigm_x[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal1_sigm_x[ipid], -1);
    gr_signal1_sigm_y[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal1_sigm_y[ipid], -1);
  }

  // Fit (x,y): 2x2D 3gaus
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fit (x,y): 2x2D 3gaus..." << endl;
  TF2 *f2_gaus3_xy[NptBins];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    f2_gaus3_xy[ipt] = new TF2(Form("f2_gaus3_xy_pt%i", ipt), "xygaus(0)+xygaus(5) + xygaus(10)+xygaus(15) + xygaus(20)+xygaus(25)", -0.5, 1.5, -0.5, 0.5);

    f2_gaus3_xy[ipt]->SetParName(0, "Peak1SglConstant");
    f2_gaus3_xy[ipt]->SetParName(1, "Peak1SglMeanX");
    f2_gaus3_xy[ipt]->SetParName(2, "Peak1SglSigmaX");
    f2_gaus3_xy[ipt]->SetParName(3, "Peak1SglMeanY");
    f2_gaus3_xy[ipt]->SetParName(4, "Peak1SglSigmaY");
    f2_gaus3_xy[ipt]->SetParName(5, "Peak1BgdConstant");
    f2_gaus3_xy[ipt]->SetParName(6, "Peak1BgdMeanX");
    f2_gaus3_xy[ipt]->SetParName(7, "Peak1BgdSigmaX");
    f2_gaus3_xy[ipt]->SetParName(8, "Peak1BgdMeanY");
    f2_gaus3_xy[ipt]->SetParName(9, "Peak1BgdSigmaY");
    f2_gaus3_xy[ipt]->SetParName(10, "Peak2SglConstant");
    f2_gaus3_xy[ipt]->SetParName(11, "Peak2SglMeanX");
    f2_gaus3_xy[ipt]->SetParName(12, "Peak2SglSigmaX");
    f2_gaus3_xy[ipt]->SetParName(13, "Peak2SglMeanY");
    f2_gaus3_xy[ipt]->SetParName(14, "Peak2SglSigmaY");
    f2_gaus3_xy[ipt]->SetParName(15, "Peak2BgdConstant");
    f2_gaus3_xy[ipt]->SetParName(16, "Peak2BgdMeanX");
    f2_gaus3_xy[ipt]->SetParName(17, "Peak2BgdSigmaX");
    f2_gaus3_xy[ipt]->SetParName(18, "Peak2BgdMeanY");
    f2_gaus3_xy[ipt]->SetParName(19, "Peak2BgdSigmaY");
    f2_gaus3_xy[ipt]->SetParName(20, "Peak3SglConstant");
    f2_gaus3_xy[ipt]->SetParName(21, "Peak3SglMeanX");
    f2_gaus3_xy[ipt]->SetParName(22, "Peak3SglSigmaX");
    f2_gaus3_xy[ipt]->SetParName(23, "Peak3SglMeanY");
    f2_gaus3_xy[ipt]->SetParName(24, "Peak3SglSigmaY");
    f2_gaus3_xy[ipt]->SetParName(25, "Peak3BgdConstant");
    f2_gaus3_xy[ipt]->SetParName(26, "Peak3BgdMeanX");
    f2_gaus3_xy[ipt]->SetParName(27, "Peak3BgdSigmaX");
    f2_gaus3_xy[ipt]->SetParName(28, "Peak3BgdMeanY");
    f2_gaus3_xy[ipt]->SetParName(29, "Peak3BgdSigmaY");

    f2_gaus3_xy[ipt]->SetParameter(0, f2_gaus2_xy[ipt][0]->GetParameter(0));
    f2_gaus3_xy[ipt]->FixParameter(1, 0.); 
    // f2_gaus3_xy[ipt]->SetParameter(1, f_signal1_mean_m2[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(2, gr_signal1_sigm_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->FixParameter(3, 0.); 
    // f2_gaus3_xy[ipt]->SetParameter(3, gr_signal1_mean_y[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(4, gr_signal1_sigm_y[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(5, f2_gaus2_xy[ipt][0]->GetParameter(5));
    f2_gaus3_xy[ipt]->FixParameter(6, 0.); 
    // f2_gaus3_xy[ipt]->SetParameter(6, gr_signal1_mean_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(7, f2_gaus2_xy[ipt][0]->GetParameter(7));
    f2_gaus3_xy[ipt]->FixParameter(8, 0.); 
    // f2_gaus3_xy[ipt]->SetParameter(8, f2_gaus2_xy[ipt][0]->GetParameter(8));
    f2_gaus3_xy[ipt]->SetParameter(9, f2_gaus2_xy[ipt][0]->GetParameter(9));
    f2_gaus3_xy[ipt]->SetParameter(10, f2_gaus2_xy[ipt][1]->GetParameter(0));
    // f2_gaus3_xy[ipt]->FixParameter(11, pid_Msqr.at(1)); 
    f2_gaus3_xy[ipt]->SetParameter(11, gr_signal1_mean_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(12, gr_signal1_sigm_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(13, 0.);
    // f2_gaus3_xy[ipt]->SetParameter(13, gr_signal1_mean_y[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(14, gr_signal1_sigm_y[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(15, f2_gaus2_xy[ipt][1]->GetParameter(5));
    // f2_gaus3_xy[ipt]->FixParameter(16, pid_Msqr.at(1)); 
    f2_gaus3_xy[ipt]->SetParameter(16, gr_signal1_mean_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(17, f2_gaus2_xy[ipt][1]->GetParameter(7));
    f2_gaus3_xy[ipt]->SetParameter(18, 0.);
    // f2_gaus3_xy[ipt]->SetParameter(18, gr_signal1_mean_y[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(19, f2_gaus2_xy[ipt][1]->GetParameter(9));
    f2_gaus3_xy[ipt]->SetParameter(20, f2_gaus2_xy[ipt][2]->GetParameter(0));
    // f2_gaus3_xy[ipt]->FixParameter(21, pid_Msqr.at(2)); 
    f2_gaus3_xy[ipt]->SetParameter(21, gr_signal1_mean_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(22, gr_signal1_sigm_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(23, f2_gaus2_xy[ipt][2]->GetParameter(3));
    f2_gaus3_xy[ipt]->SetParameter(24, f2_gaus2_xy[ipt][2]->GetParameter(4));
    f2_gaus3_xy[ipt]->SetParameter(25, f2_gaus2_xy[ipt][2]->GetParameter(5));
    // f2_gaus3_xy[ipt]->FixParameter(26, pid_Msqr.at(2)); 
    f2_gaus3_xy[ipt]->SetParameter(26, gr_signal1_mean_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f2_gaus3_xy[ipt]->SetParameter(27, f2_gaus2_xy[ipt][2]->GetParameter(7));
    f2_gaus3_xy[ipt]->SetParameter(28, f2_gaus2_xy[ipt][2]->GetParameter(8));
    f2_gaus3_xy[ipt]->SetParameter(29, f2_gaus2_xy[ipt][2]->GetParameter(9));

    f2_gaus3_xy[ipt]->SetParLimits(0, f2_gaus3_xy[ipt]->GetParameter(0)*0.5, f2_gaus3_xy[ipt]->GetParameter(0)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(1, f2_gaus3_xy[ipt]->GetParameter(1)*0.5, f2_gaus3_xy[ipt]->GetParameter(1)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(2, f2_gaus3_xy[ipt]->GetParameter(2)*0.5, f2_gaus3_xy[ipt]->GetParameter(2)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(3, f2_gaus3_xy[ipt]->GetParameter(3)*0.5, f2_gaus3_xy[ipt]->GetParameter(3)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(4, f2_gaus3_xy[ipt]->GetParameter(4)*0.5, f2_gaus3_xy[ipt]->GetParameter(4)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(5, f2_gaus3_xy[ipt]->GetParameter(5)*0.5, f2_gaus3_xy[ipt]->GetParameter(5)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(6, f2_gaus3_xy[ipt]->GetParameter(6)*0.5, f2_gaus3_xy[ipt]->GetParameter(6)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(7, f2_gaus3_xy[ipt]->GetParameter(7)*0.5, f2_gaus3_xy[ipt]->GetParameter(7)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(8, f2_gaus3_xy[ipt]->GetParameter(8)*0.5, f2_gaus3_xy[ipt]->GetParameter(8)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(9, f2_gaus3_xy[ipt]->GetParameter(9)*0.5, f2_gaus3_xy[ipt]->GetParameter(9)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(10, f2_gaus3_xy[ipt]->GetParameter(10)*0.5, f2_gaus3_xy[ipt]->GetParameter(10)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(11, f2_gaus3_xy[ipt]->GetParameter(11)*0.5, f2_gaus3_xy[ipt]->GetParameter(11)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(12, f2_gaus3_xy[ipt]->GetParameter(12)*0.5, f2_gaus3_xy[ipt]->GetParameter(12)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(13, f2_gaus3_xy[ipt]->GetParameter(13)*0.5, f2_gaus3_xy[ipt]->GetParameter(13)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(14, f2_gaus3_xy[ipt]->GetParameter(14)*0.5, f2_gaus3_xy[ipt]->GetParameter(14)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(15, f2_gaus3_xy[ipt]->GetParameter(15)*0.5, f2_gaus3_xy[ipt]->GetParameter(15)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(16, f2_gaus3_xy[ipt]->GetParameter(16)*0.5, f2_gaus3_xy[ipt]->GetParameter(16)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(17, f2_gaus3_xy[ipt]->GetParameter(17)*0.5, f2_gaus3_xy[ipt]->GetParameter(17)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(18, f2_gaus3_xy[ipt]->GetParameter(18)*0.5, f2_gaus3_xy[ipt]->GetParameter(18)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(19, f2_gaus3_xy[ipt]->GetParameter(19)*0.5, f2_gaus3_xy[ipt]->GetParameter(19)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(20, f2_gaus3_xy[ipt]->GetParameter(20)*0.5, f2_gaus3_xy[ipt]->GetParameter(20)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(21, f2_gaus3_xy[ipt]->GetParameter(21)*0.5, f2_gaus3_xy[ipt]->GetParameter(21)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(22, f2_gaus3_xy[ipt]->GetParameter(22)*0.5, f2_gaus3_xy[ipt]->GetParameter(22)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(23, f2_gaus3_xy[ipt]->GetParameter(23)*0.5, f2_gaus3_xy[ipt]->GetParameter(23)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(24, f2_gaus3_xy[ipt]->GetParameter(24)*0.5, f2_gaus3_xy[ipt]->GetParameter(24)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(25, f2_gaus3_xy[ipt]->GetParameter(25)*0.5, f2_gaus3_xy[ipt]->GetParameter(25)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(26, f2_gaus3_xy[ipt]->GetParameter(26)*0.5, f2_gaus3_xy[ipt]->GetParameter(26)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(27, f2_gaus3_xy[ipt]->GetParameter(27)*0.5, f2_gaus3_xy[ipt]->GetParameter(27)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(28, f2_gaus3_xy[ipt]->GetParameter(28)*0.5, f2_gaus3_xy[ipt]->GetParameter(28)*1.5);
    f2_gaus3_xy[ipt]->SetParLimits(29, f2_gaus3_xy[ipt]->GetParameter(29)*0.5, f2_gaus3_xy[ipt]->GetParameter(29)*1.5);
  }
  for (int ipt = 0; ipt < NptBins; ipt++) {
    cout << "\tPt bin " << ipt <<endl;
    for (int iter=0; iter<Niter2D; iter++) {
        if (!fast_check) h_pidXY[ipt]->Fit(f2_gaus3_xy[ipt],Form("%s B",fitOptions.c_str()));
    }
    timer.Print();
    timer.Continue();
  }

  // Do polinomial fits for parameter smoothening for pion, kaon selection - x,y
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fit (x,y): Do polinomial fits for pion, kaon selection..." << endl;
  TGraphErrors *gr_signal2_mean_x[Npid];
  TGraphErrors *gr_signal2_mean_y[Npid];
  TGraphErrors *gr_signal2_sigm_x[Npid];
  TGraphErrors *gr_signal2_sigm_y[Npid];
  TGraphErrors *gr_orig_signal2_mean_x[Npid];
  TGraphErrors *gr_orig_signal2_mean_y[Npid];
  TGraphErrors *gr_orig_signal2_sigm_x[Npid];
  TGraphErrors *gr_orig_signal2_sigm_y[Npid];

  for (int ipid=0; ipid<Npid; ipid++) {
    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus3_xy[ipt]->GetParameter(1+10*ipid));
      vey_tmp.push_back(f2_gaus3_xy[ipt]->GetParError(1+10*ipid));
    }
    gr_orig_signal2_mean_x[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal2_mean_x[ipid]->SetName(Form("gr_orig_signal2_mean_x_pid%i", ipid));
    gr_orig_signal2_mean_x[ipid]->SetTitle(Form("gr_orig_signal2_mean_x_pid%i", ipid));
    gr_orig_signal2_mean_x[ipid]->SetMarkerStyle(kFullSquare);

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus3_xy[ipt]->GetParameter(3+10*ipid));
      vey_tmp.push_back(f2_gaus3_xy[ipt]->GetParError(3+10*ipid));
    }
    gr_orig_signal2_mean_y[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal2_mean_y[ipid]->SetName(Form("gr_orig_signal2_mean_y_pid%i", ipid));
    gr_orig_signal2_mean_y[ipid]->SetTitle(Form("gr_orig_signal2_mean_y_pid%i", ipid));
    gr_orig_signal2_mean_y[ipid]->SetMarkerStyle(kFullSquare);

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus3_xy[ipt]->GetParameter(2+10*ipid));
      vey_tmp.push_back(f2_gaus3_xy[ipt]->GetParError(2+10*ipid));
    }
    gr_orig_signal2_sigm_x[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal2_sigm_x[ipid]->SetName(Form("gr_orig_signal2_sigm_x_pid%i", ipid));
    gr_orig_signal2_sigm_x[ipid]->SetTitle(Form("gr_orig_signal2_sigm_x_pid%i", ipid));
    gr_orig_signal2_sigm_x[ipid]->SetMarkerStyle(kFullSquare);

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      vy_tmp.push_back(f2_gaus3_xy[ipt]->GetParameter(4+10*ipid));
      vey_tmp.push_back(f2_gaus3_xy[ipt]->GetParError(4+10*ipid));
    }
    gr_orig_signal2_sigm_y[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal2_sigm_y[ipid]->SetName(Form("gr_orig_signal2_sigm_y_pid%i", ipid));
    gr_orig_signal2_sigm_y[ipid]->SetTitle(Form("gr_orig_signal2_sigm_y_pid%i", ipid));
    gr_orig_signal2_sigm_y[ipid]->SetMarkerStyle(kFullSquare);
  }

  // Smoothen the pT-dependencies
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_signal2_mean_x[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal2_mean_x[ipid], -1);
    gr_signal2_mean_y[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal2_mean_y[ipid], -1);
    gr_signal2_sigm_x[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal2_sigm_x[ipid], -1);
    gr_signal2_sigm_y[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal2_sigm_y[ipid], -1);
  }

  // Cleanup (x,y) of pions and kaons via resampling
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Cleanup (x,y) of pions and kaons via resampling..." << endl;
  TH2D *h_pika1XY[NptBins];
  double x_coord, y_coord;
  for (int ipt = 0; ipt < NptBins; ipt++) {
    cout << "\tPt bin " << ipt <<endl;

    h_pika1XY[ipt] = new TH2D(Form("h_pika1XY_pt%i", ipt), Form("h_pika1XY for %1.1f < p_{T} < %1.1f GeV/c (#pi, K only);X;Y", ptBins[ipt], ptBins[ipt+1]), 6000, -1.5, 1.5, 6000, -1.5, 1.5);

    long nentries = (long)h_pikaXY[ipt]->GetEntries();
    for (long ientry=0; ientry<nentries; ientry++) {
      h_pikaXY[ipt]->GetRandom2(x_coord,y_coord);
      
      double sigmx_pi = gr_signal2_sigm_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double sigmx_ka = gr_signal2_sigm_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double sigmx_pr = gr_signal2_sigm_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double meanx_pi = gr_signal2_mean_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double meanx_ka = gr_signal2_mean_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double meanx_pr = gr_signal2_mean_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double sigmy_pi = gr_signal2_sigm_y[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double sigmy_ka = gr_signal2_sigm_y[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double sigmy_pr = gr_signal2_sigm_y[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double meany_pi = gr_signal2_mean_y[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double meany_ka = gr_signal2_mean_y[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);
      double meany_pr = gr_signal2_mean_y[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.);

      // Cut everything outside 5*sigma for pions and kaons + cut 3*sigma protons
      if ( (pow(x_coord - meanx_pi, 2)/pow(5.*sigmx_pi, 2) + pow(y_coord - meany_pi, 2)/pow(5.*sigmy_pi, 2) < 1. ||
            pow(x_coord - meanx_ka, 2)/pow(5.*sigmx_ka, 2) + pow(y_coord - meany_ka, 2)/pow(5.*sigmy_ka, 2) < 1. ) &&
            pow(x_coord - meanx_pr, 2)/pow(3.*sigmx_pr, 2) + pow(y_coord - meany_pr, 2)/pow(3.*sigmy_pr, 2) > 1.) {
        h_pika1XY[ipt]->Fill(x_coord, y_coord);
      }
    }

    timer.Print();
    timer.Continue();
  }

  // Fit x: 2x1D gaus for pions and kaons
  // --- ToDo
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fit (x,y): 3 2x1D gaus..." << endl;
  TH1D *h_pika_x[NptBins];
  TF1 *f1_gaus_pika_x[NptBins];
  TF1 *f1_gaus2_pika_x[NptBins];
  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_pika_x[ipt] = (TH1D*) h_pika1XY[ipt]->ProjectionX(Form("h_pika_x_pt%i", ipt), 
    h_pika1XY[ipt]->GetYaxis()->FindBin(gr_signal2_mean_y[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.) - 2.5*gr_signal2_sigm_y[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.)), 
    h_pika1XY[ipt]->GetYaxis()->FindBin(gr_signal2_mean_y[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.) + 2.5*gr_signal2_sigm_y[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.)));

    f1_gaus_pika_x[ipt] = new TF1(Form("f1_gaus_pika_x_pt%i", ipt), "gaus(0)+gaus(3)", 
    gr_signal2_mean_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.) - 2.5*gr_signal2_sigm_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.),
    gr_signal2_mean_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.) + 2.5*gr_signal2_sigm_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f1_gaus2_pika_x[ipt] = new TF1(Form("f1_gaus2_pika_x_pt%i", ipt), "gaus(0)+gaus(3)+gaus(6)+gaus(9)", 
    gr_signal2_mean_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.) - 2.5*gr_signal2_sigm_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.),
    gr_signal2_mean_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.) + 2.5*gr_signal2_sigm_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));

    f1_gaus_pika_x[ipt]->SetParName(0, "Peak1Constant");
    f1_gaus_pika_x[ipt]->SetParName(1, "Peak1Mean");
    f1_gaus_pika_x[ipt]->SetParName(2, "Peak1Sigma");
    f1_gaus_pika_x[ipt]->SetParName(3, "Peak2Constant");
    f1_gaus_pika_x[ipt]->SetParName(4, "Peak2Mean");
    f1_gaus_pika_x[ipt]->SetParName(5, "Peak2Sigma");

    f1_gaus2_pika_x[ipt]->SetParName(0, "Peak1Constant");
    f1_gaus2_pika_x[ipt]->SetParName(1, "Peak1Mean");
    f1_gaus2_pika_x[ipt]->SetParName(2, "Peak1Sigma");
    f1_gaus2_pika_x[ipt]->SetParName(3, "Peak2Constant");
    f1_gaus2_pika_x[ipt]->SetParName(4, "Peak2Mean");
    f1_gaus2_pika_x[ipt]->SetParName(5, "Peak2Sigma");
    f1_gaus2_pika_x[ipt]->SetParName(6, "Bgd1Constant");
    f1_gaus2_pika_x[ipt]->SetParName(7, "Bgd1Mean");
    f1_gaus2_pika_x[ipt]->SetParName(8, "Bgd1Sigma");
    f1_gaus2_pika_x[ipt]->SetParName(9, "Bgd2Constant");
    f1_gaus2_pika_x[ipt]->SetParName(10, "Bgd2Mean");
    f1_gaus2_pika_x[ipt]->SetParName(11, "Bgd2Sigma");

    f1_gaus_pika_x[ipt]->SetParameter(0, f1_gaus2_x[ipt][0]->GetParameter(0));
    f1_gaus_pika_x[ipt]->FixParameter(1, 0.);
    // f1_gaus_pika_x[ipt]->SetParameter(1, 0.);
    f1_gaus_pika_x[ipt]->SetParameter(2, gr_signal1_sigm_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f1_gaus_pika_x[ipt]->SetParameter(3, f1_gaus2_x[ipt][1]->GetParameter(0));
    f1_gaus_pika_x[ipt]->SetParameter(4, 0.25);
    // f1_gaus_pika_x[ipt]->SetParameter(4, gr_signal1_mean_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f1_gaus_pika_x[ipt]->SetParameter(5, gr_signal1_sigm_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));

    f1_gaus2_pika_x[ipt]->SetParameter(0, f1_gaus2_x[ipt][0]->GetParameter(0));
    f1_gaus2_pika_x[ipt]->FixParameter(1, 0.);
    // f1_gaus2_pika_x[ipt]->SetParameter(1, 0.);
    f1_gaus2_pika_x[ipt]->SetParameter(2, gr_signal1_sigm_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f1_gaus2_pika_x[ipt]->SetParameter(3, f1_gaus2_x[ipt][1]->GetParameter(0));
    f1_gaus2_pika_x[ipt]->SetParameter(4, 0.25);
    // f1_gaus2_pika_x[ipt]->SetParameter(4, gr_signal1_mean_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f1_gaus2_pika_x[ipt]->SetParameter(5, gr_signal1_sigm_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f1_gaus2_pika_x[ipt]->SetParameter(6, f1_gaus2_x[ipt][0]->GetParameter(0)/10.);
    f1_gaus2_pika_x[ipt]->SetParameter(7, 0.);
    f1_gaus2_pika_x[ipt]->SetParameter(8, 10.*gr_signal1_sigm_x[0]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f1_gaus2_pika_x[ipt]->SetParameter(9, f1_gaus2_x[ipt][1]->GetParameter(0)/10.);
    f1_gaus2_pika_x[ipt]->SetParameter(10, 0.25);
    // f1_gaus2_pika_x[ipt]->SetParameter(10, gr_signal1_mean_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    f1_gaus2_pika_x[ipt]->SetParameter(11, 5.*gr_signal1_sigm_x[1]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
    

    // f1_gaus_pika_x[ipt]->SetParLimits(0, 0., f1_gaus_pika_x[ipt]->GetParameter(0)*100.);
    // f1_gaus_pika_x[ipt]->SetParLimits(1, -0.05, 0.05);
    // f1_gaus_pika_x[ipt]->SetParLimits(2, f1_gaus_pika_x[ipt]->GetParameter(2)*0.5, f1_gaus_pika_x[ipt]->GetParameter(2)*1.5);
    // f1_gaus_pika_x[ipt]->SetParLimits(3, 0., f1_gaus_pika_x[ipt]->GetParameter(3)*100.);
    // f1_gaus_pika_x[ipt]->SetParLimits(4, f1_gaus_pika_x[ipt]->GetParameter(4)*0.5, f1_gaus_pika_x[ipt]->GetParameter(4)*1.5);
    // f1_gaus_pika_x[ipt]->SetParLimits(5, f1_gaus_pika_x[ipt]->GetParameter(5)*0.5, f1_gaus_pika_x[ipt]->GetParameter(5)*1.5);

    for (int iter=0; iter<Niter1D; iter++) {
      h_pika_x[ipt]->Fit(f1_gaus_pika_x[ipt], "MRNQ");
      h_pika_x[ipt]->Fit(f1_gaus2_pika_x[ipt], "MRNQ");
    }
  }

  // Do final polinomial fits for parameter smoothening for pion, kaon selection - x,y
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Fit (x,y): Do final polinomial fits for pion, kaon selection..." << endl;
  TGraphErrors *gr_signal3_mean_x[Npid];
  TGraphErrors *gr_signal3_sigm_x[Npid];
  TGraphErrors *gr_orig_signal3_mean_x[Npid];
  TGraphErrors *gr_orig_signal3_sigm_x[Npid];

  for (int ipid=0; ipid<Npid; ipid++) {
    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      if (ipid < 2) {
        vy_tmp.push_back(f1_gaus2_pika_x[ipt]->GetParameter(1+3*ipid));
        vey_tmp.push_back(f1_gaus2_pika_x[ipt]->GetParError(1+3*ipid));
      } else {
        vy_tmp.push_back(gr_signal1_mean_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
        vey_tmp.push_back(gr_signal1_mean_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
      }
    }
    gr_orig_signal3_mean_x[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal3_mean_x[ipid]->SetName(Form("gr_orig_signal3_mean_x_pid%i", ipid));
    gr_orig_signal3_mean_x[ipid]->SetTitle(Form("gr_orig_signal3_mean_x_pid%i", ipid));
    gr_orig_signal3_mean_x[ipid]->SetMarkerStyle(kFullSquare);

    vx_tmp.clear(); vy_tmp.clear(); vex_tmp.clear(); vey_tmp.clear();
    for (int ipt=0; ipt<NptBins; ipt++) {
      vx_tmp.push_back((ptBins[ipt+1]+ptBins[ipt])/2.);
      vex_tmp.push_back(0.);
      if (ipid < 2) {
        vy_tmp.push_back(f1_gaus2_pika_x[ipt]->GetParameter(2+3*ipid));
        vey_tmp.push_back(f1_gaus2_pika_x[ipt]->GetParError(2+3*ipid));
      } else {
        vy_tmp.push_back(gr_signal1_sigm_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
        vey_tmp.push_back(gr_signal1_sigm_x[2]->Eval((ptBins[ipt+1]+ptBins[ipt])/2.));
      }
    }
    gr_orig_signal3_sigm_x[ipid] = new TGraphErrors(vx_tmp.size(), &vx_tmp[0], &vy_tmp[0], &vex_tmp[0], &vey_tmp[0]);
    gr_orig_signal3_sigm_x[ipid]->SetName(Form("gr_orig_signal3_sigm_x_pid%i", ipid));
    gr_orig_signal3_sigm_x[ipid]->SetTitle(Form("gr_orig_signal3_sigm_x_pid%i", ipid));
    gr_orig_signal3_sigm_x[ipid]->SetMarkerStyle(kFullSquare);
  }

  // Smoothen the pT-dependencies
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_signal3_mean_x[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal3_mean_x[ipid], -1);
    gr_signal3_sigm_x[ipid] = (TGraphErrors*) DoLaplaceSmooth(gr_orig_signal3_sigm_x[ipid], -1);
  }

  // Writing output
  cout << "////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "Writing output..." << endl;
  fo->mkdir("nsig_m2_1D");
  fo->cd("nsig_m2_1D");
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
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_orig_signal1_mean_m2[ipid]->Write();
    gr_orig_signal1_mean_ns[ipid]->Write();
    gr_orig_signal1_sigm_m2[ipid]->Write();
    gr_orig_signal1_sigm_ns[ipid]->Write();

    gr_signal1_mean_m2[ipid]->Write();
    gr_signal1_mean_ns[ipid]->Write();
    gr_signal1_sigm_m2[ipid]->Write();
    gr_signal1_sigm_ns[ipid]->Write();

    f_signal1_mean_m2[ipid]->Write();
    f_signal1_mean_ns[ipid]->Write();
    f_signal1_sigm_m2[ipid]->Write();
    f_signal1_sigm_ns[ipid]->Write();
  }
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_orig_signal2_mean_m2[ipid]->Write();
    gr_orig_signal2_mean_ns[ipid]->Write();
    gr_orig_signal2_sigm_m2[ipid]->Write();
    gr_orig_signal2_sigm_ns[ipid]->Write();

    gr_signal2_mean_m2[ipid]->Write();
    gr_signal2_mean_ns[ipid]->Write();
    gr_signal2_sigm_m2[ipid]->Write();
    gr_signal2_sigm_ns[ipid]->Write();

    f_signal2_mean_m2[ipid]->Write();
    f_signal2_mean_ns[ipid]->Write();
    f_signal2_sigm_m2[ipid]->Write();
    f_signal2_sigm_ns[ipid]->Write();
  }
  fo->mkdir("nsig_m2_2D");
  fo->cd("nsig_m2_2D");
  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_NsigPiMsqr[ipt]->Write();
    f2_gaus3_Msqr[ipt]->Write();
    for (int ipid=0; ipid<Npid; ipid++) {
      // f2_gaus1_Msqr[ipt][ipid]->Write();
      f2_gaus2_Msqr[ipt][ipid]->Write();
    }
  }

  fo->mkdir("x_y_1D");
  fo->cd("x_y_1D");
  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_x[ipt]->Write();
    h_pika_x[ipt]->Write();
    f1_gaus_pika_x[ipt]->Write();
    f1_gaus2_pika_x[ipt]->Write();
    for (int ipid=0; ipid<Npid; ipid++) {
      h_y[ipt][ipid]->Write();
      f1_gaus1_x[ipt][ipid]->Write();
      f1_gaus2_x[ipt][ipid]->Write();
      f1_gaus1_y[ipt][ipid]->Write();
      f1_gaus2_y[ipt][ipid]->Write();
    }
  }
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_orig_signal1_mean_x[ipid]->Write();
    gr_orig_signal1_mean_y[ipid]->Write();
    gr_orig_signal1_sigm_x[ipid]->Write();
    gr_orig_signal1_sigm_y[ipid]->Write();

    gr_signal1_mean_x[ipid]->Write();
    gr_signal1_mean_y[ipid]->Write();
    gr_signal1_sigm_x[ipid]->Write();
    gr_signal1_sigm_y[ipid]->Write();
  }
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_orig_signal2_mean_x[ipid]->Write();
    gr_orig_signal2_mean_y[ipid]->Write();
    gr_orig_signal2_sigm_x[ipid]->Write();
    gr_orig_signal2_sigm_y[ipid]->Write();

    gr_signal2_mean_x[ipid]->Write();
    gr_signal2_mean_y[ipid]->Write();
    gr_signal2_sigm_x[ipid]->Write();
    gr_signal2_sigm_y[ipid]->Write();
  }
  for (int ipid=0; ipid<Npid; ipid++) {
    gr_orig_signal3_mean_x[ipid]->Write();
    gr_orig_signal3_sigm_x[ipid]->Write();

    gr_signal3_mean_x[ipid]->Write();
    gr_signal3_sigm_x[ipid]->Write();
  }

  fo->mkdir("x_y_2D");
  fo->cd("x_y_2D");
  for (int ipt = 0; ipt < NptBins; ipt++) {
    h_pidXY[ipt]->Write();
    h_pikaXY[ipt]->Write();
    h_pika1XY[ipt]->Write();
    f2_gaus3_xy[ipt]->Write();
    for (int ipid=0; ipid<Npid; ipid++) {
      f2_gaus2_xy[ipt][ipid]->Write();
    }
  }

  timer.Stop();
  timer.Print();
}








pair<double, double> GetXY(double _val_m2, double _val_sig, 
  double _fit_mean_m2_pi, double _fit_mean_m2_ka,
  double _fit_mean_sig_pi, double _fit_mean_sig_ka,
  double _fit_sig_m2_pi,
  double _fit_sig_sig_pi)
{
  double f = _fit_sig_sig_pi / _fit_sig_m2_pi;
  // double alpha = -1.*tanh( (_fit_mean_m2_ka - _fit_mean_m2_pi)/((_fit_mean_sig_ka - _fit_mean_sig_pi)/f) );
  // double alpha = -1.*atan2( (_fit_mean_m2_ka - _fit_mean_m2_pi), ((_fit_mean_sig_ka - _fit_mean_sig_pi)/f) );
  double alpha = -1.*atan2( (_fit_mean_m2_ka - _fit_mean_m2_pi), ((_fit_mean_sig_ka - _fit_mean_sig_pi)/f) );

  double x0 = (_val_sig - _fit_mean_sig_pi)/f;
  double y0 = _val_m2 - _fit_mean_m2_pi;
  // return {x0, y0};

  double x = cos(alpha)*x0 - sin(alpha)*y0;
  double y = sin(alpha)*x0 + cos(alpha)*y0;
  return {x,y};
}

TGraphErrors *DoLaplaceSmooth(TGraphErrors *const& gr, int niter=-1)
{
  TGraphErrors *result;

  vector<double> vx, vy, vex, vey, weight;

  for (int i=0; i<gr->GetN(); i++) {
    vx.push_back(gr->GetPointX(i));
    vy.push_back(gr->GetPointY(i));
    vex.push_back(gr->GetErrorX(i));
    vey.push_back(gr->GetErrorY(i));
    weight.push_back(1./vey.at(i));
  }

  int nsteps = (niter <= 0) ? gr->GetN()/2 : niter;
  for (int istep=0; istep<nsteps; istep++) {
    DoLaplaceSmooth(vy, istep);
    // DoLaplaceSmooth(vy, istep, weight);
  }

  result = new TGraphErrors(vx.size(), &vx[0], &vy[0], &vex[0], &vey[0]);
  result->SetMarkerStyle(gr->GetMarkerStyle());
  result->SetMarkerSize(gr->GetMarkerSize());
  result->SetMarkerColor(gr->GetMarkerColor());
  result->SetLineColor(gr->GetLineColor());
  result->SetLineStyle(gr->GetLineStyle());
  result->SetLineWidth(gr->GetLineWidth());
  result->SetName(Form("%s_smoothed",gr->GetName()));
  result->SetTitle(Form("smoothed %s",gr->GetTitle()));

  return (TGraphErrors*)result->Clone();
}

void DoLaplaceSmooth(vector<double> &vect, int istep, vector<double> weight)
{
  int ibegin = istep;
  int iend = vect.size()-istep;
  for (int i = ibegin; i < iend; i++) {
    if (i == ibegin) continue;
    if (i == iend-1) continue;
    double weight1 = 1.;
    double weight2 = 1.;
    if (weight.size() == vect.size())
    {
      weight1 = weight.at(i-1);
      weight2 = weight.at(i+1);
    }
    vect.at(i) = (weight1*vect.at(i-1)+weight2*vect.at(i+1))/(weight1+weight2);
  }
}