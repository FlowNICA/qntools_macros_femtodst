#include "/mnt/pool/nica/7/parfenovpeter/Soft/QnTools/install/include/QnTools/QnDataFrame.hpp"
//#include "QnDataFrame.hpp"
#include "utils.h"

vector <vector<string>> QvQv=
{
  {"Q_TPC_L_ch_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"Q_TPC_L_ch_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"Q_TPC_L_ch_TWIST", "Q_TPC_R_ch_TWIST"},
  {"Q_TPC_L_ch_RESCALED", "Q_TPC_R_ch_RESCALED"}
};

vector<vector<string>> Qv=
{
  {"Q_TPC_F_ch_RESCALED"}
};

vector <vector<string>> uvQv=
{
  {"u_TPC_L_ch_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_ch_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_ch_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_ch_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_ch_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_ch_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_ch_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_ch_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_pion_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_pion_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_pion_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_pion_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_pion_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_pion_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_pion_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_pion_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_kaon_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_kaon_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_kaon_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_kaon_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_kaon_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_kaon_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_kaon_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_kaon_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_proton_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_proton_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_proton_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_proton_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_proton_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_proton_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_proton_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_proton_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_tofpion_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_tofpion_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_tofpion_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_tofpion_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_tofpion_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_tofpion_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_tofpion_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_tofpion_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_tofkaon_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_tofkaon_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_tofkaon_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_tofkaon_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_tofkaon_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_tofkaon_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_tofkaon_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_tofkaon_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_tofproton_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_tofproton_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_tofproton_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_tofproton_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_tofproton_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_tofproton_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_tofproton_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_tofproton_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_newpion_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_newpion_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_newpion_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_newpion_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_newpion_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_newpion_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_newpion_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_newpion_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_newkaon_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_newkaon_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_newkaon_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_newkaon_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_newkaon_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_newkaon_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_newkaon_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_newkaon_RESCALED", "Q_TPC_L_ch_RESCALED"},
  {"u_TPC_L_newproton_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_newproton_PLAIN", "Q_TPC_L_ch_PLAIN"},
  {"u_TPC_L_newproton_RECENTERED", "Q_TPC_R_ch_RECENTERED"},
  {"u_TPC_R_newproton_RECENTERED", "Q_TPC_L_ch_RECENTERED"},
  {"u_TPC_L_newproton_TWIST", "Q_TPC_R_ch_TWIST"},
  {"u_TPC_R_newproton_TWIST", "Q_TPC_L_ch_TWIST"},
  {"u_TPC_L_newproton_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_newproton_RESCALED", "Q_TPC_L_ch_RESCALED"}
};

vector <string> uv= 
  {"u_TPC_L_ch_RESCALED", "u_TPC_R_ch_RESCALED"};

vector<vector<string>> uvQvqv=
{
  {"u_TPC_F_ch_RESCALED", "Q_TPC_F_ch_RESCALED", "u_TPC_F_ch_RESCALED"}
};

void correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
	TStopwatch timer;
	timer.Start();
  int nSamples = 50;
  Qn::AxisD centAxis({"evCent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.}});
  Qn::AxisD vtxAxis({"evVtxZ", {-30.,-24.,-18.,-12.,-6.,0.,6.,12.,18.,24.,30.}});
  auto axes_correlation = Qn::MakeAxes(centAxis);//, vtxAxis);
  TChain *c=makeChain(inputFiles, "tree");
  ROOT::RDataFrame d(*c);
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  namespace P4 = Qn::Correlation::FourParticle;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wUnity1P = [](const Qn::QVector &a) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };
  auto wSumWu1P = [](const Qn::QVector &a) { return a.sumweights(); };
  auto wSumWu3P = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return a.sumweights(); };
  auto wSumWuEP = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights()/b.sumweights(); };
  auto wDenomASumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return 1./a.sumweights(); };
  auto wDenomBSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return 1./b.sumweights(); };
  auto wDenomABSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return 1./(a.sumweights()*b.sumweights()); };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};
	for (auto &corr:Qv)
	{
	  std::array<std::string, 1> qn{corr};
	  string corrName=corr;
	  corrBuilder.AddCorrelationWithInternalReader(corrName+"_v2c2_", P2::c2(2), wUnity1P, wn, qn, qn);
	  corrBuilder.AddCorrelationWithInternalReader(corrName+"_v2c4_", P4::c4(2), wUnity1P, wn, qn, qn);
	}
  for (auto &corr:QvQv)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1);

    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2X2_", P2::xx(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2Y2_", P2::yy(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2Y2_", P2::xy(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2X2_", P2::yx(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP2_", P2::ScalarProduct(2, 2), wUnity, wn, qn, qn);
  
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3X3_", P2::xx(3, 3), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3Y3_", P2::yy(3, 3), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3Y3_", P2::xy(3, 3), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3X3_", P2::yx(3, 3), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP3_", P2::ScalarProduct(3, 3), wUnity, wn, qn, qn);
  }
  for (auto &corr:uvQv)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1); 

    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2X2_", P2::xx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2Y2_", P2::yy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2Y2_", P2::xy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2X2_", P2::yx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP2_", P2::ScalarProduct(2, 2), wSumWu, wy, qn, qn);
  
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3X3_", P2::xx(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3Y3_", P2::yy(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3Y3_", P2::xy(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3X3_", P2::yx(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP3_", P2::ScalarProduct(3, 3), wSumWu, wy, qn, qn);
  }

  for (auto &corr:uvQvqv)
  {
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"_"+corr.at(1)+"_"+corr.at(2);

    corrBulder.AddCorrelationWithInternalReader(corrName+"_v2d2_", P2::d2(2), wSumWu3P, wy, qn, qn);
    corrBulder.AddCorrelationWithInternalReader(corrName+"_v2d4_", P2::d4(2), wSumWu3P, wy, qn, qn);
  }

  // ---------------- //
  // saving to output //
  // ---------------- //
  auto corrFile = TFile::Open(outputFile.c_str(), "RECREATE");
  corrFile->cd();
  auto results = corrBuilder.GetResults();
  for (auto &res : results) {
    res->Write();
  }
  corrFile->Close();
	timer.Stop();
	timer.Print();
}
