#include "/mnt/pool/nica/7/parfenovpeter/Soft/QnTools/install/include/QnTools/QnDataFrame.hpp"
//#include "QnDataFrame.hpp"
#include "utils.h"

vector <vector<string>> QvQv=
{
  {"Q_TPC_L_ch_PLAIN", "Q_TPC_R_ch_PLAIN"}
};

vector <string> Qv=
  {"Q_TPC_L_ch_PLAIN", "Q_TPC_R_ch_PLAIN", "Q_TPC_F_ch_PLAIN"};

vector <vector<string>> uvQv=
{
  {"u_TPC_L_ch_PLAIN", "Q_TPC_R_ch_PLAIN"},
  {"u_TPC_R_ch_PLAIN", "Q_TPC_L_ch_PLAIN"}
};

vector <string> uv= 
  {"u_TPC_L_ch_PLAIN", "u_TPC_R_ch_PLAIN", "u_TPC_F_ch_PLAIN"};

void correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
	TStopwatch timer;
	timer.Start();
  int nSamples = 50;
  Qn::AxisD centAxis({"evCent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.}});
  Qn::AxisD vtxAxis({"evVtxZ", 7, -70., 70.});
  auto axes_correlation = Qn::MakeAxes(centAxis, vtxAxis);
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

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};
	for (auto &corr:Qv)
	{
	  std::array<std::string, 1> qn{corr};
	  string corrName=corr;
	  corrBuilder.AddCorrelationWithInternalReader(corrName+"_c2_", P2::c2(2), wUnity1P, wn, qn, qn);
	  corrBuilder.AddCorrelationWithInternalReader(corrName+"_c4_", P4::c4(2), wUnity1P, wn, qn, qn);
	}
  for (auto &corr:QvQv)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XX_", P2::xx(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YY_", P2::yy(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XY_", P2::xy(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YX_", P2::yx(2, 2), wUnity, wn, qn, qn);
  }
  for (auto &corr:uvQv)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1); 
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XX_", P2::xx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YY_", P2::yy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XY_", P2::xy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YX_", P2::yx(2, 2), wSumWu, wy, qn, qn);
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
