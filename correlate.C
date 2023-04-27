#include "/mnt/pool/nica/7/parfenovpeter/Soft/QnTools/install/include/QnTools/QnDataFrame.hpp"
//#include "QnDataFrame.hpp"
#include "utils.h"

vector <vector<string>> QvQv=
{
  {"Q_TPC_L_ch_RESCALED", "Q_TPC_R_ch_RESCALED"}
};

vector <string> Qv=
  {"Q_TPC_L_ch_RESCALED", "Q_TPC_R_ch_RESCALED"};

vector <vector<string>> uvQv=
{
  {"u_TPC_L_ch_RESCALED", "Q_TPC_R_ch_RESCALED"},
  {"u_TPC_R_ch_RESCALED", "Q_TPC_L_ch_RESCALED"},
};

vector <string> uv= 
  {"u_TPC_L_ch_PLAIN", "u_TPC_R_ch_PLAIN"};

void correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
	TStopwatch timer;
	timer.Start();
  int nSamples = 50;
  Qn::AxisD centAxis({"evCent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.}});
  //Qn::AxisD vtxAxis({"evVtxZ", 7, -70., 70.});
	Qn::AxisD ptAxis({"trPt", {0., 0.2, 0.4, 0.6, 0.8, 1., 1.5, 2., 3.}});
	Qn::AxisD etaAxis({"trEta", 5, -1., 1.});
  auto axes_correlation = Qn::MakeAxes(centAxis); //, vtxAxis);
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
	  //corrBuilder.AddCorrelationWithInternalReader(corrName+"_c2_", P2::c2(2), wUnity1P, wn, qn, qn);
	  //corrBuilder.AddCorrelationWithInternalReader(corrName+"_c2nom_", P2::c2nom(2), wUnity1P, wn, qn, qn);
	  //corrBuilder.AddCorrelationWithInternalReader(corrName+"_c4_", P4::c4(2), wUnity1P, wn, qn, qn);
	  //corrBuilder.AddCorrelationWithInternalReader(corrName+"_c4nom_", P4::c4nom(2), wUnity1P, wn, qn, qn);
	  //corrBuilder.AddCorrelationWithInternalReader(corrName+"_n2_", P2::n2(), wUnity1P, wn, qn, qn);
	  //corrBuilder.AddCorrelationWithInternalReader(corrName+"_n4_", P4::n4(), wUnity1P, wn, qn, qn);
	}
  for (auto &corr:QvQv)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XX_", P2::xx(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YY_", P2::yy(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XY_", P2::xy(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YX_", P2::yx(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP_", P2::ScalarProduct(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_coscos_", P2::xx(2, 2), wDenomABSumWu, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_sinsin_", P2::yy(2, 2), wDenomABSumWu, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_sinsin_", P2::xy(2, 2), wDenomABSumWu, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_coscos_", P2::yx(2, 2), wDenomABSumWu, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_EP_", P2::ScalarProduct(2, 2), wDenomABSumWu, wn, qn, qn);
  }
  for (auto &corr:uvQv)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1); 
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XX_", P2::xx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YY_", P2::yy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_XY_", P2::xy(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_YX_", P2::yx(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP_", P2::ScalarProduct(2, 2), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Xcos_", P2::xx(2, 2), wDenomBSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Ysin_", P2::yy(2, 2), wDenomBSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Xsin_", P2::xy(2, 2), wDenomBSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Ycos_", P2::yx(2, 2), wDenomBSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_EP_", P2::ScalarProduct(2, 2), wDenomBSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_d2_", P2::d2(2), wSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_d4_", P4::d4(2), wSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_d2_", P2::d2(2), wUnity, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_d2nom_", P2::d2nom(2), wUnity, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_d4_", P4::d4(2), wUnity, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_d4nom_", P4::d4nom(2), wUnity, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_nd2_", P2::nd2(), wSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_nd4_", P4::nd4(), wSumWu, wy, qn, qn);
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
