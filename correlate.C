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

vector <string> Qv=
  {"Q_TPC_L_ch_RESCALED", "Q_TPC_R_ch_RESCALED"};

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
  {"u_TPC_R_proton_RESCALED", "Q_TPC_L_ch_RESCALED"}
};

vector <string> uv= 
  {"u_TPC_L_ch_RESCALED", "u_TPC_R_ch_RESCALED"};

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

    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2X2_", P2::xx(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2Y2_", P2::yy(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2Y2_", P2::xy(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2X2_", P2::yx(2, 2), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP2_", P2::ScalarProduct(2, 2), wUnity, wn, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_cos2cos2_", P2::xx(2, 2), wDenomABSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_sin2sin2_", P2::yy(2, 2), wDenomABSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_sin2sin2_", P2::xy(2, 2), wDenomABSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_cos2cos2_", P2::yx(2, 2), wDenomABSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_EP2_", P2::ScalarProduct(2, 2), wDenomABSumWu, wy, qn, qn);

    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3X3_", P2::xx(3, 3), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3Y3_", P2::yy(3, 3), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3Y3_", P2::xy(3, 3), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3X3_", P2::yx(3, 3), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP3_", P2::ScalarProduct(3, 3), wUnity, wn, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_cos3cos3_", P2::xx(3, 3), wDenomABSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_sin3sin3_", P2::yy(3, 3), wDenomABSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_sin3sin3_", P2::xy(3, 3), wDenomABSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_cos3cos3_", P2::yx(3, 3), wDenomABSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_EP3_", P2::ScalarProduct(3, 3), wDenomABSumWu, wy, qn, qn);

    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_X4X4_", P2::xx(4, 4), wUnity, wn, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y4Y4_", P2::yy(4, 4), wUnity, wn, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_X4Y4_", P2::xy(4, 4), wUnity, wn, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y4X4_", P2::yx(4, 4), wUnity, wn, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP4_", P2::ScalarProduct(4, 4), wUnity, wn, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_cos4cos4_", P2::xx(4, 4), wDenomABSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_sin4sin4_", P2::yy(4, 4), wDenomABSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_sin4sin4_", P2::xy(4, 4), wDenomABSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_cos4cos4_", P2::yx(4, 4), wDenomABSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_EP4_", P2::ScalarProduct(4, 4), wDenomABSumWu, wy, qn, qn);
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
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2cos2_", P2::xx(2, 2), wDenomBSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2sin2_", P2::yy(2, 2), wDenomBSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2sin2_", P2::xy(2, 2), wDenomBSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2cos2_", P2::yx(2, 2), wDenomBSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_EP2_", P2::ScalarProduct(2, 2), wDenomBSumWu, wy, qn, qn);

    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3X3_", P2::xx(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3Y3_", P2::yy(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3Y3_", P2::xy(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3X3_", P2::yx(3, 3), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP3_", P2::ScalarProduct(3, 3), wSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3cos3_", P2::xx(3, 3), wDenomBSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3sin3_", P2::yy(3, 3), wDenomBSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_X3sin3_", P2::xy(3, 3), wDenomBSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y3cos3_", P2::yx(3, 3), wDenomBSumWu, wy, qn, qn);
    // corrBuilder.AddCorrelationWithInternalReader(corrName+"_EP3_", P2::ScalarProduct(3, 3), wDenomBSumWu, wy, qn, qn);

    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_X4X4_", P2::xx(4, 4), wSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y4Y4_", P2::yy(4, 4), wSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_X4Y4_", P2::xy(4, 4), wSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y4X4_", P2::yx(4, 4), wSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP4_", P2::ScalarProduct(4, 4), wSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_X4cos4_", P2::xx(4, 4), wDenomBSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y4sin4_", P2::yy(4, 4), wDenomBSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_X4sin4_", P2::xy(4, 4), wDenomBSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y4cos4_", P2::yx(4, 4), wDenomBSumWu, wy, qn, qn);
    //corrBuilder.AddCorrelationWithInternalReader(corrName+"_EP4_", P2::ScalarProduct(4, 4), wDenomBSumWu, wy, qn, qn);

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
