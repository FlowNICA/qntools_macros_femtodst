#include "utils.h"
#include "makeQvectors.h"

filteredDF defineVariables(definedDF &d);
void setupQvectors();

void makeQvectors(string inputFiles="/eos/nica/mpd/sim/taranenko/jam109/auau_4.5gev_rqmdrmf_MD2/part1/jam_auau_4.5gev_rqmdrmf_MD2_3836479_481.mcpico.root", string calibFilePath="qa.root", string outFilePath="qn.root")
{
  TStopwatch timer;
  timer.Start();

  ROOT::RDataFrame d(*makeChain(inputFiles, "tStarData"));
  auto dd=defineVariables(d);
  init(dd, outFilePath, calibFilePath);
  setupQvectors(); 
  timer.Stop();
  timer.Print();
  timer.Start();
  run(dd);
  cout << "Done!\n";
  timer.Stop();
  timer.Print();
}

filteredDF defineVariables(definedDF &d)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Define function to get pt from px, py
  auto getPt = [](const ROOT::VecOps::RVec<Double_t> &px, const ROOT::VecOps::RVec<Double_t> &py) {
    auto pt = ROOT::VecOps::Map(px, py, [](float x, float y){ return sqrt(x*x + y*y); });
    return pt;
  };

  // Define function to get eta from px, py, pz
  auto getEta = [](const ROOT::VecOps::RVec<Double_t> &px, const ROOT::VecOps::RVec<Double_t> &py, const ROOT::VecOps::RVec<Double_t> &pz) {
    auto eta = ROOT::VecOps::Map(px, py, pz, [](float x, float y, float z){ 
      auto mod = sqrt(x*x + y*y + z*z);
      if (mod == z) return -999.;
      return 0.5 * log( (mod + z)/(mod - z) );
    });
    return eta;
  };

  // Define function to get azimuthal angle from px, py
  auto getPhi = [](const ROOT::VecOps::RVec<Double_t> &px, const ROOT::VecOps::RVec<Double_t> &py) {
    auto phi = ROOT::VecOps::Map(px, py, [](float x, float y){ return atan2(y,x); });
    return phi;
  };
	/////////////////////////////////////////////////////////////////////////////////////////////////

  auto dd=d
		.Define("evCent","ev_cent")
		.Define("evVtxX","ev_vtxX")
		.Define("evVtxY","ev_vtxY")
		.Define("evVtxZ","ev_vtxZ")
		.Define("evVpdZ","ev_vpdZ")
    .Define("trPt",  getPt,  {"tr_px", "tr_py"})
    .Define("trEta", getEta, {"tr_px", "tr_py", "tr_pz"})
    .Define("trPhi", getPhi, {"tr_px", "tr_py"})
    .Define("trCh",  "tr_ch")
    .Define("trNhits", "tr_Nhits")
    .Define("trNhitsFit", "tr_NhitsFit")
    .Define("trNhitsPoss", "tr_NhitsPoss")
    .Define("trDca", "tr_dca")
    .Filter("evVtxZ>-70 && evVtxZ<70") // At least one filter should be present (even if it always returns true)!!!
  ;
  
  varPatterns=
  {
    "ev(Cent|VtxX|VtxY|VtxZ|VpdZ)", // kEvent
    "",                         // kChannel 
    "tr(Pt|Eta|Phi|Ch|Nhits|NhitsFit|NhitsPoss|Dca)",  // kRecParticle  
    ""     // kSimParticle  
  };

  return dd; 
}

void setupQvectors()
{
  vector<Qn::AxisD> corrAxesEvent=
  {
    {"evCent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.}},
  };

  vector<Qn::AxisD> corrAxesParticle=
  {
    {"trPt",30,0,3},
    {"trEta",20,-1.,1.},
  };

  for (auto &axis:corrAxesEvent)
    man.AddCorrectionAxis(axis);

  Qn::Recentering recentering;
  recentering.SetApplyWidthEqualization(false);
  Qn::TwistAndRescale twistRescale;
  twistRescale.SetApplyRescale(true);
  twistRescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::DOUBLE_HARMONIC);
  
  auto sumW=Qn::QVector::Normalization::M;
  auto track=Qn::DetectorType::TRACK;
  auto channel=Qn::DetectorType::CHANNEL;
  auto plain=Qn::QVector::CorrectionStep::PLAIN;
  auto recentered=Qn::QVector::CorrectionStep::RECENTERED;
  auto twisted=Qn::QVector::CorrectionStep::TWIST;
  auto rescaled=Qn::QVector::CorrectionStep::RESCALED;
  
  std::string name;

  name = "u_TPC_F_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return abs(eta)<1.;}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return pt<3.;}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
//  man.AddHisto2D("tr", {{"trEta", 100, 0., 6.}, {"trPt",  100, 0., 3.}}, "Ones");

  name = "u_TPC_L_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>-1. && eta<-0.05);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return pt<3.;}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "u_TPC_R_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>0.05 && eta<1.);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return pt<3.;}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_F_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return abs(eta)<1.;}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<2. && pt>0.2);}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_L_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>-1. && eta<-0.05);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<2. && pt>0.2);}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_R_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>0.05 && eta<1.);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<2. && pt>0.2);}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
}
