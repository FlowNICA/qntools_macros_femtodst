#include "utils.h"
#include "makeQvectors.h"

filteredDF defineVariables(definedDF &d);
void setupQvectors();

void makeQvectors(string inputFiles="/eos/nica/mpd/sim/taranenko/jam109/auau_4.5gev_rqmdrmf_MD2/part1/jam_auau_4.5gev_rqmdrmf_MD2_3836479_481.mcpico.root", string calibFilePath="qa.root", string outFilePath="qn.root")
{
  TStopwatch timer;
  timer.Start();
  std::string treename = "mctree";
  ROOT::RDataFrame d(*makeChain(inputFiles, treename.c_str()));
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
  // Define function for a centrality
  auto getCentB = [](float b) {
    // Hard coded centrality defenition
    // based on the impact parameter
    float fcent;
    if      (b < 2.91)  fcent = 2.5; // 0-5%
    else if (b < 4.18)  fcent = 7.5; // 5-10%
    else if (b < 6.01)  fcent = 15.; // 10-20%
    else if (b < 7.37)  fcent = 25.; // 20-30%
    else if (b < 8.52)  fcent = 35.; // 30-40%
    else if (b < 9.57)  fcent = 45.; // 40-50%
    else if (b < 10.55) fcent = 55.; // 50-60%
    else if (b < 11.46) fcent = 65.; // 60-70%
    else if (b < 12.31) fcent = 75.; // 70-80%
    else                fcent = -1;
    return fcent;
  };

  // Define function to get pt from px, py
  auto getPt = [](const ROOT::VecOps::RVec<Float_t> &px, const ROOT::VecOps::RVec<Float_t> &py) {
    auto pt = ROOT::VecOps::Map(px, py, [](float x, float y){ return sqrt(x*x + y*y); });
    return pt;
  };

  // Define function to get eta from px, py, pz
  auto getEta = [](const ROOT::VecOps::RVec<Float_t> &px, const ROOT::VecOps::RVec<Float_t> &py, const ROOT::VecOps::RVec<Float_t> &pz) {
    auto eta = ROOT::VecOps::Map(px, py, pz, [](float x, float y, float z){ 
      auto mod = sqrt(x*x + y*y + z*z);
      if (mod == z) return -999.;
      return 0.5 * log( (mod + z)/(mod - z) );
    });
    return eta;
  };

  // Define function to get azimuthal angle from px, py
  auto getPhi = [](const ROOT::VecOps::RVec<Float_t> &px, const ROOT::VecOps::RVec<Float_t> &py) {
    auto phi = ROOT::VecOps::Map(px, py, [](float x, float y){ return atan2(y,x); });
    return phi;
  };

  // Define function to do a simple 1 = 1 translation into a newly defined column for array of int
  auto getArrayValueI = [](const ROOT::VecOps::RVec<int> &val) {
    auto result = ROOT::VecOps::Map(val, [](int x){ return (double)x; });
    return result;
  };

  // Define function to do a simple 1 = 1 translation into a newly defined column for array of short
  auto getArrayValueSh = [](const ROOT::VecOps::RVec<short> &val) {
    auto result = ROOT::VecOps::Map(val, [](short x){ return (double)x; });
    return result;
  };
  /////////////////////////////////////////////////////////////////////////////////////////////////

  auto dd=d
    .Define("evCent",getCentB, {"bimp"})
    .Define("trPt",  getPt,  {"momx", "momy"})
    .Define("trEta", getEta, {"momx", "momy", "momz"})
    .Define("trPhi", getPhi, {"momx", "momy"})
    .Define("trPdg", getArrayValueI, {"pdg"})
    .Define("trCh",  getArrayValueSh, {"charge"})
    .Filter("evCent!=-1.") // At least one filter should be present (even if it always returns true)!!!
  ;
  
  varPatterns=
  {
    "evCent",                   // kEvent
    "",                         // kChannel 
    "",                         // kRecParticle  
    "tr(Pt|Eta|Phi|Pdg|Ch)"     // kSimParticle  
  };

  return dd; 
}

void setupQvectors()
{
  vector<Qn::AxisD> corrAxesEvent=
  {
    {"evCent", {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.}},
  };

  vector<Qn::AxisD> corrAxesParticle=
  {
    {"trPt",30,0,3},
    {"trEta",30,-1.5,1.5},
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
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return abs(eta)<1.5;}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return pt<3.;}, "pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"},  [](float ch){return (ch!=0. && ch!=-999.);}, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
//  man.AddHisto2D("tr", {{"trEta", 100, 0., 6.}, {"trPt",  100, 0., 3.}}, "Ones");
  
  name = "u_TPC_L_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>-1.5 && eta<-0.05);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"},  [](float ch){return (ch!=0. && ch!=-999.);}, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "u_TPC_R_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta<1.5 && eta>0.05);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"},  [](float ch){return (ch!=0. && ch!=-999.);}, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "u_TPC_F_prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return abs(eta)<1.5;}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return pt<3.;}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
//  man.AddHisto2D("tr", {{"trEta", 100, 0., 6.}, {"trPt",  100, 0., 3.}}, "Ones");
  
  name = "u_TPC_L_prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>-1.5 && eta<-0.05);}, "eta_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "u_TPC_R_prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta<1.5 && eta>0.05);}, "eta_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_F_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return abs(eta)<1.5;}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<3. && pt>0.2);}, "pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"},  [](float ch){return (ch!=0. && ch!=-999.);}, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_L_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>-1.5 && eta<-0.05);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<3. && pt>0.2);}, "pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"},  [](float ch){return (ch!=0. && ch!=-999.);}, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_R_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta<1.5 && eta>0.05);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<3. && pt>0.2);}, "pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"},  [](float ch){return (ch!=0. && ch!=-999.);}, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_F_prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return abs(eta)<1.5;}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<3. && pt>0.2);}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_L_prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>-1.5 && eta<-0.05);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<3. && pt>0.2);}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "Q_TPC_R_prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta<1.5 && eta>0.05);}, "eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"},  [](float pt){return (pt<3. && pt>0.2);}, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
}
