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

  // Define function to get pdg code
  auto getPdg = [](const ROOT::VecOps::RVec<int> &code) {
    auto pid = ROOT::VecOps::Map(code, [](int x){ return (double)x; });
    return pid;
  };

  // Define charge from pdg code
  auto getCharge = [](const ROOT::VecOps::EVec<int> &code) {
    auto pid = ROOT::VecOps::Map(code, [](int x){ return (double)x; });
    TParticlePDG *particle = (TParticlePDG*) TDatabasePDG::Instance()->GetParticle(pid);
    if (!particle) return -999.;
    auto charge = particle->Charge()/3.;
    delete particle;
    return charge;
  };
  /////////////////////////////////////////////////////////////////////////////////////////////////

  auto dd=d
    .Define("evB", [](float b){ return (float)b;}, {"bimp"})
    .Define("trPt",  getPt,  {"momx", "momy"})
    .Define("trEta", getEta, {"momx", "momy", "momz"})
    .Define("trPhi", getPhi, {"momx", "momy"})
    .Define("trPdg", getPdg, {"pdg"})
    .Define("trCharge", getCharge, {"pdg"})
    .Filter("evB<=16.")
    .Filter("trCharge!=0 && trCharge!=-999.")
  ;
  
  varPatterns=
  {
    "evB",               // kEvent
    "",                  // kChannel 
    "",                  // kRecParticle  
    "tr(Pt|Eta|Phi|Pdg)" // kSimParticle  
  };

  return dd; 
}

void setupQvectors()
{
  vector<Qn::AxisD> corrAxesEvent=
  {
    {"evB", 4,0,10},
  };

  vector<Qn::AxisD> corrAxesParticle=
  {
    {"trPt",4,0,2},
    {"trEta",4,0,2},
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

  name = "tr_TPC_F_u_Prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return abs(eta)<1.5;}, "eta_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
//  man.AddHisto2D("tr", {{"trEta", 100, 0., 6.}, {"trPt",  100, 0., 3.}}, "Ones");
  
  name = "tr_TPC_L_u_Prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>-1.5 && eta<-0.05);}, "eta_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "tr_TPC_R_u_Prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta<1.5 && eta>0.05);}, "eta_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "tr_TPC_F_Q_Prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesEvent, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return abs(eta)<1.5;}, "eta_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "tr_TPC_L_Q_Prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesEvent, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta>-1.5 && eta<-0.05);}, "eta_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");

  name = "tr_TPC_R_Q_Prot";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesEvent, {1,2,3,4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"trPdg"}, [](float pdg){return int(pdg)==2212;}, "proton");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){return (eta<1.5 && eta>0.05);}, "eta_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
}
