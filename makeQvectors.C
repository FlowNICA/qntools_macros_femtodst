#include "utils.h"
#include "makeQvectors.h"
#include <random>

filteredDF defineVariables(definedDF &d);
void setupQvectors();

void makeQvectors(string inputFiles = "/eos/nica/mpd/sim/taranenko/jam109/auau_4.5gev_rqmdrmf_MD2/part1/jam_auau_4.5gev_rqmdrmf_MD2_3836479_481.mcpico.root", string calibFilePath = "qa.root", string outFilePath = "qn.root")
{
  TStopwatch timer;
  timer.Start();

  ROOT::RDataFrame d(*makeChain(inputFiles, "tStarData"));
  auto dd = defineVariables(d);
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
  // Define function to get p from px, py, pz
  auto getP = [](const RVecD &px, const RVecD &py, const RVecD &pz){
    auto p = Map(px, py, pz, [](double x, double y, double z){
        return (double)sqrt(x * x + y * y + z * z); });
    return p;
  };

  // Define function to get pt from px, py
  auto getPt = [](const RVecD &px, const RVecD &py){
    auto pt = Map(px, py, [](double x, double y){
        return (double)sqrt(x * x + y * y); });
    return pt;
  };

  // Define function to get eta from px, py, pz
  auto getEta = [](const RVecD &px, const RVecD &py, const RVecD &pz){
    auto eta = Map(px, py, pz, [](double x, double y, double z){ 
        auto mod = (double)sqrt(x*x + y*y + z*z);
        if (mod == z) return -999.;
        return (double)(0.5 * log( (mod + z)/(mod - z) )); });
    return eta;
  };

  // Define function to get azimuthal angle from px, py
  auto getPhi = [](const RVecD &px, const RVecD &py){
    auto phi = Map(px, py, [](double x, double y){
        return (double)atan2(y, x); });
    return phi;
  };

  // Assign a random number from 1 to 4
  auto rndNum = [](const RVecD &px){
    auto result = Map(px, [](double x){ 
        // Prepare variables for a random sub-event determination
        std::random_device rd; 
        std::mt19937 gen(rd()); 
        std::uniform_int_distribution<> distr(1,4); // 4 sub-events division
        return (double) distr(gen); });
    return result;
  };

  // PID functions
  //  is it pion?
  auto pidPion = [](const RVecI &isTof, const RVecD &p, const RVecF &m2,
                    const RVecF &nsigPi, const RVecF &nsigKa){
    auto result = Map(isTof, p, m2, nsigPi, nsigKa, [](int _isTof, double _p, float _m2, float _nsigPi, float _nsigKa){
        int ispion = 0;
        if (_isTof!=1 && abs(_nsigPi) < 2. && _p > 0.2 && _p < 0.6)
          ispion = 1;
        if (_isTof!=1 && abs(_nsigPi) < 2. && abs(_nsigKa) > 2. && _p > 0.6 && _p < 0.7)
          ispion = 1;
        if (_isTof==1 && abs(_nsigPi) < 3. && _m2 > -0.15 && _m2 < 0.1 && _p > 0.2 && _p < 3.4)
          ispion = 1;
        return ispion;	
      });
    return result;
  };

  //  is it kaon?
  auto pidKaon = [](const RVecI &isTof, const RVecD &p, const RVecF &m2,
                    const RVecF &nsigPi, const RVecF &nsigKa){
    auto result = Map(isTof, p, m2, nsigPi, nsigKa, [](int _isTof, double _p, float _m2, float _nsigPi, float _nsigKa){
        int iskaon = 0;
        if (_isTof!=1 && abs(_nsigKa) < 2 && _p > 0.2 && _p < 0.5)
          iskaon = 1;
        if (_isTof!=1 && abs(_nsigKa) < 2. && abs(_nsigPi) > 3. && _p > 0.5 && _p < 0.7)
          iskaon = 1;
        if (_isTof==1 && abs(_nsigKa) < 3. && _m2 > 0.2 && _m2 < 0.32 && _p > 0.2 && _p < 3.4)
          iskaon = 1;
        return iskaon;	
      });
    return result;
  };

  //  is it proton?
  auto pidProton = [](const RVecI &isTof, const RVecD &p, const RVecF &m2,
                      const RVecF &nsigPi, const RVecF &nsigPr){
    auto result = Map(isTof, p, m2, nsigPi, nsigPr, [](int _isTof, double _p, float _m2, float _nsigPi, float _nsigPr){
        int isproton = 0;
        if (_isTof!=1 && abs(_nsigPr) < 2 && _p > 0.4 && _p < 0.9)
          isproton = 1;
        if (_isTof!=1 && abs(_nsigPr) < 2. && abs(_nsigPi) > 3. && _p > 0.9 && _p < 1.2)
          isproton = 1;
        if (_isTof==1 && abs(_nsigPr) < 3. && _m2 > 0.75 && _m2 < 1.2 && _p > 0.2 && _p < 3.4)
          isproton = 1; 
        return isproton;
      });
    return result;
  };

  // NewPID (STAR-like)
  // is it pion?
  auto newpidPion = [](const RVecF &pidx, const RVecF &pidy, const RVecI &isTof, const RVecF &m2,
                       const RVecF &meanXpi, const RVecF &meanYpi, const RVecF &sigmXpi, const RVecF &sigmYpi,
                       const RVecF &meanXka, const RVecF &meanYka, const RVecF &sigmXka, const RVecF &sigmYka,
                       const RVecF &meanXpr, const RVecF &meanYpr, const RVecF &sigmXpr, const RVecF &sigmYpr){
    auto result = Map(pidx, pidy, isTof, m2, meanXpi, meanYpi, sigmXpi, sigmYpi, meanXka, meanYka, sigmXka, sigmYka, meanXpr, meanYpr, sigmXpr, sigmYpr,
                      [](float _pidx, float _pidy, int _isTof, float _m2, float _meanXpi, float _meanYpi, float _sigmXpi, float _sigmYpi, float _meanXka, float _meanYka, float _sigmXka, float _sigmYka, float _meanXpr, float _meanYpr, float _sigmXpr, float _sigmYpr){
                        int ispion = 0;
                        if ( pow(_pidx - _meanXpi,2)/pow(3.*_sigmXpi,2)+pow(_pidy - _meanYpi,2)/pow(3.*_sigmYpi,2) < 1. && _isTof==1 && _m2 > -0.15 && _m2 < 0.1)
                        /*if ( pow(_pidx - _meanXpi,2)/pow(3.*_sigmXpi,2)+pow(_pidy - _meanYpi,2)/pow(3.*_sigmYpi,2) < 1. &&
                             pow(_pidx - _meanXka,2)/pow(3.*_sigmXka,2)+pow(_pidy - _meanYka,2)/pow(3.*_sigmYka,2) > 1. &&
                             pow(_pidx - _meanXpr,2)/pow(3.*_sigmXpr,2)+pow(_pidy - _meanYpr,2)/pow(3.*_sigmYpr,2) > 1.)*/
                          ispion = 1;
                        return ispion;
                      });
    return result;
  };

  // is it kaon?
  auto newpidKaon = [](const RVecF &pidx, const RVecF &pidy, const RVecI &isTof, const RVecF &m2,
                       const RVecF &meanXpi, const RVecF &meanYpi, const RVecF &sigmXpi, const RVecF &sigmYpi,
                       const RVecF &meanXka, const RVecF &meanYka, const RVecF &sigmXka, const RVecF &sigmYka,
                       const RVecF &meanXpr, const RVecF &meanYpr, const RVecF &sigmXpr, const RVecF &sigmYpr){
    auto result = Map(pidx, pidy, isTof, m2, meanXpi, meanYpi, sigmXpi, sigmYpi, meanXka, meanYka, sigmXka, sigmYka, meanXpr, meanYpr, sigmXpr, sigmYpr,
                      [](float _pidx, float _pidy, int _isTof, float _m2, float _meanXpi, float _meanYpi, float _sigmXpi, float _sigmYpi, float _meanXka, float _meanYka, float _sigmXka, float _sigmYka, float _meanXpr, float _meanYpr, float _sigmXpr, float _sigmYpr){
                        int iskaon = 0;
                        if ( pow(_pidx - _meanXka,2)/pow(3.*_sigmXka,2)+pow(_pidy - _meanYka,2)/pow(3.*_sigmYka,2) < 1. && _isTof==1 && _m2 > 0.2 && _m2 < 0.32)
                        /*if ( pow(_pidx - _meanXpi,2)/pow(3.*_sigmXpi,2)+pow(_pidy - _meanYpi,2)/pow(3.*_sigmYpi,2) > 1. &&
                             pow(_pidx - _meanXka,2)/pow(3.*_sigmXka,2)+pow(_pidy - _meanYka,2)/pow(3.*_sigmYka,2) < 1. &&
                             pow(_pidx - _meanXpr,2)/pow(3.*_sigmXpr,2)+pow(_pidy - _meanYpr,2)/pow(3.*_sigmYpr,2) > 1.)*/
                          iskaon = 1;
                        return iskaon;
                      });
    return result;
  };

  // is it proton?
  auto newpidProton = [](const RVecF &pidx, const RVecF &pidy, const RVecI &isTof, const RVecF &m2,
                         const RVecF &meanXpi, const RVecF &meanYpi, const RVecF &sigmXpi, const RVecF &sigmYpi,
                         const RVecF &meanXka, const RVecF &meanYka, const RVecF &sigmXka, const RVecF &sigmYka,
                         const RVecF &meanXpr, const RVecF &meanYpr, const RVecF &sigmXpr, const RVecF &sigmYpr){
    auto result = Map(pidx, pidy, isTof, m2, meanXpi, meanYpi, sigmXpi, sigmYpi, meanXka, meanYka, sigmXka, sigmYka, meanXpr, meanYpr, sigmXpr, sigmYpr,
                      [](float _pidx, float _pidy, int _isTof, float _m2, float _meanXpi, float _meanYpi, float _sigmXpi, float _sigmYpi, float _meanXka, float _meanYka, float _sigmXka, float _sigmYka, float _meanXpr, float _meanYpr, float _sigmXpr, float _sigmYpr){
                        int isproton = 0;
                        if ( pow(_pidx - _meanXpr,2)/pow(3.*_sigmXpr,2)+pow(_pidy - _meanYpr,2)/pow(3.*_sigmYpr,2) < 1. && _isTof==1 && _m2 > 0.75 && _m2 < 1.2)
                        /*if ( pow(_pidx - _meanXpi,2)/pow(3.*_sigmXpi,2)+pow(_pidy - _meanYpi,2)/pow(3.*_sigmYpi,2) > 1. &&
                             pow(_pidx - _meanXka,2)/pow(3.*_sigmXka,2)+pow(_pidy - _meanYka,2)/pow(3.*_sigmYka,2) > 1. &&
                             pow(_pidx - _meanXpr,2)/pow(3.*_sigmXpr,2)+pow(_pidy - _meanYpr,2)/pow(3.*_sigmYpr,2) < 1.)*/
                          isproton = 1;
                        return isproton;
                      });
    return result;
  };
  /////////////////////////////////////////////////////////////////////////////////////////////////

  auto dd = d.Define("evCent", "ev_cent")
    .Define("evVtxX", "ev_vtxX")
    .Define("evVtxY", "ev_vtxY")
    .Define("evVtxZ", "ev_vtxZ")
    .Define("evVpdZ", "ev_vpdZ")
    .Define("trPt", getPt, {"tr_px", "tr_py"})
    .Define("trP", getP, {"tr_px", "tr_py", "tr_pz"})
    .Define("trEta", getEta, {"tr_px", "tr_py", "tr_pz"})
    .Define("trPhi", getPhi, {"tr_px", "tr_py"})
    .Define("trCh", "tr_ch")
    .Define("trNhits", "tr_Nhits")
    .Define("trNhitsFit", "tr_NhitsFit")
    .Define("trNhitsPoss", "tr_NhitsPoss")
    .Define("trDca", "tr_dca")
    .Define("trRnd4Sub", rndNum, {"tr_px"})
    .Define("trIsTof", "tr_isTofTrack")
    .Define("trM2", "tr_tofm2")
    .Define("trDedx", "tr_dedx")
    .Define("trNsigEl", "tr_nSigEl")
    .Define("trNsigPi", "tr_nSigPi")
    .Define("trNsigKa", "tr_nSigKa")
    .Define("trNsigPr", "tr_nSigPr")
    .Define("trPidPi", pidPion, {"trIsTof", "trP", "trM2", "trNsigPi", "trNsigKa"})
    .Define("trPidKa", pidKaon, {"trIsTof", "trP", "trM2", "trNsigPi", "trNsigKa"})
    .Define("trPidPr", pidProton, {"trIsTof", "trP", "trM2", "trNsigPi", "trNsigPr"})
    .Define("trNewPidX", "tr_pidX")
    .Define("trNewPidY", "tr_pidY")
    .Define("trNewPidMeanXPi", "tr_pidMeanXPi")
    .Define("trNewPidMeanYPi", "tr_pidMeanYPi")
    .Define("trNewPidSigmXPi", "tr_pidSigmXPi")
    .Define("trNewPidSigmYPi", "tr_pidSigmYPi")
    .Define("trNewPidMeanXKa", "tr_pidMeanXKa")
    .Define("trNewPidMeanYKa", "tr_pidMeanYKa")
    .Define("trNewPidSigmXKa", "tr_pidSigmXKa")
    .Define("trNewPidSigmYKa", "tr_pidSigmYKa")
    .Define("trNewPidMeanXPr", "tr_pidMeanXPr")
    .Define("trNewPidMeanYPr", "tr_pidMeanYPr")
    .Define("trNewPidSigmXPr", "tr_pidSigmXPr")
    .Define("trNewPidSigmYPr", "tr_pidSigmYPr")
    .Define("trNewPidPi", newpidPion, {"trNewPidX", "trNewPidY", "trIsTof", "trM2", "trNewPidMeanXPi", "trNewPidMeanYPi", "trNewPidSigmXPi", "trNewPidSigmYPi", "trNewPidMeanXKa", "trNewPidMeanYKa", "trNewPidSigmXKa", "trNewPidSigmYKa", "trNewPidMeanXPr", "trNewPidMeanYPr", "trNewPidSigmXPr", "trNewPidSigmYPr"})
    .Define("trNewPidKa", newpidKaon, {"trNewPidX", "trNewPidY", "trIsTof", "trM2", "trNewPidMeanXPi", "trNewPidMeanYPi", "trNewPidSigmXPi", "trNewPidSigmYPi", "trNewPidMeanXKa", "trNewPidMeanYKa", "trNewPidSigmXKa", "trNewPidSigmYKa", "trNewPidMeanXPr", "trNewPidMeanYPr", "trNewPidSigmXPr", "trNewPidSigmYPr"})
    .Define("trNewPidPr", newpidProton, {"trNewPidX", "trNewPidY", "trIsTof", "trM2", "trNewPidMeanXPi", "trNewPidMeanYPi", "trNewPidSigmXPi", "trNewPidSigmYPi", "trNewPidMeanXKa", "trNewPidMeanYKa", "trNewPidSigmXKa", "trNewPidSigmYKa", "trNewPidMeanXPr", "trNewPidMeanYPr", "trNewPidSigmXPr", "trNewPidSigmYPr"})
    .Filter("abs(evVtxZ)<30.")                      // At least one filter should be present (even if it always returns true)!!!
    .Filter("sqrt(evVtxX*evVtxX+evVtxY*evVtxY)<2.") // At least one filter should be present (even if it always returns true)!!!
    .Filter("abs(evVtxZ-evVpdZ)<3.")                // At least one filter should be present (even if it always returns true)!!!
    ;

  varPatterns =
    {
      "ev(Cent|VtxX|VtxY|VtxZ|VpdZ)",                                               // kEvent
      "",                                                                           // kChannel
      "tr(P|Pt|Eta|Phi|Ch|Nhits|NhitsFit|NhitsPoss|Dca|Rnd4Sub|PidPi|PidKa|PidPr|NewPidPi|NewPidKa|NewPidPr|M2|NsigPi|NsigKa|NsigPr|Dedx|NewPidX|NewPidY)", // kRecParticle
      ""                                                                            // kSimParticle
    };

  return dd;
}

void setupQvectors()
{
  vector<Qn::AxisD> corrAxesEvent =
    {
      {"evCent", {0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80.}},
      {"evVtxZ", {-30., -24., -18., -12., -6., 0., 6., 12., 18., 24., 30.}},
    };

  vector<Qn::AxisD> corrAxesParticle =
    {
      //{"trPt",50,0,5},
      {"trPt", {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.25, 3.5, 3.75, 4., 4.5, 5.}},
      //{"trEta", 20, -1., 1.},
      //{"trEta",4,-1.,1.},
    };

  for (auto &axis : corrAxesEvent)
    man.AddCorrectionAxis(axis);

  Qn::GainEqualization gain;
  gain.SetEqualizationMethod(Qn::GainEqualization::Method::AVERAGE);
  gain.SetNoOfEntriesThreshold(1);
  gain.SetUseChannelGroupsWeights(true);

  Qn::Recentering recentering;
  recentering.SetApplyWidthEqualization(false);
  Qn::TwistAndRescale twistRescale;
  twistRescale.SetApplyRescale(true);
  twistRescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::DOUBLE_HARMONIC);

  auto sumW = Qn::QVector::Normalization::M;
  auto track = Qn::DetectorType::TRACK;
  auto channel = Qn::DetectorType::CHANNEL;
  auto plain = Qn::QVector::CorrectionStep::PLAIN;
  auto recentered = Qn::QVector::CorrectionStep::RECENTERED;
  auto twisted = Qn::QVector::CorrectionStep::TWIST;
  auto rescaled = Qn::QVector::CorrectionStep::RESCALED;

  std::string name;

  name = "u_TPC_L_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > -1. && eta < -0.05); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 3.); },"dca_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");

  name = "u_TPC_L_pion";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > -1. && eta < -0.05); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trPidPi"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_L_kaon";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > -1. && eta < -0.05); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trPidKa"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_L_proton";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > -1. && eta < -0.05); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trPidPr"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_L_newpion";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > -1. && eta < -0.05); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trNewPidPi"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNewPidX", 520, -0.8, 1.8}, {"trNewPidY", 520, -1.8, 0.8}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_L_newkaon";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > -1. && eta < -0.05); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trNewPidKa"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNewPidX", 520, -0.8, 1.8}, {"trNewPidY", 520, -1.8, 0.8}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_L_newproton";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > -1. && eta < -0.05); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trNewPidPr"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNewPidX", 520, -0.8, 1.8}, {"trNewPidY", 520, -1.8, 0.8}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_R_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > 0.05 && eta < 1.); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((float)nhits / (float)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 3.); },"dca_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");

  name = "u_TPC_R_pion";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > 0.05 && eta < 1.); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((float)nhits / (float)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trPidPi"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_R_kaon";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > 0.05 && eta < 1.); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((float)nhits / (float)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trPidKa"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_R_proton";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > 0.05 && eta < 1.); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((float)nhits / (float)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trPidPr"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_R_newpion";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > 0.05 && eta < -1.); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trNewPidPi"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNewPidX", 520, -0.8, 1.8}, {"trNewPidY", 520, -1.8, 0.8}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_R_newkaon";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > 0.05 && eta < -1.); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trNewPidKa"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNewPidX", 520, -0.8, 1.8}, {"trNewPidY", 520, -1.8, 0.8}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");

  name = "u_TPC_R_newproton";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", corrAxesParticle, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > 0.05 && eta < -1.); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 5. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((double)nhits / (double)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 1.); },"dca_cut");
  man.AddCutOnDetector(name.c_str(), {"trNewPidPr"}, [](int pid){ return (pid == 1); },"pid_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNewPidX", 520, -0.8, 1.8}, {"trNewPidY", 520, -1.8, 0.8}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trM2", 170, -0.2, 1.5}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trP", 175, 0., 3.5}, {"trDedx", 325, 0., 6.5e-5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPi", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigKa", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trNsigPr", 330, -3.3, 3.3}, {"trP", 175, 0., 3.5}}, "Ones");
  
  name = "Q_TPC_L_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", {}, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > -1. && eta < -0.05); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 2. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((float)nhits / (float)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 3.); },"dca_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");

  name = "Q_TPC_R_ch";
  man.AddDetector(name.c_str(), track, "trPhi", "trPt", {}, {1, 2, 3, 4}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trEta"}, [](float eta){ return (eta > 0.05 && eta < 1.); },"eta_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt < 2. && pt > 0.15); },"pt_cut");
  man.AddCutOnDetector(name.c_str(), {"trP"}, [](float p){ return (p < 10.); },"p_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit"}, [](float nhits){ return (nhits > 15); },"nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsPoss"}, [](float nhitsposs){ return (nhitsposs > 0); },"nhitstrNhitsPoss_cut");
  man.AddCutOnDetector(name.c_str(), {"trNhitsFit", "trNhitsPoss"}, [](float nhits, float nhitsposs){ return ((float)nhits / (float)nhitsposs > 0.51); },"nhitsRatio_cut");
  man.AddCutOnDetector(name.c_str(), {"trDca"}, [](float dca){ return (abs(dca) < 3.); },"dca_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"evCent", 16, 0., 80.}, "Ones");
  man.AddHisto1D(name.c_str(), {"evVtxZ", 60, -30., 30.}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto1D(name.c_str(), {"trPt", 5000, 0., 5.}, "Ones");

  // man.AddCorrectionOnInputData(name.c_str(), gain);
}
