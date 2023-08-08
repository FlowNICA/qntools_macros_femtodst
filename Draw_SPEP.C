#include "Functions.C"

void Draw_SPEP(TString inFileName = "", TString outFileName = "./test_graphs.root", bool isV2 = true, bool isV3 = true)
{
  TFile *fiCorr = new TFile(inFileName.Data(), "read");
  if (inFileName == "" || !fiCorr)
  {
    std::cerr << "No input file was provided!" << std::endl;
    return;
  }
  else
  {
    std::cout << "Opened file: " << inFileName.Data() << std::endl;
  }

  // Set centrality selection
  const std::vector<std::pair<float, float>> cent_ranges_v2 = {{0., 5.}, {5., 10.}, {10., 20.}, {20., 30.}, {30., 40.}, {40., 60.}, {10., 40.}, {0., 30.}, {30., 80.}};
  const std::vector<std::pair<float, float>> cent_ranges_v3 = {{0., 10.}, {10., 20.}, {20., 30.}, {30., 40.}, {40., 50.}, {50., 60.}, {10., 40.}, {0., 30.}, {30., 80.}};

  Qn::AxisD axis_pt_v2_ch("trPt",     {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.8, 3.25, 3.75, 4.5});
  Qn::AxisD axis_pt_v2_pion("trPt",   {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.25});
  Qn::AxisD axis_pt_v2_kaon("trPt",   {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.25});
  Qn::AxisD axis_pt_v2_proton("trPt", {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.25});
  Qn::AxisD axis_pt_v3_ch("trPt",     {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.6, 2.8, 3., 3.25});
  Qn::AxisD axis_pt_v3_pion("trPt",   {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.6, 2.8, 3., 3.25});
  Qn::AxisD axis_pt_v3_kaon("trPt",   {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.6, 2.8, 3., 3.25});
  Qn::AxisD axis_pt_v3_proton("trPt", {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.6, 2.8, 3., 3.25});

  std::vector<std::string> pid_names = {"ch", "pion", "kaon", "proton"};

  // Names for <uQ> and <QQ> correlations for vn: {{uQ, QQ}, ...}
  std::vector<std::pair<std::string, std::string>> corr_names_v2 = {
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_L_pion_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_R_pion_RESCALED_Q_TPC_L_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_L_pion_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_R_pion_RESCALED_Q_TPC_L_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_L_pion_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_R_pion_RESCALED_Q_TPC_L_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_L_kaon_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_R_kaon_RESCALED_Q_TPC_L_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_L_kaon_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_R_kaon_RESCALED_Q_TPC_L_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_L_kaon_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_R_kaon_RESCALED_Q_TPC_L_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_L_proton_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_R_proton_RESCALED_Q_TPC_L_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_L_proton_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_R_proton_RESCALED_Q_TPC_L_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_L_proton_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_R_proton_RESCALED_Q_TPC_L_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_L_newpion_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_R_newpion_RESCALED_Q_TPC_L_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_L_newpion_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_R_newpion_RESCALED_Q_TPC_L_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_L_newpion_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_R_newpion_RESCALED_Q_TPC_L_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_L_newkaon_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_R_newkaon_RESCALED_Q_TPC_L_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_L_newkaon_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_R_newkaon_RESCALED_Q_TPC_L_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_L_newkaon_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_R_newkaon_RESCALED_Q_TPC_L_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_L_newproton_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_R_newproton_RESCALED_Q_TPC_L_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_L_newproton_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_R_newproton_RESCALED_Q_TPC_L_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_L_newproton_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_R_newproton_RESCALED_Q_TPC_L_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"}
  };
  std::vector<std::pair<std::string, std::string>> corr_names_v3 = {
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_L_pion_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_R_pion_RESCALED_Q_TPC_L_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_L_pion_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_R_pion_RESCALED_Q_TPC_L_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_L_pion_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_R_pion_RESCALED_Q_TPC_L_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_L_kaon_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_R_kaon_RESCALED_Q_TPC_L_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_L_kaon_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_R_kaon_RESCALED_Q_TPC_L_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_L_kaon_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_R_kaon_RESCALED_Q_TPC_L_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_L_proton_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_R_proton_RESCALED_Q_TPC_L_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_L_proton_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_R_proton_RESCALED_Q_TPC_L_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_L_proton_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_R_proton_RESCALED_Q_TPC_L_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_L_newpion_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_R_newpion_RESCALED_Q_TPC_L_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_L_newpion_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_R_newpion_RESCALED_Q_TPC_L_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_L_newpion_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_R_newpion_RESCALED_Q_TPC_L_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_L_newkaon_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_R_newkaon_RESCALED_Q_TPC_L_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_L_newkaon_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_R_newkaon_RESCALED_Q_TPC_L_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_L_newkaon_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_R_newkaon_RESCALED_Q_TPC_L_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_L_newproton_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_R_newproton_RESCALED_Q_TPC_L_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_L_newproton_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_R_newproton_RESCALED_Q_TPC_L_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_L_newproton_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_R_newproton_RESCALED_Q_TPC_L_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"}
  };

  // Indices for <uLQR> and <uRQL>: {{i_uLQR, i_iRQL}, ...}
  std::vector<std::pair<int, int>> LR_indices = {{0,1},{2,3},{4,5},
                                                {6,7},{8,9},{10,11},
                                                {12,13},{14,15},{16,17},
                                                {18,19},{20,21},{22,23},
                                                {24,25},{26,27},{28,29},
                                                {30,31},{32,33},{34,35},
                                                {36,37},{38,39},{40,41}
                                                };

  // Prepare names for vn
  std::vector<std::string> v2_tpc_names, v3_tpc_names, res2_tpc_names, res3_tpc_names;
  for (auto &name : corr_names_v2)
  {
    res2_tpc_names.push_back(Form("res_%s", name.second.c_str()));
    for (auto &cent : cent_ranges_v2)
    {
      v2_tpc_names.push_back(Form("vn_%s_cent%s%s", name.first.c_str(), std::to_string((int)cent.first).c_str(), std::to_string((int)cent.second).c_str()));
    }
  }
  for (auto &name : corr_names_v3)
  {
    res3_tpc_names.push_back(Form("res_%s", name.second.c_str()));
    for (auto &cent : cent_ranges_v3)
    {
      v3_tpc_names.push_back(Form("vn_%s_cent%s%s", name.first.c_str(), std::to_string((int)cent.first).c_str(), std::to_string((int)cent.second).c_str()));
    }
  }
  // Add merged <uLQR>+<uRQL> for vn
  for (auto &sides : LR_indices)
  {
    for (auto &cent : cent_ranges_v2)
    {
      v2_tpc_names.push_back(Form("vn_%s_%s_cent%s%s", corr_names_v2.at(sides.first).first.c_str(), corr_names_v2.at(sides.second).first.c_str(), std::to_string((int)cent
      .first).c_str(), std::to_string((int)cent
      .second).
      c_str()));
    }
    for (auto &cent : cent_ranges_v3)
    {
      v3_tpc_names.push_back(Form("vn_%s_%s_cent%s%s", corr_names_v3.at(sides.first).first.c_str(), corr_names_v3.at(sides.second).first.c_str(), std::to_string((int)cent
      .first).c_str(), std::to_string((int)cent
      .second).c_str()));
    }
  }
  
  TFile *foGraphs = new TFile(outFileName.Data(), "recreate");

  std::vector<std::pair< std::pair<Qn::DataContainerStatCalculate, std::string>, std::pair<Qn::DataContainerStatCalculate, std::string> >> corr2_tpc, corr3_tpc; // TPC-based: {{{uQ, name}, {QQ, name}}, ...}

  std::pair<std::pair<Qn::DataContainerStatCalculate, std::string>, std::pair<Qn::DataContainerStatCalculate, std::string>> dcsc_tmp;
  TGraphErrors *graph;

  if (isV2)
  {
    for (auto &name : corr_names_v2)
    {
      // Get <uQ> and <QQ> from input file
      if (GetDCStatCalculate(fiCorr, name.first, dcsc_tmp.first.first) && GetDCStatCalculate(fiCorr, name.second, dcsc_tmp.second.first))
      {
        std::cout << "Reading: " << name.first.c_str() << ", " << name.second.c_str() << std::endl;
        dcsc_tmp.first.second = name.first;
        dcsc_tmp.second.second = name.second;
        corr2_tpc.push_back(dcsc_tmp);
      }
    }
  }
  if (isV3)
  {
    for (auto &name : corr_names_v3)
    {
      // Get <uQ> and <QQ> from input file
      if (GetDCStatCalculate(fiCorr, name.first, dcsc_tmp.first.first) && GetDCStatCalculate(fiCorr, name.second, dcsc_tmp.second.first))
      {
        std::cout << "Reading: " << name.first.c_str() << ", " << name.second.c_str() << std::endl;
        dcsc_tmp.first.second = name.first;
        dcsc_tmp.second.second = name.second;
        corr3_tpc.push_back(dcsc_tmp);
      }
    }
  }

  // Calculate resolution and flow
  std::vector<Qn::DataContainerStatCalculate> v2_tpc, v3_tpc, res2_tpc, res3_tpc;
  if (isV2)
  {
    for (auto &corr : corr2_tpc)
    {
      auto res = Sqrt(corr.second.first);
      auto res_cent = res.Projection({"evCent"});
      res_cent.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
      res2_tpc.push_back(res_cent);

      auto v2_general = corr.first.first / res_cent;
      // Loop over centralities
      for (auto &cent : cent_ranges_v2)
      {
        auto v2_cent = v2_general.Rebin({"evCent", 1, cent.first, cent.second});
        auto v2_pT = v2_cent.Projection({"trPt"});
        if (corr.first.second.find(pid_names.at(1)) != std::string::npos)
          v2_pT = v2_pT.Rebin(axis_pt_v2_pion);
        else if (corr.first.second.find(pid_names.at(2)) != std::string::npos)
          v2_pT = v2_pT.Rebin(axis_pt_v2_kaon);
        else if (corr.first.second.find(pid_names.at(3)) != std::string::npos)
          v2_pT = v2_pT.Rebin(axis_pt_v2_proton);
        else
          v2_pT = v2_pT.Rebin(axis_pt_v2_ch);
        v2_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
        v2_tpc.push_back(v2_pT);
      }
    }
  }
  if (isV3)
  {
    for (auto &corr : corr3_tpc)
    {
      auto res = Sqrt(corr.second.first);
      auto res_cent = res.Projection({"evCent"});
      res_cent.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
      res3_tpc.push_back(res_cent);

      auto v3_general = corr.first.first / res_cent;
      // Loop over centralities
      for (auto &cent : cent_ranges_v3)
      {
        auto v3_cent = v3_general.Rebin({"evCent", 1, cent.first, cent.second});
        auto v3_pT = v3_cent.Projection({"trPt"});
        v3_pT = v3_pT.Rebin(axis_pt_v3_ch);
        if (corr.first.second.find(pid_names.at(1)) != std::string::npos)
          v3_pT = v3_pT.Rebin(axis_pt_v3_pion);
        else if (corr.first.second.find(pid_names.at(2)) != std::string::npos)
          v3_pT = v3_pT.Rebin(axis_pt_v3_kaon);
        else if (corr.first.second.find(pid_names.at(3)) != std::string::npos)
          v3_pT = v3_pT.Rebin(axis_pt_v3_proton);
        else
          v3_pT = v3_pT.Rebin(axis_pt_v3_ch);
        v3_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
        v3_tpc.push_back(v3_pT);
      }
    }
  }

  // Fill full vn (<uLQR> + <uRQL>)/2
  for (auto &sides : LR_indices)
  {
    if (isV2)
    {
      for (auto &cent : cent_ranges_v2)
      {
        auto v2_general = Merge(corr2_tpc.at(sides.first).first.first, corr2_tpc.at(sides.second).first.first) / res2_tpc.at(sides.first);
        auto v2_cent = v2_general.Rebin({"evCent", 1, cent.first, cent.second});
        auto v2_pT = v2_cent.Projection({"trPt"});
        if (corr2_tpc.at(sides.first).first.second.find(pid_names.at(1)) != std::string::npos)
          v2_pT = v2_pT.Rebin(axis_pt_v2_pion);
        else if (corr2_tpc.at(sides.first).first.second.find(pid_names.at(2)) != std::string::npos)
          v2_pT = v2_pT.Rebin(axis_pt_v2_kaon);
        else if (corr2_tpc.at(sides.first).first.second.find(pid_names.at(3)) != std::string::npos)
          v2_pT = v2_pT.Rebin(axis_pt_v2_proton);
        else
          v2_pT = v2_pT.Rebin(axis_pt_v2_ch);
        v2_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
        v2_tpc.push_back(v2_pT);
      }
    }
    if (isV3)
    {
      for (auto &cent : cent_ranges_v3)
      {
        auto v3_general = Merge(corr3_tpc.at(sides.first).first.first, corr3_tpc.at(sides.second).first.first) / res3_tpc.at(sides.first);
        auto v3_cent = v3_general.Rebin({"evCent", 1, cent.first, cent.second});
        auto v3_pT = v3_cent.Projection({"trPt"});
        if (corr3_tpc.at(sides.first).first.second.find(pid_names.at(1)) != std::string::npos)
          v3_pT = v3_pT.Rebin(axis_pt_v3_pion);
        else if (corr3_tpc.at(sides.first).first.second.find(pid_names.at(2)) != std::string::npos)
          v3_pT = v3_pT.Rebin(axis_pt_v3_kaon);
        else if (corr3_tpc.at(sides.first).first.second.find(pid_names.at(3)) != std::string::npos)
          v3_pT = v3_pT.Rebin(axis_pt_v3_proton);
        else
          v3_pT = v3_pT.Rebin(axis_pt_v3_ch);
        v3_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
        v3_tpc.push_back(v3_pT);
      }
    }
  }

  std::vector<std::string> v_written_res;
  std::string str_tmp;

  // Write resolution
  std::cout << "Written output resolutions:" << std::endl;
  if (isV2)
  {
    v_written_res.clear();
    foGraphs->mkdir("res2_tpc");
    foGraphs->cd("res2_tpc");
    for (int i = 0; i < res2_tpc.size(); i++)
    {
      str_tmp = res2_tpc_names.at(i);
      if (std::any_of(v_written_res.begin(), v_written_res.end(), [&str_tmp](std::string str){ return str == str_tmp; }) )
        continue;
      graph = Qn::DataContainerHelper::ToTGraph(res2_tpc.at(i));
      graph->SetName(Form("res_%s", res2_tpc_names.at(i).c_str()));
      graph->SetTitle(Form("res_%s;%s", res2_tpc_names.at(i).c_str(), "Centrality, %;R"));
      graph->Write();
      v_written_res.push_back(res2_tpc_names.at(i));
      std::cout << "\t" << graph->GetName() << std::endl;
    }
  }
  if (isV3)
  {
    foGraphs->mkdir("res3_tpc");
    foGraphs->cd("res3_tpc");
    v_written_res.clear();
    for (int i = 0; i < res3_tpc.size(); i++)
    {
      str_tmp = res3_tpc_names.at(i);
      if (std::any_of(v_written_res.begin(), v_written_res.end(), [&str_tmp](std::string str){ return str == str_tmp; }) )
        continue;
      graph = Qn::DataContainerHelper::ToTGraph(res3_tpc.at(i));
      graph->SetName(Form("res_%s", res3_tpc_names.at(i).c_str()));
      graph->SetTitle(Form("res_%s;%s", res3_tpc_names.at(i).c_str(), "Centrality, %;R"));
      graph->Write();
      v_written_res.push_back(res3_tpc_names.at(i));
      std::cout << "\t" << graph->GetName() << std::endl;
    }
  }

  // Write flow
  std::cout << "Writing output flow:" << std::endl;
  if (isV2)
  {
    foGraphs->mkdir("v2_tpc");
    foGraphs->cd("v2_tpc");
    for (int i = 0; i < v2_tpc.size(); i++)
    {
      graph = Qn::DataContainerHelper::ToTGraph(v2_tpc.at(i));
      graph->SetName(Form("%s", v2_tpc_names.at(i).c_str()));
      graph->SetTitle(Form("%s;%s", v2_tpc_names.at(i).c_str(), "p_{T}, GeV/c;v_{2}"));
      graph->Write();
      std::cout << "\t" << graph->GetName() << std::endl;
    }
  }
  if (isV3)
  {
    foGraphs->mkdir("v3_tpc");
    foGraphs->cd("v3_tpc");
    for (int i = 0; i < v3_tpc.size(); i++)
    {
      graph = Qn::DataContainerHelper::ToTGraph(v3_tpc.at(i));
      graph->SetName(Form("%s", v3_tpc_names.at(i).c_str()));
      graph->SetTitle(Form("%s;%s", v3_tpc_names.at(i).c_str(), "p_{T}, GeV/c;v_{3}"));
      graph->Write();
      std::cout << "\t" << graph->GetName() << std::endl;
    }
  }

  foGraphs->Close();
}
