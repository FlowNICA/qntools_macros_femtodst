#include "Functions.C"

void Draw_SPEP(TString inFileName = "", TString outFileName = "./test_graphs.root")
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
  const std::vector<std::pair<float, float>> cent_ranges_v2 = {{0., 5.}, {5., 10.}, {10., 20.}, {20., 30.}, {30., 40.}, {40., 60.}, {10., 40.}};
  const std::vector<std::pair<float, float>> cent_ranges_v3 = {{0., 10.}, {10., 20.}, {20., 30.}, {30., 40.}, {40., 50.}, {50., 60.}, {10., 40.}};

  Qn::AxisD axis_pt_v2("trPt", {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.8, 3.25, 3.75, 4.5});
  Qn::AxisD axis_pt_v3("trPt", {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.8, 3.25, 3.75, 4.5});

  // Names for <uQ> and <QQ> correlations for vn: {{uQ, QQ}, ...}
  std::vector<std::pair<std::string, std::string>> corr_names_v2 = {
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_SP2_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP2_evCent"},
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_X2X2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X2X2_evCent"},
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_Y2Y2_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y2Y2_evCent"}
  };
  std::vector<std::pair<std::string, std::string>> corr_names_v3 = {
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_SP3_evCent",  "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_SP3_evCent"},
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_X3X3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_X3X3_evCent"},
    {"u_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"},
    {"u_TPC_R_ch_RESCALED_Q_TPC_L_ch_RESCALED_Y3Y3_evCent", "Q_TPC_L_ch_RESCALED_Q_TPC_R_ch_RESCALED_Y3Y3_evCent"}
  };

  // Indices for <uLQR> and <uRQL>: {{i_uLQR, i_iRQL}, ...}
  std::vector<std::pair<int, int>> LR_indices = {{0,1},{2,3},{4,5}};

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

  // <QQ> - {coscos, sinsin, EP, XX, YY, SP}
  std::vector<std::pair<Qn::DataContainerStatCalculate, Qn::DataContainerStatCalculate>> corr2_tpc, corr3_tpc; // TPC-based: {{uQ, QQ}, ...}

  std::pair<Qn::DataContainerStatCalculate, Qn::DataContainerStatCalculate> dcsc_tmp;
  TGraphErrors *graph;

  for (auto &name : corr_names_v2)
  {
    // Get <uQ> and <QQ> from input file
    if (GetDCStatCalculate(fiCorr, name.first, dcsc_tmp.first) && GetDCStatCalculate(fiCorr, name.second, dcsc_tmp.second))
    {
      std::cout << "Reading: " << name.first.c_str() << ", " << name.second.c_str() << std::endl;
      corr2_tpc.push_back(dcsc_tmp);
    }
  }
  for (auto &name : corr_names_v3)
  {
    // Get <uQ> and <QQ> from input file
    if (GetDCStatCalculate(fiCorr, name.first, dcsc_tmp.first) && GetDCStatCalculate(fiCorr, name.second, dcsc_tmp.second))
    {
      std::cout << "Reading: " << name.first.c_str() << ", " << name.second.c_str() << std::endl;
      corr3_tpc.push_back(dcsc_tmp);
    }
  }

  // Calculate resolution and flow
  std::vector<Qn::DataContainerStatCalculate> v2_tpc, v3_tpc, res2_tpc, res3_tpc;
  for (auto &corr : corr2_tpc)
  {
    auto res = Sqrt(corr.second);
    auto res_cent = res.Projection({"evCent"});
    res_cent.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
    res2_tpc.push_back(res_cent);

    auto v2_general = corr.first / res_cent;
    // Loop over centralities
    for (auto &cent : cent_ranges_v2)
    {
      auto v2_cent = v2_general.Rebin({"evCent", 1, cent.first, cent.second});
      auto v2_pT = v2_cent.Projection({"trPt"});
      v2_pT = v2_pT.Rebin(axis_pt_v2);
      v2_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
      v2_tpc.push_back(v2_pT);
    }
  }
  for (auto &corr : corr3_tpc)
  {
    auto res = Sqrt(corr.second);
    auto res_cent = res.Projection({"evCent"});
    res_cent.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
    res3_tpc.push_back(res_cent);

    auto v3_general = corr.first / res_cent;
    // Loop over centralities
    for (auto &cent : cent_ranges_v3)
    {
      auto v3_cent = v3_general.Rebin({"evCent", 1, cent.first, cent.second});
      auto v3_pT = v3_cent.Projection({"trPt"});
      v3_pT = v3_pT.Rebin(axis_pt_v3);
      v3_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
      v3_tpc.push_back(v3_pT);
    }
  }

  // Fill full vn (<uLQR> + <uRQL>)/2
  for (auto &sides : LR_indices)
  {
    for (auto &cent : cent_ranges_v2)
    {
      auto v2_general = Merge(corr2_tpc.at(sides.first).first, corr2_tpc.at(sides.second).first) / res2_tpc.at(sides.first);
      auto v2_cent = v2_general.Rebin({"evCent", 1, cent.first, cent.second});
      auto v2_pT = v2_cent.Projection({"trPt"});
      v2_pT = v2_pT.Rebin(axis_pt_v2);
      v2_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
      v2_tpc.push_back(v2_pT);
    }
    for (auto &cent : cent_ranges_v3)
    {
      auto v3_general = Merge(corr3_tpc.at(sides.first).first, corr3_tpc.at(sides.second).first) / res3_tpc.at(sides.first);
      auto v3_cent = v3_general.Rebin({"evCent", 1, cent.first, cent.second});
      auto v3_pT = v3_cent.Projection({"trPt"});
      v3_pT = v3_pT.Rebin(axis_pt_v3);
      v3_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
      v3_tpc.push_back(v3_pT);
    }
  }

  // Write resolution
  foGraphs->mkdir("res2_tpc");
  foGraphs->cd("res2_tpc");
  for (int i = 0; i < res2_tpc.size(); i++)
  {
    graph = Qn::DataContainerHelper::ToTGraph(res2_tpc.at(i));
    graph->SetName(Form("res_%s", res2_tpc_names.at(i).c_str()));
    graph->SetTitle(Form("res_%s;%s", res2_tpc_names.at(i).c_str(), "Centrality, %;R"));
    graph->Write();
  }
  foGraphs->mkdir("res3_tpc");
  foGraphs->cd("res3_tpc");
  for (int i = 0; i < res3_tpc.size(); i++)
  {
    graph = Qn::DataContainerHelper::ToTGraph(res3_tpc.at(i));
    graph->SetName(Form("res_%s", res3_tpc_names.at(i).c_str()));
    graph->SetTitle(Form("res_%s;%s", res3_tpc_names.at(i).c_str(), "Centrality, %;R"));
    graph->Write();
  }

  // Write flow
  foGraphs->mkdir("v2_tpc");
  foGraphs->cd("v2_tpc");
  for (int i = 0; i < v2_tpc.size(); i++)
  {
    graph = Qn::DataContainerHelper::ToTGraph(v2_tpc.at(i));
    graph->SetName(Form("%s", v2_tpc_names.at(i).c_str()));
    graph->SetTitle(Form("%s;%s", v2_tpc_names.at(i).c_str(), "p_{T}, GeV/c;v_{2}"));
    graph->Write();
  }
  foGraphs->mkdir("v3_tpc");
  foGraphs->cd("v3_tpc");
  for (int i = 0; i < v3_tpc.size(); i++)
  {
    graph = Qn::DataContainerHelper::ToTGraph(v3_tpc.at(i));
    graph->SetName(Form("%s", v3_tpc_names.at(i).c_str()));
    graph->SetTitle(Form("%s;%s", v3_tpc_names.at(i).c_str(), "p_{T}, GeV/c;v_{2}"));
    graph->Write();
  }

  std::cout << "Written output resolutions:" << std::endl;
  for (int i = 0; i < res2_tpc.size(); i++)
    std::cout << "\t" << Form("res_%s", res2_tpc_names.at(i).c_str()) << std::endl;
  for (int i = 0; i < res3_tpc.size(); i++)
    std::cout << "\t" << Form("res_%s", res3_tpc_names.at(i).c_str()) << std::endl;

  std::cout << "Writing output flow:" << std::endl;
  for (int i = 0; i < v2_tpc.size(); i++)
    std::cout << "\t" << Form("%s", v2_tpc_names.at(i).c_str()) << std::endl;
  for (int i = 0; i < v3_tpc.size(); i++)
    std::cout << "\t" << Form("%s", v3_tpc_names.at(i).c_str()) << std::endl;

  foGraphs->Close();
}