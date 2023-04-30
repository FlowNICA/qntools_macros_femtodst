#include "Functions.C"

void Draw_SPEP(TString inFileName="", TString outFileName="./test_graphs.root")
{
    TFile *fiCorr = new TFile(inFileName.Data(),"read");
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
    const std::vector<std::pair<float, float>> cent_ranges = {{0., 5.}, {5., 10.}, {10., 20.}, {20., 30.}, {30., 40.}, {40., 60.}, {10., 40.}};

    // Set PID namings
    const int npid = 1; // Total 6 types: h+-, proton, K+, K-, pi+, pi-
    const std::vector<std::string> pidnames = {"ch"}; // Total: {"ch", "protons", "kaons", "akaons", "pions", "apions"};

    // Set correction and methods naming
    const std::vector<std::string> corrnames = {"RESCALED"}; //Total: {"PLAIN", "RECENTERED", "TWIST", "RESCALED"};
    const std::vector<std::string> QQ_methods = {"XX", "YY", "SP", "coscos", "sinsin", "EP"};
    const std::vector<std::string> uQ_methods = {"XX", "YY", "SP", "Xcos", "Ysin", "EP"};

    // Set up u-, Q-vectors naming
    const std::vector<std::string> uv_tpc_names = {"u_TPC_L", "u_TPC_R"};
    const std::vector<std::string> Qv_tpc_names = {"Q_TPC_L", "Q_TPC_R"};

    // Set up event-wise axis name
    const std::string axisName="evCent";
    
    // Construct <QQ> correlation names
    std::vector<std::string> QQ_tpc_names_test, uQ_tpc_names_test; 
    std::vector<std::string> QQ_tpc_names, uQ_tpc_names; 
    for (int iQ1=0; iQ1 < Qv_tpc_names.size(); iQ1++)
    {
      for (int iQ2=0; iQ2 < Qv_tpc_names.size(); iQ2++)
      {
        if (iQ1 == iQ2) continue;
        for (int iCorr=0; iCorr < corrnames.size(); iCorr++)
        {
          for (int iMethod=0; iMethod < QQ_methods.size(); iMethod++)
          {
            for (int ipid=0; ipid < pidnames.size(); ipid++)
            {
              QQ_tpc_names_test.push_back({Qv_tpc_names.at(iQ1)+"_"+pidnames.at(ipid)+"_"+corrnames.at(iCorr)+"_"+Qv_tpc_names.at(iQ2)+"_"+pidnames.at(0)+"_"+corrnames.at(iCorr)+"_"+QQ_methods.at(iMethod)+"_"+axisName});
            }
          }
        }
      }
    }
    // Construct <uQ> correlation names
    for (int iu=0; iu < uv_tpc_names.size(); iu++)
    {
      for (int iQ=0; iQ < Qv_tpc_names.size(); iQ++)
      {
        if (iu == iQ) continue;
        for (int iCorr=0; iCorr < corrnames.size(); iCorr++)
        {
          for (int iMethod=0; iMethod < uQ_methods.size(); iMethod++)
          {
            for (int ipid=0; ipid < pidnames.size(); ipid++)
            {
              uQ_tpc_names_test.push_back({uv_tpc_names.at(iu)+"_"+pidnames.at(ipid)+"_"+corrnames.at(iCorr)+"_"+Qv_tpc_names.at(iQ)+"_"+pidnames.at(0)+"_"+corrnames.at(iCorr)+"_"+uQ_methods.at(iMethod)+"_"+axisName});
            }
          }
        }
      }
    }

    TFile *foGraphs = new TFile(outFileName.Data(),"recreate");

    // <QQ> - {coscos, sinsin, EP, XX, YY, SP}
    std::vector<Qn::DataContainerStatCalculate> QQ_tpc; // TPC-based
    // <uQ> - {Xcos, Ysin, EP, XX, YY, SP}
    std::vector<Qn::DataContainerStatCalculate> uQ_tpc; // TPC-based

    Qn::DataContainerStatCalculate dcsc_tmp;
    TGraphErrors *graph;

    // Getting correlations from the TFile
    // <QQ>
    for (auto &name : QQ_tpc_names_test)
    {
      if (GetDCStatCalculate(fiCorr, name, dcsc_tmp))
      {
        QQ_tpc_names.push_back(name);
        QQ_tpc.push_back(dcsc_tmp);
      }
    }
    // <uQ>
    for (auto &name : uQ_tpc_names_test)
    {
      if (GetDCStatCalculate(fiCorr, name, dcsc_tmp))
      {
        uQ_tpc_names.push_back(name);
        uQ_tpc.push_back(dcsc_tmp);
      }
    }

    std::cout << "QQ correlations found:" << std::endl;
    for (auto &element : QQ_tpc_names)
    {
      std::cout << "\t" << element << std::endl;
    }
    std::cout << "uQ correlations found:" << std::endl;
    for (auto &element : uQ_tpc_names)
    {
      std::cout << "\t" << element << std::endl;
    }

    // Calculate resolution: R = sqrt(<QQ>)
    std::vector<Qn::DataContainerStatCalculate> res2_tpc;
    for (auto &qq : QQ_tpc)
    {
      auto res = Sqrt( qq );
      auto res_cent = res.Projection({"evCent"});
      res_cent.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
      res2_tpc.push_back(res_cent);
    }

    // Calculate flow: vn = <uQ>/R
    std::vector<Qn::DataContainerStatCalculate> v2_tpc;
    std::vector<std::string> v2_tpc_names;
    int ires=0;
    for (int i = 0; i < (int)((float)uQ_tpc.size()/2.); i++)
    {
      auto v2_obs = Merge(uQ_tpc.at(i), uQ_tpc.at((int)((float)uQ_tpc.size()/2.)+i)); // Merge <uLQR> and <uRQL>
      Qn::DataContainerStatCalculate v2_general;
      for (int j=0; j < uQ_methods.size(); j++)
      {
        if (uQ_tpc_names.at(i).find(uQ_methods.at(j)) == std::string::npos) continue; // methods for <uQ> and <QQ> are not equal - skip
        if (uQ_tpc_names.at(i).find("_SP_") != std::string::npos || 
          uQ_tpc_names.at(i).find("_EP_") != std::string::npos)
        {
          v2_general = v2_obs / res2_tpc.at(j); // full EP or SP methods
        }
        else
        {
          v2_general = sqrt(2.) * v2_obs / res2_tpc.at(j); // components: XX, YY, ...
        }
        ires = j;
        break;
      }
      // Loop over centralities
      for (auto &cent : cent_ranges)
      {
        auto v2_cent = v2_general.Rebin({"evCent", 1, cent.first, cent.second});
        auto v2_pT = v2_cent.Projection({"trPt"});
        v2_pT = v2_pT.Rebin({"trPt", {0.2, 0.4, 0.6, 0.8, 1., 1.5, 2., 3.}});
        v2_pT.SetErrors(Qn::Stat::ErrorType::BOOTSTRAP);
        v2_tpc_names.push_back("vn_TPC_ch_"+uQ_methods.at(ires)+"_"+"cent"+to_string((int)cent.first)+to_string((int)cent.second));
        v2_tpc.push_back(v2_pT);
      }
    }

    // Write resolution
    foGraphs->mkdir("res2_tpc");
    foGraphs->cd("res2_tpc");
    for (int i=0; i < res2_tpc.size(); i++)
    {
      graph = Qn::DataContainerHelper::ToTGraph(res2_tpc.at(i));
      graph->SetName(Form("res_%s", QQ_tpc_names.at(i).c_str()));
      graph->SetTitle(Form("res_%s;%s", QQ_tpc_names.at(i).c_str(), "Centrality, %;R"));
      graph->Write();
    }

    // Write flow
    foGraphs->mkdir("v2_tpc");
    foGraphs->cd("v2_tpc");
    for (int i=0; i < v2_tpc.size(); i++)
    {
      graph = Qn::DataContainerHelper::ToTGraph(v2_tpc.at(i));
      graph->SetName(Form("%s", v2_tpc_names.at(i).c_str()));
      graph->SetTitle(Form("%s;%s", v2_tpc_names.at(i).c_str(), "p_{T}, GeV/c;v_{2}"));
      graph->Write();
    }

    foGraphs->Close();
}