#include "Functions.C"

void Draw_QCumulants(TString iFileName, TString oFileName)
{
	TFile *fiCorr = new TFile(iFileName.Data(), "read");
	TFile *fo = new TFile(oFileName.Data(), "recreate");
	if (iFileName == "" || !fiCorr)
	{
		std::cerr << "No input file was provided!" << std::endl;
		return;
	}

	const int nsides = 3;		// 0 = L and 1 = R and 2 = F
	const int nmethods = 2; // 2; // 0 = v2{2} and 1 = v2{4}
	Qn::DataContainerStatCalculate qcorr_TPC[nsides][nmethods];
	Qn::DataContainerStatCalculate qc2_TPC[nsides][nmethods];
	Qn::DataContainerStatCalculate qv2_TPC[nsides][nmethods];
	TGraphErrors *graph;
	std::string name;
	std::string qcorr_names[nsides][nmethods] = {{"Q_TPC_L_ch_PLAIN_c2_evCentevVtxZ",
																								"Q_TPC_L_ch_PLAIN_c4_evCentevVtxZ"},
																							 {"Q_TPC_R_ch_PLAIN_c2_evCentevVtxZ",
																								"Q_TPC_R_ch_PLAIN_c4_evCentevVtxZ"},
																							 {"Q_TPC_F_ch_PLAIN_c2_evCentevVtxZ",
																								"Q_TPC_F_ch_PLAIN_c4_evCentevVtxZ"}};

	bool isCorrThere = false;

	for (int iside = 0; iside < nsides; iside++)
	{
		for (int imethod = 0; imethod < nmethods; imethod++)
		{
			name = qcorr_names[iside][imethod];
			isCorrThere = GetDCStatCalculate(fiCorr, name, qcorr_TPC[iside][imethod]);
			if (!isCorrThere)
				continue;
		}
	}

	for (int iside = 0; iside < nsides; iside++)
	{
		for (int imethod = 0; imethod < nmethods; imethod++)
		{
			qcorr_TPC[iside][imethod] = qcorr_TPC[iside][imethod].Rebin({"evCent", 14, 0., 70.});
			qcorr_TPC[iside][imethod] = qcorr_TPC[iside][imethod].Projection({"evCent"});
		}
	}

	for (int iside = 0; iside < nsides; iside++)
	{
		qc2_TPC[iside][0] = qcorr_TPC[iside][0];
		qv2_TPC[iside][0] = Sqrt(qc2_TPC[iside][0]);

		qc2_TPC[iside][1] = qcorr_TPC[iside][1] - 2. * qcorr_TPC[iside][0] * qcorr_TPC[iside][0];
		qv2_TPC[iside][1] = -1. * qc2_TPC[iside][1];
		// qv2_TPC[iside][1] = qc2_TPC[iside][1] * qc2_TPC[iside][1];
		// qv2_TPC[iside][1] = Sqrt(qv2_TPC[iside][1]);
		qv2_TPC[iside][1] = Sqrt(Sqrt(qv2_TPC[iside][1]));
	}

	// for (int iside = 0; iside < nsides; iside++)
	// {
	// 	for (int imethod = 0; imethod < nmethods; imethod++)
	// 	{
	// 		qc2_TPC[iside][imethod] = qc2_TPC[iside][imethod].Rebin({"evCent", 14, 0., 70.});
	// 		qc2_TPC[iside][imethod] = qc2_TPC[iside][imethod].Projection({"evCent"});
	// 		qv2_TPC[iside][imethod] = qv2_TPC[iside][imethod].Rebin({"evCent", 14, 0., 70.});
	// 		qv2_TPC[iside][imethod] = qv2_TPC[iside][imethod].Projection({"evCent"});
	// 	}
	// }

	fo->cd();

	for (int iside = 0; iside < nsides; iside++)
	{
		for (int imethod = 0; imethod < nmethods; imethod++)
		{
			graph = Qn::DataContainerHelper::ToTGraph(qc2_TPC[iside][imethod]);
			name = "gr_c2n_" + qcorr_names[iside][imethod];
			graph->SetName(name.c_str());
			if (imethod == 0)
				name += ";centrality, %;c_{2}{2}";
			if (imethod == 1)
				name += ";centrality, %;c_{2}{4}";
			graph->SetTitle(name.c_str());
			graph->Write();

			graph = Qn::DataContainerHelper::ToTGraph(qv2_TPC[iside][imethod]);
			name = "gr_v2n_" + qcorr_names[iside][imethod];
			graph->SetName(name.c_str());
			if (imethod == 0)
				name += ";centrality, %;v_{2}{2}";
			if (imethod == 1)
				name += ";centrality, %;v_{2}{4}";
			graph->SetTitle(name.c_str());
			graph->Write();
		}
	}
}
