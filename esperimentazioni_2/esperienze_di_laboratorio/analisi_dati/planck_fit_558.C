{
    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Planck", 400, 400);
    c1->SetGrid();
    

    // DATA /////////////////////////////////////////////////////////////////////////////
    const int n = 35;
	double x_v[n] = { 0.80,29.7,60.8,91.1,121.5,150.4,179.4,211.3,239.6,272.0,301.0,331.5,361.9,390.2,421.4,451.9,480.4,510.1,554.0,601.3,699.1,802.6,903.4,1003.0,1104.0,1202.0,1300.0,1406.0,1503.0,1602.0,1702.0,2006.0,2399.0,2800.0,2990.0};
        
    double x_v_err[n] = {  1,1,1,1,2,2,2,3,3,3,4,4,4,4,5,5,5,6,6,7,7,9,10,15,16,17,18,19,20,21,22,25,29,33,35};

	double y_i[n] = {  0.00287,0.00223,0.00167,0.00119,0.00078,0.00047,0.00021,-0.00001,-0.00013,-0.00024,-0.00032,-0.00035,-0.00039,-0.00040,-0.00041,-0.00043,-0.00044,-0.00046,-0.00049,-0.00051,-0.00054,-0.00057,-0.00059,-0.00061,-0.00061,-0.00062,-0.00063,-0.00064,-0.00066,-0.00065,-0.00067,-0.00070,-0.00075,-0.00075,-0.00079};

	double y_i_err[n] = {0.00006,0.00002,0.00002,0.00008,0.00004,0.00020,0.00060,0.00002,0.00002,0.00002,0.00006,0.00006,0.00006,0.00004,0.00012,0.00004,0.00008,0.00004,0.00004,0.00004,0.00002,0.00002,0.00004,0.00004,0.00008,0.00006,0.00004,0.00004,0.00004,0.00004,0.00010,0.00004,0.00004,0.00008,0.00008};


     // GRAPH /////////////////////////////////////////////////////////////////////////////
    TGraphErrors *gr = new TGraphErrors(n, x_v, y_i, x_v_err, y_i_err); // x_freq_err=0
    gr->SetTitle("I(V) (558)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(1);
	gr->SetFillColor(kBlue-10);
    gr->SetFillStyle(3001); // semi-transparent fill
    gr->Draw("ALP");

 
    // FIT  /////////////////////////////////////////////////////////////////////////////
	double min = 1000;
	double max = x_v[n-1];
	const double alpha = 1.5;  

	TF1 *ff1 = new TF1("ff1", "pol1");
	ff1->SetParNames("m", "q");
	ff1->SetLineColor(kRed);
	ff1->SetLineWidth(2);

	gr->Fit("ff1", "SR+", "Q", min, max);  // R=range, +=improvements
	double m = ff1->GetParameter(0);
	double m_err = ff1->GetParError(0);

	while (fabs(m/m_err) > alpha && min < x_v[n-3]) {
    	gr->Fit("ff1", "SR+", "Q", min, max);
    	m = ff1->GetParameter(0);
    	m_err = ff1->GetParError(0);
    
    	cout << "Fit range: [" << min << ", " << max << "]" << endl;
    	cout << "m/m_err = " << (m/m_err) << endl;
    
    	min += 1;  // Move the lower bound up
	}

	if (fabs(m/m_err) < alpha) {
		cout << "Significance dropped below threshold (|m/m_err| < " << alpha << ")" << endl;
	} else {
    	cout << "Reached end of data range" << endl;
	}

	
    // LEGENDA  /////////////////////////////////////////////////////////////////////////
	TLegend *leg = new TLegend(0.60, 0.60, 0.90, 0.90);
	leg->SetBorderSize(1);
	leg->SetTextSize(0.03);
	leg->SetHeader("Risultati del Fit", "C");  // Titolo centrato
	
	leg->AddEntry(ff1, "fit orizzontale", "l");  
	leg->AddEntry((TObject*)0, TString::Format("Range fit: [%.1f, %.qf]", min, max), "");
	leg->AddEntry((TObject*)0, TString::Format("m = %.3f #pm %.3f", m, m_err), "");
	leg->AddEntry((TObject*)0, TString::Format("m/#sigma_{m} = %.2f", fabs(m/m_err)), "");
	leg->AddEntry((TObject*)0, TString::Format("#alpha = %.1f (soglia)", alpha), "");
	
	leg->Draw();
}

