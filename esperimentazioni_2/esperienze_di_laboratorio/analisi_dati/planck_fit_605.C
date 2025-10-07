{
    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Planck", 400, 400);
    c1->SetGrid();
    

    // DATA /////////////////////////////////////////////////////////////////////////////
    const int n = 27;
	double x_v[n] = { 0.90,19.8,39.9,61.6,80.7,99.9,120.3,140.8,151.4,160.6,169.8,179.8,199.5,236.2,278.6,320.4,356.6,401.2,499.0,599.8,700.4,858.3,949.8,1102.0,1251.0,2007.0,2511.0 };
        
    double x_v_err[n] = { 1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,4,4,5,5,6,8,9,14,16,18,25,30};

	double y_i[n] = { 0.0028,0.0021,0.0016,0.0011,0.0007,0.0005,0.0003,0.0001,0.0000,-0.0001,-0.0001,-0.0002,-0.0003,-0.0004,-0.0005,-0.0005,-0.0006,-0.0006,-0.0007,-0.0008,-0.0009,-0.0009,-0.0010,-0.0010,-0.0010,-0.0013,-0.0013};

	double y_i_err[n] = { 0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002};


     // GRAPH /////////////////////////////////////////////////////////////////////////////
    TGraphErrors *gr = new TGraphErrors(n, x_v, y_i, x_v_err, y_i_err); // x_freq_err=0
    gr->SetTitle("I(V) (605)");
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

	while (fabs(m/m_err) > alpha && min < x_v[n-2]) {
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
	leg->AddEntry((TObject*)0, TString::Format("Range fit: [%.1f, %.f]", min, max), "");
	leg->AddEntry((TObject*)0, TString::Format("m = %.3f #pm %.3f", m, m_err), "");
	leg->AddEntry((TObject*)0, TString::Format("m/#sigma_{m} = %.2f", fabs(m/m_err)), "");
	leg->AddEntry((TObject*)0, TString::Format("#alpha = %.1f (soglia)", alpha), "");
	
	leg->Draw();
}

