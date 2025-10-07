{
    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Planck", 400, 400);
    c1->SetGrid();
    

    // DATA /////////////////////////////////////////////////////////////////////////////
    const int n = 35;
	double x_v[n] = { 0.50,102.5,201.9,300.9,400.1,500.0,600.5,700.4,721.0,740.5,759.9,
		781.1,801.0,819.7,841.0,859.9,910.4,960.6,1009.0,1062.0,1108.0,1162.0,1209.0,1261.0,
		1311.0,1358.0,1409.0,1460.0,1509.0,1560.0,1611.0,1708.0,1861.0,2011.0,2505.0 };
        
    double x_v_err[n] = { 0.5, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10,
		15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 24, 25, 30 };

	double y_i[n] = {0.0319,0.0249,0.0185,0.0131,0.0087,0.0053,0.0028,0.0012,0.0010,0.0007,0.0006,0.0004,0.0002,0.0001,0.0000,-0.0001,-0.0004,-0.0006,-0.0008,-0.0009,-0.0009,-0.0008,-0.0010,-0.0011,-0.0011,-0.0011,-0.0012,-0.0012,-0.0012,-0.0012,-0.0012,-0.0012,-0.0013,-0.0013,-0.0014 };


	double y_i_err[n] = { 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
		0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
		0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
		0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
		0.0001 };


     // GRAPH /////////////////////////////////////////////////////////////////////////////
    TGraphErrors *gr = new TGraphErrors(n, x_v, y_i, x_v_err, y_i_err); // x_freq_err=0
    gr->SetTitle("I(V) (430)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(1);
	gr->SetFillColor(kBlue-10);
    gr->SetFillStyle(3001); // semi-transparent fill
    gr->Draw("AP3");

 
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
	leg->AddEntry((TObject*)0, TString::Format("m = %.3f #pm %.5f", m, m_err), "");
	leg->AddEntry((TObject*)0, TString::Format("m/#sigma_{m} = %.2f", fabs(m/m_err)), "");
	leg->AddEntry((TObject*)0, TString::Format("#alpha = %.1f (soglia)", alpha), "");
	
	leg->Draw();
}

