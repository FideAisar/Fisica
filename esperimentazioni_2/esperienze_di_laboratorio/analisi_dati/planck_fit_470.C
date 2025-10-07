{
    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Planck", 400, 400);
    c1->SetGrid();
    

    // DATA /////////////////////////////////////////////////////////////////////////////
    const int n = 28;
	double x_v[n] = { 0.70,99.6,199.2,299.7,398.7,500.9,553.7,601.2,611.4,
		620.3,630.3,640.8,651.3,701.7,751.7,799.5,849.8,900.9,999.7,1102.0,
		1201.0,1299.0,1403.0,1499.0,1753.0,2004.0,2504.0,2990.0 };
        
    double x_v_err[n] = { 1,1,2,3,4,6,6,7,7,7,7,7,7,8,8,8,9,10,10,12,13,13,15,15,18,21,26,30 };

	double y_i[n] = { 0.0298,0.0209,0.0136,0.0080,0.0041,0.0015,0.0007,0.0001,0.0000,-0.0001,-0.0002,-0.0003,-0.0004,-0.0007,-0.0009,-0.0010,-0.0011,-0.0012,-0.0013,-0.0014,-0.0014,-0.0015,-0.0014,-0.0015,-0.0015,-0.0015,-0.0017,-0.0017 };

	double y_i_err[n] = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,
		0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,
		0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001 };


    // GRAPH /////////////////////////////////////////////////////////////////////////////
    TGraphErrors *gr = new TGraphErrors(n, x_v, y_i, x_v_err, y_i_err); // x_freq_err=0
    gr->SetTitle("I(V) (470)");
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
	leg->AddEntry((TObject*)0, TString::Format("Range fit: [%.1f, %.qf]", min, max), "");
	leg->AddEntry((TObject*)0, TString::Format("m = %.3f #pm %.3f", m, m_err), "");
	leg->AddEntry((TObject*)0, TString::Format("m/#sigma_{m} = %.2f", fabs(m/m_err)), "");
	leg->AddEntry((TObject*)0, TString::Format("#alpha = %.1f (soglia)", alpha), "");
	
	leg->Draw();
}

