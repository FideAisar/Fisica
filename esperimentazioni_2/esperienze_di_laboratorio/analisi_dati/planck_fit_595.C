{
    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Planck", 400, 400);
    c1->SetGrid();
    

    // DATA /////////////////////////////////////////////////////////////////////////////
    const int n = 29;
	double x_v[n] = { 0.01,25.2,50.4,76.0,101.5,124.9,148.9,175.0,201.2,211.2,223.8,252.2,302.4,350.7,401.6,451.2,502.2,551.0,605.7,702.1,801.4,904.3,1003.0,1235.0,1499.0,1752.0,2001.0,2501.0,2990.0 };
        
    double x_v_err[n] = {  1,1,1,1,2,2,2,2,3,3,3,3,4,4,5,5,6,6,7,8,9,10,15,17,20,23,25,30,35};

	double y_i[n] = {  0.0296,0.0252,0.0188,0.0134,0.0090,0.0059,0.0032,0.0010,0.0003,-0.0009,-0.0014,-0.0023,-0.0033,-0.0040,-0.0045,-0.0050,-0.0055,-0.0059,-0.0062,-0.0067,-0.0073,-0.0076,-0.0079,-0.0084,-0.0087,-0.0089,-0.0090,-0.0092,-0.0094};

	double y_i_err[n] = {0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002,0.0002};


     // GRAPH /////////////////////////////////////////////////////////////////////////////
    TGraphErrors *gr = new TGraphErrors(n, x_v, y_i, x_v_err, y_i_err); // x_freq_err=0
    gr->SetTitle("I(V) (595)");
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
	leg->AddEntry((TObject*)0, TString::Format("Range fit: [%.1f, %.qf]", min, max), "");
	leg->AddEntry((TObject*)0, TString::Format("m = %.3f #pm %.5f", m, m_err), "");
	leg->AddEntry((TObject*)0, TString::Format("m/#sigma_{m} = %.2f", fabs(m/m_err)), "");
	leg->AddEntry((TObject*)0, TString::Format("#alpha = %.1f (soglia)", alpha), "");
	
	leg->Draw();
}

