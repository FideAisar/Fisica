{
    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Planck", 400, 400);
    c1->SetGrid();
    

    // DATA /////////////////////////////////////////////////////////////////////////////
    const int n = 32;
	double x_v[n] = { 0.70,20.9,41.2,60.0,82.6,101.0,121.5,143.4,161.0,182.2,202.8,224.8,241.1,263.5,284.1,301.8,349.8,402.3,454.7,501.7,551.7,602.1,702.0,804.1,900.6,1003.0,1101.0,1200.0,1501.0,2001.0,2540.0,2991.0 };
        
    double x_v_err[n] = { 1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,5,5,6,6,7,12,13,14,15,16,17,20,25,30,35 };

	double y_i[n] = {  0.00305,0.00208,0.00136,0.00085,0.00046,0.00013,-0.00009,-0.00025,-0.00035,-0.00049,-0.00059,-0.00069,-0.00077,-0.00085,-0.00089,-0.00097,-0.00125,-0.00145,-0.00161,-0.00173,-0.00183,-0.00193,-0.00208,-0.00223,-0.00234,-0.00240,-0.00248,-0.00255,-0.00267,-0.00290,-0.00300,-0.00320};

	double y_i_err[n] = { 0.00030,0.00010,0.00010,0.00040,0.00020,0.00010,0.00030,0.00010,0.00010,0.00010,0.00030,0.00030,0.00030,0.00002,0.00006,0.00002,0.00004,0.00002,0.00002,0.00002,0.00001,0.00001,0.00002,0.00002,0.00004,0.00003,0.00002,0.00002,0.00002,0.00020,0.00050,0.00050 };


     // GRAPH /////////////////////////////////////////////////////////////////////////////
    TGraphErrors *gr = new TGraphErrors(n, x_v, y_i, x_v_err, y_i_err); // x_freq_err=0
    gr->SetTitle("I(V) (644)");
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
	leg->AddEntry((TObject*)0, TString::Format("Range fit: [%.1f, %.f]", min, max), "");
	leg->AddEntry((TObject*)0, TString::Format("m = %.3f #pm %.3f", m, m_err), "");
	leg->AddEntry((TObject*)0, TString::Format("m/#sigma_{m} = %.2f", fabs(m/m_err)), "");
	leg->AddEntry((TObject*)0, TString::Format("#alpha = %.1f (soglia)", alpha), "");
	
	leg->Draw();
}

