{
    // Canvas
    TCanvas *c1 = new TCanvas("c1", "Planck", 400, 400);
    c1->SetGrid();

    // DATA /////////////////////////////////////////////////////////////////////////////
    const int n = 3;

	double x_f[n] = {59.34, 100.59, 8.8328};    // frequenza in Hz
    double x_f_err[n] = {2.2114, 4.083,  4.4313};  // errore Hz

    double y_v[n] = {230, 401,  50};         // tensione in mV
    double y_v_err[n] = {16, 9, 5};            // errore mV

	// GRAPH /////////////////////////////////////////////////////////////////////////////
    TGraphErrors *gr = new TGraphErrors(n, x_f, y_v, x_f_err, y_v_err);
    gr->SetTitle("Photoelectric Effect;Frequency (Hz);Voltage (mV)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(1);
    gr->SetFillColor(kBlue - 10);
    gr->SetFillStyle(3001); // semi-transparent fill
    gr->Draw("AP");

    // FIT /////////////////////////////////////////////////////////////////////////////
    TF1 *ff1 = new TF1("ff1", "pol1", 5, 410);
    ff1->SetLineColor(kRed);
    ff1->SetLineWidth(2);

    gr->Fit("ff1", "SR+", "Q");  // S = silent, R = range, + = improved fit

    double intercetta = ff1->GetParameter(0);  // valore della frequenza a V=0
    double pendenza = ff1->GetParameter(1);    // coefficiente angolare (Hz/mV)

    cout << "Intercetta: " << intercetta << " Hz" << endl;
    cout << "Pendenza: " << pendenza << " Hz/mV" << endl;

    // Calcolo costante di Planck
    double e = 1.602e-19;  // carica elettrone in C
    double pendenza_SI = pendenza ;  // conversione mV → V
    double h = pendenza * e;           // h = e / m

    cout << "Stima di h = " << h*pow(10,-15) << " J·s" << endl;
}

