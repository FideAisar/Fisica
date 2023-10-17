{
// STANDARD CODE
// HISTOGRAMS TO BE READ FROM THE FILES
TH1D *histo1_a = new TH1D;
TH1D *histo1_b = new TH1D;
TH1D *histo2_a = new TH1D;
TH1D *histo2_b = new TH1D;
TH1D *histo3_a = new TH1D;
TH1D *histo3_b = new TH1D;
TH1D *histo4 = new TH1D;

// HISTOGRAMS TO BE PRODUCED
TH1D *histo_out_a = new TH1D("hew_"+S_plot+"_200",S_title,nbins,h_low,h_high);
TH1D *histo_out_b = new TH1D("hew_"+S_plot+"_noh",S_title,nbins,h_low,h_high);


// READ
histo1_a = (TH1D*)f_ew6_qcd0_200->Get(S_dir+"hew_"+S_plot);
histo1_b = (TH1D*)f_ew6_qcd0_noh->Get(S_dir+"hew_"+S_plot);
histo2_a = (TH1D*)f_ew4_qcd2_notop_200->Get(S_dir+"hew_"+S_plot);
histo2_b = (TH1D*)f_ew4_qcd2_notop_noh->Get(S_dir+"hew_"+S_plot);
histo3_a = (TH1D*)f_ew4_qcd2_top_200->Get(S_dir+"hew_"+S_plot);
histo3_b = (TH1D*)f_ew4_qcd2_top_noh->Get(S_dir+"hew_"+S_plot);
histo4 = (TH1D*)f_ew2_qcd4->Get(S_dir+"hew_"+S_plot);

// OPEN SUBDIR
  f->mkdir(S_plot+"_dir");
  f->cd(S_plot+"_dir");

// SUM THE VARIOUS CONTRIBUTIONS

Double_t dummyval;
for(Int_t i=0; i<=histo_out_a->GetNbinsX()+1; i++) {
   dummyval = 
    histo1_a->GetBinContent(i)
    + histo2_a->GetBinContent(i)
    + histo3_a->GetBinContent(i)
    + histo4->GetBinContent(i);
  histo_out_a->SetBinContent(i,dummyval);
}
for(Int_t i=0; i<=histo_out_b->GetNbinsX()+1; i++) {
   dummyval = 
    histo1_b->GetBinContent(i)
    + histo2_b->GetBinContent(i)
    + histo3_b->GetBinContent(i)
    + histo4->GetBinContent(i);
  histo_out_b->SetBinContent(i,dummyval);
}


// COLOURS 0:white 1:black
//         2:red 3:bright green 4:bright blue 5:yellow 6:bright pink
//         7:light blue 7:dull green 8:dull blue 13:dark grey 50:maroon 
//ASSIGN STANDARD COLOURS AND NAMES

histo1_a->SetLineColor(kRed);
histo1_a->SetLineWidth(3);
histo1_a->SetName("hew_"+S_plot+"_ew6_qcd0_200");

histo1_b->SetLineColor(kRed);
histo1_b->SetLineWidth(3);
histo1_b->SetName("hew_"+S_plot+"_ew6_qcd0_noh");


histo3_a->SetLineColor(kBlue);
histo3_a->SetLineWidth(3);
histo3_a->SetName("hew_"+S_plot+"_ew4_qcd2_top_200");

histo3_b->SetLineColor(kBlue);
histo3_b->SetLineWidth(3);
histo3_b->SetName("hew_"+S_plot+"_ew4_qcd2_top_noh");

histo2_a->SetLineColor(phaviolet);
histo2_a->SetLineWidth(3);
histo2_a->SetName("hew_"+S_plot+"_ew4_qcd2_notop_200");

histo2_b->SetLineColor(phaviolet);
histo2_b->SetLineWidth(3);
histo2_b->SetName("hew_"+S_plot+"_ew4_qcd2_notop_noh");


histo4->SetLineColor(13);
histo4->SetLineWidth(3);
histo4->SetName("hew_"+S_plot+"_ew2_qcd4");


histo_out_a->SetLineColor(kBlack);
histo_out_a->SetLineWidth(3);
histo_out_a->SetName("hew_"+S_plot+"_200");

histo_out_b->SetLineColor(kBlack);
histo_out_b->SetLineWidth(3);
histo_out_b->SetName("hew_"+S_plot+"_noh");


//SAVE HISTOGRAMS ON FILE HISTO.ROOT

histo1_a->Write();
histo1_b->Write();
histo3_a->Write();
histo3_b->Write();
histo2_a->Write();
histo2_b->Write();
histo4->Write();
histo_out_a->Write();
histo_out_b->Write();

  f->cd();

//OPTIONALLY COMPUTE THE INTEGRAL BETWEEN range _low and range_high
if(i_integral == "y"){
  TCanvas *c_xsec = new TCanvas("total_cross_section");

  TPaveText *text_xsec = new TPaveText(0.1,0.1,0.9,0.9);
  text_xsec.SetLabel("Integrated cross sections");
  s = "Total cross section for "+S_title+" in (";
  s += integral _low;
  s += ",";
  s += integral _high;
  s+=")";
  text_xsec.AddText(s.Data());
  text_xsec.AddText(" ");
  text_xsec.AddText("                mh200                   noh ");
  s= "ew         ";
  s += integro(histo1_a,integral _low,integral_high);
  s += "   ";
  s +=  integro(histo1_b,integral _low,integral_high);  
  text_xsec.AddText(s.Data());
  s= "QCD notop    ";
  s += integro(histo2_a,integral _low,integral_high);
  s += "   ";
  s += integro(histo2_b,integral _low,integral_high);  
  text_xsec.AddText(s.Data());
  s= "QCD top  ";
  s += integro(histo3_a,integral _low,integral_high);
  s += "   ";
  s += integro(histo3_b,integral _low,integral_high);  
  text_xsec.AddText(s.Data());
  s= "QCD V4j                 ";
  s += integro(histo4,integral _low,integral_high);  
  s += "               ";
  text_xsec.AddText(s.Data());

  text_xsec.Draw();
  c_xsec.Write();

}

//OPTIONALLY REBIN AND RESCALE HISTORAMS FOR BETTER PRESENTATION
if(rebin_scalefactor > 1){
histo1_a->Rebin(rebin_scalefactor);
histo1_a->Scale(1.0/rebin_scalefactor);
histo1_b->Rebin(rebin_scalefactor);
histo1_b->Scale(1.0/rebin_scalefactor);
histo2_a->Rebin(rebin_scalefactor);
histo2_a->Scale(1.0/rebin_scalefactor);
histo2_b->Rebin(rebin_scalefactor);
histo2_b->Scale(1.0/rebin_scalefactor);
histo3_a->Rebin(rebin_scalefactor);
histo3_a->Scale(1.0/rebin_scalefactor);
histo3_b->Rebin(rebin_scalefactor);
histo3_b->Scale(1.0/rebin_scalefactor);
histo4->Rebin(rebin_scalefactor);
histo4->Scale(1.0/rebin_scalefactor);
histo_out_a->Rebin(rebin_scalefactor);
histo_out_a->Scale(1.0/rebin_scalefactor);
histo_out_b->Rebin(rebin_scalefactor);
histo_out_b->Scale(1.0/rebin_scalefactor);
}

//OPTIONALLY NORMALIZE TO UNITY
if(inorm == "y"){

TH1D *histo1_a_norm = (TH1D*)histo1_a->Clone();
TH1D *histo1_b_norm = (TH1D*)histo1_b->Clone();
TH1D *histo2_a_norm = (TH1D*)histo2_a->Clone();
TH1D *histo2_b_norm = (TH1D*)histo2_b->Clone();
TH1D *histo3_a_norm = (TH1D*)histo3_a->Clone();
TH1D *histo3_b_norm = (TH1D*)histo3_b->Clone();
TH1D *histo4_norm = (TH1D*)histo4->Clone();

histo1_a_norm->SetTitle("Same distributions, normalized to unit area");
histo1_b_norm->SetTitle("Same distributions, normalized to unit area");
histo2_a_norm->SetTitle("Same distributions, normalized to unit area");
histo2_b_norm->SetTitle("Same distributions, normalized to unit area");
histo3_a_norm->SetTitle("Same distributions, normalized to unit area");
histo3_b_norm->SetTitle("Same distributions, normalized to unit area");
histo4_norm->SetTitle("Same distributions, normalized to unit area");

Double_t dummy;
dummy = integro(histo1_a);
if (dummy <= 0.0){
cout <<"integral of "<<histo1_a.GetName()<<" = 0 !"<<endl;
}else{
histo1_a_norm->Scale(1.0/dummy);}
dummy = integro(histo1_b);
if (dummy <= 0.0){
cout <<"integral of "<<histo1_b.GetName()<<" = 0 !"<<endl;
}else{
histo1_b_norm->Scale(1.0/dummy);}
dummy = integro(histo2_a);
if (dummy <= 0.0){
cout <<"integral of "<<histo2_a.GetName()<<" = 0 !"<<endl;
}else{
histo2_a_norm->Scale(1.0/dummy);}
dummy = integro(histo2_b);
if (dummy <= 0.0){
cout <<"integral of "<<histo2_b.GetName()<<" = 0 !"<<endl;
}else{
histo2_b_norm->Scale(1.0/dummy);}
dummy = integro(histo3_a);
if (dummy <= 0.0){
cout <<"integral of "<<histo3_a.GetName()<<" = 0 !"<<endl;
}else{
histo3_a_norm->Scale(1.0/dummy);}
dummy = integro(histo3_b);
if (dummy <= 0.0){
cout <<"integral of "<<histo3_b.GetName()<<" = 0 !"<<endl;
}else{
histo3_b_norm->Scale(1.0/dummy);}
dummy = integro(histo4);
if (dummy <= 0.0){
cout <<"integral of "<<histo4.GetName()<<" = 0 !"<<endl;
}else{
histo4_norm->Scale(1.0/dummy);}

TCanvas *c_a = new TCanvas("c_"+S_plot+"_200",S_title,10,10,700,1000);
c_a->SetBorderSize(2);
c_a->SetFrameFillColor(0);
c_a->SetFrameLineWidth(3);
c_a.Divide(1,2);
c_a.cd(1);


TLegend leg1(0.6,0.7,0.89,0.89);
leg1.SetMargin(0.5);
leg1.AddEntry(histo_out_a,"Total","l");
leg1.AddEntry(histo1_a,"#alpha_{ew}^{6}","l");
leg1.AddEntry(histo3_a,"#alpha_{ew}^{4} #alpha_{s}^{2} top","l");
leg1.AddEntry(histo2_a,"#alpha_{ew}^{4} #alpha_{s}^{2} notop","l");
leg1.AddEntry(histo4,"#alpha_{ew}^{2} #alpha_{s}^{4}","l");

histo_out_a->SetStats(kFALSE);
histo_out_a->GetXaxis()->SetTitle(S_Xaxis_title);
if(irange == "y"){histo_out_a->GetXaxis()->SetRangeUser(range_low,range_high);}
histo_out_a->GetYaxis()->SetTitle(S_Yaxis_title);

histo_out_a->Draw();

histo1_a->Draw("same");
histo3_a->Draw("same");
histo2_a->Draw("same");
histo4->Draw("same");

// REDRAW IN ORDER TO HAVE IT ON TOP
histo_out_a->Draw("same");

leg1.Draw();
c_a.Draw();

c_a.cd(2);  // WRITE NORMALIZED HISTOGRAMS

TLegend leg1_norm(0.6,0.7,0.89,0.89);
leg1_norm.SetMargin(0.5);
leg1_norm.AddEntry(histo1_a_norm,"#alpha_{ew}^{6}","l");
leg1_norm.AddEntry(histo3_a_norm,"#alpha_{ew}^{4} #alpha_{s}^{2} top","l");
leg1_norm.AddEntry(histo2_a_norm,"#alpha_{ew}^{4} #alpha_{s}^{2} notop","l");
leg1_norm.AddEntry(histo4_norm,"#alpha_{ew}^{2} #alpha_{s}^{4}","l");

histo4_norm->SetStats(kFALSE);
histo4_norm->GetXaxis()->SetTitle(S_Xaxis_title);
if(irange == "y"){histo4_norm->GetXaxis()->SetRangeUser(range_low,range_high);}
histo4_norm->GetYaxis()->SetTitle(S_Yaxis_title);

histo4_norm->Draw();
histo1_a_norm->Draw("same");
histo3_a_norm->Draw("same");
histo2_a_norm->Draw("same");

// REDRAW IN ORDER TO HAVE IT ON TOP
histo4_norm->Draw("same");

leg1_norm.Draw();
c_a.Draw();

c_a.Write();

// NOH+MH200_EW CANVAS 
TCanvas *c_b = new TCanvas("c_"+S_plot+"_noh",S_title,10,10,700,1000);
c_b->SetBorderSize(2);
c_b->SetFrameFillColor(0);
c_b->SetFrameLineWidth(3);
c_b.Divide(1,2);
c_b.cd(1);

TLegend leg1(0.6,0.7,0.89,0.89);
leg1.SetMargin(0.5);
leg1.AddEntry(histo_out_b,"Total","l");
leg1.AddEntry(histo1_b,"#alpha_{ew}^{6}","l");
histo1_a->SetLineColor(6);
leg1.AddEntry(histo3_b,"#alpha_{ew}^{4} #alpha_{s}^{2} top","l");
leg1.AddEntry(histo2_b,"#alpha_{ew}^{4} #alpha_{s}^{2} notop","l");
leg1.AddEntry(histo4,"#alpha_{ew}^{2} #alpha_{s}^{4}","l");
leg1.AddEntry(histo1_a,"#alpha_{ew}^{6} mh200","l");

histo_out_b->SetStats(kFALSE);
histo_out_b->GetXaxis()->SetTitle(S_Xaxis_title);
if(irange == "y"){histo_out_b->GetXaxis()->SetRangeUser(range_low,range_high);}
histo_out_b->GetYaxis()->SetTitle(S_Yaxis_title);

histo_out_b->Draw();

histo1_b->Draw("same");
histo3_b->Draw("same");
histo2_b->Draw("same");
histo4->Draw("same");
histo1_a->Draw("same");

// REDRAW IN ORDER TO HAVE IT ON TOP
histo_out_b->Draw("same");

leg1.Draw();
c_b.Draw();

leg1.Draw();
c_b.Draw();

c_b.cd(2);  // WRITE NORMALIZED HISTOGRAMS

TLegend leg1_norm(0.6,0.7,0.89,0.89);
leg1_norm.SetMargin(0.5);
leg1_norm.AddEntry(histo1_b_norm,"#alpha_{ew}^{6}","l");
histo1_a_norm->SetLineColor(6);
leg1_norm.AddEntry(histo3_b_norm,"#alpha_{ew}^{4} #alpha_{s}^{2} top","l");
leg1_norm.AddEntry(histo2_b_norm,"#alpha_{ew}^{4} #alpha_{s}^{2} notop","l");
leg1_norm.AddEntry(histo4_norm,"#alpha_{ew}^{2} #alpha_{s}^{4}","l");
leg1_norm.AddEntry(histo1_a_norm,"#alpha_{ew}^{6} mh200","l");

histo4_norm->SetStats(kFALSE);
histo4_norm->GetXaxis()->SetTitle(S_Xaxis_title);
if(irange == "y"){histo4_norm->GetXaxis()->SetRangeUser(range_low,range_high);}
histo4_norm->GetYaxis()->SetTitle(S_Yaxis_title);

histo4_norm->Draw();
histo1_b_norm->Draw("same");
histo3_b_norm->Draw("same");
histo2_b_norm->Draw("same");
histo1_a_norm->Draw("same");

// REDRAW IN ORDER TO HAVE IT ON TOP
histo4_norm->Draw("same");

leg1_norm.Draw();
c_b.Draw();

c_b.Write();

}
else{  // DO NOT NORMALIZE

TCanvas *c_a = new TCanvas("c_"+S_plot+"_200",S_title,10,10,700,600);
c_a->SetBorderSize(2);
c_a->SetFrameFillColor(0);
c_a->SetFrameLineWidth(3);

TLegend leg1(0.6,0.7,0.89,0.89);
leg1.SetMargin(0.5);
leg1.AddEntry(histo_out_a,"Total","l");
leg1.AddEntry(histo1_a,"#alpha_{ew}^{6}","l");
leg1.AddEntry(histo3_a,"#alpha_{ew}^{4} #alpha_{s}^{2} top","l");
leg1.AddEntry(histo2_a,"#alpha_{ew}^{4} #alpha_{s}^{2} notop","l");
leg1.AddEntry(histo4,"#alpha_{ew}^{2} #alpha_{s}^{4}","l");

histo_out_a->SetStats(kFALSE);
histo_out_a->GetXaxis()->SetTitle(S_Xaxis_title);
if(irange == "y"){histo_out_a->GetXaxis()->SetRangeUser(range_low,range_high);}
histo_out_a->GetYaxis()->SetTitle(S_Yaxis_title);

histo_out_a->Draw();

histo1_a->Draw("same");
histo3_a->Draw("same");
histo2_a->Draw("same");
histo4->Draw("same");

// REDRAW IN ORDER TO HAVE IT ON TOP
histo_out_a->Draw("same");

leg1.Draw();
c_a.Draw();

c_a.Write();


TCanvas *c_b = new TCanvas("c_"+S_plot+"_noh",S_title,10,10,700,600);
c_b->SetBorderSize(2);
c_b->SetFrameFillColor(0);
c_b->SetFrameLineWidth(3);

TLegend leg1(0.6,0.7,0.89,0.89);
leg1.SetMargin(0.5);
leg1.AddEntry(histo_out_b,"Total","l");
histo1_a->SetLineColor(6);
leg1.AddEntry(histo1_a,"#alpha_{ew}^{6} mh200","l");
leg1.AddEntry(histo1_b,"#alpha_{ew}^{6}","l");
leg1.AddEntry(histo3_b,"#alpha_{ew}^{4} #alpha_{s}^{2} top","l");
leg1.AddEntry(histo2_b,"#alpha_{ew}^{4} #alpha_{s}^{2} notop","l");
leg1.AddEntry(histo4,"#alpha_{ew}^{2} #alpha_{s}^{4}","l");

histo_out_b->SetStats(kFALSE);
histo_out_b->GetXaxis()->SetTitle(S_Xaxis_title);
if(irange == "y"){histo_out_b->GetXaxis()->SetRangeUser(range_low,range_high);}
histo_out_b->GetYaxis()->SetTitle(S_Yaxis_title);

histo_out_b->Draw();

histo1_b->Draw("same");
histo3_b->Draw("same");
histo2_b->Draw("same");
histo4->Draw("same");
histo1_a->Draw("same");

// REDRAW IN ORDER TO HAVE IT ON TOP
histo_out_b->Draw("same");

leg1.Draw();
c_b.Draw();

c_b.Write();
}
// END STANDARD CODE

}
