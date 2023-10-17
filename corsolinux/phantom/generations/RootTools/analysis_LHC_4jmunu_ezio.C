/*
 *analisys_LHC_4jmunu.C
 *
 *Last update: Aug 31, 2007
 *
 ******************************************************************************/

#if !defined(_CINT_) || defined(_MAKECINT_)
#include <Riostream.h>
#include <stdio.h>
#include "TROOT.h"
#include "TString.h"
#include "TObjString.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TMath.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TList.h"
#include <TAxis.h>
#endif
using namespace std;


//void analysis_LHC_4jmunu(TString, TString);
void analysis_LHC_4jmunu(TString, Double_t);
double integro(TH1D*);
double integro(TH1D*,double,double);


//void analysis_LHC_4jmunu(TString treefilename, TString opt) {
void analysis_LHC_4jmunu(TString opt, Double_t rmh) {


  gSystem->Load("libPhysics.so");
  gROOT->Reset();


  //total cross section for histogram normalization
  Double_t xsec = 0.0;
  Double_t xsec_tot = 0.0;

  //parameters for histogram rescaling
  Double_t scale_Mlv, scale_Mjj, scale_Mjjlv, scale_eta, scale_pT, scale_phi;;
  static Double_t A_Mlv, A_Mjj, A_Mjjlv, A_eta, A_pT, A_phi;

  //event kinematics
  Double_t px[8], py[8], pz[8], E[8];
  Int_t idup[8],istup[8];
  TLorentzVector particle[8], jet[8], clep[8], nu[8];
  TLorentzVector jf, jb, jc[2], jj, jjj, jlv, jfjb, V[2], Vlep_reco[2], 
    Vlep_reco_jf[2], Vlep_reco_jb[2], Vhad_jf, Vhad_jb, Vlep_jf, Vlep_jb, VV;
  Int_t n_j, n_clep, n_nu, index_jet[8], index_b[8];

  //parameters for histogram settings
  Double_t Mlv_high, Mlv_low;
  Double_t Mjj_high, Mjj_low;
  Double_t Mjjlv_high, Mjjlv_low;
  Double_t eta_high, eta_low;
  Double_t pT_high, pT_low; 
  Double_t phi_high, phi_low;
  Int_t nentries_tree=0;   //total number of entries of the selected TREE(S)
  Int_t nbins_M;           //number of bins for invariant mass histograms
  Int_t nbins_eta;         //number of bins for pseudorapidity histograms
  Int_t nbins_pT;          //number of bins for pT histograms
  Int_t nbins_phi;         //number of bins for phi histograms

  //parameters for neutrino reconstruction
  Double_t pL_nu1st=0.0;
  Double_t pL_nu2nd=0.0;
  Double_t E_nu1st=0.0;
  Double_t E_nu2nd=0.0;
  Double_t auxa, auxb, auxdelta;
  Double_t Mw2;

  TLorentzVector p_bb;        //momentum of bb system
  TLorentzVector p_bj;        //momentum of bj system
  Int_t is_a_b;
  Double_t M_bj[3];        //invariant mass of single tagged b with another jet
                           //system (3 in total)
  Double_t pt_bj[3];       //pt single tagged b with another jet system 
  Double_t DR;
  
  //user-defined selection cuts
  TString setcuts_E_min_j, setcuts_E_min_clep, setcuts_pT_min_j, 
    setcuts_pT_min_clep, setcuts_eta_max_j, setcuts_eta_max_clep, 
    setcuts_Mjj, setcuts_Mjcjc, setcuts_deltaetajfjb, setcuts_Wlv, 
    setBveto, setBsingletag, setBdoubletag, setTopcut_jjj, setTopcut_jlv_rec,
    setcuts_DRbbj;
  Double_t cut_E_min_j, cut_E_min_clep, cut_pT_min_j, cut_pT_min_clep, 
    cut_eta_max_j, cut_eta_max_clep, centralvalue_Mjcjc, cut_Mjj, cut_Mjcjc, 
    cut_deltaetajfjb, cut_pT_Wlv, cut_eta_Wlv_rec, cut_M_lvj_min_rec, 
    Btagging_region, Btagging_efficiency, 
    centralvalue_Mjjj, cut_Mjjj, centralvalue_Mjlv_rec, cut_Mjlv_rec, Mcut,
    cut_DRbbj;
  //random number generator and the corresponding random variable
  TRandom3 irand;

  static Int_t ipasscuts;
  static Int_t ib, ibtag;

  Int_t iwarning=0;


/*******************************************************************************
 * Beginning of USER-DEFINED SELECTION CUTS ("y"=cut applied; "n"=cut not 
 * applied)
 ******************************************************************************/

  //Minimum energy for jets
  cin>>setcuts_E_min_j;
//  setcuts_E_min_j="n";
  if(setcuts_E_min_j=="y"){
    cin>>cut_E_min_j;
//  cut_E_min_j = 20.0;
  }
  //Minimum energy for charged leptons
  setcuts_E_min_clep="n";
  cut_E_min_clep = 20.0;

  //Minimum pT for jets
  setcuts_pT_min_j="y";
  cut_pT_min_j = 30.0;

  //Minimum pT for charged leptons
  setcuts_pT_min_clep="y";
  cut_pT_min_clep = 20.0;

  //Maximum pseudorapidity for jets
  setcuts_eta_max_j="y";
  cut_eta_max_j=6.5;

  //Maximum pseudorapidity for charged  leptons
  setcuts_eta_max_clep="y";
  cut_eta_max_clep=3.0;

  //Mimimum invariant mass of any two jets
  setcuts_Mjj="y";
  cut_Mjj = 60.0;

  //Cut on invariant mass of the two central jets
  setcuts_Mjcjc="n";
  centralvalue_Mjcjc = 85.0; //GeV
  cut_Mjcjc = 15.0; // => window = centralMjj +/- cutMjj

  //Cut on the difference in pseudorapidity between forward and backward tag jet
  setcuts_deltaetajfjb="n";
  cut_deltaetajfjb = 4.0;

  //Cuts on the leptonic W
  setcuts_Wlv="n";
  cut_pT_Wlv = 100.0;
  cut_eta_Wlv_rec = 2.0;
  cut_M_lvj_min_rec = 250.0;

  //Cut on minimum DeltaR of bb pair with a jet 
  setcuts_DRbbj="y";
  cut_DRbbj = 2.0; 

  //Top cuts:
  //1.Cut on invariant mass of jjj
  setTopcut_jjj="y";
  centralvalue_Mjjj = 175.0; //GeV
  cut_Mjjj = 15.0; // => window = centralMjjj +/- cutMjjj

  //2.Cut on invariant mass of jlv_rec
  setTopcut_jlv_rec="y";
  centralvalue_Mjlv_rec = 175.0; //GeV
  cut_Mjlv_rec = 15.0; // => window = centralMjlv +/- cutMjlv

  Mcut = 0.0;
  
  //Mimic b-quark veto or b-tagging
  Btagging_region = 1.5;  // B-tagging active on |eta| < Btagging_region
  Btagging_efficiency=0.8;
  setBveto="n";           // rejects b-tagged events
  setBsingletag="y";      // accepts single b-tagged events only
  setBdoubletag="y";      // accepts double b-tagged events only

/*******************************************************************************
 * End of USER-DEFINED SELECTION CUTS
 ******************************************************************************/



  //Information about the applied selection cuts will be stored in a canvas
  TPaveText text_cuts(0.1,0.1,0.9,0.9);
  text_cuts.SetLabel("List of selection cuts applied");
  text_cuts.AddText(" ");
  if(setcuts_E_min_j=="y"){
    TString s="minimum energy of jets:";
    //s.Form("%f",cut_E_min_j);
    s+=cut_E_min_j;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setcuts_E_min_clep=="y"){
    TString s="minimum energy of charged leptons:";
    s+=cut_E_min_clep;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setcuts_pT_min_j=="y"){
    TString s="minimum p_{T} of jets:";
    s+=cut_pT_min_j;
    s+=" GeV";  
    text_cuts.AddText(s.Data());
  }
  if(setcuts_pT_min_clep=="y"){
    TString s="minimum p_{T} of charged leptons:";
    s+=cut_pT_min_clep;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setcuts_eta_max_j=="y"){
    TString s="maximum |#eta| of jets:";
    s+=cut_eta_max_j;
    text_cuts.AddText(s.Data());
  }
  if(setcuts_eta_max_clep=="y"){
    TString s="maximum |#eta| of charged leptons:";
    s+=cut_eta_max_clep;
    text_cuts.AddText(s.Data());
  }
  if(setcuts_Mjj=="y"){
    TString s="minimum invariant mass of any jet pair:";
    s+=cut_Mjj;
    s+=" GeV";    
    text_cuts.AddText(s.Data());
  }
  if(setcuts_Mjcjc=="y"){
    TString s="inv.mass window for the 2 central jets:";
    s+=centralvalue_Mjcjc;
    s+="      +/-";
    s+=cut_Mjcjc;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setcuts_deltaetajfjb=="y"){
    TString s="minimum |#Delta#eta| between tag jets:";
    s+=cut_deltaetajfjb;
    text_cuts.AddText(s.Data());
  }
  if(setcuts_Wlv=="y"){
    TString s="minimum p_{T} of leptonic W (reconstructed):";
    s+=cut_pT_Wlv;
    s+=" GeV";
    text_cuts.AddText(s.Data());
    s="maximum |#eta| of leptonic W (reconstructed):";
    s+=cut_eta_Wlv_rec;
    text_cuts.AddText(s.Data());
    s="minimum inv.mass M_{jl#nu} (reconstructed):";
    s+=cut_M_lvj_min_rec;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setBveto=="y"){
    TString s="b-quark veto active for |#eta|    <";
    s+=Btagging_region;
    s+="     with efficiency:";
    s+=Btagging_efficiency;
    text_cuts.AddText(s.Data());
  }
  if(setBsingletag=="y"){
    TString s="single b-quark tagging active for |#eta|    <";
    s+=Btagging_region;
    s+="     with efficiency:";
    s+=Btagging_efficiency;
    text_cuts.AddText(s.Data());
  }
  if(setBdoubletag=="y"){
    TString s="double b-quark tagging active for |#eta|    <";
    s+=Btagging_region;
    s+="     with efficiency:";
    s+=Btagging_efficiency;
    text_cuts.AddText(s.Data());
  }
  if(setcuts_DRbbj=="y"){
    TString s="minimum DelatR between bb system and jet:";
    s+=cut_DRbbj;
    text_cuts.AddText(s.Data());
  }
  if(setTopcut_jjj=="y"){
    TString s="top-veto n.1:     |M_{jjj} - M_{top}|     >";
    s+=cut_Mjjj;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setTopcut_jlv_rec=="y"){
    TString s="top-veto n.2:     |M_{jl#nu}(reco) - M_{top}|     >";
    s+=cut_Mjlv_rec;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  text_cuts.AddText(" ");





  cout<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<"                  -- analysis_LHC_4jmunu --                   "<<endl;
  cout<<"            selection procedure for pp -> muvm_4j             "<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<endl;

  cout<<"Enter the number of tree files to be considered (cross sections will be summed): "<<endl;
  Int_t ntr=0;
  cin>>ntr;
  TString treefilename[ntr];
  for(Int_t ki=0;ki<ntr;ki++){
    cout<<"Enter name of tree file no."<<ki+1<<": "<<endl;
    cin>>treefilename[ki];
  }

  //User can choose between two possible schemes for identifying TAG JETS:
  // "eta" scheme  --  tag jets = most forward/backward jets
  // "pT"  scheme  --  tag jets = jets with highest pT
  //
  Int_t itag_scheme = 0;
  cout<<endl<<"SELECT SCHEME FOR IDENTIFYING TAG JETS"<<endl;
  cout<<" eta scheme: tag jets = most forward/backward jets (enter 1)"<<endl;
  cout<<" pT scheme : tag jets = jets with highest pT       (enter 2)"<<endl;
  cout<<" >> ";
  cin>>itag_scheme;
  cout<<endl;
  cout<<"Processing. Please wait..."<<endl;
  cout<<endl;
  
  if(itag_scheme==1){
    text_cuts.AddText(" ");
    text_cuts.AddText("Tag jets identified through eta-scheme");
    text_cuts.AddText(" ");
  }
  else if(itag_scheme==2){
    text_cuts.AddText(" ");
    text_cuts.AddText("Tag jets identified through pT-scheme");
    text_cuts.AddText(" ");
  }




  // Create output file plots.root_________________________________
  TFile *f = new TFile("plots.root","RECREATE");
  //create subdirectories for single histograms
  f->mkdir("MVV_allevents");
  f->mkdir("MVV_cuts");
  f->mkdir("otherhistos_cuts");
  f->cd();
  f->Close();

  TString title;
  char app[50];
  Int_t val = (Int_t)Mcut;
  sprintf(app,"(for M_{VV} > %d GeV) ",val);



  // Set the number of bins and extrema for each histogram____________________
  //cout<<"Enter number of bins for invariant mass histograms: ";
  //cin>>nbins_M;
  nbins_M = 14000;
  nbins_eta=100;
  nbins_pT=1400;
  nbins_phi=nbins_eta;


  //set extrema for histograms________________________________________________
  Mlv_low = 0;
  Mlv_high = 14000;
  Mjj_low = 0;
  Mjj_high = 14000;

  //cout<<"Enter Mjjlv_low: ";
  //cin>>Mjjlv_low;
  //cout<<endl;
  //cout<<"Enter Mjjlv_high: ";
  //cin>>Mjjlv_high;
  Mjjlv_low=5.0;
  Mjjlv_high=14005.0;

  eta_low = -7;
  eta_high = 7;
  pT_low = 5;
  pT_high = 14005;
  phi_low = 0.0;
  phi_high = TMath::Pi();


  // Declaration of histograms__________________________________________________

  // Before cuts
  //-------------
  TH1D hew_Mlv_allproc("hew_Mlv_allproc",";M_{l#nu} (GeV)",nbins_M,Mlv_low,Mlv_high);
  hew_Mlv_allproc.SetTitle("Invariant mass of l#nu system with no selection cut");

  TH1D hew_Mjj_allproc("hew_Mjj_allproc",";M_{jj} (GeV)",nbins_M,Mjj_low,Mjj_high);
  hew_Mjj_allproc.SetTitle("Invariant mass of j_{c}j_{c} system with no selection cut");

  TH1D hew_Mjjlv_allproc("hew_Mjjlv_allproc","M_{jjl#nu}",nbins_M,Mjjlv_low,Mjjlv_high);
  hew_Mjjlv_allproc.SetTitle("M_{jjl#nu} with no selection cut;M_{jjl#nu} (GeV);");

  TH1D hew_Mbb_allproc("hew_Mbb_allproc",";M_{bb}",nbins_M,Mjj_low,Mjj_high);
  hew_Mbb_allproc.SetTitle("Invariant mass of bb system (pb/GeV)");

  // At the Higgs peak
  //-------------
  // delta-eta between tag jets
  TH1D hew_DEta_jfjb_peak("hew_DEta_jfjb_peak","|#eta|",nbins_eta,0,2*eta_high);
  title="Difference in pseudorapidity between tag jets at the Higgs peak";
  hew_DEta_jfjb_peak.SetTitle(title.Data());

  TH1D hew_pTlv_peak("hew_pTlv_peak","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT of leptonic W at the Higgs peak (pb/GeV)";
  hew_pTlv_peak.SetTitle(title.Data());

  // After cuts
  //-------------
  TH1D hew_Mlv_cuts("hew_Mlv_cuts",";M_{l#nu} (GeV)",nbins_M,Mlv_low,Mlv_high);
  hew_Mlv_cuts.SetTitle("Invariant mass of l#nu system with selection cuts");

  TH1D hew_Mjj_cuts("hew_Mjj_cuts",";M_{jj} (GeV)",nbins_M,Mjj_low,Mjj_high);
  hew_Mjj_cuts.SetTitle("Invariant mass of j{_c}j{_c} system with selection cuts");

  TH1D hew_Mjjlv_cuts("hew_Mjjlv_cuts","M_{jjl#nu}",nbins_M,Mjjlv_low,Mjjlv_high);
  hew_Mjjlv_cuts.SetTitle("M_{jjl#nu} with selection cuts;M_{jjl#nu} (GeV);");

  TH1D hew_Mbb_cuts("hew_Mbb_cuts",";M_{bb}",nbins_M,Mjj_low,Mjj_high);
  hew_Mbb_cuts.SetTitle("Invariant mass of bb system  with selection cuts (pb/GeV)");

  TH1D hew_Mbj1_cuts("hew_Mbj1_cuts",";M_{bj1}",nbins_M,Mjj_low,Mjj_high);
  hew_Mbj1_cuts.SetTitle("Invariant mass of bj1 system with selection cuts (pb/GeV)");

  TH1D hew_Mbj2_cuts("hew_Mbj2_cuts",";M_{bj2}",nbins_M,Mjj_low,Mjj_high);
  hew_Mbj2_cuts.SetTitle("Invariant mass of bj2 system with selection cuts (pb/GeV)");

  TH1D hew_Mbj3_cuts("hew_Mbj3_cuts",";M_{bj3}",nbins_M,Mjj_low,Mjj_high);
  hew_Mbj3_cuts.SetTitle("Invariant mass of bj3 system with selection cuts (pb/GeV)");
  //|Delta R| between b's

 // At the Higgs peak
  //-------------
  TH1D hew_pTbb_peak("hew_pTbb_peak","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT of bb system at the Higgs peak with selection cuts (pb/GeV)";
  hew_pTbb_peak.SetTitle(title.Data());

  TH1D hew_pTbj_peak("hew_pTbj_peak","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT of minimum mass bj system at the Higgs peak with selection cuts (pb/GeV)";
  hew_pTbj_peak.SetTitle(title.Data());

  TH1D hew_DR_bb_peak("hew_DR_bb_peak","#Delta #R ",nbins_eta,0,2*eta_high);
  title="#Delta #R between b's  at the Higgs peak with selection cuts (GeV^{-1})"+TString(app);
  hew_DR_bb_peak.SetTitle(title.Data());
 
  TH1D hew_DRmin_bbj_peak("hew_DRmin_bbj_peak","#Delta #R ",nbins_eta,0,2*eta_high);
  title="Minumum #Delta R between bb system and a jet at the Higgs peak with selection cuts (GeV^{-1})"+TString(app);
  hew_DRmin_bbj_peak.SetTitle(title.Data());
 
//  TH1D hew_DR_bj_peak("hew_DR_bj_peak","#Delta #R ",nbins_eta,0,2*eta_high);
//  title="#Delta #R between b and jet at the Higgs peak with selection cuts (GeV^{-1}"+TString(app);
//  hew_DR_bj_peak.SetTitle(title.Data());
 

  // Control plots
  //-------------
  //1. minimum pT of tag jets
  TH1D hew_pTmin_jfjb("hew_pTmin_jfjb","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="Minimum pT of tag jets "+TString(app);
  hew_pTmin_jfjb.SetTitle(title.Data());

  //2. minimum pT of all jets
  TH1D hew_pTmin_alljets("hew_pTmin_alljets","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="Minimum pT of all jets "+TString(app);
  hew_pTmin_alljets.SetTitle(title.Data());

  //3. pT of heavy boson decaying hadronically (associated to central jets)
  TH1D hew_pT_jcjc("hew_pT_jcjc","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT of the two central jets "+TString(app);
  hew_pT_jcjc.SetTitle(title.Data());

  //4. eta of the heavy boson decaying hadronically (associated to central jets)
  TH1D hew_eta_jcjc("hew_eta_jcjc","#eta",nbins_eta,eta_low,eta_high);
  title="Pseudorapidity of the two central jets "+TString(app);
  hew_eta_jcjc.SetTitle(title.Data());

  //5. invariant mass of the two tag jets
  TH1D hew_M_jfjb("hew_M_jfjb","M_{jj} (GeV)",nbins_M,Mjjlv_low,Mjjlv_high);
  title="Invariant mass of the two tag jets "+TString(app);
  hew_M_jfjb.SetTitle(title.Data());

  //6. min. inv. mass of the heavy boson decaying hadronically and one tag jet
  TH1D hew_Mmin_jcjcjfb("hew_Mmin_jcjcjfb","M_{jjj} (GeV)",nbins_M,Mjjlv_low,Mjjlv_high);
  title="Minimum inv. mass of the two central jets and one tag jet "+TString(app);
  hew_Mmin_jcjcjfb.SetTitle(title.Data());

  //7. |Delta eta| between tag jets
  TH1D hew_DeltaEta_jfjb("hew_DeltaEta_jfjb","|#eta|",nbins_eta,0,2*eta_high);
  title="|#Delta #eta| between tag jets "+TString(app);
  hew_DeltaEta_jfjb.SetTitle(title.Data());

  //8. |Delta eta| between reconstructed heavy bosons
  TH1D hew_DeltaEta_VV("hew_DeltaEta_VV","|#eta|",nbins_eta,0,2*eta_high);
  title="|#Delta #eta| between reconstructed heavy bosons "+TString(app);
  hew_DeltaEta_VV.SetTitle(title.Data());

  //9. maximum pT of tag jets
  TH1D hew_MaxPT_tagj("hew_MaxPT_tagj","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="Maximum pT of tag jets "+TString(app);
  hew_MaxPT_tagj.SetTitle(title.Data());

  //10. minimum inv. mass of tag jet + vector boson
  TH1D hew_MinM_Vj("hew_MinM_Vj","M_{jV} (GeV)",nbins_M,Mjjlv_low,Mjjlv_high);
  title="Minimum inv. mass of tag jet + vector boson "+TString(app);
  hew_MinM_Vj.SetTitle(title.Data()); 

  //10bis. complementary inv. mass of tag jet + vector boson
  //    [i.e.: data la minima massa tag-jet/bosone, la massa dell'altro tag-jet 
  //     con l'altro bosone]
  TH1D hew_comp_MinM_Vj("hew_comp_MinM_Vj","M_{jV} (GeV)",nbins_M,Mjjlv_low,Mjjlv_high);
  title="Complementary inv. mass of tag jet + vector boson "+TString(app);
  hew_comp_MinM_Vj.SetTitle(title.Data()); 

  //11. minimum delta-eta between tag jet and the vector boson closest to it
  TH1D hew_MinDEta_Vtag("hew_MinDEta_Vtag","|#Delta #eta|",nbins_eta,0,2*eta_high);
  title="Minimum |#Delta #eta| between tag jet and closest vector boson "+TString(app);
  hew_MinDEta_Vtag.SetTitle(title.Data());
  

  //11bis. complementary delta-eta between the tag jet and closest vector boson
  TH1D hew_comp_MinDEta_Vtag("hew_comp_MinDEta_Vtag","|#Delta #eta|",nbins_eta,0,2*eta_high);
  title="Complementary |#Delta #eta| between tag jet and closest vector boson "+TString(app);
  hew_comp_MinDEta_Vtag.SetTitle(title.Data());
  

  //12. |Delta phi| between tag jets
  TH1D hew_DPhi_tag("hew_DPhi_tag","|#Delta #phi| ",nbins_phi,phi_low,phi_high);
  title="|#Delta #phi| between tag jets "+TString(app);
  hew_DPhi_tag.SetTitle(title.Data());
  

  //13. |Delta phi| between reconstructed vector bosons
  TH1D hew_DPhi_V("hew_DPhi_V","|#Delta #phi| ",nbins_phi,phi_low,phi_high);
  title="|#Delta #phi| between reconstructed vector bosons "+TString(app);
  hew_DPhi_V.SetTitle(title.Data());


  hew_Mlv_allproc.SetLineColor(kBlue);
  hew_Mjj_allproc.SetLineColor(kBlue);
  hew_Mbb_allproc.SetLineColor(kBlue);
  hew_Mjjlv_allproc.SetLineColor(kBlue);
  hew_Mlv_cuts.SetLineColor(kBlue);
  hew_Mjj_cuts.SetLineColor(kBlue);
  hew_Mjjlv_cuts.SetLineColor(kBlue);
  hew_Mjjlv_allproc.SetLineColor(kBlue);
  hew_Mbb_cuts.SetLineColor(kBlue);





  //Define "buffer" histograms (temporary objects)
  TH1D *hew_Mlv_allproc_buf = (TH1D*) hew_Mlv_allproc.Clone();
  TH1D *hew_Mjj_allproc_buf = (TH1D*) hew_Mjj_allproc.Clone();
  TH1D *hew_Mjjlv_allproc_buf = (TH1D*) hew_Mjjlv_allproc.Clone();
  TH1D *hew_Mbb_allproc_buf = (TH1D*) hew_Mbb_allproc.Clone();
  TH1D *hew_Mlv_cuts_buf = (TH1D*) hew_Mlv_cuts.Clone();
  TH1D *hew_Mjj_cuts_buf = (TH1D*) hew_Mjj_cuts.Clone();
  TH1D *hew_Mjjlv_cuts_buf = (TH1D*) hew_Mjjlv_cuts.Clone();
  TH1D *hew_Mbb_cuts_buf = (TH1D*) hew_Mbb_cuts.Clone();
  TH1D *hew_Mbj1_cuts_buf = (TH1D*) hew_Mbj1_cuts.Clone();
  TH1D *hew_Mbj2_cuts_buf = (TH1D*) hew_Mbj2_cuts.Clone();
  TH1D *hew_Mbj3_cuts_buf = (TH1D*) hew_Mbj3_cuts.Clone();
  TH1D *hew_pTbb_peak_buf = (TH1D*) hew_pTbb_peak.Clone();
  TH1D *hew_pTbj_peak_buf = (TH1D*) hew_pTbj_peak.Clone();
  TH1D *hew_pTlv_peak_buf = (TH1D*) hew_pTlv_peak.Clone();
  TH1D *hew_DR_bb_peak_buf = (TH1D*) hew_DR_bb_peak.Clone();
  TH1D *hew_DRmin_bbj_peak_buf = (TH1D*) hew_DRmin_bbj_peak.Clone();
  TH1D *hew_DEta_jfjb_peak_buf = (TH1D*) hew_DEta_jfjb_peak.Clone();
  TH1D *hew_pTmin_jfjb_buf = (TH1D*) hew_pTmin_jfjb.Clone();
  TH1D *hew_pTmin_alljets_buf = (TH1D*) hew_pTmin_alljets.Clone();
  TH1D *hew_pT_jcjc_buf = (TH1D*) hew_pT_jcjc.Clone();
  TH1D *hew_eta_jcjc_buf = (TH1D*) hew_eta_jcjc.Clone();
  TH1D *hew_M_jfjb_buf = (TH1D*) hew_M_jfjb.Clone();
  TH1D *hew_Mmin_jcjcjfb_buf = (TH1D*) hew_Mmin_jcjcjfb.Clone();
  TH1D *hew_DeltaEta_jfjb_buf = (TH1D*) hew_DeltaEta_jfjb.Clone();
  TH1D *hew_DeltaEta_VV_buf = (TH1D*) hew_DeltaEta_VV.Clone();
  TH1D *hew_MaxPT_tagj_buf = (TH1D*) hew_MaxPT_tagj.Clone();
  TH1D *hew_MinM_Vj_buf = (TH1D*) hew_MinM_Vj.Clone();
  TH1D *hew_comp_MinM_Vj_buf = (TH1D*) hew_comp_MinM_Vj.Clone();
  TH1D *hew_MinDEta_Vtag_buf = (TH1D*) hew_MinDEta_Vtag.Clone();
  TH1D *hew_comp_MinDEta_Vtag_buf = (TH1D*) hew_comp_MinDEta_Vtag.Clone();
  TH1D *hew_DPhi_tag_buf = (TH1D*) hew_DPhi_tag.Clone();
  TH1D *hew_DPhi_V_buf = (TH1D*) hew_DPhi_V.Clone();








  for(Int_t ktr=0;ktr<ntr;ktr++){

  //Open the tree file and get the tree
  TFile treefile(treefilename[ktr].Data());
  TTree *tree = (TTree*) treefile.Get("events");


  // Take info from tree to fill histograms_____________________________________
  tree->SetBranchAddress("px",px);
  tree->SetBranchAddress("py",py);
  tree->SetBranchAddress("pz",pz);
  tree->SetBranchAddress("E",E);
  tree->SetBranchAddress("IDUP",idup);
  tree->SetBranchAddress("ISTUP",istup);
  if(opt=="xs") {
    // Take the value of total cross section (stored inside the tree as User
    // Information)
    TObjString *os_xsec = (TObjString*)tree->GetUserInfo()->At(0);
    TString s_xsec = os_xsec->GetString();
    xsec = s_xsec.Atof();
    xsec_tot += xsec;
    
    cout<<endl<<ktr+1<<"  >>>>>"<<treefilename[ktr].Data()<<"<<<<<"<<":"<<endl;
    cout<<"Read total cross section: "<<xsec<<" pb"<<endl;
  }
  nentries_tree = (Int_t)tree->GetEntries();





  // Start analysis_____________________________________________________________
  iwarning=0;

  // Loop over events stored inside the tree
  for(Int_t k_n=0;k_n<nentries_tree;k_n++) {

    tree->GetEntry(k_n);

    n_j=0;
    n_clep=0;
    n_nu=0;

    for(Int_t k=0;k<8;k++){

      //Make sure that the first two particles are ingoing (initial state)
      if(istup[0] != -1 || istup[1] != -1){
	cout<<"***ERROR: the first 2 particles are NOT ingoing!!!"<<endl;
	cout<<"          Please check the event file."<<endl<<endl;
	cout<<"Execution stopped"<<endl;
	exit(0);
      }

      particle[k] = TLorentzVector(px[k],py[k],pz[k],E[k]);

      if(istup[k]==1){   //final-state particles
	if(TMath::Abs(idup[k]) <= 5 || idup[k]==21){   //jet
	  jet[n_j] = TLorentzVector(px[k],py[k],pz[k],E[k]);
          index_jet[n_j]=k;
	  n_j++;
	}
 	else if(TMath::Abs(idup[k])==11 || TMath::Abs(idup[k])==13 || 
		TMath::Abs(idup[k])==15){   //charged lepton
	  clep[n_clep] = TLorentzVector(px[k],py[k],pz[k],E[k]);
	  n_clep++;
	}
	else if(TMath::Abs(idup[k])==12 || TMath::Abs(idup[k])==14 || 
		TMath::Abs(idup[k])==16){   //neutrino
	  nu[n_nu] = TLorentzVector(px[k],py[k],pz[k],E[k]);
	  n_nu++;
	}
	else{
	  cout<<"***ERROR: particle "<<idup[k]
	      <<" is not a final-state fermion or gluon!"<<endl;
	  cout<<"          Please check the event file."<<endl<<endl;
	  cout<<"Execution stopped"<<endl;
	  exit(0);
	}
      }
    }


    //Identify TAG JETS (jf,jb) and complementary ones (jc)
    if(itag_scheme==1){

      Int_t index_sort[n_j];
      Double_t eta[n_j];
      for(Int_t k=0;k<n_j;k++){
	eta[k]=jet[k].Eta();
      }
      TMath::Sort(n_j,eta,index_sort,kFALSE); // kFALSE => increasing order
      jf = jet[index_sort[n_j-1]];
      jb = jet[index_sort[0]];
      if((n_j % 2) != 0){
	cout<<"***ERROR: odd number of jets!!!"<<endl;
	exit(0);
      }
      Int_t nc = (Int_t) n_j/2;
      jc[0]=jet[index_sort[nc-1]];
      jc[1]=jet[index_sort[nc]];
    }

    else if(itag_scheme==2){
      Int_t index_sort[n_j];
      Double_t pT[n_j];
      for(Int_t k=0;k<n_j;k++){
	pT[k]=jet[k].Pt();
      }
      TMath::Sort(n_j,pT,index_sort,kFALSE); // kFALSE => increasing order
      jf = jet[index_sort[n_j-1]];
      jb = jet[index_sort[n_j-2]];
      if((n_j % 2) != 0){
	cout<<"***ERROR: odd number of jets!!!"<<endl;
	exit(0);
      }
      jc[0]=jet[index_sort[0]];
      jc[1]=jet[index_sort[1]];
    }

    else{
      cout<<"ERROR"<<endl;
      exit(0);
    }

    //Identify heavy bosons: by convention V[0] decays hadronically, V[1] 
    //leptonically.
    V[0]=jc[0]+jc[1];
    V[1]=clep[0]+nu[0];
    VV=V[0]+V[1];



    // Start selection procedure_______________________________________________

    ib=0;
    ibtag=0;
    ipasscuts=1;

    //Cut on energy, pT and pseudorapidity of jets
    for(Int_t k=0;k<n_j;k++){
      if(setcuts_E_min_j=="y"){
	if(jet[k].E() < cut_E_min_j){
	  ipasscuts=0;
	}
      }
      if(setcuts_pT_min_j=="y"){
	if(jet[k].Pt() < cut_pT_min_j){
	  ipasscuts=0;

	}
      }
      if(setcuts_eta_max_j=="y"){
	if(TMath::Abs(jet[k].Eta()) > cut_eta_max_j){
	  ipasscuts=0;
	}
      }
    }
    //Cut on energy, pT and pseudorapidity of charged leptons
    for(Int_t k=0;k<n_clep;k++){
      if(setcuts_E_min_clep=="y"){
	if(clep[k].E() < cut_E_min_clep){
	  ipasscuts=0;
	}
      }
      if(setcuts_pT_min_clep=="y"){
	if(clep[k].Pt() < cut_pT_min_clep){
	  ipasscuts=0;
	}
      }
      if(setcuts_eta_max_clep=="y"){
	if(TMath::Abs(clep[k].Eta()) > cut_eta_max_clep){
	  ipasscuts=0;
	}
      }
    }


    //Mimic b-tagging
    if(setBveto=="y" && (setBsingletag=="y" || setBdoubletag=="y")){
      cout<<"***Conflict in user-defined cuts. Check b-tagging setup"<<endl;
      exit(0);
    }

    //b-quark veto
    if(setBveto=="y"){
      for(Int_t k=2;k<8;k++) {      // for OUTGOING particles only
	if(TMath::Abs(idup[k])==5) {
	  if( (TMath::Abs(particle[k].Eta()) < Btagging_region)
	      && (irand.Uniform() > (1.0-Btagging_efficiency)) ) {
	    ipasscuts=0;
	  }
	}
        else if(TMath::Abs(idup[k])==4){
	  if( (TMath::Abs(particle[k].Eta()) < Btagging_region)
	      && (irand.Uniform() > (1.0-Cfake_efficiency)) ) {
	    ipasscuts=0;
	  }
	}
        else if(TMath::Abs(idup[k])<4 || TMath::Abs(idup[k])==21){
	  if( (TMath::Abs(particle[k].Eta()) < Btagging_region)
	      && (irand.Uniform() > (1.0-Qfake_efficiency)) ) {
	    ipasscuts=0;
	  }
	}
        
      }
    }

    //b-quark tagging
    if(setBsingletag=="y" || setBdoubletag=="y"){
      for(Int_t k=2;k<8;k++) {  //OUTGOING particles only
	if(TMath::Abs(idup[k])==5) {
	  if( (TMath::Abs(particle[k].Eta()) < Btagging_region) 
	      && (irand.Uniform() > (1.0-Btagging_efficiency)) ) {
	    ibtag++;
            is_a_b=k;   // useful for single tag
            index_b[ibtag-1]=k;
            if(ibtag==1){
              p_bb=particle[k];}
            else{
              p_bb+=particle[k];}
	  }
	}
        else if(TMath::Abs(idup[k])==4){
	  if( (TMath::Abs(particle[k].Eta()) < Btagging_region)
	      && (irand.Uniform() > (1.0-Cfake_efficiency)) ) {
	    ibtag++;
            is_a_b=k;   // useful for single tag
            index_b[ibtag-1]=k;
            if(ibtag==1){
              p_bb=particle[k];}
            else{
              p_bb+=particle[k];}
	  }
	}
        else if(TMath::Abs(idup[k])<4 || TMath::Abs(idup[k])==21){
	  if( (TMath::Abs(particle[k].Eta()) < Btagging_region)
	      && (irand.Uniform() > (1.0-Qfake_efficiency)) ) {
	    ibtag++;
            is_a_b=k;   // useful for single tag
            index_b[ibtag-1]=k;
            if(ibtag==1){
              p_bb=particle[k];}
            else{
              p_bb+=particle[k];}
	  }
	}
        
      }
      if(ibtag==0){
	ipasscuts=0;
      }
    }

 // Minimum DeltaR between bb system and one other jet
    DR=1.e4;
    for(Int_t k=0;k<n_j;k++) {
      if(index_jet[k]!= index_b[0]  && index_jet[k]!=index_b[1]){ 
        auxa=p_bb.Eta()-particle[index_jet[k]].Eta();
        auxb=p_bb.DeltaPhi(particle[index_jet[k]]);
        auxdelta = TMath::Sqrt(auxa*auxa+auxb*auxb);
        if(auxdelta < DR){
           DR = auxdelta;
        }
      }
    }

  if(setcuts_DRbbj=="y"){
      if(DR < cut_DRbbj){
	ipasscuts=0;
      }
  }
    
    //Cut on the invariant mass of any two jets
    if(ipasscuts==1){
      if(setcuts_Mjj=="y"){
	for(Int_t ki=0;ki<n_j;ki++){
	  for(Int_t kj=ki+1;kj<n_j;kj++){
	    jj = jet[ki]+jet[kj];
	    if(jj.M() < cut_Mjj) {
	      ipasscuts=0;
	    }
	  }
	}
      }
    }

    if(ipasscuts==1) {
      //RECONSTRUCT LONGITUDINAL MOMENTUM (AND THUS ENERGY) OF NEUTRINO
      //Notice: in general two solutions exist. One must take the both into
      //account (pL_nu1st, pl_nu2nd). Moreover this strategy does NOT allow
      //to cut on the invariant mass of the leptonic W!

      Mw2=80.40*80.40;  // (W mass)^2

      auxb = Mw2*clep[0].Pz() + 2.0*clep[0].Px()*nu[0].Px()*clep[0].Pz() 
	+ 2.0*clep[0].Py()*nu[0].Py()*clep[0].Pz();

      auxdelta = (clep[0].E()*clep[0].E())*
	( Mw2*Mw2                                              +
	  4.0*Mw2*clep[0].Px()*nu[0].Px()                      -
	  4.0*clep[0].E()*clep[0].E()*nu[0].Px()*nu[0].Px()    +
	  4.0*clep[0].Px()*clep[0].Px()*nu[0].Px()*nu[0].Px()  +
	  4.0*Mw2*clep[0].Py()*nu[0].Py()                      +
	  8.0*clep[0].Px()*nu[0].Px()*clep[0].Py()*nu[0].Py()  -
	  4.0*clep[0].E()*clep[0].E()*nu[0].Py()*nu[0].Py()    +
	  4.0*clep[0].Py()*clep[0].Py()*nu[0].Py()*nu[0].Py()  +
	  4.0*nu[0].Px()*nu[0].Px()*clep[0].Pz()*clep[0].Pz()  +
	  4.0*nu[0].Py()*nu[0].Py()*clep[0].Pz()*clep[0].Pz() );
      auxa = 2.0*(clep[0].E()*clep[0].E()-clep[0].Pz()*clep[0].Pz());

      pL_nu1st = (auxb - TMath::Sqrt(auxdelta))/(auxa);
      E_nu1st = TMath::Sqrt(nu[0].Px()*nu[0].Px() + nu[0].Py()*nu[0].Py() + 
			    pL_nu1st*pL_nu1st);

      pL_nu2nd = (auxb + TMath::Sqrt(auxdelta))/(auxa);
      E_nu2nd = TMath::Sqrt(nu[0].Px()*nu[0].Px() + nu[0].Py()*nu[0].Py() + 
			    pL_nu2nd*pL_nu2nd);

      if(auxdelta < 0.0) {  //no solutions
	//cout<<"WARNING: NO SOLUTIONS for neutrino longitudinal momentum!!!"
	//    <<endl;
	//cout<<"auxdelta = "<<auxdelta<<endl;
	//cout<<"k_n = "<<k_n<<endl;
	//cout<<"pL_nu1st = "<<pL_nu1st<<endl;
	//cout<<"pL_nu2nd = "<<pL_nu2nd<<endl;
	ipasscuts=0;
	iwarning++;
      }


      //pass cuts on INVARIANT MASS OF THE TWO CENTRAL JETS
      if(setcuts_Mjcjc=="y"){
	if(TMath::Abs(V[0].M()-centralvalue_Mjcjc) > cut_Mjcjc){
	  ipasscuts=0;
	}
	//reject if also the 2 tag jets fall into the invariant mass window
	jfjb=jf+jb;
	if(TMath::Abs(jfjb.M()-centralvalue_Mjcjc) < cut_Mjcjc){
	  ipasscuts=0;
	}
      }
    } //end if(ipasscuts==1) {



    if(ipasscuts==1) {
      //pass cuts on |Delta eta| between FORWARD/BACKWARD JETS
      if(setcuts_deltaetajfjb=="y"){
	if(TMath::Abs(jf.Eta()-jb.Eta()) < cut_deltaetajfjb){
	  ipasscuts=0;
	}
      }
    }

    if(ipasscuts==1) {

      //RECONSTRUCT leptonic W (2 general solutions)
      Vlep_reco[0] = TLorentzVector(V[1].Px(),V[1].Py(),clep[0].Pz()+pL_nu1st,
				    clep[0].E()+E_nu1st);
      Vlep_reco[1] = TLorentzVector(V[1].Px(),V[1].Py(),clep[0].Pz()+pL_nu2nd,
				    clep[0].E()+E_nu2nd);
      //choose minimum pseudorapidity
      Double_t eta_Wlv_rec = TMath::Min(Vlep_reco[0].Eta(),Vlep_reco[1].Eta());


      //RECONSTRUCT invariant masses of {leptonic W + jf/jb} and find minimum
      //Notice: results do NOT coincide with the "MC truth" because the
      //RECONSTRUCTED longitudinal momentum of the neutrino is used. As in 
      //general two solutions exist for pL_nu, the both are considered in 
      //finding the minimum.
      Vlep_reco_jf[0] = Vlep_reco[0]+jf;
      Vlep_reco_jf[1] = Vlep_reco[1]+jf;
      Vlep_reco_jb[0] = Vlep_reco[0]+jb;
      Vlep_reco_jb[1] = Vlep_reco[1]+jb;
      Double_t M_lvjtag_rec[4];
      Double_t M_lvj_min_rec;
      M_lvjtag_rec[0] = Vlep_reco_jf[0].M();
      M_lvjtag_rec[1] = Vlep_reco_jf[1].M();
      M_lvjtag_rec[2] = Vlep_reco_jb[0].M();
      M_lvjtag_rec[3] = Vlep_reco_jb[1].M();
      Int_t index_sort[4];
      TMath::Sort(4,M_lvjtag_rec,index_sort,kFALSE); //increasing order
      M_lvj_min_rec = M_lvjtag_rec[index_sort[0]];

      if(setcuts_Wlv=="y"){
	//further cuts on W reconstructed from leptons
	if(V[1].Pt() < cut_pT_Wlv || TMath::Abs(eta_Wlv_rec) > cut_eta_Wlv_rec
	   || M_lvj_min_rec < cut_M_lvj_min_rec){
	  ipasscuts=0;
	}
      }

    }  //end if(ipasscuts==1) {


    if(ipasscuts==1) {
      if(setTopcut_jjj=="y"){
	//cut on invariant mass of jet triplets {jjj} for top subtraction
	for(Int_t ki=0;ki<n_j;ki++){
	  for(Int_t kj=ki+1;kj<n_j;kj++){
	    for(Int_t kk=kj+1;kk<n_j;kk++){
	      jjj=jet[ki]+jet[kj]+jet[kk];
	      if(TMath::Abs(jjj.M()-centralvalue_Mjjj) < cut_Mjjj){
		ipasscuts=0;
	      }
	    }
	  }
	}
      }
    } // end if(ipasscuts==1) {


    if(ipasscuts==1) {
      if(setTopcut_jlv_rec=="y"){
	//cut on invariant mass of {jlv} triplets (2 solutions for each j) for
	//top subtraction
	for(Int_t ki=0;ki<n_j;ki++){
	  for(Int_t kj=0;kj<2;kj++){
	    jlv=jet[ki]+Vlep_reco[kj];
	    if(TMath::Abs(jlv.M()-centralvalue_Mjlv_rec) < cut_Mjlv_rec){
	      ipasscuts=0;
	    }
	  }
	}
      }
    } //end if(ipasscuts==1) {



    // Fill histograms

    hew_Mlv_allproc_buf->Fill(V[1].M());
    hew_Mjj_allproc_buf->Fill(V[0].M());
    hew_Mjjlv_allproc_buf->Fill(VV.M());
//Fill Mbb histogram before cuts
    if(ibtag==2 && setBdoubletag=="y"){
          hew_Mbb_allproc_buf->Fill(p_bb.M());
    }
    
    Double_t DeltaEta_jfjb = TMath::Abs(jf.Eta()-jb.Eta());
    if(TMath::Abs(V[0].M()-rmh)<5.0){
      hew_pTlv_peak_buf->Fill(V[1].Pt());
      hew_DEta_jfjb_peak_buf->Fill(DeltaEta_jfjb);
      if(ibtag==2 && setBdoubletag=="y"){
        auxa=particle[index_b[0]].Eta()-particle[index_b[1]].Eta();
        auxb=particle[index_b[0]].DeltaPhi(particle[index_b[1]]);
        DR = TMath::Sqrt(auxa*auxa+auxb*auxb);
        hew_DR_bb_peak_buf->Fill(DR);
      }
    }

    if(ipasscuts==1) {
      hew_Mlv_cuts_buf->Fill(V[1].M());
      hew_Mjj_cuts_buf->Fill(V[0].M());
      hew_Mjjlv_cuts_buf->Fill(VV.M());

      if(ibtag==2 && setBdoubletag=="y"){
        hew_Mbb_cuts_buf->Fill(p_bb.M());
        if(TMath::Abs(p_bb.M()-rmh)<5.0){
          hew_pTbb_peak_buf->Fill(p_bb.Pt());
 // Minimum DeltaR between bb system and one other jet
          hew_DRmin_bbj_peak_buf->Fill(DR);
        }
      }

      if(ibtag==1 && setBsingletag=="y") {
        Int_t ind=0;
        for(Int_t k=0;k<n_j;k++) {  //OUTGOING jets only
          if(index_jet[k]!=is_a_b) {
            p_bj = jet[k]+particle[is_a_b];
            M_bj[ind] = p_bj.M();
            pt_bj[ind] = p_bj.Pt();
            ind++;
          }
        }

//Order masses
        Int_t index_sort[ind];
        TMath::Sort(ind,M_bj,index_sort,kFALSE); // kFALSE => increasing order
        hew_Mbj1_cuts_buf->Fill(M_bj[index_sort[0]]);
        hew_Mbj2_cuts_buf->Fill(M_bj[index_sort[1]]);
        hew_Mbj3_cuts_buf->Fill(M_bj[index_sort[2]]);
        if(TMath::Abs(M_bj[index_sort[0]]-rmh)<5.0){
          hew_pTbj_peak_buf->Fill(pt_bj[index_sort[0]]);
        } 
      }
    }


    if(ipasscuts==1){

      //FILL OTHER HISTOGRAMS___________________________________________________
      if(VV.M() > Mcut){

	//1. minimum pT of tag jets
	Double_t pTmin_jfjb = TMath::Min(jf.Pt(),jb.Pt());
	hew_pTmin_jfjb_buf->Fill(pTmin_jfjb);

	//2. minimum pT of all jets
	Int_t index_sort[n_j];
	Double_t pT_j[n_j];
	for(Int_t k=0;k<n_j;k++){
	  pT_j[k]=jet[k].Pt();
	}
	TMath::Sort(n_j,pT_j,index_sort,kFALSE); // kFALSE => increasing order
	hew_pTmin_alljets_buf->Fill(pT_j[index_sort[0]]);

	//3. pT of heavy boson decaying hadronically (=> two central jets)
	hew_pT_jcjc_buf->Fill(V[0].Pt());

       	//4. eta of the heavy boson decaying hadronically (=> two central jets)
	hew_eta_jcjc_buf->Fill(V[0].Eta());

      	//5. invariant mass of the two tag jets
	hew_M_jfjb_buf->Fill(jfjb.M());

	//6. min. inv. mass of the heavy boson decaying hadronically and one tag
	//   jet
	Vhad_jf = V[0]+jf;
	Vhad_jb = V[0]+jb;
	Double_t M_jcjcjfb_min = TMath::Min(Vhad_jf.M(),Vhad_jb.M());
	hew_Mmin_jcjcjfb_buf->Fill(M_jcjcjfb_min);

	//7. |Delta eta| between tag jets
	hew_DeltaEta_jfjb_buf->Fill(DeltaEta_jfjb);

	//8. |Delta eta| between reconstructed heavy bosons
	Double_t DeltaEta_VV = TMath::Abs(V[0].Eta()-V[1].Eta());
	hew_DeltaEta_VV_buf->Fill(DeltaEta_VV);

	//9. maximum pT of tag jets
	Double_t pTmax_jfjb = TMath::Max(jf.Pt(),jb.Pt());
	hew_MaxPT_tagj_buf->Fill(pTmax_jfjb);

	//10. minimum inv. mass of tag jet + vector boson
	//10bis. complementary inv. mass of tag jet + vector boson
	//[i.e.: data la minima massa tag-jet/bosone, la massa dell'altro 
	//tag-jet con l'altro bosone]
	Double_t MinM_Vj=0.0;
	Vlep_jf = V[1]+jf;
	Vlep_jb = V[1]+jb;
	MinM_Vj = TMath::Min(Vhad_jf.M(),Vhad_jb.M());
	MinM_Vj = TMath::Min(MinM_Vj,Vlep_jf.M());
	MinM_Vj = TMath::Min(MinM_Vj,Vlep_jb.M());
//	if(MinM_Vj == Vhad_jf.M()){
	if(TMath::Abs(MinM_Vj - Vhad_jf.M()) < 1.0E-10){
	  hew_MinM_Vj_buf->Fill(Vhad_jf.M());
	  hew_comp_MinM_Vj_buf->Fill(Vlep_jb.M());
	}
//	else if(MinM_Vj == Vhad_jb.M()){
	else if(TMath::Abs(MinM_Vj - Vhad_jb.M()) < 1.0E-10){
	  hew_MinM_Vj_buf->Fill(Vhad_jb.M());
	  hew_comp_MinM_Vj_buf->Fill(Vlep_jf.M());
	}
//	else if(MinM_Vj == Vlep_jf.M()){
	else if(TMath::Abs(MinM_Vj - Vlep_jf.M()) < 1.0E-10){
	  hew_MinM_Vj_buf->Fill(Vlep_jf.M());
	  hew_comp_MinM_Vj_buf->Fill(Vhad_jb.M());
	}
//	else if(MinM_Vj == Vlep_jb.M()){
	else if(TMath::Abs(MinM_Vj - Vlep_jb.M()) < 1.0E-10){
	  hew_MinM_Vj_buf->Fill(Vlep_jb.M());
	  hew_comp_MinM_Vj_buf->Fill(Vhad_jf.M());
	}
	else{
	  cout<<"***ERROR !"<<endl;
	  exit(0);
	}


	//11. minimum delta-eta between tag jet and the vector boson closest 
	//    to it
	//11bis. complementary delta-eta between the tag jet and closest 
	//       vector boson
	Double_t MinDEta_Vtag=0.0;
	Double_t DEta_jfjcjc = TMath::Abs(jf.Eta() - V[0].Eta());
	Double_t DEta_jbjcjc = TMath::Abs(jb.Eta() - V[0].Eta());
	Double_t DEta_jflv = TMath::Abs(jf.Eta() - V[1].Eta());
	Double_t DEta_jblv = TMath::Abs(jb.Eta() - V[1].Eta());
	MinDEta_Vtag=TMath::Min(DEta_jfjcjc,DEta_jbjcjc);
	MinDEta_Vtag=TMath::Min(MinDEta_Vtag,DEta_jflv);
	MinDEta_Vtag=TMath::Min(MinDEta_Vtag,DEta_jblv);
//	if(MinDEta_Vtag == DEta_jfjcjc){
	if(TMath::Abs(MinDEta_Vtag - DEta_jfjcjc) < 1.0E-10){	
	  hew_MinDEta_Vtag_buf->Fill(DEta_jfjcjc);
	  hew_comp_MinDEta_Vtag_buf->Fill(DEta_jblv);
	}
//	else if(MinDEta_Vtag == DEta_jbjcjc){
	else if(TMath::Abs(MinDEta_Vtag - DEta_jbjcjc) < 1.0E-10){
	  hew_MinDEta_Vtag_buf->Fill(DEta_jbjcjc);
	  hew_comp_MinDEta_Vtag_buf->Fill(DEta_jflv);
	}
//	else if(MinDEta_Vtag == DEta_jflv){
	else if(TMath::Abs(MinDEta_Vtag - DEta_jflv) < 1.0E-10){
	  hew_MinDEta_Vtag_buf->Fill(DEta_jflv);
	  hew_comp_MinDEta_Vtag_buf->Fill(DEta_jbjcjc);
	}
//	else if(MinDEta_Vtag == DEta_jblv){
	else if(TMath::Abs(MinDEta_Vtag - DEta_jblv) < 1.0E-10){
	  hew_MinDEta_Vtag_buf->Fill(DEta_jblv);
	  hew_comp_MinDEta_Vtag_buf->Fill(DEta_jfjcjc);
	}
	else{
	  cout<<"***ERROR !!";
	  exit(0);
	}

	//12. |Delta phi| between tag jets
	Double_t DPhi_jfjb = TMath::Abs(jf.DeltaPhi(jb));
	hew_DPhi_tag_buf->Fill(DPhi_jfjb);

	//13. |Delta phi| between reconstructed vector bosons
	Double_t DPhi_VV = TMath::Abs(V[0].DeltaPhi(V[1]));
	hew_DPhi_V_buf->Fill(DPhi_VV);
      }
    }


  } //end for(Int_t k_n=0;k_n<nentries_tree;k_n++)


  cout<<"WARNING: discarded "<<iwarning<<" events with no real solutions for reconstructed pz of neutrino in file "<<treefilename[ktr].Data()<<endl;

  treefile.Close();


  // Normalize histograms to the total cross section____________________________
  if(opt=="xs") {
    // Scale factors for histograms; having been multiplied by them, 
    // histograms show the DIFFERENTIAL CROSS SECTION.
    // [scale = sigma_tot/(A*nentries)]
    A_Mlv = (Mlv_high - Mlv_low)/nbins_M;
    scale_Mlv = xsec/(A_Mlv*(nentries_tree));
    A_Mjj = (Mjj_high - Mjj_low)/nbins_M;
    scale_Mjj = xsec/(A_Mjj*(nentries_tree));
    A_Mjjlv = (Mjjlv_high - Mjjlv_low)/nbins_M;
    scale_Mjjlv = xsec/(A_Mjjlv*(nentries_tree));
    
    A_eta = (eta_high - eta_low)/nbins_eta;
    scale_eta = xsec/(A_eta*(nentries_tree));
    A_pT = (pT_high - pT_low)/nbins_pT;
    scale_pT = xsec/(A_pT*(nentries_tree));
    A_phi = (phi_high - phi_low)/nbins_phi;
    scale_phi = xsec/(A_phi*(nentries_tree));

    hew_Mlv_allproc_buf->Scale(scale_Mlv);
    hew_Mlv_cuts_buf->Scale(scale_Mlv);
    hew_Mjj_allproc_buf->Scale(scale_Mjj);
    hew_Mjj_cuts_buf->Scale(scale_Mjj);
    hew_Mbb_allproc_buf->Scale(scale_Mjj);
    hew_Mbb_cuts_buf->Scale(scale_Mjj);
    hew_Mbj1_cuts_buf->Scale(scale_Mjj);
    hew_Mbj2_cuts_buf->Scale(scale_Mjj);
    hew_Mbj3_cuts_buf->Scale(scale_Mjj);
    hew_pTbb_peak_buf->Scale(scale_pT);
    hew_pTbj_peak_buf->Scale(scale_pT);
    hew_pTlv_peak_buf->Scale(scale_pT);
    hew_DEta_jfjb_peak_buf->Scale(scale_eta);
    hew_DR_bb_peak_buf->Scale(scale_eta);
    hew_DRmin_bbj_peak_buf->Scale(scale_eta);
    hew_Mjjlv_allproc_buf->Scale(scale_Mjjlv);
    hew_Mjjlv_cuts_buf->Scale(scale_Mjjlv);
    hew_pTmin_jfjb_buf->Scale(scale_pT);
    hew_pTmin_alljets_buf->Scale(scale_pT);
    hew_pT_jcjc_buf->Scale(scale_pT);
    hew_eta_jcjc_buf->Scale(scale_eta);
    hew_M_jfjb_buf->Scale(scale_Mjjlv);
    hew_Mmin_jcjcjfb_buf->Scale(scale_Mjjlv);
    hew_DeltaEta_jfjb_buf->Scale(xsec/(((2*eta_high)/nbins_eta)*(nentries_tree)));
    hew_DeltaEta_VV_buf->Scale(xsec/(((2*eta_high)/nbins_eta)*(nentries_tree)));
    hew_MaxPT_tagj_buf->Scale(scale_pT);
    hew_MinM_Vj_buf->Scale(scale_Mjjlv);
    hew_comp_MinM_Vj_buf->Scale(scale_Mjjlv);
    hew_MinDEta_Vtag_buf->Scale(scale_eta);
    hew_comp_MinDEta_Vtag_buf->Scale(scale_eta);
    hew_DPhi_tag_buf->Scale(scale_phi);
    hew_DPhi_V_buf->Scale(scale_phi);
  }


  //Unbuffer temporary histograms (clones)
  hew_Mlv_allproc.Add(hew_Mlv_allproc_buf);
  hew_Mlv_cuts.Add(hew_Mlv_cuts_buf);
  hew_Mjj_allproc.Add(hew_Mjj_allproc_buf);
  hew_Mjj_cuts.Add(hew_Mjj_cuts_buf);
  hew_Mbb_allproc.Add(hew_Mbb_allproc_buf);
  hew_Mbb_cuts.Add(hew_Mbb_cuts_buf);
  hew_Mbj1_cuts.Add(hew_Mbj1_cuts_buf);
  hew_Mbj2_cuts.Add(hew_Mbj2_cuts_buf);
  hew_Mbj3_cuts.Add(hew_Mbj3_cuts_buf);
  hew_pTbb_peak.Add(hew_pTbb_peak_buf);
  hew_pTbj_peak.Add(hew_pTbj_peak_buf);
  hew_pTlv_peak.Add(hew_pTlv_peak_buf);
  hew_DEta_jfjb_peak.Add(hew_DEta_jfjb_peak_buf);
  hew_DR_bb_peak.Add(hew_DR_bb_peak_buf);
  hew_DRmin_bbj_peak.Add(hew_DRmin_bbj_peak_buf);
  hew_Mjjlv_allproc.Add(hew_Mjjlv_allproc_buf);
  hew_Mjjlv_cuts.Add(hew_Mjjlv_cuts_buf);
  hew_pTmin_jfjb.Add(hew_pTmin_jfjb_buf);
  hew_pTmin_alljets.Add(hew_pTmin_alljets_buf);
  hew_pT_jcjc.Add(hew_pT_jcjc_buf);
  hew_eta_jcjc.Add(hew_eta_jcjc_buf);
  hew_M_jfjb.Add(hew_M_jfjb_buf);
  hew_Mmin_jcjcjfb.Add(hew_Mmin_jcjcjfb_buf);
  hew_DeltaEta_jfjb.Add(hew_DeltaEta_jfjb_buf);
  hew_DeltaEta_VV.Add(hew_DeltaEta_VV_buf);
  hew_MaxPT_tagj.Add(hew_MaxPT_tagj_buf);
  hew_MinM_Vj.Add(hew_MinM_Vj_buf);
  hew_comp_MinM_Vj.Add(hew_comp_MinM_Vj_buf);
  hew_MinDEta_Vtag.Add(hew_MinDEta_Vtag_buf);
  hew_comp_MinDEta_Vtag.Add(hew_comp_MinDEta_Vtag_buf);
  hew_DPhi_tag.Add(hew_DPhi_tag_buf);
  hew_DPhi_V.Add(hew_DPhi_V_buf);

  hew_Mlv_allproc_buf->Reset();
  hew_Mlv_cuts_buf->Reset();
  hew_Mjj_allproc_buf->Reset();
  hew_Mjj_cuts_buf->Reset();
  hew_Mbb_allproc_buf->Reset();
  hew_Mbb_cuts_buf->Reset();
  hew_Mbj1_cuts_buf->Reset();
  hew_Mbj2_cuts_buf->Reset();
  hew_Mbj3_cuts_buf->Reset();
  hew_pTbb_peak_buf->Reset();
  hew_pTbj_peak_buf->Reset();
  hew_pTlv_peak_buf->Reset();
  hew_DEta_jfjb_peak_buf->Reset();
  hew_DR_bb_peak_buf->Reset();
  hew_DRmin_bbj_peak_buf->Reset();
  hew_Mjjlv_allproc_buf->Reset();
  hew_Mjjlv_cuts_buf->Reset();
  hew_pTmin_jfjb_buf->Reset();
  hew_pTmin_alljets_buf->Reset();
  hew_pT_jcjc_buf->Reset();
  hew_eta_jcjc_buf->Reset();
  hew_M_jfjb_buf->Reset();
  hew_Mmin_jcjcjfb_buf->Reset();
  hew_DeltaEta_jfjb_buf->Reset();
  hew_DeltaEta_VV_buf->Reset();
  hew_MaxPT_tagj_buf->Reset();
  hew_MinM_Vj_buf->Reset();
  hew_comp_MinM_Vj_buf->Reset();
  hew_MinDEta_Vtag_buf->Reset();
  hew_comp_MinDEta_Vtag_buf->Reset();
  hew_DPhi_tag_buf->Reset();
  hew_DPhi_V_buf->Reset();

  } //end for(Int_t ktr=0;ktr<ntr;ktr++){


  cout<<endl<<"Total cross section: "<<xsec_tot<<" pb"<<endl;

  


  //Information about the integrated cross sections will be stored in a canvas
  TPaveText text_xsec(0.1,0.1,0.9,0.9);
  text_xsec.SetLabel("Integrated cross sections");
  TH1D *h_ref = new TH1D(hew_Mjjlv_cuts);
  TString s;
  s="total cross section before cuts";
  s+="    ";
  s+=xsec_tot;
  s+=" pb";
  text_xsec.AddText(s.Data());
  s="total cross section after cuts";
  s+="    ";
  s+=integro(h_ref);
  s+=" pb";
  text_xsec.AddText(s.Data());
  text_xsec.AddText(" ");
  text_xsec.AddText(" M_{cut} (GeV)                       #sigma (pb)");
  for(Int_t kl=0;kl<8;kl++){
    Double_t mc=600.0+200*kl;
    s="";
    s+=mc;
    s+="                   ";
    s+=integro(h_ref,mc,14000.0);
    s+=" pb";
    text_xsec.AddText(s.Data());
  }
  text_xsec.AddText(" ");




  // Print results on file plots.root________________________________

  TFile *fo = new TFile("plots.root","update");

  //  TCanvas *c_cuts = new TCanvas("list_of_selection_cuts");
  TCanvas c_cuts("list_of_selection_cuts");
  text_cuts.Draw();
  c_cuts.Write();

  //  TCanvas *c_xsec = new TCanvas("total_cross_section");
  TCanvas c_xsec("total_cross_section");
  text_xsec.Draw();
  c_xsec.Write();

  fo->cd("MVV_allevents");
  hew_Mlv_allproc.Write();
  hew_Mjj_allproc.Write();
  hew_Mbb_allproc.Write();
  hew_Mjjlv_allproc.Write();
  hew_pTlv_peak.Write();
  hew_DEta_jfjb_peak.Write();
  hew_DR_bb_peak.Write();
  fo->cd();
  fo->cd("MVV_cuts");
  hew_Mlv_cuts.Write();
  hew_Mjj_cuts.Write();
  hew_Mbb_cuts.Write();
  hew_Mbj1_cuts.Write();
  hew_Mbj2_cuts.Write();
  hew_Mbj3_cuts.Write();
  hew_pTbb_peak.Write();
  hew_pTbj_peak.Write();
  hew_Mjjlv_cuts.Write();
  hew_DRmin_bbj_peak.Write();
  fo->cd();
  fo->cd("otherhistos_cuts");
  hew_pTmin_jfjb.Write();
  hew_pTmin_alljets.Write();
  hew_pT_jcjc.Write();
  hew_eta_jcjc.Write();
  hew_M_jfjb.Write();
  hew_Mmin_jcjcjfb.Write();
  hew_DeltaEta_jfjb.Write();
  hew_DeltaEta_VV.Write();
  hew_MaxPT_tagj.Write();
  hew_MinM_Vj.Write();
  hew_comp_MinM_Vj.Write();
  hew_MinDEta_Vtag.Write();
  hew_comp_MinDEta_Vtag.Write();
  hew_DPhi_tag.Write();
  hew_DPhi_V.Write();
  fo->cd();

  //Save file on disk
  fo->Write();

  cout<<endl;
  cout<<"created file plots.root"<<endl;
  cout<<endl;
  cout<<"Execution successful"<<endl;
  cout<<endl;

}








// Adapted from integro.c (by S.Bolognesi)

double integro(TH1D *histo)
{
  int nbins;
  double width,integral;
  width = histo->GetBinWidth(1);
  //integral = (histo->Integral())*width;
  nbins=histo->GetNbinsX();
  integral = (histo->Integral(1,nbins,"width"));
  return integral;
}

double integro(TH1D *histo,double low,double high)
{
  double width,min,max,integral;
  int highBin,lowBin;
  width = histo->GetBinWidth(1);
  
  TAxis *xaxis = histo->GetXaxis();
  min=xaxis->GetXmin();
  max=xaxis->GetXmax();
  
  if(low>=min) lowBin=(int)TMath::Ceil((low-min)/width);
  else
    {
      cout<<"*** integro.c ERROR: lower integ. limit .LT.< histogram range"<<endl;
//      return 0;
      exit(0);
    }
  
  if(high<=max) highBin=(int)TMath::Ceil((high-min)/width);
  else 
    {
//      cout<<"limite sup. > del range"<<endl;
      highBin=(int)TMath::Ceil((max-min)/width);
    }

  low = xaxis->GetBinLowEdge(lowBin);
  high = ((xaxis->GetBinLowEdge(highBin))+width);

//  cout<<"integrale"<<endl<<"da "<<low<<" (bin "<<lowBin<<")"<<endl
//      <<"a "<<high<<" (bin "<<highBin<<")"<<endl; 
  integral= (histo->Integral(lowBin,highBin))*width;
  return integral;
}
