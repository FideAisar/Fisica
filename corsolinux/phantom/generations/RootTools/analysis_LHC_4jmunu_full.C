/*
 *analysis_LHC_4jmunu_full.C
 *
 *Last update: April 12, 2008
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
#include "TRegexp.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TList.h"
#include <TAxis.h>
#endif
using namespace std;


void analysis_LHC_4jmunu_full(void);
double integro(TH1D*);
double integro(TH1D*,double,double);


void analysis_LHC_4jmunu_full() {


  gSystem->Load("libPhysics.so");
  gROOT->Reset();


  //total cross section for histogram normalization
  Double_t xsec = 0.0;
  Double_t xsec_tot = 0.0;

  //parameters for histogram rescaling
  Double_t scale_Mlv, scale_Mjj, scale_Mjjlv, scale_eta, scale_pT, scale_phi;
  static Double_t A_Mlv, A_Mjj, A_Mjjlv, A_eta, A_pT, A_phi;

  //event kinematics
  Double_t px[8], py[8], pz[8], E[8];
  Int_t idup[8],istup[8];
  TLorentzVector particle[8], jet[8], clep[8], nu[8];
  TLorentzVector jf, jb, jc[2], jj, jjj, tot_vis, jlv, jfjb, V[2], Vlep_reco[2], 
    Vlep_reco_jf[2], Vlep_reco_jb[2], Vhad_jf, Vhad_jb, Vlep_jf, Vlep_jb, VV;
  Int_t n_j, n_clep, n_nu, index_jet[8];

  //parameters for histogram settings
  Double_t Mlv_high, Mlv_low;
  Double_t Mjj_high, Mjj_low;
  Double_t Mjjlv_high, Mjjlv_low;
  Double_t eta_high, eta_low;
  Double_t pT_high, pT_low; 
  Double_t phi_high, phi_low;
  Int_t nentries_tree=0;   //total number of entries of the selected TREE(S)
  Int_t nbins_M;           //number of bins for invariant mass histograms
  Int_t nbins_MVV;         //number of bins for VV invariant mass histograms
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

  Double_t DeltaEta_VV;
  
  //user-defined selection cuts
  TString setcuts_pT_min_j,setcuts_pT_min_clep,
    setcuts_eta_max_j, setcuts_eta_max_clep, 
    setcuts_Mjj, setcuts_Mjcjc, setcuts_deltaetajfjb, setcuts_Wlv, 
    setBveto, setBsingletag, setBdoubletag, setTopcut_jjj, setTopcut_jlv_rec,
    setcuts_deltaetaVV,setcuts_TotVisM,setcuts_minpT_miss,setcuts_minDeltaRjj,
    setcuts_minpT_jc1,setcuts_minpT_jc2,setcuts_eta_max_jc1,setcuts_eta_max_jc2,
    setcuts_Mjfjb, setcuts_minpT_jcjc,setcuts_deltaetaVj;
  Double_t cut_pT_min_j, cut_pT_min_clep, 
    cut_eta_max_j, cut_eta_max_clep, centralvalue_Mjcjc, cut_Mjcjc,
    cut_Mjj, cut_deltaetajfjb, cut_pT_Wlv, cut_eta_Wlv_rec, cut_M_lvj_min_rec, 
    Btagging_region, Btagging_efficiency, Cfake_efficiency, Qfake_efficiency, 
    centralvalue_Mjjj, cut_Mjjj, centralvalue_Mjlv_rec, cut_Mjlv_rec, Mcut,
    cut_deltaetaVV,cut_TotVisM,cut_minpT_miss,cut_minDeltaRjj,
    cut_minpT_jc1,cut_minpT_jc2,cut_eta_max_jc1,cut_eta_max_jc2,
    cut_Mjfjb, cut_minpT_jcjc,cut_deltaetaVj;
  //random number generator and the corresponding random variable
  TRandom3 irand;

  Double_t DR,DRmin,Mmin_jj,eta_jc1,eta_jc2,pT_jc1,pT_jc2,pT_jcjc;   

  static Int_t ipasscuts;

  Int_t iwarning=0;

  TString dummy; 
/*******************************************************************************
 * Beginning of USER-DEFINED SELECTION CUTS ("y"=cut applied; "n"=cut not 
 * applied)
 ******************************************************************************/

  //Minimum pT for jets
  cin>>setcuts_pT_min_j>>dummy;
  if(setcuts_pT_min_j=="y"){
    cin>>cut_pT_min_j>>dummy;
  }

  //Minimum pT for charged leptons
  cin>>setcuts_pT_min_clep>>dummy;
  if(setcuts_pT_min_clep=="y"){
    cin>>cut_pT_min_clep>>dummy;
  }

  //Maximum pseudorapidity for jets
  cin>>setcuts_eta_max_j>>dummy;
  if(setcuts_eta_max_j=="y"){
    cin>>cut_eta_max_j>>dummy;
  }

  //Maximum pseudorapidity for charged  leptons
  cin>>setcuts_eta_max_clep>>dummy;
  if(setcuts_eta_max_clep=="y"){
    cin>>cut_eta_max_clep>>dummy;
  }

  //Mimimum invariant mass of any two jets
  cin>>setcuts_Mjj>>dummy;
  if(setcuts_Mjj=="y"){
    cin>>cut_Mjj>>dummy;
  }

  //Cut on invariant mass of the two central jets
  cin>>setcuts_Mjcjc>>dummy;
  if(setcuts_Mjcjc=="y"){
    cin>>centralvalue_Mjcjc>>dummy; //GeV
    cin>>cut_Mjcjc>>dummy;  // => window = centralMjj +/- cutMjj
  }

  //Cut on the difference in pseudorapidity between forward and backward tag jet
  cin>>setcuts_deltaetajfjb>>dummy;
  if(setcuts_deltaetajfjb=="y"){
    cin>>cut_deltaetajfjb >>dummy;
  }

  //Cut on the difference in pseudorapidity between jcjc system and leptonic W
  cin>>setcuts_deltaetaVV>>dummy;
  if(setcuts_deltaetaVV=="y"){
    cin>>cut_deltaetaVV >>dummy;
  }

  //Cut on the total visible mass (4j+l)
  cin>>setcuts_TotVisM>>dummy;
  if(setcuts_TotVisM=="y"){
    cin>>cut_TotVisM >>dummy;
  }

  //Cuts on the leptonic W
  cin>>setcuts_Wlv>>dummy;
  if(setcuts_Wlv=="y"){
    cin>>cut_pT_Wlv>>dummy;
    cin>>cut_eta_Wlv_rec>>dummy;
    cin>>cut_M_lvj_min_rec>>dummy;
  }

  //Top cuts:
  //2.Cut on invariant mass of jjj system
  cin>>setTopcut_jjj>>dummy;
  if(setTopcut_jjj=="y"){
    cin>>centralvalue_Mjjj>>dummy; //GeV
    cin>>cut_Mjjj>>dummy; // => window = centralMjjj +/- cutMjjj
  }

  //2.Cut on invariant mass of jlv_rec
  cin>>setTopcut_jlv_rec>>dummy;
  if(setTopcut_jlv_rec=="y"){
    cin>>centralvalue_Mjlv_rec>>dummy; //GeV
    cin>>cut_Mjlv_rec>>dummy;  // => window = centralMjlv +/- cutMjlv
  }

  //Cut on minimum |Delta R| of any two jets
  cin>>setcuts_minDeltaRjj>>dummy;
  if(setcuts_minDeltaRjj=="y"){
    cin>>cut_minDeltaRjj>>dummy;
  }

  //Cuts on missing pT
  cin>>setcuts_minpT_miss>>dummy;
  if(setcuts_minpT_miss=="y"){
    cin>>cut_minpT_miss>>dummy;
  }

  //Cut on pT_jc1
  cin>>setcuts_minpT_jc1>>dummy;
  if(setcuts_minpT_jc1=="y"){
    cin>>cut_minpT_jc1>>dummy;
  }

  //Cut on pT_jc2
  cin>>setcuts_minpT_jc2>>dummy;
  if(setcuts_minpT_jc2=="y"){
    cin>>cut_minpT_jc2>>dummy;
  }

  //Maximum pseudorapidity for jc1
  cin>>setcuts_eta_max_jc1>>dummy;
  if(setcuts_eta_max_jc1=="y"){
    cin>>cut_eta_max_jc1>>dummy;
  }

  //Maximum pseudorapidity for jc2
  cin>>setcuts_eta_max_jc2>>dummy;
  if(setcuts_eta_max_jc2=="y"){
    cin>>cut_eta_max_jc2>>dummy;
  }

  //Cut on pT_jcjc
  cin>>setcuts_minpT_jcjc>>dummy;
  if(setcuts_minpT_jcjc=="y"){
    cin>>cut_minpT_jcjc>>dummy;
  }

  //Cut on invariant mass of the fb jets
  cin>>setcuts_Mjfjb>>dummy;
  if(setcuts_Mjfjb=="y"){
    cin>>cut_Mjfjb>>dummy;
  }

  //Cut on the difference in pseudorapidity vector bosons and 
  // forward/backward tag jet
  cin>>setcuts_deltaetaVj>>dummy;
  if(setcuts_deltaetaVj=="y"){
    cin>>cut_deltaetaVj >>dummy;
  }

  cin>>Mcut>>dummy;     //M_min_VVsytem

  //Mimic b-quark veto or b-tagging
  Btagging_region = 1.5;  // B-tagging active on |eta| < Btagging_region
//  Btagging_efficiency=0.5;
//  Cfake_efficiency=0.1;
//  Qfake_efficiency=0.01;
  Btagging_efficiency=0.8;
  Cfake_efficiency=0.0;
  Qfake_efficiency=0.0;
  cin>>setBveto>>dummy;           // rejects b-tagged events

/*******************************************************************************
 * End of USER-DEFINED SELECTION CUTS
 ******************************************************************************/



  //Information about the applied selection cuts will be stored in a canvas
  TPaveText text_cuts(0.1,0.1,0.9,0.9);
  text_cuts.SetLabel("List of selection cuts applied");
  text_cuts.AddText(" ");
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
    TString s="inv.mass window for the 2 central jets and jf-jb pair:";
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
  if(setcuts_deltaetaVV=="y"){
    TString s="minimum |#Delta#eta| between jcjc and reconstructed W:";
    s+=cut_deltaetaVV;
    text_cuts.AddText(s.Data());
  }
  if(setcuts_TotVisM=="y"){
    TString s="minimum total visible mass 4j+l:";
    s+=cut_TotVisM;
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
  if(setcuts_minDeltaRjj=="y"){
    TString s="minimum DR between jets:";
    s+=cut_minDeltaRjj;
    text_cuts.AddText(s.Data());
  }
  if(setcuts_minpT_miss=="y"){
    TString s="minimum missing pT:";
    s+=cut_minpT_miss;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setcuts_minpT_jc1=="y"){
    TString s="minimum pT_jc1:";
    s+=cut_minpT_jc1;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setcuts_minpT_jc2=="y"){
    TString s="minimum pT_jc2:";
    s+=cut_minpT_jc2;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setcuts_eta_max_jc1=="y"){
    TString s="maximum |#eta| of jc1:";
    s+=cut_eta_max_jc1;
    text_cuts.AddText(s.Data());
  }
  if(setcuts_eta_max_jc2=="y"){
    TString s="maximum |#eta| of jc2:";
    s+=cut_eta_max_jc2;
    text_cuts.AddText(s.Data());
  }
  if(setcuts_minpT_jcjc=="y"){
    TString s="minimum pT_jcjc:";
    s+=cut_minpT_jcjc;
    s+=" GeV";
    text_cuts.AddText(s.Data());
  }
  if(setcuts_Mjfjb=="y"){
    TString s="minimum invariant mass of fb jet pair:";
    s+=cut_Mjfjb;
    s+=" GeV";    
    text_cuts.AddText(s.Data());
  }
  if(setcuts_deltaetaVj=="y"){
    TString s=
    "minimum |#Delta#eta| between vector bosons and tag jets:";
    s+=cut_deltaetaVj;
    text_cuts.AddText(s.Data());
  }

  if(setBveto=="y"){
    TString s="b-quark veto active for |#eta|    <";
    s+=Btagging_region;
    s+="     with efficiency:";
    s+=Btagging_efficiency;
    text_cuts.AddText(s.Data());
    s="c-quark fake efficiency:";
    s+=Cfake_efficiency;
    s+="    lightquark/g fake efficiency:";
    s+=Qfake_efficiency;
    text_cuts.AddText(s.Data());
  }

  TString s="M_min of VVsytem: ";
  s+=Mcut;
  s+=" GeV";
  text_cuts.AddText(s.Data());
  
  text_cuts.AddText(" ");





  cout<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<"               -- analysis_LHC_4jmunu_full --                 "<<endl;
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
  f->mkdir("Basic_cuts");
  f->mkdir("otherhistos_peak");
  f->cd();
  f->Close();

  TString title;
  char app[50];
  Int_t val = (Int_t)Mcut;
  sprintf(app,"(for M_{VV} > %d GeV) ",val);



  // Set the number of bins and extrema for each histogram____________________
  nbins_M = 1400;
  nbins_MVV = 14000;
  nbins_eta=100;
  nbins_pT=1400;
  nbins_phi=nbins_eta;


  //set extrema for histograms________________________________________________
  Mlv_low = 5.0;
  Mlv_high = 14005.0;
  Mjj_low = 5.0;
  Mjj_high = 14005.0;
  Mjjlv_low=5.0;
  Mjjlv_high=14005.0;

// CHECK:nbins_pT too low same range as masses but much smaller numer of bins 
  eta_low = -7;
  eta_high = 7;
  pT_low = 5;
  pT_high = 14005;
  phi_low = 0.0;
  phi_high = TMath::Pi();


  // Declaration of histograms__________________________________________________

  // After cuts
  //-------------
  TH1D hew_Mlv_cuts("hew_Mlv_cuts",";M_{l#nu} (GeV)",nbins_M,Mlv_low,Mlv_high);
  hew_Mlv_cuts.SetTitle("Invariant mass of l#nu system with selection cuts");

  TH1D hew_Mjj_cuts("hew_Mjj_cuts",";M_{jj} (GeV)",nbins_MVV,Mjj_low,Mjj_high);
  hew_Mjj_cuts.SetTitle("Invariant mass of j_{c}j_{c} system with selection cuts");

  TH1D hew_Mjjlv_cuts("hew_Mjjlv_cuts","M_{jjl#nu}",nbins_MVV,Mjjlv_low,Mjjlv_high);
  hew_Mjjlv_cuts.SetTitle("M_{jjl#nu} with selection cuts;M_{jjl#nu} (GeV);");

 // At the jcjc peak
  //-------------
  TH1D hew_Mjjlv_peak("hew_Mjjlv_peak","M_{jjl#nu}",nbins_MVV,Mjjlv_low,Mjjlv_high);
  hew_Mjjlv_peak.SetTitle("M_{jjl#nu} at the jcjc peak;M_{jjl#nu} (GeV);");

  TH1D hew_pTlv_peak("hew_pTlv_peak","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT of leptonic W at the jcjc peak (pb/GeV)";
  hew_pTlv_peak.SetTitle(title.Data());

  TH1D hew_MinM_Vhadj("hew_MinM_Vhadj","M_{jV} (GeV)",nbins_M,Mjjlv_low,Mjjlv_high);
  title="Minimum inv. mass of tag jet + had vector boson at the jcjc peak"+TString(app);
  hew_MinM_Vhadj.SetTitle(title.Data()); 

  TH1D hew_Mvis("hew_Mvis",";M_{jj} (GeV)",nbins_M,Mjj_low,Mjj_high);
  hew_Mvis.SetTitle("Total visible Invariant mass at the jcjc peak");


  // Control plots
  //--------------
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

  //14. pT of the charged lepton
  TH1D hew_pT_mu("hew_pT_mu","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT of the charged lepton "+TString(app);
  hew_pT_mu.SetTitle(title.Data());

  //15. eta of the charged lepton
  TH1D hew_eta_mu("hew_eta_mu","#eta",nbins_eta,eta_low,eta_high);
  title="Pseudorapidity of the charged lepton "+TString(app);
  hew_eta_mu.SetTitle(title.Data());

  //16. pT_miss
  TH1D hew_pT_miss("hew_pT_miss","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT_miss "+TString(app);
  hew_pT_miss.SetTitle(title.Data());

  //17. DeltaRjj_min
  TH1D hew_DRjj_min("hew_DRjj_min","p_{T} (GeV)",nbins_eta,0,2*eta_high);
  title="DRjj_min "+TString(app);
  hew_DRjj_min.SetTitle(title.Data());

  //18. pT_jc1
  TH1D hew_pT_jc1("hew_pT_jc1","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT_jc1 "+TString(app);
  hew_pT_jc1.SetTitle(title.Data());

  //19. pT_jc2
  TH1D hew_pT_jc2("hew_pT_jc2","p_{T} (GeV)",nbins_pT,pT_low,pT_high);
  title="pT_jc2 "+TString(app);
  hew_pT_jc2.SetTitle(title.Data());

  //20. minimum inv. mass of a jet pair
  TH1D hew_MinM_jj("hew_MinM_jj","Mmin_{jj} (GeV)",nbins_M,Mjjlv_low,Mjjlv_high);
  title="Minimum inv. mass of jj pairs "+TString(app);
  hew_MinM_jj.SetTitle(title.Data()); 

  //21. eta_jc1
  TH1D hew_eta_jc1("hew_eta_jc1","#eta",nbins_eta,eta_low,eta_high);
  title="eta_jc1 "+TString(app);
  hew_eta_jc1.SetTitle(title.Data());

  //22. eta_jc2
  TH1D hew_eta_jc2("hew_eta_jc2","#eta",nbins_eta,eta_low,eta_high);
  title="eta_jc2 "+TString(app);
  hew_eta_jc2.SetTitle(title.Data());


  //Define "buffer" histograms (temporary objects)
  TH1D *hew_Mlv_cuts_buf = (TH1D*) hew_Mlv_cuts.Clone();
  TH1D *hew_Mjj_cuts_buf = (TH1D*) hew_Mjj_cuts.Clone();
  TH1D *hew_Mjjlv_cuts_buf = (TH1D*) hew_Mjjlv_cuts.Clone();
  TH1D *hew_Mjjlv_peak_buf = (TH1D*) hew_Mjjlv_peak.Clone();
  TH1D *hew_pTlv_peak_buf = (TH1D*) hew_pTlv_peak.Clone();
  TH1D *hew_Mvis_buf = (TH1D*) hew_Mvis.Clone();
  TH1D *hew_MinM_Vhadj_buf = (TH1D*) hew_MinM_Vhadj.Clone();
  TH1D *hew_pTmin_jfjb_buf = (TH1D*) hew_pTmin_jfjb.Clone();
  TH1D *hew_pTmin_alljets_buf = (TH1D*) hew_pTmin_alljets.Clone();
  TH1D *hew_pT_jcjc_buf = (TH1D*) hew_pT_jcjc.Clone();
  TH1D *hew_pT_mu_buf = (TH1D*) hew_pT_mu.Clone();
  TH1D *hew_pT_miss_buf = (TH1D*) hew_pT_miss.Clone();
  TH1D *hew_pT_jc1_buf = (TH1D*) hew_pT_jc1.Clone();
  TH1D *hew_pT_jc2_buf = (TH1D*) hew_pT_jc2.Clone();
  TH1D *hew_DRjj_min_buf = (TH1D*) hew_DRjj_min.Clone();
  TH1D *hew_eta_jcjc_buf = (TH1D*) hew_eta_jcjc.Clone();
  TH1D *hew_eta_mu_buf = (TH1D*) hew_eta_mu.Clone();
  TH1D *hew_eta_jc1_buf = (TH1D*) hew_eta_jc1.Clone();
  TH1D *hew_eta_jc2_buf = (TH1D*) hew_eta_jc2.Clone();
  TH1D *hew_M_jfjb_buf = (TH1D*) hew_M_jfjb.Clone();
  TH1D *hew_Mmin_jcjcjfb_buf = (TH1D*) hew_Mmin_jcjcjfb.Clone();
  TH1D *hew_DeltaEta_jfjb_buf = (TH1D*) hew_DeltaEta_jfjb.Clone();
  TH1D *hew_DeltaEta_VV_buf = (TH1D*) hew_DeltaEta_VV.Clone();
  TH1D *hew_MaxPT_tagj_buf = (TH1D*) hew_MaxPT_tagj.Clone();
  TH1D *hew_MinM_Vj_buf = (TH1D*) hew_MinM_Vj.Clone();
  TH1D *hew_comp_MinM_Vj_buf = (TH1D*) hew_comp_MinM_Vj.Clone();
  TH1D *hew_MinM_jj_buf = (TH1D*) hew_MinM_jj.Clone();
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
  // Take the value of total cross section (stored inside the tree as User
  // Information)
  TObjString *os_xsec = (TObjString*)tree->GetUserInfo()->At(0);
  TString s_xsec = os_xsec->GetString();
  xsec = s_xsec.Atof();
  xsec_tot += xsec;
  
  cout<<endl<<ktr+1<<"  >>>>>"<<treefilename[ktr].Data()<<"<<<<<"<<":"<<endl;
  cout<<"Read total cross section: "<<xsec<<" pb"<<endl;
  nentries_tree = (Int_t)tree->GetEntries();





  // Start analysis_____________________________________________________________
  iwarning=0;

  // Loop over events stored inside the tree
  for(Int_t k_n=0;k_n<nentries_tree;k_n++) {

    tree->GetEntry(k_n);

    n_j=0;
    n_clep=0;
    n_nu=0;
    tot_vis=TLorentzVector(0.0,0.0,0.0,0.0);

   //Make sure that the first two particles are ingoing (initial state)
   if(istup[0] != -1 || istup[1] != -1){
     cout<<"***ERROR: the first 2 particles are NOT ingoing!!!"<<endl;
     cout<<"          Please check the event file."<<endl<<endl;
     cout<<"Execution stopped"<<endl;
     exit(0);
   }

    for(Int_t k=0;k<8;k++){

      particle[k] = TLorentzVector(px[k],py[k],pz[k],E[k]);

      if(istup[k]==1){   //final-state particles
	if(TMath::Abs(idup[k]) <= 5 || idup[k]==21){   //jet
	  jet[n_j] = TLorentzVector(px[k],py[k],pz[k],E[k]);
          index_jet[n_j]=k;
          tot_vis+=jet[n_j];
	  n_j++;
	}
 	else if(TMath::Abs(idup[k])==11 || TMath::Abs(idup[k])==13 || 
		TMath::Abs(idup[k])==15){   //charged lepton
	  clep[n_clep] = TLorentzVector(px[k],py[k],pz[k],E[k]);
          tot_vis+=clep[n_clep];
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
    }  // End of loop for(Int_t k=0;k<8;k++)


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
      cout<<"ERROR:WRONG itag_scheme"<<endl;
      exit(0);
    }
    //END of Identify TAG JETS (jf,jb) and complementary ones (jc)
    
    jfjb=jf+jb;

    //Identify heavy bosons: by convention V[0] decays hadronically, V[1] 
    //leptonically.
    V[0]=jc[0]+jc[1];
    V[1]=clep[0]+nu[0];
    VV=V[0]+V[1];
    Vhad_jf = V[0]+jf;
    Vhad_jb = V[0]+jb;



    // Start selection procedure_______________________________________________

    ipasscuts=1;

    //Cut on pT and pseudorapidity of jets
    for(Int_t k=0;k<n_j;k++){
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
    //Cut on pT and pseudorapidity of charged leptons
    for(Int_t k=0;k<n_clep;k++){
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

    if(ipasscuts==1){
    // Cenral jets are always classified according to decreasing pT 
       if(jc[0].Pt()>jc[1].Pt()){
         pT_jc1=jc[0].Pt();
         pT_jc2=jc[1].Pt();
         eta_jc1=jc[0].Eta();
         eta_jc2=jc[1].Eta();
         }else{
         pT_jc1=jc[1].Pt();
         pT_jc2=jc[0].Pt();
         eta_jc1=jc[1].Eta();
         eta_jc2=jc[0].Eta();
         }
    //Cut on pT and pseudorapidity of central jets
    // ordered according to their pT jc1.Pt > jc2.Pt
      if(setcuts_minpT_jc1=="y"){
        if(pT_jc1< cut_minpT_jc1){
          ipasscuts=0;
        }
      }
      if(setcuts_minpT_jc2=="y"){
        if(pT_jc2< cut_minpT_jc2){
          ipasscuts=0;
        }
      }
      if(setcuts_eta_max_jc1=="y"){
        if( eta_jc1> cut_eta_max_jc1){
          ipasscuts=0;
        }
      }
      if(setcuts_eta_max_jc2=="y"){
        if( eta_jc2> cut_eta_max_jc2){
          ipasscuts=0;
        }
      }
    //Cut on pT of central jet pair
      pT_jcjc=V[0].Pt();
      if(setcuts_minpT_jcjc=="y"){
        if(pT_jcjc< cut_minpT_jcjc){
          ipasscuts=0;
        }
      }
    //Cut on the invariant mass of any two jets
      Mmin_jj=1.0e9;
      for(Int_t ki=0;ki<n_j;ki++){
        for(Int_t kj=ki+1;kj<n_j;kj++){
          jj = jet[ki]+jet[kj];
          Mmin_jj=TMath::Min(Mmin_jj,jj.M());
        }
      }
      if(setcuts_Mjj=="y"){
        if(Mmin_jj < cut_Mjj) {
          ipasscuts=0;
        }
      }
      //Cut on minimum |Delta R| of any two jets
      DRmin=1.e6;
      for(Int_t ki=0;ki<n_j;ki++){
	for(Int_t kj=ki+1;kj<n_j;kj++){
          DR=TMath::Abs(jet[ki].DeltaR(jet[kj]));
	  if(DR < DRmin){
	    DRmin=DR;
	  }
	}
      }
      if(setcuts_minDeltaRjj=="y"){
	    if(DRmin < cut_minDeltaRjj){
	      ipasscuts=0;
	    }
      }
//    }
    //Cut on the invariant mass of VV pair
      if(VV.M() < Mcut){
	      ipasscuts=0;
	    }

    //Cut on mass window for jf-jb pair
      if(setcuts_Mjcjc=="y"){
	if(TMath::Abs(jfjb.M()-centralvalue_Mjcjc) < cut_Mjcjc){
	  ipasscuts=0;
	}
      }
     }


    Double_t DeltaEta_jfjb = TMath::Abs(jf.Eta()-jb.Eta());
    Double_t M_jfjb = jfjb.M();
    if(ipasscuts==1) {
      //pass cuts on |Delta eta| between FORWARD/BACKWARD JETS
      if(setcuts_deltaetajfjb=="y"){
	if(DeltaEta_jfjb < cut_deltaetajfjb){
	  ipasscuts=0;
	}
      }
      if(setcuts_Mjfjb=="y"){
	if(M_jfjb < cut_Mjfjb){
	  ipasscuts=0;
	}
      }
    }

    DeltaEta_VV = TMath::Abs(V[0].Eta()-V[1].Eta());
    if(ipasscuts==1) {
      if(setcuts_deltaetaVV=="y"){
         if(DeltaEta_VV < cut_deltaetaVV){
	    ipasscuts=0;
	}
      }
      if(setcuts_TotVisM=="y"){
         if(tot_vis.M()<cut_TotVisM){
	    ipasscuts=0;
	}
      }
    }

    if(ipasscuts==1){
      //Cuts on missing pT
      if(setcuts_minpT_miss=="y"){
	if(nu[0].Pt() < cut_minpT_miss){
	  ipasscuts=0;
	}
      }
    }

    //Cut on DEltaEta between tag jets and vector bosons
    Double_t MinDEta_Vtag=0.0;
    Double_t DEta_jfjcjc = TMath::Abs(jf.Eta() - V[0].Eta());
    Double_t DEta_jbjcjc = TMath::Abs(jb.Eta() - V[0].Eta());
    Double_t DEta_jflv = TMath::Abs(jf.Eta() - V[1].Eta());
    Double_t DEta_jblv = TMath::Abs(jb.Eta() - V[1].Eta());
    MinDEta_Vtag=TMath::Min(DEta_jfjcjc,DEta_jbjcjc);
    MinDEta_Vtag=TMath::Min(MinDEta_Vtag,DEta_jflv);
    MinDEta_Vtag=TMath::Min(MinDEta_Vtag,DEta_jblv);
    if(ipasscuts==1) {
      if(setcuts_deltaetaVj=="y"){
        if(MinDEta_Vtag < cut_deltaetaVj){
	    ipasscuts=0;
	}
      }
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

    } //end if(ipasscuts==1) {



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

    if(ipasscuts==1) {
      hew_Mlv_cuts_buf->Fill(V[1].M());
      hew_Mjj_cuts_buf->Fill(V[0].M());
      hew_Mjjlv_cuts_buf->Fill(VV.M());
            }

    if(ipasscuts==1){

      //FILL OTHER HISTOGRAMS AT THE jcjc PEAK_______________________________
      if(TMath::Abs(V[0].M()-centralvalue_Mjcjc) < cut_Mjcjc
        || setcuts_Mjcjc=="n"){

        hew_Mjjlv_peak_buf->Fill(VV.M());
        hew_pTlv_peak_buf->Fill(V[1].Pt());
        hew_Mvis_buf->Fill(tot_vis.M());
        hew_MinM_Vhadj_buf->Fill(TMath::Min(Vhad_jf.M(),Vhad_jb.M()));
        hew_pT_mu_buf->Fill(clep[0].Pt());
        hew_eta_mu_buf->Fill(clep[0].Eta());
        hew_pT_miss_buf->Fill(nu[0].Pt());
        hew_DRjj_min_buf->Fill(DRmin);
        hew_MinM_jj_buf->Fill(Mmin_jj);
        
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
	hew_pT_jcjc_buf->Fill(pT_jcjc);

       	//4. eta of the heavy boson decaying hadronically (=> two central jets)
	hew_eta_jcjc_buf->Fill(V[0].Eta());

      	//5. invariant mass of the two tag jets
	hew_M_jfjb_buf->Fill(M_jfjb);

	//6. min. inv. mass of the heavy boson decaying hadronically and one tag
	//   jet
	Double_t M_jcjcjfb_min = TMath::Min(Vhad_jf.M(),Vhad_jb.M());
	hew_Mmin_jcjcjfb_buf->Fill(M_jcjcjfb_min);

	//7. |Delta eta| between tag jets
	hew_DeltaEta_jfjb_buf->Fill(DeltaEta_jfjb);

	//8. |Delta eta| between reconstructed heavy bosons
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
        
	//18. Largest pT of the two central jets
        hew_pT_jc1_buf->Fill(pT_jc1);
        
	//19. Lowest pT of the two central jets
        hew_pT_jc2_buf->Fill(pT_jc2);
        
	//21. eta_jc1
        hew_eta_jc1_buf->Fill(eta_jc1);
        
	//22. eta_jc2
        hew_eta_jc2_buf->Fill(eta_jc2);
        
        
      }     // end of  if(TMath::Abs(V[0].M()-centralvalue_Mjcjc) < cut_Mjcjc)
    }       // end of if(ipasscuts==1)


  } //end for(Int_t k_n=0;k_n<nentries_tree;k_n++)


  cout<<"WARNING: discarded "<<iwarning<<" events with no real solutions for reconstructed pz of neutrino in file "<<treefilename[ktr].Data()<<endl;

  treefile.Close();


  // Normalize histograms to the total cross section____________________________
  // Scale factors for histograms; having been multiplied by them, 
  // histograms show the DIFFERENTIAL CROSS SECTION.
  // [scale = sigma_tot/(A*nentries)]
  A_Mlv = (Mlv_high - Mlv_low)/nbins_M;
  scale_Mlv = xsec/(A_Mlv*(nentries_tree));
  A_Mjj = (Mjj_high - Mjj_low)/nbins_M;
  scale_Mjj = xsec/(A_Mjj*(nentries_tree));
  A_Mjjlv = (Mjjlv_high - Mjjlv_low)/nbins_MVV;
  scale_Mjjlv = xsec/(A_Mjjlv*(nentries_tree));
  
  A_eta = (eta_high - eta_low)/nbins_eta;
  scale_eta = xsec/(A_eta*(nentries_tree));
  A_pT = (pT_high - pT_low)/nbins_pT;
  scale_pT = xsec/(A_pT*(nentries_tree));
  A_phi = (phi_high - phi_low)/nbins_phi;
  scale_phi = xsec/(A_phi*(nentries_tree));

  hew_Mlv_cuts_buf->Scale(scale_Mlv);
  hew_Mjj_cuts_buf->Scale(scale_Mjjlv);
  hew_pTlv_peak_buf->Scale(scale_pT);
  hew_MinM_Vhadj_buf->Scale(scale_Mjjlv);
  hew_Mvis_buf->Scale(scale_Mjjlv);
  hew_Mjjlv_cuts_buf->Scale(scale_Mjjlv);
  hew_Mjjlv_peak_buf->Scale(scale_Mjjlv);
  hew_pTmin_jfjb_buf->Scale(scale_pT);
  hew_pTmin_alljets_buf->Scale(scale_pT);
  hew_pT_jcjc_buf->Scale(scale_pT);
  hew_pT_mu_buf->Scale(scale_pT);
  hew_pT_miss_buf->Scale(scale_pT);
  hew_pT_jc1_buf->Scale(scale_pT);
  hew_pT_jc2_buf->Scale(scale_pT);
  hew_DRjj_min_buf->Scale(scale_eta);
  hew_eta_jcjc_buf->Scale(scale_eta);
  hew_eta_mu_buf->Scale(scale_eta);
  hew_eta_jc1_buf->Scale(scale_eta);
  hew_eta_jc2_buf->Scale(scale_eta);
  hew_M_jfjb_buf->Scale(scale_Mjj);
  hew_Mmin_jcjcjfb_buf->Scale(scale_Mjj);
  hew_DeltaEta_jfjb_buf->Scale(xsec/(((2*eta_high)/nbins_eta)*(nentries_tree)));
  hew_DeltaEta_VV_buf->Scale(xsec/(((2*eta_high)/nbins_eta)*(nentries_tree)));
  hew_MaxPT_tagj_buf->Scale(scale_pT);
  hew_MinM_Vj_buf->Scale(scale_Mjj);
  hew_comp_MinM_Vj_buf->Scale(scale_Mjj);
  hew_MinM_jj_buf->Scale(scale_Mjj);
  hew_MinDEta_Vtag_buf->Scale(scale_eta);
  hew_comp_MinDEta_Vtag_buf->Scale(scale_eta);
  hew_DPhi_tag_buf->Scale(scale_phi);
  hew_DPhi_V_buf->Scale(scale_phi);


  //Unbuffer temporary histograms (clones)
  hew_Mlv_cuts.Add(hew_Mlv_cuts_buf);
  hew_Mjj_cuts.Add(hew_Mjj_cuts_buf);
  hew_pTlv_peak.Add(hew_pTlv_peak_buf);
  hew_MinM_Vhadj.Add(hew_MinM_Vhadj_buf);
  hew_Mvis.Add(hew_Mvis_buf);
  hew_Mjjlv_cuts.Add(hew_Mjjlv_cuts_buf);
  hew_Mjjlv_peak.Add(hew_Mjjlv_peak_buf);
  hew_pTmin_jfjb.Add(hew_pTmin_jfjb_buf);
  hew_pTmin_alljets.Add(hew_pTmin_alljets_buf);
  hew_pT_jcjc.Add(hew_pT_jcjc_buf);
  hew_pT_mu.Add(hew_pT_mu_buf);
  hew_pT_miss.Add(hew_pT_miss_buf);
  hew_pT_jc1.Add(hew_pT_jc1_buf);
  hew_pT_jc2.Add(hew_pT_jc2_buf);
  hew_DRjj_min.Add(hew_DRjj_min_buf);
  hew_pT_mu.Add(hew_pT_mu_buf);
  hew_eta_jcjc.Add(hew_eta_jcjc_buf);
  hew_eta_mu.Add(hew_eta_mu_buf);
  hew_eta_jc1.Add(hew_eta_jc1_buf);
  hew_eta_jc2.Add(hew_eta_jc2_buf);
  hew_M_jfjb.Add(hew_M_jfjb_buf);
  hew_Mmin_jcjcjfb.Add(hew_Mmin_jcjcjfb_buf);
  hew_DeltaEta_jfjb.Add(hew_DeltaEta_jfjb_buf);
  hew_DeltaEta_VV.Add(hew_DeltaEta_VV_buf);
  hew_MaxPT_tagj.Add(hew_MaxPT_tagj_buf);
  hew_MinM_Vj.Add(hew_MinM_Vj_buf);
  hew_comp_MinM_Vj.Add(hew_comp_MinM_Vj_buf);
  hew_MinM_jj.Add(hew_MinM_jj_buf);
  hew_MinDEta_Vtag.Add(hew_MinDEta_Vtag_buf);
  hew_comp_MinDEta_Vtag.Add(hew_comp_MinDEta_Vtag_buf);
  hew_DPhi_tag.Add(hew_DPhi_tag_buf);
  hew_DPhi_V.Add(hew_DPhi_V_buf);

  hew_Mlv_cuts_buf->Reset();
  hew_Mjj_cuts_buf->Reset();
  hew_pTlv_peak_buf->Reset();
  hew_MinM_Vhadj_buf->Reset();
  hew_Mvis_buf->Reset();
  hew_Mjjlv_cuts_buf->Reset();
  hew_Mjjlv_peak_buf->Reset();
  hew_pTmin_jfjb_buf->Reset();
  hew_pTmin_alljets_buf->Reset();
  hew_pT_jcjc_buf->Reset();
  hew_pT_mu_buf->Reset();
  hew_pT_miss_buf->Reset();
  hew_pT_jc1_buf->Reset();
  hew_pT_jc2_buf->Reset();
  hew_DRjj_min_buf->Reset();
  hew_eta_jcjc_buf->Reset();
  hew_eta_mu_buf->Reset();
  hew_eta_jc1_buf->Reset();
  hew_eta_jc2_buf->Reset();
  hew_M_jfjb_buf->Reset();
  hew_Mmin_jcjcjfb_buf->Reset();
  hew_DeltaEta_jfjb_buf->Reset();
  hew_DeltaEta_VV_buf->Reset();
  hew_MaxPT_tagj_buf->Reset();
  hew_MinM_Vj_buf->Reset();
  hew_comp_MinM_Vj_buf->Reset();
  hew_MinM_jj_buf->Reset();
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
//  TString s;
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

  fo->cd("Basic_cuts");
  hew_Mlv_cuts.Write();
  hew_Mjj_cuts.Write();
  hew_Mjjlv_cuts.Write();
  fo->cd();
  fo->cd("otherhistos_peak");
  hew_pTlv_peak.Write();
  hew_MinM_Vhadj.Write();
  hew_Mvis.Write();
  hew_pTmin_jfjb.Write();
  hew_pTmin_alljets.Write();
  hew_pT_jcjc.Write();
  hew_pT_mu.Write();
  hew_pT_miss.Write();
  hew_pT_jc1.Write();
  hew_pT_jc2.Write();
  hew_DRjj_min.Write();
  hew_eta_jcjc.Write();
  hew_eta_mu.Write();
  hew_eta_jc1.Write();
  hew_eta_jc2.Write();
  hew_M_jfjb.Write();
  hew_Mmin_jcjcjfb.Write();
  hew_DeltaEta_jfjb.Write();
  hew_DeltaEta_VV.Write();
  hew_MaxPT_tagj.Write();
  hew_MinM_Vj.Write();
  hew_comp_MinM_Vj.Write();
  hew_MinM_jj.Write();
  hew_MinDEta_Vtag.Write();
  hew_comp_MinDEta_Vtag.Write();
  hew_DPhi_tag.Write();
  hew_DPhi_V.Write();
  hew_Mjjlv_peak.Write();
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

  cout<<"integrale"<<endl<<"da "<<low<<" (bin "<<lowBin<<")"<<endl
      <<"a "<<high<<" (bin "<<highBin<<")"<<endl; 
  integral= (histo->Integral(lowBin,highBin))*width;
  return integral;
}
