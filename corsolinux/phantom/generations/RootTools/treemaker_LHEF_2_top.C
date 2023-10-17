/*
 * treemaker_LHEF_2_top.C
 *
 * Last update: May 26, 2008
 *
 * This routine stores the information from a Les Houches Event File (LHEF) in a
 * ROOT file tree.root.
 * The file consists of a ROOT tree which provides, for each EXTERNAL
 * (i.e. not intermediate) particle, the following information:
 *  - IDUP
 *  - ISTUP
 *  - px
 *  - py
 *  - pz
 *  - E
 * Moreover the tree stores the value of the total cross section, which is
 * required to be entered as input.
 *
 * This particular version keeps only the events with at least two b-quarks
 * in the final state and at least one top in order to combine this sample
 * with the enriched sample H-->bb in which top's are vetoed 
 ******************************************************************************/


#if !defined(_CINT_) || defined(_MAKECINT_)
#include <Riostream.h>
#include <stdio.h>
#include "TString.h"
#include "TObjString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#endif



void treemaker()
{

  gSystem->Load("libPhysics.so");
  gROOT->Reset();

  
  Double_t px[8], py[8], pz[8], E[8];
  Int_t idup[8];             //Les Houches IDUP (particle ID)
  Int_t istup[8];            //Les Houches ISTUP
  Int_t imotherstatus;       //Les Houches mother status:
                             // if abs(imotherstatus)!=1, particle is a mother
                             // i.e. it does not belong to the effective final 
                             // state
  Double_t xsec,xsec_top;             //total cross section
  TLorentzVector jet[8], clep[8], nu[8], Vlep, jjj, jlv;
  Int_t n_j, n_clep, n_nu, n_top, n_b, ipasscuts;
  Double_t centralvalue_Mjjj, cut_Mjjj, centralvalue_Mjlv_rec, cut_Mjlv_rec;
  
  TString LHEFfilename;
  char line_aux[300];
  char str[300];
  TString line;
  Double_t dummy_double[7];

  Int_t nline, nline_local, ninitlines, nevts, nparticle, nevts_accepted;
  nline = 0;
  ninitlines = 0;
  nevts = 0;
  nevts_accepted = 0;
  
  centralvalue_Mjjj = 175.0;
  cut_Mjjj = 15.0;
  centralvalue_Mjlv_rec = 175.0;
  cut_Mjlv_rec = 15.0;
  
  cout.precision(10);

  cout<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<"                    -- treemaker_LHEF --                      "<<endl;
  cout<<"   version adapted to Les Houches Event File (LHEF) format    "<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<endl;

  cout<<"Enter the number of LHEF files to be read: "<<endl;
  Int_t nlhef=0;
  cin>>nlhef;
  cout<<"Enter the  absolute path of the "<< nlhef <<" LHEF files: "<<endl;
  TString gen[nlhef];
  for(Int_t kfile=0;kfile<nlhef;kfile++){
    cin>>gen[kfile];
  }


  TString rootfilename = "tree.root";
  TFile *rootfile = new TFile(rootfilename.Data(),"RECREATE");
  cout<<"created ROOT file "<<rootfilename.Data()<<endl;

  // create tree and set branches
  TTree *tree = new TTree("events","Kinematics of the event (four-momenta and particle ID's)");
  tree->Branch("px",px,"px[8]/D");
  tree->Branch("py",py,"py[8]/D");
  tree->Branch("pz",pz,"pz[8]/D");
  tree->Branch("E",E,"E[8]/D");
  tree->Branch("IDUP",idup,"idup[8]/I");
  tree->Branch("ISTUP",istup,"istup[8]/I");


  // store information about the total cross section
  cout<<"Enter the value of total cross section: ";
  cin>>xsec;
  cout.precision(10);
  cout<<endl;
  cout<<"read total cross section: xs = "<<xsec<<" pb"<<endl;
  cout<<endl;
  cout<<endl;
  cout<<"Processing. Please wait..."<<endl;




  // Start analysis_____________________________________________________________


  for(Int_t kfile=0;kfile<nlhef;kfile++){

    LHEFfilename = gen[kfile];

    // Check if the input file is in Les Houches (LHEF) format. LHEF requires 
    // the first line to look like: "<LesHouchesEvents version="X.X">"
    FILE *fp0 = fopen(LHEFfilename.Data(),"r");
    fgets(&line_aux[0],300,fp0);
    sscanf(&line_aux[0], "%s",&str[0]);
    line = str;
    if(!line.Contains("<LesHouchesEvents")) {
      cout<<"***ERROR: event file is not in standard LHEF format.!!!"<<endl;
      exit(0);
    }
    fclose(fp0);


  // Scan LHEF file____________________________________________________________

    FILE *fp = fopen(LHEFfilename.Data(),"r");

    ninitlines=0;
    line="";
    while(line != "<event>") {
      fgets(&line_aux[0],300,fp);
      sscanf(&line_aux[0], "%s",&str[0]);
      line = str;
      ninitlines++;
    }

    //At this point line contains the first "<event>" encountered: this is the
    //starting point for data acquisition. Everything must be initialized because the line <event> will not be read again
    nline = 0;
    nline_local = 1;
    nparticle = 0;
    n_j=0;
    n_clep=0;
    n_nu=0;
    n_b=0;
    
    while (line != "</LesHouchesEvents>") {

      fgets(&line_aux[0],300,fp);
      sscanf(&line_aux[0], "%s",&str[0]);
      line = str;
      nline++;

      if(line.Contains("<event>")){
	nline_local=1;
      n_j=0;
      n_clep=0;
      n_nu=0;
      n_b=0;
	continue;  //returns at beginning of the WHILE loop
      }
      else if(line.BeginsWith("#")){
	continue;
      }
      else if(nline_local==1){ // line contains common event information
	nline_local++;
	continue;
      }
      else if(line.Contains("</event>")){

        if(n_b>1){
// Check for top
        Vlep=clep[0]+nu[0];
	  ipasscuts=0;
        n_top=0;
	  // invariant mass of jet triplets {jjj}
	  for(Int_t ki=0;ki<n_j;ki++){
	    for(Int_t kj=ki+1;kj<n_j;kj++){
	      for(Int_t kk=kj+1;kk<n_j;kk++){
	        jjj=jet[ki]+jet[kj]+jet[kk];
	        if(TMath::Abs(jjj.M()-centralvalue_Mjjj) < cut_Mjjj){
	      	ipasscuts=1;
	        }
	      }
	    }
	  }
        if(ipasscuts==1){
          n_top++;
          ipasscuts=0;
        }
	  // invariant mass of Vlep+jet
	  for(Int_t ki=0;ki<n_j;ki++){
	    jlv=jet[ki]+Vlep;
	    if(TMath::Abs(jlv.M()-centralvalue_Mjlv_rec) < cut_Mjlv_rec){
	      ipasscuts=1;
	    }
	  }
        if(ipasscuts==1){
          n_top++;
        }
        }

        if(n_b<=1 || n_top>0){
	  tree->Fill();
	  nevts_accepted++;
        }
	nevts++;
	nparticle=0;
	nline_local=0;
      }
      else{  // analyse line content
	nline_local++;
	sscanf(&line_aux[0], "%*i %d %*i %*i %*i %*i %le %le %le %le %le %le %le",&imotherstatus,&dummy_double[0],&dummy_double[1],&dummy_double[2],&dummy_double[3],&dummy_double[4],&dummy_double[5],&dummy_double[6]);

	if(TMath::Abs(imotherstatus)==1){
	  sscanf(&line_aux[0], "%i %d %*i %*i %*i %*i %le %le %le %le %le %le %le",&idup[nparticle],&istup[nparticle],&px[nparticle],&py[nparticle],&pz[nparticle],&E[nparticle],&dummy_double[0],&dummy_double[1],&dummy_double[2]);


         if(istup[nparticle]==1){   //final-state particles
	     if(TMath::Abs(idup[nparticle]) < 5 || idup[nparticle]==21){   //jet
	       jet[n_j] = TLorentzVector(px[nparticle],py[nparticle],pz[nparticle],E[nparticle]);
	       n_j++;
	     }
	     if(TMath::Abs(idup[nparticle]) == 5) {   //b
	       jet[n_j] = TLorentzVector(px[nparticle],py[nparticle],pz[nparticle],E[nparticle]);
	       n_j++;
             n_b++;
	     }
 	     else if(TMath::Abs(idup[nparticle])==11 || TMath::Abs(idup[nparticle])==13 || 
		  TMath::Abs(idup[nparticle])==15){   //charged lepton
	       clep[n_clep] = TLorentzVector(px[nparticle],py[nparticle],pz[nparticle],E[nparticle]);
	       n_clep++;
	     }
	     else if(TMath::Abs(idup[nparticle])==12 || TMath::Abs(idup[nparticle])==14 || 
	    	TMath::Abs(idup[nparticle])==16){   //neutrino
	        nu[n_nu] = TLorentzVector(px[nparticle],py[nparticle],pz[nparticle],E[nparticle]);
	        n_nu++;
	      }
	   }  // End if(istup[nparticle]==1)


	  nparticle++;
      } // End if(TMath::Abs(imotherstatus)==1)

     }

    }  // end  while (line != "</LesHouchesEvents>") {


    fclose(fp);

    cout<<endl<<kfile+1<<"  >>>>>"<<gen[kfile].Data()<<"<<<<<"<<":"<<endl;
    cout<<nline+ninitlines<<" lines read in "<<LHEFfilename.Data()<<endl;
    cout<<nevts<<" events processed"<<endl;
    cout<<endl;
  
  } //end for(Int_t kfile=0;kfile<nlhef;kfile++){


   cout<<"Total cross section "<<xsec<<" pb"<<endl;
// Rescale cross section
  xsec_top=xsec*nevts_accepted/nevts;
   cout<<"Total cross section with top(s) "<<xsec_top<<" pb"<<endl;
  TString sxsec;
  sxsec.Form("%f",xsec_top);
  TObjString *osxsec = new TObjString(sxsec.Data());
  tree->GetUserInfo()->Add(osxsec);
  
  rootfile->Write();

  cout<<endl;
  cout<<"Execution successful."<<endl;
  cout<<endl;


}
