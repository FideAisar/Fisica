/*
 *treemaker_LHEF_3.C
 * 
 *
 *Last update: March 1, 2011
 *
 *This routine stores the information from a Les Houches Event File (LHEF) in a
 *ROOT file tree.root.
 *The file consists of a ROOT tree which provides, for each EXTERNAL
 *(i.e. not intermediate) particle, the following information:
 *  - IDUP
 *  - ISTUP
 *  - px
 *  - py
 *  - pz
 *  - E
 *Moreover the tree stores the value of the total cross section, which is
 *required to be entered as input.
 *
 * Max number of particles increased to 12
 ******************************************************************************/

#if !defined(_CINT_) || defined(_MAKECINT_)
#include "TROOT.h"
#include <Riostream.h>
#include <stdio.h>
#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObjString.h"
#include "TPaveText.h"
#include "TCanvas.h"
#endif

void treemaker()
{

//  gROOT->Reset();

  
  Int_t nmax_partons=12;
  Double_t px[nmax_partons], py[nmax_partons], pz[nmax_partons], E[nmax_partons];
  Int_t idup[nmax_partons];             //Les Houches IDUP (particle ID)
  Int_t istup[nmax_partons];            //Les Houches ISTUP
  Int_t imotherstatus;       //Les Houches mother status:
                             // if abs(imotherstatus)!=1, particle is a mother
                             // i.e. it does not belong to the effective final 
                             // state
  Double_t xsec;             //total cross section

  TString LHEFfilename;
  char line_aux[300];
  char str[300];
  TString line;
  Double_t dummy_double[7];
  TString dummy_string; 
  Int_t nmax_lhef=20;
  TString gen[20];

  Int_t nline, nline_local, ninitlines, nevts, nparticle;
  nline = 0;
  ninitlines = 0;
  nevts = 0; 

  cout.precision(10);

  cout<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<"                    -- treemaker_LHEF --                      "<<endl;
  cout<<"   version adapted to Les Houches Event File (LHEF) format    "<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<endl;

  cout<<"Enter the number of LHEF files to be read (<"<<nmax_lhef<<": "<<endl;
  Int_t nlhef=0;
  cin>>nlhef>>dummy_string;
  if(nlhef > nmax_lhef){
  cout<<"Max Number of LHEF files ("<<nmax_lhef<<") exceeded! Aborting!"<<endl;
  return;
  }
  else{
   cout<<"Number of LHEF files to be read is: "<<nlhef<<endl;
    cout<<"Enter the  absolute path of the "<< nlhef <<" LHEF files: "<<endl;
    for(Int_t ki=0;ki<nlhef;ki++){
      cin>>gen[ki];
    }
  }

  TString rootfilename = "tree.root";
  TFile *rootfile = new TFile(rootfilename.Data(),"RECREATE");
  cout<<"created ROOT file "<<rootfilename.Data()<<endl;

  cout<<"Enter the value of total cross section: ";
  cin>>xsec>>dummy_string;
  cout.precision(10);
  cout<<endl;
  cout<<"read total cross section: xs = "<<xsec<<" pb"<<endl;
  cout<<endl;

  // create tree and set branches
  TTree *tree = new TTree("events","Kinematics of the event (four-momenta and particle ID's)");
  tree->Branch("px",px,"px[12]/D");
  tree->Branch("py",py,"py[12]/D");
  tree->Branch("pz",pz,"pz[12]/D");
  tree->Branch("E",E,"E[12]/D");
  tree->Branch("IDUP",idup,"idup[12]/I");
  tree->Branch("ISTUP",istup,"istup[12]/I");


  // store information about the total cross section
  TString sxsec;
  sxsec.Form("%f",xsec);
  TObjString *osxsec = new TObjString(sxsec.Data());
  tree->GetUserInfo()->Add(osxsec);
  
  //Information about the integrated cross sections will be stored in a canvas
  TPaveText text_xsec(0.1,0.1,0.9,0.9);
  text_xsec.SetLabel("Integrated cross sections");
  TString s;
  s="Total cross section:";
  s+="    ";
  s+=xsec;
  s+=" pb";
  text_xsec.AddText(s.Data());
  text_xsec.AddText(" ");

  TCanvas c_xsec("total_cross_section");
  text_xsec.Draw();
  c_xsec.Write();


  // Start analysis_____________________________________________________________
  cout<<"Starting analysis "<<endl;


  for(Int_t ki=0;ki<nlhef;ki++){

    LHEFfilename = gen[ki];

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
    //starting point for data acquisition.
    nline = 0;
    nline_local = 1;
    nevts = 0;
    nparticle = 0;

    while (line != "</LesHouchesEvents>") {

      fgets(&line_aux[0],300,fp);
      sscanf(&line_aux[0], "%s",&str[0]);
      line = str;
      nline++;

      if(line.Contains("<event>")){
	nline_local=1;
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
	tree->Fill();
	nparticle=0;
	nevts++;
	nline_local=0;
      }
      else{  // analyse line content
	nline_local++;
	sscanf(&line_aux[0], "%*i %d %*i %*i %*i %*i %le %le %le %le %le %le %le",&imotherstatus,&dummy_double[0],&dummy_double[1],&dummy_double[2],&dummy_double[3],&dummy_double[4],&dummy_double[5],&dummy_double[6]);

	if(TMath::Abs(imotherstatus)==1){
	  sscanf(&line_aux[0], "%i %d %*i %*i %*i %*i %le %le %le %le %le %le %le",&idup[nparticle],&istup[nparticle],&px[nparticle],&py[nparticle],&pz[nparticle],&E[nparticle],&dummy_double[0],&dummy_double[1],&dummy_double[2]);
	  nparticle++;
	}
      }


    }  // end  while (line != "</LesHouchesEvents>") {


    fclose(fp);

    cout<<endl<<ki+1<<"  >>>>>"<<gen[ki].Data()<<"<<<<<"<<":"<<endl;
    cout<<nline+ninitlines<<" lines read in "<<LHEFfilename.Data()<<endl;
    cout<<nevts<<" events processed"<<endl;
    cout<<endl;
  
  } //end for(Int_t ki=0;ki<nlhef;ki++){



  rootfile->Write();

  cout<<endl;
  cout<<"Execution successful."<<endl;
  cout<<endl;


}
