/*
 *treemaker_LHEF.C
 *
 *Last update: Aug 29, 2007
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
 ******************************************************************************/

#if !defined(_CINT_) || defined(_MAKECINT_)
#include <Riostream.h>
#include <stdio.h>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TObjString.h"
#endif

void treemaker()
{

  gROOT->Reset();

  
  Double_t px[8], py[8], pz[8], E[8];
  Int_t idup[8];             //Les Houches IDUP (particle ID)
  Int_t istup[8];            //Les Houches ISTUP
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

  cout<<"Enter absolute path of directory containing LHEF files: "<<endl;
  cout<<">> ";
  TString basedir;
  cin>>basedir;
  cout<<"Enter the number of LHEF files to be read (starting by gen1/): "<<endl;
  Int_t nlhef=0;
  cin>>nlhef;
  TString gen[nlhef];
  for(Int_t ki=0;ki<nlhef;ki++){
    gen[ki] = "gen";
    gen[ki]+=(ki+1);
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
  TString sxsec;
  sxsec.Form("%f",xsec);
  TObjString *osxsec = new TObjString(sxsec.Data());
  tree->GetUserInfo()->Add(osxsec);
  
  cout<<endl;
  cout<<"Processing. Please wait..."<<endl;




  // Start analysis_____________________________________________________________


  for(Int_t ki=0;ki<nlhef;ki++){

    LHEFfilename = basedir+"/"+gen[ki]+"/phamom.dat";

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
