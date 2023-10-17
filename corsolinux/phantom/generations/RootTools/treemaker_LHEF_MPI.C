/*
 *treemaker_LHEF_2.C
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
 *This particular version is geared to deal with multiple particle interactions.
 *It combines two separate generations event by event superimposing the two sets
 *of particles. This is the most  naive approach to MPI.
 *It is assumed that the numer of LH files and of events is the same in the two
 *generations.
 ******************************************************************************/

#if !defined(_CINT_) || defined(_MAKECINT_)
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

void treemakerMPI()
{

  gROOT->Reset();

  
  Double_t px[10], py[10], pz[10], E[10];
  Int_t idup[10];             //Les Houches IDUP (particle ID)
  Int_t istup[10];            //Les Houches ISTUP
  Int_t imotherstatus;       //Les Houches mother status:
                             // if abs(imotherstatus)!=1, particle is a mother
                             // i.e. it does not belong to the effective final 
                             // state
  Double_t xsec,xsec1,xsec2,xsec_eff;             //total cross section

  TString LHEFfilename1,LHEFfilename2;
  char line_aux[300];
  char str[300];
  TString line1,line2;
  Double_t dummy_double[7];

  xsec_eff=14.5e9;  // In pb (14.5 mb)
  
  Int_t nline1, nline_local1, nline2, nline_local2, ninitlines1, ninitlines2, nevts, nparticle;
  nline1 = 0;
  ninitlines1 = 0;
  nline2 = 0;
  ninitlines2 = 0;
  nevts = 0;

  cout.precision(10);

  cout<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<"                    -- treemaker_LHEF --                      "<<endl;
  cout<<"   version adapted to Les Houches Event File (LHEF) format    "<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<endl;

  cout<<"Enter the number of LHEF files of type 1 to be read: "<<endl;
  Int_t nlhef1=0;
  cin>>nlhef1;
  cout<<"Enter the  absolute path of the "<< nlhef1 <<" LHEF files: "<<endl;
  TString gen1[nlhef1];
  for(Int_t ki=0;ki<nlhef1;ki++){
    cin>>gen1[ki];
  }
  cout<<"Enter the value of total cross section1: ";
  cin>>xsec1;
  cout.precision(10);
  cout<<endl;
  cout<<"read total cross section: xs1 = "<<xsec1<<" pb"<<endl;
  cout<<endl;

  cout<<"Enter the number of LHEF files of type 2 to be read: "<<endl;
  Int_t nlhef2=0;
  cin>>nlhef2;
  cout<<"Enter the  absolute path of the "<< nlhef2 <<" LHEF files: "<<endl;
  TString gen2[nlhef2];
  for(Int_t ki=0;ki<nlhef2;ki++){
    cin>>gen2[ki];
  }
  cout<<"Enter the value of total cross section2: ";
  cin>>xsec2;
  cout.precision(10);
  cout<<endl;
  cout<<"read total cross section: xs2 = "<<xsec2<<" pb"<<endl;
  cout<<endl;

  xsec=xsec1*xsec2/xsec_eff;
  cout<<"total cross section: xsec=xsec1*xsec2/xsec_eff = "<<xsec<<" pb"<<endl;
  cout<<endl;
  cout<<"Processing. Please wait..."<<endl;

  if(nlhef1 != nlhef2){
    cout<<"Different number of files in the two generations!!!!"<<endl;
    }

  TString rootfilename = "tree.root";
  TFile *rootfile = new TFile(rootfilename.Data(),"RECREATE");
  cout<<"created ROOT file "<<rootfilename.Data()<<endl;

  // create tree and set branches
  TTree *tree = new TTree("events","Kinematics of the event (four-momenta and particle ID's)");
  tree->Branch("px",px,"px[10]/D");
  tree->Branch("py",py,"py[10]/D");
  tree->Branch("pz",pz,"pz[10]/D");
  tree->Branch("E",E,"E[10]/D");
  tree->Branch("IDUP",idup,"idup[10]/I");
  tree->Branch("ISTUP",istup,"istup[10]/I");


  // store information about the total cross section
  TString sxsec;
  sxsec.Form("%f",xsec);
  TObjString *osxsec = new TObjString(sxsec.Data());
  tree->GetUserInfo()->Add(osxsec);
  
  //Information about the integrated cross sections will be stored in a canvas
  TPaveText text_xsec(0.1,0.1,0.9,0.9);
  text_xsec.SetLabel("Integrated cross sections");
  TString s;
  s="total cross section:";
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


  for(Int_t ki=0;ki<nlhef1;ki++){

    LHEFfilename1 = gen1[ki];

    // Check if the input file is in Les Houches (LHEF) format. LHEF requires 
    // the first line to look like: "<LesHouchesEvents version="X.X">"
    FILE *fp01 = fopen(LHEFfilename1.Data(),"r");
    fgets(&line_aux[0],300,fp01);
    sscanf(&line_aux[0], "%s",&str[0]);
    line1 = str;
    if(!line1.Contains("<LesHouchesEvents")) {
      cout<<"***ERROR: event file "<< LHEFfilename1 <<" is not in standard LHEF format.!!!"<<endl;
      exit(0);
    }
    fclose(fp01);

  cout<<"Processing file "<<LHEFfilename1<<endl;

    LHEFfilename2 = gen2[ki];

    // Check if the input file is in Les Houches (LHEF) format. LHEF requires 
    // the first line to look like: "<LesHouchesEvents version="X.X">"
    FILE *fp02 = fopen(LHEFfilename2.Data(),"r");
    fgets(&line_aux[0],300,fp02);
    sscanf(&line_aux[0], "%s",&str[0]);
    line2 = str;
    if(!line2.Contains("<LesHouchesEvents")) {
      cout<<"***ERROR: event file"<< LHEFfilename2 <<" is not in standard LHEF format.!!!"<<endl;
      exit(0);
    }
    fclose(fp02);

  cout<<"Processing file "<<LHEFfilename2<<endl;

  // Scan LHEF files____________________________________________________________

    FILE *fp1 = fopen(LHEFfilename1.Data(),"r");
    FILE *fp2 = fopen(LHEFfilename2.Data(),"r");


    ninitlines1=0;
    line1="";
    while(line1 != "<event>") {
      fgets(&line_aux[0],300,fp1);
      sscanf(&line_aux[0], "%s",&str[0]);
      line1 = str;
/*
      if(line1.Contains(EOF)){
          cout<<"File "<<LHEFfilename1.Data()<<"does not contain events!!!"<<endl;
          exit(1);
          }
*/
      ninitlines1++;
    }

    ninitlines2=0;
    line2="";
    while(line2 != "<event>") {
      fgets(&line_aux[0],300,fp2);
      sscanf(&line_aux[0], "%s",&str[0]);
      line2 = str;
      ninitlines2++;
    }

    //At this point line1,2 contain the first "<event>" encountered: this is the
    //starting point for data acquisition.
    nline1 = 0;
    nline_local1 = 1;
    nline2 = 0;
    nline_local2 = 1;
    nevts = 0;
    nparticle = 0;

    while (line1 != "</LesHouchesEvents>") {

      fgets(&line_aux[0],300,fp1);
      sscanf(&line_aux[0], "%s",&str[0]);
      line1 = str;
      nline1++;
      if(line1.Contains("<event>")){
	nline_local1=1;
	continue;  //returns at beginning of the WHILE loop
      }
      else if(line1.BeginsWith("#")){
	continue;
      }
      else if(nline_local1==1){ // line contains common event information
	nline_local1++;
	continue;
      }
      else if(line1.Contains("</event>")){



        // Scan second file: always stop on "</event>"
        // line2 is cleared at the end in order to have a clean start
        while (1){
          fgets(&line_aux[0],300,fp2);
          sscanf(&line_aux[0], "%s",&str[0]);
          line2 = str;
          nline2++;
          if(line2.Contains("<event>")){
	    nline_local2=1;
	    continue;  //returns at beginning of the WHILE loop
          }
          else if(line2.BeginsWith("#")){
          continue;
          }
          else if(nline_local2==1){ // line contains common event information
	    nline_local2++;
      	continue;
          }
          else if(line2.Contains("</event>")){
	    break;  // stop loop
          }
          else{  // analyse line content
        	nline_local2++;
        	sscanf(&line_aux[0], "%*i %d %*i %*i %*i %*i %le %le %le %le %le %le %le",&imotherstatus,&dummy_double[0],&dummy_double[1],&dummy_double[2],&dummy_double[3],&dummy_double[4],&dummy_double[5],&dummy_double[6]);

        	if(TMath::Abs(imotherstatus)==1){
	      sscanf(&line_aux[0], "%i %d %*i %*i %*i %*i %le %le %le %le %le %le %le",&idup[nparticle],&istup[nparticle],&px[nparticle],&py[nparticle],&pz[nparticle],&E[nparticle],&dummy_double[0],&dummy_double[1],&dummy_double[2]);
 	      nparticle++;
	    }
        }  
      } // End while (line2 != "</event>")




/*
// Test
//  if(nevts<=10 || nevts>=49990){
   cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
   cout<<"file:"<<ki<<endl;
   cout<<"event"<<nevts<<endl;
  for(Int_t kj=0;kj<nparticle;kj++){
    cout<<idup[kj]<<" "<<istup[kj]<<" "<<px[kj]<<" "<<py[kj]<<" "<<pz[kj]<<" "<<E[kj]<<endl;
     }
//  }
// End Test
*/

    	    tree->Fill();
	    nparticle=0;
	    nevts++;
          if(nevts-(nevts/5000)*5000 == 0){
            cout<<nevts<<endl;
          }
	    nline_local1=0;
	    nline_local2=0;
          line2="";

      }
      else{  // analyse line content
	nline_local1++;
	sscanf(&line_aux[0], "%*i %d %*i %*i %*i %*i %le %le %le %le %le %le %le",&imotherstatus,&dummy_double[0],&dummy_double[1],&dummy_double[2],&dummy_double[3],&dummy_double[4],&dummy_double[5],&dummy_double[6]);

	if(TMath::Abs(imotherstatus)==1){
	  sscanf(&line_aux[0], "%i %d %*i %*i %*i %*i %le %le %le %le %le %le %le",&idup[nparticle],&istup[nparticle],&px[nparticle],&py[nparticle],&pz[nparticle],&E[nparticle],&dummy_double[0],&dummy_double[1],&dummy_double[2]);
	  nparticle++;
	}
      } // End while (line1 != "</LesHouchesEvents>")


    }  // end  while (line != "</LesHouchesEvents>") {


    fclose(fp1);
    fclose(fp2);

    cout<<endl<<ki+1<<"  >>>>>"<<gen1[ki].Data()<<"<<<<<"<<":"<<endl;
    cout<<nline1+ninitlines1<<" lines read in "<<LHEFfilename1.Data()<<endl;
    cout<<endl<<ki+1<<"  >>>>>"<<gen2[ki].Data()<<"<<<<<"<<":"<<endl;
    cout<<nline2+ninitlines2<<" lines read in "<<LHEFfilename2.Data()<<endl;
    cout<<nevts<<" events processed"<<endl;
    cout<<endl;
  
  } //end for(Int_t ki=0;ki<nlhef;ki++){



  rootfile->Write();

  cout<<endl;
  cout<<"Execution successful."<<endl;
  cout<<endl;


}
