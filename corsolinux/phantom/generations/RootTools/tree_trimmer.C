/*

  TO BE FINISHED
  
  
  
 *tree_trimmer.C
 *
 *Last update: Apr 16, 2008
 *
 *This routine selects events in a tree.root file according to a set of cuts
 *and writes them down into a new root file
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
 * NOTE: It maight be necessary in the future to store also the "mothers"
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
  cout<<"                    -- tree_trimmer --                        "<<endl;
  cout<<"**************************************************************"<<endl;
  cout<<endl;

  TString treefilename;
  cout<<"Enter name of tree file "<<endl;
  cin>>treefilename;


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



  //Open the tree file and get the tree
  TFile treefile(treefilename.Data());
  TTree *tree_in = (TTree*) treefile.Get("events");


  // Take info from tree_in to fill histograms_____________________________________
  tree_in->SetBranchAddress("px",px);
  tree_in->SetBranchAddress("py",py);
  tree_in->SetBranchAddress("pz",pz);
  tree_in->SetBranchAddress("E",E);
  tree_in->SetBranchAddress("IDUP",idup);
  tree_in->SetBranchAddress("ISTUP",istup);
  // Take the value of total cross section (stored inside the tree as User
  // Information)
  TObjString *os_xsec = (TObjString*)tree_in->GetUserInfo()->At(0);
  TString s_xsec = os_xsec->GetString();
  xsec = s_xsec.Atof();
  xsec_tot += xsec;
  
  cout<<"Read total cross section: "<<xsec<<" pb"<<endl;
  nentries_tree = (Int_t)tree->GetEntries();



  // Start analysis_____________________________________________________________

  // Loop over events stored inside the tree
  for(Int_t k_n=0;k_n<nentries_tree;k_n++) {

    tree->GetEntry(k_n);

    for(Int_t k=0;k<8;k++){

      //Make sure that the first two particles are ingoing (initial state)
      if(istup[0] != -1 || istup[1] != -1){
	cout<<"***ERROR: the first 2 particles are NOT ingoing!!!"<<endl;
	cout<<"          Please check the event file."<<endl<<endl;
	cout<<"Execution stopped"<<endl;
	exit(0);
      }

      particle[k] = TLorentzVector(px[k],py[k],pz[k],E[k]);

    }





// TILL HERE








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
