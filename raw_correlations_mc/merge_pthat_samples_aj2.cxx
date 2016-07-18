/// reco PbPb
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>
#include <vector>

using namespace std;

#define nCBins 4
#define nPtBins 1
#define nTrkPtBins 6
#define nAjBins 2

float trkPtCut=1;

int parti = -999;
bool is_data = false;


enum enum_dataset_types {e_Data2011,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170, e_n_dataset_types};
int dataset_type_code = -999;

TString dataset_type_strs[e_n_dataset_types] = {"Data2011","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170"};
//Hydjet80 = 5
//Pythia80 = 14

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220};
TString dataset_pthat_strings[e_n_dataset_types+1] = {"pthat0", "pthat0","pthat15","pthat30","pthat50","pthat80","pthat120","pthat170","pthat220","pthat280","pthat370","pthat15","pthat30","pthat50","pthat80","pthat120","pthat170"};

double cross_section[e_n_dataset_types+1] = {0,0,2.034E-01,1.075E-02,1.025E-03,9.865E-05,1.129E-05,1.465E-06,2.837E-07 ,5.323E-08 ,5.934E-09,2.034E-01, 1.075E-02,1.025E-03,9.865E-05,1.129E-05,1.465E-06 };

double total_events[e_n_dataset_types] =  {0,  0,   333000 ,250000,   395000,   368000,   367000,   392000,    181000,     50000,   14000,100000, 100000,123000, 104500,97500,69000};


enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,e_n_data_mc_types};

TString data_mc_type_strs[e_n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};


float PtBins[nPtBins+1] = {100, 300};
TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

float CBins[nCBins+1] = {0, 20, 60, 100, 200};
TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};

float TrkPtBins[nTrkPtBins+1] = {0.5, 1, 2, 3, 4, 8, 300};
TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt300" };

float AjBins[nAjBins+1] = {0,0.22,0.75};
TString AjBin_strs[nAjBins+2] = {"Aj0","Aj22","Aj75","Dummy"};

double weight[e_n_dataset_types];
double inclusive_jets[e_n_data_mc_types][e_n_dataset_types][nCBins];
double dijets_lead[e_n_data_mc_types][e_n_dataset_types][nCBins];
double dijets_sub[e_n_data_mc_types][e_n_dataset_types][nCBins];
double dijets_lead_Aj[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins];
double dijets_sub_Aj[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins];
double nonleading_jets[e_n_data_mc_types][e_n_dataset_types][nCBins];

TH1F* vz[e_n_data_mc_types][e_n_dataset_types];
TH1F* vz_reweighted[e_n_data_mc_types][e_n_dataset_types];
TH1F* cent[e_n_data_mc_types][e_n_dataset_types];
TH1F* cent_reweighted[e_n_data_mc_types][e_n_dataset_types];

TH1F* all_jets_corrpT[e_n_data_mc_types][e_n_dataset_types][nCBins];
TH1F* all_jets_eta[e_n_data_mc_types][e_n_dataset_types][nCBins];
TH1F* all_jets_phi[e_n_data_mc_types][e_n_dataset_types][nCBins];

TH1F* only_leadingjets_corrpT[e_n_data_mc_types][e_n_dataset_types][nCBins];
TH1F* only_subleadingjets_corrpT[e_n_data_mc_types][e_n_dataset_types][nCBins];
TH1F* only_leadingjets_phi[e_n_data_mc_types][e_n_dataset_types][nCBins];
TH1F* only_subleadingjets_phi[e_n_data_mc_types][e_n_dataset_types][nCBins];
TH1F* only_leadingjets_eta[e_n_data_mc_types][e_n_dataset_types][nCBins];
TH1F* only_subleadingjets_eta[e_n_data_mc_types][e_n_dataset_types][nCBins];


TH1F* only_leadingjets_Aj_corrpT[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins];
TH1F* only_subleadingjets_Aj_corrpT[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins];
TH1F* only_leadingjets_Aj_phi[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins];
TH1F* only_subleadingjets_Aj_phi[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins];
TH1F* only_leadingjets_Aj_eta[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins];
TH1F* only_subleadingjets_Aj_eta[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins];


TH2D* hJetTrackMELeading[e_n_data_mc_types][e_n_dataset_types][nCBins][nTrkPtBins];
TH2D* hJetTrackMESubLeading[e_n_data_mc_types][e_n_dataset_types][nCBins][nTrkPtBins];
TH2D* hJetTrackMELeading_pTweighted[e_n_data_mc_types][e_n_dataset_types][nCBins][nTrkPtBins];
TH2D* hJetTrackMESubLeading_pTweighted[e_n_data_mc_types][e_n_dataset_types][nCBins][nTrkPtBins];

TH2D* hJetTrackMELeading_Aj[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins][nTrkPtBins];
TH2D* hJetTrackMESubLeading_Aj[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins][nTrkPtBins];
TH2D* hJetTrackMELeading_Aj_pTweighted[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins][nTrkPtBins];
TH2D* hJetTrackMESubLeading_Aj_pTweighted[e_n_data_mc_types][e_n_dataset_types][nCBins][nAjBins][nTrkPtBins];


TFile *in_file[e_n_data_mc_types][e_n_dataset_types];
TFile *in_file_jets[e_n_data_mc_types][e_n_dataset_types];

Int_t merge_pthat_samples_aj2(bool is_pythia=kFALSE, int first_mc_type = 1, int n_mc_types = 5){

  int min_file, max_file;

  if(is_pythia){
    min_file = 14;
    max_file = 17;
  }else{
    min_file = 5;
    max_file = 10;

  }

  int lbin, rbin;
 
  double summed_cross_section = cross_section[min_file]-cross_section[min_file+1];
  double summed_total_events = total_events[min_file];
  double summed_weight = (cross_section[min_file]-cross_section[min_file+1])/total_events[min_file];

  for(int file_i = min_file+1; file_i<max_file; file_i++){
    summed_cross_section += (cross_section[file_i]-cross_section[file_i+1]);
    summed_total_events += total_events[file_i];
    summed_weight +=  (cross_section[file_i]-cross_section[file_i+1])/total_events[file_i];
    cout<<summed_cross_section<<" "<<summed_total_events<<" "<<endl;
  }

  //First we get some nice merged spectra!

  for(int mc_type_code = first_mc_type; mc_type_code <n_mc_types; mc_type_code++){

    TString infile_stem = "unmerged_parts_reco_gen_me/";
  
    if(is_pythia){
      infile_stem = "unmerged_parts_reco_gen_me/";
    }

    int mc_type_code_jets = mc_type_code;

    //there are some special hist mc types that only differ in terms of track info--for these we need to refer back to the mc type that contains jet info

    if(mc_type_code==1||mc_type_code==5||mc_type_code==6) mc_type_code_jets = 2;
    else if(mc_type_code==3||mc_type_code==7||mc_type_code==8) mc_type_code_jets = 4;
    else if(mc_type_code>8&&mc_type_code%2==0) mc_type_code_jets = mc_type_code-1;

    for(int file_i = min_file; file_i<max_file; file_i++){
     
      in_file[mc_type_code][file_i] = new TFile((TString)(infile_stem+dataset_type_strs[file_i]+"_"+data_mc_type_strs[mc_type_code]+"_Merged.root"),"READ");
      in_file_jets[mc_type_code][file_i] = new TFile((TString)(infile_stem+dataset_type_strs[file_i]+"_"+data_mc_type_strs[mc_type_code_jets]+"_Merged.root"),"READ");
    
      cout<<in_file[mc_type_code][file_i]->GetName()<<endl;

      weight[file_i] = (cross_section[file_i]-cross_section[file_i+1])/total_events[file_i]/summed_weight;
      
      if(file_i==min_file)cout<<"this is min_file"<<endl;

      //weight[file_i] = 1.;

     

      cout<<"File "<<in_file[mc_type_code][file_i]->GetName()<<" Cross section: "<<cross_section[file_i]<<" N triggers: "<<total_events[file_i]<<" Weight: "<<weight[file_i]<<endl;
     
      if(mc_type_code == 2){
	vz[file_i][2]=(TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code] + "_Vz"))->Clone((TString)(data_mc_type_strs[mc_type_code] + "_Vz"));
	vz_reweighted[file_i][2] = (TH1F*) in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code] + "_Vz_Reweighted"))->Clone((TString)(data_mc_type_strs[mc_type_code] + "_Vz_Reweighted"));
	cent[file_i][2]= (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code] + "_Centrality"))->Clone((TString)(data_mc_type_strs[mc_type_code] + "_Centrality"));
	cent_reweighted[file_i][2]= (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code] + "_Centrality_Reweighted"))->Clone((TString)(data_mc_type_strs[mc_type_code] + "_Centrality_Reweighted"));

      }

      for (int ibin=0; ibin<nCBins; ibin ++){

	  all_jets_corrpT[mc_type_code][file_i][ibin] = (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets] + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);

	  all_jets_phi[mc_type_code][file_i][ibin] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_all_jets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_all_jets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);

	  all_jets_eta[mc_type_code][file_i][ibin] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_all_jets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_all_jets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);

	  //// leading jet histograms
	  only_leadingjets_corrpT[mc_type_code][file_i][ibin] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_leadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_leadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
	  only_leadingjets_phi[mc_type_code][file_i][ibin] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_leadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_leadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
	  only_leadingjets_eta[mc_type_code][file_i][ibin] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_leadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_leadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
	  //// subleading jet histograms
	  only_subleadingjets_corrpT[mc_type_code][file_i][ibin] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_subleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_subleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
	  only_subleadingjets_phi[mc_type_code][file_i][ibin] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_subleadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_subleadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
	  only_subleadingjets_eta[mc_type_code][file_i][ibin] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_subleadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_subleadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
	  

	  for(int ibin4 = 0; ibin4<nAjBins; ibin4++){


	    //// leading jet histograms
	    only_leadingjets_Aj_corrpT[mc_type_code][file_i][ibin][ibin4] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_leadingjets_Aj_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_leadingjets_Aj_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);
	    only_leadingjets_Aj_phi[mc_type_code][file_i][ibin][ibin4] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_leadingjets_Aj_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_leadingjets_Aj_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);
	    only_leadingjets_Aj_eta[mc_type_code][file_i][ibin][ibin4] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_leadingjets_Aj_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_leadingjets_Aj_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);

	    //// subleading jet histograms
	    only_subleadingjets_Aj_corrpT[mc_type_code][file_i][ibin][ibin4] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_subleadingjets_Aj_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_subleadingjets_Aj_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);
	    only_subleadingjets_Aj_phi[mc_type_code][file_i][ibin][ibin4] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_subleadingjets_Aj_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_subleadingjets_Aj_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);
	    only_subleadingjets_Aj_eta[mc_type_code][file_i][ibin][ibin4] =  (TH1F*)in_file_jets[mc_type_code][file_i]->Get((TString)(data_mc_type_strs[mc_type_code_jets]+ "_only_subleadingjets_Aj_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]))->Clone(data_mc_type_strs[mc_type_code]+"_"+dataset_type_strs[file_i] + "_only_subleadingjets_Aj_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);

	  }//Aj


      } // cent loop
    
    } //file loop

    if(mc_type_code == 2){
      vz[0][2]= (TH1F*)vz[min_file][2]->Clone((TString)(data_mc_type_strs[mc_type_code] + "_Vz_Merged"));
      vz_reweighted[0][2]= (TH1F*)vz_reweighted[min_file][2]->Clone((TString)(data_mc_type_strs[mc_type_code] + "_Vz_Reweighted_Merged"));
      cent[0][2]= (TH1F*)cent[min_file][2]->Clone((TString)(data_mc_type_strs[mc_type_code] + "_Centrality_Merged"));
      cent_reweighted[0][2]= (TH1F*)cent_reweighted[min_file][2]->Clone((TString)(data_mc_type_strs[mc_type_code] + "_Centrality_Reweighted_Merged"));
  
    }


    for (int ibin=0; ibin<nCBins; ibin ++){

      all_jets_corrpT[mc_type_code][0][ibin] = (TH1F*)all_jets_corrpT[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
      all_jets_phi[mc_type_code][0][ibin] = (TH1F*)all_jets_phi[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"all_jets_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
      all_jets_eta[mc_type_code][0][ibin] = (TH1F*)all_jets_eta[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"all_jets_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);

      only_leadingjets_corrpT[mc_type_code][0][ibin] = (TH1F*)only_leadingjets_corrpT[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"only_leadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
      only_leadingjets_phi[mc_type_code][0][ibin] = (TH1F*)only_leadingjets_phi[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"only_leadingjets_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
      only_leadingjets_eta[mc_type_code][0][ibin] = (TH1F*)only_leadingjets_eta[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"only_leadingjets_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);


      only_subleadingjets_corrpT[mc_type_code][0][ibin] = (TH1F*)only_subleadingjets_corrpT[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"only_subleadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
      only_subleadingjets_phi[mc_type_code][0][ibin] = (TH1F*)only_subleadingjets_phi[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"only_subleadingjets_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);
      only_subleadingjets_eta[mc_type_code][0][ibin] = (TH1F*)only_subleadingjets_eta[mc_type_code][min_file][ibin]->Clone(data_mc_type_strs[mc_type_code]+"only_subleadingjets_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]);

      all_jets_corrpT[mc_type_code][0][ibin]->Scale(weight[min_file]);
      all_jets_phi[mc_type_code][0][ibin]->Scale(weight[min_file]);
      all_jets_eta[mc_type_code][0][ibin]->Scale(weight[min_file]);

      only_leadingjets_corrpT[mc_type_code][0][ibin]->Scale(weight[min_file]);
      only_leadingjets_phi[mc_type_code][0][ibin]->Scale(weight[min_file]);
      only_leadingjets_eta[mc_type_code][0][ibin]->Scale(weight[min_file]);

      only_subleadingjets_corrpT[mc_type_code][0][ibin]->Scale(weight[min_file]);
      only_subleadingjets_phi[mc_type_code][0][ibin]->Scale(weight[min_file]);
      only_subleadingjets_eta[mc_type_code][0][ibin]->Scale(weight[min_file]);
 
      if(mc_type_code == 2 &&ibin ==0){
	vz[0][2]->Scale(weight[min_file]);
	vz_reweighted[0][2]->Scale(weight[min_file]);
	cent[0][2]->Scale(weight[min_file]);
	cent_reweighted[0][2]->Scale(weight[min_file]);
      }

      
      for(int file_i = min_file+1; file_i<max_file; file_i++){

	if(mc_type_code == 2 &&ibin ==0){
	  vz[0][2]->Add(vz[file_i][2],weight[file_i]);
	  vz_reweighted[0][2]->Add(vz_reweighted[file_i][2],weight[file_i]);
	  cent[0][2]->Add(cent[file_i][2],weight[file_i]);
	  cent_reweighted[0][2]->Add(cent_reweighted[file_i][2],weight[file_i]);
	    
	}
	all_jets_corrpT[mc_type_code][0][ibin]->Add(all_jets_corrpT[mc_type_code][file_i][ibin],weight[file_i]);
	all_jets_phi[mc_type_code][0][ibin]->Add(all_jets_phi[mc_type_code][file_i][ibin],weight[file_i]);
	all_jets_eta[mc_type_code][0][ibin]->Add(all_jets_eta[mc_type_code][file_i][ibin],weight[file_i]);

	only_leadingjets_corrpT[mc_type_code][0][ibin]->Add(only_leadingjets_corrpT[mc_type_code][file_i][ibin],weight[file_i]);
	only_leadingjets_phi[mc_type_code][0][ibin]->Add(only_leadingjets_phi[mc_type_code][file_i][ibin],weight[file_i]);
	only_leadingjets_eta[mc_type_code][0][ibin]->Add(only_leadingjets_eta[mc_type_code][file_i][ibin],weight[file_i]);

	only_subleadingjets_corrpT[mc_type_code][0][ibin]->Add(only_subleadingjets_corrpT[mc_type_code][file_i][ibin],weight[file_i]);
	only_subleadingjets_phi[mc_type_code][0][ibin]->Add(only_subleadingjets_phi[mc_type_code][file_i][ibin],weight[file_i]);
	only_subleadingjets_eta[mc_type_code][0][ibin]->Add(only_subleadingjets_eta[mc_type_code][file_i][ibin],weight[file_i]);

      }


      for(int ibin4 = 0; ibin4<nAjBins; ibin4++){

	only_leadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4] = (TH1F*)only_leadingjets_Aj_corrpT[mc_type_code][min_file][ibin][ibin4]->Clone(data_mc_type_strs[mc_type_code]+"only_leadingjets_Aj_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);
	only_leadingjets_Aj_phi[mc_type_code][0][ibin][ibin4] = (TH1F*)only_leadingjets_Aj_phi[mc_type_code][min_file][ibin][ibin4]->Clone(data_mc_type_strs[mc_type_code]+"only_leadingjets_Aj_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);
	only_leadingjets_Aj_eta[mc_type_code][0][ibin][ibin4] = (TH1F*)only_leadingjets_Aj_eta[mc_type_code][min_file][ibin][ibin4]->Clone(data_mc_type_strs[mc_type_code]+"only_leadingjets_Aj_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);


	only_subleadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4] = (TH1F*)only_subleadingjets_Aj_corrpT[mc_type_code][min_file][ibin][ibin4]->Clone(data_mc_type_strs[mc_type_code]+"only_subleadingjets_Aj_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);
	only_subleadingjets_Aj_phi[mc_type_code][0][ibin][ibin4] = (TH1F*)only_subleadingjets_Aj_phi[mc_type_code][min_file][ibin][ibin4]->Clone(data_mc_type_strs[mc_type_code]+"only_subleadingjets_Aj_phi_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);
	only_subleadingjets_Aj_eta[mc_type_code][0][ibin][ibin4] = (TH1F*)only_subleadingjets_Aj_eta[mc_type_code][min_file][ibin][ibin4]->Clone(data_mc_type_strs[mc_type_code]+"only_subleadingjets_Aj_eta_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]);

	only_leadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4]->Scale(weight[min_file]);
	only_leadingjets_Aj_phi[mc_type_code][0][ibin][ibin4]->Scale(weight[min_file]);
	only_leadingjets_Aj_eta[mc_type_code][0][ibin][ibin4]->Scale(weight[min_file]);

	only_subleadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4]->Scale(weight[min_file]);
	only_subleadingjets_Aj_phi[mc_type_code][0][ibin][ibin4]->Scale(weight[min_file]);
	only_subleadingjets_Aj_eta[mc_type_code][0][ibin][ibin4]->Scale(weight[min_file]);
 
      
	for(int file_i = min_file+1; file_i<max_file; file_i++){

	  only_leadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4]->Add(only_leadingjets_Aj_corrpT[mc_type_code][file_i][ibin][ibin4],weight[file_i]);
	  only_leadingjets_Aj_phi[mc_type_code][0][ibin][ibin4]->Add(only_leadingjets_Aj_phi[mc_type_code][file_i][ibin][ibin4],weight[file_i]);
	  only_leadingjets_Aj_eta[mc_type_code][0][ibin][ibin4]->Add(only_leadingjets_Aj_eta[mc_type_code][file_i][ibin][ibin4],weight[file_i]);

	  only_subleadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4]->Add(only_subleadingjets_Aj_corrpT[mc_type_code][file_i][ibin][ibin4],weight[file_i]);
	  only_subleadingjets_Aj_phi[mc_type_code][0][ibin][ibin4]->Add(only_subleadingjets_Aj_phi[mc_type_code][file_i][ibin][ibin4],weight[file_i]);
	  only_subleadingjets_Aj_eta[mc_type_code][0][ibin][ibin4]->Add(only_subleadingjets_Aj_eta[mc_type_code][file_i][ibin][ibin4],weight[file_i]);


	}
      }
      
    }
 
  }


  cout<<" NOW WE HAVE MERGED SPECTRA AND WE'RE READY FOR CORRELATIONS! "<<endl;

  for(int mc_type_code = first_mc_type; mc_type_code <n_mc_types; mc_type_code++){

    for(int file_i = min_file; file_i<max_file; file_i++){
    
      for (int ibin=0; ibin<nCBins; ibin ++){

	if(is_pythia&&ibin>0)continue;

	for(int ibin3 = 0; ibin3< nTrkPtBins; ibin3++){

	  if(mc_type_code==1||mc_type_code==3){

	    hJetTrackMELeading[mc_type_code][file_i][ibin][ibin3]= (TH2D*)in_file[mc_type_code][file_i]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	    hJetTrackMESubLeading[mc_type_code][file_i][ibin][ibin3]= (TH2D*)in_file[mc_type_code][file_i]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  

	  }else{

	    hJetTrackMELeading[mc_type_code][file_i][ibin][ibin3]= (TH2D*)in_file[mc_type_code][file_i]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	    hJetTrackMESubLeading[mc_type_code][file_i][ibin][ibin3]= (TH2D*)in_file[mc_type_code][file_i]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));


	  }

	  for(int ibin4=0; ibin4<nAjBins; ibin4++){

	    if(mc_type_code==1||mc_type_code==3){
	      hJetTrackMELeading_Aj[mc_type_code][file_i][ibin][ibin4][ibin3]= (TH2D*)in_file[mc_type_code][file_i]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	      hJetTrackMESubLeading_Aj[mc_type_code][file_i][ibin][ibin4][ibin3]= (TH2D*)in_file[mc_type_code][file_i]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	    }else{

	 
	      hJetTrackMELeading_Aj[mc_type_code][file_i][ibin][ibin4][ibin3]= (TH2D*)in_file[mc_type_code][file_i]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Aj_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	      hJetTrackMESubLeading_Aj[mc_type_code][file_i][ibin][ibin4][ibin3]= (TH2D*)in_file[mc_type_code][file_i]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Aj_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_Pt100_Pt300"+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	    }



	  }//Aj 

	}//TrkPt loop
      
      } // cent loop
    
    } //file loop
    



    cout<<"got all correlations"<<endl;

    for (int ibin=0; ibin<nCBins; ibin ++){

      if(is_pythia&&ibin>0)continue;
  
      for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	hJetTrackMELeading[mc_type_code][0][ibin][ibin3]= (TH2D*)hJetTrackMELeading[mc_type_code][min_file][ibin][ibin3]->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	hJetTrackMESubLeading[mc_type_code][0][ibin][ibin3]= (TH2D*)hJetTrackMESubLeading[mc_type_code][min_file][ibin][ibin3]->Clone((TString)(data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	hJetTrackMELeading[mc_type_code][0][ibin][ibin3]->Scale(1.*weight[min_file]);
	hJetTrackMESubLeading[mc_type_code][0][ibin][ibin3]->Scale(1.*weight[min_file]);


      }
    
 
      for(int file_i = min_file+1; file_i<max_file; file_i++){
	
	for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	  hJetTrackMELeading[mc_type_code][0][ibin][ibin3]->Add(hJetTrackMELeading[mc_type_code][file_i][ibin][ibin3], weight[file_i]);
	  hJetTrackMESubLeading[mc_type_code][0][ibin][ibin3]->Add(hJetTrackMESubLeading[mc_type_code][file_i][ibin][ibin3], weight[file_i]);


	
	}
      }


      cout<<"ready to normalize"<<endl;

      inclusive_jets[mc_type_code][0][ibin] = all_jets_corrpT[mc_type_code][0][ibin]->Integral();
      dijets_lead[mc_type_code][0][ibin] = only_leadingjets_corrpT[mc_type_code][0][ibin]->Integral();
      dijets_sub[mc_type_code][0][ibin] = only_subleadingjets_corrpT[mc_type_code][0][ibin]->Integral();

      cout<<"Njets for "<<mc_type_code<<": "<<inclusive_jets[mc_type_code][0][ibin]<<endl;
      cout<<"Leading jets for "<<mc_type_code<<": "<<dijets_lead[mc_type_code][0][ibin]<<endl;
      cout<<"SubLeading jets for "<<mc_type_code<<": "<<dijets_sub[mc_type_code][0][ibin]<<endl;
      
      //      dijets_sub[mc_type_code][0][ibin] = dijets_lead[mc_type_code][0][ibin];

      for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	hJetTrackMELeading[mc_type_code][0][ibin][ibin3]->Scale(1./dijets_lead[mc_type_code][0][ibin]);
	hJetTrackMESubLeading[mc_type_code][0][ibin][ibin3]->Scale(1./dijets_sub[mc_type_code][0][ibin]);
     
      }

   

      for(int ibin4 = 0; ibin4<nAjBins; ibin4++){

	for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	  hJetTrackMELeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]= (TH2D*)hJetTrackMELeading_Aj[mc_type_code][min_file][ibin][ibin4][ibin3]->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Aj_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	  hJetTrackMESubLeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]= (TH2D*)hJetTrackMESubLeading_Aj[mc_type_code][min_file][ibin][ibin4][ibin3]->Clone((TString)(data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Aj_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	
	  hJetTrackMELeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]->Scale(weight[min_file]);
	  hJetTrackMESubLeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]->Scale(weight[min_file]);

	}
    

	for(int file_i = min_file+1; file_i<max_file; file_i++){
	
	  for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	    hJetTrackMELeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]->Add(hJetTrackMELeading_Aj[mc_type_code][file_i][ibin][ibin4][ibin3], weight[file_i]);
	    hJetTrackMESubLeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]->Add(hJetTrackMESubLeading_Aj[mc_type_code][file_i][ibin][ibin4][ibin3], weight[file_i]);

	
	  }
	}
      }

      cout<<"ready to normalize"<<endl;
	//-------------------------
	// Combine Aj Bins & Normalize!
	//------------------------
   
   
      for(int ibin4 = 0; ibin4<nAjBins; ibin4++){


	dijets_lead_Aj[mc_type_code][0][ibin][ibin4] = only_leadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4]->Integral();
	dijets_sub_Aj[mc_type_code][0][ibin][ibin4] = only_subleadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4]->Integral();


	cout<<"Leading_Aj jets for "<<mc_type_code<<": "<<dijets_lead_Aj[mc_type_code][0][ibin][ibin4]<<endl;
	cout<<"SubLeading_Aj jets for "<<mc_type_code<<": "<<dijets_sub_Aj[mc_type_code][0][ibin][ibin4]<<endl;
      
	//      dijets_sub[mc_type_code][0][ibin] = dijets_lead[mc_type_code][0][ibin];

	for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	  hJetTrackMELeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]->Scale(1./dijets_lead_Aj[mc_type_code][0][ibin][ibin4]);
	  hJetTrackMESubLeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]->Scale(1./dijets_sub_Aj[mc_type_code][0][ibin][ibin4]);

	 
	}

      }//Aj
    }
  }
    cout<<"ready to write"<<endl;

  
  for(int mc_type_code = first_mc_type; mc_type_code<n_mc_types; mc_type_code ++){

    TFile *out_file;
    if(is_pythia){
      out_file = new TFile("Pythia_RecoGen_ME_Merged_"+data_mc_type_strs[mc_type_code]+".root","RECREATE"); 
    }else{
      out_file = new TFile("HydJet_RecoGen_ME_Merged_"+data_mc_type_strs[mc_type_code]+".root","RECREATE");
    }

    cout<<"here"<<endl;
    if(mc_type_code%2 ==0){

      if(mc_type_code == 2 ){
	  vz[0][2]->Write();
	  vz_reweighted[0][2]->Write();
	  cent[0][2]->Write();
	  cent_reweighted[0][2]->Write();
	}
    
	for (int ibin=0; ibin<nCBins; ibin ++){

	  all_jets_corrpT[mc_type_code][0][ibin]->Write();
	  all_jets_phi[mc_type_code][0][ibin]->Write();
	  all_jets_eta[mc_type_code][0][ibin]->Write();

	  only_leadingjets_corrpT[mc_type_code][0][ibin]->Write();
	  only_leadingjets_phi[mc_type_code][0][ibin]->Write();
	  only_leadingjets_eta[mc_type_code][0][ibin]->Write();

	  only_subleadingjets_corrpT[mc_type_code][0][ibin]->Write();
	  only_subleadingjets_phi[mc_type_code][0][ibin]->Write();
	  only_subleadingjets_eta[mc_type_code][0][ibin]->Write();

	  for(int ibin4 = 0; ibin4<nAjBins; ibin4++){
	    
	    only_leadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4]->Write();
	    only_leadingjets_Aj_phi[mc_type_code][0][ibin][ibin4]->Write();
	    only_leadingjets_Aj_eta[mc_type_code][0][ibin][ibin4]->Write();

	    only_subleadingjets_Aj_corrpT[mc_type_code][0][ibin][ibin4]->Write();
	    only_subleadingjets_Aj_phi[mc_type_code][0][ibin][ibin4]->Write();
	    only_subleadingjets_Aj_eta[mc_type_code][0][ibin][ibin4]->Write();

	  }
	}

      }

    for(int ibin = 0; ibin<nCBins; ibin++){

	if(is_pythia&&ibin>0)continue;

	for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){
	  hJetTrackMELeading[mc_type_code][0][ibin][ibin3]->Write();
	  hJetTrackMESubLeading[mc_type_code][0][ibin][ibin3]->Write();
 
	

	  for(int ibin4 = 0; ibin4<nAjBins; ibin4++){

	    hJetTrackMELeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]->Write();
	    hJetTrackMESubLeading_Aj[mc_type_code][0][ibin][ibin4][ibin3]->Write();
 

	  }
	}


      }
     cout<<"wrote MC type "<<mc_type_code<<endl;
  }
  return 0;
}
