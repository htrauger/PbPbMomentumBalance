#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLatex.h"


#include <iostream>
#include <vector>
#include <fstream>


Int_t me_correct_aj(bool is_pp = kFALSE){

  //-------------Input files & set eta-----------------

  TString datalabel;

  TLatex *centtex, *pttex,*ajtex;

  TF1 *pol0 = new TF1("pol0","[0]+x-x",-.3,.3);

  TLegend *l;

  TFile *fin, *fin_mb, *fout;

  float norm_temp, width_temp, width_temp_x, width_temp_y, width_temp_ref, max_bin, max_cont,bc,err, me00_range, me00_sub, me00_lead, me00_err;
    
  TString me00_range_string;
  

  if(is_pp){
    fin = new TFile("../raw_correlations/Data_pp_Aj_v5_noweight.root","READ");
    fout = new TFile((TString)("pp_Correlations.root"),"RECREATE");
       
    me00_range = 0.2;
    datalabel = "pp";
    
  }else{
    //  fin = new TFile("../raw_correlations/Data_PbPb_Aj_v5_jettrig.root","READ");
    fin = new TFile("../raw_correlations/Data_PbPb_Aj_v5_mb.root","READ");
    fin_mb = new TFile("../raw_correlations/Data_PbPb_Aj_v5_mb.root","READ");
    fout = new TFile((TString)("PbPb_Correlations.root"),"RECREATE");
       
    me00_range = 0.2;
    datalabel = "PbPb";
  }

  //----------------------------------------------------

  TCanvas *me_proj_leading_canvas = new TCanvas("me_proj_leading_canvas","",0,0,1500,2400);
  me_proj_leading_canvas->Divide(4,6,0.0000,0.0000);

  TCanvas *me_proj_sub_canvas = new TCanvas("me_proj_sub_canvas","",0,0,1500,2400);
  me_proj_sub_canvas->Divide(4,6,0.0000,0.0000);

  TCanvas *result_proj_canvas = new TCanvas("result_proj_canvas","",0,0,1500,2400);
  result_proj_canvas->Divide(4,6,0.0001,0.0001);

  TCanvas *result_proj_canvas_Aj0_Aj22 = new TCanvas("result_proj_canvas_Aj0_Aj22","",0,0,1500,2400);
  result_proj_canvas_Aj0_Aj22->Divide(4,6,0.0001,0.001);

  TCanvas *result_proj_canvas_Aj22_Aj75 = new TCanvas("result_proj_canvas_Aj22_Aj75","",0,0,1500,2400);
  result_proj_canvas_Aj22_Aj75->Divide(4,6,0.0001,0.0001);




  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 6;
  const int nAjBins = 4;

  float me00;
  int me00binl, me00binr, n_me00bins;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0, 10, 30, 50, 100};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 0-30%","Cent. 30-50%","Cent. 50-100%"};

 
  float TrkPtBins[nTrkPtBins+1] = {0.5,1, 2, 3, 4, 8, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};


  float AjBins[nAjBins+1] = {0,0.11,0.22,0.33,0.75};
  TString AjBin_strs[nAjBins+2] = {"Aj0","Aj11","Aj22","Aj33","Aj75","Dummy"};
 
   
  TString ref_name;
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  
  TH1F* only_leadingjets_corrpT[nCBins][nPtBins];
  TH1F* only_subleadingjets_corrpT[nCBins][nPtBins];
   
  TH2D* hJetTrackSignalBackgroundLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundSubLeading[nCBins][nPtBins][nTrkPtBins];
 
  TH2D* hJetTrackSignalBackgroundLeading_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundSubLeading_pTweighted[nCBins][nPtBins][nTrkPtBins];

  TH2D* hJetTrackMELeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackMESubLeading[nCBins][nPtBins][nTrkPtBins];

  TH2D *yield_lead[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_lead_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_sub[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_sub_pTweighted[nCBins][nPtBins][nTrkPtBins];

  TH1D *yield_lead_proj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH1D *yield_sub_proj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  
  TH1D *me_proj_lead[nCBins][nPtBins][nTrkPtBins];
  TH1D *me_proj_sub[nCBins][nPtBins][nTrkPtBins];
  
  TH1F* only_leadingjets_Aj_corrpT[nCBins][nPtBins][nAjBins];
  TH1F* only_subleadingjets_Aj_corrpT[nCBins][nPtBins][nAjBins];
   
  TH2D* hJetTrackSignalBackgroundLeading_Aj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundSubLeading_Aj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundLeading_Aj_pTweighted[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[nCBins][nPtBins][nAjBins][nTrkPtBins];
    
  TH2D* hJetTrackMELeading_Aj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH2D* hJetTrackMESubLeading_Aj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  
  TH2D *yield_lead_Aj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH2D *yield_lead_Aj_pTweighted[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH2D *yield_sub_Aj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH2D *yield_sub_Aj_pTweighted[nCBins][nPtBins][nAjBins][nTrkPtBins];
  
  TH1D *me_proj_lead_Aj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  TH1D *me_proj_sub_Aj[nCBins][nPtBins][nAjBins][nTrkPtBins];
  
  TString desc = "Data";

  TCanvas *dummy = new TCanvas("dummy");

  //-----------------------
  // Start getting histos
  //-----------------------

  for (int ibin=0;ibin<nCBins;ibin++){
    
 
    dummy->cd();

    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 


      only_leadingjets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_only_leadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("only_leadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
 
      only_subleadingjets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_only_subleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("only_subleadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
 
      if(ibin==1){
	only_leadingjets_corrpT[ibin][ibin2]->Add(	only_leadingjets_corrpT[ibin-1][ibin2]);
 	only_subleadingjets_corrpT[ibin][ibin2]->Add(only_subleadingjets_corrpT[ibin-1][ibin2]);

      }

      fout->cd();
      only_leadingjets_corrpT[ibin][ibin2]->Write();
      only_subleadingjets_corrpT[ibin][ibin2]->Write();
             
      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){


	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	
	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundSubLeading"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	

	hJetTrackSignalBackgroundLeading_pTweighted[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundLeading_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	

	hJetTrackSignalBackgroundSubLeading_pTweighted[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundSubLeading_pTweighted"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  hJetTrackMESubLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackMESubLeading"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	  if(ibin==1){
	    hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Add(hJetTrackSignalBackgroundLeading[ibin-1][ibin2][ibin3]);
	    hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Add(hJetTrackSignalBackgroundSubLeading[ibin-1][ibin2][ibin3]);
	    hJetTrackSignalBackgroundLeading_pTweighted[ibin][ibin2][ibin3]->Add(hJetTrackSignalBackgroundLeading_pTweighted[ibin-1][ibin2][ibin3]);
	    hJetTrackSignalBackgroundSubLeading_pTweighted[ibin][ibin2][ibin3]->Add(hJetTrackSignalBackgroundSubLeading_pTweighted[ibin-1][ibin2][ibin3]);
	    hJetTrackMELeading[ibin][ibin2][ibin3]->Add(hJetTrackMELeading[ibin-1][ibin2][ibin3]);
	    hJetTrackMESubLeading[ibin][ibin2][ibin3]->Add( hJetTrackMESubLeading[ibin-1][ibin2][ibin3]);
	  }
      }

      for(int ibin4 = 0; ibin4 <nAjBins; ibin4++){

	only_leadingjets_Aj_corrpT[ibin][ibin2][ibin4] = (TH1F*)fin->Get((TString) (desc + "_only_leadingjets_Aj_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]))->Clone((TString) ("only_leadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]));
 
	only_subleadingjets_Aj_corrpT[ibin][ibin2][ibin4] = (TH1F*)fin->Get((TString) (desc + "_only_subleadingjets_Aj_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]))->Clone((TString) ("only_subleadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]));
 

	if(ibin==1){


	  only_leadingjets_Aj_corrpT[ibin][ibin2][ibin4]->Add(only_leadingjets_Aj_corrpT[ibin-1][ibin2][ibin4]);
	  only_subleadingjets_Aj_corrpT[ibin][ibin2][ibin4]->Add(only_subleadingjets_Aj_corrpT[ibin-1][ibin2][ibin4]);

	}

	for(int ibin3 = 0; ibin3 < nTrkPtBins; ibin3++){

	  hJetTrackSignalBackgroundLeading_Aj[ibin][ibin2][ibin4][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundLeading_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	  //

	  hJetTrackSignalBackgroundSubLeading_Aj[ibin][ibin2][ibin4][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundSubLeading_Aj"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  


	  hJetTrackSignalBackgroundLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundLeading_Aj_pTweighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	  

	  hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundSubLeading_Aj_pTweighted"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_pTweighted_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	  
	  hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackMELeading_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	  hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackMESubLeading_Aj"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
	  
	  if(ibin==1){

	    hJetTrackSignalBackgroundLeading_Aj[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackSignalBackgroundLeading_Aj[ibin-1][ibin2][ibin4][ibin3]);
	    hJetTrackSignalBackgroundSubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackSignalBackgroundSubLeading_Aj[ibin-1][ibin2][ibin4][ibin3]);
	    hJetTrackSignalBackgroundLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackSignalBackgroundLeading_Aj_pTweighted[ibin-1][ibin2][ibin4][ibin3]);
	
	    hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[ibin-1][ibin2][ibin4][ibin3]);
	    hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->Add( hJetTrackMELeading_Aj[ibin-1][ibin2][ibin4][ibin3]);
	    hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackMESubLeading_Aj[ibin-1][ibin2][ibin4][ibin3]);
	

	  }



	} //ibin4


      } /// ibin3

    } // ibin2

  }//ibin

  //-------------------------
  // Normalize all!
  //------------------------
  
  for(int ibin =0; ibin<nCBins; ibin++){ 
    for (int ibin2=0;ibin2<nPtBins;ibin2++){

      for(int ibin4 = 0; ibin4<nAjBins; ibin4++){
	if(ibin4==1||ibin4==3){ 
	  continue; 
	}else{
	  only_leadingjets_Aj_corrpT[ibin][ibin2][ibin4]->Add(only_leadingjets_Aj_corrpT[ibin][ibin2][ibin4+1]);
	  only_subleadingjets_Aj_corrpT[ibin][ibin2][ibin4]->Add(only_subleadingjets_Aj_corrpT[ibin][ibin2][ibin4+1]);

	  only_leadingjets_Aj_corrpT[ibin][ibin2][ibin4]->Write();
	  only_subleadingjets_Aj_corrpT[ibin][ibin2][ibin4]->Write();
	}

      }
      for(int ibin3 = 0; ibin3<nTrkPtBins; ibin3++){

	width_temp_x = 1.;
	width_temp_y = 1.;
  
	norm_temp = only_leadingjets_corrpT[ibin][ibin2]->GetEntries();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackSignalBackgroundLeading_pTweighted[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);

	norm_temp = only_subleadingjets_corrpT[ibin][ibin2]->GetEntries();
	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackSignalBackgroundSubLeading_pTweighted[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackMESubLeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);

	for(int ibin4 = 0; ibin4< nAjBins; ibin4++){

	  if(ibin4==1||ibin4==3){ 
	    continue; 
	  }else{
	    hJetTrackSignalBackgroundLeading_Aj[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackSignalBackgroundLeading_Aj[ibin][ibin2][ibin4+1][ibin3]);
	    hJetTrackSignalBackgroundLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackSignalBackgroundLeading_Aj_pTweighted[ibin][ibin2][ibin4+1][ibin3]);
	    hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackMELeading_Aj[ibin][ibin2][ibin4+1][ibin3]);

	    hJetTrackSignalBackgroundSubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackSignalBackgroundSubLeading_Aj[ibin][ibin2][ibin4+1][ibin3]);
	    hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[ibin][ibin2][ibin4+1][ibin3]);
	    hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Add(hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4+1][ibin3]);
	 
	  }
	  


	  norm_temp = only_leadingjets_Aj_corrpT[ibin][ibin2][ibin4]->GetEntries();
	  hJetTrackSignalBackgroundLeading_Aj[ibin][ibin2][ibin4][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	  hJetTrackSignalBackgroundLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	  hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);
   
	  norm_temp = only_subleadingjets_Aj_corrpT[ibin][ibin2][ibin4]->GetEntries();
	  hJetTrackSignalBackgroundSubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	  hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	  hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);

	}
      
	//---------------------------------------
	//    ME correct
	//-------------------------------------

	//  LEADING



	yield_lead[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_Leading_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	yield_lead_pTweighted[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading_pTweighted[ibin][ibin2][ibin3]->Clone((TString)("Yield_pTweighted_Leading_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));


	me_proj_lead[ibin][ibin2][ibin3] = (TH1D*)hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_lead_temp_%d%d%d",ibin,ibin2,ibin3),1,200);

	me_proj_lead[ibin][ibin2][ibin3]->Scale (1./200);

	me_proj_lead[ibin][ibin2][ibin3] -> Fit("pol0","q","",-me00_range,me00_range);

	me00_lead = me_proj_lead[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);

	//   SUBLEADING

	me_proj_sub_canvas->cd(4*(ibin3+1)-ibin);

	yield_sub[ibin][ibin2][ibin3] = (TH2D*) hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_SubLeading_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	yield_sub_pTweighted[ibin][ibin2][ibin3] = (TH2D*) hJetTrackSignalBackgroundSubLeading_pTweighted[ibin][ibin2][ibin3]->Clone((TString)("Yield_pTweighted_SubLeading_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	me_proj_sub[ibin][ibin2][ibin3] = (TH1D*)hJetTrackMESubLeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_sub_%d%d%d",ibin,ibin2,ibin3),1,200);

	me_proj_sub[ibin][ibin2][ibin3]->Scale (1./200);
	me_proj_sub[ibin][ibin2][ibin3] -> Fit("pol0","q","",-me00_range,me00_range);

	me00_sub = me_proj_sub[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);

	me00 = (me00_sub+me00_lead)/2.;

	me00_err = TMath::Abs(me00_sub-me00_lead)/me00;

	me_proj_sub[ibin][ibin2][ibin3]->SetAxisRange(-.299,.299);
	me_proj_sub[ibin][ibin2][ibin3]->Scale(1./me00);
	me_proj_sub[ibin][ibin2][ibin3]->Draw();


	if(ibin<3){
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	if(ibin==3){
	  TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}


	
	me_proj_leading_canvas->cd(4*(ibin3+1)-ibin);

	me_proj_lead[ibin][ibin2][ibin3]->Scale(1./me00);
	me_proj_lead[ibin][ibin2][ibin3]->SetAxisRange(-.299,.299);

	me_proj_lead[ibin][ibin2][ibin3]->Draw();



	if(ibin<3){
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	if(ibin==3){
	  TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}

	dummy->cd();


	if(ibin3==0){
	  for(int k = 1; k<hJetTrackMELeading[ibin][ibin2][ibin3]->GetNbinsX()/2.+1; k++){
	  
	    bc = (me_proj_lead[ibin][ibin2][ibin3]->GetBinContent(k)+me_proj_lead[ibin][ibin2][ibin3]->GetBinContent(me_proj_lead[ibin][ibin2][ibin3]->GetNbinsX()+1-k))/2.;
	    err = TMath::Sqrt(me_proj_lead[ibin][ibin2][ibin3]->GetBinError(k)*me_proj_lead[ibin][ibin2][ibin3]->GetBinError(k)+me_proj_lead[ibin][ibin2][ibin3]->GetBinError(me_proj_lead[ibin][ibin2][ibin3]->GetNbinsX()+1-k)*me_proj_lead[ibin][ibin2][ibin3]->GetBinError(me_proj_lead[ibin][ibin2][ibin3]->GetNbinsX()+1-k))/2.;
	  
	    if(TMath::Abs(hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(k))<0.3)bc = 1.;
	    

	    for(int m = 1; m<hJetTrackMELeading[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){

	      hJetTrackMELeading[ibin][ibin2][ibin3]->SetBinContent(k,m,bc);
	      hJetTrackMELeading[ibin][ibin2][ibin3]->SetBinError(k,m,err);
	      hJetTrackMELeading[ibin][ibin2][ibin3]->SetBinContent(me_proj_lead[ibin][ibin2][ibin3]->GetNbinsX()+1-k,m,bc);
	      hJetTrackMELeading[ibin][ibin2][ibin3]->SetBinError(me_proj_lead[ibin][ibin2][ibin3]->GetNbinsX()+1-k,m,err);

	    }
	  }

	  for(int k = 1; k<hJetTrackMESubLeading[ibin][ibin2][ibin3]->GetNbinsX()/2.+1; k++){
	  
	    bc = (me_proj_sub[ibin][ibin2][ibin3]->GetBinContent(k)+me_proj_sub[ibin][ibin2][ibin3]->GetBinContent(me_proj_sub[ibin][ibin2][ibin3]->GetNbinsX()+1-k))/2.;
	    err = TMath::Sqrt(me_proj_sub[ibin][ibin2][ibin3]->GetBinError(k)*me_proj_sub[ibin][ibin2][ibin3]->GetBinError(k)+me_proj_sub[ibin][ibin2][ibin3]->GetBinError(me_proj_sub[ibin][ibin2][ibin3]->GetNbinsX()+1-k)*me_proj_sub[ibin][ibin2][ibin3]->GetBinError(me_proj_sub[ibin][ibin2][ibin3]->GetNbinsX()+1-k))/2.;
	  
	    if(TMath::Abs(hJetTrackMESubLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(k))<0.3)bc = 1.;

	    for(int m = 1; m<hJetTrackMESubLeading[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){

	      hJetTrackMESubLeading[ibin][ibin2][ibin3]->SetBinContent(k,m,bc);
	      hJetTrackMESubLeading[ibin][ibin2][ibin3]->SetBinError(k,m,err);
	      hJetTrackMESubLeading[ibin][ibin2][ibin3]->SetBinContent(me_proj_lead[ibin][ibin2][ibin3]->GetNbinsX()+1-k,m,bc);
	      hJetTrackMESubLeading[ibin][ibin2][ibin3]->SetBinError(me_proj_lead[ibin][ibin2][ibin3]->GetNbinsX()+1-k,m,err);

	    }
	  }
	
	}else{
	  for(int k = 1; k<hJetTrackMELeading[ibin][ibin2][ibin3]->GetNbinsX()+1; k++){
	  
	    bc = me_proj_lead[ibin][ibin2][ibin3]->GetBinContent(k);
	    err = me_proj_lead[ibin][ibin2][ibin3]->GetBinError(k);

	    if(TMath::Abs(hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(k))<0.3)bc = 1.;	  

	    for(int m = 1; m<hJetTrackMELeading[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){
	      hJetTrackMELeading[ibin][ibin2][ibin3]->SetBinContent(k,m,bc);
	      hJetTrackMELeading[ibin][ibin2][ibin3]->SetBinError(k,m,err);
	    }
	  }

	  for(int k = 1; k<hJetTrackMESubLeading[ibin][ibin2][ibin3]->GetNbinsX()+1; k++){
	  
	    bc = me_proj_sub[ibin][ibin2][ibin3]->GetBinContent(k);
	    err = me_proj_sub[ibin][ibin2][ibin3]->GetBinError(k);
	  
	    if(TMath::Abs(hJetTrackMESubLeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(k))<0.3)bc = 1.;

	    for(int m = 1; m<hJetTrackMESubLeading[ibin][ibin2][ibin3]->GetNbinsY()+1; m++){

	      hJetTrackMESubLeading[ibin][ibin2][ibin3]->SetBinContent(k,m,bc);
	      hJetTrackMESubLeading[ibin][ibin2][ibin3]->SetBinError(k,m,err);
	    }
	  }

	}
      

	if(ibin==3){
	    yield_lead[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin-1][ibin2][ibin3]);
	    yield_lead_pTweighted[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin-1][ibin2][ibin3]);
	    yield_sub[ibin][ibin2][ibin3]->Divide(hJetTrackMESubLeading[ibin-1][ibin2][ibin3]);
	    yield_sub_pTweighted[ibin][ibin2][ibin3]->Divide(hJetTrackMESubLeading[ibin-1][ibin2][ibin3]);
	    

	}else{

	    yield_lead[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
	    yield_lead_pTweighted[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
	    yield_sub[ibin][ibin2][ibin3]->Divide(hJetTrackMESubLeading[ibin][ibin2][ibin3]);
	    yield_sub_pTweighted[ibin][ibin2][ibin3]->Divide(hJetTrackMESubLeading[ibin][ibin2][ibin3]);
	}
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Write();
	yield_lead[ibin][ibin2][ibin3] ->Write();
	hJetTrackSignalBackgroundLeading_pTweighted[ibin][ibin2][ibin3]->Write();
	yield_lead_pTweighted[ibin][ibin2][ibin3] ->Write();
	hJetTrackMELeading[ibin][ibin2][ibin3]->Write();


	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Write();
	yield_sub[ibin][ibin2][ibin3]->Write();
	hJetTrackSignalBackgroundSubLeading_pTweighted[ibin][ibin2][ibin3]->Write();
	yield_sub_pTweighted[ibin][ibin2][ibin3]->Write();
	hJetTrackMESubLeading[ibin][ibin2][ibin3]->Write();

	for(int ibin4 = 0; ibin4< nAjBins; ibin4++){

	  if(ibin4==1||ibin4==3)continue;



	  yield_lead_Aj[ibin][ibin2][ibin4][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading_Aj[ibin][ibin2][ibin4][ibin3]->Clone((TString)("Yield_Leading_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	  yield_lead_Aj_pTweighted[ibin][ibin2][ibin4][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Clone((TString)("Yield_pTweighted_Leading_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));


	  me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3] = hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->ProjectionX(Form("me_proj_lead_temp_%d%d%d%d",ibin,ibin2,ibin3,ibin4),1,200);

	  me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->Scale (1./200);
	  me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3] -> Fit("pol0","q","",-me00_range,me00_range);

	  me00_lead = me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetFunction("pol0")->GetParameter(0);




	  yield_sub_Aj[ibin][ibin2][ibin4][ibin3] = (TH2D*) hJetTrackSignalBackgroundSubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Clone((TString)("Yield_SubLeading_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	  yield_sub_Aj_pTweighted[ibin][ibin2][ibin4][ibin3] = (TH2D*) hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Clone((TString)("Yield_pTweighted_SubLeading_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+AjBin_strs[ibin4]+"_"+AjBin_strs[ibin4+2]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	  me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3] = hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->ProjectionX(Form("me_proj_sub_%d%d%d%d",ibin,ibin2,ibin3,ibin4),1,200);

	  me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->Scale (1./200);
	  me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3] -> Fit("pol0","q","",-me00_range,me00_range);

	  me00_sub = me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetFunction("pol0")->GetParameter(0);

	  me00 = (me00_sub+me00_lead)/2.;

	  me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->Scale(1./me00);

	  me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->Scale(1./me00);


	
	  me00_err = TMath::Abs(me00_sub-me00_lead)/me00;

	  if(ibin3==0){
	for(int k = 1; k<hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()/2+1; k++){
	  
	  bc = (me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetBinContent(k)+me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetBinContent(me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k))/2.;
	  err = TMath::Sqrt(me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetBinError(k)*me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetBinError(k)+me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetBinError(me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k)*me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetBinError(me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k))/2.;
	  
	  if(TMath::Abs(hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->GetXaxis()->GetBinCenter(k))<0.3)bc = 1.;

	  for(int m = 1; m<hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsY()+1; m++){

	    hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinContent(k,m,bc);
	    hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinError(k,m,err);
	    hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinContent(me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k,m,bc);
	    hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinError(me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k,m,err);


	  }
	}

	for(int k = 1; k<hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()/2.+1; k++){
	  
	  bc = (me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetBinContent(k)+me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetBinContent(me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k))/2.;
	  err = TMath::Sqrt(me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetBinError(k)*me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetBinError(me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k))/2.;
	  
	  if(TMath::Abs(hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->GetXaxis()->GetBinCenter(k))<0.3)bc = 1.;

	  for(int m = 1; m<hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsY()+1; m++){

	    hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinContent(k,m,bc);
	    hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinError(k,m,err);
	    hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinContent(me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k,m,bc);
	    hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinError(me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1-k,m,err);


	  }
	}
	  }else{
	for(int k = 1; k<hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1; k++){
	  
	  bc = me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetBinContent(k);
	  err =me_proj_lead_Aj[ibin][ibin2][ibin4][ibin3]->GetBinError(k);
	  
	  if(TMath::Abs(hJetTrackMELeading[ibin][ibin2][ibin3]->GetXaxis()->GetBinCenter(k))<0.3)bc = 1.;

	  for(int m = 1; m<hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsY()+1; m++){

	    hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinContent(k,m,bc);
	    hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinError(k,m,err);
	  }
	}

	for(int k = 1; k<hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsX()+1; k++){
	  
	  bc = me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetBinContent(k);
	  err = me_proj_sub_Aj[ibin][ibin2][ibin4][ibin3]->GetBinError(k);
	  
	  if(TMath::Abs(hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->GetXaxis()->GetBinCenter(k))<0.3)bc = 1.;

	  for(int m = 1; m<hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->GetNbinsY()+1; m++){

	    hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinContent(k,m,bc);
	    hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->SetBinError(k,m,err);
	  }
	}

	  }

	  //	  if(ibin3<5){
	  
	  if(ibin3==5){
	    yield_lead_Aj[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
	    yield_lead_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);
	    yield_sub_Aj[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMESubLeading[ibin][ibin2][ibin3]);
	    yield_sub_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMESubLeading[ibin][ibin2][ibin3]);
	    /*	      
	  }else if(ibin==3){
	    yield_lead_Aj[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMELeading_Aj[ibin-1][ibin2][ibin4][ibin3]);
	    yield_lead_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMELeading_Aj[ibin-1][ibin2][ibin4][ibin3]);
	    yield_sub_Aj[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMESubLeading_Aj[ibin-1][ibin2][ibin4][ibin3]);
	    yield_sub_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMESubLeading_Aj[ibin-1][ibin2][ibin4][ibin3]);
	    */
	  }else{
	    yield_lead_Aj[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]);
	    yield_lead_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]);
	    yield_sub_Aj[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]);
	    yield_sub_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Divide(hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]);
	  }

	    //	  }
	  hJetTrackSignalBackgroundSubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Write();
	  yield_sub_Aj[ibin][ibin2][ibin4][ibin3]->Write();
	  hJetTrackSignalBackgroundSubLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Write();
	  yield_sub_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Write();
	  hJetTrackMESubLeading_Aj[ibin][ibin2][ibin4][ibin3]->Write();
	  hJetTrackSignalBackgroundLeading_Aj[ibin][ibin2][ibin4][ibin3]->Write();
	  yield_lead_Aj[ibin][ibin2][ibin4][ibin3] ->Write();

	  hJetTrackSignalBackgroundLeading_Aj_pTweighted[ibin][ibin2][ibin4][ibin3]->Write();
	  yield_lead_Aj_pTweighted[ibin][ibin2][ibin4][ibin3] ->Write();

	  hJetTrackMELeading_Aj[ibin][ibin2][ibin4][ibin3]->Write();


	} //ibin4

       
	//--------------------
	//Draw result histos
	//---------------------

	result_proj_canvas->cd(4*(ibin3+1)-ibin);



	int phil = yield_lead[ibin][ibin2][ibin3]->GetYaxis()->FindBin(-TMath::Pi()/2.+0.01);
	int phir = yield_lead[ibin][ibin2][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2-0.01);
	
	yield_lead_proj[ibin][ibin2][4][ibin3]= (TH1D*)yield_lead_pTweighted[ibin][ibin2][ibin3]->ProjectionX((TString)("yield_lead_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir); 


	yield_lead_proj[ibin][ibin2][4][ibin3]->Rebin(10);

	float 	width = yield_lead_proj[ibin][ibin2][4][ibin3]->GetBinWidth(1);
	yield_lead_proj[ibin][ibin2][4][ibin3]->Scale(1./width);



	//	yield_lead_proj[ibin][ibin2][4][ibin3]->SetMaximum(	yield_lead_proj[ibin][ibin2][4][ibin3]->GetBinContent(5)+20);
	yield_lead_proj[ibin][ibin2][4][ibin3]->SetMarkerStyle(10);
	yield_lead_proj[ibin][ibin2][4][ibin3]->SetMarkerSize(1);

	yield_lead_proj[ibin][ibin2][4][ibin3]->SetMarkerColor(kOrange-2);
	yield_lead_proj[ibin][ibin2][4][ibin3]->SetLineColor(kOrange-2);



	float bc, err;
	for(int k = 1; k<yield_lead_proj[ibin][ibin2][4][ibin3]->GetNbinsX()/2+1; k++){
	  bc = 	(yield_lead_proj[ibin][ibin2][4][ibin3]->GetBinContent(k)+yield_lead_proj[ibin][ibin2][4][ibin3]->GetBinContent(yield_lead_proj[ibin][ibin2][4][ibin3]->GetNbinsX()+1-k))/2.;
	  err = TMath::Sqrt(yield_lead_proj[ibin][ibin2][4][ibin3]->GetBinError(k)*yield_lead_proj[ibin][ibin2][4][ibin3]->GetBinError(k)+yield_lead_proj[ibin][ibin2][4][ibin3]->GetBinError(yield_lead_proj[ibin][ibin2][4][ibin3]->GetNbinsX()/2+1-k)*yield_lead_proj[ibin][ibin2][4][ibin3]->GetBinError(yield_lead_proj[ibin][ibin2][4][ibin3]->GetNbinsX()+1-k))/2.;

	  yield_lead_proj[ibin][ibin2][4][ibin3]->SetBinContent(k,bc);
	  yield_lead_proj[ibin][ibin2][4][ibin3]->SetBinError(k,err);

	  yield_lead_proj[ibin][ibin2][4][ibin3]->SetBinContent(yield_lead_proj[ibin][ibin2][4][ibin3]->GetNbinsX()+1-k,bc);
	  yield_lead_proj[ibin][ibin2][4][ibin3]->SetBinError(yield_lead_proj[ibin][ibin2][4][ibin3]->GetNbinsX()+1-k,err);

	}



	yield_sub_proj[ibin][ibin2][4][ibin3]= (TH1D*)yield_sub_pTweighted[ibin][ibin2][ibin3]->ProjectionX((TString)("yield_sub_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);

	yield_sub_proj[ibin][ibin2][4][ibin3]->Rebin(10);
	yield_sub_proj[ibin][ibin2][4][ibin3]->Scale(1./width);


	//	yield_sub_proj[ibin][ibin2][4][ibin3]->SetMaximum(	yield_sub_proj[ibin][ibin2][4][ibin3]->GetBinContent(5)+20);
	yield_sub_proj[ibin][ibin2][4][ibin3]->SetMarkerStyle(10);
	yield_sub_proj[ibin][ibin2][4][ibin3]->SetMarkerSize(1);
	yield_sub_proj[ibin][ibin2][4][ibin3]->SetMarkerColor(kGreen+3);
	yield_sub_proj[ibin][ibin2][4][ibin3]->SetLineColor(kGreen+3);


	for(int k = 1; k<yield_sub_proj[ibin][ibin2][4][ibin3]->GetNbinsX()/2+1; k++){
	  bc = 	(yield_sub_proj[ibin][ibin2][4][ibin3]->GetBinContent(k)+yield_sub_proj[ibin][ibin2][4][ibin3]->GetBinContent(yield_sub_proj[ibin][ibin2][4][ibin3]->GetNbinsX()+1-k))/2.;
	  err = TMath::Sqrt(yield_sub_proj[ibin][ibin2][4][ibin3]->GetBinError(k)*yield_sub_proj[ibin][ibin2][4][ibin3]->GetBinError(k)+yield_sub_proj[ibin][ibin2][4][ibin3]->GetBinError(yield_sub_proj[ibin][ibin2][4][ibin3]->GetNbinsX()/2+1-k)*yield_sub_proj[ibin][ibin2][4][ibin3]->GetBinError(yield_sub_proj[ibin][ibin2][4][ibin3]->GetNbinsX()+1-k))/2.;

	  yield_sub_proj[ibin][ibin2][4][ibin3]->SetBinContent(k,bc);
	  yield_sub_proj[ibin][ibin2][4][ibin3]->SetBinError(k,err);

	  yield_sub_proj[ibin][ibin2][4][ibin3]->SetBinContent(yield_sub_proj[ibin][ibin2][4][ibin3]->GetNbinsX()+1-k,bc);
	  yield_sub_proj[ibin][ibin2][4][ibin3]->SetBinError(yield_sub_proj[ibin][ibin2][4][ibin3]->GetNbinsX()+1-k,err);

	  
	}


	yield_lead_proj[ibin][ibin2][4][ibin3]->SetAxisRange(-2.49,2.49);
	yield_sub_proj[ibin][ibin2][4][ibin3]->SetAxisRange(-2.49,2.49);
	float ymax = 1.1*TMath::Max(yield_lead_proj[ibin][ibin2][4][ibin3]->GetMaximum(),yield_sub_proj[ibin][ibin2][4][ibin3]->GetMaximum())+5.;
	yield_lead_proj[ibin][ibin2][4][ibin3]->SetMaximum(ymax);

	//	if(ibin3==5)ymax = 10.;

	yield_lead_proj[ibin][ibin2][4][ibin3]->SetAxisRange(-2.49,2.49);
	yield_lead_proj[ibin][ibin2][4][ibin3]->Draw();
	yield_sub_proj[ibin][ibin2][4][ibin3]->Draw("same");


	if(ibin==0&&ibin3==0){
	  l = new TLegend(0.45,0.65,0.85,0.85);
	  l->AddEntry(yield_lead_proj[ibin][ibin2][4][ibin3],"Leading #Delta#eta");
	  l->AddEntry(yield_sub_proj[ibin][ibin2][4][ibin3],"SubLeading #Delta#eta");
	  l->SetLineColor(kWhite);
	  l->SetTextSize(0.05);
	  l->Draw();
	}

	centtex = new TLatex(0.15,0.8,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.15,0.7,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
	ajtex = new TLatex(0.15,0.6,"A_{J} Inclusive");
	ajtex->SetNDC();
	ajtex->Draw();


	dummy->cd();

	

	result_proj_canvas_Aj0_Aj22->cd(4*(ibin3+1)-ibin);


      	
	yield_lead_proj[ibin][ibin2][0][ibin3]= (TH1D*)yield_lead_Aj_pTweighted[ibin][ibin2][0][ibin3]->ProjectionX((TString)("yield_lead_proj_Aj0_Aj22"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);

	yield_lead_proj[ibin][ibin2][0][ibin3]->Rebin(10);

	yield_lead_proj[ibin][ibin2][0][ibin3]->Scale(1./width);

	//	yield_lead_proj[ibin][ibin2][0][ibin3]->SetMaximum(	yield_lead_proj[ibin][ibin2][0][ibin3]->GetBinContent(5)+20.);
	yield_lead_proj[ibin][ibin2][0][ibin3]->SetMarkerStyle(10);
	yield_lead_proj[ibin][ibin2][0][ibin3]->SetMarkerSize(1);
	yield_lead_proj[ibin][ibin2][0][ibin3]->SetMarkerColor(kOrange-2);
	yield_lead_proj[ibin][ibin2][0][ibin3]->SetLineColor(kOrange-2);

	for(int k = 1; k<yield_lead_proj[ibin][ibin2][0][ibin3]->GetNbinsX()/2+1; k++){
	  bc = 	(yield_lead_proj[ibin][ibin2][0][ibin3]->GetBinContent(k)+yield_lead_proj[ibin][ibin2][0][ibin3]->GetBinContent(yield_lead_proj[ibin][ibin2][0][ibin3]->GetNbinsX()+1-k))/2.;
	  err = TMath::Sqrt(yield_lead_proj[ibin][ibin2][0][ibin3]->GetBinError(k)*yield_lead_proj[ibin][ibin2][0][ibin3]->GetBinError(k)+yield_lead_proj[ibin][ibin2][0][ibin3]->GetBinError(yield_lead_proj[ibin][ibin2][0][ibin3]->GetNbinsX()/2+1-k)*yield_lead_proj[ibin][ibin2][0][ibin3]->GetBinError(yield_lead_proj[ibin][ibin2][0][ibin3]->GetNbinsX()+1-k))/2.;

	  yield_lead_proj[ibin][ibin2][0][ibin3]->SetBinContent(k,bc);
	  yield_lead_proj[ibin][ibin2][0][ibin3]->SetBinError(k,err);

	  yield_lead_proj[ibin][ibin2][0][ibin3]->SetBinContent(yield_lead_proj[ibin][ibin2][0][ibin3]->GetNbinsX()+1-k,bc);
	  yield_lead_proj[ibin][ibin2][0][ibin3]->SetBinError(yield_lead_proj[ibin][ibin2][0][ibin3]->GetNbinsX()+1-k,err);

	}

	yield_sub_proj[ibin][ibin2][0][ibin3]= (TH1D*)yield_sub_Aj_pTweighted[ibin][ibin2][0][ibin3]->ProjectionX((TString)("yield_sub_proj_Aj0_Aj22"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);

	yield_sub_proj[ibin][ibin2][0][ibin3]->Rebin(10);
	yield_sub_proj[ibin][ibin2][0][ibin3]->Scale(1./width);

	//	yield_sub_proj[ibin][ibin2][0][ibin3]->SetMaximum(	yield_sub_proj[ibin][ibin2][0][ibin3]->GetBinContent(5)+20.);
	yield_sub_proj[ibin][ibin2][0][ibin3]->SetMarkerStyle(10);
	yield_sub_proj[ibin][ibin2][0][ibin3]->SetMarkerSize(1);
	yield_sub_proj[ibin][ibin2][0][ibin3]->SetMarkerColor(kGreen+3);
	yield_sub_proj[ibin][ibin2][0][ibin3]->SetLineColor(kGreen+3);


	for(int k = 1; k<yield_sub_proj[ibin][ibin2][0][ibin3]->GetNbinsX()/2+1; k++){
	  bc = 	(yield_sub_proj[ibin][ibin2][0][ibin3]->GetBinContent(k)+yield_sub_proj[ibin][ibin2][0][ibin3]->GetBinContent(yield_sub_proj[ibin][ibin2][0][ibin3]->GetNbinsX()+1-k))/2.;
	  err = TMath::Sqrt(yield_sub_proj[ibin][ibin2][0][ibin3]->GetBinError(k)*yield_sub_proj[ibin][ibin2][0][ibin3]->GetBinError(k)+yield_sub_proj[ibin][ibin2][0][ibin3]->GetBinError(yield_sub_proj[ibin][ibin2][0][ibin3]->GetNbinsX()/2+1-k)*yield_sub_proj[ibin][ibin2][0][ibin3]->GetBinError(yield_sub_proj[ibin][ibin2][0][ibin3]->GetNbinsX()+1-k))/2.;

	  yield_sub_proj[ibin][ibin2][0][ibin3]->SetBinContent(k,bc);
	  yield_sub_proj[ibin][ibin2][0][ibin3]->SetBinError(k,err);

	  yield_sub_proj[ibin][ibin2][0][ibin3]->SetBinContent(yield_sub_proj[ibin][ibin2][0][ibin3]->GetNbinsX()+1-k,bc);
	  yield_sub_proj[ibin][ibin2][0][ibin3]->SetBinError(yield_sub_proj[ibin][ibin2][0][ibin3]->GetNbinsX()+1-k,err);

	  
	}

	yield_lead_proj[ibin][ibin2][0][ibin3]->SetAxisRange(-2.49,2.49);
	yield_sub_proj[ibin][ibin2][0][ibin3]->SetAxisRange(-2.49,2.49);
	ymax = 1.1*TMath::Max(yield_lead_proj[ibin][ibin2][0][ibin3]->GetMaximum(),yield_sub_proj[ibin][ibin2][0][ibin3]->GetMaximum())+5.;
	
	//	if(ibin3==5)ymax = 10.;

	yield_lead_proj[ibin][ibin2][0][ibin3]->SetMaximum(ymax);

	yield_lead_proj[ibin][ibin2][0][ibin3]->SetAxisRange(-2.49,2.49);
	yield_lead_proj[ibin][ibin2][0][ibin3]->Draw();
	yield_sub_proj[ibin][ibin2][0][ibin3]->Draw("same");


	if(ibin==0&&ibin3==0){

	  l->Draw();
	}

	centtex = new TLatex(0.15,0.8,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.15,0.7,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
	ajtex = new TLatex(0.15,0.6,"A_{J}<0.22");
	ajtex->SetNDC();
	ajtex->Draw();


	dummy->cd();
	

	

	result_proj_canvas_Aj22_Aj75->cd(4*(ibin3+1)-ibin);


      	
	yield_lead_proj[ibin][ibin2][2][ibin3]= (TH1D*)yield_lead_Aj_pTweighted[ibin][ibin2][2][ibin3]->ProjectionX((TString)("yield_lead_proj_Aj22_Aj75"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);
	yield_lead_proj[ibin][ibin2][2][ibin3]->Rebin(10);
	yield_lead_proj[ibin][ibin2][2][ibin3]->Scale(1./width);
	//	yield_lead_proj[ibin][ibin2][2][ibin3]->SetMaximum(	yield_lead_proj[ibin][ibin2][2][ibin3]->GetBinContent(5)+20.);
	yield_lead_proj[ibin][ibin2][2][ibin3]->SetMarkerStyle(10);
	yield_lead_proj[ibin][ibin2][2][ibin3]->SetMarkerSize(1);
	yield_lead_proj[ibin][ibin2][2][ibin3]->SetMarkerColor(kOrange-2);
	yield_lead_proj[ibin][ibin2][2][ibin3]->SetLineColor(kOrange-2);


	for(int k = 1; k<yield_lead_proj[ibin][ibin2][2][ibin3]->GetNbinsX()/2+1; k++){
	  bc = 	(yield_lead_proj[ibin][ibin2][2][ibin3]->GetBinContent(k)+yield_lead_proj[ibin][ibin2][2][ibin3]->GetBinContent(yield_lead_proj[ibin][ibin2][2][ibin3]->GetNbinsX()+1-k))/2.;
	  err = TMath::Sqrt(yield_lead_proj[ibin][ibin2][2][ibin3]->GetBinError(k)*yield_lead_proj[ibin][ibin2][2][ibin3]->GetBinError(k)+yield_lead_proj[ibin][ibin2][2][ibin3]->GetBinError(yield_lead_proj[ibin][ibin2][2][ibin3]->GetNbinsX()/2+1-k)*yield_lead_proj[ibin][ibin2][2][ibin3]->GetBinError(yield_lead_proj[ibin][ibin2][2][ibin3]->GetNbinsX()+1-k))/2.;

	  yield_lead_proj[ibin][ibin2][2][ibin3]->SetBinContent(k,bc);
	  yield_lead_proj[ibin][ibin2][2][ibin3]->SetBinError(k,err);

	  yield_lead_proj[ibin][ibin2][2][ibin3]->SetBinContent(yield_lead_proj[ibin][ibin2][2][ibin3]->GetNbinsX()+1-k,bc);
	  yield_lead_proj[ibin][ibin2][2][ibin3]->SetBinError(yield_lead_proj[ibin][ibin2][2][ibin3]->GetNbinsX()+1-k,err);

	}



	yield_sub_proj[ibin][ibin2][2][ibin3]= (TH1D*)yield_sub_Aj_pTweighted[ibin][ibin2][2][ibin3]->ProjectionX((TString)("yield_sub_proj_Aj22_Aj75"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);
	yield_sub_proj[ibin][ibin2][2][ibin3]->Rebin(10);
	yield_sub_proj[ibin][ibin2][2][ibin3]->Scale(1./width);


	yield_sub_proj[ibin][ibin2][2][ibin3]->SetMarkerStyle(10);
	yield_sub_proj[ibin][ibin2][2][ibin3]->SetMarkerSize(1);
	yield_sub_proj[ibin][ibin2][2][ibin3]->SetMarkerColor(kGreen+3);
	yield_sub_proj[ibin][ibin2][2][ibin3]->SetLineColor(kGreen+3);




	for(int k = 1; k<yield_sub_proj[ibin][ibin2][2][ibin3]->GetNbinsX()/2+1; k++){
	  bc = 	(yield_sub_proj[ibin][ibin2][2][ibin3]->GetBinContent(k)+yield_sub_proj[ibin][ibin2][2][ibin3]->GetBinContent(yield_sub_proj[ibin][ibin2][2][ibin3]->GetNbinsX()+1-k))/2.;
	  err = TMath::Sqrt(yield_sub_proj[ibin][ibin2][2][ibin3]->GetBinError(k)*yield_sub_proj[ibin][ibin2][2][ibin3]->GetBinError(k)+yield_sub_proj[ibin][ibin2][2][ibin3]->GetBinError(yield_sub_proj[ibin][ibin2][2][ibin3]->GetNbinsX()/2+1-k)*yield_sub_proj[ibin][ibin2][2][ibin3]->GetBinError(yield_sub_proj[ibin][ibin2][2][ibin3]->GetNbinsX()+1-k))/2.;

	  yield_sub_proj[ibin][ibin2][2][ibin3]->SetBinContent(k,bc);
	  yield_sub_proj[ibin][ibin2][2][ibin3]->SetBinError(k,err);

	  yield_sub_proj[ibin][ibin2][2][ibin3]->SetBinContent(yield_sub_proj[ibin][ibin2][2][ibin3]->GetNbinsX()+1-k,bc);
	  yield_sub_proj[ibin][ibin2][2][ibin3]->SetBinError(yield_sub_proj[ibin][ibin2][2][ibin3]->GetNbinsX()+1-k,err);

	  
	}

	yield_lead_proj[ibin][ibin2][2][ibin3]->SetAxisRange(-2.49,2.49);
	yield_sub_proj[ibin][ibin2][2][ibin3]->SetAxisRange(-2.49,2.49);
	ymax = 1.1*TMath::Max(yield_lead_proj[ibin][ibin2][2][ibin3]->GetMaximum(),yield_sub_proj[ibin][ibin2][2][ibin3]->GetMaximum())+5.;
	//	if(ibin3==5)ymax = 10.;

	yield_lead_proj[ibin][ibin2][2][ibin3]->SetMaximum(ymax);

	yield_lead_proj[ibin][ibin2][2][ibin3]->SetAxisRange(-2.49,2.49);
	yield_lead_proj[ibin][ibin2][2][ibin3]->Draw();
	yield_sub_proj[ibin][ibin2][2][ibin3]->Draw("same");


	if(ibin==0&&ibin3==0){

	  l->Draw();
	}
	centtex = new TLatex(0.15,0.8,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.15,0.7,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
	ajtex = new TLatex(0.15,0.6,"A_{J}>0.22");
	ajtex->SetNDC();
	ajtex->Draw();


	dummy->cd();
	
	}//ibin3
	dummy->cd();
      
    } // ibin2

  } // ibin ( centrality ) loop

  
  me_proj_leading_canvas->SaveAs((TString)("All_ME_Projections_Leading_"+datalabel+".png"));
  me_proj_sub_canvas->SaveAs((TString)("All_ME_Projections_SubLeading_"+datalabel+".png"));
  
  result_proj_canvas->SaveAs((TString)("Corrected_dEta_Projections_"+datalabel+".png"));
  result_proj_canvas_Aj0_Aj22->SaveAs((TString)("Corrected_dEta_Projections_"+datalabel+"_Aj0_Aj22.png"));
  result_proj_canvas_Aj22_Aj75->SaveAs((TString)("Corrected_dEta_Projections_"+datalabel+"_Aj22_Aj75.png"));
 
 result_proj_canvas->SaveAs((TString)("Corrected_dEta_Projections_"+datalabel+".pdf"));
  result_proj_canvas_Aj0_Aj22->SaveAs((TString)("Corrected_dEta_Projections_"+datalabel+"_Aj0_Aj22.pdf"));
  result_proj_canvas_Aj22_Aj75->SaveAs((TString)("Corrected_dEta_Projections_"+datalabel+"_Aj22_Aj75.pdf"));
  
  return 0;
} // main loop






  
     
