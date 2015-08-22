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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLatex.h"
#include "THStack.h"

#include <iostream>
#include <vector>
#include <fstream>

#include "../JetTrack2015_functions.h"

//******************************************************
// VERSION TO DO JET-SHAPES FOR PAS
//******************************************************

Int_t jet_shapes_result(bool is_subleading = kFALSE, bool use_highpT_bin = kTRUE){

  TFile *fMC[6], *f_ref[7];

  TCanvas *c_jetshape[5];
   
  TString jetetacut, etalabel,centlabel,pTlabel,Ajlabel;
  float eta_ymax;

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins, limR;
  float deta, dphi, r, bc, bg_err,temp1,temp2, rbin, temperr, err, width_temp_x, width_temp_y, width_temp, norm_temp, zerobin, temp, norm_tot, err_temp, cont;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 7;

  float RBins[17] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.25,1.5};

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  
  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {0.5,1, 2, 3, 4, 8, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<p_{T}^{assoc.}<1 GeV/c","1<p_{T}^{assoc.}<2 GeV/c","2<p_{T}^{assoc.}<3 GeV/c","3<p_{T}^{assoc.}<4 GeV/c","4<p_{T}^{assoc.}<8 GeV/c","p_{T}^{assoc.}>8 GeV/c"};
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
 
 
  TH2D* result[6][nCBins][nPtBins][nTrkPtBins];

  TH2D* resultMC_gen[6][nCBins][nPtBins][nTrkPtBins];
  TH2D* resultMC_reco[6][nCBins][nPtBins][nTrkPtBins];

  TH2D* resultMC2[6][nCBins][nPtBins][nTrkPtBins];

  TH2D* resultMC[6][nCBins][nPtBins][nTrkPtBins];

  TF1 *gaus_eta[4][1][4];
  TF1 *gaus_phi[4][1][4];
 
  float par0, par1, par2, par3, par4;

  TH1D *background_syst_rebin[12][6][4][5];

   
  TH1D* JetShapeMC[6][nCBins][nTrkPtBins];

  TH1D* JetShape[6][nCBins][nTrkPtBins];
  TH1D* JetShape2[6][nCBins][nTrkPtBins];
  TH1D* JetShape_noerr_up[6][nCBins][nTrkPtBins];
  TH1D* JetShape_noerr_down[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ref[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ref_diff[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff[6][nCBins][nTrkPtBins];

  TH1D* JetShape_ref_ratio[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ratio[6][nCBins][nTrkPtBins];
  TH1D* JetShape_ratio2[6][nCBins][nTrkPtBins];

  TH1D* JetShape_diff2[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff_noerr_up[6][nCBins][nTrkPtBins];
  TH1D* JetShape_diff_noerr_down[6][nCBins][nTrkPtBins];
  TH1D* JetShape2_geo[6][nCBins][nTrkPtBins];

  TH1D* JetShape_syst[6][nCBins][nTrkPtBins];
  TGraphAsymmErrors* JetShape_graph[6][nCBins][nTrkPtBins];



  THStack *JetShape_Stack_Up[6][nCBins];
  THStack *JetShape_Diff_Stack_Up[6][nCBins];

  THStack *JetShape_Stack_Down[6][nCBins];
  THStack *JetShape_Diff_Stack_Down[6][nCBins];

 
  double temp_cont, nextr, nextl, cos_weight,me00_range, mc_error;
  
  TString stem, datalabel,me00_range_string,stem_mc;

 
 
  float norm;

  TFile *fin = new TFile("../bg_fit/Dijet_Correlations.root", "READ");
  TFile *fout  = new TFile("Jet_Shapes.root","RECREATE"); 
  TFile *f_spillover = new TFile("../spill_over/Dijet_SpillOvers.root");

  TFile *f_jff_gen2 = new TFile("../me_correct_mc/Pythia_GenJet_GenTrack_Dijet_Correlations.root");  
  TFile *f_jff_reco2 = new TFile("../me_correct_mc/Pythia_RecoJet_RecoTrack_Dijet_Correlations.root"); 
 
  TFile *f_jff_reco_reco = new TFile("../me_correct_mc/HydJet_RecoJet_RecoTrack_Dijet_Correlations.root"); 
  TFile *f_jff_reco = new TFile("../me_correct_mc/HydJet_RecoJet_GenTrack_Sube0_Dijet_Correlations.root"); 
  TFile *f_jff_gen = new TFile("../me_correct_mc/HydJet_GenJet_GenTrack_Sube0_Dijet_Correlations.root");  



  // TFile *f_jff_reco = new TFile("../me_correct_mc/Pythia_RecoJet_GenTrack_Dijet_Correlations.root"); 
 

  

  TF2* ClosureFit = new TF2("ClosureFit", "[0]/2/TMath::Pi()/[1]/[2]*TMath::Exp(-1.*(x*x/[1]/[1]/2))*TMath::Exp(-1.*(y*y/[2]/[2]/2))",-5.,5.,-TMath::Pi()/2.,3*TMath::Pi()/2.);


  //-----------------------
  // Start getting histos
  //-----------------------
  for(int g=0; g<2; g++){

    switch(g){
    case 0:
      stem = "Yield_BkgSub_PbPb_";
      datalabel = "Leading";
      if(is_subleading)     datalabel = "SubLeading";
      break;
    case 1:
      stem = "Yield_BkgSub_pp_";
      datalabel = "Leading";
      if(is_subleading)     datalabel = "SubLeading";
      break;
    }

    for (int ibin=0;ibin<nCBins;ibin++){
  
      for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
    
	for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){

	  result[g][ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString)(stem + CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (stem + datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  if(is_subleading){

	    for(int k = 1; k<result[g][ibin][ibin2][ibin3]->GetNbinsX(); k++){
	      for(int m = 1; m<result[g][ibin][ibin2][ibin3]->GetNbinsY(); m++){
		result[g][ibin][ibin2][ibin3]->SetBinContent(k,m,result[g][ibin][ibin2][ibin3]->GetBinContent(k,m+100));
		result[g][ibin][ibin2][ibin3]->SetBinError(k,m,result[g][ibin][ibin2][ibin3]->GetBinError(k,m+100));


	      }
	    }
	  }


	  width_temp_x = result[g][ibin][ibin2][ibin3]->GetXaxis()->GetBinWidth(1);
	  width_temp_y = result[g][ibin][ibin2][ibin3]->GetYaxis()->GetBinWidth(1);
	  
	  result[g][ibin][ibin2][ibin3]->Scale(1./width_temp_x/width_temp_y);

	  if(g==0)	  resultMC2[g][ibin][0][ibin3] = new TH2D((TString)("Closure_2D_PbPb_"+CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
	  else	  resultMC2[g][ibin][0][ibin3] = new TH2D((TString)("Closure_2D_pp_"+CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
	

	  if(g==0){
	   
	    if(is_subleading)     gaus_eta[ibin][0][ibin3] = (TF1*)f_spillover->Get((TString)("GausFit_Eta_pTweighted_SubLeading_Hydjet_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_" +TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	    else     gaus_eta[ibin][0][ibin3] = (TF1*)f_spillover->Get((TString)("GausFit_Eta_pTweighted_Leading_Hydjet_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));



	    if(is_subleading)    gaus_phi[ibin][0][ibin3] = (TF1*)f_spillover->Get((TString)("GausFit_Phi_pTweighted_SubLeading_Hydjet_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_" +TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
	    else     gaus_phi[ibin][0][ibin3] = (TF1*)f_spillover->Get((TString)("GausFit_Phi_pTweighted_Leading_Hydjet_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] +"_"+TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));


	    par0= gaus_eta[ibin][0][ibin3]->GetParameter(1);
	    par1= gaus_eta[ibin][0][ibin3]->GetParameter(2);
	    par2= gaus_phi[ibin][0][ibin3]->GetParameter(2);

	    ClosureFit->SetParameter(0,par0);
	    ClosureFit->SetParameter(1,par1);
	    ClosureFit->SetParameter(2,par2);
	    

	    resultMC2[g][ibin][0][ibin3]->Eval(ClosureFit);
	        
	    result[g][ibin][ibin2][ibin3]->Add(resultMC2[g][ibin][0][ibin3],-1.);

	  }

	  if(ibin==0&&g==1){

	    if(is_subleading){
	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_gen2->Get((TString)("Raw_Yield_pTweighted_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_SubLeading_PythiaGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco2->Get((TString)("Raw_Yield_pTweighted_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_SubLeading_PythiaGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    }else{

	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_gen2->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco2->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaReco"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  }
	    resultMC_reco[g][ibin][0][ibin3]->Add(resultMC_gen[g][ibin][0][ibin3],-1.);
	  
	    resultMC_reco[g][ibin][0][ibin3]->Scale(1./width_temp_x/width_temp_y);

	  }else if(g==0){

	    if(is_subleading){
	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_gen->Get((TString)("Raw_Yield_pTweighted_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_SubLeading_PythiaForPbPbGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      if(ibin3>2)	      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco_reco->Get((TString)("Raw_Yield_pTweighted_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_SubLeading_PythiaForPbPbGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	      else      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco->Get((TString)("Raw_Yield_pTweighted_SubLeading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_SubLeading_PythiaForPbPbGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	    }else{

	      resultMC_gen[g][ibin][0][ibin3] = (TH2D*)f_jff_gen->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaForPbPbGen"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      if(ibin3>3)      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco_reco->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaReco"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	      else      resultMC_reco[g][ibin][0][ibin3] = (TH2D*)f_jff_reco->Get((TString)("Raw_Yield_pTweighted_Leading_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]))->Clone((TString)("Raw_Yield_pTweighted_Leading_PythiaReco"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  }
	    resultMC_reco[g][ibin][0][ibin3]->Add(resultMC_gen[g][ibin][0][ibin3],-1.);
	  
	    resultMC_reco[g][ibin][0][ibin3]->Scale(1./width_temp_x/width_temp_y);

	  }

	  if(ibin==0){
	    
	    for(int k = 1; k< resultMC_reco[g][ibin][0][ibin3]->GetNbinsX(); k++){
	      for(int m = 1; m< resultMC_reco[g][ibin][0][ibin3]->GetNbinsY(); m++){

		deta = resultMC_reco[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
		dphi = resultMC_reco[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(m);
		
		r = TMath::Sqrt(deta*deta+dphi*dphi);
	  
		if(r>0.3){
		  resultMC_reco[g][ibin][0][ibin3]->SetBinContent(k,m,0.);
		  resultMC_reco[g][ibin][0][ibin3]->SetBinError(k,m,0.);

	      }
	    }
	    }
	  }
	  
	  if(g==0) result[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][ibin][0][ibin3],-1.);
	  else 	  result[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][0][0][ibin3],-1.);
	  

	  resultMC[g][ibin][ibin2][ibin3] = (TH2D*)resultMC2[g][ibin][ibin2][ibin3]->Clone((TString)("Combined_MC_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  resultMC[g][ibin][ibin2][ibin3]->Add(resultMC_reco[g][0][0][ibin3],0.02/0.07);
	 	 

	  if(ibin3==0){
	    result[g][ibin][ibin2][6] = (TH2D*) result[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])));

	    resultMC[g][ibin][ibin2][6] = (TH2D*) result[g][ibin][ibin2][ibin3]->Clone(((TString) ("Summed_result_MC"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1])));
	  }else if(ibin3>0&&(use_highpT_bin||ibin3<5)){
	    result[g][ibin][ibin2][6]->Add(result[g][ibin][ibin2][ibin3]);
	    resultMC[g][ibin][ibin2][6]->Add(resultMC[g][ibin][ibin2][ibin3]);
	  }
	}
    
	for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	  
	  JetShape2[g][ibin][ibin3] = new TH1D((TString)("JetShape2_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",16,RBins);  

	  JetShape2_geo[g][ibin][ibin3] = new TH1D((TString)("JetShape_Geometry_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",16,RBins);  

	  

	  int kmin = result[g][ibin][0][ibin3]->GetXaxis()->FindBin(-TMath::Pi()/2.+.0001);
	  int kmax = result[g][ibin][0][ibin3]->GetXaxis()->FindBin(TMath::Pi()/2.-.0001);
      

	  int lmin = result[g][ibin][0][ibin3]->GetYaxis()->FindBin(-TMath::Pi()/2.+.0001);
	  int lmax = result[g][ibin][0][ibin3]->GetYaxis()->FindBin(TMath::Pi()/2.-.0001);

	
	  for(int k = kmin; k<kmax+1; k++){  //
	    for(int l = lmin; l<lmax+1; l++){

	      deta = result[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
	      dphi = result[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(l);
	   
	      r = TMath::Sqrt(deta*deta+dphi*dphi);
	      
	      
	      rbin = JetShape2[g][ibin][ibin3]->FindBin(r);
	   
	      //catch-all, if you want
	      /*
	      if(rbin > 14){
		rbin = 14;
		cout<<r<<" "<<rbin<<endl;
	      }
	      */
	      bc = result[g][ibin][0][ibin3]->GetBinContent(k,l);
	      err = result[g][ibin][0][ibin3]->GetBinError(k,l);

	      temp1 = JetShape2[g][ibin][ibin3]->GetBinContent(rbin);
	      temp2 = JetShape2_geo[g][ibin][ibin3]->GetBinContent(rbin);
	    
	      temperr = JetShape2[g][ibin][ibin3]->GetBinError(rbin);
	  
	      JetShape2[g][ibin][ibin3]->SetBinContent(rbin,temp1+bc);
	      JetShape2[g][ibin][ibin3]->SetBinError(rbin,TMath::Sqrt(temperr*temperr+err*err));
	      
	      JetShape2_geo[g][ibin][ibin3]->SetBinContent(rbin,temp2+1.);
	      JetShape2_geo[g][ibin][ibin3]->SetBinError(rbin,0.);
	      
	      
	  	  
	    }
	  }


	  JetShapeMC[g][ibin][ibin3] = new TH1D((TString)("JetShape_MC_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]),"",16,RBins);  

	  for(int k = kmin; k<kmax+1; k++){  //
	    for(int l = lmin; l<lmax+1; l++){

	      deta = resultMC[g][ibin][0][ibin3]->GetXaxis()->GetBinCenter(k);
	      dphi = resultMC[g][ibin][0][ibin3]->GetYaxis()->GetBinCenter(l);
	   
	      r = TMath::Sqrt(deta*deta+dphi*dphi);
	  
	      rbin = JetShapeMC[g][ibin][ibin3]->FindBin(r);

	      bc = resultMC[g][ibin][0][ibin3]->GetBinContent(k,l);
	      err = resultMC[g][ibin][0][ibin3]->GetBinError(k,l);

	      temp1 = JetShapeMC[g][ibin][ibin3]->GetBinContent(rbin);
	      temperr = JetShapeMC[g][ibin][ibin3]->GetBinError(rbin);
	 	  
	      JetShapeMC[g][ibin][ibin3]->SetBinContent(rbin,temp1+bc);
	      JetShapeMC[g][ibin][ibin3]->SetBinError(rbin,TMath::Sqrt(temperr*temperr+err*err));
	    }
	  }
	  

	
	  //  JetShape2[g][ibin][ibin3]->Divide(JetShape2_geo[g][ibin][ibin3]);

	  JetShape[g][ibin][ibin3] = (TH1D*)JetShape2[g][ibin][ibin3]->Clone((TString)("JetShape_"+stem+ datalabel+"_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  if(!is_subleading){


	    if(g==0&&ibin==0){
	      if(ibin3>0&&ibin3<5){
		f_ref[ibin3] = new TFile((TString)("JetShapesRef_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]+".root"),"READ");
	      }else if(ibin3==6){
		f_ref[ibin3] = new TFile((TString)("JetShapesRef_"+TrkPtBin_strs[1]+"_"+TrkPtBin_strs[5]+".root"),"READ");
	      }
	    }
	    if(ibin3==1){

	      JetShape_ref[0][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pbpb_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone("JetShape_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]);
	
	      JetShape_ref[1][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pp_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]))->Clone("JetShape_pp_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]);
     
	    }else if(ibin3>0&&ibin3<5){

	      JetShape_ref[0][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pbpb_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone("JetShape_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]);
	
	      JetShape_ref[1][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pp_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone("JetShape_pp_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]);
     

	    }else if(ibin3==6){

	      JetShape_ref[0][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pbpb_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone("JetShape_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[1]+"_"+TrkPtBin_strs[5]);
	      JetShape_ref[1][ibin][ibin3] = (TH1D*)f_ref[ibin3]->Get((TString)("pp_hist_hists_JetShapeDiffParticles_1D"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"))->Clone("JetShape_pp_Ref_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300_"+TrkPtBin_strs[1]+"_"+TrkPtBin_strs[5]);
     

	    }

	    if(ibin3==6){
	   
	      if(ibin==0){
		JetShape_ref[1][ibin][ibin3]->SetBinContent(1,12.485);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(2,4.683);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(3,1.5);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(4,0.722);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(5,0.396);
		JetShape_ref[1][ibin][ibin3]->SetBinContent(6,0.215);

		JetShape_ref[1][ibin][ibin3]->SetBinError(1,0.024);
		JetShape_ref[1][ibin][ibin3]->SetBinError(2,0.014);
		JetShape_ref[1][ibin][ibin3]->SetBinError(3,0.006);
		JetShape_ref[1][ibin][ibin3]->SetBinError(4,0.004);
		JetShape_ref[1][ibin][ibin3]->SetBinError(5,0.002);
		JetShape_ref[1][ibin][ibin3]->SetBinError(6,0.002);



		JetShape_ref[0][ibin][ibin3]->SetBinContent(1,12.636);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(2,4.560);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(3,1.318);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(4,0.707);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(5,0.486);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(6,0.294);

		JetShape_ref[0][ibin][ibin3]->SetBinError(1,0.05);
		JetShape_ref[0][ibin][ibin3]->SetBinError(2,0.028);
		JetShape_ref[0][ibin][ibin3]->SetBinError(3,0.012);
		JetShape_ref[0][ibin][ibin3]->SetBinError(4,0.009);
		JetShape_ref[0][ibin][ibin3]->SetBinError(5,0.009);
		JetShape_ref[0][ibin][ibin3]->SetBinError(6,0.009);

	      }


	      if(ibin==1){
		JetShape_ref[0][ibin][ibin3]->SetBinContent(1,12.93);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(2,4.36);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(3,1.28);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(4,0.685);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(5,0.455);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(6,0.273);

		JetShape_ref[0][ibin][ibin3]->SetBinError(1,0.046);
		JetShape_ref[0][ibin][ibin3]->SetBinError(2,0.025);
		JetShape_ref[0][ibin][ibin3]->SetBinError(3,0.011);
		JetShape_ref[0][ibin][ibin3]->SetBinError(4,0.007);
		JetShape_ref[0][ibin][ibin3]->SetBinError(5,0.007);
		JetShape_ref[0][ibin][ibin3]->SetBinError(6,0.007);

	      }
   
	      if(ibin==2){
		JetShape_ref[0][ibin][ibin3]->SetBinContent(1,13.168);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(2,4.281);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(3,1.278);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(4,0.640);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(5,0.399);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(6,0.233);

		JetShape_ref[0][ibin][ibin3]->SetBinError(1,0.071);
		JetShape_ref[0][ibin][ibin3]->SetBinError(2,0.038);
		JetShape_ref[0][ibin][ibin3]->SetBinError(3,0.016);
		JetShape_ref[0][ibin][ibin3]->SetBinError(4,0.010);
		JetShape_ref[0][ibin][ibin3]->SetBinError(5,0.008);
		JetShape_ref[0][ibin][ibin3]->SetBinError(6,0.007);

	      }


	      if(ibin==3){
		JetShape_ref[0][ibin][ibin3]->SetBinContent(1,(13.086+12.98)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(2,(4.356+4.322)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(3,(1.320+1.396)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(4,(0.637+0.707)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(5,(0.395+0.371)/2.);
		JetShape_ref[0][ibin][ibin3]->SetBinContent(6,(0.230+0.215)/2.);

		JetShape_ref[0][ibin][ibin3]->SetBinError(1,0.13);
		JetShape_ref[0][ibin][ibin3]->SetBinError(2,0.071);
		JetShape_ref[0][ibin][ibin3]->SetBinError(3,0.030);
		JetShape_ref[0][ibin][ibin3]->SetBinError(4,0.019);
		JetShape_ref[0][ibin][ibin3]->SetBinError(5,0.012);
		JetShape_ref[0][ibin][ibin3]->SetBinError(6,0.010);

	      }
	    }
	  

	    if(((ibin3>0&&ibin3<5)||ibin3==6)){

		JetShape_ref[0][ibin][ibin3]->SetMarkerStyle(20);
		JetShape_ref[0][ibin][ibin3]->SetMarkerSize(1);
		JetShape_ref[0][ibin][ibin3]->SetMarkerColor(kCyan);
		JetShape_ref[0][ibin][ibin3]->SetLineColor(kCyan);


		JetShape_ref[1][ibin][ibin3]->SetMarkerStyle(20);
		JetShape_ref[1][ibin][ibin3]->SetMarkerSize(1);
		JetShape_ref[1][ibin][ibin3]->SetMarkerColor(kCyan);
		JetShape_ref[1][ibin][ibin3]->SetLineColor(kCyan);
	
	    }

	  } //!is_subleading


       
	  for(int k = 1; k<JetShape2[g][ibin][ibin3]->GetNbinsX()+1; k++){
	    
	    float temp = JetShape2[g][ibin][ibin3]->GetBinContent(k)/JetShape2[g][ibin][ibin3]->GetBinWidth(k);
	    float err = JetShape2[g][ibin][ibin3]->GetBinError(k)/JetShape2[g][ibin][ibin3]->GetBinWidth(k);
	    
	    JetShape2[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShape2[g][ibin][ibin3]->SetBinError(k,err);

	    JetShape[g][ibin][ibin3]->SetBinContent(k,temp);
	    JetShape[g][ibin][ibin3]->SetBinError(k,err);

	  }	    

	  //	  JetShape2[g][ibin][ibin3]->Scale(0.05);
	  //JetShape[g][ibin][ibin3]->Scale(0.05);
	  //	  JetShape2_geo[g][ibin][ibin3]->Scale(0.05);

	  norm =  JetShape2[g][ibin][ibin3]->Integral(1,6,"width");

	  //	  norm = 160./0.05;
	  

	  JetShape2[g][ibin][ibin3]->Scale(1./norm);
	  //	  JetShape2_geo[g][ibin][ibin3]->Scale(1./norm);

      
	  if(ibin3==6){

	    cout<<"norm is "<<norm<<endl;

	    //	    if(!is_subleading&&!use_highpT_bin) norm = JetShape[g][ibin][ibin3]->GetBinContent(1)/JetShape_ref[g][ibin][ibin3]->GetBinContent(1);	 

	    JetShape[g][ibin][ibin3]->Scale(1./norm);  
	    JetShape[g][ibin][0]->Scale(1./norm);  
	    JetShape[g][ibin][1]->Scale(1./norm);  
	    JetShape[g][ibin][2]->Scale(1./norm);  
	    JetShape[g][ibin][3]->Scale(1./norm);  
	    JetShape[g][ibin][4]->Scale(1./norm);  
	    JetShape[g][ibin][5]->Scale(1./norm);  

	 
	  }
	  
	  
	  if(ibin<3){
	    JetShape2[g][ibin][ibin3]->GetYaxis()->SetLabelSize(0.);
	    JetShape2[g][ibin][ibin3]->GetYaxis()->SetTitleSize(0.);

	  }

	
	  JetShapeMC[g][ibin][ibin3]->Divide(JetShape2_geo[g][ibin][ibin3]);

	  JetShapeMC[g][ibin][ibin3]->Scale(1./norm);
	
	
	  fout->cd();
	  
	  JetShape2[g][ibin][ibin3]->SetAxisRange(0.,0.99,"x");
	  JetShape[g][ibin][ibin3]->SetAxisRange(0.,0.99,"x");
	  JetShape2[g][ibin][ibin3]->Write();

	  if(ibin3<6){
	    TString syst_name_pbpb = make_name("Syst_",g,ibin3,ibin,0,pTlabel, centlabel, Ajlabel);  
	    background_syst_rebin[g][ibin][ibin3][0] = (TH1D*)fin->Get(syst_name_pbpb)->Clone(syst_name_pbpb);
	  }else{
	    TString syst_name_pbpb = make_name("Syst_",g,0,ibin,0,pTlabel, centlabel, Ajlabel);  
	    syst_name_pbpb.ReplaceAll("Pt1","Pt300");
	    background_syst_rebin[g][ibin][ibin3][0] = (TH1D*)  background_syst_rebin[g][ibin][0][0] ->Clone(syst_name_pbpb);
	    for(int k = 0; k<ibin3; k++){
	      background_syst_rebin[g][ibin][ibin3][0]->Add(background_syst_rebin[g][ibin][k][0]);
	    }
	  }



	  JetShape_syst[g][ibin][ibin3] = (TH1D*) JetShape2[g][ibin][ibin3]->Clone((TString)("Jet_Shape_SystErr_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	

	  for(int k = 1; k<JetShape_syst[g][ibin][ibin3]->GetNbinsX()+1; k++){

	    bg_err = background_syst_rebin[g][ibin][ibin3][0]->GetBinError(1)/background_syst_rebin[g][ibin][ibin3][0]->GetBinWidth(1)/5.*JetShape2_geo[g][ibin][ibin3]->GetBinContent(k)/JetShape2_geo[g][ibin][ibin3]->GetBinWidth(k)/norm;

	    //	    bg_err = 0.;

	    bc = TMath::Abs(JetShape_syst[g][ibin][ibin3]->GetBinContent(k))+TMath::Abs(JetShape_syst[g][ibin][ibin3]->GetBinError(k));

	    
	    
	    JetShape_syst[g][ibin][ibin3]->SetBinError(k,TMath::Sqrt(bc*bc*(.05*.05+0.04*0.04+0.03*0.03+0.07*0.07)+bg_err*bg_err)); 
	 	      
	  }
	
	  JetShape_graph[g][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g][ibin][ibin3]);

	  JetShape_graph[g][ibin][ibin3]->SetName((TString)("Jet_Shape_SystErrGraph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	  JetShape_graph[g][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g][ibin][ibin3]);

	  JetShape_graph[g][ibin][ibin3]->SetName((TString)("Jet_Shape_SystErrGraph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  	  
	  for(int k = 1; k<JetShape_syst[g][ibin][ibin3]->GetNbinsX()+1; k++){


	    float evalpt = JetShape_syst[g][ibin][ibin3]->GetBinCenter(k);
 
	    mc_error = JetShapeMC[g][ibin][ibin3]->GetBinContent(k);

	    JetShape_graph[g][ibin][ibin3]->SetPointEYhigh(k-1,TMath::Sqrt(JetShape_syst[g][ibin][ibin3]->GetBinError(k)*JetShape_syst[g][ibin][ibin3]->GetBinError(k)+mc_error*mc_error/4.));

	    JetShape_graph[g][ibin][ibin3]->SetPointEYlow(k-1,TMath::Sqrt(JetShape_syst[g][ibin][ibin3]->GetBinError(k)*JetShape_syst[g][ibin][ibin3]->GetBinError(k)+mc_error*mc_error/4.));

	    JetShape_graph[g][ibin][ibin3]->SetMarkerStyle(20);
	    JetShape_graph[g][ibin][ibin3]->SetMarkerSize(1);
	    JetShape_graph[g][ibin][ibin3]->SetMarkerColor(kWhite);
	    JetShape_graph[g][ibin][ibin3]->SetFillColor(kBlack);
	    JetShape_graph[g][ibin][ibin3]->SetFillStyle(3004);


	    if(g==1){

	      JetShape_syst[g+1][ibin][ibin3] = (TH1D*) JetShape_syst[g-1][ibin][ibin3]->Clone((TString)("Jet_Shape_Syst_Diff_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	      JetShape_syst[g+1][ibin][ibin3]->Divide(JetShape_syst[g][ibin][ibin3]);


	      JetShape_graph[g+1][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g+1][ibin][ibin3]);

	      JetShape_graph[g+1][ibin][ibin3]->SetName((TString)("Jet_Shape_Diff_Graph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


	      JetShape_graph[g+1][ibin][ibin3] = new TGraphAsymmErrors(JetShape_syst[g+1][ibin][ibin3]);

	      JetShape_graph[g+1][ibin][ibin3]->SetName((TString)("Jet_Shape_Diff_Graph_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	  	  
	      for(int k = 1; k<JetShape_syst[g+1][ibin][ibin3]->GetNbinsX()+1; k++){

		mc_error = JetShapeMC[g][ibin][ibin3]->GetBinContent(k);

		float evalpt = JetShape_syst[g+1][ibin][ibin3]->GetBinCenter(k);
	
		JetShape_graph[g+1][ibin][ibin3]->SetPointEYhigh(k-1,TMath::Sqrt(JetShape_syst[g+1][ibin][ibin3]->GetBinError(k)*JetShape_syst[g+1][ibin][ibin3]->GetBinError(k)+mc_error*mc_error/4.));

		JetShape_graph[g+1][ibin][ibin3]->SetPointEYlow(k-1,TMath::Sqrt(JetShape_syst[g+1][ibin][ibin3]->GetBinError(k)*JetShape_syst[g+1][ibin][ibin3]->GetBinError(k)+mc_error*mc_error/4.));


	      }


	      JetShape_graph[g+1][ibin][ibin3]->SetMarkerStyle(20);
	      JetShape_graph[g+1][ibin][ibin3]->SetMarkerSize(1);
	      JetShape_graph[g+1][ibin][ibin3]->SetMarkerColor(kWhite);
	      JetShape_graph[g+1][ibin][ibin3]->SetFillColor(kBlack);
	      JetShape_graph[g+1][ibin][ibin3]->SetFillStyle(3004);

	    }
	  
	  }
	}
      }
    }
  }


  for(int k = 0; k<7; k++){
   
    TString canvas_name = "JetShape_Comparison_";
    
    if(is_subleading)canvas_name+="SubLeading_";
    if(k<6){ 
      canvas_name+=TrkPtBin_strs[k]; canvas_name+="_"; canvas_name+=TrkPtBin_strs[k+1];  
    } else {
      canvas_name+=TrkPtBin_strs[1]; canvas_name+="_"; canvas_name+=TrkPtBin_strs[6];  
    }
   
    c_jetshape[k] = new TCanvas(canvas_name," ",10,10,1500,800);
    c_jetshape[k]->Divide(4,2,0.,0.);


    for(int ibin = 0; ibin<4; ibin++){
  
  
      c_jetshape[k]->cd(4-ibin);

      JetShape2[0][ibin][k]->GetXaxis()->SetTitle("Radius(r)");
      JetShape2[0][ibin][k]->GetXaxis()->CenterTitle(true);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelFont(42);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelOffset(0.02);
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelSize(0.08);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleSize(0.065);
      JetShape2[0][ibin][k]->GetXaxis()->SetTickLength(0.025);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleOffset(1.3);
      JetShape2[0][ibin][k]->GetXaxis()->SetTitleFont(42);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitle("#rho(r)");
      JetShape2[0][ibin][k]->GetYaxis()->CenterTitle(true);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelFont(42);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelOffset(0.004);
      JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.075);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.09);
      JetShape2[0][ibin][k]->GetYaxis()->SetTickLength(0.025);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleOffset(.8);
     
      // JetShape2[0][ibin][k]->SetMinimum(0.008);
      //JetShape2[0][ibin][k]->SetMaximum(14.5);

      JetShape2[0][ibin][k]->SetMinimum(0.01);
      JetShape2[0][ibin][k]->SetMaximum(35.1);

      JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
      JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.08);

      gPad->SetLogy();

      JetShape2[0][ibin][k]->SetMarkerSize(1);
      JetShape2[0][ibin][k]->SetLineColor(kBlack);
      JetShape2[0][ibin][k]->SetMarkerColor(kBlack);
      JetShape2[0][ibin][k]->SetMarkerStyle(24);
     
      JetShape2[1][ibin][k]->SetMarkerSize(1);
      JetShape2[1][ibin][k]->SetLineColor(kBlack);
      JetShape2[1][ibin][k]->SetMarkerColor(kBlack);
      JetShape2[1][ibin][k]->SetMarkerStyle(24);
      

      if(ibin<3){
	JetShape2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      gPad->SetLogy();

      JetShape2[0][ibin][k]->Draw();

      
  
      JetShape2[1][0][k]->Draw("same");

      if(((k>0&&k<5)||k==6)&&!is_subleading){
	JetShape_ref[0][ibin][k]->Draw("same");  



	JetShape_ref[1][0][k]->SetMarkerStyle(20);

	JetShape_ref[1][0][k]->Draw("same");
      }


      if(ibin==3){
	TLegend *legend = new TLegend(0.4,0.55,0.95,0.85);
	if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref[0][ibin][k],"PbPb JetShapes");
	legend->AddEntry(JetShape2[0][ibin][k],"PbPb Present Study");
	if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref[1][ibin][k],"pp JetShapes");
	legend->AddEntry(JetShape2[1][ibin][k],"pp Present Study");
	legend->SetLineColor(kWhite);
	legend->SetTextSize(0.065);
	legend->Draw();
      }else if(ibin==2){
	if(k<4){
	TLatex  *pttex = new TLatex(0.5,0.81,TrkPtBin_labels[k]);
	pttex->SetNDC();
	pttex->SetTextSize(0.075);
	pttex->Draw();
	}else{
	  TLatex  *pttex = new TLatex(0.5,0.81,"1<p_{T}^{assoc.}<8 GeV/c");
	  pttex->SetNDC();
	  pttex->SetTextSize(0.075);
	  pttex->Draw();
	}
      }

    
      TLatex  *centtex ;
      centtex = new TLatex(0.55,0.9,CBin_labels[ibin]);
      centtex->SetNDC();
      centtex->SetTextSize(0.075);
      if(ibin==3){  
	centtex = new TLatex(0.6,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->SetTextSize(0.07);
      }
      
      centtex->Draw();

      c_jetshape[k]->cd(8-ibin);

   
    
      TString diff_name ="JetShape_diff_"; diff_name+=CBin_strs[ibin]; diff_name+= "_"; diff_name += CBin_strs[ibin+1]; diff_name+= k;

      TString diff2_name ="JetShape_diff2_"; diff2_name+=CBin_strs[ibin]; diff2_name+= "_"; diff2_name += CBin_strs[ibin+1]; diff2_name+= k;
      
     
      JetShape_diff[0][ibin][k] = (TH1D*)JetShape[0][ibin][k]->Clone(diff_name);
      JetShape_diff2[0][ibin][k] = (TH1D*)JetShape2[0][ibin][k]->Clone(diff2_name);

      JetShape_diff[0][ibin][k]->Add( JetShape[1][0][k],-1. );
      JetShape_diff2[0][ibin][k]->Add( JetShape2[1][0][k],-1. );
   

      JetShape_diff2[0][ibin][k]->SetMinimum(-2);
      JetShape_diff2[0][ibin][k]->SetMaximum(4.);
     

      JetShape_diff2[0][ibin][k]->SetMarkerStyle(24);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitle("#rho(r)_{PbPb} - #rho(r)_{pp}");
      JetShape_diff2[0][ibin][k]->GetXaxis()->SetNdivisions(408);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetNdivisions(612);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
      JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.08);

      if(ibin==3){
	JetShape_diff2[0][ibin][k]->GetXaxis()->SetTitleSize(0.9*JetShape_diff2[0][ibin][k]->GetXaxis()->GetTitleSize());
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.9*JetShape_diff2[0][ibin][k]->GetYaxis()->GetTitleSize());
	JetShape_diff2[0][ibin][k]->GetXaxis()->SetLabelSize(0.9*JetShape_diff2[0][ibin][k]->GetXaxis()->GetLabelSize());
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.9*JetShape_diff2[0][ibin][k]->GetYaxis()->GetLabelSize());
      }

      if(ibin<3){
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape_diff2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      JetShape_diff2[0][ibin][k]->Draw();

      if(((k>0&&k<5)||k==6)&&!is_subleading){
	diff_name ="JetShape_ref_diff_"; diff_name+=CBin_strs[ibin]; diff_name+= "_"; diff_name += CBin_strs[ibin+1]; diff_name+= k;
      
     
	JetShape_ref_diff[0][ibin][k] = (TH1D*)JetShape_ref[0][ibin][k]->Clone(diff_name);

	JetShape_ref_diff[0][ibin][k]->Add( JetShape_ref[1][0][k],-1. );

	JetShape_ref_diff[0][ibin][k]->SetMarkerStyle(20);
	JetShape_ref_diff[0][ibin][k]->Draw("same");
      }



      TString ratio_name ="JetShape_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= k;

      TString ratio2_name ="JetShape_ratio2_"; ratio2_name+=CBin_strs[ibin]; ratio2_name+= "_"; ratio2_name += CBin_strs[ibin+1]; ratio2_name+= k;
      
     
      JetShape_ratio[0][ibin][k] = (TH1D*)JetShape[0][ibin][k]->Clone(ratio_name);
      JetShape_ratio2[0][ibin][k] = (TH1D*)JetShape2[0][ibin][k]->Clone(ratio2_name);

      JetShape_ratio[0][ibin][k]->Divide( JetShape[1][0][k]);
      JetShape_ratio2[0][ibin][k]->Divide( JetShape2[1][0][k] );
   

      JetShape_ratio2[0][ibin][k]->SetMinimum(-1.);
      JetShape_ratio2[0][ibin][k]->SetMaximum(9.);
     

      JetShape_ratio2[0][ibin][k]->SetMarkerStyle(20);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitle("#rho(r)_{PbPb} - #rho(r)_{pp}");
      JetShape_ratio2[0][ibin][k]->GetXaxis()->SetNdivisions(408);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetNdivisions(612);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.07);
      JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.08);

      if(ibin==3){
	JetShape_ratio2[0][ibin][k]->GetXaxis()->SetTitleSize(0.9*JetShape_ratio2[0][ibin][k]->GetXaxis()->GetTitleSize());
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.9*JetShape_ratio2[0][ibin][k]->GetYaxis()->GetTitleSize());
	JetShape_ratio2[0][ibin][k]->GetXaxis()->SetLabelSize(0.9*JetShape_ratio2[0][ibin][k]->GetXaxis()->GetLabelSize());
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.9*JetShape_ratio2[0][ibin][k]->GetYaxis()->GetLabelSize());
      }

      if(ibin<3){
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetLabelSize(0.);
	JetShape_ratio2[0][ibin][k]->GetYaxis()->SetTitleSize(0.);
      }

      JetShape_ratio2[0][ibin][k]->Draw();

      if(((k>0&&k<5)||k==6)&&!is_subleading){
	ratio_name ="JetShape_ref_ratio_"; ratio_name+=CBin_strs[ibin]; ratio_name+= "_"; ratio_name += CBin_strs[ibin+1]; ratio_name+= k;
      
     
	JetShape_ref_ratio[0][ibin][k] = (TH1D*)JetShape_ref[0][ibin][k]->Clone(ratio_name);

	JetShape_ref_ratio[0][ibin][k]->Divide( JetShape_ref[1][0][k] );

	JetShape_ref_ratio[0][ibin][k]->SetMarkerStyle(20);
	JetShape_ref_ratio[0][ibin][k]->Draw("same");
      }
 
      JetShape2[0][ibin][k]->GetXaxis()->SetLabelSize(0.);

     
      if(ibin==3){
	TLegend *legend = new TLegend(0.2,0.7,0.4,0.85);
	if(((k>0&&k<5)||k==6)&&!is_subleading)	legend->AddEntry(JetShape_ref_ratio[0][ibin][k],"PbPb - pp JetShapes");
	legend->AddEntry(JetShape_ratio2[0][ibin][k],"PbPb - pp Present Study");
	legend->SetLineColor(kWhite);
	legend->SetTextSize(0.06);
	legend->Draw();
      }      

    }
 

    canvas_name+=".pdf";
    
    c_jetshape[k]->SaveAs(canvas_name);
    
    canvas_name.ReplaceAll(".pdf",".png");
    c_jetshape[k]->SaveAs(canvas_name);

  }

  for(int ibin = 0; ibin<4; ibin++){
 
    for(int k = 0; k<7; k++){

 
      JetShape_noerr_up[0][ibin][k] = new TH1D((TString)("JetShape_PbPb_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",16,RBins);  
      JetShape_noerr_down[0][ibin][k] = new TH1D((TString)("JetShape_PbPb_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",16,RBins);  
   
      JetShape_noerr_up[1][ibin][k] = new TH1D((TString)("JetShape_pp_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",16,RBins);  
      JetShape_noerr_down[1][ibin][k] = new TH1D((TString)("JetShape_pp_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",16,RBins);  
    
      JetShape_diff_noerr_up[0][ibin][k] = new TH1D((TString)("JetShape_Diff_NoErr_Up_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",16,RBins);  
      JetShape_diff_noerr_down[0][ibin][k] = new TH1D((TString)("JetShape_Diff_NoErr_Down_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]+"_"+TrkPtBin_strs[k]+"_" +TrkPtBin_strs[k+1]),"",16,RBins);  
   
      

      for(int l = 0; l< JetShape[0][ibin][k]->GetNbinsX()+1; l++){
      
	bc = JetShape[0][ibin][k]->GetBinContent(l);
          
	if(bc>0){ 
	  JetShape_noerr_up[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_down[0][ibin][k]->SetBinContent(l,0.);	

	}else{
	  JetShape_noerr_down[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_up[0][ibin][k]->SetBinContent(l,0.);	
	}

  
	bc = JetShape_diff[0][ibin][k]->GetBinContent(l);

    
	if(bc>0){ 
	  JetShape_diff_noerr_up[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_diff_noerr_down[0][ibin][k]->SetBinContent(l,0.);	

	}else{
	  JetShape_diff_noerr_down[0][ibin][k]->SetBinContent(l,bc);	
	  JetShape_diff_noerr_up[0][ibin][k]->SetBinContent(l,0.);	
	} 

	bc = JetShape[1][ibin][k]->GetBinContent(l);
      
	if(bc>0){ 
	  JetShape_noerr_up[1][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_down[1][ibin][k]->SetBinContent(l,0.);	

	}else{
	  JetShape_noerr_down[1][ibin][k]->SetBinContent(l,bc);	
	  JetShape_noerr_up[1][ibin][k]->SetBinContent(l,0.);	
	}
      }
     
      switch(k){
      case 0:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kBlue-9);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kBlue-9);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kBlue-9);
	break;
      case 1:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kYellow-9);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kYellow-9);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kYellow-9);
	break;
      case 2:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kOrange+1);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kOrange+1);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kOrange+1);
	break;
      case 3:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kViolet-5);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kViolet-5);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kViolet-5);
	break;
      case 4:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kGreen+3);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kGreen+3);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kGreen+3);
	break;
      case 5:
	JetShape_noerr_up[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_up[1][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_down[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_noerr_down[1][ibin][k]->SetFillColor(kRed+1);
	JetShape_diff_noerr_up[0][ibin][k]->SetFillColor(kRed+1);
	JetShape_diff_noerr_down[0][ibin][k]->SetFillColor(kRed+1);
	break;
      default:
	break;
      }

      JetShape_noerr_up[0][ibin][k]->SetFillStyle(1001);
      JetShape_noerr_up[1][ibin][k]->SetFillStyle(1001);
      JetShape_diff_noerr_up[0][ibin][k]->SetFillStyle(1001);

     	  
      JetShape_noerr_down[0][ibin][k]->SetFillStyle(1001);
      JetShape_noerr_down[1][ibin][k]->SetFillStyle(1001);
      JetShape_diff_noerr_down[0][ibin][k]->SetFillStyle(1001);


      
    } //k
    
    cout<<"ready to make stack for "<<ibin<<endl;

    JetShape_Stack_Up[0][ibin]=new THStack((TString)("JetShapeStack_PbPb_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Stack_Up[1][ibin]=new THStack((TString)("JetShapeStack_pp_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Diff_Stack_Up[0][ibin]=new THStack((TString)("JetShapeStack_Diff_Up_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");


    JetShape_Stack_Down[0][ibin]=new THStack((TString)("JetShapeStack_PbPb_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Stack_Down[1][ibin]=new THStack((TString)("JetShapeStack_pp_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");
    JetShape_Diff_Stack_Down[0][ibin]=new THStack((TString)("JetShapeStack_Diff_Down_"+CBin_strs[ibin]+CBin_strs[ibin+1]),"");

    cout<<"made stack"<<endl;

    if(!use_highpT_bin){
      for(int k = 0; k<5; k++){
	JetShape_Stack_Up[0][ibin]->Add(JetShape_noerr_up[0][ibin][4-k]);
	JetShape_Stack_Up[1][ibin]->Add(JetShape_noerr_up[1][ibin][4-k]);
	JetShape_Diff_Stack_Up[0][ibin]->Add(JetShape_diff_noerr_up[0][ibin][4-k]);

	JetShape_Stack_Down[0][ibin]->Add(JetShape_noerr_down[0][ibin][4-k]);
	JetShape_Stack_Down[1][ibin]->Add(JetShape_noerr_down[1][ibin][4-k]);
	JetShape_Diff_Stack_Down[0][ibin]->Add(JetShape_diff_noerr_down[0][ibin][4-k]);
      }

    }else{

      for(int k = 0; k<6; k++){
	JetShape_Stack_Up[0][ibin]->Add(JetShape_noerr_up[0][ibin][k]);
	JetShape_Stack_Up[1][ibin]->Add(JetShape_noerr_up[1][ibin][k]);
	JetShape_Diff_Stack_Up[0][ibin]->Add(JetShape_diff_noerr_up[0][ibin][k]);

	JetShape_Stack_Down[0][ibin]->Add(JetShape_noerr_down[0][ibin][k]);
	JetShape_Stack_Down[1][ibin]->Add(JetShape_noerr_down[1][ibin][k]);
	JetShape_Diff_Stack_Down[0][ibin]->Add(JetShape_diff_noerr_down[0][ibin][k]);
      }


    }

  
    JetShape_Stack_Up[0][ibin]->SetMaximum(14.1);
    JetShape_Stack_Up[1][ibin]->SetMaximum(14.1);
    JetShape_Stack_Up[0][ibin]->SetMinimum(-.95);
    JetShape_Stack_Up[1][ibin]->SetMinimum(-.95);
  
  
    
    if(is_subleading){
      JetShape_Diff_Stack_Up[0][ibin]->SetMinimum(-1.41);
      JetShape_Diff_Stack_Up[0][ibin]->SetMaximum(1.41);
   
    }else{
      JetShape_Diff_Stack_Up[0][ibin]->SetMinimum(-1.41);
      JetShape_Diff_Stack_Up[0][ibin]->SetMaximum(1.41);
    
    }
    

    if(use_highpT_bin){

      JetShape_Stack_Up[0][ibin]->SetMaximum(37.);
      JetShape_Stack_Up[1][ibin]->SetMaximum(37.);
      JetShape_Stack_Up[0][ibin]->SetMinimum(0.01);
      JetShape_Stack_Up[1][ibin]->SetMinimum(0.01);
  

    
      if(is_subleading){
	JetShape_Diff_Stack_Up[0][ibin]->SetMinimum(-1.1);
	JetShape_Diff_Stack_Up[0][ibin]->SetMaximum(1.45);
   
      }else{
	JetShape_Diff_Stack_Up[0][ibin]->SetMinimum(-1.1);
	JetShape_Diff_Stack_Up[0][ibin]->SetMaximum(1.45);
    
      }


    }
    
  


  } //ibin
  

  cout<<"drawing"<<endl;

  TCanvas *PAS_plot = new TCanvas("JetShape_ForPAS","",1200,800);
  PAS_plot->Divide(3,2,0.,0.);


  PAS_plot->cd(4);
  TLegend *legend = new TLegend(0.1,0.3,0.9,1.);
  legend->AddEntry(JetShape_noerr_up[0][0][0],"0.5<p_{T}^{assoc.}<1 GeV/c","f");
  legend->AddEntry(JetShape_noerr_up[0][0][1],"1<p_{T}^{assoc.}<2 GeV/c","f");
  legend->AddEntry(JetShape_noerr_up[0][0][2],"2<p_{T}^{assoc.}<3 GeV/c","f");
  legend->AddEntry(JetShape_noerr_up[0][0][3],"3<p_{T}^{assoc.}<4 GeV/c","f");
  legend->AddEntry(JetShape_noerr_up[0][0][4],"4<p_{T}^{assoc.}<8 GeV/c","f");

  if(!use_highpT_bin){
    legend->AddEntry(JetShape_graph[0][0][6],"Total 0.5<p_{T}^{assoc.}<8 GeV/c","lpfe");
  }else{
    legend->AddEntry(JetShape_noerr_up[0][0][5],"p_{T}^{assoc.}>8 GeV/c","f");
    legend->AddEntry(JetShape_graph[0][0][6],"Total p_{T}^{assoc.}>0.5 GeV/c","lpfe");
  }
  legend->SetTextSize(0.055);
  legend->SetLineColor(kWhite);
  legend->Draw();

  
  TLatex *aj_tex = new TLatex(0.05,0.2,"A_{J} Inclusive");
  aj_tex->SetTextSize(0.08);
  aj_tex->SetLineColor(kWhite);
  aj_tex->SetNDC();
  aj_tex->Draw();
    
  TLatex *type_tex;
  if(is_subleading) type_tex = new TLatex(0.05,0.1,"SubLeading Jet Shape");
  else type_tex = new TLatex(0.05,0.1,"Leading Jet Shape");
  type_tex->SetTextSize(0.08);
  type_tex->SetLineColor(kWhite);
  type_tex->SetNDC();
  type_tex->Draw();
    


  
  PAS_plot->cd(1);
  
  JetShape_Stack_Up[1][0]->Draw();
  JetShape_Stack_Up[1][0]->GetXaxis()->SetRangeUser(0.,0.99);

  gPad->SetLogy();
  JetShape_Stack_Up[1][0]->GetYaxis()->SetLabelSize(0.07);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitleSize(0.07);
  JetShape_Stack_Up[1][0]->GetYaxis()->SetTitle("#rho(r)");
  JetShape_Stack_Up[1][0]->GetYaxis()->CenterTitle();
  JetShape_Stack_Down[1][0]->Draw("same");
  JetShape_Stack_Up[1][0]->Draw("same");

 

  JetShape_graph[1][0][6]->Draw("same e2 P");
  JetShape2[1][0][6]->Draw("same");
  if(!is_subleading)JetShape_ref[1][0][6]->Draw("same");

  TLatex  *label_pp = new TLatex(0.2,0.9,"pp Refrence");
  label_pp->SetTextSize(0.08);
  label_pp->SetLineColor(kWhite);
  label_pp->SetNDC();
  label_pp->Draw();

  TLine *l_dr = new TLine(0.,0.,TMath::Pi()/2.,0.);
  l_dr->SetLineStyle(2);
  l_dr->Draw();


  TLatex *cms_tex_dphi = new TLatex(0.2,0.8,"CMS Preliminary");
  cms_tex_dphi->SetTextSize(0.08);
  cms_tex_dphi->SetLineColor(kWhite);
  cms_tex_dphi->SetNDC();
  cms_tex_dphi->Draw(); 



  PAS_plot->cd(2);


  
  JetShape_Stack_Up[0][3]->Draw();

  JetShape_Stack_Up[0][3]->GetXaxis()->SetRangeUser(0.,0.99);

  gPad->SetLogy();
  JetShape_Stack_Up[0][3]->GetYaxis()->SetLabelSize(0.);
  JetShape_Stack_Down[0][3]->Draw("same");

 

  JetShape_graph[0][3][6]->Draw("same e2 P");
  JetShape2[0][3][6]->Draw("same");

  if(!is_subleading) JetShape_ref[0][3][6]->Draw("same");

  TLatex  *label_per = new TLatex(0.05,0.9,"PbPb Cent. 50-100%");
  label_per->SetTextSize(0.08);
  label_per->SetLineColor(kWhite);
  label_per->SetNDC();
  label_per->Draw()
;


  PAS_plot->cd(3);
  
  JetShape_Stack_Up[0][0]->Draw();
  gPad->SetLogy();
  JetShape_Stack_Up[0][0]->GetYaxis()->SetLabelSize(0.);
  JetShape_Stack_Down[0][0]->Draw("same");

  JetShape_Stack_Up[0][0]->GetXaxis()->SetRangeUser(0.,0.99);

   
   JetShape_graph[0][0][6]->Draw("same e2 P");
   JetShape2[0][0][6]->Draw("same");

if(!is_subleading) JetShape_ref[0][0][6]->Draw("same");

  TLatex  *label_cent = new TLatex(0.05,0.9,"PbPb Cent. 0-10%");
  label_cent->SetTextSize(0.08);
  label_cent->SetLineColor(kWhite);
  label_cent->SetNDC();
  label_cent->Draw();


  
  PAS_plot->cd(5);


  JetShape_ratio[0][3][6]->SetMinimum(0.);
  JetShape_ratio[0][3][6]->SetMaximum(8.5);


 JetShape_ratio[0][3][6]->Draw();
 JetShape_ratio[0][3][6]->GetYaxis()->SetLabelSize(0.07); 
 JetShape_ratio[0][3][6]->GetYaxis()->CenterTitle();
 JetShape_ratio[0][3][6]->GetYaxis()->SetTitle("#rho(r)_{PbPb}/#rho(r)_{pp}");
 JetShape_ratio[0][3][6]->GetYaxis()->SetTitleSize(0.07);
 JetShape_ratio[0][3][6]->GetYaxis()->SetTitleOffset(1.2);

 JetShape_ratio[0][3][6]->GetXaxis()->SetLabelSize(0.06); 
 JetShape_ratio[0][3][6]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][3][6]->GetXaxis()->SetTitle("r");
 JetShape_ratio[0][3][6]->GetXaxis()->SetTitleSize(0.07);
 JetShape_ratio[0][3][6]->GetXaxis()->SetTitleOffset(1.);
 JetShape_ratio[0][3][6]->GetXaxis()->CenterTitle();
 

 JetShape_ratio[0][3][6]->SetMarkerStyle(20);
 JetShape_ratio[0][3][6]->SetMarkerSize(1);
 JetShape_ratio[0][3][6]->SetMarkerColor(kBlack);
  JetShape_ratio[0][3][6]->Draw();

  JetShape_graph[2][3][6]->Draw("same e2 P");
  JetShape_ratio[0][3][6]->Draw("same");
  if(!is_subleading) JetShape_ref_ratio[0][3][6]->Draw("same");


  TLatex  *label_ratio = new TLatex(0.05,0.9,"");
  label_ratio->SetTextSize(0.08);
  label_ratio->SetLineColor(kWhite);
  label_ratio->SetNDC();
  label_ratio->Draw();

  TLegend *legend_ratio = new TLegend(0.02,0.75,0.9,0.9);
  legend_ratio->AddEntry(JetShape_ratio[0][3][6],"Jet Shape Ratio PbPb/pp");
  if(!is_subleading)  legend_ratio->AddEntry(JetShape_ref_ratio[0][3][6],"CMS Published: PLB 730 (2014)");
  legend_ratio->SetLineColor(kWhite);
  legend_ratio->SetTextSize(0.055);
  legend_ratio->Draw("same");



  TLine *line = new TLine(0.,1.,1.,1.);
  line->SetLineStyle(2);
  line->Draw();

  PAS_plot->cd(6);


  JetShape_ratio[0][0][6]->SetMinimum(0.);
  JetShape_ratio[0][0][6]->SetMaximum(8.5);


 JetShape_ratio[0][0][6]->SetMarkerStyle(20);
 JetShape_ratio[0][0][6]->SetMarkerSize(1);
 JetShape_ratio[0][0][6]->SetMarkerColor(kBlack);
 

 JetShape_ratio[0][0][6]->Draw();
 JetShape_ratio[0][0][6]->GetYaxis()->SetLabelSize(0.0); 
  JetShape_ratio[0][0][6]->GetXaxis()->SetLabelSize(0.06); 
 JetShape_ratio[0][0][6]->GetXaxis()->CenterTitle();
 JetShape_ratio[0][0][6]->GetXaxis()->SetTitle("r");
 JetShape_ratio[0][0][6]->GetXaxis()->SetTitleSize(0.07);
 JetShape_ratio[0][0][6]->GetXaxis()->SetTitleOffset(1.);
 JetShape_ratio[0][0][6]->GetXaxis()->CenterTitle();

 
  JetShape_ratio[0][0][6]->Draw();
  JetShape_graph[2][0][6]->Draw("same e2 P");
  JetShape_ratio[0][0][6]->Draw("same");
  if(!is_subleading) JetShape_ref_ratio[0][0][6]->Draw("same");
 

  line->Draw();

  if(!use_highpT_bin){
    if(!is_subleading){
      PAS_plot->SaveAs("JetShapes_PAS.png");
      PAS_plot->SaveAs("JetShapes_PAS.pdf");
    }else{
      PAS_plot->SaveAs("JetShapes_SubLeading_PAS.png");
      PAS_plot->SaveAs("JetShapes_SubLeading_PAS.pdf");

    }
  

  }else{
   if(!is_subleading){
      PAS_plot->SaveAs("JetShapes_WithHighpT_PAS.png");
      PAS_plot->SaveAs("JetShapes_WithHighpT_PAS.pdf");
    }else{
      PAS_plot->SaveAs("JetShapes_SubLeading_WithHighpT_PAS.png");
      PAS_plot->SaveAs("JetShapes_SubLeading_WithHighpT_PAS.pdf");

    }
  }

  cout<<"pp:"<<endl;
  cout<<"Syst. norm uncertainty: "<<( JetShape_syst[1][3][6]->GetBinError(1)+JetShape_syst[1][3][6]->GetBinError(2)+JetShape_syst[1][3][6]->GetBinError(3)+JetShape_syst[1][3][6]->GetBinError(4)+JetShape_syst[1][3][6]->GetBinError(5)+JetShape_syst[1][3][6]->GetBinError(6))*0.05<<endl;

  cout<<"Stat. norm uncertainty: "<<( JetShape2[1][3][6]->GetBinError(1)+JetShape2[1][3][6]->GetBinError(2)+JetShape2[1][3][6]->GetBinError(3)+JetShape2[1][3][6]->GetBinError(4)+JetShape2[1][3][6]->GetBinError(5)+JetShape2[1][3][6]->GetBinError(6))*0.05<<endl;


  cout<<"PbPb Cent 50-100%:"<<endl;
  cout<<"Syst. norm uncertainty: "<<( JetShape_syst[0][3][6]->GetBinError(1)+JetShape_syst[0][3][6]->GetBinError(2)+JetShape_syst[0][3][6]->GetBinError(3)+JetShape_syst[0][3][6]->GetBinError(4)+JetShape_syst[0][3][6]->GetBinError(5)+JetShape_syst[0][3][6]->GetBinError(6))*0.05<<endl;

  cout<<"Stat. norm uncertainty: "<<( JetShape2[0][3][6]->GetBinError(1)+JetShape2[0][3][6]->GetBinError(2)+JetShape2[0][3][6]->GetBinError(3)+JetShape2[0][3][6]->GetBinError(4)+JetShape2[0][3][6]->GetBinError(5)+JetShape2[0][3][6]->GetBinError(6))*0.05<<endl;



  cout<<"PbPb Cent 0-10%:"<<endl;
  cout<<"Syst. norm uncertainty: "<<( JetShape_syst[0][0][6]->GetBinError(1)+JetShape_syst[0][0][6]->GetBinError(2)+JetShape_syst[0][0][6]->GetBinError(3)+JetShape_syst[0][0][6]->GetBinError(4)+JetShape_syst[0][0][6]->GetBinError(5)+JetShape_syst[0][0][6]->GetBinError(6))*0.05<<endl;

  cout<<"Stat. norm uncertainty: "<<( JetShape2[0][0][6]->GetBinError(1)+JetShape2[0][0][6]->GetBinError(2)+JetShape2[0][0][6]->GetBinError(3)+JetShape2[0][0][6]->GetBinError(4)+JetShape2[0][0][6]->GetBinError(5)+JetShape2[0][0][6]->GetBinError(6))*0.05<<endl;


  
 
  
  return 0;
}// main loop
