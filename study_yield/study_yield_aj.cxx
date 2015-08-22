#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLatex.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "../JetTrack2015_functions.h"


using namespace std;

Int_t study_yield_aj(bool is_number = kFALSE, bool use_highpT_bin = kFALSE){

  gROOT->ForceStyle();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);

  TCanvas *c_yields_phi[5];
  TCanvas *c_PbPb_pp_phi[5];

  TCanvas *c_yields_eta[5];
  TCanvas *c_PbPb_pp_eta[5];


  TH2D *result[12][6][4][5];

  TH1D *signal_dPhi[12][6][4][5];
  TH1D *signal_dPhi_syst[12][6][4][5];
  TH1D *signal_dPhi_rebin[12][6][4][5];
  TH1D *signal_dPhi_PbPb_pp[12][6][4][5];  

  TH1D *background_diff_rebin[12][6][4][5];
  TH1D *background_syst_rebin[12][6][4][5];
  TH1D *background_diff_PbPb_pp[12][6][4][5];

  TH1D *signal_dEta[12][6][4][5];
  TH1D *signal_dEta_syst[12][6][4][5];
  TH1D *signal_dEta_rebin[12][6][4][5];
  TH1D *signal_dEta_PbPb_pp[12][6][4][5];  

  TH1D *spill_over_dEta[12][6][4][5];
  TH1D *spill_over_dPhi[12][6][4][5];
 
  TH1D *jff_residual_dEta[12][6][4][5];
  TH1D *jff_residual_dPhi[12][6][4][5];


  TFile *f_in = new TFile("../bg_fit/Dijet_Correlations.root");
  if(is_number)f_in = new TFile("../bg_fit/Dijet_Correlations_NoPtWeight.root");

  TFile *f_spillover = new TFile("../spill_over/Dijet_SpillOvers.root");
  if(is_number)f_spillover = new TFile("../spill_over/Dijet_SpillOvers_NoPtWeight.root");

  TFile *f_jff_reco = new TFile("../jff_residual/All_JFFResiduals_RecoReco.root");
  TFile *f_jff_gen = new TFile("../jff_residual/All_JFFResiduals_RecoGen.root");

  // if(is_number)f_jff = new TFile("../jff_residual/All_JFFResiduals_NoPtWeight.root");

  TFile *f_out;
  
  if(is_number)   f_out = new TFile("Dijet_Results_NoPtWeight.root","RECREATE");
  else if(use_highpT_bin) f_out = new TFile("Dijet_Results_WithHighpTbin.root","RECREATE");
  else f_out = new TFile("Dijet_Results.root","RECREATE");

  TString in_name, pTlabel,centlabel,Ajlabel;

  TString AjBin_strs[4] = {"Aj0","Aj22","Aj75","AjInclusive"};

  int lbin, rbin;

  float bin_width_phi, bin_width_eta, diff_max, diff_min, signal_min, signal_max, bc, err;

  float rel_err =TMath::Sqrt(0.05*0.05+0.04*0.04+0.03*0.03+.02*.02);

  float me_err[12][6][4][3];
  
  for(int l = 0; l<3; l++){

    c_yields_phi[l] = new TCanvas(Form("yields_phi%d",l),"",10,10,1500,2400);
    c_yields_phi[l]->Divide(4,6,0,0);

    c_PbPb_pp_phi[l] = new TCanvas(Form("PbPb_minus_pp_phi%d",l),"",10,10,1500,2400);
    c_PbPb_pp_phi[l]->Divide(4,6,0,0);
    
    c_yields_eta[l] = new TCanvas(Form("yields_eta%d",l),"",10,10,1500,2400);
    c_yields_eta[l]->Divide(4,6,0,0);

    c_PbPb_pp_eta[l] = new TCanvas(Form("PbPb_minus_pp_eta%d",l),"",10,10,1500,2400);
    c_PbPb_pp_eta[l]->Divide(4,6,0,0);
 

  
    for(int i = 0; i<6; i++){

      for(int j = 0; j<4; j++){


	for(int g = 0; g<2; g++){


    
	  in_name = make_name("Yield_BkgSub_",g,i,j,l,pTlabel,centlabel,Ajlabel);

	  cout<<in_name<<endl;
	  result[g][i][j][l] = (TH2D*)f_in->Get(in_name)->Clone(in_name);

	  TString proj_name_phi = in_name; 
	  proj_name_phi.ReplaceAll("Yield_BkgSub","Proj_dPhi");

	  lbin = result[g][i][j][l]->GetXaxis()->FindBin(-2.5+.0001);
	  rbin = result[g][i][j][l]->GetXaxis()->FindBin(2.5-.0001);
	  
	  signal_dPhi[g][i][j][l] = (TH1D*)result[g][i][j][l]->ProjectionY(proj_name_phi,lbin,rbin);

	  TString leading_name = in_name; 
	  leading_name.ReplaceAll("Yield_BkgSub","Leading_dPhi");

	  TString subleading_name = in_name; 
	  subleading_name.ReplaceAll("Yield_BkgSub","SubLeading_dPhi");


	  int nbins_half = result[g][i][j][l]->GetNbinsY()/2;

	  signal_dPhi[g+4][i][j][l] = new TH1D(leading_name,"",nbins_half,-TMath::Pi()/2.,TMath::Pi()/2.);
	  signal_dPhi[g+2][i][j][l] = new TH1D(subleading_name,"",nbins_half,-TMath::Pi()/2.,TMath::Pi()/2.);

		  
	    
	  for(int k = 1; k<nbins_half+1; k++){

	    bc = (signal_dPhi[g][i][j][l]->GetBinContent(k)+signal_dPhi[g][i][j][l]->GetBinContent(nbins_half+1-k))/2.;
	    err = TMath::Sqrt(signal_dPhi[g][i][j][l]->GetBinError(k)*signal_dPhi[g][i][j][l]->GetBinError(k)+signal_dPhi[g][i][j][l]->GetBinError(nbins_half+1-k)*signal_dPhi[g][i][j][l]->GetBinError(nbins_half+1-k))/2.;
	    signal_dPhi[g+4][i][j][l]->SetBinContent(k,bc);
	    signal_dPhi[g+4][i][j][l]->SetBinContent(nbins_half+1-k,bc);
	    signal_dPhi[g+4][i][j][l]->SetBinError(k,err);
	    signal_dPhi[g+4][i][j][l]->SetBinError(nbins_half+1-k,err);

	    bc = (signal_dPhi[g][i][j][l]->GetBinContent(k+nbins_half)+signal_dPhi[g][i][j][l]->GetBinContent(2*nbins_half+1-k))/2.;
	    err = TMath::Sqrt(signal_dPhi[g][i][j][l]->GetBinError(k+nbins_half)*signal_dPhi[g][i][j][l]->GetBinError(k+nbins_half)+signal_dPhi[g][i][j][l]->GetBinError(2*nbins_half+1-k)*signal_dPhi[g][i][j][l]->GetBinError(2*nbins_half+1-k))/2.;
	    signal_dPhi[g+2][i][j][l]->SetBinContent(k,bc);
	    signal_dPhi[g+2][i][j][l]->SetBinContent(nbins_half+1-k,bc);
	    signal_dPhi[g+2][i][j][l]->SetBinError(k,err);
	    signal_dPhi[g+2][i][j][l]->SetBinError(nbins_half+1-k,err);
	  
	  }


	  if( signal_dPhi[g+4][i][j][l]->GetBinWidth(1)!= signal_dPhi[g][i][j][l]->GetBinWidth(1))cout<<"Widths do not match!"<<endl;
	  
	  bin_width_phi =  signal_dPhi[g+4][i][j][l]->GetBinWidth(1);

	  signal_dPhi[g+4][i][j][l]->Scale(1./bin_width_phi);
	  signal_dPhi[g+2][i][j][l]->Scale(1./bin_width_phi);

      
	  signal_dPhi_rebin[g+4][i][j][l] = Rebin_dPhi(signal_dPhi[g+4][i][j][l]);
	  signal_dPhi_rebin[g+2][i][j][l] = Rebin_dPhi(signal_dPhi[g+2][i][j][l]);

	 
	  TString proj_name_eta = in_name; 
	  proj_name_eta.ReplaceAll("Yield_BkgSub","Leading_dEta");

	  lbin = result[g][i][j][l]->GetYaxis()->FindBin(-TMath::Pi()/2.+.0001);
	  rbin = result[g][i][j][l]->GetYaxis()->FindBin(TMath::Pi()/2.-.0001);
	  
	  signal_dEta[g+4][i][j][l] = (TH1D*)result[g][i][j][l]->ProjectionX(proj_name_eta,lbin,rbin);

	  proj_name_eta.ReplaceAll("Leading","SubLeading");

	  lbin = result[g][i][j][l]->GetYaxis()->FindBin(TMath::Pi()/2.+.0001);
	  rbin = result[g][i][j][l]->GetYaxis()->FindBin(3.*TMath::Pi()/2.-.0001);
	  
	  signal_dEta[g+2][i][j][l] = (TH1D*)result[g][i][j][l]->ProjectionX(proj_name_eta,lbin,rbin);


	  bin_width_eta =  signal_dEta[g+4][i][j][l]->GetBinWidth(1);

	  signal_dEta[g+4][i][j][l]->Scale(1./bin_width_eta);
	  signal_dEta[g+2][i][j][l]->Scale(1./bin_width_eta);



	  int nbins_eta = result[g][i][j][l]->GetNbinsX();
    
	  for(int k = 1; k<nbins_eta+1; k++){

	    bc = (signal_dEta[g+4][i][j][l]->GetBinContent(k)+signal_dEta[g+4][i][j][l]->GetBinContent(nbins_eta+1-k))/2.;
	    err = TMath::Sqrt(signal_dEta[g+4][i][j][l]->GetBinError(k)*signal_dEta[g+4][i][j][l]->GetBinError(k)+signal_dEta[g+4][i][j][l]->GetBinError(nbins_eta+1-k)*signal_dEta[g+4][i][j][l]->GetBinError(nbins_eta+1-k))/2.;
	    signal_dEta[g+4][i][j][l]->SetBinContent(k,bc);
	    signal_dEta[g+4][i][j][l]->SetBinContent(nbins_eta+1-k,bc);
	    signal_dEta[g+4][i][j][l]->SetBinError(k,err);
	    signal_dEta[g+4][i][j][l]->SetBinError(nbins_eta+1-k,err);



	    bc = (signal_dEta[g+2][i][j][l]->GetBinContent(k)+signal_dEta[g+2][i][j][l]->GetBinContent(nbins_eta+1-k))/2.;
	    err = TMath::Sqrt(signal_dEta[g+2][i][j][l]->GetBinError(k)*signal_dEta[g+2][i][j][l]->GetBinError(k)+signal_dEta[g+2][i][j][l]->GetBinError(nbins_eta+1-k)*signal_dEta[g+2][i][j][l]->GetBinError(nbins_eta+1-k))/2.;
	    signal_dEta[g+2][i][j][l]->SetBinContent(k,bc);
	    signal_dEta[g+2][i][j][l]->SetBinContent(nbins_eta+1-k,bc);
	    signal_dEta[g+2][i][j][l]->SetBinError(k,err);
	    signal_dEta[g+2][i][j][l]->SetBinError(nbins_eta+1-k,err);



	  }



      
	  signal_dEta_rebin[g+4][i][j][l] = Rebin_dEta(signal_dEta[g+4][i][j][l]);
	  signal_dEta_rebin[g+2][i][j][l] = Rebin_dEta(signal_dEta[g+2][i][j][l]);


	  me_err[g+2][i][j][l] = max(  signal_dEta_rebin[g+2][i][j][l] ->GetBinContent(  signal_dEta_rebin[g+2][i][j][l]->FindBin(-1.51))- signal_dEta_rebin[g+2][i][j][l] ->GetBinContent(  signal_dEta_rebin[g+2][i][j][l]->FindBin(-2.01)),signal_dEta_rebin[g+2][i][j][l] ->GetBinContent(  signal_dEta_rebin[g+2][i][j][l]->FindBin(1.51))- signal_dEta_rebin[g+2][i][j][l] ->GetBinContent(  signal_dEta_rebin[g+2][i][j][l]->FindBin(2.01)));

	  me_err[g+4][i][j][l] = max(  signal_dEta_rebin[g+4][i][j][l] ->GetBinContent(  signal_dEta_rebin[g+4][i][j][l]->FindBin(-1.51))- signal_dEta_rebin[g+4][i][j][l] ->GetBinContent(  signal_dEta_rebin[g+4][i][j][l]->FindBin(-2.01)),signal_dEta_rebin[g+4][i][j][l] ->GetBinContent(  signal_dEta_rebin[g+4][i][j][l]->FindBin(1.51))- signal_dEta_rebin[g+4][i][j][l] ->GetBinContent(  signal_dEta_rebin[g+4][i][j][l]->FindBin(2.01)));
	
	  //------------------
	  //Apply corrections
	  //-------------------
	 
	  TString jff_name; 

	  if(g%2!=0){

	    if(j==0){
	      jff_name = make_name("JFF_Residual_Eta_",g+8,i,3,l,pTlabel,centlabel,Ajlabel);
	      jff_residual_dEta[g+2][i][j][l] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);

	      jff_name.ReplaceAll("Eta","Phi");
	      jff_residual_dPhi[g+2][i][j][l] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);
	 
	      jff_name = make_name("JFF_Residual_Eta_",g+10,i,3,l,pTlabel,centlabel,Ajlabel);
	      jff_residual_dEta[g+4][i][j][l] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);

	      jff_name.ReplaceAll("Eta","Phi");
	      jff_residual_dPhi[g+4][i][j][l] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);
	    }

	    signal_dEta_rebin[g+2][i][j][l]->Add(jff_residual_dEta[g+2][i][0][l],-1.);
	    signal_dPhi_rebin[g+2][i][j][l]->Add(jff_residual_dPhi[g+2][i][0][l],-1.);
	 
	    signal_dEta_rebin[g+4][i][j][l]->Add(jff_residual_dEta[g+4][i][0][l],-1.);
	    signal_dPhi_rebin[g+4][i][j][l]->Add(jff_residual_dPhi[g+4][i][0][l],-1.);


	  }else{

	    TString spill_name = make_name("SpillOvers_Eta_pTweighted_",g+8,i,j,2,pTlabel,centlabel,Ajlabel);
	    if(is_number) spill_name = make_name("SpillOvers_Eta_",g+8,i,j,2,pTlabel,centlabel,Ajlabel);
	    
	    spill_name.ReplaceAll("Pt100_Pt300_","");
	    spill_over_dEta[g+2][i][j][l] = (TH1D*)f_spillover->Get(spill_name)->Clone(spill_name);
	    spill_name.ReplaceAll("Eta","Phi");
	    spill_over_dPhi[g+2][i][j][l] = (TH1D*)f_spillover->Get(spill_name)->Clone(spill_name);
	 

	    spill_name = make_name("SpillOvers_Eta_pTweighted_",g+10,i,j,l,pTlabel,centlabel,Ajlabel);
	    if(is_number)spill_name = make_name("SpillOvers_Eta_",g+10,i,j,l,pTlabel,centlabel,Ajlabel);
	    spill_name.ReplaceAll("Pt100_Pt300_","");

	    spill_over_dEta[g+4][i][j][l] = (TH1D*)f_spillover->Get(spill_name)->Clone(spill_name);
	  
	    spill_name.ReplaceAll("Eta","Phi");
	    spill_over_dPhi[g+4][i][j][l] = (TH1D*)f_spillover->Get(spill_name)->Clone(spill_name);

	
	    signal_dEta_rebin[g+2][i][j][l]->Add(spill_over_dEta[g+2][i][j][l],-1.);
	    signal_dPhi_rebin[g+2][i][j][l]->Add(spill_over_dPhi[g+2][i][j][l],-1.);
	 
	    signal_dEta_rebin[g+4][i][j][l]->Add(spill_over_dEta[g+4][i][j][l],-1.);
	    signal_dPhi_rebin[g+4][i][j][l]->Add(spill_over_dPhi[g+4][i][j][l],-1.);


	    jff_name = make_name("JFF_Residual_Eta_",g+8,i,j,l,pTlabel,centlabel,Ajlabel);
	    if(i>4)	    jff_residual_dEta[g+2][i][j][l] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);
	    else   jff_residual_dEta[g+2][i][j][l] = (TH1D*)f_jff_gen->Get(jff_name)->Clone(jff_name);

	    jff_name.ReplaceAll("Eta","Phi");
	    if(i>4)    jff_residual_dPhi[g+2][i][j][l] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);
	    else     jff_residual_dPhi[g+2][i][j][l] = (TH1D*)f_jff_gen->Get(jff_name)->Clone(jff_name);
	    
	    jff_name = make_name("JFF_Residual_Eta_",g+10,i,j,l,pTlabel,centlabel,Ajlabel);
	    if(i>4)    jff_residual_dEta[g+4][i][j][l] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);
	    else jff_residual_dEta[g+4][i][j][l] = (TH1D*)f_jff_gen->Get(jff_name)->Clone(jff_name);

	    jff_name.ReplaceAll("Eta","Phi");
	    if(i>4)   jff_residual_dPhi[g+4][i][j][l] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);
	    else    jff_residual_dPhi[g+4][i][j][l] = (TH1D*)f_jff_gen->Get(jff_name)->Clone(jff_name);
	 
	    signal_dEta_rebin[g+2][i][j][l]->Add(jff_residual_dEta[g+2][i][j][l],-1.);
	    signal_dPhi_rebin[g+2][i][j][l]->Add(jff_residual_dPhi[g+2][i][j][l],-1.);
	  
	    signal_dEta_rebin[g+4][i][j][l]->Add(jff_residual_dEta[g+4][i][j][l],-1.);
	    signal_dPhi_rebin[g+4][i][j][l]->Add(jff_residual_dPhi[g+4][i][j][l],-1.);

	    /*    
	    if(j==0){
	      jff_name = make_name("JFF_Residual_Eta_",g+9,i,3,l,pTlabel,centlabel,Ajlabel);
	      jff_residual_dEta[g+2][i][j][l] = (TH1D*)f_jff->Get(jff_name)->Clone(jff_name);

	      jff_name.ReplaceAll("Eta","Phi");
	      jff_residual_dPhi[g+2][i][j][l] = (TH1D*)f_jff->Get(jff_name)->Clone(jff_name);
	 
	      jff_name = make_name("JFF_Residual_Eta_",g+11,i,3,l,pTlabel,centlabel,Ajlabel);
	      jff_residual_dEta[g+4][i][j][l] = (TH1D*)f_jff->Get(jff_name)->Clone(jff_name);

	      jff_name.ReplaceAll("Eta","Phi");
	      jff_residual_dPhi[g+4][i][j][l] = (TH1D*)f_jff->Get(jff_name)->Clone(jff_name);
	    }
	    
	    signal_dEta_rebin[g+2][i][j][l]->Add(jff_residual_dEta[g+2][i][0][l],-1.);
	    signal_dPhi_rebin[g+2][i][j][l]->Add(jff_residual_dPhi[g+2][i][0][l],-1.);
	 
	    signal_dEta_rebin[g+4][i][j][l]->Add(jff_residual_dEta[g+4][i][0][l],-1.);
	    signal_dPhi_rebin[g+4][i][j][l]->Add(jff_residual_dPhi[g+4][i][0][l],-1.);
	    */

	    
	  }

	
	  cout<<"and here"<<endl;

	  //---------------------------------------------------
	  //---------------------------------------------------

	  signal_dEta_rebin[g+4][i][j][l]->SetAxisRange(-2.49,2.49);
	  signal_dEta_rebin[g+2][i][j][l]->SetAxisRange(-2.49,2.49);
	  
	  signal_dPhi_rebin[g+4][i][j][l]->SetAxisRange(-2.49,2.49);
	  signal_dPhi_rebin[g+2][i][j][l]->SetAxisRange(-2.49,2.49);


	  signal_dEta_rebin[g+4][i][j][l]->Scale(-1.);

	  signal_dEta_rebin[g+4][i][j][l]->SetMarkerSize(1);
	  signal_dEta_rebin[g+4][i][j][l]->SetMarkerColor(kOrange-2);
	  signal_dEta_rebin[g+4][i][j][l]->SetLineColor(kOrange-2);

	  signal_dEta_rebin[g+2][i][j][l]->SetMarkerSize(1);
	  signal_dEta_rebin[g+2][i][j][l]->SetMarkerColor(kGreen+3);
	  signal_dEta_rebin[g+2][i][j][l]->SetLineColor(kGreen+3);

	  signal_dPhi_rebin[g+4][i][j][l]->Scale(-1.);

	  signal_dPhi_rebin[g+4][i][j][l]->SetMarkerSize(1);
	  signal_dPhi_rebin[g+4][i][j][l]->SetMarkerColor(kOrange-2);
	  signal_dPhi_rebin[g+4][i][j][l]->SetLineColor(kOrange-2);

	  signal_dPhi_rebin[g+2][i][j][l]->SetMarkerSize(1);
	  signal_dPhi_rebin[g+2][i][j][l]->SetMarkerColor(kGreen+3);
	  signal_dPhi_rebin[g+2][i][j][l]->SetLineColor(kGreen+3);





	}//g to get hists

	switch(i){
	case 0: 
	  signal_max = 20.;
	  signal_min = -20.;
	  diff_max = 12.;
	  diff_min = -12.;
	  break;
	case 1: 
	  signal_max = 20.;
	  signal_min = -20.;
	  diff_max = 12.;
	  diff_min = -12.;
	  break;
	case 2: 
	  signal_max = 20.;
	  signal_min = -20.;
	  diff_max = 12.;
	  diff_min = -12.;
	  break;
	case 3:
	  signal_max = 20.;
	  signal_min = -20.; 
	  diff_max = 12.;
	  diff_min = -12.;
	  break;
	case 4: 
	  signal_max = 30.;
	  signal_min = -30.;
	  diff_max = 12.;
	  diff_min = -12.;
	  break;
	case 5: 
	  signal_max = 500.;
	  signal_min = -500.;
	  diff_max = 100.;
	  diff_min = -100.;
	  break;
	default: 
	  break;
	}






	c_yields_phi[l]->cd(4*i+j+1);

	signal_dPhi_rebin[2][i][j][l]->SetMarkerStyle(20);
	signal_dPhi_rebin[3][i][j][l]->SetMarkerStyle(4);
	signal_dPhi_rebin[4][i][j][l]->SetMarkerStyle(20);
	signal_dPhi_rebin[5][i][j][l]->SetMarkerStyle(4);

	signal_dPhi_rebin[2][i][j][l]->SetMinimum(signal_min);
	signal_dPhi_rebin[2][i][j][l]->SetMaximum(signal_max);


	if(j==0){
	  signal_dPhi_rebin[2][i][j][l]->GetYaxis()->SetLabelSize(0.06);
	  signal_dPhi_rebin[2][i][j][l]->GetYaxis()->SetTitleSize(0.06);
	  signal_dPhi_rebin[2][i][j][l]->GetYaxis()->SetTitle("1/N_{evt} 1/d#Delta#phi");
	}else{
	  signal_dPhi_rebin[2][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	}

	if(i==4){
	  signal_dPhi_rebin[2][i][j][l]->GetXaxis()->SetLabelSize(0.06);
	  signal_dPhi_rebin[2][i][j][l]->GetXaxis()->SetTitleSize(0.06);
	  signal_dPhi_rebin[2][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	  signal_dPhi_rebin[2][i][j][l]->GetXaxis()->CenterTitle();

	}
	signal_dPhi_rebin[2][i][j][l]->Draw();
	signal_dPhi_rebin[3][i][j][l]->Draw("same");
	signal_dPhi_rebin[4][i][j][l]->Draw("same");
	signal_dPhi_rebin[5][i][j][l]->Draw("same");

	TString diff_name_pbpb = make_name("Diff_",0,i,j,l,pTlabel, centlabel, Ajlabel);  	diff_name_pbpb+="_Rebin";


	background_diff_rebin[0][i][j][l] = (TH1D*)f_in->Get(diff_name_pbpb)->Clone(diff_name_pbpb);

	TString diff_name_pp = make_name("Diff_",1,i,j,l,pTlabel, centlabel, Ajlabel);  	diff_name_pp+="_Rebin";

	background_diff_rebin[1][i][j][l] = (TH1D*)f_in->Get(diff_name_pp)->Clone(diff_name_pp);


	background_diff_rebin[0][i][j][l]->Scale(5./2); //projection range is 5. here, sideband range was 2. ***BE CAREFUL**
	background_diff_rebin[1][i][j][l]->Scale(5./2); //projection range is 5. here, sideband range was 2. ***BE CAREFUL**


	TString syst_name_pbpb = make_name("Syst_",0,i,j,l,pTlabel, centlabel, Ajlabel);  


	background_syst_rebin[0][i][j][l] = (TH1D*)f_in->Get(syst_name_pbpb)->Clone(syst_name_pbpb);

	TString syst_name_pp = make_name("Syst_",1,i,j,l,pTlabel, centlabel, Ajlabel);  

	background_syst_rebin[1][i][j][l] = (TH1D*)f_in->Get(syst_name_pp)->Clone(syst_name_pp);


	background_syst_rebin[0][i][j][l]->Scale(5./2); //projection range is 5. here, sideband range was 2. ***BE CAREFUL**
	background_syst_rebin[1][i][j][l]->Scale(5./2); //projection range is 5. here, sideband range was 2. ***BE CAREFUL**


	TString syst_name = make_name("SubLeading_dEta_Syst_",0,i,j,l,pTlabel,centlabel,Ajlabel);
	signal_dEta_syst[2][i][j][l] = (TH1D*)signal_dEta_rebin[2][i][j][l]->Clone(syst_name);

	syst_name.ReplaceAll("SubLeading","Leading");
	signal_dEta_syst[4][i][j][l] = (TH1D*)signal_dEta_rebin[4][i][j][l]->Clone(syst_name);
	  
	syst_name.ReplaceAll("Eta","Phi");
	signal_dPhi_syst[4][i][j][l] = (TH1D*)signal_dPhi_rebin[4][i][j][l]->Clone(syst_name);

	syst_name.ReplaceAll("Leading","SubLeading");
	signal_dPhi_syst[2][i][j][l] = (TH1D*)signal_dPhi_rebin[2][i][j][l]->Clone(syst_name);


	syst_name = make_name("SubLeading_dEta_Syst_",1,i,j,l,pTlabel,centlabel,Ajlabel);
	signal_dEta_syst[3][i][j][l] = (TH1D*)signal_dEta_rebin[3][i][j][l]->Clone(syst_name);

	syst_name.ReplaceAll("SubLeading","Leading");
	signal_dEta_syst[5][i][j][l] = (TH1D*)signal_dEta_rebin[5][i][j][l]->Clone(syst_name);
	  
	syst_name.ReplaceAll("Eta","Phi");
	signal_dPhi_syst[5][i][j][l] = (TH1D*)signal_dPhi_rebin[5][i][j][l]->Clone(syst_name);

	syst_name.ReplaceAll("Leading","SubLeading");
	signal_dPhi_syst[3][i][j][l] = (TH1D*)signal_dPhi_rebin[3][i][j][l]->Clone(syst_name);



	if(use_highpT_bin&&i==5){
	    
	  for(int k = 0; k<5; k++){
	    signal_dEta_syst[2][i][j][l]->Add(signal_dEta_syst[2][k][j][l]);
	    signal_dPhi_syst[2][i][j][l]->Add(signal_dPhi_syst[2][k][j][l]);

	    signal_dEta_syst[4][i][j][l]->Add(signal_dEta_syst[4][k][j][l]);
	    signal_dPhi_syst[4][i][j][l]->Add(signal_dPhi_syst[4][k][j][l]);

	    signal_dEta_syst[3][i][j][l]->Add(signal_dEta_syst[3][k][j][l]);
	    signal_dPhi_syst[3][i][j][l]->Add(signal_dPhi_syst[3][k][j][l]);

	    signal_dEta_syst[5][i][j][l]->Add(signal_dEta_syst[5][k][j][l]);
	    signal_dPhi_syst[5][i][j][l]->Add(signal_dPhi_syst[5][k][j][l]);


	  }
	}
	    

	if(!use_highpT_bin&&i==4){
	    
	  for(int k = 0; k<4; k++){
	    signal_dEta_syst[2][i][j][l]->Add(signal_dEta_syst[2][k][j][l]);
	    signal_dPhi_syst[2][i][j][l]->Add(signal_dPhi_syst[2][k][j][l]);

	    signal_dEta_syst[4][i][j][l]->Add(signal_dEta_syst[4][k][j][l]);
	    signal_dPhi_syst[4][i][j][l]->Add(signal_dPhi_syst[4][k][j][l]);

	    signal_dEta_syst[3][i][j][l]->Add(signal_dEta_syst[3][k][j][l]);
	    signal_dPhi_syst[3][i][j][l]->Add(signal_dPhi_syst[3][k][j][l]);

	    signal_dEta_syst[5][i][j][l]->Add(signal_dEta_syst[5][k][j][l]);
	    signal_dPhi_syst[5][i][j][l]->Add(signal_dPhi_syst[5][k][j][l]);


	  }
	    
	}
	if((use_highpT_bin&&i==5)||(!use_highpT_bin&&i==4)){
	  for(int k = 1; k<signal_dPhi_syst[2][i][j][l]->GetNbinsX()+1; k++){

	    bc =  signal_dPhi_syst[2][i][j][l]->GetBinContent(k);
	    err = TMath::Sqrt(bc*rel_err*bc*rel_err+background_syst_rebin[0][i][j][l]->GetBinError(1)*background_syst_rebin[0][i][j][l]->GetBinError(1)+me_err[2][i][j][l]*me_err[2][i][j][l]+spill_over_dPhi[2][0][j][l]->GetBinContent(k)*spill_over_dPhi[2][0][j][l]->GetBinContent(k)/4.);
	    signal_dPhi_syst[2][i][j][l]->SetBinError(k,err);

	    bc =  signal_dPhi_syst[4][i][j][l]->GetBinContent(k);
	    err = TMath::Sqrt(bc*rel_err*bc*rel_err+background_syst_rebin[0][i][j][l]->GetBinError(1)*background_syst_rebin[0][i][j][l]->GetBinError(1)+me_err[4][i][j][l]*me_err[4][i][j][l]+spill_over_dPhi[4][0][j][l]->GetBinContent(k)*spill_over_dPhi[4][0][j][l]->GetBinContent(k)/4.);
	    signal_dPhi_syst[4][i][j][l]->SetBinError(k,err);
	  }
	    
	  for(int k = 1; k<signal_dPhi_syst[3][i][j][l]->GetNbinsX()+1; k++){
	    bc =  signal_dPhi_syst[3][i][j][l]->GetBinContent(k);
	    err =TMath::Sqrt(rel_err*bc*rel_err*bc+background_syst_rebin[1][i][j][l]->GetBinError(1)*background_syst_rebin[1][i][j][l]->GetBinError(1)+me_err[3][i][j][l]*me_err[3][i][j][l]);
	    signal_dPhi_syst[3][i][j][l]->SetBinError(k,err);

	    bc =  signal_dPhi_syst[5][i][j][l]->GetBinContent(k);
	    err =TMath::Sqrt(rel_err*bc*rel_err*bc+background_syst_rebin[1][i][j][l]->GetBinError(1)*background_syst_rebin[1][i][j][l]->GetBinError(1)+me_err[5][i][j][l]*me_err[5][i][j][l]);
	    signal_dPhi_syst[5][i][j][l]->SetBinError(k,err);
	  }
	  
	  for(int k = 1; k<background_syst_rebin[0][i][j][l]->GetNbinsX()+1; k++){
	    bc = background_syst_rebin[0][i][j][l]->GetBinContent(k);
	    err =TMath::Sqrt(rel_err*bc*rel_err*bc+background_syst_rebin[0][i][j][l]->GetBinError(1)*background_syst_rebin[0][i][j][l]->GetBinError(1)+me_err[2][i][j][l]*me_err[2][i][j][l]+me_err[4][i][j][l]*me_err[4][i][j][l]);
	    background_syst_rebin[0][i][j][l]->SetBinError(k,err);
	

	    bc = background_syst_rebin[1][i][j][l]->GetBinContent(k);
	    err =TMath::Sqrt(rel_err*bc*rel_err*bc+background_syst_rebin[1][i][j][l]->GetBinError(1)*background_syst_rebin[1][i][j][l]->GetBinError(1)+me_err[3][i][j][l]*me_err[3][i][j][l]+me_err[5][i][j][l]*me_err[5][i][j][l]);
	    background_syst_rebin[1][i][j][l]->SetBinError(k,err);
	
	  }

	}
      

	//------------
	// Drawing
	//------------



	background_diff_rebin[0][i][j][l]->SetMarkerStyle(20);
	background_diff_rebin[0][i][j][l]->SetMarkerSize(1);
	background_diff_rebin[0][i][j][l]->SetLineColor(kRed);
	background_diff_rebin[0][i][j][l]->SetMarkerColor(kRed);

	background_diff_rebin[1][i][j][l]->SetMarkerStyle(4);
	background_diff_rebin[1][i][j][l]->SetMarkerSize(1);
	background_diff_rebin[1][i][j][l]->SetLineColor(kRed);
	background_diff_rebin[1][i][j][l]->SetMarkerColor(kRed);


	TLine *zero = new TLine(-1.5,0.,1.5,0.);
	zero->SetLineStyle(2);
	zero->Draw();

	TPaveText *labels;
	if(j==0){
	  labels = new TPaveText(0.18,0.75,0.45,0.95,"NDC");
	}else{
	  labels = new TPaveText(0.05,0.75,0.45,0.95,"NDC");
	}  
	labels->SetName("labels");
	labels->SetFillColor(0);
	labels->SetLineColor(0);
	labels->SetTextAlign(11);
	labels->AddText(centlabel);
	labels->AddText(Ajlabel);
	labels->AddText(pTlabel);
	labels->SetTextSize(ts2-0.01);
	labels->Draw("same");

	if(i==0&&j==0){
	  TLegend *legend = new TLegend(0.6,0.75,0.95,0.95);

	  legend->AddEntry(signal_dPhi_rebin[2][i][j][l],"PbPb Sub.");
	  legend->AddEntry(signal_dPhi_rebin[3][i][j][l],"pp Sub.");
	  legend->AddEntry(signal_dPhi_rebin[4][i][j][l],"PbPb Lead.");
	  legend->AddEntry(signal_dPhi_rebin[5][i][j][l],"pp Lead.");
	  //	legend->AddEntry(background_diff_rebin[0][i][j][l],"PbPbBkg S-L");
	  //legend->AddEntry(background_diff_rebin[1][i][j][l],"ppBkg S-L");

	  legend->SetTextSize(ts2-0.01);
	  legend->SetLineColor(0);
	  legend->Draw();
	}

	c_PbPb_pp_phi[l]->cd(4*i+j+1);

	TString sub_name_pbpb_pp_phi = make_name("Sub_PbPb_minus_pp_dPhi_",0,i,j,l,pTlabel, centlabel, Ajlabel);  
	signal_dPhi_PbPb_pp[2][i][j][l] = (TH1D*)signal_dPhi_rebin[2][i][j][l]->Clone(sub_name_pbpb_pp_phi);
	signal_dPhi_PbPb_pp[2][i][j][l]->Add(signal_dPhi_rebin[3][i][j][l],-1.);

	TString lead_name_pbpb_pp_phi = make_name("Lead_PbPb_minus_pp_dPhi",0,i,j,l,pTlabel, centlabel, Ajlabel);  
	signal_dPhi_PbPb_pp[4][i][j][l] = (TH1D*)signal_dPhi_rebin[4][i][j][l]->Clone(lead_name_pbpb_pp_phi);
	signal_dPhi_PbPb_pp[4][i][j][l]->Add(signal_dPhi_rebin[5][i][j][l],-1.);

	TString diff_name_pbpb_pp_phi = make_name("Diff_PbPb_minus_pp_dPhi",0,i,j,l,pTlabel, centlabel, Ajlabel);  
	background_diff_PbPb_pp[0][i][j][l] = (TH1D*)background_diff_rebin[0][i][j][l]->Clone(diff_name_pbpb_pp_phi);
	background_diff_PbPb_pp[0][i][j][l]->Add(background_diff_rebin[1][i][j][l],-1.);	


	signal_dPhi_PbPb_pp[2][i][j][l]->SetMaximum(diff_max);
	signal_dPhi_PbPb_pp[2][i][j][l]->SetMinimum(diff_min);

	signal_dPhi_PbPb_pp[2][i][j][l]->Draw();
	signal_dPhi_PbPb_pp[4][i][j][l]->Draw("same");


	background_diff_PbPb_pp[0][i][j][l]->Draw("same");
	zero->Draw();
	labels->SetLineColor(0);
	labels->Draw("same");

	if(j==0&&i==0){
	  TLegend *legend_diff = new TLegend(0.6,0.75,0.95,0.95);

	  legend_diff->AddEntry(signal_dPhi_PbPb_pp[2][i][j][l],"Sub. PbPb-pp");
	  legend_diff->AddEntry(signal_dPhi_PbPb_pp[4][i][j][l],"Lead. PbPb-pp");
	  legend_diff->AddEntry(background_diff_PbPb_pp[0][i][j][l],"Bkg S-L PbPb-pp");

	  legend_diff->SetLineColor(0);
	  legend_diff->SetTextSize(ts2-0.01);
	  legend_diff->Draw();

	}

	//_____________
	//  Draw dEta
	//-------------


	c_yields_eta[l]->cd(4*i+j+1);

	signal_dEta_rebin[2][i][j][l]->SetMarkerStyle(20);
	signal_dEta_rebin[3][i][j][l]->SetMarkerStyle(4);
	signal_dEta_rebin[4][i][j][l]->SetMarkerStyle(20);
	signal_dEta_rebin[5][i][j][l]->SetMarkerStyle(4);

	signal_dEta_rebin[2][i][j][l]->SetMinimum(signal_min);
	signal_dEta_rebin[2][i][j][l]->SetMaximum(signal_max);

	if(j==0){
	  signal_dEta_rebin[2][i][j][l]->GetYaxis()->SetLabelSize(0.06);
	  signal_dEta_rebin[2][i][j][l]->GetYaxis()->SetTitleSize(0.06);
	  signal_dEta_rebin[2][i][j][l]->GetYaxis()->SetTitle("1/N_{evt} 1/d#Delta#eta");
	}else{
	  signal_dEta_rebin[2][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	}

	if(i==4){
	  signal_dEta_rebin[2][i][j][l]->GetXaxis()->SetLabelSize(0.06);
	  signal_dEta_rebin[2][i][j][l]->GetXaxis()->SetTitleSize(0.06);
	  signal_dEta_rebin[2][i][j][l]->GetXaxis()->SetTitle("#Delta#eta");
	  signal_dEta_rebin[2][i][j][l]->GetXaxis()->CenterTitle();

	}
	signal_dEta_rebin[2][i][j][l]->Draw();
	signal_dEta_rebin[3][i][j][l]->Draw("same");
	signal_dEta_rebin[4][i][j][l]->Draw("same");
	signal_dEta_rebin[5][i][j][l]->Draw("same");
	

	zero->Draw();
	labels->Draw("same");

	if(i==0&&j==0){
	  TLegend *legend = new TLegend(0.6,0.75,0.95,0.95);

	  legend->AddEntry(signal_dPhi_rebin[2][i][j][l],"PbPb Sub.");
	  legend->AddEntry(signal_dPhi_rebin[3][i][j][l],"pp Sub.");
	  legend->AddEntry(signal_dPhi_rebin[4][i][j][l],"PbPb Lead.");
	  legend->AddEntry(signal_dPhi_rebin[5][i][j][l],"pp Lead.");
	  legend->SetTextSize(ts2-0.01);
	  legend->SetLineColor(0);
	  legend->Draw();
	}

	c_PbPb_pp_eta[l]->cd(4*i+j+1);

	TString sub_name_pbpb_pp_eta = make_name("Sub_PbPb_minus_pp_dEta",0,i,j,l,pTlabel, centlabel, Ajlabel);  
	signal_dEta_PbPb_pp[2][i][j][l] = (TH1D*)signal_dEta_rebin[2][i][j][l]->Clone(sub_name_pbpb_pp_eta);
	signal_dEta_PbPb_pp[2][i][j][l]->Add(signal_dEta_rebin[3][i][j][l],-1.);

	TString lead_name_pbpb_pp_eta = make_name("Lead_PbPb_minus_pp_dEta",0,i,j,l,pTlabel, centlabel, Ajlabel);  
	signal_dEta_PbPb_pp[4][i][j][l] = (TH1D*)signal_dEta_rebin[4][i][j][l]->Clone(lead_name_pbpb_pp_eta);
	signal_dEta_PbPb_pp[4][i][j][l]->Add(signal_dEta_rebin[5][i][j][l],-1.);

	signal_dEta_PbPb_pp[2][i][j][l]->SetMaximum(diff_max);
	signal_dEta_PbPb_pp[2][i][j][l]->SetMinimum(diff_min);

	signal_dEta_PbPb_pp[2][i][j][l]->Draw();
	signal_dEta_PbPb_pp[4][i][j][l]->Draw("same");

	zero->Draw();
	labels->SetLineColor(0);
	labels->Draw("same");

	if(j==0&&i==0){
	  TLegend *legend_diff = new TLegend(0.6,0.75,0.95,0.95);

	  legend_diff->AddEntry(signal_dEta_PbPb_pp[2][i][j][l],"Sub. PbPb-pp");
	  legend_diff->AddEntry(signal_dEta_PbPb_pp[4][i][j][l],"Lead. PbPb-pp");
	  legend_diff->SetLineColor(0);
	  legend_diff->SetTextSize(ts2-0.01);
	  legend_diff->Draw();


	}

	//---------------------
	//  Write everything!
	//----------------------

	f_out->cd();

	signal_dPhi_rebin[2][i][j][l]->Write();
	signal_dPhi_rebin[4][i][j][l]->Write();
	signal_dPhi_rebin[3][i][j][l]->Write();
	signal_dPhi_rebin[5][i][j][l]->Write();
	
	background_diff_rebin[0][i][j][l]->Write();
	background_diff_rebin[1][i][j][l]->Write();

	signal_dEta_rebin[2][i][j][l]->Write();
	signal_dEta_rebin[4][i][j][l]->Write();
	signal_dEta_rebin[3][i][j][l]->Write();
	signal_dEta_rebin[5][i][j][l]->Write();


	background_syst_rebin[0][i][j][l]->Write();
	background_syst_rebin[1][i][j][l]->Write();


	  
	signal_dPhi_syst[2][i][j][l]->Write();
	signal_dPhi_syst[4][i][j][l]->Write();

	signal_dPhi_syst[3][i][j][l]->Write();
	signal_dPhi_syst[5][i][j][l]->Write();

      }
    }

    if(l<2){
      c_yields_phi[l]->SaveAs((TString)("Yields_dPhi_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));
      c_PbPb_pp_phi[l]->SaveAs((TString)("Yields_dPhi_PbPb_minus_pp_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));
      c_yields_eta[l]->SaveAs((TString)("Yields_dEta_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));
      c_PbPb_pp_eta[l]->SaveAs((TString)("Yields_dEta_PbPb_minus_pp_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));
 
    }else{
      c_yields_phi[l]->SaveAs((TString)("Yields_dPhi_AjInclusive.png"));
      c_PbPb_pp_phi[l]->SaveAs((TString)("Yields_dPhi_PbPb_minus_pp_AjInclusive.png"));
      c_yields_eta[l]->SaveAs((TString)("Yields_dEta_AjInclusive.png"));
      c_PbPb_pp_eta[l]->SaveAs((TString)("Yields_dEta_PbPb_minus_pp_AjInclusive.png"));
    }
  }
  
 

  for(int l = 0; l<3; l++){
    for(int j = 0; j<4; j++){
      for(int i = 0; i<6; i++){

	cout<<l<<" "<<j<<" "<<i<<" "<<spill_over_dPhi[4][i][j][l]->Integral("width")<<" "<<signal_dEta_rebin[4][i][j][l]->Integral("width")<<endl;
	
	//	cout<<j<<" "<<l<<" "<<i<<signal_dEta_rebin[4][i][j][l]->Integral("width")<<" "<<signal_dEta_rebin[2][i][j][l]->Integral("width")<<endl;
      }
    }
  }

  if(!is_number){
    cout<<"JFF RESIDUAL"<<endl;

    for(int g = 2; g<6; g++){

      for(int l = 0; l<3; l++){
	for(int j = 0; j<4; j++){
	  for(int i = 0; i<6; i++){
       
	    if(g%2!=0&&j!=0)continue;
	    //	    cout<<l<<" "<<j<<" "<<i<<" "<<jff_residual_dPhi[g][i][j][l]->Integral("width")<<" "<<signal_dEta_rebin[g][i][j][l]->Integral("width")<<endl;
	
	    //	cout<<j<<" "<<l<<" "<<i<<signal_dEta_rebin[4][i][j][l]->Integral("width")<<" "<<signal_dEta_rebin[2][i][j][l]->Integral("width")<<endl;

	    cout<<g<<" "<<l<<" "<<j<<" "<<i<<" "<<me_err[g][i][j][l]<<endl;
	  }
	}
      }
    }
  }


  return 0;
}
