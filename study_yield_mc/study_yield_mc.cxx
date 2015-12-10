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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TExec.h"
#include "TLatex.h"
#include "THStack.h"
#include <iostream>
#include <vector>
#include <fstream>

#include "../JetTrack2015_functions.h"


using namespace std;

Int_t study_yield_mc(bool use_highpT_bin = kFALSE){

 
  gROOT->ForceStyle();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);
  gStyle->SetHatchesLineWidth(2);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);


 
  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RecoReco_Minus_GenGen, n_data_mc_types};

  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RecoReco_Minus_GenGen"};
      


  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 6;
  const int nAjBins = 3;



  float CBins[nCBins+1] = {0, 10, 30, 50, 100};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%", "Cent. 30-50%","Cent. 50-100%"};

 
  float TrkPtBins[nTrkPtBins+1] = {0.5,1, 2, 3, 4, 8, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};

  float AjBins[nAjBins+1] = {0,0.22,0.75};
  TString AjBin_labels[nAjBins+1] = {"A_{J}<0.22","A_{J}>0.22","A_{J} Inclusive"};

  TFile *fin[12][n_data_mc_types];
  TFile *fout[12][n_data_mc_types];
  TFile *f_in_syst = new TFile("../final_plots/Systematic_Error.root");  


  TFile *f_spillover = new TFile("../spill_over/Dijet_SpillOvers.root");
  TFile  *f_jff_reco = new TFile("../jff_residual/All_JFFResiduals_RecoReco.root");
  TFile  *f_jff_gen = new TFile("../jff_residual/All_JFFResiduals_RecoGen.root");
  TFile *f_err;

  TH1D *spill_over_dPhi[12][6][4][3];
  TH1D *spill_over_dEta[12][6][4][3];
  
  TH1D *jff_residual_dPhi[12][6][4][3];
  TH1D *jff_residual_dEta[12][6][4][3];
 
 
  TH2D *raw[12][6][4][3][n_data_mc_types];
  TH2D *result[12][6][4][3][n_data_mc_types];
  TH2D *result_sub_lead[12][6][4][3][n_data_mc_types];
  TH1D *result_sub_lead_proj[12][6][4][3][n_data_mc_types];
 
  TH2D* background[12][6][4][3][n_data_mc_types];
  TH1D* background_left[12][6][4][3][n_data_mc_types];
  TH1D* background_right[12][6][4][3][n_data_mc_types];
  TH1D* background_proj[12][6][4][3][n_data_mc_types];
  TH1D* background_proj_rebin[12][6][4][3][n_data_mc_types];

  TH1D *phi_proj[12][6][4][3][n_data_mc_types];
  TH1D *phi_proj_rebin[12][6][4][3][n_data_mc_types];
  TH1D *phi_proj_rebin2[12][6][4][3][n_data_mc_types];
  TH1D *eta_proj[12][6][4][3][n_data_mc_types];
  TH1D *eta_proj_rebin[12][6][4][3][n_data_mc_types];
  TH1D *eta_proj_rebin2[12][6][4][3][n_data_mc_types];
 
  TH1D *signal_dPhi_rebin[12][6][4][5][n_data_mc_types];
  TH1D *signal_dPhi_tot[12][6][4][5][n_data_mc_types];
  TH1D *signal_dPhi_noerr_up[12][6][4][5][n_data_mc_types];
  TH1D *signal_dPhi_noerr_down[12][6][4][5][n_data_mc_types];

  TH1D *sub_lead_dPhi_rebin[12][6][4][5][n_data_mc_types];
  TH1D *sub_lead_dPhi_tot[12][6][4][5][n_data_mc_types];
  TH1D *sub_lead_dPhi_noerr_up[12][6][4][5][n_data_mc_types];
  TH1D *sub_lead_dPhi_noerr_down[12][6][4][5][n_data_mc_types];

  TH1D *blank[12][6][4][5][n_data_mc_types];  

  TH1D *background_diff_rebin[12][6][4][5][n_data_mc_types];
  TH1D *background_diff_tot[12][6][4][5][n_data_mc_types];
  TH1D *background_diff_noerr_up[12][6][4][5][n_data_mc_types];
  TH1D *background_diff_noerr_down[12][6][4][5][n_data_mc_types];


  TH1D *signal_dPhi_syst[12][6][4][5];
  TH1D *sub_lead_dPhi_syst[12][6][4][5];
  TH1D *background_syst_rebin[12][6][4][5];
  

  THStack *dPhi_pTdist_up[12][4][3][n_data_mc_types];
  THStack *dPhi_Sub_Lead_NoBkg_up[12][4][3][n_data_mc_types];
  THStack *dPhi_pTdist_down[12][4][3][n_data_mc_types];
  THStack *dPhi_Sub_Lead_NoBkg_down[12][4][3][n_data_mc_types];

  THStack *Integral_up[12][4][3][n_data_mc_types];
  THStack *Integral_down[12][4][3][n_data_mc_types];


  TH1D *Integral[12][4][3][n_data_mc_types];
  TH1D *Integral_noerr_up[12][4][3][n_data_mc_types];
  TH1D *Integral_noerr_down[12][4][3][n_data_mc_types];


  TCanvas *c_jet[5][n_data_mc_types];
  TCanvas *c_nobkgsub[5][n_data_mc_types];
  TCanvas *c_longrange[5][n_data_mc_types];
  TCanvas *c_longrange_int;
  TCanvas *c_all_int;
 
  TString in_name, stack_name, tot_name, pTlabel,centlabel,Ajlabel,infile_name, datalabel,name_stem;

  TString AjBin_strs[4] = {"Aj0","Aj22","Aj75","AjInclusive"};
  TString AjName_strs[4] = {"Aj0_Aj22_","Aj22_Aj75_",""};

  int lbin, rbin;

  float bin_width_phi, bin_width_eta, diff_max, diff_min, double_diff_max, double_diff_min, no_bgsub_min, no_bgsub_max, signal_min, signal_max;
  double integral, int_err;

  float pt_bin_centers[6] = {0.75,1.5,2.5,3.5,6.,12.};
  float pt_bin_bounds[7] = {0.5,1.,2.,3.,4.,8.,16.};
  float pt_bin_errors[6] = {0.25,.5,.5,.5,2.,4.};
  TLatex *label_pp, *label_per, *label_cent;
 

  float temp1, err1,dx_eta,dx_phi,temp, check_ymin, check_ymax;

  TLine *lineEta, *linePhi;

  TLatex *jet_tex, *aj_tex;
  
  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  TCanvas *dummy = new TCanvas();
  int gstart = 2; 
  int gend = 6;

  float njets[12]= {6689, 8891,3985,1371,5624,2410,729,12313,15382,6395,2100};

 float bc_pbpb, bc_pp, err_change,new_err, old_err;

  for(int l = 1; l<6; l++){

    for(int g=gstart; g<gend; g++){

      for(int p=0; p<3; p++){
 
	if(g%2==0) infile_name = (TString)("../me_correct_mc/HydJet_"+data_mc_type_strs[l]+"_Dijet_Correlations.root");
	else infile_name = (TString)("../me_correct_mc/Pythia_"+data_mc_type_strs[l]+"_Dijet_Correlations.root");
	
	if(l<5)	fin[g][l] = new TFile(infile_name,"READ");

	
	//	if(g%2!=0&&l==4)infile_name = (TString)("../me_correct_mc/HydJet_GenJet_GenTrack_Sube0_Dijet_Correlations.root");

	//----------------------------------------------------
	//  Start of main i & j loops 
	//-----------------------------------------------------

	if(g<4)name_stem = (TString)("Yield_pTweighted_SubLeading_");
	else name_stem = (TString)("Yield_pTweighted_Leading_");

	if(g%2==0)name_stem+="Hydjet_";
	else name_stem+="Pythia_";
   
	for(int i=0; i<nTrkPtBins; i++){
	
	  for (int j=0; j<4; j++){

	    if(g%2!=0&&j<3)continue;

	    //	    if(g%2!=0&&l==4)  name_stem.ReplaceAll("Pythia","Hydjet");

	    if(l<5){
	      raw[g][i][j][p][l] = (TH2D*)fin[g][l]->Get((TString)(name_stem+CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_"+AjName_strs[p]+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]))->Clone((TString)(name_stem+"_"+data_mc_type_strs[l]+"_"+CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_"+AjName_strs[p]+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]));

	      

	    }else{
	      raw[g][i][j][p][l] = (TH2D*)raw[g][i][j][p][1]->Clone((TString)(name_stem+"_"+data_mc_type_strs[l]+"_"+CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_"+AjName_strs[p]+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]));
	      //  if(g%2==0)      raw[g][i][j][p][l]->Add(raw[g+1][i][3][p][4] ,-1.);
	      //      else    
		raw[g][i][j][p][l]->Add(raw[g][i][j][p][4] ,-1.);
	    }



	    TString result_string = "Result_";
	    if(g<4)result_string+="SubLeading_";
	    else result_string+="Leading_";
	    
	    if(g%2==0)result_string+="Hydjet_";
	    else result_string+="Pythia_";
   
	    

	    result[g][i][j][p][l] = (TH2D*)raw[g][i][j][p][l]->Clone((TString)(result_string+data_mc_type_strs[l]+"_"+CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_"+AjName_strs[p]+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]));


	    llimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(1.+0.0001);
	    rlimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(2.-0.0001);
	    
	    background_proj[g][i][j][p][l] = (TH1D*)result[g][i][j][p][l]->ProjectionY(Form("ProjectedBackground%d%d%d%d",g,i,j,p),llimiteta,rlimiteta);



	    llimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(-2.+0.0001);
	    rlimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(-1.-0.0001);
	    
	    background_left[g][i][j][p][l] = (TH1D*)result[g][i][j][p][l]->ProjectionY(Form("LeftSideBackground%d%d%d%d",g,i,j,p),llimiteta,rlimiteta);
	

	    if(l>3&&(i==0||i==2)&&j>1){


	    llimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(1.3+0.0001);
	    rlimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(2.3-0.0001);
	    
	    background_proj[g][i][j][p][l] = (TH1D*)result[g][i][j][p][l]->ProjectionY(Form("ProjectedBackground%d%d%d%d",g,i,j,p),llimiteta,rlimiteta);



	    llimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(-2.3+0.0001);
	    rlimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(-1.3-0.0001);
	    
	    background_left[g][i][j][p][l] = (TH1D*)result[g][i][j][p][l]->ProjectionY(Form("LeftSideBackground%d%d%d%d",g,i,j,p),llimiteta,rlimiteta);



	    }
	    

	    if(l>3&&(i==1)&&j>1){


	      llimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(1.4+0.0001);
	      rlimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(2.4-0.0001);
	    
	      background_proj[g][i][j][p][l] = (TH1D*)result[g][i][j][p][l]->ProjectionY(Form("ProjectedBackground%d%d%d%d",g,i,j,p),llimiteta,rlimiteta);



	      llimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(-2.4+0.0001);
	      rlimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(-1.4-0.0001);
	    
	      background_left[g][i][j][p][l] = (TH1D*)result[g][i][j][p][l]->ProjectionY(Form("LeftSideBackground%d%d%d%d",g,i,j,p),llimiteta,rlimiteta);



	    }


	    if(i>3&&l==4&&j<2&&g==4){

	
	      llimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(0.8+0.0001);
	      rlimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(1.8-0.0001);
	    
	      background_proj[g][i][j][p][l] = (TH1D*)result[g][i][j][p][l]->ProjectionY(Form("ProjectedBackground%d%d%d%d",g,i,j,p),llimiteta,rlimiteta);
	      background_left[g][i][j][p][l] = (TH1D*)result[g][i][j][p][l]->ProjectionY(Form("LeftSideBackground%d%d%d%d",g,i,j,p),llimiteta,rlimiteta);
	 
	    }


	    background_proj[g][i][j][p][l]->Add(background_left[g][i][j][p][l]);

	    dx_phi =  background_proj[g][i][j][p][l]->GetBinWidth(1);

	 
	    background_proj_rebin[g][i][j][p][l] = Rebin_dPhi(background_proj[g][i][j][p][l]);

	    background_proj_rebin[g][i][j][p][l]->Scale(5./2./dx_phi);

	    background_proj[g][i][j][p][l]->Scale(1./2/(rlimiteta-llimiteta+1));

	    for(int k = 1 ; k<background_proj_rebin[g][i][j][p][l]->GetNbinsX()/2+1; k++){
	    
	      temp = (background_proj_rebin[g][i][j][p][l]->GetBinContent(background_proj_rebin[g][i][j][p][l]->GetNbinsX()+1-k)+background_proj_rebin[g][i][j][p][l]->GetBinContent(k))/2.;
	      
	 
	      background_proj_rebin[g][i][j][p][l]->SetBinContent(background_proj_rebin[g][i][j][p][l]->GetNbinsX()+1-k, temp);
	      background_proj_rebin[g][i][j][p][l]->SetBinContent(k, temp);
	    }

	    background[g][i][j][p][l] = (TH2D*)result[g][i][j][p][l]->Clone(Form("Background%d%d%d%d",g,i,j,l));

	    for(int k = 1;  k<result[g][i][j][p][l]->GetNbinsY(); k++){
	      temp1 = background_proj[g][i][j][p][l]->GetBinContent(k);
	      err1 = background_proj[g][i][j][p][l]->GetBinError(k);
	    
	      for(int m = 1;  m<result[g][i][j][p][l]->GetNbinsX(); m++){
		background[g][i][j][p][l]->SetBinContent(m,k,temp1);
		background[g][i][j][p][l]->SetBinError(m,k,err1);
		
	      }

	    }

	    result[g][i][j][p][l]->Add(background[g][i][j][p][l],-1.);
	  }
	}
      }
    }
  
 
    for(int g=gstart; g<gend; g++){

      for(int p=0; p<3; p++){
 

	for(int i = 0; i<nTrkPtBins; i++){
	  for(int j = 0; j<nCBins; j++){
	    if(g%2!=0&&j<3)continue;

	    if(l==5){

	      TString result_string = "Result_";
	      if(g<4)result_string+="SubLeading_";
	      else result_string+="Leading_";
	    
	      if(g%2==0)result_string+="Hydjet_";
	      else result_string+="Pythia_";
   
	      result[g][i][j][p][l] = (TH2D*)result[g][i][j][p][1]->Clone((TString)(result_string+data_mc_type_strs[l]+"_"+CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_"+AjName_strs[p]+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]));

	      if(g%2==0) result[g][i][j][p][l]->Add( result[g+1][i][3][p][4],-1.);
	      else 	result[g][i][j][p][l]->Add( result[g][i][j][p][4],-1.);


	      background_proj_rebin[g][i][j][p][l] = (TH1D*)background_proj_rebin[g][i][j][p][1]->Clone((TString)("BkgProjRebin"+data_mc_type_strs[l]+"_"+CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_"+AjName_strs[p]+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]));


	      if(g%2==0)    background_proj_rebin[g][i][j][p][l]->Add( background_proj_rebin[g+1][i][3][p][4],-1.);
	      else    background_proj_rebin[g][i][j][p][l]->Add( background_proj_rebin[g][i][j][p][4],-1.);
	    }

	    //-------------------------------
	    //dEta projection
	    //------------------------

	    TString eta_proj_name=  make_name("Eta_Projections_",g+6,i,j,p,pTlabel,centlabel,Ajlabel);
	    
	    llimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(-2.5+.001);
	    rlimiteta = result[g][i][j][p][l]->GetXaxis()->FindBin(2.5-.001);

	    llimitphi = result[g][i][j][p][l]->GetYaxis()->FindBin(-philim+.001);
	    rlimitphi = result[g][i][j][p][l]->GetYaxis()->FindBin(philim-.001);
	    

	    eta_proj[g][i][j][p][l] = result[g][i][j][p][l]->ProjectionX(eta_proj_name,llimitphi,rlimitphi);
	    dx_eta = eta_proj[g][i][j][p][l]->GetBinWidth(1);
	    eta_proj[g][i][j][p][l]->Scale(1/dx_eta);

	    TString eta_proj_name_rebin = eta_proj_name;
	    eta_proj_name_rebin.ReplaceAll("Eta_Proj","Eta_Proj_Rebin");
	 
	    eta_proj_rebin[g][i][j][p][l] = (TH1D*)Rebin_dEta(eta_proj[g][i][j][p][l]);
	    eta_proj_rebin[g][i][j][p][l]->SetName(eta_proj_name_rebin);
	 

	    //-------------------------------
	    //dPhi projection
	    //------------------------

	    TString phi_proj_name= make_name("Phi_Projections_",g+6,i,j,p,pTlabel,centlabel,Ajlabel);
	 
	    phi_proj[g][i][j][p][l] = result[g][i][j][p][l]->ProjectionY(phi_proj_name,llimiteta,rlimiteta);
	    dx_phi = phi_proj[g][i][j][p][l]->GetBinWidth(1);
	    phi_proj[g][i][j][p][l]->Scale(1/dx_phi);

	    TString phi_proj_name_rebin = phi_proj_name;
	    phi_proj_name_rebin.ReplaceAll("Phi_Proj","Phi_Proj_Rebin");
	 
	    signal_dPhi_rebin[g][i][j][p][l] = (TH1D*)Rebin_dPhi(phi_proj[g][i][j][p][l]);
	 
	    for(int k = 1 ; k<eta_proj_rebin[g][i][j][p][l]->GetNbinsX()/2+1; k++){
	    
	      temp = (eta_proj_rebin[g][i][j][p][l]->GetBinContent(eta_proj_rebin[g][i][j][p][l]->GetNbinsX()+1-k)+eta_proj_rebin[g][i][j][p][l]->GetBinContent(k))/2.;
	      eta_proj_rebin[g][i][j][p][l]->SetBinContent(eta_proj_rebin[g][i][j][p][l]->GetNbinsX()+1-k, temp);
	      eta_proj_rebin[g][i][j][p][l]->SetBinContent(k, temp);
	    }

	    for(int k = 1 ; k<signal_dPhi_rebin[g][i][j][p][l]->GetNbinsX()/2+1; k++){
	    
	      temp = (signal_dPhi_rebin[g][i][j][p][l]->GetBinContent(signal_dPhi_rebin[g][i][j][p][l]->GetNbinsX()+1-k)+signal_dPhi_rebin[g][i][j][p][l]->GetBinContent(k))/2.;
	      signal_dPhi_rebin[g][i][j][p][l]->SetBinContent(signal_dPhi_rebin[g][i][j][p][l]->GetNbinsX()+1-k, temp);
	      signal_dPhi_rebin[g][i][j][p][l]->SetBinContent(k, temp);
	    }

	  }
	}

      }

    }
   
    
  }

  cout<<"now for study yield step..."<<endl;

 
  for(int g=gstart; g<gend; g++){
    
    for(int p = 0; p<3; p++){
      
      
      if(g<4)datalabel = "SubLeading";
      else datalabel = "Leading";
   
    
      for(int i = 0; i<nTrkPtBins; i++){
	for(int j = 0; j<nCBins; j++){
	  if(g%2!=0&&j<3)continue;
	  //	  if(j==1)continue;

	  for(int l = 1; l<6; l++){

	    //	    cout<<g<<" "<<i<<" "<<j<<" "<<p<<" "<<l<<endl;
	  
	    switch(i){
	    case 0: 
	      check_ymax = 11.;
	      check_ymin = -5.; 
	      break;
	    case 1: 
	      check_ymax = 11.; 
	      check_ymin = -5.; 
	      break;
	    case 2: 
	      check_ymax = 20.; 
	      check_ymin = -5.; 
	      break;
	    case 3: 
	      check_ymax = 25.; 
	      check_ymin = -5.; 
	   
	      break;
	    case 4: 
	      check_ymax = 100.;
	      check_ymin = -10.; 
	      break;
	    case 5: 
	      check_ymax = 400.;
	      check_ymin = -10.; 
	      break;
	    }


	    eta_proj_rebin[g][i][j][p][l]->SetMarkerSize(1);
	    signal_dPhi_rebin[g][i][j][p][l]->SetMarkerSize(1);


	    eta_proj_rebin[g][i][j][p][l]->SetMarkerColor(kBlack);
	    eta_proj_rebin[g][i][j][p][l]->SetLineColor(kBlack);
	  
	    signal_dPhi_rebin[g][i][j][p][l]->SetMarkerColor(kBlack);
	    signal_dPhi_rebin[g][i][j][p][l]->SetLineColor(kBlack);
	
	    //------------------
	    // dEta plotting
	    //-------------------

	    eta_proj_rebin[g][i][j][p][l]->SetMaximum(check_ymax);
	    eta_proj_rebin[g][i][j][p][l]->SetMinimum(check_ymin);

   
	    eta_proj_rebin[g][i][j][p][l]->SetMarkerColor(1);
	    eta_proj_rebin[g][i][j][p][l]->SetAxisRange(-2.5+.2,2.5-.2,"x");


	  
	    eta_proj_rebin[g][i][j][p][l]->GetYaxis()->SetLabelSize(tstitle);
	    eta_proj_rebin[g][i][j][p][l]->GetXaxis()->SetLabelSize(tstitle);
	    eta_proj_rebin[g][i][j][p][l]->GetXaxis()->SetTitle("#Delta#eta");
	    eta_proj_rebin[g][i][j][p][l]->GetXaxis()->SetTitleSize(tstitle);
	    eta_proj_rebin[g][i][j][p][l]->GetXaxis()->SetTitleOffset(xoffset);
	    eta_proj_rebin[g][i][j][p][l]->GetYaxis()->SetTitle("Y#equiv1/N_{jet} d^{2}N/(d#Delta#eta) (GeV/c)^{-1}");
	    eta_proj_rebin[g][i][j][p][l]->GetYaxis()->SetTitleOffset(yoffset);
	    eta_proj_rebin[g][i][j][p][l]->GetYaxis()->SetTitleSize(tstitle);

 
	    eta_proj_rebin[g][i][j][p][l]->GetXaxis()->CenterTitle();
	    eta_proj_rebin[g][i][j][p][l]->GetYaxis()->CenterTitle();
	    eta_proj_rebin[g][i][j][p][l]->GetXaxis()->SetLabelSize(0.);
	  
	    if(g%2==0){
	      eta_proj_rebin[g][i][j][p][l]->GetYaxis()->SetTitleSize(0.0);
	      eta_proj_rebin[g][i][j][p][l]->GetYaxis()->SetLabelSize(0.0);
	    }
	
	    signal_dPhi_rebin[g][i][j][p][l]->SetMaximum(check_ymax);
	    signal_dPhi_rebin[g][i][j][p][l]->SetMinimum(check_ymin);
		    
	    signal_dPhi_rebin[g][i][j][p][l]->GetXaxis()->CenterTitle();
	    signal_dPhi_rebin[g][i][j][p][l]->GetYaxis()->CenterTitle();
	    signal_dPhi_rebin[g][i][j][p][l]->GetYaxis()->SetLabelSize(tstitle);
	    signal_dPhi_rebin[g][i][j][p][l]->GetXaxis()->SetLabelSize(tstitle);
	    signal_dPhi_rebin[g][i][j][p][l]->GetXaxis()->SetTitle("#Delta#phi");
	    signal_dPhi_rebin[g][i][j][p][l]->GetXaxis()->SetTitleSize(tstitle);
	    signal_dPhi_rebin[g][i][j][p][l]->GetXaxis()->SetTitleOffset(xoffset);
	    signal_dPhi_rebin[g][i][j][p][l]->GetYaxis()->SetTitle("Y#equiv1/N_{jet} d^{2}N/(d#Delta#phi) (GeV/c)^{-1}");
	    signal_dPhi_rebin[g][i][j][p][l]->GetYaxis()->SetTitleOffset(yoffset);
	    signal_dPhi_rebin[g][i][j][p][l]->GetYaxis()->SetTitleSize(tstitle);
	    signal_dPhi_rebin[g][i][j][p][l]->GetXaxis()->SetLabelSize(0.);
	  
	    if(g%2==0){
	
	      signal_dPhi_rebin[g][i][j][p][l]->GetYaxis()->SetTitleSize(0.0);
	      signal_dPhi_rebin[g][i][j][p][l]->GetYaxis()->SetLabelSize(0.0);
	    }
	  } //l
	  //------------------
	  //Apply corrections
	  //-------------------

	
	  if(g%2==0){
	    TString spill_name = make_name("SpillOvers_Eta_pTweighted_",g+6,i,j,p,pTlabel,centlabel,Ajlabel);
	    spill_name.ReplaceAll("Pt100_Pt300_","");
	   
	    spill_over_dEta[g][i][j][p] = (TH1D*)f_spillover->Get(spill_name)->Clone(spill_name);
	    spill_name.ReplaceAll("Eta","Phi");
	    spill_over_dPhi[g][i][j][p] = (TH1D*)f_spillover->Get(spill_name)->Clone(spill_name);

	    eta_proj_rebin[g][i][j][p][1]->Add(spill_over_dEta[g][i][j][p],-1.);
	    signal_dPhi_rebin[g][i][j][p][1]->Add(spill_over_dPhi[g][i][j][p],-1.);


	    eta_proj_rebin[g][i][j][p][5]->Add(spill_over_dEta[g][i][j][p],-1.);
	    signal_dPhi_rebin[g][i][j][p][5]->Add(spill_over_dPhi[g][i][j][p],-1.);


	    eta_proj_rebin[g][i][j][p][2]->Add(spill_over_dEta[g][i][j][p],-1.);
	    signal_dPhi_rebin[g][i][j][p][2]->Add(spill_over_dPhi[g][i][j][p],-1.);


	  }	 
	
	
	  TString  jff_name = make_name("JFF_Residual_Eta_",g+6,i,j,p,pTlabel,centlabel,Ajlabel);

	  if(i>3) jff_residual_dEta[g][i][j][p] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);
	  else	  jff_residual_dEta[g][i][j][p] = (TH1D*)f_jff_gen->Get(jff_name)->Clone(jff_name);
	  
	  jff_name.ReplaceAll("Eta","Phi");
	  if(i>3) jff_residual_dPhi[g][i][j][p] = (TH1D*)f_jff_reco->Get(jff_name)->Clone(jff_name);
	  else  jff_residual_dPhi[g][i][j][p] = (TH1D*)f_jff_gen->Get(jff_name)->Clone(jff_name);

	  eta_proj_rebin[g][i][j][p][1]->Add(jff_residual_dEta[g][i][j][p],-1.);

	  signal_dPhi_rebin[g][i][j][p][1]->Add(jff_residual_dPhi[g][i][j][p],-1.);
	

	  eta_proj_rebin[g][i][j][p][5]->Add(jff_residual_dEta[g][i][j][p],-1.);

	  signal_dPhi_rebin[g][i][j][p][5]->Add(jff_residual_dPhi[g][i][j][p],-1.);

	  eta_proj_rebin[g][i][j][p][2]->Add(jff_residual_dEta[g][i][j][p],-1.);

	  signal_dPhi_rebin[g][i][j][p][2]->Add(jff_residual_dPhi[g][i][j][p],-1.);


	} //close j
      } // close i 
    
    } //close l
  
  } //close g
 

  cout<<"*****Switching to plotting code*********"<<endl;;

  ///------------------------------------------------
  ///------------------------------------------------
  // And now run through the full data plotting code
  ///------------------------------------------------
  ///------------------------------------------------



  for(int mc_code = 1; mc_code < 6; mc_code++){

    TString outfile_name;
	
    outfile_name = (TString)("HydJet_"+data_mc_type_strs[mc_code]+"_Projections.root");
    fout[0][mc_code] = new TFile(outfile_name,"RECREATE");
    outfile_name = (TString)("Pythia_"+data_mc_type_strs[mc_code]+"_Dijet_Projections.root");
    fout[1][mc_code] = new TFile(outfile_name,"RECREATE");
   

    for(int l = 0; l<3; l++){

      for(int i = 0; i<6; i++){

	for(int j = 0; j<4; j++){


	  if(i==0&&mc_code==5&&j==0){

	    TString syst_name_pbpb = make_name("Syst_",0,4,0,l,pTlabel, centlabel, Ajlabel);  
	    background_syst_rebin[0][0][0][l] = (TH1D*)f_in_syst->Get(syst_name_pbpb)->Clone(syst_name_pbpb);

	    cout<<"here"<<endl;
	    TString syst_name_pp = make_name("Syst_",1,4,3,l,pTlabel, centlabel, Ajlabel);  
	    background_syst_rebin[1][0][3][l] = (TH1D*)f_in_syst->Get(syst_name_pp)->Clone(syst_name_pp);

	    syst_name_pbpb.ReplaceAll("Syst","Syst_Diff");
	    background_syst_rebin[7][0][0][l] = (TH1D*)background_syst_rebin[0][0][0][l]->Clone(syst_name_pbpb);
	    background_syst_rebin[7][0][0][l]->Add(background_syst_rebin[1][0][3][l],-1.);


	    syst_name_pbpb = make_name("Syst_",0,4,2,l,pTlabel, centlabel, Ajlabel);  
	    background_syst_rebin[0][0][3][l] = (TH1D*)f_in_syst->Get(syst_name_pbpb)->Clone(syst_name_pbpb);

	    syst_name_pbpb.ReplaceAll("Syst","Syst_Diff");
	    background_syst_rebin[7][0][3][l] = (TH1D*)background_syst_rebin[0][0][0][l]->Clone(syst_name_pbpb);
	    background_syst_rebin[7][0][3][l]->Add(background_syst_rebin[1][0][3][l],-1.);



	
	    cout<<"background done"<<endl;

	    signal_dPhi_syst[2][0][0][l] = (TH1D*)f_in_syst->Get((TString)("SubLeading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"))->Clone((TString)("SubLeading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"));



	    signal_dPhi_syst[4][0][0][l] = (TH1D*)f_in_syst->Get((TString)("Leading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"))->Clone((TString)("Leading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"));

	    signal_dPhi_syst[2][0][3][l] = (TH1D*)f_in_syst->Get((TString)("SubLeading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent10_Cent30"))->Clone((TString)("SubLeading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent0_Cent10"));

	    signal_dPhi_syst[4][0][3][l] = (TH1D*)f_in_syst->Get((TString)("Leading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent10_Cent30"))->Clone((TString)("Leading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent0_Cent10"));


	    signal_dPhi_syst[3][0][3][l] = (TH1D*)f_in_syst->Get((TString)("SubLeading_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"))->Clone((TString)("SubLeading_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"));

	    signal_dPhi_syst[5][0][3][l] = (TH1D*)f_in_syst->Get((TString)("Leading_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"))->Clone((TString)("Leading_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"));
   
	    cout<<"here"<<endl;

	    sub_lead_dPhi_syst[2][0][0][l] = (TH1D*)f_in_syst->Get((TString)("Diff_NoBkgSub_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"))->Clone((TString)("Diff_NoBkgSub_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"));
  
	    sub_lead_dPhi_syst[2][0][3][l] = (TH1D*)f_in_syst->Get((TString)("Diff_NoBkgSub_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent10_Cent30"))->Clone((TString)("Diff_NoBkgSub_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent0_Cent10"));

	    sub_lead_dPhi_syst[3][0][3][l] = (TH1D*)f_in_syst->Get((TString)("Diff_NoBkgSub_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"))->Clone((TString)("Diff_NoBkgSub_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"));

	    cout<<"and here"<<endl;

    
	    in_name = make_name("PbPb_pp_WithErrors_Syst_",8,0,0,l,pTlabel,centlabel,Ajlabel);
	    signal_dPhi_syst[8][0][0][l] = (TH1D*)signal_dPhi_syst[2][0][0][l]->Clone(in_name);
	    signal_dPhi_syst[8][0][0][l]->Add( signal_dPhi_syst[3][0][3][l],-1.);

	    in_name.ReplaceAll("SubLeading","Leading");
	    signal_dPhi_syst[10][0][0][l] = (TH1D*)signal_dPhi_syst[4][0][0][l]->Clone(in_name);
	    signal_dPhi_syst[10][0][0][l]->Add( signal_dPhi_syst[5][0][3][l],-1.);


    
   
	    for(int k = 1; k< signal_dPhi_syst[8][0][0][l]->GetNbinsX()+1; k++){
      
	      bc_pbpb =  signal_dPhi_syst[2][0][0][l]->GetBinContent(k);
	      bc_pp =  signal_dPhi_syst[3][0][3][l]->GetBinContent(k);
	      err_change = TMath::Sqrt(0.04*0.04+0.05*0.05)*TMath::Abs(TMath::Sqrt(bc_pbpb*bc_pbpb+bc_pp*bc_pp)-TMath::Abs(bc_pbpb - bc_pp));
	      old_err = signal_dPhi_syst[8][0][0][l]->GetBinError(k);
	      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
     
	      signal_dPhi_syst[8][0][0][l]->SetBinError(k,new_err);
    
	      bc_pbpb =  signal_dPhi_syst[4][0][0][l]->GetBinContent(k);
	      bc_pp =  signal_dPhi_syst[5][0][3][l]->GetBinContent(k);
	      err_change = TMath::Sqrt(0.04*0.04+0.05*0.05)*TMath::Abs(TMath::Sqrt(bc_pbpb*bc_pbpb+bc_pp*bc_pp)-TMath::Abs(bc_pbpb - bc_pp));
	      old_err = signal_dPhi_syst[10][0][0][l]->GetBinError(k);
	      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
	      signal_dPhi_syst[10][0][0][l]->SetBinError(k,new_err);

    

	    }




	    in_name.ReplaceAll("SubLeading","DoubleDiff");
	    signal_dPhi_syst[6][0][0][l] = (TH1D*)signal_dPhi_syst[8][0][0][l]->Clone(in_name);
	    signal_dPhi_syst[6][0][0][l]->Add( signal_dPhi_syst[10][0][0][l]);



   
	    for(int k = 1; k< signal_dPhi_syst[8][0][0][l]->GetNbinsX()+1; k++){
      
	      bc_pbpb =  signal_dPhi_syst[2][0][0][l]->GetBinContent(k);
	      bc_pp =  signal_dPhi_syst[4][0][0][l]->GetBinContent(k);
	      err_change = TMath::Sqrt(0.04*.04+.05*0.05)*(bc_pbpb-TMath::Abs(signal_dPhi_syst[6][0][0][l]->GetBinContent(k)));

 
	      old_err = signal_dPhi_syst[6][0][0][l]->GetBinError(k);
	      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
	      signal_dPhi_syst[6][0][0][l]->SetBinError(k,new_err);

	      cout<<k<<" "<<bc_pbpb<<" "<<bc_pp<<" "<<err_change<<" "<<old_err<<" "<<new_err<<endl;
    
	    }




	    in_name = make_name("PbPb_pp_WithErrors_Syst_",8,0,0,l,pTlabel,centlabel,Ajlabel);
	    signal_dPhi_syst[8][0][3][l] = (TH1D*)signal_dPhi_syst[2][0][3][l]->Clone(in_name);
	    signal_dPhi_syst[8][0][3][l]->Add( signal_dPhi_syst[3][0][3][l],-1.);

	    in_name.ReplaceAll("SubLeading","Leading");
	    signal_dPhi_syst[10][0][3][l] = (TH1D*)signal_dPhi_syst[4][0][3][l]->Clone(in_name);
	    signal_dPhi_syst[10][0][3][l]->Add( signal_dPhi_syst[5][0][3][l],-1.);


    
   
	    for(int k = 1; k< signal_dPhi_syst[8][0][3][l]->GetNbinsX()+1; k++){
      
	      bc_pbpb =  signal_dPhi_syst[2][0][3][l]->GetBinContent(k);
	      bc_pp =  signal_dPhi_syst[3][0][3][l]->GetBinContent(k);
	      err_change = TMath::Sqrt(0.04*0.04+0.05*0.05)*TMath::Abs(TMath::Sqrt(bc_pbpb*bc_pbpb+bc_pp*bc_pp)-TMath::Abs(bc_pbpb - bc_pp));
	      old_err = signal_dPhi_syst[8][0][3][l]->GetBinError(k);
	      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
     
	      signal_dPhi_syst[8][0][3][l]->SetBinError(k,new_err);
    
	      bc_pbpb =  signal_dPhi_syst[4][0][3][l]->GetBinContent(k);
	      bc_pp =  signal_dPhi_syst[5][0][3][l]->GetBinContent(k);
	      err_change = TMath::Sqrt(0.04*0.04+0.05*0.05)*TMath::Abs(TMath::Sqrt(bc_pbpb*bc_pbpb+bc_pp*bc_pp)-TMath::Abs(bc_pbpb - bc_pp));
	      old_err = signal_dPhi_syst[10][0][3][l]->GetBinError(k);
	      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
	      signal_dPhi_syst[10][0][3][l]->SetBinError(k,new_err);

    

	    }




	    in_name.ReplaceAll("SubLeading","DoubleDiff");
	    signal_dPhi_syst[6][0][3][l] = (TH1D*)signal_dPhi_syst[8][0][3][l]->Clone(in_name);
	    signal_dPhi_syst[6][0][3][l]->Add( signal_dPhi_syst[10][0][3][l]);



   
	    for(int k = 1; k< signal_dPhi_syst[8][0][3][l]->GetNbinsX()+1; k++){
      
	      bc_pbpb =  signal_dPhi_syst[2][0][3][l]->GetBinContent(k);
	      bc_pp =  signal_dPhi_syst[4][0][3][l]->GetBinContent(k);
	      err_change = TMath::Sqrt(0.04*.04+.05*0.05)*(bc_pbpb-TMath::Abs(signal_dPhi_syst[6][0][3][l]->GetBinContent(k)));

 
	      old_err = signal_dPhi_syst[6][0][3][l]->GetBinError(k);
	      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
	      signal_dPhi_syst[6][0][3][l]->SetBinError(k,new_err);

	      cout<<k<<" "<<bc_pbpb<<" "<<bc_pp<<" "<<err_change<<" "<<old_err<<" "<<new_err<<endl;
    
	    }





	   
	    cout<<"ready"<<endl;

	    for(int k = 0; k< sub_lead_dPhi_syst[3][0][3][l]->GetNbinsX()+1; k++){

	      background_syst_rebin[0][0][0][l]->SetBinContent(k,0.000000001);
	      background_syst_rebin[7][0][0][l]->SetBinContent(k,0.000000001);
	      background_syst_rebin[0][0][3][l]->SetBinContent(k,0.000000001);
	      background_syst_rebin[1][0][3][l]->SetBinContent(k,0.000000001);
	      background_syst_rebin[7][0][3][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[2][0][0][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[4][0][0][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[2][0][3][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[4][0][3][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[6][0][0][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[8][0][0][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[10][0][0][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[6][0][3][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[8][0][3][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[10][0][3][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[3][0][3][l]->SetBinContent(k,0.000000001);
	      signal_dPhi_syst[5][0][3][l]->SetBinContent(k,0.000000001);
	      sub_lead_dPhi_syst[2][0][0][l]->SetBinContent(k,0.000000001);
	      sub_lead_dPhi_syst[2][0][3][l]->SetBinContent(k,0.000000001);
	      sub_lead_dPhi_syst[3][0][3][l]->SetBinContent(k,0.000000001);

	    }

	    cout<<"set zeros"<<endl;


	    background_syst_rebin[0][0][0][l]->SetMarkerSize(0.05);
	    background_syst_rebin[7][0][0][l]->SetMarkerSize(0.05);
	    background_syst_rebin[0][0][3][l]->SetMarkerSize(0.05);
	    background_syst_rebin[1][0][3][l]->SetMarkerSize(0.05);
	    background_syst_rebin[7][0][3][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[2][0][0][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[4][0][0][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[2][0][3][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[6][0][0][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[8][0][0][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[10][0][0][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[6][0][3][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[8][0][3][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[10][0][3][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[3][0][3][l]->SetMarkerSize(0.05);
	    signal_dPhi_syst[5][0][3][l]->SetMarkerSize(0.05);
	    sub_lead_dPhi_syst[2][0][0][l]->SetMarkerSize(0.05);
	    sub_lead_dPhi_syst[2][0][3][l]->SetMarkerSize(0.05);
	    sub_lead_dPhi_syst[3][0][3][l]->SetMarkerSize(0.05);


	    background_syst_rebin[0][0][0][l]->SetFillStyle(3004);
	    background_syst_rebin[7][0][0][l]->SetFillStyle(3004);
	    background_syst_rebin[0][0][3][l]->SetFillStyle(3004);
	    background_syst_rebin[1][0][3][l]->SetFillStyle(3004);
	    background_syst_rebin[7][0][3][l]->SetFillStyle(3004);
	    signal_dPhi_syst[2][0][0][l]->SetFillStyle(3004);
	    signal_dPhi_syst[4][0][0][l]->SetFillStyle(3004);
	    signal_dPhi_syst[2][0][3][l]->SetFillStyle(3004);
	    signal_dPhi_syst[6][0][0][l]->SetFillStyle(3004);
	    signal_dPhi_syst[8][0][0][l]->SetFillStyle(3004);
	    signal_dPhi_syst[10][0][0][l]->SetFillStyle(3004);
	    signal_dPhi_syst[6][0][3][l]->SetFillStyle(3004);
	    signal_dPhi_syst[8][0][3][l]->SetFillStyle(3004);
	    signal_dPhi_syst[10][0][3][l]->SetFillStyle(3004);
	    signal_dPhi_syst[3][0][3][l]->SetFillStyle(3004);
	    signal_dPhi_syst[5][0][3][l]->SetFillStyle(3004);
	    sub_lead_dPhi_syst[2][0][0][l]->SetFillStyle(3004);
	    sub_lead_dPhi_syst[2][0][3][l]->SetFillStyle(3004);
	    sub_lead_dPhi_syst[3][0][3][l]->SetFillStyle(3004);



	    background_syst_rebin[0][0][0][l]->SetFillColor(kBlack);
	    background_syst_rebin[7][0][0][l]->SetFillColor(kBlack);
	    background_syst_rebin[0][0][3][l]->SetFillColor(kBlack);
	    background_syst_rebin[1][0][3][l]->SetFillColor(kBlack);
	    background_syst_rebin[7][0][3][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[2][0][0][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[4][0][0][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[2][0][3][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[6][0][0][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[8][0][0][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[10][0][0][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[6][0][3][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[8][0][3][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[10][0][3][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[3][0][3][l]->SetFillColor(kBlack);
	    signal_dPhi_syst[5][0][3][l]->SetFillColor(kBlack);
	    sub_lead_dPhi_syst[2][0][0][l]->SetFillColor(kBlack);
	    sub_lead_dPhi_syst[2][0][3][l]->SetFillColor(kBlack);
	    sub_lead_dPhi_syst[3][0][3][l]->SetFillColor(kBlack);






	    in_name = make_name("Double_Diff_NoBkgSub_Syst_",6,0,0,l,pTlabel,centlabel,Ajlabel);
	    sub_lead_dPhi_syst[6][0][0][l] = (TH1D*)sub_lead_dPhi_syst[2][0][0][l]->Clone(in_name);
	    sub_lead_dPhi_syst[6][0][0][l]->Add( sub_lead_dPhi_syst[3][0][3][l],-1.);


	    in_name = make_name("Double_Diff_NoBkgSub_Syst_",6,0,2,l,pTlabel,centlabel,Ajlabel);
	    sub_lead_dPhi_syst[6][0][3][l] = (TH1D*) sub_lead_dPhi_syst[2][0][3][l]->Clone(in_name);
	    sub_lead_dPhi_syst[6][0][3][l]->Add( sub_lead_dPhi_syst[3][0][3][l],-1.);

	  }
	  cout<<"and here"<<endl;  
	    for(int g = 0; g<2; g++){

	      if(j!=3&&g==1)continue;

	    in_name = make_name("sub_lead_dPhi_NoErrUp_",g,i,j,l,pTlabel,centlabel,Ajlabel);
	    sub_lead_dPhi_noerr_up[g][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	    in_name.ReplaceAll("Up","Down");
	    sub_lead_dPhi_noerr_down[g][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);


	    in_name = make_name("Sub_Lead_2D",g,i,j,l,pTlabel,centlabel,Ajlabel);

	    result_sub_lead[g][i][j][l][mc_code] = (TH2D*) raw[g+2][i][j][l][mc_code]->Clone(in_name);
	    result_sub_lead[g][i][j][l][mc_code]->Add(raw[g+4][i][j][l][mc_code],-1.);

	    in_name.ReplaceAll("2D","dPhi");

	    lbin = result_sub_lead[g][i][j][l][mc_code]->GetXaxis()->FindBin(-2.499);
	    rbin = result_sub_lead[g][i][j][l][mc_code]->GetXaxis()->FindBin(2.499);

	    result_sub_lead_proj[g][i][j][l][mc_code]= (TH1D*)result_sub_lead[g][i][j][l][mc_code]->ProjectionY(in_name,lbin,rbin);

	    dx_phi =  result_sub_lead_proj[g][i][j][l][mc_code]->GetBinWidth(1);
	    result_sub_lead_proj[g][i][j][l][mc_code]->Scale(1./dx_phi);

	    sub_lead_dPhi_rebin[g][i][j][l][mc_code] = Rebin_dPhi(result_sub_lead_proj[g][i][j][l][mc_code]);


	    cout<<"about to correct "<<g<<" "<<i<<" "<<j<<" "<<l<<endl;
	    //------------------
	    //Spill-over and JFF correct


	    if(g%2==0&&(mc_code==1||mc_code==2||mc_code==5)){
	      sub_lead_dPhi_rebin[g][i][j][l][mc_code]->Add(spill_over_dPhi[g+2][i][j][l],-1.);
	      sub_lead_dPhi_rebin[g][i][j][l][mc_code]->Add(spill_over_dPhi[g+4][i][j][l],1.);
  
	      sub_lead_dPhi_rebin[g][i][j][l][mc_code]->Add(jff_residual_dPhi[g+2][i][j][l],-1.);
	      sub_lead_dPhi_rebin[g][i][j][l][mc_code]->Add(jff_residual_dPhi[g+4][i][j][l]);
	 
	    }else if(mc_code==1||mc_code==2||mc_code==5){
	      sub_lead_dPhi_rebin[g][i][j][l][mc_code]->Add(jff_residual_dPhi[g+2][i][j][l],-1.);
	      sub_lead_dPhi_rebin[g][i][j][l][mc_code]->Add(jff_residual_dPhi[g+4][i][j][l]);
	    }

	
	    cout<<"corrected"<<endl;

	    //--------------------
		  
	    signal_dPhi_rebin[g+4][i][j][l][mc_code]->Scale(-1.);
	 
	    for(int k = 1 ; k<sub_lead_dPhi_rebin[g][i][j][l][mc_code]->GetNbinsX()/2+1; k++){
	    
	      temp = (sub_lead_dPhi_rebin[g][i][j][l][mc_code]->GetBinContent(sub_lead_dPhi_rebin[g][i][j][l][mc_code]->GetNbinsX()+1-k)+sub_lead_dPhi_rebin[g][i][j][l][mc_code]->GetBinContent(k))/2.;
	      sub_lead_dPhi_rebin[g][i][j][l][mc_code]->SetBinContent(sub_lead_dPhi_rebin[g][i][j][l][mc_code]->GetNbinsX()+1-k, temp);
	      sub_lead_dPhi_rebin[g][i][j][l][mc_code]->SetBinContent(k, temp);
	    }


	    in_name = make_name("Diff_",g,i,j,l,pTlabel,centlabel,Ajlabel);
	    in_name+="_Rebin";

	    background_diff_rebin[g][i][j][l][mc_code] = (TH1D*)background_proj_rebin[g+2][i][j][l][mc_code]->Clone(in_name);
	    background_diff_rebin[g][i][j][l][mc_code]->Add(background_proj_rebin[g+4][i][j][l][mc_code],-1.);


	    in_name.ReplaceAll("Diff","Diff_NoErrUp");
	    background_diff_noerr_up[g][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	    in_name.ReplaceAll("Up","Down");
	    background_diff_noerr_down[g][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);

	
	    in_name = make_name("Sub_Lead_NoErrUp_",g,i,j,l,pTlabel,centlabel,Ajlabel);
	    sub_lead_dPhi_noerr_up[g][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	    in_name.ReplaceAll("Up","Down");
	    sub_lead_dPhi_noerr_down[g][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);



	    in_name = make_name("dPhi_NoErrUp_",g+2,i,j,l,pTlabel,centlabel,Ajlabel);
	    signal_dPhi_noerr_up[g+2][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	    in_name.ReplaceAll("Up","Down");
	    signal_dPhi_noerr_down[g+2][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);

	    in_name.ReplaceAll("SubLeading","Leading");
	    signal_dPhi_noerr_down[g+4][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	    in_name.ReplaceAll("Down","Up");
	    signal_dPhi_noerr_up[g+4][i][j][l][mc_code] = new TH1D(in_name,"",nbounds_phi-1 ,bin_bounds_phi);

  
	    fout[g][mc_code]->cd();

	    background_diff_rebin[g][i][j][l][mc_code]->Write();	    
	    sub_lead_dPhi_rebin[g][i][j][l][mc_code]->Write();	    
	    signal_dPhi_rebin[g+2][i][j][l][mc_code]->Write();	    
	    signal_dPhi_rebin[g+4][i][j][l][mc_code]->Write();	    

	  }

	}
    
	for(int j = 0; j<4; j++){

	  for(int g = 0; g<2; g++){
	 
 
	    if(j!=3&&g==1)continue;
	    
	    for(int k = 1; k<26; k++){
	      if(background_diff_rebin[g][i][j][l][mc_code]->GetBinContent(k)>0) background_diff_noerr_up[g][i][j][l][mc_code]->SetBinContent(k, background_diff_rebin[g][i][j][l][mc_code]->GetBinContent(k));
	      else background_diff_noerr_down[g][i][j][l][mc_code]->SetBinContent(k, background_diff_rebin[g][i][j][l][mc_code]->GetBinContent(k));

	      if( signal_dPhi_rebin[g+2][i][j][l][mc_code]->GetBinContent(k)>0) signal_dPhi_noerr_up[g+2][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[g+2][i][j][l][mc_code]->GetBinContent(k));
	      else signal_dPhi_noerr_down[g+2][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[g+2][i][j][l][mc_code]->GetBinContent(k));

	      if(signal_dPhi_rebin[g+4][i][j][l][mc_code]->GetBinContent(k)>0)   signal_dPhi_noerr_up[g+4][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[g+4][i][j][l][mc_code]->GetBinContent(k));
	      else   signal_dPhi_noerr_down[g+4][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[g+4][i][j][l][mc_code]->GetBinContent(k));
	
	      if(sub_lead_dPhi_rebin[g][i][j][l][mc_code]->GetBinContent(k)>0) sub_lead_dPhi_noerr_up[g][i][j][l][mc_code]->SetBinContent(k,sub_lead_dPhi_rebin[g][i][j][l][mc_code]->GetBinContent(k));
	      else  sub_lead_dPhi_noerr_down[g][i][j][l][mc_code]->SetBinContent(k,sub_lead_dPhi_rebin[g][i][j][l][mc_code]->GetBinContent(k));
	    }

	    switch(i){
	    case 0:
	      signal_dPhi_noerr_up[g+2][i][j][l][mc_code]->SetFillColor(kBlue-9);
	      signal_dPhi_noerr_up[g+4][i][j][l][mc_code]->SetFillColor(kBlue-9);
	      background_diff_noerr_up[g][i][j][l][mc_code]->SetFillColor(kBlue-9);
	      sub_lead_dPhi_noerr_up[g][i][j][l][mc_code]->SetFillColor(kBlue-9);

	      signal_dPhi_noerr_down[g+2][i][j][l][mc_code]->SetFillColor(kBlue-9);
	      signal_dPhi_noerr_down[g+4][i][j][l][mc_code]->SetFillColor(kBlue-9);
	      background_diff_noerr_down[g][i][j][l][mc_code]->SetFillColor(kBlue-9);
	      sub_lead_dPhi_noerr_down[g][i][j][l][mc_code]->SetFillColor(kBlue-9);

	      break;
	    case 1:
	      signal_dPhi_noerr_up[g+2][i][j][l][mc_code]->SetFillColor(kYellow-9);
	      signal_dPhi_noerr_up[g+4][i][j][l][mc_code]->SetFillColor(kYellow-9);
	      background_diff_noerr_up[g][i][j][l][mc_code]->SetFillColor(kYellow-9);
	      sub_lead_dPhi_noerr_up[g][i][j][l][mc_code]->SetFillColor(kYellow-9);

	      signal_dPhi_noerr_down[g+2][i][j][l][mc_code]->SetFillColor(kYellow-9);
	      signal_dPhi_noerr_down[g+4][i][j][l][mc_code]->SetFillColor(kYellow-9);
	      background_diff_noerr_down[g][i][j][l][mc_code]->SetFillColor(kYellow-9);
	      sub_lead_dPhi_noerr_down[g][i][j][l][mc_code]->SetFillColor(kYellow-9);
	      break;
	    case 2:
	      signal_dPhi_noerr_up[g+2][i][j][l][mc_code]->SetFillColor(kOrange+1);
	      signal_dPhi_noerr_up[g+4][i][j][l][mc_code]->SetFillColor(kOrange+1);
	      background_diff_noerr_up[g][i][j][l][mc_code]->SetFillColor(kOrange+1);
	      sub_lead_dPhi_noerr_up[g][i][j][l][mc_code]->SetFillColor(kOrange+1);

	      signal_dPhi_noerr_down[g+2][i][j][l][mc_code]->SetFillColor(kOrange+1);
	      signal_dPhi_noerr_down[g+4][i][j][l][mc_code]->SetFillColor(kOrange+1);
	      background_diff_noerr_down[g][i][j][l][mc_code]->SetFillColor(kOrange+1);
	      sub_lead_dPhi_noerr_down[g][i][j][l][mc_code]->SetFillColor(kOrange+1);
	      break;
	    case 3:
	      signal_dPhi_noerr_up[g+2][i][j][l][mc_code]->SetFillColor(kViolet-5);
	      signal_dPhi_noerr_up[g+4][i][j][l][mc_code]->SetFillColor(kViolet-5);
	      background_diff_noerr_up[g][i][j][l][mc_code]->SetFillColor(kViolet-5);
	      sub_lead_dPhi_noerr_up[g][i][j][l][mc_code]->SetFillColor(kViolet-5);

	      signal_dPhi_noerr_down[g+2][i][j][l][mc_code]->SetFillColor(kViolet-5);
	      signal_dPhi_noerr_down[g+4][i][j][l][mc_code]->SetFillColor(kViolet-5);
	      background_diff_noerr_down[g][i][j][l][mc_code]->SetFillColor(kViolet-5);
	      sub_lead_dPhi_noerr_down[g][i][j][l][mc_code]->SetFillColor(kViolet-5);
	      break;
	    case 4:
	      signal_dPhi_noerr_up[g+2][i][j][l][mc_code]->SetFillColor(kGreen+3);
	      signal_dPhi_noerr_up[g+4][i][j][l][mc_code]->SetFillColor(kGreen+3);
	      background_diff_noerr_up[g][i][j][l][mc_code]->SetFillColor(kGreen+3);
	      sub_lead_dPhi_noerr_up[g][i][j][l][mc_code]->SetFillColor(kGreen+3);

	      signal_dPhi_noerr_down[g+2][i][j][l][mc_code]->SetFillColor(kGreen+3);
	      signal_dPhi_noerr_down[g+4][i][j][l][mc_code]->SetFillColor(kGreen+3);
	      background_diff_noerr_down[g][i][j][l][mc_code]->SetFillColor(kGreen+3);
	      sub_lead_dPhi_noerr_down[g][i][j][l][mc_code]->SetFillColor(kGreen+3);
	      break;
	    case 5:
	      signal_dPhi_noerr_up[g+2][i][j][l][mc_code]->SetFillColor(kRed+1);
	      signal_dPhi_noerr_up[g+4][i][j][l][mc_code]->SetFillColor(kRed+1);
	      background_diff_noerr_up[g][i][j][l][mc_code]->SetFillColor(kRed+1);
	      sub_lead_dPhi_noerr_up[g][i][j][l][mc_code]->SetFillColor(kRed+1);

	      signal_dPhi_noerr_down[g+2][i][j][l][mc_code]->SetFillColor(kRed+1);
	      signal_dPhi_noerr_down[g+4][i][j][l][mc_code]->SetFillColor(kRed+1);
	      background_diff_noerr_down[g][i][j][l][mc_code]->SetFillColor(kRed+1);
	      sub_lead_dPhi_noerr_down[g][i][j][l][mc_code]->SetFillColor(kRed+1);
	      break;
	    default:
	      break;
	    }

	    signal_dPhi_noerr_up[g+2][i][j][l][mc_code]->SetFillStyle(1001);
	    signal_dPhi_noerr_up[g+4][i][j][l][mc_code]->SetFillStyle(1001);
	    background_diff_noerr_up[g][i][j][l][mc_code]->SetFillStyle(1001);
	    sub_lead_dPhi_noerr_up[g][i][j][l][mc_code]->SetFillStyle(1001);
	  
	    signal_dPhi_noerr_down[g+2][i][j][l][mc_code]->SetFillStyle(1001);
	    signal_dPhi_noerr_down[g+4][i][j][l][mc_code]->SetFillStyle(1001);
	    background_diff_noerr_down[g][i][j][l][mc_code]->SetFillStyle(1001);
	    sub_lead_dPhi_noerr_down[g][i][j][l][mc_code]->SetFillStyle(1001);
	  

	  }
	}
      }
    }



    cout<<"and now we start plotting...."<<endl;
    for(int l = 0; l<3; l++){
   
      c_jet[l][mc_code] = new TCanvas(Form("MomentumBalanceJet%d%d",l,mc_code),"",10,10,1200,1600);
      c_jet[l][mc_code]->Divide(3,4,0,0);

      c_longrange[l][mc_code] = new TCanvas(Form("MomentumBalanceLongRange%d%d",l,mc_code),"",10,10,1200,800);
      c_longrange[l][mc_code]->Divide(3,2,0,0);
  
      c_nobkgsub[l][mc_code] = new TCanvas(Form("MomentumBalanceNoBkgSub%d%d",l,mc_code),"",10,10,1200,800);
      c_nobkgsub[l][mc_code]->Divide(3,2,0,0);
    

      cout<<"here..."<<endl;


      for(int j = 0; j<4; j++){
	//	if(j==1)continue;
      
	tot_name = make_name("AllHists_Summed_",2,0,j,l,pTlabel,centlabel,Ajlabel);
	signal_dPhi_tot[2][0][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[2][0][j][l][mc_code]->Clone(tot_name);
    
	tot_name.ReplaceAll("SubLeading","Leading");
	signal_dPhi_tot[4][0][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[4][0][j][l][mc_code]->Clone(tot_name);

	tot_name = make_name("AllHists_BackgroundDiff_",0,0,j,l,pTlabel,centlabel,Ajlabel);
	background_diff_tot[0][0][j][l][mc_code] = (TH1D*)background_diff_rebin[0][0][j][l][mc_code]->Clone(tot_name);
   
	tot_name = make_name("AllHists_NoBkgSub_",0,0,j,l,pTlabel,centlabel,Ajlabel);
	sub_lead_dPhi_tot[0][0][j][l][mc_code] = (TH1D*)sub_lead_dPhi_rebin[0][0][j][l][mc_code]->Clone(tot_name);
    
	if(use_highpT_bin){
	  for(int k = 1; k<6; k++){
	    signal_dPhi_tot[2][0][j][l][mc_code]->Add(signal_dPhi_rebin[2][k][j][l][mc_code]);
	    signal_dPhi_tot[4][0][j][l][mc_code]->Add(signal_dPhi_rebin[4][k][j][l][mc_code]);
	    background_diff_tot[0][0][j][l][mc_code]->Add( background_diff_rebin[0][k][j][l][mc_code]);
	    sub_lead_dPhi_tot[0][0][j][l][mc_code]->Add(sub_lead_dPhi_rebin[0][k][j][l][mc_code]);
	  }
	}else{

	  for(int k = 1; k<5; k++){
	    signal_dPhi_tot[2][0][j][l][mc_code]->Add(signal_dPhi_rebin[2][k][j][l][mc_code]);
	    signal_dPhi_tot[4][0][j][l][mc_code]->Add(signal_dPhi_rebin[4][k][j][l][mc_code]);
	    background_diff_tot[0][0][j][l][mc_code]->Add( background_diff_rebin[0][k][j][l][mc_code]);
	    sub_lead_dPhi_tot[0][0][j][l][mc_code]->Add(sub_lead_dPhi_rebin[0][k][j][l][mc_code]);
	  }




	}
	signal_dPhi_tot[2][0][j][l][mc_code]->SetMarkerStyle(4);
	signal_dPhi_tot[4][0][j][l][mc_code]->SetMarkerStyle(4);
	background_diff_tot[0][0][j][l][mc_code]->SetMarkerStyle(4);
	sub_lead_dPhi_tot[0][0][j][l][mc_code]->SetMarkerStyle(4);

      }

      tot_name.ReplaceAll("PbPb","pp");
      signal_dPhi_tot[5][0][3][l][mc_code] = (TH1D*)signal_dPhi_rebin[5][0][3][l][mc_code]->Clone(tot_name);

      tot_name.ReplaceAll("Leading","SubLeading");
      signal_dPhi_tot[3][0][3][l][mc_code] = (TH1D*)signal_dPhi_rebin[3][0][3][l][mc_code]->Clone(tot_name);


      tot_name = make_name("AllHists_BackgroundDiff_",1,0,3,l,pTlabel,centlabel,Ajlabel);
      background_diff_tot[1][0][3][l][mc_code] = (TH1D*)background_diff_rebin[1][0][3][l][mc_code]->Clone(tot_name);
   
      tot_name = make_name("AllHists_NoBkgSub_",1,0,3,l,pTlabel,centlabel,Ajlabel);
      sub_lead_dPhi_tot[1][0][3][l][mc_code] = (TH1D*)sub_lead_dPhi_rebin[1][0][3][l][mc_code]->Clone(tot_name);
    
      if(use_highpT_bin){
	for(int k = 1; k<6; k++){
	  signal_dPhi_tot[3][0][3][l][mc_code]->Add(signal_dPhi_rebin[3][k][3][l][mc_code]);
	  signal_dPhi_tot[5][0][3][l][mc_code]->Add(signal_dPhi_rebin[5][k][3][l][mc_code]);
	  background_diff_tot[1][0][3][l][mc_code]->Add(background_diff_rebin[1][k][3][l][mc_code]);
	  sub_lead_dPhi_tot[1][0][3][l][mc_code]->Add( sub_lead_dPhi_rebin[1][k][3][l][mc_code]);
	}
      }else{
	for(int k = 1; k<5; k++){
	  signal_dPhi_tot[3][0][3][l][mc_code]->Add(signal_dPhi_rebin[3][k][3][l][mc_code]);
	  signal_dPhi_tot[5][0][3][l][mc_code]->Add(signal_dPhi_rebin[5][k][3][l][mc_code]);
	  background_diff_tot[1][0][3][l][mc_code]->Add(background_diff_rebin[1][k][3][l][mc_code]);
	  sub_lead_dPhi_tot[1][0][3][l][mc_code]->Add( sub_lead_dPhi_rebin[1][k][3][l][mc_code]);
	}
      }

      signal_dPhi_tot[3][0][3][l][mc_code]->SetMarkerStyle(4);
      signal_dPhi_tot[5][0][3][l][mc_code]->SetMarkerStyle(4);
      background_diff_tot[1][0][3][l][mc_code]->SetMarkerStyle(4);
      sub_lead_dPhi_tot[1][0][3][l][mc_code]->SetMarkerStyle(4);
 
      tot_name = make_name("AllHists_PbPb_pp_",2,0,0,l,pTlabel,centlabel,Ajlabel);
      signal_dPhi_tot[8][0][0][l][mc_code] = (TH1D*) signal_dPhi_tot[2][0][0][l][mc_code]->Clone(tot_name);
      signal_dPhi_tot[8][0][0][l][mc_code]->Add( signal_dPhi_tot[3][0][3][l][mc_code],-1. );

      tot_name.ReplaceAll("SubLeading","Leading");
      signal_dPhi_tot[10][0][0][l][mc_code] = (TH1D*) signal_dPhi_tot[4][0][0][l][mc_code]->Clone(tot_name);
      signal_dPhi_tot[10][0][0][l][mc_code]->Add( signal_dPhi_tot[5][0][3][l][mc_code],-1. );

      tot_name.ReplaceAll("Leading","DoubleDiff");
 
      signal_dPhi_tot[6][0][0][l][mc_code] = (TH1D*) signal_dPhi_tot[8][0][0][l][mc_code]->Clone(tot_name);
      signal_dPhi_tot[6][0][0][l][mc_code]->Add( signal_dPhi_tot[10][0][0][l][mc_code]);


      tot_name = make_name("BackgroundDiff_PbPb_pp_",0,0,0,l,pTlabel,centlabel,Ajlabel);
      background_diff_tot[7][0][0][l][mc_code] = (TH1D*) background_diff_tot[0][0][0][l][mc_code]->Clone(tot_name);
      background_diff_tot[7][0][0][l][mc_code]->Add( background_diff_tot[1][0][3][l][mc_code],-1. );



      tot_name = make_name("Sub_Lead_PbPb_pp_",2,0,0,l,pTlabel,centlabel,Ajlabel);
      sub_lead_dPhi_tot[6][0][0][l][mc_code] = (TH1D*) sub_lead_dPhi_tot[0][0][0][l][mc_code]->Clone(tot_name);
      sub_lead_dPhi_tot[6][0][0][l][mc_code]->Add( sub_lead_dPhi_tot[1][0][3][l][mc_code],-1. );

   
  
 



      tot_name = make_name("AllHists_PbPb_pp_",2,0,0,l,pTlabel,centlabel,Ajlabel);
      signal_dPhi_tot[8][0][2][l][mc_code] = (TH1D*) signal_dPhi_tot[2][0][2][l][mc_code]->Clone(tot_name);
      signal_dPhi_tot[8][0][2][l][mc_code]->Add( signal_dPhi_tot[3][0][3][l][mc_code],-1. );

      tot_name.ReplaceAll("SubLeading","Leading");
      signal_dPhi_tot[10][0][2][l][mc_code] = (TH1D*) signal_dPhi_tot[4][0][2][l][mc_code]->Clone(tot_name);
      signal_dPhi_tot[10][0][2][l][mc_code]->Add( signal_dPhi_tot[5][0][3][l][mc_code],-1. );

      tot_name.ReplaceAll("Leading","DoubleDiff");
 
      signal_dPhi_tot[6][0][2][l][mc_code] = (TH1D*) signal_dPhi_tot[8][0][2][l][mc_code]->Clone(tot_name);
      signal_dPhi_tot[6][0][2][l][mc_code]->Add( signal_dPhi_tot[10][0][2][l][mc_code] );


      tot_name = make_name("BackgroundDiff_PbPb_pp_",0,0,2,l,pTlabel,centlabel,Ajlabel);
      background_diff_tot[7][0][2][l][mc_code] = (TH1D*) background_diff_tot[0][0][2][l][mc_code]->Clone(tot_name);
      background_diff_tot[7][0][2][l][mc_code]->Add( background_diff_tot[1][0][3][l][mc_code],-1. );

   
      tot_name = make_name("Sub_Lead_PbPb_pp_",2,0,2,l,pTlabel,centlabel,Ajlabel);
      sub_lead_dPhi_tot[6][0][2][l][mc_code] = (TH1D*) sub_lead_dPhi_tot[0][0][2][l][mc_code]->Clone(tot_name);
      sub_lead_dPhi_tot[6][0][2][l][mc_code]->Add( sub_lead_dPhi_tot[1][0][3][l][mc_code],-1. );

   




     
      for(int j = 0; j<4; j++){
	//	if(j==1)continue;

	stack_name = make_name("AllHists_Up_",2,0,j,l,pTlabel,centlabel,Ajlabel);
	stack_name.ReplaceAll("TrkPt05_TrkPt1_","");

	dPhi_pTdist_up[2][j][l][mc_code] = new THStack(stack_name,"");
      
	stack_name.ReplaceAll("SubLeading","Leading");
	dPhi_pTdist_up[4][j][l][mc_code] = new THStack(stack_name,"");

	stack_name.ReplaceAll("Leading","Leading_PbPb_pp");
	dPhi_pTdist_up[10][j][l][mc_code] = new THStack(stack_name,"");

	stack_name.ReplaceAll("Leading","SubLeading");
	dPhi_pTdist_up[8][j][l][mc_code] = new THStack(stack_name,"");
      
	stack_name.ReplaceAll("SubLeading","DoubleDiff");
	dPhi_pTdist_up[6][j][l][mc_code] = new THStack(stack_name,"");

	stack_name.ReplaceAll("DoubleDiff","SideBandDiff");
	dPhi_pTdist_up[0][j][l][mc_code] = new THStack(stack_name,"");
      
	stack_name.ReplaceAll("SideBandDiff","SideBandDoubleDiff");
	dPhi_pTdist_up[7][j][l][mc_code] = new THStack(stack_name,"");

	stack_name.ReplaceAll("SideBandDoubleDiff","SubLeadingMinusLeading");
	dPhi_Sub_Lead_NoBkg_up[0][j][l][mc_code] = new THStack(stack_name,"");
    
	stack_name.ReplaceAll("SubLeadingMinusLeading","DoubleDiffNoBkgSub");
	dPhi_Sub_Lead_NoBkg_up[6][j][l][mc_code] = new THStack(stack_name,"");



	stack_name = make_name("AllHists_Down_",2,0,j,l,pTlabel,centlabel,Ajlabel);
	stack_name.ReplaceAll("TrkPt05_TrkPt1_","");

	dPhi_pTdist_down[2][j][l][mc_code] = new THStack(stack_name,"");
      
	stack_name.ReplaceAll("SubLeading","Leading");
	dPhi_pTdist_down[4][j][l][mc_code] = new THStack(stack_name,"");

	stack_name.ReplaceAll("Leading","Leading_PbPb_pp");
	dPhi_pTdist_down[10][j][l][mc_code] = new THStack(stack_name,"");

	stack_name.ReplaceAll("Leading","SubLeading");
	dPhi_pTdist_down[8][j][l][mc_code] = new THStack(stack_name,"");
      
	stack_name.ReplaceAll("SubLeading","DoubleDiff");
	dPhi_pTdist_down[6][j][l][mc_code] = new THStack(stack_name,"");

	stack_name.ReplaceAll("DoubleDiff","SideBandDiff");
	dPhi_pTdist_down[0][j][l][mc_code] = new THStack(stack_name,"");
      
	stack_name.ReplaceAll("SideBandDiff","SideBandDoubleDiff");
	dPhi_pTdist_down[7][j][l][mc_code] = new THStack(stack_name,"");

	stack_name.ReplaceAll("SideBandDoubleDiff","SubLeadingMinusLeading");
	dPhi_Sub_Lead_NoBkg_down[0][j][l][mc_code] = new THStack(stack_name,"");
    
	stack_name.ReplaceAll("SubLeadingMinusLeading","DoubleDiffNoBkgSub");
	dPhi_Sub_Lead_NoBkg_down[6][j][l][mc_code] = new THStack(stack_name,"");

    
	for(int i = 0; i<6; i++){
	  in_name = make_name("PbPb_pp_WithErrors_",8,i,j,l,pTlabel,centlabel,Ajlabel);
	  signal_dPhi_rebin[8][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[2][i][j][l][mc_code]->Clone(in_name);
	  signal_dPhi_rebin[8][i][j][l][mc_code]->Add(signal_dPhi_rebin[3][i][3][l][mc_code],-1.);

	  in_name.ReplaceAll("SubLeading","Leading");
	  signal_dPhi_rebin[10][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[4][i][j][l][mc_code]->Clone(in_name);
	  signal_dPhi_rebin[10][i][j][l][mc_code]->Add(signal_dPhi_rebin[5][i][3][l][mc_code],-1.);

	  in_name.ReplaceAll("SubLeading","DoubleDiff");
	  signal_dPhi_rebin[6][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[8][i][j][l][mc_code]->Clone(in_name);
	  signal_dPhi_rebin[6][i][j][l][mc_code]->Add(signal_dPhi_rebin[10][i][j][l][mc_code]);


	  in_name.ReplaceAll("DoubleDiff","SideBandDoubleDiff");
	  background_diff_rebin[7][i][j][l][mc_code] = (TH1D*)background_diff_rebin[0][i][j][l][mc_code]->Clone(in_name);
	  background_diff_rebin[7][i][j][l][mc_code]->Add(background_diff_rebin[1][i][3][l][mc_code],-1.);

	  in_name.ReplaceAll("SideBand","NoBkgSub");
	  sub_lead_dPhi_rebin[6][i][j][l][mc_code] = (TH1D*)sub_lead_dPhi_rebin[0][i][j][l][mc_code]->Clone(in_name);
	  sub_lead_dPhi_rebin[6][i][j][l][mc_code]->Add(sub_lead_dPhi_rebin[1][i][3][l][mc_code],-1.);
  

	  //Now make NoErr versions


	  in_name = make_name("PbPb_pp_Up_",8,i,j,l,pTlabel,centlabel,Ajlabel);
	  signal_dPhi_noerr_up[8][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[8][i][j][l][mc_code]->Clone(in_name);

	  in_name.ReplaceAll("SubLeading","Leading");
	  signal_dPhi_noerr_up[10][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[10][i][j][l][mc_code]->Clone(in_name);

	  in_name.ReplaceAll("SubLeading","DoubleDiff");
	  signal_dPhi_noerr_up[6][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[6][i][j][l][mc_code]->Clone(in_name);

	  in_name.ReplaceAll("DoubleDiff","SideBandDoubleDiff");
	  background_diff_noerr_up[7][i][j][l][mc_code] = (TH1D*)background_diff_rebin[7][i][j][l][mc_code]->Clone(in_name);

	  in_name.ReplaceAll("SideBand","NoBkgSub");
	  sub_lead_dPhi_noerr_up[6][i][j][l][mc_code] = (TH1D*)sub_lead_dPhi_rebin[6][i][j][l][mc_code]->Clone(in_name);


	  in_name = make_name("PbPb_pp_Down_",8,i,j,l,pTlabel,centlabel,Ajlabel);
	  signal_dPhi_noerr_down[8][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[8][i][j][l][mc_code]->Clone(in_name);

	  in_name.ReplaceAll("SubLeading","Leading");
	  signal_dPhi_noerr_down[10][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[10][i][j][l][mc_code]->Clone(in_name);

	  in_name.ReplaceAll("SubLeading","DoubleDiff");
	  signal_dPhi_noerr_down[6][i][j][l][mc_code] = (TH1D*)signal_dPhi_rebin[6][i][j][l][mc_code]->Clone(in_name);


	  in_name.ReplaceAll("DoubleDiff","SideBandDoubleDiff");
	  background_diff_noerr_down[7][i][j][l][mc_code] = (TH1D*)background_diff_rebin[7][i][j][l][mc_code]->Clone(in_name);


	  in_name.ReplaceAll("SideBand","NoBkgSub");
	  sub_lead_dPhi_noerr_down[6][i][j][l][mc_code] = (TH1D*)sub_lead_dPhi_rebin[6][i][j][l][mc_code]->Clone(in_name);
      
	  for(int k = 1; k<26; k++){

	    background_diff_noerr_up[7][i][j][l][mc_code]->SetBinError(k,0.);
	    background_diff_noerr_down[7][i][j][l][mc_code]->SetBinError(k,0.);
	    signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetBinError(k,0.);
	    signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetBinError(k,0.);
	    signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetBinError(k,0.);
	    signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetBinError(k,0.);
	    signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetBinError(k,0.);
	    signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetBinError(k,0.);
	    sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetBinError(k,0.);
	    sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetBinError(k,0.);


	    if(background_diff_rebin[7][i][j][l][mc_code]->GetBinContent(k)>0){
	      background_diff_noerr_up[7][i][j][l][mc_code]->SetBinContent(k, background_diff_rebin[7][i][j][l][mc_code]->GetBinContent(k));
	      background_diff_noerr_down[7][i][j][l][mc_code]->SetBinContent(k,0.);
	    } else{
	      background_diff_noerr_down[7][i][j][l][mc_code]->SetBinContent(k, background_diff_rebin[7][i][j][l][mc_code]->GetBinContent(k));
	      background_diff_noerr_up[7][i][j][l][mc_code]->SetBinContent(k, 0.);
	    }

	    if(	signal_dPhi_rebin[8][i][j][l][mc_code]->GetBinContent(k)>0){
	      signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[8][i][j][l][mc_code]->GetBinContent(k));
	      signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetBinContent(k,0.);
	    }else{
	      signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[8][i][j][l][mc_code]->GetBinContent(k));
	      signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetBinContent(k,0.);
	    }

	    if(signal_dPhi_rebin[10][i][j][l][mc_code]->GetBinContent(k)>0){
	      signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[10][i][j][l][mc_code]->GetBinContent(k));
	      signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetBinContent(k,0.);
	    }else{ 
	      signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[10][i][j][l][mc_code]->GetBinContent(k));
	      signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetBinContent(k,0.);
	    }



	    if(signal_dPhi_rebin[6][i][j][l][mc_code]->GetBinContent(k)>0){
	      signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[6][i][j][l][mc_code]->GetBinContent(k));
	      signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetBinContent(k,0.);
	    }else{ 
	      signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetBinContent(k,signal_dPhi_rebin[6][i][j][l][mc_code]->GetBinContent(k));
	      signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetBinContent(k,0.);
	    }



	    if(sub_lead_dPhi_rebin[6][i][j][l][mc_code]->GetBinContent(k)>0){
	      sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetBinContent(k,sub_lead_dPhi_rebin[6][i][j][l][mc_code]->GetBinContent(k));
	      sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetBinContent(k,0.);
	    }else{
	      sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetBinContent(k,sub_lead_dPhi_rebin[6][i][j][l][mc_code]->GetBinContent(k));
	      sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetBinContent(k,0.);
	    }
	  }


	  switch(i){
	  case 0:
	    signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kBlue-9);
	    background_diff_noerr_up[7][i][j][l][mc_code]->SetFillColor(kBlue-9);
	    sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kBlue-9);

	    signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kBlue-9);
	    background_diff_noerr_down[7][i][j][l][mc_code]->SetFillColor(kBlue-9);
	    sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kBlue-9);

	    break;
	  case 1:
	    signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    background_diff_noerr_up[7][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kYellow-9);

	    signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    background_diff_noerr_down[7][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kYellow-9);
	    break;
	  case 2:
	    signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    background_diff_noerr_up[7][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kOrange+1);

	    signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    background_diff_noerr_down[7][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kOrange+1);
	    break;
	  case 3:
	    signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    background_diff_noerr_up[7][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kViolet-5);

	    signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    background_diff_noerr_down[7][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kViolet-5);
	    break;
	  case 4:
	    signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    background_diff_noerr_up[7][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kGreen+3);

	    signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    background_diff_noerr_down[7][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kGreen+3);
	    break;
	  case 5:
	    signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kRed+1);
	    background_diff_noerr_up[7][i][j][l][mc_code]->SetFillColor(kRed+1);
	    sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillColor(kRed+1);

	    signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kRed+1);
	    background_diff_noerr_down[7][i][j][l][mc_code]->SetFillColor(kRed+1);
	    sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillColor(kRed+1);
	    break;
	  default:
	    break;
	  }

	  signal_dPhi_noerr_up[8][i][j][l][mc_code]->SetFillStyle(1001);
	  signal_dPhi_noerr_up[10][i][j][l][mc_code]->SetFillStyle(1001);
	  signal_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillStyle(1001);
	  background_diff_noerr_up[7][i][j][l][mc_code]->SetFillStyle(1001);
	  sub_lead_dPhi_noerr_up[6][i][j][l][mc_code]->SetFillStyle(1001);
	  
	  signal_dPhi_noerr_down[8][i][j][l][mc_code]->SetFillStyle(1001);
	  signal_dPhi_noerr_down[10][i][j][l][mc_code]->SetFillStyle(1001);
	  signal_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillStyle(1001);
	  background_diff_noerr_down[7][i][j][l][mc_code]->SetFillStyle(1001);
	  sub_lead_dPhi_noerr_down[6][i][j][l][mc_code]->SetFillStyle(1001);
	  
	}
      
	if(use_highpT_bin){

	  for(int i = 1; i<6; i++){

	    dPhi_pTdist_up[2][j][l][mc_code]->Add(signal_dPhi_noerr_up[2][6-i][j][l][mc_code]);
	    dPhi_pTdist_up[4][j][l][mc_code]->Add(signal_dPhi_noerr_up[4][6-i][j][l][mc_code]);
	    dPhi_pTdist_up[8][j][l][mc_code]->Add(signal_dPhi_noerr_up[8][6-i][j][l][mc_code]);
	    dPhi_pTdist_up[10][j][l][mc_code]->Add(signal_dPhi_noerr_up[10][6-i][j][l][mc_code]);  
	    dPhi_pTdist_up[6][j][l][mc_code]->Add(signal_dPhi_noerr_up[6][6-i][j][l][mc_code]);
	    dPhi_pTdist_up[0][j][l][mc_code]->Add(background_diff_noerr_up[0][6-i][j][l][mc_code]);
	    dPhi_pTdist_up[7][j][l][mc_code]->Add(background_diff_noerr_up[7][6-i][j][l][mc_code]);
	    dPhi_Sub_Lead_NoBkg_up[0][j][l][mc_code]->Add(sub_lead_dPhi_noerr_up[0][6-i][j][l][mc_code]);
	    dPhi_Sub_Lead_NoBkg_up[6][j][l][mc_code]->Add(sub_lead_dPhi_noerr_up[6][6-i][j][l][mc_code]);

	    dPhi_pTdist_down[2][j][l][mc_code]->Add(signal_dPhi_noerr_down[2][6-i][j][l][mc_code]);
	    dPhi_pTdist_down[4][j][l][mc_code]->Add(signal_dPhi_noerr_down[4][6-i][j][l][mc_code]);
	    dPhi_pTdist_down[8][j][l][mc_code]->Add(signal_dPhi_noerr_down[8][6-i][j][l][mc_code]);
	    dPhi_pTdist_down[10][j][l][mc_code]->Add(signal_dPhi_noerr_down[10][6-i][j][l][mc_code]);  
	    dPhi_pTdist_down[6][j][l][mc_code]->Add(signal_dPhi_noerr_down[6][6-i][j][l][mc_code]);
	    dPhi_pTdist_down[0][j][l][mc_code]->Add(background_diff_noerr_down[0][6-i][j][l][mc_code]);
	    dPhi_pTdist_down[7][j][l][mc_code]->Add(background_diff_noerr_down[7][6-i][j][l][mc_code]);
	    dPhi_Sub_Lead_NoBkg_down[0][j][l][mc_code]->Add(sub_lead_dPhi_noerr_down[0][6-i][j][l][mc_code]);
	    dPhi_Sub_Lead_NoBkg_down[6][j][l][mc_code]->Add(sub_lead_dPhi_noerr_down[6][6-i][j][l][mc_code]);


	  }


	}else{

	  cout<<"make one set of stacks"<<endl;
	  for(int i = 1; i<6; i++){

	    dPhi_pTdist_up[2][j][l][mc_code]->Add(signal_dPhi_noerr_up[2][5-i][j][l][mc_code]);
	    dPhi_pTdist_up[4][j][l][mc_code]->Add(signal_dPhi_noerr_up[4][5-i][j][l][mc_code]);
	    dPhi_pTdist_up[8][j][l][mc_code]->Add(signal_dPhi_noerr_up[8][5-i][j][l][mc_code]);
	    dPhi_pTdist_up[10][j][l][mc_code]->Add(signal_dPhi_noerr_up[10][5-i][j][l][mc_code]);  
	    dPhi_pTdist_up[6][j][l][mc_code]->Add(signal_dPhi_noerr_up[6][5-i][j][l][mc_code]);
	    dPhi_pTdist_up[0][j][l][mc_code]->Add(background_diff_noerr_up[0][5-i][j][l][mc_code]);
	    dPhi_pTdist_up[7][j][l][mc_code]->Add(background_diff_noerr_up[7][5-i][j][l][mc_code]);
	    dPhi_Sub_Lead_NoBkg_up[0][j][l][mc_code]->Add(sub_lead_dPhi_noerr_up[0][5-i][j][l][mc_code]);
	    dPhi_Sub_Lead_NoBkg_up[6][j][l][mc_code]->Add(sub_lead_dPhi_noerr_up[6][5-i][j][l][mc_code]);

	    dPhi_pTdist_down[2][j][l][mc_code]->Add(signal_dPhi_noerr_down[2][5-i][j][l][mc_code]);
	    dPhi_pTdist_down[4][j][l][mc_code]->Add(signal_dPhi_noerr_down[4][5-i][j][l][mc_code]);
	    dPhi_pTdist_down[8][j][l][mc_code]->Add(signal_dPhi_noerr_down[8][5-i][j][l][mc_code]);
	    dPhi_pTdist_down[10][j][l][mc_code]->Add(signal_dPhi_noerr_down[10][5-i][j][l][mc_code]);  
	    dPhi_pTdist_down[6][j][l][mc_code]->Add(signal_dPhi_noerr_down[6][5-i][j][l][mc_code]);
	    dPhi_pTdist_down[0][j][l][mc_code]->Add(background_diff_noerr_down[0][5-i][j][l][mc_code]);
	    dPhi_pTdist_down[7][j][l][mc_code]->Add(background_diff_noerr_down[7][5-i][j][l][mc_code]);
	    dPhi_Sub_Lead_NoBkg_down[0][j][l][mc_code]->Add(sub_lead_dPhi_noerr_down[0][5-i][j][l][mc_code]);
	    dPhi_Sub_Lead_NoBkg_down[6][j][l][mc_code]->Add(sub_lead_dPhi_noerr_down[6][5-i][j][l][mc_code]);


	  }
	}
       signal_max = 104.;
      signal_min = -104.;
      

      diff_max = 16.;
      diff_min = -16.;
      
      double_diff_max = 11.;
      double_diff_min = -11.; 
      
      no_bgsub_max = 25.;
      no_bgsub_min = -25.;

      if(use_highpT_bin){
	signal_max = 440.;
	signal_min = -440.;
      

	diff_max = 38.;
	diff_min = -38.;
      
	double_diff_max = 38.;
	double_diff_min = -38.; 
      
	no_bgsub_max = 191.;
	no_bgsub_min = -191.;
      }



	dPhi_pTdist_up[2][j][l][mc_code]->SetMaximum(signal_max);
	dPhi_pTdist_up[2][j][l][mc_code]->SetMinimum(signal_min);


	dPhi_pTdist_up[8][j][l][mc_code]->SetMaximum(diff_max);
	dPhi_pTdist_up[8][j][l][mc_code]->SetMinimum(diff_min);


	dPhi_pTdist_up[10][j][l][mc_code]->SetMaximum(diff_max);
	dPhi_pTdist_up[10][j][l][mc_code]->SetMinimum(diff_min);

	dPhi_Sub_Lead_NoBkg_up[0][j][l][mc_code]->SetMaximum(no_bgsub_max);
	dPhi_Sub_Lead_NoBkg_up[0][j][l][mc_code]->SetMinimum(no_bgsub_min);
  
	dPhi_pTdist_up[6][j][l][mc_code]->SetMaximum(double_diff_max);
	dPhi_pTdist_up[6][j][l][mc_code]->SetMinimum(double_diff_min);

   
	dPhi_pTdist_up[0][j][l][mc_code]->SetMaximum(diff_max);
	dPhi_pTdist_up[0][j][l][mc_code]->SetMinimum(diff_min);
     

	dPhi_pTdist_up[7][j][l][mc_code]->SetMaximum(double_diff_max);
	dPhi_pTdist_up[7][j][l][mc_code]->SetMinimum(double_diff_min);

	dPhi_Sub_Lead_NoBkg_up[6][j][l][mc_code]->SetMaximum(double_diff_max);
	dPhi_Sub_Lead_NoBkg_up[6][j][l][mc_code]->SetMinimum(double_diff_min);



	dPhi_pTdist_up[2][j][l][mc_code]->SetMaximum(signal_max);
	dPhi_pTdist_up[2][j][l][mc_code]->SetMinimum(signal_min);


	dPhi_pTdist_up[8][j][l][mc_code]->SetMaximum(diff_max);
	dPhi_pTdist_up[8][j][l][mc_code]->SetMinimum(diff_min);


	dPhi_pTdist_down[10][j][l][mc_code]->SetMaximum(diff_max);
	dPhi_pTdist_down[10][j][l][mc_code]->SetMinimum(diff_min);

	dPhi_Sub_Lead_NoBkg_down[0][j][l][mc_code]->SetMaximum(no_bgsub_max);
	dPhi_Sub_Lead_NoBkg_down[0][j][l][mc_code]->SetMinimum(no_bgsub_min);
  
	dPhi_pTdist_down[6][j][l][mc_code]->SetMaximum(double_diff_max);
	dPhi_pTdist_down[6][j][l][mc_code]->SetMinimum(double_diff_min);

   
	dPhi_pTdist_down[0][j][l][mc_code]->SetMaximum(diff_max);
	dPhi_pTdist_down[0][j][l][mc_code]->SetMinimum(diff_min);
     

	dPhi_pTdist_down[7][j][l][mc_code]->SetMaximum(diff_max);
	dPhi_pTdist_down[7][j][l][mc_code]->SetMinimum(diff_min);

	dPhi_Sub_Lead_NoBkg_down[6][j][l][mc_code]->SetMaximum(no_bgsub_min);
	dPhi_Sub_Lead_NoBkg_down[6][j][l][mc_code]->SetMinimum(no_bgsub_min);
  
  
      
      }
   

      stack_name = make_name("AllHists_Up_",5,0,3,l,pTlabel,centlabel,Ajlabel);
      dPhi_pTdist_up[5][3][l][mc_code] = new THStack(stack_name,"");
      
      stack_name.ReplaceAll("Leading","SubLeading");
      dPhi_pTdist_up[3][3][l][mc_code] = new THStack(stack_name,"");

      stack_name.ReplaceAll("SubLeading","SideBandDiff");
      dPhi_pTdist_up[1][3][l][mc_code] = new THStack(stack_name,"");

      stack_name.ReplaceAll("SideBandDiff","SubLeadingMinusLeading");
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code] = new THStack(stack_name,"");
 

      stack_name = make_name("AllHists_Down_",5,0,3,l,pTlabel,centlabel,Ajlabel);
      dPhi_pTdist_down[5][3][l][mc_code] = new THStack(stack_name,"");
      
      stack_name.ReplaceAll("Leading","SubLeading");
      dPhi_pTdist_down[3][3][l][mc_code] = new THStack(stack_name,"");

      stack_name.ReplaceAll("SubLeading","SideBandDiff");
      dPhi_pTdist_down[1][3][l][mc_code] = new THStack(stack_name,"");

      stack_name.ReplaceAll("SideBandDiff","SubLeadingMinusLeading");
      dPhi_Sub_Lead_NoBkg_down[1][3][l][mc_code] = new THStack(stack_name,"");
 
      if(use_highpT_bin){
	for(int i = 1; i<7; i++){

	  dPhi_pTdist_up[3][3][l][mc_code]->Add(signal_dPhi_noerr_up[3][6-i][3][l][mc_code]);
	  dPhi_pTdist_up[5][3][l][mc_code]->Add(signal_dPhi_noerr_up[5][6-i][3][l][mc_code]);
	  dPhi_pTdist_up[1][3][l][mc_code]->Add(background_diff_noerr_up[1][6-i][3][l][mc_code]);
	  dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->Add(sub_lead_dPhi_noerr_up[1][6-i][3][l][mc_code]);

	  dPhi_pTdist_down[3][3][l][mc_code]->Add(signal_dPhi_noerr_down[3][6-i][3][l][mc_code]);
	  dPhi_pTdist_down[5][3][l][mc_code]->Add(signal_dPhi_noerr_down[5][6-i][3][l][mc_code]);
	  dPhi_pTdist_down[1][3][l][mc_code]->Add(background_diff_noerr_down[1][6-i][3][l][mc_code]);
	  dPhi_Sub_Lead_NoBkg_down[1][3][l][mc_code]->Add(sub_lead_dPhi_noerr_down[1][6-i][3][l][mc_code]);
  
	}

      }else{

	for(int i = 1; i<6; i++){

	  dPhi_pTdist_up[3][3][l][mc_code]->Add(signal_dPhi_noerr_up[3][5-i][3][l][mc_code]);
	  dPhi_pTdist_up[5][3][l][mc_code]->Add(signal_dPhi_noerr_up[5][5-i][3][l][mc_code]);
	  dPhi_pTdist_up[1][3][l][mc_code]->Add(background_diff_noerr_up[1][5-i][3][l][mc_code]);
	  dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->Add(sub_lead_dPhi_noerr_up[1][5-i][3][l][mc_code]);

	  dPhi_pTdist_down[3][3][l][mc_code]->Add(signal_dPhi_noerr_down[3][5-i][3][l][mc_code]);
	  dPhi_pTdist_down[5][3][l][mc_code]->Add(signal_dPhi_noerr_down[5][5-i][3][l][mc_code]);
	  dPhi_pTdist_down[1][3][l][mc_code]->Add(background_diff_noerr_down[1][5-i][3][l][mc_code]);
	  dPhi_Sub_Lead_NoBkg_down[1][3][l][mc_code]->Add(sub_lead_dPhi_noerr_down[1][5-i][3][l][mc_code]);
  
	}

      }

      dPhi_pTdist_up[3][3][l][mc_code]->SetMaximum(signal_max);
      dPhi_pTdist_up[3][3][l][mc_code]->SetMinimum(signal_min);

      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->SetMinimum(no_bgsub_min);
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->SetMaximum(no_bgsub_max);

      dPhi_pTdist_up[1][3][l][mc_code]->SetMaximum(diff_max);
      dPhi_pTdist_up[1][3][l][mc_code]->SetMinimum(diff_min);

      dPhi_pTdist_down[3][3][l][mc_code]->SetMaximum(signal_max);
      dPhi_pTdist_down[3][3][l][mc_code]->SetMinimum(signal_min);

      dPhi_Sub_Lead_NoBkg_down[1][3][l][mc_code]->SetMinimum(diff_min);
      dPhi_Sub_Lead_NoBkg_down[1][3][l][mc_code]->SetMaximum(diff_max);

      dPhi_pTdist_down[1][3][l][mc_code]->SetMaximum(diff_max);
      dPhi_pTdist_down[1][3][l][mc_code]->SetMinimum(diff_min);

 

      ///////////////////////////////////////////////////////////////

      //     DRAWING STARTS HERE  

      ///////////////////////////////////////////////////////////////
    
      cout<<"Drawing now!"<<endl;

      c_jet[l][mc_code]->cd(0);

    TLatex *aj_tex = new TLatex(0.07,0.98,Ajlabel);
    aj_tex->SetTextSize(0.025);
    aj_tex->SetLineColor(kWhite);
    aj_tex->SetNDC();
    aj_tex->Draw();
    
    TLatex *type_tex = new TLatex(0.3,0.98,"Jet Peak, |#Delta#eta| < 2.5");
    type_tex->SetTextSize(0.025);
    type_tex->SetLineColor(kWhite);
    type_tex->SetNDC();
    type_tex->Draw();
   
    TLatex   *MC_tex = new TLatex(0.07,0.96,"PYTHIA 6.423 Tune Z2");
    MC_tex->SetTextFont(43);
    MC_tex->SetTextSizePixels(25);
    MC_tex->SetLineColor(kWhite);
    MC_tex->SetNDC();
    MC_tex->Draw();

    TLatex   *Hyd_tex = new TLatex(0.3,0.96,"HYDJET v1.8");
    Hyd_tex->SetTextFont(43);
    Hyd_tex->SetTextSizePixels(25);
    Hyd_tex->SetLineColor(kWhite);
    Hyd_tex->SetNDC();
    Hyd_tex->Draw();
   
    TLatex   *jet_reco_tex = new TLatex(0.605,0.98,"anti-kT R = 0.3, |#eta_{jet}| < 1.6");
    jet_reco_tex->SetTextFont(43);
    jet_reco_tex->SetTextSizePixels(25);
    jet_reco_tex->SetLineColor(kWhite);
    jet_reco_tex->SetNDC();
    jet_reco_tex->Draw();

    TLatex   *jet_cut_tex = new TLatex(0.605,0.96,"120 < p_{T,1}< 300, p_{T,2}> 50 GeV/c, #Delta#phi_{1,2}> 5#pi/6");
    jet_cut_tex->SetTextFont(43);
    jet_cut_tex->SetTextSizePixels(25);
    jet_cut_tex->SetLineColor(kWhite);
    jet_cut_tex->SetNDC();
    jet_cut_tex->Draw();
    


   

      c_jet[l][mc_code]->cd(4);


      TLegend *legend = new TLegend(0.1,0.3,0.9,1.);
      legend->AddEntry(signal_dPhi_noerr_up[2][0][0][l][mc_code],"0.5<p_{T}^{assoc.}<1 GeV/c","f");
      legend->AddEntry(signal_dPhi_noerr_up[2][1][0][l][mc_code],"1<p_{T}^{assoc.}<2 GeV/c","f");
      legend->AddEntry(signal_dPhi_noerr_up[2][2][0][l][mc_code],"2<p_{T}^{assoc.}<3 GeV/c","f");
      legend->AddEntry(signal_dPhi_noerr_up[2][3][0][l][mc_code],"3<p_{T}^{assoc.}<4 GeV/c","f");
      legend->AddEntry(signal_dPhi_noerr_up[2][4][0][l][mc_code],"4<p_{T}^{assoc.}<8 GeV/c","f");
      if(use_highpT_bin)      legend->AddEntry(signal_dPhi_noerr_up[2][5][0][l][mc_code],"p_{T}^{assoc.}>8 GeV/c","f");
      legend->SetTextSize(0.055);
      legend->SetLineColor(kWhite);
      legend->Draw();



  
      c_jet[l][mc_code]->cd(1);
 
      dPhi_pTdist_up[3][3][l][mc_code]->Draw();
      dPhi_pTdist_up[3][3][l][mc_code]->GetYaxis()->SetLabelSize(0.07); 
      dPhi_pTdist_up[3][3][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
      dPhi_pTdist_up[3][3][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_pTdist_up[3][3][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[3][3][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_pTdist_up[3][3][l][mc_code]->GetYaxis()->SetTitle("1/N_{evt} dp_{T}/d#Delta#phi (GeV/c)");
      dPhi_pTdist_up[3][3][l][mc_code]->GetYaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[3][3][l][mc_code]->GetYaxis()->CenterTitle();
      dPhi_pTdist_up[5][3][l][mc_code]->Draw("same");
 
      dPhi_pTdist_down[5][3][l][mc_code]->Draw("same");
      dPhi_pTdist_down[3][3][l][mc_code]->Draw("same");
  
  
      if(mc_code==5) signal_dPhi_syst[3][0][3][l]->Draw("same e2");
      if(mc_code==5) signal_dPhi_syst[5][0][3][l]->Draw("same e2");
      signal_dPhi_tot[3][0][3][l][mc_code]->Draw("same");
      signal_dPhi_tot[5][0][3][l][mc_code]->Draw("same");
 
      cout<<"0"<<endl;
 
      label_pp = new TLatex(0.2,0.9,"PYTHIA");
      label_pp->SetTextSize(0.07);
      label_pp->SetLineColor(kWhite);
      label_pp->SetNDC();
      label_pp->Draw();

      TLine *l_phi = new TLine(-TMath::Pi()/2.,0.,TMath::Pi()/2.,0.);
      l_phi->SetLineStyle(2);
      l_phi->Draw();

      TLatex *prelim_tex_dphi = new TLatex(0.35,0.05,"Preliminary Simulation");
      prelim_tex_dphi->SetTextFont(53);
      prelim_tex_dphi->SetTextSizePixels(25);
      prelim_tex_dphi->SetLineColor(kWhite);
      prelim_tex_dphi->SetNDC();
      prelim_tex_dphi->Draw(); 
  
  

      TLatex *cms_tex_dphi = new TLatex(0.2,0.05,"CMS");
      cms_tex_dphi->SetTextSize(0.07);
      cms_tex_dphi->SetLineColor(kWhite);
      cms_tex_dphi->SetNDC();
      cms_tex_dphi->Draw(); 



      c_jet[l][mc_code]->cd(2);

      dPhi_pTdist_up[2][0][l][mc_code]->Draw();
      dPhi_pTdist_up[2][0][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      dPhi_pTdist_up[4][0][l][mc_code]->Draw("same");
      dPhi_pTdist_down[4][0][l][mc_code]->Draw("same");
      dPhi_pTdist_down[2][0][l][mc_code]->Draw("same");

      signal_dPhi_tot[2][0][0][l][mc_code]->Draw("same");
      signal_dPhi_tot[4][0][0][l][mc_code]->Draw("same");
 
      if(mc_code==5)    signal_dPhi_syst[2][0][0][l]->Draw("same e2");
      if(mc_code==5)    signal_dPhi_syst[4][0][0][l]->Draw("same e2");
 
   

      label_per = new TLatex(0.05,0.9,"P+H Cent. 50-100%");
      label_per->SetTextSize(0.07);
      label_per->SetLineColor(kWhite);
      label_per->SetNDC();
      label_per->Draw();
      l_phi->Draw();

      c_jet[l][mc_code]->cd(3);
   
      dPhi_pTdist_up[2][2][l][mc_code]->Draw();
      dPhi_pTdist_up[2][2][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      dPhi_pTdist_up[4][2][l][mc_code]->Draw("same");
      dPhi_pTdist_down[4][2][l][mc_code]->Draw("same");
      dPhi_pTdist_down[2][2][l][mc_code]->Draw("same");

    

      if(mc_code==5)   signal_dPhi_syst[2][0][3][l]->Draw("same e2");
      if(mc_code==5)  signal_dPhi_syst[4][0][3][l]->Draw("same e2");
 
      cout<<"1"<<endl;

      signal_dPhi_tot[2][0][2][l][mc_code]->Draw("same");
      signal_dPhi_tot[4][0][2][l][mc_code]->Draw("same");
 
  
      label_cent = new TLatex(0.05,0.9,"P+H Cent. 0-30%");
      label_cent->SetTextSize(0.07);
      label_cent->SetLineColor(kWhite);
      label_cent->SetNDC();
      label_cent->Draw();

      TLatex *orientation = new TLatex(.97,0.,"      <--Leading     SubLeading-->");
      orientation->SetTextSize(0.07);
      orientation->SetLineColor(kWhite);
      orientation->SetNDC();
      orientation->SetTextAngle(90.);
      orientation->Draw();
 
      l_phi->Draw();


      c_jet[l][mc_code]->cd(5);
      dPhi_pTdist_up[8][0][l][mc_code]->Draw();
      dPhi_pTdist_up[8][0][l][mc_code]->GetYaxis()->SetLabelSize(0.07); 
      dPhi_pTdist_up[8][0][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_pTdist_up[8][0][l][mc_code]->GetYaxis()->SetTitle("1/N_{evt} dp_{T}/d#Delta#phi (GeV/c)");
      dPhi_pTdist_up[8][0][l][mc_code]->GetYaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[8][0][l][mc_code]->GetYaxis()->SetTitleOffset(1.2);
      dPhi_pTdist_up[8][0][l][mc_code]->GetYaxis()->CenterTitle();
      dPhi_pTdist_down[8][0][l][mc_code]->Draw("same");


      if(mc_code==5)  signal_dPhi_syst[8][0][0][l]->Draw("same e2");
      signal_dPhi_tot[8][0][0][l][mc_code]->Draw("same");
    
      cout<<"2"<<endl;

      TLatex *label = new TLatex(0.05,0.9,"Subleading (P+H)-PYTHIA");
      label->SetTextSize(0.07);
      label->SetLineColor(kWhite);
      label->SetNDC();
      label->Draw();

      l_phi->Draw();


      c_jet[l][mc_code]->cd(6);

      dPhi_pTdist_up[8][2][l][mc_code]->Draw();
      dPhi_pTdist_up[8][2][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      l_phi->Draw();
      dPhi_pTdist_down[8][2][l][mc_code]->Draw("same");

      if(mc_code==5)     signal_dPhi_syst[8][0][3][l]->Draw("same e2");


      signal_dPhi_tot[8][0][2][l][mc_code]->Draw("same");

   
      c_jet[l][mc_code]->cd(8);
      dPhi_pTdist_up[10][0][l][mc_code]->Draw();
      dPhi_pTdist_up[10][0][l][mc_code]->GetYaxis()->SetLabelSize(0.07); 
      dPhi_pTdist_up[10][0][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_pTdist_up[10][0][l][mc_code]->GetYaxis()->SetTitle("1/N_{evt} dp_{T}/d#Delta#phi (GeV/c)");
      dPhi_pTdist_up[10][0][l][mc_code]->GetYaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[10][0][l][mc_code]->GetYaxis()->SetTitleOffset(1.2);
      dPhi_pTdist_up[10][0][l][mc_code]->GetYaxis()->CenterTitle();
      dPhi_pTdist_down[10][0][l][mc_code]->Draw("same");

      if(mc_code==5)    signal_dPhi_syst[10][0][0][l]->Draw("same e2");
      signal_dPhi_tot[10][0][0][l][mc_code]->Draw("same");
      cout<<"3"<<endl;
 
 
      label = new TLatex(0.05,0.9,"Leading (P+H)-PYTHIA (#times -1)");
      label->SetTextSize(0.07);
      label->SetLineColor(kWhite);
      label->SetNDC();
      label->Draw();

      l_phi->Draw();

      c_jet[l][mc_code]->cd(9);
      dPhi_pTdist_up[10][2][l][mc_code]->Draw();
      dPhi_pTdist_up[10][2][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      dPhi_pTdist_down[10][2][l][mc_code]->Draw("same");
      //  orientation->Draw();

      if(mc_code==5)  signal_dPhi_syst[10][0][3][l]->Draw("same e2");
      signal_dPhi_tot[10][0][2][l][mc_code]->Draw("same");
 
      cout<<"4"<<endl;

      l_phi->Draw();

      c_jet[l][mc_code]->cd(11);
      dPhi_pTdist_up[6][0][l][mc_code]->Draw();
      dPhi_pTdist_up[6][0][l][mc_code]->GetYaxis()->SetLabelSize(0.07); 
 
      dPhi_pTdist_up[6][0][l][mc_code]->GetYaxis()->SetTitle("1/N_{evt} dp_{T}/d#Delta#phi (GeV/c)");
      dPhi_pTdist_up[6][0][l][mc_code]->GetYaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[6][0][l][mc_code]->GetYaxis()->SetTitleOffset(1.2);
      dPhi_pTdist_up[6][0][l][mc_code]->GetYaxis()->CenterTitle();

      dPhi_pTdist_up[6][0][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_pTdist_up[6][0][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[6][0][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_pTdist_up[6][0][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
      dPhi_pTdist_down[6][0][l][mc_code]->Draw("same");

      if(mc_code==5)    signal_dPhi_syst[6][0][3][l]->Draw("same e2");
      signal_dPhi_tot[6][0][0][l][mc_code]->Draw("same");
 
 
      label = new TLatex(0.05,0.9,"Subleading-Leading");
      label->SetTextSize(0.07);
      label->SetLineColor(kWhite);
      label->SetNDC();
      label->Draw();
      TLatex *label2 = new TLatex(0.05,0.8,"((P+H)-PYTHIA)");
      label2->SetTextSize(0.07);
      label2->SetLineColor(kWhite);
      label2->SetNDC();
      label2->Draw();

      l_phi->Draw();

      c_jet[l][mc_code]->cd(12);
      dPhi_pTdist_up[6][2][l][mc_code]->Draw();
      dPhi_pTdist_up[6][2][l][mc_code]->GetYaxis()->SetLabelSize(0.);

      dPhi_pTdist_up[6][2][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_pTdist_up[6][2][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[6][2][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_pTdist_up[6][2][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
      dPhi_pTdist_down[6][2][l][mc_code]->Draw("same");

      if(mc_code==5)   signal_dPhi_syst[6][0][3][l]->Draw("same e2");
      signal_dPhi_tot[6][0][2][l][mc_code]->Draw("same");

   
 

      l_phi->Draw();
      //  orientation->Draw();


      c_longrange[l][mc_code]->cd(0);


      aj_tex = new TLatex(0.07,0.96,Ajlabel);
      aj_tex->SetTextSize(0.03);
      aj_tex->SetLineColor(kWhite);
      aj_tex->SetNDC();
      aj_tex->Draw();
      aj_tex->Draw();
   
      type_tex = new TLatex(0.3,0.96,"Long Range, |#Delta#eta|< 2.5");
      type_tex->SetTextSize(0.03);
      type_tex->SetLineColor(kWhite);
      type_tex->SetNDC();
      type_tex->Draw();
   
      TLatex  *luminosity_tex_pp = new TLatex(0.07,0.915,"PYTHIA 6.423 Tune Z2");
      luminosity_tex_pp->SetTextFont(43);
      luminosity_tex_pp->SetTextSizePixels(25);
      luminosity_tex_pp->SetLineColor(kWhite);
      luminosity_tex_pp->SetNDC();
      luminosity_tex_pp->Draw();
   
      TLatex  *luminosity_tex_PbPb = new TLatex(0.3,0.915,"HYDJET v1.8");
      luminosity_tex_PbPb->SetTextFont(43);
      luminosity_tex_PbPb->SetTextSizePixels(25);
      luminosity_tex_PbPb->SetLineColor(kWhite);
      luminosity_tex_PbPb->SetNDC();
      luminosity_tex_PbPb->Draw();
    

      jet_reco_tex = new TLatex(0.605,0.96,"anti-kT R = 0.3, |#eta_{jet}| < 1.6");
      jet_reco_tex->SetTextFont(43);
      jet_reco_tex->SetTextSizePixels(25);
      jet_reco_tex->SetLineColor(kWhite);
      jet_reco_tex->SetNDC();
      jet_reco_tex->Draw();

      jet_cut_tex = new TLatex(0.605,0.915,"120 < p_{T,1}< 300, p_{T,2}> 50 GeV/c, #Delta#phi_{1,2}> 5#pi/6");
      jet_cut_tex->SetTextFont(43);
      jet_cut_tex->SetTextSizePixels(25);
      jet_cut_tex->SetLineColor(kWhite);
      jet_cut_tex->SetNDC();
      jet_cut_tex->Draw();

      c_longrange[l][mc_code]->cd(4);
      legend->Draw();
     
      c_longrange[l][mc_code]->cd(1);
      dPhi_pTdist_up[1][3][l][mc_code]->Draw();
      label_pp->Draw();
 
      dPhi_pTdist_up[1][3][l][mc_code]->GetYaxis()->SetTitle("1/N_{evt} dp_{T}/d#Delta#phi (GeV/c)");
      dPhi_pTdist_up[1][3][l][mc_code]->GetYaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[1][3][l][mc_code]->GetYaxis()->CenterTitle();
      dPhi_pTdist_up[1][3][l][mc_code]->GetYaxis()->SetLabelSize(0.07); 


      dPhi_pTdist_up[1][3][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_pTdist_up[1][3][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[1][3][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_pTdist_up[1][3][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
      dPhi_pTdist_down[1][3][l][mc_code]->Draw("same");

      if(mc_code==5)   background_syst_rebin[1][0][3][l]->Draw("same e2");
      background_diff_tot[1][0][3][l][mc_code]->Draw("same");
  
      cout<<"5"<<endl;
      l_phi->Draw();

      cms_tex_dphi->Draw(); 
      prelim_tex_dphi->Draw(); 
      
      c_longrange[l][mc_code]->cd(2);

      dPhi_pTdist_up[0][0][l][mc_code]->Draw();
      dPhi_pTdist_up[0][0][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      dPhi_pTdist_down[0][0][l][mc_code]->Draw("same");

      if(mc_code==5)  background_syst_rebin[0][0][0][l]->Draw("same e2");
      background_diff_tot[0][0][0][l][mc_code]->Draw("same");
 
      label_per->Draw();
      l_phi->Draw();
      
      c_longrange[l][mc_code]->cd(3);
 
      dPhi_pTdist_up[0][2][l][mc_code]->Draw();
      dPhi_pTdist_up[0][2][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      dPhi_pTdist_down[0][2][l][mc_code]->Draw("same");
      label_cent->Draw();
      orientation->Draw();
   
      if(mc_code==5)    background_syst_rebin[0][0][3][l]->Draw("same e2");
      background_diff_tot[0][0][2][l][mc_code]->Draw("same");

      l_phi->Draw();
  
      c_longrange[l][mc_code]->cd(5);
      dPhi_pTdist_up[7][0][l][mc_code]->Draw();
  

      dPhi_pTdist_up[7][0][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_pTdist_up[7][0][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[7][0][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_pTdist_up[7][0][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
   
      dPhi_pTdist_up[7][0][l][mc_code]->GetYaxis()->SetTitle("1/N_{evt} dp_{T}/d#Delta#phi (GeV/c)");
      dPhi_pTdist_up[7][0][l][mc_code]->GetYaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[7][0][l][mc_code]->GetYaxis()->CenterTitle();
      dPhi_pTdist_up[7][0][l][mc_code]->GetYaxis()->SetLabelSize(0.07); 
      l_phi->Draw();
      dPhi_pTdist_down[7][0][l][mc_code]->Draw("same");

      if(mc_code==5)  background_syst_rebin[7][0][0][l]->Draw("same e2");
      background_diff_tot[7][0][0][l][mc_code]->Draw("same");

      cout<<"6"<<endl;

      TLatex *label_per2 = new TLatex(0.05,0.9,"(P+H)-PYTHIA 50-100% Cent.");
      label_per2->SetTextSize(0.07);
      label_per2->SetLineColor(kWhite);
      label_per2->SetNDC();
      label_per2->Draw();
  
     
      c_longrange[l][mc_code]->cd(6);
      dPhi_pTdist_up[7][2][l][mc_code]->Draw();
      dPhi_pTdist_up[7][2][l][mc_code]->GetYaxis()->SetLabelSize(0.);
   

      dPhi_pTdist_up[7][2][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_pTdist_up[7][2][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_pTdist_up[7][2][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_pTdist_up[7][2][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
      dPhi_pTdist_down[7][2][l][mc_code]->Draw("same");

     if(mc_code==5)  background_syst_rebin[7][0][3][l]->Draw("same e2");
      background_diff_tot[7][0][2][l][mc_code]->Draw("same");

      cout<<"7"<<endl;

      l_phi->Draw();
    

      c_nobkgsub[l][mc_code]->cd(0);


      aj_tex->Draw();
    
      type_tex = new TLatex(0.3,0.96,"Hemisphere Difference, |#Delta#eta|< 2.5");
      type_tex->SetTextSize(0.03);
      type_tex->SetLineColor(kWhite);
      type_tex->SetNDC();
      type_tex->Draw();

      luminosity_tex_pp->Draw();
      luminosity_tex_PbPb->Draw();
      jet_reco_tex->Draw();
      jet_cut_tex->Draw();
    



      c_nobkgsub[l][mc_code]->cd(4);





      legend->Draw();
      
    
      c_nobkgsub[l][mc_code]->cd(1);
 
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->Draw();
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->GetYaxis()->SetLabelSize(0.07); 
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->GetYaxis()->SetTitle("1/N_{evt} dp_{T}/d#Delta#phi (GeV/c)");
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->GetYaxis()->SetTitleSize(0.07);
      dPhi_Sub_Lead_NoBkg_up[1][3][l][mc_code]->GetYaxis()->CenterTitle();
      dPhi_Sub_Lead_NoBkg_down[1][3][l][mc_code]->Draw("same");   


      if(mc_code==5)     sub_lead_dPhi_syst[3][0][3][l]->Draw("same e2");
      sub_lead_dPhi_tot[1][0][3][l][mc_code]->Draw("same");

      label_pp->Draw();

      l_phi->Draw();

      cms_tex_dphi->Draw(); 
      prelim_tex_dphi->Draw(); 

      c_nobkgsub[l][mc_code]->cd(2);


      dPhi_Sub_Lead_NoBkg_up[0][0][l][mc_code]->Draw();
      dPhi_Sub_Lead_NoBkg_up[0][0][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      dPhi_Sub_Lead_NoBkg_down[0][0][l][mc_code]->Draw("same");

      if(mc_code==5)     sub_lead_dPhi_syst[2][0][0][l]->Draw("same e2");
      sub_lead_dPhi_tot[0][0][0][l][mc_code]->Draw("same");
      label_per->Draw();
      l_phi->Draw();


      c_nobkgsub[l][mc_code]->cd(3);
   
      dPhi_Sub_Lead_NoBkg_up[0][2][l][mc_code]->Draw();
      dPhi_Sub_Lead_NoBkg_up[0][2][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      dPhi_Sub_Lead_NoBkg_down[0][2][l][mc_code]->Draw("same");

      if(mc_code==5)     sub_lead_dPhi_syst[2][0][3][l]->Draw("same e2");
      sub_lead_dPhi_tot[0][0][2][l][mc_code]->Draw("same");
      label_cent->Draw();

      orientation->Draw();
 
      l_phi->Draw();


      c_nobkgsub[l][mc_code]->cd(5);
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->Draw();
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetYaxis()->SetLabelSize(0.07); 
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetYaxis()->CenterTitle();
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetYaxis()->SetTitle("1/N_{evt} dp_{T}/d#Delta#phi (GeV/c)");
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetYaxis()->SetTitleSize(0.07);
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetYaxis()->SetTitleOffset(1.2);

      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetXaxis()->SetTitleOffset(1.2);
      dPhi_Sub_Lead_NoBkg_up[6][0][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_Sub_Lead_NoBkg_down[6][0][l][mc_code]->Draw("same");

      if(mc_code==5) sub_lead_dPhi_syst[6][0][0][l]->Draw("same e2");
      sub_lead_dPhi_tot[6][0][0][l][mc_code]->Draw("same");

      label = new TLatex(0.05,0.9,"(P+H)-PYTHIA");
      label->SetTextSize(0.07);
      label->SetLineColor(kWhite);
      label->SetNDC();
      label->Draw();

      l_phi->Draw();


      c_nobkgsub[l][mc_code]->cd(6);

      dPhi_Sub_Lead_NoBkg_up[6][2][l][mc_code]->Draw();
      dPhi_Sub_Lead_NoBkg_up[6][2][l][mc_code]->GetYaxis()->SetLabelSize(0.);
      dPhi_Sub_Lead_NoBkg_up[6][2][l][mc_code]->GetXaxis()->SetLabelSize(0.07); 
      dPhi_Sub_Lead_NoBkg_up[6][2][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_Sub_Lead_NoBkg_up[6][2][l][mc_code]->GetXaxis()->SetTitle("#Delta#phi");
      dPhi_Sub_Lead_NoBkg_up[6][2][l][mc_code]->GetXaxis()->SetTitleSize(0.07);
      dPhi_Sub_Lead_NoBkg_up[6][2][l][mc_code]->GetXaxis()->SetTitleOffset(1.2);
      dPhi_Sub_Lead_NoBkg_up[6][2][l][mc_code]->GetXaxis()->CenterTitle();
      dPhi_Sub_Lead_NoBkg_down[6][2][l][mc_code]->Draw("same");

      if(mc_code==5)   sub_lead_dPhi_syst[6][0][3][l]->Draw("same e2");
      sub_lead_dPhi_tot[6][0][2][l][mc_code]->Draw("same");

      l_phi->Draw();
 

      if(l<2){


	if(use_highpT_bin){

	  c_longrange[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_LongRange_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+"_WithHighpT.pdf"));
	  c_longrange[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_LongRange_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+"_WithHighpT.png"));


	  c_nobkgsub[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_NoBkgSub_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+"_WithHighpT.png"));
	  c_nobkgsub[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_NoBkgSub_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+"_WithHighpT.pdf"));
 

	
	  c_jet[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_JetRelated_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+"_WithHighpT.png"));
	  c_jet[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_JetRelated_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+"_WithHighpT.pdf"));



	}else{
	  c_longrange[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_LongRange_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf"));
	  c_longrange[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_LongRange_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));


	  c_nobkgsub[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_NoBkgSub_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));
	  c_nobkgsub[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_NoBkgSub_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf"));
 

	
	  c_jet[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_JetRelated_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));
	  c_jet[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_JetRelated_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf"));

	}

      }else{

	if(use_highpT_bin){
	  c_jet[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_JetRelated_AjInclusive_WithHighpT.png"));
	  c_jet[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_JetRelated_AjInclusive_WithHighpT.pdf"));

	  c_longrange[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_LongRange_AjInclusive_WithHighpT.png"));
	  c_longrange[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_LongRange_AjInclusive_WithHighpT.pdf"));

	  c_nobkgsub[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_NoBkgSub_AjInclusive_WithHighpT.png"));
	  c_nobkgsub[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_NoBkgSub_AjInclusive_WithHighpT.pdf"));
	}else{

	  c_jet[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_JetRelated_AjInclusive.png"));
	  c_jet[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_JetRelated_AjInclusive.pdf"));

	  c_longrange[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_LongRange_AjInclusive.png"));
	  c_longrange[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_LongRange_AjInclusive.pdf"));

	  c_nobkgsub[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_NoBkgSub_AjInclusive.png"));
	  c_nobkgsub[l][mc_code]->SaveAs((TString)(data_mc_type_strs[mc_code]+"_MpT_NoBkgSub_AjInclusive.pdf"));



	}

      }
    
				    
      /*
      //-----------------------------
      //  Make Integrals
      //-----------------------------

      Integral[1][3][l][mc_code] = new TH1D(Form("GenGen_Integral_SideBand_pp_Aj%d",l),"",6,pt_bin_bounds));
      Integral[3][3][l][mc_code] = new TH1D(Form("GenGen_Integral_SubLeading_pp_Aj%d",l),"",6,pt_bin_bounds));
      Integral[5][3][l][mc_code] = new TH1D(Form("GenGen_Integral_Leading_pp_Aj%d",l),"",6,pt_bin_bounds);



      Integral_noerr_up[1][3][l][mc_code] = new TH1D(Form("GenGen_Integral_SideBand_pp_Aj%d_Noerr_Up",l),"",6,pt_bin_bounds);
      Integral_noerr_up[3][3][l][mc_code] = new TH1D(Form("GenGen_Integral_SubLeading_pp_Aj%d_Noerr_Up",l),"",6,pt_bin_bounds);
      Integral_noerr_up[5][3][l][mc_code] = new TH1D(Form("GenGen_Integral_Leading_pp_Aj%d_Noerr_Up",l),"",6,pt_bin_bounds);


      Integral_noerr_down[1][3][l][mc_code] = new TH1D(Form("GenGen_Integral_SideBand_pp_Aj%d_Noerr_Down",l),"",6,pt_bin_bounds);
      Integral_noerr_down[3][3][l][mc_code] = new TH1D(Form("GenGen_Integral_SubLeading_pp_Aj%d_Noerr_Down",l),"",6,pt_bin_bounds);
      Integral_noerr_down[5][3][l][mc_code] = new TH1D(Form("GenGen_Integral_Leading_pp_Aj%d_Noerr_Down",l),"",6,pt_bin_bounds);

      for(int i = 0; i<6; i++){

      integral = background_diff_rebin[1][i][0][l][mc_code]->IntegralAndError(1,background_diff_rebin[1][i][0][l][mc_code]->GetNbinsX(),int_err,"width");
      Integral[1][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral[1][3][l][mc_code]->SetBinError(i+1,int_err);
   
      if(integral>0){
      Integral_noerr_up[1][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_down[1][3][l][mc_code]->SetBinContent(i+1,0.);
      }else{
      Integral_noerr_down[1][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_up[1][3][l][mc_code]->SetBinContent(i+1,0.);
      }

      integral = signal_dPhi_rebin[3][i][0][l][mc_code]->IntegralAndError(1,signal_dPhi_rebin[3][i][0][l][mc_code]->GetNbinsX(),int_err,"width");
      Integral[3][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral[3][3][l][mc_code]->SetBinError(i+1,int_err);
 
      if(integral>0){
      Integral_noerr_up[3][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_down[3][3][l][mc_code]->SetBinContent(i+1,0.);
      }else{
      Integral_noerr_down[3][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_up[3][3][l][mc_code]->SetBinContent(i+1,0.);
      }


      integral = signal_dPhi_rebin[5][i][0][l][mc_code]->IntegralAndError(1,signal_dPhi_rebin[5][i][0][l][mc_code]->GetNbinsX(),int_err,"width");
      Integral[5][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral[5][3][l][mc_code]->SetBinError(i+1,int_err);

      if(integral>0){
      Integral_noerr_up[5][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_down[5][3][l][mc_code]->SetBinContent(i+1,0.);
      }else{
      Integral_noerr_down[5][3][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_up[5][3][l][mc_code]->SetBinContent(i+1,0.);
      }

      }

  
      for(int j = 0; j<4; j++){
 
      Integral[0][j][l][mc_code] = new TH1D(Form("GenGen_Integral_SideBand_PbPb_Cent%d_Aj%d",j,l),"",6,pt_bin_bounds);
      Integral[2][j][l][mc_code] = new TH1D(Form("GenGen_Integral_SubLeading_PbPb_Cent%d_Aj%d",j,l),"",6,pt_bin_bounds);
      Integral[4][j][l][mc_code] = new TH1D(Form("GenGen_Integral_Leading_PbPb_Cent%d_Aj%d",j,l),"",6,pt_bin_bounds);

      Integral_noerr_up[0][j][l][mc_code] = new TH1D(Form("GenGen_Integral_SideBand_PbPb_Cent%d_Aj%d_NoerrUp",j,l),"",6,pt_bin_bounds);
      Integral_noerr_up[2][j][l][mc_code] = new TH1D(Form("GenGen_Integral_SubLeading_PbPb_Cent%d_Aj%d_NoerrUp",j,l),"",6,pt_bin_bounds);
      Integral_noerr_up[4][j][l][mc_code] = new TH1D(Form("GenGen_Integral_Leading_PbPb_Cent%d_Aj%d_NoerrUp",j,l),"",6,pt_bin_bounds);


      Integral_noerr_down[0][j][l][mc_code] = new TH1D(Form("GenGen_Integral_SideBand_PbPb_Cent%d_Aj%d_NoerrDown",j,l),"",6,pt_bin_bounds);
      Integral_noerr_down[2][j][l][mc_code] = new TH1D(Form("GenGen_Integral_SubLeading_PbPb_Cent%d_Aj%d_NoerrDown",j,l),"",6,pt_bin_bounds);
      Integral_noerr_down[4][j][l][mc_code] = new TH1D(Form("GenGen_Integral_Leading_PbPb_Cent%d_Aj%d_NoerrDown",j,l),"",6,pt_bin_bounds);
  
   
      for(int i = 0; i<6; i++){

      integral = background_diff_rebin[0][i][j][l][mc_code]->IntegralAndError(1,background_diff_rebin[0][i][j][l][mc_code]->GetNbinsX(),int_err,"width");
      Integral[0][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral[0][j][l][mc_code]->SetBinError(i+1,int_err);

      if(integral>0){
      Integral_noerr_up[0][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_down[0][j][l][mc_code]->SetBinContent(i+1,0.);
      }else{
      Integral_noerr_down[0][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_up[0][j][l][mc_code]->SetBinContent(i+1,0.);
      }

      integral = signal_dPhi_rebin[2][i][j][l][mc_code]->IntegralAndError(1,signal_dPhi_rebin[2][i][j][l][mc_code]->GetNbinsX(),int_err,"width");
      Integral[2][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral[2][j][l][mc_code]->SetBinError(i+1,int_err);


      if(integral>0){
      Integral_noerr_up[2][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_down[2][j][l][mc_code]->SetBinContent(i+1,0.);
      }else{
      Integral_noerr_down[2][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_up[2][j][l][mc_code]->SetBinContent(i+1,0.);
      }

      integral= signal_dPhi_rebin[4][i][j][l][mc_code]->IntegralAndError(1,signal_dPhi_rebin[4][i][j][l][mc_code]->GetNbinsX(),int_err,"width");
      Integral[4][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral[4][j][l][mc_code]->SetBinError(i+1,int_err);


      if(integral>0){
      Integral_noerr_up[4][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_down[4][j][l][mc_code]->SetBinContent(i+1,0.);
      }else{
      Integral_noerr_down[4][j][l][mc_code]->SetBinContent(i+1,integral);
      Integral_noerr_up[4][j][l][mc_code]->SetBinContent(i+1,0.);
      }

      }
   
     
      Integral[6][j][l][mc_code] = (TH1D*)Integral[0][j][l][mc_code]->Clone(Form("GenGen_Integral_SideBand_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral[6][j][l][mc_code]->Add(Integral[1][3][l][mc_code],-1.);
 
      Integral[8][j][l][mc_code] = (TH1D*)Integral[2][j][l][mc_code]->Clone(Form("GenGen_Integral_SubLeading_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral[8][j][l][mc_code]->Add(Integral[3][3][l][mc_code],-1.);
 
      Integral[10][j][l][mc_code] = (TH1D*)Integral[4][j][l][mc_code]->Clone(Form("GenGen_Integral_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral[10][j][l][mc_code]->Add(Integral[5][3][l][mc_code], -1.);

      Integral[9][j][l][mc_code] = (TH1D*)Integral[8][j][l][mc_code]->Clone(Form("GenGen_Integral_Sub_minus_leading_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral[9][j][l][mc_code]->Add(Integral[10][j][l][mc_code]) ;
 


      Integral[0][j][l][mc_code]->SetMarkerSize(1.5);
      Integral[0][j][l][mc_code]->SetMarkerColor(kBlack);
      Integral[0][j][l][mc_code]->SetLineColor(kBlack);
      Integral[0][j][l][mc_code]->SetMarkerStyle(4);
      Integral[2][j][l][mc_code]->SetMarkerSize(1.5);
      Integral[2][j][l][mc_code]->SetMarkerColor(kBlack);
      Integral[2][j][l][mc_code]->SetLineColor(kBlack);
      Integral[2][j][l][mc_code]->SetMarkerStyle(4);
      Integral[4][j][l][mc_code]->SetMarkerSize(1.5);
      Integral[4][j][l][mc_code]->SetMarkerColor(kBlack);
      Integral[4][j][l][mc_code]->SetLineColor(kBlack);
      Integral[4][j][l][mc_code]->SetMarkerStyle(4);
      Integral[6][j][l][mc_code]->SetMarkerSize(1.5);
      Integral[6][j][l][mc_code]->SetMarkerColor(kBlack);
      Integral[6][j][l][mc_code]->SetLineColor(kBlack);
      Integral[6][j][l][mc_code]->SetMarkerStyle(4);
      Integral[8][j][l][mc_code]->SetMarkerSize(1.5);
      Integral[8][j][l][mc_code]->SetMarkerColor(kBlack);
      Integral[8][j][l][mc_code]->SetLineColor(kBlack);
      Integral[8][j][l][mc_code]->SetMarkerStyle(4);
      Integral[9][j][l][mc_code]->SetMarkerSize(1.5);
      Integral[9][j][l][mc_code]->SetMarkerColor(kBlack);
      Integral[9][j][l][mc_code]->SetLineColor(kBlack);
      Integral[9][j][l][mc_code]->SetMarkerStyle(4);
      Integral[10][j][l][mc_code]->SetMarkerSize(1.5);
      Integral[10][j][l][mc_code]->SetMarkerColor(kBlack); 
      Integral[10][j][l][mc_code]->SetLineColor(kBlack);
      Integral[10][j][l][mc_code]->SetMarkerStyle(4);
 



 
      Integral_noerr_up[6][j][l][mc_code] = (TH1D*)Integral[6][j][l][mc_code]->Clone(Form("GenGen_Integral_noerr_up_SideBand_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral_noerr_up[8][j][l][mc_code] = (TH1D*)Integral[8][j][l][mc_code]->Clone(Form("GenGen_Integral_noerr_up_SubLeading_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral_noerr_up[9][j][l][mc_code] = (TH1D*)Integral[9][j][l][mc_code]->Clone(Form("GenGen_Integral_noerr_up_SubLeading_minus_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral_noerr_up[10][j][l][mc_code] = (TH1D*)Integral[10][j][l][mc_code]->Clone(Form("GenGen_Integral_noerr_up_Leading_minus_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));

      Integral_noerr_down[6][j][l][mc_code] = (TH1D*)Integral[6][j][l][mc_code]->Clone(Form("GenGen_Integral_noerr_down_SideBand_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral_noerr_down[8][j][l][mc_code] = (TH1D*)Integral[8][j][l][mc_code]->Clone(Form("GenGen_Integral_noerr_down_SubLeading_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral_noerr_down[9][j][l][mc_code] = (TH1D*)Integral[9][j][l][mc_code]->Clone(Form("GenGen_Integral_noerr_down_SubLeading_minus_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));
      Integral_noerr_down[10][j][l][mc_code] =  (TH1D*)Integral[10][j][l][mc_code]->Clone(Form("GenGen_Integral_noerr_down_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));




      for(int k = 1; k< Integral_noerr_up[6][j][l][mc_code]->GetNbinsX()+1; k++){
      Integral_noerr_up[6][j][l][mc_code]->SetBinError(k,0.);
      Integral_noerr_up[8][j][l][mc_code]->SetBinError(k,0.);
      Integral_noerr_up[9][j][l][mc_code]->SetBinError(k,0.);
      Integral_noerr_up[10][j][l][mc_code]->SetBinError(k,0.);
   
      Integral_noerr_down[6][j][l][mc_code]->SetBinError(k,0.);
      Integral_noerr_down[8][j][l][mc_code]->SetBinError(k,0.);
      Integral_noerr_down[9][j][l][mc_code]->SetBinError(k,0.);
      Integral_noerr_down[10][j][l][mc_code]->SetBinError(k,0.);
     

      if(  Integral[6][j][l][mc_code]->GetBinContent(k)>0){
      Integral_noerr_down[6][j][l][mc_code]->SetBinContent(k,0.);
      }else{
      Integral_noerr_up[6][j][l][mc_code]->SetBinContent(k,0.);
      }


      if(  Integral[8][j][l][mc_code]->GetBinContent(k)>0){
      Integral_noerr_down[8][j][l][mc_code]->SetBinContent(k,0.);
      }else{
      Integral_noerr_up[8][j][l][mc_code]->SetBinContent(k,0.);
      }
     

      if(  Integral[9][j][l][mc_code]->GetBinContent(k)>0){
      Integral_noerr_down[9][j][l][mc_code]->SetBinContent(k,0.);
      }else{
      Integral_noerr_up[9][j][l][mc_code]->SetBinContent(k,0.);
      }


      if(  Integral[10][j][l][mc_code]->GetBinContent(k)>0){
      Integral_noerr_down[10][j][l][mc_code]->SetBinContent(k,0.);
      }else{
      Integral_noerr_up[10][j][l][mc_code]->SetBinContent(k,0.);
      }


      }



      }

  


      }//l

      c_longrange_int = new TCanvas("LongrangeIntegrals","",10,10,800,400);
      c_longrange_int->Divide(2,1,0,0);
   
      float int_max = 16.;
      float int_min = -4.;
  
      c_longrange_int->cd(1);
 
      Integral[1][3][0][4]->SetMarkerStyle(20); 
      Integral[1][3][0][4]->SetMarkerSize(2);
      Integral[1][3][0][4]->SetMarkerColor(kBlack);
      Integral[1][3][0][4]->SetLineColor(kBlack);
  
 
      Integral[1][3][1][4]->SetMarkerStyle(20); 
      Integral[1][3][1][4]->SetMarkerSize(2);
      Integral[1][3][1][4]->SetMarkerColor(kBlack);
      Integral[1][3][1][4]->SetLineColor(kBlack);

      Integral[0][0][0][4]->SetMarkerStyle(24); 
      Integral[0][0][0][4]->SetMarkerSize(2.5);
      Integral[0][0][0][4]->SetMarkerColor(kRed+2);
      Integral[0][0][0][4]->SetLineColor(kRed+2);



      Integral[0][0][1][4]->SetMarkerStyle(24); 
      Integral[0][0][1][4]->SetMarkerSize(2.5);
      Integral[0][0][1][4]->SetMarkerColor(kRed+2);
      Integral[0][0][1][4]->SetLineColor(kRed+2);



      Integral[0][3][0][4]->SetMarkerStyle(21); 
      Integral[0][3][0][4]->SetMarkerSize(2);
      Integral[0][3][0][4]->SetMarkerColor(kRed+2);
      Integral[0][3][0][4]->SetLineColor(kRed+2);
  
 
      Integral[0][3][1][4]->SetMarkerStyle(21); 
      Integral[0][3][1][4]->SetMarkerSize(2);
      Integral[0][3][1][4]->SetMarkerColor(kRed+2);
      Integral[0][3][1][4]->SetLineColor(kRed+2);
    
 
      Integral[0][0][0][4]->SetMaximum(int_max);
      Integral[0][0][0][4]->SetMinimum(int_min);
  
 
      Integral[0][0][0][4]->Draw();
      Integral[0][0][0][4]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/c)");
      Integral[0][0][0][4]->GetXaxis()->SetTitleSize(0.07);
      Integral[0][0][0][4]->GetXaxis()->CenterTitle();
      Integral[0][0][0][4]->GetXaxis()->SetLabelSize(0.07); 
      Integral[0][0][0][4]->GetXaxis()->SetRangeUser(0.501,7.99);
 

      Integral[0][0][0][4]->GetYaxis()->SetTitle("1/N_{evt} #Sigma p_{T} (GeV/c)");
      Integral[0][0][0][4]->GetYaxis()->SetTitleSize(0.07);
      Integral[0][0][0][4]->GetYaxis()->CenterTitle();
      Integral[0][0][0][4]->GetYaxis()->SetTitleOffset(0.8);
      Integral[0][0][0][4]->GetYaxis()->SetLabelSize(0.07); 
 


      Integral[0][0][0][4]->Draw();
      Integral[0][0][0][4]->Draw("same");
      Integral[0][3][0][4]->Draw("same");


      aj_tex = new TLatex(0.2,0.9,"A_{J}<0.22");
      aj_tex->SetTextSize(0.07);
      aj_tex->SetLineColor(kWhite);
      aj_tex->SetNDC();
      aj_tex->Draw(); 
 


  
      TLine *l_int = new TLine(0.5,0.,8.0,0.);
      l_int->SetLineStyle(2);
      l_int->Draw();


      TLegend *legend_int = new TLegend(0.2,0.6,0.97,0.85);

      Integral[0][3][0][4]->SetFillColor(kBlack);
      Integral[0][3][0][4]->SetFillStyle(3004);
 
      legend_int->AddEntry(Integral[0][0][0][4],"PbPb 50-100% 'Ridge' Asymmetry");
      legend_int->AddEntry(Integral[0][3][0][4],"PbPb 0-10% 'Ridge' Asymmetry","lpfe");
  
      legend_int->SetLineColor(kWhite);
      legend_int->SetTextSize(0.05);
      legend_int->Draw();

 

      c_longrange_int->cd(2);
 
      Integral[0][0][1][4]->SetMaximum(int_max);
      Integral[0][0][1][4]->SetMinimum(int_min);
 


      Integral[0][0][1][4]->Draw();
      Integral[0][0][1][4]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/c)");
      Integral[0][0][1][4]->GetXaxis()->SetTitleSize(0.07);
      Integral[0][0][1][4]->GetXaxis()->CenterTitle();
      Integral[0][0][1][4]->GetXaxis()->SetLabelSize(0.07); 
      Integral[0][0][1][4]->GetXaxis()->SetRangeUser(0.501,7.99);
      Integral[0][0][1][4]->GetYaxis()->SetTitleSize(0.0);
      Integral[0][0][1][4]->GetYaxis()->SetLabelSize(0.0);
  
      Integral[0][0][1][4]->Draw();
 
      Integral[0][0][1][4]->Draw("same");
      Integral[0][3][1][4]->Draw("same");
  
      TLatex *longrange_tex = new TLatex(0.05,0.8,"Long range 'ridge' asymmetry |#Delta#eta|<2.5");
      longrange_tex->SetTextSize(0.05);
      longrange_tex->SetLineColor(kWhite);
      longrange_tex->SetNDC();
      longrange_tex->Draw(); 
 


      aj_tex = new TLatex(0.05,0.9,"A_{J}>0.22");
      aj_tex->SetTextSize(0.07);
      aj_tex->SetLineColor(kWhite);
      aj_tex->SetNDC();
      aj_tex->Draw(); 
 

      TLatex *cms_tex = new TLatex(0.4,0.9,"CMS Preliminary");
      cms_tex->SetTextSize(0.07);
      cms_tex->SetLineColor(kWhite);
      cms_tex->SetNDC();
      cms_tex->Draw(); 





      l_int->Draw();
      c_longrange_int->SaveAs("GenGen_Integral_Longrange.png");
      c_longrange_int->SaveAs("GenGen_Integral_Longrange.pdf");

      //---------------------
      //Integral summary
      //---------------------


 


      c_all_int = new TCanvas("GenGen_IntegralSummary","",10,10,1000,800);
      c_all_int->Divide(2,2,0,0);
  
      c_all_int->cd(1);

 
    
      Integral_noerr_up[10][3][0][4]->SetMarkerSize(0);
      Integral_noerr_up[10][3][0][4]->SetFillColor(kYellow-9);
      Integral_noerr_up[10][3][1][4]->SetMarkerSize(0);
      Integral_noerr_up[10][3][1][4]->SetFillColor(kYellow-9);
  
      Integral_noerr_down[10][3][0][4]->SetMarkerSize(0);
      Integral_noerr_down[10][3][0][4]->SetFillColor(kYellow-9);
      Integral_noerr_down[10][3][1][4]->SetMarkerSize(0);
      Integral_noerr_down[10][3][1][4]->SetFillColor(kYellow-9);
 
      Integral_noerr_up[8][3][0][4]->SetMarkerSize(0);
      Integral_noerr_up[8][3][1][4]->SetMarkerSize(0);
      Integral_noerr_up[8][3][0][4]->SetFillColor(kGreen+3);
      Integral_noerr_up[8][3][1][4]->SetFillColor(kGreen+3);

      Integral_noerr_down[8][3][0][4]->SetMarkerSize(0);
      Integral_noerr_down[8][3][1][4]->SetMarkerSize(0);
      Integral_noerr_down[8][3][0][4]->SetFillColor(kGreen+3);
      Integral_noerr_down[8][3][1][4]->SetFillColor(kGreen+3);



 
      Integral_noerr_up[9][3][0][4]->SetMarkerSize(0);
      Integral_noerr_up[9][3][0][4]->SetFillColor(kBlue-9);
      Integral_noerr_up[9][3][1][4]->SetMarkerSize(0);
      Integral_noerr_up[9][3][1][4]->SetFillColor(kBlue-9);
 
      Integral_noerr_up[6][3][0][4]->SetMarkerSize(0);
      Integral_noerr_up[6][3][1][4]->SetMarkerSize(0);
      Integral_noerr_up[6][3][0][4]->SetFillColor(kRed+2);
      Integral_noerr_up[6][3][1][4]->SetFillColor(kRed+2);




 
      Integral_noerr_down[9][3][0][4]->SetMarkerSize(0);
      Integral_noerr_down[9][3][0][4]->SetFillColor(kBlue-9);
      Integral_noerr_down[9][3][1][4]->SetMarkerSize(0);
      Integral_noerr_down[9][3][1][4]->SetFillColor(kBlue-9);
 
      Integral_noerr_down[6][3][0][4]->SetMarkerSize(0);
      Integral_noerr_down[6][3][1][4]->SetMarkerSize(0);
      Integral_noerr_down[6][3][0][4]->SetFillColor(kRed+2);
      Integral_noerr_down[6][3][1][4]->SetFillColor(kRed+2);

 
 
      TString int_name = "int_jet_up_Aj0_Aj22";
 
      Integral_up[8][3][0][4] = new THStack(int_name, "");
      Integral_up[8][3][0][4]->Add(Integral_noerr_up[8][3][0][4]);
      Integral_up[8][3][0][4]->Add(Integral_noerr_up[10][3][0][4]);
 
  
      int_name.ReplaceAll("up","down");
      Integral_down[8][3][0][4] = new THStack(int_name, "");
      Integral_down[8][3][0][4]->Add(Integral_noerr_down[8][3][0][4]);
      Integral_down[8][3][0][4]->Add(Integral_noerr_down[10][3][0][4]);
   


      int_name = "int_jet_up_Aj22_Aj75";

      Integral_up[8][3][1][4] = new THStack(int_name, "");
      Integral_up[8][3][1][4]->Add(Integral_noerr_up[8][3][1][4]);
      Integral_up[8][3][1][4]->Add(Integral_noerr_up[10][3][1][4]);
  
  
      int_name.ReplaceAll("up","down");
      Integral_down[8][3][1][4] = new THStack(int_name, "");
      Integral_down[8][3][1][4]->Add(Integral_noerr_down[8][3][1][4]);
      Integral_down[8][3][1][4]->Add(Integral_noerr_down[10][3][1][4]);
  

  
      int_name = "int_diff_up_Aj0_Aj22";
 
      Integral_up[9][3][0][4] = new THStack(int_name, "");
      Integral_up[9][3][0][4]->Add(Integral_noerr_up[9][3][0][4]);
      Integral_up[9][3][0][4]->Add(Integral_noerr_up[6][3][0][4]);
 
  
      int_name.ReplaceAll("up","down");
      Integral_down[9][3][0][4] = new THStack(int_name, "");
      Integral_down[9][3][0][4]->Add(Integral_noerr_down[9][3][0][4]);
      Integral_down[9][3][0][4]->Add(Integral_noerr_down[6][3][0][4]);
   


      int_name = "int_diff_up_Aj22_Aj75";

      Integral_up[9][3][1][4] = new THStack(int_name, "");
      Integral_up[9][3][1][4]->Add(Integral_noerr_up[9][3][1][4]);
      Integral_up[9][3][1][4]->Add(Integral_noerr_up[6][3][1][4]);
  
  
      int_name.ReplaceAll("up","down");
      Integral_down[9][3][1][4] = new THStack(int_name, "");
      Integral_down[9][3][1][4]->Add(Integral_noerr_down[9][3][1][4]);
      Integral_down[9][3][1][4]->Add(Integral_noerr_down[6][3][1][4]);
  




      ///DRAWING//



      float summary_int_max = 14.;
      float summary_int_min = -14.;
 
      float summary_int_max2 = 12.;
      float summary_int_min2 = -13.;
 

      Integral_up[8][3][0][4]->Draw();
  

      Integral_up[8][3][0][4]->SetMaximum(summary_int_max);
      Integral_up[8][3][0][4]->SetMinimum(summary_int_min);
      Integral_up[8][3][0][4]->GetYaxis()->SetLabelSize(0.07);
      Integral_up[8][3][0][4]->GetYaxis()->SetTitleSize(0.065);
      Integral_up[8][3][0][4]->GetYaxis()->SetTitleOffset(2.5);
      Integral_up[8][3][0][4]->GetYaxis()->SetTitle("(1/N_{evt}#Sigmap_{T})_{PbPb}-(1/N_{evt}#Sigmap_{T})_{pp} (GeV/c)");
      Integral_up[8][3][0][4]->GetYaxis()->CenterTitle();
      Integral_up[8][3][0][4]->GetYaxis()->SetTitleOffset(0.8);

      Integral_up[8][3][0][4]->GetXaxis()->SetLabelSize(0.0); 
      Integral_up[8][3][0][4]->GetXaxis()->SetRangeUser(0.501,7.99);
  


      Integral_down[8][3][0][4]->Draw("same ");


      Integral[8][3][0][4]->Add(   Integral[10][3][0][4]);
      Integral[8][3][0][4]->Draw("same");
   

      l_int->Draw();



      TLegend *legend_int_lead = new TLegend(0.42,0.7,0.95,0.95);
      legend_int_lead->AddEntry(Integral_noerr_up[8][3][0][4],"SubLeading jet modification","f");
      legend_int_lead->AddEntry(Integral_noerr_up[10][3][0][4],"Leading jet modification","f");
      legend_int_lead->SetTextSize(0.055);
      legend_int_lead->SetLineColor(kWhite);
      legend_int_lead->Draw();
   

   
      aj_tex = new TLatex(0.2,0.9,"A_{J}<0.22");
      aj_tex->SetTextSize(0.07);
      aj_tex->SetLineColor(kWhite);
      aj_tex->SetNDC();
      aj_tex->Draw(); 
 

      TLine *l_int2 = new TLine(0.5,0.,8.0,0.);
      l_int2->SetLineStyle(2);
      l_int2->Draw();



      c_all_int->cd(2);

      Integral_up[8][3][1][4]->Draw();
  

      Integral_up[8][3][1][4]->SetMaximum(summary_int_max);
      Integral_up[8][3][1][4]->SetMinimum(summary_int_min);
      Integral_up[8][3][1][4]->GetYaxis()->SetLabelSize(0.0);
      Integral_up[8][3][1][4]->GetYaxis()->SetTitleSize(0.0);
 

      Integral_up[8][3][1][4]->GetXaxis()->SetLabelSize(0.0); 
      Integral_up[8][3][1][4]->GetXaxis()->SetRangeUser(0.501,7.99);
  
      Integral_down[8][3][1][4]->Draw("same ");

      Integral[8][3][1][4]->Add( Integral[10][3][1][4]);
      Integral[8][3][1][4]->Draw("same");
  

      l_int->Draw();

      aj_tex = new TLatex(0.05,0.9,"A_{J}>0.22");
      aj_tex->SetTextSize(0.07);
      aj_tex->SetLineColor(kWhite);
      aj_tex->SetNDC();
      aj_tex->Draw(); 
  
  
      cms_tex->Draw(); 




      TLatex *cent_tex = new TLatex(0.4,0.7,"0-10% Central");
      cent_tex->SetTextSize(0.07);
      cent_tex->SetLineColor(kWhite);
      cent_tex->SetNDC();
      cent_tex->Draw(); 
  

      TLatex *data_tex = new TLatex(0.4,0.8,"PbPb - pp");
      data_tex->SetTextSize(0.07);
      data_tex->SetLineColor(kWhite);
      data_tex->SetNDC();
      data_tex->Draw(); 


      l_int2->Draw();
      c_all_int->cd(3);

      cout<<"here"<<endl;
      Integral_up[9][3][0][4]->Draw();


      Integral_up[9][3][0][4]->SetMaximum(summary_int_max2);
      Integral_up[9][3][0][4]->SetMinimum(summary_int_min2);
      Integral_up[9][3][0][4]->GetYaxis()->SetLabelSize(0.06);
      Integral_up[9][3][0][4]->GetYaxis()->SetTitleSize(0.055);
      Integral_up[9][3][0][4]->GetYaxis()->SetTitleOffset(2.5);
      Integral_up[9][3][0][4]->GetYaxis()->SetTitle("(1/N_{evt}#Sigmap_{T})_{PbPb}-(1/N_{evt}#Sigmap_{T})_{pp} (GeV/c)");
      Integral_up[9][3][0][4]->GetYaxis()->CenterTitle();
      Integral_up[9][3][0][4]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/c)");
      Integral_up[9][3][0][4]->GetXaxis()->SetTitleSize(0.07);
      Integral_up[9][3][0][4]->GetXaxis()->CenterTitle();
      Integral_up[9][3][0][4]->GetXaxis()->SetLabelSize(0.07); 
      Integral_up[9][3][0][4]->GetYaxis()->SetTitleOffset(0.8);
      Integral_up[9][3][0][4]->GetXaxis()->SetRangeUser(.501,7.99);
   
      Integral_down[9][3][0][4]->Draw("same");

      Integral[9][3][0][4]->Add(Integral[6][3][0][4]);
      Integral[9][3][0][4]->Draw("same");
     


      l_int->Draw();

      TLegend *legend_int_diff = new TLegend(0.42,0.7,0.95,0.95);
      legend_int_diff->AddEntry(Integral_noerr_up[9][3][0][4],"Dijet modification","f");
      legend_int_diff->AddEntry(Integral_noerr_up[6][3][0][4],"Long range difference","f");
      legend_int_diff->SetTextSize(0.055);
      legend_int_diff->SetLineColor(kWhite);
      legend_int_diff->Draw();

      l_int2->Draw();

      c_all_int->cd(4);
 
      Integral_up[9][3][1][4]->Draw();


      Integral_up[9][3][1][4]->SetMaximum(summary_int_max2);
      Integral_up[9][3][1][4]->SetMinimum(summary_int_min2);
      Integral_up[9][3][1][4]->GetYaxis()->SetLabelSize(0.0);
      Integral_up[9][3][1][4]->GetYaxis()->SetTitleSize(0.0);
      Integral_up[9][3][1][4]->GetXaxis()->SetTitle("p_{T}^{assoc} (GeV/c)");
      Integral_up[9][3][1][4]->GetXaxis()->SetTitleSize(0.07);
      Integral_up[9][3][1][4]->GetXaxis()->CenterTitle();
      Integral_up[9][3][1][4]->GetXaxis()->SetLabelSize(0.07); 
      Integral_up[9][3][1][4]->GetYaxis()->SetTitleOffset(0.8);
      Integral_up[9][3][1][4]->GetXaxis()->SetRangeUser(.501,7.99);
   
      Integral_down[9][3][1][4]->Draw("same");

      Integral[9][3][1][4]->Add( Integral[6][3][1][4]);
      Integral[9][3][1][4]->Draw("same");
     



      l_int2->Draw();

      c_all_int->SaveAs("GenGen_Integral_PbPb_pp_Summary.png");
      c_all_int->SaveAs("GenGen_Integral_PbPb_pp_Summary.pdf");
      */
    }
  }

  return 0;

}  //and we're done.
 
