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
#include "TGraphErrors.h"


#include <iostream>
#include <vector>
#include <fstream>

#include "../JetTrack2015_functions.h"


Int_t jff_residual_aj(bool is_recogen = kFALSE, bool is_number = kFALSE){


  gROOT->ForceStyle();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);



  const int nCBins= 4;
  const int nPtBins=1;
  const int nTrkPtBins=6;
  const int nAjBins = 3;


 
  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};

  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};

  int data_mc_type_code = -999;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};

 
  float TrkPtBins[nTrkPtBins+1] = {0.5,1, 2, 3, 4, 8, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};


  float AjBins[nAjBins+1] = {0,0.22,0.75};
  TString AjBin_strs[nAjBins+2] = {"Aj0","Aj22","Aj75"};
 
 

  float x, offset, value;
  TF1 *do_offset = new TF1("do_offset","-1.*[0]+x-x",-3.,3.);



  Double_t xAxis[6] = {-100,-50,-30,-10,0}; 
  TH1D* int_cent[12][6][3];
  TH1D* blank[6][3];
  TH1D* blank2[6][3];

   
  TFile *fin[12];
  
  TH2D *result[12][6][4][3];
  TH2D *result2[12][6][4][3];


  TH2D* background[12][6][4][3];
  TH1D* background_left[12][6][4][3];
  TH1D* background_right[12][6][4][3];
  TH1D* background_proj[12][6][4][3];


  TH1D *phi_proj[12][6][4][3];
  TH1D *phi_proj_rebin[12][6][4][3];
  TH1D *phi_proj_rebin2[12][6][4][3];
  TH1D *eta_proj[12][6][4][3];
  TH1D *eta_proj_rebin[12][6][4][3];
  TH1D *eta_proj_rebin2[12][6][4][3];

  TH1D *eta_proj_ref[12][6][4][3];
  TH1D *phi_proj_ref[12][6][4][3];


  TH1D *jff_residual_eta[12][6][4][3];
  TH1D *jff_residual_phi[12][6][4][3];


  TH1D *jff_percent_eta[12][6][4][3];
  TH1D *jff_percent_phi[12][6][4][3];

 


  TCanvas *corr_canvas_eta[12][6][3];
  TCanvas *corr_canvas_phi[12][6][3];



  vector<float> pTbin_centers;
  pTbin_centers.push_back(0.75);
  pTbin_centers.push_back(1.5);
  pTbin_centers.push_back(2.5);
  pTbin_centers.push_back(3.5);
  pTbin_centers.push_back(6.0);
  pTbin_centers.push_back(12.0);
 
  vector<float> pTbin_errors;
  pTbin_errors.push_back(.25);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(2.);
  pTbin_errors.push_back(4.);

  vector<float> Closure_integral_eta0;
  vector<float> Closure_integral_phi0;
  vector<float> Closure_integral_eta1;
  vector<float> Closure_integral_phi1;
  vector<float> Closure_integral_eta2;
  vector<float> Closure_integral_phi2;
  vector<float> Closure_integral_eta3;
  vector<float> Closure_integral_phi3;
 
 
  TGraphErrors *Closure_integral_eta_pT[12][4][3];
  TGraphErrors *Closure_integral_phi_pT[12][4][3];


  TGraphErrors *Closure_integral_eta_pT2[12][4][3];
 

  vector<float> closure_integral_values, closure_integral_errors;

  TH1D *Closure_integral_eta_cent[12][5][3];
  TH1D *Closure_integral_eta_cent2[12][5][3];

  TLine *lineCent, *linePt;


  TCanvas *cintegral_eta_pT[12][3];
  TCanvas *cintegral_phi_pT[12][3];
   
  TCanvas *cintegral_eta_cent[12][3];
  TCanvas *cintegral_phi_cent[12][3];
  
  TString in_name, plotname, outname, funcname, centlabel, datalabel,pTlabel,Ajlabel;

  Double_t check_ymax, check_ymin, residual_ymin, residual_ymax, dx_eta, dx_phi, bc, err, evalpt, temp1, err1, r, deta, dphi;

  TLegend *lcheck, *leta, *lHminusP;
  
  TLine *linePhi, *lineEta;
 
  TLegend *l40,*l41,*l42;

  int llimiteta1, rlimiteta1,llimiteta2, rlimiteta2; 

  /////////////////////////

  etalim = 1.;
  philim = 1.;

  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  int gstart = 2; 
  int gend = 12;


  TFile *fout;
  if(is_number)  fout = new TFile("All_JFFResiduals_NoPtWeight.root","RECREATE");
  else if(is_recogen) fout = new TFile("All_JFFResiduals_RecoGen.root","RECREATE");
  else fout = new TFile("All_JFFResiduals_RecoReco.root","RECREATE");

  for(int g=gstart; g<gend; g++){

    if(g==6||g==7)continue; 
       
    switch(g){
    case 2:
      fin[g] = new TFile("../me_correct_mc/Pythia_RecoJet_RecoTrack_Dijet_Correlations.root","READ");
      if(is_recogen)   fin[g] = new TFile("../me_correct_mc/Pythia_RecoJet_GenTrack_Dijet_Correlations.root","READ");
      break;
    case 3:
      fin[g] = new TFile("../me_correct_mc/Pythia_GenJet_GenTrack_Dijet_Correlations.root","READ");
      break;
    case 4:
      break;
    case 5:
      break;
    case 8:
      fin[g] = new TFile("../me_correct_mc/HydJet_RecoJet_RecoTrack_Dijet_Correlations.root","READ");
      if(is_recogen) fin[g] = new TFile("../me_correct_mc/HydJet_RecoJet_GenTrack_Sube0_Dijet_Correlations.root","READ");
      break;
    case 9:
      fin[g] = new TFile("../me_correct_mc/HydJet_GenJet_GenTrack_Dijet_Correlations.root","READ");
       if(is_recogen) fin[g] = new TFile("../me_correct_mc/HydJet_GenJet_GenTrack_Sube0_Dijet_Correlations.root","READ");
      break;
    case 10:
      break;
    case 11:
      break;
    }
    
  
    
    //----------------------------------------------------
    //  Start of main i & j loops 
    //-----------------------------------------------------

    for(int l = 0; l<nAjBins; l++){
      for(int i=0; i<nTrkPtBins; i++){
	
	if(g>5){

	  corr_canvas_eta[g][i][l] = new TCanvas(Form("CorrCanvasEta%d%d%d",g,i,l)," ",10,10,1500,800);
	  corr_canvas_eta[g][i][l]->Divide(4,2,0.,0.);
  
	  corr_canvas_phi[g][i][l] = new TCanvas(Form("CorrCanvasPhi%d%d%d",g,i,l)," ",10,10,1500,800);
	  corr_canvas_phi[g][i][l]->Divide(4,2,0.,0.);
	}
    

	for(int j=0; j<4; j++){
	
	  if(g<6&&j<3)continue;

	  TString in_name = make_name("Raw_Yield_pTweighted_",g,i,j,l,centlabel,pTlabel,Ajlabel);
	  if(is_number)  in_name = make_name("Raw_Yield_",g,i,j,l,centlabel,pTlabel,Ajlabel);

	  in_name.ReplaceAll("Pt100_Pt300_","");

	  /*	  in_name.ReplaceAll("PbPb","Pythia");
	  in_name.ReplaceAll("pp","Pythia");
	  //	  in_name.ReplaceAll("Hydjet","");
	  if(g>5) in_name.ReplaceAll("Pythia","Hydjet");
	  */

	  
	  in_name.ReplaceAll("PbPb_","");
	  in_name.ReplaceAll("pp_","");
	  in_name.ReplaceAll("Pythia_","");
	  in_name.ReplaceAll("Hydjet_","");

	  cout<<g<<" "<<in_name<<endl;

	  if(g<4||(g>5&&g<10)){
	    result[g][i][j][l] = (TH2D*)fin[g]->Get(in_name)->Clone(in_name);
	    }else{
	      result[g][i][j][l] = (TH2D*)fin[g-2]->Get(in_name)->Clone(in_name);
	    }

	  cout<<"got histo"<<endl;

	  for(int k = 1; k< result[g][i][j][l]->GetNbinsX(); k++){
	    for(int m = 1; m< result[g][i][j][l]->GetNbinsY(); m++){

	      deta = result[g][i][j][l]->GetXaxis()->GetBinCenter(k);
	      dphi = result[g][i][j][l]->GetYaxis()->GetBinCenter(m);
		
	      r = TMath::Sqrt(deta*deta+dphi*dphi);
	  
	      if(r>0.3){
		result[g][i][j][l]->SetBinContent(k,m,0.);
		result[g][i][j][l]->SetBinError(k,m,0.);

	      }
	    }
	  }

	  if(g%2==0){
	    llimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-2.+0.0001);
	    rlimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-1.0-0.0001);
	  	    
	    llimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(1.+0.0001);
	    rlimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(2.-0.0001);

	  }else{
	    llimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-2.+0.0001);
	    rlimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-1.0-0.0001);
	    
	    llimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(1.+0.0001);
	    rlimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(2.-0.0001);
	  }
	    
	  //-------------------------------
	  //dEta projection
	  //------------------------

	  TString eta_proj_name= in_name;
	  eta_proj_name.ReplaceAll("Raw_Yield","Eta_Proj");
	  eta_proj_name.ReplaceAll("Yield","Eta_Proj");
	    
	  llimiteta = result[g][i][j][l]->GetXaxis()->FindBin(-etalim+.001);
	  rlimiteta = result[g][i][j][l]->GetXaxis()->FindBin(etalim-.001);

	  llimitphi = result[g][i][j][l]->GetYaxis()->FindBin(-philim+.001);
	  rlimitphi = result[g][i][j][l]->GetYaxis()->FindBin(philim-.001);
	    

	  eta_proj[g][i][j][l] = result[g][i][j][l]->ProjectionX(eta_proj_name,llimitphi,rlimitphi);
	  dx_eta = eta_proj[g][i][j][l]->GetBinWidth(1);
	  eta_proj[g][i][j][l]->Scale(1/dx_eta);

	  for(int k = 1; k<eta_proj[g][i][j][l]->GetNbinsX()/2+1; k++){
	    bc =   (eta_proj[g][i][j][l]->GetBinContent(k)+eta_proj[g][i][j][l]->GetBinContent(eta_proj[g][i][j][l]->GetNbinsX()+1-k))/2.;
	    err =   TMath::Sqrt(eta_proj[g][i][j][l]->GetBinError(k)*eta_proj[g][i][j][l]->GetBinError(k)+eta_proj[g][i][j][l]->GetBinError(eta_proj[g][i][j][l]->GetNbinsX()+1-k*eta_proj[g][i][j][l]->GetBinError(eta_proj[g][i][j][l]->GetNbinsX()+1-k)))/2.;
	    eta_proj[g][i][j][l]->SetBinContent(k,bc);
	    eta_proj[g][i][j][l]->SetBinContent(eta_proj[g][i][j][l]->GetNbinsX()+1-k,bc);

	    eta_proj[g][i][j][l]->SetBinError(k,err);
	    eta_proj[g][i][j][l]->SetBinError(eta_proj[g][i][j][l]->GetNbinsX()+1-k,err);
	  }	    



	  TString eta_proj_name_rebin = eta_proj_name;
	  eta_proj_name_rebin.ReplaceAll("Eta_Proj","Eta_Proj_Rebin");
	 
	  eta_proj_rebin[g][i][j][l] = (TH1D*)Rebin_dEta(eta_proj[g][i][j][l]);
	  eta_proj_rebin[g][i][j][l]->SetName(eta_proj_name_rebin);


	  //-------------------------------
	  //dPhi projection
	  //------------------------

	  TString phi_proj_name= in_name;
	  phi_proj_name.ReplaceAll("Raw_Yield","Phi_Proj");
	  phi_proj_name.ReplaceAll("Yield","Phi_Proj");

	  phi_proj[g][i][j][l] = result[g][i][j][l]->ProjectionY(phi_proj_name,llimiteta,rlimiteta);
	  dx_phi = phi_proj[g][i][j][l]->GetBinWidth(1);
	  phi_proj[g][i][j][l]->Scale(1/dx_phi);



	  int nbins_half = phi_proj[g][i][j][l]->GetNbinsX()/2;

	  for(int k = 1; k<nbins_half/2+1; k++){
	    bc =   (phi_proj[g][i][j][l]->GetBinContent(k)+phi_proj[g][i][j][l]->GetBinContent(nbins_half+1-k))/2.;
	    err =   TMath::Sqrt(phi_proj[g][i][j][l]->GetBinError(k)*phi_proj[g][i][j][l]->GetBinError(k)+phi_proj[g][i][j][l]->GetBinError(nbins_half+1-k)*phi_proj[g][i][j][l]->GetBinError(nbins_half+1-k))/2.;
	    phi_proj[g][i][j][l]->SetBinContent(k,bc);
	    phi_proj[g][i][j][l]->SetBinContent(nbins_half+1-k,bc);

	    phi_proj[g][i][j][l]->SetBinError(k,err);
	    phi_proj[g][i][j][l]->SetBinError(nbins_half+1-k,err);
	  }	    

	  TString phi_proj_name_rebin = phi_proj_name;
	  phi_proj_name_rebin.ReplaceAll("Phi_Proj","Phi_Proj_Rebin");
  

	  phi_proj_rebin[g][i][j][l] = (TH1D*)Rebin_dPhi(phi_proj[g][i][j][l]);
	  phi_proj_rebin[g][i][j][l]->SetName(phi_proj_name_rebin);
	
	  switch(i){
	  case 0: 
	    check_ymax = 5.;
	    check_ymin = -1.; 
	    residual_ymax = 1.;
	    residual_ymin = -1.; 

	    break;
	  case 1: 
	    check_ymax = 11.; 
	    check_ymin = -2.; 
	    residual_ymax = 2.;
	    residual_ymin = -2.; 

	    break;
	  case 2: 
	    check_ymax = 17.; 
	    check_ymin = -3.; 
	    residual_ymax = 3.;
	    residual_ymin = -3.; 

	    break;
	  case 3: 
	    check_ymax = 20.; 
	    check_ymin = -5.; 
	    residual_ymax = 5;
	    residual_ymin = -5; 

	   
	    break;
	  case 4: 
	    check_ymax = 100.;
	    check_ymin = -15.; 
	    residual_ymax = 15.;
	    residual_ymin = -15.; 

	    break;
	  case 5: 
	    check_ymax = 500.;
	    check_ymin = -50.; 
	    residual_ymax = 50;
	    residual_ymin = -50; 

	    break;
	  }


	  eta_proj_rebin[g][i][j][l]->SetLineColor(kBlack);
	  eta_proj_rebin[g][i][j][l]->SetMarkerColor(kBlack);
	  eta_proj_rebin[g][i][j][l]->SetMarkerStyle(20);
	  eta_proj_rebin[g][i][j][l]->SetMarkerSize(1);
	  eta_proj_rebin[g][i][j][l]->SetMinimum(check_ymin);
	  eta_proj_rebin[g][i][j][l]->SetMaximum(check_ymax);

	  phi_proj_rebin[g][i][j][l]->SetLineColor(kBlack);
	  phi_proj_rebin[g][i][j][l]->SetMarkerColor(kBlack);
	  phi_proj_rebin[g][i][j][l]->SetMarkerStyle(20);
	  phi_proj_rebin[g][i][j][l]->SetMarkerSize(1);
	  phi_proj_rebin[g][i][j][l]->SetMinimum(check_ymin);
	  phi_proj_rebin[g][i][j][l]->SetMaximum(check_ymax);

	  if(g<6){
	    phi_proj_rebin[g][i][j][l]->SetLineColor(kRed);
	    phi_proj_rebin[g][i][j][l]->SetMarkerColor(kRed);
	    eta_proj_rebin[g][i][j][l]->SetLineColor(kRed);
	    eta_proj_rebin[g][i][j][l]->SetMarkerColor(kRed);
	  }


	  //-------------------

	  //   Saving & Plotting!

	  //--------------------
	
	  if(g%2==0)continue;

	  TString residual_name_eta =  make_name("JFF_Residual_Eta_",g,i,j,l,centlabel,pTlabel,Ajlabel);
	  if(g>6)residual_name_eta.ReplaceAll("Pythia","Hydjet");
	  else residual_name_eta.ReplaceAll("pp","Pythia");

	  jff_residual_eta[g][i][j][l] = (TH1D*)eta_proj_rebin[g-1][i][j][l]->Clone(residual_name_eta);
	  jff_residual_eta[g][i][j][l]->Add(eta_proj_rebin[g][i][j][l],-1.);
	  jff_residual_eta[g][i][j][l]->SetMinimum(residual_ymin);
	  jff_residual_eta[g][i][j][l]->SetMaximum(residual_ymax);
	 
	  fout->cd();

	  jff_residual_eta[g][i][j][l]->Write();

	  TString residual_name_phi =  make_name("JFF_Residual_Phi_",g,i,j,l,centlabel,pTlabel,Ajlabel);
	  if(g>6)residual_name_phi.ReplaceAll("Pythia","Hydjet");
	  else residual_name_phi.ReplaceAll("pp","Pythia");

	  jff_residual_phi[g][i][j][l] = (TH1D*)phi_proj_rebin[g-1][i][j][l]->Clone(residual_name_phi);
	  jff_residual_phi[g][i][j][l]->Add(phi_proj_rebin[g][i][j][l],-1.);
	  jff_residual_phi[g][i][j][l]->SetMinimum(residual_ymin);
	  jff_residual_phi[g][i][j][l]->SetMaximum(residual_ymax);

	  jff_residual_phi[g][i][j][l]->Write();



	  /*
	  TString percent_name_eta =  make_name("JFF_Percent_Eta_",g,i,j,l,centlabel,pTlabel,Ajlabel);
	  if(g>6)percent_name_eta.ReplaceAll("Pythia","Hydjet");
	  else percent_name_eta.ReplaceAll("pp","Pythia");

	  jff_percent_eta[g][i][j][l] = (TH1D*)jff_residual_eta[g][i][j][l]->Clone(percent_name_eta);
	  jff_percent_eta[g][i][j][l]->Divide(eta_proj_rebin[g][i][j][l]);
	  jff_percent_eta[g][i][j][l]->SetMinimum(residual_ymin);
	  jff_percent_eta[g][i][j][l]->SetMaximum(residual_ymax);


	  TString percent_name_phi =  make_name("JFF_Percent_Phi_",g,i,j,l,centlabel,pTlabel,Ajlabel);
	  if(g>6)percent_name_phi.ReplaceAll("Pythia","Hydjet");
	  else percent_name_phi.ReplaceAll("pp","Pythia");

	  jff_percent_phi[g][i][j][l] = (TH1D*)jff_residual_phi[g][i][j][l]->Clone(percent_name_phi);
	  jff_percent_phi[g][i][j][l]->Divide(phi_proj_rebin[g][i][j][l]);
	  jff_percent_phi[g][i][j][l]->SetMinimum(residual_ymin);
	  jff_percent_phi[g][i][j][l]->SetMaximum(residual_ymax);

	  */ 
	  if((i==0&&j==0)||(g%2!=0&&i==0)){
	    Closure_integral_eta0.clear();
	    Closure_integral_phi0.clear();
	    Closure_integral_eta1.clear();
	    Closure_integral_phi1.clear();
	    Closure_integral_eta2.clear();
	    Closure_integral_phi2.clear();
	    Closure_integral_eta3.clear();
	    Closure_integral_phi3.clear();
	  }

	  llimiteta = jff_residual_eta[g][i][j][l]->GetXaxis()->FindBin(-1.0+.0001);
	  rlimiteta = jff_residual_eta[g][i][j][l]->GetXaxis()->FindBin(1.0-.0001);
	     
	  double Yield_eta = jff_residual_eta[g][i][j][l]->Integral(llimiteta,rlimiteta,"width")/ abs(eta_proj_rebin[g][i][j][l]->Integral(llimiteta,rlimiteta,"width"));	      



	  llimitphi = jff_residual_phi[g][i][j][l]->GetXaxis()->FindBin(-1.0+.0001);
	  rlimitphi = jff_residual_phi[g][i][j][l]->GetXaxis()->FindBin(1.0-.0001);

	  double Yield_phi = jff_residual_phi[g][i][j][l]->Integral(llimitphi,rlimitphi,"width")/abs(phi_proj_rebin[g][i][j][l]->Integral(llimitphi,rlimitphi,"width"));	      

	  switch(j){
	  case 0:
	    Closure_integral_eta0.push_back(Yield_eta);
	    Closure_integral_phi0.push_back(Yield_eta);
	    break;
	  case 1:
	    Closure_integral_eta1.push_back(Yield_eta);
	    Closure_integral_phi1.push_back(Yield_eta);
	    break;
	  case 2:
	    Closure_integral_eta2.push_back(Yield_eta);
	    Closure_integral_phi2.push_back(Yield_eta);
	    break;
	  case 3:
	    Closure_integral_eta3.push_back(Yield_eta);
	    Closure_integral_phi3.push_back(Yield_eta);
	    break;
	  }
      

	  if(g>5){
	
	    corr_canvas_eta[g][i][l]->cd(j+1);


	    eta_proj_rebin[g][i][j][l]->SetMarkerStyle(4);
	    if(j==0) eta_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.06);
	    else eta_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    eta_proj_rebin[g][i][j][l]->Draw();


	  
	    eta_proj_rebin[g-1][i][j][l]->Draw("same");

	    eta_proj_rebin[g-7][i][3][l]->Draw("same");
	    eta_proj_rebin[g-6][i][3][l]->SetMarkerStyle(4);
	    eta_proj_rebin[g-6][i][3][l]->Draw("same");

	    drawlabels(g,i,j);

	    TLine *l_eta = new TLine(-1.5,0.,1.5,0.);
	    l_eta->SetLineStyle(2);
	    l_eta->Draw();
	  
	    corr_canvas_eta[g][i][l]->cd(j+5);

	    if(j==0) jff_residual_eta[g][i][j][l]->GetYaxis()->SetLabelSize(0.06);
	    else jff_residual_eta[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    jff_residual_eta[g][i][j][l]->GetXaxis()->SetLabelSize(0.06);
	    jff_residual_eta[g][i][j][l]->GetXaxis()->SetTitleSize(0.06);
	    jff_residual_eta[g][i][j][l]->GetXaxis()->SetTitle("#Delta#eta");
	    jff_residual_eta[g][i][j][l]->GetXaxis()->CenterTitle();
	  
	    jff_residual_eta[g][i][j][l]->Draw();
	    jff_residual_eta[g-6][i][3][l]->Draw("same");
	    l_eta->Draw();
	  
	    TLegend *legend = new TLegend(0.2,0.75,0.9,0.95);
	    legend->AddEntry( eta_proj_rebin[g-1][i][j][l],"P+H Sube0 RecoGen");
	    legend->AddEntry( eta_proj_rebin[g][i][j][l],"P+H Sube0 GenGen");
	    legend->AddEntry( eta_proj_rebin[g-7][i][3][l],"Pythia RecoGen");
	    legend->AddEntry( eta_proj_rebin[g-6][i][3][l],"Pythia GenGen");

	    legend->SetTextSize(0.05);
	    legend->SetLineColor(kWhite);

	    if(j==0)legend->Draw();

					


	    corr_canvas_phi[g][i][l]->cd(j+1);

	    phi_proj_rebin[g][i][j][l]->SetMarkerStyle(4);
	    if(j==0) phi_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.06);
	    else phi_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    phi_proj_rebin[g][i][j][l]->Draw();


	  
	    phi_proj_rebin[g-1][i][j][l]->Draw("same");

	    phi_proj_rebin[g-7][i][3][l]->Draw("same");
	    phi_proj_rebin[g-6][i][3][l]->SetMarkerStyle(4);
	    phi_proj_rebin[g-6][i][3][l]->Draw("same");

	    drawlabels(g,i,j);

	    TLine *l_phi = new TLine(-1.5,0.,1.5,0.);
	    l_phi->SetLineStyle(2);
	    l_phi->Draw();
	  
	    corr_canvas_phi[g][i][l]->cd(j+5);

	    if(j==0) jff_residual_phi[g][i][j][l]->GetYaxis()->SetLabelSize(0.06);
	    else jff_residual_phi[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    jff_residual_phi[g][i][j][l]->GetXaxis()->SetLabelSize(0.06);
	    jff_residual_phi[g][i][j][l]->GetXaxis()->SetTitleSize(0.06);
	    jff_residual_phi[g][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	    jff_residual_phi[g][i][j][l]->GetXaxis()->CenterTitle();
	  
	    jff_residual_phi[g][i][j][l]->Draw();
	    jff_residual_phi[g-6][i][3][l]->Draw("same");
	    l_phi->Draw();
	
	  }
      
	    if(g%2!=0&&g>5){

	      TString save_name_eta = "JFF_Residual_Corrections_Eta_";
	      if(g==9) save_name_eta+="SubLeading_";
	      if(g==11) save_name_eta+="Leading_";
	      save_name_eta+=TrkPtBin_strs[i]; save_name_eta+="_"; save_name_eta+=TrkPtBin_strs[i+1];
	      if(l==2){ save_name_eta+="_AjInclusive";
	      }else{  save_name_eta+="_";  save_name_eta+=AjBin_strs[l]; save_name_eta+="_"; save_name_eta+=AjBin_strs[l+1]; }
	      if(is_recogen) save_name_eta+="_RecoGen";
	      else save_name_eta+="_RecoReco";
	      save_name_eta+=".png";
	      corr_canvas_eta[g][i][l]->SaveAs(save_name_eta);
	    
	      save_name_eta.ReplaceAll(".png",".pdf");
	      corr_canvas_eta[g][i][l]->SaveAs(save_name_eta);

	      TString save_name_phi = "JFF_Residual_Corrections_Phi_";
	      if(g==9) save_name_phi+="SubLeading_";
	      if(g==11) save_name_phi+="Leading_";
	      save_name_phi+=TrkPtBin_strs[i]; save_name_phi+="_"; save_name_phi+=TrkPtBin_strs[i+1]; 
	      if(l==2){	save_name_phi+="_AjInclusive";
	      }else{ save_name_phi+="_";  save_name_phi+=AjBin_strs[l]; save_name_phi+="_"; save_name_phi+=AjBin_strs[l+1]; }
	     
	      if(is_recogen) save_name_phi+="_RecoGen";
	      else save_name_phi+="_RecoReco";
	  
	      save_name_phi+=".png";
	      corr_canvas_phi[g][i][l]->SaveAs(save_name_phi);
	      save_name_phi.ReplaceAll(".png",".pdf");
	      corr_canvas_phi[g][i][l]->SaveAs(save_name_phi);
	    }
	
	    } //j

	}//i
    

      cout<<"starting integrals"<<endl;
      TString integral_eta_pT_name = "integral_eta_pT";
      integral_eta_pT_name+=g;
      integral_eta_pT_name+=l;

      cintegral_eta_pT[g][l] = new TCanvas(integral_eta_pT_name,"",10,10,1500,500);
      cintegral_eta_pT[g][l]->Divide(4,1,0.,0.);


      TString integral_eta_cent_name = "integral_eta_cent";
      integral_eta_cent_name+=g;
      integral_eta_cent_name+=l;

      cintegral_eta_cent[g][l] = new TCanvas(integral_eta_cent_name,"",10,10,2000,500);
      cintegral_eta_cent[g][l]->Divide(5,1,0.,0.);



      for(int j = 0; j<4; j++){

	if(g<6&&j<3)continue;

	in_name = make_name("Result_",g,0,j,l,centlabel,pTlabel,Ajlabel);


	cintegral_eta_pT[g][l]->cd(j+1);

	TString ClosureIntegralEtaPt_name = in_name;
	ClosureIntegralEtaPt_name.ReplaceAll("Result","Closure_Integral_Eta");
	ClosureIntegralEtaPt_name.ReplaceAll("_TrkPt05_TrkPt1","");

	switch(j){
	case 0:
	  Closure_integral_eta_pT[g][j][l] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta0[0],&pTbin_errors[0],&closure_integral_errors[0]);
	  break;
	case 1:
	  Closure_integral_eta_pT[g][j][l] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta1[0],&pTbin_errors[0],&closure_integral_errors[0]);
	  break;
	case 2:
	  Closure_integral_eta_pT[g][j][l] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta2[0],&pTbin_errors[0],&closure_integral_errors[0]);
	  break;
	case 3:
	  Closure_integral_eta_pT[g][j][l] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&Closure_integral_eta3[0],&pTbin_errors[0],&closure_integral_errors[0]);
	  break;

	}

	Closure_integral_eta_pT[g][j][l]->SetName(ClosureIntegralEtaPt_name);
	cout<<g<<ClosureIntegralEtaPt_name<<endl;
      

	Closure_integral_eta_pT[g][j][l]->SetMarkerColor(1);
	Closure_integral_eta_pT[g][j][l]->SetMarkerSize(1);
	Closure_integral_eta_pT[g][j][l]->SetLineColor(1);

	switch(g){
	case 3:
	  Closure_integral_eta_pT[g][j][l]->SetMarkerStyle(34);
	  break;
	case 5: 
	  Closure_integral_eta_pT[g][j][l]->SetMarkerStyle(21);
	  break;
	case 9:
	  Closure_integral_eta_pT[g][j][l]->SetMarkerStyle(34);
	  break;
	case 11: 
	  Closure_integral_eta_pT[g][j][l]->SetMarkerStyle(21);
	  break;
	default:
	  Closure_integral_eta_pT[g][j][l]->SetMarkerStyle(10);
	  break;
	}

	Closure_integral_eta_pT[g][j][l]->SetMinimum(-.5);
	Closure_integral_eta_pT[g][j][l]->SetMaximum(1.5);
	    
	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetRangeUser(.501,7.99);
	Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetNdivisions(306);
	Closure_integral_eta_pT[g][j][l]->Draw("p X A");
	 

	Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetLabelSize(ts);
	   


	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");
	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetTitleSize(ts2);
	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetTitleOffset(xoffset+0.2);
	Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetTitle("(Reco - Gen)/Gen");

	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetNdivisions(8);
   
	
	Closure_integral_eta_pT[g][j][l]->GetXaxis()->CenterTitle();
	Closure_integral_eta_pT[g][j][l]->GetYaxis()->CenterTitle();
	   
	if(j>0){
	  Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetTitleSize(0.0);
	  Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetLabelSize(0.0);
	  Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetTitleSize(ts);
	  Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetLabelSize(ts);
	  Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetTitleOffset(xoffset+0.15);
	}else{
	  Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetLabelSize(ts3);
	  Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetLabelOffset(0.015);
	  Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetTitleOffset(1.);
	  Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetTitleSize(ts2);
	  Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetLabelSize(ts2);
	}



	Closure_integral_eta_pT[g][j][l]->SetMarkerSize(2);



	linePt = new TLine(1.,0,16.,0);
	linePt->SetLineStyle(2);
	linePt->SetLineWidth(1);
	linePt->Draw("same");

      
	if(g!=11)continue;

	
	closure_integral_values.clear();
	closure_integral_errors.clear();
	
	for(int k = 0; k<6; k++){
	  double pt_val, x_val;
	  
	  Closure_integral_eta_pT[g][j][l]->GetPoint(k,x_val,pt_val);
	  closure_integral_values.push_back(pt_val);
	  closure_integral_errors.push_back(pt_val/2.);

	  cout<<pt_val<<endl;
	  
	}


      
	 	  
	Closure_integral_eta_pT2[g][j][l] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&closure_integral_values[0],&pTbin_errors[0],&closure_integral_errors[0]);


	Closure_integral_eta_pT2[11][j][l]->SetFillColor(kOrange-2);
	Closure_integral_eta_pT2[11][j][l]->Draw("same e2");
	  
    

	Closure_integral_eta_pT[5][3][l]->SetMarkerSize(2.5);
	Closure_integral_eta_pT[5][3][l]->SetMarkerColor(kBlue);
	Closure_integral_eta_pT[5][3][l]->SetLineColor(kBlue);
	Closure_integral_eta_pT[5][3][l]->Draw("same p X");


	Closure_integral_eta_pT[3][3][l]->SetMarkerColor(kCyan);
	Closure_integral_eta_pT[3][3][l]->SetLineColor(kCyan);
	Closure_integral_eta_pT[3][3][l]->Draw("same p X");

	Closure_integral_eta_pT[11][j][l]->Draw("same p X");

	Closure_integral_eta_pT[9][j][l]->SetMarkerColor(kGreen+1);
	Closure_integral_eta_pT[9][j][l]->SetLineColor(kGreen+1);
	Closure_integral_eta_pT[9][j][l]->Draw("same p X");


	if(g==11&&j==0){ 
	  l40 = new TLegend(textalign2,texty1-.05,0.8,texty4-.1);
	  l40->SetName("l40");
	  l40->SetTextFont(43);
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);

	  l40->AddEntry(Closure_integral_eta_pT[11][j][l],"Leading P+H","p");
	  l40->AddEntry(Closure_integral_eta_pT[9][j][l],"Subleading P+H","p");


	  l40->AddEntry(Closure_integral_eta_pT[5][3][l],"Leading PYTHIA","p");
	  l40->AddEntry(Closure_integral_eta_pT[3][3][l],"Subleading PYTHIA","p");
	
	  l40->Draw("same");

	}
      
	drawlabels_int_pt2(g,j);
      
      }

      if(g%2==0)continue;
  
      for(int i = 0; i<nTrkPtBins-1; i++){
	
	TString ClosureIntegralEtaCent_name = "ClosureIntegralEtaCent";
	ClosureIntegralEtaCent_name+=g;
	ClosureIntegralEtaCent_name+=i;
	ClosureIntegralEtaCent_name+=l;
	Closure_integral_eta_cent[g][i][l] = new TH1D(ClosureIntegralEtaCent_name,"",4,xAxis);

	for(int k=0; k<4; k++){
	  evalpt = pTbin_centers.at(i);
	  if(g<6) value = Closure_integral_eta_pT[g][3][l]->Eval(evalpt);
	  else	value = Closure_integral_eta_pT[g][k][l]->Eval(evalpt);
	  Closure_integral_eta_cent[g][i][l]->SetBinContent(k+1,value);
	}
  
	switch(g){
	case 1: 
	  Closure_integral_eta_cent[g][i][l]->SetMarkerStyle(10);
	  break;
	case 3:
	  Closure_integral_eta_cent[g][i][l]->SetMarkerStyle(34);
	  break;
	case 5: 
	  Closure_integral_eta_cent[g][i][l]->SetMarkerStyle(21);
	  break;
	case 7: 
	  Closure_integral_eta_cent[g][i][l]->SetMarkerStyle(10);
	  break;
	case 9:
	  Closure_integral_eta_cent[g][i][l]->SetMarkerStyle(34);
	  break;
	case 11: 
	  Closure_integral_eta_cent[g][i][l]->SetMarkerStyle(21);
	  break;
	default:
	  Closure_integral_eta_cent[g][i][l]->SetMarkerStyle(10);
	  break;
	}


	Closure_integral_eta_cent[g][i][l]->SetMarkerSize(2);
	Closure_integral_eta_cent[g][i][l]->SetLineColor(kBlack);

	Closure_integral_eta_cent[g][i][l]->SetLineColor(kBlack);
	Closure_integral_eta_cent[g][i][l]->SetLineColor(kBlack);
	Closure_integral_eta_cent[g][i][l]->GetYaxis()->SetNdivisions(306);

   
	if(g==11){


	  cintegral_eta_cent[g][l]->cd(i+1);


	  TString histnameblank = "blank_hist";
	  histnameblank+=g;
	  histnameblank+=i;
	  histnameblank+=l;

	  blank[i][l] = new TH1D(histnameblank,"",4,xAxis);

	
	  TString histnameblank2 = "blank_hist2";
	  histnameblank2+=g;
	  histnameblank2+=i;

	
	  blank[i][l]->SetMinimum(-.3);
	  blank[i][l]->SetMaximum(.5);
	  blank[i][l]->GetXaxis()->SetTitle("Centrality (%)");
	  blank[i][l]->GetXaxis()->SetTitleOffset(1.1);
	  blank[i][l]->GetXaxis()->CenterTitle(true);
	  blank[i][l]->GetXaxis()->SetTitleSize(ts);

	  blank[i][l]->GetYaxis()->SetTitle("(dN/dp_{T})_{PbPb}- (dN/dp_{T})_{pp} (GeV/c)^{-1}");
	  blank[i][l]->GetYaxis()->SetTitleSize(0.);
	  blank[i][l]->GetYaxis()->CenterTitle(true);
	  //	blank[i][l]->GetYaxis()->SetLabelOffset(yoffset);
	  //	blank[i][l]->GetYaxis()->SetLabelSize(0.);
   
	  blank[i][l]->GetYaxis()->SetTickLength(0.025);

	  blank[i][l]->GetXaxis()->SetBinLabel(1,"50-100");
	  blank[i][l]->GetXaxis()->SetBinLabel(2,"30-50");
	  blank[i][l]->GetXaxis()->SetBinLabel(3,"10-30");
	  blank[i][l]->GetXaxis()->SetBinLabel(4," 0-10");
    
	  blank[i][l]->GetXaxis()->SetLabelSize(0.08);
	  blank[i][l]->GetXaxis()->SetLabelOffset(0.015);
	
	  blank[i][l]->GetXaxis()->LabelsOption("h");
	  blank[i][l]->GetXaxis()->SetTickLength(0.0);



	  switch(i){
	  case 0: 
	    //  gPad->SetLeftMargin(0.2);
	    blank[i][l]->GetYaxis()->SetTitleSize(ts);
	    blank[i][l]->GetXaxis()->SetTitleOffset(1.1);
	    blank[i][l]->GetXaxis()->SetTitleSize(0.07);
	    blank[i][l]->SetLabelSize(0.95*blank[i][l]->GetXaxis()->GetLabelSize());
	    blank[i][l]->GetYaxis()->SetLabelSize(ts2);
	    break;
	  case 3:
	    // gPad->SetRightMargin(0.02);
	    break;
	  default:
	    break;
	  }

	  //----------------------------------

	  blank[i][l]->GetXaxis()->SetTitle("Centrality (%)");
	  blank[i][l]->GetXaxis()->SetTitleOffset(1.1);
	  blank[i][l]->GetXaxis()->CenterTitle(true);
	  blank[i][l]->GetXaxis()->SetTitleSize(ts);

	  blank[i][l]->GetYaxis()->SetTitle("(dN/dp_{T})_{PbPb}- (dN/dp_{T})_{pp} (GeV/c)^{-1}");
	  blank[i][l]->GetYaxis()->SetTitleSize(0.);
	  blank[i][l]->GetYaxis()->CenterTitle(true);
	  //	blank[i][l]->GetYaxis()->SetLabelOffset(yoffset);
	  //	blank[i][l]->GetYaxis()->SetLabelSize(0.);
   
	  blank[i][l]->GetYaxis()->SetTickLength(0.025);

	  blank[i][l]->GetXaxis()->SetBinLabel(1,"50-100");
	  blank[i][l]->GetXaxis()->SetBinLabel(2,"30-50");
	  blank[i][l]->GetXaxis()->SetBinLabel(3,"10-30");
	  blank[i][l]->GetXaxis()->SetBinLabel(4," 0-10");
    
	  blank[i][l]->GetXaxis()->SetLabelSize(0.08);
	  blank[i][l]->GetXaxis()->SetLabelOffset(0.015);
	
	  blank[i][l]->GetXaxis()->LabelsOption("h");
	  blank[i][l]->GetXaxis()->SetTickLength(0.0);

	  switch(i){
	  case 0: 
	    //  gPad->SetLeftMargin(0.2);
	    blank[i][l]->GetYaxis()->SetTitleSize(ts);
	    blank[i][l]->GetXaxis()->SetTitleOffset(1.1);
	    blank[i][l]->GetXaxis()->SetTitleSize(0.07);
	    blank[i][l]->SetLabelSize(0.95*blank[i][l]->GetXaxis()->GetLabelSize());
	    blank[i][l]->GetYaxis()->SetLabelSize(ts2);
	    break;
	  default:
	    blank[i][l]->GetYaxis()->SetLabelSize(0.);
	    break;
	  }



	  blank[i][l]->Draw();


    

	  Closure_integral_eta_cent[5][i][l]->SetMarkerSize(2.5);
	  Closure_integral_eta_cent[5][i][l]->SetMarkerColor(kBlue);
	  Closure_integral_eta_cent[5][i][l]->Draw("same p X");

	  Closure_integral_eta_cent[3][i][l]->SetMarkerColor(kCyan);
	  Closure_integral_eta_cent[3][i][l]->Draw("same p X");

	  Closure_integral_eta_cent[11][i][l]->Draw("same p X");

	  Closure_integral_eta_cent[9][i][l]->SetMarkerColor(kGreen+1);
	  Closure_integral_eta_cent[9][i][l]->Draw("same p X");




	  TLatex *pt_label = new TLatex(0.2,0.85, TrkPtBin_labels[i]);
	  pt_label->SetNDC();
	  pt_label->SetTextSizePixels(tspixels);
	  pt_label->Draw();

	  TLatex *aj_label = new TLatex(0.2,0.8, Ajlabel);
	  aj_label->SetNDC();
	  aj_label->SetTextSizePixels(tspixels);
	  aj_label->Draw();
	}

      }
    
      if(g!=11)continue;			
			       
      cintegral_eta_pT[g][l]->cd(0);
								      
      TLatex *canvas_title = new TLatex(0.06,0.9,"CMS Preliminary Simulation");
      canvas_title->SetTextSizePixels(tspixels);
      canvas_title->SetTextFont(63);
      canvas_title->Draw();

      TLatex *canvas_title2 = new TLatex(0.295,0.9,"PYTHIA+HYDJET");
      canvas_title2->SetTextSizePixels(tspixels);
      canvas_title2->Draw();

      if(l==2){

	cintegral_eta_pT[g][l]->SaveAs("Integral_Closure_pT_Leading_AjInclusive.pdf");
	cintegral_eta_pT[g][l]->SaveAs("Integral_Closure_pT_Leading_AjInclusive.png");

	cintegral_eta_cent[g][l]->SaveAs("Integral_Closure_Cent_Leading_AjInclusive.pdf");
	cintegral_eta_cent[g][l]->SaveAs("Integral_Closure_Cent_Leading_AjInclusive.png");
   
      }else{
	cintegral_eta_pT[g][l]->SaveAs("Integral_Closure_pT_Leading_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf");
	cintegral_eta_pT[g][l]->SaveAs("Integral_Closure_pT_Leading_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png");
	cintegral_eta_cent[g][l]->SaveAs("Integral_Closure_Cent_Leading_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf");
	cintegral_eta_cent[g][l]->SaveAs("Integral_Closure_Cent_Leading_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png");
   
      }

    
    }//l
  }//g
  
  return 0;
}
