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
#include "TGaxis.h"
#include "TStyle.h"
#include "TLatex.h"
#include "THStack.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "../JetTrack2015_functions.h"


using namespace std;

Int_t plot_results(){

  gROOT->ForceStyle();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);
  gStyle->SetHatchesLineWidth(2);

  // gStyle->SetTextFont(63);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.25);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  gStyle->SetNdivisions(510);

 
  TH1D *signal_dPhi_rebin[12][6][4][5];
  TH1D *signal_dPhi_tot[12][6][4][5];
  TH1D *signal_dPhi_syst[12][6][4][5];
  TH1D *signal_dPhi_noerr_up[12][6][4][5];
  TH1D *signal_dPhi_noerr_down[12][6][4][5];

  TH1D *sub_lead_dPhi_rebin[12][6][4][5];
  TH1D *sub_lead_dPhi_syst[12][6][4][5];
  TH1D *sub_lead_dPhi_tot[12][6][4][5];
  TH1D *sub_lead_dPhi_noerr_up[12][6][4][5];
  TH1D *sub_lead_dPhi_noerr_down[12][6][4][5];

  TH1D *blank[12][6][4][5];  

  TH1D *background_diff_rebin[12][6][4][5];
  TH1D *background_syst_rebin[12][6][4][5];
  TH1D *background_diff_tot[12][6][4][5];
  TH1D *background_diff_noerr_up[12][6][4][5];
  TH1D *background_diff_noerr_down[12][6][4][5];

  TH1D *signal_dEta_rebin[12][6][4][5];

  THStack *dPhi_pTdist_up[12][4][3];
  THStack *dPhi_Sub_Lead_NoBkg_up[12][4][3];
  THStack *dPhi_pTdist_down[12][4][3];
  THStack *dPhi_Sub_Lead_NoBkg_down[12][4][3];

  THStack *Integral_up[12][4][3];
  THStack *Integral_down[12][4][3];


  TH1D *Integral[12][4][3];
  TH1D *Integral_syst[12][4][3];
  TH1D *Integral_noerr_up[12][4][3];
  TH1D *Integral_noerr_down[12][4][3];


  TCanvas *c_jet[5];
  TCanvas *c_nobkgsub[5];
  TCanvas *c_longrange[5];
  TCanvas *c_longrange_int;
  TCanvas *c_all_int;



  TFile *f_in = new TFile("../study_yield/Dijet_Results.root");
  TFile *f_in_nobkgsub = new TFile("../study_yield/Dijet_Results_NoBkgSub.root");

  TFile *f_out = new TFile("Systematic_Error.root","RECREATE");
 
  TString in_name, stack_name, tot_name, pTlabel,centlabel,Ajlabel;

  TString AjBin_strs[4] = {"Aj0","Aj22","Aj75","AjInclusive"};
  TString AjName_strs[4] = {"Aj0_Aj22_","Aj22_Aj75_",""};

  int lbin, rbin;

  float bin_width_phi, bin_width_eta, diff_max, diff_min, double_diff_max, double_diff_min, no_bgsub_min, no_bgsub_max, signal_min, signal_max;
  double integral, int_err;

  float pt_bin_centers[6] = {0.75,1.5,2.5,3.5,6.,12.};
  float pt_bin_bounds[7] = {0.5,1.,2.,3.,4.,8.,16.};
  float pt_bin_errors[6] = {0.25,.5,.5,.5,2.,4.};

  TLatex *label_pp, *label_per, *label_cent;
  float bc_pbpb, bc_pp, err_change,new_err, old_err;

  //First, get all histos
  
  for(int l = 0; l<3; l++){

    for(int i = 0; i<6; i++){

      for(int j = 0; j<4; j++){


	TString syst_name_pbpb = make_name("Syst_",0,i,j,l,pTlabel, centlabel, Ajlabel);  
	background_syst_rebin[0][i][j][l] = (TH1D*)f_in->Get(syst_name_pbpb)->Clone(syst_name_pbpb);

	
	TString syst_name_pp = make_name("Syst_",1,i,3,l,pTlabel, centlabel, Ajlabel);  
	background_syst_rebin[1][i][j][l] = (TH1D*)f_in->Get(syst_name_pp)->Clone(syst_name_pp);


	if(i==4){
	  for(int k = 0; k<4; k++){
	    background_syst_rebin[0][i][j][l]->Add( background_syst_rebin[0][k][j][l]);
	    background_syst_rebin[1][i][j][l]->Add( background_syst_rebin[1][k][j][l]);

	  }

	  syst_name_pbpb.ReplaceAll("Syst","Syst_Diff");
	  background_syst_rebin[7][i][j][l] = (TH1D*)background_syst_rebin[0][i][j][l]->Clone(syst_name_pbpb);
	  background_syst_rebin[7][i][j][l]->Add(background_syst_rebin[1][i][j][l],-1.);
	}



	for(int g = 0; g<2; g++){
	  
	  //  if(j!=3&&g==1)continue;

	
	  in_name = make_name("Diff_",g,i,j,l,pTlabel,centlabel,Ajlabel);
	  in_name+="_Rebin";
	  background_diff_rebin[g][i][j][l] = (TH1D*)f_in->Get(in_name)->Clone(in_name);

	  in_name.ReplaceAll("Diff","Diff_NoErrUp");
	  background_diff_noerr_up[g][i][j][l] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	  in_name.ReplaceAll("Up","Down");
	  background_diff_noerr_down[g][i][j][l] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);


	  in_name = make_name("SubLeading_dPhi_",g,i,j,l,pTlabel,centlabel,Ajlabel);
	  in_name+="_Rebin";

	  signal_dPhi_rebin[g+2][i][j][l] = (TH1D*)f_in->Get(in_name)->Clone(in_name);

	  in_name.ReplaceAll("SubLeading","Leading");
	  signal_dPhi_rebin[g+4][i][j][l] = (TH1D*)f_in->Get(in_name)->Clone(in_name);
	  
	  in_name.ReplaceAll("Phi","Eta");
	  signal_dEta_rebin[g+4][i][j][l] = (TH1D*)f_in->Get(in_name)->Clone(in_name);

	  in_name.ReplaceAll("Leading","SubLeading");
	  signal_dEta_rebin[g+2][i][j][l] = (TH1D*)f_in->Get(in_name)->Clone(in_name);

	  in_name = make_name("Proj_dPhi_",g,i,j,l,pTlabel,centlabel,Ajlabel);
	  in_name+="_Rebin";
	  sub_lead_dPhi_rebin[g][i][j][l] = (TH1D*)f_in_nobkgsub->Get(in_name)->Clone(in_name);
	  
	  in_name.ReplaceAll("dPhi","dPhi_NoErrUp");
	  sub_lead_dPhi_noerr_up[g][i][j][l] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	  in_name.ReplaceAll("Up","Down");
	  sub_lead_dPhi_noerr_down[g][i][j][l] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);


	  in_name = make_name("dPhi_NoErrUp_",g+2,i,j,l,pTlabel,centlabel,Ajlabel);
	  signal_dPhi_noerr_up[g+2][i][j][l] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	  in_name.ReplaceAll("Up","Down");
	  signal_dPhi_noerr_down[g+2][i][j][l] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);

	  in_name.ReplaceAll("SubLeading","Leading");
	  signal_dPhi_noerr_down[g+4][i][j][l] = new TH1D(in_name,"",nbounds_phi-1,bin_bounds_phi);
	  in_name.ReplaceAll("Down","Up");
	  signal_dPhi_noerr_up[g+4][i][j][l] = new TH1D(in_name,"",nbounds_phi-1 ,bin_bounds_phi);

	  for(int k = 1; k<26; k++){
	    if(background_diff_rebin[g][i][j][l]->GetBinContent(k)>0) background_diff_noerr_up[g][i][j][l]->SetBinContent(k, background_diff_rebin[g][i][j][l]->GetBinContent(k));
	    else background_diff_noerr_down[g][i][j][l]->SetBinContent(k, background_diff_rebin[g][i][j][l]->GetBinContent(k));

	    if(	signal_dPhi_rebin[g+2][i][j][l]->GetBinContent(k)>0) signal_dPhi_noerr_up[g+2][i][j][l]->SetBinContent(k,signal_dPhi_rebin[g+2][i][j][l]->GetBinContent(k));
	    else signal_dPhi_noerr_down[g+2][i][j][l]->SetBinContent(k,signal_dPhi_rebin[g+2][i][j][l]->GetBinContent(k));

	    if(signal_dPhi_rebin[g+4][i][j][l]->GetBinContent(k)>0)   signal_dPhi_noerr_up[g+4][i][j][l]->SetBinContent(k,signal_dPhi_rebin[g+4][i][j][l]->GetBinContent(k));
	    else   signal_dPhi_noerr_down[g+4][i][j][l]->SetBinContent(k,signal_dPhi_rebin[g+4][i][j][l]->GetBinContent(k));
	
	    if(sub_lead_dPhi_rebin[g][i][j][l]->GetBinContent(k)>0) sub_lead_dPhi_noerr_up[g][i][j][l]->SetBinContent(k,sub_lead_dPhi_rebin[g][i][j][l]->GetBinContent(k));
	    else  sub_lead_dPhi_noerr_down[g][i][j][l]->SetBinContent(k,sub_lead_dPhi_rebin[g][i][j][l]->GetBinContent(k));
	  }
	

	  switch(i){
	  case 0:
	    signal_dPhi_noerr_up[g+2][i][j][l]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_up[g+4][i][j][l]->SetFillColor(kBlue-9);
	    background_diff_noerr_up[g][i][j][l]->SetFillColor(kBlue-9);
	    sub_lead_dPhi_noerr_up[g][i][j][l]->SetFillColor(kBlue-9);

	    signal_dPhi_noerr_down[g+2][i][j][l]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_down[g+4][i][j][l]->SetFillColor(kBlue-9);
	    background_diff_noerr_down[g][i][j][l]->SetFillColor(kBlue-9);
	    sub_lead_dPhi_noerr_down[g][i][j][l]->SetFillColor(kBlue-9);

	    break;
	  case 1:
	    signal_dPhi_noerr_up[g+2][i][j][l]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_up[g+4][i][j][l]->SetFillColor(kYellow-9);
	    background_diff_noerr_up[g][i][j][l]->SetFillColor(kYellow-9);
	    sub_lead_dPhi_noerr_up[g][i][j][l]->SetFillColor(kYellow-9);

	    signal_dPhi_noerr_down[g+2][i][j][l]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_down[g+4][i][j][l]->SetFillColor(kYellow-9);
	    background_diff_noerr_down[g][i][j][l]->SetFillColor(kYellow-9);
	    sub_lead_dPhi_noerr_down[g][i][j][l]->SetFillColor(kYellow-9);
	    break;
	  case 2:
	    signal_dPhi_noerr_up[g+2][i][j][l]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_up[g+4][i][j][l]->SetFillColor(kOrange+1);
	    background_diff_noerr_up[g][i][j][l]->SetFillColor(kOrange+1);
	    sub_lead_dPhi_noerr_up[g][i][j][l]->SetFillColor(kOrange+1);

	    signal_dPhi_noerr_down[g+2][i][j][l]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_down[g+4][i][j][l]->SetFillColor(kOrange+1);
	    background_diff_noerr_down[g][i][j][l]->SetFillColor(kOrange+1);
	    sub_lead_dPhi_noerr_down[g][i][j][l]->SetFillColor(kOrange+1);
	    break;
	  case 3:
	    signal_dPhi_noerr_up[g+2][i][j][l]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_up[g+4][i][j][l]->SetFillColor(kViolet-5);
	    background_diff_noerr_up[g][i][j][l]->SetFillColor(kViolet-5);
	    sub_lead_dPhi_noerr_up[g][i][j][l]->SetFillColor(kViolet-5);

	    signal_dPhi_noerr_down[g+2][i][j][l]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_down[g+4][i][j][l]->SetFillColor(kViolet-5);
	    background_diff_noerr_down[g][i][j][l]->SetFillColor(kViolet-5);
	    sub_lead_dPhi_noerr_down[g][i][j][l]->SetFillColor(kViolet-5);
	    break;
	  case 4:
	    signal_dPhi_noerr_up[g+2][i][j][l]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_up[g+4][i][j][l]->SetFillColor(kGreen+3);
	    background_diff_noerr_up[g][i][j][l]->SetFillColor(kGreen+3);
	    sub_lead_dPhi_noerr_up[g][i][j][l]->SetFillColor(kGreen+3);

	    signal_dPhi_noerr_down[g+2][i][j][l]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_down[g+4][i][j][l]->SetFillColor(kGreen+3);
	    background_diff_noerr_down[g][i][j][l]->SetFillColor(kGreen+3);
	    sub_lead_dPhi_noerr_down[g][i][j][l]->SetFillColor(kGreen+3);
	    break;
	  case 5:
	    signal_dPhi_noerr_up[g+2][i][j][l]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_up[g+4][i][j][l]->SetFillColor(kRed+1);
	    background_diff_noerr_up[g][i][j][l]->SetFillColor(kRed+1);
	    sub_lead_dPhi_noerr_up[g][i][j][l]->SetFillColor(kRed+1);

	    signal_dPhi_noerr_down[g+2][i][j][l]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_down[g+4][i][j][l]->SetFillColor(kRed+1);
	    background_diff_noerr_down[g][i][j][l]->SetFillColor(kRed+1);
	    sub_lead_dPhi_noerr_down[g][i][j][l]->SetFillColor(kRed+1);
	    break;
	  default:
	    break;
	  }

	  signal_dPhi_noerr_up[g+2][i][j][l]->SetFillStyle(1001);
	  signal_dPhi_noerr_up[g+4][i][j][l]->SetFillStyle(1001);
	  background_diff_noerr_up[g][i][j][l]->SetFillStyle(1001);
	  sub_lead_dPhi_noerr_up[g][i][j][l]->SetFillStyle(1001);
	  
	  signal_dPhi_noerr_down[g+2][i][j][l]->SetFillStyle(1001);
	  signal_dPhi_noerr_down[g+4][i][j][l]->SetFillStyle(1001);
	  background_diff_noerr_down[g][i][j][l]->SetFillStyle(1001);
	  sub_lead_dPhi_noerr_down[g][i][j][l]->SetFillStyle(1001);
	  

    	}
      }
    }
  }


  for(int l = 0; l<3; l++){
   
    c_jet[l] = new TCanvas(Form("MomentumBalanceJet%d",l),"",10,10,1200,1600);
    c_jet[l]->Divide(3,4,0,0);

    c_longrange[l] = new TCanvas(Form("MomentumBalanceLongRange%d",l),"",10,10,1200,800);
    c_longrange[l]->Divide(3,2,0,0);
  
    c_nobkgsub[l] = new TCanvas(Form("MomentumBalanceNoBkgSub%d",l),"",10,10,1200,800);
    c_nobkgsub[l]->Divide(3,2,0,0);
    

    cout<<"here..."<<endl;


    for(int j = 0; j<4; j++){
      
      tot_name = make_name("AllHists_Summed_",2,0,j,l,pTlabel,centlabel,Ajlabel);
      signal_dPhi_tot[2][0][j][l] = (TH1D*)signal_dPhi_rebin[2][0][j][l]->Clone(tot_name);
    
      tot_name.ReplaceAll("SubLeading","Leading");
      signal_dPhi_tot[4][0][j][l] = (TH1D*)signal_dPhi_rebin[4][0][j][l]->Clone(tot_name);

      tot_name = make_name("AllHists_BackgroundDiff_",0,0,j,l,pTlabel,centlabel,Ajlabel);
      background_diff_tot[0][0][j][l] = (TH1D*)background_diff_rebin[0][0][j][l]->Clone(tot_name);
   
      tot_name = make_name("AllHists_NoBkgSub_",0,0,j,l,pTlabel,centlabel,Ajlabel);
      sub_lead_dPhi_tot[0][0][j][l] = (TH1D*)sub_lead_dPhi_rebin[0][0][j][l]->Clone(tot_name);
    
   
      for(int k = 1; k<5; k++){
	signal_dPhi_tot[2][0][j][l]->Add(signal_dPhi_rebin[2][k][j][l]);
	signal_dPhi_tot[4][0][j][l]->Add(signal_dPhi_rebin[4][k][j][l]);
	background_diff_tot[0][0][j][l]->Add( background_diff_rebin[0][k][j][l]);
	sub_lead_dPhi_tot[0][0][j][l]->Add(sub_lead_dPhi_rebin[0][k][j][l]);
      }

      signal_dPhi_tot[2][0][j][l]->SetMarkerStyle(4);
      signal_dPhi_tot[4][0][j][l]->SetMarkerStyle(4);
      background_diff_tot[0][0][j][l]->SetMarkerStyle(4);
      sub_lead_dPhi_tot[0][0][j][l]->SetMarkerStyle(4);

    }

    tot_name.ReplaceAll("PbPb","pp");
    signal_dPhi_tot[5][0][3][l] = (TH1D*)signal_dPhi_rebin[5][0][3][l]->Clone(tot_name);

    tot_name.ReplaceAll("Leading","SubLeading");
    signal_dPhi_tot[3][0][3][l] = (TH1D*)signal_dPhi_rebin[3][0][3][l]->Clone(tot_name);


    tot_name = make_name("AllHists_BackgroundDiff_",1,0,3,l,pTlabel,centlabel,Ajlabel);
    background_diff_tot[1][0][3][l] = (TH1D*)background_diff_rebin[1][0][3][l]->Clone(tot_name);
   
    tot_name = make_name("AllHists_NoBkgSub_",1,0,3,l,pTlabel,centlabel,Ajlabel);
    sub_lead_dPhi_tot[1][0][3][l] = (TH1D*)sub_lead_dPhi_rebin[1][0][3][l]->Clone(tot_name);
    
   
    for(int k = 1; k<5; k++){
      signal_dPhi_tot[3][0][3][l]->Add(signal_dPhi_rebin[3][k][3][l]);
      signal_dPhi_tot[5][0][3][l]->Add(signal_dPhi_rebin[5][k][3][l]);
      background_diff_tot[1][0][3][l]->Add(background_diff_rebin[1][k][3][l]);
      sub_lead_dPhi_tot[1][0][3][l]->Add( sub_lead_dPhi_rebin[1][k][3][l]);
    }


    signal_dPhi_tot[3][0][3][l]->SetMarkerStyle(4);
    signal_dPhi_tot[5][0][3][l]->SetMarkerStyle(4);
    background_diff_tot[1][0][3][l]->SetMarkerStyle(4);
    sub_lead_dPhi_tot[1][0][3][l]->SetMarkerStyle(4);
 
 

    signal_dPhi_syst[2][0][0][l] = (TH1D*)f_in->Get((TString)("SubLeading_dPhi_Syst_PbPb_Cent50_Cent100_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("SubLeading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"));
  

    signal_dPhi_syst[4][0][0][l] = (TH1D*)f_in->Get((TString)("Leading_dPhi_Syst_PbPb_Cent50_Cent100_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("Leading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"));

    signal_dPhi_syst[2][0][2][l] = (TH1D*)f_in->Get((TString)("SubLeading_dPhi_Syst_PbPb_Cent10_Cent30_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("SubLeading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent10_Cent30"));

    signal_dPhi_syst[4][0][2][l] = (TH1D*)f_in->Get((TString)("Leading_dPhi_Syst_PbPb_Cent10_Cent30_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("Leading_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent10_Cent30"));


     signal_dPhi_syst[3][0][3][l] = (TH1D*)f_in->Get((TString)("SubLeading_dPhi_Syst_pp_Cent0_Cent10_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("SubLeading_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"));
     signal_dPhi_syst[5][0][3][l] = (TH1D*)f_in->Get((TString)("Leading_dPhi_Syst_pp_Cent0_Cent10_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("Leading_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"));
   


    sub_lead_dPhi_syst[2][0][0][l] = (TH1D*)f_in_nobkgsub->Get((TString)("SubLeading_dPhi_Syst_PbPb_Cent50_Cent100_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("Diff_NoBkgSub_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent50_Cent100"));
  
    sub_lead_dPhi_syst[2][0][2][l] = (TH1D*)f_in_nobkgsub->Get((TString)("SubLeading_dPhi_Syst_PbPb_Cent10_Cent30_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("Diff_NoBkgSub_dPhi_Syst_PbPb_"+AjName_strs[l]+"Cent10_Cent30"));

     sub_lead_dPhi_syst[3][0][3][l] = (TH1D*)f_in_nobkgsub->Get((TString)("SubLeading_dPhi_Syst_pp_Cent0_Cent10_Pt100_Pt300_"+AjName_strs[l]+"TrkPt4_TrkPt8"))->Clone((TString)("Diff_NoBkgSub_dPhi_Syst_pp_"+AjName_strs[l]+"Cent0_Cent10"));

  
    signal_dPhi_syst[2][0][0][l]->SetMarkerStyle(20);
    signal_dPhi_syst[4][0][0][l]->SetMarkerStyle(20);
    signal_dPhi_syst[2][0][2][l]->SetMarkerStyle(20);
    signal_dPhi_syst[4][0][2][l]->SetMarkerStyle(20);
    signal_dPhi_syst[3][0][3][l]->SetMarkerStyle(20);
    signal_dPhi_syst[5][0][3][l]->SetMarkerStyle(20);

    cout<<"and here"<<endl;

    background_syst_rebin[0][4][2][l]->SetMarkerStyle(20);
    background_syst_rebin[0][4][0][l]->SetMarkerStyle(20);
    background_syst_rebin[1][4][3][l]->SetMarkerStyle(20);
    background_syst_rebin[7][4][2][l]->SetMarkerStyle(20);
    background_syst_rebin[7][4][0][l]->SetMarkerStyle(20);
  
    cout<<"also here"<<endl;
    
    sub_lead_dPhi_syst[2][0][0][l]->SetMarkerStyle(20);
    sub_lead_dPhi_syst[2][0][2][l]->SetMarkerStyle(20);
    sub_lead_dPhi_syst[3][0][3][l]->SetMarkerStyle(20);
 

    signal_dPhi_syst[2][0][0][l]->SetMarkerColor(kWhite);
    signal_dPhi_syst[4][0][0][l]->SetMarkerColor(kWhite);
    signal_dPhi_syst[2][0][2][l]->SetMarkerColor(kWhite);
    signal_dPhi_syst[4][0][2][l]->SetMarkerColor(kWhite);
    signal_dPhi_syst[3][0][3][l]->SetMarkerColor(kWhite);
    signal_dPhi_syst[5][0][3][l]->SetMarkerColor(kWhite);


    background_syst_rebin[0][4][2][l]->SetMarkerColor(kWhite);
    background_syst_rebin[0][4][0][l]->SetMarkerColor(kWhite);
    background_syst_rebin[1][4][3][l]->SetMarkerColor(kWhite);
    background_syst_rebin[7][4][2][l]->SetMarkerColor(kWhite);
    background_syst_rebin[7][4][0][l]->SetMarkerColor(kWhite);


    sub_lead_dPhi_syst[2][0][0][l]->SetMarkerColor(kWhite);
    sub_lead_dPhi_syst[2][0][2][l]->SetMarkerColor(kWhite);
    sub_lead_dPhi_syst[3][0][3][l]->SetMarkerColor(kWhite);

  
  
    signal_dPhi_syst[2][0][0][l]->SetFillStyle(3004);
    signal_dPhi_syst[4][0][0][l]->SetFillStyle(3004);
    signal_dPhi_syst[2][0][2][l]->SetFillStyle(3004);
    signal_dPhi_syst[4][0][2][l]->SetFillStyle(3004);
    signal_dPhi_syst[3][0][3][l]->SetFillStyle(3004);
    signal_dPhi_syst[5][0][3][l]->SetFillStyle(3004);
    background_syst_rebin[0][4][2][l]->SetFillStyle(3004);
    background_syst_rebin[0][4][0][l]->SetFillStyle(3004);
    background_syst_rebin[1][4][3][l]->SetFillStyle(3004);
    background_syst_rebin[7][4][2][l]->SetFillStyle(3004);
    background_syst_rebin[7][4][0][l]->SetFillStyle(3004);


    sub_lead_dPhi_syst[2][0][0][l]->SetFillStyle(3004);
    sub_lead_dPhi_syst[2][0][2][l]->SetFillStyle(3004);
    sub_lead_dPhi_syst[3][0][3][l]->SetFillStyle(3004);

  
 
    signal_dPhi_syst[2][0][0][l]->SetFillColor(kBlack);
    signal_dPhi_syst[4][0][0][l]->SetFillColor(kBlack);
    signal_dPhi_syst[2][0][2][l]->SetFillColor(kBlack);
    signal_dPhi_syst[4][0][2][l]->SetFillColor(kBlack);
    signal_dPhi_syst[3][0][3][l]->SetFillColor(kBlack);
    signal_dPhi_syst[5][0][3][l]->SetFillColor(kBlack);
    background_syst_rebin[0][4][0][l]->SetFillColor(kBlack);
    background_syst_rebin[0][4][2][l]->SetFillColor(kBlack);
    background_syst_rebin[1][4][3][l]->SetFillColor(kBlack);
    background_syst_rebin[7][4][2][l]->SetFillColor(kBlack);
    background_syst_rebin[7][4][0][l]->SetFillColor(kBlack);


    sub_lead_dPhi_syst[2][0][0][l]->SetFillColor(kBlack);
    sub_lead_dPhi_syst[2][0][2][l]->SetFillColor(kBlack);
    sub_lead_dPhi_syst[3][0][3][l]->SetFillColor(kBlack);



    
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





    in_name = make_name("PbPb_pp_WithErrors_Syst_",8,0,3,l,pTlabel,centlabel,Ajlabel);
    signal_dPhi_syst[8][0][2][l] = (TH1D*)signal_dPhi_syst[2][0][2][l]->Clone(in_name);
    signal_dPhi_syst[8][0][2][l]->Add( signal_dPhi_syst[3][0][3][l],-1.);

    in_name.ReplaceAll("SubLeading","Leading");
    signal_dPhi_syst[10][0][2][l] = (TH1D*)signal_dPhi_syst[4][0][2][l]->Clone(in_name);
    signal_dPhi_syst[10][0][2][l]->Add( signal_dPhi_syst[5][0][3][l],-1.);



  
   
    for(int k = 1; k< signal_dPhi_syst[8][0][2][l]->GetNbinsX()+1; k++){
      
      bc_pbpb =  signal_dPhi_syst[2][0][2][l]->GetBinContent(k);
      bc_pp =  signal_dPhi_syst[3][0][3][l]->GetBinContent(k);
      err_change = TMath::Sqrt(0.04*0.04+0.05*0.05)*(TMath::Sqrt(bc_pbpb*bc_pbpb+bc_pp*bc_pp)-TMath::Abs(bc_pbpb - bc_pp));
      old_err = signal_dPhi_syst[8][0][2][l]->GetBinError(k);
      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
      signal_dPhi_syst[8][0][2][l]->SetBinError(k,new_err);
    
      bc_pbpb =  signal_dPhi_syst[4][0][2][l]->GetBinContent(k);
      bc_pp =  signal_dPhi_syst[5][0][3][l]->GetBinContent(k);
      err_change = TMath::Sqrt(0.04*0.04+0.05*0.05)*(TMath::Sqrt(bc_pbpb*bc_pbpb+bc_pp*bc_pp)-TMath::Abs(bc_pbpb - bc_pp)); 
      old_err = signal_dPhi_syst[10][0][2][l]->GetBinError(k);
      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
      signal_dPhi_syst[10][0][2][l]->SetBinError(k,new_err);

    }





    in_name.ReplaceAll("SubLeading","DoubleDiff");
    signal_dPhi_syst[6][0][2][l] = (TH1D*)signal_dPhi_syst[8][0][2][l]->Clone(in_name);
    signal_dPhi_syst[6][0][2][l]->Add( signal_dPhi_syst[10][0][2][l]);
 


   
    for(int k = 1; k< signal_dPhi_syst[8][0][2][l]->GetNbinsX()+1; k++){
      
      bc_pbpb =  signal_dPhi_syst[2][0][2][l]->GetBinContent(k);
      bc_pp =  signal_dPhi_syst[4][0][2][l]->GetBinContent(k);
      err_change = TMath::Sqrt(0.04*0.04+0.05*0.05)*(bc_pbpb-TMath::Abs(signal_dPhi_syst[6][0][2][l]->GetBinContent(k)));
    
      old_err = signal_dPhi_syst[6][0][2][l]->GetBinError(k);
      new_err = TMath::Sqrt(old_err*old_err-err_change*err_change);
      signal_dPhi_syst[6][0][2][l]->SetBinError(k,new_err);
    
    }




    in_name = make_name("Double_Diff_NoBkgSub_Syst_",8,0,0,l,pTlabel,centlabel,Ajlabel);
    sub_lead_dPhi_syst[8][0][0][l] = (TH1D*)sub_lead_dPhi_syst[2][0][0][l]->Clone(in_name);
    sub_lead_dPhi_syst[8][0][0][l]->Add( sub_lead_dPhi_syst[3][0][3][l],-1.);


    in_name = make_name("Double_Diff_NoBkgSub_Syst_",8,0,3,l,pTlabel,centlabel,Ajlabel);
    sub_lead_dPhi_syst[8][0][2][l] = (TH1D*) sub_lead_dPhi_syst[2][0][2][l]->Clone(in_name);
    sub_lead_dPhi_syst[8][0][2][l]->Add( sub_lead_dPhi_syst[3][0][3][l],-1.);



    tot_name = make_name("AllHists_PbPb_pp_",2,0,0,l,pTlabel,centlabel,Ajlabel);
    signal_dPhi_tot[8][0][0][l] = (TH1D*) signal_dPhi_tot[2][0][0][l]->Clone(tot_name);
    signal_dPhi_tot[8][0][0][l]->Add( signal_dPhi_tot[3][0][3][l],-1. );

    tot_name.ReplaceAll("SubLeading","Leading");
    signal_dPhi_tot[10][0][0][l] = (TH1D*) signal_dPhi_tot[4][0][0][l]->Clone(tot_name);
    signal_dPhi_tot[10][0][0][l]->Add( signal_dPhi_tot[5][0][3][l],-1. );

    tot_name.ReplaceAll("Leading","DoubleDiff");
 
    signal_dPhi_tot[6][0][0][l] = (TH1D*) signal_dPhi_tot[8][0][0][l]->Clone(tot_name);
    signal_dPhi_tot[6][0][0][l]->Add( signal_dPhi_tot[10][0][0][l] );


    tot_name = make_name("BackgroundDiff_PbPb_pp_",2,0,0,l,pTlabel,centlabel,Ajlabel);
    background_diff_tot[7][0][0][l] = (TH1D*) background_diff_tot[0][0][0][l]->Clone(tot_name);
    background_diff_tot[7][0][0][l]->Add( background_diff_tot[1][0][3][l],-1. );



    tot_name = make_name("Sub_Lead_PbPb_pp_",2,0,0,l,pTlabel,centlabel,Ajlabel);
    sub_lead_dPhi_tot[6][0][0][l] = (TH1D*) sub_lead_dPhi_tot[0][0][0][l]->Clone(tot_name);
    sub_lead_dPhi_tot[6][0][0][l]->Add( sub_lead_dPhi_tot[1][0][3][l],-1. );

   
   



    tot_name = make_name("AllHists_PbPb_pp_",2,0,0,l,pTlabel,centlabel,Ajlabel);
    signal_dPhi_tot[8][0][2][l] = (TH1D*) signal_dPhi_tot[2][0][2][l]->Clone(tot_name);
    signal_dPhi_tot[8][0][2][l]->Add( signal_dPhi_tot[3][0][3][l],-1. );

    tot_name.ReplaceAll("SubLeading","Leading");
    signal_dPhi_tot[10][0][2][l] = (TH1D*) signal_dPhi_tot[4][0][2][l]->Clone(tot_name);
    signal_dPhi_tot[10][0][2][l]->Add( signal_dPhi_tot[5][0][3][l],-1. );

    tot_name.ReplaceAll("Leading","DoubleDiff");
 
    signal_dPhi_tot[6][0][2][l] = (TH1D*) signal_dPhi_tot[8][0][2][l]->Clone(tot_name);
    signal_dPhi_tot[6][0][2][l]->Add( signal_dPhi_tot[10][0][2][l] );


    tot_name = make_name("BackgroundDiff_PbPb_pp_",2,0,3,l,pTlabel,centlabel,Ajlabel);
    background_diff_tot[7][0][2][l] = (TH1D*) background_diff_tot[0][0][2][l]->Clone(tot_name);
    background_diff_tot[7][0][2][l]->Add( background_diff_tot[1][0][3][l],-1. );

   
    tot_name = make_name("Sub_Lead_PbPb_pp_",2,0,3,l,pTlabel,centlabel,Ajlabel);
    sub_lead_dPhi_tot[6][0][2][l] = (TH1D*) sub_lead_dPhi_tot[0][0][2][l]->Clone(tot_name);
    sub_lead_dPhi_tot[6][0][2][l]->Add( sub_lead_dPhi_tot[1][0][3][l],-1. );

   




     
    for(int j = 0; j<4; j++){

      stack_name = make_name("AllHists_Up_",2,0,j,l,pTlabel,centlabel,Ajlabel);
      stack_name.ReplaceAll("TrkPt05_TrkPt1_","");

      dPhi_pTdist_up[2][j][l] = new THStack(stack_name,"");
      
      stack_name.ReplaceAll("SubLeading","Leading");
      dPhi_pTdist_up[4][j][l] = new THStack(stack_name,"");

      stack_name.ReplaceAll("Leading","Leading_PbPb_pp");
      dPhi_pTdist_up[10][j][l] = new THStack(stack_name,"");

      stack_name.ReplaceAll("Leading","SubLeading");
      dPhi_pTdist_up[8][j][l] = new THStack(stack_name,"");
      
      stack_name.ReplaceAll("SubLeading","DoubleDiff");
      dPhi_pTdist_up[6][j][l] = new THStack(stack_name,"");

      stack_name.ReplaceAll("DoubleDiff","SideBandDiff");
      dPhi_pTdist_up[0][j][l] = new THStack(stack_name,"");
      
      stack_name.ReplaceAll("SideBandDiff","SideBandDoubleDiff");
      dPhi_pTdist_up[7][j][l] = new THStack(stack_name,"");

      stack_name.ReplaceAll("SideBandDoubleDiff","SubLeadingMinusLeading");
      dPhi_Sub_Lead_NoBkg_up[0][j][l] = new THStack(stack_name,"");
    
      stack_name.ReplaceAll("SubLeadingMinusLeading","DoubleDiffNoBkgSub");
      dPhi_Sub_Lead_NoBkg_up[6][j][l] = new THStack(stack_name,"");



      stack_name = make_name("AllHists_Down_",2,0,j,l,pTlabel,centlabel,Ajlabel);
      stack_name.ReplaceAll("TrkPt05_TrkPt1_","");

      dPhi_pTdist_down[2][j][l] = new THStack(stack_name,"");
      
      stack_name.ReplaceAll("SubLeading","Leading");
      dPhi_pTdist_down[4][j][l] = new THStack(stack_name,"");

      stack_name.ReplaceAll("Leading","Leading_PbPb_pp");
      dPhi_pTdist_down[10][j][l] = new THStack(stack_name,"");

      stack_name.ReplaceAll("Leading","SubLeading");
      dPhi_pTdist_down[8][j][l] = new THStack(stack_name,"");
      
      stack_name.ReplaceAll("SubLeading","DoubleDiff");
      dPhi_pTdist_down[6][j][l] = new THStack(stack_name,"");

      stack_name.ReplaceAll("DoubleDiff","SideBandDiff");
      dPhi_pTdist_down[0][j][l] = new THStack(stack_name,"");
      
      stack_name.ReplaceAll("SideBandDiff","SideBandDoubleDiff");
      dPhi_pTdist_down[7][j][l] = new THStack(stack_name,"");

      stack_name.ReplaceAll("SideBandDoubleDiff","SubLeadingMinusLeading");
      dPhi_Sub_Lead_NoBkg_down[0][j][l] = new THStack(stack_name,"");
    
      stack_name.ReplaceAll("SubLeadingMinusLeading","DoubleDiffNoBkgSub");
      dPhi_Sub_Lead_NoBkg_down[6][j][l] = new THStack(stack_name,"");

    
      for(int i = 0; i<6; i++){
	in_name = make_name("PbPb_pp_WithErrors_",8,i,j,l,pTlabel,centlabel,Ajlabel);
	signal_dPhi_rebin[8][i][j][l] = (TH1D*)signal_dPhi_rebin[2][i][j][l]->Clone(in_name);
	signal_dPhi_rebin[8][i][j][l]->Add(signal_dPhi_rebin[3][i][3][l],-1.);

	in_name.ReplaceAll("SubLeading","Leading");
	signal_dPhi_rebin[10][i][j][l] = (TH1D*)signal_dPhi_rebin[4][i][j][l]->Clone(in_name);
	signal_dPhi_rebin[10][i][j][l]->Add(signal_dPhi_rebin[5][i][3][l],-1.);

	in_name.ReplaceAll("SubLeading","DoubleDiff");
	signal_dPhi_rebin[6][i][j][l] = (TH1D*)signal_dPhi_rebin[8][i][j][l]->Clone(in_name);
	signal_dPhi_rebin[6][i][j][l]->Add(signal_dPhi_rebin[10][i][j][l]);


	in_name.ReplaceAll("DoubleDiff","SideBandDoubleDiff");
	background_diff_rebin[7][i][j][l] = (TH1D*)background_diff_rebin[0][i][j][l]->Clone(in_name);
	background_diff_rebin[7][i][j][l]->Add(background_diff_rebin[1][i][3][l],-1.);

	in_name.ReplaceAll("SideBand","NoBkgSub");
	sub_lead_dPhi_rebin[6][i][j][l] = (TH1D*)sub_lead_dPhi_rebin[0][i][j][l]->Clone(in_name);
	sub_lead_dPhi_rebin[6][i][j][l]->Add(sub_lead_dPhi_rebin[1][i][3][l],-1.);
  

	//Now make NoErr versions


	in_name = make_name("PbPb_pp_Up_",8,i,j,l,pTlabel,centlabel,Ajlabel);
	signal_dPhi_noerr_up[8][i][j][l] = (TH1D*)signal_dPhi_rebin[8][i][j][l]->Clone(in_name);

	in_name.ReplaceAll("SubLeading","Leading");
	signal_dPhi_noerr_up[10][i][j][l] = (TH1D*)signal_dPhi_rebin[10][i][j][l]->Clone(in_name);

	in_name.ReplaceAll("SubLeading","DoubleDiff");
	signal_dPhi_noerr_up[6][i][j][l] = (TH1D*)signal_dPhi_rebin[6][i][j][l]->Clone(in_name);

	in_name.ReplaceAll("DoubleDiff","SideBandDoubleDiff");
	background_diff_noerr_up[7][i][j][l] = (TH1D*)background_diff_rebin[7][i][j][l]->Clone(in_name);

	in_name.ReplaceAll("SideBand","NoBkgSub");
	sub_lead_dPhi_noerr_up[6][i][j][l] = (TH1D*)sub_lead_dPhi_rebin[6][i][j][l]->Clone(in_name);


	in_name = make_name("PbPb_pp_Down_",8,i,j,l,pTlabel,centlabel,Ajlabel);
	signal_dPhi_noerr_down[8][i][j][l] = (TH1D*)signal_dPhi_rebin[8][i][j][l]->Clone(in_name);

	in_name.ReplaceAll("SubLeading","Leading");
	signal_dPhi_noerr_down[10][i][j][l] = (TH1D*)signal_dPhi_rebin[10][i][j][l]->Clone(in_name);

	in_name.ReplaceAll("SubLeading","DoubleDiff");
	signal_dPhi_noerr_down[6][i][j][l] = (TH1D*)signal_dPhi_rebin[6][i][j][l]->Clone(in_name);


	in_name.ReplaceAll("DoubleDiff","SideBandDoubleDiff");
	background_diff_noerr_down[7][i][j][l] = (TH1D*)background_diff_rebin[7][i][j][l]->Clone(in_name);


	in_name.ReplaceAll("SideBand","NoBkgSub");
	sub_lead_dPhi_noerr_down[6][i][j][l] = (TH1D*)sub_lead_dPhi_rebin[6][i][j][l]->Clone(in_name);
      
	for(int k = 1; k<26; k++){

	  background_diff_noerr_up[7][i][j][l]->SetBinError(k,0.);
	  background_diff_noerr_down[7][i][j][l]->SetBinError(k,0.);
	  signal_dPhi_noerr_up[8][i][j][l]->SetBinError(k,0.);
	  signal_dPhi_noerr_down[8][i][j][l]->SetBinError(k,0.);
	  signal_dPhi_noerr_up[10][i][j][l]->SetBinError(k,0.);
	  signal_dPhi_noerr_down[10][i][j][l]->SetBinError(k,0.);
	  signal_dPhi_noerr_up[6][i][j][l]->SetBinError(k,0.);
	  signal_dPhi_noerr_down[6][i][j][l]->SetBinError(k,0.);
	  sub_lead_dPhi_noerr_up[6][i][j][l]->SetBinError(k,0.);
	  sub_lead_dPhi_noerr_down[6][i][j][l]->SetBinError(k,0.);


	  if(background_diff_rebin[7][i][j][l]->GetBinContent(k)>0){
	    background_diff_noerr_up[7][i][j][l]->SetBinContent(k, background_diff_rebin[7][i][j][l]->GetBinContent(k));
	    background_diff_noerr_down[7][i][j][l]->SetBinContent(k,0.);
	  } else{
	    background_diff_noerr_down[7][i][j][l]->SetBinContent(k, background_diff_rebin[7][i][j][l]->GetBinContent(k));
	    background_diff_noerr_up[7][i][j][l]->SetBinContent(k, 0.);
	  }

	  if(	signal_dPhi_rebin[8][i][j][l]->GetBinContent(k)>0){
	    signal_dPhi_noerr_up[8][i][j][l]->SetBinContent(k,signal_dPhi_rebin[8][i][j][l]->GetBinContent(k));
	    signal_dPhi_noerr_down[8][i][j][l]->SetBinContent(k,0.);
	  }else{
	    signal_dPhi_noerr_down[8][i][j][l]->SetBinContent(k,signal_dPhi_rebin[8][i][j][l]->GetBinContent(k));
	    signal_dPhi_noerr_up[8][i][j][l]->SetBinContent(k,0.);
	  }

	  if(signal_dPhi_rebin[10][i][j][l]->GetBinContent(k)>0){
	    signal_dPhi_noerr_up[10][i][j][l]->SetBinContent(k,signal_dPhi_rebin[10][i][j][l]->GetBinContent(k));
	    signal_dPhi_noerr_down[10][i][j][l]->SetBinContent(k,0.);
	  }else{ 
	    signal_dPhi_noerr_down[10][i][j][l]->SetBinContent(k,signal_dPhi_rebin[10][i][j][l]->GetBinContent(k));
	    signal_dPhi_noerr_up[10][i][j][l]->SetBinContent(k,0.);
	  }

	  if(signal_dPhi_rebin[6][i][j][l]->GetBinContent(k)>0){
	    signal_dPhi_noerr_up[6][i][j][l]->SetBinContent(k,signal_dPhi_rebin[6][i][j][l]->GetBinContent(k));
	    signal_dPhi_noerr_down[6][i][j][l]->SetBinContent(k,0.);
	  }else{ 
	    signal_dPhi_noerr_down[6][i][j][l]->SetBinContent(k,signal_dPhi_rebin[6][i][j][l]->GetBinContent(k));
	    signal_dPhi_noerr_up[6][i][j][l]->SetBinContent(k,0.);
	  }



	  if(sub_lead_dPhi_rebin[6][i][j][l]->GetBinContent(k)>0){
	    sub_lead_dPhi_noerr_up[6][i][j][l]->SetBinContent(k,sub_lead_dPhi_rebin[6][i][j][l]->GetBinContent(k));
	    sub_lead_dPhi_noerr_down[6][i][j][l]->SetBinContent(k,0.);
	  }else{
	    sub_lead_dPhi_noerr_down[6][i][j][l]->SetBinContent(k,sub_lead_dPhi_rebin[6][i][j][l]->GetBinContent(k));
	    sub_lead_dPhi_noerr_up[6][i][j][l]->SetBinContent(k,0.);
	  }
	}


	  switch(i){
	  case 0:
	    signal_dPhi_noerr_up[8][i][j][l]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_up[10][i][j][l]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_up[6][i][j][l]->SetFillColor(kBlue-9);
	    background_diff_noerr_up[7][i][j][l]->SetFillColor(kBlue-9);
	    sub_lead_dPhi_noerr_up[6][i][j][l]->SetFillColor(kBlue-9);

	    signal_dPhi_noerr_down[8][i][j][l]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_down[10][i][j][l]->SetFillColor(kBlue-9);
	    signal_dPhi_noerr_down[6][i][j][l]->SetFillColor(kBlue-9);
	    background_diff_noerr_down[7][i][j][l]->SetFillColor(kBlue-9);
	    sub_lead_dPhi_noerr_down[6][i][j][l]->SetFillColor(kBlue-9);

	    break;
	  case 1:
	    signal_dPhi_noerr_up[8][i][j][l]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_up[10][i][j][l]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_up[6][i][j][l]->SetFillColor(kYellow-9);
	    background_diff_noerr_up[7][i][j][l]->SetFillColor(kYellow-9);
	    sub_lead_dPhi_noerr_up[6][i][j][l]->SetFillColor(kYellow-9);

	    signal_dPhi_noerr_down[8][i][j][l]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_down[10][i][j][l]->SetFillColor(kYellow-9);
	    signal_dPhi_noerr_down[6][i][j][l]->SetFillColor(kYellow-9);
	    background_diff_noerr_down[7][i][j][l]->SetFillColor(kYellow-9);
	    sub_lead_dPhi_noerr_down[6][i][j][l]->SetFillColor(kYellow-9);
	    break;
	  case 2:
	    signal_dPhi_noerr_up[8][i][j][l]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_up[10][i][j][l]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_up[6][i][j][l]->SetFillColor(kOrange+1);
	    background_diff_noerr_up[7][i][j][l]->SetFillColor(kOrange+1);
	    sub_lead_dPhi_noerr_up[6][i][j][l]->SetFillColor(kOrange+1);

	    signal_dPhi_noerr_down[8][i][j][l]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_down[10][i][j][l]->SetFillColor(kOrange+1);
	    signal_dPhi_noerr_down[6][i][j][l]->SetFillColor(kOrange+1);
	    background_diff_noerr_down[7][i][j][l]->SetFillColor(kOrange+1);
	    sub_lead_dPhi_noerr_down[6][i][j][l]->SetFillColor(kOrange+1);
	    break;
	  case 3:
	    signal_dPhi_noerr_up[8][i][j][l]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_up[10][i][j][l]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_up[6][i][j][l]->SetFillColor(kViolet-5);
	    background_diff_noerr_up[7][i][j][l]->SetFillColor(kViolet-5);
	    sub_lead_dPhi_noerr_up[6][i][j][l]->SetFillColor(kViolet-5);

	    signal_dPhi_noerr_down[8][i][j][l]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_down[10][i][j][l]->SetFillColor(kViolet-5);
	    signal_dPhi_noerr_down[6][i][j][l]->SetFillColor(kViolet-5);
	    background_diff_noerr_down[7][i][j][l]->SetFillColor(kViolet-5);
	    sub_lead_dPhi_noerr_down[6][i][j][l]->SetFillColor(kViolet-5);
	    break;
	  case 4:
	    signal_dPhi_noerr_up[8][i][j][l]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_up[10][i][j][l]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_up[6][i][j][l]->SetFillColor(kGreen+3);
	    background_diff_noerr_up[7][i][j][l]->SetFillColor(kGreen+3);
	    sub_lead_dPhi_noerr_up[6][i][j][l]->SetFillColor(kGreen+3);

	    signal_dPhi_noerr_down[8][i][j][l]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_down[10][i][j][l]->SetFillColor(kGreen+3);
	    signal_dPhi_noerr_down[6][i][j][l]->SetFillColor(kGreen+3);
	    background_diff_noerr_down[7][i][j][l]->SetFillColor(kGreen+3);
	    sub_lead_dPhi_noerr_down[6][i][j][l]->SetFillColor(kGreen+3);
	    break;
	  case 5:
	    signal_dPhi_noerr_up[8][i][j][l]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_up[10][i][j][l]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_up[6][i][j][l]->SetFillColor(kRed+1);
	    background_diff_noerr_up[7][i][j][l]->SetFillColor(kRed+1);
	    sub_lead_dPhi_noerr_up[6][i][j][l]->SetFillColor(kRed+1);

	    signal_dPhi_noerr_down[8][i][j][l]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_down[10][i][j][l]->SetFillColor(kRed+1);
	    signal_dPhi_noerr_down[6][i][j][l]->SetFillColor(kRed+1);
	    background_diff_noerr_down[7][i][j][l]->SetFillColor(kRed+1);
	    sub_lead_dPhi_noerr_down[6][i][j][l]->SetFillColor(kRed+1);
	    break;
	  default:
	    break;
	  }

	  signal_dPhi_noerr_up[8][i][j][l]->SetFillStyle(1001);
	  signal_dPhi_noerr_up[10][i][j][l]->SetFillStyle(1001);
	  signal_dPhi_noerr_up[6][i][j][l]->SetFillStyle(1001);
	  background_diff_noerr_up[7][i][j][l]->SetFillStyle(1001);
	  sub_lead_dPhi_noerr_up[6][i][j][l]->SetFillStyle(1001);
	  
	  signal_dPhi_noerr_down[8][i][j][l]->SetFillStyle(1001);
	  signal_dPhi_noerr_down[10][i][j][l]->SetFillStyle(1001);
	  signal_dPhi_noerr_down[6][i][j][l]->SetFillStyle(1001);
	  background_diff_noerr_down[7][i][j][l]->SetFillStyle(1001);
	  sub_lead_dPhi_noerr_down[6][i][j][l]->SetFillStyle(1001);
	  
      }


      for(int i = 1; i<6; i++){

	dPhi_pTdist_up[2][j][l]->Add(signal_dPhi_noerr_up[2][5-i][j][l]);
	dPhi_pTdist_up[4][j][l]->Add(signal_dPhi_noerr_up[4][5-i][j][l]);
	dPhi_pTdist_up[8][j][l]->Add(signal_dPhi_noerr_up[8][5-i][j][l]);
	dPhi_pTdist_up[10][j][l]->Add(signal_dPhi_noerr_up[10][5-i][j][l]);  
	dPhi_pTdist_up[6][j][l]->Add(signal_dPhi_noerr_up[6][5-i][j][l]);
	dPhi_pTdist_up[0][j][l]->Add(background_diff_noerr_up[0][5-i][j][l]);
	dPhi_pTdist_up[7][j][l]->Add(background_diff_noerr_up[7][5-i][j][l]);
	dPhi_Sub_Lead_NoBkg_up[0][j][l]->Add(sub_lead_dPhi_noerr_up[0][5-i][j][l]);
	dPhi_Sub_Lead_NoBkg_up[6][j][l]->Add(sub_lead_dPhi_noerr_up[6][5-i][j][l]);

	dPhi_pTdist_down[2][j][l]->Add(signal_dPhi_noerr_down[2][5-i][j][l]);
	dPhi_pTdist_down[4][j][l]->Add(signal_dPhi_noerr_down[4][5-i][j][l]);
	dPhi_pTdist_down[8][j][l]->Add(signal_dPhi_noerr_down[8][5-i][j][l]);
	dPhi_pTdist_down[10][j][l]->Add(signal_dPhi_noerr_down[10][5-i][j][l]);  
	dPhi_pTdist_down[6][j][l]->Add(signal_dPhi_noerr_down[6][5-i][j][l]);
	dPhi_pTdist_down[0][j][l]->Add(background_diff_noerr_down[0][5-i][j][l]);
	dPhi_pTdist_down[7][j][l]->Add(background_diff_noerr_down[7][5-i][j][l]);
	dPhi_Sub_Lead_NoBkg_down[0][j][l]->Add(sub_lead_dPhi_noerr_down[0][5-i][j][l]);
	dPhi_Sub_Lead_NoBkg_down[6][j][l]->Add(sub_lead_dPhi_noerr_down[6][5-i][j][l]);


      }
      signal_max = 104.;
      signal_min = -106.;
      

      diff_max = 16.;
      diff_min = -16.;
      
      double_diff_max = 11.;
      double_diff_min = -11.; 
      
      no_bgsub_max = 25.;
      no_bgsub_min = -25.;


      dPhi_pTdist_up[2][j][l]->SetMaximum(signal_max);
      dPhi_pTdist_up[2][j][l]->SetMinimum(signal_min);


      dPhi_pTdist_up[8][j][l]->SetMaximum(diff_max);
      dPhi_pTdist_up[8][j][l]->SetMinimum(diff_min);


      dPhi_pTdist_up[10][j][l]->SetMaximum(diff_max);
      dPhi_pTdist_up[10][j][l]->SetMinimum(diff_min);

      dPhi_Sub_Lead_NoBkg_up[0][j][l]->SetMaximum(no_bgsub_max);
      dPhi_Sub_Lead_NoBkg_up[0][j][l]->SetMinimum(no_bgsub_min);
  
      dPhi_pTdist_up[6][j][l]->SetMaximum(double_diff_max);
      dPhi_pTdist_up[6][j][l]->SetMinimum(double_diff_min);

   
      dPhi_pTdist_up[0][j][l]->SetMaximum(diff_max);
      dPhi_pTdist_up[0][j][l]->SetMinimum(diff_min);
     

      dPhi_pTdist_up[7][j][l]->SetMaximum(double_diff_max);
      dPhi_pTdist_up[7][j][l]->SetMinimum(double_diff_min);

      dPhi_Sub_Lead_NoBkg_up[6][j][l]->SetMaximum(double_diff_max);
      dPhi_Sub_Lead_NoBkg_up[6][j][l]->SetMinimum(double_diff_min);



      dPhi_pTdist_up[2][j][l]->SetMaximum(signal_max);
      dPhi_pTdist_up[2][j][l]->SetMinimum(signal_min);


      dPhi_pTdist_up[8][j][l]->SetMaximum(diff_max);
      dPhi_pTdist_up[8][j][l]->SetMinimum(diff_min);


      dPhi_pTdist_down[10][j][l]->SetMaximum(diff_max);
      dPhi_pTdist_down[10][j][l]->SetMinimum(diff_min);

      dPhi_Sub_Lead_NoBkg_down[0][j][l]->SetMaximum(no_bgsub_max);
      dPhi_Sub_Lead_NoBkg_down[0][j][l]->SetMinimum(no_bgsub_min);
  
      dPhi_pTdist_down[6][j][l]->SetMaximum(double_diff_max);
      dPhi_pTdist_down[6][j][l]->SetMinimum(double_diff_min);

   
      dPhi_pTdist_down[0][j][l]->SetMaximum(diff_max);
      dPhi_pTdist_down[0][j][l]->SetMinimum(diff_min);
     

      dPhi_pTdist_down[7][j][l]->SetMaximum(diff_max);
      dPhi_pTdist_down[7][j][l]->SetMinimum(diff_min);

      dPhi_Sub_Lead_NoBkg_down[6][j][l]->SetMaximum(no_bgsub_min);
      dPhi_Sub_Lead_NoBkg_down[6][j][l]->SetMinimum(no_bgsub_min);
  
  
      
    }
   

    stack_name = make_name("AllHists_Up_",5,0,3,l,pTlabel,centlabel,Ajlabel);
    dPhi_pTdist_up[5][3][l] = new THStack(stack_name,"");
      
    stack_name.ReplaceAll("Leading","Subleading");
    dPhi_pTdist_up[3][3][l] = new THStack(stack_name,"");

    stack_name.ReplaceAll("Subleading","SideBandDiff");
    dPhi_pTdist_up[1][3][l] = new THStack(stack_name,"");

    stack_name.ReplaceAll("SideBandDiff","SubLeadingMinusLeading");
    dPhi_Sub_Lead_NoBkg_up[1][3][l] = new THStack(stack_name,"");
 

    stack_name = make_name("AllHists_Down_",5,0,3,l,pTlabel,centlabel,Ajlabel);
    dPhi_pTdist_down[5][3][l] = new THStack(stack_name,"");
      
    stack_name.ReplaceAll("Leading","SubLeading");
    dPhi_pTdist_down[3][3][l] = new THStack(stack_name,"");

    stack_name.ReplaceAll("SubLeading","SideBandDiff");
    dPhi_pTdist_down[1][3][l] = new THStack(stack_name,"");

    stack_name.ReplaceAll("SideBandDiff","SubLeadingMinusLeading");
    dPhi_Sub_Lead_NoBkg_down[1][3][l] = new THStack(stack_name,"");
 

    for(int i = 1; i<6; i++){

      dPhi_pTdist_up[3][3][l]->Add(signal_dPhi_noerr_up[3][5-i][3][l]);
      dPhi_pTdist_up[5][3][l]->Add(signal_dPhi_noerr_up[5][5-i][3][l]);
      dPhi_pTdist_up[1][3][l]->Add(background_diff_noerr_up[1][5-i][3][l]);
      dPhi_Sub_Lead_NoBkg_up[1][3][l]->Add(sub_lead_dPhi_noerr_up[1][5-i][3][l]);

      dPhi_pTdist_down[3][3][l]->Add(signal_dPhi_noerr_down[3][5-i][3][l]);
      dPhi_pTdist_down[5][3][l]->Add(signal_dPhi_noerr_down[5][5-i][3][l]);
      dPhi_pTdist_down[1][3][l]->Add(background_diff_noerr_down[1][5-i][3][l]);
      dPhi_Sub_Lead_NoBkg_down[1][3][l]->Add(sub_lead_dPhi_noerr_down[1][5-i][3][l]);
  
    }
    dPhi_pTdist_up[3][3][l]->SetMaximum(signal_max);
    dPhi_pTdist_up[3][3][l]->SetMinimum(signal_min);

    dPhi_Sub_Lead_NoBkg_up[1][3][l]->SetMinimum(no_bgsub_min);
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->SetMaximum(no_bgsub_max);

    dPhi_pTdist_up[1][3][l]->SetMaximum(diff_max);
    dPhi_pTdist_up[1][3][l]->SetMinimum(diff_min);

    dPhi_pTdist_down[3][3][l]->SetMaximum(signal_max);
    dPhi_pTdist_down[3][3][l]->SetMinimum(signal_min);

    dPhi_Sub_Lead_NoBkg_down[1][3][l]->SetMinimum(diff_min);
    dPhi_Sub_Lead_NoBkg_down[1][3][l]->SetMaximum(diff_max);

    dPhi_pTdist_down[1][3][l]->SetMaximum(diff_max);
    dPhi_pTdist_down[1][3][l]->SetMinimum(diff_min);

 

    ///////////////////////////////////////////////////////////////

    //     DRAWING STARTS HERE  

    ///////////////////////////////////////////////////////////////
    
    c_jet[l]->cd(0);

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
   
    TLatex   *luminosity_tex_pp = new TLatex(0.07,0.96,"pp 5.3 pb^{-1} (2.76 TeV)");
    luminosity_tex_pp->SetTextFont(43);
    luminosity_tex_pp->SetTextSizePixels(25);
    luminosity_tex_pp->SetLineColor(kWhite);
    luminosity_tex_pp->SetNDC();
    luminosity_tex_pp->Draw();
   
    TLatex   *luminosity_tex_PbPb = new TLatex(0.3,0.96,"PbPb 166 #mub^{-1} (2.76 TeV)");
    luminosity_tex_PbPb->SetTextFont(43);
    luminosity_tex_PbPb->SetTextSizePixels(25);
    luminosity_tex_PbPb->SetLineColor(kWhite);
    luminosity_tex_PbPb->SetNDC();
    luminosity_tex_PbPb->Draw();
    

    TLatex   *jet_reco_tex = new TLatex(0.605,0.98,"anti-k_{T} R = 0.3, |#eta_{jet}| < 1.6");
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
    


   

    c_jet[l]->cd(4);

    signal_dPhi_syst[2][0][0][l]->SetLineColor(kBlack);

    TLegend *legend = new TLegend(0.1,0.15,0.8,.85);
    legend->AddEntry(signal_dPhi_noerr_up[2][0][0][l],"0.5 < p_{T}^{assoc.}< 1 GeV/c","f");
    legend->AddEntry(signal_dPhi_noerr_up[2][1][0][l],"1 < p_{T}^{assoc.}< 2 GeV/c","f");
    legend->AddEntry(signal_dPhi_noerr_up[2][2][0][l],"2 < p_{T}^{assoc.}< 3 GeV/c","f");
    legend->AddEntry(signal_dPhi_noerr_up[2][3][0][l],"3 < p_{T}^{assoc.}< 4 GeV/c","f");
    legend->AddEntry(signal_dPhi_noerr_up[2][4][0][l],"4 < p_{T}^{assoc.}< 8 GeV/c","f");
    legend->AddEntry(signal_dPhi_syst[2][0][0][l],"0.5 < p_{T}^{assoc.}< 8 GeV/c","lpfe");
    legend->AddEntry(signal_dPhi_syst[2][0][0][l],"|#eta_{track}| < 2.4 ","");
    legend->SetTextSize(0.06);
    legend->SetLineColor(kWhite);
    legend->Draw();

  
 
  
    c_jet[l]->cd(1);
 
    dPhi_pTdist_up[3][3][l]->Draw();
    dPhi_pTdist_up[3][3][l]->GetYaxis()->SetLabelSize(0.08); 
    dPhi_pTdist_up[3][3][l]->GetXaxis()->SetLabelSize(0.08); 
    dPhi_pTdist_up[3][3][l]->GetXaxis()->SetTitle("#Delta#phi");
    dPhi_pTdist_up[3][3][l]->GetXaxis()->SetTitleOffset(0.7);
    dPhi_pTdist_up[3][3][l]->GetYaxis()->SetTitleOffset(.6);
    dPhi_pTdist_up[3][3][l]->GetYaxis()->SetNdivisions(510);
    dPhi_pTdist_up[3][3][l]->GetXaxis()->SetTitleSize(0.09);
    dPhi_pTdist_up[3][3][l]->GetXaxis()->CenterTitle();
    dPhi_pTdist_up[3][3][l]->GetYaxis()->SetTitle("P  (GeV/c)");
    dPhi_pTdist_up[3][3][l]->GetYaxis()->SetTitleSize(0.09);
    dPhi_pTdist_up[3][3][l]->GetYaxis()->CenterTitle();
    dPhi_pTdist_up[5][3][l]->Draw("same");
    dPhi_pTdist_down[5][3][l]->Draw("same");
    dPhi_pTdist_down[3][3][l]->Draw("same");
  
   
    signal_dPhi_syst[3][0][3][l]->Draw("same e2");
    signal_dPhi_syst[5][0][3][l]->Draw("same e2");
  
    signal_dPhi_tot[3][0][3][l]->Draw("same");
    signal_dPhi_tot[5][0][3][l]->Draw("same");
 
 
    label_pp = new TLatex(0.2,0.9,"pp Reference");
    label_pp->SetTextFont(63);
    label_pp->SetTextSizePixels(30);
    label_pp->SetLineColor(kWhite);
    label_pp->SetNDC();
    label_pp->Draw();

    TLine *l_phi = new TLine(-TMath::Pi()/2.,0.,TMath::Pi()/2.,0.);
    l_phi->SetLineStyle(2);
    l_phi->Draw();


    TLatex *cms_tex_dphi = new TLatex(0.2,0.05,"CMS");
    cms_tex_dphi->SetTextFont(63);
    cms_tex_dphi->SetTextSizePixels(30);
    cms_tex_dphi->SetLineColor(kWhite);
    cms_tex_dphi->SetNDC();
    cms_tex_dphi->Draw(); 

    TLatex *prelim_tex_dphi = new TLatex(0.38,0.05,"Preliminary");
    prelim_tex_dphi->SetTextFont(53);
    prelim_tex_dphi->SetTextSizePixels(30);
    prelim_tex_dphi->SetLineColor(kWhite);
    prelim_tex_dphi->SetNDC();
    prelim_tex_dphi->Draw(); 
  
    
    gPad->RedrawAxis();

   c_jet[l]->cd(2);

   dPhi_pTdist_up[2][0][l]->Draw();
   dPhi_pTdist_up[2][0][l]->GetYaxis()->SetLabelSize(0.);
   dPhi_pTdist_up[2][0][l]->GetYaxis()->SetNdivisions(306);
   dPhi_pTdist_up[4][0][l]->Draw("same");
   dPhi_pTdist_down[4][0][l]->Draw("same");
   dPhi_pTdist_down[2][0][l]->Draw("same");


   signal_dPhi_syst[2][0][0][l]->Draw("same e2");
   signal_dPhi_syst[4][0][0][l]->Draw("same e2");
 


   signal_dPhi_tot[2][0][0][l]->Draw("same");
   signal_dPhi_tot[4][0][0][l]->Draw("same");
 
   label_per = new TLatex(0.05,0.9,"PbPb Cent. 50-100%");
   label_per->SetTextFont(63);
   label_per->SetTextSizePixels(30);
   label_per->SetLineColor(kWhite);
   label_per->SetNDC();
   label_per->Draw();
   l_phi->Draw();

   c_jet[l]->cd(3);
   
   dPhi_pTdist_up[2][2][l]->Draw();
   dPhi_pTdist_up[2][2][l]->GetYaxis()->SetLabelSize(0.);
   dPhi_pTdist_up[2][2][l]->GetYaxis()->SetNdivisions(306);
   dPhi_pTdist_up[4][2][l]->Draw("same");
   dPhi_pTdist_down[4][2][l]->Draw("same");
   dPhi_pTdist_down[2][2][l]->Draw("same");

 signal_dPhi_syst[2][0][2][l]->Draw("same e2");
   signal_dPhi_syst[4][0][2][l]->Draw("same e2");
 


   signal_dPhi_tot[2][0][2][l]->Draw("same");
   signal_dPhi_tot[4][0][2][l]->Draw("same");
 
  
   label_cent = new TLatex(0.05,0.9,"PbPb Cent. 0-30%");
   label_cent->SetTextFont(63);
   label_cent->SetTextSizePixels(30);
   label_cent->SetLineColor(kWhite);
   label_cent->SetNDC();
   label_cent->Draw();

   TLatex *orientation = new TLatex(.97,0.,"      <-Leading          Subleading->");
   orientation->SetTextSize(0.065);
   orientation->SetLineColor(kWhite);
   orientation->SetNDC();
   orientation->SetTextAngle(90.);
   orientation->Draw();
 
   l_phi->Draw();


   c_jet[l]->cd(5);
   dPhi_pTdist_up[8][0][l]->Draw();
 dPhi_pTdist_up[8][0][l]->GetYaxis()->SetLabelSize(0.0); 
   dPhi_pTdist_up[8][0][l]->GetXaxis()->CenterTitle();
   dPhi_pTdist_up[8][0][l]->GetYaxis()->SetTitle("P_{PbPb}- P_{pp}  (GeV/c)");
   dPhi_pTdist_up[8][0][l]->GetYaxis()->SetTitleSize(0.0);
   dPhi_pTdist_up[8][0][l]->GetYaxis()->SetTitleOffset(0.8);
   dPhi_pTdist_up[8][0][l]->GetYaxis()->CenterTitle();
   dPhi_pTdist_down[8][0][l]->Draw("same");

   signal_dPhi_syst[8][0][0][l]->Draw("same e2");
   signal_dPhi_tot[8][0][0][l]->Draw("same");
 
    

   TLatex *label = new TLatex(0.05,0.9,"Subleading PbPb - pp");
   label->SetTextFont(63);
   label->SetTextSizePixels(30);
   label->SetLineColor(kWhite);
   label->SetNDC();
   label->Draw();

   l_phi->Draw();

   gPad->RedrawAxis();

   c_jet[l]->cd(6);

   dPhi_pTdist_up[8][2][l]->Draw();
   dPhi_pTdist_up[8][2][l]->GetYaxis()->SetLabelSize(0.);
   l_phi->Draw();
   dPhi_pTdist_down[8][2][l]->Draw("same");

   signal_dPhi_syst[8][0][2][l]->Draw("same e2");
   signal_dPhi_tot[8][0][2][l]->Draw("same");
  
   
   c_jet[l]->cd(8);
   dPhi_pTdist_up[10][0][l]->Draw();
   dPhi_pTdist_up[10][0][l]->GetYaxis()->SetLabelSize(0.0); 
   dPhi_pTdist_up[10][0][l]->GetXaxis()->CenterTitle();
   dPhi_pTdist_up[10][0][l]->GetYaxis()->SetTitle("P_{PbPb}- P_{pp}  (GeV/c)");
   dPhi_pTdist_up[10][0][l]->GetYaxis()->SetTitleSize(0.0);
   dPhi_pTdist_up[10][0][l]->GetYaxis()->SetTitleOffset(0.8);
   dPhi_pTdist_up[10][0][l]->GetYaxis()->CenterTitle();
   dPhi_pTdist_down[10][0][l]->Draw("same");

   signal_dPhi_syst[10][0][0][l]->Draw("same e2");
   signal_dPhi_tot[10][0][0][l]->Draw("same");
 
  

   label = new TLatex(0.05,0.9,"Leading PbPb - pp (#times -1)");
     label->SetTextFont(63);
   label->SetTextSizePixels(30);
   label->SetLineColor(kWhite);
   label->SetNDC();
   label->Draw();

   l_phi->Draw();

   gPad->RedrawAxis();

  
   c_jet[l]->cd(9);
   dPhi_pTdist_up[10][2][l]->Draw();
   dPhi_pTdist_up[10][2][l]->GetYaxis()->SetLabelSize(0.);
   dPhi_pTdist_down[10][2][l]->Draw("same");
 

   signal_dPhi_syst[10][0][2][l]->Draw("same e2");
   signal_dPhi_tot[10][0][2][l]->Draw("same");
 
 

   l_phi->Draw();

   c_jet[l]->cd(11);
   dPhi_pTdist_up[6][0][l]->Draw();
   dPhi_pTdist_up[6][0][l]->GetYaxis()->SetLabelSize(0.0); 
 
   dPhi_pTdist_up[6][0][l]->GetYaxis()->SetTitle("#DeltaP_{PbPb}- #DeltaP_{pp}  (GeV/c)");
   dPhi_pTdist_up[6][0][l]->GetYaxis()->SetTitleSize(0.0);
   dPhi_pTdist_up[6][0][l]->GetYaxis()->SetTitleOffset(0.8);
   dPhi_pTdist_up[6][0][l]->GetYaxis()->CenterTitle();

   dPhi_pTdist_up[6][0][l]->GetXaxis()->SetTitle("#Delta#phi");
   dPhi_pTdist_up[6][0][l]->GetXaxis()->SetTitleOffset(0.7);
   dPhi_pTdist_up[6][0][l]->GetXaxis()->SetTitleSize(0.09);
   dPhi_pTdist_up[6][0][l]->GetXaxis()->CenterTitle();
   dPhi_pTdist_up[6][0][l]->GetXaxis()->SetLabelSize(0.08); 
   dPhi_pTdist_down[6][0][l]->Draw("same");

   signal_dPhi_syst[6][0][0][l]->Draw("same e2");
   signal_dPhi_tot[6][0][0][l]->Draw("same");
 
 
   label = new TLatex(0.05,0.9,"Subleading - Leading");
     label->SetTextFont(63);
   label->SetTextSizePixels(30);
   label->SetLineColor(kWhite);
   label->SetNDC();
   label->Draw();
   TLatex *label2 = new TLatex(0.05,0.8,"(PbPb - pp)");
   label2->SetTextFont(63);
   label2->SetTextSizePixels(30);
   label2->SetLineColor(kWhite);
   label2->SetNDC();
   label2->Draw();

   l_phi->Draw();

   TPave *cover_x_l = new TPave(-.08,0.,0.1,0.13);
   cover_x_l->SetFillColor(kWhite);
   cover_x_l->SetLineColor(kWhite);
   cover_x_l->SetOption("NDC NB");
   cover_x_l->Draw();

   TPave *cover_x_r = new TPave(0.83,0.,1.05,0.13);
   cover_x_r->SetFillColor(kWhite);
   cover_x_r->SetLineColor(kWhite);
   cover_x_r->SetOption("NDC NB");
   cover_x_r->Draw();

 
   c_jet[l]->cd(12);
   dPhi_pTdist_up[6][2][l]->Draw();
   dPhi_pTdist_up[6][2][l]->GetYaxis()->SetLabelSize(0.);

   dPhi_pTdist_up[6][2][l]->GetXaxis()->SetTitle("#Delta#phi");
   dPhi_pTdist_up[6][2][l]->GetXaxis()->SetTitleOffset(0.7);
   dPhi_pTdist_up[6][2][l]->GetXaxis()->SetTitleSize(0.09);
   dPhi_pTdist_up[6][2][l]->GetXaxis()->CenterTitle();
   dPhi_pTdist_up[6][2][l]->GetXaxis()->SetLabelSize(0.08); 
   dPhi_pTdist_down[6][2][l]->Draw("same");

   signal_dPhi_syst[6][0][2][l]->Draw("same e2");
   signal_dPhi_tot[6][0][2][l]->Draw("same");

 

   l_phi->Draw();
   //  orientation->Draw();

   cover_x_l->Draw();
   cover_x_r->Draw();


   c_jet[l]->cd(4);
   
   TGaxis *dummy_axis_diff = new TGaxis(1.,0.,1.0,.975,diff_min,diff_max);

   dummy_axis_diff->ImportAxisAttributes( dPhi_pTdist_up[3][3][l]->GetYaxis());
   dummy_axis_diff->SetTitleOffset(.9);
   dummy_axis_diff->SetTitleSize(0.09);
   dummy_axis_diff->SetLabelSize(0.08);
   dummy_axis_diff->SetTickSize(0.);
   dummy_axis_diff->SetTitle("P_{PbPb}- P_{pp}  (GeV/c)");
   dummy_axis_diff->Draw();

  
   c_jet[l]->cd(7);
  dummy_axis_diff->Draw();

  c_jet[l]->cd(10);
  TGaxis *dummy_axis_double_diff = new TGaxis(1.,0.13,1.0,.975,double_diff_min,double_diff_max);

  dummy_axis_double_diff->ImportAxisAttributes( dPhi_pTdist_up[3][3][l]->GetYaxis());

  dummy_axis_double_diff->SetTitleOffset(0.9);
  dummy_axis_double_diff->SetTitleSize(0.085);
  dummy_axis_double_diff->SetLabelSize(0.075);
  dummy_axis_double_diff->SetTickSize(0.);
  dummy_axis_double_diff->SetTitle("#DeltaP_{PbPb}- #DeltaP_{pp}  (GeV/c)   ");
  dummy_axis_double_diff->Draw();
 

  c_jet[l]->cd(4);
  TGaxis *dummy_axis_phi = new TGaxis(0.165,1.,1.0,1.,-TMath::Pi()/2.,TMath::Pi()/2.);

  dummy_axis_phi->ImportAxisAttributes( dPhi_pTdist_up[3][3][l]->GetXaxis());
  dummy_axis_phi->SetTitleOffset(0.8);
  dummy_axis_phi->SetTickSize(0.);
  dummy_axis_phi->Draw();
 
  TPave *cover_x_l2 = new TPave(0.9,0.9,1.05,0.99);
  cover_x_l2->SetFillColor(kWhite);
  cover_x_l2->SetLineColor(kWhite);
  cover_x_l2->SetOption("NDC NB");
  cover_x_l2->Draw();

     
 
  TPave *cover_x_r2 = new TPave(-.05,0.9,0.5,0.99);
  cover_x_r2->SetFillColor(kWhite);
  cover_x_r2->SetLineColor(kWhite);
  cover_x_r2->SetOption("NDC NB");
  cover_x_r2->Draw();

     



   c_longrange[l]->cd(0);

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
   
   luminosity_tex_pp = new TLatex(0.07,0.915,"pp 5.3 pb^{-1} (2.76 TeV)");
   luminosity_tex_pp->SetTextFont(43);
   luminosity_tex_pp->SetTextSizePixels(25);
   luminosity_tex_pp->SetLineColor(kWhite);
   luminosity_tex_pp->SetNDC();
   luminosity_tex_pp->Draw();
   
   luminosity_tex_PbPb = new TLatex(0.3,0.915,"PbPb 166 #mub^{-1} (2.76 TeV)");
   luminosity_tex_PbPb->SetTextFont(43);
   luminosity_tex_PbPb->SetTextSizePixels(25);
   luminosity_tex_PbPb->SetLineColor(kWhite);
   luminosity_tex_PbPb->SetNDC();
   luminosity_tex_PbPb->Draw();
    

   jet_reco_tex = new TLatex(0.605,0.96,"anti-k_{T} R = 0.3, |#eta_{jet}| < 1.6");
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

   
   c_longrange[l]->cd(4);
   legend->Draw();
  
   dummy_axis_phi->Draw();
   cover_x_l2->Draw();

    
   c_longrange[l]->cd(1);
   dPhi_pTdist_up[1][3][l]->Draw();
   label_pp->Draw();
 
   dPhi_pTdist_up[1][3][l]->GetYaxis()->SetTitle("#DeltaP  (GeV/c)");
   dPhi_pTdist_up[1][3][l]->GetYaxis()->SetTitleSize(0.09);
   dPhi_pTdist_up[1][3][l]->GetYaxis()->CenterTitle();
   dPhi_pTdist_up[1][3][l]->GetYaxis()->SetLabelSize(0.08); 
   dPhi_pTdist_up[1][3][l]->GetYaxis()->SetTitleOffset(0.7);

   dPhi_pTdist_up[1][3][l]->GetXaxis()->SetTitle("#Delta#phi");
   dPhi_pTdist_up[1][3][l]->GetXaxis()->SetTitleOffset(0.7);
   dPhi_pTdist_up[1][3][l]->GetXaxis()->SetTitleSize(0.09);
   dPhi_pTdist_up[1][3][l]->GetXaxis()->CenterTitle();
   dPhi_pTdist_up[1][3][l]->GetXaxis()->SetLabelSize(0.08); 
   dPhi_pTdist_down[1][3][l]->Draw("same");

   background_syst_rebin[1][4][3][l]->Draw("same e2");
   background_diff_tot[1][0][3][l]->Draw("same");
  
  l_phi->Draw();

  cms_tex_dphi->Draw(); 
  prelim_tex_dphi->Draw(); 
   c_longrange[l]->cd(2);

   dPhi_pTdist_up[0][0][l]->Draw();
   dPhi_pTdist_up[0][0][l]->GetYaxis()->SetLabelSize(0.);
   dPhi_pTdist_down[0][0][l]->Draw("same");


   background_syst_rebin[0][4][0][l]->Draw("same e2");
   background_diff_tot[0][0][0][l]->Draw("same");
 
   label_per->Draw();
   l_phi->Draw();

   c_longrange[l]->cd(3);
 
   dPhi_pTdist_up[0][2][l]->Draw();
   dPhi_pTdist_up[0][2][l]->GetYaxis()->SetLabelSize(0.);
   dPhi_pTdist_down[0][2][l]->Draw("same");
   label_cent->Draw();
   orientation->Draw();
   
   background_syst_rebin[0][4][2][l]->Draw("same e2");
   background_diff_tot[0][0][2][l]->Draw("same");

   l_phi->Draw();

   c_longrange[l]->cd(5);
   dPhi_pTdist_up[7][0][l]->Draw();
  

   dPhi_pTdist_up[7][0][l]->GetXaxis()->SetTitle("#Delta#phi");
   dPhi_pTdist_up[7][0][l]->GetXaxis()->SetTitleOffset(0.7);
   dPhi_pTdist_up[7][0][l]->GetXaxis()->SetTitleSize(0.09);
   dPhi_pTdist_up[7][0][l]->GetXaxis()->CenterTitle();
   dPhi_pTdist_up[7][0][l]->GetXaxis()->SetLabelSize(0.08); 
   
   dPhi_pTdist_up[7][0][l]->GetYaxis()->SetTitle("#DeltaP_{PbPb}- #DeltaP_{pp}  (GeV/c)");
   dPhi_pTdist_up[7][0][l]->GetYaxis()->SetTitleSize(0.0);
   dPhi_pTdist_up[7][0][l]->GetYaxis()->CenterTitle();
   dPhi_pTdist_up[7][0][l]->GetYaxis()->SetLabelSize(0.0); 
   dPhi_pTdist_up[7][0][l]->GetYaxis()->SetTitleOffset(0.5);
   l_phi->Draw();
   dPhi_pTdist_down[7][0][l]->Draw("same");

   background_syst_rebin[7][4][0][l]->Draw("same e2");
   background_diff_tot[7][0][0][l]->Draw("same");

  TLatex *label_per2 = new TLatex(0.05,0.9,"PbPb - pp");
  label_per2->SetTextFont(63);
  label_per2->SetTextSizePixels(30);
   label_per2->SetLineColor(kWhite);
   label_per2->SetNDC();
   label_per2->Draw();
  
   cover_x_r = new TPave(-0.05,0.,.1,0.13);
   cover_x_r->SetFillColor(kWhite);
   cover_x_r->SetLineColor(kWhite);
   cover_x_r->SetOption("NDC NB");
   cover_x_r->Draw();
   cover_x_r->Draw();

   cover_x_l = new TPave(0.85,0.,1.05,0.13);
   cover_x_l->SetFillColor(kWhite);
   cover_x_l->SetLineColor(kWhite);
   cover_x_l->SetOption("NDC NB");
   cover_x_l->Draw();
   cover_x_l->Draw();


   c_longrange[l]->cd(6);
   dPhi_pTdist_up[7][2][l]->Draw();
   dPhi_pTdist_up[7][2][l]->GetYaxis()->SetLabelSize(0.);
   

   dPhi_pTdist_up[7][2][l]->GetXaxis()->SetTitle("#Delta#phi");
   dPhi_pTdist_up[7][2][l]->GetXaxis()->SetTitleOffset(0.7);
   dPhi_pTdist_up[7][2][l]->GetXaxis()->SetTitleSize(0.09);
   dPhi_pTdist_up[7][2][l]->GetXaxis()->CenterTitle();
   dPhi_pTdist_up[7][2][l]->GetXaxis()->SetLabelSize(0.08); 
   dPhi_pTdist_down[7][2][l]->Draw("same");

   background_syst_rebin[7][4][2][l]->Draw("same e2");
   background_diff_tot[7][0][2][l]->Draw("same");
   l_phi->Draw();

   cover_x_r->Draw();
   cover_x_l->Draw();

   c_longrange[l]->cd(4);
   dummy_axis_double_diff->Draw();
   
   c_nobkgsub[l]->cd(0);
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
    


   c_nobkgsub[l]->cd(4);
   legend->Draw();
   dummy_axis_phi->Draw();
   cover_x_l2->Draw();
   

  
   
    c_nobkgsub[l]->cd(1);
 
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->Draw();
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetYaxis()->SetLabelSize(0.08); 
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetXaxis()->SetLabelSize(0.08); 
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetXaxis()->SetTitle("#Delta#phi");
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetXaxis()->SetTitleOffset(0.7);
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetYaxis()->SetTitleOffset(0.7);
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetXaxis()->SetTitleSize(0.09);
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetXaxis()->CenterTitle();
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetYaxis()->SetTitle("#DeltaP  (GeV/c)");
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetYaxis()->SetTitleSize(0.09);
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetYaxis()->CenterTitle();
    dPhi_Sub_Lead_NoBkg_up[1][3][l]->GetYaxis()->SetNdivisions(408);
    dPhi_Sub_Lead_NoBkg_down[1][3][l]->Draw("same");   


    sub_lead_dPhi_syst[3][0][3][l]->Draw("same e2");
    sub_lead_dPhi_tot[1][0][3][l]->Draw("same");

    label_pp->Draw();

    l_phi->Draw();

     cms_tex_dphi->Draw(); 
  prelim_tex_dphi->Draw();  
    c_nobkgsub[l]->cd(2);


    dPhi_Sub_Lead_NoBkg_up[0][0][l]->Draw();
    dPhi_Sub_Lead_NoBkg_up[0][0][l]->GetYaxis()->SetLabelSize(0.);
    dPhi_Sub_Lead_NoBkg_up[0][0][l]->GetYaxis()->SetNdivisions(408);
    dPhi_Sub_Lead_NoBkg_down[0][0][l]->Draw("same");

    sub_lead_dPhi_syst[2][0][0][l]->Draw("same e2");
    sub_lead_dPhi_tot[0][0][0][l]->Draw("same");
    label_per->Draw();
    l_phi->Draw();

    c_nobkgsub[l]->cd(3);
   
    dPhi_Sub_Lead_NoBkg_up[0][2][l]->Draw();
    dPhi_Sub_Lead_NoBkg_up[0][2][l]->GetYaxis()->SetLabelSize(0.);
    dPhi_Sub_Lead_NoBkg_up[0][2][l]->GetYaxis()->SetNdivisions(408);
    dPhi_Sub_Lead_NoBkg_down[0][2][l]->Draw("same");

    sub_lead_dPhi_syst[2][0][2][l]->Draw("same e2");
    sub_lead_dPhi_tot[0][0][2][l]->Draw("same");
    label_cent->Draw();

    orientation->Draw();
 
    l_phi->Draw();


    c_nobkgsub[l]->cd(5);
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->Draw();
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetYaxis()->SetLabelSize(0.0); 
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetYaxis()->CenterTitle();
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetYaxis()->SetTitle("#DeltaP_{PbPb}- #DeltaP_{pp}  (GeV/c)");
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetYaxis()->SetTitleSize(0.0);
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetYaxis()->SetTitleOffset(0.5);

    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetXaxis()->SetLabelSize(0.08); 
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetXaxis()->CenterTitle();
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetXaxis()->SetTitle("#Delta#phi");
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetXaxis()->SetTitleSize(0.09);
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetXaxis()->SetTitleOffset(0.7);
    dPhi_Sub_Lead_NoBkg_up[6][0][l]->GetXaxis()->CenterTitle();
    dPhi_Sub_Lead_NoBkg_down[6][0][l]->Draw("same");

    sub_lead_dPhi_syst[8][0][0][l]->Draw("same e2");
    sub_lead_dPhi_tot[6][0][0][l]->Draw("same");

    label = new TLatex(0.05,0.9,"PbPb - pp");
    label->SetTextSize(0.08);
    label->SetLineColor(kWhite);
    label->SetNDC();
    label->Draw();

    l_phi->Draw();
    cover_x_r->Draw();
    cover_x_l->Draw();

    c_nobkgsub[l]->cd(6);

    dPhi_Sub_Lead_NoBkg_up[6][2][l]->Draw();
    dPhi_Sub_Lead_NoBkg_up[6][2][l]->GetYaxis()->SetLabelSize(0.);
    dPhi_Sub_Lead_NoBkg_up[6][2][l]->GetXaxis()->SetLabelSize(0.08); 
    dPhi_Sub_Lead_NoBkg_up[6][2][l]->GetXaxis()->CenterTitle();
    dPhi_Sub_Lead_NoBkg_up[6][2][l]->GetXaxis()->SetTitle("#Delta#phi");
    dPhi_Sub_Lead_NoBkg_up[6][2][l]->GetXaxis()->SetTitleSize(0.09);
    dPhi_Sub_Lead_NoBkg_up[6][2][l]->GetXaxis()->SetTitleOffset(0.7);
    dPhi_Sub_Lead_NoBkg_up[6][2][l]->GetXaxis()->CenterTitle();
    dPhi_Sub_Lead_NoBkg_down[6][2][l]->Draw("same");

    sub_lead_dPhi_syst[8][0][2][l]->Draw("same e2");
    sub_lead_dPhi_tot[6][0][2][l]->Draw("same");

    l_phi->Draw();
    
    cover_x_r->Draw();
    cover_x_l->Draw();
  
    c_nobkgsub[l]->cd(4);
   dummy_axis_double_diff->Draw();
   
 
   if(l<2){
     c_jet[l]->Print("Missing_pT_JetRelated_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png");
     c_jet[l]->Print("Missing_pT_JetRelated_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf");

     c_longrange[l]->Print("Missing_pT_LongRange_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png");
     c_longrange[l]->Print("Missing_pT_LongRange_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf");

     c_nobkgsub[l]->Print("Missing_pT_NoBkgSub_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png");
     c_nobkgsub[l]->Print("Missing_pT_NoBkgSub_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf");
 
   }else{
     c_jet[l]->Print("Missing_pT_JetRelated_AjInclusive.png");
     c_jet[l]->Print("Missing_pT_JetRelated_AjInclusive.pdf");

     c_longrange[l]->Print("Missing_pT_LongRange_AjInclusive.png");
     c_longrange[l]->Print("Missing_pT_LongRange_AjInclusive.pdf");

     c_nobkgsub[l]->Print("Missing_pT_NoBkgSub_AjInclusive.png");
     c_nobkgsub[l]->Print("Missing_pT_NoBkgSub_AjInclusive.pdf");


   }
  

   f_out->cd();

 
   background_syst_rebin[0][4][2][l]->Write();
   background_syst_rebin[1][4][3][l]->Write();
   background_syst_rebin[7][4][2][l]->Write();
   background_syst_rebin[0][4][0][l]->Write();
   background_syst_rebin[7][4][0][l]->Write();
   signal_dPhi_syst[4][0][0][l]->Write();
   signal_dPhi_syst[2][0][0][l]->Write();
   signal_dPhi_syst[2][0][2][l]->Write();
   signal_dPhi_syst[4][0][2][l]->Write();
   signal_dPhi_syst[3][0][3][l]->Write();
   signal_dPhi_syst[5][0][3][l]->Write();
   sub_lead_dPhi_syst[2][0][0][l]->Write();
   sub_lead_dPhi_syst[2][0][2][l]->Write();
   sub_lead_dPhi_syst[3][0][3][l]->Write();
   //-----------------------------
   //  Make Integrals
   //-----------------------------

   Integral[1][3][l] = new TH1D(Form("Integral_SideBand_pp_Aj%d",l),"",6,pt_bin_bounds);
   Integral[3][3][l] = new TH1D(Form("Integral_SubLeading_pp_Aj%d",l),"",6,pt_bin_bounds);
   Integral[5][3][l] = new TH1D(Form("Integral_Leading_pp_Aj%d",l),"",6,pt_bin_bounds);



   Integral_noerr_up[1][3][l] = new TH1D(Form("Integral_SideBand_pp_Aj%d_Noerr_Up",l),"",6,pt_bin_bounds);
   Integral_noerr_up[3][3][l] = new TH1D(Form("Integral_SubLeading_pp_Aj%d_Noerr_Up",l),"",6,pt_bin_bounds);
   Integral_noerr_up[5][3][l] = new TH1D(Form("Integral_Leading_pp_Aj%d_Noerr_Up",l),"",6,pt_bin_bounds);


   Integral_noerr_down[1][3][l] = new TH1D(Form("Integral_SideBand_pp_Aj%d_Noerr_Down",l),"",6,pt_bin_bounds);
   Integral_noerr_down[3][3][l] = new TH1D(Form("Integral_SubLeading_pp_Aj%d_Noerr_Down",l),"",6,pt_bin_bounds);
   Integral_noerr_down[5][3][l] = new TH1D(Form("Integral_Leading_pp_Aj%d_Noerr_Down",l),"",6,pt_bin_bounds);



   Integral_syst[1][3][l] = (TH1D*) Integral[1][3][l]->Clone(Form("Integral_Syst_SideBand_pp_Aj%d",l));
   Integral_syst[3][3][l] = (TH1D*) Integral[1][3][l]->Clone(Form("Integral_Syst_SubLeading_pp_Aj%d",l));
   Integral_syst[5][3][l] = (TH1D*) Integral[1][3][l]->Clone(Form("Integral_Syst_Leading_pp_Aj%d",l));
 
  



   for(int i = 0; i<6; i++){

     integral = background_diff_rebin[1][i][0][l]->IntegralAndError(1,background_diff_rebin[1][i][0][l]->GetNbinsX(),int_err,"width");
     Integral[1][3][l]->SetBinContent(i+1,integral);
     Integral[1][3][l]->SetBinError(i+1,int_err);
   
     if(integral>0){
       Integral_noerr_up[1][3][l]->SetBinContent(i+1,integral);
       Integral_noerr_down[1][3][l]->SetBinContent(i+1,0.);
     }else{
       Integral_noerr_down[1][3][l]->SetBinContent(i+1,integral);
       Integral_noerr_up[1][3][l]->SetBinContent(i+1,0.);
     }

   

     Integral_syst[1][3][l]->SetBinContent(i+1,integral);
     Integral_syst[1][3][l]->SetBinError(i+1,TMath::Sqrt(integral*integral*(0.05*0.05+0.04*0.04+0.03*0.03+0.02*0.02)+0.06*0.06+0.22*0.22));
    

     integral = signal_dPhi_rebin[3][i][0][l]->IntegralAndError(1,signal_dPhi_rebin[3][i][0][l]->GetNbinsX(),int_err,"width");
     Integral[3][3][l]->SetBinContent(i+1,integral);
     Integral[3][3][l]->SetBinError(i+1,int_err);
 
     if(integral>0){
       Integral_noerr_up[3][3][l]->SetBinContent(i+1,integral);
       Integral_noerr_down[3][3][l]->SetBinContent(i+1,0.);
     }else{
       Integral_noerr_down[3][3][l]->SetBinContent(i+1,integral);
       Integral_noerr_up[3][3][l]->SetBinContent(i+1,0.);
     }


     Integral_syst[3][3][l]->SetBinContent(i+1,integral);
     Integral_syst[3][3][l]->SetBinError(i+1,TMath::Sqrt(integral*integral*(0.05*0.05+0.04*0.04+0.03*0.03+0.02*0.02)+.06*.06+0.22*0.22));

     integral = signal_dPhi_rebin[5][i][0][l]->IntegralAndError(1,signal_dPhi_rebin[5][i][0][l]->GetNbinsX(),int_err,"width");
     Integral[5][3][l]->SetBinContent(i+1,integral);
     Integral[5][3][l]->SetBinError(i+1,int_err);

     if(integral>0){
       Integral_noerr_up[5][3][l]->SetBinContent(i+1,integral);
       Integral_noerr_down[5][3][l]->SetBinContent(i+1,0.);
     }else{
       Integral_noerr_down[5][3][l]->SetBinContent(i+1,integral);
       Integral_noerr_up[5][3][l]->SetBinContent(i+1,0.);
     }


     Integral_syst[5][3][l]->SetBinContent(i+1,integral);
     Integral_syst[5][3][l]->SetBinError(i+1,TMath::Sqrt(integral*integral*(0.05*0.05+0.04*0.04+0.03*0.03+0.02*0.02)+0.06*0.06+0.22*0.22));


   }

   Integral_syst[1][3][l]->SetMarkerSize(0.);
   Integral_syst[1][3][l]->SetFillColor(kGray+2);
   Integral_syst[3][3][l]->SetMarkerSize(0.);
   Integral_syst[3][3][l]->SetFillColor(kGray+2);
   Integral_syst[5][3][l]->SetMarkerSize(0.);
   Integral_syst[5][3][l]->SetFillColor(kGray+2);
  



   //Integral[5][3][l]->Scale(-1.);	
   //Integral_syst[5][3][l]->Scale(-1.);

   for(int j = 0; j<4; j++){
 
     Integral[0][j][l] = new TH1D(Form("Integral_SideBand_PbPb_Cent%d_Aj%d",j,l),"",6,pt_bin_bounds);
     Integral[2][j][l] = new TH1D(Form("Integral_SubLeading_PbPb_Cent%d_Aj%d",j,l),"",6,pt_bin_bounds);
     Integral[4][j][l] = new TH1D(Form("Integral_Leading_PbPb_Cent%d_Aj%d",j,l),"",6,pt_bin_bounds);

     Integral_noerr_up[0][j][l] = new TH1D(Form("Integral_SideBand_PbPb_Cent%d_Aj%d_NoerrUp",j,l),"",6,pt_bin_bounds);
     Integral_noerr_up[2][j][l] = new TH1D(Form("Integral_SubLeading_PbPb_Cent%d_Aj%d_NoerrUp",j,l),"",6,pt_bin_bounds);
     Integral_noerr_up[4][j][l] = new TH1D(Form("Integral_Leading_PbPb_Cent%d_Aj%d_NoerrUp",j,l),"",6,pt_bin_bounds);


     Integral_noerr_down[0][j][l] = new TH1D(Form("Integral_SideBand_PbPb_Cent%d_Aj%d_NoerrDown",j,l),"",6,pt_bin_bounds);
     Integral_noerr_down[2][j][l] = new TH1D(Form("Integral_SubLeading_PbPb_Cent%d_Aj%d_NoerrDown",j,l),"",6,pt_bin_bounds);
     Integral_noerr_down[4][j][l] = new TH1D(Form("Integral_Leading_PbPb_Cent%d_Aj%d_NoerrDown",j,l),"",6,pt_bin_bounds);



     Integral_syst[0][j][l] = (TH1D*) Integral[0][j][l]->Clone(Form("Integral_Syst_SideBand_PbPb_Cent%d_Aj%d",j,l));
     Integral_syst[2][j][l] = (TH1D*) Integral[2][j][l]->Clone(Form("Integral_Syst_SubLeading_PbPb_Cent%d_Aj%d",j,l));
     Integral_syst[4][j][l] = (TH1D*) Integral[4][j][l]->Clone(Form("Integral_Syst_Leading_PbPb_Cent%d_Aj%d",j,l));
 
  
   
     for(int i = 0; i<6; i++){

       integral = background_diff_rebin[0][i][j][l]->IntegralAndError(1,background_diff_rebin[0][i][j][l]->GetNbinsX(),int_err,"width");
       Integral[0][j][l]->SetBinContent(i+1,integral);
       Integral[0][j][l]->SetBinError(i+1,int_err);

       if(integral>0){
	 Integral_noerr_up[0][j][l]->SetBinContent(i+1,integral);
	 Integral_noerr_down[0][j][l]->SetBinContent(i+1,0.);
       }else{
	 Integral_noerr_down[0][j][l]->SetBinContent(i+1,integral);
	 Integral_noerr_up[0][j][l]->SetBinContent(i+1,0.);
       }

       Integral_syst[0][j][l]->SetBinContent(i+1,integral);
       Integral_syst[0][j][l]->SetBinError(i+1,TMath::Sqrt(integral*integral*(0.05*0.05+0.04*0.04+0.03*0.03+0.02*0.02)+0.36*0.36+1.05*1.05));
   
       integral = signal_dPhi_rebin[2][i][j][l]->IntegralAndError(1,signal_dPhi_rebin[2][i][j][l]->GetNbinsX(),int_err,"width");
       Integral[2][j][l]->SetBinContent(i+1,integral);
       Integral[2][j][l]->SetBinError(i+1,int_err);


       if(integral>0){
	 Integral_noerr_up[2][j][l]->SetBinContent(i+1,integral);
	 Integral_noerr_down[2][j][l]->SetBinContent(i+1,0.);
       }else{
	 Integral_noerr_down[2][j][l]->SetBinContent(i+1,integral);
	 Integral_noerr_up[2][j][l]->SetBinContent(i+1,0.);
       }


       Integral_syst[2][j][l]->SetBinContent(i+1,integral);
       Integral_syst[2][j][l]->SetBinError(i+1,TMath::Sqrt(integral*integral*(0.05*0.05+0.04*0.04+0.03*0.03+0.02*0.02)+0.36*0.36+1.05*1.05));

     

       integral= signal_dPhi_rebin[4][i][j][l]->IntegralAndError(1,signal_dPhi_rebin[4][i][j][l]->GetNbinsX(),int_err,"width");
       Integral[4][j][l]->SetBinContent(i+1,integral);
       Integral[4][j][l]->SetBinError(i+1,int_err);


       if(integral>0){
	 Integral_noerr_up[4][j][l]->SetBinContent(i+1,integral);
	 Integral_noerr_down[4][j][l]->SetBinContent(i+1,0.);
       }else{
	 Integral_noerr_down[4][j][l]->SetBinContent(i+1,integral);
	 Integral_noerr_up[4][j][l]->SetBinContent(i+1,0.);
       }

       Integral_syst[4][j][l]->SetBinContent(i+1,integral);
       Integral_syst[4][j][l]->SetBinError(i+1,TMath::Sqrt(integral*integral*(0.05*0.05+0.04*0.04+0.03*0.03+0.02*0.02)+0.36*0.36+1.05*1.05));
    

     }
   
     
     Integral[6][j][l] = (TH1D*)Integral[0][j][l]->Clone(Form("Integral_SideBand_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral[6][j][l]->Add(Integral[1][3][l],-1.);
 
     Integral[8][j][l] = (TH1D*)Integral[2][j][l]->Clone(Form("Integral_SubLeading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral[8][j][l]->Add(Integral[3][3][l],-1.);
 
     Integral[10][j][l] = (TH1D*)Integral[4][j][l]->Clone(Form("Integral_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral[10][j][l]->Add(Integral[5][3][l],-1.);

     Integral[9][j][l] = (TH1D*)Integral[8][j][l]->Clone(Form("Integral_Sub_minus_leading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral[9][j][l]->Add(Integral[10][j][l]);
 


     Integral_syst[6][j][l] = (TH1D*)Integral_syst[0][j][l]->Clone(Form("Integral_syst_SideBand_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_syst[6][j][l]->Add(Integral_syst[1][3][l],-1.);
 
     Integral_syst[8][j][l] = (TH1D*)Integral_syst[2][j][l]->Clone(Form("Integral_syst_SubLeading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_syst[8][j][l]->Add(Integral_syst[3][3][l],-1.);
 
     Integral_syst[10][j][l] = (TH1D*)Integral_syst[4][j][l]->Clone(Form("Integral_syst_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_syst[10][j][l]->Add(Integral_syst[5][3][l],-1.);


     Integral_syst[9][j][l] = (TH1D*)Integral_syst[8][j][l]->Clone(Form("Integral_syst_SubLeading_minus_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_syst[9][j][l]->Add(Integral_syst[10][j][l]);
 



     Integral_syst[0][j][l]->SetMarkerSize(1.5);
     Integral_syst[0][j][l]->SetMarkerStyle(20);
     Integral_syst[0][j][l]->SetMarkerColor(kWhite);
     Integral_syst[0][j][l]->SetFillColor(kBlack);
     Integral_syst[0][j][l]->SetFillStyle(3004);
     Integral_syst[2][j][l]->SetMarkerSize(1.5);
     Integral_syst[2][j][l]->SetMarkerStyle(20);
     Integral_syst[2][j][l]->SetMarkerColor(kWhite);
     Integral_syst[2][j][l]->SetFillColor(kBlack);
     Integral_syst[2][j][l]->SetFillStyle(3004);
     Integral_syst[4][j][l]->SetMarkerSize(1.5);
     Integral_syst[4][j][l]->SetMarkerStyle(20);
     Integral_syst[4][j][l]->SetMarkerColor(kWhite);
     Integral_syst[4][j][l]->SetFillColor(kBlack);
     Integral_syst[4][j][l]->SetFillStyle(3004);
     Integral_syst[6][j][l]->SetFillColor(kBlack);
     Integral_syst[6][j][l]->SetFillStyle(3004);
     Integral_syst[6][j][l]->SetMarkerSize(1.5);
     Integral_syst[6][j][l]->SetMarkerStyle(20);
     Integral_syst[6][j][l]->SetMarkerColor(kWhite);
     Integral_syst[8][j][l]->SetMarkerSize(1.5);
     Integral_syst[8][j][l]->SetMarkerStyle(20);
     Integral_syst[8][j][l]->SetMarkerColor(kBlue-9);
     Integral_syst[8][j][l]->SetFillColor(kBlack);
     Integral_syst[8][j][l]->SetFillStyle(3004);
     Integral_syst[9][j][l]->SetFillColor(kBlack);
     Integral_syst[9][j][l]->SetFillStyle(3004);
     Integral_syst[9][j][l]->SetMarkerSize(1.5);
     Integral_syst[9][j][l]->SetMarkerStyle(20);
     Integral_syst[9][j][l]->SetMarkerColor(kWhite);
     Integral_syst[10][j][l]->SetMarkerSize(1.5);
     Integral_syst[10][j][l]->SetMarkerStyle(20);
     Integral_syst[10][j][l]->SetMarkerColor(kWhite);
     Integral_syst[10][j][l]->SetFillColor(kBlack);
     Integral_syst[10][j][l]->SetFillStyle(3004);
 


     Integral[0][j][l]->SetMarkerSize(1.5);
     Integral[0][j][l]->SetMarkerColor(kBlack);
     Integral[0][j][l]->SetLineColor(kBlack);
     Integral[0][j][l]->SetMarkerStyle(4);
     Integral[2][j][l]->SetMarkerSize(1.5);
     Integral[2][j][l]->SetMarkerColor(kBlack);
     Integral[2][j][l]->SetLineColor(kBlack);
     Integral[2][j][l]->SetMarkerStyle(4);
     Integral[4][j][l]->SetMarkerSize(1.5);
     Integral[4][j][l]->SetMarkerColor(kBlack);
     Integral[4][j][l]->SetLineColor(kBlack);
     Integral[4][j][l]->SetMarkerStyle(4);
     Integral[6][j][l]->SetMarkerSize(1.5);
     Integral[6][j][l]->SetMarkerColor(kBlack);
     Integral[6][j][l]->SetLineColor(kBlack);
     Integral[6][j][l]->SetMarkerStyle(4);
     Integral[8][j][l]->SetMarkerSize(1.5);
     Integral[8][j][l]->SetMarkerColor(kBlack);
     Integral[8][j][l]->SetLineColor(kBlack);
     Integral[8][j][l]->SetMarkerStyle(4);
     Integral[9][j][l]->SetMarkerSize(1.5);
     Integral[9][j][l]->SetMarkerColor(kBlack);
     Integral[9][j][l]->SetLineColor(kBlack);
     Integral[9][j][l]->SetMarkerStyle(4);
     Integral[10][j][l]->SetMarkerSize(1.5);
     Integral[10][j][l]->SetMarkerColor(kBlack); 
     Integral[10][j][l]->SetLineColor(kBlack);
     Integral[10][j][l]->SetMarkerStyle(4);
 



 
     Integral_noerr_up[6][j][l] = (TH1D*)Integral[6][j][l]->Clone(Form("Integral_noerr_up_SideBand_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_noerr_up[8][j][l] = (TH1D*)Integral[8][j][l]->Clone(Form("Integral_noerr_up_SubLeading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_noerr_up[9][j][l] = (TH1D*)Integral[9][j][l]->Clone(Form("Integral_noerr_up_SubLeading_minus_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_noerr_up[10][j][l] = (TH1D*)Integral[10][j][l]->Clone(Form("Integral_noerr_up_Leading_minus_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));

     Integral_noerr_down[6][j][l] = (TH1D*)Integral[6][j][l]->Clone(Form("Integral_noerr_down_SideBand_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_noerr_down[8][j][l] = (TH1D*)Integral[8][j][l]->Clone(Form("Integral_noerr_down_SubLeading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_noerr_down[9][j][l] = (TH1D*)Integral[9][j][l]->Clone(Form("Integral_noerr_down_SubLeading_minus_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));
     Integral_noerr_down[10][j][l] = (TH1D*)Integral[10][j][l]->Clone(Form("Integral_noerr_down_Leading_PbPb_minus_ppCent%d_Aj%d",j,l));




     for(int k = 1; k< Integral_noerr_up[6][j][l]->GetNbinsX()+1; k++){
       Integral_noerr_up[6][j][l]->SetBinError(k,0.);
       Integral_noerr_up[8][j][l]->SetBinError(k,0.);
       Integral_noerr_up[9][j][l]->SetBinError(k,0.);
       Integral_noerr_up[10][j][l]->SetBinError(k,0.);
   
       Integral_noerr_down[6][j][l]->SetBinError(k,0.);
       Integral_noerr_down[8][j][l]->SetBinError(k,0.);
       Integral_noerr_down[9][j][l]->SetBinError(k,0.);
       Integral_noerr_down[10][j][l]->SetBinError(k,0.);
     

       if(  Integral[6][j][l]->GetBinContent(k)>0){
	 Integral_noerr_down[6][j][l]->SetBinContent(k,0.);
       }else{
	 Integral_noerr_up[6][j][l]->SetBinContent(k,0.);
       }


        if(  Integral[8][j][l]->GetBinContent(k)>0){
	 Integral_noerr_down[8][j][l]->SetBinContent(k,0.);
       }else{
	 Integral_noerr_up[8][j][l]->SetBinContent(k,0.);
       }
     

	if(  Integral[9][j][l]->GetBinContent(k)>0){
	 Integral_noerr_down[9][j][l]->SetBinContent(k,0.);
       }else{
	 Integral_noerr_up[9][j][l]->SetBinContent(k,0.);
       }


	if(  Integral[10][j][l]->GetBinContent(k)>0){
	  Integral_noerr_down[10][j][l]->SetBinContent(k,0.);
	}else{
	  Integral_noerr_up[10][j][l]->SetBinContent(k,0.);
	}


     }



  }

  


  }//l

  gStyle->SetPadBottomMargin(0.2);

   c_longrange_int = new TCanvas("LongrangeIntegrals","",10,10,800,400);
   c_longrange_int->Divide(2,1,0,0);
   


   float int_max = 14.6;
   float int_min = -2.;
  
   c_longrange_int->cd(1);
 
   Integral[1][3][0]->SetMarkerStyle(20); 
   Integral[1][3][0]->SetMarkerSize(1);
   Integral[1][3][0]->SetMarkerColor(kBlack);
   Integral[1][3][0]->SetLineColor(kBlack);
  
 
   Integral[1][3][1]->SetMarkerStyle(20); 
   Integral[1][3][1]->SetMarkerSize(1);
   Integral[1][3][1]->SetMarkerColor(kBlack);
   Integral[1][3][1]->SetLineColor(kBlack);

   Integral[0][0][0]->SetMarkerStyle(24); 
   Integral[0][0][0]->SetMarkerSize(1.5);
   Integral[0][0][0]->SetMarkerColor(kRed+2);
   Integral[0][0][0]->SetLineColor(kRed+2);



   Integral[0][0][1]->SetMarkerStyle(24); 
   Integral[0][0][1]->SetMarkerSize(1.5);
   Integral[0][0][1]->SetMarkerColor(kRed+2);
   Integral[0][0][1]->SetLineColor(kRed+2);



   Integral[0][2][0]->SetMarkerStyle(21); 
   Integral[0][2][0]->SetMarkerSize(1);
   Integral[0][2][0]->SetMarkerColor(kRed+2);
   Integral[0][2][0]->SetLineColor(kRed+2);
  
 
   Integral[0][2][1]->SetMarkerStyle(21); 
   Integral[0][2][1]->SetMarkerSize(1);
   Integral[0][2][1]->SetMarkerColor(kRed+2);
   Integral[0][2][1]->SetLineColor(kRed+2);
    
 
   Integral[0][0][0]->SetMaximum(int_max);
   Integral[0][0][0]->SetMinimum(int_min);
  
 
   Integral[0][0][0]->Draw();
   Integral[0][0][0]->GetXaxis()->SetTitle("p_{T}^{assoc}   (GeV/c)");
   Integral[0][0][0]->GetXaxis()->SetTitleSize(0.08);
   Integral[0][0][0]->GetXaxis()->SetTitleOffset(1.);
   Integral[0][0][0]->GetXaxis()->CenterTitle();
   Integral[0][0][0]->GetXaxis()->SetLabelSize(0.08); 
   Integral[0][0][0]->GetXaxis()->SetRangeUser(0.501,7.6);
 

   Integral[0][0][0]->GetYaxis()->SetTitle("#Sigma #DeltaP   (GeV/c)");
   Integral[0][0][0]->GetYaxis()->SetTitleSize(0.08);
   Integral[0][0][0]->GetYaxis()->CenterTitle();
   Integral[0][0][0]->GetYaxis()->SetTitleOffset(.7);
   Integral[0][0][0]->GetYaxis()->SetLabelSize(0.08); 
 


   Integral[0][0][0]->Draw();
 
   //   Integral_syst[0][0][0]->Draw("same e2");   
   //  Integral_syst[0][3][0]->SetFillStyle(0);
   Integral_syst[0][2][0]->Draw("same e2");   
  
   Integral_syst[1][3][0]->Draw("same e2");
   Integral[0][0][0]->Draw("same");
   Integral[0][2][0]->Draw("same");


   TLatex *aj_tex = new TLatex(0.2,0.9,"A_{J} < 0.22");
   aj_tex->SetTextSize(0.08);
   aj_tex->SetLineColor(kWhite);
   aj_tex->SetNDC();
   aj_tex->Draw(); 
 


  
   TLine *l_int = new TLine(0.5,0.,8.0,0.);
   l_int->SetLineStyle(2);
   l_int->Draw();
  

   TLegend *legend_int = new TLegend(0.17,0.6,0.97,0.85);

   Integral[0][2][0]->SetFillColor(kBlack);
   Integral[0][2][0]->SetFillStyle(3004);
 
   legend_int->AddEntry(Integral[0][0][0],"PbPb 50-100% ");
   legend_int->AddEntry(Integral[0][2][0],"PbPb 0-30%","lpfe");
   legend_int->AddEntry(Integral_syst[1][3][0],"pp Reference","f");
 
   legend_int->SetLineColor(kWhite);
   legend_int->SetTextSize(0.08);
   legend_int->Draw();

 
   TPave *cover_x_int = new TPave(0.95,-0.5,1.05,0.18);
   cover_x_int->SetFillColor(kWhite);
   cover_x_int->SetLineColor(kWhite);
   cover_x_int->SetOption("NDC NB");
   cover_x_int->Draw();

 

   c_longrange_int->cd(2);
 
   Integral[0][0][1]->SetMaximum(int_max);
   Integral[0][0][1]->SetMinimum(int_min);
 


   Integral[0][0][1]->Draw();
   Integral[0][0][1]->GetXaxis()->SetTitle("p_{T}^{assoc}   (GeV/c)");
   Integral[0][0][1]->GetXaxis()->SetTitleSize(0.08);
   Integral[0][0][1]->GetXaxis()->SetTitleOffset(1.);
   Integral[0][0][1]->GetXaxis()->CenterTitle();
   Integral[0][0][1]->GetXaxis()->SetLabelSize(0.08); 
   Integral[0][0][1]->GetXaxis()->SetRangeUser(0.501,7.99);
   Integral[0][0][1]->GetYaxis()->SetTitleSize(0.0);
   Integral[0][0][1]->GetYaxis()->SetTitleOffset(.8);
   Integral[0][0][1]->GetYaxis()->SetLabelSize(0.0);
  
   Integral[0][0][1]->Draw();
 
   //Integral_syst[0][0][1]->Draw("same e2");   
   Integral_syst[0][2][1]->Draw("same e2");   

   
   Integral_syst[1][3][1]->Draw("same e2");
   Integral[0][0][1]->Draw("same");
   Integral[0][2][1]->Draw("same");
 

   aj_tex = new TLatex(0.05,0.9,"A_{J} > 0.22");
   aj_tex->SetTextSize(0.08);
   aj_tex->SetLineColor(kWhite);
   aj_tex->SetNDC();
   aj_tex->Draw(); 
 

   TLatex *cms_tex = new TLatex(0.4,0.9,"CMS");
   cms_tex->SetTextSize(0.08);
   cms_tex->SetLineColor(kWhite);
   cms_tex->SetNDC();
   cms_tex->Draw(); 


    TLatex *prelim_tex = new TLatex(0.55,0.9,"Preliminary");
    prelim_tex->SetTextFont(53);
    prelim_tex->SetTextSizePixels(25);
    prelim_tex->SetLineColor(kWhite);
    prelim_tex->SetNDC();
    prelim_tex->Draw(); 
  

   l_int->Draw();

   cover_x_int->Draw();
   c_longrange_int->cd(0);

    
    TLatex *type_tex_int = new TLatex(0.1,0.93,"Long range asymmetry |#Delta#eta|<2.5");
    type_tex_int->SetTextSize(0.06);
    type_tex_int->SetLineColor(kWhite);
    type_tex_int->SetNDC();
    type_tex_int->Draw();
   
    TLatex   *luminosity_tex_int_pp = new TLatex(0.1,0.87,"pp 5.3 pb^{-1} (2.76 TeV)");
    luminosity_tex_int_pp->SetTextFont(43);
    luminosity_tex_int_pp->SetTextSizePixels(18);
    luminosity_tex_int_pp->SetLineColor(kWhite);
    luminosity_tex_int_pp->SetNDC();
    luminosity_tex_int_pp->Draw();
   
    TLatex   *luminosity_tex_int_PbPb = new TLatex(0.1,0.81,"PbPb 166 #mub^{-1} (2.76 TeV)");
    luminosity_tex_int_PbPb->SetTextFont(43);
    luminosity_tex_int_PbPb->SetTextSizePixels(18);
    luminosity_tex_int_PbPb->SetLineColor(kWhite);
    luminosity_tex_int_PbPb->SetNDC();
    luminosity_tex_int_PbPb->Draw();
    

    TLatex   *jet_reco_tex_int = new TLatex(0.55,0.87,"anti-k_{T} R = 0.3, |#eta_{jet}| < 1.6");
    jet_reco_tex_int->SetTextFont(43);
    jet_reco_tex_int->SetTextSizePixels(18);
    jet_reco_tex_int->SetLineColor(kWhite);
    jet_reco_tex_int->SetNDC();
    jet_reco_tex_int->Draw();

    TLatex   *jet_cut_tex_int = new TLatex(0.55,0.81,"120 < p_{T,1}< 300, p_{T,2}> 50 GeV/c, #Delta#phi_{1,2}> 5#pi/6");
    jet_cut_tex_int->SetTextFont(43);
    jet_cut_tex_int->SetTextSizePixels(18);
    jet_cut_tex_int->SetLineColor(kWhite);
    jet_cut_tex_int->SetNDC();
    jet_cut_tex_int->Draw();
    
  
   c_longrange_int->Print("Integral_Longrange.png");
   c_longrange_int->Print("Integral_Longrange.pdf");


   //---------------------
   //Integral summary
   //---------------------


 


  c_all_int = new TCanvas("IntegralSummary","",10,10,1000,800);
  c_all_int->Divide(2,2,0,0);
  
  c_all_int->cd(1);

 
    
  Integral_noerr_up[10][2][0]->SetMarkerSize(0);
  Integral_noerr_up[10][2][0]->SetFillColor(kYellow-9);
  Integral_noerr_up[10][2][1]->SetMarkerSize(0);
  Integral_noerr_up[10][2][1]->SetFillColor(kYellow-9);
  
  Integral_noerr_down[10][2][0]->SetMarkerSize(0);
  Integral_noerr_down[10][2][0]->SetFillColor(kYellow-9);
  Integral_noerr_down[10][2][1]->SetMarkerSize(0);
  Integral_noerr_down[10][2][1]->SetFillColor(kYellow-9);
 
  Integral_noerr_up[8][2][0]->SetMarkerSize(0);
  Integral_noerr_up[8][2][1]->SetMarkerSize(0);
  Integral_noerr_up[8][2][0]->SetFillColor(kGreen+3);
  Integral_noerr_up[8][2][1]->SetFillColor(kGreen+3);

  Integral_noerr_down[8][2][0]->SetMarkerSize(0);
  Integral_noerr_down[8][2][1]->SetMarkerSize(0);
  Integral_noerr_down[8][2][0]->SetFillColor(kGreen+3);
  Integral_noerr_down[8][2][1]->SetFillColor(kGreen+3);



 
  Integral_noerr_up[9][2][0]->SetMarkerSize(0);
  Integral_noerr_up[9][2][0]->SetFillColor(kBlue-9);
  Integral_noerr_up[9][2][1]->SetMarkerSize(0);
  Integral_noerr_up[9][2][1]->SetFillColor(kBlue-9);
 
  Integral_noerr_up[6][2][0]->SetMarkerSize(0);
  Integral_noerr_up[6][2][1]->SetMarkerSize(0);
  Integral_noerr_up[6][2][0]->SetFillColor(kRed+2);
  Integral_noerr_up[6][2][1]->SetFillColor(kRed+2);




 
  Integral_noerr_down[9][2][0]->SetMarkerSize(0);
  Integral_noerr_down[9][2][0]->SetFillColor(kBlue-9);
  Integral_noerr_down[9][2][1]->SetMarkerSize(0);
  Integral_noerr_down[9][2][1]->SetFillColor(kBlue-9);
 
  Integral_noerr_down[6][2][0]->SetMarkerSize(0);
  Integral_noerr_down[6][2][1]->SetMarkerSize(0);
  Integral_noerr_down[6][2][0]->SetFillColor(kRed+2);
  Integral_noerr_down[6][2][1]->SetFillColor(kRed+2);

 
 
   TString int_name = "int_jet_up_Aj0_Aj22";
 
   Integral_up[8][2][0] = new THStack(int_name, "");
   Integral_up[8][2][0]->Add(Integral_noerr_up[8][2][0]);
   Integral_up[8][2][0]->Add(Integral_noerr_up[10][2][0]);
 
  
   int_name.ReplaceAll("up","down");
   Integral_down[8][2][0] = new THStack(int_name, "");
   Integral_down[8][2][0]->Add(Integral_noerr_down[8][2][0]);
   Integral_down[8][2][0]->Add(Integral_noerr_down[10][2][0]);
   


   int_name = "int_jet_up_Aj22_Aj75";

   Integral_up[8][2][1] = new THStack(int_name, "");
   Integral_up[8][2][1]->Add(Integral_noerr_up[8][2][1]);
   Integral_up[8][2][1]->Add(Integral_noerr_up[10][2][1]);
  
  
   int_name.ReplaceAll("up","down");
   Integral_down[8][2][1] = new THStack(int_name, "");
   Integral_down[8][2][1]->Add(Integral_noerr_down[8][2][1]);
   Integral_down[8][2][1]->Add(Integral_noerr_down[10][2][1]);
  

  
   int_name = "int_diff_up_Aj0_Aj22";
 
   Integral_up[9][2][0] = new THStack(int_name, "");
   Integral_up[9][2][0]->Add(Integral_noerr_up[9][2][0]);
   Integral_up[9][2][0]->Add(Integral_noerr_up[6][2][0]);
 
  
   int_name.ReplaceAll("up","down");
   Integral_down[9][2][0] = new THStack(int_name, "");
   Integral_down[9][2][0]->Add(Integral_noerr_down[9][2][0]);
   Integral_down[9][2][0]->Add(Integral_noerr_down[6][2][0]);
   


   int_name = "int_diff_up_Aj22_Aj75";

   Integral_up[9][2][1] = new THStack(int_name, "");
   Integral_up[9][2][1]->Add(Integral_noerr_up[9][2][1]);
   Integral_up[9][2][1]->Add(Integral_noerr_up[6][2][1]);
  
  
   int_name.ReplaceAll("up","down");
   Integral_down[9][2][1] = new THStack(int_name, "");
   Integral_down[9][2][1]->Add(Integral_noerr_down[9][2][1]);
   Integral_down[9][2][1]->Add(Integral_noerr_down[6][2][1]);
  


   //peripheral


 
    
  Integral_noerr_up[10][0][0]->SetMarkerSize(0);
  Integral_noerr_up[10][0][0]->SetFillColor(kYellow-9);
  Integral_noerr_up[10][0][1]->SetMarkerSize(0);
  Integral_noerr_up[10][0][1]->SetFillColor(kYellow-9);
  
  Integral_noerr_down[10][0][0]->SetMarkerSize(0);
  Integral_noerr_down[10][0][0]->SetFillColor(kYellow-9);
  Integral_noerr_down[10][0][1]->SetMarkerSize(0);
  Integral_noerr_down[10][0][1]->SetFillColor(kYellow-9);
 
  Integral_noerr_up[8][0][0]->SetMarkerSize(0);
  Integral_noerr_up[8][0][1]->SetMarkerSize(0);
  Integral_noerr_up[8][0][0]->SetFillColor(kGreen+3);
  Integral_noerr_up[8][0][1]->SetFillColor(kGreen+3);

  Integral_noerr_down[8][0][0]->SetMarkerSize(0);
  Integral_noerr_down[8][0][1]->SetMarkerSize(0);
  Integral_noerr_down[8][0][0]->SetFillColor(kGreen+3);
  Integral_noerr_down[8][0][1]->SetFillColor(kGreen+3);



 
  Integral_noerr_up[9][0][0]->SetMarkerSize(0);
  Integral_noerr_up[9][0][0]->SetFillColor(kBlue-9);
  Integral_noerr_up[9][0][1]->SetMarkerSize(0);
  Integral_noerr_up[9][0][1]->SetFillColor(kBlue-9);
 
  Integral_noerr_up[6][0][0]->SetMarkerSize(0);
  Integral_noerr_up[6][0][1]->SetMarkerSize(0);
  Integral_noerr_up[6][0][0]->SetFillColor(kRed+2);
  Integral_noerr_up[6][0][1]->SetFillColor(kRed+2);




 
  Integral_noerr_down[9][0][0]->SetMarkerSize(0);
  Integral_noerr_down[9][0][0]->SetFillColor(kBlue-9);
  Integral_noerr_down[9][0][1]->SetMarkerSize(0);
  Integral_noerr_down[9][0][1]->SetFillColor(kBlue-9);
 
  Integral_noerr_down[6][0][0]->SetMarkerSize(0);
  Integral_noerr_down[6][0][1]->SetMarkerSize(0);
  Integral_noerr_down[6][0][0]->SetFillColor(kRed+2);
  Integral_noerr_down[6][0][1]->SetFillColor(kRed+2);

 
 
  int_name = "int_jet_up_Aj0_Aj22";
 
   Integral_up[8][0][0] = new THStack(int_name, "");
   Integral_up[8][0][0]->Add(Integral_noerr_up[8][0][0]);
   Integral_up[8][0][0]->Add(Integral_noerr_up[10][0][0]);
 
  
   int_name.ReplaceAll("up","down");
   Integral_down[8][0][0] = new THStack(int_name, "");
   Integral_down[8][0][0]->Add(Integral_noerr_down[8][0][0]);
   Integral_down[8][0][0]->Add(Integral_noerr_down[10][0][0]);
   


   int_name = "int_jet_up_Aj22_Aj75";

   Integral_up[8][0][1] = new THStack(int_name, "");
   Integral_up[8][0][1]->Add(Integral_noerr_up[8][0][1]);
   Integral_up[8][0][1]->Add(Integral_noerr_up[10][0][1]);
  
  
   int_name.ReplaceAll("up","down");
   Integral_down[8][0][1] = new THStack(int_name, "");
   Integral_down[8][0][1]->Add(Integral_noerr_down[8][0][1]);
   Integral_down[8][0][1]->Add(Integral_noerr_down[10][0][1]);
  

  
   int_name = "int_diff_up_Aj0_Aj22";
 
   Integral_up[9][0][0] = new THStack(int_name, "");
   Integral_up[9][0][0]->Add(Integral_noerr_up[9][0][0]);
   Integral_up[9][0][0]->Add(Integral_noerr_up[6][0][0]);
 
  
   int_name.ReplaceAll("up","down");
   Integral_down[9][0][0] = new THStack(int_name, "");
   Integral_down[9][0][0]->Add(Integral_noerr_down[9][0][0]);
   Integral_down[9][0][0]->Add(Integral_noerr_down[6][0][0]);
   


   int_name = "int_diff_up_Aj22_Aj75";

   Integral_up[9][0][1] = new THStack(int_name, "");
   Integral_up[9][0][1]->Add(Integral_noerr_up[9][0][1]);
   Integral_up[9][0][1]->Add(Integral_noerr_up[6][0][1]);
  
  
   int_name.ReplaceAll("up","down");
   Integral_down[9][0][1] = new THStack(int_name, "");
   Integral_down[9][0][1]->Add(Integral_noerr_down[9][0][1]);
   Integral_down[9][0][1]->Add(Integral_noerr_down[6][0][1]);
  





   ///DRAWING//
   c_all_int->cd(1);


   float summary_int_max = 16.;
   float summary_int_min = -13.;
 
   float summary_int_max2 = 16.;
   float summary_int_min2 = -13.;
 

   Integral_up[8][2][0]->Draw();
  

   Integral_up[8][2][0]->SetMaximum(summary_int_max);
   Integral_up[8][2][0]->SetMinimum(summary_int_min);
   Integral_up[8][2][0]->GetYaxis()->SetLabelSize(0.085);
   Integral_up[8][2][0]->GetYaxis()->SetTitleSize(0.085);
   Integral_up[8][2][0]->GetYaxis()->SetTitleOffset(0.7);
   Integral_up[8][2][0]->GetYaxis()->SetTitle("#Sigma P_{PbPb}- #Sigma P_{pp}   (GeV/c)");
    Integral_up[8][2][0]->GetYaxis()->CenterTitle();
   Integral_up[8][2][0]->GetYaxis()->SetTitleOffset(0.8);

   Integral_up[8][2][0]->GetXaxis()->SetLabelSize(0.0); 
   Integral_up[8][2][0]->GetXaxis()->SetRangeUser(0.501,7.99);
  


   Integral_down[8][2][0]->Draw("same ");


   Integral_syst[8][2][0]->Add( Integral_syst[10][2][0]);
   Integral_syst[8][2][0]->Draw("same e2");


   Integral[8][2][0]->Add(   Integral[10][2][0]);
   Integral[8][2][0]->Draw("same");
   

   l_int->Draw();



   TLegend *legend_int_lead = new TLegend(0.17,0.65,0.95,0.85);
   legend_int_lead->AddEntry(Integral_noerr_up[8][2][0],"Subleading jet ","f");
   legend_int_lead->AddEntry(Integral_noerr_up[10][2][0],"Leading jet ","f");
   legend_int_lead->AddEntry(Integral_syst[8][2][0],"Subleading+leading jet ","lpfe");
   legend_int_lead->SetTextSize(0.07);
   legend_int_lead->SetLineColor(kWhite);
   legend_int_lead->Draw();
   

   
   aj_tex = new TLatex(0.2,0.9,"A_{J}<0.22");
   aj_tex->SetTextSize(0.08);
   aj_tex->SetLineColor(kWhite);
   aj_tex->SetNDC();
   aj_tex->Draw(); 
 

   TLine *l_int2 = new TLine(0.5,0.,8.0,0.);
   l_int2->SetLineStyle(2);
   l_int2->Draw();



   c_all_int->cd(2);

   Integral_up[8][2][1]->Draw();
  

   Integral_up[8][2][1]->SetMaximum(summary_int_max);
   Integral_up[8][2][1]->SetMinimum(summary_int_min);
   Integral_up[8][2][1]->GetYaxis()->SetLabelSize(0.0);
   Integral_up[8][2][1]->GetYaxis()->SetTitleSize(0.0);
 

   Integral_up[8][2][1]->GetXaxis()->SetLabelSize(0.0); 
   Integral_up[8][2][1]->GetXaxis()->SetRangeUser(0.501,7.99);
  
   Integral_down[8][2][1]->Draw("same ");


   Integral_syst[8][2][1]->Add(Integral_syst[10][2][1]);
   Integral_syst[8][2][1]->Draw("same e2");


   Integral[8][2][1]->Add( Integral[10][2][1]);
   Integral[8][2][1]->Draw("same");
  

   l_int->Draw();

   aj_tex = new TLatex(0.05,0.9,"A_{J}>0.22");
   aj_tex->SetTextSize(0.08);
   aj_tex->SetLineColor(kWhite);
   aj_tex->SetNDC();
   aj_tex->Draw(); 
  
  
   cms_tex->Draw(); 
   prelim_tex->Draw(); 




   TLatex *cent_tex = new TLatex(0.4,0.7,"0-30% Central");
   cent_tex->SetTextSize(0.08);
   cent_tex->SetLineColor(kWhite);
   cent_tex->SetNDC();
   cent_tex->Draw(); 
  

   TLatex *data_tex = new TLatex(0.4,0.8,"PbPb - pp");
   data_tex->SetTextSize(0.1);
   data_tex->SetLineColor(kWhite);
   data_tex->SetNDC();
   data_tex->Draw(); 


   l_int2->Draw();
    c_all_int->cd(3);

   cout<<"here"<<endl;
   Integral_up[9][2][0]->Draw();


   Integral_up[9][2][0]->SetMaximum(summary_int_max2);
   Integral_up[9][2][0]->SetMinimum(summary_int_min2);
   Integral_up[9][2][0]->GetYaxis()->SetLabelSize(0.07);
   Integral_up[9][2][0]->GetYaxis()->SetTitleSize(0.07);
   Integral_up[9][2][0]->GetYaxis()->SetTitleOffset(.9);
   Integral_up[9][2][0]->GetYaxis()->SetTitle("#Sigma #DeltaP_{PbPb}- #Sigma #DeltaP_{pp}   (GeV/c)");
   Integral_up[9][2][0]->GetYaxis()->CenterTitle();
   Integral_up[9][2][0]->GetXaxis()->SetTitle("p_{T}^{assoc}   (GeV/c)");
   Integral_up[9][2][0]->GetXaxis()->SetTitleSize(0.08);
   Integral_up[9][2][0]->GetXaxis()->CenterTitle();
   Integral_up[9][2][0]->GetXaxis()->SetLabelSize(0.08); 
   Integral_up[9][2][0]->GetXaxis()->SetTitleOffset(1.);
   Integral_up[9][2][0]->GetXaxis()->SetRangeUser(.501,7.99);
   
   Integral_down[9][2][0]->Draw("same");


   Integral_syst[9][2][0]->Add(Integral_syst[6][2][0]);
   Integral_syst[9][2][0]->Draw("same e2");

   Integral[9][2][0]->Add(Integral[6][2][0]);
   Integral[9][2][0]->Draw("same");
     


   l_int->Draw();

   TLegend *legend_int_diff = new TLegend(0.17,0.75,0.95,0.95);
   legend_int_diff->AddEntry(Integral_noerr_up[9][2][0],"Subleading+leading jet ","f");
   legend_int_diff->AddEntry(Integral_noerr_up[6][2][0],"Long range","f");
   legend_int_diff->AddEntry(Integral_syst[9][2][0],"Total","lpfe");
   legend_int_diff->SetTextSize(0.062);
   legend_int_diff->SetLineColor(kWhite);
   legend_int_diff->Draw();

   l_int2->Draw();
 cover_x_int->Draw();


 c_all_int->cd(4);
 
    Integral_up[9][2][1]->Draw();


   Integral_up[9][2][1]->SetMaximum(summary_int_max2);
   Integral_up[9][2][1]->SetMinimum(summary_int_min2);
   Integral_up[9][2][1]->GetYaxis()->SetLabelSize(0.0);
   Integral_up[9][2][1]->GetYaxis()->SetTitleSize(0.0);
    Integral_up[9][2][1]->GetXaxis()->SetTitle("p_{T}^{assoc}   (GeV/c)");
   Integral_up[9][2][1]->GetXaxis()->SetTitleSize(0.08);
   Integral_up[9][2][1]->GetXaxis()->CenterTitle();
   Integral_up[9][2][1]->GetXaxis()->SetLabelSize(0.08); 
   Integral_up[9][2][1]->GetYaxis()->SetTitleOffset(0.8);
   Integral_up[9][2][1]->GetXaxis()->SetRangeUser(.501,7.99);
   
   Integral_down[9][2][1]->Draw("same");


   Integral_syst[9][2][1]->Add(Integral_syst[6][2][1]);
   Integral_syst[9][2][1]->Draw("same e2");


   Integral[9][2][1]->Add( Integral[6][2][1]);
   Integral[9][2][1]->Draw("same");
     



   l_int2->Draw();
   cover_x_int->Draw();


   c_all_int->cd(0);

   
    TLatex   *luminosity_tex_all_pp = new TLatex(0.1,0.96,"pp 5.3 pb^{-1} (2.76 TeV)");
    luminosity_tex_all_pp->SetTextFont(43);
    luminosity_tex_all_pp->SetTextSizePixels(25);
    luminosity_tex_all_pp->SetLineColor(kWhite);
    luminosity_tex_all_pp->SetNDC();
    luminosity_tex_all_pp->Draw();
   
    TLatex   *luminosity_tex_all_PbPb = new TLatex(0.1,0.92,"PbPb 166 #mub^{-1} (2.76 TeV)");
    luminosity_tex_all_PbPb->SetTextFont(43);
    luminosity_tex_all_PbPb->SetTextSizePixels(25);
    luminosity_tex_all_PbPb->SetLineColor(kWhite);
    luminosity_tex_all_PbPb->SetNDC();
    luminosity_tex_all_PbPb->Draw();
    

    TLatex   *jet_reco_tex_all = new TLatex(0.53,0.96,"anti-k_{T} R = 0.3, |#eta_{jet}| < 1.6");
    jet_reco_tex_all->SetTextFont(43);
    jet_reco_tex_all->SetTextSizePixels(25);
    jet_reco_tex_all->SetLineColor(kWhite);
    jet_reco_tex_all->SetNDC();
    jet_reco_tex_all->Draw();

    TLatex   *jet_cut_tex_all = new TLatex(0.53,0.92,"120 < p_{T,1}< 300, p_{T,2}> 50 GeV/c, #Delta#phi_{1,2}> 5#pi/6");
    jet_cut_tex_all->SetTextFont(43);
    jet_cut_tex_all->SetTextSizePixels(25);
    jet_cut_tex_all->SetLineColor(kWhite);
    jet_cut_tex_all->SetNDC();
    jet_cut_tex_all->Draw();
    


   c_all_int->Print("Integral_PbPb_pp_Summary.png");
   c_all_int->Print("Integral_PbPb_pp_Summary.pdf");



   //peripheral

   c_all_int->cd(1);


   Integral_up[8][0][0]->Draw();
  

  
   Integral_up[8][0][0]->SetMaximum(summary_int_max);
   Integral_up[8][0][0]->SetMinimum(summary_int_min);
   Integral_up[8][0][0]->GetYaxis()->SetLabelSize(0.085);
   Integral_up[8][0][0]->GetYaxis()->SetTitleSize(0.085);
   Integral_up[8][0][0]->GetYaxis()->SetTitleOffset(0.7);
   Integral_up[8][0][0]->GetYaxis()->SetTitle("#Sigma P_{PbPb}- #Sigma P_{pp}   (GeV/c)");
    Integral_up[8][0][0]->GetYaxis()->CenterTitle();
   Integral_up[8][0][0]->GetYaxis()->SetTitleOffset(0.8);

   Integral_up[8][0][0]->GetXaxis()->SetLabelSize(0.0); 
   Integral_up[8][0][0]->GetXaxis()->SetRangeUser(0.501,7.99);
  
  

   Integral_down[8][0][0]->Draw("same ");


   Integral_syst[8][0][0]->Add( Integral_syst[10][0][0]);
   Integral_syst[8][0][0]->Draw("same e2");


   Integral[8][0][0]->Add(   Integral[10][0][0]);
   Integral[8][0][0]->Draw("same");
   

   l_int->Draw();

   legend_int_lead->Draw();
   

   
   aj_tex = new TLatex(0.2,0.9,"A_{J}<0.22");
   aj_tex->SetTextSize(0.08);
   aj_tex->SetLineColor(kWhite);
   aj_tex->SetNDC();
   aj_tex->Draw(); 
 

    l_int2->Draw();



   c_all_int->cd(2);

   Integral_up[8][0][1]->Draw();
  

   Integral_up[8][0][1]->SetMaximum(summary_int_max);
   Integral_up[8][0][1]->SetMinimum(summary_int_min);
   Integral_up[8][0][1]->GetYaxis()->SetLabelSize(0.0);
   Integral_up[8][0][1]->GetYaxis()->SetTitleSize(0.0);
 

   Integral_up[8][0][1]->GetXaxis()->SetLabelSize(0.0); 
   Integral_up[8][0][1]->GetXaxis()->SetRangeUser(0.501,7.99);
  
   Integral_down[8][0][1]->Draw("same ");


   Integral_syst[8][0][1]->Add(Integral_syst[10][0][1]);
   Integral_syst[8][0][1]->Draw("same e2");


   Integral[8][0][1]->Add( Integral[10][0][1]);
   Integral[8][0][1]->Draw("same");
  

   l_int->Draw();

   aj_tex = new TLatex(0.05,0.9,"A_{J}>0.22");
   aj_tex->SetTextSize(0.08);
   aj_tex->SetLineColor(kWhite);
   aj_tex->SetNDC();
   aj_tex->Draw(); 
  
  
   cms_tex->Draw(); 
   prelim_tex->Draw(); 




   cent_tex = new TLatex(0.4,0.7,"50-100% Peripheral");
   cent_tex->SetTextSize(0.08);
   cent_tex->SetLineColor(kWhite);
   cent_tex->SetNDC();
   cent_tex->Draw(); 
  

   data_tex->Draw(); 


   l_int2->Draw();
    c_all_int->cd(3);

   cout<<"here"<<endl;
   Integral_up[9][0][0]->Draw();



   Integral_up[9][0][0]->SetMaximum(summary_int_max2);
   Integral_up[9][0][0]->SetMinimum(summary_int_min2);
   Integral_up[9][0][0]->GetYaxis()->SetLabelSize(0.07);
   Integral_up[9][0][0]->GetYaxis()->SetTitleSize(0.07);
   Integral_up[9][0][0]->GetYaxis()->SetTitleOffset(.9);
   Integral_up[9][0][0]->GetYaxis()->SetTitle("#Sigma #DeltaP_{PbPb}- #Sigma #DeltaP_{pp}   (GeV/c)");
   Integral_up[9][0][0]->GetYaxis()->CenterTitle();
   Integral_up[9][0][0]->GetXaxis()->SetTitle("p_{T}^{assoc}   (GeV/c)");
   Integral_up[9][0][0]->GetXaxis()->SetTitleSize(0.08);
   Integral_up[9][0][0]->GetXaxis()->CenterTitle();
   Integral_up[9][0][0]->GetXaxis()->SetLabelSize(0.08); 
   Integral_up[9][0][0]->GetXaxis()->SetTitleOffset(1.);
   Integral_up[9][0][0]->GetXaxis()->SetRangeUser(.501,7.99);
   
   Integral_down[9][0][0]->Draw("same");


   Integral_syst[9][0][0]->Add(Integral_syst[6][0][0]);
   Integral_syst[9][0][0]->Draw("same e2");

   Integral[9][0][0]->Add(Integral[6][0][0]);
   Integral[9][0][0]->Draw("same");
     


   l_int->Draw();
   
   legend_int_diff->Draw();

   l_int2->Draw();
   cover_x_int->Draw();

    c_all_int->cd(4);
 
    Integral_up[9][0][1]->Draw();

   Integral_up[9][0][1]->SetMaximum(summary_int_max2);
   Integral_up[9][0][1]->SetMinimum(summary_int_min2);
   Integral_up[9][0][1]->GetYaxis()->SetLabelSize(0.0);
   Integral_up[9][0][1]->GetYaxis()->SetTitleSize(0.0);
    Integral_up[9][0][1]->GetXaxis()->SetTitle("p_{T}^{assoc}   (GeV/c)");
   Integral_up[9][0][1]->GetXaxis()->SetTitleSize(0.08);
   Integral_up[9][0][1]->GetXaxis()->CenterTitle();
   Integral_up[9][0][1]->GetXaxis()->SetLabelSize(0.08); 
   Integral_up[9][0][1]->GetYaxis()->SetTitleOffset(0.8);
   Integral_up[9][0][1]->GetXaxis()->SetRangeUser(.501,7.99);   Integral_up[9][0][1]->GetXaxis()->SetRangeUser(.501,7.99);
   
   Integral_down[9][0][1]->Draw("same");


   Integral_syst[9][0][1]->Add(Integral_syst[6][0][1]);
   Integral_syst[9][0][1]->Draw("same e2");


   Integral[9][0][1]->Add( Integral[6][0][1]);
   Integral[9][0][1]->Draw("same");
     



   l_int2->Draw();
   cover_x_int->Draw();


   c_all_int->cd(0);

   luminosity_tex_all_pp->Draw();
   
   luminosity_tex_all_PbPb->Draw();
    

     jet_reco_tex_all->Draw();

    jet_cut_tex_all->Draw();
    


   c_all_int->Print("Integral_PbPb_pp_Summary_Peripheral.png");
   c_all_int->Print("Integral_PbPb_pp_Summary_Peripheral.pdf");


   return 0;
  
}
