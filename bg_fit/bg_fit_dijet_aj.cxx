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

Int_t bg_fit_dijet_aj(bool is_number = kFALSE)
{

  TCanvas *dummy = new TCanvas("dummy");

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


  double phimin = -1.5;
  double phimax = 4.5;

  float x, y;
  TF1 *fit0 = new TF1("fit0","[0]",-3.,3.);

  Float_t v1_tot[12][6][4][5];
  Float_t v2_tot[12][6][4][5];
  Float_t v3_tot[12][6][4][5];

  TH1D *raw_eta[12][6][4][5];

  TH2D *yield[12][6][4][5];
  TH1D *yield_proj[12][6][4][5];
 
  TH1D *yield_rebin[12][6][4][5];
 
  TH1D *in_hist[12][6][4][5];
  TH1D *in_hist_rebin[12][6][4][5];
  TH1D *in_hist_ref[12][6][4][5];

  TH2D *in_hist2[12][6][4][5];
  TH2D *out_hist[12][6][4][5];

  TH1D *background_diff[12][6][4][5];
  TH1D *background_diff_rebin[12][6][4][5];
  TH1D *background_syst[12][6][4][5];

  TH2D *result[12][6][4][5];
  TH2D *result2[12][6][4][5];

  double error[12][6][4][5];
  double error_temp[12][6][4][5][8];

  TCanvas *cfit[12][5];
  TCanvas *cdiff_phi[12][5];
 	
  float yleft[12][6][4][5];
  float yright[12][6][4][5];
	

  int bin; 

  float dx_phi,nbins_sideband, nbins_peak;
  //----------------------------------------
  //   v2_jet[g][i][l] (for all centralities)
  //-------------------------------------------
  

  double bglevel, A_AS, V_2, V_1, V_3,alpha, beta, bglevelmax, bglevelmin, A_ASmax, A_ASmin, V_2max, V_1max, V_1min, V_3max, V_3min, V_2min, fitfunc0,fitmin0,fitmax0,histvalue, evalpt, error2, bglevel_err, V_1_err, V_2_err, V_3_err, A_AS_err, alpha_err, beta_err,bc, err, temp_bc, temp,temp_err, ymax, ymin, etamin1_val, etamax1_val, etamin2_val, etamax2_val, refvalue, referror, projmax,projmin, diff_min, diff_max;

  TString in_name, plotname, outname,bkgname2, resultname, centlabel, pTlabel,Ajlabel,checkcanvasnameEta,checkcanvasnamePhi,fitcanvasname, diff_etacanvasname, diff_phicanvasname,zoomedcanvasname, alphaCentcanvasname,alphapTcanvasname, betapTcanvasname, betaCentcanvasname;

  TString AjBin_strs[4] = {"Aj0","Aj22","Aj75","AjInclusive"};
  
  
  TLegend *lfit;

  TPaveText *pave;

  //------------------------------------------
  //Fit function and generic parameter sets
  //------------------------------------------

  TF1 *fourier = new TF1("fourier",
			 "[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))");
 
  fourier->SetParName(0,"Bkg level");
  fourier->SetParName(1,"V_{1}");
  fourier->SetParName(2,"V_{2}");
  fourier->SetParName(3,"V_{3}");
 
  TF1 *fourier_fit[12][6][4][3];
  TF1 *fourier_up[12][6][4][3];
  TF1 *fourier_down[12][6][4][3];
  TF1 *fourier_fit_diff[12][6][4][3];

 
  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
  
  TFile *f_PbPb = new TFile("../me_correct/PbPb_Correlations.root","READ");
  TFile *f_pp = new TFile("../me_correct/pp_Correlations.root","READ");

  TFile *f_out;
  if(is_number) f_out = new TFile("Dijet_Correlations_NoPtWeight.root","RECREATE");
  else f_out = new TFile("Dijet_Correlations.root","RECREATE");

  //---------------------------------
  //     Start i and l loops -- set up canvases etc.
  //------------------------------------------------


  for(int g = 2; g<6; g++){
    for(int l = 0; l<3; l++){

      fitcanvasname = "FitCanvas";
      fitcanvasname+= g;
      fitcanvasname+= l;
    
      cfit[g][l] = new TCanvas(fitcanvasname,"",0,0,1500,2000);
      cfit[g][l]->Divide(4,6,0.000001,0.000001);
    

      diff_phicanvasname= "DiffCanvas";
      diff_phicanvasname+= g;
      diff_phicanvasname+= l;
    
      cdiff_phi[g][l] = new TCanvas(diff_phicanvasname,"",0,0,1500,2000);
      cdiff_phi[g][l]->Divide(4,6,0.000,0.000);
    


      for(int i=0; i<6; i++){

	//----------------------------------------------------
	//  Start of j loop -- individual histo analysis begins
	//-----------------------------------------------------
	
	for (int j=0; j<4; j++){

	  cfit[g][l]->cd(4*i+j+1);
	  


	  //******************************************//

	  //  SET BACKGROUND SAMPLING REGION HERE//

	  //******************************************//

	  projmax = 2.5;
	  projmin = 1.5;


	  if(i==5){
	    projmax = 2.;
	    projmin = 1.;
	  }



	  etamin1_val = -projmax+0.001;
	  etamax1_val = -projmin-0.001;
	  etamin2_val =  projmin+0.001;
	  etamax2_val =  projmax -0.001;
	 
	  //******************************************//
	
	  //-------------------------------------
	  //  Get and assign all histograms 
	  //-------------------------------------

	  in_name = make_name("Yield_pTweighted_",g,i,j,l,pTlabel, centlabel, Ajlabel);
	  if(is_number) in_name = make_name("Yield_",g,i,j,l,pTlabel, centlabel, Ajlabel);

	  cout<<in_name<<endl;
	  if(g%2==0){
	    yield[g][i][j][l] = (TH2D*)f_PbPb->Get(in_name)->Clone(in_name);
	  }else{
	    yield[g][i][j][l] = (TH2D*)f_pp->Get(in_name)->Clone(in_name);
	  }

	  etalim = 1.5;

	  llimiteta = yield[g][i][j][l]->GetXaxis()->FindBin(-etalim+0.001);
	  rlimiteta = yield[g][i][j][l]->GetXaxis()->FindBin(etalim-0.001);
	
	  nbins_peak = (rlimiteta-llimiteta+1);


	  TString yieldprojname = in_name;
	  yieldprojname.ReplaceAll("Yield","Yield_proj");
	  yield_proj[g][i][j][l] = (TH1D*)yield[g][i][j][l]->ProjectionY(yieldprojname,llimiteta,rlimiteta);
	  yield_proj[g][i][j][l]->Scale(1./nbins_peak);
	
	  TString bkgname = in_name;
	  bkgname.ReplaceAll("Yield","Summed_bkg");

	  TString tempname = "temp_name";
	  tempname+=g;
	  tempname+=i;
	  tempname+=j;
	  tempname+=l;


	  Int_t etamin1 = yield[g][i][j][l]->GetXaxis()->FindBin(etamin1_val);
	  Int_t etamax1 = yield[g][i][j][l]->GetXaxis()->FindBin(etamax1_val);
	  Int_t etamin2 = yield[g][i][j][l]->GetXaxis()->FindBin(etamin2_val);
	  Int_t etamax2 = yield[g][i][j][l]->GetXaxis()->FindBin(etamax2_val);
 
	  nbins_sideband = etamax1-etamin1+etamax2-etamin2+2;
	
	  in_hist[g][i][j][l] = yield[g][i][j][l]->ProjectionY(bkgname,etamin1,etamax1);
	  in_hist[g][i][j][l]->Add(yield[g][i][j][l]->ProjectionY(tempname,etamin2,etamax2));
	  in_hist[g][i][j][l]->Scale(1./nbins_sideband);


	  // Deal properly with bad bins:

	  if(i<3){

	    for(int k = 1; k<in_hist[g][i][j][l]->GetNbinsX()+1; k++){
	    
	      temp_bc = in_hist[g][i][j][l]->GetBinContent(k+1);
	      temp_err = in_hist[g][i][j][l]->GetBinError(k+1);
	   
	   
	      if(in_hist[g][i][j][l]->GetBinError(k)<0.2*temp_err||in_hist[g][i][j][l]->GetBinError(k)>2.*temp_err){
		cout<<"******Replacing bin content**********"<<endl;
		cout<<g<<" "<<i<<" "<<j<<" "<<l<<" "<<in_hist[g][i][j][l]->GetBinContent(k)<<" "<<temp_bc<<" "<<in_hist[g][i][j][l]->GetBinError(k)<<" "<<temp_err<<endl;
		in_hist[g][i][j][l]->SetBinContent(k,temp_bc);
		in_hist[g][i][j][l]->SetBinError(k,temp_err);
	    
	      }
	
	   	   
	    }

	    if(g%2==0&&i==2&&j==0){

	  
	    for(int k = 1; k<in_hist[g][i][j][l]->GetNbinsX()/4+1; k++){
	    
	      temp_bc = in_hist[g][i][j][l]->GetBinContent(in_hist[g][i][j][l]->GetNbinsX()/2+1-k);
	      in_hist[g][i][j][l]->SetBinContent(k,temp_bc);
		    
	      }
	
	   	   
	    }



	  }

	    
	}
      }
    }

  }
  
  cout<<"ready to start fitting"<<endl;


  int g = 4;
  //-------------------------------------
  //      PbPb Fits
  //-------------------------------------
  for(int l = 0; l<3; l++){
      
    for(int j = 0; j<4; j++){

      for(int i = 0; i<6; i++){

	in_name = make_name("Yield_",g,i,j,l,pTlabel, centlabel, Ajlabel);

	in_name.ReplaceAll("Leading_","");

	cout<<"------------------------------------------------------------------------------------"<<endl;
	cout<<in_name<<endl;
	cout<<"------------------------------------------------------------------------------------"<<endl;
	




	TString diff_name = make_name("Diff_",0,i,j,l,pTlabel, centlabel, Ajlabel);  


	background_diff[0][i][j][l] = (TH1D*)in_hist[2][i][j][l]->Clone(diff_name);

	
	background_diff[0][i][j][l]->Add(in_hist[4][i][j][l],-1.);


	background_diff[0][i][j][l]->Scale(nbins_sideband);
	
	dx_phi = background_diff[0][i][j][l]->GetBinWidth(1);
	background_diff[0][i][j][l]->Scale(1./dx_phi);




	background_diff_rebin[0][i][j][l] = Rebin_dPhi(	background_diff[0][i][j][l]);

	background_diff_rebin[0][i][j][l]->SetAxisRange(-1.5,1.5);



	
	for(int k = 1; k<background_diff_rebin[0][i][j][l]->FindBin(TMath::Pi()/2.-.001)/2.+1; k++){

	  temp = (background_diff_rebin[0][i][j][l]->GetBinContent(k)+ background_diff_rebin[0][i][j][l]->GetBinContent(background_diff_rebin[0][i][j][l]->FindBin(TMath::Pi()/2.-.002)+1-k))/2.;

	  err = (background_diff_rebin[0][i][j][l]->GetBinError(k)+ background_diff_rebin[0][i][j][l]->GetBinError(background_diff_rebin[0][i][j][l]->FindBin(TMath::Pi()/2.-.002)+1-k))/2.;

	  background_diff_rebin[0][i][j][l]->SetBinContent(k,temp);
	  background_diff_rebin[0][i][j][l]->SetBinContent(background_diff_rebin[0][i][j][l]->FindBin(TMath::Pi()/2.-.002)+1-k,temp);
	
	  background_diff_rebin[0][i][j][l]->SetBinError(k,err);
	  background_diff_rebin[0][i][j][l]->SetBinError(background_diff_rebin[0][i][j][l]->FindBin(TMath::Pi()/2.-.002)+1-k,err);


	}



	background_diff_rebin[0][i][j][l]->SetMarkerSize(1);
	background_diff_rebin[0][i][j][l]->SetMarkerStyle(20);
	background_diff_rebin[0][i][j][l]->SetMarkerColor(kRed);
	background_diff_rebin[0][i][j][l]->SetLineColor(kRed);

	/*
	switch(i){
	case 0: 
	  diff_max = 8.;
	  diff_min = -3.;
	  break;
	case 1: 
	  diff_max = 2.;
	  diff_min = -1.;
	  break;
	case 2:
	  diff_max = 1.;
	  diff_min = -.5;
	  break;
	case 3: 
	  diff_max = 1.;
	  diff_min = -.5;
	  break;
	case 4: 
	  diff_max = 1.;
	  diff_min = -.5;
	  break;
	}
	*/

	diff_min = -3.;
	diff_max = 8.;

	background_diff_rebin[0][i][j][l]->SetMinimum(diff_min);
	background_diff_rebin[0][i][j][l]->SetMaximum(diff_max);

      
	if(j==0){
	  background_diff_rebin[0][i][j][l]->GetYaxis()->SetLabelSize(0.06);
	  background_diff_rebin[0][i][j][l]->GetYaxis()->SetTitleSize(0.06);
	  background_diff_rebin[0][i][j][l]->GetYaxis()->SetTitle("1/N_{evt} 1/d#Delta#phi");
	}else{
	  background_diff_rebin[0][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	}

	if(i==5){
	  background_diff_rebin[0][i][j][l]->GetXaxis()->SetLabelSize(0.06);
	  background_diff_rebin[0][i][j][l]->GetXaxis()->SetTitleSize(0.06);
	  background_diff_rebin[0][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	  background_diff_rebin[0][i][j][l]->GetXaxis()->CenterTitle();

	}	cout<<"and here 2"<<endl;

	TString diff_name2 = make_name("Diff_",1,i,j,l,pTlabel, centlabel, Ajlabel);

	background_diff[1][i][j][l] = (TH1D*)in_hist[3][i][j][l]->Clone(diff_name2);
	
	background_diff[1][i][j][l]->Add(in_hist[5][i][j][l],-1.);


	background_diff[1][i][j][l]->Scale(nbins_sideband);
	


	dx_phi = background_diff[1][i][j][l]->GetBinWidth(1);
	background_diff[1][i][j][l]->Scale(1./dx_phi);

	background_diff_rebin[1][i][j][l] = Rebin_dPhi(background_diff[1][i][j][l]);
	

	background_diff_rebin[1][i][j][l]->SetAxisRange(-1.5,1.5);

	
	for(int k = 1; k<background_diff_rebin[1][i][j][l]->FindBin(TMath::Pi()/2.-.001)/2.+1; k++){

	  temp = (background_diff_rebin[1][i][j][l]->GetBinContent(k)+ background_diff_rebin[1][i][j][l]->GetBinContent(background_diff_rebin[1][i][j][l]->FindBin(TMath::Pi()/2.-.002)+1-k))/2.;

	  err = (background_diff_rebin[1][i][j][l]->GetBinError(k)+ background_diff_rebin[1][i][j][l]->GetBinError(background_diff_rebin[1][i][j][l]->FindBin(TMath::Pi()/2.-.002)+1-k))/2.;

	  background_diff_rebin[1][i][j][l]->SetBinContent(k,temp);
	  background_diff_rebin[1][i][j][l]->SetBinContent(background_diff_rebin[1][i][j][l]->FindBin(TMath::Pi()/2.-.002)+1-k,temp);
	
	  background_diff_rebin[1][i][j][l]->SetBinError(k,err);
	  background_diff_rebin[1][i][j][l]->SetBinError(background_diff_rebin[1][i][j][l]->FindBin(TMath::Pi()/2.-.002)+1-k,err);


	}



	background_diff_rebin[1][i][j][l]->SetMarkerSize(1);
	background_diff_rebin[1][i][j][l]->SetMarkerStyle(10);
	background_diff_rebin[1][i][j][l]->SetMarkerColor(kBlack);
	background_diff_rebin[1][i][j][l]->SetLineColor(kBlack);

	nbins =  in_hist[4][i][j][l]->GetNbinsX();


	for(int k = 1; k<nbins/2+1; k++){

	  in_hist[4][i][j][l]->SetBinContent(k+nbins/2,in_hist[2][i][j][l]->GetBinContent(k));
	  in_hist[4][i][j][l]->SetBinError(k+nbins/2,in_hist[2][i][j][l]->GetBinError(k));
	 
	}

	cfit[4][l]->cd(4*i+j+1);
		
	fourier->ReleaseParameter(3);
	
	bglevel= in_hist[4][i][j][l]->GetBinContent(in_hist[4][i][j][l]->FindBin(TMath::Pi()/2));
	fourier->SetParameter(0,0.);
	fourier->SetParameter(1,bglevel);
	fourier->SetParameter(2,0.03);
	fourier->SetParameter(3,0.);
	fourier->SetParLimits(2,0.,1.0);
	fourier->SetParLimits(1,-1.,1.0);
	fourier->SetParLimits(3,0.,1.0);


	V_3 = 0.;
	V_2 = 0.;
	V_1 = 0.;
	
	if(i>3)fourier->FixParameter(3,0.);

	in_hist[4][i][j][l]->Fit("fourier","","",-TMath::Pi()/2.,3.*TMath::Pi()/2.);


	bglevel = fourier->GetParameter(0);
	V_1     = fourier->GetParameter(1);
	V_2     = fourier->GetParameter(2);
	V_3     = fourier->GetParameter(3);

	if(abs(V_1)<1e-6)V_1=0.;
	if(abs(V_2)<1e-6)V_2=0.;
	if(abs(V_3)<1e-6)V_3=0.;
	if(abs(V_1_err)<1e-6)V_1_err=0.;
	if(abs(V_2_err)<1e-6)V_2_err=0.;
	if(abs(V_3_err)<1e-6)V_3_err=0.;


	bglevel_err = fourier->GetParError(0);
	V_1_err = fourier->GetParError(1);
	V_2_err = fourier->GetParError(2);
	V_3_err = fourier->GetParError(3);



	v1_tot[4][i][j][l] = V_1;
	v2_tot[4][i][j][l] = V_2;
	v3_tot[4][i][j][l] = V_3;



	fourier_fit_diff[4][i][j][l] = new TF1(Form("FourierFitDiffPbPb_lead%d%d%d",i,j,l),"[0]*(2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))",-TMath::Pi()/2.,TMath::Pi()/2.);

	fourier_fit[4][i][j][l] = new TF1(Form("FourierFitPbPb_lead%d%d%d",i,j,l),"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))",-TMath::Pi()/2.,3.*TMath::Pi()/2.);

	fourier_up[4][i][j][l] = new TF1(Form("FourierFitPbPb_up%d%d%d",i,j,l),"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))",-TMath::Pi()/2.,3.*TMath::Pi()/2.);

	fourier_down[4][i][j][l] = new TF1(Form("FourierFitPbPb_down%d%d%d",i,j,l),"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))",-TMath::Pi()/2.,3.*TMath::Pi()/2.);

	fourier_fit[4][i][j][l]->SetParameter(1,V_1);
	fourier_fit[4][i][j][l]->SetParameter(3,V_3);
	fourier_fit[4][i][j][l]->SetParameter(2,V_2);
	fourier_fit[4][i][j][l]->SetParameter(0,bglevel*nbins_peak/dx_phi);
	fourier_fit[4][i][j][l]->SetLineColor(kRed);
	fourier_fit[4][i][j][l]->SetLineWidth(2);

	float	temp_error=0.;	

	
	fourier->SetParameter(0,bglevel*nbins_peak/dx_phi);

	fourier->SetParameter(1,V_1+V_1_err);
	temp_error+=(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.));

	fourier->SetParameter(1,V_1-V_1_err);
	temp_error+=(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.));

	fourier->SetParameter(1,V_1);

	fourier->SetParameter(2,V_2+V_2_err);
	temp_error+=(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.));

	fourier->SetParameter(2,V_2-V_2_err);
	temp_error+=(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.));

	fourier->SetParameter(2,V_2);

	fourier->SetParameter(3,V_3+V_3_err);
	temp_error+=(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.));

	fourier->SetParameter(3,V_3-V_3_err);
	temp_error+=(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.));

	fourier->SetParameter(3,V_3);

	fourier->SetParameter(0,(bglevel+bglevel_err)*nbins_peak/dx_phi);
	temp_error+=(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.));

	fourier->SetParameter(0,(bglevel-bglevel_err)*nbins_peak/dx_phi);
	temp_error+=(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[4][i][j][l]->Eval(0.));

	fourier->SetParameter(0,bglevel);


	error[4][i][j][l]=TMath::Sqrt(temp_error/8.);

	
	fourier_up[4][i][j][l]->SetLineColor(kBlue);
	fourier_up[4][i][j][l]->SetLineWidth(4);
	fourier_down[4][i][j][l]->SetLineColor(kBlue);
	fourier_down[4][i][j][l]->SetLineWidth(4);
	


	
	fourier_up[4][i][j][l]->SetParameter(1,V_1);
	fourier_up[4][i][j][l]->SetParameter(3,V_3);
	fourier_up[4][i][j][l]->SetParameter(2,V_2);
	fourier_up[4][i][j][l]->SetParameter(0,bglevel*nbins_peak/dx_phi+error[4][i][j][l]);
	
	fourier_down[4][i][j][l]->SetParameter(1,V_1);
	fourier_down[4][i][j][l]->SetParameter(3,V_3);
	fourier_down[4][i][j][l]->SetParameter(2,V_2);
	fourier_down[4][i][j][l]->SetParameter(0,bglevel*nbins_peak/dx_phi-error[4][i][j][l]);

	TString syst_name = make_name("Syst_",0,i,j,l,pTlabel, centlabel, Ajlabel);  

	background_syst[0][i][j][l] = (TH1D*)background_diff_rebin[0][i][j][l]->Clone(syst_name);

	for(int k = 0; k<background_syst[0][i][j][l]->GetNbinsX()+1; k++){
	  err = error[4][i][j][l];
	  background_syst[0][i][j][l]->SetBinError(k,err);
	}

	//---------------------------------------------
	//    Fill PbPb new background
	//---------------------------------------------

	outname= in_name;
	outname.ReplaceAll("Yield","Fit_bkg");    
	out_hist[4][i][j][l] = (TH2D*)yield[g][i][j][l]->Clone(outname);
	
	bkgname2 = in_name;
	bkgname2.ReplaceAll("Yield","Summed2D_bkg");
	in_hist2[4][i][j][l] = (TH2D*)yield[g][i][j][l]->Clone(bkgname2);

	cout<<"we're going to replace "<<out_hist[4][i][j][l]->GetNbinsX()<<" by "<<out_hist[4][i][j][l]->GetNbinsY()<<" bins"<<endl;
	

	nbins = out_hist[4][i][j][l]->GetNbinsY();

	for(int k = 1; k<out_hist[4][i][j][l]->GetNbinsY()+1;k++){
   
	  evalpt = in_hist[4][i][j][l]->GetBinCenter(k);
	  histvalue = fourier->Eval(evalpt);
	  error2 = error[4][i][j][l];

	  refvalue = in_hist[4][i][j][l]->GetBinContent(k);
	  referror = in_hist[4][i][j][l]->GetBinError(k);

	  for(int m = 1; m<out_hist[4][i][j][l]->GetNbinsX()+1; m++){
   
	    in_hist2[4][i][j][l]->SetBinContent(m,k,refvalue);
	    in_hist2[4][i][j][l]->SetBinError(m,k,referror*sqrt(nbins_sideband));
	 
	    out_hist[4][i][j][l]->SetBinContent(m,k,histvalue);
	    out_hist[4][i][j][l]->SetBinError(m,k,referror*sqrt(nbins_sideband));

	    if(k<nbins/2){
	      yield[4][i][j][l]->SetBinContent(m,k+nbins/2,yield[2][i][j][l]->GetBinContent(m,k));
	      yield[4][i][j][l]->SetBinError(m,k+nbins/2,yield[2][i][j][l]->GetBinError(m,k));
	    }
	  }
	}

	resultname = in_name;
	resultname.ReplaceAll("Yield","Yield_BkgSub");

	result[4][i][j][l] = (TH2D*)yield[4][i][j][l]->Clone(resultname);
	/*
	if(i<4){
	  result[4][i][j][l]->Add(out_hist[4][i][j][l],-1.);
	}else{
	*/
	result[4][i][j][l]->Add(in_hist2[4][i][j][l],-1.);
	//}


	//--------------------------------------------
	//   Draw fit plots
	//---------------------------------------------
      

	if(i<6){

	  cfit[4][l]->cd(4*i+j+1);

	  TString outhistproj_name = in_name;
	  outhistproj_name.ReplaceAll("Yield","New1d_bkg");

	  for(int k = 1; k<nbins/2+1; k++){

	    yield_proj[4][i][j][l]->SetBinContent(k+nbins/2,yield_proj[2][i][j][l]->GetBinContent(k));
	    yield_proj[4][i][j][l]->SetBinError(k+nbins/2,yield_proj[2][i][j][l]->GetBinError(k));
	 
	  }
	
	  in_hist_rebin[4][i][j][l] = Rebin_dPhi_full(in_hist[4][i][j][l]);
          
	  in_hist_rebin[4][i][j][l]->Scale(1./dx_phi*nbins_peak);
	  in_hist_rebin[4][i][j][l]->SetMarkerStyle(20);
	  in_hist_rebin[4][i][j][l]->SetMarkerSize(1);
	  in_hist_rebin[4][i][j][l]->SetMarkerColor(kBlack);
	  in_hist_rebin[4][i][j][l]->SetLineColor(kBlack);


	  //	in_hist_rebin[4][i][j][l]->SetAxisRange(-1.5,4.8);

	  in_hist_rebin[4][i][j][l]->SetMaximum(ymax);
	  in_hist_rebin[4][i][j][l]->SetMinimum(ymin);
	  in_hist_rebin[4][i][j][l]->Draw();
	
	
	  yield_rebin[4][i][j][l]= Rebin_dPhi_full(yield_proj[4][i][j][l]);
   

	  TString yield_rebin_name = in_name;
	  yield_rebin_name += "_1D_rebin";
	  yield_rebin_name +=i;
	  yield_rebin_name+=j;
	  yield_rebin[4][i][j][l]->SetName(yield_rebin_name);
	


	  yield_rebin[4][i][j][l]->SetMarkerColor(kViolet-5);

	  yield_rebin[4][i][j][l]->SetLineColor(kViolet-5);
	  yield_rebin[4][i][j][l]->SetMarkerStyle(10);

	  yield_rebin[4][i][j][l]->Scale(1.*nbins_peak/dx_phi);

	  ymin =yield_rebin[4][i][j][l]->GetMinimum()-1;
	  ymax =yield_rebin[4][i][j][l]->GetMaximum()+2.;

	  if((in_hist_rebin[4][i][j][l]->GetMaximum()+.004)>ymax){ymax = in_hist_rebin[4][i][j][l]->GetMaximum()+.004;}

	  if(i==0)ymax=ymin+9.;

	  yield_rebin[4][i][j][l]->SetMaximum(ymax);
	  yield_rebin[4][i][j][l]->SetMinimum(ymin);


	  yield_rebin[4][i][j][l]->GetXaxis()->SetLabelSize(0.07);
	  yield_rebin[4][i][j][l]->GetXaxis()->SetTitleSize(0.07);
	  yield_rebin[4][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	  yield_rebin[4][i][j][l]->GetXaxis()->CenterTitle();
	  yield_rebin[4][i][j][l]->GetYaxis()->SetLabelSize(0.07);
	  yield_rebin[4][i][j][l]->GetYaxis()->SetTitleSize(0.07);
	  yield_rebin[4][i][j][l]->GetYaxis()->SetTitle("1/N_{evt} d p_{T}/d#Delta#phi");
	  yield_rebin[4][i][j][l]->GetYaxis()->CenterTitle();

	  yield_rebin[4][i][j][l]->Draw();

	  in_hist_rebin[4][i][j][l]->Draw("same");
	  fourier_fit[4][i][j][l]->Draw("same");
	
	  fourier_up[4][i][j][l]->Draw("same");
	  fourier_down[4][i][j][l]->Draw("same");
	
	  lfit = new TLegend(0.65,0.8,0.88,0.94);
	  lfit->SetFillColor(kWhite);
	  lfit->SetLineColor(kWhite);
	  lfit->AddEntry(yield_rebin[4][i][j][l],"Yield");
	  lfit->AddEntry(in_hist_rebin[4][i][j][l],"Sum. Bkg");
	  lfit->SetTextSize(ts2-0.01);
	  lfit->Draw("same");


	  lfit->Draw("same");


	  pave = new TPaveText(0.18,0.75,0.45,0.9,"NDC");

	  pave->SetName("pave");
	  pave->SetFillColor(0);
	  pave->SetLineColor(0);
	  pave->SetTextAlign(11);
	  pave->AddText(centlabel);
	  pave->AddText(pTlabel);
	  pave->AddText(Ajlabel);
	  pave->SetTextSize(ts2-0.01);
	  pave->Draw("same");

	  TLine *pi_over_2 = new TLine(TMath::Pi()/2.,ymin, TMath::Pi()/2.,ymax);
	  pi_over_2->SetLineStyle(2);
	  pi_over_2->Draw();

	}

	//-------------------
  
	//  Now for pp fits

	//-------------------

	in_name.ReplaceAll("PbPb","pp");

	for(int k = 1; k<nbins/2+1; k++){

	  in_hist[5][i][j][l]->SetBinContent(k+nbins/2,in_hist[3][i][j][l]->GetBinContent(k));
	  in_hist[5][i][j][l]->SetBinError(k+nbins/2,in_hist[3][i][j][l]->GetBinError(k));
	 
	}
      


	cfit[5][l]->cd(4*i+j+1);
		

	bglevel= in_hist[5][i][j][l]->GetBinContent(in_hist[5][i][j][l]->FindBin(TMath::Pi()/2));

	fourier->SetParameter(0,bglevel);
	fourier->SetParLimits(2,0.,1.0);
	fourier->SetParLimits(1,-1.,1.0);
	fourier->SetParLimits(3,0.,1.0);



	V_3= 0.;
	V_2 = 0.;
	V_1 = 0.;
	
	in_hist[5][i][j][l]->Fit("fourier","","",-TMath::Pi()/2.,3.*TMath::Pi()/2.);

	bglevel = fourier->GetParameter(0);
	V_1     = fourier->GetParameter(1);
	V_2     = fourier->GetParameter(2);
	V_3     = fourier->GetParameter(3);


	bglevel_err = fourier->GetParError(0);
	V_1_err = fourier->GetParError(1);
	V_2_err = fourier->GetParError(2);
	V_3_err = fourier->GetParError(3);


	/*

	if(abs(V_1)<1e-6)V_1=0.;
	if(abs(V_2)<1e-6)V_2=0.;
	if(abs(V_3)<1e-6)V_3=0.;

	if(abs(V_1_err)<1e-6)V_1_err=0.;
	if(abs(V_2_err)<1e-6)V_2_err=0.;
	if(abs(V_3_err)<1e-6)V_3_err=0.;
	*/

	v1_tot[5][i][j][l] = V_1;
	v2_tot[5][i][j][l] = V_2;
	v3_tot[5][i][j][l] = V_3;


	fourier_fit[5][i][j][l] = new TF1(Form("FourierFitpp_lead%d%d%d",i,j,l),"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))",-TMath::Pi()/2.,3.*TMath::Pi()/2.);

	fourier_up[5][i][j][l] = new TF1(Form("FourierFitPbPb_up%d%d%d",i,j,l),"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))",-TMath::Pi()/2.,3.*TMath::Pi()/2.);

	fourier_down[5][i][j][l] = new TF1(Form("FourierFitPbPb_down%d%d%d",i,j,l),"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))",-TMath::Pi()/2.,3.*TMath::Pi()/2.);

	fourier_fit[5][i][j][l]->SetParameter(1,V_1);
	fourier_fit[5][i][j][l]->SetParameter(3,V_3);
	fourier_fit[5][i][j][l]->SetParameter(2,V_2);
	fourier_fit[5][i][j][l]->SetParameter(0,bglevel*nbins_peak/dx_phi);
	fourier_fit[5][i][j][l]->SetLineColor(kRed);
	fourier_fit[5][i][j][l]->SetLineWidth(2);

	error[5][i][j][l] = 0.;


	fourier->SetParameter(0,bglevel*nbins_peak/dx_phi);

	fourier->SetParameter(1,V_1+V_1_err);
	error[5][i][j][l]+=(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.));

	fourier->SetParameter(1,V_1-V_1_err);
	error[5][i][j][l]+=(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.));

	fourier->SetParameter(1,V_1);

	fourier->SetParameter(2,V_2+V_2_err);
	error[5][i][j][l]+=(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.));

	fourier->SetParameter(2,V_2-V_2_err);
	error[5][i][j][l]+=(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.));

	fourier->SetParameter(2,V_2);

	fourier->SetParameter(3,V_3+V_3_err);
	error[5][i][j][l]+=(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.));

	fourier->SetParameter(3,V_3-V_3_err);
	error[5][i][j][l]+=(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.));

	fourier->SetParameter(3,V_3);

	fourier->SetParameter(0,(bglevel+bglevel_err)*nbins_peak/dx_phi);
	error[5][i][j][l]+=(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.));

	fourier->SetParameter(0,(bglevel-bglevel_err)*nbins_peak/dx_phi);
	error[5][i][j][l]+=(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.))*(fourier->Eval(0.)-fourier_fit[5][i][j][l]->Eval(0.));

	fourier->SetParameter(0,bglevel);

	error[5][i][j][l]=TMath::Sqrt(error[5][i][j][l]/8.);
	
	fourier_up[5][i][j][l]->SetLineColor(kBlue);
	fourier_up[5][i][j][l]->SetLineWidth(4);
	fourier_down[5][i][j][l]->SetLineColor(kBlue);
	fourier_down[5][i][j][l]->SetLineWidth(4);
	

	fourier_up[5][i][j][l]->SetParameter(1,V_1);
	fourier_up[5][i][j][l]->SetParameter(3,V_3);
	fourier_up[5][i][j][l]->SetParameter(2,V_2);
	fourier_up[5][i][j][l]->SetParameter(0,bglevel*nbins_peak/dx_phi+error[5][i][j][l]);
	
	fourier_down[5][i][j][l]->SetParameter(1,V_1);
	fourier_down[5][i][j][l]->SetParameter(3,V_3);
	fourier_down[5][i][j][l]->SetParameter(2,V_2);
	fourier_down[5][i][j][l]->SetParameter(0,bglevel*nbins_peak/dx_phi-error[5][i][j][l]);



	
	 syst_name = make_name("Syst_",1,i,j,l,pTlabel, centlabel, Ajlabel);  

	background_syst[1][i][j][l] = (TH1D*)background_diff_rebin[1][i][j][l]->Clone(syst_name);


	for(int k = 0; k<background_syst[1][i][j][l]->GetNbinsX()+1; k++){
	  err = error[5][i][j][l];
	  background_syst[1][i][j][l]->SetBinError(k,err);
	}





	//---------------------------------------------
	//    Fill and subtract pp new background
	//---------------------------------------------

	outname= in_name;
	outname.ReplaceAll("Yield","Fit_bkg");    
	out_hist[5][i][j][l] = (TH2D*)yield[5][i][j][l]->Clone(outname);
	
	bkgname2 = in_name;
	bkgname2.ReplaceAll("Yield","Summed2D_bkg");
	in_hist2[5][i][j][l] = (TH2D*)yield[5][i][j][l]->Clone(bkgname2);


	for(int k = 1; k<out_hist[5][i][j][l]->GetNbinsY()+1;k++){
   
	  evalpt = in_hist[5][i][j][l]->GetBinCenter(k);
	  histvalue = fourier->Eval(evalpt);
	  error2 = error[5][i][j][l];

	  refvalue = in_hist[5][i][j][l]->GetBinContent(k);
	  referror = in_hist[5][i][j][l]->GetBinError(k);
	

	  for(int m = 1; m<out_hist[5][i][j][l]->GetNbinsX(); m++){
   
	    in_hist2[5][i][j][l]->SetBinContent(m,k,refvalue);
	    in_hist2[5][i][j][l]->SetBinError(m,k,referror*sqrt(nbins_sideband));
	 
	    out_hist[5][i][j][l]->SetBinContent(m,k,histvalue);
	    out_hist[5][i][j][l]->SetBinError(m,k,referror*sqrt(nbins_sideband));

	    if(k<nbins/2){
	      yield[5][i][j][l]->SetBinContent(m,k+nbins/2,yield[3][i][j][l]->GetBinContent(m,k));
	      yield[5][i][j][l]->SetBinError(m,k+nbins/2,yield[3][i][j][l]->GetBinError(m,k));
	    }

	  }
	}

	resultname = in_name;
	resultname.ReplaceAll("Yield","Yield_BkgSub");

	result[5][i][j][l] = (TH2D*)yield[5][i][j][l]->Clone(resultname);
	/*
	if(i<4){
	  result[5][i][j][l]->Add(out_hist[5][i][j][l],-1.);
	}else{
	*/
	result[5][i][j][l]->Add(in_hist2[5][i][j][l],-1.);
	//}
     
	cfit[5][l]->cd(4*i+j+1);



	TString outhistproj_name = in_name;
	outhistproj_name.ReplaceAll("Yield","New1d_bkg");

	for(int k = 1; k<nbins/2+1; k++){

	  yield_proj[5][i][j][l]->SetBinContent(k+nbins/2,yield_proj[3][i][j][l]->GetBinContent(k));
	  yield_proj[5][i][j][l]->SetBinError(k+nbins/2,yield_proj[3][i][j][l]->GetBinError(k));
	 
	}
	

	in_hist_rebin[5][i][j][l] = Rebin_dPhi_full(in_hist[5][i][j][l]);
          
	in_hist_rebin[5][i][j][l]->Scale(nbins_peak/dx_phi);
	in_hist_rebin[5][i][j][l]->SetMarkerStyle(20);
	in_hist_rebin[5][i][j][l]->SetMarkerSize(1);
	in_hist_rebin[5][i][j][l]->SetMarkerColor(kBlack);
	in_hist_rebin[5][i][j][l]->SetLineColor(kBlack);



	yield_rebin[5][i][j][l]= Rebin_dPhi_full(yield_proj[5][i][j][l]);
   

	TString yield_rebin_name = in_name;
	yield_rebin_name += "_1D_rebin";
	yield_rebin_name +=i;
	yield_rebin_name+=j;
	yield_rebin[5][i][j][l]->SetName(yield_rebin_name);
	


	yield_rebin[5][i][j][l]->SetMarkerColor(kViolet-5);

	yield_rebin[5][i][j][l]->SetLineColor(kViolet-5);
	yield_rebin[5][i][j][l]->SetMarkerStyle(10);




	yield_rebin[5][i][j][l]->Scale(1.*nbins_peak/dx_phi);


	ymin =yield_rebin[5][i][j][l]->GetMinimum()-1.;
	ymax =yield_rebin[5][i][j][l]->GetMaximum()+2.;
	if((in_hist_rebin[5][i][j][l]->GetMaximum()+.004)>ymax){ymax = in_hist_rebin[5][i][j][l]->GetMaximum()+.004;}

	yield_rebin[5][i][j][l]->SetMaximum(ymax);
	yield_rebin[5][i][j][l]->SetMinimum(ymin);

	yield_rebin[5][i][j][l]->GetXaxis()->SetLabelSize(0.07);
	yield_rebin[5][i][j][l]->GetXaxis()->SetTitleSize(0.07);
	yield_rebin[5][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	yield_rebin[5][i][j][l]->GetXaxis()->CenterTitle();
	yield_rebin[5][i][j][l]->GetYaxis()->SetLabelSize(0.07);
	yield_rebin[5][i][j][l]->GetYaxis()->SetTitleSize(0.07);
	yield_rebin[5][i][j][l]->GetYaxis()->SetTitle("1/N_{evt} d p_{T}/d#Delta#phi");
	yield_rebin[5][i][j][l]->GetYaxis()->CenterTitle();


	yield_rebin[5][i][j][l]->Draw();
	in_hist_rebin[5][i][j][l]->Draw("same");
	/*
	fourier_up[5][i][j][l]->Draw("same");
	fourier_down[5][i][j][l]->Draw("same");
	*/
	fourier_fit[5][i][j][l]->Draw("same");
	
	lfit->Draw("same");
	pave->Draw("same");

	TLine *pi_over_2 = new TLine(TMath::Pi()/2.,ymin, TMath::Pi()/2.,ymax);
	pi_over_2->SetLineStyle(2);
	pi_over_2->Draw();
	pi_over_2->Draw();

	cdiff_phi[g][l]->cd(4*i+j+1);
	
	background_diff_rebin[0][i][j][l]->Draw();
	background_diff_rebin[1][i][j][l]->Draw("same");

	fourier_fit_diff[4][i][j][l]->Draw("same");
	

	TPaveText *pave2;
	if(j==0){
	  pave2= new TPaveText(0.18,0.75,0.45,0.9,"NDC");
	}else{
	    pave2= new TPaveText(0.05,0.75,0.45,0.9,"NDC");
	}
	pave2->SetName("pave2");
	pave2->SetFillColor(0);
	pave2->SetLineColor(0);
	pave2->SetTextAlign(11);
	pave2->AddText(centlabel);
	pave2->AddText(pTlabel);
	pave2->AddText(Ajlabel);
	pave2->SetTextSize(ts2-0.01);
	pave2->Draw("same");
      
	TLine *lphi = new TLine(-TMath::Pi()/2.,0.,TMath::Pi()/2.,0.);
	lphi->SetLineStyle(2);
	lphi->Draw();

 


	float diff_val = background_diff_rebin[0][i][j][l]->Integral(1,background_diff_rebin[0][i][j][l]->GetNbinsX(),"width");
    
	TString int_label = Form("PbPb Int. =%.2f",diff_val);  
 
	TLatex *int_tex = new TLatex(0.55,0.9,int_label);
	int_tex->SetTextSize(ts2-0.01);
	int_tex->SetNDC();
	int_tex->SetTextColor(kRed);
	int_tex->Draw();


	TString int_label2 = Form("Over |#Delta#eta|<2.5: %.2f",diff_val*5./2.);

	TLatex *int_tex2 = new TLatex(0.55,0.85,int_label2);
	int_tex2->SetTextSize(ts2-0.01);
	int_tex2->SetNDC();
	int_tex2->SetTextColor(kRed);
	int_tex2->Draw();


      

	float diff_val3 = background_diff_rebin[1][i][j][l]->Integral(1,background_diff_rebin[0][i][j][l]->GetNbinsX(),"width");
   
	TString int_label3 = Form("pp Int. = %.2f",diff_val3);

	TLatex *int_tex3 = new TLatex(0.55,0.75,int_label3);
	int_tex3->SetTextSize(ts2-0.01);
	int_tex3->SetNDC();
	int_tex3->Draw();


	TString int_label4 = Form("Over |#Delta#eta|<2.5: %.2f",diff_val3*5/2.);


	TLatex *int_tex4 = new TLatex(0.55,0.7,int_label4);
	int_tex4->SetTextSize(ts2-0.01);
	int_tex4->SetNDC();
	int_tex4->Draw();

	dummy->cd();
      

	yield[4][i][j][l]->Write();
	out_hist[4][i][j][l]->Write();
	in_hist2[4][i][j][l]->Write();
	result[4][i][j][l]->Write();
	background_diff_rebin[0][i][j][l]->Write();
	background_syst[0][i][j][l]->Write();

	yield[5][i][j][l]->Write();
	out_hist[5][i][j][l]->Write();
	in_hist2[5][i][j][l]->Write();
	result[5][i][j][l]->Write();
	background_diff_rebin[1][i][j][l]->Write();
	background_syst[1][i][j][l]->Write();


      } //closes the [i] (centrality) loop

    

    } //closes the [j] (pT) loop
    cout<<"about to save"<<endl;

    if(l<2){
    cfit[4][l]->SaveAs((TString)("DijetFreeFits_PbPb_GluedBackgroundV1V2V3"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf"));
    cfit[4][l]->SaveAs((TString)("DijetFreeFits_PbPb_GluedBackgroundV1V2V3"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));
    cfit[5][l]->SaveAs((TString)("DijetFreeFits_pp_GluedBackgroundV1V2V3"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf"));
    cfit[5][l]->SaveAs((TString)("DijetFreeFits_pp_GluedBackgroundV1V2V3"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));

    cdiff_phi[4][l]->SaveAs((TString)("BackgroundDifference_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png"));
    cdiff_phi[4][l]->SaveAs((TString)("BackgroundDifference_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf"));
   
    }else{
      cfit[4][l]->SaveAs((TString)("DijetFreeFits_PbPb_GluedBackgroundV1V2V3_AjInclusive.pdf"));
      cfit[4][l]->SaveAs((TString)("DijetFreeFits_PbPb_GluedBackgroundV1V2V3_AjInclusive.png"));
      
      cfit[5][l]->SaveAs((TString)("DijetFreeFits_pp_GluedBackgroundV1V2V3_AjInclusive.pdf"));
      cfit[5][l]->SaveAs((TString)("DijetFreeFits_pp_GluedBackgroundV1V2V3_AjInclusive.png"));

      cdiff_phi[4][l]->SaveAs("BackgroundDifference_AjInclusive.png");
      cdiff_phi[4][l]->SaveAs("BackgroundDifference_AjInclusive.pdf");
    }
    

  }

  TCanvas *c_demo = new TCanvas("Decomposition Display","",10,10,1200,800);
  c_demo->Divide(3,2,0.0001,0.0001);

  c_demo->cd(1);
  
  yield_rebin[5][1][3][1]->SetMinimum(0.1);
  yield_rebin[5][1][3][1]->SetMaximum(22.);

  yield_rebin[5][1][3][1]->GetYaxis()->SetTitleOffset(0.8);
  yield_rebin[5][1][3][1]->Draw();
  in_hist_rebin[5][1][3][1]->Draw("same");
  in_hist_rebin[5][1][3][1]->Rebin(2);
  in_hist_rebin[5][1][3][1]->Scale(1./2);
  in_hist_rebin[5][1][3][1]->Draw("same");
 
  /* fourier_up[5][1][3][1]->Draw("same");
  fourier_down[5][1][3][1]->Draw("same");
  fourier_fit[5][1][3][1]->Draw("same");
  */
  yield_rebin[5][1][3][1]->Draw("axis same");

  TLatex *label_pp = new TLatex(0.2,0.88,"pp Reference");
  label_pp->SetTextSize(0.07);
  label_pp->SetLineColor(kWhite);
  label_pp->SetNDC();
  label_pp->Draw();

  TLine *pi_over_2 = new TLine(TMath::Pi()/2.,0.1, TMath::Pi()/2.,22.);
  pi_over_2->SetLineStyle(2);
  pi_over_2->Draw();

  TLegend *legend = new TLegend(0.17,0.73,0.9,0.83);
  legend->AddEntry(yield_rebin[5][1][3][1],"Dijet Peaks, |#Delta#eta|<1.5");
  legend->AddEntry(in_hist_rebin[5][1][3][1],"Long Range Dist., |#Delta#eta|<2.5");
  //legend->AddEntry(fourier_fit[5][1][3][1],"Fit to Long Range Dist.","l");
  legend->SetTextSize(0.05);
  legend->SetLineColor(kWhite);
  legend->Draw();


  c_demo->cd(2);

  yield_rebin[4][1][0][1]->SetMinimum(7.);
  yield_rebin[4][1][0][1]->SetMaximum(29.);
 

  yield_rebin[4][1][0][1]->Draw();
  in_hist_rebin[4][1][0][1]->Rebin(2);
  in_hist_rebin[4][1][0][1]->Scale(1./2);
 in_hist_rebin[4][1][0][1]->Draw("same");
 /*
  fourier_up[4][1][0][1]->Draw("same");
  fourier_down[4][1][0][1]->Draw("same");
  fourier_fit[4][1][0][1]->Draw("same");
 */
  TLatex *label_per = new TLatex(0.2,0.88,"PbPb Cent. 50-100%");
  label_per->SetTextSize(0.07);
  label_per->SetLineColor(kWhite);
  label_per->SetNDC();
  label_per->Draw();

  pi_over_2 = new TLine(TMath::Pi()/2.,7., TMath::Pi()/2.,29.);
  pi_over_2->SetLineStyle(2);
  pi_over_2->Draw();


  TLatex *label_pT = new TLatex(0.2,0.8,"1<p_{T}^{assoc.}<2 GeV/c");
  label_pT->SetTextSize(0.07);
  label_pT->SetLineColor(kWhite);
  label_pT->SetNDC();
  label_pT->Draw();



  c_demo->cd(3);

  yield_rebin[4][1][3][1]->SetMinimum(140.);
  yield_rebin[4][1][3][1]->SetMaximum(162.);

  yield_rebin[4][1][3][1]->Draw();
  in_hist_rebin[4][1][3][1]->Rebin(2);
  in_hist_rebin[4][1][3][1]->Scale(1./2);
  in_hist_rebin[4][1][3][1]->Draw("same");
  /*
  fourier_up[4][1][3][1]->Draw("same");
  fourier_down[4][1][3][1]->Draw("same");
  fourier_fit[4][1][3][1]->Draw("same");
  */  
  TLatex *label_cent = new TLatex(0.2,0.88,"PbPb Cent. 0-10%");
  label_cent->SetTextSize(0.07);
  label_cent->SetLineColor(kWhite);
  label_cent->SetNDC();
  label_cent->Draw();

  pi_over_2 = new TLine(TMath::Pi()/2.,140., TMath::Pi()/2.,162.);
  pi_over_2->SetLineStyle(2);
  pi_over_2->Draw();

  TLatex *label_aj = new TLatex(0.2,0.78,"A_{J}>0.22");
  label_aj->SetTextSize(0.07);
  label_aj->SetLineColor(kWhite);
  label_aj->SetNDC();
  label_aj->Draw();




  c_demo->cd(4);
  
  yield_rebin[5][4][3][1]->SetMinimum(0.1);
  yield_rebin[5][4][3][1]->SetMaximum(47.);

  yield_rebin[5][4][3][1]->GetYaxis()->SetTitleOffset(0.8);
  yield_rebin[5][4][3][1]->Draw();
  in_hist_rebin[5][4][3][1]->Draw("same");
  in_hist_rebin[5][4][3][1]->Rebin(2);
  in_hist_rebin[5][4][3][1]->Scale(1./2);
  in_hist_rebin[5][4][3][1]->Draw("same");
 
  /* fourier_up[5][1][3][1]->Draw("same");
  fourier_down[5][1][3][1]->Draw("same");
  fourier_fit[5][1][3][1]->Draw("same");
  */
  yield_rebin[5][4][3][1]->Draw("axis same");

  
  pi_over_2 = new TLine(TMath::Pi()/2.,0.1, TMath::Pi()/2.,47.);
  pi_over_2->SetLineStyle(2);
  pi_over_2->Draw();

  label_pp->Draw();

  c_demo->cd(5);

  yield_rebin[4][4][0][1]->SetMinimum(0.1);
  yield_rebin[4][4][0][1]->SetMaximum(47.);
 

  yield_rebin[4][4][0][1]->Draw();
  in_hist_rebin[4][4][0][1]->Rebin(2);
  in_hist_rebin[4][4][0][1]->Scale(1./2);
 in_hist_rebin[4][4][0][1]->Draw("same");

  pi_over_2->Draw();
  label_per->Draw();

  label_pT = new TLatex(0.2,0.8,"4<p_{T}^{assoc.}<8 GeV/c");
  label_pT->SetTextSize(0.07);
  label_pT->SetLineColor(kWhite);
  label_pT->SetNDC();
  label_pT->Draw();

  c_demo->cd(6);

  yield_rebin[4][4][3][1]->SetMinimum(0.1);
  yield_rebin[4][4][3][1]->SetMaximum(47.);

  yield_rebin[4][4][3][1]->Draw();
  in_hist_rebin[4][4][3][1]->Rebin(2);
  in_hist_rebin[4][4][3][1]->Scale(1./2);
  in_hist_rebin[4][4][3][1]->Draw("same");
 
  label_cent->Draw();
  pi_over_2->Draw();

  c_demo->SaveAs("SideBandDemonstrationPlot.png");
  c_demo->SaveAs("SideBandDemonstrationPlot.pdf");
 
  for(int j = 0; j<4; j++){
    for(int l = 0; l<3; l++){
      for(int i = 0; i<5; i++){

	//	if(i==0&&j==0&&l==0)   cout<<"CentIndex AjIndex TrkPtIndex V1 V2 V3"<<endl;
	cout<<"PbPb "<<j<<" "<<l<<" "<<i<<" "<<v1_tot[5][i][j][l]<<" "<<v2_tot[5][i][j][l]<<" "<<v3_tot[5][i][j][l]<<" "<<error[4][i][j][l]<<endl;

	//	cout<<"PbPb "<<j<<" "<<l<<" "<<i<<" "<<error[4][i][j][l]<<endl;
      }
    }
  }

  for(int j = 0; j<4; j++){
    for(int l = 0; l<3; l++){
      for(int i = 0; i<5; i++){

	cout<<"pp "<<j<<" "<<l<<" "<<i<<" "<<v1_tot[5][i][j][l]<<" "<<v2_tot[5][i][j][l]<<" "<<v3_tot[5][i][j][l]<<" "<<error[5][i][j][l]<<endl;

      }
    }
  }
  
  return 0;
}
