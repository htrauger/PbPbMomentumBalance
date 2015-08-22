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

#include <iostream>
#include <vector>
#include <fstream>

#include "../JetTrack2015_functions.h"


using namespace std;

Int_t spill_over_aj(bool is_number = kFALSE){

  //Set Style:
 
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.16);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  gStyle->SetTextFont(43);
  gStyle->SetCanvasBorderMode(0);
 
 

  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 6;
  const int nAjBins = 3;

 
enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};

TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};

  float CBins[nCBins+1] = {0, 10, 30, 50, 100};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%", "Cent. 30-50%","Cent. 50-100%"};

 
  float TrkPtBins[nTrkPtBins+1] = {0.5,1, 2, 3, 4, 8, 300};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt05","TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt300" };
  TString TrkPtBin_labels[nTrkPtBins] = {"0.5<pT<1","1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};



  float AjBins[nAjBins+1] = {0,0.22,0.75};
  TString AjBin_strs[nAjBins+1] = {"Aj0","Aj22","Aj75"};
 
  TF1 *fit0 = new TF1("fit0","[0]",-3.,3.);
 
  float x;
  TF1 *do_offset = new TF1("do_offset","-1.*[0]+x-x",-3.,3.);
  float offset;
   
  TFile *fin, *fin2, *fout;

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
 
  TH1D *closure_eta[12][6][4][3];
  TH1D *closure_phi[12][6][4][3];

  TH2D *closure_2D_direct[12][6][4][3];
  TH2D *closure_2D_direct_gen[12][6][4][3];
  TH1D *closure_eta_direct[12][6][4][3];
  TH1D *closure_phi_direct[12][6][4][3];
  TH1D *closure_eta_direct_rebin[12][6][4][3];
  TH1D *closure_phi_direct_rebin[12][6][4][3];

  TString in_name, plotname, outname, funcname, centlabel,Ajlabel, datalabel, jettype,jettype2, pTlabel,checkcanvasnameEta,checkcanvasnamePhi,rawetacanvasname,HminusPphicanvasname, HminusPetacanvasname,PbPbminuscanvas_phi,PbPbminuscanvas_eta, checkcanvasnameEta_wide, EtaClosureName;

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

  TH1D *Closure_integral_eta_cent[12][4][3];
  TH1D *Closure_integral_eta_cent2[12][4][3];

  TLine *lineCent, *linePt;


  TCanvas *cintegral_eta_pT[12][3];
  TCanvas *cintegral_phi_pT[12][3];
    
  TCanvas *cintegral_eta_cent[12][3];
  TCanvas *cintegral_phi_cent[12][3];



  Double_t xAxis[5] = {-100,-50,-30,-10,0}; 
  TH1D* int_cent[12][4][3];
  TH1D* blank[4][3];
  TH1D* blank2[4][3];

   

  Double_t check_ymax, check_ymin, dx_eta, dx_phi, bc, err, evalpt, temp1, err1;

  TF1 *gaus1d = new TF1("gaus1d","[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))");

  TF1 *gaus_phi[12][5][4][3];
  TF1 *gaus_eta[12][5][4][3];
 
  TLegend *lcheck, *leta, *lHminusP;
  
  TLine *linePhi, *lineEta;
 
  int mc_type_code;

  int llimiteta1, rlimiteta1,llimiteta2, rlimiteta2; 

  /////////////////////////

  etalim = 1.;
  philim = 1.;

  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  fin = new TFile("../me_correct_mc/HydJet_RecoJet_GenTrack_NoSube0_Dijet_Correlations.root","READ");
  fin2 = new TFile("../me_correct_mc/HydJet_GenJet_GenTrack_NoSube0_Dijet_Correlations.root","READ");
 
  if(is_number) fout = new TFile("Dijet_SpillOvers_NoPtWeight.root","RECREATE");
  else  fout = new TFile("Dijet_SpillOvers.root","RECREATE");

  for(int g=8; g<12; g++){

    //  There will only be one big "g-loop".
    if(g%2!=0)continue;
  
    //----------------------------------------------
    //     Start i and l loops -- set up canvases
    //------------------------------------------------

    for(int l = 0; l<3; l++){
     
      //----------------------------------------------------
      //  Start of main i & j loops 
      //-----------------------------------------------------

      for(int i=0; i<6; i++){
	
	for (int j=0; j<4; j++){

	  in_name = make_name("Yield_pTweighted_",g,i,j,l,centlabel, pTlabel, Ajlabel);
	  if(is_number) in_name = make_name("Yield_",g,i,j,l,centlabel, pTlabel, Ajlabel);


	  in_name.ReplaceAll("Pt100_Pt300_","");

	  cout<<in_name<<endl;

	  result[g][i][j][l] = (TH2D*)fin->Get(in_name)->Clone(in_name);
	
	  llimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-2.5+0.0001);
	  rlimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-1.-0.0001);

	  	    
	  llimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(1.+0.0001);
	  rlimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(2.5-0.0001);

	  background_left[g][i][j][l] = (TH1D*)result[g][i][j][l]->ProjectionY(Form("LeftSideBackground%d%d%d%d",g,i,j,l),llimiteta1,rlimiteta1);
	    
	  background_proj[g][i][j][l] = (TH1D*)result[g][i][j][l]->ProjectionY(Form("ProjectedBackground%d%d%d%d",g,i,j,l),llimiteta2,rlimiteta2);

	  background_proj[g][i][j][l]->Add(background_left[g][i][j][l]);

	  background_proj[g][i][j][l]->Scale(1./2/(rlimiteta1-llimiteta1+1));
	    
	  background[g][i][j][l] = (TH2D*)result[g][i][j][l]->Clone(Form("Background%d%d%d%d",g,i,j,l));

	  for(int k = 1;  k<result[g][i][j][l]->GetNbinsY(); k++){
	    temp1 = background_proj[g][i][j][l]->GetBinContent(k);
	    err1 = background_proj[g][i][j][l]->GetBinError(k);
	    
	    for(int m = 1;  m<result[g][i][j][l]->GetNbinsX(); m++){
	      background[g][i][j][l]->SetBinContent(m,k,temp1);
	      background[g][i][j][l]->SetBinError(m,k,err1);
		
	    }

	  }

	  result[g][i][j][l]->Add(background[g][i][j][l],-1.);

	  //-------------------------------
	  //dEta projection
	  //------------------------

	  TString eta_proj_name= in_name;
	  eta_proj_name.ReplaceAll("Yield","Eta_Proj");
	    
	  llimiteta = result[g][i][j][l]->GetXaxis()->FindBin(-etalim+.001);
	  rlimiteta = result[g][i][j][l]->GetXaxis()->FindBin(etalim-.001);

	  llimitphi = result[g][i][j][l]->GetYaxis()->FindBin(-philim+.001);
	  rlimitphi = result[g][i][j][l]->GetYaxis()->FindBin(philim-.001);
	    

	  eta_proj[g][i][j][l] = result[g][i][j][l]->ProjectionX(eta_proj_name,llimitphi,rlimitphi);
	  dx_eta = eta_proj[g][i][j][l]->GetBinWidth(1);
	  eta_proj[g][i][j][l]->Scale(1/dx_eta);

	  //	  eta_proj[g][i][j][l]->Scale(0.8);
	  

	  TString eta_proj_name_rebin = eta_proj_name;
	  eta_proj_name_rebin.ReplaceAll("Eta_Proj","Eta_Proj_Rebin");
	 
	  eta_proj_rebin[g][i][j][l] = (TH1D*)Rebin_dEta(eta_proj[g][i][j][l]);
	  eta_proj_rebin[g][i][j][l]->SetName(eta_proj_name_rebin);
	 
	  //-------------------------------
	  //dPhi projection
	  //------------------------

	  TString phi_proj_name= in_name;
	  phi_proj_name.ReplaceAll("Yield","Phi_Proj");

	  phi_proj[g][i][j][l] = result[g][i][j][l]->ProjectionY(phi_proj_name,llimiteta,rlimiteta);
	  dx_phi = phi_proj[g][i][j][l]->GetBinWidth(1);
	  phi_proj[g][i][j][l]->Scale(1/dx_phi);

	  TString phi_proj_name_rebin = phi_proj_name;
	  phi_proj_name_rebin.ReplaceAll("Phi_Proj","Phi_Proj_Rebin");
	 
	  phi_proj_rebin[g][i][j][l] = (TH1D*)Rebin_dPhi(phi_proj[g][i][j][l]);
	  phi_proj_rebin[g][i][j][l]->SetName(phi_proj_name_rebin);

	  float totbins = eta_proj_rebin[g][i][j][l]->GetNbinsX();

	
	    offset = (eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(1.001))+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(.99))+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(-.99))+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(-1.01))+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(1.51))+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(-1.51)))/6.;


  if(g==8&&i==1){
	  offset = (eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(1.001))+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(-1.01))+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(0.99))+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->FindBin(-0.99)))/4.;
	  }


	    do_offset->SetParameter(0, offset);

	  eta_proj_rebin[g][i][j][l]->Add(do_offset);
	  phi_proj_rebin[g][i][j][l]->Add(do_offset);
	  
	  for(int k = 1; k<eta_proj_rebin[g][i][j][l]->GetNbinsX()/2+1; k++){
	    bc = (eta_proj_rebin[g][i][j][l]->GetBinContent(k)+eta_proj_rebin[g][i][j][l]->GetBinContent(eta_proj_rebin[g][i][j][l]->GetNbinsX()+1-k))/2.;
	    eta_proj_rebin[g][i][j][l]->SetBinContent(k,bc);
	    eta_proj_rebin[g][i][j][l]->SetBinContent(eta_proj_rebin[g][i][j][l]->GetNbinsX()+1-k,bc);
	  }
	  

	  for(int k = 1; k<phi_proj_rebin[g][i][j][l]->GetNbinsX()/2+1; k++){
	    bc = (phi_proj_rebin[g][i][j][l]->GetBinContent(k)+phi_proj_rebin[g][i][j][l]->GetBinContent(phi_proj_rebin[g][i][j][l]->GetNbinsX()+1-k))/2.;
	    phi_proj_rebin[g][i][j][l]->SetBinContent(k,bc);
	    phi_proj_rebin[g][i][j][l]->SetBinContent(phi_proj_rebin[g][i][j][l]->GetNbinsX()+1-k,bc);
	  }
	  
	  if(g==2&&i==1&&j==0&&l==1){
	    eta_proj_rebin[g][i][j][l] = eta_proj_rebin[g][i][j][0];
	    phi_proj_rebin[g][i][j][l] = phi_proj_rebin[g][i][j][0];

	  }

	  //-----------------
	  //MAKE SPILLOVERS
	  //-----------------

	     
	  TString gaus_eta_name = in_name;
	  gaus_eta_name.ReplaceAll("Yield","GausFit_Eta");

	  EtaClosureName = in_name;
	  EtaClosureName.ReplaceAll("Pt100_Pt300_","");
	  EtaClosureName.ReplaceAll("Yield","SpillOvers_Eta");
	   	     
	  TString gaus_phi_name = gaus_eta_name;
	  gaus_phi_name.ReplaceAll("Eta","Phi");


	  TString PhiClosureName = EtaClosureName;
	  PhiClosureName.ReplaceAll("Eta","Phi");


	  check_ymin = -2.;
	  check_ymax = 11.;

	  

	  eta_proj_rebin[g][i][j][l]->SetLineColor(1);
	  eta_proj_rebin[g][i][j][l]->SetMarkerStyle(10);
	  eta_proj_rebin[g][i][j][l]->SetMarkerColor(1);
	  eta_proj_rebin[g][i][j][l]->SetMarkerSize(1);

	  eta_proj_rebin[g][i][j][l]->SetMaximum(check_ymax);
	  eta_proj_rebin[g][i][j][l]->SetMinimum(check_ymin);
	  eta_proj_rebin[g][i][j][l]->Draw();     


	  eta_proj_rebin[g][i][j][l]->GetXaxis()->CenterTitle();
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->CenterTitle();
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetTitle("#Delta#eta");
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitle("1/N_{jet} d^{2}N/(d#Delta#phi dp_{T})");
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);
	  if(i<3){
	    eta_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	  }
	  if(j>0){
	    eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	    eta_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	  }
	      
	  llimiteta = eta_proj_rebin[g][i][j][l]->GetXaxis()->FindBin(-0.5+.0001);
	  rlimiteta = eta_proj_rebin[g][i][j][l]->GetXaxis()->FindBin(0.5-.0001);
	     
	  double Yield_eta = eta_proj_rebin[g][i][j][l]->Integral(llimiteta,rlimiteta,"width");	      
	      
	     	     
	    
	  if(Yield_eta<0.){Yield_eta=0.;}

	      
	  gaus1d->FixParameter(0,0.);
	  gaus1d->FixParameter(1,Yield_eta);
	      
	  gaus1d->ReleaseParameter(2);
	  gaus1d->SetParameter(2,0.2);
	    
	  gaus1d->SetParLimits(2,0.2,1.5);

	   	    	      
	  eta_proj_rebin[g][i][j][l]->Fit("gaus1d","","",-1.5,1.5);
	      
	  
	  gaus_eta[g][i][j][l] = new TF1(gaus_eta_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

	  for(int k= 0; k<3; k++){
	
	    double temp = gaus1d->GetParameter(k);
	    gaus_eta[g][i][j][l]->SetParameter(k, temp);

	  }

	  double err_temp= gaus1d->GetParError(2);

	
	  eta_proj_rebin[g][i][j][l]->SetMaximum(check_ymax);
	  eta_proj_rebin[g][i][j][l]->SetMinimum(check_ymin);
	      
	  //    eta_proj_rebin[g][i][j][l]->SetMaximum(1.5);
	  eta_proj_rebin[g][i][j][l]->Draw();     

	 

	  closure_eta[g][i][j][l] = (TH1D*)eta_proj_rebin[g][i][j][l]->Clone(EtaClosureName);
   
	  for(int k=1; k< eta_proj_rebin[g][i][j][l]->GetNbinsX()+1; k++){
	    evalpt = eta_proj_rebin[g][i][j][l]->GetBinCenter(k);
	    bc = gaus1d->Eval(evalpt);
	    closure_eta[g][i][j][l]->SetBinContent(k,bc);
	    // closure_eta[g][i][j][l]->SetBinError(k,err_temp);
	    closure_eta[g][i][j][l]->SetBinError(k,0.);
	  }

	  gaus_eta[g][i][j][l]->SetLineColor(kBlue);
	  gaus_eta[g][i][j][l]->Draw("same");

	  closure_eta[g][i][j][l]->SetMarkerColor(kBlue);
	  closure_eta[g][i][j][l]->SetLineColor(kBlue);
	  closure_eta[g][i][j][l]->SetMarkerStyle(10);
	
	  
	  if(j>0){ lHminusP= new TLegend(legendoffset,texty3,texty1,0.95);
	    lHminusP->SetTextSize(ts);
	  }
	  if(j==0){ 
	      lHminusP = new TLegend(legendoffset2-.04,texty3,texty1,0.95);
	      lHminusP->SetTextSize(ts2);
	    }

	    lHminusP->SetFillColor(kWhite);
	    lHminusP->SetLineColor(kWhite);

	    lHminusP->AddEntry(eta_proj_rebin[g][i][j][l],"HYD.-PYTH.","lpfe");
	     	
	    if(i==3){lHminusP->SetTextSize(ts2);}
	    lHminusP->Draw("same");


	    phi_proj_rebin[g][i][j][l]->SetLineColor(1);
	    phi_proj_rebin[g][i][j][l]->SetMarkerStyle(10);
	    phi_proj_rebin[g][i][j][l]->SetMarkerColor(1);
	    phi_proj_rebin[g][i][j][l]->SetMarkerSize(1);


	    phi_proj_rebin[g][i][j][l]->SetMaximum(check_ymax);
	    phi_proj_rebin[g][i][j][l]->SetMinimum(check_ymin);
	    phi_proj_rebin[g][i][j][l]->Draw();     
      
 
	    phi_proj_rebin[g][i][j][l]->GetXaxis()->CenterTitle();
	    phi_proj_rebin[g][i][j][l]->GetYaxis()->CenterTitle();
	    phi_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	    phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	    phi_proj_rebin[g][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	    phi_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	    phi_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	    phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitle("1/N_{jet} d^{2}N/(d#Delta#phi dp_{T})");
	    phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	    phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);
	    if(i<3){
	      phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	    }
	    if(j>0){
	      phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	      phi_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    }

	    llimitphi = phi_proj_rebin[g][i][j][l]->GetXaxis()->FindBin(-1.0+.01);
	    rlimitphi = phi_proj_rebin[g][i][j][l]->GetXaxis()->FindBin(1.0-.001);
	   	      
	 
	    double Yield_phi = phi_proj_rebin[g][i][j][l]->Integral(llimitphi,rlimitphi,"width");	      
	      
	    Yield_phi = Yield_eta;
	 	      
	    //	    if(Yield_phi<0){Yield_phi=0;}
	    

	    gaus1d->SetParLimits(2,0.2,1.5);
	     	      
	    phi_proj_rebin[g][i][j][l]->Fit("gaus1d");
	    
	    gaus_phi[g][i][j][l] = new TF1(gaus_phi_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

	    for(int k= 0; k<3; k++){
	      double temp = gaus1d->GetParameter(k);
	      gaus_phi[g][i][j][l]->SetParameter(k, temp);
	    }

	    err_temp = gaus1d->GetParError(2);
	 
	    phi_proj_rebin[g][i][j][l]->SetMaximum(check_ymax);
	    phi_proj_rebin[g][i][j][l]->SetMinimum(check_ymin);
	    phi_proj_rebin[g][i][j][l]->Draw();    

	    gaus_phi[g][i][j][l]->SetLineColor(kBlue);
	    gaus_phi[g][i][j][l]->Draw("same"); 


	    closure_phi[g][i][j][l] = (TH1D*)phi_proj_rebin[g][i][j][l]->Clone(PhiClosureName);
      
	    for(int k=0; k< phi_proj_rebin[g][i][j][l]->GetNbinsX()+1; k++){
	      evalpt = phi_proj_rebin[g][i][j][l]->GetBinCenter(k);
	      bc = gaus1d->Eval(evalpt);
	      closure_phi[g][i][j][l]->SetBinContent(k,bc);
	      //  closure_phi[g][i][j][l]->SetBinError(k,err_temp);
	      closure_phi[g][i][j][l]->SetBinError(k,0.);
	    }

	    gaus_phi[g][i][j][l]->SetLineColor(kBlue);

	    closure_phi[g][i][j][l]->SetMarkerColor(kBlue);
	    closure_phi[g][i][j][l]->SetLineColor(kBlue);
	    closure_phi[g][i][j][l]->SetMarkerStyle(10);
	    closure_phi[g][i][j][l]->Draw("same");

	    cout<<"Here"<<endl;


 
	    eta_proj_rebin[g][i][j][l]->Write();
	    closure_eta[g][i][j][l]->Write();
	    gaus_eta[g][i][j][l]->Write();

	    phi_proj_rebin[g][i][j][l]->Write();
	    closure_phi[g][i][j][l]->Write();
	    gaus_phi[g][i][j][l]->Write();
	  

	  	
	} //close j
      } // close i
     
  
   
    }  // Closes the [l] (gen vs. reco) loop

 
  }// closes g loop



  cout<<"starting plots for AN"<<endl;
  //**********************
  //  CLOSURE PLOTS FOR AN
  //***********************

  TCanvas *cClosuresAN_eta[12][4][3],*cClosuresAN_phi[12][4][3];
  float result_min, result_max,  val_l,val_r,err_l,err_r, value, pt_val;
  TLegend *l40,*l41,*l42;

 for(int g=8;g<12;g++){
      
    switch(g){
    case 8: 
      datalabel = "Subleading"; break;
    case 10:
      datalabel = "Leading"; break;
    default:
      continue; 
      break;
    }

    for(int l = 0; l<3; l++){

    for(int i = 0; i<6; i++){
      TString pTrange;
      switch(i){
      case 0: 
	result_max = 11.;
	result_min = -1.;
	break;
      case 1: 
	result_max = 11.;
	result_min = -1.;
	break;
      case 2: 
	result_max = 8.5;
	result_min = -1.; 
	break;
      case 3: 
	result_max = 6.2;
	result_min = -.45;
	break;
      case 4: 
	result_max = 11.;
	result_min = -.1;
	break;
      case 5: 
	result_max = 11.;
	result_min = -1.;
	break;
      }
   
      in_name = make_name("Yield_",g,i,0,l,centlabel, pTlabel, Ajlabel);

      TString ClosureANeta_name = in_name;
      ClosureANeta_name.ReplaceAll("Yield_","SpillOver_Eta_");
      ClosureANeta_name.ReplaceAll("Cent50_Cent100_Pt100_Pt300_","");
      cClosuresAN_eta[g][i][l] = new TCanvas(ClosureANeta_name," ",10,10,1500,400);
      cClosuresAN_eta[g][i][l]->Divide(4,1,0.,0.);
   


      TString ClosureANphi_name = in_name;
      ClosureANphi_name.ReplaceAll("Yield_","SpillOver_Phi_");
      ClosureANphi_name.ReplaceAll("Cent50_Cent100_Pt100_Pt300_","");
      cClosuresAN_phi[g][i][l] = new TCanvas(ClosureANphi_name," ",10,10,1500,400);      cClosuresAN_phi[g][i][l]->Divide(4,1,0.,0.);
	

      for(int j=0; j<4; j++){

	in_name = make_name("Yield_",g,i,j,l,centlabel, pTlabel, Ajlabel);

	cClosuresAN_eta[g][i][l]->cd(j+1);
	eta_proj_rebin[g][i][j][l]->SetMinimum(result_min);
	eta_proj_rebin[g][i][j][l]->SetMaximum(result_max);
	eta_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(ts2);
	eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle2);
	eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleOffset(1.);
	eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitle("");
	eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(0);
	
	if(j==0){ eta_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(ts2-0.01);
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleSize(ts);
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(ts2);
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelOffset(0.013);
	}
	eta_proj_rebin[g][i][j][l]->SetAxisRange(-1.4,1.41,"x");
	eta_proj_rebin[g][i][j][l]->Draw();
	gaus_eta[g][i][j][l]->Draw("same");

	//	closure_eta_ref[g][i][j][l]->Draw("same");
	eta_proj_rebin[g][i][j][l]->Draw("same");

	drawlabels(g,i,j);

	if(j==3){
	  TLatex *tex24eta = new TLatex(textalign,texty2,phirangelabel);
	  tex24eta->SetName("tex24eta");
	  tex24eta->SetNDC();
	  tex24eta->SetTextSizePixels(tspixels);
	  tex24eta->Draw();
	}

	lineEta = new TLine(-1.5,0,1.5,0); 
	lineEta->Draw("same");
	
	


	lineEta = new TLine(-1.5,0,1.5,0);
	lineEta->SetLineStyle(2);
	lineEta->Draw("same");


	lineEta->Draw("same");

	 	 
	//----------------------------------
	// Closure Plots for AN (dPhi)
	//----------------------------------

	cClosuresAN_phi[g][i][l]->cd(j+1);
	  
	phi_proj_rebin[g][i][j][l]->SetMinimum(-1.);
	phi_proj_rebin[g][i][j][l]->SetMaximum(result_max);
	phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(ts2);
	phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle2);
	phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleOffset(1.);
	phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitle("");
	phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(0.);
	  

	if(j==0){ phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(ts2-0.01);
	  phi_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleSize(ts);
	  phi_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(ts2);
	  phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelOffset(0.013);
	}

	phi_proj_rebin[g][i][j][l]->SetAxisRange(-1.4,1.41,"x");
	phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(ts2);
	if(j==0){phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(ts2-0.01);}
	phi_proj_rebin[g][i][j][l]->Draw();
	gaus_phi[g][i][j][l]->Draw("same");

	//	closure_phi_ref[g][i][j][l]->Draw("same");
	phi_proj_rebin[g][i][j][l]->Draw("same");
	 
	drawlabels(g,i,j);

	if(j==3){
	  TLatex *tex24phi = new TLatex(textalign,texty2,etarangelabel);
	  tex24phi->SetName("tex24phi");
	  tex24phi->SetNDC();
	  tex24phi->SetTextSizePixels(tspixels);
	  tex24phi->Draw();
	}



	linePhi = new TLine(-1.51,0,1.51,0);
	  linePhi->SetLineStyle(2);
	  linePhi->Draw("same");


	linePhi->Draw("same");
	    
	    if(i==0&&j==0){
	      Closure_integral_eta0.clear();
	      Closure_integral_phi0.clear();
	      Closure_integral_eta1.clear();
	      Closure_integral_phi1.clear();
	      Closure_integral_eta2.clear();
	      Closure_integral_phi2.clear();
	      Closure_integral_eta3.clear();
	      Closure_integral_phi3.clear();
	    }

	    llimiteta = eta_proj_rebin[g][i][j][l]->GetXaxis()->FindBin(-1.0+.0001);
	    rlimiteta = eta_proj_rebin[g][i][j][l]->GetXaxis()->FindBin(1.0-.0001);
	     
	    double Yield_eta = eta_proj_rebin[g][i][j][l]->Integral(llimiteta,rlimiteta,"width");	      

	    if(Yield_eta<0.){Yield_eta=0.;}

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
	 
      }//close j


      cClosuresAN_eta[g][i][l]->cd(0);
 
      TLatex *canvas_title = new TLatex(0.06,0.92,"CMS Preliminary Simulation");
      canvas_title->SetTextSizePixels(tspixels);
      canvas_title->SetTextFont(63);
      canvas_title->Draw();

      TLatex *canvas_title2 = new TLatex(0.295,0.92,"PYTHIA+HYDJET");
      canvas_title2->SetTextSizePixels(tspixels);
      canvas_title2->Draw();

      ClosureANeta_name+=".pdf";
      cClosuresAN_eta[g][i][l]->SaveAs(ClosureANeta_name);
      ClosureANeta_name.ReplaceAll("pdf","png");
      cClosuresAN_eta[g][i][l]->SaveAs(ClosureANeta_name);


      cClosuresAN_phi[g][i][l]->cd(0);
 
      canvas_title->Draw();
      canvas_title2->Draw();

      ClosureANphi_name+=".pdf";
      cClosuresAN_phi[g][i][l]->SaveAs(ClosureANphi_name);
      ClosureANphi_name.ReplaceAll("pdf","png");
      cClosuresAN_phi[g][i][l]->SaveAs(ClosureANphi_name);
    

    }//close i 
 
        
      TString integral_eta_pT_name = "integral_eta_pT";
      integral_eta_pT_name+=g;    integral_eta_pT_name+=l;

      cintegral_eta_pT[g][l] = new TCanvas(integral_eta_pT_name,"",10,10,1500,500);
      cintegral_eta_pT[g][l]->Divide(4,1,0.,0.);



      for(int j = 0; j<4; j++){


	in_name = make_name("Result_",g,3,j,0,centlabel,pTlabel,Ajlabel);


	cintegral_eta_pT[g][l]->cd(j+1);

	TString ClosureIntegralEtaPt_name = in_name;
	ClosureIntegralEtaPt_name.ReplaceAll("Result","Closure_Integral_Eta");
	ClosureIntegralEtaPt_name.ReplaceAll("_TrkPt8_TrkPt300","");
	ClosureIntegralEtaPt_name.ReplaceAll("_Pt100_Pt300","");

	cout<< ClosureIntegralEtaPt_name<<endl;	

	cout<<"here"<<endl;

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

	cout<<"and here"<<endl;
	Closure_integral_eta_pT[g][j][l]->SetName(ClosureIntegralEtaPt_name);

	Closure_integral_eta_pT[g][j][l]->SetMarkerColor(1);
	Closure_integral_eta_pT[g][j][l]->SetLineColor(1);

	switch(g){
	case 8:
	  Closure_integral_eta_pT[g][j][l]->SetMarkerStyle(34);
	  break;
	case 10: 
	  Closure_integral_eta_pT[g][j][l]->SetMarkerStyle(21);
	  break;
	default:
	  Closure_integral_eta_pT[g][j][l]->SetMarkerStyle(10);
	  break;
	}

	Closure_integral_eta_pT[g][j][l]->SetMinimum(-1.);
	Closure_integral_eta_pT[g][j][l]->SetMaximum(3.9);
	    
	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetRangeUser(.501,7.9);
	Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetNdivisions(306);
	Closure_integral_eta_pT[g][j][l]->Draw("p X A");
	 

	Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetLabelSize(ts);
	   


	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");
	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetTitleSize(ts2);
	Closure_integral_eta_pT[g][j][l]->GetXaxis()->SetTitleOffset(xoffset+0.2);
	Closure_integral_eta_pT[g][j][l]->GetYaxis()->SetTitle("(dN/dp_{T})_{P+H} - (dN/dp_{T})_{PYTH} (GeV/c)^{-1}");

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



	linePt = new TLine(1.,0,8.,0);
	linePt->SetLineStyle(2);
	linePt->SetLineWidth(1);
	linePt->Draw("same");

	cout<<"here and g = "<<g<<endl;
   
	if(g!=10)continue;

	  closure_integral_values.clear();
	  closure_integral_errors.clear();
	
	  for(int k = 0; k<6; k++){
	    double pt_val, x_val;
	  
	    Closure_integral_eta_pT[g][j][l]->GetPoint(k,x_val,pt_val);
	    closure_integral_values.push_back(pt_val);
	    closure_integral_errors.push_back(pt_val/2.);
	  
	  }
	 	  
	  Closure_integral_eta_pT2[g][j][l] = new TGraphErrors(pTbin_centers.size(),&pTbin_centers[0],&closure_integral_values[0],&pTbin_errors[0],&closure_integral_errors[0]);

	  cout<<"and even here"<<endl;

	  Closure_integral_eta_pT2[10][j][l]->SetFillColor(kOrange-2);
	  Closure_integral_eta_pT2[10][j][l]->Draw("same e2");
	  
	  
	  Closure_integral_eta_pT[10][j][l]->Draw("same p X");

	  Closure_integral_eta_pT[8][j][l]->SetMarkerColor(kCyan+3);
	  Closure_integral_eta_pT[8][j][l]->SetLineColor(kCyan+3);
	  Closure_integral_eta_pT[8][j][l]->Draw("same p X");
	
	  if(j==0){ 
	    l40 = new TLegend(textalign2,texty1-.05,0.8,texty4-.1);
	    l40->SetName("l40");
	    l40->SetTextFont(43);
	    l40->SetTextSizePixels(tspixels);
	    l40->SetFillColor(kWhite);
	    l40->SetLineColor(kWhite);
	    if(g==10){
	      l40->AddEntry(Closure_integral_eta_pT[10][j][l],"Leading Spill-Over","p");
	      l40->AddEntry(Closure_integral_eta_pT[8][j][l],"Subleading Spill-Over","p");
	      l40->AddEntry(Closure_integral_eta_pT2[10][j][l],"Uncertainty Assigned for","f");
	      l40->AddEntry(Closure_integral_eta_pT2[10][j][l],"Leading Jet Spill-Over","");
	    }
	   
	    l40->Draw("same");

	  
	  }
	  drawlabels_int_pt2(g,j);
	 
	  if(j==3){
	    TLatex *aj_tex = new TLatex(0.05, 0.8,Ajlabel);
	    aj_tex->SetTextSizePixels(tspixels);
	    aj_tex->SetNDC();
	    aj_tex->Draw();
	  }

      }
				       
      cintegral_eta_pT[g][l]->cd(0);
								      
      TLatex *canvas_title = new TLatex(0.06,0.9,"CMS Preliminary Simulation");
      canvas_title->SetTextSizePixels(tspixels);
      canvas_title->SetTextFont(63);
      canvas_title->Draw();

      TLatex *canvas_title2 = new TLatex(0.295,0.9,"PYTHIA+HYDJET");
      canvas_title2->SetTextSizePixels(tspixels);
      canvas_title2->Draw();


      if(l==2){

	cintegral_eta_pT[g][l]->SaveAs("Integral_SpillOver_pT_Leading_AjInclusive.pdf");
	cintegral_eta_pT[g][l]->SaveAs("Integral_SpillOver_pT_Leading_AjInclusive.png");
      }else{
	cintegral_eta_pT[g][l]->SaveAs("Integral_SpillOver_pT_Leading_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".pdf");
	cintegral_eta_pT[g][l]->SaveAs("Integral_SpillOver_pT_Leading_"+AjBin_strs[l]+"_"+AjBin_strs[l+1]+".png");
       }

    }//close Aj for plots


 } //close g for plots
 

 
 return 0;


}  //and we're done.
 
