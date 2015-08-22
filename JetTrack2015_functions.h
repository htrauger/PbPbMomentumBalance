
#include "JetTrack2015_universals.h"


TString make_name(TString stem, int g, int i, int j, int l, TString &centlabel, TString &pTlabel, TString &Ajlabel){

  TString in_name = stem;

  switch(g){
  case 0:
    in_name+= "PbPb_";
    break;
  case 1:
    in_name+= "pp_";
    break;
  case 2: 
    in_name+= "SubLeading_PbPb_";
    break;
  case 3:
    in_name+= "SubLeading_pp_";
    break;
  case 4: 
    in_name+= "Leading_PbPb_";
    break;
  case 5:
    in_name+= "Leading_pp_";
    break;
  case 6:
    in_name+= "Hydjet_";
    break;
  case 7:
    in_name+= "Pythia_";
    break;
  case 8: 
    in_name+= "SubLeading_Hydjet_";
    break;
  case 9:
    in_name+= "SubLeading_Pythia_";
    break;
  case 10: 
    in_name+= "Leading_Hydjet_";
    break;
  case 11:
    in_name+= "Leading_Pythia_";
    break;
  }
  

  switch(j){
  case 0: 
   in_name+= "Cent50_Cent100_Pt100_Pt300_";
    centlabel = "Cent. 50-100%";
    break;

  case 1:
    in_name+= "Cent30_Cent50_Pt100_Pt300_";
    centlabel = "Cent. 30-50%";
    break;

  case 2: 
    in_name+= "Cent10_Cent30_Pt100_Pt300_";
    centlabel = "Cent. 10-30%";
    break;

  case 3:
    in_name+= "Cent0_Cent10_Pt100_Pt300_";
    centlabel = "Cent. 0-10%";
    break;

  default:
    break;
  }

  switch(l){
  case 0:
    Ajlabel = "A_{J}<0.22";
    in_name+= "Aj0_Aj22_";
    break;
 case 1:
    Ajlabel = "A_{J}>0.22";
    in_name+= "Aj22_Aj75_";
    break;
  case 2:
    Ajlabel = "A_{J} Inclusive";
    in_name+= "";
    break;
  }

  switch(i){
 case 0:
    pTlabel = "0.5<p_{T}<1";
    in_name+= "TrkPt05_TrkPt1";
    break;
  case 1:
    pTlabel = "1<p_{T}<2";
    in_name+= "TrkPt1_TrkPt2";
    break;
  case 2:
    pTlabel = "2<p_{T}<3";
    in_name+= "TrkPt2_TrkPt3";
    break;
  case 3:
    pTlabel = "3<p_{T}<4";
    in_name+= "TrkPt3_TrkPt4";
    break;
  case 4:
    pTlabel = "4<p_{T}<8";
    in_name+= "TrkPt4_TrkPt8";
    break;
  case 5:
    pTlabel = "p_{T}>8";
    in_name+= "TrkPt8_TrkPt300";
    break;
  }

  return in_name;
}






TH1D *Rebin_dPhi(TH1D* hold){
  // New bin boundaries
  // Double_t bin_bounds_phi[]   =  {-1.50796,-1.256637, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531, 1.256637,1.50796};
 
  Int_t new_bins = sizeof(bin_bounds_phi) / sizeof(bin_bounds_phi[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_phi[i] - hold->GetBinLowEdge(j) )<=1e-3 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_phi[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }
    
  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_phi );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}




TH1D *Rebin_dEta(TH1D* hold){
  // New bin boundaries
  Double_t bin_bounds_eta[]
    //  =   {-2.5,-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.,2.5};
    =   {-3.5,-3.,-2.5,-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.,2.5,3.,3.5};
 
  Int_t new_bins = sizeof(bin_bounds_eta) / sizeof(bin_bounds_eta[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_eta[i] - hold->GetBinLowEdge(j) )<1e-4 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_eta[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }

  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_eta );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}




TH1D *Rebin_dEta2(TH1D* hold){
  // New bin boundaries
  Double_t bin_bounds_eta[]
    //  =   {-2.5,-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.,2.5};
    =   {-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.};
 
  Int_t new_bins = sizeof(bin_bounds_eta) / sizeof(bin_bounds_eta[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_eta[i] - hold->GetBinLowEdge(j) )<1e-4 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_eta[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }

  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_eta );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}






TH1D *Rebin_dEta3(TH1D* hold){
  // New bin boundaries
  Double_t bin_bounds_eta[]
    =   {-3.,-2.5,-2.,-1.5,  -1,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1, 1.5,2.,2.5,3.};
 
  Int_t new_bins = sizeof(bin_bounds_eta) / sizeof(bin_bounds_eta[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_eta[i] - hold->GetBinLowEdge(j) )<1e-4 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_eta[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }

  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_eta );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}


 



TH1D *Rebin_dPhi_full(TH1D* hold){
  // New bin boundaries
  Double_t bin_bounds_phi[]=  {-1.507968,-1.382304,-1.25664,-1.130976,-1.005312,-0.879648,-0.753984,-0.62832,-0.502656,-0.376992,-0.251328,-0.125664,0.125664,0.251328,0.376992,0.502656,0.62832,0.753984,0.879648,1.005312,1.130976,1.25664,1.382304,1.507968,1.633632,1.82212,2.01062,2.19911,2.38761,2.57611,2.7646,2.9531,3.14159,3.33009,3.51858,3.70708,3.89557,4.08407,4.27257,4.46106,4.58673};

  Int_t new_bins = sizeof(bin_bounds_phi) / sizeof(bin_bounds_phi[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_phi[i] - hold->GetBinLowEdge(j) )<=1e-3 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_phi[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }
    
  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_phi );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}


TH1D *Rebin_dPhi_full_old(TH1D* hold){
  // New bin boundaries
  
 Double_t bin_bounds_phi[]=  {-1.507968,-1.382304,-1.25664,-1.130976,-1.005312,-0.879648,-0.753984,-0.62832,-0.502656,-0.376992,-0.251328,-0.125664,0.125664,0.251328,0.376992,0.502656,0.62832,0.753984,0.879648,1.005312,1.130976,1.25664,1.382304,1.507968,1.633632,1.759296,1.88496,2.010624,2.136288,2.261952,2.387616,2.51328,2.638944,2.764608,2.890272,3.015936,3.1416,3.267264,3.392928,3.518592,3.644256,3.76992,3.895584,4.021248,4.146912,4.272576,4.39824,4.523904,4.649568,4.775232,4.900896,5.02656};


  Int_t new_bins = sizeof(bin_bounds_phi) / sizeof(bin_bounds_phi[0])-1; // Should be 35 or so now
  

  // important cross check. These need to be binboundaries of the old one!
  for (int i=0; i<new_bins+1 ; ++i ){
    bool binok=false;
    for (int j=1; j<=hold->GetNbinsX()+1 ; ++j ){
      if ( fabs( bin_bounds_phi[i] - hold->GetBinLowEdge(j) )<=1e-3 ) {
	binok=true;
	break;
      }
    }
    if (!binok ){
      cerr << "Bin edge " << bin_bounds_phi[i] << " is NOT a bin boundary" << endl;
      return 0;
    }
  }
    
  TString new_name = TString( hold->GetName() ) + "_Rebin" ;
  TH1D* hnew = (TH1D*) hold->Rebin ( new_bins, new_name, bin_bounds_phi );
  
  // Now renormalize per dPhi
  Double_t dx = hold->GetXaxis()->GetBinWidth(1);
  for (int j=1; j<=hnew->GetNbinsX()+1 ; ++j ){
    hnew->SetBinContent( j, hnew->GetBinContent( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
    hnew->SetBinError(   j, hnew->GetBinError( j ) * dx / hnew->GetXaxis()->GetBinWidth(j) );
  }
      

  return hnew;
}




void drawlabels(int g, int i, int j){
  
  TLatex *tex00, *tex01, *tex02,*tex03,*tex10, *tex11, *tex12, *tex13, *tex20, *tex21, *tex22, *tex23, *tex30, *tex31, *tex32, *tex33;

  gStyle->SetTextFont(43);
  
  float texty0=1.05;

  switch(j){
  case 0:{
	  
    tex13 = new TLatex(textalign2,texty0,"CMS Preliminary");
    tex13->SetName("tex11");
    tex13->SetNDC();
    tex13->SetTextSizePixels(tspixels);
    tex13->SetTextFont(63);
    //tex13->Draw();   	  
 
    tex00 = new TLatex(textalign2,texty1,"Centrality 50-100%");
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->Draw("same");
	    
    switch(i){
    case 0:
      tex01 = new TLatex(textalign2,texty2,"0.5<p_{T}^{assoc}<1GeV/c");
      break;  
    case 1:
      tex01 = new TLatex(textalign2,texty2,"1<p_{T}^{assoc}<2GeV/c");
      break;
    case 2:
      tex01 = new TLatex(textalign2,texty2,"2<p_{T}^{assoc}<3GeV/c");
      break;
    case 3:
      tex01 = new TLatex(textalign2,texty2,"3<p_{T}^{assoc}<4GeV/c");
      break;
    case 4:
      tex01 = new TLatex(textalign2,texty2,"4<p_{T}^{assoc}<8GeV/c");
      break;
    case 5:
      tex01 = new TLatex(textalign2,texty2,"p_{T}^{assoc}>8GeV/c");
      break;
    }
	      
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->SetLineWidth(2);
    tex00->Draw("same");
    tex01->SetName("tex01");
    tex01->SetNDC();
    tex01->SetTextSizePixels(tspixels);
    tex01->SetLineWidth(2);
    tex01->Draw("same");
  }
    break;
  case 1:   
    {
      tex10 = new TLatex(textalign,texty1,"Centrality 30-50%");
      tex10->SetName("tex10");
      tex10->SetNDC();
      tex10->SetTextSizePixels(tspixels);
      tex10->Draw("same");
    
      if(g<6){
	tex11 = new TLatex(textalign,texty0,"PbPb 166 #mub^{-1} (2.76 TeV)");
	tex11->SetName("tex11");
	tex11->SetNDC();
	tex11->SetTextSizePixels(tspixels);
	//	tex11->SetTextFont(63);
	//	tex11->Draw();
     
      }

     if(g==0||g==1||g==6||g==7){
	tex31 = new TLatex(textalign,texty2,"p_{T,jet}>120 GeV/c");
	tex31->SetNDC();
	tex31->SetTextSizePixels(tspixels);
	tex31->SetLineWidth(2);
	tex31->Draw();
     }
     if((g>1&&g<6)||g>7){
       tex31 = new TLatex(textalign,texty2,"p_{T, jet1}> 120 GeV/c");
       tex31->SetName("tex31");
       tex31->SetNDC();
       tex31->SetTextSizePixels(tspixels);
       tex31->SetLineWidth(2);
       tex31->Draw();
       tex33 = new TLatex(textalign,texty3,"p_{T, jet2}> 50 GeV/c");
       tex33->SetName("tex33");
       tex33->SetNDC();
       tex33->SetTextSizePixels(tspixels);
       tex33->SetLineWidth(2);
       tex33->Draw();
       tex32 = new TLatex(textalign,texty4,"#Delta#phi_{1,2}> 5#pi/6");
       tex32->SetName("tex32");
       tex32->SetNDC();
       tex32->SetTextSizePixels(tspixels);
       tex32->SetLineWidth(2);
       tex32->Draw();
     }




    }
    break;

  case 2:
    {
      tex20 = new TLatex(textalign,texty1,"Centrality 10-30%");
      tex20->SetName("tex02");
      tex20->SetNDC();
      tex20->SetTextSizePixels(tspixels);
      tex20->SetLineWidth(2);
      tex20->Draw();
	      	   
   
      if(g<6){

     	tex12 = new TLatex(textalign,texty0,"pp 5.3 pb^{-1} (2.76 TeV)");
	tex12->SetName("tex12");
	tex12->SetNDC();
	tex12->SetTextSizePixels(tspixels);
	//	tex12->SetTextFont(63);
	//	tex12->Draw();
	
      }




   tex21 = new TLatex(textalign,texty2,"anti-k_{T} jets, R=0.3");
      tex21->SetName("tex2");
      tex21->SetNDC();
      tex21->SetTextSizePixels(tspixels);
	   	    
      tex22 = new TLatex(textalign,texty3,"|#eta_{jet}|<1.6");
      tex22->SetName("tex22");
      tex22->SetNDC();
      tex22->SetTextSizePixels(tspixels);
      tex21->Draw();
      tex22->Draw();

    }
    break; 
	 
  case 3:
    {
      tex30 = new TLatex(textalign,texty1,"Centrality 0-10%");
      tex30->SetName("tex30");
      tex30->SetNDC();
      tex30->SetTextSizePixels(tspixels);
      tex30->SetLineWidth(2);
      tex30->Draw("same");

    }
    break;
  }

}




void drawlabels_closures(int g, int i, int j){
  
  TLatex *tex00, *tex01, *tex02,*tex03,*tex10, *tex11, *tex12, *tex13, *tex20, *tex21, *tex22, *tex23, *tex30, *tex31, *tex32, *tex33;

  gStyle->SetTextFont(43);
  
  float texty0=1.05;

  switch(j){
  case 0:{
	  
    tex13 = new TLatex(textalign,texty0,"CMS Preliminary");
    tex13->SetName("tex11");
    tex13->SetNDC();
    tex13->SetTextSizePixels(tspixels);
    tex13->SetTextFont(63);
    //tex13->Draw();   	  
 
    tex00 = new TLatex(textalign,texty1,"Centrality 50-100%");
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->Draw("same");
	    
    switch(i){
    case 0:
      tex01 = new TLatex(textalign,texty2,"0.5<p_{T}^{assoc}<1GeV/c");
      break;
    case 1:
      tex01 = new TLatex(textalign,texty2,"1<p_{T}^{assoc}<2GeV/c");
      break;
    case 2:
      tex01 = new TLatex(textalign,texty2,"2<p_{T}^{assoc}<3GeV/c");
      break;
    case 3:
      tex01 = new TLatex(textalign,texty2,"3<p_{T}^{assoc}<4GeV/c");
      break;
    case 4:
      tex01 = new TLatex(textalign,texty2,"4<p_{T}^{assoc}<8GeV/c");
      break;
    case 5:
      tex01 = new TLatex(textalign,texty2,"p_{T}^{assoc}>8GeV/c");
      break;
    }
	      
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->SetLineWidth(2);
    tex00->Draw("same");
    tex01->SetName("tex01");
    tex01->SetNDC();
    tex01->SetTextSizePixels(tspixels);
    tex01->SetLineWidth(2);
    tex01->Draw("same");
  }
    break;
  case 1:   
    {
      tex10 = new TLatex(textalign,texty1,"Centrality 30-50%");
      tex10->SetName("tex10");
      tex10->SetNDC();
      tex10->SetTextSizePixels(tspixels);
      tex10->Draw("same");
    
      if(g<6){
	tex11 = new TLatex(textalign,texty0,"PbPb 166 #mub^{-1} (2.76 TeV)");
	tex11->SetName("tex11");
	tex11->SetNDC();
	tex11->SetTextSizePixels(tspixels);
	//	tex11->SetTextFont(63);
	//	tex11->Draw();
      }


	if(g==0||g==1||g==6||g==7){
	  tex31 = new TLatex(textalign,texty2,"p_{T,jet}>120 GeV/c");
	  tex31->SetNDC();
	  tex31->SetTextSizePixels(tspixels);
	  tex31->SetLineWidth(2);
	  tex31->Draw();
	}
	if((g>1&&g<6)||g>7){
	  tex31 = new TLatex(textalign,texty2,"p_{T, jet1}> 120 GeV/c");
	  tex31->SetName("tex31");
	  tex31->SetNDC();
	  tex31->SetTextSizePixels(tspixels);
	  tex31->SetLineWidth(2);
	  tex31->Draw();
	  tex33 = new TLatex(textalign,texty3,"p_{T, jet2}> 50 GeV/c");
	  tex33->SetName("tex33");
	  tex33->SetNDC();
	  tex33->SetTextSizePixels(tspixels);
	  tex33->SetLineWidth(2);
	  tex33->Draw();
	  tex32 = new TLatex(textalign,texty4,"#Delta#phi_{1,2}> 5#pi/6");
	  tex32->SetName("tex32");
	  tex32->SetNDC();
	  tex32->SetTextSizePixels(tspixels);
	  tex32->SetLineWidth(2);
	  tex32->Draw();
	}
    
    }
    break;

  case 2:
    {
      tex20 = new TLatex(textalign,texty1,"Centrality 0-30%");
      tex20->SetName("tex02");
      tex20->SetNDC();
      tex20->SetTextSizePixels(tspixels);
      tex20->SetLineWidth(2);
      tex20->Draw();
	      	   
   
      if(g<6){

     	tex12 = new TLatex(textalign,texty0,"pp 5.3 pb^{-1} (2.76 TeV)");
	tex12->SetName("tex12");
	tex12->SetNDC();
	tex12->SetTextSizePixels(tspixels);
	//	tex12->SetTextFont(63);
	//	tex12->Draw();
	
      }




   tex21 = new TLatex(textalign,texty2,"anti-k_{T} jets, R=0.3");
      tex21->SetName("tex2");
      tex21->SetNDC();
      tex21->SetTextSizePixels(tspixels);
	   	    
      tex22 = new TLatex(textalign,texty3,"|#eta_{jet}|<1.6");
      tex22->SetName("tex22");
      tex22->SetNDC();
      tex22->SetTextSizePixels(tspixels);
      tex21->Draw();
      tex22->Draw();

    }
    break; 
	 
  case 3:
    {
      tex30 = new TLatex(textalign,texty1,"Centrality 0-10%");
      tex30->SetName("tex30");
      tex30->SetNDC();
      tex30->SetTextSizePixels(tspixels);
      tex30->SetLineWidth(2);
      tex30->Draw("same");


    }
    break;
  }

}



void drawlabels_int_cent2(int g, int i){
  
  TLatex *tex00, *tex01, *tex02,*tex03,*tex10, *tex11, *tex12, *tex13, *tex20, *tex21, *tex22, *tex23, *tex30, *tex31, *tex32, *tex33;
  

  textalign = 0.18;
  

  gStyle->SetTextFont(43);
  
  float texty0=1.;


  texty0+=-0.11;
  texty1+=-0.11;
  texty2+=-0.11;
  texty3+=-0.11;
  texty4+=-0.11;
  


  switch(i){
  case 0:{
	  
    tex13 = new TLatex(textalign,texty0,"CMS Preliminary");
    tex13->SetName("tex13");
    tex13->SetNDC();
    tex13->SetTextSizePixels(tspixels);
    tex13->SetTextFont(63);
    // tex13->Draw();   	  
 
    tex00 = new TLatex(textalign,texty1,"1<p_{T}^{assoc}<2GeV/c");
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->Draw("same");

  }
    break;
  case 1:   
    {
      tex10 = new TLatex(textalign,texty1,"2<p_{T}^{assoc}<3 GeV/c");
      tex10->SetName("tex10");
      tex10->SetNDC();
      tex10->SetTextSizePixels(tspixels);
      tex10->Draw("same");
    
      if(g<6){
	tex11 = new TLatex(textalign,texty0,"PbPb 166 #mub^{-1} (2.76 TeV)");
	tex11->SetName("tex11");
	tex11->SetNDC();
	tex11->SetTextSizePixels(tspixels);
	//	tex11->SetTextFont(63);
	//	tex11->Draw();
     
      }

     if(g==0||g==1||g==6||g==7){
	tex31 = new TLatex(textalign,texty2,"p_{T,jet}>120 GeV/c");
	tex31->SetNDC();
	tex31->SetTextSizePixels(tspixels);
	tex31->SetLineWidth(2);
	tex31->Draw();
     }
     if((g>1&&g<6)||g>7){
       tex31 = new TLatex(textalign,texty2,"p_{T, jet1}> 120 GeV/c");
       tex31->SetName("tex31");
       tex31->SetNDC();
       tex31->SetTextSizePixels(tspixels);
       tex31->SetLineWidth(2);
       tex31->Draw();
       tex33 = new TLatex(textalign,texty3,"p_{T, jet2}> 50 GeV/c");
       tex33->SetName("tex33");
       tex33->SetNDC();
       tex33->SetTextSizePixels(tspixels);
       tex33->SetLineWidth(2);
       tex33->Draw();
       tex32 = new TLatex(textalign,texty4,"#Delta#phi_{1,2}> 5#pi/6");
       tex32->SetName("tex32");
       tex32->SetNDC();
       tex32->SetTextSizePixels(tspixels);
       tex32->SetLineWidth(2);
       tex32->Draw();
     }




    }
    break;

  case 2:
    {
      tex20 = new TLatex(textalign,texty1,"3<p_{T}^{assoc}<4 GeV/c");
      tex20->SetName("tex02");
      tex20->SetNDC();
      tex20->SetTextSizePixels(tspixels);
      tex20->SetLineWidth(2);
      tex20->Draw();
	      	   
   
      if(g<6){

     	tex12 = new TLatex(textalign,texty0,"pp 5.3 pb^{-1} (2.76 TeV)");
	tex12->SetName("tex12");
	tex12->SetNDC();
	tex12->SetTextSizePixels(tspixels);
	//	tex12->SetTextFont(63);
	//	tex12->Draw();
	
      }




   tex21 = new TLatex(textalign,texty2,"anti-k_{T} jets, R=0.3");
      tex21->SetName("tex2");
      tex21->SetNDC();
      tex21->SetTextSizePixels(tspixels);
	   	    
      tex22 = new TLatex(textalign,texty3,"|#eta_{jet}|<1.6");
      tex22->SetName("tex22");
      tex22->SetNDC();
      tex22->SetTextSizePixels(tspixels);
      tex21->Draw();
      tex22->Draw();

    }
    break; 
	 
  case 3:
    {
      tex30 = new TLatex(textalign,texty1,"4<p_{T}^{assoc}<8 GeV/c");
      tex30->SetName("tex30");
      tex30->SetNDC();
      tex30->SetTextSizePixels(tspixels);
      tex30->SetLineWidth(2);
      tex30->Draw("same");

    }
    break;
  }

}





void drawlabels_int_pt2(int g, int j){
  
  TLatex *tex00, *tex01, *tex02,*tex03,*tex10, *tex11, *tex12, *tex13, *tex20, *tex21, *tex22, *tex23, *tex30, *tex31, *tex32, *tex33;

  gStyle->SetTextFont(43);
  
  float texty0=1.05;

  switch(j){
  case 0:{
	  
    if(g<6){  tex13 = new TLatex(textalign2,texty0,"CMS Preliminary");
    }else{  tex13 = new TLatex(textalign2,texty0,"CMS Preliminary Simulation PYTHIA+HYDJET"); }
    tex13->SetName("tex11");
    tex13->SetNDC();
    tex13->SetTextSizePixels(tspixels);
    tex13->SetTextFont(63);
    //   tex13->Draw();   	  
 
    tex00 = new TLatex(textalign2,texty1,"Centrality 50-100%");
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->Draw();
  }
   
    break;
 
  case 1:   
    {
      tex10 = new TLatex(textalign,texty1,"Centrality 30-50%");
      tex10->SetName("tex10");
      tex10->SetNDC();
      tex10->SetTextSizePixels(tspixels);
      tex10->Draw("same");
    
      if(g<6){
	tex11 = new TLatex(textalign,texty0,"PbPb 166 #mub^{-1} (2.76 TeV)");
	tex11->SetName("tex11");
	tex11->SetNDC();
	tex11->SetTextSizePixels(tspixels);
	//	tex11->SetTextFont(63);
	//	tex11->Draw();
     
      }

     if(g==0||g==1||g==6||g==7){
	tex31 = new TLatex(textalign,texty2,"p_{T,jet}>120 GeV/c");
	tex31->SetNDC();
	tex31->SetTextSizePixels(tspixels);
	tex31->SetLineWidth(2);
	tex31->Draw();
     }
     if((g>1&&g<6)||g>7){
       tex31 = new TLatex(textalign,texty2,"p_{T, jet1}> 120 GeV/c");
       tex31->SetName("tex31");
       tex31->SetNDC();
       tex31->SetTextSizePixels(tspixels);
       tex31->SetLineWidth(2);
       tex31->Draw();
       tex33 = new TLatex(textalign,texty3,"p_{T, jet2}> 50 GeV/c");
       tex33->SetName("tex33");
       tex33->SetNDC();
       tex33->SetTextSizePixels(tspixels);
       tex33->SetLineWidth(2);
       tex33->Draw();
       tex32 = new TLatex(textalign,texty4,"#Delta#phi_{1,2}> 5#pi/6");
       tex32->SetName("tex32");
       tex32->SetNDC();
       tex32->SetTextSizePixels(tspixels);
       tex32->SetLineWidth(2);
       tex32->Draw();
     }




    }
    break;

  case 2:
    {
      tex20 = new TLatex(textalign,texty1,"Centrality 10-30%");
      tex20->SetName("tex02");
      tex20->SetNDC();
      tex20->SetTextSizePixels(tspixels);
      tex20->SetLineWidth(2);
      tex20->Draw();
	      	   
   
      if(g<6){

     	tex12 = new TLatex(textalign,texty0,"pp 5.3 pb^{-1} (2.76 TeV)");
	tex12->SetName("tex12");
	tex12->SetNDC();
	tex12->SetTextSizePixels(tspixels);
	//	tex12->SetTextFont(63);
	//	tex12->Draw();
	
      }




   tex21 = new TLatex(textalign,texty2,"anti-k_{T} jets, R=0.3");
      tex21->SetName("tex2");
      tex21->SetNDC();
      tex21->SetTextSizePixels(tspixels);
	   	    
      tex22 = new TLatex(textalign,texty3,"|#eta_{jet}|<1.6");
      tex22->SetName("tex22");
      tex22->SetNDC();
      tex22->SetTextSizePixels(tspixels);
      tex21->Draw();
      tex22->Draw();

    }
    break; 
	 
  case 3:
    {
      tex30 = new TLatex(textalign,texty1,"Centrality 0-10%");
      tex30->SetName("tex30");
      tex30->SetNDC();
      tex30->SetTextSizePixels(tspixels);
      tex30->SetLineWidth(2);
      tex30->Draw("same");

    }
    break;
  }

}





void drawlabels_int_pt(int g, int j){
  
  TLatex *tex00, *tex01, *tex02,*tex03,*tex10, *tex11, *tex12, *tex13, *tex20, *tex21, *tex22, *tex23, *tex30, *tex31, *tex32, *tex33;

 
  switch(j){
  case 0:{
	     	   
    tex00 = new TLatex(textalign2,texty1,"Centrality 50-100%");
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->SetLineWidth(2);
    tex00->Draw("same");

    tex01 = new TLatex(textalign2,texty2,"CMS Preliminary");
    tex01->SetName("tex01");
    tex01->SetNDC();
    tex01->SetTextSizePixels(tspixels);
    tex01->SetLineWidth(2);
    tex01->Draw();
  }
    break;
  case 1:   
    {
      tex10 = new TLatex(textalign,texty1,"Centrality 30-50%");
      tex10->SetName("tex10");
      tex10->SetNDC();
      tex10->SetTextSizePixels(tspixels);
      tex10->SetLineWidth(2);
      tex10->Draw("same");
      tex11 = new TLatex(textalign,texty2,"#sqrt{S_{NN}} = 2.76 TeV");
      tex11->SetName("tex11");
      tex11->SetNDC();
      tex11->SetTextSizePixels(tspixels);
      tex11->SetLineWidth(2);
      tex11->Draw();
      tex12 = new TLatex(textalign,texty3,"PbPb #int L dt = 150 #mub^{-1}");
      tex12->SetName("tex12");
      tex12->SetNDC();
      tex12->SetTextSizePixels(tspixels);
      tex12->SetLineWidth(2);
      tex12->Draw();
      tex13 = new TLatex(textalign,texty4,"pp #int L dt = 5.3 pb^{-1}");
      tex13->SetName("tex11");
      tex13->SetNDC();
      tex13->SetTextSizePixels(tspixels);
      tex13->SetLineWidth(2);
      tex13->Draw();
    }
    break;

  case 2:
    {
      tex20 = new TLatex(textalign,texty1,"Centrality 10-30%");
      tex20->SetName("tex02");
      tex20->SetNDC();
      tex20->SetTextSizePixels(tspixels);
      tex20->SetLineWidth(2);
	      	   
      tex21 = new TLatex(textalign,texty2,"anti-k_{T} jets, R=0.3");
      tex21->SetName("tex2");
      tex21->SetNDC();
      tex21->SetTextSizePixels(tspixels);
	   	    
      tex22 = new TLatex(textalign,texty3,"|#eta|<1.6");
      tex22->SetName("tex22");
      tex22->SetNDC();
      tex22->SetTextSizePixels(tspixels);
      tex20->Draw();
      tex21->Draw();
      tex22->Draw();
    }
    break; 
	 
  case 3:
    {
      tex30 = new TLatex(textalign,texty1,"Centrality 0-10%");
      tex30->SetName("tex30");
      tex30->SetNDC();
      tex30->SetTextSizePixels(tspixels);
      tex30->SetLineWidth(2);
      tex30->Draw("same");
      if(g==1){
	tex31 = new TLatex(textalign,texty2,"p_{T,jet}>120 GeV/c");
	tex31->SetNDC();
	tex31->SetTextSizePixels(tspixels);
	tex31->SetLineWidth(2);
	tex31->Draw();
      }
      if(g==5||g==3){
	tex31 = new TLatex(textalign,texty2,"p_{T,1}>120 GeV/c, p_{T,2}> 50GeV/c");
	tex31->SetName("tex31");
	tex31->SetNDC();
	tex31->SetTextSizePixels(tspixels);
	tex31->SetLineWidth(2);
	tex31->Draw();
	tex32 = new TLatex(textalign,texty3,"#Delta#phi_{1,2}> 5#pi/6");
	tex32->SetName("tex32");
	tex32->SetNDC();
	tex32->SetTextSizePixels(tspixels);
	tex32->SetLineWidth(2);
	tex32->Draw();
      }
    }
    break;
  }
}




void drawlabels_int_pt_mc(int g, int j){
  
  TLatex *tex00, *tex01, *tex02,*tex03,*tex10, *tex11, *tex12, *tex13, *tex20, *tex21, *tex22, *tex23, *tex30, *tex31, *tex32, *tex33;
  gStyle->SetTextFont(423);

 if(g>6){cout<<"monte carlo"<<endl; }

  switch(j){
  case 0:{
	     	   
    tex00 = new TLatex(textalign2,texty1,"Centrality 50-100%");
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->SetLineWidth(2);
    tex00->Draw("same");

    tex01 = new TLatex(textalign2,texty2,"CMS Preliminary");
    tex01->SetName("tex01");
    tex01->SetNDC();
    tex01->SetTextSizePixels(tspixels);
    tex01->SetLineWidth(2);
    tex01->Draw();
  }
    break;
  case 1:   
    {
      tex10 = new TLatex(textalign,texty1,"Centrality 30-50%");
      tex10->SetName("tex10");
      tex10->SetNDC();
      tex10->SetTextSizePixels(tspixels);
      tex10->SetLineWidth(2);
      tex10->Draw("same");
      tex11 = new TLatex(textalign,texty2,"#sqrt{S_{NN}} = 2.76 TeV");
      tex11->SetName("tex11");
      tex11->SetNDC();
      tex11->SetTextSizePixels(tspixels);
      tex11->SetLineWidth(2);
      tex11->Draw();
      tex12 = new TLatex(textalign,texty3,"(PYTHIA+HYDJET) - PYTHIA");
      tex12->SetName("tex12");
      tex12->SetNDC();
      tex12->SetTextSizePixels(tspixels);
      tex12->SetLineWidth(2);
      tex12->Draw();
     
    }
    break;

  case 2:
    {
      tex20 = new TLatex(textalign,texty1,"Centrality 10-30%");
      tex20->SetName("tex02");
      tex20->SetNDC();
      tex20->SetTextSizePixels(tspixels);
      tex20->SetLineWidth(2);
	      	   
      tex21 = new TLatex(textalign,texty2,"anti-k_{T} jets, R=0.3");
      tex21->SetName("tex2");
      tex21->SetNDC();
      tex21->SetTextSizePixels(tspixels);
	   	    
      tex22 = new TLatex(textalign,texty3,"|#eta|<1.6");
      tex22->SetName("tex22");
      tex22->SetNDC();
      tex22->SetTextSizePixels(tspixels);
      tex20->Draw();
      tex21->Draw();
      tex22->Draw();
    }
    break; 
	 
  case 3:
    {
      tex30 = new TLatex(textalign,texty1,"Centrality 0-10%");
      tex30->SetName("tex30");
      tex30->SetNDC();
      tex30->SetTextSizePixels(tspixels);
      tex30->SetLineWidth(2);
      tex30->Draw("same");
      
      tex31 = new TLatex(textalign,texty2,"p_{T,jet}>120 GeV/c");
      tex31->SetNDC();
      tex31->SetTextSizePixels(tspixels);
      tex31->SetLineWidth(2);
      tex31->Draw();
    
    }
    break;
  }
}




void drawlabels_int_cent(int g, int i){

  gStyle->SetTextFont(423);
  
  TLatex *tex00, *tex01, *tex02,*tex03,*tex10, *tex11, *tex12, *tex13, *tex20, *tex21, *tex22, *tex23, *tex30, *tex31, *tex32, *tex33;
 
  texty1+=-0.02;
  texty2+=-0.02;
  texty3+=-0.02;
  texty4+=-0.02;
  
  

  textalign = 0.18;
  
  switch(i){
  case 0:
    tex00 = new TLatex(textalign,texty1,"1<p_{T}^{assoc}<2GeV/c");	      
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->SetLineWidth(2);
    tex00->Draw("same");
    tex01 = new TLatex(textalign,texty2,"CMS Preliminary");
    tex01->SetName("tex01");
    tex01->SetNDC();
    tex01->SetTextSizePixels(tspixels);
    tex01->SetLineWidth(2);
    tex01->Draw();
    break;
  case 1:   
    {
      tex10 = new TLatex(textalign,texty1,"2<p_{T}^{assoc}<3GeV/c");
      tex10->SetName("tex10");
      tex10->SetNDC();
      tex10->SetTextSizePixels(tspixels);
      tex10->SetLineWidth(2);
      tex10->Draw("same");
      tex11 = new TLatex(textalign,texty2,"#sqrt{S_{NN}} = 2.76 TeV");
      tex11->SetName("tex11");
      tex11->SetNDC();
      tex11->SetTextSizePixels(tspixels);
      tex11->SetLineWidth(2);
      tex11->Draw();
      tex12 = new TLatex(textalign,texty3,"PbPb #int L dt = 150 #mub^{-1}");
      tex12->SetName("tex12");
      tex12->SetNDC();
      tex12->SetTextSizePixels(tspixels);
      tex12->SetLineWidth(2);
      tex12->Draw();
      tex13 = new TLatex(textalign,texty4,"pp #int L dt = 5.3 pb^{-1}");
      tex13->SetName("tex11");
      tex13->SetNDC();
      tex13->SetTextSizePixels(tspixels);
      tex13->SetLineWidth(2);
      tex13->Draw();
    }
    break;
    
  case 2:
    {
      tex20 = new TLatex(textalign,texty1,"3<p_{T}^{assoc}<4GeV/c");
      tex20->SetName("tex02");
      tex20->SetNDC();
      tex20->SetTextSizePixels(tspixels);
      tex20->SetLineWidth(2);
	      	   
      tex21 = new TLatex(textalign,texty2,"anti-k_{T} jets, R=0.3");
      tex21->SetName("tex2");
      tex21->SetNDC();
      tex21->SetTextSizePixels(tspixels);
	   	    
      tex22 = new TLatex(textalign,texty3,"|#eta|<1.6");
      tex22->SetName("tex22");
      tex22->SetNDC();
      tex22->SetTextSizePixels(tspixels);
      tex20->Draw();
      tex21->Draw();
      tex22->Draw();
    }
    break; 
	 
  case 3:
    {
      tex30 = new TLatex(textalign,texty1,"4<p_{T}^{assoc}<8GeV/c");
      tex30->SetName("tex30");
      tex30->SetNDC();
      tex30->SetTextSizePixels(tspixels);
      tex30->SetLineWidth(2);
      tex30->Draw("same");
      if(g==1){
	tex31 = new TLatex(textalign,texty2,"p_{T,jet}>120 GeV/c");
	tex31->SetNDC();
	tex31->SetTextSizePixels(tspixels);
	tex31->SetLineWidth(2);
	tex31->Draw();
      }
      if(g==3||g==5){
	tex31 = new TLatex(textalign,texty2,"p_{T,1}>120 GeV/c, p_{T,2}>50 GeV/c");
	tex31->SetName("tex31");
	tex31->SetNDC();
	tex31->SetTextSizePixels(tspixels);
	tex31->SetLineWidth(2);
	tex31->Draw();
	tex32 = new TLatex(textalign,texty3,"#Delta#phi_{1,2}> 5#pi/6");
	tex32->SetName("tex32");
	tex32->SetNDC();
	tex32->SetTextSizePixels(tspixels);
	tex32->SetLineWidth(2);
	tex32->Draw();
      }
    }
    break;
  }

}


void drawlabels_4by4(int g, int i, int j,TString datalabel){

  TLatex *tex00, *tex01, *tex02,*tex03,*tex10, *tex11, *tex12, *tex13, *tex20, *tex21, *tex22, *tex23, *tex30, *tex31, *tex32, *tex33;
 
   
  switch(j){
  case 0:{
	     	   
    tex00 = new TLatex(textalign2,texty1,"Centrality 50-100%");
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->SetLineWidth(2);
   
    tex00->Draw("same");
	    
   
    switch(i){
    case 0:
      tex01 = new TLatex(textalign2,texty2,"0.5<p_{T}^{assoc}<1GeV/c");
      break;
    case 1:
      tex01 = new TLatex(textalign2,texty2,"1<p_{T}^{assoc}<2GeV/c");
      break;
    case 2:
      tex01 = new TLatex(textalign2,texty2,"2<p_{T}^{assoc}<3GeV/c");
      break;
    case 3:
      tex01 = new TLatex(textalign2,texty2,"3<p_{T}^{assoc}<4GeV/c");
      break;
    case 4:
      tex01 = new TLatex(textalign2,texty2,"4<p_{T}^{assoc}<8GeV/c");
      break;
    case 5:
      tex01 = new TLatex(textalign2,texty2,"p_{T}^{assoc}>8GeV/c");
      break;

    }
	      
    tex00->SetName("tex00");
    tex00->SetNDC();
    tex00->SetTextSizePixels(tspixels);
    tex00->SetLineWidth(2);
    tex00->Draw("same");
    tex01->SetName("tex01");
    tex01->SetNDC();
    tex01->SetTextSizePixels(tspixels);
    tex01->SetLineWidth(2);
    tex01->Draw("same");

  
    TString datalabeltext = datalabel;
    datalabeltext+= " Jets";
    if(g==12){datalabeltext = datalabel;}
    tex02 = new TLatex(textalign2,texty3,datalabeltext);
    tex02->SetName("tex02");
    tex02->SetNDC();
    tex02->SetTextSizePixels(tspixels);
    tex02->SetLineWidth(2);
    tex02->Draw();
    if(i==3){
      tex00->SetTextSize(ts3);
      tex01->SetTextSize(ts3);
      tex02->SetTextSize(ts3);
    }
  }

  
    break;
  case 1:   
    {
      tex10 = new TLatex(textalign,texty1,"Centrality 30-50%");
      tex10->SetName("tex10");
      tex10->SetNDC();
      tex10->SetTextSizePixels(tspixels);
      tex10->SetLineWidth(2);
      tex10->Draw("same");
      tex11 = new TLatex(textalign,texty2,"anti-k_{T} jets, R=0.3");
      tex11->SetName("tex11");
      tex11->SetNDC();
      tex11->SetTextSizePixels(tspixels);
      tex11->SetLineWidth(2);
      tex11->Draw();
      tex13 = new TLatex(textalign,texty3,"|#eta| < 1.6");
      tex13->SetName("tex13");
      tex13->SetNDC();
      tex13->SetTextSizePixels(tspixels);
      tex13->SetLineWidth(2);
      tex13->Draw();
      if(i==3){
	tex10->SetTextSizePixels(tspixels);
	tex11->SetTextSizePixels(tspixels);
	tex13->SetTextSizePixels(tspixels);
      }
    }
    break;

  case 2:
    {
      tex20 = new TLatex(textalign,texty1,"Centrality 10-30%");
      //	  if((g==6||g==10)&&i==0){tex20 = new TLatex(textalign,texty1,"Centrality 0-30%");}
      tex20->SetName("tex02");
      tex20->SetNDC();
      tex20->SetTextSizePixels(tspixels);
      tex20->SetLineWidth(2);
      tex20->Draw("same");
	   
      if(g<2||(g>5&&g<8)||g==12){  
	tex21 = new TLatex(textalign,texty2,"p_{T,jet} > 120 GeV/c");
	tex21->SetName("tex2");
	tex21->SetNDC();
	tex21->SetTextSizePixels(tspixels);
	tex21->SetLineWidth(2);
	tex21->Draw();
	if(i==3){
	  tex20->SetTextSizePixels(tspixels);
	  tex21->SetTextSizePixels(tspixels);
	}
      }
      if((g>1&&g<6)||(g>7&&g<12)){  
	tex21 = new TLatex(textalign,texty2,"p_{T,1} > 120 GeV/c");
	tex22 = new TLatex(textalign,texty3,"p_{T,2} > 50 GeV/c");
	tex23 = new TLatex(textalign,texty4,"|#Delta#phi_{1,2}| > 5#pi/6");
	tex21->SetNDC();
	tex21->SetTextSizePixels(tspixels);
	tex21->SetLineWidth(2);
	tex21->Draw();	   
	tex22->SetName("tex22");
	tex22->SetNDC();
	tex22->SetTextSizePixels(tspixels);
	tex22->SetLineWidth(2);
	tex22->Draw();
	tex23->SetName("tex23");
	tex23->SetNDC();
	tex23->SetTextSizePixels(tspixels);
	tex23->SetLineWidth(2);
	tex23->Draw();
	if(i==3){
	  tex20->SetTextSizePixels(tspixels);
	  tex21->SetTextSizePixels(tspixels);
	  tex22->SetTextSizePixels(tspixels);
	  tex23->SetTextSizePixels(tspixels);
	}
      }
    }
    break; 
  case 3:
    {
      tex30 = new TLatex(textalign,texty1,"Centrality 0-10%");
      //if((g==6||g==10)&&i==0){tex30 = new TLatex(textalign,texty1,"Centrality 0-30%");}
      tex30->SetName("tex30");
      tex30->SetNDC();
      tex30->SetTextSizePixels(tspixels);
      tex30->SetLineWidth(2);
      tex30->Draw("same");
      tex31 = new TLatex(textalign,texty2,"CMS Preliminary");
      tex31->SetName("tex31");
      tex31->SetNDC();
      tex31->SetTextSizePixels(tspixels);
      tex31->SetLineWidth(2);
      tex31->Draw();
      tex32 = new TLatex(textalign,texty3,etarangelabel);
      tex32->SetName("tex32");
      tex32->SetNDC();
      tex32->SetTextSizePixels(tspixels);
      tex32->SetLineWidth(2);
      tex32->Draw();
      if(i==3){
	tex30->SetTextSizePixels(tspixels);
	tex31->SetTextSizePixels(tspixels);
	tex32->SetTextSizePixels(tspixels);
      }
    }
    break;
  }
}





void drawlabels_PAS_bkg(int i){
  
  TLatex *tex00, *tex01, *tex02,*tex03,*tex04;

  gStyle->SetTextFont(43);

  tex04 = new TLatex(0.12,0.87,"PbPb 166 #mub^{-1} (2.76 TeV)");
  tex04->SetName("tex02");
  tex04->SetNDC();
  // tex04->SetTextFont(63);
  tex04->SetTextSizePixels(tspixels);
  tex04->Draw();

 
  tex01 = new TLatex(0.12,0.92,"CMS Preliminary");
  tex01->SetName("tex02");
  tex01->SetNDC();
  tex01->SetTextFont(63);
  tex01->SetTextSizePixels(tspixels);
    

  tex03 = new TLatex(0.12,0.82,"p_{T,jet}>120 GeV/c");
  tex03->SetName("tex02");
  tex03->SetNDC();
  tex03->SetTextSizePixels(tspixels);
    



  switch(i){
  case 0: tex00 = new TLatex(0.64,0.82,"0.5<p_{T}^{assoc}<1 GeV/c");break;
  case 1: tex00 = new TLatex(0.64,0.82,"1<p_{T}^{assoc}<2 GeV/c");break;
  case 2: tex00 = new TLatex(0.64,0.82,"2<p_{T}^{assoc}<3 GeV/c");break;
  case 3: tex00 = new TLatex(0.64,0.82,"3<p_{T}^{assoc}<4 GeV/c");break;
  case 4: tex00 = new TLatex(0.64,0.82,"4<p_{T}^{assoc}<8 GeV/c");break;
  case 5: tex00 = new TLatex(0.64,0.82,"p_{T}^{assoc}>8 GeV/c");break;
  }

  tex00->SetName("tex00");
  tex00->SetNDC();
  tex00->SetTextSizePixels(tspixels);
 
  tex02 = new TLatex(0.66,0.87,"Centrality 0-10%");
  tex02->SetName("tex02");
  tex02->SetNDC();
  tex02->SetTextSizePixels(tspixels);
   
  
  tex00->Draw();
  tex01->Draw();
  tex02->Draw();
  tex03->Draw();

}

