

  //*********************************************************
  // SET INTEGRAL LIMITS AND LABELS
  //*************************************************

  double etalim = 1.0;
  double philim = 1.0;

  TString etarangelabel = "Projected |#Delta#eta|<1.0";
  TString phirangelabel = "Projected |#Delta#phi|<1.0";

  float llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  //*********************************************************

Int_t nbounds_phi = 16;
Int_t nbounds_eta = 20;
Int_t nbounds_eta2 = 18;
 
Double_t bin_bounds_phi[16]   ={-1.5708, -1.22522  , -1.03673  ,-0.84823,  -0.659734, -0.471239, -0.282743, -0.0942478, 0.0942478,  0.282743,  0.471239,  0.659734,  0.84823,  1.03673, 1.22522,1.5708};

//Double_t bin_bounds_phi[26]   =  {-1.5708 ,-1.44513 ,-1.31947 ,-1.19381 ,-1.06814 ,-0.942478 ,-0.816814 ,-0.69115 ,-0.565487 ,-0.439823 ,-0.314159 ,-0.188496 ,-0.0628319 ,0.0628319 ,0.188496 ,0.314159 ,0.439823 ,0.565487 ,0.69115 ,0.816814 ,0.942478 ,1.06814 ,1.19381 ,1.31947 ,1.44513 ,1.5708};


Double_t bin_bounds_phi_wide[12]   =  {-1.5708 ,-1.19381 ,-0.942478  ,-0.69115 ,-0.439823  ,-0.188496,0.188496 ,0.439823,0.69115 ,0.942478  ,1.19381,1.5708};

//Double_t bin_bounds_phi_wide[12]   =  {-1.5708 ,-1.31947 ,-1.06814 ,-0.816814 ,-0.565487  ,-0.314159 , 0.314159  ,0.565487  ,0.816814,1.06814  ,1.31947 ,1.5708};

Double_t bin_bounds_eta[20]
    =   {-2.5,-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.,2.5};

 Double_t bin_bounds_eta2[18]
    =   {-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.};

   //------------------------
  //Canvas formatting specs
  //-------------------------

  float textalign = 0.05;
  float textalign2 = 0.23;
  float ts = 0.07;
  float ts2 = 0.065;
  float ts3 = 0.057;
  float tstitle = 0.08;
float tstitle2 = 0.08;
  float xlabeloffset2 = 0.5;
  float xoffset = 0.85;
  float xoffset2 = 1.0;
  float yoffset = 0.85;
  float texty1 = 0.9;
  float texty2 = 0.8;
  float texty3 = 0.7;
  float texty4 = 0.6;
  float legendoffset = 0.65;
  float legendoffset2 = 0.7;

float tspixels = 25;

  //------------------------------
  



 

