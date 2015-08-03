//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar 15 23:00:38 2014 by ROOT version 5.32/00
// from TTree mixing_tree/
// found on file: /afs/cern.ch/work/p/pkurt/private/datasets/new_reco_minintuples/pbpb_data/NewReco_PbPbData_mini_ntuple_v1_p0.root
//////////////////////////////////////////////////////////

#ifndef mixing_tree_h
#define mixing_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream>
using namespace std;

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class mixing_tree {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Int_t           mult;
  vector<float>   *pt;
  vector<float>   *phi;
  vector<float>   *eta;
  vector<int>     *chg;
  vector<int>     *sube;
  Int_t           nTrk;
  Int_t           nParticle;
  Int_t           nPFpart;
  vector<float>   *pPt;
  vector<float>   *pPhi;
  vector<float>   *pEta;
  vector<float>   *trkEta;
  vector<float>   *trkPhi;
  vector<float>   *trkPt;
  vector<float>   *trkAlgo;
  vector<float>   *highPurity;
  vector<float>   *vz;
  Int_t           pHBHENoiseFilter;
  Int_t           pcollisionEventSelection;
  Int_t           HLT_HIJet80_v1;
  Int_t           HLT_HIJet80_v7;
  Int_t           HLT_HIMinBiasHfOrBSC_v1;
  Int_t           pPAcollisionEventSelectionPA;
  Int_t           HLT_PAJet80_NoJetID_v1;
  Int_t           hiBin;
  vector<float>   *jteta;
  vector<float>   *jtphi;
  vector<float>   *jtpt;
  vector<float>   *corrpt;
  Float_t         pthat;
  vector<float>   *trackMax;
  vector<float>   *trkDxy1;
  vector<float>   *trkDxyError1;
  vector<float>   *trkDz1;
  vector<float>   *trkDzError1;
  vector<float>   *trkPtError;
  vector<float>   *geneta;
  vector<float>   *genphi;
  vector<float>   *genpt;
 
  vector<Int_t> *pfId;
  vector<float> *pfPt;
  vector<float> *pfVsPt;
  vector<float> *pfEta;
  vector<float> *pfPhi;
  vector<float> *sumpt;

  // List of branches
  TBranch        *b_mult;   //!
  TBranch        *b_pt;   //!
  TBranch        *b_phi;   //!
  TBranch        *b_eta;   //!
  TBranch        *b_chg;   //!
  TBranch        *b_sube;   //!
  TBranch        *b_nTrk;   //!
  TBranch        *b_nPFpart;   //!
  TBranch        *b_nParticle;   //!
  TBranch        *b_pPt;   //!
  TBranch        *b_pPhi;   //!
  TBranch        *b_pEta;   //!
  TBranch        *b_trkEta;   //!
  TBranch        *b_trkPhi;   //!
  TBranch        *b_trkPt;   //!
  TBranch        *b_trkAlgo;   //!
  TBranch        *b_highPurity;   //!
  TBranch        *b_vz;   //!
  TBranch        *b_pHBHENoiseFilter;   //!
  TBranch        *b_pcollisionEventSelection;   //!
  TBranch        *b_HLT_HIJet80_v1;   //!
  TBranch        *b_HLT_HIJet80_v7;   //!
  TBranch        *b_HLT_HIMinBiasHfOrBSC_v1;   //!
  TBranch        *b_pPAcollisionEventSelectionPA;
  TBranch        *b_HLT_PAJet80_NoJetID_v1;
  TBranch        *b_hiBin;   //!
  TBranch        *b_jteta;   //!
  TBranch        *b_jtphi;   //!
  TBranch        *b_jtpt;   //!
  TBranch        *b_corrpt;   //!
  TBranch        *b_pthat;   //!
  TBranch        *b_trackMax;   //!
  TBranch        *b_trkDxy1;   //!
  TBranch        *b_trkDxyError1;   //!
  TBranch        *b_trkDz1;   //!
  TBranch        *b_trkDzError1;   //!
  TBranch        *b_trkPtError;   //!
  TBranch        *b_geneta;   //!
  TBranch        *b_genphi;   //!
  TBranch        *b_genpt;   //!
 

  TBranch        *b_pfId;   //!
  TBranch        *b_pfPt;   //!
  TBranch        *b_pfVsPt;   //!
  TBranch        *b_pfEta;   //!
  TBranch        *b_pfPhi;   //!
  TBranch        *b_sumpt;   //!


  mixing_tree(TTree *tree=0);
  virtual ~mixing_tree();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef mixing_tree_cxx
mixing_tree::mixing_tree(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/work/p/pkurt/private/datasets/new_reco_minintuples/pbpb_data/NewReco_PbPbData_mini_ntuple_v1_p0.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/afs/cern.ch/work/p/pkurt/private/datasets/new_reco_minintuples/pbpb_data/NewReco_PbPbData_mini_ntuple_v1_p0.root");
    }
    f->GetObject("mixing_tree",tree);

  }
  Init(tree);
}

mixing_tree::~mixing_tree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t mixing_tree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t mixing_tree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void mixing_tree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  pt = 0;
  phi = 0;
  eta = 0;
  chg = 0;
  sube = 0;
  pPt = 0;
  pPhi = 0;
  pEta = 0;
  trkEta = 0;
  trkPhi = 0;
  trkPt = 0;
  trkAlgo = 0;
  highPurity = 0;
  vz = 0;
  jteta = 0;
  jtphi = 0;
  jtpt = 0;
  corrpt = 0;
  trackMax = 0;
  trkDxy1 = 0;
  trkDxyError1 = 0;
  trkDz1 = 0;
  trkDzError1 = 0;
  trkPtError = 0;
  geneta = 0;
  genphi = 0;
  genpt = 0;
 
  pfId = 0;
  pfPt = 0;
  pfVsPt = 0;
  pfEta = 0;
  pfPhi = 0;
  sumpt = 0;

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("mult", &mult, &b_mult);
  fChain->SetBranchAddress("pt", &pt, &b_pt);
  fChain->SetBranchAddress("phi", &phi, &b_phi);
  fChain->SetBranchAddress("eta", &eta, &b_eta);
  fChain->SetBranchAddress("chg", &chg, &b_chg);
  fChain->SetBranchAddress("sube", &sube, &b_sube);
  fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
  fChain->SetBranchAddress("nPFpart", &nPFpart, &b_nPFpart);
  fChain->SetBranchAddress("nParticle", &nParticle, &b_nParticle);
  fChain->SetBranchAddress("pPt", &pPt, &b_pPt);
  fChain->SetBranchAddress("pPhi", &pPhi, &b_pPhi);
  fChain->SetBranchAddress("pEta", &pEta, &b_pEta);
  fChain->SetBranchAddress("trkEta", &trkEta, &b_trkEta);
  fChain->SetBranchAddress("trkPhi", &trkPhi, &b_trkPhi);
  fChain->SetBranchAddress("trkPt", &trkPt, &b_trkPt);
  fChain->SetBranchAddress("trkAlgo", &trkAlgo, &b_trkAlgo);
  fChain->SetBranchAddress("highPurity", &highPurity, &b_highPurity);
  fChain->SetBranchAddress("vz", &vz, &b_vz);
  fChain->SetBranchAddress("pHBHENoiseFilter", &pHBHENoiseFilter, &b_pHBHENoiseFilter);
  fChain->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection, &b_pcollisionEventSelection);
  fChain->SetBranchAddress("HLT_HIJet80_v1", &HLT_HIJet80_v1, &b_HLT_HIJet80_v1);
  fChain->SetBranchAddress("HLT_HIJet80_v7", &HLT_HIJet80_v7, &b_HLT_HIJet80_v7);
  fChain->SetBranchAddress("HLT_HIMinBiasHfOrBSC_v1", &HLT_HIMinBiasHfOrBSC_v1, &b_HLT_HIMinBiasHfOrBSC_v1);
  fChain->SetBranchAddress("pPAcollisionEventSelectionPA", &pPAcollisionEventSelectionPA, &b_pPAcollisionEventSelectionPA);
  fChain->SetBranchAddress("HLT_PAJet80_NoJetID_v1", &HLT_PAJet80_NoJetID_v1, &b_HLT_PAJet80_NoJetID_v1);
  fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
  fChain->SetBranchAddress("jteta", &jteta, &b_jteta);
  fChain->SetBranchAddress("jtphi", &jtphi, &b_jtphi);
  fChain->SetBranchAddress("jtpt", &jtpt, &b_jtpt);
  fChain->SetBranchAddress("corrpt", &corrpt, &b_corrpt);
  fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
  fChain->SetBranchAddress("trackMax", &trackMax, &b_trackMax);
  fChain->SetBranchAddress("trkDxy1", &trkDxy1, &b_trkDxy1);
  fChain->SetBranchAddress("trkDxyError1", &trkDxyError1, &b_trkDxyError1);
  fChain->SetBranchAddress("trkDz1", &trkDz1, &b_trkDz1);
  fChain->SetBranchAddress("trkDzError1", &trkDzError1, &b_trkDzError1);
  fChain->SetBranchAddress("trkPtError", &trkPtError, &b_trkPtError);
  fChain->SetBranchAddress("geneta", &geneta, &b_geneta);
  fChain->SetBranchAddress("genphi", &genphi, &b_genphi);
  fChain->SetBranchAddress("genpt", &genpt, &b_genpt);
 
  fChain->SetBranchAddress("pfId", &pfId, &b_pfId);
  fChain->SetBranchAddress("pfPt", &pfPt, &b_pfPt);
  fChain->SetBranchAddress("pfVsPt", &pfVsPt, &b_pfVsPt);
  fChain->SetBranchAddress("pfEta", &pfEta, &b_pfEta);
  fChain->SetBranchAddress("pfPhi", &pfPhi, &b_pfPhi);
  fChain->SetBranchAddress("sumpt", &sumpt, &b_sumpt);



  Notify();
}

Bool_t mixing_tree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void mixing_tree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t mixing_tree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef mixing_tree_cxx
