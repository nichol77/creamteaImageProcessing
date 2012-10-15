//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep  7 11:26:17 2009 by ROOT version 5.20/00
// from TTree Absorbed/Absorbed
// found on file: /unix/anita1/creamtea/minerva/fakecontainer_10cmtarget/pca/pca_fakecontainer_10cmtarget_million_1.root
//////////////////////////////////////////////////////////

#ifndef AbsorbedLooper_h
#define AbsorbedLooper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class AbsorbedLooper {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        xGrad;
   Double_t        xCut;
   Double_t        yGrad;
   Double_t        yCut;
   Double_t        xyzFitQual;
   Double_t        xGradReco;
   Double_t        xCutReco;
   Double_t        yGradReco;
   Double_t        yCutReco;
   Double_t        xyzFitQualReco;

   // List of branches
   TBranch        *b_xGrad;   //!
   TBranch        *b_xCut;   //!
   TBranch        *b_yGrad;   //!
   TBranch        *b_yCut;   //!
   TBranch        *b_xyzFitQual;   //!
   TBranch        *b_xGradReco;   //!
   TBranch        *b_xCutReco;   //!
   TBranch        *b_yGradReco;   //!
   TBranch        *b_yCutReco;   //!
   TBranch        *b_xyzFitQualReco;   //!

   AbsorbedLooper(TTree *tree=0);
   virtual ~AbsorbedLooper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    CutReco(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     MakeSliceHists(int numBins=100);
   virtual void     MakeSliceHistsIteratively(int binWidth=100);
   virtual void     MakeSliceHistsIterativelyReco(int binWidth);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AbsorbedLooper_cxx
AbsorbedLooper::AbsorbedLooper(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/unix/creamtea/sfayer/1mdetector_10cmtarget/pca/pca_small1mdetector_10cmtarget_million_1.root");
      if (!f) {
         f = new TFile("/unix/creamtea/sfayer/1mdetector_10cmtarget/pca/pca_small1mdetector_10cmtarget_million_1.root");
      }
      tree = (TTree*)gDirectory->Get("Absorbed");

   }
   Init(tree);
}

AbsorbedLooper::~AbsorbedLooper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AbsorbedLooper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AbsorbedLooper::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AbsorbedLooper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("xGrad", &xGrad, &b_xGrad);
   fChain->SetBranchAddress("xCut", &xCut, &b_xCut);
   fChain->SetBranchAddress("yGrad", &yGrad, &b_yGrad);
   fChain->SetBranchAddress("yCut", &yCut, &b_yCut);
   fChain->SetBranchAddress("xyzFitQual", &xyzFitQual, &b_xyzFitQual);
   fChain->SetBranchAddress("xGradReco", &xGradReco, &b_xGradReco);
   fChain->SetBranchAddress("xCutReco", &xCutReco, &b_xCutReco);
   fChain->SetBranchAddress("yGradReco", &yGradReco, &b_yGradReco);
   fChain->SetBranchAddress("yCutReco", &yCutReco, &b_yCutReco);
   fChain->SetBranchAddress("xyzFitQualReco", &xyzFitQualReco, &b_xyzFitQualReco);
   Notify();
}

Bool_t AbsorbedLooper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AbsorbedLooper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AbsorbedLooper::Cut(Long64_t /*entry*/)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  Int_t detLevel=-800;
  Int_t cutAt=500;
  if(TMath::Abs((xGrad*detLevel)+xCut)>cutAt) return -1;
  if(TMath::Abs((yGrad*detLevel)+yCut)>cutAt) return -1;
  if(xyzFitQual>1) return -1;
  return 1;
}

Int_t AbsorbedLooper::CutReco(Long64_t /*entry*/)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  Int_t detLevel=-800;
  Int_t cutAt=500;
  if(TMath::Abs((xGradReco*detLevel)+xCutReco)>cutAt) return -1;
  if(TMath::Abs((yGradReco*detLevel)+yCutReco)>cutAt) return -1;
  if(xyzFitQual>1e3) return -1; //This cut needs to be examined
  return 1;
}
#endif // #ifdef AbsorbedLooper_cxx
