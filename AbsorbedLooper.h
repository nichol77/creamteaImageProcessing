//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 31 11:02:09 2009 by ROOT version 5.20/00
// from TTree Absorbed/Absorbed
// found on file: ../data/strips_650/fakecontainer_10cmtargetat_0p5_1_0p5_steelboxat_m0p5_3_m0p5/pca_fakecontainer_10cmtarget_steelbox_million_1.root
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

   // List of branches
   TBranch        *b_gradX;   //!
   TBranch        *b_cutX;   //!
   TBranch        *b_gradY;   //!
   TBranch        *b_cutY;   //!
   TBranch        *b_xyzFitQual;   //!

   AbsorbedLooper(TTree *tree=0);
   virtual ~AbsorbedLooper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     MakeSliceHists(int numBins=20);
   virtual void     MakeSliceHistsIteratively(int binWidth=1000);
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../data/strips_650/fakecontainer_10cmtargetat_0p5_1_0p5_steelboxat_m0p5_3_m0p5/pca_fakecontainer_10cmtarget_steelbox_million_1.root");
      if (!f) {
         f = new TFile("../data/strips_650/fakecontainer_10cmtargetat_0p5_1_0p5_steelboxat_m0p5_3_m0p5/pca_fakecontainer_10cmtarget_steelbox_million_1.root");
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

   fChain->SetBranchAddress("xGrad", &xGrad, &b_gradX);
   fChain->SetBranchAddress("xCut", &xCut, &b_cutX);
   fChain->SetBranchAddress("yGrad", &yGrad, &b_gradY);
   fChain->SetBranchAddress("yCut", &yCut, &b_cutY);
   fChain->SetBranchAddress("xyzFitQual", &xyzFitQual, &b_xyzFitQual);
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
  Int_t detLevel=-7000;
  Int_t cutAt=5000;
  if(TMath::Abs((xGrad*detLevel)+xCut)>cutAt) return -1;
  if(TMath::Abs((yGrad*detLevel)+yCut)>cutAt) return -1;
  if(xyzFitQual>1) return -1;
  return 1;
}
#endif // #ifdef AbsorbedLooper_cxx
