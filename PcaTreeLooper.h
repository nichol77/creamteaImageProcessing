//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  3 07:14:30 2009 by ROOT version 5.20/00
// from TTree pcaTree/pcaTree
// found on file: /unix/anita1/creamtea/strips_650/fakecontainer_10cmtargetandwaterboxat_0p5_1_0p5_hollowsteelboxat_m0p5_3_m0p5/pca/pca_fakecontainer_10cmtargetandwaterbox_hollowsteelbox_million_1.root
//////////////////////////////////////////////////////////

#ifndef WillTreeLooper_h
#define WillTreeLooper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3.h>
#include <map>

class WillTreeLooper {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::map <int, std::map<int, double> > L;
   std::map <int, double> Lambda; ///< Does it make sense to switch these to arrays?
   //   Int_t fVoxelCountl
   // double *Lambda;
   std::map <int, double> Gradient; ///< Does it make sense to switch these to arrays?
   //double *Gradient;
   std::map <int, double> Sigma; ///< Does it make sense to switch these to arrays?
   std::map <int, double> S; ///< Does it make sense to switch these to arrays?


   // Declaration of leaf types
   Double_t        xPosTrue;
   Double_t        yPosTrue;
   Double_t        zPosTrue;
   Double_t        thetaTrue;
   Double_t        thetaxzTrue;
   Double_t        thetayzTrue;
   Double_t        xzGradTrue[2];
   Double_t        yzGradTrue[2];
   Double_t        xzCutTrue[2];
   Double_t        yzCutTrue[2];
   Double_t        xyzFitQualTrue[2];
   Int_t           gotRecoPCA;
   Double_t        xPosReco;
   Double_t        yPosReco;
   Double_t        zPosReco;
   Double_t        thetaReco;
   Double_t        thetaxzReco;
   Double_t        thetayzReco;
   Double_t        xzGradReco[2];
   Double_t        yzGradReco[2];

   // List of branches
   TBranch        *b_xPosTrue;   //!
   TBranch        *b_yPosTrue;   //!
   TBranch        *b_zPosTrue;   //!
   TBranch        *b_thetaTrue;   //!
   TBranch        *b_thetaxzTrue;   //!
   TBranch        *b_thetayzTrue;   //!
   TBranch        *b_xzGradTrue;   //!
   TBranch        *b_yzGradTrue;   //!
   TBranch        *b_xzCutTrue;   //!
   TBranch        *b_yzCutTrue;   //!
   TBranch        *b_xyzFitQualTrue;   //!
   TBranch        *b_gotRecoPCA;   //!
   TBranch        *b_xPosReco;   //!
   TBranch        *b_yPosReco;   //!
   TBranch        *b_zPosReco;   //!
   TBranch        *b_thetaReco;   //!
   TBranch        *b_thetaxzReco;   //!
   TBranch        *b_thetayzReco;   //!
   TBranch        *b_xzGradReco;   //!
   TBranch        *b_yzGradReco;   //!

   WillTreeLooper(TTree *tree=0);
   virtual ~WillTreeLooper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void SLFill(int first, int last, int Nx, int Ny, int Nz);
   virtual void LambdaFill(int VoxelCount);
   virtual void GradientFill();
   virtual void LambdaAlpha(double Alpha);
   virtual void SigmaFill();
   virtual double Cost(double Alpha, int first, int last);
   virtual void     FillPosHist(TH3F *histPos, Double_t thetaCut=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void DrawSlices(int topSlice, int bottomSlice, int Nx, int Ny, int Nz,char* fileNameLambda);
};

#endif

#ifdef WillTreeLooper_cxx
WillTreeLooper::WillTreeLooper(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/unix/anita1/creamtea/strips_650/fakecontainer_10cmtargetandwaterboxat_0p5_1_0p5_hollowsteelboxat_m0p5_3_m0p5/pca/pca_fakecontainer_10cmtargetandwaterbox_hollowsteelbox_million_1.root");
      if (!f) {
         f = new TFile("/unix/anita1/creamtea/strips_650/fakecontainer_10cmtargetandwaterboxat_0p5_1_0p5_hollowsteelboxat_m0p5_3_m0p5/pca/pca_fakecontainer_10cmtargetandwaterbox_hollowsteelbox_million_1.root");
      }
      tree = (TTree*)gDirectory->Get("pcaTree");

   }
   Init(tree);
}

WillTreeLooper::~WillTreeLooper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t WillTreeLooper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t WillTreeLooper::LoadTree(Long64_t entry)
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

void WillTreeLooper::Init(TTree *tree)
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

   fChain->SetBranchAddress("xPosTrue", &xPosTrue, &b_xPosTrue);
   fChain->SetBranchAddress("yPosTrue", &yPosTrue, &b_yPosTrue);
   fChain->SetBranchAddress("zPosTrue", &zPosTrue, &b_zPosTrue);
   fChain->SetBranchAddress("thetaTrue", &thetaTrue, &b_thetaTrue);
   fChain->SetBranchAddress("thetaxzTrue", &thetaxzTrue, &b_thetaxzTrue);
   fChain->SetBranchAddress("thetayzTrue", &thetayzTrue, &b_thetayzTrue);
   fChain->SetBranchAddress("xzCutTrue", xzCutTrue, &b_xzCutTrue);
   fChain->SetBranchAddress("yzCutTrue", yzCutTrue, &b_yzCutTrue);
   fChain->SetBranchAddress("xzGradTrue", xzGradTrue, &b_xzGradTrue);
   fChain->SetBranchAddress("yzGradTrue", yzGradTrue, &b_yzGradTrue);
   fChain->SetBranchAddress("xyzFitQualTrue", xyzFitQualTrue, &b_xyzFitQualTrue);
   fChain->SetBranchAddress("gotRecoPCA", &gotRecoPCA, &b_gotRecoPCA);
   fChain->SetBranchAddress("xPosReco", &xPosReco, &b_xPosReco);
   fChain->SetBranchAddress("yPosReco", &yPosReco, &b_yPosReco);
   fChain->SetBranchAddress("zPosReco", &zPosReco, &b_zPosReco);
   fChain->SetBranchAddress("thetaReco", &thetaReco, &b_thetaReco);
   fChain->SetBranchAddress("thetaxzReco", &thetaxzReco, &b_thetaxzReco);
   fChain->SetBranchAddress("thetayzReco", &thetayzReco, &b_thetayzReco);
   fChain->SetBranchAddress("xzGradReco", xzGradReco, &b_xzGradReco);
   fChain->SetBranchAddress("yzGradReco", yzGradReco, &b_yzGradReco);
   Notify();
}

Bool_t WillTreeLooper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void WillTreeLooper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t WillTreeLooper::Cut(Long64_t /*entry*/)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef WillTreeLooper_cxx


