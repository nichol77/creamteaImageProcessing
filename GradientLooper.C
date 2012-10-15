#define GradientLooper_cxx
#include "GradientLooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void GradientLooper::Loop(int sliceNo)
{
//   In a ROOT session, you can do:
//      Root > .L GradientLooper.C
//      Root > GradientLooper t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

      Long64_t ientry = LoadTree(0);
      if (ientry < 0) break;
      nb = fChain->GetEntry(0);
      // if (Cut(ientry) < 0) continue;


      //Work out Voxel Start and End from Slice Number. 0 = top, 100 = bottom.
      int SliceNo = 50;
      int Nx, Ny, Nz;
      Nx = Ny = Nz = 100;
      int VoxStart = SliceNo * Nx * Ny;
      int VoxEnd = VoxStart + Nx * Ny;

      //Make TH2F
      TH2F XYSlice = new TH2F("XYSlice", "Gradient map of a slice of the Volume", Nx, 130.0, Ny, 130.0);
      
      for(int cX = 0; cX < Nx; cX++){
	for(int cY = 0; cY < Ny; cY++){
	  XYSlice.SetBinContent(cX,cY,Array.second[VoxStart + cY * Nx + cX]);
	}
      }
      




      //TH2F

}
