#define PcaTreeLooper_cxx
#include "PcaTreeLooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <map>
#include <TMath.h>
//For Calculation Time
#include <time.h>
#include <stdio.h>
#include <vector>
using namespace std;

void PcaTreeLooper::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L PcaTreeLooper.C
//      Root > PcaTreeLooper t
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

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
  }
}

void PcaTreeLooper::FillPosHist(TH3F *histPos, Double_t thetaCut)
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  
  Int_t numStars=100;
  Long64_t starEvery=nentries/numStars;
  if(starEvery==0) starEvery++;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%starEvery==0) std::cerr << "*";
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(thetaTrue<thetaCut) continue;
    if(xyzFitQualTrue[0]>1 || xyzFitQualTrue[1]>1) continue;
    
    histPos->Fill(xPosTrue,yPosTrue,zPosTrue);    
  }
  std::cerr << "\n";
}





///////////////////////////////////////////////////////////////////////Will's Functions
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//Define RootMap
ClassImp(RootMap);

RootMap::RootMap(){};
RootMap::~RootMap(){};





//Fill Matrix L with muons from FIRST to LAST, and make the S vector. Ni is the number of voxels in the i direction.
void PcaTreeLooper::SLFill(int first, int last, int Nx, int Ny, int Nz){

  if(fChain==0) return;
  
  //Calculation Time
  clock_t LijStart = clock();
  
  double S;

  char RootFileOut[FILENAME_MAX];
  sprintf(RootFileOut, "L_%d_%d.root",first,last);
  TFile *fpout = new TFile(RootFileOut,"RECREATE");
  TTree *L_Tree = new TTree("L","Matrix L");
  RootMap *PtrLj = new RootMap();
  L_Tree->Branch("Voxels","RootMap",&PtrLj);
  L_Tree->Branch("ThetaSquared",&S,"ThetaSquared/D");




  
  //Size of cuboid
  double MaxX = 6500;
  double MinX = -6500;
  double MaxY = 6500;
  double MinY = -6500;
  double MaxZ = 6500;
  double MinZ = -6500;



  
  //Voxel dimensions
  double Vx = (MaxX-MinX)/Nx;
  double Vy = (MaxY-MinY)/Ny;
  double Vz = (MaxZ-MinZ)/Nz;


  double VoxelDiag = pow(pow(Vx,2)+pow(Vy,2)+pow(Vz,2),0.5);





  //Loop over Muon Branches of trees in Chain
  for(int jentry = first; jentry < last; jentry++){
    Long64_t ientry = LoadTree((Long64_t)jentry);
    if(ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    special = 0;
    
    //Only use muon if there is a PCA
    if(gotRecoPCA != 1) break;
    
    //Scattering Angle Cuttof
    if(thetaTrue > 2) break;
    
    
    
    
    
    
    
    
    PtrLj->Array.clear();
    for(std::map<int,double>::iterator iter = L[jentry].begin(); iter != L[jentry].end(); ++iter){
      PtrLj->Array[iter->first] = iter->second;
    }
    S_array[jentry] = S = thetaTrue*thetaTrue;
    L_Tree->Fill();
    
    
    
  }//Close loop over muons
  
  
  
  
  
  
  L_Tree->AutoSave();
  
  fpout->Close();
  
  
  printf("Lij & S Fill: %f\n", ((double)(clock() - LijStart) / CLOCKS_PER_SEC));
  
  
  
  
  
  
}//End LSFill


