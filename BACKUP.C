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




 //Initialise L
  typedef map<int, map<int, double> > MapMap; 
  MapMap L;
  
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
   
  //Initialise Variables
  int EndX, EndY, EndZ, EndID, CurrentX, CurrentY, CurrentZ, vID, xzCase, yzCase, xSign, ySign, errorcounter, special;
  errorcounter = 0;
  double xEntry, yEntry, zEntry, xExit, yExit, zExit, zExitX, zExitY, Length;
  Long64_t nbytes = 0;
  Long64_t nb = 0;

  //error margin for float inequalities
  float Epsilon = 0.0000001;
  xSign=ySign=xzCase=yzCase=0;
  xExit=yExit=zExit=0;
  for(int jentry = first; jentry < last; jentry++){
    Long64_t ientry = LoadTree((Long64_t)jentry);
    if(ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    special = 0;
    
    //Only use muon if there is a PCA
    if(gotRecoPCA != 1) break;

    //Scattering Angle Cuttof
    if(thetaTrue > 2) break;

    //Must start and end in the box.
    double xEntryTest, yEntryTest, xExitTest, yExitTest;
    xEntryTest = xzGradTrue[0] * MaxZ + xzCutTrue[0];
    yEntryTest = yzGradTrue[0] * MaxZ + yzCutTrue[0];
    xExitTest = xzGradTrue[1] * MinZ + xzCutTrue[1];
    yExitTest = yzGradTrue[1] * MinZ + yzCutTrue[1];
    if(xEntryTest < MaxX && xEntryTest > MinX && yEntryTest < MaxY && yEntryTest > MinY && xExitTest > MinX && xExitTest < MaxX && yExitTest > MinY && yExitTest < MaxY){
      
      //Is PCA in box?
      if((xPosTrue >= MinX) && (xPosTrue <= MaxX) && (yPosTrue >= MinY) && (yPosTrue <= MaxY) && (zPosTrue >= MinZ) && (zPosTrue <= MaxZ)){
	
	
	for(int part=0; part < 2; part++){
	  
	  
	  if(part == 0){
	    xEntry = xEntryTest;
	    yEntry = yEntryTest;
	    zEntry = MaxZ;
	    EndX = (int) floor((xPosTrue-MinX)/Vx);
	    EndY = (int) floor((yPosTrue-MinY)/Vy);
	    EndZ = (int) floor((zPosTrue-MinZ)/Vz);
	  }
	  else{
	    EndX = (int) floor((xExitTest-MinX)/Vx);
	    EndY = (int) floor((yExitTest-MinY)/Vy);
	    EndZ = 0;
	  }
	  
	  //What Voxel Does it Start in?
	  CurrentX = (int) floor((xEntry-MinX)/Vx);
	  CurrentY = (int) floor((yEntry-MinY)/Vy);
	  CurrentZ = (int) floor((zEntry-MinZ)/Vz);
	  vID = CurrentX + CurrentY*Nx + CurrentZ*Nx*Ny;
	  
	  //Which Voxel is this muon track's last?
	  EndID = EndX + EndY*Nz + EndZ*Nx*Ny;
	  
	  
	  if(xzGradTrue[part] > 0){
	    xzCase=0;//If gradient is positive, then x will leave via the CurrentX*Vx Wall.
	    xSign=-1;
	  }else if(xzGradTrue[part] < 0){
	    xzCase=1;//Same as above, but (CurrentX+1)*Vx Wall.
	    xSign=1;
	  }
	  
	  if(yzGradTrue[part] > 0){
	    yzCase=0;
	    ySign=-1;
	  }else if(yzGradTrue[part] < 0){
	    yzCase=1;
	    ySign=1;
	  }
	  
	  //Loop over voxels for this muon track.
	  while((vID != EndID) && (special < 2*Nz)){
	    zExitX = ((CurrentX+xzCase)*Vx+MinX-xzCutTrue[part])/xzGradTrue[part];
	    zExitY = ((CurrentY+yzCase)*Vy+MinY-yzCutTrue[part])/yzGradTrue[part];
	    
	    zExit = max(max(zExitX,zExitY),CurrentZ*Vz+MinZ);
	    xExit = xzGradTrue[part]*zExit+xzCutTrue[part];
	    yExit = yzGradTrue[part]*zExit+yzCutTrue[part];
	    
	    Length = pow(pow(xExit-xEntry,2)+pow(yExit-yEntry,2)+pow(zExit-zEntry,2),0.5);
	    if(Length < VoxelDiag){L[jentry][vID] = Length;}
	    else{errorcounter++;}
	    
	    if((zExit <= zExitX+Epsilon) && (zExit >= zExitX-Epsilon)){CurrentX += xSign;}
	    if((zExit <= zExitY+Epsilon) && (zExit >= zExitY-Epsilon)){CurrentY += ySign;}
	    if((zExit <= CurrentZ*Vz+MinZ+ Epsilon) && (zExit >= CurrentZ*Vz+MinZ-Epsilon)){CurrentZ -= 1;}
	    
	    xEntry = xExit;
	    yEntry = yExit;
	    zEntry = zExit;
	    vID = CurrentX + CurrentY*Nx + CurrentZ*Nx*Ny;
	  if(vID < 0) cout << vID << endl;
	    special++;
	    
	  }
	  
	  
	  if(vID == EndID && part == 0){ //PCA Voxel. Two parts, up to PCA, then up to exit.
	    
	    //Up to PCA
	    Length = pow(pow(xExit-xEntry,2)+pow(yExit-yEntry,2)+pow(zExit-zEntry,2),0.5);
	    if(Length < VoxelDiag){L[jentry][vID] = Length;}
	    else{errorcounter++;}
	    
	    //Other side of PCA
	    zExitX = (xPosTrue-xzCutTrue[part])/xzGradTrue[part];
	    zExitY = (yPosTrue-yzCutTrue[part])/yzGradTrue[part];
	    
	    zExit = max(max(zExitX,zExitY),CurrentZ*Vz+MinZ);
	    xExit = xzGradTrue[part]*zExit+xzCutTrue[part];
	    yExit = yzGradTrue[part]*zExit+yzCutTrue[part];
	    
	    //Add on second part of length in PCA Voxel
	    Length = pow(pow(xExit-xEntry,2)+pow(yExit-yEntry,2)+pow(zExit-zEntry,2),0.5);
	    if(Length < VoxelDiag){L[jentry][vID] += Length;}
	    else{errorcounter++;}
	    
	    if((zExit <= zExitX+Epsilon) && (zExit >= zExitX-Epsilon)){CurrentX += xSign;}
	    if((zExit <= zExitY+Epsilon) && (zExit >= zExitY-Epsilon)){CurrentY += ySign;}
	    if((zExit <= CurrentZ*Vz+MinZ+ Epsilon) && (zExit >= CurrentZ*Vz+MinZ-Epsilon)){CurrentZ -= 1;}
	    
	    xEntry = xExit;
	    yEntry = yExit;
	    zEntry = zExit;
	    vID = CurrentX + CurrentY*Nx + CurrentZ*Nx*Ny;
	    
	  }else if(vID == EndID && part == 1){
	    //Clause for final voxel. Do it here rather than in the main loop so the PCA voxel is easier to do.
	    zExit = MinZ;
	    xExit = xzGradTrue[part]*zExit+xzCutTrue[part];
	    yExit = yzGradTrue[part]*zExit+yzCutTrue[part];
	    Length = pow(pow(xExit-xEntry,2)+pow(yExit-yEntry,2)+pow(zExit-zEntry,2),0.5);
	    if(Length < VoxelDiag){L[jentry][vID] += Length;}
	    else{errorcounter++;}	  
	  }
	  
	  
	  
	  
	  
	}//Clost loop over two parts of muon track
	
	PtrLj->Array.clear();
	for(std::map<int,double>::iterator iter = L[jentry].begin(); iter != L[jentry].end(); ++iter){
	  PtrLj->Array[iter->first] = iter->second;
	}
	S_array[jentry] = S = thetaTrue*thetaTrue;
	L_Tree->Fill();
      }//Close if for is in PCA box? 
    }
    
    
    
  }//Close loop over muons

  




  L_Tree->AutoSave();

  fpout->Close();


  printf("Lij & S Fill: %f\n", ((double)(clock() - LijStart) / CLOCKS_PER_SEC));





  
}//End LFill




void PcaTreeLooper::LambdaFill(int VoxelCount){
  clock_t LambdaStart = clock();
  
  typedef map<int, map<int, double> > MapMap; 


  for(int i=0; i < VoxelCount; i++){
    Lambda[i] = 0.000003067;
  }


  printf("Lambda Fill: %f\n", ((double)(clock() - LambdaStart) / CLOCKS_PER_SEC));
}



//new Lambda Vector
void PcaTreeLooper::LambdaAlpha(double Alpha){
  clock_t LambdaStart = clock();
  /*
  //Get Old Lambda
  //Get Gradient
  //Get Alpha*
  */

  double AlphaStar;
  double newLambda;


  for(std::map<int,double>::iterator iter = Lambda.begin(); iter != Lambda.end(); ++iter){
    newLambda = iter->second - Alpha*Gradient[iter->first];
    if(newLambda < 0.000003067) newLambda = 0.000003067;
    Lambda[iter->first] = newLambda;
  }



  printf("Lambda Alpha Fill: %f\n", ((double)(clock() - LambdaStart) / CLOCKS_PER_SEC));

}






//Gradient Vector.
void PcaTreeLooper::GradientFill(int first, int last){
  clock_t GradientStart = clock();
  
  char RootFileOut[FILENAME_MAX];
  sprintf(RootFileOut, "L_%d_%d.root",first,last);
  TFile *fpin = new TFile(RootFileOut,"OLD");
  TTree *L_Tree = (TTree*) fpin->Get("L");
  RootMap *PtrLj = 0;
  L_Tree->SetBranchAddress("Voxels",&PtrLj);
  double S;
  L_Tree->SetBranchAddress("ThetaSquared",&S);


  double sum = 0;
  Gradient.clear();
  Long64_t total = L_Tree->GetEntries();
  for(Long64_t i = 0; i < total; ++i){
    L_Tree->GetEntry(i);
    for(std::map<int,double>::iterator iter = PtrLj->Array.begin(); iter != PtrLj->Array.end(); ++iter){
      Gradient[iter->first] += iter->second * (1 + (S_array[iter->first]/Sigma[iter->first]));
      sum = iter->second * (1 + (S_array[iter->first]/Sigma[iter->first]));  
    }
    
    
  }
  printf("Gradient Fill: %f\n", ((double)(clock() - GradientStart) / CLOCKS_PER_SEC));
}



void PcaTreeLooper::SigmaFill(int first, int last){
  
  clock_t SigmaStart = clock();

  char RootFileOut[FILENAME_MAX];
  sprintf(RootFileOut, "L_%d_%d.root",first,last);
  TFile *fpin = new TFile(RootFileOut,"OLD");
  TTree *L_Tree = (TTree*) fpin->Get("L");
  RootMap *PtrLj = 0;
  L_Tree->SetBranchAddress("Voxels",&PtrLj);
  double newLambda;

  double prefactor = 13.6*13.6/1000000;
  Long64_t total = L_Tree->GetEntries();
  for(Long64_t i = 0; i < total; ++i){ //loop over muons
    L_Tree->GetEntry(i);
    Sigma[i] = 0;
    for(std::map<int,double>::iterator iter = PtrLj->Array.begin(); iter != PtrLj->Array.end(); ++iter){
      Sigma[i] += iter->second * Lambda[iter->first];
    
    }
    Sigma[i] *= prefactor;
  }
  
  fpin->Close();
  printf("Sigma Fill: %f\n", ((double)(clock() - SigmaStart) / CLOCKS_PER_SEC));
}



//cost Function
double PcaTreeLooper::Cost(double Alpha){
  double cost=0;
  char RootFileOut[FILENAME_MAX];
  sprintf(RootFileOut, "L_0_10000.root");
  TFile *fpin = new TFile(RootFileOut,"OLD");
  TTree *L_Tree = (TTree*) fpin->Get("L");
  RootMap *PtrLj = 0;
  L_Tree->SetBranchAddress("Voxels",&PtrLj);
  Long64_t total = L_Tree->GetEntries();

  map<int,double> SigmaA;
  double newLambda;
  double prefactor =  13.6*13.6/1000000;

  
  //SigmaA[]
  if(fabs(Alpha) > 0.0000000001){
    for(Long64_t i = 0; i < total; ++i){ //loop over muons
      L_Tree->GetEntry(i);
      for(std::map<int,double>::iterator iter = PtrLj->Array.begin(); iter != PtrLj->Array.end(); ++iter){
	newLambda = Lambda[iter->first] - Alpha*Gradient[iter->first];
	if(newLambda < 0.000003067) newLambda = 0.000003067;
	SigmaA[i] += iter->second * newLambda;
      }
      SigmaA[i] *= prefactor;
    }
  }else{
    for(std::map<int,double>::iterator iter = Sigma.begin(); iter != Sigma.end(); ++iter){
      SigmaA[iter->first] = iter->second;
    }
  }


  for(Long64_t i = 0; i < total; ++i){ //loop over muons
    if(Sigma[i] > 0){
      cost += S_array[i]/SigmaA[i] + TMath::Log(SigmaA[i]);
    }
  }

  return cost;
}
