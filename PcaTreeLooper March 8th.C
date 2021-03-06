
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


double prefactor = 13.6*13.6/1000;


//Define RootMap
ClassImp(RootMap);

RootMap::RootMap(){};
RootMap::~RootMap(){};





//Fill Matrix L with muons from FIRST to LAST, and make the S vector. Ni is the number of voxels in the i direction.
void PcaTreeLooper::SLFill(int first, int last, int Nx, int Ny, int Nz, char * Lijname){

  if(fChain==0) return;
  
  //Calculation Time
  clock_t LijStart = clock();
  
  double S;

  //char RootFileOut[FILENAME_MAX];
  //sprintf(RootFileOut, "L_%d_%d.root",first,last);
  TFile *fpout = new TFile(Lijname,"RECREATE");
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

  //FitQuality Cutoff
  double FitQual = 0.5;

  
  //Voxel dimensions
  double Vx = (MaxX-MinX)/Nx;
  double Vy = (MaxY-MinY)/Ny;
  double Vz = (MaxZ-MinZ)/Nz;

  //In the event of a muon "missing" it's PCA, ignore it. this is rare, and down to floating point errors.
  int PCAmiss;
  int MissCount = 0;



  Long64_t nb = 0;
  Long64_t nbytes = 0;
  
  map<int,double> L;
  int count = 0;
  //Loop over Muon Branches of trees in Chain
  for(int jentry = first; jentry < last; jentry++){
    Long64_t ientry = LoadTree((Long64_t)jentry);
    if(ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    
    //Reset L Array
    L.clear();

    //Reset PCAmiss flag
    PCAmiss = 0;

    //Only use muon if there is a PCA
    if(gotRecoPCA != 1) break;
    
    //Scattering Angle Cuttof --- Cuts out about 5% of muons for a value of 1.
    if(thetaTrue < 0.2){

      //Fit Quality Cuttof --- Cuts out additional ~5% of muons for a value of 0.5.
      if(xyzFitQualTrue[0] < FitQual && xyzFitQualTrue[1] < FitQual){

	//PCA in Box defined above?  --- Cuts out additional ~20% of muons.
	if(xPosTrue > MinX && xPosTrue < MaxX && yPosTrue < MaxY && yPosTrue > MinY && zPosTrue < MaxZ && zPosTrue > MinZ){
	
	  //Must Enter through top plane, and exit through bottom plane. --- Cuts out additional ~0.5% of muons.
	  double xEntry = xzGradTrue[0] * MaxZ + xzCutTrue[0];
	  double xExit = xzGradTrue[1] * MinZ + xzCutTrue[1];
	  double yEntry = yzGradTrue[0] * MaxZ + yzCutTrue[0];
	  double yExit = yzGradTrue[1] * MinZ + yzCutTrue[1];
	  if(xEntry > MinX && xEntry < MaxX && xExit > MinX && xExit < MaxX && yEntry > MinY && yEntry < MaxY && yExit > MinY && yExit < MaxY){

	    xzGradTrue[0] = (xEntry-xPosTrue)/(MaxZ-zPosTrue);
	    xzGradTrue[1] = (xExit-xPosTrue)/(MinZ-zPosTrue);
	    yzGradTrue[0] = (yEntry-yPosTrue)/(MaxZ-zPosTrue);
	    yzGradTrue[1] = (yExit-yPosTrue)/(MinZ-zPosTrue);

	    xzCutTrue[0] = xPosTrue-xzGradTrue[0]*zPosTrue;
	    xzCutTrue[1] = xPosTrue-xzGradTrue[1]*zPosTrue;
	    yzCutTrue[0] = yPosTrue-yzGradTrue[0]*zPosTrue;
	    yzCutTrue[1] = yPosTrue-yzGradTrue[1]*zPosTrue;


	    
	    //PCA Voxel ID
	    int PCAxID = (int) floor((xPosTrue - MinX) / Vx);
	    int PCAyID = (int) floor((yPosTrue - MinY) / Vy);
	    int PCAzID = (int) floor((zPosTrue - MinZ) / Vz);

	    int pcaID = PCAxID + Nx * PCAyID + Nx * Ny * PCAzID;
	    	      
	    double xCurrent, yCurrent, zCurrent, ZVoxelExitX, ZVoxelExitY, ZVoxelSideX, ZVoxelSideY, Zbottom, Xleave, Yleave;
	    int VoxelEnd;
	    map<int,int> xCase, yCase, xSign, ySign;
	    int xCurID, yCurID, zCurID, CurID;
	    int part = 0;

	    xCurrent = xEntry;
	    yCurrent = yEntry;
	    zCurrent = MaxZ;
	    VoxelEnd = pcaID;
	    
	    //If the Gradient is positive, then the voxel number decreases, if it is negative, it increases.
	    for(int a = 0; a < 2; a++){
	      xCase[a] = (xzGradTrue[a] > 0) ? 0 : 1;
	      yCase[a] = (yzGradTrue[a] > 0) ? 0 : 1;
	      xSign[a] = (xzGradTrue[a] > 0) ? -1 : 1;
	      ySign[a] = (yzGradTrue[a] > 0) ? -1 : 1;
	    }
	    
	    //Current Voxel ID and stuff
	    xCurID = (int) floor((xCurrent-MinX)/Vx);
	    yCurID = (int) floor((yCurrent-MinY)/Vy);
	    zCurID = (int) floor((zCurrent-MinZ)/Vz);
	    CurID = xCurID + Nx * yCurID + Nx*Ny * zCurID;
	    
	    //While loop until it hits the PCA. If it misses due to floating point errors, it is caught and the loop is exited and the second loop is never run.
	    while(part == 0){
	      
	      Zbottom = zCurID*Vz + MinZ;
	      Xleave = xzGradTrue[part] * Zbottom + xzCutTrue[part];
	      Yleave = yzGradTrue[part] * Zbottom + yzCutTrue[part];
	      
	      //Check if it leaves through the bottom plane of the voxel.
	      if(Xleave < (xCurID+1)*Vx + MinX && Xleave > xCurID*Vx + MinX && Yleave < (yCurID+1)*Vy + MinY && Yleave > yCurID*Vy + MinY){
		L[CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
		xCurrent = Xleave;
		yCurrent = Yleave;
		zCurrent = Zbottom;
		zCurID -= 1;
		CurID -= Nx*Ny;
		
	      }else{
		
		//Work out Exit coords for X Y
		ZVoxelSideX = (xCurID+xCase[part])*Vx + MinX;
		ZVoxelSideY = (yCurID+yCase[part])*Vy + MinY;
		
		//Work out Z coords for these exits in X&Y
		ZVoxelExitX = (ZVoxelSideX - xzCutTrue[part])/xzGradTrue[part];
		ZVoxelExitY = (ZVoxelSideY - yzCutTrue[part])/yzGradTrue[part];
		
		//which is higher?
		double ZExit = max(ZVoxelExitX, ZVoxelExitY);
		if(ZVoxelExitX > ZVoxelExitY){
		  xCurID += xSign[part];
		  CurID += xSign[part];
		}else{
		  yCurID += ySign[part];
		  CurID += Nx * ySign[part];
		}
		
		//New Current Coords
		Xleave = xzGradTrue[part] * ZExit + xzCutTrue[part];
		Yleave = yzGradTrue[part] * ZExit + yzCutTrue[part];
		L[CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
		xCurrent = Xleave;
		yCurrent = Yleave;
		zCurrent = ZExit;
		
	      }
	      
	      //Get out of loop
	      if(CurID == VoxelEnd) part = 1;
	      
	      //Get out of loop AND prevent next loop from being ru; floating point error somewhere.
	      if(CurID < 0 || zCurrent < MinZ || xCurrent < MinX || xCurrent > MaxX || yCurrent < MinY || yCurrent > MaxY) part = 2, PCAmiss = 1;
	      
	    }//Close While part=0 loop.
	    
	    //We should now be at the PCA, do this voxel.
	    
	    //Work out Extra PCA distance, will be added on later.
	    int PCAID = CurID;
	    double ExtraPCA = sqrt((xCurrent-xPosTrue)*(xCurrent-xPosTrue)+(yCurrent-yPosTrue)*(yCurrent-yPosTrue)+(zCurrent-zPosTrue)*(zCurrent-zPosTrue));
	    
	    xCurrent = xPosTrue;
	    yCurrent = yPosTrue;
	    zCurrent = zPosTrue;
	    
	    

	    VoxelEnd = (int) floor((xExit-MinX)/Vx) + Nx * (int) floor((yExit-MinY)/Vy);
	    int check = 0;
	    

	    while(part == 1){
	      check++;
	      
	      Zbottom = zCurID*Vz + MinZ;
	      Xleave = xzGradTrue[part] * Zbottom + xzCutTrue[part];
	      Yleave = yzGradTrue[part] * Zbottom + yzCutTrue[part];
	      

	      //if(jentry == 477) printf("%f \t %f \t %f \t %f \t %f \t %f %d \n",xCurrent,Xleave,yCurrent,Yleave,Zbottom, zCurrent,check);

	      //Check if it leaves through the bottom plane of the voxel.
	      if(Xleave < (xCurID+1)*Vx + MinX && Xleave > xCurID*Vx + MinX && Yleave < (yCurID+1)*Vy + MinY && Yleave > yCurID*Vy + MinY){
		L[CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
		xCurrent = Xleave;
		yCurrent = Yleave;
		zCurrent = Zbottom;
		zCurID -= 1;
		CurID -= Nx*Ny;
	      }else{
		
		//Work out Exit coords for X Y
		ZVoxelSideX = (xCurID+xCase[part])*Vx + MinX;
		ZVoxelSideY = (yCurID+yCase[part])*Vy + MinY;
		
		//Work out Z coords for these exits in X&Y
		ZVoxelExitX = (ZVoxelSideX - xzCutTrue[part])/xzGradTrue[part];
		ZVoxelExitY = (ZVoxelSideY - yzCutTrue[part])/yzGradTrue[part];

		
		//if(jentry == 477) printf("\t\t\t\t\t\t\t\t\t\t\t%f \t %f \t %f \t %f \t  %f \t %d \n",ZVoxelSideX,ZVoxelSideY,zCurrent,ZVoxelExitX,ZVoxelExitY,check);
		
		//which is higher?
		double ZExit = max(ZVoxelExitX, ZVoxelExitY);
		
		//New Current Coords
		Xleave = xzGradTrue[part] * ZExit + xzCutTrue[part];
		Yleave = yzGradTrue[part] * ZExit + yzCutTrue[part];
		L[CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+(ZExit-zCurrent)*(ZExit-zCurrent),0.5);
		xCurrent = Xleave;
		yCurrent = Yleave;
		zCurrent = ZExit;
		if(ZVoxelExitX > ZVoxelExitY){
		  xCurID += xSign[part];
		  CurID += xSign[part];
		}else{
		  yCurID += ySign[part];
		  CurID += Nx * ySign[part];
		}

	      }
	     
	      //Get out of loop
	      if(CurID == VoxelEnd) part = 2;
	      
	      //Get out of loop AND prevent next loop from being run, as something has gone wrong, provided it hasn't just hit the last voxel.
	      if(part == 1 && (CurID < 0 || zCurrent < MinZ || xCurrent < MinX || xCurrent > MaxX || yCurrent < MinY || yCurrent > MaxY)) part = 2, PCAmiss = 1;
	      
	    }//Close While part=1 loop.
	    
	    
	    count++;
	    

	    //only fill if the muon hit it's PCA -- cuts out about 3% more muons.
	    if(PCAmiss == 0){

	      //Add in the extra bit of length for the PCA that's not accounted for by the loop.
	      L[PCAID] += ExtraPCA;
	
	      //Fill SL Root File.
	      PtrLj->Array.clear();
	      for(std::map<int,double>::iterator iter = L.begin(); iter != L.end(); ++iter){
		PtrLj->Array[iter->first] = iter->second;
	      }
	      //Fill Scattering Angle Array, so that the root file needn't be needlessly looped through when S is already in memory.
	      S_array[jentry] = S = thetaTrue*thetaTrue;

	      L_Tree->Fill();
	      

	    }else{//Close if for "if missed PCA don't fill L Matrix". else to count misses.
	      MissCount ++;
	    }

	  

	  }//Close check that the muon enters and exits through the scintillator planes.

	}//Close "Is PCA in box" If.
	
      }//Close Fit Quality Cuttof If.

    }//Close Scattering Angle Cuttoff If.
    
  }//Close loop over muons
  cout << count << "\t" << MissCount << endl;


  //Save Root File, then Close it.
  L_Tree->AutoSave();
  fpout->Close();
  
  printf("Lij & S Fill: %f\n", ((double)(clock() - LijStart) / CLOCKS_PER_SEC));
    
}//End LSFill





//First instance of Lambda being Filled.
void PcaTreeLooper::LambdaFill(int VoxelCount, char* Lambdaname){
  clock_t LambdaStart = clock();

  //Save Lambda as a root file so it can be looked at later.
  //char RootFileOut[FILENAME_MAX];
  //sprintf(RootFileOut, "Lambda.root");
  TFile *fpout = new TFile(Lambdaname,"RECREATE");
  TTree *Lambda_Tree = new TTree("Lambda","Lambda");
  RootMap *Lambda_in_tree = new RootMap();
  Lambda_Tree->Branch("Voxels","RootMap",&Lambda_in_tree);
  
  for(int i=0; i < VoxelCount; i++){
    Lambda[i] = 0.000003291;
    Lambda_in_tree->Array[i] = 0.000003291;
  }
  Lambda_Tree->Fill();
  Lambda_Tree->AutoSave();
  printf("Lambda Fill: %f\n", ((double)(clock() - LambdaStart) / CLOCKS_PER_SEC));
  fpout->Close();
}



//Update Lambda by subtracting Alpha*Grad(Psi) from it.
void PcaTreeLooper::LambdaAlpha(double Alpha, char * Lambdaname){
  clock_t LambdaStart = clock();

  //Save Lambda as a root file so it can be looked at later.
  //char RootFileOut[FILENAME_MAX];
  //sprintf(RootFileOut, "Lambda.root");
  TFile *fpout = new TFile(Lambdaname,"RECREATE");
  TTree *Lambda_Tree = new TTree("Lambda","Lambda");
  RootMap *Lambda_in_tree = new RootMap();
  Lambda_Tree->Branch("Voxels","RootMap",&Lambda_in_tree);

  double newLambda;

  for(std::map<int,double>::iterator iter = Lambda.begin(); iter != Lambda.end(); ++iter){
    newLambda = iter->second - Alpha*Gradient[iter->first];
    //      if(Gradient[iter->first] > 10000000) printf("%d =>  %f\n",iter->first,Gradient[iter->first]);
    if(newLambda < 0.000003067) newLambda = 0.000003067;
    Lambda[iter->first] = newLambda;
    Lambda_in_tree->Array[iter->first] = newLambda;
  }
  printf("Lambda Alpha Fill: %f\n", ((double)(clock() - LambdaStart) / CLOCKS_PER_SEC));
  Lambda_Tree->Fill();
  Lambda_Tree->AutoSave();
  fpout->Close();
}



void PcaTreeLooper::SigmaFill(int first, int last, char* Lijname, char* Lambdaname, char* SigSname){
  
  clock_t SigmaStart = clock();

  //Open SL Root file.
  //char RootFileOut[FILENAME_MAX];
  //sprintf(RootFileOut, "L_%d_%d.root",first,last);
  TFile *fpinLij = new TFile(Lijname,"OLD");
  TTree *L_Tree = (TTree*) fpinLij->Get("L");
  RootMap *PtrLj = 0;
  L_Tree->SetBranchAddress("Voxels",&PtrLj);

  //Open Lambda File
  //Open Gradient Root File
  TFile *fpinLambda = new TFile(Lambdaname,"OLD");
  TTree *Lambda_Tree = (TTree*) fpinLambda->Get("Lambda");
  RootMap *PtrLambda = 0;
  Lambda_Tree->SetBranchAddress("Voxels",&PtrLambda);
  Lambda_Tree->GetEntry(0);


  double newLambda;

  double S;
  double SigmaForTree;
  //sprintf(RootFileOut, "SigS_%d_%d.root",first,last);
  TFile *fpoutSigS = new TFile(SigSname,"RECREATE");
  TTree *SigS_Tree = new TTree("SS","Sigma and S");

  SigS_Tree->Branch("SigmaSquared",&SigmaForTree,"SigmaSquared/D");
  SigS_Tree->Branch("ThetaSquared",&S,"ThetaSquared/D");


  Long64_t total = L_Tree->GetEntries();
  for(Long64_t i = 0; i < total; ++i){ //loop over muons
    L_Tree->GetEntry(i);
    Sigma[i] = 0;
    for(std::map<int,double>::iterator iter = PtrLj->Array.begin(); iter != PtrLj->Array.end(); ++iter){
      Sigma[i] += iter->second * PtrLambda->Array[iter->first];
    }
    S = S_array[i];
    Sigma[i] *= prefactor;
    SigmaForTree = Sigma[i];
    SigS_Tree->Fill();
  }
  
  SigS_Tree->AutoSave();

  printf("Sigma Fill: %f\n", ((double)(clock() - SigmaStart) / CLOCKS_PER_SEC));
  fpinLij->Close();
  fpinLambda->Close();
  fpoutSigS->Close();


}





//Gradient Vector.
void PcaTreeLooper::GradientFill(int first, int last, char* Lijname, char* Gname){
  clock_t GradientStart = clock();
  
  //open LS root file.
  //char RootFileOut[FILENAME_MAX];
  //sprintf(RootFileOut, "L_%d_%d.root",first,last);
  TFile *fpin = new TFile(Lijname,"OLD");
  TTree *L_Tree = (TTree*) fpin->Get("L");
  RootMap *PtrLj = 0;
  L_Tree->SetBranchAddress("Voxels",&PtrLj);
  double S;
  L_Tree->SetBranchAddress("ThetaSquared",&S);

  //Save gradient as a root file so it can be looked at later.
  //sprintf(RootFileOut,"G.root");
  TFile *fpout = new TFile(Gname,"RECREATE");
  TTree *G_Tree = new TTree("G","Gradient");
  RootMap *Grad = new RootMap();
  G_Tree->Branch("Voxels","RootMap",&Grad);
  
  

  Gradient.clear();
  Grad->Array.clear();

  //Fill the array with zeros before filling it with other numbers. This is to prevent calls to undefeined elements.
  for(int i = 0; i < 999999; i++){
    Grad->Array[i] = 0;
    Gradient[i] = 0;
  }

  Long64_t total = L_Tree->GetEntries();
  for(Long64_t i = 0; i < total; ++i){
    L_Tree->GetEntry(i);
    for(std::map<int,double>::iterator iter = PtrLj->Array.begin(); iter != PtrLj->Array.end(); ++iter){
      Grad->Array[iter->first] += iter->second * ((1 - (S_array[(int)i]/Sigma[(int)i])))/Sigma[i];
      Gradient[iter->first] += iter->second * ((1 - (S_array[(int)i]/Sigma[(int)i])))/Sigma[i];
    }
    
    
  }

  G_Tree->Fill();
  
  G_Tree->AutoSave();
  fpin->Close();
  fpout->Close();
  printf("Gradient Fill: %f\n", ((double)(clock() - GradientStart) / CLOCKS_PER_SEC));
}






//cost Function
double PcaTreeLooper::Cost(double Alpha){
  double cost=0;
  char RootFileOut[FILENAME_MAX];
  sprintf(RootFileOut, "L_0_100000.root");
  TFile *fpin = new TFile(RootFileOut,"OLD");
  TTree *L_Tree = (TTree*) fpin->Get("L");
  RootMap *PtrLj = 0;
  L_Tree->SetBranchAddress("Voxels",&PtrLj);
  Long64_t total = L_Tree->GetEntries();

  map<int,double> SigmaA;
  double newLambda;

  
  //SigmaA[]
  if(fabs(Alpha) > 0.0000000001){
    for(Long64_t i = 0; i < total; ++i){ //loop over muons
      L_Tree->GetEntry(i);
      SigmaA[i] = 0;
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
    cost += S_array[i]/SigmaA[i] + TMath::Log(SigmaA[i]);
  }

  return cost;

  fpin->Close();
}




//Draw Gradient Slices
void PcaTreeLooper::DrawGradientSlices(int topSlice, int bottomSlice, int Nx, int Ny, int Nz){
  if(topSlice > bottomSlice && topSlice <= Nz && bottomSlice >= 0){
    //Remember - 0 is the bottom of the box
    //Open Gradient Root File
    char GradFile[80];
    sprintf(GradFile,"Gradient.root");
    TFile *fpin = new TFile(GradFile,"OLD");
    TTree *G_Tree = (TTree*) fpin->Get("G");
    RootMap *PtrG = 0;
    G_Tree->SetBranchAddress("Voxels",&PtrG);
    G_Tree->GetEntry(0);

    cout << "1" << endl;

    //Open Lambda File
    //Open Gradient Root File
    char LambdaFile[80];
    sprintf(LambdaFile,"Lambda.root");
    TFile *fpinLambda = new TFile(LambdaFile,"OLD");
    TTree *Lambda_Tree = (TTree*) fpinLambda->Get("Lambda");
    RootMap *PtrLambda = 0;
    Lambda_Tree->SetBranchAddress("Voxels",&PtrLambda);
    Lambda_Tree->GetEntry(0);

    cout << "2" << endl;
    

    TFile *fpOu1 = new TFile("GradientSlices.root","RECREATE");
    fpOu1->Write();
    TCanvas *canProj = new TCanvas();
    canProj->Divide(3,3);
    cout << "3" << endl;

    TFile *fpOu2 = new TFile("LambdaSlices.root","RECREATE");
    fpOu2->Write();
    TCanvas *canProj2 = new TCanvas();
    canProj2->Divide(3,3);
    cout << "4" << endl;
    
    int z = 50;
    char HistName[80];
    for(int SliceNo = 1; SliceNo < 10; SliceNo++){
      sprintf(HistName,"Slice_%d",SliceNo);
      TH2F *CurrentSlice = new TH2F(HistName,HistName, Nx, 0.0, (double) Nx, Ny, 0.0, (double) Ny);
      TH2F *CurrentLSlice = new TH2F(HistName,HistName, Nx, 0.0, (double) Nx, Ny, 0.0, (double) Ny);

      for(int x = 1; x <= 100; x++){
	for(int y = 1; y <= 100; y++){
	  CurrentSlice->SetBinContent(x,y,PtrG->Array[(z+SliceNo)*Nx*Ny + y*Nx + x]);
	  CurrentLSlice->SetBinContent(x,y,PtrLambda->Array[(z+SliceNo)*Nx*Ny + y*Nx + x]);
	}
      }


    canProj->cd(SliceNo);
    CurrentSlice->Draw("colz");
    CurrentSlice->SetStats(0);
    CurrentSlice->GetXaxis()->SetTitle("Voxel ID");

    canProj2->cd(SliceNo);
    CurrentLSlice->Draw("colz");
    CurrentLSlice->SetStats(0);
    CurrentLSlice->GetXaxis()->SetTitle("Voxel ID");
      
    }
      
      
    
    fpinLambda->Close();
    fpin->Close();
    
  }
}
