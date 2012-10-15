#define LambdaPcaTreeLooperMLSD_cxx
#include "LambdaPcaTreeLooperMLSD.h"
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

#define MAKE_DEBUG_TREE 1

//#define MAKE_GRAD_TREE 1

#define LAMBDA_AIR (1./3.039e5)

//Matrix multiplication
void LambdaPcaTreeLooperMLSD::mm_mul (double A[2][2], double B[2][2], double C[2][2]){ 
 
          int i, j, k;
           double sum;
   
           for (i = 0; i < 2; i++) {
                   for (j = 0; j < 2; j++) {
                           sum = 0;
                           for (k = 0; k < 2; k++) {
                                   sum += A[i][k] * B[k][j];
                           }
                           C[i][j] = sum;
                   }
           }
  }

// Evaluate the product DAD
double LambdaPcaTreeLooperMLSD::DAD (double D[2], double A[2][2]){

	double dad = 0.;
	for(int i=0; i<2; i++){
		for(int j=0; j<2; j++){

		dad+= A[i][j]*D[i]*D[j];
		}
	}

	return dad;

}

void LambdaPcaTreeLooperMLSD::Loop()
{
   //   In a ROOT session, you can do:
   //      Root > .L LambdaPcaTreeLooperMLSD.C
   //      Root > LambdaPcaTreeLooperMLSD t
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

void LambdaPcaTreeLooperMLSD::FillPosHist(TH3F *histPos, Double_t thetaCut)
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


//Fill Matrix L with muons from FIRST to LAST, and make the dx and dtheta vectors. Ni is the number of voxels in the i direction.
void LambdaPcaTreeLooperMLSD::SLFill(int first, int last, int Nx, int Ny, int Nz){
   //At some point might update this so that one can specify the first and last and only have array that size
   fNumMuons=last;

   S = new double[fNumMuons];
   dx = new double[fNumMuons];
   dy = new double[fNumMuons];
   dThetax = new double[fNumMuons];
   dThetay = new double[fNumMuons];
   PreFactorEng = new double[fNumMuons];
   memset(S,0,fNumMuons*sizeof(double));
   memset(PreFactorEng,0,fNumMuons*sizeof(double));
   memset(dx,0,fNumMuons*sizeof(double));
   memset(dy,0,fNumMuons*sizeof(double));
   memset(dThetax,0,fNumMuons*sizeof(double));
   memset(dThetay,0,fNumMuons*sizeof(double));

#ifdef MAKE_DEBUG_TREE
   Int_t muonNum,voxelNum;
   Double_t lij;
   Double_t tij;
Double_t aij;
   Double_t bij;
   Double_t cij;
   Double_t dxNum;
   Double_t dxthetaNum;
   Double_t dthetaNum;
   TFile *fpLij = new TFile("/unix/creamtea/sfayer/temp/lij.root","RECREATE");
   TTree *lijTree = new TTree("lijTree","lijTree");
   lijTree->Branch("muonNum",&muonNum,"muonNum/I");
   lijTree->Branch("voxelNum",&voxelNum,"voxelNum/I");
   lijTree->Branch("lij",&lij,"lij/D");
   lijTree->Branch("tij",&tij,"tij/D");
lijTree->Branch("aij",&aij,"aij/D");
   lijTree->Branch("bij",&bij,"bij/D");
   lijTree->Branch("cij",&cij,"cij/D");
   lijTree->Branch("dxNum",&dxNum,"dxNum/D");
   lijTree->Branch("dxthetaNum",&dxthetaNum,"dxthetaNum/D");
   lijTree->Branch("dthetaNum",&dthetaNum,"dthetaNum/D");
#endif


   if(fChain==0) return;
  
   //Calculation Time
   clock_t LijStart = clock();

   //Size of cuboid
   double MaxX = 500;
   double MinX = -500;
   double MaxY = 500;
   double MinY = -500;
   double MaxZ = 500;
   double MinZ = -500;

//   double deltax = 0.;
//   double deltay = 0.;
//   double Lxy = 0.;

   //FitQuality Cutoff
   double FitQual = 0.5;

   //Voxel dimensions
   double Vx = (MaxX-MinX)/Nx;
   double Vy = (MaxY-MinY)/Ny;
   double Vz = (MaxZ-MinZ)/Nz;

   //In the event of a muon "missing" it's PCA, ignore it. this is rare, and down to doubleing point errors.
   int PCAmiss;
   int MissCount = 0;

   Long64_t nb = 0;
   Long64_t nbytes = 0;

   int count = 0;
   //Loop over Muon Branches of trees in Chain
   for(int jentry = first; jentry < last; jentry++){
      std::map<int, double> tempMap;
      std::map<int, double> tempMapT;
      Long64_t ientry = LoadTree((Long64_t)jentry);
      if(ientry < 0) break;
      nb = fChain->GetEntry(jentry); nbytes += nb;

      //Reset PCAmiss flag
      PCAmiss = 0;

      //Only use muon if there is a PCA
      if(gotRecoPCA != 1) continue;
    
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
		  if(pcaID<0 || pcaID>=Nx*Ny*Nz) {
		     std::cout << pcaID <<"\n";
		  }


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
		  
		  //Hack because z=500 is outside the box
		  if(zCurID==Nz) zCurID--;
		  if(yCurID==Ny) yCurID--;
		  if(xCurID==Nx) xCurID--;

		  CurID = xCurID + Nx * yCurID + Nx*Ny * zCurID;

		  if(CurID<0 || CurID>=Nx*Ny*Nz) {
		     std::cout << "CurID " << CurID << "\t" << xCurrent << "\t" << yCurrent << "\t" << zCurrent << "\t" << zCurID <<  "\n";
		  }
	    
		  //While loop until it hits the PCA. If it misses due to doubleing point errors, it is caught and the loop is exited and the second loop is never run.
		  while(part == 0){
		     Zbottom = zCurID*Vz + MinZ;
		     Xleave = xzGradTrue[part] * Zbottom + xzCutTrue[part];
		     Yleave = yzGradTrue[part] * Zbottom + yzCutTrue[part];
	      
		     //Check if it leaves through the bottom plane of the voxel.
		     if(Xleave < (xCurID+1)*Vx + MinX && Xleave > xCurID*Vx + MinX && Yleave < (yCurID+1)*Vy + MinY && Yleave > yCurID*Vy + MinY){
			//			L[jentry][CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
			tempMap[CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
			tempMapT[CurID] = Zbottom + MaxZ;

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
			//			L[jentry][CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
			tempMap[CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
			tempMapT[CurID] = MaxZ - Zbottom ;
			xCurrent = Xleave;
			yCurrent = Yleave;
			zCurrent = ZExit;
		
		     }
	      
		     //Get out of loop
		     if(CurID == VoxelEnd) part = 1;
      
		     //Get out of loop AND prevent next loop from being run; doubleing point error somewhere.
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
	      
		     //Check if it leaves through the bottom plane of the voxel.
		     if(Xleave < (xCurID+1)*Vx + MinX && Xleave > xCurID*Vx + MinX && Yleave < (yCurID+1)*Vy + MinY && Yleave > yCurID*Vy + MinY){
			//L[jentry][CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
			tempMap[CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+Vz*Vz,0.5);
			tempMapT[CurID] = Zbottom + MaxZ;
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
		
			//New Current Coords
			Xleave = xzGradTrue[part] * ZExit + xzCutTrue[part];
			Yleave = yzGradTrue[part] * ZExit + yzCutTrue[part];
			//			L[jentry][CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+(ZExit-zCurrent)*(ZExit-zCurrent),0.5);
			tempMap[CurID] = pow((xCurrent - Xleave)*(xCurrent-Xleave)+(yCurrent - Yleave)*(yCurrent-Yleave)+(ZExit-zCurrent)*(ZExit-zCurrent),0.5);
			tempMapT[CurID] = MaxZ - Zbottom;
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
		     tempMap[PCAID] += ExtraPCA;
	
		     //Fill Scattering Angle Array, so that the root file needn't be needlessly looped through when S is already in memory.
		     S[jentry] = thetaTrue*thetaTrue;
		     dThetax[jentry] = thetaxzTrue;
		     dThetay[jentry] = thetayzTrue;
		     dx[jentry] = dxTrue; //in mm
		     dy[jentry] = dyTrue; //in mm

		     double Lxy = TMath::Sqrt(1. + (xzGradTrue[0]*xzGradTrue[0]) + (yzGradTrue[0]*yzGradTrue[0]));

			if (intEng<2000.) PreFactorEng[jentry]=13.6*13.6/(1000.*1000.);
		     else PreFactorEng[jentry]=13.6*13.6/(4000.*4000.);
		     //Now copy tempMap to the big L map of maps
		     //		     std::cout << tempMap.size() << "\n";
#ifdef MAKE_DEBUG_TREE
		     muonNum=jentry;
#endif
		     
		     for(std::map<int, double>::iterator tempIt=tempMap.begin();
			 tempIt!=tempMap.end();tempIt++) {
			L[jentry][tempIt->first]=tempMap[tempIt->first]*Lxy;//tempIt->second;
			T[jentry][tempIt->first]= tempMapT[tempIt->first]*Lxy;
			if (T[jentry][tempIt->first]>((MaxZ-MinZ)*TMath::Sqrt(3.))) T[jentry][tempIt->first]=(MaxZ-MinZ)*TMath::Sqrt(3.);//assumes a cube
			A[jentry][tempIt->first]=L[jentry][tempIt->first];
			B[jentry][tempIt->first]=(0.5*L[jentry][tempIt->first]*L[jentry][tempIt->first])+(L[jentry][tempIt->first]*T[jentry][tempIt->first]);
			C[jentry][tempIt->first]=(L[jentry][tempIt->first]*L[jentry][tempIt->first]*L[jentry][tempIt->first]/3.)+(L[jentry][tempIt->first]*L[jentry][tempIt->first]*T[jentry][tempIt->first])+(L[jentry][tempIt->first]*T[jentry][tempIt->first]*T[jentry][tempIt->first]);
#ifdef MAKE_DEBUG_TREE
			lij=tempIt->second;
			voxelNum=tempIt->first;
			tij = T[jentry][tempIt->first];
			aij = A[jentry][tempIt->first];
			bij = B[jentry][tempIt->first];
			cij = C[jentry][tempIt->first];
			lijTree->Fill();
#endif
		     }
		     

#ifdef MAKE_DEBUG_TREE
		     dxNum = sqrt((dx[jentry]*dx[jentry]+dy[jentry]*dy[jentry])/2);
		     dxthetaNum = (dx[jentry]*dThetax[jentry]+dy[jentry]*dThetay[jentry])/2;
		     dthetaNum = sqrt((dThetax[jentry]*dThetax[jentry]+dThetay[jentry]*dThetay[jentry])/2);
		     lijTree->Fill();
#endif

		     if(jentry == 146891) {
			std::cout << jentry << "\t" << S[jentry] << "\t" << L[jentry][PCAID] << "\n";
		     }

		  }else{//Close if for "if missed PCA don't fill L Matrix". else to count misses.
		     MissCount ++;
		  }

	       }//Close check that the muon enters and exits through the scintillator planes.

	    }//Close "Is PCA in box" If.

	 }//Close Fit Quality Cuttof If.

      }//Close Scattering Angle Cuttoff If.
    
   }//Close loop over muons
   cout << count << "\t" << MissCount << "\t" << L.size() <<  endl;
#ifdef MAKE_DEBUG_TREE
   lijTree->AutoSave();
   fpLij->Close();
#endif

   printf("Lij & S Fill: %f\n", ((double)(clock() - LijStart) / CLOCKS_PER_SEC));
        
}//End LSFill


//First instance of Lambda being Filled.
void LambdaPcaTreeLooperMLSD::LambdaFill(int VoxelCount){
   fVoxelCount=VoxelCount;
   clock_t LambdaStart = clock();

   Lambda = new double [fVoxelCount];


   for(int i=0; i < fVoxelCount; i++){
      //      Lambda[i] = 0.000003291;
      Lambda[i]=LAMBDA_AIR;
   }
   std::cout << "Lambda[0] " << Lambda[0] << "\n";
   printf("Lambda Fill: %f\n", ((double)(clock() - LambdaStart) / CLOCKS_PER_SEC));

}


//Update Lambda by subtracting Alpha*Grad(Psi) from it.
void LambdaPcaTreeLooperMLSD::LambdaAlpha(double Alpha){
   clock_t LambdaStart = clock();

   double newLambda;

   for(int voxId=0;voxId<fVoxelCount;voxId++) {
      newLambda = Lambda[voxId] - Alpha*Gradient[voxId];
      if(newLambda < LAMBDA_AIR) newLambda = LAMBDA_AIR;
      Lambda[voxId] = newLambda;
   }
   printf("Lambda Alpha Fill: %f\n", ((double)(clock() - LambdaStart) / CLOCKS_PER_SEC));
}
//Update Lambda with MLSD formula
void LambdaPcaTreeLooperMLSD::LambdaNew(){
   clock_t LambdaStart = clock();

   double newLambda;

   for(int voxId=0;voxId<fVoxelCount;voxId++) {
      newLambda = Sj[voxId]/(2.*M[voxId]);
      if(newLambda < LAMBDA_AIR) newLambda = LAMBDA_AIR;
      if(isnan(newLambda)) newLambda = LAMBDA_AIR;
      if(isinf(newLambda)) newLambda = LAMBDA_AIR;
      if(M[voxId]<4) newLambda = LAMBDA_AIR;
      /*if(newLambda>10) {
std::cout << "voxID: " << voxId << "    Lambda: " << newLambda << "    S: " << Sj[voxId] << "     M: " << M[voxId] << endl;
	//newLambda=LAMBDA_AIR;
}*/

Lambda[voxId] = newLambda;
   }
   printf("Lambda Alpha Fill: %f\n", ((double)(clock() - LambdaStart) / CLOCKS_PER_SEC));
}

//This function fills D^{-1} matrix
void LambdaPcaTreeLooperMLSD::SigmaFill(){
  
   clock_t SigmaStart = clock();
#ifdef MAKE_DEBUG_TREE
   static TFile *fpSigma = new TFile("/unix/creamtea/sfayer/temp/sigma.root","RECREATE");
   static TTree *sigmaTree = (TTree*) fpSigma->Get("sigmaTree");
   static Double_t sigmax=0;
   static Double_t sigmaxtheta=0;
   static Double_t sigmatheta=0;
#endif
   static int doneInit=0;
   static Int_t muonNum=0;
   static Double_t determinant=0;
   static Double_t prefactorNum=0;
   static Double_t aDNum=0;
   static Double_t bDNum=0;
   static Double_t cDNum=0;
   static Double_t thetaNum=0;
   static Double_t dxNum=0;
   static Double_t LambdaNum=0;
   if(!doneInit) {
      aD = new double[fNumMuons];
      bD = new double[fNumMuons];
      cD = new double[fNumMuons];
      DET = new double[fNumMuons];
#ifdef MAKE_DEBUG_TREE
      if(!sigmaTree) {
	 sigmaTree = new TTree("sigmaTree","sigmaTree");
	 sigmaTree->Branch("muonNum",&muonNum,"muonNum/I");
	 sigmaTree->Branch("loopNum",&doneInit,"loopNum/I");
	 sigmaTree->Branch("prefactorNum",&prefactorNum,"prefactorNum/D");
	 sigmaTree->Branch("aDNum",&aDNum, "aDNum/D");
	 sigmaTree->Branch("bDNum",&bDNum, "bDNum/D");
	 sigmaTree->Branch("cDNum",&cDNum, "cDNum/D");
	 sigmaTree->Branch("sigmax",&sigmax,"sigmax/D");
	 sigmaTree->Branch("sigmaxtheta",&sigmaxtheta,"sigmaxtheta/D");
	 sigmaTree->Branch("sigmatheta",&sigmatheta,"sigmatheta/D");
	 sigmaTree->Branch("thetaNum",&thetaNum, "thetaNum/D");
	 sigmaTree->Branch("dxNum",&dxNum, "dxNum/D");
	 sigmaTree->Branch("determinant",&determinant, "determinant/D");
      }
      else {
	 sigmaTree->SetBranchAddress("muonNum",&muonNum);
	 sigmaTree->SetBranchAddress("loopNum",&doneInit);
	 sigmaTree->SetBranchAddress("aDNum",&aDNum);
	 sigmaTree->SetBranchAddress("bDNum",&bDNum);
	 sigmaTree->SetBranchAddress("cDNum",&cDNum);
	 sigmaTree->SetBranchAddress("sigmax",&sigmax);
	 sigmaTree->SetBranchAddress("sigmaxtheta",&sigmaxtheta);
	 sigmaTree->SetBranchAddress("sigmatheta",&sigmatheta);
	 sigmaTree->SetBranchAddress("thetaNum",&thetaNum);
	 sigmaTree->SetBranchAddress("dxNum",&dxNum);
	 sigmaTree->SetBranchAddress("determinant",&determinant);
      }
#endif
   }
   doneInit++;
   memset(aD,0,fNumMuons*sizeof(double));
   memset(bD,0,fNumMuons*sizeof(double));
   memset(cD,0,fNumMuons*sizeof(double));

   for(std::map<int,map<int,double> >::iterator iter1 = L.begin(); iter1 != L.end(); ++iter1){ //loop over muons
      aD[iter1->first] = 0.;
      bD[iter1->first] = 0.;
      cD[iter1->first] = 0.;
      determinant = 0.;

      for(std::map<int,double>::iterator iter2 = (iter1->second).begin(); iter2 != (iter1->second).end(); ++iter2){ //loop over voxels
	if (iter2->second==0.) {
	  std::cout << "iter2->second==0" << endl;
	  continue;
      }
	 cD[iter1->first] += A[iter1->first][iter2->first]*Lambda[iter2->first];
	 bD[iter1->first] += B[iter1->first][iter2->first]*Lambda[iter2->first];
	 aD[iter1->first] += C[iter1->first][iter2->first]*Lambda[iter2->first];
      }//end loop over voxels
      if(iter1->first==0) {
	 std::cout << "Before factor: " << cD[iter1->first] << "\n";
      }
	cD[iter1->first] *= PreFactorEng[iter1->first];
	bD[iter1->first] *= PreFactorEng[iter1->first];
	aD[iter1->first] *= PreFactorEng[iter1->first];
#ifdef MAKE_DEBUG_TREE
	sigmax = sqrt(aD[iter1->first]);
	sigmaxtheta = bD[iter1->first];
	sigmatheta = sqrt(cD[iter1->first]);
#endif

      cD[iter1->first] += 0.000049;
//      bD[iter1->first] -= 0.000049;
//      cD[iter1->first] += 0.000049;

	determinant = (aD[iter1->first]*cD[iter1->first]) - TMath::Power(bD[iter1->first], 2.);
	DET[iter1->first]=determinant;	

	aD[iter1->first]/=determinant;
	bD[iter1->first]/=(determinant*(-1.));
	cD[iter1->first]/=determinant;

      if(iter1->first==0) {
	 std::cout << "After factor: " << cD[iter1->first] << "\n";
     }

#ifdef MAKE_DEBUG_TREE
      prefactorNum=PreFactorEng[iter1->first];
      muonNum=iter1->first;   
	aDNum = aD[iter1->first];
	bDNum = bD[iter1->first];
	cDNum = cD[iter1->first];
	thetaNum = dThetax[iter1->first];
	dxNum = dx[iter1->first];
      sigmaTree->Fill();
#endif
   } //end loop over muons
#ifdef MAKE_DEBUG_TREE
   sigmaTree->AutoSave();
#endif

   printf("Sigma Fill: %f\n", ((double)(clock() - SigmaStart) / CLOCKS_PER_SEC));
}


//Gradient Vector.
void LambdaPcaTreeLooperMLSD::GradientFill(){
   clock_t GradientStart = clock();
#ifdef MAKE_GRAD_TREE
   static TFile *fpGrad = new TFile("/unix/creamtea/sfayer/temp/grad.root","RECREATE");
   static TTree *gradTree = (TTree*) fpGrad->Get("gradTree");
#endif
   static Int_t voxelNum=0;
   static Double_t lambdaNum=0;
   static Double_t gradientNum=0;
   static int madeGradient=0;
   static double mat_W[2][2];
   static double mat_D[2][2]; //SigmaD
   static double GRADmat_D[2][2]; // GRAD SigmaD
   static double mat_R[2][2]; // SigmaD*W
   static double mat_R2[2][2]; // SigmaD*W
   static double mat_Tot[2][2]; // SigmaD*W*SigmaD
   static double mat_Tot2[2][2]; // SigmaD*W*SigmaD
   static double vector[2];
   static double determinant;
   static double GRADdeterminant;
   static double GRADa;
   static double GRADb;
   static double GRADc;
   static double GRADtrace;
   static double GRADproduct;
   static double GRADproductY;
   static double trace;
   static double product;
   static double productY;

static double Lij;
static double Bij;
static double Cij;
static double ai;
static double bi;
static double ci;


   if(!madeGradient) {
      Gradient = new double [fVoxelCount];
#ifdef MAKE_GRAD_TREE
      if(!gradTree) {
	 gradTree= new TTree("gradTree","Tree of voxel gradients");
	 gradTree->Branch("voxelNum",&voxelNum,"voxelNum/I");
	 gradTree->Branch("loopNum",&madeGradient,"loopNum/I");
	 gradTree->Branch("lambdaNum",&lambdaNum,"lambdaNum/D");
	 gradTree->Branch("gradientNum",&gradientNum,"gradientNum/D");
	 gradTree->Branch("lambdaNum",&lambdaNum,"lambdaNum/D");
      }      
#endif
   }
   madeGradient++;
   memset(Gradient,0,fVoxelCount*sizeof(double)); 
   //  Gradient.clear();

   for(std::map<int,map<int,double> >::iterator iter1 = L.begin(); iter1 != L.end(); ++iter1){ 

	determinant = DET[iter1->first];
	 mat_D[0][0] = aD[iter1->first];
	 mat_D[0][1] = bD[iter1->first];
	 mat_D[1][0] = bD[iter1->first];
	 mat_D[1][1] = cD[iter1->first];

	 ai = aD[iter1->first];
	 bi = bD[iter1->first];
	 ci = cD[iter1->first];
       for(std::map<int,double>::iterator iter = (iter1->second).begin(); iter != (iter1->second).end(); ++iter){

	Lij = L[iter1->first][iter->first];
	Bij = B[iter1->first][iter->first];
	Cij = C[iter1->first][iter->first];

	GRADdeterminant = PreFactorEng[iter1->first]*determinant*((Lij*ai) + (Cij*ci) + (2*Bij*bi));
	GRADa = ((PreFactorEng[iter1->first]*Cij) - (GRADdeterminant*ai))/determinant;
	GRADb = -((PreFactorEng[iter1->first]*Bij) - (GRADdeterminant*bi))/determinant;
	GRADc = ((PreFactorEng[iter1->first]*Lij) - (GRADdeterminant*ci))/determinant;
	GRADtrace = (GRADa*Lij) + (2*GRADb*Bij) + (GRADc*Cij);


	GRADproduct = 0.5*( (2*ai*GRADa*dThetax[iter1->first]*dThetax[iter1->first]*Lij)  +  (2*(GRADa*bi + GRADb*ai)*(dThetax[iter1->first]*dThetax[iter1->first]*Bij  +  dThetax[iter1->first]*dx[iter1->first]*Lij)) + (2*bi*GRADb*(2*dThetax[iter1->first]*dx[iter1->first]*Bij + dx[iter1->first]*dx[iter1->first]*Lij + dThetax[iter1->first]*dThetax[iter1->first]*Cij )) + (2*ci*GRADc*dx[iter1->first]*dx[iter1->first]*Cij) + (2*(GRADb*ci + GRADc*bi)*(dx[iter1->first]*dx[iter1->first]*Bij  +  dThetax[iter1->first]*dx[iter1->first]*Cij)) + (2*dThetax[iter1->first]*dx[iter1->first]*Bij*(GRADa*ci + GRADc*ai)) );

	GRADproduct += 0.5*( (2*ai*GRADa*dThetay[iter1->first]*dThetay[iter1->first]*Lij)  +  (2*(GRADa*bi + GRADb*ai)*(dThetay[iter1->first]*dThetay[iter1->first]*Bij  +  dThetay[iter1->first]*dy[iter1->first]*Lij)) + (2*bi*GRADb*(2*dThetay[iter1->first]*dy[iter1->first]*Bij + dy[iter1->first]*dy[iter1->first]*Lij + dThetay[iter1->first]*dThetay[iter1->first]*Cij )) + (2*ci*GRADc*dy[iter1->first]*dy[iter1->first]*Cij) + (2*(GRADb*ci + GRADc*bi)*(dy[iter1->first]*dy[iter1->first]*Bij  +  dThetay[iter1->first]*dy[iter1->first]*Cij)) + (2*dThetay[iter1->first]*dy[iter1->first]*Bij*(GRADa*ci + GRADc*ai)) );


	  Gradient[iter->first] += ( (1./Lambda[iter->first]) + (0.5*PreFactorEng[iter1->first]*(product - trace + Lambda[iter->first]*(GRADproduct - GRADtrace)))  );


      }
   }
   for(voxelNum=0;voxelNum<fVoxelCount;voxelNum++) {
      gradientNum=Gradient[voxelNum];
      lambdaNum=Lambda[voxelNum];
#ifdef MAKE_GRAD_TREE
      gradTree->Fill();
#endif
   }
#ifdef MAKE_GRAD_TREE
   gradTree->AutoSave();
#endif
   printf("Gradient Fill: %f\n", ((double)(clock() - GradientStart) / CLOCKS_PER_SEC));
}


//cost Function
double LambdaPcaTreeLooperMLSD::Cost(double Alpha, int first, int last, int iteration){

   double cost=0;
   static double mat_W[2][2];
   static double mat_D[2][2]; //SigmaD
   static double mat_R[2][2]; // SigmaD*W
   static double mat_Tot[2][2]; // SigmaD*W*SigmaD
   static double trace;
   static double product;

  static int loopCount=0;
   static double myAlpha=Alpha;
   static double newLambda;
   static int muonNum;
   static double muonCost;
   static double costVoxel;

#ifdef MAKE_DEBUG_TREE
   static double sNum;
   static double mNum;
   static double W00voxel;
   static double W01voxel;
   static double W11voxel;
   static double thetaNum;
   static double dxNum;
   static int voxelNum;
   static double gradientNum;
   static double lambdaNum;
static double voxel1;
   static TFile *fpCost = new TFile("/unix/creamtea/sfayer/temp/cost.root","RECREATE");
   static TTree *costTree = 0;//(TTree*) fpCost->Get("costTree");
   static TTree *iterationTree = 0;
   if(!costTree) {
      costTree= new TTree("costTree","costTree");
      costTree->Branch("alpha",&myAlpha,"alpha/D");
      costTree->Branch("loopCount",&loopCount,"loopCount/I");
      costTree->Branch("muonNum",&muonNum,"muonNum/I");
      costTree->Branch("sNum",&sNum,"sNum/D");
      costTree->Branch("mNum",&mNum,"mNum/D");
      costTree->Branch("lambdaNum",&lambdaNum,"lambdaNum/D");
      costTree->Branch("gradientNum",&gradientNum,"gradientNum/D");
      costTree->Branch("W00voxel",&W00voxel,"W00voxel/D");
      costTree->Branch("W01voxel",&W01voxel,"W01voxel/D");
      costTree->Branch("W11voxel",&W11voxel,"W11voxel/D");
      costTree->Branch("thetaNum",&thetaNum,"thetaNum/D");
      costTree->Branch("dxNum",&dxNum,"dxNum/D");
      costTree->Branch("costVoxel",&costVoxel,"costVoxel/D");
      costTree->Branch("trace",&trace,"trace/D");
      costTree->Branch("product",&product,"product/D");
      costTree->Branch("voxelNum",&voxelNum,"voxelNum/I");
      }
   if(!iterationTree) {
     iterationTree= new TTree("iterationTree","iterationTree");
	iterationTree->Branch("voxel1",&voxel1,"voxel1/D");
}
#endif

   myAlpha=Alpha;
   static int madeSM=0;
if(!madeSM){
      Sj = new double[fVoxelCount];
      M = new double[fVoxelCount];
 avetheta = new double[fVoxelCount];
 mintheta = new double[fVoxelCount];
 maxtheta = new double[fVoxelCount];
 avedist = new double[fVoxelCount];
 mindist = new double[fVoxelCount];
 maxdist = new double[fVoxelCount];
 rmstheta = new double[fVoxelCount];
 rmsdist = new double[fVoxelCount];
 nomuons = new int[fVoxelCount];
memset(avetheta,0,fVoxelCount*sizeof(double)); 
memset(avedist,0,fVoxelCount*sizeof(double)); 
memset(mintheta,10000,fVoxelCount*sizeof(double)); 
memset(mindist,10000,fVoxelCount*sizeof(double)); 
memset(maxtheta,0,fVoxelCount*sizeof(double)); 
memset(maxdist,0,fVoxelCount*sizeof(double));
memset(rmstheta,0,fVoxelCount*sizeof(double)); 
memset(rmsdist,0,fVoxelCount*sizeof(double)); 
memset(nomuons,0,fVoxelCount*sizeof(int)); 
	madeSM++;
 }

   memset(Sj,0,fVoxelCount*sizeof(double)); 
   memset(M,0,fVoxelCount*sizeof(double)); 
   for(std::map <int, std::map<int, double> >::iterator muonIter=L.begin();
       muonIter!=L.end();
       muonIter++) {
      muonNum=muonIter->first;

	 mat_D[0][0] = aD[muonNum];
	 mat_D[0][1] = bD[muonNum];
	 mat_D[1][0] = bD[muonNum];
	 mat_D[1][1] = cD[muonNum];

      for(std::map<int,double>::iterator iter = (muonIter->second).begin(); iter != (muonIter->second).end(); ++iter){
if (iteration==0){
if(abs(dThetax[muonNum])<mintheta[iter->first])mintheta[iter->first]=abs(dThetax[muonNum]);
if(abs(dThetay[muonNum])<mintheta[iter->first])mintheta[iter->first]=abs(dThetay[muonNum]);
if(abs(dThetax[muonNum])>maxtheta[iter->first])maxtheta[iter->first]=abs(dThetax[muonNum]);
if(abs(dThetay[muonNum])>maxtheta[iter->first])maxtheta[iter->first]=abs(dThetay[muonNum]);
if(abs(dx[muonNum])<mindist[iter->first])mindist[iter->first]=abs(dx[muonNum]);
if(abs(dy[muonNum])<mindist[iter->first])mindist[iter->first]=abs(dy[muonNum]); 
if(abs(dx[muonNum])>maxdist[iter->first])maxdist[iter->first]=abs(dx[muonNum]);
if(abs(dy[muonNum])>maxdist[iter->first])maxdist[iter->first]=abs(dy[muonNum]); 
avedist[iter->first]+=((dx[muonNum]+dy[muonNum])/2);
avetheta[iter->first]+=((dThetax[muonNum]+dThetay[muonNum])/2);
rmsdist[iter->first]+=(((dx[muonNum]*dx[muonNum])+(dy[muonNum]*dy[muonNum]))/2);
rmstheta[iter->first]+=(((dThetax[muonNum]*dThetax[muonNum])+(dThetay[muonNum]*dThetay[muonNum]))/2);
nomuons[iter->first]+=1;
}
	 newLambda = Lambda[iter->first];// - Alpha*Gradient[iter->first];
	 mat_W[0][0] = A[muonNum][iter->first];
	 mat_W[0][1] = B[muonNum][iter->first];
	 mat_W[1][0] = B[muonNum][iter->first];
	 mat_W[1][1] = C[muonNum][iter->first];
	  LambdaPcaTreeLooperMLSD::mm_mul (mat_D, mat_W, mat_R);
	 trace = mat_R[0][0] + mat_R[1][1];
	  LambdaPcaTreeLooperMLSD::mm_mul (mat_R, mat_D, mat_Tot);
	 product = (mat_Tot[0][0]*dThetax[muonNum]*dThetax[muonNum]) + (mat_Tot[1][0]*dThetax[muonNum]*dx[muonNum]) + (mat_Tot[0][1]*dThetax[muonNum]*dx[muonNum]) + (mat_Tot[1][1]*dx[muonNum]*dx[muonNum]);
	 
	 if(newLambda<LAMBDA_AIR) newLambda=LAMBDA_AIR;
	 if(isnan(newLambda)) newLambda=LAMBDA_AIR;
	 
 Sj[iter->first] += 0.5*(2*newLambda + (newLambda*newLambda*PreFactorEng[muonNum]*(product-trace)));

//now add the Y component
	 product = (mat_Tot[0][0]*dThetay[muonNum]*dThetay[muonNum]) + (mat_Tot[1][0]*dThetay[muonNum]*dy[muonNum]) + (mat_Tot[0][1]*dThetay[muonNum]*dy[muonNum]) + (mat_Tot[1][1]*dy[muonNum]*dy[muonNum]);

 Sj[iter->first] += 0.5*(2*newLambda + (newLambda*newLambda*PreFactorEng[muonNum]*(product-trace)));
	 M[iter->first] +=1.;

if(newLambda>1 || (iter->first >= 400&&iter->first<450&&M[iter->first]>0)) {
  std::cout << "voxelnum: " << iter->first << "    newLambda: " << newLambda << "    product: " << product << "    trace: " << trace << "    S: " << Sj[iter->first] << " M:   " << M[iter->first] << " mat_D[0][0] = " << mat_D[0][0] << " mat_D[0][1] = " << mat_D[0][1] << "   mat_D[1][0] = " << mat_D[1][0] << "   mat_D[1][1] = " << mat_D[1][1] << "   mat_W[0][0] = " << mat_W[0][0] << "   mat_W[0][1] = " << mat_W[0][1] << "   mat_W[1][0] = " << mat_W[1][0] << "   mat_W[1][1] = " << mat_W[1][1] <<" dThetay[muonNum]: " << dThetay[muonNum] <<   "dy[muonNum]: " << dy[muonNum] << " dThetax[muonNum]: " << dThetax[muonNum] <<   "dx[muonNum]: " << dx[muonNum] << " Prefactor: " << PreFactorEng[muonNum] << endl;
	   //newLambda=LAMBDA_AIR;
	 }

#ifdef MAKE_DEBUG_TREE
  W00voxel =  mat_W[0][0];
   W01voxel =  mat_W[1][0];
   W11voxel =  mat_W[1][1];
   thetaNum = dThetax[muonNum];
   dxNum = dx[muonNum];
   voxelNum = iter->first;
   lambdaNum = newLambda;
   sNum = Sj[iter->first];
   mNum = M[iter->first];
if (iter->first==998062){
voxel1 = newLambda;
iterationTree->Fill();
}
   costTree->Fill();
#endif	

	}

   }

  for(int voxId=0;voxId<fVoxelCount;voxId++) {
	costVoxel = M[voxId]*TMath::Log(Lambda[voxId]) + (Sj[voxId]*0.5/Lambda[voxId]);
	cost+= costVoxel;
}

   loopCount++;

#ifdef MAKE_DEBUG_TREE
   iterationTree->AutoSave();
   costTree->AutoSave();
#endif
   return cost;
}


//Draw Slices
void LambdaPcaTreeLooperMLSD::DrawSlices(int topSlice, int bottomSlice, int Nx, int Ny, int Nz, char* FileNameLambda){
   if(topSlice > bottomSlice && topSlice <= Nz && bottomSlice >= 0){

      TFile *fpOut = new TFile(FileNameLambda,"RECREATE");
      fpOut->Write();
      TCanvas *canProj = new TCanvas();
      canProj->Divide(3,3);
    
      //Slice start height.
//      int z = 50;
      int z = 0;
      char HistName[80];
      for(int SliceNo = 1; SliceNo < Nz; SliceNo++){
	 sprintf(HistName,"Slice_%d",SliceNo);
	 TH2F *CurrentLSlice = new TH2F(HistName,HistName, Nx, 0.0, (double) Nx, Ny, 0.0, (double) Ny);

	 for(int x = 1; x <= Nx; x++){
	    for(int y = 1; y <= Ny; y++){
	       CurrentLSlice->SetBinContent(x,y,Lambda[(z+(SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
	    }
	 }

	 canProj->cd(SliceNo);
	 CurrentLSlice->Draw("colz");
	 CurrentLSlice->SetStats(0);
	 CurrentLSlice->GetXaxis()->SetTitle("X ID");
	 CurrentLSlice->GetYaxis()->SetTitle("Y ID");
	 CurrentLSlice->Write();
 
      }
      fpOut->Close();
    
   }
  
}

//Draw Muons and Input data
void LambdaPcaTreeLooperMLSD::DrawMuons(int outputSlice, int Nx, int Ny, int Nz, char* FileNameMuon) {
      TFile *muOut = new TFile(FileNameMuon,"RECREATE");
      muOut->Write();

      TCanvas *canProj = new TCanvas();
      canProj->Divide(3,3);
      char* data[9];
      data[0] = "average_theta";
	data[1] = "min_theta";
	data[2] = "max_theta";
      data[3] = "average_displacement";
	data[4] = "min_displacement";
	data[5] = "max_displacement";
	data[6] = "number_of_muons";
	data[7] = "rms_theta";
	data[8] = "rms_displacement";

      char HistName[80];
      int SliceNo = outputSlice;
for (int d=0;d<9;d++){
	 sprintf(HistName,"%s_%d",data[d],SliceNo);
	 TH2F *CurrentMSlice = new TH2F(HistName,HistName, Nx, 0.0, (double) Nx, Ny, 0.0, (double) Ny);
	 for(int x = 0; x <= Nx; x++){
	    for(int y = 0; y <= Ny; y++){
	       if(d==0) {
avetheta[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]/=nomuons[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x];
CurrentMSlice->SetBinContent(x,y,avetheta[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
}
		if(d==1) CurrentMSlice->SetBinContent(x,y,mintheta[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
		if(d==2) CurrentMSlice->SetBinContent(x,y,maxtheta[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
 	       if(d==3) {
avedist[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]/=nomuons[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x];
CurrentMSlice->SetBinContent(x,y,avedist[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
}
		if(d==4) CurrentMSlice->SetBinContent(x,y,mindist[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
		if(d==5) CurrentMSlice->SetBinContent(x,y,maxdist[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
 	       if(d==6) CurrentMSlice->SetBinContent(x,y,nomuons[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
if(d==7) {
rmstheta[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]/=nomuons[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x];
CurrentMSlice->SetBinContent(x,y,rmstheta[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
}
if(d==8) {
rmsdist[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]/=nomuons[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x];
CurrentMSlice->SetBinContent(x,y,rmsdist[((SliceNo-1))*Nx*Ny + (y-1)*Nx + x]);
}
	    }
	 }
	 canProj->cd(d);
	 CurrentMSlice->Draw("colz");
	 CurrentMSlice->SetStats(0);
	 CurrentMSlice->GetXaxis()->SetTitle("X ID");
	 CurrentMSlice->GetYaxis()->SetTitle("Y ID");
	 CurrentMSlice->Write();
 
      }
      muOut->Close();
  
}
