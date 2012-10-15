#define LambdaPcaTreeLooper_cxx
#include "LambdaPcaTreeLooper.h"
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

#define LAMBDA_AIR (1./3.039e4)

void LambdaPcaTreeLooper::Loop()
{
   //   In a ROOT session, you can do:
   //      Root > .L LambdaPcaTreeLooper.C
   //      Root > LambdaPcaTreeLooper t
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

void LambdaPcaTreeLooper::FillPosHist(TH3F *histPos, Double_t thetaCut)
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





/////////////////////////Will's Functions////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


double prefactor = 13.6*13.6/1000000;
//double prefactor = 13.6*13.6;



//Fill Matrix L with muons from FIRST to LAST, and make the S vector. Ni is the number of voxels in the i direction.
void LambdaPcaTreeLooper::SLFill(int first, int last, int Nx, int Ny, int Nz){
   //At some point might update this so that one can specify the first and last and only have array that size
   fNumMuons=last;

   S = new double[fNumMuons];
   PreFactorEng = new double[fNumMuons];
   memset(S,0,fNumMuons*sizeof(double));
   memset(PreFactorEng,0,fNumMuons*sizeof(double));

#ifdef MAKE_DEBUG_TREE
   Int_t muonNum,voxelNum;
   Double_t lij;
   TFile *fpLij = new TFile("/unix/creamtea/sfayer/temp/lijOut.root","RECREATE");
   TTree *lijTree = new TTree("lijTree","lijTree");
   lijTree->Branch("muonNum",&muonNum,"muonNum/I");
   lijTree->Branch("voxelNum",&voxelNum,"voxelNum/I");
   lijTree->Branch("lij",&lij,"lij/D");
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
		  
		  //Hack because z=6500 is outside the box
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
			xCurrent = Xleave;
			yCurrent = Yleave;
			zCurrent = ZExit;
		
		     }
	      
		     //Get out of loop
		     if(CurID == VoxelEnd) part = 1;
      
		     //Get out of loop AND prevent next loop from being ru; doubleing point error somewhere.
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
		     // Fill in the PreFactorEng array with true energy of this particular event
		     //		     PreFactorEng[jentry]=13.6*13.6/(intEng*intEng);
		     PreFactorEng[jentry]=13.6*13.6/(1000*1000);
		     
		     //Now copy tempMap to the big L map of maps
		     //		     std::cout << tempMap.size() << "\n";
#ifdef MAKE_DEBUG_TREE
		     muonNum=jentry;
#endif		     

		     for(std::map<int, double>::iterator tempIt=tempMap.begin();
			 tempIt!=tempMap.end();tempIt++) {
			L[jentry][tempIt->first]=tempIt->second;
#ifdef MAKE_DEBUG_TREE
			lij=tempIt->second;
			voxelNum=tempIt->first;
			lijTree->Fill();
#endif
		     }
		     



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
void LambdaPcaTreeLooper::LambdaFill(int VoxelCount){
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
void LambdaPcaTreeLooper::LambdaAlpha(double Alpha){
   clock_t LambdaStart = clock();

   double newLambda;

   //  for(std::map<int,double>::iterator iter = Lambda.begin(); iter != Lambda.end(); ++iter){
   for(int voxId=0;voxId<fVoxelCount;voxId++) {
      newLambda = Lambda[voxId] - Alpha*Gradient[voxId];
      //      if(newLambda < 0.000003067) newLambda = 0.000003067;
      if(newLambda < LAMBDA_AIR) newLambda = LAMBDA_AIR;
      Lambda[voxId] = newLambda;
   }
   printf("Lambda Alpha Fill: %f\n", ((double)(clock() - LambdaStart) / CLOCKS_PER_SEC));
}


void LambdaPcaTreeLooper::SigmaFill(){
  
   clock_t SigmaStart = clock();
#ifdef MAKE_DEBUG_TREE
   static TFile *fpSigma = new TFile("/unix/creamtea/sfayer/temp/sigmaOut.root","RECREATE");
   static TTree *sigmaTree = 0;// (TTree*) fpSigma->Get("sigmaTree");
#endif
   static int doneInit=0;
   static Int_t muonNum=0;
   static Double_t sigmaNum=0;
   static Double_t prefactorNum=0;
   static Double_t sNum=0;
   if(!doneInit) {
      Sigma = new double[fNumMuons];
#ifdef MAKE_DEBUG_TREE
      if(!sigmaTree) {
	 sigmaTree = new TTree("sigmaTree","sigmaTree");
	 sigmaTree->Branch("muonNum",&muonNum,"muonNum/I");
	 sigmaTree->Branch("loopNum",&doneInit,"loopNum/I");
	 sigmaTree->Branch("prefactorNum",&prefactorNum,"prefactorNum/D");
	 sigmaTree->Branch("sigmaNum",&sigmaNum,"sigmaNum/D");
	 sigmaTree->Branch("sNum",&sNum,"sNum/D");
      }
      else {
	 sigmaTree->SetBranchAddress("muonNum",&muonNum);
	 sigmaTree->SetBranchAddress("loopNum",&doneInit);
	 sigmaTree->SetBranchAddress("sigmaNum",&sigmaNum);
	 sigmaTree->SetBranchAddress("sNum",&sNum);
      }
#endif
   }
   doneInit++;
   //   Sigma.clear();
   memset(Sigma,0,fNumMuons*sizeof(double));

      

   for(std::map<int,map<int,double> >::iterator iter1 = L.begin(); iter1 != L.end(); ++iter1){ //loop over muons
      Sigma[iter1->first] = 0;
      for(std::map<int,double>::iterator iter2 = (iter1->second).begin(); iter2 != (iter1->second).end(); ++iter2){
	 Sigma[iter1->first] += iter2->second * Lambda[iter2->first];
      }
      //      Sigma[iter1->first] *= prefactor;///(intEng*intEng);    
      //    Sigma[iter1->first] *= prefactor;
      if(iter1->first==0) {
	 std::cout << "Before factor: " << Sigma[iter1->first] << "\n";
      }
      Sigma[iter1->first] *= PreFactorEng[iter1->first];
      prefactorNum=PreFactorEng[iter1->first];
      Sigma[iter1->first] += 0.000049;
      muonNum=iter1->first;   
      sigmaNum=Sigma[iter1->first];
      sNum=S[iter1->first];

      if(iter1->first==0) {
	 std::cout << "After factor: " << Sigma[iter1->first] << "\n";
      }

#ifdef MAKE_DEBUG_TREE
      sigmaTree->Fill();
#endif
   }
#ifdef MAKE_DEBUG_TREE
   sigmaTree->AutoSave();
#endif

   printf("Sigma Fill: %f\n", ((double)(clock() - SigmaStart) / CLOCKS_PER_SEC));
}


//Gradient Vector.
void LambdaPcaTreeLooper::GradientFill(){
   clock_t GradientStart = clock();
#ifdef MAKE_DEBUG_TREE
   static TFile *fpGrad = new TFile("/unix/creamtea/sfayer/temp/gradOut.root","RECREATE");
   static TTree *gradTree = (TTree*) fpGrad->Get("gradTree");
#endif
   static Int_t voxelNum=0;
   static Double_t lambdaNum=0;
   static Double_t gradientNum=0;
   static int madeGradient=0;
   if(!madeGradient) {
      Gradient = new double [fVoxelCount];
#ifdef MAKE_DEBUG_TREE
      if(!gradTree) {
	 gradTree= new TTree("gradTree","Tree of voxel gradients");
	 gradTree->Branch("voxelNum",&voxelNum,"voxelNum/I");
	 gradTree->Branch("loopNum",&madeGradient,"loopNum/I");
	 gradTree->Branch("lambdaNum",&lambdaNum,"lambdaNum/D");
	 gradTree->Branch("gradientNum",&gradientNum,"gradientNum/D");
      }      
#endif
   }
   madeGradient++;
   
   memset(Gradient,0,fVoxelCount*sizeof(double)); 
   //  Gradient.clear();

   for(std::map<int,map<int,double> >::iterator iter1 = L.begin(); iter1 != L.end(); ++iter1){
       for(std::map<int,double>::iterator iter = (iter1->second).begin(); iter != (iter1->second).end(); ++iter){
	  Gradient[iter->first] += iter->second * ((1 - (S[iter1->first]/Sigma[iter1->first])))/Sigma[iter1->first];
      }
   }
   for(voxelNum=0;voxelNum<fVoxelCount;voxelNum++) {
      gradientNum=Gradient[voxelNum];
      lambdaNum=Lambda[voxelNum];
#ifdef MAKE_DEBUG_TREE
      gradTree->Fill();
#endif
   }
#ifdef MAKE_DEBUG_TREE
   gradTree->AutoSave();
#endif
   printf("Gradient Fill: %f\n", ((double)(clock() - GradientStart) / CLOCKS_PER_SEC));
}


//cost Function
double LambdaPcaTreeLooper::Cost(double Alpha, int first, int last){

   double cost=0;
  
   static int loopCount=0;
   static double myAlpha=Alpha;
   static double newLambda;
   static double SigmaA;
   static double muonCost;
   static int muonNum;
   static double sNum;
#ifdef MAKE_DEBUG_TREE
   static TFile *fpCost = new TFile("/unix/creamtea/sfayer/temp/costOut.root","RECREATE");
   static TTree *costTree =0;
   
   if(!costTree) {
      costTree= new TTree("costTree","costTree");
      costTree->Branch("alpha",&myAlpha,"alpha/D");
      costTree->Branch("sigmaA",&SigmaA,"sigmaA/D");
      costTree->Branch("loopCount",&loopCount,"loopCount/I");
      costTree->Branch("muonNum",&muonNum,"muonNum/I");
      costTree->Branch("muonCost",&muonCost,"muonCost/D");
      costTree->Branch("sNum",&sNum,"sNum/D");
   }
#endif
   myAlpha=Alpha;


   //  map<int,double> SigmaA; 
    
   

   //SigmaA[]
   //   for(int i = first; i < last; ++i){ //loop over muons
   for(std::map <int, std::map<int, double> >::iterator muonIter=L.begin();
       muonIter!=L.end();
       muonIter++) {
      muonNum=muonIter->first;

      SigmaA=0;
      for(std::map<int,double>::iterator iter = (muonIter->second).begin(); iter != (muonIter->second).end(); ++iter){
	 //	 std::cout << i << "\t" << iter->first << "\n";
	 newLambda = Lambda[iter->first] - Alpha*Gradient[iter->first];
	 //	 if(newLambda < 0.000003067) newLambda = 0.000003067;
	 if(newLambda<LAMBDA_AIR) newLambda=LAMBDA_AIR;
	 //      SigmaA[muonNum] += iter->second * newLambda;
	 SigmaA += iter->second * newLambda;
      }
      //    SigmaA[muonNum] *= prefactor;
//      SigmaA *= prefactor; // What is the prefactor
//      SigmaA *= prefactor; ///(intEng*intEng); // What is the prefactor
      SigmaA *= PreFactorEng[muonNum];
      SigmaA += 0.000049;  //What is this?
      muonCost=S[muonNum]/SigmaA + TMath::Log(SigmaA);
      sNum=S[muonNum];
      //    cost += S[muonNum]/SigmaA[muonNum] + TMath::Log(SigmaA[muonNum]);
      cost += muonCost;
#ifdef MAKE_DEBUG_TREE
      costTree->Fill();
#endif
   }
   loopCount++;
#ifdef MAKE_DEBUG_TREE
   costTree->AutoSave();
#endif
   return cost;
}


//Draw Slices
void LambdaPcaTreeLooper::DrawSlices(int topSlice, int bottomSlice, int Nx, int Ny, int Nz, char* FileNameLambda){
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

	 for(int x = 0; x <= Nx; x++){
	    for(int y = 0; y <= Ny; y++){
	       CurrentLSlice->SetBinContent(x,y,Lambda[(z+SliceNo)*Nx*Ny + y*Nx + x]);
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
