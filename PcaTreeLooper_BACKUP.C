#define PcaTreeLooper_cxx
#include "PcaTreeLooper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <map>
//For Calculation Time
#include <time.h>
#include <stdio.h>
#include <vector>

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

//void PcaTreeLooper::frprmn(float p[], int n, float ftol, int *iter, float *fret, float (*func)(float []), void (*dfunc)(float[], float []))
void PcaTreeLooper::frprmn()
{
  
  //To Make Density Array And L Matrix -> 400x400x400 with 100000 muons = 280 seconds.
  //The Rate seems linear with the number of muons used.
  //The effect on computation time due to voxel count is seemingly not simple.
  
  
  //For Calculation Time
  clock_t Lijstart = clock();
  
  //Part from the skeleton code for a loop, no point doing anything if it's not going to work.
  if (fChain == 0) return;
  
  
  
  //Size of cuboid
  double MaxX = 6500;
  double MinX = -6500;
  double MaxY = 6500;
  double MinY = -6500;
  double MaxZ = 6500;
  double MinZ = -6500;
  
  
  //Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nentries = 100000; //For testing, so we don't take all the muons
  
  
  int Nx = 15; //Number of voxels in x direction
  int Ny = 15; //Number of voxels in y direction
  int Nz = 15; //Number of voxels in z direction
  
  //Voxel dimensions (multiplication faster)
  double Vx = (MaxX-MinX)/Nx;
  double Vy = (MaxY-MinY)/Ny;
  double Vz = (MaxZ-MinZ)/Nz;
  
  //Populate Voxel scattering densities & an empty gradient thing.
  map<int,double> Lambda;
  map<int,double> gradient;
  //Starting Position for Voxel Densities.
  for(int AA=0;AA<Nx;AA++){
    for(int BB=0;BB<Ny;BB++){
      for(int CC=0;CC<Nz;CC++){
	int VoxelID = AA + Nx*BB + Nx*Ny*CC;
	Lambda[VoxelID] = 0.000003067;
	gradient[VoxelID] = 1;
	}
    }
  }
  
  //Voxel ID: 
  // k = Floor(ID/Nz)
  // j = Floor((ID - k * Nz)/Ny)
  // i = (ID - k*Nz - j*Ny)/Nx
  
  //////////////////populate Lij. i = muon ID, j = voxel ID
  
  
  
  
  //Scattering Angle per muon
  map<int,double> S;
  
  
  //...................................................................................... -> Make Lij Matrix <- .............
  
  
  //Define Lij Matrix kind of thing
  //Reference -> http://groups.google.com/group/comp.lang.c++.moderated/browse_thread/thread/bb499f84b93c4343?pli=1
  typedef map<int, map<int, double> > MapMap; 
  MapMap L;
  // N * L[i][j] <- Et Voila. Lij.
  
 
  
  //Initialise some variables so we don't get loads of errors.
  int EndX, EndY, EndZ, EndID, CurrentX, CurrentY, CurrentZ, vID, special, xzCase, yzCase, xSign, ySign, errorcounter;
  errorcounter = 0;
  double xEntry, yEntry, zEntry, xExit, yExit, zExit, zExitX, zExitY, Length;
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  //error marjin for float inequalities.
  float Epsilon = 0.0000001;
  xSign=ySign=xzCase=yzCase=0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //If PCA is in box. If not, we don't care.
    if(gotRecoPCA=1){
      //Is PCA inside Box?
      if(xPosTrue >= MinX && xPosTrue <= MaxX && yPosTrue >= MinY && yPosTrue <= MaxY && zPosTrue >= MinZ && zPosTrue <= MaxZ){
	//Scattering Angle of muon.
	S[(int) jentry] = thetaTrue;
	for(int part=0; part < 2; part++){
	  //part==0 for before pca, part==1 after pca.
	  if(part == 0){
	    xEntry = xzGradTrue[0] * MaxZ + xzCutTrue[0];
	    yEntry = yzGradTrue[0] * MaxZ + yzCutTrue[0];
	    zEntry = MaxZ;
	    EndX = (int) floor((xPosTrue-MinX)/Vx);
	    EndY = (int) floor((yPosTrue-MinY)/Vy);
	    EndZ = (int) floor((zPosTrue-MinZ)/Vz);
	  }
	  else{
	    xEntry = xPosTrue;
	    yEntry = yPosTrue;
	    zEntry = zPosTrue;
	    EndX = (int) floor(((xzGradTrue[1] * MinZ + xzCutTrue[1])-MinX)/Vx);
	    EndY = (int) floor(((yzGradTrue[1] * MinZ + yzCutTrue[1])-MinY)/Vy);
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
	  
	  //stupid "may not be initialised" errors.
	  xExit = yExit = zExit = 0;
	  

	  //Loop over voxels for this muon track.
	  special=0;
	  //cout << "********************************" << xzGradTrue[part] << " : " << yzGradTrue[part] << " ||| " << EndID << endl;
	  while(vID != EndID && special < 30){
	    //cout << CurrentX << " " << CurrentY << " " << CurrentZ << endl;
	    zExitX = ((CurrentX+xzCase)*Vx+MinX-xzCutTrue[part])/xzGradTrue[part];
	    zExitY = ((CurrentY+yzCase)*Vy+MinY-yzCutTrue[part])/yzGradTrue[part];

	    zExit = max(max(zExitX,zExitY),CurrentZ*Vz+MinZ);
	    xExit = xzGradTrue[part]*zExit+xzCutTrue[part];
	    yExit = yzGradTrue[part]*zExit+yzCutTrue[part];

	    Length = pow(pow(xExit-xEntry,2)+pow(yExit-yEntry,2)+pow(zExit-zEntry,2),0.5);
	    if(Length < pow(pow(Vx,2)+pow(Vy,2)+pow(Vz,2),0.5)){L[jentry][vID] = Length;}
	    else{errorcounter++;}

	    if(zExit <= zExitX+Epsilon && zExit >= zExitX-Epsilon){CurrentX += xSign;}
	    if(zExit <= zExitY+Epsilon && zExit >= zExitY-Epsilon){CurrentY += ySign;}
	    if(zExit <= CurrentZ*Vz+MinZ+ Epsilon&& zExit >= CurrentZ*Vz+MinZ-Epsilon){CurrentZ -= 1;}

	    xEntry = xExit;
	    yEntry = yExit;
	    zEntry = zExit;
	    vID = CurrentX + CurrentY*Nx + CurrentZ*Nx*Ny;

	    special++;
	  }
	  if(vID == EndID && part == 0){
	    //cout << CurrentX << " " << CurrentY << " " << CurrentZ << endl;
	    //Loop for PCA Voxel

	    //Up to PCA
	    Length = pow(pow(xExit-xEntry,2)+pow(yExit-yEntry,2)+pow(zExit-zEntry,2),0.5);
	    if(Length < pow(pow(Vx,2)+pow(Vy,2)+pow(Vz,2),0.5)){L[jentry][vID] = Length;}
	    else{errorcounter++;}

	    //Other side of PCA
	    zExitX = (xPosTrue-xzCutTrue[part])/xzGradTrue[part];
	    zExitY = (yPosTrue-yzCutTrue[part])/yzGradTrue[part];

	    zExit = max(max(zExitX,zExitY),CurrentZ*Vz+MinZ);
	    xExit = xzGradTrue[part]*zExit+xzCutTrue[part];
	    yExit = yzGradTrue[part]*zExit+yzCutTrue[part];

	    //Add on second part of length in PCA Voxel
	    Length = pow(pow(xExit-xEntry,2)+pow(yExit-yEntry,2)+pow(zExit-zEntry,2),0.5);
	    if(Length < pow(pow(Vx,2)+pow(Vy,2)+pow(Vz,2),0.5)){L[jentry][vID] += Length;}
	    else{errorcounter++;}
	    
	    if(zExit <= zExitX+Epsilon && zExit >= zExitX-Epsilon){CurrentX += xSign;}
	    if(zExit <= zExitY+Epsilon && zExit >= zExitY-Epsilon){CurrentY += ySign;}
	    if(zExit <= CurrentZ*Vz+MinZ+ Epsilon&& zExit >= CurrentZ*Vz+MinZ-Epsilon){CurrentZ -= 1;}

	    xEntry = xExit;
	    yEntry = yExit;
	    zEntry = zExit;
	    vID = CurrentX + CurrentY*Nx + CurrentZ*Nx*Ny;


	  }else if(vID == EndID && part == 1){
	    //cout << CurrentX << " " << CurrentY << " " << CurrentZ << endl;
	    //Loop for final voxel. Do it here rather than in the main loop so the PCA voxel is easier to do.
	    zExit = MinZ;
	    xExit = xzGradTrue[part]*zExit+xzCutTrue[part];
	    yExit = yzGradTrue[part]*zExit+yzCutTrue[part];
	    Length = pow(pow(xExit-xEntry,2)+pow(yExit-yEntry,2)+pow(zExit-zEntry,2),0.5);
	    if(Length < pow(pow(Vx,2)+pow(Vy,2)+pow(Vz,2),0.5)){L[jentry][vID] += Length;}
	    else{errorcounter++;}	  }
	}
	//Move onto next muon.
	
      }
    }
 
  
      
  }


//...................................................................................... -> Lij Matrix Made <- .............


 printf("Lij Fill: %f\n", ((double)clock() - Lijstart) / CLOCKS_PER_SEC);


 //For Calculation Time
 clock_t Gradientstart = clock();
 
 double Gradient[Nx*Ny*Nz];

//Gradient Vector
 int minVox = 0;
 int maxVox = Nx * Ny * Nz;
 TH1F *MuonsPerVoxel = new TH1F("MuonsPerVoxel","Muons Per Voxel",maxVox-minVox,minVox,maxVox);
 TH1F *TotalLengthPerVoxel = new TH1F("TotalLengthPerVoxel","Total Pathlength Per Voxel",maxVox-minVox,minVox,maxVox);
 TH1F *VoxelGradient = new TH1F("VoxelGradient","Gradient of each Voxel",maxVox-minVox,minVox,maxVox);

 map<int,double> Sigma;
 for (int iter1 = 0; iter1 < nentries; ++iter1) {
   for(std::map<int,double>::iterator iter2 = L[iter1].begin(); iter2 != L[iter1].end(); ++iter2){
     Sigma[iter1] += iter2->second * (Lambda[iter2->first]);
   }
 }

 for(int iter1 = 0; iter1 < (int)nentries; ++iter1){
   for(std::map<int,double>::iterator iter2 = L[iter1].begin(); iter2 != L[iter1].end(); ++iter2){
     int VoxelID = iter2->first;
     double PathLength = iter2->second;
     if(VoxelID >= minVox && VoxelID <= maxVox){
       Gradient[VoxelID] = PathLength *(1-S[iter1]*S[iter1]/Sigma[iter1]);
       VoxelGradient->Fill(VoxelID,Gradient[VoxelID]);
       MuonsPerVoxel->Fill(VoxelID,1);
       TotalLengthPerVoxel->Fill(VoxelID,PathLength);
     }
   }
 }

TH1F *CostAlpha = new TH1F("CostAlpha","Psi(Lambda(0)-Alpha*Grad(Lambda(0)))",maxVox-minVox,minVox,maxVox);
 //Psi(Lambda(0) - Alpha*Grad(Lambda(0)))
 for(int b=0;b<10;b++){
   float Alpha = b;
   
   //Work out sigma for this alpha.
   for (int iter1 = 0; iter1 < nentries; ++iter1) {
     for(std::map<int,double>::iterator iter2 = L[iter1].begin(); iter2 != L[iter1].end(); ++iter2){
       double newLambda =Lambda[iter2->first]-Alpha*Gradient[iter2->first;
       if(newLambda < 0.000003067 ){newLambda = 0.000003067;}
       Sigma[iter1] += iter2->second * (newLambda);
     }
   }
   for(int VoxelID=minVox;VoxelID<maxVox;VoxelID++){
     CostAlpha->Fill(VoxelID,S[VoxelID]*S[VoxelID]/Sigma[VoxelID] + log(Sigma[VoxelID]));
   }

   


 }
 
 


 printf("Functions: %f\n", ((double)clock() - Gradientstart) / CLOCKS_PER_SEC);
 cout << errorcounter << endl;


  TFile *fpOu1 = new TFile("MuonPath1.root","RECREATE");
  fpOu1->Write();
  TCanvas *canProj = new TCanvas();
  canProj->Divide(1,3);
  canProj->cd(1);
  CostAlpha->Draw();
  CostAlpha->SetStats(0);
  CostAlpha->GetXaxis()->SetTitle("Voxel ID");
  canProj->cd(2);
  TotalLengthPerVoxel->Draw();
  TotalLengthPerVoxel->SetStats(0);
  TotalLengthPerVoxel->GetXaxis()->SetTitle("Voxel ID");
  canProj->cd(3);
  VoxelGradient->Draw();
  VoxelGradient->SetStats(0);
  VoxelGradient->GetXaxis()->SetTitle("Voxel ID");


  //  TH1F *VoxelGradient = new TH1F("VoxelGradient","Gradient at each voxel",maxVox,0,maxVox);
  //  for(int j = minVox; j <= maxVox; j++){
  //  VoxelGradient->SetBinContent(j,gradient[j]);
  //  }

 
	
	//Before the Fletcher-Reeves-Polak-Ribiere minimisations is performed, a few vectors / matricies need defining.
	
	
	
	
		//n - number of voxels
		//n = Ni*Nj*Nk
	
		//*func - This will be our Cost Function input is a vector (the vector which contains the scattering density of each voxel)
	
		//*dfunc - This will be our Gradient. it has two arguements, and no output. 
	
		//Start Algo.
	
		//p[] = voxel density map
	//ID of Voxel = I + J*Ni + K*Ni*Nj
	//n = size of map p[] - number of voxels = Nx*Ny*Nz
		//ftol = convergence tolerance
		//fret
	
	//Stuff for Algorithm from Recipes in C
	//	#include <math.h>
//	#include "nrutil.h"
//	#define ITMAX 200 //Max iterations
//	#define EPS 1.0e-10 //Small number to rectify the case of converging exactly to zero
//	#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);
	
	
	/*
	  Given a starting point, p[1..n] (in our case the vector Lambda), Fletcher-Reeves-Polak-Ribiere minimisation
	is performed on a function "func", using it's gradient as calculated by dfunc. The convergence tolerance on
	the function value is input as "ftol". Returned quantities are "p[1..n]" (the Value of Lambda for the minimum
	cost), "iter" (the number of iterations that were performed), and "fret" (the minimum value of the cost function).
	The routine linmin is called to perform line minimizations.	
	*/
	
	/*		void linmin(float_t p[], float_t xi[], int_t n, float_t *fret, float_t (*func)(float []));
		int_t j,its;
		float_t gg,gam,fp,dgg;
		float_t *g, *h, *xi;
		
		g=vector(1,n);
		h=vector(1,n);
		xi=vector(1,n);
		fp=(*func)(p); //Need to define Cost Function
		(*dfunc)(p,xi); //Need to define gradient of cost function
		for (j=1;j<=n;j++){
			g[j]=-xi[j];
			xi[j]=h[j]=g[j];
		}
		for(its=1;its<=ITMAX;its++){
			*iter=its;
			linmin(p, xi, n, fret, func); // Need to Not use this, and write my own linmin function.
			if(2.0*fabs(*fret-fp)<=ftol*(fabs(*fret)+fabs(fp)+EPS)){
				FREEALL
				return;
			}
			fp=*fret;
			(*dfunc)(p,xi);
			dgg = gg = 0.0;
			for (j=1;j<=n;j++){
				gg += g[j]*g[j];      
				 dgg += xi[j]*xi[j];  // Fletcher-Reeves.
				dgg += (xi[j]+g[j])*xi[j]; // Polak-Ribiere
			}
			if(gg=0.0){ //unlikely - If gradient is exactly zero, then we are already done
				FREEALL
				return;
			}
			gam = dgg/gg;
			for (j=1; j<=n ; j++){
				g[j] = -xi[j];
				xi[j] = h[j] = g[j] + gam * h[j];
			}
		}
		nrerror("Too many iterations in frprmn");
		
		*/
		
}
