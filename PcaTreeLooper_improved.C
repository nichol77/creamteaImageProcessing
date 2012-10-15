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
clock_t start = clock();

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
Long64_t nentries = 5000; //For testing, so we don't take all the muons


int Nx = 10; //Number of voxels in x direction
int Ny = 10; //Number of voxels in y direction
int Nz = 10; //Number of voxels in z direction

//Recipricol Voxel dimensions (multiplication faster)
double Vx = Nx/(MaxX-MinX);
double Vy = Ny/(MaxY-MinY);
double Vz = Nz/(MaxZ-MinZ);

//Populate Voxel scattering densities
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
vector<double> S ((int)nentries);


//...................................................................................... -> Make Lij Matrix <- .............	

//Define Lij Matrix kind of thing
//Reference -> http://groups.google.com/group/comp.lang.c++.moderated/browse_thread/thread/bb499f84b93c4343?pli=1
typedef map<int, map<int, int> > MapMap; 
MapMap L;
// N * L[i][j] <- Et Voila. Lij.


	
//Initialise some variables so we don't get loads of errors.
double xEntry, yEntry, xExit, yExit, TopVx, TopVy, TopVz, BottomVx, BottomVy, BottomVz, RLength, CurrentX, CurrentY, CurrentZ;
int VoxelI, VoxelJ, VoxelK, ThisI, ThisJ, ThisK, i, ThisVoxelID;
Long64_t nbytes = 0;
Long64_t nb = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   nbytes += nb;
     
 
  S[(int) jentry] = thetaTrue;

	for(int part=0; part < 2; part++){
//Current (i.e. first) Voxel
int CurrentX = floor((xzGradTrue[0] * MaxZ + xzCutTrue[0])/Vx);
int CurrentY = floor((yzGradTrue[0] * MaxZ + yzCutTrue[0])/Vy);
int CurrentZ = Nz-1;
int vID = CurrentX + Nx * CurentY + Nx * Ny * CurrentZ;

//End Voxel
if(part == 0){
int EndX = floor(xPosTrue/Vx);
int EndY = floor(yPosTrue/Vy);
int EndZ = floor(zPosTrue/Vz);
int EndVoxelID = PcaX + Nx * PcaY + Nx * Ny * PcaZ;
}else if(part == 1){
int EndX = floor(xPosTrue/Vx);
int EndY = floor(yPosTrue/Vy);
int EndZ = floor(zPosTrue/Vz);
int EndVoxelID = PcaX + Nx * PcaY + Nx * Ny * PcaZ;
}

//Entry & Exit coords. Leave Exit blank so it's done in the while loop.
double xEntry = xzGradTrue[part] * MaxZ + xzCutTrue[part];
double yEntry = yzGradTrue[part] * MaxZ + yzCutTrue[part];
double zEntry = MaxZ;
double xExit, yExit, zExit;

int go = 1;
while(go != 0){
if(xzGradTrue[part] > 0){zExit = zEntry - (xEntry - CurrentX * Vx)/xzGradTrue[part];}
if(xzGradTrue[part] < 0){zExit = zEntry - (xEntry - (CurrentX+1) * Vx)/xzGradTrue[part];}
if(xzGradTrue[part] = 0){zExit = zEntry - Vz;}
if(yzGradTrue[part] > 0){zExit = zEntry - (yEntry - CurrentY * Vy)/yzGradTrue[part];}
if(yzGradTrue[part] < 0){zExit = zEntry - (yEntry - (CurrentY+1) * Vx)/yzGradTrue[part];}
if(yzGradTrue[part] = 0){zExit = zEntry – Vz;} //!!

xExit = xzGradTrue[part] * zExit + xzCutTrue[part];
yExit = yzGradTrue[part] * zExit + yzCutTrue[part];
//Length
L[i][vID] = pow(pow(xEntry-xExit,2) + pow(yEntry-yExit,2) + pow(zEntry-zExit,2),0.5);

//Set the new Entry Points to be the old Exit Points
xEntry = xExit;
yEntry = yExit;
zEntry = zExit;

//Which Voxel is it in now?
CurrentX = round(xEntry/Vx);
CurrentY = round(yEntry/Vy);
CurrentZ = round(zEntry/Vz);
vID = CurrentX + Nx * CurentY + Nx * Ny * CurrentZ;

if(EndVoxelID == vID){go = 0;}
}
}  


  }
}

 
//...................................................................................... -> Lij Matrix Made <- .............




//Gradient Vector
 double temp = 0;
 int minVox = 0;
 int maxVox = Nx*Ny*Nz;
 TH1F *MuonsPerVoxel = new TH1F("MuonsPerVoxel","Muons Per Voxel",maxVox-minVox,minVox,maxVox);
 TH1F *TotalLengthPerVoxel = new TH1F("TotalLengthPerVoxel","Total Pathlength Per Voxel",maxVox-minVox,minVox,maxVox);


 map<int,double> Sigma;
 for (int iter1 = 0; iter1 < nentries; ++iter1) {
   for(std::map<int,int>::iterator iter2 = L[iter1].begin(); iter2 != L[iter1].end(); ++iter2){
     Sigma[iter1] += iter2->second * N * (Lambda[iter2->first]);
   }
 }
 for(int iter1 = 0; iter1 < (int)nentries; ++iter1){
   for(std::map<int,int>::iterator iter2 = L[iter1].begin(); iter2 != L[iter1].end(); ++iter2){
     int VoxelID = iter2->first;
     double PathLength = iter2->second * N;
     if(VoxelID >= 0 && VoxelID < Nx*Ny*Nz && PathLength > 0){
       temp = PathLength * (1 - S[iter1]/Sigma[iter1]);
       if(temp != temp){temp = 0; cout << VoxelID << endl;}
       if(fabs(temp)==DBL_MAX|| temp > 1.7 * pow(10,250)){temp = 1.7 * pow(10,250);}
       gradient[VoxelID]+=temp;
     }
     if(VoxelID >= minVox && VoxelID <= maxVox){
       MuonsPerVoxel->Fill(VoxelID,1);
       TotalLengthPerVoxel->Fill(VoxelID,PathLength);
     }
   }
  }
 double b = 0;
 double Psi = 0;
 map<int,double> Lambda2 = Lambda;
 TH1F *CostFunctionAlpha = new TH1F("CostFunctionAlpha","Psi(Lambda - alpha*Gradient)",400,-200,200);
 for(int a = -200; a < 200 ;a++){
   Psi = 0;
   b = 0.0000000001 * (double) a;

   for(int v = 0; v < Nx*Ny*Nz; ++v){
     Lambda2[v] = Lambda[v] + b*gradient[v];
     if(Lambda2[v] < 0.000003067 ){
       Lambda2[v] = 0.000003067;
     }     
   }

   for (int iter1 = 0; iter1 < nentries; ++iter1) {
     for(std::map<int,int>::iterator iter2 = L[iter1].begin(); iter2 != L[iter1].end(); ++iter2){
       Sigma[iter1] += iter2->second * N * (Lambda2[iter2->first]);
     }
   }
   
   for(int iter1 = 0; iter1 < (int)nentries; ++iter1){
     if(Sigma[iter1] > 0){
       Psi += S[iter1]*S[iter1]/Sigma[iter1] + log(Sigma[iter1]);
     }
   }
   //cout << Psi << endl;
   CostFunctionAlpha->Fill(a,Psi);
   Sigma.clear();
 }
 
 







  TH1F *VoxelGradient = new TH1F("VoxelGradient","Gradient Of Each Voxel)",maxVox-minVox,minVox,maxVox);
  for(int i = minVox; i < maxVox; i++){
    VoxelGradient->Fill(i,gradient[i]);
  }

  TFile *fpOu1 = new TFile("MuonPath1.root","RECREATE");
  fpOu1->Write();
  TCanvas *canProj = new TCanvas();
  canProj->Divide(1,3);
  canProj->cd(1);
  MuonsPerVoxel->Draw();
  MuonsPerVoxel->SetStats(0);
  MuonsPerVoxel->GetXaxis()->SetTitle("Voxel ID");
  canProj->cd(2);
  VoxelGradient->Draw();
  VoxelGradient->SetStats(0);
  VoxelGradient->GetXaxis()->SetTitle("Voxel ID");
  canProj->cd(3);
  CostFunctionAlpha->Draw();
  CostFunctionAlpha->SetStats(0);
  CostFunctionAlpha->GetXaxis()->SetTitle("a times a million or something");
 


  //  TH1F *VoxelGradient = new TH1F("VoxelGradient","Gradient at each voxel",maxVox,0,maxVox);
  //  for(int j = minVox; j <= maxVox; j++){
  //  VoxelGradient->SetBinContent(j,gradient[j]);
  //  }

 
 printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
	
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
