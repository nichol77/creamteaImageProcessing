#include <TFile.h>
#include <TH3.h>
#include <TTree.h>
#include <TROOT.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include "DetectorDefs.hh"

double PI = 3.141592654;

void histoConvert(double xPos, double yPos, double zPos, double sideLength, double sphereRadius, char* inputHisto,  TFile* inFile);


void imageQuality(double nt, double nb, double Nt, double Nb);


int main() {
  
    char* inputHisto = "bkgndP";
  // char* inputHisto = "perfect";

  TFile* in = new TFile("./analysed.root","READ");
  //  TFile* in = new TFile("./images/finalNewData/7800_20000_r0.3_g2_7_2_0.08_0.07.root","READ");
  //TFile* in = new TFile("./testHistograms.root","READ");

    //std::cout << "Enter name of histogram: ";
    //std::cin >> inputHisto;
    //std::cout << "\n";

  histoConvert(SPHERE_X_M, SPHERE_Y_M, SPHERE_Z_M, SIDELENGTH, SPHERE_RADIUS_CM/100., inputHisto, in);

  return 0;

}


  void histoConvert(double xPos, double yPos, double zPos, double sideLength, double sphereRadius, char* inputHisto,  TFile* inFile) {
  
    TH3D* in = (TH3D*) inFile->Get(inputHisto); // loads the input histo

    double Rnorm = sphereRadius / sideLength;    // needs to be so that if Rnorm = 1/2, spheres stretched from one side of room to other, if centered
        
    // information to clone the axis to the out histogram
    int noBinsX = in->GetNbinsX();
    int noBinsY = in->GetNbinsY();
    int noBinsZ = in->GetNbinsZ();

    double Rx = Rnorm * noBinsX;
    double Ry = Rnorm * noBinsY;
    double Rz = Rnorm * noBinsZ;
  
    double nt = 0; // No. of bins inside target filled
    double nb = 0; // No. of bins outside target filled
    double Nt = (4 * PI * (pow(Rx,3)))/3; // Total No. of bins inside target
    double Nb = (noBinsX * noBinsY * noBinsZ) - Nt; // Total No. of bins outside target

    // normalise positions so that eg top corner is 1, 1, 1

    double xNormSphere = ( xPos + (sideLength / 2) ) / sideLength;
    double yNormSphere = ( yPos + (sideLength / 2) ) / sideLength;
    double zNormSphere = ( zPos + (sideLength / 2) ) / sideLength;
 


    double xAxMin = in->GetXaxis()->GetXmin();
    double xAxMax = in->GetXaxis()->GetXmax();

    double yAxMin = in->GetYaxis()->GetXmin();
    double yAxMax = in->GetYaxis()->GetXmax();

    double zAxMin = in->GetZaxis()->GetXmin();
    double zAxMax = in->GetZaxis()->GetXmax();

    // Sphere centre bin position

    double x0 = xNormSphere * noBinsX; 
    double y0 = yNormSphere * noBinsY;
    double z0 = zNormSphere * noBinsZ;  
  


    int ni = 0;


     for (int i = 0; i <= noBinsX; i++) {

       for (int j = 0; j<= noBinsY; j++) {

         for (int k = 0; k <= noBinsZ; k++) {
	    
    	   double binContent = in->GetBinContent(i,j,k);
          
         if(binContent > 0) {
           
          int xi = i;
          int yj = j;
          int zk = k;
	    
        ni++;
	
	double riSquared = pow((xi-x0),2) + pow((yj-y0),2) + pow((zk-z0),2);
	double ri = sqrt(riSquared);
	double RxSquared = pow(Rx,2);

        std::cout << xi << ", " << yj << ", " << zk << "    ni = " << ni << "     ri = " << ri;
        

	  
	if(riSquared < RxSquared) {
    	  nt++;
	  std::cout << ", nt = " << nt << " \n";
	  } else {
	  nb++;
	  std::cout << ", nb = " << nb << " \n";
	  } 

        // if(( xi < x0 + Rx ) && ( xi > x0 - Rx ) && ( yj < y0 + Ry ) && ( yj > y0 - Ry ) && ( zk < z0 + Rz ) && ( zk > z0 - Rz) ) {
	// nt++;
	//std::cout << ", nt = " << nt << " \n";
	// } else {
	//  nb++;
	//  std::cout << ", nb = " << nb << " \n";
	//  } 
	      
        }
        } 
      }
    }
    
  
  
  double iQ = ((nt / Nt) - (nb / Nb))*100;  
  
    std::cout << "no x bins = " << noBinsX << " no y bins = " << noBinsY << " no z bins = " << noBinsZ << "\n";
    
    double nEntries = in->GetEntries();

    std::cout << "no entries = " << nEntries << "\n";
      
    std::cout << "Sphere Radius = " << sphereRadius << ", sideLength = " << sideLength << ", Rnorm = sphere radius over sidelength = " << Rnorm << "\n";

    std::cout << "Rx = Rnorm * no bins x = " << Rx << ", Ry = " << Ry << ", Rz = " << Rz << "\n";

    std::cout << "Nt = (4/3) * PI * (sqrt(Rx^2 + Ry^2 + Rz^2))^3 = " << Nt << ", Nb = (no bins x * no bins y * no bins z) - Nt = " << Nb << "\n";

    std::cout << "xNormSphere = (xPos+(sideLength/2))/sideLength = " << xNormSphere << ", yNormSphere = " << yNormSphere << ", zNormSphere = " << zNormSphere << "\n";

    std::cout << "xAxMax, xAxMin = " << xAxMax << ", " << xAxMin << " y = " << yAxMax << ", " << yAxMin << ", z = " << zAxMax << ", " << zAxMin << "\n";

    std::cout << "x0 = xNormSphere * noBinsX = " << x0 << ", y0 = " << y0 << ", z0 = " << z0 << "\n";


  std::cout << "Image Quality for sphere positioned at: (" << xPos << ", " << yPos << ", " << zPos << "), of size " << sphereRadius << "m," << " with nt = " << nt << " nb = " << nb << ", Nt = " << Nt << ", and Nb = " << Nb << ", iQ = " << iQ << "% \n";
 
   
}
