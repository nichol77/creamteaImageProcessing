#include <TFile.h>
#include <TH3.h>
#include <TTree.h>
#include <TROOT.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include "../global_vars.hh"

double PI = 3.141592654;

void testHistoCreate(double xPos, double yPos, double zPos, double sideLength, double sphereRadius, char* inputHisto, TFile* inFile, TFile* outFile);

int main(){

  char* inputHisto = "bkgnd";
  
  TFile* in = new TFile("./analysed.root","READ");

// Open output file (apparently needs to be done before booking)

  TFile* out = new TFile("./testHistograms.root","RECREATE");

  testHistoCreate(XPOS, YPOS, ZPOS, SIDELENGTH, SPHERERADIUS, inputHisto, in, out);

  out->Write();

  in->Close();
  out->Close();

  return 0;

}

// Book histograms

void testHistoCreate(double xPos, double yPos, double zPos, double sideLength, double sphereRadius, char* inputHisto, TFile* inFile, TFile* outFile) {

  TH3D* in = (TH3D*) inFile->Get(inputHisto); // loads the input histo

  double Rnorm = sphereRadius / sideLength;    // needs to be so that if Rnorm = 1/2, spheres stretched from one side of room to other, if centered
        
  // information to clone the axis to the out histogram
  
  int noBinsX = in->GetNbinsX();
  int noBinsY = in->GetNbinsY();
  int noBinsZ = in->GetNbinsZ();

  double Rx = Rnorm * noBinsX;
  //double Ry = Rnorm * noBinsY;
  //double Rz = Rnorm * noBinsZ;
  
  //double nt = 0; // No. of bins inside target filled
  //double nb = 0; // No. of bins outside target filled
  //double Nt = (4 * PI * (pow(Rx,3)))/3; // Total No. of bins inside target
  //double Nb = (noBinsX * noBinsY * noBinsZ) - Nt; // Total No. of bins outside target

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
 
  TH3D* perfectHisto = new TH3D("perfect", "perfect image", noBinsX,xAxMin,xAxMax, noBinsY,yAxMin,yAxMax, noBinsZ,zAxMin,zAxMax); 

//  TH3D* badHisto = new TH3D("bad", "bad image",  noBinsX,xAxMin,xAxMax, noBinsY,yAxMin,yAxMax, noBinsZ,zAxMin,zAxMax); 

//  TH3D* perfectBadHisto = new TH3D("perfectBad", "perfect bad image", noBinsX,xAxMin,xAxMax, noBinsY,yAxMin,yAxMax, noBinsZ,zAxMin,zAxMax); 

  int ni = 0;
  
  for (int i = 0; i <= noBinsX; i++) {
    
    for (int j = 0; j <= noBinsY; j++) {
     
      for (int k = 0; k <= noBinsZ; k++) {
	    
  	int xi = i;
        int yj = j;
        int zk = k;
	
	ni++;
	
        double riSquared = pow((xi-x0),2) + pow((yj-y0),2) + pow((zk-z0),2);
        double ri = sqrt(riSquared);
        double RxSquared = pow(Rx,2);
      	
        if(riSquared < RxSquared) { 

	  std::cout << xi << ", " << yj << ", " << zk << "    ni = " << ni << "     ri = " << ri << "\n";

          perfectHisto->SetBinContent (i,j,k,1);

        } 

	//else if (riSquared > RxSquared) {
	 
	//perfectBadHisto->SetBinContent (i,j,k,1);

  	//}

      }

    }
  
  }
  
  outFile->Add(perfectHisto);
//  outFile->Add(badHisto);
//  outFile->Add(perfectBadHisto);

 

}
