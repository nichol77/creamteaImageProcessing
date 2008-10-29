//  *************************************************  //
//  *                                               *  //
//  * This creates two sets of histograms in outF:  *  //
//  *                                               *  //
//  * bkgnd (plain, N and P) which show the source  *  //
//  *       scattering histogram with a background  *  //
//  *       subtracted, N and P correspond to the   *  //
//  *       indvidual exess and dearth views of     *  //
//  *       this subtraction.                       *  //
//  *                                               *  //
//  * absorbed (plain, B and S for each slice)      *  //
//  *       which show the absorbed muons projected *  //
//  *       onto a plane in y with x and z          *  //
//  *       calculated with the gradient and cut    *  //
//  *       from StepThrough S and B correspond     *  //
//  *       to the individual source and background *  //
//  *       views                                   *  //
//  *                                               *  //
//  * Sam 18-09-08                                  *  //
//  *                                               *  //
//  *************************************************  //
  

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

// ******************** Variables ********************

char sourceF = "Source3DGS.dat"; // .dat file used to import test histo , extent and modifier values into loader

char outF    = "analysed.root";  // .root file name to save out put to

float length = SIDELENGTH*1000;  // length in mm of the side of the middle detector (layer 1) SIDELENGTH defined in global_var.hh

// **************** Preload functions ****************

// bkgndSubtraction subtracts a background histogram from the source histogram and then compares the value
// if it is larger than the error on the background file then a point is plotted.

void bkgndSubtraction    (char* inHistoName, char* outHistName, char* backgroundHisto, // in, out and background histogram names
			  int extent, float modifier,                                  // extent of the cube to test density over (cube side length (2*extent + 1)
			  TFile* inFile, TFile* outFile, TFile* backgroundFile);       // in, out and background files

// absorbedSubtraction takes <numberOfSlices> slices through y and at each slice calculates the projected position of
// the absorbed muons at that slice for both a background and a source file. These are then subtracted 
// and if the difference is larger than the background error a point is plotted

void absorbedSubtraction (char* bkgndTreeName, char* absTreeName, char* outHistName,   // tree and histo. names
			  int numberOfSlices, float sliceWidth ,                       // width and number of slices to take
			  int extent, float modifier,                                  // extent to test over and modifier
			  TFile* backgroundFile, TFile* inFile, TFile* outFile,        // in, out and background files
			  float startY = 0, int nBinsX = 100, int nBinsZ = 100,        // default values, starting position of y(mm) and number of bins to use in x and z
			  float xMax = length/2, float xMin = -length/2,               // range of x and z to use (mm)
			  float zMax = length/2, float zMin = -length/2);            

// Loader is used to strip variables from a .dat file so that 3DgradScan doesn't have to be 
// recompiled every time threshold modifiers, test extents or input/background files are changed
// the values are read from sourceF

int loader               (int* extent, int* extent2,                                   // two integer extents
			  float* modifier1, float* modifier2,                          // modifiers
			  char* inputFile, char* bkgndFile);                           // input and background files

// ********************** Main **********************

int main ()
{
  char inputFile[256]; // variables to take values from loader
  char bkgndFile[256];
  int extent = 1;
  int extent2 = 2;
  float modifier1 = 0;
  float modifier2 = 0;

  if ( !loader(&extent,&extent2, &modifier1,&modifier2, inputFile, bkgndFile) ) // if the loading fails it returns 0 
    {
      return 0;
    }
   
  TFile* out = new TFile(outF,"RECREATE");

  TFile* in = new TFile(inputFile,"READ");

  TFile* background = new TFile(bkgndFile,"READ");

  //******************

  bkgndSubtraction ("PCAh","bkgnd","PCAh",extent, modifier1, in, out, background);

  absorbedSubtraction ("Absorbed","Absorbed","absorbed" ,9,100 ,extent2,modifier2 ,background, in, out,500);

  out->Write();

  in->Close();
  out->Close();
  background->Close();

  std::cout<<std::endl;
  std::cout<<"done"<<std::endl;

  return 0;
}


// ************* Function definitions **************


void bkgndSubtraction  (char* inHistoName, char* outHistName,char* backgroundHisto, int extent, float modifier, TFile* inFile, TFile* outFile, TFile* backgroundFile)
{
  TH3D* in = (TH3D*) inFile->Get(inHistoName); // loads the input histo

  TH3D* bkgnd = (TH3D*) backgroundFile->Get(backgroundHisto); // loads the background histo

  // information to clone the axis to the out histogram
  int nBinsX = in->GetNbinsX();
  int nBinsY = in->GetNbinsY();
  int nBinsZ = in->GetNbinsZ();

  float xAxisMin = in->GetXaxis()->GetXmin();
  float xAxisMax = in->GetXaxis()->GetXmax();

  float yAxisMin = in->GetYaxis()->GetXmin();
  float yAxisMax = in->GetYaxis()->GetXmax();

  float zAxisMin = in->GetZaxis()->GetXmin();
  float zAxisMax = in->GetZaxis()->GetXmax();

  //******************** create output histograms

  char outHistNameN [25];
  char outHistNameP [25];

  strcpy (outHistNameN, outHistName);
  strcpy (outHistNameP, outHistName);

  strcat (outHistNameN, "N");
  strcat (outHistNameP, "P");

  TH3D* out  = new TH3D(outHistName,  outHistName,  nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax); // total output
  TH3D* outN = new TH3D(outHistNameN, outHistNameN, nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax); // where source < background
  TH3D* outP = new TH3D(outHistNameP, outHistNameP, nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax); // where source > background

  // calculations for tests

  double nEntriesS = in->GetEntries();    //number of source entries
  double nEntriesB = bkgnd->GetEntries();  // number of background entries

  double vol = std::pow ( (float)(2*extent + 1), 3) ; // volume of the test region

  double volTotal = (nBinsX * nBinsY * nBinsZ);

  double threshold = modifier * std::sqrt ((nEntriesB+nEntriesS)/volTotal);  // <modifier> number of std deviations on the density of the test region 

  double scalingConst = nEntriesS / nEntriesB; // used to account for different numbers of muons

  //******************** run the test on all points

  for (int i = 0; i <= nBinsX; i++)
    {
      for (int j = 0; j<= nBinsY; j++)
	{
	  for (int k = 0; k <= nBinsZ; k++)
	    {

	      double sourceVal = 0;
	      double bkgndVal  = 0;

	      for (int x = -extent; x <= extent; x++)
		{
		  for (int y = -extent; y <= extent; y++)
		    {
		      for (int z = -extent; z <= extent; z++)
			{
			  sourceVal += in->GetBinContent ((i+x),(j+y),(k+z));
			  bkgndVal  += bkgnd->GetBinContent ((i+x),(j+y),(k+z));
			} 
		    }		  
		}
	      sourceVal = sourceVal /(vol*scalingConst);  // density to be tested
	      bkgndVal  = (bkgndVal*scalingConst)  /vol;  // expected density

	      if( std::fabs (sourceVal - bkgndVal) > threshold )// if there is a significant difference in densities 
		{
		  out->SetBinContent (i,j,k,1);
		  
		  if (sourceVal > bkgndVal)
		    {
		     outP->SetBinContent (i,j,k,1);
		    }
		  else if (sourceVal < bkgndVal)
		    {
		      outN->SetBinContent (i,j,k,1);
		    }		    
		}

	    }
	}
      if (!(i%10)) std::cout<<i<< " of " << nBinsX <<std::endl;
    }
  outFile->Add(out);
  outFile->Add(outP);
  outFile->Add(outN);

  std::cout<<(threshold)<<" threshold "<<std::endl;
}

void absorbedSubtraction (char* bkgndTreeName, char* absTreeName, char* outHistName, int numberOfSlices,float sliceWidth ,int extent, float modifier,TFile* backgroundFile, TFile* inFile, TFile* outFile,float startY, int nBinsX, int nBinsZ,float xMax, float xMin,float zMax, float zMin) 
{
  // ******************** absorbed muon tree
  double xCut, xGrad, zCut, zGrad;

  TTree* abs = (TTree*) inFile->Get(absTreeName);
  
  abs->SetBranchAddress("xCut",&xCut);
  abs->SetBranchAddress("zCut",&zCut);

  abs->SetBranchAddress("xGrad",&xGrad);
  abs->SetBranchAddress("zGrad",&zGrad);

  int nEntriesS = abs->GetEntries();

  // ******************** background absorbed muon tree
  double xCutB, xGradB, zCutB, zGradB;

  TTree* absB = (TTree*) backgroundFile->Get(bkgndTreeName);
  
  absB->SetBranchAddress("xCut",&xCutB);
  absB->SetBranchAddress("zCut",&zCutB);

  absB->SetBranchAddress("xGrad",&xGradB);
  absB->SetBranchAddress("zGrad",&zGradB);

  int nEntriesB = absB->GetEntries();

  // ******************** calculate the threshold and similar

  float vol = std::pow ((float) (2*extent+1),3);

  float threshold = modifier; 

  float scalingConstant = ((float) nEntriesS) / ((float)nEntriesB); 
  
  // ******************** run for each slice 

  for (int j = 0; j < numberOfSlices; j++)
    {
      double y = (j*sliceWidth + startY); // the 'real' value of y
      char name  [25];
      char nameS [25];
      char nameB [25];
      char yVal  [5];
	  
      sprintf (yVal, "%d",(int)y);

      strcpy (name,outHistName);
      strcat (name, yVal);
      
      strcpy (nameB, outHistName);
      strcat (nameB, "B");
      strcat (nameB, yVal);
      
      strcpy (nameS, outHistName);
      strcat (nameS, "S");
      strcat (nameS, yVal);

      TH2D* out  = new TH2D (name, name,  nBinsX,xMin,xMax, nBinsZ,zMin,zMax);
      TH2D* outS = new TH2D (nameS,nameS, nBinsX,xMin,xMax, nBinsZ,zMin,zMax);
      TH2D* outB = new TH2D (nameB,nameB, nBinsX,xMin,xMax, nBinsZ,zMin,zMax);
      
      // ******************** create the source histo
      for (int entry = 0; entry < nEntriesS; entry++)
	{
	  abs->GetEntry(entry);
	  
	  if(xGrad == 0 || zGrad == 0) continue; // if there is null gradient ignore and continue      

	  double x = (y - xCut)/ xGrad;  // calculate x and z positions using the y value (x = (y-cut)/grad)
	  double z = (y - zCut)/ zGrad;  
	  
	  outS->Fill (x,z,1);
	}
  
      outFile->Add(outS);
      // ******************** create background histo 

      for (int entryB = 0; entryB < nEntriesB; entryB++)
	{
	  absB->GetEntry(entryB);
	  
	  if(xGradB == 0 || zGradB == 0) continue; // if there is null gradient ignore and continue      

	  double x = (y - xCutB)/ xGradB;  // calculate x and z positions using the y value (x = (y-cut)/grad)
	  double z = (y - zCutB)/ zGradB;
	  
	  outB->Fill (x,z,1);
	}
       outFile->Add(outB);

      // ******************** test the source density against background density

      for (int i =0; i < nBinsX; i++)
	{
	  for (int k = 0; k < nBinsZ; k++)
	    {
	      float bkgnd = 0;
	      float test  = 0;

	      for (int io = -extent; io <= extent; io++)
		{
		  for (int ko = -extent; ko <= extent; ko++)
		    {
		      test  += outS->GetBinContent ((i+io),(k+ko));
		      bkgnd += outB->GetBinContent ((i+io),(k+ko));
		    }
		}

	      test  = (test/vol) /scalingConstant; // ensures that the densities are comparable
	      bkgnd = (bkgnd/vol)*scalingConstant;

	      if (std::fabs (test-bkgnd) > threshold)
		{
		  out->SetBinContent (i,k,1);  
		}    
	    }
	  if (!(i%10)) std::cout<<(j*nBinsX)+ i<<" of "<<(numberOfSlices*nBinsX)<<std::endl;
	}
      outFile->Add(out);
    }
  std::cout<<vol<<" "<<threshold<<std::endl;
}

int loader (int* extent, int* extent2, float* modifier1, float* modifier2, char* inputFile, char* bkgndFile)
{
  char v_inputFile[256];
      char v_bkgndFile[256];

  // default values
      int  v_extent = 1;
      int v_extent2 = 1;
      float v_modifier1 = 0;
      float v_modifier2 =0;

      ifstream in;
      in.open (sourceF);
      in>>v_extent>>v_extent2>>v_modifier1>>v_modifier2>>v_inputFile>>v_bkgndFile;
      in.close();
  
      if(!v_inputFile[0])
  {
  std::cout<<"no file loaded. Quiting"<<std::endl;
      
          return 0;      
    }
      else if(!v_bkgndFile[0])
  {
  std::cout<<"no histogram specified. Quitting"<<std::endl;

          return 0;
    }
      else
  {      
  *extent = v_extent;
          *extent2 = v_extent2;
          *modifier1 = v_modifier1;
          *modifier2 = v_modifier2;

          std::strcpy(inputFile,v_inputFile);
          std::strcpy(bkgndFile, v_bkgndFile);

          return 1;
    }
}
