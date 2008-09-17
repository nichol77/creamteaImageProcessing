/*  
 * backgroundTest_2 works the rest is still needing development
 *
 *
 * The idea of this is that unlike the edge detection program
 * that works by finding areas with a higher than expected
 * density of scattering points (it should be renamed)
 * this will find points based on the 'gradient' of scattering 
 * points surrounding it. 
 * 
 * ie a point with 10 hits on one side and 0 on the other will
 * have a high gradient.
 * 
 * this scan will pass along each axis individually, first x then y 
 * and finally z.  
 * 
 * At a later point a method of 'smoothing' may be used initially to
 * reduce background noise.
 *
 * NB 'v_' prefix indicates values taken from input streams that
 *         are therefore non const - therefore they get passed to
 *         new variables.
 * 
 * Sam 5 - 08 - 08
 */

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

// file used to import test histo , extent and modifier values
#define SOURCE_F "Source3DGS.dat"

// file to load background data and the name of the relevant histogram

//using namespace std;

float length = SIDELENGTH*1000;

void backgroundTest_2(char* inHistoName, char* outHistName, char* backgroundHisto, int extent, float modifier, TFile* inFile, TFile* outFile, TFile* backgroundFile);

void backgroundTest  (char* inHistoName, char* outHistName, char* backgroundHisto, float modifier, TFile* inFile, TFile* outFile, TFile* backgroundFile);

void averager        (char* inHistoName, char* outHistName, int extent, TFile* inFile, TFile* outFile);

void densityTest     (char* inHistoName, char* outHistName, int extent, float modifier, TFile* inFile, TFile* outFile);

void densityGradient (char* inHistoName, char* outHistName, int extent, float modifier, TFile* inFile, TFile* outFile);

void linearGradient  (char* inHistoName, char* outHistName, int extent, float modifier, TFile* inFile, TFile* outFile);

void absorbed        (char* inHistoName, char* absTreeName, char* outHistName,int extent, float modifier, TFile* inFile, TFile* absorbedFile, TFile* outFile);

void absorbedSlices  (char* inHistoName, char* absTreeName, char* outHistName, int numberOfSlices,float modifier,TFile* inFile, TFile* absorbedFile, TFile* outFile);

void absorbedSubtraction (char* bkgndTreeName, char* absTreeName, char* outHistName, int numberOfSlices,float sliceWidth , int extent, float modifier, TFile* backgroundFile, TFile* absorbedFile, TFile* outFile, float startY = 0, int nBinsX = 100, int nBinsZ = 100, float xMax = length/2, float xMin = -length/2, float zMax = length/2, float zMin = -length/2);

int loader           (int* extent, int* sliceNumber, float* dtModifier, float* gradModifier, char* inputFile, char* histoName);

double permute       (TH3D* in, int index, int x, int y, int z, char option);

//********************

int main ()
{
  char inputFile[256];
  char histoName[256];
  int extent = 1;
  int sliceNumber = 2;
  float dtModifier = 0;
  float gradModifier = 0;

  if ( !loader(&extent,&sliceNumber, &dtModifier,&gradModifier, inputFile, histoName) ) // if the loading fails it returns 0
    {
      return 0;
    }
   
  TFile* out = new TFile("analysed.root","RECREATE");

  TFile* in = new TFile(inputFile,"READ");

  TFile* background = new TFile(histoName,"READ");

  //******************

  backgroundTest_2 ("PCAh","bkgnd","PCAh",extent, dtModifier, in, out, background);

  absorbedSubtraction ("Absorbed","Absorbed","absorbed" ,9,100 ,sliceNumber,gradModifier ,background, in, out,500);

  out->Write();

  in->Close();
  out->Close();
  background->Close();

  std::cout<<std::endl;
  std::cout<<"done"<<std::endl;

  return 0;
}


//********************

// tests the deviation of the source files from the background file by comparing their densities 
void backgroundTest_2  (char* inHistoName, char* outHistName,char* backgroundHisto, int extent, float modifier, TFile* inFile, TFile* outFile, TFile* backgroundFile)
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

// subtracts the background value of each point from the input-Histogram, background values taken from an 'empty' run of the area
void backgroundTest  (char* inHistoName, char* outHistName,char* backgroundHisto, float modifier, TFile* inFile, TFile* outFile, TFile* backgroundFile)
{

  TH3D* in = (TH3D*) inFile->Get(inHistoName); // loads the input histo

      TH3D* bkgnd = (TH3D*) backgroundFile->Get(backgroundHisto); // loads the background histo

      int nBinsX = in->GetNbinsX();
      int nBinsY = in->GetNbinsY();
      int nBinsZ = in->GetNbinsZ();

      float xAxisMin = in->GetXaxis()->GetXmin();
      float xAxisMax = in->GetXaxis()->GetXmax();

      float yAxisMin = in->GetYaxis()->GetXmin();
      float yAxisMax = in->GetYaxis()->GetXmax();

      float zAxisMin = in->GetZaxis()->GetXmin();
      float zAxisMax = in->GetZaxis()->GetXmax();

  //double nEntriesS = in->GetEntries();    //number of source entries
  //double nEntriesB = bkgnd->GetEntries();  // number of background entries

      double threshold = modifier; //*std::sqrt ((nEntriesB + nEntriesS)/vol); // modifier is number of standard deviations from zero to consider significant

      TH3D* out = new TH3D(outHistName, outHistName, nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax);

      for (int i = 0; i <= nBinsX; i++)
  {
  for (int j = 0; j<= nBinsY; j++)
  {
  for (int k = 0; k <= nBinsZ; k++)
  {
  double sourceVal = in->GetBinContent (i,j,k);
	          double bkgndVal = bkgnd->GetBinContent (i,j,k);

	      //std::cout<<sourceVal<<" "<< bkgndVal<<std::endl;
	      
	          if (std::fabs (sourceVal - bkgndVal)> threshold ) // if the difference is greater than <modifier> number of std deviations consider it significant
  {
  out->SetBinContent (i,j,k,(sourceVal - bkgndVal) );
		  //std::cout<<"here i am" <<std::endl;
		}
	    }
	}
          std::cout<<i<< " of " << nBinsX <<std::endl;
    }

      outFile->Add(out);

      std::cout<<threshold<<" threshold"<<std::endl;
}

// uses a mean average to place a value at each point
void averager (char* inHistoName, char* outHistName, int extent, TFile* inFile, TFile* outFile)
{
  TH3D* in = (TH3D*) inFile->Get(inHistoName);

      int nBinsX = in->GetNbinsX();
      int nBinsY = in->GetNbinsY();
      int nBinsZ = in->GetNbinsZ();

      float xAxisMin = in->GetXaxis()->GetXmin();
      float xAxisMax = in->GetXaxis()->GetXmax();

      float yAxisMin = in->GetYaxis()->GetXmin();
      float yAxisMax = in->GetYaxis()->GetXmax();

      float zAxisMin = in->GetZaxis()->GetXmin();
      float zAxisMax = in->GetZaxis()->GetXmax();

      double vol = std::pow( ((double) (extent*2 + 1)) ,3); 

      TH3D* out = new TH3D(outHistName, outHistName, nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax);

      for (int i = 0; i <= nBinsX; i++)
  {
  for (int j = 0; j<= nBinsY; j++)
  {
  for (int k = 0; k <= nBinsZ; k++)
  {
	      
  double test =0;

	          for (int x = -extent; x <= extent; x++)
  {
  for (int y = -extent; y <= extent; y++)
  {
  for (int z = -extent; z <= extent; z++)
  {
  test += in->GetBinContent ( (i+x),(j+y),(k+z) );
			} 
		    }		  
		}
	          test = test/vol;
	      
	          if (!out->GetBinContent(i,j,k)) 
  {
  out->SetBinContent (i,j,k,(int) test);
		}
	    }
	}
          std::cout<<i<< " of " << nBinsX <<std::endl;
    }

      outFile->Add(out);
}

// tests the density of a region side length 2*extent+1 against the expected uniform density of the same region
void densityTest (char* inHistoName, char* outHistName, int extent, float modifier, TFile* inFile, TFile* outFile)
{
  TH3D* in = (TH3D*) inFile->Get(inHistoName);

      int nBinsX = in->GetNbinsX();
      int nBinsY = in->GetNbinsY();
      int nBinsZ = in->GetNbinsZ();

      float xAxisMin = in->GetXaxis()->GetXmin();
      float xAxisMax = in->GetXaxis()->GetXmax();

      float yAxisMin = in->GetYaxis()->GetXmin();
      float yAxisMax = in->GetYaxis()->GetXmax();

      float zAxisMin = in->GetZaxis()->GetXmin();
      float zAxisMax = in->GetZaxis()->GetXmax();

      double nEntries = in->GetEntries();

      TH3D* out = new TH3D(outHistName, outHistName, nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax);

      double threshold = ((extent*2 + 1)*(extent*2 + 1)*(extent*2 + 1)) * nEntries / (nBinsX * nBinsY * nBinsZ); // predicted density for the test region

      double error = modifier; // TODO set this so its the std dev on the number of PCA within the region

      std::cout<<  threshold <<std::endl;

      for (int i = 0; i <= nBinsX; i++)
  {
  for (int j = 0; j<= nBinsY; j++)
  {
  for (int k = 0; k <= nBinsZ; k++)
  {
	      
  double test =0;

	          for (int x = -extent; x <= extent; x++)
  {
  for (int y = -extent; y <= extent; y++)
  {
  for (int z = -extent; z <= extent; z++)
  {
  test += in->GetBinContent ( (i+x),(j+y),(k+z) );
			} 
		    }		  
		}


	          if ( (test < (threshold - error)) || ((threshold + error) < test) )
  {
  out->SetBinContent(i,j,k,1);
		  // std::cout<<test<<std::endl;
		  // std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		}
	    }
	}
          std::cout<<i<< " of " << nBinsX <<std::endl;
    }

      outFile->Add(out);

      std::cout<<threshold<<" threshold and error: "<<error<<std::endl;
  
  //outFile->Write();
}

// gradient edge test that uses the gradient across faces of the region
void densityGradient (char* inHistoName, char* outHistName, int extent, float modifier, TFile* inFile, TFile* outFile) 
{
  TH3D* in = (TH3D*) inFile->Get(inHistoName);

  // get the dimensions and range of the input histogram and make the output the same.

      int nBinsX = in->GetNbinsX();
      int nBinsY = in->GetNbinsY();
      int nBinsZ = in->GetNbinsZ();

      float xAxisMin = in->GetXaxis()->GetXmin();
      float xAxisMax = in->GetXaxis()->GetXmax();

      float yAxisMin = in->GetYaxis()->GetXmin();
      float yAxisMax = in->GetYaxis()->GetXmax();

      float zAxisMin = in->GetZaxis()->GetXmin();
      float zAxisMax = in->GetZaxis()->GetXmax();

      TH3D* out = new TH3D(outHistName, outHistName, nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax);

  // these arrays should allow cycling through the histogram in all 3 directions the <ijk>OP correspond to the 
  // offset (xyz) co-ordinates and the <ijk>P values are the source point co-ordinates

  // threshold calculates as (mean error on the number of hits per cell) * modifier 
      float threshold = modifier;// * std::sqrt(nEntries) / (nBinsZ * nBinsY * nBinsX);


      for(int perm = 0; perm < 3; perm++) // loop through all permutations
  {
  // ********************

  for (int i = 0; i < nBinsX; i++) // loop through the entire histogram
  {
  for (int j = 0; j < nBinsY; j++)
  {
  for (int k = 0; k < nBinsZ; k++)
  {
  // ********************

  float test = 0;
		  
		      for(int x = -extent; x <= extent; x++) // loop through each point surronding the current
  {
  for(int y = -extent; y <= extent; y++)
  {
  for(int z = -extent; z <= extent; z++)
  {
  //test += weight[std::abs(x)][std::abs(y)][std::abs(z)]* permute(out,perm,(i+x),(j+y),(k+z),'G');
  if(z<0)
  {
  test += permute(in,perm,(i+x),(j+y),(k+z),'G');
				}
			          else if(z>0)
  {
  test -= permute(in,perm,(i+x),(j+y),(k+z),'G');
				}
			          else
  {
  continue;
				}
			    } 
			}		  
		    }

		      if ( (std::fabs(test) > threshold) && !(permute(out,perm,i,j,k,'G')) )
  {
  permute(out,perm,i,j,k,'F');
		    }
		  // ********************
		}
	    }

	      std::cout<< i <<" of "<<nBinsX<<" on permutation set "<<perm<<std::endl;
	}
      // ********************
    }

      outFile->Add(out);
}

// gradient edge test that uses the gradient along a line through the source
void linearGradient (char* inHistoName, char* outHistName, int extent, float modifier, TFile* inFile, TFile* outFile) 
{
  TH3D* in = (TH3D*) inFile->Get(inHistoName);

      int nBinsX = in->GetNbinsX();
      int nBinsY = in->GetNbinsY();
      int nBinsZ = in->GetNbinsZ();

      float xAxisMin = in->GetXaxis()->GetXmin();
      float xAxisMax = in->GetXaxis()->GetXmax();

      float yAxisMin = in->GetYaxis()->GetXmin();
      float yAxisMax = in->GetYaxis()->GetXmax();

      float zAxisMin = in->GetZaxis()->GetXmin();
      float zAxisMax = in->GetZaxis()->GetXmax();

      float threshold = modifier;// * std::sqrt(nEntries) / (nBinsZ * nBinsY * nBinsX);

      TH3D* out = new TH3D(outHistName, outHistName, nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax);

      for (int i = 0; i <= nBinsX; i++)
  {
  for (int j = 0; j<= nBinsY; j++)
  {
  for (int k = 0; k <= nBinsZ; k++)
  {

  double test[3] = {0,0,0};
	          for (int r = -extent; r <= extent; r++)
  {
  if (r==0)continue; // ignore the source position
		  
		      for (int perm = 0; perm < 3; perm ++)
  {
  test[0] += (1/r) * permute (in, 0, i, j, (k+r), 'G' );
	 	    }
		  //   std::cout<<r<<" "<<test[0];
		}
	      //std::cout<<std::endl;
	      // std::cout<<test[0]<<std::endl;

	          if ( (test[0]>threshold) ||  (test[1]>threshold) ||  (test[2]>threshold)  )
  {
  permute(out,0,i,j,k,'F');
		}
	      
	    }
	}
          std::cout<<i<< " of " << nBinsX <<std::endl;
    }
      outFile->Add(out);
   
  //outFile->Write();
}

// tests whether it is likely a muon has been absorbed
void absorbed  (char* inHistoName, char* absTreeName, char* outHistName,int extent, float modifier, TFile* inFile, TFile* absorbedFile, TFile* outFile) 
{
  // scattered muon histo
  TH3D* in = (TH3D*) inFile->Get(inHistoName);

      int nBinsX = in->GetNbinsX();
      int nBinsY = in->GetNbinsY();
      int nBinsZ = in->GetNbinsZ();

      float xAxisMin = in->GetXaxis()->GetXmin();
      float xAxisMax = in->GetXaxis()->GetXmax();

      float yAxisMin = in->GetYaxis()->GetXmin();
      float yAxisMax = in->GetYaxis()->GetXmax();

      float zAxisMin = in->GetZaxis()->GetXmin();
      float zAxisMax = in->GetZaxis()->GetXmax();

      float binWidth = in->GetXaxis()->GetBinWidth(1);
  
      float threshold = modifier;
  
  //output histo
      TH3D* out = new TH3D(outHistName, outHistName, nBinsX,xAxisMin,xAxisMax, nBinsY,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax);

  // absorbed muon histo
      double xCut, xGrad, zCut, zGrad;

      TTree* abs = (TTree*) absorbedFile->Get(absTreeName);
  
      abs->SetBranchAddress("xCut",&xCut);
      abs->SetBranchAddress("zCut",&zCut);

      abs->SetBranchAddress("xGrad",&xGrad);
      abs->SetBranchAddress("zGrad",&zGrad);

      int nEntries = abs->GetEntries();

      for(int entry = 0; entry < nEntries; entry++)
  {
  abs->GetEntry(entry);

          if(xGrad == 0 || zGrad == 0) continue; // if there is null gradient ignore and continue      

          int counter = (2*extent + 1); // tracks and 'catterpillars' the face
          double face[counter]; // will take the value for each new face of the surronding cube
      
          for (int n = 0; n < counter; n++) //initialise the array
  {
  face[n] = 0;
	}

          for (int j = 0; j < nBinsY; j++) // tracks the absorbed muon through the box
  {
  if(counter >= (2*extent + 1))
  {
  counter = 0; // causes catterpillar effect
	    }
	  
	      face[counter] = 0; // reset the array 
	 
	  // calculate position by bin
	      double y = (j*binWidth) + yAxisMin; // the 'real' value of y
	      double x = (y - xCut)/ xGrad;  // calculate x and z positions using the y value (x = (y-cut)/grad)
	      double z = (y - zCut)/ zGrad;  

	      int i = (int) ((x - xAxisMin) / binWidth ); // calculate the bin number for x and z
	      int k = (int) ((z - zAxisMin) / binWidth );
	  
	      for (int iO = -extent; iO <= extent; iO++) // cycle through the surronding values
  {
  for (int kO = -extent; kO <= extent; kO++)
  {
  face[counter] += in->GetBinContent((i+iO), (j + extent),(k+kO)); // cycle through the face to be added
		}
	    }

	      double testVal = 0; // the value that will take the density of points
	  
	      for (int faceN = 0; faceN < (2*extent + 1); faceN ++)
  {
  testVal += face[faceN]; 
	    }

	      if ((testVal > threshold) && !(out->GetBinContent (i,j,k)) )
  {
  out->SetBinContent (i,j,k,1);
	    }
	      counter++;
	}

          if (!(entry%100)) std::cout<<entry<<" of "<<nEntries<<std::endl; 
    }
      outFile->Add(out);
}

// this will create a series of 2D histograms showing slices (in y) of the absorbed muons and where they pass through that slice
void absorbedSlices (char* inHistoName, char* absTreeName, char* outHistName, int numberOfSlices,float modifier,TFile* inFile, TFile* absorbedFile, TFile* outFile) 
{
  // input histogram

  TH3D* in = (TH3D*) inFile->Get(inHistoName);

      int nBinsX =1000;// in->GetNbinsX();
      int nBinsZ =1000;// in->GetNbinsZ();

      float xAxisMin = in->GetXaxis()->GetXmin();
      float xAxisMax = in->GetXaxis()->GetXmax();

      float yAxisMin = in->GetYaxis()->GetXmin();
      float yAxisMax = in->GetYaxis()->GetXmax();

      float zAxisMin = in->GetZaxis()->GetXmin();
      float zAxisMax = in->GetZaxis()->GetXmax();

      float sliceWidth = (in->GetYaxis()->GetXmax())/ numberOfSlices;

  // absorbed muon histo
      double xCut, xGrad, zCut, zGrad;

      TTree* abs = (TTree*) absorbedFile->Get(absTreeName);
  
      abs->SetBranchAddress("xCut",&xCut);
      abs->SetBranchAddress("zCut",&zCut);

      abs->SetBranchAddress("xGrad",&xGrad);
      abs->SetBranchAddress("zGrad",&zGrad);

      int nEntries = abs->GetEntries();

      float threshold = modifier + std::sqrt(nEntries / (nBinsZ * nBinsX));
      float density = nEntries / (nBinsZ * nBinsX);

      TH3D* out2 = new TH3D (outHistName, outHistName, nBinsX,xAxisMin,xAxisMax, numberOfSlices,yAxisMin,yAxisMax, nBinsZ,zAxisMin,zAxisMax); 

      for (int j = 0; j < numberOfSlices; j++)
  {
  double y = (j*sliceWidth) + yAxisMin; // the 'real' value of y
          char name [25];
          char yVal [3];
	  
          sprintf (yVal, "%d",j);

      //    std::cout<<name<<" "<<yVal<<std::endl;
      
          strcpy (name,outHistName);
          strcat (name, yVal);
      
      //    std::cout<<name<<std::endl;

          TH2D* out = new TH2D (name,name, nBinsX,xAxisMin,xAxisMax, nBinsZ,zAxisMin,zAxisMax);

          for (int entry = 0; entry < nEntries; entry++)
  {
  abs->GetEntry(entry);
	  
	      if(xGrad == 0 || zGrad == 0) continue; // if there is null gradient ignore and continue      

	      double x = (y - xCut)/ xGrad;  // calculate x and z positions using the y value (x = (y-cut)/grad)
	      double z = (y - zCut)/ zGrad;  

	      out->Fill (x,z);

	      if (!(entry%5000)) std::cout<<(j*nEntries + entry)<<" of "<<(nEntries*numberOfSlices)<<" j is "<< j<<std::endl;
	}


      // this section finds the regions of interest by removing the noise,
      // if the point has an above mean (by <modifier> std dev) then it 
      // is marked otherwise it is removed.

          for (int i =0; i < nBinsX; i++)
  {
  for (int k = 0; k < nBinsZ; k++)
  {
  float test = out->GetBinContent (i,k);
	          if (std::fabs (test - density) > threshold)
  {
  out->SetBinContent (i,k,1);
		      out2->SetBinContent (i,j,k,1);
		}
	          else
  {
  out->SetBinContent (i,k,0);
		}
	    }
	}
      
          outFile->Add(out2);
          outFile->Add(out);
    }
}

void absorbedSubtraction (char* bkgndTreeName, char* absTreeName, char* outHistName,  // names of the background and input tress as well as the name to give the output histo
			  int numberOfSlices,float sliceWidth ,                       // number of slices and their width in y
			  int extent, float modifier,                                 // extent to compare background to source over and modifier for threshold
			  TFile* backgroundFile, TFile* absorbedFile, TFile* outFile, // files to use
			  float startY, int nBinsX, int nBinsZ,                      // optional starting position for the y slice and number of bins to use in x and z
			  float xMax, float xMin,float zMax, float zMin)              // values for the range of x and z)     
{
  // ******************** absorbed muon tree
  double xCut, xGrad, zCut, zGrad;

  TTree* abs = (TTree*) absorbedFile->Get(absTreeName);
  
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


// basic load function to import threshold and range values as well as the name of the root file to use and the histogram to use
int loader (int* extent, int* sliceNumber, float* dtModifier, float* gradModifier, char* inputFile, char* histoName)
  //returns 0 for error 1 for anything else
{
  char v_inputFile[256];
      char v_histoName[256];

  // default values
      int  v_extent = 1;
      int v_sliceNumber = 1;
      float v_dtModifier = 0;
      float v_gradModifier =0;

      ifstream in;
      in.open (SOURCE_F);
      in>>v_extent>>v_sliceNumber>>v_dtModifier>>v_gradModifier>>v_inputFile>>v_histoName;
      in.close();
  
      if(!v_inputFile[0])
  {
  std::cout<<"no file loaded. Quiting"<<std::endl;
      
          return 0;      
    }
      else if(!v_histoName[0])
  {
  std::cout<<"no histogram specified. Quitting"<<std::endl;

          return 0;
    }
      else
  {      
  *extent = v_extent;
          *sliceNumber = v_sliceNumber;
          *dtModifier = v_dtModifier;
          *gradModifier = v_gradModifier;

          std::strcpy(inputFile,v_inputFile);
          std::strcpy(histoName, v_histoName);

          return 1;
    }
}

// small function to allow all three directions to be tested
double permute(TH3D* in, int index, int x, int y, int z, char option)
  // This is designed so i can cycle through the various permutations
  // of x y and z
  // the x y and z are the 'native values' while the index dictates
  // the order in which to use them in the 'in' histogram
  //
  // index 1=(x,y,z) 2=(y,z,x) 3=(z,x,y)
  //
  // option defines either fill (F or f) or getBinContent (G or g) 
{
  if(option =='G'||option =='g')
  {
  switch (index)
  {
	    case 0: 
  return (in->GetBinContent (x,y,z));
	      break;
	    case 1:
  return (in->GetBinContent (y,z,x));
	      break;
	    case 2: 
  return (in->GetBinContent (z,x,y));
	      break;
	    default:
  return 0;
	      break;
	}
    }
      else if(option =='F'||option =='f')
  {
  switch (index)
  {
	    case 0: 
  in->SetBinContent (x,y,z,1);
	      break;
	    case 1:
  in->SetBinContent (y,z,x,1);
	      break;
	    case 2: 
  in->SetBinContent (z,x,y,1);
	      break;
	    default:
  return 0;
	      break;
	}
    }

      return 0;
}








