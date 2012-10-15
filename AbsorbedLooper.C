#define AbsorbedLooper_cxx
#include "AbsorbedLooper.h"
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include <list>
#define NUM_ITS 1

void AbsorbedLooper::MakeSliceHists(int numBins)
{
//   In a ROOT session, you can do:
//      Root > .L AbsorbedLooper.C
//      Root > AbsorbedLooper t
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

   const int dz=10; //mm
   const int numSlices=2*(500/dz)+1;
   TH2F *histAbsorbed[numSlices]={0};
   Int_t zSlice[numSlices]={0};
   char histName[180];
   char histTitle[180];
   for(int i=0;i<numSlices;i++) {
     int z=-500+i*dz;
     zSlice[i]=z;
     sprintf(histName,"histSlice_%d",z);			
     sprintf(histTitle,"Slice at Z=%d mm",z);
     histAbsorbed[i]=new TH2F(histName,histTitle,numBins,-500,500,
			      numBins,-500,500);
     histAbsorbed[i]->SetXTitle("x-position (mm)");
     histAbsorbed[i]->SetYTitle("y-position (mm)");
   }
   

   
   TH1D *histxcut =  new TH1D("histxcut","histxcut",2000,-1000,1000);
   
   Long64_t nentries = fChain->GetEntries();

   std::cout << "Entries: " << nentries << "\n";
   Long64_t nbytes = 0, nb = 0;
   Int_t countEntries=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) < 0) continue;
      countEntries++;
      histxcut->Fill(xCut);

      for(int slice=0;slice<numSlices;slice++) {
	Double_t xPos=xGrad*zSlice[slice]+xCut;
	Double_t yPos=yGrad*zSlice[slice]+yCut;
	histAbsorbed[slice]->Fill(xPos,yPos);
      }
   }
   std::cout << "Found " << countEntries << " stoppers.\n";

   
  TCanvas *canNeg = new TCanvas("canNeg","canNeg");
  canNeg->Divide(4,4);
  for(int slice=0;slice<50;slice++) {    
    canNeg->cd(slice+1);
    histAbsorbed[slice]->Draw("colz");
  }

  TCanvas *canPos = new TCanvas("canPos","canPos");
  canPos->Divide(4,4);
  for(int slice=50;slice<100;slice++) {    
    canPos->cd(slice-49);
    histAbsorbed[slice]->Draw("colz");
  }
  
  TCanvas *canTarget = new TCanvas("canTarget","canTarget");
  histAbsorbed[50]->Draw("colz");

  TCanvas *xcutcanvas = new TCanvas();
  histxcut->Draw();
  
}



void AbsorbedLooper::MakeSliceHistsIteratively(int binWidth)
{
//   In a ROOT session, you can do:
//      Root > .L AbsorbedLooper.C
//      Root > AbsorbedLooper t
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

   int dz=binWidth;
   const int numSlices=2*(500/dz)+1;
   std::cout << "Using " << numSlices << " slices\n";
   char outName[FILENAME_MAX];
   sprintf(outName,"/unix/creamtea/sfayer/absorbed/absorbedOut_%d.root",binWidth);
   TFile *fp = new TFile (outName,"RECREATE");
   TH2F *histRawAbsorbed[numSlices];
   TH2F *histWeightedAbsorbed[numSlices];
   TH2F *histWeightedAbsorbedIter[NUM_ITS][numSlices];
   TH2F *histGaussianWeightedAbsorbed[numSlices];
   Int_t zSlice[numSlices];
   Double_t dzSlice[numSlices];
   Int_t numBins=numSlices;
   char histName[180];
   char histTitle[180];
   Double_t maxXY=500+binWidth/2;
   for(int i=0;i<numSlices;i++) {
     int z=-500+i*dz;
     zSlice[i]=z;
     dzSlice[i]=z;
     sprintf(histName,"histRawSlice_%d",(z+500));			
     sprintf(histTitle,"Raw Absorptions, Slice at Z=%d mm",(z+500));     
     histRawAbsorbed[i]=new TH2F(histName,histTitle,numBins,-maxXY,maxXY,
				 numBins,-maxXY,maxXY);
     histRawAbsorbed[i]->SetXTitle("x-position (mm)");
     histRawAbsorbed[i]->SetYTitle("y-position (mm)");
     sprintf(histName,"histWeightedSlice_%d",(z+500));		
     sprintf(histTitle,"Weighted Absorptions, Slice at Z=%d mm",(z+500));     	
     histWeightedAbsorbed[i]=new TH2F(histName,histTitle,numBins,
					  -maxXY,maxXY,
					  numBins,-maxXY,maxXY);
     histWeightedAbsorbed[i]->SetXTitle("x-position (mm)");
     histWeightedAbsorbed[i]->SetYTitle("y-position (mm)");
     sprintf(histName,"histGaussianWeightedSlice_%d",(z+500));		
     sprintf(histTitle,"Gaussian Weighted Absorptions, Slice at Z=%d mm",(z+500));     	
     histGaussianWeightedAbsorbed[i]=new TH2F(histName,histTitle,numBins,
					  -maxXY,maxXY,
					  numBins,-maxXY,maxXY);
     histGaussianWeightedAbsorbed[i]->SetXTitle("x-position (mm)");
     histGaussianWeightedAbsorbed[i]->SetYTitle("y-position (mm)");

//      for(int it=0;it<NUM_ITS;it++) {
//        sprintf(histName,"histWeightedSliceIt%d_%d",it,z);		
//        sprintf(histTitle,"Weighted Absorptions (Iter %d), Slice at Z=%d mm",it,z); 
//        histWeightedAbsorbedIter[it][i]=new TH2F(histName,histTitle,numBins,
// 						-maxXY,maxXY,
// 						numBins,-maxXY,maxXY);
//        histWeightedAbsorbedIter[it][i]->SetXTitle("x-position (mm)");
//        histWeightedAbsorbedIter[it][i]->SetYTitle("y-position (mm)");
//      }
   }

   std::list<Long64_t> validStopperEntryList;
   
   
   Long64_t nentries = fChain->GetEntries();

   std::cout << "Entries: " << nentries << "\n";
   Long64_t nbytes = 0, nb = 0;
   Int_t countEntries=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (Cut(ientry) < 0) continue;
      countEntries++;
      validStopperEntryList.push_back(jentry);
      for(int slice=0;slice<numSlices;slice++) {
	Double_t xPos=xGrad*zSlice[slice]+xCut;
	Double_t yPos=yGrad*zSlice[slice]+yCut;
	histRawAbsorbed[slice]->Fill(xPos,yPos);
      }
   }
   std::cout << "Found " << countEntries << " stoppers.\n";
   
   std::list<Long64_t>::iterator validIt;     
   Int_t binx,biny,bin;
   Double_t arbThresh=1./numSlices;
   arbThresh=0;
   
   char graphName[180];
   //Now loop over the valid entries
   Int_t muonCount=0;
   TF1 *fitty = new TF1("fitty","gaus",-500,500);

   for(validIt=validStopperEntryList.begin(); 
       validIt != validStopperEntryList.end() ; 
       validIt++) {
     Long64_t entryNum=(*validIt);     
     Long64_t ientry = LoadTree(entryNum);
     if (ientry < 0) break;
     nb = fChain->GetEntry(entryNum);   nbytes += nb;
     
     //Now for the clever bit, well it is not really all that clever but it is
     // moderately clever. We find the bin in each slice that the muon passes 
     // through and calculate a weight (for each slice on each track) that is
     // proportional to the number of coincident tracks in each voxel.
     Double_t sumNxyz=0;
     Double_t weights[numSlices];
     for(int slice=0;slice<numSlices;slice++) {
       Double_t xPos=xGrad*zSlice[slice]+xCut;
       Double_t yPos=yGrad*zSlice[slice]+yCut;
       binx=histRawAbsorbed[slice]->GetXaxis()->FindBin(xPos);
       biny=histRawAbsorbed[slice]->GetXaxis()->FindBin(yPos);
       weights[slice]=histRawAbsorbed[slice]->GetBinContent(binx,biny)-1;
       sumNxyz+=weights[slice];
     }
     
     
     for(int slice=0;slice<numSlices;slice++) {
       weights[slice]/=sumNxyz;
       //       std::cout << slice << "\t" << weights[slice] << "\t" << sumNxyz << "\n";
       Double_t xPos=xGrad*zSlice[slice]+xCut;
       Double_t yPos=yGrad*zSlice[slice]+yCut;
       if(weights[slice]>arbThresh)
	 histWeightedAbsorbed[slice]->Fill(xPos,yPos,weights[slice]);
     }          
     sprintf(graphName,"grMuon%d",muonCount);
     fitty->SetParameters(0,0,1);
     TGraph *gr = new TGraph(numSlices,dzSlice,weights);
     gr->SetName(graphName);
     gr->Fit("fitty","QR");     
     gr->Write();
     muonCount++;


     if(fitty->GetParameter(0)>0.05 && TMath::Abs(fitty->GetParameter(1))<maxXY) {
       //       std::cout << muonCount-1 << endl;
       for(int slice=0;slice<numSlices;slice++) {
	 Double_t xPos=xGrad*zSlice[slice]+xCut;
	 Double_t yPos=yGrad*zSlice[slice]+yCut;
	 histGaussianWeightedAbsorbed[slice]->Fill(xPos,yPos,fitty->Eval(dzSlice[slice]));
       }
     }
						 
     
   }
   //   Same again but using the weighted histograms as input


   
 //   for(int it=1;it<NUM_ITS;it++) {
//      muonCount=0;
//      std::cout << "Starting iteration number " << it << "\n";
//      //Now loop over the valid entries
//      for(validIt=validStopperEntryList.begin(); 
// 	 validIt != validStopperEntryList.end() ; 
// 	 validIt++) {
//        Long64_t entryNum=(*validIt);     
//        Long64_t ientry = LoadTree(entryNum);
//        if (ientry < 0) break;
//        nb = fChain->GetEntry(entryNum);   nbytes += nb;
       
//        //Now for the clever bit, well it is not really all that clever but it is
//        // moderately clever. We find the bin in each slice that the muon passes 
//        // through and calculate a weight (for each slice on each track) that is
//        // proportional to the number of coincident tracks in each voxel.
//        Double_t sumNxyz=0;
//        Double_t weights[numSlices];
//        for(int slice=0;slice<numSlices;slice++) {
// 	 Double_t xPos=xGrad*zSlice[slice]+xCut;
// 	 Double_t yPos=yGrad*zSlice[slice]+yCut;
// 	 if(it==1) {
// 	   binx=histWeightedAbsorbed[slice]->GetXaxis()->FindBin(xPos);
// 	   biny=histWeightedAbsorbed[slice]->GetXaxis()->FindBin(yPos);
// 	   weights[slice]=histWeightedAbsorbed[slice]->GetBinContent(binx,biny);
// 	 }
// 	 else {
// 	   binx=histWeightedAbsorbedIter[it-1][slice]->GetXaxis()->FindBin(xPos);
// 	   biny=histWeightedAbsorbedIter[it-1][slice]->GetXaxis()->FindBin(yPos);
// 	   weights[slice]=histWeightedAbsorbedIter[it-1][slice]->GetBinContent(binx,biny);
// 	 }


// 	 sumNxyz+=weights[slice];
//        }
       
//        for(int slice=0;slice<numSlices;slice++) {
// 	 weights[slice]/=sumNxyz;
// 	 Double_t xPos=xGrad*zSlice[slice]+xCut;
// 	 Double_t yPos=yGrad*zSlice[slice]+yCut;
// 	 histWeightedAbsorbedIter[it][slice]->Fill(xPos,yPos,weights[slice]);
//        }          
       
//        sprintf(graphName,"grMuon%d_%d",muonCount,it);
//        fitty->SetParameters(0,0,1);
//        TGraph *gr = new TGraph(numSlices,dzSlice,weights);
//        gr->SetName(graphName);
//        gr->Fit("fitty","QR");     
//        gr->Write();
//        muonCount++;
//      }
//    }

   Double_t maxRaw=0;
   Double_t maxWeighted=0;
   Double_t maxGaussian=0;
   for(int slice=0;slice<numSlices;slice++) {    
     if(histRawAbsorbed[slice]->GetMaximum()>maxRaw)
       maxRaw=histRawAbsorbed[slice]->GetMaximum();
     if(histWeightedAbsorbed[slice]->GetMaximum()>maxWeighted)
       maxWeighted=histWeightedAbsorbed[slice]->GetMaximum();
     if(histGaussianWeightedAbsorbed[slice]->GetMaximum())
       maxGaussian=histGaussianWeightedAbsorbed[slice]->GetMaximum();
   }

  TCanvas *canRaw = new TCanvas("canRaw","canRaw");
  Int_t numCols=TMath::Sqrt(numSlices);
  Int_t numRows=numSlices/numCols;  
  if(numSlices%numCols) numRows++;

  canRaw->Divide(numCols,numRows,0,0);
  for(int slice=0;slice<numSlices;slice++) {    
    canRaw->cd(slice+1);
    gPad->SetLogz();
    histRawAbsorbed[slice]->SetStats(0);
    histRawAbsorbed[slice]->SetMinimum(1);
    histRawAbsorbed[slice]->SetMaximum(maxRaw);
    histRawAbsorbed[slice]->Draw("colz");
  }
  char canName[180];
  sprintf(canName,"canWeighted");
  TCanvas *canWeighted = new TCanvas(canName,canName);
  canWeighted->Divide(numCols,numRows,0,0);
  for(int slice=0;slice<numSlices;slice++) {    
    canWeighted->cd(slice+1);
    gPad->SetLogz();
    histWeightedAbsorbed[slice]->SetStats(0);
    histWeightedAbsorbed[slice]->SetMinimum(0.0001);
    histWeightedAbsorbed[slice]->SetMaximum(maxWeighted);
    histWeightedAbsorbed[slice]->Draw("colz");
  }


  sprintf(canName,"canGaussianWeighted");
  TCanvas *canGaussianWeighted = new TCanvas(canName,canName);
  canGaussianWeighted->Divide(numCols,numRows,0,0);
  for(int slice=0;slice<numSlices;slice++) {    
    canGaussianWeighted->cd(slice+1);
    gPad->SetLogz();
    histGaussianWeightedAbsorbed[slice]->SetStats(0);
    histGaussianWeightedAbsorbed[slice]->SetMinimum(0.0001);
    histGaussianWeightedAbsorbed[slice]->SetMaximum(maxGaussian);
    histGaussianWeightedAbsorbed[slice]->Draw("colz");
  }

  
  for(int slice=(numSlices/2)-6;slice<=(numSlices/2)+6;slice++) {
    sprintf(canName,"canSlice%d",slice);
    TCanvas *canSlice = new TCanvas(canName,canName);
    canSlice->Divide(1,3);
    canSlice->cd(1);
    gPad->SetLogz(1);
    histRawAbsorbed[slice]->SetMinimum(1);
    histRawAbsorbed[slice]->Draw("colz");

    canSlice->cd(2);
    gPad->SetLogz(1);
    histWeightedAbsorbed[slice]->SetMinimum(0.05);
    histWeightedAbsorbed[slice]->Draw("colz");
    canSlice->cd(3);
    gPad->SetLogz(1);
    histGaussianWeightedAbsorbed[slice]->SetMinimum(0.05);
    histGaussianWeightedAbsorbed[slice]->Draw("colz");
  }
  
  fp->Write();
}



void AbsorbedLooper::MakeSliceHistsIterativelyReco(int binWidth)
{
//   In a ROOT session, you can do:
//      Root > .L AbsorbedLooper.C
//      Root > AbsorbedLooper t
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

   int dz=binWidth;
   const int numSlices=2*(500/dz)+1;
   std::cout << "Using " << numSlices << " slices\n";
   TFile *fp = new TFile ("/unix/creamtea/sfayer/absorbed/absorbedRecoOut.root","RECREATE");
   TH2F *histRawAbsorbed[numSlices];
   TH2F *histWeightedAbsorbed[numSlices];
   TH2F *histGaussianWeightedAbsorbed[numSlices];
   Int_t zSlice[numSlices];
   Double_t dzSlice[numSlices];
   Int_t numBins=numSlices;
   char histName[180];
   char histTitle[180];
   Double_t maxXY=500+binWidth/2;
   for(int i=0;i<numSlices;i++) {
     int z=-500+i*dz;
     zSlice[i]=z;
     dzSlice[i]=z;
     sprintf(histName,"histRecoRawSlice_%d",z);			
     sprintf(histTitle,"Raw Absorptions, Slice at Z=%d mm",z);     
     histRawAbsorbed[i]=new TH2F(histName,histTitle,numBins,-maxXY,maxXY,
				 numBins,-maxXY,maxXY);
     histRawAbsorbed[i]->SetXTitle("x-position (mm)");
     histRawAbsorbed[i]->SetYTitle("y-position (mm)");
     sprintf(histName,"histRecoWeightedSlice_%d",z);		
     sprintf(histTitle,"Weighted Absorptions, Slice at Z=%d mm",z);     	
     histWeightedAbsorbed[i]=new TH2F(histName,histTitle,numBins,
					  -maxXY,maxXY,
					  numBins,-maxXY,maxXY);
     histWeightedAbsorbed[i]->SetXTitle("x-position (mm)");
     histWeightedAbsorbed[i]->SetYTitle("y-position (mm)");
     sprintf(histName,"histRecoGaussianWeightedSlice_%d",z);		
     sprintf(histTitle,"Gaussian Weighted Absorptions, Slice at Z=%d mm",z);     	
     histGaussianWeightedAbsorbed[i]=new TH2F(histName,histTitle,numBins,
					  -maxXY,maxXY,
					  numBins,-maxXY,maxXY);
     histGaussianWeightedAbsorbed[i]->SetXTitle("x-position (mm)");
     histGaussianWeightedAbsorbed[i]->SetYTitle("y-position (mm)");
   }

   std::list<Long64_t> validStopperEntryList;
   
   
   Long64_t nentries = fChain->GetEntries();

   std::cout << "Entries: " << nentries << "\n";
   Long64_t nbytes = 0, nb = 0;
   Int_t countEntries=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (CutReco(ientry) < 0) continue;
      countEntries++;
      validStopperEntryList.push_back(jentry);
      for(int slice=0;slice<numSlices;slice++) {
	Double_t xPos=xGradReco*zSlice[slice]+xCutReco;
	Double_t yPos=yGradReco*zSlice[slice]+yCutReco;
	histRawAbsorbed[slice]->Fill(xPos,yPos);
      }
   }
   std::cout << "Found " << countEntries << " stoppers.\n";
   
   std::list<Long64_t>::iterator validIt;     
   Int_t binx,biny,bin;
   Double_t arbThresh=1./numSlices;
   arbThresh=0;

   char graphName[180];
   //Now loop over the valid entries
   Int_t muonCount=0;
   TF1 *fitty = new TF1("fitty","gaus",-500,500);

   for(validIt=validStopperEntryList.begin(); 
       validIt != validStopperEntryList.end() ; 
       validIt++) {
     Long64_t entryNum=(*validIt);     
     Long64_t ientry = LoadTree(entryNum);
     if (ientry < 0) break;
     nb = fChain->GetEntry(entryNum);   nbytes += nb;
     
     //Now for the clever bit, well it is not really all that clever but it is
     // moderately clever. We find the bin in each slice that the muon passes 
     // through and calculate a weight (for each slice on each track) that is
     // proportional to the number of coincident tracks in each voxel.
     Double_t sumNxyz=0;
     Double_t weights[numSlices];
     for(int slice=0;slice<numSlices;slice++) {
       Double_t xPos=xGrad*zSlice[slice]+xCut;
       Double_t yPos=yGrad*zSlice[slice]+yCut;
       binx=histRawAbsorbed[slice]->GetXaxis()->FindBin(xPos);
       biny=histRawAbsorbed[slice]->GetXaxis()->FindBin(yPos);
       weights[slice]=histRawAbsorbed[slice]->GetBinContent(binx,biny)-1;
       sumNxyz+=weights[slice];
     }
     
     
     for(int slice=0;slice<numSlices;slice++) {
       weights[slice]/=sumNxyz;
       //       std::cout << slice << "\t" << weights[slice] << "\t" << sumNxyz << "\n";
       Double_t xPos=xGrad*zSlice[slice]+xCut;
       Double_t yPos=yGrad*zSlice[slice]+yCut;
       if(weights[slice]>arbThresh)
	 histWeightedAbsorbed[slice]->Fill(xPos,yPos,weights[slice]);
     }          
     sprintf(graphName,"grMuon%d",muonCount);
     fitty->SetParameters(0,0,1);
     TGraph *gr = new TGraph(numSlices,dzSlice,weights);
     gr->SetName(graphName);
     gr->Fit("fitty","QR");     
     gr->Write();
     muonCount++;


     if(fitty->GetParameter(0)>0.05 && TMath::Abs(fitty->GetParameter(1))<maxXY) {
       //       std::cout << muonCount-1 << endl;
       for(int slice=0;slice<numSlices;slice++) {
	 Double_t xPos=xGrad*zSlice[slice]+xCut;
	 Double_t yPos=yGrad*zSlice[slice]+yCut;
	 histGaussianWeightedAbsorbed[slice]->Fill(xPos,yPos,fitty->Eval(dzSlice[slice]));
       }
     }
						 
     
   }
   //Same again but using the weighted histograms as input
//    for(int it=1;it<5;it++) {
//      //Now loop over the valid entries
//      for(validIt=validStopperEntryList.begin(); 
// 	 validIt != validStopperEntryList.end() ; 
// 	 validIt++) {
//        Long64_t entryNum=(*validIt);     
//        Long64_t ientry = LoadTree(entryNum);
//        if (ientry < 0) break;
//        nb = fChain->GetEntry(entryNum);   nbytes += nb;
       
//        //Now for the clever bit, well it is not really all that clever but it is
//        // moderately clever. We find the bin in each slice that the muon passes 
//        // through and calculate a weight (for each slice on each track) that is
//        // proportional to the number of coincident tracks in each voxel.
//      Double_t sumNxyz=0;
//      Double_t weights[numSlices];
//      for(int slice=0;slice<numSlices;slice++) {
//        Double_t xPos=xGrad*zSlice[slice]+xCut;
//        Double_t yPos=yGrad*zSlice[slice]+yCut;
//        binx=histWeightedAbsorbed[it-1][slice]->GetXaxis()->FindBin(xPos);
//        biny=histWeightedAbsorbed[it-1][slice]->GetXaxis()->FindBin(yPos);
//        weights[slice]=histWeightedAbsorbed[it-1][slice]->GetBinContent(binx,biny);
//        sumNxyz+=weights[slice];
//      }

//      for(int slice=0;slice<numSlices;slice++) {
//        weights[slice]/=sumNxyz;
//        Double_t xPos=xGrad*zSlice[slice]+xCut;
//        Double_t yPos=yGrad*zSlice[slice]+yCut;
//        histWeightedAbsorbed[slice]->Fill(xPos,yPos,weights[slice]);
//      }          
//      }
//    }

   Double_t maxRaw=0;
   Double_t maxWeighted=0;
   Double_t maxGaussian=0;
   for(int slice=0;slice<numSlices;slice++) {    
     if(histRawAbsorbed[slice]->GetMaximum()>maxRaw)
       maxRaw=histRawAbsorbed[slice]->GetMaximum();
     if(histWeightedAbsorbed[slice]->GetMaximum()>maxWeighted)
       maxWeighted=histWeightedAbsorbed[slice]->GetMaximum();
     if(histGaussianWeightedAbsorbed[slice]->GetMaximum())
       maxGaussian=histGaussianWeightedAbsorbed[slice]->GetMaximum();
   }

  TCanvas *canRaw = new TCanvas("canRaw","canRaw");
  Int_t numCols=TMath::Sqrt(numSlices);
  Int_t numRows=numSlices/numCols;  
  if(numSlices%numCols) numRows++;

  canRaw->Divide(numCols,numRows,0,0);
  for(int slice=0;slice<numSlices;slice++) {    
    canRaw->cd(slice+1);
    gPad->SetLogz();
    histRawAbsorbed[slice]->SetStats(0);
    histRawAbsorbed[slice]->SetMinimum(1);
    histRawAbsorbed[slice]->SetMaximum(maxRaw);
    histRawAbsorbed[slice]->Draw("colz");
  }
  char canName[180];
  sprintf(canName,"canWeighted");
  TCanvas *canWeighted = new TCanvas(canName,canName);
  canWeighted->Divide(numCols,numRows,0,0);
  for(int slice=0;slice<numSlices;slice++) {    
    canWeighted->cd(slice+1);
    gPad->SetLogz();
    histWeightedAbsorbed[slice]->SetStats(0);
    histWeightedAbsorbed[slice]->SetMinimum(0.0001);
    histWeightedAbsorbed[slice]->SetMaximum(maxWeighted);
    histWeightedAbsorbed[slice]->Draw("colz");
  }


  sprintf(canName,"canGaussianWeighted");
  TCanvas *canGaussianWeighted = new TCanvas(canName,canName);
  canGaussianWeighted->Divide(numCols,numRows,0,0);
  for(int slice=0;slice<numSlices;slice++) {    
    canGaussianWeighted->cd(slice+1);
    gPad->SetLogz();
    histGaussianWeightedAbsorbed[slice]->SetStats(0);
    histGaussianWeightedAbsorbed[slice]->SetMinimum(0.0001);
    histGaussianWeightedAbsorbed[slice]->SetMaximum(maxGaussian);
    histGaussianWeightedAbsorbed[slice]->Draw("colz");
  }

  
  for(int slice=(numSlices/2)-6;slice<=(numSlices/2)+6;slice++) {
    sprintf(canName,"canSlice%d",slice);
    TCanvas *canSlice = new TCanvas(canName,canName);
    canSlice->Divide(1,3);
    canSlice->cd(1);
    gPad->SetLogz(1);
    histRawAbsorbed[slice]->SetMinimum(1);
    histRawAbsorbed[slice]->Draw("colz");

    canSlice->cd(2);
    gPad->SetLogz(1);
    histWeightedAbsorbed[slice]->SetMinimum(0.05);
    histWeightedAbsorbed[slice]->Draw("colz");
    canSlice->cd(3);
    gPad->SetLogz(1);
    histGaussianWeightedAbsorbed[slice]->SetMinimum(0.05);
    histGaussianWeightedAbsorbed[slice]->Draw("colz");
  }
  
  fp->Write();
}
