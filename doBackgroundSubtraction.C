


void doBackgroundSubtraction()
{
  char *backgroundBase="/unix/anita1/creamtea/strips_650/nocontainer/pca/pca_nocontainer";
  //  char *targetBase="/unix/anita1/creamtea/strips_650/container_30cmtungsten/pca/pca_container_tungsten30cm";
  char *targetBase="/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm";
  //  char *targetBase="/unix/anita1/creamtea/strips_650/container/pca/pca_container";

  Int_t numBackground=10000;
  Int_t numTarget=10000;
  
  char fileName[180];

  TChain *backTree = new TChain("pcaTree");
  for(Int_t startFile=1;startFile<numBackground;startFile+=1000) {
    Int_t endFile=startFile+999;
    sprintf(fileName,"%s_%d_%d.root",backgroundBase,startFile,endFile);
    backTree->Add(fileName);
    //    cout << fileName << endl;
  }


  TChain *targetTree = new TChain("pcaTree");
  for(Int_t startFile=1;startFile<numTarget;startFile+=1000) {
    Int_t endFile=startFile+999;
    sprintf(fileName,"%s_%d_%d.root",targetBase,startFile,endFile);
    targetTree->Add(fileName);
    //    cout << fileName << endl;
  }

  gSystem->CompileMacro("PcaTreeLooper.C","k");
  PcaTreeLooper backLooper(backTree);
  PcaTreeLooper targetLooper(targetTree);


  std::cout << backTree->GetEntries() << "\t" <<  targetTree->GetEntries() << std::endl;

  const Int_t numBins=26;
  TFile *fpOut = new TFile("outHists.root","RECREATE");
  TH3F *histBack = new TH3F("histBack","Background",numBins,-6500,6500,numBins,-6500,6500,numBins,-6500,6500);
  TH3F *histTarget = new TH3F("histTarget","Target",numBins,-6500,6500,numBins,-6500,6500,numBins,-6500,6500);
  TH3F *histDiff = new TH3F("histDiff","Difference",numBins,-6500,6500,numBins,-6500,6500,numBins,-6500,6500);
  TH3F *histDiffPos = new TH3F("histDiffPos","Positive Difference",numBins,-6500,6500,numBins,-6500,6500,numBins,-6500,6500);
  TH3F *histDiffNeg = new TH3F("histDiffNeg","Negative Difference",numBins,-6500,6500,numBins,-6500,6500,numBins,-6500,6500);
  

  TH2F *histDiffPosXY = new TH2F("histDiffPosXY","Positive Difference XY",numBins,-6500,6500,numBins,-6500,6500);
  TH2F *histDiffPosXZ = new TH2F("histDiffPosXZ","Positive Difference XZ",numBins,-6500,6500,numBins,-6500,6500);
  TH2F *histDiffPosYZ = new TH2F("histDiffPosYZ","Positive Difference YZ",numBins,-6500,6500,numBins,-6500,6500);

  TH2F *histDiffPosXYSlice[numBins];
  char histName[180];
  char histTitle[180];
  for(int i=0;i<numBins;i++) {
    sprintf(histName,"histDiffPosXYSlice%d",i+1);
    sprintf(histTitle,"Positive Difference XY (Slice %d)",i+1);
    histDiffPosXYSlice[i] = new TH2F(histName,histTitle,numBins,-6500,6500,numBins,-6500,6500);    

  }

  //  char plotCond[180];
  //  sprintf(plotCond,"thetaTrue>%f",0);
  //  backTree->Project("histBack","xPosTrue:yPosTrue:zPosTrue",plotCond);
  //  targetTree->Project("histTarget","xPosTrue:yPosTrue:zPosTrue",plotCond);


  backLooper.FillPosHist(histBack,0.2);
  targetLooper.FillPosHist(histTarget,0.2);

  
  Int_t numSigma=4;

  for(int binx=1;binx<=histBack->GetNbinsX();binx++) {
    for(int biny=1;biny<=histBack->GetNbinsY();biny++) {
      for(int binz=1;binz<=histBack->GetNbinsZ();binz++) {
	Double_t targetVal=histTarget->GetBinContent(binx,biny,binz);
	Double_t backVal=histBack->GetBinContent(binx,biny,binz);
	Double_t errorVal=sqrt(backVal);
	backVal*=numTarget;
	backVal/=numBackground;
	errorVal*=numTarget;
	errorVal/=numBackground;
	


	histDiff->SetBinContent(binx,biny,binz,targetVal-backVal);

	  
	
	if((targetVal-backVal)>numSigma*errorVal && backVal>0) {
	  histDiffPos->SetBinContent(binx,biny,binz,(targetVal-backVal)/errorVal);
	  histDiffPosXY->SetBinContent(binx,biny,histDiffPosXY->GetBinContent(binx,biny)+(targetVal-backVal)/errorVal);
	  histDiffPosXZ->SetBinContent(binx,binz,histDiffPosXZ->GetBinContent(binx,binz)+(targetVal-backVal)/errorVal);
	  histDiffPosYZ->SetBinContent(biny,binz,histDiffPosYZ->GetBinContent(biny,binz)+(targetVal-backVal)/errorVal);
	  
	  histDiffPosXYSlice[binz-1]->SetBinContent(binx,biny,histDiffPosXYSlice[binz-1]->GetBinContent(binx,biny)+(targetVal-backVal)/errorVal);

	  std::cout << targetVal << "\t" << backVal << "\t"
		    << (targetVal-backVal)/errorVal << "\n";
	}
	else if((backVal-targetVal)>numSigma*errorVal) {
	  histDiffNeg->SetBinContent(binx,biny,binz,-1*(targetVal-backVal));
	}


      }
    }
  }

  fpOut->Write();
//   TCanvas *can = new TCanvas();
//   histDiff->Draw();
//   TCanvas *can = new TCanvas();
//   histDiffPos->Draw();
  TCanvas *canProj = new TCanvas();
  canProj->Divide(1,3);
  canProj->cd(1);
  histDiffPosXY->Draw("colz");
  histDiffPosXY->SetStats(0);
  histDiffPosXY->GetXaxis()->SetTitle("X Position (mm)");
  histDiffPosXY->GetYaxis()->SetTitle("Y Position (mm)");
  canProj->cd(2);
  histDiffPosXZ->Draw("colz");
  histDiffPosXZ->SetStats(0);
  histDiffPosXZ->GetXaxis()->SetTitle("X Position (mm)");
  histDiffPosXZ->GetYaxis()->SetTitle("Z Position (mm)");
  canProj->cd(3);
  histDiffPosYZ->Draw("colz");
  histDiffPosYZ->SetStats(0);
  histDiffPosYZ->GetXaxis()->SetTitle("Y Position (mm)");
  histDiffPosYZ->GetYaxis()->SetTitle("Z Position (mm)");

//   TCanvas *can = new TCanvas();
//   histDiffNeg->Draw();
	

//   TCanvas *canSlice = new TCanvas("canSlice","canSlice");
//   canSlice->Divide(5,6);
//   for(int i=0;i<numBins;i++) {
//     canSlice->cd(i+1);
//     histDiffPosXYSlice[i]->Draw("colz");
//   }

}
