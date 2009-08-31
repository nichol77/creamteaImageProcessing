

void quickMakeHisto() {

  
  new TCanvas();
  TChain *nocontainerTree = new TChain("pcaTree");
  nocontainerTree->Add("/home/rjn/creamtea/data/nocontainer/pca_nocontainer_0_1000files.root");
  nocontainerTree->Add("/home/rjn/creamtea/data/nocontainer/pca_nocontainer_1000_2000files.root");
  nocontainerTree->Add("/home/rjn/creamtea/data/nocontainer/pca_nocontainer_3000_4000files.root");

  TFile *fpNocontainer = new TFile("nocontainer.root","RECREATE");
  TH3F *histNocontainer = new TH3F("histNoContainer","histNoContainer",100,-6000,+6000,100,-6000,+6000,100,-6000,+6000);
  
  nocontainerTree->Draw("xPosTrue:yPosTrue:zPosTrue>>histNoContainer","abs(zPosTrue)<3000 && abs(xPosTrue)<6000 && abs(yPosTrue)<6000 && thetaTrue>0.00","");
  histNocontainer->Write();


  new TCanvas();
  TChain *containerTree = new TChain("pcaTree");
  containerTree->Add("/home/rjn/creamtea/data/container13m/pca_container_1_1000files.root");
  containerTree->Add("/home/rjn/creamtea/data/container13m/pca_container_1001_2000files.root");
  containerTree->Add("/home/rjn/creamtea/data/container13m/pca_container_2001_3000files.root");

  TFile *fpcontainer = new TFile("container.root","RECREATE");
  TH3F *histcontainer = new TH3F("histContainer","histContainer",100,-6000,+6000,100,-6000,+6000,100,-6000,+6000);
  
  containerTree->Draw("xPosTrue:yPosTrue:zPosTrue>>histContainer","abs(zPosTrue)<3000 && abs(xPosTrue)<6000 && abs(yPosTrue)<6000 && thetaTrue>0.00","");
  histcontainer->Write();


  new TCanvas();
  TChain *targetTree = new TChain("pcaTree");
  targetTree->Add("/home/rjn/creamtea/data/target/pca_container_1_1000files.root");
  targetTree->Add("/home/rjn/creamtea/data/target/pca_container_1001_2000files.root");
  targetTree->Add("/home/rjn/creamtea/data/target/pca_container_2001_3000files.root");

  TFile *fptarget = new TFile("target.root","RECREATE");
  TH3F *histtarget = new TH3F("histTarget","histTarget",100,-6000,+6000,100,-6000,+6000,100,-6000,+6000);
  
  targetTree->Draw("xPosTrue:yPosTrue:zPosTrue>>histTarget","abs(zPosTrue)<3000 && abs(xPosTrue)<6000 && abs(yPosTrue)<6000 && thetaTrue>0.00","");
  histtarget->Write();



}
