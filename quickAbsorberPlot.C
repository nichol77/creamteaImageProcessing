

void quickAbsorberPlot() 
{
  TChain *Absorbed = new TChain("Absorbed");
  //  Absorbed->Add("/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm_1_1000.root");
  //  Absorbed->Add("/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm_1001_2000.root");
  //  Absorbed->Add("/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm_2001_3000.root");
  
  char inputName[180];
  for(int i=0;i<=9000;i+=1000) {
    //    sprintf(inputName,"/unix/anita1/creamtea/strips_650/fakecontainer/pca/pca_fakecontainer_%d_%d.root",i+1,i+1000);
    //    sprintf(inputName,"/unix/anita1/creamtea/strips_650/container/pca/pca_container_%d_%d.root",i+1,i+1000);
    //    sprintf(inputName,"/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm_%d_%d.root",i+1,i+1000);
    sprintf(inputName,"/unix/anita1/creamtea/strips_650/fakecontainer_10cmtargetat_1_3_1/pca/pca_fakecontainer_10cmtarget_%d_%d.root",i+1,i+1000);
    Absorbed->Add(inputName);
  }
  
  
  TCanvas *can = new TCanvas("can","can");
  can->Divide(2,2);
  can->cd(1);
  Absorbed->Draw("yGrad*0000+yCut:xGrad*0000+xCut","abs(xGrad*-7000+xCut)<5000 && abs(yGrad*-7000+yCut)<5000","colz");
  can->cd(2);
  Absorbed->Draw("yGrad*1000+yCut:xGrad*1000+xCut","abs(xGrad*-7000+xCut)<5000 && abs(yGrad*-7000+yCut)<5000","colz");
  can->cd(3);
  Absorbed->Draw("yGrad*2000+yCut:xGrad*2000+xCut","abs(xGrad*-7000+xCut)<5000 && abs(yGrad*-7000+yCut)<5000","colz");
  can->cd(4);
  Absorbed->Draw("yGrad*3000+yCut:xGrad*3000+xCut","abs(xGrad*-7000+xCut)<5000 && abs(yGrad*-7000+yCut)<5000","colz");
  

  TCanvas *can2 = new TCanvas("can2","can2");
  can2->Divide(2,2);
  can2->cd(1);
  Absorbed->Draw("yGrad*0000+yCut:xGrad*0000+xCut","abs(xGrad*-7000+xCut)<5000 && abs(yGrad*-7000+yCut)<5000","colz");
  can2->cd(2);
  Absorbed->Draw("yGrad*-1000+yCut:xGrad*-1000+xCut","abs(xGrad*-7000+xCut)<5000 && abs(yGrad*-7000+yCut)<5000","colz");
  can2->cd(3);
  Absorbed->Draw("yGrad*-2000+yCut:xGrad*-2000+xCut","abs(xGrad*-7000+xCut)<5000 && abs(yGrad*-7000+yCut)<5000","colz");
  can2->cd(4);
  Absorbed->Draw("yGrad*-3000+yCut:xGrad*-3000+xCut","abs(xGrad*-7000+xCut)<5000 && abs(yGrad*-7000+yCut)<5000","colz");
  


}
