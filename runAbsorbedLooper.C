

void runAbsorbedLooper() 
{
  gSystem->CompileMacro("AbsorbedLooper.C","k");

  TChain *Absorbed = new TChain("Absorbed");
  //  Absorbed->Add("/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm_1_1000.root");
  //  Absorbed->Add("/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm_1001_2000.root");
  //  Absorbed->Add("/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm_2001_3000.root");
  
  char inputName[180];
  //Local file old name system
  //  for(int i=0;i<=1000;i+=1000) {
    //    sprintf(inputName,"/unix/anita1/creamtea/strips_650/fakecontainer/pca/pca_fakecontainer_%d_%d.root",i+1,i+1000);
    //    sprintf(inputName,"/unix/anita1/creamtea/strips_650/container/pca/pca_container_%d_%d.root",i+1,i+1000);
    //    sprintf(inputName,"/unix/anita1/creamtea/strips_650/container_10cmtungsten/pca/pca_tungsten10cm_%d_%d.root",i+1,i+1000);
    //    sprintf(inputName,"/unix/anita1/creamtea/strips_650/fakecontainer_10cmtargetat_1_3_1/pca/pca_fakecontainer_10cmtarget_%d_%d.root",i+1,i+1000);

  //  }
//Downloaded from UCL
for(int tag=1;tag<=1;tag++) {
  //  sprintf(inputName,"http://www.hep.ucl.ac.uk/~rjn/creamtea/pcaFiles/fakecontainer_10cmtargetat_0p5_1_0p5_steelboxat_m0p5_3_m0p5/pca_fakecontainer_10cmtarget_steelbox_million_%d.root",tag);
  //  sprintf(inputName,"/home/rjn/creamtea/data/strips_650/fakecontainer_10cmtargetat_0p5_1_0p5_steelboxat_m0p5_3_m0p5/pca_fakecontainer_10cmtarget_steelbox_million_%d.root",tag);
  //  sprintf(inputName,"/unix/anita1/creamtea/strips_650/fakecontainer_10cmtargetat_0p5_1_0p5_hollowsteelboxat_m0p5_3_m0p5/pca/pca_fakecontainer_10cmtarget_hollowsteelbox_million_%d.root",tag);
   sprintf(inputName,"/unix/anita1/creamtea/minerva/fakecontainer_10cmtarget/pca/pca_fakecontainer_10cmtarget_million_%d.root",tag);
   //   sprintf(inputName,"/unix/anita1/creamtea/minerva/fakecontainer_5cmtarget/pca/pca_fakecontainer_5cmtarget_million_%d.root",tag);
   //  sprintf(inputName,"/unix/anita1/creamtea/minerva/fakecontainer_notarget/pca/pca_fakecontainer_notarget_million_%d.root",tag);
    Absorbed->Add(inputName);
 }
 
 AbsorbedLooper theLooper(Absorbed);
 // theLooper.MakeSliceHists(480);
 // theLooper.MakeSliceHistsIteratively(100);
 theLooper.MakeSliceHistsIterativelyReco(100);
 
  

}
