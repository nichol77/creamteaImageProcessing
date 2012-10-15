#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "/unix/creamtea/sfayer/1mdetector_5cmtarget/include/DetectorDefs.hh"

using namespace std;

double getMean (TFile* in_t,int nobinsx,int nobinsy,int nobinsz);

double getDeviation (TFile* in_t, double mean_t,int nobinsx,int nobinsy,int nobinsz);

double getSignificance (TFile* in_t, TFile* out_t, double mean_t, double deviation_t,int nobinsx,int nobinsy,int nobinsz);

double targetfind (char *BaseFile,int seconds, TFile* in_t,int sig_t,int noBinsX,int noBinsY,int noBinsZ);

double backgroundSubtract (TFile* in_t, int seconds_t,int noBinsx,int noBinsy,int noBinsz) {
char BackName[180];
char SliceName[20];
char *backname="~/imageProcessing/trunk/iterations/Lambda_40iterations_blank_target";
sprintf(BackName,"%s_%d_seconds.root",backname,seconds_t);
TFile* background = new TFile(BackName);
double lambdaT, lambdaB;
double meanB, deviationB;
double meanT, deviationT;
TFile *flambdas = new TFile("~/imageProcessing/trunk/tmp/Lambdas.root","RECREATE");
   TTree *lambdasTree = new TTree("lambdasTree","lambdasTree");
   lambdasTree->Branch("lambdaT",&lambdaT,"lambdaT/I");
lambdasTree->Branch("lambdaB",&lambdaB,"lambdaB/I");
for(int SliceNo = 1; SliceNo < noBinsz; SliceNo++) {
sprintf(SliceName,"Slice_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) in_t->Get(SliceName);
 TH2D* BackgroundSlice = (TH2D*) background->Get(SliceName);
for (int x = 1; x < noBinsx; x++) {
	 for (int y = 1; y < noBinsy; y++) {
	   lambdaT = CurrentSlice->GetBinContent(x,y);
	lambdaB = BackgroundSlice->GetBinContent(x,y);
	lambdasTree->Fill();
	 }
 }
 }
lambdasTree->AutoSave();
std::cout << "Background at " << seconds_t << ": " << endl;
meanB = getMean(background,noBinsx,noBinsy,noBinsz);
deviationB = getDeviation(background,meanB,noBinsx,noBinsy,noBinsz);
std::cout << "Whole image at " << seconds_t << ": " << endl;
meanT = getMean(in_t,noBinsx,noBinsy,noBinsz);
deviationT = getDeviation(in_t,meanT,noBinsx,noBinsy,noBinsz);
flambdas->Close();
return ((deviationT-deviationB)/deviationB*100);
}

double pictureMap (char *BaseFile,int seconds,TFile* in_t,int noBinsx,int noBinsy,int noBinsz){
	TObjArray *HlistPic = new TObjArray(100);
	char HistName[80];
	char HistoPic[20];
	char HistoPicAbove[20];
	char HistoPicTwoAbove[20];
TFile* out_t = new TFile("~/imageProcessing/trunk/tmp/image_Seen.root", "RECREATE");
out_t->Write();
TCanvas *canProj = new TCanvas();
      canProj->Divide(3,3);
 double currentLambda = 0;
double significance;
double imagemean;
 double deviation;
imagemean = getMean(in_t,noBinsx,noBinsy,noBinsz);
deviation = getDeviation(in_t,imagemean,noBinsx,noBinsy,noBinsz);
std::cout<<"Looking for targets"<<endl;
//for(int sig=5;sig>0;sig--){
int sig=1;
for(int SliceNo = 1; SliceNo < noBinsz; SliceNo++){
 sprintf(HistoPic, "%s_%d", "Image_Seen", SliceNo);
 TH2F *Pic = new TH2F(HistoPic, HistoPic,noBinsx,0.,1., noBinsy,0.,1.);
 HlistPic->Add(Pic);
 sprintf(HistName,"Slice_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) in_t->Get(HistName);
 for (int x = 1; x < noBinsx; x++) {
	 for (int y = 1; y < noBinsy; y++) {
	 currentLambda = CurrentSlice->GetBinContent(x,y);
	 if ((currentLambda >= (0.03*sig))&&(Pic->GetBinContent(x,y)<1)){ 
		 Pic->SetBinContent(x,y,100.);
	 } else {
if(Pic->GetBinContent(x,y)<1) Pic->SetBinContent(x,y,0.);
}
}
}
canProj->cd(SliceNo);
CurrentSlice->Draw("colz");
	 CurrentSlice->SetStats(0);
	 CurrentSlice->GetXaxis()->SetTitle("X ID");
	 CurrentSlice->GetYaxis()->SetTitle("Y ID");
}
 HlistPic->Write();

// remove isolated peaks
for(int SliceNo = 1; SliceNo <= (noBinsz-3); SliceNo++){
	sprintf(HistoPic, "%s_%d", "Image_Seen", SliceNo);
TH2D* ImageSlice = (TH2D*) out_t->Get(HistoPic);
	sprintf(HistoPicAbove, "%s_%d", "Image_Seen", (SliceNo+1));
TH2D* ImageSliceAbove = (TH2D*) out_t->Get(HistoPicAbove);
sprintf(HistoPicTwoAbove, "%s_%d", "Image_Seen", (SliceNo+2));
TH2D* ImageSliceTwoAbove = (TH2D*) out_t->Get(HistoPicTwoAbove);
for (int x = 1; x <= (noBinsx-3); x++) {
  for (int y = 1; y <= (noBinsy-3); y++) {
		 if (ImageSliceAbove->GetBinContent((x+1),(y+1))==100.) {
		 if ((ImageSlice->GetBinContent((x+1),(y+1))==0.)&&(ImageSliceTwoAbove->GetBinContent((x+1),(y+1))==0.)&&(ImageSliceAbove->GetBinContent((x+1),y)==0.)&&(ImageSliceAbove->GetBinContent((x+1),(y+2))==0.)&&(ImageSliceAbove->GetBinContent(x,(y+1))==0.)&&(ImageSliceAbove->GetBinContent((x+2),(y+1))==0.)){
			ImageSliceAbove->SetBinContent((x+1),(y+1), 0.);
		 }
	 }
	 }
}
}

//now calculate signifcance of highlighted voxels
significance = getSignificance(in_t,out_t,imagemean,deviation,noBinsx,noBinsy,noBinsz);
//group peaks to targets
targetfind(BaseFile,seconds,out_t,sig,noBinsx,noBinsy,noBinsz);
HlistPic->Write();
//}
out_t->Close();
return significance;
}

double getMean (TFile* in_t,int nobinsx,int nobinsy,int nobinsz) {
double mean;
char SliceName[20];
for(int SliceNo = 1; SliceNo < nobinsz; SliceNo++) {
sprintf(SliceName,"Slice_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) in_t->Get(SliceName);
for (int x = 1; x < nobinsx; x++) {
	 for (int y = 1; y < nobinsy; y++) {
	   mean += CurrentSlice->GetBinContent(x,y);
	 }
	 mean /= nobinsy;
 }
 mean /= nobinsx;
 }
mean /= nobinsz;
std::cout << "Mean of image: " << mean << endl;
return mean;
}

double getDeviation (TFile* in_t, double mean_t,int nobinsx,int nobinsy,int nobinsz){
double stddeviation;
char SliceName[20];
for(int SliceNo = 1; SliceNo < nobinsz; SliceNo++){
sprintf(SliceName,"Slice_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) in_t->Get(SliceName);
 for (int x = 1; x < nobinsx; x++) {
	 for (int y = 1; y < nobinsy; y++) {
		 stddeviation += ((CurrentSlice->GetBinContent(x,y))-mean_t)*((CurrentSlice->GetBinContent(x,y))-mean_t);
	 }
 }
}
stddeviation /= (nobinsx*nobinsy*nobinsz);
stddeviation = sqrt (stddeviation);
std::cout << "Standard deviation: " << stddeviation << endl;
return stddeviation;
}

double getSignificance (TFile* in_t, TFile* out_t, double mean_t, double deviation_t, int nobinsx, int nobinsy, int nobinsz){
double calcsignificance = 0;
int voxels = 1;
char HistoPic[20];
char SliceName[20];
for(int SliceNo = 1; SliceNo < nobinsz; SliceNo++){
	sprintf(HistoPic, "%s_%d", "Image_Seen", SliceNo);
	TH2D* ImageSlice = (TH2D*) out_t->Get(HistoPic);
	sprintf(SliceName,"Slice_%d",SliceNo);
	TH2D* CurrentSlice = (TH2D*) in_t->Get(SliceName);
	for (int x = 1; x < nobinsx; x++) {
	 for (int y = 1; y < nobinsy; y++) {
		 if (ImageSlice->GetBinContent(x,y)==1.){
			 calcsignificance = (CurrentSlice->GetBinContent(x,y))-mean_t;
			 voxels++;
		 }
	 }
	}
}
	calcsignificance /= voxels;
	calcsignificance /= deviation_t;
return calcsignificance;
}

// now attempt to improve image quality
// plug holes
int plugholes(int flag,TFile* out_t,int noBinsx,int noBinsy,int noBinsz) {
int added=0;
char HistoPic[20];
char HistoPicAbove[20];
char HistoPicTwoAbove[20];
for(int SliceNo = 1; SliceNo <= (noBinsz-1); SliceNo++){
sprintf(HistoPic, "%s_%d", "Image_Seen", SliceNo);
TH2D* ImageSlice = (TH2D*) out_t->Get(HistoPic);
 for (int x = 1; x < noBinsx; x++) {
   for (int y = 1; y < (noBinsy-2); y++) {
	   if ((ImageSlice->GetBinContent(x,y)==(flag))&&(ImageSlice->GetBinContent(x,(y+2))==(flag))) {
			 ImageSlice->SetBinContent(x,(y+1),(flag));
			added++;
		 }
	 }
 }
 for (int y = 0; y < noBinsy; y++) {
	 for (int x = 0; x < (noBinsx-2); x++) {
	   if ((ImageSlice->GetBinContent(x,y)==(flag))&&(ImageSlice->GetBinContent((x+2),y)==(flag))) {
			 ImageSlice->SetBinContent((x+1),y,(flag));
			added++;
		 }
	 }
 }
 if (SliceNo <= (noBinsz-3)){
sprintf(HistoPicAbove, "%s_%d", "Image_Seen", (SliceNo+1));
TH2D* ImageSliceAbove = (TH2D*) out_t->Get(HistoPicAbove);
sprintf(HistoPicTwoAbove, "%s_%d", "Image_Seen", (SliceNo+2));
TH2D* ImageSliceTwoAbove = (TH2D*) out_t->Get(HistoPicTwoAbove);
for (int y = 1; y < noBinsy; y++) {
	 for (int x = 1; x < noBinsx; x++) {
	   if ((ImageSlice->GetBinContent(x,y)==(flag))&&(ImageSliceTwoAbove->GetBinContent(x,y)==(flag))) {
			 ImageSliceAbove->SetBinContent(x,y,(flag));
			added++;
		 }
	 }
}
}
}
return added;
}

int expandtarget (int flag,TFile* out_t,int noBinsx,int noBinsy,int noBinsz){
//std::cout<<"finding target " << flag << " bounds"<<endl;
  int added=1;
  int total=0;
char SliceName[20];
char HistoPicAbove[20];
char HistoPicBelow[20];
  for(int iter=0;iter<=100;iter++){
for(int SliceNo = 2; SliceNo < noBinsz; SliceNo++){
sprintf(SliceName,"Image_Seen_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) out_t->Get(SliceName);
sprintf(HistoPicAbove, "%s_%d", "Image_Seen", (SliceNo+1));
TH2D* ImageSliceAbove = (TH2D*) out_t->Get(HistoPicAbove);
sprintf(HistoPicBelow, "%s_%d", "Image_Seen", (SliceNo-1));
TH2D* ImageSliceBelow = (TH2D*) out_t->Get(HistoPicBelow);
 for (int x = 2; x < noBinsx; x++) {
	 for (int y = 2; y < noBinsy; y++) {
	   if (CurrentSlice->GetBinContent(x,y)==(flag)) {
	     if(CurrentSlice->GetBinContent((x+1),y)==(100.)) {
	       CurrentSlice->SetBinContent((x+1),y,(flag));
					   added++;
					   }
	     if(CurrentSlice->GetBinContent((x-1),y)==(100.)) {
	       CurrentSlice->SetBinContent((x-1),y,(flag));
					   added++;
					   }
		 if(CurrentSlice->GetBinContent(x,(y+1))==(100.)) {
		   CurrentSlice->SetBinContent(x,(y+1),(flag));
			added++;
		 }
		 if(CurrentSlice->GetBinContent(x,(y-1))==(100.)) {
		   CurrentSlice->SetBinContent(x,(y-1),(flag));
			added++;
		 }
	     if(ImageSliceAbove->GetBinContent(x,y)==(100.)) {
		   ImageSliceAbove->SetBinContent(x,y,(flag));
			added++;
		 }
	     if(ImageSliceBelow->GetBinContent(x,y)==(100.)) {
	       ImageSliceBelow->SetBinContent(x,y,(flag));
			added++;
		 }
	   }
	 }
 }
	     }
 total+=added;
//cout<<"total" <<total<<endl;
 if(added==0) break;
 added=0;
  }
	     return total;
 }

double targetfind (char *BaseFile,int seconds, TFile* out,int sig_t,int noBinsX,int noBinsY,int noBinsZ){
double avelambda; 
int flags[100]= {0,0,0,0,0,0,0,0,0,0};
int a = 1;
   int surrounding;
   char SliceName[20];
char HistoPicAbove[20];
char HistoPicBelow[20];
for(int SliceNo = 2; SliceNo < (noBinsZ-1); SliceNo++){
sprintf(SliceName,"Image_Seen_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) out->Get(SliceName);
sprintf(HistoPicAbove, "%s_%d", "Image_Seen", (SliceNo+1));
TH2D* ImageSliceAbove = (TH2D*) out->Get(HistoPicAbove);
sprintf(HistoPicBelow, "%s_%d", "Image_Seen", (SliceNo-1));
TH2D* ImageSliceBelow = (TH2D*) out->Get(HistoPicBelow);
for (int x = 2; x < (noBinsX-1); x++) {
	 for (int y = 2; y < (noBinsY-1); y++) {
	   if(CurrentSlice->GetBinContent(x,y)==(100.)){
	     surrounding=0;
	     if(CurrentSlice->GetBinContent((x+1),y)==(100.)) surrounding++;
	     if(CurrentSlice->GetBinContent((x-1),y)==(100.)) surrounding++;
	     if(CurrentSlice->GetBinContent(x,(y+1))==(100.)) surrounding++;
	     if(CurrentSlice->GetBinContent(x,(y-1))==(100.)) surrounding++;
	     if(ImageSliceAbove->GetBinContent(x,y)==(100.)) surrounding++;
	     if(ImageSliceBelow->GetBinContent(x,y)==(100.)) surrounding++;
	     if (surrounding >= 6) {
		//std::cout<<"target found!"<<endl;
	       CurrentSlice->SetBinContent(x,y,((a+((sig_t-1)*20.))));
	       flags[(a+((sig_t-1)*20))]=expandtarget((a+((sig_t-1)*20)),out,(noBinsX-1),(noBinsY-1),(noBinsZ-1));
		flags[(a+((sig_t-1)*20))]+=plugholes((a+((sig_t-1)*20)),out,noBinsX,noBinsY,noBinsZ);
	       a++;
		if(a>=20) goto output;
	   }
	 }
	if(SliceNo>=(noBinsZ-2)&&y>=(noBinsY-2)&&x>=(noBinsX-2)){
std::cout<<"finished"<< endl;
goto output;
 }
 }
 }
   }
output:
for (int b=(((sig_t-1)*20)+1);b<(sig_t*20);b++){
if(flags[b]>0) {
std::cout << "Target " << b << " has size " << flags[b] << " hits" << endl;
}
}
if(sig_t==1){
char TargetFile[180];
double content;
sprintf(TargetFile,"~/imageProcessing/trunk/targets/%s_%d_seconds.root",BaseFile,seconds);
TFile* targets = new TFile(TargetFile, "RECREATE");
TObjArray *HistoList = new TObjArray(100);
	char HistoTarget[20];
for(int SliceNo = 1; SliceNo < noBinsZ; SliceNo++){
 sprintf(HistoTarget, "%s_%d", "Section", SliceNo);
 TH2F *target = new TH2F(HistoTarget, HistoTarget,noBinsX,0.,1., noBinsY,0.,1.);
TH2D* Section = (TH2D*) out->Get(HistoTarget);
for (int x = 1; x < noBinsX; x++) {
	 for (int y = 1; y < noBinsY; y++) {
content = Section->GetBinContent(x,y);
target->SetBinContent(x,y,content);
}
}
 HistoList->Add(target);
target->Draw("colz");
	 target->SetStats(0);
	 target->GetXaxis()->SetTitle("X ID");
	 target->GetYaxis()->SetTitle("Y ID");
}
HistoList->Write();
targets->Close();
}
avelambda=0.;
return avelambda;
}

void makethree(char *BaseFile,int seconds, int noBinsX,int noBinsY,int noBinsZ){
char TargetFile[180];
double dummy;
sprintf(TargetFile,"~/imageProcessing/trunk/targets/%s_%d_seconds.root",BaseFile,seconds);
TFile* targets = new TFile(TargetFile);
TFile* image3D = new TFile("~/imageProcessing/trunk/3Dimage.root","RECREATE");
char *threeN = "Targets_3D";
TH3F *three = new TH3F(threeN,threeN,noBinsX,0.,1., noBinsY,0.,1.,noBinsZ,0.,1.);
char HistoTarget[20];
for(int SliceNo = 1; SliceNo < noBinsZ; SliceNo++){
 sprintf(HistoTarget, "%s_%d", "Image_Seen", SliceNo);
TH2D* CurrentSlice = (TH2D*) targets->Get(HistoTarget);
for (int x = 1; x < noBinsX; x++) {
	 for (int y = 1; y < noBinsY; y++) {
	dummy=CurrentSlice->GetBinContent(x,y);
if((dummy>1)&&(dummy!=10))three->SetBinContent(x,y,SliceNo,dummy);
}
}
}
three->Write();
targets->Close();
}

int main (int argc, char**argv){
char inputFile[180];
char BaseFile[180];
int seconds = 60;
int start, finish;
//char *input="~/imageProcessing/trunk/iterations/Lambda_40iterations_10cmtarget";
if(argc<3) {
      std::cerr << "Usage:\t" << "<input file location> <base name> <seconds(optional)>\n";
      return -1;
   }
strncpy(inputFile,argv[1],180);
strncpy(BaseFile,argv[2],180);
if(argc==4){
start=atoi(argv[3]);
finish=atoi(argv[3]);
}
if(argc>4){
start=atoi(argv[3]);
finish=atoi(argv[4]);
}
for(seconds=start;seconds<=finish;seconds++){
if(seconds==1||seconds%5==0){
char FileName[180];
int noBinsX = NX;
int noBinsY = NY;
int noBinsZ = NZ;
sprintf(FileName,"%s%s_%d_seconds.root",inputFile, BaseFile,seconds);
TFile* in = new TFile(FileName);
double y;
bool bo;
y = backgroundSubtract(in,seconds,noBinsX,noBinsY,noBinsZ);
std::cout << "Percentage difference in deviations is " << y << "" << endl;
/*std::cout << "Would you like to map image? (1/0)" << endl;
std::cin >> bo;
if(bo==1){*/
double z;
std::cout << "Now mapping significant lambdas" << endl;
  z = pictureMap(BaseFile,seconds,in,noBinsX,noBinsY,noBinsZ);
std::cout << "Image significance compared to image average is " << z << " sigma" << endl;
/*std::cout << "Would you like to make 3D image? (1/0)" << endl;
std::cin >> bo;
if(bo==1){
makethree(BaseFile,seconds,noBinsX,noBinsY,noBinsZ);
}
}*/
}
}
return 0;
}


