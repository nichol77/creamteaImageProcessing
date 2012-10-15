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

#define LAMBDA_AIR (1./3.039e4)

using namespace std;

double getMean (TFile* in_t,int nobinsx,int nobinsy,int nobinsz);

double getDeviation (TFile* in_t, double mean_t,int nobinsx,int nobinsy,int nobinsz);

double getSignificance (TFile* in_t, TFile* out_t, double mean_t, double deviation_t,int nobinsx,int nobinsy,int nobinsz);

double targetfind (char *BaseFile,char *OutFile,int seconds, TFile* in, TFile* out,int noBinsX,int noBinsY,int noBinsZ);

double backgroundSubtract (TFile* in_t, int seconds_t,int noBinsx,int noBinsy,int noBinsz) {
char BackName[180];
char SliceName[20];
char *backname="/unix/creamtea/sfayer/blank_target/lambdamap/Lambda_40iterations_blank_target";
sprintf(BackName,"%s_%d_seconds.root",backname,seconds_t);
TFile* background = new TFile(BackName);
double lambdaT, lambdaB;
double meanB, deviationB;
double meanT, deviationT;
TFile *flambdas = new TFile("/unix/creamtea/sfayer/temp/Lambdas.root","RECREATE");
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
deviationB = getDeviation(background,LAMBDA_AIR,noBinsx,noBinsy,noBinsz);
std::cout << "Whole image at " << seconds_t << ": " << endl;
meanT = getMean(in_t,noBinsx,noBinsy,noBinsz);
deviationT = getDeviation(in_t,LAMBDA_AIR,noBinsx,noBinsy,noBinsz);
flambdas->Close();
background->Close();
return ((deviationT-deviationB)/deviationB*100);
}

// remove isolated peaks
void removeiso (TFile* out_t,int noBinsx,int noBinsy,int noBinsz){
	char HistoPic[20];
	char HistoPicAbove[20];
	char HistoPicTwoAbove[20];
for(int SliceNo = 1; SliceNo <= (noBinsz-3); SliceNo++){
	sprintf(HistoPic, "%s_%d", "Image_Seen", SliceNo);
TH2D* ImageSlice = (TH2D*) out_t->Get(HistoPic);
	sprintf(HistoPicAbove, "%s_%d", "Image_Seen", (SliceNo+1));
TH2D* ImageSliceAbove = (TH2D*) out_t->Get(HistoPicAbove);
sprintf(HistoPicTwoAbove, "%s_%d", "Image_Seen", (SliceNo+2));
TH2D* ImageSliceTwoAbove = (TH2D*) out_t->Get(HistoPicTwoAbove);
for (int x = 1; x <= (noBinsx-3); x++) {
  for (int y = 1; y <= (noBinsy-3); y++) {
		 if (ImageSliceAbove->GetBinContent((x+1),(y+1))==-1.) {
		 if ((ImageSlice->GetBinContent((x+1),(y+1))==0.)&&(ImageSliceTwoAbove->GetBinContent((x+1),(y+1))==0.)&&(ImageSliceAbove->GetBinContent((x+1),y)==0.)&&(ImageSliceAbove->GetBinContent((x+1),(y+2))==0.)&&(ImageSliceAbove->GetBinContent(x,(y+1))==0.)&&(ImageSliceAbove->GetBinContent((x+2),(y+1))==0.)){
			ImageSliceAbove->SetBinContent((x+1),(y+1), 0.);
		 }
	 }
	 }
}
}
}

double pictureMap (char *BaseFile,char *OutFile,int seconds,TFile* in_t,int noBinsx,int noBinsy,int noBinsz){
	TObjArray *HlistPic = new TObjArray(100);
	char HistName[80];
	char HistoPic[20];
	char HistoPicAbove[20];
	char HistoPicTwoAbove[20];
TFile* out_t = new TFile("/unix/creamtea/sfayer/temp/image_Seen.root", "RECREATE");
out_t->Write();
TCanvas *canProj = new TCanvas();
      canProj->Divide(3,3);
 double currentLambda = 0;
double significance;
double imagemean;
 double deviation;
imagemean = getMean(in_t,noBinsx,noBinsy,noBinsz);
deviation = getDeviation(in_t,LAMBDA_AIR,noBinsx,noBinsy,noBinsz);
std::cout<< "Looking for lambda greater than " << (LAMBDA_AIR+(2.5*deviation)) <<endl;
std::cout<<"Looking for targets"<<endl;
for(int SliceNo = 1; SliceNo < noBinsz; SliceNo++){
 sprintf(HistoPic, "%s_%d", "Image_Seen", SliceNo);
 TH2F *Pic = new TH2F(HistoPic, HistoPic,noBinsx,0.,1., noBinsy,0.,1.);
 HlistPic->Add(Pic);
 sprintf(HistName,"Slice_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) in_t->Get(HistName);
 for (int x = 1; x < noBinsx; x++) {
	 for (int y = 1; y < noBinsy; y++) {
	 currentLambda = CurrentSlice->GetBinContent(x,y);
	 /*if ((LAMBDA_AIR+(2.5*deviation))>0.04){
	if((currentLambda >= (LAMBDA_AIR+(2.5*deviation)))&&(Pic->GetBinContent(x,y)<1)){ 
		 Pic->SetBinContent(x,y,-1.);
	 } else {
if(Pic->GetBinContent(x,y)<1) Pic->SetBinContent(x,y,0.);
}
} else {*/
if((currentLambda >= (0.008))&&(Pic->GetBinContent(x,y)<1)){ 
		 Pic->SetBinContent(x,y,-1.);
	 } else {
if(Pic->GetBinContent(x,y)<1) Pic->SetBinContent(x,y,0.);
}
//}
}
}
canProj->cd(SliceNo);
CurrentSlice->Draw("colz");
	 CurrentSlice->SetStats(0);
	 CurrentSlice->GetXaxis()->SetTitle("X ID");
	 CurrentSlice->GetYaxis()->SetTitle("Y ID");
}

removeiso(out_t,noBinsx,noBinsy,noBinsz);

HlistPic->Write();

//now calculate signifcance of highlighted voxels
significance = getSignificance(in_t,out_t,LAMBDA_AIR,deviation,noBinsx,noBinsy,noBinsz);
//group peaks to targets
targetfind(BaseFile,OutFile,seconds,in_t,out_t,noBinsx,noBinsy,noBinsz);
//HlistPic->Write();
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
		 if (ImageSlice->GetBinContent(x,y)==-1.){
			 calcsignificance += (CurrentSlice->GetBinContent(x,y))-mean_t;
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
int plugholes(int flag,double lambdaVal,TFile* in_t,TFile* out_t,int noBinsx,int noBinsy,int noBinsz) {
int added=0;
char HistoPic[20];
char LambdaSlice[20];
char HistoPicAbove[20];
char HistoPicTwoAbove[20];
for(int SliceNo = 1; SliceNo <= (noBinsz-1); SliceNo++){
sprintf(LambdaSlice,"Slice_%d",SliceNo);
 TH2D* SliceLambdas = (TH2D*) in_t->Get(LambdaSlice);
sprintf(HistoPic, "%s_%d", "Image_Seen", SliceNo);
TH2D* ImageSlice = (TH2D*) out_t->Get(HistoPic);
 for (int x = 1; x < noBinsx; x++) {
   for (int y = 1; y < (noBinsy-2); y++) {
	   if ((ImageSlice->GetBinContent(x,y)==(flag))&&(ImageSlice->GetBinContent(x,(y+2))==(flag))) {
			 ImageSlice->SetBinContent(x,(y+1),(flag));
			lambdaVal+=SliceLambdas->GetBinContent(x,y);
			added++;
		 }
	 }
 }
 for (int y = 0; y < noBinsy; y++) {
	 for (int x = 0; x < (noBinsx-2); x++) {
	   if ((ImageSlice->GetBinContent(x,y)==(flag))&&(ImageSlice->GetBinContent((x+2),y)==(flag))) {
			 ImageSlice->SetBinContent((x+1),y,(flag));
			lambdaVal+=SliceLambdas->GetBinContent(x,y);
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
			lambdaVal+=SliceLambdas->GetBinContent(x,y);
			added++;
		 }
	 }
}
}
}
return added;
}

int expandtarget (int flag,double lambdaVal,TFile* in_t,TFile* out_t,int noBinsx,int noBinsy,int noBinsz){
//std::cout<<"finding target " << flag << " bounds"<<endl;
  int added=1;
  int total=0;
char SliceName[20];
char LambdaSlice[20];
char HistoPicAbove[20];
char HistoPicBelow[20];
  for(int iter=0;iter<=100;iter++){
for(int SliceNo = 2; SliceNo < noBinsz; SliceNo++){
sprintf(LambdaSlice,"Slice_%d",SliceNo);
 TH2D* SliceLambdas = (TH2D*) in_t->Get(LambdaSlice);
sprintf(SliceName,"Image_Seen_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) out_t->Get(SliceName);
sprintf(HistoPicAbove, "%s_%d", "Image_Seen", (SliceNo+1));
TH2D* ImageSliceAbove = (TH2D*) out_t->Get(HistoPicAbove);
sprintf(HistoPicBelow, "%s_%d", "Image_Seen", (SliceNo-1));
TH2D* ImageSliceBelow = (TH2D*) out_t->Get(HistoPicBelow);
 for (int x = 2; x < noBinsx; x++) {
	 for (int y = 2; y < noBinsy; y++) {
	   if (CurrentSlice->GetBinContent(x,y)==(flag)) {
	     if(CurrentSlice->GetBinContent((x+1),y)==(-1.)) {
	       CurrentSlice->SetBinContent((x+1),y,(flag));
		lambdaVal+=SliceLambdas->GetBinContent(x,y);
					   added++;
					   }
	     if(CurrentSlice->GetBinContent((x-1),y)==(-1.)) {
	       CurrentSlice->SetBinContent((x-1),y,(flag));
		lambdaVal+=SliceLambdas->GetBinContent(x,y);
					   added++;
					   }
		 if(CurrentSlice->GetBinContent(x,(y+1))==(-1.)) {
		   CurrentSlice->SetBinContent(x,(y+1),(flag));
		lambdaVal+=SliceLambdas->GetBinContent(x,y);
			added++;
		 }
		 if(CurrentSlice->GetBinContent(x,(y-1))==(-1.)) {
		   CurrentSlice->SetBinContent(x,(y-1),(flag));
		lambdaVal+=SliceLambdas->GetBinContent(x,y);
			added++;
		 }
	     if(ImageSliceAbove->GetBinContent(x,y)==(-1.)) {
		   ImageSliceAbove->SetBinContent(x,y,(flag));
		lambdaVal+=SliceLambdas->GetBinContent(x,y);
			added++;
		 }
	     if(ImageSliceBelow->GetBinContent(x,y)==(-1.)) {
	       ImageSliceBelow->SetBinContent(x,y,(flag));
		lambdaVal+=SliceLambdas->GetBinContent(x,y);
			added++;
		 }
	   }
	 }
 }
	     }
 total+=added;
//cout<<"total" <<total<<endl;
 if(added==0) {break;}
 added=0;
  }
	     return total;
 }

double targetfind (char *BaseFile,char *OutFile,int seconds, TFile* in, TFile* out,int noBinsX,int noBinsY,int noBinsZ){
double avelambda; 
int numhits[100];
double lamhits[100];
int totaltarg;
int a = 1;
for (int initial=1;initial<100;initial++){
numhits[(initial)]=0;
lamhits[(initial)]=0;
}
   int surrounding;
double pixelLambda=0;
   char SliceName[20];
char LambdaSlice[20];
char HistoPicAbove[20];
char HistoPicBelow[20];
for(int SliceNo = 2; SliceNo < (noBinsZ-1); SliceNo++){
sprintf(LambdaSlice,"Slice_%d",SliceNo);
 TH2D* SliceLambdas = (TH2D*) in->Get(LambdaSlice);
sprintf(SliceName,"Image_Seen_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) out->Get(SliceName);
sprintf(HistoPicAbove, "%s_%d", "Image_Seen", (SliceNo+1));
TH2D* ImageSliceAbove = (TH2D*) out->Get(HistoPicAbove);
sprintf(HistoPicBelow, "%s_%d", "Image_Seen", (SliceNo-1));
TH2D* ImageSliceBelow = (TH2D*) out->Get(HistoPicBelow);
for (int x = 2; x < (noBinsX-1); x++) {
	 for (int y = 2; y < (noBinsY-1); y++) {
	   if(CurrentSlice->GetBinContent(x,y)==(-1.)){
	     surrounding=0;
	     if(CurrentSlice->GetBinContent((x+1),y)==(-1.)) surrounding++;
	     if(CurrentSlice->GetBinContent((x+1),(y+1))==(-1.)) surrounding++;
	     if(CurrentSlice->GetBinContent((x-1),(y-1))==(-1.)) surrounding++;
	     if(CurrentSlice->GetBinContent((x+1),(y-1))==(-1.)) surrounding++;
	     if(CurrentSlice->GetBinContent((x-1),(y+1))==(-1.)) surrounding++;
	     if(CurrentSlice->GetBinContent((x-1),y)==(-1.)) surrounding++;
	     if(CurrentSlice->GetBinContent(x,(y+1))==(-1.)) surrounding++;
	     if(CurrentSlice->GetBinContent(x,(y-1))==(-1.)) surrounding++;
	     if(ImageSliceAbove->GetBinContent(x,y)==(-1.)) surrounding++;
	     if(ImageSliceAbove->GetBinContent((x+1),y)==(-1.)) surrounding++;
	     if(ImageSliceAbove->GetBinContent((x-1),y)==(-1.)) surrounding++;
	     if(ImageSliceAbove->GetBinContent(x,(y+1))==(-1.)) surrounding++;
	     if(ImageSliceAbove->GetBinContent(x,(y-1))==(-1.)) surrounding++;
	     if(ImageSliceBelow->GetBinContent(x,y)==(-1.)) surrounding++;
	     if(ImageSliceBelow->GetBinContent((x+1),y)==(-1.)) surrounding++;
	     if(ImageSliceBelow->GetBinContent((x-1),y)==(-1.)) surrounding++;
	     if(ImageSliceBelow->GetBinContent(x,(y+1))==(-1.)) surrounding++;
	     if(ImageSliceBelow->GetBinContent(x,(y-1))==(-1.)) surrounding++;
	     if (surrounding >= 14) {
		//std::cout<<"target found!"<<endl;
		pixelLambda=SliceLambdas->GetBinContent(x,y);
	       CurrentSlice->SetBinContent(x,y,(a));
	       numhits[(a)]=expandtarget((a),pixelLambda,in,out,(noBinsX-1),(noBinsY-1),(noBinsZ-1)); 
		lamhits[(a)]=pixelLambda/numhits[(a)];
		numhits[(a)]+=plugholes((a),pixelLambda,in,out,noBinsX,noBinsY,noBinsZ);
	       a++;
		pixelLambda=0;
		if(a>=100) goto output;
	   }
	 }
	if(SliceNo>=(noBinsZ-2)&&y>=(noBinsY-2)&&x>=(noBinsX-2)){
totaltarg=a-1;
std::cout<<"finished"<< endl;
goto output;
 }
 }
 }
   }
output:
std::cout << "Number of targets found = " << totaltarg << endl;
for (int b=1;b<(100);b++){
if(isnan(numhits[b])==1){
numhits[(b)]=0;
lamhits[(b)]=0;
}
if(numhits[b]>0) {
avelambda+=lamhits[b];
std::cout << "Target " << b << " has size " << numhits[b] << " hits, with mean lambda " << lamhits[b] << endl;
}
}
avelambda/=totaltarg;
char TargetFile[180];
double content;
sprintf(TargetFile,"/unix/creamtea/sfayer/%s/targets/%s_%d_seconds.root",OutFile,BaseFile,seconds);
TFile* targets = new TFile(TargetFile, "RECREATE");
TObjArray *HistoList = new TObjArray(100);
	char HistoTarget[20];
	char HistoRef[20];
for(int SliceNo = 1; SliceNo < noBinsZ; SliceNo++){
 sprintf(HistoTarget, "%s_%d", "Section", SliceNo);
sprintf(HistoRef, "%s_%d", "Image_Seen", SliceNo);
 TH2F *target = new TH2F(HistoTarget, HistoTarget,noBinsX,0.,1., noBinsY,0.,1.);
TH2D* CurrentSlice = (TH2D*) out->Get(HistoRef);
for (int x = 1; x < noBinsX; x++) {
	 for (int y = 1; y < noBinsY; y++) {
content = CurrentSlice->GetBinContent(x,y);
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
return avelambda;
}

void makethree(int targnum,char *BaseFile,char *OutFile,int seconds, int noBinsX,int noBinsY,int noBinsZ){
char TargetFile[180];
char ThreeFile[180];
double dummy;
sprintf(TargetFile,"/unix/creamtea/sfayer/%s/targets/%s_%d_seconds.root",OutFile,BaseFile,seconds);
sprintf(ThreeFile,"~/imageProcessing/trunk/3Dimages/%s_%d_seconds.root",BaseFile,seconds);
TFile* targets = new TFile(TargetFile);
TFile* image3D = new TFile(ThreeFile,"RECREATE");
TObjArray *ThreeList = new TObjArray();
char *threeN = "Targets_3D";
TH3F *three = new TH3F(threeN,threeN,noBinsX,0.,1., noBinsY,0.,1.,noBinsZ,0.,1.);
char HistoTarget[20];
for(int SliceNo = 1; SliceNo < noBinsZ; SliceNo++){
 sprintf(HistoTarget, "%s_%d", "Section", SliceNo);
TH2D* CurrentSlice = (TH2D*) targets->Get(HistoTarget);
for (int x = 1; x < noBinsX; x++) {
	 for (int y = 1; y < noBinsY; y++) {
	dummy=CurrentSlice->GetBinContent(x,y);
if((dummy>0)||(dummy==-1)){three->SetBinContent(x,y,SliceNo,1);}
}
}
}
char Name[20];
for (int imagenum = 1;imagenum<=targnum;imagenum++){
sprintf(Name,"Target_%d",imagenum);
TH3F *curtarg = new TH3F(Name,Name,noBinsX,0.,1., noBinsY,0.,1.,noBinsZ,0.,1.);
char HistoTarget[20];
for(int SliceNo = 1; SliceNo < noBinsZ; SliceNo++){
 sprintf(HistoTarget, "%s_%d", "Section", SliceNo);
TH2D* CurrentSlice = (TH2D*) targets->Get(HistoTarget);
for (int x = 1; x < noBinsX; x++) {
	 for (int y = 1; y < noBinsY; y++) {
	dummy=CurrentSlice->GetBinContent(x,y);
if(dummy==imagenum){curtarg->SetBinContent(x,y,SliceNo,1);}
}
}
}
ThreeList->Add(curtarg);
curtarg->GetXaxis()->SetTitle("X ID");
curtarg->GetYaxis()->SetTitle("Y ID");
curtarg->GetZaxis()->SetTitle("Z ID");
}
three->Write();
ThreeList->Write();
image3D->Close();
targets->Close();
}

int main (int argc, char**argv){
char inputFile[180];
char BaseFile[180];
char OutFile[30];
int seconds = 60;
int start, finish;
//char *input="~/imageProcessing/trunk/iterations/Lambda_40iterations_10cmtarget";
if(argc<3) {
      std::cerr << "Usage:\t" << "<input file location> <base name> <output location> <seconds(optional)>\n";
      return -1;
   }
strncpy(inputFile,argv[1],180);
strncpy(BaseFile,argv[2],180);
strncpy(OutFile,argv[3],30);
if(argc==5){
start=atoi(argv[4]);
finish=atoi(argv[4]);
}
if(argc>5){
start=atoi(argv[4]);
finish=atoi(argv[5]);
}
for(seconds=start;seconds<=finish;seconds++){
if(seconds==1||seconds%5==0){
char FileName[180];
int noBinsX = 100;
int noBinsY = 100;
int noBinsZ = 100;
sprintf(FileName,"%s%s_%d_seconds.root",inputFile,BaseFile,seconds);
TFile* in = new TFile(FileName);
std::cout << FileName << endl;
double y;
bool bo;
//y = backgroundSubtract(in,seconds,noBinsX,noBinsY,noBinsZ);
//std::cout << "Percentage difference in deviations is " << y << "" << endl;
/*std::cout << "Would you like to map image? (1/0)" << endl;
std::cin >> bo;
if(bo==1){*/
double z;
std::cout << "Now mapping significant lambdas" << endl;
  z = pictureMap(BaseFile,OutFile,seconds,in,noBinsX,noBinsY,noBinsZ);
std::cout << "Image significance compared to image average is " << z << " sigma" << endl;
int x=1;
/*std::cout << "Map 3D image for how many targets? (0 is no mapping)" << endl;
std::cin >> x;
if(x>0){*/
makethree(x,BaseFile,OutFile,seconds,noBinsX,noBinsY,noBinsZ);
//}
//}
in->Close();
}
}
return 0;
}


