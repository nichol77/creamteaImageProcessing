#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TGraph.h>
#include <TTree.h>
#include <TROOT.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

void imageQuality () {
TFile* out = new TFile("~/imageProcessing/trunk/FoMTargetout.root", "RECREATE");
TObjArray *HlistSignal = new TObjArray(35);
TObjArray *HlistBackground = new TObjArray(35);
TObjArray *ListFoM = new TObjArray(20);
TObjArray *ListTot = new TObjArray(20);
TObjArray *ListPer = new TObjArray(20);
TObjArray *ListGhost = new TObjArray(20);
char HistName[80];
 char HistTarget[80];
char FileName[180];
 char FileTarget[180];
char HistoSign[80];
char HistoBack[80];
char FoMTitle[80];
char TotTarTitle[80];
char GhostTitle[80];
char PerTitle[80];
double sphere_x,sphere_y,sphere_z;
for(int TarSize =0;TarSize<=10;TarSize++) {
TGraph *FoM_graph = new TGraph();
TGraph *percentages = new TGraph();
TGraph *ghosting = new TGraph();
TGraph *total_targets = new TGraph();
  std::cout<<"Now calculationing target size: " << TarSize << endl;
if(TarSize==1){//begin position allocation
sphere_x=0.2;
sphere_y=-0.05;
sphere_z=0;
} else {
if(TarSize==2){
sphere_x=0.05;
sphere_y=0.2;
sphere_z=0;
} else {
if(TarSize==3){
sphere_x=0.1;
sphere_y=0;
sphere_z=0;
} else {
if(TarSize==4){
sphere_x=-0.1;
sphere_y=0;
sphere_z=0;
} else {
if(TarSize==5){
sphere_x=0;
sphere_y=0;
sphere_z=0;
} else {
if(TarSize==6){
sphere_x=0;
sphere_y=0.1;
sphere_z=0;
} else {
if(TarSize==7){
sphere_x=0;
sphere_y=-0.1;
sphere_z=0;
} else {
if(TarSize==8){
sphere_x=0.1;
sphere_y=0.1;
sphere_z=0;
} else {
if(TarSize==9){
sphere_x=-0.1;
sphere_y=-0.1;
sphere_z=0;
} else {
if(TarSize==10){
sphere_x=0.25;
sphere_y=0.3;
sphere_z=0;
}else{
sphere_x=0;
sphere_y=0;
sphere_z=0;
} } } } } } } } } } //end position allocation
char *TarType = "";//including underscore before
int noBinsX = 100;
 int noBinsY = 100;
int noBinsZ = 100;
double side_length = 1;//1m target volume
double dx = side_length/noBinsX;
double dy = side_length/noBinsY;
double dz = side_length/noBinsZ;
double sphere_radius = TarSize/100.; // everything in m
// half of container dimensions in m
double container_x = (0.99/2);
double container_y = (0.99/2);
double container_z = (0.99/2);
// Number of bins inside the container
int in_x = (noBinsX/2) - int(container_x/dx) ;
int fin_x = (noBinsX/2) + int(container_x/dx) ;
int in_y = (noBinsY/2) - int(container_y/dy) ;
int fin_y = (noBinsY/2) + int(container_y/dy) ;
int in_z = (noBinsZ/2) - int(container_z/dz) ;
int fin_z = (noBinsZ/2) + int(container_z/dz) ;
//Normalise the position of the sphere to the slice coordinate;
sphere_x = (sphere_x + (side_length/2.));
sphere_y = (sphere_y + (side_length/2.));
sphere_z = (sphere_z + (side_length/2.));
sphere_radius += sqrt((dx*dx)+(dy*dy)+(dz*dz))/2.;
double x_pos = dx/2.;
double y_pos = dy/2.;
double z_pos = dz/2.;
double distance = 0.;
 double target;
 int point, intarget, notreg, overreg, outtarget, total;
 if(TarSize==0) { point = 1;}
for (int NFile = 1; NFile <= 90 ; NFile ++) {
  if (NFile == 1 ||(NFile%5 == 0)) {
if(TarSize==0){
sprintf(FileName,"/unix/creamtea/sfayer/blank_target/lambdamap/Lambda_40iterations_blank_target_%d_seconds.root",NFile);
sprintf(FileTarget,"/unix/creamtea/sfayer/blank_target/targets/Lambda_40iterations_blank_target_%d_seconds.root",NFile);
}else{
  sprintf(FileName,"/unix/creamtea/sfayer/1mdetector_%dcmtarget/lambdamap/Lambda_40iterations_%dcmtarget%s_%d_seconds.root",TarSize,TarSize,TarType,NFile);
  sprintf(FileTarget,"/unix/creamtea/sfayer/1mdetector_%dcmtarget/targets/Lambda_40iterations_%dcmtarget%s_%d_seconds.root",TarSize,TarSize,TarType,NFile);
}
TFile* in = new TFile(FileName);
 TFile* Ref = new TFile(FileTarget);
out->cd();
 sprintf(HistoSign, "%s_%d", "Signal", NFile);
sprintf(HistoBack, "%s_%d", "Background", NFile);
//Histo for Current Signal and Current Background
TH1F *CS = new TH1F(HistoSign, HistoSign,100,0.,1.);
TH1F *CB = new TH1F(HistoBack, HistoBack,100,0.,1.);
HlistSignal->Add(CS);
HlistBackground->Add(CB);
intarget=0;
outtarget=0;
notreg=0;
overreg=0;
for(int SliceNo = in_z; SliceNo <= fin_z; SliceNo++){
sprintf(HistName,"Slice_%d",SliceNo);
sprintf(HistTarget,"Section_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) in->Get(HistName);
 TH2F* Targetref = (TH2F*) Ref->Get(HistTarget);
z_pos = dz*(SliceNo*1. + 0.5);
for(int x = in_x; x <= fin_x; x++){
x_pos = dx*(x*1. + 0.5);
for(int y = in_y; y <= fin_y; y++){
y_pos = dy*(y*1. + 0.5);
distance = sqrt(pow((x_pos - sphere_x),2.) + pow((y_pos - sphere_y),2.) + pow((z_pos - sphere_z),2.) ) ;
 double target = Targetref->GetBinContent(x,y);
 if (distance <= sphere_radius) {
   if ((target)>=1) {
     CS->Fill(CurrentSlice->GetBinContent(x,y));
     intarget++;
     } else {
     CB->Fill(CurrentSlice->GetBinContent(x,y));
notreg++;
}
} else {
   CB->Fill(CurrentSlice->GetBinContent(x,y));
   if ((target)>=1){
     overreg++;
       } else outtarget++;
      }
} // end for y
} // end for x
} // end slices z
 total = intarget+overreg+outtarget+notreg;
 double reg = ((intarget)*100.);
reg/=(overreg+intarget);
if (isnan(reg)) reg = 0;
 double tarpercent = (intarget*100.);
tarpercent/=(notreg+intarget);
 double correct = ((intarget+outtarget)*100.);
correct/=total;
double FoM = (CS->GetMean() - CB->GetMean())/ CB->GetRMS();
double termSign = pow(CS->GetMeanError()/CB->GetRMS(), 2.);
double termBack = pow(CB->GetMeanError()/CB->GetRMS(), 2.);
double termSigmaBack = pow(CB->GetRMSError()*FoM/CB->GetRMS(), 2.);
double FoMerror = sqrt(termSign + termBack + termSigmaBack);
FoM_graph->SetPoint(point, NFile*1., FoM);
FoM_graph->SetMarkerColor(TarSize);
total_targets->SetPoint(point, NFile*1.,(intarget+overreg));
 percentages->SetPoint(point, NFile*1., tarpercent);
ghosting->SetPoint(point, NFile*1., reg);
//FoM_graph->SetPointError(point, 0., FoMerror);
//std::cout << intarget << " " << outtarget << " " << overreg << " " << notreg << " " << total << " " << tarpercent << " " << reg << endl;
point ++;
in->Close();
Ref->Close();
} }
sprintf(FoMTitle,"FoM_%d",TarSize);
FoM_graph->SetTitle(FoMTitle);
FoM_graph->SetName(FoMTitle);
FoM_graph->GetXaxis()->SetTitle("Seconds");
FoM_graph->GetYaxis()->SetTitle("FoM");
FoM_graph->SetMarkerColor(kRed+TarSize);
FoM_graph->Draw("alp");
ListFoM->Add(FoM_graph);
sprintf(TotTarTitle,"Total_%d",TarSize);
total_targets->SetTitle(TotTarTitle);
 total_targets->SetName(TotTarTitle);
total_targets->GetXaxis()->SetTitle("Seconds");
total_targets->GetYaxis()->SetTitle("Voxels of Target found");
total_targets->SetMarkerColor(kRed+TarSize);
total_targets->Draw("same");
ListTot->Add(total_targets);
sprintf(PerTitle,"Percentages_%d",TarSize);
 percentages->SetTitle(PerTitle);
 percentages->SetName(PerTitle);
percentages->GetXaxis()->SetTitle("Seconds");
percentages->GetYaxis()->SetTitle("Percentage");
 ListPer->Add(percentages);
sprintf(GhostTitle,"Ghosting_%d",TarSize);
ghosting->SetTitle(GhostTitle);
 ghosting->SetName(GhostTitle);
ghosting->GetXaxis()->SetTitle("Seconds");
ghosting->GetYaxis()->SetTitle("Percentage");
 ListGhost->Add(ghosting);
}
out->cd();
HlistSignal->Write();
HlistBackground->Write();
ListTot->Write();
ListFoM->Write();
ListPer->Write();
ListGhost->Write();
out->Close();
}

int main(){
imageQuality();
return 0;
}
