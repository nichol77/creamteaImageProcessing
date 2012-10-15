#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TGraph.h>
#include <TTree.h>
#include <TROOT.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "/unix/creamtea/sfayer/1mdetector_5cmtarget/include/DetectorDefs.hh"
#include <cmath>
#include <iostream>

using namespace std;

void imageQuality () {
TObjArray *HlistSignal = new TObjArray(35);
TObjArray *HlistBackground = new TObjArray(35);
TGraph *FoM_graph = new TGraph();
char HistName[80];
char FileName[180];
char HistoSign[80];
char HistoBack[80];
int noBinsX = NX;
 int noBinsY = NY;
int noBinsZ = NZ;
double side_length = SIDELENGTH;
double dx = side_length/noBinsX;
double dy = side_length/noBinsY;
double dz = side_length/noBinsZ;
double sphere_radius = SPHERE_RADIUS_CM/100.; // everything in m
double sphere_x = SPHERE_X_M ;
double sphere_y = SPHERE_Y_M ;
double sphere_z = SPHERE_Z_M ;
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
int point = 1 ;
char *input="~/imageProcessing/trunk/iterations/Lambda_40iterations_5cmTarget";
for (int NFile = 1; NFile < 201 ; NFile ++) {
  if (NFile == 1 ||( (NFile<131)&&(NFile%5 == 0))||((NFile>=131)&&(NFile%10 == 0)) ) {
sprintf(FileName,"%s_%d_seconds.root",input,NFile);
TFile* in = new TFile(FileName);
 sprintf(HistoSign, "%s_%d", "Signal", NFile);
sprintf(HistoBack, "%s_%d", "Background", NFile);
//Histo for Current Signal and Current Background
TH1F *CS = new TH1F(HistoSign, HistoSign,100,0.,1.);
TH1F *CB = new TH1F(HistoBack, HistoBack,100,0.,1.);
HlistSignal->Add(CS);
HlistBackground->Add(CB);
for(int SliceNo = in_z; SliceNo <= fin_z; SliceNo++){
sprintf(HistName,"Slice_%d",SliceNo);
 TH2D* CurrentSlice = (TH2D*) in->Get(HistName);
z_pos = dz*(SliceNo*1. + 0.5);
for(int x = in_x; x <= fin_x; x++){
x_pos = dx*(x*1. + 0.5);
for(int y = in_y; y <= fin_y; y++){
y_pos = dy*(y*1. + 0.5);
distance = sqrt(pow((x_pos - sphere_x),2.) + pow((y_pos - sphere_y),2.) + pow((z_pos - sphere_z),2.) ) ;
if (distance <= sphere_radius){
CS->Fill(CurrentSlice->GetBinContent(x,y));
} else CB->Fill(CurrentSlice->GetBinContent(x,y));
} // end for y
} // end for x
} // end slices z
double FoM = (CS->GetMean() - CB->GetMean())/ CB->GetRMS();
double termSign = pow(CS->GetMeanError()/CB->GetRMS(), 2.);
double termBack = pow(CB->GetMeanError()/CB->GetRMS(), 2.);
double termSigmaBack = pow(CB->GetRMSError()*FoM/CB->GetRMS(), 2.);
double FoMerror = sqrt(termSign + termBack + termSigmaBack);
FoM_graph->SetPoint(point, NFile*1., FoM);
//FoM_graph->SetPointError(point, 0., FoMerror);
point ++;
} }
TFile* out = new TFile("~/imageProcessing/trunk/FoMout.root", "RECREATE");
 HlistSignal->Write();
HlistBackground->Write();
FoM_graph->SetTitle("Figure_of_Merit");
FoM_graph->SetName("FoM_graph");
FoM_graph->GetXaxis()->SetTitle("Seconds");
FoM_graph->GetYaxis()->SetTitle("FoM");
FoM_graph->Write();
out->Close();
}

int main(){
imageQuality();
return 0;
}
