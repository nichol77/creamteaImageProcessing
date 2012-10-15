#include <TFile.h>
#include <TH3.h>
#include <TTree.h>
#include <TROOT.h>
#include <stdlib.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <cstring>
#include "/unix/creamtea/sfayer/1mdetector_10cmtarget/include/DetectorDefs.hh"
#include "LambdaPcaTreeLooperMLSD.C"

using namespace std;

/*
Parabolic Interpolation and Brent's method in 1D
from "Numerical recipes in C++", pag 402
*/

int main(int argc, char** argv){

  int seconds = 1;  // Or you can give a different number of seconds from terminal
  if(argc > 1){
    seconds = atoi(argv[1]);
  } 

  int muons = seconds*10000/60; //28K = ~1 second for 13 m detector, 1 (muon/cm)/minute. max 315381

  cout << "Number of muons: " <<muons << endl;

  //char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_30cmtarget/pca/pca_fakecontainer_30cmtarget_million";
//  char *targetBase="/unix/creamtea/minerva/fakecontainer_10cmtarget/pcaMLSD/pca_fakecontainer_10cmtarget_million";
  char *targetBase="/unix/creamtea/sfayer/1mdetector_10cmtarget/pca/pca_small1mdetector_10cmtarget_million";
  //char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_5cmtarget/pca/pca_fakecontainer_5cmtarget_million";
  //char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_notarget/pca/pca_fakecontainer_notarget_million";

  Int_t numTargetMillions=2;
  Int_t numTarget=numTargetMillions*1000;
  char fileName[180];

  TChain *targetTree = new TChain("pcaTree");
  for(Int_t startFile=1;startFile<=numTargetMillions;startFile++) {
    sprintf(fileName,"%s_%d.root",targetBase,startFile);
    targetTree->Add(fileName);
    cout << fileName << endl;
  }

  int Nx, Ny, Nz;
  Nx = 10;
    Ny = 10;
    Nz = 10;
  char fileNameLambda[180];

  int iterations = 50;

  sprintf(fileNameLambda,"~/imageProcessing/trunk/iterations/Lambda_%diterations_10cmtarget_%d_secondsXYTHRESH2GeV.root",iterations,seconds);
char *fileNameMuons="~/imageProcessing/trunk/tmp/inputMuons.root";

  //gSystem->CompileMacro("LambdaPcaTreeLooperMLSD.C","k");

  LambdaPcaTreeLooperMLSD targetLooperMLSD(targetTree);
  targetLooperMLSD.SLFill(0,muons,Nx,Ny,Nz);
  targetLooperMLSD.LambdaFill(Nx*Ny*Nz);
  targetLooperMLSD.SigmaFill();
//  targetLooperMLSD.GradientFill();
double fx;
  for(int iii = 0; iii < iterations; iii++){
/*
    double ax = -1.0e-6;
    double bx = 0.0;
    double cx = 1.0e-6;
    double tol = 3.0e-8; // this tol has been chosen for double. it is 2*sqrt(epsilon_double)
    double xmin = 0;
    int test = 0;
 
    int ITMAX = 100; //ITMAX is the maximum allowed number of iterations
    double CGOLD = 0.3819660; // Golden Ratio
    double ZEPS = 1.0e-10; //ZEPS small number that protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero
// ***********************************************************
//Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
//between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
//the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
//the minimum is returned as xmin, and the minimum function value is returned as brent, the
//returned function value.
// ************************************************************
  int iter;
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) // if b>=0. return fabs(a) otherwise return -fabs(a)
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
    double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    double e=0.0; //this will be the distance moved the step before last.

    a=(ax < cx ? ax : cx); //a and b must be in ascending order
    b=(ax > cx ? ax : cx); //but input abscissas need not be.
    x=w=v=bx; //initialisations...
    fw=fv=fx=targetLooperMLSD.Cost(x,0,muons);
   
    for(iter=1;iter<=ITMAX;iter++){ //main program loop.
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      cout << iii << "." << iter << "\t" << x << ":\t" <<  fx << endl;
      if(fabs(x-xm) <= (tol2-0.5*(b-a))){ //Test for done here.
	xmin=x;
	cout << iii << "." << iter <<  "\tReturn: " << x << ":\t" <<  fx << endl;
	test = 1;
	break;
      }
      if(fabs(e) > tol1){ //Construct a trial parabolic fit.
	r=(x-w)*(fx-fv);
	q=(x-v)*(fx-fw);
	p=(x-v)*q-(x-w)*r;
	q=2.0*(q-r);
	if(q > 0.0) p = -p;
	q=fabs(q);
	etemp=e;
	e=d;
	if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	  d=CGOLD*(e=(x >= xm ? a-x : b-x));
	//The above conditions determine the acceptability of the parabolic fit. Here we take the golden section step into the larger of the two segments.
	else {
	  d=p/q; //Take the parabolic step.
	  u=x+d;
	  if (u-a < tol2 || b-u < tol2)
	    d=SIGN(tol1,xm-x);
	}
      }else{
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      }
      u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));

    targetLooperMLSD.LambdaAlpha(u);      
  targetLooperMLSD.SigmaFill();
      fu=targetLooperMLSD.Cost(u,0,muons); //This is the one function evaluation per iteration.
    targetLooperMLSD.LambdaAlpha(-u);
     
      if(fu <= fx){ //now decide what to do with our evaluation.
	if(u >= x) a=x; else b=x;
	SHFT(v,w,x,u) //Housekeeping follows:
	  SHFT(fv,fw,fx,fu)
	  }else{
	if(u < x) a=u; else b=u;
	if(fu <= fw || w == x){
	  v=w;
	  w=u;
	  fv=fw;
	  fw=fu;
	}else if (fu <= fv || v == x || v == w){
	  v=u;
	  fv=fu;
	}
      } //done with housekeeping, back for another iteration
    } //end loop iter
    if(test == 0){
      xmin=x; //Never get here
      cout << iii << "." << iter <<  "\tReturn: " << x << ":\t" <<  fx << endl;
    }*/
    cout << "**********************" << iii+1 << "***************************" << endl;


    fx=targetLooperMLSD.Cost(0.1,0,muons,iii);
    targetLooperMLSD.LambdaNew();
//    targetLooperMLSD.LambdaAlpha(xmin);
    targetLooperMLSD.SigmaFill();
//    targetLooperMLSD.GradientFill();
  } //end loop iii

  targetLooperMLSD.DrawSlices(Nz,1,Nx,Ny,Nz,fileNameLambda);
  targetLooperMLSD.DrawMuons(50,Nx,Ny,Nz,fileNameMuons);
 
}
