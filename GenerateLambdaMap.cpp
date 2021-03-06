#include <TFile.h>
#include <TH3.h>
#include <TTree.h>
#include <TROOT.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <map>
#include <cstring>
#include "/unix/creamtea/sfayer/10by10by10_10cm/include/DetectorDefs.hh"
#include "LambdaPcaTreeLooper.C"

using namespace std;

int main(int argc, char** argv){
char iterationsFile[180];
char targetBase[180];
  int seconds = 60;
if(argc<3) {
      std::cerr << "Usage:\t" << "<input file location and base name> <output name> <seconds(optional, default=60)>\n";
return -1;
   }
strncpy(targetBase,argv[1],180);
strncpy(iterationsFile,argv[2],180);
  if(argc>3) seconds=atoi(argv[3]);
  int muons = seconds*10000/60; //28K = ~1 second.
  cout << "Number of muons: " << muons << endl;

  //char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_30cmtarget/pca/pca_fakecontainer_30cmtarget_million";
  //char *targetBase="/unix/creamtea/sfayer/1mdetector_10cmtarget_iron/pca/pca_small1mdetector_10cmtarget_iron_million";
  //char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_5cmtarget/pca/pca_fakecontainer_5cmtarget_million";
  //char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_notarget/pca/pca_fakecontainer_notarget_million";

  Int_t numTargetMillions=1;
  Int_t numTarget=numTargetMillions*1000;
  char fileName[180];

  TChain *targetTree = new TChain("pcaTree");
  for(Int_t startFile=1;startFile<=numTargetMillions;startFile++) {
    sprintf(fileName,"%s_%d.root", targetBase,startFile);
    targetTree->Add(fileName);
    cout << fileName << endl;
  }
  int Nx, Ny, Nz;
  Nx = 400;
  Ny = 400;
  Nz = 400;
  char fileNameLambda[180];

  int iterations = 40;

  sprintf(fileNameLambda,"~/imageProcessing/trunk/iterations/pic400Lambda_%diterations_%s_%d_seconds.root",iterations, iterationsFile,seconds);
  //gSystem->CompileMacro("LambdaPcaTreeLooper.C","k");
  cout << fileNameLambda << endl;

  LambdaPcaTreeLooper targetLooper(targetTree);
  targetLooper.SLFill(0,muons,Nx,Ny,Nz);
  targetLooper.LambdaFill(Nx*Ny*Nz);
  targetLooper.SigmaFill();
  targetLooper.GradientFill();
  
  for(int iii = 0; iii < iterations; iii++){

    double ax = -1.0e-5;
    double bx = 0.0;
    double cx = 1.0e-5;
    double tol = 3.0e-8;
    double xmin = 0;
    int test = 0;
    
    int ITMAX = 100;
    double CGOLD = 0.3819660;
    double ZEPS = 1.0e-10;
    int iter;
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
    double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    double e=0.0; //this will be the distance moved the step before last.
    a=(ax < cx ? ax : cx); //a and b must be in ascending order
    b=(ax > cx ? ax : cx); //but input abscissas need not be.
    x=w=v=bx; //initialisations...
    fw=fv=fx=targetLooper.Cost(x,0,muons);
    
    for(iter=1;iter<=ITMAX;iter++){ //main program loop.
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      // cout << iii << "." << iter << "\t" << x << ":\t" <<  fx << endl;
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
      
      fu=targetLooper.Cost(u,0,muons); //This is the one function evaluation per iteration.
      
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
    }
    if(test == 0){
      xmin=x; //Never get here
      cout << iii << "." << iter <<  "\tReturn: " << x << ":\t" <<  fx << endl;
    }
    cout << "**********************" << iii+1 << "***************************" << endl;

    targetLooper.LambdaAlpha(xmin);
    targetLooper.SigmaFill();
    targetLooper.GradientFill();
  }

  targetLooper.DrawSlices(50,48,Nx,Ny,Nz,fileNameLambda);
}
