void compilething()
{
  char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_30cmtarget/pca/pca_fakecontainer_30cmtarget_million";
  //char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_10cmtarget/pca/pca_fakecontainer_10cmtarget_million";
  // char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_5cmtarget/pca/pca_fakecontainer_5cmtarget_million";
  //char *targetBase="/unix/anita1/creamtea/minerva/fakecontainer_notarget/pca/pca_fakecontainer_notarget_million";


  Int_t numTargetMillions=10;
  Int_t numTarget=numTargetMillions*1000;
  char fileName[180];

  TChain *targetTree = new TChain("pcaTree");
  for(Int_t startFile=1;startFile<numTargetMillions;startFile++) {
    sprintf(fileName,"%s_%d.root",targetBase,startFile);
    targetTree->Add(fileName);
    //    cout << fileName << endl;
  }

  int muons = 112000; //28K = 1 second.
  int Nx, Ny, Nz;
  Nx = Ny = Nz = 100;


  gSystem->CompileMacro("PcaTreeLooper.C","k");






  PcaTreeLooper targetLooper(targetTree);
  targetLooper.SLFill(0,muons,Nx,Ny,Nz);
  targetLooper.LambdaFill(Nx*Ny*Nz);
  targetLooper.SigmaFill(0,muons);
  targetLooper.GradientFill(0,muons);

  
  for(int iii = 0; iii < 3; iii++){

    double ax = -1.0e-5;
    double bx = 0.0;
    double cx = 1.0e-5;
    double tol = 1.0e-13;
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
    fw=fv=fx=targetLooper.Cost(x);
    
    for(iter=1;iter<=ITMAX;iter++){ //main program loop.
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      cout << iii << "\t" << x << ":\t" <<  fx << endl;
      if(fabs(x-xm) <= (tol2-0.5*(b-a))){ //Test for done here.
	xmin=x;
	cout << iii <<  "\tReturn: " << x << ":\t" <<  fx << endl;
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
      
      fu=targetLooper.Cost(u); //This is the one function evaluation per iteration.
      
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
      cout << iii <<  "\tReturn: " << x << ":\t" <<  fx << endl;
    }
  cout << "**********************" << iii+1 << "***************************" << endl;

  targetLooper.LambdaAlpha(xmin);
  targetLooper.SigmaFill(0,muons);
  targetLooper.GradientFill(0,muons);
  }

  targetLooper.DrawGradientSlices(50,48,100,100,100);
  

}

