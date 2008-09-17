{
TFile* in = new TFile ("analysed.root","READ");
TCanvas* out = new TCanvas("out","out");
out->Divide (3,3,0.001,0.001);

for (int i = 0; i<9; i++)
{
  out->cd (i+1);
  
  char name [20];
  char number [5];
  
  strcpy (name,"absorbed");
  sprintf (number,"%d",(i*100+500) );
  strcat (name,number);
  
  TH2D* oa = (TH2D*)in->Get(name);
  oa->DrawClone("same");

  int top = (i*5) + 25;
  int bottom = ((i+1) * 5) + 25;

  bkgndP->GetYaxis ()->SetRange (top,bottom); /* blue */
  TH2D* op = bkgndP->Project3D ("xz");
  op->SetMarkerColor (4);
  op->DrawClone("same");

  bkgndN->GetYaxis()->SetRange (top,bottom); /* red */
  TH2D* on = bkgndN->Project3D ("xz");
  on->SetMarkerColor (2);
  on->DrawClone("same");
}



}
