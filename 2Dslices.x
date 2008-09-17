/* This should take a series of slices through the z direction showing the 2D histograms xy
 *
 * Sam 8-08-08
 *
 */

{
  TFile* in = new TFile("analysed.root","READ");

  TH3D* bkgnd = (TH3D*) in->Get ("bkgnd");

  TCanvas* out = new TCanvas("out","out");
  
  out->Divide(4,5,0.001,0.001,0);

  for (int i = 0; i < 20; i++)
    {
      out->cd(20-i);

      bkgnd->GetXaxis()->SetRange(i*5,(i+1)*5);
      
      char name[10];
      char s[3];
      sprintf(s,"%d\n",i+1);
      strcpy(name, "yz" );
      strcat(name,s);

      TH2D* yz = bkgnd->Project3D(name);

      yz->Draw();

      bkgnd->GetYaxis()->SetRange(0,100);
    }

  return 0;
}


