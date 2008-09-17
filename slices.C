{
TFile* in = new TFile("analysed.root");

//TCanvas* out[20];

TCanvas* out = new TCanvas("out","out");

out->Divide(5,4,0.001,0.001,0);

//for (int j = 0; j < 5; j++)
//{
  
  //TH1* n[20];
  
  for ( int i = 0; i <20; i++)
    {
      out->cd(20-i);
      char name[10];
      char s[3];
      sprintf(s,"%d\n",i+1);
      strcpy(name, "_pz" );
      strcat(name,s);

      //bkgnd->GetZaxis()->SetRange(40,70);
     
      TH1* n = bkgnd->ProjectionZ(name,25,75, i*5, (i+1)*5 );

      n->Draw();
      //  }

    }
  
return 0;
}


/* char cName[4];
   char jNo[2];
  
   sprintf(jNo,"%d",j);
   strcpy(cName, "c" );
   strcat(cName,jNo);

   out[j] = new TCanvas(cName,cName, 10,10,1000,790);
  
   out[j]->Divide(4,4,0.001,0.001,0);*/
