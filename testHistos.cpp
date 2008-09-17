/*
 * This should create a simple 3D histo gram that can be used to 
 * check edge detection systems
 *
 * Initially it will have a single cube in the centre with a high
 * density of points and nothing surronding it. 
 * 
 * Sam 06-08-08
 */

# include <TFile.h>
# include <TH3.h>
# include <TROOT.h>
# include <stdlib.h>
# include <iostream>


int main ()
{
  TFile* testFile = new TFile("test.root", "RECREATE");
  
  TH3D* test = new TH3D("test","test", 50,-50,50, 50,-50,50, 50,-50,50);
  
  for (int i = -10; i <= 10; i++)
    {
      for (int j = -2; j <= 2; j++)
	{
	  for (int k = -2; k <= 2; k++)
	    {
	      test->Fill (i, j, k, 5);
	    }
	}
    }
  
  testFile->Write();

  return 0;
}
