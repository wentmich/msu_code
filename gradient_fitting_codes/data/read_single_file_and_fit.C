
//Global variables
double Leff, r0, Iamps;
char fname[]="EngeProd_0.10.txt";

double func_EngeProd(double *z, double *p)
{
  double zi =  z[0] - Leff/2.0;
  double ze = -z[0] - Leff/2.0;
  zi /= 2*r0; ze /= 2*r0;
  double fi =  1./(1.+exp(p[0]+zi*(p[1]+zi*(p[2]+zi*(p[3]+zi*(p[4]+zi*p[5]))))));
  double fe =  1./(1.+exp(p[6]+ze*(p[7]+ze*(p[8]+ze*(p[9]+ze*(p[10]+ze*p[11]))))));
  //printf(" z= %lg    f= %lg\n", z, p[1]);
  //printf("   p = %lg %lg %lg %lg\n", p[0], p[1], p[2], p[3]);
  return p[12]*fi*fe;
}

void read_single_file_and_fit()
{
  
	
	
	
	
	
	
	
	
	
	
	
	
  //Examples using class EngeProdObj and TF1 style funcEngeProd
  // see section 5.3.3 of ROOTUsersGuideLetter2014.pdf 
  const int Np = 13;  //number of parameter
  
  //Initiallize values
  double p[Np] =
  { 0., 3., 0., 0., 0., 0.    //0-5, Enge incident side
   ,0., 4., 0., 0., 0., 0.    //6-11, exit
   ,0.3                                //12, B0
  };
  //Initialize
  r0 = 0.2;         //reference radius
  double ylim = 2.2;
  const int NLINE=256;
  FILE *fp;         //file pointer
  char line[NLINE]; //buffer to read lines
  
  double zmax=1., zmin=-zmax, delz = 0.005;   int Nz;
  double *za, *fa;
  double zz[1], ftmp;

  Nz = 1 + (int)( (zmax-zmin)/delz );
  printf(" points per excitation, Nz= %i \n", Nz);
  za = (double*) malloc( Nz * sizeof(double) );
  fa = (double*) malloc( Nz * sizeof(double) );

  //Read field gradient data from file. Generate scatter plot.
      fp = fopen(fname, "r");
      fgets(line, NLINE, fp);   sscanf(line, " %lg ", &Iamps);
      fgets(line, NLINE, fp);   sscanf(line, " %lg ", &r0);
      fgets(line, NLINE, fp);   sscanf(line, " %i ", &Nz);
      fgets(line, NLINE, fp);
      printf(" Iamps= %lg, r0= %lg, Nz= %i \n", Iamps, r0, Nz);
      for( int k=0; k<Nz; k++)
      {
        fgets(line, NLINE, fp);
        sscanf(line," %lg %lg\n", &za[k], &fa[k]);
        //printf(" %lg %lg\n", za[k], fa[k]);
      }
      fclose(fp);
      
      //-------- create graph and do Fit and other here
      TGraph * data = new TGraph(Nz, za, fa);
      
      //Plot data in TCanvas
      if( (TCanvas*)gROOT->FindObject("cEP") ) {} else
      TCanvas * cEP = new TCanvas("cEP","cEP", 40,40, 700,450);
         data->SetTitle("field; z; f");
         data->SetLineColor(0);
         //data->SetMaximum(ylim);
         data->SetMarkerStyle(1);
         data->Draw();
   
  // find effective length using integral and b(z=0)
  int z_zero_index;
  for (int i = 0; i < Nz; i++)
  {
  	   if (fabs(za[i]) < 1.0e-10)
  	   	   z_zero_index = i;
  }
  double Bzero = fa[z_zero_index];
  double integral = data -> Integral();
  Leff = integral/Bzero;
  printf("Leff = %lg \n",Leff);
   
  //Define the fuction to fit
  TF1 *func = new TF1( "func", func_EngeProd, zmin,zmax, Np);
  for( int i=0; i<Np; i++) func->SetParameter(i, p[i]);
     
  //Optional: Fit the parameters to data data 
  bool fit_TF1=1;  //Set to 1 to allow fitting
  if(fit_TF1)
  {
    data->Fit("func");
  }//fit_TF1
  func->Draw("same");
  
  //See also https://root.cern.ch/root/HowtoFit.html
  //printf(" p[0] = %lg   \n", p[0]);
  
  double coefficients[Np];
  for (int i = 0; i < Np; ++i)
  	  coefficients[i] = func -> GetParameter(i);

  
  return;
}