// Generate simulated B field vs z and current for some arbitrary magnet
// Write B_vs_z for varying I/Imax to series of files
//  2018-06
class EngeProdObj {
  public:
    EngeProdObj();
    ~EngeProdObj();
    double Evaluate(double x) const;
    double Derivative(double x) const;
    double Derivative2(double x) const;
    double Derivative3(double x) const;
  private:
    double L;
};

//Global variables
double Leff, r0;

//.C, =======================================
EngeProdObj::EngeProdObj(){ }
EngeProdObj::~EngeProdObj(){ }

double EngeProdObj::EngeProd(double z, double *p)
{
  double zi = z - Leff/2.0;
  double ze = -z - Leff/2.0;
  zi /= 2*r0; ze /= 2*r0;
  double fi =  1./(1.+exp(p[0]+zi*(p[1]+zi*(p[2]+zi*(p[3]+zi*(p[4]+zi*p[5]))))));
  double fe =  1./(1.+exp(p[6]+ze*(p[7]+ze*(p[8]+ze*(p[9]+ze*(p[10]+ze*p[11]))))));
  //printf(" z= %lg    f= %lg\n", z, p[1]);
  //printf("   p = %lg %lg %lg %lg\n", p[0], p[1], p[2], p[3]);
  return p[51]*fi*fe;
}
double func_EngeProd(double *z, double *p)
{
  double zi = z[0] - Leff/2.0;
  double ze = -z[0] - Leff/2.0;
  zi /= 2*r0; ze /= 2*r0;
  double fi =  1./(1.+exp(p[0]+zi*(p[1]+zi*(p[2]+zi*(p[3]+zi*(p[4]+zi*p[5]))))));
  double fe =  1./(1.+exp(p[6]+ze*(p[7]+ze*(p[8]+ze*(p[9]+ze*(p[10]+ze*p[11]))))));
  //printf(" z= %lg    f= %lg\n", z, p[1]);
  //printf("   p = %lg %lg %lg %lg\n", p[0], p[1], p[2], p[3]);
  return p[51]*fi*fe;
}
double f_c0i(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[12]+p[13]*aI+p[14]*aI*aI; }
double f_c0e(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[15]+p[16]*aI+p[17]*aI*aI; }
double f_c1i(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[18]+p[19]*aI+p[20]*aI*aI; }
double f_c1e(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[21]+p[22]*aI+p[23]*aI*aI; }
double f_c2i(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[24]+p[25]*aI+p[26]*aI*aI; }
double f_c2e(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[27]+p[28]*aI+p[29]*aI*aI; }
double f_c3i(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[30]+p[31]*aI+p[32]*aI*aI; }
double f_c3e(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[33]+p[34]*aI+p[35]*aI*aI; }
double f_c4i(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[36]+p[37]*aI+p[38]*aI*aI; }
double f_c4e(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[39]+p[40]*aI+p[41]*aI*aI; }
double f_c5i(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[42]+p[43]*aI+p[44]*aI*aI; }
double f_c5e(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[45]+p[46]*aI+p[47]*aI*aI; }
//Field parameter B0 vs I
double f_B0(double *I, double *p){ 
  double aI = fabs(I[0]);  return p[48]+p[49]*aI+p[50]*aI*aI; }

double EngeProdObj::Evaluate(double x, double *p){ return EngeProd( x, p); }

double EngeProdObj::Derivative(double x, double *p, double dx=0.001){
  return (EngeProd(x+dx,p)-EngeProd(x-dx,p))/(2.*dx);
}

double EngeProdObj::Derivative2(double x, double *p, double dx=0.001){
  return ( EngeProd(x+dx,p) -2*EngeProd(x,p) + EngeProd(x-dx,p) ) / (dx*dx);
}

double EngeProdObj::Derivative3(double x, double *p, double dx=0.001){
  return ( EngeProd(x+2*dx,p) -2*EngeProd(x+dx,p)
    +2*EngeProd(x-dx,p) -EngeProd(x-2*dx,p) ) / (2.*dx*dx*dx);
}


void generate_field_vs_z_and_curr()
{
  //Examples using class EngeProdObj and TF1 style funcEngeProd
  // see section 5.3.3 of ROOTUsersGuideLetter2014.pdf 
  const int Np = 52;
  double p[Np] =
  { 0., 3., 0., 0., 0., 0.    //0-5, Enge incident side
   ,0., 4., 0., 0., 0., 0.    //6-11, exit
   //Enge coef vs current based on a1900 dipole study, params_30d-mode.xlsx
   ,3.95e-1, 3.61e-3, -1.51e-5     //12-14, c0i dependence on I
   ,3.95e-1, 3.61e-3, -1.51e-5          //15-17, c0e dependence on I
   ,4., 4.18e-4, -1.45e-5   //18-20, c1i dependence on I
   ,4., 4.18e-4, -1.45e-5   //21-23, c1e dependence on I
   ,-1.34, 9.51e-3, -1.67e-5     //24-26, c2i dependence on I
   ,-1.34, 9.51e-3, -1.67e-5     //27-29, c2e dependence on I
   ,0.535,-4.01e-3, 1.02e-5      //30-32, c3i dependence on I
   ,0.535,-4.01e-3, 1.02e-5      //33-35, c3e dependence on I
   ,0., 0.0, 0.0          //36-38, c4i dependence on I
   ,0., 0.0, 0.0          //39-41, c4e dependence on I
   ,0., 0.0, 0.0          //42-44, c5i dependence on I
   ,0., 0.0, 0.0          //45-47, c5e dependence on I
   ,-0.114, 2.15e-2, -5.26e-5          //48-50, B0 dependence on I
   ,2.0                                //51, B0
   }; 
  Leff = 0.7; r0 = 0.2;
  double ylim = 2.2;
    double Ilim = 180;
  const int NLINE=256;
  FILE *fp, *fp2; char line[NLINE];
  
  double zmax=1., zmin=-zmax, delz = 0.005;   int Nz;
  double *za, *fa;
  double zz[1], ftmp;
  
  Nz = 1 + (int)( (zmax-zmin)/delz );
  printf(" points per excitation, Nz= %i \n", Nz);
  za = (double*) malloc( Nz * sizeof(double) );
  fa = (double*) malloc( Nz * sizeof(double) );
    
  double frac, frac0=0.1, dfrac = 0.1;
  double *Ia, *cf, Iamps[1];
  int Nfrac;
  Nfrac = 1 + (int)( (1.0-frac0)/dfrac );
  printf(" points of excitation, Nfrac= %i \n", Nfrac);
  Ia = (double*) malloc( Nfrac*sizeof(double) );
  cf = (double*) malloc( Nfrac*sizeof(double) );
    for( int n=0; n<Nfrac; n++){
      frac = n*dfrac+frac0;
      Ia[n] = frac*Ilim;
      printf(" n=%i  frac=%lg  Iamps=%lg \n", n, frac, Ia[n]);
    }

  //Implement Enge coef vs I dependence
  bool gen_EngeProd_Z_I=1;
  if(gen_EngeProd_Z_I)
  {    
    //Evaluate func_EngeProd over range of excitation
    TF1 *func = new TF1( "func", func_EngeProd, zmin,zmax, Np);
    for( int i=0; i<Np; i++) func->SetParameter(i, p[i]);
    
    //Write Enge coefs vs Iamps
    sprintf(line, "EngeCoefs.csv");
    fp2 = fopen(line, "w");
    fprintf(fp2,
 " Iamps, B0, c0i, c1i, c2i, c3i, c4i, c5i, c0e, c1e, c2e, c3e, c4e, c5e\n");
    for( int n=0; n<Nfrac; n++){
      Iamps[0] = Ia[n];
      //Set Enge coefs at Iamps
      p[0] = f_c0i( Iamps, p);
       p[1] = f_c1i( Iamps, p);
       p[2] = f_c2i( Iamps, p);
       p[3] = f_c3i( Iamps, p);
       p[4] = f_c4i( Iamps, p);
       p[5] = f_c5i( Iamps, p);
      p[6] = f_c0e( Iamps, p);
       p[7] = f_c1e( Iamps, p);
       p[8] = f_c2e( Iamps, p);
       p[9] = f_c3e( Iamps, p);
       p[10] = f_c4e( Iamps, p);
       p[11] = f_c5e( Iamps, p);
      //Set B0 at Iamps
      p[51] = f_B0( Iamps, p);
      //
      fprintf(fp2," %lg,",Iamps[0]);
      fprintf(fp2," %lg,",p[51]);
      for( int j=0; j<11; j++) fprintf(fp2," %lg,",p[j]);
      fprintf(fp2," %lg \n",p[11]);
      //Write EngeProd vs z to files for each Iamps
      sprintf(line, "EngeProd_%2.2f.txt", Ia[n]/Ilim);
      fp = fopen(line, "w");
      fprintf(fp, " %lg   Amps (Ilim=%lg) \n", Iamps[0], Ilim);
      fprintf(fp, " %lg   r0[m] \n", r0);
      fprintf(fp, " %i   Nz data points \n", Nz);
      fprintf(fp, " z[m]  EngeProd \n");
      printf(" I/Ilim= %04.2e  c0i= %06.2e  c0e= %06.2e, fp= %s\n",
        Ia[n]/Ilim, p[0],p[6], line);
      for( int i=0; i<Np; i++) func->SetParameter(i, p[i]);
      for( double z=zmin; z<=zmax; z+=delz )
      {
        zz[0] = z;
        ftmp = func_EngeProd(zz,p);
        fprintf(fp, " %lg   %lg \n", z, ftmp);
        //printf(" z= %06.2e  EngeProd= %06.2e \n", z, ftmp);
      }
      fclose(fp);
    }    
    fclose(fp2);
  }//gen_EngeProd_Z_I
  
  return;
}