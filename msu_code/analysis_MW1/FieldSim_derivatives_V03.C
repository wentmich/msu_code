/*
- Simulate multipole fields vs z and theta using COSY for magnets.
- Decompose harmonics at all z
- Evaluate at misalignment values dx,dy to estimate derivative of dipole
  components. These derivatives can be used to estimate misalignment of 
  magnets
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>
//# include "plot_and_fit_gradient.h"

typedef struct {
  int status;
  double r, t, z;
  double Br, Bt, Bz;
  double xc, yc; //shift of magnet center (-xc,-yc)
} SimData;

//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Measured
// Default values
#define N_HARMONIC (11)
#define N_AXIAL (801)
#define N_ANGLE (72)
#define N_TEMP (8)
#define MAX_FILENAME_LENGTH (512)
#define LLEN (512)

typedef struct {                                                        
  double T[N_TEMP];
  double z;
} zTemperatures;

typedef struct {
  int status; //0 if not filled properly
  double x, y, z, th; //expected values
  double xm, ym, zm, rm, thm; //actual values at measurement
  double vBr; //measured Hall probe voltage
  double Br; //evaluated field based on calibration factors
  //double radius, normradius, theta;
} MeasData;

  zTemperatures *zTemps;
  double *B_r_n, *A_r_n; //Nz,Nn for B=normal A=skew
  double *IFTfRe, *IFTfIm; //Nz
  double *b_n0_Re, *b_n0_Im, *a_n0_Re, *a_n0_Im; //Nz,Nn
  double *b_n0_Re_integ, *a_n0_Re_integ; //Nn
  double *Br_vs_th, *th_arr; //Nt
  double *z_arr, *fz_arr, *fz_arr_B, *fz_arr_A; //Nz
  double *k_arr, *FTfRe, *FTfIm; //Nk
  //MMMMMMMMMMMMMMMMMMM END

  //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv                               
  //Global variables
  int i, j, k, c, ee, n, m, index;
  double PI=3.141592654, PI2=2*PI;
  double DEGRAD=3.141592654/180.;

  // ** Field parameter settings
  const int NC  = 12;    //max number of field coefficients for each 2n-pole
  const int NEE = 2;     //entrance(0) exit(1)
  const int NNC = 11;    //max n harmonics for 2n poles for field params
  const double e = 2.71828;
  const char sep[2]=" ,";
  double LEFF = 0.782 ;  //Effective length (overriden by file value)
  int OrdCOSY=6;             //Order of COSY calculations

  double r_probes = 0.1873 ;  //radius at probes  (FSQ5, FST1)
//  double r_probes = 0.1173 ;  //radius at probes (FSQ2, FST2)
//  double r_probes = 0.0973 ;  //radius at probes (FSQ1)
//  double r_probes = 0.2 ;  //test
  // NEEded for simulation analysis
  double R0 = 0.2; //*0.8 ;   //reference radius of model (pole-tip radius)
  double R_AN = R0*0.8;       //radius of integrated harmonic analysis

  //Randomization parameters
  //Base error based on design characteristics
  // ..\Q_MapperProbeErrrors.xlsx
  int Nevents=5;                  //number of randomized events
  double zrms = 117.0e-6 ;          // z-error  
  double rrms =  62.0e-6 ;          // radial error
  double trms =  55.0e-6/r_probes ; // angle about ring
  double nrms =   3.1e-3 ;          // angle orientation
  double Br_rms = 2.0e-4;           // probe error in Tesla

  //double z_min = -0.9, z_max = 0.9, z_step = 0.005 ;
//TEMP--
  double z_min = -1.0, z_max = 1.0, z_step = 0.005;
  //double z_min =-0.1, z_max=0.1, z_step = 0.1;

  double th_min = 0., th_max = 2.*PI, th_step = 5.*DEGRAD;

  int time0, time1;
//  time0 = time( NULL );
//  printf( " Nz*Nt= %i  evaluated in (sec.): %i \n", Nz*Nt, (time1-time0));
//  printf( " points/sec = %lf \n", (double)(Nz*Nt)/(double)(time1-time0) );

  //Field array dependent
  double *fcoef;       //array with field coefficients
 // double *z_arr;
  const int Nz = 1+(z_max-z_min)/z_step ;
  int Nk = Nz;
  const int Nt = (int)( (th_max-th_min)/th_step)+1;
  const int Nn = 11;
  const int Nm = 11;
  //Arrays to keep generalized gradients of the field model
  //(Enge func. and such)
  double *b_nm, *a_nm;
  double *d1_b_nm, *d1_a_nm;
  //Nn, B-normal and A-skew fields for each multiple n as function of current I
  double *Bn, *An;
  double *d_b_nm, *d2_b_nm, *b_nm_1, *b_nm_tmp;
  double *d_a_nm, *d2_a_nm, *a_nm_1, *a_nm_tmp;


//MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM Measured
void Analyze_Field_Initialize()
{
  //initialize arrays 
  //if(( data=(MeasData *) calloc( Nt*Nz, sizeof(MeasData)) )==NULL){
  //  printf("ERROR in calloc\n"); exit(-1);}
  if(( zTemps=(zTemperatures *) calloc( Nz, sizeof(zTemperatures)) )==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  //Temporary 1D arrays used for integration
  if(( Br_vs_th=(double *) calloc( Nt, sizeof(double)) )==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if(( th_arr=(double *) calloc( Nt, sizeof(double)) )==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if(( z_arr=(double *) calloc( Nz, sizeof(double)) )==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((fz_arr=(double *) calloc(Nz,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((fz_arr_B=(double *) calloc(Nz,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((fz_arr_A=(double *) calloc(Nz,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}

  if((IFTfRe=(double *) calloc(Nz,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((IFTfIm=(double *) calloc(Nz,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}

  if((b_n0_Re=(double *) calloc(Nz*Nn,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((b_n0_Im=(double *) calloc(Nz*Nn,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((a_n0_Re=(double *) calloc(Nz*Nn,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((a_n0_Im=(double *) calloc(Nz*Nn,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}

  //For decomposed n harmonics from th dependent data
  if(( B_r_n=(double *) calloc( Nn*Nz, sizeof(double)) )==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if(( A_r_n=(double *) calloc( Nn*Nz, sizeof(double)) )==NULL){
    printf("ERROR in calloc\n"); exit(-1);}

  //For integrated generalized gradients
  if((b_n0_Re_integ=(double *) calloc(Nn,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((a_n0_Re_integ=(double *) calloc(Nn,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}

  //Initialize arrays that remain constant for analysis
  for(k=0; k<Nz; k++) z_arr[k] = (double)k * z_step + z_min;
  for(j=0; j<Nt; j++) th_arr[j] = (double)j * th_step/DEGRAD + th_min/DEGRAD;
  //printf(" th_step=%lf  th_min=%lf \n", th_step/DEGRAD, th_min/DEGRAD);
  //exit(0);
  return;
}

void harmonics_B_r_n__A_r_n_MeasData( MeasData *data ){
  //B_r_n, A_r_n: Decomposed data about th for each harmonic n up to (Nn-1)
  int index2;
  double th_min_rad = th_arr[0]*DEGRAD;
  double th_max_rad = th_arr[Nt-1]*DEGRAD;  
  for(n=0; n<Nn; n++){
   for(k=0; k<Nz; k++)
   //for(k=400; k<401; k++)
   {
    for(j=0; j<Nt; j++)
    {
      index = j + Nt*k;
      Br_vs_th[j] = data[index].Br;
      //if(n==0) printf(" j=%i, th=%lf  Br=%lf\n", j, th_arr[j], Br_vs_th[j]);
    }//j
    index2 = k + n*Nz;
    B_r_n[index2] = filon_tab_sin(Nt, Br_vs_th, th_min_rad, th_max_rad, n)/PI;
    A_r_n[index2] = filon_tab_cos(Nt, Br_vs_th, th_min_rad, th_max_rad, n)/PI;
    //if( n==3 ){
    // printf(" B filon k=%i, index2=%i  %lf\n", k, index2, B_r_n[index2]);
     //printf(" A filon k=%i, index2=%i  %lf\n", k, index2, A_r_n[index2]);
    //}
   }//k
  }//n
  //Optional: Write certain data to file for testing.
  /*
  FILE *fp;
  fp = fopen("B_r_n_test.txt","w");
  n=2;
  fprintf(fp,"n = %i\n", n);
  fprintf(fp,"  z   B_r_n \n");
   for(k=0; k<Nz; k++)
   {
    index2 = k + n*Nz;
    fprintf(fp," %i  %lf \n", z_arr[k], B_r_n[index2]);
   }//k
  fclose(fp);
  */
  return;
}

void harmonics_B_r_n__A_r_n_SimData( SimData *data,
      double *SB_r_n, double *SA_r_n ){
  //B_r_n, A_r_n: Decomposed data about th for each harmonic n up to (Nn-1)
  int index2;
  double th_min_rad = th_arr[0]*DEGRAD;
  double th_max_rad = th_arr[Nt-1]*DEGRAD;  
  for(n=0; n<Nn; n++){
   for(k=0; k<Nz; k++)
   //for(k=400; k<401; k++)
   {
    for(j=0; j<Nt; j++)
    {
      index = j + Nt*k;
      Br_vs_th[j] = data[index].Br;
      //if(n==2) printf(" j=%i, th=%lf  Br=%lf\n", j, th_arr[j], Br_vs_th[j]);
    }//j
    index2 = k + n*Nz;
    SB_r_n[index2] = filon_tab_sin(Nt, Br_vs_th, th_min_rad, th_max_rad, n)/PI;
    SA_r_n[index2] = filon_tab_cos(Nt, Br_vs_th, th_min_rad, th_max_rad, n)/PI;
    //if( n==3 ){
    // printf(" B filon k=%i, index2=%i  %lf\n", k, index2, B_r_n[index2]);
     //printf(" A filon k=%i, index2=%i  %lf\n", k, index2, A_r_n[index2]);
    //}
   }//k
  }//n
  /* Optional: Write certain data to file for testing.
  FILE *fp;
  fp = fopen("B_r_n_test.txt","w");
  n=1;
  fprintf(fp,"n = %i\n", n);
  fprintf(fp," k   z   B_r_n  A_r_n \n");
  printf(" k   z   B_r_n  A_r_n \n");
   for(k=0; k<Nz; k++)
   {
    index2 = k + n*Nz;
    fprintf(fp," %i  %lf  %le  %le \n", k, z_arr[k], SB_r_n[index2], SA_r_n[index2]);
    printf(" %i  %lf  %le  %le \n", k, z_arr[k], SB_r_n[index2], SA_r_n[index2]);
   }//k
  fclose(fp); */
  return;
}

void harmonics_BA_r_n_write( SimData *data,
 double *SB_r_n, double *SA_r_n, const char *fname)
{
// Write B_r_n and A_r_n decomposed harmonics to file
  int index2;
  FILE *fp;
  fp = fopen( fname, "w");
  fprintf(fp, " z, n, B_r_n, A_r_n \n");
  for(n=0; n<Nn; n++){
  fprintf(fp," \n");
   for(k=0; k<Nz; k++)
   {
    index = 0 + Nt*k; //to pic z value from first element of theta array
    index2 = k + n*Nz;
    if( n==1 || n==2 ){
     fprintf(fp, " %lf, %i, %lf, %lf \n",
       data[index].z, n, SB_r_n[index2], SA_r_n[index2] );
     //printf(" B filon k=%i, index2=%i  %lf\n", k, index2, B_r_n[index2]);
     //printf(" A filon k=%i, index2=%i  %lf\n", k, index2, A_r_n[index2]);
    }
   }//k
  }//n

  fclose(fp);
  return;
}

void derivative_BA_r_1_WRT_dx_dy( const char *fname, int n, 
  double dx, double dy,
  double *B_0_0, double *A_0_0,
  double *B_dx_0, double *A_dx_0,
  double *B_0_dy, double *A_0_dy )
// Given Normal & Skew (B & A) harmonic components evaluated at misaligned 
//  quad at dx and dy, evaluate derivatives at all z.
//  Writes B_r_2 since value are to be applied for quad misalignment.
{
  FILE *fp;
  fp = fopen( fname, "w");
  // B,A 2D arrays (n,z)
  int i1, i2;
  fprintf(fp, " z,  dB_r_1/dx,  dA_r_1/dy, B_r_2 \n");
  for(k=0; k<Nz; k++){
    i1 = k + 2*Nz; //index for n=2 to write B_r_2
    i2 = k + 1*Nz; //index for n=1 to write dB/dx and dA/dy
    fprintf(fp, " %lf, %lf, %lf, %lf \n", z_arr[k],
    (B_dx_0[i2]-B_0_0[i2])/dx,
    (A_0_dy[i2]-A_0_0[i2])/dy,
    B_0_0[i1] );
  }
  fclose(fp);
  return;
}

void zFourierTransfAnalysis(double *B_r_n, double *A_r_n)
{
  // Perform z-Fourier transform analysis and save results to corresponding 
  //  directory for data provided at each radius (Note: may be only one)
  int m_max, kmax, kmin, dk, n_harm;
  if( Nz==1 )
  {//11111111111111111111111111111111111111111
    printf(" zFourierTransfAnalysis: Nz=1. No zFourier to do.\n");
    printf("   Assuming no z-dependence; hence, only m=0 exists.\n");
  for(n=0; n<Nn; n++){
   n_harm = n+1; //n=0 not allowed in p_factor_radial.
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    fz_arr[k] = B_r_n[index];
   }//k
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    b_n0_Re[index] = B_r_n[index]/pow((r_probes/R0), n-1); //r_probes, R0
    b_n0_Im[index] = 0.;
    a_n0_Re[index] = A_r_n[index]/pow((r_probes/R0), n-1);
    a_n0_Im[index] = 0.;
   }//k
  }//n
    return;
  }//11111111111111111111111111111111111111111 end
  
  //Assuming Nz>=3
  double W=10; //Factor on width for optimal accuracy
  m_max = 20;
  kmax = 40*(W/z_max);
  kmin = -kmax;
  dk = (kmax-kmin)/(double)(Nk-1);
  // Arrays in k-space 
  if((k_arr=(double *) calloc(Nk,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((FTfRe=(double *) calloc(Nk,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  if((FTfIm=(double *) calloc(Nk,sizeof(double)))==NULL){
    printf("ERROR in calloc\n"); exit(-1);}
  for ( i=0; i<Nk; i++ ){
    k_arr[i] = ( (double)(Nk-i-1)*kmin+(double)(i)*kmax ) / (double)(Nk-1);
  }
  // ------   Br_B (Normal component) ---------------------------------
  for(n=0; n<Nn; n++){
  //for(n=2; n<3; n++){  printf(" **WARNING** Only n=2 being extracted.\n");
   n_harm = n+1; //n=0 not allowed in p_factor_radial./
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    fz_arr[k] = B_r_n[index];
   }//k
      //determine k-space solutions  Br,n(k)
      FTfRe = FourierTransCos( Nk, k_arr, Nz, z_arr, fz_arr);
      FTfIm = FourierTransSin( Nk, k_arr, Nz, z_arr, fz_arr);
      //determine k-space solutions  b_n,0(k)
      FTfRe = p_factor_radial( Nk, k_arr, n_harm, m_max, FTfRe, r_probes, R0);
      FTfIm = p_factor_radial( Nk, k_arr, n_harm, m_max, FTfIm, r_probes, R0);
      //determine z-space solutions b_n,0(z)
      IFTfRe = InvFourierTransRe(Nk,k_arr,FTfRe,FTfIm,Nz,z_arr);
      IFTfIm = InvFourierTransIm(Nk,k_arr,FTfRe,FTfIm,Nz,z_arr);
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    b_n0_Re[index] = IFTfRe[k];
    b_n0_Im[index] = IFTfIm[k];
   }//k
  }//n
  // ------   Br_A (Skew component) ---------------------------------
  for(n=0; n<Nn; n++){
  //for(n=2; n<3; n++){
   n_harm = n+1;
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    fz_arr[k] = A_r_n[index];
   }//k
      //determine k-space solutions  Ar,n(k)
      FTfRe = FourierTransCos(Nk,k_arr,Nz,z_arr,fz_arr);
      FTfIm = FourierTransSin(Nk,k_arr,Nz,z_arr,fz_arr);
      //determine k-space solutions  a_n,0(k)
      FTfRe = p_factor_radial( Nk, k_arr, n_harm, m_max, FTfRe, r_probes, R0);
      FTfIm = p_factor_radial( Nk, k_arr, n_harm, m_max, FTfIm, r_probes, R0);
      //determine z-space solutions a_n,0(z)
      IFTfRe = InvFourierTransRe(Nk,k_arr,FTfRe,FTfIm,Nz,z_arr);
      IFTfIm = InvFourierTransIm(Nk,k_arr,FTfRe,FTfIm,Nz,z_arr);
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    a_n0_Re[index] = IFTfRe[k];
    a_n0_Im[index] = IFTfIm[k];
   }//k
  }//n
  return;
}


void write_b_n0_a_n0_Re_and_Im( )
{
/* Michael understands this function. It creates eleven new files for the normal
and skew components of the on-axis gradient. These are created one for each
current I think. At zero current, the gradient should be zero everywhere. YES!
That was right, good job. */
  printf(",\n");
//wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww Re, Im
// Write Re and Im components to file
  FILE *outptr;
  char outname[512]; // buffer for file name
  int n_harm;

  //Write b_n,0(z) values to file
  for(n=0; n<Nn; n++){ // Nn is 11. Number of output files? yes I do believe it is. 
	  n_harm = n;
	  sprintf(outname,"b_n%02d_m0.txt", (int)n_harm); // puts file name in buffer
	  if((outptr=fopen(outname,"w"))==NULL){
		 printf("ERROR: can't open output file %s for write\n",outname);
		 exit(-1);
	  }
	  fprintf(outptr,"      z,       b_n0_Re,       b_n0_Im\n");
	   for(k=0; k<Nz; k++)
	   {
		index = k + n*Nz;
		fprintf(outptr," %14.6lg, %14.6lg, %14.6lg \n",
		  z_arr[k], b_n0_Re[index], b_n0_Im[index]);
	   }//k
	  fclose(outptr);  
  }//n
  //Write a_n,0(z) values to file
  for(n=0; n<Nn; n++){
  n_harm = n+1;
  sprintf(outname,"a_n%02d_m0.txt", (int)n_harm);
  if((outptr=fopen(outname,"w"))==NULL){
     printf("ERROR: can't open output file %s for write\n",outname);
     exit(-1);
  }
  fprintf(outptr,"      z,       a_n0_Re,       a_n0_Im\n");
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    fprintf(outptr," %14.6lg %14.6lg %14.6lg \n",
      z_arr[k], a_n0_Re[index], a_n0_Im[index]);
   }//k
  fclose(outptr);  
  }//n
  return;
//wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww Re, Im
}

void rescale_b_n0_a_n0( double refR_old, double refR_new,
       double r_old, double r_new)
// refR = reference radius (usually stays constant)
// r_ = radius at which field is desired. r_new = r where field are desired at
{
  for(n=0; n<Nn; n++){
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    b_n0_Re[index] *= pow((r_old/r_new)*(refR_new/refR_old),(n-1));
   }//k
  }//n
  return;
}

void integrate_gen_gradients( char *fname_integrals )
//Integrate generalized gradients for each harmonic.
// This is the result of integrating at r=R0. 
{
  for(n=0; n<Nn; n++){
   for(k=0; k<Nz; k++)
   {
    index = (int) k + n*Nz;
    fz_arr_B[k] = b_n0_Re[index];
    fz_arr_A[k] = a_n0_Re[index];
   }//k
   if( Nz==0 ) {
     printf(" integrate_gen_gradients: Error. Nz=0\n"); exit(0);
   } else if( Nz==1 ) {
     printf(" integrate_gen_gradients: Error. Nz=1\n"); exit(0);
   } else {
     printf(" Nz= %i \n",Nz); exit(0);
     b_n0_Re_integ[n] = integd_trapez( Nz, z_arr, fz_arr_B);
     a_n0_Re_integ[n] = integd_trapez( Nz, z_arr, fz_arr_A);
   }
  }//n
  if( strlen(fname_integrals)>0 )
  {//---------
   FILE *outptr;
   printf("  write integrated harmonics to file: %s\n", fname_integrals);
   if((outptr=fopen(fname_integrals,"w"))==NULL){
      printf("ERROR: can't open output file %s for write\n",fname_integrals);
      exit(-1);
   }
   printf(        " n       b_n0_Re_int    a_n0_Re_int \n");
   fprintf(outptr," n       b_n0_Re_int    a_n0_Re_int \n");
   for(n=0; n<Nn; n++){
     fprintf(outptr," %i %14.6lg %14.6lg \n",
         n, b_n0_Re_integ[n], a_n0_Re_integ[n]);
   printf(        " %i %14.6lg %14.6lg \n",
       n, b_n0_Re_integ[n], a_n0_Re_integ[n]);
   }
   fclose(outptr);
  } else {
   printf(        " n       b_n0_Re_int    a_n0_Re_int \n");
   for(n=0; n<Nn; n++){
    printf(        " %i %14.6lg %14.6lg \n",
       n, b_n0_Re_integ[n], a_n0_Re_integ[n]);
   }
  }//strlen--------
  return;
}

void integrate_gradients_r_factor( double r_factor, char *fname_integrals )
//Integrates B_r_n, A_r_n about (r/R0)=r_factor
// For example, if integral factor at 80% of R0 is needed 
{
  for(n=0; n<Nn; n++){
   for(k=0; k<Nz; k++)
   {
    index = k + n*Nz;
    fz_arr_B[k] = b_n0_Re[index]*pow( r_factor, (n-1) );
    fz_arr_A[k] = a_n0_Re[index]*pow( r_factor, (n-1) );;
   }//k
   if( Nz==0 ) {
     printf(" integrate_gradients_r_factor: Error. Nz=0\n"); exit(0);
   } else if( Nz==1 ) {
     //printf(" integrate_gradients_r_factor: Nz=1. No integration needed.\n");
     b_n0_Re_integ[n] = fz_arr_B[0];
     a_n0_Re_integ[n] = fz_arr_A[0];
   } else {
     printf(" Nz= %i \n",Nz);
     b_n0_Re_integ[n] = integd_trapez( Nz, z_arr, fz_arr_B);
     a_n0_Re_integ[n] = integd_trapez( Nz, z_arr, fz_arr_A);
   }
  }//n
  if( strlen(fname_integrals)>0 )
  {//---------
   FILE *outptr;
   printf("  write integrated harmonics to file: %s\n", fname_integrals);
   if((outptr=fopen(fname_integrals,"w"))==NULL){
      printf("ERROR: can't open output file %s for write\n",fname_integrals);
      exit(-1);
   }
   printf(        " n       b_n0_Re_int    a_n0_Re_int \n");
   fprintf(outptr," n       b_n0_Re_int    a_n0_Re_int \n");
   for(n=0; n<Nn; n++){
     fprintf(outptr," %i %14.6lg %14.6lg \n",
         n, b_n0_Re_integ[n], a_n0_Re_integ[n]);
   printf(        " %i %14.6lg %14.6lg \n",
       n, b_n0_Re_integ[n], a_n0_Re_integ[n]);
   }
   fclose(outptr);
  } else {
   printf(        " n       b_n0_Re_int    a_n0_Re_int \n");
   for(n=0; n<Nn; n++){
    printf(        " %i %14.6lg %14.6lg \n",
       n, b_n0_Re_integ[n], a_n0_Re_integ[n]);
   }
  }//strlen--------
  return;
}

void read_integ_file(char *filename, int Nharm, int Nfile, 
  double b_arr[], double a_arr[] )
{
  int n;
  double d1, d2;
  printf("read_integ_file: %s\n",filename);
  FILE *fp; char buffer[512]; int j;
  fp = fopen( filename, "r");
  if (fp==NULL) { printf("no file, %s. Aborting.\n", filename); exit(-1); }
  fgets(buffer, 512, fp); //cout << buffer;
  for( j=0; j<Nharm; j++ ){
    fgets(buffer, 512, fp); //cout << buffer;
    index = j + Nfile*Nharm;
    sscanf(buffer," %i %lf %lf\n", &n, &b_arr[index], &a_arr[index]);
    //sscanf(buffer," %i %lf %lf\n", &n, &d1, &d2);
    if( j != n ) { printf("read_integ_file: n mismatch !\n"); exit(-1); }
    //b_arr[index] = d1; a_arr[index] = d2;
  //  printf("  %i %i %i %lf %lf\n", index, j, n, b_arr[index], a_arr[index]);
  }
  fclose( fp );
  return;
}
//MMMMMMMMMMMMMMMMMMM END

void ROOT_plot_contour( TH2D *Hist2D, TCanvas *pCanv ){
  //Plot them
   gStyle->SetPalette(1);
   gStyle->SetOptStat(0);
   gStyle->SetTitleW(0.99);
   gStyle->SetTitleH(0.08);
   const Int_t NCont=20;
   Double_t contours[NCont];
   Double_t cont_max = Hist2D->GetMaximum() ;
   Double_t cont_min = Hist2D->GetMinimum() ;
   Double_t dcont = (cont_max-cont_min)/(NCont-1);
   for( Int_t i = 0; i<NCont; i++ ){
     contours[i] = cont_min + dcont*(Double_t)i ;
   }
   Hist2D->SetContour(NCont, contours);
   // Draw contours as filled regions, and Save points
   pCanv->SetRightMargin(0.15);
   pCanv->SetTopMargin(0.15);
   Hist2D->Draw("CONT Z LIST");
   pCanv->Update();
   pCanv->SaveAs("ROOT_plot_contour.png");
  return;
}

void Plot_Br_Bt_Bz_VS_z_th( SimData *data )
{
  //Optional:
  //Plot histograms
  if( (TH2D*)gROOT->FindObject("pHistBr") ) delete pHistBr;
  TH2D *HistBr = new TH2D("pHistBr","Br(r,th)",
    Nz, z_min, z_max, Nt, th_min/DEGRAD, th_max/DEGRAD);
  if( (TH2D*)gROOT->FindObject("pHistBt") ) delete pHistBt;
  TH2D *HistBt = new TH2D("pHistBt","Bt(r,th)",
    Nz, z_min, z_max, Nt, th_min/DEGRAD, th_max/DEGRAD);
  if( (TH2D*)gROOT->FindObject("pHistBz") ) delete pHistBz;
  TH2D *HistBz = new TH2D("pHistBz","Bz(r,th)",
    Nz, z_min, z_max, Nt, th_min/DEGRAD, th_max/DEGRAD);
  //Fill TH2D's
  for( k=0; k<Nz; k++ ){
   for( j=0; j<Nt; j++ ){
     index = j + Nt*k ;
     HistBr->SetBinContent( k, j, data[index].Br);
     HistBt->SetBinContent( k, j, data[index].Bt);
     HistBz->SetBinContent( k, j, data[index].Bz);
   }
  }
  if( (TCanvas*)gROOT->FindObject("pCanvBr") ) delete pCanvBr;
  TCanvas* pCanvBr = new TCanvas("pCanvBr","Br(z,th)",0,0,600,600);
  ROOT_plot_contour( HistBr, pCanvBr );
  /*
  if( (TCanvas*)gROOT->FindObject("pCanvBz") ) delete pCanvBz ;
  TCanvas* pCanvBz = new TCanvas("pCanvBz","Bz(z,th)",100,100,600,600);
  ROOT_plot_contour( HistBz, pCanvBz );
  if( (TCanvas*)gROOT->FindObject("pCanvBt") ) delete pCanvBt;
  TCanvas* pCanvBt = new TCanvas("pCanvBt","Bth(z,th)",50,50,600,600);
  ROOT_plot_contour( HistBt, pCanvBt );
  */
  return;
}

void Plot_Br_VS_z_th( SimData *data )
{
  //Optional:
  //Plot histograms
  if( (TH2D*)gROOT->FindObject("pHistBr") ) delete pHistBr;
  TH2D *HistBr = new TH2D("pHistBr","Br(r,th)",
    Nz, z_min, z_max, Nt, th_min/DEGRAD, th_max/DEGRAD);
  //Fill TH2D's
  for( k=0; k<Nz; k++ ){
   for( j=0; j<Nt; j++ ){
     index = j + Nt*k ;
     HistBr->SetBinContent( k, j, data[index].Br);
   }
  }
  if( (TCanvas*)gROOT->FindObject("pCanvBr") ) delete pCanvBr;
  TCanvas* pCanvBr = new TCanvas("pCanvBr","Br(z,th)",0,0,600,600);
  ROOT_plot_contour( HistBr, pCanvBr );
  return;
}


void ROOT_TGraph_n_y1_y2(Int_t Nh, double *y1, char *title1,
  double *y2, char *title2 )
{
// create or regenerate graph, depending on existence
  double ymax = y1[0], ymin = ymax;
  for(Int_t n=0; n<Nh; n++ ){
    if( ymin > y1[n] ) ymin = y1[n];
    if( ymax < y1[n] ) ymax = y1[n];
    if( ymin > y2[n] ) ymin = y2[n];
    if( ymax < y2[n] ) ymax = y2[n];
  }
  ymax *= 1.1;
  const Int_t NNh=Nh;
  double x[NNh];
  for( Int_t i=0; i<Nh; i++ ) x[i] = (double)i;
  if( (TCanvas*)gROOT->FindObject("c1") ){ //regenerate
    c1->Delete();
  	TCanvas *c1 = new TCanvas("c1","Graph Draw Options",200,10,1200,600);
  } else { //create
  	TCanvas *c1 = new TCanvas("c1","Graph Draw Options",200,10,1200,600);
  }
  c1->Divide(2,1);
    gStyle->SetOptStat(0);
    gPad->SetGrid();
  c1->cd(1);
  TGraph *gr1 = new TGraph(Nh, x, y1);
  gr1->Draw("APB");
    gr1->SetTitle(title1);
    gr1->SetMarkerStyle(0);
    gr1->GetYaxis()->SetRangeUser(ymin,ymax);
    gr1->SetMarkerColor(4);
    gr1->SetFillColor(4);
    gr1->SetMarkerSize(1.5);
  c1->cd(2);
  TGraph *gr2 = new TGraph(Nh, x, y2);
  gr2->Draw("APB");
    gr2->SetTitle(title2);
    gr2->SetMarkerStyle(0);
    gr2->GetYaxis()->SetRangeUser(ymin,ymax);
     gr2->SetMarkerColor(4);
     gr2->SetFillColor(4);
     gr2->SetMarkerSize(1.5);
  c1->SaveAs("Plot_Harmonics.png");
  return;
}

void  harmonics_BA_r_n_Plot_at_z( const double z_select,
  double *B_r_n, double *A_r_n, const char *fname)
{
  FILE *fp;
  fp = fopen( fname, "w");
  int k_select = 1+(z_select-z_min)/z_step ;
  printf(" plot of harmonics at z=%lf \n", z_select);
  printf(" Nz=%i  k=%i \n", Nz, k);
  printf(" See file: %s\n", fname);
  double *B_r_n__z; B_r_n__z = (double*) malloc( Nn*sizeof(double) );
  double *A_r_n__z; A_r_n__z = (double*) malloc( Nn*sizeof(double) );
  fprintf(fp, " plot of harmonics at z=%lf \n", z_select);
  fprintf(fp, " n,  B_r_n,  A_r_n \n");
  for(n=0; n<Nn; n++){
   for(k=0; k<Nz; k++){
    index = k + n*Nz;
    if( k==k_select ){ 
      B_r_n__z[n] = B_r_n[index];
      A_r_n__z[n] = A_r_n[index];
      fprintf(fp, " %i, %lf, %lf \n", n, B_r_n__z[n], A_r_n__z[n]);
    }
   }
  }
  fclose( fp );
  ROOT_TGraph_n_y1_y2( Nn, 
    B_r_n__z, "Normal components; n; B_r_n",
    A_r_n__z, "  Skew components; n; A_r_n");
  free( B_r_n__z );
  free( A_r_n__z );
  return;
}

void Field_R0_r_t_z__Br_Bt_Bz(
  double R0, double r, double th, double z,
  double &Br, double &Bt, double &Bz )
{
  //th in radian units
  double r_R0=r/R0, nn, mm, tmp;
  int ii, denom, kk;
  kk = (z-z_min)/z_step ;
  Br=0.; Bt=0.; Bz=0.;
  for( n=1; n<Nn; n++) { //Nz
    for( m=0; m<Nm; m++) { //Nm
      ii = kk + Nz*( m + Nm*n ) ;
      nn = (double)n;
      mm = (double)m;
      denom = (n+2*m);
      tmp = (nn+2.0*mm);
      Br += ( a_nm[ii]*cos(nn*th) + b_nm[ii]*sin(nn*th) )*pow(r_R0,n+2*m-1);
      if( denom!=0 ){
        //printf(" n%i,m%i  %lf \n", n,m, nn/tmp);
        Bt += (-a_nm[ii]*sin(nn*th) + b_nm[ii]*cos(nn*th) )*pow(r_R0,n+2*m-1)*nn/tmp;
        Bz += ( d1_a_nm[ii]*cos(nn*th) + d1_b_nm[ii]*sin(nn*th) )*pow(r_R0,n+2*m)*R0/tmp;
      } else printf("Error: (n+2m)=0  n=%i m=%i \n", n, m);
    }
  }
  //printf("  z, kk = %lf %i \n",z,kk);
  return;
}

void read_M5_params_NORM( char *fname, double *Bn_)
{
// Michael understands this section of the code. It will read in all of the
// initial values from fname and will also get the reference radius, the effective
// length, and the field at R0 for each multipole component.
  //Read fields at R0 for each NORMAL multipole component (no skews for M5)
  //  scan line of data from line into array  printf(" read_M5_params_NORM: %s \n", fname);
  int i;
  FILE *fp; char line[512]; char *pch; //const char sep[2]=" ,";
  // creates a dummy file, a line buffer, a token buffer, and a delimiter variable
  fp = fopen(fname, "r");
  //opens the file fname to read it and places it in the dummy file
  //Read header comments
  sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
  //R0: 
  //fgets(line,512,fp);
  sscanf(line," %lf\n", &R0);
  // gets the reference radius
  //LEFF: 
  fgets(line,512,fp);
  sscanf(line," %lf\n", &LEFF);
  // gets the effective length
  sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
  // makes first character of "line" a pound sign. Then says, while the first 
  //character of "line" is a pound sign, get the next line of the file. 
  // in the end, "line will be the first line of the file w/o a "#" at the beggining
    i=2;
    pch = strtok(line, sep); sscanf(pch,"%lf", &Bn_[i]);
    for( i=3; i<7 ; i++ ){
      pch = strtok(NULL, sep); sscanf(pch,"%lf", &Bn_[i]);
    }
  // this part above just reads in the fields at R0 for each multipole from the file
  // it looks complicated due to the syntax of the file in which the parameters are kept
  //Read Enge coefs ( see read_fcoef() 
  sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
  while( !feof(fp) ){
  	// while the end of file has not been reached  
    //scan line of data into array
    pch = strtok(line, sep); sscanf(pch,"%i", &n);
    // gets first column
    pch = strtok(NULL, sep); sscanf(pch,"%i", &ee);
    // gets second column
    for( c=0; c<NC; c++){
      // gets all twelve paramters from the line its referencing
      index = c + NC*( NEE*n + ee );
      // makes sure parameter is going to the right place
      pch = strtok(NULL, sep); sscanf(pch,"%lf", &fcoef[index]);
      //printf(" %i \n", index);
    }
    fgets(line,512,fp);
  }
  // and the part above just reads in all the data from the Enge coefficients. 
  // again, complications arise from the syntax of the file
    //scan LAST line of data into array
    pch = strtok(line, sep); sscanf(pch,"%i", &n);
    pch = strtok(NULL, sep); sscanf(pch,"%i", &ee);
    for( c=0; c<NC; c++){
      // this loop will read in each one of the parameters from the file
      // the total is 168 parameters indexed 0-167. 
      index = c + NC*( NEE*n + ee );
      pch = strtok(NULL, sep); sscanf(pch,"%lf", &fcoef[index]);
      //printf(" %i \n", index);
    }
  /* Optional: Inspect coefs
  printf(" n  ee  fcoef[c]... \n");
  for( n=0; n<NNC; n++){
    for( ee=0; ee<NEE; ee++){
      printf(" %i  %i ", n, ee);
      for( c=0; c<NC; c++){
        index = c + NC*( NEE*n + ee );
        if(c<(NC-1)) { printf(" %lf", fcoef[index]); }
        else { printf(" %lf\n", fcoef[index]); }
      }//c
    }//ee
  }//n
  exit(0);
  */
  
  return;
}

void read_multipole_fields( char *fname, int Nn_, double Bn_[], double An_[] )
{
  printf(" read_multipole_fields: %s \n", fname);
  int i;
  FILE *fp; char line[512]; char *pch; const char sep[2]=" ,";
  fp = fopen(fname, "r");
  //Read header comments
  sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
    //Normal components: scan line of data from line into array
    i=0;
    pch = strtok(line, sep); sscanf(pch,"%lf", &Bn_[i]);
    for( i=1; i<Nn_ ; i++ ){
      pch = strtok(NULL, sep); sscanf(pch,"%lf", &Bn_[i]);
    }
    //Skew components: 
    fgets(line,512,fp);
    i=0;
    pch = strtok(line, sep); sscanf(pch,"%lf", &An_[i]);
    for( i=1; i<Nn_ ; i++ ){
      pch = strtok(NULL, sep); sscanf(pch,"%lf", &An_[i]);
    }
    //LEFF: 
    fgets(line,512,fp);
    sscanf(line," %lf\n", &LEFF);
  return;
}

void read_fcoef( char *fname )
/* Read field parameters, such as Enge coefficients, that describe the 
     field of multipole for each 
       n = harmonic number (i.e. 2n poles)
       ee = 0-entrance side  1-exit side
       Below is an example for reading n=1 and n=2 with 12 coefficients.
       The number of coefficients NC must be the same for each line.
       Note that NC=12 for this example. 0-5 are for the Enge coefficients.
       The rest can be used for the so-called extended functions
       (e.g. Gaussians). Lines at the beginning that begin with a # character
       are treated as comments.
# n  ee  c0 c1 c2 c3 c4 c5  c6 c7 c8  c9 c10 c11
  1  0   0. 2.8 0. 0. 0. 0.   0 0 0  0 0 0
  1  1   0. 2.8 0. 0. 0. 0.   0 0 0  0 0 0
  2  0   0. 3. 0. 0. 0. 0.   0 0 0  0 0 0
  2  1   0. 3. 0. 0. 0. 0.   0 0 0  0 0 0
*/
{
  printf(" read_fcoef: %s \n", fname);
  FILE *fp; char line[512]; char *pch; const char sep[2]=" ,";
  fp = fopen(fname, "r");
  //Read header comments
  sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
  while( !feof(fp) ){
    //scan line of data into array
    pch = strtok(line, sep); sscanf(pch,"%i", &n);
    pch = strtok(NULL, sep); sscanf(pch,"%i", &ee);
    for( c=0; c<NC; c++){
      index = c + NC*( NEE*n + ee );
      pch = strtok(NULL, sep); sscanf(pch,"%lf", &fcoef[index]);
      //printf(" %i \n", index);
    }
    fgets(line,512,fp);
  }
    //scan LAST line of data into array
    pch = strtok(line, sep); sscanf(pch,"%i", &n);
    pch = strtok(NULL, sep); sscanf(pch,"%i", &ee);
    for( c=0; c<NC; c++){
      index = c + NC*( NEE*n + ee );
      pch = strtok(NULL, sep); sscanf(pch,"%lf", &fcoef[index]);
      //printf(" %i \n", index);
    }
  /*Optional:
  printf(" n  ee  fcoef[c]... \n");
  for( n=0; n<NNC; n++){
    for( ee=0; ee<NEE; ee++){
      printf(" %i  %i ", n, ee);
      for( c=0; c<NC; c++){
        index = c + NC*( NEE*n + ee );
        if(c<(NC-1)) { printf(" %lf", fcoef[index]); }
        else { printf(" %lf\n", fcoef[index]); }
      }//c
    }//ee
  }//n
  */
  return;
}

double fcoef_func( int n_, int ee_, int c_ )
{
  double coef; int i;
  i = c_ + NC*( NEE*n_ + ee_ );
  coef = fcoef[i];
  return coef;
}

void Field_Simulation_via_COSY_M5( SimData *data )
{
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  //  Field evaluation from Enge & other ceofficients via COSY
  //  data    = array with (x y z) points
  //  FS_fname = temp output file for COSY output (x y z Bx By Bz)
  //
  char FS_fname[50]="TMP_COSY_out.txt" ;
  double xx, yy, zz, Bxx, Byy, Bzz, Brr, Btt, tt;
  //
  bool call_COSY=1;
  if( call_COSY ){
  //-- First write parameters to file for COSY to read
  FILE *fpC;
  fpC = fopen("TMP_Enge_params_for_COSY.txt", "w");
  fprintf(fpC, " %d\n", OrdCOSY);
  fprintf(fpC, " %lf\n", LEFF);
  fprintf(fpC, " %lf\n", R0);
  for( n=1 ; n<Nn; n++ ){
   for( ee=0; ee<NEE; ee++ ){
    fprintf(fpC, " %i\n", n);
    fprintf(fpC, " %i\n", ee);
    for( c=0; c<NC; c++ ){
      index = c + NC*( NEE*n + ee );
      fprintf(fpC, " %lf\n", fcoef[index]);
    }
   }
  }
  fprintf(fpC,"Normal field components\n");
  for( i=0; i<Nn ; i++ ) fprintf(fpC," %lf \n", Bn[i]);
  fprintf(fpC,"Skew field components\n");
  for( i=0; i<Nn ; i++ ) fprintf(fpC," %lf \n", An[i]);
  fprintf(fpC,"Points for field evaluation\n");
  fprintf(fpC,"%s\n",FS_fname);
  //Write positions to evaluate fields at
  for( k=0; k<Nz; k++ ){ //Nz
   for( j=0; j<Nt; j++ ){   //Nt
//  for( k=450; k<451; k++ ){ //Nz
//   for( j=0; j<5; j++ ){   //Nt
     index = j + Nt*k ;
     xx = -data[index].xc + data[index].r * cos( data[index].t );
//printf("%lf\n",data[index].r); exit(0);
     yy = -data[index].yc + data[index].r * sin( data[index].t );
     fprintf(fpC, " %lf,%lf,%lf\n", xx, yy, data[index].z);
   }
  }
  fprintf(fpC,"END\n");
  fclose(fpC);

  //-- Next, call COSY to read parameters and positions
  //   COSY will evaluate field at each point using GENERATE_M5_FIELD_EXT
  time0 = time( NULL );
  system("runcosy91.sh COSY_field_Gen.fox");
  time1 = time( NULL );
  printf( " Nz*Nt= %i  evaluated in (sec.): %i \n", Nz*Nt, (time1-time0));
  printf( " points/sec = %lf \n", (double)(Nz*Nt)/(double)(time1-time0) );
  printf( " sec/point = %lf \n", (double)(time1-time0)/(double)(Nz*Nt) );
  }

  //Read back data from COSY output and put into data
  FILE *fpCO;
  char line[512]; int Nscan;
  fpCO = fopen(FS_fname, "r");
  if (fpCO==NULL) { printf("File not found: %s\n",FS_fname); exit(0); }  
  for( k=0; k<Nz; k++ ){ //Nz
   for( j=0; j<Nt; j++ ){   //Nt
//  for( k=450; k<451; k++ ){ //Nz
//   for( j=0; j<5; j++ ){   //Nt
     index = j + Nt*k ;
     /*
     xx = data[index].r * cos( data[index].t );
     yy = data[index].r * sin( data[index].t );
     fprintf(fpC, " %lf,%lf,%lf\n", xx, yy, data[index].z);
     */
     fgets( line, 512, fpCO );
     Nscan = sscanf( line, " %lf %lf %lf %lf %lf %lf\n",
      &xx,&yy,&zz, &Bxx,&Byy,&Bzz);
     //printf( " %lf %lf %lf %lf %lf %lf\n", xx,yy,zz, Bxx,Byy,Bzz);
     if( Nscan!=6 ){
       printf(" Error reading point: k%i j%i \n",k,j);
       exit(0);
     }
     //SimData: r,t,z,Br,Bt,Bz
     //tt = atan2(yy,xx);
     data[index].Br = Bxx*cos(data[index].t) + Byy*sin(data[index].t);
     data[index].Bt =-Bxx*sin(data[index].t) + Byy*cos(data[index].t);
     data[index].Bz = Bzz;
   }
  }
  fclose( fpCO );
  return;
  //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc END
}

void File_Save_format_z_Br_MAP_V03( char *fname, SimData *data )
{
  //Emulates mapper format.
  //Write data to file in z_Br format (i.e. format expected from Br vs z mapper)
  FILE *fzBr;
  fzBr = fopen(fname,"w");
  fprintf(fzBr,"# ---------- Simulated data ------------- \n");
  fprintf(fzBr,"# theta scan settings (must have constant th_step)\n");
  fprintf(fzBr," %lf  %lf \n", th_min/DEGRAD, th_step/DEGRAD);
  fprintf(fzBr,"# z-position scan settings (must have constant zstep)\n");
  fprintf(fzBr," %lf  %lf  %lf\n", z_min, z_max, z_step);

  fprintf(fzBr,"# File with probe calibration factors and other data\n");
  fprintf(fzBr,"NO_FILE_FOR_PROBE_DATA\n");
  fprintf(fzBr,"# File laser scanner data\n");
  fprintf(fzBr,"NO_FILE_FOR_LASER_SCANNER_DATA\n");
  fprintf(fzBr,"# radii: R0, r_probes, R_AN\n");
  fprintf(fzBr," %lf \n", r_probes);  // probe radius
  double Vm;
  for( k=0; k<Nz; k++ ){
  fprintf(fzBr," %lf \n", time(NULL) );  // timestamp
  fprintf(fzBr," 25.1 25.2 25.3 25.4 25.5 25.6 25.7 25.8 \n");
  index = 0 + Nt*k ;
  fprintf(fzBr," %lf \n", data[index].z);
   for( j=0; j<(Nt-1); j++ ){
     index = j + Nt*k ;
     Vm = 0.2*data[index].Br ; //Assuming simple linear relation
     //Vm = (double)j*th_step/DEGRAD; //for test
     fprintf(fzBr," %lf  %lf \n", Vm, data[index].Br);
   }
  }
  fprintf(fzBr,"END\n");
  // Other parameters that are needed by FourierBR01.cpp
  fprintf(fzBr," %lf \n", R0);        //reference radius of magnet
  fprintf(fzBr," %lf \n", R_AN);      //reference radius for analysis  
  fclose(fzBr);
  return;
}

void read_z_Br_MAP_V3( char *fdata, MeasData *data)
{
  // Data structure as of 2018-06-27
  // Assumes data is in form where only z value and measurements at each probe 
  //  about the circle are measured.
  printf("    read data points ................ \n ");
  //Variables that are global in V2
  int Nt_data, Nt, Nz;
  double z_min, z_max, z_step, th_min, th_max, th_step;
  double timestamp_d;        //time stamp value
  char timestamp_c[200];     //time stamp string

  FILE *fp;
  int TheEnd, m, test, NL=0, valid_points;
  char file_probe_data[LLEN];
  char file_laser_scan_data[LLEN];
  double tmp[20];
  char buffer1[LLEN];
  fp = fopen(fdata, "r");
  printf("read file: %s\n",fdata);
  if (fp==NULL) { printf("no file. Aborting.\n"); exit(0); }
  //fgets(buffer1,LLEN,fp); printf("%s",buffer1); //title line

  //___ theta properties
  sprintf(buffer1,"#");
  while(strncmp(buffer1,"#",1)==0){fgets(buffer1,LLEN,fp); NL++;}
    sscanf(buffer1," %lf %lf\n",&th_min, &th_step);
  //Note: Nt must be odd for integration purposes. N_data can be one less than
  //  Nt since first point of circle does not have to be measured twice.
  //  The missing point is to be filled later after reading data.
  th_max = th_min + 360.;
  Nt_data = (th_max-th_min)/th_step ;
  Nt = Nt_data+1 ;
  printf(" Nt_data=%i, th_min=%lf, th_max=%lf, th_step=%lf (deg.) \n",
    Nt_data, th_min, th_max, th_step);
  printf(" Nt (must be odd) = %i \n", Nt);

  //___ z properties
  sprintf(buffer1,"#");
  while(strncmp(buffer1,"#",1)==0){fgets(buffer1,LLEN,fp); NL++;}
    sscanf(buffer1," %lf %lf %lf\n",&z_min, &z_max, &z_step);
  Nz = 1+(z_max-z_min)/z_step ;
  printf(" Nz=%i, z_min=%lf, z_max=%lf, z_step=%lf \n",
    Nz, z_min, z_max, z_step);

  //___ file name of file with for probe arrangement data
  sprintf(buffer1,"#");
  while(strncmp(buffer1,"#",1)==0){fgets(buffer1,LLEN,fp); NL++;}
    sscanf(buffer1,"%s\n", &file_probe_data);
    printf("  probe properties: %s\n",file_probe_data);

  //___ file name of file with laser scanner data
  sprintf(buffer1,"#");
  while(strncmp(buffer1,"#",1)==0){fgets(buffer1,LLEN,fp); NL++;}
    sscanf(buffer1,"%s\n", &file_laser_scan_data);
    printf("  laser scan data: %s\n",file_laser_scan_data);

  //___ radius of probes 
  sprintf(buffer1,"#");
  while(strncmp(buffer1,"#",1)==0){fgets(buffer1,LLEN,fp); NL++;}
   sscanf(buffer1,"%lf\n",&r_probes);

  int index0, indexNt, cnt;
  double x, y, z, r, th, Bx, By, Bz, Br, Bt;
  double xm, ym, zm, rm, thm;

  // Initialize data points as invalid (.status=0)  
  cnt=0;
    for(k=0; k<Nz; k++){
      z = (double)k * z_step + z_min;                                            
      for(j=0; j<Nt; j++){
        th = (double)j * th_step;
        index = j + Nt*k;
        x = r_probes*cos(th*DEGRAD);
        y = r_probes*sin(th*DEGRAD);
	      data[index].x=x;
   	    data[index].y=y;
  	    data[index].z=z;
   	    data[index].th=th;
	      data[index].status=0; // make point invalid
	      cnt++;
      }// k
    }// j
  printf("    total data points initialized, cnt = %i \n\n", cnt);
  // Initialize values of zTemps(z).T(i)
  for(k=0; k<Nz; k++){
    for( i=0; i<N_TEMP; i++){
      zTemps[k].T[i] = 0.;
      zTemps[k].z = (double)k * z_step + z_min;
    }
    //printf("%f %f %f\n",zTemps[k].z, zTemps[k].T[3],zTemps[k].T[7]);
  }

  //Read the rest of the data in the file  
  while(!TheEnd)
  {
    fgets(buffer1,LLEN,fp); NL++;
    if( feof(fp) || strncmp(buffer1,"END",3)==0 )
    {
      TheEnd=1; printf(" _______ Encounted END as expected in file.\n");
      continue ;
    }
    else
    {
      //2a___ time stamp
      fgets(buffer1,LLEN,fp); NL++;
      if( (test=sscanf(buffer1," %lf\n", &timestamp_d))!=1)
      { printf("%i: Error in reading timestamp_d.\n",NL);
        cout << buffer1 ; exit(-1);
      }
      //2b___ read zTemps
      if( (test=sscanf(buffer1," %lf %lf %lf %lf %lf %lf %lf %lf\n",
        &tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5],&tmp[6],&tmp[7]))!=8)
      { printf("%i: Error in reading temperatures.\n",NL);
        cout << buffer1 ; exit(-1);
      }
      //2c___ read z position
      fgets(buffer1,LLEN,fp); NL++;
      if( (test=sscanf(buffer1," %lf\n",&zm))!=1)
      { printf("%i: Error in reading zm.\n",NL);
        cout << buffer1 ; exit(-1);
      }
      //define z index
      k = (int)floor( 0.5+ (zm-z_min)/z_step );
      //printf(" zm, k = %lf  %i \n", zm, k);
      for( m=0; m<N_TEMP; m++) zTemps[k].T[m]=tmp[m];
      //check: for( i=0; i<N_TEMP; i++) cout << zTemps[k].T[i] <<endl;
      //read value at each th value
      for( j=0; j<Nt_data; j++)
      {
        fgets(buffer1,LLEN,fp); NL++;
        //printf("j,k = %i,%i %s",j,k,buffer1);
        if( (test=sscanf(buffer1," %lf %lf\n",&tmp[0],&tmp[1]))!=2)
        { printf("%i: Error in reading v and Br.\n",NL);
          cout << buffer1 ; exit(-1);
        }
        index = j + Nt*k;
        data[index].vBr = tmp[0];
        data[index].Br = tmp[1];
  	    data[index].status=1; // make point valid
  	    //Fill last point with value in first to complete full circle
        if( j==0 )
        {
          indexNt = (Nt-1) + Nt*k;
          data[indexNt].vBr = tmp[0];
          data[indexNt].Br = tmp[1];
    	    data[indexNt].status=1; // make point valid
    	    valid_points++;
        }
  	    valid_points++;
      }//j
    }//z
  }

  //Other parameters that are needed from simulation runs.
  // These are not really needed for raw probe data.
  fgets(buffer1,LLEN,fp); NL++;
    sscanf(buffer1,"%lf\n",&R0);
     printf("  R0 = %lf \n", R0);
  fgets(buffer1,LLEN,fp); NL++;
    sscanf(buffer1,"%lf\n",&R_AN);
    printf("  R_AN = %lf \n", R_AN);

  fclose(fp);
  printf(" _______ lines read NL= %i\n",NL);
  //To debug 
  /*k=0;
    for(j=0; j<Nt; j++){
    index = j + Nt*k;
      printf(" %i %1f  %1f  %lf\n", 
        index, data[index].th, data[index].vBr, data[index].Br );
    }*/
  if((-valid_points)>0) { 
    printf("    missing points = %i \n", Nt*Nz-valid_points);
  } else printf("    Read expected number of points, %i\n", valid_points);
  return;
}


void File_Save_format_GEN_V03( char *fname, SimData *data )
{
  //Generalized format (simplified for analysis V03)
  //Write data to file in simplified format
  FILE *fzBr;
  fzBr = fopen(fname,"w");
  fprintf(fzBr,"# ---------- Simulated data ------------- \n");
  fprintf(fzBr,"# model R0= %lf \n", R0);
  fprintf(fzBr,"# --------------------------------------- \n");
  fprintf(fzBr,"# Br(z,th) data\n");
  fprintf(fzBr,"# z followed by pairs of (th,Br) \n");
  fprintf(fzBr," %lf  # r_probes\n", r_probes );
  fprintf(fzBr," %lf  # th_step(deg)\n", th_step/DEGRAD );
  fprintf(fzBr," %lf  # z_step \n", z_step );
  fprintf(fzBr,"  %i  # Nz (number of z points) \n", Nz );
  fprintf(fzBr,"DATA_Z_TH_BR\n");

  for( k=0; k<Nz; k++ ){
  index = 0 + Nt*k ;
  fprintf(fzBr,"  %lf \n", data[index].z); //z position
   for( j=0; j<Nt; j++ ){
     index = j + Nt*k ;
     fprintf(fzBr,"%lf %lf\n", data[index].t/DEGRAD, data[index].Br );
   }
  }
  fprintf(fzBr,"END\n");
  fclose(fzBr);
  return;
}

double Enge( double x, double d, int n, int ee, int extF )
// x = z/R0    d = 2*R0 = full aperture
// n = pole/2 number   ee = 1/2 => entrance/exit
// extF = flag for extended function addition
// assumes x division by aperture has already been done
{
  double value;
  double x_d;
  x_d = x/d;
//  value = coef(n,ee,0) + x_d*( coef(n,ee,1) + x_d*( coef(n,ee,2)
//     + x_d*( coef(n,ee,3) + x_d*( coef(n,ee,4) + x_d*coef(n,ee,5) )))) ;
  value =    fcoef_func(n,ee,0) + x_d*( fcoef_func(n,ee,1) 
     + x_d*( fcoef_func(n,ee,2) + x_d*( fcoef_func(n,ee,3)
     + x_d*( fcoef_func(n,ee,4) + x_d*fcoef_func(n,ee,5) )))) ;
  if( value<-60. ){
    value=-60.;
  } else if( value>30. ){
    value=30.;
  }
  value = 1./(1.+exp(value)) ;
  //check: for( int i=0; i<6; i++) printf(" %lf \n", fcoef_func(n,ee,i) );
  return value;
}

double Enge_product( double x, double d, double LEFF, int n, int extF )
// Product of entrance and exit "Enge" functions about +-LEFF/2
{
  double Enge1, Enge2;
  double x1, x2;
  x1 = (-0.5*LEFF-x)/d ;
  x2 = (x-0.5*LEFF)/d ;
  Enge1 = fcoef_func(n,0,0) + x1*( fcoef_func(n,0,1)
   + x1*( fcoef_func(n,0,2) + x1*( fcoef_func(n,0,3)
   + x1*( fcoef_func(n,0,4) +   x1*fcoef_func(n,0,5) )))) ;
  if( Enge1<-60. ){
    Enge1=-60.;
  } else if( Enge1>30. ){
    Enge1=30.;
  }
  Enge1 = 1./(1.+exp(Enge1)) ;
  //
  Enge2 = fcoef_func(n,1,0) + x2*( fcoef_func(n,1,1)
   + x2*( fcoef_func(n,1,2) + x2*( fcoef_func(n,1,3)
   + x2*( fcoef_func(n,1,4) +   x2*fcoef_func(n,1,5) )))) ;
  if( Enge2<-60. ){
    Enge2=-60.;
  } else if( Enge2>30. ){
    Enge2=30.;
  }
  Enge2 = 1./(1.+exp(Enge2)) ;
  //for( int i=0; i<6; i++) printf(" %lf \n", fcoef_func(n,0,i) );
  return Enge1*Enge2 ;
}

double *Enge_product_array( int NN, double x[], double d, double LEFF, 
  int n, int extF )
// Product of entrance and exit "Enge" functions about +-LEFF/2
{
  double *array;
  int i;
  array = (double*) malloc( NN * sizeof(double) );
  for( i=0; i<NN; i++ ){
    array[i] = Enge_product( x[i], d, LEFF, n, extF );
  }
  free( array );
  return array;
}

void eval_gradients_in_expansion()
{
  // Evaluate higher order b_n,m and a_n,m from m=0 term about z(k)
  // Put results into arrays b_nm and a_nm, index = k + Nz*( m + Nm*n ) 
  double coef;
  for( n=0; n<Nn; n++)
  {
    //Impose Enge function onto m=0 term
    m=0;
    b_nm_1 = Enge_product_array( Nz, z_arr, R0, LEFF, n, 0);
    a_nm_1 = Enge_product_array( Nz, z_arr, R0, LEFF, n, 0);
    //Impose field factor and save to b_nm, a_nm array
    for( k=0; k<Nz; k++ ){
        index = k + Nz*( m + Nm*n ) ;
        //printf(" n,m,k = %i %i %i   index= %i \n", n,m,k, index);
        b_nm[index] = Bn[n]*b_nm_1[k];
        a_nm[index] = An[n]*a_nm_1[k];
    }
    //Evalute higher order terms from derivatives
    for( m=1; m<Nm; m++)
    {
      d_b_nm = derivative_d( Nz, z_arr, b_nm_1);
      d2_b_nm = derivative_d( Nz, z_arr, d_b_nm);  
      d_a_nm = derivative_d( Nz, z_arr, a_nm_1);
      d2_a_nm = derivative_d( Nz, z_arr, d_a_nm);  
      for( k=0; k<Nz; k++ ){
        coef = (double)( 4.*m*(n+m)*(n+2.*m-2.) ) ;
        if(coef>0.)
        {
          coef = (double)(n+2.*m)*R0*R0/coef ;
        } else {
          //Denominator is zero for (n,m)=(0,1); hence, for n=0 there is no 
          // reason to try to evaluate n=0 terms.
          //printf("ERROR:  n%i,m%i  %lf\n",n,m,coef);
          coef = 0. ;
        }
        b_nm_tmp[k] = d2_b_nm[k] * coef ;
        b_nm_1[k] = b_nm_tmp[k];  //for next iteration
        a_nm_tmp[k] = d2_a_nm[k] * coef ;
        a_nm_1[k] = a_nm_tmp[k];
        index = k + Nz*( m + Nm*n ) ;
        //printf(" n,m,k = %i %i %i   index= %i \n", n,m,k, index);
        b_nm[index] = b_nm_1[k];
        a_nm[index] = a_nm_1[k];
        //keep first derivative for evaluation of Bz
        d1_b_nm[index] = d_b_nm[k];
        d1_a_nm[index] = d_a_nm[k];        
      }
    }
  }
  /*
  //Optional:   //Write b_2,0 and b_2,1 to file
  FILE *fp2;
  fp2 = fopen("test_b_2-0_Nn.txt", "w");
  fprintf(fp2,"  z   b_n0 b_n1 b_n2 b_n3 b_n4 b_n5 \n");
  n=2;
  for( k=0; k<Nz; k++ ){
    fprintf(fp2," %lf  %lf %lf %lf %lf %lf %lf \n",
      z_arr[k],
      b_nm[k + Nz*( 0 + Nm*n )],
      b_nm[k + Nz*( 1 + Nm*n )],
      b_nm[k + Nz*( 2 + Nm*n )],
      b_nm[k + Nz*( 3 + Nm*n )],
      b_nm[k + Nz*( 4 + Nm*n )],
      b_nm[k + Nz*( 5 + Nm*n )] );
    //     k + Nz*( m + Nm*n ) ;
  }
  fclose(fp2);
  //
  fp2 = fopen("test_b_6-0_Nn.txt", "w");
  fprintf(fp2,"  z   b_n0 b_n1 b_n2 b_n3 b_n4 b_n5 \n");
  n=6;
  for( k=0; k<Nz; k++ ){
    fprintf(fp2," %lf  %lf %lf %lf %lf %lf %lf \n",
      z_arr[k],
      b_nm[k + Nz*( 0 + Nm*n )],
      b_nm[k + Nz*( 1 + Nm*n )],
      b_nm[k + Nz*( 2 + Nm*n )],
      b_nm[k + Nz*( 3 + Nm*n )],
      b_nm[k + Nz*( 4 + Nm*n )],
      b_nm[k + Nz*( 5 + Nm*n )] );
    //     k + Nz*( m + Nm*n ) ;
  }
  fclose(fp2);
  */
  return;
}

void Eval_Fields_via_Expansion(SimData *data)
{
  //EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
  //  Field evaluation from Enge & other ceofficients vie EXPANSION method
  // NOTE: Abandoned the use of this method until some discrepancies can be 
  //       resolved. Use COSY method instead.
 bool Exp_method=1;
 if( Exp_method )
 {
  eval_gradients_in_expansion();

  double z, r=r_probes, t;
  for( k=0; k<Nz; k++ ){
   z = k*z_step + z_min ; //printf(" k=%i  z= %lf\n",k,z);
   for( j=0; j<Nt; j++ ){
     index = j + Nt*k ;
     data[index].z = z ;
     data[index].t = (double)j * th_step;
     data[index].r = r ;
     //printf(" j=%i  th= %lf\n",j,data[index].th/DEGRAD);
   }
  }
   printf(" Nt=%i, th_step= %lf\n", Nt, th_step/DEGRAD);

  //Evaluate fields at ideal position values
  time0 = time( NULL );
  for( k=0; k<Nz; k++ ){
   for( j=0; j<Nt; j++ ){
     index = j + Nt*k ;
     Field_R0_r_t_z__Br_Bt_Bz(R0,
       data[index].r, data[index].t, data[index].z,
       data[index].Br, data[index].Bt, data[index].Bz);
     /*printf(" %lf %lf %lf %lf %lf %lf\n",
       data[index].r, data[index].t/DEGRAD, data[index].z,
       data[index].Br, data[index].Bt, data[index].Bz);*/
   }
  }
  time1 = time( NULL );
  printf( " Nz*Nt= %i  evaluated in (sec.): %i \n", Nz*Nt, (time1-time0));
  printf( " points/sec = %lf \n", (double)(Nz*Nt)/(double)(time1-time0) );
  // NOTE:  
  // Nn=11, Nm=11, Nz*Nt= 58473  evaluated in (sec.): 18, points/sec = 3248.5
  // Nn=11, Nm=11, points/sec = 11694.6 (3.6 times faster)
  
  //Write data to file in z_Br format
  // Use format expected from Br vs z mapper
  FILE *fzBr;
  fzBr = fopen("FldSim.table","w");
  fprintf(fzBr,"# Simulated data\n");
  fprintf(fzBr,"NO_FILE_FOR_HALL_PROBE_ARRANGEMENT\n");
  fprintf(fzBr," %lf  %lf \n", th_min/DEGRAD, th_step/DEGRAD);
  fprintf(fzBr," %lf  %lf  %lf\n", z_min, z_max, z_step);
  fprintf(fzBr," %lf \n", r_probes);
  double Vm;
  for( k=0; k<Nz; k++ ){
  fprintf(fzBr," 25.1 25.2 25.3 25.4 25.5 25.6 25.7 25.8 \n");
  index = 0 + Nt*k ;
  fprintf(fzBr," %lf \n", data[index].z);
   for( j=0; j<(Nt-1); j++ ){
     index = j + Nt*k ;
     //Vm = 0.2*data[index].Br ; //Assuming simple linear relation
     Vm = (double)j*th_step/DEGRAD;
     fprintf(fzBr," %lf  %lf \n", Vm, data[index].Br);
   }
  }
  fprintf(fzBr,"END\n");
  fclose(fzBr);
  
  //EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE END
 } //Exp_method
}

/*
double rand_gaus_rms_trunc( double rms, double rms_lim )
{
  // Gaussian distribution with rms and truncated at |rms_lim|
  double rnd1, rnd2, val=(rms_lim+1.);
  while( fabs(val)>rms_lim ){
    rnd1 = (double)rand() / double(RAND_MAX);
    rnd2 = (double)rand() / double(RAND_MAX);
    val = rms*sqrt(2.*log(1/rnd1))*cos(6.283185307179586*rnd2);
  }
  return val;
}
*/

void FieldSim_derivatives_V03()
{
  //gROOT->Reset();
  gROOT->ProcessLine(".L filon_MP.c");
  gROOT->ProcessLine(".L plot_and_fit_gradient.C");
  // loads the file above so that its functions can be used
  // loads another file to use its functions
  printf("     FieldSimulationV03: ------------------->\n\n");
  // print line at beginning of code
  char f_dir[512]="/projects/aris/MP/field_params_DEV/param_files_2016/FSQ5/";
  char f_params[128] = "M5_params_FSQ5_9.4Tpm_TEST.txt";
  // creates two names for files or folders to be used lated in code
  //initialize field parameters (e.g. Enge) 
  fcoef = (double*) malloc( NC*NEE*NNC * sizeof(double) );
  z_arr = (double*) malloc( Nz * sizeof(double) );
  Bn = (double*) malloc( Nn * sizeof(double) );
  An = (double*) malloc( Nn * sizeof(double) );
  // allocates memory to arrays that are going to be used in the program
  
  //Initialize z position array using even step size z_step from z_min
  for( i=0; i<Nz; i++ )
  {
    z_arr[i] = z_min + z_step*(double)i ;
  }
  // puts in the z array to be used. the entire array is now made
  //read field description parameters (ref. radius, Enge coefs, etc.)
  //  from file name { f_dir + f_params }
  char fn[640];
  sprintf(fn, "%s%s", f_dir, f_params);
  // makes the file name where the parameters will be read from   
  for( i=0; i<Nn; i++ ) { Bn[i]=0.; An[i]=0.; } //initialize all gradients
  // puts zeros everywhere in the two written arrays. 
  // now everything excpet fcoef is full. 
  read_M5_params_NORM( fn, Bn);
  // fills up Bn array with values from file fn. It also fills up the array fcoef
  // with all of the coefficients from the file. 
  for( i=2; i<7; i++ ) { printf(" Bn[%i] = %lf \n", i, Bn[i]); }

  //Define data with no randomization
  SimData *data0;   // no randomization
  SimData *data_dx;
  SimData *data_dy;
  data0=  (SimData*) malloc( Nt*Nz*sizeof(SimData) );
  data_dx=(SimData*) malloc( Nt*Nz*sizeof(SimData) );
  data_dy=(SimData*) malloc( Nt*Nz*sizeof(SimData) );
  //Initialize arrays for analysis
  //Creates some new data structures of type SimData which is a bunch of doubles and and integer
  Analyze_Field_Initialize();
  
  //Evaluate field distributions using COSY for misaligned magnet
  //Only transverse shift by dx and dy are considered
  double dx=0.01, dy=0.01;

  double zz, rr=r_probes, tt;
  for( k=0; k<Nz; k++ ){
  zz = k*z_step + z_min ; //printf(" k=%i  z= %lf\n",k,z);
  for( j=0; j<Nt; j++ )
  {
     index = j + Nt*k ;
     tt = (double)j * th_step ;
     data0[index].z = zz ; data0[index].t = tt;  data0[index].r = rr ;
     data0[index].xc = 0. ;  data0[index].yc = 0.;
     //printf(" j=%i  th= %lf\n",j,FieldPoint[index].th/DEGRAD);
     data_dx[index].z = zz ; data_dx[index].t = tt;  data_dx[index].r = rr ;
     data_dx[index].xc = dx ;  data_dx[index].yc = 0.;
     //
     data_dy[index].z = zz ; data_dy[index].t = tt;  data_dy[index].r = rr ;
     data_dy[index].xc = 0. ;  data_dy[index].yc = dy;
  }
  }
  
  Field_Simulation_via_COSY_M5( data0 );
  // Optional: Plot field distribution
  Plot_Br_Bt_Bz_VS_z_th( data0 );
  Field_Simulation_via_COSY_M5( data_dx );
  Field_Simulation_via_COSY_M5( data_dy );
  // Optional: Save field data to file 
  File_Save_format_GEN_V03("COSY_0_0.txt", data0);
  //File_Save_format_GEN_V03("COSY_dx_0.txt", data_dx);
  //File_Save_format_GEN_V03("COSY_0_dy.txt", data_dy);

  //Carry out harmonic analysis on each
  double *B_r_n_0, *A_r_n_0; //Nz,Nn for B=normal A=skew
  double *B_r_n_dx, *A_r_n_dx; //Nz,Nn for B=normal A=skew
  double *B_r_n_dy, *A_r_n_dy; //Nz,Nn for B=normal A=skew
  B_r_n_0 = (double*) malloc( Nn*Nz * sizeof(double) );
  A_r_n_0 = (double*) malloc( Nn*Nz * sizeof(double) );
  B_r_n_dx = (double*) malloc( Nn*Nz * sizeof(double) );
  A_r_n_dx = (double*) malloc( Nn*Nz * sizeof(double) );
  B_r_n_dy = (double*) malloc( Nn*Nz * sizeof(double) );
  A_r_n_dy = (double*) malloc( Nn*Nz * sizeof(double) );

  harmonics_B_r_n__A_r_n_SimData( data0, B_r_n_0, A_r_n_0 );
  harmonics_B_r_n__A_r_n_SimData( data_dx, B_r_n_dx, A_r_n_dx );
  harmonics_B_r_n__A_r_n_SimData( data_dy, B_r_n_dy, A_r_n_dy );
  
  //Optional: Evaluate on-axis gradients 
  zFourierTransfAnalysis( B_r_n_0, A_r_n_0);
  write_b_n0_a_n0_Re_and_Im();

  /*
  Optional: Print to file,  harmonic terms for all z
  harmonics_BA_r_n_write(data0, B_r_n_0, A_r_n_0, "BA_r_n_0_0.txt" );
  harmonics_BA_r_n_write(data_dx, B_r_n_dx, A_r_n_dx, "BA_r_n_dx_0.txt" );
  harmonics_BA_r_n_write(data_dy, B_r_n_dy, A_r_n_dy, "BA_r_n_0_dy.txt" );
  */
  
  //Optional: Print to file,  harmonic terms for only one z
  //And plot decomposed harmonics
  //                          z , B_r_n  ,  A_r_n  ,  file name
  harmonics_BA_r_n_Plot_at_z( 0., B_r_n_0,  A_r_n_0,  "z_BA_r_n_0_0.txt" );
  harmonics_BA_r_n_Plot_at_z( 0., B_r_n_dx, A_r_n_dx, "z_BA_r_n_dx_0.txt" );
  harmonics_BA_r_n_Plot_at_z( 0., B_r_n_dy, A_r_n_dy, "z_BA_r_n_0_dy.txt" );

  //Evaluate dx & dy derivatives of B_r_1 & A_r_1 given overall shift in  
  //Magnet alignment by dx,dy  
  derivative_BA_r_1_WRT_dx_dy( "BA_sens_dx_dy.txt", n, dx, dy, 
    B_r_n_0, A_r_n_0, 
    B_r_n_dx, A_r_n_dx,
    B_r_n_dy, A_r_n_dy  );
  
  //Optionals:
  //Move output files to directory called 'output'
  //system("if [ ! -d output ]; then mkdir output; fi"); //create directory
  //system("mv -f *.out COSY*.txt BA*.txt z_BA*.txt ROOT*.png Plot*.png  output/");
  //Delete temporary COSY files that were generated
  system("rm -f RKLOG.DAT TMP_COSY_out.txt *.lis foxyinp.dat TMP_Enge_params_for_COSY.txt");
  //--------------------------------------------
  free( fcoef );
  free( Bn ); free( An );
  free( z_arr );
  free( data0 ); free( data_dx ); free( data_dy );
  free( B_r_n_0 ); free( A_r_n_0 );
  free( B_r_n_dx ); free( A_r_n_dx );
  free( B_r_n_dy ); free( A_r_n_dy );
  //--------------------------------------------
  printf(" ..END.\n");
  return;
}