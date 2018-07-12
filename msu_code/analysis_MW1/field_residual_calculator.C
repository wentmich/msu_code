/* 
- read in two data structures from files containing the field distributions
- calculate the residuals between the fields in each data structure
- return a chisquared value between the two field distributions
*/

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>
# include <sstream>
# include <TCanvas.h>
# include <TF1.h>
# include <TGraph.h>
# include <fstream>

/* DECLARATION OF GLOBAL VARIBALES:
    sep     = delimiters used by "strtok" funciton
    PI      = constant
    DEGRAD  = degrees to radians conversion
    LEFF    = effective length of magnet
    R0      = reference radius of magnet measurements
    e       = constant
    z_min   = minimum z value in measurements
    z_max   = maximum ----------""-----------
    z_step  = step size for z measurements
    th_min  = minimum theta value in radians
    th_max  = maximum ----------""-----------
    th_step = step size for theta measurements
    fcoef   = pointer to array of Enge coefficients
    Bn      = pointer to array of maximum normal B components for multipoles
    An      = --------------""----------- skew   -----------""--------------
    Nz      = number of distinct z data points
    Nk      = number of distinct z data points
    Nt      = number of distinct theta data points
    Nn      = number of multipole components
    NEE     = number of Enge functions in model
    Nm      = number of multipole components
    NC      = number of coefficients
    index   = index used in read_M5_params_NORM()  
    NP      = number of total parameters in Enge model function  */
    
const char sep[2]=" ,";
double PI=3.141592654;
double DEGRAD=3.141592654/180.;
double LEFF = 0.782 ;
double R0 = 0.2;
double e  = 2.718281828459045;
double z_min = -1.0, z_max = 1.0, z_step = 0.005;
double th_min = 0., th_max = 2.*PI, th_step = 5.*DEGRAD;
double *fcoef;
double *Bn;
double *An;
const int Nz = 1+(z_max-z_min)/z_step ;
int Nk = Nz;
const int Nt = (int)( (th_max-th_min)/th_step)+1;
const int Nn = 11;
const int NNC = 11;
const int NEE = 2;
const int Nm = 11;
const int NC = 12;
int NP = NC + 1;
int index;
int OrdCOSY = 6;

/* DECLARATION OF DATA STRUCTURES:
	FldData = structure used to store arrays of field data   */

typedef struct
{
	bool status; // flag, if==1 data point is valid, ==0 invalid
	double r, t, z; // t(radians)  r,z(usually meters)
	double Br, Bt, Bz; // B(usually Tesla)
	double xc, yc; //shift of magnet center (-xc,-yc)
} FldData;

/* DECLARATION OF FUNCTIONS: */

FldData * File_Read_zScan_FldData( char *fname)
{
/* Function reads in the data from the file "fname" and places it into a
"FldData" type data structure called "data" which it returns. */
	FldData *data
	data=  (FldData*) malloc(Nt*Nz*sizeof(FldData));
	printf("\tFile_Read_zScan_FldData: %s\n", fname);
	char format_flag[] = "DATA_Z_TH_BR_FldData";
	FILE * field_file;
	char buffer1[LLEN];
	field_file = fopen(fname,"r");
	
	if (field_file==NULL)
	{
		printf("no file, %s. Aborting.\n", fname);
		exit(-1);
	}
	
	sprintf(buffer1,"#");
	
	while(strncmp(buffer1,"#",1)==0)
	{
		fgets(buffer1,LLEN,fin);
		printf(buffer1);
	}
	
	sscanf(buffer1," %lf \n", &r_probes );
	fgets(buffer1,LLEN,fin);
	sscanf(buffer1," %lf \n", &th_step );
	th_step *= DEGRAD;
	fgets(buffer1,LLEN,fin);
	sscanf(buffer1," %lf \n", &z_step );
	fgets(buffer1,LLEN,fin);
	sscanf(buffer1," %i \n", &Nz );
	fgets(buffer1,LLEN,fin);
	sscanf(buffer1," %i \n", &Nt );
	fgets(buffer1,LLEN,fin);
	printf(buffer1);
	if( strncmp(buffer1,format_flag,strlen(format_flag) )!=0 )
	{
		printf(" Unexpected header: %s\n Continue?(y/n)", buffer1);
		cin >> buffer1; if( strncmp(buffer1,"y",1)!=0 ) exit(0);
	}
	int kk, ntest;
	for( k=0; k<Nz; k++ )
	{
		fgets(buffer1,LLEN,fin);
		sscanf(buffer1," %i \n", &kk );
		for( j=0; j<Nt; j++ )
		{
			index = j + Nt*k ;
			fgets(buffer1,LLEN,fin); //cout<<buffer1;
			ntest = sscanf(buffer1," %lf %lf %lf\n", &data[index].z, &data[index].t, &data[index].Br );
			data[index].t *= DEGRAD; //convert degrees to radians
			if( ntest!=3 )
			{
				printf(" Read unexpected input: %s\n  Aborting.", buffer1);
				exit(0);
			}
		}
	}
	fgets(buffer1,LLEN,fin); //cout<<buffer1;
	if( strncmp(buffer1,"END",3)==0 ) printf(" ...END. Read normally.\n");
	fclose(fin);
	return data;
}

void read_M5_params_NORM( char *fname)
{
/* reads in the Enge coefficients and field parameters from the file denoted by
"fname." It then returns two arrays, one for the field parameters (Bn_) and one
for the Enge coefficients (fcoef). */
	double * Bn_;
	double * fcoef;
	fcoef = (double*) malloc( NC*NEE*NNC * sizeof(double) );
	Bn = (double*) malloc( Nn * sizeof(double) );
	FILE *fp; char line[512]; char *pch;
	const char sep[2]=" ,";
	fp = fopen(fname, "r");
	sprintf(line,"#");
	while(strncmp(line,"#",1)==0)
		fgets(line,512,fp);
	sscanf(line," %lf\n", &R0);
	fgets(line,512,fp);
	sscanf(line," %lf\n", &LEFF);
	sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
	int i=2;
	pch = strtok(line, sep); sscanf(pch,"%lf", &Bn_[i]);
	for( i=3; i<7 ; i++ )
	{
		pch = strtok(NULL, sep);
		sscanf(pch,"%lf", &Bn_[i]);
	} 
	sprintf(line,"#");
	while(strncmp(line,"#",1)==0)
		fgets(line,512,fp);
	while( !feof(fp) )
	{
		pch = strtok(line, sep);
		sscanf(pch,"%i", &n);
		pch = strtok(NULL, sep);
		sscanf(pch,"%i", &ee);
		for( c=0; c<NC; c++)
		{
			index = c + NC*( NEE*n + ee );
			pch = strtok(NULL, sep);
			sscanf(pch,"%lf", &fcoef[index]);
		}
		fgets(line,512,fp);
	}
	pch = strtok(line, sep);
	sscanf(pch,"%i", &n);
	pch = strtok(NULL, sep);
	sscanf(pch,"%i", &ee);
	for( c=0; c<NC; c++)
	{
		index = c + NC*( NEE*n + ee );
		pch = strtok(NULL, sep);
		sscanf(pch,"%lf", &fcoef[index]);
	}
	return;
}

void Field_Simulation_via_COSY_M5_FldData(FldData *data, double * fcoef, double * Bn)
{
/* takes the r, t, z elements in "data" and uses them to create x, y, z points
for the simulation. It creates a file for Enge coefficients and field paramters
to be read in by the COSY code. The field simulation is done, and the final
field result is read into "data." */
	char FS_fname[50]="TMP_COSY_out.txt" ;
	double xx, yy, zz, Bxx, Byy, Bzz, Brr, Btt, tt;
	bool call_COSY=1;
	if( call_COSY )
	{
		FILE *fpC;
		fpC = fopen("TMP_Enge_params_for_COSY.txt", "w");
		fprintf(fpC, " %d\n", OrdCOSY);
		fprintf(fpC, " %lf\n", LEFF);
		fprintf(fpC, " %lf\n", R0);
		for( n=1 ; n<Nn; n++ )
		{
			for( ee=0; ee<NEE; ee++ )
			{
				fprintf(fpC, " %i\n", n);
				fprintf(fpC, " %i\n", ee);
				for( c=0; c<NC; c++ )
				{
					index = c + NC*( NEE*n + ee );
					fprintf(fpC, " %lf\n", fcoef[index]);
				}
			}
		}
		fprintf(fpC,"Normal field components\n");
		for( i=0; i<Nn ; i++ )
			fprintf(fpC," %lf \n", Bn[i]);
		fprintf(fpC,"Skew field components\n");
		for( i=0; i<Nn ; i++ )
			fprintf(fpC," %lf \n", An[i]);
		fprintf(fpC,"Points for field evaluation\n");
		fprintf(fpC,"%s\n",FS_fname);
		for( k=0; k<Nz; k++ )
		{
			for( j=0; j<Nt; j++ )
			{
				index = j + Nt*k ;
				xx = -data[index].xc + data[index].r * cos( data[index].t );
				yy = -data[index].yc + data[index].r * sin( data[index].t );
				fprintf(fpC, " %lf,%lf,%lf\n", xx, yy, data[index].z);
			}
		}
		fprintf(fpC,"END\n");
		fclose(fpC);
		time0 = time( NULL );
		system("runcosy91.sh COSY_field_Gen.fox");
		time1 = time( NULL );
		printf( " Nz*Nt= %i  evaluated in (sec.): %i \n", Nz*Nt, (time1-time0));
		printf( " points/sec = %lf \n", (double)(Nz*Nt)/(double)(time1-time0) );
		printf( " sec/point = %lf \n", (double)(time1-time0)/(double)(Nz*Nt) );
	}
	FILE *fpCO;
	char line[512]; 
	int Nscan;
	fpCO = fopen(FS_fname, "r");
	if (fpCO==NULL)
	{
		printf("File not found: %s\n",FS_fname);
		exit(0);
	}  
	for( k=0; k<Nz; k++ )
	{
		for( j=0; j<Nt; j++ )
		{
			index = j + Nt*k ;
			fgets( line, 512, fpCO );
			Nscan = sscanf( line, " %lf %lf %lf %lf %lf %lf\n", &xx, &yy, &zz, &Bxx,&Byy,&Bzz);
			if( Nscan!=6 )
			{
				printf(" Error reading point: k%i j%i \n",k,j);
				exit(0);
			}
			data[index].Br = Bxx*cos(data[index].t) + Byy*sin(data[index].t);
			data[index].Bt =-Bxx*sin(data[index].t) + Byy*cos(data[index].t);
			data[index].Bz = Bzz;
		}
	}
	fclose( fpCO );
	return;
}

double get_chi_squared(FldData * exp_data, FldData * model_data)
{
/* finds the sum of the squared residuals between two data sets and returns that
value as "chi_squared." */
	double add_on = 0;
	double exp_b = 0;
	double model_b = 0;
	double chi_squared = 0;
	for (int i = 0; i < (Nt*Nz); ++i)
	{
		exp_b = exp_data[i].Br;
		model_b = model_data[i].Br;
		add_on = pow(2,(exp_b - model_b));
		ch_squared = chi_squared + add_on;
	}
	return chi_squared;
}

double get_chi_squared_from_function(FldData * exp_data, double * fcoef, double * Bn)
{
	