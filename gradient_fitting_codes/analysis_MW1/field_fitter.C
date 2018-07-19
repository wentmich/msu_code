/* This program executes a full field fit on the magnetic field data of an SQ
using MINUIT and the minimizer MIGRAD. The program minimizes the sum of the
squared residuals as calculated by the functions in the header file 
"field_residual_calculator.h." The initial parameters are read in from the file
"fname_init." Intermediate and final parameters are read into the file
"fname_tmp." The field is fitted to the field described by the data in 
"field_file." The COSY code only uses n = {2,...,6}, so all other parameters are
held fixed in the fit. */

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
# include <TGraph2D.h>
# include "Math/Minimizer.h"
# include "Math/Factory.h"
# include "Math/Functor.h"
# include "field_residual_calculator.h"

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
    NP      = number of total parameters in Enge model function
    OrdCOSY = number of poles analyzed by COSY 
    r_probes= radius of the modeled probes
    NPOLES  = number of multipole components accounted for in fit
    ------- = integers used as indices in loops 
    fname** = initial parameters and intermediate paramters for the fit   */
    
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
double r_probes = 0.1873;
int NPOLES = 11;
int nParams = NP*NPOLES;
int i, j, k, c, ee, n, m, index;

char * fname_model = (char*) "M5_params_FSQ5_9.4Tpm_TEST.txt";
char * fname_init  = (char*) "M5_init_pars.txt";
char * fname_tmp   = (char*) "M5_tmp_fit_file.txt";
char * field_file  = (char*) "COSY_0_0.txt";

double * fcoef = (double*) malloc( NC*NEE*NNC * sizeof(double) );
double * Bn    = (double*) malloc( Nn * sizeof(double) );
double * Bn_tmp = (double*) malloc(sizeof(double)*NPOLES);
double * fcoef_tmp = (double*) malloc( NC*NEE*NNC * sizeof(double) );
/* DECLARATION OF DATA STRUCTURES:
	FldData = structure used to store arrays of field data   */

typedef struct
{
	bool status;
	double r, t, z;
	double Br, Bt, Bz;
	double xc, yc;
} FldData;

/* DECLARATION OF FUNCTIONS:   */

double target_function(double * pars)
{
	TStopwatch timer;
	timer.Start();
	// parameters are read in correctly. They are the initialized variables
	for (int i = 0; i < NPOLES; ++i)
	{
		Bn_tmp[i] = pars[i*NP];
		for (int j = 0; j < NC/2; ++j)
			fcoef_tmp[i*NC*NEE + j] = pars[i*NP + j +1];
		for (int k = NC/2; k < NC; ++k)
			fcoef_tmp[i*NC*NEE + k] = 0;
		for (int m = 0; m < NC/2; ++m)
			fcoef_tmp[i*NC*NEE + m + NC] = pars[i*NP + m + NC/2 +1];
		for (int n = NC/2; n < NC; ++n)
			fcoef_tmp[i*NC*NEE + n + NC] = 0;
	}
	
	write_M5_params_NORM(Bn_tmp, fcoef_tmp, fname_tmp);
	// function definitely works up to here.
	// system("vi M5_tmp_fit_file.txt");
	double chi = field_residual_calculator_from_COSY_field(field_file, fname_tmp);
	printf("chi^2 = %lg\n", chi);
	//gSystem -> Unload("field_residual_calculator.C");
	//field_residual_calculator();
	double clock = timer.RealTime();
	cout << "Time for iteration: " << clock << endl;
	timer.Reset();
	return chi;
}

string to_string(int number)
{
	stringstream ss;
	ss << number;
	string str = ss.str();
	return str;
}

int field_fitter(const char * min_name = "Minuit", const char * fnc_name = "Migrad")
{	
	ROOT::Math::Minimizer * min = 
	ROOT::Math::Factory::CreateMinimizer(min_name, fnc_name);
	
	min -> SetMaxFunctionCalls(100000000);
	min -> SetMaxIterations(10000);
	min -> SetTolerance(0.001);
	min -> SetPrintLevel(1);
	
	ROOT::Math::Functor func(&target_function, nParams);
	min -> SetFunction(func);
	
	double * variables = (double*) malloc(sizeof(double)*nParams);
	double * step      = (double*) malloc(sizeof(double)*nParams);
	double * minimum   = (double*) malloc(sizeof(double)*nParams);
	double * maximum   = (double*) malloc(sizeof(double)*nParams);
	
	//read in initial parameters into the fcoef and Bn array
	read_M5_params_NORM(fname_init);

	for (int i = 0; i < NPOLES; ++i)
	{
		variables[i*NP] = Bn[i];
		for (int j = 0; j < NC/2; ++j)
			variables[i*NP + j + 1] = fcoef[NC*NEE*i + j];
		for (int j = 0; j < NC/2; ++j)
			variables[i*NP + j + 1 + NC/2] = fcoef[NC*NEE*i + j + NC];
	}

	for (int i = 0; i < NP*NPOLES; ++i)
	{
		step[i] = 0.01;
		minimum[i]  = -1.0e2;
		maximum[i]  =  1.0e2;
	}

	for (int i = 0; i < (NP*2); ++i)
		min -> SetFixedVariable(i, to_string(i), variables[i]);
	for (int i = (NP*2); i < (NP*7); ++i)
		min -> SetLimitedVariable(i, to_string(i), variables[i], step[i], minimum[i], maximum[i]);
	for (int i = (NP*7); i < (NP*11); ++i)
		min -> SetFixedVariable(i, to_string(i), variables[i]);
	
	
	//program good until here at least...
	min -> Minimize();
	double fit_min_value = min -> MinValue();
	
	cout << "Minimization Achieved. Minimum Chi^2: " << fit_min_value << endl;
	return 0; 
}