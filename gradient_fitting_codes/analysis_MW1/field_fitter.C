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
int i, j, k, c, ee, n, m, index;

char * fname_model = (char*) "M5_params_FSQ5_9.4Tpm_TEST.txt";
char * fname_fit   = (char*) "M5_sim_fit.txt";
char * fname_init  = (char*) "M5_params_FSQ5_9.4Tpm_TEST.txt";

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

double target_function(double * parameters)
{
	gROOT -> ProcessLine(".L field_residual_calculator.C");
	double * Bn_tmp = (double*) malloc(sizeof(double)*NPOLES);
	double * fcoef_tmp = (double*) malloc(sizeof(double)*NPOLES*NC);
	for (int i = 0; i < NPOLES; ++i)
	{
		Bn_tmp[i] = parameters[i*NP];
		for (int j = 0; j < NC; ++j)
			fcoef_tmp[j*i] = parameters[i*NP + j + 1];
	}
	write_M5_params_NORM(Bn_tmp, fcoef_tmp);
	double chi_squared = field_residual_calculator();
	return chi_squared;
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
	gROOT -> ProcessLine(".L field_residual_calculator.C");
	
	double * fcoef = (double*) malloc( NC*NEE*NNC * sizeof(double) );
	double * Bn    = (double*) malloc( Nn * sizeof(double) );
	
	ROOT::Math::Minimizer * min = 
	ROOT::Math::Factory::CreateMinimizer(min_name, fnc_name);
	
	min -> SetMaxFunctionCalls(100000000);
	min -> SetMaxIterations(3);
	min -> SetTolerance(0.1);
	min -> SetPrintLevel(1);
	
	int nParams = NP*NPOLES;
	
	ROOT::Math::Functor func(&target_function, nParams);
	min -> SetFunction(func);
	
	double * variables = (double*) malloc(sizeof(double)*nParams);
	double * step      = (double*) malloc(sizeof(double)*nParams);
	double * minimum   = (double*) malloc(sizeof(double)*nParams);
	double * maximum   = (double*) malloc(sizeof(double)*nParams);
	
	read_M5_params_NORM(fname_init);
	
	for (int i = 0; i < NPOLES; ++i)
	{
		variables[i*NP] = Bn[i];
		cout << "!!!!!!!!!!!!!!!!!!!!" << endl;
		for (int j = 0; j < NC; ++j)
			variables[i*NP + j + 1] = fcoef[j*i];
	}

	for (int i = 0; i < NP*NPOLES; ++i)
	{
		step[i] = 0.00001;
		minimum[i]  = -1.0e3;
		maximum[i]  =  1.0e3;
	}

	for (int i = 0; i < NP*NPOLES; ++i)
		min -> SetLimitedVariable(i, to_string(i), variables[i], step[i], minimum[i], maximum[i]);
	
	min -> Minimize();
	fit_min_value = min -> MinValue();
	
	cout << "Minimization Achieved. Minimum Value: " << fit_min_value << endl;
	return 0; 
}