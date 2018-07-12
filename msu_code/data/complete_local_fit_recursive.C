/* The program reads in data from a text file and creates a scatter plot. It 
calculates the effective length of the magnet which the data pertains to. Then
it fits an Enge function to the scattered data and returns the fitted parameters
along with a plot of the fitted function. These fitted parameters are taken for 
each of ten "NF" currents and are then fitted individually (coefficient as a
function of current). The final output is an array of fitted parameters for the 
Enge coefficient polynomials. The function is used for local fits with both
variable and constant effective length. It also outputs the residual sum for the 
2D Enge function which contains the coefficient polynomials and their associated
parameters. 

NOTE: This works for files of the following form...

 90   Amps (Ilim=180) 
 0.2   r0[m] 
 401   NZ data points 
 z[m]  EngeProd 
 -1   0.00223685 
 -0.995   0.00234938 
 -0.99   0.00246708 
 -0.985   0.00259017 
 -0.98   0.00271888 
 -0.975   0.00285345 
 -0.97   0.00299411 
 -0.965   0.00314113 
 -0.96   0.00329477 
 -0.955   0.00345531           */

#include <fstream>
#include <iostream>
#include <string>
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include <iomanip> 
#include <stdio.h>
#include "TCanvas.h"
#include "TFrame.h"
#include "TBenchmark.h"
#include <TGraph.h>
#include "TF1.h"
#include "TH1.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TPaveText.h"
#include <sstream>


using namespace std; 


//double LEFF, R0, IAMPS, ZMAX, ZMIN, DELZ, zz[1], b_zero, integral,
//chi_leff_calc, chi_leff_const, e;
double IAMPS;
const int NP = 13; const int NZ = 401; const int NF = 10; const int NPP = 3;
double ZMAX = 1.; double ZMIN = -ZMAX; double DELZ = 0.005; double LEFF = 0.7;
double R0 = 0.2; double chi_leff_calc = 0;
double chi_leff_const = 0; double e = 2.718281828459045;
double CURRENTS_ARRAY[NF] = {18.,36.,54.,72.,90.,108.,126.,144.,162.,180.};

double * leff_values = (double*) malloc(sizeof(double)*NF);
    
double * length_params = (double*) malloc(sizeof(double)*NPP);

double * all_coefficients = (double*) malloc(sizeof(double)*NF*NP);

double * fitted_parameters = (double*) malloc(sizeof(double)*(NPP+1));
    
double * current_coefficient_set = (double*)malloc(sizeof(double)*NF);
    
double * fitted_parameters_calc_leff = (double*) malloc(sizeof(double)*NPP*NP);
    
double * fitted_parameters_constant_leff = (double*) malloc(sizeof(double)*NPP*NP);
    
double * temp_array = (double*) malloc(sizeof(double)*(NPP));
    
struct ReadData
{
	double * za;
	double * fa;
	ReadData()
	{
		za = (double*) malloc(sizeof(double)*NZ*NF);
		fa = (double*) malloc(sizeof(double)*NZ*NF);
	}
}read_data,all_final_data;


//double all_final_data.za = (double*) malloc(sizeof(double)*NZ*NF);
//double all_final_data.fa = (double*) malloc(sizeof(double)*NZ*NF);
    
/* Defines all the global variables to be used in the program:
LEFF = effective length (calculated for half of program)
R0   = reference radius
ZMAX = range for z position
DELZ = step size in data collection
CURRENTS_ARRAY = list of currents used in data collection
NP   = number of parameters
NZ   = number of data points
NF   = number of currents
NPP  = degree of polynomial + 1
e    = defines the constant e to use in place of exp function
read_data.za   = temporary array for z position points
read_data.fa   = temporary array for field points 
b_zero = changing variable for b_field at z = 0. used in calc_LEFF()
integral = changing variable for value of integral used in calc_LEFF()
leff_values = array for effective lenths (NF elements) 
length_params = array to temporarily hold values of fitted LEFF parameters
all_coefficients = array for all polynomial coefficients
fitted_parameters = temporary array for fitted polynomial coefficients
current_coefficient_set = temporary array for fitting polynomials
fitted_parameters****** = final p_n parameters 
temp_array = array to hold polynomial coefficients and chi_squared
all_final_data = array to hold every (z,b) point for all currents */

double func_second_poly(double *I, double *params)
{
// defines and returns function with three paramters (a second order polynomial)
    double func = (params[0] + params[1]*I[0] + params[2]*pow(I[0],2));
    return func;
}

double * read_initial_parameters(const char * FileName = "")
{
/* this function reads in data from a file to create the initial parameter array
It outputs an array of "NP" parameters which can then be used to set the initial
conditions on the fit. */
	double * initial_parameters;
	initial_parameters = (double*) malloc(sizeof(double)*NP);
	const int Nline = 256;
	char string_location[Nline];
	FILE * file;
	file = fopen(FileName, "r");
	for (int i = 0; i < NP; ++i)
	{
		fgets(string_location, Nline, file);
		sscanf(string_location, "%lg \n", &initial_parameters[i]);
	}
	return initial_parameters;
}

void write_final_parameters(const char * FileName, double * final_parameters)
{
	ofstream file;
	file.open(FileName);
	for (int i = 0; i < NP; ++i)
		file << final_parameters[i] << '\n';
	file.close();
}

double * scatter_data_from_array(double * a, int j)
{
/* Plots data as a scatter plot (CURRENTS_ARRAY = x_axis, a = y_axis). The input
array should have NF elements. It then fits a second degree polynomial to the 
data and returns an array of the fitted parameters with the residual sum. */
    for (int i = 0; i < NF; ++i)
        current_coefficient_set[i] = a[j + i*NP];
    
    TGraph * data = new TGraph(NF,CURRENTS_ARRAY,current_coefficient_set); 
    // fits data with a second degree polynomial as defined above   
    TF1 * func = new TF1("func", func_second_poly, 0, 200, NPP);
    for ( int i = 0; i < NPP; ++i)
        func -> SetParameter(i, 1.);
    data -> Fit("func");
    for (int i = 0; i < NPP; ++i)
        fitted_parameters[i] = func -> GetParameter(i);
    return fitted_parameters;
}

double func_EngeProd(double *z, double *p)
{
/* creates a function to fit to the scattered data. It is a one dimensional 
Enge function (a function of z position only). In the case that LEFF is 
calculated, a value for it can be substituted into the function calling this
function creator. It also only has 13 parameters for the c_i coefficients. */
  double zi =  z[0] - LEFF/2.0;             // defines z variable
  double ze = -z[0] - LEFF/2.0;             // defines z variable
  zi /= 2*R0; ze /= 2*R0;
  double fi=1./(1.+exp(p[0]+zi*(p[1]+zi*(p[2]+zi*(p[3]+zi*(p[4]+zi*p[5]))))));
  double fe=1./(1.+exp(p[6]+ze*(p[7]+ze*(p[8]+ze*(p[9]+ze*(p[10]+ze*p[11]))))));
  return p[12]*fi*fe;
}

double gen_Enge_func(double * z, double * p)
{
/* Creates the function to be fitted to the input data. Here Leff is a funciton
of current. It can be held constant by setting "p[39]" to LEFF and "p[40]/p[41]"
to zero. The function returns a double object with "NP" parameters. */
	double zi =  z[0] - (p[39]+p[40]*p[42]+p[41]*p[42]*p[42])/2.0; // defines the z variable
	double ze = -z[0] - (p[39]+p[40]*p[42]+p[41]*p[42]*p[42])/2.0; // defines the z variable
	zi /= 2*R0; ze /= 2*R0;
	double fi=1./(1.+pow(e,((p[0]+p[1]*p[42]+p[2]*p[42]*p[42]) + zi*((p[3]+p[4]*p[42]+p[5]*p[42]*p[42])
		+ zi*((p[6]+p[7]*p[42]+p[8]*p[42]*p[42]) + zi*((p[9]+p[10]*p[42]+p[11]*p[42]*p[42])
		+ zi*((p[12]+p[13]*p[42]+p[14]*p[42]*p[42]) + zi*(p[15]+p[16]*p[42]+p[17]*p[42]*p[42]))))))));
	double fe=1./(1.+pow(e,((p[18]+p[19]*p[42]+p[20]*p[42]*p[42]) + ze*((p[21]+p[22]*
		p[42]+p[23]*p[42]*p[42]) + ze*((p[24]+p[25]*p[42]+p[26]*p[42]*p[42]) + ze*((p[27]+p[28]*p[42]+p[29]
		*p[42]*p[42])+ze*((p[30]+p[31]*p[42]+p[32]*p[42]*p[42])+ze*(p[33]+p[34]*p[42]+p[35]*p[42]*p[42]))))))));
	double B = p[36] + p[37]*p[42] + p[38]*p[42]*p[42];
	double FinalFunction = B*fi*fe;
	return FinalFunction;
}

double get_chi_squared(double * parameters, double * length_params, double * z, double * b)
{
/* this function takes the array parameters as an argument. It creates a
function from "gen_Enge_func" and sets all the parameters to the proper settings
from the array. It then evaluates the total residual for each of the currents
and sums them together. */
	double chi_squared = 0;
	TF1 * fitted_function = new TF1("fitted_function",gen_Enge_func,ZMIN,ZMAX,(NPP*NP + NPP + 1));
	for (int i = 0; i < (NP*NPP); ++i)
		fitted_function -> SetParameter(i,parameters[i]);
	for (int i = (NP*NPP); i < (NP*NPP + NPP); ++i)
		fitted_function -> SetParameter(i,length_params[i - (NP*NPP)]);
/* has now created a function of the same form as the complete fitted function 
a current parameter to be specified. */
	for (int j = 0; j < NF; ++j)
	{
		fitted_function -> SetParameter(NP*NPP + NPP, CURRENTS_ARRAY[j]);
		for (int i = 0; i < NZ; ++i)
		{
			double func_value = fitted_function -> Eval(z[/*j*NZ + */i]);
			double residual   = pow(func_value - b[j*NZ + i],2);
			chi_squared = chi_squared + residual;
		}
	}
	return chi_squared;
}
	

char * create_name_from_index(double index = 0.10)
{
/* Function creates a files name using the index number given as an argument. It
returns a constant character array pointer by creating a string and converting
it to a pointer. */
    int k = 2; // sets the number of decimal precision to be read into string
    stringstream stream;
    stream << fixed << setprecision(k) << index; // makes "index" a string
    string mystring = stream.str();  
    
    string front ("EngeProd_");
    string back  (".txt");
    string FileName = front + mystring + back;  // creates file name string
    
    const int n = FileName.length();
    char * fname;
    fname = (char*) malloc((n+1) * sizeof(char)); // space for string + null
    strcpy(fname, FileName.c_str());
    return fname;
}

double calc_LEFF(double * z_data, double * b_field)
{
/* function to determine the effective length using the definition 
LEFF = integral(b_field)/(b_field at z=0). The input arrays should each have NZ 
elements. */
	double b_zero;
	double integral;
    int z_zero_index = 0;
    for (int i = 0; i < NZ; i++){
    // finds index of z = 0
        if (fabs(z_data[i]) < 1.0e-10)
        z_zero_index = i;}
    TGraph * data = new TGraph(NZ, z_data, b_field); // plots arrays
    b_zero = b_field[z_zero_index];
    integral = data -> Integral();
    LEFF = integral/b_zero;
    return LEFF;
}

ReadData read_single_file_make_ReadData(const char * f_name)
{
/* function reads in data from file "fname" and creates two arrays inside of a 
structure ReadData. It then returns these two arrays. */
	ReadData temp_data_holder;
	const int NLINE=256;
    FILE *fp;                   //file pointer
    char line[NLINE];           //buffer to read lines
    printf(" points per excitation, NZ= %i \n", NZ);
    temp_data_holder.za=(double*)malloc(NZ*sizeof(double));//creates z array
    temp_data_holder.fa=(double*)malloc(NZ*sizeof(double));//creates field array

    //Read field gradient data from file. Generate scatter plot.
    fp = fopen(f_name, "r");
    fgets(line, NLINE, fp);   sscanf(line, " %lg ", &IAMPS);
    fgets(line, NLINE, fp);   sscanf(line, " %lg ", &R0);
    fgets(line, NLINE, fp);   sscanf(line, " %i ", &NZ);
    fgets(line, NLINE, fp);
    printf(" IAMPS= %lg, R0= %lg, NZ= %i \n", IAMPS, R0, NZ);
    for( int k=0; k<NZ; k++)
    {
      fgets(line, NLINE, fp);                       // reads data to z,f arrays
      sscanf(line," %lg %lg\n", &temp_data_holder.za[k], &temp_data_holder.fa[k]);
    }
    fclose(fp);
    
    return temp_data_holder;
}

void read_single_file_and_fit(double index = 0.10,int j=1, int i=1, double * p)
{
/* function reads in a file of the form shown above. It then scatters the data
and fits a function to the data. The file it reads in is "EngeProd_*.**.txt".
For an explanation of the fitting process and the TF1 class, see 
5.3.3 of the ROOTUsersGuideLetter2014.pdf. See also 
https://root.cern.ch/root/HowtoFit.html*/
	
    char * fname = create_name_from_index(index); // creates file name
    read_data = read_single_file_make_ReadData(fname);
    for (int k = 0; k < NZ; ++k)
    {
    	all_final_data.za[k + j*NZ] = read_data.za[k];
    	all_final_data.fa[k + j*NZ] = read_data.fa[k];
    }
    // creates scatter plot of z vs field
    TGraph * data = new TGraph(NZ, read_data.za, read_data.fa);
    // find effective length using integral and b(z=0)
    if (i == 0){
    	LEFF = calc_LEFF(read_data.za,read_data.fa);
    	leff_values[j] = LEFF;}
    else
    	LEFF = 0.7;
    // Define the fuction to fit
    TF1 *func = new TF1( "func", func_EngeProd, ZMIN,ZMAX, NP);
    for( int i=0; i<NP; i++)
        func->SetParameter(i, p[i]); // initialize the parameters
     
    // Optional: Fit the parameters to data 
    bool fit_TF1=1;  //Set to 1 to allow fitting
    if(fit_TF1)
        data->Fit("func");
    
    // double coefficients[NP]; // creates an array of the fitted parameters
    for (int i = 0; i < NP; ++i)
        all_coefficients[i+(j*NP)] = func -> GetParameter(i);
}

void wipe_file(char * wiped_file)
{
    fstream temp_file;
    temp_file.open(wiped_file,ios::out | ios::trunc);
    temp_file.close();
}

int main()
{
/* main function reads, scatters, and fits the data for all the current files in
the list "list_of_files". This list gives the numbers to be input into the file
name EngeProd_*.**.txt. In the "read_single_file_and_fit" function, set the 
third parameters to zero to fit with a calculated effective length and set it to 
one to fit to fit with a constant effective legnth. */
	//double all_final_data.za = (double*) malloc(sizeof(double)*NZ*NF);
	//double all_final_data.fa = (double*) malloc(sizeof(double)*NZ*NF);
    double list_of_files[]={0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00};
    double chi_calc = 1.; double chi_const = 1.; double chi_acceptable = 0.1;
    int iterations_calc = 0; int interations_const = 0;
    double * fin_params = (double*) malloc(sizeof(double)*NP);

	while (chi_calc > chi_acceptable)
	{
		double * pars = read_initial_parameters("initial_parameters_local.txt");
		for (int i = 0; i < NP; ++i)
			cout << pars[i] << endl;
		for (int i = 0; i < 10; ++i)
			read_single_file_and_fit(list_of_files[i],i,0,pars);
		for (int j = 0; j < NP; ++j){
			temp_array = scatter_data_from_array(all_coefficients, j);
			for (int i = 0; i < NPP; ++i)
				fitted_parameters_calc_leff[i+(j*NPP)] = temp_array[i];}
		TF1 * leff_poly = new TF1("leff_poly",func_second_poly,0,180,NPP);
		for (int i = 0; i < NPP; ++i)
			leff_poly -> SetParameter(i,1.);
		TGraph * leff_vs_curr = new TGraph(NF, CURRENTS_ARRAY, leff_values);
		
		leff_vs_curr -> Fit("leff_poly");
		for (int i = 0; i < NPP; ++i)
			length_params[i] = leff_poly -> GetParameter(i);
		chi_calc = get_chi_squared(fitted_parameters_calc_leff, length_params, all_final_data.za, all_final_data.fa);
		
		for (int i = 0; i < NP; ++i)
			fin_params[i] = fitted_parameters_calc_leff[NPP*i];
	
		wipe_file("initial_parameters_local.txt");
		write_final_parameters("initial_parameters_local.txt",fin_params);
		++iterations_calc;
		for (int i = 0; i < NP; ++i)
			cout << pars[i] << endl;
		cout << chi_calc << endl;
	}
/*
    for (int i = 0; i < 10; ++i)
        read_single_file_and_fit(list_of_files[i],i,1,pars);
    for (int j = 0; j < NP; ++j){
        temp_array = scatter_data_from_array(all_coefficients, j);
        for (int i = 0; i < NPP; ++i)
            fitted_parameters_constant_leff[i+(j*NPP)] = temp_array[i];}
    length_params[0] = 0.7; length_params[1] = 0; length_params[2] = 0;
    chi_leff_const = get_chi_squared(fitted_parameters_constant_leff, length_params, all_final_data.za, all_final_data.fa);
  */  
    cout << "Calculated LEFF" << '\t' << "Fixed LEFF" << endl;
    for (int i = 0; i < (NP*NPP); ++i)
        cout << fitted_parameters_calc_leff[i] << '\t' << endl;
        //<< fitted_parameters_constant_leff[i] << endl;
    cout << "chi leff constant = " << chi_const << endl;
    cout << "chi leff calc     = " << chi_calc  << endl;
    cout << "interations for calc fit = " << iterations_calc << endl;
    return 0;
}
