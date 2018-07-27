/* This program starts by reading in the data from "NF" files containing the 
z_position and the magnetic field for a series of points within the magnet
at a fixed radius "r0." The data from the "NF" files is then placed into a 
single large structure of three arrays, InputData. These arrays are plotted into
a three dimensional scatter plot. The function fits a 2D Enge function to the 
scatter plot. The Enge function is a function of z_position and current. The 
effective length is calculated here. It reads the initial parameters from the
file "initial_parameters_global.txt". 
NOTE: the file reads in data from files of the following form...

 90   Amps (Ilim=180) 
 0.2   r0[m] 
 401   Nz data points 
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
 -0.955   0.00345531 
 
As shown, the first column is the z_position, and the second column is the 
magnetic field. The program will run recursively for a certain number of
iterations and then read final parameters to "final_parameters_global.txt" */

#include <iostream>
#include <string>
#include "TF1.h"
#include <TF2.h>
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
#include <TGraph2D.h>
#include <fstream>
#include "TStopwatch.h"
#include "TObject.h"
using namespace std; 

double LEFF = 0.7; double R0 = 0.2; double ZMIN = -1.; double ZMAX = -ZMIN;
double IMAX = 200.; double IMIN = 0.; 
double e = 2.718281828459045;

const int NZ = 401; const int NF = 10; 
int TOTAL_DATA_POINTS = NZ * NF;
const int I_LIM = 180; const int NP = 42;

double currents_array[NF] = {18.,36.,54.,72.,90.,108.,126.,144.,162.,180.};

// Declaration of structure
struct input_datas
{
	double * Iamps;
	double * z_pos;
	double * b_field;
	input_datas()
	{
		Iamps = (double*) malloc(sizeof(double)*NZ);
		z_pos = (double*) malloc(sizeof(double)*NZ);
		b_field = (double*) malloc(sizeof(double)*NZ);
	}
};

double max_iterations = 30;

input_datas read_file_make_input_datas(const char * FileName)
{
/* The function reads in data from a file in the form shown above and creates a
structure "new_data" in the form of "input_datas." Each column of the file is
read to an array "z_pos" or "b_field," and "Iamps" is filled with the same value
of the current. */
	const int Nline = 256;               // maximum characters read at a time
	char string_location[Nline];         // place to store characters
	double current[1];                   // value of the current for the file
	input_datas new_data;                 // structure to be returned
	FILE * file;
	file = fopen(FileName, "r");         // opens the file
	new_data.Iamps   = (double*) malloc(NZ * sizeof(double));
	new_data.z_pos   = (double*) malloc(NZ * sizeof(double));
	new_data.b_field = (double*) malloc(NZ * sizeof(double));
	fgets(string_location, Nline, file); 
	sscanf(string_location, "%lg", &current[0]);
	fgets(string_location, Nline, file); // opens and the lets go of lines 2-4
	fgets(string_location, Nline, file);
	fgets(string_location, Nline, file);
	for ( int i = 0; i < NZ; ++i)
	{
		// opens each line and reads in the two values for z_pos and b_field
		fgets(string_location, Nline, file);
		sscanf(string_location, "%lg %lg\n",
			&new_data.z_pos[i], &new_data.b_field[i]);
		new_data.Iamps[i] = current[0];//sets each current entry to file current
	}
	return new_data;
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

double * combine_arrays(double * array1, double * array2, int k)
{
/* Takes two arrays of the length specified by "NZ" and replaces 
them into a single array. Returns a pointer to the combined array." */
	double * NewArray;
	NewArray = (double*) malloc(sizeof(double)*(NZ*k) 
		+ sizeof(double)*NZ);
	int Nelements1 = NZ*k;
	int Nelements2 = NZ;
	for (int i = 0; i < Nelements1; ++i)
		NewArray[i] = array1[i];
	for (int i = Nelements1; i < (Nelements1 + Nelements2); ++i)
		NewArray[i] = array2[i - Nelements1];
	return NewArray;
}

TGraph2D * scatter(double * array1, double * array2, double * array3)
{
/* Creates a TGraph2D object of a scatter plot of the data of three arrays. The 
order of the array arguments determines the axes of the plot: array1 = x, array2
= y, array3 = z. */
	if (sizeof(array1) != sizeof(array2))
		cout << "Function not compatible; arrays have different size." << endl;
	else 
	{
		if (sizeof(array2) != sizeof(array3))
			cout << "Function not compatible; arrays have different size." << endl;
		else
		{
			int Nelements = TOTAL_DATA_POINTS;
			TGraph2D * ScatterPlot = new TGraph2D(Nelements, array1, array2, array3);
			return ScatterPlot;
		}
	}
}

double func_EngeProd_calc_leff(double * z, double * p)
{
/* Creates the function to be fitted to the input data. Here Leff is held
constant. The function returns a double object with "NP" parameters. */
	double I = z[1];			   // defines the current variable
	double zi =  z[0] - (p[39]+p[40]*I+p[41]*I*I)/2.0;  // defines the z variable
	double ze = -z[0] - (p[39]+p[40]*I+p[41]*I*I)/2.0;  // defines the z variable
	zi /= 2*R0; ze /= 2*R0;
	double fi=1./(1.+pow(e,((p[0]+p[1]*I+p[2]*I*I) + zi*((p[3]+p[4]*I+p[5]*I*I)
		+ zi*((p[6]+p[7]*I+p[8]*I*I) + zi*((p[9]+p[10]*I+p[11]*I*I)
		+ zi*((p[12]+p[13]*I+p[14]*I*I) + zi*(p[15]+p[16]*I+p[17]*I*I))))))));
	double fe=1./(1.+pow(e,((p[18]+p[19]*I+p[20]*I*I) + ze*((p[21]+p[22]*
		I+p[23]*I*I) + ze*((p[24]+p[25]*I+p[26]*I*I) + ze*((p[27]+p[28]*I+p[29]
		*I*I)+ze*((p[30]+p[31]*I+p[32]*I*I)+ze*(p[33]+p[34]*I+p[35]*I*I))))))));
	double B = p[36] + p[37]*I + p[38]*I*I;
	double FinalFunction = B*fi*fe;
	return FinalFunction;
}

double calc_Leff(double * z_data, double * b_field, TGraph * data)
{
/* function to determine the effective length using the definition 
Leff = integral(b_field)/(b_field at z=0). */
	int z_zero_index = -5;
    for (int i = 0; i < NZ; i++)
    {
    	if (fabs(z_data[i]) < 1.0e-10)
        z_zero_index = i;
    }
    double Bzero = b_field[z_zero_index];
    double integral = data -> Integral();
    double Leff = integral/Bzero;
    if (z_zero_index < 0)
    	cout << "invalid calc_Leff function" << endl;
    return Leff;
}

double * read_initial_parameters(const char * FileName)
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

void write_final_parameters_float(const char * FileName, float * final_parameters)
{
	ofstream file;
	file.open(FileName);
	for (int i = 0; i < NP; ++i)
		file << final_parameters[i] << '\n';
	file.close();
}

void wipe_file(const char * wiped_file)
{
    fstream temp_file;
    temp_file.open(wiped_file,ios::out | ios::trunc);
    temp_file.close();
}

int main()
{
/* The main function reads in the data from all txt files of the form 
"EngeProd_*.**.txt" and creates three arrays for current, b_field, and z_pos. 
It creates a scatter plot of the data and fits a function to that plot. It then
prints the fitted parameters, and it saves a plot of the fit to the folder this
code is kept in. */
	input_datas WholeData;	  // creates a space for the data to be read to
	WholeData.Iamps   = (double*) malloc(NZ * sizeof(double));
	WholeData.z_pos   = (double*) malloc(NZ * sizeof(double));
	WholeData.b_field = (double*) malloc(NZ * sizeof(double));
	double FitParameters[NP];
	double FileIndices[] = {0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00};
	double Effective_Lengths[NF];
    float * chi_squared_array = (float*) malloc(sizeof(float)*max_iterations);
    float * del_chi_over_cpu_time = (float*) malloc(sizeof(float)*max_iterations);
	int j = 0;
	while (j < NF)
	{
	/* The while loop reads in all the data from the files and puts it into 
	three double arrays. */
		const char * fname = create_name_from_index(FileIndices[j]);
		float curr = FileIndices[j] * I_LIM;
		input_datas FileData = read_file_make_input_datas(fname);
		TGraph * graph = new TGraph(NZ, 
			FileData.z_pos, FileData.b_field);
		Effective_Lengths[j] = calc_Leff(FileData.z_pos,FileData.b_field,graph);
		if (j == 0)
		{
				WholeData.Iamps = FileData.Iamps;
				WholeData.z_pos = FileData.z_pos;
				WholeData.b_field = FileData.b_field;
		}
		else 
		{
			WholeData.Iamps = combine_arrays(WholeData.Iamps, FileData.Iamps,j);
			WholeData.z_pos = combine_arrays(WholeData.z_pos, FileData.z_pos,j);
			WholeData.b_field =
				combine_arrays(WholeData.b_field, FileData.b_field,j);
		}
		j = j + 1;
	}
	
	TGraph * leff_vs_curr = new TGraph(NF,currents_array,Effective_Lengths);
	
	TF1 * leff_poly = new TF1("leff_poly","[0]+[1]*x+[2]*x*x",0,200);
	double leff_pars[3] = {1.,1.,1.};
	for (int h = 0; h < 3; ++h)
		leff_poly -> SetParameter(h,leff_pars[h]);
	leff_vs_curr -> Fit("leff_poly");
	
	TGraph2D * ScatteredData = scatter(WholeData.z_pos,
		WholeData.Iamps,WholeData.b_field);
	double * p = (double*) malloc(sizeof(double)*NP);
	int iteration = 0;
	while (iteration < max_iterations)
	{
		TStopwatch watch;
		watch.Start();
		if (iteration == 0)
			p = read_initial_parameters("initial_parameters_global.txt");
		else
		    p = read_initial_parameters("parameters.txt");
		p[39] = leff_poly -> GetParameter(0);
		p[40] = leff_poly -> GetParameter(1);
		p[41] = leff_poly -> GetParameter(2);
		for (int i = 0; i < NP; ++i)
			cout << p[i] << endl;
		TF2 * func = new TF2("func", func_EngeProd_calc_leff, ZMIN, ZMAX, IMIN, IMAX,NP);
		for (int i = 0; i < NP; ++i)
			func -> SetParameter(i,p[i]);
		for (int i = (NP-3); i < NP; ++i)
			func -> SetParLimits(i,p[i],p[i]);
		for (int i = 0; i < (NP-3); ++i)
			func -> SetParLimits(i,-10,10);
		ScatteredData -> Fit("func");
		for (int i = 0; i < (NP); ++i)
			FitParameters[i] = func -> GetParameter(i);
        double timer = watch.CpuTime();
        if (func -> GetChisquare() > chi_squared_array[iteration-1])
        {
        	if (iteration == 0)
        		chi_squared_array[iteration] = func -> GetChisquare();
        	else
        	{
        		chi_squared_array[iteration] = chi_squared_array[iteration-1];
        		cout << "chi squared kept the same" << endl;
        	}
        }
        else 
        	chi_squared_array[iteration] = func -> GetChisquare();
        del_chi_over_cpu_time[iteration] = (chi_squared_array[iteration-1]-chi_squared_array[iteration])/timer;
		wipe_file("parameters.txt");
		write_final_parameters("parameters.txt", FitParameters);
		cout << iteration << endl;
		watch.Reset();
		iteration = iteration + 1;		
	}
	write_final_parameters("final_parameters_global_leff_calc.txt",FitParameters);
    wipe_file("chi_squared_vs_iterations.txt");
    write_final_parameters_float("chi_squared_vs_iterations.txt", chi_squared_array);
	write_final_parameters_float("del_chi_over_cpu_time.txt", del_chi_over_cpu_time);
    for (int i = 0; i < NP; ++i)
		cout << FitParameters[i] << endl;
	return 0;
}