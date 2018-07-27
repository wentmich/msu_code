/* This function takes in the values for "NP" Enge coefficients and fits a 
second order polynomial to the data set. It can then return the fitted 
parameters. 
NOTE: the program works for a text file of the following form...

18	0.47884	0.0017
36	0.24747	0.00029
54	0.25407	0.00022
72	0.25474	0.00026
90	0.2678	0.00021
108	0.26025	0.00029
126	0.26746	0.00029
144	0.26469	0.00033
162	0.259	0.00037
180	0.25611	0.00034

The first column is the current, the second is the Enge coefficient value, and 
the third is the error on the Enge coefficient value.
*/

// Global varibales
char fname[] = "EngeCoefs_p0.txt";
int number_of_lines = 10;
const int NP = 3;

struct file_data
{
// define structure to be used as output for read_file_to_array()
	double * Iamps;    // list of currents
	double * coef;     // list of the Enge coefficient values
	double * error;    // list of the errors in value of the Enge coefficient
};

double func_second_poly(double *I, double *params)
{
// defines and returns function with three paramters (a second order polynomial)
	func = (params[0] + params[1]*I[0] + params[2]*pow(I[0],2));
	return func;
}

file_data read_file_to_array(char name[] = "EngeCoefs_p0.txt")
{
/* Reads in a file and creates three arrays from the three columns of the file.
The arrays are for current, enge coefficient, and error. These arrays are 
returned in the form of a structure "file_data." */
	const int Nline = 256;            // max number of characters to be read in
	char string_location[Nline];      // creates a location to read file line
	file_data file_arrays;
	FILE * file;
	file = fopen(name, "r");
	file_arrays.Iamps = (double*) malloc(number_of_lines * sizeof(double));
	file_arrays.coef  = (double*) malloc(number_of_lines * sizeof(double));
	file_arrays.error = (double*) malloc(number_of_lines * sizeof(double));
	for ( int i = 0; i < number_of_lines; ++i)
	{
		// reads in data from the file line by line into the file_arrays struct
		fgets(string_location, Nline, file);
		sscanf(string_location, " %lg %lg %lg\n", &file_arrays.Iamps[i], 
			&file_arrays.coef[i], &file_arrays.error[i]);
	}
	fclose(file);
	return file_arrays;
}

double * scatter_data_from_array(double * a, double * b, double * c, double * d)
{
/* Plots data as a scatter plot (a = x_axis, b = y_axis, c = x_error,
d = y_error). The data is taken from the first two input arrays, and errors are 
taken from the second two. It then returns an array of the fitted parameters. */
	double * fitted_parameters;
	fitted_parameters = (double*) malloc(sizeof(double)*NP);
	TGraphErrors * data = new TGraphErrors(number_of_lines, a, b, c, d); 
	data -> SetTitle("Enge Coefficient Plot;I (amps);Enge Coefficient");
	data -> SetLineColor(0);      // line will be white
	data -> SetMarkerStyle(8);    // markers are circles
	data -> Draw();
	// fits data with a second degree polynomial as defined above	
	double params[NP] = {1.,1.,1.};
	TF1 * func = new TF1("func", func_second_poly, 0, 200, NP);
	for ( int i = 0; i < NP; ++i)
		func -> SetParameter(i, params[i]);
	data -> Fit("func");
	func -> Draw("same");
	for (int i = 0; i < NP; ++i)
		fitted_parameters[i] = func -> GetParameter(i);
	return fitted_parameters;
}

// main function calling on previously defined functions
void main(char name[] = fname)
{
/* Reads in data for the specified function, plots it, and then
fits a polynomial to it. */
	file_data file_arrays;
	file_arrays = read_file_to_array(name);
	const int var = number_of_lines;
	double error_x[var];
	double * parameters;
	parameters = (double*) malloc(sizeof(double)*NP);
	parameters = scatter_data_from_array(file_arrays.Iamps, file_arrays.coef,
		error_x, file_arrays.error);
}