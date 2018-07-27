// Fit a second order polynomial to each Enge coefficient as a function of current

// Global varibales
char fname[] = "EngeCoefs_p0.txt";
int number_of_lines = 10;
const int Np = 4;

struct arrays
	{
		double * Iamps;
		double * coef;
		double * error;
	};
// defines a 1D function with three paramters (a second order polynomial)

double func_second_poly(double *I, double *params)
{
	func = (params[0] + params[1]*I[0] + params[2]*pow(I[0],2) + params[3]*pow(I[0],3));
	return func;
}

// reads file and creates two arrays from the two file columns for current and associated Enge coefficient
arrays read_file_to_array(char name[] = "EngeCoefs_p0.txt")
{
	const int Nline = 256;
	char string_location[Nline];
	arrays file_arrays;
	FILE * file;
	file = fopen(name, "r");
	file_arrays.Iamps = (double*) malloc(number_of_lines * sizeof(double));
	file_arrays.coef  = (double*) malloc(number_of_lines * sizeof(double));
	file_arrays.error = (double*) malloc(number_of_lines * sizeof(double));
	for ( int i = 0; i < number_of_lines; ++i)
	{
		fgets(string_location, Nline, file);
		sscanf(string_location, " %lg %lg %lg\n", &file_arrays.Iamps[i], &file_arrays.coef[i], &file_arrays.error[i]);
	}
	fclose(file);
	return file_arrays;
}
// now we have the data needed in two separate arrays
// plotting the data as a scatter plot using the two arrays
void scatter_data_from_array(double * a, double * b, double * c, double * d)
{
	//TGraph * data = new TGraph(number_of_lines, a, b);
	
	TGraphErrors * data = new TGraphErrors(number_of_lines, a, b, c, d);
	if ((TCanvas*)gROOT -> FindObject("cCOEF")){}
	else
	{
		TCanvas * cEFIT = new TCanvas("cEFIT","cEFIT",40,40,700,450);
	}
	data -> SetTitle("Enge Coefficient Plot (Calculated Leff);I (amps);Enge Coefficient");
	//data -> SetTitle("Calculated Effective Length vs. Current; I (amps); Effective Length (m)");
	data -> SetLineColor(1);
	data -> SetMarkerStyle(8);
	data -> Draw("ap");
	
	cEFIT -> Update();
	
	double params[Np] = {1.,1.,1.,1.};
	TF1 * func = new TF1("func", func_second_poly, 0, 200, Np);
	for ( int i = 0; i < Np; ++i)
	{
		func -> SetParameter(i, params[i]);
	}
	data -> Fit("func");
	func -> Draw("same");
	
	
}
// main function calling on previously defined functions
void main(char name[] = "EngeCoefs_p0.txt")
{
	arrays file_arrays;
	file_arrays = read_file_to_array(name);
	double * Iamps, * coef, * error;
	Iamps = file_arrays.Iamps;
	coef  = file_arrays.coef;
	error = file_arrays.error;
	const int var = number_of_lines - 1;
	double error_x[var];
	scatter_data_from_array(Iamps,coef,error_x,error);
}