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

//int Nz = 401;

struct data
{
	double * z_pos;
	double * grad;
	data()
	{
		z_pos = (double*) malloc(sizeof(double)*Nz);
		grad  = (double*) malloc(sizeof(double)*Nz);
	}
}read_data;
char * f_name = (char *)"/projects/fribmappers_data/analysis_MW1/output/b_n06_m0.out";
char * file_pars = (char *)"/projects/aris/MP/field_params_DEV/param_files_2016/FSQ5/M5_params_FSQ5_9.4Tpm_TEST.txt";
char * split_string;
//const char * delim = " ,";
//const int NP = 13;
//double R0 = 1.;
//double LEFF = 1.;
//double e = 2.71828;

char * make_file_name(int file_number)
{
	if (file_number == 0)
		return ((char *)"/projects/fribmappers_data/analysis_MW1/b_n00_m0.txt");
	char * out_name = (char *) malloc(sizeof(char)*60);
	sprintf(out_name, (char *) "/projects/fribmappers_data/analysis_MW1/b_n%02d_m0.txt", file_number);
	return out_name;
}

double * read_initial_pars(char * name = (char *)"")
{
	double * pars = (double*) malloc(sizeof(double)*(NC + 3));
	double * dump = (double*) malloc(sizeof(double)*1);
	const int Nline = 600;
	char line[Nline];
	FILE * file;
	file = fopen(name, "r");
	for (int i = 0; i < 4; ++i)
		fgets(line, Nline, file);
	fgets(line, Nline, file); sscanf(line, "%lg", &pars[0]);
	fgets(line, Nline, file); sscanf(line, "%lg", &pars[1]);
	for (int i = 0; i < 3; ++i)
		fgets(line, Nline, file);
	sscanf(line, "%lg", &pars[14]);
	for (int i = 0; i < 6; ++i)
		fgets(line, Nline, file);
	fgets(line, Nline, file);
	sscanf(line, "%lg %lg %lg %lg %lg %lg %lg %lg", &dump[0], &dump[0],
		&pars[2], &pars[3], &pars[4], &pars[5], &pars[6], &pars[7]);
	fgets(line, Nline, file);
	sscanf(line, "%lg %lg %lg %lg %lg %lg %lg %lg", &dump[0], &dump[0],
		&pars[8], &pars[9], &pars[10], &pars[11], &pars[12], &pars[13]);
	return pars;
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
  double fi=1./(1.+pow(e,(p[0]+zi*(p[1]+zi*(p[2]+zi*(p[3]+zi*(p[4]+zi*p[5])))))));
  double fe=1./(1.+pow(e,(p[6]+ze*(p[7]+ze*(p[8]+ze*(p[9]+ze*(p[10]+ze*p[11])))))));
  return p[12]*fi*fe;
}

data read_gradient_data(char * input_name)
{
	const char * fname = input_name;
	//printf(fname);
	data read_data;
	const int Nline = 100;
	char line[Nline];
	FILE * file;
	file = fopen(fname, "r");
	fgets(line,Nline,file);
    for (int i = 0; i < Nz; ++i)
    {
    	fgets(line,Nline,file);
    	split_string = strtok(line,sep);
    	sscanf(split_string,"%lg",&read_data.z_pos[i]);
    	split_string = strtok(NULL,sep);
    	sscanf(split_string,"%lg",&read_data.grad[i]);
    		
    	//sscanf(line,"%lg %lg \n", &read_data.z_pos[i], &read_data.grad[i]);
    }
    for (int i = 0; i < Nz; ++i)
    {
    	read_data.grad[i] = read_data.grad[i];
    	//cout << read_data.z_pos[i] << endl;
    }
    return read_data;
}

TGraph * graph_data(data read_data)
{
	TGraph * graph = new TGraph(Nz, read_data.z_pos, read_data.grad);
    graph -> SetMarkerStyle(0);
    return graph;
}

TF1 * fit_graph(TGraph * graph, double * params)
{
	R0 = params[0];
    LEFF = params[1];
    TF1 * func = new TF1("func", func_EngeProd, -1,1,NP);
    for (int i = 2; i < NC + 3; ++i)
    {
    	// make the fitted Enge function product third degree
    	func -> SetParameter(i-2, params[i]);
    	if (i-2==4 || i-2==5 || i-2==10 || i-2==11)
    		func -> FixParameter(i-2,0);
    }
    graph -> Fit("func");
    return func;
}

double * combine_arrays(double * array1, double * array2, int k)
{
/* Takes two arrays of the length specified by "Nz" and replaces 
them into a single array. Returns a pointer to the combined array." */
	double * NewArray;
	NewArray = (double*) malloc(sizeof(double)*((NC + 3)*k) 
		+ sizeof(double)*(NC + 3));
	int Nelements1 = Nz*k;
	int Nelements2 = Nz;
	for (int i = 0; i < Nelements1; ++i)
		NewArray[i] = array1[i];
	for (int i = Nelements1; i < (Nelements1 + Nelements2); ++i)
		NewArray[i] = array2[i - Nelements1];
	return NewArray;
}

int write_M5_params_NORM(double * parameters)
{
    ofstream file;
    file.open("M5_sim_fit.txt");
    file << "#Field excitation parameters for each multipole n \r\n";
    file << "# row1: reference radius r0 [m] \r\n";
    file << "# row2: LEFF[m] \r\n";
    file << "# row3: Field excitation on multipole \r\n";
    file << '\t' << parameters[0] << "\r\n";
    file << '\t' << parameters[1];
    file << "\r\n#Fields at r0 for each multipole component \r\n";
    file << "#n= \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 \t 10 \r\n";
    for (int i = 2; i < 11; ++i)
    	file << '\t' << parameters[((NC + 3)*i)+14];
    file << "\r\n";
    file << "# Enge coefs \r\n";
    file << "# n \t ee \t c0 \t c1 \t c2 \t c3 \t c4 \t c5 \t c6 \t c7 \t c8 \t c9 \t c10 \t c11 \r\n";
    for (int i = 0; i < 11; ++i)
    {
    	file << i << '\t' << "0" << '\t';
    	for (int k = 0; k < ((NP-1)/2); ++k)
    		file << parameters[(i*(NC + 3))+(k+2)] << '\t';
    	for (int k = 0; k < (NP-1)/2; ++k)
    		file << "0 \t";
    	file << "\r\n";
    	
    	file << i << '\t' << "1" << '\t';
    	for (int k = (int)(NP-1)/2; k < (NP-1); ++k)
    		file << parameters[(i*(NC + 3))+(k+2)] << '\t';
    	for (int k = 0; k < (NP-1)/2; ++k)
    		file << "0 \t";
    	file << "\r\n";
    }    	
	
    file.close();
    return 0;
}

void plot_and_fit_gradient(int num, double * all)
{
	TCanvas * c1 = new TCanvas();
	f_name = make_file_name(num);
	double * params = (double*) malloc(sizeof(double)*(NC + 3));
	params = read_initial_pars(file_pars);	
	read_data = read_gradient_data(f_name);
	TGraph * graph = graph_data(read_data);
	TF1 * func = fit_graph(graph, params);
	graph -> Draw("p0");
	func -> SetLineColor(5);
	func -> Draw("same");
	//c1 -> SaveAs("c%i.png",num);
	all[num*(NC + 3)] = R0;
	all[num*(NC + 3) + 1] = LEFF;
	for (int i = 2; i < (NC + 3); ++i)
		all[num*(NC + 3) + i] = func -> GetParameter(i-2);
}

double * fit()
{	
	int Nn = 11;
	double * all_params = (double*) malloc(sizeof(double)*(NC + 3)*Nn);
	for(int n=0; n<Nn; n++)
	{
  	    plot_and_fit_gradient(n, all_params);
  	    cout << "completed cycle: " << n << endl;
  	} 
  	return all_params;
}