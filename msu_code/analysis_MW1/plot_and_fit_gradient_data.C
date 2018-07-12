# include <stdio.h>
# include <math.h>
# include <string.h>
# include <time.h>
# include <sstream>
# include <TCanvas.h>
# include <TF1.h>
# include <TGraph.h>
# include <fstream>

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
int index;

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

char * make_file_name(int file_number)
{
	if (file_number == 0)
		return ((char *)"/projects/fribmappers_data/analysis_MW1/b_n00_m0.txt");
	char * out_name = (char *) malloc(sizeof(char)*60);
	sprintf(out_name, (char *) "/projects/fribmappers_data/analysis_MW1/b_n%02d_m0.txt", file_number);
	return out_name;
}

TGraph * graph_data(data read_data)
{
	TGraph * graph = new TGraph(Nz, read_data.z_pos, read_data.grad);
	graph -> SetMarkerStyle(0);
	return graph;
}

double calc_LEFF(double * z_data, double * b_field)
{
/* function to determine the effective length using the definition 
LEFF = integral(b_field)/(b_field at z=0). The input arrays should each have NZ 
elements. */
	double b_zero;
	double integral;
    int z_zero_index = 0;
    for (int i = 0; i < Nz; i++){
    // finds index of z = 0
        if (fabs(z_data[i]) < 1.0e-10)
        z_zero_index = i;}
    TGraph * data = new TGraph(Nz, z_data, b_field); // plots arrays
    b_zero = b_field[z_zero_index];
    integral = data -> Integral();
    LEFF = integral/b_zero;
    return LEFF;
}

void * fit_graph(TGraph * graph, double * coefficients, double * bfield, int pole)
{
	TCanvas * c1 = new TCanvas();
	int NP = NC + 1;
	TF1 * func = new TF1("func", func_EngeProd, -1, 1, NP);
	for (int i = 0; i < (NP - 1)/2; ++i)
	{
		func -> SetParameter(i, coefficients[(NC*NEE*pole) + i]);
		if (i == 4 || i == 5)
			func -> FixParameter(i,0);
	}
	for (int i = (NP - 1)/2; i < (NP - 1); ++i)
	{
		func -> SetParameter(i, coefficients[(NC*NEE*pole) + NC + i]);
		if (i == 10 || i == 11)
			func -> FixParameter(i,0);
	}
	func -> SetParameter(NP-1, bfield[pole]);
	graph -> Fit("func");
	for (int i = 0; i < (NP - 1)/2; ++i)
		coefficients[(NC*NEE*pole) + i] = func -> GetParameter(i);
	for (int i = 0; i < (NP - 1)/2; ++i)
		coefficients[(NC*NEE*pole) + NC + i] = func -> GetParameter(i + (NP - 1)/2);
	bfield[pole] = func -> GetParameter(NP-1);
	graph -> Draw();
	func -> Draw("same");
	return;
}

void write_M5_params_NORM(double * bfield, double * parameters)
{
    ofstream file;
    file.open("M5_sim_fit.txt");
    file << "#Field excitation parameters for each multipole n \r\n";
    file << "# row1: reference radius r0 [m] \r\n";
    file << "# row2: LEFF[m] \r\n";
    file << "# row3: Field excitation on multipole \r\n";
    file << '\t' << R0 << "\r\n";
    file << '\t' << LEFF;
    file << "\r\n#Fields at r0 for each multipole component \r\n";
    file << "#n= \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 \t 10 \r\n";
    for (int i = 2; i < Nn; ++i)
    	file << '\t' << bfield[i];
    file << "\r\n";
    file << "# Enge coefs \r\n";
    file << "# n \t ee \t c0 \t c1 \t c2 \t c3 \t c4 \t c5 \t c6 \t c7 \t c8 \t c9 \t c10 \t c11 \r\n";
    for (int i = 0; i < Nn; ++i)
    {
    	file << i << '\t' << "0" << '\t';
    	for (int k = 0; k < NC; ++k)
    		file << parameters[NEE*i*NC + k] << '\t';
    	file << "\r\n";
    	
    	file << i << '\t' << "1" << '\t';
    	for (int k = NC; k < NC*2; ++k)
    		file << parameters[NEE*i*NC + k] << '\t';
    	file << "\r\n";
    }    	
	file.close();
    return;
}

data read_gradient_data(char * input_name)
{
	const char * fname = input_name;
	data read_data;
	const int Nline = 100;
	char line[Nline];
	FILE * file;
	file = fopen(fname, "r");
	fgets(line,Nline,file);
	char * split_string;
    for (int i = 0; i < Nz; ++i)
    {
    	fgets(line,Nline,file);
    	split_string = strtok(line,sep);
    	sscanf(split_string,"%lg",&read_data.z_pos[i]);
    	split_string = strtok(NULL,sep);
    	sscanf(split_string,"%lg",&read_data.grad[i]);
    }
    for (int i = 0; i < Nz; ++i)
    {
    	read_data.grad[i] = read_data.grad[i];
    }
    return read_data;
}

void read_M5_params_NORM( char *fname, double *Bn_)
{
	int n, ee;
	int i;
	FILE *fp; char line[512]; char *pch; //const char sep[2]=" ,";
	fp = fopen(fname, "r");
	sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
	sscanf(line," %lf\n", &R0);
	fgets(line,512,fp);
	sscanf(line," %lf\n", &LEFF);
	cout << R0 << '\t' << LEFF << endl;
	sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
	i=2;
	pch = strtok(line, sep);
	cout << pch << endl;
	cout << "!!!!!!!!!!!!!" << endl;
	sscanf(pch,"%lf", &Bn_[i]);
	for( i=3; i<7 ; i++ ){
	  pch = strtok(NULL, sep); sscanf(pch,"%lf", &Bn_[i]);
	}
	sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
	cout << "!!!!!!!!!!!!!" << endl;
	while( !feof(fp) ){
	pch = strtok(line, sep); sscanf(pch,"%i", &n);
	pch = strtok(NULL, sep); sscanf(pch,"%i", &ee);
	for(int c=0; c<NC; c++){
	  index = c + NC*( NEE*n + ee );
	  pch = strtok(NULL, sep); sscanf(pch,"%lf", &fcoef[index]);
	}
	fgets(line,512,fp);
	}
	pch = strtok(line, sep); sscanf(pch,"%i", &n);
	pch = strtok(NULL, sep); sscanf(pch,"%i", &ee);
	for( c=0; c<NC; c++){
	  index = c + NC*( NEE*n + ee );
	  pch = strtok(NULL, sep); sscanf(pch,"%lf", &fcoef[index]);
	}
	return;
}

int main()
{
	char f_dir[512]="/projects/aris/MP/field_params_DEV/param_files_2016/FSQ5/";
	char f_params[128] = "M5_params_FSQ5_9.4Tpm_TEST.txt";
	char file_pars[640];
	sprintf(file_pars, "%s%s", f_dir, f_params);
	Bn = (double*) malloc( Nn * sizeof(double) );
	fcoef = (double*) malloc( NC*NEE*NNC * sizeof(double) );
	An = (double*) malloc( Nn * sizeof(double) );
	read_M5_params_NORM(file_pars, Bn);
	for (int num = 0; num < Nn; ++num)
	{
		read_data = read_gradient_data(make_file_name(num));
		LEFF = calc_LEFF(read_data.z_pos, read_data.grad);
		TGraph * graph = graph_data(read_data);
		fit_graph(graph, fcoef, Bn, num);
		cout << LEFF << endl;
	}
	for (int i = 0; i < 256; ++i)
		cout << fcoef[i] << endl;
	write_M5_params_NORM(Bn, fcoef);
	return 0;
}