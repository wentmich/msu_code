/* This program reads in on-axis gradient data created by the program 
"FieldSim_derivatives_V03.C." and fits Enge product functions to the data. The
files it reads from are "/projects/fribmappers_data/analysis_MW1/b_n**d_m0.txt,"
and are written by the above program. */

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
    index   = index used in read_M5_params_NORM()  */
    
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

/* DECLARATION OF DATA STRUCTURES:
    data     = structure used to read in the data from the text files */

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

/* DECLARATION OF FUNCTIONS: */

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
/* creates a file name of the form 
/projects/fribmappers_data/analysis_MW1/b_n**d_m0.txt. This is used to read in 
the on-axis gradient data. */
    if (file_number == 0)
        return ((char *)"/projects/fribmappers_data/analysis_MW1/b_n00_m0.txt");
    char * out_name = (char *) malloc(sizeof(char)*60);
    sprintf(out_name, (char *) "/projects/fribmappers_data/analysis_MW1/b_n%02d_m0.txt", file_number);
    return out_name;
}

TGraph * graph_data(data read_data)
{
/* creates a graph object from the two arrays comprising the "data" structure 
type. It plots data.z_pos on the x-axis and data.grad on the y-axis. */
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
    for (int i = 0; i < Nz; i++)
    {
        if (fabs(z_data[i]) < 1.0e-10)
            z_zero_index = i;
    }
    TGraph * data = new TGraph(Nz, z_data, b_field); // plots arrays
    b_zero = b_field[z_zero_index];
    integral = data -> Integral();
    LEFF = integral/b_zero;
    return LEFF;
}

void fit_graph(TGraph *graph, double *coefficients, double *bfield, int pole, bool fit)
{
/* fits an Enge product model to a graph. The initial values of the Enge model
are taken from "coefficients" and "bfield". The Enge coefficients come from 
"coefficients," and the b_field constant comes from "bfield." The integer "pole"
tells the fitter which coefficients and which bfield value to take. Although the
model function is actually fifth order, here it has been set to third order by 
fixing the fourth and fifth order parameters to zero. */
    TCanvas * c1 = new TCanvas();
    int NP = NC + 1;
    TF1 * func = new TF1("func", func_EngeProd, -1, 1, NP);
    for (int i = 0; i < (NP - 1)/2; ++i)
    {
        func -> SetParameter(i, coefficients[(NC*NEE*pole) + i]);
        if (i == 4 || i == 5)
            func -> FixParameter(i,0);
        //cout << coefficients[(NC*NEE*pole) + i] << endl;
    }
    for (int i = (NP - 1)/2; i < (NP - 1); ++i)
    {
        func -> SetParameter(i, coefficients[(NC*NEE*pole) + (NC/2) + i]);
        if (i == 10 || i == 11)
            func -> FixParameter(i,0);
        //cout << coefficients[(NC*NEE*pole) + (NC/2) + i] << endl;
    }
    func -> SetParameter(NP-1, bfield[pole]);
    if (fit)
    {
		graph -> Fit("func");
		for (int i = 0; i < (NP - 1)/2; ++i)
			coefficients[(NC*NEE*pole) + i] = func -> GetParameter(i);
		for (int i = 0; i < (NP - 1)/2; ++i)
			coefficients[(NC*NEE*pole) + NC + i] = func -> GetParameter(i+(NP-1)/2);
		bfield[pole] = func -> GetParameter(NP-1);
	}
    graph -> Draw();
    func -> Draw("same");
    return;
}

void write_M5_params_NORM(double * bfield, double * parameters)
{
/* writes the final fitted coefficients and bfield values to a file named
"M5_sim_fit.txt." The parameters are written to copy the format of the file from
which the initial parameters are read. */
    ofstream file;
    file.open("M5_sim_fit.txt");
    file << "#Field excitation parameters for each multipole n \r\n";
    file << "# row1: reference radius r0 [m] \r\n";
    file << "# row2: LEFF[m] \r\n";
    file << "# row3: Field excitation on multipole \r\n";
    file << "   " << R0 << "\r\n";
    file << "   " << LEFF;
    file << "\r\n#Fields at r0 for each multipole component \r\n";
    file << "#n=  2        3     4     5     6     7     8     9     10\r\n";
    for (int i = 2; i < Nn; ++i)
        file << "     " << bfield[i];
    file << "\r\n";
    file << "# Enge coefs \r\n";
    file << "# n  ee  c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11\r\n";
    for (int i = 0; i < Nn; ++i)
    {
        file << "  " << i << "  0  ";
        for (int k = 0; k < NC; ++k)
            file << parameters[NEE*i*NC + k] << " ";
        file << "\r\n";
        
        file << "  " << i << "  1  ";
        for (int k = NC; k < NC*2; ++k)
            file << parameters[NEE*i*NC + k] << " ";
        file << "\r\n";
    }       
    file.close();
    return;
}

data read_gradient_data(char * input_name)
{
/* reads the on-axis gradient data from the file specified in the main funciton.
The function reads it into a "data" type data structure. NOTE: the files from 
which it is reading contain both real and imaginary components of the gradient. 
In this function, the imaginary components are ignored. */
    const char * fname = (const char *) input_name;
    data file_data;
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
        sscanf(split_string,"%lg",&file_data.z_pos[i]);
        split_string = strtok(NULL,sep);
        sscanf(split_string,"%lg",&file_data.grad[i]);
    }
    /*for (int i = 0; i < Nz; ++i)
    {
        file_data.grad[i] = file_data.grad[i]*(10e13);
    }*/
    return file_data;
}

void read_M5_params_NORM( char *fname, double *Bn_)
{
/* reads in the initial parameters values for the fit from the file 
"M5_params_FSQ5_9.4Tpm_TEST.txt." These parameters are read in the same way they
are read in "FieldSim_derivatives_V03.C." The Enge coefficients are read into 
the array "fcoef," and the field parameters are read into the array specified as 
the "Bn_" argument. */
    int n, ee;
    int i;
    FILE *fp; char line[512]; char *pch; //const char sep[2]=" ,";
    fp = fopen(fname, "r");
    sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
    sscanf(line," %lf\n", &R0);
    fgets(line,512,fp);
    sscanf(line," %lg\n", &LEFF);
    sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
    i=2;
    pch = strtok(line, sep);
    sscanf(pch,"%lf", &Bn_[i]);
    for( i=3; i<7 ; i++ ){
      pch = strtok(NULL, sep); sscanf(pch,"%lf", &Bn_[i]);
    }
    sprintf(line,"#"); while(strncmp(line,"#",1)==0) fgets(line,512,fp);
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

/* DECLARATION OF MAIN FUNCTION TO RUN THE FITTING PROCEDURE: */

int plot_and_fit_gradient_data()
{
/* main function to run from root. It declares the file name from which the
initial parameters should be read. It then reads the paramters, creates graphs
for each of the multipole components, and fits Enge product models to them. A 
flag can be set to plot the function with initial parameters instead of fitted
parameters. */
    char f_dir[512]="/projects/aris/MP/field_params_DEV/param_files_2016/FSQ5/";
    char f_params[128] = "M5_params_FSQ5_9.4Tpm_TEST.txt";
    char file_pars[640];
    sprintf(file_pars, "%s%s", f_dir, f_params);
    Bn = (double*) malloc( Nn * sizeof(double) );
    fcoef = (double*) malloc( NC*NEE*NNC * sizeof(double) );
    An = (double*) malloc( Nn * sizeof(double) );
    read_M5_params_NORM(file_pars, Bn);
    double LEFF_const = LEFF; // initialize
    // Evaluate LEFF from dominant multipole term, nLEFF
    // using integration on the gradient
    int nLEFF = 2;
    read_data = read_gradient_data(make_file_name(nLEFF));
    LEFF = calc_LEFF(read_data.z_pos, read_data.grad);
    // Loop through all gradients and fit field (e.g. Enge) parameters on each
    bool fit = 1;
    // fitting flag. if true, complete fits. if false, only plot the initial parameters
    
    for (int pole = 1; pole < Nn; ++pole)
    {
        read_data = read_gradient_data(make_file_name(pole));
        TGraph * graph = graph_data(read_data);
        fit_graph(graph, fcoef, Bn, pole, fit);
        printf(" pole = %i,  LEFF = %lf\n", pole, LEFF);
    }
    write_M5_params_NORM(Bn, fcoef);
    return 0;
}