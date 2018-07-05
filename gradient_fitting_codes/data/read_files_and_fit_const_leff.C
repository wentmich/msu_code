/* The program reads in data from a text file and creates a scatter plot. It 
calculates the effective length of the magnet which the data pertains to. Then
it fits an Enge function to the scattered data and returns the fitted parameters
along with a plot of the fitted function. The function is used for local fits
with a constant effective length.
Other than effective length calculations, this program is the same as 
"read_files_and_fit_calc_leff.C".
NOTE: This works for files of the following form...

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
 
The function performs these tasks for one file at a time. It is a local fit. */

/* Defines all the global variables to be used in the program */
double Leff, R0, Iamps;
const int NP = 13;
R0 = 0.2;              
double ylim = 2.2;
double zmax=1., zmin=-zmax, delz = 0.005;   int Nz;
double *za, *fa;
double zz[1], ftmp;
Nz = 1 + (int)( (zmax-zmin)/delz );

double func_EngeProd(double *z, double *p)
{
/* creates a function to fit to the scattered data. */
  double zi =  z[0] - Leff/2.0;				// defines z variable
  double ze = -z[0] - Leff/2.0;				// defines z variable
  zi /= 2*R0; ze /= 2*R0;
  double fi=1./(1.+exp(p[0]+zi*(p[1]+zi*(p[2]+zi*(p[3]+zi*(p[4]+zi*p[5]))))));
  double fe=1./(1.+exp(p[6]+ze*(p[7]+ze*(p[8]+ze*(p[9]+ze*(p[10]+ze*p[11]))))));
  return p[12]*fi*fe;
}

const char * create_name_from_index(double index = 0.10)
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
	const char * fname;
	fname = (char*) malloc((n+1) * sizeof(char)); // space for string + null
	strcpy(fname, FileName.c_str());
	
	return fname;
}

void read_single_file_and_fit(double index = 0.10)
{
/* function reads in a file of the form shown above. It then scatters the data
and fits a function to the data. The file it reads in is "EngeProd_*.**.txt".
For an explanation of the fitting process and the TF1 class, see 
5.3.3 of the ROOTUsersGuideLetter2014.pdf. See also 
https://root.cern.ch/root/HowtoFit.html*/
	const char * fname = create_name_from_index(index); // creates file name
	const int NLINE=256;
	FILE *fp;                   //file pointer
	char line[NLINE];           //buffer to read lines
	printf(" points per excitation, Nz= %i \n", Nz);
	za = (double*) malloc( Nz * sizeof(double) );	// creates z array
	fa = (double*) malloc( Nz * sizeof(double) );	// creates field array

	//Read field gradient data from file. Generate scatter plot.
    fp = fopen(fname, "r");
    fgets(line, NLINE, fp);   sscanf(line, " %lg ", &Iamps);
    fgets(line, NLINE, fp);   sscanf(line, " %lg ", &R0);
    fgets(line, NLINE, fp);   sscanf(line, " %i ", &Nz);
    fgets(line, NLINE, fp);
    printf(" Iamps= %lg, R0= %lg, Nz= %i \n", Iamps, R0, Nz);
    for( int k=0; k<Nz; k++)
    {
      fgets(line, NLINE, fp);						// reads data to z,f arrays
      sscanf(line," %lg %lg\n", &za[k], &fa[k]);
    }
    fclose(fp);

    TGraph * data = new TGraph(Nz, za, fa);// creates scatter plot of z vs field
    data->SetTitle("field; z; f");	
    data->SetLineColor(0);
    data->SetMarkerStyle(1);
    data->Draw();					// plots data on canvas
   
    // find effective length using integral and b(z=0)
    Leff = 0.7;
   
    // Define the fuction to fit
    double p[NP] =
		{0., 3., 0., 0., 0., 0.    //0-5, Enge incident side
		,0., 4., 0., 0., 0., 0.    //6-11, exit
		,0.3};                     //12, B0
    TF1 *func = new TF1( "func", func_EngeProd, zmin,zmax, NP);
    for( int i=0; i<NP; i++)
    	func->SetParameter(i, p[i]); // initialize the parameters
     
    // Optional: Fit the parameters to data 
    bool fit_TF1=1;  //Set to 1 to allow fitting
    if(fit_TF1)
    	data->Fit("func");
    func->Draw("same"); // plot the fitted function with the scatter plot 
  
    int k = 2; // sets the number of decimal precision to be read into string
	stringstream stream;
	stream << fixed << setprecision(k) << index; // makes "index" a string
	string mystring = stream.str();  
	
    double coefficients[NP]; // creates an array of the fitted parameters
    double errors[NP];	     // creates an array of the errors
    for (int i = 0; i < NP; ++i)
    	coefficients[i] = func -> GetParameter(i);
    for (int i = 0; i < NP; ++i)
    	errors[i] = func -> GetParError(i);
    // print fitted parameters to text file with effective length
    ofstream myfile;
    string myfile_front ( "Current_" );
    string myfile_back  ( "_coefficients.txt" );
    string temp_name = myfile_front + mystring + myfile_back;
    const int m = temp_name.length();
    char myfile_name[m+1];
    strcpy(myfile_name, temp_name.c_str());
    myfile.open(myfile_name);
    myfile << "Iamps = " << 180 * index << "," << '\t' << "R0 = " 
    	<< R0 << "," << '\t' << "Nz = " << Nz << "," << '\t' 
    	<< "Leff = " << Leff << '\r';
    for (int i = 0; i < NP; ++i)
    	myfile << fixed << setprecision(5) << "p" << i << '\t' 
    		<< coefficients[i] << '\t' << errors[i] << '\r';
    myfile.close();
}

void main()
{
/* main function reads, scatters, and fits the data for all the current files in
the list "list_of_files". */
    double list_of_files[] = {0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00};
    for (int i = 0; i < 10; ++i)
        read_single_file_and_fit(list_of_files[i]);
}