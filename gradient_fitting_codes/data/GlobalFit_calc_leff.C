/* This program starts by reading in the data from "NF" files containing the 
z_position and the magnetic field for a series of points within the magnet
at a fixed radius "r0." The data from the "NF" files is then placed into a 
single large structure of three arrays, InputData. These arrays are plotted into
a three dimensional scatter plot. The function fits a 2D Enge function to the 
scatter plot. The Enge function is a function of z_position and current. The 
effective length is kept constant here. 
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
magnetic field. */

// Global variables
const double LEFF, R0, I_AMPS, ZMAX, ZMIN, IMAX, IMIN, e;
const int NUMBER_OF_POINTS, NF, TOTAL_DATA_POINTS, I_LIM, NP, NZ;

LEFF = 0.7; R0 = 0.2; ZMIN = -1.; ZMAX = -ZMIN; IMAX = 200.; IMIN = 0.; 
e = 2.718281828459045;

NUMBER_OF_POINTS = 401; NF = 10; TOTAL_DATA_POINTS = NUMBER_OF_POINTS * NF;
I_LIM = 180; NP = 42; NZ = 401;

double currents_array[NF] = {18.,36.,54.,72.,90.,108.,126.,144.,162.,180.};

// Declaration of structure
struct input_data
{
	double * Iamps;
	double * z_pos;
	double * b_field;
};

input_data read_file_make_input_data(const char * FileName)
{
/* The function reads in data from a file in the form shown above and creates a
structure "new_data" in the form of "input_data." Each column of the file is
read to an array "z_pos" or "b_field," and "Iamps" is filled with the same value
of the current. */
	const int Nline = 256;               // maximum characters read at a time
	char string_location[Nline];         // place to store characters
	double current[1];                   // value of the current for the file
	input_data new_data;                 // structure to be returned
	FILE * file;
	file = fopen(FileName, "r");         // opens the file
	new_data.Iamps   = (double*) malloc(NUMBER_OF_POINTS * sizeof(double));
	new_data.z_pos   = (double*) malloc(NUMBER_OF_POINTS * sizeof(double));
	new_data.b_field = (double*) malloc(NUMBER_OF_POINTS * sizeof(double));
	fgets(string_location, Nline, file); 
	sscanf(string_location, "%lg", &current[0]);
	fgets(string_location, Nline, file); // opens and the lets go of lines 2-4
	fgets(string_location, Nline, file);
	fgets(string_location, Nline, file);
	for ( int i = 0; i < NUMBER_OF_POINTS; ++i)
	{
		// opens each line and reads in the two values for z_pos and b_field
		fgets(string_location, Nline, file);
		sscanf(string_location, "%lg %lg\n",
			&new_data.z_pos[i], &new_data.b_field[i]);
		new_data.Iamps[i] = current[0];//sets each current entry to file current
	}
	return new_data;
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

double * combine_arrays(double * array1, double * array2, int k)
{
/* Takes two arrays of the length specified by "NUMBER_OF_POINTS" and replaces 
them into a single array. Returns a pointer to the combined array." */
	double * NewArray;
	NewArray = (double*) malloc(sizeof(double)*(NUMBER_OF_POINTS*k) 
		+ sizeof(double)*NUMBER_OF_POINTS);
	int Nelements1 = NUMBER_OF_POINTS*k;
	int Nelements2 = NUMBER_OF_POINTS;
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
	}
	else
	{
		int Nelements = TOTAL_DATA_POINTS;
		TGraph2D * ScatterPlot = new TGraph2D(Nelements, array1, array2, array3);
		return ScatterPlot;
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
	int z_zero_index;
    for (int i = 0; i < NZ; i++)
    {
    	if (fabs(z_data[i]) < 1.0e-10)
        z_zero_index = i;
    }
    double Bzero = b_field[z_zero_index];
    double integral = data -> Integral();
    double Leff = integral/Bzero;
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

// main function to read in all files, create a scatter plot and fit
TCanvas * main()
{
/* The main function reads in the data from all txt files of the form 
"EngeProd_*.**.txt" and creates three arrays for current, b_field, and z_pos. 
It creates a scatter plot of the data and fits a function to that plot. It then
prints the fitted parameters, and it saves a plot of the fit to the folder this
code is kept in. */
	TCanvas * ScatterCanvas = // creates a canvas to plot the data and fit
		new TCanvas("ScatterCanvas","Three Dimensional B Field Data",900,600);
	input_data WholeData;	  // creates a space for the data to be read to
	WholeData.Iamps   = (double*) malloc(NUMBER_OF_POINTS * sizeof(double));
	WholeData.z_pos   = (double*) malloc(NUMBER_OF_POINTS * sizeof(double));
	WholeData.b_field = (double*) malloc(NUMBER_OF_POINTS * sizeof(double));
	double FitParameters[NP];
	double FileIndices[] = {0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00};
	double Effective_Lengths[NF];
	int j = 0;
	while (j < NF)
	{
	/* The while loop reads in all the data from the files and puts it into 
	three double arrays. */
		const char * fname = create_name_from_index(FileIndices[j]);
		float curr = FileIndices[j] * I_LIM;
		input_data FileData = read_file_make_input_data(fname);
		
		TGraph * graph = new TGraph(NUMBER_OF_POINTS, 
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
/*	
	double p[NP] =
	{3.95e-1, 3.61e-3, -1.51e-5     			// c0i dependence on I
	,4., 4.18e-4, -1.45e-5   					// c1i dependence on I
    ,-1.34, 9.51e-3, -1.67e-5     				// c2i dependence on I
    ,0.535,-4.01e-3, 1.02e-5      				// c3i dependence on I
    ,0., 0.0, 0.0          						// c4i dependence on I
    ,0., 0.0, 0.0          						// c5i dependence on I
    ,3.95e-1, 3.61e-3, -1.51e-5      			// c0e dependence on I
    ,4., 4.18e-4, -1.45e-5   					// c1e dependence on I
    ,-1.34, 9.51e-3, -1.67e-5     				// c2e dependence on I
    ,0.535,-4.01e-3, 1.02e-5      				// c3e dependence on I
    ,0., 0.0, 0.0          						// c4e dependence on I
    ,0., 0.0, 0.0 								// c5e dependence on I
    ,-0.114, 2.15e-2, -5.26e-5,0.,0.,0.};		// B0  dependence on I
*/
	double p[NP] =
	{0., 0.0, 0.0     				// c0i dependence on I
	,3., 0.0, 0.0					// c1i dependence on I
    ,0., 0.0, 0.0    				// c2i dependence on I
    ,0., 0.0, 0.0      				// c3i dependence on I
    ,0., 0.0, 0.0          			// c4i dependence on I
    ,0., 0.0, 0.0          			// c5i dependence on I
    ,0., 0.0, 0.0      				// c0e dependence on I
    ,4., 0.0, 0.0  					// c1e dependence on I
    ,0., 0.0, 0.0     				// c2e dependence on I
    ,0., 0.0, 0.0     				// c3e dependence on I
    ,0., 0.0, 0.0          			// c4e dependence on I
    ,0., 0.0, 0.0 					// c5e dependence on I
    ,0.3, 0., 0.,0.,0.,0.};			// B0  dependence on I
    
    p[39] = leff_poly -> GetParameter(0);
    p[40] = leff_poly -> GetParameter(1);
    p[41] = leff_poly -> GetParameter(2);

    
	TF2 * func = new TF2("func", func_EngeProd_calc_leff, ZMIN, ZMAX, IMIN, IMAX,NP);
	for (int i = 0; i < NP; ++i)
		func -> SetParameter(i,p[i]);
	for (int i = (NP-3); i < NP; ++i)
		func -> SetParLimits(i,p[i],p[i]);
	for (int i = 0; i < (NP-3); ++i)
		func -> SetParLimits(i,-10,10);
	ScatteredData -> Fit("func");
	for (int i = 0; i < (NP-3); ++i)
		FitParameters[i] = func -> GetParameter(i);
	
	func -> SetTitle("Magnetic Field Global Fit; ; ;Magnetic Field (T)");
	func -> SetLineColor(1);
	func -> SetLineWidth(1);
	func -> Draw("surf1z");
	ScatteredData -> Draw("same");
    
    for (int i = 0; i < NP; ++i)
		cout << FitParameters[i] << endl;
	//ScatterCanvas -> SaveAs("GlobalFieldFit.png");
	return ScatterCanvas;
}