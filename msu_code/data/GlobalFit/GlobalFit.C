// Fit data to an Enge function of z position and current to map B field
// Global variables
double Leff, r0, Iamps;
r0 = 0.2; Leff = 0.7;
int NumberOfPoints = 401;
int TotalDataPoints = NumberOfPoints * 10;
const int Np = 39;
const double e = 2.718281828459045;
double zmax = 1.;
double zmin = -zmax;
double Imin = 0.;
double Imax = 200.;
int Ilim = 180;
struct InputData
{
	double * Iamps;
	double * z_pos;
	double * b_field;
};

// Function to read in files and create a structure of arrays
// Files are sorted based on current, so they must be read in based on that criteria
InputData read_file_make_InputData(const char * FileName, float current)
{
	const int Nline = 256;
	char string_location[Nline];
	InputData NewData;
	FILE * file;
	file = fopen(FileName, "r");
	NewData.Iamps   = (double*) malloc(NumberOfPoints * sizeof(double));
	NewData.z_pos   = (double*) malloc(NumberOfPoints * sizeof(double));
	NewData.b_field = (double*) malloc(NumberOfPoints * sizeof(double));
	fgets(string_location, Nline, file);
	fgets(string_location, Nline, file);
	fgets(string_location, Nline, file);
	fgets(string_location, Nline, file);
	for ( int i = 0; i < NumberOfPoints; ++i)
	{
		fgets(string_location, Nline, file);
		sscanf(string_location, "%lg %lg\n",
			&NewData.z_pos[i], &NewData.b_field[i]);
		NewData.Iamps[i] = current;
	}
	return NewData;
}

// Create a file name given the double index (i.e. 0.10, 0.20, ..., 1.00)
const char * create_name_from_index(double index = 0.10)
{
	int k = 2;
	stringstream stream;
	stream << fixed << setprecision(k) << index;
	string mystring = stream.str();
	
	string front ("EngeProd_");
	string back  (".txt");
	string FileName = front + mystring + back;
	
	const int n = FileName.length();
	const char * fname;
	fname = (char*) malloc((n+1) * sizeof(char));
	strcpy(fname, FileName.c_str());
	
	return fname;
}

// Compare two arrays for equality
int compare_arrays(double * array1, double * array2)
{
	int n = sizeof(array1)/sizeof(double);
	if (sizeof(array1) != size0f(array2))
		return 0;
	else
	{
		for (int i = 0; i < n; ++i)
		{
			if (array1[i] != array2[i])
				return 0;
		}
		return 1;
	}
}

// Function to combine multiple arrays into one larger array
double * combine_arrays(double * array1, double * array2, int k)
{

	double * NewArray;
	NewArray = (double*) malloc(sizeof(double)*(NumberOfPoints*k) + sizeof(double)*NumberOfPoints);
	int Nelements1 = NumberOfPoints*k;
	int Nelements2 = NumberOfPoints;
	for (int i = 0; i < Nelements1; ++i)
		NewArray[i] = array1[i];
	for (int i = Nelements1; i < (Nelements1 + Nelements2); ++i)
		NewArray[i] = array2[i - Nelements1];
	return NewArray;
}

// Plot data from three arrays into a 2D scatter plot
TGraph2D * scatter(double * array1, double * array2, double * array3)
{
	if (sizeof(array1) != sizeof(array2))
		cout << "Function not compatible; arrays have different size." << endl;
	else 
	{
		if (sizeof(array2) != sizeof(array3))
			cout << "Function not compatible; arrays have different size." << endl;
	}
	else
	{
		int Nelements = TotalDataPoints;
		TGraph2D * ScatterPlot = new TGraph2D(Nelements, array1, array2, array3);
		return ScatterPlot;
	}
}

// Create a product of Enge functions to be fitted
double func_EngeProd(double * z, double * p)
{
	double zi =  z[0] - Leff/2.0;
	double ze = -z[0] - Leff/2.0;
	zi /= 2*r0; ze /= 2*r0;
	double I = z[1];
	double fi = 1./(1.+pow(e,((p[0]+p[1]*I+p[2]*I*I) + zi*((p[3]+p[4]*I+p[5]*I*I)
		+ zi*((p[6]+p[7]*I+p[8]*I*I) + zi*((p[9]+p[10]*I+p[11]*I*I)
			+ zi*((p[12]+p[13]*I+p[14]*I*I) + zi*(p[15]+p[16]*I+p[17]*I*I))))))));
	double fe = 1./(1.+pow(e,((p[18]+p[19]*I+p[20]*I*I) + ze*((p[21]+p[22]*I+p[23]*I*I)
		+ ze*((p[24]+p[25]*I+p[26]*I*I) + ze*((p[27]+p[28]*I+p[29]*I*I)
			+ ze*((p[30]+p[31]*I+p[32]*I*I) + ze*(p[33]+p[34]*I+p[35]*I*I))))))));
	double B  = p[36] + p[37]*I + p[38]*I*I;
	double FinalFunction = B*fi*fe;
	return FinalFunction;
}

double 1Dfunc_EngeProd(double * z, double * p)
{
	double zi =  z[0] - Leff/2.0;
	double ze = -z[0] - Leff/2.0;
	zi /= 2*r0; ze /= 2*r0;
	double fi = p[36]/(1.+pow(e,((p[0]+p[1]*p[37]+p[2]*p[37]*p[37]) + zi*((p[3]+p[4]*p[37]+p[5]*p[37]*p[37])
		+ zi*((p[6]+p[7]*p[37]+p[8]*p[37]*p[37]) + zi*((p[9]+p[10]*p[37]+p[11]*p[37]*p[37])
			+ zi*((p[12]+p[13]*p[37]+p[14]*p[37]*p[37]) + zi*(p[15]+p[16]*p[37]+p[17]*p[37]*p[37]))))))));
	double fe = 1./(1.+pow(e,((p[18]+p[19]*p[37]+p[20]*p[37]*p[37]) + ze*((p[21]+p[22]*p[37]+p[23]*p[37]*p[37])
		+ ze*((p[24]+p[25]*p[37]+p[26]*p[37]*p[37]) + ze*((p[27]+p[28]*p[37]+p[29]*p[37]*p[37])
			+ ze*((p[30]+p[31]*p[37]+p[32]*p[37]*p[37]) + ze*(p[33]+p[34]*p[37]+p[35]*p[37]*p[37]))))))));
	double FinalFunction = fi*fe;
	return FinalFunction;
}

// main function to read in all files, create a scatter plot and fit
TCanvas * main()
{
	TCanvas * ScatterCanvas = new TCanvas("ScatterCanvas","Three Dimensional B Field Data",800,600);
	//ScatterCanvas -> Divide(2,7);
	InputData WholeData;
	WholeData.Iamps   = (double*) malloc(NumberOfPoints * sizeof(double));
	WholeData.z_pos   = (double*) malloc(NumberOfPoints * sizeof(double));
	WholeData.b_field = (double*) malloc(NumberOfPoints * sizeof(double));
	double FitParameters[Np-1];
	double FileIndices[] = {0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00};
	int j = 0;
	while (j < 10)
	{
		const char * fname = create_name_from_index(FileIndices[j]);
		float curr = FileIndices[j] * Ilim;
		InputData FileData = read_file_make_InputData(fname, curr);
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
				WholeData.b_field = combine_arrays(WholeData.b_field, FileData.b_field,j);
		}
		j = j + 1;
	}
	
	TGraph2D * ScatteredData = scatter(WholeData.z_pos,WholeData.Iamps,WholeData.b_field);
	
	double p[Np] =
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
    ,-0.114, 2.15e-2, -5.26e-5};				// B0 dependence on I
	
	TF2 * func = new TF2("func", func_EngeProd, zmin, zmax, Imin, Imax,Np);
	for (int i = 0; i < Np; ++i)
		func -> SetParameter(i,p[i]);
	for (int i = 0; i < Np; ++i)
		func -> SetParLimits(i,-10,10);
	ScatteredData -> Fit("func");
	for (int i = 0; i < Np; ++i)
		FitParameters[i] = func -> GetParameter(i);
	
	func -> SetTitle("Magnetic Field Global Fit;Z Position (m);Current (A);Magnetic Field (T)");
	func -> Draw("surf1");
	ScatteredData -> Draw("same p0");
	auto legend = new TLegend(0.1,0.7,0.35,0.82);
    legend->AddEntry("ScatteredData","Input Data Points","l");
    legend->AddEntry("func","Globally Fitted Function","l");
    legend->Draw();
    //ScatterCanvas -> GetXAxis() -> SetTitleOffset(1.5);
    //ScatterCanvas -> GetYAxis() -> SetTitleOffset(1.5);
    //ScatterCanvas -> Update();
    
    for (int i = 0; i < Np; ++i)
		cout << p[i] << endl;
	ScatterCanvas -> SaveAs("I:/projects/fribmappers/tools/fit_test_data_vs_curr/data/EngeCoefFit/GlobalFieldFit.png")
	return ScatterCanvas;
}