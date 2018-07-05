// Read in coefficient data and plot expected vs. fitted behavior
//Global Structure
struct coefficients
{
	double * expected; double * fitted;
};
const int Np = 39;
int StartPoints[12] = {0,3,6,9,12,15,18,21,24,27,30,33};
int Ks[12] = {0,1,2,3,4,5,6,7,8,9,10,11};

const char * create_name_from_index(int index = 0)
{
	int k = 2;
	stringstream stream;
	stream << index;
	string mystring = stream.str();
	
	string front ("c");
	string back  ("_ExpVsGlobalFit.png");
	string FileName = front + mystring + back;
	
	const int n = FileName.length();
	const char * fname;
	fname = (char*) malloc((n+1) * sizeof(char));
	strcpy(fname, FileName.c_str());
	
	return fname;
}

TCanvas * make_plot()
{	
	TCanvas * c1 = new TCanvas("c1","polynomialgraph",800,600);
	for (int j = 0; j < 12; ++j)
	{
		int StartPoint = StartPoints[j];
		int k = j;
	//c1 -> DrawFrame(0,0,180,4);
	const int Nline = 256;
	char string_location[Nline];
	coefficients Data;
	Data.expected = (double*) malloc(Np*sizeof(double));
	Data.fitted   = (double*) malloc(Np*sizeof(double));
	FILE * file;
	file = fopen("GlobalFit_ExpectedData.txt","r");
	for (int i = 0; i < Np; ++i)
	{
		fgets(string_location, Nline, file);
		sscanf(string_location, "%lg %lg \n", &Data.fitted[i], &Data.expected[i]);
	}
	// Now we have two arrays with all the stuff in them
	TF1 * ExpectedPoly = new TF1("ExpectedPoly","[0]+[1]*x+[2]*x*x",0,180);
	ExpectedPoly -> SetParameters(Data.expected[StartPoint],Data.expected[StartPoint+1],Data.expected[StartPoint+2]);
	TF1 * FittedPoly   = new TF1("FittedPoly","[0]+[1]*x+[2]*x*x",0,180);
	FittedPoly   -> SetParameters(Data.fitted[StartPoint],Data.fitted[StartPoint+1],Data.fitted[StartPoint+2]);
	ExpectedPoly -> SetTitle("Expected Polynomial vs. Globally Fitted Polynomial;Current (A);Enge Coefficient");
	ExpectedPoly -> SetLineColor(1);
	ExpectedPoly -> SetLineStyle(1);
	FittedPoly   -> SetLineColor(4);
	FittedPoly   -> SetLineStyle(2);
	FittedPoly   -> Draw();
	ExpectedPoly -> Draw("same");
	
	/*
	double ExpMax = ExpectedPoly -> GetMaximum(0,180);
	double ExpMin = ExpectedPoly -> GetMinimum(0,180);
	double FitMax = FittedPoly   -> GetMaximum(0,180);
	double FitMin = FittedPoly   -> GetMinimum(0,180);
	if (ExpMax > FitMax)
	{
		if (ExpMin < FitMin)
			FittedPoly -> GetYaxis() -> SetRangeUser(ExpMin,ExpMax);
		else
			FittedPoly -> GetYaxis() -> SetRangeUser(FitMin,ExpMax);
	}
	else 
	{
		if (ExpMin < FitMin)
			FittedPoly -> GetYaxis() -> SetRangeUser(ExpMin,FitMax);
		else
			FittedPoly -> GetYaxis() -> SetRangeUser(FitMin,FitMax);
	}
	c1 -> Modified();
	c1 -> Update();
	*/
	auto legend = new TLegend(0.1,0.7,0.35,0.83);
    legend->AddEntry("ExpectedPoly","Expected","l");
    legend->AddEntry("FittedPoly"  ,"Fitted"  ,"l");
    legend->Draw();
    
    const char* name = create_name_from_index(k);
    c1 -> SaveAs(name);
    //cout << j << '\t' << ExpMax << endl;
    }
}