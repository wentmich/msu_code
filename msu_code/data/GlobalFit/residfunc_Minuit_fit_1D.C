#include <iostream>
#include <fstream>
#include <cmath>
//#include <ctype>
//#include "TFitter.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TROOT.h"
#include <string.h>
#include "TGraph.h"
#include "TGraph2D.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
#include "TArc.h"
#include "TF1.h"
#include "TBox.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TLegend.h"
using namespace std;

// --------------------------------Define GLobal Variables --------------------
// number of data points for "residfunc"
Long64_t NEntries;
Float_t FitMinValue;
// parameters
Float_t A0, xc0, sx0; // single gaussian
// global data arrays to hold actual field data to fit the theoretical model to
Float_t *dataX,*dataY;

//Define functions to use as theoretical model
Float_t one_gaus(Float_t X){
  Float_t value;
  value = (X-xc0)/sx0 ;
  value = A0*exp(-0.5*value*value)/(2.50662827463100024*sx0);
  return value;
}

//Define the generic rate equation for arrays of parameters
Float_t yeqn( Int_t n, Float_t t){
  //  source[n]   decay[n]    decay source[n-1]
  return R[n] - a[n]*y0[n] + a[n-1]*y0[n-1];
}


// ------------------------------- Evaluation of Fit --------------------------
// Residual function that evaluates how good a parameter set describes 
// the field described in the global arrays dataX, dataY, ...
// Minimizer supplies the variable parameters into the array "par"
// unction takes the constant magnet parameters from the global structure
// then it calculates the theoretical field at each data point defined by the
// based on http://root.cern.ch/root/html/tutorials/fit/NumericalMinimization.C.html
double residfunc(const double *par){
  double sum=0.0, weight, theo, tmp;
  Int_t ParIndex; 
  //--------------------------------------------FIT
  //Define the parameters allowed for fitting
  ParIndex = 0;
  A0 = par[ParIndex++];
  xc0 = par[ParIndex++];
  sx0 = par[ParIndex++];
  //Weighted residual sum between theoretical function and data
  for(Long64_t i=0;i<NEntries;i++){
    theo = one_gaus(dataX[i]);
    weight = 1.0 ; // default
    tmp = theo - dataY[i];
    sum += tmp*tmp*weight;
  }
  return (double)sum;
}


void FitToModel(){
  // definitions for minimizer, should try to move to subroutine
  const char * minName = "Minuit"; //"Minuit", "Minuit2", ...
  const char * algoName = "" ; // "", "Combined",...
  cout << "   minName = " << minName << endl;

  //first create an instance of the minimizer
  ROOT::Math::Minimizer* obj =
  ROOT::Math::Factory::CreateMinimizer(minName, algoName);

   // set tolerance , etc... for the minimizer
   obj->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
   obj->SetMaxIterations(10000);  // for GSL 
   obj->SetTolerance(0.001);
   obj->SetPrintLevel(1);

   // determine number of parameters
   const int parnum=3; //cout << " parnum = " << parnum << endl;
    
   // create function wrapper for minmizer
   ROOT::Math::Functor f(&residfunc,parnum); 
   // define the function to be used by minimizer (using the wrapper)
   obj->SetFunction(f);

   // define arrays with parameters and steps
   double variable[parnum];
   double step[parnum] ;
   double min[parnum] ;
   double max[parnum] ;
   //Initialize parameters to be fitted
   A0=1.;
   xc0=-1.;
   sx0=1.;
   // fill arrays with variables and steps for use by minimizer
   Int_t n =0 ;
   variable[n]=A0;    step[n] =  0.1;  min[n]= -1.0e4; max[n]= +1.0e4; n++;
   variable[n]=xc0;   step[n] =  0.01; min[n]= -1.0e2; max[n]= +1.0e2; n++;
   variable[n]=sx0;   step[n] =  0.01; min[n]=     0.; max[n]= +1.0e2; n++;
   cout << "variables set, n = " << n << endl;
   n=0;
   //obj->SetLimitedVariable(n, name, variable, step, min, max);
   obj->SetLimitedVariable(n, "A0",variable[n],step[n],min[n],max[n]); n++;
   obj->SetLimitedVariable(n,"xc0",variable[n],step[n],min[n],max[n]); n++;
   obj->SetLimitedVariable(n,"sx0",variable[n],step[n],min[n],max[n]); n++;

   // do the minimization
   obj->Minimize(); 
   FitMinValue = obj->MinValue();
   cout << "Fit minimum value found = " << FitMinValue << endl ;
}

void load_data(){
  // allocate arrays to make data readable from "residfunc"
  NEntries = 5;
  Float_t *dataXarray = new Float_t[NEntries];
  Float_t *dataYarray = new Float_t[NEntries];
  // fill arrays
  Int_t i;
  i=0; dataXarray[i] = -2.0; dataYarray[i] = 0.1;
  i=1; dataXarray[i] = -1.5; dataYarray[i] = 0.5;
  i=2; dataXarray[i] = -1.0; dataYarray[i] = 1.0;
  i=3; dataXarray[i] = -0.5; dataYarray[i] = 0.4;
  i=4; dataXarray[i] = -0.0; dataYarray[i] = 0.2;
  // set the global pointers for access from "residfunc"
  dataX = &(dataXarray[0]);
  dataY = &(dataYarray[0]);
}

void CanvPlot1D_data_n_fit(){
  printf("CanvPlot1D_data_n_fit: \n");
 //Determine standard deviation of residuals between data and model
  const Long64_t N=NEntries;
  Float_t stdev=0., minX=dataX[0], maxX=dataX[0];
  for (Long64_t k=0; k<N; k++){
  	stdev += pow( one_gaus(dataX[k]) - dataY[k], 2);
  	if (dataX[k]>maxX) maxX=dataX[k] ;
  	if (dataX[k]<minX) minX=dataX[k] ;
  	//printf("%f\t %f\n", one_gaus(dataX[k]), dataY[k]);
  }
  stdev = sqrt(stdev/(N-1));
  printf(" stdev = %f \n", stdev);
  printf(" minX, maxX = %f, %f \n", minX, maxX);
 //Determine errors at +-stdev (i.e. 68% confidence level)
  Float_t eX[N], eY[N];
  for (Long64_t k=0; k<N; k++){
  	eX[k]=0.;
  	eY[k]=stdev;
  }
 //Evaluate theoretical curve over range of data for plot
  const Int_t Nth=10;
  Float_t delX, extX, Xth[Nth], Yth[Nth];
  extX = 0.1*fabs(maxX-minX) ; //factor of range to extend plot by
  delX = (maxX+extX - (minX-extX))/(Float_t)(Nth-1) ;
  for (Int_t i=0; i<Nth; i++){
  	Xth[i] = minX-extX + delX*(Float_t)i ;
  	Yth[i] = one_gaus(Xth[i]) ;
  	printf("%f\t %f\n", Xth[i], Yth[i]);
  }
  printf(" delX = %f \n", delX);

  TCanvas *c1D = new TCanvas("c1D","plots",  50,50,  700,470);
  //
  TGraph *grth = new TGraph(Nth, Xth, Yth);
  grth->Draw("AC");
      grth->SetTitle("Data and fit; x; y");
      grth->SetLineWidth(2.);
      grth->SetLineColor(kGreen+1);
      grth->SetFillColor(0);//avoids black background in Legend
  //
  TGraphErrors *gr = new TGraphErrors(N, dataX, dataY, eX, eY);
  gr->Draw("*");
      gr->SetMarkerColor(kBlack);
      //gr->SetMarkerStyle(20);
      gr->SetMarkerSize(0.8);
      gr->SetFillStyle(0);
      gr->SetFillColor(0);//avoids black background on symbol
      //error bar lines
      gr->SetLineWidth(2.);
      gr->SetLineColor(kGray+3);
      //gr->SetFillColor(0);//avoids black background in Legend
  TLegend *leg = new TLegend(0.8,0.85, 0.95,0.99) ; //x1,y1, x2,y2
      leg->AddEntry(gr,"data","P");
      leg->AddEntry(grth,"model","L");
      leg->Draw();
  c1D->SaveAs("residfunc_Minuit_fit_1D.gif");
  return;
}

// main function:  $root -l residfunc_Minuit_fit_1D.C+()
int residfunc_Minuit_fit_1D(){
  load_data();
  FitToModel();
  CanvPlot1D_data_n_fit();
  cout << "...end."<<endl;
  return 0;
}
