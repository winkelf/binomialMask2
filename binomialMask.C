#include "helperFunctions.C"
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include <cmath>
#include "TCanvas.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "Math/ProbFunc.h"
#include <fitsio.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <TH2F.h>
#include <TCanvas.h>
#include <nlohmann/json.hpp>
#include <chrono>

using namespace chrono;
using namespace std;


vector<float> analyze(float threshold, string fits_filename_, string fileName_, int ohdu) {

  int   n = 60;
  float r = n/2;
  
  ////////////////////////////////////////////////////////////////////
  ////////////////////// Get previous masks (only for data) //////////
  ////////////////////////////////////////////////////////////////////

  vector<int> maskList = {
	                  1,     // Neighbor
			  4,     // Bleed
			  8,     // Halo
			  16,    // Crosstalk
			  64,    // Edge
			  128,   // Serial register
			  512,   // Hot pixel
			  1024,  // Hot column
			  4096,  // Extended bleed
			  //8192,  // LEC
			  //8194,  // LEC + event?
			  16384, // Full well mask
			  32768, // Cluster shape
                         };

  // Read the FITS file and load image data
  const char* fits_filename = "";
  const char* fileName      = "";

  fits_filename = fits_filename_.c_str();
  fileName      = fileName_.c_str();

  vector<vector<double>> image_data;
  vector<vector<double>> previousMask;

  try {
    previousMask = read_fits_image(fits_filename, image_data, ohdu);
  } catch (const runtime_error& e) {
    cerr << e.what() << endl;
  }

  // Get image from rootfile's TTree and fill a matrix
  vector<vector<int>> matrix_data = loadImage(fileName, ohdu);

  //plotMatrices(matrix_data, previousMask, maskList, true);
  //exit(1);

  ////////////////////////////////////////////////////////////////////
  ////////////////////// Scan grid ///////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  int Nx = 3200;
  int Ny = 520; 

  for (int xcoord = 0; xcoord <= Nx - n; xcoord = xcoord + n/2) {	
    for (int ycoord = 0; ycoord <= Ny - n; ycoord = ycoord + n/2) {	    
  
      double masked  = 0;
      int nevents    = 0;
      
      int n1e_events = 0;
      int n2e_events = 0;

      // Coordinates of the circle center
      int x_center = xcoord + r;
      int y_center = ycoord + r;

      int trials = 0;

      // Count hits in sliding window
      for (int i = ycoord; i < ycoord + n; ++i) {
        for (int j = xcoord; j < xcoord + n; ++j) {

          // Check if point is inside circle		
          bool isInside = isInsideCircle(r, x_center, y_center, j, i); 
          if (!isInside) continue;	

          // Get mask value in pixel, if any 
          int maskInPixel = int(previousMask.at(i).at(j));
	  if (maskInPixel==8194) continue;

          int isMasked = 0;

          for (int m=0; m<maskList.size(); m++) {
            if ((maskInPixel&maskList.at(m))==maskList.at(m) ) {
              isMasked += 1;
            }
          }

          if (isMasked!=0){ continue; }

          // Count number of binomial trials and surviving events 
          trials+=1;
	  if (matrix_data.at(i).at(j)==1)  { n1e_events += 1; }
	  if (matrix_data.at(i).at(j)==2)  { n2e_events += 1; }
        }
      }

      // If a window has 0 trials, is because is already fully masked. Skip it.
      if (trials==0) continue;
      if (n1e_events==0) continue;

      // This number should be the same for all the windows except if a mask is being applied
      int binomN = trials;
      
      long double P = static_cast<float>(n1e_events) / static_cast<float>(binomN);

      int effCheck = binomN - 1;

      long double P_eff = static_cast<float>(n1e_events) / static_cast<float>(effCheck);
   
      if (P>threshold){

        // Set mask in mask matrix
        for (int i = ycoord; i < ycoord + n; ++i) {
          for (int j = xcoord; j < xcoord + n; ++j) {
            // Check if point is inside circle          
            bool isInside = isInsideCircle(r, x_center, y_center, j, i);
            if (!isInside) continue;
            previousMask.at(i).at(j) = 65536;
          }
        }

        //plotMatrices(matrix_data, previousMask, maskList, true);
        //exit(1);

      }//pval
    }// ycoord
  }// xcoord

  maskList.push_back(65536);

  //plotMatrices(matrix_data, previousMask, maskList, true);
  //exit(1);

  // count unmasked pixels
  int nUnmaskedPixels = countUnmaskedPixels(previousMask, maskList, true);

  /////////////////////////////////////////////////////////////////////
  ///////// Counting n-electron events that survived mask /////////////
  /////////////////////////////////////////////////////////////////////

  //cout << "  Number of unmasked pixels: "<< nUnmaskedPixels << endl;

  int n1eEventsUnmasked = 0;
  int n2eEventsUnmasked = 0;

  for (int x = 0; x < 520; ++x) {
    for (int y = 0; y < 3200; ++y) {
      int event       = matrix_data.at(x).at(y);
      int maskInPixel = previousMask.at(x).at(y);
      if (maskInPixel==8194) continue;

      int isMasked = 0;

      for (int m=0; m<maskList.size(); m++) {
        if ((maskInPixel&maskList.at(m))==maskList.at(m) ) {
          isMasked += 1;
        }
      }

      if (isMasked!=0){ continue; }

      if (event==1) n1eEventsUnmasked+=1;
      if (event==2) n2eEventsUnmasked+=1;
    }
  }

  float fl_nUnmaskedPixels   = static_cast<float>( nUnmaskedPixels   );
  float fl_n1eEventsUnmasked = static_cast<float>( n1eEventsUnmasked );
  float fl_n2eEventsUnmasked = static_cast<float>( n2eEventsUnmasked );
  
  vector<float> rVector = { fl_n1eEventsUnmasked, fl_n2eEventsUnmasked, fl_nUnmaskedPixels};

  return rVector;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////                Main           /////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int binomialMask(){

  // Make list of thresholds/cut values to loop over	
  vector<float> thresholds = {};
  for (float i=0; i<0.02; i+=0.0005) thresholds.push_back(static_cast<float>(i));

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetCanvasColor(kBlack);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetPadColor(kBlack);
  gStyle->SetStatColor(kBlack);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(kBlack);
  gStyle->SetTitleTextColor(kWhite);
  gStyle->SetLabelColor(kWhite, "XYZ");
  gStyle->SetTitleColor(kWhite, "XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetHistFillColor(kBlack);
  gStyle->SetHistLineColor(kWhite);
  gStyle->SetFuncColor(kWhite);
  gStyle->SetFrameLineColor(kWhite);
  gStyle->SetTitleFontSize(0.025);

  TFile *f = new TFile("mc_q_distribution_gauss.root", "OPEN");

  // Histo
  TH1D* h2 = new TH1D("","Distribution of q = -2ln #left[ #frac{B(k_{1},p_{1},N) #times B(k_{2},p_{2},N)}{B(#hat{k_{1}},p_{1},N) #times B(#hat{k_{2}},p_{2},N)} #right]",50,0,.1);

  // Open the config file
  std::ifstream configFile("../git/configs/config_multiImage.json");

  if (!configFile.is_open()) {
    cerr << "Could not open config file." << std::endl;
    exit(1);
  }

  string fits_filename_;
  string fileName_;

  const char* fits_filename = "";
  const char* fileName      = "";

  // Parse the JSON
  nlohmann::json config;
 
  // Clock 
  auto a = high_resolution_clock::now();

  vector<float> n1ElectronEventsVec = {};
  vector<float> n2ElectronEventsVec = {};
  vector<float> unmaskedPixelsVec   = {};
  vector<float> oneElectronRateVec  = {};
  vector<float> twoElectronRateVec  = {};

  vector<int> ohdus = {
	                1,
	                2,
	                3,
	                4	
                      };
  
  try {
    configFile >> config;
    for (auto ii : thresholds) {
      cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
      cout<<"Threshold = "<<ii <<endl;
  
      int nImages = 0;
      int nImagesInConfig = config["filePairs"].size();

      float n1ElectronEvents = 0;
      float n2ElectronEvents = 0;
      float unmaskedPixels   = 0;


      for (const auto& pair : config["filePairs"]) {

        string fits_filename_ = pair["fitsFile"];
        string fileName_      = pair["rootFile"];

        fits_filename         = fits_filename_.c_str();
        fileName              = fileName_.c_str();
     
        string fileNameNoPath = fileName_.substr(fileName_.find_last_of("/\\") + 1);
        
	nImages+=1;

        //if (nImages<90) continue;

	cout<<"                                                    "<<endl;
	cout<<"  READING FILE: "<< nImages << "/"<< nImagesInConfig <<endl;
	cout<<"                "<< fileNameNoPath                   <<endl;

	for (int oh=0; oh<ohdus.size(); oh++){     

          int ohdu = ohdus.at(oh);

	  cout<<"                   *ohdu "<< ohdu << endl;

          vector<float> qs  = analyze(ii, fits_filename, fileName, ohdu);
	  
	  cout<<"                      - Number of 1 electron events "<< qs.at(0) << endl;
	  cout<<"                      - Number of 2 electron events "<< qs.at(1) << endl;
	  cout<<"                      - Number of unmasked pixels   "<< qs.at(2) << endl;

          n1ElectronEvents += qs.at(0);
          n2ElectronEvents += qs.at(1);
          unmaskedPixels   += qs.at(2);

        }
        if (nImages>120) break;
      } 

      cout<<"                      - Number of 1 electron events "<< n1ElectronEvents << endl;
      cout<<"                      - Number of 2 electron events "<< n2ElectronEvents << endl;
      cout<<"                      - Number of unmasked pixels   "<< unmaskedPixels   << endl;
      //exit(1);

      n1ElectronEventsVec.push_back( n1ElectronEvents                );
      n2ElectronEventsVec.push_back( n2ElectronEvents                );
      unmaskedPixelsVec  .push_back( unmaskedPixels                  );
      oneElectronRateVec .push_back( n1ElectronEvents/unmaskedPixels  );
      twoElectronRateVec .push_back( n2ElectronEvents/unmaskedPixels  );
    }
  } catch (const std::exception& e) {
    cerr << "Error parsing config file: " << e.what() << endl;
    exit(1);
  }


  // output txt to plot with python
  std::ofstream outFile0("more_survElec.txt");
  std::ofstream outFile1("more_survPix.txt");
  std::ofstream outFile2("more_rates.txt");

  if (!outFile0) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile1) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile2) { cerr << "Error opening file for writing!" << endl; }


  for (size_t t = 0; t < thresholds.size(); ++t) {
    outFile0 << thresholds.at(t) <<"\t" << n1ElectronEventsVec.at(t) << "\t" << n2ElectronEventsVec.at(t) <<endl;
    outFile1 << thresholds.at(t) <<"\t" << unmaskedPixelsVec.  at(t) << "\t" << unmaskedPixelsVec.  at(t) <<endl;
    outFile2 << thresholds.at(t) <<"\t" << oneElectronRateVec. at(t) << "\t" << twoElectronRateVec. at(t) <<endl;
  }

  outFile0.close();
  outFile1.close();
  outFile2.close();

  // chrono
  auto b = high_resolution_clock::now();
  cout << "took " << duration_cast<seconds>(b - a).count() << " seconds" <<  endl;

  return 0;
}
