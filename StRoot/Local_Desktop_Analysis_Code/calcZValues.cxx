/////ROOT LIBRARIES///////////////
#include "TROOT.h"
#include "TSystem.h" 
#include "TApplication.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <iomanip>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TPad.h>
#include <TMath.h>
#include <TF1.h>
#include <THnSparse.h>
#include <TLeaf.h>
#include <TLatex.h>
#include "TError.h"
gErrorIgnoreLevel = 4000;  //default is kInfo

///////C++ Libraries/////////////////
#include <iostream>
#include <iomanip>
#include <math.h>
#include <map>
#include <vector>
#include <utility>
#include <climits>
#include <sstream>
#include <string>
#include <fstream>
#include <TAttMarker.h>
#include <algorithm>
#include <memory>
#include <iterator>
#include <ctype.h>
#include <bitset>
#include <stdio.h>
#include <stdlib.h>



using namespace std;




int calcZValues(){

    //starting values
    int NEtaBins = 9;
    int NPhiBins = 24;
    
    
    int NAbsEtaBins;
    int NAbsPhiBins;
    int checkSum;
    
    NAbsEtaBins = ((NEtaBins-1)/2) + 1;
    NAbsPhiBins = (NPhiBins/4) + 1;
    
    cout << "N_Abs_Eta_bins = " << NAbsEtaBins << "    " << "N_Abs_Phi_bins = " << NAbsPhiBins << endl << endl;

    checkSum = (NEtaBins*NPhiBins)*(NEtaBins*NPhiBins);
    
    cout << "checkSum:   " << checkSum << endl << endl;
    
    int maxIndex = (NEtaBins*NPhiBins);
    
    double sum = 0;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    }
