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

#include "C:/Users/ajentsch/desktop/corrHistogramMaker.h"
//void formatCorrHist(TH2D* hist);
//void formatCorrHist(TH2D* hist, TString title);

Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if ( (reject && (x[0] > 1.81 && x[0] < 1.91)) || (reject && (x[0] > 1.6 && x[0] < 1.7)) ){     //|| (reject && (x[0] > 2.0 && x[0] < 2.1))
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}

double calcScaleFactorforFit(double sigma)
{
    double preFactor = TMath::Sqrt(TMath::PiOver2());
    
    return preFactor*sigma*.5;
    
}
        

int corrHistogramMaker(){
  

    

//-------------------------------------------------------------
//BEGIN VARIABLES SECTION
//-------------------------------------------------------------  
    int NUM_PHI_BINS = 12;
    int NUM_ETA_BINS = 9;
    
    int NUM_CENT_BINS = 16;
    int NUM_VZ_BINS   = 10;
    int NUM_PT_BINS   = 4;
    
    double NUM_FILES_COMBINED = 1;
    
    double ETA_RANGE = 1.55;
    
    TString numPhi = "12";
    TString numEta = "9";
    
    TString inputRootFolder = "D0_Hadron_Correlation_Root_Files/";    
    TString subFolder1 = "corrHistograms_";
    TString subFolder2 = "_Eta_";
    TString subFolder3 = "_Phi_bins";
    TString rootFile   = ".root";
    TString outputRoot = "corrHistogramMaker";
    TString slash      = "/";
    //-----------------------
    //file naming and input/output
    //-----------------------
    
    bool PRINT_SUB_HISTOS_RAW = false;
    bool PRINT_SUB_HISTOS_SCALED = false;
    bool PRINT_SUB_HISTOS_DEL_RHO = false;
    bool PRINT_DELRHO_OVER_REF_HISTOS = false;
    
    
    TString mainPath       = "C:/Users/ajentsch/Desktop/";    
    TString mainFolder     = "D0 Correlation Output Histograms/";
    TString binSubFolder   = subFolder1 + numEta + subFolder2 + numPhi + subFolder3;
    
    TString inputFileName  = mainPath + inputRootFolder + numEta + subFolder2 + numPhi + subFolder3 + rootFile;
     
    
    TString path           = mainPath + mainFolder + binSubFolder + slash;
    TString outputFileName = path + outputRoot + rootFile;
    
    TString fileType    = ".png"; //file type for output histogram pictures
    TString fileTypeEps = ".eps"; //file type for progess report
    
    TString binLabelVz[10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    TString binLabelCent[16] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
    TString binLabelPt[4]    = {"0", "1", "2", ""};
    
    TString PtBinLabel[4]   = {"_PtBin_", "_PtBin_", "_PtBin_", ""};
    
    TString outputFolders[9] = {"SibCorr/", "MixedCorr/", "ScaledMixedCorr/", "SibMinusScaledMixCorr/", 
                                 "delRho_over_rhoRefCorr/", "FullySubtractedCorr/", "USMinusLS/", 
                                 "SideBandSubtraction/", "Centralities/"};
                                
    TString invariantMassOutputFolders[4] = {"PtBin0/", "PtBin1/", "PtBin2/", "PtIntegrated/"};
    TString invariantMassHistogramsFolder = "invariantMassHistograms/";
        
                                
    TString VzSubFolders[11] = {"VzBin0/", "VzBin1/", "VzBin2/", "VzBin3/", "VzBin4/", "VzBin5/", "VzBin6/", "VzBin7/", "VzBin8/", "VzBin9/", "VzIntegrated/"};
    TString PtSubFolders[4] = {"PtBin0/", "PtBin1/", "PtBin2/", "PtIntegrated/"};
    TString bandFolders[5]   = {"LeftSideBand/","UnlikeSign/", "RightSideBand/", "LikeSign/", "SideBandAverage/"};    
    
    TFile *file = new TFile(inputFileName); //Root file with raw histograms to use for calculations
    
    TFile *output = new TFile(outputFileName, "RECREATE"); //root file to store calculated and scaled histograms
    
    
   
    TH1D*   oneDUSHistos[4][4];     
    TString oneDUSStrings[4] = {"D0_US_invMass_Pt_Bin_0", "D0_US_invMass_Pt_Bin_1", "D0_US_invMass_Pt_Bin_2", "unlikeSign"};
                                              
    TH1D*   oneDLSHistos[4][4];     
    TString oneDLSStrings[4] = {"LS_invMass_Pt_Bin_0", "LS_invMass_Pt_Bin_1", "LS_invMass_Pt_Bin_2", "LikeSignBG"};
                                              
    TH1D*   oneDScaledLSHistos[4][4];     
    TString oneDScaledLSStrings[4] = {"ScaledLS_PtBin_0", "ScaledLS_PtBin_1", "ScaledLS_PtBin_2", "ScaledLS"};
                                                  
    TH1D*   oneDSubtractedInvMassHistos[4][4];                                              
    TString oneDSubtractedInvMassStrings[4] = {"D0_Minus_Scaled_LS_BG_PtBin_0", "D0_Minus_Scaled_LS_BG_PtBin_1", "D0_Minus_Scaled_LS_BG_PtBin_2", "D0_Minus_Scaled_LS_BG"}; 
    
    TH1D*   oneDSubtractedInvMassHistosFunction[4][4];                                              
    TString oneDSubtractedInvMassStringsFunction[4] = {"US_Minus_Expo_Fit_BG_PtBin_0", "US_Minus_Expo_Fit_BG_PtBin_1", "US_Minus_Expo_Fit_BG_PtBin_2", "US_Minus_Expo_Fit_BG"}; 
                                                                                                                          
    TString binLabelCentralityClass[3] = {"_Peripheral", "_MidCentral", "_Central"};
    
    
    double integralSideBandCounts[2];
    double integralCounts[2];
   
    //-----------------------------------
    //Correlation variables and arrays
    //-----------------------------------
    
    // US histograms
    TH2D* sibCorrBin[4][4][10][16];    //Raw US sibling histograms binned by centrality
    TH2D* mixCorrBin[4][4][10][16];    //Raw US mixed histograms binned by centrality
    TH2D* scaledMixCorrBin[4][4][10][16];  //storage for the eventual scaled mixed US histograms  
    TH2D* sibMinusScaledMix[4][4][10][16]; // Raw US Sibling minus scaled mixed                   
    TH2D* delRhoOverRhoRefBin[4][4][10][16];
    
    //ERRORS STORED HERE--ALL ERRORS STORED AS VARIANCES -- SQUARE ROOT MUST BE TAKEN FOR ACTUAL ERROR
    
    double sibCorrBinStatErrors[4][4][10][16][9][12];    
    double mixCorrBinStatErrors[4][4][10][16][9][12];
    double scaledMixCorrBinStatErrors[4][4][10][16][9][12];
    double sibMinusScaledMixStatErrors[4][4][10][16][9][12];
    double delRhoOverRhoRefBinStatErrors[4][4][10][16][9][12];
    //------------------------------------------------------------------
    
    TH2D* sibCorrBinVzInt[4][4][16];
    TH2D* mixCorrBinVzInt[4][4][16];
    TH2D* scaledMixCorrBinVzInt[4][4][16];  
    TH2D* sibMinusScaledMixVzInt[4][4][16]; 
    TH2D* delRhoOverRhoRefBinVzInt[4][4][16];
    
    TH1D* sibCorrBinPhiProj[4][4][10][16];    //Raw US sibling histograms binned by centrality
    TH1D* mixCorrBinPhiProj[4][4][10][16];    //Raw US mixed histograms binned by centrality
    TH1D* scaledMixCorrBinPhiProj[4][4][10][16];  //storage for the eventual scaled mixed US histograms
    TH1D* sibMinusScaledMixPhiProj[4][4][10][16]; // Raw US Sibling minus scaled mixed
    TH1D* delRhoOverRhoRefBinPhiProj[4][4][10][16];
    
    TH1D* sibCorrBinPhiProjVzInt[4][4][16];
    TH1D* mixCorrBinPhiProjVzInt[4][4][16];
    TH1D* scaledMixCorrBinPhiProjVzInt[4][4][16];  
    TH1D* sibMinusScaledMixPhiProjVzInt[4][4][16]; 
    TH1D* delRhoOverRhoRefBinPhiProjVzInt[4][4][16];
  
    
    TH1D* sibCorrUSPhiProjMinusLSPhiProjVzInt[3][16];
   
   
    TH2D* sideBandAverage[4][16];
    TH1D* sideBandAveragePhiProj[4][16];    
    TH2D* fullySubtractedCorrSideBandCent[4][16];
    TH1D* fullySubtractedCorrSideBandCentPhiProj[4][16];
    
    TH2D* delRhoOverRhoRefBinVzIntCentInt[4][4][3];
    TH1D* delRhoOverRhoRefBinVzIntCentIntPhiProj[4][4][3];
    
    TH2D* fullySubtractedCorrLSCent[4][16];
    TH1D* fullySubtractedCorrLSCentPhiProj[4][16];
    
    //TH2D* mixCorrUSVzInt[3][16];
    //TH2D* mixCorrLSVzInt[3][16];
    
    //Correlations with integrated centralities
    TH2D* fullySubtractedCorrCustomCentralityUS[4][3];
    TH2D* fullySubtractedCorrCustomCentralityLS[4][3];
    TH1D* fullySubtractedCorrCustomCentralityLSPhiProj[4][3];
    TH2D* fullySubtractedCorrCustomCentralitySideBand[4][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandPhiProj[4][3];
    double VzWeightFactorCustom[3][10];
    double VzWeightFactorCustomTotal;
    
    TString fullySubtractedLabelCustomCentrality = "FullSubtractedCorr_";
    TString centralityBin[3] = {"Peripheral", "MidCentral", "Central"};
    
    TString SibCorrLabels[4]          = {"Sibling_SideBandLeft_correlation", "Sibling_US_correlation", "Sibling_SideBandRight_correlation", "Sibling_LS_correlation"};
    TString MixCorrLabels[4]          = {"Mixed_SideBandLeft_correlation", "Mixed_US_correlation", "Mixed_SideBandRight_correlation", "Mixed_LS_correlation"};
    TString ScaledMixCorrLabels[4]     = {"Scaled_Mixed_SideBandLeft_correlation", "Scaled_Mixed_US_correlation", "Scaled_Mixed_SideBandRight_correlation", "Scaled_Mixed_LS_correlation"};
    TString SibMinusScaledMixLabels[4] = {"Sib_Minus_Scaled_Mix_SideBandLeft_correlation", "Sib_Minus_Scaled_Mix_US_correlation", "Sib_Minus_Scaled_Mix_SideBandRight_correlation", "Sib_Minus_Scaled_Mix_LS_correlation"};
    TString delRhoOverRhoRefLabels[4]  = {"delRho_over_RhoRef_SideBandLeft_correlation", "delRho_over_RhoRef_US_correlation", "delRho_over_RhoRef_SideBandRight_correlation", "delRho_over_RhoRef_LS_correlation"};
   
    TString subtractedCorrLabels[2] = {"US_Minus_SidebandAverage_Corr", "US_Minus_LS_Corr"};
    
    TString sideBandAverageLabel = "sideBandAverage";
    
    TString phiProj = "Phi_projection_";
    
    TString USMinusLSSubtractedLabel1 = "US_Minus_LS_Subtracted_CentBin_";
    
    TString withFitLabel = "_With_Fit";
    
    //TString delRhoOverRhoRefUSMinusLSLabel1 = "delRho_over_RhoRef_US_Minus_Scaled_LS_";
    
    TString centBinLabel = "_CentBin_";
    TString VzBinLabel = "_VzBin_";
    TString VzIntLabel = "_Vz_Integrated_";
    TString bandLabel[4] = {"SideBandLeft", "US", "SideBandRight", "LS"};

    double numSib = 0;
    double numMix = 0;
    double numMixUS = 0;
    double numMixLS = 0;
    
    double integralRawMixHistos[4][4][10][16];
    double integralRawMixHistosVzInt[4][4][16];
    //double integralRawMixHistos[4][4][10][16];
    double integralRawMixHistosVzIntCentInt[4][4][3];
   
    double scaleFactorVzCent = 0;
    
    double scaleFactorVz = 0;
    double scaleFactorLS = 0;
    
    
    double LSScaleFactor[4][4] = {{1.0, 1.0, 1.0, 1.0}, {1.0,1.0,1.0,1.0}};
    double ScalingRatio[3][16];
    
    double BOverSPlusBLS[4][4];
    double BOverSPlusBFit[4][4];
    double SPlusBOverSLS[4][4];
    double SPlusBOverSFit[4][4];
    
    double totalPairsSibling = 0;
    double integralError[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    //fully subtraced US minus LS
    
    TString delRhoOverRhoRefUSMinusLSLabel1 = "US_delRhoOverRho_Minus_LS_delRhoOverRho_";
    
    //----------------------------
    //General use variables
    //----------------------------
   
    TString str1;
    TString str2;
    TString str3;
    TString str4;
    

//-------------------------------------------------------------
//END VARIABLES SECTION
//-------------------------------------------------------------    
    
    TCanvas *c = new TCanvas("c2", "Histograms", 1100, 850);

    double countsPeakRegionUS = 0;
    double countsPeakRegionLS = 0;
    double countsSideBandUS   = 0;
    double countsSideBandLS   = 0;
    double sideBandLow      = 2.0;
    double sideBandHigh     = 2.1;
    double massLow          = 1.82;
    double massHigh         = 1.90;
    
    double integralFit      = 0;
    
    double integral         = 0;
    //double LSScaleFactor    = 0;
    
    Int_t binMassLowUS;
    Int_t binMassHighUS;
    Int_t binMassLowLS;
    Int_t binMassHighLS;
        
    TString title;
    TString tmp;
        
    Double_t par[5];
    Double_t parBG[2];
        
    TF1 * g1 = new TF1("m1", "gaus", 1.8, 1.9);
    TF1 * e1 = new TF1("m2", "expo", 1.6, 1.78);
    
    //TF1 * eL = new TF1("m3", "expo", 1.69, 1.81);
    //TF1 * eR = new TF1("m4", "expo", 1.92, 2.0);
    
    TF1 * BGFunc = new TF1("BGfunc", "expo", 1.6, 2.1);
    
    TF1 * fun = new TF1("fun","gaus(0)+expo(3)",1.6,2.1);
    
    TF1 * finalGaussian = new TF1("final Guassian", "gaus", 1.82, 1.9);
    
    double par0;
    double par1;
    double par2;
    double par3;

    TAxis* xAxisUS;
    TAxis* xAxisLS;
    
    int LabelVersion = 0;
    
//------------------------------------------------------------------
//Code section to write information about this dataset to file
//------------------------------------------------------------------    
    /* histOfCuts->SetBinContent(1, D0InvMassLow); //D0InvMassLow
    histOfCuts->SetBinContent(2, D0InvMassHigh); //D0InvMassHigh
    histOfCuts->SetBinContent(3, USSideBandLeftLow); //USSideBandLeftLow
    histOfCuts->SetBinContent(4, USSideBandLeftHigh); //USSideBandLeftHigh
    histOfCuts->SetBinContent(5, USSideBandRightLow); //USSideBandRightLow
    histOfCuts->SetBinContent(6, USSideBandRightHigh); //USSideBandRightHigh
    histOfCuts->SetBinContent(7, d0PtLow); //d0PtLow 
    histOfCuts->SetBinContent(8, d0PtHigh); //d0PtHigh
    histOfCuts->SetBinContent(9, d0DecayLengthMin); //d0DecayLengthMin
    histOfCuts->SetBinContent(10, d0DecayLengthMax); //d0DecayLengthMax
    histOfCuts->SetBinContent(11, daughterDCA); //daughterDCA
    histOfCuts->SetBinContent(12, d0DaughterPionPtMin); //d0DaughterPionPtMin
    histOfCuts->SetBinContent(13, d0DaughterKaonPtMin); //d0DaughterKaonPtMin
    histOfCuts->SetBinContent(14, kaonDCA); //kaonDCA
    histOfCuts->SetBinContent(15, pionDCA); //pionDCA
    histOfCuts->SetBinContent(16, d0DCAtoPV); //d0DCAtoPV
    histOfCuts->SetBinContent(17, hadronPt); //hadronPt
    histOfCuts->SetBinContent(18, BUFFER_SIZE); //BUFFER_SIZE
    histOfCuts->SetBinContent(19, NUM_PHI_BINS); //NUM_PHI_BINS
    histOfCuts->SetBinContent(20, NUM_ETA_BINS); //NUM_ETA_BINS
    
    histOfCuts->SetBinContent(17, hadronPtMin); //hadronPtMin
    histOfCuts->SetBinContent(18, hadronPtMax); //hadronPtMax
    histOfCuts->SetBinContent(19, BUFFER_SIZE); //BUFFER_SIZE
    histOfCuts->SetBinContent(20, NUM_PHI_BINS); //NUM_PHI_BINS
    histOfCuts->SetBinContent(21, NUM_ETA_BINS); //NUM_ETA_BINS
    histOfCuts->SetBinContent(22, 1); // Require HFT
    histOfCuts->SetBinContent(23, trackChi2Max); // chi2 cut
    histOfCuts->SetBinContent(24, trackDCAtoPvtx); // DCA cut
    histOfCuts->SetBinContent(29, 1); //Scale factor counter to produce text file*/
    
    
    double numEvents = 0;
    
    numEvents = (((TH1I*) file->Get("number of events used"))->GetBinContent(2))/(1000000);
    
    ofstream cutFile;
    TString cutFileName = "Cuts_and_data_information.txt";
    TString cutOutputFile = path + cutFileName;
    cutFile.open (cutOutputFile);
    
    cutFile << "Important information about this data run" << endl << endl;
    cutFile << "Number of Events: " << numEvents << "M Events" << endl << endl;
    
    NUM_FILES_COMBINED = ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(29);
    
    ((TH1D*) file->Get("HistOfCuts"))->Scale(1/NUM_FILES_COMBINED);
    
    cutFile << "Trigger D0 Cuts" << endl << endl;
    cutFile << "D0 InvMass Signal Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(1) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(2) << endl;
    cutFile << "D0 SideBandLeft   Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(3) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(4) << endl;
    cutFile << "D0 SideBandRight  Band -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(5) << "\t\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(6) << endl;
    cutFile << "D0 Pt Cuts             -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(7) << "\t\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(8) << endl;
    cutFile << "D0 DecayLengthCut      -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(9) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(10) << endl;
    cutFile << "D0 DaughterPtCut       -\t" << "Pi : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(12) << "\t" << "K  : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(13) << endl;
    cutFile << "D0 K/Pi DCA to PV      -\t" << "Pi : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(15) << "\t" << "K  : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(14) << endl;
    cutFile << "DaughterDCA            -\t" << "Less Than : " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(11) << endl;
    cutFile << "D0 DCA to PV           -\t" << "Greater Than: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(16) << endl;
    
    cutFile << "Associated hadron Cuts" << endl << endl;
    cutFile << "associated hadron Pt   -\t" << "Low: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(17) << "\t" << "High: " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(18) << endl;
    cutFile << "Num of events used to mix  -\t" << ": " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(19) << endl;
    cutFile << "HFT Tracks only??  -\t" << "YES(1), NO(0): " << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(22) << endl;
    cutFile << "Track Chi2 Max:  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(23) << endl;
    cutFile << "Track DCA Max:  -\t" << ((TH1D*) file->Get("HistOfCuts"))->GetBinContent(24) << endl;
    
    cutFile.close();

//---------------------------------------------    
//BEGIN INVARIANT MASS HISTOGRAM INFORMATION   
//---------------------------------------------
    
    ofstream invariantMassInfo;
    TString invariantMassInfoName = "Signal_and_background_information.txt";
    TString invariantMassInfoFile = path + invariantMassInfoName;
    invariantMassInfo.open (invariantMassInfoFile);
    
    
    for(int k = 0; k < NUM_PT_BINS; k++){  //Initialize the histograms 
        for(int i = 0; i < 3; i++){
        
            //if(k == 3) { LabelVersion = 1; }
       
            tmp = oneDUSStrings[k] + binLabelCentralityClass[i];
            oneDUSHistos[k][i] = (TH1D*)file->Get(tmp);     //US spectrum
            tmp = oneDLSStrings[k] + binLabelCentralityClass[i];
            oneDLSHistos[k][i] = (TH1D*)file->Get(tmp);     //LS spectrum
        
            oneDScaledLSHistos[k][i] = (TH1D*)oneDLSHistos[k][i]->Clone();  //clone of LS spectrum to be scaled   
            tmp = oneDScaledLSStrings[k] + binLabelCentralityClass[i];
            oneDScaledLSHistos[k][i]->SetTitle(oneDScaledLSStrings[k]);        
            tmp = oneDSubtractedInvMassStrings[k] + binLabelCentralityClass[i];
            oneDSubtractedInvMassHistos[k][i] = new TH1D(tmp, tmp, 50, 1.6, 2.1); //final spectrum -- US minus LS (scaled)
       
            oneDUSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
            oneDUSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDUSHistos[k][i]->GetYaxis()->SetTitleOffset(1.2);
            oneDLSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
            oneDLSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDLSHistos[k][i]->GetYaxis()->SetTitleOffset(1.2);
            oneDScaledLSHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
            oneDScaledLSHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDScaledLSHistos[k][i]->GetYaxis()->SetTitleOffset(1.2);
            oneDSubtractedInvMassHistos[k][i]->GetXaxis()->SetTitle("Invariant Mass GeV/c^{2}");
            oneDSubtractedInvMassHistos[k][i]->GetYaxis()->SetTitle("counts");
            oneDSubtractedInvMassHistos[k][i]->GetYaxis()->SetTitleOffset(1.2);
        
            xAxisUS = oneDUSHistos[k][i]->GetXaxis();
            xAxisLS = oneDLSHistos[k][i]->GetXaxis();
        
            binMassLowUS  = oneDUSHistos[k][i]->GetXaxis()->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
            binMassHighUS = oneDUSHistos[k][i]->GetXaxis()->FindBin(sideBandHigh);
            binMassLowLS = oneDLSHistos[k][i]->GetXaxis()->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
            binMassHighLS = oneDLSHistos[k][i]->GetXaxis()->FindBin(sideBandHigh);
        
            countsSideBandUS = oneDUSHistos[k][i]->Integral(binMassLowUS, binMassHighUS);
            countsSideBandLS = oneDLSHistos[k][i]->Integral(binMassLowLS, binMassHighLS);
        
            //Calculate peak region values
        
            binMassLowUS  = xAxisUS->FindBin(massLow);            
            binMassHighUS = xAxisUS->FindBin(massHigh);
            binMassLowLS  = xAxisLS->FindBin(massLow);            
            binMassHighLS = xAxisLS->FindBin(massHigh);
        
            oneDUSHistos[k][i]->Sumw2();
            oneDUSHistos[k][i]->SetMarkerStyle(20);
            oneDLSHistos[k][i]->Sumw2();
            oneDLSHistos[k][i]->SetMarkerStyle(20);
            oneDSubtractedInvMassHistos[k][i]->Sumw2();
            oneDSubtractedInvMassHistos[k][i]->SetMarkerStyle(20);
        
            oneDUSHistos[k][i]->Draw();
                
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDUSStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
        
            //oneDUSHistos[i]->SetTitle("");
            //oneDUSHistos[i]->Draw();
            //tmp = path + oneDUSStrings[i] + fileTypeEps;
            //c->SaveAs(tmp);
        
            LSScaleFactor[k][i] = countsSideBandUS/countsSideBandLS;
        
            oneDScaledLSHistos[k][i]->Scale(LSScaleFactor[k][i]);            //normalize LS spectrum here
        
            countsPeakRegionUS = oneDUSHistos[k][i]->Integral(binMassLowUS, binMassHighUS);
            countsPeakRegionLS = oneDScaledLSHistos[k][i]->Integral(binMassLowLS, binMassHighLS);  //this gets an estimate for B from the LS scaled to the BG in the US, using a single sideband
        
            invariantMassInfo << "____________________________________________________________________________________________________" << endl;
            invariantMassInfo << "PtBin (0: 0-1, 1: 1-4, 2: 4-20, 3: all): " << k << endl;
            invariantMassInfo << "CentBin (0: per, 1: mid, 2: cent, 3: all): " << i << endl;
            //cout << "LS Scale Factor: " << LSScaleFactor[i] << endl;
            invariantMassInfo << "S+B (from integral in mass window): " << countsPeakRegionUS << endl;
        
            oneDSubtractedInvMassHistos[k][i]->Add(oneDUSHistos[k][i], oneDScaledLSHistos[k][i], 1, -1); // form subtracted spectrum here
        
            oneDSubtractedInvMassHistosFunction[k][i] = (TH1D*)oneDSubtractedInvMassHistos[k][i]->Clone();
            tmp = oneDSubtractedInvMassStringsFunction[k] + binLabelCentralityClass[i];
            oneDSubtractedInvMassHistosFunction[k][i]->SetTitle(tmp);
        
            oneDSubtractedInvMassHistos[k][i]->Fit(g1, "qR0");
            oneDSubtractedInvMassHistos[k][i]->Fit(e1, "qR0+");
       
            g1->GetParameters(&par[0]);
            e1->GetParameters(&par[3]);
        
            fun->SetParameters(par);
            oneDSubtractedInvMassHistos[k][i]->Fit(fun, "qR+");
            oneDSubtractedInvMassHistos[k][i]->SetMarkerColor(2);
        
            BGFunc->SetParameter(0, fun->GetParameter(3));
            BGFunc->SetParameter(1, fun->GetParameter(4));
        
            //------------------------This is where things are controversial -- how do I fit the residual background????-------------------
        
        
            TF1 *fl = new TF1("fl",fline,1.6,2.1,2);
            fl->SetParameters(2,-1);
            //fit only the linear background excluding the signal area
            reject = kTRUE;
            oneDSubtractedInvMassHistos[k][i]->Fit(fl,"qR");
            reject = kFALSE;
        
            oneDSubtractedInvMassHistosFunction[k][i]->Add(fl, -1);
        
        
            //-----------------------------------------------------------------------------------------------------------------------------
        
            oneDSubtractedInvMassHistosFunction[k][i]->Fit(finalGaussian, "qR");
        
            integralFit = oneDSubtractedInvMassHistosFunction[k][i]->Integral(binMassLowUS, binMassHighUS);
        
            //integralFit is the estimate for "S"
            //countsPeakRegionUS is S+B
            //B = countsPeakRegionUS - integralFit
        
            BOverSPlusBFit[k][i] = (countsPeakRegionUS - integralFit)/countsPeakRegionUS;    //  B/S+B scale factor from background subtraction from fit
            SPlusBOverSFit[k][i] = countsPeakRegionUS/integralFit;                           //  S+B/S used as a final scale factor on the final correlation quantity -- from fit
        
            invariantMassInfo << "S (from fit estimate): " << integralFit << endl;
            invariantMassInfo << "B (from fit estimate): " << (countsPeakRegionUS-integralFit) << endl;
            invariantMassInfo << "S+B/S from LS subtraction + expo Fit:  " << SPlusBOverSFit[k][i] << endl;
            invariantMassInfo << "B/S+B from LS subtraction + expo Fit: " << BOverSPlusBFit[k][i] << endl << endl;
        
            oneDUSHistos[k][i]->Draw();
            oneDScaledLSHistos[k][i]->Draw("SAME");
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDUSStrings[k] + binLabelCentralityClass[i] + withFitLabel + fileType;
            c->SaveAs(tmp);
            oneDLSHistos[k][i]->Draw();
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDLSStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
            oneDSubtractedInvMassHistos[k][i]->Draw();
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDSubtractedInvMassStrings[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
            oneDSubtractedInvMassHistosFunction[k][i]->Draw();
            tmp = path + invariantMassHistogramsFolder + invariantMassOutputFolders[k] + oneDSubtractedInvMassStringsFunction[k] + binLabelCentralityClass[i] + fileType;
            c->SaveAs(tmp);
        }
    }
    
    invariantMassInfo.close();

    cout << endl << "Invariant Mass histograms process and S & B calculated for each Pt/Centrality bin" << endl << endl;
    
//------------------------------------------
//BEGIN CORRELATION HISTOGRAM INFORMATION
//------------------------------------------

    /**********************************************************************************
        EXTRACT SIBLING AND MIXED HISTOGRAMS FROM RAW DATA FILE HERE
    **********************************************************************************/
    for(int band = 0; band < 4; band++){ //begin loop to get raw sibling and mixed histograms in both US and LS
        for(int i = 0; i < NUM_CENT_BINS; i++){
            for(int k = 0; k < NUM_PT_BINS; k++){
    
            //Sibling histogram information
        
            str1 = SibCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
            sibCorrBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
        
            str1 = SibCorrLabels[band] + phiProj + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
            sibCorrBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            formatCorrHist(sibCorrBinVzInt[band][k][i]);
        
            //mixed histogram information
            
            str1 = MixCorrLabels[band]  + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
            mixCorrBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
        
            str1 = MixCorrLabels[band]  + phiProj + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
            mixCorrBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            formatCorrHist(mixCorrBinPhiProjVzInt[band][k][i]);
        
        
            for(int j = 0; j < 10; j++){
            
                str1 = SibCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                str2 = MixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                
                sibCorrBin[band][k][j][i] = (TH2D*)file->Get(str1);
                mixCorrBin[band][k][j][i] = (TH2D*)file->Get(str2);
                
                sibCorrBin[band][k][j][i]->Sumw2();
                mixCorrBin[band][k][j][i]->Sumw2();
                
                //This gets the stat errors for the raw sibling and mixed histograms
                for(int etaBin = 0; etaBin < NUM_ETA_BINS; etaBin++){
                    for(int phiBin = 0; phiBin < NUM_PHI_BINS; phiBin++){
                        
                        sibCorrBinStatErrors[band][k][j][i][etaBin][phiBin] = sibCorrBin[band][k][j][i]->GetBinContent(etaBin, phiBin);    
                        mixCorrBinStatErrors[band][k][j][i][etaBin][phiBin] = mixCorrBin[band][k][j][i]->GetBinContent(etaBin, phiBin);
                    }
                }//Stat errors
                
                
                if(k==3){totalPairsSibling = totalPairsSibling + sibCorrBin[band][k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);}
                
                sibCorrBinVzInt[band][k][i]->Add(sibCorrBin[band][k][j][i]);
                formatCorrHist(sibCorrBin[band][k][j][i]);  //histogram formatting for axes
                mixCorrBinVzInt[band][k][i]->Add(mixCorrBin[band][k][j][i]);
                formatCorrHist(mixCorrBin[band][k][j][i]);  //histogram formatting for axes
                
            
                //Sibling US
                if(PRINT_SUB_HISTOS_RAW){
                    sibCorrBin[band][k][j][i]->Draw("SURF1");                              //sibling 2D US histogram
                    str1 = path + outputFolders[0] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + SibCorrLabels[band] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;
                    c->SaveAs(str1);
                }    
        
                str1 = phiProj + SibCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];         //sibling phi projection US
                sibCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());  
                sibCorrBinPhiProj[band][k][j][i] = (TH1D*)sibCorrBin[band][k][j][i]->ProjectionY();  
                if(PRINT_SUB_HISTOS_RAW){
                    str1 = path + outputFolders[0] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + SibCorrLabels[band] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;        
                    sibCorrBinPhiProj[band][k][j][i]->Draw();
                    c->SaveAs(str1);
                }
                
                sibCorrBinPhiProjVzInt[band][k][i]->Add(sibCorrBinPhiProj[band][k][j][i]);
                         
                //Mixed 
                if(PRINT_SUB_HISTOS_RAW){
                    mixCorrBin[band][k][j][i]->Draw("SURF1");                              //mixling 2D US histogram
                    str1 = path + outputFolders[1] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + MixCorrLabels[band] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;
                    c->SaveAs(str1);
                }
                
                str1 = phiProj + MixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];     //mixling phi projection US
                mixCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());  
                mixCorrBinPhiProj[band][k][j][i] = (TH1D*)mixCorrBin[band][k][j][i]->ProjectionY();  
                
                if(PRINT_SUB_HISTOS_RAW){
                    str1 = path + outputFolders[1] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + MixCorrLabels[band] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;        
                    mixCorrBinPhiProj[band][k][j][i]->Draw();
                    c->SaveAs(str1);
                }
        
                mixCorrBinPhiProjVzInt[band][k][i]->Add(mixCorrBinPhiProj[band][k][j][i]);
                         
            }
        
            if(PRINT_SUB_HISTOS_RAW){
                
                //sibling
                sibCorrBinVzInt[band][k][i]->Draw("SURF1");
                str1 = path + outputFolders[0] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + SibCorrLabels[band] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                c->SaveAs(str1);
        
                sibCorrBinPhiProjVzInt[band][k][i]->Draw();
                str1 = path + outputFolders[0] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + SibCorrLabels[band] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                c->SaveAs(str1);
        
                //mixed
                mixCorrBinVzInt[band][k][i]->Draw("SURF1");
                str1 = path + outputFolders[1] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + MixCorrLabels[band] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                c->SaveAs(str1);
        
                mixCorrBinPhiProjVzInt[band][k][i]->Draw();
                str1 = path + outputFolders[1] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + MixCorrLabels[band] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                c->SaveAs(str1);
        
           }
       }
    }//end loop to get raw sibling and mixed histograms in both US and LS
    
    cout << "Total sibling pairs for band " << band << " : " << totalPairsSibling << endl;
    
   }

   cout << endl << "Sibling and mixed histograms extracted" << endl;
   
    /**********************************************************************************
        NORMALIZE THE MIXED HISTOGRAMS HERE
    **********************************************************************************/
   
  for(int band = 0; band < 4; band++){ // begin loop to make scaled mixed histograms
    for(int i = 0; i < NUM_CENT_BINS; i++){ 
        for(int k = 0; k < NUM_PT_BINS; k++){
            
            str1 = ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
            scaledMixCorrBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            
            str1 = phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
            scaledMixCorrBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                  
            for(int j = 0; j < 10; j++){
        
                scaledMixCorrBin[band][k][j][i] = (TH2D*) mixCorrBin[band][k][j][i]->Clone();
                
                //US histograms
                numSib = sibCorrBin[band][k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                numMixUS = mixCorrBin[band][k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                integralRawMixHistos[band][k][j][i] = numMixUS;
                scaledMixCorrBin[band][k][j][i]->Scale(numSib/numMixUS);
                //VzScaleFactorUS[j] = numMixUS;
                
                for(int etaBin = 0; etaBin < NUM_ETA_BINS; etaBin++){
                    for(int phiBin = 0; phiBin < NUM_PHI_BINS; phiBin++){
                        
                        scaledMixCorrBinStatErrors[band][k][j][i][etaBin][phiBin] = mixCorrBinStatErrors[band][k][j][i][etaBin][phiBin]/((numSib/numMixUS)*(numSib/numMixUS));    
                        
                    }
                }//Stat errors
         
         
                //US scaled histograms
                
                str1 = ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
        
                formatCorrHist(scaledMixCorrBin[band][k][j][i], str1);  //formatting
        
                if(PRINT_SUB_HISTOS_SCALED){
                    scaledMixCorrBin[band][k][j][i]->Draw("SURF1");
                    str1 = path + outputFolders[2] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;
                    c->SaveAs(str1);
                }
                
                str1 = phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];       //scaled mix phi projection US 
                scaledMixCorrBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());  
                scaledMixCorrBinPhiProj[band][k][j][i] = (TH1D*)scaledMixCorrBin[band][k][j][i]->ProjectionY();  
                
                if(PRINT_SUB_HISTOS_SCALED){
                    str1 = path + outputFolders[2] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;        
                    scaledMixCorrBinPhiProj[band][k][j][i]->Draw();
                    c->SaveAs(str1);
                }
            }
                
            if(PRINT_SUB_HISTOS_SCALED){
            
                str1 = path + outputFolders[2] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;        
                scaledMixCorrBinVzInt[band][k][i]->Draw("SURF1");
                c->SaveAs(str1);
        
                str1 = path + outputFolders[2] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + ScaledMixCorrLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;        
                scaledMixCorrBinPhiProjVzInt[band][k][i]->Draw();
                c->SaveAs(str1);
        
           }
        }
      }  //end loop to make scaled mixed histograms  */
   }
   //--------------------------------------------------------calculate some useful quantities here--------------------------------------
   
   cout << endl << "Mixed Histograms Normalized." << endl << endl;
   
   for(int band = 0; band < 4; band++){                               //build arrays to store information for adding the various alpha's (Vz, centrality, etc.)
       for(int k = 0; k < NUM_PT_BINS; k++){
           
           integralRawMixHistosVzIntCentInt[band][k][0] = 0;
           integralRawMixHistosVzIntCentInt[band][k][1] = 0;
           integralRawMixHistosVzIntCentInt[band][k][2] = 0;
   
           for(int i = 0; i < NUM_CENT_BINS; i++){
            
                integralRawMixHistosVzInt[band][k][i] = 0;
                //cout << band << "   " << k << "   "<< i << endl;
               
                for(int j = 0; j < 10; j++){
            
                    integralRawMixHistosVzInt[band][k][i] = integralRawMixHistosVzInt[band][k][i] + integralRawMixHistos[band][k][j][i];
                    
                }   
                
                        
           }
       }
    }   
 //------------------------------------------------------------------------------------------------------------------------------------
 
   cout << "Weight Factors Calculated." << endl << endl;
    /**********************************************************************************
        CALCULATE DELTA RHO HERE
    **********************************************************************************/
 
    for(int band = 0; band < 4; band++){ // begin loop to make sib minus scaledMixed histos (delRho)
        for(int i = 0; i < NUM_CENT_BINS; i++){ 
            for(int k = 0; k < NUM_PT_BINS; k++){
            
                str1 = SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                sibMinusScaledMixVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                str1 = phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                sibMinusScaledMixPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           
                for(int j = 0; j < 10; j++){
        
                
                    str1 = SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                
                    sibMinusScaledMix[band][k][j][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                    sibMinusScaledMix[band][k][j][i]->Add(sibCorrBin[band][k][j][i], scaledMixCorrBin[band][k][j][i], 1, -1);
                    //integralError = 0.0;
                    for(int etaBin = 0; etaBin < NUM_ETA_BINS; etaBin++){
                        for(int phiBin = 0; phiBin < NUM_PHI_BINS; phiBin++){
                        
                            sibMinusScaledMixStatErrors[band][k][j][i][etaBin][phiBin] = sibCorrBinStatErrors[band][k][j][i][etaBin][phiBin] + scaledMixCorrBinStatErrors[band][k][j][i][etaBin][phiBin];    
                            integralError[phiBin] = integralError[phiBin] + sibMinusScaledMixStatErrors[band][k][j][i][etaBin][phiBin];
                        }
                    }//Stat errors
                
                      
                
                    formatCorrHist(sibMinusScaledMix[band][k][j][i]); //formatting
        
                    if(PRINT_SUB_HISTOS_DEL_RHO){
                        sibMinusScaledMix[band][k][j][i]->Draw("SURF1");
                        str1 = path + outputFolders[3] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + SibMinusScaledMixLabels[band] + binLabelCent[i] + fileType;
                        c->SaveAs(str1);
                   
                    }
                   
                    str1 = phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    sibMinusScaledMixPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());  
                    sibMinusScaledMixPhiProj[band][k][j][i] = (TH1D*)sibMinusScaledMix[band][k][j][i]->ProjectionY();  
                
                    if(PRINT_SUB_HISTOS_DEL_RHO){
                        str1 = path + outputFolders[3] +  bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;        
                        
                        for(int phiBin = 0; phiBin < NUM_PHI_BINS; phiBin++){
                                                
                            sibMinusScaledMixPhiProj[band][k][j][i]->SetBinError(phiBin+1, integralError[phiBin]);
                        }
                        
                        sibMinusScaledMixPhiProj[band][k][j][i]->Draw("E0");
                        c->SaveAs(str1);
                    }
                
                    sibMinusScaledMixVzInt[band][k][i]->Add(sibMinusScaledMix[band][k][j][i]);
                    sibMinusScaledMixPhiProjVzInt[band][k][i]->Add(sibMinusScaledMixPhiProj[band][k][j][i]);
                }
     
                if(PRINT_SUB_HISTOS_DEL_RHO){ 
                    str1 = path+outputFolders[3]+bandFolders[band]+PtSubFolders[k]+VzSubFolders[10]+SibMinusScaledMixLabels[band]+PtBinLabel[k]+binLabelPt[k]+VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                    sibMinusScaledMixVzInt[band][k][i]->Draw("SURF1");
                    c->SaveAs(str1);
             
                    str1 = path + outputFolders[3] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + SibMinusScaledMixLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                    sibMinusScaledMixPhiProjVzInt[band][k][i]->Draw();
                    c->SaveAs(str1);
                }                
            }
        }// end loop to make sibling minus scaled mixed histograms (delRho)
    }
  
   cout << "deltaRho for each band (US, LS, SBL, SBR) calculated." << endl << endl;
  
   //-----------------------------------------------------------------------------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------------------------------------------------------------------------
   //-----------------------------------------BEGIN SECTION TO CALCULATE ACTUAL CORRELATIONS--------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------------------------------------------------------------------------
   
   
    for(int band = 0; band < 4; band++){     // begin loop to make delRho/RhoRef US histos -- this will calculate the delRho/Rho histograms for all 4 bands -- LS, US, SBR, SBL
        for(int i = 0; i < NUM_CENT_BINS; i++){ 
            for(int k = 0; k < NUM_PT_BINS; k++){
            
                str1 = delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                delRhoOverRhoRefBinVzInt[band][k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                str1 = phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i];
                delRhoOverRhoRefBinPhiProjVzInt[band][k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            
                for(int j = 0; j < 10; j++){

                    str1 = delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];
                    delRhoOverRhoRefBin[band][k][j][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                    
                    if(scaledMixCorrBin[band][k][j][i]->GetEntries() > 0){
                    
                        delRhoOverRhoRefBin[band][k][j][i]->Divide(sibMinusScaledMix[band][k][j][i], scaledMixCorrBin[band][k][j][i], 1, 1);
                    }
                    
                    formatCorrHist(delRhoOverRhoRefBin[band][k][j][i]);
                    //VzScaleFactorUS[j] = delRhoOverRhoRefUSBin[k][j][i]->Integral(1, NUM_ETA_BINS, 1, NUM_PHI_BINS);
                    delRhoOverRhoRefBin[band][k][j][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                
                    if(PRINT_DELRHO_OVER_REF_HISTOS){ 
                        delRhoOverRhoRefBin[band][k][j][i]->Draw("SURF1");
                        str1 = path + outputFolders[4] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;
                        c->SaveAs(str1);
                    }    
        
                    str1 = phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i];        //sib - scaled mix phi projection US 
                    delRhoOverRhoRefBinPhiProj[band][k][j][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());  
                    delRhoOverRhoRefBinPhiProj[band][k][j][i] = (TH1D*)delRhoOverRhoRefBin[band][k][j][i]->ProjectionY();  
                
                    if(PRINT_DELRHO_OVER_REF_HISTOS){               
                        str1 = path + outputFolders[4] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[j] + phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzBinLabel + binLabelVz[j] + centBinLabel + binLabelCent[i] + fileType;      
                        delRhoOverRhoRefBinPhiProj[band][k][j][i]->Draw();
                        c->SaveAs(str1);
                    }
                
                    scaleFactorVz = integralRawMixHistos[band][k][j][i]/integralRawMixHistosVzInt[band][k][i];
                    delRhoOverRhoRefBinVzInt[band][k][i]->Add(delRhoOverRhoRefBin[band][k][j][i], scaleFactorVz);
                
                }
 
                str1 = path + outputFolders[4] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType; 
                delRhoOverRhoRefBinVzInt[band][k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                if(PRINT_DELRHO_OVER_REF_HISTOS){ 
                    delRhoOverRhoRefBinVzInt[band][k][i]->Draw("SURF1");
                    c->SaveAs(str1);
                    
                }
            
                delRhoOverRhoRefBinPhiProjVzInt[band][k][i] = (TH1D*)delRhoOverRhoRefBinVzInt[band][k][i]->ProjectionY();
            
                str1 = path + outputFolders[4] + bandFolders[band] + PtSubFolders[k] + VzSubFolders[10] + phiProj + delRhoOverRhoRefLabels[band] + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType;
                if(PRINT_DELRHO_OVER_REF_HISTOS){ 
                    delRhoOverRhoRefBinPhiProjVzInt[band][k][i]->Draw();
                    c->SaveAs(str1);
                }    
            }
        } //end loop to make delRho/Rho_ref histograms*/
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------
    //---------------------------THIS SECTION DOES THE FINAL SUBTRACTION IN THE 16 MULTIPLICITY BINS------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------
    
        for(int i = 0; i < NUM_CENT_BINS; i++){// begin side band subtraction here -- this will take the US delRho/Rho and subtract the SB average delRho/Rho from it
            for(int k = 0; k < NUM_PT_BINS; k++){
           
                str1 = sideBandAverageLabel + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                sideBandAverage[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                str1 = phiProj + sideBandAverageLabel + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                sideBandAveragePhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                
                sideBandAverage[k][i]->Add(delRhoOverRhoRefBinVzInt[0][k][i], delRhoOverRhoRefBinVzInt[2][k][i], .5, .5);
                sideBandAveragePhiProj[k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[0][k][i], delRhoOverRhoRefBinPhiProjVzInt[2][k][i], .5, .5);
                
                str1 = subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrSideBandCent[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                str1 = phiProj + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrSideBandCentPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           
                fullySubtractedCorrSideBandCent[k][i]->Add(delRhoOverRhoRefBinVzInt[1][k][i], sideBandAverage[k][i], 1, -BOverSPlusBFit[k][3]);
               
                str1 = path + outputFolders[4] + bandFolders[4] + PtSubFolders[k] + VzSubFolders[10] + sideBandAverageLabel + PtBinLabel[k] + binLabelPt[k] + VzIntLabel + centBinLabel + binLabelCent[i] + fileType; 
                sideBandAverage[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                sideBandAverage[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                str1 = path + outputFolders[5] + outputFolders[7] + PtSubFolders[k] + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrSideBandCent[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrSideBandCent[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                fullySubtractedCorrSideBandCent[k][i]->Write();
                
                str1 = path + outputFolders[5] + outputFolders[7] + PtSubFolders[k] + phiProj + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrSideBandCentPhiProj[k][i] = (TH1D*)fullySubtractedCorrSideBandCent[k][i]->ProjectionY();
                fullySubtractedCorrSideBandCentPhiProj[k][i]->Draw();
                c->SaveAs(str1);
           
           
            }    
        }// end side band subtraction here
        
        
        for(int i = 0; i < NUM_CENT_BINS; i++){// begin LS subtraction here  -- this will take the US delRho/Rho and subtract the LS delRho/Rho from it
            for(int k = 0; k < NUM_PT_BINS; k++){
           
                str1 = subtractedCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrLSCent[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                str1 = phiProj + subtractedCorrLabels[1] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i];
                fullySubtractedCorrLSCentPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           
                fullySubtractedCorrLSCent[k][i]->Add(delRhoOverRhoRefBinVzInt[1][k][i], delRhoOverRhoRefBinVzInt[3][k][i], 1, -BOverSPlusBFit[k][3]);
                
                //fullySubtractedCorrLSCentPhiProj[k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[1][k][i], delRhoOverRhoRefBinPhiProjVzInt[3][k][i], 1, -BOverSPlusB[k]);
            
                str1 = path + outputFolders[5] + outputFolders[6] + PtSubFolders[k] + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrLSCent[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrLSCent[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                str1 = path + outputFolders[5] + outputFolders[6] + PtSubFolders[k] + phiProj + subtractedCorrLabels[0] + PtBinLabel[k] + binLabelPt[k] + centBinLabel + binLabelCent[i] + fileType; 
                fullySubtractedCorrLSCentPhiProj[k][i] = (TH1D*) fullySubtractedCorrLSCent[k][i]->ProjectionY();
                fullySubtractedCorrLSCentPhiProj[k][i]->Draw();
                c->SaveAs(str1);
           
           
            }    
        }// end LS subtraction here
          

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------      
    //---------------------------THIS SECTION DOES THE FINAL SUBTRACTION IN THE 3 centrality BINS----------------------------------------------------------------      
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------

    
          
      
    
    TString bandSubtract[3]          = {"LS_","SideBand_", "US"};
    int temp = 0;
    int centralityCombinations[6] = {0, 2, 3, 8, 9, 15}; //This was calculated for combining multiplicity bins into the proper centrality bins
    int numBinsInCentrality = 0;
   
    for(int k = 0; k < NUM_PT_BINS; k++){
        for(int i = 0; i < 3; i++){
           
           
           str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[2] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralityUS[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           
           str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralityLS[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           str1 =  fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralitySideBand[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           str1 =  phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralityLSPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           str1 =  phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i];
           fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           
           for(int band = 0; band < 4; band++){
           
                delRhoOverRhoRefBinVzIntCentInt[band][k][i] = new TH2D("","",  NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                delRhoOverRhoRefBinVzIntCentIntPhiProj[band][k][i] = new TH1D("","", NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
           
           
           }
           
           numBinsInCentrality = centralityCombinations[(2*i)+1] - centralityCombinations[(2*i)];
           
           for(int band = 0; band < 4; band++){
                for(int bin = 0; bin < numBinsInCentrality+1; bin++){
               
                    temp = centralityCombinations[(2*i)]+bin;
                    integralRawMixHistosVzIntCentInt[band][k][i] = integralRawMixHistosVzIntCentInt[band][k][i] + integralRawMixHistosVzInt[band][k][temp];
           
                }
           }
       }    
       
       cout << endl << endl;
        
       for(int i = 0; i < 3; i++){        
           
           numBinsInCentrality = centralityCombinations[(2*i)+1] - centralityCombinations[(2*i)];
           
           for(int bin = 0; bin < numBinsInCentrality+1; bin++){
           
                temp = centralityCombinations[(2*i)]+bin;
                
                scaleFactorVzCent = integralRawMixHistosVzInt[1][k][temp]/integralRawMixHistosVzIntCentInt[1][k][i];                 //US
                delRhoOverRhoRefBinVzIntCentInt[1][k][i]->Add(delRhoOverRhoRefBinVzInt[1][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[1][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[1][k][temp], scaleFactorVzCent);
                
                cout << "US cent Scale Factor (cent bin: " << i << ") :" << scaleFactorVzCent << endl;
                
                scaleFactorVzCent = integralRawMixHistosVzInt[3][k][temp]/integralRawMixHistosVzIntCentInt[3][k][i];                 //LS
                delRhoOverRhoRefBinVzIntCentInt[3][k][i]->Add(delRhoOverRhoRefBinVzInt[3][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[3][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[3][k][temp], scaleFactorVzCent);
                
                scaleFactorVzCent = integralRawMixHistosVzInt[2][k][temp]/integralRawMixHistosVzIntCentInt[2][k][i];                 //Right SB
                delRhoOverRhoRefBinVzIntCentInt[2][k][i]->Add(delRhoOverRhoRefBinVzInt[2][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[2][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[2][k][temp], scaleFactorVzCent);
                
                scaleFactorVzCent = integralRawMixHistosVzInt[0][k][temp]/integralRawMixHistosVzIntCentInt[0][k][i];                 //Left SB
                delRhoOverRhoRefBinVzIntCentInt[0][k][i]->Add(delRhoOverRhoRefBinVzInt[0][k][temp], scaleFactorVzCent);
                delRhoOverRhoRefBinVzIntCentIntPhiProj[0][k][i]->Add(delRhoOverRhoRefBinPhiProjVzInt[0][k][temp], scaleFactorVzCent);
            
           }
       }    
    }   
      
       double phiProjectionScaleFactor = 1.0/9.0;
       
        for(int k = 0; k < 4; k++){
            for(int i = 0; i < 3; i++){
      
                cout << "k, i: " << k << " , " << i << endl;
               
                fullySubtractedCorrCustomCentralityLS[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[1][k][i], delRhoOverRhoRefBinVzIntCentInt[3][k][i], 1, -BOverSPlusBFit[k][i]);
                fullySubtractedCorrCustomCentralityLS[k][i]->Scale(SPlusBOverSFit[k][i]);
                
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[1][k][i]);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[0][k][i], -.5*BOverSPlusBFit[k][i]);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Add(delRhoOverRhoRefBinVzIntCentInt[2][k][i], -.5*BOverSPlusBFit[k][i]);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Scale(SPlusBOverSFit[k][i]);
                
                fullySubtractedCorrCustomCentralityLSPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralityLS[k][i]->ProjectionY();
                fullySubtractedCorrCustomCentralityLSPhiProj[k][i]->Scale(phiProjectionScaleFactor);
                
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i] = (TH1D*) fullySubtractedCorrCustomCentralitySideBand[k][i]->ProjectionY();
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Scale(phiProjectionScaleFactor);
       
                str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
                fullySubtractedCorrCustomCentralityLS[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrCustomCentralityLS[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                fullySubtractedCorrCustomCentralityLS[k][i]->Write();
           
                str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
                fullySubtractedCorrCustomCentralitySideBand[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
                c->SaveAs(str1);
                
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Write();
           
                str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[0] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
                fullySubtractedCorrCustomCentralityLSPhiProj[k][i]->Draw();
                c->SaveAs(str1);
           
                str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fileType;
                fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw();
                c->SaveAs(str1);
            }    
        }
 
 
 
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------FITTING AND PARAMETER EXTRACTION DONE HERE---------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    TCanvas *c = new TCanvas("c3", "Histograms", 1100, 850);
 
    TF2 *myfit   = new TF2("parameterFitFunction","[0] + [1]*cos(y) + [2]*2*cos(2*y) + [3]*(exp(-0.5*(y*y)/([4]*[4]))+exp(-0.5*((y-6.28)*(y-6.28))/([4]*[4])))*exp(-0.5*(x*x)/([5]*[5]))", -ETA_RANGE, ETA_RANGE,-1.57,4.71);
    
    TF1 *myfit1D = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x)" , -1.57, 4.71);
    
    TF1 *fitProjection = new TF1("parameterFitFunction","[0] + [1]*cos(x) + [2]*2*cos(2*x) + [3]*[4]*(exp(-0.5*(x*x)/([5]*[5]))+exp(-0.5*((x-6.28)*(x-6.28))/([5]*[5])))" , -1.57, 4.71);
    
    TF1 *quadrupole1    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    TF1 *quadrupole2    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    TF1 *quadrupole3    = new TF1("Quadrupole", "[0]*2*cos(2*x)", -1.57, 4.71);
    
    TF1 *pileupLine1    = new TF1("PileupLine1", "[0] + [1]*x", -2.0, -1.0);
    TF1 *pileupLine2    = new TF1("PileupLine2", "[0] + [1]*x", -1.0, 0);
    TF1 *pileupLine3    = new TF1("PileupLine3", "[0] + [1]*x", 0, 1.0);
    TF1 *pileupLine4    = new TF1("PileupLine4", "[0] + [1]*x", 1.0, 2.0);
    
    TH1D *testPileupFunc = new TH1D("test", "test", 100, -2, 2);
    
    //testPileupFunc->SetMaximum(1);
    //testPileupFunc->SetMinimum(-1);
    
    pileupLine1->FixParameter(0, -2);
    pileupLine1->FixParameter(1, -1);
    pileupLine2->FixParameter(0, 1);
    pileupLine2->FixParameter(1, 2);
    pileupLine3->FixParameter(0, 1);
    pileupLine3->FixParameter(1, -2);
    pileupLine4->FixParameter(0, -2);
    pileupLine4->FixParameter(1, 1);
    
    str1 = path + outputFolders[5] + outputFolders[8] + fileType;
    testPileupFunc->Draw();
    pileupLine1->Draw("SAME");
    pileupLine2->Draw("SAME");
    pileupLine3->Draw("SAME");
    pileupLine4->Draw("SAME");
    
    c->SaveAs(str1); 
 
    TH2D* fullySubtractedCorrCustomCentralitySideBandFit[4][3];
    TH2D* fullySubtractedCorrCustomCentralitySideBandRes[4][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[4][3];
    TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[4][3];
    
    TH2D* fullySubtractedCorrCustomCentralitySideBandPileupCorrected[4][3];
    //TH2D* fullySubtractedCorrCustomCentralitySideBandRes[4][3];
    //TH1D* fullySubtractedCorrCustomCentralitySideBandFitPhiProj[4][3];
    //TH1D* fullySubtractedCorrCustomCentralitySideBandResPhiProj[4][3];
    
    TString fitLabel = "_Fit";
    TString resLabel = "_Residual";
    TString pileupLabel = "PileupCorrected";
    
    ofstream fitParameterFile2D;
    TString fitParameterFileName2D = "Extracted_Parameters_from_fitting_2D.txt";
    TString fitParamterOutputFile2D = path + fitParameterFileName2D;
    fitParameterFile2D.open (fitParamterOutputFile2D);
    
    ofstream fitParameterFile1D;
    TString fitParameterFileName1D = "Extracted_Parameters_from_fitting_1D.txt";
    TString fitParamterOutputFile1D = path + fitParameterFileName1D;
    fitParameterFile1D.open (fitParamterOutputFile1D);
    
    double v2_2D;
    double v2_2D_Error;
    double v2_1D;
    double v2_1D_Error;
    
    int bin1;
    int bin2;
    int bin3;
    int bin4;
    
    double GaussIntegralParameter;
 
    for(int k = 0; k < 4; k++){
        for(int i = 0; i < 3; i++){ 
 
           
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFit[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            
            str1 = phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel;
            fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
            str1 = phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = new TH1D(str1, str1, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
 
            str1 = fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + pileupLabel;
            fullySubtractedCorrCustomCentralitySideBandPileupCorrected[k][i] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
 
            fitParameterFile2D << "_________________________________________________________________________________________________________" << endl;
            fitParameterFile1D << "_________________________________________________________________________________________________________" << endl;
 
           
            fullySubtractedCorrCustomCentralitySideBandPileupCorrected[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone();

              
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->GetXaxis()->SetRangeUser(-ETA_RANGE, ETA_RANGE);

            myfit1D->SetParLimits(0, -.1, .1); //offset A0
            myfit1D->SetParLimits(1, -.1, .1); //v1     AD
            myfit1D->SetParLimits(2, -.1, .1); //v2     AQ
            
            myfit->SetParLimits(0, -.1, .1);    //offset A0
            myfit->SetParLimits(1, -.1,  .1);    //v1     AD
            myfit->SetParLimits(2, -.1,  .1);    //v2     AQ
            myfit->SetParLimits(3,  0,  1);    //Jetamp A1
            myfit->SetParLimits(4, 0.1, 3);     //sigPhi
            myfit->SetParLimits(5, 0.1, 3);    //sigEta
            
            //////////////////////2d fitting///////////////////////////////////////
            
            fullySubtractedCorrCustomCentralitySideBand[k][i]->Fit(myfit, "qR0E");
            //fullySubtractedCorrCustomCentralitySideBandFit[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone("hfit");
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Eval(myfit);
            
            GaussIntegralParameter = calcScaleFactorforFit(myfit->GetParameter(5));
            
            cout << "Pt bin: " << k << "   Cent Bin: " << i << "   ProjectionNumber: " << GaussIntegralParameter << endl;
            
            fitParameterFile2D << endl << "pt bin (3 is integrated): " << k << "    Centrality Bin (0 - per, 1 - mid cent, 2 - cent): " << i << endl;
            fitParameterFile2D << endl << "A0" << "\t\t" << "AD" << "\t\t" <<"AQ" << "\t\t" <<"Jet Amp A1" << "\t\t" << "SigEta" << "\t\t" <<"SigPhi" << "\t\t" << endl;
            
            fitParameterFile2D << endl << myfit->GetParameter(0) << "\t" << myfit->GetParameter(1) << "\t" 
                                       << myfit->GetParameter(2) << "\t" << myfit->GetParameter(3) << "\t" 
                                       << myfit->GetParameter(5) << "\t" << myfit->GetParameter(4) << endl; 
            
            fitParameterFile2D << myfit->GetParError(0) << "\t" << myfit->GetParError(1) << "\t" 
                               << myfit->GetParError(2) << "\t" << myfit->GetParError(3) << "\t" 
                               << myfit->GetParError(5) << "\t" << myfit->GetParError(4) << endl;            
            
            fullySubtractedCorrCustomCentralitySideBandRes[k][i] = (TH2D*)fullySubtractedCorrCustomCentralitySideBand[k][i]->Clone("hres");
            //resHistUS[i]->Add(fitHistUS[i], angCorrUS[i], 1, -1);
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Add(fullySubtractedCorrCustomCentralitySideBandFit[k][i], -1);
            
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Draw("SURF1");
            c->SaveAs(str1);
 
            fullySubtractedCorrCustomCentralitySideBandFit[k][i]->Write();
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandRes[k][i]->Draw("SURF1");
            c->SaveAs(str1); 
            
            /////////////////////////////1d fitting///////////////////////////////////////////
    
            
    
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Fit(myfit1D, "qR0E");
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Clone("hfit");
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Eval(myfit1D);
 
            fitParameterFile1D << endl << "pt bin (3 is integrated): " << k << "    Centrality Bin (0 - per, 1 - mid cent, 2 - cent): " << i << endl;
            fitParameterFile1D << endl << "A0" << "\t\t" << "AD" << "\t\t" <<"AQ" << endl; // "\t\t" << "Jet Amp A1" << "\t\t" << "SigPhi" << endl;
            
            fitParameterFile1D << endl << myfit1D->GetParameter(0) << "\t" << myfit1D->GetParameter(1) << "\t" << myfit1D->GetParameter(2) << endl; //<< "\t" << myfit1D->GetParameter(3) << "\t" << myfit1D->GetParameter(5) << endl; 
            fitParameterFile1D << myfit1D->GetParError(0) << "\t" << myfit1D->GetParError(1) << "\t" << myfit1D->GetParError(2) << endl;//<< "\t" << myfit1D->GetParError(3) << "\t" << myfit1D->GetParError(5) << endl; 
            
            v2_2D = TMath::Sqrt(myfit->GetParameter(2));
            v2_2D_Error = v2_2D*(.5)*(myfit->GetParError(2)/myfit->GetParameter(2));
            
            v2_1D = TMath::Sqrt(myfit1D->GetParameter(2));
            v2_1D_Error = v2_2D*(.5)*(myfit1D->GetParError(2)/myfit1D->GetParameter(2));
            
            fitParameterFile2D << "V2 = " << v2_2D << " +/- " << v2_2D_Error << endl;
            fitParameterFile1D << "V2 = " << v2_1D << " +/- " << v2_1D_Error << endl;
            
            fitProjection->SetParameter(0, myfit->GetParameter(0));
            fitProjection->SetParameter(1, myfit->GetParameter(1));
            fitProjection->SetParameter(2, myfit->GetParameter(2));
            fitProjection->SetParameter(3, myfit->GetParameter(3));
            fitProjection->SetParameter(5, myfit->GetParameter(4));
            fitProjection->SetParameter(4, GaussIntegralParameter);
            
            quadrupole1->SetParameter(0, myfit->GetParameter(2));
            quadrupole2->SetParameter(0, myfit1D->GetParameter(2));
            quadrupole3->SetParameter(0, fitProjection->GetParameter(2));
            
            fitProjection->SetLineColor(2); //red
            quadrupole1->SetLineColor(3); // Green, from 2D fit
            quadrupole2->SetLineColor(4); //blue, from 1D fit
            myfit1D->SetLineColor(7); // light blue
            quadrupole3->SetLineColor(41);
            
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i] = (TH1D*)fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Clone("hres");
            //resHistUS[i]->Add(fitHistUS[i], angCorrUS[i], 1, -1);
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Add(myfit1D, -1);
            
            
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + fitLabel + fileType;
            //fullySubtractedCorrCustomCentralitySideBandFitPhiProj[k][i]->Draw();
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->Draw("EX0");
            fullySubtractedCorrCustomCentralitySideBandPhiProj[k][i]->SetMarkerStyle(20);
            
            fitProjection->Draw("SAME");
            quadrupole1->Draw("SAME");
            quadrupole2->Draw("SAME");
            //quadrupole3->Draw("SAME");
            myfit1D->Draw("SAME");
            c->SaveAs(str1);
 
            str1 = path + outputFolders[5] + outputFolders[8] + PtSubFolders[k] + phiProj + fullySubtractedLabelCustomCentrality + bandSubtract[1] + PtBinLabel[k] + binLabelPt[k] + centralityBin[i] + resLabel + fileType;
            fullySubtractedCorrCustomCentralitySideBandResPhiProj[k][i]->Draw();
            c->SaveAs(str1); 
 
 
        }
    }    
 
    fitParameterFile2D.close();
    fitParameterFile1D.close();
 
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //--------------------------------------Histogram formatting and output for talks----------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
 
    TCanvas *talkCanvas = new TCanvas("c3", "Histograms", 3300, 1700);
    talkCanvas->Divide(3,2);
    
    int padNumber = 1;
    TString talkOutput = "Histograms_for_talk";
    TString ptRanges[2] = {" 1 < p_{t} < 4 GeV ", " 4 < p_{t} < 20 GeV "};
    TString centRanges[3] = {" 50-100% ", " 20-50% ", " 0-20% "};
    
    double xValueforPave[3] = {1000, 2100, 3200};
    double yValueforPave[2] = {300, 1150};
    
    TString firstLine;
    TString secondLine;
    
    //PaveText *textBox = new TPaveText(
    
    for(int k = 1; k < 3; k++){
            for(int i = 0; i < 3; i++){
 
                
                formatCorrHist(fullySubtractedCorrCustomCentralitySideBand[k][i]);
                
                talkCanvas->cd(padNumber);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->SetStats(false);
                fullySubtractedCorrCustomCentralitySideBand[k][i]->SetTitle();
                fullySubtractedCorrCustomCentralitySideBand[k][i]->Draw("SURF1");
                
                TPaveText *textBox = new TPaveText(0.8, 0.8, 1.0, 1.0, "NDC");
                
                firstLine = ptRanges[k-1];
                secondLine = centRanges[i];
                
                textBox->AddText(firstLine.Data());
                textBox->AddText(secondLine.Data());
                textBox->Draw();
 
 
                padNumber++;
        }    
    }
 
    str1 = path + talkOutput + fileType;
    talkCanvas->SaveAs(str1);




















    file->Close();
    output->Close();
   
}

