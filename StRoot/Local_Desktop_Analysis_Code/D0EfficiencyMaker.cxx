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

Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if ( (reject && (x[0] > 1.81 && x[0] < 1.91)) || (reject && (x[0] > 1.6 && x[0] < 1.7)) ){     //|| (reject && (x[0] > 2.0 && x[0] < 2.1))
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}

int D0EfficiencyMaker(){
    
    TString path = "C:/Users/ajentsch/desktop/";
    TString rootFileName = "invMassHistograms.root";
    TString rootFileInput = path + rootFileName;
    
    TString invMassHistogramsFolder = "invMassHistograms/";

    TString fine = "FINE_";
    TString binLabelPtBin[11] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
    TString binLabelCentralityClass[3] = {"_Peripheral", "_MidCentral", "_Central"};
    TString bandBin[2] = {"D0_US_invMass_", "LS_invMass_"};
    TString ptLabel = "Pt_Bin_";
    TString final = "FullySubtractedD0InvMass_";
    TString finalHisto = "signalVsPt";
    TString LSSubtracted = "LS_Subtracted";
    TString USRaw = "US_Raw";
    TString SignalVsPt = "Signal/N_{ev} Vs. P_{t}";
    TString finalHistogramsLabel;
    
    TString ptRangeForPave;
    TString dash = "-";
    TString gev = " GeV/c";
    TString signal;
    TString counts = "S/N_{ev} = ";
    
    TString totalSignalForPave;
    
    TString outputFileName;
    TString pngExtenstion = ".png";

    int NUM_CENT_BINS = 3;
    int NUM_PT_BINS = 10;
    
    double sideBandLow      = 2.0;
    double sideBandHigh     = 2.1;
    double massLow          = 1.82;
    double massHigh         = 1.90;
    double countsPeakRegionUS = 0;
    double countsPeakRegionLS = 0;
    double countsSideBandUS   = 0;
    double countsSideBandLS   = 0;
    double sideBandScaleFactor = 1;
    double finalSignalCounts = 0;
    
    Int_t binMassLowUS;
    Int_t binMassHighUS;
    Int_t binMassLowLS;
    Int_t binMassHighLS;
    
    int numPerEvents = 0;
    int numMidCentEvents = 0;
    int numCentEvents = 0;
    
    double numEventsCent[3] = {0,0,0};
    
    double signalVsPtByCent[3][10];
    
    TAxis* xAxisUS;
    TAxis* xAxisLS;
    
    
    TString USHistString;
    TString LSHistString;
    
    ofstream yieldFile;
    TString yieldFileName = "D0_Yield_data.txt";
    TString yieldFileOutputFile = path + yieldFileName;
    yieldFile.open (yieldFileOutputFile);
    
    TH1D* USHistograms[3][10];
    TH1D* LSHistograms[3][10];
    TH1D* scaledLSHistograms[3][10];
    TH1D* LSSubtractedHistograms[3][10];
    TH1D* finalD0Histogams[3][10];
    
    TH1D* SVsPt[3];
    
    for(int i = 0; i < 3; i++){
    
        str1 = SignalVsPt + binLabelCentralityClass[i];   
        SVsPt[i] = new TH1D(str1, str1, 10, 0, 10);
        SVsPt[i]->SetMarkerStyle(20);
        
    }    
    
    TString str = "event counts per Vz/Cent bin";
    
    TFile *inputRootFile = new TFile(rootFileInput);
    
    TKey *key = inputRootFile->GetKey(str);
    if (key) {
            TH1F *EventCounter = key->ReadObjectAny(TH2I::Class());
    }
    
    //inputRootFile.ls();
    
    for(int mult = 1; mult < 3; mult++){
        for(int vz = 1; vz < 11; vz++){
            
            numEventsCent[0] = numEventsCent[0] + EventCounter->GetBinContent(mult, vz);
        }
    }

    for(int mult = 4; mult < 9; mult++){
        for(int vz = 1; vz < 11; vz++){
            
            numEventsCent[1] = numEventsCent[1] + EventCounter->GetBinContent(mult, vz);
        }
    }
    
    for(int mult = 10; mult < 17; mult++){
        for(int vz = 1; vz < 11; vz++){
            
            numEventsCent[2] = numEventsCent[2] + EventCounter->GetBinContent(mult, vz);
        }
    }
    
    cout << endl;
    cout << "Peripheral Events: " << numEventsCent[0] << endl;
    cout << "Mid-Central Events: " << numEventsCent[1] << endl;
    cout << "Central Events: " << numEventsCent[2] << endl;
    
    int TotalEvents = numEventsCent[0] + numEventsCent[1] + numEventsCent[2];
    
    cout << "Total : " << TotalEvents << endl << endl;
    
    TCanvas *outputCanvas = new TCanvas("c1", "c1", 3500, 1400);
    outputCanvas->Divide(5,2);
    
    for(int i = 0; i < NUM_CENT_BINS; i++){
        
        //yieldFile << "Centrality Bin: " << binLabelCentralityClass[i] << endl << endl;
         
        for(int k = 0; k < NUM_PT_BINS; k++){

            outputCanvas->cd(k+1);
        
            TPaveText *textBox = new TPaveText(0.6, 0.75, 0.9, 0.9, "NDC");
        
            //yieldFile << "Pt Bin: " << binLabelPtBin[k] << "----->  ";
            //yieldFile << "Pt Bin: " << binLabelPtBin[k] << "  S/N_ev = ";
            
            USHistString = fine + bandBin[0] + ptLabel + binLabelPtBin[k] + binLabelCentralityClass[i];
            LSHistString = fine + bandBin[1] + ptLabel + binLabelPtBin[k] + binLabelCentralityClass[i];
            finalHistogramsLabel = final + ptLabel + binLabelPtBin[k] + binLabelCentralityClass[i];
            
            //yieldFile << USHistString << "\t" << LSHistString << endl;
            
            USHistograms[i][k] = new TH1D("","", 50, 1.6, 2.1);
            LSHistograms[i][k] = new TH1D("","", 50, 1.6, 2.1);
            scaledLSHistograms[i][k] = new TH1D("","", 50, 1.6, 2.1);
            LSSubtractedHistograms[i][k] = new TH1D("","", 50, 1.6, 2.1);
            finalD0Histogams[i][k] = new TH1D("","", 50, 1.6, 2.1);
            
            USHistograms[i][k] = (TH1D*) inputRootFile->Get(USHistString);
            LSHistograms[i][k] = (TH1D*) inputRootFile->Get(LSHistString);
            USHistograms[i][k]->Sumw2();
            USHistograms[i][k]->SetMarkerStyle(20);
            LSHistograms[i][k]->Sumw2();
            LSHistograms[i][k]->SetMarkerStyle(20);
            
            scaledLSHistograms[i][k] = (TH1D*)LSHistograms[i][k]->Clone();
            
            xAxisUS = USHistograms[i][k]->GetXaxis();
            xAxisLS = LSHistograms[i][k]->GetXaxis();
        
            binMassLowUS  = USHistograms[i][k]->GetXaxis()->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
            binMassHighUS = USHistograms[i][k]->GetXaxis()->FindBin(sideBandHigh);
            binMassLowLS = LSHistograms[i][k]->GetXaxis()->FindBin(sideBandLow);            //get normalization factors to scale the LS distribution
            binMassHighLS = LSHistograms[i][k]->GetXaxis()->FindBin(sideBandHigh);
            
            countsSideBandUS = USHistograms[i][k]->Integral(binMassLowUS, binMassHighUS);
            countsSideBandLS = LSHistograms[i][k]->Integral(binMassLowLS, binMassHighLS);
            
            sideBandScaleFactor = countsSideBandUS/countsSideBandLS;
            
            binMassLowUS  = xAxisUS->FindBin(massLow);            
            binMassHighUS = xAxisUS->FindBin(massHigh);
            binMassLowLS  = xAxisLS->FindBin(massLow);            
            binMassHighLS = xAxisLS->FindBin(massHigh);
            
            
            
            scaledLSHistograms[i][k]->Scale(sideBandScaleFactor);
            
            countsPeakRegionUS = USHistograms[i][k]->Integral(binMassLowUS, binMassHighUS);      //S+B
            
            
            
            LSSubtractedHistograms[i][k]->Add(USHistograms[i][k], scaledLSHistograms[i][k], 1, -1);
            
            
            
            finalD0Histogams[i][k] = (TH1D*)LSSubtractedHistograms[i][k]->Clone();
            
            TF1 *fl = new TF1("fl",fline,1.6,2.1,2);
            fl->SetParameters(2,-1);
            //fit only the linear background excluding the signal area
            reject = kTRUE;
            LSSubtractedHistograms[i][k]->Fit(fl,"qR");
            reject = kFALSE;
        
            finalD0Histogams[i][k]->Add(fl, -1);
            
            finalSignalCounts = finalD0Histogams[i][k]->Integral(binMassLowUS, binMassHighUS);
            finalD0Histogams[i][k]->SetStats(0);
            finalD0Histogams[i][k]->SetTitle(finalHistogramsLabel);
            finalD0Histogams[i][k]->SetMarkerStyle(20);
            finalD0Histogams[i][k]->Draw();
            
            ptRangeForPave = binLabelPtBin[k] + dash + binLabelPtBin[k+1] + gev;
            finalSignalCounts = finalSignalCounts/numEventsCent[i];
            
            if(finalSignalCounts > 0) { 
            
                SVsPt[i]->SetBinContent(k+1, finalSignalCounts);
                signalVsPtByCent[i][k] = finalSignalCounts;
            }
            
            else signalVsPtByCent[i][k] = 0;
            
            cout << finalSignalCounts << endl;
            signal.Form("%.9f", finalSignalCounts);
   
            //yieldFile << signal << endl;
            
            totalSignalForPave = counts + signal;
            
            textBox->AddText(ptRangeForPave);
            textBox->AddText(totalSignalForPave);
            textBox->GetLine(0)->SetTextSize(.055);
            textBox->GetLine(1)->SetTextSize(.055);
            
            textBox->Draw("SAME");
            
        }
        
        yieldFile << endl;
        outputFileName = path + invMassHistogramsFolder + final + binLabelCentralityClass[i] + pngExtenstion;
        outputCanvas->SaveAs(outputFileName);
        
       
    }
    
    /*for(int i = 0; i < NUM_CENT_BINS; i++){
        for(int k = 0; k < NUM_PT_BINS; k++){
            
            outputCanvas->cd(k+1);
            
            TPaveText *textBox = new TPaveText(0.6, 0.75, 0.9, 0.9, "NDC");
            
            LSSubtractedHistograms[i][k]->SetStats(0);
            LSSubtractedHistograms[i][k]->Draw();
            
            
            ptRangeForPave = binLabelPtBin[k] + dash + binLabelPtBin[k+1] + gev;
            textBox->AddText(ptRangeForPave);
            textBox->GetLine(0)->SetTextSize(.055);
            
            
            textBox->Draw("SAME");
            outputFileName = path + invMassHistogramsFolder + LSSubtracted + binLabelCentralityClass[i] + pngExtenstion;
            outputCanvas->SaveAs(outputFileName);
            
        }
    }*/
    
    /*for(int i = 0; i < NUM_CENT_BINS; i++){
        for(int k = 0; k < NUM_PT_BINS; k++){
            
            outputCanvas->cd(k+1);
            
            TPaveText *textBox = new TPaveText(0.6, 0.75, 0.9, 0.9, "NDC");
            
            USHistograms[i][k]->SetStats(0);
            USHistograms[i][k]->Draw();
            
            
            ptRangeForPave = binLabelPtBin[k] + dash + binLabelPtBin[k+1] + gev;
            textBox->AddText(ptRangeForPave);
            textBox->GetLine(0)->SetTextSize(.055);
            
            
            textBox->Draw("SAME");
            outputFileName = path + invMassHistogramsFolder + USRaw + binLabelCentralityClass[i] + pngExtenstion;
            outputCanvas->SaveAs(outputFileName);
            
        }
    }*/
    
    TCanvas *finalCanvas = new TCanvas("c2", "c2", 3300, 850);
    finalCanvas->Divide(3,1);
    
    for(int i = 0; i < 3; i++){
    
        finalCanvas->cd(i+1); 
        //finalCanvas->SetGrid();        
        SVsPt[i]->Draw("P");
        //SVsPt[i]->GetXaxis()->SetAxisColor(17);
        //SVsPt[i]->GetYaxis()->SetAxisColor(17);
        //finalCanvas->RedrawAxis();
        
    }
        outputFileName = path + invMassHistogramsFolder + finalHisto + pngExtenstion;
        finalCanvas->SaveAs(outputFileName);
        
        
    //////////////////////NOW WE WILL CALCULATE EFFICIENCIES/////////////////////////////////////////

    double m0 = 1.86; //GeV
    
    TF1 *levyPeripheral = new TF1("LevyPer", "((x*x)/(TMath::Sqrt((x*x)+(1.86*1.86))))*([0]/((1+((TMath::Sqrt((x*x)+(1.86*1.86))-1.86)/([1]*[2])))**[1]))", 0, 10);
    TF1 *levyMidCentral = new TF1("LevyMidCent", "((x*x)/(TMath::Sqrt((x*x)+(1.86*1.86))))*([0]/((1+((TMath::Sqrt((x*x)+(1.86*1.86))-1.86)/([1]*[2])))**[1]))", 0, 10);
    TF1 *levyCentral = new TF1("LevyCent", "((x*x)/(TMath::Sqrt((x*x)+(1.86*1.86))))*([0]/((1+((TMath::Sqrt((x*x)+(1.86*1.86))-1.86)/([1]*[2])))**[1]))", 0, 10); 

    
    
    levyPeripheral->SetParameter(0, .0214399993); // A  
    levyPeripheral->SetParameter(1, 12.99); // n
    levyPeripheral->SetParameter(2, .3256); //T
    levyMidCentral->SetParameter(0, .15630001); // A  
    levyMidCentral->SetParameter(1, 212.0); // n
    levyMidCentral->SetParameter(2, .37059999); //T
    levyCentral->SetParameter(0, .36400002); // A  
    levyCentral->SetParameter(1, 15.91); // n
    levyCentral->SetParameter(2, .30699998); //T
    
    /*TCanvas *functionCanvas = new TCanvas("f1", "f1", 1100, 850);
    
    levyPeripheral->Draw();
    
    TString func = "func";
    
    str1 = path + invMassHistogramsFolder + func + pngExtenstion;
    functionCanvas->SaveAs(str1);*/
    
    double levyPerCentrality[3];
    //double levyPerEval;
    //double levyMidEval;
    //double levyCentEval;
    double efficiency;
    double pt;
    
    TH1D* effPlot[3];
    
    TString effLabel = "Efficiency";
    
    for(int i = 0; i < 3; i++){
        
        str1 = effLabel + binLabelCentralityClass[i];
        effPlot[i] = new TH1D(str1, str1, 10, 0, 10);
        effPlot[i]->SetMarkerStyle(20);
    }
    
    
    
    yieldFile << "Centrality Bin: " << binLabelCentralityClass[0] << endl << endl;
    for(int k = 0; k < 10; k++){
    
        pt = k + 0.5;
        
        levyPerCentrality[0] = levyPeripheral->Eval(pt, 0, 0, 0);
        efficiency = signalVsPtByCent[0][k]/levyPerCentrality[0];
        
        yieldFile << "Pt Bin: " << binLabelPtBin[k] << "  Eff = " << efficiency << endl;
        
        if(efficiency > 0) { effPlot[0]->SetBinContent(k+1, efficiency); }
    }
    
    yieldFile << "Centrality Bin: " << binLabelCentralityClass[1] << endl << endl;
         
    for(int k = 0; k < 10; k++){
    
        pt = k + 0.5;
       
        levyPerCentrality[1] = levyMidCentral->Eval(pt, 0, 0, 0);          
        efficiency = signalVsPtByCent[1][k]/levyPerCentrality[1];
        
        yieldFile << "Pt Bin: " << binLabelPtBin[k] << "  Eff = " << efficiency << endl;
        
        if(efficiency > 0) { effPlot[1]->SetBinContent(k+1, efficiency); }
    }   
    
    yieldFile << "Centrality Bin: " << binLabelCentralityClass[2] << endl << endl;
    
    for(int k = 0; k < 10; k++){
    
        pt = k + 0.5;
       
        levyPerCentrality[2] = levyCentral->Eval(pt, 0, 0, 0);          
        efficiency = signalVsPtByCent[2][k]/levyPerCentrality[2];
        
        yieldFile << "Pt Bin: " << binLabelPtBin[k] << "  Eff = " << efficiency << endl;
        
        if(efficiency > 0) { effPlot[2]->SetBinContent(k+1, efficiency); }
    }   
    
    TCanvas *effCanvas = new TCanvas("c2", "c2", 2200, 1700);
    //effCanvas->Divide(2,2);
    
    
    for(int i = 0; i < 3; i++){
    
        //effCanvas->cd(i+1);
        effCanvas->SetLogy();
        effPlot[i]->Draw("P");
        
        str1 = path + invMassHistogramsFolder + effLabel + binLabelCentralityClass[i] + pngExtenstion;
        effCanvas->SaveAs(str1);
        
    }
    
    //str1 = path + invMassHistogramsFolder + effLabel + pngExtenstion;
    //effCanvas->SaveAs(str1);
    
    
    yieldFile.close();
    
}