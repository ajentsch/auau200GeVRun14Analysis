#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"

#include "StPicoCharmContainers/StPicoD0Event.h"
#include "StPicoCharmContainers/StKaonPion.h"

#include "StMixedEventBuffer/StMixedEventBuffer.h"
#include "StMixedEventBuffer/StMixerEvent.h"
#include "StMixedEventBuffer/StMixerTrack.h"

#include "StPicoD0AnaMaker.h"
#include "StPicoHFMaker/StHFCuts.h"


/****

Author: Alex Jentsch

Current version of AnaMaker (by date): 2/22/2016

Description of current functionality:

1. Can read in both Trees and PicoDst -> Stores all basic QA info for PicoDsts
2. Can make invMass Histos for various pt bins
3. Can produce sibling histograms for various pt bins


Still Needs:

1. Various BG production methods (side band, mixing, etc.)


Update (3/10/2016)

Currently adding event mixing to code. Using the StPicoMixedEventMaker as a baseline, but changing its implementation significantly.
I want to be able to handle all of the mixing in my AnaMaker with all of the mixing functions being used as if from any other normal class.


Update (3/31/2016)

Buffering code built and seems to function. Lots of testing left to do. Need to make sure the stored event information can be retrieved 
and can be binned into histograms.


Update (4/19/2016)

Event mixing is working and the buffering is not causing any problems. May be changing the event mixing algorithm to better the statistics.

Update (6/3/2016)

Code is fully functional for both non-identified and D0-Hadron correlations. Still working on getting all of the cuts right.
Code has been updated for the new code structure with the updated production.

Update (6/28/2016)

Hadron-Hadron Correlation code has been removed from this version. The event mixer has been simplified and completely redone. Lots of variables were renamed.

****/

ClassImp(StPicoD0AnaMaker)

StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
                                   char const * outName, StPicoDstMaker* picoDstMaker): 
                            StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL),
                            mOutFileName(outName), mInputFileList(inputFilesList),
                            mOutputFile(NULL), mChain(NULL), mEventCounter(0), mHFCuts(NULL)
{}

//---------------------------------------Initialization-----------------------------------------
Int_t StPicoD0AnaMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFileList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoD0AnaMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
   mOutputFile->cd();

   if (!mHFCuts)
    mHFCuts = new StHFCuts;   
   mHFCuts->init();
   
   
   //-----------------------FLAGS--------------------------------//
   DEBUG               = false;
   DEBUG_MIX_BUFFER    = false;
   USE_CENT_BINS       = true;
   USE_VZ_BINS         = true;
   USE_PT_BINS         = true;
   D0_HADRON_CORR      = true;
   EVENT_MIXING        = true;
   //----------------------------------------------------------//
   
   //---------------------Important constants-----------------//
   BUFFER_SIZE   = 5;
   NUM_PHI_BINS  = 12;
   NUM_ETA_BINS  = 9;
   NUM_VZ_BINS   = 10;
   NUM_CENT_BINS = 15; 
   NUM_PT_BINS   = 3;
   
   
   // --------------------Event Mixer Buffer-------------------------------------
   
    if(USE_VZ_BINS) { nVzBins = NUM_VZ_BINS; }  //flag to set binning on Vz
    else nVzBins = 1;
    
    if(USE_CENT_BINS) { nCentBins = NUM_CENT_BINS; }
    else nCentBins = 1;
    
    //********************D0_Hadron section************************//
    
    if(D0_HADRON_CORR){
    
        for(int i = 0; i < nCentBins; i++){
            for( int j = 0; j < nVzBins; j++){
        
                eventBufferPicoEvent[j][i] = new StMixedEventBuffer();
                eventBufferPicoEvent[j][i]->setBufferSize(BUFFER_SIZE);    //set buffer size here -- the amount of events in each 2d bin
                eventBufferPicoEvent[j][i]->setBufferCounter(0);
            }
        }
    }        
   
    if(D0_HADRON_CORR){
   
        for(int i = 0; i < nCentBins; i++){
            for( int j = 0; j < nVzBins; j++){
        
                eventBufferD0Candidate[j][i] = new StMixedEventBuffer();
                eventBufferD0Candidate[j][i]->setBufferSize(BUFFER_SIZE);    //set buffer size here -- the amount of events in each 2d bin
                eventBufferD0Candidate[j][i]->setBufferCounter(0);
            }
        }
    }       
   
    //--------------------CUTS------------------------------------------
   
    // D0 Cuts
   
    kaonPtCut           = .15;  //daughter pt cut
    pionPtCut           = .15;  //daughter pt cut
    
    D0InvMassLow        = 1.82; //signal region US invariant mass
    D0InvMassHigh       = 1.90; //signal region US invariant mass
    
    USSideBandLeftLow   = 1.62;
    USSideBandLeftHigh  = 1.70;
    USSideBandRightLow  = 2.0;
    USSideBandRightHigh = 2.1;
    
    d0PtLow             = 0.0;
    d0PtHigh            = 20.0;
    d0DecayLengthMin    = .0220;
    d0DecayLengthMax    = 999999.0;
    daughterDCA         = .0055;
    d0DaughterPionPtMin = .15;
    d0DaughterKaonPtMin = .15;
    kaonDCA             = .008;
    pionDCA             = .008;
    d0DCAtoPV           = .0065;
    
    //hadron cuts
    
    hadronPtMin         = .15;
    hadronPtMax         = 15.45;
    
    trackChi2max        = 3.0;
    trackDCAtoPvtx      = 3.0;
    
	//Run 14 MinBias triggers
    trigger[0] = 450050;
    trigger[1] = 450060;
    trigger[2] = 450005;
    trigger[3] = 450015;
    trigger[4] = 450025;
    //-----------------------------------------------------------------
    
    eventNumber = 1;
    
    TString VzBinLabel   = "_VzBin_";
    TString CentBinLabel = "_CentBin_";
    TString PtBinLabel   = "_PtBin_";
    
    //Labels for D0-hadron US and LS corr histograms
    
    TString SibCorrLabels[4] = {"Sibling_SideBandLeft_correlation", "Sibling_US_correlation", "Sibling_SideBandRight_correlation", "Sibling_LS_correlation"};
    TString MixCorrLabels[4] = {"Mixed_SideBandLeft_correlation", "Mixed_US_correlation", "Mixed_SideBandRight_correlation", "Mixed_LS_correlation"};
    
    //other labels
   
    TString eventCounterLabel = "Event Count Vz ";
    TString etaLabel          = "Inclusive Hadron Eta";
    TString phiLabel          = "Inclusive Hadron Phi";
    TString etaPhiLabel       = "Inclusive 2D Hadron Eta/Phi";
    
    TString invMassD0PtBin    = "D0_US_invMass_Pt_Bin_";
    TString invMassLSPtBin    = "LS_invMass_Pt_Bin_";
    
    TString str1;
    TString str2;
    TString str3;
    
    TString binLabelVz[10]   = {"0", "1", "2", "3", "4", "5", "6", "7", "8","9"};
    TString binLabelCent[15] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"};
    TString binLabelPt[6]    = {"0", "1", "2", "3", "4", "5"};
   
   // --------------------Begin User Variables-----------------------------------
   
   //ptDist          = new TH1D("Pt Distribution", "Pt Distribution", 1000, 0, 10);              
   invMass         = new TH1D("unlikeSign", "unlikeSign", 50, 1.6, 2.1);
   likeSignBG      = new TH1D("LikeSignBG", "LikeSignBG", 50, 1.6, 2.1);
   
   
   if(USE_PT_BINS){
        for(int k = 0; k < NUM_PT_BINS; k++){
        
            str1 =  invMassD0PtBin + binLabelPt[k];       
            D0InvMassPtBin[k] = new TH1D(str1, str1, 50, 1.6, 2.1);
            str1 =  invMassLSPtBin + binLabelPt[k];       
            LSInvMassPtBin[k] = new TH1D(str1, str1, 50, 1.6, 2.1);
            
        }
    }        
   
   
   if(D0_HADRON_CORR){//begin D0-Hadron Correlation Histograms
        for(int band = 0; band < 4; band++){
            for(int i = 0; i < nVzBins; i++){ //Initialize all of the histograms for storing the sibling and mixed information
                for(int j = 0; j < nCentBins; j++){
        
                    str1 = SibCorrLabels[band] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                    str2 = MixCorrLabels[band] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                
                    sibCorrBin[band][i][j]             = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                    mixCorrBin[band][i][j]             = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
               
                }
            }
        } 
    }        
    
    if(D0_HADRON_CORR && USE_PT_BINS){//begin D0-Hadron Correlation Histograms -- USE PT BINS HERE
        for(int band = 0; band < 4; band++){
            for(int i = 0; i < nVzBins; i++){ 
                for(int j = 0; j < nCentBins; j++){
                    for(int k = 0; k < NUM_PT_BINS; k++){
                
                
                        str1 = SibCorrLabels[band] + PtBinLabel + binLabelPt[k] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                        str2 = MixCorrLabels[band] + PtBinLabel + binLabelPt[k] + VzBinLabel + binLabelVz[i] + CentBinLabel + binLabelCent[j];
                    
                        sibCorrBinPt[band][k][i][j] = new TH2D(str1, str1, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                        mixCorrBinPt[band][k][i][j] = new TH2D(str2, str2, NUM_ETA_BINS, -2, 2, NUM_PHI_BINS, -TMath::PiOver2(), 3*TMath::PiOver2());
                    
                    }
                }            
            }
        }    
    }//end D0-Hadron Correlation Histograms -- USE PT BINS HERE    
   
    if(USE_VZ_BINS){
        for(int i = 0; i < nVzBins; i++){ 
            
            str1 = phiLabel + VzBinLabel + binLabelVz[i];
            str2 = etaLabel + VzBinLabel + binLabelVz[i];
            str3 = etaPhiLabel + VzBinLabel + binLabelVz[i];
            phiDistVz[i]    = new TH1D(str1, str1, 48, -TMath::Pi(), TMath::Pi());
            etaDistVz[i]    = new TH1D(str2, str2, 50, -1, 1);
            etaPhiDistVz[i] = new TH2D(str3, str3, 25, -1, 1, 24, -TMath::Pi(), TMath::Pi());
            etaPhiDistVz[i]->GetXaxis()->SetTitle("#eta");
            etaPhiDistVz[i]->GetYaxis()->SetTitle("#phi");
        }    
    }
 
   
   //QA Histograms
   eventCounter    = new TH1I("number of events used", "number of events used", 4, 0, 4);
   trackCounter    = new TH1I("number of tracks per event", "number of tracks per event", 1500, 0, 1499);
   usedTracks      = new TH1I("Used tracks", "Used tracks", 1500, 0, 1499);
   kaonEtaDist     = new TH1D("Kaon Eta Distribution", "Kaon Eta Distribution", 250, -1, 1);
   pionEtaDist     = new TH1D("Pion Eta Distribution", "Pion Eta Distribution", 250, -1, 1);
   kaonPhiDist     = new TH1D("Kaon Phi Distribution", "Kaon Phi Distribution", 250, -2*TMath::Pi(), 2*TMath::Pi());
   pionPhiDist     = new TH1D("Pion Phi Distribution", "Pion Phi Distribution", 250, -2*TMath::Pi(), 2*TMath::Pi());
   DCAtoPrimaryVertex    = new TH1D("track DCA to PV", "track DCA to PV", 500, 0.0, 10.0);
   DCAtoPrimaryVertexCut = new TH1D("track DCA to PV cut check", "track DCA to PV cut check", 500, 0.0, 10.0);
   hadronPtDist    = new TH1D("Inclusive Hadron pt", "Inclusive Hadron pt", 250, 0, 10);
   hadronPhiDist   = new TH1D("Inclusive Hadron Phi", "Inclusive Hadron Phi", 250, -TMath::Pi(), TMath::Pi());
   hadronEtaDist   = new TH1D("Inclusvie Hadron Eta", "Inclusive Hadron Eta", 250, -1, 1);
   hadronChi2      = new TH1D("Chi2 for hadron tracks", "Chi2 for hadron tracks", 500, 0, 10);
   pVtxX           = new TH1D("X position of pVtx", "X position of pVtx", 500, -10, 10);
   pVtxY           = new TH1D("Y position of pVtx", "Y position of pVtx", 500, -10, 10);
   pVtxZ           = new TH1D("Z position of pVtx", "Z position of pVtx", 500, -7, 7);
   dEdxVsPt        = new TH2D("dEdx_vs_P", "dEdx_vs_P", 250, 0, 10, 250, 0, 10);
   invBetaVsPt     = new TH2D("#Beta^{-1} Vs. P", "#Beta^{-1} Vs. P", 250, 0, 10, 250, 0, 4);
   vZandCentBinPerEvent = new TH2I("event counts per Vz/Cent bin", "event counts per Vz/Cent bin", 15, 0, 15, 10, 0, 10);
   vZandCentBinPerEvent->GetXaxis()->SetTitle("Centrality Bin");
   vZandCentBinPerEvent->GetYaxis()->SetTitle("Vz Bin");
   //QA for mass-cut D0  
   D0ptDist        = new TH1D("D0 Candidate pt Dist (mass cut)", "D0 Candidate pt Dist (mass cut)", 250, 0, 10);
   D0EtaDist       = new TH1D("D0 Eta Dist", "D0 #eta Dist. (mass cut)", 250, -1, 1);
   D0PhiDist       = new TH1D("D0 Phi Dist", "D0 #phi Dist. (mass cut)", 250, -2*TMath::Pi(), 2*TMath::Pi());
   d0CountPerEvent = new TH1I("number of D0 candidates per event", "number of D0 candidates per event", 50, 0, 50);
   histOfCuts      = new TH1D("HistOfCuts", "HistOfCuts", 30, 1, 30);  
   
   
   //Histogram formatting

   eventCounter->GetXaxis()->SetBinLabel(1,"minBias+D0 Cand.");
   eventCounter->GetXaxis()->SetBinLabel(2,"minBias");
   eventCounter->GetXaxis()->SetBinLabel(3,"total");
   eventCounter->GetXaxis()->SetBinLabel(4,"events from bad runs");
   
//----------------------End User Variables------------------------------------
  
   //Fill Cuts histogram here for producing output text file
    histOfCuts->SetBinContent(1, D0InvMassLow); //D0InvMassLow
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
    histOfCuts->SetBinContent(17, hadronPtMin); //hadronPtMin
    histOfCuts->SetBinContent(18, hadronPtMax); //hadronPtMax
    histOfCuts->SetBinContent(19, BUFFER_SIZE); //BUFFER_SIZE
    histOfCuts->SetBinContent(20, NUM_PHI_BINS); //NUM_PHI_BINS
    histOfCuts->SetBinContent(21, NUM_ETA_BINS); //NUM_ETA_BINS
    histOfCuts->SetBinContent(22, 0); // Require HFT
    histOfCuts->SetBinContent(23, trackChi2max); // chi2 cut
    histOfCuts->SetBinContent(24, trackDCAtoPvtx); // DCA cut
    histOfCuts->SetBinContent(29, 1); //Scale factor counter to produce text file
    
    
   
   return kStOK;
}

//-----------------------------------------------------------------------------

//------------------------------------Destructor------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
   /*  */
}
//-----------------------------------------------------------------------------

//------------------------------------Finish-----------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
   LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
   mOutputFile->cd();
   // save user variables here
   
   if(D0_HADRON_CORR){
        for(int band = 0; band < 4; band++){
            for(int i = 0; i < nVzBins; i++){
                for(int j = 0; j < nCentBins; j++){
        
                    sibCorrBin[band][i][j]->Write();
                    mixCorrBin[band][i][j]->Write();
           
                    for(int k = 0; k < NUM_PT_BINS; k++){
            
                        sibCorrBinPt[band][k][i][j]->Write();
                        mixCorrBinPt[band][k][i][j]->Write();
                    }
                }    
            }
        }
    }        
    if(D0_HADRON_CORR){
        for(int i = 0; i < nVzBins; i++){
            for(int j = 0; j < nCentBins; j++){
                cout << "deleting buffer[" << i << "][" << j << "]" << endl;
                delete eventBufferPicoEvent[i][j];
                delete eventBufferD0Candidate[i][j];
            }
        }
    }
    
    if(USE_VZ_BINS){
        for(int i = 0; i < nVzBins; i++){
           
                etaDistVz[i]->Write();
                phiDistVz[i]->Write();
                etaPhiDistVz[i]->Write();
        }
    }      


   if(USE_PT_BINS){
        for(int i = 0; i < NUM_PT_BINS; i++){
            
            D0InvMassPtBin[i]->Write();
            LSInvMassPtBin[i]->Write();
        }
   }        
  
   
   invMass->Write();
   likeSignBG->Write();
   D0EtaDist->Write();
   D0PhiDist->Write();
   D0ptDist->Write(); 
   eventCounter->Write();
   kaonEtaDist->Write();
   pionEtaDist->Write();
   kaonPhiDist->Write();
   pionPhiDist->Write();
   hadronPtDist->Write();
   hadronPhiDist->Write();
   hadronEtaDist->Write();
   trackCounter->Write();
   dEdxVsPt->Write();
   invBetaVsPt->Write();
   DCAtoPrimaryVertex->Write();
   DCAtoPrimaryVertexCut->Write();
   d0CountPerEvent->Write();
   vZandCentBinPerEvent->Write();
   usedTracks->Write();
   histOfCuts->Write();
   pVtxX->Write();
   pVtxY->Write();
   pVtxZ->Write();
   hadronChi2->Write();
   
   mOutputFile->Close();
  

   return kStOK;
}
//-----------------------------------------------------------------------------

//----------------------------------Make---------------------------------------
Int_t StPicoD0AnaMaker::Make(){ //begin Make member function
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();


   if (!picoDst)
   {
      LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if(mPicoD0Event->runId() != picoDst->event()->runId() ||
       mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
     LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
     LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
     exit(1);
   }

   //-------------------Begin User Analysis---------------------------
 
/****************************************************************************************************************/
/****************************GLOBAL AND LOCAL VARIABLES SET AND INITIALIZED BEGIN*********************************/  
/***************************************************************************************************************/ 
    
    
    //-------------- Various global cuts and pt ranges---------------- 
  
    trackCount = 0;
    centralityBin = 0;  //This number will be between 0 and 8 -- 9 bins total
    
    VzBin = 0;
    bandBin = -1; 
	
    //--------------Local Variables used in the code-------------------
    double delPhi         = 0;
    double delPhiCp       = 0;
    double delEta         = 0;
    double pt             = 0;
    double phi            = 0;
    double eta            = 0;
    double dEdx           = 0;
    double beta           = 0;
    int    PIDflag        = 0;
    int    numD0s         = 0;
    bool   minBiasFlag    = false;
    int    realTracks     = 0;
    
    StPicoTrack* trk;
    //StPicoTrack* trk2;
    StThreeVectorF trackMom;
    StThreeVectorF trackMom2;
    double bField         = picoDst->event()->bField();
    StThreeVectorF pVtx   = picoDst->event()->primaryVertex();
    StThreeVectorF kaonPionMom;
    
    std::vector <StThreeVectorF> mAssociatedHadronList;
    
    bool storeEventToMix        = true;
    bool storeKaonPionEvent     = false;
    int d0Counter = 0;
    
/****************************************************************************************************************/
/****************************GLOBAL AND LOCAL VARIABLES SET AND INITIALIZED END********************************/  
/***************************************************************************************************************/   
    
 
/****************************************************************************************************************/
/****************************PRELIMINARY EVENT CUTS BEGIN********************************************************/  
/***************************************************************************************************************/    
 
    if(!mHFCuts->isGoodRun(picoDst->event())){
      
        eventCounter->Fill(2);
        eventCounter->Fill(3);
        return kStOk; //makes sure the event comes from good run
    }

    eventCounter->Fill(2);
    
    for(unsigned int i = 0; i < 5; i++){
    
        if(picoDst->event()->isTrigger(trigger[i])){
        
            minBiasFlag = true;
            break;
        }
    }        
    
    if(!minBiasFlag) {return kStOK;}
    
    eventCounter->Fill(1);

    TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();    //The Kaon-Pion list is generated here

    trackCount = picoDst->numberOfTracks();            
    trackCounter->Fill(picoDst->numberOfTracks());
   
/****************************************************************************************************************/
/****************************PRELIMINARY EVENT CUTS END********************************************************/  
/***************************************************************************************************************/ 
  
   
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and tracks BEGIN**********/
/****************************************************************************************************************/   
   
    for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){ // Begin loop to fill basic track information and make nega/posi list
                                                                 //gets pt, eta, phi, dEdx and beta from TOF
        trk = picoDst->track(i);
        trackMom = trk->gMom(pVtx, bField);
        
        if(!mHFCuts->isGoodTrack(trk)) { continue; }    //Checks the HFT requirement for a track
        
        trackDCA = ((trk->helix().origin())-pVtx).mag();
        DCAtoPrimaryVertex->Fill(trackDCA);            //QA plot to see full DCA distribution
        
        if(trk->chi2() > trackChi2max) { continue; }  //check chi2 cut
        if(!checkDCAtoPV(trackDCA))  { continue; }   // track quality cut for DCA to PV
        
        pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
        if(pt < hadronPtMin || pt > hadronPtMax){continue;}              //check pt cut
        
        if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }  //remove electron contamination
        if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
        
        eta  = trackMom.pseudoRapidity();
		if(eta > 1 || eta < -1) { continue; }
		
        realTracks++;                        //after a track passes all cuts, it is considered a "real track"
        
        mAssociatedHadronList.push_back(trackMom); //store the associated tracks to a list for quicker pairing later
        
        phi  = TMath::ATan2(trackMom.y(),trackMom.x());  
        dEdx = trk->dEdx();
        hadronChi2->Fill(trk->chi2());   
        DCAtoPrimaryVertexCut->Fill(trackDCA);        
            
        if(mHFCuts->hasTofPid(trk)){    
         
            beta = mHFCuts->getTofBeta(trk);                                                
            invBetaVsPt->Fill(trackMom.mag(), (1/beta));        
        }
        
        dEdxVsPt->Fill(trackMom.mag(), dEdx);
        hadronPtDist->Fill(pt);                                                
        hadronEtaDist->Fill(eta);
        etaDistVz[VzBin]->Fill(eta);
        phiDistVz[VzBin]->Fill(phi);
        etaPhiDistVz[VzBin]->Fill(eta,phi);
        
        if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaKaon()) < 2.0)){       //ONLY check nSigma for TPC track -- used to be mHFCuts->isTPCKaon(trk), but this include pt cut
                                                                                  //need to fix this. Need to figure out how to access the nSigma from the cuts
              kaonEtaDist->Fill(eta);
              kaonPhiDist->Fill(phi);
              continue;
        }

        if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaPion()) < 3.0)){     //ONLY check nSigma for TPC track
               
              pionEtaDist->Fill(eta);
              pionPhiDist->Fill(phi);
              continue;
        }
        
    } //End loop to fill basic track information and make nega/posi list   
       
    usedTracks->Fill(realTracks);   
    
    centralityBin = getCentralityBin(realTracks);  //get centrality bin -- These mult bins are still rough -- pending open questions
   
    if(USE_VZ_BINS){
        Vz = picoDst->event()->primaryVertex().z();
        VzBin = getVzBin(Vz);                  //get Vz bin
    }
   
    else VzBin = 0;
    
    if(centralityBin == -1 || VzBin == -1) { return kStOk; }
    
    pVtxX->Fill(picoDst->event()->primaryVertex().x());
    pVtxY->Fill(picoDst->event()->primaryVertex().y());
    pVtxZ->Fill(picoDst->event()->primaryVertex().z()); 
   
    vZandCentBinPerEvent->Fill(centralityBin+1, VzBin+1);   
   
    if(DEBUG){ 
        
                cout << endl << endl;
                cout << "*********************EVENT START******************" << endl;
                cout << "We are on event # " << eventNumber << endl;
                cout << "This event has " << realTracks << " tracks." << endl;
                if(USE_VZ_BINS){ cout << "Vz: " << Vz << "   VzBin: " << VzBin << "    Centrality Bin: " << centralityBin << endl; }
                else cout << "Centrality Bin: " << centralityBin << endl;
                cout << endl;
                cout << endl;
                eventNumber++;
    }         
           
        
/****************************************************************************************************************/  
/********************QA Loop for storing all of the basic information about the events and tracks END**********/
/****************************************************************************************************************/         
 

//***************************************************************************************
//BEGIN D0 BLOCK OF CODE HERE -- MAKE SURE THE NON_IDENTIFIED FLAG IS SET TO FALSE!!!!
//***************************************************************************************    

//Still need to add checks for TOF information-----   

/****************************************************************************************************************/  
/***************************************LOOP TO STORE EVENTS IN MIXER BEGIN***************************************/
/****************************************************************************************************************/   
   
    for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){//begin loop to check if event can be used (check if a D0 exists)
   
        StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
        StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
        StPicoTrack const* pion = picoDst->track(kp->pionIdx());
                    
                    //  ptmin ptmax   decayLenMin&Max   daughterDCA kaon/pion pt kaon/pion DCA  DCA to PV
        if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh && isGoodPair(kp) &&                              //Check for candidate in peak range (LS & US)
            cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
            d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) {

                     storeEventToMix = false;
                     storeKaonPionEvent = true;                     
                     if(DEBUG){
                         cout << "D0 Candidate information" << endl << endl;
                         cout << "Mass: "<< kp->m() << "       US (-1) or LS (1): " << kaon->charge()*pion->charge() << endl;
                         cout << endl;
                     }
                     
                     numD0s++;
        }     

       else if(kp->m() > USSideBandLeftLow && kp->m() < USSideBandLeftHigh && isGoodPair(kp) && kaon->charge()*pion->charge() < 0 &&             //Check for candidate in left side band (US)
            cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
            d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) {

                     if(DEBUG){
                         cout << "Left side band Candidate information" << endl << endl;
                         cout << "Mass: "<< kp->m() << endl;
                         cout << endl;
                     } 
            
                     storeEventToMix = false;
                     storeKaonPionEvent = true;                     
            }               
            
      else if(kp->m() > USSideBandRightLow && kp->m() < USSideBandRightHigh && isGoodPair(kp) && kaon->charge()*pion->charge() < 0 &&          //Check for candidate in right side band (US)
            cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
            d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) {

                     if(DEBUG){
                         cout << "Right side band Candidate information" << endl << endl;
                         cout << "Mass: "<< kp->m() << endl;
                         cout << endl;
                     } 
            
                     storeEventToMix = false;
                     storeKaonPionEvent = true;                     
            }                     
        
    }//end loop to check if event can be used (check if a D0 exists)        
    
    if(DEBUG) {cout << "Number of D0 candidates in this event: " << numD0s << endl << endl; }
    
    if(EVENT_MIXING){  // begin event mixing on/off switch
        
        if(eventBufferPicoEvent[VzBin][centralityBin]->getBufferSize() == BUFFER_SIZE) { storeEventToMix = false; }
        
        if(storeEventToMix){ //Begin block to store non-d0 containing events 
    
            eventBufferPicoEvent[VzBin][centralityBin]->addEvent(picoDst);
               
            for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){//begin loop to add picoDST tracks to buffer
                         
                        trk = picoDst->track(i);                                                             
                        trackMom = trk->gMom(pVtx, bField); 
                        if(!mHFCuts->isGoodTrack(trk)) { continue; }
                        
                        trackDCA = ((trk->helix().origin())-pVtx).mag();
        
                        if(trk->chi2() > trackChi2max) { continue; }  //check chi2 cut
                        if(!checkDCAtoPV(trackDCA))    { continue; }   // track quality cut for DCA to PV
                        
                        pt   = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y()));
                        if(pt < hadronPtMin || pt > hadronPtMax){continue;}   
                        
                        if(trackMom.mag() > .20 && trackMom.mag() < .45 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
                        if(trackMom.mag() > .70 && trackMom.mag() < .80 && trk->nSigmaElectron() > -1.5 && trk->nSigmaElectron() < 1.5) { continue; }
                        
                        eta  = trackMom.pseudoRapidity();
                        if(eta > 1 || eta < -1) { continue; }
                        
                        if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaKaon()) < 2.0)){ PIDflag = 1; }
                        else if(mHFCuts->isGoodTrack(trk) && (fabs(trk->nSigmaPion()) < 3.0)) { PIDflag = 2; }
                        else PIDflag = 0;
            
                        eventBufferPicoEvent[VzBin][centralityBin]->addTrackToEvent(eventBufferPicoEvent[VzBin][centralityBin]->getBufferIndex()-1, trackMom, trk->charge(), PIDflag);   //0 -- any hadron, 1 -- kaon, 2 -- pion
                  
                    }//end loop to add picoDST tracks to buffer
                    
                    if(DEBUG){
                         
                        int bufIndex = 0;
                        bufIndex = eventBufferPicoEvent[VzBin][centralityBin]->getBufferIndex();                        
                        cout << "Event stored in buffer." << "Num of tracks = " << eventBufferPicoEvent[VzBin][centralityBin]->getEvent(bufIndex-1)->getNoTracks() << endl;
                   }
                   
        } // end conditional for the event being stored in the mixer
        
        if(storeKaonPionEvent){//begin flag to store D0  
        
            eventBufferD0Candidate[VzBin][centralityBin]->addEvent(picoDst);
            
            for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){//begin loop to store D0s to buffer
   
                StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
                StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
                StPicoTrack const* pion = picoDst->track(kp->pionIdx());
                
                if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh && isGoodPair(kp) &&                                 //begin conditional to store D0 candidate to buffer
                    cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
                    d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) {          
                    
                    kaonPionMom.set(kp->lorentzVector().px(),kp->lorentzVector().py(),kp->lorentzVector().pz());
                   
                    eventBufferD0Candidate[VzBin][centralityBin]->addKaonPionToEvent(eventBufferD0Candidate[VzBin][centralityBin]->getBufferIndex()-1, kaonPionMom, kp->m(), 
                                                                                     kp->kaonIdx(), kp->pionIdx(), kaon->charge()*pion->charge());
                    
               }//end conditional to store D0 candidate to buffer
               
               if(kp->m() > USSideBandLeftLow && kp->m() < USSideBandLeftHigh && isGoodPair(kp) && kaon->charge()*pion->charge() < 0 &&  //begin conditional to store US left sideband
                    cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
                    d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) {          
                    
                    kaonPionMom.set(kp->lorentzVector().px(),kp->lorentzVector().py(),kp->lorentzVector().pz());
                    
                    eventBufferD0Candidate[VzBin][centralityBin]->addKaonPionToEvent(eventBufferD0Candidate[VzBin][centralityBin]->getBufferIndex()-1, kaonPionMom, kp->m(), 
                                                                                     kp->kaonIdx(), kp->pionIdx(), kaon->charge()*pion->charge());
                    
               }//end conditional to store US left sideband
               
               if(kp->m() > USSideBandRightLow && kp->m() < USSideBandRightHigh && isGoodPair(kp) && kaon->charge()*pion->charge() < 0 &&  //begin conditional to store US right sideband
                    cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
                    d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) {          
                    
                    kaonPionMom.set(kp->lorentzVector().px(),kp->lorentzVector().py(),kp->lorentzVector().pz());
                    
                    eventBufferD0Candidate[VzBin][centralityBin]->addKaonPionToEvent(eventBufferD0Candidate[VzBin][centralityBin]->getBufferIndex()-1, kaonPionMom, kp->m(), 
                                                                                     kp->kaonIdx(), kp->pionIdx(), kaon->charge()*pion->charge());
                    
               }//end conditional to store US right sideband
           }//end loop to store D0s to buffer
        }//end flag to store D0   
     }//end event mixing on/off switch    
/****************************************************************************************************************/  
/***************************************LOOP TO STORE EVENTS IN MIXER END***************************************/
/****************************************************************************************************************/   
   
/****************************************************************************************************************/   
/*****************************************BEGIN MAIN SIBLING LOOP************************************************/
/****************************************************************************************************************/
    if(DEBUG) { cout << "   Begin sibling pair code block" << endl << endl; }
   
    for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx){ // begin main sibling loop
     
        StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    
        StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
        StPicoTrack const* pion = picoDst->track(kp->pionIdx());
      
        if(!isGoodPair(kp)) continue;   
        if(!cutCheck(kp, d0PtLow,  d0PtHigh,  d0DecayLengthMin,  d0DecayLengthMax,  daughterDCA,
                     d0DaughterPionPtMin, d0DaughterKaonPtMin, kaonDCA,  pionDCA,  d0DCAtoPV)) { continue; }    
      
        ptBin = getPtBin(kp->pt());
        
        if(kaon->charge()*pion->charge() < 0){// begin Unlike-sign conditional 
	      
            invMass->Fill(kp->m());     
            if(ptBin > -1) { D0InvMassPtBin[ptBin]->Fill(kp->m()); }           
            
	    }//end Unlike-sign conditional 
      
	    if(kaon->charge()*pion->charge() > 0){//begin Like-sign conditional 
          
            likeSignBG->Fill(kp->m());   
            if(ptBin > -1) { LSInvMassPtBin[ptBin]->Fill(kp->m()); }
          
        }//end Like-sign conditional 
 
     
    /****************************************************************************************************************/
    /****************************SIBLING EVENT PAIRS FORMATIION BEGINS HERE*****************************************/  
    /***************************************************************************************************************/     
        
        bandBin = -1;  //set default case value -- no pair formation
        
        if(kaon->charge()*pion->charge() < 0){// begin US conditional
                
            if(kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh){ // begin loop over both unlike-sign AND LIKE SIGN D0 candidates 

                if(DEBUG) { cout << "   D0 Candidate Mass: " << kp->m() << "    PtBin : " << ptBin << endl; }
            
                d0Counter = d0Counter + 1;
                D0ptDist->Fill(kp->pt()); 
                D0EtaDist->Fill(kp->eta());
                D0PhiDist->Fill(kp->phi());
                    
                bandBin = 1;
            }
                
            else if(kp->m() > USSideBandLeftLow &&  kp->m() < USSideBandLeftHigh){ //Side Band Left
                
                if(DEBUG) { cout << "   SideBandLeft: " << kp->m() << "    PtBin : " << ptBin << endl;}    
                bandBin = 0; 
            }
            
            else if(kp->m() > USSideBandRightLow &&  kp->m() < USSideBandRightHigh){ //Side Band Right
            
                if(DEBUG) { cout << "   SideBandRight: " << kp->m() << "    PtBin : " << ptBin << endl;}    
                bandBin = 2; 
            }
                
        }  //end US conditional  
            
        else if(kaon->charge()*pion->charge() > 0 && kp->m() > D0InvMassLow && kp->m() < D0InvMassHigh) {// like sign conditional 

            if(DEBUG) { cout << "   LS Mass: " << kp->m() << "    PtBin : " << ptBin << endl; }
            bandBin = 3; 
        }          
            
        
            
            if(bandBin > -1 && D0_HADRON_CORR){// begin sibling pair formation with whatever "band" we are in, if we have a candidate in the band -- NEED TO MAKE A TRACK LIST TO SPEED THIS UP!!!
                
                realTracks = 0;
                
                for(unsigned int i = 0; i < picoDst->numberOfTracks(); ++i){ // begin picoDST loop for d0-hadron correlations
           
                    if(i == kp->kaonIdx() || i == kp->pionIdx()) { continue; }    // Need to check this -- should avoid doing correlations with a D0 candidate daughter
                    
                    
                    trackMom = mAssociatedHadronList[i];
        
                    realTracks++;
                    
                    delPhi = kp->phi()-phi;
                    
                    if(delPhi < -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                    else if(delPhi > TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                    
                    delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                    if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                    delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                    if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                    delEta = TMath::Abs(kp->eta()-eta);
                   
                    sibCorrBin[bandBin][VzBin][centralityBin]->Fill(delEta, delPhi);
                    sibCorrBin[bandBin][VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                    sibCorrBin[bandBin][VzBin][centralityBin]->Fill(-delEta, delPhi);
                    sibCorrBin[bandBin][VzBin][centralityBin]->Fill(delEta, delPhiCp);
                    if(ptBin > -1) { 
                    
                        sibCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(delEta, delPhi);
                        sibCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(-delEta, delPhiCp); 
                        sibCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(-delEta, delPhi);
                        sibCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(delEta, delPhiCp); 
                    }
                
                }// end picoDST loop for d0-hadron correlations
                
                if(DEBUG) { cout << endl << "   Number of pairs with this candidate: " << realTracks << endl; }
        
                
           }// end sibling pair formation with whatever "band" we are in, if we have a candidate in the band
    }// End main sibling loop

    
    
    /****************************************************************************************************************/
    /****************************SIBLING EVENT PAIRS FORMATIION ENDS HERE********************************************/  
    /****************************************************************************************************************/
    
/****************************************************************************************************************/   
/*****************************************END MAIN SIBLING LOOP**************************************************/
/****************************************************************************************************************/    
    
    
//------------------------------------------------------------------------------------------------------------------------//   
    
/***************************************************************************************************************/
/*************************************EVENT MIXING BEINGS HERE**************************************************/
/***************************************************************************************************************/
    
    if(EVENT_MIXING){// begin switch to turn on event mixing
        if((eventBufferPicoEvent[VzBin][centralityBin]->getBufferMaxSize() == eventBufferPicoEvent[VzBin][centralityBin]->getBufferSize()) //begin buffer-full mixing conditional
            && eventBufferD0Candidate[VzBin][centralityBin]->getBufferSize() >= 1){ 
    
            if(DEBUG_MIX_BUFFER){cout << "Mixing events in Vz/Centrality Bin " << VzBin<< "/" << centralityBin << endl << endl;}
            
            StMixerEvent* eventWithD0;
            StMixerEvent* hadronEvent;
            StMixerTrack kaonPionTrack1;
            StMixerTrack hadronTrack2;
            int nD0Tracks = 0;
            int nHadronTracks = 0;
            // D0Mom = 0;
            double D0Pt = 0;
        
            
            
            if(DEBUG_MIX_BUFFER) { cout << "Code gets here" << endl; }
        
            for(int k = 0; k <  eventBufferD0Candidate[VzBin][centralityBin]->getBufferSize(); k++){ // begin d0 event buffer loop
            
                eventWithD0 = eventBufferD0Candidate[VzBin][centralityBin]->getEvent(k); //This event in the buffer will be the trigger D0 event
                nD0Tracks = eventWithD0->getKaonPionListSize();
        
                for(int i = 0; i < eventBufferPicoEvent[VzBin][centralityBin]->getBufferSize(); i++){ //begin associated event buffer loop
            
                    hadronEvent = eventBufferPicoEvent[VzBin][centralityBin]->getEvent(i); 
                    nHadronTracks = hadronEvent->getNoTracks();
                    if(DEBUG_MIX_BUFFER) { cout << "Number of hadron track: " << nHadronTracks << endl; }   
                    if(DEBUG_MIX_BUFFER) { cout << "buffer position for hadron event: " << i << endl; }
            
                    for (int idx = 0; idx < nD0Tracks; idx++){//begin loop over event D0 kaon-pion list
                
                        kaonPionTrack1 = eventWithD0->getKaonPionAt(idx);
                        if(DEBUG_MIX_BUFFER) { cout << "index position for D0 track: " << idx << endl; }
                        kaonPionMom = kaonPionTrack1.gMom();
                        D0Pt = TMath::Sqrt((kaonPionMom.x()*kaonPionMom.x())+(kaonPionMom.y()*kaonPionMom.y()));
                    
                        ptBin = getPtBin(D0Pt);
                    
                        if(kaonPionTrack1.charge() < 0 && kaonPionTrack1.mass() > D0InvMassLow && kaonPionTrack1.mass() < D0InvMassHigh)                 { bandBin = 1; }    //US center band
                        else if(kaonPionTrack1.charge() < 0 && kaonPionTrack1.mass() > USSideBandLeftLow &&  kaonPionTrack1.mass() < USSideBandLeftHigh) { bandBin = 0; }    //US left band
                        else if(kaonPionTrack1.charge() < 0 && kaonPionTrack1.mass() > USSideBandRightLow && kaonPionTrack1.mass() < USSideBandRightHigh){ bandBin = 2; }    //US right band
                        else if(kaonPionTrack1.charge() > 0 && kaonPionTrack1.mass() > D0InvMassLow && kaonPionTrack1.mass() < D0InvMassHigh)            { bandBin = 3; }    //LS conditional
                    
                            if(DEBUG_MIX_BUFFER) { cout << "Mixing with mother in bin: " << bandBin << endl; }
                    
                            for(int j = 0; j < nHadronTracks; j++){//begin loop over hadron event tracks
                        
                                hadronTrack2 = hadronEvent->getTrack(j);  
                                             
                                trackMom = hadronTrack2.gMom();
               
                                pt = TMath::Sqrt((trackMom.x()*trackMom.x())+(trackMom.y()*trackMom.y())); 
                                phi = TMath::ATan2(trackMom.y(),trackMom.x());  
                                eta = trackMom.pseudoRapidity();
                 
                                delPhi = kaonPionTrack1.gMom().phi()-phi;
                                if(delPhi <= -TMath::Pi()) { delPhi = delPhi + 2*TMath::Pi(); }     //shift [-2Pi, 2Pi] -> [-Pi,Pi]
                                else if(delPhi >= TMath::Pi()) { delPhi = delPhi - 2*TMath::Pi(); }
                 
                                delPhi = TMath::Abs(delPhi);//Gets absolute value of angle 
                                if(delPhi >= 3*TMath::PiOver2()){ delPhi = delPhi - 2*TMath::Pi(); }//....and shifts it
                    
                                delPhiCp = -delPhi;//Gets negative copy of absolute value of angle.....
                                if(delPhiCp < -TMath::PiOver2()){ delPhiCp = delPhiCp + 2*TMath::Pi(); }//....and shifts it
                
                
                                delEta = TMath::Abs(kaonPionTrack1.gMom().pseudoRapidity()-eta);
                
                                mixCorrBin[bandBin][VzBin][centralityBin]->Fill(delEta, delPhi); //delEta & delPhi stored in US mixed correlation histogram
                                mixCorrBin[bandBin][VzBin][centralityBin]->Fill(-delEta, delPhiCp); //delEta & delPhi stored in US mixed correlation histogram
                                mixCorrBin[bandBin][VzBin][centralityBin]->Fill(-delEta, delPhi); //delEta & delPhi stored in US mixed correlation histogram
                                mixCorrBin[bandBin][VzBin][centralityBin]->Fill(delEta, delPhiCp); //delEta & delPhi stored in US mixed correlation histogram
                                if(ptBin > -1) { 
                                
                                    mixCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(delEta, delPhi);
                                    mixCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(-delEta, delPhiCp);
                                    mixCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(-delEta, delPhi);
                                    mixCorrBinPt[bandBin][ptBin][VzBin][centralityBin]->Fill(delEta, delPhiCp);
                                }
                            }// end loop over hadron event tracks
                        }// end loop over event D0 kaon-pion list
                    }// end associated event buffer loop
                }// end d0 event buffer loop
            
            eventBufferPicoEvent[VzBin][centralityBin]->removeFirstEvent(); //MAY NEED TO FIX THIS ALGORITHM TO AVOID A BIAS
            
            for(int numEvents = 0; numEvents < eventBufferD0Candidate[VzBin][centralityBin]->getBufferSize(); numEvents++){ //This will remove all of the events in the D0 buffer, since they are now used.
            
                eventBufferD0Candidate[VzBin][centralityBin]->removeFirstEvent(); 
            
            }
      
            if(DEBUG) { 
        
                cout << endl << endl;
                cout << "First Event in US buffer for Vz/centrality bin "<< VzBin << "/"<< centralityBin<< " cleared" << endl;
                cout << endl;
            }     
        
        }// end buffer-full US mixing conditional        
    }// end switch to turn on event mixer  
    
/***************************************************************************************************************/
/*************************************EVENT MIXING ENDS HERE****************************************************/
/***************************************************************************************************************/   
       
       
       
       d0CountPerEvent->Fill(d0Counter);
       //invMassMinusBG->Add(invMass, likeSignBG, 1, -1);
       //D0PeakMinusBG->Add(D0PeakPlusBG, D0LikeSignBG, 1, -1);
    
   //if(DEBUG) { cout << endl << "***************EVENT END****************************" << endl; }
    
   //-------------------End Current Event Analysis--------------------------------
   
   return kStOK;

}//end Make member function

//---------------------User Functions-------------------------------------

bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  //  To be replaced by mHFCuts->isGoodSecondaryVertexPair(kp))
  bool pairCuts = kp->m() > mHFCuts->cutSecondaryPairMassMin() && 
    kp->m() < mHFCuts->cutSecondaryPairMassMax() &&
    std::cos(kp->pointingAngle()) > mHFCuts->cutSecondaryPairCosThetaMin() &&
    kp->decayLength()  > mHFCuts->cutSecondaryPairDecayLengthMin() && 
    kp->decayLength()  < mHFCuts->cutSecondaryPairDecayLengthMax() &&
    kp->dcaDaughters() < mHFCuts->cutSecondaryPairDcaDaughtersMax();

  return (mHFCuts->isGoodTrack(kaon) && mHFCuts->isGoodTrack(pion) &&
	  mHFCuts->isTPCKaon(kaon) && mHFCuts->isTPCPion(pion) && 
	  pairCuts);

}


bool StPicoD0AnaMaker::cutCheck(StKaonPion const* const kp, double ptMin, double ptMax, double decayLengthMin, double decayLengthMax, 
                                                            double dcaDaughters, double kaonPtCut, double pionPtCut, 
                                                            double dcaKaon, double dcaPion, double dcaV0toPV) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  
      bool truthCuts = kp->pt() > ptMin && kp->pt() < ptMax &&
                       kp->decayLength() > decayLengthMin && 
                       kp->decayLength() < decayLengthMax &&
                       kp->dcaDaughters() < dcaDaughters  &&
                       kaon->gPt() > kaonPtCut && pion->gPt() > pionPtCut &&
                       kp->kaonDca() > dcaKaon && kp->pionDca() > dcaPion &&
                       kp->perpDcaToVtx() < dcaV0toPV;
  
      return truthCuts;

}


int StPicoD0AnaMaker::getCentralityBin(int nTracks){

    /*if(nTracks >= 2   && nTracks < 16)  { return 0;  }
    if(nTracks >= 16  && nTracks < 31)  { return 1;  }
    if(nTracks >= 31  && nTracks < 51)  { return 2;  }
    if(nTracks >= 51  && nTracks < 79)  { return 3;  }
    if(nTracks >= 79  && nTracks < 116) { return 4;  }
    if(nTracks >= 116 && nTracks < 164) { return 5;  }
    if(nTracks >= 164 && nTracks < 224) { return 6;  }
    if(nTracks >= 224 && nTracks < 295) { return 7;  }
    if(nTracks >= 295 && nTracks < 335) { return 8;  }
    if(nTracks >= 335 && nTracks < 376) { return 9;  }
    if(nTracks >= 376 && nTracks < 420) { return 10; }
    if(nTracks >= 420 && nTracks < 470) { return 11; }
    if(nTracks >= 470 && nTracks < 550) { return 12; }
    if(nTracks >= 550 && nTracks < 620) { return 13; }
    if(nTracks >= 620)                  { return 14; }*/
    
    if(nTracks >= 2   && nTracks < 14)  { return 0;  }
    if(nTracks >= 14  && nTracks < 32)  { return 1;  }
    if(nTracks >= 34  && nTracks < 67)  { return 2;  }
    if(nTracks >= 67  && nTracks < 115) { return 3;  }
    if(nTracks >= 115 && nTracks < 183) { return 4;  }
    if(nTracks >= 183 && nTracks < 243) { return 5;  }
    if(nTracks >= 243 && nTracks < 300) { return 6;  }
    if(nTracks >= 300 && nTracks < 370) { return 7;  }
    if(nTracks >= 370 && nTracks < 450) { return 8;  }
    if(nTracks >= 450 && nTracks < 520) { return 9;  }
    if(nTracks >= 520 && nTracks < 580) { return 10;  }
    if(nTracks >= 580 && nTracks < 650) { return 11;  }
    if(nTracks >= 650 && nTracks < 710) { return 12;  }
    if(nTracks >= 710 && nTracks < 770) { return 13;  }
    if(nTracks >= 770)                  { return 14; }

    else return -1;

}    

int StPicoD0AnaMaker::getVzBin(double Vz){

    if(Vz >= -6.0  && Vz < -4.8)  { return 0;  }   //bin 0: -6 to -4.8
    if(Vz >= -4.8  && Vz < -3.6)  { return 1;  }      //bin 1: -4.8 to -3.6
    if(Vz >= -3.6  && Vz < -2.4)  { return 2;  }          //bin 2: -3.6 to -2.4
    if(Vz >= -2.4  && Vz < -1.2)  { return 3;  }      //bin 3: -2.4 to -1.2
    if(Vz >= -1.2  && Vz < 0)     { return 4;  }      //bin 4: -1.2 to 0
    if(Vz >= 0     && Vz < 1.2)   { return 5;  }
    if(Vz >= 1.2   && Vz < 2.4)   { return 6;  }
    if(Vz >= 2.4   && Vz < 3.6)   { return 7;  }
    if(Vz >= 3.6   && Vz < 4.8)   { return 8;  }
    if(Vz >= 4.8   && Vz < 6.0)   { return 9;  }
    

    else return -1;

}    

int StPicoD0AnaMaker::getPtBin(double pt){

    if(pt >=  0.0  && pt <  1.0)    { return 0;  }   
    if(pt >=  1.0  && pt <  4.0)    { return 1;  }
    if(pt >=  4.0  && pt <  20.0)   { return 2;  }        
    //if(pt >=  3.0  && pt <  4.0)   { return 3;  }
    //if(pt >=  4.0  && pt <  5.0)   { return 4;  }
    //if(pt >=  5.0  && pt <  10.0)  { return 5;  }
    
    else return -1;

}    

bool StPicoD0AnaMaker::checkDCAtoPV(float trackDCA){

     return (trackDCA <= trackDCAtoPvtx);
     
}

