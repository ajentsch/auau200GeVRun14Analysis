#ifndef StPicoHHCorrMaker_h
#define StPicoHHCorrMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis. 
 *
 *  Please write your analysis in the ::Make() function.
 *
 *  Authors:  **Mustafa Mustafa (mmustafa@lbl.gov)
 *
 *  **Code Maintainer
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "TH2.h"
#include "StMixedEventBuffer/StMixedEventBuffer.h"

class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StHFCuts;
class StMixedEventBuffer;


class StPicoHHCorrMaker : public StMaker
{
  public:
    StPicoHHCorrMaker(char const * name, char const * inputFilesList, 
                     char const * outName, StPicoDstMaker* picoDstMaker);
    virtual ~StPicoHHCorrMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);   

      

  private:
    StPicoHHCorrMaker() {}
    void readNextEvent();

    bool isGoodPair(StKaonPion const*) const;
    bool cutCheck(StKaonPion const*, double, double, double, double, double, double, double, double, double, double) const;
    int  getCentralityBin(int nTracks);
    int  getVzBin(double Vz);

    int NUM_PHI_BINS;
    int NUM_ETA_BINS;
    int NUM_VZ_BINS;
    int NUM_CENT_BINS;
    
    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
  
    StMixedEventBuffer* eventBuffer[10][11];    //To use for unidentified hadron
    
    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    StHFCuts* mHFCuts;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    
    int eventNumber;
    int trackCount;
    double Vz;
    int centralityBin;
    int VzBin;
    
    int nVzBins;
    int nCentBins;
    int BUFFER_SIZE;
    
             //  ptmin ptmax   decayLenMin&Max   daughterDCA kaon/pion pt kaon/pion DCA  DCA to PV
        //if(!isGoodPair(kp)) continue;             
        //if(!cutCheck(kp, 0.15,  20.0,  .0200,  999999.0,  .0055,  1.2,  1.2,  .008,  .008,  .0065)) { continue; }
    
    double kaonPtCut;
    double pionPtCut;
    double hadronPtCutLow;
    double hadronPtCutHigh;
    
    //double d0PtLow;
    //double d0PtHigh;
    //double d0DecayLengthMin;
    //double d0DecayLengthMax;
    //double daughterDCA;
    //double d0DaughterPionPtMin;
    //double d0DaughterKaonPtMin;
    double kaonDCA;
    double pionDCA;
    //double d0DCAtoPV;
    
    //double D0InvMassLow;
    //double D0InvMassHigh;
    
    bool DEBUG;                //important flags for debugging and for switching binning on and off
    bool DEBUG_MIX_BUFFER;
    bool USE_CENT_BINS;
    bool USE_VZ_BINS;
    //bool D0_HADRON_CORR;
    bool NON_IDENTIFIED_CORR;
                                           // vpdmb-5-p-nobsmd-hlt // vpdmb-5-p-nobsmd-hlt // vpdmb-5-p-nobsmd // vpdmb-5-p-nobsmd // vpdmb-5-p-nobsmd 
    unsigned int trigger[5];    
    
    
    TH2D* sibCorrPlusPlus[10][11];  //for non-identified particles
    TH2D* sibCorrPlusMinus[10][11];
    TH2D* sibCorrMinusPlus[10][11];
    TH2D* sibCorrMinusMinus[10][11];
    
    TH2D* mixCorrPlusPlus[10][11];  //for non-identified particles
    TH2D* mixCorrPlusMinus[10][11];
    TH2D* mixCorrMinusPlus[10][11];
    TH2D* mixCorrMinusMinus[10][11];
    
    TH1D* etaDistVz[10];
    TH1D* phiDistVz[10];
    
    TH1D* eventCategoryCounter[10]; //organized by vz bin, each one binned by centrality
    
    
    TH1D* ptDist;
    TH1D* invMass;
    //TH1D* USsiblingAngCorrPhi;
    //TH1D* USsiblingAngCorrEta;
    TH1D* kaonDist;
    TH1D* pionDist;
    TH1D* likeSignBG;
    TH1D* invMassMinusBG;
    //TH1D* D0EtaDist;
    //TH1D* D0PhiDist;
   
    TH1D* eventCounter;
    
    //TH1D* kaonDCAfromD0;
    //TH1D* pionDCAfromD0;
    //TH1D* decayLengthQA;
    //TH1D* pointingAngleQA;
    //TH1D* daughterDCAQA;
    //TH1D* D0ptDist;
    TH1D* kaonPtDist;
    TH1D* pionPtDist;
    TH1D* kaonEtaDist;
    TH1D* pionEtaDist;
    TH1D* kaonPhiDist;
    TH1D* pionPhiDist;
    TH1D* kaonDCAprimary;
    TH1D* pionDCAprimary;
    /*TH1D* invMassBin1;	
    TH1D* invMassBin2;	
    TH1D* invMassBin3;	
    TH1D* invMassBin4;	
    TH1D* invMassBin5;
    TH1D* angCorrPhiBin1;
    TH1D* angCorrEtaBin1;
    TH2D* angCorr2DBin1;
    TH1D* angCorrPhiBin2;
    TH1D* angCorrEtaBin2;
    TH2D* angCorr2DBin2;
    TH1D* angCorrPhiBin3;
    TH1D* angCorrEtaBin3;
    TH2D* angCorr2DBin3;
    TH1D* angCorrPhiBin4;
    TH1D* angCorrEtaBin4;
    TH2D* angCorr2DBin4;
    TH1D* angCorrPhiBin5;
    TH1D* angCorrEtaBin5;
    TH2D* angCorr2DBin5;

    TH1D* likeSignBin1;	
    TH1D* likeSignBin2;	
    TH1D* likeSignBin3;	
    TH1D* likeSignBin4;	
    TH1D* likeSignBin5;	 
    TH1D* USminusLSBin1;
    TH1D* USminusLSBin2;
    TH1D* USminusLSBin3;
    TH1D* USminusLSBin4;
    TH1D* USminusLSBin5;*/
    
    //TH1D* D0PeakPlusBG;
    //TH1D* D0LikeSignBG;
    //TH1D* D0PeakMinusBG;
    TH1D* trackCounter;    

    TH1D* hadronPtDist;
    TH1D* hadronPhiDist;
    TH1D* hadronEtaDist;
    
    TH1D* mixedEventTest;
    TH1D* mixedHadronPtDist;
    TH1D* mixedHadronPhiDist;
    TH1D* mixedHadronEtaDist;
    
    TH1D* mixedEventKPInvMass;
    TH2D* dEdxVsPt;
    TH2D* invBetaVsPt;
    //TH2D* mixedEventAngCorr2DAll;
    //TH2D* fullCorr;
    //
    //TH2D* USsibling2DAngCorrAll;
    //TH2D* LSsibling2DAngCorrAll;
    
   // TH2D* USmixedEventAngCorr2DAll;
   // TH2D* LSmixedEventAngCorr2DAll;
    
    TH1I* d0CountPerEvent;

    ClassDef(StPicoHHCorrMaker, 0)
};

inline int StPicoHHCorrMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoHHCorrMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

inline void StPicoHHCorrMaker::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
