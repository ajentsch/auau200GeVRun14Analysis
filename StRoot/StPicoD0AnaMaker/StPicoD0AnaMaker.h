#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

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


class StPicoD0AnaMaker : public StMaker
{
  public:
    StPicoD0AnaMaker(char const * name, char const * inputFilesList, 
                     char const * outName, StPicoDstMaker* picoDstMaker);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;

    void setHFCuts(StHFCuts* cuts);   

      

  private:
    StPicoD0AnaMaker() {}
    void readNextEvent();

    bool isGoodPair(StKaonPion const*) const;
    bool cutCheck(StKaonPion const*, double, double, double, double, double, double, double, double, double, double) const;
    int  getCentralityBin(int nTracks);
    int  getVzBin(double Vz);
    int  getPtBin(double pt);
    bool checkDCAtoPV(float trackDCA);

    int NUM_PHI_BINS;
    int NUM_ETA_BINS;
    int NUM_VZ_BINS;
    int NUM_CENT_BINS;
    int NUM_PT_BINS;
    
    StPicoDstMaker* mPicoDstMaker;
    StPicoD0Event* mPicoD0Event;
    
    StMixedEventBuffer* eventBufferPicoEvent[10][16];      //To use for the actual PicoEvent
    StMixedEventBuffer* eventBufferD0Candidate[10][16];      //To use for D0 candidate
    //StKaonPion*         d0CandidateBuffer[10][11];  //stores the candidate D0s use for eventual mixing
    
    TString mOutFileName;
    TString mInputFileList;
    TFile* mOutputFile;
    TChain* mChain;
    int mEventCounter;

    StHFCuts* mHFCuts;

    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate
    
    /*double ptRange[6] = {0.0, 1.0, 2.0, 3.0, 5.0, 10.0};
    double decayLengthCuts[5] = {.0145, .0181, .0212, .0247, .0259};
    double daughterDCACuts[5] = {.0084, .0066, .0057, .0050, .0060};
    double dcaKaonPV[5] = {.0103, .0091, .0095, .0079, .0058};
    double dcaPionPV[5] = {.0110, .0111, .0086, .0081, .0062};
    double dcaV0toPV[5] = {.0061, .0049, .0038, .0038, .0040};*/
    
    int eventNumber;
    int trackCount;
    double Vz;
    int centralityBin;
    int VzBin;
    int ptBin;
    int bandBin;
    
    int nVzBins;
    int nCentBins;
    int BUFFER_SIZE;
    
             //  ptmin ptmax   decayLenMin&Max   daughterDCA kaon/pion pt kaon/pion DCA  DCA to PV
        //if(!isGoodPair(kp)) continue;             
        //if(!cutCheck(kp, 0.15,  20.0,  .0200,  999999.0,  .0055,  1.2,  1.2,  .008,  .008,  .0065)) { continue; }
    
    double trackDCA;
    double kaonPtCut;
    double pionPtCut;
    double hadronPtMin;
    double hadronPtMax;
    float trackChi2max;
    float trackDCAtoPvtx;
    double d0PtLow;
    double d0PtHigh;
    double d0DecayLengthMin;
    double d0DecayLengthMax;
    double daughterDCA;
    double d0DaughterPionPtMin;
    double d0DaughterKaonPtMin;
    double kaonDCA;
    double pionDCA;
    double d0DCAtoPV;
    double D0InvMassLow;
    double D0InvMassHigh;
    double USSideBandLeftLow;
    double USSideBandLeftHigh;
    double USSideBandRightLow;
    double USSideBandRightHigh;
    
    bool DEBUG;                //important flags for debugging and for switching binning on and off
    bool DEBUG_MIX_BUFFER;
    bool USE_CENT_BINS;
    bool USE_VZ_BINS;
    bool D0_HADRON_CORR;
    bool EVENT_MIXING;
    bool USE_PT_BINS;
                                           // vpdmb-5-p-nobsmd-hlt // vpdmb-5-p-nobsmd-hlt // vpdmb-5-p-nobsmd // vpdmb-5-p-nobsmd // vpdmb-5-p-nobsmd 
    unsigned int trigger[5];    
    
    TH2D* sibCorrBin[4][10][16];
    TH2D* mixCorrBin[4][10][16];
    TH2D* sibCorrBinPt[4][6][10][16];
    TH2D* mixCorrBinPt[4][6][10][16];
    TH1D* etaDistVz[10];
    TH1D* phiDistVz[10];
    TH2D* etaPhiDistVz[10];
    TH1D* D0InvMassPtBin[6];
    TH1D* LSInvMassPtBin[6];
    TH1D* ptDist;
    TH1D* invMass;
    TH1D* likeSignBG;
    TH1D* D0EtaDist;
    TH1D* D0PhiDist;
    TH1I* eventCounter;
    TH1D* D0ptDist;
    TH1D* kaonPtDist;
    TH1D* pionPtDist;
    TH1D* kaonEtaDist;
    TH1D* pionEtaDist;
    TH1D* kaonPhiDist;
    TH1D* pionPhiDist;
    TH1I* trackCounter;    
    TH1D* hadronPtDist;
    TH1D* hadronPhiDist;
    TH1D* hadronEtaDist;
    TH2D* dEdxVsPt;
    TH2D* invBetaVsPt;
    TH1I* usedTracks;
    TH1I* d0CountPerEvent;
    TH2I* vZandCentBinPerEvent;
    TH1D* histOfCuts;
    TH1D* hadronChi2;
    TH1D* pVtxX;
    TH1D* pVtxY;
    TH1D* pVtxZ;
    TH1D* DCAtoPrimaryVertex;
    TH1D* DCAtoPrimaryVertexCut;
    TH2D* phiD0vsPhiH[10][16];     
    TH2D* etaD0vsEtaH[10][16];
    TH2D* phiD0vsEtaD0[10][16];
    TH2D* phiHvsEtaH[10][16];

    ClassDef(StPicoD0AnaMaker, 0)
};

inline int StPicoD0AnaMaker::getEntries() const 
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoD0AnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}

inline void StPicoD0AnaMaker::setHFCuts(StHFCuts* cuts)   
{ 
  mHFCuts = cuts; 
}

#endif
