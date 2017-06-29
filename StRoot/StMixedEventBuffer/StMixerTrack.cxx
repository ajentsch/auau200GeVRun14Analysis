#include <limits>

#include "StMixerTrack.h"
#include "StPicoDstMaker/StPicoTrack.h"

StMixerTrack::StMixerTrack() : mMom(StThreeVectorF()), mPIDFlag(-1), mCharge(0), mMass(-1), mKaonIdx(-1), mPionIdx(-1), mKaonMom(StThreeVectorF()), mPionMom(StThreeVectorF())
{
}

StMixerTrack::StMixerTrack(StMixerTrack const * trk) : mMom(trk->mMom), mPIDFlag(trk->mPIDFlag), mCharge(trk->mCharge),
                                                       mMass(trk->mMass), mKaonIdx(trk->mKaonIdx), mPionIdx(trk->mPionIdx)
{
} 

StMixerTrack::StMixerTrack(StThreeVectorF trackMomentum, int charge, int PID){

    mMom = trackMomentum;
    mCharge = charge;
    mPIDFlag = PID;

}    

StMixerTrack::StMixerTrack(StThreeVectorF trackMomentum, double mass, int kaonIdx, int pionIdx, int charge, StThreeVectorF kaonMom, StThreeVectorF pionMom){

    mMom = trackMomentum;
    mMass = mass;
    mKaonIdx = kaonIdx;
    mPionIdx = pionIdx;
    mCharge = charge;  //charge = -1 for D0 and +1 for a LS KPi pair
    mKaonMom = kaonMom;
    mPionMom = pionMom;
}    