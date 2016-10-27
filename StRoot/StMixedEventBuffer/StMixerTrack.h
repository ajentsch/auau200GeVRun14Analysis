#ifndef StMixerTrack_hh
#define StMixerTrack_hh


/* **************************************************
 *
 * Track class used for mixed event buffer
 *
 *
 * ************************************************/

#include "StarClassLibrary/StThreeVectorF.hh"

class StPicoTrack;

class StMixerTrack{
 public:
  StMixerTrack();
  StMixerTrack(StMixerTrack const * trk);
  ~StMixerTrack(){/*cout << "track destructor called" << endl*/;};
  StMixerTrack(StThreeVectorF trackMomentum, int charge, int PID);
  StMixerTrack(StThreeVectorF trackMomentum, double mass, int kaonIdx, int pionIdx, int charge);
  
  StThreeVectorF const gMom() const;
  int particleType();
  int charge();
  double mass();         
  int kaonIdx();      
  int pionIdx();      
  
  
 private:
  
  StThreeVectorF mMom;
  int mPIDFlag;
  int mCharge;
  double mMass;
  int mKaonIdx;
  int mPionIdx;
  
};

inline StThreeVectorF const StMixerTrack::gMom() const { return(mMom); }
inline int StMixerTrack::particleType() { return mPIDFlag; }
inline int StMixerTrack::charge()       { return mCharge;  }
inline double StMixerTrack::mass()         { return mMass;    }
inline int StMixerTrack::kaonIdx()      { return mKaonIdx;    }
inline int StMixerTrack::pionIdx()      { return mPionIdx;    }

#endif
