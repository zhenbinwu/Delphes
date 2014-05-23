#ifndef IsoTrack_h
#define IsoTrack_h

/** \class IsoTrack
 *
 *  Sums transverse momenta of IsoTrack objects (tracks, calorimeter towers, etc)
 *  within a DeltaR cone around a candidate and calculates fraction of this sum
 *  to the candidate's transverse momentum. outputs candidates that have
 *  the transverse momenta fraction within (PTRatioMin, PTRatioMax].
 *
 *  $Date: 2013-11-04 13:14:33 +0100 (Mon, 04 Nov 2013) $
 *  $Revision: 1317 $
 *
 *
 *  \author Z. Wu - Baylor
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;

class ExRootFilter;
class IsoTrackClassifier;

class IsoTrack: public DelphesModule
{
public:

  IsoTrack();
  ~IsoTrack();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fDeltaRMax;

  Double_t fPTRatioMax;

  Double_t fPTSumMax;
  
  Double_t fIsoTrackPTMin;
  Double_t fIsoTrackEtaMax;

  Bool_t fUsePTSum;

  IsoTrackClassifier *fClassifier; //!

  ExRootFilter *fIsolationFilter;
  ExRootFilter *fCandidateFilter;

  TIterator *fItIsoTrackInputArray; //!

  TIterator *fItCandidateInputArray; //!

  TIterator *fItRhoInputArray; //!

  const TObjArray *fIsoTrackInputArray; //!

  const TObjArray *fCandidateInputArray; //!

  const TObjArray *fRhoInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(IsoTrack, 1)
};

#endif
