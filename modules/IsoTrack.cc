
/** \class IsoTrack
 *
 *  Sums transverse momenta of Track objects within a DeltaR cone around a
 *  candidate and calculates fraction of this sum to the candidate's
 *  transverse momentum. outputs candidates that have the transverse momenta
 *  fraction within (PTRatioMin, PTRatioMax].
 *
 *  $Date: 2013-11-04 13:14:33 +0100 (Mon, 04 Nov 2013) $
 *  $Revision: 1317 $
 *
 *
 *  \author  Z. Wu - Baylor
 *
 */

#include "modules/IsoTrack.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

class IsoTrackClassifier : public ExRootClassifier
{
public:

  IsoTrackClassifier() {}

  Int_t GetCategory(TObject *object);

  Double_t fPTMin;

};

//------------------------------------------------------------------------------

Int_t IsoTrackClassifier::GetCategory(TObject *object)
{
  Candidate *track = static_cast<Candidate*>(object);
  const TLorentzVector &momentum = track->Momentum;

  if(momentum.Pt() < fPTMin) return -1;

  return 0;
}

//------------------------------------------------------------------------------

IsoTrack::IsoTrack() :
  fClassifier(0), 
  fIsolationFilter(0), fCandidateFilter(0),
  fItIsoTrackInputArray(0), fItCandidateInputArray(0)
{
  fClassifier = new IsoTrackClassifier;
}

//------------------------------------------------------------------------------

IsoTrack::~IsoTrack()
{
}

//------------------------------------------------------------------------------

void IsoTrack::Init()
{

  fDeltaRMax = GetDouble("DeltaRMax", 0.5);

  fPTRatioMax = GetDouble("PTRatioMax", 0.1);

  fPTSumMax = GetDouble("PTSumMax", 5.0);

  fUsePTSum = GetBool("UsePTSum", false);


  fIsoTrackPTMin = GetDouble("IsoTrackPTMin", 10);
  fIsoTrackEtaMax = GetDouble("IsoTrackEtaMax", 2.4);

  fClassifier->fPTMin = GetDouble("PTMin", 0.5);

  // import input array(s)

  fIsoTrackInputArray = ImportArray(GetString("IsolationInputArray", "TrackMerger/tracks"));
  fItIsoTrackInputArray = fIsoTrackInputArray->MakeIterator();

  fIsolationFilter = new ExRootFilter(fIsoTrackInputArray);

  fCandidateInputArray = ImportArray(GetString("CandidateInputArray", "Calorimeter/electrons"));
  fItCandidateInputArray = fCandidateInputArray->MakeIterator();

  fCandidateFilter = new ExRootFilter(fCandidateInputArray);

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "electrons"));
}

//------------------------------------------------------------------------------

void IsoTrack::Finish()
{
  if(fIsolationFilter) delete fIsolationFilter;
  if(fCandidateFilter) delete fCandidateFilter;
  if(fItCandidateInputArray) delete fItCandidateInputArray;
  if(fItIsoTrackInputArray) delete fItIsoTrackInputArray;
}

//------------------------------------------------------------------------------

void IsoTrack::Process()
{
  Candidate *candidate, *IsoTrack;
  TObjArray *IsoTrackArray;
  TObjArray *CanTrackArray;
  Double_t sum, ratio;
  Int_t counter;
  Double_t eta = 0.0;

  // select IsoTrack objects
  fIsolationFilter->Reset();
  IsoTrackArray = fIsolationFilter->GetSubArray(fClassifier, 0);

  fCandidateFilter->Reset();
  CanTrackArray = fCandidateFilter->GetSubArray(fClassifier, 0);

  if(IsoTrackArray == 0) return;

  TIter itIsoTrackArray(IsoTrackArray);
  TIter itCanTrackArray(CanTrackArray);

  // loop over all input jets
  //fItCandidateInputArray->Reset();
  itCanTrackArray.Reset();
  while((candidate = static_cast<Candidate*>(itCanTrackArray.Next())))
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = TMath::Abs(candidateMomentum.Eta());

    // loop over all input tracks
    sum = 0.0;
    counter = 0;
    itIsoTrackArray.Reset();
    while((IsoTrack = static_cast<Candidate*>(itIsoTrackArray.Next())))
    {
      const TLorentzVector &IsoTrackMomentum = IsoTrack->Momentum;

      if(candidateMomentum.DeltaR(IsoTrackMomentum) <= fDeltaRMax &&
         !candidate->Overlaps(IsoTrack))
      {
        sum += IsoTrackMomentum.Pt();
        ++counter;
      }
    }

    ratio = sum/candidateMomentum.Pt();

    candidate->IsolationVar = ratio;

    if((fUsePTSum && sum > fPTSumMax) || ratio > fPTRatioMax) continue;

    if (candidate->Momentum.Pt() > fIsoTrackPTMin && 
        std::fabs(candidate->Momentum.Eta() < fIsoTrackEtaMax))
    {
      fOutputArray->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------
