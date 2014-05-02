
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
 *  \author Z. Wu
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
  //Double_t fZVertexResolution;

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
  fItIsoTrackInputArray(0), fItCandidateInputArray(0),
  fItRhoInputArray(0)
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
  const char *rhoInputArrayName;

  fDeltaRMax = GetDouble("DeltaRMax", 0.5);

  fPTRatioMax = GetDouble("PTRatioMax", 0.1);

  fPTSumMax = GetDouble("PTSumMax", 5.0);

  fUsePTSum = GetBool("UsePTSum", false);


  fIsoTrackPTMin = GetDouble("IsoTrackPTMin", 10);
  fIsoTrackEtaMax = GetDouble("IsoTrackEtaMax", 2.4);

  fClassifier->fPTMin = GetDouble("PTMin", 0.5);
  //fClassifier->fZVertexResolution  = GetDouble("ZVertexResolution", 0.005)*1.0E3;

  // import input array(s)

  fIsoTrackInputArray = ImportArray(GetString("IsolationInputArray", "Delphes/partons"));
  fItIsoTrackInputArray = fIsoTrackInputArray->MakeIterator();

  fIsolationFilter = new ExRootFilter(fIsoTrackInputArray);

  fCandidateInputArray = ImportArray(GetString("CandidateInputArray", "Calorimeter/electrons"));
  fItCandidateInputArray = fCandidateInputArray->MakeIterator();

  fCandidateFilter = new ExRootFilter(fCandidateInputArray);

  rhoInputArrayName = GetString("RhoInputArray", "");
  if(rhoInputArrayName[0] != '\0')
  {
    fRhoInputArray = ImportArray(rhoInputArrayName);
    fItRhoInputArray = fRhoInputArray->MakeIterator();
  }
  else
  {
    fRhoInputArray = 0;
  }

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "electrons"));
}

//------------------------------------------------------------------------------

void IsoTrack::Finish()
{
  if(fItRhoInputArray) delete fItRhoInputArray;
  if(fIsolationFilter) delete fIsolationFilter;
  if(fCandidateFilter) delete fCandidateFilter;
  if(fItCandidateInputArray) delete fItCandidateInputArray;
  if(fItIsoTrackInputArray) delete fItIsoTrackInputArray;
}

//------------------------------------------------------------------------------

void IsoTrack::Process()
{
  Candidate *candidate, *IsoTrack, *object;
  TObjArray *IsoTrackArray;
  TObjArray *CanTrackArray;
  Double_t sum, ratio;
  Int_t counter;
  Double_t eta = 0.0;
  Double_t rho = 0.0;

  if(fRhoInputArray && fRhoInputArray->GetEntriesFast() > 0)
  {
    candidate = static_cast<Candidate*>(fRhoInputArray->At(0));
    rho = candidate->Momentum.Pt();
  }

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

    // find rho
    rho = 0.0;
    if(fRhoInputArray)
    {
      fItRhoInputArray->Reset();
      while((object = static_cast<Candidate*>(fItRhoInputArray->Next())))
      {
        if(eta >= object->Edges[0] && eta < object->Edges[1])
        {
          rho = object->Momentum.Pt();
        }
      }
    }

    // correct sum for pile-up contamination
    sum = sum - rho*fDeltaRMax*fDeltaRMax*TMath::Pi();

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
