
/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables
 *
 *  \author S. Zenz
 *
 */

#include "modules/PileUpJetID.h"

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

PileUpJetID::PileUpJetID() :
  fItJetInputArray(0),fTrackInputArray(0),fNeutralInputArray(0)
{

}

//------------------------------------------------------------------------------

PileUpJetID::~PileUpJetID()
{

}

//------------------------------------------------------------------------------

void PileUpJetID::Init()
{
  fJetPTMin = GetDouble("JetPTMin", 20.0);
  fParameterR = GetDouble("ParameterR", 0.5);
  fUseConstituents = GetInt("UseConstituents", 0);

  fAverageEachTower = false; // for timing

  // import input array(s)

  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();


  //  cout << "BeforE SCZ additions in init" << endl;
  //  cout << GetString("TrackInputArray", "ParticlePropagator/tracks") << endl;
  //  cout << GetString("EFlowTrackInputArray", "ParticlePropagator/tracks") << endl;

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "ParticlePropagator/tracks"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();


  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

  //  cout << " end of INIT " << endl;

}

//------------------------------------------------------------------------------

void PileUpJetID::Finish()
{
  //  cout << "In finish" << endl;

  if(fItJetInputArray) delete fItJetInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;

}

//------------------------------------------------------------------------------

void PileUpJetID::Process()
{
  //  cout << "start of process" << endl;

  Candidate *candidate, *constituent;
  TLorentzVector momentum, area;

  //  cout << "BeforE SCZ additions in process" << endl;

  // SCZ
  Candidate *trk;

  // loop over all input candidates
  fItJetInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    momentum = candidate->Momentum;
    area = candidate->Area;

    float sumT0 = 0.;
    float sumT1 = 0.;
    float sumT10 = 0.;
    float sumT20 = 0.;
    float sumT30 = 0.;
    float sumT40 = 0.;
    float sumWeightsForT = 0.;
    candidate->nTimes = 0;

    float sumpt = 0.;
    float sumptch = 0.;
    float sumptchpv = 0.;
    float sumptchpu = 0.;
    float sumdrsqptsq = 0.;
    float sumptsq = 0.;
    int nc = 0;
    int nn = 0;
    float pt_ann[5];

    for (int i = 0 ; i < 5 ; i++) {
      pt_ann[i] = 0.;
    }

    if (fUseConstituents) {
      TIter itConstituents(candidate->GetCandidates());
      while((constituent = static_cast<Candidate*>(itConstituents.Next()))) {
        float pt = constituent->Momentum.Pt();
        float dr = candidate->Momentum.DeltaR(constituent->Momentum);
	//	cout << " There exists a constituent with dr=" << dr << endl;
	sumpt += pt;
	sumdrsqptsq += dr*dr*pt*pt;
	sumptsq += pt*pt;
	if (constituent->Charge == 0) {
	  nn++;
	} else {
	  if (constituent->IsRecoPU) {
	    sumptchpu += pt;
	  } else {
	    sumptchpv += pt;
	  }
	  sumptch += pt;
	  nc++;
	}
	for (int i = 0 ; i < 5 ; i++) {
	  if (dr > 0.1*i && dr < 0.1*(i+1)) {
	    pt_ann[i] += pt;
	  }
	}
	float tow_sumT = 0;
	float tow_sumW = 0;
	for (int i = 0 ; i < constituent->ecal_E_t.size() ; i++) {
	  float w = TMath::Sqrt(constituent->ecal_E_t[i].first);
	  if (fAverageEachTower) {
            tow_sumT += w*constituent->ecal_E_t[i].second;
            tow_sumW += w;
	  } else {
	    sumT0 += w*constituent->ecal_E_t[i].second;
	    sumT1 += w*gRandom->Gaus(constituent->ecal_E_t[i].second,0.001);
	    sumT10 += w*gRandom->Gaus(constituent->ecal_E_t[i].second,0.010);
	    sumT20 += w*gRandom->Gaus(constituent->ecal_E_t[i].second,0.020);
	    sumT30 += w*gRandom->Gaus(constituent->ecal_E_t[i].second,0.030);
	    sumT40 += w*gRandom->Gaus(constituent->ecal_E_t[i].second,0.040);
	    sumWeightsForT += w;
	    candidate->nTimes++;
	  }
	}
	if (fAverageEachTower && tow_sumW > 0.) {
	  sumT0 += tow_sumT;
	  sumT1 += tow_sumW*gRandom->Gaus(tow_sumT/tow_sumW,0.001);
          sumT10 += tow_sumW*gRandom->Gaus(tow_sumT/tow_sumW,0.0010);
          sumT20 += tow_sumW*gRandom->Gaus(tow_sumT/tow_sumW,0.0020);
          sumT30 += tow_sumW*gRandom->Gaus(tow_sumT/tow_sumW,0.0030);
          sumT40 += tow_sumW*gRandom->Gaus(tow_sumT/tow_sumW,0.0040);
	  sumWeightsForT += tow_sumW;
	  candidate->nTimes++;
	}
      }
    } else {
      // Not using constituents, using dr
      fItTrackInputArray->Reset();
      while ((trk = static_cast<Candidate*>(fItTrackInputArray->Next()))) {
	if (trk->Momentum.DeltaR(candidate->Momentum) < fParameterR) {
	  float pt = trk->Momentum.Pt();
	  sumpt += pt;
	  sumptch += pt;
	  if (trk->IsRecoPU) {
	    sumptchpu += pt;
	  } else {
	    sumptchpv += pt;
	  }
	  float dr = candidate->Momentum.DeltaR(trk->Momentum);
	  sumdrsqptsq += dr*dr*pt*pt;
	  sumptsq += pt*pt;
	  nc++;
	  for (int i = 0 ; i < 5 ; i++) {
	    if (dr > 0.1*i && dr < 0.1*(i+1)) {
	      pt_ann[i] += pt;
	    }
	  }
	}
      }
      fItNeutralInputArray->Reset();
      while ((constituent = static_cast<Candidate*>(fItNeutralInputArray->Next()))) {
	if (constituent->Momentum.DeltaR(candidate->Momentum) < fParameterR) {
	  float pt = constituent->Momentum.Pt();
	  sumpt += pt;
	  float dr = candidate->Momentum.DeltaR(constituent->Momentum);
	  sumdrsqptsq += dr*dr*pt*pt;
	  sumptsq += pt*pt;
	  nn++;
	  for (int i = 0 ; i < 5 ; i++) {
	    if (dr > 0.1*i && dr < 0.1*(i+1)) {
            pt_ann[i] += pt;
	    }
	  }
	}
      }
    }

    if (sumptch > 0.) {
      candidate->Beta = sumptchpv/sumptch;
      candidate->BetaStar = sumptchpu/sumptch;
    } else {
      candidate->Beta = -999.;
      candidate->BetaStar = -999.;
    }
    if (sumptsq > 0.) {
      candidate->MeanSqDeltaR = sumdrsqptsq/sumptsq;
    } else {
      candidate->MeanSqDeltaR = -999.;
    }
    candidate->NCharged = nc;
    candidate->NNeutrals = nn;
    if (sumpt > 0.) {
      candidate->PTD = TMath::Sqrt(sumptsq) / sumpt;
      for (int i = 0 ; i < 5 ; i++) {
        candidate->FracPt[i] = pt_ann[i]/sumpt;
      }
    } else {
      candidate->PTD = -999.;
      for (int i = 0 ; i < 5 ; i++) {
        candidate->FracPt[i] = -999.;
      }
    }

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
