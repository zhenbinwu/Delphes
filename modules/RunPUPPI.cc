/** \class RunPUPPI
 *
 *  Finds jets using FastJet library.
 *
 *  $Date: 2013-05-24 01:26:33 +0200 (Fri, 24 May 2013) $
 *  $Revision: 1122 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *  \with modifications by J. Stupak and J. Dolen
 *
 */

#include "PUPPI/puppiCleanContainer.hh"
#include "PUPPI/RecoObj.hh"

#include "modules/RunPUPPI.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

//#include "ExRootAnalysis/ExRootResult.h"
//#include "ExRootAnalysis/ExRootFilter.h"
//#include "ExRootAnalysis/ExRootClassifier.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

//------------------------------------------------------------------------------

RunPUPPI::RunPUPPI() :
  fItTrackInputArray(0), fItNeutralInputArray(0)
{

}

//------------------------------------------------------------------------------

RunPUPPI::~RunPUPPI()
{

}

//------------------------------------------------------------------------------

void RunPUPPI::Init()
{
  // define algorithm

  fTrackerEta = GetDouble("TrackerEta", 2.5);

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "Calorimeter/towers"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();
  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "Calorimeter/towers"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();


  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "weightedparticles"));
}

//------------------------------------------------------------------------------

void RunPUPPI::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;
}

//------------------------------------------------------------------------------

void RunPUPPI::Process()
{
  //  Candidate *candidate, *constituent;
  Candidate *candidate;
  TLorentzVector momentum;
  //  Double_t deta, dphi, detaMax, dphiMax;
  //  Int_t number;
  //  Double_t rho = 0;
  //  PseudoJet jet, area;

  DelphesFactory *factory = GetFactory();

  // loop over input objects
  fItTrackInputArray->Reset();
  fItNeutralInputArray->Reset();

  vector<RecoObj> puppiInputVector;



  while((candidate = static_cast<Candidate*>(fItTrackInputArray->Next())))
    {   
      momentum = candidate->Momentum;
      RecoObj curPseudoJet;
      curPseudoJet.pt  = momentum.M();
      curPseudoJet.eta = momentum.Eta();
      curPseudoJet.phi = momentum.Phi();
      curPseudoJet.m   = momentum.M();
      if (candidate->IsRecoPU) {
	curPseudoJet.id  = 3;
      } else {
	curPseudoJet.id  = 2;
      }
      puppiInputVector.push_back(curPseudoJet);
    }
  while((candidate = static_cast<Candidate*>(fItNeutralInputArray->Next())))
    {
      momentum = candidate->Momentum;
      RecoObj curPseudoJet;
      curPseudoJet.pt  = momentum.M();
      curPseudoJet.eta = momentum.Eta();
      curPseudoJet.phi = momentum.Phi();
      curPseudoJet.m   = momentum.M();
      curPseudoJet.id  = 1;
      puppiInputVector.push_back(curPseudoJet);
    }

  puppiCleanContainer curEvent(puppiInputVector,fTrackerEta);
  std::vector < fastjet::PseudoJet > puppiParticles = curEvent.puppiEvent(7,0.5);

  for (std::vector<fastjet::PseudoJet>::iterator it = puppiParticles.begin() ; it != puppiParticles.end() ; it++) {
    candidate = factory->NewCandidate();
    candidate->Momentum.SetXYZT(it->px(),it->py(),it->pz(),it->e());
    fOutputArray->Add(candidate);
  }
}
