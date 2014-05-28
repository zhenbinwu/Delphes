
/** \class TrackPileUpSubtractor
 *
 *  Subtract pile-up contribution from tracks.
 *
 *  $Date: 2013-03-24 21:25:01 +0100 (Sun, 24 Mar 2013) $
 *  $Revision: 1073 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/TrackPileUpSubtractor.h"

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

TrackPileUpSubtractor::TrackPileUpSubtractor():fZVertexResolutionFormula(0)
{
  fZVertexResolutionFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

TrackPileUpSubtractor::~TrackPileUpSubtractor()
{
  if(fZVertexResolutionFormula) delete fZVertexResolutionFormula;
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Init()
{
  fZVertexResolution  = GetDouble("ZVertexResolution", 0.005)*1.0E3;
  fZVertexResolutionFormula->Compile(GetString("ZVertexResolutionFormula", "1.0"));

  // import arrays with output from other modules

  ExRootConfParam param = GetParam("InputArray");
  Long_t i, size;
  const TObjArray *array;
  TIterator *iterator;

  size = param.GetSize();
  for(i = 0; i < size/2; ++i)
  {
    array = ImportArray(param[i*2].GetString());
    iterator = array->MakeIterator();

    fInputMap[iterator] = ExportArray(param[i*2 + 1].GetString());;
  }

  fPVInputArray = ImportArray(GetString("PVInputArray", "PV"));
  fPVItInputArray = fPVInputArray->MakeIterator();

}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Finish()
{
  map< TIterator *, TObjArray * >::iterator itInputMap;
  TIterator *iterator;

  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;

    if(iterator) delete iterator;
  }
}

//------------------------------------------------------------------------------

void TrackPileUpSubtractor::Process()
{
  Candidate *candidate, *particle;
  map< TIterator *, TObjArray * >::iterator itInputMap;
  TIterator *iterator;
  TObjArray *array;
  Double_t z;
  Double_t PVZ = 0.;

  fPVItInputArray->Reset();
  Candidate *pv = static_cast<Candidate*>(fPVItInputArray->Next());
  if (pv) PVZ = pv->Position.Z();
  //  cout << "SCZ TrackPileUpSubtractor Debug PVZ=" << PVZ << endl;

  // loop over all input arrays
  for(itInputMap = fInputMap.begin(); itInputMap != fInputMap.end(); ++itInputMap)
  {
    iterator = itInputMap->first;
    array = itInputMap->second;

    // loop over all candidates
    iterator->Reset();
    while((candidate = static_cast<Candidate*>(iterator->Next())))
    {
      particle = static_cast<Candidate*>(candidate->GetCandidates()->Last());
      z = particle->Position.Z();

      // apply pile-up subtraction
      // assume perfect pile-up subtraction for tracks outside fZVertexResolution
      //if(candidate->IsPU && TMath::Abs(z-PVZ) > fZVertexResolution) {
      if(TMath::Abs(z-PVZ) > fZVertexResolutionFormula->Eval(0, particle->Position.Eta())) {
        candidate->IsRecoPU = 1;
      } else {
        candidate->IsRecoPU = 0;
        array->Add(candidate);
      }
    }
  }
}

