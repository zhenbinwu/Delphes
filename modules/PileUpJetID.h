#ifndef PileUpJetID_h
#define PileUpJetID_h

/** \class PileUpJetID
 *
 *  Subtract pile-up contribution from jets using the fastjet area method
 *
 *  $Date: 2012-11-18 15:57:08 +0100 (Sun, 18 Nov 2012) $
 *  $Revision: 814 $
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;
class DelphesFormula;

class PileUpJetID: public DelphesModule
{
public:

  PileUpJetID();
  ~PileUpJetID();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fJetPTMin;
  Double_t fParameterR;

  // If set to true, may have weird results for PFCHS
  // If set to false, uses everything within dR < fParameterR even if in other jets &c.
  // Results should be very similar for PF
  Int_t fUseConstituents; 

  TIterator *fItJetInputArray; //!

  const TObjArray *fJetInputArray; //!

  const TObjArray *fTrackInputArray; // SCZ
  const TObjArray *fNeutralInputArray; 

  TIterator *fItTrackInputArray; // SCZ
  TIterator *fItNeutralInputArray; // SCZ

  TObjArray *fOutputArray; //!

  ClassDef(PileUpJetID, 1)
};

#endif
