#ifndef PileUpJetID_h
#define PileUpJetID_h

/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables
 *
 *  \author S. Zenz
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
