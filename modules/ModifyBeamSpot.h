#ifndef ModifyBeamSpot_h
#define ModifyBeamSpot_h

/** \class ModifyBeamSpot
 *
 *  \author S. Zenz
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class ModifyBeamSpot: public DelphesModule
{
public:

  ModifyBeamSpot();
  ~ModifyBeamSpot();

  void Init();
  void Process();
  void Finish();

private:

  DelphesFormula *fFormula; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  
  TObjArray *fOutputArray; //!

  Double_t fZVertexSpread;
  Double_t currentZ, currentT;
  Double_t currentPU;

  // Store Z of PV
  TObjArray *fPVOutputArray; //!

  ClassDef(ModifyBeamSpot, 1)
};

#endif
