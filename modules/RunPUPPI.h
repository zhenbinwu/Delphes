#ifndef RunPUPPI_h
#define RunPUPPI_h

/** \class RunPUPPI
 *
 * Runs PUPPI!
 *
 * SZ
 *
 */

#include "classes/DelphesModule.h"
#include <vector>

class TObjArray;
class TIterator;


class RunPUPPI: public DelphesModule
{
public:

  RunPUPPI();
  ~RunPUPPI();

  void Init();
  void Process();
  void Finish();
  
private:

  Double_t fTrackerEta;

  TIterator *fItTrackInputArray;
  TIterator *fItNeutralInputArray; //!

  const TObjArray *fTrackInputArray;
  const TObjArray *fNeutralInputArray; //!

  TObjArray *fOutputArray;

  ClassDef(RunPUPPI, 1)
};

#endif
