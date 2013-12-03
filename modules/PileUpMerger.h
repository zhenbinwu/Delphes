#ifndef PileUpMerger_h
#define PileUpMerger_h

/** \class PileUpMerger
 *
 *  Merges particles from pile-up sample into event
 *
 *
 *  $Date: 2013-02-12 15:13:59 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 907 $
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

namespace CLHEP {
  class RandGaussQ;
  class HepRandomEngine;
}

class TObjArray;
class DelphesPileUpReader;

class PileUpMerger: public DelphesModule
{
public:

  PileUpMerger();
  ~PileUpMerger();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fMeanPileUp;
  Double_t fZVertexSpread;

  /*
048     fX0 =        p.getParameter<double>("X0")*cm;
049     fY0 =        p.getParameter<double>("Y0")*cm;
050     fZ0 =        p.getParameter<double>("Z0")*cm;
051     fSigmaZ =    p.getParameter<double>("SigmaZ")*cm;
052     alpha_ =     p.getParameter<double>("Alpha")*radian;
053     phi_ =       p.getParameter<double>("Phi")*radian;
054     fbetastar =  p.getParameter<double>("BetaStar")*cm;
055     femittance = p.getParameter<double>("Emittance")*cm; // this is not the normalized emittance
056     fTimeOffset = p.getParameter<double>("TimeOffset")*ns*c_light; // HepMC time units are mm
  */

  Double_t phi_;
  Double_t fbetastar;
  Double_t fX0,fY0,fZ0;
  Double_t fSigmaZ;
  Double_t alpha_;
  Double_t femittance;
  Double_t fTimeOffset;
  CLHEP::RandGaussQ*  fRandom ;
  CLHEP::HepRandomEngine* fEngine;

  DelphesPileUpReader *fReader;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  double BetaFunction(double z, double z0);

  TObjArray *fNPUOutputArray; //!                                                                                                                                                    
  ClassDef(PileUpMerger, 2)
};

#endif
