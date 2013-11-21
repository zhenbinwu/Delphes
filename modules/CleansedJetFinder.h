#ifndef CleansedJetFinder_h
#define CleansedJetFinder_h

/** \class CleansedJetFinder
 *
 *  Finds jets using FastJet library.
 *
 *  $Date: 2013-03-06 18:06:55 +0100 (Wed, 06 Mar 2013) $
 *  $Revision: 1031 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <vector>

class TObjArray;
class TIterator;

namespace fastjet {
  class JetDefinition;
  class AreaDefinition;
  class Selector;
  
  namespace contrib {
    class JetCleanser;
  }
}

class CleansedJetFinder: public DelphesModule
{
public:

  CleansedJetFinder();
  ~CleansedJetFinder();

  void Init();
  void Process();
  void Finish();
  
private:

  void *fPlugin; //!
  fastjet::JetDefinition *fDefinition; //!
   
  Int_t fJetAlgorithm;
  Double_t fParameterR;
  Double_t fJetPTMin;
  Double_t fConeRadius;
  Double_t fSeedThreshold;
  Double_t fConeAreaFraction;
  Int_t fMaxIterations;
  Int_t fMaxPairSize;
  Int_t fIratch;
  Double_t fAdjacencyCut;
  Double_t fOverlapThreshold;
  
  // --- FastJet Area method --------
  
  fastjet::AreaDefinition *fAreaDefinition;
  Int_t fAreaAlgorithm;  
  //  Bool_t  fComputeRho; 
  //  Double_t fRhoEtaMax;
  
  // -- ghost based areas --
  Double_t fGhostEtaMax;
  Int_t fRepeat;
  Double_t fGhostArea;
  Double_t fGridScatter;
  Double_t fPtScatter;
  Double_t fMeanGhostPt; 
  
  // -- voronoi areas --
  Double_t fEffectiveRfact;

  // -- cleansing parameters --
  Double_t fLinearGamma;
  Double_t fTrimmingParameter;
  Int_t fFilteringParameter;
  Int_t fSubjetAlgorithm;
  Double_t fSubjetParameterR;

  TIterator *fItChargedPrimaryInputArray; //!
  TIterator *fItChargedPileUpInputArray; //!
  TIterator *fItNeutralInputArray; //!

  const TObjArray *fChargedPrimaryInputArray; //!
  const TObjArray *fChargedPileUpInputArray; //!
  const TObjArray *fNeutralInputArray; //!

  TObjArray *fPlainOutputArray; //!
  TObjArray *fJVFOutputArray; // SCZ
  TObjArray *fLinearOutputArray; // SCZ
  TObjArray *fLinearTrimmedOutputArray; // SCZ
  TObjArray *fLinearFilteredOutputArray; // SCZ
  TObjArray *fGaussianOutputArray; // SCZ

  Bool_t fSavePlain;
  Bool_t fSaveJVF;
  Bool_t fSaveLinear;
  Bool_t fSaveLinearTrimmed;
  Bool_t fSaveLinearFiltered;
  Bool_t fSaveGaussian;

  fastjet::JetDefinition *fSubjetDefinition; //!

  fastjet::contrib::JetCleanser *jvf_cleanser_B;
  fastjet::contrib::JetCleanser *linear_cleanser_B;
  fastjet::contrib::JetCleanser *linear_cleanser_B_trimmed;
  fastjet::contrib::JetCleanser *linear_cleanser_B_filtered;
  fastjet::contrib::JetCleanser *gaussian_cleanser_B;
 

  /*
  modules/CleansedJetFinder.cc:305:38: error: 'jvf_cleanser_B' was not declared in this scope
    modules/CleansedJetFinder.cc:316:38: error: 'linear_cleanser_B' was not declared in this scope
    modules/CleansedJetFinder.cc:327:38: error: 'gaussian_cleanser_B' was not declared in this scope
  */

  ClassDef(CleansedJetFinder, 1)
};

#endif
