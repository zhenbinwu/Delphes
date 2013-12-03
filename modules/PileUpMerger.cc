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

#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Random/JamesRandom.h"

#include "modules/PileUpMerger.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

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

PileUpMerger::PileUpMerger() :
  fReader(0), fItInputArray(0)
{
}

//------------------------------------------------------------------------------

PileUpMerger::~PileUpMerger()
{
}

//------------------------------------------------------------------------------

double PileUpMerger::BetaFunction(double z, double z0)
{
  return sqrt(femittance*(fbetastar+(((z-z0)*(z-z0))/fbetastar)));
}

void PileUpMerger::Init()
{
  const char *fileName;

  fMeanPileUp  = GetDouble("MeanPileUp", 10);
  fZVertexSpread = GetDouble("ZVertexSpread", 0.05)*1.0E3;

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

343 NominalCollisionVtxSmearingParameters = cms.PSet(
344     Phi = cms.double(0.000142),
345     BetaStar = cms.double(55.0),
346     Emittance = cms.double(1.006e-07),
347     Alpha = cms.double(0.0),
348     SigmaZ = cms.double(5.3),
349     TimeOffset = cms.double(0.0),
350     Y0 = cms.double(0.0),
351     X0 = cms.double(0.05),
352     Z0 = cms.double(0.0)
353 )
  */

  // All the units need to be checked!!!!
  fTimeOffset = 0.0;
  phi_ = 0.000142;
  fbetastar = 55.;
  alpha_ = 0.0;
  fSigmaZ = 5.3;
  fY0 = 0.;
  fX0 = 0.05;
  fZ0 = 0.0;
  femittance = 1.006e-07;

  fEngine = new CLHEP::HepJamesRandom();
  fRandom = new CLHEP::RandGaussQ(fEngine); // (getEngine()); // effect on random seeding?

  fileName = GetString("PileUpFile", "MinBias.pileup");
  fReader = new DelphesPileUpReader(fileName);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
  fNPUOutputArray = ExportArray(GetString("NPUOutputArray", "NPU"));
}

//------------------------------------------------------------------------------

void PileUpMerger::Finish()
{
  if(fReader) delete fReader;
}

//------------------------------------------------------------------------------

void PileUpMerger::Process()
{
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle;
  Int_t pid;
  Float_t x, y, z, t;
  Float_t px, py, pz, e;
  Double_t dz, dphi;
  Int_t poisson, event;
  Long64_t allEntries, entry;
  Candidate *candidate;
  DelphesFactory *factory;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    fOutputArray->Add(candidate);
  }

  factory = GetFactory();

  poisson = gRandom->Poisson(fMeanPileUp);

  allEntries = fReader->GetEntries();

  for(event = 0; event < poisson; ++event)
  {
    do
    {
      entry = TMath::Nint(gRandom->Rndm()*allEntries);
    }
    while(entry >= allEntries);

    fReader->ReadEntry(entry);

    //    dz = gRandom->Gaus(0.0, fZVertexSpread);
    //    dphi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
    dz = 0.;
    dphi = 0.;

    double X,Y,Z;
   
    double tmp_sigz = fRandom->fire(0., fSigmaZ);
    Z = tmp_sigz + fZ0;

    double tmp_sigx = BetaFunction(Z,fZ0); 
    // need sqrt(2) for beamspot width relative to single beam width
    tmp_sigx /= sqrt(2.0);
    X = fRandom->fire(0.,tmp_sigx) + fX0; // + Z*fdxdz ;
    
    double tmp_sigy = BetaFunction(Z,fZ0);
    // need sqrt(2) for beamspot width relative to single beam width
    tmp_sigy /= sqrt(2.0);
    Y = fRandom->fire(0.,tmp_sigy) + fY0; // + Z*fdydz;
    
    double tmp_sigt = fRandom->fire(0., fSigmaZ);
    double T = tmp_sigt + fTimeOffset; 
    
    while(fReader->ReadParticle(pid, x, y, z, t, px, py, pz, e))
    {
      candidate = factory->NewCandidate();

      // Get rid of BS position in PU
      // Hard-coded to deal with http://red-gridftp11.unl.edu/Snowmass/MinBias100K_14TeV.pileup
      x = x - 2.44; 
      y = y - 3.93;

      candidate->PID = pid;

      candidate->Status = 1;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      candidate->IsPU = event+1; // might as well store which PU vertex this comes from so they can be separated

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum.RotateZ(dphi);

      //      candidate->Position.SetXYZT(x, y, z + dz, t);
      //      candidate->Position.SetXYZT(X,Y,Z,T);
      candidate->Position.SetXYZT(x,y,Z,T); // Use CMSSW Z,T - keep old 
      candidate->Position.RotateZ(dphi);

      // 0.24, 0.39
      //      candidate->Position.

      fOutputArray->Add(candidate);
    }
  }
  candidate = factory->NewCandidate();
  candidate->Momentum.SetPtEtaPhiE((float)poisson, 0.0, 0.0, (float)poisson); // cheating and storing NPU as a float
  fNPUOutputArray->Add(candidate);

}

//------------------------------------------------------------------------------

