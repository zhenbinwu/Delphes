
/**
 *
 *  Definition of classes to be stored in the root tree.
 *  Function CompareXYZ sorts objects by the variable XYZ that MUST be
 *  present in the data members of the root tree class of the branch.
 *
 *  $Date: 2008-06-04 13:57:24 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"

#include "classes/DelphesFactory.h"
#include "classes/SortableObject.h"

CompBase *GenParticle::fgCompare = 0;
CompBase *Photon::fgCompare = CompPT<Photon>::Instance();
CompBase *Electron::fgCompare = CompPT<Electron>::Instance();
CompBase *Muon::fgCompare = CompPT<Muon>::Instance();
CompBase *IsoTrack::fgCompare = CompPT<IsoTrack>::Instance();
CompBase *Jet::fgCompare = CompPT<Jet>::Instance();
CompBase *Track::fgCompare = CompPT<Track>::Instance();
CompBase *Tower::fgCompare = CompE<Tower>::Instance();
CompBase *Candidate::fgCompare = CompMomentumPt<Candidate>::Instance();

//------------------------------------------------------------------------------

TLorentzVector GenParticle::P4()
{
  TLorentzVector vec;
  vec.SetPxPyPzE(Px, Py, Pz, E);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Photon::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Electron::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Muon::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector IsoTrack::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Jet::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

TLorentzVector Jet::AreaP4()
{
  TLorentzVector vec(AreaX,AreaY,AreaZ,AreaT);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Track::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Tower::P4()
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(ET, Eta, Phi, 0.0);
  return vec;
}


//------------------------------------------------------------------------------

Candidate::Candidate() :
  PID(0), Status(0), M1(-1), M2(-1), D1(-1), D2(-1),
  IsolationVar(0),
  Charge(0), Mass(0.0),
  IsPU(0), IsRecoPU(0), IsEMCand(0), IsConstituent(0),
  BTag(0), TauTag(0), Eem(0.0), Ehad(0.0),
  WTag(0), TopTag(0), HTag(0),
  Tau1(-999), Tau2(-999), Tau3(-999),
  NSubJets(-999), MassDrop(-999), TrimmedMass(-999),
  Beta(-999.), BetaStar(-999.), 
  MeanSqDeltaR(-999.), PTD(-999.),
  NCharged(-999), NNeutrals(-999),
  DeltaEta(0.0), DeltaPhi(0.0),
  Momentum(0.0, 0.0, 0.0, 0.0),
  Position(0.0, 0.0, 0.0, 0.0),
  Area(0.0, 0.0, 0.0, 0.0),
  //  t0(-999999.),
  //  t1(-999999.),
  //  t10(-999999.),
  //  t20(-999999.),
  //  t30(-999999.),
  //  t40(-999999.),
  //  nTimes(-1),
  fFactory(0),
  fArray(0)
{
  Edges[0] = 0.0;
  Edges[1] = 0.0;
  Edges[2] = 0.0;
  Edges[3] = 0.0;
  for (int i = 0 ; i < 5 ; i++) {
    FracPt[i] = -999.;
  }
}

//------------------------------------------------------------------------------

void Candidate::AddCandidate(Candidate *object)
{
  if(!fArray) fArray = fFactory->NewArray();
  fArray->Add(object);
}

//------------------------------------------------------------------------------

TObjArray *Candidate::GetCandidates()
{
  if(!fArray) fArray = fFactory->NewArray();
  return fArray;
}

//------------------------------------------------------------------------------

Bool_t Candidate::Overlaps(const Candidate *object) const
{
  const Candidate *candidate;

  if(object->GetUniqueID() == GetUniqueID()) return kTRUE;

  if(fArray)
  {
    TIter it(fArray);
    while((candidate = static_cast<Candidate *>(it.Next())))
    {
      if(candidate->Overlaps(object)) return kTRUE;
    }
  }

  if(object->fArray)
  {
    TIter it(object->fArray);
    while((candidate = static_cast<Candidate *>(it.Next())))
    {
      if(candidate->Overlaps(this)) return kTRUE;
    }
  }

  return kFALSE;
}


//------------------------------------------------------------------------------

TObject *Candidate::Clone(const char *newname) const
{
  Candidate *object = fFactory->NewCandidate();
  Copy(*object);
  return object;
}

//------------------------------------------------------------------------------

void Candidate::Copy(TObject &obj) const
{
  Candidate &object = static_cast<Candidate &>(obj);
  Candidate *candidate;

  object.PID = PID;
  object.Status = Status;
  object.IsolationVar = IsolationVar;
  object.M1 = M1;
  object.M2 = M2;
  object.D1 = D1;
  object.D2 = D2;
  object.Charge = Charge;
  object.Mass = Mass;
  object.IsPU = IsPU;
  object.IsRecoPU = IsRecoPU;
  object.IsEMCand = IsEMCand;
  object.IsConstituent = IsConstituent;
  object.BTag = BTag;
  object.TauTag = TauTag;
  object.WTag = WTag;
  object.TopTag = TopTag;
  object.HTag = HTag;
  object.Tau1=Tau1; 
  object.Tau2=Tau2; 
  object.Tau3=Tau3;
  object.NSubJets=NSubJets;
  object.MassDrop=MassDrop;
  object.TrimmedMass=TrimmedMass;
  object.Eem = Eem;
  object.Ehad = Ehad;
  object.Edges[0] = Edges[0];
  object.Edges[1] = Edges[1];
  object.Edges[2] = Edges[2];
  object.Edges[3] = Edges[3];
  object.DeltaEta = DeltaEta;
  object.DeltaPhi = DeltaPhi;
  object.Momentum = Momentum;
  object.Position = Position;
  object.Area = Area;
  object.Beta = Beta;
  object.BetaStar = BetaStar;
  object.MeanSqDeltaR = MeanSqDeltaR;
  object.PTD = PTD;
  object.NCharged = NCharged; 
  object.NNeutrals = NNeutrals;
  for (int i = 0 ; i < 5 ; i++) {
    object.FracPt[i] = FracPt[i];
  }

  // Copy jet timing info
  //  object.t0 = t0;
  //  object.t1 = t1;
  // object.t10 = t10;
  // object.t20 = t20;
  // object.t30 = t30;
  // object.t40 = t40;
  // object.nTimes = nTimes;

  // Copy cluster timing info
  copy(ecal_E_t.begin(),ecal_E_t.end(),back_inserter(object.ecal_E_t));

  object.fFactory = fFactory;
  object.fArray = 0;

  if(fArray && fArray->GetEntriesFast() > 0)
  {
    TIter itArray(fArray);
    TObjArray *array = object.GetCandidates();
    while((candidate = static_cast<Candidate *>(itArray.Next())))
    {
      array->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------

void Candidate::Clear(Option_t* option)
{
  SetUniqueID(0);
  ResetBit(kIsReferenced);
  PID = 0;
  Status = 0;
  IsolationVar =0.;
  M1 = -1; M2 = -1; D1 = -1; D2 = -1;
  Charge = 0;
  Mass = 0.0;
  IsPU = 0;
  IsRecoPU = 0;
  IsConstituent = 0;
  IsEMCand = 0;
  BTag = 0;
  TauTag = 0;
  WTag = 0;
  TopTag = 0;
  HTag = 0;
  Tau1 = -999;
  Tau2 = -999;
  Tau3 = -999;
  NSubJets=-999;
  MassDrop = -999;
  TrimmedMass = -999;
  Eem = 0.0;
  Ehad = 0.0;
  Edges[0] = 0.0;
  Edges[1] = 0.0;
  Edges[2] = 0.0;
  Edges[3] = 0.0;
  DeltaEta = 0.0;
  DeltaPhi = 0.0;
  Momentum.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Area.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Beta = -999.;
  BetaStar = -999.;
  MeanSqDeltaR = -999.;
  PTD = -999.;
  NCharged = -1;
  NNeutrals = -1;
  fArray = 0;

  ecal_E_t.clear();
}
