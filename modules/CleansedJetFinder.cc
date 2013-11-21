/** \class CleansedJetFinder
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

#include "fastjet/plugins/JetCleanser/JetCleanser.hh"

#include "modules/CleansedJetFinder.h"

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
#include <vector>

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Filter.hh"

#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"

// For some reason, including this here and in FastJetFinder breaks everything
//#include "fastjet/tools/Njettiness.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//------------------------------------------------------------------------------

CleansedJetFinder::CleansedJetFinder() :
  fPlugin(0), fDefinition(0), fAreaDefinition(0), fItChargedPrimaryInputArray(0), fItChargedPileUpInputArray(0), fItNeutralInputArray(0), fSubjetDefinition(0)
{

}

//------------------------------------------------------------------------------

CleansedJetFinder::~CleansedJetFinder()
{

}

//------------------------------------------------------------------------------

void CleansedJetFinder::Init()
{
  JetDefinition::Plugin *plugin = NULL;

  // define algorithm

  fJetAlgorithm = GetInt("JetAlgorithm", 6);
  fParameterR = GetDouble("ParameterR", 0.5);

  fConeRadius = GetDouble("ConeRadius", 0.5);
  fSeedThreshold = GetDouble("SeedThreshold", 1.0);
  fConeAreaFraction = GetDouble("ConeAreaFraction", 1.0);
  fMaxIterations = GetInt("MaxIterations", 100);
  fMaxPairSize = GetInt("MaxPairSize", 2);
  fIratch = GetInt("Iratch", 1);
  fAdjacencyCut = GetDouble("AdjacencyCut", 2.0);
  fOverlapThreshold = GetDouble("OverlapThreshold", 0.75);

  fJetPTMin = GetDouble("JetPTMin", 10.0);

  // ---  Jet Area Parameters ---
  fAreaAlgorithm = GetInt("AreaAlgorithm", 0);
  //  fComputeRho = GetBool("ComputeRho", false);
  //  fRhoEtaMax = GetDouble("RhoEtaMax", 5.0);
  // - ghost based areas -
  fGhostEtaMax = GetDouble("GhostEtaMax", 5.0);
  fRepeat = GetInt("Repeat", 1);
  fGhostArea = GetDouble("GhostArea", 0.01);
  fGridScatter = GetDouble("GridScatter", 1.0);
  fPtScatter = GetDouble("PtScatter", 0.1);
  fMeanGhostPt = GetDouble("MeanGhostPt", 1.0E-100);
  // - voronoi based areas -
  fEffectiveRfact = GetDouble("EffectiveRfact", 1.0);

  // - cleansing parameters -
  fSavePlain = GetBool("SavePlain",true);
  fSaveJVF = GetBool("SaveJVF",true);
  fSaveLinear = GetBool("SaveLinear",true);
  fSaveLinearTrimmed = GetBool("SaveLinearTrimmed",true);
  fSaveLinearFiltered = GetBool("SaveLinearFiltered",true);
  fSaveGaussian = GetBool("SaveGaussian",false);

  fLinearGamma = GetDouble("LinearGamma",0.55);
  fTrimmingParameter = GetDouble("TrimmingParameter",0.05);
  fSubjetParameterR = GetDouble("SubjetRParameter",0.3);
  fFilteringParameter = GetInt("FilteringParameter",2);
  fSubjetAlgorithm = GetInt("SubjetAlgorithm",4);
  

  /*
  // -- cleansing parameters --
  Double_t fLinearGamma;
  Double_t fTrimmingParam;
  Int_t fFilteringParam;
  Int_t fSubjetAlgorithm;
  Double_t fSubjetParameterR;


  Bool_t fSavePlain;
  Bool_t fSaveJVF;
  Bool_t fSaveLinear;
  Bool_t fSaveGaussian;
  */


  switch(fAreaAlgorithm)
  {
    case 1:
      fAreaDefinition = new fastjet::AreaDefinition(active_area_explicit_ghosts, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
    case 2:
      fAreaDefinition = new fastjet::AreaDefinition(one_ghost_passive_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
    case 3:
      fAreaDefinition = new fastjet::AreaDefinition(passive_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
    case 4:
      fAreaDefinition = new fastjet::AreaDefinition(VoronoiAreaSpec(fEffectiveRfact));
      break;
    case 5:
      fAreaDefinition = new fastjet::AreaDefinition(active_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
    default:
    case 0:
      fAreaDefinition = 0;
      break;
  }

  switch(fJetAlgorithm)
  {
    case 1: 
      plugin = new fastjet::CDFJetCluPlugin(fSeedThreshold, fConeRadius, fAdjacencyCut, fMaxIterations, fIratch, fOverlapThreshold);
      fDefinition = new fastjet::JetDefinition(plugin);
      break;
    case 2:
      plugin = new fastjet::CDFMidPointPlugin(fSeedThreshold, fConeRadius, fConeAreaFraction, fMaxPairSize, fMaxIterations, fOverlapThreshold);
      fDefinition = new fastjet::JetDefinition(plugin);
      break;
    case 3:
      plugin = new fastjet::SISConePlugin(fConeRadius, fOverlapThreshold, fMaxIterations, fJetPTMin);
      fDefinition = new fastjet::JetDefinition(plugin);
      break;
    case 4:
      fDefinition = new fastjet::JetDefinition(fastjet::kt_algorithm, fParameterR);
      break;
    case 5:
      fDefinition = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fParameterR);
      break;
    default:
    case 6:
      fDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm, fParameterR);
      break;
  }

  switch(fSubjetAlgorithm) 
  {
    default:
    case 4:
      fSubjetDefinition = new fastjet::JetDefinition(fastjet::kt_algorithm, fSubjetParameterR);
      break;
    case 5:
      fSubjetDefinition = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fSubjetParameterR);
      break;
    case 6:
      fSubjetDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm, fSubjetParameterR);
      break;
  }


  jvf_cleanser_B = new JetCleanser(*fSubjetDefinition, JetCleanser::jvf_cleansing, JetCleanser::input_nc_separate);

  linear_cleanser_B = new JetCleanser(*fSubjetDefinition, JetCleanser::linear_cleansing, JetCleanser::input_nc_separate);
  linear_cleanser_B->SetLinearParameters(fLinearGamma);

  linear_cleanser_B_trimmed = new JetCleanser(*linear_cleanser_B);
  linear_cleanser_B_trimmed->SetTrimming(fTrimmingParameter);

  linear_cleanser_B_filtered = new JetCleanser(*linear_cleanser_B);
  linear_cleanser_B_filtered->SetFiltering(fFilteringParameter);

  gaussian_cleanser_B = new JetCleanser(*fSubjetDefinition, JetCleanser::gaussian_cleansing, JetCleanser::input_nc_separate);
  gaussian_cleanser_B->SetGaussianParameters(0.67,0.62,0.20,0.25);

  
  fPlugin = plugin;

  ClusterSequence::print_banner();
  
  // import input array

  //  set ChargedPrimaryInputArray EFlowMergerCharged/eflow
  //  set ChargedPileUpInputArray TrackPileUpSubtractorRejected/eflowTracks
  //  set NeutralInputArray Calorimeter/eflowTowers

  fChargedPrimaryInputArray = ImportArray(GetString("ChargedPrimaryInputArray", "Calorimeter/towers"));
  fItChargedPrimaryInputArray = fChargedPrimaryInputArray->MakeIterator();
  fChargedPileUpInputArray = ImportArray(GetString("ChargedPileUpInputArray", "Calorimeter/towers"));
  fItChargedPileUpInputArray = fChargedPileUpInputArray->MakeIterator();
  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "Calorimeter/towers"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();


  // create output arrays
  fPlainOutputArray = ExportArray(GetString("PlainOutputArray", "jets"));
  fJVFOutputArray = ExportArray(GetString("JVFOutputArray", "jets"));
  fLinearOutputArray = ExportArray(GetString("LinearOutputArray", "jets"));
  fLinearTrimmedOutputArray = ExportArray(GetString("LinearTrimmedOutputArray", "jets"));
  fLinearFilteredOutputArray = ExportArray(GetString("LinearFilteredOutputArray", "jets"));
  fGaussianOutputArray = ExportArray(GetString("GaussianOutputArray", "jets"));

}

//------------------------------------------------------------------------------

void CleansedJetFinder::Finish()
{
  if(fItChargedPrimaryInputArray) delete fItChargedPrimaryInputArray;
  if(fItChargedPileUpInputArray) delete fItChargedPileUpInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;
  if(jvf_cleanser_B) delete jvf_cleanser_B;
  if(linear_cleanser_B) delete linear_cleanser_B;
  if(linear_cleanser_B_trimmed) delete linear_cleanser_B_trimmed;
  if(linear_cleanser_B_filtered) delete linear_cleanser_B_filtered;
  if(gaussian_cleanser_B) delete gaussian_cleanser_B;
  if(fSubjetDefinition) delete fSubjetDefinition;
  if(fDefinition) delete fDefinition;
  if(fAreaDefinition) delete fAreaDefinition;
  if(fPlugin) delete static_cast<JetDefinition::Plugin*>(fPlugin);
}

//------------------------------------------------------------------------------

void CleansedJetFinder::Process()
{
  Candidate *candidate, *constituent;
  TLorentzVector momentum;
  Double_t deta, dphi, detaMax, dphiMax;
  Int_t number;
  Double_t rho = 0;
  PseudoJet jet, area;
  vector<PseudoJet> hard_event_charged, pileup_charged;
  vector<PseudoJet> full_event;  // via "calo cells"
  vector<PseudoJet> full_event_neutral;      // via "particle flow"
  //  ClusterSequence *sequence;

  DelphesFactory *factory = GetFactory();


  // loop over input objects
  fItChargedPrimaryInputArray->Reset();
  fItChargedPileUpInputArray->Reset();
  fItNeutralInputArray->Reset();
  number = 0;
  while((candidate = static_cast<Candidate*>(fItChargedPrimaryInputArray->Next())))
  {   
    momentum = candidate->Momentum;
    jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
    jet.set_user_index(number);
    full_event.push_back(jet);
    hard_event_charged.push_back(jet);
    //    cout << "pushed back to inputlist pt eta phi m: " << momentum.Pt() << " " << momentum.Eta() << " " << momentum.Phi() << " " << momentum.M() << endl;
    ++number;
  }
  while((candidate = static_cast<Candidate*>(fItChargedPileUpInputArray->Next())))
    {
      momentum = candidate->Momentum;
      jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
      jet.set_user_index(number);
      full_event.push_back(jet);
      pileup_charged.push_back(jet);
      //    cout << "pushed back to inputlist pt eta phi m: " << momentum.Pt() << " " << momentum.Eta() << " " << momentum.Phi() << " " << momentum.M() << endl;
      ++number;
    }
  while((candidate = static_cast<Candidate*>(fItNeutralInputArray->Next())))
    {
      momentum = candidate->Momentum;
      jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
      jet.set_user_index(number);
      full_event.push_back(jet);
      full_event_neutral.push_back(jet);
      //    cout << "pushed back to inputlist pt eta phi m: " << momentum.Pt() << " " << momentum.Eta() << " " << momentum.Phi() << " " << momentum.M() << endl;
      ++number;
    }


  vector< vector<PseudoJet> > sets;
  sets.push_back( full_event );           // calorimeter cells
  sets.push_back( hard_event_charged );   // tracks from primary interaction
  sets.push_back( pileup_charged );       // tracks from pileup
  sets.push_back( full_event_neutral );   // neutral particles

  // collect jets
                                                                                                                                                                 
  vector< vector<PseudoJet> > jet_sets = ClusterSets(*fDefinition, full_event, sets, fJetPTMin);
  vector<PseudoJet> jets_plain     = jet_sets[0];
  vector<PseudoJet> jets_tracks_LV = jet_sets[1];
  vector<PseudoJet> jets_tracks_PU = jet_sets[2];
  vector<PseudoJet> jets_neutrals  = jet_sets[3];
  //  vector<PseudoJet> jets_plain;
  //  vector<PseudoJet> jets_tracks_LV;
  //  vector<PseudoJet> jets_tracks_PU;
  //  vector<PseudoJet> jets_neutrals;

  //  PseudoJet plain_dijet, truth_dijet, jvf_cleansed_dijet, lin_cleansed_dijet, gau_cleansed_dijet;
  //  unsigned n_jets;

  // ---------
  // Cleansing
  // ---------

  for (unsigned i=0; i<jets_plain.size(); i++) {
    PseudoJet plain_jet = jets_plain[i];

    if (fSavePlain) {
      if (plain_jet.pt() > fJetPTMin) {
	momentum.SetPxPyPzE(plain_jet.px(), plain_jet.py(), plain_jet.pz(), plain_jet.E());
	candidate = factory->NewCandidate();
	candidate->Momentum = momentum;
	candidate->Area.SetPxPyPzE(0.,0.,0.,0.);
	fPlainOutputArray->Add(candidate);
      }
    }
    if (fSaveJVF) {
      PseudoJet jvf_cleansed_jet = (*jvf_cleanser_B)( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(),
						   jets_tracks_PU[i].constituents() );
      if (jvf_cleansed_jet.pt() > fJetPTMin) {
	momentum.SetPxPyPzE(jvf_cleansed_jet.px(), jvf_cleansed_jet.py(), jvf_cleansed_jet.pz(), jvf_cleansed_jet.E());
	candidate = factory->NewCandidate();
	candidate->Momentum = momentum;
	candidate->Area.SetPxPyPzE(0.,0.,0.,0.);
	fJVFOutputArray->Add(candidate);
      }
    }
    if (fSaveLinear) {
      PseudoJet lin_cleansed_jet = (*linear_cleanser_B)( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(),
						      jets_tracks_PU[i].constituents() );
      if (lin_cleansed_jet.pt() > fJetPTMin) {
	momentum.SetPxPyPzE(lin_cleansed_jet.px(), lin_cleansed_jet.py(), lin_cleansed_jet.pz(), lin_cleansed_jet.E());
	candidate = factory->NewCandidate();
	candidate->Momentum = momentum;
	candidate->Area.SetPxPyPzE(0.,0.,0.,0.);
	fLinearOutputArray->Add(candidate);
      }
    }
    /*
      if (lin_cleansed_jet.pt() > 200. && fParameterR > 0.65) {
	cout << "SCZ No grooming pt eta phi m " << lin_cleansed_jet.pt() << " " << lin_cleansed_jet.eta() << " " << lin_cleansed_jet.phi() << " " << lin_cleansed_jet.m() << endl;
	lin_cleansed_jet = (*linear_cleanser_B_trimmed)( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(),
							 jets_tracks_PU[i].constituents() );
        cout << "SCZ With trimming pt eta phi m " << lin_cleansed_jet.pt() << " " << lin_cleansed_jet.eta() << " " << lin_cleansed_jet.phi() << " " << lin_cleansed_jet.m() << endl;
        lin_cleansed_jet = (*linear_cleanser_B_filtered)( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(),
                                                         jets_tracks_PU[i].constituents() );
        cout << "SCZ With filtering pt eta phi m " << lin_cleansed_jet.pt() << " " << lin_cleansed_jet.eta() << " " << lin_cleansed_jet.phi() << " " << lin_cleansed_jet.m() << endl;
      }
      */
    if (fSaveLinearTrimmed) {
      PseudoJet lin_cleansed_jet = (*linear_cleanser_B_trimmed)( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(),
						       jets_tracks_PU[i].constituents() );
      if (lin_cleansed_jet.pt() > fJetPTMin) {
        momentum.SetPxPyPzE(lin_cleansed_jet.px(), lin_cleansed_jet.py(), lin_cleansed_jet.pz(), lin_cleansed_jet.E());
        candidate = factory->NewCandidate();
        candidate->Momentum = momentum;
        candidate->Area.SetPxPyPzE(0.,0.,0.,0.);
        fLinearTrimmedOutputArray->Add(candidate);
      }
    }

    if (fSaveLinearFiltered) {
      PseudoJet lin_cleansed_jet = (*linear_cleanser_B_filtered)( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(),
                                                       jets_tracks_PU[i].constituents() );
      if (lin_cleansed_jet.pt() > fJetPTMin) {
        momentum.SetPxPyPzE(lin_cleansed_jet.px(), lin_cleansed_jet.py(), lin_cleansed_jet.pz(), lin_cleansed_jet.E());
        candidate = factory->NewCandidate();
        candidate->Momentum = momentum;
        candidate->Area.SetPxPyPzE(0.,0.,0.,0.);
        fLinearFilteredOutputArray->Add(candidate);
      }
    }

    if (fSaveGaussian) {
      PseudoJet gau_cleansed_jet = (*gaussian_cleanser_B)( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(),
							jets_tracks_PU[i].constituents() );
      if (gau_cleansed_jet.pt() > fJetPTMin) {
	momentum.SetPxPyPzE(gau_cleansed_jet.px(), gau_cleansed_jet.py(), gau_cleansed_jet.pz(), gau_cleansed_jet.E());
	candidate = factory->NewCandidate();
	candidate->Momentum = momentum;
	candidate->Area.SetPxPyPzE(0.,0.,0.,0.);
	fGaussianOutputArray->Add(candidate);
      }
    }
    //  delete sequence;
  }
  hard_event_charged.clear();
  pileup_charged.clear();
  full_event.clear();
  full_event_neutral.clear();
  jets_plain.clear();
  jets_neutrals.clear();
  jets_tracks_LV.clear();
  jets_tracks_PU.clear();
  //  for (int i = 0 ; i < jet_sets.size() ; i++) jet_sets[i].clear();
  //  jet_sets.clear();
}
