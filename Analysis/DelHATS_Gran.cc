// ===========================================================================
// 
//       Filename:  DelHATS.cc
// 
//    Description:  An example code for Delphes exercise for lepton efficiency
//    comparison between PhaseI and PhaseII 0PU samples
// 
//        Version:  1.0
//        Created:  02/11/2014 01:44:39 PM
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu, John Stupak
//        Company:  HATS@LPC
// 
// ===========================================================================

// Classes from STL
#include <cstdlib>
#include <iostream>  
#include <string>
#include <vector>
#include <cassert>

// Classes from ROOT
#include "TH1.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TApplication.h"

// Classes from Delphes
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "classes/DelphesClasses.h"


  template <class T>
std::map<int, int> MatchingLepton(TClonesArray *branchParticle, TClonesArray *branchLep, int PID);


bool PassSelection(TClonesArray *branchJet, TClonesArray *branchMet);
double CalcHT(TClonesArray* branchJet);
int CalcNJets(TClonesArray* branchJet);
bool ElectronVeto(TClonesArray *branchElectron);
bool MuonVeto(TClonesArray *branchMuon);
bool TauVeto(TClonesArray* branchJet);
bool IsoTrackVeto(TClonesArray* branchIsoTrk);
int FindTauDecays(int Wtau, TClonesArray *branchParticle );
int EventCategory(TClonesArray* branchParticle);
bool IsoTrackVeto(const IsoTrack* isk);
std::vector<int> GetFinalHad(std::vector<int> VGenHad, TClonesArray *branchParticle);


std::map<int, int> MatchingMuon(TClonesArray *branchParticle, TClonesArray *branchMuon);
std::map<int, int> MatchingElectron(TClonesArray *branchParticle, TClonesArray *branchElectron);
// ===  FUNCTION  ============================================================
//         Name:  main
//  Description:  Wrap up of Example1.C
// ===========================================================================
int main ( int argc, char *argv[] )
{
  if(argc != 3)
  {
    std::cout << " Usage: DelHATS input_file output_file" << std::endl;
    std::cout << " input_file - input file in ROOT format ('Delphes' tree)," << std::endl;
    std::cout << " output_file - output file in ROOT format" << std::endl;
    return EXIT_FAILURE;
  }

  // Getting the input filename
  const std::string inputFile_name  = argv[1];
  const std::string outputFile_name = argv[2];

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile_name.c_str());

  // Create the output file
  TFile outputfile(outputFile_name.c_str(), "RECREATE");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet        = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet     = treeReader->UseBranch("GenJet");
  TClonesArray *branchElectron   = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon       = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton     = treeReader->UseBranch("Photon");
  TClonesArray *branchMet        = treeReader->UseBranch("MissingET");
  TClonesArray *branchParticle   = treeReader->UseBranch("Particle");
  TClonesArray *branchIsoTrk     = treeReader->UseBranch("IsoTrack");
  
  // Book histograms
//----------------------------------------------------------------------------
//  Example 
//---------------------------------------------------------------------------
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 500, 0.0, 1000);

  int lepveto = 0;
  int isoveto = 0;
//----------------------------------------------------------------------------
//  Lepton Efficiency exercise
//----------------------------------------------------------------------------
  TH1 *histGenLepEta      = new TH1F("histGenLepEta ", "GenLepton Eta", 120, -6, 6);
  TH1 *histMatchGenLepEta = new TH1F("histMatchGenLepEta", "Matched GenLepton Eta", 120, -6, 6);
  TH1 *histLepEffEta      = new TH1F("histLepEffEta", "Lepton Efficiency", 120, -6, 6);


  TH1 *histGenEleEta      = new TH1F("histGenEleEta ", "GenEletron Eta", 120, -6, 6);
  TH1 *histMatchGenEleEta = new TH1F("histMatchGenEleEta", "Matched GenEletron Eta", 120, -6, 6);
  TH1 *histEleISKEffEta      = new TH1F("histEleISKEffEta", "Electron Isotrak Efficiency", 120, -6, 6);
  TH1 *histGenElePt      = new TH1F("histGenElePt ", "GenEletron Pt", 100, 0, 200);
  TH1 *histMatchGenElePt = new TH1F("histMatchGenElePt", "Matched GenEletron Pt", 100, 0, 200);
  TH1 *histEleISKEffPt      = new TH1F("histEleISKEffPt", "Electron Isotrak Efficiency", 100, 0, 200);

  TH1 *histGenMuonEta      = new TH1F("histGenMuonEta ", "GenMuontron Eta", 120, -6, 6);
  TH1 *histMatchGenMuonEta = new TH1F("histMatchGenMuonEta", "Matched GenMuontron Eta", 120, -6, 6);
  TH1 *histMuonISKEffEta      = new TH1F("histMuonISKEffEta", "Muonctron Isotrak Efficiency", 120, -6, 6);
  TH1 *histGenMuonPt      = new TH1F("histGenMuonPt ", "GenMuontron Pt", 100, 0, 200);
  TH1 *histMatchGenMuonPt = new TH1F("histMatchGenMuonPt", "Matched GenMuontron Pt", 100, 0, 200);
  TH1 *histMuonISKEffPt      = new TH1F("histMuonISKEffPt", "Muonctron Isotrak Efficiency", 100, 0, 200);

  TH1 *histGenTauEta      = new TH1F("histGenTauEta ", "GenTautron Eta", 120, -6, 6);
  TH1 *histMatchGenTauEta = new TH1F("histMatchGenTauEta", "Matched GenTautron Eta", 120, -6, 6);
  TH1 *histTauISKEffEta      = new TH1F("histTauISKEffEta", "Tauctron Isotrak Efficiency", 120, -6, 6);
  TH1 *histGenTauPt      = new TH1F("histGenTauPt ", "GenTautron Pt", 100, 0, 200);
  TH1 *histMatchGenTauPt = new TH1F("histMatchGenTauPt", "Matched GenTautron Pt", 100, 0, 200);
  TH1 *histTauISKEffPt      = new TH1F("histTauISKEffPt", "Tauctron Isotrak Efficiency", 100, 0, 200);

  TH1 *histMissIsKEta      = new TH1F("histMissISKEta ", "MissISK Eta", 120, -6, 6);
  TH1 *histMissIsKPT      = new TH1F("histMissISKPT", "MissISK PT", 100, 0, 200);

  TH1 *histLostEleEta      = new TH1F("histLostEleEta ", "Lost Electron Eta", 120, -6, 6);
  TH1 *histLostElePT      = new TH1F("histLostElePT", "Lost Electron PT", 100, 0, 200);
  TH1 *histLostMuonEta      = new TH1F("histLostMuonEta ", "Lost Muon Eta", 120, -6, 6);
  TH1 *histLostMuonPT      = new TH1F("histLostMuonPT", "Lost Muon PT", 100, 0, 200);

  TH1 *histLostEleISKEta      = new TH1F("histLostEleISKEta ", "Lost Electron Eta", 120, -6, 6);
  TH1 *histLostEleISKPT      = new TH1F("histLostEleISKPT", "Lost Electron PT", 100, 0, 200);
  TH1 *histLostMuonISKEta      = new TH1F("histLostMuonISKEta ", "Lost Muon Eta", 120, -6, 6);
  TH1 *histLostMuonISKPT      = new TH1F("histLostMuonISKPT", "Lost Muon PT", 100, 0, 200);
//----------------------------------------------------------------------------
//  JEtMET execise
//----------------------------------------------------------------------------
  //TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 500, 0.0, 1000);
  //TH1 *histJetEta = new TH1F("jet_eta", "jet Eta", 120, -6, 6);

  ////  MET and MET resolution
  //TH1 *histMET = new TH1F("MET", "MET", 100, 0.0, 1000);
  //TH1 *histMET_X = new TH1F("MET_X", "MET_X", 300, -300.0, 300.0);
  //TH1 *histMET_Y = new TH1F("MET_Y", "MET_Y", 300, -300.0, 300.0);

  //// MHT and MHT resolution
  //TH1 *histMHT = new TH1F("MHT", "MHT", 100, 0.0, 1000);
  //TH1 *histMHT_X = new TH1F("MHT_X", "MHT_X", 300, -300.0, 300.0);
  //TH1 *histMHT_Y = new TH1F("MHT_Y", "MHT_Y", 300, -300.0, 300.0);

  //// Jet response: approximate by TProfile
  //TProfile *histJetResPT = new TProfile("histJetResPT ", 
      //"Jet response as a fucntion of GenJet Pt", 100, 0, 500);
  //TProfile *histJetResEta = new TProfile("histJetResEta ", 
      //"Jet response as a function of GenJet Eta", 100, -5, 5);

  //// Jet resolution for Eta (0, 2.5, 4, 5) and Pt (30, 50);
  //TH1 *histJetEta1 = new TH1F("histJetEta1 ", "Jet Resolution with eta (0, 2.5)", 200, 4, 4);
  //TH1 *histJetEta2 = new TH1F("histJetEta2 ", "Jet Resolution with eta (2.5, 4)", 200, 4, 4);
  //TH1 *histJetEta3 = new TH1F("histJetEta3 ", "Jet Resolution with eta (4, 5)", 200, 4, 4);


//----------------------------------------------------------------------------
//  1st index: had, electron, muon, leptonic tau, 1-prong tau, 3-prong tau
//  2nd index: original, +lepton veto, +tau veto, +isotrack veto
//----------------------------------------------------------------------------
  double EventCount[6][4];
  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      EventCount[i][j] = 0;
    }
  }

  int lepevent = 0;
//----------------------------------------------------------------------------
//   Loop over all events
//----------------------------------------------------------------------------
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    //if (entry > 500 ) break;
    if (entry % 500 == 0)
      std::cout << "--------------------" << entry << std::endl;

    //----------------------------------------------------------------------------
    //  Event selections
    //----------------------------------------------------------------------------

    //----------------------------------------------------------------------------
    //  Lepton Efficiency Exercise
    //----------------------------------------------------------------------------
    std::map<int, int> MatchIdxe = MatchingElectron(branchParticle, branchElectron);
    std::map<int, int> MatchIdxm = MatchingMuon(branchParticle, branchMuon);
    std::map<int, int> MatchIske = MatchingLepton<IsoTrack>(branchParticle, branchIsoTrk, 11);
    std::map<int, int> MatchIskm = MatchingLepton<IsoTrack>(branchParticle, branchIsoTrk, 13);
    std::map<int, int> MatchIskt = MatchingLepton<IsoTrack>(branchParticle, branchIsoTrk, 15);


    if (MatchIdxe.size() + MatchIdxm.size() != 0) lepevent++;

    //----------------------------------------------------------------------------
    //  Plot lost lepton
    //----------------------------------------------------------------------------

    for(std::map<int, int>::iterator it=MatchIdxe.begin();
        it!=MatchIdxe.end(); it++)
    {
      assert(MatchIske.find(it->first) != MatchIske.end());

      GenParticle *genall = (GenParticle*) branchParticle->At(it->first);
      histGenEleEta->Fill(genall->Eta);
      histGenElePt->Fill(genall->PT);

      if (it->second == -1)
      {
        GenParticle *gen = (GenParticle*) branchParticle->At(it->first);
        histLostEleEta->Fill(gen->Eta);
        histLostElePT->Fill(gen->PT);

        if (MatchIske[it->first] == -1)
        {
          histLostEleISKEta->Fill(gen->Eta);
          histLostEleISKPT->Fill(gen->PT);
        }
      } else {
        histMatchGenEleEta->Fill(genall->Eta);
        histMatchGenElePt->Fill(genall->PT);
      }

    }


    for(std::map<int, int>::iterator it=MatchIdxm.begin();
        it!=MatchIdxm.end(); it++)
    {
      assert(MatchIskm.find(it->first) != MatchIskm.end());

      GenParticle *genall = (GenParticle*) branchParticle->At(it->first);
      histGenMuonEta->Fill(genall->Eta);
      histGenMuonPt->Fill(genall->PT);

      if (it->second == -1 )
      {
        GenParticle *gen = (GenParticle*) branchParticle->At(it->first);
        histLostMuonEta->Fill(gen->Eta);
        histLostMuonPT->Fill(gen->PT);

        if (MatchIskm[it->first] == -1)
        {

          histLostMuonISKEta->Fill(gen->Eta);
          histLostMuonISKPT->Fill(gen->PT);
        }
      }else{
        histMatchGenMuonEta->Fill(genall->Eta);
        histMatchGenMuonPt->Fill(genall->PT);
      }
    }


    if (! PassSelection(branchJet, branchMet)) continue;

    //----------------------------------------------------------------------------
    //  Filling in the number of events after vetoes
    //----------------------------------------------------------------------------
    int lepcount = MatchIske.size() + MatchIskm.size() + MatchIskt.size();
    //assert(lepcount < 3);

    // Ignore dilepton events
    if (lepcount >= 2) continue;


    //int category =1;
    int category = EventCategory( branchParticle);
    //std::cout << "category : " << category << std::endl;
    if (category == -1) continue;
    // Only consider single lepton events?

    // Original count
    EventCount[category][0]++;

    // Lepton Veto
    if (ElectronVeto(branchElectron)) continue;
    if (MuonVeto(branchMuon)) continue;
    EventCount[category][1]++;

    // Tau Veto
    if (TauVeto(branchJet)) continue;
    EventCount[category][2]++;

    if (IsoTrackVeto(branchIsoTrk)) continue;
    EventCount[category][3]++;


  } // End of looping events


  std::cout << " total lep vevent : " << lepevent << std::endl;

  std::cout << "[table border=\"1\"]" << std::endl;

  std::cout << "[b][center]Cut[/center][/b]            | [b][center]Hadronic[/center][/b] | [b][center]Electron[/center][/b] "
    << " | [b][center]Muon[/center][/b] | [b][center]leptonic tau[/center][/b] "
    << " | [b][center]1-prong Tau[/center][/b] | [b][center]2-prong tau[/center][/b] "
    <<" |- " << std::endl;

  std::cout << "[b][center]Original[/center][/b]       | " << EventCount[0][0] << "       | " << EventCount[1][0] << "    | " << EventCount[2][0] << "    | " << EventCount[3][0] << "  | " << EventCount[4][0] << "  | " << EventCount[5][0] << "  |- "<< std::endl;
  std::cout << "[b][center]+Lepton veto[/center][/b]   | " << EventCount[0][1] << "       | " << EventCount[1][1] << "    | " << EventCount[2][1] << "    | " << EventCount[3][1] << "  | " << EventCount[4][1] << "  | " << EventCount[5][1] << "  |- "<< std::endl;
  std::cout << "[b][center]+Tau veto[/center][/b]      | " << EventCount[0][2] << "       | " << EventCount[1][2] << "    | " << EventCount[2][2] << "    | " << EventCount[3][2] << "  | " << EventCount[4][2] << "  | " << EventCount[5][2] << "  |- "<< std::endl;
  std::cout << "[b][center]+Isotrack veto[/center][/b] | " << EventCount[0][3] << "       | " << EventCount[1][3] << "    | " << EventCount[2][3] << "    | " << EventCount[3][3] << "  | " << EventCount[4][3] << "  | " << EventCount[5][3] << "  |- "<< std::endl;

  std::cout << "[/table]" << std::endl;

 
  std::cout << " Total number of Events: " << numberOfEntries  << std::endl;
  histLepEffEta = (TH1*)histMatchGenLepEta->Clone("histLepEffEta");
  histLepEffEta->SetTitle("Lepton Efficiency");
  histLepEffEta->Divide(histGenLepEta);

  histEleISKEffEta = (TH1*)histMatchGenEleEta->Clone("histEleISKEffEta");
  histEleISKEffEta->SetTitle("Eleton Efficiency");
  histEleISKEffEta->Divide(histGenEleEta);


  histMuonISKEffEta = (TH1*)histMatchGenMuonEta->Clone("histMuonISKEffEta");
  histMuonISKEffEta->SetTitle("Muonton Efficiency");
  histMuonISKEffEta->Divide(histGenMuonEta);

  histTauISKEffEta = (TH1*)histMatchGenTauEta->Clone("histTauISKEffEta");
  histTauISKEffEta->SetTitle("Tauton Efficiency");
  histTauISKEffEta->Divide(histGenTauEta);

  histEleISKEffPt = (TH1*)histMatchGenElePt->Clone("histEleISKEffPt");
  histEleISKEffPt->SetTitle("Eleton Efficiency");
  histEleISKEffPt->Divide(histGenElePt);


  histMuonISKEffPt = (TH1*)histMatchGenMuonPt->Clone("histMuonISKEffPt");
  histMuonISKEffPt->SetTitle("Muonton Efficiency");
  histMuonISKEffPt->Divide(histGenMuonPt);

  histTauISKEffPt = (TH1*)histMatchGenTauPt->Clone("histTauISKEffPt");
  histTauISKEffPt->SetTitle("Tauton Efficiency");
  histTauISKEffPt->Divide(histGenTauPt);


  // Saving resulting histograms
  histJetPT->Write();
  histGenLepEta->Write();
  histMatchGenLepEta->Write();
  histLepEffEta->Write();

  histGenEleEta->Write();
  histGenElePt->Write();
  histMatchGenEleEta->Write();
  histMatchGenElePt->Write();
  histEleISKEffEta->Write();
  histEleISKEffPt->Write();

  histGenMuonEta->Write();
  histGenMuonPt->Write();
  histMatchGenMuonEta->Write();
  histMatchGenMuonPt->Write();
  histMuonISKEffEta->Write();
  histMuonISKEffPt->Write();


  histGenTauEta->Write();
  histGenTauPt->Write();
  histMatchGenTauEta->Write();
  histMatchGenTauPt->Write();
  histTauISKEffEta->Write();
  histTauISKEffPt->Write();

  histMissIsKEta->Write();
  histMissIsKPT->Write();

  histLostEleEta->Write();
  histLostElePT->Write();
  histLostMuonEta->Write();
  histLostMuonPT->Write();

  histLostEleISKEta->Write();
  histLostEleISKPT->Write();
  histLostMuonISKEta->Write();
  histLostMuonISKPT->Write();

  //histJetEta->Write();
  //histMET->Write();
  //histMET_X->Write();
  //histMET_Y->Write();
  //histMHT->Write();
  //histMHT_X->Write();
  //histMHT_Y->Write();
  //histJetResPT->Write();
  //histJetResEta->Write();
  //histJetEta1->Write();
  //histJetEta2->Write();
  //histJetEta3->Write();

  outputfile.Close();
  return EXIT_SUCCESS;
}				// ----------  end of function main  ----------

// ===  FUNCTION  ============================================================
//         Name:  MatchingLepton
//  Description:  A template fucntion for matching reco lepton with the
//  GenParticles specified by the PID (ele 11, muon 13, tau 15). Return a
//  mapping of the index of GenParticle to reco lepton
// ===========================================================================
  template <class T>
std::map<int, int> MatchingLepton(TClonesArray *branchParticle, TClonesArray *branchLep, int PID)
{
  //Mapping the GenParticle index with Lepton index
  std::map<int, int> MatchIdx;

//----------------------------------------------------------------------------
//  Getting the Gen Lepton
//----------------------------------------------------------------------------
  int GenSize = branchParticle->GetEntries(); 
  for (int i = 0; i < GenSize; ++i)
  {
    GenParticle *p = (GenParticle*) branchParticle->At(i);
    //std::cout << "i " << i << " p " << p->PID << " status " << p->Status  << std::endl;
    if (p->Status != 3 ) //Only select stable particle
      continue;
    if ( (p->M1 != -1 && fabs(((GenParticle*)branchParticle->At(p->M1))->PID) != 24) && 
        (p->M2 != -1 && fabs(((GenParticle*)branchParticle->At(p->M2))->PID) != 24 ))
        continue;  //Making sure the lepton from W decay 
    if (std::fabs(p->PID) == PID) //Matched to the wanted lepton
    {
      //std::cout << p->PID << std::endl;
      MatchIdx[i] = -1;
    }
  }

//----------------------------------------------------------------------------
//  Getting the matched lepton, simply take deltaR < 0.3 as matched
//----------------------------------------------------------------------------
  for (int i = 0; i < branchLep->GetEntries(); ++i)
  {
    T *lep = (T*)branchLep->At(i);

    if (IsoTrackVeto(lep))
    {
      for(std::map<int, int>::iterator it=MatchIdx.begin();
          it!=MatchIdx.end(); it++)
      {
        GenParticle *p = (GenParticle*) branchParticle->At(it->first);
        if (p->P4().DeltaR(lep->P4())<0.3)
        {
          it->second = i;
          break;
        }
      }

    }
  }

  return MatchIdx;
}       // -----  end of function MatchingLepton  -----

// ===  FUNCTION  ============================================================
//         Name:  PassSelection
//  Description:  
// ===========================================================================
bool PassSelection(TClonesArray *branchJet, TClonesArray *branchMet)
  //TClonesArray *branchElectron, TClonesArray *branchMuon,  TClonesArray *branchIsoTrk)
{
  assert(branchMet->GetEntries() == 1);

  //MET > 200GeV
  MissingET* met = (MissingET*) branchMet->At(0);
  if (met->MET < 200) return false;

  //HT > 450
  if (CalcHT(branchJet) < 450) return false;

  
  //> 1 pfjet
  if (CalcNJets(branchJet) <= 1) return false;

  return true;
}       // -----  end of function PassSelection  -----

// ===  FUNCTION  ============================================================
//         Name:  ElectronVeto
//  Description:  
// ===========================================================================
bool ElectronVeto(TClonesArray *branchElectron)
{
  int elecount = 0;
  
  for (int i = 0; i < branchElectron->GetEntries(); ++i)
  {
    Electron *ele = (Electron*)branchElectron->At(i);
    //if (std::fabs(ele) < 2.4) //Excluding 1.442 < |eta|< 1.566 
    if (std::fabs(ele->Eta) < 1.442 || 1.566 < std::fabs(ele->Eta) < 2.4)
    {
      if (ele->PT > 10 && ele->IsolationVar < 0.15 && ele->EhadOverEem < 0.15)
      {
        elecount++;
      }
    }
  }

  return elecount > 0 ? true : false;
}       // -----  end of function ElectronVeto  -----


// ===  FUNCTION  ============================================================
//         Name:  MuonVeto
//  Description:  
// ===========================================================================
bool MuonVeto(TClonesArray *branchMuon)
{
  
  int muoncount = 0;
  
  for (int i = 0; i < branchMuon->GetEntries(); ++i)
  {
    Muon *muon = (Muon*)branchMuon->At(i);
    if (std::fabs(muon->Eta) < 2.4) //Excluding 1.442 < |eta|< 1.566 
    {
      if (muon->PT > 10 && muon->IsolationVar < 0.5 )
      {
        muoncount++;
      }
    }
  }


  return muoncount > 0 ? true : false;
}       // -----  end of function MuonVeto  -----


// ===  FUNCTION  ============================================================
//         Name:  TauVeto
//  Description:  
// ===========================================================================
bool TauVeto(TClonesArray* branchJet)
{
  int taucount = 0;
  
  for (int i = 0; i < branchJet->GetEntries(); ++i)
  {
    Jet *jet = (Jet*)branchJet->At(i);
    if (std::fabs(jet->Eta) < 2.3 && jet->TauTag) 
    {
      if (jet->PT > 20)
      {
        taucount++;
      }
    }
  }

  return taucount > 0 ? true : false;
}       // -----  end of function TauVeto  -----

// ===  FUNCTION  ============================================================
//         Name:  CalcHT
//  Description:  
// ===========================================================================
double CalcHT(TClonesArray* branchJet)
{
  double HT = 0.0;

  for (int i = 0; i < branchJet->GetEntries(); ++i)
  {
    Jet *jet = (Jet*)branchJet->At(i);
    if (std::fabs(jet->Eta) < 3.0 && jet->PT > 50) 
    {
      HT += jet->PT;
    }
  }

  return HT;
}       // -----  end of function CalcHT  -----

// ===  FUNCTION  ============================================================
//         Name:  CalcNJets
//  Description: i
// ===========================================================================
int CalcNJets(TClonesArray* branchJet)
{

  int njets = 0;

  for (int i = 0; i < branchJet->GetEntries(); ++i)
  {
    Jet *jet = (Jet*)branchJet->At(i);
    if (std::fabs(jet->Eta) < 2.4 && jet->PT > 100) 
    {
      njets++;
    }
  }

  return njets;
}       // -----  end of function CalcNJets  -----

// ===  FUNCTION  ============================================================
//         Name:  IsoTrackVeto
//  Description:  
// ===========================================================================
bool IsoTrackVeto(TClonesArray* branchIsoTrk)
{
  int iskcount = 0;
  
  for (int i = 0; i < branchIsoTrk->GetEntries(); ++i)
  {
    IsoTrack *isk = (IsoTrack*)branchIsoTrk->At(i);
    if(IsoTrackVeto(isk)) iskcount++;
  }

  return iskcount > 0 ? true : false;
  
}       // -----  end of function IsoTrackVeto  -----


// ===  FUNCTION  ============================================================
//         Name:  IsoTrackVeto
//  Description:  
// ===========================================================================
bool IsoTrackVeto(const IsoTrack* isk)
{
  if (isk->IsEMCand)
  {
    
    if (std::fabs(isk->Eta) < 2.5 && isk->PT > 5
        && isk->IsolationVar < 0.2) 
      return true;
  } else {
    if (std::fabs(isk->Eta) < 2.5 && isk->PT > 10
        && isk->IsolationVar < 0.1) 
      return true;
  }

  return false;
  
}       // -----  end of function IsoTrackVeto  -----


// ===  FUNCTION  ============================================================
//         Name:  EventCategory
//  Description:  
// ===========================================================================
int EventCategory(TClonesArray* branchParticle)
{
  
  int cat = -1;
  int lepcount = 0;

  for (int i = 0; i < branchParticle->GetEntries(); ++i)
  {
    GenParticle *p = (GenParticle*) branchParticle->At(i);


    if (p->Status != 3 ) //Only select stable particle
      continue;
    if ( (p->M1 != -1 && fabs(((GenParticle*)branchParticle->At(p->M1))->PID) != 24) && 
        (p->M2 != -1 && fabs(((GenParticle*)branchParticle->At(p->M2))->PID) != 24 ))
        continue;  //Making sure the lepton from W decay 
    if (std::fabs(p->PID) == 11) // Electron 
    {
      lepcount ++;
      cat = 1;
    }

    if (std::fabs(p->PID) == 13) // Muon
    {
      lepcount ++;
      cat = 2;
    }

    if (std::fabs(p->PID) == 15) // tau
    {
      lepcount ++;
      cat = FindTauDecays(i, branchParticle);
    }
  }

  if (lepcount == 0) cat = 0;

  return cat;
}       // -----  end of function EventCategory  -----


// ===  FUNCTION  ============================================================
//         Name:  FindTauDecays
//  Description:  
// ===========================================================================
int FindTauDecays(int Wtau, TClonesArray *branchParticle )
{

  // Tau will radiate photon
  // Looping until we reach the last tau
  std::map<int, int> TempDau;
  bool finaltau = false;

  while (!finaltau)
  {
    TempDau.clear();
    finaltau = true; //Assume finaltau is true
    for (int j = Wtau+1; j < branchParticle->GetEntries(); ++j)
    {
      GenParticle *p = (GenParticle*) branchParticle->At(j);
      if (p->M1 == Wtau || p->M2 == Wtau)
      {
        TempDau[j] = p->Status;
        if (std::fabs(p->PID) == 15)
        {
          Wtau= j;
          finaltau = false; // find a tau
          break;
        }
      }
    }
  }

  // Get the final tau, now look into its decay products
  TempDau.clear();
  for (int j = Wtau+1; j < branchParticle->GetEntries(); ++j)
  {
    GenParticle *p = (GenParticle*) branchParticle->At(j);
    if (p->M1 == Wtau || p->M2 == Wtau)
    {
      //std::cout << "M 1 " << p->M1 <<"  M2 " << p->M2 << std::endl;
      TempDau[j] = p->PID;
      //std::cout << "i " << j << " p " << p->PID << " status " << p->Status
      //<< " m1 " << p->M1 << " m2 " << p->M2 <<" eta " << p->Eta << " PT " << p->PT<< std::endl;
    }
  }

  //----------------------------------------------------------------------------
  //  Only interested in pion (211), rho(213),  a1(20213)
  //----------------------------------------------------------------------------
  std::vector<int> VHad;

  for(std::map<int, int>::iterator fit=TempDau.begin();
      fit!=TempDau.end(); fit++)
  {
    GenParticle *p = (GenParticle*) branchParticle->At(fit->first);

    if (std::fabs(p->PID) == 211  ||
        std::fabs(p->PID) == 213 || std::fabs(p->PID) == 20213 )
    {
      VHad.push_back(fit->first);
    }
    //BUG:: yes, I countting W as leptonic decay directlyy.. really not that
    //important for this study anyway
    if (std::fabs(p->PID) == 11 || std::fabs(p->PID) == 13 || std::fabs(p->PID) == 24)
      return 3;

  }

//----------------------------------------------------------------------------
//  For  hadronic , find out how many prongs?
//----------------------------------------------------------------------------
  std::vector<int> Nhads =  GetFinalHad(VHad, branchParticle);

  if (Nhads.size() == 1) return 4;
  if (Nhads.size() == 3) return 5;

  return -1;
}       // -----  end of function FindTauDecays  -----


// ===  FUNCTION  ============================================================
//         Name:  GetFinalHad
//  Description:  
// ===========================================================================
std::vector<int> GetFinalHad(std::vector<int> VGenHad, TClonesArray *branchParticle)
{

  std::vector<int> VFinalHad;
  std::vector<int> VNotFinalHad;

  for (int i = 0; i < VGenHad.size(); ++i)
  {
    GenParticle *p1 = (GenParticle*) branchParticle->At(VGenHad.at(i));
    for (int j = VGenHad.at(i); j < branchParticle->GetEntries(); ++j)
    {
      GenParticle *p = (GenParticle*) branchParticle->At(j);
      if (p->M1 == VGenHad.at(i) || p->M2 == VGenHad.at(i))
      {
        if (p->Status == 1 && p->PID != 22)
          VFinalHad.push_back(p->PID);
        else
          VNotFinalHad.push_back(j);
      }
    }
  }

  // Now, loop until all finals
  while (VNotFinalHad.size() != 0)
  {
    std::vector<int> VTempNotFinalHad = VNotFinalHad;
    VNotFinalHad.clear();

    for (int i = 0; i < VTempNotFinalHad.size(); ++i)
    {
      for (int j = VTempNotFinalHad.at(i); j < branchParticle->GetEntries(); ++j)
      {
        GenParticle *p = (GenParticle*) branchParticle->At(j);
        if (p->M1 == VTempNotFinalHad.at(i) || p->M2 == VTempNotFinalHad.at(i))
        {
          if (p->Status == 1 && p->PID != 22)
            VFinalHad.push_back(p->PID);
          else
            VNotFinalHad.push_back(j);
        }
      }
    }

  }

//----------------------------------------------------------------------------
//  Print out the prongs of tau decays
//----------------------------------------------------------------------------
/*
 *  std::cout << " Final tau: ";
 *  for (int i = 0; i < VFinalHad.size(); ++i)
 *  {
 *    
 *    std::cout << " " << VFinalHad.at(i);
 *  }
 *
 *  std::cout  << std::endl;
 */


  return  VFinalHad;
}       // -----  end of function GetFinalHad  -----


std::map<int, int> MatchingElectron(TClonesArray *branchParticle, TClonesArray *branchElectron)
{
  //Mapping the GenParticle index with Lepton index
  std::map<int, int> MatchIdx;

//----------------------------------------------------------------------------
//  Getting the Gen Lepton
//----------------------------------------------------------------------------
  int GenSize = branchParticle->GetEntries(); 
  for (int i = 0; i < GenSize; ++i)
  {
    GenParticle *p = (GenParticle*) branchParticle->At(i);
    //std::cout << "i " << i << " p " << p->PID << " status " << p->Status  << std::endl;
    if (p->Status != 3 ) //Only select stable particle
      continue;
    if ( (p->M1 != -1 && fabs(((GenParticle*)branchParticle->At(p->M1))->PID) != 24) && 
        (p->M2 != -1 && fabs(((GenParticle*)branchParticle->At(p->M2))->PID) != 24 ))
        continue;  //Making sure the lepton from W decay 
    if (std::fabs(p->PID) == 11) //Matched to the wanted lepton
    {
      //std::cout << p->PID << std::endl;
      MatchIdx[i] = -1;
    }
  }

//----------------------------------------------------------------------------
//  Getting the matched lepton, simply take deltaR < 0.3 as matched
//----------------------------------------------------------------------------
  for (int i = 0; i < branchElectron->GetEntries(); ++i)
  {
    Electron *ele = (Electron*)branchElectron->At(i);
    if (std::fabs(ele->Eta) < 1.442 || 1.566 < std::fabs(ele->Eta) < 2.4)
    {
      if (ele->PT > 10 && ele->IsolationVar < 0.15 && ele->EhadOverEem < 0.15)
      {
        for(std::map<int, int>::iterator it=MatchIdx.begin();
            it!=MatchIdx.end(); it++)
        {
          GenParticle *p = (GenParticle*) branchParticle->At(it->first);
          if (p->P4().DeltaR(ele->P4())<0.3)
          {
            it->second = i;
            break;
          }
        }
      }
    }
  }

  return MatchIdx;
}       // -----  end of function MatchingElectron  -----

std::map<int, int> MatchingMuon(TClonesArray *branchParticle, TClonesArray *branchMuon)
{
  //Mapping the GenParticle index with Lepton index
  std::map<int, int> MatchIdx;

//----------------------------------------------------------------------------
//  Getting the Gen Lepton
//----------------------------------------------------------------------------
  int GenSize = branchParticle->GetEntries(); 
  for (int i = 0; i < GenSize; ++i)
  {
    GenParticle *p = (GenParticle*) branchParticle->At(i);
    //std::cout << "i " << i << " p " << p->PID << " status " << p->Status  << std::endl;
    if (p->Status != 3 ) //Only select stable particle
      continue;
    if ( (p->M1 != -1 && fabs(((GenParticle*)branchParticle->At(p->M1))->PID) != 24) && 
        (p->M2 != -1 && fabs(((GenParticle*)branchParticle->At(p->M2))->PID) != 24 ))
        continue;  //Making sure the lepton from W decay 
    if (std::fabs(p->PID) == 13) //Matched to the wanted lepton
    {
      //std::cout << p->PID << std::endl;
      MatchIdx[i] = -1;
    }
  }

//----------------------------------------------------------------------------
//  Getting the matched lepton, simply take deltaR < 0.3 as matched
//----------------------------------------------------------------------------
  for (int i = 0; i < branchMuon->GetEntries(); ++i)
  {
    Muon *muon = (Muon*)branchMuon->At(i);

    if (std::fabs(muon->Eta) < 2.4) //Excluding 1.442 < |eta|< 1.566 
    {
      if (muon->PT > 10 && muon->IsolationVar < 0.5 )
      {

        for(std::map<int, int>::iterator it=MatchIdx.begin();
            it!=MatchIdx.end(); it++)
        {
          GenParticle *p = (GenParticle*) branchParticle->At(it->first);
          if (p->P4().DeltaR(muon->P4())<0.3)
          {
            it->second = i;
            break;
          }
        }
      }
    }
  }

  return MatchIdx;
}       // -----  end of function MatchingMuon  -----
