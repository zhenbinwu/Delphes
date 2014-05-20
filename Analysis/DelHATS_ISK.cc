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
  TClonesArray *branchIsoTrk   = treeReader->UseBranch("IsoTrack");
  
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
//  1st index: had, lep, 1lep, 2lep
//  2nd index: original, +lepton veto, +tau veto, +isotrack veto
//----------------------------------------------------------------------------
  double EventCount[4][4];
  for (int i = 0; i < 4; ++i)
  {
    for (int j = 0; j < 4; ++j)
    {
      EventCount[i][j] = 0;
    }
  }
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
//  Lepton Efficiency Exercise
//----------------------------------------------------------------------------
    std::map<int, int> MatchIdxe = MatchingLepton<Electron>(branchParticle, branchElectron, 11);
    std::map<int, int> MatchIdxm = MatchingLepton<Muon>(branchParticle, branchMuon, 13);
    std::map<int, int> MatchIske = MatchingLepton<Muon>(branchParticle, branchIsoTrk, 11);
    std::map<int, int> MatchIskm = MatchingLepton<Muon>(branchParticle, branchIsoTrk, 13);
    std::map<int, int> MatchIskt = MatchingLepton<Muon>(branchParticle, branchIsoTrk, 15);



//----------------------------------------------------------------------------
//  Filling in the number of events after vetoes
//----------------------------------------------------------------------------
    int lepcount = MatchIske.size() + MatchIskm.size() + MatchIskt.size();

    assert(lepcount < 3);
    //if (lepcount >= 3)

    // Original count
    EventCount[lepcount][0]++;
    if (lepcount > 0) EventCount[3][0]++;
    
    int vetocount[3] = {0, 0, 0};
    // after lepton veto
    if (branchElectron->GetEntries() + branchMuon->GetEntries() == 0 )
      vetocount[0]=1;
    
    if (vetocount[0])
    {
      bool hastau = false;
      for (int jit = 0; jit < branchJet->GetEntries(); ++jit)
      {
        Jet *jet = (Jet*) branchJet->At(jit);
        if (jet->TauTag) 
        {
          hastau = true;
          break;
        }
      }

      if (hastau == false) vetocount[1] = 1;
    }

    if (vetocount[1])
    {
      if (branchIsoTrk->GetEntries() == 0)
        vetocount[2] = 1;
    }

    // Fill in the matrix

    for (int i = 0; i < 3; ++i)
    {
      EventCount[lepcount][i+1] += vetocount[i];
      if (lepcount > 0) EventCount[3][i+1] += vetocount[i];
    }


//----------------------------------------------------------------------------
//  Plot lost lepton
//----------------------------------------------------------------------------

    for(std::map<int, int>::iterator it=MatchIdxe.begin();
      it!=MatchIdxe.end(); it++)
    {
      assert(MatchIske.find(it->first) != MatchIske.end());
      if (it->second == -1 && MatchIske[it->first] == -1)
      //if (it->second == -1)
      {
        GenParticle *gen = (GenParticle*) branchParticle->At(it->first);
        histLostEleEta->Fill(gen->Eta);
        histLostElePT->Fill(gen->PT);
      }
    }


    for(std::map<int, int>::iterator it=MatchIdxm.begin();
      it!=MatchIdxm.end(); it++)
    {
      assert(MatchIskm.find(it->first) != MatchIskm.end());
      if (it->second == -1 && MatchIskm[it->first] == -1)
      //if (it->second == -1 )
      {
        GenParticle *gen = (GenParticle*) branchParticle->At(it->first);
        histLostMuonEta->Fill(gen->Eta);
        histLostMuonPT->Fill(gen->PT);
      }
    }






    std::vector<int> EIsk;
    std::vector<int> MIsk;
    std::vector<int> TIsk;

    
    for(std::map<int, int>::iterator it=MatchIske.begin();
      it!=MatchIske.end(); it++)
    {
      EIsk.push_back(it->second);
      GenParticle *gen = (GenParticle*) branchParticle->At(it->first);
      histGenEleEta->Fill(gen->Eta);
      histGenElePt->Fill(gen->PT);
      if (it->second != -1)
      {
        histMatchGenEleEta->Fill(gen->Eta);
        histMatchGenElePt->Fill(gen->PT);
      }
    }

    for(std::map<int, int>::iterator it=MatchIskm.begin();
      it!=MatchIskm.end(); it++)
    {
      MIsk.push_back(it->second);
      GenParticle *gen = (GenParticle*) branchParticle->At(it->first);
      histGenMuonEta->Fill(gen->Eta);
      histGenMuonPt->Fill(gen->PT);
      if (it->second != -1)
      {
        histMatchGenMuonEta->Fill(gen->Eta);
        histMatchGenMuonPt->Fill(gen->PT);
      }
    }


    for(std::map<int, int>::iterator it=MatchIskt.begin();
      it!=MatchIskt.end(); it++)
    {
      TIsk.push_back(it->second);
      GenParticle *gen = (GenParticle*) branchParticle->At(it->first);
      histGenTauEta->Fill(gen->Eta);
      histGenTauPt->Fill(gen->PT);
      if (it->second != -1)
      {
        histMatchGenTauEta->Fill(gen->Eta);
        histMatchGenTauPt->Fill(gen->PT);
      }
    }



//----------------------------------------------------------------------------
//  Loop over track collection for matching
//----------------------------------------------------------------------------
    for (int i = 0; i < branchIsoTrk->GetEntries(); ++i)
    {
      int found = 0;
      if (std::find(EIsk.begin(), EIsk.end(), i) != EIsk.end()) found++;
      if (std::find(MIsk.begin(), MIsk.end(), i) != MIsk.end()) found++;
      if (std::find(TIsk.begin(), TIsk.end(), i) != TIsk.end()) found++;
      if (found > 1)
      {
        std::cout<<"Run to \033[0;31m"<<__func__<<"\033[0m at \033[1;36m"<< __FILE__<<"\033[0m, line \033[0;34m"<< __LINE__<<"\033[0m"<< std::endl; 
      }
      if (found == 0)
      {
        Muon *trk = (Muon*) branchIsoTrk->At(i);
        histMissIsKEta->Fill(trk->Eta);
        histMissIsKPT->Fill(trk->PT);
      }
    }

    //if (TIsk.size() != 0)
    //{
      //for (int i = 0; i < branchParticle->GetEntries(); ++i)
      //{
        //GenParticle *p = (GenParticle*) branchParticle->At(i);

        //std::cout << "i " << i << " p " << p->PID << " status " << p->Status
          //<< " m1 " << p->M1 << " m2 " << p->M2 <<" eta " << p->Eta << " PT " << p->PT<< std::endl;
      //}
    //}

  } // End of looping events



  std::cout << "[table border=\"1\"]" << std::endl;

  std::cout << "[b][center]Cut[/center][/b]            | [b][center]Hadronic[/center][/b] | [b][center]1 lepton[/center][/b] | [b][center]2 lepton[/center][/b] | [b][center]leptonic[/center][/b] |- " << std::endl;

  std::cout << "[b][center]Original[/center][/b]       | " << EventCount[0][0] << "       | " << EventCount[1][0] << "    | " << EventCount[2][0] << "    | " << EventCount[3][0] << "    |- "<< std::endl;
  std::cout << "[b][center]+Lepton veto[/center][/b]   | " << EventCount[0][1] << "       | " << EventCount[1][1] << "    | " << EventCount[2][1] << "    | " << EventCount[3][1] << "    |- "<< std::endl;
  std::cout << "[b][center]+Tau veto[/center][/b]      | " << EventCount[0][2] << "       | " << EventCount[1][2] << "    | " << EventCount[2][2] << "    | " << EventCount[3][2] << "    |- "<< std::endl;
  std::cout << "[b][center]+Isotrack veto[/center][/b] | " << EventCount[0][3] << "       | " << EventCount[1][3] << "    | " << EventCount[2][3] << "    | " << EventCount[3][3] << "    |- "<< std::endl;

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
    for(std::map<int, int>::iterator it=MatchIdx.begin();
      it!=MatchIdx.end(); it++)
    {
      GenParticle *p = (GenParticle*) branchParticle->At(it->first);
      if (p->P4().DeltaR(lep->P4())<0.5)
      {
        it->second = i;
        break;
      }
    }
  }

  return MatchIdx;
}       // -----  end of function MatchingLepton  -----
