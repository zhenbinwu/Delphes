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
//   Loop over all events
//----------------------------------------------------------------------------
  int lept=0;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
      // Take first jet
      Jet *jet = (Jet*) branchJet->At(0);
      
      // Plot jet transverse momentum
      histJetPT->Fill(jet->PT);
      
      // Print jet transverse momentum
      //std::cout << jet->PT << std::endl;
    }

    if (entry % 500 == 0)
      std::cout << "--------------------" << entry << std::endl;
//----------------------------------------------------------------------------
//  Lepton Efficiency Exercise
//----------------------------------------------------------------------------
    std::map<int, int> MatchIdxe = MatchingLepton<Electron>(branchParticle, branchElectron, 11);
    std::map<int, int> MatchIdxm = MatchingLepton<Muon>(branchParticle, branchMuon, 13);
    if (MatchIdxe.size() + MatchIdxm.size() == 0) continue;
    lept++;
    std::map<int, int> MatchIske = MatchingLepton<Electron>(branchParticle, branchIsoTrk, 11);
    std::map<int, int> MatchIskm = MatchingLepton<Muon>(branchParticle, branchIsoTrk, 13);
    //for(std::map<int, int>::iterator it=MatchIdx.begin();
      //it!=MatchIdx.end(); it++)
    //{
      //GenParticle *gen = (GenParticle*) branchParticle->At(it->first);
      //histGenLepEta->Fill(gen->Eta);
      //if (it->second != -1)
      //{
        //histMatchGenLepEta->Fill(gen->Eta);
      //}
    //}


    if ((branchElectron->GetEntries() + branchMuon->GetEntries()) == 0) 
    {
      lepveto++;
      int NISK = 0;
      for (int k = 0; k < branchIsoTrk->GetEntries(); ++k)
      {
        Muon *trk = (Muon*) branchIsoTrk->At(k);
        if (trk->PT > 10) NISK++;
      }
      if (NISK == 0) isoveto++;

    }

  } // End of looping events

  std::cout <<"Total " <<  lept << "Lepton Veto " << lepveto << " with ISOTRK "  << isoveto<< std::endl;
 
  histLepEffEta = (TH1*)histMatchGenLepEta->Clone("histLepEffEta");
  histLepEffEta->SetTitle("Lepton Efficiency");
  histLepEffEta->Divide(histGenLepEta);

  // Saving resulting histograms
  histJetPT->Write();
  histGenLepEta->Write();
  histMatchGenLepEta->Write();
  histLepEffEta->Write();
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
    std::cout << "i " << i << " p " << p->PID << " status " << p->Status  << std::endl;
    //if (p->Status != 1 ) //Only select stable particle
      //continue;
    if ( (p->M1 != -1 && fabs(((GenParticle*)branchParticle->At(p->M1))->PID) != 24) && 
        (p->M2 != -1 && fabs(((GenParticle*)branchParticle->At(p->M2))->PID) != 24 ))
        continue;  //Making sure the lepton from W decay 
    if (std::fabs(p->PID) == PID) //Matched to the wanted lepton
    {
      std::cout << p->PID << std::endl;
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
      if (p->P4().DeltaR(lep->P4())<0.3)
      {
        it->second = i;
        break;
      }
    }
  }

  return MatchIdx;
}       // -----  end of function MatchingLepton  -----
