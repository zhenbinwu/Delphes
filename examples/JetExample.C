/*
root -l examples/Example1.C\(\"delphes_output.root\"\)

or just:

./JetExample delphes_output.root

N.B. you must touch JetExample.cpp in order to convince Make that there is anything to do to build this!
*/

//------------------------------------------------------------------------------

void JetExample(const char *inputFile)
{
  //  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchRawJetNoPU = treeReader->UseBranch("RawJetNoPU");
  TClonesArray *branchRawJet = treeReader->UseBranch("RawJet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchRho = treeReader->UseBranch("Rho");
  TClonesArray *branchNPU = treeReader->UseBranch("NPU");

  // Constituents will be 0 otherwise
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
  TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
 

  // Book histograms

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    cout << "Event " << entry << endl;
    
    for (int i = 0 ; i < branchRho->GetEntries() ; i++) {
      Rho *rho = (Rho*) branchRho->At(i);
      cout << "  Rho (" << rho->Edges[0] << "-" << rho->Edges[1] << "): " << rho->Rho << endl;
    }

    // I have cheated and recorder the true number of pileup vertices in a "ScalarHT" object!
    ScalarHT *NPU = (ScalarHT*) branchNPU->At(0);
    int nPUvertices_true = (int)NPU->HT;
    cout << "  Number of true pileup vertices: " << nPUvertices_true << endl;

    // Loop over jets
    for (int i = 0 ; i < branchJet->GetEntries() ; i++) {
      Jet *jet = (Jet*) branchJet->At(i);
      if (jet->PT > 100.) {
        cout << "  Jet " << i << endl;
        cout << "    pT: " << jet->PT << endl;
        cout << "    Eta: " << jet->Eta << endl;
        cout << "    Area: " << jet->AreaP4().Pt() << endl;
        cout << "    Time at calo: " << jet->t0 << endl;
	cout << "    Time at calo smeared by 1ps 10ps 20ps 30ps 40ps: " << jet->t1 << " " << jet->t10 << " " << jet->t20 << " " << jet->t30 << " " << jet->t40 << endl;
	cout << "    Beta*:" << jet->BetaStar << endl;
	cout << "    NCharged: " << jet->NCharged << endl;
	cout << "    Number of constituents: " << jet->Constituents.GetEntries() << endl;
	if (branchEFlowTrack && branchEFlowTower && branchEFlowMuon) {
	  for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
	    TObject *obj = jet->Constituents[j];
	    if (obj->IsA() == Tower::Class()) {
	      Tower *tow = static_cast<Tower *> ( obj ) ;
	      cout << "     Jet constituent Et Eta Phi Time (at calo) " << tow->ET << " " << tow->Eta << " " << tow->Phi << " " << tow->t0 << endl;
	    } else {
	      cout << "  not a tower - could check if it's a track instead (cf. Example3.C)" << endl;
	    }
	  }
	}
      }
    }

    // Loop over jets
    for (int i = 0 ; i < branchRawJet->GetEntries() ; i++) {
      Jet *jet = (Jet*) branchRawJet->At(i);
      if (jet->PT > 100.) {
	cout << "  Raw Jet " << i << endl;
	cout << "    pT: " << jet->PT << endl;
	cout << "    Eta: " << jet->Eta << endl;
	cout << "    Area: " << jet->AreaP4().Pt() << endl;
      }
    }
    for (int i = 0 ; i < branchRawJetNoPU->GetEntries() ; i++) {
      Jet *jet = (Jet*) branchRawJetNoPU->At(i);
      if (jet->PT > 30.) {
        cout << "  Raw Jet (no PU mixing) " << i << endl;
        cout << "    pT: " << jet->PT << endl;
        cout << "    Eta: " << jet->Eta << endl;
        cout << "    Area: " << jet->AreaP4().Pt() << endl;
      }
    }

    for (int i = 0 ; i < branchGenJet->GetEntries() ; i++) {
      Jet *jet = (Jet*) branchGenJet->At(i);
      if (jet->PT > 30.) {
	cout << "  Gen Jet " << i << endl;
	cout << "    pT: " << jet->PT << endl;
	cout << "    Eta: " << jet->Eta << endl;
        cout << "    Area: " << jet->AreaP4().Pt() << endl;

      }
    }


  }

}
