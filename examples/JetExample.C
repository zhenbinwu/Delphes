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
  TClonesArray *branchRho = treeReader->UseBranch("Rho");
  TClonesArray *branchNPU = treeReader->UseBranch("NPU");
  
  // Book histograms

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    cout << "Event " << entry << endl;
    ScalarHT *rho = (ScalarHT*) branchRho->At(0);
    cout << "  Rho: " << rho->HT << endl;

    // I have cheated and recorder the true number of pileup vertices in a "ScalarHT" object!
    ScalarHT *NPU = (ScalarHT*) branchNPU->At(0);
    int nPUvertices_true = (int)NPU->HT;
    cout << "  Number of true pileup vertices: " << nPUvertices_true << endl;

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
