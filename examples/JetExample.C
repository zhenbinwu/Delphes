/*
root -l examples/Example1.C\(\"delphes_output.root\"\)

or

./JetExample delphes_output.root
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
  TClonesArray *branchRawJet = treeReader->UseBranch("RawJet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchRho = treeReader->UseBranch("Rho");
  
  // Book histograms

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    cout << "Event " << entry << endl;
    ScalarHT *rho = (ScalarHT*) branchRho->At(0);
    cout << "  Rho: " << rho->HT << endl;

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

    for (int i = 0 ; i < branchGenJet->GetEntries() ; i++) {
      Jet *jet = (Jet*) branchGenJet->At(i);
      if (jet->PT > 50.) {
	cout << "  Gen Jet " << i << endl;
	cout << "    pT: " << jet->PT << endl;
	cout << "    Eta: " << jet->Eta << endl;
      }
    }


  }

}
