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

  TClonesArray *branchGenJetWithPU = treeReader->UseBranch("GenJetWithPU");


  // Constituents will be 0 otherwise
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
  TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchBeamSpotParticle = treeReader->UseBranch("BeamSpotParticle");
  TClonesArray *branchGenParticleWithPU = treeReader->UseBranch("ParticleWithPU");

  bool verbose = true;
  bool doResTree = false;

  float geneta, genpt, recpt, dr;
  TFile *f;
  TTree *tree;

  if (doResTree) {
    f = new TFile("out.root","RECREATE");
    tree = new TTree("tree","tree");
    tree->Branch("geneta",&geneta,"geneta/F");
    tree->Branch("genpt",&genpt,"genpt/F");
    tree->Branch("recpt",&recpt,"recpt/F");
    tree->Branch("dr",&dr,"dr/F");
  }


  // Book histograms

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    if (verbose||entry%5000==0) cout << "Event " << entry << " / " << numberOfEntries << endl;
    
    for (int i = 0 ; i < branchRho->GetEntries() ; i++) {
      Rho *rho = (Rho*) branchRho->At(i);
      if (verbose) cout << "  Rho (" << rho->Edges[0] << "-" << rho->Edges[1] << "): " << rho->Rho << endl;
    }


    cout << "before scalarHT" << endl;
    // I have cheated and recorded the true number of pileup vertices in a "ScalarHT" object!
    ScalarHT *NPU = (ScalarHT*) branchNPU->At(0);
    int nPUvertices_true = (int)NPU->HT;
    if (verbose) cout << "  Number of true pileup vertices: " << nPUvertices_true << endl;

    if (doResTree) {
      // some very simple variables in a very simple tree
      for (int i = 0 ; i < branchGenJet->GetEntries() ; i++) {
	Jet *genjet = (Jet*) branchGenJet->At(i);
	TLorentzVector genp4 = genjet->P4();
	genpt = genjet->PT;
	geneta = genjet->Eta;
	dr = 9999.;
	recpt = -1.;
	for (int j = 0 ; j < branchJet->GetEntries() ; j++ ) {
	  Jet *jet = (Jet*) branchJet->At(j);
	  TLorentzVector p4 = jet->P4();
	  float tempdr = p4.DeltaR(genp4);
	  if (tempdr < dr) {
	    dr = tempdr;
	    recpt = jet->PT;
	  }
	}
	tree->Fill(); // per genjet
      }
    } // doResTree


    cout << "before branchBeamSpotParticle" << endl;
    // One particle from primary vertex
    if (verbose && branchBeamSpotParticle) {
      cout << branchBeamSpotParticle->GetEntries() << endl;
      GenParticle *part = (GenParticle*) branchBeamSpotParticle->At(0);
      cout <<  "  True primary vertex X Y Z T: " << part->X << " " << part->Y << " " << part->Z << " " << part->T << endl;
      //      cout << "     Generator particle PID Pt Eta Phi X Y Z T (at origin) "  << part->PID << " "
      //	   << part->PT << " " << part->Eta << " " << part->Phi << " " << part->X << " " << part->Y << " " << part->Z << " " << part->T << endl;
    }
    cout << "after branchBeamSpotParticle" << endl;


    // Status code 3 particle collection
    if (verbose && branchGenParticle) {
      for (int i = 0 ; i < branchGenParticle->GetEntries() ; i++ ) {
	GenParticle *part = (GenParticle*) branchGenParticle->At(i);
	cout << "     Status code 3 generator particle PID Pt Eta Phi Z T (at origin) "  << part->PID << " "
	     << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
      }
    }
      
    if (verbose) {
      // Loop over jets
      for (int i = 0 ; i < branchJet->GetEntries() ; i++) {
	//	cout << "  Jet " << i << "/" << branchJet->GetEntries() << endl;
	Jet *jet = (Jet*) branchJet->At(i);
	//	cout << jet->PT << endl;
	//	continue;
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
		cout << "     Jet constituent Et Eta Phi Time (at calo) " << tow->ET << " " << tow->Eta << " " << tow->Phi << " " << tow->TOuter << endl;
		if (jet->P4().DeltaR(tow->P4()) > 0.5) cout << "        !!!! dr=" << jet->P4().DeltaR(tow->P4()) << endl;
	      } else {
		cout << "  not a tower - could check if it's a track instead (cf. Example3.C)" << endl;
	      }
	    }
	  }
	}
      }
      
      // Loop over jets
      if (branchRawJet) {
	for (int i = 0 ; i < branchRawJet->GetEntries() ; i++) {
	  Jet *jet = (Jet*) branchRawJet->At(i);
	  if (jet->PT > 100.) {
	    cout << "  Raw Jet " << i << endl;
	    cout << "    pT: " << jet->PT << endl;
	    cout << "    Eta: " << jet->Eta << endl;
	    cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  }
	}
      }
      if (branchRawJetNoPU) {
	for (int i = 0 ; i < branchRawJetNoPU->GetEntries() ; i++) {
	  Jet *jet = (Jet*) branchRawJetNoPU->At(i);
	  if (jet->PT > 30.) {
	    cout << "  Raw Jet (no PU mixing) " << i << endl;
	    cout << "    pT: " << jet->PT << endl;
	    cout << "    Eta: " << jet->Eta << endl;
	    cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  }
	}
      }
	
      for (int i = 0 ; i < branchGenJet->GetEntries() ; i++) {
	Jet *jet = (Jet*) branchGenJet->At(i);
	if (jet->PT > 30.) {
	  cout << "  Gen Jet " << i << endl;
	  cout << "    pT: " << jet->PT << endl;
	  cout << "    Eta: " << jet->Eta << endl;
	  cout << "    Area: " << jet->AreaP4().Pt() << endl;
	  cout << "    Constituents: " << jet->Constituents.GetEntries() << endl;
	  if (branchGenParticleWithPU) { // now same colection
	    for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
	      TObject *obj = jet->Constituents[j];
	      if (obj->IsA() == GenParticle::Class()) {
		GenParticle *part = static_cast<GenParticle *> ( obj ) ;
		cout << "     Jet constituent Pt Eta Phi Z T (at origin) " << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
	      } else {
		cout << "     Jet constituent is not a particle (?)" << endl;
	      }
	    }
	  }
	}
      }

      if (branchGenJetWithPU) {
	for (int i = 0 ; i < branchGenJetWithPU->GetEntries() ; i++) {
	  Jet *jet = (Jet*) branchGenJetWithPU->At(i);
	  if (jet->PT > 100.) {
	    cout << "  Gen Jet (with PU included) " << i << endl;
	    cout << "    pT: " << jet->PT << endl;
	    cout << "    Eta: " << jet->Eta << endl;
	    cout << "    Area: " << jet->AreaP4().Pt() << endl;
	    cout << "    Constituents: " << jet->Constituents.GetEntries() << endl;
	    if (branchGenParticleWithPU) {
	      for (int j = 0 ; j < jet->Constituents.GetEntries() ; j++) {
		TObject *obj = jet->Constituents[j];
		if (obj->IsA() == GenParticle::Class()) {
		  GenParticle *part = static_cast<GenParticle *> ( obj ) ;
		  cout << "     Jet constituent Pt Eta Phi Z T (at origin) " << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
		} else {
		  cout << "     Jet constituent is not a particle (?)" << endl;
		}
	      }
	    }
	  }
	}
      }


    } // verbose 
  } // event

  if (doResTree) {
    f->cd();
    tree->Write();
    f->Close();
  }

}
