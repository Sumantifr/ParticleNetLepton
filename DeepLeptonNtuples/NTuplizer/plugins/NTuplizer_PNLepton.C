// user include files
#include "NTuplizer_PNLepton.h"
 
//
// member functions
//

// ------------ method called for each event  ------------
void
PNLepton::analyze(const edm::Event& iEvent, const edm::EventSetup& pset) {
  
  using namespace edm;
  Nevt++;
  
  irun = iEvent.id().run();
  ilumi = iEvent.luminosityBlock();
  ievt = iEvent.id().event();
  
  if (Nevt%100==1)cout <<"PNLepton::analyze "<<Nevt<<" "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
    
  // First store all MC information //
    
  if(isMC){
	
	// MC weights //
	  
    edm::Handle<GenEventInfoProduct>eventinfo ;  
    iEvent.getByToken(tok_wt_,eventinfo) ;
        
    if (eventinfo.isValid()){	Generator_weight = eventinfo->weight();  }
    else {  Generator_weight =  -10000;  }
   
    // Gen particles //
    
    edm::Handle<std::vector<reco::GenParticle>> genparticles;  
	iEvent.getByToken(tok_genparticles_,genparticles);
    
    GenPartons.clear();
    
	if(genparticles.isValid()){
	
		for(unsigned ig=0; ig<(genparticles->size()); ig++){
			
			if(!(((*genparticles)[ig].status()==1)||(abs((*genparticles)[ig].status())==22)||((*genparticles)[ig].status()==23))) continue;
				  
			//if(!((abs((*genparticles)[ig].pdgId())>=1 && abs((*genparticles)[ig].pdgId())<=6) || (abs((*genparticles)[ig].pdgId())>=11 && abs((*genparticles)[ig].pdgId())<=16) || (abs((*genparticles)[ig].pdgId())>=21 && abs((*genparticles)[ig].pdgId())<=25))) continue;
			if(!(abs((*genparticles)[ig].pdgId())==11 || abs((*genparticles)[ig].pdgId())==13 || abs((*genparticles)[ig].pdgId())==15) ) continue;
			// important condition on pdg id -> May be changed in other analyses //
	  
			GenParton parton; 
	  
			parton.pt = (*genparticles)[ig].pt();
			parton.eta = (*genparticles)[ig].eta();
			parton.phi = (*genparticles)[ig].phi();
			parton.mass = (*genparticles)[ig].mass();
			parton.p4.SetPtEtaPhiM(parton.pt,parton.eta,parton.phi,parton.mass);
    
			parton.status = (*genparticles)[ig].status();
			parton.pdgId = (*genparticles)[ig].pdgId();
	
			parton.fromhard = (*genparticles)[ig].isHardProcess();
			parton.fromhardbFSR = (*genparticles)[ig].fromHardProcessBeforeFSR();
			parton.isPromptFinalState = (*genparticles)[ig].isPromptFinalState();
			parton.isLastCopyBeforeFSR = (*genparticles)[ig].isLastCopyBeforeFSR();
			parton.isDirectPromptTauDecayProductFinalState = (*genparticles)[ig].isDirectPromptTauDecayProductFinalState();
		
			int mompdg, momstatus, grmompdg;
	  
			// mother pdg id & status //
			const Candidate * mom = (*genparticles)[ig].mother();
			mompdg = mom->pdgId();
			momstatus = mom->status();
			const Candidate * momtmp = (*genparticles)[ig].mother();
			while(mompdg == parton.pdgId)
			{
				mompdg = momtmp->mother()->pdgId();
				momstatus = momtmp->mother()->status();
				momtmp = momtmp->mother();
			}
	  
			// grand-mother pdg id //
			bool found_grmom = false;
			if(mom->numberOfMothers()>0){
				const Candidate * grmom  = mom->mother();
				for(int iter=0; iter<10; iter++){
					if(grmom->pdgId() != mom->pdgId()){
						grmompdg  = grmom->pdgId();
						found_grmom = true;
						break;
					}else{
						if(grmom->numberOfMothers()>0){
							grmom = grmom->mother();
						}
						else{ break; }
					}
				}
			}
			if(!found_grmom){
				grmompdg  = -10000000;
			}
			
			parton.mompdgId = mompdg;
			parton.momstatus = momstatus;
			parton.grmompdgId = grmompdg; 
						
			GenPartons.push_back(parton);
    
			if(int(GenPartons.size())>=npartmx) break;
	
		}
    }
    
    // Pileup information //
    
    Pileup_nPU = 0;
	Pileup_nTrueInt = 0;
    
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(pileup_, PupInfo);
    if (PupInfo.isValid()) {
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
		if (PVI->getBunchCrossing()==0) {
			Pileup_nPU = PVI->getPU_NumInteractions();
			Pileup_nTrueInt = PVI->getTrueNumInteractions();
			break;
			}
		}
    }

  }//isMC
  
  // Primary vertex info //
  
  Handle<VertexCollection> primaryVertices;
  iEvent.getByToken(tok_primaryVertices_, primaryVertices);
  reco::Vertex PV_vertex;
  
  if (primaryVertices.isValid()) {
	  
	if(primaryVertices->size() > 0){  
		PV_vertex = primaryVertices->at(0); 
	} 

    int nprimi_org = 0;
    
    for (reco::VertexCollection::const_iterator vert=primaryVertices->begin(); vert<primaryVertices->end(); vert++) {
      if (vert->isValid() && !vert->isFake()) {
		if (vert->ndof() > 4 && fabs(vert->position().z()) <= 24 && fabs(vert->position().Rho()) <= 2) {
			nprimi_org++;
			}
		}
    }
   
    PV_npvsGood = nprimi_org;
    
  } else { 
    PV_npvsGood = -100;
  }
  
  // Energy density info //

  edm::Handle<double> Rho_PF;
  iEvent.getByToken(tok_Rho_,Rho_PF);
  Rho = *Rho_PF;

  
  // ====== RECO-objects now  ==========//
  
  edm::Handle<edm::View<pat::PackedCandidate>> pfcands;                                                                                                          
  iEvent.getByToken(tok_pfcands_, pfcands);      
  
  edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVertices;
  iEvent.getByToken(tok_sv,secondaryVertices);
                                                                                                         
  
  edm::Handle<edm::View<pat::Muon>> muons;                                                                                                          
  iEvent.getByToken(tok_muons_, muons);       
  
  edm::Handle<edm::View<pat::Jet>> pfjetAK4s;
  iEvent.getByToken(tok_pfjetAK4s_, pfjetAK4s);                                                                                   
  
  std::vector<TLorentzVector> tlep;
  
  // Muons //
    
  nLepton = 0; 
  
  vector<Lepton> leptons;      
                                                                                                                                   
  if(muons.isValid() && muons->size()>0) {                                                                                                           
    
	edm::View<pat::Muon>::const_iterator muon1;                                                                                                      

    for( muon1 = muons->begin(); muon1 < muons->end(); muon1++ ) {                                                                                   

		if (StoreMuon(*muon1,minmuPt,maxEta)) {
			
			Lepton lepton;  
			Initialize(lepton);                                                              
			
			TrackRef trktrk = muon1->innerTrack();   
			
			lepton.pt = muon1->pt();                                                                                                                                                                             
			lepton.p = trktrk->p();                                                                                                          
			lepton.eta = muon1->eta();                                                                                                              
			lepton.phi = muon1->phi();  
			lepton.mass = muon1->mass();  
			lepton.charge = muon1->charge();                                                                                                                                                                                                                          
            lepton.pdgId = muon1->pdgId();
            lepton.tightcharge = (muon1->muonBestTrack()->ptError()/muon1->muonBestTrack()->pt() < 0.2)?2:0;
                                                                                                                                                       
			//MiniIsolation //     
			                                                                                 
			vector<float> isovalues;
			Read_MiniIsolation(muon1,Rho,isovalues);
			lepton.minisoall = isovalues[0];
			lepton.minchiso = isovalues[1];
			lepton.minnhiso = isovalues[2];
			lepton.minphiso = isovalues[3];
			                                         
			// Basic id variables //    
			                                  
			lepton.isPFCand = muon1->isPFMuon();                                                                                                        
			lepton.isGlobal = muon1->isGlobalMuon();                                                                                                    
			lepton.isTracker = muon1->isTrackerMuon();                                                                                                  
			lepton.isLoose = (muon::isLooseMuon(*muon1));                                                                                           
			lepton.isMedium = (muon::isMediumMuon(*muon1));                                                                                            
			lepton.isMedPr = false;                                                                          
			if(muon::isMediumMuon(*muon1)) {                                                                                                             
				if ((std::abs(muon1->muonBestTrack()->dz(PV_vertex.position())) < 0.1) && (std::abs(muon1->muonBestTrack()->dxy(PV_vertex.position())) < 0.02)){                                                                                                                  
					lepton.isMedPr = true;                                                                                                              
				}                                                                                                                                          
			}                                                                                                                                      
			lepton.isGoodGlobal = (muon1->isGlobalMuon() && muon1->globalTrack()->normalizedChi2() < 3 && muon1->combinedQuality().chi2LocalPosition < 12 && muon1->combinedQuality().trkKink < 20 && (muon::segmentCompatibility(*muon1)) > 0.303);                     
			lepton.isTight = (muon::isTightMuon(*muon1,PV_vertex));                                                                                    
			lepton.isHighPt = (muon::isHighPtMuon(*muon1,PV_vertex));                                                                                  
			lepton.isHighPttrk = (muon::isTrackerHighPtMuon(*muon1,PV_vertex));   
			
			// Displacement //
			
			lepton.dxy = muon1->muonBestTrack()->dxy(PV_vertex.position());                                                                         
			lepton.dz = muon1->muonBestTrack()->dz(PV_vertex.position());  
			lepton.dxyError = muon1->edB(pat::Muon::PV2D);   
			lepton.ip3d =  muon1->dB(pat::Muon::PV3D);    
			lepton.sip3d =  muon1->dB(pat::Muon::PV3D)/muon1->edB(pat::Muon::PV3D);     
			
			// Displacement w.r.t secondary vertex //
			 
			lepton.dxy_sv = DistanceFromSV_Muon(muon1->muonBestTrack(),secondaryVertices);  
			
			// GEN particle matching //
   
			if(isMC){
				lepton.genPartFlav = getGenPartonFlavor(lepton.eta,lepton.phi,GenPartons,13);
			}else{
				lepton.genPartFlav = -100;
			}
    
			// Energy info //
			                                                   
			lepton.e_ECAL = (muon1->calEnergy()).em * 1./muon1->energy();                                                                                                  
			lepton.e_HCAL = (muon1->calEnergy()).had* 1./muon1->energy();
			
			// Track info //
			                                                                         
			lepton.posmatch = muon1->combinedQuality().chi2LocalPosition;                                                                           
			lepton.trkKink = muon1->combinedQuality().trkKink;                                                                                       
			lepton.segmentComp = muon::segmentCompatibility(*muon1);                                                                                     
			                                                                                                                                                                                                            
			TrackRef trkglb =muon1->globalTrack();                                                                                                       
			if ((!muon1->isGlobalMuon())) {                                                                                                              
				if (muon1->isTrackerMuon()) {                                                                                                              
					trkglb =muon1->innerTrack();                                                                                                             
				} else {                                                                                                                                   
					trkglb =muon1->outerTrack();                                                                                                             
				}                                                                                                                                          
			}
			                                                                                                                                            
			lepton.chi2 = trkglb->normalizedChi2();                                                                                                  
			lepton.ndof = (int)trkglb->ndof();                                                                                                       
			lepton.hit = (int)trkglb->hitPattern().numberOfValidMuonHits();                                                                              
			lepton.nStations = (int)muon1->numberOfMatchedStations();                                                                                          
			lepton.pixhit = (int)trktrk->hitPattern().numberOfValidPixelHits();                                                                          
			lepton.nTrackerLayers = (int)trktrk->hitPattern().trackerLayersWithMeasurement();                                                                    
			lepton.valfrac = trktrk->validFraction();   
			lepton.lostHits = (int)trktrk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS); 
			
			// Isolation variables //
			                                                    
			lepton.pfRelIso04_drcor = (muon1->pfIsolationR04().sumChargedHadronPt + max(0., muon1->pfIsolationR04().sumNeutralHadronEt + muon1->pfIsolationR04().sumPhotonEt - 0.5*muon1->pfIsolationR04().sumPUPt))/muon1->pt();                                               
			lepton.pfRelIso04_ChargedHadron = muon1->pfIsolationR04().sumChargedHadronPt*1./muon1->pt();     
			lepton.pfRelIso04_NeutralHadron = muon1->pfIsolationR04().sumNeutralHadronEt*1./muon1->pt();     
			lepton.pfRelIso04_Photon = muon1->pfIsolationR04().sumPhotonEt*1./muon1->pt();     
			lepton.pfRelIso04_PileUp = muon1->pfIsolationR04().sumPUPt*1./muon1->pt();   
			
			lepton.pfRelIso03_drcor = (muon1->pfIsolationR03().sumChargedHadronPt + max(0., muon1->pfIsolationR03().sumNeutralHadronEt + muon1->pfIsolationR03().sumPhotonEt - 0.5*muon1->pfIsolationR03().sumPUPt))/muon1->pt();                                               
			lepton.pfRelIso03_ChargedHadron = muon1->pfIsolationR03().sumChargedHadronPt*1./muon1->pt();     
			lepton.pfRelIso03_NeutralHadron = muon1->pfIsolationR03().sumNeutralHadronEt*1./muon1->pt();     
			lepton.pfRelIso03_Photon = muon1->pfIsolationR03().sumPhotonEt*1./muon1->pt();     
			lepton.pfRelIso03_PileUp = muon1->pfIsolationR03().sumPUPt*1./muon1->pt();  
			lepton.tkRelIso = muon1->isolationR03().sumPt/muon1->tunePMuonBestTrack()->pt();
			
			// store 4-vector to vector of leptons //
			TLorentzVector p4;
			p4.SetPtEtaPhiM(lepton.pt, lepton.eta, lepton.phi, lepton.mass);
			tlep.push_back(p4);
			
			// now push lepton to vector of leptons //
			leptons.push_back(lepton);
			
			if (++nLepton>=njetmx) break;                                                                                                                 
		
		}                                                                                                                                              
      }                                                                                                                                               
  }// muon loop 
  
  // Electrons //
    
  for(const auto& electron1 : iEvent.get(tok_electrons_) ) {                                                                                          
 
    if (!StoreElectron(electron1,minePt,maxEta)) continue;
                    
    Lepton lepton;
    Initialize(lepton);                     
                                                   
	GsfTrackRef gsftrk1 = electron1.gsfTrack();   																														
    TrackRef ctftrk = electron1.closestCtfTrackRef();    
    
    // Basic kinematic variables //
    
    lepton.pt = electron1.pt();   	                                                                                 
    lepton.eta = electron1.eta();                                                                                                                 
    lepton.phi = electron1.phi();                                                                                                                 
    lepton.mass = electron1.mass();                                                                                              
    lepton.p = electron1.trackMomentumAtVtx().R();
    lepton.charge = electron1.charge();     
    lepton.pdgId = electron1.pdgId(); 
    lepton.tightcharge = electron1.isGsfCtfScPixChargeConsistent() + electron1.isGsfScPixChargeConsistent();
    
    // Energy variables //
    
    lepton.e_ECAL = electron1.ecalEnergy()*1./electron1.energy();   
    lepton.e_HCAL = (1.-electron1.ecalEnergy())*1./electron1.energy();  
    
    // MVA id //
    
    lepton.mvaid_Fallv2WP90 = electron1.electronID(melectronID_isowp90);                                                                                 
    lepton.mvaid_Fallv2WP90_noIso = electron1.electronID(melectronID_noisowp90);                                                                             
    lepton.mvaid_Fallv2WP80 = electron1.electronID(melectronID_isowp80);                                                                                 
    lepton.mvaid_Fallv2WP80_noIso = electron1.electronID(melectronID_noisowp80);   
    
    // displacement //
                                                                                 
    lepton.dxy = gsftrk1->dxy(PV_vertex.position());  
    lepton.dxyError = electron1.edB(pat::Electron::PV2D);                                                                                           
    lepton.dz = gsftrk1->dz(PV_vertex.position()); 
    lepton.dzError = electron1.edB(pat::Electron::PVDZ);    
    lepton.ip3d =  electron1.dB(pat::Electron::PV3D); 
    lepton.sip3d =  electron1.dB(pat::Electron::PV3D)/electron1.edB(pat::Electron::PV3D);    
    
    // Displacement w.r.t secondary vertex //
                                                                                                                                                   
    lepton.dxy_sv = DistanceFromSV_Electron(gsftrk1,secondaryVertices);                                                                                          
    
    // GEN particle matching //
    
    if(isMC){
		lepton.genPartFlav = getGenPartonFlavor(lepton.eta,lepton.phi,GenPartons,11);
	}else{
		lepton.genPartFlav = -100;
	}
    
    // supercluste info //
    
    lepton.supcl_energy = electron1.superCluster()->energy();  
    lepton.supcl_eta = electron1.superCluster()->eta();                                                                                           
    lepton.supcl_phi = electron1.superCluster()->phi();                                                                                           
    
    // shape of energy deposition //             
                                                                                                                                                                                                                                              
	lepton.sigmaietaieta = electron1.full5x5_sigmaIetaIeta();                                                                                         
    lepton.sigmaiphiiphi = electron1.full5x5_sigmaIphiIphi();  
    lepton.hcaloverecal = electron1.full5x5_hcalOverEcal();                                                                                         
    lepton.r9full = electron1.full5x5_r9(); 
    lepton.e1x5bye5x5 = 1.-electron1.full5x5_e1x5()/electron1.full5x5_e5x5();   
                                                                                                                                                                                                                                                  
    lepton.eoverp = 1./electron1.superCluster()->energy();                                                                                             
    lepton.hoe = electron1.hadronicOverEm();                                                                                                   
    lepton.ecloverpout = electron1.eEleClusterOverPout();  
    lepton.convVeto = electron1.passConversionVeto();
    
    lepton.eInvMinusPInv =  (1-electron1.eSuperClusterOverP())/electron1.ecalEnergy(); // OR 1.0/(electron1.ecalEnergy())-1.0/(electron1.trackMomentumAtVtx().R())
	lepton.etain = electron1.deltaEtaSuperClusterTrackAtVtx();                                                                                    
	lepton.phiin = electron1.deltaPhiSuperClusterTrackAtVtx();                                                                                     
	lepton.supcl_preshvsrawe = electron1.superCluster()->preshowerEnergy()/electron1.superCluster()->rawEnergy();                                                                                                                                            
	lepton.supcl_etaWidth = electron1.superCluster()->etaWidth();                                                                                     
	lepton.supcl_phiWidth = electron1.superCluster()->phiWidth(); 
	lepton.deltaetacltrkcalo = electron1.deltaEtaSeedClusterTrackAtCalo();   
	                                                                                                        
	lepton.convtxprob = electron1.convVtxFitProb();   
	lepton.fbrem = electron1.fbrem();                         
	                                                                    
    // isolation variables //                              
                                                                                                                                                                                       
    const reco::GsfElectron::PflowIsolationVariables& pfIso = electron1.pfIsolationVariables();  
    lepton.pfRelIso03_ChargedHadron = pfIso.sumChargedHadronPt*1./electron1.pt();                                                                
    lepton.pfRelIso03_NeutralHadron = pfIso.sumNeutralHadronEt*1./electron1.pt(); 
    lepton.pfRelIso03_Photon = pfIso.sumPhotonEt*1./electron1.pt();                                                                                                                         
    lepton.pfRelIso03_drcor = (pfIso.sumChargedHadronPt + max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5*pfIso.sumPUPt))*1./electron1.pt();   
  
    lepton.dr03EcalRecHitSumEt_Rel = electron1.dr03EcalRecHitSumEt()*1./electron1.pt(); 
    lepton.dr03HcalDepth1TowerSumEt_Rel = electron1.dr03HcalDepth1TowerSumEt()*1./electron1.pt(); 
    lepton.dr03HcalDepth2TowerSumEt_Rel = electron1.dr03HcalDepth2TowerSumEt()*1./electron1.pt();
    lepton.dr03TkSumPt_Rel =  electron1.dr03TkSumPt()*1./electron1.pt();
    lepton.dr03TkSumPtHEEP_Rel =  electron1.dr03TkSumPtHEEP()*1./electron1.pt();
    /*
    const reco::GsfElectron::PflowIsolationVariables& PfIso = electron1.isolationVariables04();//dr04IsolationVariables(); //pfIsolationVariables();  
    lepton.pfRelIso04_ChargedHadron = PfIso.sumChargedHadronPt*1./electron1.pt();                                                                
    lepton.pfRelIso04_NeutralHadron = PfIso.sumNeutralHadronEt*1./electron1.pt(); 
    lepton.pfRelIso04_Photon = PfIso.sumPhotonEt*1./electron1.pt();                                                                                                                         
    lepton.pfRelIso04_drcor = (PfIso.sumChargedHadronPt + max(0., PfIso.sumNeutralHadronEt + PfIso.sumPhotonEt - 0.5*PfIso.sumPUPt))*1./electron1.pt();   
    */
    vector<float> pfisovalues;                                                                                     
    Read_ElePFIsolation(&electron1,Rho,pfisovalues);
    lepton.pfRelIso03_eacor = pfisovalues[0];
    lepton.pfRelIso04_eacor = pfisovalues[1];
    
    // track info //
    
    lepton.chi2 = gsftrk1->normalizedChi2();                                                                                                                 
    lepton.ndof = (int)gsftrk1->ndof();                                                                                                            
	lepton.lostHits = (int)gsftrk1->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
	lepton.hit = (int)gsftrk1->hitPattern().numberOfValidMuonHits();                                                                              
	//lepton.nStations = (int)gsftrk1->numberOfMatchedStations();                                                                                          
	lepton.pixhit = (int)gsftrk1->hitPattern().numberOfValidPixelHits();                                                                          
	lepton.nTrackerLayers = (int)gsftrk1->hitPattern().trackerLayersWithMeasurement();                                                                    
	lepton.valfrac = gsftrk1->validFraction();   
	lepton.closeTrackNLayers = (int)electron1.closestCtfTrackNLayers();      
    lepton.closeTrackNormChi2 = electron1.closestCtfTrackNormChi2();   
    
    //MiniIsolation: begin//                                                                                      
	
	vector<float> isovalues;
	Read_ElePFIsolation(&electron1,Rho,isovalues);
	lepton.minisoall = isovalues[0];
	lepton.minchiso = isovalues[1];
	lepton.minnhiso = isovalues[2];
	lepton.minphiso = isovalues[3];
	
	// store 4-vector to vector of leptons //
	TLorentzVector p4;
	p4.SetPtEtaPhiM(lepton.pt, lepton.eta, lepton.phi, lepton.mass);
	tlep.push_back(p4);
	
	// now push lepton to vector of leptons //
	leptons.push_back(lepton);
	
    if(++nLepton>=njetmx) break;                                                                                                                      
  }
  
 // AK4 Jet //
  
  nPFJetAK4 = 0;
  
  for (unsigned jet = 0; jet< pfjetAK4s->size(); jet++) {
      
	const auto &ak4jet = (*pfjetAK4s)[jet];
    TLorentzVector pfjetAK4_4v(ak4jet.correctedP4("Uncorrected").px(),ak4jet.correctedP4("Uncorrected").py(),ak4jet.correctedP4("Uncorrected").pz(), ak4jet.correctedP4("Uncorrected").energy());
    
    double tmprecpt = pfjetAK4_4v.Pt();
   
    double total_cor =1;
    Read_JEC(total_cor,tmprecpt,pfjetAK4_4v.Eta(),Rho,isData,ak4jet,jecL1FastAK4,jecL2RelativeAK4,jecL3AbsoluteAK4,jecL2L3ResidualAK4);  
    PFJetAK4_JEC[nPFJetAK4] = total_cor;
    
    tmprecpt = pfjetAK4_4v.Pt();
    if(tmprecpt<minjPt) continue;
    if(abs(pfjetAK4_4v.Rapidity())>maxEta) continue;
      
    PFJetAK4_pt[nPFJetAK4] = 	tmprecpt;
    PFJetAK4_eta[nPFJetAK4] = 	pfjetAK4_4v.Eta();
    PFJetAK4_y[nPFJetAK4] = pfjetAK4_4v.Rapidity();
    PFJetAK4_phi[nPFJetAK4] = pfjetAK4_4v.Phi();
    PFJetAK4_mass[nPFJetAK4] = pfjetAK4_4v.M(); 
    
    // Jet id //
      
    JetIDVars AK4idvars;
      
    AK4idvars.NHF = ak4jet.neutralHadronEnergyFraction();
    AK4idvars.NEMF = ak4jet.neutralEmEnergyFraction();
    AK4idvars.MUF = ak4jet.muonEnergyFraction();
    AK4idvars.CHF = ak4jet.chargedHadronEnergyFraction();
    AK4idvars.CEMF = ak4jet.chargedEmEnergyFraction();
    AK4idvars.NumConst = (ak4jet.chargedMultiplicity()+ak4jet.neutralMultiplicity());
    AK4idvars.NumNeutralParticle = ak4jet.neutralMultiplicity();
    AK4idvars.CHM = ak4jet.chargedHadronMultiplicity();
     
    PFJetAK4_jetID[nPFJetAK4] = getJetID(AK4idvars,"CHS",year,PFJetAK4_eta[nPFJetAK4],false,isUltraLegacy);
    PFJetAK4_jetID_tightlepveto[nPFJetAK4] = getJetID(AK4idvars,"CHS",year,PFJetAK4_eta[nPFJetAK4],true,isUltraLegacy);
      
    PFJetAK4_hadronflav[nPFJetAK4] = ak4jet.hadronFlavour();
    PFJetAK4_partonflav[nPFJetAK4] = ak4jet.partonFlavour();
      
    PFJetAK4_qgl[nPFJetAK4] = ak4jet.userFloat("QGTagger:qgLikelihood");
    PFJetAK4_PUID[nPFJetAK4] = ak4jet.userFloat("pileupJetId:fullDiscriminant");
    
    // B tagging stuffs //
    
    PFJetAK4_btag_DeepCSV[nPFJetAK4] = ak4jet.bDiscriminator("pfDeepCSVJetTags:probb")+ak4jet.bDiscriminator("pfDeepCSVJetTags:probbb");
    PFJetAK4_btag_DeepFlav[nPFJetAK4] = ak4jet.bDiscriminator("pfDeepFlavourJetTags:probb") + ak4jet.bDiscriminator("pfDeepFlavourJetTags:probbb")+ak4jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
     
    nPFJetAK4++;	
    if(nPFJetAK4 >= njetmx) { break;}
    
  }
    
  // Skimming condition //
  
  if((nLepton)>=1){
	  
	for(unsigned ilep=0; ilep<leptons.size(); ilep++){
	
		if(save_only_muons) { 
			if (abs(leptons[ilep].pdgId) != 13) continue;
		}
		else if(save_only_electrons) { 
			if (abs(leptons[ilep].pdgId) != 11) continue;
		}
	
		lepton_genPartFlav = leptons[ilep].genPartFlav;
		lepton_pdgId = leptons[ilep].pdgId;
		
		// create generator labels //
		
		label_Muon_Prompt = label_Muon_fromTau = label_Muon_fromHadron = label_Muon_fromPhoton = label_Muon_unknown = false;
		label_Electron_Prompt = label_Electron_fromTau = label_Electron_fromHadron = label_Electron_fromPhoton = label_Electron_unknown = false;
		
		if(abs(leptons[ilep].pdgId)==13){
			if(leptons[ilep].genPartFlav==1) { label_Muon_Prompt = true; }
			else if (leptons[ilep].genPartFlav==15) { label_Muon_fromTau = true; }
			else if (leptons[ilep].genPartFlav==5 || leptons[ilep].genPartFlav==4) { label_Muon_fromHadron = true; }
			else if (leptons[ilep].genPartFlav==22) { label_Muon_fromPhoton = true; }
			else { label_Muon_unknown = true; }
		}
		
		if(abs(leptons[ilep].pdgId)==11){
			if(leptons[ilep].genPartFlav==1) { label_Electron_Prompt = true; }
			else if (leptons[ilep].genPartFlav==15) { label_Electron_fromTau = true; }
			else if (leptons[ilep].genPartFlav==5 || leptons[ilep].genPartFlav==4) { label_Electron_fromHadron = true; }
			else if (leptons[ilep].genPartFlav==22) { label_Electron_fromPhoton = true; }
			else { label_Electron_unknown = true; }
		}
		
		// created generator labels //
	
		lepton_pt = leptons[ilep].pt;
		lepton_eta = leptons[ilep].eta;
		lepton_phi = leptons[ilep].phi;
		lepton_mass = leptons[ilep].mass;		
		lepton_p = leptons[ilep].p ;
	
		lepton_charge = leptons[ilep].charge;
		lepton_tightcharge = leptons[ilep].tightcharge;
		
		lepton_dxy = leptons[ilep].dxy;
		lepton_dz = leptons[ilep].dz;
		lepton_dxyError = leptons[ilep].dxyError;
		lepton_dzError = leptons[ilep].dzError;
		lepton_ip3d = leptons[ilep].ip3d;
		lepton_sip3d = leptons[ilep].sip3d;
		lepton_dxy_sv = leptons[ilep].dxy_sv;
  
		lepton_chi2 = leptons[ilep].chi2;
		lepton_ndof = leptons[ilep].ndof;
		lepton_trkKink = leptons[ilep].trkKink;
		lepton_hit = leptons[ilep].hit;
		lepton_pixhit = leptons[ilep].pixhit;
		lepton_nTrackerLayers = leptons[ilep].nTrackerLayers;
		lepton_lostHits = leptons[ilep].lostHits;
    
		lepton_e_ECAL = leptons[ilep].e_ECAL;
		lepton_e_HCAL = leptons[ilep].e_HCAL;
		lepton_hoe = leptons[ilep].hoe;
		  
		lepton_minchiso = leptons[ilep].minchiso;
		lepton_minnhiso = leptons[ilep].minnhiso;
		lepton_minphiso = leptons[ilep].minphiso;
		lepton_minisoall = leptons[ilep].minisoall;
		
		lepton_pfRelIso03_drcor = leptons[ilep].pfRelIso03_drcor;
		lepton_pfRelIso03_ChargedHadron = leptons[ilep].pfRelIso03_ChargedHadron;
		lepton_pfRelIso03_NeutralHadron = leptons[ilep].pfRelIso03_NeutralHadron;
		lepton_pfRelIso03_Photon = leptons[ilep].pfRelIso03_Photon;
		lepton_pfRelIso03_PileUp = leptons[ilep].pfRelIso03_PileUp;
		lepton_tkRelIso = leptons[ilep].tkRelIso;
		
		lepton_nStations = leptons[ilep].nStations;
		lepton_segmentComp = leptons[ilep].segmentComp;
    
		lepton_isPFCand = leptons[ilep].isPFCand;
		lepton_isGlobal = leptons[ilep].isGlobal;
		lepton_isTracker = leptons[ilep].isTracker;
		lepton_isLoose = leptons[ilep].isLoose;
		lepton_isGoodGlobal = leptons[ilep].isGoodGlobal;
		lepton_isMedium = leptons[ilep].isMedium;
		lepton_isMedPr = leptons[ilep].isMedPr;
		lepton_isTight = leptons[ilep].isTight;
		lepton_isHighPt = leptons[ilep].isHighPt;
		lepton_isHighPttrk = leptons[ilep].isHighPttrk;

		lepton_posmatch = leptons[ilep].posmatch;
    
		lepton_mvaid_Fallv2WP90 = leptons[ilep].mvaid_Fallv2WP90;
		lepton_mvaid_Fallv2WP90_noIso = leptons[ilep].mvaid_Fallv2WP90_noIso;
		lepton_mvaid_Fallv2WP80 = leptons[ilep].mvaid_Fallv2WP80;
		lepton_mvaid_Fallv2WP80_noIso = leptons[ilep].mvaid_Fallv2WP80_noIso;
  
		lepton_eoverp = leptons[ilep].eoverp;
 
		lepton_pfRelIso03_eacor = leptons[ilep].pfRelIso03_eacor ;
		lepton_pfRelIso04_eacor = leptons[ilep].pfRelIso04_eacor;
		lepton_dr03EcalRecHitSumEt_Rel = leptons[ilep].dr03EcalRecHitSumEt_Rel;
		lepton_dr03HcalDepth1TowerSumEt_Rel = leptons[ilep].dr03HcalDepth1TowerSumEt_Rel;
		lepton_dr03HcalDepth2TowerSumEt_Rel = leptons[ilep].dr03HcalDepth2TowerSumEt_Rel;
		lepton_dr03TkSumPt_Rel = leptons[ilep].dr03TkSumPt_Rel;
		lepton_dr03TkSumPtHEEP_Rel = leptons[ilep].dr03TkSumPtHEEP_Rel;
  
		lepton_eInvMinusPInv = leptons[ilep].eInvMinusPInv;
		lepton_supcl_eta = leptons[ilep].supcl_eta;
		lepton_supcl_phi = leptons[ilep].supcl_phi;
		lepton_supcl_energy = leptons[ilep].supcl_energy;
		lepton_sigmaietaieta = leptons[ilep].sigmaietaieta;
		lepton_sigmaiphiiphi = leptons[ilep].sigmaiphiiphi;
		lepton_r9full = leptons[ilep].r9full;
		lepton_hcaloverecal = leptons[ilep].hcaloverecal;
		lepton_ecloverpout = leptons[ilep].ecloverpout;
		lepton_convVeto = leptons[ilep].convVeto;

		lepton_etain = leptons[ilep].etain;
		lepton_phiin = leptons[ilep].phiin;
		lepton_fbrem = leptons[ilep].fbrem;
		lepton_supcl_etaWidth = leptons[ilep].supcl_etaWidth;
		lepton_supcl_phiWidth = leptons[ilep].supcl_phiWidth;
  
		lepton_e1x5bye5x5 = leptons[ilep].e1x5bye5x5;
		lepton_convtxprob = leptons[ilep].convtxprob;
		lepton_deltaetacltrkcalo = leptons[ilep].deltaetacltrkcalo;
		lepton_supcl_preshvsrawe = leptons[ilep].supcl_preshvsrawe;
  
		lepton_closeTrackNLayers = leptons[ilep].closeTrackNLayers;
		lepton_closeTrackNormChi2 = leptons[ilep].closeTrackNormChi2;
	
		// PF candidates within dR<0.5 around lepton //
		
		nPFCand = 0;  
		
		vector<PFlowCandidate> pfcandidates;                                                                                                                                      
  
		if(pfcands.isValid() && pfcands->size()>0) {                                                                                                           
    
			edm::View<pat::PackedCandidate>::const_iterator pfcand1;     
			
			float dRmax = 0.5;  
			float dRmin = -1.e-3;                                                                                               

			for( pfcand1 = pfcands->begin(); pfcand1 < pfcands->end(); pfcand1++ ) {    
		
				if(delta2R(leptons[ilep].eta,leptons[ilep].phi,pfcand1->eta(),pfcand1->phi())<dRmax && delta2R(leptons[ilep].eta,leptons[ilep].phi,pfcand1->eta(),pfcand1->phi())>dRmin){
				
					PFlowCandidate pfcand;
					Initialize_PFlowCandidate(pfcand);
		
					pfcand.pt = pfcand1->pt();
					pfcand.eta = pfcand1->eta();
					pfcand.phi = pfcand1->phi();
					pfcand.mass = pfcand1->mass();
					pfcand.phiAtVtx = pfcand1->phiAtVtx();
					pfcand.pdgId = pfcand1->pdgId();
					pfcand.puppiWeight = pfcand1->puppiWeight();
					pfcand.puppiWeightNoLep = pfcand1->puppiWeightNoLep();
		
					pfcand.status = pfcand1->status();
					pfcand.caloFraction = pfcand1->caloFraction();
					pfcand.hcalFraction = pfcand1->hcalFraction();
		
					if(pfcand1->hasTrackDetails()){
			
						pfcand.time = pfcand1->time();
						pfcand.charge = pfcand1->charge();
						pfcand.dz = pfcand1->dz();
						pfcand.dzError = pfcand1->dzError();
						pfcand.dxy = pfcand1->dxy();
						pfcand.dxyError = pfcand1->dxyError();
						pfcand.vertexChi2 = pfcand1->vertexChi2();
						pfcand.lostInnerHits = pfcand1->lostInnerHits();
						pfcand.pvAssocQuality = pfcand1->pvAssociationQuality();
						pfcand.trackHighPurity = pfcand1->trackHighPurity();
						pfcand.pixelhits = pfcand1->numberOfPixelHits();
			
						const reco::Track *trktrk = pfcand1->bestTrack();
						pfcand.trkChi2 = trktrk->normalizedChi2();
						pfcand.nTrackerLayers = trktrk->hitPattern().trackerLayersWithMeasurement();
	
						//cout<<pfcand1->ptTrk()<<" "<<<endl;
						//cout<<pfcand1->vertexNormalizedChi2()<<"\t"
						//<<pfcand1->vertexNdof()<<"\t"
						//<<pfcand1->rawCaloFraction()<<"\t"
						//<<pfcand1->isIsolatedChargedHadron()<<"\t"
	
					}
					
					pfcandidates.push_back(pfcand);
			
					if (++nPFCand>=nconsmax) break;  
				
				} // DR cond
		
		
			}// loop over pfcands
	
		}
		
		ChargePFCand_pt_rel.clear();
		ChargePFCand_deta.clear();
		ChargePFCand_dphi.clear();
		ChargePFCand_dphiAtVtx.clear();
		ChargePFCand_puppiWeight.clear();
		ChargePFCand_puppiWeightNoLep.clear();
		ChargePFCand_caloFraction.clear();
		ChargePFCand_hcalFraction.clear();
		ChargePFCand_pdgId.clear();
		
		ChargePFCand_dz.clear();
		ChargePFCand_dzError.clear();
		ChargePFCand_dxy.clear();
		ChargePFCand_dxyError.clear();
		ChargePFCand_trkChi2.clear();
		ChargePFCand_vertexChi2.clear();
		ChargePFCand_charge.clear();
		ChargePFCand_pvAssocQuality.clear();
		ChargePFCand_trkQuality.clear();
		ChargePFCand_nTrackerLayers.clear();
		ChargePFCand_pixelhits.clear();
		ChargePFCand_status.clear();
		ChargePFCand_time.clear();
		ChargePFCand_trackHighPurity.clear();
		
		NeutralPFCand_pt_rel.clear();
		NeutralPFCand_deta.clear();
		NeutralPFCand_dphi.clear();
		NeutralPFCand_dphiAtVtx.clear();
		NeutralPFCand_puppiWeight.clear();
		NeutralPFCand_puppiWeightNoLep.clear();
		NeutralPFCand_caloFraction.clear();
		NeutralPFCand_hcalFraction.clear();
		NeutralPFCand_pdgId.clear();
		
		for(unsigned ipf=0; ipf<pfcandidates.size(); ipf++){
			
			if(abs(pfcandidates[ipf].charge)>0){
			
				ChargePFCand_pt_rel.push_back(pfcandidates[ipf].pt/lepton_pt);
				ChargePFCand_deta.push_back(pfcandidates[ipf].eta - lepton_eta);
				ChargePFCand_dphi.push_back(PhiInRange(pfcandidates[ipf].phi - lepton_phi));
				ChargePFCand_dphiAtVtx.push_back(PhiInRange(pfcandidates[ipf].phiAtVtx - lepton_phi));;
				ChargePFCand_puppiWeight.push_back(pfcandidates[ipf].puppiWeight);
				ChargePFCand_puppiWeightNoLep.push_back(pfcandidates[ipf].puppiWeightNoLep);
				ChargePFCand_caloFraction.push_back(pfcandidates[ipf].caloFraction);
				ChargePFCand_hcalFraction.push_back(pfcandidates[ipf].hcalFraction);
				ChargePFCand_pdgId.push_back(pfcandidates[ipf].pdgId);
  
				ChargePFCand_dz.push_back(pfcandidates[ipf].dz);
				ChargePFCand_dzError.push_back(pfcandidates[ipf].dzError);
				ChargePFCand_dxy.push_back(pfcandidates[ipf].dxy);
				ChargePFCand_dxyError.push_back(pfcandidates[ipf].dxyError);
				ChargePFCand_trkChi2.push_back(pfcandidates[ipf].trkChi2);
				ChargePFCand_vertexChi2.push_back(pfcandidates[ipf].vertexChi2);
				ChargePFCand_charge.push_back(pfcandidates[ipf].charge);
				ChargePFCand_lostInnerHits.push_back(pfcandidates[ipf].lostInnerHits);
				ChargePFCand_pvAssocQuality.push_back(pfcandidates[ipf].pvAssocQuality);
				ChargePFCand_trkQuality.push_back(pfcandidates[ipf].trkQuality);
				ChargePFCand_nTrackerLayers.push_back(pfcandidates[ipf].nTrackerLayers);
				ChargePFCand_pixelhits.push_back(pfcandidates[ipf].pixelhits);
				ChargePFCand_status.push_back(pfcandidates[ipf].status);
				ChargePFCand_time.push_back(pfcandidates[ipf].time);
				ChargePFCand_trackHighPurity.push_back(pfcandidates[ipf].trackHighPurity);
			
			}
			
			else{
			
				NeutralPFCand_pt_rel.push_back(pfcandidates[ipf].pt/lepton_pt);
				NeutralPFCand_deta.push_back(pfcandidates[ipf].eta - lepton_eta);
				NeutralPFCand_dphi.push_back(PhiInRange(pfcandidates[ipf].phi - lepton_phi));
				NeutralPFCand_dphiAtVtx.push_back(PhiInRange(pfcandidates[ipf].phiAtVtx - lepton_phi));;
				NeutralPFCand_puppiWeight.push_back(pfcandidates[ipf].puppiWeight);
				NeutralPFCand_puppiWeightNoLep.push_back(pfcandidates[ipf].puppiWeightNoLep);
				NeutralPFCand_caloFraction.push_back(pfcandidates[ipf].caloFraction);
				NeutralPFCand_hcalFraction.push_back(pfcandidates[ipf].hcalFraction);
				NeutralPFCand_pdgId.push_back(pfcandidates[ipf].pdgId);
				
			}
  
		}
		
		nChargePFCand = (int)ChargePFCand_deta.size();
		nNeutralPFCand = (int)NeutralPFCand_deta.size();
		
		nSV = 0;
		
		SV_x.clear(); SV_y.clear(); SV_z.clear(); 
		SV_ndof.clear(); SV_chi2.clear();
		SV_mass.clear(), SV_pt_rel.clear();
		SV_dlen.clear(), SV_dlenSig.clear();
		SV_dxy.clear(); SV_dxySig.clear();
		SV_pAngle.clear();
		
		for(unsigned int isv=0; isv<(secondaryVertices->size()); isv++){                                                                                        
	
			const auto &sv = (*secondaryVertices)[isv];                                                                                                           
			reco::TrackBase::Point svpoint(sv.vx(),sv.vy(),sv.vz());
	
			float dRmax = 0.5;
		
			if(delta2R(leptons[ilep].eta,leptons[ilep].phi,sv.eta(),sv.phi())<dRmax){
			
				SV_x.push_back(sv.vx());
				SV_y.push_back(sv.vy());
				SV_z.push_back(sv.vz());
				SV_ndof.push_back(sv.vertexNdof());
				SV_chi2.push_back(sv.vertexNormalizedChi2());
				SV_mass.push_back(sv.mass());
				SV_pt_rel.push_back(sv.pt()/lepton_pt);
			
				VertexDistance3D vdist;
				Measurement1D dl = vdist.distance(PV_vertex, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
				SV_dlen.push_back(dl.value());
				SV_dlenSig.push_back(dl.significance());
			
				VertexDistanceXY vdistXY;
				Measurement1D d2d = vdistXY.distance(PV_vertex, VertexState(RecoVertex::convertPos(sv.position()), RecoVertex::convertError(sv.error())));
				SV_dxy.push_back(d2d.value());
				SV_dxySig.push_back(d2d.significance());
			
				double dx = (PV_vertex.x() - sv.vx()), dy = (PV_vertex.y() - sv.vy()), dz = (PV_vertex.z() - sv.vz());
				double pdotv = (dx * sv.px() + dy * sv.py() + dz * sv.pz()) / sv.p() / sqrt(dx * dx + dy * dy + dz * dz);
				SV_pAngle.push_back( acos(pdotv));
		
				if (++nSV>=nconsmax) break;
				  
			} // DR condition
			
		}// loop over secondary vertices	
		
		float dR_lj_min = 0.4;
		int i_nearjet = -1;
		
		for(int ijet=0; ijet<nPFJetAK4; ijet++){
			if(delta2R(lepton_eta,lepton_phi,PFJetAK4_eta[ijet],PFJetAK4_phi[ijet])<dR_lj_min){
				dR_lj_min = delta2R(lepton_eta,lepton_phi,PFJetAK4_eta[ijet],PFJetAK4_phi[ijet]);
				i_nearjet = ijet;
			}
		}
		
		if(i_nearjet>=0){
		
			TLorentzVector jet_p4;
			jet_p4.SetPtEtaPhiM(PFJetAK4_pt[i_nearjet],PFJetAK4_eta[i_nearjet],PFJetAK4_phi[i_nearjet],PFJetAK4_mass[i_nearjet]);
			TLorentzVector lep_p4;
			lep_p4.SetPtEtaPhiM(lepton_pt,lepton_eta,lepton_phi,lepton_mass);
			jet_p4 -= lep_p4;
			
			lepton_jetPtRelv2 = (lep_p4.Vect().Perp(jet_p4.Vect()))*1./lepton_pt;
			lepton_jetRelIso = PFJetAK4_pt[i_nearjet]*1./lepton_pt;
			lepton_jetbtag = PFJetAK4_btag_DeepFlav[i_nearjet];
			
		}
		else{
			
			lepton_jetPtRelv2 = -100;
			lepton_jetRelIso = -100;
			lepton_jetbtag = -100;
			
			}
  	
		T1->Fill();
		
	} // ilep

  }
  // End of skimming 
  
}
