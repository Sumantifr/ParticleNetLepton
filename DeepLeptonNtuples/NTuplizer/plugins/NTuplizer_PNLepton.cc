// -*- C++ -*-
//
// Package:    Run2_2016/TopplusB
// Class:      TopplusB
// 
/**\class NTuplizer_XYH NTuplizer_XYH.cc 
   
   Description: [one line class summary]
   
   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Suman Chatterjee
//         Created:  Fri, 1 Oct 2021 16:22:44 GMT
//

//Twikis used:
/*
Prefiring weights: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe#2018_UL 
Electron MVA ID: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
Photon MVA ID: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2 
Main EGamma: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#MVA_based_electron_Identificatio
JEC: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources
DeepAKX: https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepAKXTagging
Btag SF (recipe): https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration
Btag SF (2018UL): https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
Rochester correction: https://gitlab.cern.ch/akhukhun/roccor
*/

// system include files
#include <memory>
// user include files
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TAxis.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TRandom.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "GeneratorInterface/Pythia8Interface/plugins/ReweightUserHooks.h"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include <string>
#include <iostream>
#include <fstream>
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include  "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "RecoEgamma/ElectronIdentification/interface/ElectronMVAEstimatorRun2.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReportEntry.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReport.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include "FWCore/Utilities/interface/typelookup.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/GhostedAreaSpec.hh>
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/SecondaryVertexConverter.h"

// Rochester correction for muons //
#include "RoccoR.h"

//for storing vectors in tree//

# include <vector>


#ifdef __MAKECINT__
    
    #pragma link C++ class std::vector+;
    #pragma link C++ class std::vector<float>+;
    #pragma link C++ class std::vector<int>+;
    #pragma link C++ class std::vector<bool>+;
    #pragma link C++ class std::vector<std::vector<float> >+;
    
#endif


using namespace std;
using namespace edm;
using namespace reco;  
using namespace CLHEP;
using namespace trigger;
using namespace math;
using namespace fastjet;
using namespace fastjet::contrib;
using namespace btagbtvdeep;

const float mu_mass = 0.105658;
const float el_mass = 0.000511;
const float pival = acos(-1.);

static const int nsrc = 24;
const char* jecsrcnames[nsrc] = {
	 "AbsoluteStat", "AbsoluteScale","AbsoluteMPFBias", 
	 "FlavorQCD", "Fragmentation", 
	 "PileUpDataMC",  "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", //"PileUpPtHF",
	 "PileUpPtRef",
	 "RelativeFSR", "RelativeJEREC1", "RelativeJEREC2", //"RelativeJERHF",
	 "RelativePtBB", "RelativePtEC1", "RelativePtEC2", //"RelativePtHF", 
	 "RelativeBal", "RelativeSample", "RelativeStatEC", "RelativeStatFSR", //"RelativeStatHF", 
	 "SinglePionECAL", "SinglePionHCAL","TimePtEta",
	 "Total"
	};
const int njecmcmx = 2*nsrc + 1 ;

struct triggervar{
  TLorentzVector  trg4v;
  bool		  both;
  bool            level1;
  bool            highl;
  int             ihlt;
  int             prescl;
  int             pdgId;
  int			  type;
};

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double theta_to_eta(double theta) { return -log(tan(theta/2.)); }

double PhiInRange(const double& phi) {
  double phiout = phi;
  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;
  return phiout;
}

double delta2R(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1 - eta2,2) +pow(PhiInRange(phi1 - phi2),2));
}

class GenParton{

 public:

  float  eta;
  float  mass;
  float  phi;
  float  pt;

  int status;
  int pdgId;
  int mompdgId;
  int momstatus;
  int grmompdgId;
  
  bool fromhard;
  bool fromhardbFSR;
  bool isPromptFinalState;
  bool isLastCopyBeforeFSR;
  bool isDirectPromptTauDecayProductFinalState;
  
  TLorentzVector p4;

} ;

class Lepton{

	public:
	
		float pt;
		float eta;
		float phi;
		float mass;
		float p;
		int pdgId;
		int charge;
		int tightcharge;
		//distance from PV & SV//
		float dxy;
		float dxyError;
		float dz;
		float dzError;
		float ip3d;
		float sip3d;
		float dxy_sv;
		//Track information //
		float chi2;
		float valfrac;
		int ndof;
		int trkKink;
		int hit;
		int pixhit;
		int nTrackerLayers;
		int lostHits;
		//Energy information //
		float e_ECAL;
		float e_HCAL;
		float hoe;
		//Mini-isolation variables //
		float minchiso;
		float minnhiso;
		float minphiso;
		float minisoall; 
		//PF-isolation variables //
		float pfRelIso03_drcor;
		float pfRelIso03_ChargedHadron;
		float pfRelIso03_NeutralHadron;
		float pfRelIso03_Photon;
		float pfRelIso03_PileUp;
		float tkRelIso;
		float pfRelIso04_drcor;
		float pfRelIso04_ChargedHadron;
		float pfRelIso04_NeutralHadron;
		float pfRelIso04_Photon;
		float pfRelIso04_PileUp;
		// GEN matching //
		int genPartFlav;
  
		// Electron-specific variables //
  
		// ID booleans //
		bool mvaid_Fallv2WP90_noIso;
		bool mvaid_Fallv2WP90;
		bool mvaid_Fallv2WP80_noIso;
		bool mvaid_Fallv2WP80;
		// isolation //
		float pfRelIso03_eacor;
		float pfRelIso04_eacor;
		// energy cluster //
		float dr03EcalRecHitSumEt_Rel;
		float dr03HcalDepth1TowerSumEt_Rel;
		float dr03HcalDepth2TowerSumEt_Rel;
		float dr03TkSumPt_Rel;
		float dr03TkSumPtHEEP_Rel;
		float eoverp;
		float etain;
		float phiin;
		float fbrem; 
		// super-clyster// 
		float supcl_eta;
		float supcl_phi;
		float supcl_energy; 
		float r9full;
		float supcl_etaWidth;
		float supcl_phiWidth;
		float hcaloverecal;
		int  closeTrackNLayers;
		float closeTrackNormChi2;
		float e1x5bye5x5;
		float normchi2;
		float convtxprob;
		float ecloverpout;
		float deltaetacltrkcalo;
		float supcl_preshvsrawe;
		float sigmaietaieta;
		float sigmaiphiiphi;
		float eInvMinusPInv;
		bool convVeto;
  
		//Muon-specific variables //
  
		// ID booleans //
		bool isPFCand;
		bool isGlobal;
		bool isTracker;
		bool isGoodGlobal;
		bool isTight;
		bool isHighPt;
		bool isHighPttrk;
		bool isMedium;
		bool isMedPr;
		bool isLoose;
		//track//
		float posmatch;
		float segmentComp;
		int nStations;
  
};

void Initialize(Lepton &lepton){
	
	lepton.pt = -100 ;
	lepton.eta = -100 ;
	lepton.phi = -100 ;
	lepton.mass = -100 ;
	lepton.p = -100 ;
	
	lepton.pdgId = -100 ;
	lepton.charge = -100 ;
	lepton.tightcharge = -100 ;
		
	lepton.dxy = -100;
	lepton.dz = -100;
	lepton.dxyError = -100;
	lepton.dzError = -100;
	lepton.ip3d = -100;
	lepton.sip3d = -100;
	lepton.dxy_sv = -100;
		
	lepton.genPartFlav = -100;
  
  	lepton.chi2 = -100;
	lepton.ndof = -100;
	lepton.trkKink = -100;
	lepton.hit = -100;
	lepton.pixhit = -100;
	lepton.nTrackerLayers = -100;
	lepton.lostHits = -100;
    
	lepton.e_ECAL = -100;
	lepton.e_HCAL = -100;
	lepton.hoe = -100;
		  
	lepton.minchiso = -100;
	lepton.minnhiso = -100;
	lepton.minphiso = -100;
	lepton.minisoall = -100;
		
	lepton.pfRelIso03_drcor = -100;
	lepton.pfRelIso03_ChargedHadron = -100;
	lepton.pfRelIso03_NeutralHadron = -100;
	lepton.pfRelIso03_Photon = -100;
	lepton.pfRelIso03_PileUp = -100;
	lepton.tkRelIso = -100;
   
	lepton.mvaid_Fallv2WP90 = false;
	lepton.mvaid_Fallv2WP90_noIso = false;
	lepton.mvaid_Fallv2WP80 = false;
	lepton.mvaid_Fallv2WP80_noIso = false;
   
	lepton.eoverp = -100;
 
	lepton.pfRelIso03_eacor = -100;
	lepton.pfRelIso04_eacor = -100;
	lepton.dr03EcalRecHitSumEt_Rel = -100;
	lepton.dr03HcalDepth1TowerSumEt_Rel = -100;
	lepton.dr03HcalDepth2TowerSumEt_Rel = -100;
	lepton.dr03TkSumPt_Rel = -100;
	lepton.dr03TkSumPtHEEP_Rel = -100;
    
	lepton.eInvMinusPInv = -100;
	lepton.supcl_eta = -100;
	lepton.supcl_phi = -100;
	lepton.supcl_energy = -100;
	lepton.sigmaietaieta = -100;
	lepton.sigmaiphiiphi = -100;
	lepton.r9full = -100;
	lepton.hcaloverecal = -100;
	lepton.ecloverpout = -100;
	lepton.convVeto = -100;
	lepton.etain = -100;
	lepton.phiin = -100;
	lepton.fbrem = -100;
	lepton.supcl_etaWidth = -100;
	lepton.supcl_phiWidth = -100;
 
	lepton.e1x5bye5x5 = -100;
	lepton.convtxprob = -100;
	lepton.deltaetacltrkcalo = -100;
	lepton.supcl_preshvsrawe = -100;
 
	lepton.closeTrackNLayers = -100;
	lepton.closeTrackNormChi2 = -100;
	
	lepton.isPFCand = false;
	lepton.isGlobal = false;
	lepton.isTracker = false;
	lepton.isLoose = false;
	lepton.isGoodGlobal = false;
	lepton.isMedium = false;
	lepton.isMedPr = false;
	lepton.isTight = false;
	lepton.isHighPt = false;
	lepton.isHighPttrk = false;
  
	lepton.posmatch = -100;
	lepton.nStations = -100;
	lepton.segmentComp = -100;
 
};

class PFlowCandidate{
	
	public: 
		float pt;
		float eta;
		float phi;
		float mass;
		float phiAtVtx;
		float dz;
		float dzError;
		float dxy;
		float dxyError;
		float puppiWeight;
		float puppiWeightNoLep;
		float trkChi2;
		float vertexChi2;
		float caloFraction;
		float hcalFraction;
		int pdgId;
		int charge;
		int lostInnerHits;
		float pvAssocQuality;
		float trkQuality;
		int nTrackerLayers;
		int pixelhits;
		int status;
		int time;
		bool trackHighPurity;
 
};

void Initialize_PFlowCandidate(PFlowCandidate &pfcand){

		pfcand.pt = -100;
		pfcand.eta = -100;
		pfcand.phi = -100;
		pfcand.mass = -100;
		pfcand.phiAtVtx = -100;
		pfcand.dz = -100;
		pfcand.dzError = -100;
		pfcand.dxy = -100;
		pfcand.dxyError = -100;
		pfcand.puppiWeight = -100;
		pfcand.puppiWeightNoLep = -100;
		pfcand.trkChi2 = -100;
		pfcand.vertexChi2 = -100;
		pfcand.caloFraction = -100;
		pfcand.hcalFraction = -100;
		pfcand.pdgId = -100;
		pfcand.charge = -100;
		pfcand.lostInnerHits = -100;
		pfcand.pvAssocQuality = -100;
		pfcand.trkQuality = -100;
		pfcand.nTrackerLayers = -100;
		pfcand.pixelhits = -100;
		pfcand.status = -100;
		pfcand.time = -100;
		pfcand.trackHighPurity = false;
};

int getGenPartonFlavor(float lep_eta, float lep_phi, vector<GenParton>GenParts, int lep_pdgId, float dRcut=0.4)
{
	bool gen_match = false;
	unsigned i_gen_match;
	
	for(unsigned igen=0; igen<GenParts.size(); igen++){
		if(delta2R(lep_eta,lep_phi,GenParts[igen].eta,GenParts[igen].phi)<dRcut && abs(GenParts[igen].pdgId)==lep_pdgId){
			if(GenParts[igen].status==1){
				gen_match = true;
				i_gen_match = igen;
			}
		}
	}
	
	int mom_pdgId = (GenParts[i_gen_match].momstatus == 2) ? (abs(GenParts[i_gen_match].grmompdgId)) : (abs(GenParts[i_gen_match].mompdgId));
	int genPartFlav = -100;
	
	if(gen_match){	
		//if (mom_pdgId==24 || mom_pdgId==23 || mom_pdgId==25 || 
		if (GenParts[i_gen_match].isPromptFinalState){
			if (mom_pdgId==22){
				genPartFlav = 22;
			}
			else{
				genPartFlav = 1;
			}
		}
		//else if (mom_pdgId==15 ||
		else if (GenParts[i_gen_match].isDirectPromptTauDecayProductFinalState){
			genPartFlav = 15;
		} 
		else if (mom_pdgId==22){
			genPartFlav = 22;
		} 
		else if ( (mom_pdgId/1000 == 5) || (mom_pdgId/100 == 5) || (mom_pdgId == 5)){
			genPartFlav = 5;
		} 
		else if ( (mom_pdgId/1000 == 4) || (mom_pdgId/100 == 4) || (mom_pdgId == 4)){
			genPartFlav = 4;
		} 
		else {
			genPartFlav = 3;
		}
	}
	else{
			genPartFlav = -100;
	}
	
	return 	genPartFlav;
}

struct JetIDVars
{
  float NHF, NEMF, MUF, CHF, CEMF;
  int NumConst, NumNeutralParticle, CHM;
};

bool getJetID(JetIDVars vars, string jettype="CHS", int year=2018, double eta=0, bool tightLepVeto=true, bool UltraLegacy=false){
  
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
  
  if (jettype!="CHS" && jettype!="PUPPI"){
    cout<<"Don't know your jet type! I know only CHS & PUPPI :D"<<endl;
    return false;
  }
  
  float NHF, NEMF, MUF, CHF, CEMF;
  int NumConst, NumNeutralParticle, CHM;
  
  NHF = vars.NHF; 
  NEMF = vars.NEMF;
  MUF = vars.MUF;
  CHF = vars.CHF;
  CEMF = vars.CEMF;
  NumConst = vars.NumConst;
  NumNeutralParticle = vars.NumNeutralParticle;
  CHM = vars.CHM;
  
  bool JetID = false;
  
  if(!UltraLegacy){
    
    if(year==2018 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10));
    }
    
    if(year==2018 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)<=3.0 && NHF<0.99) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>2 && NumNeutralParticle<15));
    }
    
    if(year==2017 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>10));
    }
    
    if(year==2017 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) ||
 (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NHF<0.99) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>2 && NumNeutralParticle<15));
    }

    if(year==2016 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.90 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CEMF<0.99 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9  && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.01 && NHF<0.98 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>10));
	}
    
    if(year==2016 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.9 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CEMF<0.99 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ));
      if(fabs(eta)>2.7) { JetID = false; }
	}
  }
  
  else {
    
    if(year==2017||year==2018){
      
      if(jettype=="CHS"){
	
	JetID = ( fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.01 && NEMF<0.99 && NumNeutralParticle>1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10) ;
      }
      
      if(jettype=="PUPPI"){
	
	JetID =  ( fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NHF<0.9999 ) ||( fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>2 ) ;
      }
      // there is a inconsistency between table & lines in https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
      // table is chosen as it is consistent with the slides https://indico.cern.ch/event/937597/contributions/3940302/attachments/2073315/3481068/ULJetID_UL17_UL18_AK4PUPPI.pdf 
    }
    
    if(year==2016){
      
      if(jettype=="CHS"){
	
	JetID =  ( fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.4 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticle>1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10) ;
	
      }
      
      if(jettype=="PUPPI"){
	
	JetID = ( fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.4 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.98 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NumNeutralParticle>=1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>2  ) ;
      }
    }	
  }
  
  return JetID;
  
}

bool Muon_Tight_ID(bool muonisGL,bool muonisPF, float muonchi, float muonhit, float muonmst, float muontrkvtx, float muondz, float muonpixhit, float muontrklay){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Tight_Muon
	bool tightid = false;
	if(muonisGL && muonisPF){
		if(muonchi<10 && muonhit>0 && muonmst>1){
			if(fabs(muontrkvtx)<0.2 && fabs(muondz)<0.5){
				if(muonpixhit>0 && muontrklay>5){
					tightid = true;
				}
			}
		}
	}
	return tightid;
}

bool StoreMuon(pat::Muon muon1, float ptcut, float etacut){
	
	if (((muon1.isTrackerMuon() || muon1.isGlobalMuon()) && (muon1.isPFMuon())) && (muon1.pt()>=ptcut) && (fabs(muon1.eta())<=etacut)) {                                                                
			return true;
	}
	else{
			return false;
		}
}

bool StoreElectron(pat::Electron electron1, float ptcut, float etacut){
	
	GsfTrackRef gsftrk1 = electron1.gsfTrack();                                                                                                      
    if ((!gsftrk1.isNull()) && (electron1.pt()>=ptcut) && (fabs(electron1.eta())<=etacut) && (gsftrk1->ndof()>=9)) {
			return true;
		}
    else{
			return false;
		}
}

void Read_JEC(double &total_JEC,  double &tmprecpt, 
			  double jeteta, double Rho, bool isData,
			  pat::Jet jet,
			  FactorizedJetCorrector *jecL1Fast, FactorizedJetCorrector *jecL2Relative, FactorizedJetCorrector *jecL3Absolute, FactorizedJetCorrector*jecL2L3Residual)
{
	
    double total_cor =1;
      
    jecL1Fast->setJetPt(tmprecpt); jecL1Fast->setJetA(jet.jetArea()); jecL1Fast->setRho(Rho);jecL1Fast->setJetEta(jeteta);
    double corFactorL1Fast = jecL1Fast->getCorrection();
    total_cor *= corFactorL1Fast;
    tmprecpt = tmprecpt * corFactorL1Fast;
      
    jecL2Relative->setJetPt(tmprecpt); jecL2Relative->setJetEta(jeteta);
    double corFactorL2Relative = jecL2Relative->getCorrection();
    total_cor *= corFactorL2Relative ;
    tmprecpt = tmprecpt * corFactorL2Relative;
      
    jecL3Absolute->setJetPt(tmprecpt); jecL3Absolute->setJetEta(jeteta);
    double corFactorL3Absolute = jecL3Absolute->getCorrection();
    total_cor *= corFactorL3Absolute ;
    tmprecpt = tmprecpt * corFactorL3Absolute;
      
    double corFactorL2L3Residual=1.;
      
    if(isData){
		jecL2L3Residual->setJetPt(tmprecpt); jecL2L3Residual->setJetEta(jeteta);
		corFactorL2L3Residual = jecL2L3Residual->getCorrection();
		total_cor*= corFactorL2L3Residual;
		tmprecpt *=corFactorL2L3Residual;
	}
	
	total_JEC = total_cor;
	
	return;     
}

void Read_JER(std::string mPtResoFile, std::string mPtSFFile, double tmprecpt, TLorentzVector pfjet4v, double Rho, edm::Handle<reco::GenJetCollection>  genjets, double dRcut, vector<double> &SFs)
{
 
	JME::JetResolution resolution;
	resolution = JME::JetResolution(mPtResoFile.c_str());
	JME::JetResolutionScaleFactor res_sf;
	res_sf = JME::JetResolutionScaleFactor(mPtSFFile.c_str());
	
	JME::JetParameters parameters_5 = {{JME::Binning::JetPt, tmprecpt}, {JME::Binning::JetEta, pfjet4v.Eta()}, {JME::Binning::Rho, Rho}};
	double rp = resolution.getResolution(parameters_5);
	double gaus_rp = gRandom->Gaus(0.,rp);
	double sf = res_sf.getScaleFactor(parameters_5, Variation::NOMINAL);
	double sf_up = res_sf.getScaleFactor(parameters_5, Variation::UP);
	double sf_dn = res_sf.getScaleFactor(parameters_5, Variation::DOWN);
	
	bool match = false;
	int match_gen = -1;
		
	for (unsigned get = 0; get<(genjets->size()); get++) {
		TLorentzVector genjet4v((*genjets)[get].px(),(*genjets)[get].py(),(*genjets)[get].pz(), (*genjets)[get].energy());
		if((delta2R(pfjet4v.Rapidity(),pfjet4v.Phi(),genjet4v.Rapidity(),genjet4v.Phi()) < (dRcut)) &&(fabs(tmprecpt-genjet4v.Pt())<(3*fabs(rp)*tmprecpt))){
			match = true;
			match_gen = get;
			break;
		}
	}
		
	if(match && (match_gen>=0)){
	  
		SFs.push_back((sf-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
		SFs.push_back((sf_up-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
		SFs.push_back((sf_dn-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
	  
	}else{
	  
		SFs.push_back(sqrt(max(0.,(sf*sf-1))) * gaus_rp);
		SFs.push_back(sqrt(max(0.,(sf_up*sf_up-1))) * gaus_rp);
		SFs.push_back(sqrt(max(0.,(sf_dn*sf_dn-1))) * gaus_rp);
	}
      	
}

float getEtaForEA(auto obj){
	float eta;
	if(abs(obj->pdgId())==11||abs(obj->pdgId())==22) { eta = obj->superCluster()->eta(); }     
	else { eta = obj->eta(); }
	return eta;    
}

std::unique_ptr<EffectiveAreas> ea_mu_miniiso_, ea_el_miniiso_;

void Read_MiniIsolation(auto obj, double Rho, vector<float> &isovalues)
{
	pat::PFIsolation iso = obj->miniPFIsolation();                                                                                                                                                                                                   
	float chg = iso.chargedHadronIso();                                                                                                                     
	float neu = iso.neutralHadronIso();                                                                                                                     
	float pho = iso.photonIso();                                                                                       
	                                                                                   
	float ea;
	if(abs(obj->pdgId())==13) { ea = ea_mu_miniiso_->getEffectiveArea(fabs(getEtaForEA(obj))); }
	else { ea = ea_el_miniiso_->getEffectiveArea(fabs(getEtaForEA(obj))); }  
	                                                                                    
	float R = 10.0/std::min(std::max(obj->pt(), 50.0),200.0);                                                                      
	ea *= std::pow(R / 0.3, 2);                                                                                                                  	
	float tot = (chg+std::max(0.0,neu+pho-(Rho)*ea));
	
	isovalues.push_back(tot);
	isovalues.push_back(chg);
	isovalues.push_back(neu);
	isovalues.push_back(pho);	
	
	for(unsigned ij=0; ij<isovalues.size(); ij++){
		isovalues[ij] *= 1./obj->pt();
	}
}

std::unique_ptr<EffectiveAreas> ea_el_pfiso_;

void Read_ElePFIsolation(auto obj, double Rho, vector<float> &isovalues)
{
	auto iso = obj->pfIsolationVariables();   
	auto  ea = ea_el_pfiso_->getEffectiveArea(fabs(getEtaForEA(obj)));                                                    
    float val = iso.sumChargedHadronPt + max(0., iso.sumNeutralHadronEt + iso.sumPhotonEt - (Rho)*ea); 
    float val04 = (obj->chargedHadronIso()+std::max(0.0,obj->neutralHadronIso()+obj->photonIso()-(Rho)*ea*16./9.));
    isovalues.push_back(val);
    isovalues.push_back(val04);
    
    for(unsigned ij=0; ij<isovalues.size(); ij++){
		isovalues[ij] *= 1./obj->pt();
	}    
}

float DistanceFromSV_Muon(TrackRef track, edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVertices){

	float dzmin = 1000;                                                                                                                              
    float dxymin = 1000;
    
    if(secondaryVertices.isValid()){                                                                                                                              
		for(unsigned int isv=0; isv<(secondaryVertices->size()); isv++){                                                                                            
		const auto &sv = (*secondaryVertices)[isv];                                                                                                               
     	  reco::TrackBase::Point svpoint(sv.vx(),sv.vy(),sv.vz());
		  float dztmp =fabs(track->dz(svpoint));
		  if(dztmp < dzmin){                                                                                                      
			dzmin = dztmp;                                                                                                        
			dxymin = track->dxy(svpoint);    
			}                                                                                                                                            
		}                                                                                                                                              
    }
    
    return  dxymin;	
}

float DistanceFromSV_Electron(GsfTrackRef track, edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVertices){

	float dzmin = 1000;                                                                                                                              
    float dxymin = 1000;
    
    if(secondaryVertices.isValid()){                                                                                                                              
		for(unsigned int isv=0; isv<(secondaryVertices->size()); isv++){                                                                                            
		const auto &sv = (*secondaryVertices)[isv];                                                                                                               
     	  reco::TrackBase::Point svpoint(sv.vx(),sv.vy(),sv.vz());
		  float dztmp =fabs(track->dz(svpoint));
		  if(dztmp < dzmin){                                                                                                      
			dzmin = dztmp;                                                                                                        
			dxymin = track->dxy(svpoint);    
			}                                                                                                                                            
		}                                                                                                                                              
    }
    
    return  dxymin;	
} 

//class declaration
//
class PNLepton : public edm::EDAnalyzer {
public:
  explicit PNLepton(const edm::ParameterSet&);
  ~PNLepton();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  void fillmetarray();
  // ----------member data ---------------------------
  int Nevt;
  bool isData;
  bool isMC;
  bool isFastSIM;
  int year;
  bool isUltraLegacy;
  bool store_electron_scalnsmear, store_electron_addvariabs;
  bool store_fatjet_constituents;
  bool read_btagSF;
  bool subtractLepton_fromAK4, subtractLepton_fromAK8;
  bool save_only_muons, save_only_electrons;
  
  uint nPDFsets;
  
  std::string theHLTTag;
  
  edm::EDGetTokenT<reco::VertexCollection> tok_primaryVertices_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> tok_sv;
  edm::EDGetTokenT<edm::View<pat::Jet>>tok_pfjetAK4s_;
  edm::EDGetTokenT<edm::View<pat::Muon>> tok_muons_;
  edm::EDGetTokenT<edm::View<pat::Electron>> tok_electrons_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tok_pfcands_;
  
  edm::EDGetTokenT<std::vector<reco::GenParticle>>tok_genparticles_;
  
  edm::EDGetTokenT<HepMCProduct> tok_HepMC ;
  edm::EDGetTokenT<GenEventInfoProduct> tok_wt_;
  
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileup_;
  edm::EDGetTokenT<double> tok_Rho_;
  
  //edm::EDGetTokenT <edm::ValueMap <bool> > tok_mvaPhoID_FallV2_WP90;
  //edm::EDGetTokenT <edm::ValueMap <bool> > tok_mvaPhoID_FallV2_WP80;
  edm::EDGetTokenT <edm::ValueMap <float> > tok_mvaPhoID_FallV2_raw;

  // object cuts //
  int iTag;
  int iTagMET;
  
  double minjPt;
  double maxEta;
  double minmuPt;
  double minePt;
  
  // Root file & tree //
  
  TFile* theFile;
  
  TTree* T1;
    
  unsigned ievt;
  
  static const int njetmx = 20; 
  static const int npartmx = 50; 
  static const int nconsmax = 100; 
  
  int irunold;
  int irun, ilumi, ifltr, ibrnch;
  
  int PV_npvsGood;
  
  double Generator_weight;
  double weights[njetmx];
  
  double Rho ;
  
 /* 
  int PFJetAK8_ncons[njetmxAK8];
  std::vector <std::vector <float> > PFJetAK8_cons_pt, PFJetAK8_cons_eta, PFJetAK8_cons_phi, PFJetAK8_cons_mass;
  std::vector <std::vector <int> > PFJetAK8_cons_pdgId;
  //float PFJetAK8_cons_pt[njetmxAK8][nconsmax], PFJetAK8_cons_eta[njetmxAK8][nconsmax], PFJetAK8_cons_phi[njetmxAK8][nconsmax], PFJetAK8_cons_mass[njetmxAK8][nconsmax];
*/  

  int nPFJetAK4;
  float PFJetAK4_pt[njetmx], PFJetAK4_eta[njetmx], PFJetAK4_y[njetmx], PFJetAK4_phi[njetmx], PFJetAK4_mass[njetmx];
  float PFJetAK4_btag_DeepCSV[njetmx], PFJetAK4_btag_DeepFlav[njetmx]; 
  bool PFJetAK4_jetID[njetmx], PFJetAK4_jetID_tightlepveto[njetmx];
  float PFJetAK4_JEC[njetmx];
  int PFJetAK4_hadronflav[njetmx], PFJetAK4_partonflav[njetmx];
  int PFJetAK4_Ncons[njetmx];
  float PFJetAK4_qgl[njetmx], PFJetAK4_PUID[njetmx];
  
  // Lepton variables //
  
  int nLepton;
  
  // scalar of all variables //
  
  // Common variables between muon & electrons //
  
  //basic kinematic variables //
  float lepton_pt, lepton_eta, lepton_phi, lepton_mass, lepton_p;
  int lepton_pdgId, lepton_charge, lepton_tightcharge;
  //distance from PV & SV//
  float lepton_dxy, lepton_dxyError, lepton_dz, lepton_dzError, lepton_ip3d, lepton_sip3d, lepton_dxy_sv;
  //Track information //
  float lepton_chi2,  Lepton_valfrac;
  int lepton_ndof, lepton_trkKink, lepton_hit, lepton_pixhit, lepton_nTrackerLayers, lepton_lostHits;
  //Energy information //
  float lepton_e_ECAL, lepton_e_HCAL, lepton_hoe;
  //Mini-isolation variables //
  float lepton_minchiso, lepton_minnhiso, lepton_minphiso, lepton_minisoall; 
  //PF-isolation variables //
  float lepton_pfRelIso03_drcor;
  float lepton_pfRelIso03_ChargedHadron, lepton_pfRelIso03_NeutralHadron, lepton_pfRelIso03_Photon, lepton_pfRelIso03_PileUp, lepton_tkRelIso;
  float lepton_pfRelIso04_drcor;
  float lepton_pfRelIso04_ChargedHadron, lepton_pfRelIso04_NeutralHadron, lepton_pfRelIso04_Photon, lepton_pfRelIso04_PileUp;
  
  // Electron-specific variables //
  
  // ID booleans //
  bool lepton_mvaid_Fallv2WP90_noIso, lepton_mvaid_Fallv2WP90, lepton_mvaid_Fallv2WP80_noIso, lepton_mvaid_Fallv2WP80;
  // isolation //
  float lepton_pfRelIso03_eacor, lepton_pfRelIso04_eacor;
  // energy cluster //
  float lepton_dr03EcalRecHitSumEt_Rel, lepton_dr03HcalDepth1TowerSumEt_Rel, lepton_dr03HcalDepth2TowerSumEt_Rel, lepton_dr03TkSumPt_Rel, lepton_dr03TkSumPtHEEP_Rel;
  float lepton_eoverp, lepton_etain, lepton_phiin, lepton_fbrem; 
  // super-clyster// 
  float lepton_supcl_eta, lepton_supcl_phi, lepton_supcl_energy; 
  float lepton_r9full;
  float lepton_supcl_etaWidth;
  float lepton_supcl_phiWidth;
  float lepton_hcaloverecal;
  int  lepton_closeTrackNLayers;
  float lepton_closeTrackNormChi2;
  float lepton_e1x5bye5x5;
  float lepton_normchi2;
  float lepton_convtxprob;
  float lepton_ecloverpout;
  float lepton_deltaetacltrkcalo;
  float lepton_supcl_preshvsrawe;
  float lepton_sigmaietaieta, lepton_sigmaiphiiphi;
  float lepton_eInvMinusPInv;
  bool lepton_convVeto; 
  int lepton_genPartFlav;
  //track//
  float lepton_posmatch, lepton_segmentComp;
  int lepton_nStations;
  // ID booleans //
  bool lepton_isPFCand, lepton_isGlobal, lepton_isTracker;
  bool lepton_isGoodGlobal, lepton_isTight, lepton_isHighPt, lepton_isHighPttrk, lepton_isMedium, lepton_isMedPr, lepton_isLoose;
  
  float lepton_jetPtRelv2, lepton_jetRelIso, lepton_jetbtag;
  
  bool label_Muon_Prompt, label_Muon_fromTau, label_Muon_fromHadron, label_Muon_fromPhoton, label_Muon_unknown;
  bool label_Electron_Prompt, label_Electron_fromTau, label_Electron_fromHadron, label_Electron_fromPhoton, label_Electron_unknown;
  
  // scalar end
  
  // PF candidates //
  
  int nPFCand;
  
  int nChargePFCand;
  vector<float> ChargePFCand_pt_rel;
  vector<float> ChargePFCand_deta;
  vector<float> ChargePFCand_dphi;
  vector<float> ChargePFCand_dphiAtVtx;
  vector<float> ChargePFCand_puppiWeight;
  vector<float> ChargePFCand_puppiWeightNoLep;
  vector<float> ChargePFCand_caloFraction;
  vector<float> ChargePFCand_hcalFraction;
  vector<int> ChargePFCand_pdgId;
  
  vector<float> ChargePFCand_dz;
  vector<float> ChargePFCand_dzError;
  vector<float> ChargePFCand_dxy;
  vector<float> ChargePFCand_dxyError;
  vector<float> ChargePFCand_trkChi2;
  vector<float> ChargePFCand_vertexChi2;
  vector<int> ChargePFCand_charge;
  vector<int> ChargePFCand_lostInnerHits;
  vector<float> ChargePFCand_pvAssocQuality;
  vector<float> ChargePFCand_trkQuality;
  vector<int> ChargePFCand_nTrackerLayers;
  vector<int> ChargePFCand_pixelhits;
  vector<float> ChargePFCand_status;
  vector<float> ChargePFCand_time;
  vector<bool> ChargePFCand_trackHighPurity;
  
  int nNeutralPFCand;
  vector<float> NeutralPFCand_pt_rel;
  vector<float> NeutralPFCand_deta;
  vector<float> NeutralPFCand_dphi;
  vector<float> NeutralPFCand_dphiAtVtx;
  vector<float> NeutralPFCand_puppiWeight;
  vector<float> NeutralPFCand_puppiWeightNoLep;
  vector<float> NeutralPFCand_caloFraction;
  vector<float> NeutralPFCand_hcalFraction;
  vector<int> NeutralPFCand_pdgId;
 
  // secondary vertices //
  int nSV;
  vector<float> SV_x;
  vector<float> SV_y;
  vector<float> SV_z;
  vector<float> SV_pt_rel;
  vector<float> SV_mass;
  vector<float> SV_chi2;
  vector<int> SV_ndof;
  vector<float> SV_dxy;
  vector<float> SV_dxySig;
  vector<float> SV_dlen;
  vector<float> SV_dlenSig;
  vector<float> SV_pAngle;
  
  // Collision Info //
 
  int Pileup_nPU;
  int Pileup_nTrueInt;
  
  HLTPrescaleProvider hltPrescaleProvider_;
    
  // ---- Jet Corrector Parameter ---- //
  
  JetCorrectorParameters *L1FastAK4, *L2RelativeAK4, *L3AbsoluteAK4, *L2L3ResidualAK4;
  vector<JetCorrectorParameters> vecL1FastAK4, vecL2RelativeAK4, vecL3AbsoluteAK4, vecL2L3ResidualAK4;
  FactorizedJetCorrector *jecL1FastAK4, *jecL2RelativeAK4, *jecL3AbsoluteAK4, *jecL2L3ResidualAK4;
  
  // ---- Jet Corrector Parameter End---- //

  // BTagCalibration Begin //

  BTagCalibration calib_deepcsv, calib_deepflav;
  BTagCalibrationReader reader_deepcsv, reader_deepflav;
  
  // BTagCalibration End //

  // ---- Jet Resolution Parameter ---- //
  
  std::string mJECL1FastFileAK4, mJECL2RelativeFileAK4, mJECL3AbsoluteFileAK4, mJECL2L3ResidualFileAK4;
  std::string mPtResoFileAK4, mPtSFFileAK4;
  
  // ---- Jet Resolution Parameter End---- //
  
  // ---- B tagging scale factor files --- //
  
  std::string mBtagSF_DeepCSV;
  std::string mBtagSF_DeepFlav;
  
  // ---- B tagging scale factor files End --- //
  
  // ---- Rochester correction files --- //
  
  std::string mRochcorFolder;
  
  // Electron MVA ID //
  
  std::string melectronID_isowp90, melectronID_noisowp90;
  std::string melectronID_isowp80, melectronID_noisowp80;
  
  // Photon MVA ID //
  
  std::string mPhoID_FallV2_WP90, mPhoID_FallV2_WP80;
  std::string mPhoID_SpringV1_WP90, mPhoID_SpringV1_WP80;
  
  // Rochester correction for muons//
  
  RoccoR roch_cor; 
  
  // GEN particles //
  vector<GenParton> GenPartons;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

PNLepton::PNLepton(const edm::ParameterSet& pset):
  hltPrescaleProvider_(pset, consumesCollector(), *this)  
{
  //now do what ever initialization is needed
  
  edm::Service<TFileService> fs;
  
  isData    = pset.getUntrackedParameter<bool>("Data",false);
  isMC      = pset.getUntrackedParameter<bool>("MonteCarlo", false);
  isFastSIM      = pset.getUntrackedParameter<bool>("FastSIM", false);
  year		= pset.getUntrackedParameter<int>("YEAR", 2018);
  isUltraLegacy = pset.getUntrackedParameter<bool>("UltraLegacy", false);
  theHLTTag = pset.getUntrackedParameter<string>("HLTTag", "HLT");
  store_electron_scalnsmear = pset.getUntrackedParameter<bool>("store_electron_scalnsmear", false);
  store_electron_addvariabs = pset.getUntrackedParameter<bool>("store_electron_addvariabs", false);
  read_btagSF = pset.getUntrackedParameter<bool>("Read_btagging_SF", false);
  save_only_muons = pset.getUntrackedParameter<bool>("Save_only_Muons", false);
  save_only_electrons = pset.getUntrackedParameter<bool>("Save_only_Electrons", false);

  minjPt = pset.getUntrackedParameter<double>("minjPt",25.);
  minmuPt = pset.getUntrackedParameter<double>("minmuPt",10.);
  minePt = pset.getUntrackedParameter<double>("minePt",10.);
 
  maxEta = pset.getUntrackedParameter<double>("maxEta",3.);
 
  ea_mu_miniiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_MuonMiniIso")).fullPath()));
  ea_el_miniiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_EleMiniIso")).fullPath()));
  ea_el_pfiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_ElePFIso")).fullPath()));
 
  tok_primaryVertices_ =consumes<reco::VertexCollection>( pset.getParameter<edm::InputTag>("PrimaryVertices"));
  tok_sv =consumes<reco::VertexCompositePtrCandidateCollection>( pset.getParameter<edm::InputTag>("SecondaryVertices"));
  
  tok_Rho_ = consumes<double>(pset.getParameter<edm::InputTag>("PFRho"));
    
  tok_muons_ = consumes<edm::View<pat::Muon>> ( pset.getParameter<edm::InputTag>("Muons"));
  tok_electrons_ = consumes<edm::View<pat::Electron>> ( pset.getParameter<edm::InputTag>("Electrons"));
  tok_pfjetAK4s_= consumes<edm::View<pat::Jet>>( pset.getParameter<edm::InputTag>("PFJetsAK4"));
  tok_pfcands_= consumes<edm::View<pat::PackedCandidate>>( pset.getParameter<edm::InputTag>("pfCands"));
  
  
  if(isMC){
	  
    tok_genparticles_ = consumes<std::vector<reco::GenParticle>>( pset.getParameter<edm::InputTag>("GenParticles"));
    
    tok_HepMC = consumes<HepMCProduct>(pset.getParameter<edm::InputTag>("Generator"));
    tok_wt_ = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("Generator")) ;
    pileup_ = consumes<std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("slimmedAddPileupInfo"));
    
  }
  
  melectronID_isowp90       = pset.getParameter<std::string>("electronID_isowp90");
  melectronID_noisowp90     = pset.getParameter<std::string>("electronID_noisowp90");
  melectronID_isowp80       = pset.getParameter<std::string>("electronID_isowp80");
  melectronID_noisowp80     = pset.getParameter<std::string>("electronID_noisowp80");
  
  mJECL1FastFileAK4         = pset.getParameter<std::string>("jecL1FastFileAK4");
  mJECL2RelativeFileAK4     = pset.getParameter<std::string>("jecL2RelativeFileAK4");
  mJECL3AbsoluteFileAK4     = pset.getParameter<std::string>("jecL3AbsoluteFileAK4");
  mJECL2L3ResidualFileAK4   = pset.getParameter<std::string> ("jecL2L3ResidualFileAK4");
  
  mPtResoFileAK4  = pset.getParameter<std::string>("PtResoFileAK4");
  mPtSFFileAK4  = pset.getParameter<std::string>("PtSFFileAK4");
    
  mBtagSF_DeepCSV = pset.getParameter<std::string>("BtagSFFile_DeepCSV");
  mBtagSF_DeepFlav = pset.getParameter<std::string>("BtagSFFile_DeepFlav");
  
  mRochcorFolder = pset.getParameter<std::string>("RochcorFolder");
  
  //T1 = new TTree("Events", "ParticleNet");
  T1 = fs->make<TTree>("tree","ParticleNet") ;
 
  T1->Branch("irun", &irun, "irun/I");  
  T1->Branch("ilumi", &ilumi, "ilumi/I");  
  T1->Branch("ievt", &ievt, "ievt/i");
  
  // primary vertices //
  
  T1->Branch("PV_npvsGood", &PV_npvsGood, "PV_npvsGood/I");
 
  // Lepton info //
  
  //T1->Branch("nLepton",&nLepton,"nLepton/I");
  
  // Labels //
  
  T1->Branch("label_Muon_Prompt", &label_Muon_Prompt, "label_Muon_Prompt/O");
  T1->Branch("label_Muon_fromHadron", &label_Muon_fromHadron, "label_Muon_fromHadron/O");
  T1->Branch("label_Muon_fromTau", &label_Muon_fromTau, "label_Muon_fromTau/O");
  T1->Branch("label_Muon_fromPhoton", &label_Muon_fromPhoton, "label_Muon_fromPhoton/O");
  T1->Branch("label_Muon_unknown", &label_Muon_unknown, "label_Muon_unknown/O");
  T1->Branch("label_Electron_Prompt", &label_Electron_Prompt, "label_Electron_Prompt/O");
  T1->Branch("label_Electron_fromTau", &label_Electron_fromTau, "label_Electron_fromTau/O");
  T1->Branch("label_Electron_fromHadron", &label_Electron_fromHadron, "label_Electron_fromHadron/O");
  T1->Branch("label_Electron_fromPhoton", &label_Electron_fromPhoton, "label_Electron_fromPhoton/O");
  T1->Branch("label_Electron_unknown", &label_Electron_unknown, "label_Electron_unknown/O");
  
  // common kinematic variables //
  
  T1->Branch("lepton_pt", &lepton_pt,"lepton_pt/F");
  T1->Branch("lepton_p", &lepton_p,"lepton_p/F");
  T1->Branch("lepton_eta", &lepton_eta,"lepton_eta/F");
  T1->Branch("lepton_phi", &lepton_phi,"lepton_phi/F");
  T1->Branch("lepton_mass", &lepton_mass,"lepton_mass/F");
  
  T1->Branch("lepton_charge", &lepton_charge,"lepton_charge/I");
  T1->Branch("lepton_tightcharge", &lepton_tightcharge,"lepton_tightcharge/I");
  T1->Branch("lepton_pdgId", &lepton_pdgId,"lepton_pdgId/I");
  
  // common distance from PV & SV //
  
  T1->Branch("lepton_dxy", &lepton_dxy,"lepton_dxy/F");
  T1->Branch("lepton_dz", &lepton_dz,"lepton_dz/F");
  T1->Branch("lepton_dxyError", &lepton_dxyError,"lepton_dxyError/F");
  T1->Branch("lepton_dzError", &lepton_dzError,"lepton_dzError/F");
  T1->Branch("lepton_ip3d", &lepton_ip3d,"lepton_ip3d/F");
  T1->Branch("lepton_sip3d", &lepton_sip3d,"lepton_sip3d/F");
  T1->Branch("lepton_dxy_sv", &lepton_dxy_sv,"lepton_dxy_sv/F");
  
  // Generator-matching //
  
  T1->Branch("lepton_genPartFlav", &lepton_genPartFlav,"lepton_genPartFlav/I");
  
  // common Track information //
  
  T1->Branch("lepton_chi2", &lepton_chi2,"lepton_chi2/F");
  T1->Branch("lepton_ndof", &lepton_ndof,"lepton_ndof/I");
  T1->Branch("lepton_trkKink", &lepton_trkKink,"lepton_trkKink/F");
  T1->Branch("lepton_hit", &lepton_hit,"lepton_hit/I");
  T1->Branch("lepton_pixhit", &lepton_pixhit,"lepton_pixhit/I");
  T1->Branch("lepton_nTrackerLayers", &lepton_nTrackerLayers,"lepton_nTrackerLayers/I"); 
  T1->Branch("lepton_lostHits", &lepton_lostHits,"lepton_lostHits/I");
  
  // common calorimeter energy info //
  
  T1->Branch("lepton_e_ECAL", &lepton_e_ECAL,"lepton_e_ECAL/F");
  T1->Branch("lepton_e_HCAL", &lepton_e_HCAL,"lepton_e_HCAL/F");
  T1->Branch("lepton_hoe", &lepton_hoe,"lepton_hoe/F");
  
  // common mini-isolation variables //
  
  T1->Branch("lepton_minisoch", &lepton_minchiso, "lepton_minchiso/F");
  T1->Branch("lepton_minisonh", &lepton_minnhiso, "lepton_minnhiso/F");
  T1->Branch("lepton_minisoph", &lepton_minphiso, "lepton_minphiso/F");
  T1->Branch("lepton_minisoall", &lepton_minisoall, "lepton_minisoall/F");
  
  // common isolation variables //
  
  T1->Branch("lepton_pfRelIso03_drcor", &lepton_pfRelIso03_drcor,"lepton_pfRelIso03_drcor/F");
  T1->Branch("lepton_pfRelIso03_ChargedHadron", &lepton_pfRelIso03_ChargedHadron,"lepton_pfRelIso03_ChargedHadron/F");
  T1->Branch("lepton_pfRelIso03_NeutralHadron", &lepton_pfRelIso03_NeutralHadron,"lepton_pfRelIso03_NeutralHadron/F");
  T1->Branch("lepton_pfRelIso03_Photon", &lepton_pfRelIso03_Photon,"lepton_pfRelIso03_Photon/F");
  T1->Branch("lepton_pfRelIso03_PileUp", &lepton_pfRelIso03_PileUp,"lepton_pfRelIso03_PileUp/F");
  T1->Branch("lepton_tkRelIso", &lepton_tkRelIso,"lepton_tkRelIso/F");
  
  // Muon-specfic variables //
  
  T1->Branch("lepton_nStations", &lepton_nStations,"lepton_nStations/I");
  T1->Branch("lepton_segmentComp", &lepton_segmentComp,"lepton_segmentComp/F");
  
  // Muon ID booleans //
  
  T1->Branch("lepton_isPFCand", &lepton_isPFCand,"lepton_isPFCand/O");
  T1->Branch("lepton_isGlobal", &lepton_isGlobal,"lepton_isGlobal/O");
  T1->Branch("lepton_isTracker", &lepton_isTracker,"lepton_isTracker/O");
  T1->Branch("lepton_isLoose", &lepton_isLoose,"lepton_isLoose/O");
  T1->Branch("lepton_isGoodGlobal", &lepton_isGoodGlobal,"lepton_isGoodGlobal/O");
  T1->Branch("lepton_isMedium", &lepton_isMedium,"lepton_isMedium/O");
  T1->Branch("lepton_isMedPr", &lepton_isMedPr,"lepton_isMedPr/O");
  T1->Branch("lepton_isTight", &lepton_isTight,"lepton_isTight/O");
  T1->Branch("lepton_isHighPt", &lepton_isHighPt,"lepton_isHighPt/O"); 
  T1->Branch("lepton_isHighPttrk", &lepton_isHighPttrk,"lepton_isHighPttrk/O");
  
  T1->Branch("lepton_posmatch", &lepton_posmatch,"lepton_posmatch/F");
  
  // Electron ID booleans //
  
  T1->Branch("lepton_mvaid_Fallv2WP90", &lepton_mvaid_Fallv2WP90,"lepton_mvaid_Fallv2WP90/O");
  T1->Branch("lepton_mvaid_Fallv2WP90_noIso", &lepton_mvaid_Fallv2WP90_noIso,"lepton_mvaid_Fallv2WP90_noIso/O");
  T1->Branch("lepton_mvaid_Fallv2WP80", &lepton_mvaid_Fallv2WP80,"lepton_mvaid_Fallv2WP80/O");
  T1->Branch("lepton_mvaid_Fallv2WP80_noIso", &lepton_mvaid_Fallv2WP80_noIso,"lepton_mvaid_Fallv2WP80_noIso/O");
  
  // Electron energy cluster shapes //
  
  T1->Branch("lepton_eoverp", &lepton_eoverp,"lepton_eoverp/F");
  
  // Electron-specific isolation variables //
 
  T1->Branch("lepton_pfRelIso03_eacor", &lepton_pfRelIso03_eacor,"lepton_pfRelIso03_eacor/F");
  T1->Branch("lepton_pfReliso04_eacor", &lepton_pfRelIso04_eacor,"lepton_pfRelIso04_eacor/F");
  
  T1->Branch("lepton_dr03EcalRecHitSumEt_Rel", &lepton_dr03EcalRecHitSumEt_Rel,"lepton_dr03EcalRecHitSumEt_Rel/F");
  T1->Branch("lepton_dr03HcalDepth1TowerSumEt_Rel", &lepton_dr03HcalDepth1TowerSumEt_Rel,"lepton_dr03HcalDepth1TowerSumEt_Rel/F");
  T1->Branch("lepton_dr03HcalDepth2TowerSumEt_Rel", &lepton_dr03HcalDepth2TowerSumEt_Rel,"lepton_dr03HcalDepth2TowerSumEt_Rel/F");
  T1->Branch("lepton_dr03TkSumPt_Rel", &lepton_dr03TkSumPt_Rel,"lepton_dr03TkSumPt_Rel/F");
  T1->Branch("lepton_dr03TkSumPtHEEP_Rel", &lepton_dr03TkSumPtHEEP_Rel,"lepton_dr03TkSumPtHEEP_Rel/F");
  
  // Electron supercluster information //
  
  T1->Branch("lepton_eInvMinusPInv", &lepton_eInvMinusPInv,"lepton_eInvMinusPInv/F");
  T1->Branch("lepton_supcl_eta", &lepton_supcl_eta,"lepton_supcl_eta/F");
  T1->Branch("lepton_supcl_phi", &lepton_supcl_phi,"lepton_supcl_phi/F");
  T1->Branch("lepton_supcl_energy", &lepton_supcl_energy,"lepton_supcl_energy/F");
  T1->Branch("lepton_sigmaietaieta", &lepton_sigmaietaieta, "lepton_sigmaietaieta/F");
  T1->Branch("lepton_sigmaiphiiphi", &lepton_sigmaiphiiphi, "lepton_sigmaiphiiphi/F");
  T1->Branch("lepton_r9full", &lepton_r9full, "lepton_r9full/F");
  
  T1->Branch("lepton_hcaloverecal", &lepton_hcaloverecal, "lepton_hcaloverecal/F");
  T1->Branch("lepton_ecloverpout", &lepton_ecloverpout, "lepton_ecloverpout/F"); 
  T1->Branch("lepton_convVeto", &lepton_convVeto, "lepton_convVeto/O");

  T1->Branch("lepton_etain", &lepton_etain,"lepton_etain/F");
  T1->Branch("lepton_phiin", &lepton_phiin,"lepton_phiin/F");
  T1->Branch("lepton_fbrem", &lepton_fbrem,"lepton_fbrem/F");
  T1->Branch("lepton_supcl_etaWidth", &lepton_supcl_etaWidth, "lepton_supcl_etaWidth/F");
  T1->Branch("lepton_supcl_phiWidth", &lepton_supcl_phiWidth, "lepton_supcl_phiWidth/F");
  
  T1->Branch("lepton_e1x5bye5x5", &lepton_e1x5bye5x5, "lepton_e1x5bye5x5/F");
  T1->Branch("lepton_convtxprob", &lepton_convtxprob, "lepton_convtxprob/F");
  T1->Branch("lepton_deltaetacltrkcalo", &lepton_deltaetacltrkcalo, "lepton_deltaetacltrkcalo/F");
  T1->Branch("lepton_supcl_preshvsrawe", &lepton_supcl_preshvsrawe, "lepton_supcl_preshvsrawe/F");
  
  T1->Branch("lepton_closeTrackNLayers", &lepton_closeTrackNLayers, "lepton_closeTrackNLayers/I");
  T1->Branch("lepton_closeTrackNormChi2", &lepton_closeTrackNormChi2, "lepton_closeTrackNormChi2/F");
  
  // Nearest jet features //
  
  T1->Branch("lepton_jetPtRelv2", &lepton_jetPtRelv2, "lepton_closeTrackNLayers/F");
  T1->Branch("lepton_jetRelIso", &lepton_jetRelIso, "lepton_jetRelIso/F");
  T1->Branch("lepton_jetbtag", &lepton_jetbtag, "lepton_jetbtag/F");
  
  // PF candidates //
  
  T1->Branch("nPFCand",&nPFCand,"nPFCand/I");
  
  T1->Branch("ChargePFCand_pt_rel",&ChargePFCand_pt_rel);
  T1->Branch("ChargePFCand_deta",&ChargePFCand_deta);
  T1->Branch("ChargePFCand_dphi",&ChargePFCand_dphi);
  T1->Branch("ChargePFCand_dphiAtVtx",&ChargePFCand_dphiAtVtx);
  //T1->Branch("ChargePFCand_mass",&ChargePFCand_mass);
  T1->Branch("ChargePFCand_pdgId",&ChargePFCand_pdgId);
  T1->Branch("ChargePFCand_caloFraction",&ChargePFCand_caloFraction);
  T1->Branch("ChargePFCand_hcalFraction",&ChargePFCand_hcalFraction);
  T1->Branch("ChargePFCand_puppiWeight",&ChargePFCand_puppiWeight);
  T1->Branch("ChargePFCand_puppiWeightNoLep",&ChargePFCand_puppiWeightNoLep);
  
  T1->Branch("ChargePFCand_dz",&ChargePFCand_dz);
  T1->Branch("ChargePFCand_dzError",&ChargePFCand_dzError);
  T1->Branch("ChargePFCand_dxy",&ChargePFCand_dxy);
  T1->Branch("ChargePFCand_dxyError",&ChargePFCand_dxyError);
  T1->Branch("ChargePFCand_vertexChi2",&ChargePFCand_vertexChi2);
  T1->Branch("ChargePFCand_time",&ChargePFCand_time);
  T1->Branch("ChargePFCand_charge",&ChargePFCand_charge);
  T1->Branch("ChargePFCand_lostInnerHits",&ChargePFCand_lostInnerHits);
  T1->Branch("ChargePFCand_pvAssocQuality",&ChargePFCand_pvAssocQuality);
  T1->Branch("ChargePFCand_status",&ChargePFCand_status);
  T1->Branch("ChargePFCand_pixelhits",&ChargePFCand_pixelhits);
  T1->Branch("ChargePFCand_nTrackerLayers",&ChargePFCand_nTrackerLayers);
  T1->Branch("ChargePFCand_trkChi2",&ChargePFCand_trkChi2);
  T1->Branch("ChargePFCand_trackHighPurity",&ChargePFCand_trackHighPurity);
  
  T1->Branch("NeutralPFCand_pt_rel",&NeutralPFCand_pt_rel);
  T1->Branch("NeutralPFCand_deta",&NeutralPFCand_deta);
  T1->Branch("NeutralPFCand_dphi",&NeutralPFCand_dphi);
  T1->Branch("NeutralPFCand_dphiAtVtx",&NeutralPFCand_dphiAtVtx);
  //T1->Branch("NeutralPFCand_mass",&NeutralPFCand_mass);
  T1->Branch("NeutralPFCand_pdgId",&NeutralPFCand_pdgId);
  T1->Branch("NeutralPFCand_caloFraction",&NeutralPFCand_caloFraction);
  T1->Branch("NeutralPFCand_hcalFraction",&NeutralPFCand_hcalFraction);
  T1->Branch("NeutralPFCand_puppiWeight",&NeutralPFCand_puppiWeight);
  T1->Branch("NeutralPFCand_puppiWeightNoLep",&NeutralPFCand_puppiWeightNoLep);
  
  // Secondary vertices //
  
  T1->Branch("nSV",&nSV,"nSV/I");
  T1->Branch("SV_x",&SV_x);
  T1->Branch("SV_y",&SV_y);
  T1->Branch("SV_z",&SV_z);
  T1->Branch("SV_chi2",&SV_chi2);
  T1->Branch("SV_ndof",&SV_ndof);
  T1->Branch("SV_pt",&SV_pt_rel);
  T1->Branch("SV_mass",&SV_mass);
  T1->Branch("SV_dxy",&SV_dxy);
  T1->Branch("SV_dxySig",&SV_dxySig);
  T1->Branch("SV_dlen",&SV_dlen);
  T1->Branch("SV_dlenSig",&SV_dlenSig);
  T1->Branch("SV_pAngle",&SV_pAngle);
  
  // MC Info //
  
  if(isMC){
	 
	T1->Branch("Generator_weight", &Generator_weight, "Generator_weight/D") ;
	T1->Branch("Pileup_nPU",&Pileup_nPU,"Pileup_nPU/I");
	T1->Branch("Pileup_nTrueInt",&Pileup_nTrueInt,"Pileup_nTrueInt/I");
  
  } //isMC
 
  
  Nevt=0;
}


PNLepton::~PNLepton()
{
 
  // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


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
			//if(!((*genparticles)[ig].isHardProcess())) continue;
	  
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


// ------------ method called once each job just before starting event loop  ------------
void 
PNLepton::beginJob()
{
  
  Nevt = 0;
  
  ////JEC /////
  
  L1FastAK4       = new JetCorrectorParameters(mJECL1FastFileAK4.c_str());
  L2RelativeAK4   = new JetCorrectorParameters(mJECL2RelativeFileAK4.c_str());
  L3AbsoluteAK4   = new JetCorrectorParameters(mJECL3AbsoluteFileAK4.c_str());
  L2L3ResidualAK4 = new JetCorrectorParameters(mJECL2L3ResidualFileAK4.c_str());
  
  vecL1FastAK4.push_back(*L1FastAK4);
  vecL2RelativeAK4.push_back(*L2RelativeAK4);
  vecL3AbsoluteAK4.push_back(*L3AbsoluteAK4);
  vecL2L3ResidualAK4.push_back(*L2L3ResidualAK4);
  
  jecL1FastAK4       = new FactorizedJetCorrector(vecL1FastAK4);
  jecL2RelativeAK4   = new FactorizedJetCorrector(vecL2RelativeAK4);
  jecL3AbsoluteAK4   = new FactorizedJetCorrector(vecL3AbsoluteAK4);
  jecL2L3ResidualAK4 = new FactorizedJetCorrector(vecL2L3ResidualAK4);
  
  if(read_btagSF){
	calib_deepcsv = BTagCalibration("DeepCSV", mBtagSF_DeepCSV.c_str());
	reader_deepcsv = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}); 
	reader_deepcsv.load(calib_deepcsv, BTagEntry::FLAV_B, "comb");
	reader_deepcsv.load(calib_deepcsv, BTagEntry::FLAV_C, "comb");
	reader_deepcsv.load(calib_deepcsv, BTagEntry::FLAV_UDSG, "incl");
  
	calib_deepflav = BTagCalibration("DeepJet", mBtagSF_DeepFlav.c_str());
	reader_deepflav = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}); 
	reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
	reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
	reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");
  }
  
  if(isUltraLegacy)
  {
	if(year==2018){
		roch_cor.init((mRochcorFolder+"RoccoR2018UL.txt").c_str()); 
	}
	if(year==2017){
		roch_cor.init((mRochcorFolder+"RoccoR2017UL.txt").c_str()); 
	}
	if(year==2016){
		roch_cor.init((mRochcorFolder+"RoccoR2016aUL.txt").c_str()); 
	}
  }
  else{
		roch_cor.init((mRochcorFolder+"RoccoR2017.txt").c_str()); 
	  }
  
  //**Important**//
  //For precision top physics, change "comb" to "mujets" in BTagCalibrationReader above //
  //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18#Additional_information
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PNLepton::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
void 
PNLepton::beginRun(edm::Run const& iRun, edm::EventSetup const& pset)
{	
  bool changed(true);
  if(!isFastSIM){
	hltPrescaleProvider_.init(iRun,pset,theHLTTag,changed);
	HLTConfigProvider const&  hltConfig_ = hltPrescaleProvider_.hltConfigProvider();
  }
  //hltConfig_.dump("Triggers");
  //hltConfig_.dump("PrescaleTable");
}

// ------------ method called when ending the processing of a run  ------------
void 
PNLepton::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PNLepton::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PNLepton::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PNLepton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PNLepton);
