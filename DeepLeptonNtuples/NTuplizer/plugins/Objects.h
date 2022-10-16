#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>

using namespace std;

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
		float dxySig;
		float dzSig;
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
		int genPartIdx;
		int genPartFlav;
		
		// Electron-specific variables //
  
		// ID booleans //
		bool mvaid_Fallv2WP90_noIso;
		bool mvaid_Fallv2WP90;
		bool mvaid_Fallv2WP80_noIso;
		bool mvaid_Fallv2WP80;
		// ID floats //
		float mvaFall17V1Iso;
		float mvaFall17V1noIso;
		float mvaFall17V2Iso;
		float mvaFall17V2noIso;
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
		float dEtaInSeed;
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
  
		pat::Muon MuonPtr;
		pat::Electron ElectronPtr;
};

void Initialize(Lepton &lepton){
	
	lepton.pt = 0 ;
	lepton.eta = -100 ;
	lepton.phi = -100 ;
	lepton.mass = 0 ;
	lepton.p = 0 ;
	
	lepton.pdgId = 0 ;
	lepton.charge = 0 ;
	lepton.tightcharge = 0 ;
		
	lepton.dxy = 0;
	lepton.dz = 0;
	lepton.dxyError = 0;
	lepton.dzError = 0;
	lepton.dxySig = 0;
	lepton.dzSig = 0;
	lepton.ip3d = 0;
	lepton.sip3d = 0;
	lepton.dxy_sv = 0;
		
	lepton.genPartFlav = -1;
	lepton.genPartIdx = 0;
  
  	lepton.chi2 = 999;
	lepton.ndof = 0;
	lepton.trkKink = 0;
	lepton.hit = 0;
	lepton.pixhit = 0;
	lepton.nTrackerLayers = 0;
	lepton.lostHits = -1;
    
	lepton.e_ECAL = 0;
	lepton.e_HCAL = 0;
	lepton.hoe = 0;
		  
	lepton.minchiso = 0;
	lepton.minnhiso = 0;
	lepton.minphiso = 0;
	lepton.minisoall = 0;
		
	lepton.pfRelIso03_drcor = 0;
	lepton.pfRelIso03_ChargedHadron = 0;
	lepton.pfRelIso03_NeutralHadron = 0;
	lepton.pfRelIso03_Photon = 0;
	lepton.pfRelIso03_PileUp = 0;
	lepton.tkRelIso = 0;
   
	lepton.mvaid_Fallv2WP90 = false;
	lepton.mvaid_Fallv2WP90_noIso = false;
	lepton.mvaid_Fallv2WP80 = false;
	lepton.mvaid_Fallv2WP80_noIso = false;
	
	lepton.mvaFall17V1Iso = -1;
	lepton.mvaFall17V1noIso = -1;
	lepton.mvaFall17V2Iso = -1;
	lepton.mvaFall17V2noIso = -1;
   
	lepton.eoverp = 0;
 
	lepton.pfRelIso03_eacor = 0;
	lepton.pfRelIso04_eacor = 0;
	lepton.dr03EcalRecHitSumEt_Rel = 0;
	lepton.dr03HcalDepth1TowerSumEt_Rel = 0;
	lepton.dr03HcalDepth2TowerSumEt_Rel = 0;
	lepton.dr03TkSumPt_Rel = 0;
	lepton.dr03TkSumPtHEEP_Rel = 0;
    
	lepton.eInvMinusPInv = 0;
	lepton.supcl_eta = -100;
	lepton.supcl_phi = -100;
	lepton.supcl_energy = 999;
	lepton.sigmaietaieta = 999;
	lepton.sigmaiphiiphi = 999;
	lepton.r9full = 999;
	lepton.hcaloverecal = 0;
	lepton.ecloverpout = 0;
	lepton.convVeto = 0;
	lepton.etain = 999;
	lepton.dEtaInSeed = 999;
	lepton.phiin = 999;
	lepton.fbrem = 0;
	lepton.supcl_etaWidth = 999;
	lepton.supcl_phiWidth = 999;
 
	lepton.e1x5bye5x5 = 0;
	lepton.convtxprob = 0;
	lepton.deltaetacltrkcalo = 0;
	lepton.supcl_preshvsrawe = 0;
 
	lepton.closeTrackNLayers = 0;
	lepton.closeTrackNormChi2 = 999;
	
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
  
	lepton.posmatch = 0;
	lepton.nStations = 0;
	lepton.segmentComp = 0;
 
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
		float dzSig;
		float dxy;
		float dxyError;
		float dxySig;
		float puppiWeight;
		float puppiWeightNoLep;
		float trkChi2;
		float vertexChi2;
		float caloFraction;
		float hcalFraction;
		float hcalFractionCalib;
		int pdgId;
		int charge;
		int lostInnerHits;
		float pvAssocQuality;
		int nTrackerLayers;
		int pixelhits;
		int status;
		int time;
		bool trackHighPurity;
		bool isElectron;
		bool isMuon;
		bool isPhoton;
		bool isChargedHadron;
		bool isNeutralHadron;
		bool fromPV;
		bool hasTrackDetails;
};

void Initialize_PFlowCandidate(PFlowCandidate &pfcand){

		pfcand.pt = -100;
		pfcand.eta = -100;
		pfcand.phi = -100;
		pfcand.mass = -100;
		pfcand.phiAtVtx = -100;
		pfcand.dz = -0;
		pfcand.dzError = 0;
		pfcand.dzSig = 0;
		pfcand.dxy = 0;
		pfcand.dxyError = 0;
		pfcand.dxySig = 0;
		pfcand.puppiWeight = -100;
		pfcand.puppiWeightNoLep = -100;
		pfcand.trkChi2 = 999;
		pfcand.vertexChi2 = 999;
		pfcand.caloFraction = 0;
		pfcand.hcalFraction = 0;
		pfcand.hcalFractionCalib = -100;
		pfcand.pdgId =  0;
		pfcand.charge = 0;
		pfcand.lostInnerHits = -1;
		pfcand.pvAssocQuality = 0;
		pfcand.nTrackerLayers = 0;
		pfcand.pixelhits = 0;
		pfcand.status = 0;
		pfcand.time = -100;
		pfcand.trackHighPurity = false;
		pfcand.isElectron = false;
		pfcand.isMuon = false;
		pfcand.isPhoton = false;
		pfcand.isChargedHadron = false;
		pfcand.isNeutralHadron = false;
		pfcand.fromPV = false;
		pfcand.hasTrackDetails = false;
};
