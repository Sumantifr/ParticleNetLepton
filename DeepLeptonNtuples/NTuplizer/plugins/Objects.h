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
