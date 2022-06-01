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

// Including my object file //
# include "Objects.h"

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

int getGenPartonIndex(float lep_eta, float lep_phi, vector<GenParton>GenParts, int lep_pdgId, float dRcut=0.4)
{
	int i_gen_match = -1;
	
	for(unsigned igen=0; igen<GenParts.size(); igen++){
		if(GenParts[igen].status==1 && abs(GenParts[igen].pdgId)==lep_pdgId){
			if(delta2R(lep_eta,lep_phi,GenParts[igen].eta,GenParts[igen].phi)<dRcut){
				dRcut = delta2R(lep_eta,lep_phi,GenParts[igen].eta,GenParts[igen].phi);
				i_gen_match = (int)igen;
			}
		}
	}
	
	return i_gen_match;	
}

//int getGenPartonFlavor(float lep_eta, float lep_phi, vector<GenParton>GenParts, int lep_pdgId, float dRcut=0.4)
int getGenPartonFlavor(vector<GenParton>GenParts, int i_gen_match)
{	
	int genPartFlav = 0;
	
	if(i_gen_match>=0){	
		
		int mom_pdgId = (GenParts[i_gen_match].momstatus == 2) ? (abs(GenParts[i_gen_match].grmompdgId)) : (abs(GenParts[i_gen_match].mompdgId));
		
		if (GenParts[i_gen_match].isPromptFinalState){
			if (mom_pdgId==22){
				genPartFlav = 22;
			}
			else{
				genPartFlav = 1;
			}
		}
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
			genPartFlav = 0;
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
  
  double DR_PFCand_Maximum;
  double DR_PFCand_Minimum;
  double DR_SV_Maximum;
  double DR_AK4Jet_Maximum;
  
  // Root file & tree //
  
  TFile* theFile;
  
  TTree* T1;
    
  unsigned ievt;
  
  static const int njetmx = 20; 
  static const int npartmx = 50; 
  static const int nconsmax = 100; 
  static const int nsvmax = 10;
  
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
  int lepton_pdgId;
  float lepton_charge, lepton_tightcharge;
  //distance from PV & SV//
  float lepton_dxy, lepton_dxyError, lepton_dxySig, lepton_dz, lepton_dzError, lepton_dzSig, lepton_ip3d, lepton_sip3d, lepton_dxy_sv;
  //Track information //
  float lepton_chi2,  Lepton_valfrac;
  float lepton_ndof, lepton_trkKink, lepton_hit, lepton_pixhit, lepton_nTrackerLayers, lepton_lostHits;
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
  float lepton_eoverp, lepton_etain, lepton_dEtaInSeed, lepton_phiin, lepton_fbrem; 
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
  float lepton_isPFCand, lepton_isGlobal, lepton_isTracker;
  float lepton_isGoodGlobal, lepton_isTight, lepton_isHighPt, lepton_isHighPttrk, lepton_isMedium, lepton_isMedPr, lepton_isLoose;
  
  float lepton_jetPtRelv2, lepton_jetPtRelv2_log, lepton_jetRelIso, lepton_jetbtag;
  
  int label_Muon_Prompt, label_Muon_fromTau, label_Muon_fromHFHadron, label_Muon_fromLFHadron, label_Muon_fromPhoton, label_Muon_unknown, label_Muon_fromPhotonORunknown;
  int label_Electron_Prompt, label_Electron_fromTau, label_Electron_fromHFHadron, label_Electron_fromLFHadron, label_Electron_fromPhoton, label_Electron_unknown, label_Electron_fromPhotonORunknown;
  int label_unknown;
  bool label_Muon, label_Electron;
  bool label_fromTop, label_fromW, label_fromZ, label_fromH, label_fromNP, label_fromQCD, label_fromQCD_b, label_fromQCD_c, label_fromQCD_l, label_others, label_noGenMatch;
  
  // scalar end
  
  // GEN particle variables //
  
  int nGenPart;
  vector<float> GenPart_pt;
  vector<float> GenPart_eta;
  vector<float> GenPart_phi;
  vector<float> GenPart_mass;
  vector<float> GenPart_pdgId;
  vector<float> GenPart_mompdgId;
  vector<float> GenPart_grmompdgId;
  
  // PF candidates //
  
  int nPFCand;
  
  int nChargePFCand;
  vector<float> PFCand_pt_rel;
  vector<float> PFCand_pt_rel_log;
  vector<float> PFCand_eta_rel;
  vector<float> PFCand_phi_rel;
  vector<float> PFCand_phiAtVtx_rel;
  vector<float> PFCand_deltaR;
  vector<float> PFCand_puppiWeight;
  vector<float> PFCand_puppiWeightNoLep;
  vector<float> PFCand_caloFraction;
  vector<float> PFCand_hcalFraction;
  vector<float> PFCand_hcalFractionCalib;
  vector<float> PFCand_pdgId;
  vector<float> PFCand_energy_log;
  
  vector<float> PFCand_dz;
  vector<float> PFCand_dzError;
  vector<float> PFCand_dzSig;
  vector<float> PFCand_dxy;
  vector<float> PFCand_dxyError;
  vector<float> PFCand_dxySig;
  vector<float> PFCand_trkChi2;
  vector<float> PFCand_vertexChi2;
  vector<float> PFCand_charge;
  vector<float> PFCand_lostInnerHits;
  vector<float> PFCand_pvAssocQuality;
  vector<float> PFCand_nTrackerLayers;
  vector<float> PFCand_pixelhits;
  vector<float> PFCand_status;
  //vector<float> PFCand_time;
  vector<float> PFCand_trackHighPurity;
  vector<float> PFCand_isElectron;
  vector<float> PFCand_isMuon;
  vector<float> PFCand_isChargedHadron;
  vector<float> PFCand_fromPV;
  
  int nNeutralPFCand;
  /*
  vector<float> NeutralPFCand_pt_rel;
  vector<float> NeutralPFCand_eta_rel;
  vector<float> NeutralPFCand_deltaR;
  vector<float> NeutralPFCand_phi_rel;
  vector<float> NeutralPFCand_phiAtVtx_rel;
  vector<float> NeutralPFCand_puppiWeight;
  vector<float> NeutralPFCand_puppiWeightNoLep;
  vector<float> NeutralPFCand_caloFraction;
  vector<float> NeutralPFCand_hcalFraction;
  vector<float> NeutralPFCand_hcalFractionCalib;
  vector<int> NeutralPFCand_pdgId;
  vector<bool> NeutralPFCand_isPhoton;
  vector<bool> NeutralPFCand_isNeutralHadron;
  */
  // secondary vertices //
  int nSV;
  vector<float> SV_pt_rel;
  vector<float> SV_pt_rel_log;
  vector<float> SV_eta_rel;
  vector<float> SV_phi_rel;
  vector<float> SV_deltaR;
  vector<float> SV_mass;
  vector<float> SV_chi2;
  vector<float> SV_ntracks;
  vector<float> SV_ndof;
  vector<float> SV_dxy;
  vector<float> SV_dxyError;
  vector<float> SV_dxySig;
  vector<float> SV_d3d;
  vector<float> SV_d3dError;
  vector<float> SV_d3dSig;
  vector<float> SV_pAngle;
  vector<float> SV_cospAngle;
  
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
  
  DR_PFCand_Maximum = pset.getUntrackedParameter<double>("DR_PFCand_Maximum",0.5);
  DR_PFCand_Minimum = pset.getUntrackedParameter<double>("DR_PFCand_Minimum",5.e-3);
  DR_SV_Maximum = pset.getUntrackedParameter<double>("DR_SV_Maximum",0.5);
  DR_AK4Jet_Maximum = pset.getUntrackedParameter<double>("DR_AK4Jet_Maximum",0.4);
 
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
  
  T1->Branch("label_Muon_Prompt", &label_Muon_Prompt, "label_Muon_Prompt/I");
  T1->Branch("label_Muon_fromHFHadron", &label_Muon_fromHFHadron, "label_Muon_fromHFHadron/I");
  T1->Branch("label_Muon_fromLFHadron", &label_Muon_fromLFHadron, "label_Muon_fromLFHadron/I");
  T1->Branch("label_Muon_fromTau", &label_Muon_fromTau, "label_Muon_fromTau/I");
  T1->Branch("label_Muon_fromPhoton", &label_Muon_fromPhoton, "label_Muon_fromPhoton/I");
  T1->Branch("label_Muon_unknown", &label_Muon_unknown, "label_Muon_unknown/I");
  T1->Branch("label_Muon_fromPhotonORunknown", &label_Muon_fromPhotonORunknown, "label_Muon_fromPhotonORunknown/I");
  T1->Branch("label_Electron_Prompt", &label_Electron_Prompt, "label_Electron_Prompt/I");
  T1->Branch("label_Electron_fromTau", &label_Electron_fromTau, "label_Electron_fromTau/I");
  T1->Branch("label_Electron_fromHFHadron", &label_Electron_fromHFHadron, "label_Electron_fromHFHadron/I");
  T1->Branch("label_Electron_fromLFHadron", &label_Electron_fromLFHadron, "label_Electron_fromLFHadron/I");
  T1->Branch("label_Electron_fromPhoton", &label_Electron_fromPhoton, "label_Electron_fromPhoton/I");
  T1->Branch("label_Electron_unknown", &label_Electron_unknown, "label_Electron_unknown/I");
  T1->Branch("label_Electron_fromPhotonORunknown", &label_Electron_fromPhotonORunknown, "label_Electron_fromPhotonORunknown/I");
  T1->Branch("label_unknown", &label_unknown, "label_unknown/I");
  T1->Branch("label_Muon", &label_Muon, "label_Muon/O");
  T1->Branch("label_Electron", &label_Electron, "label_Electron/O");
  T1->Branch("label_fromTop", &label_fromTop, "label_fromTop/O");
  T1->Branch("label_fromW", &label_fromW, "label_fromW/O");
  T1->Branch("label_fromZ", &label_fromZ, "label_fromZ/O");
  T1->Branch("label_fromH", &label_fromH, "label_fromH/O");
  T1->Branch("label_fromNP", &label_fromNP, "label_fromNP/O");
  T1->Branch("label_fromQCD", &label_fromQCD, "label_fromQCD/O");
  T1->Branch("label_fromQCD_b", &label_fromQCD_b, "label_fromQCD_b/O");
  T1->Branch("label_fromQCD_c", &label_fromQCD_c, "label_fromQCD_c/O");
  T1->Branch("label_fromQCD_l", &label_fromQCD_l, "label_fromQCD_l/O");
  T1->Branch("label_others", &label_others, "label_others/O");
  T1->Branch("label_noGenMatch", &label_noGenMatch, "label_noGenMatch/O");
 
  // common kinematic variables //
  
  T1->Branch("lepton_pt", &lepton_pt,"lepton_pt/F");
  T1->Branch("lepton_p", &lepton_p,"lepton_p/F");
  T1->Branch("lepton_eta", &lepton_eta,"lepton_eta/F");
  T1->Branch("lepton_phi", &lepton_phi,"lepton_phi/F");
  T1->Branch("lepton_mass", &lepton_mass,"lepton_mass/F");
  
  T1->Branch("lepton_charge", &lepton_charge,"lepton_charge/F");
  T1->Branch("lepton_tightcharge", &lepton_tightcharge,"lepton_tightcharge/F");
  T1->Branch("lepton_pdgId", &lepton_pdgId,"lepton_pdgId/I");
  
  // common distance from PV & SV //
  
  T1->Branch("lepton_dxy", &lepton_dxy,"lepton_dxy/F");
  T1->Branch("lepton_dz", &lepton_dz,"lepton_dz/F");
  T1->Branch("lepton_dxyError", &lepton_dxyError,"lepton_dxyError/F");
  T1->Branch("lepton_dzError", &lepton_dzError,"lepton_dzError/F");
  T1->Branch("lepton_dxySig", &lepton_dxySig,"lepton_dxySig/F");
  T1->Branch("lepton_dzSig", &lepton_dzSig,"lepton_dzSig/F");
  T1->Branch("lepton_ip3d", &lepton_ip3d,"lepton_ip3d/F");
  T1->Branch("lepton_sip3d", &lepton_sip3d,"lepton_sip3d/F");
  T1->Branch("lepton_dxy_sv", &lepton_dxy_sv,"lepton_dxy_sv/F");
  
  // Generator-matching //
  
  T1->Branch("lepton_genPartFlav", &lepton_genPartFlav,"lepton_genPartFlav/I");
  
  // common Track information //
  
  T1->Branch("lepton_chi2", &lepton_chi2,"lepton_chi2/F");
  T1->Branch("lepton_ndof", &lepton_ndof,"lepton_ndof/F");
  T1->Branch("lepton_hit", &lepton_hit,"lepton_hit/F");
  T1->Branch("lepton_pixhit", &lepton_pixhit,"lepton_pixhit/F");
  T1->Branch("lepton_nTrackerLayers", &lepton_nTrackerLayers,"lepton_nTrackerLayers/F"); 
  T1->Branch("lepton_lostHits", &lepton_lostHits,"lepton_lostHits/F");
  
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
  
  if(save_only_muons && !save_only_electrons){
  
	// Muon-specfic variables //
  
	T1->Branch("lepton_trkKink", &lepton_trkKink,"lepton_trkKink/F");
	T1->Branch("lepton_nStations", &lepton_nStations,"lepton_nStations/I");
	T1->Branch("lepton_segmentComp", &lepton_segmentComp,"lepton_segmentComp/F");
	T1->Branch("lepton_posmatch", &lepton_posmatch,"lepton_posmatch/F");
  
	// Muon ID booleans //
  
	T1->Branch("lepton_isPFCand", &lepton_isPFCand,"lepton_isPFCand/F");
	T1->Branch("lepton_isGlobal", &lepton_isGlobal,"lepton_isGlobal/F");
	T1->Branch("lepton_isTracker", &lepton_isTracker,"lepton_isTracker/F");
	T1->Branch("lepton_isLoose", &lepton_isLoose,"lepton_isLoose/F");
	T1->Branch("lepton_isGoodGlobal", &lepton_isGoodGlobal,"lepton_isGoodGlobal/F");
	T1->Branch("lepton_isMedium", &lepton_isMedium,"lepton_isMedium/F");
	T1->Branch("lepton_isMedPr", &lepton_isMedPr,"lepton_isMedPr/F");
	T1->Branch("lepton_isTight", &lepton_isTight,"lepton_isTight/F");
	T1->Branch("lepton_isHighPt", &lepton_isHighPt,"lepton_isHighPt/F"); 
	T1->Branch("lepton_isHighPttrk", &lepton_isHighPttrk,"lepton_isHighPttrk/F");
  
  }
  
  if(!save_only_muons && save_only_electrons){
  
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
	T1->Branch("lepton_dEtaInSeed", &lepton_dEtaInSeed,"lepton_dEtaInSeed/F");
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
  
  }
  
  // Nearest jet features //
  
  T1->Branch("lepton_jetPtRelv2", &lepton_jetPtRelv2, "lepton_jetPtRelv2/F");
  T1->Branch("lepton_jetPtRelv2_log", &lepton_jetPtRelv2_log, "lepton_jetPtRelv2_log/F");
  T1->Branch("lepton_jetRelIso", &lepton_jetRelIso, "lepton_jetRelIso/F");
  T1->Branch("lepton_jetbtag", &lepton_jetbtag, "lepton_jetbtag/F");
  
  // PF candidates //
  
  T1->Branch("nPFCand",&nPFCand,"nPFCand/I");
  
  T1->Branch("PFCand_pt_rel",&PFCand_pt_rel);
  T1->Branch("PFCand_pt_rel_log",&PFCand_pt_rel_log);
  T1->Branch("PFCand_eta_rel",&PFCand_eta_rel);
  T1->Branch("PFCand_phi_rel",&PFCand_phi_rel);
  T1->Branch("PFCand_phiAtVtx_rel",&PFCand_phiAtVtx_rel);
  T1->Branch("PFCand_deltaR",&PFCand_deltaR);
  //T1->Branch("PFCand_mass",&PFCand_mass);
  T1->Branch("PFCand_pdgId",&PFCand_pdgId);
  T1->Branch("PFCand_caloFraction",&PFCand_caloFraction);
  T1->Branch("PFCand_hcalFraction",&PFCand_hcalFraction);
  T1->Branch("PFCand_hcalFractionCalib",&PFCand_hcalFractionCalib);
  T1->Branch("PFCand_puppiWeight",&PFCand_puppiWeight);
  T1->Branch("PFCand_puppiWeightNoLep",&PFCand_puppiWeightNoLep);
  T1->Branch("PFCand_energy_log",&PFCand_energy_log);
  
  T1->Branch("PFCand_dz",&PFCand_dz);
  T1->Branch("PFCand_dzError",&PFCand_dzError);
  T1->Branch("PFCand_dzSig",&PFCand_dzSig);
  T1->Branch("PFCand_dxy",&PFCand_dxy);
  T1->Branch("PFCand_dxyError",&PFCand_dxyError);
  T1->Branch("PFCand_dxySig",&PFCand_dxySig);
  //T1->Branch("PFCand_vertexChi2",&PFCand_vertexChi2);
  //T1->Branch("PFCand_time",&PFCand_time);
  T1->Branch("PFCand_charge",&PFCand_charge);
  T1->Branch("PFCand_lostInnerHits",&PFCand_lostInnerHits);
  T1->Branch("PFCand_pvAssocQuality",&PFCand_pvAssocQuality);
  T1->Branch("PFCand_status",&PFCand_status);
  T1->Branch("PFCand_pixelhits",&PFCand_pixelhits);
  T1->Branch("PFCand_nTrackerLayers",&PFCand_nTrackerLayers);
  T1->Branch("PFCand_trkChi2",&PFCand_trkChi2);
  T1->Branch("PFCand_trackHighPurity",&PFCand_trackHighPurity);
  T1->Branch("PFCand_isElectron",&PFCand_isElectron);
  T1->Branch("PFCand_isMuon",&PFCand_isMuon);
  T1->Branch("PFCand_isChargedHadron",&PFCand_isChargedHadron);
  T1->Branch("PFCand_fromPV",&PFCand_fromPV);
  /*
  T1->Branch("NeutralPFCand_pt_rel",&NeutralPFCand_pt_rel);
  T1->Branch("NeutralPFCand_eta_rel",&NeutralPFCand_eta_rel);
  T1->Branch("NeutralPFCand_phi_rel",&NeutralPFCand_phi_rel);
  T1->Branch("NeutralPFCand_phiAtVtx_rel",&NeutralPFCand_phiAtVtx_rel);
  T1->Branch("NeutralPFCand_deltaR",&NeutralPFCand_deltaR);
  //T1->Branch("NeutralPFCand_mass",&NeutralPFCand_mass);
  T1->Branch("NeutralPFCand_pdgId",&NeutralPFCand_pdgId);
  T1->Branch("NeutralPFCand_caloFraction",&NeutralPFCand_caloFraction);
  T1->Branch("NeutralPFCand_hcalFraction",&NeutralPFCand_hcalFraction);
  T1->Branch("NeutralPFCand_hcalFractionCalib",&NeutralPFCand_hcalFractionCalib);
  T1->Branch("NeutralPFCand_puppiWeight",&NeutralPFCand_puppiWeight);
  T1->Branch("NeutralPFCand_puppiWeightNoLep",&NeutralPFCand_puppiWeightNoLep);
  T1->Branch("NeutralPFCand_isPhoton",&NeutralPFCand_isPhoton);
  T1->Branch("NeutralPFCand_isNeutralHadron",&NeutralPFCand_isNeutralHadron);
  */
  // Secondary vertices //
  
  T1->Branch("nSV",&nSV,"nSV/I");
  T1->Branch("SV_pt_rel",&SV_pt_rel);
  T1->Branch("SV_pt_rel_log",&SV_pt_rel_log);
  T1->Branch("SV_mass",&SV_mass);
  T1->Branch("SV_eta_rel",&SV_eta_rel);
  T1->Branch("SV_phi_rel",&SV_phi_rel);
  T1->Branch("SV_deltaR",&SV_deltaR);
  T1->Branch("SV_chi2",&SV_chi2);
  T1->Branch("SV_ntracks",&SV_ntracks);
  T1->Branch("SV_ndof",&SV_ndof);
  T1->Branch("SV_dxy",&SV_dxy);
  T1->Branch("SV_dxyError",&SV_dxyError);
  T1->Branch("SV_dxySig",&SV_dxySig);
  T1->Branch("SV_d3d",&SV_d3d);
  T1->Branch("SV_d3dError",&SV_d3dError);
  T1->Branch("SV_d3dSig",&SV_d3dSig);
  T1->Branch("SV_pAngle",&SV_pAngle);
  T1->Branch("SV_cospAngle",&SV_cospAngle);
  
  // MC Info //
  
  if(isMC){
	 
	T1->Branch("Generator_weight", &Generator_weight, "Generator_weight/D") ;
	T1->Branch("Pileup_nPU",&Pileup_nPU,"Pileup_nPU/I");
	T1->Branch("Pileup_nTrueInt",&Pileup_nTrueInt,"Pileup_nTrueInt/I");
	
	T1->Branch("nGenPart",&nGenPart,"nGenPart/I");
	T1->Branch("GenPart_pt",&GenPart_pt);
	T1->Branch("GenPart_eta",&GenPart_eta);
	T1->Branch("GenPart_phi",&GenPart_phi);
	T1->Branch("GenPart_mass",&GenPart_mass);
	T1->Branch("GenPart_pdgId",&GenPart_pdgId);
	T1->Branch("GenPart_mompdgId",&GenPart_mompdgId);
	T1->Branch("GenPart_grmompdgId",&GenPart_grmompdgId);
	T1->Branch("PFCand_hcalFraction",&PFCand_hcalFraction);
	T1->Branch("PFCand_hcalFractionCalib",&PFCand_hcalFractionCalib);
	T1->Branch("PFCand_puppiWeight",&PFCand_puppiWeight);
	T1->Branch("PFCand_puppiWeightNoLep",&PFCand_puppiWeightNoLep);
  
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
