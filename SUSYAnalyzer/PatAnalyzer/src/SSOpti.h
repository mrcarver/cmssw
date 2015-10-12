#ifndef _L1ITMu_SSOpti_h_
#define _L1ITMu_SSOpti_h_

#include <memory>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <TStyle.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>

#include "L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#define MAXOBJECT 30

typedef edm::ParameterSet PSet;


class SSOpti : public edm::EDAnalyzer {

//class L1TMuonTextDumper : public edm::EDProducer {
public:
  SSOpti(const PSet&);
  ~SSOpti() {}
	private:
 void analyze(const edm::Event&, const edm::EventSetup&); 
 //void produce(edm::Event&, const edm::EventSetup&); 
  void beginJob();
  void endJob();
  float CalcdR(float eta1, float phi1, float eta2, float phi2);
  int FindFirstMother(reco::GenParticle gp);
  
  FILE *write, *dphi, *tptest;
  
  TFile* fout;
  edm::Service<TFileService> fs;
  TNtuple* theNtuple;
  
  unsigned int Sample;
  //char Sample;
  
  ///////////////////////
  //// Ntuple Objects ///
  ///////////////////////

  int ngenJets, ngenLep, nHLTLep, nL1Lep, DblMuTrig, DblMuTrigNS, DblElTrig, DblElTrigNS4, DblElTrigNS6, EMuTrig, EMuTrigNS, L1M6HT150, L1M8HT125, L1HT175, L1DEG6HT150, L1DEG4HT150;
  int genLepId, HLTLepId, L1LepId, genLepOrigin, GENHT, HLTHT, L1HT;
  
  float genJetPt[MAXOBJECT], genJetEta[MAXOBJECT], genJetPhi[MAXOBJECT], genLepPt[MAXOBJECT], genLepEta[MAXOBJECT], genLepPhi[MAXOBJECT], genLepIso[MAXOBJECT], HLTLepPt[MAXOBJECT], HLTLepEta[MAXOBJECT], HLTLepPhi[MAXOBJECT], L1LepPt[MAXOBJECT], L1LepEta[MAXOBJECT], L1LepPhi[MAXOBJECT];
  
  edm::ParameterSet LUTparam;
 // CSCSectorReceiverLUT* srLUTs_[5][2];
  const float ptscale[33] = { 
  	-1.,   0.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,
    4.5,   5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,  
    16.0,  18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0, 
    50.0,  60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0, 1.E6 };
    

   
   
   ////////////////////////////
   //// Variables For Tree ////
   ////////////////////////////
   TTree* outputTree;
   
   int _nL1CJets;
   int _nL1Muons;
   int _nL1EGs;
   int _nGenJets;
   int _n_HLT_PFJets;
   int _nGenMuons;
   int _nGenElectrons;
   
   double _L1CJet_eta[50], _L1CJet_phi[50], _L1CJet_pt[50];
   double _L1Muon_eta[50], _L1Muon_phi[50], _L1Muon_pt[50];
   double _L1EG_eta[50], _L1EG_phi[50], _L1EG_pt[50];
   double _GenJet_eta[50], _GenJet_phi[50], _GenJet_pt[50];
   double _PFJet_eta[50], _PFJet_phi[50], _PFJet_pt[50];
   double _GenMuon_eta[50], _GenMuon_phi[50], _GenMuon_pt[50];
   double _GenElectron_eta[50], _GenElectron_phi[50], _GenElectron_pt[50]; 
   
   
   double _L1_HT, _HLT_HT, _Gen_HT, _Gen_HT2;
   double GPLepton[3][4], RSLeptonTest[4][3][301];//, RSLepton[3][301];
   double GPLepton2[2][3][4], GPLepton2HB[2][3][4], RSLeptonTest2[2][4][3][301], ASLeptonTest2[2][4][3][301], RSLeptonTest2HB[2][4][3][301], ASLeptonTest2HB[2][4][3][301];
   	///[ptrange][prompt/nonpropmt]  then[ptrange][isocut]
   
   TH1F* GLEta, *RecoLeptonPt[4];
   
   ///// Gen Plots /////
   
   TH1F* dRQL[2][3], *dRJL[2][3], *dRQM[3], *dRJM[3], *dRQE[3], *dRJE[3];
   TH1F* GenLIso[2][3], *GenMIso[3], *GenEIso[3];
   
   
   ///// RECO Plots /////
   
   TH1F* RecodRLRJ[4][3];
   TH1F* RecodRLGJ[3];
   TH1F* RecodRLGQ[4][3];
   
   TH1F* RecoLIsoPrompt[2][3], *RecoLIsoNotPrompt[3];
   TH1F* RecoLAIsoPrompt[2][3], *RecoLAIsoNotPrompt[3];//finished
   
   TH1F* RecoLIsoPromptHB[2][3], *RecoLIsoNotPromptHB[3];//finished
   TH1F* RecoLAIsoPromptHB[2][3], *RecoLAIsoNotPromptHB[3];//finished
   
   TH1F* RecoLIsoPrompt2[2][2][3], *RecoLIsoNotPrompt2[2][3];
   TH1F* JetCSV[2][3][2], *LeptonJetPtRatio[2][3][2];//[e/mu][ptrange][prompt/nonprompt]
   
   
   
   ///// Efficiency Plots /////
   
   TH1F* GenPromptLepton, *EffvsRelIso[3][4], *GEffvsRelIso[3][4];///[ptrange][prompt/nonpropmt]
   TH1F* EffvsRelIso2[2][3][4], *GEffvsRelIso2[2][3][4];///[e/mu][ptrange][prompt/nonprompt]
   
   TH1F* EffvsAIso2[2][3][4], *GEffvsAIso2[2][3][4];///[e/mu][ptrange][prompt/nonprompt]//finished
   
   TH1F* EffvsRelIso2HB[2][3][4], *GEffvsRelIso2HB[2][3][4];///[e/mu][ptrange][prompt/nonprompt]
   TH1F* EffvsAIso2HB[2][3][4], *GEffvsAIso2HB[2][3][4];///[e/mu][ptrange][prompt/nonprompt]
   

  bool _dogen;
  edm::InputTag _geninput;
  std::vector<edm::InputTag> _tpinputs, _convTrkInputs;
  edm::Service<TFileService> histofile;
  
  
};




#endif
