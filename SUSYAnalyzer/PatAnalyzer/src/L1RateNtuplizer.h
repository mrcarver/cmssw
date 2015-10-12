#ifndef L1RateNtuplizer_H
#define L1RateNtuplizer_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "SUSYAnalyzer/PatAnalyzer/interface/GenParticleManager.h"
#include "SUSYAnalyzer/PatAnalyzer/interface/Statistics.h"
#include "SUSYAnalyzer/PatAnalyzer/interface/Tools.h"
#include "SUSYAnalyzer/PatAnalyzer/interface/OnTheFlyCorrections.hh"


#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"




//Root Classes

#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TClonesArray.h"

//Standard C++ classes
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <ostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <memory>
#include <iomanip>

using namespace std;

const int nLeptonsMax = 8;

class L1RateNtuplizer : public edm::EDAnalyzer {
public:
    
    explicit L1RateNtuplizer(const edm::ParameterSet & iConfig);
	
    ~L1RateNtuplizer(){};
    
private:
    
    //virtual void analyze(edm::Event & iEvent, const edm::EventSetup & iSetup);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob(void);
    
   // void fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu);
    //void fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* el);
    
    std::string Sample;
    edm::InputTag IT_muon;
    edm::InputTag IT_electron;
    edm::InputTag IT_tau;
    edm::InputTag IT_tauDiscriminator;
    edm::InputTag IT_jet;
    edm::InputTag IT_pfmet;
    edm::InputTag IT_beamspot;
    edm::InputTag IT_hltresults;
    edm::InputTag IT_METFilters;
	//edm::InputTag electronsCollection_;
	edm::InputTag MVAidCollection_;
	
	EGammaMvaEleEstimatorCSA14* myMVATrig;
    
    edm::Service<TFileService> fs;
    FILE *outfile;
    
    TH1F *Nvtx;
    
    //desired output variables
    TTree* outputTree;
    
    string _corrLevel;
    
    
    double _relIsoCutE;
    double _relIsoCutMu;
    double _relIsoCutEloose;
    double _relIsoCutMuloose;
    
    bool _chargeConsistency;
    
    double _minPt0;
    double _minPt1;
    double _tightD0Mu;
    double _tightD0E;
    double _looseD0Mu;
    double _looseD0E;
	double _weight;
    
    double _jetPtCut;
    double _jetEtaCut;
    
    double _tauPt;
    double _tauEta;
    
    bool _regression;
    
    
    //genlevel particles
    GenParticleManager GPM;
    OnTheFlyCorrections* fMetCorrector;
    
    int _n_bJets;
    int _n_Jets;
    
    double _jetEta[20];
    double _jetPhi[20];
    double _jetPt[20];
    bool _bTagged[20];
    double _csv[20];
    
    int _n_bJetsAll;
    int _n_JetsAll;
    
    double _jetEtaAll[100];
    double _jetPhiAll[100];
    double _jetPtAll[100];
    bool _bTaggedAll[100];
    double _csvAll[100];
    int _leptonIndex;
    int _closeIndex[nLeptonsMax];
    
    
    TClonesArray* _leptonP4;
    TClonesArray* _jetP4;
    TClonesArray* _jetAllP4;
	
	double _L1Mu_Pt[40];
	double _L1Mu_Eta[40];
	double _L1Mu_Phi[40];
	int _nL1Mus;
	
	
	TClonesArray* _Mu17Mu8Objs;
	TClonesArray* _Mu17TkMu8Objs;
	TClonesArray* _Ele23Ele12Objs;
	TClonesArray* _Mu23Ele12Objs;
	TClonesArray* _Ele23Mu8Objs;
	
	TClonesArray* _DblMuHTObjs;
	TClonesArray* _DblEleHTObjs;
	TClonesArray* _MuEleHTObjs;
	
	TClonesArray* _MuBtagObjs;
	TClonesArray* _EleBtagObjs;
	
	TClonesArray* _MuControlObjs[4];
	TClonesArray* _IsoMuControlObjs[4];
	TClonesArray* _EleControlObjs[5];
	TClonesArray* _IsoEleControlObjs[4];
	
	int _nMu17Mu8Objs;
	int _nMu17TkMu8Objs;
	int _nEle23Ele12Objs;
	int _nMu23Ele12Objs;
	int _nEle23Mu8Objs;
	
	int _nDblMuHTObjs;
	int _nDblEleHTObjs;
	int _nMuEleHTObjs;
	int _nMuBtagObjs;
	int _nEleBtagObjs;
	
	int _nMuControlObjs[4];
	int _nIsoMuControlObjs[4];
	int _nEleControlObjs[5];
	int _nIsoEleControlObjs[4];
    
    int _nLeptons;
    
    int _eventType; //ee,mm,em
    bool _sb;
    bool _doubleF;
    int _index1 = -1;
    int _index2 = -1;

    int _PFHT800;
	double _hltHT;
	double _L1HT;
	
	int _Mu17Mu8;
	int _Mu17TkMu8;
	int _Ele23Ele12;
	int _Mu23Ele12;
	int _Ele23Mu8;
	
	int _DblMuHT;
	int _DblEleHT;
	int _MuEleHT;
	
	int _MuBtag;
	int _EleBtag;
	
	int _MuControl[4];
	int _IsoMuControl[4];
	int _EleControl[5];
	int _IsoEleControl[4];
	
    int _indeces[nLeptonsMax];
    int _flavors[nLeptonsMax];
	int _pdgids[nLeptonsMax];
    double _charges[nLeptonsMax];
    double _isolation[nLeptonsMax];
	double _miniIsolation[nLeptonsMax];
	double _miniIsolationEA[nLeptonsMax];
	double _MVAVal[nLeptonsMax];
	double _ptrel[nLeptonsMax];
    double _isolationComponents[nLeptonsMax][4];
	int _inheritance[nLeptonsMax][30];
    double _isolationMC[nLeptonsMax][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    
    int _origin[nLeptonsMax];
    int _originReduced[nLeptonsMax];
    double _PVchi2;
    double _PVerr[3];
    double _ipPV[nLeptonsMax];
    double _ipPVerr[nLeptonsMax];
    double _ipPVmc[nLeptonsMax];
    
    double _ipZPV[nLeptonsMax];
    double _ipZPVerr[nLeptonsMax];
    
    double _3dIP[nLeptonsMax];
    double _3dIPerr[nLeptonsMax];
    double _3dIPsig[nLeptonsMax];
	
	int _missingHits[nLeptonsMax];
	double _dphi[nLeptonsMax];
	double _deta[nLeptonsMax]; 
	double _sigIeta[nLeptonsMax];
	double _HoE[nLeptonsMax];
	double _pMe[nLeptonsMax];
	
    double _mt[nLeptonsMax];
    
    double _closeJetPt[nLeptonsMax];
    double _closeJetPtAll[nLeptonsMax];
    
    double _closeJetAngAll[nLeptonsMax];
    double _ptRelAll[nLeptonsMax];
    double _closeJetPtAllMC[nLeptonsMax];
    double _closeJetPtAllstatus[nLeptonsMax];
    int _partonIdMatched[nLeptonsMax];
    bool _sameParton[nLeptonsMax];
    
    bool _isloose[nLeptonsMax];
    bool _istight[nLeptonsMax];
	bool _istightNoIso[nLeptonsMax];
	bool _istightNoIsoSIP[nLeptonsMax];
	bool _islooseMVA[nLeptonsMax];
	bool _istightMVANoIsoSIP[nLeptonsMax];
	bool _istightMVANoIsoSIP_LMVA[nLeptonsMax];
	bool _istightMVA[nLeptonsMax];
	bool _chargeConsistent[nLeptonsMax];
	
	int _genCharge[nLeptonsMax];
	double _genHT;
    
    int _n_PV;
    
    int _n_electrons;
    int _n_muons;
    int _n_taus;
    
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;
    

    double _mompt[nLeptonsMax];
    double _momphi[nLeptonsMax];
    double _mometa[nLeptonsMax];
    int _mompdg[nLeptonsMax];
    
    
	double _rho;
    double _met;
    double _met_phi;
    double HT;
    
    long _nEventsTotal;
    long _nEventsTotalCounted;
    long _nEventsFiltered;
    
    TH1D* _hCounter;
	TH1D* MvaHist;
    
    double _regVars[15];
    double hJet_ptRaw;
    double hJet_genPt;
    double hJet_pt;
    double hJet_phi;
    double hJet_eta;
    double hJet_e;
    
    double hJet_ptLeadTrack;
    
    double hJet_vtx3dL;
    double hJet_vtx3deL;
    double hJet_vtxMass;
    double hJet_vtxPt;
    
    double hJet_cef;
    double hJet_nconstituents;
    double hJet_JECUnc;
    
    double hJet_SoftLeptptRel;
    double hJet_SoftLeptPt;
    double hJet_SoftLeptdR;
    
    double hJet_SoftLeptIdlooseMu;
    double hJet_SoftLeptId95;
	
	TH1F *Triggers;
	
	int _DiLepHTTrigs[3], _nDiLepHTObjs[3];
	TClonesArray* _DiLepHTObjs[3];
	TString DiLepHTTrigNames[3] = {"HLT_DoubleMu8_Mass8_HTT300","HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300","HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300"};
	
	int _DiLepTrigs[5], _nDiLepObjs[5];
	TClonesArray* _DiLepObjs[5];
	TString DiLepTrigNames[5] = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
								 "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
								 "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
								 "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
								 "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"};
					//[flav][iso?][thresh]
	int _ControlTrigs[2][2][5], _nControlObjs[2][2][5];
	TClonesArray* _ControlObjs[2][2][5];
	TString ControlTrigNames[2][2][5] = {{{"HLT_Ele8_CaloIdM_TrackIdM_PFJet30","HLT_Ele12_CaloIdM_TrackIdM_PFJet30","HLT_Ele18_CaloIdM_TrackIdM_PFJet30","HLT_Ele23_CaloIdM_TrackIdM_PFJet30","HLT_Ele33_CaloIdM_TrackIdM_PFJet30"},
										  {"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele34_CaloIdL_TrackIdL_IsoVL_PFJet30"}},
										 {{"HLT_Mu8_v","HLT_Mu17_v","HLT_Mu24_v","HLT_Mu34_v","HLT_Mu35_v"},
										  {"HLT_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_v","HLT_Mu24_TrkIsoVVL_v","HLT_Mu34_TrkIsoVVL_v","HLT_Mu35_TrkIsoVVL_v"}}};
	
	int _BControlTrigs[2], _nBControlObjs[2];
	TClonesArray* _BControlObjs[2];
	TString BControlTrigNames[2] = {"HLT_Mu10_CentralPFJet30_BTagCSV0p5","HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p5"};
	
};

#endif
