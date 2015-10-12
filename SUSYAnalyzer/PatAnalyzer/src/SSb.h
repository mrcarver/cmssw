#ifndef SSb_H
#define SSb_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

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


#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"


#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"


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

class SSb : public edm::EDAnalyzer {
public:
    
    explicit SSb(const edm::ParameterSet & iConfig);
    ~SSb(){};
    
private:
    
    //virtual void analyze(edm::Event & iEvent, const edm::EventSetup & iSetup);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob(void);
    
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
    
    edm::Service<TFileService> fs;

    double _relIsoCutE;
    double _relIsoCutMu;
    bool _chargeConsistency;
    
    double _minPt0;
    double _minPt1;
    double _tightD0e;
    double _tightD0mu;
    double _looseD0;
    
    double _jetPtCut;
    double _jetEtaCut;
    
    double _tauPt;
    double _tauEta;

    
   //genlevel particles
    GenParticleManager GPM;
    JetCorrectionUncertainty* _totalJEC;
    
    //OutputFile
    TString outputFile;
    
    //desired output variables
    TTree* outputTree;
    TH1D* _hCounter;
    long _nEventsTotal;
    long _nEventsTotalCounted;
    long _nEventsFiltered;
    
    TClonesArray* _leptonP4;
    TClonesArray* _jetP4;
    TClonesArray* _jetAllP4;
    TClonesArray* _closeJetP4;
    
    double _parton20HT;

    int _indeces[10];
    int _flavors[10];
    double _charges[10];
    double _isolation[10];
    int _noniso;
    int _iso;
    int _origin0[10];
    int _originReduced0[10];
    bool _fromTop[10];
    bool _fromBjet[10];
    
    bool _chargeMisId0[10];
    double _ipPV[10];
    double _ipPVerr[10];

    double _ipZPV[10];
    double _ipZPVerr[10];
    
    double _3dIP[10];
    double _3dIPerr[10];
    double _3dIPsig[10];
    
    double _PVerror[3];
    
    double _closeJetPtAll[10];
    double _closeJetAngAll[10];
    double _ptRelAll[10];
    
    int _eventType; //ee,mm,em
    bool _sb;
    bool _doubleF;
    double _mll;
    double _mtHard;
    double _ht;
    double _chargePair;
    int _index1;
    int _index2;
    
    
    int _n_b_MC_pt30[3];
    int _n_b_MC[3];
    
    int _n_bJets;
    int _n_Jets;
    
    double _jetEta[50];
    double _jetPhi[50];
    double _jetPt[50];
    bool _bTagged[50];
    
    int _n_bJetsAll;
    int _n_JetsAll;
    int _n_JetsnoTau;
    int _nLeptons;
    
    int _n_bJetsAll30;
    int _n_JetsAll30;
    
    int _mcFlavor[50];
    
    double _jetEtaAll[50];
    double _jetPhiAll[50];
    double _jetPtAll[50];
    bool _bTaggedAll[50];

    bool _bTaggedAllSSVHE[50];
    bool _bTaggedAllSSVHP[50];
    bool _bTaggedAllTCHE[50];
    bool _bTaggedAllTCHP[50];

    int _n_PV;
    
    
    int _n_electrons;
    int _n_muons;
    int _n_taus;
    int _n_tausPreSel;
    
    long _eventNb;
    long _runNb;
    long _lumiBlock;
    
    bool _trigDiEle;
    bool _trigDiMu;
    bool _trigEMu;
    bool _trigMuE;
    
    bool passedEMu_HT;
    bool passedDiMu_HT;
    bool passedDiEle_HT;
    bool passedEMu_MET;
    bool passedDiMu_MET;
    bool passedDiEle_MET;
    
    int _sameSign[2][2];
    
    
    double pfMET;
    double pfMET_Phi;
    double HT;
    
    double pfMETtype1;
    double pfMETtype1_Phi;
    
    
    float _mgluino;
    float _mchi0;
    float _xparameter;
    
    long nEventsTotal;
    long nEventsFiltered;
    
    TH2F* hElectronsMuonsPre;
    TH2F* hElectronsMuons;
    TH2F* hTauEMu;
    
    int _origin[5][10];
    int _originReduced[5][10];
    int _pdgId[5][10];
    int _originMatched[5][10];
    bool _chargeMisId[5][10];
    double _chargeMC[5][10];

    //HLT paths
    std::vector< std::string> vMuMu_HLTs;
    std::vector< std::string> vElEl_HLTs;
    std::vector< std::string> vElMu_HLTs;
    bool requireHLT;
    
    
};

#endif
