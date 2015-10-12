#ifndef FakeMuons_H
#define FakeMuons_H

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

class FakeMuons : public edm::EDAnalyzer {
public:
    
    explicit FakeMuons(const edm::ParameterSet & iConfig);
    ~FakeMuons(){};
    
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
    
    FILE *outfile;

    string _corrLevel;
    
//    double _relIsoCut;
    double _relIsoCutE;
    double _relIsoCutEloose;
    double _relIsoCutMu;
    double _relIsoCutMuloose;

    bool _chargeConsistency;
    
    bool _NonIsoMuonTrig;
    bool _NonIsoElectronTrig;

    
    double _minPt0;
    double _minPt1;
    double _tightD0Mu;
    double _tightD0E;
    double _looseD0Mu;
    double _looseD0E;
    
    double _jetPtCut;
    double _jetEtaCut;
    
    double _tauPt;
    double _tauEta;
    
    //genlevel particles
    GenParticleManager GPM;
    
    OnTheFlyCorrections* fMetCorrector;
    
    reco::tau::RecoTauVertexAssociator* tauVtx;

    //OutputFile
    TString outputFile;

    //desired output variables
    TTree* outputTree;
    
    TH1D* _triggers;
    
    int _eventType; //ee,mm,em
    
    int _n_bJets;
    int _n_Jets;
    
    
    double _jetEta[100];
    double _jetPhi[100];
    double _jetPt[100];
    bool _bTagged[100];

    double _csv[100];
    double _csvAll[100];

    
    int _n_bJetsAll;
    int _n_JetsAll;
    int _nLeptons;
    
    double _jetEtaAll[100];
    double _jetPhiAll[100];
    double _jetPtAll[100];
    bool _bTaggedAll[100];

    TClonesArray* _leptonP4;
    TClonesArray* _jetP4;
    TClonesArray* _jetAllP4;

    bool _csExercise;
    int _csJetIndex;
    
    int _nonisoOverlap;
    
    bool _bOverlap[3];
    int _indeces[3];
    int _flavors[3];
    double _charges[3];
    double _isolation[3];
    double _isolationMC[3][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    int _noniso;
    int _origin0[3];
    int _originReduced0[3];
    bool _chargeMisId0[3];
    double _ipPV[3];
    double _mt0[3];
    double _angleWl[3];
    double _mll[3];
    
    double _closeJetPt[3];
    double _closeJetPtAll[3];
    double _closeJetAng[3];
    double _closeJetAngAll[3];
    double _ptRel[3];
    double _ptRelAll[3];
    
    double _closeJetPtAllMC[3];
    int _partonIdMatched[10][2];
    bool _sameParton[3];
    
    double _closeJetPtAllNotItself[3];
    double _closeJetAngAllNotItself[3];
    double _ptRelAllNotItself[3];
    
    int _n_bJetsAll30;
    int _n_JetsAll30;
   
    bool _isloose[3];
    bool _istight[3];
    
    int _n_PV;
    
    
    int _n_electrons;
    int _n_muons;
    int _n_taus;
    int _n_tausPreSel;
    
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;

    bool _trigDiEle;
    bool _trigEleHad;
    
    bool _trigDiEle_CS;
    bool _trigDiEle_CS_2;
    bool _trigEleHad_CS;
    
    bool _trigDiMu;
    bool _trigEMu;
    bool _trigMuE;
    
    bool _trigSingMu;
    
    bool _trgLepNjets[3][5];
    bool _trgSingLep[2][4];
    
    int _sameSign[3][2];
    
    double _charge[3][3];
    double _pt[3][3];
    
    double _chargeMC[5][3];
    int _origin[5][3];
    int _originReduced[5][3];
    int _pdgId[5][3];
    int _originMatched[5][3];
    bool _chargeMisId[5][3];
    
    double _d0_tau[3];
    double _dz_tau[3];
    double _vtx_tau_comty[3];
    
    
    double _ip[3][3];
    double _ipSig[3][3];
    double _eta[3][3];
    double _phi[3][3];
    bool _remove[3][3];

    
    double _iso_trk[2][3];
    double _iso_hcal[2][3];
    double _iso_ecal[2][3];
    double _iso_det[2][3];
    double _iso_pf[2][3];
    double _iso_pf_parts[5][3];
    
    double _mt[3];
    double _Et[3][3];
    double _px[3][3];
    double _py[3][3];
    double _pz[3][3];
    double _E[3][3];
    double _mInvPair[6][6];
        

    int _nHard;
    double _hardIpt[8];
    double _hardIphi[8];
    double _hardIeta[8];
    int _hardIpdg[8];
    
    double _mompt;
    double _momphi;
    double _mometa;
    int _mompdg;

    double _leppt[8][4];
    double _lepphi[8][4];
    double _lepeta[8][4];
    double _lepiso[8][4];
    int _leppdg[8][4];
    int _lepNum[8];
    
    double _met;
    double _met_phi;
    double HT;
    
    long _nEventsTotal;
    long _nEventsTotalCounted;
    long _nEventsFiltered;
    
    TH1D* _hCounter;
    
    TH1F* hElectronOrigin;
    TH1F* hMuonOrigin;

    TH2F* hRecoElectronOrigin;
    TH2F* hRecoMuonOrigin;
    
    TH2F* hElectronsMuonsPre;
    TH2F* hElectronsMuons;
    TH2F* hTauEMu;

    
    //HLT paths
    std::vector< std::string> vMuMu_HLTs;
    std::vector< std::string> vElEl_HLTs;
    std::vector< std::string> vElMu_HLTs;
    bool requireHLT;
    
    //MET cuts
    //  double value_pfmet;
    
    //jet cuts
    double pt_jet;
    double eta_jet;
    bool   jetLeptonCleaning;
    
    //muon cuts
    double pt_mu;
    double d0_mu;
    double reliso_mu;
    bool usePFIso_Mu;
    
    //electron cuts
    double pt_el;
    double d0_el;
    double reliso_el;
    double usePFIso_El;
    bool chargeConsistency;
    
    edm::Service<TFileService> fs;
    
    float Iso_Cut[7];
    float Pt_Cut[3];
    float Pt_Cut_J[7];
    float ptb[14];
    

    
    
    //histograms
    TH1F *hMET, *hHT, *hNJets, *hMuonPt;
    TH1F *MET_DoubleMu_Den, *MET_DoubleMu_Num, *MET_DoubleEl_Den, *MET_DoubleEl_Num, *MET_DoubleEM_Den, *MET_DoubleEM_Num;
    TH1F *HT_DoubleMu_Den, *HT_DoubleMu_Num, *HT_DoubleEl_Den, *HT_DoubleEl_Num, *HT_DoubleEM_Den, *HT_DoubleEM_Num;
    
    TH1F *Nevt_jet;
    TH1F *Nvtx;
    
    TH1F *Iso_Lep_Btg[7][9][5];
    TH1F *Iso_Lep_Btg_Det[7][9][5];
    TH1F *dxy_Lep_Btg[7][9][5];
	TH1F *pt_Lep_Btg[7][9];
	TH1F *pt_Lep_Btg_Iso[7][9];
    TH1F *Eta_Lep_Btg[7][9];
    TH1F *EtaIso_Lep_Btg[7][9];
    
    TH1F *pt_Jet1_Btg[10][6];
    TH1F *pt_Jet2_Btg[10][6];
    TH1F *dPhi_Jet1_Lep[10][6];
    TH1F *dPhi_Jet2_Lep[10][6];
    TH1F *dR_Jet1_Lep[10][6];
    TH1F *dR_Jet2_Lep[10][6];
    
    TH1F* iso_trk[3];
    TH1F* iso_hcal[3];
    TH1F* iso_ecal[3];
    TH1F* iso_det[3];
    
    TH2F* iso_ecal_hcal[3];
    TH2F* iso_ecal_trk[3];
    TH2F* iso_trk_hcal[3];
    
    TH1F* iso_pf[3];
    TH2F* iso_pf_det[3];
    
    TH1F* pt1[3];
    TH1F* pt2[3];
    TH1F* ptBoth[3];
    
    TH2F* iso_DET_pt1[3];
    TH2F* iso_DET_pt2[3];
    TH2F* iso_DET_ptBoth[3];

    TH2F* iso_PF_pt1[3];
    TH2F* iso_PF_pt2[3];
    TH2F* iso_PF_ptBoth[3];
    
    /*
     TH1F *pt_Btg_Jet1_Flvr_b;
     TH1F *pt_Btg_Jet1_Flvr_c;
     TH1F *pt_Btg_Jet1_Flvr_uds;
     TH1F *pt_Btg_Jet2_Flvr_b;
     TH1F *pt_Btg_Jet2_Flvr_c;
     TH1F *pt_Btg_Jet2_Flvr_uds;
     
     TH1F *Btg_TrkVtxCounting_Jet1;
     TH1F *Btg_TrkVtxCounting_Jet1_b;
     TH1F *Btg_TrkVtxCounting_Jet1_c;
     TH1F *Btg_TrkVtxCounting_Jet1_uds;
     
     TH1F *Btg_TrkVtxCounting_Jet2;
     TH1F *Btg_TrkVtxCounting_Jet2_b;
     TH1F *Btg_TrkVtxCounting_Jet2_c;
     TH1F *Btg_TrkVtxCounting_Jet2_uds;
     
     TH1F *Iso_LFromW[5][9][5];
     TH1F *Iso_LFromTauW[5][9][5];
     TH1F *Iso_LFromB[5][9][5];
     TH1F *Iso_LFromCB[5][9][5];
     TH1F *Iso_LFromCC[5][9][5];
     TH1F *Iso_LFromTauB[5][9][5];
     TH1F *Iso_LFromCW[5][9][5];
     TH1F *Iso_LFromOth[5][9][5];
     TH1F *Iso_LFromFake[5][9][5];
     */
};

#endif
