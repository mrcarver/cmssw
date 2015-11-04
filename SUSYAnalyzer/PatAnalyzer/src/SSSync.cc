#include "SSSync.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/CandAlgos/interface/CandMatcher.h"
#include "PhysicsTools/HepMCCandAlgos/interface/MCTruthPairSelector.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/MergeableCounter.h"


#include "PhysicsTools/CandUtils/interface/CandMatcherNew.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"


using namespace std;
using namespace edm;
using namespace reco;
using namespace tools;
using namespace math;
using namespace reco::tau;

SSSync::SSSync(const edm::ParameterSet & iConfig) :
_relIsoCutE(0.09),
_relIsoCutMu(0.10),
_relIsoCutEloose(999.), //0.6
_relIsoCutMuloose(999.), //1.0
_chargeConsistency(true),
_minPt0(10.),
_minPt1(20.),
_tightD0Mu(0.05),
_tightD0E(0.02),
_looseD0Mu(0.05),
_looseD0E(0.05),
//_looseD0Mu(0.2),
//_looseD0E(9999999.),
//_jetPtCut(40.),
_jetPtCut(18.),
_jetEtaCut(2.4),
_tauPt(20),
_tauEta(2.3),
_regression(false)
{
    Sample              = iConfig.getUntrackedParameter<std::string>("SampleLabel") ;
	doFR 				= iConfig.getUntrackedParameter<bool>("doFR");
	Dov1v3Corrections 	= iConfig.getUntrackedParameter<bool>("Dov1v3Corrections");
    IT_muon             = iConfig.getParameter<edm::InputTag>("MuonLabel") ;
    IT_electron         = iConfig.getParameter<edm::InputTag>("ElectronLabel") ;
    IT_tau              = iConfig.getParameter<edm::InputTag>("TauLabel") ;
    //IT_tauDiscriminator = iConfig.getParameter<edm::InputTag>("TauDiscriminatorLabel") ;
    IT_jet              = iConfig.getParameter<edm::InputTag>("JetLabel");
    IT_pfmet            = iConfig.getParameter<edm::InputTag>("METLabel")  ;
    IT_beamspot         = iConfig.getParameter<edm::InputTag>("BeamSpotLabel");
    IT_hltresults       = iConfig.getParameter<edm::InputTag>("HLTResultsLabel");
    //IT_METFilters       = iConfig.getParameter<edm::InputTag>("METFilter");
	
	MVAidCollection_      = iConfig.getParameter<edm::InputTag>("MVAId");
	
	id_      = iConfig.getParameter<edm::InputTag>("id");
	
	
	consumes<edm::ValueMap<bool> >(id_);
    consumes<edm::ValueMap<unsigned> >(id_);
	consumes<edm::ValueMap<float> >(MVAidCollection_);
    consumes<std::string>(id_);
	

}




void SSSync::beginJob()
{
	
	
	
	
	Triggers = fs->make<TH1F>("Triggers","Triggers",450,0,450);//

    Nvtx           = fs->make<TH1F>("N_{vtx}"        , "Number of vertices;N_{vtx};events / 1"  ,    40, 0., 40.);
    
    _hCounter = fs->make<TH1D>("hCounter", "Events counter", 5,0,5);
	
	const char* flav[2] = {"Ele","Mu"};
	const char* iso[2] = {"","_Iso"};
	const char* types[3] = {"","_FO","_Tight"};
	const char* MTcut[2] = {"","_Inverted"};
	const char* pres[2] = {"","_Prescaled"};
	
	float pbs[6] = {10.0,15.0,25.0,35.0,50.0,70.0};

	

    outputTree = new TTree("SSSyncTree","SSSyncTree");
    _leptonP4 = new TClonesArray("TLorentzVector", 8);
    for (int i=0; i!=8; ++i) {
        new ( (*_leptonP4)[i] ) TLorentzVector();
    }
    //outputTree->Branch("_leptonP4", "TClonesArray", &_leptonP4, 32000, 0);


	
	
	outputTree->Branch("_nL1Mus", &_nL1Mus, "_nL1Mus/I");
	outputTree->Branch("_L1HT", &_L1HT, "_L1HT/D");
	outputTree->Branch("_L1Mu_Pt", &_L1Mu_Pt, "_L1Mu_Pt[40]/D");
	outputTree->Branch("_L1Mu_Eta", &_L1Mu_Eta, "_L1Mu_Eta[40]/D");
	outputTree->Branch("_L1Mu_Phi", &_L1Mu_Phi, "_L1Mu_Phi[40]/D");

	outputTree->Branch("_DiLepHTTrigs", &_DiLepHTTrigs, "_DiLepHTTrigs[3]/I");
	outputTree->Branch("_DiLepTrigs", &_DiLepTrigs, "_DiLepTrigs[5]/I");
	outputTree->Branch("_ControlTrigs", &_ControlTrigs, "_ControlTrigs[2][2][5]/I");
	outputTree->Branch("_ControlTrigPrescale", &_ControlTrigPrescale, "_ControlTrigPrescale[2][2][5]/I");
	outputTree->Branch("_BControlTrigs", &_BControlTrigs, "_BControlTrigs[2]/I");
	outputTree->Branch("_BControlTrigPrescale", &_BControlTrigPrescale, "_BControlTrigPrescale[2]/I");
	
	for(int t=0;t<5;t++){
	
		
		_DiLepObjs[t] = new TClonesArray("TLorentzVector",40);
		for(int i=0;i<40;i++)
			new ( (*_DiLepObjs[t])[i] ) TLorentzVector();
		
		//outputTree->Branch(Form("_DiLepObjs%i",t),"TClonesArray", &_DiLepObjs[t], 32000, 0);
		
		//outputTree->Branch(Form("_nDiLepObjs%i",t), &_nDiLepObjs[t], "_nDiLepObjs[5]/I");
		if(t<3){
			
			_DiLepHTObjs[t] = new TClonesArray("TLorentzVector",40);
			for(int i=0;i<40;i++)
				new ( (*_DiLepHTObjs[t])[i] ) TLorentzVector();
			
			//outputTree->Branch(Form("_DiLepHTObjs%i",t),"TClonesArray", &_DiLepHTObjs[t], 32000, 0);
			
			//outputTree->Branch(Form("_nDiLepHTObjs%i",t), &_nDiLepHTObjs[t], "_nDiLepHTObjs[3]/I");
			if(t<2){
			
				_BControlObjs[t] = new TClonesArray("TLorentzVector",40);
				for(int i=0;i<40;i++)
					new ( (*_BControlObjs[t])[i] ) TLorentzVector();
					
				//outputTree->Branch(Form("_BControlObjs%i",t),"TClonesArray", &_BControlObjs[t], 32000, 0);	
				//outputTree->Branch(Form("_BControlTrigs%i",t), &_BControlTrigs[t], "_BControlTrigs[2]/I");//
				//outputTree->Branch(Form("_nBControlObjs%i",t), &_nBControlObjs[t], "_nBControlObjs[2]/I");//
			}
		}
		
		for(int x=0;x<2;x++){
			for(int y=0;y<2;y++){
				
				_ControlObjs[x][y][t] = new TClonesArray("TLorentzVector",40);
				for(int i=0;i<40;i++)
					new ( (*_ControlObjs[x][y][t])[i] ) TLorentzVector();
				
				//outputTree->Branch(Form("_ControlObjs%i%i%i",x,y,t),"TClonesArray", &_ControlObjs[x][y][t], 32000, 0);
				//outputTree->Branch(Form("_ControlTrigs%i%i%i",x,y,t), &_ControlTrigs[x][y][t], "_ControlTrigs[2][2][5]/I");
				//outputTree->Branch(Form("_nControlObjs%i%i%i",x,y,t), &_nControlObjs[x][y][t], "_nControlObjs[2][2][5]/I");
			}
		}
		

	}	
	

    
    outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
    outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
    outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
	
	
	
	outputTree->Branch("_weight", &_weight , "_weight/D");
	
	
    outputTree->Branch("_nLeptons", &_nLeptons, "_nLeptons/I");///
    outputTree->Branch("_flavors", &_flavors, "_flavors[8]/I");
	outputTree->Branch("_pdgids", &_pdgids, "_pdgids[8]/I");
    outputTree->Branch("_charges", &_charges, "_charges[8]/D");
    outputTree->Branch("_isolation", &_isolation, "_isolation[8]/D");
	outputTree->Branch("_ptrel", &_ptrel, "_ptrel[8]/D");
	outputTree->Branch("_EcalPFIso", &_EcalPFIso, "_EcalPFIso[8]/D");
	outputTree->Branch("_HcalPFIso", &_HcalPFIso, "_HcalPFIso[8]/D");
	outputTree->Branch("_TrackIso", &_TrackIso, "_TrackIso[8]/D");
	outputTree->Branch("_miniIsolation", &_miniIsolation, "_miniIsolation[8]/D");
	outputTree->Branch("_miniIsolationEA", &_miniIsolationEA, "_miniIsolationEA[8]/D");
	outputTree->Branch("_MVAVal", &_MVAVal, "_MVAVal[8]/D");
	outputTree->Branch("_selectedLeptonPt", &_selectedLeptonPt, "_selectedLeptonPt[8]/D");
	outputTree->Branch("_selectedLeptonEta", &_selectedLeptonEta, "_selectedLeptonEta[8]/D");
	outputTree->Branch("_selectedLeptonPhi", &_selectedLeptonPhi, "_selectedLeptonPhi[8]/D");
    outputTree->Branch("_isolationComponents", &_isolationComponents, "_isolationComponents[8][4]/D");
	outputTree->Branch("_inheritance", &_inheritance, "_inheritance[8][30]/I");
    outputTree->Branch("_isolationMC", &_isolationMC, "_isolationMC[8][4]/D");//
    
    outputTree->Branch("_index1", &_index1, "_index1/I");
    outputTree->Branch("_index2", &_index2, "_index2/I");

    outputTree->Branch("_sb", &_sb, "_sb/O");
    outputTree->Branch("_doubleF", &_doubleF, "_doubleF/O");
    
    outputTree->Branch("_origin", &_origin, "_origin[8]/I");
    outputTree->Branch("_originReduced", &_originReduced, "_originReduced[8]/I");
    
    outputTree->Branch("_PVchi2", &_PVchi2, "_PVchi2/D");
    outputTree->Branch("_PVerr", &_PVerr, "_PVerr[3]/D");
    
    outputTree->Branch("_ipPV", &_ipPV, "_ipPV[8]/D");
    outputTree->Branch("_ipPVerr", &_ipPVerr, "_ipPVerr[8]/D");
    outputTree->Branch("_ipZPV", &_ipZPV, "_ipZPV[8]/D");
    outputTree->Branch("_ipZPVerr", &_ipZPVerr, "_ipZPVerr[8]/D");
    
    outputTree->Branch("_ipPVmc", &_ipPVmc, "_ipPVmc[4]/D");
    
    outputTree->Branch("_3dIP", &_3dIP, "_3dIP[8]/D");
    outputTree->Branch("_3dIPerr", &_3dIPerr, "_3dIPerr[8]/D");
    outputTree->Branch("_3dIPsig", &_3dIPsig, "_3dIPsig[8]/D");
    
	outputTree->Branch("_missingHits", &_missingHits, "_missingHits[8]/I");
    outputTree->Branch("_dphi", &_dphi, "_dphi[8]/D");
	outputTree->Branch("_deta", &_deta, "_deta[8]/D");
	outputTree->Branch("_sigIeta", &_sigIeta, "_sigIeta[8]/D");
	outputTree->Branch("_HoE", &_HoE, "_HoE[8]/D");
	outputTree->Branch("_pMe", &_pMe, "_pMe[8]/D");
	
    outputTree->Branch("_mt", &_mt, "_mt[8]/D");
    outputTree->Branch("_isloose", &_isloose, "_isloose[8]/O");
	outputTree->Branch("_isFO", &_isFO, "_isFO[8]/O");
	outputTree->Branch("_isTight", &_isTight, "_isTight[8]/O");
    outputTree->Branch("_istight", &_istight, "_istight[8]/O");
	outputTree->Branch("_istightNoIso", &_istightNoIso, "_istightNoIso[8]/O");
	outputTree->Branch("_istightNoIsoSIP", &_istightNoIsoSIP, "_istightNoIsoSIP[8]/O");
	outputTree->Branch("_islooseMVA", &_islooseMVA, "_islooseMVA[8]/O");
	outputTree->Branch("_istightMVA", &_istightMVA, "_istightMVA[8]/O");
	outputTree->Branch("_istightMVANoIsoSIP", &_istightMVANoIsoSIP, "_istightMVANoIsoSIP[8]/O");
	outputTree->Branch("_istightMVANoIsoSIP_LMVA", &_istightMVANoIsoSIP_LMVA, "_istightMVANoIsoSIP_LMVA[8]/O");
	outputTree->Branch("_chargeConsistent", &_chargeConsistent, "_chargeConsistent[8]/O");
	
    outputTree->Branch("_genHT",&_genHT,"_genHT/D");
	outputTree->Branch("_genCharge",&_genCharge,"_genCharge[8]/I");
    
    outputTree->Branch("_closeJetPtAll", &_closeJetPtAll, "_closeJetPtAll[8]/D");
    outputTree->Branch("_closeJetAngAll", &_closeJetAngAll, "_closeJetAngAll[8]/D");
    outputTree->Branch("_ptRelAll", &_ptRelAll, "_ptRelAll[8]/D");
	outputTree->Branch("_ptRatioAll", &_ptRatioAll, "_ptRatioAll[8]/D");
	outputTree->Branch("_closeIndex", &_closeIndex, "_closeIndex[8]/I");
    
    outputTree->Branch("_closeJetPtAllMC", &_closeJetPtAllMC, "_closeJetPtAllMC[4]/D");
    outputTree->Branch("_closeJetPtAllstatus", &_closeJetPtAllstatus, "_closeJetPtAllstatus[4]/D");
    outputTree->Branch("_partonIdMatched", &_partonIdMatched, "_partonIdMatched[4]/I");
    outputTree->Branch("_sameParton", &_sameParton, "_sameParton[4]/O");
    
    if (_regression) {
        outputTree->Branch("_regVars", &_regVars, "_regVars[15]/D");
        
        outputTree->Branch("hJet_ptRaw", &hJet_ptRaw, "hJet_ptRaw/D");
        outputTree->Branch("hJet_genPt", &hJet_genPt, "hJet_genPt/D");
        outputTree->Branch("hJet_pt", &hJet_pt, "hJet_pt/D");
        outputTree->Branch("hJet_phi", &hJet_phi, "hJet_phi/D");
        outputTree->Branch("hJet_eta", &hJet_eta, "hJet_eta/D");
        outputTree->Branch("hJet_e", &hJet_e, "hJet_e/D");
        
        outputTree->Branch("hJet_ptLeadTrack", &hJet_ptLeadTrack, "hJet_ptLeadTrack/D");
        
        outputTree->Branch("hJet_vtx3dL", &hJet_vtx3dL, "hJet_vtx3dL/D");
        outputTree->Branch("hJet_vtx3deL", &hJet_vtx3deL, "hJet_vtx3deL/D");
        outputTree->Branch("hJet_vtxMass", &hJet_vtxMass, "hJet_vtxMass/D");
        outputTree->Branch("hJet_vtxPt", &hJet_vtxPt, "hJet_vtxPt/D");
        
        outputTree->Branch("hJet_cef", &hJet_cef, "hJet_cef/D");
        
        outputTree->Branch("hJet_nconstituents", &hJet_nconstituents, "hJet_nconstituents/D");
        outputTree->Branch("hJet_JECUnc", &hJet_JECUnc, "hJet_JECUnc/D");
        
        outputTree->Branch("hJet_SoftLeptptRel", &hJet_SoftLeptptRel, "hJet_SoftLeptptRel/D");
        outputTree->Branch("hJet_SoftLeptPt", &hJet_SoftLeptPt, "hJet_SoftLeptPt/D");
        outputTree->Branch("hJet_SoftLeptdR", &hJet_SoftLeptdR, "hJet_SoftLeptdR/D");
        
        outputTree->Branch("hJet_SoftLeptIdlooseMu", &hJet_SoftLeptIdlooseMu, "hJet_SoftLeptIdlooseMu/D");
        outputTree->Branch("hJet_SoftLeptId95", &hJet_SoftLeptId95, "hJet_SoftLeptId95/D");
    }
    
    /*****    outputTree->Branch("_closeJetPt", &_closeJetPt, "_closeJetPt[3]/D");
     outputTree->Branch("_closeJetAng", &_closeJetAng, "_closeJetAng[3]/D");
     outputTree->Branch("_ptRel", &_ptRel, "_ptRel[3]/D");
     
     outputTree->Branch("_ptRelAllNotItself", &_ptRelAllNotItself, "_ptRelAllNotItself[3]/D");
     outputTree->Branch("_closeJetPtAllNotItself", &_closeJetPtAllNotItself, "_closeJetPtAllNotItself[3]/D");
     outputTree->Branch("_closeJetAngAllNotItself", &_closeJetAngAllNotItself, "_closeJetAngAllNotItself[3]/D");
     */
   
   
   
     outputTree->Branch("_n_bJetsAll", &_n_bJetsAll, "_n_bJetsAll/I");
     outputTree->Branch("_n_JetsAll", &_n_JetsAll, "_n_JetsAll/I");
     outputTree->Branch("_bTaggedAll", &_bTaggedAll, "_bTaggedAll[100]/O");
     outputTree->Branch("_jetEtaAll", &_jetEtaAll, "_jetEtaAll[100]/D");
     outputTree->Branch("_jetPhiAll", &_jetPhiAll, "_jetPhiAll[100]/D");
     outputTree->Branch("_jetPtAll", &_jetPtAll, "_jetPtAll[100]/D");
     outputTree->Branch("_csvAll", &_csvAll, "_csvAll[100]/D");
     
    outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
    
	outputTree->Branch("_rho", &_rho, "_rho/D");
	outputTree->Branch("_rhoNC", &_rhoNC, "_rhoNC/D");
    outputTree->Branch("_met", &_met, "_met/D");
    outputTree->Branch("_met_phi", &_met_phi, "_met_phi/D");
    outputTree->Branch("HT", &HT, "HT/D");
    
    
    outputTree->Branch("_mompt", &_mompt, "_mompt[4]/D");
    outputTree->Branch("_momphi", &_momphi, "_momphi[4]/D");
    outputTree->Branch("_mometa", &_mometa, "_mometa[4]/D");
    outputTree->Branch("_mompdg", &_mompdg, "_mompdg[4]/I");
    
    outputTree->Branch("_n_bJets", &_n_bJets, "_n_bJets/I");
    outputTree->Branch("_n_Jets", &_n_Jets, "_n_Jets/I");
    outputTree->Branch("_bTagged", &_bTagged, "_bTagged[20]/O");
    outputTree->Branch("_jetEta", &_jetEta, "_jetEta[20]/D");
    outputTree->Branch("_jetPhi", &_jetPhi, "_jetPhi[20]/D");
    outputTree->Branch("_jetPt", &_jetPt, "_jetPt[20]/D");
    outputTree->Branch("_csv", &_csv, "_csv[20]/D");
    

    GPM = GenParticleManager();
	
	
    
    bool isData = !(Sample=="ElectronsMC");
	
	std::string GT = "Summer15_25nsV2_MC";
	if(isData)
		GT = "Summer15_25nsV5_DATA";
	
    fMetCorrector = new OnTheFlyCorrections(GT, isData); //isData = true
	
    _corrLevel = "L1FastJet";
	_corrLevelAbs = "L3Absolute";
	_corrLevelRes = "L2L3Residual";
	_corrLevelFull = "L3Absolute";
	
    if (isData) _corrLevelFull = "L2L3Residual";
    
    _nEventsTotal = 0;
    _nEventsFiltered = 0;
    _nEventsTotalCounted = 0;
	
	
	std::cout<<"begin end\n";
	
	calib = new BTagCalibration("csvv2","/lfs/scratch/mrcarver/Fall15AnalysisCode/CMSSW_7_4_14/src/SUSYAnalyzer/PatAnalyzer/test/CSVv2.csv");
	reader = new BTagCalibrationReader(calib,     // calibration instance
                             BTagEntry::OP_MEDIUM,  // operating point
                             "mujets",               // measurement type
                             "central");//
							 
	
	//TFile *file = new TFile("/lfs/scratch/mrcarver/Fall15AnalysisCode/CMSSW_7_4_14/src/SUSYAnalyzer/PatAnalyzer/test/btageff__ttbar_powheg_pythia8_25ns.root","read");
	//btagEff = new TH2F("btagEff","",15,0,15,7,0,2.6);
	//btagEff->Read("h2_BTaggingEff_csv_med_Eff_c");
	//file->Close();
	
}

void SSSync::endJob() {
    //outputTree -> Write();
    // store nEventsTotal and nEventsFiltered in preferred way
    //std::cout<<_nEventsTotalCounted<<std::endl;
    //std::cout<<_nEventsFiltered<<std::endl;
    
    delete fMetCorrector;
    
    
}

void SSSync::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
  using namespace edm;

  bool changed(false);
  //edm::InputTag *triggerEventTag_  = new edm::InputTag("TriggerResults");
  //if(!hltConfig_.init(iRun,iSetup,triggerEventTag_->process(),changed) ){
  if(!hltConfig_.init(iRun,iSetup,"HLT",changed) ){
	edm::LogError( "CandidateTriggerObjectProducer" ) <<
	  "Error! Can't initialize HLTConfigProvider";
	throw cms::Exception("HLTConfigProvider::init() returned non 0");
  }

  return;

}

void SSSync::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{

	
	//std::cout<<"ana1\n";

	bool debug = false;

	if(debug) std::cout<<"0\n";
	
	
	

	
	edm::Handle<edm::TriggerResults> trigResults;
	edm::InputTag trigResultsTag("TriggerResults","","HLT");
	iEvent.getByLabel(trigResultsTag,trigResults);
	/*
	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	edm::InputTag triggerObjects_("selectedPatTrigger","","PAT");//PAT for MC and RECO for data
	edm::InputTag triggerObjectsData_("selectedPatTrigger","","RECO");//PAT for MC and RECO for data
	if(Sample=="ElectronsMC")
		iEvent.getByLabel(triggerObjects_, triggerObjects);
	else
		iEvent.getByLabel(triggerObjectsData_, triggerObjects);
	*/
	const edm::TriggerNames trigNames = iEvent.triggerNames(*trigResults);
	
	
	
	
	 //size_t filterIndex = triggerSummary->filterIndex( triggerFilterEle_ );
	// trigger::TriggerObjectCollection triggerObjects = triggerSummary->getObjects();
	 //if( !(filterIndex >= triggerSummary->sizeFilters()) ){
	//     const trigger::Keys& keys = triggerSummary->filterKeys( filterIndex );
	//     for( size_t j = 0; j < keys.size(); ++j ){
	//         trigger::TriggerObject foundObject = triggerObjects[keys[j]];
	//         if(fabs(foundObject.id()) == 11){ //It's an electron

	//}
	//}
	//}       
	
	//std::cout<<"trigresults size = "<<trigResults->size()<<"\n";//
	//unsigned int tsize = trigResults->size();
	//tsize--;
	//tsize--;
	//std::cout << "\n === TRIGGER PATHS === " << std::endl;
	
	if(debug) std::cout<<"0.1\n";
	
	for(int z=0;z<3;z++){
		_DiLepHTTrigs[z] = 0;
		//std::cout<<"dilepHT names "<<z<<" = "<<DiLepHTTrigNames[z]<<"\n";	
	}
	
	for(int z=0;z<2;z++){
		_BControlTrigs[z] = 0;
		_BControlTrigPrescale[z] = 1;
	}
	
	for(int z=0;z<5;z++)
		_DiLepTrigs[z] = 0;
	
	for(int z=0;z<2;z++){
		for(int x=0;x<2;x++){
			for(int y=0;y<5;y++){
				_ControlTrigs[z][x][y] = 0;
				_ControlTrigPrescale[z][x][y] = 1;
			}
		}
	}
	
	if(debug) std::cout<<"0.2\n";
	
	//HLTConfigProvider *hltconfig = new HLTConfigProvider;////
	//hltconfig->init(
	
	//edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesH_; 
	//edm::InputTag triggerPrescaleToken("patTrigger","","PAT");
  	//iEvent.getByLabel( triggerPrescaleToken, triggerPrescalesH_);
	
    for (unsigned int i = 0, n = trigResults->size(); i < n; ++i) {
	//for (unsigned int i = 0; i< 300; i++) {
     
	
	// const char *path = trigNames.triggerName(i).c_str();
	//  Triggers->GetXaxis()->SetBinLabel(i+1,path);

	  if(trigResults->accept(i)){
	  
	  	TString PathName = trigNames.triggerName(i);
		//std::cout << "\n\nTrigger " << trigNames.triggerName(i) << 
    	//  ": " << (trigResults->accept(i) ? "PASS" : "fail (or not run)") 
    	 // << std::endl;
	  		/*
	    //get prescale info from hltConfig_
		//std::vector<int> *prescales, *l1prescales;
		if(debug) std::cout<<"0.2.0.0\n";
		std::pair<std::vector<std::pair<std::string,int> >,int> detailedPrescaleInfo = hltConfig_.prescaleValuesInDetail(iEvent, iEventSetup, trigNames.triggerName(i));   
		if(debug) std::cout<<"0.2.0\n";
		//prescales->push_back( triggerPrescalesH_.isValid() ? detailedPrescaleInfo.second : -1 );
		//prescales->push_back( detailedPrescaleInfo.second  );
	  //	std::cout<<" HLT prescale = "<<detailedPrescaleInfo.second<<"\n";
	  
		// save l1 prescale values in standalone vector
		std::vector <int> l1prescalevals;
		for( size_t varind = 0; varind < detailedPrescaleInfo.first.size(); varind++ ){
		  l1prescalevals.push_back(detailedPrescaleInfo.first.at(varind).second);
		}
		if(debug) std::cout<<"0.2.0.1.1\n";
		// find and save minimum l1 prescale of any ORed L1 that seeds the HLT
		std::vector<int>::iterator result = std::min_element(std::begin(l1prescalevals), std::end(l1prescalevals));
		size_t minind = std::distance(std::begin(l1prescalevals), result);
		// sometimes there are no L1s associated with a HLT. In that case, this branch stores -1 for the l1prescale
		//l1prescales->push_back( minind < l1prescalevals.size() ? l1prescalevals.at(minind) : -1 );
		*/
		//std::cout<<" L1 prescale = "<<l1prescalevals.at(minind)<<"\n";
	      
	  if(debug) std::cout<<PathName<<"\n";
	  	//int Prescale = 1;
		//if(triggerPrescalesH_.isValid() && minind)
		//int	Prescale = detailedPrescaleInfo.second * l1prescalevals.at(minind);
		//if(debug) std::cout<<"0.2.0.2\n";
		
		//std::cout<<"Prescale = "<<Prescale<<"\n";
		//std::cout<<PathName<<" passed!\n";
		//std::pair<int,int> prescale = hltConfig_.prescaleValues(iEvent,iEventSetup,trigNames.triggerName(i));
		//std::cout<<PathName<<" L1 prescale = "<<prescale.first<<" and HLT = "<<prescale.second<<"\n";
	  	for(int z=0;z<3;z++){
			
			
			if(PathName.Contains(DiLepHTTrigNames[z])){
				_DiLepHTTrigs[z] = 1;
				//std::cout<<DiLepHTTrigNames[z]<<" Passed!\n";
			}
				
		}
		for(int z=0;z<2;z++){
			if(PathName.Contains(BControlTrigNames[z])){
				_BControlTrigs[z] = 1;
				_BControlTrigPrescale[z] = 1;//detailedPrescaleInfo.second * l1prescalevals.at(minind);//Prescale;
			}
		}
		for(int z=0;z<5;z++){
			if(PathName.Contains(DiLepTrigNames[z]))
				_DiLepTrigs[z] = 1;
		}
		for(int z=0;z<2;z++){
			for(int x=0;x<2;x++){
				for(int y=0;y<5;y++){
					if(PathName.Contains(ControlTrigNames[z][x][y])){
						_ControlTrigs[z][x][y] = 1;
						_ControlTrigPrescale[z][x][y] = 1;//detailedPrescaleInfo.second * l1prescalevals.at(minind);//Prescale;
						
					}
				}
			}
		}
	  
	//  Triggers->Fill(i);
	  
	  }
	
    }
	
	
	/*
	std::vector<TLorentzVector> MBmuon, EBeles;
	
	edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigObjs;
	iEvent.getByLabel("selectedPatTrigger","",trigObjs);
	for(std::vector<pat::TriggerObjectStandAlone>::const_iterator to = trigObjs->begin();to != trigObjs->end();to++){
		
		TLorentzVector obj;
		obj.SetPtEtaPhiM(to->pt(),to->eta(),to->phi(),to->energy());
		if(to->hasFilterLabel("hltL3fL1sMu0L1f0L2f3QL3Filtered10Q"))
			MBmuon.push_back(obj);
		
		
		if(to->hasFilterLabel("hltSingleEle10CaloIdMTrackIdMDphiFilter"))
			EBeles.push_back(obj);
			
		//if(_BControlTrigs[1]){
		//std::vector<std::string> st = to->filterLabels();
		//for(std::vector<std::string>::iterator sti = st.begin();sti != st.end();sti++)
		///	std::cout<<*sti<<"\n";
		//}
	}
	
	*/
	
	if(debug) std::cout<<"0.4\n";
	//============ Primary vertices ============
    edm::InputTag IT_goodVtx = edm::InputTag("goodOfflinePrimaryVertices");
    edm::Handle<std::vector<Vertex> > theVertices;
    iEvent.getByLabel( "goodOfflinePrimaryVertices", theVertices) ;
    if( ! theVertices.isValid() ) ERR(IT_goodVtx ) ;
    int nvertex = theVertices->size();
    
    _n_PV = nvertex;
    
    Nvtx->Fill(TMath::Min(nvertex,39));//
    if(! nvertex ){
        cout << "[WARNING]: No candidate primary vertices passed the quality cuts, so skipping event" << endl;
        return ;
    }
	

	_PFHT800 = 0;
	
	
	/*
	std::string ht800 = (Sample == "ElectronsMC" ? "HLT_PFHT900_v1" : "HLT_PFHT800_v1");
	bool tht800 = trigResults->accept(trigNames.triggerIndex(ht800));
	if(tht800)
		_PFHT800 = 1;
	*/

	if(debug) std::cout<<"1.1\n";
	
	int triglepcounter[8] = {0,0,0,0,0,0,0,0};
	int mbtlc = 0, ebtlc = 0;
	int ConLepCounter[4][5] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
	
	for(int i=0;i<5;i++){
	
		_nDiLepObjs[i] = 0;
		if(i<3)
			_nDiLepHTObjs[i] = 0;
		
		if(i<2)
			_nBControlObjs[i] = 0;
	
		for(int j=0;j<2;j++){
			for(int k=0;k<2;k++)
				_nControlObjs[j][k][i] = 0;
			
		}
	}
	
	for(int i=0;i<8;i++){
		for(int k=0;k<30;k++){
			_inheritance[i][k] = 0;
		}
	}

	/*
	if(debug) std::cout<<"1.2\n";
	for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      
	    obj.unpackPathNames(trigNames);
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);
   
		for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
		
			//std::cout<<"pathNamesAll "<<h<<" "<<pathNamesAll[h]<<"\n";
           // bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
            bool isL3   = obj.hasPathName( pathNamesAll[h], true, true ); 
            //bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
            //bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
			
			TString pn = pathNamesAll[h];
			//if(pn) std::cout<<"pathNamesAll "<<h<<" "<<pn<<"\n";
			for(int z=0;z<5;z++){
			
				if(pn.Contains(DiLepTrigNames[z]) && _DiLepTrigs[z] && isL3 && _nDiLepObjs[z] < 40){
					_nDiLepObjs[z]++;
					((TLorentzVector *)_DiLepObjs[z]->At(_nDiLepObjs[z]))->SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
				} 
				
				if(z<3){
				
					if(pn.Contains(DiLepHTTrigNames[z]) && _DiLepHTTrigs[z] && isL3 && _nDiLepHTObjs[z] < 40){
						_nDiLepHTObjs[z]++;
						((TLorentzVector *)_DiLepHTObjs[z]->At(_nDiLepHTObjs[z]))->SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
					} 
					
					if(z<2){
						if(pn.Contains(BControlTrigNames[z]) && _BControlTrigs[z] && isL3 && _nBControlObjs[z] < 40){
							_nBControlObjs[z]++;
							((TLorentzVector *)_BControlObjs[z]->At(_nBControlObjs[z]))->SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
						} 
					}
				}
				
				for(int x=0;x<2;x++){
					for(int y=0;y<2;y++){
						if(pn.Contains(ControlTrigNames[x][y][z]) && _ControlTrigs[x][y][z] && isL3 && _nControlObjs[x][y][z] < 40){
							_nControlObjs[x][y][z]++;
							((TLorentzVector *)_ControlObjs[x][y][z]->At(_nControlObjs[x][y][z]))->SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
						} 
					}
				}
			}
			
			
			
			if(pathNamesAll[h] == ht800 && tht800 && isL3)
				_hltHT = obj.pt();
		
			
			
        }
	}
	*/

	float genqpt = 0.0;
	
    edm::Handle<GenParticleCollection> TheGenParticles;
    //bool islepton;
    if (Sample=="ElectronsMC") {
        //******************************************************************************************************************
        // Gen level particles                  ****************************************************************************
        //******************************************************************************************************************
        iEvent.getByLabel("prunedGenParticles", TheGenParticles);
        std::vector<const GenParticle*> vGenElectrons, vGenMuons, vGenNPElectrons, vGenNPMuons, vGenW;
        if( TheGenParticles.isValid() )
        {
            GPM.SetCollection(TheGenParticles);
            GPM.Classify();
            vGenMuons = GPM.filterByStatus(GPM.getPromptMuons(),1);
            vGenElectrons = GPM.filterByStatus(GPM.getPromptElectrons(),1);
            vGenNPMuons = GPM.filterByStatus(GPM.getNonPromptMuons(),1);
            vGenNPElectrons = GPM.filterByStatus(GPM.getNonPromptElectrons(),1);
            //std::cout<<"*************"<<std::endl;
			
			
			
			for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ )
            {
			
				int id = TMath::Abs(p->pdgId());
				if ((id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 21 || id == 22 ) && (p->status() == 23)){
                      genqpt += p->pt();
                }
			
			}
        }
        //******************************************************************************************************************
        //******************************************************************************************************************
        //******************************************************************************************************************
        
        //**************************************************************************************
        // MC//
        //**************************************************************************************
    }
	
	_genHT = genqpt;
    
    _runNb = iEvent.id().run();
    _eventNb = iEvent.id().event();
    _lumiBlock = iEvent.luminosityBlock();
	

   if(debug) std::cout<<"1.4\n";
	
	
    //============ Counter done ============
    
    //============ Beamspot ============
    edm::Handle< reco::BeamSpot > theBeamSpot;
    iEvent.getByLabel( IT_beamspot, theBeamSpot );
    if( ! theBeamSpot.isValid() ) ERR( IT_beamspot ) ;
    BeamSpot::Point  BS= theBeamSpot->position();
    //==================================
    
	/*
	//=========== L1 Particles =============
	edm::Handle<std::vector<l1extra::L1EtMissParticle>> L1HTContainer;
	iEvent.getByLabel("l1extraParticles","MHT",L1HTContainer);
	
	edm::Handle<std::vector<l1extra::L1MuonParticle>> L1MuContainer;
	iEvent.getByLabel("l1extraParticles","",L1MuContainer);
		
	for(std::vector<l1extra::L1MuonParticle>::const_iterator nm = L1MuContainer->begin();nm != L1MuContainer->end();nm++){
	
		
		if(_nL1Mus > 40) continue;
		
		
		_L1Mu_Pt[_nL1Mus] = nm->pt();
		_L1Mu_Eta[_nL1Mus] = nm->eta();
		_L1Mu_Phi[_nL1Mus] = nm->phi();
		
	
	
		_nL1Mus++;
	}
	//======================================
   */
   if(debug) std::cout<<"1.5\n";
    
    Vertex::Point PV = theVertices->begin()->position();
    const Vertex* PVtx = &((*theVertices)[0]);
    _PVchi2 = PVtx->chi2();
    _PVerr[0] = PVtx->xError();
    _PVerr[1] = PVtx->yError();
    _PVerr[2] = PVtx->zError();
    //==================================
    
    //============ Pat MET ============
    edm::Handle< vector<pat::MET> > ThePFMET;
    iEvent.getByLabel(IT_pfmet, ThePFMET);
    if( ! ThePFMET.isValid() ) ERR( IT_pfmet );
    const vector<pat::MET> *pfmetcol = ThePFMET.product();
    const pat::MET *pfmet;
    pfmet = &(pfmetcol->front());
	
    _met = pfmet->pt();
    _met_phi = pfmet->phi();
	//double rawmetX = pfmet->px();
    //double rawmetY = pfmet->py();
    //==================================
    
    //============ Pat Muons ============
    edm::Handle< std::vector<pat::Muon> > thePatMuons;
    iEvent.getByLabel( IT_muon, thePatMuons );
    if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
    //==================================
    
    //============ Pat Electrons ============
    edm::Handle< std::vector<pat::Electron> > thePatElectrons;
    iEvent.getByLabel( IT_electron, thePatElectrons );
    if( ! thePatElectrons.isValid() ) ERR( IT_electron );
    //==================================
    
    //============ Conversions ============
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    iEvent.getByLabel("reducedEgamma","reducedConversions", theConversions);
    //==================================
    if(_eventNb == 110365493)  std::cout<<"6\n";
    
    //============ Pat Jets ============
    edm::Handle< std::vector< pat::Jet> > thePatJets;
    iEvent.getByLabel(IT_jet , thePatJets );
    if( ! thePatJets.isValid() ) ERR(IT_jet);
    //==================================
    
   
    //============ Rho ============
    edm::Handle<double> rhoJets;
    iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetAll","") , rhoJets);//kt6PFJets////
    double myRhoJets = *rhoJets;
	_rho = myRhoJets;
	
	edm::Handle<double> rhoJetsNC;
    iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetCentralNeutral","") , rhoJetsNC);//kt6PFJets////fixedGridRhoFastjetAll
    double myRhoJetsNC = *rhoJetsNC;
	_rhoNC = myRhoJetsNC;
    //==================================
	
	bool FUCK = false;
	bool PD = false;
	if(Dov1v3Corrections){
	
	
	//in miniAOD v1
    //double rawmetX = pfmet->shiftedP2_74x(pat::MET::METUncertainty(12),pat::MET::Raw).px;
    //double rawmetY = pfmet->shiftedP2_74x(pat::MET::METUncertainty(12),pat::MET::Raw).py;

    //in miniAOD v2
    double rawmetX = pfmet->uncorPx(); 
    double rawmetY = pfmet->uncorPy();

    
    double corrMetx = rawmetX;
    double corrMety = rawmetY;

    
    for( std::vector<pat::Jet>::const_iterator jet = (*thePatJets).begin(); jet != (*thePatJets).end(); jet++ ) {
        TLorentzVector pJet;
        pJet.SetPtEtaPhiE(((&*jet)->correctedP4("Uncorrected")).Pt(), ((&*jet)->correctedP4("Uncorrected")).Eta(),((&*jet)->correctedP4("Uncorrected")).Phi(), ((&*jet)->correctedP4("Uncorrected")).E());
        for(int unsigned it=0; it!=(&*jet)->numberOfDaughters(); ++it) {
            const pat::PackedCandidate * icand = dynamic_cast<const pat::PackedCandidate *> ((&*jet)->daughter(it));
            if(  fabs(icand->pdgId()) == 13 ) {
                double mupt = icand->pt();
                double muphi = icand->phi();
                for( std::vector<pat::Muon>::const_iterator mu = (*thePatMuons).begin() ; mu != (*thePatMuons).end() ; mu++ ) {
                    if ( (fabs(mu->pt() - mupt)/mupt < 0.1) && (fabs(mu->phi() - muphi) < 0.01)  ) {
                        if (mu->isGlobalMuon() || mu->isStandAloneMuon()) {
                            TLorentzVector pMu;
                            pMu.SetPtEtaPhiE(mu->pt(),mu->eta(),mu->phi(),mu->energy());
                            //pMu.SetPtEtaPhiE(icand->pt(),icand->eta(),icand->phi(),icand->energy());
                            pJet-=pMu;
                        }
                    }
                }
            }
        }
        std::pair <double, double> corr = fMetCorrector->getCorrections(
                                                                      ((&*jet)->correctedP4("Uncorrected")).Pt(),
                                                                      ((&*jet)->correctedP4("Uncorrected")).Eta(),
                                                                      pJet.Pt(),
                                                                      pJet.Phi(),
                                                                      (&*jet)->neutralEmEnergyFraction()+(&*jet)->chargedEmEnergyFraction(),
                                                                      myRhoJets,
                                                                      (&*jet)->jetArea());
        corrMetx += corr.first ;
        corrMety += corr.second;
    }
    double newmet    = sqrt(corrMetx*corrMetx + corrMety*corrMety);
    double newmetphi = atan2(corrMety, corrMetx);

    
    //std::cout<<newmet<<" "<<_met<<"; "<<newmetphi<<" "<<_met_phi<<" "<<std::endl;
    //std::cout<<std::endl;

    
    _met = newmet;
    _met_phi = newmetphi;
	
	
	
	
	}

    //============= 3D IP ==============
    //ESHandle<TransientTrackBuilder> theTTBuilder;
    //iEventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);//
    //==================================
	
	//============Packed PF Cands for MiniIso=========
	edm::Handle<pat::PackedCandidateCollection> pfcands;
    iEvent.getByLabel("packedPFCandidates", pfcands);
	//================================================
	
	
	//============= MVA Electron ID =================
	edm::Handle<edm::ValueMap<float> > valuesMap;
    iEvent.getByLabel(MVAidCollection_,valuesMap);
	//===============================================
	
	
	_weight = 1;
	if(Sample=="ElectronsMC"){//============Individual Event Weight=============
	edm::Handle<GenEventInfoProduct> pdfvariables;
	iEvent.getByLabel("generator",pdfvariables);
	_weight = pdfvariables->weight();
	
	if(_weight > 0)
		_weight = 1;
	else
		_weight = -1;
		
	int aindex = _n_PV;
	if(_n_PV > 59)
		aindex = 59;
		
	_weight *= PUarray[aindex];
    
	}
	
	
	 //============ Total number of events is the sum of the events ============
    //============ in each of these luminosity blocks ============
    if(Sample=="ElectronsMC"){
		_nEventsTotalCounted += _weight;
		_hCounter->Fill(0.,_weight);
	}
	else{
		_nEventsTotalCounted++;
    	_hCounter->Fill(0);
	}

    

	std::vector<const pat::Muon* > sMu = MVALooseMuonSelector( *thePatMuons, _minPt0, PV, _looseD0Mu);
  
    std::vector<const pat::Electron* > sEl = MVALooseElectronSelector( thePatElectrons, _minPt0, PV, _looseD0E, _chargeConsistency, theConversions, BS,*valuesMap);//, myMVATrig);
	
	//for(auto ele=thePatElectrons->begin();ele!=thePatElectrons->end();++ele){
   //   const edm::Ptr<pat::Electron> elePtr(thePatElectrons,ele-thePatElectrons->begin()); //value map is keyed of edm::Ptrs so we need to make one
	//  float value = (*valuesMap)[elePtr];
	//  std::cout<<"MVA output = "<<value<<"\n";
   // }

    
    std::vector<const pat::Jet* > SelectedJetsAll = JetSelectorAll(*thePatJets, 5., 3.0);
    
    std::vector<const pat::Jet* > SelectedJets = JetSelector(*thePatJets, _jetPtCut, _jetEtaCut);//_jetEtaCut

	if(debug) std::cout<<"2\n";


    if (sEl.size() + sMu.size() < 2 && !doFR) return;//two lepton selection
	
	if(sEl.size() + sMu.size() < 1 ) return;//FR selection//
    
    HT = 0.;
    std::vector< const pat::Jet* > Bjets;
	
    _n_bJetsAll = 0;
    
    int n_bJetsAll30 = 0;
    
	
    _n_Jets = 0;
    _n_bJets = 0;
   // std::cout<<"before jets\n";
   if(debug) std::cout<<"3\n";
 
    for(unsigned int i = 0 ; i < SelectedJetsAll.size() ;i++ ){
        
        
        
		//std::cout<<"before L1\n";
		
       
		
		//std::cout<<"before L3abs\n";
		
		
            
        _jetEtaAll[i] = SelectedJetsAll[i]->eta();//uncEta;
        _jetPhiAll[i] = SelectedJetsAll[i]->phi();//uncPhi;
        _jetPtAll[i] = SelectedJetsAll[i]->pt();//uncPt;//*corr;
		
		if(debug) std::cout<<"3.1\n";
		if(Dov1v3Corrections){
		
			double uncPt = (SelectedJetsAll[i]->correctedP4("Uncorrected")).Pt();
        	double uncEta = (SelectedJetsAll[i]->correctedP4("Uncorrected")).Eta();
        	//double uncPhi = (SelectedJetsAll[i]->correctedP4("Uncorrected")).Phi();
			
			double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[i]->jetArea(),_corrLevel);
			double L2L3corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[i]->jetArea(),_corrLevelFull);
			_jetPtAll[i] = uncPt*L2L3corr;
		
		}
		if(debug) std::cout<<"3.2\n";
		
		//if(PD) std::cout<<" jet "<<i<<" raw pT = "<<(SelectedJetsAll[i]->correctedP4("Uncorrected")).Pt()<<", corr pT = "<<_jetPtAll[i]<<", eta = "<<_jetEtaAll[i]<<", phi = "<<_jetPhiAll[i]<<"\n";//, L1corr = "<<corr<<", L2L3 = "<<Lcorr<<", rho = "<<myRhoJets<<", area = "<<SelectedJetsAll[i]->jetArea()<<"\n";
        
        //((TLorentzVector *)_jetAllP4->At(i))->SetPtEtaPhiM( _jetPtAll[i], _jetEtaAll[i], _jetPhiAll[i], 0 );
        
        _csvAll[i] = SelectedJetsAll[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        
        if(SelectedJetsAll[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.814) {
            _bTaggedAll[i] = true;
            _n_bJetsAll++;
            if (_jetPtAll[i] > _jetPtCut) {
                n_bJetsAll30++;
            }
        } else _bTaggedAll[i] = false;
        
        
    }
    _n_JetsAll = SelectedJetsAll.size();
	
  
    int _sameSign[2][2] = {{0, 0}, {0, 0}};
    int leptonCounter = 0;
	
	if(debug) std::cout<<"3.9\n";
    for(unsigned int i = 0 ; i < sMu.size() ;i++ ){
        
        
        const pat::Muon *iM = sMu[i];
        _leptonIndex = i;
        
        if (leptonCounter == 8) continue;

        _flavors[leptonCounter] = 1;
		_pdgids[leptonCounter] = iM->pdgId();
        _charges[leptonCounter] = iM->charge();
        _isolation[leptonCounter] = pfRelIso15(iM,myRhoJets);
        _miniIsolation[leptonCounter] = getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, false,myRhoJetsNC);
		_TrackIso[leptonCounter] = 0.0;
		_EcalPFIso[leptonCounter] = 0.0;
		_HcalPFIso[leptonCounter] = 0.0;
		
        _sameSign[int(_isolation[leptonCounter] > 0.1)][int((_charges[leptonCounter]+1)/2)]++;

        
        double chargedHadronIso = iM->pfIsolationR03().sumChargedHadronPt;
        double neutralHadronIso = iM->pfIsolationR03().sumNeutralHadronEt;
        double photonIso = iM->pfIsolationR03().sumPhotonEt;
        
		if(debug) std::cout<<"3.9.1\n";
        double Aeff[5] = { 0.0913, 0.0765, 0.0546, 0.0728, 0.1177 };
		double Aeff_Fall15Anal[5] = { 0.0735, 0.0619, 0.0465, 0.0433, 0.0577 };
    	double CorrectedTerm=0.0;
		double miniR = TMath::Max(0.05,TMath::Min(0.2,10.0/iM->pt()));
		double mR2 = miniR*miniR;
    
    	if( TMath::Abs( iM->eta() ) < 0.8 ) CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 0 ]*(mR2/0.09);
    	else if( TMath::Abs( iM->eta() ) > 0.8 && TMath::Abs( iM->eta() ) < 1.3  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 1 ]*(mR2/0.09);
    	else if( TMath::Abs( iM->eta() ) > 1.3 && TMath::Abs( iM->eta() ) < 2.0  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 2 ]*(mR2/0.09);//
   	    else if( TMath::Abs( iM->eta() ) > 2.0 && TMath::Abs( iM->eta() ) < 2.2  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 3 ]*(mR2/0.09);
   	    else if( TMath::Abs( iM->eta() ) > 2.2 && TMath::Abs( iM->eta() ) < 2.5  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 4 ]*(mR2/0.09);
      
        _isolationComponents[leptonCounter][0] = chargedHadronIso/iM->pt();
        _isolationComponents[leptonCounter][1] = neutralHadronIso/iM->pt();
        _isolationComponents[leptonCounter][2] = photonIso/iM->pt();
        _isolationComponents[leptonCounter][3] = CorrectedTerm/iM->pt();
		
		_miniIsolationEA[leptonCounter] = TMath::Max(0.0,_miniIsolation[leptonCounter] - CorrectedTerm/iM->pt());
		_MVAVal[leptonCounter] = 10.0;
		_chargeConsistent[leptonCounter] = true;
		
		if(iM->genParticle())
			_genCharge[leptonCounter] = iM->genParticle()->pdgId();
		else
			_genCharge[leptonCounter] = 0;
	
		std::vector<const reco::Candidate*> moms;
		if(iM->genParticle() && iM->genParticle()->mother(0)){
		
		
		
		moms.push_back(iM->genParticle()->mother(0));
		for(int i=0;i<30;i++){
			if(!moms[i] || !moms[i]->mother(0))
			break;
			
			moms.push_back(	moms[i]->mother(0) );
		}
		}
		for(unsigned int j=0;j<moms.size();j++){
			_inheritance[leptonCounter][j] = moms[j]->pdgId();
		}
		
		if(debug) std::cout<<"3.9.2\n";
		
        _ipPV[leptonCounter] = TMath::Abs(iM->innerTrack()->dxy(PV));
		if(debug) std::cout<<"3.9.2.1\n";
        _ipPVerr[leptonCounter] = iM->innerTrack()->dxyError();
        if(debug) std::cout<<"3.9.2.2\n";
        _ipZPV[leptonCounter] = iM->innerTrack()->dz(PV);
		if(debug) std::cout<<"3.9.2.3\n";
        _ipZPVerr[leptonCounter] = iM->innerTrack()->dzError();
		if(debug) std::cout<<"3.9.2.4\n";
		
		const reco::GenParticle *mc = iM->genParticle();
		if(debug) std::cout<<"3.9.2.5\n";
		_origin[leptonCounter] = -1;
		//if(mc && iM->genParticle())
		//	_origin[leptonCounter] = GPM.origin(mc);
		
		if(debug) std::cout<<"3.9.2\n";
			
		_originReduced[leptonCounter] = GPM.originReduced(_origin[leptonCounter]);
        
  		if(debug) std::cout<<"3.9.3\n";
		_missingHits[leptonCounter] = -1;
		_dphi[leptonCounter] = -1.;
		_deta[leptonCounter] = -1.;
		_sigIeta[leptonCounter] = -1.;
		_HoE[leptonCounter] = -1.;
		_pMe[leptonCounter] = -1.;


            _3dIP[leptonCounter]    = iM->dB(pat::Muon::PV3D);// ip3D.value();
            _3dIPerr[leptonCounter] = iM->edB(pat::Muon::PV3D);//ip3D.error();
            _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);//ip3D.significance();
   
   
		_isloose[leptonCounter] = (_isolation[leptonCounter]< 0.5 && _ipPV[leptonCounter] < 0.05 && _ipZPV[leptonCounter] < 0.1);//
		_islooseMVA[leptonCounter] = (_isolation[leptonCounter]< 0.5 && _ipPV[leptonCounter] < 0.05 && _ipZPV[leptonCounter] < 0.1);//
		
		
	
		
        if (_isloose[leptonCounter])
        	_istight[leptonCounter] = ( tools::SSSyncIsTight(*iM, _minPt0, *theVertices->begin(), _looseD0Mu, _3dIPsig[leptonCounter]) ) && (_isolation[leptonCounter]< 0.1);
        else _istight[leptonCounter] = false;
		
		if( tools::SSSyncIsTight(*iM, _minPt0, *theVertices->begin(), _looseD0Mu, _3dIPsig[leptonCounter]) ) 
        	_istightNoIso[leptonCounter] = true;
        else _istightNoIso[leptonCounter] = false;
		
		if( tools::SSSyncIsTightNoSIP(*iM, _minPt0, *theVertices->begin(), _looseD0Mu) ) 
        	_istightNoIsoSIP[leptonCounter] = true;
        else _istightNoIsoSIP[leptonCounter] = false;
		
		if( tools::MVAIsTightNoIsoSIP(*iM, _minPt0, *theVertices->begin(), _looseD0Mu)  ) 
        	_istightMVANoIsoSIP[leptonCounter] = true;
        else _istightMVANoIsoSIP[leptonCounter] = false;
        
		
		_istightMVANoIsoSIP_LMVA[leptonCounter] = false;
		
		
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iM->pt(), iM->eta(), iM->phi(), iM->energy());
		_selectedLeptonPt[leptonCounter] = iM->pt();
		_selectedLeptonEta[leptonCounter] = iM->eta();
		_selectedLeptonPhi[leptonCounter] = iM->phi();

        if(debug) std::cout<<"3.9.4\n";
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);    
		
		
		if(_3dIPsig[leptonCounter] < 4 && _chargeConsistent[leptonCounter] && _missingHits[leptonCounter] < 1 && _miniIsolation[leptonCounter] < 0.4  && _ipPV[leptonCounter] < 0.05 
			&& fabs(_ipZPV[leptonCounter]) < 0.1  && ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt() > 10.0 && _istightMVANoIsoSIP[leptonCounter]){
			_isFO[leptonCounter] = true;
		}
		else
			_isFO[leptonCounter] = false;
	
	
        _closeJetPtAll[leptonCounter] = 0;
        _closeJetAngAll[leptonCounter] = 10000;
        _ptRelAll[leptonCounter] = 0;
		_ptRatioAll[leptonCounter] = 0.0001;
        
        _closeIndex[leptonCounter] = 0;
		if(debug) std::cout<<"3.9.5\n";
		
		//std::cout<<"mu pt = "<<iM->pt()<<"\n";
        for(unsigned int k = 0 ; k < SelectedJetsAll.size() ;k++ ){
			
			double uncPt = (SelectedJetsAll[k]->correctedP4("Uncorrected")).Pt();
        	double uncEta = (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta();
			double L1corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[k]->jetArea(),_corrLevel);
			double L2L3corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[k]->jetArea(),_corrLevelAbs);
			
			
			
			TLorentzVector pJet; pJet.SetPtEtaPhiE( uncPt*L1corr, _jetEtaAll[k], _jetPhiAll[k], L1corr*(SelectedJetsAll[k]->correctedP4("Uncorrected")).E() );
			
			pJet -= *((TLorentzVector *)_leptonP4->At(leptonCounter));
			
			pJet *= L2L3corr/L1corr;
			
			pJet += *((TLorentzVector *)_leptonP4->At(leptonCounter));
			
			
			
			/*
			double PtUse = uncPt*L1corr;
			PtUse -= iM->pt();
			PtUse *= L2L3corr/L1corr;
			PtUse += iM->pt();
			
			TLorentzVector pJet; pJet.SetPtEtaPhiE( PtUse, _jetEtaAll[k], _jetPhiAll[k], L1corr*(SelectedJetsAll[k]->correctedP4("Uncorrected")).E() );
			*/
			//double LepAwarePt = (uncPt*L1corr - iM->pt())*(L2L3corr/L1corr) + iM->pt();
			//double LepAwareE = (L2L3corr/L1corr)*((SelectedJetsAll[k]->correctedP4("Uncorrected")).E()*L1corr - iM->energy()) + iM->energy();
			
			//TLorentzVector pJet; pJet.SetPtEtaPhiE( LepAwarePt, _jetEtaAll[k], _jetPhiAll[k], LepAwareE );
			
			//pJet-=*((TLorentzVector *)_leptonP4->At(leptonCounter));//need still?yes
            
			double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pJet );
			
			//std::cout<<"\nuncorr pt = "<<uncPt<<", corrected = "<<pJet.Pt()<<", unceta = "<<uncEta<<", reg eta = "<<pJet.Eta()<<", uncphi = "<<(SelectedJetsAll[k]->correctedP4("Uncorrected")).Phi()<<", regphi = "<<pJet.Phi()<<", ang = "<<ang<<"\n";
			//std::cout<<"L1 corr = "<<L1corr<<", L2L3 = "<<L2L3corr/L1corr<<", full = "<<L2L3corr<<"\n";
			
            if (ang < _closeJetAngAll[leptonCounter]) {
                _closeJetAngAll[leptonCounter] = ang;
                _closeJetPtAll[leptonCounter] = pJet.Pt();//LepAwarePt;//_jetPtAll[k]
				//_closeJetEAll[leptonCounter] = corr2/corr*((SelectedJetsAll[k]->correctedP4("Uncorrected")).E()*corr - ((TLorentzVector *)leptonP4->At(leptonCounter))->E()) + ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
                _ptRelAll[leptonCounter] =  ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect() - ((TLorentzVector *)_leptonP4->At(leptonCounter))->Vect()); 
										  //((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect());//
                _ptRatioAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()/pJet.Pt();
				_closeIndex[leptonCounter] = k;
				
				
				if(PD) std::cout<<"Closest jet "<<k<<" muon "<<((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()<<",: pT = "<<_jetPtAll[k]<<
						", corr pT = "<<_closeJetPtAll[leptonCounter]<<", eta = "<<_jetEtaAll[k]<<", phi = "<<_jetPhiAll[k]<<", dR = "<<ang<<", L1corr = "<<L1corr<<
						", L2L3 = "<<L2L3corr<<", ptrel = "<<_ptRelAll[leptonCounter]<<", ptratio = "<<_ptRatioAll[leptonCounter]<<"\n";
            }
        }
       
        
        leptonCounter++;

        
    }
  
  
    
	if(debug) std::cout<<"4\n";
	
	int EleFOIndex = -1;
	
    //std::cout<<"before ele\n";
    for(unsigned int i = 0 ; i < sEl.size() ;i++ ){
        
        const pat::Electron *iE = sEl[i];
		//const pat::Electron 
        _leptonIndex = i;
        if (leptonCounter == 8) continue;
        _flavors[leptonCounter] = 0;
		_pdgids[leptonCounter] = iE->pdgId();
        _charges[leptonCounter] = iE->charge();
        _miniIsolation[leptonCounter] = getMiniIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iE), 0.05, 0.2, 10., false, false,myRhoJetsNC);
		//std::cout<<"miniIsoE = "<<_miniIsolation[leptonCounter]<<"\n";
        _isolation[leptonCounter] = pfRelIso15(iE, myRhoJets);
        _ipPV[leptonCounter] = TMath::Abs(iE->gsfTrack()->dxy(PV));
		_ipZPV[leptonCounter] = TMath::Abs(iE->gsfTrack()->dz(PV));
        
		_TrackIso[leptonCounter] = iE->dr03TkSumPt()/iE->pt();
		_EcalPFIso[leptonCounter] = iE->ecalPFClusterIso()/iE->pt();
		_HcalPFIso[leptonCounter] = iE->hcalPFClusterIso()/iE->pt();
		
		
        _sameSign[int(_isolation[leptonCounter] > 0.1)][int((_charges[leptonCounter]+1)/2)]++;
		
		if(debug) std::cout<<"4.1\n";
        
        double Aeff[5] = { 0.1013, 0.0988, 0.0572, 0.0842, 0.1530 };
		double Aeff_Fall15Anal[7] = { 0.1752, 0.1862, 0.1411, 0.1534, 0.1903 , 0.2243, 0.2687 };
        double CorrectedTerm=0.0;
		double miniR = TMath::Max(0.05,TMath::Min(0.2,10.0/iE->pt()));
		double mR2 = miniR*miniR;
        if( TMath::Abs( iE->superCluster()->eta() ) < 1.0 ) CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 0 ]*(mR2/0.09);
        else if( TMath::Abs( iE->superCluster()->eta() ) > 1.0 && TMath::Abs( iE->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 1 ]*(mR2/0.09);
        else if( TMath::Abs( iE->superCluster()->eta() ) > 1.479 && TMath::Abs( iE->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 2 ]*(mR2/0.09);
        else if( TMath::Abs( iE->superCluster()->eta() ) > 2.0 && TMath::Abs( iE->superCluster()->eta() ) < 2.2  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 3 ]*(mR2/0.09);
        else if( TMath::Abs( iE->superCluster()->eta() ) > 2.2 && TMath::Abs( iE->superCluster()->eta() ) < 2.3  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 4 ]*(mR2/0.09);
		else if( TMath::Abs( iE->superCluster()->eta() ) > 2.3 && TMath::Abs( iE->superCluster()->eta() ) < 2.4  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 5 ]*(mR2/0.09);
        else if( TMath::Abs( iE->superCluster()->eta() ) > 2.4 && TMath::Abs( iE->superCluster()->eta() ) < 2.5  )   CorrectedTerm = myRhoJets * Aeff_Fall15Anal[ 6 ]*(mR2/0.09);
        
		if(debug) std::cout<<"4.2\n";
		
        _isolationComponents[leptonCounter][0] = iE->pfIsolationVariables().sumChargedHadronPt/iE->pt();
        _isolationComponents[leptonCounter][1] = iE->pfIsolationVariables().sumNeutralHadronEt/iE->pt();
        _isolationComponents[leptonCounter][2] = iE->pfIsolationVariables().sumPhotonEt/iE->pt();
        _isolationComponents[leptonCounter][3] = CorrectedTerm/iE->pt();// CorrectedTerm/iE->pt() 
		
		_miniIsolationEA[leptonCounter] = TMath::Max(0.0,_miniIsolation[leptonCounter] - CorrectedTerm/iE->pt());
        
		
		if(debug) std::cout<<"4.3\n";
		
		if(iE->genParticle())
			_genCharge[leptonCounter] = iE->genParticle()->pdgId();
		else
			_genCharge[leptonCounter] = 0;
		
		std::vector<const reco::Candidate*> moms;
		if(iE->genParticle() && iE->genParticle()->mother(0)){
		moms.push_back(iE->genParticle()->mother(0));
		for(int i=0;i<30;i++){
			if(!moms[i] || !moms[i]->mother(0))
			break;
		
			moms.push_back(	moms[i]->mother(0) );
		}
		}
		for(unsigned int j=0;j<moms.size();j++){
			_inheritance[leptonCounter][j] = moms[j]->pdgId();
		}
		
		if(debug) std::cout<<"4.4\n";
		
		float value = -10.0;
		float mvaVal = -10.0;//myMVATrig->mvaValue(*iE,false);//
		for(auto ele=thePatElectrons->begin();ele!=thePatElectrons->end();++ele){
		const edm::Ptr<pat::Electron> elePtr(thePatElectrons,ele-thePatElectrons->begin()); //value map is keyed of edm::Ptrs so we need to make one
			if(iE->gsfTrack() == ele->gsfTrack())
				mvaVal = (*valuesMap)[elePtr];
		}
		
		
		double mvaThresh = 1.0, LmvaThresh = 1.0;
		_MVAVal[leptonCounter] = mvaVal;
		if( TMath::Abs(iE->eta()) < 0.8){
			mvaThresh = 0.87;
			LmvaThresh = -0.11;//-0.32
		}
		else if(TMath::Abs(iE->eta()) >= 0.8 && TMath::Abs(iE->eta()) < 1.479){
			mvaThresh = 0.60;
			LmvaThresh = -0.35;//-0.50
		}
		else{
			mvaThresh = 0.17;
			LmvaThresh = -0.55;//-0.70
		}
		
		_missingHits[leptonCounter] = iE->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
	 	
		if(mvaVal > mvaThresh && iE->isGsfCtfScPixChargeConsistent() && _missingHits[leptonCounter] < 1)
			_istightMVANoIsoSIP[leptonCounter] = true;
        else _istightMVANoIsoSIP[leptonCounter] = false;
		
		if(mvaVal > LmvaThresh && iE->isGsfCtfScPixChargeConsistent() && _missingHits[leptonCounter] < 1)
			_istightMVANoIsoSIP_LMVA[leptonCounter] = true;
        else _istightMVANoIsoSIP_LMVA[leptonCounter] = false;
		
		
		
        _3dIP[leptonCounter]    = iE->dB(pat::Electron::PV3D);//ip3D.value();
        _3dIPerr[leptonCounter] = iE->edB(pat::Electron::PV3D);//ip3D.error();
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);//ip3D.significance();
		
		_chargeConsistent[leptonCounter] = iE->isGsfCtfScPixChargeConsistent();//
		
		if(debug) std::cout<<"e mc1\n";
		const reco::GenParticle *mc = iE->genParticle();
		_origin[leptonCounter] = -1;
		//if(mc)
		//	_origin[leptonCounter] = GPM.origin(mc);
			
		if(debug) std::cout<<"e mc2\n";
		_originReduced[leptonCounter] = GPM.originReduced(_origin[leptonCounter]);
        
		
		
		_dphi[leptonCounter] = TMath::Abs(iE->deltaPhiSuperClusterTrackAtVtx());
		_deta[leptonCounter] = TMath::Abs(iE->deltaEtaSuperClusterTrackAtVtx());
		_sigIeta[leptonCounter] = TMath::Abs(iE->scSigmaIEtaIEta());
		_HoE[leptonCounter] = TMath::Abs(iE->hadronicOverEm());
		_pMe[leptonCounter] = TMath::Abs((1.0/iE->ecalEnergy()) - (1.0 / iE->p()));//TMath::Abs(1.0/iE->ecalEnergy() - iE->eSuperClusterOverP()/iE->ecalEnergy());
		
        _closeJetPtAll[leptonCounter] = 0;
        _closeJetAngAll[leptonCounter] = 10000;
        _ptRelAll[leptonCounter] = 0;
		_ptRatioAll[leptonCounter] = 0.001;
        
        
        _isloose[leptonCounter] = (_isolation[leptonCounter]<0.5 && _ipPV[leptonCounter] < 0.05 && _ipZPV[leptonCounter] < 0.01);
		
	
		
		if (_isloose[leptonCounter])
        	_istight[leptonCounter] = ( tools::SSSyncIsTight(*iE, _minPt0, *theVertices->begin(), _looseD0E, _chargeConsistency, theConversions, BS, _3dIPsig[leptonCounter]) ) && (_isolation[leptonCounter]<0.1);
        else _istight[leptonCounter] = false;
		
		if( tools::SSSyncIsTight(*iE, _minPt0, *theVertices->begin(), _looseD0E, _chargeConsistency, theConversions, BS, _3dIPsig[leptonCounter]) ) 
        	_istightNoIso[leptonCounter] = true;
        else 
		    _istightNoIso[leptonCounter] = false;
			
			
		if( tools::SSSyncIsTightNoSIP(*iE, _minPt0, *theVertices->begin(), _looseD0E, _chargeConsistency, theConversions, BS) ) 
        	_istightNoIsoSIP[leptonCounter] = true;
        else 
		    _istightNoIsoSIP[leptonCounter] = false;
		
		
		
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iE->pt(), iE->eta(), iE->phi(), iE->energy());
		_selectedLeptonPt[leptonCounter] = iE->pt();
		_selectedLeptonEta[leptonCounter] = iE->eta();
		_selectedLeptonPhi[leptonCounter] = iE->phi();
 
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
      
	  	if(_3dIPsig[leptonCounter] < 4 && _chargeConsistent[leptonCounter] && _missingHits[leptonCounter] < 1 && _miniIsolation[leptonCounter] < 0.4  && _ipPV[leptonCounter] < 0.05 
			&& fabs(_ipZPV[leptonCounter]) < 0.1  && ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt() > 10.0 && _istightMVANoIsoSIP[leptonCounter]){
			_isFO[leptonCounter] = true;
		}
		else
			_isFO[leptonCounter] = false;
        
        
        _closeIndex[leptonCounter] = 0;
        for(unsigned int k = 0 ; k < SelectedJetsAll.size() ;k++ ){
		
		
			double uncPt = (SelectedJetsAll[k]->correctedP4("Uncorrected")).Pt();
        	double uncEta = (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta();
			double L1corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[k]->jetArea(),_corrLevel);
			double L2L3corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[k]->jetArea(),_corrLevelAbs);
			
			TLorentzVector pJet; pJet.SetPtEtaPhiE( uncPt*L1corr, _jetEtaAll[k], _jetPhiAll[k], L1corr*(SelectedJetsAll[k]->correctedP4("Uncorrected")).E() );
			
			pJet -= *((TLorentzVector *)_leptonP4->At(leptonCounter));
			
			pJet *= L2L3corr/L1corr;
			
			pJet += *((TLorentzVector *)_leptonP4->At(leptonCounter));
			/*
			double PtUse = uncPt*L1corr;
			PtUse -= iE->pt();;
			PtUse *= L2L3corr/L1corr;
			PtUse += iE->pt();
			
			TLorentzVector pJet; pJet.SetPtEtaPhiE( PtUse, _jetEtaAll[k], _jetPhiAll[k], L1corr*(SelectedJetsAll[k]->correctedP4("Uncorrected")).E() );
			*/
			//double LepAwarePt = (uncPt*L1corr - iE->pt())*(L2L3corr/L1corr) + iE->pt();
			//double LepAwareE = (L2L3corr/L1corr)*((SelectedJetsAll[k]->correctedP4("Uncorrected")).E()*L1corr - iE->energy()) + iE->energy();
			
            double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pJet );
            if (ang < _closeJetAngAll[leptonCounter]) {
                _closeJetAngAll[leptonCounter] = ang;
                _closeJetPtAll[leptonCounter] = pJet.Pt();//LepAwarePt;
				//_closeJetEAll[leptonCounter] = corr2/corr*((SelectedJetsAll[k]->correctedP4("Uncorrected")).E()*corr - ((TLorentzVector *)leptonP4->At(leptonCounter))->E()) + ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
                _ptRelAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect() - ((TLorentzVector *)_leptonP4->At(leptonCounter))->Vect());
										 //((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect());
				
				_ptRatioAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()/pJet.Pt();
                _closeIndex[leptonCounter] = k;
				
				//std::cout<<"Closest jet ele: pT = "<<_jetPtAll[k]<<", corr pT = "<<_closeJetPtAll[leptonCounter]<<", eta = "<<_jetEtaAll[k]<<", phi = "<<_jetPhiAll[k]<<", dR = "<<ang<<", L1corr = "<<L1corr<<", L2L3 = "<<L2L3corr<<"\n";
            }
        }
		
      
        leptonCounter++;
        
    }
	

	
	if(debug) std::cout<<"5\n";
	
	
	//if( doFR && leptonCounter < 1) return;	
    if (leptonCounter < 2 && !doFR) return;
	
	_n_Jets = 0;
    _n_bJets = 0;
    HT = 0;
	//std::cout<<"jet cleaning\n";
	double PMC[2] = {1,1}, PData[2] = {1,1};
	double ptbins[19] = {20.,25.,30.,35.,40.,45.,50.,60.,70.,80.,90.,100.,120.,150.,200.,300.,400.,600.,1000000.};
	
	double bTageffs[7][17] = {{0.498440,0.564021,0.610909,0.632182,0.628327,0.674033,0.662796,0.664053,0.675163,0.681013,0.676228,0.678714,0.687691,0.695168,0.667931,0.605706,0.530881},
							  {0.518034,0.582882,0.628564,0.652494,0.651326,0.667326,0.683174,0.687201,0.697793,0.701271,0.698518,0.700552,0.707651,0.709330,0.682399,0.625348,0.553085},
							  {0.488457,0.555124,0.601580,0.625985,0.625929,0.641381,0.657951,0.667800,0.676599,0.681218,0.681516,0.682655,0.682393,0.682125,0.646767,0.596776,0.520092},
							  {0.452275,0.551838,0.566909,0.597384,0.587755,0.599446,0.619942,0.623788,0.630893,0.639091,0.634120,0.629618,0.633678,0.605852,0.564228,0.511104,0.422896},
							  {0.425729,0.495264,0.542952,0.576627,0.576082,0.576421,0.596198,0.602617,0.606897,0.614185,0.608507,0.604401,0.603599,0.568326,0.525934,0.484535,0.409251},
							  {0.330925,0.399433,0.448948,0.480190,0.477676,0.485342,0.504022,0.501825,0.507183,0.511612,0.501227,0.495817,0.504509,0.468249,0.422269,0.382096,0.345576},
							  {0.061321,0.081894,0.095642,0.105105,0.105836,0.110704,0.115588,0.115834,0.119451,0.123472,0.117514,0.122229,0.129296,0.120606,0.124337,0.121655,0.200000}};
	
	
    for(unsigned int i = 0 ; i < SelectedJets.size() ;i++ ){
        _jetEta[_n_Jets] = SelectedJets[i]->eta();
        _jetPhi[_n_Jets] = SelectedJets[i]->phi();
		
        _jetPt[_n_Jets] = SelectedJets[i]->pt();
		
		
		if(Dov1v3Corrections){
		
		double uncPt = (SelectedJets[i]->correctedP4("Uncorrected")).Pt();
        double uncEta = (SelectedJets[i]->correctedP4("Uncorrected")).Eta();
		double L1corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[i]->jetArea(),_corrLevel);
		double L2L3corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJetsAll[i]->jetArea(),_corrLevelFull);
		_jetPt[_n_Jets] = uncPt*L2L3corr;
       
	   
	   }
	   
	   int ebin = 0, pbin = 0;
	   
	   for(int b=1;b<8;b++){
	   
	   		if(fabs(_jetEta[_n_Jets]) < b*0.4 && fabs(_jetEta[_n_Jets]) >= (b-1)*0.4)
				ebin = b-1;
	   
	   }
	   
	   for(int b=1;b<18;b++){
	   	
	   		if(_jetPt[_n_Jets] < ptbins[b] && _jetPt[_n_Jets] >= ptbins[b-1])
				 pbin = b-1;
	   
	   }
	   
	  
	   double JetSF = reader->eval(BTagEntry::FLAV_B, _jetEta[_n_Jets], _jetPt[_n_Jets]);
	
        //TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[_n_Jets],_jetEta[_n_Jets],_jetPhi[_n_Jets],0);
       // double dR1 = ((TLorentzVector *)_leptonP4->At(_index1))->DeltaR( jt );
       // double dR2 = ((TLorentzVector *)_leptonP4->At(_index2))->DeltaR( jt );
        
        //if (dR1 < 0.4 || dR2 < 0.4) continue;
        
        _csv[_n_Jets] = SelectedJets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		
		if(_csv[_n_Jets] > 0.89){
			PMC[0] *= bTageffs[ebin][pbin];
			PData[0] *= JetSF*bTageffs[ebin][pbin];
		}
		else{
			PMC[1] *= (1-bTageffs[ebin][pbin]);
			PData[1] *= (1-(JetSF*bTageffs[ebin][pbin]));
		}
		
		//if(SelectedJets[i]->pt() > 40)
		//	std::cout<<"csv = "<<SelectedJets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")<<"\n";
        
        if(_csv[_n_Jets] > 0.814) {
            _bTagged[_n_Jets] = true;
            _n_bJets++;
        } else _bTagged[_n_Jets] = false;
		
		//std::cout<<"jet pT = "<<_jetPt[_n_Jets]<<", eta = "<<_jetEta[_n_Jets]<<", phi = "<<_jetPhi[_n_Jets]<<"\n";
        
        HT+= _jetPt[_n_Jets];
        _n_Jets++;
    }
	
	
	_weight *= (PData[0]*PData[1])/(PMC[0]*PMC[1]);
	
	if(debug) std::cout<<"8";
	
    _nLeptons = leptonCounter;
   
   // if(PD) std::cout<<"before 1 Lep\n";
   
	
	if( doFR && leptonCounter < 1) return;	
	
	bool FillFR = false;
	
	//if(PD) std::cout<<" 1 Lep\n";
	
	if( doFR){
	std::vector<int> FOIndex;
	for(int j=0;j<_nLeptons;j++){
	
		((TLorentzVector*)_leptonP4->At(j))->SetPtEtaPhiM(_selectedLeptonPt[j],_selectedLeptonEta[j],_selectedLeptonPhi[j],0);
	
		if(_3dIPsig[j] < 4 && _chargeConsistent[j] && _missingHits[j] < 1 && _miniIsolation[j] < 0.4  && _ipPV[j] < 0.05 
				&& fabs(_ipZPV[j]) < 0.1  && ((TLorentzVector *)_leptonP4->At(j))->Pt() > 10.0 
				){//FO Lep Req
				
				
			
				//if(PD) std::cout<<"lepton"<<j<<" pt/eta/phi = "<<_selectedLeptonPt[j]<<"/"<<_selectedLeptonEta[j]<<"/"<<_selectedLeptonPhi[j]<<"\n";
				
				if(_flavors[j] && _istightMVANoIsoSIP[j])
					FOIndex.push_back(j);
				else if(!_flavors[j])
					FOIndex.push_back(j);
			
		}
	}
	
	
	
	bool awayJet = false;
	
	for(unsigned int fos=0;fos<FOIndex.size();fos++){
	
		for(int i = 0 ; i < _n_Jets ;i++ ){
		
			TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],0);
			
			
			if( _jetPt[i] > 40 && ((TLorentzVector*)_leptonP4->At(FOIndex[fos]))->DeltaR(jt) > 1){
				awayJet = true;
				
				//if(PD) std::cout<<"Jet"<<i<<" pt/eta/phi = "<<_jetPt[i]<<"/"<<_jetEta[i]<<"/"<<_jetPhi[i]<<"\n";
			}
		
		}
	
		if(awayJet )//&& _met < 20 && _mt[FOIndex[0]] < 20)
			FillFR = true;
	
	
	}
	
	/*if(FOIndex.size() == 1){
	
	if(PD) std::cout<<"1 FO\n";
	
		for(int i = 0 ; i < _n_Jets ;i++ ){
		
			TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],0);
			
			
			if( _jetPt[i] > 40 && ((TLorentzVector*)_leptonP4->At(FOIndex[0]))->DeltaR(jt) > 1){
				awayJet = true;
				
				if(PD) std::cout<<"Jet"<<i<<" pt/eta/phi = "<<_jetPt[i]<<"/"<<_jetEta[i]<<"/"<<_jetPhi[i]<<"\n";
			}
		
		}
	
		if(awayJet )//&& _met < 20 && _mt[FOIndex[0]] < 20)
			FillFR = true;
	
	}*/
	if(PD){
	if(awayJet)
		std::cout<<"away Jet \n";
	
	if(_met < 20)
		std::cout<<"MET < 20 and = "<<_met<<", phi = "<<_met_phi<<"\n";
	else
		std::cout<<"MET > 20 and = "<<_met<<", phi = "<<_met_phi<<"\n";
	
	if(_mt[FOIndex[0]] < 20)
		std::cout<<"MT < 20 and = "<<_mt[FOIndex[0]]<<"\n";
	else
		std::cout<<"MT > 20 and = "<<_mt[FOIndex[0]]<<"\n";
	
	}
	}
	
	if(doFR && !FillFR) return;
	
	
    outputTree->Fill();

}




DEFINE_FWK_MODULE(SSSync);
