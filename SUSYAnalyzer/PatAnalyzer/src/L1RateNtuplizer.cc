#include "L1RateNtuplizer.h"
#include "DataFormats/Common/interface/Handle.h"
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

L1RateNtuplizer::L1RateNtuplizer(const edm::ParameterSet & iConfig) :
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
_jetPtCut(25.),
_jetEtaCut(2.4),
_tauPt(20),
_tauEta(2.3),
_regression(false)
{
    Sample              = iConfig.getUntrackedParameter<std::string>("SampleLabel") ;
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
	std::vector<std::string> myManualCatWeigths;
	
	
	myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml");
	myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml");
	myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml");
	myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml");
	myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml");
	myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml");
	
	
	//myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_50ns_EB_BDT.weights.xml");
   // myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_50ns_EE_BDT.weights.xml");
	
	vector<string> myManualCatWeigthsTrig;
	string the_path;
    for (unsigned i  = 0 ; i < myManualCatWeigths.size() ; i++){
        the_path = edm::FileInPath ( myManualCatWeigths[i] ).fullPath();
        myManualCatWeigthsTrig.push_back(the_path);
    }
	
	myMVATrig = new EGammaMvaEleEstimatorCSA14();
    myMVATrig->initialize("BDT",
                          EGammaMvaEleEstimatorCSA14::kNonTrigPhys14,
						  //EGammaMvaEleEstimatorCSA14::kTrig,
                          true,
                          myManualCatWeigthsTrig);
    
    //outfile = fopen("FakeSync.txt", "w");
}




void L1RateNtuplizer::beginJob()
{
	
	
	
	
	Triggers = fs->make<TH1F>("Triggers","Triggers",450,0,450);

    Nvtx           = fs->make<TH1F>("N_{vtx}"        , "Number of vertices;N_{vtx};events / 1"  ,    40, 0., 40.);
    
    _hCounter = fs->make<TH1D>("hCounter", "Events counter", 5,0,5);
	MvaHist = fs->make<TH1D>("MvaHist","MvaHist",50,-1,1);
	
	std::cout<<"treeb\n";
    outputTree = new TTree("L1RateNtuplizerTree","L1RateNtuplizerTree");
    std::cout<<"treea\n";
   
	outputTree->Branch("_nL1Mus", &_nL1Mus, "_nL1Mus/I");
	outputTree->Branch("_L1HT", &_L1HT, "_L1HT/D");
	outputTree->Branch("_L1Mu_Pt", &_L1Mu_Pt, "_L1Mu_Pt[40]/D");
	outputTree->Branch("_L1Mu_Eta", &_L1Mu_Eta, "_L1Mu_Eta[40]/D");
	outputTree->Branch("_L1Mu_Phi", &_L1Mu_Phi, "_L1Mu_Phi[40]/D");

	

    
    outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
    outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
    outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
	
	
	
	outputTree->Branch("_weight", &_weight , "_weight/D");
	

    GPM = GenParticleManager();
	
	
    _nEventsTotal = 0;
    _nEventsFiltered = 0;
    _nEventsTotalCounted = 0;
	
	
	
    
}

void L1RateNtuplizer::endJob() {
    //outputTree -> Write();
    // store nEventsTotal and nEventsFiltered in preferred way
    std::cout<<_nEventsTotalCounted<<std::endl;
    std::cout<<_nEventsFiltered<<std::endl;
    
   // delete fMetCorrector;
    
    
}

void L1RateNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{

	bool debug = true;
	edm::Handle<edm::TriggerResults> trigResults;
	edm::InputTag trigResultsTag("TriggerResults","","HLT");
	iEvent.getByLabel(trigResultsTag,trigResults);
	
	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	edm::InputTag triggerObjects_("selectedPatTrigger","","PAT");//PAT for MC and RECO for data
	iEvent.getByLabel(triggerObjects_, triggerObjects);
	
	//============ Primary vertices ============
    edm::InputTag IT_goodVtx = edm::InputTag("goodOfflinePrimaryVertices");
    edm::Handle<std::vector<Vertex> > theVertices;
    iEvent.getByLabel( "goodOfflinePrimaryVertices", theVertices) ;
    if( ! theVertices.isValid() ) ERR(IT_goodVtx ) ;
    int nvertex = theVertices->size();
    /*
    _n_PV = nvertex;
    
    Nvtx->Fill(TMath::Min(nvertex,39));//
    if(! nvertex ){
        cout << "[WARNING]: No candidate primary vertices passed the quality cuts, so skipping event" << endl;
        return ;
    }
	*/
	
	/*
	
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
    }*/

    _runNb = iEvent.id().run();
    _eventNb = iEvent.id().event();
    _lumiBlock = iEvent.luminosityBlock();
	
		
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
	
	
    //============ Counter done ============
    
    //============ Beamspot ============
    edm::Handle< reco::BeamSpot > theBeamSpot;
    iEvent.getByLabel( IT_beamspot, theBeamSpot );
    if( ! theBeamSpot.isValid() ) ERR( IT_beamspot ) ;
    BeamSpot::Point  BS= theBeamSpot->position();
    //==================================
    
	
	//=========== L1 Particles =============
	edm::Handle<std::vector<l1extra::L1EtMissParticle>> L1HTContainer;
	iEvent.getByLabel("l1extraParticles","MHT",L1HTContainer);
	
	_L1HT = L1HTContainer->begin()->etTotal();
	
	edm::Handle<std::vector<l1extra::L1MuonParticle>> L1MuContainer;
	iEvent.getByLabel("l1extraParticles","",L1MuContainer);
	
	
	for(int i=0;i<40;i++){
		
		
		_L1Mu_Pt[i] = -999;
		_L1Mu_Eta[i] = -999;
		_L1Mu_Phi[i] = -999;
	
	
	}
	
	_nL1Mus = 0;
	if(_nL1Mus > 40)
		_nL1Mus = 40;
		
	for(std::vector<l1extra::L1MuonParticle>::const_iterator nm = L1MuContainer->begin();nm != L1MuContainer->end();nm++){
	
		
		if(_nL1Mus > 40) continue;
		
		
		_L1Mu_Pt[_nL1Mus] = nm->pt();
		_L1Mu_Eta[_nL1Mus] = nm->eta();
		_L1Mu_Phi[_nL1Mus] = nm->phi();
		
	
	
		_nL1Mus++;
	}
	//======================================
    
    
 
	_weight = 1;
	if(Sample=="ElectronsMC"){//============Individual Event Weight=============
	edm::Handle<GenEventInfoProduct> pdfvariables;
	iEvent.getByLabel("generator",pdfvariables);
	_weight = pdfvariables->weight();
	
	if(_weight > 0)
		_weight = 1;
	else
		_weight = -1;
    
	}
	


    outputTree->Fill();

    
}




DEFINE_FWK_MODULE(L1RateNtuplizer);
