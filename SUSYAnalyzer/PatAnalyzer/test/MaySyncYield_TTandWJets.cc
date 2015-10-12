#include <iostream>
#include <sstream>
#include <TStyle.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TGaxis.h>
#include <THStack.h>
#include <TLegend.h>
#include <TObject.h>
#include <TF1.h>
#include <TF2.h>


int SR(int nbjets, float met, int njets, float HT, int analysis){///analysis = 0 for low pt and 1 for high pt

	///sr[nbjets][met][njets][HT]
	int sr[4][3][3][3] = {{{{-1,-1,-1},{0,0,0},{0,0,0}},{{-1,-1,-1},{0,1,2},{0,3,4}},{{-1,-1,-1},{0,5,6},{0,7,8}}},
						  {{{-1,-1,-1},{10,10,10},{10,10,10}},{{-1,-1,-1},{10,11,12},{10,13,14}},{{-1,-1,-1},{10,15,16},{10,17,18}}},
					 	  {{{-1,-1,-1},{20,20,20},{20,20,20}},{{-1,-1,-1},{20,21,22},{20,23,24}},{{-1,-1,-1},{20,25,26},{20,27,28}}},
						  {{{-1,-1,-1},{30,30,30},{30,30,30}},{{-1,-1,-1},{30,31,32},{30,33,34}},{{-1,-1,-1},{30,35,36},{30,37,38}}}};


	int nj = 0, nbj = 0, m = 0, h = 0, aht[2] = {250,200};
	if(njets == 2 || njets == 3)
		nj = 1;
	else if(njets >= 4)
		nj = 2;
		
	if(met > 50 && met <= 120)
		m = 1;
	else if(met > 120)
		m = 2;
		
	if(HT > aht[analysis] && HT <= 400)
		h = 1;
	else if(HT > 400)
		h = 2;
		
	if(nbjets == 1)
		nbj = 1;
	else if(nbjets == 2)
		nbj = 2;
	else if(nbjets >= 3)
		nbj = 3;
		
	return sr[nbj][m][nj][h];
}

int GetNewSignalRegion(int nbjets, float met, int njets, float HT, float mTmin, int ept){

	if(ept == 0){
	
	int sr = 1;
	
	if(HT < 300)
		return 0;
	
	if(mTmin > 120)
		return 8;
		
	if(nbjets >= 3)
		return 7;
		
	sr += nbjets*2;
	if(met > 200)
		sr++;
		


	return sr;
	
	}
	
	if(ept == 1){
	
	int First6 [2][2][2] = {{{1,2},{3,4}},{{3,5},{3,6}}};//[MET][jets][HT]
	
	int Last4[2][2] = {{1,2},{3,4}};//[MET][HT]
	
	
	int sr = 0, nj = 0, m = 0, ht = 0;
	
	if(HT > 1600 && met <= 500)
		return 26;
		
	if(met > 500)
		return 25;
		
	if(mTmin > 120 && met <= 500)
		return (HT > 300 ? 24 : 23);
		
	if(nbjets > 3)
		nbjets = 3;
	
	sr += nbjets*6;
	
	if(njets > 4) nj = 1;
	
	if(met > 500)
		m = 2;
	else if(met > 200)
		m = 1;
		
	if(HT > 1600)
		ht = 2;
	else if(HT > 300)
		ht = 1;
	
	
	if(mTmin > 120)
		return sr+= Last4[m][ht];
	else
		return sr+= First6[m][nj][ht];
	
	}
	
	
	if(ept == 2){
	
	int First8 [2][2][2][2] = {{{{1,2},{3,4}},{{3,5},{3,6}}},{{{3,7},{3,8}},{{3,8},{3,8}}}};//[mTmin][MET][jets][HT]
	
	int Last6 [2][2][2] = {{{1,2},{3,4}},{{5,6},{5,6}}};//[mTmin][MET][HT]
	
	if(met < 50)
		return 0;
	
	int sr = 0, nj = 0, mt = 0, m = 0, ht = 0;
	
	if(HT > 1600 && met <= 500)
		return 32;
		
	if(met > 500 && HT > 300)
		return 31;
		
	if(nbjets > 3)
		nbjets = 3;
	
	sr += nbjets*8;
	
	if(njets > 4) nj = 1;
	
	if(met > 500)
		m = 2;
	else if(met > 200)
		m = 1;
		
	if(HT > 1600)
		ht = 2;
	else if(HT > 300)
		ht = 1;
		
	if(mTmin > 120)
		mt = 1;
	
	if(nbjets > 2)
		return sr+= Last6[mt][m][ht];
	else
		return sr+= First8[mt][m][nj][ht];
	
	
	}
	
	return 0;

}


int GetHHSignalRegion(int nbjets, float met, int njets, float HT, float mTmin){
	
	int First8 [2][2][2][2] = {{{{1,2},{3,4}},{{3,5},{3,6}}},{{{3,7},{3,8}},{{3,8},{3,8}}}};//[mTmin][MET][jets][HT]
	
	int Last6 [2][2][2] = {{{1,2},{3,4}},{{5,6},{5,6}}};//[mTmin][MET][HT]
	
	if(met < 50)
		return 0;
	
	int sr = 0, nj = 0, mt = 0, m = 0, ht = 0;
	
	if(HT > 1600 && met <= 500)
		return 32;
		
	if(met > 500 && HT > 300)
		return 31;
		
	if(nbjets > 3)
		nbjets = 3;
	
	sr += nbjets*8;
	
	if(njets > 4) nj = 1;
	
	if(met > 500)
		m = 2;
	else if(met > 200)
		m = 1;
		
	if(HT > 1600)
		ht = 2;
	else if(HT > 300)
		ht = 1;
		
	if(mTmin > 120)
		mt = 1;
	
	if(nbjets > 2)
		return sr+= Last6[mt][m][ht];
	else
		return sr+= First8[mt][m][nj][ht];
		

}

int GetHLSignalRegion(int nbjets, float met, int njets, float HT, float mTmin){
	
	int First6 [2][2][2] = {{{1,2},{3,4}},{{3,5},{3,6}}};//[MET][jets][HT]
	
	int Last4[2][2] = {{1,2},{3,4}};//[MET][HT]
	
	
	int sr = 0, nj = 0, m = 0, ht = 0;
	
	if(HT > 1600 && met <= 500)
		return 26;
		
	if(met > 500)
		return 25;
		
	if(mTmin > 120 && met <= 500)
		return (HT > 300 ? 24 : 23);
		
	if(nbjets > 3)
		nbjets = 3;
	
	sr += nbjets*6;
	
	if(njets > 4) nj = 1;
	
	if(met > 500)
		m = 2;
	else if(met > 200)
		m = 1;
		
	if(HT > 1600)
		ht = 2;
	else if(HT > 300)
		ht = 1;
	
	
	if(mTmin > 120)
		return sr+= Last4[m][ht];
	else
		return sr+= First6[m][nj][ht];
		

}

int GetLLSignalRegion(int nbjets, float met, float mTmin, float HT){	
	
	int sr = 1;
	
	if(HT < 300)
		return 0;
	
	if(mTmin > 120)
		return 8;
		
	if(nbjets >= 3)
		return 7;
		
	sr += nbjets*2;
	if(met > 200)
		sr++;
		


	return sr;
		

}


enum CType{BJets,CJets,Light,Tau};

double deltaPhi(TLorentzVector* lepton, double met, double metPhi) {
    double W_px = lepton->Px() + met*TMath::Cos(metPhi);
    double W_py = lepton->Py() + met*TMath::Sin(metPhi);
    double W_phi = TMath::ATan2(W_py,W_px);
    return(fabs(W_phi -lepton->Phi()));
}

void MaySyncYield_TTandWJets()
{
    //int leptonStudy = 1; //1=muon, 0-electron

    
	gStyle->SetOptFit(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetPadColor(kWhite);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetNdivisions(505,"XY");
	
	gStyle->SetLabelFont(132,"XYZ");
	gStyle->SetLabelSize(0.06,"XYZ");
	gStyle->SetLabelOffset(0.001,"X");
	
    gStyle->SetTitleFont(132,"");
	gStyle->SetTitleFont(132,"XYZ");
	gStyle->SetTitleFontSize(0.09);
	gStyle->SetTitleSize(0.07,"XY");
	gStyle->SetTitleXOffset(0.8);
	gStyle->SetTitleYOffset(1.0);
    
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    //gStyle->SetNdivisions(510, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
	
    gStyle->SetPadBottomMargin(0.18);
	gStyle->SetPadRightMargin(0.04);
	gStyle->SetPadTopMargin(0.04);
	gStyle->SetPadLeftMargin(0.15);
    
    gStyle->SetStatFont(132);
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatColor(10);
    gStyle->SetStatFontSize(0.08);
    //gStyle->SetTitleBorderSize(1);
    gStyle->SetTitleFont(132);
    gStyle->SetTitleFontSize(0.08);
	
	gStyle->SetMarkerSize(1.0);
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerColor(gStyle->GetHistLineColor());
	
	gStyle->SetPalette(1,0);
    TGaxis::SetMaxDigits(3);
	gStyle->SetPadGridX(true);
	gStyle->SetPadGridY(true);
	
	gStyle->SetFuncColor(kRed);
    
    int _n_bJets;
    int _n_Jets;
    
    double _jetEta[20];
    double _jetPhi[20];
    double _jetPt[20];
    bool _bTagged[20];
    double _csv[20];
    
    //int _n_bJetsAll;
    int _n_JetsAll;
    
    double _jetEtaAll[100];
    double _jetPhiAll[100];
    double _jetPtAll[100];
    //bool _bTaggedAll[100];
    //double _csvAll[100];
    //int _leptonIndex;
   	int _closeIndex[8];
	
	int _nLeptons;
    
    //int _eventType; //ee,mm,em
    bool _sb;
    bool _doubleF;
    int _index1 = -1;
    int _index2 = -1;

    
    //int _indeces[4];
    int _flavors[8];
    double _charges[8];
	int _pdgids[8];
	double _MVAVal[8];
    double _isolation[8];
	double _miniIsolation[8];
	double _miniIsolationEA[8];
    double _isolationComponents[8][4];
    double _isolationMC[8][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    
    int _origin[8];
    int _originReduced[8];
    double _PVchi2;
    double _PVerr[3];
    double _ipPV[8];
    double _ipPVerr[8];
    double _ipPVmc[8];
    
    double _ipZPV[8];
    double _ipZPVerr[8];
    
    double _3dIP[8];
    double _3dIPerr[8];
    double _3dIPsig[8];
	
	int _missingHits[8];
	double _dphi[8];
	double _deta[8]; 
	double _sigIeta[8];
	double _HoE[8];
	double _pMe[8];
    
    double _mt[8];
    
    //double _closeJetPt[4];
    double _closeJetPtAll[8];
    double _closeJetAngAll[8];
    double _ptRelAll[8];
	
    double _closeJetPtAllMC[4];
    double _closeJetPtAllstatus[4];
    int _partonIdMatched[4];
    bool _sameParton[4];
    
    bool _isloose[8];
    bool _istight[8];
	bool _istightNoIso[8];
	bool _istightNoIsoSIP[8];
	bool _islooseMVA[8];
	bool _istightMVANoIsoSIP[8];
	bool _istightMVANoIsoSIP_LMVA[8];
	bool _istightMVA[8];
    
    int _n_PV;
    
    //int _n_electrons;
    //int _n_muons;
   	// int _n_taus;
    
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;
	
	double _mompt[4];
    double _momphi[4];
    double _mometa[4];
    int _mompdg[4];
    
    
    double _met;
    double _met_phi;
    double HT;
    
    TClonesArray* _leptonP4 = new TClonesArray("TLorentzVector", 8);
    for (int i=0; i!=8; ++i) {
        new ( (*_leptonP4)[i] ) TLorentzVector();
    }
    
    TClonesArray* _jetP4 = new TClonesArray("TLorentzVector", 50);
    for (int i=0; i!=50; ++i) {
        new ( (*_jetP4)[i] ) TLorentzVector();
    }
    
    TClonesArray* _jetAllP4 = new TClonesArray("TLorentzVector", 50);
    for (int i=0; i!=50; ++i) {
        new ( (*_jetAllP4)[i] ) TLorentzVector();
    }
    
	////////////////////////////////
	//// Read in FR from Histos ////
	////////////////////////////////
	
	float fptbins[6] = {10.0,15.0,25.0,35.0,50.0,70.0};
	float etabins[4] = {0.0,1.0,2.0,2.4};
	const char *fv[2] = {"Ele","Mu"};
	const char *whichFR[5] = {"","CorrectedPt","Jet","JetLow","JetHigh"};
	const char *whichFO[5] = {"Tight","FO1","FO2","FO3","FO4"};
	const char *whichSample[2] = {"QCD","InSitu"};
	const char *EventPt[3] = {"LL","HL","HH"};
	
	TH2F *FakeRates[2][2][5][4];//[QCD or InSitu][flavor][pt,corrected pt, jet pt,jet pt low lep pt, jet pt high lep pt][fakeable object]
	TFile *FRfiles[2][2];//[QCD or InSitu][Flavor]
	
	const char* fileNames[2][2] = {{"FakeRateHistosMay28/QCD_Ele_FR.root","FakeRateHistosMay28/QCD_Mu_FR.root"},
								   {"FakeRateHistosMay28/InSitu_Ele_FR_FO4Loose.root","FakeRateHistosMay28/InSitu_Mu_FR.root"}};
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			FRfiles[i][j] = new TFile(fileNames[i][j],"READ");
			for(int k=0;k<5;k++){
				for(int l=0;l<4;l++){
	
					std::cout<<Form("FR2DMeasure%s_%s_%s",whichFR[k],fv[j],whichFO[l+1])<<"\n";
					FakeRates[i][j][k][l] = new TH2F(Form("FR2DMeasure%s_%s_%s",whichFR[k],fv[j],whichFO[l+1]),Form("FR2DMeasure%s_%s_%s",fv[i],whichFR[k],whichFO[l+1]),5,fptbins,3,etabins);
					FakeRates[i][j][k][l]->Read(Form("FR2DMeasure%s_%s_%s",whichFR[k],fv[j],whichFO[l+1]));
					
				}
			}
		}
	}
  
    
   	// TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    double lumi = 10000;
	double ZZto4L = 15.4;
	
	double WJetsToLNu = 20508.9;
	double TTJets = 809.1;//2 files
	
	
	//double SampleXSection[9] = {TTWJets,TTZJets,WZJets,ZZto4L,DYJets,TTJets,WJetsToLNu,T14t12,T14t15};
	//TString Sample = "WJets";
	TString Sample = "TTJets";
	
	TH1F *SRYields[3];
	
	SRYields[0] = new TH1F("LL_Yields","Low-Low Signal Region Yields",8,1,9);
	SRYields[1] = new TH1F("HL_Yields","High-Low Signal Region Yields",25,1,26);
	SRYields[2] = new TH1F("HH_Yields","High-High Signal Region Yields",32,1,33);
    
    
	//int NumEvents[9] = {0};
	//Double_t SampleScales[9] = {0};
	//int WhichSample = 5;
	FILE *EvtYield = fopen("NewestSRYields.txt","w");
	
	int evs = 0;
    	std::cout<<"\n\n";

	
	TH1F *TotalEvents = new TH1F("TotalEvents","Total Events",5,0,5);
	
	//TFile *hfiles = new TFile("file:MaySyncOutput.root","READ");
	TFile *hfiles = new TFile("/afs/cern.ch/work/m/mcarver/CMSSW_7_2_1_patch1/src/TTJetsNtuples/TTJets.root","READ");//
	hfiles->cd("SyncExercise");
	TH1D* _hCounters = new TH1D("hCounter", "Events counter", 5,0,5);
	_hCounters->Read("hCounter");
	evs += _hCounters->GetBinContent(1);
	
	Double_t scale = lumi*TTJets/evs;
	
	std::cout<<"scale = "<<scale<<"\n";
	
	
	TotalEvents->SetBinContent(1,evs);
	
    TChain *outputTree=new TChain("SyncExercise/SSSyncTree");
	outputTree->Add("/afs/cern.ch/work/m/mcarver/CMSSW_7_2_1_patch1/src/TTJetsNtuples/TTJets.root");
	
	outputTree->SetBranchAddress("_leptonP4", &_leptonP4);
	
	outputTree->SetBranchAddress("_eventNb",   &_eventNb); 
    outputTree->SetBranchAddress("_runNb",     &_runNb);    
    outputTree->SetBranchAddress("_lumiBlock", &_lumiBlock);
    
    outputTree->SetBranchAddress("_nLeptons", &_nLeptons);
    outputTree->SetBranchAddress("_flavors", &_flavors);
    outputTree->SetBranchAddress("_charges", &_charges);
	outputTree->SetBranchAddress("_pdgids", &_pdgids);
   	outputTree->SetBranchAddress("_isolation", &_isolation);
	outputTree->SetBranchAddress("_MVAVal", &_MVAVal);
	outputTree->SetBranchAddress("_miniIsolation", &_miniIsolation);
	outputTree->SetBranchAddress("_miniIsolationEA", &_miniIsolationEA);
   	outputTree->SetBranchAddress("_isolationComponents", &_isolationComponents);
    outputTree->SetBranchAddress("_isolationMC", &_isolationMC);
    
    outputTree->SetBranchAddress("_index1", &_index1);
    outputTree->SetBranchAddress("_index2", &_index2);

    outputTree->SetBranchAddress("_sb", &_sb);
    outputTree->SetBranchAddress("_doubleF", &_doubleF);
    
    outputTree->SetBranchAddress("_origin", &_origin);
    outputTree->SetBranchAddress("_originReduced", &_originReduced); 
    
    outputTree->SetBranchAddress("_PVchi2", &_PVchi2);
    outputTree->SetBranchAddress("_PVerr", &_PVerr);
    
    outputTree->SetBranchAddress("_ipPV", &_ipPV);
    outputTree->SetBranchAddress("_ipPVerr", &_ipPVerr);
    outputTree->SetBranchAddress("_ipZPV", &_ipZPV);
    outputTree->SetBranchAddress("_ipZPVerr", &_ipZPVerr);
    
    outputTree->SetBranchAddress("_ipPVmc", &_ipPVmc);
    
    outputTree->SetBranchAddress("_3dIP", &_3dIP);
    outputTree->SetBranchAddress("_3dIPerr", &_3dIPerr);
    outputTree->SetBranchAddress("_3dIPsig", &_3dIPsig);
	
	outputTree->SetBranchAddress("_missingHits", &_missingHits);
	outputTree->SetBranchAddress("_dphi", &_dphi);
	outputTree->SetBranchAddress("_deta", &_deta);
	outputTree->SetBranchAddress("_sigIeta", &_sigIeta);
	outputTree->SetBranchAddress("_HoE", &_HoE);
	outputTree->SetBranchAddress("_pMe", &_pMe);
    
    
    outputTree->SetBranchAddress("_mt", &_mt);
    outputTree->SetBranchAddress("_isloose", &_isloose);
	outputTree->SetBranchAddress("_islooseMVA", &_islooseMVA);
    outputTree->SetBranchAddress("_istight", &_istight);
	outputTree->SetBranchAddress("_istightMVA", &_istightMVA);
	outputTree->SetBranchAddress("_istightNoIso", &_istightNoIso);
	outputTree->SetBranchAddress("_istightNoIsoSIP", &_istightNoIsoSIP);
	outputTree->SetBranchAddress("_istightMVANoIsoSIP", &_istightMVANoIsoSIP);
    
    
    outputTree->SetBranchAddress("_closeJetPtAll", &_closeJetPtAll); 
    outputTree->SetBranchAddress("_closeJetAngAll", &_closeJetAngAll);
    outputTree->SetBranchAddress("_ptRelAll", &_ptRelAll);
    
    outputTree->SetBranchAddress("_closeJetPtAllMC", &_closeJetPtAllMC);
    outputTree->SetBranchAddress("_closeJetPtAllstatus", &_closeJetPtAllstatus);
    outputTree->SetBranchAddress("_partonIdMatched", &_partonIdMatched);
    outputTree->SetBranchAddress("_sameParton", &_sameParton);
	
	outputTree->SetBranchAddress("_n_PV", &_n_PV);
    
    outputTree->SetBranchAddress("_met", &_met);
    outputTree->SetBranchAddress("_met_phi", &_met_phi);
    outputTree->SetBranchAddress("HT", &HT);
    
    
    outputTree->SetBranchAddress("_mompt", &_mompt);
    outputTree->SetBranchAddress("_momphi", &_momphi);
    outputTree->SetBranchAddress("_mometa", &_mometa);
    outputTree->SetBranchAddress("_mompdg", &_mompdg);
    
    outputTree->SetBranchAddress("_n_bJets", &_n_bJets);
    outputTree->SetBranchAddress("_n_Jets", &_n_Jets);
	outputTree->SetBranchAddress("_n_JetsAll", &_n_JetsAll);
    outputTree->SetBranchAddress("_bTagged", &_bTagged);
    outputTree->SetBranchAddress("_jetEta", &_jetEta);
    outputTree->SetBranchAddress("_jetPhi", &_jetPhi);
    outputTree->SetBranchAddress("_jetPt", &_jetPt);
	outputTree->SetBranchAddress("_jetEtaAll", &_jetEtaAll);
    outputTree->SetBranchAddress("_jetPhiAll", &_jetPhiAll);
    outputTree->SetBranchAddress("_jetPtAll", &_jetPtAll);
    outputTree->SetBranchAddress("_csv", &_csv);
	outputTree->SetBranchAddress("_closeIndex", &_closeIndex);
    
    
    long nEntries = outputTree->GetEntries();
    
    //std::cout<<"Entries "<<nEntries<<std::endl;
	//if(nEntries > 50000)  
	//   nEntries = 50000;
		
	//FILE *dump = fopen("output.txt","w");
	
	
	////////////////
	int SRs[4][3][3] = {{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}}};//[sigReg][lepton event type][lepton pt type]
	
	bool verbose = false;
    	for (long it=0; it!=nEntries; ++it) {
    
		
        outputTree->GetEntry(it);
		
		bool PassSingleLep = false;
        
        if (it%10000 == 0)
            std::cout<<'.'<<std::flush;//std::cout<<((TLorentzVector*)_leptonP4->At(0))->Pt()<<std::endl;
		
		
			
		

		//////////////////////////////////////////////
		//// Select leptons by type(loose, tight) ////
		//////////////////////////////////////////////
		
		std::vector<int> tight, looseout;
		for(int i=0;i<_nLeptons;i++){
		
			if(_miniIsolation[i] < 0.4  && _ipPV[i] < 0.05 && fabs(_ipZPV[i]) < 0.1 && _originReduced[i] != 0){//remove lmva for sync stuff
				
				
				if(_eventNb == 110738 && _lumiBlock == 1108)
					std::cout<<"\nloose lepton "<<_pdgids[i]<<", pt = "<<((TLorentzVector*)_leptonP4->At(i))->Pt()<<",eta = "<<((TLorentzVector*)_leptonP4->At(i))->Eta()<<
					           ", phi = "<<((TLorentzVector*)_leptonP4->At(i))->Phi()<<", mIso = "<<_miniIsolation[i]<<
							   ", isoEA = "<<_miniIsolation[i]<<", d0 = "<<_ipPV[i]<<
							   ", dz = "<<_ipZPV[i]<<", SIP = "<<_3dIPsig[i]<<", MVAVal = "<<_MVAVal[i]<<
							   ", missingHits = "<<_missingHits[i]<<"\n";
				
				bool mvaT = false;
							   
				if(fabs(((TLorentzVector*)_leptonP4->At(i))->Eta()) < 0.8){
					if(_MVAVal[i] > -0.11)
						mvaT = true;
				}
				else if(fabs(((TLorentzVector*)_leptonP4->At(i))->Eta()) < 1.479){
					if(_MVAVal[i] > -0.35)
						mvaT = true;
				}
				else{
					if(_MVAVal[i] > -0.55)
						mvaT = true;
				}
							   
				
				if(mvaT)			
					looseout.push_back(i);
				
			}
			
			if(_istightMVANoIsoSIP[i] && ((TLorentzVector*)_leptonP4->At(i))->Pt() > 10 
					&& _3dIPsig[i] < 4 && _originReduced[i] == 0
					&& _ipPV[i] < 0.05 && fabs(_ipZPV[i]) < 0.1){
					
					
					
				if(_eventNb == 110738 && _lumiBlock == 1108)
					std::cout<<"\nTight lepton "<<_pdgids[i]<<", pt = "<<((TLorentzVector*)_leptonP4->At(i))->Pt()<<",eta = "<<((TLorentzVector*)_leptonP4->At(i))->Eta()<<
					           ", phi = "<<((TLorentzVector*)_leptonP4->At(i))->Phi()<<", mIso = "<<_miniIsolation[i]<<
							   ", isoEA = "<<_miniIsolation[i]<<", d0 = "<<_ipPV[i]<<
							   ", dz = "<<_ipZPV[i]<<", SIP = "<<_3dIPsig[i]<<", pTratio = "<<(((TLorentzVector*)_leptonP4->At(i))->Pt())/_closeJetPtAll[i]<<
							   ", pTrel = "<<_ptRelAll[i]<<", closest jet angle = "<<_closeJetAngAll[i]<<", pt = "<<_closeJetPtAll[i]<<"\n";
				
				
				float relCut = 7.0;
				float ratioCut = 0.7;
				float isoCut = 0.1;
				
				if(fabs(_pdgids[i]) == 13){
					ratioCut = 0.68;
					isoCut = 0.14;
					relCut = 6.7;
				}
				
				bool RatioRel = false, IsoCut = false;
				
				if(_closeJetAngAll[i] > 0.4)
					RatioRel = true;
				else if( ( (((TLorentzVector*)_leptonP4->At(i))->Pt())/_closeJetPtAll[i] > ratioCut || _ptRelAll[i] > relCut ) && _closeJetAngAll[i] < 0.4 )
					RatioRel = true;
					
				if( _miniIsolation[i] < isoCut )
					IsoCut = true;
				
				
				if(RatioRel && IsoCut){
					tight.push_back(i);
					if(_eventNb == 110738 && _lumiBlock == 1108)
						std::cout<<"Passes Tight\n\n";
					
				}
			
			
			}
				
		
		
		}
		if(_eventNb == 110738 && _lumiBlock == 1108)
			std::cout<<"loose size = "<<looseout.size()<<" and tight size = "<<tight.size()<<"\n";
	
		if(tight.size() < 2) continue;// or < 1 for SB selection
		
		///////////////////////////
		//// Make charge pairs ////
		///////////////////////////
		
		std::vector<std::pair<int,int>> Pairs;
		
		for(unsigned int i=0;i<tight.size();i++){
			for(unsigned int j=0;j<looseout.size();j++){
		
				std::pair<int,int> tmp;
				if(_charges[tight[i]] == _charges[looseout[j]] && tight[i] != looseout[j]){
				
					tmp.first = tight[i];
					tmp.second = looseout[j];
					Pairs.push_back(tmp);
				
				}
			}
		}
		
		if(verbose) std::cout<<"before pairs.size()\n";
		
		if(!Pairs.size())
		      continue;
			 
		if(_eventNb == 110738 && _lumiBlock == 1108)
			std::cout<<"At least one SS pair\n";
		
		/////////////////////////////////
		//// Low Mass Resonance Veto ////
		/////////////////////////////////
		std::vector<int> GPairIndex, GPairType;
		
		
		for(unsigned int i=0;i<Pairs.size();i++){
		
			bool LMV = false;
			bool GV = false;
			bool ZV = false;
			
			TLorentzVector both = *((TLorentzVector*)_leptonP4->At(Pairs[i].first)) + *((TLorentzVector*)_leptonP4->At(Pairs[i].second));
			if(both.M() < 8) 
				LMV = true;
		
			for(unsigned int j=0;j<looseout.size();j++){
			
				
				if(looseout[j] == Pairs[i].first || looseout[j] == Pairs[i].second)
					continue;
				
				TLorentzVector b1 = *((TLorentzVector*)_leptonP4->At(looseout[j])) + *((TLorentzVector*)_leptonP4->At(Pairs[i].first));	
				TLorentzVector b2 = *((TLorentzVector*)_leptonP4->At(looseout[j])) + *((TLorentzVector*)_leptonP4->At(Pairs[i].second));	
			
				bool SS = true;
				if(_charges[Pairs[i].first]*_charges[looseout[j]] < 1)
					SS = false;
					
				bool SF[2] = {false,false};
				if(_flavors[looseout[j]] == _flavors[Pairs[i].first])
					SF[0] = true;
				
				if(_flavors[looseout[j]] == _flavors[Pairs[i].second])
					SF[1] = true;
			
					
					
				if(!SS && SF[0] && b1.M() > 76 && b1.M() < 106){
					ZV = true;
					if(_eventNb == 110738 && _lumiBlock == 1108 && ZV){
						std::cout<<"Fails ZV\n";
						std::cout<<"Tight pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].first))->Pt()<<", loose pT = "
								 <<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", InvM = "<<b1.M()<<"\n";
					}
				}
				
				if(!SS && SF[1] && b2.M() > 76 && b2.M() < 106){
					ZV = true;
					if(_eventNb == 110738 && _lumiBlock == 1108 && ZV){
						std::cout<<"Fails ZV\n";
						std::cout<<"Tight pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].second))->Pt()<<", loose pT = "
								 <<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", InvM = "<<b2.M()<<"\n";
					}
				}
				
				if(!SS && SF[0] && b1.M() < 12){
					GV = true;
					if(_eventNb == 110738 && _lumiBlock == 1108 && GV){
						std::cout<<"Fails GV\n";
						std::cout<<"Tight pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].first))->Pt()<<", loose pT = "
								 <<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", InvM = "<<b1.M()<<"\n";
					}
				}
				
				if(!SS && SF[1] && b2.M() < 12){
					GV = true;
					if(_eventNb == 110738 && _lumiBlock == 1108 && GV){
						std::cout<<"Fails GV\n";
						std::cout<<"Tight pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].second))->Pt()<<", loose pT = "
								 <<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", InvM = "<<b2.M()<<"\n";
					}
				}
			
		
			}
			
			if(_eventNb == 110738 && _lumiBlock == 1108 && LMV){
				std::cout<<"Fails LMV\n";
				std::cout<<"Pt1 = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].first))->Pt()<<", Pt2 = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].second))->Pt()<<
						   ", InvM = "<<both.M()<<"\n";
			}
				
			
		
			
			
			if(!LMV && !ZV && !GV){
				GPairIndex.push_back(i);
				GPairType.push_back(_flavors[Pairs[i].first] + _flavors[Pairs[i].second]);
			}
				
		
		}
		
		
		
		
		if(!GPairIndex.size())
		      continue;
		
		if(_eventNb == 110738 && _lumiBlock == 1108)
			std::cout<<"Low Mass Veto Pass\n";
			
		if(verbose) std::cout<<"At Least One Good Pair\n";
		
		/////////////////////////////////////////////
		//// Select best pair by type and pT sum ////
		/////////////////////////////////////////////
		
		int finalPair = 0;
		float hpairpt = 0, hpairtype = -1;
		for(unsigned int i=0;i<GPairIndex.size();i++){
		
			float tsumpt = ((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[i]].first))->Pt() + ((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[i]].second))->Pt();
		
			if(GPairType[i] > hpairtype){
				finalPair = i;
				hpairpt = tsumpt;
				hpairtype = GPairType[i];
			}
			else if((GPairType[i] == hpairtype) && (tsumpt > hpairpt)){
					finalPair = i;
					hpairpt = tsumpt;
					hpairtype = GPairType[i];
			}
			
		}
		
		/////////////////////////////
		//// Now Do Jet Cleaning ////
		/////////////////////////////
		
		if(_eventNb == 110738 && _lumiBlock == 1108)
				std::cout<<"Start Jet cleaning with "<<_n_Jets<<" Jets\n\n\n\n";
        
		float CleanedHT = 0;
		int nCJets = 0, nCBJets = 0;
		for(int i = 0 ; i < _n_Jets ;i++ ){
		
			TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],0);
			
			if(_eventNb == 110738 && _lumiBlock == 1108)
				std::cout<<"Jet "<<i<<"pt = "<<_jetPt[i]<<", eta = "<<_jetEta[i]<<", phi = "<<_jetPhi[i]<<", bjet = "<<(_csv[i] > 0.814)<<"\n";
			
			
			bool dr = false;
			for(unsigned int j=0;j<looseout.size();j++){
	
				if(((TLorentzVector *)_leptonP4->At(looseout[j]))->DeltaR( jt ) < 0.4 && ((TLorentzVector *)_leptonP4->At(looseout[j]))->Pt() > 10.0 && jt.Pt() > ( (_csv[i] > 0.814) ? 25 : 40)){
					dr = true;
					if(_eventNb == 110738 && _lumiBlock == 1108)
						std::cout<<"Fails cleaning with lepton pT "<<((TLorentzVector *)_leptonP4->At(looseout[j]))->Pt()<<", and dR = "<<((TLorentzVector *)_leptonP4->At(looseout[j]))->DeltaR( jt )<<"\n\n\n";
				}
			
			}
 

			if(dr) continue;
			
			if(_csv[i] > 0.814)
				nCBJets++;
				
			if(jt.Pt() < 40) continue;
 
        	CleanedHT += _jetPt[i];
        	nCJets++;
			
		}
			
		
		if(_eventNb == 110738 && _lumiBlock == 1108)
			std::cout<<"Njets = "<<nCJets<<", Nbjets = "<<nCBJets<<", CleanedHT = "<<CleanedHT<<", MET = "<<_met<<"\n";

		
		///////////////////////////////
		//// Set the Signal Region ////
		///////////////////////////////
		
		int SigReg = 0;
		bool B0 = false;
		
		if(CleanedHT > 80 && nCJets > 1){
		
			if(CleanedHT < 500 && _met > 30){
			
				if(nCBJets >= 0)
					B0 = true;
				
				if(nCBJets == 1)
					SigReg = 1;
				else if(nCBJets == 2)
					SigReg = 2;
				else if(nCBJets >= 3)
					SigReg = 3;
					
			}//
			else if(CleanedHT >= 500){
			
				if(nCBJets >= 0)
					B0 = true;
				
				if(nCBJets == 1)
					SigReg = 1;
				else if(nCBJets == 2)
					SigReg = 2;
				else if(nCBJets >= 3)
					SigReg = 3;
			}
		}
		
		
		
		
		if(_eventNb == 110738 && _lumiBlock == 1108 && B0)
			std::cout<<"Baseline Region Reached\n";
		
		if(B0){//got to baseline
		
			//////////////////////////////////
			//// Get FR Scale for Closure ////float fptbins[6] = {10.0,15.0,25.0,35.0,50.0,70.0};
			//////////////////////////////////float etabins[4] = {0.0,1.0,2.0,2.4};
			
			int flav = 0;
			if(fabs(_pdgids[Pairs[GPairIndex[finalPair]].second]) == 13)
				flav = 1;
			
			
			float isoCut[2] = {0.1,0.14};
			float relCut[2] = {7,6.7};
			float ratCut[2] = {0.7,0.68};
			
			float CorrectedFakePt = ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt();
		
			if(_ptRelAll[Pairs[GPairIndex[finalPair]].second] > relCut[flav])
				CorrectedFakePt = ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()*(1 + TMath::Max(_miniIsolation[Pairs[GPairIndex[finalPair]].second] - isoCut[flav],0.0));
			else
				CorrectedFakePt = TMath::Max(ratCut[flav]*_closeJetPtAll[Pairs[GPairIndex[finalPair]].second],((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt());
			
			
		
			float FRScales[2][2][5][4] = {{{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}},
									  {{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}}};//[QCD or InSitu][Flavor][pt,corrected,jet][FOX]
		
			int ebin[5] = {0,0,0,0,0};
			int pbin[5] = {0,0,0,0,0};
		
			float etas[5] = {fabs(((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Eta()),
					    	 fabs(((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Eta()),
						 	 fabs(_jetEtaAll[_closeIndex[Pairs[GPairIndex[finalPair]].second]]),
						 	 fabs(_jetEtaAll[_closeIndex[Pairs[GPairIndex[finalPair]].second]]),
						 	 fabs(_jetEtaAll[_closeIndex[Pairs[GPairIndex[finalPair]].second]])};
			float pts[5]  = {fabs(((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()),
						 	CorrectedFakePt,
						 	_closeJetPtAll[Pairs[GPairIndex[finalPair]].second],
						 	_closeJetPtAll[Pairs[GPairIndex[finalPair]].second],
						 	_closeJetPtAll[Pairs[GPairIndex[finalPair]].second]};
						 
			for(int b=0;b<5;b++){
		
		
				if(etas[b] < etabins[1])
					ebin[b] = 1;
				else if(etas[b] < etabins[2])
					ebin[b] = 2;
				else
					ebin[b] = 3;
				
			
				if(pts[b] >= fptbins[0] && pts[b] < fptbins[1])
					pbin[b] = 1;
				else if(pts[b] < fptbins[2])
					pbin[b] = 2;
				else if(pts[b] < fptbins[3])
					pbin[b] = 3;
				else if(pts[b] < fptbins[4])
					pbin[b] = 4;
				else
					pbin[b] = 5;
			}
		
		
			for(int i=0;i<2;i++){
				for(int j=0;j<2;j++){
					for(int k=0;k<5;k++){
						for(int l=0;l<4;l++){
					
							FRScales[i][j][k][l] = (FakeRates[i][j][k][l]->GetBinContent(pbin[k],ebin[k]))/(1 - FakeRates[i][j][k][l]->GetBinContent(pbin[k],ebin[k]));
							if(FakeRates[i][j][k][l]->GetBinContent(pbin[k],ebin[k]) == 1)
								FRScales[i][j][k][l] = 1;
							
					
						}
					}
				}
			}
		
			
			int ept = 0, etype = 0;
			
			
			float mTmin = _mt[Pairs[GPairIndex[finalPair]].first];
			if(_mt[Pairs[GPairIndex[finalPair]].second] < _mt[Pairs[GPairIndex[finalPair]].first])
				mTmin = _mt[Pairs[GPairIndex[finalPair]].second];
				
			
			
			int SignalRegion = GetHHSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin);
			
			if(_eventNb == 110738 && _lumiBlock == 1108)
				std::cout<<"Final Pair pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt()<<" and "<<((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()<<"\n";
			
			if(((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt() > 25)
				ept++;
			
			if(((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt() > 25)
				ept++;
		
			SRs[0][GPairType[finalPair]][ept]++;
			
			
			if(_miniIsolation[Pairs[GPairIndex[finalPair]].second] < 0.4 
					&& _istightMVANoIsoSIP[Pairs[GPairIndex[finalPair]].second]){//FO1
			
				if(!(_miniIsolation[Pairs[GPairIndex[finalPair]].second] < isoCut[flav]
						 && (_ptRelAll[Pairs[GPairIndex[finalPair]].second] > relCut[flav] 
						    ||  ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()/_closeJetPtAll[Pairs[GPairIndex[finalPair]].second] > ratCut[flav] ) ) ){///loose but not tight
			
					
					if(flav == 1)
						SRYields[ept]->Fill(GetNewSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin,ept),FRScales[0][1][1][0]);//[QCD][flavor][corrected pT][FO1]
					
			
				}
			}
			
			if(_miniIsolation[Pairs[GPairIndex[finalPair]].second] < 0.4 
					&& _istightMVANoIsoSIP_LMVA[Pairs[GPairIndex[finalPair]].second]){//FO2
			
				if(!(_miniIsolation[Pairs[GPairIndex[finalPair]].second] < isoCut[flav]
						 && (_ptRelAll[Pairs[GPairIndex[finalPair]].second] > relCut[flav] 
						    ||  ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()/_closeJetPtAll[Pairs[GPairIndex[finalPair]].second] > ratCut[flav] ) ) ){///loose but not tight
			
					
					if(flav == 0)
						SRYields[ept]->Fill(GetNewSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin,ept),FRScales[0][0][1][2]);//[QCD][flavor][corrected pT][FO2]
			
				}
			}
			
			
			/*if(_miniIsolation[Pairs[GPairIndex[finalPair]].second] < 0.4 
					&& _istightMVANoIsoSIP[Pairs[GPairIndex[finalPair]].second]){//FO3
				
				bool trial = false;
				if(_ptRelAll[Pairs[GPairIndex[finalPair]].second] > relCut[flav] ||  CorrectedFakePt/_closeJetPtAll[Pairs[GPairIndex[finalPair]].second] > ratCut[flav] )
					trial = true;
				
				
				if(trial && !(_miniIsolation[Pairs[GPairIndex[finalPair]].second] < isoCut[flav]
						 && (_ptRelAll[Pairs[GPairIndex[finalPair]].second] > relCut[flav] 
						    ||  ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()/_closeJetPtAll[Pairs[GPairIndex[finalPair]].second] > ratCut[flav] ) ) ){///loose but not tight
			
			
			
				}
			}
			
			if(_miniIsolation[Pairs[GPairIndex[finalPair]].second] < 0.4 
					&& _istightMVANoIsoSIP[Pairs[GPairIndex[finalPair]].second]){//FO4
				
				bool FO4 = false;
				if(_ptRelAll[Pairs[GPairIndex[finalPair]].second] > relCut[flav] || (_closeJetPtAll[Pairs[GPairIndex[finalPair]].second]/((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()) < ((1/ratCut[flav]) + _miniIsolation[Pairs[GPairIndex[finalPair]].second] ))
					FO4 = true;
				
				
				if(FO4 && !(_miniIsolation[Pairs[GPairIndex[finalPair]].second] < isoCut[flav]
						 && (_ptRelAll[Pairs[GPairIndex[finalPair]].second] > relCut[flav] 
						    ||  ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()/_closeJetPtAll[Pairs[GPairIndex[finalPair]].second] > ratCut[flav] ) ) ){///loose but not tight
			
			
			
				}
			}*/
			
			//SRYields[ept]->Fill(GetNewSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin,ept),scale);
			
			if(ept > 1) fprintf(EvtYield,Form("%1d %9d %12d\t%2d\t%+2d %5.1f\t%+2d %5.1f\t%d\t%2d\t%5.1f\t%6.1f\t%2d\n", _runNb, _lumiBlock, _eventNb, looseout.size(), 
																											 _pdgids[Pairs[GPairIndex[finalPair]].first], ((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt(), 
																											 _pdgids[Pairs[GPairIndex[finalPair]].second],((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt(), 
																											 nCJets, nCBJets, _met, CleanedHT, SignalRegion) );
																											 
			if(SigReg)
				SRs[SigReg][GPairType[finalPair]][ept]++;
				
		}
		
		
		
		
    }//event loop
	
	const char* type[3] = {"ee","em","mm"};
	const char*	eventpt[3] = {"LL","HL","HH"};
	int nbins[3] = {8,25,32};
	
	TString fileOutName = "FakeOutputs/Predicted/Sample+"_Yields.root";
	TFile *out = new TFile(fileOutName,"RECREATE");
	
	for(int i=0;i<3;i++){
		for(int j=1;j<nbins[i]+1;j++){
			
			SRYields[i]->SetBinContent(j,(SRYields[i]->GetBinContent(j)*scale));
	
		}
	}
	
	out->cd();
	SRYields[0]->Write();
    SRYields[1]->Write();
	SRYields[2]->Write();
	
	for(int i=0;i<4;i++){
		std::cout<<"SR "<<i<<"0\n";
		for(int j=0;j<3;j++){
			std::cout<<"Event Type "<<type[j]<<"\n";
			for(int k=0;k<3;k++){
				
				std::cout<<eventpt[k]<<" :::: "<<SRs[i][j][k]<<"\n";
			
			}
		}
	}
	
	
	
	

	
	


}
