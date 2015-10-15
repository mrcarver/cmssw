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
	
	if(HT > 1125 && met <= 300)
		return 26;
		
	if(met > 300)
		return 25;
		
	if(mTmin > 120 && met <= 300)
		return (HT > 300 ? 24 : 23);
		
	if(nbjets > 3)
		nbjets = 3;
	
	sr += nbjets*6;
	
	if(njets > 4) nj = 1;
	
	if(met > 300)
		m = 2;
	else if(met > 200)
		m = 1;
		
	if(HT > 1125)
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
	
	if(HT > 1125 && met <= 300)
		return 32;
		
	if(met > 300 && HT > 300)
		return 31;
		
	if(nbjets > 3)
		nbjets = 3;
	
	sr += nbjets*8;
	
	if(njets > 4) nj = 1;
	
	if(met > 300)
		m = 2;
	else if(met > 200)
		m = 1;
		
	if(HT > 1125)
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

void MaySyncYield()
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
   	//int _closeIndex[4];
	
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
	int _genCharge[8];
    double _PVchi2;
    double _PVerr[3];
    double _ipPV[8];
    double _ipPVerr[8];
    double _ipPVmc[8];
	
	double _selectedLeptonPt[8];
	double _selectedLeptonEta[8];
	double _selectedLeptonPhi[8];
    
    double _ipZPV[8];
    double _ipZPVerr[8];
    
    double _3dIP[8];
    double _3dIPerr[8];
    double _3dIPsig[8];
	
	int _missingHits[8];
	bool _chargeConsistent[8];
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
	
	double _EcalPFIso[8];
	double _HcalPFIso[8];
	double _TrackIso[8];
    
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
    
    double _genHT;
    double _met;
    double _met_phi;
    double HT;
	
	int _DiLepTrigs[5], _DiLepHTTrigs[3];
    
    TClonesArray* _leptonP4 = new TClonesArray("TLorentzVector", 8);
    for (int i=0; i!=8; ++i) {
        new ( (*_leptonP4)[i] ) TLorentzVector();
    }
    
    //TClonesArray* _jetP4 = new TClonesArray("TLorentzVector", 50);
    //for (int i=0; i!=50; ++i) {
    //    new ( (*_jetP4)[i] ) TLorentzVector();
    //}
    
    //TClonesArray* _jetAllP4 = new TClonesArray("TLorentzVector", 50);
    //for (int i=0; i!=50; ++i) {
    //    new ( (*_jetAllP4)[i] ) TLorentzVector();
    //}
    
  
    
   	// TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    double lumi = 10000;
	double ZZto4L = 15.4;
	
	double WZJets = 2.165;//1 file
	double WJetsToLNu = 20508.9;//5 files
	double DYJets = 6024;//2 files
	double DYHT[4] = {246.761,66.3448,8.31342,2.76733};//2 files per HT bin
	double WJetsHT[4] = {2234.91,580.068,68.4003,23.1363};//2 files per HT bin
	double T14t12 = 0.0856418;
	double T14t15 = 0.0141903;
	double TTWJets = 0.6647;//2 files
	double TTZJets = 0.8565;//2 files
	double TTJets = 809.1;//6 files
	double TTbarH = 0.5085;
	double T5Deg_mG1k_mS300_mCh285_mChi280 = 0.325388;
	double T5Deg_mG1k_mS300_mChi280_4body = 0.325388;
	double T5ttttDeg_G13k_S300_Ch285_23Body = 0.0460525;
	double T6ttWW_600_425_50 = 0.174599;
	double T6ttWW_650_150_50 = 0.107045;
	
	
	
	
	double SampleXSection[9] = {TTWJets,TTZJets,WZJets,ZZto4L,DYJets,TTJets,WJetsToLNu,T14t12,T14t15};
	TString Sample = "DYJets_HT600toInf";
	
	
	
	
	TH1F *SRYields[3];
	
	SRYields[0] = new TH1F("LL_Yields","Low-Low Signal Region Yields",8,1,9);
	SRYields[1] = new TH1F("HL_Yields","High-Low Signal Region Yields",25,1,26);
	SRYields[2] = new TH1F("HH_Yields","High-High Signal Region Yields",32,1,33);
	

	//int NumEvents[9] = {0};
	//Double_t SampleScales[9] = {0};
	//int WhichSample = 5;
	FILE *EvtYield = fopen("FuckYieldstexts.txt","w");
	FILE *EvtYieldHL = fopen("EvtYieldHL.txt","w");
	
	
	int evs = 0;
    	std::cout<<"\n\n";

	
	TH1F *TotalEvents = new TH1F("TotalEvents","Total Events",5,0,5);
	
	//TFile *hfiles = new TFile("file:MaySyncOutput.root","READ");
	TFile *hfiles = new TFile("file:TTW_sync.root","READ");//
	//TFile *hfiles = new TFile("/cms/data/store/user/t2/users/mrcarver/Fall15AnalysisTuples/MC/TTJets/Full.root","READ");//
	hfiles->cd("SyncExercise");
	TH1D* _hCounters = new TH1D("hCounter", "Events counter", 5,0,5);
	_hCounters->Read("hCounter");
	evs += _hCounters->GetBinContent(1);
	
	
	Double_t scale = lumi*DYHT[3]/evs;
	
	std::cout<<"Scale = "<<scale<<"\n";
	
	TotalEvents->SetBinContent(1,evs);
	
    TChain *outputTree=new TChain("SyncExercise/SSSyncTree");
	outputTree->Add("file:TTW_sync.root");
	//outputTree->Add("/cms/data/store/user/t2/users/mrcarver/Fall15AnalysisTuples/MC/TTJets/Full.root");
	
	outputTree->SetBranchAddress("_selectedLeptonPt", &_selectedLeptonPt);
	outputTree->SetBranchAddress("_selectedLeptonEta", &_selectedLeptonEta);
	outputTree->SetBranchAddress("_selectedLeptonPhi", &_selectedLeptonPhi);
	
	outputTree->SetBranchAddress("_DiLepHTTrigs", &_DiLepHTTrigs);
	outputTree->SetBranchAddress("_DiLepTrigs"    , &_DiLepTrigs);
	
	outputTree->SetBranchAddress("_eventNb",   &_eventNb); 
    outputTree->SetBranchAddress("_runNb",     &_runNb);    
    outputTree->SetBranchAddress("_lumiBlock", &_lumiBlock);
	outputTree->SetBranchAddress("_genHT", &_genHT);
    
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
	
	outputTree->SetBranchAddress("_EcalPFIso", &_EcalPFIso);
	outputTree->SetBranchAddress("_HcalPFIso", &_HcalPFIso);
	outputTree->SetBranchAddress("_TrackIso", &_TrackIso  );
    
    outputTree->SetBranchAddress("_index1", &_index1);
    outputTree->SetBranchAddress("_index2", &_index2);

    outputTree->SetBranchAddress("_sb", &_sb);
    outputTree->SetBranchAddress("_doubleF", &_doubleF);
    
    outputTree->SetBranchAddress("_origin", &_origin);
	outputTree->SetBranchAddress("_genCharge", &_genCharge);
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
	outputTree->SetBranchAddress("_chargeConsistent", &_chargeConsistent);
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
	outputTree->SetBranchAddress("_istightMVANoIsoSIP_LMVA", &_istightMVANoIsoSIP_LMVA);
    
    
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
    
    
    long nEntries = outputTree->GetEntries();
    
    std::cout<<"Entries "<<nEntries<<std::endl;
	//if(nEntries > 50000)  
	//   nEntries = 1000;
		
	//FILE *dump = fopen("output.txt","w");
	
	
	////////////////
	int SRs[4][3][3] = {{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}}};//[sigReg][lepton event type][lepton pt type]
	
	bool verbose = false;
    	for (long it=0; it!=nEntries; ++it) {
    
		
        outputTree->GetEntry(it);
		
		bool PassSingleLep = false;
        
        if (it%10000 == 0)
            std::cout<<'.'<<std::flush;//std::cout<<((TLorentzVector*)_leptonP4->At(0))->Pt()<<std::endl;
		
		
		//if(_genHT > 100)
		//	continue;
		/////////////////////////////
		//// Now Do Jet Cleaning ////
		/////////////////////////////used to be just before setting signal region. So far matches perfect so think its good. Leave for now though
		
		if(_eventNb == 381733 && _lumiBlock == 1153) 
				std::cout<<"Start Jet cleaning with "<<_n_Jets<<" Jets\n\n\n\n";
        
		float CleanedHT = 0;
		int nCJets = 0, nCBJets = 0;
		int LepCleanIndex = -1;
		for(int i = 0 ; i < _n_Jets ;i++ ){
		
			TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],0);
			
			if(_eventNb == 381733 && _lumiBlock == 1153)
				std::cout<<"Jet "<<i<<"pt = "<<_jetPt[i]<<", eta = "<<_jetEta[i]<<", phi = "<<_jetPhi[i]<<", bjet = "<<(_csv[i] > 0.814)<<"\n";
			
			
			bool dr = false;
			/*for(unsigned int j=0;j<looseout.size();j++){part of what was used when just before setting signal region
	
				if(((TLorentzVector *)_leptonP4->At(looseout[j]))->DeltaR( jt ) < 0.4 && ((TLorentzVector *)_leptonP4->At(looseout[j]))->Pt() > 10.0 && jt.Pt() > ( (_csv[i] > 0.814) ? 25 : 40)){
					dr = true;
					if(_eventNb == 381733 && _lumiBlock == 1153)
						std::cout<<"Fails cleaning with lepton pT "<<((TLorentzVector *)_leptonP4->At(looseout[j]))->Pt()<<", and dR = "<<((TLorentzVector *)_leptonP4->At(looseout[j]))->DeltaR( jt )<<"\n\n\n";
				}
			
			}*/
			
			for(unsigned int j=0;j<_nLeptons;j++){
				if(!i) ((TLorentzVector*)_leptonP4->At(j))->SetPtEtaPhiM(_selectedLeptonPt[j],_selectedLeptonEta[j],_selectedLeptonPhi[j],0);
	
				if(LepCleanIndex != 99999 && _3dIPsig[j] < 4 && _chargeConsistent[j] && _missingHits[j] < 1 && _miniIsolation[j] < 0.4  && _ipPV[j] < 0.05 
						&& fabs(_ipZPV[j]) < 0.1 && ((TLorentzVector *)_leptonP4->At(j))->DeltaR( jt ) < 0.4 && ((TLorentzVector *)_leptonP4->At(j))->Pt() > 10.0 
						&& jt.Pt() > ( (_csv[i] > 0.814) ? 25 : 40)){
					dr = true;
					LepCleanIndex = j;
					if(_eventNb == 381733 && _lumiBlock == 1153)
						std::cout<<"Fails cleaning with lepton pT "<<((TLorentzVector *)_leptonP4->At(j))->Pt()<<", and dR = "<<((TLorentzVector *)_leptonP4->At(j))->DeltaR( jt )<<"\n\n\n";
				}
			
			}
 

			if(dr) continue;
			
			if(_csv[i] > 0.814)
				nCBJets++;
				
			if(jt.Pt() < 40) continue;
 
        	CleanedHT += _jetPt[i];
        	nCJets++;
			
		}
			
		
		if(_eventNb == 381733 && _lumiBlock == 1153)
			std::cout<<"Njets = "<<nCJets<<", Nbjets = "<<nCBJets<<", CleanedHT = "<<CleanedHT<<", MET = "<<_met<<"\n";

		//////////////////////////////////////////////
		//// Select leptons by type(loose, tight) ////
		//////////////////////////////////////////////
		
		std::vector<int> tight, looseout;
		for(int i=0;i<_nLeptons;i++){
		
			//std::cout<<"in lep\n";
		
			//((TLorentzVector*)_leptonP4->At(i))->SetPtEtaPhiM(_selectedLeptonPt[i],_selectedLeptonEta[i],_selectedLeptonPhi[i],0);
			
			//std::cout<<"pt = "<<((TLorentzVector*)_leptonP4->At(i))->Pt()<<"\n";
		
			if(_miniIsolation[i] < 0.4  && _ipPV[i] < 0.05 && fabs(_ipZPV[i]) < 0.1){//remove looseMVA for sync running
				
				
				if(_eventNb == 381733 && _lumiBlock == 1153)
					std::cout<<"\nloose lepton "<<_pdgids[i]<<", pt = "<<((TLorentzVector*)_leptonP4->At(i))->Pt()<<",eta = "<<((TLorentzVector*)_leptonP4->At(i))->Eta()<<
					           ", phi = "<<((TLorentzVector*)_leptonP4->At(i))->Phi()<<", mIso = "<<_miniIsolation[i]<<
							   ", isoEA = "<<_miniIsolation[i]<<", d0 = "<<_ipPV[i]<<
							   ", dz = "<<_ipZPV[i]<<", SIP = "<<_3dIPsig[i]<<", MVAVal = "<<_MVAVal[i]<<
							   ", missingHits = "<<_missingHits[i]<<"\n";
				
				
				looseout.push_back(i);
				
			}
			
			if(_istightMVANoIsoSIP[i] && ((TLorentzVector*)_leptonP4->At(i))->Pt() > 10 && _3dIPsig[i] < 4  //containts info about MVA tight cuts
					//&& _originReduced[i] == 0 //&& _pdgids[i] == _genCharge[i]//prompt and same charge as gen particle
					&& _ipPV[i] < 0.05 && fabs(_ipZPV[i]) < 0.1){
					
					
					
				if(_eventNb == 381733 && _lumiBlock == 1153)
					std::cout<<"\nTight lepton "<<_pdgids[i]<<", pt = "<<((TLorentzVector*)_leptonP4->At(i))->Pt()<<",eta = "<<((TLorentzVector*)_leptonP4->At(i))->Eta()<<
					           ", phi = "<<((TLorentzVector*)_leptonP4->At(i))->Phi()<<", mIso = "<<_miniIsolation[i]<<
							   ", isoEA = "<<_miniIsolation[i]<<", d0 = "<<_ipPV[i]<<
							   ", dz = "<<_ipZPV[i]<<", SIP = "<<_3dIPsig[i]<<", pTratio = "<<(((TLorentzVector*)_leptonP4->At(i))->Pt())/_closeJetPtAll[i]<<
							   ", pTrel = "<<_ptRelAll[i]<<", closest jet angle = "<<_closeJetAngAll[i]<<", pt = "<<_closeJetPtAll[i]<<"\n";
				
				
				float relCut = 7.2;
				float ratioCut = 0.8;
				float isoCut = 0.12;
				
				if(fabs(_pdgids[i]) == 13){
					ratioCut = 0.76;
					isoCut = 0.16;
					relCut = 7.2;
				}
				
				bool RatioRel = false, IsoCut = false, IsoEmu = true;
				
				if(CleanedHT < 300)
					IsoEmu = false;
				
				if(_closeJetAngAll[i] > 0.4)
					RatioRel = true;
				else if( ( (((TLorentzVector*)_leptonP4->At(i))->Pt())/_closeJetPtAll[i] > ratioCut || _ptRelAll[i] > relCut ) && _closeJetAngAll[i] < 0.4 )
					RatioRel = true;
					
				if( _miniIsolation[i] < isoCut )
					IsoCut = true;
				
				
				if(IsoCut && _eventNb == 381733 && _lumiBlock == 1153)
					std::cout<<"passes Iso\n";
					
				if(RatioRel && _eventNb == 381733 && _lumiBlock == 1153)
					std::cout<<"passes RatioRel\n";
				
				if(_EcalPFIso[i] < 0.45 &&  _HcalPFIso[i] < 0.25 && _TrackIso[i] < 0.2)
					IsoEmu = true;
				
				if(RatioRel && IsoCut && IsoEmu){
					tight.push_back(i);
					if(_eventNb == 381733 && _lumiBlock == 1153)
						std::cout<<"Passes Tight\n\n";
					
				}
			
			
			}
				
		
		
		}
		if(_eventNb == 381733 && _lumiBlock == 1153)
			std::cout<<"loose size = "<<looseout.size()<<" and tight size = "<<tight.size()<<"\n";
	
		if(tight.size() < 2) continue;// or < 1 for SB selection
		
		///////////////////////////
		//// Make charge pairs ////
		///////////////////////////
		
		std::vector<int> SSL[2];
		for(unsigned int i=0;i<tight.size();i++){
				
				if(_charges[tight[i]] < 0)
					SSL[0].push_back(tight[i]);
				else if(_charges[tight[i]] > 0)
					SSL[1].push_back(tight[i]);	
				
		}
			
		std::vector<std::pair<int,int>> Pairs, SBPairs;
		
		if(_eventNb == 381733 && _lumiBlock == 1153) std::cout<<"made SS Pairs\n";
		
		for(unsigned int i=0;i<SSL[0].size();i++){
		  for(unsigned int j=i+1;j<SSL[0].size();j++){
		      
			  std::pair<int,int> tmp;
			  
			  
			  if(  ((TLorentzVector*)_leptonP4->At(SSL[0][i]))->Pt() >  ((TLorentzVector*)_leptonP4->At(SSL[0][j]))->Pt()   ){
			      tmp.first = SSL[0][i];
				  tmp.second = SSL[0][j];
			  }
			  else{
			  	  tmp.first = SSL[0][j];
				  tmp.second = SSL[0][i];
			  }
				
			  Pairs.push_back(tmp);
			  
		  }
		}
		
		for(unsigned int i=0;i<SSL[1].size();i++){
		  for(unsigned int j=i+1;j<SSL[1].size();j++){
		      
			  std::pair<int,int> tmp;
			  
			 
			  if(  ((TLorentzVector*)_leptonP4->At(SSL[1][i]))->Pt() >  ((TLorentzVector*)_leptonP4->At(SSL[1][j]))->Pt()   ){
			      tmp.first = SSL[1][i];
				  tmp.second = SSL[1][j];
			  }
			  else{
			  	  tmp.first = SSL[1][j];
				  tmp.second = SSL[1][i];
			  }
			 
			  Pairs.push_back(tmp);
			  
		  }
		}
		
		if(_eventNb == 381733 && _lumiBlock == 1153) std::cout<<"before pairs.size()\n";
		
		if(!Pairs.size())
		      continue;
			 
		//if(_eventNb == 381733 && _lumiBlock == 1153)
		if(_eventNb == 381733 && _lumiBlock == 1153)
			std::cout<<"At least one SS pair\n";
		
		/////////////////////////////////
		//// Low Mass Resonance Veto ////
		/////////////////////////////////
		std::vector<int> GPairIndex, GPairType;
		
		
		for(unsigned int i=0;i<Pairs.size();i++){
		
			bool LMV = false;
			bool GV = false;
			bool ZV = false;
			bool OneFake = false;
			bool ChargeFlip = false;
			
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
					if(_eventNb == 381733 && _lumiBlock == 1153 && ZV){
						std::cout<<"Fails ZV\n";
						std::cout<<"Tight pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].first))->Pt()<<", loose pT = "
								 <<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", InvM = "<<b1.M()<<"\n";
					}
				}
				
				if(!SS && SF[1] && b2.M() > 76 && b2.M() < 106){
					ZV = true;
					if(_eventNb == 381733 && _lumiBlock == 1153 && ZV){
						std::cout<<"Fails ZV\n";
						std::cout<<"Tight pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].second))->Pt()<<", loose pT = "
								 <<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", InvM = "<<b2.M()<<"\n";
					}
				}
				
				if(!SS && SF[0] && b1.M() < 12){
					GV = true;
					if(_eventNb == 381733 && _lumiBlock == 1153 && GV){
						std::cout<<"Fails GV\n";
						std::cout<<"Tight pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].first))->Pt()<<", loose pT = "
								 <<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", InvM = "<<b1.M()<<"\n";
					}
				}
				
				if(!SS && SF[1] && b2.M() < 12){
					GV = true;
					if(_eventNb == 381733 && _lumiBlock == 1153 && GV){
						std::cout<<"Fails GV\n";
						std::cout<<"Tight pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].second))->Pt()<<", loose pT = "
								 <<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", InvM = "<<b2.M()<<"\n";
					}
				}
			
		
			}
			
			if(_eventNb == 381733 && _lumiBlock == 1153 && LMV){
				std::cout<<"Fails LMV\n";
				std::cout<<"Pt1 = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].first))->Pt()<<", Pt2 = "<<((TLorentzVector*)_leptonP4->At(Pairs[i].second))->Pt()<<
						   ", InvM = "<<both.M()<<"\n";
			}
				
			
			if(_originReduced[Pairs[i].first] == 0 && _originReduced[Pairs[i].second] != 0)
				OneFake = true;
			
			if(_originReduced[Pairs[i].first] != 0 && _originReduced[Pairs[i].second] == 0)
				OneFake = true;
				
			if(_pdgids[Pairs[i].first] == _genCharge[Pairs[i].first] && _pdgids[Pairs[i].second] != _genCharge[Pairs[i].second])
				ChargeFlip = true;
				
			if(_pdgids[Pairs[i].first] != _genCharge[Pairs[i].first] && _pdgids[Pairs[i].second] == _genCharge[Pairs[i].second])
				ChargeFlip = true;
			
			
			if(!LMV && !ZV && !GV){// && OneFake && ChargeFlip){
				GPairIndex.push_back(i);
				GPairType.push_back(_flavors[Pairs[i].first] + _flavors[Pairs[i].second]);
			}
				
		
		}
		
		
		if(_eventNb == 381733 && _lumiBlock == 1153) std::cout<<"fuck\n";
		
		if(!GPairIndex.size())
		      continue;
		
		if(_eventNb == 381733 && _lumiBlock == 1153)
			std::cout<<"Low Mass Veto Pass\n";
			
		if(_eventNb == 381733 && _lumiBlock == 1153) std::cout<<"At Least One Good Pair\n";
		
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
		
		
		
		///////////////////////////////////////////
		//// Apply Trigger Flavor Matching Req ////
		///////////////////////////////////////////
		
		
		//Final pair type DE - Cross - DM
		//GPairType[finalPair]
		//FPT is needed because this is backwards in ntuplizer. Should fix but will do later
		int FPT = GPairType[finalPair];
		if(GPairType[finalPair] == 0)
			FPT = 2;
			
		if(GPairType[finalPair] == 2)
			FPT = 0;
		//
		//
		//
		//
		
		
		if(_eventNb == 381733 && _lumiBlock == 1153)
			std::cout<<"Slected best pair type of "<<GPairType[finalPair]<<" and FPT = "<<FPT<<"\n";
		
		
		//DM 0-1 //Cross 2-3 //DE 4
		//_DiLepTrigs
		
		//DM - Cross - DE
		//_DiLepHTTrigs
		
		bool PassFlavMatch = false;
		
		if(CleanedHT < 300){
		
			if(FPT < 1 && (_DiLepTrigs[0] || _DiLepTrigs[1]))
				PassFlavMatch = true;
				
			if(FPT == 1 && (_DiLepTrigs[2] || _DiLepTrigs[3]))
				PassFlavMatch = true;
				
			if(FPT > 1 && _DiLepTrigs[4])
				PassFlavMatch = true;
		
	
		}
		else{
		
			if(_DiLepHTTrigs[FPT])
				PassFlavMatch = true;
		
		}
		
		if(!PassFlavMatch)
			continue;
			
		if(_eventNb == 381733 && _lumiBlock == 1153)
			std::cout<<"Passed Flavor match\n";
		
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
		
		
		
		
		if(_eventNb == 381733 && _lumiBlock == 1153 && B0)
			std::cout<<"Baseline Region Reached\n";
		
		//int Analysis = 1;///analysis = 0 for low pt and 1 for high pt
		//int SignalRegion = SR(nCBJets,_met,nCJets,CleanedHT,Analysis);
		
		if(B0){//got to baseline
		
			
			int ept = 0, etype = 0;
			
			float mTmin = _mt[Pairs[GPairIndex[finalPair]].first];
			if(_mt[Pairs[GPairIndex[finalPair]].second] < _mt[Pairs[GPairIndex[finalPair]].first])
				mTmin = _mt[Pairs[GPairIndex[finalPair]].second];
			
			int SignalRegion = GetHHSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin);
			
			if(_eventNb == 381733 && _lumiBlock == 1153)
				std::cout<<"Final Pair pT = "<<((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt()<<" and "<<((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()<<"\n";
			
			if(((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt() > 25)
				ept++;
			
			if(((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt() > 25)
				ept++;
		
			SRs[0][GPairType[finalPair]][ept]++;
			
			//if(SRs[0][GPairType[finalPair]][ept])
			//	std::cout<<"filled something\n";
			
			
			//if(ept == 0)
			//	SRYields[ept]->Fill(GetLLSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin));
			//else if(ept == 1)
			//	SRYields[ept]->Fill(GetHLSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin));
			//else if(ept == 2)
			//	SRYields[ept]->Fill(GetHHSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin));
			
			if(_eventNb == 381733 && _lumiBlock == 1153)
				std::cout<<"Ok here and all vars = "<<nCBJets<<","<<_met<<","<<nCJets<<","<<CleanedHT<<","<<mTmin<<","<<ept<<"\n";
				
			SRYields[ept]->Fill(GetNewSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin,ept));
			
			int HLSR = GetNewSignalRegion(nCBJets, _met, nCJets, CleanedHT, mTmin,ept);
			
			
			if(_eventNb == 381733 && _lumiBlock == 1153)
				std::cout<<"Ok past the fill\n";
			
			if(ept > 1) fprintf(EvtYield,Form("%1d %9d %12d\t%2d\t%+2d %5.1f\t%+2d %5.1f\t%d\t%2d\t%5.1f\t%6.1f\t%2d\n", _runNb, _lumiBlock, _eventNb, looseout.size(), 
																											 _pdgids[Pairs[GPairIndex[finalPair]].first], ((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt(), 
																											 _pdgids[Pairs[GPairIndex[finalPair]].second],((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt(), 
																											 nCJets, nCBJets, _met, CleanedHT, SignalRegion) );
																											 
			if(ept == 1) fprintf(EvtYieldHL,Form("%1d %9d %12d\t%2d\t%+2d %5.1f\t%+2d %5.1f\t%d\t%2d\t%5.1f\t%6.1f\t%2d\n", _runNb, _lumiBlock, _eventNb, looseout.size(), 
																											 _pdgids[Pairs[GPairIndex[finalPair]].first], ((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt(), 
																											 _pdgids[Pairs[GPairIndex[finalPair]].second],((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt(), 
																											 nCJets, nCBJets, _met, CleanedHT, HLSR) );
			
			if(SigReg)
				SRs[SigReg][GPairType[finalPair]][ept]++;
				
				
			
			
		}
		
		
		
		
    }//event loop
	
	const char* type[3] = {"ee","em","mm"};
	const char*	eventpt[3] = {"LL","HL","HH"};
	int nbins[3] = {8,25,32};
	
	
    
	TString fileOutName = "Fuck_Yields.root";
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
