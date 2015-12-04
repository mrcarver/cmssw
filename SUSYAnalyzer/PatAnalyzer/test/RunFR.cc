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

	if(met < 50)
		return 0;

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

void RunFR()
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
	
    gStyle->SetPadBottomMargin(0.14);
	gStyle->SetPadRightMargin(0.04);
	gStyle->SetPadTopMargin(0.1);
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
	gStyle->SetPaintTextFormat("4.3f");
	
	TH1F *HT_Yield = new TH1F("SS_HT",";H_{T} [GeV];Entries / 40 GeV",20,0,800);
	TH1F *MET_Yield = new TH1F("SS_MET",";MET [GeV];Entries / 20 GeV",35,0,700);
	TH1F *MTmin_Yield = new TH1F("SS_MTmin",";M_{T} min [GeV];Entries / 10 GeV",20,0,200);
	TH1F *NJet_Yield = new TH1F("SS_NJ",";N jets;Entries / 1 jet",10,0,10);
	TH1F *NBJet_Yield = new TH1F("SS_NBJ",";N b-jets;Entries / 1 b-jet",5,0,5);
    
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
	double _ptRatioAll[8];
	
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
	double _weight;
    
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
	
	int _DiLepTrigs[5], _DiLepHTTrigs[3], _ControlTrigs[2][2][5];
	int _ControlTrigPrescale[2][2][5];
    
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
    double lumi = 2110.0;
	double ZZto4L = 15.4;
	
	double WZJets = 2.165;//1 file
	double WJetsToLNu = 20508.9*3;//5 files
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
	//FILE *EvtYield = fopen("FuckYieldstexts.txt","w");
	//FILE *EvtYieldHL = fopen("EvtYieldHL.txt","w");
	FILE *EventDump = fopen("ObsEventList.txt","w");
	
	TString dataSamples[3] = {"DM","DE","ME"};
	TString FRSamples[4] = {"MC/WJets_G1FO_Total.root","MC/DYJets_G1FO_Total.root","Data/DM_Total_O5v4_2p11fb_v6JEC.root","Data/DE_Total_O5v4_2p11fb_v6JEC.root"};
	TString sampleOutName[4] = {"WJets","DYJets","DoubleMu","DoubleEle"};
	double FRxsec[4] = {WJetsToLNu,DYJets,1.0,1.0};
	TString Flavs[2] = {"Ele","Mu"};
	TString Isos[2]  = {"","_Iso"};
	TString ND[2]    = {"_Denom","_Numer"};
	float pbs[6]     = {10.0,15.0,25.0,35.0,50.0,70.0};
	float ebs[2][4]  = {{0.0,0.8,1.479,2.5},{0.0,1.2,2.1,2.4}};
	double IsoCuts[2][3] = {{7.2,0.8,0.12},{7.2,0.76,0.16}};
	
	FILE *EleDump = fopen("EleDump.txt","w");
	FILE *MuDump = fopen("MuDump.txt","w");
	
	bool MakeNumDen = false;
	
	if(MakeNumDen){
	for(int samp=0;samp<4;samp++){
	
	//if(samp != 2)
	//  continue;
	
	TString fileToUse = "/cms/data/store/user/t2/users/mrcarver/Fall15AnalysisTuples/FRTuples/"+FRSamples[samp];
	bool onMC = true;
	if(samp > 1)
		onMC = false;
	
	int evs = 0;
   	std::cout<<"\n\n";

	
	TH1F *TotalEvents = new TH1F("TotalEvents","Total Events",5,0,5);
	
	TFile *hfiles = new TFile(fileToUse,"READ");//
	

	
	hfiles->cd("SyncExercise");
	TH1D* _hCounters = new TH1D("hCounter", "Events counter", 5,0,5);
	_hCounters->Read("hCounter");
	evs += _hCounters->GetBinContent(1);
	
	
	Double_t scale = lumi*FRxsec[samp]/evs;
	
	std::cout<<"Scale = "<<scale<<"\n";
	
	TotalEvents->SetBinContent(1,evs);
	
    TChain *outputTree=new TChain("SyncExercise/SSSyncTree");
	
	outputTree->Add(fileToUse);
	
	
	
	
	outputTree->SetBranchAddress("_selectedLeptonPt", &_selectedLeptonPt);
	outputTree->SetBranchAddress("_selectedLeptonEta", &_selectedLeptonEta);
	outputTree->SetBranchAddress("_selectedLeptonPhi", &_selectedLeptonPhi);
	
	outputTree->SetBranchAddress("_weight", &_weight);
	
	outputTree->SetBranchAddress("_DiLepHTTrigs", &_DiLepHTTrigs);
	outputTree->SetBranchAddress("_DiLepTrigs"    , &_DiLepTrigs);
	outputTree->SetBranchAddress("_ControlTrigs"    , &_ControlTrigs);
	
	outputTree->SetBranchAddress("_ControlTrigPrescale", &_ControlTrigPrescale);
	
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
	outputTree->SetBranchAddress("_ptRatioAll", &_ptRatioAll);
    
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
	
	
	
	TH1F *tightMTDist[2][2],*FRvsEta[2][2][2], *FRvsPt[2][2][2];
	TH2F *FRMap[2][2][2];
	
	
	
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			tightMTDist[i][j] = new TH1F(Flavs[i]+"_EWK_MT_Dist"+Isos[j],"",15,0,150);
			
			for(int k=0;k<2;k++){
			
				FRvsEta[i][j][k] = new TH1F(Flavs[i]+"_FRvsEta"+Isos[j]+ND[k],"",3,ebs[i]);
				FRvsPt[i][j][k]  = new TH1F(Flavs[i]+"_FRvsPt" +Isos[j]+ND[k],"",5,pbs);
				FRMap[i][j][k]   = new TH2F(Flavs[i]+"_FRMap"  +Isos[j]+ND[k],"",5,pbs,3,ebs[i]);
			
			}
		}
	}
    
    
    long nEntries = outputTree->GetEntries();
    
    std::cout<<"Entries "<<nEntries<<std::endl;
	//if(nEntries > 50000)  
	//   nEntries = 1000;
		
	//FILE *dump = fopen("output.txt","w");
	
	
	////////////////
	double SRs[4][3][3] = {{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}}};//[sigReg][lepton event type][lepton pt type]
	
	float mvaEtaCuts[2][3] = {{-0.7,-0.83,-0.92},{-0.155,-0.56,-0.76}};
	
	
	bool verbose = false;
    	for (long it=0; it!=nEntries; ++it) {
		
		
        outputTree->GetEntry(it);
        
        if (it%10000 == 0)
            std::cout<<'.'<<std::flush;//std::cout<<((TLorentzVector*)_leptonP4->At(0))->Pt()<<std::endl;

		//////////////////////////////////////////////
		//// Select FO Leptons ///////////////////////
		//////////////////////////////////////////////
		
		std::vector<int> tight[2], looseout, FOindex[2];//[non-iso - 0, iso - 1]
		for(int i=0;i<_nLeptons;i++){
		
			((TLorentzVector*)_leptonP4->At(i))->SetPtEtaPhiM(_selectedLeptonPt[i],_selectedLeptonEta[i],_selectedLeptonPhi[i],0);
			
			
			
			if(_eventNb == 232506989){
							   
				std::cout<<"\nTight lepton "<<_pdgids[i]<<", pt = "<<((TLorentzVector*)_leptonP4->At(i))->Pt()<<",eta = "<<((TLorentzVector*)_leptonP4->At(i))->Eta()<<
					           ", phi = "<<((TLorentzVector*)_leptonP4->At(i))->Phi()<<", mIso = "<<_miniIsolation[i]<<
							   ", isoEA = "<<_miniIsolation[i]<<", d0 = "<<_ipPV[i]<<", MVAVal = "<<_MVAVal[i]<<
							   ", dz = "<<_ipZPV[i]<<", SIP = "<<_3dIPsig[i]<<", pTratio = "<<(((TLorentzVector*)_leptonP4->At(i))->Pt())/_closeJetPtAll[i]<<", ptratio 2 = "<<_ptRatioAll[i]<<
							   ", pTrel = "<<_ptRelAll[i]<<", closest jet angle = "<<_closeJetAngAll[i]<<", pt = "<<_closeJetPtAll[i]<<", met = "<<_met<<", mt = "<<_mt[i]<<"\n";
				
			}
			
			
			if(_miniIsolation[i] < 0.4  && _ipPV[i] < 0.05 && fabs(_ipZPV[i]) < 0.1 && ((TLorentzVector*)_leptonP4->At(i))->Pt() > 10.0 
									    && _3dIPsig[i] < 4){//should be all of tight minus full iso and mva
				
				
				if(_eventNb == 232506989) std::cout<<"inside nominal FO condition\n";
				
				
				float relCut = 7.2;
				float ratioCut = 0.8;
				float isoCut = 0.12;
				
				int ebin = -1;
				if(fabs(((TLorentzVector*)_leptonP4->At(i))->Eta()) < 0.8)
					ebin = 0;
				else if(fabs(((TLorentzVector*)_leptonP4->At(i))->Eta()) < 1.479)
					ebin = 1;
				else if(fabs(((TLorentzVector*)_leptonP4->At(i))->Eta()) < 2.5)	
					ebin = 2;
				
				if(fabs(_pdgids[i]) == 13){
					ratioCut = 0.76;
					isoCut = 0.16;
					relCut = 7.2;
				}
				
				bool RatioRel = false, IsoCut = false, IsoEmu = false;
				
				if(_closeJetAngAll[i] > 0.4)
					RatioRel = true;
				else if( ( (((TLorentzVector*)_leptonP4->At(i))->Pt())/_closeJetPtAll[i] > ratioCut || _ptRelAll[i] > relCut ) && _closeJetAngAll[i] < 0.4 )
					RatioRel = true;
					
				if( _miniIsolation[i] < isoCut )
					IsoCut = true;
				
				if(_EcalPFIso[i] < 0.45 &&  _HcalPFIso[i] < 0.25 && _TrackIso[i] < 0.2)
					IsoEmu = true;
				
				
				
				if(IsoEmu && ebin > -1 && _MVAVal[i] > mvaEtaCuts[1][ebin]){
					if(_flavors[i] && _istightMVANoIsoSIP[i] && fabs(((TLorentzVector*)_leptonP4->At(i))->Eta()) < 2.4)
						FOindex[1].push_back(i);
					else if(!_flavors[i])
						FOindex[1].push_back(i);
					
				}
				
				if(_eventNb == 232506989){
				
					std::cout<<"\nebin = "<<ebin<<", mvaetacuts[0][ebin] = "<<mvaEtaCuts[0][ebin]<<"\n";
				
				}
				
				if(ebin > -1 && _MVAVal[i] > mvaEtaCuts[0][ebin]){
					
					if(_flavors[i] && _istightMVANoIsoSIP[i] && fabs(((TLorentzVector*)_leptonP4->At(i))->Eta()) < 2.4)
						FOindex[0].push_back(i);
					else if(!_flavors[i])
						FOindex[0].push_back(i);
						
						
					if(_eventNb == 232506989) std::cout<<"in FO 0 here\n";
					
				}
				
				//bool pTreq = true;
				//if(!_flavors[i] && ((TLorentzVector*)_leptonP4->At(i))->Pt() <= 15.0)
				//	pTreq = false;
				
				if(RatioRel && IsoCut && _istightMVANoIsoSIP[i]){// && pTreq){
					
					tight[0].push_back(i);
					
					if(IsoEmu)
						tight[1].push_back(i);
						
					
					if(_eventNb == 232506989) std::cout<<"in FO 1\n";
				
				}
			
			}
				
		
		
		}
	
		if(FOindex[0].size() < 1 && FOindex[1].size() < 1) continue;// or < 1 for SB selection
		
		
		if(_eventNb == 232506989) std::cout<<"met = "<<_met<<"\n";
		
		///////////////////////////////////////////////
		//// Get Tight MT Dist For EWK Subtraction ////
		///////////////////////////////////////////////
		int eleTrigIndex[2] = {1,0};
		
		for(int iso=0;iso<2;iso++){
		
			if(tight[iso].size() == 1 && FOindex[0].size() == 1){
			
				bool awayJet = false;
				for(int i = 0 ; i < _n_Jets ;i++ ){
					TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],0);
					if( _jetPt[i] > 40 && ((TLorentzVector*)_leptonP4->At(tight[iso][0]))->DeltaR(jt) > 1)
						awayJet = true;
				}		
					
				if(awayJet && _met > 30){	
					
					
					if(_flavors[tight[iso][0]] == 0 && _ControlTrigs[0][iso][eleTrigIndex[iso]] && samp != 2){
						if(onMC)
							tightMTDist[0][iso]->Fill(_mt[tight[iso][0]],scale*_weight);
						else
							tightMTDist[0][iso]->Fill(_mt[tight[iso][0]],_ControlTrigPrescale[0][iso][eleTrigIndex[iso]]);
					}
					
					
					if(_flavors[tight[iso][0]] == 1 && ((TLorentzVector*)_leptonP4->At(tight[iso][0]))->Pt() < 25 && _ControlTrigs[1][iso][0] && samp != 3){
						if(onMC)
							tightMTDist[1][iso]->Fill(_mt[tight[iso][0]],scale*_weight);
						else
							tightMTDist[1][iso]->Fill(_mt[tight[iso][0]],_ControlTrigPrescale[1][iso][0]);
					}
					
					
					if(_flavors[tight[iso][0]] == 1 && ((TLorentzVector*)_leptonP4->At(tight[iso][0]))->Pt() >= 25 && _ControlTrigs[1][iso][1] && samp != 3){
						if(onMC)
							tightMTDist[1][iso]->Fill(_mt[tight[iso][0]],scale*_weight);
						else
							tightMTDist[1][iso]->Fill(_mt[tight[iso][0]],_ControlTrigPrescale[1][iso][1]);
					}
					
				}
		
			}
		
		}
		
		////////////////////////////////////////////////////////
		//// Fill Denom & Numer of FR (1D and 2D) //////////////
		////////////////////////////////////////////////////////
		
		
		
		for(int iso=0;iso<2;iso++){
		
			if(_eventNb == 232506989) std::cout<<"iso = "<<iso<<"\n\n";
			
		
			if(FOindex[iso].size() == 1 && FOindex[0].size() == 1){
		
				bool awayJet = false;
				for(int i = 0 ; i < _n_Jets ;i++ ){
					TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],0);
					
					if(_eventNb == 232506989) std::cout<<"jet"<<i<<" pt = "<<_jetPt[i]<<", delta R = "<<((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->DeltaR(jt)<<"\n";
					
					if( _jetPt[i] > 40 && ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->DeltaR(jt) > 1)
						awayJet = true;
				}
				
				
				////Calculated Corrected Lepton Pt////
				double CorrectedPt = ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt();//IsoCuts[ele,mu][rel,rat,iso]
				if(_ptRelAll[FOindex[iso][0]] > IsoCuts[_flavors[FOindex[iso][0]]][0])
					CorrectedPt = CorrectedPt*(1+TMath::Max(0.0,_miniIsolation[FOindex[iso][0]] - IsoCuts[_flavors[FOindex[iso][0]]][2])); 
				else if(_closeJetAngAll[FOindex[iso][0]] < 0.4)
					CorrectedPt = TMath::Max(CorrectedPt,_closeJetPtAll[FOindex[iso][0]]*IsoCuts[_flavors[FOindex[iso][0]]][1]);
				else 
					CorrectedPt *= 1;
					
				if(CorrectedPt >= 70)
					CorrectedPt = 69;
					
				////Decide if lepton is in numerator////
				bool isNumer = false;
				for(unsigned int t=0;t<tight[iso].size();t++){
					if(FOindex[iso][0] == tight[iso][t])
						isNumer = true;
				}
				
				if(_eventNb == 232506989) std::cout<<"mt = "<<_mt[FOindex[iso][0]]<<", CorrectedPt = "<<CorrectedPt<<"\n";
		
				if(awayJet && _met < 20 && _mt[FOindex[iso][0]] < 20){////FR event selection satisfied 
				
				
					if(_eventNb == 232506989) std::cout<<"FR event selection satisfied\n";
					
					if(_eventNb == 232506989) std::cout<<"controltrigs[0]["<<iso<<"]["<<eleTrigIndex[iso]<<" = "<<_ControlTrigs[0][iso][eleTrigIndex[iso]]<<"\n";
		
					if(_flavors[FOindex[iso][0]] == 0 && _ControlTrigs[0][iso][eleTrigIndex[iso]] && samp != 2){
							
							if(_eventNb == 232506989) std::cout<<"electron with "<<iso<<" trigger passing\n";
						
						if(onMC){
							FRMap[0][iso][0]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),scale*_weight);
							FRvsPt[0][iso][0]->Fill(CorrectedPt														   ,scale*_weight);
							FRvsEta[0][iso][0]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,scale*_weight);
							if(isNumer){
								FRMap[0][iso][1]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),scale*_weight);
								FRvsPt[0][iso][1]->Fill(CorrectedPt														   ,scale*_weight);
								FRvsEta[0][iso][1]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,scale*_weight);
							}
						}
						else{
						
							if(!iso && CorrectedPt > 35){
							
								if(_eventNb == 232506989) std::cout<<"inside print function\n";
							
								fprintf(EleDump,Form("%12lu %6.2f %6.2f %5.2f %.2f %.2f %5.2f %5.2f %1i %7i\n",_eventNb, ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt()
																												, CorrectedPt, ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()
																												, _miniIsolation[FOindex[iso][0]]
																												, (((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt())/_closeJetPtAll[FOindex[iso][0]]
																												, _ptRelAll[FOindex[iso][0]]
																												, _MVAVal[FOindex[iso][0]], isNumer
																												, _ControlTrigPrescale[0][iso][eleTrigIndex[iso]] ));
								
							
							}
						
						
							FRMap[0][iso][0]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),_ControlTrigPrescale[0][iso][eleTrigIndex[iso]]);
							FRvsPt[0][iso][0]->Fill(CorrectedPt														   ,_ControlTrigPrescale[0][iso][eleTrigIndex[iso]]);
							FRvsEta[0][iso][0]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,_ControlTrigPrescale[0][iso][eleTrigIndex[iso]]);
							if(isNumer){
								FRMap[0][iso][1]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),_ControlTrigPrescale[0][iso][eleTrigIndex[iso]]);
								FRvsPt[0][iso][1]->Fill(CorrectedPt														   ,_ControlTrigPrescale[0][iso][eleTrigIndex[iso]]);
								FRvsEta[0][iso][1]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,_ControlTrigPrescale[0][iso][eleTrigIndex[iso]]);
							}
						}
					}
					
					
					if(_flavors[FOindex[iso][0]] == 1 && ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt() < 25 && _ControlTrigs[1][iso][0] && samp != 3){
					
						if(_eventNb == 232506989) std::cout<<"muon with "<<iso<<" trigger passing\n";
						
						if(onMC){
							FRMap[1][iso][0]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),scale*_weight);
							FRvsPt[1][iso][0]->Fill(CorrectedPt												           ,scale*_weight);
							FRvsEta[1][iso][0]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,scale*_weight);
							if(isNumer){
								FRMap[1][iso][1]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),scale*_weight);
								FRvsPt[1][iso][1]->Fill(CorrectedPt												           ,scale*_weight);
								FRvsEta[1][iso][1]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,scale*_weight);
							}
						}
						else{
						
						
							if(iso && CorrectedPt > 35){
							
								if(_eventNb == 232506989) std::cout<<"inside print function\n";
							
								fprintf(MuDump,Form("%12lu %6.2f %6.2f %5.2f %.2f %.2f %5.2f %5.2f %1i %7i\n",_eventNb, ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt()
																												, CorrectedPt, ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()
																												, _miniIsolation[FOindex[iso][0]]
																												, (((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt())/_closeJetPtAll[FOindex[iso][0]]
																												, _ptRelAll[FOindex[iso][0]]
																												, _MVAVal[FOindex[iso][0]], isNumer
																												, _ControlTrigPrescale[1][iso][0] ));
								
							
							}
						
						
							FRMap[1][iso][0]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),_ControlTrigPrescale[1][iso][0]);
							FRvsPt[1][iso][0]->Fill(CorrectedPt														   ,_ControlTrigPrescale[1][iso][0]);
							FRvsEta[1][iso][0]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,_ControlTrigPrescale[1][iso][0]);
							if(isNumer){
								FRMap[1][iso][1]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),_ControlTrigPrescale[1][iso][0]);
								FRvsPt[1][iso][1]->Fill(CorrectedPt														   ,_ControlTrigPrescale[1][iso][0]);
								FRvsEta[1][iso][1]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,_ControlTrigPrescale[1][iso][0]);
							}
						}
					}
					
					
					if(_flavors[FOindex[iso][0]] == 1 && ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt() >= 25 && _ControlTrigs[1][iso][1] && samp != 3){
					
						if(_eventNb == 232506989) std::cout<<"muon with "<<iso<<" trigger passing\n";
					
						if(onMC){
							FRMap[1][iso][0]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),scale*_weight);
							FRvsPt[1][iso][0]->Fill(CorrectedPt												           ,scale*_weight);
							FRvsEta[1][iso][0]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,scale*_weight);
							if(isNumer){
								FRMap[1][iso][1]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),scale*_weight);
								FRvsPt[1][iso][1]->Fill(CorrectedPt												           ,scale*_weight);
								FRvsEta[1][iso][1]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,scale*_weight);
							}
						}
						else{
						
							if(iso && CorrectedPt > 35){
							
								if(_eventNb == 232506989) std::cout<<"inside print function\n";
							
								fprintf(MuDump,Form("%12lu %6.2f %6.2f %5.2f %.2f %.2f %5.2f %5.2f %1i %7i\n",_eventNb, ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt()
																												, CorrectedPt, ((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()
																												, _miniIsolation[FOindex[iso][0]]
																												, (((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Pt())/_closeJetPtAll[FOindex[iso][0]]
																												, _ptRelAll[FOindex[iso][0]]
																												, _MVAVal[FOindex[iso][0]], isNumer
																												, _ControlTrigPrescale[1][iso][1] ));
								
							
							}
						
						
							FRMap[1][iso][0]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),_ControlTrigPrescale[1][iso][1]);
							FRvsPt[1][iso][0]->Fill(CorrectedPt														   ,_ControlTrigPrescale[1][iso][1]);
							FRvsEta[1][iso][0]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,_ControlTrigPrescale[1][iso][1]);
							if(isNumer){
								FRMap[1][iso][1]->Fill(CorrectedPt,fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta()),_ControlTrigPrescale[1][iso][1]);
								FRvsPt[1][iso][1]->Fill(CorrectedPt														   ,_ControlTrigPrescale[1][iso][1]);
								FRvsEta[1][iso][1]->Fill(fabs(((TLorentzVector*)_leptonP4->At(FOindex[iso][0]))->Eta())		   ,_ControlTrigPrescale[1][iso][1]);
							}
						}
					}
		
				}
		
			}
		
		}
		
		
		
    }//event loop
	
	
	
	TString fileOutName = sampleOutName[samp]+"_FROut_2p11_Test.root";
	
	TFile *out = new TFile(fileOutName,"RECREATE");
	out->cd();
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			tightMTDist[i][j]->Write();
			for(int k=0;k<2;k++){
				FRMap[i][j][k]->Write();
				FRvsPt[i][j][k]->Write();
				FRvsEta[i][j][k]->Write();
			}
		}
	}
	out->Close();
	
	
	}//sample loop
	}//make numden if
	
	
	bool doPlots = true;
	if(doPlots){
	
	float nFactor[2][2] = {{1,1},{1,1}};
	float numEve[4][2][2] = {{{0,0},{0,0}},{{0,0},{0,0}},{{0,0},{0,0}},{{0,0},{0,0}}};
	TH1F *mtdist[4][2][2];
	THStack *stack[2][2];
	TLegend *legend = new TLegend(0.3,0.5,0.5,0.7);
	TCanvas *norms[2][2];
	
	TH2F *FRMap[4][2][2][2];
	TH1F *FRvsPt[4][2][2][2], *FRvsEta[4][2][2][2];
	
	for(int i=0;i<4;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<2;k++){
				mtdist[i][j][k] = new TH1F(sampleOutName[i]+Flavs[j]+"_EWK_MT_Dist"+Isos[k],"",15,0,150);
				for(int nd=0;nd<2;nd++){
					FRvsEta[i][j][k][nd] = new TH1F(sampleOutName[i]+Flavs[j]+"_FRvsEta"+Isos[k]+ND[nd],"",3,ebs[j]);
					FRvsPt[i][j][k][nd]  = new TH1F(sampleOutName[i]+Flavs[j]+"_FRvsPt" +Isos[k]+ND[nd],"",5,pbs);
					FRMap[i][j][k][nd]   = new TH2F(sampleOutName[i]+Flavs[j]+"_FRMap"  +Isos[k]+ND[nd],"",5,pbs,3,ebs[j]);
				}
			}
		}
	}
	
	for(int samp=0;samp<4;samp++){
	
		TString fileOutName = sampleOutName[samp]+"_FROut_2p11_Test.root";
		TFile *file1 = new TFile(sampleOutName[samp]+"_FROut_2p11_Test.root","READ");
		
		for(int f=0;f<2;f++){
			for(int iso=0;iso<2;iso++){
		
				mtdist[samp][f][iso]->Read(Flavs[f]+"_EWK_MT_Dist"+Isos[iso]);
				for(int nd=0;nd<2;nd++){
					FRvsEta[samp][f][iso][nd]->Read(Flavs[f]+"_FRvsEta"+Isos[iso]+ND[nd]);
					FRvsPt[samp][f][iso][nd]->Read(Flavs[f]+"_FRvsPt" +Isos[iso]+ND[nd]);
					FRMap[samp][f][iso][nd]->Read(Flavs[f]+"_FRMap"  +Isos[iso]+ND[nd]);
				}
				mtdist[samp][f][iso]->GetXaxis()->SetTitle("M_{T} [GeV]");
				mtdist[samp][f][iso]->GetYaxis()->SetTitle("Events / 10 GeV");
				if(samp < 2){
					mtdist[samp][f][iso]->SetLineColor(samp+2);
					mtdist[samp][f][iso]->SetFillColor(samp+2);
				}
				else{
					mtdist[samp][f][iso]->SetLineColor(1);
					mtdist[samp][f][iso]->SetMarkerColor(1);
				}
				
				numEve[samp][f][iso] = mtdist[samp][f][iso]->Integral(8,11);
				
		
			}
		}
	}
	
	TCanvas *FinalFRMap[2][2];
	//TH2F *FinalFRs[2][2];
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			nFactor[i][j] = (numEve[2][i][j]+numEve[3][i][j])/(numEve[0][i][j]+numEve[1][i][j]);
			
			std::cout<<"SF "<<i<<j<<" = "<<nFactor[i][j]<<"\n";
			/*
			for(int b=1;b<16;b++){
				mtdist[0][i][j]->SetBinContent(b, mtdist[0][i][j]->GetBinContent(b)*nFactor[i][j]);
				mtdist[1][i][j]->SetBinContent(b, mtdist[1][i][j]->GetBinContent(b)*nFactor[i][j]);
			
			}
			
			stack[i][j] = new THStack();
			stack[i][j]->Add(mtdist[0][i][j],"hist");
			stack[i][j]->Add(mtdist[1][i][j],"hist");
			
			TCanvas *norms = new TCanvas(Flavs[i]+Isos[j]+"Norms","",2);
			norms->cd();
			
			mtdist[2][i][j]->GetYaxis()->SetRangeUser(0,(mtdist[2][i][j]->GetBinContent(7)+mtdist[3][i][j]->GetBinContent(7))*1.4);
			mtdist[2][i][j]->Draw();
			stack[i][j]->Draw("same");
			mtdist[2][i][j]->Draw("same");
			mtdist[3][i][j]->Draw("same");
			
			*/
			//////////////////////////
			////Do EWK Subtraction////
			//////////////////////////
			for(int nd=0;nd<2;nd++){
			
				for(int xb=1;xb<6;xb++){
					for(int yb=1;yb<4;yb++){
			
						FRMap[2][i][j][nd]->SetBinContent(xb,yb, (FRMap[2][i][j][nd]->GetBinContent(xb,yb) - (FRMap[0][i][j][1]->GetBinContent(xb,yb) + FRMap[1][i][j][1]->GetBinContent(xb,yb))*nFactor[i][j]  ) );
						FRMap[3][i][j][nd]->SetBinContent(xb,yb, (FRMap[3][i][j][nd]->GetBinContent(xb,yb) - (FRMap[0][i][j][1]->GetBinContent(xb,yb) + FRMap[1][i][j][1]->GetBinContent(xb,yb))*nFactor[i][j]  ) );
					}
				
					//do vs pT
				
					//if(xb < 4){//do vs eta
					//}
				}
				
				
			
			}
			
			//FinalFRs[i][j] = new TH2F(Flavs[i]+Isos[j],"",5,pbs,3,ebs[i]);
			FRMap[2][i][j][1]->Divide(FRMap[2][i][j][1],FRMap[2][i][j][0],1.,1.,"cp");
			FRMap[3][i][j][1]->Divide(FRMap[3][i][j][1],FRMap[3][i][j][0],1.,1.,"cp");
			
			FinalFRMap[i][j] = new TCanvas(Flavs[i]+Isos[j],"",2);
			FinalFRMap[i][j]->cd();
			FinalFRMap[i][j]->SetGridy(false);
			FinalFRMap[i][j]->SetGridx(false);
			FRMap[2][i][j][1]->SetMarkerColor(1);
			FRMap[2][i][j][1]->SetMarkerSize(2);
			FRMap[3][i][j][1]->SetMarkerColor(1);
			FRMap[3][i][j][1]->SetMarkerSize(2);
			FRMap[2][i][j][1]->GetZaxis()->SetRangeUser(0,0.5);
			FRMap[3][i][j][1]->GetZaxis()->SetRangeUser(0,0.5);
			if(!i)
				FRMap[3][i][j][1]->Draw("colz text e");
			else
				FRMap[2][i][j][1]->Draw("colz text e");
 			
			
		}
	}
	
	
	
	
	
	
	
	}
	
	const char* type[3] = {"ee","em","mm"};
	const char*	eventpt[3] = {"LL","HL","HH"};
	int nbins[3] = {8,26,32};

	


}
