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
#include <TLatex.h>
#include <TPaveStats.h>
#include <TString.h>

int SR(int nbjets, float met, int njets, float HT, int analysis){///analysis = 0 for low pt and 1 for high pt

	///sr[nbjets][met][njets][HT]
	int sr[3][3][3][3] = {{{{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},{{-1,-1,-1},{-1,1,2},{-1,3,4}},{{-1,-1,-1},{-1,5,6},{-1,7,8}}},
						  {{{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},{{-1,-1,-1},{-1,11,12},{-1,13,14}},{{-1,-1,-1},{-1,15,16},{-1,17,18}}},
					 	  {{{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},{{-1,-1,-1},{-1,21,22},{-1,23,24}},{{-1,-1,-1},{-1,25,26},{-1,27,28}}}};


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
	else if(nbjets >= 2)
		nbj = 2;
		
	return sr[nbj][m][nj][h];
}

int pdgid(int flavor, int charge){

	int out = 0;
	if(!flavor)
		out = 11;
	else
		out = 13;

	return out*charge*-1;
}


double deltaPhi(TLorentzVector* lepton, double met, double metPhi) {
    double W_px = lepton->Px() + met*TMath::Cos(metPhi);
    double W_py = lepton->Py() + met*TMath::Sin(metPhi);
    double W_phi = TMath::ATan2(W_py,W_px);
    return(fabs(W_phi -lepton->Phi()));
}

void PureDileptonEffs()
{
    //int leptonStudy = 1; //1=muon, 0-electron

    gStyle->SetPaintTextFormat("4.2f");
	gStyle->SetOptFit(1);
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
	gStyle->SetPadRightMargin(0.16);
	gStyle->SetPadTopMargin(0.04);
	gStyle->SetPadLeftMargin(0.11);
    
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
    //int _n_JetsAll;
    
    //double _jetEtaAll[100];
    //double _jetPhiAll[100];
    //double _jetPtAll[100];
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
    double _isolation[8];
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
    double _closeJetPtAll[4];
    
    double _closeJetAngAll[4];
    double _ptRelAll[4];
    double _closeJetPtAllMC[4];
    double _closeJetPtAllstatus[4];
    int _partonIdMatched[4];
    bool _sameParton[4];
    
    bool _isloose[8];
    bool _istight[8];
	bool _istightNoIso[8];
    
    int _n_PV;
    
    //int _n_electrons;
    //int _n_muons;
   // int _n_taus;
    
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;
	
	
	int _Mu17Mu8;
	int _Mu17TkMu8;
	int _Ele23Ele12;
	int _Mu23Ele12;
	int _Ele23Mu8;
	
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
	
	
	int _nMu17Mu8Objs;
	int _nMu17TkMu8Objs;
	int _nEle23Ele12Objs;
	int _nMu23Ele12Objs;
	int _nEle23Mu8Objs;
	
	TClonesArray* _Mu17Mu8Objs = new TClonesArray("TLorentzVector", 30);
    for (int i=0; i!=30; ++i) {
        new ( (*_Mu17Mu8Objs)[i] ) TLorentzVector();
    }
	
	TClonesArray* _Mu17TkMu8Objs = new TClonesArray("TLorentzVector", 30);
    for (int i=0; i!=30; ++i) {
        new ( (*_Mu17TkMu8Objs)[i] ) TLorentzVector();
    }
	
	TClonesArray* _Ele23Ele12Objs = new TClonesArray("TLorentzVector", 30);
    for (int i=0; i!=30; ++i) {
        new ( (*_Ele23Ele12Objs)[i] ) TLorentzVector();
    }
	
	TClonesArray* _Mu23Ele12Objs = new TClonesArray("TLorentzVector", 30);
    for (int i=0; i!=30; ++i) {
        new ( (*_Mu23Ele12Objs)[i] ) TLorentzVector();
    }
	
	TClonesArray* _Ele23Mu8Objs = new TClonesArray("TLorentzVector", 30);
    for (int i=0; i!=30; ++i) {
        new ( (*_Mu17Mu8Objs)[i] ) TLorentzVector();
    }
    
    TClonesArray* _jetP4 = new TClonesArray("TLorentzVector", 50);
    for (int i=0; i!=50; ++i) {
        new ( (*_jetP4)[i] ) TLorentzVector();
    }
    
    TClonesArray* _jetAllP4 = new TClonesArray("TLorentzVector", 50);
    for (int i=0; i!=50; ++i) {
        new ( (*_jetAllP4)[i] ) TLorentzVector();
    }

	const char* firse[2] = {"Zero/","One/"};
	const char* tt25file = "/lfs/scratch/mrcarver/TT25_Ntuples_TriggerObjects/";
	
	enum Triggers{  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1,//0
					HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1,//1
					HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1,//2
					HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1,//3
					HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1};//4
					
	const char* PureTrigs[5] = {  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1",//0
							   "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1",//1
							   "HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1",//2
							   "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1",//3
							   "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1"};
							   
	TString PureTrigsFnames[5] = {  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1",//0
							   "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1",//1
							   "HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1",//2
							   "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1",//3
							   "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1"};
							   
	const char* PureTrigs1d[7] = {   "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1",//0
							   		 "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1",//1
							  		 "HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1",//2
							  		 "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1_Mu",//3
							   		 "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1_Mu",
									 "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1_Ele",//3
							   		 "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1_Ele"};
									 
	TString PureTrigs1dFnames[7] = {   "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1",//0
							   		 "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1",//1
							  		 "HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1",//2
							  		 "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1_Mu",//3
							   		 "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1_Mu",
									 "HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1_Ele",//3
							   		 "HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1_Ele"};
	
	const char* titles[5] = { "M17M8","M17TM8","E23E12","M23E12","E23M8"}; 
	const char* frac[2] = {"Denom","Numer"};
	const char* axis[5][2] = {{"Trailing","Leading"},{"Trailing","Leading"},{"Trailing","Leading"},{"Electron","Muon"},{"Electron","Muon"}};
	const char* titles1d[7] = {"M17M81d","M17TM81d","E23E121d","M23E12_Mu1d","E23M8_Mu1d","M23E12_Ele1d","E23M8_Ele1d"};
	const char* axis1d[7] = {"Sub Leading Muon","Sub Leading Muon","Sub Leading Electron","Muon","Muon","Electron","Electron"};

	TH2F *PureDiLepton[5][2];
	TH1F *PureDiLepton1D[7][2], *PureDiLepton1DIso[7][2];
	
	for(int i=0;i<2;i++){
		for(int j=0;j<5;j++){
	
			PureDiLepton[j][i] = new TH2F(Form("%s_%s",titles[j],frac[i]),Form("%s;%s p_{T} [GeV];%s p_{T} [GeV];efficiency",PureTrigs[j],axis[j][1],axis[j][0]),10,0,50,10,0,50);
	
	
		}
		
		for(int j=0;j<7;j++){
	
			PureDiLepton1D[j][i] = new TH1F(Form("%s_%s",titles1d[j],frac[i]),Form("%s;%s p_{T} [GeV];efficiency",PureTrigs1d[j],axis1d[j]),35,0,70);
			PureDiLepton1DIso[j][i] = new TH1F(Form("%s_%s_Iso",titles1d[j],frac[i]),Form("%s;%s Relative Iso;efficiency",PureTrigs1d[j],axis1d[j]),40,0,2);
		
	
		}
	}
	

	
   
	for(int t=0;t<2;t++){
	for(int f=0;f<3;f++){
    
    std::stringstream ss;
	ss << tt25file << firse[t] <<"Job_"<< f+1 << ".root";
	
	std::cout<<ss.str().c_str()<<"\n";
    TChain *outputTree = new TChain("SyncExercise/SSSyncTree");
	
	outputTree->Add(ss.str().c_str());
	
	outputTree->SetBranchAddress("_leptonP4", &_leptonP4);
	outputTree->SetBranchAddress("_Mu17Mu8Objs", &_Mu17Mu8Objs);
	outputTree->SetBranchAddress("_Mu17TkMu8Objs", &_Mu17TkMu8Objs);
	outputTree->SetBranchAddress("_Ele23Ele12Objs", &_Ele23Ele12Objs);
	outputTree->SetBranchAddress("_Mu23Ele12Objs", &_Mu23Ele12Objs);
	outputTree->SetBranchAddress("_Ele23Mu8Objs", &_Ele23Mu8Objs);
	

	outputTree->SetBranchAddress("_nMu17Mu8Objs", &_nMu17Mu8Objs);
	outputTree->SetBranchAddress("_nMu17TkMu8Objs", &_nMu17TkMu8Objs);
	outputTree->SetBranchAddress("_nEle23Ele12Objs", &_nEle23Ele12Objs);
	outputTree->SetBranchAddress("_nMu23Ele12Objs", &_nMu23Ele12Objs);
	outputTree->SetBranchAddress("_nEle23Mu8Objs", &_nEle23Mu8Objs);
	
	outputTree->SetBranchAddress("_eventNb",   &_eventNb); 
    outputTree->SetBranchAddress("_runNb",     &_runNb);    
    outputTree->SetBranchAddress("_lumiBlock", &_lumiBlock);
	
	outputTree->SetBranchAddress("_Mu17Mu8",   &_Mu17Mu8); 
    outputTree->SetBranchAddress("_Mu17TkMu8",     &_Mu17TkMu8);    
    outputTree->SetBranchAddress("_Ele23Ele12", &_Ele23Ele12);
	outputTree->SetBranchAddress("_Mu23Ele12",     &_Mu23Ele12);    
    outputTree->SetBranchAddress("_Ele23Mu8", &_Ele23Mu8);
    
    outputTree->SetBranchAddress("_nLeptons", &_nLeptons);
    outputTree->SetBranchAddress("_flavors", &_flavors);
    outputTree->SetBranchAddress("_charges", &_charges);
	outputTree->SetBranchAddress("_pdgids", &_pdgids);
    outputTree->SetBranchAddress("_isolation", &_isolation);
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
    outputTree->SetBranchAddress("_istight", &_istight);
	outputTree->SetBranchAddress("_istightNoIso", &_istightNoIso);
    
    
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
    outputTree->SetBranchAddress("_bTagged", &_bTagged);
    outputTree->SetBranchAddress("_jetEta", &_jetEta);
    outputTree->SetBranchAddress("_jetPhi", &_jetPhi);
    outputTree->SetBranchAddress("_jetPt", &_jetPt);
    outputTree->SetBranchAddress("_csv", &_csv);
    
    
    long nEntries = outputTree->GetEntries();
    
    std::cout<<"Entries "<<nEntries<<std::endl;
	//if(nEntries > 50000)  
	//	nEntries = 50000;
	
	
    for (long it=0; it!=nEntries; ++it) {
    
		
        outputTree->GetEntry(it);
        
        if (it%10000 == 0)
            std::cout<<'.'<<std::flush;
		
		
		/////////////////////////////////////
		//// Select Best Pair of Leptons ////
		/////////////////////////////////////
		std::vector<int> Smus, Sels, Smus_ni, Sels_ni;
		for(int i=0;i<_nLeptons;i++){
		

			if(_istight[i] && fabs(_pdgids[i]) == 13)
				Smus.push_back(i);
				
			if(_istight[i] && fabs(_pdgids[i]) == 11)
				Sels.push_back(i);
				
			if(_istightNoIso[i] && fabs(_pdgids[i]) == 13)
				Smus_ni.push_back(i);
				
			if(_istightNoIso[i] && fabs(_pdgids[i]) == 11)
				Sels_ni.push_back(i);
				
		
		}
		
		
		if(Smus.size() == 2){
		
			     //first muon[two triggers],  secondmuon[two triggers]
			bool match1[2] = {false,false}, match2[2] = {false,false};
			
			for(int m1=0;m1<_nMu17Mu8Objs;m1++){
			
			
				if(((TLorentzVector*)_leptonP4->At(Smus[0]))->DeltaR(*((TLorentzVector*)_Mu17Mu8Objs->At(m1))) < 0.4 && ((TLorentzVector*)_Mu17Mu8Objs->At(m1))->Pt() > 6)
					match1[0] = true;
					
				if(((TLorentzVector*)_leptonP4->At(Smus[1]))->DeltaR(*((TLorentzVector*)_Mu17Mu8Objs->At(m1))) < 0.4 && ((TLorentzVector*)_Mu17Mu8Objs->At(m1))->Pt() > 6)
					match2[0] = true;
		
			}
			
			for(int m1=0;m1<_nMu17TkMu8Objs;m1++){
			
			
				if(((TLorentzVector*)_leptonP4->At(Smus[0]))->DeltaR(*((TLorentzVector*)_Mu17TkMu8Objs->At(m1))) < 0.4 && ((TLorentzVector*)_Mu17TkMu8Objs->At(m1))->Pt() > 6)
					match1[1] = true;
					
				if(((TLorentzVector*)_leptonP4->At(Smus[1]))->DeltaR(*((TLorentzVector*)_Mu17TkMu8Objs->At(m1))) < 0.4 && ((TLorentzVector*)_Mu17TkMu8Objs->At(m1))->Pt() > 6)
					match2[1] = true;
		
			}
			
		
		
			int low=0, high = 1;
			if(((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt() > ((TLorentzVector*)_leptonP4->At(Smus[1]))->Pt()){
				low = 1;
				high = 0;
			}
			
			PureDiLepton[Triggers::HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1][0]->Fill(((TLorentzVector*)_leptonP4->At(Smus[high]))->Pt(),((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt());
			if(_Mu17Mu8 && match1[0] && match2[0])
				PureDiLepton[Triggers::HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1][1]->Fill(((TLorentzVector*)_leptonP4->At(Smus[high]))->Pt(),((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt());
			
			PureDiLepton[Triggers::HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1][0]->Fill(((TLorentzVector*)_leptonP4->At(Smus[high]))->Pt(),((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt());
			if(_Mu17TkMu8 && match1[1] && match2[1])
				PureDiLepton[Triggers::HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1][1]->Fill(((TLorentzVector*)_leptonP4->At(Smus[high]))->Pt(),((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt());
			
		
			
			if(((TLorentzVector*)_leptonP4->At(Smus[high]))->Pt() > 35){
				
				PureDiLepton1D[0][0]->Fill(((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt());
				PureDiLepton1D[1][0]->Fill(((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt());
				
				if(_Mu17Mu8 && match1[0] && match2[0]){
					PureDiLepton1D[0][1]->Fill(((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt());
					
					//if(((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt() < 3)
						//std::cout<<"muon < 3 GeV event "<<_eventNb<<"  run "<<_runNb<<"\n";
					
				}
			
				
				if(_Mu17TkMu8 && match1[1] && match2[1])
					PureDiLepton1D[1][1]->Fill(((TLorentzVector*)_leptonP4->At(Smus[low]))->Pt());
			}
		
		}
		
		
		if(Smus.size() == 1 && Sels.size() == 1){
		
			bool matchm[2] = {false,false}, matche[2] = {false,false};
			
			
			for(int t1=0;t1<_nMu23Ele12Objs;t1++){
			
			
				if(((TLorentzVector*)_leptonP4->At(Smus[0]))->DeltaR(*((TLorentzVector*)_Mu23Ele12Objs->At(t1))) < 0.4 && ((TLorentzVector*)_Mu23Ele12Objs->At(t1))->Pt() > 6)
					matchm[0] = true;
					
				if(((TLorentzVector*)_leptonP4->At(Sels[0]))->DeltaR(*((TLorentzVector*)_Mu23Ele12Objs->At(t1))) < 0.4 && ((TLorentzVector*)_Mu23Ele12Objs->At(t1))->Pt() > 6)
					matche[0] = true;
			
			
			}
			
			for(int t1=0;t1<_nEle23Mu8Objs;t1++){
			
			
				if(((TLorentzVector*)_leptonP4->At(Smus[0]))->DeltaR(*((TLorentzVector*)_Ele23Mu8Objs->At(t1))) < 0.4 && ((TLorentzVector*)_Ele23Mu8Objs->At(t1))->Pt() > 6)
					matchm[1] = true;
					
				if(((TLorentzVector*)_leptonP4->At(Sels[0]))->DeltaR(*((TLorentzVector*)_Ele23Mu8Objs->At(t1))) < 0.4 && ((TLorentzVector*)_Ele23Mu8Objs->At(t1))->Pt() > 6)
					matche[1] = true;
			
			
			}
			
			PureDiLepton[Triggers::HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1][0]->Fill(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt(),((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt());
			if(_Mu23Ele12 && matchm[0] && matche[0])
				PureDiLepton[Triggers::HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1][1]->Fill(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt(),((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt());
			
			PureDiLepton[Triggers::HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1][0]->Fill(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt(),((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt());
			if(_Ele23Mu8 && matchm[1] && matche[1])
				PureDiLepton[Triggers::HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1][1]->Fill(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt(),((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt());
			
			if(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt() > 35){
				PureDiLepton1D[3][0]->Fill(((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt());
				PureDiLepton1D[4][0]->Fill(((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt());
			
				if(_Mu23Ele12 && matchm[0] && matche[0])
					PureDiLepton1D[3][1]->Fill(((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt());
					
				if(_Ele23Mu8 && matchm[1] && matche[1])
					PureDiLepton1D[4][1]->Fill(((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt());
			
			}
			
			if(((TLorentzVector*)_leptonP4->At(Smus[0]))->Pt() > 35){
				PureDiLepton1D[5][0]->Fill(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt());
				PureDiLepton1D[6][0]->Fill(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt());
			
				if(_Mu23Ele12 && matchm[0] && matche[0])
					PureDiLepton1D[5][1]->Fill(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt());
				
				if(_Ele23Mu8 && matchm[1] && matche[1])
					PureDiLepton1D[6][1]->Fill(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt());
			
			}

		
		}
		
		if(Sels.size() == 2){
		
		
			bool match1 = false, match2 = false;
			
			
			for(int e1=0;e1<_nEle23Ele12Objs;e1++){
			
				
				if(((TLorentzVector*)_leptonP4->At(Sels[0]))->DeltaR(*((TLorentzVector*)_Ele23Ele12Objs->At(e1))) < 0.4 && ((TLorentzVector*)_Ele23Ele12Objs->At(e1))->Pt() > 6)
					match1 = true;
					
					
				if(((TLorentzVector*)_leptonP4->At(Sels[1]))->DeltaR(*((TLorentzVector*)_Ele23Ele12Objs->At(e1))) < 0.4 && ((TLorentzVector*)_Ele23Ele12Objs->At(e1))->Pt() > 6)
					match2 = true;
			
			
			
			}
		
			int low = 0, high = 1;
			if(((TLorentzVector*)_leptonP4->At(Sels[0]))->Pt() > ((TLorentzVector*)_leptonP4->At(Sels[1]))->Pt()){
				low = 1;
				high = 0;
			}
			
			PureDiLepton[Triggers::HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1][0]->Fill(((TLorentzVector*)_leptonP4->At(Sels[high]))->Pt(),((TLorentzVector*)_leptonP4->At(Sels[low]))->Pt());
			if(_Ele23Ele12 && match1 && match2)
				PureDiLepton[Triggers::HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1][1]->Fill(((TLorentzVector*)_leptonP4->At(Sels[high]))->Pt(),((TLorentzVector*)_leptonP4->At(Sels[low]))->Pt());
			
			
			if(((TLorentzVector*)_leptonP4->At(Sels[high]))->Pt() > 35){
				PureDiLepton1D[2][0]->Fill(((TLorentzVector*)_leptonP4->At(Sels[low]))->Pt());
				if(_Ele23Ele12 && match1 && match2){
					PureDiLepton1D[2][1]->Fill(((TLorentzVector*)_leptonP4->At(Sels[low]))->Pt());
					
					//if(((TLorentzVector*)_leptonP4->At(Sels[low]))->Pt() < 3)
					//	std::cout<<"ele < 3 GeV event "<<_eventNb<<" run "<<_runNb<<"\n";
					
				}
					
			}
			
		
		}
		
		
		if(Smus_ni.size() == 2){
		
		
			int low = 0, high = 1;
			if(((TLorentzVector*)_leptonP4->At(Smus_ni[0]))->Pt() > ((TLorentzVector*)_leptonP4->At(Smus_ni[1]))->Pt()){
				low = 1;
				high = 0;
			}
			
			
			if(((TLorentzVector*)_leptonP4->At(Smus_ni[low]))->Pt() > 35 && ((TLorentzVector*)_leptonP4->At(Smus_ni[high]))->Pt() > 35 && _istight[Smus_ni[high]]){
				
				PureDiLepton1DIso[0][0]->Fill(_isolation[Smus_ni[low]]);
				PureDiLepton1DIso[1][0]->Fill(_isolation[Smus_ni[low]]);
				
				if(_Mu17Mu8)
					PureDiLepton1DIso[0][1]->Fill(_isolation[Smus_ni[low]]);
			
				
				if(_Mu17TkMu8)
					PureDiLepton1DIso[1][1]->Fill(_isolation[Smus_ni[low]]);
			}
		
		
		}
		
		if(Sels_ni.size() == 2){
		
			int low = 0, high = 1;
			if(((TLorentzVector*)_leptonP4->At(Sels_ni[0]))->Pt() > ((TLorentzVector*)_leptonP4->At(Sels_ni[1]))->Pt()){
				low = 1;
				high = 0;
			}
			
			
			if(((TLorentzVector*)_leptonP4->At(Sels_ni[high]))->Pt() > 35 && ((TLorentzVector*)_leptonP4->At(Sels_ni[low]))->Pt() > 35 && _istight[Sels_ni[high]]){
				
				
				PureDiLepton1DIso[2][0]->Fill(_isolation[Sels_ni[low]]);
				if(_Ele23Ele12)
					PureDiLepton1DIso[2][1]->Fill(_isolation[Sels_ni[low]]);
					
			}
			
		
		}
		
		if(Smus_ni.size() == 1 && Sels_ni.size() == 1){
		
	
			if(((TLorentzVector*)_leptonP4->At(Sels_ni[0]))->Pt() > 35 && _istight[Sels_ni[0]] && ((TLorentzVector*)_leptonP4->At(Smus_ni[0]))->Pt() > 35){
				PureDiLepton1DIso[3][0]->Fill(_isolation[Smus_ni[0]]);
				PureDiLepton1DIso[4][0]->Fill(_isolation[Smus_ni[0]]);
			
				if(_Mu23Ele12)
					PureDiLepton1DIso[3][1]->Fill(_isolation[Smus_ni[0]]);
					
				if(_Ele23Mu8)
					PureDiLepton1DIso[4][1]->Fill(_isolation[Smus_ni[0]]);
			
			}
			
			if(((TLorentzVector*)_leptonP4->At(Smus_ni[0]))->Pt() > 35 && _istight[Smus_ni[0]] && ((TLorentzVector*)_leptonP4->At(Sels_ni[0]))->Pt() > 35){
				PureDiLepton1DIso[5][0]->Fill(_isolation[Sels_ni[0]]);
				PureDiLepton1DIso[6][0]->Fill(_isolation[Sels_ni[0]]);
			
				if(_Mu23Ele12)
					PureDiLepton1DIso[5][1]->Fill(_isolation[Sels_ni[0]]);
				
				if(_Ele23Mu8)
					PureDiLepton1DIso[6][1]->Fill(_isolation[Sels_ni[0]]);
			
			}

		
		}
		
		
    }//event loop
    
	
	
	}
	}//loops over different dirs and ntpules per sample
	
	
	
	
	for(int z=0;z<2;z++){
		for(int y=0;y<5;y++){
		
			PureDiLepton[y][z]->Sumw2();
		
		}
		for(int j=0;j<7;j++){
			PureDiLepton1D[j][z]->Sumw2();
			PureDiLepton1DIso[j][z]->Sumw2();
		}
	}
	
	
	for(int j=0;j<5;j++){
		
		PureDiLepton[j][1]->Divide(PureDiLepton[j][1],PureDiLepton[j][0],1.,1.,"B");
		
	}
	
	for(int j=0;j<7;j++){
		
		PureDiLepton1D[j][1]->Divide(PureDiLepton1D[j][1],PureDiLepton1D[j][0],1.,1.,"B");
		PureDiLepton1DIso[j][1]->Divide(PureDiLepton1DIso[j][1],PureDiLepton1DIso[j][0],1.,1.,"B");
		
	}
	
	
	
	
	//DM[1]->Draw("COLZ");
	//DM[1]->Draw("e text same");
	
	TFile* out = new TFile("PDL.root","RECREATE");
	out->cd();
	
	TCanvas *Out2D[5], *Out1D[7], *Out1DIso[7];
	
	TLine *ptgrid[4];
	
	ptgrid[0] = new TLine(10,0,10,50);ptgrid[0]->SetLineWidth(5);ptgrid[0]->SetLineColor(5);
	ptgrid[1] = new TLine(0,10,50,10);ptgrid[1]->SetLineWidth(5);ptgrid[1]->SetLineColor(5);
	ptgrid[2] = new TLine(25,0,25,50);ptgrid[2]->SetLineWidth(5);ptgrid[2]->SetLineColor(5);
	ptgrid[3] = new TLine(0,25,50,25);ptgrid[3]->SetLineWidth(5);ptgrid[3]->SetLineColor(5);

	/*for(int j=0;j<5;j++){
		
		
		Out2D[j] = new TCanvas(PureTrigs[j],PureTrigs[j],600,500);
		Out2D[j]->cd();
		
		
		PureDiLepton[j][1]->GetYaxis()->SetTitleOffset(0.65);
		PureDiLepton[j][1]->GetZaxis()->SetTitleSize(0.07);
		PureDiLepton[j][1]->GetZaxis()->SetTitleOffset(0.75);
		PureDiLepton[j][1]->Draw("COLZ");
		PureDiLepton[j][1]->Draw("e text same");
		
		for(int g=0;g<4;g++)
		  ptgrid[g]->Draw("same");
		  
		Out2D[j]->SaveAs("SavePDL6/"+PureTrigsFnames[j]+"_2D.pdf");
		
	}*/
	
	
	TF1* fiFunf = new TF1("effs","0.5*[0]*(TMath::Erf((x-[1])/[2])+1)",0,1000);
	fiFunf->SetParNames("#varepsilon_{#infty}", "x_{1/2} (GeV)", "#sigma (GeV)");
	TLine *lines[7];
	TLatex* threshold[7];
	
	double x12[7] = {8.0,8.0,14.0,23.0,9.0,14.0,25.0};
	double wid[7] = {1.5,1.0,15.0,6.0,3.0,15.0,15.0};
	
	double x12i[7] = {0.6,0.6,0.5,0.7,0.7,0.3,0.3};
	double widi[7] = {1,1,1.0,1,1,0.4,0.4};
	
	for(int j=0;j<7;j++){
	
		if(j != 1) continue;
	
		fiFunf->SetParameters(0.5, x12[j], wid[j]);
		
		Out1D[j] = new TCanvas(Form("%s_1D",PureTrigs1d[j]),Form("%s_1D",PureTrigs1d[j]),2);
		Out1D[j]->cd();
		
		PureDiLepton1D[j][1]->GetYaxis()->SetRangeUser(0,1.1);
		PureDiLepton1D[j][1]->GetYaxis()->SetTitleOffset(0.60);
		PureDiLepton1D[j][1]->Fit("effs");
		
		
		
		Out1D[j]->Update();
		//PureDiLepton1D[j][1]->Draw();
		
		TPaveStats *ps = (TPaveStats*)Out1D[j]->GetPrimitive("stats");
        ps->SetX1NDC(.5);
        ps->SetX2NDC(.9);
        ps->SetY1NDC(.25);
        ps->SetY2NDC(.55);
        gPad->Modified();
		
		double plat = fiFunf->GetParameter(1)+1.5*fiFunf->GetParameter(2);
        lines[j] = new TLine(plat,0,plat,1);
        lines[j]->SetLineColor(kGreen);
        lines[j]->SetLineWidth(3);
        lines[j]->Draw("same");
		
		TString tst = Form("%3.1f GeV", plat);
        threshold[j] = new TLatex(plat+10, 1.0,tst);
        threshold[j]->Draw("same");
		
		//Out1D[j]->SaveAs("SavePDL6/"+PureTrigs1dFnames[j]+"_1D.pdf");
		
		
		
		
		/*fiFunf->SetParameters(0.5,x12i[j], widi[j]);
		
		Out1DIso[j] = new TCanvas(Form("%s_1D_Iso",PureTrigs1d[j]),Form("%s_1D_Iso",PureTrigs1d[j]),2);
		Out1DIso[j]->cd();
		
		PureDiLepton1DIso[j][1]->GetYaxis()->SetRangeUser(0,1.1);
		PureDiLepton1DIso[j][1]->GetYaxis()->SetTitleOffset(0.60);
		PureDiLepton1DIso[j][1]->Fit("effs");
		
		
		
		Out1DIso[j]->Update();
		
		TPaveStats *ps = (TPaveStats*)Out1DIso[j]->GetPrimitive("stats");
        ps->SetX1NDC(.5);
        ps->SetX2NDC(.9);
        ps->SetY1NDC(.55);
        ps->SetY2NDC(.85);
        gPad->Modified();*/
		
		/*Out1DIso[j] = new TCanvas(Form("%s_1D_Iso",PureTrigs1d[j]),Form("%s_1D_Iso",PureTrigs1d[j]),2);
		Out1DIso[j]->cd();
		
		PureDiLepton1DIso[j][1]->GetYaxis()->SetRangeUser(0,1.0);
		PureDiLepton1DIso[j][1]->GetYaxis()->SetTitleOffset(0.8);
		PureDiLepton1DIso[j][1]->Draw();
		
		Out1DIso[j]->SaveAs("SavePDL/"+PureTrigs1dFnames[j]+"_Iso_1D.pdf");*/
		
	}


}
