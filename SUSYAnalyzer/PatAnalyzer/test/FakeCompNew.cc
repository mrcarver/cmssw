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

enum CType{BJets,CJets,Light,Tau};

double deltaPhi(TLorentzVector* lepton, double met, double metPhi) {
    double W_px = lepton->Px() + met*TMath::Cos(metPhi);
    double W_py = lepton->Py() + met*TMath::Sin(metPhi);
    double W_phi = TMath::ATan2(W_py,W_px);
    return(fabs(W_phi -lepton->Phi()));
}

void FakeCompNew()
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
    	double _isolation[8];
		double _miniIsolation[8];
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
		bool _istightNoIsoSIP[8];
		bool _islooseMVA[8];
		bool _istightMVANoIsoSIP[8];
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
    
  
    
   	// TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    double lumi = 10000;
	double xSection = 0.82705;//ttW
	double ZZto4L = 15.4;
	double TTWJets = 0.82705;//
	double TTZJets = 0.62705;
	double WZJets = 19.83;
	double DYJets = 6025.2;
	double TTJets = 809.1;
	double WJetsToLNu = 61524;//?
	double T14t12 = 0.0856418;
	double T14t15 = 0.0141903;
	double SampleXSection[9] = {TTWJets,TTZJets,WZJets,ZZto4L,DYJets,TTJets,WJetsToLNu,T14t12,T14t15};
    
    
	int SignalRegionYields[9][30] = {0};
	int SR0SampleYields[9][4][3][3] = {0};

	//int NumEvents[9] = {0};
	//Double_t SampleScales[9] = {0};
	//int WhichSample = 5;
	
	
	int evs = 0;
    	std::cout<<"\n\n";
	//const char* file = "/cms/data/store/user/t2/users/mrcarver/FRCompTuples/TTJets3/Zero/";
	const char* file = "/cms/data/store/user/t2/users/mrcarver/FRCompTuples/TTJetsMVA/Zero/";
	
	
	//const char* file = "/cms/data/store/user/t2/users/mrcarver/FRCompTuples/DYJetsHT/OneH/Zero/";
	//const char* file = "/cms/data/store/user/t2/users/mrcarver/FRCompTuples/WJetsHT/OneH/Zero/";
	
	
	for(int f=0;f<6;f++){
	
	std::stringstream ss;
	ss << file  <<"Job_"<< f+1 << ".root";
	
	std::cout<<"File name is "<<ss.str().c_str()<<"\n";
	
	TFile *hfiles = new TFile(ss.str().c_str(),"READ");
	hfiles->cd("SyncExercise");
	TH1D* _hCounters = new TH1D("hCounter", "Events counter", 5,0,5);
	_hCounters->Read("hCounter");
	evs += _hCounters->GetBinContent(1);
	
	
	}
	
	
	
	Double_t scale = lumi*SampleXSection[5]/evs;
	//SampleScales[5] = scale;
	//NumEvents[5] = evs;
	std::cout<<"Scale = "<<scale<<" and total entries in all files = "<<evs<<"\n";
	
	//[0] = b, [1] = c, [2] = light, [3] = tau, first index is on fake pT 10-25 or > 25
    	int SRTotal[2][2][4] = {{{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0}}};
	int SRComp[2][2][4][4] = {{{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}},{{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}}};
	
	
	float SBTotal[2][2][4][20], SBComp[2][2][4][20][4];
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<4;k++){
				for(int l=0;l<20;l++){
					SBTotal[i][j][k][l] = 0;
					for(int m=0;m<4;m++)
						SBComp[i][j][k][l][m] = 0;
				}
			}
		}
	}
	std::cout<<"here1\n";
	int NIso[2][2][2][11], NIsoJ[2][2][2][2][20], NIsoJAll[2][2][20];
	int NonIso[2][2][2][3][11], NonIsoJ[2][2][2][2][3][20], NonIsoJAll[2][2][3][20];
	
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<2;k++){
				for(int z=0;z<2;z++){
				for(int l=0;l<20;l++){
				
					NIsoJ[i][j][k][z][l] = 0;
					if(l < 11 && !z)
						NIso[i][j][k][l] = 0;
						
					if(!j && !z)
						NIsoJAll[i][k][l] = 0;
					
					for(int m=0;m<3;m++){
					
						NonIsoJ[i][j][k][m][z][l] = 0;
						if(!j && !z)
							NonIsoJAll[i][k][m][l] = 0;
						
						if(l < 11 && !z)
							NonIso[i][j][k][m][l] = 0;
					}
				}
				}
			}
		}
	}
	//std::cout<<"nisoja[0][0][0] = "<<NIsoJAll[0][0][0]<<"\n";
	
	
	TH1F* dRJ = new TH1F("dRJ","dRJ",50,0,5);
	
	std::cout<<"here2\n";
	for(int f=0;f<6;f++){
    
    	std::stringstream ss;
	ss << file <<"Job_"<< f+1 << ".root";
    	TChain *outputTree=new TChain("SyncExercise/SSSyncTree");
		outputTree->Add(ss.str().c_str());
	
		outputTree->SetBranchAddress("_leptonP4", &_leptonP4);
	
		outputTree->SetBranchAddress("_eventNb",   &_eventNb); 
    	outputTree->SetBranchAddress("_runNb",     &_runNb);    
    	outputTree->SetBranchAddress("_lumiBlock", &_lumiBlock);
    
    	outputTree->SetBranchAddress("_nLeptons", &_nLeptons);
    	outputTree->SetBranchAddress("_flavors", &_flavors);
    	outputTree->SetBranchAddress("_charges", &_charges);
		outputTree->SetBranchAddress("_pdgids", &_pdgids);
   		outputTree->SetBranchAddress("_isolation", &_isolation);
		outputTree->SetBranchAddress("_miniIsolation", &_miniIsolation);
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
	//outputTree->SetBranchAddress("_n_JetsAll", &_n_JetsAll);
    	outputTree->SetBranchAddress("_bTagged", &_bTagged);
    	outputTree->SetBranchAddress("_jetEta", &_jetEta);
    	outputTree->SetBranchAddress("_jetPhi", &_jetPhi);
    	outputTree->SetBranchAddress("_jetPt", &_jetPt);
	//outputTree->SetBranchAddress("_jetEtaAll", &_jetEtaAll);
    	//outputTree->SetBranchAddress("_jetPhiAll", &_jetPhiAll);
    	//outputTree->SetBranchAddress("_jetPtAll", &_jetPtAll);
    	outputTree->SetBranchAddress("_csv", &_csv);
    
    
    	long nEntries = outputTree->GetEntries();
    
    	//std::cout<<"Entries "<<nEntries<<std::endl;
		//if(nEntries > 10000)  
	    //   nEntries = 10000;
		
	FILE *dump = fopen("output.txt","w");
	int SRYields[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	
	////////////////
	int SRs[4][3][3] = {{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}}};
	
	bool verbose = false;
    	for (long it=0; it!=nEntries; ++it) {
    
		
        outputTree->GetEntry(it);
        
        if (it%10000 == 0)
            std::cout<<'.'<<std::flush;//std::cout<<((TLorentzVector*)_leptonP4->At(0))->Pt()<<std::endl;
		
		

		//////////////////////////////////////////////
		//// Select leptons by type(loose, tight) ////
		//////////////////////////////////////////////
		
		std::vector<int> tight,fakeSB,fakeNoIsoSIP, loose, looseout, all;
		for(int i=0;i<_nLeptons;i++){
		
			//if(_miniIsolation[i] < 0.6)
				all.push_back(i);
				
			if(_isloose[i] && !_istight[i])
				loose.push_back(i);
			
			if(_islooseMVA[i] && _miniIsolation[i] < 0.3 && _3dIPsig[i] < 6)//and isolation < 0.5
				looseout.push_back(i);
			
			if(_istightMVANoIsoSIP[i] && ((TLorentzVector*)_leptonP4->At(i))->Pt() > 10 && _miniIsolation[i] < 0.05 && _3dIPsig[i] < 4 && !_originReduced[i])
				tight.push_back(i);
				
			if(_istightNoIso[i] && _originReduced[i] > 0 && ((TLorentzVector*)_leptonP4->At(i))->Pt() > 10)
				fakeSB.push_back(i);
				
			if(_istightMVANoIsoSIP[i] && _originReduced[i] > 0 && ((TLorentzVector*)_leptonP4->At(i))->Pt() > 10)
				fakeNoIsoSIP.push_back(i);
		
		}
	
		if(!fakeSB.size() && !fakeNoIsoSIP.size()) continue;// or < 1 for SB selection
		
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
		
		if(verbose) std::cout<<"made SS Pairs\n";
		
		for(unsigned int i=0;i<SSL[0].size();i++){
		  for(unsigned int j=i+1;j<SSL[0].size();j++){
		      
			  std::pair<int,int> tmp;
			  tmp.first = SSL[0][i];
			  tmp.second = SSL[0][j];
			  Pairs.push_back(tmp);
			  
		  }
		}
		
		for(unsigned int i=0;i<SSL[1].size();i++){
		  for(unsigned int j=i+1;j<SSL[1].size();j++){
		      
			  std::pair<int,int> tmp;
			  tmp.first = SSL[1][i];
			  tmp.second = SSL[1][j];
			  Pairs.push_back(tmp);
			  
		  }
		}
		
		if(verbose) std::cout<<"before fakeNoIsoSIP\n";
		
		for(unsigned int i=0;i<tight.size();i++){
		  for(unsigned int j=0;j<fakeNoIsoSIP.size();j++){
		  	std::pair<int,int> tmp;
			if(_charges[tight[i]] == _charges[fakeNoIsoSIP[j]]){
			  tmp.first = tight[i];
		  	  tmp.second = fakeNoIsoSIP[j];
			  SBPairs.push_back(tmp);
			 if(verbose) std::cout<<"makes SB pair\n";
			}
		  }
		}
		
		if(verbose) std::cout<<"after fakeNoIsoSIP\n";
		
		if(SSL[0].size() < 2 && SSL[1].size() < 2 && !SBPairs.size())//&& SBPairs.size() < 1 for SB selection
		      continue;
		
		/////////////////////////////////
		//// Low Mass Resonance Veto ////
		/////////////////////////////////
		std::vector<int> GPairIndex, GPairType, GPairIndexSB, GPairTypeSB;
		for(unsigned int i=0;i<Pairs.size();i++){
		
			TLorentzVector both = *((TLorentzVector*)_leptonP4->At(Pairs[i].first)) + *((TLorentzVector*)_leptonP4->At(Pairs[i].second));
			if(both.M() > 8){
				GPairIndex.push_back(i);
				GPairType.push_back(_flavors[Pairs[i].first] + _flavors[Pairs[i].second]);
			}
		
		}
		
		if(verbose) std::cout<<"LMV, SBPairs.size() = "<<SBPairs.size()<<"\n";
		
		for(unsigned int i=0;i<SBPairs.size();i++){
		
			if(verbose) std::cout<<"inside LMV loop\n";
			TLorentzVector both = *((TLorentzVector*)_leptonP4->At(SBPairs[i].first)) + *((TLorentzVector*)_leptonP4->At(SBPairs[i].second));
			if(verbose) std::cout<<"inside LMV loop2\n";
			if(both.M() > 8){
				if(verbose) std::cout<<"> 8 cond\n";
				GPairIndexSB.push_back(i);
				GPairTypeSB.push_back(_flavors[SBPairs[i].first] + _flavors[SBPairs[i].second]);
				
				if(verbose) std::cout<<"One Good SB Pair\n";
				
			}
		
		}
		
		if(!GPairIndex.size() && !GPairIndexSB.size())
		      continue;
			
		if(verbose) std::cout<<"At Least One Good Pair\n";
		
		/////////////////////////////////////////////
		//// Select best pair by type and pT sum ////
		/////////////////////////////////////////////
		
		int finalPair = 0, finalPairSB = 0;
		float hpairpt = 0, hpairtype = -1, hpairptSB = 0, hpairtypeSB = -1;
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
		
		
		for(unsigned int i=0;i<GPairIndexSB.size();i++){
		
			float tsumpt = ((TLorentzVector*)_leptonP4->At(SBPairs[GPairIndexSB[i]].first))->Pt() + ((TLorentzVector*)_leptonP4->At(SBPairs[GPairIndexSB[i]].second))->Pt();
		
			if(GPairTypeSB[i] > hpairtypeSB){
				finalPairSB = i;
				hpairptSB = tsumpt;
				hpairtypeSB = GPairTypeSB[i];
			}
			else if((GPairTypeSB[i] == hpairtypeSB) && (tsumpt > hpairptSB)){
					finalPairSB = i;
					hpairptSB = tsumpt;
					hpairtypeSB = GPairTypeSB[i];
			}
			
		}
		
		/////////////////////////////
		//// Now Do Jet Cleaning ////
		/////////////////////////////
        
		float CleanedHT = 0;
		int nCJets = 0, nCBJets = 0;
		for(int i = 0 ; i < _n_Jets ;i++ ){
		
		
			TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],0);
			bool dr = false;
			for(unsigned int j=0;j<looseout.size();j++){
	
				if(((TLorentzVector *)_leptonP4->At(looseout[j]))->DeltaR( jt ) < 0.4 && ((TLorentzVector *)_leptonP4->At(looseout[j]))->Pt() > 10.0)
					dr = true;
			
			}
 
			bool bjet = false;
			if(_csv[i] > 0.814)
			   bjet = true;


			if(dr) continue;
 
        		CleanedHT += _jetPt[i];
        		nCJets++;
			if(_csv[i] > 0.814)
				nCBJets++;
			
    		}
		

		///////////////////////
		//// Z/Gamma* Veto ////
		///////////////////////
		bool zveto = false, gveto = false, lmveto = false;
		
		for(unsigned int i=0;i<tight.size();i++){
			for(unsigned int j=0;j<looseout.size();j++){
			
			
				TLorentzVector b1 = *((TLorentzVector*)_leptonP4->At(looseout[j])) + *((TLorentzVector*)_leptonP4->At(tight[i]));	
				
				if(b1.M() > 76 && b1.M() < 106 && _flavors[looseout[j]] == _flavors[tight[i]] && _charges[looseout[j]]*_charges[tight[i]] < 0)
					zveto = true;
				
				if(b1.M() < 12 && _charges[looseout[j]]*_charges[tight[i]] < 0)
					gveto = true;
					
				if(_charges[looseout[j]]*_charges[tight[i]] > 0 && b1.M() < 8 && i != j)
					lmveto = true;
				
			}
		}
		
		
		if(zveto) continue;
		
		if(verbose) std::cout<<"Z Veto\n";

		if(gveto) continue;
		
		if(verbose) std::cout<<"G Veto\n";
		
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
		
		if(B0 && verbose)
			std::cout<<"Baseline Region Reached\n";
		
		//int Analysis = 1;///analysis = 0 for low pt and 1 for high pt
		//int SignalRegion = SR(nCBJets,_met,nCJets,CleanedHT,Analysis);
		
		/////////////////////////////////
		//// First Do SR Composition ////
		/////////////////////////////////
		
		if(GPairIndex.size() && B0){
			
		/////////////////////////////////////////////
		//// Make sure one fake and one not fake ////
		/////////////////////////////////////////////

		int nfake = 0, forigin = -1, foriginR = -1, fpt = 0, flav = 0;
		
		
		
		if(_originReduced[Pairs[GPairIndex[finalPair]].first] > 0){
			nfake++;
			forigin = _origin[Pairs[GPairIndex[finalPair]].first];
			foriginR = _originReduced[Pairs[GPairIndex[finalPair]].first];
			if(((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt() > 25)
				fpt = 1;
				
			if(fabs(_pdgids[Pairs[GPairIndex[finalPair]].first]) == 13)
				flav = 1;
		}
			
		if(_originReduced[Pairs[GPairIndex[finalPair]].second] > 0){
			nfake++;
			forigin = _origin[Pairs[GPairIndex[finalPair]].second];
			foriginR = _originReduced[Pairs[GPairIndex[finalPair]].second];
			if(((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt() > 25)
				fpt = 1;
				
			if(fabs(_pdgids[Pairs[GPairIndex[finalPair]].second]) == 13)
				flav = 1;
		}
		
		if(verbose) std::cout<<"Before Fake Cut\n";

		if(nfake != 1) continue;
			
		SRTotal[flav][fpt][SigReg]++;
		
		//if(forigin == 4 || forigin == 5 || forigin == 7 || forigin == 10 || forigin == 11 || forigin == 13)
		//	SRComp[flav][fpt][SigReg][CType::Tau]++;
		if(foriginR == 1)
			SRComp[flav][fpt][SigReg][CType::BJets]++;
		else if(foriginR == 2)
			SRComp[flav][fpt][SigReg][CType::CJets]++;
		else if(foriginR == 3)
			SRComp[flav][fpt][SigReg][CType::Light]++;
		
		
		}
		
		///////////////////////////////
		//// Now do SB Composition ////
		///////////////////////////////
		
		if(GPairIndexSB.size() && B0){
		
		//std::cout<<"inside SB Comp\n";
		
		int HTbin = 0;
		if(CleanedHT > 400)
			HTbin = 1;
		

		int forigin = _origin[SBPairs[GPairIndexSB[finalPairSB]].second], foriginR = _originReduced[SBPairs[GPairIndexSB[finalPairSB]].second], fpt = 0, flav = 0;
		if(((TLorentzVector *)_leptonP4->At(SBPairs[GPairIndexSB[finalPairSB]].second))->Pt() > 25)
			fpt = 1;
		
		if(fabs(_pdgids[SBPairs[GPairIndexSB[finalPairSB]].second]) == 13)
			flav = 1;
		
		int iso = 0;
		for(int i=1;i<20;i++){
		
			if(_miniIsolation[SBPairs[GPairIndexSB[finalPairSB]].second] > 0.1*(i+1) && _miniIsolation[SBPairs[GPairIndexSB[finalPairSB]].second] <= 0.1*(i+2))
				iso = i;
			
			if(i == 19 && _miniIsolation[SBPairs[GPairIndexSB[finalPairSB]].second] > 2.0)
				iso = 19;
		
		}
			
		SBTotal[flav][fpt][SigReg][iso]++;
		
		//if(forigin == 4 || forigin == 5 || forigin == 7 || forigin == 10 || forigin == 11 || forigin == 13)
		//	SBComp[flav][fpt][SigReg][iso][CType::Tau]++;
		if(foriginR == 1)
			SBComp[flav][fpt][SigReg][iso][CType::BJets]++;
		else if(foriginR == 2)
			SBComp[flav][fpt][SigReg][iso][CType::CJets]++;
		else if(foriginR == 3)
			SBComp[flav][fpt][SigReg][iso][CType::Light]++;
			
		int ptb = 0, sip = 0;
		for(int i=3;i<13;i++){
			if(((TLorentzVector*)_leptonP4->At(SBPairs[GPairIndexSB[finalPairSB]].second))->Pt() > 5*i && ((TLorentzVector*)_leptonP4->At(SBPairs[GPairIndexSB[finalPairSB]].second))->Pt() <= 5*(i+1))
				ptb = i-2;
				
			if(((TLorentzVector*)_leptonP4->At(SBPairs[GPairIndexSB[finalPairSB]].second))->Pt() > 65)
				ptb = 10;
		}
		
		float cjpT = -1, cjdR = 999;
		int jptb = 0;
		/*for(int i = 0 ; i < _n_JetsAll ;i++ ){
		
		
			TLorentzVector jt; jt.SetPtEtaPhiM(_jetPtAll[i],_jetEtaAll[i],_jetPhiAll[i],0);
			if(  ((TLorentzVector *)_leptonP4->At(SBPairs[GPairIndexSB[finalPairSB]].second))->DeltaR( jt ) < cjdR && jt.Pt() > 10){
				cjdR = ((TLorentzVector *)_leptonP4->At(SBPairs[GPairIndexSB[finalPairSB]].second))->DeltaR( jt );
				cjpT = jt.Pt();
			}
			
		}*/
		
		dRJ->Fill(_closeJetAngAll[SBPairs[GPairIndexSB[finalPairSB]].second]);
		
		for(int i=1;i<21;i++){
			if(_closeJetPtAll[SBPairs[GPairIndexSB[finalPairSB]].second] > 5*i && _closeJetPtAll[SBPairs[GPairIndexSB[finalPairSB]].second] <= 5*(i+1))
				jptb = i-1;
				
			if(_closeJetPtAll[SBPairs[GPairIndexSB[finalPairSB]].second] > 105)
				jptb = 19;
		}
			
		if(_3dIPsig[SBPairs[GPairIndexSB[finalPairSB]].second] > 4)	
			sip = 1;
			
		
		
		if(_miniIsolation[SBPairs[GPairIndexSB[finalPairSB]].second] < 0.05){
			NIso[flav][HTbin][sip][ptb]++;
			if(jptb > -1) NIsoJ[flav][HTbin][sip][fpt][jptb]++;
			if(jptb > -1) NIsoJAll[flav][sip][jptb]++;
		}
		else if(_miniIsolation[SBPairs[GPairIndexSB[finalPairSB]].second] > 0.05){
		
		if(_miniIsolation[SBPairs[GPairIndexSB[finalPairSB]].second] < 0.5){
			NonIso[flav][HTbin][sip][0][ptb]++;
			if(jptb > -1) NonIsoJ[flav][HTbin][sip][fpt][0][jptb]++;
			if(jptb > -1) NonIsoJAll[flav][sip][0][jptb]++;
		}
		
		if(_miniIsolation[SBPairs[GPairIndexSB[finalPairSB]].second] < 1.0){
			NonIso[flav][HTbin][sip][1][ptb]++;
			if(jptb > -1) NonIsoJ[flav][HTbin][sip][fpt][1][jptb]++;
			if(jptb > -1) NonIsoJAll[flav][sip][1][jptb]++;
		}
		
		//if(_miniIsolation[SBPairs[GPairIndexSB[finalPairSB]].second] > 0.1){
			NonIso[flav][HTbin][sip][2][ptb]++;
			if(jptb > -1) NonIsoJ[flav][HTbin][sip][fpt][2][jptb]++;
			if(jptb > -1) NonIsoJAll[flav][sip][2][jptb]++;
		//}
		
		}
		
		}
		
		//for(unsigned int i=0;i<fakeSB.size();i++){
		/*if(fakeSB.size() == 1){
			
			
			
			int iso = 0, flav = 0, fpt = 0;
			for(int i=1;i<20;i++){
		
				if(_isolation[fakeSB[0]] > 0.1*(i+1) && _isolation[fakeSB[0]] <= 0.1*(i+2))
					iso = i;
			
				if(i == 19 && _isolation[fakeSB[0]] > 2.0)
					iso = 19;
		
			}
			
			if(fabs(_pdgids[fakeSB[0]]) == 13)
				flav = 1;
				
			if(((TLorentzVector *)_leptonP4->At(fakeSB[0]))->Pt() > 25)
				fpt = 1;
			
			
			if(_isolation[fakeSB[0]] > 0.1){
			
			
			
				SBTotal[flav][fpt][SigReg][iso]++;
		
				if(_originReduced[fakeSB[0]] == 1)
					SBComp[flav][fpt][SigReg][iso][CType::BJets]++;
				else if(_originReduced[fakeSB[0]] == 2)
					SBComp[flav][fpt][SigReg][iso][CType::CJets]++;
				else if(_originReduced[fakeSB[0]] == 3)
					SBComp[flav][fpt][SigReg][iso][CType::Light]++;
			
			}
			
			if(_isolation[fakeSB[0]] < 0.1){
			
				SRTotal[flav][fpt][SigReg]++;
		
				if(_originReduced[fakeSB[0]] == 1)
					SRComp[flav][fpt][SigReg][CType::BJets]++;
				else if(_originReduced[fakeSB[0]] == 2)
					SRComp[flav][fpt][SigReg][CType::CJets]++;
				else if(_originReduced[fakeSB[0]] == 3)
					SRComp[flav][fpt][SigReg][CType::Light]++;
			
			
			}
		
		
		}
		
		if(fakeNoIsoSIP.size() == 1){
			//std::cout<<"here\n";
		
			//if(_3dIPsig[fakeNoIsoSIP[0]] > 4){
			
			//std::cout<<"here2\n";
				
				int ptb = 0,flav = 0, sip = 0;
				for(int i=3;i<13;i++){
					if(((TLorentzVector*)_leptonP4->At(fakeNoIsoSIP[0]))->Pt() > 5*i && ((TLorentzVector*)_leptonP4->At(fakeNoIsoSIP[0]))->Pt() <= 5*(i+1))
						ptb = i-2;
						
					if(((TLorentzVector*)_leptonP4->At(fakeNoIsoSIP[0]))->Pt() > 65)
						ptb = 10;
				}
				
				if(fabs(_pdgids[fakeNoIsoSIP[0]]) == 13)
					flav = 1;
					
				if(_3dIPsig[fakeNoIsoSIP[0]] > 4)	
					sip = 1;
					
			
				
				if(_isolation[fakeNoIsoSIP[0]] < 0.1)
					NIso[flav][sip][ptb]++;
				else if(_isolation[fakeNoIsoSIP[0]] < 0.5)
					NonIso[flav][sip][0][ptb]++;
				else if(_isolation[fakeNoIsoSIP[0]] < 1.0)
					NonIso[flav][sip][1][ptb]++;
				else
					NonIso[flav][sip][2][ptb]++;
				
		
			//}
		}*/
		
		
    }//event loop
    
	
	}
	


	
	float FCompPerctent[4] = {0,0,0,0};
	std::string oris[4] = {"b-jets","c-jets","uds","taus"};
	std::string prange[2] = {"10 - 25","> 25"};
	std::string flavor[2] = {"Electrons","Muons"};
	
	char* orisC[4] = {"b-jets","c-jets","uds   ","taus  "};
	char* prangeC[2] = {"10 - 25","> 25"};
	char* flavorC[2] = {"Electrons","Muons"};
	
	FILE *SRoutput = fopen("SRFakeOrigin_DYHT_MVASelection_RegIso.txt","w");
	
	std::cout<<"\nSRComp[0][0][0][0] = "<<SRComp[0][0][0][0]<<", SRTotal[0][0][0] = "<<SRTotal[0][0][0]<<"\n\n";
	
	for(int i=0;i<2;i++){
		//std::cout<<"Fake "<<flavor[i]<<"\n\n";
		fprintf(SRoutput,"Fake %s\n\n",flavorC[i]);
		for(int j=0;j<2;j++){
			//std::cout<<"Fake pT "<<prange[j]<<" GeV\n\n";
			fprintf(SRoutput,"Fake pT %s GeV\n\n",prangeC[j]);
			for(int k=0;k<4;k++){
				for(int l=0;l<4;l++){
			
					if(SRTotal[i][j][k] == 0){
						//std::cout<<"SR"<<k<<"0 "<<oris[l]<<" % = 0 +/- 0\n";
						fprintf(SRoutput,"SR %d0 %s % = 0 +/- 0\n",k,orisC[l]);
					}
					else{
				
						float counts = SRComp[i][j][k][l], total = SRTotal[i][j][k];
				
						//std::cout<<"SR"<<k<<"0 "<<oris[l]<<" % = "<<counts/total<<" +/- "<<sqrt(counts)/total<<"\n";
						fprintf(SRoutput,"SR %d0 %s % = %1.4f +/- %1.4f\n",k,orisC[l],counts/total,sqrt(counts)/total);
					
					}
		
				}
				//std::cout<<"\n\n\n";
				fprintf(SRoutput,"\n\n\n");
			}
		}
	}


	TH1F *SBCompvsIso[2][2][4][4];
	const char *pr[2] = {"Low","High"};
	const char *SRs[4] = {"00","10","20","30"};
	const char *ori[4] = {"Bs ","Cs ","UDS","Ts "};
	const char *fla[2] = {"Ele","Mu"};
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<4;k++){
				for(int l=0;l<4;l++){
			
					SBCompvsIso[i][j][k][l] = new TH1F(Form("SBComp_%s_%spT_SB%s_%s",fla[i],pr[j],SRs[k],ori[l]),"Fake Composition vs RelIso;RelIso;Event Percentage",20,0.1,2.1);
			
				}
			}
		}
	}
	
	TH1F *FRs[2][4][2][3], *FRJs[2][2][2][2][3], *FRJAll[2][2][3];
	const char* fl[2] = {"E","M"};
	const char* si[2] = {"L","G"};
	const char* st[3] = {"0.5","1.0","Inf"};
	
	for(int i=0;i<2;i++){//flavor
		for(int j=0;j<2;j++){//SR or HT bin
			for(int k=0;k<2;k++){//L or G
				for(int z=0;z<2;z++){//fpt
				for(int l=0;l<3;l++){//Iso SB definition
				
				
					if(!z)
						FRs[i][j][k][l] = new TH1F(Form("%sFR%d_S%s4_B%d0",fl[i],l,si[k],j),Form("N_{iso}/N_{non-iso} SB < %s;Fake p_{T};Tight/Loose",st[l]),11,10,65);
					FRJs[i][j][k][z][l] = new TH1F(Form("%sFRJ%d_S%s4_B%d0_FpT%d",fl[i],l,si[k],j,z),Form("N_{iso}/N_{non-iso} SB < %s;Closest Jet p_{T};Tight/Loose",st[l]),20,5,105);
					if(!j && !z)
						FRJAll[i][k][l] = new TH1F(Form("%sFRJ%d_S%s4_BAll",fl[i],l,si[k]),Form("N_{iso}/N_{non-iso} SB < %s;Closest Jet p_{T};Tight/Loose",st[l]),20,5,105);
				}
				}
			}
		}
	}
	
	
	
	TFile* out = new TFile("SBCompOut_ttbar_MVASelection_3_MiniIso.root","RECREATE");
	out->cd();
	
	dRJ->Write();
	for(int z=0;z<2;z++){
	for(int k=0;k<2;k++){
	for(int l=0;l<2;l++){
	for(int j=0;j<2;j++){
	for(int i=0;i<20;i++){
	
		float fracj[3] = {NonIsoJ[j][l][k][z][0][i],NonIsoJ[j][l][k][z][1][i],NonIsoJ[j][l][k][z][2][i]};
		float nisoj = NIsoJ[j][l][k][z][i];
		
		//if(nisoj) std::cout<<"nisoj = "<<nisoj<<"\n";
		//if(fracj[0]) std::cout<<"fracj[0] = "<<fracj[0]<<"\n";
		
		float entj[3] = {nisoj/(fracj[0]+nisoj),nisoj/(fracj[1]+nisoj),nisoj/(fracj[2]+nisoj)};
		float entEj[3] = {sqrt(nisoj)/(fracj[0]+nisoj),sqrt(nisoj)/(fracj[1]+nisoj),sqrt(nisoj)/(fracj[2]+nisoj)};
		
		//if(entj[0] > 0) std::cout<<"entj[0] = "<<entj[0]<<"\n";
		//if(entj[1] > 0) std::cout<<"entj[1] = "<<entj[1]<<"\n";
		
		if(entj[0]> 0) FRJs[j][l][k][z][0]->SetBinContent(i+1,entj[0]);
		if(entj[0]> 0) FRJs[j][l][k][z][1]->SetBinContent(i+1,entj[1]);
		if(entj[0]> 0) FRJs[j][l][k][z][2]->SetBinContent(i+1,entj[2]);
		
		if(entEj[0]> 0) FRJs[j][l][k][z][0]->SetBinError(i+1,entEj[0]);
		if(entEj[0]> 0) FRJs[j][l][k][z][1]->SetBinError(i+1,entEj[1]);
		if(entEj[0]> 0) FRJs[j][l][k][z][2]->SetBinError(i+1,entEj[2]);
		
		
		if(!l && !z){
			
			float fracja[3] = {NonIsoJAll[j][k][0][i],NonIsoJAll[j][k][1][i],NonIsoJAll[j][k][2][i]};
			float nisoja = NIsoJAll[j][k][i];
			//std::cout<<"Nisoja = "<<NIsoJAll[j][k][i]<<"\n";
			
			//std::cout<<k<<"_"<<l<<"_"<<j<<"_"<<i<<"  nonisoja0 = "<<NonIsoJAll[j][k][0][i]<<", 1 = "<<NonIsoJAll[j][k][1][i]<<", 2 = "<<NonIsoJAll[j][k][2][i]<<"\n";
			float entja[3] = {0,0,0};
			//std::cout<<k<<"_"<<l<<"_"<<j<<"_"<<i<<"  nonisoja0 = "<<entja[0]<<", 1 = "<<entja[1]<<", 2 = "<<entja[2]<<"\n";
			float entEja[3] = {0,0,0};
			//std::cout<<k<<"_"<<l<<"_"<<j<<"_"<<i<<"  nonisoja0 = "<<entEja[0]<<", 1 = "<<entEja[1]<<", 2 = "<<entEja[2]<<"\n\n";
			
			for(int x=0;x<3;x++){
				
				if(fracja[x]+nisoja){
					entja[x] = nisoja/(fracja[x]+nisoja);
					entEja[x] = sqrt(nisoja)/(fracja[x]+nisoja);
				}
			}
	
			FRJAll[j][k][0]->SetBinContent(i+1,entja[0]);
			FRJAll[j][k][1]->SetBinContent(i+1,entja[1]);
			FRJAll[j][k][2]->SetBinContent(i+1,entja[2]);
		
			FRJAll[j][k][0]->SetBinError(i+1,entEja[0]);
			FRJAll[j][k][1]->SetBinError(i+1,entEja[1]);
			FRJAll[j][k][2]->SetBinError(i+1,entEja[2]);
	
		}
		
		if(i<11 && !z){
		float frac[3] = {NonIso[j][l][k][0][i],NonIso[j][l][k][1][i],NonIso[j][l][k][2][i]};
		float niso = NIso[j][l][k][i];
		float ent[3] = {niso/(frac[0]+niso),niso/(frac[1]+niso),niso/(frac[2]+niso)};
		float entE[3] = {sqrt(niso)/(frac[0]+niso),sqrt(niso)/(frac[1]+niso),sqrt(niso)/(frac[2]+niso)};
		FRs[j][l][k][0]->SetBinContent(i+1,ent[0]);
		FRs[j][l][k][1]->SetBinContent(i+1,ent[1]);
		FRs[j][l][k][2]->SetBinContent(i+1,ent[2]);
		
		FRs[j][l][k][0]->SetBinError(i+1,entE[0]);
		FRs[j][l][k][1]->SetBinError(i+1,entE[1]);
		FRs[j][l][k][2]->SetBinError(i+1,entE[2]);
		}
		
		
	
	}
	
		FRs[j][l][k][0]->SetMarkerStyle(20);
		FRs[j][l][k][1]->SetMarkerStyle(20);
		FRs[j][l][k][2]->SetMarkerStyle(20);
		FRs[j][l][k][0]->SetMarkerColor(1);
		FRs[j][l][k][1]->SetMarkerColor(2);
		FRs[j][l][k][2]->SetMarkerColor(3);
		FRs[j][l][k][0]->SetLineColor(1);
		FRs[j][l][k][1]->SetLineColor(2);
		FRs[j][l][k][2]->SetLineColor(3);
		FRs[j][l][k][0]->Write();
		FRs[j][l][k][1]->Write();
		FRs[j][l][k][2]->Write();
		
		FRJs[j][l][k][z][0]->SetMarkerStyle(20);
		FRJs[j][l][k][z][1]->SetMarkerStyle(20);
		FRJs[j][l][k][z][2]->SetMarkerStyle(20);
		FRJs[j][l][k][z][0]->SetMarkerColor(1);
		FRJs[j][l][k][z][1]->SetMarkerColor(2);
		FRJs[j][l][k][z][2]->SetMarkerColor(3);
		FRJs[j][l][k][z][0]->SetLineColor(1);
		FRJs[j][l][k][z][1]->SetLineColor(2);
		FRJs[j][l][k][z][2]->SetLineColor(3);
		FRJs[j][l][k][z][0]->Write();
		FRJs[j][l][k][z][1]->Write();
		FRJs[j][l][k][z][2]->Write();
		
		if(!l && !z){
			
			FRJAll[j][k][0]->SetMarkerStyle(20);
			FRJAll[j][k][1]->SetMarkerStyle(20);
			FRJAll[j][k][2]->SetMarkerStyle(20);
			FRJAll[j][k][0]->SetMarkerColor(1);
			FRJAll[j][k][1]->SetMarkerColor(2);
			FRJAll[j][k][2]->SetMarkerColor(3);
			FRJAll[j][k][0]->SetLineColor(1);
			FRJAll[j][k][1]->SetLineColor(2);
			FRJAll[j][k][2]->SetLineColor(3);
			FRJAll[j][k][0]->Write();
			FRJAll[j][k][1]->Write();
			FRJAll[j][k][2]->Write();
		
		}
	}
	}
	}
	}
	
	
	
	
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			for(int k=0;k<4;k++){
				for(int l=0;l<4;l++){
					for(int m=0;m<20;m++){
				
						float counts = SBComp[i][j][k][m][l], total = SBTotal[i][j][k][m];
					
						if(total == 0){
							SBCompvsIso[i][j][k][l]->SetBinContent(m+1,0);
						}
						else{
							SBCompvsIso[i][j][k][l]->SetBinContent(m+1,counts/total);
							SBCompvsIso[i][j][k][l]->SetBinError(m+1,sqrt(counts)/total);
						}
					
					
					}
					SBCompvsIso[i][j][k][l]->Write();
				}
			}
		}
	}















}
