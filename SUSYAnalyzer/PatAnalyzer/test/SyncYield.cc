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


//int Reduce(int origin) {
    /*    W_L,  // 0
     W_T_L, // 1
     W_B_L, // 2
     W_B_D_L, //3
     W_B_D_T_L, // 4
     W_B_T_L, // 5
     W_D_L, // 6
     W_D_T_L, //7
     B_L, // 8
     B_D_L, //9
     B_D_T_L, //10
     B_T_L,  // 11
     D_L, //12
     D_T_L, //13
     B_Baryon, // 14
     C_Baryon, //15
     pi_0, //16
     photon_, //17
     F_L, //18
     N_U_L_L // 19
     *//*
    int originR;
    if (origin < 2 )
        originR = 0;
    else if ( (origin == 2) || (origin == 13398 && _lumiBlock == 134) || (origin == 8) || (origin==11) || (origin==14) )
        originR = 1;
    else if ((origin == 3 ) || (origin == 4 ) || (origin == 6 ) || (origin==7) || (origin == 9 ) || (origin == 10 ) || (origin==12) || (origin==13) || (origin==15))
        originR = 2;
    else originR = origin - 13;
    
    if ((originR == 4) || (originR == 13398 && _lumiBlock == 134))
        originR = 3;
    else if (originR == 6) {
        originR = 4;
    }
    return originR;
    /*
     0    "Prompt",
     1    "b-jets",
     2    "c-jets",
     3    "uds",
     4    "Unknown"
     */
//}*/

double deltaPhi(TLorentzVector* lepton, double met, double metPhi) {
    double W_px = lepton->Px() + met*TMath::Cos(metPhi);
    double W_py = lepton->Py() + met*TMath::Sin(metPhi);
    double W_phi = TMath::ATan2(W_py,W_px);
    return(fabs(W_phi -lepton->Phi()));
}

void SyncYield()
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
    //double xSection = 102.5;
    //double xSection = 24.56;
    //double xSection = 37509;
    //double xSection = 6642;
	//double xSection = 424.5;//ttJets
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
	const char* firse[2] = {"First/","Second/"};
	const char* dirs[9] = {"TTWJetsNtuple/","TTZJetsNtuple/","WZJetsTo3LNuNtuple/","ZZto4LNtuple/","DYJetsM50ToLLNtuple/","TTJetsNtuple/","WJetsToLNuNtuple/","T1tttt12Ntuple/","T1tttt15Ntuple/"};
	int TwoDirs[9] = {2,2,1,2,2,2,2,2,1};
	int NFiles[9][2] = {{1,1},{1,1},{1,0},{2,2},{2,3},{19,19},{6,9},{1,1},{1,0}};
	int NumEvents[9] = {0};
	Double_t SampleScales[9] = {0};
	int WhichSample = 7;
	
	for(int sample=0;sample<9;sample++){
		
	if(sample != WhichSample) continue;
	
	int evs = 0;
    std::cout<<"\n\n";
	//const char* file = "SyncOutput/SSSync_Interactive_Full_PHYS14_25_V2_fixed3";//Sync File
	const char* file = "/scratch/osg/mrcarver/";//TTWJetsNtuple/Second/Job_1";
	
	for(int i=0;i<TwoDirs[sample];i++){
	for(int f=0;f<NFiles[sample][i];f++){
	
	std::stringstream ss;
	ss << file << dirs[sample] << firse[i] <<"Job_"<< f+1 << ".root";
	
	std::cout<<"File name is "<<ss.str().c_str()<<"\n";
	
	TFile *hfiles = new TFile(ss.str().c_str(),"READ");
	hfiles->cd("SyncExercise");
	TH1D* _hCounters = new TH1D("hCounter", "Events counter", 5,0,5);
	_hCounters->Read("hCounter");
	evs += _hCounters->GetBinContent(1);
	
	
	}
	}
	
	
	Double_t scale = lumi*SampleXSection[sample]/evs;
	SampleScales[sample] = scale;
	NumEvents[sample] = evs;
	std::cout<<"Scale = "<<scale<<" and total entries in all files = "<<evs<<"\n";
	
    
	
	
	for(int t=0;t<TwoDirs[sample];t++){
	for(int f=0;f<NFiles[sample][t];f++){
    
    std::stringstream ss;
	ss << file << dirs[sample] << firse[t] <<"Job_"<< f+1 << ".root";
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
    
    //std::cout<<"Entries "<<nEntries<<std::endl;
	//if(nEntries > 50)	
		//nEntries = 50;
		
	FILE *dump = fopen("output.txt","w");
	int SRYields[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	////My Plots/////
	//TH1F* SearchRegions = new TH1F("SearchRegions","Events per Search Region",29,0,29);SearchRegions->SetLineWidth(2);
	
	//TH1F* nSSL = new TH1F("nSSL","Number of SS Leptons",10,0,10);
	
	////////////////
	int SRs[4][3][3] = {{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}}};
	
	
    for (long it=0; it!=nEntries; ++it) {
    
		
        outputTree->GetEntry(it);
        
        if (it%10000 == 0)
            std::cout<<'.'<<std::flush;//std::cout<<((TLorentzVector*)_leptonP4->At(0))->Pt()<<std::endl;
		
		
		//int signalregion = SR(_n_bJets, _met,_n_Jets,HT,1);
		/////////////////////////////////////
		//// Select Best Pair of Leptons ////
		/////////////////////////////////////
		std::vector<int> li, hi, tight, loose, looseout, all;
		for(int i=0;i<_nLeptons;i++){
		
			
			if(_eventNb == 13398 && _lumiBlock == 134){
			  std::cout<<"\n\n"<<i<<",evt "<<_eventNb<<"::pT = "<<((TLorentzVector*)_leptonP4->At(i))->Pt()<<", eta = "<<((TLorentzVector*)_leptonP4->At(i))->Eta()<<", iso = "<<_isolation[i]<<" and SIP = "<<_3dIPsig[i]<<", pdgid = "<<pdgid(_flavors[i],_charges[i])<<
			  ", dxy = "<<_ipPV[i]<<", dz = "<<_ipZPV[i]<<", _dphi = "<<_dphi[i]<<", deta = "<<_deta[i]<<", SigIetaIeta = "<<_sigIeta[i]<<", HoE = "<<_HoE[i]<<" and pMe = "<<_pMe[i]<<",missing inner hits = "<<_missingHits[i]<<", isloose = "<<_isloose[i]<<", istight = "<<_istight[i]<<"\n";
			}
		
			//if(_isolation[i] < 0.6)
				all.push_back(i);
		
			if(  ((TLorentzVector*)_leptonP4->At(i))->Pt() <= 25 && ((TLorentzVector*)_leptonP4->At(i))->Pt() > 10  && _istight[i])
				li.push_back(i);
			else if(((TLorentzVector*)_leptonP4->At(i))->Pt() > 25 && _istight[i])
				hi.push_back(i);
				
			if(_isloose[i] && !_istight[i])
				loose.push_back(i);
			
			if(_isloose[i])
				looseout.push_back(i);
			
			if(_istight[i])
				tight.push_back(i);
				
			
			for(int j=1+1;j<_nLeptons;j++){

				TLorentzVector b1 = *((TLorentzVector*)_leptonP4->At(i)) + *((TLorentzVector*)_leptonP4->At(j));	
			
				if(_eventNb == 13398 && _lumiBlock == 134)
					std::cout<<"i = "<<i<<", j = "<<j<<" and Mll = "<<b1.M()<<"\n";

			}
		
		}
		//std::cout<<_eventNb<<"\n";;
		
		if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"\nloose.size = "<<looseout.size()<<"\nTight sizse = "<<tight.size()<<"\n\n";
		
		if(tight.size() < 2) continue;
		
		
		
		std::vector<int> SSL[2];
		for(unsigned int i=0;i<tight.size();i++){
			
				//if(_eventNb == 13398 && _lumiBlock == 1347){
					//std::cout<<"\npdg = "<<pdgid(_flavors[loose[i]],_charges[loose[i]])<<" and pt"<<i<<" = "<<((TLorentzVector*)_leptonP4->At(loose[i]))->Pt()<<"\n";
				//}
				if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"tc"<<i<<" = "<<_charges[tight[i]]<<" and pt = "<<((TLorentzVector*)_leptonP4->At(tight[i]))->Pt()<<"\n";
				
				if(_charges[tight[i]] < 0)
					SSL[0].push_back(tight[i]);
				else if(_charges[tight[i]] > 0)
					SSL[1].push_back(tight[i]);
					
				
		}
		
		if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"Two or more tight. SSL0 size = "<<SSL[0].size()<<", SSL1 size = "<<SSL[1].size()<<"\n";
		
		if(SSL[0].size() < 2 && SSL[1].size() < 2)
			continue;
			
		if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"two SS\n";
		
		std::vector<std::pair<int,int>> Pairs;
		
		//std::cout<<"here2\n";
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
		/////////////////////////////////
		//// Low Mass Resonance Veto ////
		/////////////////////////////////
		std::vector<int> GPairIndex, GPairType;
		for(unsigned int i=0;i<Pairs.size();i++){
		
			TLorentzVector both = *((TLorentzVector*)_leptonP4->At(Pairs[i].first)) + *((TLorentzVector*)_leptonP4->At(Pairs[i].second));
			if(both.M() > 8){
				GPairIndex.push_back(i);
				GPairType.push_back(_flavors[Pairs[i].first] + _flavors[Pairs[i].second]);
			}
		
		}
		
		if(!GPairIndex.size())
			continue;
		
		
			
		if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"low mass\n";
			
		if(GPairIndex.size() > 2){ //find out what to do when two of same type, only 1 event . Take highest pT pair
			//nSSL->Fill(0);
			//continue;
		}
		
		
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
		
		if(_eventNb == 13398 && _lumiBlock == 134){
			std::cout<<"Final pair indexs = "<<Pairs[GPairIndex[finalPair]].first<<" and "<<Pairs[GPairIndex[finalPair]].second<<"\n";
		}
		
		
		
		float htp[2] = {((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt(),((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt()};
		//int pdgi[2] = {pdgid(_flavors[Pairs[GPairIndex[finalPair]].first],_charges[Pairs[GPairIndex[finalPair]].first]),pdgid(_flavors[Pairs[GPairIndex[finalPair]].second],_charges[Pairs[GPairIndex[finalPair]].second])};
		int pdgi[2] = {_pdgids[Pairs[GPairIndex[finalPair]].first],_pdgids[Pairs[GPairIndex[finalPair]].second]};
		if(htp[1] > htp[0]){
			htp[0] = ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt();
			htp[1] = ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt();
			//pdgi[0] = pdgid(_flavors[Pairs[GPairIndex[finalPair]].second],_charges[Pairs[GPairIndex[finalPair]].second]);
			//pdgi[1] = pdgid(_flavors[Pairs[GPairIndex[finalPair]].first],_charges[Pairs[GPairIndex[finalPair]].first]);
			
			pdgi[0] = _pdgids[Pairs[GPairIndex[finalPair]].second];
			pdgi[1] = _pdgids[Pairs[GPairIndex[finalPair]].first];
		}
		
		
		/////////////////////////////
		//// Now Do Jet Cleaning ////
		/////////////////////////////
        
		float CleanedHT = 0;
		int nCJets = 0, nCBJets = 0;
		for(int i = 0 ; i < _n_Jets ;i++ ){
		
		
			TLorentzVector jt; jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],0);
			bool dr = false;
			//if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"\nNew Jet\n";
			for(unsigned int j=0;j<looseout.size();j++){
			
				//double dR = ((TLorentzVector *)_leptonP4->At(looseout[j]))->DeltaR( jt );
				if(((TLorentzVector *)_leptonP4->At(looseout[j]))->DeltaR( jt ) < 0.4 && ((TLorentzVector *)_leptonP4->At(looseout[j]))->Pt() > 10.0)
					dr = true;
					
				if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"J"<<i<<",L"<<j<<"::dR = "<<((TLorentzVector *)_leptonP4->At(looseout[j]))->DeltaR( jt )<<"\n";
			
			}
       
        	
        	//double dR1 = ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->DeltaR( jt );
        	//double dR2 = ((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->DeltaR( jt );
			
			bool bjet = false;
			if(_csv[i] > 0.814)
			   bjet = true;
			 
			
			if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"jet "<<i<<" pT = "<<jt.Pt()<<", eta = "<<jt.Eta()<<", phi = "<<jt.Phi()<<", bjet = "<<bjet<<"\n";//, dR = "<<dr<<"\n";
        
        	//if (dR1 < 0.4 || dR2 < 0.4) continue;
			if(dr) continue;
        
        	
        
        	CleanedHT += _jetPt[i];
        	nCJets++;
			if(_csv[i] > 0.814)
				nCBJets++;
			
    	}
		
		
		
		//std::cout<<"\n"<<Form("%1d %9d %12d\t%2d\t%+2d %5.1f\t%+2d %5.1f\t%d\t%2d\t%5.1f\t%6.1f",_runNb,_lumiBlock,_eventNb,loose.size(),pdgid(_flavors[Pairs[GPairIndex[finalPair]].first],_charges[Pairs[GPairIndex[finalPair]].first]),
			  //((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt(),pdgid(_flavors[Pairs[GPairIndex[finalPair]].second],_charges[Pairs[GPairIndex[finalPair]].second]),
			//((TLorentzVector *)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt(),
			//nCJets,nCBJets,_met,CleanedHT);
		
		///////////////////////
		//// Z/Gamma* Veto ////
		///////////////////////
		bool zveto = false, gveto = false, lmveto = false;
		
		for(unsigned int i=0;i<tight.size();i++){
			for(unsigned int j=0;j<looseout.size();j++){
			
			
				TLorentzVector b1 = *((TLorentzVector*)_leptonP4->At(looseout[j])) + *((TLorentzVector*)_leptonP4->At(tight[i]));	
				
				//std::cout<<"tight"<<i<<", loose"<<j<<"
				
				if(b1.M() > 76 && b1.M() < 106 && _flavors[looseout[j]] == _flavors[tight[i]] && _charges[looseout[j]]*_charges[tight[i]] < 0){
					zveto = true;
					if(_eventNb == 13398 && _lumiBlock == 134) std::cout<<"\nzveto:::tight"<<i<<" "<<((TLorentzVector*)_leptonP4->At(tight[i]))->Pt()<<" loose"<<j<<" "<<((TLorentzVector*)_leptonP4->At(looseout[j]))->Pt()<<", Mll = "<<b1.M()<<"\n";
				}
				
				if(b1.M() < 12 && _charges[looseout[j]]*_charges[tight[i]] < 0)
					gveto = true;
					
				if(_charges[looseout[j]]*_charges[tight[i]] > 0 && b1.M() < 8 && i != j)
					lmveto = true;
				
			}
			
			//for(int j=0;j<tight.size();j++){
			
				//TLorentzVector b1 = *((TLorentzVector*)_leptonP4->At(tight[j])) + *((TLorentzVector*)_leptonP4->At(tight[i]));	
			
				//if(_eventNb == 13398 && _lumiBlock == 134)
				//	std::cout<<"i = "<<i<<", j = "<<j<<" and Mll = "<<b1.M()<<"\n";
			
			//}
		}
		
		if(_eventNb == 13398 && _lumiBlock == 134)
			std::cout<<"before z veto\n";
		
		if(zveto) continue;
		if(_eventNb == 13398 && _lumiBlock == 134)
			std::cout<<"z veto\n";
		if(gveto) continue;
		if(_eventNb == 13398 && _lumiBlock == 134)
			std::cout<<"g veto\n";
			
		//if(lmveto)
		//	continue;
		if(_eventNb == 13398 && _lumiBlock == 134)
			std::cout<<"low mass veto\n";
		
		
		
				
				
		
		
		
		//if(!_istight[Pairs[GPairIndex[finalPair]].first]) continue;
		//if(!_istight[Pairs[GPairIndex[finalPair]].second]) continue;
		
		
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
					
			}
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
		
		
		
		int SignalRegion = SR(nCBJets,_met,nCJets,CleanedHT,1);///analysis = 0 for low pt and 1 for high pt
		
		
		int highlow = 0;//need to fix this if allow loose leptons < 10 GeV. Actualy not because tight requirement is still > 10 GeV
		if(((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].first))->Pt() > 25)
			highlow++;
			
		if(((TLorentzVector*)_leptonP4->At(Pairs[GPairIndex[finalPair]].second))->Pt() > 25)
			highlow++;
			
		if(highlow == 2 && SignalRegion >= 0){
		
			SRYields[SignalRegion]++;
			SignalRegionYields[sample][SignalRegion]++;
		
		}
		
		
		if(B0){
			SRs[0][highlow][GPairType[finalPair]]++;
			SR0SampleYields[sample][0][highlow][GPairType[finalPair]]++;
			fprintf(dump,Form("%1d %9d %12d\t%2d\t%+2d %5.1f\t%+2d %5.1f\t%d\t%2d\t%5.1f\t%6.1f \n",_runNb,_lumiBlock,_eventNb,looseout.size(),pdgi[0],
		 		 htp[0],pdgi[1],htp[1],nCJets,nCBJets,_met,CleanedHT));
		}
		
		if(SigReg){
			SRs[SigReg][highlow][GPairType[finalPair]]++;
			SR0SampleYields[sample][SigReg][highlow][GPairType[finalPair]]++;
		}
		
		
    }//event loop
    
	
	}
	}//loops over different dirs and ntpules per sample
	
	
	}//loop over samples
	
	//TFile* out = new TFile("syout.root","RECREATE");
	//out->cd();
	//nSSL->Write();
	//SearchRegions->Write();
	//out->Close();
	/*commented out when running over all the samples
	std::cout<<"\n\n\n";
	//float his[29] = {67.25,3.74, 0.91, 1.20, 2.01, 2.41, 1.57, 0.85, 2.72, 0, 31.12, 4.27, 1.12,1.98,3.80, 2.80, 2.03, 1.03,4.50,0,12.73, 1.62,0.50, 1.02,1.96,1.06,0.41, 0.41,2.55};
	
	float ECO[4][2][3] = {{{52,80,124},{168,196,278}},{{12,22,35},{51,65,80}},{{21,29,45},{52,63,104}},{{13,23,33},{50,51,73}}};
	float UCS[4][2][3] = {{{52,186,140},{178,418,295}},{{13,45,39},{56,130,77}},{{19,73,50},{59,133,117}},{{16,59,40},{48,121,78}}};
	*/
	std::cout<<"\n\n";
	std::string hls[3] = {"Low-Low","High- Low","High-High"};
	std::string nums[4] = {"0 ","10","20","30"};
	for(int i=0;i<4;i++){
		for(int j=2;j>0;j--){
			
			std::cout<<"SR "<<nums[i]<<"::::"<<hls[j]<<"::::Event Yield = "<<SR0SampleYields[WhichSample][i][j][0]*SampleScales[WhichSample]<<"    "<<SR0SampleYields[WhichSample][i][j][1]*SampleScales[WhichSample]<<"     "<<SR0SampleYields[WhichSample][i][j][2]*SampleScales[WhichSample]<<"\n";
			//std::cout<<"SR "<<nums[i]<<"::::"<<hls[j]<<"::::Event Yield = "<<SRs[i][j][0]<<"-> "<<SRs[i][j][0]/ECO[i][j-1][0]<<"    "<<SRs[i][j][1]<<"-> "<<SRs[i][j][1]/ECO[i][j-1][1]<<"     "<<SRs[i][j][2]<<"-> "<<SRs[i][j][2]/ECO[i][j-1][2]<<"\n";
			//std::cout<<"SR "<<nums[i]<<"::::"<<hls[j]<<"::::Event Yield = "<<SRs[i][j][0]<<"-> "<<SRs[i][j][0]/UCS[i][j-1][0]<<"    "<<SRs[i][j][1]<<"-> "<<SRs[i][j][1]/UCS[i][j-1][1]<<"     "<<SRs[i][j][2]<<"-> "<<SRs[i][j][2]/UCS[i][j-1][2]<<"\n";
		}
	}
	
	std::cout<<"\n\n";
	
	FILE *ychart = fopen("TexYields.txt","w");
	//fprintf(ychart,"\ hline \n\hline");
	
	for(int i=0;i<29;i++){
	
		//if(i == 9 || i == 19) continue;
		fprintf(ychart,Form("\n%1d",i));
		//std::cout<<"\n\n"<<i;
		for(int s=0;s<9;s++){
		
			fprintf(ychart,Form("\t& \t%5.3f pm \t%5.3f",SignalRegionYields[s][i]*SampleScales[s],sqrt(SignalRegionYields[s][i])*SampleScales[s]));
			//std::cout<<"  &  "<<SignalRegionYields[s][i]*SampleScales[s]<<" +/- "<<sqrt(SignalRegionYields[s][i])*SampleScales[s];//
		
		}
		
		//std::cout<<"SR"<<i<<":: Event Yield = "<<SR0SampleYields[0][i]*SampleScales[0]<<" +/- "<<sqrt(SR0SampleYields[0][i])*SampleScales[0]<<"\n";
		
		//std::cout<<"SR"<<i<<":: Event Yield = "<<SearchRegions->GetBinContent(i+1)*scale<<" +/- "<<sqrt(SearchRegions->GetBinContent(i+1))*scale<<" ::::  "<<SearchRegions->GetBinContent(i+1)*scale/(his[i])<<"\n";
		//std::cout<<"SR"<<i<<":: Event Yield = "<<SearchRegions->GetBinContent(i+1)<<" +/- "<<sqrt(SearchRegions->GetBinContent(i+1))<<"\n";
	
	}

}
