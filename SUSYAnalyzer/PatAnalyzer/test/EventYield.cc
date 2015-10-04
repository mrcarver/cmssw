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



int Reduce(int origin) {
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
     */
    int originR;
    if (origin < 2 )
        originR = 0;
    else if ( (origin == 2) || (origin == 5) || (origin == 8) || (origin==11) || (origin==14) )
        originR = 1;
    else if ((origin == 3 ) || (origin == 4 ) || (origin == 6 ) || (origin==7) || (origin == 9 ) || (origin == 10 ) || (origin==12) || (origin==13) || (origin==15))
        originR = 2;
    else originR = origin - 13;
    
    if ((originR == 4) || (originR == 5))
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
}

double deltaPhi(TLorentzVector* lepton, double met, double metPhi) {
    double W_px = lepton->Px() + met*TMath::Cos(metPhi);
    double W_py = lepton->Py() + met*TMath::Sin(metPhi);
    double W_phi = TMath::ATan2(W_py,W_px);
    return(fabs(W_phi -lepton->Phi()));
}

void EventYield()
{
    int leptonStudy = 1; //1=muon, 0-electron

    
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
    
    int _n_bJetsAll;
    int _n_JetsAll;
    
    double _jetEtaAll[100];
    double _jetPhiAll[100];
    double _jetPtAll[100];
    bool _bTaggedAll[100];
    double _csvAll[100];
    int _leptonIndex;
    int _closeIndex[4];
	
	int _nLeptons;
    
    int _eventType; //ee,mm,em
    bool _sb;
    bool _doubleF;
    int _index1 = -1;
    int _index2 = -1;

    
    int _indeces[4];
    int _flavors[4];
    double _charges[4];
    double _isolation[4];
    double _isolationComponents[4][4];
    double _isolationMC[4][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    
    int _origin[4];
    int _originReduced[4];
    double _PVchi2;
    double _PVerr[3];
    double _ipPV[4];
    double _ipPVerr[4];
    double _ipPVmc[4];
    
    double _ipZPV[4];
    double _ipZPVerr[4];
    
    double _3dIP[4];
    double _3dIPerr[4];
    double _3dIPsig[4];
    
    double _mt[4];
    
    double _closeJetPt[4];
    double _closeJetPtAll[4];
    
    double _closeJetAngAll[4];
    double _ptRelAll[4];
    double _closeJetPtAllMC[4];
    double _closeJetPtAllstatus[4];
    int _partonIdMatched[4];
    bool _sameParton[4];
    
    bool _isloose[4];
    bool _istight[4];
    
    int _n_PV;
    
    int _n_electrons;
    int _n_muons;
    int _n_taus;
    
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
    
   
    
  
    
    TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    double lumi = 10000;
    //double xSection = 102.5;
    //double xSection = 24.56;
    //double xSection = 37509;
    //double xSection = 6642;
	//double xSection = 424.5;//ttJets
	double xSection = 0.82705;//ttW
    
    
    
    //TFile *hfile=new TFile("ttwOutput/Job_1.root","READ");
    
   // hfile->cd("FakeElectrons");
   // _hCounter->Read("hCounter");
   // Double_t scale = lumi*xSection/(_hCounter->GetBinContent(1));
   // hfile->Close();
	
	int evs = 0;
    
	for(int i=1;i<22;i++){
		
		const char* file = "ttwOutput3/Job_";
		std::stringstream ss;
		ss << file << i << ".root";
		
		TFile *hfiles = new TFile(ss.str().c_str(),"READ");
		hfiles->cd("FakeElectrons");
		TH1D* _hCounters = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounters->Read("hCounter");
		evs += _hCounters->GetBinContent(1);
	}
	
	Double_t scale = lumi*xSection/evs;
	std::cout<<"Scale = "<<scale<<" and total entries in all files = "<<evs<<"\n";
	
    TClonesArray* _leptonP4 = new TClonesArray("TLorentzVector", 6);
    for (int i=0; i!=6; ++i) {
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
    
	
	////My Plots/////
	TH1F* SearchRegions = new TH1F("SearchRegions","Events per Search Region",29,0,29);SearchRegions->SetLineWidth(2);
	////////////////
	
    for(int files=1;files<22;files++){
	
	
	const char* file = "ttwOutput3/Job_";
	std::stringstream ss;
	ss << file << files << ".root";
    
    TChain *outputTree=new TChain("FakeElectrons/fakeTree");
	//outputTree->Add("ttwOutput/Job_1.root");
	outputTree->Add(ss.str().c_str());
	
	outputTree->SetBranchAddress("_leptonP4", &_leptonP4);
	
	outputTree->SetBranchAddress("_eventNb",   &_eventNb); 
    outputTree->SetBranchAddress("_runNb",     &_runNb);    
    outputTree->SetBranchAddress("_lumiBlock", &_lumiBlock);
    
    outputTree->SetBranchAddress("_nLeptons", &_nLeptons);
    outputTree->SetBranchAddress("_flavors", &_flavors);
    outputTree->SetBranchAddress("_charges", &_charges);
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

	
	
	
	
	
    for (long it=0; it!=nEntries; ++it) {
    
        outputTree->GetEntry(it);
        
        if (it%200 == 0)
            std::cout<<'.'<<std::flush;//std::cout<<((TLorentzVector*)_leptonP4->At(0))->Pt()<<std::endl;
			
			
		int signalregion = SR(_n_bJets, _met,_n_Jets,HT,1);
		
		if(_nLeptons != 2) continue;
		//if(_nLeptons == 2 && (_charges[0] != _charges[1])) continue;
        
		bool passLep = false;
		
		int l20 = 0;
		
		//std::vector<std::pair<int,int>> SelLep;
		
		for(int i=0;i<_nLeptons;i++){
		
			if(((TLorentzVector*)_leptonP4->At(i))->Pt() > 20)
				l20++;
				
			for(int j=i+1;j<_nLeptons;j++){
		
				if(_charges[i] != _charges[j]) continue;
				
				if(((TLorentzVector*)_leptonP4->At(i))->Pt() <= 20 || ((TLorentzVector*)_leptonP4->At(j))->Pt() <= 20) continue;
				
				if(((TLorentzVector*)_leptonP4->At(i))->Eta() >= 2.4 || ((TLorentzVector*)_leptonP4->At(j))->Eta() >= 2.4) continue;
				
				TLorentzVector both = *((TLorentzVector*)_leptonP4->At(i)) + *((TLorentzVector*)_leptonP4->At(j));
				
				if(both.M() <= 8) continue;
				
				if(_flavors[i] == 1 && _isolation[i] >= 0.1) continue;
				if(_flavors[j] == 1 && _isolation[j] >= 0.1) continue;
				if(_flavors[i] == 0 && _isolation[i] >= 0.09) continue;
				if(_flavors[j] == 0 && _isolation[j] >= 0.09) continue;
				
				if(_flavors[i] == 0 && (((TLorentzVector*)_leptonP4->At(i))->Eta() >= 1.442 && ((TLorentzVector*)_leptonP4->At(i))->Eta() <= 1.566 )) continue;
				if(_flavors[j] == 0 && (((TLorentzVector*)_leptonP4->At(j))->Eta() >= 1.442 && ((TLorentzVector*)_leptonP4->At(j))->Eta() <= 1.566 )) continue;
				
				passLep = true;
				//std::pair tp (i,j);
				//SelLep.push_back(tp);
			}
		}
		
		if(passLep && signalregion > 0)
			SearchRegions->Fill(signalregion);
			
		if(HT > 80 && passLep && _n_Jets >= 2){
		
			int bf = 0;
			
			if(_n_bJets == 1)
				bf = 1;
				
			if(_n_bJets >= 2)
				bf = 2;
			
			if(HT < 500 && _met > 30)
				SearchRegions->Fill(bf*10);
			else if(HT >= 500)
				SearchRegions->Fill(bf*10);
				
			
			
		}
        
    }
    
    //std::cout<<" processed all events"<<std::endl;
	}

	
	TFile* out = new TFile("hope_miss.root","RECREATE");
	out->cd();
	SearchRegions->Write();
	out->Close();
	
	std::cout<<"\n";
	float his[29] = {67.25,3.74, 0.91, 1.20, 2.01, 2.41, 1.57, 0.85, 2.72, 0, 31.12, 4.27, 1.12,1.98,3.80, 2.80, 2.03, 1.03,4.50,0,12.73, 1.62,0.50, 1.02,1.96,1.06,0.41, 0.41,2.55};
	
	
	for(int i=0;i<29;i++){
	
		if(i == 9 || i == 19) continue;
		
		std::cout<<"SR"<<i<<":: Event Yield = "<<SearchRegions->GetBinContent(i+1)*scale<<" +/- "<<sqrt(SearchRegions->GetBinContent(i+1))*scale<<" ::::  "<<SearchRegions->GetBinContent(i+1)*scale/(his[i])<<"\n";
	
	}

}
