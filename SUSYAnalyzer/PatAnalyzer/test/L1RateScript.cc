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

float DiEleScale(float pt1, float pt2){

	//return 0.5*0.9367*(TMath::Erf((pt - 9.463)/10.41) + 1);//TT_TuneZStar
	//return 0.5*0.9458*(TMath::Erf((pt - 8.962)/11.57) + 1);//TTJets
	//return 0.5*0.9446*(TMath::Erf((pt - 12.71)/10.9) + 1);//TTJets no iso selection

	float val1 = 0.5*0.9437*(TMath::Erf((pt1 - 11.65)/11.13) + 1);
	float val2 = 0.5*0.9437*(TMath::Erf((pt2 - 11.65)/11.13) + 1);
	float plat = 0.9437;
	
	return val1*val2/plat;

}

float MuEleScaleEle(float pt){

	//return 0.5*0.944*(TMath::Erf((pt - 9.376)/8.745) + 1);//TT_TuneZStar
	return 0.5*0.9198*(TMath::Erf((pt - 8.687)/10.65) + 1);//TTJets

}

float MuEleScale(float pte, float ptm, int flav){

	//return 0.5*0.944*(TMath::Erf((pt - 9.376)/8.745) + 1);//TT_TuneZStar
	//return 0.5*0.9198*(TMath::Erf((pt - 8.687)/10.65) + 1);//TTJets
	
	float vale = 0.5*0.9415*(TMath::Erf((pte - 11.75)/11.45) + 1);
	float valm = 0.5*0.9528*(TMath::Erf((ptm - 5.149)/5.514) + 1);
	float plat[2] = {0.9415,0.9528};
	
	return vale*valm/plat[flav];

}

float MuEleScaleMu(float pt){

	//return 0.5*0.9231*(TMath::Erf((pt - 10.15)/1.354) + 1);//TT_TuneZStar
	return 0.5*0.9221*(TMath::Erf((pt - 5.003)/4.847) + 1);//TTJets

}

float DiMuonScale(float pt1, float pt2){

	//return 0.5*0.9432*(TMath::Erf((pt - 10.04)/1.374) + 1);//TT_TuneZStar
	//return 0.5*0.9397*(TMath::Erf((pt - 4.999)/3.533) + 1);//TTJets
	
	float val1 = 0.5*0.94*(TMath::Erf((pt1 - 4.932)/1.935) + 1);
	float val2 = 0.5*0.94*(TMath::Erf((pt2 - 4.932)/1.935) + 1);
	float plat = 0.94;
	
	return val1*val2/plat;

}

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
	
	//if(met > 500)
	//	m = 2;
	if(met > 200)
		m = 1;
		
	//if(HT > 1600)
	//	ht = 2;
	if(HT > 300)
		ht = 1;
	
	
	if(nbjets > 2)
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

void L1RateScript()
{
    //int leptonStudy = 1; //1=muon, 0-electron

    
	gStyle->SetOptFit(1);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetPadColor(kWhite);
	gStyle->SetOptStat(1);
	gStyle->SetOptTitle(1);
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
    
 
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;
	
	
	
	double _L1Mu_Pt[40];
	double _L1Mu_Eta[40];
	double _L1Mu_Phi[40];
	int _nL1Mus;
	double _L1HT;
	
	TH1F *L1_HTT_Rate = new TH1F("L1_HTT_Rate",";L1_HTT;Rate [Hz]",100,0,1000);L1_HTT_Rate->SetLineColor(2);L1_HTT_Rate->SetMarkerColor(2);L1_HTT_Rate->SetLineWidth(2);
	TH1F *L1_Mu8_HTT_Rate = new TH1F("L1_Mu8_HTT_Rate",";L1_HTT;Rate [Hz]",100,0,1000);L1_Mu8_HTT_Rate->SetLineColor(3);L1_Mu8_HTT_Rate->SetMarkerColor(3);L1_Mu8_HTT_Rate->SetLineWidth(2);
	TH1F *L1_Mu_Rate = new TH1F("L1_Mu_Rate",";L1_Mu_Rate;Rate [Hz]",50,0,100);L1_Mu_Rate->SetLineColor(4);L1_Mu_Rate->SetMarkerColor(4);L1_Mu_Rate->SetLineWidth(2);

	double evs = 0.0;
    	std::cout<<"\n\n";

	
	TH1F *TotalEvents = new TH1F("TotalEvents","Total Events",5,0,5);
	TFile *hfiles = new TFile("file:L1ntuplizer.root","READ");//
	hfiles->cd("SyncExercise");
	TH1D* _hCounters = new TH1D("hCounter", "Events counter", 5,0,5);
	_hCounters->Read("hCounter");
	evs += _hCounters->GetBinContent(1);
	
	
	Double_t scale = (2508*11246.0)/(evs);
	
	std::cout<<"events = "<<evs<<" and Scale = "<<scale<<"\n";
	
	
	
	TotalEvents->SetBinContent(1,evs);
	
    TChain *outputTree=new TChain("SyncExercise/L1RateNtuplizerTree");
	outputTree->Add("file:L1ntuplizer.root");
	
	outputTree->SetBranchAddress("_nL1Mus", &_nL1Mus);
	outputTree->SetBranchAddress("_L1HT", &_L1HT);
	outputTree->SetBranchAddress("_L1Mu_Pt", &_L1Mu_Pt);
	outputTree->SetBranchAddress("_L1Mu_Eta", &_L1Mu_Eta);
	outputTree->SetBranchAddress("_L1Mu_Phi", &_L1Mu_Phi);
	
	outputTree->SetBranchAddress("_eventNb",   &_eventNb); 
    outputTree->SetBranchAddress("_runNb",     &_runNb);    
    outputTree->SetBranchAddress("_lumiBlock", &_lumiBlock);
	
    
	TLorentzVector jt;
	
    long nEntries = outputTree->GetEntries();
    
    //std::cout<<"Entries "<<nEntries<<std::endl;
    	//if(nEntries > 100000)  
    	//   nEntries = 100000;
		
	//FILE *dump = fopen("output.txt","w");
	
	
	////////////////
	int SRs[4][3][3] = {{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0}}};//[sigReg][lepton event type][lepton pt type]
	
	bool verbose = false;
    	for (long it=0; it!=nEntries; ++it) {
    
		
        outputTree->GetEntry(it);
		
        
        if (it%10000 == 0)
            std::cout<<'.'<<std::flush;//std::cout<<((TLorentzVector*)_leptonP4->At(0))->Pt()<<std::endl;
			
		bool Mu8 = false;
		bool Mu[50] = {false};
		for(int i=0;i<_nL1Mus;i++){
		
		
			if(_L1Mu_Pt[i] >= 8)
				Mu8 = true;
				
			for(int j=0;j<50;j++){
			
				Mu[j] = false;
				if(_L1Mu_Pt[i] >= j*2)
					Mu[j] = true;
			
			}
		
		}
			
		for(int i=0;i<100;i++){
		
			if(_L1HT >= i*10){
				L1_HTT_Rate->SetBinContent(i+1,(L1_HTT_Rate->GetBinContent(i+1) + 1));
				
				if(Mu8)
					L1_Mu8_HTT_Rate->SetBinContent(i+1,(L1_Mu8_HTT_Rate->GetBinContent(i+1) + 1));
				
			}
			
			if(i<50){
				if(Mu[i])
					L1_Mu_Rate->SetBinContent(i+1,(L1_Mu_Rate->GetBinContent(i+1) + 1));
		
			}
		}
		
		
		
    }//event loop
	
	for(int i=0;i<100;i++){
	
		L1_HTT_Rate->SetBinContent(i+1,(L1_HTT_Rate->GetBinContent(i+1))*scale);
		L1_Mu8_HTT_Rate->SetBinContent(i+1,(L1_Mu8_HTT_Rate->GetBinContent(i+1))*scale);
		L1_Mu_Rate->SetBinContent(i+1,(L1_Mu_Rate->GetBinContent(i+1))*scale);
		
		L1_HTT_Rate->SetBinError(i+1,sqrt(L1_HTT_Rate->GetBinContent(i+1))*scale);
		L1_Mu8_HTT_Rate->SetBinError(i+1,sqrt(L1_Mu8_HTT_Rate->GetBinContent(i+1))*scale);
		L1_Mu_Rate->SetBinError(i+1,sqrt(L1_Mu_Rate->GetBinContent(i+1))*scale);
	
	
	}
	

	TString fileOutName = "L1Rate_Output_new.root";
	TFile *out = new TFile(fileOutName,"RECREATE");
	
	out->cd();
	L1_HTT_Rate->Write();
	L1_Mu8_HTT_Rate->Write();
	L1_Mu_Rate->Write();


}
