#include <iostream>
#include <algorithm>
#include <vector>
#include <TStyle.h>



using namespace std;

void DataFRPlots()
{

    int fontToUse = 132;
    gStyle->SetOptFit(1);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetPadColor(kWhite);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(1);
	gStyle->SetNdivisions(505,"XY");
    
    gStyle->SetAxisColor(1, "XYZ");
    gStyle->SetStripDecimals(kTRUE);
    gStyle->SetTickLength(0.03, "XYZ");
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);
	
	gStyle->SetLabelFont(fontToUse,"XYZ");
	gStyle->SetLabelSize(0.05,"XYZ");
	gStyle->SetLabelOffset(0.001,"X");
	
    gStyle->SetTitleFont(fontToUse,"");
	gStyle->SetTitleFont(fontToUse,"XYZ");
	gStyle->SetTitleFontSize(0.06);
	gStyle->SetTitleSize(0.06,"XY");
	//gStyle->SetTitleXOffset(1.0);
	gStyle->SetTitleXOffset(0.9);
	gStyle->SetTitleYOffset(1.25);
    
    gStyle->SetErrorX(0.);
    
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadLeftMargin(0.16);
    
    gStyle->SetStatFont(fontToUse);
    gStyle->SetStatColor(10);
    gStyle->SetStatFontSize(0.08);
    gStyle->SetTitleFont(fontToUse);
    gStyle->SetTitleFontSize(0.08);
	
	gStyle->SetMarkerSize(1.);
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerColor(gStyle->GetHistLineColor());
	
	gStyle->SetPalette(1,0);
	
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);
	gStyle->SetFuncColor(kRed);
	gStyle->SetPadRightMargin(0.07);
	gStyle->SetPadBottomMargin(0.16);

	TString flavs[2] = {"Ele","Mu"};
	TString pts[2] = {"Low","High"};
	TString sbs[4] = {"00","10","20","30"};
	
	float pbs[6] = {10.0,15.0,25.0,35.0,50.0,70.0};
	
	TString PNames[21] = {"HLT_Ele8_CaloIdM_TrackIdM_PFJet30_PtDist","HLT_Ele12_CaloIdM_TrackIdM_PFJet30_PtDist","HLT_Ele18_CaloIdM_TrackIdM_PFJet30_PtDist","HLT_Ele23_CaloIdM_TrackIdM_PFJet30_PtDist","HLT_Ele33_CaloIdM_TrackIdM_PFJet30_PtDist",
						  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_PtDist","HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_PtDist","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_PtDist","HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_PtDist",
						  "HLT_Mu8_v_PtDist","HLT_Mu17_v_PtDist","HLT_Mu24_v_PtDist","HLT_Mu34_v_PtDist",
						  "HLT_Mu8_TrkIsoVVL_v_PtDist","HLT_Mu17_TrkIsoVVL_v_PtDist","HLT_Mu24_TrkIsoVVL_v_PtDist","HLT_Mu34_TrkIsoVVL_v_PtDist",
						  "EleBtag_PtDist","EleBtag_PtDist_Mod","MuBtag_PtDist","MuBtag_PtDist_Mod"};

	TString nameFile = "ModBtag_Prescaled_FR";
	TFile* file = new TFile("/cms/data/store/user/t2/users/mrcarver/Run2Tuples/Data/"+nameFile+".root","read");
 	file->cd("Run2Ntuplizer");
    
	TH1F *TotalTight = new TH1F("TotalTight","TotalTight",5,pbs), *TotalFO = new TH1F("TotalFO","TotalFO",5,pbs);
	
	for(int f=17;f<19;f++){
		
		TH1F* FO = new TH1F("FO","FO",5,0,5), *Tight = new TH1F("Tight","Tight",5,0,5);
		FO->Read(PNames[f]);
		Tight->Read(PNames[f]+"_Tight");
	
		*TotalTight = *TotalTight + *Tight;
		*TotalFO = *TotalFO + *FO;
	
	}
	
	//TGraphAsymmErrors *FR = new TGraphAsymmErrors(TotalTight,TotalFO,"cp");
	//FR->SetLineColor(1);
	//FR->SetMarkerColor(1);
	//FR->SetLineWidth(2);
	//FR->Draw();

    

    
}
