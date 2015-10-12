#include <iostream>
#include <algorithm>
#include <vector>
#include <TStyle.h>



using namespace std;

void PromptContam()
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
	
	double lumi = 10000;//10 1/fb
	double lumi50ns = 40.2;
	double xSection = 0.82705;//ttW
	double ZZto4L = 15.4;
	double TTWJets = 0.82705;//
	double TTZJets = 0.62705;
	double WZJets = 19.83;
	double DYJets = 6025.2;
	double TTJets = 809.1;
	double WJetsToLNu = 61524.;//?
	double T14t12 = 0.0856418;
	double T14t15 = 0.0141903;

	TString nameFile = "/cms/data/store/user/t2/users/mrcarver/Run2Tuples/FRPromptContam/";

    TString files[3] = {"WJetsAll_Histos_FR","ZJetsAll_Histos_FR","TTJetsAll_Histos_FR"};
	double xsecs[3] = {WJetsToLNu,DYJets,TTJets};
	
	TString PNames[21] = {"HLT_Ele8_CaloIdM_TrackIdM_PFJet30","HLT_Ele12_CaloIdM_TrackIdM_PFJet30","HLT_Ele18_CaloIdM_TrackIdM_PFJet30","HLT_Ele23_CaloIdM_TrackIdM_PFJet30","HLT_Ele33_CaloIdM_TrackIdM_PFJet30",
						  "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30","HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30",
						  "HLT_Mu8_v","HLT_Mu17_v","HLT_Mu24_v","HLT_Mu34_v",
						  "HLT_Mu8_TrkIsoVVL_v","HLT_Mu17_TrkIsoVVL_v","HLT_Mu24_TrkIsoVVL_v","HLT_Mu34_TrkIsoVVL_v",
						  "EleBtag","EleBtag","MuBtag","MuBtag"};
   TCanvas *plots[21];
   TLegend *legend = new TLegend(0.75,0.7,0.92,0.89);
   legend->SetBorderSize(0);
   legend->SetFillColor(0);
   
   TString btag = "_Mt", btagm = "_Mt_Mod", reg = "MTDist";
   TString ptdist = "_PtDist", ptdistpre = "_PtDist_Prescaled";
   
   for(int p=3;p<4;p++){
   
   
	THStack* stack =  new THStack();
	
	TH1F *histos[3];
	
	for(int f=0;f<2;f++){
	
		TFile* file = new TFile(nameFile+files[f]+".root","read");
		file->cd("SyncExercise");
		
		double evs = 0;
		TH1D* _hCounters = new TH1D("hCounter", "Events counter", 5,0,5);
		_hCounters->Read("hCounter");
		evs += _hCounters->GetBinContent(1);
		
		double scale = lumi50ns*xsecs[f]/evs;
		
		
		TH1F *hist = new TH1F("hist","hist",5,0,5);
		histos[f] = new TH1F("histos","histos",5,0,5);
		
		/*if(p < 17){
			hist->Read(PNames[p]+reg);
			histos[f]->Read(PNames[p]+reg);
			//data->Read(PNames[p]+reg);
		}
		else if(p == 18 || p == 21){
			hist->Read(PNames[p]+btagm);
			histos[f]->Read(PNames[p]+btagm);
			//data->Read(PNames[p]+btagm);
		}
		else{
			hist->Read(PNames[p]+btag);
			histos[f]->Read(PNames[p]+btag);
			//data->Read(PNames[p]+btag);
		}*/
		
		
		hist->Read(PNames[p]+"MTDist");
		histos[f]->Read(PNames[p]+"MTDist");
		
		for(int b=0;b<5;b++){
		
			hist->SetBinContent(b+1,hist->GetBinContent(b+1)*scale);
			histos[f]->SetBinContent(b+1,histos[f]->GetBinContent(b+1)*scale);
		
		}
		
		hist->SetLineColor(f+2);
		hist->SetFillColor(f+2);
		histos[f]->SetLineColor(f+2);
		histos[f]->SetFillColor(f+2);
		
		stack->Add(hist,"hist");
		
	}
	
	double totalMC = histos[0]->Integral() + histos[1]->Integral();// + histos[2]->Integral();
	for(int f=0;f<3;f++){
		for(int b=0;b<5;b++){
		
			//histos[f]->SetBinContent(b+1,histos[f]->GetBinContent(b+1)/totalMC);
			
		}
		
		//stack->Add(histos[f]);
	}
	
	if(!p){
		legend->AddEntry(histos[0],"WJets");
		legend->AddEntry(histos[1],"DYJets");
		//legend->AddEntry(histos[2],"TTJets");
	}
	
	TFile* file = new TFile("/cms/data/store/user/t2/users/mrcarver/Run2Tuples/Data/ModBtag_Prescaled_FR_gjson.root","read");
	file->cd("Run2Ntuplizer");
	

	
	TH1F *data = new TH1F("data","data",5,0,5);
	
	
	
	
	if(p < 17)
		data->Read(PNames[p]+reg);
	else if(p == 18 || p == 21)
		data->Read(PNames[p]+btagm);
	else
		data->Read(PNames[p]+btag);
	
	data->Sumw2();
	data->SetMarkerColor(1);
	data->SetMarkerStyle(20);
	data->SetMarkerSize(1.5);
	data->SetLineColor(1);
	data->SetLineWidth(2);
	data->SetTitle("");
	
	if(!p)
		legend->AddEntry(data,"Data");
	
	//setting prescale
	//for(int b=0;b<5;b++)
	//	data->SetBinContent(b+1,data->GetBinContent(b+1)*80);
	
 
 	TString title = PNames[p]+"Mt";
	if(p == 18 || p == 20)
		title = PNames[p]+"Mt_Mod";
 
 	plots[p] = new TCanvas(title,title,2);
	plots[p]->cd();
	
	
	/*if(data->Integral() > totalMC){
		data->Draw("same");
		stack->Draw("same");
		data->Draw("same");
		legend->Draw("same");
	}
	else{
		stack->Draw();
		data->Draw("same");
		legend->Draw("same");
	}
	
	//plots[p]->SaveAs("PromptPlots_NN/"+PNames[p]+".png");
	*/
	
	
	
	
	
	//TCanvas *plot = new TCanvas("plot","plot",2);
	//plots[p]->cd();

    TPad *p1 = new TPad("p1","p1",0,0.3,1,1.0), *p2 = new TPad("p2","p2",0,0.05,1,0.3);

    p1->SetBottomMargin(0);
    p1->SetGridx();
    p1->SetGridy();
	p1->SetLogy();
    p1->Draw();
    p1->cd();
    data->SetStats(0);
	data->SetTitle("");
	data->GetYaxis()->SetLabelSize(0.07);
	
	int largestBin = 0, lowestBin = 99999999;
	for(int i=1;i<6;i++){
		if(data->GetBinContent(i) > largestBin)
			largestBin = data->GetBinContent(i);
			
		if(histos[0]->GetBinContent(i)+histos[1]->GetBinContent(i)/*+histos[2]->GetBinContent(i)*/ < lowestBin)
			lowestBin = histos[0]->GetBinContent(i)+histos[1]->GetBinContent(i);//+histos[2]->GetBinContent(i);
	}
	data->GetYaxis()->SetRangeUser(lowestBin*0.8,largestBin*1.5);
	
	
	
	if(data->Integral() > totalMC){
		data->Draw("same");
		stack->Draw("same");
		data->Draw("same");
		legend->Draw("same");
	}
	else{
		stack->Draw();
		data->Draw("same");
		legend->Draw("same");
	}
	
    //data->Draw();
	//stack->Draw("same");
	//data->Draw("same");
	
	
	
	float pbs[6] = {10.0,15.0,25.0,35.0,50.0,70.0};
	
	
	
	
    plots[p]->cd();
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.4);
    p2->SetGridx();
    p2->SetGridy();
    p2->Draw();
    p2->cd();
	
	
	//TH1F* h3 = new TH1F("h3","h3",5,pbs);
	TH1F* h3 = new TH1F("h3","h3",15,0,150);
	for(int b=0;b<5;b++){
	
		//std::cout<<"data"<<b<<" = "<<data->GetBinContent(b+1)<<"\n";
		//std::cout<<"stack"<<b<<" = "<<histos[0]->GetBinContent(b+1)+histos[1]->GetBinContent(b+1)+histos[2]->GetBinContent(b+1)<<"\n";
	
		h3->SetBinContent(b+1,histos[0]->GetBinContent(b+1)+histos[1]->GetBinContent(b+1));//+histos[2]->GetBinContent(b+1));
		
	}
    
    h3->SetLineColor(kBlack);
	h3->Divide(data);
	
	float hb = 0.0;
	for(int b=1;b<6;b++){
		if(h3->GetBinContent(b) > hb)
			hb = h3->GetBinContent(b);
	}
	
	h3->SetLineColor(kBlack);
    h3->SetMinimum(0.0);
    h3->SetMaximum(hb*1.3);
	h3->SetMarkerStyle(20);
    h3->SetMarkerSize(1.0);
    h3->SetTitle("");
    h3->GetYaxis()->SetTitle("MC/Data");
	h3->GetYaxis()->SetTitleSize(0.15);
	h3->GetYaxis()->SetTitleOffset(0.4);
	h3->GetYaxis()->SetLabelSize(0.18);
    h3->GetXaxis()->SetTitle("Leptonp_{T}");
    h3->GetXaxis()->SetTitleSize(0.2);
    h3->GetXaxis()->SetTitleSize(0.2);
    h3->GetXaxis()->SetLabelSize(0.2);
	h3->SetStats(0);
    h3->Draw("ep");
	
	//plots[p]->SaveAs("SeptemberFRPlots/"+title+"_PtDist.png");
	
	}	
     

}
