#include <iostream>
#include <algorithm>
#include <vector>
#include <TStyle.h>



using namespace std;

void MakePlots()
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

    TH1F* CompPlot[2][2][4][4];
	
	TString flavs[2] = {"Ele","Mu"};
	TString pts[2] = {"Low","High"};
	TString sbs[4] = {"00","10","20","30"};
	TString ori[4] = {"Bs","Cs","UDS","Ts"};

	TString nameFile = "SBCompOut_3sr";
	
	float SRValues[2][2][4][4][2] = {{{{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}},
									  {{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}}},
									 
									 
									 {{{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}},
									  {{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{0,0},{0,0}}}}};
    
    TFile* file = new TFile(nameFile+".root","read");
    TCanvas* plot[2][2][4];
 
    for(int i=0;i<2;++i) {
		//if(i) continue;
		for(int j=0;j<2;j++){
			//if(j) continue;
			for(int k=0;k<4;k++){
			//	if(k) continue;
		
        		TString tmp = "SBComp_"+flavs[i]+"_"+pts[j]+"pT_SB"+sbs[k];
        		plot[i][j][k] = new TCanvas(tmp,tmp,2);
				plot[i][j][k]->cd();
				
				TLegend *leg = new TLegend(0.4,0.4,0.5,0.6);
				leg->SetBorderSize(0);
				leg->SetFillColor(0);
				
				for(int l=0;l<3;l++){
					
					TString ts = "SBComp_"+flavs[i]+"_"+pts[j]+"pT_SB"+sbs[k]+"_"+ori[l];
				
        			CompPlot[i][j][k][l] = new TH1F(ts,ts,500, 0,1000);
        			CompPlot[i][j][k][l]->Read(ts);
					
					
					//TH1F *tmpPlot = new TH1F("tmpPlot","Fake Composition vs RelIso;RelIso;Event Percentage",21,0.0,2.1);
					
					//tmpPlot->SetTitle("");
					//tmpPlot->SetMarkerColor(l+1);
					//tmpPlot->SetLineColor(l+1);
        			//tmpPlot->GetYaxis()->SetRangeUser(0.001,1.1);
					
					//tmpPlot->SetBinContent(1,
					
					
					CompPlot[i][j][k][l]->SetTitle("");
					CompPlot[i][j][k][l]->SetMarkerColor(l+1);
					CompPlot[i][j][k][l]->SetLineColor(l+1);
        			CompPlot[i][j][k][l]->GetYaxis()->SetRangeUser(0.001,1.1);
					leg->AddEntry(CompPlot[i][j][k][l],ori[l]);
					
					if(!l)
					    CompPlot[i][j][k][l]->Draw();
					else
					    CompPlot[i][j][k][l]->Draw("same");
        	
				}
				
				//plot[i][j][k]->SetLogy();
				leg->Draw("same");
				
				
        		plot[i][j][k]->SaveAs("SBCompPlots/"+tmp+".png");
			}
        }
    }
    
    

    
}
