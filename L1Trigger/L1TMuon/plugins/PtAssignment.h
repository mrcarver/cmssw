////Function to calculte pT for a given InternalTrack
////
////2985826856 old checksum11 in DataFormats/L1TMuon/src/classes_def.xml
////1494215132 12

#include "EmulatorClasses.h"
#include "L1TMuonTextDumper.h"
#include "Forest.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void CalculatePt(L1TMuon::InternalTrack track){

	bool verbose = false;

	int dphi[6] = {-999,-999,-999,-999,-999,-999}, deta[6] = {-999,-999,-999,-999,-999,-999};
	int clct[4] = {-999,-999,-999,-999}, cscid[4] = {-999,-999,-999,-999};
	int phis[4] = {-999,-999,-999,-999}, etas[4] = {-999,-999,-999,-999}, mode = 0;;
	
	float theta_angle = ((track.theta)*0.2874016 + 8.5)*(3.14159265359/180);
	float eta = (-1)*log(tan(theta_angle/2));

	const TriggerPrimitiveStationMap stubs = track.getStubs();
		
	if(verbose) std::cout<<"Track eta = "<<eta<<" and has hits in stations ";//
	
	int x=0;
	for(unsigned int s=8;s<12;s++){
		if((stubs.find(s)->second).size() == 1){
			
			if(verbose) std::cout<<(stubs.find(s)->second)[0]->detId<CSCDetId>().station()<<" ";
			
			etas[s-8] = (fabs((stubs.find(s)->second)[0]->getCMSGlobalEta()) + 0.9)/(0.0125);
			phis[s-8] = track.phis[x];//(stubs.find(s)->second)[0]->getCMSGlobalPhi();//
			clct[s-8] = (stubs.find(s)->second)[0]->getPattern();
			cscid[s-8] = (stubs.find(s)->second)[0]->Id();
			
			switch(s-7){
				case 1: mode |= 1;break;
				case 2: mode |= 2;break;
				case 3: mode |= 4;break;
				case 4: mode |= 8;break;
				default: mode |= 0;
			}
			x++;
		}
	}
	
	if(verbose) std::cout<<"\nMode = "<<mode<<std::endl; 
	
    if(phis[0] > 0 && phis[1] > 0){
		dphi[0] = phis[1] - phis[0];
		deta[0] = etas[1] - etas[0];
	}
	if(phis[0] > 0 && phis[2] > 0){
		dphi[1] = phis[2] - phis[0];
		deta[1] = etas[2] - etas[0];
	}
	if(phis[0] > 0 && phis[3] > 0){
		dphi[2] = phis[3] - phis[0];
		deta[2] = etas[3] - etas[0];
	}
	if(phis[1] > 0 && phis[2] > 0){
		dphi[3] = phis[2] - phis[1];
		deta[3] = etas[2] - etas[1];
	}
	if(phis[1] > 0 && phis[3] > 0){
		dphi[4] = phis[3] - phis[1];
		deta[4] = etas[3] - etas[1];
	}
	if(phis[2] > 0 && phis[3] > 0){
		dphi[5] = phis[3] - phis[2];
		deta[5] = etas[3] - etas[2];
	}
	
	
	if(verbose){
		for(int f=0;f<4;f++){
			std::cout<<"\nphis["<<f<<"] = "<<phis[f]<<" and etas = "<<etas[f]<<std::endl;
			std::cout<<"\nclct["<<f<<"] = "<<clct[f]<<" and cscid = "<<cscid[f]<<std::endl;
		}
	
		for(int u=0;u<6;u++)
			std::cout<<"\ndphi["<<u<<"] = "<<dphi[u]<<" and deta = "<<deta[u]<<std::endl;
	}
	
	
	//if(mode == 7){
	
		//Forest* forest = new Forest();
		
		//forest->loadForestFromXML("ModeVariables/trees/7",64);
		
		//std::cout<<"forest loaded\n\n";
	
	//}
	
	
	
	
	
	
	
	
	//std::cout<<"\n\n";
}
