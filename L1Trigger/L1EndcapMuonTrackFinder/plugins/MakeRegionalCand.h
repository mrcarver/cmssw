//////Takes in the values of track quantites and returns an 
//////L1TReginalMuonCandidate with proper values and legths
//////
//////
//////
//////


#include "DataFormats/L1TMuon/interface/L1TRegionalMuonCandidate.h"
#include "DataFormats/L1TMuon/interface/L1TRegionalMuonCandidateFwd.h"


l1t::L1TRegionalMuonCandidate MakeRegionalCand(float pt, int phi, int theta, 
											   int sign, int quality, 
											   int trackaddress, int link){
											   
	l1t::L1TRegionalMuonCandidate Cand;									   

	float theta_angle = (theta*0.2874016 + 8.5)*(3.14159265359/180);
	float eta = (-1)*log(tan(theta_angle/2));
	
	int iEta = 0;
	if(eta < 0) iEta & 256;
	

	Cand.setHwPt(pt);
	Cand.setHwEta(iEta);
  	Cand.setHwPhi(phi/4);//this is relative phi not global. Needs to be decided on still
  	Cand.setHwSign(1);
  	Cand.setHwQual(quality);
  	Cand.setHwTrackAddress(1);
  	Cand.setLink(1);

	return Cand;

}
