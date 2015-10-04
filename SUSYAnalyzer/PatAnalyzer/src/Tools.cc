#include "SUSYAnalyzer/PatAnalyzer/interface/Tools.h"


float tools::mass(float pt1 , float pt2, float eta1 , float eta2, float phi1, float phi2)
{
    TLorentzVector v1,v2;
    v1.SetPtEtaPhiM(pt1, eta1, phi1, 0.);
    v2.SetPtEtaPhiM(pt2,eta2, phi2, 0.);
    v1 = v1 + v2;
    return v1.Mag();
}

void tools::ERR( edm::InputTag& IT )
{
    cerr << "[ERROR] : " << IT << " is not a valid input label for this event.  SKIPPING EVENT " << endl;
}

//Muon pfRelIso
double tools::pfRelIso15(const pat::Muon *mu, double myRho)
{
    double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
    double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
    double photonIso = mu->pfIsolationR03().sumPhotonEt;
    double beta = mu->pfIsolationR03().sumPUPt;
	double Aeff[5] = { 0.0913, 0.0765, 0.0546, 0.0728, 0.1177 };
    double CorrectedTerm=0.0;
    
    if( TMath::Abs( mu->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( mu->eta() ) > 0.8 && TMath::Abs( mu->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( mu->eta() ) > 1.3 && TMath::Abs( mu->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( mu->eta() ) > 2.0 && TMath::Abs( mu->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
    else if( TMath::Abs( mu->eta() ) > 2.2 && TMath::Abs( mu->eta() ) < 2.5  )     CorrectedTerm = myRho * Aeff[ 4 ];
	
	
    double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - CorrectedTerm) )/mu->pt() ;
	
    return pfRelIsoMu;
}

double tools::pfRelIso(const pat::Muon *mu)
{
    double chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
    double neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
    double photonIso = mu->pfIsolationR03().sumPhotonEt;
    double beta = mu->pfIsolationR03().sumPUPt;
    double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
    return pfRelIsoMu;
}

//Electron pfRelIso
double tools::pfRelIso(const pat::Electron *iE, double myRho)
{
    //double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13  };
    //double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14  };
	double Aeff[5] = { 0.1013, 0.0988, 0.0572, 0.0842, 0.1530 };
    double CorrectedTerm=0.0;
    /*if( TMath::Abs( iE->eta() ) < 1.0                                            )   CorrectedTerm = myRho * Aeff[ 0 ];
     else if( TMath::Abs( iE->eta() ) > 1.0   && TMath::Abs( iE->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
     else if( TMath::Abs( iE->eta() ) > 1.479 && TMath::Abs( iE->eta() ) < 2.0    )   CorrectedTerm = myRho * Aeff[ 2 ];
     else if( TMath::Abs( iE->eta() ) > 2.0   && TMath::Abs( iE->eta() ) < 2.2    )   CorrectedTerm = myRho * Aeff[ 3 ];
     else if( TMath::Abs( iE->eta() ) > 2.2   && TMath::Abs( iE->eta() ) < 2.3    )   CorrectedTerm = myRho * Aeff[ 4 ];
     else if( TMath::Abs( iE->eta() ) > 2.3   && TMath::Abs( iE->eta() ) < 2.4    )   CorrectedTerm = myRho * Aeff[ 5 ];*/
    
    if( TMath::Abs( iE->superCluster()->eta() ) < 1.0 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 1.0 && TMath::Abs( iE->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 1.479 && TMath::Abs( iE->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 2.0 && TMath::Abs( iE->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 2.2 && TMath::Abs( iE->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 2.3 && TMath::Abs( iE->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
    else CorrectedTerm = myRho * Aeff[ 6 ];
    
    //double pfRelIsoE = (iE->chargedHadronIso() + TMath::Max(0.0, iE->neutralHadronIso() + iE->photonIso() - CorrectedTerm ) ) /iE->pt() ;
	double pfRelIsoE = (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - 0.5*iE->pfIsolationVariables().sumPUPt ) ) /iE->pt() ;

    return pfRelIsoE;
}

double tools::pfRelIso15(const pat::Electron *iE, double myRho)
{

	double Aeff[5] = { 0.1013, 0.0988, 0.0572, 0.0842, 0.1530 };
    double CorrectedTerm=0.0;
    
    if( TMath::Abs( iE->superCluster()->eta() ) < 0.8 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 0.8 && TMath::Abs( iE->superCluster()->eta() ) < 1.3  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 1.3 && TMath::Abs( iE->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 2.0 && TMath::Abs( iE->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
    else if( TMath::Abs( iE->superCluster()->eta() ) > 2.2 && TMath::Abs( iE->superCluster()->eta() ) < 2.5  )     CorrectedTerm = myRho * Aeff[ 4 ];
    
    double pfRelIsoE = (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - CorrectedTerm ) ) /iE->pt() ;
	//double pfRelIsoE = (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - 0.5*iE->pfIsolationVariables().sumPUPt ) ) /iE->pt() ;

    return pfRelIsoE;
}

double tools::pfRelIso(const pat::Electron *iE)
{
   
    //double pfRelIsoE = (iE->chargedHadronIso() + TMath::Max(0.0, iE->neutralHadronIso() + iE->photonIso() - CorrectedTerm ) ) /iE->pt() ;
	double pfRelIsoE = (iE->pfIsolationVariables().sumChargedHadronPt + TMath::Max(0.0, iE->pfIsolationVariables().sumNeutralHadronEt + iE->pfIsolationVariables().sumPhotonEt - 0.5*iE->pfIsolationVariables().sumPUPt ) ) /iE->pt() ;

    return pfRelIsoE;
}

double tools::getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                        const reco::Candidate* ptcl,  
                        double r_iso_min, double r_iso_max, double kt_scale,
                        bool use_pfweight, bool charged_only, double rho) {

    if (ptcl->pt()<5.) return 99999.;

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.); 
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;

      double dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;
      
      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          double wpf(1.);
          /*if (use_pfweight){
            double wpv(0.), wpu(0.);
            for (const pat::PackedCandidate &jpfc : *pfcands) {
              double jdr = deltaR(pfc, jpfc);
              if (pfc.charge()!=0 || jdr<0.00001) continue;
              double jpt = jpfc.pt();
              if (pfc.fromPV()>1) wpv *= jpt/jdr;
              else wpu *= jpt/jdr;
            }
            wpv = log(wpv);
            wpu = log(wpu);
            wpf = wpv/(wpv+wpu);
          }*/
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += wpf*pfc.pt();
          /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += wpf*pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    double iso(0.);
    //if (charged_only){
    //  iso = iso_ch;
    //} else {
    //  iso = iso_ph + iso_nh;
    //  if (!use_pfweight) iso -= 0.5*iso_pu;
    //  if (iso>0) iso += iso_ch;
    //  else iso = iso_ch;
    //}
	
	int em = 0;
	if(ptcl->isMuon())
		em = 1;
		
	double Aeff[2][5] = {{ 0.1013, 0.0988, 0.0572, 0.0842, 0.1530 },{ 0.0913, 0.0765, 0.0546, 0.0728, 0.1177 }};
	
	
	double CorrectedTerm=0.0;
	double riso2 = r_iso*r_iso;
    if( TMath::Abs( ptcl->eta() ) < 0.8 ) CorrectedTerm = rho * Aeff[em][ 0 ]*(riso2/0.09);
    else if( TMath::Abs( ptcl->eta() ) > 0.8 && TMath::Abs( ptcl->eta() ) < 1.3  )   CorrectedTerm = rho * Aeff[em][ 1 ]*(riso2/0.09);
    else if( TMath::Abs( ptcl->eta() ) > 1.3 && TMath::Abs( ptcl->eta() ) < 2.0  )   CorrectedTerm = rho * Aeff[em][ 2 ]*(riso2/0.09);
    else if( TMath::Abs( ptcl->eta() ) > 2.0 && TMath::Abs( ptcl->eta() ) < 2.2  )   CorrectedTerm = rho * Aeff[em][ 3 ]*(riso2/0.09);
    else if( TMath::Abs( ptcl->eta() ) > 2.2 && TMath::Abs( ptcl->eta() ) < 2.5  )   CorrectedTerm = rho * Aeff[em][ 4 ]*(riso2/0.09);
    
	//std::cout<<"riso = "<<r_iso<<", iso_nh = "<<iso_nh<<", iso_ch = "<<iso_ch<<", iso_ph = "<<iso_ph<<",rho = "<<rho<<" adn CTerm = "<<CorrectedTerm<<"\n";  
	 
	iso = iso_ch + TMath::Max(0.0, iso_ph + iso_nh - CorrectedTerm );
	
	
    iso = iso/ptcl->pt();

    return iso;
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Matt Additions //////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<const pat::Electron* > tools::MVALooseElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS, EGammaMvaEleEstimatorCSA14* myMVATrig){
	
	std::vector<const pat::Electron* > vElectrons;
	for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
	const reco::GsfTrackRef gsfTrack = el->gsfTrack();
    if (!gsfTrack.isNonnull()) {
        continue;
    }

	if( TMath::Abs(el->eta()) > 2.5 ) continue;
	if( el->pt() < 7 ) continue;
	if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;
	bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
    if( vtxFitConversion )  continue;
	//if( TMath::Abs(gsfTrack->dxy(PV)) > 0.05  )  continue;
   // if( TMath::Abs(gsfTrack->dz(PV)) > 0.1  ) continue;
	
	double mvaThresh = 1.0;
	if( TMath::Abs(el->eta()) < 0.8)
		mvaThresh = -0.32;//-0.11;
	else if(TMath::Abs(el->eta()) >= 0.8 && TMath::Abs(el->eta()) < 1.479)
		mvaThresh = -0.50;//-0.35;
	else
		mvaThresh = -0.70;//-0.55;
		
	double mvaVal = myMVATrig->mvaValue(*el,false);
	if(mvaVal < mvaThresh) continue;
	
	vElectrons.push_back(&*el );

	}
	
	return vElectrons;

}

/*bool tools::isMVATightNoIsoSIP(const pat::Electron  thePatElectron, EGammaMvaEleEstimatorCSA14* myMVATrig)
{//

        if( bool_electron_chargeConsistency && !thePatElectron.isGsfCtfScPixChargeConsistent() )  return false;
     
        double mvaVal = myMVATrig->mvaValue(thePatElectron,false);
		
		double mvaThresh = 1.0;
		if( TMath::Abs(thePatElectron.eta()) < 0.8)
			mvaThresh = 0.73;
		else if(TMath::Abs(thePatElectron.eta()) >= 0.8 && TMath::Abs(thePatElectron.eta()) < 1.479)
			mvaThresh = 0.57;
		else
			mvaThresh = 0.05;
		
		if(mvaVal < mvaThresh) return false;
	
    	return true;
}*/

std::vector<const pat::Electron* > tools::SSSyncElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
	double v_electron_dxy = 0.05;
	bool verbose = false;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
		//std::cout<<"\nNew Electron\n";
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        if(verbose) std::cout<<"1\n";
        //if( el->pt() < 7 ) continue;
		if(verbose) std::cout<<"2 - pt = "<<el->pt()<<"\n";
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
		if(verbose) std::cout<<"3 - eta = "<<el->eta()<<"\n";
        //if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        if(verbose) std::cout<<"4\n";
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;
		if(verbose) std::cout<<"5\n";
		if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_dxy  )  continue;
		if(verbose) std::cout<<"6\n";
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        if(verbose) std::cout<<"7\n";
        if( TMath::Abs( el->eta()) < 1.479  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.8  ) continue;
			if(verbose) std::cout<<"8_1\n";
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
			//if(verbose) std::cout<<"9_1 - sigmaieta = "<<TMath::Abs(el->scSigmaIEtaIEta())<<"\n";
            //if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
			if(verbose) std::cout<<"9_1 - sigmaieta = "<<TMath::Abs(el->full5x5_sigmaIetaIeta())<<"\n";
            if( TMath::Abs(el->full5x5_sigmaIetaIeta()) > 0.01 ) continue;
			if(verbose) std::cout<<"10_1\n";
            if( TMath::Abs(el->hadronicOverEm())  > 0.15  ) continue;  //recommended is 0.12 but HLT applies 0.1
			//if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.053 ) continue;
			if(verbose) std::cout<<"11_1\n";
	    }
        else if( TMath::Abs( el->eta() ) < 2.5 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.7 ) continue;
			if(verbose) std::cout<<"8_2\n";
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) continue;
			if(verbose) std::cout<<"9_2\n";
            if( TMath::Abs(el->full5x5_sigmaIetaIeta()) > 0.03 ) continue;
			if(verbose) std::cout<<"10_2\n";
            //if( TMath::Abs(el->hadronicOverEm()) > 0.099 ) continue;   /// at the HLT 0.075  recommended is 0.1
			//if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.11 ) continue;
		}
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}


bool tools::SSSyncIsTight(const pat::Electron  thePatElectron,
                                                              double v_electron_pt,
                                                              reco::Vertex PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS,float sip)
															  //reco::TransientTrack theTT)
{//
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
	double v_electron_dxy = 0.05;
	bool verbose = false;
    
        const reco::GsfTrackRef gsfTrack = thePatElectron.gsfTrack();
        if (!gsfTrack.isNonnull()) {
            return false;
        }
        
		if(verbose) std::cout<<"\n\n1: pt = "<<thePatElectron.pt()<<"\n";
		
        //if( thePatElectron.pt() < v_electron_pt ) return false;
		if(verbose) std::cout<<"2\n";
        if( TMath::Abs(thePatElectron.eta()) > v_electron_eta ) return false;
		if(verbose) std::cout<<"3\n";
        if( bool_electron_chargeConsistency && !thePatElectron.isGsfCtfScPixChargeConsistent() )  return false;
        if(verbose) std::cout<<"4\n";
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (thePatElectron), theConversions, BS);
        if( vtxFitConversion )  return false;
        if(verbose) std::cout<<"5 :::: "<<gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<<"\n";
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) return false;
		if(verbose) std::cout<<"6\n";
		if( TMath::Abs(gsfTrack->dxy(PV.position())) > v_electron_dxy  )  return false;
		if(verbose) std::cout<<"7\n";
        if( TMath::Abs(gsfTrack->dz(PV.position())) > v_electron_dz  ) return false;
        if(verbose) std::cout<<"8\n";
        if( TMath::Abs( thePatElectron.eta()) < 1.479  ) {
            if( TMath::Abs(thePatElectron.deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) return false;
			if(verbose) std::cout<<"9_1\n";
            if( TMath::Abs(thePatElectron.deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) return false;
			if(verbose) std::cout<<"10_1\n";
            if( TMath::Abs(thePatElectron.full5x5_sigmaIetaIeta()) > 0.01 ) return false;
			if(verbose) std::cout<<"11_1\n";
            if( TMath::Abs(thePatElectron.hadronicOverEm())  > 0.12  ) return false;  //recommended is 0.12 but HLT applies 0.1
			if(verbose) std::cout<<"12_1\n";
			if( TMath::Abs(1.0/thePatElectron.ecalEnergy() - thePatElectron.eSuperClusterOverP()/thePatElectron.ecalEnergy()) > 0.05 ) return false;
			
			
	    }
        else if( TMath::Abs( thePatElectron.eta() ) < 2.5 ) {
            if( TMath::Abs(thePatElectron.deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) return false;
			if(verbose) std::cout<<"9_2\n";
            if( TMath::Abs(thePatElectron.deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) return false;
			if(verbose) std::cout<<"10_2\n";
            if( TMath::Abs(thePatElectron.full5x5_sigmaIetaIeta()) > 0.03 ) return false;
			if(verbose) std::cout<<"11_2\n";
            if( TMath::Abs(thePatElectron.hadronicOverEm()) > 0.1 ) return false;   /// at the HLT 0.075  recommended is 0.1
			if(verbose) std::cout<<"12_2\n";
			if( TMath::Abs(1.0/thePatElectron.ecalEnergy() - thePatElectron.eSuperClusterOverP()/thePatElectron.ecalEnergy()) > 0.05 ) return false;
			
		}
		if(verbose) std::cout<<"13\n";
		//reco::TransientTrack tt  = theTTBuilder.build(gsfTrack);
       // Measurement1D ip3D = IPTools::absoluteImpactParameter3D(theTT, PV).second;
		//if(ip3D.significance() >= 4) return false;
		if(sip >= 4) return false;
		
   
    return true;
}


bool tools::SSSyncIsTightNoSIP(const pat::Electron  thePatElectron,
                                                              double v_electron_pt,
                                                              reco::Vertex PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
															  //reco::TransientTrack theTT)
{//
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
	double v_electron_dxy = 0.05;
	bool verbose = false;
    
        const reco::GsfTrackRef gsfTrack = thePatElectron.gsfTrack();
        if (!gsfTrack.isNonnull()) {
            return false;
        }
        
		if(verbose) std::cout<<"\n\n1: pt = "<<thePatElectron.pt()<<"\n";
		
        //if( thePatElectron.pt() < v_electron_pt ) return false;
		if(verbose) std::cout<<"2\n";
        if( TMath::Abs(thePatElectron.eta()) > v_electron_eta ) return false;
		if(verbose) std::cout<<"3\n";
        if( bool_electron_chargeConsistency && !thePatElectron.isGsfCtfScPixChargeConsistent() )  return false;
        if(verbose) std::cout<<"4\n";
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (thePatElectron), theConversions, BS);
        if( vtxFitConversion )  return false;
        if(verbose) std::cout<<"5 :::: "<<gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<<"\n";
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) return false;
		if(verbose) std::cout<<"6\n";
		if( TMath::Abs(gsfTrack->dxy(PV.position())) > v_electron_dxy  )  return false;
		if(verbose) std::cout<<"7\n";
        if( TMath::Abs(gsfTrack->dz(PV.position())) > v_electron_dz  ) return false;
        if(verbose) std::cout<<"8\n";
        if( TMath::Abs( thePatElectron.eta()) < 1.479  ) {
            if( TMath::Abs(thePatElectron.deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) return false;
			if(verbose) std::cout<<"9_1\n";
            if( TMath::Abs(thePatElectron.deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) return false;
			if(verbose) std::cout<<"10_1\n";
            if( TMath::Abs(thePatElectron.full5x5_sigmaIetaIeta()) > 0.01 ) return false;
			if(verbose) std::cout<<"11_1\n";
            if( TMath::Abs(thePatElectron.hadronicOverEm())  > 0.12  ) return false;  //recommended is 0.12 but HLT applies 0.1
			if(verbose) std::cout<<"12_1\n";
			if( TMath::Abs(1.0/thePatElectron.ecalEnergy() - thePatElectron.eSuperClusterOverP()/thePatElectron.ecalEnergy()) > 0.05 ) return false;
			
			
	    }
        else if( TMath::Abs( thePatElectron.eta() ) < 2.5 ) {
            if( TMath::Abs(thePatElectron.deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) return false;
			if(verbose) std::cout<<"9_2\n";
            if( TMath::Abs(thePatElectron.deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) return false;
			if(verbose) std::cout<<"10_2\n";
            if( TMath::Abs(thePatElectron.full5x5_sigmaIetaIeta()) > 0.03 ) return false;
			if(verbose) std::cout<<"11_2\n";
            if( TMath::Abs(thePatElectron.hadronicOverEm()) > 0.1 ) return false;   /// at the HLT 0.075  recommended is 0.1
			if(verbose) std::cout<<"12_2\n";
			if( TMath::Abs(1.0/thePatElectron.ecalEnergy() - thePatElectron.eSuperClusterOverP()/thePatElectron.ecalEnergy()) > 0.05 ) return false;
			
		}
		if(verbose) std::cout<<"13\n";
		//reco::TransientTrack tt  = theTTBuilder.build(gsfTrack);
        //Measurement1D ip3D = IPTools::absoluteImpactParameter3D(theTT, PV).second;
		//if(ip3D.significance() >= 4) return false;
		
   
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> tools::SSSyncElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
    double v_electron_dxy = 0.05;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
		index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
       // if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
		
        //if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;
		if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_dxy  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( TMath::Abs( el->eta()) < 1.479  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.8  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.15  ) continue;  //recommended is 0.12 but HLT applies 0.1
			//if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.053 ) continue;
			
	    }
        else if( TMath::Abs( el->eta() ) < 2.5 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.7 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            //if( TMath::Abs(el->hadronicOverEm()) > 0.099 ) continue;   /// at the HLT 0.075  recommended is 0.1
			//if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.11 ) continue;
		}
        
        vElectrons.push_back(index);
	}
    return vElectrons;
}

std::vector<const pat::Muon* > tools::MVALooseMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;

    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        if ( mu->pt()  < 5 ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        bool got = (mu->isTrackerMuon() || mu->isGlobalMuon());
        if ( !got  ) continue;
        if ( !mu->isPFMuon() ) continue;
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        
        //if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        //if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

std::vector<const pat::Muon* > tools::SSSyncMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;

    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
	bool verbose = false;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
		if(verbose) std::cout<<"\n\nNew Muon pt = "<<mu->pt()<<"\n";
        //if ( mu->pt()  < 5 ) continue;
		if(verbose) std::cout<<"eta = "<<mu->eta()<<"\n";
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        bool got = (mu->isTrackerMuon() || mu->isGlobalMuon());
		if(verbose) std::cout<<"tracker or global\n";
        if ( !got  ) continue;
		if(verbose) std::cout<<"pf muon\n";
        if ( !mu->isPFMuon() ) continue;
       // if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
       // if( innerTrack.isNull() ) continue;//
       /// if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
       // if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
		//if( 
        
        //const reco::TrackRef globalTrack = mu->globalTrack() ;
        //if( globalTrack.isNull() ) continue;
       // if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
       // if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        if(verbose) std::cout<<"d0 = "<<TMath::Abs(innerTrack->dxy(PV))<<"\n";
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
		if(verbose) std::cout<<"dz = "<<TMath::Abs(innerTrack->dz(PV))<<"\n";
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

bool tools::MVAIsTightNoIsoSIP(const pat::Muon thePatMuon,
                                                      double v_muon_pt,
                                                      reco::Vertex PV,
                                                      double v_muon_d0)
													  //reco::TransientTrack theTT)//
{
  
    	//if( thePatMuon.pt() < 10) return false;
	
        //if ( TMath::Abs( thePatMuon.eta() ) > v_muon_eta ) return false;
       // if ( !thePatMuon.isGlobalMuon()  ) return false;
        //if ( !thePatMuon.isPFMuon() ) return false;
		//if( thePatMuon.globalTrack()->normalizedChi2() > 2 ) return false;
		//if(thePatMuon.combinedQuality().chi2LocalPosition > 11) return false;
		//if(thePatMuon.combinedQuality().trkKink > 19) return false;
		//if(thePatMuon.innerTrack()->validFraction() < 0.8) return false;
		//if((thePatMuon.innerTrack()->ptError())/(thePatMuon.innerTrack()->pt()) > 0.2) return false;
		//if(thePatMuon.segmentCompatibility() < 0.451) return false;
		//reco::TransientTrack tt  = theTTBuilder.build(innerTrack);
       // Measurement1D ip3D = IPTools::absoluteImpactParameter3D(theTT, PV).second;
		//if(ip3D.significance() >= 4) return false;
		
		if((thePatMuon.innerTrack()->ptError())/(thePatMuon.innerTrack()->pt()) > 0.2) return false;
		
		bool goodGlb = thePatMuon.isGlobalMuon() && thePatMuon.globalTrack()->normalizedChi2() < 3 
												 && thePatMuon.combinedQuality().chi2LocalPosition < 12 
												 && thePatMuon.combinedQuality().trkKink < 20;
		
		bool isTight = thePatMuon.innerTrack()->validFraction() >= 0.8 
						&& thePatMuon.segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
		
		
        
    return isTight;
}

bool tools::SSSyncIsTight(const pat::Muon thePatMuon,
                                                      double v_muon_pt,
                                                      reco::Vertex PV,
                                                      double v_muon_d0,float sip)
													  //reco::TransientTrack theTT)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;

    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
        //if ( thePatMuon.pt()  < v_muon_pt ) return false;
        if ( TMath::Abs( thePatMuon.eta() ) > v_muon_eta ) return false;
        if ( !thePatMuon.isGlobalMuon()  ) return false;
        if ( !thePatMuon.isPFMuon() ) return false;
        if ( thePatMuon.numberOfMatchedStations() < v_muon_numberOfMatchedStations ) return false;   //we should add this to skim
        
        const reco::TrackRef innerTrack = thePatMuon.innerTrack();
        if( innerTrack.isNull() ) return false;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) return false;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) return false;
		//if( 
        
        const reco::TrackRef globalTrack = thePatMuon.globalTrack() ;
        if( globalTrack.isNull() ) return false;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) return false;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) return false;
        
        if(TMath::Abs(innerTrack->dxy(PV.position())) > v_muon_d0  ) return false;
        if(TMath::Abs(innerTrack->dz(PV.position()))   > v_muon_dz  ) return false;
		
		//reco::TransientTrack tt  = theTTBuilder.build(innerTrack);
       // Measurement1D ip3D = IPTools::absoluteImpactParameter3D(theTT, PV).second;
		//if(ip3D.significance() >= 4) return false;
		if(sip >= 4) return false;
        
    return true;
}

bool tools::SSSyncIsTightNoSIP(const pat::Muon thePatMuon,
                                                      double v_muon_pt,
                                                      reco::Vertex PV,
                                                      double v_muon_d0)
													 // reco::TransientTrack theTT)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;

    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
        //if ( thePatMuon.pt()  < v_muon_pt ) return false;
        if ( TMath::Abs( thePatMuon.eta() ) > v_muon_eta ) return false;
        if ( !thePatMuon.isGlobalMuon()  ) return false;
        if ( !thePatMuon.isPFMuon() ) return false;
        if ( thePatMuon.numberOfMatchedStations() < v_muon_numberOfMatchedStations ) return false;   //we should add this to skim
        
        const reco::TrackRef innerTrack = thePatMuon.innerTrack();
        if( innerTrack.isNull() ) return false;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) return false;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) return false;
		//if( 
        
        const reco::TrackRef globalTrack = thePatMuon.globalTrack() ;
        if( globalTrack.isNull() ) return false;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) return false;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) return false;
        
        if(TMath::Abs(innerTrack->dxy(PV.position())) > v_muon_d0  ) return false;
        if(TMath::Abs(innerTrack->dz(PV.position()))   > v_muon_dz  ) return false;
		
		//reco::TransientTrack tt  = theTTBuilder.build(innerTrack);
       // Measurement1D ip3D = IPTools::absoluteImpactParameter3D(theTT, PV).second;
		//if(ip3D.significance() >= 4) return false;
        
    return true;
}

std::vector<int > tools::SSSyncMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<int> vMuons;
    int counter = -1;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        counter++;
        //if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
		bool got = (mu->isTrackerMuon() || mu->isGlobalMuon());
        if ( !got  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
       // if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
       // if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;
        
        vMuons.push_back(counter);
    }
    
    return vMuons;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const reco::Muon* > tools::EffsMuonSelector(const std::vector<reco::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;

    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<const reco::Muon* > vMuons;
    for( std::vector<reco::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
		//if( 
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        //if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
        //if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

std::vector<const pat::Muon* > tools::ssbMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;

    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
		//if( 
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        //if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
        //if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int > tools::ssbMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<int> vMuons;
    int counter = -1;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        counter++;
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
       // if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
       // if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;
        
        vMuons.push_back(counter);
    }
    
    return vMuons;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<const pat::Electron* > tools::SSOptimizerElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
    //bool bool_electron_ecalDriven = true;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        //if (!el->trackerDrivenSeed() ) continue;
        
        //if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        //if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;

        /*
        if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
        //if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) >= 0.053 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
        
        if( TMath::Abs( el->eta()) < 1.479  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
			if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
			if( TMath::Abs(gsfTrack->dxy(PV)) > 0.02  )  continue;
        	if( TMath::Abs(gsfTrack->dz(PV)) > 0.1  ) continue;
	    }
        else if( TMath::Abs( el->eta() ) < 2.5 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.030 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;   /// at the HLT 0.075  recommended is 0.1
			if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
			if( TMath::Abs(gsfTrack->dxy(PV)) > 0.02  )  continue;
       	 	if( TMath::Abs(gsfTrack->dz(PV)) > 0.1  ) continue;
		}
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}

std::vector<int> tools::SSOptimizerElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.03;
    //bool bool_electron_ecalDriven = true;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        //if (!el->trackerDrivenSeed() ) continue;

        
        //if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        //if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        /*if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
       // if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
        
        if( TMath::Abs( el->eta()) < 1.479  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
			if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
			if( TMath::Abs(gsfTrack->dxy(PV)) > 0.02  )  continue;
        	if( TMath::Abs(gsfTrack->dz(PV)) > 0.1  ) continue;
	    }
        else if( TMath::Abs( el->eta() ) < 2.5 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.030 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.1 ) continue;   /// at the HLT 0.075  recommended is 0.1
			if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
			if( TMath::Abs(gsfTrack->dxy(PV)) > 0.02  )  continue;
       	 	if( TMath::Abs(gsfTrack->dz(PV)) > 0.1  ) continue;
		}
        
        vElectrons.push_back(index);
	}
    return vElectrons;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Electron* > tools::ssbElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.03;
    //bool bool_electron_ecalDriven = true;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        //if (!el->trackerDrivenSeed() ) continue;
        
        //if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        //if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;

        /*
        if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
        //if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) >= 0.053 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;
        
        if( TMath::Abs( el->eta()) < 1.479  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.051  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.015 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.10  ) continue;  //recommended is 0.12 but HLT applies 0.1
			if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.053 ) continue;
			if( TMath::Abs(gsfTrack->dxy(PV)) > 0.012  )  continue;
        	if( TMath::Abs(gsfTrack->dz(PV)) > 0.030  ) continue;
	    }
        else if( TMath::Abs( el->eta() ) < 2.5 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.056 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.023 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.030 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.099 ) continue;   /// at the HLT 0.075  recommended is 0.1
			if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.11 ) continue;
			if( TMath::Abs(gsfTrack->dxy(PV)) > 0.068  )  continue;
       	 	if( TMath::Abs(gsfTrack->dz(PV)) > 0.78  ) continue;
		}
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> tools::ssbElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.03;
    //bool bool_electron_ecalDriven = true;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        //if (!el->trackerDrivenSeed() ) continue;

        
        //if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        //if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        /*if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
       // if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;
        
        if( TMath::Abs( el->eta()) < 1.479  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.051  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.015 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.10  ) continue;  //recommended is 0.12 but HLT applies 0.1
			if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.053 ) continue;
			if( TMath::Abs(gsfTrack->dxy(PV)) > 0.012  )  continue;
        	if( TMath::Abs(gsfTrack->dz(PV)) > 0.030  ) continue;
	    }
        else if( TMath::Abs( el->eta() ) < 2.5 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.056 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.023 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.030 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.099 ) continue;   /// at the HLT 0.075  recommended is 0.1
			if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.11 ) continue;
			if( TMath::Abs(gsfTrack->dxy(PV)) > 0.068  )  continue;
       	 	if( TMath::Abs(gsfTrack->dz(PV)) > 0.78  ) continue;
	    }
        
        vElectrons.push_back(index);
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Muon* > tools::ssbMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0,
                                                      std::vector<const pat::Jet* > SelectedBJets)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
        if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;

        TLorentzVector pLep; pLep.SetPtEtaPhiE(mu->pt(), mu->eta(), mu->phi(), mu->energy());
        bool RemoveL = false;
        for (unsigned int nj = 0; nj!=SelectedBJets.size(); ++nj) {
            TLorentzVector pJet; pJet.SetPtEtaPhiM( SelectedBJets[nj]->pt(), SelectedBJets[nj]->eta(), SelectedBJets[nj]->phi(), 0 );
            double ang = pLep.DeltaR( pJet );
            if (ang < 0.4) {
                RemoveL = true;
                break;
            }
        }
        if (!RemoveL)
            vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int > tools::ssbMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                              double v_muon_pt,
                                              reco::Vertex::Point PV,
                                              double v_muon_d0,
                                              std::vector<const pat::Jet* > SelectedBJets)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.1;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    double v_muon_emVetoEt=4.;
    double v_muon_hadVetoEt= 6.;
    
    std::vector<int> vMuons;
    int counter = -1;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ )
	{
        counter++;
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > v_muon_dz  ) continue;
        
        if( mu->isolationR03().emVetoEt > v_muon_emVetoEt ) continue;
        if( mu->isolationR03().hadVetoEt> v_muon_hadVetoEt ) continue;
        
        TLorentzVector pLep; pLep.SetPtEtaPhiE(mu->pt(), mu->eta(), mu->phi(), mu->energy());
        bool RemoveL = false;
        for (unsigned int nj = 0; nj!=SelectedBJets.size(); ++nj) {
            TLorentzVector pJet; pJet.SetPtEtaPhiM( SelectedBJets[nj]->pt(), SelectedBJets[nj]->eta(), SelectedBJets[nj]->phi(), 0 );
            double ang = pLep.DeltaR( pJet );
            if (ang < 0.4) {
                RemoveL = true;
                break;
            }
        }
        if (!RemoveL)
        vMuons.push_back(counter);
    }
    
    return vMuons;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Electron* > tools::ssbElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS,
                                                              std::vector<const pat::Jet* > SelectedBJets)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
    //bool bool_electron_ecalDriven = true;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        //if (!el->trackerDrivenSeed() ) continue;
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        /*if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
        
        if( TMath::Abs( el->eta()) < 1.5  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.1  ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.075 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }
        
        TLorentzVector pLep; pLep.SetPtEtaPhiE(el->pt(), el->eta(), el->phi(), el->energy());
        bool RemoveL = false;
        for (unsigned int nj = 0; nj!=SelectedBJets.size(); ++nj) {
            TLorentzVector pJet; pJet.SetPtEtaPhiM( SelectedBJets[nj]->pt(), SelectedBJets[nj]->eta(), SelectedBJets[nj]->phi(), 0 );
            double ang = pLep.DeltaR( pJet );
            if (ang < 0.4) {
                RemoveL = true;
                break;
            }
        }
        if (!RemoveL)
            vElectrons.push_back(&*el );
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<int> tools::ssbElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                 double v_electron_pt,
                                                 reco::Vertex::Point PV,
                                                 double v_electron_d0,
                                                 bool bool_electron_chargeConsistency,
                                                 edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                 reco::BeamSpot::Point BS,
                                                 std::vector<const pat::Jet* > SelectedBJets)
{
    double v_electron_eta = 2.4;
    double v_electron_dz = 0.1;
    //bool bool_electron_ecalDriven = true;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) > 1.4442 && TMath::Abs(el->superCluster()->eta()) < 1.566 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        
        if (!el->trackerDrivenSeed() ) continue;
        
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        /*if( el->pt() < 20. ){
            if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
        }*/
        
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
		
        if( TMath::Abs( el->eta()) < 1.5  ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.1  ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.075 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }
        
        TLorentzVector pLep; pLep.SetPtEtaPhiE(el->pt(), el->eta(), el->phi(), el->energy());
        bool RemoveL = false;
        for (unsigned int nj = 0; nj!=SelectedBJets.size(); ++nj) {
            TLorentzVector pJet; pJet.SetPtEtaPhiM( SelectedBJets[nj]->pt(), SelectedBJets[nj]->eta(), SelectedBJets[nj]->phi(), 0 );
            double ang = pLep.DeltaR( pJet );
            if (ang < 0.4) {
                RemoveL = true;
                break;
            }
        }
        if (!RemoveL)
            vElectrons.push_back(index);
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Muon Selector

std::vector<const pat::Muon* > tools::fakeMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.2;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}


std::vector<int> tools::fakeMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                                       double v_muon_pt,
                                                       reco::Vertex::Point PV,
                                                       double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.2;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    std::vector<int > vMuons;
    int index = -1;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        index++;
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(index);
    }
    
    return vMuons;
}


std::vector<const pat::Muon* > tools::ewkMuonSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.2;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
        
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;

        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(&*mu);
    }
    
    return vMuons;
}

std::vector<int> tools::ewkMuonSelectorIndex(const std::vector<pat::Muon>  & thePatMuons,
                                                      double v_muon_pt,
                                                      reco::Vertex::Point PV,
                                                      double v_muon_d0)
{
    double v_muon_eta = 2.4;
    double v_muon_dz = 0.5;
    
    int v_muon_numberOfMatchedStations = 2;
    int v_muon_nPixelValidHits = 1;
    int v_muon_numberOfValidMuonHits = 1;
    int v_muon_nValidHits = 6;
    double v_muon_chi2Norm = 10.;
    
    int index = -1;
    std::vector<int > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        index++;
        if ( mu->pt() < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        //if ( !mu->isTrackerMuon() ) continue;
        if ( !mu->isGlobalMuon()  ) continue;
        if ( !mu->isPFMuon() ) continue;
        if ( mu->numberOfMatchedStations() < v_muon_numberOfMatchedStations ) continue;   //we should add this to skim
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        if( innerTrack->hitPattern().trackerLayersWithMeasurement() < v_muon_nValidHits ) continue;
        if( innerTrack->hitPattern().numberOfValidPixelHits() < v_muon_nPixelValidHits  ) continue;
        
        const reco::TrackRef globalTrack = mu->globalTrack() ;
        if( globalTrack.isNull() ) continue;
        if( globalTrack->normalizedChi2() > v_muon_chi2Norm ) continue;
        if( globalTrack->hitPattern().numberOfValidMuonHits() < v_muon_numberOfValidMuonHits ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > v_muon_d0  ) continue;
        if(TMath::Abs(innerTrack->dz(PV)) > v_muon_dz  ) continue;
        
        vMuons.push_back(index);
	}
    return vMuons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<const pat::Electron* > tools::ewkElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    double v_electron_eta=2.4;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        } 
     /*   std::cout<<TMath::Abs(el->eta())<<std::endl;
        std::cout<<TMath::Abs(el->superCluster()->eta())<<std::endl;
        std::cout<<el->pt()<<std::endl;
        std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
        std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;
        std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;
        std::cout<<"Conversion  "<<ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS)<<std::endl;
        std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
        
        std::cout<<"********"<<std::endl;*/
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        //if (!el->trackerDrivenSeed() ) continue;

        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;

        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
    
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
                
        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;
        
        if( TMath::Abs( el->eta()) < 1.5  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else if( TMath::Abs( el->eta() ) < 2.4 )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
        }
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<int> tools::ewkElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                           double v_electron_pt,
                                                           reco::Vertex::Point PV,
                                                           double v_electron_d0,
                                                           bool bool_electron_chargeConsistency,
                                                           edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                           reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    double v_electron_eta=2.4;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;
        //if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        //if (!el->trackerDrivenSeed() ) continue;

        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 1 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;
        
        if( TMath::Abs( el->eta()) < 1.5  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.15  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
	    }
        else if( TMath::Abs( el->eta() ) < 2.4 )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.10 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.009 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
	    }
        
        vElectrons.push_back(index);
	}
    return vElectrons;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<const pat::Electron* > tools::fakeElectronSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                              double v_electron_pt,
                                                              reco::Vertex::Point PV,
                                                              double v_electron_d0,
                                                              bool bool_electron_chargeConsistency,
                                                              edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                              reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.1;
    double v_electron_eta=2.5;
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        /*   std::cout<<TMath::Abs(el->eta())<<std::endl;
         std::cout<<TMath::Abs(el->superCluster()->eta())<<std::endl;
         std::cout<<el->pt()<<std::endl;
         std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
         std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;
         std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;
         std::cout<<"Conversion  "<<ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS)<<std::endl;
         std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
         
         std::cout<<"********"<<std::endl;*/
        
        if( el->pt() < v_electron_pt ) continue;
        
        //std::cout<<TMath::Abs(el->eta())<<std::endl;
        
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;

        //std::cout<<TMath::Abs(gsfTrack->dxy(PV))<<std::endl;
        //std::cout<<TMath::Abs(gsfTrack->dz(PV))<<std::endl;

        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;

        //std::cout<<TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy())<<std::endl;

        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        //std::cout<<vtxFitConversion<<std::endl;
        if( vtxFitConversion )  continue;
        //std::cout<<gsfTrack->trackerExpectedHitsInner().numberOfHits()<<std::endl;
        
        //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
        
        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
        }
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Electron Selector
std::vector<int> tools::fakeElectronSelectorIndex(const std::vector<pat::Electron>  & thePatElectrons,
                                                 double v_electron_pt,
                                                 reco::Vertex::Point PV,
                                                 double v_electron_d0,
                                                 bool bool_electron_chargeConsistency,
                                                 edm::Handle< std::vector<reco::Conversion> > &theConversions,
                                                 reco::BeamSpot::Point BS)
{
    double v_electron_dz = 0.1;
    double v_electron_eta=2.5;
    
    int index = -1;
    std::vector<int > vElectrons;
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        index++;
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        //if( TMath::Abs(el->superCluster()->eta()) < 1.566 &&  TMath::Abs(el->superCluster()->eta()) > 1.4442 ) continue;
        
        
        if( TMath::Abs(gsfTrack->dxy(PV)) > v_electron_d0  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV))  > v_electron_dz  ) continue;
        
        if( bool_electron_chargeConsistency && !el->isGsfCtfScPixChargeConsistent() )  continue;
        
        //if( el->pt() < 20. )
	    //{
        //    if( ! ( el->fbrem() > 0.15 || ( TMath::Abs( el->superCluster()->eta() ) < 1.0 && el->eSuperClusterOverP() > 0.95 ) ) ) continue;
	    //}
        
        
        //if( TMath::Fabs(( 1./ el->trackMomentumAtVtx().p() ) - ( 1./ el->ecalEnergy() ) ) > 0.05 ) continue;
        if( TMath::Abs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy()) > 0.05 ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        
	    //if( gsfTrack->trackerExpectedHitsInner().numberOfHits() > 0 ) continue;
		if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0 ) continue;
        
        if( TMath::Abs( el->superCluster()->eta()) < 1.479  )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.06  ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.004 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm())  > 0.12  ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else //if( TMath::Abs( el->superCluster()->eta() ) < v_electron_eta )
	    {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.03 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.10 ) continue;   /// at the HLT 0.075  recommended is 0.1
        }
        
        vElectrons.push_back(index);
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<const pat::Jet* > tools::JetSelectorAll(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta)
{
    bool    bool_jet_id= false;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
       /* if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }*/
        vJets.push_back( &*jet );
        
        /*unsigned int nConst = jet->getPFConstituents().size();
         std::cout<<"Number of constituents "<<nConst<<std::endl;
         for (unsigned int i=0; i!=nConst; ++i) {
         std::cout<<jet->getPFConstituent(i)->reco::LeafCandidate::vz()<<std::endl;
         }*/
        
        
    }
    return vJets;
}

std::vector<int> tools::JetSelectorIndexAll(const std::vector<pat::Jet>  & thePatJets,
                                         double  v_jet_pt,
                                         double  v_jet_eta)
{
    bool    bool_jet_id= false;
    
    std::vector<int> vJets;
    
    int i = -1;
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        i++;
       /* if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }*/
        vJets.push_back( i );
        
    }
    return vJets;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta)
{
    bool    bool_jet_id= true;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
		//std::cout<<"jet pt/eta/phi = "<<jet->pt()<<"/"<<jet->eta()<<"/"<<jet->phi()<<"\n";
        if( jet->pt() < v_jet_pt )continue;
		//std::cout<<"1\n";
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
			//std::cout<<"2\n";
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
			//std::cout<<"3\n";
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
			//std::cout<<"4\n";
            if( ( jet->neutralMultiplicity() + jet->chargedMultiplicity() ) < 2 ) continue;
			//std::cout<<"5\n";
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() < 0.000001 ) continue;
				//std::cout<<"5\n";
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
				//std::cout<<"6\n";
                if( jet->chargedMultiplicity() < 1 ) continue;
				//std::cout<<"7\n";
            }
	    }
        vJets.push_back( &*jet );
        
        /*unsigned int nConst = jet->getPFConstituents().size();
        std::cout<<"Number of constituents "<<nConst<<std::endl;
        for (unsigned int i=0; i!=nConst; ++i) {
            std::cout<<jet->getPFConstituent(i)->reco::LeafCandidate::vz()<<std::endl;
        }*/
        
        
    }
    return vJets;
}

std::vector<int> tools::JetSelectorIndex(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta)
{
    bool    bool_jet_id= true;
    
    std::vector<int> vJets;
    
    int i = -1;
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        i++;
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralMultiplicity() + jet->chargedMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0.000001 ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() < 1 ) continue;
            }
	    }
        vJets.push_back( i );
        
    }
    return vJets;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons)
{
    bool    bool_jet_id= true;
    double  v_jet_leptonVetoDR=0.4;
    //double  v_jet_leptonVetoDR = -1.;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_jet_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
	    
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}

std::vector<int> tools::JetSelectorIndex(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons)
{
    bool    bool_jet_id= true;
    double  v_jet_leptonVetoDR=0.4;
    //double  v_jet_leptonVetoDR = -1.;
    
    std::vector<int> vJets;
    
    int index = -1;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        index++;
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_jet_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
	    
        if( vetoJet ) continue;
	    
        
        vJets.push_back( index );
    }
    return vJets;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Jet* > tools::JetSelector(const std::vector<pat::Jet>  & thePatJets,
                                                 double  v_jet_pt,
                                                 double  v_jet_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons,
                                                 std::map<const reco::PFTau*, int> SelectedTaus)
{
    bool    bool_jet_id= true;
    double  v_jet_leptonVetoDR=0.4;
    
    std::vector< const pat::Jet* > vJets;
    
    for( std::vector<pat::Jet>::const_iterator jet = thePatJets.begin(); jet != thePatJets.end(); jet++ )
	{
        if( jet->pt() < v_jet_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_jet_eta) continue;
        if( bool_jet_id )
	    {
            if( jet->neutralHadronEnergyFraction() >= 0.99 ) continue;
            if( jet->neutralEmEnergyFraction() >= 0.99 ) continue;
            if( ( jet->neutralHadronMultiplicity() + jet->chargedHadronMultiplicity() ) < 2 ) continue;
            if( TMath::Abs( jet->eta() ) < 2.4 )
            {
                if( jet->chargedHadronEnergyFraction() == 0. ) continue;
                if( jet->chargedEmEnergyFraction() >= 0.99 ) continue;
                if( jet->chargedMultiplicity() == 0 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_jet_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }

        for(std::map<const reco::PFTau*, int >::iterator it = SelectedTaus.begin() ; it != SelectedTaus.end() ;it++ ){
            
            const reco::PFTau *itau = it->first;

            float dphi = TMath::ACos( TMath::Cos( itau->phi()-jet->phi() ) );
            float deta = itau->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_jet_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Photon* > tools::PhotonSelector(const std::vector<pat::Photon>  & thePatPhotons,
                                                 double  v_photon_pt,
                                                 double  v_photon_eta,
                                                 std::vector<const pat::Electron*> vElectrons,
                                                 std::vector<const pat::Muon*> vMuons,
                                                 std::map<const reco::PFTau*, int> SelectedTaus)
{
    bool    bool_photon_id= true;
    double  v_photon_leptonVetoDR=0.2;
    
    std::vector< const pat::Photon* > vJets;
    
    for( std::vector<pat::Photon>::const_iterator jet = thePatPhotons.begin(); jet != thePatPhotons.end(); jet++ )
	{
        //std::cout<<jet->pt()<<" "<<jet->eta()<<" "<<jet->hadTowOverEm()<<" "<<jet->sigmaIetaIeta()<<std::endl;
        if( jet->pt() < v_photon_pt )continue;
        if( TMath::Abs( jet->eta() ) > v_photon_eta) continue;
        if( bool_photon_id )
	    {
            if( jet->hadTowOverEm() >= 0.05 ) continue;
            if( TMath::Abs(jet->sigmaIetaIeta()) > 0.012 ) continue;
            if( TMath::Abs( jet->eta() ) < 1.479 )
            {
                if( jet->hadTowOverEm() >= 0.05 ) continue;
                if( TMath::Abs(jet->sigmaIetaIeta()) > 0.034 ) continue;
            }
	    }
        
        bool vetoJet = false;
        for(unsigned int i = 0 ; i < vMuons.size() ;i++ )
        {
            const pat::Muon *mu = vMuons[i];
            float dphi = TMath::ACos( TMath::Cos( mu->phi()-jet->phi() ));
            float deta = mu->eta()-jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta) ;
            
            if(dr < v_photon_leptonVetoDR )
            {
                vetoJet = true;
                break;
            }
	    }
        if( vetoJet ) continue;
        
        for(unsigned int i = 0 ; i < vElectrons.size() ;i++ )
        {
            const pat::Electron *el = vElectrons[i];
            float dphi = TMath::ACos( TMath::Cos( el->phi()-jet->phi() ) );
            float deta = el->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_photon_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        for(std::map<const reco::PFTau*, int >::iterator it = SelectedTaus.begin() ; it != SelectedTaus.end() ;it++ ){
            
            const reco::PFTau *itau = it->first;
            
            float dphi = TMath::ACos( TMath::Cos( itau->phi()-jet->phi() ) );
            float deta = itau->eta() - jet->eta();
            float dr = TMath::Sqrt( dphi*dphi + deta*deta );
            
            if(dr < v_photon_leptonVetoDR)
            {
                vetoJet = true;
                break;
            }
	    }
        
        if( vetoJet ) continue;
	    
        
        vJets.push_back( &*jet );
    }
    return vJets;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Isolated Muon Selector
std::vector<const pat::Muon* > tools::MuonSelector_Iso(const std::vector< const pat::Muon *>  & thePatMuons, double v_muon_reliso, bool usePFiso, int Iso)
{
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<const pat::Muon *>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        
        float pfRelIsoMu = pfRelIso(*mu);
        float detRelIso = ((*mu)->isolationR03().emEt + (*mu)->isolationR03().hadEt + (*mu)->isolationR03().sumPt) / (*mu)->pt() ;
        float RelIso = detRelIso;
        if(usePFiso) RelIso = pfRelIsoMu;
        if( RelIso  < v_muon_reliso && Iso==0) continue;  //NonIso
        if( RelIso  > v_muon_reliso && Iso==1) continue;  //Iso
        
        vMuons.push_back(*mu);
    }
    
    return vMuons;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Isolated Electron Selector

std::vector<const pat::Electron* > tools::ElectronSelector_Iso(const std::vector<const pat::Electron *>  & thePatElectrons,
                                                               double v_electron_reliso, bool usePFiso, double myRho, int Iso)
{
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<const pat::Electron *>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        
        float ecalIso = TMath::Abs((*el)->eta()) > 1.47 ? (*el)->dr03EcalRecHitSumEt() : TMath::Max((*el)->dr03EcalRecHitSumEt()-1.,0.);
        float detRelIso = ((*el)->dr03TkSumPt() + (*el)->dr03HcalTowerSumEt() + ecalIso ) / (*el)->pt() ;
        float pfRelIsoEl = pfRelIso(*el, myRho);
        float RelIso = detRelIso;
        if(usePFiso) RelIso = pfRelIsoEl;
        if( RelIso  < v_electron_reliso && Iso==0) continue;  //NonIso
        if( RelIso  > v_electron_reliso && Iso==1) continue;  //Iso
        
        vElectrons.push_back(*el );
	}
    return vElectrons;
}

// clean vs selected muons
std::vector<const pat::Electron* > tools::ElectronSelector_Iso(const std::vector<const pat::Electron *>  & thePatElectrons,
                                                               double v_electron_reliso, bool usePFiso, double myRho, int Iso,
                                                               std::vector<const pat::Muon* > & thePatMuons)
{
    
    std::vector<const pat::Electron* > vElectrons;
    for( std::vector<const pat::Electron *>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ )
	{
        
        float ecalIso = TMath::Abs((*el)->eta()) > 1.47 ? (*el)->dr03EcalRecHitSumEt() : TMath::Max((*el)->dr03EcalRecHitSumEt()-1.,0.);
        float detRelIso = ((*el)->dr03TkSumPt() + (*el)->dr03HcalTowerSumEt() + ecalIso ) / (*el)->pt() ;
        float pfRelIsoEl = pfRelIso(*el, myRho);
        float RelIso = detRelIso;
        if(usePFiso) RelIso = pfRelIsoEl;
        if( RelIso  < v_electron_reliso && Iso==0) continue;  //NonIso
        if( RelIso  > v_electron_reliso && Iso==1) continue;  //Iso
        
        bool Remove = false;
        TLorentzVector Ele; Ele.SetPtEtaPhiE( (*el)->pt(), (*el)->eta(), (*el)->phi(), (*el)->energy() );
        for( std::vector<const pat::Muon *>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
            TLorentzVector Mu; Mu.SetPtEtaPhiE( (*mu)->pt(), (*mu)->eta(), (*mu)->phi(), (*mu)->energy() );
            float dR = Mu.DeltaR( Ele );
            if( dR < 0.1 ) {
                Remove = true;
                break;
            }
        }
        if (!Remove)
            vElectrons.push_back(*el );
	}
    return vElectrons;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tools::MT_calc(TLorentzVector Vect, double MET, double MET_Phi){
    
    double MT=sqrt(2* Vect.Pt() * MET * ( 1 - (TMath::Cos(Vect.Phi() - MET_Phi )) ) );
    
    return MT;
}

double tools::Mll_calc(TLorentzVector Vect1, TLorentzVector Vect2){
    return (Vect1 + Vect2).Mag();
}

//ID=#MET+5*(#MT+3*(#Mll+3*#category))
int tools::srID(double met, double mt, double mll, double channel) {
    if ((channel > 0) && (channel!=4))
        return TMath::Min(int(met/50),4) + 5*(int(mt>120)+int(mt>160) + 3*(2*int(mll>100) + 3*channel));
    else 
        return TMath::Min(int(met/50),4) + 5*(int(mt>120)+int(mt>160) + 3*(int(mll>75) + int(mll>105) + 3*channel));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tools::RelIso_El(const pat::Electron *el, int usePFiso, double myRho){
    double reliso= 999.;
    //double  Aeff[ 7 ] = { 0.10, 0.12, 0.085, 0.11, 0.12, 0.12, 0.13  };
    double  Aeff[ 7 ] = { 0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14  };
    
    float ecalIso = TMath::Abs(el->eta()) > 1.47 ? el->dr03EcalRecHitSumEt() : TMath::Max(el->dr03EcalRecHitSumEt()-1.,0.);
    float detRelIso = (el->dr03TkSumPt() + el->dr03HcalTowerSumEt() + ecalIso ) / el->pt() ;
    
    // PF Isolation
    double CorrectedTerm=0.0;
    if( TMath::Abs( el->superCluster()->eta() ) < 1.0 ) CorrectedTerm = myRho * Aeff[ 0 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 1.0 && TMath::Abs( el->superCluster()->eta() ) < 1.479  )   CorrectedTerm = myRho * Aeff[ 1 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 1.479 && TMath::Abs( el->superCluster()->eta() ) < 2.0  )   CorrectedTerm = myRho * Aeff[ 2 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 2.0 && TMath::Abs( el->superCluster()->eta() ) < 2.2  )     CorrectedTerm = myRho * Aeff[ 3 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 2.2 && TMath::Abs( el->superCluster()->eta() ) < 2.3  )     CorrectedTerm = myRho * Aeff[ 4 ];
    else if( TMath::Abs( el->superCluster()->eta() ) > 2.3 && TMath::Abs( el->superCluster()->eta() ) < 2.4  )     CorrectedTerm = myRho * Aeff[ 5 ];
    else CorrectedTerm = myRho * Aeff[ 6 ];
   
    float pfRelIso = (el->chargedHadronIso() + TMath::Max(0.0, el->neutralHadronIso() + el->photonIso() - CorrectedTerm ) ) /el->pt() ;
    
    if(usePFiso) reliso = pfRelIso;
    if(!usePFiso) reliso = detRelIso;
    
    
    return reliso;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double tools::RelIso_Mu(const pat::Muon *mu, int usePFiso){
    double reliso= 999.;
    
    float beta_i    = mu->pfIsolationR03().sumPUPt;
    float pfRelIso  = (mu->pfIsolationR03().sumChargedHadronPt +
                       TMath::Max ( 0.0 ,(mu->pfIsolationR03().sumNeutralHadronEt +
                                          mu->pfIsolationR03().sumPhotonEt - 0.5 *beta_i )) ) /mu->pt() ;
    
    float detRelIso = (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt) / mu->pt() ;
    
    
    
    if(usePFiso) reliso = pfRelIso;
    if(!usePFiso) reliso = detRelIso;
    
    
    return reliso;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int tools::DY(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons){
    TLorentzVector lep1, lep2, dilep, lep3;  // lepton 4-vectors
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
        const pat::Electron *ele1 = thePatElectrons[i];
        for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
            const pat::Electron *ele2 = thePatElectrons[j];
            if(ele2->charge()==ele1->charge()) continue;
            lep1.SetPtEtaPhiE(ele1->pt(),ele1->eta(),ele1->phi(),ele1->energy());
            lep2.SetPtEtaPhiE(ele2->pt(),ele2->eta(),ele2->phi(),ele2->energy());
            dilep=lep1+lep2;
            if(dilep.M()<12) return 1;
        }
    }
	
    for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
        const pat::Muon *mu1 = thePatMuons[i];
        for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
            const pat::Muon *mu2 = thePatMuons[j];
            if(mu2->charge()==mu1->charge()) continue;
            lep1.SetPtEtaPhiE(mu1->pt(),mu1->eta(),mu1->phi(),mu1->energy());
            lep2.SetPtEtaPhiE(mu2->pt(),mu2->eta(),mu2->phi(),mu2->energy());
            dilep=lep1+lep2;
            if(dilep.M()<12) return 1;
        }
    }
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int tools::OSSF(std::vector<const pat::Electron*>  & thePatElectrons, std::vector<const pat::Muon*>  & thePatMuons){
    int ossf= 0;
    for(unsigned int i = 0   ; i < thePatElectrons.size() ;i++ ) {
        for(unsigned int j = i+1 ; j < thePatElectrons.size() ;j++ ) {
            const pat::Electron *ele1 = thePatElectrons[i];
            const pat::Electron *ele2 = thePatElectrons[j];
            
            if( ( ele1->charge() + ele2->charge() ) ==0 ) ossf = 1;
        }}
    
    for(unsigned int i = 0   ; i < thePatMuons.size() ;i++ ) {
        for(unsigned int j = i+1 ; j < thePatMuons.size() ;j++ ) {
            const pat::Muon *mu1 = thePatMuons[i];
            const pat::Muon *mu2 = thePatMuons[j];
            
            if( ( mu1->charge() + mu2->charge() ) ==0 ) ossf = 1;
        }}
    
    return ossf;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool tools::PassEventAcceptance(std::vector<const pat::Electron*>  & thePatElectrons,std::vector<const pat::Muon*>  & thePatMuons){
    int n20=0;
    int n10=0;
    for(unsigned int i = 0 ; i < thePatElectrons.size() ;i++ ) {
        const pat::Electron *ele1 = thePatElectrons[i];
        if(ele1->pt()>20){
            n20++;
            n10++;
        }
        else if(ele1->pt()>10) n10++;
    }
    for(unsigned int i = 0 ; i < thePatMuons.size() ;i++ ) {
        const pat::Muon *mu1 = thePatMuons[i];
        if(mu1->pt()>20){
            n20++;
            n10++;
        }
        else if(mu1->pt()>10) n10++;
    }
    if(n20>0 && n10==3) return kTRUE;
    else return kFALSE;
}

//*****************************************************************************************************************************
//**** Tau Selector ***********************************************************************************************************
//*****************************************************************************************************************************
std::map<const reco::PFTau*, int > tools::TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
                                                    double v_tau_pt,
                                                    double v_tau_eta,
                                                    edm::Handle<reco::PFTauDiscriminator> & electron,
                                                    edm::Handle<reco::PFTauDiscriminator> & muon,
                                                    edm::Handle<reco::PFTauDiscriminator> & iso,
                                                    edm::Handle<reco::PFTauDiscriminator> & decay){
    std::map<const reco::PFTau*, int> vTaus;
    
    for( unsigned  i=0; i<PFTaus->size(); i++ ) {
        
        //std::cout<<"Tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
        //std::cout<<"pt threshold"<<std::endl;
        
        if((*PFTaus)[i].pt()<v_tau_pt) continue;
        
        //std::cout<<"eta cut"<<std::endl;
        
        if(TMath::Abs((*PFTaus)[i].eta())>v_tau_eta) continue;
        
        reco::PFTauRef tauCandidate(PFTaus, i);
        //std::cout<<"electron discr."<<std::endl;
        if( (*electron)[tauCandidate] < 0.5 ) continue;
        
        //std::cout<<"electron discr."<<std::endl;
        if( (*muon)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"muon discr."<<std::endl;
        if( (*iso)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"iso discr."<<std::endl;
        if( (*decay)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"tau pass"<<std::endl;
        vTaus.insert(std::pair<const reco::PFTau*, int >(&((*PFTaus)[i]), i));
    }
    return vTaus;
}


std::map<const reco::PFTau*, int > tools::TauSelector(edm::Handle<reco::PFTauCollection>  & PFTaus,
                                                      double v_tau_pt,
                                                      double v_tau_eta,
                                                      edm::Handle<reco::PFTauDiscriminator> & electron,
                                                      edm::Handle<reco::PFTauDiscriminator> & muon,
                                                      edm::Handle<reco::PFTauDiscriminator> & iso,
                                                      edm::Handle<reco::PFTauDiscriminator> & decay,
                                                      std::vector<const pat::Muon*> & thePatMuons,
                                                      const std::vector<const pat::Electron*>  & thePatElectrons){
    std::map<const reco::PFTau*, int> vTaus;
    
    for( unsigned  i=0; i<PFTaus->size(); i++ ) {
        
        //std::cout<<"Tau "<<(*PFTaus)[i].pt()<<" "<< (*PFTaus)[i].eta()<<" "<< (*PFTaus)[i].phi()<<std::endl;
        //std::cout<<"pt threshold"<<std::endl;
        
        if((*PFTaus)[i].pt()<v_tau_pt) continue;
        
        //std::cout<<"eta cut"<<std::endl;
        
        if(TMath::Abs((*PFTaus)[i].eta())>v_tau_eta) continue;
        
        reco::PFTauRef tauCandidate(PFTaus, i);
        //std::cout<<"electron discr."<<std::endl;
        if( (*electron)[tauCandidate] < 0.5 ) continue;
        
        //std::cout<<"electron discr."<<std::endl;
        if( (*muon)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"muon discr."<<std::endl;
        if( (*iso)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"iso discr."<<std::endl;
        if( (*decay)[tauCandidate] < 0.5 ) continue;
        //std::cout<<"tau pass"<<std::endl;
        
        bool Remove = false;
        TLorentzVector Tau; Tau.SetPtEtaPhiE( (*PFTaus)[i].pt(), (*PFTaus)[i].eta(), (*PFTaus)[i].phi(), (*PFTaus)[i].energy() );
        for( std::vector<const pat::Muon *>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
            TLorentzVector Mu; Mu.SetPtEtaPhiE( (*mu)->pt(), (*mu)->eta(), (*mu)->phi(), (*mu)->energy() );
            float dR = Mu.DeltaR( Tau );
            if( dR < 0.1 ) {
                Remove = true;
                break;
            }
        }
        if (!Remove) {
            for( std::vector<const pat::Electron *>::const_iterator mu = thePatElectrons.begin() ; mu != thePatElectrons.end() ; mu++ ) {
                TLorentzVector Mu; Mu.SetPtEtaPhiE( (*mu)->pt(), (*mu)->eta(), (*mu)->phi(), (*mu)->energy() );
                float dR = Mu.DeltaR( Tau );
                if( dR < 0.1 ) {
                    Remove = true;
                    break;
                }
            }
        }
        if (!Remove)
            vTaus.insert(std::pair<const reco::PFTau*, int >(&((*PFTaus)[i]), i));
    }
    return vTaus;
}


std::vector<const pat::Electron* > tools::ssbElectronVetoSelector(const std::vector<pat::Electron>  & thePatElectrons,
                                                                  double Rho,
                                                                  reco::Vertex::Point PV,
                                                                  const char* cutName)
{
    double v_electron_pt = 0 ;
    TString CutName = cutName;
    if( CutName == "gammastar" ) v_electron_pt = 5;
    else if( CutName == "Z"    ) v_electron_pt = 10;
    
    double v_electron_eta = 2.4 ;
    double v_electron_d0 = 0.04;
    double v_electron_reliso = 0.2;
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    
    std::vector<const pat::Electron* > vElectrons;
    
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        
        //    if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
        //    if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        if(TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
        if(TMath::Abs(el->gsfTrack()->dz(PV)) > v_electron_dz  ) continue;
        
        if( TMath::Abs( el->eta() ) < 1.4442 ){
            //  if( el->isEB() ){
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.8 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.15 ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.7 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
        }
        
        float pfRelIsoLocal = pfRelIso(&*el, Rho) ;
        if( pfRelIsoLocal > v_electron_reliso ) continue;
        
        vElectrons.push_back(&*el );
	}
    return vElectrons;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<const pat::Muon* > tools::ssbMuonVetoSelector(const std::vector<pat::Muon>  & thePatMuons,
                                                          const char* cutName)
{
    //bool bool_pfIsolation = true;
    float v_muon_pt = 0;
    float v_muon_eta = 2.4;
    //  float v_muon_dz  = 999.;
    float v_muon_iso = 0.2;
    //  bool ZMass = false;
    //double ZMass = 0;
    
    TString CutName = cutName;
    
    if( CutName == "gammastar" ) v_muon_pt = 5;
    else if( CutName == "Z"    ) v_muon_pt = 10;
    
    std::vector<const pat::Muon* > vMuons;
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        if ( mu->isTrackerMuon() || mu->isGlobalMuon() ) {
            if ( !mu->isPFMuon() ) continue;            
            // PF Isolation
            float pfRelIsoLocal = pfRelIso(&*mu) ;
            if( pfRelIsoLocal > v_muon_iso ) continue;
            vMuons.push_back(&*mu);
        }
    }
    return vMuons;
}

bool tools::cleanUp(const pat::Muon* testMu,
                   std::vector<const pat::Muon* > allMu,
                   const char* cutName)
{
    
    TLorentzVector TLVec; TLVec.SetPtEtaPhiE( testMu->pt(), testMu->eta(), testMu->phi(), testMu->energy() );
    TString CutName = cutName;

    for (unsigned int i=0; i!=allMu.size(); ++i) {
        if (testMu->charge() + allMu[i]->charge() != 0) continue;
    
        TLorentzVector P1; P1.SetPtEtaPhiE( allMu[i]->pt(), allMu[i]->eta(), allMu[i]->phi(), allMu[i]->energy() );
        
        float Mass = ( P1 + TLVec ).M();
        if( CutName == "Z" ){
            if( Mass < 106. && Mass > 76. )  {
                return true;
            }
        }
        else if( CutName == "gammastar" ){
            if( Mass < 12. )  {
                return true;
            }
        }
    }
    
    return false;
}

bool tools::cleanUp(const pat::Electron* testMu,
                   std::vector<const pat::Electron* > allMu,
                   const char* cutName)
{
    
    TLorentzVector TLVec; TLVec.SetPtEtaPhiE( testMu->pt(), testMu->eta(), testMu->phi(), testMu->energy() );
    TString CutName = cutName;

    for (unsigned int i=0; i!=allMu.size(); ++i) {
        if (testMu->charge() + allMu[i]->charge() != 0) continue;
        
        TLorentzVector P1; P1.SetPtEtaPhiE( allMu[i]->pt(), allMu[i]->eta(), allMu[i]->phi(), allMu[i]->energy() );
        
        float Mass = ( P1 + TLVec ).M();
        if( CutName == "Z" ){
            if( Mass < 106. && Mass > 76. )  {
                return true;
            }
        }
        else if( CutName == "gammastar" ){
            if( Mass < 12. )  {
                return true;
            }
        }
    }
    
    return false;
}


/*void tools::removeZCandForEE(const std::vector<pat::Electron>  & thePatElectrons,
                             edm::Handle< std::vector<reco::Conversion> > &theConversions,
                             reco::BeamSpot::Point BS,
                             reco::Vertex::Point PV,
                             double  Rho,
                             std::vector< Leptons > & LepCand,
                             bool signalLepOnly,
                             const char* cutName ) {
    
    
    TString CutName = cutName;
    
    double ZMass = 0;
    
    
    std::pair<unsigned int , TString  > myPair;
    std::vector< pair<unsigned int,TString >> VecOfPair;
 
    double v_electron_pt = 0 ;
    
    
    if( CutName == "gammastar" ) v_electron_pt = 5;
    else if( CutName == "Z"    ) v_electron_pt = 10;
    
    
    double v_electron_eta = 2.4 ;
    double v_electron_d0 = 0.04;
    double v_electron_reliso = 0.2;
    double v_electron_dz = 0.2;
    //bool bool_electron_ecalDriven = true;
    
    for( std::vector<pat::Electron>::const_iterator el = thePatElectrons.begin() ; el != thePatElectrons.end() ; el++ ) {
        
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        
        if( el->pt() < v_electron_pt ) continue;
        if( TMath::Abs(el->eta()) > v_electron_eta ) continue;
        
        //    if( TMath::Abs( el->superCluster()->eta() ) < 1.566 &&  TMath::Abs( el->superCluster()->eta() ) > 1.4442 ) continue;
        //    if( bool_electron_ecalDriven && !el->ecalDrivenSeed() ) continue;
        if(TMath::Abs(el->gsfTrack()->dxy(PV)) > v_electron_d0  )  continue;
        if(TMath::Abs(el->gsfTrack()->dz(PV)) > v_electron_dz  ) continue;
        
        if( TMath::Abs( el->eta() ) < 1.4442 ){
            //  if( el->isEB() ){
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.8 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.007 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.01 ) continue;
            if( TMath::Abs(el->hadronicOverEm()) > 0.15 ) continue;  //recommended is 0.12 but HLT applies 0.1
        }
        else if( TMath::Abs( el->eta() ) < 2.4 ) {
            if( TMath::Abs(el->deltaPhiSuperClusterTrackAtVtx()) > 0.7 ) continue;
            if( TMath::Abs(el->deltaEtaSuperClusterTrackAtVtx()) > 0.01 ) continue;
            if( TMath::Abs(el->scSigmaIEtaIEta()) > 0.03 ) continue;
        }
        
                
        float pfRelIsoLocal = pfRelIso(el, Rho) ;
        if( pfRelIsoLocal > v_electron_reliso ) continue;
        
        std::vector< Leptons >::iterator iter;
        for( iter = LepCand.begin() ; iter < LepCand.end(); iter++ ){
            
            TString Type         = iter->type;
            TString Flavor       = iter->flavor;
            int Charge           = iter->charge;
            TLorentzVector TLVec = iter->tlvector;
            
            if( Flavor == "m" ) continue;
            if( signalLepOnly && Type == "sideband" ) continue;
            
            if( el->charge() == Charge  ) continue;
            
            TLorentzVector P1; P1.SetPtEtaPhiE( el->pt(), el->eta(), el->phi(), el->energy() );
            
            float Mass = ( P1 + TLVec ).M();
            
            if( CutName == "Z" ){
                if( Mass < 106. && Mass > 76. )  {
                    ZMass = Mass;
                    myPair.first  = iter->index;
                    myPair.second = iter->flavor;
                    VecOfPair.push_back( myPair );
                }
            }
            else if( CutName == "gammastar" ){
                if( Mass < 12. )  {
                    ZMass = Mass;
                    myPair.first  = iter->index;
                    myPair.second = iter->flavor;
                    VecOfPair.push_back( myPair );
                }
            }
        }
    }
    
    std::vector< Leptons >::iterator iter = LepCand.begin();
    while ( iter != LepCand.end() ) {
        
        TString Type    = iter->type;
        TString Flavor  = iter->flavor;
        unsigned int Index       = iter->index;
        
        //       if( Flavor == "m" ) continue;
        //       if( signalLepOnly && Type == "sideband" ) continue;
        
        bool remove = false;
        std::vector< pair<unsigned int,TString >>::iterator tmp;
        for( tmp = VecOfPair.begin(); tmp < VecOfPair.end(); tmp++ ){
            
            if( (*tmp).first == Index &&  (*tmp).second == Flavor ) remove = true;
        }
        
        if( remove ) iter = LepCand.erase( iter );
        else
            ++iter;
        
    }
}




void tools::removeZCandForMM(const std::vector<pat::Muon>  & thePatMuons,
                             reco::Vertex::Point PV,
                             std::vector< Leptons > & LepCand,
                             bool signalLepOnly,
                             const char* cutName ){
    
    
    //bool bool_pfIsolation = true;
    float v_muon_pt = 0;
    float v_muon_eta = 2.4;
    //  float v_muon_dz  = 999.;
    float v_muon_iso = 0.2;
    //  bool ZMass = false;
    double ZMass = 0;
    
    TString CutName = cutName;
    
    
    if( CutName == "gammastar" ) v_muon_pt = 5;
    else if( CutName == "Z"    ) v_muon_pt = 10;
    
    std::pair<unsigned int , TString  > myPair;
    std::vector< pair< unsigned int,TString >> VecOfPair;
    
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons.begin() ; mu != thePatMuons.end() ; mu++ ) {
        
        if ( mu->pt()  < v_muon_pt ) continue;
        if ( TMath::Abs( mu->eta() ) > v_muon_eta ) continue;
        if ( mu->isTrackerMuon() || mu->isGlobalMuon() ) {
            if ( !mu->isPFMuon() ) continue;
            
            //    TVector3 momPerp(0,0,0);
            //    momPerp.SetPtEtaPhi(mu->pt(),mu->eta(),mu->phi());
            //    TVector3 posPerp(mu->vx()-PV.x(), mu->vy() - PV.y(), 0);
            //    float dzcorr = mu->vz() - PV.z() - posPerp.Dot(momPerp)/mu->pt() * (mu->pz()/mu->pt());
            //   if(TMath::Abs(dzcorr ) > v_muon_dz  ) continue;
            
            // PF Isolation
            float chargedHadronIso = mu->pfIsolationR03().sumChargedHadronPt;
            float neutralHadronIso = mu->pfIsolationR03().sumNeutralHadronEt;
            float photonIso = mu->pfIsolationR03().sumPhotonEt;
            float beta = mu->pfIsolationR03().sumPUPt;
            float pfRelIso  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/mu->pt() ;
            
            // Det Isolation
            float detRelIso = (mu->isolationR03().emEt + mu->isolationR03().hadEt + mu->isolationR03().sumPt) / mu->pt() ;
            
            if( bool_pfIsolation ){
                
                if( pfRelIso > v_muon_iso ) continue;
            }
            else{
                
                if( detRelIso > v_muon_iso ) continue;
            }
            
            
            std::vector< Leptons >::iterator iter;
            for( iter = LepCand.begin() ; iter < LepCand.end(); iter++ ){
                
                TString Type         = iter->type;
                TString Flavor       = iter->flavor;
                int Charge           = iter->charge;
                TLorentzVector TLVec = iter->tlvector;
                
                if( signalLepOnly && Type == "sideband" ) continue;
                
                if( Flavor == "e" ) continue;
                if( mu->charge() == Charge  ) continue;
                
                TLorentzVector P1; P1.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
                
                float Mass = ( P1 + TLVec ).M();
                
                
                if( CutName == "Z" ){
                    if( Mass < 106. && Mass > 76. )  {
                        ZMass = Mass;
                        myPair.first  = iter->index;
                        myPair.second = iter->flavor;
                        VecOfPair.push_back( myPair );
                    }
                }
                else if( CutName == "gammastar" ){
                    if( Mass < 12. )  {
                        ZMass = Mass;
                        myPair.first  = iter->index;
                        myPair.second = iter->flavor;
                        VecOfPair.push_back( myPair );
                    }
                }
            }
        }
    }
    
    std::vector< Leptons >::iterator iter = LepCand.begin();
    while ( iter != LepCand.end() ) {
        
        TString Type    = iter->type;
        TString Flavor  = iter->flavor;
        unsigned int Index       = iter->index;
        
        //   if( Flavor == "e" ) continue;
        //   if( signalLepOnly && Type == "sideband" ) continue;
        
        bool remove = false;
        std::vector< pair< unsigned int,TString >>::iterator tmp;
        for( tmp = VecOfPair.begin(); tmp < VecOfPair.end(); tmp++ ){
            
            if( (*tmp).first == Index &&  (*tmp).second == Flavor ) remove = true;
        }
        
        if( remove ) iter = LepCand.erase( iter );
        else
            ++iter;
        
    }
    
}

*/

double tools::JER (double eta) {

    double feta = fabs(eta);
    double scf = 1;
    if (feta < 0.5)
        scf = 1.052;
    else if (feta < 1.1)
        scf = 1.057;
    else if (feta < 1.7)
        scf = 1.096;
    else if (feta < 2.3)
        scf = 1.134;
    else scf = 1.288;
    
    return scf;
}

float tools::quadsum(float a, float b) {
    return sqrt(a*a+b*b);
}

//______________________________________________________________________________
#define JERSUMMER11
float tools::smear_pt_res(float pt, float genpt, float eta)
{
    eta = fabs(eta);
    if (genpt>15. && (fabs(pt - genpt) / pt)<0.5) {  // limit the effect to the core
        double res    = 1.0;
        //double resErr = 0.0;
#ifdef JERSUMMER11
        // from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        if (eta <= 0.5) {
            res    = 1.052;
            //resErr = quadsum(0.012, 0.062);
        } else if (0.5 < eta && eta <= 1.1) {
            res    = 1.057;
            //resErr = quadsum(0.012, 0.062);
        } else if (1.1 < eta && eta <= 1.7) {
            res    = 1.096;
            //resErr = quadsum(0.017, 0.063);
        } else if (1.7 < eta && eta <= 2.3) {
            res    = 1.134;
            //resErr = quadsum(0.035, 0.087);
        } else {
            res    = 1.288;
            //resErr = quadsum(0.127, 0.155);
        }
#else
        // from VHbb analysis
        if (eta <= 1.1) {
            res    = 1.05;
            //resErr = 0.05;
        } else if (1.1 < eta && eta <= 2.5) {
            res    = 1.10;
            //resErr = 0.10;
        } else {
            res    = 1.30;
            //resErr = 0.20;
        }
#endif
        float deltapt = (pt - genpt) * res;
        return TMath::Max(float(0.), genpt + deltapt);
    }
    return pt;
}

bool  tools::srSSbID(int nJets, int nbJets, double met, double Ht, int id) {

    if ((nbJets >= 2 ) && (nJets >= 2) && (Ht > 80) && (met > 30) && (id == 1)) return true;
    if ((nbJets >= 2 ) && (nJets >= 2) && (Ht > 80) && (met > 30) && (id == 2)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 200) && (met > 120) && (id == 3)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 200) && (met > 50) && (id == 4)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 320) && (met > 50) && (id == 5)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 320) && (met > 120) && (id == 6)) return true;
    if ((nbJets >= 3 ) && (nJets >= 3) && (Ht > 200) && (met > 50) && (id == 7)) return true;
    if ((nbJets >= 2 ) && (nJets >= 4) && (Ht > 320) && (met > 0) && (id == 8)) return true;

    return false;
}

bool  tools::srSSbID28(int nJets, int nbJets, double met, double Ht, int id) {
    
    if ((id == 0) || (id == 9) || (id == 10) || (id == 19) || (id == 20)) return false;
    
    if (nJets < 2) return false;
    if (met < 50) return false;
    if (Ht < 200) return false;
    
    if (!((id < 10) || ((nbJets == 1) && (id/10 == 1)) || ((nbJets > 1) && (id/10 == 2)))) return false;
    if (!(((met < 120) && (id%10 < 5)) || ((met > 120) && (id%10 >=5)))) return false;
    if (!(((Ht < 400) && (id%2 == 1)) || ((Ht > 400) && (id%2 != 1)))) return false;
    if (!(((nJets < 4) && ((((id-1)%10)%4)/2 == 0)) || ((nJets >= 4) && ((((id-1)%10)%4)/2 != 0)))) return false;
    
    return true;

}

//______________________________________________________________________________
#include "TLorentzVector.h"
double tools::evalEt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Et();
}

double tools::evalMt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Mt();
}

std::vector<double> tools::RegressionVars(const pat::Jet *jet, float genpt, const pat::Muon* mu) {
    std::vector<double> vars;
    
    vars.push_back(smear_pt_res((jet->correctedP4("Uncorrected")).Pt(), genpt, jet->eta()));
    vars.push_back(jet->pt());
    vars.push_back(evalEt(jet->pt(), jet->eta(), jet->phi(), jet->energy()));
    vars.push_back(evalMt(jet->pt(), jet->eta(), jet->phi(), jet->energy()));
    
    double hJet_ptLeadTrack = 0;
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
    }
    vars.push_back(hJet_ptLeadTrack);
    
    double hJet_vtx3dL = 0;
    double hJet_vtx3deL = 0;
    double hJet_vtxMass = 0;
    double hJet_vtxPt = 0;
    
    const reco::SecondaryVertexTagInfo* scdVtx = jet->tagInfoSecondaryVertex("secondaryVertex");
    
    if (scdVtx) {
        //std::cout<<"Vertetx info: "<<scdVtx->nVertices()<<std::endl;
        if (scdVtx->nVertices()) {
            const reco::Vertex &sv1 = scdVtx->secondaryVertex(0);
            if (!sv1.isFake()) {
                Measurement1D distance1 = scdVtx->flightDistance(0, true);
                hJet_vtx3dL = distance1.value();
                hJet_vtx3deL = distance1.error();
                
                math::XYZTLorentzVectorD p4vtx = sv1.p4();
                hJet_vtxMass = p4vtx.M();
                hJet_vtxPt = p4vtx.Pt();
            }
        }
    } //else
      //  std::cout<<"no info"<<std::endl;
    
    vars.push_back(TMath::Max(0.,hJet_vtx3dL));
    vars.push_back(TMath::Max(0.,hJet_vtx3deL));
    vars.push_back(TMath::Max(0.,hJet_vtxMass));
    vars.push_back(TMath::Max(0.,hJet_vtxPt));
    
    vars.push_back(jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction());
    vars.push_back(jet->getPFConstituents().size());
    
    vars.push_back(0.); //JECUnc in the main file
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    vars.push_back(pLep.Perp(pJet.Vect()));
    vars.push_back(mu->pt());
    vars.push_back(pLep.DeltaR(pJet));
    
    /*values[0] = "breg_rawptJER := smear_pt_res(hJet_ptRaw, hJet_genPt, hJet_eta)";
    values[1] = "breg_pt := hJet_pt";
    values[2] = "breg_et := evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)";
    values[3] = "breg_mt := evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)";
    values[4] = "breg_leadtrackpt := hJet_ptLeadTrack";
    values[5] = "breg_vtx3dL := max(0,hJet_vtx3dL)";
    values[6] = "breg_vtx3deL := max(0,hJet_vtx3deL)";
    values[7] = "breg_vtxMass := max(0,hJet_vtxMass)";
    values[8] = "breg_vtxPt := max(0,hJet_vtxPt)";
    values[9] = "breg_cef := hJet_cef";
    values[10] = "breg_ntot := hJet_nconstituents";
    values[11] = "breg_eJEC := hJet_JECUnc";
    values[12] = "breg_softlepptrel := max(0,hJet_SoftLeptptRel*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[13] = "breg_softleppt := max(0,hJet_SoftLeptPt*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[14] = "breg_softlepdR := max(0,hJet_SoftLeptdR*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[15] = "breg_evt_rho25 := rho25";*/
    
    return vars;
    
}

