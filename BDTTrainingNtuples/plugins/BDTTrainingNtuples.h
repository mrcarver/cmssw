#ifndef _BDTTrainingNtuples_h_
#define _BDTTrainingNtuples_h_
//asd
#include <memory>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "L1Trigger/L1TMuon/interface/deprecate/SubsystemCollectorFactory.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "L1Trigger/L1TMuon/interface/deprecate/GeometryTranslator.h"
#include "L1Trigger/L1TMuon/interface/deprecate/MuonTriggerPrimitive.h"

#include "L1Trigger/L1TMuonEndCap/interface/MuonInternalTrack.h"
#include "L1Trigger/L1TMuonEndCap/interface/MuonInternalTrackFwd.h"


#include "L1Trigger/L1TMuonEndCap/interface/PhiMemoryImage.h"
#include "L1Trigger/L1TMuonEndCap/interface/EmulatorClasses.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCTFPtLUT.h"
#include "L1Trigger/CSCTrackFinder/interface/CSCSectorReceiverLUT.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"

#include "L1Trigger/L1TMuonEndCap/interface/PrimitiveConverter.h"
#include "L1Trigger/L1TMuonEndCap/interface/PtAssignment.h"

typedef edm::ParameterSet PSet;

class BDTTrainingNtuples : public edm::EDProducer {
public:
  BDTTrainingNtuples(const PSet&);
  BDTTrainingNtuples() {}
	
 void produce(edm::Event&, const edm::EventSetup&); 
  void beginJob();
  void endJob(); 

 private:
  
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> inputTokenCSC;
  int bxShiftCSC = 0;
  PrimitiveConverter primConv_;
  l1t::EmtfPtAssignment ptAssignment_;
  
};


#endif
