import FWCore.ParameterSet.Config as cms

###EMTF emulator configuration
simEmtfDigis = cms.EDProducer("L1TMuonEndCapTrackProducer",
                              #CSCInput = cms.InputTag('simCscTriggerPrimitiveDigis','MPCSORTED')
							  CSCInput = cms.InputTag('csctfDigis',''),
							  GMTInput = cms.InputTag('gtDigis','','L1TMuonEmulation'),
							  LegInput = cms.InputTag('csctfDigis','','L1TMuonEmulation')
                              )


