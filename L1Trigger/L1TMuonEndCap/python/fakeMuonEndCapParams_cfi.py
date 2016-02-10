import FWCore.ParameterSet.Config as cms


# NOT YET IMPLEMENTED ...

# This file is included in the expected place, so standard sequences will find it as expected.

emtfParamsSource = cms.ESSource(
   "EmptyESSource",
   recordName = cms.string('L1TMuonEndcapParamsRcd'),
   iovIsRunNotTime = cms.bool(True),
   firstValid = cms.vuint32(1)
)

###EMTF ESProducer. Fills CondFormats from XML files.
emtfParams = cms.ESProducer(
   "L1TMuonEndcapParamsESProducer",
   PT_assignment_version = cms.uint32(1),
   NumPhiBits = cms.int32(7)
)




