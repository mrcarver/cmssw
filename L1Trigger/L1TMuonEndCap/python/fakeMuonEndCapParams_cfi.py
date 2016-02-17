import FWCore.ParameterSet.Config as cms

emtfParamsSource = cms.ESSource(
	"EmptyESSource",
	recordName = cms.string('L1TMuonEndcapParamsRcd'),
	iovIsRunNotTime = cms.bool(True),
	firstValid = cms.vuint32(1)
)

##EMTF ESProducer. Fills CondFormats from XML files.
emtfParams = cms.ESProducer(
	"L1TMuonEndcapParamsESProducer",
   PT_assignment_version = cms.uint32(1),
   St1PhiMatchWindow = cms.int32(15),
   St2PhiMatchWindow = cms.int32(15),
   St3PhiMatchWindow = cms.int32(7),
   St4PhiMatchWindow = cms.int32(7),
   XmlPtDir = cms.string('ModeVariables_v1_dTheta')
)




