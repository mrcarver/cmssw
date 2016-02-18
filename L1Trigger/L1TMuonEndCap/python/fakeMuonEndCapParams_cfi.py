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
   configXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/emtf_config/EMTF_XML_Config.xml"),
)




