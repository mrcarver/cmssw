import FWCore.ParameterSet.Config as cms

process = cms.Process('L1EMTF')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1))
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.L1TMuonTrackFinderEndCap.L1TMuonTriggerPrimitiveProducer_cfi')
#process.load('Configuration.Geometry.GeometryExtendedPostLS1Reco_cff')
#process.load('Configuration.Geometry.GeometryExtendedPostLS1_cff')

process.load('Configuration.Geometry.GeometryExtended2015devReco_cff')
process.load('Configuration.Geometry.GeometryExtended2015dev_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.L1TMuonEndcapTrackFinder = cms.EDProducer(
    'L1TMuonUpgradedTrackFinder',
    
	
    primitiveSrcs = cms.VInputTag(
    cms.InputTag('L1TMuonTriggerPrimitives','CSC'),
    ),
   
)


process.L1TMuonTriggerPrimitives.CSC.src = cms.InputTag('csctfDigis')

del process.L1TMuonTriggerPrimitives.RPC  
del process.L1TMuonTriggerPrimitives.DT 

process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")



process.content = cms.EDAnalyzer("EventContentAnalyzer")

infile = [

		'/store/data/Run2015C/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/254/380/00000/7A04A4B0-5E46-E511-9A5C-02163E013891.root',
		]


process.source = cms.Source(
    'PoolSource',
    fileNames = cms.untracked.vstring(infile)
    
    )



outCommands = cms.untracked.vstring('keep *')

process.FEVTDEBUGoutput = cms.OutputModule(
    "PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = outCommands,
    fileName = cms.untracked.string('Emulator_EDM_Out.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    )
)

process.L1TMuonSequence = cms.Sequence( process.csctfDigis +
										process.L1TMuonTriggerPrimitives + 
									process.L1TMuonEndcapTrackFinder)

process.L1TMuonPath = cms.Path(process.L1TMuonSequence)

process.outPath = cms.EndPath(process.FEVTDEBUGoutput)

process.schedule = cms.Schedule(process.L1TMuonPath,process.outPath)

