import FWCore.ParameterSet.Config as cms

process = cms.Process("L1MicroGMTEmulator")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')


process.source = cms.Source("PoolSource",
    # ttbar:
     fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_7_2_2_patch1/RelValTTbarLepton_13/GEN-SIM-RECO/MCRUN2_72_V3_71XGENSIM-v2/00000/18295FD0-1B74-E411-893D-0025905A6094.root",
        "/store/relval/CMSSW_7_2_2_patch1/RelValTTbarLepton_13/GEN-SIM-RECO/MCRUN2_72_V3_71XGENSIM-v2/00000/92B21B80-1074-E411-8367-0025905B85A2.root",
        "/store/relval/CMSSW_7_2_2_patch1/RelValTTbarLepton_13/GEN-SIM-RECO/MCRUN2_72_V3_71XGENSIM-v2/00000/D25833CC-1B74-E411-AC75-0025905B8592.root",
     ),
    #fileNames = cms.untracked.vstring (
    #    "file:/afs/cern.ch/work/j/jlingema/private/scratch0/L1TDev/CMSSW_7_3_0_pre1/src/L1Trigger/J_psi_to_mumu_high_boost.root"
    #),
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS1', '')

process.uGMTInputProducer = cms.EDProducer("l1t::uGMTInputProducerFromGen",
)



# process.load("L1Trigger.L1TGlobalMuon.microgmtemulator_cfi")


# process.microGMTEmulator.barrelTFInput = cms.InputTag("uGMTInputProducer", "BarrelTFMuons")
# process.microGMTEmulator.overlapTFInput = cms.InputTag("uGMTInputProducer", "OverlapTFMuons")
# process.microGMTEmulator.forwardTFInput = cms.InputTag("uGMTInputProducer", "ForwardTFMuons")
# process.microGMTEmulator.triggerTowerInput = cms.InputTag("uGMTInputProducer", "TriggerTowerSums")

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'keep *_*_*_L1MicroGMTEmulator',
    ),
    fileName = cms.untracked.string('ttbar_large_sample.root')
)


#process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
    process.uGMTInputProducer
    #*process.microGMTEmulator
    )

process.e = cms.EndPath(process.out)
