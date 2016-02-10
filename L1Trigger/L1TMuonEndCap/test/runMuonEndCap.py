# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
import commands


inputFile=""
outputFile=""

def getVal(arg):
    i=0
    while i < len(arg) and arg[i] != "=": i+=1
    return arg[i+1:]
	
## loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i] 
    if "isMC" in sys.argv[i] :
        isMC=False if getVal(sys.argv[i]) == "False" else True
    elif "isMSUGRA" in sys.argv[i]:
        isMSUGRA=True
    elif "isSMS" in sys.argv[i]:
        isSMS=True
    elif "skim" in sys.argv[i]:
        skim=getVal(sys.argv[i])
    elif "output" in sys.argv[i]:
        outputFile=getVal(sys.argv[i])
    elif "input" in sys.argv[i] :
        inputFile=getVal(sys.argv[i])


process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))



process.source = cms.Source('PoolSource',
 fileNames = cms.untracked.vstring(
			#'/cms/data/store/data/Run2015D/ZeroBias/RAW/v1/000/259/721/00000/86FD2518-7A78-E511-B5B7-02163E0145C1.root',
			#'/cms/data/store/data/Run2015D/ZeroBias/RAW/v1/000/259/721/00000/0A2BA199-5878-E511-B77E-02163E01190D.root',
			#'/cms/data/store/mc/Fall13dr/MuPlus_Pt-1to150_PositiveEndcap-gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/20000/04643A9D-F382-E311-9D3C-0017A4770C34.root'
								)
	                    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

# PostLS1 geometry used
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
############################
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


####Event Setup Producer
process.load('L1Trigger.L1TMuonEndCap.fakeMuonEndCapParams_cfi')
process.esProd = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet = cms.VPSet(
      cms.PSet(record = cms.string('L1TMuonEndcapParamsRcd'),
               data = cms.vstring('L1TMuonEndcapParams'))
                   ),
   verbose = cms.untracked.bool(True)
)

process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load('Configuration.StandardSequences.L1Emulator_cff')

#process.simCscTriggerPrimitiveDigis.CSCWireDigiProducer = cms.InputTag( 'MuonCSCDigis', 'MuonCSCWireDigi' )

process.L1ReEmulSeq = cms.Sequence(#process.L1Emulator
                                    process.ecalDigis
                                   + process.hcalDigis
                                   + process.gtDigis
                                   + process.gtEvmDigis
                                   + process.csctfDigis
                                   + process.dttfDigis
								   + process.RawToDigi
								   # + process.simL1Emulator
                                   )

####OMTF Emulator
process.load('L1Trigger.L1TMuonEndCap.simMuonEndCapDigis_cfi')

process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

process.L1TMuonSeq = cms.Sequence( process.esProd          
                                   + process.simEmtfDigis 
#                                   + process.dumpED
#                                   + process.dumpES
)

process.L1TMuonPath = cms.Path(process.L1ReEmulSeq+process.L1TMuonSeq)

for i in inputFile.split(","):
   print "Adding: ", i
   process.source.fileNames.append(i)


process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("FUCK.root") )

#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/mrcarver/CondFormatsCrap/CMSSW_8_0_0_pre5/src/L1Trigger/L1TMuonEndCap/test/JSONfile.txt').getVLuminosityBlockRange() 


outCommands = cms.untracked.vstring('keep *')

process.out = cms.OutputModule("PoolOutputModule", 
   fileName = cms.untracked.string("l1tomtf_superprimitives1.root"),
   outputCommands = outCommands,
   
)

#process.output_step = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.L1TMuonPath)
#process.schedule.extend([process.output_step])

from SLHCUpgradeSimulations.Configuration.muonCustoms import customise_csc_PostLS1
process = customise_csc_PostLS1(process)
