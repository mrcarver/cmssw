import sys
import FWCore.ParameterSet.Config as cms

process = cms.Process('TEXTDUMP')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1))
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")



isMC     = bool(True)
isMSUGRA = bool(False)
isSMS    = bool(False)
useCrab  = bool(False)
doSusyTopProjection = bool(False)
inputFile=""
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

#process.load("CommonTools.ParticleFlow.pfIsolation_cfg")

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.GlobalTag.globaltag="POSTLS170_V6::All" # tag for 50ns good alca
#process.GlobalTag.globaltag="PHYS14_25_V1::All" #phys14 25ns tag

#process.L1TMuonText = cms.EDProducer(
process.Eff = cms.EDAnalyzer(
    'SSOpti',
    doGen = cms.untracked.bool(True),
	sample = cms.uint32(5),#1-T1tttt1200::2-T1tttt1500::3-T5ttttD
    genSrc = cms.untracked.InputTag("genParticles"),
    lutParam = cms.PSet(
    isBeamStartConf = cms.untracked.bool(True),
    ReadPtLUT = cms.bool(False),
    PtMethod = cms.untracked.uint32(32)
    )
)


process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" # checks for fake PVs automatically
	, filterParams =cms.PSet(
	 		  	  minNdof = cms.double( 4. )
				, maxZ    = cms.double( 24. )
				, maxRho  = cms.double( 2. ) )		
				, filter       = cms.bool( False ) # use only as producer
				, src          = cms.InputTag( 'offlineSlimmedPrimaryVertices' )
		
)

#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso, setupPFPhotonIso

process.content = cms.EDAnalyzer("EventContentAnalyzer")


#fileOutName = "Tester0.root"
#fileOutName = "AllSamples/ttbar.root"
fileOutName = "T5_1300_C_Out.root"
#fileOutName = "AllSamples/T1_1200_Out.root"

#for i in inputFile.split(","):
#    print "Adding: ", i
#    process.source.fileNames.append(i)

process.source = cms.Source(
    'PoolSource',
	duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
	
	
		#TTbar :::: Sample 7
		#'/cms/data/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/007B37D4-8B70-E411-BC2D-0025905A6066.root',
		#'/cms/data/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/06843FC5-8370-E411-9B8C-0025905A60AA.root',
		#'/cms/data/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0A867F71-8C70-E411-9CC9-0025905A48D6.root',
		#'/cms/data/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0C1D0A70-8870-E411-BAB1-0025905A612C.root',
	
		#T1tttt Gluino 1200 LSP 800 ::: Sample 1
		#'/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0CD15D7F-4E6B-E411-AEB4-002590DB9216.root',
		#'/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/4CB22889-4E6B-E411-B0AA-002481E0D144.root',
		#'/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/261DD410-5E6B-E411-933B-002590D57E3C.root',
		#'/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/D8C225F4-5D6B-E411-A778-0025901D42C0.root',#problem with this file...
	
		#T1tttt Gluino 1500 LSP 100 ::: Sample 2
		#'/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/4836AA58-5A6B-E411-8865-20CF305B053E.root',
		#'/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/829D372D-7F6B-E411-81B1-0025907B5048.root',
		#'/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/C43D68C5-8D6B-E411-B030-0025907750A0.root',
		
		#T5tttt Gluino 800 Stop 300 Chargino 285 Chi 280 2-3 Body decay ::: Sample 3
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.00.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.01.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.02.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.03.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.04.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.05.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.06.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.07.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.08.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.09.root',
	   ##'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.10.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.11.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.12.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.13.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.14.root',
		#'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo800_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.15.root',
		
		#T5tttt Gluino 800 Stop 300 Chi 280 4 Body offshell decay ::: Sample 4
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.00.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.01.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.02.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.03.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.04.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.05.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.06.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.07.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.08.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.09.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.10.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.11.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.12.root',
	   #'
	   #'
	   #'

	   
	   
		#T5tttt Gluino 1300 Stop 300 Chargino 285 Chi 280 2-3 Body decay ::: Sample 5
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.00.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.01.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.02.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.03.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.04.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.05.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.06.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.07.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.08.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.09.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.10.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.11.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.12.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.13.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.14.root',
		'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec/T1tttt_2J_mGo1300_mStop300_mCh285_mChi280_pythia8-23bodydec.MINIAODSIM.15.root',
		
		#T5tttt Gluino 1300 Stop 300 Chi 280 4 Body offshell decay ::: Sample 6
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.00.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.01.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.02.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.03.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.04.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.05.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.06.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.07.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.08.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.09.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.10.root',
	   #'/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo1300_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.11.root',
	   #'
	   #'
	   #'
	   #'
	  
		
	
	
	)
    
    #eventsToProcess = cms.untracked.VEventRange('1:590','1:602','1:640','1:906','1:1029','1:1318','1:1622','1:1765','1:1817','1:2205','1:2206','1:2494','1:2541')
    #eventsToProcess = cms.untracked.VEventRange('1:1925-1:1927','1:2769-1:2771','1:3669-1:3671'),
    #eventsToSkip= cms.untracked.VEventRange('1:896')
    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(
	fileOutName
))

process.EffSequence = cms.Sequence(process.goodOfflinePrimaryVertices*process.Eff)
#process.L1TMuonSequence = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence + process.phoIsoSequence)

process.EffPath = cms.Path(process.EffSequence)

process.schedule = cms.Schedule(process.EffPath)

