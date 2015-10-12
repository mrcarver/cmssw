import sys
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts

skim="NOSkim"
isMC     = bool(True)
isMSUGRA = bool(False)
isSMS    = bool(False)
useCrab  = bool(False)
doSusyTopProjection = bool(False)
inputFile=""
outputFile="ssb"
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

        
if skim=="" : print "WARNING: No Skim Conditions have been provided \n"

#process = cms.Process("FakeElectrons")
process = cms.Process("pippo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'WARNING' # Options: INFO, WARNING, ERROR
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

if isMC:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    #process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco
	process.GlobalTag.globaltag="POSTLS170_V6::All" # tag for 50ns good alca
else:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    #process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco
	process.GlobalTag.globaltag="POSTLS170_V6::All" # tag for 50ns good alca

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),                    
                             fileNames = cms.untracked.vstring(),      
                             )


#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

#if isData:
#    process.source.fileNames = cms.untracked.vstring( 
#    #'file:../../../CMSSW_5_2_3_PatTuple_test_PFiso_El_newMETCor.root'
#    #'file:DoubleElectron_Run2012B_TwoLepSkim_SSDiLepSkim_193752_194076.root'
#    )
#    runOnData(process)
#       #'/store/user/ddidar/CRAB_output/DoubleEle/CMSSW_5_2_3_PATtupple_2Lep2JetSkim_83_1_sOE.root'
#       #                                               )

#else:
#    process.source.fileNames = cms.untracked.vstring( )
#    runOnData(process)
#              #'file:/raid/raid10/ddidar/CMSSW_5_2_3/src/CMSSW_5_2_3_PatTuple.root'
#	      #'file:/raid/raid10/ddidar/CMSSW_5_2_3/src/CMSSW_5_2_3_TTJets_2Lep2JetSkim.root'
#              #'/store/user/ddidar/CRAB_output/MC/TTJets/CMSSW_5_2_3_TTJets_2Lep2JetSkim_103_1_7vU.root'
#	      #)
   

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )
## input files
for i in inputFile.split(","):
	print "Adding: ", i
	process.source.fileNames.append(i)
	
#process.source = cms.Source(
#	'PoolSource',
#	fileNames = cms.untracked.vstring('/cms/data/store/user/lshchuts/MCSignalSamples/CSA14/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec/T1tttt_2J_mGo800_mStop300_mChi280_pythia8-4bodydec.MINIAODSIM.00.root')
#		#'/scratch/lfs/lesya/TTW/miniAOD-prod_PAT_1.root')
#)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#FakeElectrons
process.FakeElectrons = cms.EDAnalyzer("SSb13",#"FakeMuonsYa",
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       #TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
                                       JetLabel = cms.InputTag("slimmedJets"),
                                       BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                       HLTResultsLabel = cms.InputTag("TriggerResults"),
                                       METLabel = cms.InputTag("slimmedMETs"),
                                       #METFilter = cms.InputTag("TriggerResults::PAT"),
                                       qualityCuts = PFTauQualityCuts,
                                       SampleLabel = cms.untracked.string("ElectronsMC") # defines a piece of code to run; helps to avoid code recompilation
                                       )

process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" # checks for fake PVs automatically
	, filterParams =cms.PSet(
	 		  	  minNdof = cms.double( 4. )
				, maxZ    = cms.double( 24. )
				, maxRho  = cms.double( 2. ) )		
				, filter       = cms.bool( False ) # use only as producer
				, src          = cms.InputTag( 'offlineSlimmedPrimaryVertices' )
)


process.electronMatch = cms.EDProducer("MCMatcher",
                                       src         = cms.InputTag("slimmedElectrons"),          # RECO objects to match
                                       matched     = cms.InputTag("prunedGenParticles"),          # mc-truth particle collection
                                       mcPdgId     = cms.vint32(),                          # one or more PDG ID (15 = tau); absolute values (see below)
                                       checkCharge = cms.bool(False),                       # True = require RECO and MC objects to have the same charge
                                       mcStatus    = cms.vint32(),                          # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                       # NOTE that Taus can only be status 3 or 2, never 1!
                                       maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed 
                                       maxDPtRel   = cms.double(10.),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                                       resolveAmbiguities    = cms.bool(False),                  # Forbid two RECO objects to match to the same GEN object
                                       resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
                                       )

process.muonMatch = cms.EDProducer("MCMatcher",
                                       src         = cms.InputTag("slimmedMuons"),              # RECO objects to match
                                       matched     = cms.InputTag("prunedGenParticles"),          # mc-truth particle collection
                                       mcPdgId     = cms.vint32(),                          # one or more PDG ID (15 = tau); absolute values (see below)
                                       checkCharge = cms.bool(False),                       # True = require RECO and MC objects to have the same charge
                                       mcStatus    = cms.vint32(),                          # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                       # NOTE that Taus can only be status 3 or 2, never 1!
                                       maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed
                                       maxDPtRel   = cms.double(10.),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                                       resolveAmbiguities    = cms.bool(False),                  # Forbid two RECO objects to match to the same GEN object
                                       resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
                                       )

## reco-generator(parton) matching for jets
process.jetPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, -1, -2, -3, -4, -5, 21),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(3.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(False),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.jetPartonMatch2 = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, -1, -2, -3, -4, -5, 21),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(2),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(3.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(False),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.electronPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                             src = cms.InputTag("slimmedElectrons"),      # RECO objects to match
                                             matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                             mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                             mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                             checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                             maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                             maxDPtRel = cms.double(300.0),                # Minimum deltaPt/Pt for the match
                                             resolveAmbiguities = cms.bool(False),        # Forbid two RECO objects to match to the same GEN object
                                             resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                             )

process.muonPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                         src = cms.InputTag("slimmedMuons"),      # RECO objects to match
                                         matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                         mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                         mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                         checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                         maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                         maxDPtRel = cms.double(300.0),                # Minimum deltaPt/Pt for the match
                                         resolveAmbiguities = cms.bool(False),        # Forbid two RECO objects to match to the same GEN object
                                         resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                         )
# A set of quality cuts used for the PFTaus.  Note that the quality cuts are
# different for the signal and isolation regions.  (Currently, only in Nhits)

if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
   # *process.tauMatch
#    *process.muonMatch
#    *process.electronMatch
#    *process.jetPartonMatch
#    *process.jetPartonMatch2
#    *process.electronPartonMatch
#    *process.muonPartonMatch
    *process.FakeElectrons
    )
else: 
    import PhysicsTools.PythonAnalysis.LumiList as LumiList
#    process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/hpc/lesya/CMSSW_5_3_6_patch1/src/SUSYAnalyzer/PatAnalyzer/test/JSON/Cert_190456-207898_8TeV_PromptReco_Collisions12_JSON.txt').getVLuminosityBlockRange()
    process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test/fullJSON_SUSYFakes.json').getVLuminosityBlockRange()
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
       *process.FakeElectrons
    )
    
