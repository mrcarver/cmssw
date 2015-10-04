import sys
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts


skim     ="NOSkim"
isMC     = bool(True)
isMSUGRA = bool(False)
isSMS    = bool(False)
useCrab  = bool(False)
doSusyTopProjection = bool(False)
inputFile=""
outputFile="results/FakeMuons_Interactive.root"
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

if isMC:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco
else:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),                    
                             fileNames = cms.untracked.vstring(),      
                             )

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

if isMC:
    process.source.fileNames = cms.untracked.vstring(               
#    '/store/user/ddidar/MC_JEC53X/WZJetsTo3LNu/PatTuple_WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_NOSkim_93_1_zD0.root'
#    '/cms/data/store/user/ddidar/MC_JEC53X/TTTo2L2Nu2B/PatTuple_TTTo2L2Nu2B_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_NOSkim_89_0_e88.root'
#    '/cms/data/store/user/ddidar/MC_JEC53X/BJetsNew/BJets_HT-100To250/PatTuple_BJets_HT-100To250_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_OneLepOneBJetSkim_105_1_Fsk.root'
#    '/cms/data/store/user/ddidar/MC_JEC53X/BJets_HT-100To250/PatTuple_BJets_HT-100To250_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_OneLepOneBJetSkim_105_1_Vhy.root'
#    '/cms/data/store/user/peveraer/SbottomTop/PatTuple_SbottomTop_100_1_3Iq.root'
#    '/cms/data/store/user/skhi/SignalSamples/mSbottom340-480_mChargino140-300_GoodSample_v4/PatTuple_SbottomTopv6_macneill-SbottomTopv6-ebe1b6443ab75fb2cd06c2581e9a7621_USER_NOSkim_10_1_RsK.root'
#    '/cms/data/store/user/ddidar/MC_JEC53X/WJetsToLNu_OneLepOneBJetSkim/PatTuple_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_OneLepOneBJetSkim_997_1_goM.root'
#     '/scratch/lfs/lesya/FakeSync.root'
#     '/scratch/lfs/lesya/FakeSync_jets_test_correctGT.root'
#     '/cms/data/store/user/ddidar/MC_JEC53X/WJetsToLNu_NOSkim//PatTuple_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_NOSkim_1045_3_3cj.root'
#      '/cms/data/store/user/ddidar/MC_JEC53X/QCD_Pt_20_MuEnrichedPt_15_CorrectGT/PatTuple_QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3_AODSIM_NOSkim_970_1_GY6.root'
#     '/scratch/lfs/lesya/FakeSync_WJets.root'
#     '/scratch/lfs/lesya/FakeSync_TTWJets.root'
#     '/cms/data/store/user/ddidar/MC_JEC53X/WJetsToLNu_NOSkim/PatTuple_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_NOSkim_983_3_gek.root'
     '/cms/data/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/003E832C-8AFC-E311-B7AA-002590596490.root'
    )

else:
   import PhysicsTools.PythonAnalysis.LumiList as LumiList
#   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test/fullJSON_SUSYFakes.json').getVLuminosityBlockRange()

   process.source.fileNames = cms.untracked.vstring( 
#    #'file:../../../CMSSW_5_2_3_PatTuple_test_PFiso_El_newMETCor.root'
#    'file:/raid/raid10/ddidar/CMSSW_5_2_3/src/CMSSW_5_2_3_PatTuple.root'
#    '/store/user/ddidar/DATA_52X/DoubleElectron_Run2012B-PromptReco-v1_AOD/TwoLepSkim/193754_196509/PatTuple_DoubleElectron_Run2012B-PromptReco-v1_AOD_TwoLepSkim_193754_196509_1813_1_Il0.root'
#     '/store/user/ddidar/DATA_52X/ElectronHad_Run2012B-PromptReco-v1_AOD/TwoLepSkim/193754_196509/PatTuple_ElectronHad_Run2012B-PromptReco-v1_AOD_TwoLepSkim_193754_196509_999_1_PPc.root'
#    '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/MuHad_Run2012D-PromptReco-v1_AOD/OneLepOneBJetSkim/203768_207925/PatTuple_MuHad_Run2012D-PromptReco-v1_AOD_OneLepOneBJetSkim_203768_207925_214_1_YFH.root'
#    '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/ElectronHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193621/PatTuple_ElectronHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193621_110_1_24w.root'
#    '/cms/data/store/user/ddidar/DATA_53X/ElectronHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193754/PatTuple_ElectronHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193754_101_1_aAA.root',
#   '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/MuHad_Run2012B-13Jul2012-v1_AOD/OneLepOneBJetSkim/193833_196531/PatTuple_MuHad_Run2012B-13Jul2012-v1_AOD_OneLepOneBJetSkim_193833_196531_451_1_huk.root'
#    '/cms/data/store/user/ddidar/DATA_53X/MuHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193754/PatTuple_MuHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193754_332_1_G9q.root' 
#    '/scratch/lfs/lesya/FakeSync_jets_test_data_3rdShit.root'
#     '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/DoubleMu_Run2012D-PromptReco-v1_AOD/NOSkim/203767_207926/PatTuple_DoubleMu_Run2012D-PromptReco-v1_AOD_NOSkim_203767_207926_85_1_CCE.root'
    '/scratch/lfs/lesya/FakeSync_jets_test_data_Electrons.root'
    )
#    runOnData(process)
   

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )

## input files
#for i in inputFile.split(","):
#    print "Adding: ", i
#    process.source.fileNames.append(i)

#FakeElectrons
process.FakeElectrons = cms.EDAnalyzer("FakeMuonsSB",
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

#process.tauMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
#                                     src = cms.InputTag("hpsPFTauProducer"),
#                                     matched = cms.InputTag("prunedGenParticles"),
#                                     distMin = cms.double(0.15),
#                                     matchPDGId = cms.vint32()
#                                     )

#process.tauMCMatch = cms.EDProducer("MCMatcher",
#    				src     = cms.InputTag("hpsPFTauProducer"),
#    				matched = cms.InputTag("prunedGenParticles"),
#                   mcPdgId     = cms.vint32(),
#                   checkCharge = cms.bool(False),
#                   mcStatus = cms.vint32(),
#                   maxDeltaR = cms.double(0.15),
#                   maxDPtRel = cms.double(3.0),
#                   resolveAmbiguities = cms.bool(True),
#                   resolveByMatchQuality = cms.bool(False),
#)


#
# Example for a configuration of the MC match
# for taus (cuts are NOT tuned)
# (using old values from TQAF, january 2008)
#
process.tauMatch = cms.EDProducer("MCMatcher",
                          src         = cms.InputTag("hpsPFTauProducer"),          # RECO objects to match
                          matched     = cms.InputTag("prunedGenParticles"),              # mc-truth particle collection
                          mcPdgId     = cms.vint32(15),                            # one or more PDG ID (15 = tau); absolute values (see below)
                          checkCharge = cms.bool(True),                            # True = require RECO and MC objects to have the same charge
                          mcStatus    = cms.vint32(2),                             # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                          # NOTE that Taus can only be status 3 or 2, never 1!
                          maxDeltaR   = cms.double(999.9),                         # Minimum deltaR for the match.     By default any deltaR is allowed (why??)
                          maxDPtRel   = cms.double(999.9),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                          resolveAmbiguities    = cms.bool(True),                  # Forbid two RECO objects to match to the same GEN object
                          resolveByMatchQuality = cms.bool(False),                 # False = just match input in order; True = pick lowest deltaR pair first
                          )

process.tauGenJetMatch = cms.EDProducer("GenJetMatcher",             # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                src         = cms.InputTag("hpsPFTauProducer"),          # RECO jets (any View<Jet> is ok)
                                matched     = cms.InputTag("tauGenJets"),  # GEN jets  (must be GenJetCollection) "tauGenJetsSelectorAllHadrons"
                                mcPdgId     = cms.vint32(),                              # n/a
                                mcStatus    = cms.vint32(),                              # n/a
                                checkCharge = cms.bool(False),                           # n/a
                                maxDeltaR   = cms.double(0.1),                           # Minimum deltaR for the match
                                maxDPtRel   = cms.double(3.0),                           # Minimum deltaPt/Pt for the match
                                resolveAmbiguities    = cms.bool(True),                  # Forbid two RECO objects to match to the same GEN object
                                resolveByMatchQuality = cms.bool(False),                 # False = just match input in order; True = pick lowest deltaR pair first
                                )
process.tauGenJets = cms.EDProducer(
                            "TauGenJetProducer",
                            GenParticles =  cms.InputTag('prunedGenParticles'),
                            includeNeutrinos = cms.bool( False ),
                            verbose = cms.untracked.bool( False )
                            )

process.patMCTruth_Tau =  cms.Sequence ( process.tauMatch+
                                process.tauGenJets*
                                process.tauGenJetMatch )


process.electronMatch = cms.EDProducer("MCMatcher",
                                  src         = cms.InputTag("slimmedElectrons"),          # RECO objects to match
                                  matched     = cms.InputTag("prunedGenParticles"),              # mc-truth particle collection
                                  mcPdgId     = cms.vint32(),                            # one or more PDG ID (15 = tau); absolute values (see below)
                                  checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
                                  mcStatus    = cms.vint32(),                             # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                  # NOTE that Taus can only be status 3 or 2, never 1!
                                  maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed (why??)
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
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.electronPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                             src = cms.InputTag("slimmedElectrons"),      # RECO objects to match
                                             matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                             mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                             mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                             checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                             maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                             maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                             resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                             resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                             )

process.muonPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                         src = cms.InputTag("slimmedMuons"),      # RECO objects to match
                                         matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                         mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                         mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                         checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                         maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                         maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                         resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                         resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                         )

process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" # checks for fake PVs automatically
                                                  , filterParams =cms.PSet(
                                                                           minNdof = cms.double( 4. )
                                                                           , maxZ    = cms.double( 24. )
                                                                           , maxRho  = cms.double( 2. ) )		
                                                  , filter       = cms.bool( False ) # use only as producer
                                                  , src          = cms.InputTag( 'offlineSlimmedPrimaryVertices' )
)

if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
#   *process.tauMatch
   *process.muonMatch
   *process.electronMatch
   *process.jetPartonMatch
   *process.electronPartonMatch
   *process.muonPartonMatch
    *process.FakeElectrons
        )
else:
   import PhysicsTools.PythonAnalysis.LumiList as LumiList
#   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/hpc/lesya/CMSSW_5_3_6_patch1/src/SUSYAnalyzer/PatAnalyzer/test/JSON/Merged_190456-208686_8TeV_Collision12.json').getVLuminosityBlockRange()
   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test/fullJSON_SUSYFakes.json').getVLuminosityBlockRange() 
   process.p = cms.Path(
	process.goodOfflinePrimaryVertices
       *process.FakeElectrons
   )
    

