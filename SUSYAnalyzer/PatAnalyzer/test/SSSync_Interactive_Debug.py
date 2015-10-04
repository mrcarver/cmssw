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
#outputFile=""
#outputFile="SyncOutput/TriggerNamesTest.root"
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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if isMC:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    #process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco
	#process.GlobalTag.globaltag="PHYS14_25_V2::All" # tag fo PHYS14 25ns samples
	process.GlobalTag.globaltag="MCRUN2_74_V9::All" # tag for Run2 MC samples
else:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'), 
							 #eventsToProcess = cms.untracked.VEventRange('1:1557092937','1:1624289770','1:1137407433','1:1658883658'),                   
                             fileNames = cms.untracked.vstring(),      
                             )

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

if isMC:
    process.source.fileNames = cms.untracked.vstring(   
	
	 #     '/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0CD15D7F-4E6B-E411-AEB4-002590DB9216.root'
	 #     '/cms/data/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU40bx25_tsg_PHYS14_25_V1-v1/00000/06E41ADB-7870-E411-8850-0025905A605E.root'
     #      '/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
	  #     '/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0EA5C14E-BC76-E411-8BBA-0025907DC9D0.root',
	  
	#  '/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/641A7B60-BC76-E411-8221-0025904B578E.root',
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/64E57633-6F75-E411-AE68-00266CFFA604.root',
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66203DDB-3775-E411-9E9E-002590AC4C9E.root',
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66689298-3874-E411-BBFB-0025904B1340.root',
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/66A6025E-D275-E411-B485-002590DB9188.root,'
	
	
	
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/68548FD6-7C75-E411-8EC4-002590AC4B3E.root',
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/68B439F0-3174-E411-A50C-0025904B144E.root',
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6A99352C-9F74-E411-8244-002590A2CD44.root',
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6C725C73-0A75-E411-91E3-003048D4364C.root',
	#'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/6E47D20A-0276-E411-9773-00215AD4D6E2.root',

	 #'/cms/data/store/relval/CMSSW_7_4_3/RelValTTbar_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/14354550-6F0C-E511-9702-0025905964BA.root',
	 #'/cms/data/store/relval/CMSSW_7_4_3/RelValTTbar_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/A6642653-6F0C-E511-8A33-003048FFD796.root',
	 #'/cms/data/store/relval/CMSSW_7_4_3/RelValTTbar_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/FA7DF3CE-6D0C-E511-AAFA-002618943870.root',
	 
	 #'/cms/data/store/relval/CMSSW_7_4_3/RelValTTbar_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/00FF2F84-860C-E511-9EE9-0025905A6136.root',
	 #'/cms/data/store/relval/CMSSW_7_4_3/RelValTTbar_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/36D3684E-D30C-E511-8384-003048FFD754.root',
	 #'/cms/data/store/relval/CMSSW_7_4_3/RelValTTbar_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/AC892B8D-730C-E511-A640-0025905964B2.root',
	 #'/cms/data/store/relval/CMSSW_7_4_3/RelValTTbar_13/MINIAODSIM/PU25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/D0C9CA51-D30C-E511-AD93-0025905964B4.root'
	 
	 #'/store/relval/CMSSW_7_4_3/RelValSMS-T1tttt_mGl-1500_mLSP-100_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/4C8CCB35-740C-E511-9CCE-0025905964B4.root',
     # '/store/relval/CMSSW_7_4_3/RelValSMS-T1tttt_mGl-1500_mLSP-100_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/BA6F98D8-710C-E511-84EC-003048FF86CA.root'
	  
	  	#'/store/relval/CMSSW_7_4_3/RelValZMM_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/00932919-730C-E511-988E-002618943981.root',
       #'/store/relval/CMSSW_7_4_3/RelValZMM_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/0C97DD1D-730C-E511-9BB9-002618FDA287.root'
	   
	   #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/0A6800BE-720C-E511-846C-0025905B85EE.root',
       #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/8CE318BC-720C-E511-895B-0025905AA9CC.root',
       #'/store/relval/CMSSW_7_4_3/RelValZEE_13/MINIAODSIM/PUpmx25ns_MCRUN2_74_V9_fixMem_FastSim-v1/00000/EAA3A6BF-720C-E511-ADDA-0026189438CE.root'
	
    #'/cms/data/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/022B08C4-C702-E511-9995-D4856459AC30.root',
	 '/store/relval/CMSSW_7_4_6_patch6/RelValSMS-T1tttt_13/GEN-SIM-DIGI-RAW-HLTDEBUG/74X_mcRun2_asymptotic_realisticBS_v0_2015Jul24PU-v1/00000/10243BAB-4932-E511-A27D-0025905A60D2.root',
	 
	 '/cms/data/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/A0523330-7214-E511-BB04-00259073E506.root',
	 '/cms/data/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/A0707FB3-32FD-E411-810D-00259074AEDE.root',
	  '/cms/data/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/4874FDF3-9F02-E511-B0B0-008CFA0A59E0.root',
	)

else:
   import PhysicsTools.PythonAnalysis.LumiList as LumiList
#   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test/fullJSON_SUSYFakes.json').getVLuminosityBlockRange()

   process.source.fileNames = cms.untracked.vstring( 
#    #'file:../../../CMSSW_5_2_3_PatTuple_test_PFiso_El_newMETCor.root'
#    'file:/raid/raid10/ddidar/CMSSW_5_2_3/src/CMSSW_5_2_3_PatTuple.root'
#    '/store/user/ddidar/DATA_52X/DoubleElectron_Run2012B-PromptReco-v1_AOD/TwoLepSkim/193754_196509/PatTuple_DoubleElectron_Run2012B-PromptReco-v1_AOD_TwoLepSkim_193754_196509_1813_1_Il0.root'
#     '/store/user/ddidar/DATA_52X/ElectronHad_Run2012B-PromptReco-v1_AOD/TwoLepSkim/193754_196509/PatTuple_ElectronHad_Run2012B-PromptReco-v1_AOD_TwoLepSkim_193754_196509_999_1_PPc.root'
#         '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/MuHad_Run2012D-PromptReco-v1_AOD/OneLepOneBJetSkim/203768_207925/PatTuple_MuHad_Run2012D-PromptReco-v1_AOD_OneLepOneBJetSkim_203768_207925_214_1_YFH.root'
#         '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/ElectronHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193621/PatTuple_ElectronHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193621_110_1_24w.root'
#         '/cms/data/store/user/ddidar/DATA_53X/ElectronHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193754/PatTuple_ElectronHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193754_101_1_aAA.root',
#        '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/MuHad_Run2012B-13Jul2012-v1_AOD/OneLepOneBJetSkim/193833_196531/PatTuple_MuHad_Run2012B-13Jul2012-v1_AOD_OneLepOneBJetSkim_193833_196531_451_1_huk.root'
#         '/cms/data/store/user/ddidar/DATA_53X/MuHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193754/PatTuple_MuHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193754_332_1_G9q.root' 
#    '/scratch/lfs/lesya/FakeSync_jets_test_data_3rdShit.root'
#          '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/DoubleMu_Run2012D-PromptReco-v1_AOD/NOSkim/203767_207926/PatTuple_DoubleMu_Run2012D-PromptReco-v1_AOD_NOSkim_203767_207926_85_1_CCE.root'
    #'/scratch/lfs/lesya/FakeSync_jets_test_data_Electrons.root'
    )
#    runOnData(process)
   

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
									#SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

#process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )
process.TFileService = cms.Service("TFileService", fileName = cms.string("DMPathTest.root") )

## input files
#for i in inputFile.split(","):
#   print "Adding: ", i
#   process.source.fileNames.append(i)

process.SyncExercise = cms.EDAnalyzer("SSSync",#"FakeMuonsYa",
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       #TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
                                       JetLabel = cms.InputTag("slimmedJets"),
                                       BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                       HLTResultsLabel = cms.InputTag("TriggerResults"),
									   MVAId = cms.InputTag("mvaTrigV050nsCSA14","","addMVAid"), #not used when running on PAT
                                       METLabel = cms.InputTag("slimmedMETs"),
                                       #METFilter = cms.InputTag("TriggerResults::PAT"),
                                       qualityCuts = PFTauQualityCuts,
                                       SampleLabel = cms.untracked.string("ElectronsMC"), # defines a piece of code to run; helps to avoid code recompilation
									   doFR = cms.untracked.bool(True)
                                       )

#FakeElectrons
#process.FakeElectrons = cms.EDAnalyzer("SSSync",#"FakeMuonsYa",
#                                       MuonLabel = cms.InputTag("slimmedMuons"),
#                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
#                                       TauLabel = cms.InputTag("slimmedTaus"),
#                                       #TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
#                                       JetLabel = cms.InputTag("slimmedJets"),
#                                      BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
#                                      HLTResultsLabel = cms.InputTag("TriggerResults"),
#                                      METLabel = cms.InputTag("slimmedMETs"),
#                                      #METFilter = cms.InputTag("TriggerResults::PAT"),
#                                      qualityCuts = PFTauQualityCuts,
#                                       SampleLabel = cms.untracked.string("ElectronsMC") # defines a piece of code to run; helps to avoid code recompilation
#                                       )

#process.tauMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
#                                     src = cms.InputTag("slimmedTaus"),
#                                     matched = cms.InputTag("prunedGenParticles"),
#                                     distMin = cms.double(0.15),
#                                     matchPDGId = cms.vint32()
#                                     )

#process.tauMCMatch = cms.EDProducer("MCMatcher",
#    				src     = cms.InputTag("slimmedTaus"),
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
                          src         = cms.InputTag("slimmedTaus"),          # RECO objects to match
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
                                src         = cms.InputTag("slimmedTaus"),          # RECO jets (any View<Jet> is ok)
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
                            prunedGenParticles =  cms.InputTag('prunedGenParticles'),
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
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, -1, -2, -3, -4, -5, 21),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
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

process.jetPartonMatchJets = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("ak5GenJets"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(3.0),                # Minimum deltaPt/Pt for the match
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
#   *process.muonMatch
#   *process.electronMatch
#   *process.jetPartonMatch
#   *process.jetPartonMatch2
#   *process.electronPartonMatch
#   *process.muonPartonMatch
 #   *process.FakeElectrons
 	#*process.mvaNonTrigV025nsPHYS14
	*process.SyncExercise
        )
else:
   import PhysicsTools.PythonAnalysis.LumiList as LumiList
#   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/hpc/lesya/CMSSW_5_3_6_patch1/src/SUSYAnalyzer/PatAnalyzer/test/JSON/Merged_190456-208686_8TeV_Collision12.json').getVLuminosityBlockRange()
   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test/fullJSON_SUSYFakes.json').getVLuminosityBlockRange() 
   process.p = cms.Path(
	process.goodOfflinePrimaryVertices
       *process.FakeElectrons
   )
    

