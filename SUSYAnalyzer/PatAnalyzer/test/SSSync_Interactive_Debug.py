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
outputFile=""
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
process.MessageLogger.cerr.FwkReport.reportEvery = 100

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

# START ELECTRON ID SECTION
#
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#

# Load tools and function definitions
#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

#process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

#from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
#process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)


#add in the heep ID to the producer. You can run with other IDs but heep ID must be loaded with setupVIDSelection, not setupAllVIDSelection as heep works differently because mini-aod and aod are defined in the same file to ensure consistancy (you cant change cuts of aod without changing miniaod
#process.load('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff')
#setupVIDSelection(process.egmGsfElectronIDs,process.heepElectronID_HEEPV51)

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

# Do not forget to add the egmGsfElectronIDSequence to the path,
# as in the example below!

#
# END ELECTRON ID SECTION
#


process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'), 
							 #eventsToProcess = cms.untracked.VEventRange('1:40603'),                   
                             fileNames = cms.untracked.vstring(),      
                             )

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

if isMC:
    process.source.fileNames = cms.untracked.vstring(   


    #'/cms/data/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/022B08C4-C702-E511-9995-D4856459AC30.root',
	 #'/store/relval/CMSSW_7_4_6_patch6/RelValSMS-T1tttt_13/GEN-SIM-DIGI-RAW-HLTDEBUG/74X_mcRun2_asymptotic_realisticBS_v0_2015Jul24PU-v1/00000/10243BAB-4932-E511-A27D-0025905A60D2.root',
	 
	 '/cms/data/store/mc/RunIISpring15DR74/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/30000/60087A61-9134-E511-B0C6-0025905B855E.root',
	 
	 #'/cms/data/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/A0523330-7214-E511-BB04-00259073E506.root',
	 #'/cms/data/store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/A0707FB3-32FD-E411-810D-00259074AEDE.root',
	 # '/cms/data/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/4874FDF3-9F02-E511-B0B0-008CFA0A59E0.root',
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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

#process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )
process.TFileService = cms.Service("TFileService", fileName = cms.string("TTW_sync_oldpTrel.root") )

## input files
#for i in inputFile.split(","):
#   print "Adding: ", i
#   process.source.fileNames.append(i)

process.SyncExercise = cms.EDAnalyzer("SSSync",#"FakeMuonsYa",
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       JetLabel = cms.InputTag("slimmedJets"),
                                       BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                       HLTResultsLabel = cms.InputTag("TriggerResults"),
									   MVAId = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"), #not used when running on PAT
                                       METLabel = cms.InputTag("slimmedMETs"),
                                       #METFilter = cms.InputTag("TriggerResults::PAT"),
                                       qualityCuts = PFTauQualityCuts,
                                       SampleLabel = cms.untracked.string("ElectronsMC"), # defines a piece of code to run; helps to avoid code recompilation
									   doFR = cms.untracked.bool(False),
									   id = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51")
                                       )



process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" # checks for fake PVs automatically
                                                  , filterParams =cms.PSet(
                                                                           minNdof = cms.double( 4. )
                                                                           , maxZ    = cms.double( 24. )
                                                                           , maxRho  = cms.double( 2. ) )		
                                                  , filter       = cms.bool( False ) # use only as producer
                                                  , src          = cms.InputTag( 'offlineSlimmedPrimaryVertices' )
)


process.writeedmfile = cms.OutputModule("PoolOutputModule",
                               outputCommands=cms.untracked.vstring(
                                   'keep *',
                                   ),
                               fileName=cms.untracked.string("MVAoutput.root"),
                               )
							   
#process.outp = cms.EndPath(process.writeedmfile)

if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
	#*process.egmGsfElectronIDSequence
	*process.electronMVAValueMapProducer
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
    

sched = cms.Schedule(process.p)#, process.outp)
