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
process = cms.Process("AnalysisTuple")

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
							 #eventsToProcess = cms.untracked.VEventRange('1:85776299','1:99480374'),                   
                             fileNames = cms.untracked.vstring(),      
                             )

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

if isMC:
    process.source.fileNames = cms.untracked.vstring(            
	
	)

else:
   import PhysicsTools.PythonAnalysis.LumiList as LumiList
#   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test/fullJSON_SUSYFakes.json').getVLuminosityBlockRange()

   process.source.fileNames = cms.untracked.vstring( 
   
    )
#    runOnData(process)
   

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("newTrig.root") )

## input files
for i in inputFile.split(","):
   print "Adding: ", i
   process.source.fileNames.append(i)

process.SyncExercise = cms.EDAnalyzer("SSSync",#"FakeMuonsYa",
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       #TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
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

if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
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
    

