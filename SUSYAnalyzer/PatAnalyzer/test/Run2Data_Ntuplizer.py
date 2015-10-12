import sys
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts


skim     ="NOSkim"
isMC     = bool(False)
isMSUGRA = bool(False)
isSMS    = bool(False)
useCrab  = bool(False)
doSusyTopProjection = bool(False)
#inputFile=""
#outputFile=""
outputFile="Purities_fuck.root"
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
    #process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco
	process.GlobalTag.globaltag="GR_P_V56::All" # tag fo 74X 2015B July2015 

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')


process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             #duplicateCheckMode = cms.untracked.string('noDuplicateCheck'), 
							 #eventsToProcess = cms.untracked.VEventRange('251643:59426492'),        
							 #lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_7_4_7/src/SUSYAnalyzer/PatAnalyzer/test/JSON/json_DCSONLY_Run2015B.txt').getVLuminosityBlockRange()           
                             fileNames = cms.untracked.vstring(),      
                             )

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

if isMC:
    process.source.fileNames = cms.untracked.vstring(               

	 #'/cms/data/store/mc/Phys14DR/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0CD15D7F-4E6B-E411-AEB4-002590DB9216.root'
	 #'/cms/data/store/mc/Phys14DR/TT_Tune4C_13TeV-pythia8-tauola/MINIAODSIM/PU40bx25_tsg_PHYS14_25_V1-v1/00000/06E41ADB-7870-E411-8850-0025905A605E.root'
      #'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root',
	  #'/cms/data/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0EA5C14E-BC76-E411-8BBA-0025907DC9D0.root',
	  #'/cms/data/store/mc/RunIISpring15DR74/TT_TuneZ2star_13TeV-powheg-pythia6-tauola/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/002BE1CB-4307-E511-A582-008CFA197DB4.root'
	)

else:
   import FWCore.PythonUtilities.LumiList as LumiList
   #process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_7_4_7/src/SUSYAnalyzer/PatAnalyzer/test/JSON/json_DCSONLY_Run2015B.txt').getVLuminosityBlockRange()
   #process.source.lumisToProcess = LumiList.LumiList(filename = '/lfs/scratch/mrcarver/CMSSW_7_4_6/src/SUSYAnalyzer/PatAnalyzer/test/GoldenJSON.txt').getVLuminosityBlockRange()

   process.source.fileNames = cms.untracked.vstring( 
    #'/cms/data/store/data/Run2015B/MuonEG/MINIAOD/PromptReco-v1//000/251/643/00000/6EE03EBE-BE2C-E511-9EE7-02163E0138B3.root',
	'/cms/data/store/data/Run2015C/DoubleEG/MINIAOD/PromptReco-v1/000/254/199/00000/1AA8B4B1-C845-E511-BCE2-02163E012AC7.root',
	#'/cms/data/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/162/00000/BE906A2C-4327-E511-8014-02163E0121CC.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/162/00000/12284DB9-4227-E511-A438-02163E013674.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/164/00000/402F0995-A326-E511-86BB-02163E013948.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/167/00000/70C4A781-A826-E511-95B4-02163E013414.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/168/00000/627E9C65-DD26-E511-87FB-02163E013576.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/244/00000/E42FEF61-6E27-E511-B93A-02163E0143C0.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/251/00000/0292F6F9-8A27-E511-9074-02163E012213.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/252/00000/9ADEE140-9C27-E511-919A-02163E011D23.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/493/00000/323EBCB2-D428-E511-86E5-02163E01463E.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/497/00000/826586F3-FC28-E511-89EE-02163E01340A.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/498/00000/20B6D9C3-FE28-E511-BD2D-02163E0117FF.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/499/00000/EA52602B-0929-E511-AB5B-02163E011955.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/500/00000/4897C658-4129-E511-AA73-02163E0135F3.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/521/00000/B4B3A942-6429-E511-8A0C-02163E01413E.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/522/00000/E817E288-5C29-E511-9EDC-02163E011A5A.root',
    #'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/560/00000/DA0DAC6D-DE29-E511-B863-02163E0133B6.root',
	#'/cms/data/store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/561/00000/5AC326A0-142A-E511-B70C-02163E0128E3.root',
	
	)
	#runOnData(process)
   

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("fuck.root") )

## input files
#for i in inputFile.split(","):
#   print "Adding: ", i
#   process.source.fileNames.append(i)

process.Run2Ntuplizer = cms.EDAnalyzer("SSSync",#"FakeMuonsYa",
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
                                       SampleLabel = cms.untracked.string("ElectronsData"), # defines a piece of code to run; helps to avoid code recompilation
									   doFR = cms.untracked.bool(True)
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
	*process.Run2Ntuplizer
        )
else:
   import FWCore.PythonUtilities.LumiList as LumiList
   #process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_7_4_7/src/SUSYAnalyzer/PatAnalyzer/test/JSON/json_DCSONLY_Run2015B.txt').getVLuminosityBlockRange() 
   #process.source.lumisToProcess = LumiList.LumiList(filename = '/lfs/scratch/mrcarver/CMSSW_7_4_6/src/SUSYAnalyzer/PatAnalyzer/test/GoldenJSON.txt').getVLuminosityBlockRange()
   process.p = cms.Path(
	process.goodOfflinePrimaryVertices
       *process.Run2Ntuplizer
   )
    

