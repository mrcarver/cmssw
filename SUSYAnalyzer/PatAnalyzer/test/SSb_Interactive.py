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
outputFile="results/SSb_Interactive.root"
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

#process = cms.Process("SSb")
process = cms.Process("pippo")
process.load('Configuration.StandardSequences.Reconstruction_cff')

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

if isMC:
    process.source.fileNames = cms.untracked.vstring(               
     '/cms/data/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/003E832C-8AFC-E311-B7AA-002590596490.root'
)



else:
    process.source.fileNames = cms.untracked.vstring( 
#     '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/ElectronHad_Run2012A-13Jul2012-v1_AOD/TwoLepSkim/190456_193621/PatTuple_ElectronHad_Run2012A-13Jul2012-v1_AOD_TwoLepSkim_190456_193621_368_2_D2g.root'
#      '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/MuEG_Run2012D-PromptReco-v1_AOD/TwoLepSkim/203768_207925/PatTuple_MuEG_Run2012D-PromptReco-v1_AOD_TwoLepSkim_203768_207925_1047_1_ovT.root'
     '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/DoubleMu_Run2012D-PromptReco-v1_AOD/TwoLepSkim/207925_209152/PatTuple_DoubleMu_Run2012D-PromptReco-v1_AOD_TwoLepSkim_207925_209152_195_2_xZU.root'
    )
#    runOnData(process)
   

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )

## input files
#for i in inputFile.split(","):
#    print "Adding: ", i
#    process.source.fileNames.append(i)

#SSb
process.SSb = cms.EDAnalyzer("SSb",
                                    MuonLabel = cms.InputTag("patMuons"),
                                    ElectronLabel = cms.InputTag("patElectrons"),
                                    TauLabel = cms.InputTag("hpsPFTauProducer"),
                                    TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
                                    JetLabel = cms.InputTag("patJetsAK5PF"),
                                    BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                    HLTResultsLabel = cms.InputTag("TriggerResults::HLT"),
                                    METLabel = cms.InputTag("patMETsPF"),
                                    METFilter = cms.InputTag("TriggerResults::PatAnalyzer"),
                                    qualityCuts = PFTauQualityCuts,
                                    SampleLabel = cms.untracked.string("MCsample") # defines a piece of code to run; helps to avoid code recompilation
                                    )

#process.tauMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
#                                     src = cms.InputTag("hpsPFTauProducer"),
#                                     matched = cms.InputTag("genParticles"),
#                                     distMin = cms.double(0.15),
#                                     matchPDGId = cms.vint32()
#                                     )

#process.tauMCMatch = cms.EDProducer("MCMatcher",
#    				src     = cms.InputTag("hpsPFTauProducer"),
#    				matched = cms.InputTag("genParticles"),
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
                          matched     = cms.InputTag("genParticles"),              # mc-truth particle collection
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
                            GenParticles =  cms.InputTag('genParticles'),
                            includeNeutrinos = cms.bool( False ),
                            verbose = cms.untracked.bool( False )
                            )

process.patMCTruth_Tau =  cms.Sequence ( process.tauMatch+
                                process.tauGenJets*
                                process.tauGenJetMatch )

process.muonMatch = cms.EDProducer("MCMatcher",
                                   src         = cms.InputTag("patMuons"),              # RECO objects to match
                                   matched     = cms.InputTag("genParticles"),          # mc-truth particle collection
                                   mcPdgId     = cms.vint32(),                          # one or more PDG ID (15 = tau); absolute values (see below)
                                   checkCharge = cms.bool(False),                       # True = require RECO and MC objects to have the same charge
                                   mcStatus    = cms.vint32(),                          # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                   # NOTE that Taus can only be status 3 or 2, never 1!
                                   maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed
                                   maxDPtRel   = cms.double(10.),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                                   resolveAmbiguities    = cms.bool(False),                  # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
                                   )
process.electronMatch = cms.EDProducer("MCMatcher",
                                  src         = cms.InputTag("patElectrons"),          # RECO objects to match
                                  matched     = cms.InputTag("genParticles"),              # mc-truth particle collection
                                  mcPdgId     = cms.vint32(),                            # one or more PDG ID (15 = tau); absolute values (see below)
                                  checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
                                  mcStatus    = cms.vint32(),                             # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                  # NOTE that Taus can only be status 3 or 2, never 1!
                                  maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed (why??)
                                  maxDPtRel   = cms.double(10.),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                                  resolveAmbiguities    = cms.bool(False),                  # Forbid two RECO objects to match to the same GEN object
                                  resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
                                  )

## reco-generator(parton) matching for jets
process.jetPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("patJetsAK5PF"),      # RECO objects to match
                                        matched = cms.InputTag("genParticles"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 6, 21),# one or more PDG ID (quarks except top; gluons)
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
                                                  , src          = cms.InputTag( 'offlinePrimaryVertices' )
)

if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
    #   *process.tauMatch
    *process.muonMatch
    *process.electronMatch
    *process.jetPartonMatch
#    *process.patMCTruth_Tau
    *process.SSb
        )
else: 
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
       *process.SSb
    )
    

