import FWCore.ParameterSet.Config as cms

HighPtDielectronPSet = cms.PSet(
    hltPathsToCheck = cms.vstring(
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v", # Run2
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v" # Run1
        ),
    recElecLabel  = cms.InputTag("gedGsfElectrons"),
    # -- Analysis specific cuts
    minCandidates = cms.uint32(2),
    # -- Analysis specific binnings
    parametersTurnOn = cms.vdouble( 0, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000,
                                    1100, 1200, 1500
                                   ),
    )
