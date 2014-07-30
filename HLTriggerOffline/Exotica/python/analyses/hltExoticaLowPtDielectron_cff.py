import FWCore.ParameterSet.Config as cms

LowPtDielectronPSet = cms.PSet(
    hltPathsToCheck = cms.vstring(
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v", # Run2
        "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v" # Run1
        ),
    recElecLabel  = cms.InputTag("gedGsfElectrons"),
    # -- Analysis specific cuts
    minCandidates = cms.uint32(2),
    # -- Analysis specific binnings

    parametersTurnOn = cms.vdouble( 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                                    60, 70, 80, 100
                                   ),
    )
