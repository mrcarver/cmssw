#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "uGMTtokens.h"

namespace l1t {
   UGMTtokens::UGMTtokens(const edm::ParameterSet& cfg, edm::ConsumesCollector& cc) : PackerTokens(cfg, cc)
   {
      auto tag = cfg.getParameter<edm::InputTag>("InputLabel");

      regionalMuonCandTokenBMTF_ = cc.consumes<RegionalMuonCandBxCollection>(tag);
      regionalMuonCandTokenOMTF_ = cc.consumes<RegionalMuonCandBxCollection>(tag);
      regionalMuonCandTokenEMTF_ = cc.consumes<RegionalMuonCandBxCollection>(tag);
      muonToken_ = cc.consumes<MuonBxCollection>(tag);
   }
}
