#ifndef __l1microgmtconfiguration_h
#define __l1microgmtconfiguration_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1TMuon/interface/L1TRegionalMuonCandidateFwd.h"
#include "DataFormats/L1TMuon/interface/L1TGMTInternalMuonFwd.h"
#include "DataFormats/L1TMuon/interface/L1TGMTInputCaloSumFwd.h"

#include <map>
#include <utility>

namespace l1t {
  class MicroGMTConfiguration {
    public:
      // All possible inputs for LUTs
      enum input_t { 
        PT, PT_COARSE, PHI, ETA, ETA_COARSE, QUALITY, DELTA_ETA_RED, DELTA_PHI_RED
      };

      // All possible input muon types
      enum muon_t {
        BARRELTF = 0, OVERLAPTF_NEG = 1, OVERLAPTF_POS = 2, FORWARDTF_NEG = 3, FORWARDTF_POS = 4, UNSET
      };
      typedef std::pair<input_t, int> PortType; 
      typedef L1TRegionalMuonCandidateCollection InputCollection;
      typedef MuonBxCollection OutputCollection;
      typedef Muon OutMuon;
      typedef L1TGMTInternalMuon InterMuon;
      typedef L1TGMTInternalMuonCollection InterMuonCollection;
      typedef L1TGMTInternalMuonList InterMuonList;
      typedef L1TGMTInputCaloSum CaloInput;
      typedef L1TGMTInputCaloSumCollection CaloInputCollection;
      // Two's complement for a given bit-length
      static unsigned getTwosComp(const int signed_int, const int width);

      
  };
}
#endif /* defined (__l1microgmtconfiguration_h) */ 