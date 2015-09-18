#ifndef uGMTcollections_h
#define uGMTcollections_h

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "EventFilter/L1TRawToDigi/interface/UnpackerCollections.h"

namespace l1t {
   class UGMTcollections : public UnpackerCollections {
      public:
         UGMTcollections(edm::Event& e) :
            UnpackerCollections(e),
            regionalMuonCandsBMTF_(new RegionalMuonCandBxCollection()),
            regionalMuonCandsOMTF_(new RegionalMuonCandBxCollection()),
            regionalMuonCandsEMTF_(new RegionalMuonCandBxCollection()),
            muons_(new MuonBxCollection()) {};

         virtual ~UGMTcollections();

         inline RegionalMuonCandBxCollection* getRegionalMuonCandsBMTF() { return regionalMuonCandsBMTF_.get(); };
         inline RegionalMuonCandBxCollection* getRegionalMuonCandsOMTF() { return regionalMuonCandsOMTF_.get(); };
         inline RegionalMuonCandBxCollection* getRegionalMuonCandsEMTF() { return regionalMuonCandsEMTF_.get(); };
         inline MuonBxCollection* getMuons() { return muons_.get(); };

      private:
         std::auto_ptr<RegionalMuonCandBxCollection> regionalMuonCandsBMTF_;
         std::auto_ptr<RegionalMuonCandBxCollection> regionalMuonCandsOMTF_;
         std::auto_ptr<RegionalMuonCandBxCollection> regionalMuonCandsEMTF_;
         std::auto_ptr<MuonBxCollection> muons_;
   };
}

#endif
