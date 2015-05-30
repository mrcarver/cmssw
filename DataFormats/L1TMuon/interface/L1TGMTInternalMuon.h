#ifndef __l1t_gmt_internal_muon_h__
#define __l1t_gmt_internal_muon_h__

#include "DataFormats/L1TMuon/interface/L1TRegionalMuonCandidate.h"

namespace l1t {
class L1TGMTInternalMuon : public L1TRegionalMuonCandidate {
  public:
    L1TGMTInternalMuon() : 
      L1TRegionalMuonCandidate(), type_(MicroGMTConfiguration::muon_t::UNSET), hwRank_(0), hwCancelBit_(0), hwWins_(0), 
        hwIsoSum_(0), hwDeltaEta_(0), hwDeltaPhi_(0), hwAbsIso_(0), hwRelIso_(0), hwCaloIndex_(-1,-1)
      {};

    L1TGMTInternalMuon(int pt, int phi, int eta, int sign, int signvalid, int quality, int rank, int cancelbit) : 
      L1TRegionalMuonCandidate(pt, phi, eta, sign, signvalid, quality), hwRank_(rank), hwCancelBit_(cancelbit), hwIsoSum_(0), 
        hwDeltaEta_(0), hwDeltaPhi_(0), hwAbsIso_(0), hwRelIso_(0), hwCaloIndex_(-1,-1)
      {};

    L1TGMTInternalMuon(const L1TRegionalMuonCandidate& other) : L1TRegionalMuonCandidate(other), type_(MicroGMTConfiguration::muon_t::UNSET), hwRank_(0), hwCancelBit_(0), hwIsoSum_(0), hwDeltaEta_(0), hwDeltaPhi_(0), hwAbsIso_(-1), hwRelIso_(-1)
    {};

    virtual ~L1TGMTInternalMuon() {};

    void setHwCancelBit(int bit) { hwCancelBit_ = bit; };
    void setHwRank(int bits) { hwRank_ = bits; };
    void setHwWins(int wins) { hwWins_ = wins; };
    void increaseWins() { hwWins_++; };
    void setHwIsoSum(int isosum) { hwIsoSum_ = isosum; };
    void setHwAbsIso(int iso) { hwAbsIso_ = iso; };
    void setHwRelIso(int iso) { hwRelIso_ = iso; };
    void setType(MicroGMTConfiguration::muon_t type) { type_ = type; };
    void setExtrapolation(int deta, int dphi) { hwDeltaEta_ = deta; hwDeltaPhi_ = dphi; };
    void setHwCaloEta(int idx) { hwCaloIndex_.second = idx; }
    void setHwCaloPhi(int idx) { hwCaloIndex_.first = idx; }

    const int hwCancelBit() const { return hwCancelBit_; };
    const int hwRank() const { return hwRank_; };
    const int hwWins() const { return hwWins_; };
    const int hwIsoSum() const { return hwIsoSum_; };
    const int hwDEta() const { return hwDeltaEta_; };
    const int hwDPhi() const { return hwDeltaPhi_; };
    const int hwAbsIso() const { return hwAbsIso_; };
    const int hwRelIso() const { return hwRelIso_; };
    const MicroGMTConfiguration::muon_t type() const { return type_; };
    const int hwCaloEta() const { return hwCaloIndex_.second; }
    const int hwCaloPhi() const { return hwCaloIndex_.first; }
    
  private:
    MicroGMTConfiguration::muon_t type_;
    int hwRank_;
    int hwCancelBit_;
    int hwWins_;
    int hwIsoSum_;
    int hwDeltaEta_;
    int hwDeltaPhi_;
    int hwAbsIso_;
    int hwRelIso_;
    std::pair<int, int> hwCaloIndex_;
};

}

#endif /* define __l1t_gmt_internal_muon_h__ */