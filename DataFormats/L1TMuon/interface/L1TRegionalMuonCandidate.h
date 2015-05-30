#ifndef __l1t_regional_muon_candidate_h__
#define __l1t_regional_muon_candidate_h__

namespace l1t {
class L1TRegionalMuonCandidate {
  public:
    L1TRegionalMuonCandidate() : 
      hwPt_(0), hwPhi_(0), hwEta_(0), hwSignBit(0), hwSignValid_(0), hwQuality_(0), hwTrackAddress_(0)
      {};

    L1TRegionalMuonCandidate(int pt, int phi, int eta, int sign, int signvalid, int quality) : 
      hwPt_(pt), hwPhi_(phi), hwEta_(eta), hwSignBit(sign), hwSignValid_(signvalid), hwQuality_(quality), hwTrackAddress_(0)
      {};

    virtual ~L1TRegionalMuonCandidate() {};

    void setHwPt(int bits) { hwPt_ = bits; };
    void setHwPhi(int bits) { hwPhi_ = bits; };
    void setHwEta(int bits) { hwEta_ = bits; };
    void setHwSign(int bits) { hwSignBit = bits; };
    void setHwSignValid(int bits) { hwSignValid_ = bits; };
    void setHwQual(int bits) { hwQuality_ = bits; };
    void setHwTrackAddress(int bits) { hwTrackAddress_ = bits; };
    void setLink(int link) { link_ = link; };

    const int hwPt() const { return hwPt_; };
    const int hwPhi() const { return hwPhi_; };
    const int hwEta() const { return hwEta_; };
    const int hwSign() const { return hwSignBit; };
    const int hwSignValid() const { return hwSignValid_; };
    const int hwQual() const { return hwQuality_; };
    const int hwTrackAddress() const { return hwTrackAddress_; };
    const int link() const { return link_; }
  private:
    int hwPt_;
    int hwPhi_;
    int hwEta_;
    int hwSignBit;
    int hwSignValid_;
    int hwQuality_;
    int hwTrackAddress_;
    int link_;
};

}

#endif /* define __l1t_regional_muon_candidate_h__ */