#ifndef __l1microgmtcanceloutunit_h
#define __l1microgmtcanceloutunit_h

#include "MicroGMTConfiguration.h"
#include "MicroGMTMatchQualLUT.h"

namespace l1t {
  class MicroGMTCancelOutUnit {
    public: 
      explicit MicroGMTCancelOutUnit (const edm::ParameterSet&);
      virtual ~MicroGMTCancelOutUnit ();
      // Cancel-Out is set to 1 for the lower quality muon, if a match is found according to match LUTs
      void setCancelOutBits(MicroGMTConfiguration::InterMuonList&);
    private:
      // This goes through two neighboring sections and checks for matches
      void getCancelOutBits( std::vector<MicroGMTConfiguration::InterMuonList::iterator> &, std::vector<MicroGMTConfiguration::InterMuonList::iterator> &);

      MicroGMTMatchQualLUT m_boPosMatchQualLUT;
      MicroGMTMatchQualLUT m_boNegMatchQualLUT;
      MicroGMTMatchQualLUT m_foPosMatchQualLUT;
      MicroGMTMatchQualLUT m_foNegMatchQualLUT;
      MicroGMTMatchQualLUT m_brlSingleMatchQualLUT;
      MicroGMTMatchQualLUT m_ovlPosSingleMatchQualLUT;
      MicroGMTMatchQualLUT m_ovlNegSingleMatchQualLUT;
      MicroGMTMatchQualLUT m_fwdPosSingleMatchQualLUT;
      MicroGMTMatchQualLUT m_fwdNegSingleMatchQualLUT;
      std::map<int, MicroGMTMatchQualLUT*> m_lutDict;
  };
}
#endif /* defined(__l1microgmtcanceloutunit_h) */