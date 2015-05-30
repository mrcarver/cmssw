#ifndef __l1microgmtmatchquallut_h
#define __l1microgmtmatchquallut_h

#include "MicroGMTLUT.h"
#include "MicroGMTConfiguration.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace l1t {
    class MicroGMTMatchQualLUT : public MicroGMTLUT {
      public:
        MicroGMTMatchQualLUT () : m_dEtaRedMask(0), m_dPhiRedMask(0), m_dEtaRedInWidth(-1), m_dPhiRedInWidth(-1) {};
        explicit MicroGMTMatchQualLUT (const edm::ParameterSet&, std::string);
        virtual ~MicroGMTMatchQualLUT ();

        int lookup(int dEta, int dPhi) const;

        int hashInput(int dEta, int dPhi) const;
        void unHashInput(int input, int& dEta, int& dPhi) const;

        int getDeltaEtaWidth() const { return m_dEtaRedInWidth; }
        int getDeltaPhiWidth() const { return m_dPhiRedInWidth; }
      private:
        int m_dEtaRedMask; 
        int m_dPhiRedMask; 
        int m_dEtaRedInWidth;
        int m_dPhiRedInWidth;
    };
}
#endif /* defined(__l1microgmtmatchquallut_h) */