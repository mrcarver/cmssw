#include "../interface/MicroGMTIsolationUnit.h"

#include "DataFormats/L1TMuon/interface/L1TGMTInputCaloSum.h"
#include "DataFormats/L1TMuon/interface/L1TGMTInternalMuon.h"
#include "DataFormats/L1TMuon/interface/L1TRegionalMuonCandidate.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


l1t::MicroGMTIsolationUnit::MicroGMTIsolationUnit (const edm::ParameterSet& iConfig) :
  m_BEtaExtrapolation(iConfig, "BEtaExtrapolationLUTSettings", 0), m_BPhiExtrapolation(iConfig, "BPhiExtrapolationLUTSettings", 1), m_OEtaExtrapolation(iConfig, "OEtaExtrapolationLUTSettings", 0),
  m_OPhiExtrapolation(iConfig, "OPhiExtrapolationLUTSettings", 1), m_FEtaExtrapolation(iConfig, "FEtaExtrapolationLUTSettings", 0), m_FPhiExtrapolation(iConfig, "FPhiExtrapolationLUTSettings", 1),
  m_IdxSelMemEta(iConfig, "IdxSelMemEtaLUTSettings", 0), m_IdxSelMemPhi(iConfig, "IdxSelMemPhiLUTSettings", 1), m_RelIsoCheckMem(iConfig, "RelIsoCheckMemLUTSettings"),  
  m_AbsIsoCheckMem(iConfig, "AbsIsoCheckMemLUTSettings"), m_initialSums(false)
{
  m_etaExtrapolationLUTs[MicroGMTConfiguration::muon_t::BARRELTF] = &m_BEtaExtrapolation;
  m_phiExtrapolationLUTs[MicroGMTConfiguration::muon_t::BARRELTF] = &m_BPhiExtrapolation;
  m_etaExtrapolationLUTs[MicroGMTConfiguration::muon_t::OVERLAPTF_POS] = &m_OEtaExtrapolation;
  m_etaExtrapolationLUTs[MicroGMTConfiguration::muon_t::OVERLAPTF_NEG] = &m_OEtaExtrapolation;
  m_phiExtrapolationLUTs[MicroGMTConfiguration::muon_t::OVERLAPTF_POS] = &m_OPhiExtrapolation;
  m_phiExtrapolationLUTs[MicroGMTConfiguration::muon_t::OVERLAPTF_NEG] = &m_OPhiExtrapolation;
  m_etaExtrapolationLUTs[MicroGMTConfiguration::muon_t::FORWARDTF_POS] = &m_FEtaExtrapolation;
  m_etaExtrapolationLUTs[MicroGMTConfiguration::muon_t::FORWARDTF_NEG] = &m_FEtaExtrapolation;
  m_phiExtrapolationLUTs[MicroGMTConfiguration::muon_t::FORWARDTF_POS] = &m_FPhiExtrapolation;
  m_phiExtrapolationLUTs[MicroGMTConfiguration::muon_t::FORWARDTF_NEG] = &m_FPhiExtrapolation;
}

l1t::MicroGMTIsolationUnit::~MicroGMTIsolationUnit ()
{
}

int 
l1t::MicroGMTIsolationUnit::getCaloIndex(MicroGMTConfiguration::InterMuon& mu) const 
{
  // handle the wrap-around of phi:
  int phi = (mu.hwPhi() + mu.hwDPhi())%576;
  if (phi < 0) {
    phi = 576+phi;
  }

  int phiIndex = m_IdxSelMemPhi.lookup(phi);
  int eta = mu.hwEta()+mu.hwDEta();
  eta = MicroGMTConfiguration::getTwosComp(eta, 9);
  int etaIndex = m_IdxSelMemEta.lookup(eta);
  mu.setHwCaloEta(etaIndex);
  mu.setHwCaloPhi(phiIndex);

  return phiIndex + etaIndex*36; 
}

void 
l1t::MicroGMTIsolationUnit::extrapolateMuons(MicroGMTConfiguration::InterMuonList& inputmuons) const {
  MicroGMTConfiguration::InterMuonList::iterator mu;
  for (mu = inputmuons.begin(); mu != inputmuons.end(); ++mu) {   
    // only use 6 LSBs of pt:
    int ptRed = mu->hwPt() & 0b111111;
    // here we drop the two LSBs and masking the MSB
    int etaAbsRed = (std::abs(mu->hwEta()) >> 2) & ((1 << 6) - 1);
    
    int deltaPhi = 0;
    int deltaEta = 0;

    if (mu->hwPt() < 64) { // extrapolation only for "low" pT muons
      int sign = 1;
      if (mu->hwSign() == 0) {
        sign = -1;
      }
      deltaPhi = (m_phiExtrapolationLUTs.at(mu->type())->lookup(etaAbsRed, ptRed) << 3) * sign;
      deltaEta = (m_etaExtrapolationLUTs.at(mu->type())->lookup(etaAbsRed, ptRed) << 3);
    }

    mu->setExtrapolation(deltaEta, deltaPhi);
  }
}

void
l1t::MicroGMTIsolationUnit::calculate5by1Sums(const MicroGMTConfiguration::CaloInputCollection& inputs) 
{
  m_5by1TowerSums.clear();
  if (inputs.size() == 0) return;

  for (int iphi = 0; iphi < 36; ++iphi) {
    int iphiIndexOffset = iphi*28;
    m_5by1TowerSums.push_back(inputs[iphiIndexOffset].etBits()+inputs[iphiIndexOffset+1].etBits()+inputs[iphiIndexOffset+2].etBits());//ieta = 0 (tower -28)
    m_5by1TowerSums.push_back(inputs[iphiIndexOffset-1].etBits()+inputs[iphiIndexOffset].etBits()+inputs[iphiIndexOffset+1].etBits()+inputs[iphiIndexOffset+2].etBits()); // 
    for (int ieta = 2; ieta < 26; ++ieta) {
      int sum = 0;
      for (int dIEta = -2; dIEta <= 2; ++dIEta) {
        sum += inputs[iphiIndexOffset+dIEta].etBits();
      }
      m_5by1TowerSums.push_back(sum);
    }
    m_5by1TowerSums.push_back(inputs[iphiIndexOffset+1].etBits()+inputs[iphiIndexOffset].etBits()+inputs[iphiIndexOffset-1].etBits()+inputs[iphiIndexOffset-2].etBits());
    m_5by1TowerSums.push_back(inputs[iphiIndexOffset].etBits()+inputs[iphiIndexOffset-1].etBits()+inputs[iphiIndexOffset-2].etBits());//ieta = 0 (tower 28)
  }

  m_initialSums = true;
}


int 
l1t::MicroGMTIsolationUnit::calculate5by5Sum(unsigned index) const
{
  if (index > m_5by1TowerSums.size()) {
    edm::LogWarning("energysum out of bounds!");
    return 0;
  }
  // phi wrap around:
  int returnSum = 0; 
  for (int dIPhi = -2; dIPhi <= 2; ++dIPhi) {
    int currIndex = (index + dIPhi*28)%1008; // wrap-around at top
    if (currIndex < 0) currIndex = 1008+currIndex;
    if ((unsigned)currIndex < m_5by1TowerSums.size()) {
      returnSum += m_5by1TowerSums[currIndex];
    } else {
      edm::LogWarning("energysum out of bounds!");
    }
  } 
  return std::min(31, returnSum);
}

void 
l1t::MicroGMTIsolationUnit::isolate(MicroGMTConfiguration::InterMuonList& muons) const
{
  MicroGMTConfiguration::InterMuonList::iterator muIt;
  for (muIt = muons.begin(); muIt != muons.end(); ++muIt) {
    int caloIndex = getCaloIndex(*muIt);
    int energySum = calculate5by5Sum(caloIndex);
    muIt->setHwIsoSum(energySum);

    int absIso = m_AbsIsoCheckMem.lookup(energySum);
    int relIso = m_RelIsoCheckMem.lookup(energySum, muIt->hwPt());

    muIt->setHwRelIso(relIso);
    muIt->setHwAbsIso(absIso);
  }
}

void l1t::MicroGMTIsolationUnit::setTowerSums(const MicroGMTConfiguration::CaloInputCollection& inputs) {
  m_towerEnergies.clear();
  if (inputs.size() == 0) return;
  for (const auto& input:inputs) {
    if ( input.etBits() != 0 ) {
      m_towerEnergies[input.hwEta()*36+input.hwPhi()] = input.etBits();
    }
  }

  m_initialSums = true;

} 

void l1t::MicroGMTIsolationUnit::isolatePreSummed(MicroGMTConfiguration::InterMuonList& muons) const 
{
  MicroGMTConfiguration::InterMuonList::iterator muIt;
  for (muIt = muons.begin(); muIt != muons.end(); ++muIt) {
    int caloIndex = getCaloIndex(*muIt);
    int energySum = 0;
    if (m_towerEnergies.count(caloIndex) == 1) {
      energySum = m_towerEnergies.at(caloIndex);
    }

    muIt->setHwIsoSum(energySum);

    int absIso = m_AbsIsoCheckMem.lookup(energySum);
    int relIso = m_RelIsoCheckMem.lookup(energySum, muIt->hwPt());

    muIt->setHwRelIso(relIso);
    muIt->setHwAbsIso(absIso);
  }
  
}
