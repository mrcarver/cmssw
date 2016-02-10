#ifndef L1MuEMTFParametersRCD_H
#define L1MuEMTFParametersRCD_H

#include "boost/mpl/vector.hpp"

#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "CondFormats/DataRecord/interface/L1TriggerKeyListRcd.h"
#include "CondFormats/DataRecord/interface/L1TriggerKeyRcd.h"

class L1MuEMTFParametersRcd : public
edm::eventsetup::DependentRecordImplementation<L1MuEMTFParametersRcd, boost::mpl::vector<L1TriggerKeyListRcd,L1TriggerKeyRcd> > {};

#endif
