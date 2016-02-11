// -*- C++ -*-
//
// Class:      L1TMuonEndcapParamsESProducer
//
// Original Author:  Matthew Carver UF
//         Created:
//
//


// system include files
#include <memory>
#include "boost/shared_ptr.hpp"

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"

#include "CondFormats/L1TObjects/interface/L1TMuonEndcapParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonEndcapParamsRcd.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

// class declaration

class L1TMuonEndcapParamsESProducer : public edm::ESProducer {
   public:
      L1TMuonEndcapParamsESProducer(const edm::ParameterSet&);
      ~L1TMuonEndcapParamsESProducer();
   
      typedef boost::shared_ptr<L1TMuonEndcapParams> ReturnType;

      ReturnType produce(const L1TMuonEndcapParamsRcd&);
   private:
      L1TMuonEndcapParams m_params;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1TMuonEndcapParamsESProducer::L1TMuonEndcapParamsESProducer(const edm::ParameterSet& iConfig)
{
   //the following line is needed to tell the framework what
   // data is being produced
   setWhatProduced(this);
   // Firmware version
   unsigned fwVersion = iConfig.getParameter<unsigned>("PT_assignment_version");
   m_params.SetVersion(fwVersion);
   
   int PBits = iConfig.getParameter<int>("NumPhiBits");
   m_params.SetNumPhiBits(PBits);


}


L1TMuonEndcapParamsESProducer::~L1TMuonEndcapParamsESProducer()
{
}



//
// member functions
//

// ------------ method called to produce the data  ------------
L1TMuonEndcapParamsESProducer::ReturnType
L1TMuonEndcapParamsESProducer::produce(const L1TMuonEndcapParamsRcd& iRecord)
{
   using namespace edm::es;
   boost::shared_ptr<L1TMuonEndcapParams> pEMTFParams;

   pEMTFParams = boost::shared_ptr<L1TMuonEndcapParams>(new L1TMuonEndcapParams(m_params));
   return pEMTFParams;
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(L1TMuonEndcapParamsESProducer);
