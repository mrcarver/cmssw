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
   unsigned Version = iConfig.getParameter<unsigned>("PT_assignment_version");
   m_params.SetPtAssignVersion(Version);
   
   int S1MatchWindow = iConfig.getParameter<int>("St1PhiMatchWindow");
   m_params.SetSt1PhiMatchWindow(S1MatchWindow);
   
   int S2MatchWindow = iConfig.getParameter<int>("St2PhiMatchWindow");
   m_params.SetSt2PhiMatchWindow(S2MatchWindow);
   
   int S3MatchWindow = iConfig.getParameter<int>("St3PhiMatchWindow");
   m_params.SetSt3PhiMatchWindow(S3MatchWindow);
   
   int S4MatchWindow = iConfig.getParameter<int>("St4PhiMatchWindow");
   m_params.SetSt4PhiMatchWindow(S4MatchWindow);
   
   std::string XmlPtDir = iConfig.getParameter<std::string>("XmlPtDir");
   m_params.SetXmlPtTreeDir(XmlPtDir);


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
