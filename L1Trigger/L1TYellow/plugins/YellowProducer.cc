// -*- C++ -*-
//
// Package:    YellowProducer
// Class:      YellowProducer
// 
/**\class YellowProducer YellowProducer.cc L1Trigger/L1TYellow/src/YellowProducer.cc

 Description: Emulation of Fictitious Level-1 Yellow Trigger for demonstration purposes

 Implementation:
     [Notes on implementation]
*/
//


// system include files
#include <boost/shared_ptr.hpp>

// user include files

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/L1TYellowParamsRcd.h"
#include "CondFormats/L1TYellow/interface/YellowParams.h"
#include "DataFormats/L1TYellow/interface/YellowDigi.h"
#include "DataFormats/L1TYellow/interface/YellowOutput.h"
#include "L1Trigger/L1TYellow/interface/YellowFirmware.h"
#include "L1Trigger/L1TYellow/interface/YellowFirmwareFactory.h"



using namespace std;
using namespace edm;

namespace l1t {

//
// class declaration
//
  
  class YellowProducer : public EDProducer {
  public:
    explicit YellowProducer(const ParameterSet&);
    ~YellowProducer();
    
    static void fillDescriptions(ConfigurationDescriptions& descriptions);
    
  private:
    virtual void produce(Event&, EventSetup const&);
    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(Run const&iR, EventSetup const&iE);
    virtual void endRun(Run const& iR, EventSetup const& iE);
  
    // ----------member data ---------------------------
    unsigned long long m_yellowParamsCacheId; // Cache-ID from current parameters, to check if needs to be updated.
    boost::shared_ptr<const YellowParams>  m_dbpars; // Database parameters for the trigger, to be updated as needed.
    boost::shared_ptr<YellowFirmware> m_fw; // Algorithm to run per event, depends on database parameters.

    YellowFirmwareFactory m_factory; // Factory to produce algorithms based on DB parameters

    EDGetToken yellowDigisToken;
  };
  
  //
  // constructors and destructor
  //
  YellowProducer::YellowProducer(const ParameterSet& iConfig)
  {
    // register what you produce
    produces<YellowOutputCollection>();
    
    // register what you consume and keep token for later access:
    yellowDigisToken = consumes<YellowDigiCollection>(iConfig.getParameter<InputTag>("fakeRawToDigi"));
    
    // set cache id to zero, will be set at first beginRun:
    m_yellowParamsCacheId = 0;
  }


  YellowProducer::~YellowProducer()
  {
  }


//
// member functions
//

// ------------ method called to produce the data  ------------
void
YellowProducer::produce(Event& iEvent, const EventSetup& iSetup)
{

  LogInfo("l1t|yellow") << "YellowProducer::produce function called...\n";

  cout << "cout version: YellowProducer::produce function called...\n";

  
  Handle<YellowDigiCollection> inputDigis;
  iEvent.getByToken(yellowDigisToken,inputDigis);
  
  std::auto_ptr<YellowOutputCollection> outColl (new YellowOutputCollection);
  YellowOutput iout;

  if (inputDigis->size()){
    if (m_fw) {
      m_fw->processEvent(*inputDigis, *outColl);
    } else {
      cout << "firmware is invalid...\n";
    }
    // already complained in beginRun, doing nothing now will send empty collection to event, as desired.
  } else {
    LogError("l1t|yellow") << "YellowProducer: input Digis have zero size.\n";
  }

  //iout.setRawData(20);  outColl->push_back(iout);
  //iout.setRawData(18);  outColl->push_back(iout);
  //iout.setRawData(30);  outColl->push_back(iout);

  iEvent.put(outColl);
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
YellowProducer::beginJob()
{
  cout << "YellowProducer begin JOB being called...\n";

}

// ------------ method called once each job just after ending the event loop  ------------
void 
YellowProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------

void YellowProducer::beginRun(Run const&iR, EventSetup const&iE){

  cout << "YellowProducer Begin Run Called!\n";

  unsigned long long id = iE.get<L1TYellowParamsRcd>().cacheIdentifier();

  if (id != m_yellowParamsCacheId){ // Need to update:
    m_yellowParamsCacheId = id;

    // Retrieve the  yellow parameters from the event setup:
    ESHandle<YellowParams> yParameters;
    iE.get<L1TYellowParamsRcd>().get(yParameters);

    m_dbpars = boost::shared_ptr<const YellowParams>(yParameters.product());
  
    if (m_dbpars){
      cout << "dbpars is non-null...\n";
    } else {
      cout << "dbpars is null...\n";
    }

    // Set the current algorithm version based on DB pars from database:
    m_fw = m_factory.create(*m_dbpars);

    if (! m_fw) {
      // we complain here once per run
      LogError("l1t|yellow") << "YellowProducer:  could not retreive DB params from Event Setup\n";
    }
  }



}

// ------------ method called when ending the processing of a run  ------------
void YellowProducer::endRun(Run const& iR, EventSetup const& iE){
  
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
YellowProducer::fillDescriptions(ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

} // namespace

//define this as a plug-in
DEFINE_FWK_MODULE(l1t::YellowProducer);


// OR?:
////define this as a plug-in
//DEFINE_FWK_MODULE(YellowProducer);
//} // namespace

