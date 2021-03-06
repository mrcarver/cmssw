#ifndef SimG4Core_OscarMTProducer_H
#define SimG4Core_OscarMTProducer_H

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"

#include "SimG4Core/Application/interface/OscarMTMasterThread.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <memory>

class SimProducer;
class RunManagerMTWorker;

class OscarMTProducer : public edm::stream::EDProducer<
  edm::GlobalCache<OscarMTMasterThread>,
  edm::RunCache<int> // for some reason void doesn't compile
>
{
public:
  typedef std::vector<boost::shared_ptr<SimProducer> > Producers;

  explicit OscarMTProducer(edm::ParameterSet const & p, const OscarMTMasterThread *);
  virtual ~OscarMTProducer();

  static std::unique_ptr<OscarMTMasterThread> initializeGlobalCache(const edm::ParameterSet& iConfig);
  static std::shared_ptr<int> globalBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup, const OscarMTMasterThread *masterThread);
  static void globalEndRun(const edm::Run& iRun, const edm::EventSetup& iSetup, const RunContext *iContext);
  static void globalEndJob(OscarMTMasterThread *masterThread);


  virtual void endRun(const edm::Run & r,const edm::EventSetup& c) override;
  virtual void produce(edm::Event & e, const edm::EventSetup& c) override;

private:
  Producers     m_producers;
  std::unique_ptr<RunManagerMTWorker> m_runManagerWorker;
  //edm::EDGetTokenT<edm::HepMCProduct> m_HepMC;
};

#endif
