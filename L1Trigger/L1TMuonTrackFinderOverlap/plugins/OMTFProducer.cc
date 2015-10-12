#include <iostream>
#include <strstream>
#include <vector>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"

#include "L1Trigger/L1TMuonTrackFinderOverlap/plugins/OMTFProducer.h"
#include "L1Trigger/L1TMuonTrackFinderOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonTrackFinderOverlap/interface/OMTFinputMaker.h"
#include "L1Trigger/L1TMuonTrackFinderOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonTrackFinderOverlap/interface/OMTFSorter.h"
#include "L1Trigger/L1TMuonTrackFinderOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonTrackFinderOverlap/interface/XMLConfigWriter.h"

using namespace L1TMuon;

OMTFProducer::OMTFProducer(const edm::ParameterSet& cfg):
  theConfig(cfg),
  trigPrimSrc(cfg.getParameter<edm::InputTag>("TriggerPrimitiveSrc")){

  produces<l1t::RegionalMuonCandBxCollection >("OMTF");

  inputToken = consumes<TriggerPrimitiveCollection>(trigPrimSrc);

  if(!theConfig.exists("omtf")){
    edm::LogError("OMTFProducer")<<"omtf configuration not found in cfg.py";
  }

  myInputMaker = new OMTFinputMaker();
  mySorter = new OMTFSorter();
  myWriter = 0;

  dumpResultToXML = theConfig.getParameter<bool>("dumpResultToXML");
  dumpDetailedResultToXML = theConfig.getParameter<bool>("dumpDetailedResultToXML");
  dumpGPToXML = theConfig.getParameter<bool>("dumpGPToXML");
  theConfig.getParameter<std::string>("XMLDumpFileName");

  if(dumpResultToXML || dumpGPToXML){
    myWriter = new XMLConfigWriter();
    std::string fName = "OMTF_Events";
    myWriter->initialiseXMLDocument(fName);
  }

  myOMTFConfig = 0;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
OMTFProducer::~OMTFProducer(){

  delete myOMTFConfig;
  delete myOMTF;

  delete myInputMaker;
  delete mySorter;

  if (myWriter) delete myWriter;

}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFProducer::beginJob(){

  if(theConfig.exists("omtf")){
    myOMTFConfig = new OMTFConfiguration(theConfig.getParameter<edm::ParameterSet>("omtf"));
    myOMTF = new OMTFProcessor(theConfig.getParameter<edm::ParameterSet>("omtf"));
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFProducer::endJob(){

  if(dumpResultToXML && !dumpGPToXML){
    std::string fName = theConfig.getParameter<std::string>("XMLDumpFileName");
    myWriter->finaliseXMLDocument(fName);
  }

  if(dumpGPToXML && !dumpResultToXML){

    GoldenPattern *dummy = new GoldenPattern(Key(0,0,0));
    dummy->reset();

    std::string fName = "OMTF";
    myWriter->initialiseXMLDocument(fName);
    const std::map<Key,GoldenPattern*> & myGPmap = myOMTF->getPatterns();
    for(auto itGP: myGPmap){
      //std::cout<<*itGP.second<<std::endl;
      //myWriter->writeGPData(*itGP.second);
      if(itGP.second->key().thePtCode>5) myWriter->writeGPData(*itGP.second,*dummy, *dummy, *dummy);
    }
    fName = "GPs.xml";
    myWriter->finaliseXMLDocument(fName);
    ///Write GPs merged by 4 above iPt19, and by 2 below//
    //////////////////////////////////////////////////////
    ///4x merging
    fName = "OMTF";
    myWriter->initialiseXMLDocument(fName);
    myOMTF->averagePatterns(1);
    myOMTF->averagePatterns(-1);

    writeMergedGPs();

    fName = "GPs_4x.xml";
    myWriter->finaliseXMLDocument(fName);
  }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFProducer::writeMergedGPs(){
  
  const std::map<Key,GoldenPattern*> & myGPmap = myOMTF->getPatterns();

  GoldenPattern *dummy = new GoldenPattern(Key(0,0,0));
  dummy->reset();

  unsigned int iPtMin = 6;
  Key aKey = Key(1, iPtMin,-1);
  while(myGPmap.find(aKey)!=myGPmap.end()){

    GoldenPattern *aGP1 = myGPmap.find(aKey)->second;
    GoldenPattern *aGP2 = dummy;
    GoldenPattern *aGP3 = dummy;
    GoldenPattern *aGP4 = dummy;

    std::cout<<"A "<<aGP1->key()<<std::endl;

    ++aKey.thePtCode;
    if(aKey.thePtCode<=31 && myGPmap.find(aKey)!=myGPmap.end()) aGP2 =  myGPmap.find(aKey)->second;

    if(aKey.thePtCode>19){
      ++aKey.thePtCode;
      if(aKey.thePtCode<=31 && myGPmap.find(aKey)!=myGPmap.end()) aGP3 =  myGPmap.find(aKey)->second;

      ++aKey.thePtCode;
      if(aKey.thePtCode<=31 && myGPmap.find(aKey)!=myGPmap.end()) aGP4 =  myGPmap.find(aKey)->second;
    }
    ++aKey.thePtCode;
    myWriter->writeGPData(*aGP1,*aGP2, *aGP3, *aGP4);

    ///Write the opposite charge.
    Key aTmpKey = aGP1->key();
    std::cout<<aTmpKey<<std::endl;
    aTmpKey.theCharge = 1;
    if(myGPmap.find(aTmpKey)!=myGPmap.end()) aGP1 =  myGPmap.find(aTmpKey)->second;
    else aGP1 = dummy;

    aTmpKey = aGP2->key();
    aTmpKey.theCharge = 1;
    if(myGPmap.find(aTmpKey)!=myGPmap.end()) aGP2 =  myGPmap.find(aTmpKey)->second;
    else aGP2 = dummy;

    aTmpKey = aGP3->key();
    aTmpKey.theCharge = 1;
    if(myGPmap.find(aTmpKey)!=myGPmap.end()) aGP3 =  myGPmap.find(aTmpKey)->second;
    else aGP3 = dummy;

    aTmpKey = aGP4->key();
    aTmpKey.theCharge = 1;
    if(myGPmap.find(aTmpKey)!=myGPmap.end()) aGP4 =  myGPmap.find(aTmpKey)->second;
    else aGP4 = dummy;
    
    myWriter->writeGPData(*aGP1,*aGP2, *aGP3, *aGP4);
    }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFProducer::produce(edm::Event& iEvent, const edm::EventSetup& evSetup){

  std::ostringstream myStr;

  myInputMaker->initialize(evSetup);

  edm::Handle<TriggerPrimitiveCollection> trigPrimitives;
  iEvent.getByToken(inputToken, trigPrimitives);

  ///Filter digis by dropping digis from selected (by cfg.py) subsystems
  const L1TMuon::TriggerPrimitiveCollection filteredDigis = filterDigis(*trigPrimitives);

  std::auto_ptr<l1t::RegionalMuonCandBxCollection > myCands(new l1t::RegionalMuonCandBxCollection);

  if(dumpResultToXML) aTopElement = myWriter->writeEventHeader(iEvent.id().event());

  //l1t::tftype mtfType = l1t::tftype::bmtf;
  //l1t::tftype mtfType = l1t::tftype::omtf_pos;
  //l1t::tftype mtfType = l1t::tftype::omtf_neg;
  //l1t::tftype mtfType = l1t::tftype::emtf_pos;

  // NOTE: assuming all is for bx 0
  int bx = 0;

  ///Loop over all processors, each covering 60 deg in phi
  for(unsigned int iProcessor=0;iProcessor<6;++iProcessor){

    myStr<<" iProcessor: "<<iProcessor;

    ///Input data with phi ranges shifted for each processor, so it fits 11 bits range
    const OMTFinput *myInputPos = myInputMaker->buildInputForProcessor(filteredDigis,iProcessor, l1t::tftype::omtf_pos);
    OMTFinput myShiftedInputPos =  myOMTF->shiftInput(iProcessor,*myInputPos);

    const OMTFinput *myInputNeg = myInputMaker->buildInputForProcessor(filteredDigis,iProcessor, l1t::tftype::omtf_neg);
    OMTFinput myShiftedInputNeg =  myOMTF->shiftInput(iProcessor,*myInputNeg);

    l1t::RegionalMuonCandBxCollection myOTFCandidatesPos, myOTFCandidatesNeg;
    ///Results for each GP in each logic region of given processor
    //Retreive all candidates returned by sorter: upto 3 non empty ones with different phi or charge
    const std::vector<OMTFProcessor::resultsMap> & myResultsNeg = myOMTF->processInput(iProcessor,myShiftedInputNeg);
    mySorter->sortProcessor(myResultsNeg, myOTFCandidatesNeg, bx);

    const std::vector<OMTFProcessor::resultsMap> & myResultsPos = myOMTF->processInput(iProcessor,myShiftedInputPos);
    mySorter->sortProcessor(myResultsPos, myOTFCandidatesPos, bx);

    ///Shift phi scales, and put MicroGMT candidates into myCands collection
    processCandidates(iProcessor, bx, myCands, myOTFCandidatesPos, l1t::tftype::omtf_pos);
    processCandidates(iProcessor, bx, myCands, myOTFCandidatesNeg, l1t::tftype::omtf_neg);

    ///Write data to XML file
    if(dumpResultToXML){
      xercesc::DOMElement * aProcElement = myWriter->writeEventData(aTopElement,iProcessor,myShiftedInputPos);
      for(unsigned int iRefHit=0;iRefHit<OMTFConfiguration::nTestRefHits;++iRefHit){
	///Dump only regions, where a candidate was found
	InternalObj myCand = mySorter->sortRefHitResults(myResultsPos[iRefHit],0);//charge=0 means ignore charge
	if(myCand.pt){
	  myWriter->writeCandidateData(aProcElement,iRefHit,myCand);
	  if(dumpDetailedResultToXML){
	    for(auto & itKey: myResultsNeg[iRefHit]) myWriter->writeResultsData(aProcElement,
										iRefHit,
										itKey.first,itKey.second);
	  }
	}
      }
    }
  }

  //dumpResultToXML = true;
  myStr<<" Number of candidates: "<<myCands->size(bx);
  edm::LogInfo("OMTFOMTFProducer")<<myStr.str();

  iEvent.put(myCands, "OMTF");
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void OMTFProducer::processCandidates(unsigned int iProcessor, int bx,
				     std::auto_ptr<l1t::RegionalMuonCandBxCollection > & myCands,
				     l1t::RegionalMuonCandBxCollection & myOTFCandidates,
				     l1t::tftype mtfType){

  ////Switch from internal processor n bit scale to global one
  int procOffset = OMTFConfiguration::globalPhiStart(iProcessor);
  if(procOffset<0) procOffset+=(int)OMTFConfiguration::nPhiBins;
  ///Set local 0 at iProcessor x 15 deg
  procOffset-=(15+iProcessor*60)/360.0*OMTFConfiguration::nPhiBins;
  int lowScaleEnd = pow(2,OMTFConfiguration::nPhiBits-1);

    for(unsigned int iCand=0; iCand<myOTFCandidates.size(bx); ++iCand){
      // shift phi from processor to global coordinates
      l1t::RegionalMuonCand cand = myOTFCandidates.at(bx, iCand);
      int phiValue = (cand.hwPhi()+procOffset+lowScaleEnd);
      if(phiValue>=(int)OMTFConfiguration::nPhiBins) phiValue-=OMTFConfiguration::nPhiBins;
      phiValue/=10; //MicroGMT has 10x coarser scale than OMTF

      cand.setHwPhi(phiValue);
      cand.setTFIdentifiers(iProcessor,mtfType);
      // store candidate
      if(cand.hwPt()) myCands->push_back(bx, cand);
    }
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
const L1TMuon::TriggerPrimitiveCollection OMTFProducer::filterDigis(const L1TMuon::TriggerPrimitiveCollection & vDigi){

  if(!theConfig.getParameter<bool>("dropRPCPrimitives") &&
     !theConfig.getParameter<bool>("dropDTPrimitives") &&
     !theConfig.getParameter<bool>("dropCSCPrimitives")) return vDigi;

  L1TMuon::TriggerPrimitiveCollection filteredDigis;
  for(auto it:vDigi){
    switch (it.subsystem()) {
    case L1TMuon::TriggerPrimitive::kRPC: {
      if(!theConfig.getParameter<bool>("dropRPCPrimitives")) filteredDigis.push_back(it);
      break;
    }
    case L1TMuon::TriggerPrimitive::kDT: {
      if(!theConfig.getParameter<bool>("dropDTPrimitives")) filteredDigis.push_back(it);
      break;
    }
    case L1TMuon::TriggerPrimitive::kCSC: {
      if(!theConfig.getParameter<bool>("dropCSCPrimitives")) filteredDigis.push_back(it);
      break;
    }
    case L1TMuon::TriggerPrimitive::kNSubsystems: {break;}
    }
  }
  return filteredDigis;
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
