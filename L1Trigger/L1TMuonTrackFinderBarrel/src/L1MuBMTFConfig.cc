//-------------------------------------------------
//
//   Class: L1MuBMTFConfig
//
//   Description: BMTrackFinder parameters for L1MuBMTrackFinder
//
//
//
//   Author :
//   N. Neumeister            CERN EP
//   J. Troconiz              UAM Madrid
//
//--------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------

#include "L1Trigger/L1TMuonTrackFinderBarrel/src/L1MuBMTFConfig.h"

//---------------
// C++ Headers --
//---------------

#include <iostream>
#include <string>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

using namespace std;

// --------------------------------
//       class L1MuBMTFConfig
//---------------------------------

//----------------
// Constructors --
//----------------

L1MuBMTFConfig::L1MuBMTFConfig(const edm::ParameterSet & ps) {

    m_ps = &ps;
    setDefaults();

}


//--------------
// Destructor --
//--------------
L1MuBMTFConfig::~L1MuBMTFConfig() {}


//--------------
// Operations --
//--------------

void L1MuBMTFConfig::setDefaults() {

  m_BMDigiInputTag = m_ps->getParameter<edm::InputTag>("DTDigi_Source");
  //m_CSCTrSInputTag = m_ps->getParameter<edm::InputTag>("CSCStub_Source");

  m_debug = true;
  m_dbgLevel = m_ps->getUntrackedParameter<int>("Debug",0);

  //m_overlap = m_ps->getUntrackedParameter<bool>("Overlap",true);

  // set min and max bunch crossing
  m_BxMin = m_ps->getUntrackedParameter<int>("BX_min",-9);
  m_BxMax = m_ps->getUntrackedParameter<int>("BX_max", 7);

  // set Filter for Extrapolator
  m_extTSFilter = m_ps->getUntrackedParameter<int>("Extrapolation_Filter",1);

  // set switch for open LUTs usage
  m_openLUTs = m_ps->getUntrackedParameter<bool>("Open_LUTs",false);

  // set switch for EX21 usage
  m_useEX21 = m_ps->getUntrackedParameter<bool>("Extrapolation_21",false);

  // set switch for eta track finder usage
  m_etaTF = m_ps->getUntrackedParameter<bool>("EtaTrackFinder",true);

  // set switch for etaFlag cancellation of CSC segments
  //m_etacanc = m_ps->getUntrackedParameter<bool>("CSC_Eta_Cancellation",false);

  // set Filter for Out-of-time Track Segments
  m_TSOutOfTimeFilter = m_ps->getUntrackedParameter<bool>("OutOfTime_Filter",false);
  m_TSOutOfTimeWindow = m_ps->getUntrackedParameter<int>("OutOfTime_Filter_Window",1);

  // set precision for extrapolation
  m_NbitsExtPhi  = m_ps->getUntrackedParameter<int>("Extrapolation_nbits_Phi", 8);
  m_NbitsExtPhib = m_ps->getUntrackedParameter<int>("Extrapolation_nbits_PhiB",8);

  // set precision for pt-assignment
  m_NbitsPtaPhi  = m_ps->getUntrackedParameter<int>("PT_Assignment_nbits_Phi", 12);
  m_NbitsPtaPhib = m_ps->getUntrackedParameter<int>("PT_Assignment_nbits_PhiB",10);

  // set precision for phi-assignment look-up tables
  m_NbitsPhiPhi  = m_ps->getUntrackedParameter<int>("PHI_Assignment_nbits_Phi", 10);
  m_NbitsPhiPhib = m_ps->getUntrackedParameter<int>("PHI_Assignment_nbits_PhiB",10);

  if ( Debug(1) ) cout << endl;
  if ( Debug(1) ) cout << "*******************************************" << endl;
  if ( Debug(1) ) cout << "**** L1 barrel Track Finder settings : ****" << endl;
  if ( Debug(1) ) cout << "*******************************************" << endl;
  if ( Debug(1) ) cout << endl;

  if ( Debug(1) ) cout << "L1 barrel Track Finder : BM Digi Source:  " <<  m_BMDigiInputTag << endl;
  //if ( Debug(1) ) cout << "L1 barrel Track Finder : CSC Stub Source: " <<  m_CSCTrSInputTag << endl;
  if ( Debug(1) ) cout << endl;

  if ( Debug(1) ) cout << "L1 barrel Track Finder : debug level: " << m_dbgLevel << endl;

  if ( Debug(1) ) cout << "L1 barrel Track Finder : minimal bunch-crossing : " << m_BxMin << endl;
  if ( Debug(1) ) cout << "L1 barrel Track Finder : maximal bunch-crossing : " << m_BxMax << endl;

  if ( Debug(1) ) cout << "L1 barrel Track Finder : Extrapolation Filter : " << m_extTSFilter << endl;

  if ( Debug(1) && m_openLUTs) {
    cout << "L1 barrel Track Finder : use open LUTs : on" << endl;
  }
  if ( Debug(1) && !m_openLUTs) {
    cout << "L1 barrel Track Finder : use open LUTs : off" << endl;
  }

  if ( Debug(1) && m_useEX21 ) {
    cout << "L1 barrel Track Finder : use EX21 extrapolations : on" << endl;
  }
  if ( Debug(1) && !m_useEX21 ) {
    cout << "L1 barrel Track Finder : use EX21 extrapolations : off" << endl;
  }

  if ( Debug(1) && m_etaTF ) {
    cout << "L1 barrel Track Finder : Eta Track Finder : on" << endl;
  }
  if ( Debug(1) && !m_etaTF ) {
    cout << "L1 barrel Track Finder : Eta Track Finder : off" << endl;
  }

  //if ( Debug(1) && m_etacanc ) {
  //  cout << "L1 barrel Track Finder : CSC etaFlag cancellation : on" << endl;
  //}
  //if ( Debug(1) && !m_etacanc ) {
  //  cout << "L1 barrel Track Finder : CSC etaFlag cancellation : off" << endl;
  //}

  if ( Debug(1) && m_TSOutOfTimeFilter ) {
    cout << "L1 barrel Track Finder : out-of-time TS filter : on" << endl;
    cout << "L1 barrel Track Finder : out-of-time TS filter window : " << m_TSOutOfTimeWindow << endl;
  }
  if ( Debug(1) && !m_TSOutOfTimeFilter ) {
    cout << "L1 barrel Track Finder : out-of-time TS filter : off" << endl;
  }

  if ( Debug(1) ) cout << "L1 barrel Track Finder : # of bits used for phi  (extrapolation)  : " << m_NbitsExtPhi << endl;
  if ( Debug(1) ) cout << "L1 barrel Track Finder : # of bits used for phib (extrapolation)  : " << m_NbitsExtPhib << endl;
  if ( Debug(1) ) cout << "L1 barrel Track Finder : # of bits used for phi  (pt-assignment)  : " << m_NbitsPtaPhi << endl;
  if ( Debug(1) ) cout << "L1 barrel Track Finder : # of bits used for phib (pt-assignment)  : " << m_NbitsPtaPhib << endl;
  if ( Debug(1) ) cout << "L1 barrel Track Finder : # of bits used for phi  (phi-assignment) : " << m_NbitsPhiPhi << endl;
  if ( Debug(1) ) cout << "L1 barrel Track Finder : # of bits used for phib (phi-assignment) : " << m_NbitsPhiPhib << endl;

}


// static data members

edm::InputTag L1MuBMTFConfig::m_BMDigiInputTag = edm::InputTag();
//edm::InputTag L1MuBMTFConfig::m_CSCTrSInputTag = edm::InputTag();

bool L1MuBMTFConfig::m_debug = false;
int  L1MuBMTFConfig::m_dbgLevel = -1;
int  L1MuBMTFConfig::m_BxMin = -9;
int  L1MuBMTFConfig::m_BxMax =  7;
int  L1MuBMTFConfig::m_extTSFilter  = 1;
bool L1MuBMTFConfig::m_openLUTs  = false;
bool L1MuBMTFConfig::m_useEX21 = false;
bool L1MuBMTFConfig::m_etaTF = true;
//bool L1MuBMTFConfig::m_etacanc = false;
bool L1MuBMTFConfig::m_TSOutOfTimeFilter = false;
int  L1MuBMTFConfig::m_TSOutOfTimeWindow = 1;
int  L1MuBMTFConfig::m_NbitsExtPhi  = 8;
int  L1MuBMTFConfig::m_NbitsExtPhib = 8;
int  L1MuBMTFConfig::m_NbitsPtaPhi  = 12;
int  L1MuBMTFConfig::m_NbitsPtaPhib = 10;
int  L1MuBMTFConfig::m_NbitsPhiPhi  = 10;
int  L1MuBMTFConfig::m_NbitsPhiPhib = 10;
