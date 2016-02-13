///
/// \class L1TMuonEndcapParams
///
/// Description: Placeholder for EMTF parameters
///
///
/// \author: Matthew Carver
///

#ifndef L1TEMTFParams_h
#define L1TEMTFParams_h

#include <memory>
#include <iostream>
#include <vector>

#include "CondFormats/Serialization/interface/Serializable.h"
#include "CondFormats/L1TObjects/interface/LUT.h"

class L1TMuonEndcapParams {

public:

 	void SetNumPhiBits(int phiBits){NumPhiBits_ = phiBits;};
	void SetVersion(unsigned version){version_ = version;};
	
	int GetPhiBits() const {return NumPhiBits_;};
	unsigned GetVersion(){return version_;};
	
	
	L1TMuonEndcapParams() { version_=1; }
   ~L1TMuonEndcapParams() {}

	void print(std::ostream&) const;

private:
	
  unsigned version_;
  int NumPhiBits_;

  COND_SERIALIZABLE;
};
#endif
