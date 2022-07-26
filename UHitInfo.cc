//____________________________________________________________________
//
// UHitInfo is a data class for storing geant4 hit information
//
//
#include "UHitInfo.h"
#include "TMath.h"
#ifndef __IOSTREAM__
#include <iostream>
#endif
using std::endl;

//____________________________________________________________________
// ClassImp(UHitInfo);

//____________________________________________________________________
UHitInfo::UHitInfo() {
  // we assume the there is no default value allowed for proper operation
  fEntries = 0;
}
//____________________________________________________________________
UHitInfo::~UHitInfo() {}

//____________________________________________________________________
ostream& operator<<(ostream& os, UHitInfo* hitinfo) {
  os << " Entries=" << hitinfo->GetEntries() << endl;
  return os;
}
