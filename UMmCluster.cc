//____________________________________________________________________
//
// UMmCluster is a data class for storing geant4 hit information
//
//
#include "UMmCluster.h"
#include "TMath.h"
#ifndef __IOSTREAM__
#include <iostream>
#endif
using std::endl;

//____________________________________________________________________
// ClassImp(UMmCluster);

//____________________________________________________________________
UMmCluster::UMmCluster() {
  // we assume the there is no default value allowed for proper operation
  fNstrips = 0;
}
//____________________________________________________________________
UMmCluster::~UMmCluster() {}

//____________________________________________________________________
ostream& operator<<(ostream& os, UMmCluster* clust) {
  os << " MaxNStrip=" << clust->GetNstrips() << endl;
  return os;
}
