//____________________________________________________________________
//
// SensorHit is a data class for storing mapped data for
// one QDC detector
//

//
// $Id: SensorHit.cpp,v 1.1 2012/02/09 23:51:09 dc Exp $
//
#include "USensorHit.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif
using std::endl;
using std::cout;

//____________________________________________________________________
// ClassImp(USensorHit);

//____________________________________________________________________
USensorHit::USensorHit() {
  fClusterSize = 0;
  fDevX = 0;
  fDevY = 0;
  fDevZ = 0;
  fBelongsToPrimaryTrack = false;
  fUsed = false;
}

//____________________________________________________________________
USensorHit::USensorHit(const double x, const double y, const double z, const int clustsize) {
  fX = x;
  fY = y;
  fZ = z;
  fClusterSize = clustsize;
  fBelongsToPrimaryTrack = false;
  fUsed = false;

  fLocalX = x;
  fLocalY = y;
}

//____________________________________________________________________
USensorHit::~USensorHit() {}

//____________________________________________________________________
ostream& operator<<(ostream& os, USensorHit* hit) {
  os << " Hit(x,y,z):"
     << ", " << hit->GetX() << ", " << hit->GetY() << ", " << hit->GetZ() << endl;
  return os;
}
