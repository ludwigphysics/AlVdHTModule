//____________________________________________________________________
//
// SensorPixel is a data class for storing mapped data for
// one QDC detector
//

//
// $Id: SensorPixel.cpp,v 1.1 2012/02/09 23:51:09 dc Exp $
//
#include "USensorPixel.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif
using std::endl;
using std::cout;

//____________________________________________________________________
USensorPixel::USensorPixel() {
  fLine = 0;
  fColumn = 0;
  fPitch = 18.4 * 1e-3;  // in mm
  fUsed = false;
}

//____________________________________________________________________
USensorPixel::USensorPixel(const int iL, const int iC) {
  fLine = iL;
  fColumn = iC;
  fPitch = 18.4 * 1e-3;  // in mm
  fUsed = false;
}

//____________________________________________________________________
USensorPixel::~USensorPixel() {}
