//____________________________________________________________________
//
// SensorPixel is a data class for storing mapped data for
// one QDC detector
//

//
// $Id: SensorPixel.cpp,v 1.1 2012/02/09 23:51:09 dc Exp $
//
#include "UChipPixel.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif
using std::endl;
using std::cout;

//____________________________________________________________________
UChipPixel::UChipPixel() {
  fLine = 0;
  fColumn = 0;
  fPitchY       = 29.24*1e-3; // in mm
  fPitchX       = 26.88*1e-3; // in mm
  fOffX = 0.6;
  fUsed = false;
}

//____________________________________________________________________
UChipPixel::UChipPixel(const int iL, const int iC) {
  fLine = iL;
  fColumn = iC;
  fPitchY       = 29.24*1e-3; // in mm
  fPitchX       = 26.88*1e-3; // in mm
  fOffX = 0.6;
  fUsed = false;
}

//____________________________________________________________________
UChipPixel::~UChipPixel() {}
