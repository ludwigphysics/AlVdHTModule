//____________________________________________________________________
//
// UG4Hit is a data class for storing geant4 hit information
//
//
#include "UG4Hit.h"
#ifndef __IOSTREAM__
#include <iostream>
#endif
using std::endl;

//____________________________________________________________________
// ClassImp(UG4Hit);

//____________________________________________________________________
UG4Hit::UG4Hit() {
  fTotalDistance = 0;
  fDeviation = 11111;
  fIsFake = false;
  fFrontVtpc = 0;
  // we assume the there is no default value allowed for proper operation
}
//____________________________________________________________________
UG4Hit::~UG4Hit() {}

//____________________________________________________________________
double UG4Hit::FindDistance(UG4Hit* hit) { return TMath::Sqrt((hit->GetX() - fPosIn[0]) * (hit->GetX() - fPosIn[0]) + (hit->GetY() - fPosIn[1]) * (hit->GetY() - fPosIn[1]) + (hit->GetZ() - fPosIn[2]) * (hit->GetZ() - fPosIn[2])); }

//____________________________________________________________________
double UG4Hit::FindDistance() { return TMath::Sqrt((fPosOut[0] - fPosIn[0]) * (fPosOut[0] - fPosIn[0]) + (fPosOut[1] - fPosIn[1]) * (fPosOut[1] - fPosIn[1]) + (fPosOut[2] - fPosIn[2]) * (fPosOut[2] - fPosIn[2])); }

//____________________________________________________________________
ostream& operator<<(ostream& os, UG4Hit* hit) {
  os << " TrackID=" << hit->GetTrackID() << endl;
  return os;
}
