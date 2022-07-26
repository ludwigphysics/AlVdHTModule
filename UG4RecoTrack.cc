
#include "UG4RecoTrack.h"
#include "TLorentzVector.h"
#include "ULine3D.h"

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef UTIL_UVector3D
#include "UVector3D.h"
#endif

#ifndef __IOSTREAM__
#include <iostream>
#endif

using std::endl;

// ClassImp(UG4RecoTrack)

UG4RecoTrack::UG4RecoTrack() {
  fline = new Line3D();
  // fline->GetOrigin();
  // fline->Print();
  fTLV = new TLorentzVector();
  fTotalVtpcDistance = 0;
}

UG4RecoTrack::UG4RecoTrack(const double px_, const double py_, const double pz_, const double e, Vector3D& /*origin*/) {
  fTLV = new TLorentzVector(px_, py_, pz_, e);
  fTotalVtpcDistance = 0;
}

UG4RecoTrack::~UG4RecoTrack() {
  if (fline != 0) delete fline;
  if (fTLV != 0) delete fTLV;
}

//__________________________________________________________________________

double UG4RecoTrack::FindDistCA(UG4RecoTrack* recotrack) {
  Line3D line(recotrack->Getline()->GetOrigin(), recotrack->Getline()->GetDirection());
  double Z = fline->GetClosestProximityPoint(line).Z();
  return Z;
}

//_______________________________________________________________________

ostream& operator<<(ostream& os, UG4RecoTrack* track) {
  os << " TrackID==" << track->GetTrackID() << endl;

  return os;
}
