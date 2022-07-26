#include "UVdTrack.h"
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

// ClassImp(UVdTrack)

UVdTrack::UVdTrack() {
  fline = new Line3D();
  flinef = new Line3D();
  flineb = new Line3D();
  fTLV = new TLorentzVector();
  fTLV_kf = new TLorentzVector();
  fArrayIndex[0] = -1;
  fArrayIndex[1] = -1;
  fArrayIndex[2] = -1;
  fArrayIndex[3] = -1;
  fTabArrayIndex[0] = -1;
  fTabArrayIndex[1] = -1;
  fTabArrayIndex[2] = -1;
  fTabArrayIndex[3] = -1;
  fTpcMatchingFlag = 0;
  fCharge = 11111;
  fCombMeth = 0;
  fTagForVtx = 0;
  fRemoveTrack = false;
  fFitKf = false;
}

UVdTrack::UVdTrack(Vector3D& origin, Vector3D& direction) {
  fline = new Line3D(origin, direction);
  flinef = new Line3D(origin, direction);
  flineb = new Line3D(origin, direction);
  fTLV = new TLorentzVector();
  fTLV_kf = new TLorentzVector();
  fRemoveTrack = false;
  fDirectionX = direction.X();
  fDirectionY = direction.Y();
  fDirectionZ = direction.Z();
  fOriginX = origin.X();
  fOriginY = origin.Y();
  fOriginZ = origin.Z();
  fArrayIndex[0] = -1;
  fArrayIndex[1] = -1;
  fArrayIndex[2] = -1;
  fArrayIndex[3] = -1;
  fTabArrayIndex[0] = -1;
  fTabArrayIndex[1] = -1;
  fTabArrayIndex[2] = -1;
  fTabArrayIndex[3] = -1;
  fTpcMatchingFlag = 0;
  fCharge = 11111;
  fCombMeth = 0;
  fTagForVtx = 0;
  fFitKf = false;
}

UVdTrack::UVdTrack(const Vector3D& origin, const Vector3D& direction) {
  fline = new Line3D(origin, direction);
  flinef = new Line3D(origin, direction);
  flineb = new Line3D(origin, direction);
  fTLV = new TLorentzVector();
  fTLV_kf = new TLorentzVector();
  fRemoveTrack = false;
  fDirectionX = direction.X();
  fDirectionY = direction.Y();
  fDirectionZ = direction.Z();
  fOriginX = origin.X();
  fOriginY = origin.Y();
  fOriginZ = origin.Z();
  fArrayIndex[0] = -1;
  fArrayIndex[1] = -1;
  fArrayIndex[2] = -1;
  fArrayIndex[3] = -1;
  fTpcMatchingFlag = 0;
  fCharge = 11111;
  fCombMeth = 0;
  fTagForVtx = 0;
  fFitKf = false;
}

//__________________________________________________________________________
void UVdTrack::Activate() {
  fline->SetOrigin(fOriginX, fOriginY, fOriginZ);
  fline->SetDirection(fDirectionX, fDirectionY, fDirectionZ);

  flinef->SetOrigin(fOriginX_f, fOriginY_f, fOriginZ_f);
  flinef->SetDirection(fDirectionX_f, fDirectionY_f, fDirectionZ_f);

  flineb->SetOrigin(fOriginX_b, fOriginY_b, fOriginZ_b);
  flineb->SetDirection(fDirectionX_b, fDirectionY_b, fDirectionZ_b);
}

//__________________________________________________________________________
UVdTrack::UVdTrack(UG4RecoTrack* recotrack) {
  fline = new Line3D(recotrack->Getline()->GetOrigin(), recotrack->Getline()->GetDirection());
  flinef = new Line3D(recotrack->Getline()->GetOrigin(), recotrack->Getline()->GetDirection());
  flineb = new Line3D(recotrack->Getline()->GetOrigin(), recotrack->Getline()->GetDirection());
  this->SetTrackID(recotrack->GetTrackID());
  this->SetPdgID(recotrack->GetPdgID());
  this->SetParentPdgID(recotrack->GetParentPdgID());
  this->SetParentTrackID(recotrack->GetParentTrackID());
  this->SetVtpcHitInd(recotrack->GetVtpcHitInd());
  this->SetVtpcName(recotrack->GetVtpcName());
  fTLV = new TLorentzVector();
  fTLV->SetPxPyPzE(recotrack->GetPX(), recotrack->GetPY(), recotrack->GetPZ(), recotrack->GetEnergy());
  /// maybe setting of fTLV_kf will be needed at some point
  fRemoveTrack = false;
  fArrayIndex[0] = -1;
  fArrayIndex[1] = -1;
  fArrayIndex[2] = -1;
  fArrayIndex[3] = -1;
  fTpcMatchingFlag = 0;
  fCharge = 11111;
  fCombMeth = 0;
  fTagForVtx = 0;
}

//__________________________________________________________________________

UVdTrack::UVdTrack(UVdTrack* track)
{
  fline = new Line3D((track->Getline())->GetOrigin(),(track->Getline())->GetDirection());
  flinef = new Line3D((track->Getlinef())->GetOrigin(),(track->Getlinef())->GetDirection());
  flineb = new Line3D((track->Getlineb())->GetOrigin(),(track->Getlineb())->GetDirection());
  fTLV =  new TLorentzVector();
  fTLV_kf =  new TLorentzVector();
  fRemoveTrack =  track->IsMarkedForRemoval();

  SetLineParams();
  SetLinefParams();
  SetLinebParams();

  SetVdsHitIDs(track->GetHitIdAtStation(0),
               track->GetHitIdAtStation(1),
               track->GetHitIdAtStation(2),
               track->GetHitIdAtStation(3));

  SetHitArrayIndexes(track->GetHitIndexOnStation(0),
                     track->GetHitIndexOnStation(1),
                     track->GetHitIndexOnStation(2),
                     track->GetHitIndexOnStation(3));

  SetTabArrayIndexes(track->GetTabIndexOnStation(0),
                     track->GetTabIndexOnStation(1),
                     track->GetTabIndexOnStation(2),
                     track->GetTabIndexOnStation(3));

  SetMomentum(track->GetPx(),track->GetPy(),track->GetPz());
  SetKfMomentum(track->GetPx_kf(),track->GetPy_kf(),track->GetPz_kf());

  fTpcMatchingFlag = 0;
  fCharge = track->GetCharge();
  fTpcCharge = track->GetTpcCharge();
  fCombMeth = track->GetCombMeth();
  fTagForVtx = track->GetTagForVtx();

  SetFlag(track->GetFlag());
  SetVtpcName(track->GetVtpcName());
  SetVtpcHitInd(track->GetVtpcHitInd());
  SetTrackID(track->GetTrackID()*100); // it is easy to connect track IDs
  SetChi2Ndf(track->GetChi2Ndf());
  SetSlopeChange(track->GetSlopeChange());
  SetCurvature(track->GetCurvature());

}

//__________________________________________________________________________
UVdTrack::~UVdTrack() {
  if (fline != 0) delete fline;
  if (flinef != 0) delete flinef;
  if (flineb != 0) delete flineb;
  if (fTLV != 0) delete fTLV;
  if (fTLV_kf != 0) delete fTLV_kf;
}

//__________________________________________________________________________

double UVdTrack::FindDistCA(UVdTrack* recotrack) {
  Line3D line(recotrack->Getline()->GetOrigin(), recotrack->Getline()->GetDirection());
  double Z = fline->GetClosestProximityPoint(line).Z();
  return Z;
}

//__________________________________________________________________________

ostream& operator<<(ostream& os, UVdTrack* track) {
  os << " TrackID==" << track->GetTrackID();

  return os;
}
