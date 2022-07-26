//-*-mode:c++-*-
//
#ifndef UTIL_UG4RecoTrack
#define UTIL_UG4RecoTrack

#ifndef UTIL_UDataObject
#include "UDataObject.h"
#endif
#ifndef UTIL_UG4Hit
#include "UG4Hit.h"
#endif

#include "TLorentzVector.h"
#include "ULine3D.h"

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef UTIL_UVector3D
#include "UVector3D.h"
#endif
#ifndef WIN32
#include <iostream>
#else
#include <iostream.h>
#endif
#include <cassert>
#include "TMath.h"

class UG4RecoTrack : public UDataObject {
 public:
  UG4RecoTrack();
  UG4RecoTrack(const double px_, const double py_, const double pz_, const double e, Vector3D& origin);

  virtual ~UG4RecoTrack();

  int GetTrackID() { return fTrackID; }
  int GetPdgID() { return fPdgId; }
  int GetParentPdgID() { return fParentPdgId; }
  int GetParentTrackID() { return fParentTrackID; }
  double GetX() { return fline->GetOrigin().X(); }
  double GetY() { return fline->GetOrigin().Y(); }
  double GetZ() { return fline->GetOrigin().Z(); }

  double GetXatZ(double z) {
    double xatz = fline->GetOrigin().X() + (z - fline->GetOrigin().Z()) * fline->GetDirection().X() / fline->GetDirection().Z();
    return xatz;
  }

  double GetYatZ(double z) {
    double yatz = fline->GetOrigin().Y() + (z - fline->GetOrigin().Z()) * fline->GetDirection().Y() / fline->GetDirection().Z();
    return yatz;
  }

  double GetPX() { return fTLV->Px(); }
  double GetPY() { return fTLV->Py(); }
  double GetPZ() { return fTLV->Pz(); }
  double GetEnergy() { return fTLV->Energy(); }
  double GetMass() { return fTLV->M(); }
  double GetP() { return TMath::Sqrt(fTLV->Px() * fTLV->Px() + fTLV->Py() * fTLV->Py() + fTLV->Pz() * fTLV->Pz()); }
  TLorentzVector* GetTLV() { return fTLV; }
  Line3D* Getline() { return fline; }
  double GetTotalVtpcDistance() { return fTotalVtpcDistance; }
  int GetVtpcHitInd() { return fVtpcHitInd; }
  TString GetVtpcName() { return fVtpcName; }
  UG4Hit* GetFirstVdHit() { return fFirstVdHit; }
  double GetChi2Ndf() { return fChi2Ndf; }
  int GetWrongVdVtpcMatching() { return fWrongVdVtpcMatching; }
  double GetIsPion() { return fIsPion; }
  double GetIsKaon() { return fIsKaon; }
  double GetIsProt() { return fIsProt; }

  void SetTrackID(const int i) { fTrackID = i; }
  void SetPdgID(const int i) { fPdgId = i; }
  void SetParentPdgID(const int i) { fParentPdgId = i; }
  void SetParentTrackID(const int i) { fParentTrackID = i; }
  // void   SetPosition(const double x_,const double y_, const double z_)
  //{ x = x_;  y= y_; z = z_;}
  void SetEnergy(const double e) { fEnergy = e; }
  void SetMomentumAndEnergy(const double px_, const double py_, const double pz_, const double e) { fTLV->SetPxPyPzE(px_, py_, pz_, e); }

  void SetTotalVtpcDistance(double v) { fTotalVtpcDistance = v; }
  void SetVtpcName(const char* name) { fVtpcName = name; }
  void SetVtpcHitInd(int v) { fVtpcHitInd = v; }
  void SetFirstVdHit(UG4Hit* hit) { fFirstVdHit = hit; }
  void SetChi2Ndf(double v) { fChi2Ndf = v; }
  void SetWrongVdVtpcMatching(int v) { fWrongVdVtpcMatching = v; }

  void SetIsPion(double v) { fIsPion = v; }
  void SetIsKaon(double v) { fIsKaon = v; }
  void SetIsProt(double v) { fIsProt = v; }

  double FindDistCA(UG4RecoTrack* recotrack);

 private:
  // int fAtEntry;
  int fTrackID;
  int fPdgId;
  int fParentPdgId;
  int fParentTrackID;
  // double x,y,z;
  double fEnergy;
  double fMass;
  TLorentzVector* fTLV;
  Line3D* fline;
  // Line3D line;
  TString fVtpcName;
  int fVtpcHitInd;
  UG4Hit* fFirstVdHit;       //! not persistent
  double fChi2Ndf;           //! not persistent
  int fWrongVdVtpcMatching;  //! not persistent

  double fTotalVtpcDistance;

  double fIsPion;
  double fIsKaon;
  double fIsProt;

 public:
  friend ostream& operator<<(ostream& os, UG4RecoTrack* stepinfo);
  //  ClassDef(UG4RecoTrack,8)  //  DC raw data class
};

#endif

// $Log:$
