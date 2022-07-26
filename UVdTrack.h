//-*-mode:c++-*-
//
#ifndef UTIL_UVdTrack
#define UTIL_UVdTrack

#ifndef UTIL_UDataObject
#include "UDataObject.h"
#endif
#ifndef UTIL_UG4RecoTrack
#include "UG4RecoTrack.h"
#endif

#include "TLorentzVector.h"
#include "ULine3D.h"

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TF1
#include "TF1.h"
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
#include <utl/Vector.h>
#include <utl/Point.h>

class UVdTrack : public UDataObject {
 public:
  UVdTrack();
  UVdTrack(Vector3D& origin, Vector3D& direction);
  UVdTrack(const Vector3D& origin, const Vector3D& direction);
  UVdTrack(UG4RecoTrack* recotrack);
  UVdTrack(UVdTrack* vdtrack);

  virtual ~UVdTrack();

  void Activate();
  void ActivateTLV(double scaleFact = 1.0) {
    fTLV->SetPxPyPzE(fPx * scaleFact, fPy * scaleFact, fPz * scaleFact, 999);
    fMomentum = fTLV->P();
  }
  void ActivateTLV_Kf(double scaleFact = 1.0) {
    fTLV_kf->SetPxPyPzE(fPx_kf * scaleFact, fPy_kf * scaleFact, fPz_kf * scaleFact, 999);
    fKfMomentum = fTLV_kf->P();
  }

  unsigned long long GetTrackID() { return fTrackID; }
  int GetPdgID() { return fPdgId; }
  int GetParentPdgID() { return fParentPdgId; }
  int GetParentTrackID() { return fParentTrackID; }
  double GetX() { return fline->GetOrigin().X(); }
  double GetY() { return fline->GetOrigin().Y(); }
  double GetZ() { return fline->GetOrigin().Z(); }
  double GetDX() { return fline->GetDirection().X(); }
  double GetDY() { return fline->GetDirection().Y(); }
  double GetDZ() { return fline->GetDirection().Z(); }

  // front line
  double GetX_f() { return flinef->GetOrigin().X(); }
  double GetY_f() { return flinef->GetOrigin().Y(); }
  double GetZ_f() { return flinef->GetOrigin().Z(); }
  double GetDX_f() { return flinef->GetDirection().X(); }
  double GetDY_f() { return flinef->GetDirection().Y(); }
  double GetDZ_f() { return flinef->GetDirection().Z(); }

  // back line
  double GetX_b() { return flineb->GetOrigin().X(); }
  double GetY_b() { return flineb->GetOrigin().Y(); }
  double GetZ_b() { return flineb->GetOrigin().Z(); }
  double GetDX_b() { return flineb->GetDirection().X(); }
  double GetDY_b() { return flineb->GetDirection().Y(); }
  double GetDZ_b() { return flineb->GetDirection().Z(); }

  double GetXatZ(double z) {
    double xatz = fline->GetOrigin().X() + (z - fline->GetOrigin().Z()) * fline->GetDirection().X() / fline->GetDirection().Z();
    return xatz;
  }
  double GetYatZ(double z) {
    double yatz = fline->GetOrigin().Y() + (z - fline->GetOrigin().Z()) * fline->GetDirection().Y() / fline->GetDirection().Z();
    return yatz;
  }
  // front line
  double GetXatZ_f(double z) {
    double xatz = flinef->GetOrigin().X() + (z - flinef->GetOrigin().Z()) * flinef->GetDirection().X() / flinef->GetDirection().Z();
    return xatz;
  }
  double GetYatZ_f(double z) {
    double yatz = flinef->GetOrigin().Y() + (z - flinef->GetOrigin().Z()) * flinef->GetDirection().Y() / flinef->GetDirection().Z();
    return yatz;
  }

  // back line
  double GetXatZ_b(double z) {
    double xatz = flineb->GetOrigin().X() + (z - flineb->GetOrigin().Z()) * flineb->GetDirection().X() / flineb->GetDirection().Z();
    return xatz;
  }

  double GetYatZ_b(double z) {
    double yatz = flineb->GetOrigin().Y() + (z - flineb->GetOrigin().Z()) * flineb->GetDirection().Y() / flineb->GetDirection().Z();
    return yatz;
  }

  double GetXatZ_pol2(double z) {
    double xatz = ffPol2.Eval(z);
    return xatz;
  }

  Line3D* Getline() { return fline; }
  Line3D* Getlinef() { return flinef; }  // front line
  Line3D* Getlineb() { return flineb; }  // back line
  int GetVtpcHitInd() { return fVtpcHitInd; }
  TString GetVtpcName() { return fVtpcName; }

  // Vd tracks do not contain momentum information. The fTLV was only added for debugging purpose
  double GetPX() { return fTLV->Px(); }
  double GetPY() { return fTLV->Py(); }
  double GetPZ() { return fTLV->Pz(); }
  double GetEnergy() { return fTLV->Energy(); }
  double GetMass() { return fTLV->M(); }
  TLorentzVector* GetTLV() { return fTLV; }
  TLorentzVector* GetTLV_kf() { return fTLV_kf; }
  double GetChi2Ndf() { return fChi2Ndf; }
  TObject* GetHitIdAtStation(int i) { return fVdsHitIDs[i]; }
  int GetHitIndexOnStation(int i) { return fArrayIndex[i]; }
  int GetTabIndexOnStation(int i) { return fTabArrayIndex[i]; }
  int GetFlag() { return fTrackFlag; }
  int GetTpcMatchingFlag() { return fTpcMatchingFlag; }
  int GetNumberOfTpcClusters() { return fNumberOfTpcClusters; }
  int GetCharge() { return fCharge; }
  int GetTpcCharge() { return fTpcCharge; }
  double GetSlopeChange() { return fSlopeChange; }
  double GetCurvature() { return fCurvature; }
  double GetMomentum() { return fMomentum; }
  double GetKfMomentum() { return fKfMomentum; }
  int GetCombMeth() { return fCombMeth; }
  int GetTagForVtx() { return fTagForVtx; }
  int GetTpcTrackIndex() { return fTpcTrackIndex; }
  utl::Point GetPositionKf() { return fPositionKf; }
  utl::Vector GetMomentumKf() { return fMomentumKf; }
  int GetChargeKf() { return fChargeKf; }
  bool GetFitKf() { return fFitKf; }
  double GetDedx() { return fDedx; }

  Float_t GetPx(){return fPx;}
  Float_t GetPy(){return fPy;}
  Float_t GetPz(){return fPz;}
  Float_t GetPx_kf(){return fPx_kf;}
  Float_t GetPy_kf(){return fPy_kf;}
  Float_t GetPz_kf(){return fPz_kf;}

  void SetMomentum(double px, double py, double pz) {
    fPx = px;
    fPy = py;
    fPz = pz;
  }
  void SetKfMomentum(double px, double py, double pz) {
    fPx_kf = px;
    fPy_kf = py;
    fPz_kf = pz;
  }

  void SetPol2(TF1 pol2) { ffPol2 = pol2; }
  void SetFlag(int i) { fTrackFlag = i; }
  void SetVtpcName(const char* name) { fVtpcName = name; }
  void SetVtpcHitInd(int i) { fVtpcHitInd = i; }
  void SetTrackID(unsigned long long i) { fTrackID = i; }
  void SetPdgID(const int i) { fPdgId = i; }
  void SetParentPdgID(const int i) { fParentPdgId = i; }
  void SetParentTrackID(const int i) { fParentTrackID = i; }
  void SetChi2Ndf(double v) { fChi2Ndf = v; }
  void SetTpcMatchingFlag(int v) { fTpcMatchingFlag = v; }
  void SetCharge(int v) { fCharge = v; }
  void SetTpcCharge(int v) { fTpcCharge = v; }
  void SetSlopeChange(double v) { fSlopeChange = v; }
  void SetCurvature(double v) { fCurvature = v; }
  void SetMomentum(double v) { fMomentum = v; }
  void SetKfMomentum(double v) { fKfMomentum = v; }
  void SetCombMeth(int v) { fCombMeth = v; }
  void SetTagForVtx(int v) { fTagForVtx = v; }
  void SetTpcTrackIndex(int v) { fTpcTrackIndex = v; }
  void SetNumberOfTpcClusters(int v) { fNumberOfTpcClusters = v; }
  void SetPositionKf(const utl::Point& v) { fPositionKf = v; }
  void SetMomentumKf(const utl::Vector& v) { fMomentumKf = v; }
  void SetChargeKf(int v) { fChargeKf = v; }
  void SetFitKf(bool v) { fFitKf = v; }
  void SetDedx(double v) { fDedx = v; }

  void SetVdsHitIDs(TObject* hit1, TObject* hit2, TObject* hit3, TObject* hit4) {
    // hit id are hit object pointers converted to long. They are unique.
    fVdsHitIDs[0] = (TObject*)hit1;
    fVdsHitIDs[1] = (TObject*)hit2;
    fVdsHitIDs[2] = (TObject*)hit3;
    fVdsHitIDs[3] = (TObject*)hit4;
  }

  void SetHitArrayIndexes(int ai1, int ai2, int ai3, int ai4) {
    // hit id are hit object pointers converted to long. They are unique.
    fArrayIndex[0] = ai1;
    fArrayIndex[1] = ai2;
    fArrayIndex[2] = ai3;
    fArrayIndex[3] = ai4;
  }

  void SetTabArrayIndexes(int ai1, int ai2, int ai3, int ai4) {
    // hit id are hit object pointers converted to long. They are unique.
    fTabArrayIndex[0] = ai1;
    fTabArrayIndex[1] = ai2;
    fTabArrayIndex[2] = ai3;
    fTabArrayIndex[3] = ai4;
  }

  double GetDirectionX() { return fDirectionX; }

  double FindDistCA(UVdTrack* vdtrack);
  void MarkForRemoval(bool v) { fRemoveTrack = v; }
  bool IsMarkedForRemoval() { return fRemoveTrack; }

 void SetLineParams(){
    fOriginX = fline->GetOrigin().X();
    fOriginY = fline->GetOrigin().Y();
    fOriginZ = fline->GetOrigin().Z();
    fDirectionX = fline->GetDirection().X();
    fDirectionY = fline->GetDirection().Y();
    fDirectionZ = fline->GetDirection().Z();
  }
  void SetLinefParams() {
    fOriginX_f = flinef->GetOrigin().X();
    fOriginY_f = flinef->GetOrigin().Y();
    fOriginZ_f = flinef->GetOrigin().Z();
    fDirectionX_f = flinef->GetDirection().X();
    fDirectionY_f = flinef->GetDirection().Y();
    fDirectionZ_f = flinef->GetDirection().Z();
  }
  void SetLinebParams() {
    fOriginX_b = flineb->GetOrigin().X();
    fOriginY_b = flineb->GetOrigin().Y();
    fOriginZ_b = flineb->GetOrigin().Z();
    fDirectionX_b = flineb->GetDirection().X();
    fDirectionY_b = flineb->GetDirection().Y();
    fDirectionZ_b = flineb->GetDirection().Z();
  }

 private:
  // int fAtEntry;
  unsigned long long fTrackID;
  int fPdgId;
  int fParentPdgId;
  int fParentTrackID;
  double fPx;
  double fPy;
  double fPz;
  double fPx_kf;
  double fPy_kf;
  double fPz_kf;
  double fMomentum;         // ! same information is contaided in fPx, fPy, fPz
  double fKfMomentum;       // ! momentum from Kalman Fitter
  TLorentzVector* fTLV;     // ! to allow simplify the kinamatic algebra
  TLorentzVector* fTLV_kf;  // ! same as fTLV but contains from Kalman Fitter
  Line3D* fline;            // !
  Line3D* flinef;           // ! front track line for field runs
  Line3D* flineb;           // ! back track line for field runs
  TF1 ffPol2;               // ! pol2 line fitted in back mode

  double fDirectionX;
  double fDirectionY;
  double fDirectionZ;
  double fOriginX;
  double fOriginY;
  double fOriginZ;

  double fDirectionX_f;
  double fDirectionY_f;
  double fDirectionZ_f;
  double fOriginX_f;
  double fOriginY_f;
  double fOriginZ_f;

  double fDirectionX_b;
  double fDirectionY_b;
  double fDirectionZ_b;
  double fOriginX_b;
  double fOriginY_b;
  double fOriginZ_b;
  //
  int fCombMeth;
  int fTagForVtx;

  int fVtpcHitInd;
  TString fVtpcName;
  double fChi2Ndf;
  TObject* fVdsHitIDs[4];  // !
  bool fRemoveTrack;        // !
  int fTrackFlag;           // !
  int fArrayIndex[4];       //   to keep information about hit position in the table.
  int fTabArrayIndex[4];    //   to keep information about hit position in the table.
  int fTpcMatchingFlag;
  int fNumberOfTpcClusters;
  int fCharge;     //   track charge in units of elementary charge
  int fTpcCharge;  //   track charge in units of elementary charge
  double fSlopeChange;
  double fCurvature;
  int fTpcTrackIndex;  // !
  utl::Point fPositionKf;
  utl::Vector fMomentumKf;
  double fChargeKf;
  bool fFitKf;
  double fDedx;

 public:
  friend ostream& operator<<(ostream& os, UVdTrack* stepinfo);
  //  ClassDef(UVdTrack,21)  //  DC raw data class
};

#endif

// $Log:$
