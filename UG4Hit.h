//-*-mode:c++-*-
//
#ifndef UTIL_UG4Hit
#define UTIL_UG4Hit

#ifndef UTIL_UDataObject
#include "UDataObject.h"
#endif

#include "TMath.h"

class UG4Hit : public UDataObject {
 public:
  UG4Hit();
  virtual ~UG4Hit();

  int GetTrackID() { return fTrackID; }
  int GetPdgID() { return fPdgId; }
  int GetParentPdgID() { return fParentPdgId; }
  int GetParentTrackID() { return fParentTrackId; }
  int GetEntryFlag() { return fAtEntry; }
  int GetExitFlag() { return fAtExit; }
  int GetPointId() { return fPointId; }
  int GetPoint() { return fPoint; }
  int GetFrontVtpc() { return fFrontVtpc; }
  // int  GetAddPoint()                { return fAddPoint; }

  double GetX() { return fPosIn[0]; }
  double GetY() { return fPosIn[1]; }
  double GetZ() { return fPosIn[2]; }
  double GetPX() { return fMomIn[0]; }
  double GetPY() { return fMomIn[1]; }
  double GetPZ() { return fMomIn[2]; }
  double GetEnergy() { return fEnergyIn; }

  double GetXOut() { return fPosOut[0]; }
  double GetYOut() { return fPosOut[1]; }
  double GetZOut() { return fPosOut[2]; }
  double GetPXOut() { return fMomOut[0]; }
  double GetPYOut() { return fMomOut[1]; }
  double GetPZOut() { return fMomOut[2]; }
  double GetEnergyOut() { return fEnergyOut; }
  double GetTotalDistance() { return fTotalDistance; }
  double GetDeviation() { return fDeviation; }
  int GetStationID() { return fStationID; }
  double GetP() { return TMath::Sqrt(GetPX() * GetPX() + GetPY() * GetPY() + GetPZ() * GetPZ()); }

  void SetEntryFlag(const int i) { fAtEntry = i; }
  void SetExitFlag(const int i) { fAtExit = i; }
  void SetTrackID(const int i) { fTrackID = i; }
  void SetPdgID(const int i) { fPdgId = i; }
  void SetParentPdgID(const int i) { fParentPdgId = i; }
  void SetParentTrackID(const int i) { fParentTrackId = i; }
  void SetPointId(const int i) { fPointId = i; }
  void SetPoint(const int i) { fPoint = i; }
  void SetFrontVtpc(const int i) { fFrontVtpc = i; }
  // void   SetAddPoint(const int i)  { fAddPoint  = i; }

  void SetPosition(const double x, const double y, const double z) {
    fPosIn[0] = x;
    fPosIn[1] = y;
    fPosIn[2] = z;
  }
  void SetMomentum(const double px, const double py, const double pz) {
    fMomIn[0] = px;
    fMomIn[1] = py;
    fMomIn[2] = pz;
  }
  void SetEnergy(const double e) { fEnergyIn = e; }

  void SetPositionOut(const double x, const double y, const double z) {
    fPosOut[0] = x;
    fPosOut[1] = y;
    fPosOut[2] = z;
  }
  void SetMomentumOut(const double px, const double py, const double pz) {
    fMomOut[0] = px;
    fMomOut[1] = py;
    fMomOut[2] = pz;
  }
  void SetEnergyOut(const double e) { fEnergyOut = e; }
  void SetTotalDistance(const double v) { fTotalDistance = v; }
  void SetDeviation(const double v) { fDeviation = v; }
  void SetStationID(int v) { fStationID = v; }

  void MarkAsFake() { fIsFake = true; }
  bool IsFake() { return fIsFake; }

  double FindDistance(UG4Hit* hit);
  double FindDistance();

 private:
  double fTotalDistance;
  double fDeviation;
  // double fX_out,  fY_out,  fZ_out;
  // double fPx_out, fPy_out, fPz_out;
  int fFrontVtpc;  //!
  int fAtEntry;
  int fAtExit;
  int fTrackID;
  int fPdgId;
  int fParentPdgId;
  int fParentTrackId;
  int fPointId;
  int fPoint;
  int fStationID;  //!
  // int fAddPoint;
  double fPosIn[3];
  double fMomIn[3];
  double fEnergyIn;
  double fPosOut[3];
  double fMomOut[3];
  double fEnergyOut;
  bool fIsFake;

 public:
  friend ostream& operator<<(ostream& os, UG4Hit* stepinfo);
  //  ClassDef(UG4Hit,7)  //  DC raw data class
};

#endif

// $Log:$
