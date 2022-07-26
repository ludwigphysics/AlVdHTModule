//-*-mode:c++-*-

#ifndef U_USensorHit
#define U_USensorHit

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif

class USensorHit : public TObject {
 public:
  USensorHit();
  USensorHit(const double x, const double y, const double z, const int clustsize);
  virtual ~USensorHit();

  void SetLocalX(double v) { fLocalX = v; }
  void SetLocalY(double v) { fLocalY = v; }

  void SetX(double v) { fX = v; }
  void SetY(double v) { fY = v; }
  void SetZ(double v) { fZ = v; }

  void SetDevX(double v) { fDevX = v; }
  void SetDevY(double v) { fDevY = v; }
  void SetDevZ(double v) { fDevZ = v; }

  void SetClusterSize(int v) { fClusterSize = v; }

  void SetClusterLine(int v) { fClusterLine = v; }
  void SetClusterColumn(int v) { fClusterColumn = v; }

  void SetStationID(int v) { fStationID = v; }
  void SetMatchedToTPCTrack(int v) { fMatchedToTPCTrack = v; }
  void SetSensorName(const char* name) { fSensorName = name; }
  void SetIndex(int v) { fIndex = v; }
  void SetTrackIndex(int v) { fTrackIndex = v; }

  void SetBelongsToPrimaryTrack(bool bv) { fBelongsToPrimaryTrack = bv; }
  void SetUsed(bool bv) { fUsed = bv; }

  double GetLocalX() { return fLocalX; }
  double GetLocalY() { return fLocalY; }

  double GetX() { return fX; }
  double GetY() { return fY; }
  double GetZ() { return fZ; }

  double GetDevX() { return fDevX; }
  double GetDevY() { return fDevY; }
  double GetDevZ() { return fDevZ; }

  int GetClusterSize() { return fClusterSize; }
  int GetClusterLine() { return fClusterLine; }
  int GetClusterColumn() { return fClusterColumn; }

  int GetMatchedToTPCTrack() { return fMatchedToTPCTrack; }
  int GetStationID() { return fStationID; }
  bool GetBelongsToPrimaryTrack() { return fBelongsToPrimaryTrack; }
  bool GetUsed() { return fUsed; }

  TString GetSensorName() { return fSensorName; }
  int GetIndex() { return fIndex; }
  int GetTrackIndex() { return fTrackIndex; }

 private:
  double fLocalX;  // to be used for geometry tunning with minuit
  double fLocalY;  //

  double fX;  //
  double fY;  //
  double fZ;  //

  double fDevX;  // !
  double fDevY;  // !
  double fDevZ;  // !

  int fClusterSize;
  int fClusterLine;
  int fClusterColumn;

  int fTrackIndex;
  int fMatchedToTPCTrack;
  int fStationID;  //!
  int fIndex;      //!

  TString fSensorName;

  bool fBelongsToPrimaryTrack;  //!
  bool fUsed;                   //!

 public:
  friend ostream& operator<<(ostream& os, USensorHit* pixel);

  //  ClassDef(USensorHit,8)  //  SensorHit data class
};

#endif
