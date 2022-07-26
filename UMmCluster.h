//-*-mode:c++-*-
//
#ifndef UTIL_UMmCluster
#define UTIL_UMmCluster

#ifndef UTIL_UDataObject
#include "UDataObject.h"
#endif

class UMmCluster : public UDataObject {
 public:
  UMmCluster();
  virtual ~UMmCluster();

  void UpdateCluster(int stripNum, double pos, double charge) {
    fStripNumber[fNstrips] = stripNum;
    fStripPos[fNstrips] = pos;
    fStripCharge[fNstrips] = charge;
    fNstrips++;
  }

  void SetHitPos(double v) { fHitPos = v; }
  void SetClosestStrip(int v) { fClosestStrip = v; }
  void SetMean(double v) { fMean = v; }
  void SetSigma(double v) { fSigma = v; }
  void SetOverlapParam(double v) { fOverlapParameter = v; }

  double GetHitPos() { return fHitPos; }
  int GetClosestStrip() { return fClosestStrip; }
  double GetSigma() { return fSigma; }
  double GetMean() { return fMean; }
  int GetNstrips() { return fNstrips; }
  double GetOverlapParam() { return fOverlapParameter; }

  int GetStripNumer(int is) { return fStripNumber[is]; }
  double GetStripPos(int is) { return fStripPos[is]; }
  double GetStripCharge(int is) { return fStripCharge[is]; }

 private:
  double fHitPos;
  int fClosestStrip;
  int fStripNumber[20];
  double fStripPos[20];
  double fStripCharge[20];
  int fNstrips;
  double fMean;
  double fSigma;
  double fOverlapParameter;

 public:
  friend ostream& operator<<(ostream& os, UMmCluster* stepinfo);
  //  ClassDef(UMmCluster,0)  //  DC raw data class
};

#endif

// $Log:$
