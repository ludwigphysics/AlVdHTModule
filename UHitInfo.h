//-*-mode:c++-*-
//
#ifndef UTIL_UHitInfo
#define UTIL_UHitInfo

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TError
#include "TError.h"
#endif

class UHitInfo : public TObject {
 public:
  UHitInfo();
  virtual ~UHitInfo();

  void AddHitInfo(int hitId, int tabId) {
    fHitId[fEntries] = hitId;
    fTabId[fEntries] = tabId;
    fEntries++;
    if (fEntries > 199) {
      Error("AddHitInfo", "To many hits at this (ax,ay) location");
    }
  }

  void SetVds1Vds2Hit() { fVds1Vds2Hit = true; }

  bool ContainsVds1Vds2Hit() { return fVds1Vds2Hit; }
  int GetClusterStatus() { return fClusterStatus; }
  int GetIndex(int i) { return fHitId[i]; }
  int GetTabId(int i) { return fTabId[i]; }

  int GetEntries() { return fEntries; }
  int GetIndx() { return fIndx; }
  int GetIndy() { return fIndy; }

  void SetClusterStatus(const int i) { fClusterStatus = i; }
  void SetIJ(int i, int j) {
    fIndx = i;
    fIndy = j;
  }
  void Reset() { fEntries = 0; }

 private:
  int fEntries;       // number of hits
  int fHitId[200];    // hit location in table (assumed up to 10 hits in one HT space point)
  int fTabId[200];    // table id
  bool fVds1Vds2Hit;  // informs whether Vds1 or Vds2 has hit
  int fClusterStatus;
  int fIndx;  // line of hit info
  int fIndy;  // column of hit info

 public:
  friend ostream& operator<<(ostream& os, UHitInfo* stepinfo);
  //  ClassDef(UHitInfo,1)  //  DC raw data class
};

#endif

// $Log:$
