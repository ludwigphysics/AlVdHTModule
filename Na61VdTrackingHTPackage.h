// -*- mode: c++ -*-
//
// Package build to all for simple traking for secondary tracks
// using HT method,
//
// Author: Paweł Staszel
//
//
#ifndef NA61_Na61VdTrackingHTPackage
#define NA61_Na61VdTrackingHTPackage

//#ifndef NA61_Na61Module
//#include "Na61Module.h"
//#endif

#ifndef NA61_Na61VdTrackingHTModule
#include "Na61VdTrackingHTModule.h"
#endif

#ifndef NA61_Na61VdParameters
#include "Na61VdParameters.h"
#endif

#ifndef UTIL_UEventNode
#include "UEventNode.h"
#endif
#ifndef UTIL_UVdEvent
#include "UVdEvent.h"
#endif
#ifndef UTIL_USensorHit
#include "USensorHit.h"
#endif
#ifndef UTIL_UHitInfo
#include "UHitInfo.h"
#endif
#ifndef UTIL_UVdTrack
#include "UVdTrack.h"
#endif
#ifndef ROOT_TH1F
#include "TH1F.h"
#endif
#ifndef ROOT_TH2F
#include "TH2F.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif
#ifndef ROOT_TChain
#include "TChain.h"
#endif
#ifndef ROOT_TFile
#include "TFile.h"
#endif
#ifndef ROOT_TObjArray
#include "TObjArray.h"
#endif

#ifndef AccSteppingAction_h
#define AccSteppingAction_h
#endif
#ifndef ROOT_TGraph
#include "TGraph.h"
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include "cmath"

using std::cout;
using std::endl;
using std::ios;
using std::setprecision;
using std::setw;
using std::ifstream;
using std::cin;

class Na61VdTrackingHTPackage : public Na61VdTrackingModule {
 public:
  Na61VdTrackingHTPackage();
  Na61VdTrackingHTPackage(const char* name, const char* title);
  virtual ~Na61VdTrackingHTPackage() {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option = "B") const;  // *MtENU*

  void SetDistCut(double v) { fDistCut = v; }
  void SetVtpcSetup(int v) { fOldSetup = v; }
  void SetFillHitOccupancy(int v) { fFillHitOccupancy = v; }
  void SetAddInGrid2Tracks(int v) { fAddInGrid2Tracks = v; }
  void SetFiducialCut(double v) { fFiduCut = v; }
  void SetPosRes(double v) { fPosRes = v; }
  void Setbw(double vx, double vy) {
    fbwx = vx;
    fbwy = vy;
  }
  void SetProduction(int v) { fProduction = v; }
  void SetMakeClusters(int v) { fMakeClusters = v; }
  void SetChi2Cut(double v) { fChi2Cut = v; }
  void SetVzOffset(double v) { fVzOffset = v; }
  void SetStep(double v) { fStep = v; }
  void SetRangeToCover(double v) { fRangeToCover = v; }

  int GetAllTracks(int i) { return fAllTracks[i]; }

 private:
  void CleanTracks(UEventNode* out);
  void RemoveMultiplyTracks(TString devname1, TString devname2, UEventNode* out);
  bool TracksMatch(UVdTrack* track1, UVdTrack* track2);

 private:
  Na61VdTrackingHTModule* fTrackingHTModules[20];

  double fStep;          // step between ranges (typicaly 5 mm)
  double fRangeToCover;  // range of secondary vertexes to cover
  int fNModules;

  int fAllTracks[20];

  int fNRes;

  int fProduction;

  int findx;
  int findy;
  int fRingLevel;

  double fChi2Cut;

  double fPosRes;
  double fbwx;
  double fbwy;
  double fFiduCut;
  double fDistCut;
  double fVzOffset;
  int fOldSetup;
  int fFillHitOccupancy;
  int fAddInGrid2Tracks;
  int fMakeClusters;

 public:
  //  ClassDef(Na61VdTrackingHTPackage,0) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
