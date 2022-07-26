// -*- mode: c++ -*-
//
// This class reconstructs G4Na61 simulated hits using
// GEANT track information. This information is not available in the real data,
// that is why I put "bypass" in the name of this module.
// However, the reconstructed tracks can be later selected according to the
// detection system acceptance and smeared out accounting for the limited
// position/momentum resolution.
// In this way we can account for the general features of the measurement process
// avoiding complicated (and not yet developed) digitization and regular reconstruction
// procedures
//
// Author: Paweł Staszel
//
//
#ifndef NA61_Na61FrameMergingModule
#define NA61_Na61FrameMergingModule

#ifndef NA61_Na61Module
#include "Na61Module.h"
#endif
#ifndef UTIL_UEventNode
#include "UEventNode.h"
#endif
#ifndef UTIL_UEvent
#include "UEvent.h"
#endif
#ifndef UTIL_UG4Hit
#include "UG4Hit.h"
#endif

#include "TLorentzVector.h"
#include "ULine3D.h"

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
#include "TSystem.h"
#include "TDirectory.h"
#include "Math/QuantFuncMathCore.h"
#include <time.h>

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

class Na61FrameMergingModule : public Na61Module {
 public:
  Na61FrameMergingModule();
  Na61FrameMergingModule(const char* name, const char* title);
  virtual ~Na61FrameMergingModule() {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNode, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option = "B") const;  // *MtENU*

  // void SetSpreadMode(int v){fSpreadMode=v;
  // cout<<"fSpreadMode="<<fSpreadMode<<endl;
  //}

  void SetXRowOffset(int v) { fXRowOffset = v; }
  void SetRowOverlap(int v) { fRowOverlap = v; }
  void SetNewMerge(int v) { fNewMerge = v; }
  void SetMinTimer(double v) { fMinTimer = v; }
  void SetMaxTimer(double v) { fMaxTimer = v; }
  void SetFrameTimers(int i, double v0, double v1, double v2 = 999999, double v3 = 999999, double v4 = 999999) {
    // the timers should be given in micro-seconds
    // at lease v0 and v1 should be specified explicitely
    fTimers[i][0] = v0;
    fTimers[i][1] = v1;
    fTimers[i][2] = v2;
    fTimers[i][3] = v3;
    fTimers[i][4] = v4;
  }

  void SetRunId(int v) { fRunId = v; }
  void SetCutOnNoisyPixels(double v) { fPixelNoiseCut = v; }

  void SetToJuraArm(bool v) { fJura = v; }

  void SetHomePath(TString v) { fHomePath = v; }

 private:
  void MergeFrames(UDataTable** tabs_merge, UDataTable** tabs_f0, UDataTable** tabs_f1, UDataTable** tabs_f2, UDataTable** tabs_f3, UDataTable** tabs_f4);
  void MergeFrames_SIM(UDataTable** tabs_merge, UDataTable** tabs_f0);
  void MergeFrames_new(UDataTable** tabs_merge, UDataTable** tabs_f1, UDataTable** tabs_f2);

  void AddPixels(UDataTable** tabs_merge, UDataTable** tabs);
  void AddPixels_new(int frameId, int row_x, UDataTable** tabs_merge, UDataTable** tabs);
  void AddPixels_SIM(int frameId, UDataTable** tabs_merge, UDataTable** tabs);
  void ReadNoisyPixelInfo();
  void ResetFreqArray();

 private:
  int fCounter_f1_1;
  int fCounter_f1_2;
  int fCounter_f2_1;
  int fCounter_f2_2;
  int fRowOverlap;
  int fXRowOffset;

  bool fJura;

  int fNewMerge;
  double fM26CycleTime;
  bool fNoisyPixelsDefined;
  int fRunId;
  static const int findx = 576;
  static const int findy = 1152;
  TString fHomePath;

  double fPixelNoiseCut;
  int fEvents;

  double fMinTimer;
  double fMaxTimer;
  double fTimers[3][5];  // first index refers to peripheral fpga

  double fNVds;
  double fVdsZ[4];

  bool faxayFreq[8][576 + 1][1152 + 1];

  TH1F* fhPixels_f0[8];
  TH1F* fhPixels_f1[8];
  TH1F* fhPixels_f2[8];
  TH1F* fhPixels_f3[8];
  TH1F* fhPixelsMerged[8];

  TH2F* fhPixelMap_f34[8];
  TH1F* fhPixelFrequency[8];
  TH2F* fhNoisyPixelMap[8];

  TH1F* fhTimer_f0[3];
  TH1F* fhTimer_f1[3];
  TH1F* fhTimer_f2[3];
  TH1F* fhTimer_f3[3];
  TH1F* fhTimer_f4[3];

  TH1F* fhTimer_acce_f0[3];
  TH1F* fhTimer_acce_f1[3];
  TH1F* fhTimer_acce_f2[3];
  TH1F* fhTimer_acce_f3[3];
  TH1F* fhTimer_acce_f4[3];

  TH1F* fhPixelsVsTimer[8];
  TH1F* fhPixelsVsTimer_ref[8];

  TH2F* fhTimer_1vs2;
  TH2F* fhTimer_1vs3;
  TH2F* fhTimer_2vs3;

 public:
  //  ClassDef(Na61FrameMergingModule,0) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
