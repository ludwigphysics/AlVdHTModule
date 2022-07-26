// -*- mode: c++ -*-
//
// Module for primary track reconstruction using both VD arms
//
// Author: Paweł Staszel
//
//
#ifndef NA61_Na61AlPrimaryVertexRecoModule
#define NA61_Na61AlPrimaryVertexRecoModule

#ifndef NA61_Na61PrimaryVertexRecoModule
#include "Na61PrimaryVertexRecoModule.h"
#endif

#ifndef NA61_Na61VdParameters
#include "Na61VdParameters.h"
#endif

#ifndef NA61_Na61ArmParameters
#include "Na61ArmParameters.h"
#endif

#ifndef UTIL_UVdEvent
#include "UVdEvent.h"
#endif
#ifndef UTIL_USensorHit
#include "USensorHit.h"
#endif
#ifndef UTIL_UVdTrack
#include "UVdTrack.h"
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

class Na61AlPrimaryVertexRecoModule : public Na61PrimaryVertexRecoModule {
 public:
  Na61AlPrimaryVertexRecoModule();
  Na61AlPrimaryVertexRecoModule(const char* namec, const char* title);
  virtual ~Na61AlPrimaryVertexRecoModule() {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* outNodeJ, UEventNode* outNodeS);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option = "B") const;  // *MENU*


 private:

 public:

  Na61VdParameters* fParams;
  Na61ArmParameters* fArmParamsJ;
  Na61ArmParameters* fArmParamsS;

  int fRunId;

 private:

 public:
  //    ClassDef(Na61PrimaryVertexRecoModule,1) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
