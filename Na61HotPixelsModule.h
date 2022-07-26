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
#ifndef NA61_Na61HotPixelsModule
#define NA61_Na61HotPixelsModule

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
#include "TRandom3.h"
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

class Na61HotPixelsModule : public Na61Module
{ 
  
public:
  
  Na61HotPixelsModule();
  Na61HotPixelsModule(const Char_t* name, const Char_t* title);
  virtual ~Na61HotPixelsModule () {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNode, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option="B") const; // *MtENU*

  //void SetSpreadMode(Int_t v){fSpreadMode=v;
  // cout<<"fSpreadMode="<<fSpreadMode<<endl;
  //} 

  void SetToJuraArm(bool v){fJura=v;}
  void SetRunId(Int_t v){fRunId=v;}
  void SetCutOnNoisyPixels(Float_t v){fPixelNoiseCut=v;}

  void SetHomePath(TString v){fHomePath = v;}

private:

  void AddPixels(UDataTable** tabs_out, UDataTable** tabs);
  void ReadNoisyPixelInfo();
  void ResetFreqArray();

private:

  bool fJura;
  Bool_t fNoisyPixelsDefined;
  Int_t fRunId;
  // MIMOSA
  //static const Int_t findx = 576;
  //static const Int_t findy = 1152;
  // ALPIDE
  static const Int_t findx = 512;
  static const Int_t findy = 1024;


  TString fHomePath;

  Float_t fPixelNoiseCut;
  Int_t fEvents;


  Float_t fNVds;
  Float_t fVdsZ[4];

  Bool_t faxayFreq[34][512+1][1024+1];

  TH1F* fhPixels[34];
  TH1F* fhCleanPixels[34];

  TH2F* fhPixelMap[34];
  TH1F* fhPixelFrequency[34];
  TH2F* fhNoisyPixelMap[34];


public:

  //ClassDef(Na61HotPixelsModule,0) // DOCUMENT ME
  
};

#endif
//____________________________________________________________________
//
