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
#ifndef NA61_Na61VdTrackingPostModule
#define NA61_Na61VdTrackingPostModule

//#ifndef NA61_Na61Module
//#include "Na61Module.h"
//#endif

#ifndef NA61_Na61VdTrackingModule
#include "Na61VdTrackingModule.h"
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

class Na61VdTrackingPostModule : public Na61VdTrackingModule
{ 
  
public:
  
  Na61VdTrackingPostModule();
  Na61VdTrackingPostModule(const Char_t* name, const Char_t* title);
  virtual ~Na61VdTrackingPostModule () {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNode, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option="B") const; // *MENU*



  void SetPrimaryVertex(Float_t vx, Float_t vy, Float_t vz){
    fPrimaryVertex -> SetX(vx); 
    fPrimaryVertex -> SetY(vy); 
    fPrimaryVertex -> SetZ(vz);
  }
  

private:


TH1F* fhFullTracksPost;

public:
  
  //ClassDef(Na61VdTrackingPostModule,0) // DOCUMENT ME
  
};

#endif
//____________________________________________________________________
//
