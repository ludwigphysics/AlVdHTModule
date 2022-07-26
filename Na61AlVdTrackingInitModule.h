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
#ifndef NA61_Na61AlVdTrackingInitModule
#define NA61_Na61AlVdTrackingInitModule

#ifndef NA61_Na61Module
#include "Na61Module.h"
#endif
#ifndef NA61_Na61AlVdTrackingModule
#include "Na61AlVdTrackingModule.h"
#endif
#ifndef NA61_Na61VdParameters
#include "Na61VdParameters.h"
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


class Na61AlVdTrackingInitModule : public Na61AlVdTrackingModule
{ 
  
public:
  
  Na61AlVdTrackingInitModule();
  Na61AlVdTrackingInitModule(const Char_t* name, const Char_t* title);
  virtual ~Na61AlVdTrackingInitModule () {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNode, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option="B") const; // *MENU*


private:


  //void Make3HitTracks(TString devname, UDataTable** hits1, UDataTable** hits2, UDataTable** hits3,UEventNode* out, Int_t* TabIndArray);

  Int_t FindDevId(TString devname){ 
    for(Int_t i=0;i<100;i++) if(fDevName[i].EqualTo(devname))return i;
    
    cout<<"devname not found on the list - check it out"<<endl;
    return -1111; 
  }

  
private:

  Float_t fNsig_dev;  
  Int_t fEvents;

  
public:
  
  //  ClassDef(Na61AlVdTrackingInitModule,0) // DOCUMENT ME
  
};

#endif
//____________________________________________________________________
//
