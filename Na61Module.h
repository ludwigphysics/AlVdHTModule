// -*- mode: c++ -*-
// base module class use in vd-standalone reconstruction
#ifndef NA61_Na61Module
#define NA61_Na61Module

#include "UEvent.h"

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

#ifndef ROOT_TFile
#include "TFile.h"
#endif

// Terror is needed due to kWarning
#ifndef ROOT_TError
#include "TError.h"
#endif

#ifndef ROOT_TStopwatch
#include "TStopwatch.h"
#endif

#ifndef __IOSTREAM__
#include <iostream>
#endif
#ifndef __IOMANIP__
#include <iomanip>
#endif

//#define UEvent UEvent_ // hack due to clash with UniGen

using std::cout;
using std::endl;
using std::flush;
using std::setw;
using std::setprecision;

class TList;
class UEventNode;
//class UEvent;

// static TFile* gHistFile;

class Na61Module : public TNamed {
 public:
  // State and Status types
  enum ENa61ModuleState {
    kSetup = 1,  // At creation/setup time
    kInit,       // At Init   time
    kBegin,      // At Begin  time
    kEvent,      // At Event  time
    kEnd,        // At End    time
    kFinish      // At Finish time
  };
  enum ENa61ModuleStatus {
    kOk = ::kWarning - 100,     // Module is doing fine
    kStop = ::kWarning + 100,   // Module wants pipeline stoped
    kFailure = ::kError + 100,  // Module failed
    kAbort = ::kError + 200     // Module failed miserably
  };
  enum {
    kShowFailures = 3,  // Show all failures
    kShowStop = 10,     // Show all stop conditions
  };

  TString fSensorNames[8];
  TString fAlSensorNames[34];
  TString fDevNames[20];
  TString fMatchStr[18];

  double fSensorXwidth;
  double fSensorYwidth;

  double fElecM;
  double fPionM;
  double fKaonM;
  double fProtM;

  double fK0sM;
  double fLambM;
  double fDPlusM;
  double fD0M;


  Int_t GetAlSensorNumberByName(TString chip_str){
    for(Int_t i=0;i<34;i++) if(fAlSensorNames[i].EqualTo(chip_str))return i;
    
    cout<<"chip not found on the list - check it out"<<endl;
    return -1111; 
  }

  int GetIdev(TString devname);
  int GetMatchIdev(TString matchname);
  double GetSensorXwidth() { return fSensorXwidth; }
  double GetSensorYwidth() { return fSensorYwidth; }

 private:
  // Level for amount of debug information
  // output during execution of the module entries
  // The application programmer should use the
  // member functions DebugLevel() for decision to
  // output information
  int fDebugLevel;       //! Debug level
  mutable int fVerbose;  //! Option flag

  // Statistics
  TStopwatch fTimer;          //! Timer for stats
  double fCpuTime;            //! Culminative CPU time spend
  int fEventCalled;           //! Counter times Event is called.
  mutable int fStopCount;     //! Counter for Stop Output
  mutable int fFailureCount;  //! Counter for Failure Output
  int fStopLimit;             //! Limit for Stop Count
  int fFailureLimit;          //! Limit for Failure Count

  // Histograms and required tables
  TList* fHistograms;         //! Histogram list
  TList* fRequiredTableList;  //! List of needed tables (not used)

 protected:
  // State and status of module
  mutable ENa61ModuleState fState;  //! What state is the module in
  mutable int fStatus;              //! How fares the module

  // Histogram flags
  bool fHistOn;      //! Histogramming enabled
  bool fHistBooked;  //! Histogrammed Booked

  int fDetectorSetupOption;

  virtual void SetState(ENa61ModuleState state) { fState = state; }
  virtual void SetStatus(int status) { fStatus = status; }

 public:
  Na61Module();
  Na61Module(const char* name, const char* title);
  virtual ~Na61Module();

  virtual void Init();     // Called once per job
  virtual void Begin() { SetState(kBegin); }  // Called once per run
  virtual void Event();                       // Called once per event
  virtual void Event(UEventNode* inev, UEventNode* outev);
  virtual void End() { SetState(kEnd); }        // Called once per run
  virtual void Finish() { SetState(kFinish); }  // Called once per job

  virtual void SetInputEvent(UEvent* /*Ev*/){};   // abstract class, should be defined in IO module
  virtual void SetOutputEvent(UEvent* /*Ev*/){};  // abstract class, should be defined in IO module
  // Histograms
  virtual void Book();
  void DisableHistograms();
  void EnableHistograms();
  bool HistOn() const { return fHistOn; }
  bool HistBooked() const { return fHistBooked; }
  void SetHistOn() { fHistOn = true; }
  void SetHistOff() { fHistOn = false; }
  void SetHistBooked() { fHistBooked = true; }
  virtual void DefineHistograms(){};
  TList* HistogramList() const { return fHistograms; }

  // Required tables (not used)
  virtual void AddRequiredTable(TObject* tablename);
  virtual TList* GetRequiredTableList() { return fRequiredTableList; }
  void SetRequiredData(char* name);

  // Statistics
  virtual void EventStatisticsStart();  // To be called by asp
  virtual void EventStatisticsEnd();    // To be called by asp
  virtual void ListEventStatistics();   // To be called by asp

  virtual void SetDetectorSetupOption(int v) { fDetectorSetupOption = v; }

  // Error, State and Status handling
  void Stop(const char* location, const char* msgfmt, ...) const;
  void Failure(const char* location, const char* msgfmt, ...) const;
  void Abort(const char* location, const char* msgfmt, ...) const;
  void Info(int lvl, const char* location, const char* msgfmt, ...) const;
  void Debug(int lvl, const char* location, const char* msgfmt, ...) const;
  virtual void SetStopLimit(int limit = 10) { fStopLimit = limit; }
  virtual void SetFailureLimit(int limit = 10) { fFailureLimit = limit; }
  virtual int GetStopLimit() const { return fStopLimit; }
  virtual int GetFailureLimit() const { return fFailureLimit; }
  virtual int GetStopCount() const { return fStopCount; }
  virtual int GetFailureCount() const { return fFailureCount; }
  virtual void Reset() { SetStatus(kOk); }
  ENa61ModuleState GetState() const { return fState; }
  int GetStatus() const { return fStatus; }

  // Information and debug
  virtual void SetDebugLevel(const int level) { fDebugLevel = level; }
  int DebugLevel() { return fDebugLevel; }  // Debuglevel value
  virtual void SetVerbose(const int verbose) { fVerbose = verbose; }
  int Verbose() { return fVerbose; }
  virtual void Print(Option_t* option = "DB");
  virtual void Print(Option_t* option = "DB") const;

  // Miscelaneous information
  virtual int GetEventNumber(UEvent* ev) const;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Int_t fNDev;
  TString fDevName[100];     // array size should be updated later 

  Int_t   fNDevComb;
  TString fDevCombName[200]; // array size should be updated later 
  Int_t fMatchId1[200]; // array size should be updated later 
  Int_t fMatchId2[200]; // array size should be updated later 

  Int_t fDevi[200];
  Int_t fDevj[200];

  //// arrays to keep "string name" to "sensor position" information
  Int_t fVdsArr[100][3];
  Int_t fColArr[100][3];
  Int_t fSenArr[100][3];

  void AnaSensorString(TString stringName,Int_t* ic,Int_t* is);

  //  ClassDef(Na61Module,1)  // BRAHMS Module definitions
};

#endif
