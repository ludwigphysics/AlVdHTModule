// -*- mode: c++ -*-
//
///////////////////////////////////////////////////////////////////////
//                                                                   //
//    UVdEvent                                                         //
//                                                                   //
//    G4NA61 event class                                             //
//                                                                   //
//    Na61Event manages access to the raw and reconstructed data     //
//    for one event in the G4NA61.                                   //
//    environment                                                    //
//                                                                   //
//    Original Author: Kris Hagel from a model by Gunther Roland     //
//    Created :                                                      //
//    Version : 1.0                                                  //
//    Changed : July 12 1997
//    Modified by Paweł Staszel (for NA61 application)               //
//                                                                   //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef UTIL_UVdEvent
#define UTIL_UVdEvent
#include "TNamed.h"
#include "THashTable.h"
#include "TBrowser.h"
#include "UEventNode.h"
#include "UVector3D.h"

class UVdEvent : public UEventNode {
 public:
  UVdEvent();
  UVdEvent(const char* Name, int run, int event);
  virtual ~UVdEvent();

  virtual unsigned int GetSensorCounter(int i) { return fSensorCounter[i]; }
  virtual bool GetEventIsCorrupted() { return fCorruptedEvent; }
  virtual bool GetTriggerFlagSet() { return fTriggerFlagSet; }
  virtual int GetPrimaryVertexStatus() { return fPrimaryVertexStatus; }
  virtual int GetRunNumber() { return fRunNumber; }
  virtual int GetEventNumber() { return fEventNumber; }
  virtual int GetNFrames() { return fNFrames; }
  virtual int GetTimer(int iff, int i) {  // iff=fpga#, i=frame#
    if (i > fNFrames) Error("GetTimer1", "Wrong wrong frame number, chech it out");
    return fTimer[iff][i];
  }
  virtual int GetJuraTimer(int i) {
    if (i > fNFrames) Error("GetJuraTimer", "Wrong frame number, chech it out");
    return fJuraTimer[i];
  }
  virtual int GetSaleveTimer(int i) {
    if (i > fNFrames) Error("GetSaleveTimer", "Wrong frame number, chech it out");
    return fSaleveTimer[i];
  }

  virtual double GetPrimaryVertexX() { return fPrimaryVertexX; }
  virtual double GetPrimaryVertexY() { return fPrimaryVertexY; }
  virtual double GetPrimaryVertexZ() { return fPrimaryVertexZ; }

  virtual int GetPrimaryTracks() { return fPrimaryTracks; }
  virtual int GetAdded4HitTracks() { return fAdded4HitTracks; }

  void SetSensorCounter(int i, unsigned int v) { fSensorCounter[i] = v; }
  void SetTriggerFlagSet(bool v) { fTriggerFlagSet = v; }
  void SetPrimaryVertexStatus(int v) { fPrimaryVertexStatus = v; }
  void SetRunNumber(int v) { fRunNumber = v; }
  void SetEventNumber(int v) { fEventNumber = v; }
  void SetNFrames(int v) { fNFrames = v; }
  void SetTimers(int i, int v1, int v2, int v3) {
    fTimer[0][i] = v1;
    fTimer[1][i] = v2;
    fTimer[2][i] = v3;
  }
  void SetJuraTimer(int i, int v) { fJuraTimer[i] = v; }
  void SetSaleveTimer(int i, int v) { fSaleveTimer[i] = v; }
  void SetEventIsCorrupted(bool v) { fCorruptedEvent = v; }
  void SetPrimaryTracks(int v) { fPrimaryTracks = v; }      // it means number of primary tracks
  void SetAdded4HitTracks(int v) { fAdded4HitTracks = v; }  // it means number of 4hit added tracks

  void SetPrimaryVertex(double vx, double vy, double vz) {
    fPrimaryVertexX = vx;
    fPrimaryVertexY = vy;
    fPrimaryVertexZ = vz;
  }

  virtual void Print(Option_t* opt = "");  //*MENU*
  virtual void Copy(UVdEvent& event);
  virtual void Browse(TBrowser* b);

 private:
  int fPrimaryTracks;
  int fAdded4HitTracks;
  bool fTriggerFlagSet;
  bool fCorruptedEvent;
  int fRunNumber;
  int fEventNumber;
  double fPrimaryVertexX;
  double fPrimaryVertexY;
  double fPrimaryVertexZ;
  int fPrimaryVertexStatus;

  int fNFrames;
  int fTimer[3][10];

  int fJuraTimer[10];
  int fSaleveTimer[10];

  unsigned int fSensorCounter[8];

 public:
  //ClassDef(UVdEvent,9) // G4NA61 event data class
};

#endif
