// -*- mode: c++ -*-
//
///////////////////////////////////////////////////////////////////////
//                                                                   //
//    UEvent                                                         //
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
//    Modified by Paweł Staszel (for NA61 apptication)               //
//                                                                   //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef UTIL_UEvent
#define UTIL_UEvent
#include "TNamed.h"
#include "THashTable.h"
#include "TBrowser.h"
#include "UEventNode.h"
#include "UVector3D.h"

#define UEvent UEvent_ // hack due to clash with UniGen

class UEvent : public UEventNode {
 public:
  UEvent();
  UEvent(const char* Name, int run, int event);
  virtual ~UEvent();

  virtual int GetRunNumber() { return fRunNumber; }
  virtual int GetEventNumber() { return fEventNumber; }
  virtual double GetEventVertexZ() { return fEventVertexZ; }

  virtual double GetPrimaryVertexX() { return fPrimaryVertexX; }
  virtual double GetPrimaryVertexY() { return fPrimaryVertexY; }
  virtual double GetPrimaryVertexZ() { return fPrimaryVertexZ; }

  virtual double GetPrimaryRecoVertexX() { return fPrimaryRecoVertexX; };
  virtual double GetPrimaryRecoVertexY() { return fPrimaryRecoVertexY; };
  virtual double GetPrimaryRecoVertexZ() { return fPrimaryRecoVertexZ; };

  bool GetPrimaryRecoVertexSet() { return fPrimaryRecoVertexSet; }

  void SetPrimaryVertex(double vx, double vy, double vz) {
    fPrimaryVertexX = vx;
    fPrimaryVertexY = vy;
    fPrimaryVertexZ = vz;
  }

  void SetPrimaryRecoVertex(double vx, double vy, double vz) {
    fPrimaryRecoVertexX = vx;
    fPrimaryRecoVertexY = vy;
    fPrimaryRecoVertexZ = vz;
    fPrimaryRecoVertexSet = true;
  }

  void ResetPrimaryRecoVertex() { fPrimaryRecoVertexSet = false; }

  void SetEventVertexZ(double v) { fEventVertexZ = v; }

  virtual void Print(Option_t* opt = "");  //*MENU*
  virtual void Copy(UEvent& event);
  virtual void Browse(TBrowser* b);

 private:
  double fEventVertexZ;
  int fRunNumber;    //! non persistent!
  int fEventNumber;  //! non persistent!
  double fPrimaryVertexX;
  double fPrimaryVertexY;
  double fPrimaryVertexZ;
  double fPrimaryRecoVertexX;
  double fPrimaryRecoVertexY;
  double fPrimaryRecoVertexZ;
  bool fPrimaryRecoVertexSet;  //! non persistent!
 public:
  //  ClassDef(UEvent,4) // G4NA61 event data class
};

#endif
