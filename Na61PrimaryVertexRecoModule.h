// -*- mode: c++ -*-
//
// Module for primary track reconstruction using both VD arms
//
// Author: Paweł Staszel
//
//
#ifndef NA61_Na61PrimaryVertexRecoModule
#define NA61_Na61PrimaryVertexRecoModule

#ifndef NA61_Na61Module
#include "Na61Module.h"
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

class Na61PrimaryVertexRecoModule : public Na61Module {
 public:
  Na61PrimaryVertexRecoModule();
  Na61PrimaryVertexRecoModule(const char* namec, const char* title);
  virtual ~Na61PrimaryVertexRecoModule() {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* outNodeJ, UEventNode* outNodeS);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option = "B") const;  // *MENU*

  void SetRunId(int v) { fRunId = v; }
  void SetPrimaryVertex(Vector3D v) { fPrimaryVertex = v; }
  void SetPrimaryVertexStatus(int v) { fPrimaryVertexStatus = v; }
  void SetZPrim(double v){fZprim = v;}

  Vector3D GetPrimaryVertex() { return fPrimaryVertex; }
  int GetPrimaryVertexStatus() { return fPrimaryVertexStatus; }

  virtual int GetAllTracksJura() { return fAllTracksJ; }
  virtual int GetAllTracksSaleve() { return fAllTracksS; }
  virtual int GetAllTracks() { return fAllTracksJ + fAllTracksS; }

  /// localy used function

  void ClassifyTracks(UEventNode* outJ, UEventNode* outS);
  void ClassifyTracks_pPb(UEventNode* outJ, UEventNode* outS);
  void ClassifyTracksAl(UEventNode* outJ, UEventNode* outS);

  void FindPrimaryVertex(int flag, UEventNode* outJ, UEventNode* outS);
  void FindPrimaryVertexAl(int flag, UEventNode* outJ, UEventNode* outS);
  void PrimaryVertexWithWeigths(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex);

  // void SetupLocalToGlobalShift(); // moved to params
  void LocalToGlobal(int iarm, UEventNode* out);
  void LocalToGlobalAl(int iarm, UEventNode* out);
  void TransformHits(UEventNode* out, double dx, double dy, double dz, double sb, double cb, double sg, double cg, double sa, double ca);
  void TransformHitsAl(UEventNode* out, double dx, double dy, double dz, double sb, double cb, double sg, double cg, double sa, double ca);

  void FillClusterCorrelation(UEventNode* inJ, UEventNode* inS);
  void FillClusterCorrelationAl(UEventNode* inJ, UEventNode* inS);

  void FillHitPositions(Int_t arm, UEventNode* out);


 public:
  double fZprim;

  int fPrimaryVertexStatus;

  Vector3D fPrimaryVertex;
  double fNsigJ;
  double fNsigS;

  double fRotX_J;
  double fRotY_J;
  double fRotX_S;
  double fRotY_S;

  double fOffsetX_J;
  double fOffsetY_J;
  double fOffsetZ_J;

  double fOffsetX_S;
  double fOffsetY_S;
  double fOffsetZ_S;

  double fOffAx_J;
  double fSigAx_J;
  double fOffAy_J;
  double fSigAy_J;
  double fOffAx_S;
  double fSigAx_S;
  double fOffAy_S;
  double fSigAy_S;

  // to allow acces in Na61PrimaryVertexRecoPostModule
  Na61VdParameters* fParams;
  Na61ArmParameters* fArmParamsJ;
  Na61ArmParameters* fArmParamsS;

  int fRunId;


  // localy used members
  int fJuraProdCounter;
  int fSaleveProdCounter;

  int fAllTracksJ;
  int fAllTracksS;

  TH2F* fhZX_fine[2];
  TH2F* fhXY_fine[2];
  TH1F* fhX_Al1[2];
  TH1F* fhX_Al2[2];
  TH1F* fhX_Al3[2];
  TH1F* fhX_Al4[2];
  TH1F* fhY_Al1[2];
  TH1F* fhY_Al2[2];
  TH1F* fhY_Al3[2];
  TH1F* fhY_Al4[2];

  TH1F* fhCurvature_back;

  TH1F* fhVtxDx;
  TH1F* fhVtxDy;
  TH1F* fhVtxDz;

  TH1F* fhVtxDx_Glob;
  TH1F* fhVtxDy_Glob;
  TH1F* fhVtxDz_Glob;

  TH2F* fhVtxStatus_JuraVsSaleve;
  TH2F* fhFullTracks_JvsS;
  TH1F* fhFullTracks_JS;
  TH1F* fhFullTracks_J;
  TH1F* fhFullTracks_S;
  TH1F* fhFullTracks4h_JS;
  TH1F* fhFullTracks4h_J;
  TH1F* fhFullTracks4h_S;

  TH1F* fhRecoVertexZ[3];
  TH1F* fhRecoVertexZ_fine[3];
  TH2F* fhRecoVertexXY[3];

  TH1F* fhRecoVertexZ_flag1[3];
  TH1F* fhRecoVertexZ_fine_flag1[3];
  TH2F* fhRecoVertexXY_flag1[3];

  TH1F* fhAxJ_flag02;
  TH1F* fhAyJ_flag02;

  TH1F* fhAxJ_flag0;
  TH1F* fhAyJ_flag0;
  TH1F* fhAxJ_flag1;
  TH1F* fhAxJ_flag2;
  TH1F* fhAyJ_flag2;

  TH1F* fhAxS_flag02;
  TH1F* fhAyS_flag02;

  TH1F* fhAxS_flag0;
  TH1F* fhAyS_flag0;
  TH1F* fhAxS_flag1;
  TH1F* fhAxS_flag2;
  TH1F* fhAyS_flag2;

  TH1F* fhAxJ;
  TH1F* fhAyJ;
  TH1F* fhAxS;
  TH1F* fhAyS;

  TH2F* fhClusters_JvsS[8];
  TH2F* fhClusters_JmS[8];

  TH1F* fhVtxX;
  TH1F* fhVtxY;
  TH1F* fhTracksVsVtxX_J;
  TH1F* fhTracksVsVtxY_J;
  TH1F* fhTracksVsVtxX_S;
  TH1F* fhTracksVsVtxY_S;

 public:
  //  ClassDef(Na61PrimaryVertexRecoModule,1) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
