// -*- mode: c++ -*-
//
// Module for primary track reconstruction using both VD arms
//
// Author: Paweł Staszel
//
//
#ifndef NA61_Na61PrimaryVertexRecoHTModule
#define NA61_Na61PrimaryVertexRecoHTModule

#ifndef NA61_Na61PrimaryVertexRecoModule
#include "Na61PrimaryVertexRecoModule.h"
#endif

#ifndef NA61_Na61Module
#include "Na61Module.h"
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

class Na61PrimaryVertexRecoHTModule : public Na61PrimaryVertexRecoModule {
 public:
  Na61PrimaryVertexRecoHTModule();
  Na61PrimaryVertexRecoHTModule(const char* namec, const char* title);
  virtual ~Na61PrimaryVertexRecoHTModule() {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode);
  virtual void Event(UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option = "B") const;  // *MENU*

  void SetTrackInput(int v) { fTrackInput = v; }
  void SetZprim(double v) { fZprim = v; }
  void SetZCut(double v) { fZCut = v; }

 private:
  void ClassifyTracksHT(UEventNode* out);
  void FindPrimaryVertexPostHT(int stage, int flag, UEventNode* out);
  void PrimaryVertexWithWeigths_Y(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex);
  void PrimaryVertexWithWeigths_X(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex);
  void PrimaryVertexWithWeigths_XY(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex, double sigx, double sigy);
  void PrimaryVertexWithWeigths_XY2(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex, double sigx, double sigy);
  void CheckPrimaryVertexResolutionHT(int stage, UEventNode* out);

  void MethodHT(UEventNode* outNode);
  void TagPrimaryTracks(int stage, UEventNode* outNode);

  void Find2TrackVertexes(UEventNode* out);

 private:
  int fTrackInput;
  int fStage;
  
  double fSigx;  // parameters used to weight x vs y
  double fSigy;


  double fZprim;
  double fZCut;

  // Vector3D fPrimaryVertex; // three stages
  Vector3D fPrimaryVertex[3];  // three stages
  Vector3D fPrimaryVertexJ;
  Vector3D fPrimaryVertexS;

  TH1F* fhVtxDx;
  TH1F* fhVtxDy;
  TH1F* fhVtxDz;

  TH1F* fhVtxDx_Glob;
  TH1F* fhVtxDy_Glob;
  TH1F* fhVtxDz_Glob;

  TH1F* fhVtxDx_fine_Glob;
  TH1F* fhVtxDy_fine_Glob;
  TH1F* fhVtxDz_fine_Glob;

  TH1F* fhVtxDx_Glob2[3];
  TH1F* fhVtxDy_Glob2[3];
  TH1F* fhVtxDz_Glob2[3];

  TH1F* fhVtxDx_fine_Glob2[3];
  TH1F* fhVtxDy_fine_Glob2[3];
  TH1F* fhVtxDz_fine_Glob2[3];

  TH2F* fhVtxDx_fine_vs_multi;
  TH2F* fhVtxDy_fine_vs_multi;
  TH2F* fhVtxDz_fine_vs_multi;

  TH2F* fhVtxDx_fine_vs_multi_La[3];
  TH2F* fhVtxDy_fine_vs_multi_La[3];
  TH2F* fhVtxDz_fine_vs_multi_La[3];

  TH2F* fhDx_2TrackVertexRes;
  TH2F* fhDy_2TrackVertexRes;
  TH2F* fhDz_2TrackVertexRes;

  TH2F* fhVtxStatus_JuraVsSaleve;
  TH2F* fhFullTracks_JvsS;
  TH1F* fhFullTracks_JS;
  TH1F* fhFullTracks_J;
  TH1F* fhFullTracks_S;
  TH1F* fhFullTracks4h_JS;
  TH1F* fhFullTracks4h_J;
  TH1F* fhFullTracks4h_S;

  TH1F* fhRecoVertexZ;  // after last stage
  TH1F* fhRecoVertexZ_fine;
  TH1F* fhRecoVertexZ_xfine;
  TH2F* fhRecoVertexXY;

  TH1F* fhRecoVertexZ_J;  // after last stage
  TH1F* fhRecoVertexZ_fine_J;
  TH2F* fhRecoVertexXY_J;

  TH1F* fhRecoVertexZ_S;  // after last stage
  TH1F* fhRecoVertexZ_fine_S;
  TH2F* fhRecoVertexXY_S;

  TH1F* fhAxJ;
  TH1F* fhAyJ;
  TH1F* fhAxS;
  TH1F* fhAyS;

  TH1F* fhVtxX;
  TH1F* fhVtxY;
  TH1F* fhTracksVsVtxX_J;
  TH1F* fhTracksVsVtxY_J;
  TH1F* fhTracksVsVtxX_S;
  TH1F* fhTracksVsVtxY_S;

  TH1F* fhTrackDx[3];
  TH1F* fhTrackDy[3];
  TH1F* fhTrackDz[3];

  TH1F* fhTrackDx_tag[3];
  TH1F* fhTrackDy_tag[3];
  TH1F* fhTrackDz_tag[3];

 public:
  //  ClassDef(Na61PrimaryVertexRecoHTModule,1) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
