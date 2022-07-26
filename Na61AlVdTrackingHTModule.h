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
// Some new changes: Mateusz Bajda (start in 03_05_22)
//
#ifndef NA61_Na61AlVdTrackingHTModule
#define NA61_Na61AlVdTrackingHTModule

//#ifndef NA61_Na61Module
//#include "Na61Module.h"
//#endif

#ifndef NA61_Na61AlVdTrackingModule
#include "Na61AlVdTrackingModule.h"
#endif

#ifndef NA61_Na61VdParameters
#include "Na61VdParameters.h"
#endif

#ifndef UTIL_UEventNode
#include "UEventNode.h"
#endif
#ifndef UTIL_UVdEvent
#include "UVdEvent.h"
#endif
#ifndef UTIL_USensorHit
#include "USensorHit.h"
#endif
#ifndef UTIL_UHitInfo
#include "UHitInfo.h"
#endif
#ifndef UTIL_UVdTrack
#include "UVdTrack.h"
#endif
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
#ifndef ROOT_TObjArray
#include "TObjArray.h"
#endif

#ifndef AccSteppingAction_h
#define AccSteppingAction_h
#endif
#ifndef ROOT_TGraph
#include "TGraph.h"
#endif

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

class Na61AlVdTrackingHTModule : public Na61AlVdTrackingModule {
 public:
  Na61AlVdTrackingHTModule();
  Na61AlVdTrackingHTModule(const char* name, const char* title);
  virtual ~Na61AlVdTrackingHTModule() {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option = "B") const;  // *MtENU*
  void SetDistCut(double v) { fDistCut = v; }
  void SetVtpcSetup(int v) { fOldSetup = v; }
  void SetFillHitOccupancy(int v) { fFillHitOccupancy = v; }
  void SetAddInGrid2Tracks(int v) { fAddInGrid2Tracks = v; }
  void SetFiducialCut(double v) { fFiduCut = v; }
  void SetPosRes(double v) { fPosRes = v; }
  void Setbw(double vx, double vy) {
    fbwx = vx;
    fbwy = vy;
  }
  void SetProduction(int v) { fProduction = v; }
  void SetMakeClusters(int v) { fMakeClusters = v; }
  void SetChi2Cut(double v) { fChi2Cut = v; }
  void SetVzOffset(double v) { fVzOffset = v; }
  void SetAdd4HitTracks(bool v) { fAdd4HitTracks = v; }

  int GetAllTracks() { return fAllTracks; }
  int GetAllHtTracks() { return fAllHtTracks; }

 private:
  TObjArray* GetHits(TString vdsStr, UEventNode* inJ, UEventNode* inS);

  void CalculateDeviations(USensorHit* hit1, USensorHit* hit2, USensorHit* hit3, USensorHit* hit4, double& devx, double& devy);
  void CalculateResDeviations(USensorHit** hits, double* sig, double* resx, double* resy);
  void Add4HitTracks(UEventNode* inNodeS, UEventNode* inNodeJ, UEventNode* outNode);
  void SetTrackHits(UVdTrack* track, UEventNode* inNode);
  void FillCombinatorialCurvature(UEventNode* inNode);

  void CheckRing(int i, int j);
  void CheckRingNew(int i, int j, int NRes);
  void RecoWithSimultaneousHT(UEventNode* inNodeS, UEventNode* inNodeJ, UEventNode* outNode);
  void FillHSaxayWithClusters();

  bool CheckAcceptance(USensorHit* hit, UDataTable* tabVTPC1, UDataTable* tabVTPC2);
  bool CheckAcceptance_Step1(USensorHit* hit1_vtpc1L, UDataTable* tabVTPC1_L, UDataTable* tabVTPC1_R, UDataTable* tabVTPC2_L, UDataTable* tabVTPC2_R);
  bool CheckAcceptance_Step1New(USensorHit* hit1_vtpc1L, UDataTable* tabVTPC1_R, UDataTable* tabVTPC2_L, UDataTable* tabVTPC2_R);
  bool CheckAcceptance_Step2(USensorHit* hit1_vtpc1R, UDataTable* tabVTPC1_R, UDataTable* tabVTPC2_L, UDataTable* tabVTPC2_R);
  bool CheckAcceptance_Step2New(USensorHit* hit1_vtpc1R, UDataTable* tabVTPC2_L, UDataTable* tabVTPC2_R);
  bool CheckAcceptance_Step3(USensorHit* hit1_vtpc2L, UDataTable* tabVTPC2_L, UDataTable* tabVTPC2_R);
  bool CheckAcceptance_Step3New(USensorHit* hit1_vtpc2L, UDataTable* tabVTPC2_R);
  bool CheckAcceptance_Step4(USensorHit* hit1_vtpc2R, UDataTable* tabVTPC2_R);
  bool CheckAcceptance_Step4New(USensorHit* hit1_vtpc2R);

  // void CreateAndDefineRecoTrack(USensorHit* hit,UDataTable* tabInGrid1,UDataTable* tabInGrid2,UDataTable* table);
  void CheckforClusterization();
  void MakeClusters(int Mini, UDataTable** tabArr);
  void MakeClustersNew();

  void CheckTrackMatching(UVdTrack* track1, UVdTrack* track2);
  bool CheckTrackMatching2(UVdTrack* track1, UVdTrack* track2);

  void AnaClustersAndCreateTracks(UDataTable* table, UDataTable** tabArr);
  void AnaClustersAndCreateTracksNew(UDataTable* table, TObjArray** tabArr);

  void ExtractTracksFromCluster(TObjArray& hitsVds1, TObjArray& hitsVds2, TObjArray& hitsVds3, TObjArray& hitsVds4, UDataTable* table);

  void FitAndCreateTrack(USensorHit** hits, UDataTable* table);
  int ClearHits4(USensorHit** hits, int* mVds);


 private:
  Na61VdParameters* fVdParams;

  bool fAdd4HitTracks;  // Addeded to handle TrackingHTPackage, should be set to TRUE only for the initial module (--vzoff=0)
  int fAllTracks;
  int fAllHtTracks;
  int fSignalTracks;  // signal tracks counter
  int fMakeClusters;

  TObjArray fCluster;
  int fResIndx[20];
  int fResIndy[20];
  int fNRes;

  TObjArray* fHitInfoArray;
  unsigned int faxayToInd[1000][1000];

  int fProduction;
  int fclustercheck;

  int fMaxXInd;
  int fMaxYInd;
  int findx;
  int findy;
  int fi, fj;  // indices showig location of the seed cluster
  int fRingLevel;

  double fChi2Cut;
  double fVzOffset;

  double fPosRes;
  double fbwx;
  double fbwy;
  double fFiduCut;
  double fDistCut;
  int fOldSetup;
  int fFillHitOccupancy;
  int fAddInGrid2Tracks;

  TH1F* fhDevX[4];
  TH1F* fhDevY[4];
  TH1F* fhDevX_acce[4];
  TH1F* fhDevY_acce[4];

  TH1F* fhResX[4];
  TH1F* fhResY[4];
  TH1F* fhResX_acce[4];
  TH1F* fhResY_acce[4];

  TH1F* fhCurvatureComb_back;
  TH1F* fhPrimaryVertexZ;
  TH1F* fhPrimaryVertexZ_acce;

  TH1F* fhZminN4;
  TH1F* fhZminN4Match;
  TH1F* fhZminN3;
  TH1F* fhZminN3Match;

  TH1F* fhPosSpreadX;
  TH1F* fhPosSpreadY;

  TH1F* fhHitsVds1;
  TH1F* fhHitsVds2;
  TH1F* fhHitsVds3;
  TH1F* fhHitsVds4;

  TH1F* fhAdded4hitTracks;

  TH1F* fhChi2xNdf_N3;
  TH1F* fhChi2yNdf_N3;
  TH1F* fhChi2xNdf_N4;
  TH1F* fhChi2yNdf_N4;

  TH1F* fhChi2Ndf;
  TH1F* fhChi2xNdf;
  TH1F* fhChi2yNdf;

  TH1F* fhChi2Ndf_Acce;
  TH1F* fhChi2xNdf_Acce;
  TH1F* fhChi2yNdf_Acce;

  TH2D* fhHSaxay;
  TH2D* fhHSaxay_clust;
  TH2F* fhYPt;

  TH1F* fhRecoTracks;
  TH1F* fhRecoHtTracks;

  TGraph* fGr_axay_vdtrack;

  // TH1F* fhPz2;

 public:
  //  ClassDef(Na61aLVdTrackingHTModule,0) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
