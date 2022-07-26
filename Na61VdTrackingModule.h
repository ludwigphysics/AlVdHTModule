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
#ifndef NA61_Na61VdTrackingModule
#define NA61_Na61VdTrackingModule

#ifndef NA61_Na61Module
#include "Na61Module.h"
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

class Na61VdTrackingModule : public Na61Module {
 public:
  Na61VdTrackingModule();
  Na61VdTrackingModule(const char* name, const char* title);
  virtual ~Na61VdTrackingModule() {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNode, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option = "B") const;  // *MENU*

  void SetOutTableName(TString v) { fOutTableName = v; }
  void SetHistDirName(TString v) { fHistDirName = v; }
  virtual void SetProductionMode(int v) { fProductionMode = v; }
  void SetRunId(int v) { fRunId = v; }
  void SetCutId(int v) { fCutId = v; }
  void SetNsig_dev(double v) { fNsig_dev = v; }
  void SetNsig_pvert(double v) { fNsig_pvert = v; }

  int GetAll3StationTracks(int i) { return fAll3StationTracks[i]; }
  int GetAllFullTracks(int i) { return fhDx_match_acce[i]->GetEntries(); }
  Vector3D* GetPrimaryVertex() { return fPrimaryVertex; }
  // TH1F* GetfhFullTracksAll(){return fhFullTracksAll;}
  bool GetFieldRun() { return fFieldRun; }

  void SetToJuraArm(bool v) { fJura = v; }

 public:
  void FitLine_w2(USensorHit** hits, double* sig, double& ax, double& ay, double& bx, double& by, double& chi2x, double& chi2y, double& N, double& zmin);

  void FitLineY_w2(USensorHit** hits, double* sig, double& ay, double& by, double& N, double& zmin);

  void CombineTracks(TString matchname, int idev1, int idev2, UEventNode* outNode);
  void CombineWithPrimaryVertex(TString matchstr, UEventNode* outNode);

  void CreateFullTrack(int imatch, UVdTrack* track1, UVdTrack* track2, UDataTable* tracktab);
  void CreateFullTrack_new(int imatch, UVdTrack* track1, UVdTrack* track2, UDataTable* tracktab);
  void CreateFullTrackx_new(int imatch, UVdTrack* track, UDataTable* tracktab);

  void FindPrimaryVertex(UEventNode* outNode);
  void PrimaryVertexWithWeigths(TObjArray& tracks, double zprim, Vector3D& pvertex);

  void FindPrimaryVertexPost(UEventNode* outNode);

  void Make3HitTracks(TString devname, int csize, UDataTable* hits1, UDataTable* hits2, UDataTable* hits3, UEventNode* outNode, int* tabInd);
  void Create3HitTrack(int idev, USensorHit* hit1, USensorHit* hit2, USensorHit* hit3, int* IndArray, int* TabIndArray, UDataTable* tracktab);
  void CheckTrackMatching(UVdTrack* track1, UVdTrack* track2);
  void CheckTrackMatching2(UVdTrack* track1, UVdTrack* track2);
  void RemoveGhostTracks(TString devname1, TString devname2, UEventNode* out);

  void Delete3HitTracks(UEventNode* outNode);
  int GetAllTracks(UEventNode* out);

  void RefitTracks(UEventNode* out);
  void RefitTracks(int Ntab, UEventNode* out);

 private:
  void CalculateSums(double& sum_up_x, double& sum_do_x, double& sum_up_y, double& sum_do_y, TObjArray tracks);
  void RefitTrack_front(UVdTrack* track, int iflag = 0);
  void RefitTrack_back(UVdTrack* track);

  void FillDeviations(int imatch, USensorHit* hit1, USensorHit* hit2, USensorHit* hit3);

 public:
  TString fHistDirName;
  TString fOutTableName;

  UDataTable* fHitsTabs[8];

  int fProductionMode;
  bool fPrimaryVertexDefined;
  Vector3D* fPrimaryVertex;
  Na61ArmParameters* fParams;

  TF1* fFGauss;

  bool fLargeClusters;
  bool fJura;
  int fCutId;
  int fAll3StationTracks[20];

  int fRunId;

  // public histos
  TH1F* fhFullTracks;
  TH1F* fhFullTracksAll;

 private:
  double fZprim;
  bool fFieldRun;

  int fTrackID;

  double fNsig_dev;
  double fNsig_pvert;

  TH1F* fhCurvature_front;
  TH1F* fhCurvature_back;
  TH1F* fhSlopeChange;
  TH1F* fdx3hit[4];
  TH1F* fdx4hit_front[4];
  TH1F* fdx4hit_back[4];

  TH1F* fhDZmin;

  TH1F* fhChi2Ndf_x;
  TH1F* fhChi2Ndf_y;

  TH2F* fh_devx_devy[20];

  TH1F* fhX_dev[20];
  TH1F* fhY_dev[20];
  TH1F* fhX_dev_cuty[20];
  TH1F* fhY_dev_cutx[20];

  TH1F* fhX_dev_acce[20];
  TH1F* fhY_dev_acce[20];

  TH1F* fhChi2X[20];
  TH1F* fhChi2Y[20];

  TH1F* fhAx[20];
  TH1F* fhAy[20];

  TH2F* fhNsVsNe[20];
  TH1F* fhTracks[20];

  TH1F* fhX_dev_lc[6];
  TH1F* fhY_dev_lc[6];
  TH1F* fhX_dev_cuty_lc[6];
  TH1F* fhY_dev_cutx_lc[6];

  TH1F* fhClusterSize[6];
  TH1F* fhClusterSizeInPick[6];

  TH1F* fhAx_lc[6];
  TH1F* fhAy_lc[6];

  TH1F* fhDx_match[4];
  TH1F* fhDy_match[4];
  TH1F* fhDax_match[4];
  TH1F* fhDay_match[4];
  TH1F* fhDx_match_cut1[4];
  TH1F* fhDy_match_cut1[4];
  TH1F* fhDax_match_cut1[4];
  TH1F* fhDay_match_cut1[4];
  TH1F* fhDx_match_acce[4];
  TH1F* fhDy_match_acce[4];
  TH1F* fhDax_match_acce[4];
  TH1F* fhDay_match_acce[4];

  TH1F* fhAx_full[18];
  TH1F* fhAy_full[18];
  TH1F* fhAy_full_prod[18];

  TH1F* fhX_dev_full[18];
  TH1F* fhY_dev_full[18];

  TH1F* fhX_dev_Vds1[18];
  TH1F* fhX_dev_Vds2[18];
  TH1F* fhX_dev_Vds3[18];
  TH1F* fhX_dev_Vds4[18];

  TH1F* fhY_dev_Vds1[18];
  TH1F* fhY_dev_Vds2[18];
  TH1F* fhY_dev_Vds3[18];
  TH1F* fhY_dev_Vds4[18];

  TH1F* fhDAx[18];
  TH1F* fhDAy[18];
  TH1F* fhDAx_cuty[18];
  TH1F* fhDAy_cutx[18];

  TH1F* fhDx[18];
  TH1F* fhDy[18];
  TH1F* fhDx_acce[18];
  TH1F* fhDy_acce[18];

  TH1F* fhRecoVertexZ;
  TH1F* fhRecoVertexZ_fine;
  TH2F* fhRecoVertexXY;

  TH1F* fhRecoVertexZ_fine_x;
  TH1F* fhRecoVertexZ_fine_y;

  TH1F* fhd2;

  TH1F* fhRecoVertexZ_w2;
  TH1F* fhRecoVertexZ_fine_w2;
  TH2F* fhRecoVertexXY_w2;

  TH1F* fhRecoVertexZ_w2_Post;
  TH1F* fhRecoVertexZ_fine_w2_Post;
  TH2F* fhRecoVertexXY_w2_Post;

 public:
  //  ClassDef(Na61VdTrackingModule,0) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
