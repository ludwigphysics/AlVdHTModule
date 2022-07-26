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
#ifndef NA61_Na61AlVdTrackingModule
#define NA61_Na61AlVdTrackingModule

#ifndef NA61_Na61Module
#include "Na61Module.h"
#endif
#ifndef NA61_Na61AlVdArmParameters
#include "Na61AlVdArmParameters.h"
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

class Na61AlVdTrackingModule : public Na61Module
{ 
  
public:
  
  Na61AlVdTrackingModule();
  Na61AlVdTrackingModule(const Char_t* name, const Char_t* title);
  virtual ~Na61AlVdTrackingModule () {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNode, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option="B") const; // *MENU*


  void SetOutTableName(TString v){fOutTableName = v;}
  void SetHistDirName(TString v){fHistDirName = v;}
  virtual void SetProductionMode(Int_t v){fProductionMode=v;}
  void SetRunId(Int_t v){fRunId=v;}
  void SetCutId(Int_t v){fCutId = v;}
  void SetNsig_dev(Float_t v) {fNsig_dev = v;}
  void SetNsig_pvert(Float_t v) {fNsig_pvert = v;}

  Int_t GetAll3StationTracks(Int_t i){ return fAll3StationTracks[i];}
  Int_t GetAllFullTracks(Int_t i){ return fhDx_match_acce[i]->GetEntries();}
  Vector3D* GetPrimaryVertex(){return fPrimaryVertex;}
  //TH1F* GetfhFullTracksAll(){return fhFullTracksAll;}
  Bool_t GetFieldRun() {return fFieldRun;}


  void SetToJuraArm(Bool_t v){fJura = v;}

public:

  //void FitLine_w2(USensorHit** hits,Float_t* sig, Float_t& ax,Float_t& ay,Float_t& bx,Float_t& by,
		 // Float_t& chi2x,Float_t& chi2y,Float_t& N, Float_t& zmin);

  void FitLine_w2(USensorHit** hits, double* sig, double& ax, double& ay, double& bx, double& by, 
  double& chi2x, double& chi2y, double& N, double& zmin);

  void FitLineY_w2(USensorHit** hits,Float_t* sig, Float_t& ay,Float_t& by,Float_t& N, Float_t& zmin);

  void CombineTracks(Int_t imatch, Int_t idev1, Int_t idev2, UEventNode* outNode);
  //void CombineWithPrimaryVertex(TString matchstr, UEventNode* outNode);

  //void CreateFullTrack(Int_t imatch,UVdTrack* track1, UVdTrack* track2, UDataTable* tracktab);
  void CreateFullTrack_new(Int_t imatch,UVdTrack* track1, UVdTrack* track2, UDataTable* tracktab);
  //void CreateFullTrackx_new(Int_t imatch, UVdTrack* track, UDataTable* tracktab);

  //void FindPrimaryVertexAtd(UEventNode* out, TString* devname);
  void FindPrimaryVertex(UEventNode* outNode);
  void PrimaryVertexWithWeigths(TObjArray& tracks,Float_t zprim,Vector3D& pvertex);

  void FindPrimaryVertexPost(UEventNode* outNode);

  void Make3HitTracks(Int_t idev, UEventNode* out, Int_t* TabIndArray);
  void Make3HitTracks(TString devname, Int_t csize, UDataTable* hits1, UDataTable* hits2, UDataTable* hits3, UEventNode* outNode, Int_t* tabInd);
  void Create3HitTrack(Int_t idev,USensorHit* hit1, USensorHit* hit2, USensorHit* hit3, Int_t* IndArray, Int_t* TabIndArray, UDataTable* tracktab);
  //void Create3HitTrackAl(Int_t idev,TString devname, USensorHit* hit1, USensorHit* hit2, USensorHit* hit3, Int_t* IndArray, Int_t* TabIndArray, UDataTable* tracktab);
  void CheckTrackMatching(UVdTrack* track1, UVdTrack* track2);
  void CheckTrackMatching2(UVdTrack* track1, UVdTrack* track2);
  void RemoveGhostTracks(TString devname1, TString devname2, UEventNode* out);


  void Delete3HitTracks(UEventNode* outNode);
  Int_t GetAllTracks(UEventNode* out);
  void FillHitPositions(UEventNode* out);


  void RefitTracks(UEventNode* out);
  void RefitTracks(Int_t Ntab,UEventNode* out);

  //_____________________________________________________________________________
  Int_t GetAlSensorNumber(Int_t station, Int_t column, Int_t number_in_column)
  {
    Int_t station_fact = 0;
    Int_t column_fact  = 0;
    if(station==2) {station_fact = 3;  column_fact = 3;} 
    if(station==3) {station_fact = 9;  column_fact = 5;}
    if(station==4) {station_fact = 19; column_fact = 5;}
    
    return station_fact + column*column_fact + number_in_column;
  }
  //_____________________________________________________________________________
  Int_t GetAlSensorNumber2(Int_t station, Int_t column, Int_t number_in_column)
  {
    Int_t station_fact = 0;
    Int_t column_fact  = 0;
    Int_t offset=1;
    if(station==2) {station_fact = 3;  column_fact = 3; offset=1;} 
    if(station==3) {station_fact = 9;  column_fact = 5; offset=2;}
    if(station==4) {station_fact = 19; column_fact = 5; offset=2;}
    
    return station_fact + column*column_fact + number_in_column + offset;

  }

  //_______________________________________________________________________________
  Int_t GetAlSensorNumber(Int_t station, Int_t number_in_station)
  {
    Int_t station_fact = 0;
    if(station==2) {station_fact = 3;} 
    if(station==3) {station_fact = 9;}
    if(station==4) {station_fact = 19;}
    
    return station_fact + number_in_station;
  }
  
  
private:

  void CalculateSums(Float_t& sum_up_x, Float_t& sum_do_x, Float_t& sum_up_y, Float_t& sum_do_y, TObjArray tracks);
  void RefitTrack_front(UVdTrack* track);
  void RefitTrack_back(UVdTrack* track);

  void FillDeviations(Int_t imatch, USensorHit* hit1, USensorHit* hit2, USensorHit* hit3);

public:

  TString fHistDirName;
  TString fOutTableName;

  UDataTable* fHitsTabs[34];

  Int_t fProductionMode;
  Bool_t fPrimaryVertexDefined;
  Vector3D* fPrimaryVertex;
  Na61AlVdArmParameters* fParams;

  TF1* fFGauss;

  Bool_t fLargeClusters;
  Bool_t fJura;
  Int_t fCutId;
  Int_t fAll3StationTracks[20];

  Int_t fRunId;

  // public histos
  TH1F* fhFullTracks;
  TH1F* fhFullTracksAll;


private:

  Float_t fZprim;
  Bool_t fFieldRun;

  Int_t fTrackID;

  Float_t fNsig_dev;
  Float_t fNsig_pvert;

  TH2F* fhZX_fine;
  TH1F* fhX_Al1;
  TH1F* fhX_Al2;
  TH1F* fhX_Al3;
  TH1F* fhX_Al4;
  TH1F* fhY_Al1;
  TH1F* fhY_Al2;
  TH1F* fhY_Al3;
  TH1F* fhY_Al4;


  TH1F* fhCurvature_front;  
  TH1F* fhCurvature_back;    
  TH1F* fhCurvature_back_prod;    
  TH1F* fhSlopeChange;  
  TH1F* fdx3hit[4];
  TH1F* fdx4hit_front[4];
  TH1F* fdx4hit_back[4];


  TH1F* fhDZmin;
  TH1F* fhTracks_prod;
  TH1F* fhTracks2_prod;
  TH1F* fhChi2Ndf_x;
  TH1F* fhChi2Ndf_y;



  TH1F* fhX_dev[100];
  TH1F* fhY_dev[100];
  TH1F* fhX_dev_cuty[100];
  TH1F* fhY_dev_cutx[100];

  TH1F* fhX_dev_acce[100];
  TH1F* fhY_dev_acce[100];
  TH2F* fh_devx_devy[100];


  TH1F* fhChi2X[100];
  TH1F* fhChi2Y[100];

  TH1F* fhAxAtd;
  TH1F* fhAyAtd;
  TH1F* fhAxAtd_prod;
  TH1F* fhAyAtd_prod;

  TH1F* fhAx[20];
  TH1F* fhAy[20];

  TH2F* fhNsVsNe[20];
  TH1F* fhTracks[20];


  TH1F* fhClusterSize;

  TH1F* fhDx_match[200];
  TH1F* fhDy_match[200];
  TH1F* fhDax_match[200];
  TH1F* fhDay_match[200];
  TH1F* fhDx_match_cut1[200];
  TH1F* fhDy_match_cut1[200];
  TH1F* fhDax_match_cut1[200];
  TH1F* fhDay_match_cut1[200];
  TH1F* fhDx_match_acce[200];
  TH1F* fhDy_match_acce[200];
  TH1F* fhDax_match_acce[200];
  TH1F* fhDay_match_acce[200];


  TH1F* fhAx_full[200];
  TH1F* fhAy_full[200];
  TH1F* fhAy_full_prod[200];

  TH1F* fhX_dev_full[200];
  TH1F* fhY_dev_full[200];

  TH1F* fhX_dev_Vds1[200];
  TH1F* fhX_dev_Vds2[200];
  TH1F* fhX_dev_Vds3[200];
  TH1F* fhX_dev_Vds4[200];

  TH1F* fhY_dev_Vds1[200];
  TH1F* fhY_dev_Vds2[200];
  TH1F* fhY_dev_Vds3[200];
  TH1F* fhY_dev_Vds4[200];


  TH1F* fhRecoVertexZ;
  TH1F* fhRecoVertexZ_fine;

  TH1F* fhRecoVertexZ_x;
  TH1F* fhRecoVertexZ_y;
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
  
  //  ClassDef(Na61AlVdTrackingModule,0) // DOCUMENT ME
  
};

#endif
//____________________________________________________________________
//
