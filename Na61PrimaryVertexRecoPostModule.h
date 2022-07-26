// -*- mode: c++ -*-
//
// Module for primary track reconstruction using both VD arms
//
// Author: Paweł Staszel
//
//
#ifndef NA61_Na61PrimaryVertexRecoPostModule
#define NA61_Na61PrimaryVertexRecoPostModule

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

class Na61PrimaryVertexRecoPostModule : public  Na61PrimaryVertexRecoModule
{ 

public:

  Na61PrimaryVertexRecoPostModule();
  Na61PrimaryVertexRecoPostModule(const Char_t* namec, const Char_t* title);
  virtual ~Na61PrimaryVertexRecoPostModule() { }
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option="B") const; // *MENU*

  virtual Int_t GetAllTracksJura() {return fAllTracksJ;}
  virtual Int_t GetAllTracksSaleve() {return fAllTracksS;}
  virtual Int_t GetAllTracks() {return fAllTracksJ + fAllTracksS;}

private:

  void ClassifyTracks(UEventNode* outJ, UEventNode* outS);
  void ClassifyTracks_pPb(UEventNode* outJ, UEventNode* outS);
  void FindPrimaryVertexPost(Int_t flag, UEventNode* outJ, UEventNode* outS);
  void PrimaryVertexWithWeigths_Y(Int_t flag,TObjArray& tracks,Float_t zprim,Vector3D& pvertex);
  void PrimaryVertexWithWeigths_X(Int_t flag,TObjArray& tracks,Float_t zprim,Vector3D& pvertex);
  void PrimaryVertexWithWeigths_XY(Int_t flag,TObjArray& tracks,Float_t zprim,Vector3D& pvertex,Float_t sigx,Float_t sigy);
  void PrimaryVertexWithWeigths_XY2(Int_t flag,TObjArray& tracks,Float_t zprim,Vector3D& pvertex,Float_t sigx,Float_t sigy);
  void CheckPrimaryVertexResolution(Int_t flag, UEventNode* outJ, UEventNode* outS);  
  
  void MethodComb(UEventNode* inNodeJ,UEventNode* inNodeS);

private:

  Bool_t fGoodVertex;

  Float_t fZprim;

  Int_t fJuraProdCounter;
  Int_t fSaleveProdCounter;

  Vector3D fPrimaryVertexPost;
  Vector3D fPrimaryVertexJ;
  Vector3D fPrimaryVertexS;

  Int_t fCombReco;

  Int_t fAllTracksJ;    
  Int_t fAllTracksS;    

  TH1F* fhVtxDx;
  TH1F* fhVtxDy;
  TH1F* fhVtxDz;

  TH1F* fhVtxDx_Glob;
  TH1F* fhVtxDy_Glob;
  TH1F* fhVtxDz_Glob;

  TH1F* fhVtxDx_GlobAcce;
  TH1F* fhVtxDy_GlobAcce;
  TH1F* fhVtxDz_GlobAcce;

  TH1F* fhVtxDx_Glob2;
  TH1F* fhVtxDy_Glob2;
  TH1F* fhVtxDz_Glob2;

  TH1F* fhVtxDx_Glob2Acce;
  TH1F* fhVtxDy_Glob2Acce;
  TH1F* fhVtxDz_Glob2Acce;

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

  TH1F* fhVtxX;
  TH1F* fhVtxY;
  TH1F* fhTracksVsVtxX_J;
  TH1F* fhTracksVsVtxY_J;
  TH1F* fhTracksVsVtxX_S;
  TH1F* fhTracksVsVtxY_S;

  TH1F* fhTrackDx;
  TH1F* fhTrackDy;
  TH1F* fhTrackDz;

public:
  
  //ClassDef(Na61PrimaryVertexRecoPostModule,1) // DOCUMENT ME
  
};

#endif
//____________________________________________________________________
//
