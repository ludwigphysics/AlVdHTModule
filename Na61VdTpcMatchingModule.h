// -*- mode: c++ -*-
//
// This module make hit map which is used to select only hits from the
// vicinity of the track extrapolation to the given VD station 
//
// Author: Paweł Staszel
//
//
#ifndef NA61_Na61VdTpcMatchingModule
#define NA61_Na61VdTpcMatchingModule

#include <evt/Event.h>
#include <evt/RecEvent.h>
#include <utl/ShineUnits.h>
#include <utl/PDGParticleIds.h>
#include <utl/DatabasePDG.h>
#include <io/EventFileChain.h>

#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include <map>
#include <set>
#include <utility>
#include <det/Silicon.h>

#include <evt/BOSBank.h>
#include <evt/BOSRecord.h>
#include <det/Detector.h>
#include <evt/raw/Silicon.h>
#include <det/Silicon.h>
#include <fwk/CentralConfig.h>
#include <det/MagneticFieldTracker.h>
#include <utl/UTCDateTime.h>
//#include <utl/WithUnit.h>
#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>
#include <list>
#include <sstream>

#include <modutils/KalmanFitterAK.h>
#include <modutils/KalmanFilterWB.h>


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
#ifndef UTIL_UEventNode
#include "UEventNode.h"
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

using namespace std;
//using namespace io;
using namespace utl;
using namespace evt;
using namespace evt::rec;
using namespace utl::ParticleConst;

using std::cout;
using std::endl;
using std::ios;
using std::setprecision;
using std::setw;
using std::ifstream;
using std::cin;

class Na61VdTpcMatchingModule : public Na61Module
{ 

public:

  Na61VdTpcMatchingModule();
  Na61VdTpcMatchingModule(const Char_t* name, const Char_t* title);
  virtual ~Na61VdTpcMatchingModule() { }
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(const det::MagneticFieldTracker& tra, evt::Event& event, UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option="B") const; // *MENU*


  void SetRunId(Int_t v){fRunId=v;}
  void SetTracker(const det::MagneticFieldTracker& /*v*/){}
  void SetNSig(double v) {fNSig = v;}

  double GetTpcTracks(){return fTpcTracks;}
  Int_t   GetNumberOfVdTpcTracks(){return fNumberOfVdTpcTracks;}

private:

  void CorrectForTpcOffset(RecEvent& recEvent, const Vertex& mainVertex);
  void TransformToNA61Frame_HT(UEventNode* node, TObjArray& vdtracks, UEventNode* inJ, UEventNode* inS);
  void TransformHits(UEventNode* node, double dx, double dy, double dz);
  void MatchTracks(const det::MagneticFieldTracker& tra, RecEvent& recEvent,const Vertex& mainVertex,TObjArray& vdtracks);
  UVdTrack* MakeTpcTrackAsVd(const det::MagneticFieldTracker& tra, const evt::rec::Track& track);

private:

	UInt_t fMatchArr[5000]; 
  //static det::MagneticFieldTracker& ftracker;

  Na61VdParameters* fVdParams;
  //Na61ArmParameters* fParams;
  double fNSig;
  
  double fMomDiffVtx;
  double fMomDiffVtx2;

  double fTpcTracks;
  Int_t fNumberOfVdTpcTracks;

  double fOffsetX;
  double fOffsetY;
  double fOffsetZ;
  
  Int_t fRunId;


  TH1F* hMomAllTpc;
  TH1F* hMomAllPos;
  TH1F* hMomAllNeg;
  TH1F* hPtAllTpc;
  

  TH2F* fhdEdx_all;
  TH2F* fhdEdxPos_all;
  TH2F* fhdEdxNeg_all;

  TH2F* fhdEdx_mtpc;
  TH2F* fhdEdxPos_mtpc;
  TH2F* fhdEdxNeg_mtpc;

  TH2F* fhdEdxVsClusters_all;
  TH2F* fhdEdxVsClustersCut_all;
  TH2F* fhdEdxVsClusters_mtpc;
  TH2F* fhdEdxVsClustersCut_mtpc;
  TH2F* fhClustersVsMom_all;

  TH2F* fhMomPoints_ZX;
  TH2F* fhMomPoints_ZX_fine;
 
  TH1F* fhTanY;
  TH1F* fhTanY_acce;
  TH2F* fhTanYvsY;
  TH2F* fhTanYvsY_acce;
  
  TH2F* fhY_JJ_vsY_cutx;
  TH2F* fhY_SS_vsY_cutx;
  TH2F* fhY_JS_vsY_cutx;
  TH2F* fhY_SJ_vsY_cutx;
  
  TH2F* fhY_JJ_vsY_cutx_r1;
  TH2F* fhY_SS_vsY_cutx_r1;
  TH2F* fhY_JS_vsY_cutx_r1;
  TH2F* fhY_SJ_vsY_cutx_r1;
  
  TH2F* fhY_JJ_vsY_cutx_r2;
  TH2F* fhY_SS_vsY_cutx_r2;
  TH2F* fhY_JS_vsY_cutx_r2;
  TH2F* fhY_SJ_vsY_cutx_r2;
  
  TH2F* fhX_JJ_vsZ_cuty;
  TH2F* fhX_SS_vsZ_cuty;
  TH2F* fhX_JS_vsZ_cuty;
  TH2F* fhX_SJ_vsZ_cuty;
  
  TH2F* fhXvsPz_J;
  TH2F* fhXvsPz_S;
  TH2F* fhXvsPx_J;
  TH2F* fhXvsPx_S;
  
  TH1F* fhPosMatch[3];
  TH1F* fhPosMatch_best[3];
  TH1F* fhMomMatch[3];
  TH1F* fhMomMatch_cutx[3];
  TH1F* fhMomMatch_cuty[3];
  TH1F* fhMomMatch_cutxy[3]; 
  TH1F* fhMomMatch_best[3];
  
  TH2F* fhMatch_TanYvsY_acce;
  TH2F* fhMatch_TanYvsY_best;

  TH2F* fhChargeCorr;
  TH2F* fhMomVsKfMom;
  TH2F* fhMomVsKfMom_kf;
  TH2F* fhMomVsKfMom_acce;
  TH1F* fhMomDiff;
  TH1F* fhMomDiff_kf;
  
  TH2F* fhTdistVsMomDiff;
  TH2F* fhXY_AtTarget;
  TH2F* fhXY_AtTarget_clean;

  TH1F* fhMainVertexVsVdVertex_X;
  TH1F* fhMainVertexVsVdVertex_Y;
  TH1F* fhMainVertexVsVdVertex_Z;

  TH2F* fhTracksVsTracksWithMatch;

  TH2F* fhXY_Vds4_All;
  TH2F* fhXY_Vds4_WithMatch;
  TH2F* fhXY_VTpc1_All;
  TH2F* fhXY_VTpc1_WithMatch;

  TH2F* fhXY_Vds4_NoMatch;
  TH2F* fhXY_fTPC_NoMatch;
  TH2F* fhXY_Vds4_Match;
  TH2F* fhXY_fTPC_Match;

  TH1F* fhCurvature_All;
  TH1F* fhCurvature_WithMatch;

  TH2F* fhCurvature_m_vs_m;
  TH2F* fhCurvature_m_vs_nm;

  TH2F* fhXOnVds4_m_vs_m;
  TH2F* fhXOnVds4_m_vs_nm;
  TH2F* fhYOnVds4_m_vs_m;
  TH2F* fhYOnVds4_m_vs_nm;

  TH1F* fhClusters_All;
  TH1F* fhClusters_All_Acce;
  TH1F* fhClusters_WithMatch;
  
  TH1F* fhPt;

public:
  
  //ClassDef(Na61VdTpcMatchingModule,0) // DOCUMENT ME
  
};

#endif
//____________________________________________________________________
//
