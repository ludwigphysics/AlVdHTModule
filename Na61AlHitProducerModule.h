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
#ifndef NA61_Na61AlHitProducerModule
#define NA61_Na61AlHitProducerModule

#ifndef NA61_Na61Module
#include "Na61Module.h"
#endif
#ifndef UTIL_UEventNode
#include "UEventNode.h"
#endif
#ifndef UTIL_USensorHit
#include "USensorHit.h"
#endif

#ifndef NA61_Na61AlVdArmParameters
#include "Na61AlVdArmParameters.h"
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

class Na61AlHitProducerModule : public Na61Module {
 public:
  Na61AlHitProducerModule();
  Na61AlHitProducerModule(const char* name, const char* title);
  virtual ~Na61AlHitProducerModule() {}
  virtual void DefineHistograms();
  virtual void Init();
  virtual void Begin();
  virtual void Event(UEventNode* inNode, UEventNode* outNode);
  virtual void End();
  virtual void Finish();
  virtual void Print(Option_t* option = "B") const;  // *MtENU*

  void SetProductionMode(int v) { fProductionMode = v; }
  void SetRunId(int v) { fRunId = v; }
  void SetToJuraArm(bool v) { fJura = v; }

  void SetDrotX(int sensor_id, double rotang) {
    fSensorId = sensor_id;
    fDrotX = rotang;
  }
  void SetDrotY(int sensor_id, double rotang) {
    fSensorId = sensor_id;
    fDrotY = rotang;
  }
  void SetDrotZ(int sensor_id, double rotang) {
    fSensorId = sensor_id;
    fDrotZ = rotang;
  }
  void SetVolumeDz(int sensor_id, double dz) {
    fSensorId = sensor_id;
    fDz = dz;
  }

  void SetRunInfo(TFile* file);  // call after init of IOMModule

  int GetNumberOfClusters(int i) { return fClusters[i]; }          // number of clusters in i-th sensor
  int GetNumberOfRealClusters(int i) { return fRealClusters[i]; }  // number of real clusters (size>1)in i-th sensor

 private:
  void MakeClusters(TString sensorname, UDataTable* pixels, UDataTable* hits);
  void CheckRing2(int i, int j, UDataTable* pixels, TObjArray* cluster);
  void SetaxayToIndArray(UDataTable* pixels);
  void FillHSaxayWithPixels(UDataTable* pixs, TH2F* hh_pix,TH2F* hh_pixpos);
  void FillHSaxayWithClusters(UDataTable* hits, TH2F* hh_clust, TH2F* hh_clustpos, int ii);
  double FillClusterDistr(UDataTable* hits, TH1F* h_multi, TH1F* h_multiLt1, TH1F* h_size);

  void FillHorizontalDistributions(UEventNode* outNode);

  void LocalToGlobal(int ii, UDataTable* tab);

 private:
  Na61AlVdArmParameters* fParams;
  double fDrotX;
  double fDrotY;
  double fDrotZ;
  double fDz;
  int fSensorId;

  bool fJura;

  int fNorm;
  int fRunId;
  int fProductionMode;
  TObjArray* fPixelInfoArray;
  TObjArray* fClusterArray;
  int faxayToInd[512 + 1][1024 + 1];
  static const int findx = 512;
  static const int findy = 1024;
  int fNRes;
  int fResIndx[2500];
  int fResIndy[2500];
  TObjArray fCluster;

  TH2F* fhHitsXVsZ;
  TH2F* fhHitsYVsZ;

  TH2F* fhPixelMap[34];
  TH2F* fhPixelPosMap[34];
  TH2F* fhClusterMap[34];
  TH2F* fhClusterPosMap[34];
  TH2F* fhLargeClusterPosMap[34];
  TH1F* fhClusters[34];
  TH1F* fhClustersLt1[34];

  TH1F* fhClusterSize[34];

  int fClusters[34];      // all clusters
  int fRealClusters[34];  // size>1

  TH1F* fhClusterVerticalDist_Vds1;
  TH1F* fhClusterVerticalDist_Vds2;
  TH1F* fhClusterVerticalDist_Vds3;
  TH1F* fhClusterVerticalDist_Vds4;
  TH1F* fhClusterHorizontalDist_Vds1;
  TH1F* fhClusterHorizontalDist_Vds2;
  TH1F* fhClusterHorizontalDist_Vds3;
  TH1F* fhClusterHorizontalDist_Vds4;


  TH1F* fhSignalToNoise;

 public:
  //  ClassDef(Na61HitProducerModule,0) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
