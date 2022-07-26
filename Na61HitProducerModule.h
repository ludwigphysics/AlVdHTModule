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
#ifndef NA61_Na61HitProducerModule
#define NA61_Na61HitProducerModule

#ifndef NA61_Na61Module
#include "Na61Module.h"
#endif
#ifndef UTIL_UEventNode
#include "UEventNode.h"
#endif
#ifndef UTIL_USensorHit
#include "USensorHit.h"
#endif

#ifndef NA61_Na61ArmParameters
#include "Na61ArmParameters.h"
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

class Na61HitProducerModule : public Na61Module {
 public:
  Na61HitProducerModule();
  Na61HitProducerModule(const char* name, const char* title);
  virtual ~Na61HitProducerModule() {}
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
  void FillHSaxayWithPixels(UDataTable* pixs, TH2F* hh_pixs);
  void FillHSaxayWithClusters(UDataTable* hits, TH2F* hh_clust, TH2F* hh_clustpos, TH2F* hh_largeclustpos, int ii);
  double FillClusterDistr(UDataTable* hits, TH1F* h_multi, TH1F* h_multiLt1, TH1F* h_size);

  void FillHorizontalDistributions(UEventNode* outNode);

  void LocalToGlobal(int ii, UDataTable* tab);

 private:
  Na61ArmParameters* fParams;
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
  int faxayToInd[576 + 1][1152 + 1];
  static const int findx = 576;
  static const int findy = 1152;
  int fNRes;
  int fResIndx[2500];
  int fResIndy[2500];
  TObjArray fCluster;

  TH2F* fhHitsXVsZ;
  TH2F* fhHitsYVsZ;

  TH2F* fhPixelMap[8];
  TH2F* fhClusterMap[8];
  TH2F* fhClusterPosMap[8];
  TH2F* fhLargeClusterPosMap[8];
  TH1F* fhClusters[8];
  TH1F* fhClustersLt1[8];

  TH1F* fhClusterSize[8];

  int fClusters[8];      // all clusters
  int fRealClusters[8];  // size>1

  TH1F* fhClusterVerticalDist_Vds1;
  TH1F* fhClusterVerticalDist_Vds2;
  TH1F* fhClusterVerticalDist_Vds3;
  TH1F* fhClusterVerticalDist_Vds4a;
  TH1F* fhClusterVerticalDist_Vds4b;

  TH1F* fhSignalToNoise;

 public:
  //  ClassDef(Na61HitProducerModule,0) // DOCUMENT ME
};

#endif
//____________________________________________________________________
//
