//-*-mode:c++-*-

#ifndef Na61_Na61VdParameters
#define Na61_Na61VdParameters

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TString
#include "TString.h"
#endif

#ifndef ROOT_TFile
#include "TFile.h"
#endif

#ifndef ROOT_TF1
#include "TF1.h"
#endif

#ifndef UTIL_UVector3D
#include "UVector3D.h"
#endif

class Na61VdParameters : public TObject {
 public:
  Na61VdParameters(TString armname);
  virtual ~Na61VdParameters();

  void Init();

  // void FindShineRun(int runid,int& runid_shine,int& chunks);
  //  void FindShineRun(int runid, int& calrunid, int& runid_shine,int& chunks);

  void SetRunId(int run_id) { fRunId = run_id; }
  void SetMatchParamsFile(TFile* file) { fMatchParamsFile = file; }
  void SetN(int v) { fN = v; }
  void PrintSetupInfo();

  double GetOffsetToTpcX() { return fOffsetToTpcX; }
  double GetOffsetToTpcY() { return fOffsetToTpcY; }
  double GetOffsetToTpcZ() { return fOffsetToTpcZ; }

  double GetJuraArmOffAx() { return fOffAx_J; }
  double GetJuraArmSigAx() { return fSigAx_J; }
  double GetJuraArmOffAy() { return fOffAy_J; }
  double GetJuraArmSigAy() { return fSigAy_J; }

  double GetSaleveArmOffAx() { return fOffAx_S; }
  double GetSaleveArmSigAx() { return fSigAx_S; }
  double GetSaleveArmOffAy() { return fOffAy_S; }
  double GetSaleveArmSigAy() { return fSigAy_S; }

  double GetJuraArmRotX() { return fRotX_J; }
  double GetJuraArmRotY() { return fRotY_J; }
  double GetJuraArmRotZ() { return fRotZ_J; }
  double GetSaleveArmRotX() { return fRotX_S; }
  double GetSaleveArmRotY() { return fRotY_S; }
  double GetSaleveArmRotZ() { return fRotZ_S; }
  Vector3D GetJuraArmOffset() { return fJuraArmOffset; }
  Vector3D GetSaleveArmOffset() { return fSaleveArmOffset; }

  double GetBeamSpotOffsetX() { return fBeamSpotOffsetX; }
  double GetBeamSpotSigmaX() { return fBeamSpotSigmaX; }
  double GetBeamSpotOffsetY() { return fBeamSpotOffsetY; }
  double GetBeamSpotSigmaY() { return fBeamSpotSigmaY; }

  double GetVtxDzSigma() { return fVtxDzSigma; }

  double GetDevOffx(int i) { return fDevOffx[i]; }
  double GetDevSigx(int i) { return fDevSigx[i]; }
  double GetDevOffy(int i) { return fDevOffy[i]; }
  double GetDevSigy(int i) { return fDevSigy[i]; }

  double GetResOffx(int i) { return fResOffx[i]; }
  double GetResSigx(int i) { return fResSigx[i]; }
  double GetResOffy(int i) { return fResOffy[i]; }
  double GetResSigy(int i) { return fResSigy[i]; }

  int GetDxMatchParams(double z, double x, double xvd, double& mean, double& sigma);
  int GetDyMatchParams(double z, double y, double x, double xvd, double& mean, double& sigma);
  bool GetInit() { return fInit; }

  void SetdOffZ(double v) { fdOffZ = v; }
  static unsigned int FindCalibRuns(const unsigned int run);

  // private:
  // private methods

  void MakeCompansationForRotation();

  void SetupCommonVdSystem();
  void SetupCommonVdSystem_pPb();
  void SetupCommonVdSystem_xela150();
  void SetupCommonVdSystem_xela75();
  void SetupCommonVdSystem_xela40();
  void SetupCommonVdSystem_pPb2018();
  void SetupCommonVdSystem_pPb2022();
  void SetupCommonVdSystem_pbpb();

  
  void ReadMatchParams();
  // void SetupFunctions(TString histName,TF1* mean, TF1* sigma);
  void SetupFunctions(TString histName, int i);

  void SetupCommonNA61System();

 private:
  bool fCommonVdSystemSetup;
  bool fInit;
  TFile* fMatchParamsFile;
  int fRunId;

  double fdOffZ;

  double fOffsetToTpcX;
  double fOffsetToTpcY;
  double fOffsetToTpcZ;

  double fDevOffx[4];
  double fDevSigx[4];
  double fDevOffy[4];
  double fDevSigy[4];

  double fResOffx[4];
  double fResSigx[4];
  double fResOffy[4];
  double fResSigy[4];

  int fN;

  double fBeamSpotOffsetX;
  double fBeamSpotSigmaX;
  double fBeamSpotOffsetY;
  double fBeamSpotSigmaY;

  double fVtxDzSigma;

  double fOffAx_J;
  double fSigAx_J;
  double fOffAy_J;
  double fSigAy_J;
  double fOffAx_S;
  double fSigAx_S;
  double fOffAy_S;
  double fSigAy_S;

  // parmeters to define common system
  double fRotX_J;
  double fRotY_J;
  double fRotZ_J;
  double fRotX_S;
  double fRotY_S;
  double fRotZ_S;

  Vector3D fJuraArmOffset;
  Vector3D fSaleveArmOffset;

  TF1* fMean_dx[8];
  TF1* fSigma_dx[8];

  TF1* fMean_dy[8];
  TF1* fSigma_dy[8];

  /*
    TF1* fMean_dx_JJ_1;
    TF1* fMean_dx_JS_1;
    TF1* fMean_dx_SS_1;
    TF1* fMean_dx_SJ_1;
    TF1* fSigma_dx_JJ_1;
    TF1* fSigma_dx_JS_1;
    TF1* fSigma_dx_SS_1;
    TF1* fSigma_dx_SJ_1;

    TF1* fMean_dx_JJ_2;
    TF1* fMean_dx_JS_2;
    TF1* fMean_dx_SS_2;
    TF1* fMean_dx_SJ_2;
    TF1* fSigma_dx_JJ_2;
    TF1* fSigma_dx_JS_2;
    TF1* fSigma_dx_SS_2;
    TF1* fSigma_dx_SJ_2;

    TF1* fMean_dy_JJ;
    TF1* fMean_dy_JS;
    TF1* fMean_dy_SS;
    TF1* fMean_dy_SJ;
    TF1* fSigma_dy_JJ;
    TF1* fSigma_dy_JS;
    TF1* fSigma_dy_SS;
    TF1* fSigma_dy_SJ;
  */

 public:
  friend ostream& operator<<(ostream& os, Na61VdParameters* pixel);

  //  ClassDef(Na61VdParameters,1)  //  QDC data class
};

#endif
