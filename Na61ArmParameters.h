//-*-mode:c++-*-

#ifndef Na61_Na61ArmParameters
#define Na61_Na61ArmParameters

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif

#ifndef UTIL_UVector3D
#include "UVector3D.h"
#endif

class Na61ArmParameters : public TObject {
 public:
  Na61ArmParameters(TString armname);
  virtual ~Na61ArmParameters();

  void Init();

  void SetRunId(int run_id) { fRunId = run_id; }
  void SetDrotX(int sensor_id, double rot) {
    fSensorId = sensor_id;
    fDrotX = rot;
  }
  void SetDrotY(int sensor_id, double rot) {
    fSensorId = sensor_id;
    fDrotY = rot;
  }
  void SetDrotZ(int sensor_id, double rot) {
    fSensorId = sensor_id;
    fDrotZ = rot;
  }
  void SetVolumeDz(int sensor_id, double dz) { fSensorId = sensor_id, fDz = dz; }

  void SetRotX(int sensor_id, double v) { fRotX[sensor_id] = v; }  // rotation angle around x axis
  void SetRotY(int sensor_id, double v) { fRotY[sensor_id] = v; }  // rotation angle around y axis
  void SetRotZ(int sensor_id, double v) { fRotZ[sensor_id] = v; }  // rotation angle around z axis

  void SetRotGlobX(int sensor_id, double v) { fRotGlobX[sensor_id] = v; }  // rotation angle around x axis
  void SetRotGlobY(int sensor_id, double v) { fRotGlobY[sensor_id] = v; }  // rotation angle around y axis
  void SetRotGlobZ(int sensor_id, double v) { fRotGlobZ[sensor_id] = v; }  // rotation angle around z axis

  // positions of sensors volumes
  void SetVolumeX(int sensor_id, double v) { fVolumeX[sensor_id] = v; }
  void SetVolumeY(int sensor_id, double v) { fVolumeY[sensor_id] = v; }
  void SetVolumeZ(int sensor_id, double v) { fVolumeZ[sensor_id] = v; }

  void SetVolumeGlobX(int sensor_id, double v) { fVolumeGlobX[sensor_id] = v; }
  void SetVolumeGlobY(int sensor_id, double v) { fVolumeGlobY[sensor_id] = v; }
  void SetVolumeGlobZ(int sensor_id, double v) { fVolumeGlobZ[sensor_id] = v; }

  double GetRotX(int sensor_id) { return fRotX[sensor_id]; }  // rotation angle around x axis
  double GetRotY(int sensor_id) { return fRotY[sensor_id]; }  // rotation angle around y axis
  double GetRotZ(int sensor_id) { return fRotZ[sensor_id]; }  // rotation angle around z axis

  double GetRotGlobX(int sensor_id) { return fRotGlobX[sensor_id]; }  // rotation angle around x axis
  double GetRotGlobY(int sensor_id) { return fRotGlobY[sensor_id]; }  // rotation angle around y axis
  double GetRotGlobZ(int sensor_id) { return fRotGlobZ[sensor_id]; }  // rotation angle around z axis

  // positions of sensors volumes
  double GetVolumeX(int sensor_id) { return fVolumeX[sensor_id]; }
  double GetVolumeY(int sensor_id) { return fVolumeY[sensor_id]; }
  double GetVolumeZ(int sensor_id) { return fVolumeZ[sensor_id]; }

  double GetVolumeGlobX(int sensor_id) { return fVolumeGlobX[sensor_id]; }
  double GetVolumeGlobY(int sensor_id) { return fVolumeGlobY[sensor_id]; }
  double GetVolumeGlobZ(int sensor_id) { return fVolumeGlobZ[sensor_id]; }

  double GetDevOffsetX(int idev) { return fOffsetx[idev]; }
  double GetDevOffsetY(int idev) { return fOffsety[idev]; }
  double GetDevSigmaX(int idev) { return fSigmax[idev]; }
  double GetDevSigmaY(int idev) { return fSigmay[idev]; }

  double GetCombineOffsetX(int idev) { return fOffx[idev]; }
  double GetCombineOffsetY(int idev) { return fOffy[idev]; }
  double GetCombineSigmaX(int idev) { return fSigx[idev]; }
  double GetCombineSigmaY(int idev) { return fSigy[idev]; }

  double GetDxOffset(int imatch) { return fOff_dx[imatch]; }
  double GetDyOffset(int imatch) { return fOff_dy[imatch]; }
  double GetDaxOffset(int imatch) { return fOff_dax[imatch]; }
  double GetDayOffset(int imatch) { return fOff_day[imatch]; }

  double GetDxSigma(int imatch) { return fSig_dx[imatch]; }
  double GetDySigma(int imatch) { return fSig_dy[imatch]; }
  double GetDaxSigma(int imatch) { return fSig_dax[imatch]; }
  double GetDaySigma(int imatch) { return fSig_day[imatch]; }

  double GetMeanLocX(int i, int j) { return fMeanLocX[i][j]; }  // "i" denotes GlobalCandidatesArray and "j" denotes sensor
  double GetMeanLocY(int i, int j) { return fMeanLocY[i][j]; }
  double GetSigmaLocX(int i, int j) { return fSigmaLocX[i][j]; }
  double GetSigmaLocY(int i, int j) { return fSigmaLocY[i][j]; }

  double GetMeanX(int i, int j) { return fMeanX[i][j]; }  // "i" denotes GlobalCandidatesArray and "j" denotes sensor
  double GetMeanY(int i, int j) { return fMeanY[i][j]; }
  double GetSigmaX(int i, int j) { return fSigmaX[i][j]; }
  double GetSigmaY(int i, int j) { return fSigmaY[i][j]; }

  double GetAxCut() { return fAxCut; }

  bool GetInit() { return fInit; }

  const char* GetArmName() { return fArmName.Data(); }

  void TransformVolumesToGlobal(double dx, double dy, double dz, double rotx, double roty, double rotz);
  void SetVolumesInGlobal(double dx, double dy, double dz, double rotx, double roty, double rotz);

 private:
  // private methods

  // geometry parameters and matching parameters
  /////// JURA
  void SetupSensorGeometry_JuraDec2016_mv1(int run_id);
  void SetupSensorRotation_JuraDec2016_mv1(int run_id);
  void SetupSensorGeometry_JuraDec2016(int run_id);
  void SetupSensorRotation_JuraDec2016(int run_id);
  void SetupSensorGeometry_JuraNov2018(int run_id);
  void SetupSensorRotation_JuraNov2018(int run_id);

  void SetupOffsetsAndSigmas_JuraDec2016_nofield(int run_id);
  void SetupOffsetsAndSigmas_JuraDec2016_field(int run_id);

  void SetupOffsetsAndSigmas_JuraOct2017_nofield(int run_id);
  void SetupOffsetsAndSigmas_JuraOct2017_field(int run_id);
  void SetupOffsetsAndSigmas_JuraOct2017_field_xela75(int run_id);
  void SetupOffsetsAndSigmas_JuraOct2017_field_xela40(int run_id);
  void SetupOffsetsAndSigmas_JuraNov2018_field(int run_id);
  

  void CellurarAutomatonParams_Jura(int run_id);
  void CellurarAutomatonParams_Jura_XeLa(int run_id);

  /////// SALEVE
  void SetupSensorGeometry_SaleveDec2016_mv1(int run_id);
  void SetupSensorRotation_SaleveDec2016_mv1(int run_id);
  void SetupSensorGeometry_SaleveDec2016(int run_id);
  void SetupSensorRotation_SaleveDec2016(int run_id);
  void SetupSensorGeometry_SaleveNov2018(int run_id);
  void SetupSensorRotation_SaleveNov2018(int run_id);

  void SetupOffsetsAndSigmas_SaleveDec2016_nofield(int run_id);
  void SetupOffsetsAndSigmas_SaleveDec2016_field(int run_id);

  void SetupOffsetsAndSigmas_SaleveOct2017_nofield(int run_id);
  void SetupOffsetsAndSigmas_SaleveOct2017_field(int run_id);
  void SetupOffsetsAndSigmas_SaleveOct2017_field_xela75(int run_id);
  void SetupOffsetsAndSigmas_SaleveOct2017_field_xela40(int run_id);
  void SetupOffsetsAndSigmas_SaleveNov2018_field(int run_id);
  
   // old paramater from July 2016 test
  void SetupSensorGeometry_July2016(int run_id);
  void SetupSensorRotation_July2016();
  void SetupOffsetsAndSigmas_July2016(int run_id);

  void CellurarAutomatonParams_Saleve(int run_id);
  void CellurarAutomatonParams_Saleve_XeLa(int run_id);

  void SetCAMembers(double* meanLocX_Vds1, double* sigmaLocX_Vds1, double* meanLocY_Vds1, double* sigmaLocY_Vds1, double* meanLocX_Vds2, double* sigmaLocX_Vds2, double* meanLocY_Vds2, double* sigmaLocY_Vds2, double* meanLocX_Vds3, double* sigmaLocX_Vds3, double* meanLocY_Vds3, double* sigmaLocY_Vds3, double* meanLocX_Vds4, double* sigmaLocX_Vds4, double* meanLocY_Vds4, double* sigmaLocY_Vds4, double* meanX_Vds1, double* sigmaX_Vds1, double* meanY_Vds1, double* sigmaY_Vds1, double* meanX_Vds2, double* sigmaX_Vds2, double* meanY_Vds2, double* sigmaY_Vds2, double* meanX_Vds3, double* sigmaX_Vds3, double* meanY_Vds3, double* sigmaY_Vds3, double* meanX_Vds4, double* sigmaX_Vds4, double* meanY_Vds4, double* sigmaY_Vds4);

 private:
  int fSensorId;
  TString fArmName;

  bool fJura;

  int fRunId;

  bool fInit;

  double fAxCut;  // cut used is VdTrackingModule to select target production

  // parameters used in CellurarAutomaton
  double fMeanLocX[4][8];
  double fMeanLocY[4][8];
  double fSigmaLocX[4][8];
  double fSigmaLocY[4][8];
  double fMeanX[4][8];
  double fMeanY[4][8];
  double fSigmaX[4][8];
  double fSigmaY[4][8];

  // variables used for geo tunning
  double fDrotX;
  double fDrotY;
  double fDrotZ;
  double fDz;

  double fRotX[8];
  double fRotY[8];
  double fRotZ[8];
  double fVolumeX[8];
  double fVolumeY[8];
  double fVolumeZ[8];

  double fRotGlobX[8];
  double fRotGlobY[8];
  double fRotGlobZ[8];
  double fVolumeGlobX[8];
  double fVolumeGlobY[8];
  double fVolumeGlobZ[8];

  // 7,8-th are extra 3sensor combination
  double fOffsetx[20];
  double fOffsety[20];
  double fSigmax[20];
  double fSigmay[20];

  double fOffx[18];
  double fOffy[18];
  double fSigx[18];
  double fSigy[18];

  double fOff_dx[4];
  double fOff_dy[4];
  double fOff_dax[4];
  double fOff_day[4];

  double fSig_dx[4];
  double fSig_dy[4];
  double fSig_dax[4];
  double fSig_day[4];

 public:
  // friend ostream& operator<< (ostream& os,Na61ArmParameters* pixel);

  //  ClassDef(Na61ArmParameters,2)  //  QDC data class
};

#endif
