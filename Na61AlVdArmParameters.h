//-*-mode:c++-*-

#ifndef Na61_Na61AlVdArmParameters
#define Na61_Na61AlVdArmParameters


#ifndef ROOT_TObject
#include "TObject.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif

#ifndef UTIL_UVector3D
#include "UVector3D.h"
#endif



class Na61AlVdArmParameters: public TObject 
{

public:
  Na61AlVdArmParameters(TString armname); 
  virtual ~Na61AlVdArmParameters();  

  void Init();

  void SetRunId(Int_t run_id){fRunId = run_id;}
  void SetDrotX(Int_t sensor_id, Float_t rot){fSensorId = sensor_id; fDrotX = rot;}
  void SetDrotY(Int_t sensor_id, Float_t rot){fSensorId = sensor_id; fDrotY = rot;}
  void SetDrotZ(Int_t sensor_id, Float_t rot){fSensorId = sensor_id; fDrotZ = rot;}
  void SetVolumeDz(Int_t sensor_id, Float_t dz){fSensorId = sensor_id, fDz = dz;}

  void SetRotX(Int_t sensor_id, Float_t v){fRotX[sensor_id] = v;}    // rotation angle around x axis
  void SetRotY(Int_t sensor_id, Float_t v){fRotY[sensor_id] = v;}    // rotation angle around y axis
  void SetRotZ(Int_t sensor_id, Float_t v){fRotZ[sensor_id] = v;}    // rotation angle around z axis

  void SetRotGlobX(Int_t sensor_id, Float_t v){fRotGlobX[sensor_id] = v;}    // rotation angle around x axis
  void SetRotGlobY(Int_t sensor_id, Float_t v){fRotGlobY[sensor_id] = v;}    // rotation angle around y axis
  void SetRotGlobZ(Int_t sensor_id, Float_t v){fRotGlobZ[sensor_id] = v;}    // rotation angle around z axis

  //positions of sensors volumes
  void SetVolumeX(Int_t sensor_id, Float_t v) {fVolumeX[sensor_id] = v;} 
  void SetVolumeY(Int_t sensor_id, Float_t v) {fVolumeY[sensor_id] = v;}
  void SetVolumeZ(Int_t sensor_id, Float_t v) {fVolumeZ[sensor_id] = v;}

  void SetVolumeGlobX(Int_t sensor_id, Float_t v) {fVolumeGlobX[sensor_id] = v;} 
  void SetVolumeGlobY(Int_t sensor_id, Float_t v) {fVolumeGlobY[sensor_id] = v;}
  void SetVolumeGlobZ(Int_t sensor_id, Float_t v) {fVolumeGlobZ[sensor_id] = v;}

  Float_t GetRotX(Int_t sensor_id){return fRotX[sensor_id];}    // rotation angle around x axis
  Float_t GetRotY(Int_t sensor_id){return fRotY[sensor_id];}    // rotation angle around y axis
  Float_t GetRotZ(Int_t sensor_id){return fRotZ[sensor_id];}    // rotation angle around z axis

  Float_t GetRotGlobX(Int_t sensor_id){return fRotGlobX[sensor_id];}    // rotation angle around x axis
  Float_t GetRotGlobY(Int_t sensor_id){return fRotGlobY[sensor_id];}    // rotation angle around y axis
  Float_t GetRotGlobZ(Int_t sensor_id){return fRotGlobZ[sensor_id];}    // rotation angle around z axis

  //positions of sensors volumes
  Float_t GetVolumeX(Int_t sensor_id) {return fVolumeX[sensor_id];} 
  Float_t GetVolumeY(Int_t sensor_id) {return fVolumeY[sensor_id];}
  Float_t GetVolumeZ(Int_t sensor_id) {return fVolumeZ[sensor_id];}

  Float_t GetVolumeGlobX(Int_t sensor_id) {return fVolumeGlobX[sensor_id];} 
  Float_t GetVolumeGlobY(Int_t sensor_id) {return fVolumeGlobY[sensor_id];}
  Float_t GetVolumeGlobZ(Int_t sensor_id) {return fVolumeGlobZ[sensor_id];}
  
  Float_t GetDevOffsetX(Int_t idev){return fOffsetx[idev];}
  Float_t GetDevOffsetY(Int_t idev){return fOffsety[idev];}
  Float_t GetDevSigmaX(Int_t idev){return fSigmax[idev];}
  Float_t GetDevSigmaY(Int_t idev){return fSigmay[idev];}
  
  Float_t GetCombineOffsetX(Int_t idev){return fOffx[idev];}
  Float_t GetCombineOffsetY(Int_t idev){return fOffy[idev];}
  Float_t GetCombineSigmaX(Int_t idev){return fSigx[idev];}
  Float_t GetCombineSigmaY(Int_t idev){return fSigy[idev];}
  
  Float_t GetDxOffset(Int_t imatch){return fOff_dx[imatch];}
  Float_t GetDyOffset(Int_t imatch){return fOff_dy[imatch];}
  Float_t GetDaxOffset(Int_t imatch){return fOff_dax[imatch];}
  Float_t GetDayOffset(Int_t imatch){return fOff_day[imatch];}
  
  Float_t GetDxSigma(Int_t imatch){return fSig_dx[imatch];}
  Float_t GetDySigma(Int_t imatch){return fSig_dy[imatch];}
  Float_t GetDaxSigma(Int_t imatch){return fSig_dax[imatch];}
  Float_t GetDaySigma(Int_t imatch){return fSig_day[imatch];}

  Float_t GetMeanLocX(Int_t i, Int_t j){return fMeanLocX[i][j];} // "i" denotes GlobalCandidatesArray and "j" denotes sensor 
  Float_t GetMeanLocY(Int_t i, Int_t j){return fMeanLocY[i][j];} 
  Float_t GetSigmaLocX(Int_t i, Int_t j){return fSigmaLocX[i][j];} 
  Float_t GetSigmaLocY(Int_t i, Int_t j){return fSigmaLocY[i][j];} 

  Float_t GetMeanX(Int_t i, Int_t j){return fMeanX[i][j];} // "i" denotes GlobalCandidatesArray and "j" denotes sensor 
  Float_t GetMeanY(Int_t i, Int_t j){return fMeanY[i][j];} 
  Float_t GetSigmaX(Int_t i, Int_t j){return fSigmaX[i][j];} 
  Float_t GetSigmaY(Int_t i, Int_t j){return fSigmaY[i][j];} 

  Float_t GetAxCut(){return fAxCut;}

  Bool_t GetInit(){return fInit;}  
  

  void TransformVolumesToGlobal(Float_t dx, Float_t dy, Float_t dz,Float_t rotx, Float_t roty, Float_t rotz);
  void SetVolumesInGlobal(Float_t dx, Float_t dy, Float_t dz,Float_t rotx, Float_t roty, Float_t rotz);

  void LocalToGlobal(Int_t ii, TObjArray* hits);

  
private:
  // private methods

  void SetupSensorGeometry_JuraNominal(Int_t run_id);
  void SetupSensorGeometry_SaleveNominal(Int_t run_id);

  void SetupSensorRotation_JuraNominal(Int_t run_id);
  void SetupSensorRotation_SaleveNominal(Int_t run_id);

  void SetupSensorGeometry_Jura_pPb2022(Int_t run_id);
  void SetupSensorGeometry_Saleve_pPb2022(Int_t run_id);

  void SetupSensorRotation_Jura_pPb2022(Int_t run_id);
  void SetupSensorRotation_Saleve_pPb2022(Int_t run_id);

  void SetupOffsetsAndSigmas_JuraNominal(Int_t run_id);
  void SetupOffsetsAndSigmas_SaleveNominal(Int_t run_id);

  void SetupOffsetsAndSigmas_Jura_pPb2022(Int_t run_id);
  void SetupOffsetsAndSigmas_Saleve_pPb2022(Int_t run_id);

  
private:

  Int_t fSensorId;
  TString fArmName;

  Bool_t fJura;

  Int_t fRunId;  

  Bool_t fInit;

  Float_t fAxCut; // cut used is VdTrackingModule to select target production

  Int_t fNDev;

  // parameters used in CellurarAutomaton
  Float_t fMeanLocX[4][8];
  Float_t fMeanLocY[4][8];
  Float_t fSigmaLocX[4][8];
  Float_t fSigmaLocY[4][8];
  Float_t fMeanX[4][8];
  Float_t fMeanY[4][8];
  Float_t fSigmaX[4][8];
  Float_t fSigmaY[4][8];

  // variables used for geo tunning
  Float_t fDrotX;
  Float_t fDrotY;
  Float_t fDrotZ;
  Float_t fDz;

  Float_t fRotX[34];
  Float_t fRotY[34];
  Float_t fRotZ[34];
  Float_t fVolumeX[34];
  Float_t fVolumeY[34];
  Float_t fVolumeZ[34];

  Float_t fRotGlobX[34];
  Float_t fRotGlobY[34];
  Float_t fRotGlobZ[34];
  Float_t fVolumeGlobX[34];
  Float_t fVolumeGlobY[34];
  Float_t fVolumeGlobZ[34];

  // 3sensor combination used in Make3HitTracks
  Float_t fOffsetx[150];
  Float_t fOffsety[150];
  Float_t fSigmax[150];
  Float_t fSigmay[150];

  Float_t fOffx[18];
  Float_t fOffy[18];
  Float_t fSigx[18];
  Float_t fSigy[18];

  Float_t fOff_dx[4];
  Float_t fOff_dy[4];
  Float_t fOff_dax[4];
  Float_t fOff_day[4];

  Float_t fSig_dx[4];
  Float_t fSig_dy[4];
  Float_t fSig_dax[4];
  Float_t fSig_day[4];


public:
  //friend ostream& operator<< (ostream& os,Na61AlVdArmParameters* pixel);
  
  //  ClassDef(Na61AlVdArmParameters,0)  //  QDC data class
};

#endif

