//____________________________________________________________________
//
// Interfase to set and acces arm geometry and matching
// parameters

#include "Na61AlVdArmParameters.h"
#include "USensorHit.h"

#ifndef ROOT_TMath
#include "TMath.h"
#endif
#ifndef ROOT_TObjArray
#include "TObjArray.h"
#endif

#ifndef __IOSTREAM__
#include <iostream>
#endif
using std::endl;
using std::cout;

//____________________________________________________________________
//ClassImp(Na61AlVdArmParameters);

//____________________________________________________________________
Na61AlVdArmParameters::Na61AlVdArmParameters(TString armname) 
{
  // normal constructor 
  fInit = kFALSE;
  fArmName = armname;
  fJura = kFALSE;
  if(armname.Contains("Jura"))fJura = kTRUE;

  fNDev = 150;
}

//____________________________________________________________________
Na61AlVdArmParameters::~Na61AlVdArmParameters() 
{
}

//____________________________________________________________________
void Na61AlVdArmParameters::Init() 
{

  SetDrotX(0,0); 
  SetDrotY(0,0); 
  SetDrotZ(0,0); 
  SetVolumeDz(0,0);

  fInit=kTRUE;

  if(fRunId<1191){

    if(fJura) SetupSensorRotation_JuraNominal(fRunId);  // used for physics reco of Dec2016 data
    else SetupSensorRotation_SaleveNominal(fRunId);  // used for physics reco of Dec2016 data
    
    if(fJura) SetupSensorGeometry_JuraNominal(fRunId);
    else SetupSensorGeometry_SaleveNominal(fRunId);

  }else if(fRunId<10000){ /////// pPb2022 

    if(fJura) SetupSensorRotation_Jura_pPb2022(fRunId);  // used for physics reco of Dec2016 data
    else SetupSensorRotation_Saleve_pPb2022(fRunId);  // used for physics reco of Dec2016 data
    
    if(fJura) SetupSensorGeometry_Jura_pPb2022(fRunId);
    else SetupSensorGeometry_Saleve_pPb2022(fRunId);
    
  } else {
    if(fJura) SetupSensorRotation_JuraNominal(fRunId);  // used for physics reco of Dec2016 data
    else SetupSensorRotation_SaleveNominal(fRunId);  // used for physics reco of Dec2016 data
    
    if(fJura) SetupSensorGeometry_JuraNominal(fRunId);
    else SetupSensorGeometry_SaleveNominal(fRunId);

    cout<<"no sensor geometry defined for this run number: fRunId="<<fRunId<<" thus the nominal geometry was assumed"<<endl;
  }


  ///////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// offsets for 3hit tracks ////////////////////////////////////
  if(fRunId<1191){

    if(fJura)SetupOffsetsAndSigmas_JuraNominal(fRunId);
    else SetupOffsetsAndSigmas_SaleveNominal(fRunId);

  }else if(fRunId<10000){ /////// pPb2022 

    if(fJura)SetupOffsetsAndSigmas_Jura_pPb2022(fRunId);
    else SetupOffsetsAndSigmas_Saleve_pPb2022(fRunId);
    
  } else {
    if(fJura)SetupOffsetsAndSigmas_JuraNominal(fRunId);
    else SetupOffsetsAndSigmas_SaleveNominal(fRunId);
    cout<<"no sensor offsets defined for this run number: fRunId="<<fRunId<<" thus the nominal offsets were assumed"<<endl;
  }  

}

//_____________________________________________________________
void Na61AlVdArmParameters::TransformVolumesToGlobal(Float_t dx, Float_t dy, Float_t dz, Float_t rotx, Float_t roty, Float_t rotz)
{
  cout<<" Na61AlVdArmParameters::TransformVolumesToGlobal: Transforming Volumes, arm offx="<<dx<<"  rotx="<<rotx<<" roty="<<roty<<" rotz="<<rotz<<endl;
  
  Double_t sa = TMath::Sin(rotz); 
  Double_t ca = TMath::Cos(rotz);
  Double_t sb = TMath::Sin(roty); 
  Double_t cb = TMath::Cos(roty);
  Double_t sg = TMath::Sin(rotx);
  Double_t cg = TMath::Cos(rotx);
  
  
  for(Int_t i=0;i<8;i++){
    
    Float_t x0 = GetVolumeX(i);
    Float_t y0 = GetVolumeY(i);
    Float_t z0 = GetVolumeZ(i);
    
    z0 = z0 - 75.0;
    
    //Float_t x1 =      cb * x0                +     sb * z0;
    //Float_t y1 =  -sg*sb * x0   +   cg * y0  +  sg*cb * z0;
    //Float_t z1 =  -cg*sb * x0   -   sg * y0  +  cg*cb * z0;
    
    //ZXY
    Float_t x1 =  (ca*cb-sa*sg*sb) * x0   +   sa*cg * y0  +    (ca*sb+sa*sg*cb) * z0;
    Float_t y1 = -(sa*cb+ca*sg*sb) * x0   +   ca*cg * y0  +   (-sa*sb+ca*sg*cb) * z0;
    Float_t z1 =            -cg*sb * x0   -      sg * y0  +              cg*cb  * z0;


    z1 = z1 + 75.0;
    
    //cout<<"z0="<<z0<<" z1="<<z1<<endl;
    
    SetVolumeX(i,x1+dx);
    SetVolumeY(i,y1+dy);
    SetVolumeZ(i,z1+dz);
    SetRotX(i,GetRotX(i)+rotx);
    SetRotY(i,GetRotY(i)+roty);
    
  }

  cout<<"after transformation volumeX_Al1_0="<<fVolumeX[0]<<"   volumeX_Al3_0="<<fVolumeX[2]<<"  volumeX_Al3_1="<<fVolumeX[3]<<endl; 


}

//_____________________________________________________________
void Na61AlVdArmParameters::SetVolumesInGlobal(Float_t dx, Float_t dy, Float_t dz, Float_t rotx, Float_t roty, Float_t rotz)
{
  cout<<" Na61AlVdArmParameters::SetVolumesInGlobal: Setting Volumes, arm offx="<<dx<<"  rotx="<<rotx<<" roty="<<roty<<" rotz="<<rotz<<endl;
  
  Double_t sa = TMath::Sin(rotz); 
  Double_t ca = TMath::Cos(rotz);
  Double_t sb = TMath::Sin(roty); 
  Double_t cb = TMath::Cos(roty);
  Double_t sg = TMath::Sin(rotx);
  Double_t cg = TMath::Cos(rotx);
  
  
  for(Int_t i=0;i<20;i++){
    
    Float_t x0 = GetVolumeX(i);
    Float_t y0 = GetVolumeY(i);
    Float_t z0 = GetVolumeZ(i);
    
    z0 = z0 - 75.0;
    
    //ZXY
    Float_t x1 =  (ca*cb-sa*sg*sb) * x0   +   sa*cg * y0  +    (ca*sb+sa*sg*cb) * z0;
    Float_t y1 = -(sa*cb+ca*sg*sb) * x0   +   ca*cg * y0  +   (-sa*sb+ca*sg*cb) * z0;
    Float_t z1 =            -cg*sb * x0   -      sg * y0  +              cg*cb  * z0;


    z1 = z1 + 75.0;
    
    //cout<<"z0="<<z0<<" z1="<<z1<<endl;
    
    SetVolumeGlobX(i,x1+dx);
    SetVolumeGlobY(i,y1+dy);
    SetVolumeGlobZ(i,z1+dz);
    //SetVolumeGlobX(i,x0+dx);
    //SetVolumeGlobY(i,y0+dy);
    //SetVolumeGlobZ(i,z0+dz);
    SetRotGlobX(i,GetRotX(i)+rotx);
    SetRotGlobY(i,GetRotY(i)+roty);
    SetRotGlobZ(i,GetRotZ(i)+rotz);
    
  }

  cout<<"Volumes in global frame: volumeX_Al1_0="<<fVolumeGlobX[0]<<"   volumeX_Al3_0="<<fVolumeGlobX[2]<<"  volumeX_Al3_1="<<fVolumeGlobX[3]<<endl; 

}


//____________________________________________________________________________________________
void Na61AlVdArmParameters::SetupSensorRotation_JuraNominal(Int_t run_id)
{
  cout<<"Na61AlVdArmParameters::SetupSensorRotation_JuraNominal: Setting sensor rotation for AlVd for run: "<<run_id<<endl;


  //X rotation
  double rotX_Al1_0 = 0;
  double rotX_Al1_1 = 0;
  double rotX_Al1_2 = 0;

  double rotX_Al2_0 = 0;
  double rotX_Al2_1 = 0;
  double rotX_Al2_2 = 0;
  double rotX_Al2_3 = 0;
  double rotX_Al2_4 = 0;
  double rotX_Al2_5 = 0;

  double rotX_Al3_0 = 0;
  double rotX_Al3_1 = 0;
  double rotX_Al3_2 = 0;
  double rotX_Al3_3 = 0;
  double rotX_Al3_4 = 0;
  double rotX_Al3_5 = 0;
  double rotX_Al3_6 = 0;
  double rotX_Al3_7 = 0;
  double rotX_Al3_8 = 0;
  double rotX_Al3_9 = 0;

  double rotX_Al4_0 = 0; 
  double rotX_Al4_1 = 0;
  double rotX_Al4_2 = 0;
  double rotX_Al4_3 = 0;
  double rotX_Al4_4 = 0;
  double rotX_Al4_5 = 0; 
  double rotX_Al4_6 = 0;
  double rotX_Al4_7 = 0;
  double rotX_Al4_8 = 0;
  double rotX_Al4_9 = 0;
  double rotX_Al4_10 = 0; 
  double rotX_Al4_11 = 0;
  double rotX_Al4_12 = 0;
  double rotX_Al4_13 = 0;
  double rotX_Al4_14 = 0;

  // Y rotation
  double rotY_Al1_0 = 0;
  double rotY_Al1_1 = 0;
  double rotY_Al1_2 = 0;

  double rotY_Al2_0 = 0;
  double rotY_Al2_1 = 0;
  double rotY_Al2_2 = 0;
  double rotY_Al2_3 = 0.1047;
  double rotY_Al2_4 = 0.1047;
  double rotY_Al2_5 = 0.1047;

  double rotY_Al3_0 = 0;
  double rotY_Al3_1 = 0;
  double rotY_Al3_2 = 0;
  double rotY_Al3_3 = 0;
  double rotY_Al3_4 = 0;
  double rotY_Al3_5 = 0.1047;
  double rotY_Al3_6 = 0.1047;
  double rotY_Al3_7 = 0.1047;
  double rotY_Al3_8 = 0.1047;
  double rotY_Al3_9 = 0.1047;


  double rotY_Al4_0 = 0;
  double rotY_Al4_1 = 0;
  double rotY_Al4_2 = 0;
  double rotY_Al4_3 = 0;
  double rotY_Al4_4 = 0;
  double rotY_Al4_5 = 0.1047;
  double rotY_Al4_6 = 0.1047;
  double rotY_Al4_7 = 0.1047;
  double rotY_Al4_8 = 0.1047;
  double rotY_Al4_9 = 0.1047;
  double rotY_Al4_10 = 0.1396;
  double rotY_Al4_11 = 0.1396;
  double rotY_Al4_12 = 0.1396;
  double rotY_Al4_13 = 0.1396;
  double rotY_Al4_14 = 0.1396;

  // Z rotation
  double rotZ_Al1_0 = 0;
  double rotZ_Al1_1 = 0;
  double rotZ_Al1_2 = 0;

  double rotZ_Al2_0 = 0;
  double rotZ_Al2_1 = 0;
  double rotZ_Al2_2 = 0;
  double rotZ_Al2_3 = 0;
  double rotZ_Al2_4 = 0;
  double rotZ_Al2_5 = 0;

  double rotZ_Al3_0 = 0;
  double rotZ_Al3_1 = 0;
  double rotZ_Al3_2 = 0;
  double rotZ_Al3_3 = 0;
  double rotZ_Al3_4 = 0;
  double rotZ_Al3_5 = 0;
  double rotZ_Al3_6 = 0;
  double rotZ_Al3_7 = 0;
  double rotZ_Al3_8 = 0;
  double rotZ_Al3_9 = 0;

  double rotZ_Al4_0 = 0;
  double rotZ_Al4_1 = 0;
  double rotZ_Al4_2 = 0;
  double rotZ_Al4_3 = 0;
  double rotZ_Al4_4 = 0;
  double rotZ_Al4_5 = 0;
  double rotZ_Al4_6 = 0;
  double rotZ_Al4_7 = 0;
  double rotZ_Al4_8 = 0;
  double rotZ_Al4_9 = 0;
  double rotZ_Al4_10 = 0;
  double rotZ_Al4_11 = 0;
  double rotZ_Al4_12 = 0;
  double rotZ_Al4_13 = 0;
  double rotZ_Al4_14 = 0;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  fRotX[0] = rotX_Al1_0;
  fRotX[1] = rotX_Al1_1;
  fRotX[2] = rotX_Al1_2;

  fRotX[3] = rotX_Al2_0;
  fRotX[4] = rotX_Al2_1; 
  fRotX[5] = rotX_Al2_2;
  fRotX[6] = rotX_Al2_3;
  fRotX[7] = rotX_Al2_4;
  fRotX[8] = rotX_Al2_5;

  fRotX[9]  = rotX_Al3_0;
  fRotX[10] = rotX_Al3_1;
  fRotX[11] = rotX_Al3_2;
  fRotX[12] = rotX_Al3_3;
  fRotX[13] = rotX_Al3_4;
  fRotX[14] = rotX_Al3_5;
  fRotX[15] = rotX_Al3_6;
  fRotX[16] = rotX_Al3_7;
  fRotX[17] = rotX_Al3_8;
  fRotX[18] = rotX_Al3_9;

  fRotX[19] = rotX_Al4_0;
  fRotX[20] = rotX_Al4_1;
  fRotX[21] = rotX_Al4_2;
  fRotX[22] = rotX_Al4_3;
  fRotX[23] = rotX_Al4_4;
  fRotX[24] = rotX_Al4_5;
  fRotX[25] = rotX_Al4_6;
  fRotX[26] = rotX_Al4_7;
  fRotX[27] = rotX_Al4_8;
  fRotX[28] = rotX_Al4_9;
  fRotX[29] = rotX_Al4_10;
  fRotX[30] = rotX_Al4_11;
  fRotX[31] = rotX_Al4_12;
  fRotX[32] = rotX_Al4_13;
  fRotX[33] = rotX_Al4_14;

  ///////////////////////////////////////

  fRotY[0] = rotY_Al1_0;
  fRotY[1] = rotY_Al1_1;
  fRotY[2] = rotY_Al1_2;

  fRotY[3] = rotY_Al2_0;
  fRotY[4] = rotY_Al2_1; 
  fRotY[5] = rotY_Al2_2;
  fRotY[6] = rotY_Al2_3;
  fRotY[7] = rotY_Al2_4;
  fRotY[8] = rotY_Al2_5;

  fRotY[9]  = rotY_Al3_0;
  fRotY[10] = rotY_Al3_1;
  fRotY[11] = rotY_Al3_2;
  fRotY[12] = rotY_Al3_3;
  fRotY[13] = rotY_Al3_4;
  fRotY[14] = rotY_Al3_5;
  fRotY[15] = rotY_Al3_6;
  fRotY[16] = rotY_Al3_7;
  fRotY[17] = rotY_Al3_8;
  fRotY[18] = rotY_Al3_9;

  fRotY[19] = rotY_Al4_0;
  fRotY[20] = rotY_Al4_1;
  fRotY[21] = rotY_Al4_2;
  fRotY[22] = rotY_Al4_3;
  fRotY[23] = rotY_Al4_4;
  fRotY[24] = rotY_Al4_5;
  fRotY[25] = rotY_Al4_6;
  fRotY[26] = rotY_Al4_7;
  fRotY[27] = rotY_Al4_8;
  fRotY[28] = rotY_Al4_9;
  fRotY[29] = rotY_Al4_10;
  fRotY[30] = rotY_Al4_11;
  fRotY[31] = rotY_Al4_12;
  fRotY[32] = rotY_Al4_13;
  fRotY[33] = rotY_Al4_14;

  ///////////////////////////////////////////

  fRotZ[0] = rotZ_Al1_0;
  fRotZ[1] = rotZ_Al1_1;
  fRotZ[2] = rotZ_Al1_2;

  fRotZ[3] = rotZ_Al2_0;
  fRotZ[4] = rotZ_Al2_1; 
  fRotZ[5] = rotZ_Al2_2;
  fRotZ[6] = rotZ_Al2_3;
  fRotZ[7] = rotZ_Al2_4;
  fRotZ[8] = rotZ_Al2_5;

  fRotZ[9]  = rotZ_Al3_0;
  fRotZ[10] = rotZ_Al3_1;
  fRotZ[11] = rotZ_Al3_2;
  fRotZ[12] = rotZ_Al3_3;
  fRotZ[13] = rotZ_Al3_4;
  fRotZ[14] = rotZ_Al3_5;
  fRotZ[15] = rotZ_Al3_6;
  fRotZ[16] = rotZ_Al3_7;
  fRotZ[17] = rotZ_Al3_8;
  fRotZ[18] = rotZ_Al3_9;

  fRotZ[19] = rotZ_Al4_0;
  fRotZ[20] = rotZ_Al4_1;
  fRotZ[21] = rotZ_Al4_2;
  fRotZ[22] = rotZ_Al4_3;
  fRotZ[23] = rotZ_Al4_4;
  fRotZ[24] = rotZ_Al4_5;
  fRotZ[25] = rotZ_Al4_6;
  fRotZ[26] = rotZ_Al4_7;
  fRotZ[27] = rotZ_Al4_8;
  fRotZ[28] = rotZ_Al4_9;
  fRotZ[29] = rotZ_Al4_10;
  fRotZ[30] = rotZ_Al4_11;
  fRotZ[31] = rotZ_Al4_12;
  fRotZ[32] = rotZ_Al4_13;
  fRotZ[33] = rotZ_Al4_14;



  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;
  
  cout<<"Atd: fDrotX = "<< fDrotX<<"   fDrotY = "<< fDrotY<<"    fDrotZ = "<<fDrotZ<<endl;


}
//____________________________________________________________________________________________
void Na61AlVdArmParameters::SetupSensorRotation_Jura_pPb2022(Int_t run_id)
{
  cout<<"Na61AlVdArmParameters::SetupSensorRotation_Jura_pPb2022: Setting sensor rotation for AlVd run: "<<run_id<<endl;

double rotX_Al1_0 = 0;
double rotX_Al1_1 = 0;
double rotX_Al1_2 = 0;
double rotX_Al2_0 = -0.00214203;
double rotX_Al2_1 = 4.0665e-05;
double rotX_Al2_2 = 0.00290969;
double rotX_Al2_3 = -0.00124912;
double rotX_Al2_4 = -0.000160175;
double rotX_Al2_5 = -0.00200341;
double rotX_Al3_0 = -0.00142411;
double rotX_Al3_1 = -0.00431986;
double rotX_Al3_2 = -0.00174657;
double rotX_Al3_3 = 0.00199561;
double rotX_Al3_4 = -2.73469e-05;
double rotX_Al3_5 = -0.00450169;
double rotX_Al3_6 = -0.00152413;
double rotX_Al3_7 = -0.00111412;
double rotX_Al3_8 = -0.00329143;
double rotX_Al3_9 = -0.00251563;
double rotX_Al4_0 = -0.00809654;
double rotX_Al4_1 = -0.0099946;
double rotX_Al4_2 = 0.00168576;
double rotX_Al4_3 = 0.0222536;
double rotX_Al4_4 = 0.000939226;
double rotX_Al4_5 = -0.000262824;
double rotX_Al4_6 = 0.00249684;
double rotX_Al4_7 = -0.000347481;
double rotX_Al4_8 = 0.00199988;
double rotX_Al4_9 = -0.0102828;
double rotX_Al4_10 = -0.00370239;
double rotX_Al4_11 = 0.00147103;
double rotX_Al4_12 = 0.00250779;
double rotX_Al4_13 = 0.00319486;
double rotX_Al4_14 = 0.00148538;


///////////////////////////////////////////
double rotY_Al1_0 = 0;
double rotY_Al1_1 = 0;
double rotY_Al1_2 = 0;
double rotY_Al2_0 = -0.00548024;
double rotY_Al2_1 = -0.00414973;
double rotY_Al2_2 = -0.000800975;
double rotY_Al2_3 = 0.110442;
double rotY_Al2_4 = 0.106138;
double rotY_Al2_5 = 0.103616;
double rotY_Al3_0 = 0.00630733;
double rotY_Al3_1 = 0.00368437;
double rotY_Al3_2 = 0.00160947;
double rotY_Al3_3 = 0.00805897;
double rotY_Al3_4 = 0.00149505;
double rotY_Al3_5 = 0.106872;
double rotY_Al3_6 = 0.111179;
double rotY_Al3_7 = 0.110705;
double rotY_Al3_8 = 0.106026;
double rotY_Al3_9 = 0.114431;
double rotY_Al4_0 = 0.0286247;
double rotY_Al4_1 = 0.00709539;
double rotY_Al4_2 = 0.000947338;
double rotY_Al4_3 = 0.016937;
double rotY_Al4_4 = -0.00244546;
double rotY_Al4_5 = 0.103365;
double rotY_Al4_6 = 0.121426;
double rotY_Al4_7 = 0.113687;
double rotY_Al4_8 = 0.118807;
double rotY_Al4_9 = 0.112718;
double rotY_Al4_10 = 0.141422;
double rotY_Al4_11 = 0.145772;
double rotY_Al4_12 = 0.143593;
double rotY_Al4_13 = 0.138806;
double rotY_Al4_14 = 0.167619;


///////////////////////////////////////////////////////////
double rotZ_Al1_0 = 0;
double rotZ_Al1_1 = 0;
double rotZ_Al1_2 = 0;
double rotZ_Al2_0 = 0.000750624;
double rotZ_Al2_1 = 6.07035e-05;
double rotZ_Al2_2 = -0.000104254;
double rotZ_Al2_3 = 0.00055113;
double rotZ_Al2_4 = 0.000151553;
double rotZ_Al2_5 = -0.000586879;
double rotZ_Al3_0 = -0.000214632;
double rotZ_Al3_1 = -0.00106656;
double rotZ_Al3_2 = -0.000508824;
double rotZ_Al3_3 = -0.000717756;
double rotZ_Al3_4 = -0.00227059;
double rotZ_Al3_5 = 0.00101441;
double rotZ_Al3_6 = 0.000310183;
double rotZ_Al3_7 = -3.2358e-05;
double rotZ_Al3_8 = -0.000554956;
double rotZ_Al3_9 = 0.000236719;
double rotZ_Al4_0 = -0.00265538;
double rotZ_Al4_1 = -0.000157674;
double rotZ_Al4_2 = 0.000517842;
double rotZ_Al4_3 = 0.000244059;
double rotZ_Al4_4 = -0.00178256;
double rotZ_Al4_5 = 0.00341459;
double rotZ_Al4_6 = -0.000674913;
double rotZ_Al4_7 = 0.000314333;
double rotZ_Al4_8 = 0.000158358;
double rotZ_Al4_9 = -0.0016615;
double rotZ_Al4_10 = 0.000327343;
double rotZ_Al4_11 = 0.000954545;
double rotZ_Al4_12 = 0.000184516;
double rotZ_Al4_13 = -0.00142507;
double rotZ_Al4_14 = 0.00391568;


////////////////////////////////////////////////////////////////////////////////////////////////

  fRotX[0] = rotX_Al1_0;
  fRotX[1] = rotX_Al1_1;
  fRotX[2] = rotX_Al1_2;

  fRotX[3] = rotX_Al2_0;
  fRotX[4] = rotX_Al2_1; 
  fRotX[5] = rotX_Al2_2;
  fRotX[6] = rotX_Al2_3;
  fRotX[7] = rotX_Al2_4;
  fRotX[8] = rotX_Al2_5;

  fRotX[9]  = rotX_Al3_0;
  fRotX[10] = rotX_Al3_1;
  fRotX[11] = rotX_Al3_2;
  fRotX[12] = rotX_Al3_3;
  fRotX[13] = rotX_Al3_4;
  fRotX[14] = rotX_Al3_5;
  fRotX[15] = rotX_Al3_6;
  fRotX[16] = rotX_Al3_7;
  fRotX[17] = rotX_Al3_8;
  fRotX[18] = rotX_Al3_9;

  fRotX[19] = rotX_Al4_0;
  fRotX[20] = rotX_Al4_1;
  fRotX[21] = rotX_Al4_2;
  fRotX[22] = rotX_Al4_3;
  fRotX[23] = rotX_Al4_4;
  fRotX[24] = rotX_Al4_5;
  fRotX[25] = rotX_Al4_6;
  fRotX[26] = rotX_Al4_7;
  fRotX[27] = rotX_Al4_8;
  fRotX[28] = rotX_Al4_9;
  fRotX[29] = rotX_Al4_10;
  fRotX[30] = rotX_Al4_11;
  fRotX[31] = rotX_Al4_12;
  fRotX[32] = rotX_Al4_13;
  fRotX[33] = rotX_Al4_14;

  ///////////////////////////////////////

  fRotY[0] = rotY_Al1_0;
  fRotY[1] = rotY_Al1_1;
  fRotY[2] = rotY_Al1_2;

  fRotY[3] = rotY_Al2_0;
  fRotY[4] = rotY_Al2_1; 
  fRotY[5] = rotY_Al2_2;
  fRotY[6] = rotY_Al2_3;
  fRotY[7] = rotY_Al2_4;
  fRotY[8] = rotY_Al2_5;

  fRotY[9]  = rotY_Al3_0;
  fRotY[10] = rotY_Al3_1;
  fRotY[11] = rotY_Al3_2;
  fRotY[12] = rotY_Al3_3;
  fRotY[13] = rotY_Al3_4;
  fRotY[14] = rotY_Al3_5;
  fRotY[15] = rotY_Al3_6;
  fRotY[16] = rotY_Al3_7;
  fRotY[17] = rotY_Al3_8;
  fRotY[18] = rotY_Al3_9;

  fRotY[19] = rotY_Al4_0;
  fRotY[20] = rotY_Al4_1;
  fRotY[21] = rotY_Al4_2;
  fRotY[22] = rotY_Al4_3;
  fRotY[23] = rotY_Al4_4;
  fRotY[24] = rotY_Al4_5;
  fRotY[25] = rotY_Al4_6;
  fRotY[26] = rotY_Al4_7;
  fRotY[27] = rotY_Al4_8;
  fRotY[28] = rotY_Al4_9;
  fRotY[29] = rotY_Al4_10;
  fRotY[30] = rotY_Al4_11;
  fRotY[31] = rotY_Al4_12;
  fRotY[32] = rotY_Al4_13;
  fRotY[33] = rotY_Al4_14;

  ///////////////////////////////////////////

  fRotZ[0] = rotZ_Al1_0;
  fRotZ[1] = rotZ_Al1_1;
  fRotZ[2] = rotZ_Al1_2;

  fRotZ[3] = rotZ_Al2_0;
  fRotZ[4] = rotZ_Al2_1; 
  fRotZ[5] = rotZ_Al2_2;
  fRotZ[6] = rotZ_Al2_3;
  fRotZ[7] = rotZ_Al2_4;
  fRotZ[8] = rotZ_Al2_5;

  fRotZ[9]  = rotZ_Al3_0;
  fRotZ[10] = rotZ_Al3_1;
  fRotZ[11] = rotZ_Al3_2;
  fRotZ[12] = rotZ_Al3_3;
  fRotZ[13] = rotZ_Al3_4;
  fRotZ[14] = rotZ_Al3_5;
  fRotZ[15] = rotZ_Al3_6;
  fRotZ[16] = rotZ_Al3_7;
  fRotZ[17] = rotZ_Al3_8;
  fRotZ[18] = rotZ_Al3_9;

  fRotZ[19] = rotZ_Al4_0;
  fRotZ[20] = rotZ_Al4_1;
  fRotZ[21] = rotZ_Al4_2;
  fRotZ[22] = rotZ_Al4_3;
  fRotZ[23] = rotZ_Al4_4;
  fRotZ[24] = rotZ_Al4_5;
  fRotZ[25] = rotZ_Al4_6;
  fRotZ[26] = rotZ_Al4_7;
  fRotZ[27] = rotZ_Al4_8;
  fRotZ[28] = rotZ_Al4_9;
  fRotZ[29] = rotZ_Al4_10;
  fRotZ[30] = rotZ_Al4_11;
  fRotZ[31] = rotZ_Al4_12;
  fRotZ[32] = rotZ_Al4_13;
  fRotZ[33] = rotZ_Al4_14;



  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;
  
  cout<<"Atd: fDrotX = "<< fDrotX<<"   fDrotY = "<< fDrotY<<"    fDrotZ = "<<fDrotZ<<endl;


}

//____________________________________________________________________________________________
void Na61AlVdArmParameters::SetupSensorRotation_SaleveNominal(Int_t run_id)
{
  cout<<"Na61AlVdArmParameters::SetupSensorRotation_SaleveNominal: Setting sensor rotation for AlVd for run: "<<run_id<<endl;


  //X rotation
  double rotX_Al1_0 = 0;
  double rotX_Al1_1 = 0;
  double rotX_Al1_2 = 0;

  double rotX_Al2_0 = 0;
  double rotX_Al2_1 = 0;
  double rotX_Al2_2 = 0;
  double rotX_Al2_3 = 0;
  double rotX_Al2_4 = 0;
  double rotX_Al2_5 = 0;

  double rotX_Al3_0 = 0;
  double rotX_Al3_1 = 0;
  double rotX_Al3_2 = 0;
  double rotX_Al3_3 = 0;
  double rotX_Al3_4 = 0;
  double rotX_Al3_5 = 0;
  double rotX_Al3_6 = 0;
  double rotX_Al3_7 = 0;
  double rotX_Al3_8 = 0;
  double rotX_Al3_9 = 0;

  double rotX_Al4_0 = 0; 
  double rotX_Al4_1 = 0;
  double rotX_Al4_2 = 0;
  double rotX_Al4_3 = 0;
  double rotX_Al4_4 = 0;
  double rotX_Al4_5 = 0; 
  double rotX_Al4_6 = 0;
  double rotX_Al4_7 = 0;
  double rotX_Al4_8 = 0;
  double rotX_Al4_9 = 0;
  double rotX_Al4_10 = 0; 
  double rotX_Al4_11 = 0;
  double rotX_Al4_12 = 0;
  double rotX_Al4_13 = 0;
  double rotX_Al4_14 = 0;

  // Y rotation
  double rotY_Al1_0 = 0;
  double rotY_Al1_1 = 0;
  double rotY_Al1_2 = 0;

  double rotY_Al2_0 = 0;
  double rotY_Al2_1 = 0;
  double rotY_Al2_2 = 0;
  double rotY_Al2_3 = -0.1047;
  double rotY_Al2_4 = -0.1047;
  double rotY_Al2_5 = -0.1047;

  double rotY_Al3_0 = 0;
  double rotY_Al3_1 = 0;
  double rotY_Al3_2 = 0;
  double rotY_Al3_3 = 0;
  double rotY_Al3_4 = 0;
  double rotY_Al3_5 = -0.1047;
  double rotY_Al3_6 = -0.1047;
  double rotY_Al3_7 = -0.1047;
  double rotY_Al3_8 = -0.1047;
  double rotY_Al3_9 = -0.1047;


  double rotY_Al4_0 = 0;
  double rotY_Al4_1 = 0;
  double rotY_Al4_2 = 0;
  double rotY_Al4_3 = 0;
  double rotY_Al4_4 = 0;
  double rotY_Al4_5 = -0.1047;
  double rotY_Al4_6 = -0.1047;
  double rotY_Al4_7 = -0.1047;
  double rotY_Al4_8 = -0.1047;
  double rotY_Al4_9 = -0.1047;
  double rotY_Al4_10 = -0.1396;
  double rotY_Al4_11 = -0.1396;
  double rotY_Al4_12 = -0.1396;
  double rotY_Al4_13 = -0.1396;
  double rotY_Al4_14 = -0.1396;

  // Z rotation
  double rotZ_Al1_0 = 0;
  double rotZ_Al1_1 = 0;
  double rotZ_Al1_2 = 0;

  double rotZ_Al2_0 = 0;
  double rotZ_Al2_1 = 0;
  double rotZ_Al2_2 = 0;
  double rotZ_Al2_3 = 0;
  double rotZ_Al2_4 = 0;
  double rotZ_Al2_5 = 0;

  double rotZ_Al3_0 = 0;
  double rotZ_Al3_1 = 0;
  double rotZ_Al3_2 = 0;
  double rotZ_Al3_3 = 0;
  double rotZ_Al3_4 = 0;
  double rotZ_Al3_5 = 0;
  double rotZ_Al3_6 = 0;
  double rotZ_Al3_7 = 0;
  double rotZ_Al3_8 = 0;
  double rotZ_Al3_9 = 0;

  double rotZ_Al4_0 = 0;
  double rotZ_Al4_1 = 0;
  double rotZ_Al4_2 = 0;
  double rotZ_Al4_3 = 0;
  double rotZ_Al4_4 = 0;
  double rotZ_Al4_5 = 0;
  double rotZ_Al4_6 = 0;
  double rotZ_Al4_7 = 0;
  double rotZ_Al4_8 = 0;
  double rotZ_Al4_9 = 0;
  double rotZ_Al4_10 = 0;
  double rotZ_Al4_11 = 0;
  double rotZ_Al4_12 = 0;
  double rotZ_Al4_13 = 0;
  double rotZ_Al4_14 = 0;

  fRotX[0] = rotX_Al1_0;
  fRotX[1] = rotX_Al1_1;
  fRotX[2] = rotX_Al1_2;

  fRotX[3] = rotX_Al2_0;
  fRotX[4] = rotX_Al2_1; 
  fRotX[5] = rotX_Al2_2;
  fRotX[6] = rotX_Al2_3;
  fRotX[7] = rotX_Al2_4;
  fRotX[8] = rotX_Al2_5;

  fRotX[9]  = rotX_Al3_0;
  fRotX[10] = rotX_Al3_1;
  fRotX[11] = rotX_Al3_2;
  fRotX[12] = rotX_Al3_3;
  fRotX[13] = rotX_Al3_4;
  fRotX[14] = rotX_Al3_5;
  fRotX[15] = rotX_Al3_6;
  fRotX[16] = rotX_Al3_7;
  fRotX[17] = rotX_Al3_8;
  fRotX[18] = rotX_Al3_9;

  fRotX[19] = rotX_Al4_0;
  fRotX[20] = rotX_Al4_1;
  fRotX[21] = rotX_Al4_2;
  fRotX[22] = rotX_Al4_3;
  fRotX[23] = rotX_Al4_4;
  fRotX[24] = rotX_Al4_5;
  fRotX[25] = rotX_Al4_6;
  fRotX[26] = rotX_Al4_7;
  fRotX[27] = rotX_Al4_8;
  fRotX[28] = rotX_Al4_9;
  fRotX[29] = rotX_Al4_10;
  fRotX[30] = rotX_Al4_11;
  fRotX[31] = rotX_Al4_12;
  fRotX[32] = rotX_Al4_13;
  fRotX[33] = rotX_Al4_14;

  ///////////////////////////////////////

  fRotY[0] = rotY_Al1_0;
  fRotY[1] = rotY_Al1_1;
  fRotY[2] = rotY_Al1_2;

  fRotY[3] = rotY_Al2_0;
  fRotY[4] = rotY_Al2_1; 
  fRotY[5] = rotY_Al2_2;
  fRotY[6] = rotY_Al2_3;
  fRotY[7] = rotY_Al2_4;
  fRotY[8] = rotY_Al2_5;

  fRotY[9]  = rotY_Al3_0;
  fRotY[10] = rotY_Al3_1;
  fRotY[11] = rotY_Al3_2;
  fRotY[12] = rotY_Al3_3;
  fRotY[13] = rotY_Al3_4;
  fRotY[14] = rotY_Al3_5;
  fRotY[15] = rotY_Al3_6;
  fRotY[16] = rotY_Al3_7;
  fRotY[17] = rotY_Al3_8;
  fRotY[18] = rotY_Al3_9;

  fRotY[19] = rotY_Al4_0;
  fRotY[20] = rotY_Al4_1;
  fRotY[21] = rotY_Al4_2;
  fRotY[22] = rotY_Al4_3;
  fRotY[23] = rotY_Al4_4;
  fRotY[24] = rotY_Al4_5;
  fRotY[25] = rotY_Al4_6;
  fRotY[26] = rotY_Al4_7;
  fRotY[27] = rotY_Al4_8;
  fRotY[28] = rotY_Al4_9;
  fRotY[29] = rotY_Al4_10;
  fRotY[30] = rotY_Al4_11;
  fRotY[31] = rotY_Al4_12;
  fRotY[32] = rotY_Al4_13;
  fRotY[33] = rotY_Al4_14;

  ///////////////////////////////////////////

  fRotZ[0] = rotZ_Al1_0;
  fRotZ[1] = rotZ_Al1_1;
  fRotZ[2] = rotZ_Al1_2;

  fRotZ[3] = rotZ_Al2_0;
  fRotZ[4] = rotZ_Al2_1; 
  fRotZ[5] = rotZ_Al2_2;
  fRotZ[6] = rotZ_Al2_3;
  fRotZ[7] = rotZ_Al2_4;
  fRotZ[8] = rotZ_Al2_5;

  fRotZ[9]  = rotZ_Al3_0;
  fRotZ[10] = rotZ_Al3_1;
  fRotZ[11] = rotZ_Al3_2;
  fRotZ[12] = rotZ_Al3_3;
  fRotZ[13] = rotZ_Al3_4;
  fRotZ[14] = rotZ_Al3_5;
  fRotZ[15] = rotZ_Al3_6;
  fRotZ[16] = rotZ_Al3_7;
  fRotZ[17] = rotZ_Al3_8;
  fRotZ[18] = rotZ_Al3_9;

  fRotZ[19] = rotZ_Al4_0;
  fRotZ[20] = rotZ_Al4_1;
  fRotZ[21] = rotZ_Al4_2;
  fRotZ[22] = rotZ_Al4_3;
  fRotZ[23] = rotZ_Al4_4;
  fRotZ[24] = rotZ_Al4_5;
  fRotZ[25] = rotZ_Al4_6;
  fRotZ[26] = rotZ_Al4_7;
  fRotZ[27] = rotZ_Al4_8;
  fRotZ[28] = rotZ_Al4_9;
  fRotZ[29] = rotZ_Al4_10;
  fRotZ[30] = rotZ_Al4_11;
  fRotZ[31] = rotZ_Al4_12;
  fRotZ[32] = rotZ_Al4_13;
  fRotZ[33] = rotZ_Al4_14;



  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;
  
  cout<<"AlVd: fDrotX = "<< fDrotX<<"   fDrotY = "<< fDrotY<<"    fDrotZ = "<<fDrotZ<<endl;


}

//____________________________________________________________________________________________
void Na61AlVdArmParameters::SetupSensorRotation_Saleve_pPb2022(Int_t run_id)
{
  cout<<"Na61AlVdArmParameters::SetupSensorRotation_Saleve_pPb2022: Setting sensor rotation for AlVd for run: "<<run_id<<endl;


  //X rotation
double rotX_Al1_0 = 0;
double rotX_Al1_1 = 0;
double rotX_Al1_2 = 0;
double rotX_Al2_0 = -0.00180032;
double rotX_Al2_1 = 0.000100265;
double rotX_Al2_2 = -4.23804e-06;
double rotX_Al2_3 = -0.0133984;
double rotX_Al2_4 = -0.00141736;
double rotX_Al2_5 = 3.65723e-06;
double rotX_Al3_0 = 0;
double rotX_Al3_1 = -0.000119797;
double rotX_Al3_2 = -0.000201875;
double rotX_Al3_3 = -0.0032288;
double rotX_Al3_4 = -0.0032;
double rotX_Al3_5 = 0.000177;
double rotX_Al3_6 = 0.000177412;
double rotX_Al3_7 = -0.000491506;
double rotX_Al3_8 = 0.000882609;
double rotX_Al3_9 = 0.00089;
double rotX_Al4_0 = 0.0100599;
double rotX_Al4_1 = -0.00192405;
double rotX_Al4_2 = 0.000118691;
double rotX_Al4_3 = 0.00175432;
double rotX_Al4_4 = 0.0222436;
double rotX_Al4_5 = 0.000284655;
double rotX_Al4_6 = 0.000339353;
double rotX_Al4_7 = 0.00190646;
double rotX_Al4_8 = 0.00282779;
double rotX_Al4_9 = 0.0161384;
double rotX_Al4_10 = -0.00166019;
double rotX_Al4_11 = 0.00428771;
double rotX_Al4_12 = -2.80949e-05;
double rotX_Al4_13 = 0.00202596;
double rotX_Al4_14 = 0.0097448;


  // Y rotation
double rotY_Al1_0 = 0;
double rotY_Al1_1 = 0;
double rotY_Al1_2 = 0;
double rotY_Al2_0 = -0.00110435;
double rotY_Al2_1 = -0.00335627;
double rotY_Al2_2 = -0.00377789;
double rotY_Al2_3 = -0.105646;
double rotY_Al2_4 = -0.107099;
double rotY_Al2_5 = -0.1149;
double rotY_Al3_0 = 0.00109271;
double rotY_Al3_1 = 0.00109271;
double rotY_Al3_2 = 0.00246709;
double rotY_Al3_3 = 0.00131136;
double rotY_Al3_4 = 0.001;
double rotY_Al3_5 = -0.1047;
double rotY_Al3_6 = -0.103912;
double rotY_Al3_7 = -0.100242;
double rotY_Al3_8 = -0.103738;
double rotY_Al3_9 = -0.1047;
double rotY_Al4_0 = -0.00183913;
double rotY_Al4_1 = 0.00660215;
double rotY_Al4_2 = 0.00143964;
double rotY_Al4_3 = -3.52733e-05;
double rotY_Al4_4 = 0.000116461;
double rotY_Al4_5 = -0.10525;
double rotY_Al4_6 = -0.0961704;
double rotY_Al4_7 = -0.0988247;
double rotY_Al4_8 = -0.0964409;
double rotY_Al4_9 = -0.0857792;
double rotY_Al4_10 = -0.159838;
double rotY_Al4_11 = -0.132392;
double rotY_Al4_12 = -0.125027;
double rotY_Al4_13 = -0.125834;
double rotY_Al4_14 = -0.135841;

  // Z rotation
double rotZ_Al1_0 = 0;
double rotZ_Al1_1 = 0;
double rotZ_Al1_2 = 0;
double rotZ_Al2_0 = -0.000202337;
double rotZ_Al2_1 = 0.000367774;
double rotZ_Al2_2 = -0.00022168;
double rotZ_Al2_3 = -0.000907295;
double rotZ_Al2_4 = 0.000261368;
double rotZ_Al2_5 = 0.000388155;
double rotZ_Al3_0 = -0.00024;
double rotZ_Al3_1 = -0.00024817;
double rotZ_Al3_2 = -0.000449916;
double rotZ_Al3_3 = -0.00117867;
double rotZ_Al3_4 = -0.0013;
double rotZ_Al3_5 = 0;
double rotZ_Al3_6 = -1.56415e-05;
double rotZ_Al3_7 = -0.000141338;
double rotZ_Al3_8 = -0.000929353;
double rotZ_Al3_9 = -0.00093;
double rotZ_Al4_0 = 0.00404804;
double rotZ_Al4_1 = -0.000214447;
double rotZ_Al4_2 = 2.09493e-05;
double rotZ_Al4_3 = -0.00111186;
double rotZ_Al4_4 = 0.0027438;
double rotZ_Al4_5 = -0.00178866;
double rotZ_Al4_6 = -0.000699514;
double rotZ_Al4_7 = -0.000249681;
double rotZ_Al4_8 = -0.000484944;
double rotZ_Al4_9 = 0.00330408;
double rotZ_Al4_10 = 0.000337487;
double rotZ_Al4_11 = 0.000725047;
double rotZ_Al4_12 = 4.06263e-05;
double rotZ_Al4_13 = 3.46138e-05;
double rotZ_Al4_14 = -0.00479974;

  fRotX[0] = rotX_Al1_0;
  fRotX[1] = rotX_Al1_1;
  fRotX[2] = rotX_Al1_2;

  fRotX[3] = rotX_Al2_0;
  fRotX[4] = rotX_Al2_1; 
  fRotX[5] = rotX_Al2_2;
  fRotX[6] = rotX_Al2_3;
  fRotX[7] = rotX_Al2_4;
  fRotX[8] = rotX_Al2_5;

  fRotX[9]  = rotX_Al3_0;
  fRotX[10] = rotX_Al3_1;
  fRotX[11] = rotX_Al3_2;
  fRotX[12] = rotX_Al3_3;
  fRotX[13] = rotX_Al3_4;
  fRotX[14] = rotX_Al3_5;
  fRotX[15] = rotX_Al3_6;
  fRotX[16] = rotX_Al3_7;
  fRotX[17] = rotX_Al3_8;
  fRotX[18] = rotX_Al3_9;

  fRotX[19] = rotX_Al4_0;
  fRotX[20] = rotX_Al4_1;
  fRotX[21] = rotX_Al4_2;
  fRotX[22] = rotX_Al4_3;
  fRotX[23] = rotX_Al4_4;
  fRotX[24] = rotX_Al4_5;
  fRotX[25] = rotX_Al4_6;
  fRotX[26] = rotX_Al4_7;
  fRotX[27] = rotX_Al4_8;
  fRotX[28] = rotX_Al4_9;
  fRotX[29] = rotX_Al4_10;
  fRotX[30] = rotX_Al4_11;
  fRotX[31] = rotX_Al4_12;
  fRotX[32] = rotX_Al4_13;
  fRotX[33] = rotX_Al4_14;

  ///////////////////////////////////////

  fRotY[0] = rotY_Al1_0;
  fRotY[1] = rotY_Al1_1;
  fRotY[2] = rotY_Al1_2;

  fRotY[3] = rotY_Al2_0;
  fRotY[4] = rotY_Al2_1; 
  fRotY[5] = rotY_Al2_2;
  fRotY[6] = rotY_Al2_3;
  fRotY[7] = rotY_Al2_4;
  fRotY[8] = rotY_Al2_5;

  fRotY[9]  = rotY_Al3_0;
  fRotY[10] = rotY_Al3_1;
  fRotY[11] = rotY_Al3_2;
  fRotY[12] = rotY_Al3_3;
  fRotY[13] = rotY_Al3_4;
  fRotY[14] = rotY_Al3_5;
  fRotY[15] = rotY_Al3_6;
  fRotY[16] = rotY_Al3_7;
  fRotY[17] = rotY_Al3_8;
  fRotY[18] = rotY_Al3_9;

  fRotY[19] = rotY_Al4_0;
  fRotY[20] = rotY_Al4_1;
  fRotY[21] = rotY_Al4_2;
  fRotY[22] = rotY_Al4_3;
  fRotY[23] = rotY_Al4_4;
  fRotY[24] = rotY_Al4_5;
  fRotY[25] = rotY_Al4_6;
  fRotY[26] = rotY_Al4_7;
  fRotY[27] = rotY_Al4_8;
  fRotY[28] = rotY_Al4_9;
  fRotY[29] = rotY_Al4_10;
  fRotY[30] = rotY_Al4_11;
  fRotY[31] = rotY_Al4_12;
  fRotY[32] = rotY_Al4_13;
  fRotY[33] = rotY_Al4_14;

  ///////////////////////////////////////////

  fRotZ[0] = rotZ_Al1_0;
  fRotZ[1] = rotZ_Al1_1;
  fRotZ[2] = rotZ_Al1_2;

  fRotZ[3] = rotZ_Al2_0;
  fRotZ[4] = rotZ_Al2_1; 
  fRotZ[5] = rotZ_Al2_2;
  fRotZ[6] = rotZ_Al2_3;
  fRotZ[7] = rotZ_Al2_4;
  fRotZ[8] = rotZ_Al2_5;

  fRotZ[9]  = rotZ_Al3_0;
  fRotZ[10] = rotZ_Al3_1;
  fRotZ[11] = rotZ_Al3_2;
  fRotZ[12] = rotZ_Al3_3;
  fRotZ[13] = rotZ_Al3_4;
  fRotZ[14] = rotZ_Al3_5;
  fRotZ[15] = rotZ_Al3_6;
  fRotZ[16] = rotZ_Al3_7;
  fRotZ[17] = rotZ_Al3_8;
  fRotZ[18] = rotZ_Al3_9;

  fRotZ[19] = rotZ_Al4_0;
  fRotZ[20] = rotZ_Al4_1;
  fRotZ[21] = rotZ_Al4_2;
  fRotZ[22] = rotZ_Al4_3;
  fRotZ[23] = rotZ_Al4_4;
  fRotZ[24] = rotZ_Al4_5;
  fRotZ[25] = rotZ_Al4_6;
  fRotZ[26] = rotZ_Al4_7;
  fRotZ[27] = rotZ_Al4_8;
  fRotZ[28] = rotZ_Al4_9;
  fRotZ[29] = rotZ_Al4_10;
  fRotZ[30] = rotZ_Al4_11;
  fRotZ[31] = rotZ_Al4_12;
  fRotZ[32] = rotZ_Al4_13;
  fRotZ[33] = rotZ_Al4_14;

  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;
  
  cout<<"AlVd: fDrotX = "<< fDrotX<<"   fDrotY = "<< fDrotY<<"    fDrotZ = "<<fDrotZ<<endl;

}


//____________________________________________________________________________________________
void Na61AlVdArmParameters::SetupSensorGeometry_JuraNominal(Int_t run_id)
{
  cout<<"Na61AlVdArmParameters::SetupSensorGeometry_JuraNominal: Setting sensor geometry for AlVd detector for run: "<<run_id<<endl;

  ////////////// X coordinate

  double VolumeX_Al1_0 = 0;
  double VolumeX_Al1_1 = 0;
  double VolumeX_Al1_2 = 0;

  double VolumeX_Al2_0 = 0;
  double VolumeX_Al2_1 = 0;
  double VolumeX_Al2_2 = 0;
  double VolumeX_Al2_3 = 12.56;
  double VolumeX_Al2_4 = 12.56;
  double VolumeX_Al2_5 = 12.56;

  double VolumeX_Al3_0 = 0 + 0.625;
  double VolumeX_Al3_1 = 0 + 0.625;
  double VolumeX_Al3_2 = 0 + 0.625;
  double VolumeX_Al3_3 = 0 + 0.625;
  double VolumeX_Al3_4 = 0 + 0.625;
  double VolumeX_Al3_5 = 12.56 + 0.625;
  double VolumeX_Al3_6 = 12.56 + 0.625;
  double VolumeX_Al3_7 = 12.56 + 0.625;
  double VolumeX_Al3_8 = 12.56 + 0.625;
  double VolumeX_Al3_9 = 12.56 + 0.625;

  double VolumeX_Al4_0 = 1.5 + 0.5818;
  double VolumeX_Al4_1 = 1.5 + 0.5818;
  double VolumeX_Al4_2 = 1.5 + 0.5818;
  double VolumeX_Al4_3 = 1.5 + 0.5818;
  double VolumeX_Al4_4 = 1.5 + 0.5818;
  double VolumeX_Al4_5 = 14.06 + 0.5818;
  double VolumeX_Al4_6 = 14.06 + 0.5818;
  double VolumeX_Al4_7 = 14.06 + 0.5818;
  double VolumeX_Al4_8 = 14.06 + 0.5818;
  double VolumeX_Al4_9 = 14.06 + 0.5818;
  double VolumeX_Al4_10 = 26.45 + 0.5818;
  double VolumeX_Al4_11 = 26.45 + 0.5818;
  double VolumeX_Al4_12 = 26.45 + 0.5818;
  double VolumeX_Al4_13 = 26.45 + 0.5818;
  double VolumeX_Al4_14 = 26.45 + 0.5818;


  // Y coordinate

  double VolumeY_Al1_0 = -30;
  double VolumeY_Al1_1 =   0;
  double VolumeY_Al1_2 =  30;

  double VolumeY_Al2_0 = -30;
  double VolumeY_Al2_1 =   0;
  double VolumeY_Al2_2 =  30;
  double VolumeY_Al2_3 = -30;
  double VolumeY_Al2_4 =   0;
  double VolumeY_Al2_5 =  30;

  double VolumeY_Al3_0 = -60;
  double VolumeY_Al3_1 = -30;
  double VolumeY_Al3_2 =   0;
  double VolumeY_Al3_3 =  30;
  double VolumeY_Al3_4 =  60;
  double VolumeY_Al3_5 = -60;
  double VolumeY_Al3_6 = -30;
  double VolumeY_Al3_7 =   0;
  double VolumeY_Al3_8 =  30;
  double VolumeY_Al3_9 =  60;

  double VolumeY_Al4_0 = -60;
  double VolumeY_Al4_1 = -30;
  double VolumeY_Al4_2 =   0;
  double VolumeY_Al4_3 =  30;
  double VolumeY_Al4_4 =  60;
  double VolumeY_Al4_5 = -60;
  double VolumeY_Al4_6 = -30;
  double VolumeY_Al4_7 =   0;
  double VolumeY_Al4_8 =  30;
  double VolumeY_Al4_9 =  60;
  double VolumeY_Al4_10 = -60;
  double VolumeY_Al4_11 = -30;
  double VolumeY_Al4_12 =   0;
  double VolumeY_Al4_13 =  30;
  double VolumeY_Al4_14 =  60;

  ///////////////////////////////////////////////
  double VolumeZ_Al1_0 = 0;
  double VolumeZ_Al1_1 = 0;
  double VolumeZ_Al1_2 = 0;

  double VolumeZ_Al2_0 = 50.;
  double VolumeZ_Al2_1 = 50.;
  double VolumeZ_Al2_2 = 50.;
  double VolumeZ_Al2_3 = 50.-4.43;
  double VolumeZ_Al2_4 = 50.-4.43;
  double VolumeZ_Al2_5 = 50.-4.43;

  double VolumeZ_Al3_0 = 100.;
  double VolumeZ_Al3_1 = 100.;
  double VolumeZ_Al3_2 = 100.;
  double VolumeZ_Al3_3 = 100.;
  double VolumeZ_Al3_4 = 100.;
  double VolumeZ_Al3_5 = 100.-4.43;
  double VolumeZ_Al3_6 = 100.-4.43;
  double VolumeZ_Al3_7 = 100.-4.43;
  double VolumeZ_Al3_8 = 100.-4.43;
  double VolumeZ_Al3_9 = 100.-4.43;


  double VolumeZ_Al4_0 = 150.0;
  double VolumeZ_Al4_1 = 150.0;
  double VolumeZ_Al4_2 = 150.0;
  double VolumeZ_Al4_3 = 150.0;
  double VolumeZ_Al4_4 = 150.0;
  double VolumeZ_Al4_5 = 150.0-4.43;
  double VolumeZ_Al4_6 = 150.0-4.43;
  double VolumeZ_Al4_7 = 150.0-4.43;
  double VolumeZ_Al4_8 = 150.0-4.43;
  double VolumeZ_Al4_9 = 150.0-4.43;
  double VolumeZ_Al4_10 = 150.0-9.44;
  double VolumeZ_Al4_11 = 150.0-9.44;
  double VolumeZ_Al4_12 = 150.0-9.44;
  double VolumeZ_Al4_13 = 150.0-9.44;
  double VolumeZ_Al4_14 = 150.0-9.44;

  ////////////////////// Some diagnostics: //////////////////////////
  
  cout<<"VolumeY_Al3_3-VolumeY_Al3_2 = "<<VolumeY_Al3_3-VolumeY_Al3_2<<endl;
  cout<<"VolumeY_Al3_2-VolumeY_Al3_1 = "<<VolumeY_Al3_2-VolumeY_Al3_1<<endl;
  cout<<"VolumeY_Al4_3-VolumeY_Al4_2 = "<<VolumeY_Al4_3-VolumeY_Al4_2<<endl;
  cout<<"VolumeY_Al4_2-VolumeY_Al4_1 = "<<VolumeY_Al4_2-VolumeY_Al4_1<<endl;

  ///////////////////////////////////////////////////////////////////////////  

  fVolumeX[0] = VolumeX_Al1_0;
  fVolumeX[1] = VolumeX_Al1_1;
  fVolumeX[2] = VolumeX_Al1_2;

  fVolumeX[3] = VolumeX_Al2_0;
  fVolumeX[4] = VolumeX_Al2_1;
  fVolumeX[5] = VolumeX_Al2_2;
  fVolumeX[6] = VolumeX_Al2_3;
  fVolumeX[7] = VolumeX_Al2_4;
  fVolumeX[8] = VolumeX_Al2_5;

  fVolumeX[9] =  VolumeX_Al3_0;
  fVolumeX[10] = VolumeX_Al3_1;
  fVolumeX[11] = VolumeX_Al3_2;
  fVolumeX[12] = VolumeX_Al3_3;
  fVolumeX[13] = VolumeX_Al3_4;
  fVolumeX[14] = VolumeX_Al3_5;
  fVolumeX[15] = VolumeX_Al3_6;
  fVolumeX[16] = VolumeX_Al3_7;
  fVolumeX[17] = VolumeX_Al3_8;
  fVolumeX[18] = VolumeX_Al3_9;

  fVolumeX[19] = VolumeX_Al4_0;
  fVolumeX[20] = VolumeX_Al4_1;
  fVolumeX[21] = VolumeX_Al4_2;
  fVolumeX[22] = VolumeX_Al4_3;
  fVolumeX[23] = VolumeX_Al4_4;
  fVolumeX[24] = VolumeX_Al4_5;
  fVolumeX[25] = VolumeX_Al4_6;
  fVolumeX[26] = VolumeX_Al4_7;
  fVolumeX[27] = VolumeX_Al4_8;
  fVolumeX[28] = VolumeX_Al4_9;
  fVolumeX[29] = VolumeX_Al4_10;
  fVolumeX[30] = VolumeX_Al4_11;
  fVolumeX[31] = VolumeX_Al4_12;
  fVolumeX[32] = VolumeX_Al4_13;
  fVolumeX[33] = VolumeX_Al4_14;

  ///////////////////////////////////////

  fVolumeY[0] = VolumeY_Al1_0;
  fVolumeY[1] = VolumeY_Al1_1;
  fVolumeY[2] = VolumeY_Al1_2;

  fVolumeY[3] = VolumeY_Al2_0;
  fVolumeY[4] = VolumeY_Al2_1;
  fVolumeY[5] = VolumeY_Al2_2;
  fVolumeY[6] = VolumeY_Al2_3;
  fVolumeY[7] = VolumeY_Al2_4;
  fVolumeY[8] = VolumeY_Al2_5;

  fVolumeY[9] =  VolumeY_Al3_0;
  fVolumeY[10] = VolumeY_Al3_1;
  fVolumeY[11] = VolumeY_Al3_2;
  fVolumeY[12] = VolumeY_Al3_3;
  fVolumeY[13] = VolumeY_Al3_4;
  fVolumeY[14] = VolumeY_Al3_5;
  fVolumeY[15] = VolumeY_Al3_6;
  fVolumeY[16] = VolumeY_Al3_7;
  fVolumeY[17] = VolumeY_Al3_8;
  fVolumeY[18] = VolumeY_Al3_9;

  fVolumeY[19] = VolumeY_Al4_0;
  fVolumeY[20] = VolumeY_Al4_1;
  fVolumeY[21] = VolumeY_Al4_2;
  fVolumeY[22] = VolumeY_Al4_3;
  fVolumeY[23] = VolumeY_Al4_4;
  fVolumeY[24] = VolumeY_Al4_5;
  fVolumeY[25] = VolumeY_Al4_6;
  fVolumeY[26] = VolumeY_Al4_7;
  fVolumeY[27] = VolumeY_Al4_8;
  fVolumeY[28] = VolumeY_Al4_9;
  fVolumeY[29] = VolumeY_Al4_10;
  fVolumeY[30] = VolumeY_Al4_11;
  fVolumeY[31] = VolumeY_Al4_12;
  fVolumeY[32] = VolumeY_Al4_13;
  fVolumeY[33] = VolumeY_Al4_14;

  ///////////////////////////////////////////

  fVolumeZ[0] = VolumeZ_Al1_0;
  fVolumeZ[1] = VolumeZ_Al1_1;
  fVolumeZ[2] = VolumeZ_Al1_2;

  fVolumeZ[3] = VolumeZ_Al2_0;
  fVolumeZ[4] = VolumeZ_Al2_1;
  fVolumeZ[5] = VolumeZ_Al2_2;
  fVolumeZ[6] = VolumeZ_Al2_3;
  fVolumeZ[7] = VolumeZ_Al2_4;
  fVolumeZ[8] = VolumeZ_Al2_5;

  fVolumeZ[9] =  VolumeZ_Al3_0;
  fVolumeZ[10] = VolumeZ_Al3_1;
  fVolumeZ[11] = VolumeZ_Al3_2;
  fVolumeZ[12] = VolumeZ_Al3_3;
  fVolumeZ[13] = VolumeZ_Al3_4;
  fVolumeZ[14] = VolumeZ_Al3_5;
  fVolumeZ[15] = VolumeZ_Al3_6;
  fVolumeZ[16] = VolumeZ_Al3_7;
  fVolumeZ[17] = VolumeZ_Al3_8;
  fVolumeZ[18] = VolumeZ_Al3_9;

  fVolumeZ[19] = VolumeZ_Al4_0;
  fVolumeZ[20] = VolumeZ_Al4_1;
  fVolumeZ[21] = VolumeZ_Al4_2;
  fVolumeZ[22] = VolumeZ_Al4_3;
  fVolumeZ[23] = VolumeZ_Al4_4;
  fVolumeZ[24] = VolumeZ_Al4_5;
  fVolumeZ[25] = VolumeZ_Al4_6;
  fVolumeZ[26] = VolumeZ_Al4_7;
  fVolumeZ[27] = VolumeZ_Al4_8;
  fVolumeZ[28] = VolumeZ_Al4_9;
  fVolumeZ[29] = VolumeZ_Al4_10;
  fVolumeZ[30] = VolumeZ_Al4_11;
  fVolumeZ[31] = VolumeZ_Al4_12;
  fVolumeZ[32] = VolumeZ_Al4_13;
  fVolumeZ[33] = VolumeZ_Al4_14;

}

//____________________________________________________________________________________________
void Na61AlVdArmParameters::SetupSensorGeometry_Jura_pPb2022(Int_t run_id)
{
  cout<<"Na61AlVdArmParameters::SetupSensorGeometry_Jura_pPb2022: Setting sensor geometry for AlVd detector for run: "<<run_id<<endl;

  ////////////// X coordinate
double VolumeX_Al1_0 = 0;
double VolumeX_Al1_1 = 0;
double VolumeX_Al1_2 = 0;
double VolumeX_Al2_0 = 0.0286684;
double VolumeX_Al2_1 = 0.029687;
double VolumeX_Al2_2 = 0.0269481;
double VolumeX_Al2_3 = 12.2698;
double VolumeX_Al2_4 = 12.2881;
double VolumeX_Al2_5 = 12.2804;
double VolumeX_Al3_0 = 0.705516;
double VolumeX_Al3_1 = 0.719004;
double VolumeX_Al3_2 = 0.671605;
double VolumeX_Al3_3 = 0.662531;
double VolumeX_Al3_4 = 0.637535;
double VolumeX_Al3_5 = 12.9479;
double VolumeX_Al3_6 = 12.955;
double VolumeX_Al3_7 = 12.957;
double VolumeX_Al3_8 = 12.911;
double VolumeX_Al3_9 = 12.9379;
double VolumeX_Al4_0 = 2.24709;
double VolumeX_Al4_1 = 2.18248;
double VolumeX_Al4_2 = 2.17962;
double VolumeX_Al4_3 = 2.20343;
double VolumeX_Al4_4 = 2.25218;
double VolumeX_Al4_5 = 14.3444;
double VolumeX_Al4_6 = 14.3674;
double VolumeX_Al4_7 = 14.3523;
double VolumeX_Al4_8 = 14.3471;
double VolumeX_Al4_9 = 14.31;
double VolumeX_Al4_10 = 26.687;
double VolumeX_Al4_11 = 26.6889;
double VolumeX_Al4_12 = 26.6813;
double VolumeX_Al4_13 = 26.6882;
double VolumeX_Al4_14 = 26.7905;


  // Y coordinate
double VolumeY_Al1_0 = -30;
double VolumeY_Al1_1 = 0;
double VolumeY_Al1_2 = 30;
double VolumeY_Al2_0 = -30.1614;
double VolumeY_Al2_1 = -0.00477168;
double VolumeY_Al2_2 = 30.1451;
double VolumeY_Al2_3 = -30.1805;
double VolumeY_Al2_4 = -0.0415238;
double VolumeY_Al2_5 = 30.1079;
double VolumeY_Al3_0 = -60.4058;
double VolumeY_Al3_1 = -30.2556;
double VolumeY_Al3_2 = -0.0738215;
double VolumeY_Al3_3 = 30.0917;
double VolumeY_Al3_4 = 60.3259;
double VolumeY_Al3_5 = -60.4099;
double VolumeY_Al3_6 = -30.2643;
double VolumeY_Al3_7 = -0.124073;
double VolumeY_Al3_8 = 30.001;
double VolumeY_Al3_9 = 60.1553;
double VolumeY_Al4_0 = -60.3803;
double VolumeY_Al4_1 = -30.1556;
double VolumeY_Al4_2 = 0.0114818;
double VolumeY_Al4_3 = 30.1726;
double VolumeY_Al4_4 = 60.5315;
double VolumeY_Al4_5 = -60.3086;
double VolumeY_Al4_6 = -30.141;
double VolumeY_Al4_7 = 0.0140532;
double VolumeY_Al4_8 = 30.1588;
double VolumeY_Al4_9 = 60.3037;
double VolumeY_Al4_10 = -60.2718;
double VolumeY_Al4_11 = -30.1199;
double VolumeY_Al4_12 = 0.0324327;
double VolumeY_Al4_13 = 30.1808;
double VolumeY_Al4_14 = 60.4196;


////// Z coordinate
double VolumeZ_Al1_0 = 0;
double VolumeZ_Al1_1 = 0;
double VolumeZ_Al1_2 = 0;
double VolumeZ_Al2_0 = 50.0746;
double VolumeZ_Al2_1 = 50.0259;
double VolumeZ_Al2_2 = 50.0007;
double VolumeZ_Al2_3 = 45.4908;
double VolumeZ_Al2_4 = 45.5634;
double VolumeZ_Al2_5 = 45.5708;
double VolumeZ_Al3_0 = 99.8579;
double VolumeZ_Al3_1 = 99.9907;
double VolumeZ_Al3_2 = 99.7474;
double VolumeZ_Al3_3 = 99.9193;
double VolumeZ_Al3_4 = 100.242;
double VolumeZ_Al3_5 = 95.609;
double VolumeZ_Al3_6 = 95.6988;
double VolumeZ_Al3_7 = 95.7421;
double VolumeZ_Al3_8 = 95.5775;
double VolumeZ_Al3_9 = 95.731;
double VolumeZ_Al4_0 = 150.946;
double VolumeZ_Al4_1 = 150.505;
double VolumeZ_Al4_2 = 150.383;
double VolumeZ_Al4_3 = 150.346;
double VolumeZ_Al4_4 = 151.267;
double VolumeZ_Al4_5 = 145.941;
double VolumeZ_Al4_6 = 145.836;
double VolumeZ_Al4_7 = 145.648;
double VolumeZ_Al4_8 = 145.544;
double VolumeZ_Al4_9 = 145.575;
double VolumeZ_Al4_10 = 140.411;
double VolumeZ_Al4_11 = 140.352;
double VolumeZ_Al4_12 = 140.188;
double VolumeZ_Al4_13 = 140.135;
double VolumeZ_Al4_14 = 140.547;


///////////////////////////////// some diagnosis //////////////////////////

  cout<<"VolumeY_Al3_4-VolumeY_Al3_3 = "<<VolumeY_Al3_4-VolumeY_Al3_3<<endl;
  cout<<"VolumeY_Al3_3-VolumeY_Al3_2 = "<<VolumeY_Al3_3-VolumeY_Al3_2<<endl;
  cout<<"VolumeY_Al3_2-VolumeY_Al3_1 = "<<VolumeY_Al3_2-VolumeY_Al3_1<<endl;
  cout<<"VolumeY_Al3_1-VolumeY_Al3_0 = "<<VolumeY_Al3_1-VolumeY_Al3_0<<endl;

  cout<<"VolumeY_Al4_4-VolumeY_Al4_3 = "<<VolumeY_Al4_4-VolumeY_Al3_3<<endl;
  cout<<"VolumeY_Al4_3-VolumeY_Al4_2 = "<<VolumeY_Al4_3-VolumeY_Al3_2<<endl;
  cout<<"VolumeY_Al4_2-VolumeY_Al4_1 = "<<VolumeY_Al4_2-VolumeY_Al3_1<<endl;
  cout<<"VolumeY_Al4_1-VolumeY_Al4_0 = "<<VolumeY_Al4_1-VolumeY_Al3_0<<endl;


  ///////////////////////////////////////////////////////////////////////////  

  fVolumeX[0] = VolumeX_Al1_0;
  fVolumeX[1] = VolumeX_Al1_1;
  fVolumeX[2] = VolumeX_Al1_2;

  fVolumeX[3] = VolumeX_Al2_0;
  fVolumeX[4] = VolumeX_Al2_1;
  fVolumeX[5] = VolumeX_Al2_2;
  fVolumeX[6] = VolumeX_Al2_3;
  fVolumeX[7] = VolumeX_Al2_4;
  fVolumeX[8] = VolumeX_Al2_5;

  fVolumeX[9] =  VolumeX_Al3_0;
  fVolumeX[10] = VolumeX_Al3_1;
  fVolumeX[11] = VolumeX_Al3_2;
  fVolumeX[12] = VolumeX_Al3_3;
  fVolumeX[13] = VolumeX_Al3_4;
  fVolumeX[14] = VolumeX_Al3_5;
  fVolumeX[15] = VolumeX_Al3_6;
  fVolumeX[16] = VolumeX_Al3_7;
  fVolumeX[17] = VolumeX_Al3_8;
  fVolumeX[18] = VolumeX_Al3_9;

  fVolumeX[19] = VolumeX_Al4_0;
  fVolumeX[20] = VolumeX_Al4_1;
  fVolumeX[21] = VolumeX_Al4_2;
  fVolumeX[22] = VolumeX_Al4_3;
  fVolumeX[23] = VolumeX_Al4_4;
  fVolumeX[24] = VolumeX_Al4_5;
  fVolumeX[25] = VolumeX_Al4_6;
  fVolumeX[26] = VolumeX_Al4_7;
  fVolumeX[27] = VolumeX_Al4_8;
  fVolumeX[28] = VolumeX_Al4_9;
  fVolumeX[29] = VolumeX_Al4_10;
  fVolumeX[30] = VolumeX_Al4_11;
  fVolumeX[31] = VolumeX_Al4_12;
  fVolumeX[32] = VolumeX_Al4_13;
  fVolumeX[33] = VolumeX_Al4_14;

  ///////////////////////////////////////

  fVolumeY[0] = VolumeY_Al1_0;
  fVolumeY[1] = VolumeY_Al1_1;
  fVolumeY[2] = VolumeY_Al1_2;

  fVolumeY[3] = VolumeY_Al2_0;
  fVolumeY[4] = VolumeY_Al2_1;
  fVolumeY[5] = VolumeY_Al2_2;
  fVolumeY[6] = VolumeY_Al2_3;
  fVolumeY[7] = VolumeY_Al2_4;
  fVolumeY[8] = VolumeY_Al2_5;

  fVolumeY[9] =  VolumeY_Al3_0;
  fVolumeY[10] = VolumeY_Al3_1;
  fVolumeY[11] = VolumeY_Al3_2;
  fVolumeY[12] = VolumeY_Al3_3;
  fVolumeY[13] = VolumeY_Al3_4;
  fVolumeY[14] = VolumeY_Al3_5;
  fVolumeY[15] = VolumeY_Al3_6;
  fVolumeY[16] = VolumeY_Al3_7;
  fVolumeY[17] = VolumeY_Al3_8;
  fVolumeY[18] = VolumeY_Al3_9;

  fVolumeY[19] = VolumeY_Al4_0;
  fVolumeY[20] = VolumeY_Al4_1;
  fVolumeY[21] = VolumeY_Al4_2;
  fVolumeY[22] = VolumeY_Al4_3;
  fVolumeY[23] = VolumeY_Al4_4;
  fVolumeY[24] = VolumeY_Al4_5;
  fVolumeY[25] = VolumeY_Al4_6;
  fVolumeY[26] = VolumeY_Al4_7;
  fVolumeY[27] = VolumeY_Al4_8;
  fVolumeY[28] = VolumeY_Al4_9;
  fVolumeY[29] = VolumeY_Al4_10;
  fVolumeY[30] = VolumeY_Al4_11;
  fVolumeY[31] = VolumeY_Al4_12;
  fVolumeY[32] = VolumeY_Al4_13;
  fVolumeY[33] = VolumeY_Al4_14;

  ///////////////////////////////////////////

  fVolumeZ[0] = VolumeZ_Al1_0;
  fVolumeZ[1] = VolumeZ_Al1_1;
  fVolumeZ[2] = VolumeZ_Al1_2;

  fVolumeZ[3] = VolumeZ_Al2_0;
  fVolumeZ[4] = VolumeZ_Al2_1;
  fVolumeZ[5] = VolumeZ_Al2_2;
  fVolumeZ[6] = VolumeZ_Al2_3;
  fVolumeZ[7] = VolumeZ_Al2_4;
  fVolumeZ[8] = VolumeZ_Al2_5;

  fVolumeZ[9] =  VolumeZ_Al3_0;
  fVolumeZ[10] = VolumeZ_Al3_1;
  fVolumeZ[11] = VolumeZ_Al3_2;
  fVolumeZ[12] = VolumeZ_Al3_3;
  fVolumeZ[13] = VolumeZ_Al3_4;
  fVolumeZ[14] = VolumeZ_Al3_5;
  fVolumeZ[15] = VolumeZ_Al3_6;
  fVolumeZ[16] = VolumeZ_Al3_7;
  fVolumeZ[17] = VolumeZ_Al3_8;
  fVolumeZ[18] = VolumeZ_Al3_9;

  fVolumeZ[19] = VolumeZ_Al4_0;
  fVolumeZ[20] = VolumeZ_Al4_1;
  fVolumeZ[21] = VolumeZ_Al4_2;
  fVolumeZ[22] = VolumeZ_Al4_3;
  fVolumeZ[23] = VolumeZ_Al4_4;
  fVolumeZ[24] = VolumeZ_Al4_5;
  fVolumeZ[25] = VolumeZ_Al4_6;
  fVolumeZ[26] = VolumeZ_Al4_7;
  fVolumeZ[27] = VolumeZ_Al4_8;
  fVolumeZ[28] = VolumeZ_Al4_9;
  fVolumeZ[29] = VolumeZ_Al4_10;
  fVolumeZ[30] = VolumeZ_Al4_11;
  fVolumeZ[31] = VolumeZ_Al4_12;
  fVolumeZ[32] = VolumeZ_Al4_13;
  fVolumeZ[33] = VolumeZ_Al4_14;


  for(Int_t i=0;i<34;i++){
    fVolumeX[i] = fVolumeX[i] - VolumeX_Al1_1; 
    fVolumeY[i] = fVolumeY[i] - VolumeY_Al1_1; 
    fVolumeZ[i] = fVolumeZ[i] - VolumeZ_Al1_1;
  } 

 
  fVolumeZ[fSensorId] = fVolumeZ[fSensorId] + fDz;


}

//____________________________________________________________________________________________
void Na61AlVdArmParameters::SetupSensorGeometry_SaleveNominal(Int_t run_id)
{
  cout<<"Na61AlVdArmParameters::SetupSensorGeometry_SaleveNominal: Setting sensor geometry for AlVd detector for run: "<<run_id<<endl;

  ////////////// X coordinate
  double VolumeX_Al1_0 = 0;
  double VolumeX_Al1_1 = 0;
  double VolumeX_Al1_2 = 0;

  double VolumeX_Al2_0 = 0;
  double VolumeX_Al2_1 = 0;
  double VolumeX_Al2_2 = 0;
  double VolumeX_Al2_3 = -12.56;
  double VolumeX_Al2_4 = -12.56;
  double VolumeX_Al2_5 = -12.56;

  double VolumeX_Al3_0 = 0 - 0.625;
  double VolumeX_Al3_1 = 0 - 0.625;
  double VolumeX_Al3_2 = 0 - 0.625;
  double VolumeX_Al3_3 = 0 - 0.625;
  double VolumeX_Al3_4 = 0 - 0.625;
  double VolumeX_Al3_5 = -12.56 - 0.625;
  double VolumeX_Al3_6 = -12.56 - 0.625;
  double VolumeX_Al3_7 = -12.56 - 0.625;
  double VolumeX_Al3_8 = -12.56 - 0.625;
  double VolumeX_Al3_9 = -12.56 - 0.625;

  double VolumeX_Al4_0 = -1.5 - 0.5818;
  double VolumeX_Al4_1 = -1.5 - 0.5818;
  double VolumeX_Al4_2 = -1.5 - 0.5818;
  double VolumeX_Al4_3 = -1.5 - 0.5818;
  double VolumeX_Al4_4 = -1.5 - 0.5818;
  double VolumeX_Al4_5 = -14.06 - 0.5818;
  double VolumeX_Al4_6 = -14.06 - 0.5818;
  double VolumeX_Al4_7 = -14.06 - 0.5818;
  double VolumeX_Al4_8 = -14.06 - 0.5818;
  double VolumeX_Al4_9 = -14.06 - 0.5818;
  double VolumeX_Al4_10 = -26.45 - 0.5818;
  double VolumeX_Al4_11 = -26.45 - 0.5818;
  double VolumeX_Al4_12 = -26.45 - 0.5818;
  double VolumeX_Al4_13 = -26.45 - 0.5818;
  double VolumeX_Al4_14 = -26.45 - 0.5818;


  // Y coordinate
  double VolumeY_Al1_0 = -30;
  double VolumeY_Al1_1 =   0;
  double VolumeY_Al1_2 =  30;

  double VolumeY_Al2_0 = -30;
  double VolumeY_Al2_1 =   0;
  double VolumeY_Al2_2 =  30;
  double VolumeY_Al2_3 = -30;
  double VolumeY_Al2_4 =   0;
  double VolumeY_Al2_5 =  30;

  double VolumeY_Al3_0 = -60;
  double VolumeY_Al3_1 = -30;
  double VolumeY_Al3_2 =   0;
  double VolumeY_Al3_3 =  30;
  double VolumeY_Al3_4 =  60;
  double VolumeY_Al3_5 = -60;
  double VolumeY_Al3_6 = -30;
  double VolumeY_Al3_7 =   0;
  double VolumeY_Al3_8 =  30;
  double VolumeY_Al3_9 =  60;

  double VolumeY_Al4_0 = -60;
  double VolumeY_Al4_1 = -30;
  double VolumeY_Al4_2 =   0;
  double VolumeY_Al4_3 =  30;
  double VolumeY_Al4_4 =  60;
  double VolumeY_Al4_5 = -60;
  double VolumeY_Al4_6 = -30;
  double VolumeY_Al4_7 =   0;
  double VolumeY_Al4_8 =  30;
  double VolumeY_Al4_9 =  60;
  double VolumeY_Al4_10 = -60;
  double VolumeY_Al4_11 = -30;
  double VolumeY_Al4_12 =   0;
  double VolumeY_Al4_13 =  30;
  double VolumeY_Al4_14 =  60;

  // Z coordinates
  double VolumeZ_Al1_0 = 0;
  double VolumeZ_Al1_1 = 0;
  double VolumeZ_Al1_2 = 0;

  double VolumeZ_Al2_0 = 50.;
  double VolumeZ_Al2_1 = 50.;
  double VolumeZ_Al2_2 = 50.;
  double VolumeZ_Al2_3 = 50.-4.43;
  double VolumeZ_Al2_4 = 50.-4.43;
  double VolumeZ_Al2_5 = 50.-4.43;

  double VolumeZ_Al3_0 = 100.;
  double VolumeZ_Al3_1 = 100.;
  double VolumeZ_Al3_2 = 100.;
  double VolumeZ_Al3_3 = 100.;
  double VolumeZ_Al3_4 = 100.;
  double VolumeZ_Al3_5 = 100.-4.43;
  double VolumeZ_Al3_6 = 100.-4.43;
  double VolumeZ_Al3_7 = 100.-4.43;
  double VolumeZ_Al3_8 = 100.-4.43;
  double VolumeZ_Al3_9 = 100.-4.43;


  double VolumeZ_Al4_0 = 150.0;
  double VolumeZ_Al4_1 = 150.0;
  double VolumeZ_Al4_2 = 150.0;
  double VolumeZ_Al4_3 = 150.0;
  double VolumeZ_Al4_4 = 150.0;
  double VolumeZ_Al4_5 = 150.0-4.43;
  double VolumeZ_Al4_6 = 150.0-4.43;
  double VolumeZ_Al4_7 = 150.0-4.43;
  double VolumeZ_Al4_8 = 150.0-4.43;
  double VolumeZ_Al4_9 = 150.0-4.43;
  double VolumeZ_Al4_10 = 150.0-9.44;
  double VolumeZ_Al4_11 = 150.0-9.44;
  double VolumeZ_Al4_12 = 150.0-9.44;
  double VolumeZ_Al4_13 = 150.0-9.44;
  double VolumeZ_Al4_14 = 150.0-9.44;


  /////////////////////////// some diagnostics: /////////////////////////////

  cout<<"VolumeY_Al3_3-VolumeY_Al3_2 = "<<VolumeY_Al3_3-VolumeY_Al3_2<<endl;
  cout<<"VolumeY_Al3_2-VolumeY_Al3_1 = "<<VolumeY_Al3_2-VolumeY_Al3_1<<endl;
  cout<<"VolumeY_Al4_3-VolumeY_Al4_2 = "<<VolumeY_Al4_3-VolumeY_Al4_2<<endl;
  cout<<"VolumeY_Al4_2-VolumeY_Al4_1 = "<<VolumeY_Al4_2-VolumeY_Al4_1<<endl;

  ///////////////////////////////////////////////////////////////////////////  

  
  fVolumeX[0] = VolumeX_Al1_0;
  fVolumeX[1] = VolumeX_Al1_1;
  fVolumeX[2] = VolumeX_Al1_2;

  fVolumeX[3] = VolumeX_Al2_0;
  fVolumeX[4] = VolumeX_Al2_1;
  fVolumeX[5] = VolumeX_Al2_2;
  fVolumeX[6] = VolumeX_Al2_3;
  fVolumeX[7] = VolumeX_Al2_4;
  fVolumeX[8] = VolumeX_Al2_5;

  fVolumeX[9] =  VolumeX_Al3_0;
  fVolumeX[10] = VolumeX_Al3_1;
  fVolumeX[11] = VolumeX_Al3_2;
  fVolumeX[12] = VolumeX_Al3_3;
  fVolumeX[13] = VolumeX_Al3_4;
  fVolumeX[14] = VolumeX_Al3_5;
  fVolumeX[15] = VolumeX_Al3_6;
  fVolumeX[16] = VolumeX_Al3_7;
  fVolumeX[17] = VolumeX_Al3_8;
  fVolumeX[18] = VolumeX_Al3_9;

  fVolumeX[19] = VolumeX_Al4_0;
  fVolumeX[20] = VolumeX_Al4_1;
  fVolumeX[21] = VolumeX_Al4_2;
  fVolumeX[22] = VolumeX_Al4_3;
  fVolumeX[23] = VolumeX_Al4_4;
  fVolumeX[24] = VolumeX_Al4_5;
  fVolumeX[25] = VolumeX_Al4_6;
  fVolumeX[26] = VolumeX_Al4_7;
  fVolumeX[27] = VolumeX_Al4_8;
  fVolumeX[28] = VolumeX_Al4_9;
  fVolumeX[29] = VolumeX_Al4_10;
  fVolumeX[30] = VolumeX_Al4_11;
  fVolumeX[31] = VolumeX_Al4_12;
  fVolumeX[32] = VolumeX_Al4_13;
  fVolumeX[33] = VolumeX_Al4_14;

  ///////////////////////////////////////

  fVolumeY[0] = VolumeY_Al1_0;
  fVolumeY[1] = VolumeY_Al1_1;
  fVolumeY[2] = VolumeY_Al1_2;

  fVolumeY[3] = VolumeY_Al2_0;
  fVolumeY[4] = VolumeY_Al2_1;
  fVolumeY[5] = VolumeY_Al2_2;
  fVolumeY[6] = VolumeY_Al2_3;
  fVolumeY[7] = VolumeY_Al2_4;
  fVolumeY[8] = VolumeY_Al2_5;

  fVolumeY[9] =  VolumeY_Al3_0;
  fVolumeY[10] = VolumeY_Al3_1;
  fVolumeY[11] = VolumeY_Al3_2;
  fVolumeY[12] = VolumeY_Al3_3;
  fVolumeY[13] = VolumeY_Al3_4;
  fVolumeY[14] = VolumeY_Al3_5;
  fVolumeY[15] = VolumeY_Al3_6;
  fVolumeY[16] = VolumeY_Al3_7;
  fVolumeY[17] = VolumeY_Al3_8;
  fVolumeY[18] = VolumeY_Al3_9;

  fVolumeY[19] = VolumeY_Al4_0;
  fVolumeY[20] = VolumeY_Al4_1;
  fVolumeY[21] = VolumeY_Al4_2;
  fVolumeY[22] = VolumeY_Al4_3;
  fVolumeY[23] = VolumeY_Al4_4;
  fVolumeY[24] = VolumeY_Al4_5;
  fVolumeY[25] = VolumeY_Al4_6;
  fVolumeY[26] = VolumeY_Al4_7;
  fVolumeY[27] = VolumeY_Al4_8;
  fVolumeY[28] = VolumeY_Al4_9;
  fVolumeY[29] = VolumeY_Al4_10;
  fVolumeY[30] = VolumeY_Al4_11;
  fVolumeY[31] = VolumeY_Al4_12;
  fVolumeY[32] = VolumeY_Al4_13;
  fVolumeY[33] = VolumeY_Al4_14;

  ///////////////////////////////////////////

  fVolumeZ[0] = VolumeZ_Al1_0;
  fVolumeZ[1] = VolumeZ_Al1_1;
  fVolumeZ[2] = VolumeZ_Al1_2;

  fVolumeZ[3] = VolumeZ_Al2_0;
  fVolumeZ[4] = VolumeZ_Al2_1;
  fVolumeZ[5] = VolumeZ_Al2_2;
  fVolumeZ[6] = VolumeZ_Al2_3;
  fVolumeZ[7] = VolumeZ_Al2_4;
  fVolumeZ[8] = VolumeZ_Al2_5;

  fVolumeZ[9] =  VolumeZ_Al3_0;
  fVolumeZ[10] = VolumeZ_Al3_1;
  fVolumeZ[11] = VolumeZ_Al3_2;
  fVolumeZ[12] = VolumeZ_Al3_3;
  fVolumeZ[13] = VolumeZ_Al3_4;
  fVolumeZ[14] = VolumeZ_Al3_5;
  fVolumeZ[15] = VolumeZ_Al3_6;
  fVolumeZ[16] = VolumeZ_Al3_7;
  fVolumeZ[17] = VolumeZ_Al3_8;
  fVolumeZ[18] = VolumeZ_Al3_9;

  fVolumeZ[19] = VolumeZ_Al4_0;
  fVolumeZ[20] = VolumeZ_Al4_1;
  fVolumeZ[21] = VolumeZ_Al4_2;
  fVolumeZ[22] = VolumeZ_Al4_3;
  fVolumeZ[23] = VolumeZ_Al4_4;
  fVolumeZ[24] = VolumeZ_Al4_5;
  fVolumeZ[25] = VolumeZ_Al4_6;
  fVolumeZ[26] = VolumeZ_Al4_7;
  fVolumeZ[27] = VolumeZ_Al4_8;
  fVolumeZ[28] = VolumeZ_Al4_9;
  fVolumeZ[29] = VolumeZ_Al4_10;
  fVolumeZ[30] = VolumeZ_Al4_11;
  fVolumeZ[31] = VolumeZ_Al4_12;
  fVolumeZ[32] = VolumeZ_Al4_13;
  fVolumeZ[33] = VolumeZ_Al4_14;


  for(Int_t i=0;i<34;i++){
    fVolumeX[i] = fVolumeX[i] - VolumeX_Al1_1; 
    fVolumeY[i] = fVolumeY[i] - VolumeY_Al1_1; 
    fVolumeZ[i] = fVolumeZ[i] - VolumeZ_Al1_1;
  } 

 
  fVolumeZ[fSensorId] = fVolumeZ[fSensorId] + fDz;


}

//____________________________________________________________________________________________
void Na61AlVdArmParameters::SetupSensorGeometry_Saleve_pPb2022(Int_t run_id)
{
  cout<<"Na61AlVdArmParameters::SetupSensorGeometry_Saleve_pPb2022: Setting sensor geometry for AlVd detector for run: "<<run_id<<endl;

  ////////////// X coordinate
double VolumeX_Al1_0 = 0;
double VolumeX_Al1_1 = 0;
double VolumeX_Al1_2 = 0;
double VolumeX_Al2_0 = 0.0152141;
double VolumeX_Al2_1 = 0.0149972;
double VolumeX_Al2_2 = 0.0223952;
double VolumeX_Al2_3 = -12.7761;
double VolumeX_Al2_4 = -12.794;
double VolumeX_Al2_5 = -12.7579;
double VolumeX_Al3_0 = -0.49;
double VolumeX_Al3_1 = -0.478075;
double VolumeX_Al3_2 = -0.504861;
double VolumeX_Al3_3 = -0.513511;
double VolumeX_Al3_4 = -0.5;
double VolumeX_Al3_5 = -13.3;
double VolumeX_Al3_6 = -13.3191;
double VolumeX_Al3_7 = -13.308;
double VolumeX_Al3_8 = -13.3423;
double VolumeX_Al3_9 = -13.3;
double VolumeX_Al4_0 = -2.06391;
double VolumeX_Al4_1 = -2.03882;
double VolumeX_Al4_2 = -2.05906;
double VolumeX_Al4_3 = -2.05589;
double VolumeX_Al4_4 = -2.04698;
double VolumeX_Al4_5 = -14.7723;
double VolumeX_Al4_6 = -14.8026;
double VolumeX_Al4_7 = -14.8083;
double VolumeX_Al4_8 = -14.8291;
double VolumeX_Al4_9 = -14.7983;
double VolumeX_Al4_10 = -27.2902;
double VolumeX_Al4_11 = -27.2118;
double VolumeX_Al4_12 = -27.1731;
double VolumeX_Al4_13 = -27.1979;
double VolumeX_Al4_14 = -27.229;

  // Y coordinate
double VolumeY_Al1_0 = -30;
double VolumeY_Al1_1 = 0;
double VolumeY_Al1_2 = 30;
double VolumeY_Al2_0 = -30.1577;
double VolumeY_Al2_1 = -0.0204851;
double VolumeY_Al2_2 = 30.1168;
double VolumeY_Al2_3 = -30.0853;
double VolumeY_Al2_4 = 0.0496958;
double VolumeY_Al2_5 = 30.1697;
double VolumeY_Al3_0 = -60.3;
double VolumeY_Al3_1 = -30.1959;
double VolumeY_Al3_2 = -0.0764792;
double VolumeY_Al3_3 = 30.0447;
double VolumeY_Al3_4 = 60.2;
double VolumeY_Al3_5 = -60.3;
double VolumeY_Al3_6 = -30.1428;
double VolumeY_Al3_7 = 0.0220242;
double VolumeY_Al3_8 = 30.1776;
double VolumeY_Al3_9 = 60.3;
double VolumeY_Al4_0 = -60.1851;
double VolumeY_Al4_1 = -30.0576;
double VolumeY_Al4_2 = 0.0503087;
double VolumeY_Al4_3 = 30.1622;
double VolumeY_Al4_4 = 60.2776;
double VolumeY_Al4_5 = -60.3085;
double VolumeY_Al4_6 = -30.1432;
double VolumeY_Al4_7 = 0.00437655;
double VolumeY_Al4_8 = 30.1535;
double VolumeY_Al4_9 = 60.2726;
double VolumeY_Al4_10 = -60.2646;
double VolumeY_Al4_11 = -30.0151;
double VolumeY_Al4_12 = 0.160529;
double VolumeY_Al4_13 = 30.3211;
double VolumeY_Al4_14 = 60.3641;

  // Z coordinates
double VolumeZ_Al1_0 = 0;
double VolumeZ_Al1_1 = 0;
double VolumeZ_Al1_2 = 0;
double VolumeZ_Al2_0 = 49.9997;
double VolumeZ_Al2_1 = 50.0755;
double VolumeZ_Al2_2 = 50.0293;
double VolumeZ_Al2_3 = 45.57;
double VolumeZ_Al2_4 = 45.7152;
double VolumeZ_Al2_5 = 45.5688;
double VolumeZ_Al3_0 = 100;
double VolumeZ_Al3_1 = 100;
double VolumeZ_Al3_2 = 100.188;
double VolumeZ_Al3_3 = 100.002;
double VolumeZ_Al3_4 = 100;
double VolumeZ_Al3_5 = 95.57;
double VolumeZ_Al3_6 = 95.5698;
double VolumeZ_Al3_7 = 95.4398;
double VolumeZ_Al3_8 = 95.5645;
double VolumeZ_Al3_9 = 95.57;
double VolumeZ_Al4_0 = 149.691;
double VolumeZ_Al4_1 = 149.707;
double VolumeZ_Al4_2 = 149.888;
double VolumeZ_Al4_3 = 149.633;
double VolumeZ_Al4_4 = 149.007;
double VolumeZ_Al4_5 = 145.591;
double VolumeZ_Al4_6 = 145.6;
double VolumeZ_Al4_7 = 145.517;
double VolumeZ_Al4_8 = 145.491;
double VolumeZ_Al4_9 = 145.184;
double VolumeZ_Al4_10 = 140.56;
double VolumeZ_Al4_11 = 140.152;
double VolumeZ_Al4_12 = 139.892;
double VolumeZ_Al4_13 = 140.024;
double VolumeZ_Al4_14 = 139.617;

///////////////////////////////// some diagnosis //////////////////////////

  cout<<"VolumeY_Al3_4-VolumeY_Al3_3 = "<<VolumeY_Al3_4-VolumeY_Al3_3<<endl;
  cout<<"VolumeY_Al3_3-VolumeY_Al3_2 = "<<VolumeY_Al3_3-VolumeY_Al3_2<<endl;
  cout<<"VolumeY_Al3_2-VolumeY_Al3_1 = "<<VolumeY_Al3_2-VolumeY_Al3_1<<endl;
  cout<<"VolumeY_Al3_1-VolumeY_Al3_0 = "<<VolumeY_Al3_1-VolumeY_Al3_0<<endl;

  cout<<"VolumeY_Al4_4-VolumeY_Al4_3 = "<<VolumeY_Al4_4-VolumeY_Al3_3<<endl;
  cout<<"VolumeY_Al4_3-VolumeY_Al4_2 = "<<VolumeY_Al4_3-VolumeY_Al3_2<<endl;
  cout<<"VolumeY_Al4_2-VolumeY_Al4_1 = "<<VolumeY_Al4_2-VolumeY_Al3_1<<endl;
  cout<<"VolumeY_Al4_1-VolumeY_Al4_0 = "<<VolumeY_Al4_1-VolumeY_Al3_0<<endl;

  ///////////////////////////////////////////////////////////////////////////  

  
  fVolumeX[0] = VolumeX_Al1_0;
  fVolumeX[1] = VolumeX_Al1_1;
  fVolumeX[2] = VolumeX_Al1_2;

  fVolumeX[3] = VolumeX_Al2_0;
  fVolumeX[4] = VolumeX_Al2_1;
  fVolumeX[5] = VolumeX_Al2_2;
  fVolumeX[6] = VolumeX_Al2_3;
  fVolumeX[7] = VolumeX_Al2_4;
  fVolumeX[8] = VolumeX_Al2_5;

  fVolumeX[9] =  VolumeX_Al3_0;
  fVolumeX[10] = VolumeX_Al3_1;
  fVolumeX[11] = VolumeX_Al3_2;
  fVolumeX[12] = VolumeX_Al3_3;
  fVolumeX[13] = VolumeX_Al3_4;
  fVolumeX[14] = VolumeX_Al3_5;
  fVolumeX[15] = VolumeX_Al3_6;
  fVolumeX[16] = VolumeX_Al3_7;
  fVolumeX[17] = VolumeX_Al3_8;
  fVolumeX[18] = VolumeX_Al3_9;

  fVolumeX[19] = VolumeX_Al4_0;
  fVolumeX[20] = VolumeX_Al4_1;
  fVolumeX[21] = VolumeX_Al4_2;
  fVolumeX[22] = VolumeX_Al4_3;
  fVolumeX[23] = VolumeX_Al4_4;
  fVolumeX[24] = VolumeX_Al4_5;
  fVolumeX[25] = VolumeX_Al4_6;
  fVolumeX[26] = VolumeX_Al4_7;
  fVolumeX[27] = VolumeX_Al4_8;
  fVolumeX[28] = VolumeX_Al4_9;
  fVolumeX[29] = VolumeX_Al4_10;
  fVolumeX[30] = VolumeX_Al4_11;
  fVolumeX[31] = VolumeX_Al4_12;
  fVolumeX[32] = VolumeX_Al4_13;
  fVolumeX[33] = VolumeX_Al4_14;

  ///////////////////////////////////////

  fVolumeY[0] = VolumeY_Al1_0;
  fVolumeY[1] = VolumeY_Al1_1;
  fVolumeY[2] = VolumeY_Al1_2;

  fVolumeY[3] = VolumeY_Al2_0;
  fVolumeY[4] = VolumeY_Al2_1;
  fVolumeY[5] = VolumeY_Al2_2;
  fVolumeY[6] = VolumeY_Al2_3;
  fVolumeY[7] = VolumeY_Al2_4;
  fVolumeY[8] = VolumeY_Al2_5;

  fVolumeY[9] =  VolumeY_Al3_0;
  fVolumeY[10] = VolumeY_Al3_1;
  fVolumeY[11] = VolumeY_Al3_2;
  fVolumeY[12] = VolumeY_Al3_3;
  fVolumeY[13] = VolumeY_Al3_4;
  fVolumeY[14] = VolumeY_Al3_5;
  fVolumeY[15] = VolumeY_Al3_6;
  fVolumeY[16] = VolumeY_Al3_7;
  fVolumeY[17] = VolumeY_Al3_8;
  fVolumeY[18] = VolumeY_Al3_9;

  fVolumeY[19] = VolumeY_Al4_0;
  fVolumeY[20] = VolumeY_Al4_1;
  fVolumeY[21] = VolumeY_Al4_2;
  fVolumeY[22] = VolumeY_Al4_3;
  fVolumeY[23] = VolumeY_Al4_4;
  fVolumeY[24] = VolumeY_Al4_5;
  fVolumeY[25] = VolumeY_Al4_6;
  fVolumeY[26] = VolumeY_Al4_7;
  fVolumeY[27] = VolumeY_Al4_8;
  fVolumeY[28] = VolumeY_Al4_9;
  fVolumeY[29] = VolumeY_Al4_10;
  fVolumeY[30] = VolumeY_Al4_11;
  fVolumeY[31] = VolumeY_Al4_12;
  fVolumeY[32] = VolumeY_Al4_13;
  fVolumeY[33] = VolumeY_Al4_14;

  ///////////////////////////////////////////

  fVolumeZ[0] = VolumeZ_Al1_0;
  fVolumeZ[1] = VolumeZ_Al1_1;
  fVolumeZ[2] = VolumeZ_Al1_2;

  fVolumeZ[3] = VolumeZ_Al2_0;
  fVolumeZ[4] = VolumeZ_Al2_1;
  fVolumeZ[5] = VolumeZ_Al2_2;
  fVolumeZ[6] = VolumeZ_Al2_3;
  fVolumeZ[7] = VolumeZ_Al2_4;
  fVolumeZ[8] = VolumeZ_Al2_5;

  fVolumeZ[9] =  VolumeZ_Al3_0;
  fVolumeZ[10] = VolumeZ_Al3_1;
  fVolumeZ[11] = VolumeZ_Al3_2;
  fVolumeZ[12] = VolumeZ_Al3_3;
  fVolumeZ[13] = VolumeZ_Al3_4;
  fVolumeZ[14] = VolumeZ_Al3_5;
  fVolumeZ[15] = VolumeZ_Al3_6;
  fVolumeZ[16] = VolumeZ_Al3_7;
  fVolumeZ[17] = VolumeZ_Al3_8;
  fVolumeZ[18] = VolumeZ_Al3_9;

  fVolumeZ[19] = VolumeZ_Al4_0;
  fVolumeZ[20] = VolumeZ_Al4_1;
  fVolumeZ[21] = VolumeZ_Al4_2;
  fVolumeZ[22] = VolumeZ_Al4_3;
  fVolumeZ[23] = VolumeZ_Al4_4;
  fVolumeZ[24] = VolumeZ_Al4_5;
  fVolumeZ[25] = VolumeZ_Al4_6;
  fVolumeZ[26] = VolumeZ_Al4_7;
  fVolumeZ[27] = VolumeZ_Al4_8;
  fVolumeZ[28] = VolumeZ_Al4_9;
  fVolumeZ[29] = VolumeZ_Al4_10;
  fVolumeZ[30] = VolumeZ_Al4_11;
  fVolumeZ[31] = VolumeZ_Al4_12;
  fVolumeZ[32] = VolumeZ_Al4_13;
  fVolumeZ[33] = VolumeZ_Al4_14;


}




//_____________________________________________________________
void Na61AlVdArmParameters::SetupOffsetsAndSigmas_JuraNominal(Int_t run_id)
{

  // dev cuts
  
  cout<<" Na61AlVdArmParameters::SetupOffsetsAndSigmas_JuraNominal: setting offsets and sigmas for run: "<<run_id<<endl;

  // offsets and sigmas related to match is deviations   
  // this parameters are used in reconstruction of simulation data
  // see definitions of DevNames in Na61AlVdTrackingInitModule::Init()
fOffsetx[0] = 0.00000;    fSigmax[0] =  0.0081830;
fOffsety[0] = 0.00000;    fSigmay[0] =  0.0088311;
fOffsetx[1] = 0.00000;    fSigmax[1] =  0.0117829;
fOffsety[1] = 0.00000;    fSigmay[1] =  0.0126218;
fOffsetx[2] = 0.00000;    fSigmax[2] =  0.0120318;
fOffsety[2] = 0.00000;    fSigmay[2] =  0.0126434;
fOffsetx[3] = 0.00000;    fSigmax[3] =  0.0127788;
fOffsety[3] = 0.00000;    fSigmay[3] =  0.0135693;
fOffsetx[4] = 0.00000;    fSigmax[4] =  0.0121559;
fOffsety[4] = 0.00000;    fSigmay[4] =  0.0132953;
fOffsetx[5] = 0.00000;    fSigmax[5] =  0.0114347;
fOffsety[5] = 0.00000;    fSigmay[5] =  0.0121771;
fOffsetx[6] = 0.00000;    fSigmax[6] =  0.0120442;
fOffsety[6] = 0.00000;    fSigmay[6] =  0.0124017;
fOffsetx[7] = 0.00000;    fSigmax[7] =  0.0124819;
fOffsety[7] = 0.00000;    fSigmay[7] =  0.0129045;
fOffsetx[8] = 0.00000;    fSigmax[8] =  0.0145894;
fOffsety[8] = 0.00000;    fSigmay[8] =  0.0148319;
fOffsetx[9] = 0.00000;    fSigmax[9] =  0.0139158;
fOffsety[9] = 0.00000;    fSigmay[9] =  0.0147842;
fOffsetx[10] = 0.00000;    fSigmax[10] =  0.0114380;
fOffsety[10] = 0.00000;    fSigmay[10] =  0.0117170;
fOffsetx[11] = 0.00000;    fSigmax[11] =  0.0146778;
fOffsety[11] = 0.00000;    fSigmay[11] =  0.0145518;
fOffsetx[12] = 0.00000;    fSigmax[12] =  0.0136237;
fOffsety[12] = 0.00000;    fSigmay[12] =  0.0151663;
fOffsetx[13] = 0.00000;    fSigmax[13] =  0.0166559;
fOffsety[13] = 0.00000;    fSigmay[13] =  0.0193654;
fOffsetx[14] = 0.00000;    fSigmax[14] =  0.0163712;
fOffsety[14] = 0.00000;    fSigmay[14] =  0.0172706;
fOffsetx[15] = 0.00000;    fSigmax[15] =  0.0071570;
fOffsety[15] = 0.00000;    fSigmay[15] =  0.0063084;
fOffsetx[16] = 0.00000;    fSigmax[16] =  0.0078875;
fOffsety[16] = 0.00000;    fSigmay[16] =  0.0082534;
fOffsetx[17] = 0.00000;    fSigmax[17] =  0.0098661;
fOffsety[17] = 0.00000;    fSigmay[17] =  0.0103857;
fOffsetx[18] = 0.00000;    fSigmax[18] =  0.0097136;
fOffsety[18] = 0.00000;    fSigmay[18] =  0.0102804;
fOffsetx[19] = 0.00000;    fSigmax[19] =  0.0094875;
fOffsety[19] = 0.00000;    fSigmay[19] =  0.0110667;
fOffsetx[20] = 0.00000;    fSigmax[20] =  0.0087451;
fOffsety[20] = 0.00000;    fSigmay[20] =  0.0109275;
fOffsetx[21] = 0.00000;    fSigmax[21] =  0.0186680;
fOffsety[21] = 0.00000;    fSigmay[21] =  0.0196205;
fOffsetx[22] = 0.00000;    fSigmax[22] =  0.0165328;
fOffsety[22] = 0.00000;    fSigmay[22] =  0.0237489;
fOffsetx[23] = 0.00000;    fSigmax[23] =  0.0246564;
fOffsety[23] = 0.00000;    fSigmay[23] =  0.0239062;
fOffsetx[24] = 0.00000;    fSigmax[24] =  0.0236755;
fOffsety[24] = 0.00000;    fSigmay[24] =  0.0408588;
fOffsetx[25] = 0.00000;    fSigmax[25] =  0.0113653;
fOffsety[25] = 0.00000;    fSigmay[25] =  0.0113806;
fOffsetx[26] = 0.00000;    fSigmax[26] =  0.0148538;
fOffsety[26] = 0.00000;    fSigmay[26] =  0.0161898;
fOffsetx[27] = 0.00000;    fSigmax[27] =  0.0152769;
fOffsety[27] = 0.00000;    fSigmay[27] =  0.0167876;
fOffsetx[28] = 0.00000;    fSigmax[28] =  0.0183926;
fOffsety[28] = 0.00000;    fSigmay[28] =  0.0209386;
fOffsetx[29] = 0.00000;    fSigmax[29] =  0.0179723;
fOffsety[29] = 0.00000;    fSigmay[29] =  0.0204218;
fOffsetx[30] = 0.00000;    fSigmax[30] =  0.0191416;
fOffsety[30] = 0.00000;    fSigmay[30] =  0.0225640;
fOffsetx[31] = 0.00000;    fSigmax[31] =  0.0220689;
fOffsety[31] = 0.00000;    fSigmay[31] =  0.0320766;
fOffsetx[32] = 0.00000;    fSigmax[32] =  0.0111282;
fOffsety[32] = 0.00000;    fSigmay[32] =  0.0119685;
fOffsetx[33] = 0.00000;    fSigmax[33] =  0.0162771;
fOffsety[33] = 0.00000;    fSigmay[33] =  0.0176218;
fOffsetx[34] = 0.00000;    fSigmax[34] =  0.0151228;
fOffsety[34] = 0.00000;    fSigmay[34] =  0.0170053;
fOffsetx[35] = 0.00000;    fSigmax[35] =  0.0183065;
fOffsety[35] = 0.00000;    fSigmay[35] =  0.0188756;
fOffsetx[36] = 0.00000;    fSigmax[36] =  0.0167719;
fOffsety[36] = 0.00000;    fSigmay[36] =  0.0218830;
fOffsetx[37] = 0.00000;    fSigmax[37] =  0.0210523;
fOffsety[37] = 0.00000;    fSigmay[37] =  0.0248084;
fOffsetx[38] = 0.00000;    fSigmax[38] =  0.0224151;
fOffsety[38] = 0.00000;    fSigmay[38] =  0.0317123;
fOffsetx[39] = 0.00000;    fSigmax[39] =  0.0120219;
fOffsety[39] = 0.00000;    fSigmay[39] =  0.0119852;
fOffsetx[40] = 0.00000;    fSigmax[40] =  0.0148893;
fOffsety[40] = 0.00000;    fSigmay[40] =  0.0162337;
fOffsetx[41] = 0.00000;    fSigmax[41] =  0.0151564;
fOffsety[41] = 0.00000;    fSigmay[41] =  0.0146136;
fOffsetx[42] = 0.00000;    fSigmax[42] =  0.0166262;
fOffsety[42] = 0.00000;    fSigmay[42] =  0.0192276;
fOffsetx[43] = 0.00000;    fSigmax[43] =  0.0229795;
fOffsety[43] = 0.00000;    fSigmay[43] =  0.0229538;
fOffsetx[44] = 0.00000;    fSigmax[44] =  0.0192721;
fOffsety[44] = 0.00000;    fSigmay[44] =  0.0278967;
fOffsetx[45] = 0.00000;    fSigmax[45] =  0.0298885;
fOffsety[45] = 0.00000;    fSigmay[45] =  0.0352376;
fOffsetx[46] = 0.00000;    fSigmax[46] =  0.0118183;
fOffsety[46] = 0.00000;    fSigmay[46] =  0.0122505;
fOffsetx[47] = 0.00000;    fSigmax[47] =  0.0118089;
fOffsety[47] = 0.00000;    fSigmay[47] =  0.0121104;
fOffsetx[48] = 0.00000;    fSigmax[48] =  0.0132968;
fOffsety[48] = 0.00000;    fSigmay[48] =  0.0128608;
fOffsetx[49] = 0.00000;    fSigmax[49] =  0.0143239;
fOffsety[49] = 0.00000;    fSigmay[49] =  0.0142901;
fOffsetx[50] = 0.00000;    fSigmax[50] =  0.0118276;
fOffsety[50] = 0.00000;    fSigmay[50] =  0.0136148;
fOffsetx[51] = 0.00000;    fSigmax[51] =  0.0127312;
fOffsety[51] = 0.00000;    fSigmay[51] =  0.0134487;
fOffsetx[52] = 0.00000;    fSigmax[52] =  0.0130201;
fOffsety[52] = 0.00000;    fSigmay[52] =  0.0135285;
fOffsetx[53] = 0.00000;    fSigmax[53] =  0.0137123;
fOffsety[53] = 0.00000;    fSigmay[53] =  0.0138798;
fOffsetx[54] = 0.00000;    fSigmax[54] =  0.0114311;
fOffsety[54] = 0.00000;    fSigmay[54] =  0.0127794;
fOffsetx[55] = 0.00000;    fSigmax[55] =  0.0122366;
fOffsety[55] = 0.00000;    fSigmay[55] =  0.0115744;
fOffsetx[56] = 0.00000;    fSigmax[56] =  0.0119618;
fOffsety[56] = 0.00000;    fSigmay[56] =  0.0122811;
fOffsetx[57] = 0.00000;    fSigmax[57] =  0.0125846;
fOffsety[57] = 0.00000;    fSigmay[57] =  0.0128048;
fOffsetx[58] = 0.00000;    fSigmax[58] =  0.0132080;
fOffsety[58] = 0.00000;    fSigmay[58] =  0.0135849;
fOffsetx[59] = 0.00000;    fSigmax[59] =  0.0133734;
fOffsety[59] = 0.00000;    fSigmay[59] =  0.0126713;
  

  fOffx[0] =  0.0;  fSigx[0] =  0.0015; 
  fOffy[0] =  0.0;  fSigy[0] =  0.0015;

  fOffx[1] =  0.0;  fSigx[1] =  0.0015; 
  fOffy[1] =  0.0;  fSigy[1] =  0.0015;

  fOffx[2] =  0.0;  fSigx[2] =  0.0015; 
  fOffy[2] =  0.0;  fSigy[2] =  0.0015;

  fOffx[3] =  0.0;  fSigx[3] =  0.0015; 
  fOffy[3] =  0.0;  fSigy[3] =  0.0015;

  // cuts on matching with primary vertex
  fOffx[4] =  0.0000002;  fSigx[4] =  0.0003999;
  fOffy[4] = -0.0000040;  fSigy[4] =  0.0000669;
  fOffx[5] =  0.0000138;  fSigx[5] =  0.0002442;
  fOffy[5] =  0.0000001;  fSigy[5] =  0.0000753;
  fOffx[6] =  0.0000059;  fSigx[6] =  0.0002075;
  fOffy[6] = -0.0000009;  fSigy[6] =  0.0000618;
  fOffx[7] =  0.0000125;  fSigx[7] =  0.0002822;
  fOffy[7] = -0.0000017;  fSigy[7] =  0.0000851;
  fOffx[8] =  0.0000020;  fSigx[8] =  0.0003823;
  fOffy[8] =  0.0000048;  fSigy[8] =  0.0000706;
  fOffx[9] =  0.000064;   fSigx[9] =  0.00029;
  fOffy[9] =  0.0000122;  fSigy[9] =  0.0000735;
  fOffx[10] = -0.0000020;  fSigx[10] =  0.0001917;
  fOffy[10] =  0.0000028;  fSigy[10] =  0.0000612;
  fOffx[11] =  0.0000073;  fSigx[11] =  0.0003012;
  fOffy[11] =  0.0000032;  fSigy[11] =  0.0000962;
  fOffx[12] =  0.0;  fSigx[12] =  0.0;
  fOffy[12] =  0.0;  fSigy[12] =  0.0;
  fOffx[13] =  0.0;  fSigx[13] =  0.0;
  fOffy[13] =  0.0;  fSigy[13] =  0.0;
  fOffx[14] =  0.0;  fSigx[14] =  0.0;
  fOffy[14] =  0.0;  fSigy[14] =  0.0;
  fOffx[15] = -0.00006;    fSigx[15] =  0.0005;
  fOffy[15] = -0.0000206;  fSigy[15] =  0.0001054;
  fOffx[16] =  0.0;  fSigx[16] =  0.0;
  fOffy[16] =  0.0;  fSigy[16] =  0.0;
  fOffx[17] = -0.0000051;  fSigx[17] =  0.0002168;
  fOffy[17] =  0.0000167;  fSigy[17] =  0.0000788;


  // extrapolated track cuts
  fOff_dx[0] =   0.0;          fSig_dx[0] = 0.008;
  fOff_dy[0] =   0.0;          fSig_dy[0] = 0.008;
  fOff_dax[0] =  0.0;          fSig_dax[0] = 0.00018;
  fOff_day[0] =  0.0;          fSig_day[0] = 0.00018;

  fOff_dx[1] =   0.0;          fSig_dx[1] = 0.008;
  fOff_dy[1] =   0.0;          fSig_dy[1] = 0.008;
  fOff_dax[1] =  0.0;          fSig_dax[1] = 0.00018;
  fOff_day[1] =  0.0;          fSig_day[1] = 0.00108;

  fOff_dx[2] =   0.0;           fSig_dx[2] = 0.008;
  fOff_dy[2] =   0.0;           fSig_dy[2] = 0.008;
  fOff_dax[2] =  0.0;           fSig_dax[2] = 0.00018;
  fOff_day[2] =  0.0;           fSig_day[2] = 0.00108;

  fOff_dx[3] =   0.0;           fSig_dx[3] = 0.008;
  fOff_dy[3] =   0.0;           fSig_dy[3] = 0.008;
  fOff_dax[3] =  0.0;           fSig_dax[3] = 0.00018;
  fOff_day[3] =  0.0;           fSig_day[3] = 0.00108;

  fAxCut = 0.046;

}

//_____________________________________________________________
void Na61AlVdArmParameters::SetupOffsetsAndSigmas_Jura_pPb2022(Int_t run_id)
{

  // dev cuts
  
  cout<<" Na61AlVdArmParameters::SetupOffsetsAndSigmas_Jura_pPb2022: setting offsets and sigmas for run: "<<run_id<<endl;

  // offsets and sigmas related to match is deviations   
  // this parameters are used in reconstruction of simulation data
  // see definitions of DevNames in Na61Module::Init()

fOffsetx[0] = -0.00061;    fSigmax[0] =  0.0067302;
fOffsety[0] = -0.00006;    fSigmay[0] =  0.0073977;
fOffsetx[1] = -0.00674;    fSigmax[1] =  0.0107603;
fOffsety[1] = 0.00322;    fSigmay[1] =  0.0105001;
fOffsetx[2] = -0.00302;    fSigmax[2] =  0.0091280;
fOffsety[2] = -0.00312;    fSigmay[2] =  0.0099346;
fOffsetx[3] = -0.00320;    fSigmax[3] =  0.0106574;
fOffsety[3] = -0.00032;    fSigmay[3] =  0.0114870;
fOffsetx[4] = -0.00284;    fSigmax[4] =  0.0098697;
fOffsety[4] = -0.00499;    fSigmay[4] =  0.0119499;
fOffsetx[5] = -0.00133;    fSigmax[5] =  0.0091687;
fOffsety[5] = 0.00008;    fSigmay[5] =  0.0096434;
fOffsetx[6] = -0.00210;    fSigmax[6] =  0.0096354;
fOffsety[6] = 0.00979;    fSigmay[6] =  0.0103535;
fOffsetx[7] = 0.00196;    fSigmax[7] =  0.0099129;
fOffsety[7] = -0.00407;    fSigmay[7] =  0.0102232;
fOffsetx[8] = -0.00010;    fSigmax[8] =  0.0108037;
fOffsety[8] = -0.00497;    fSigmay[8] =  0.0119970;
fOffsetx[9] = -0.00033;    fSigmax[9] =  0.0109378;
fOffsety[9] = 0.00068;    fSigmay[9] =  0.0117354;
fOffsetx[10] = -0.00143;    fSigmax[10] =  0.0101495;
fOffsety[10] = -0.00006;    fSigmay[10] =  0.0102334;
fOffsetx[11] = 0.00489;    fSigmax[11] =  0.0122557;
fOffsety[11] = 0.00233;    fSigmay[11] =  0.0115173;
fOffsetx[12] = 0.00129;    fSigmax[12] =  0.0112047;
fOffsety[12] = -0.00070;    fSigmay[12] =  0.0127531;
fOffsetx[13] = 0.00067;    fSigmax[13] =  0.0128065;
fOffsety[13] = -0.00473;    fSigmay[13] =  0.0128024;
fOffsetx[14] = 0.00028;    fSigmax[14] =  0.0128110;
fOffsety[14] = -0.00958;    fSigmay[14] =  0.0127404;
fOffsetx[15] = -0.00323;    fSigmax[15] =  0.0068851;
fOffsety[15] = -0.00092;    fSigmay[15] =  0.0076514;
fOffsetx[16] = 0.00289;    fSigmax[16] =  0.0085890;
fOffsety[16] = -0.00001;    fSigmay[16] =  0.0084863;
fOffsetx[17] = 0.00263;    fSigmax[17] =  0.0091430;
fOffsety[17] = -0.00366;    fSigmay[17] =  0.0105483;
fOffsetx[18] = 0.00321;    fSigmax[18] =  0.0092605;
fOffsety[18] = -0.00879;    fSigmay[18] =  0.0103721;
fOffsetx[19] = -0.01615;    fSigmax[19] =  0.0080757;
fOffsety[19] = 0.02464;    fSigmay[19] =  0.0084759;
fOffsetx[20] = -0.00425;    fSigmax[20] =  0.0072706;
fOffsety[20] = -0.01767;    fSigmay[20] =  0.0082406;
fOffsetx[21] = 0.00633;    fSigmax[21] =  0.0180825;
fOffsety[21] = -0.00036;    fSigmay[21] =  0.0194171;
fOffsetx[22] = -0.00121;    fSigmax[22] =  0.0155918;
fOffsety[22] = 0.02122;    fSigmay[22] =  0.0204143;
fOffsetx[23] = 0.07724;    fSigmax[23] =  0.0238577;
fOffsety[23] = -0.14947;    fSigmay[23] =  0.0325646;
fOffsetx[24] = -0.04474;    fSigmax[24] =  0.0213471;
fOffsety[24] = -0.00383;    fSigmay[24] =  0.0287481;
fOffsetx[25] = 0.00727;    fSigmax[25] =  0.0118295;
fOffsety[25] = -0.00012;    fSigmay[25] =  0.0111296;
fOffsetx[26] = 0.00916;    fSigmax[26] =  0.0164384;
fOffsety[26] = -0.00629;    fSigmay[26] =  0.0176689;
fOffsetx[27] = 0.00977;    fSigmax[27] =  0.0143596;
fOffsety[27] = 0.00891;    fSigmay[27] =  0.0153609;
fOffsetx[28] = -0.00494;    fSigmax[28] =  0.0212075;
fOffsety[28] = -0.01153;    fSigmay[28] =  0.0197675;
fOffsetx[29] = -0.00214;    fSigmax[29] =  0.0170820;
fOffsety[29] = 0.01406;    fSigmay[29] =  0.0213213;
fOffsetx[30] = 0.01353;    fSigmax[30] =  0.0283987;
fOffsety[30] = -0.14363;    fSigmay[30] =  0.0277767;
fOffsetx[31] = -0.04810;    fSigmax[31] =  0.0277129;
fOffsety[31] = -0.02910;    fSigmay[31] =  0.0305021;
fOffsetx[32] = 0.00046;    fSigmax[32] =  0.0111276;
fOffsety[32] = 0.00004;    fSigmay[32] =  0.0115985;
fOffsetx[33] = 0.00222;    fSigmax[33] =  0.0164753;
fOffsety[33] = 0.00321;    fSigmay[33] =  0.0179595;
fOffsetx[34] = 0.00205;    fSigmax[34] =  0.0174388;
fOffsety[34] = -0.00368;    fSigmay[34] =  0.0161512;
fOffsetx[35] = 0.00132;    fSigmax[35] =  0.0214001;
fOffsety[35] = 0.02619;    fSigmay[35] =  0.0211245;
fOffsetx[36] = 0.00256;    fSigmax[36] =  0.0185241;
fOffsety[36] = -0.01096;    fSigmay[36] =  0.0205508;
fOffsetx[37] = -0.02526;    fSigmax[37] =  0.0222440;
fOffsety[37] = -0.05995;    fSigmay[37] =  0.0239418;
fOffsetx[38] = -0.05265;    fSigmax[38] =  0.0211198;
fOffsety[38] = -0.10214;    fSigmay[38] =  0.0273995;
fOffsetx[39] = 0.00108;    fSigmax[39] =  0.0123555;
fOffsety[39] = -0.00020;    fSigmay[39] =  0.0123541;
fOffsetx[40] = 0.00105;    fSigmax[40] =  0.0142720;
fOffsety[40] = -0.00046;    fSigmay[40] =  0.0143683;
fOffsetx[41] = 0.00147;    fSigmax[41] =  0.0161543;
fOffsety[41] = 0.00615;    fSigmay[41] =  0.0169193;
fOffsetx[42] = 0.00094;    fSigmax[42] =  0.0177037;
fOffsety[42] = -0.00540;    fSigmay[42] =  0.0188059;
fOffsetx[43] = 0.00834;    fSigmax[43] =  0.0194563;
fOffsety[43] = 0.00087;    fSigmay[43] =  0.0200858;
fOffsetx[44] = 0.00729;    fSigmax[44] =  0.0195573;
fOffsety[44] = -0.12019;    fSigmay[44] =  0.0264865;
fOffsetx[45] = -0.00119;    fSigmax[45] =  0.0251066;
fOffsety[45] = -0.07308;    fSigmay[45] =  0.0228785;
fOffsetx[46] = 0.00632;    fSigmax[46] =  0.0122743;
fOffsety[46] = -0.00642;    fSigmay[46] =  0.0146689;
fOffsetx[47] = 0.00939;    fSigmax[47] =  0.0112065;
fOffsety[47] = 0.00535;    fSigmay[47] =  0.0134322;
fOffsetx[48] = 0.01560;    fSigmax[48] =  0.0141300;
fOffsety[48] = -0.00119;    fSigmay[48] =  0.0131996;
fOffsetx[49] = 0.00915;    fSigmax[49] =  0.0124050;
fOffsety[49] = 0.00449;    fSigmay[49] =  0.0123876;
fOffsetx[50] = -0.00342;    fSigmax[50] =  0.0129022;
fOffsety[50] = -0.01043;    fSigmay[50] =  0.0139941;
fOffsetx[51] = -0.00620;    fSigmax[51] =  0.0133156;
fOffsety[51] = 0.00534;    fSigmay[51] =  0.0130968;
fOffsetx[52] = -0.00126;    fSigmax[52] =  0.0134186;
fOffsety[52] = -0.00427;    fSigmay[52] =  0.0140513;
fOffsetx[53] = -0.00480;    fSigmax[53] =  0.0131771;
fOffsety[53] = -0.00082;    fSigmay[53] =  0.0138320;
fOffsetx[54] = 0.00401;    fSigmax[54] =  0.0120507;
fOffsety[54] = 0.00059;    fSigmay[54] =  0.0121625;
fOffsetx[55] = 0.00167;    fSigmax[55] =  0.0126791;
fOffsety[55] = -0.00138;    fSigmay[55] =  0.0119011;
fOffsetx[56] = -0.00524;    fSigmax[56] =  0.0112406;
fOffsety[56] = 0.00701;    fSigmay[56] =  0.0119354;
fOffsetx[57] = -0.00411;    fSigmax[57] =  0.0119507;
fOffsety[57] = -0.00807;    fSigmay[57] =  0.0126374;
fOffsetx[58] = 0.00493;    fSigmax[58] =  0.0115153;
fOffsety[58] = 0.01282;    fSigmay[58] =  0.0114676;
fOffsetx[59] = 0.00369;    fSigmax[59] =  0.0126307;
fOffsety[59] = -0.00900;    fSigmay[59] =  0.0128697;
  

  fOffx[0] =  0.0;  fSigx[0] =  0.0015; 
  fOffy[0] =  0.0;  fSigy[0] =  0.0015;

  fOffx[1] =  0.0;  fSigx[1] =  0.0015; 
  fOffy[1] =  0.0;  fSigy[1] =  0.0015;

  fOffx[2] =  0.0;  fSigx[2] =  0.0015; 
  fOffy[2] =  0.0;  fSigy[2] =  0.0015;

  fOffx[3] =  0.0;  fSigx[3] =  0.0015; 
  fOffy[3] =  0.0;  fSigy[3] =  0.0015;

  // cuts on matching with primary vertex
  fOffx[4] =  0.0000002;  fSigx[4] =  0.0003999;
  fOffy[4] = -0.0000040;  fSigy[4] =  0.0000669;
  fOffx[5] =  0.0000138;  fSigx[5] =  0.0002442;
  fOffy[5] =  0.0000001;  fSigy[5] =  0.0000753;
  fOffx[6] =  0.0000059;  fSigx[6] =  0.0002075;
  fOffy[6] = -0.0000009;  fSigy[6] =  0.0000618;
  fOffx[7] =  0.0000125;  fSigx[7] =  0.0002822;
  fOffy[7] = -0.0000017;  fSigy[7] =  0.0000851;
  fOffx[8] =  0.0000020;  fSigx[8] =  0.0003823;
  fOffy[8] =  0.0000048;  fSigy[8] =  0.0000706;
  fOffx[9] =  0.000064;   fSigx[9] =  0.00029;
  fOffy[9] =  0.0000122;  fSigy[9] =  0.0000735;
  fOffx[10] = -0.0000020;  fSigx[10] =  0.0001917;
  fOffy[10] =  0.0000028;  fSigy[10] =  0.0000612;
  fOffx[11] =  0.0000073;  fSigx[11] =  0.0003012;
  fOffy[11] =  0.0000032;  fSigy[11] =  0.0000962;
  fOffx[12] =  0.0;  fSigx[12] =  0.0;
  fOffy[12] =  0.0;  fSigy[12] =  0.0;
  fOffx[13] =  0.0;  fSigx[13] =  0.0;
  fOffy[13] =  0.0;  fSigy[13] =  0.0;
  fOffx[14] =  0.0;  fSigx[14] =  0.0;
  fOffy[14] =  0.0;  fSigy[14] =  0.0;
  fOffx[15] = -0.00006;    fSigx[15] =  0.0005;
  fOffy[15] = -0.0000206;  fSigy[15] =  0.0001054;
  fOffx[16] =  0.0;  fSigx[16] =  0.0;
  fOffy[16] =  0.0;  fSigy[16] =  0.0;
  fOffx[17] = -0.0000051;  fSigx[17] =  0.0002168;
  fOffy[17] =  0.0000167;  fSigy[17] =  0.0000788;


  // extrapolated track cuts
  fOff_dx[0] =   0.0;          fSig_dx[0] = 0.008;
  fOff_dy[0] =   0.0;          fSig_dy[0] = 0.008;
  fOff_dax[0] =  0.0;          fSig_dax[0] = 0.00018;
  fOff_day[0] =  0.0;          fSig_day[0] = 0.00018;

  fOff_dx[1] =   0.0;          fSig_dx[1] = 0.008;
  fOff_dy[1] =   0.0;          fSig_dy[1] = 0.008;
  fOff_dax[1] =  0.0;          fSig_dax[1] = 0.00018;
  fOff_day[1] =  0.0;          fSig_day[1] = 0.00108;

  fOff_dx[2] =   0.0;           fSig_dx[2] = 0.008;
  fOff_dy[2] =   0.0;           fSig_dy[2] = 0.008;
  fOff_dax[2] =  0.0;           fSig_dax[2] = 0.00018;
  fOff_day[2] =  0.0;           fSig_day[2] = 0.00108;

  fOff_dx[3] =   0.0;           fSig_dx[3] = 0.008;
  fOff_dy[3] =   0.0;           fSig_dy[3] = 0.008;
  fOff_dax[3] =  0.0;           fSig_dax[3] = 0.00018;
  fOff_day[3] =  0.0;           fSig_day[3] = 0.00108;

  fAxCut = 0.046;

}

//_____________________________________________________________
void Na61AlVdArmParameters::SetupOffsetsAndSigmas_Saleve_pPb2022(Int_t run_id)
{

  // dev cuts
  
  cout<<" Na61AlVdArmParameters::SetupOffsetsAndSigmas_Saleve_pPb2022: setting offsets and sigmas for run: "<<run_id<<endl;

  // offsets and sigmas related to match is deviations   
  // this parameters are used in reconstruction of simulation data
  // see definitions of DevNames in Na61Module::Init()

fOffsetx[0] = 0.00052;    fSigmax[0] =  0.0081478;
fOffsety[0] = 0.00194;    fSigmay[0] =  0.0087284;
fOffsetx[1] = -0.00127;    fSigmax[1] =  0.0117678;
fOffsety[1] = 0.00962;    fSigmay[1] =  0.0124911;
fOffsetx[2] = -0.00126;    fSigmax[2] =  0.0123523;
fOffsety[2] = -0.00403;    fSigmay[2] =  0.0122652;
fOffsetx[3] = 0.00048;    fSigmax[3] =  0.0129831;
fOffsety[3] = -0.00083;    fSigmay[3] =  0.0138286;
fOffsetx[4] = -0.00046;    fSigmax[4] =  0.0126441;
fOffsety[4] = 0.00147;    fSigmay[4] =  0.0136715;
fOffsetx[5] = -0.00302;    fSigmax[5] =  0.0115573;
fOffsety[5] = 0.00374;    fSigmay[5] =  0.0121640;
fOffsetx[6] = -0.00138;    fSigmax[6] =  0.0125311;
fOffsety[6] = 0.00210;    fSigmay[6] =  0.0124774;
fOffsetx[7] = -0.00030;    fSigmax[7] =  0.0126124;
fOffsety[7] = -0.00254;    fSigmay[7] =  0.0129103;
fOffsetx[8] = 0.00046;    fSigmax[8] =  0.0143662;
fOffsety[8] = 0.00080;    fSigmay[8] =  0.0155222;
fOffsetx[9] = 0.00014;    fSigmax[9] =  0.0134463;
fOffsety[9] = -0.00033;    fSigmay[9] =  0.0148247;
fOffsetx[10] = -0.00925;    fSigmax[10] =  0.0115264;
fOffsety[10] = -0.00386;    fSigmay[10] =  0.0118138;
fOffsetx[11] = -0.00434;    fSigmax[11] =  0.0147866;
fOffsety[11] = 0.00143;    fSigmay[11] =  0.0149313;
fOffsetx[12] = -0.00348;    fSigmax[12] =  0.0138537;
fOffsety[12] = -0.00994;    fSigmay[12] =  0.0153002;
fOffsetx[13] = -0.01253;    fSigmax[13] =  0.0159859;
fOffsety[13] = -0.00634;    fSigmay[13] =  0.0186210;
fOffsetx[14] = 0.00611;    fSigmax[14] =  0.0166263;
fOffsety[14] = 0.00248;    fSigmay[14] =  0.0184417;
fOffsetx[15] = -0.00778;    fSigmax[15] =  0.0073041;
fOffsety[15] = -0.00801;    fSigmay[15] =  0.0063414;
fOffsetx[16] = 0.00015;    fSigmax[16] =  0.0077808;
fOffsety[16] = -0.00193;    fSigmay[16] =  0.0081730;
fOffsetx[17] = 0.00085;    fSigmax[17] =  0.0100555;
fOffsety[17] = 0.01289;    fSigmay[17] =  0.0106754;
fOffsetx[18] = 0.00191;    fSigmax[18] =  0.0096942;
fOffsety[18] = -0.01203;    fSigmay[18] =  0.0101228;
fOffsetx[19] = -0.00230;    fSigmax[19] =  0.0095248;
fOffsety[19] = -0.00443;    fSigmay[19] =  0.0109051;
fOffsetx[20] = -0.00383;    fSigmax[20] =  0.0085403;
fOffsety[20] = 0.00532;    fSigmay[20] =  0.0105116;
fOffsetx[21] = -0.00128;    fSigmax[21] =  0.0189166;
fOffsety[21] = -0.00126;    fSigmay[21] =  0.0206265;
fOffsetx[22] = -0.00173;    fSigmax[22] =  0.0167143;
fOffsety[22] = -0.00073;    fSigmay[22] =  0.0246100;
fOffsetx[23] = -0.00795;    fSigmax[23] =  0.0231093;
fOffsety[23] = -0.03944;    fSigmay[23] =  0.0223868;
fOffsetx[24] = -0.03612;    fSigmax[24] =  0.0220605;
fOffsety[24] = 0.02232;    fSigmay[24] =  0.0386947;
fOffsetx[25] = -0.00482;    fSigmax[25] =  0.0114232;
fOffsety[25] = 0.00016;    fSigmay[25] =  0.0110900;
fOffsetx[26] = 0.00092;    fSigmax[26] =  0.0151986;
fOffsety[26] = -0.00073;    fSigmay[26] =  0.0163601;
fOffsetx[27] = 0.00271;    fSigmax[27] =  0.0150502;
fOffsety[27] = -0.00035;    fSigmay[27] =  0.0168292;
fOffsetx[28] = -0.00362;    fSigmax[28] =  0.0183752;
fOffsety[28] = 0.00144;    fSigmay[28] =  0.0220285;
fOffsetx[29] = 0.00231;    fSigmax[29] =  0.0164715;
fOffsety[29] = -0.00287;    fSigmay[29] =  0.0221297;
fOffsetx[30] = 0.01790;    fSigmax[30] =  0.0181321;
fOffsety[30] = -0.05047;    fSigmay[30] =  0.0206633;
fOffsetx[31] = -0.02958;    fSigmax[31] =  0.0211182;
fOffsety[31] = 0.00871;    fSigmay[31] =  0.0333584;
fOffsetx[32] = 0.00238;    fSigmax[32] =  0.0111659;
fOffsety[32] = -0.00182;    fSigmay[32] =  0.0117230;
fOffsetx[33] = 0.00179;    fSigmax[33] =  0.0166167;
fOffsety[33] = 0.00099;    fSigmay[33] =  0.0176649;
fOffsetx[34] = 0.00106;    fSigmax[34] =  0.0157870;
fOffsety[34] = 0.00205;    fSigmay[34] =  0.0177645;
fOffsetx[35] = 0.00093;    fSigmax[35] =  0.0192812;
fOffsety[35] = -0.02807;    fSigmay[35] =  0.0186383;
fOffsetx[36] = 0.00429;    fSigmax[36] =  0.0174071;
fOffsety[36] = 0.00845;    fSigmay[36] =  0.0220076;
fOffsetx[37] = -0.00702;    fSigmax[37] =  0.0207646;
fOffsety[37] = -0.04094;    fSigmay[37] =  0.0237881;
fOffsetx[38] = -0.05573;    fSigmax[38] =  0.0229989;
fOffsety[38] = 0.07998;    fSigmay[38] =  0.0303300;
fOffsetx[39] = 0.00479;    fSigmax[39] =  0.0120123;
fOffsety[39] = 0.00038;    fSigmay[39] =  0.0120830;
fOffsetx[40] = 0.00639;    fSigmax[40] =  0.0155377;
fOffsety[40] = 0.00596;    fSigmay[40] =  0.0155286;
fOffsetx[41] = 0.00044;    fSigmax[41] =  0.0163949;
fOffsety[41] = -0.00166;    fSigmay[41] =  0.0150694;
fOffsetx[42] = 0.00155;    fSigmax[42] =  0.0170398;
fOffsety[42] = -0.00982;    fSigmay[42] =  0.0201019;
fOffsetx[43] = -0.00972;    fSigmax[43] =  0.0234023;
fOffsety[43] = -0.00557;    fSigmay[43] =  0.0231473;
fOffsetx[44] = -0.02818;    fSigmax[44] =  0.0197605;
fOffsety[44] = -0.02112;    fSigmay[44] =  0.0271989;
fOffsetx[45] = -0.12566;    fSigmax[45] =  0.0357871;
fOffsety[45] = 0.06482;    fSigmay[45] =  0.0348985;
fOffsetx[46] = -0.00088;    fSigmax[46] =  0.0113449;
fOffsety[46] = -0.00467;    fSigmay[46] =  0.0120763;
fOffsetx[47] = 0.00013;    fSigmax[47] =  0.0118819;
fOffsety[47] = -0.00232;    fSigmay[47] =  0.0121023;
fOffsetx[48] = 0.00427;    fSigmax[48] =  0.0130047;
fOffsety[48] = -0.01007;    fSigmay[48] =  0.0124894;
fOffsetx[49] = 0.00354;    fSigmax[49] =  0.0146709;
fOffsety[49] = 0.00790;    fSigmay[49] =  0.0145554;
fOffsetx[50] = -0.00020;    fSigmax[50] =  0.0120631;
fOffsety[50] = 0.00522;    fSigmay[50] =  0.0137532;
fOffsetx[51] = -0.00357;    fSigmax[51] =  0.0125383;
fOffsety[51] = 0.00109;    fSigmay[51] =  0.0137208;
fOffsetx[52] = 0.00123;    fSigmax[52] =  0.0136812;
fOffsety[52] = 0.00115;    fSigmay[52] =  0.0129456;
fOffsetx[53] = -0.00008;    fSigmax[53] =  0.0135456;
fOffsety[53] = 0.00720;    fSigmay[53] =  0.0139243;
fOffsetx[54] = -0.00337;    fSigmax[54] =  0.0119288;
fOffsety[54] = 0.00111;    fSigmay[54] =  0.0135442;
fOffsetx[55] = -0.00566;    fSigmax[55] =  0.0118299;
fOffsety[55] = -0.00283;    fSigmay[55] =  0.0115203;
fOffsetx[56] = 0.00493;    fSigmax[56] =  0.0124754;
fOffsety[56] = 0.00005;    fSigmay[56] =  0.0130605;
fOffsetx[57] = 0.00401;    fSigmax[57] =  0.0118699;
fOffsety[57] = -0.00675;    fSigmay[57] =  0.0142717;
fOffsetx[58] = 0.00958;    fSigmax[58] =  0.0127854;
fOffsety[58] = -0.00576;    fSigmay[58] =  0.0133355;
fOffsetx[59] = 0.01117;    fSigmax[59] =  0.0138672;
fOffsety[59] = -0.00207;    fSigmay[59] =  0.0133792;


  fOffx[0] =  0.0;  fSigx[0] =  0.0015; 
  fOffy[0] =  0.0;  fSigy[0] =  0.0015;

  fOffx[1] =  0.0;  fSigx[1] =  0.0015; 
  fOffy[1] =  0.0;  fSigy[1] =  0.0015;

  fOffx[2] =  0.0;  fSigx[2] =  0.0015; 
  fOffy[2] =  0.0;  fSigy[2] =  0.0015;

  fOffx[3] =  0.0;  fSigx[3] =  0.0015; 
  fOffy[3] =  0.0;  fSigy[3] =  0.0015;

  // cuts on matching with primary vertex
  fOffx[4] =  0.0000002;  fSigx[4] =  0.0003999;
  fOffy[4] = -0.0000040;  fSigy[4] =  0.0000669;
  fOffx[5] =  0.0000138;  fSigx[5] =  0.0002442;
  fOffy[5] =  0.0000001;  fSigy[5] =  0.0000753;
  fOffx[6] =  0.0000059;  fSigx[6] =  0.0002075;
  fOffy[6] = -0.0000009;  fSigy[6] =  0.0000618;
  fOffx[7] =  0.0000125;  fSigx[7] =  0.0002822;
  fOffy[7] = -0.0000017;  fSigy[7] =  0.0000851;
  fOffx[8] =  0.0000020;  fSigx[8] =  0.0003823;
  fOffy[8] =  0.0000048;  fSigy[8] =  0.0000706;
  fOffx[9] =  0.000064;   fSigx[9] =  0.00029;
  fOffy[9] =  0.0000122;  fSigy[9] =  0.0000735;
  fOffx[10] = -0.0000020;  fSigx[10] =  0.0001917;
  fOffy[10] =  0.0000028;  fSigy[10] =  0.0000612;
  fOffx[11] =  0.0000073;  fSigx[11] =  0.0003012;
  fOffy[11] =  0.0000032;  fSigy[11] =  0.0000962;
  fOffx[12] =  0.0;  fSigx[12] =  0.0;
  fOffy[12] =  0.0;  fSigy[12] =  0.0;
  fOffx[13] =  0.0;  fSigx[13] =  0.0;
  fOffy[13] =  0.0;  fSigy[13] =  0.0;
  fOffx[14] =  0.0;  fSigx[14] =  0.0;
  fOffy[14] =  0.0;  fSigy[14] =  0.0;
  fOffx[15] = -0.00006;    fSigx[15] =  0.0005;
  fOffy[15] = -0.0000206;  fSigy[15] =  0.0001054;
  fOffx[16] =  0.0;  fSigx[16] =  0.0;
  fOffy[16] =  0.0;  fSigy[16] =  0.0;
  fOffx[17] = -0.0000051;  fSigx[17] =  0.0002168;
  fOffy[17] =  0.0000167;  fSigy[17] =  0.0000788;


  // extrapolated track cuts
  fOff_dx[0] =   0.0;          fSig_dx[0] = 0.008;
  fOff_dy[0] =   0.0;          fSig_dy[0] = 0.008;
  fOff_dax[0] =  0.0;          fSig_dax[0] = 0.00018;
  fOff_day[0] =  0.0;          fSig_day[0] = 0.00018;

  fOff_dx[1] =   0.0;          fSig_dx[1] = 0.008;
  fOff_dy[1] =   0.0;          fSig_dy[1] = 0.008;
  fOff_dax[1] =  0.0;          fSig_dax[1] = 0.00018;
  fOff_day[1] =  0.0;          fSig_day[1] = 0.00108;

  fOff_dx[2] =   0.0;           fSig_dx[2] = 0.008;
  fOff_dy[2] =   0.0;           fSig_dy[2] = 0.008;
  fOff_dax[2] =  0.0;           fSig_dax[2] = 0.00018;
  fOff_day[2] =  0.0;           fSig_day[2] = 0.00108;

  fOff_dx[3] =   0.0;           fSig_dx[3] = 0.008;
  fOff_dy[3] =   0.0;           fSig_dy[3] = 0.008;
  fOff_dax[3] =  0.0;           fSig_dax[3] = 0.00018;
  fOff_day[3] =  0.0;           fSig_day[3] = 0.00108;

  fAxCut = 0.046;

}

//_____________________________________________________________
void Na61AlVdArmParameters::SetupOffsetsAndSigmas_SaleveNominal(Int_t run_id)
{

  // dev cuts
  
  cout<<" Na61AlVdArmParameters::SetupOffsetsAndSigmas_SaleveNominal: setting offsets and sigmas for run: "<<run_id<<endl;

  // offsets and sigmas related to match is deviations   
  // this parameters are used in reconstruction of simulation data
  // see definitions of DevNames in Na61AlVdTrackingInitModule::Init()

fOffsetx[0] = 0.00000;    fSigmax[0] =  0.0081830;
fOffsety[0] = 0.00000;    fSigmay[0] =  0.0088311;
fOffsetx[1] = 0.00000;    fSigmax[1] =  0.0117829;
fOffsety[1] = 0.00000;    fSigmay[1] =  0.0126218;
fOffsetx[2] = 0.00000;    fSigmax[2] =  0.0120318;
fOffsety[2] = 0.00000;    fSigmay[2] =  0.0126434;
fOffsetx[3] = 0.00000;    fSigmax[3] =  0.0127788;
fOffsety[3] = 0.00000;    fSigmay[3] =  0.0135693;
fOffsetx[4] = 0.00000;    fSigmax[4] =  0.0121559;
fOffsety[4] = 0.00000;    fSigmay[4] =  0.0132953;
fOffsetx[5] = 0.00000;    fSigmax[5] =  0.0114347;
fOffsety[5] = 0.00000;    fSigmay[5] =  0.0121771;
fOffsetx[6] = 0.00000;    fSigmax[6] =  0.0120442;
fOffsety[6] = 0.00000;    fSigmay[6] =  0.0124017;
fOffsetx[7] = 0.00000;    fSigmax[7] =  0.0124819;
fOffsety[7] = 0.00000;    fSigmay[7] =  0.0129045;
fOffsetx[8] = 0.00000;    fSigmax[8] =  0.0145894;
fOffsety[8] = 0.00000;    fSigmay[8] =  0.0148319;
fOffsetx[9] = 0.00000;    fSigmax[9] =  0.0139158;
fOffsety[9] = 0.00000;    fSigmay[9] =  0.0147842;
fOffsetx[10] = 0.00000;    fSigmax[10] =  0.0114380;
fOffsety[10] = 0.00000;    fSigmay[10] =  0.0117170;
fOffsetx[11] = 0.00000;    fSigmax[11] =  0.0146778;
fOffsety[11] = 0.00000;    fSigmay[11] =  0.0145518;
fOffsetx[12] = 0.00000;    fSigmax[12] =  0.0136237;
fOffsety[12] = 0.00000;    fSigmay[12] =  0.0151663;
fOffsetx[13] = 0.00000;    fSigmax[13] =  0.0166559;
fOffsety[13] = 0.00000;    fSigmay[13] =  0.0193654;
fOffsetx[14] = 0.00000;    fSigmax[14] =  0.0163712;
fOffsety[14] = 0.00000;    fSigmay[14] =  0.0172706;
fOffsetx[15] = 0.00000;    fSigmax[15] =  0.0071570;
fOffsety[15] = 0.00000;    fSigmay[15] =  0.0063084;
fOffsetx[16] = 0.00000;    fSigmax[16] =  0.0078875;
fOffsety[16] = 0.00000;    fSigmay[16] =  0.0082534;
fOffsetx[17] = 0.00000;    fSigmax[17] =  0.0098661;
fOffsety[17] = 0.00000;    fSigmay[17] =  0.0103857;
fOffsetx[18] = 0.00000;    fSigmax[18] =  0.0097136;
fOffsety[18] = 0.00000;    fSigmay[18] =  0.0102804;
fOffsetx[19] = 0.00000;    fSigmax[19] =  0.0094875;
fOffsety[19] = 0.00000;    fSigmay[19] =  0.0110667;
fOffsetx[20] = 0.00000;    fSigmax[20] =  0.0087451;
fOffsety[20] = 0.00000;    fSigmay[20] =  0.0109275;
fOffsetx[21] = 0.00000;    fSigmax[21] =  0.0186680;
fOffsety[21] = 0.00000;    fSigmay[21] =  0.0196205;
fOffsetx[22] = 0.00000;    fSigmax[22] =  0.0165328;
fOffsety[22] = 0.00000;    fSigmay[22] =  0.0237489;
fOffsetx[23] = 0.00000;    fSigmax[23] =  0.0246564;
fOffsety[23] = 0.00000;    fSigmay[23] =  0.0239062;
fOffsetx[24] = 0.00000;    fSigmax[24] =  0.0236755;
fOffsety[24] = 0.00000;    fSigmay[24] =  0.0408588;
fOffsetx[25] = 0.00000;    fSigmax[25] =  0.0113653;
fOffsety[25] = 0.00000;    fSigmay[25] =  0.0113806;
fOffsetx[26] = 0.00000;    fSigmax[26] =  0.0148538;
fOffsety[26] = 0.00000;    fSigmay[26] =  0.0161898;
fOffsetx[27] = 0.00000;    fSigmax[27] =  0.0152769;
fOffsety[27] = 0.00000;    fSigmay[27] =  0.0167876;
fOffsetx[28] = 0.00000;    fSigmax[28] =  0.0183926;
fOffsety[28] = 0.00000;    fSigmay[28] =  0.0209386;
fOffsetx[29] = 0.00000;    fSigmax[29] =  0.0179723;
fOffsety[29] = 0.00000;    fSigmay[29] =  0.0204218;
fOffsetx[30] = 0.00000;    fSigmax[30] =  0.0191416;
fOffsety[30] = 0.00000;    fSigmay[30] =  0.0225640;
fOffsetx[31] = 0.00000;    fSigmax[31] =  0.0220689;
fOffsety[31] = 0.00000;    fSigmay[31] =  0.0320766;
fOffsetx[32] = 0.00000;    fSigmax[32] =  0.0111282;
fOffsety[32] = 0.00000;    fSigmay[32] =  0.0119685;
fOffsetx[33] = 0.00000;    fSigmax[33] =  0.0162771;
fOffsety[33] = 0.00000;    fSigmay[33] =  0.0176218;
fOffsetx[34] = 0.00000;    fSigmax[34] =  0.0151228;
fOffsety[34] = 0.00000;    fSigmay[34] =  0.0170053;
fOffsetx[35] = 0.00000;    fSigmax[35] =  0.0183065;
fOffsety[35] = 0.00000;    fSigmay[35] =  0.0188756;
fOffsetx[36] = 0.00000;    fSigmax[36] =  0.0167719;
fOffsety[36] = 0.00000;    fSigmay[36] =  0.0218830;
fOffsetx[37] = 0.00000;    fSigmax[37] =  0.0210523;
fOffsety[37] = 0.00000;    fSigmay[37] =  0.0248084;
fOffsetx[38] = 0.00000;    fSigmax[38] =  0.0224151;
fOffsety[38] = 0.00000;    fSigmay[38] =  0.0317123;
fOffsetx[39] = 0.00000;    fSigmax[39] =  0.0120219;
fOffsety[39] = 0.00000;    fSigmay[39] =  0.0119852;
fOffsetx[40] = 0.00000;    fSigmax[40] =  0.0148893;
fOffsety[40] = 0.00000;    fSigmay[40] =  0.0162337;
fOffsetx[41] = 0.00000;    fSigmax[41] =  0.0151564;
fOffsety[41] = 0.00000;    fSigmay[41] =  0.0146136;
fOffsetx[42] = 0.00000;    fSigmax[42] =  0.0166262;
fOffsety[42] = 0.00000;    fSigmay[42] =  0.0192276;
fOffsetx[43] = 0.00000;    fSigmax[43] =  0.0229795;
fOffsety[43] = 0.00000;    fSigmay[43] =  0.0229538;
fOffsetx[44] = 0.00000;    fSigmax[44] =  0.0192721;
fOffsety[44] = 0.00000;    fSigmay[44] =  0.0278967;
fOffsetx[45] = 0.00000;    fSigmax[45] =  0.0298885;
fOffsety[45] = 0.00000;    fSigmay[45] =  0.0352376;
fOffsetx[46] = 0.00000;    fSigmax[46] =  0.0118183;
fOffsety[46] = 0.00000;    fSigmay[46] =  0.0122505;
fOffsetx[47] = 0.00000;    fSigmax[47] =  0.0118089;
fOffsety[47] = 0.00000;    fSigmay[47] =  0.0121104;
fOffsetx[48] = 0.00000;    fSigmax[48] =  0.0132968;
fOffsety[48] = 0.00000;    fSigmay[48] =  0.0128608;
fOffsetx[49] = 0.00000;    fSigmax[49] =  0.0143239;
fOffsety[49] = 0.00000;    fSigmay[49] =  0.0142901;
fOffsetx[50] = 0.00000;    fSigmax[50] =  0.0118276;
fOffsety[50] = 0.00000;    fSigmay[50] =  0.0136148;
fOffsetx[51] = 0.00000;    fSigmax[51] =  0.0127312;
fOffsety[51] = 0.00000;    fSigmay[51] =  0.0134487;
fOffsetx[52] = 0.00000;    fSigmax[52] =  0.0130201;
fOffsety[52] = 0.00000;    fSigmay[52] =  0.0135285;
fOffsetx[53] = 0.00000;    fSigmax[53] =  0.0137123;
fOffsety[53] = 0.00000;    fSigmay[53] =  0.0138798;
fOffsetx[54] = 0.00000;    fSigmax[54] =  0.0114311;
fOffsety[54] = 0.00000;    fSigmay[54] =  0.0127794;
fOffsetx[55] = 0.00000;    fSigmax[55] =  0.0122366;
fOffsety[55] = 0.00000;    fSigmay[55] =  0.0115744;
fOffsetx[56] = 0.00000;    fSigmax[56] =  0.0119618;
fOffsety[56] = 0.00000;    fSigmay[56] =  0.0122811;
fOffsetx[57] = 0.00000;    fSigmax[57] =  0.0125846;
fOffsety[57] = 0.00000;    fSigmay[57] =  0.0128048;
fOffsetx[58] = 0.00000;    fSigmax[58] =  0.0132080;
fOffsety[58] = 0.00000;    fSigmay[58] =  0.0135849;
fOffsetx[59] = 0.00000;    fSigmax[59] =  0.0133734;
fOffsety[59] = 0.00000;    fSigmay[59] =  0.0126713;


  fOffx[0] =  0.0;  fSigx[0] =  0.0015; 
  fOffy[0] =  0.0;  fSigy[0] =  0.0015;

  fOffx[1] =  0.0;  fSigx[1] =  0.0015; 
  fOffy[1] =  0.0;  fSigy[1] =  0.0015;

  fOffx[2] =  0.0;  fSigx[2] =  0.0015; 
  fOffy[2] =  0.0;  fSigy[2] =  0.0015;

  fOffx[3] =  0.0;  fSigx[3] =  0.0015; 
  fOffy[3] =  0.0;  fSigy[3] =  0.0015;

  // cuts on matching with primary vertex
  fOffx[4] =  0.0000002;  fSigx[4] =  0.0003999;
  fOffy[4] = -0.0000040;  fSigy[4] =  0.0000669;
  fOffx[5] =  0.0000138;  fSigx[5] =  0.0002442;
  fOffy[5] =  0.0000001;  fSigy[5] =  0.0000753;
  fOffx[6] =  0.0000059;  fSigx[6] =  0.0002075;
  fOffy[6] = -0.0000009;  fSigy[6] =  0.0000618;
  fOffx[7] =  0.0000125;  fSigx[7] =  0.0002822;
  fOffy[7] = -0.0000017;  fSigy[7] =  0.0000851;
  fOffx[8] =  0.0000020;  fSigx[8] =  0.0003823;
  fOffy[8] =  0.0000048;  fSigy[8] =  0.0000706;
  fOffx[9] =  0.000064;   fSigx[9] =  0.00029;
  fOffy[9] =  0.0000122;  fSigy[9] =  0.0000735;
  fOffx[10] = -0.0000020;  fSigx[10] =  0.0001917;
  fOffy[10] =  0.0000028;  fSigy[10] =  0.0000612;
  fOffx[11] =  0.0000073;  fSigx[11] =  0.0003012;
  fOffy[11] =  0.0000032;  fSigy[11] =  0.0000962;
  fOffx[12] =  0.0;  fSigx[12] =  0.0;
  fOffy[12] =  0.0;  fSigy[12] =  0.0;
  fOffx[13] =  0.0;  fSigx[13] =  0.0;
  fOffy[13] =  0.0;  fSigy[13] =  0.0;
  fOffx[14] =  0.0;  fSigx[14] =  0.0;
  fOffy[14] =  0.0;  fSigy[14] =  0.0;
  fOffx[15] = -0.00006;    fSigx[15] =  0.0005;
  fOffy[15] = -0.0000206;  fSigy[15] =  0.0001054;
  fOffx[16] =  0.0;  fSigx[16] =  0.0;
  fOffy[16] =  0.0;  fSigy[16] =  0.0;
  fOffx[17] = -0.0000051;  fSigx[17] =  0.0002168;
  fOffy[17] =  0.0000167;  fSigy[17] =  0.0000788;


  // extrapolated track cuts
  fOff_dx[0] =   0.0;          fSig_dx[0] = 0.008;
  fOff_dy[0] =   0.0;          fSig_dy[0] = 0.008;
  fOff_dax[0] =  0.0;          fSig_dax[0] = 0.00018;
  fOff_day[0] =  0.0;          fSig_day[0] = 0.00018;

  fOff_dx[1] =   0.0;          fSig_dx[1] = 0.008;
  fOff_dy[1] =   0.0;          fSig_dy[1] = 0.008;
  fOff_dax[1] =  0.0;          fSig_dax[1] = 0.00018;
  fOff_day[1] =  0.0;          fSig_day[1] = 0.00108;

  fOff_dx[2] =   0.0;           fSig_dx[2] = 0.008;
  fOff_dy[2] =   0.0;           fSig_dy[2] = 0.008;
  fOff_dax[2] =  0.0;           fSig_dax[2] = 0.00018;
  fOff_day[2] =  0.0;           fSig_day[2] = 0.00108;

  fOff_dx[3] =   0.0;           fSig_dx[3] = 0.008;
  fOff_dy[3] =   0.0;           fSig_dy[3] = 0.008;
  fOff_dax[3] =  0.0;           fSig_dax[3] = 0.00018;
  fOff_day[3] =  0.0;           fSig_day[3] = 0.00108;

  fAxCut = 0.046;

}

//_____________________________________________________________________
void Na61AlVdArmParameters::LocalToGlobal(Int_t ii, TObjArray* hits)
{

  // ii - sensor index (0-33)


  Float_t alpha   = GetRotZ(ii);
  Float_t beta    = GetRotY(ii);
  Float_t gamma   = GetRotX(ii);
  Float_t VolumeX = GetVolumeX(ii);
  Float_t VolumeY = GetVolumeY(ii);
  Float_t VolumeZ = GetVolumeZ(ii);

  Double_t sa = TMath::Sin(alpha); 
  Double_t ca = TMath::Cos(alpha); 
  Double_t sb = TMath::Sin(beta); 
  Double_t cb = TMath::Cos(beta);
  Double_t sg = TMath::Sin(gamma);
  Double_t cg = TMath::Cos(gamma);

  //cout<<"ii="<<ii<<" "<<VolumeX<<" "<<VolumeY<<" "<<VolumeZ<<endl; 

  for(Int_t i=0;i<hits->GetEntries();i++){
    USensorHit* hit = (USensorHit*)hits->At(i);
    Float_t x1 = hit->GetLocalX();
    Float_t y1 = hit->GetLocalY();
    Float_t z1 = 0;

    Float_t x2; 
    Float_t y2; 
    Float_t z2;
    
    //ZXY
    x2 =  (ca*cb-sa*sg*sb) * x1   +   sa*cg * y1  +    (ca*sb+sa*sg*cb) * z1;
    y2 = -(sa*cb+ca*sg*sb) * x1   +   ca*cg * y1  +   (-sa*sb+ca*sg*cb) * z1;
    z2 =            -cg*sb * x1   -      sg * y1  +              cg*cb  * z1;
   
    hit->SetX(x2 + VolumeX);
    hit->SetY(y2 + VolumeY);
    hit->SetZ(z2 + VolumeZ);

  } 
  
}


/*
//____________________________________________________________________
ostream& operator<< (ostream & os,Na61AlVdArmParameters *armpars)
{
  os << " Parameters for "
     <<armpars->GetArmName()<< " arm "
     <<endl;
  return os;
}
*/
