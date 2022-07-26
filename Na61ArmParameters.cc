//____________________________________________________________________
//
// Interfase to set and acces arm geometry and matching
// parameters

#include "Na61ArmParameters.h"
#include <utl/ErrorLogger.h>

#ifndef ROOT_TMath
#include "TMath.h"
#endif

#ifndef __IOSTREAM__
#include <iostream>
#endif

#include <sstream>

using std::endl;
using std::cout;

using namespace utl;
using namespace std;

//____________________________________________________________________
// ClassImp(Na61ArmParameters);

//____________________________________________________________________
Na61ArmParameters::Na61ArmParameters(TString armname) {
  // normal constructor
  fInit = false;
  fArmName = armname;
  fJura = false;
  if (armname.Contains("Jura")) fJura = true;
  fSensorId = 0;
}

//____________________________________________________________________
Na61ArmParameters::~Na61ArmParameters() {}

//____________________________________________________________________
void Na61ArmParameters::Init() {
  SetDrotX(0, 0);
  SetDrotY(0, 0);
  SetDrotZ(0, 0);
  SetVolumeDz(0, 0);

  fInit = true;
  if (fJura) {
    SetupSensorRotation_JuraDec2016_mv1(fRunId);  // used for physics reco of Dec2016 data
    SetupSensorGeometry_JuraDec2016_mv1(fRunId);
    // SetupSensorRotation_JuraDec2016(fRunId);
    // SetupSensorGeometry_JuraDec2016(fRunId);
	if (fRunId >= 37797) {							//PbPb150 2018
        SetupSensorRotation_JuraNov2018(fRunId);
		SetupSensorGeometry_JuraNov2018(fRunId);
        
	}	
		if (fRunId < 27469) {
		  SetupOffsetsAndSigmas_JuraDec2016_field(fRunId);
		} else if (fRunId < 27479) {
		  SetupOffsetsAndSigmas_JuraDec2016_nofield(fRunId);
		} else if (fRunId < 33396) {
		  SetupOffsetsAndSigmas_JuraOct2017_field(fRunId);  // XeLa150
		} else if (fRunId < 33408) {                          // should be refine
		  SetupOffsetsAndSigmas_JuraOct2017_nofield(fRunId);  // Xela150 no field
		} else if (fRunId < 33800) {
		  SetupOffsetsAndSigmas_JuraOct2017_field_xela75(fRunId);  // XeLa75
		} else if (fRunId < 37797) {
		  SetupOffsetsAndSigmas_JuraOct2017_field_xela40(fRunId);  // XeLa40 and further
		} else {
		  SetupOffsetsAndSigmas_JuraNov2018_field(fRunId);  // PbPb 2018
		}
    CellurarAutomatonParams_Jura_XeLa(fRunId);

    // if(fRunId>182)SetupOffsetsAndSigmas_JuraDec2016_nofield(fRunId);
    // else SetupOffsetsAndSigmas_JuraDec2016_field(fRunId);

  } else {
    SetupSensorRotation_SaleveDec2016_mv1(fRunId);  // used for physics reco of Dec2016 data
    SetupSensorGeometry_SaleveDec2016_mv1(fRunId);
    // SetupSensorRotation_SaleveDec2016(fRunId);
    // SetupSensorGeometry_SaleveDec2016(fRunId);
	if (fRunId >= 37797) {							//PbPb150 2018
		SetupSensorRotation_SaleveNov2018(fRunId);
        SetupSensorGeometry_SaleveNov2018(fRunId);
	}	

		if (fRunId < 27469) {
		  SetupOffsetsAndSigmas_SaleveDec2016_field(fRunId);
		} else if (fRunId < 27497) {
		  SetupOffsetsAndSigmas_SaleveDec2016_nofield(fRunId);
		} else if (fRunId < 33396) {
		  SetupOffsetsAndSigmas_SaleveOct2017_field(fRunId);  // XeLa150
		} else if (fRunId < 33408) {
		  SetupOffsetsAndSigmas_SaleveOct2017_nofield(fRunId);  // Xela150 no field
		 } else if (fRunId < 33800) {
		  SetupOffsetsAndSigmas_SaleveOct2017_field_xela75(fRunId);  // XeLa75 
		 } else if (fRunId < 37797) {
		  SetupOffsetsAndSigmas_SaleveOct2017_field_xela40(fRunId);  // XeLa40 and further
		} else {
		  SetupOffsetsAndSigmas_SaleveNov2018_field(fRunId);  // PbPb2018
		}
    CellurarAutomatonParams_Saleve_XeLa(fRunId);

    // if(fRunId>182) SetupOffsetsAndSigmas_SaleveDec2016_nofield(fRunId);
    // else SetupOffsetsAndSigmas_SaleveDec2016_field(fRunId);
  }
}

//_____________________________________________________________
void Na61ArmParameters::TransformVolumesToGlobal(double dx, double dy, double dz, double rotx, double roty, double rotz) {
    ostringstream info;
    info << "Na61ArmParameters::TransformVolumesToGlobal: Transforming Volumes, arm offx= " << dx << "  rotx=" << rotx << " roty=" << roty << " rotz=" << rotz;
    INFO(info);

  double sa = TMath::Sin(rotz);
  double ca = TMath::Cos(rotz);
  double sb = TMath::Sin(roty);
  double cb = TMath::Cos(roty);
  double sg = TMath::Sin(rotx);
  double cg = TMath::Cos(rotx);

  for (int i = 0; i < 8; i++) {
    double x0 = GetVolumeX(i);
    double y0 = GetVolumeY(i);
    double z0 = GetVolumeZ(i);

    z0 = z0 - 75.0;

    // double x1 =      cb * x0                +     sb * z0;
    // double y1 =  -sg*sb * x0   +   cg * y0  +  sg*cb * z0;
    // double z1 =  -cg*sb * x0   -   sg * y0  +  cg*cb * z0;

    // ZXY
    double x1 = (ca * cb - sa * sg * sb) * x0 + sa * cg * y0 + (ca * sb + sa * sg * cb) * z0;
    double y1 = -(sa * cb + ca * sg * sb) * x0 + ca * cg * y0 + (-sa * sb + ca * sg * cb) * z0;
    double z1 = -cg * sb * x0 - sg * y0 + cg * cb * z0;

    z1 = z1 + 75.0;

    // cout<<"z0="<<z0<<" z1="<<z1<<endl;

    SetVolumeX(i, x1 + dx);
    SetVolumeY(i, y1 + dy);
    SetVolumeZ(i, z1 + dz);
    SetRotX(i, GetRotX(i) + rotx);
    SetRotY(i, GetRotY(i) + roty);
  }
    ostringstream info_fin;
    info_fin << "after transformation volumeX_Vds1_0=" << fVolumeX[0] << "   volumeX_Vds3_0=" << fVolumeX[2] << "  volumeX_Vds3_1=" << fVolumeX[3];
    INFO(info_fin);
}

//_____________________________________________________________
void Na61ArmParameters::SetVolumesInGlobal(double dx, double dy, double dz, double rotx, double roty, double rotz) {
    ostringstream info;
    info << "Na61ArmParameters::SetVolumesInGlobal: Setting Volumes, arm offx=" << dx << "  rotx=" << rotx << " roty=" << roty << " rotz=" << rotz;
    INFO(info);
    
  double sa = TMath::Sin(rotz);
  double ca = TMath::Cos(rotz);
  double sb = TMath::Sin(roty);
  double cb = TMath::Cos(roty);
  double sg = TMath::Sin(rotx);
  double cg = TMath::Cos(rotx);

  for (int i = 0; i < 8; i++) {
    double x0 = GetVolumeX(i);
    double y0 = GetVolumeY(i);
    double z0 = GetVolumeZ(i);

    z0 = z0 - 75.0;

    // ZXY
    double x1 = (ca * cb - sa * sg * sb) * x0 + sa * cg * y0 + (ca * sb + sa * sg * cb) * z0;
    double y1 = -(sa * cb + ca * sg * sb) * x0 + ca * cg * y0 + (-sa * sb + ca * sg * cb) * z0;
    double z1 = -cg * sb * x0 - sg * y0 + cg * cb * z0;

    z1 = z1 + 75.0;

    // cout<<"z0="<<z0<<" z1="<<z1<<endl;

    SetVolumeGlobX(i, x1 + dx);
    SetVolumeGlobY(i, y1 + dy);
    SetVolumeGlobZ(i, z1 + dz);
    // SetVolumeGlobX(i,x0+dx);
    // SetVolumeGlobY(i,y0+dy);
    // SetVolumeGlobZ(i,z0+dz);
    SetRotGlobX(i, GetRotX(i) + rotx);
    SetRotGlobY(i, GetRotY(i) + roty);
  }
    ostringstream info_fin;
    info_fin << "Volumes in global frame: volumeX_Vds1_0=" << fVolumeGlobX[0] << "   volumeX_Vds3_0=" << fVolumeGlobX[2] << "  volumeX_Vds3_1=" << fVolumeGlobX[3];
    INFO(info_fin);
}

//___________________________________________________________________________________________
void Na61ArmParameters::CellurarAutomatonParams_Jura(int /*run_id*/) {
    ostringstream info;
    info << " Na61ArmParameters::CellurarAutomatonParams_Jura: default";
    INFO(info);
    
  // Vds1 (primary tracks)
  double meanLocX_Vds1[8] = {-0.20595, -0.30500, -0.24483, -0.24645, -0.24752, -0.25225, 0.00000, -0.25558};
  double sigmaLocX_Vds1[8] = {0.04000, 0.04000, 0.01997, 0.02108, 0.01752, 0.01831, 0.00000, 0.02302};
  double meanLocY_Vds1[8] = {0.11689, 0.14484, 0.13275, 0.12756, 0.14250, 0.13666, 0.00000, 0.13788};
  double sigmaLocY_Vds1[8] = {0.02912, 0.02121, 0.01135, 0.01277, 0.01233, 0.01799, 0.00000, 0.02590};

  double meanX_Vds1[8] = {-0.00571, -0.01610, -0.00995, -0.01016, -0.01047, -0.01090, 0.00000, -0.01095};
  double sigmaX_Vds1[8] = {0.00400, 0.00398, 0.00215, 0.00233, 0.00169, 0.00173, 0.00000, 0.00240};
  double meanY_Vds1[8] = {0.00075, 0.00371, 0.00198, 0.00226, 0.00202, 0.00231, 0.00000, 0.00230};
  double sigmaY_Vds1[8] = {0.00345, 0.00256, 0.00120, 0.00128, 0.00120, 0.00178, 0.00000, 0.00238};

  // Vds2 (primary tracks)
  double meanLocX_Vds2[8] = {0.00000, -0.21659, -0.27899, -0.27847, -0.24159, -0.24187, 0.00000, -0.24111};
  double sigmaLocX_Vds2[8] = {0.00000, 0.04000, 0.02494, 0.02450, 0.01800, 0.01938, 0.00000, 0.01326};
  double meanLocY_Vds2[8] = {0.00000, 0.12577, 0.14028, 0.13934, 0.14257, 0.13799, 0.00000, 0.14495};
  double sigmaLocY_Vds2[8] = {0.00000, 0.03741, 0.01450, 0.01484, 0.01286, 0.01530, 0.00000, 0.01215};

  double meanX_Vds2[8] = {0.00000, -0.00734, -0.01348, -0.01344, -0.01014, -0.01015, 0.00000, -0.00993};
  double sigmaX_Vds2[8] = {0.00000, 0.00400, 0.00241, 0.00239, 0.00175, 0.00190, 0.00000, 0.00131};
  double meanY_Vds2[8] = {0.00000, 0.00159, 0.00276, 0.00347, 0.00211, 0.00246, 0.00000, 0.00267};
  double sigmaY_Vds2[8] = {0.00000, 0.00400, 0.00142, 0.00144, 0.00127, 0.00151, 0.00000, 0.00135};

  // Vds3 (secondary tracks)
  double meanLocX_Vds3[8] = {-0.23950, -0.21564, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds3[8] = {0.01423, 0.02227, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanLocY_Vds3[8] = {0.13073, 0.12847, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds3[8] = {0.01325, 0.02254, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  double meanX_Vds3[8] = {-0.00950, -0.00719, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds3[8] = {0.00142, 0.00226, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanY_Vds3[8] = {0.00190, 0.00144, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds3[8] = {0.00118, 0.00213, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  // Vds4 (secondary tracks)
  double meanLocX_Vds4[8] = {-0.24212, -0.24026, -0.21653, -0.21475, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds4[8] = {0.01179, 0.01195, 0.02450, 0.01485, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanLocY_Vds4[8] = {0.13311, 0.13268, 0.13420, 0.09889, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds4[8] = {0.00982, 0.01215, 0.01832, 0.02010, 0.00000, 0.00000, 0.00000, 0.00000};

  double meanX_Vds4[8] = {-0.00973, -0.00962, -0.00716, -0.00731, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds4[8] = {0.00114, 0.00128, 0.00248, 0.00195, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanY_Vds4[8] = {0.00203, 0.00185, 0.00216, -0.00036, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds4[8] = {0.00091, 0.00124, 0.00170, 0.00187, 0.00000, 0.00000, 0.00000, 0.00000};

  SetCAMembers(meanLocX_Vds1, sigmaLocX_Vds1, meanLocY_Vds1, sigmaLocY_Vds1, meanLocX_Vds2, sigmaLocX_Vds2, meanLocY_Vds2, sigmaLocY_Vds2, meanLocX_Vds3, sigmaLocX_Vds3, meanLocY_Vds3, sigmaLocY_Vds3, meanLocX_Vds4, sigmaLocX_Vds4, meanLocY_Vds4, sigmaLocY_Vds4, meanX_Vds1, sigmaX_Vds1, meanY_Vds1, sigmaY_Vds1, meanX_Vds2, sigmaX_Vds2, meanY_Vds2, sigmaY_Vds2, meanX_Vds3, sigmaX_Vds3, meanY_Vds3, sigmaY_Vds3, meanX_Vds4, sigmaX_Vds4, meanY_Vds4, sigmaY_Vds4);
}
//___________________________________________________________________________________________
void Na61ArmParameters::CellurarAutomatonParams_Jura_XeLa(int /*run_id*/) {
    ostringstream info;
    info << " Na61ArmParameters::CellurarAutomatonParams_Jura: XeLa";
    INFO(info);
    
  // Vds1 (primary tracks)
  double meanLocX_Vds1[8] = {-0.09752, -0.01796, 0.00022, 0.00175, 0.01300, -0.00500, 0.00700, 0.00139};
  double sigmaLocX_Vds1[8] = {0.06731, 0.04964, 0.03304, 0.03419, 0.06311, 0.07189, 0.07812, 0.07618};
  double meanLocY_Vds1[8] = {0.02830, 0.01365, 0.00066, 0.00144, 0.00070, -0.00306, -0.00172, 0.00629};
  double sigmaLocY_Vds1[8] = {0.02752, 0.02935, 0.01193, 0.01262, 0.01281, 0.01511, 0.02584, 0.02122};

  double meanX_Vds1[8] = {-0.00969, -0.00175, 0.00005, 0.00035, -0.00070, -0.00004, 0.00014, 0.00054};
  double sigmaX_Vds1[8] = {0.00666, 0.00498, 0.00330, 0.00346, 0.00627, 0.00707, 0.00779, 0.00771};
  double meanY_Vds1[8] = {0.00278, 0.00126, 0.00006, 0.00020, 0.00006, -0.00011, 0.00006, 0.00020};
  double sigmaY_Vds1[8] = {0.00278, 0.00294, 0.00119, 0.00124, 0.00123, 0.00151, 0.00168, 0.00180};

  // Vds2 (primary tracks)
  double meanLocX_Vds2[8] = {0.00000, -0.06820, -0.00367, -0.00334, -0.00157, -0.00342, 0.00069, -0.00292};
  double sigmaLocX_Vds2[8] = {0.00000, 0.08000, 0.03709, 0.03794, 0.03395, 0.03205, 0.02944, 0.03128};
  double meanLocY_Vds2[8] = {0.00000, 0.04867, 0.00152, 0.00570, 0.00059, 0.00118, -0.00529, 0.00980};
  double sigmaLocY_Vds2[8] = {0.00000, 0.04971, 0.01247, 0.01954, 0.01250, 0.01430, 0.01676, 0.01632};

  double meanX_Vds2[8] = {0.00000, -0.00658, -0.00035, -0.00021, -0.00004, 0.00007, 0.00012, 0.00010};
  double sigmaX_Vds2[8] = {0.00000, 0.00884, 0.00372, 0.00381, 0.00337, 0.00318, 0.00298, 0.00294};
  double meanY_Vds2[8] = {0.00000, 0.00430, 0.00015, 0.00062, 0.00007, 0.00026, 0.00002, 0.00032};
  double sigmaY_Vds2[8] = {0.00000, 0.00487, 0.00125, 0.00197, 0.00117, 0.00138, 0.00115, 0.00130};

  // Vds3 (secondary tracks)
  double meanLocX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanLocY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  double meanX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  // Vds4 (secondary tracks)
  double meanLocX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanLocY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  double meanX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  SetCAMembers(meanLocX_Vds1, sigmaLocX_Vds1, meanLocY_Vds1, sigmaLocY_Vds1, meanLocX_Vds2, sigmaLocX_Vds2, meanLocY_Vds2, sigmaLocY_Vds2, meanLocX_Vds3, sigmaLocX_Vds3, meanLocY_Vds3, sigmaLocY_Vds3, meanLocX_Vds4, sigmaLocX_Vds4, meanLocY_Vds4, sigmaLocY_Vds4, meanX_Vds1, sigmaX_Vds1, meanY_Vds1, sigmaY_Vds1, meanX_Vds2, sigmaX_Vds2, meanY_Vds2, sigmaY_Vds2, meanX_Vds3, sigmaX_Vds3, meanY_Vds3, sigmaY_Vds3, meanX_Vds4, sigmaX_Vds4, meanY_Vds4, sigmaY_Vds4);
}

//___________________________________________________________________________________________
void Na61ArmParameters::CellurarAutomatonParams_Saleve(int /*run_id*/) {
    ostringstream info;
    info << " Na61ArmParameters::CellurarAutomatonParams_Saleve: default";
    INFO(info);
  // Vds1 (primary tracks)
  double meanLocX_Vds1[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds1[8] = {0.04000, 0.04000, 0.01997, 0.02108, 0.01752, 0.01831, 0.00000, 0.02302};
  double meanLocY_Vds1[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds1[8] = {0.02912, 0.02121, 0.01135, 0.01277, 0.01233, 0.01799, 0.00000, 0.02590};

  double meanX_Vds1[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds1[8] = {0.00400, 0.00398, 0.00215, 0.00233, 0.00169, 0.00173, 0.00000, 0.00240};
  double meanY_Vds1[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds1[8] = {0.00075, 0.00371, 0.00198, 0.00226, 0.00202, 0.00231, 0.00000, 0.00230};

  // Vds2 (primary tracks)
  double meanLocX_Vds2[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds2[8] = {0.00000, 0.04000, 0.02494, 0.02450, 0.01800, 0.01938, 0.00000, 0.01326};
  double meanLocY_Vds2[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds2[8] = {0.00000, 0.03741, 0.01450, 0.01484, 0.01286, 0.01530, 0.00000, 0.01215};

  double meanX_Vds2[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds2[8] = {0.00000, 0.00400, 0.00241, 0.00239, 0.00175, 0.00190, 0.00000, 0.00131};
  double meanY_Vds2[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds2[8] = {0.00000, 0.00400, 0.00142, 0.00144, 0.00127, 0.00151, 0.00000, 0.00135};

  // Vds3 (secondary tracks)
  double meanLocX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds3[8] = {0.01423, 0.02227, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanLocY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds3[8] = {0.01325, 0.02254, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  double meanX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds3[8] = {0.00142, 0.00226, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds3[8] = {0.00118, 0.00213, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  // Vds4 (secondary tracks)
  double meanLocX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds4[8] = {0.01179, 0.01195, 0.02450, 0.01485, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanLocY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds4[8] = {0.00982, 0.01215, 0.01832, 0.02010, 0.00000, 0.00000, 0.00000, 0.00000};

  double meanX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds4[8] = {0.00114, 0.00128, 0.00248, 0.00195, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds4[8] = {0.00091, 0.00124, 0.00170, 0.00187, 0.00000, 0.00000, 0.00000, 0.00000};

  SetCAMembers(meanLocX_Vds1, sigmaLocX_Vds1, meanLocY_Vds1, sigmaLocY_Vds1, meanLocX_Vds2, sigmaLocX_Vds2, meanLocY_Vds2, sigmaLocY_Vds2, meanLocX_Vds3, sigmaLocX_Vds3, meanLocY_Vds3, sigmaLocY_Vds3, meanLocX_Vds4, sigmaLocX_Vds4, meanLocY_Vds4, sigmaLocY_Vds4, meanX_Vds1, sigmaX_Vds1, meanY_Vds1, sigmaY_Vds1, meanX_Vds2, sigmaX_Vds2, meanY_Vds2, sigmaY_Vds2, meanX_Vds3, sigmaX_Vds3, meanY_Vds3, sigmaY_Vds3, meanX_Vds4, sigmaX_Vds4, meanY_Vds4, sigmaY_Vds4);
}

//___________________________________________________________________________________________
void Na61ArmParameters::CellurarAutomatonParams_Saleve_XeLa(int /*run_id*/) {
    ostringstream info;
    info << " Na61ArmParameters::CellurarAutomatonParams_Saleve: XeLa";
    INFO(info);
   
  // Vds1 (primary tracks)
  double meanLocX_Vds1[8] = {0.03902, 0.00296, 0.00373, 0.00501, -0.00700, -0.00114, -0.00212, -0.01300};
  double sigmaLocX_Vds1[8] = {0.04970, 0.03579, 0.03301, 0.03115, 0.06000, 0.05676, 0.06000, 0.06000};
  double meanLocY_Vds1[8] = {0.01168, 0.00571, 0.00120, -0.00061, 0.00071, 0.00096, -0.00015, 0.00036};
  double sigmaLocY_Vds1[8] = {0.06000, 0.02533, 0.01181, 0.01295, 0.01452, 0.01491, 0.02056, 0.02519};

  double meanX_Vds1[8] = {0.00385, 0.00026, 0.00007, -0.00008, -0.00070, -0.00034, -0.00016, -0.00050};
  double sigmaX_Vds1[8] = {0.00523, 0.00359, 0.00331, 0.00300, 0.00600, 0.00585, 0.00600, 0.00600};
  double meanY_Vds1[8] = {0.00114, 0.00052, -0.00001, 0.00020, 0.00005, 0.00014, 0.00017, -0.00009};
  double sigmaY_Vds1[8] = {0.00600, 0.00253, 0.00118, 0.00125, 0.00144, 0.00145, 0.00191, 0.00238};

  // Vds2 (primary tracks)
  double meanLocX_Vds2[8] = {0.00000, 0.04617, 0.00398, 0.00408, -0.00099, 0.00025, -0.00073, -0.00139};
  double sigmaLocX_Vds2[8] = {0.00000, 0.06000, 0.03391, 0.03211, 0.03218, 0.03302, 0.03070, 0.02871};
  double meanLocY_Vds2[8] = {0.00000, 0.03502, -0.00012, 0.00335, -0.00048, 0.00152, -0.00184, 0.00178};
  double sigmaLocY_Vds2[8] = {0.00000, 0.06000, 0.01438, 0.01774, 0.01270, 0.01417, 0.01696, 0.01728};

  double meanX_Vds2[8] = {0.00000, 0.00460, 0.00017, -0.00003, -0.00010, -0.00022, -0.00005, -0.00020};
  double sigmaX_Vds2[8] = {0.00000, 0.00600, 0.00339, 0.00314, 0.00322, 0.00324, 0.00299, 0.00287};
  double meanY_Vds2[8] = {0.00000, 0.00410, -0.00021, 0.00063, -0.00005, 0.00025, -0.00009, 0.00011};
  double sigmaY_Vds2[8] = {0.00000, 0.00600, 0.00144, 0.00172, 0.00126, 0.00138, 0.00149, 0.00159};

  // Vds3 (secondary tracks)
  double meanLocX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanLocY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  double meanX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds3[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  // Vds4 (secondary tracks)
  double meanLocX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanLocY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaLocY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  double meanX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaX_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double meanY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};
  double sigmaY_Vds4[8] = {0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000};

  SetCAMembers(meanLocX_Vds1, sigmaLocX_Vds1, meanLocY_Vds1, sigmaLocY_Vds1, meanLocX_Vds2, sigmaLocX_Vds2, meanLocY_Vds2, sigmaLocY_Vds2, meanLocX_Vds3, sigmaLocX_Vds3, meanLocY_Vds3, sigmaLocY_Vds3, meanLocX_Vds4, sigmaLocX_Vds4, meanLocY_Vds4, sigmaLocY_Vds4, meanX_Vds1, sigmaX_Vds1, meanY_Vds1, sigmaY_Vds1, meanX_Vds2, sigmaX_Vds2, meanY_Vds2, sigmaY_Vds2, meanX_Vds3, sigmaX_Vds3, meanY_Vds3, sigmaY_Vds3, meanX_Vds4, sigmaX_Vds4, meanY_Vds4, sigmaY_Vds4);
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorRotation_JuraDec2016_mv1(int run_id) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorRotation_JuraDec2016_mv1: Setting sensor rotation for Jura Arm."<< " runid=" << run_id;
    INFO(info);

  double rotZ_Vds1_0 = -0.001077150897938;
  double rotZ_Vds2_0 = -0.001333294613345;
  double rotZ_Vds3_0 = 0.0009651053567656;
  double rotZ_Vds3_1 = 0.002172309935405;
  double rotZ_Vds4_0 = -0.03411835508146;
  double rotZ_Vds4_1 = -0.03613710979898;
  double rotZ_Vds4_2 = -0.04926660135119;
  double rotZ_Vds4_3 = -0.04926660135119;

  double rotY_Vds1_0 = -0.004228523725557;
  double rotY_Vds2_0 = 0.01313458228081;
  double rotY_Vds3_0 = -0.001163468864343;
  double rotY_Vds3_1 = 0.006851994070694;
  double rotY_Vds4_0 = 0.004850537460977;
  double rotY_Vds4_1 = 0.01226067755686;
  double rotY_Vds4_2 = 0;
  double rotY_Vds4_3 = 0.05999983341597;

  double rotX_Vds1_0 = 0.002525233766176;
  double rotX_Vds2_0 = -0.0007203204019122;
  double rotX_Vds3_0 = 0.002392604157584;
  double rotX_Vds3_1 = -0.003025546997931;
  double rotX_Vds4_0 = 0.010301824214;
  double rotX_Vds4_1 = -0.01141781094349;
  double rotX_Vds4_2 = 0;
  double rotX_Vds4_3 = -0.0112466629558;

  if (run_id > 300) {
    /* //based on run 1494 hand tunning
    rotX_Vds1_0 = 0.002525233766176;
    rotX_Vds2_0 = -0.00152000300129;
    rotX_Vds3_0 = -0.004044605362795 - 0.003;
    rotX_Vds3_1 = 0.0003659570787289;
    rotX_Vds4_0 =  -0.008989249763769;
    rotX_Vds4_1 = -0.008989249763769;
    rotX_Vds4_2 = -0.003114992646512;
    rotX_Vds4_3 = -0.01264625952445;

    rotY_Vds1_0 = -0.004228523725557;
    rotY_Vds2_0 = 0.01107346554849;
    rotY_Vds3_0 = 0.003420866321225;
    rotY_Vds3_1 = 0.01106057476881;
    rotY_Vds4_0 = -0.001875635037882;
    rotY_Vds4_1 = -0.001875635037882;
    rotY_Vds4_2 = 0.09982030596094;
    rotY_Vds4_3 = 0.09239895117879;

    rotZ_Vds1_0 = -0.001077150897938;
    rotZ_Vds2_0 = -0.001988461318324;
    rotZ_Vds3_0 = 0.0008332944005683 - 0.0004;
    rotZ_Vds3_1 = 0.001368606858636;
    rotZ_Vds4_0 = -0.005590195951921 - 0.00094;
    rotZ_Vds4_1 = -0.008904756991813;
    rotZ_Vds4_2 = -0.003868326967948 - 0.0017;
    rotZ_Vds4_3 = -0.003529297926015;
    */
    rotZ_Vds1_0 = 0;
    rotZ_Vds2_0 = -0.0007606737574483;
    rotZ_Vds3_0 = 0.001791322393959;
    rotZ_Vds3_1 = 0.003073113645296;
    rotZ_Vds4_0 = -0.004876847280671;
    rotZ_Vds4_1 = -0.006617490586695;
    rotY_Vds1_0 = 0;
    rotY_Vds2_0 = 0.007929535037832;
    rotY_Vds3_0 = 0.00678684091879;
    rotY_Vds3_1 = 0.01597817556361;
    rotY_Vds4_0 = 0.009640138452539;
    rotY_Vds4_1 = 0.01676232415152;
    rotX_Vds1_0 = 0;
    rotX_Vds2_0 = -0.00261954465003;
    rotX_Vds3_0 = 0.001883803662747;
    rotX_Vds3_1 = -0.001997533664629;
    rotX_Vds4_0 = 0.01433297906138;
    rotX_Vds4_1 = -0.01465099022064;

    rotZ_Vds4_3 = -0.001504566206239;
    rotY_Vds4_3 = 0.1067804635883;
    rotX_Vds4_3 = -0.01304942709162;
    rotZ_Vds4_2 = -0.004247184306278;
    rotY_Vds4_2 = 0.1095980053803;
    rotX_Vds4_2 = 0.006533483949134;
  }

  fRotX[0] = rotX_Vds1_0;
  fRotX[1] = rotX_Vds2_0;
  fRotX[2] = rotX_Vds3_0;
  fRotX[3] = rotX_Vds3_1;
  fRotX[4] = rotX_Vds4_0;
  fRotX[5] = rotX_Vds4_1;
  fRotX[6] = rotX_Vds4_2;
  fRotX[7] = rotX_Vds4_3;

  fRotY[0] = rotY_Vds1_0;
  fRotY[1] = rotY_Vds2_0;
  fRotY[2] = rotY_Vds3_0;
  fRotY[3] = rotY_Vds3_1;
  fRotY[4] = rotY_Vds4_0;
  fRotY[5] = rotY_Vds4_1;
  fRotY[6] = rotY_Vds4_2;
  fRotY[7] = rotY_Vds4_3;

  fRotZ[0] = rotZ_Vds1_0;
  fRotZ[1] = rotZ_Vds2_0;
  fRotZ[2] = rotZ_Vds3_0;
  fRotZ[3] = rotZ_Vds3_1;
  fRotZ[4] = rotZ_Vds4_0;
  fRotZ[5] = rotZ_Vds4_1;
  fRotZ[6] = rotZ_Vds4_2;
  fRotZ[7] = rotZ_Vds4_3;

  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;

    ostringstream info_fin;
    info_fin <<  "Jura rots: " << fRotZ[4] << " " << fRotZ[5] << "  " << fRotZ[6] << " " << fRotZ[7] << "\n"
	  << "fDrotX = " << fDrotX << "   fDrotY = " << fDrotY << "    fDrotZ = " << fDrotZ << "\n"
	  << "0: " << fRotX[0] << " " << fRotY[0] << " " << fRotZ[0] << "\n"
	  << "1: " << fRotX[1] << " " << fRotY[1] << " " << fRotZ[1] << "\n"
	  << "2: " << fRotX[2] << " " << fRotY[2] << " " << fRotZ[2] << "\n"
	  << "3: " << fRotX[3] << " " << fRotY[3] << " " << fRotZ[3] << "\n"
	  << "4: " << fRotX[4] << " " << fRotY[4] << " " << fRotZ[4] << "\n"
	  << "5: " << fRotX[5] << " " << fRotY[5] << " " << fRotZ[5] << "\n"
	  << "6: " << fRotX[6] << " " << fRotY[6] << " " << fRotZ[6] << "\n"
	  << "7: " << fRotX[7] << " " << fRotY[7] << " " << fRotZ[7] << "\n";
    INFO(info_fin);

}
//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorRotation_JuraNov2018(int run_id) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorRotation_JuraNov2018: Setting sensor rotation for Jura Arm."<< " runid=" << run_id;
    INFO(info);


  //JURA:
	double rotX_Vdj1_0 = 0;
	double rotX_Vdj2_0 = -0.003805265625976488;
	double rotX_Vdj3_0 = 0.0006856180553823504;
	double rotX_Vdj3_1 = -0.002341891159920111;
	double rotX_Vdj4_0 = 0.01290384378273127;
	double rotX_Vdj4_1 = -0.01480513314115244;
	double rotX_Vdj4_2 = 0.004198979592726484;
	double rotX_Vdj4_3 = -0.01017593697835698;
	double rotY_Vdj1_0 = 0;
	double rotY_Vdj2_0 = 0.00611108415455372;
	double rotY_Vdj3_0 = 0.005178334696321465;
	double rotY_Vdj3_1 = 0.01094299550853452;
	double rotY_Vdj4_0 = 0.00921539793445629;
	double rotY_Vdj4_1 = 0.01728122885453472;
	double rotY_Vdj4_2 = 0.1150741333633511;
	double rotY_Vdj4_3 = 0.1025287302175883;
	double rotZ_Vdj1_0 = 0;
	double rotZ_Vdj2_0 = -0.0006346247233531069;
	double rotZ_Vdj3_0 = 0.002061627840650872;
	double rotZ_Vdj3_1 = 0.003239435802170329;
	double rotZ_Vdj4_0 = -0.004520146167220872;
	double rotZ_Vdj4_1 = -0.006463341886944501;
	double rotZ_Vdj4_2 = -0.003821259077932772;
	double rotZ_Vdj4_3 = -0.001526849158316024;

  fRotX[0] = rotX_Vdj1_0;
  fRotX[1] = rotX_Vdj2_0;
  fRotX[2] = rotX_Vdj3_0;
  fRotX[3] = rotX_Vdj3_1;
  fRotX[4] = rotX_Vdj4_0;
  fRotX[5] = rotX_Vdj4_1;
  fRotX[6] = rotX_Vdj4_2;
  fRotX[7] = rotX_Vdj4_3;

  fRotY[0] = rotY_Vdj1_0;
  fRotY[1] = rotY_Vdj2_0;
  fRotY[2] = rotY_Vdj3_0;
  fRotY[3] = rotY_Vdj3_1;
  fRotY[4] = rotY_Vdj4_0;
  fRotY[5] = rotY_Vdj4_1;
  fRotY[6] = rotY_Vdj4_2;
  fRotY[7] = rotY_Vdj4_3;

  fRotZ[0] = rotZ_Vdj1_0;
  fRotZ[1] = rotZ_Vdj2_0;
  fRotZ[2] = rotZ_Vdj3_0;
  fRotZ[3] = rotZ_Vdj3_1;
  fRotZ[4] = rotZ_Vdj4_0;
  fRotZ[5] = rotZ_Vdj4_1;
  fRotZ[6] = rotZ_Vdj4_2;
  fRotZ[7] = rotZ_Vdj4_3;

  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;

   ostringstream info_fin;
    info_fin <<  "Jura rots: " << fRotZ[4] << " " << fRotZ[5] << "  " << fRotZ[6] << " " << fRotZ[7] << "\n"
	  << "fDrotX = " << fDrotX << "   fDrotY = " << fDrotY << "    fDrotZ = " << fDrotZ << "\n"
	  << "0: " << fRotX[0] << " " << fRotY[0] << " " << fRotZ[0] << "\n"
	  << "1: " << fRotX[1] << " " << fRotY[1] << " " << fRotZ[1] << "\n"
	  << "2: " << fRotX[2] << " " << fRotY[2] << " " << fRotZ[2] << "\n"
	  << "3: " << fRotX[3] << " " << fRotY[3] << " " << fRotZ[3] << "\n"
	  << "4: " << fRotX[4] << " " << fRotY[4] << " " << fRotZ[4] << "\n"
	  << "5: " << fRotX[5] << " " << fRotY[5] << " " << fRotZ[5] << "\n"
	  << "6: " << fRotX[6] << " " << fRotY[6] << " " << fRotZ[6] << "\n"
	  << "7: " << fRotX[7] << " " << fRotY[7] << " " << fRotZ[7] << "\n";
    INFO(info_fin);

}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorRotation_JuraDec2016(int /*run_id*/) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorRotation_JuraDec2016: Setting sensor rotation for Jura Arm";
    INFO(info);

  double rotZ_Vds1_0 = 0.00001689;
  double rotZ_Vds2_0 = -0.00079875;
  double rotZ_Vds3_0 = 0.0007;
  double rotZ_Vds3_1 = 0.00197;
  double rotZ_Vds4_0 = -0.0351;
  double rotZ_Vds4_1 = -0.03708;
  double rotZ_Vds4_2 = -0.0503;
  double rotZ_Vds4_3 = -0.0503;

  double rotY_Vds4_2 = 6.0 * TMath::DegToRad();
  double rotY_Vds4_3 = 6.0 * TMath::DegToRad();

  fRotX[0] = 0;
  fRotX[1] = 0;
  fRotX[2] = 0;
  fRotX[3] = 0;
  fRotX[4] = 0;
  fRotX[5] = 0;
  fRotX[6] = 0;
  fRotX[7] = 0;

  fRotY[0] = 0;
  fRotY[1] = 0;
  fRotY[2] = 0;
  fRotY[3] = 0;
  fRotY[4] = 0;
  fRotY[5] = 0;
  fRotY[6] = rotY_Vds4_2;
  fRotY[7] = rotY_Vds4_3;

  fRotZ[0] = rotZ_Vds1_0;
  fRotZ[1] = rotZ_Vds2_0;
  fRotZ[2] = rotZ_Vds3_0;
  fRotZ[3] = rotZ_Vds3_1;
  fRotZ[4] = rotZ_Vds4_0;
  fRotZ[5] = rotZ_Vds4_1;
  fRotZ[6] = rotZ_Vds4_2;
  fRotZ[7] = rotZ_Vds4_3;

  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorRotation_SaleveDec2016_mv1(int run_id) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorRotation_SaleveDec2016_mv1: Setting sensor rotation for Saleve Arm:"<< " runid=" << run_id;;
    INFO(info); 
    
  double rotZ_Vds1_0 = -0.002029551416201;
  double rotZ_Vds2_0 = -0.004397669961837;
  double rotZ_Vds3_0 = 0.0009790805104435;
  double rotZ_Vds3_1 = 0.001814672612311;
  double rotZ_Vds4_0 = 0.0112099866484;
  double rotZ_Vds4_1 = 0.01236909132521;
  double rotZ_Vds4_2 = 0.0008688995346964;
  double rotZ_Vds4_3 = 0.001755049655373;

  double rotY_Vds1_0 = 0.00976946612238;
  double rotY_Vds2_0 = 0.0006590761406203;
  double rotY_Vds3_0 = -0.002177815091619;
  double rotY_Vds3_1 = 0.009779322026311;
  double rotY_Vds4_0 = 0.001988322769912;
  double rotY_Vds4_1 = 0.01318601349458;
  double rotY_Vds4_2 = -0.05999863420543;
  double rotY_Vds4_3 = -0.05966692958544;

  double rotX_Vds1_0 = 0.00175562231091;
  double rotX_Vds2_0 = -0.0006708305935796;
  double rotX_Vds3_0 = 0.01087737560298;
  double rotX_Vds3_1 = -0.01138932540372;
  double rotX_Vds4_0 = 0.001265873671791;
  double rotX_Vds4_1 = -0.006669128962208;
  double rotX_Vds4_2 = -0.0007651277058818;
  double rotX_Vds4_3 = 0.006614750368053;

  if (run_id > 300) {
    /* based on run 1494 hand tunning
    rotX_Vds4_2 += -0.005;
    rotZ_Vds4_2 +=  0.006;
    rotZ_Vds4_3 +=  0.00487;
    */
    rotZ_Vds1_0 = 0;
    rotZ_Vds2_0 = -0.002402976454345;
    rotZ_Vds3_0 = 0.002468302147118;
    rotZ_Vds3_1 = 0.002967674352445;
    rotZ_Vds4_1 = 0.01747033607538;
    rotZ_Vds4_0 = 0.01670677748249;
    rotY_Vds1_0 = 0;
    rotY_Vds2_0 = 0.004602856356315;
    rotY_Vds3_0 = 0.008456887078114;
    rotY_Vds3_1 = 0.01034076544688;
    rotY_Vds4_1 = 0.01204902350159;
    rotY_Vds4_0 = 0.008159548813939;
    rotX_Vds1_0 = 0;
    rotX_Vds2_0 = -0.0001733583461049;
    rotX_Vds3_0 = 0.01482887765152;
    rotX_Vds3_1 = -0.01366236085178;
    rotX_Vds4_1 = -0.008605750218265;
    rotX_Vds4_0 = 0.00228327776192;
    rotZ_Vds4_2 = 0.009014005406779;
    rotY_Vds4_2 = -0.08949541204317;
    rotX_Vds4_2 = 0.001943320525718;
    rotZ_Vds4_3 = 0.008775899300644;
    rotY_Vds4_3 = -0.0800007770111;
    rotX_Vds4_3 = -0.0003027196075689;
  }

  fRotX[0] = rotX_Vds1_0;
  fRotX[1] = rotX_Vds2_0;
  fRotX[2] = rotX_Vds3_0;
  fRotX[3] = rotX_Vds3_1;
  fRotX[4] = rotX_Vds4_0;
  fRotX[5] = rotX_Vds4_1;
  fRotX[6] = rotX_Vds4_2;
  fRotX[7] = rotX_Vds4_3;

  fRotY[0] = rotY_Vds1_0;
  fRotY[1] = rotY_Vds2_0;
  fRotY[2] = rotY_Vds3_0;
  fRotY[3] = rotY_Vds3_1;
  fRotY[4] = rotY_Vds4_0;
  fRotY[5] = rotY_Vds4_1;
  fRotY[6] = rotY_Vds4_2;
  fRotY[7] = rotY_Vds4_3;

  fRotZ[0] = rotZ_Vds1_0;
  fRotZ[1] = rotZ_Vds2_0;
  fRotZ[2] = rotZ_Vds3_0;
  fRotZ[3] = rotZ_Vds3_1;
  fRotZ[4] = rotZ_Vds4_0;
  fRotZ[5] = rotZ_Vds4_1;
  fRotZ[6] = rotZ_Vds4_2;
  fRotZ[7] = rotZ_Vds4_3;

  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;

   ostringstream info_fin;
    info_fin <<  "Saleve rots: " << fRotZ[4] << " " << fRotZ[5] << "  " << fRotZ[6] << " " << fRotZ[7] << "\n"
	  << "fDrotX = " << fDrotX << "   fDrotY = " << fDrotY << "    fDrotZ = " << fDrotZ << "\n"
	  << "0: " << fRotX[0] << " " << fRotY[0] << " " << fRotZ[0] << "\n"
	  << "1: " << fRotX[1] << " " << fRotY[1] << " " << fRotZ[1] << "\n"
	  << "2: " << fRotX[2] << " " << fRotY[2] << " " << fRotZ[2] << "\n"
	  << "3: " << fRotX[3] << " " << fRotY[3] << " " << fRotZ[3] << "\n"
	  << "4: " << fRotX[4] << " " << fRotY[4] << " " << fRotZ[4] << "\n"
	  << "5: " << fRotX[5] << " " << fRotY[5] << " " << fRotZ[5] << "\n"
	  << "6: " << fRotX[6] << " " << fRotY[6] << " " << fRotZ[6] << "\n"
	  << "7: " << fRotX[7] << " " << fRotY[7] << " " << fRotZ[7] << "\n";
    INFO(info_fin);

}



//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorRotation_SaleveNov2018(int run_id) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorRotation_SaleveNov2018: Setting sensor rotation for Saleve Arm:"<< " runid=" << run_id;;
    INFO(info); 

  //SALEVE:
	double rotX_Vds1_0 = 0;
	double rotX_Vds2_0 = -0.003179498046911193;
	double rotX_Vds3_0 = 0.01237810975148561;
	double rotX_Vds3_1 = -0.01471077540540409;
	double rotX_Vds4_0 = 0.0003659713268111126;
	double rotX_Vds4_1 = -0.008722211034154562;
	double rotX_Vds4_2 = 0.002495232651873243;
	double rotX_Vds4_3 = -0.007194152710942758;
	double rotY_Vds1_0 = 0;
	double rotY_Vds2_0 = 0.0001937501299672101;
	double rotY_Vds3_0 = 0.002396727424447187;
	double rotY_Vds3_1 = 0.003797999562713775;
	double rotY_Vds4_0 = 0.003683276876527766;
	double rotY_Vds4_1 = 0.007252068662043309;
	double rotY_Vds4_2 = -0.1005492626959199;
	double rotY_Vds4_3 = -0.09030519607686563;
	double rotZ_Vds1_0 = 0;
	double rotZ_Vds2_0 = -0.002181269256232229;
	double rotZ_Vds3_0 = 0.002594467455915656;
	double rotZ_Vds3_1 = 0.003393339132357273;
	double rotZ_Vds4_0 = 0.01683265956769983;
	double rotZ_Vds4_1 = 0.01796274628785534;
	double rotZ_Vds4_2 = 0.009472102252258238;
	double rotZ_Vds4_3 = 0.008510565815155885;
  

  fRotX[0] = rotX_Vds1_0;
  fRotX[1] = rotX_Vds2_0;
  fRotX[2] = rotX_Vds3_0;
  fRotX[3] = rotX_Vds3_1;
  fRotX[4] = rotX_Vds4_0;
  fRotX[5] = rotX_Vds4_1;
  fRotX[6] = rotX_Vds4_2;
  fRotX[7] = rotX_Vds4_3;

  fRotY[0] = rotY_Vds1_0;
  fRotY[1] = rotY_Vds2_0;
  fRotY[2] = rotY_Vds3_0;
  fRotY[3] = rotY_Vds3_1;
  fRotY[4] = rotY_Vds4_0;
  fRotY[5] = rotY_Vds4_1;
  fRotY[6] = rotY_Vds4_2;
  fRotY[7] = rotY_Vds4_3;

  fRotZ[0] = rotZ_Vds1_0;
  fRotZ[1] = rotZ_Vds2_0;
  fRotZ[2] = rotZ_Vds3_0;
  fRotZ[3] = rotZ_Vds3_1;
  fRotZ[4] = rotZ_Vds4_0;
  fRotZ[5] = rotZ_Vds4_1;
  fRotZ[6] = rotZ_Vds4_2;
  fRotZ[7] = rotZ_Vds4_3;

  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;

  ostringstream info_fin;
    info_fin <<  "Saleve rots: " << fRotZ[4] << " " << fRotZ[5] << "  " << fRotZ[6] << " " << fRotZ[7] << "\n"
	  << "fDrotX = " << fDrotX << "   fDrotY = " << fDrotY << "    fDrotZ = " << fDrotZ << "\n"
	  << "0: " << fRotX[0] << " " << fRotY[0] << " " << fRotZ[0] << "\n"
	  << "1: " << fRotX[1] << " " << fRotY[1] << " " << fRotZ[1] << "\n"
	  << "2: " << fRotX[2] << " " << fRotY[2] << " " << fRotZ[2] << "\n"
	  << "3: " << fRotX[3] << " " << fRotY[3] << " " << fRotZ[3] << "\n"
	  << "4: " << fRotX[4] << " " << fRotY[4] << " " << fRotZ[4] << "\n"
	  << "5: " << fRotX[5] << " " << fRotY[5] << " " << fRotZ[5] << "\n"
	  << "6: " << fRotX[6] << " " << fRotY[6] << " " << fRotZ[6] << "\n"
	  << "7: " << fRotX[7] << " " << fRotY[7] << " " << fRotZ[7] << "\n";
    INFO(info_fin);

}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorRotation_SaleveDec2016(int /*run_id*/) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorRotation_SaleveDec2016: Setting sensor rotation for Saleve Arm";
    INFO(info); 

  // valid on 18Feb2017
  double rotZ_Vds1_0 = 0.0;
  double rotZ_Vds2_0 = -0.00365;
  double rotZ_Vds3_0 = 0.000238;
  double rotZ_Vds3_1 = 0.0007713;
  double rotZ_Vds4_0 = 0.00864511;
  double rotZ_Vds4_1 = 0.00981514;
  double rotZ_Vds4_2 = 0.0;
  double rotZ_Vds4_3 = -0.0035;

  double rotY_Vds4_2 = -6.0 * TMath::DegToRad();
  double rotY_Vds4_3 = -6.0 * TMath::DegToRad() + 0.02;

  // valid on 18Feb2017
  fRotX[0] = 0;
  fRotX[1] = 0;
  fRotX[2] = 0.008 - 0.000625;
  fRotX[3] = 0.008 - 0.000625;
  fRotX[4] = 0.006 - 0.0004;
  fRotX[5] = 0.006 - 0.0004;
  fRotX[6] = 0;
  fRotX[7] = 0.002;

  fRotY[0] = 0;
  fRotY[1] = 0;
  fRotY[2] = 0;
  fRotY[3] = 0;
  fRotY[4] = 0;
  fRotY[5] = 0;
  fRotY[6] = rotY_Vds4_2;
  fRotY[7] = rotY_Vds4_3;

  fRotZ[0] = rotZ_Vds1_0;
  fRotZ[1] = rotZ_Vds2_0;
  fRotZ[2] = rotZ_Vds3_0;
  fRotZ[3] = rotZ_Vds3_1;
  fRotZ[4] = rotZ_Vds4_0;
  fRotZ[5] = rotZ_Vds4_1;
  fRotZ[6] = rotZ_Vds4_2;
  fRotZ[7] = rotZ_Vds4_3;

  // for tunning
  fRotZ[fSensorId] = fRotZ[fSensorId] + fDrotZ;
  fRotY[fSensorId] = fRotY[fSensorId] + fDrotY;
  fRotX[fSensorId] = fRotX[fSensorId] + fDrotX;

  // if(fSensorId==2)fRotX[fSensorId+1] = fRotX[fSensorId+1] + fDrotX*0.05;
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorGeometry_JuraDec2016_mv1(int run_id) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorGeometry_JuraDec2016_mv1: Setting sensor geometry for Jura Arm";
    INFO(info); 
    
  double VolumeX_Vds1_0 = 0.1368855771607;
  double VolumeY_Vds1_0 = -0.1091340091322;
  double VolumeZ_Vds1_0 = -0.1210278592745;

  double VolumeX_Vds2_0 = -0.02447030454697;
  double VolumeY_Vds2_0 = -0.1009601314079;
  double VolumeZ_Vds2_0 = 50.0297843027;

  double VolumeX_Vds3_0 = 1.52049959432;
  double VolumeY_Vds3_0 = -10.58683310383;
  double VolumeZ_Vds3_0 = 100.0121171975;

  double VolumeX_Vds3_1 = 1.503298814002;
  double VolumeY_Vds3_1 = 11.21095151991;
  double VolumeZ_Vds3_1 = 99.95562845583;

  double VolumeX_Vds4_0 = 3.596373926511;
  double VolumeY_Vds4_0 = -11.12835412828;
  double VolumeZ_Vds4_0 = 149.9121174084;

  double VolumeX_Vds4_1 = 2.812465341623;
  double VolumeY_Vds4_1 = 10.63050089359;
  double VolumeZ_Vds4_1 = 149.973888799;

  double VolumeX_Vds4_2 = 0;
  double VolumeY_Vds4_2 = 0;
  double VolumeZ_Vds4_2 = 0;

  double VolumeX_Vds4_3 = 11.58727136231;
  double VolumeY_Vds4_3 = 10.34854015557;
  double VolumeZ_Vds4_3 = 144.5578450037;

  if (run_id > 300) {
    /* // based on run 1494 hand tuning
    VolumeX_Vds1_0 = 0.0903336613607   + 0.00140494;
    VolumeX_Vds2_0 = -0.03233480219526 + 0.000337177;
    VolumeX_Vds3_0 = 1.654018716392    - 6.43779e-05;
    VolumeX_Vds3_1 = 1.631337788366    - 0.00785947;
    VolumeX_Vds4_0 = 2.815024544775+0.1848 + 0.00389424;
    VolumeX_Vds4_1 = 2.815024544775     + 0.00228749;
    VolumeX_Vds4_2 = 12.07309367102 + 0.01748;
    VolumeX_Vds4_3 = 12.02499054956 + 0.003542;

    VolumeY_Vds1_0 = -0.1319711440422   - 0.000294733;
    VolumeY_Vds2_0 = -0.07973324647261  + 0.000393217;
    VolumeY_Vds3_0 = -10.53491732484    + 0.00167946;
    VolumeY_Vds3_1 = 11.26335154434     - 0.000481122;
    VolumeY_Vds4_0 = -10.58847523866-0.00444 + 0.00291669;
    VolumeY_Vds4_1 = 11.18658154678     - 0.00421351;
    VolumeY_Vds4_2 =-10.53121920953 - 0.003483;
    VolumeY_Vds4_3 = 11.23023408421 + 0.003767;

    VolumeZ_Vds1_0 = -0.1210278592745;
    VolumeZ_Vds2_0 = 50.04250667007;
    VolumeZ_Vds3_0 = 99.99070397969;
    VolumeZ_Vds3_1 = 99.9455097154;
    VolumeZ_Vds4_0 = 149.9152792855;
    VolumeZ_Vds4_1 = 149.9152792855;
    VolumeZ_Vds4_2 = 144.3244905877;
    VolumeZ_Vds4_3 = 144.4124135478;
    */

/*    VolumeZ_Vds1_0 = 0;
    VolumeZ_Vds2_0 = 50.0;
    VolumeZ_Vds3_0 = 100.0;
    VolumeZ_Vds3_1 = 100.0;
    VolumeZ_Vds4_1 = 149.90;
    VolumeZ_Vds4_0 = 149.90;
    VolumeY_Vds1_0 = 0;
    VolumeY_Vds2_0 = 0.0;
    VolumeY_Vds3_0 = -10.0;
    VolumeY_Vds4_0 = -10.0;
    VolumeY_Vds3_1 = 11.0;
    VolumeY_Vds4_1 = 11.0;
    VolumeX_Vds1_0 = 0;
    VolumeX_Vds2_0 = -0.0;
    VolumeX_Vds3_0 = 1.56;
    VolumeX_Vds3_1 = 1.56;
    VolumeX_Vds4_1 = 2.76;
    VolumeX_Vds4_0 = 2.90;
    VolumeZ_Vds4_3 = 144.0;
    VolumeY_Vds4_3 = 11.0;
    VolumeX_Vds4_3 = 11.0;
    VolumeZ_Vds4_2 = 144.0;
    VolumeY_Vds4_2 = -10.0;
    VolumeX_Vds4_2 = 11.0;
  */  
    VolumeZ_Vds1_0 = 0;
    VolumeZ_Vds2_0 = 50.1716392241;
    VolumeZ_Vds3_0 = 100.1441791009;
    VolumeZ_Vds3_1 = 100.0737824693;
    VolumeZ_Vds4_1 = 149.9608448119;
    VolumeZ_Vds4_0 = 149.9104174119;
    VolumeY_Vds1_0 = 0;
    VolumeY_Vds2_0 = 0.05392611472887;
    VolumeY_Vds3_0 = -10.4015557085;
    VolumeY_Vds4_0 = -10.45940058177;
    VolumeY_Vds3_1 = 11.39405339932;
    VolumeY_Vds4_1 = 11.31072888905;
    VolumeX_Vds1_0 = 0;
    VolumeX_Vds2_0 = -0.1180399455156;
    VolumeX_Vds3_0 = 1.559755433126;
    VolumeX_Vds3_1 = 1.561798380756;
    VolumeX_Vds4_1 = 2.765623280863;
    VolumeX_Vds4_0 = 2.909311492743;
    VolumeZ_Vds4_3 = 144.262296082;
    VolumeY_Vds4_3 = 11.34198036382;
    VolumeX_Vds4_3 = 11.95014448732;
    VolumeZ_Vds4_2 = 144.1991711016;
    VolumeY_Vds4_2 = -10.40721411454;
    VolumeX_Vds4_2 = 11.97852602446;
 
  }

  fVolumeX[0] = VolumeX_Vds1_0;
  fVolumeX[1] = VolumeX_Vds2_0;
  fVolumeX[2] = VolumeX_Vds3_0;
  fVolumeX[3] = VolumeX_Vds3_1;
  fVolumeX[4] = VolumeX_Vds4_0;
  fVolumeX[5] = VolumeX_Vds4_1;
  fVolumeX[6] = VolumeX_Vds4_2;
  fVolumeX[7] = VolumeX_Vds4_3;

  fVolumeY[0] = VolumeY_Vds1_0;
  fVolumeY[1] = VolumeY_Vds2_0;
  fVolumeY[2] = VolumeY_Vds3_0;
  fVolumeY[3] = VolumeY_Vds3_1;
  fVolumeY[4] = VolumeY_Vds4_0;
  fVolumeY[5] = VolumeY_Vds4_1;
  fVolumeY[6] = VolumeY_Vds4_2;
  fVolumeY[7] = VolumeY_Vds4_3;

  fVolumeZ[0] = VolumeZ_Vds1_0;
  fVolumeZ[1] = VolumeZ_Vds2_0;
  fVolumeZ[2] = VolumeZ_Vds3_0;
  fVolumeZ[3] = VolumeZ_Vds3_1;
  fVolumeZ[4] = VolumeZ_Vds4_0;
  fVolumeZ[5] = VolumeZ_Vds4_1;
  fVolumeZ[6] = VolumeZ_Vds4_2;
  fVolumeZ[7] = VolumeZ_Vds4_3;

  for (int i = 0; i < 8; i++) {
    fVolumeX[i] = fVolumeX[i] - VolumeX_Vds1_0;
    fVolumeY[i] = fVolumeY[i] - VolumeY_Vds1_0;
    fVolumeZ[i] = fVolumeZ[i] - VolumeZ_Vds1_0;
  }

    ostringstream info_fin;
    info_fin << "Jura: volumeX_Vds1_0=" << fVolumeX[0] << "   volumeX_Vds3_0=" << fVolumeX[2] << "  volumeX_Vds3_1=" << fVolumeX[3]<< "\n"
			<< "0: " << fVolumeX[0] << " " << fVolumeY[0] << " " << fVolumeZ[0] <<  "\n"
			<< "1: " << fVolumeX[1] << " " << fVolumeY[1] << " " << fVolumeZ[1] <<  "\n"
			<< "2: " << fVolumeX[2] << " " << fVolumeY[2] << " " << fVolumeZ[2] <<  "\n"
			<< "3: " << fVolumeX[3] << " " << fVolumeY[3] << " " << fVolumeZ[3] <<  "\n"
			<< "4: " << fVolumeX[4] << " " << fVolumeY[4] << " " << fVolumeZ[4] <<  "\n"
			<< "5: " << fVolumeX[5] << " " << fVolumeY[5] << " " << fVolumeZ[5] <<  "\n"
			<< "6: " << fVolumeX[6] << " " << fVolumeY[6] << " " << fVolumeZ[6] <<  "\n"
			<< "7: " << fVolumeX[7] << " " << fVolumeY[7] << " " << fVolumeZ[7] ;
    INFO(info_fin);

  fVolumeZ[fSensorId] = fVolumeZ[fSensorId] + fDz;
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorGeometry_JuraNov2018(int /* run_id */) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorGeometry_JuraNov2018: Setting sensor geometry for Jura Arm";
    INFO(info); 

	double VolumeX_Vdj1_0 = 0;
	double VolumeX_Vdj2_0 = -0.1223014141774704;
	double VolumeX_Vdj3_0 = 1.560134364320414;
	double VolumeX_Vdj3_1 = 1.567317519977404;
	double VolumeX_Vdj4_0 = 2.902218240885595;
	double VolumeX_Vdj4_1 = 2.763984128019565;
	double VolumeX_Vdj4_2 = 11.98964287150629;
	double VolumeX_Vdj4_3 = 11.97062541569299;
	double VolumeY_Vdj1_0 = 0;
	double VolumeY_Vdj2_0 = 0.05392703344366791;
	double VolumeY_Vdj3_0 = -10.40650329075579;
	double VolumeY_Vdj3_1 = 11.39073799157225;
	double VolumeY_Vdj4_0 = -10.46476699132205;
	double VolumeY_Vdj4_1 = 11.30737653434626;
	double VolumeY_Vdj4_2 = -10.41199051162143;
	double VolumeY_Vdj4_3 = 11.33928817221769;
	double VolumeZ_Vdj1_0 = 0;
	double VolumeZ_Vdj2_0 = 50.18440654199076;
	double VolumeZ_Vdj3_0 = 100.1402915611285;
	double VolumeZ_Vdj3_1 = 100.0984774543834;
	double VolumeZ_Vdj4_0 = 149.9083927471845;
	double VolumeZ_Vdj4_1 = 149.9817782505897;
	double VolumeZ_Vdj4_2 = 144.2944115854085;
	double VolumeZ_Vdj4_3 = 144.4382832887717;

  fVolumeX[0] = VolumeX_Vdj1_0;
  fVolumeX[1] = VolumeX_Vdj2_0;
  fVolumeX[2] = VolumeX_Vdj3_0;
  fVolumeX[3] = VolumeX_Vdj3_1;
  fVolumeX[4] = VolumeX_Vdj4_0;
  fVolumeX[5] = VolumeX_Vdj4_1;
  fVolumeX[6] = VolumeX_Vdj4_2;
  fVolumeX[7] = VolumeX_Vdj4_3;

  fVolumeY[0] = VolumeY_Vdj1_0;
  fVolumeY[1] = VolumeY_Vdj2_0;
  fVolumeY[2] = VolumeY_Vdj3_0;
  fVolumeY[3] = VolumeY_Vdj3_1;
  fVolumeY[4] = VolumeY_Vdj4_0;
  fVolumeY[5] = VolumeY_Vdj4_1;
  fVolumeY[6] = VolumeY_Vdj4_2;
  fVolumeY[7] = VolumeY_Vdj4_3;

  fVolumeZ[0] = VolumeZ_Vdj1_0;
  fVolumeZ[1] = VolumeZ_Vdj2_0;
  fVolumeZ[2] = VolumeZ_Vdj3_0;
  fVolumeZ[3] = VolumeZ_Vdj3_1;
  fVolumeZ[4] = VolumeZ_Vdj4_0;
  fVolumeZ[5] = VolumeZ_Vdj4_1;
  fVolumeZ[6] = VolumeZ_Vdj4_2;
  fVolumeZ[7] = VolumeZ_Vdj4_3;

  for (int i = 0; i < 8; i++) {
    fVolumeX[i] = fVolumeX[i] - VolumeX_Vdj1_0;
    fVolumeY[i] = fVolumeY[i] - VolumeY_Vdj1_0;
    fVolumeZ[i] = fVolumeZ[i] - VolumeZ_Vdj1_0;
  }

    ostringstream info_fin;
    info_fin << "Jura: volumeX_Vds1_0=" << fVolumeX[0] << "   volumeX_Vds3_0=" << fVolumeX[2] << "  volumeX_Vds3_1=" << fVolumeX[3]<< "\n"
			<< "0: " << fVolumeX[0] << " " << fVolumeY[0] << " " << fVolumeZ[0] <<  "\n"
			<< "1: " << fVolumeX[1] << " " << fVolumeY[1] << " " << fVolumeZ[1] <<  "\n"
			<< "2: " << fVolumeX[2] << " " << fVolumeY[2] << " " << fVolumeZ[2] <<  "\n"
			<< "3: " << fVolumeX[3] << " " << fVolumeY[3] << " " << fVolumeZ[3] <<  "\n"
			<< "4: " << fVolumeX[4] << " " << fVolumeY[4] << " " << fVolumeZ[4] <<  "\n"
			<< "5: " << fVolumeX[5] << " " << fVolumeY[5] << " " << fVolumeZ[5] <<  "\n"
			<< "6: " << fVolumeX[6] << " " << fVolumeY[6] << " " << fVolumeZ[6] <<  "\n"
			<< "7: " << fVolumeX[7] << " " << fVolumeY[7] << " " << fVolumeZ[7] ;
    INFO(info_fin);

  fVolumeZ[fSensorId] = fVolumeZ[fSensorId] + fDz;
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorGeometry_JuraDec2016(int run_id) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorGeometry_JuraDec2016: Setting sensor geometry for Jura Arm";
    INFO(info); 
  //
  double VolumeX_Vds1_0 = 0.0;
  double VolumeY_Vds1_0 = 0.0;
  double VolumeZ_Vds1_0 = 0.0;

  double VolumeX_Vds2_0 = 0.0;
  double VolumeY_Vds2_0 = 0.0;
  double VolumeZ_Vds2_0 = 50.0;

  double VolumeX_Vds3_0 = 1.5;
  double VolumeY_Vds3_0 = -10.6;
  double VolumeZ_Vds3_0 = 100.0;

  double VolumeX_Vds3_1 = 1.5;
  double VolumeY_Vds3_1 = 10.6;
  double VolumeZ_Vds3_1 = 100.0;

  double VolumeX_Vds4_0 = 3.0;
  double VolumeY_Vds4_0 = -10.6;
  double VolumeZ_Vds4_0 = 150.0;

  double VolumeX_Vds4_1 = 3.0;
  double VolumeY_Vds4_1 = 10.6;
  double VolumeZ_Vds4_1 = 150.0;

  double VolumeX_Vds4_2 = 12.5;
  double VolumeY_Vds4_2 = -10.6;
  double VolumeZ_Vds4_2 = 150.0 - 5.522;

  double VolumeX_Vds4_3 = 12.5 - 0.614 - 0.0134 - 0.1192 - 0.05417 - 0.05894 - 0.06311 - 0.00636;
  double VolumeY_Vds4_3 = 10.6 - 0.3831 - 0.00843 + 0.05884 + 0.02811 + 0.031 + 0.029 + 0.01365;
  double VolumeZ_Vds4_3 = 150.0 - 5.522;

  if (run_id >= 156) {
    VolumeX_Vds1_0 = VolumeX_Vds1_0 + 0.13988 + 0.00508461;
    VolumeY_Vds1_0 = VolumeY_Vds1_0 - 0.111871 - 0.00264779;
    VolumeX_Vds2_0 = VolumeX_Vds2_0 - 0.0251283 - 0.00100725;
    VolumeY_Vds2_0 = VolumeY_Vds2_0 - 0.0994618 + 0.00003;
    VolumeX_Vds3_0 = VolumeX_Vds3_0 + 0.0215244 - 0.00252343;
    VolumeY_Vds3_0 = VolumeY_Vds3_0 + 0.0191063 - 0.0056756;
    VolumeX_Vds3_1 = VolumeX_Vds3_1 + 0.00109159 + 0.00385736;
    VolumeY_Vds3_1 = VolumeY_Vds3_1 + 0.614944 + 0.00767852;
    VolumeX_Vds4_0 = VolumeX_Vds4_0 + 0.603613 + 0.000303359;
    VolumeY_Vds4_0 = VolumeY_Vds4_0 - 0.513785 - 0.0176856;
    VolumeX_Vds4_1 = VolumeX_Vds4_1 - 0.203238 + 0.0141308;
    VolumeY_Vds4_1 = VolumeY_Vds4_1 + 0.0414128 + 0.00953964;
  }

  fVolumeX[0] = VolumeX_Vds1_0;
  fVolumeX[1] = VolumeX_Vds2_0;
  fVolumeX[2] = VolumeX_Vds3_0;
  fVolumeX[3] = VolumeX_Vds3_1;
  fVolumeX[4] = VolumeX_Vds4_0;
  fVolumeX[5] = VolumeX_Vds4_1;
  fVolumeX[6] = VolumeX_Vds4_2;
  fVolumeX[7] = VolumeX_Vds4_3;

  fVolumeY[0] = VolumeY_Vds1_0;
  fVolumeY[1] = VolumeY_Vds2_0;
  fVolumeY[2] = VolumeY_Vds3_0;
  fVolumeY[3] = VolumeY_Vds3_1;
  fVolumeY[4] = VolumeY_Vds4_0;
  fVolumeY[5] = VolumeY_Vds4_1;
  fVolumeY[6] = VolumeY_Vds4_2;
  fVolumeY[7] = VolumeY_Vds4_3;

  fVolumeZ[0] = VolumeZ_Vds1_0;
  fVolumeZ[1] = VolumeZ_Vds2_0;
  fVolumeZ[2] = VolumeZ_Vds3_0;
  fVolumeZ[3] = VolumeZ_Vds3_1;
  fVolumeZ[4] = VolumeZ_Vds4_0;
  fVolumeZ[5] = VolumeZ_Vds4_1;
  fVolumeZ[6] = VolumeZ_Vds4_2;
  fVolumeZ[7] = VolumeZ_Vds4_3;

  if (fSensorId == 7) {
    fVolumeZ[fSensorId] = fVolumeZ[fSensorId] + fDz;
  }
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorGeometry_SaleveDec2016_mv1(int run_id) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorGeometry_SaleveDec2016_mv1: Setting sensor geometry for Saleve Arm";
    INFO(info); 

  double VolumeX_Vds1_0 = -0.1888393910569;
  double VolumeY_Vds1_0 = 0.01340957134613;
  double VolumeZ_Vds1_0 = 0.03055500682257;

  double VolumeX_Vds2_0 = 0.1564477117856;
  double VolumeY_Vds2_0 = 0.0646028721518;
  double VolumeZ_Vds2_0 = 49.99484950415;

  double VolumeX_Vds3_0 = -1.084603556914;
  double VolumeY_Vds3_0 = -10.97048689887;
  double VolumeZ_Vds3_0 = 98.93176829138;

  double VolumeX_Vds3_1 = -1.063238160189;
  double VolumeY_Vds3_1 = 10.81017555368;
  double VolumeZ_Vds3_1 = 98.92136116813;

  double VolumeX_Vds4_0 = -3.534082438535;
  double VolumeY_Vds4_0 = -10.84096861682;
  double VolumeZ_Vds4_0 = 149.2656718775;

  double VolumeX_Vds4_1 = -3.278034162402;
  double VolumeY_Vds4_1 = 10.92144124989;
  double VolumeZ_Vds4_1 = 149.2535000946;

  double VolumeX_Vds4_2 = -11.84097715346;
  double VolumeY_Vds4_2 = -10.95780042925;
  double VolumeZ_Vds4_2 = 143.5593651331;

  double VolumeX_Vds4_3 = -11.78112110205;
  double VolumeY_Vds4_3 = 10.85803349738;
  double VolumeZ_Vds4_3 = 143.5911373822;

  if (run_id > 300) {
    //////////  hand tuning for p+Pb
    VolumeX_Vds1_0 += -0.0373143 - 0.000415944;
    VolumeY_Vds1_0 += -0.0376786 - 0.00379813;
    VolumeX_Vds2_0 += 0.0822789 - 0.000866984;
    VolumeY_Vds2_0 += 0.0501522 + 0.00143758;
    VolumeX_Vds3_0 += 0.0625869 + 0.00155447;
    VolumeY_Vds3_0 += 0.0773752 + 0.0105226;
    VolumeX_Vds3_1 += 0.0566314 + 0.00155861;
    VolumeY_Vds3_1 += 0.075937 + 0.00854973;
    VolumeX_Vds4_0 += -0.0883276 + 4.69306e-05;
    VolumeY_Vds4_0 += -0.0817247 - 0.0071413;
    VolumeX_Vds4_1 += -0.0758553 - 0.00187708;
    VolumeY_Vds4_1 += -0.0840612 - 0.00957048;
    VolumeX_Vds4_2 += 0.7838 + 0.042 - 0.03388;
    VolumeY_Vds4_2 += -0.5138 - 0.00137 + 0.01227 + 0.0015;
    VolumeX_Vds4_3 += 0.8398 + 0.0598 + 0.01448;
    VolumeY_Vds4_3 += -0.5189 + 0.01272;

    VolumeZ_Vds4_2 += 0.2;
  }

  if (run_id > 609) {
    /*
    //////////  hand tuning for xela based on run 1494
    VolumeX_Vds4_2 += 0.02914;      VolumeY_Vds4_2 += -0.06362;
    VolumeX_Vds4_3 += 0.05371;      VolumeY_Vds4_3 += -0.0641;

    VolumeX_Vds1_0 += -0.00268863;    VolumeY_Vds1_0 += -0.00274306;
    VolumeX_Vds2_0 += -0.0112062;     VolumeY_Vds2_0 += 0.00308116;
    VolumeX_Vds3_0 += 0.0133786;      VolumeY_Vds3_0 += 0.0210239;
    VolumeX_Vds3_1 += 0.0123462;      VolumeY_Vds3_1 += -0.00929255;
    VolumeX_Vds4_0 += -0.021034;      VolumeY_Vds4_0 += 0.011134;
    VolumeX_Vds4_1 += 0.00920406;     VolumeY_Vds4_1 += -0.0232034;
    VolumeX_Vds4_2 += 0.0281203;      VolumeY_Vds4_2 += 0.0304504;
    VolumeX_Vds4_3 += 0.0278484;      VolumeY_Vds4_3 += -0.00938758;
    */

    VolumeZ_Vds1_0 = 0;
    VolumeZ_Vds2_0 = 49.82727714403;
    VolumeZ_Vds3_0 = 99.00178372593;
    VolumeZ_Vds3_1 = 98.97709945964;
    VolumeZ_Vds4_1 = 149.1208748434;
    VolumeZ_Vds4_0 = 149.1656016029;
    VolumeY_Vds1_0 = 0;
    VolumeY_Vds2_0 = 0.1527970710493;
    VolumeY_Vds3_0 = -10.83465811508;
    VolumeY_Vds3_1 = 10.94249722165;
    VolumeY_Vds4_1 = 10.8768255235;
    VolumeY_Vds4_0 = -10.88007622063;
    VolumeX_Vds1_0 = 0;
    VolumeX_Vds2_0 = 0.4743779936569;
    VolumeX_Vds3_0 = -0.7832193132034;
    VolumeX_Vds3_1 = -0.7329469858254;
    VolumeX_Vds4_1 = -3.054237537133;
    VolumeX_Vds4_0 = -3.429999374246;
    VolumeZ_Vds4_2 = 143.6028808398;
    VolumeY_Vds4_2 = -11.44307097507;
    VolumeX_Vds4_2 = -10.74352281642;
    VolumeZ_Vds4_3 = 143.5712428089;
    VolumeY_Vds4_3 = 10.36267152744;
    VolumeX_Vds4_3 = -10.51493441266;
  }

  fVolumeX[0] = VolumeX_Vds1_0;
  fVolumeX[1] = VolumeX_Vds2_0;
  fVolumeX[2] = VolumeX_Vds3_0;
  fVolumeX[3] = VolumeX_Vds3_1;
  fVolumeX[4] = VolumeX_Vds4_0;
  fVolumeX[5] = VolumeX_Vds4_1;
  fVolumeX[6] = VolumeX_Vds4_2;
  fVolumeX[7] = VolumeX_Vds4_3;

  fVolumeY[0] = VolumeY_Vds1_0;
  fVolumeY[1] = VolumeY_Vds2_0;
  fVolumeY[2] = VolumeY_Vds3_0;
  fVolumeY[3] = VolumeY_Vds3_1;
  fVolumeY[4] = VolumeY_Vds4_0;
  fVolumeY[5] = VolumeY_Vds4_1;
  fVolumeY[6] = VolumeY_Vds4_2;
  fVolumeY[7] = VolumeY_Vds4_3;

  fVolumeZ[0] = VolumeZ_Vds1_0;
  fVolumeZ[1] = VolumeZ_Vds2_0;
  fVolumeZ[2] = VolumeZ_Vds3_0;
  fVolumeZ[3] = VolumeZ_Vds3_1;
  fVolumeZ[4] = VolumeZ_Vds4_0;
  fVolumeZ[5] = VolumeZ_Vds4_1;
  fVolumeZ[6] = VolumeZ_Vds4_2;
  fVolumeZ[7] = VolumeZ_Vds4_3;

  for (int i = 0; i < 8; i++) {
    fVolumeX[i] = fVolumeX[i] - VolumeX_Vds1_0;
    fVolumeY[i] = fVolumeY[i] - VolumeY_Vds1_0;
    fVolumeZ[i] = fVolumeZ[i] - VolumeZ_Vds1_0;
  }

    ostringstream info_fin;
    info_fin << "Saleve: volumeX_Vds1_0=" << fVolumeX[0] << "   volumeX_Vds3_0=" << fVolumeX[2] << "  volumeX_Vds3_1=" << fVolumeX[3]<< "\n"
			<< "0: " << fVolumeX[0] << " " << fVolumeY[0] << " " << fVolumeZ[0] <<  "\n"
			<< "1: " << fVolumeX[1] << " " << fVolumeY[1] << " " << fVolumeZ[1] <<  "\n"
			<< "2: " << fVolumeX[2] << " " << fVolumeY[2] << " " << fVolumeZ[2] <<  "\n"
			<< "3: " << fVolumeX[3] << " " << fVolumeY[3] << " " << fVolumeZ[3] <<  "\n"
			<< "4: " << fVolumeX[4] << " " << fVolumeY[4] << " " << fVolumeZ[4] <<  "\n"
			<< "5: " << fVolumeX[5] << " " << fVolumeY[5] << " " << fVolumeZ[5] <<  "\n"
			<< "6: " << fVolumeX[6] << " " << fVolumeY[6] << " " << fVolumeZ[6] <<  "\n"
			<< "7: " << fVolumeX[7] << " " << fVolumeY[7] << " " << fVolumeZ[7] ;
    INFO(info_fin);

  fVolumeZ[fSensorId] = fVolumeZ[fSensorId] + fDz;
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorGeometry_SaleveNov2018(int /* run_id */) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorGeometry_SaleveNov2018: Setting sensor geometry for Saleve Arm";
    INFO(info); 

	double VolumeX_Vds1_0 = 0;
	double VolumeX_Vds2_0 = 0.4723372269140442;
	double VolumeX_Vds3_0 = -0.7797037182792339;
	double VolumeX_Vds3_1 = -0.7250047998755361;
	double VolumeX_Vds4_0 = -3.440362777992692;
	double VolumeX_Vds4_1 = -3.064209912433843;
	double VolumeX_Vds4_2 = -10.7647330383394;
	double VolumeX_Vds4_3 = -10.52986903239746;
	double VolumeY_Vds1_0 = 0;
	double VolumeY_Vds2_0 = 0.1558462220144408;
	double VolumeY_Vds3_0 = -10.83506345404083;
	double VolumeY_Vds3_1 = 10.94491891028382;
	double VolumeY_Vds4_0 = -10.88449758887995;
	double VolumeY_Vds4_1 = 10.87495061745295;
	double VolumeY_Vds4_2 = -11.44303992453929;
	double VolumeY_Vds4_3 = 10.36307841809384;
	double VolumeZ_Vds1_0 = 0;
	double VolumeZ_Vds2_0 = 49.81915316749762;
	double VolumeZ_Vds3_0 = 99.01074052298848;
	double VolumeZ_Vds3_1 = 99.05925423955934;
	double VolumeZ_Vds4_0 = 149.1576543247651;
	double VolumeZ_Vds4_1 = 149.2289319891831;
	double VolumeZ_Vds4_2 = 143.5867591473741;
	double VolumeZ_Vds4_3 = 143.6046453830875;

  fVolumeX[0] = VolumeX_Vds1_0;
  fVolumeX[1] = VolumeX_Vds2_0;
  fVolumeX[2] = VolumeX_Vds3_0;
  fVolumeX[3] = VolumeX_Vds3_1;
  fVolumeX[4] = VolumeX_Vds4_0;
  fVolumeX[5] = VolumeX_Vds4_1;
  fVolumeX[6] = VolumeX_Vds4_2;
  fVolumeX[7] = VolumeX_Vds4_3;

  fVolumeY[0] = VolumeY_Vds1_0;
  fVolumeY[1] = VolumeY_Vds2_0;
  fVolumeY[2] = VolumeY_Vds3_0;
  fVolumeY[3] = VolumeY_Vds3_1;
  fVolumeY[4] = VolumeY_Vds4_0;
  fVolumeY[5] = VolumeY_Vds4_1;
  fVolumeY[6] = VolumeY_Vds4_2;
  fVolumeY[7] = VolumeY_Vds4_3;

  fVolumeZ[0] = VolumeZ_Vds1_0;
  fVolumeZ[1] = VolumeZ_Vds2_0;
  fVolumeZ[2] = VolumeZ_Vds3_0;
  fVolumeZ[3] = VolumeZ_Vds3_1;
  fVolumeZ[4] = VolumeZ_Vds4_0;
  fVolumeZ[5] = VolumeZ_Vds4_1;
  fVolumeZ[6] = VolumeZ_Vds4_2;
  fVolumeZ[7] = VolumeZ_Vds4_3;

  for (int i = 0; i < 8; i++) {
    fVolumeX[i] = fVolumeX[i] - VolumeX_Vds1_0;
    fVolumeY[i] = fVolumeY[i] - VolumeY_Vds1_0;
    fVolumeZ[i] = fVolumeZ[i] - VolumeZ_Vds1_0;
  }

      ostringstream info_fin;
    info_fin << "Saleve: volumeX_Vds1_0=" << fVolumeX[0] << "   volumeX_Vds3_0=" << fVolumeX[2] << "  volumeX_Vds3_1=" << fVolumeX[3]<< "\n"
			<< "0: " << fVolumeX[0] << " " << fVolumeY[0] << " " << fVolumeZ[0] <<  "\n"
			<< "1: " << fVolumeX[1] << " " << fVolumeY[1] << " " << fVolumeZ[1] <<  "\n"
			<< "2: " << fVolumeX[2] << " " << fVolumeY[2] << " " << fVolumeZ[2] <<  "\n"
			<< "3: " << fVolumeX[3] << " " << fVolumeY[3] << " " << fVolumeZ[3] <<  "\n"
			<< "4: " << fVolumeX[4] << " " << fVolumeY[4] << " " << fVolumeZ[4] <<  "\n"
			<< "5: " << fVolumeX[5] << " " << fVolumeY[5] << " " << fVolumeZ[5] <<  "\n"
			<< "6: " << fVolumeX[6] << " " << fVolumeY[6] << " " << fVolumeZ[6] <<  "\n"
			<< "7: " << fVolumeX[7] << " " << fVolumeY[7] << " " << fVolumeZ[7] ;
    INFO(info_fin);


  fVolumeZ[fSensorId] = fVolumeZ[fSensorId] + fDz;
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorGeometry_SaleveDec2016(int run_id) {
    ostringstream info;
    info << "Na61ArmParameters::SetupSensorGeometry_SaleveDec2016: Setting sensor geometry for Saleve Arm";
    INFO(info); 

  double VolumeX_Vds1_0 = 0.0;
  double VolumeY_Vds1_0 = 0.0;
  double VolumeZ_Vds1_0 = 0.0;

  double VolumeX_Vds2_0 = 0.0;
  double VolumeY_Vds2_0 = 0.0;
  double VolumeZ_Vds2_0 = 50.0;

  double VolumeX_Vds3_0 = -1.5;
  double VolumeY_Vds3_0 = -10.6;
  double VolumeZ_Vds3_0 = 100.0;

  double VolumeX_Vds3_1 = -1.5;
  double VolumeY_Vds3_1 = 10.6;
  double VolumeZ_Vds3_1 = 100.0;

  double VolumeX_Vds4_0 = -3.0;
  double VolumeY_Vds4_0 = -10.6;
  double VolumeZ_Vds4_0 = 150.0;

  double VolumeX_Vds4_1 = -3.0;
  double VolumeY_Vds4_1 = 10.6;
  double VolumeZ_Vds4_1 = 150.0;

  double VolumeX_Vds4_2 = -12.5 + 0.563 + 0.0077 + 0.07 + 0.02438;
  double VolumeY_Vds4_2 = -10.6 - 0.3794 - 0.024 + 0.03558;
  double VolumeZ_Vds4_2 = 150.0 - 5.522 - 1.;

  double VolumeX_Vds4_3 = -12.5 + 0.624 + 0.0099 + 0.08092 - 0.003758 + 0.00592;
  double VolumeY_Vds4_3 = 10.6 + 0.3391 - 0.024 - 0.08146 - 0.003598;
  double VolumeZ_Vds4_3 = 150.0 - 5.522 - 1. - 0.18 - 0.07;

  // double VolumeX_Vds4_2 = -13.0 + 1.063 + 0.0077 + 0.07;
  // double VolumeY_Vds4_2 = -10.6 - 0.3794 - 0.024;
  // double VolumeZ_Vds4_2 =  150.0 - 5.522;

  // double VolumeX_Vds4_3 = -13.0 + 1.124 +0.0099;
  // double VolumeY_Vds4_3 =  10.6 + 0.3391-0.024;
  // double VolumeZ_Vds4_3 =  150.0 - 5.522;

  if (run_id >= 156) {
    VolumeX_Vds1_0 = VolumeX_Vds1_0 - 0.18702 + 0.00177199;
    VolumeY_Vds1_0 = VolumeY_Vds1_0 + 0.0159644 + 0.00278645;
    VolumeX_Vds2_0 = VolumeX_Vds2_0 + 0.156424 + 9.27313e-05;
    VolumeY_Vds2_0 = VolumeY_Vds2_0 + 0.0686147 - 0.00545221;
    VolumeX_Vds3_0 = VolumeX_Vds3_0 + 0.423041 - 0.00992251;
    VolumeY_Vds3_0 = VolumeY_Vds3_0 - 0.368374 - 0.00284804;
    VolumeX_Vds3_1 = VolumeX_Vds3_1 + 0.430443 + 0.000261024;
    VolumeY_Vds3_1 = VolumeY_Vds3_1 + 0.213552 - 0.00674657;
    VolumeX_Vds4_0 = VolumeX_Vds4_0 - 0.516104 - 0.0102149;
    VolumeY_Vds4_0 = VolumeY_Vds4_0 - 0.246314 + 0.0109612;
    VolumeX_Vds4_1 = VolumeX_Vds4_1 - 0.306784 + 0.0180117;
    VolumeY_Vds4_1 = VolumeY_Vds4_1 + 0.316557 + 0.00129917;
  }

  fVolumeX[0] = VolumeX_Vds1_0;
  fVolumeX[1] = VolumeX_Vds2_0;
  fVolumeX[2] = VolumeX_Vds3_0;
  fVolumeX[3] = VolumeX_Vds3_1;
  fVolumeX[4] = VolumeX_Vds4_0;
  fVolumeX[5] = VolumeX_Vds4_1;
  fVolumeX[6] = VolumeX_Vds4_2;
  fVolumeX[7] = VolumeX_Vds4_3;

  fVolumeY[0] = VolumeY_Vds1_0;
  fVolumeY[1] = VolumeY_Vds2_0;
  fVolumeY[2] = VolumeY_Vds3_0;
  fVolumeY[3] = VolumeY_Vds3_1;
  fVolumeY[4] = VolumeY_Vds4_0;
  fVolumeY[5] = VolumeY_Vds4_1;
  fVolumeY[6] = VolumeY_Vds4_2;
  fVolumeY[7] = VolumeY_Vds4_3;

  fVolumeZ[0] = VolumeZ_Vds1_0;
  fVolumeZ[1] = VolumeZ_Vds2_0;
  fVolumeZ[2] = VolumeZ_Vds3_0;
  fVolumeZ[3] = VolumeZ_Vds3_1;
  fVolumeZ[4] = VolumeZ_Vds4_0;
  fVolumeZ[5] = VolumeZ_Vds4_1;
  fVolumeZ[6] = VolumeZ_Vds4_2;
  fVolumeZ[7] = VolumeZ_Vds4_3;

  fVolumeZ[fSensorId] = fVolumeZ[fSensorId] + fDz * 0.1;

  // corelate rotation x with dz shift
  fVolumeZ[0] = fVolumeZ[0] - (154. + fVolumeY[0]) * TMath::Sin(fRotX[0]);
  fVolumeZ[1] = fVolumeZ[1] - (154. + fVolumeY[1]) * TMath::Sin(fRotX[1]);
  fVolumeZ[2] = fVolumeZ[2] - (154. + fVolumeY[2]) * TMath::Sin(fRotX[2]);
  fVolumeZ[3] = fVolumeZ[3] - (154. + fVolumeY[3]) * TMath::Sin(fRotX[3]);
  fVolumeZ[4] = fVolumeZ[4] - (154. + fVolumeY[4]) * TMath::Sin(fRotX[4]);
  fVolumeZ[5] = fVolumeZ[5] - (154. + fVolumeY[5]) * TMath::Sin(fRotX[5]);
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_JuraDec2016_field(int /* run_id */) {
  // dev cuts
    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_JuraDec2016_field: setting offsets and sigmas";
    INFO(info); 

  fOffsetx[0] = 0.00009;
  fSigmax[0] = 0.00842;
  fOffsety[0] = 0.00009;
  fSigmay[0] = 0.00507;
  fOffsetx[1] = 0.00008;
  fSigmax[1] = 0.00845;
  fOffsety[1] = 0.00016;
  fSigmay[1] = 0.00524;
  fOffsetx[2] = 0.00045;
  fSigmax[2] = 0.01040;
  fOffsety[2] = 0.00010;
  fSigmay[2] = 0.00544;
  fOffsetx[3] = 0.00056;
  fSigmax[3] = 0.01050;
  fOffsety[3] = -0.00024;
  fSigmay[3] = 0.00603;

  fOffsetx[4] = 0.0;
  fSigmax[4] = 0.0;
  fOffsety[4] = 0.0;
  fSigmay[4] = 0.0;
  fOffsetx[5] = 0.00135;
  fSigmax[5] = 0.01353;
  fOffsety[5] = 0.00057;
  fSigmay[5] = 0.00697;
  fOffsetx[6] = 0.00016;
  fSigmax[6] = 0.01241;
  fOffsety[6] = 0.00011;
  fSigmay[6] = 0.00574;
  fOffsetx[7] = 0.0;
  fSigmax[7] = 0.0;
  fOffsety[7] = 0.0;
  fSigmay[7] = 0.0;
  fOffsetx[8] = 0.00218;
  fSigmax[8] = 0.04568;
  fOffsety[8] = -0.00004;
  fSigmay[8] = 0.00540;
  fOffsetx[9] = -0.00014;
  fSigmax[9] = 0.01820;
  fOffsety[9] = -0.00022;
  fSigmay[9] = 0.00569;
  fOffsetx[10] = 0.00046;
  fSigmax[10] = 0.01210;
  fOffsety[10] = -0.00034;
  fSigmay[10] = 0.00633;
  fOffsetx[11] = 0.0;
  fSigmax[11] = 0.0;
  fOffsety[11] = 0.0;
  fSigmay[11] = 0.0;
  fOffsetx[12] = 0.00174;
  fSigmax[12] = 0.02427;
  fOffsety[12] = -0.00008;
  fSigmay[12] = 0.00553;
  fOffsetx[13] = -0.00049;
  fSigmax[13] = 0.01887;
  fOffsety[13] = 0.00014;
  fSigmay[13] = 0.00611;
  fOffsetx[14] = 0.0;
  fSigmax[14] = 0.0;
  fOffsety[14] = 0.0;
  fSigmay[14] = 0.0;
  fOffsetx[15] = 0.0;
  fSigmax[15] = 0.0;
  fOffsety[15] = 0.0;
  fSigmay[15] = 0.0;
  fOffsetx[16] = 0.0;
  fSigmax[16] = 0.0;
  fOffsety[16] = 0.0;
  fSigmay[16] = 0.0;
  fOffsetx[17] = 0.0;
  fSigmax[17] = 0.0;
  fOffsety[17] = 0.0;
  fSigmay[17] = 0.0;
  fOffsetx[18] = 0.0;
  fSigmax[18] = 0.0;
  fOffsety[18] = 0.0;
  fSigmay[18] = 0.0;
  fOffsetx[19] = 0.00095;
  fSigmax[19] = 0.02936;
  fOffsety[19] = -0.00199;
  fSigmay[19] = 0.00725;

  fOffx[0] = 0.0;
  fSigx[0] = 0.0015;
  fOffy[0] = 0.0;
  fSigy[0] = 0.0015;

  fOffx[1] = 0.0;
  fSigx[1] = 0.0015;
  fOffy[1] = 0.0;
  fSigy[1] = 0.0015;

  fOffx[2] = 0.0;
  fSigx[2] = 0.0015;
  fOffy[2] = 0.0;
  fSigy[2] = 0.0015;

  fOffx[3] = 0.0;
  fSigx[3] = 0.0015;
  fOffy[3] = 0.0;
  fSigy[3] = 0.0015;

  // cuts on matching with primary vertex
  fOffx[4] = 0.0000002;
  fSigx[4] = 0.0003999;
  fOffy[4] = -0.0000040;
  fSigy[4] = 0.0000669;
  fOffx[5] = 0.0000138;
  fSigx[5] = 0.0002442;
  fOffy[5] = 0.0000001;
  fSigy[5] = 0.0000753;
  fOffx[6] = 0.0000059;
  fSigx[6] = 0.0002075;
  fOffy[6] = -0.0000009;
  fSigy[6] = 0.0000618;
  fOffx[7] = 0.0000125;
  fSigx[7] = 0.0002822;
  fOffy[7] = -0.0000017;
  fSigy[7] = 0.0000851;
  fOffx[8] = 0.0000020;
  fSigx[8] = 0.0003823;
  fOffy[8] = 0.0000048;
  fSigy[8] = 0.0000706;
  fOffx[9] = 0.000064;
  fSigx[9] = 0.00029;
  fOffy[9] = 0.0000122;
  fSigy[9] = 0.0000735;
  fOffx[10] = -0.0000020;
  fSigx[10] = 0.0001917;
  fOffy[10] = 0.0000028;
  fSigy[10] = 0.0000612;
  fOffx[11] = 0.0000073;
  fSigx[11] = 0.0003012;
  fOffy[11] = 0.0000032;
  fSigy[11] = 0.0000962;
  fOffx[12] = 0.0;
  fSigx[12] = 0.0;
  fOffy[12] = 0.0;
  fSigy[12] = 0.0;
  fOffx[13] = 0.0;
  fSigx[13] = 0.0;
  fOffy[13] = 0.0;
  fSigy[13] = 0.0;
  fOffx[14] = 0.0;
  fSigx[14] = 0.0;
  fOffy[14] = 0.0;
  fSigy[14] = 0.0;
  fOffx[15] = -0.00006;
  fSigx[15] = 0.0005;
  fOffy[15] = -0.0000206;
  fSigy[15] = 0.0001054;
  fOffx[16] = 0.0;
  fSigx[16] = 0.0;
  fOffy[16] = 0.0;
  fSigy[16] = 0.0;
  fOffx[17] = -0.0000051;
  fSigx[17] = 0.0002168;
  fOffy[17] = 0.0000167;
  fSigy[17] = 0.0000788;

  // extrapolated track cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = 0.046;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_SaleveDec2016_field(int /*run_id*/) {
  // dev cuts
    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_SaleveDec2016_field: setting offsets and sigmas";
    INFO(info); 
    
  fOffsetx[0] = 0.00079;
  fSigmax[0] = 0.01064;
  fOffsety[0] = 0.00002;
  fSigmay[0] = 0.00646;
  fOffsetx[1] = 0.00045;
  fSigmax[1] = 0.00963;
  fOffsety[1] = -0.00005;
  fSigmay[1] = 0.00619;
  fOffsetx[2] = 0.00009;
  fSigmax[2] = 0.01301;
  fOffsety[2] = -0.00051;
  fSigmay[2] = 0.00604;
  fOffsetx[3] = 0.00117;
  fSigmax[3] = 0.01356;
  fOffsety[3] = -0.00017;
  fSigmay[3] = 0.00590;
  fOffsetx[4] = 0.00028;
  fSigmax[4] = 0.01623;
  fOffsety[4] = -0.00050;
  fSigmay[4] = 0.00669;
  fOffsetx[5] = 0.00092;
  fSigmax[5] = 0.01706;
  fOffsety[5] = 0.00058;
  fSigmay[5] = 0.00734;
  fOffsetx[6] = -0.00038;
  fSigmax[6] = 0.01426;
  fOffsety[6] = -0.00072;
  fSigmay[6] = 0.00607;
  fOffsetx[7] = 0.0;
  fSigmax[7] = 0.015;
  fOffsety[7] = -0.00058;
  fSigmay[7] = 0.00663;
  fOffsetx[8] = 0.0;
  fSigmax[8] = 0.0;
  fOffsety[8] = 0.0;
  fSigmay[8] = 0.0;
  fOffsetx[9] = 0.0;
  fSigmax[9] = 0.015;
  fOffsety[9] = -0.00049;
  fSigmay[9] = 0.00719;
  fOffsetx[10] = 0.00106;
  fSigmax[10] = 0.01528;
  fOffsety[10] = -0.00031;
  fSigmay[10] = 0.00589;
  fOffsetx[11] = 0.0;
  fSigmax[11] = 0.0;
  fOffsety[11] = 0.0;
  fSigmay[11] = 0.0;
  fOffsetx[12] = -0.00135;
  fSigmax[12] = 0.02990;
  fOffsety[12] = 0.00017;
  fSigmay[12] = 0.00698;
  fOffsetx[13] = 0.00245;
  fSigmax[13] = 0.01290;
  fOffsety[13] = -0.00104;
  fSigmay[13] = 0.00625;
  fOffsetx[14] = 0.00740;
  fSigmax[14] = 0.00000;
  fOffsety[14] = 0.97827;
  fSigmay[14] = 0.18627;
  fOffsetx[15] = 0.01102;
  fSigmax[15] = 0.03589;
  fOffsety[15] = 0.00032;
  fSigmay[15] = 0.00956;

  fOffsetx[16] = 0.00691;
  fSigmax[16] = 0.03643;
  fOffsety[16] = 0.00232;
  fSigmay[16] = 0.00771;
  fOffsetx[17] = 0.0;
  fSigmax[17] = 0.0;
  fOffsety[17] = 0.0;
  fSigmay[17] = 0.0;
  fOffsetx[18] = 0.0;
  fSigmax[18] = 0.0;
  fOffsety[18] = 0.0;
  fSigmay[18] = 0.0;
  fOffsetx[19] = 0.00576;
  fSigmax[19] = 0.03554;
  fOffsety[19] = -0.00380;
  fSigmay[19] = 0.00807;

  // cuts: 0-3 are not use
  fOffx[0] = -0.004793;
  fSigx[0] = 0.003;
  fOffy[0] = 0.00747;
  fSigy[0] = 0.0026;

  fOffx[1] = -0.00427;
  fSigx[1] = 0.002;
  fOffy[1] = 0.00186;
  fSigy[1] = 0.0016;

  fOffx[2] = -0.00211;
  fSigx[2] = 0.0018;
  fOffy[2] = 0.00057;
  fSigy[2] = 0.00167;

  fOffx[3] = -0.000357;
  fSigx[3] = 0.0015;
  fOffy[3] = -0.004726;
  fSigy[3] = 0.00165;

  // cuts on matching with primary vertex
  fOffx[4] = 0.0000164;
  fSigx[4] = 0.0004559;
  fOffy[4] = -0.0000006;
  fSigy[4] = 0.0000928;
  fOffx[5] = -0.0000028;
  fSigx[5] = 0.0001825;
  fOffy[5] = -0.0000052;
  fSigy[5] = 0.0000672;
  fOffx[6] = 0.0000128;
  fSigx[6] = 0.0002172;
  fOffy[6] = -0.0000033;
  fSigy[6] = 0.0000762;
  fOffx[7] = 0.0000220;
  fSigx[7] = 0.0003412;
  fOffy[7] = -0.0000052;
  fSigy[7] = 0.0001009;
  fOffx[8] = 0.0000299;
  fSigx[8] = 0.0004655;
  fOffy[8] = -0.0000055;
  fSigy[8] = 0.0001009;
  fOffx[9] = -0.0000224;
  fSigx[9] = 0.0001773;
  fOffy[9] = 0.0000161;
  fSigy[9] = 0.0000662;
  fOffx[10] = -0.0000194;
  fSigx[10] = 0.0002103;
  fOffy[10] = 0.0000059;
  fSigy[10] = 0.0000768;
  fOffx[11] = 0.0000144;
  fSigx[11] = 0.0003053;
  fOffy[11] = 0.0000132;
  fSigy[11] = 0.0001221;
  fOffx[12] = 0.0000487;
  fSigx[12] = 0.0009069;
  fOffy[12] = 0.0000050;
  fSigy[12] = 0.0001081;
  fOffx[13] = -0.0000075;
  fSigx[13] = 0.0002071;
  fOffy[13] = 0.0000041;
  fSigy[13] = 0.0000783;
  fOffx[14] = 0.0000060;
  fSigx[14] = 0.0002434;
  fOffy[14] = -0.0000117;
  fSigy[14] = 0.0000799;
  fOffx[15] = 0.0000177;
  fSigx[15] = 0.0008047;
  fOffy[15] = 0.0000205;
  fSigy[15] = 0.0001162;
  fOffx[16] = -0.0000223;
  fSigx[16] = 0.0002392;
  fOffy[16] = 0.0000205;
  fSigy[16] = 0.0000881;
  fOffx[17] = 0.0000246;
  fSigx[17] = 0.0002412;
  fOffy[17] = 0.0000248;
  fSigy[17] = 0.0000937;

  // extrapolated tracks cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = -0.046;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_JuraOct2017_field(int /*run_id*/) {
  // dev cuts

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_JuraOct2017_field: setting offsets and sigmas";
    INFO(info); 
 
  /*
  fOffsetx[0] = 0.00051;    fSigmax[0] =  0.0125105;
  fOffsety[0] = -0.00006;    fSigmay[0] =  0.0052764;
  fOffsetx[1] = 0.00086;    fSigmax[1] =  0.0129962;
  fOffsety[1] = 0.00007;    fSigmay[1] =  0.0054475;
  fOffsetx[2] = 0.00002;    fSigmax[2] =  0.0122530;
  fOffsety[2] = 0.00004;    fSigmay[2] =  0.0052581;
  fOffsetx[3] = 0.00020;    fSigmax[3] =  0.0129346;
  fOffsety[3] = -0.00011;    fSigmay[3] =  0.0058659;
  fOffsetx[4] = 0.00064;    fSigmax[4] =  0.0150838;
  fOffsety[4] = -0.00018;    fSigmay[4] =  0.0058707;
  fOffsetx[5] = 0.00044;    fSigmax[5] =  0.0156745;
  fOffsety[5] = -0.00010;    fSigmay[5] =  0.0066797;
  */

  fOffsetx[0] = 0.00005;
  fSigmax[0] = 0.0119579;
  fOffsety[0] = -0.00016;
  fSigmay[0] = 0.0051942;
  fOffsetx[1] = 0.00024;
  fSigmax[1] = 0.0125211;
  fOffsety[1] = 0.00008;
  fSigmay[1] = 0.0053513;
  fOffsetx[2] = 0.00000;
  fSigmax[2] = 0.0119924;
  fOffsety[2] = 0.00012;
  fSigmay[2] = 0.0052590;
  fOffsetx[3] = 0.00035;
  fSigmax[3] = 0.0126723;
  fOffsety[3] = -0.00002;
  fSigmay[3] = 0.0058513;
  fOffsetx[4] = 0.00089;
  fSigmax[4] = 0.0151112;
  fOffsety[4] = 0.00009;
  fSigmay[4] = 0.0057940;
  fOffsetx[5] = 0.00120;
  fSigmax[5] = 0.0156643;
  fOffsety[5] = -0.00000;
  fSigmay[5] = 0.0066776;

  /* used before
  fOffsetx[6] = 0.00004;     fSigmax[6] =  0.0124127;
  fOffsety[6] = 0.00001;     fSigmay[6] =  0.0054054;
  fOffsetx[7] = 0.01539;    fSigmax[7] =  0.1000000;
  fOffsety[7] = 0.00224;    fSigmay[7] =  0.0065496;
  fOffsetx[8] = 0.01539;    fSigmax[8] =  0.1000000;
  fOffsety[8] = 0.00224;    fSigmay[8] =  0.0065496;
  fOffsetx[9] = -0.00393;    fSigmax[9] =  0.0392162;
  fOffsety[9] = -0.00023;    fSigmay[9] =  0.0056589;
  fOffsetx[10] = 0.00025;    fSigmax[10] =  0.0132427;
  fOffsety[10] = -0.00007;    fSigmay[10] =  0.0060167;
  fOffsetx[11] = 0.01539;    fSigmax[11] =  0.1000000;
  fOffsety[11] = 0.00224;    fSigmay[11] =  0.0065496;
  fOffsetx[12] = 0.00073;    fSigmax[12] =  0.0115300;
  fOffsety[12] = -0.00049;    fSigmay[12] =  0.0051311;
  fOffsetx[13] = -0.00124;    fSigmax[13] =  0.0545703;
  fOffsety[13] = -0.00020;    fSigmay[13] =  0.0072022;
  fOffsetx[14] = 0.00200;    fSigmax[14] =  0.0314421;
  fOffsety[14] = -0.00026;    fSigmay[14] =  0.0077302;
  fOffsetx[15] = -0.0;       fSigmax[15] =  0.03;
  fOffsety[15] = -0.0;       fSigmay[15] =  0.007;
  fOffsetx[16] = 0.00245;    fSigmax[16] =  0.0330714;
  fOffsety[16] = -0.00026;    fSigmay[16] =  0.0069557;
  fOffsetx[17] = 0.00343;    fSigmax[17] =  0.0305652;
  fOffsety[17] = -0.00004;    fSigmay[17] =  0.0080821;
  fOffsetx[18] = 0.0    ;    fSigmax[18] =  0.03;
  fOffsety[18] = 0.0    ;    fSigmay[18] =  0.007;
  fOffsetx[19] = 0.00158;    fSigmax[19] =  0.0325604;
  fOffsety[19] = 0.00059;    fSigmay[19] =  0.0070290;
  */

  // added on 5.12.2018 (based on run 1416 with --effic=1)
  fOffsetx[6] = 0.00001;
  fSigmax[6] = 0.0122519;
  fOffsety[6] = 0.00013;
  fSigmay[6] = 0.0053732;
  fOffsetx[7] = 0.00020;
  fSigmax[7] = 0.0232666;
  fOffsety[7] = 0.00035;
  fSigmay[7] = 0.0066678;
  fOffsetx[8] = -0.00064;
  fSigmax[8] = 0.0227686;
  fOffsety[8] = -0.00016;
  fSigmay[8] = 0.0062499;
  fOffsetx[9] = 0.00013;
  fSigmax[9] = 0.0127437;
  fOffsety[9] = -0.00014;
  fSigmay[9] = 0.0055400;
  fOffsetx[10] = 0.00036;
  fSigmax[10] = 0.0128984;
  fOffsety[10] = -0.00002;
  fSigmay[10] = 0.0059792;
  fOffsetx[11] = 0.00042;
  fSigmax[11] = 0.0240284;
  fOffsety[11] = -0.00051;
  fSigmay[11] = 0.0072616;
  fOffsetx[12] = -0.00066;
  fSigmax[12] = 0.0224608;
  fOffsety[12] = -0.00011;
  fSigmay[12] = 0.0056808;
  fOffsetx[13] = 0.00021;
  fSigmax[13] = 0.0133975;
  fOffsety[13] = 0.00016;
  fSigmay[13] = 0.0057623;
  fOffsetx[14] = 0.00092;
  fSigmax[14] = 0.0154056;
  fOffsety[14] = 0.00006;
  fSigmay[14] = 0.0059983;
  fOffsetx[15] = 0.00114;
  fSigmax[15] = 0.0281642;
  fOffsety[15] = -0.00008;
  fSigmay[15] = 0.0070190;
  fOffsetx[16] = 0.00135;
  fSigmax[16] = 0.0295323;
  fOffsety[16] = -0.00003;
  fSigmay[16] = 0.0072097;
  fOffsetx[17] = 0.00127;
  fSigmax[17] = 0.0159710;
  fOffsety[17] = 0.00003;
  fSigmay[17] = 0.0068877;
  fOffsetx[18] = 0.00127;
  fSigmax[18] = 0.0284584;
  fOffsety[18] = -0.00003;
  fSigmay[18] = 0.0079667;
  fOffsetx[19] = 0.00107;
  fSigmax[19] = 0.0300292;
  fOffsety[19] = 0.00038;
  fSigmay[19] = 0.0071514;

  // cuts on matching - mostly with primary angle. 0-3 are not used (4 hit tracks)
  fOffx[0] = 0.0;
  fSigx[0] = 0.0015;
  fOffy[0] = 0.0;
  fSigy[0] = 0.0015;

  fOffx[1] = 0.0;
  fSigx[1] = 0.0015;
  fOffy[1] = 0.0;
  fSigy[1] = 0.0015;

  fOffx[2] = 0.0;
  fSigx[2] = 0.0015;
  fOffy[2] = 0.0;
  fSigy[2] = 0.0015;

  fOffx[3] = 0.0;
  fSigx[3] = 0.0015;
  fOffy[3] = 0.0;
  fSigy[3] = 0.0015;

  /*
  fOffx[4] = -0.000006;    fSigx[4] =  0.0004279;
  fOffy[4] = 0.000002;     fSigy[4] =  0.0000665;
  fOffx[5] = 0.0000245;    fSigx[5] =  0.0002937;
  fOffy[5] = -0.000014;    fSigy[5] =  0.0000798;
  fOffx[6] = 0.000035;    fSigx[6] =  0.0003745;
  fOffy[6] = 0.000000;    fSigy[6] =  0.0000622;
  fOffx[7] = 0.000039;    fSigx[7] =  0.0005336;
  fOffy[7] = 0.000001;    fSigy[7] =  0.0000859;
  fOffx[8] = -0.000001;    fSigx[8] =  0.0004522;
  fOffy[8] = 0.000004;    fSigy[8] =  0.0000696;
  fOffx[9] = 0.000053;    fSigx[9] =  0.0003072;
  fOffy[9] = 0.000034;    fSigy[9] =  0.0000872;
  fOffx[10] = 0.000023;    fSigx[10] =  0.0002674;
  fOffy[10] = 0.000002;    fSigy[10] =  0.0000555;
  fOffx[11] = 0.000035;    fSigx[11] =  0.0005618;
  fOffy[11] = 0.000004;    fSigy[11] =  0.0001057;
  fOffx[12] = -0.000100;    fSigx[12] =  0.0010234;
  fOffy[12] = -0.000004;    fSigy[12] =  0.0000900;
  fOffx[13] = 0.000019;    fSigx[13] =  0.0002918;
  fOffy[13] = -0.000026;    fSigy[13] =  0.0000801;
  fOffx[14] = -0.000005;    fSigx[14] =  0.0003049;
  fOffy[14] = -0.000004;    fSigy[14] =  0.0000689;
  fOffx[15] = -0.000137;    fSigx[15] =  0.0009280;
  fOffy[15] = -0.000000;    fSigy[15] =  0.0000863;
  fOffx[16] = -0.000002;    fSigx[16] =  0.0003396;
  fOffy[16] = 0.000009;    fSigy[16] =  0.0000919;
  fOffx[17] = -0.000009;    fSigx[17] =  0.0003040;
  fOffy[17] = -0.000008;    fSigy[17] =  0.0000701;
  */

  // added on 5.12.2018 (based on run 1416 with --effic=1)
  fOffx[4] = 0.000004;
  fSigx[4] = 0.0005086;
  fOffy[4] = 0.000002;
  fSigy[4] = 0.0000816;
  fOffx[5] = 0.000057;
  fSigx[5] = 0.0003873;
  fOffy[5] = -0.000001;
  fSigy[5] = 0.0000553;
  fOffx[6] = 0.000041;
  fSigx[6] = 0.0003862;
  fOffy[6] = -0.000001;
  fSigy[6] = 0.0000757;
  fOffx[7] = 0.000027;
  fSigx[7] = 0.0004324;
  fOffy[7] = -0.000004;
  fSigy[7] = 0.0000910;
  fOffx[8] = 0.000028;
  fSigx[8] = 0.0005404;
  fOffy[8] = 0.000005;
  fSigy[8] = 0.0000852;
  fOffx[9] = 0.000053;
  fSigx[9] = 0.0004066;
  fOffy[9] = 0.000003;
  fSigy[9] = 0.0000560;
  fOffx[10] = 0.000040;
  fSigx[10] = 0.0003586;
  fOffy[10] = 0.000003;
  fSigy[10] = 0.0000670;
  fOffx[11] = 0.000032;
  fSigx[11] = 0.0004453;
  fOffy[11] = 0.000003;
  fSigy[11] = 0.0000939;
  fOffx[12] = 0.000040;
  fSigx[12] = 0.0006226;
  fOffy[12] = -0.000000;
  fSigy[12] = 0.0001123;
  fOffx[13] = 0.000036;
  fSigx[13] = 0.0004220;
  fOffy[13] = -0.000002;
  fSigy[13] = 0.0000563;
  fOffx[14] = 0.000024;
  fSigx[14] = 0.0004388;
  fOffy[14] = -0.000004;
  fSigy[14] = 0.0000816;
  fOffx[15] = 0.000043;
  fSigx[15] = 0.0006686;
  fOffy[15] = -0.000002;
  fSigy[15] = 0.0001125;
  fOffx[16] = 0.000033;
  fSigx[16] = 0.0004312;
  fOffy[16] = 0.000001;
  fSigy[16] = 0.0000590;
  fOffx[17] = 0.000017;
  fSigx[17] = 0.0004416;
  fOffy[17] = 0.000002;
  fSigy[17] = 0.0000836;

  // extrapolated track cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = 0.046;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_JuraOct2017_field_xela75(int /*run_id*/) {
  // dev cuts

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_JuraOct2017_field_xela75: setting offsets and sigmas";
    INFO(info); 

	fOffsetx[0] = -0.00013;    fSigmax[0] =  0.0076985;
	fOffsety[0] = -0.00006;    fSigmay[0] =  0.0052700;
	fOffsetx[1] = -0.00001;    fSigmax[1] =  0.0078553;
	fOffsety[1] = -0.00002;    fSigmay[1] =  0.0053686;
	fOffsetx[2] = 0.00022;    fSigmax[2] =  0.0079318;
	fOffsety[2] = -0.00007;    fSigmay[2] =  0.0052580;
	fOffsetx[3] = 0.00032;    fSigmax[3] =  0.0080787;
	fOffsety[3] = -0.00015;    fSigmay[3] =  0.0058277;
	fOffsetx[4] = 0.00052;    fSigmax[4] =  0.0088404;
	fOffsety[4] = -0.00007;    fSigmay[4] =  0.0057859;
	fOffsetx[5] = 0.00048;    fSigmax[5] =  0.0094937;
	fOffsety[5] = -0.00001;    fSigmay[5] =  0.0066601;


	//based on 33453, 33480, 33481
	//eff=1

	fOffsetx[6] = 0.00025;    fSigmax[6] =  0.0080397;
	fOffsety[6] = -0.00006;    fSigmay[6] =  0.0053183;
	fOffsetx[7] = 0.00029;    fSigmax[7] =  0.0135319;
	fOffsety[7] = -0.00004;    fSigmay[7] =  0.0062076;
	fOffsetx[8] = 0.00045;    fSigmax[8] =  0.0130302;
	fOffsety[8] = -0.00009;    fSigmay[8] =  0.0057612;
	fOffsetx[9] = -0.00013;    fSigmax[9] =  0.0080553;
	fOffsety[9] = -0.00006;    fSigmay[9] =  0.0053464;
	fOffsetx[10] = 0.00033;    fSigmax[10] =  0.0082170;
	fOffsety[10] = -0.00015;    fSigmay[10] =  0.0059120;
	fOffsetx[11] = 0.00050;    fSigmax[11] =  0.0140662;
	fOffsety[11] = -0.00022;    fSigmay[11] =  0.0068467;
	fOffsetx[12] = 0.00020;    fSigmax[12] =  0.0131346;
	fOffsety[12] = -0.00018;    fSigmay[12] =  0.0058293;
	fOffsetx[13] = -0.00003;    fSigmax[13] =  0.0081889;
	fOffsety[13] = -0.00000;    fSigmay[13] =  0.0054522;
	fOffsetx[14] = 0.00052;    fSigmax[14] =  0.0090444;
	fOffsety[14] = -0.00006;    fSigmay[14] =  0.0057281;
	fOffsetx[15] = 0.00079;    fSigmax[15] =  0.0151246;
	fOffsety[15] = -0.00012;    fSigmay[15] =  0.0067096;
	fOffsetx[16] = 0.00059;    fSigmax[16] =  0.0161309;
	fOffsety[16] = -0.00004;    fSigmay[16] =  0.0068553;
	fOffsetx[17] = 0.00051;    fSigmax[17] =  0.0096534;
	fOffsety[17] = 0.00001;    fSigmay[17] =  0.0065872;
	fOffsetx[18] = 0.00084;    fSigmax[18] =  0.0155045;
	fOffsety[18] = 0.00008;    fSigmay[18] =  0.0077070;
	fOffsetx[19] = 0.00090;    fSigmax[19] =  0.0161542;
	fOffsety[19] = 0.00038;    fSigmay[19] =  0.0069156;



  // cuts on matching - mostly with primary angle. 0-3 are not used (4 hit tracks)
  fOffx[0] = 0.0;
  fSigx[0] = 0.0015;
  fOffy[0] = 0.0;
  fSigy[0] = 0.0015;

  fOffx[1] = 0.0;
  fSigx[1] = 0.0015;
  fOffy[1] = 0.0;
  fSigy[1] = 0.0015;

  fOffx[2] = 0.0;
  fSigx[2] = 0.0015;
  fOffy[2] = 0.0;
  fSigy[2] = 0.0015;

  fOffx[3] = 0.0;
  fSigx[3] = 0.0015;
  fOffy[3] = 0.0;
  fSigy[3] = 0.0015;

 
//based on 33453, 33480, 33481
	//eff=1
	fOffx[4] = 0.000018;    fSigx[4] =  0.0004797;
	fOffy[4] = 0.000001;    fSigy[4] =  0.0000801;
	fOffx[5] = 0.000032;    fSigx[5] =  0.0003022;
	fOffy[5] = -0.000001;    fSigy[5] =  0.0000531;
	fOffx[6] = 0.000014;    fSigx[6] =  0.0003005;
	fOffy[6] = -0.000001;    fSigy[6] =  0.0000717;
	fOffx[7] = 0.000010;    fSigx[7] =  0.0003356;
	fOffy[7] = -0.000002;    fSigy[7] =  0.0000872;
	fOffx[8] = 0.000037;    fSigx[8] =  0.0005087;
	fOffy[8] = 0.000004;    fSigy[8] =  0.0000819;
	fOffx[9] = 0.000041;    fSigx[9] =  0.0003166;
	fOffy[9] = 0.000003;    fSigy[9] =  0.0000551;
	fOffx[10] = 0.000025;    fSigx[10] =  0.0003039;
	fOffy[10] = 0.000003;    fSigy[10] =  0.0000723;
	fOffx[11] = 0.000019;    fSigx[11] =  0.0003401;
	fOffy[11] = 0.000003;    fSigy[11] =  0.0000878;
	fOffx[12] = 0.000037;    fSigx[12] =  0.0005769;
	fOffy[12] = -0.000002;    fSigy[12] =  0.0001065;
	fOffx[13] = 0.000022;    fSigx[13] =  0.0003243;
	fOffy[13] = -0.000003;    fSigy[13] =  0.0000578;
	fOffx[14] = 0.000011;    fSigx[14] =  0.0003328;
	fOffy[14] = -0.000005;    fSigy[14] =  0.0000821;
	fOffx[15] = 0.000023;    fSigx[15] =  0.0006223;
	fOffy[15] = -0.000002;    fSigy[15] =  0.0001067;
	fOffx[16] = 0.000020;    fSigx[16] =  0.0003423;
	fOffy[16] = 0.000000;    fSigy[16] =  0.0000599;
	fOffx[17] = 0.000009;    fSigx[17] =  0.0003351;
	fOffy[17] = 0.000001;    fSigy[17] =  0.0000826;



//DATA BELOW IS NOT USE 

  // extrapolated track cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = 0.046;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_JuraOct2017_field_xela40(int /*run_id*/) {
  // dev cuts

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_JuraOct2017_field_xela40: setting offsets and sigmas";
    INFO(info); 

	fOffsetx[0] = -0.00064;    fSigmax[0] =  0.0061088;
	fOffsety[0] = 0.00032;    fSigmay[0] =  0.0053099;
	fOffsetx[1] = -0.00038;    fSigmax[1] =  0.0062452;
	fOffsety[1] = 0.00010;    fSigmay[1] =  0.0054542;
	fOffsetx[2] = 0.00125;    fSigmax[2] =  0.0066680;
	fOffsety[2] = 0.00013;    fSigmay[2] =  0.0052974;
	fOffsetx[3] = 0.00110;    fSigmax[3] =  0.0067922;
	fOffsety[3] = -0.00055;    fSigmay[3] =  0.0059882;
	fOffsetx[4] = -0.00081;    fSigmax[4] =  0.0068503;
	fOffsety[4] = -0.00051;    fSigmay[4] =  0.0057768;
	fOffsetx[5] = -0.00109;    fSigmax[5] =  0.0076321;
	fOffsety[5] = -0.00020;    fSigmay[5] =  0.0066762;



	//based on 35206, 35229, 35264
	//eff=1

	fOffsetx[6] = 0.00121;    fSigmax[6] =  0.0066873;
	fOffsety[6] = 0.00013;    fSigmay[6] =  0.0053301;
	fOffsetx[7] = -0.00090;    fSigmax[7] =  0.0090248;
	fOffsety[7] = 0.00048;    fSigmay[7] =  0.0063589;
	fOffsetx[8] = 0.00027;    fSigmax[8] =  0.0080700;
	fOffsety[8] = 0.00060;    fSigmay[8] =  0.0057702;
	fOffsetx[9] = -0.00065;    fSigmax[9] =  0.0063557;
	fOffsety[9] = 0.00037;    fSigmay[9] =  0.0054542;
	fOffsetx[10] = 0.00110;    fSigmax[10] =  0.0068821;
	fOffsety[10] = -0.00056;    fSigmay[10] =  0.0060279;
	fOffsetx[11] = 0.00084;    fSigmax[11] =  0.0088966;
	fOffsety[11] = -0.00087;    fSigmay[11] =  0.0067519;
	fOffsetx[12] = 0.00015;    fSigmax[12] =  0.0075607;
	fOffsety[12] = -0.00047;    fSigmay[12] =  0.0058814;
	fOffsetx[13] = -0.00036;    fSigmax[13] =  0.0064730;
	fOffsety[13] = 0.00013;    fSigmay[13] =  0.0056036;
	fOffsetx[14] = -0.00079;    fSigmax[14] =  0.0068950;
	fOffsety[14] = -0.00050;    fSigmay[14] =  0.0057880;
	fOffsetx[15] = -0.00125;    fSigmax[15] =  0.0094172;
	fOffsety[15] = -0.00044;    fSigmay[15] =  0.0067748;
	fOffsetx[16] = -0.00121;    fSigmax[16] =  0.0095455;
	fOffsety[16] = 0.00025;    fSigmay[16] =  0.0069499;
	fOffsetx[17] = -0.00108;    fSigmax[17] =  0.0076726;
	fOffsety[17] = -0.00020;    fSigmay[17] =  0.0067047;
	fOffsetx[18] = -0.00161;    fSigmax[18] =  0.0100365;
	fOffsety[18] = -0.00014;    fSigmay[18] =  0.0078137;
	fOffsetx[19] = -0.00100;    fSigmax[19] =  0.0096044;
	fOffsety[19] = 0.00039;    fSigmay[19] =  0.0071199;





  // cuts on matching - mostly with primary angle. 0-3 are not used (4 hit tracks)
  fOffx[0] = 0.0;
  fSigx[0] = 0.0015;
  fOffy[0] = 0.0;
  fSigy[0] = 0.0015;

  fOffx[1] = 0.0;
  fSigx[1] = 0.0015;
  fOffy[1] = 0.0;
  fSigy[1] = 0.0015;

  fOffx[2] = 0.0;
  fSigx[2] = 0.0015;
  fOffy[2] = 0.0;
  fSigy[2] = 0.0015;

  fOffx[3] = 0.0;
  fSigx[3] = 0.0015;
  fOffy[3] = 0.0;
  fSigy[3] = 0.0015;

 
	//based on 35206, 35229, 35264
	//eff=1
	fOffx[4] = 0.000079;    fSigx[4] =  0.0004961;
	fOffy[4] = -0.000007;    fSigy[4] =  0.0000920;
	fOffx[5] = 0.000040;    fSigx[5] =  0.0002921;
	fOffy[5] = 0.000001;    fSigy[5] =  0.0000536;
	fOffx[6] = 0.000001;    fSigx[6] =  0.0002765;
	fOffy[6] = 0.000005;    fSigy[6] =  0.0000729;
	fOffx[7] = -0.000020;    fSigx[7] =  0.0003166;
	fOffy[7] = 0.000003;    fSigy[7] =  0.0000870;
	fOffx[8] = 0.000086;    fSigx[8] =  0.0005209;
	fOffy[8] = 0.000005;    fSigy[8] =  0.0000951;
	fOffx[9] = 0.000038;    fSigx[9] =  0.0002993;
	fOffy[9] = 0.000002;    fSigy[9] =  0.0000550;
	fOffx[10] = 0.000007;    fSigx[10] =  0.0002802;
	fOffy[10] = 0.000003;    fSigy[10] =  0.0000737;
	fOffx[11] = -0.000014;    fSigx[11] =  0.0003206;
	fOffy[11] = 0.000003;    fSigy[11] =  0.0000871;
	fOffx[12] = -0.000040;    fSigx[12] =  0.0005815;
	fOffy[12] = -0.000007;    fSigy[12] =  0.0001096;
	fOffx[13] = -0.000030;    fSigx[13] =  0.0003075;
	fOffy[13] = -0.000003;    fSigy[13] =  0.0000551;
	fOffx[14] = -0.000029;    fSigx[14] =  0.0003119;
	fOffy[14] = 0.000000;    fSigy[14] =  0.0000844;
	fOffx[15] = -0.000074;    fSigx[15] =  0.0006267;
	fOffy[15] = -0.000003;    fSigy[15] =  0.0001112;
	fOffx[16] = -0.000047;    fSigx[16] =  0.0003250;
	fOffy[16] = -0.000000;    fSigy[16] =  0.0000574;
	fOffx[17] = -0.000037;    fSigx[17] =  0.0003154;
	fOffy[17] = 0.000001;    fSigy[17] =  0.0000863;




//DATA BELOW IS NOT USE 

  // extrapolated track cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = 0.046;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_JuraNov2018_field(int /*run_id*/) {
  // dev cuts

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_JuraNov2018_field: setting offsets and sigmas";
    INFO(info); 

  
  // based on runs 38943, 38945
    fOffsetx[0] = -0.00047;    fSigmax[0] =  0.0121721;
	fOffsety[0] = -0.00007;    fSigmay[0] =  0.0051208;
	fOffsetx[1] = -0.00056;    fSigmax[1] =  0.0123533;
	fOffsety[1] = -0.00028;    fSigmay[1] =  0.0052407;
	fOffsetx[2] = 0.00072;    fSigmax[2] =  0.0120619;
	fOffsety[2] = -0.00011;    fSigmay[2] =  0.0050356;
	fOffsetx[3] = 0.00077;    fSigmax[3] =  0.0121453;
	fOffsety[3] = -0.00018;    fSigmay[3] =  0.0054774;
	fOffsetx[4] = 0.00146;    fSigmax[4] =  0.0154854;
	fOffsety[4] = -0.00012;    fSigmay[4] =  0.0056163;
	fOffsetx[5] = 0.00176;    fSigmax[5] =  0.0157114;
	fOffsety[5] = 0.00044;    fSigmay[5] =  0.0064464;



  

  // based on runs 38943, 38945
	fOffsetx[6] = 0.00071;    fSigmax[6] =  0.0124913;
	fOffsety[6] = -0.00012;    fSigmay[6] =  0.0051762;
	fOffsetx[7] = 0.00037;    fSigmax[7] =  0.0252778;
	fOffsety[7] = -0.00023;    fSigmay[7] =  0.0062707;
	fOffsetx[8] = 0.00049;    fSigmax[8] =  0.0242578;
	fOffsety[8] = -0.00006;    fSigmay[8] =  0.0059340;
	fOffsetx[9] = -0.00052;    fSigmax[9] =  0.0131763;
	fOffsety[9] = -0.00005;    fSigmay[9] =  0.0053393;
	fOffsetx[10] = 0.00078;    fSigmax[10] =  0.0126347;
	fOffsety[10] = -0.00018;    fSigmay[10] =  0.0056236;
	fOffsetx[11] = 0.00078;    fSigmax[11] =  0.0249693;
	fOffsety[11] = -0.00005;    fSigmay[11] =  0.0067428;
	fOffsetx[12] = -0.00021;    fSigmax[12] =  0.0240501;
	fOffsety[12] = -0.00026;    fSigmay[12] =  0.0060168;
	fOffsetx[13] = -0.00053;    fSigmax[13] =  0.0133533;
	fOffsety[13] = -0.00032;    fSigmay[13] =  0.0054808;
	fOffsetx[14] = 0.00149;    fSigmax[14] =  0.0157453;
	fOffsety[14] = -0.00007;    fSigmay[14] =  0.0057566;
	fOffsetx[15] = 0.00159;    fSigmax[15] =  0.0289911;
	fOffsety[15] = -0.00010;    fSigmay[15] =  0.0068231;
	fOffsetx[16] = 0.00122;    fSigmax[16] =  0.0308354;
	fOffsety[16] = -0.00029;    fSigmay[16] =  0.0070663;
	fOffsetx[17] = 0.00173;    fSigmax[17] =  0.0159722;
	fOffsety[17] = 0.00043;    fSigmay[17] =  0.0065490;
	fOffsetx[18] = 0.00189;    fSigmax[18] =  0.0288911;
	fOffsety[18] = 0.00043;    fSigmay[18] =  0.0076184;
	fOffsetx[19] = 0.00117;    fSigmax[19] =  0.0312012;
	fOffsety[19] = 0.00003;    fSigmay[19] =  0.0071397;

	
  // cuts on matching - mostly with primary angle
  fOffx[0] = 0.0;
  fSigx[0] = 0.0015;
  fOffy[0] = 0.0;
  fSigy[0] = 0.0015;

  fOffx[1] = 0.0;
  fSigx[1] = 0.0015;
  fOffy[1] = 0.0;
  fSigy[1] = 0.0015;

  fOffx[2] = 0.0;
  fSigx[2] = 0.0015;
  fOffy[2] = 0.0;
  fSigy[2] = 0.0015;

  fOffx[3] = 0.0;
  fSigx[3] = 0.0015;
  fOffy[3] = 0.0;
  fSigy[3] = 0.0015;

	fOffx[4] = 0.000053;    fSigx[4] =  0.0005136;
	fOffy[4] = -0.000005;    fSigy[4] =  0.0000844;
	fOffx[5] = 0.000053;    fSigx[5] =  0.0004032;
	fOffy[5] = -0.000005;    fSigy[5] =  0.0000587;
	fOffx[6] = 0.000013;    fSigx[6] =  0.0004141;
	fOffy[6] = -0.000005;    fSigy[6] =  0.0000772;
	fOffx[7] = 0.000004;    fSigx[7] =  0.0004488;
	fOffy[7] = -0.000007;    fSigy[7] =  0.0000941;
	fOffx[8] = 0.000054;    fSigx[8] =  0.0005418;
	fOffy[8] = 0.000007;    fSigy[8] =  0.0000860;
	fOffx[9] = 0.000051;    fSigx[9] =  0.0004109;
	fOffy[9] = 0.000007;    fSigy[9] =  0.0000597;
	fOffx[10] = 0.000005;    fSigx[10] =  0.0004172;
	fOffy[10] = 0.000004;    fSigy[10] =  0.0000777;
	fOffx[11] = -0.000010;    fSigx[11] =  0.0004530;
	fOffy[11] = 0.000005;    fSigy[11] =  0.0000945;
	fOffx[12] = 0.000084;    fSigx[12] =  0.0006259;
	fOffy[12] = -0.000004;    fSigy[12] =  0.0001129;
	fOffx[13] = 0.000047;    fSigx[13] =  0.0004359;
	fOffy[13] = -0.000006;    fSigy[13] =  0.0000639;
	fOffx[14] = 0.000017;    fSigx[14] =  0.0004662;
	fOffy[14] = -0.000008;    fSigy[14] =  0.0000885;
	fOffx[15] = 0.000092;    fSigx[15] =  0.0006726;
	fOffy[15] = 0.000007;    fSigy[15] =  0.0001139;
	fOffx[16] = 0.000036;    fSigx[16] =  0.0004480;
	fOffy[16] = 0.000007;    fSigy[16] =  0.0000656;
	fOffx[17] = -0.000005;    fSigx[17] =  0.0004687;
	fOffy[17] = 0.000004;    fSigy[17] =  0.0000889;

  // for(int i=4;i<18;i++){
  // fOffx[i] = 0.;  fSigx[i] = 0.00012;
  // fOffy[i] = 0;   fSigy[i] = 0.00012;
  //}

  // extrapolated track cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = 0.046;
}






//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_SaleveOct2017_field(int /*run_id*/) {
    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_SaleveOct2017_field: setting offsets and sigmas";
    INFO(info); 
    
  // dev cuts

  /*
  fOffsetx[0] = 0.00075;    fSigmax[0] =  0.0120180;
  fOffsety[0] = -0.00005;    fSigmay[0] =  0.0054471;
  fOffsetx[1] = 0.00080;    fSigmax[1] =  0.0122806;
  fOffsety[1] = 0.00047;    fSigmay[1] =  0.0055119;
  fOffsetx[2] = -0.00026;    fSigmax[2] =  0.0130468;
  fOffsety[2] = 0.00037;    fSigmay[2] =  0.0056783;
  fOffsetx[3] = -0.00006;    fSigmax[3] =  0.0132904;
  fOffsety[3] = 0.00013;    fSigmay[3] =  0.0057786;
  fOffsetx[4] = 0.00021;    fSigmax[4] =  0.0147948;
  fOffsety[4] = -0.00025;    fSigmay[4] =  0.0064864;
  fOffsetx[5] = 0.00015;    fSigmax[5] =  0.0148384;
  fOffsety[5] = -0.00003;    fSigmay[5] =  0.0067930;
  */
  fOffsetx[0] = 0.00083;
  fSigmax[0] = 0.0121259;
  fOffsety[0] = 0.00001;
  fSigmay[0] = 0.0053444;
  fOffsetx[1] = 0.00085;
  fSigmax[1] = 0.0121875;
  fOffsety[1] = 0.00029;
  fSigmay[1] = 0.0055506;
  fOffsetx[2] = -0.00026;
  fSigmax[2] = 0.0132520;
  fOffsety[2] = 0.00034;
  fSigmay[2] = 0.0057336;
  fOffsetx[3] = -0.00018;
  fSigmax[3] = 0.0134057;
  fOffsety[3] = 0.00009;
  fSigmay[3] = 0.0057581;
  fOffsetx[4] = 0.00024;
  fSigmax[4] = 0.0149585;
  fOffsety[4] = -0.00027;
  fSigmay[4] = 0.0065166;
  fOffsetx[5] = 0.00013;
  fSigmax[5] = 0.0152758;
  fOffsety[5] = -0.00036;
  fSigmay[5] = 0.0067251;

  /*
  fOffsetx[6] = -0.00081;    fSigmax[6] =  0.0125783;
  fOffsety[6] = 0.00056;    fSigmay[6] =  0.0056776;
  fOffsetx[7] = 0.01539;    fSigmax[7] =  0.1000000;
  fOffsety[7] = 0.00224;    fSigmay[7] =  0.0065496;
  fOffsetx[8] = 0.00327;    fSigmax[8] =  0.0213047;
  fOffsety[8] = 0.00007;    fSigmay[8] =  0.0059202;
  fOffsetx[9] = 0.00127;    fSigmax[9] =  0.0070579;
  fOffsety[9] = -0.00024;    fSigmay[9] =  0.0072694;
  fOffsetx[10] = -0.00055;    fSigmax[10] =  0.0130693;
  fOffsety[10] = 0.00030;    fSigmay[10] =  0.0058341;
  fOffsetx[11] = 0.01539;    fSigmax[11] =  0.1000000;
  fOffsety[11] = 0.00224;    fSigmay[11] =  0.0065496;
  fOffsetx[12] = 0.00318;    fSigmax[12] =  0.0226907;
  fOffsety[12] = 0.00101;    fSigmay[12] =  0.0058905;
  fOffsetx[13] = 0.00031;    fSigmax[13] =  0.0493463;
  fOffsety[13] = 0.00094;    fSigmay[13] =  0.0082138;
  fOffsetx[14] = -0.00031;    fSigmax[14] =  0.0431165;
  fOffsety[14] = -0.00050;    fSigmay[14] =  0.0094946;
  fOffsetx[15] = -0.0;       fSigmax[15] =  0.03;
  fOffsety[15] = 0.0;        fSigmay[15] =  0.007;
  fOffsetx[16] = 0.00134;    fSigmax[16] =  0.0330440;
  fOffsety[16] = 0.00074;    fSigmay[16] =  0.0074459;
  fOffsetx[17] = -0.00218;    fSigmax[17] =  0.0316778;
  fOffsety[17] = -0.00027;    fSigmay[17] =  0.0077238;
  fOffsetx[18] = -0.0;        fSigmax[18] =  0.03;
  fOffsety[18] = -0.0;        fSigmay[18] =  0.007;
  fOffsetx[19] = 0.00086;     fSigmax[19] =  0.0311639;
  fOffsety[19] = -0.00103;    fSigmay[19] =  0.0072396;
  */
  fOffsetx[6] = -0.00025;
  fSigmax[6] = 0.0134590;
  fOffsety[6] = 0.00035;
  fSigmay[6] = 0.0058272;
  fOffsetx[7] = 0.00168;
  fSigmax[7] = 0.0257715;
  fOffsety[7] = -0.00009;
  fSigmay[7] = 0.0071964;
  fOffsetx[8] = 0.00224;
  fSigmax[8] = 0.0243067;
  fOffsety[8] = -0.00005;
  fSigmay[8] = 0.0062961;
  fOffsetx[9] = 0.00076;
  fSigmax[9] = 0.0128543;
  fOffsety[9] = -0.00003;
  fSigmay[9] = 0.0052777;
  fOffsetx[10] = -0.00019;
  fSigmax[10] = 0.0136464;
  fOffsety[10] = 0.00008;
  fSigmay[10] = 0.0058561;
  fOffsetx[11] = 0.00175;
  fSigmax[11] = 0.0253464;
  fOffsety[11] = 0.00008;
  fSigmay[11] = 0.0073436;
  fOffsetx[12] = 0.00177;
  fSigmax[12] = 0.0243253;
  fOffsety[12] = 0.00048;
  fSigmay[12] = 0.0063641;
  fOffsetx[13] = 0.00083;
  fSigmax[13] = 0.0128290;
  fOffsety[13] = 0.00036;
  fSigmay[13] = 0.0058198;
  fOffsetx[14] = 0.00020;
  fSigmax[14] = 0.0154548;
  fOffsety[14] = -0.00025;
  fSigmay[14] = 0.0067878;
  fOffsetx[15] = 0.00043;
  fSigmax[15] = 0.0274352;
  fOffsety[15] = -0.00007;
  fSigmay[15] = 0.0079939;
  fOffsetx[16] = 0.00083;
  fSigmax[16] = 0.0296477;
  fOffsety[16] = 0.00030;
  fSigmay[16] = 0.0076491;
  fOffsetx[17] = 0.00013;
  fSigmax[17] = 0.0155037;
  fOffsety[17] = -0.00029;
  fSigmay[17] = 0.0069319;
  fOffsetx[18] = 0.00064;
  fSigmax[18] = 0.0269351;
  fOffsety[18] = -0.00039;
  fSigmay[18] = 0.0081094;
  fOffsetx[19] = 0.00069;
  fSigmax[19] = 0.0295223;
  fOffsety[19] = -0.00054;
  fSigmay[19] = 0.0075383;

  // cuts on matching - mostly with primary angle
  fOffx[0] = -0.004793;
  fSigx[0] = 0.003;
  fOffy[0] = 0.00747;
  fSigy[0] = 0.0026;

  fOffx[1] = -0.00427;
  fSigx[1] = 0.002;
  fOffy[1] = 0.00186;
  fSigy[1] = 0.0016;

  fOffx[2] = -0.00211;
  fSigx[2] = 0.0018;
  fOffy[2] = 0.00057;
  fSigy[2] = 0.00167;

  fOffx[3] = -0.000357;
  fSigx[3] = 0.0015;
  fOffy[3] = -0.004726;
  fSigy[3] = 0.00165;

  /*
  fOffx[4] = 0.000009;    fSigx[4] =  0.0004344;
  fOffy[4] = 0.000002;    fSigy[4] =  0.0000667;
  fOffx[5] = -0.000006;    fSigx[5] =  0.0002478;
  fOffy[5] = 0.000001;    fSigy[5] =  0.0000736;
  fOffx[6] = -0.000014;    fSigx[6] =  0.0002667;
  fOffy[6] = -0.000001;    fSigy[6] =  0.0000568;
  fOffx[7] = -0.000001;    fSigx[7] =  0.0005851;
  fOffy[7] = 0.000000;    fSigy[7] =  0.0000952;
  fOffx[8] = 0.000009;    fSigx[8] =  0.0004510;
  fOffy[8] = -0.000018;    fSigy[8] =  0.0000697;
  fOffx[9] = 0.000043;    fSigx[9] =  0.0003674;
  fOffy[9] = 0.000006;    fSigy[9] =  0.0000939;
  fOffx[10] = -0.000001;    fSigx[10] =  0.0002548;
  fOffy[10] = -0.000005;    fSigy[10] =  0.0000593;
  fOffx[11] = 0.000022;    fSigx[11] =  0.0005888;
  fOffy[11] = -0.000006;    fSigy[11] =  0.0001153;
  fOffx[12] = 0.00009;      fSigx[12] =  0.001336;
  fOffy[12] = -0.000010;    fSigy[12] =  0.0000981;
  fOffx[13] = 0.000001;    fSigx[13] =  0.0002673;
  fOffy[13] = -0.000006;    fSigy[13] =  0.0000792;
  fOffx[14] = -0.000039;    fSigx[14] =  0.0002979;
  fOffy[14] = -0.000013;    fSigy[14] =  0.0000704;
  fOffx[15] = 0.000074;    fSigx[15] =  0.0010720;
  fOffy[15] = 0.000015;    fSigy[15] =  0.0001014;
  fOffx[16] = 0.000040;    fSigx[16] =  0.0006272;
  fOffy[16] = 0.000016;    fSigy[16] =  0.0001224;
  fOffx[17] = -0.000018;    fSigx[17] =  0.0002965;
  fOffy[17] = 0.000022;    fSigy[17] =  0.0000720;
  */

  fOffx[4] = -0.000044;
  fSigx[4] = 0.0005468;
  fOffy[4] = 0.000003;
  fSigy[4] = 0.0000847;
  fOffx[5] = -0.000011;
  fSigx[5] = 0.0004095;
  fOffy[5] = 0.000000;
  fSigy[5] = 0.0000589;
  fOffx[6] = -0.000012;
  fSigx[6] = 0.0003923;
  fOffy[6] = -0.000002;
  fSigy[6] = 0.0000784;
  fOffx[7] = -0.000019;
  fSigx[7] = 0.0004249;
  fOffy[7] = -0.000002;
  fSigy[7] = 0.0000946;
  fOffx[8] = -0.000027;
  fSigx[8] = 0.0005536;
  fOffy[8] = -0.000013;
  fSigy[8] = 0.0000878;
  fOffx[9] = 0.000028;
  fSigx[9] = 0.0004148;
  fOffy[9] = -0.000000;
  fSigy[9] = 0.0000615;
  fOffx[10] = 0.000029;
  fSigx[10] = 0.0004008;
  fOffy[10] = 0.000003;
  fSigy[10] = 0.0000810;
  fOffx[11] = 0.000026;
  fSigx[11] = 0.0004671;
  fOffy[11] = 0.000005;
  fSigy[11] = 0.0001056;
  fOffx[12] = -0.000031;
  fSigx[12] = 0.0007023;
  fOffy[12] = -0.000005;
  fSigy[12] = 0.0001167;
  fOffx[13] = -0.000033;
  fSigx[13] = 0.0004444;
  fOffy[13] = -0.000003;
  fSigy[13] = 0.0000626;
  fOffx[14] = -0.000077;
  fSigx[14] = 0.0004544;
  fOffy[14] = -0.000004;
  fSigy[14] = 0.0000902;
  fOffx[15] = 0.000011;
  fSigx[15] = 0.0006968;
  fOffy[15] = 0.000006;
  fSigy[15] = 0.0001169;
  fOffx[16] = 0.000009;
  fSigx[16] = 0.0004506;
  fOffy[16] = 0.000005;
  fSigy[16] = 0.0000660;
  fOffx[17] = -0.000025;
  fSigx[17] = 0.0004606;
  fOffy[17] = 0.000007;
  fSigy[17] = 0.0000913;

  // extrapolated tracks cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = -0.03;
}


//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_SaleveOct2017_field_xela75(int /*run_id*/) {

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_SaleveOct2017_field_xela75: setting offsets and sigmas";
    INFO(info); 
    
  // dev cuts

	fOffsetx[0] = 0.00018;    fSigmax[0] =  0.0079099;
	fOffsety[0] = 0.00016;    fSigmay[0] =  0.0054965;
	fOffsetx[1] = 0.00030;    fSigmax[1] =  0.0079130;
	fOffsety[1] = 0.00048;    fSigmay[1] =  0.0055444;
	fOffsetx[2] = 0.00002;    fSigmax[2] =  0.0088651;
	fOffsety[2] = 0.00038;    fSigmay[2] =  0.0058292;
	fOffsetx[3] = 0.00015;    fSigmax[3] =  0.0087304;
	fOffsety[3] = 0.00031;    fSigmay[3] =  0.0058111;
	fOffsetx[4] = 0.00012;    fSigmax[4] =  0.0097437;
	fOffsety[4] = -0.00029;    fSigmay[4] =  0.0065144;
	fOffsetx[5] = 0.00021;    fSigmax[5] =  0.0094949;
	fOffsety[5] = 0.00016;    fSigmay[5] =  0.0065275;


	//based on 33453, 33480, 33481
	//eff=1
	
	fOffsetx[6] = 0.00006;    fSigmax[6] =  0.0090320;
	fOffsety[6] = 0.00039;    fSigmay[6] =  0.0059306;
	fOffsetx[7] = 0.00103;    fSigmax[7] =  0.0146935;
	fOffsety[7] = 0.00027;    fSigmay[7] =  0.0069087;
	fOffsetx[8] = 0.00092;    fSigmax[8] =  0.0140243;
	fOffsety[8] = 0.00019;    fSigmay[8] =  0.0062927;
	fOffsetx[9] = 0.00016;    fSigmax[9] =  0.0082113;
	fOffsety[9] = 0.00015;    fSigmay[9] =  0.0055828;
	fOffsetx[10] = 0.00017;    fSigmax[10] =  0.0088758;
	fOffsety[10] = 0.00029;    fSigmay[10] =  0.0059007;
	fOffsetx[11] = 0.00144;    fSigmax[11] =  0.0150790;
	fOffsety[11] = 0.00044;    fSigmay[11] =  0.0067341;
	fOffsetx[12] = 0.00124;    fSigmax[12] =  0.0135563;
	fOffsety[12] = 0.00081;    fSigmay[12] =  0.0061767;
	fOffsetx[13] = 0.00019;    fSigmax[13] =  0.0082137;
	fOffsety[13] = 0.00049;    fSigmay[13] =  0.0056044;
	fOffsetx[14] = 0.00012;    fSigmax[14] =  0.0099649;
	fOffsety[14] = -0.00029;    fSigmay[14] =  0.0065121;
	fOffsetx[15] = 0.00124;    fSigmax[15] =  0.0158974;
	fOffsety[15] = -0.00018;    fSigmay[15] =  0.0076435;
	fOffsetx[16] = 0.00005;    fSigmax[16] =  0.0164152;
	fOffsety[16] = 0.00047;    fSigmay[16] =  0.0073853;
	fOffsetx[17] = 0.00018;    fSigmax[17] =  0.0096090;
	fOffsety[17] = 0.00016;    fSigmay[17] =  0.0065640;
	fOffsetx[18] = 0.00003;    fSigmax[18] =  0.0155422;
	fOffsety[18] = 0.00050;    fSigmay[18] =  0.0075950;
	fOffsetx[19] = 0.00020;    fSigmax[19] =  0.0161923;
	fOffsety[19] = -0.00007;    fSigmay[19] =  0.0073794;


  // cuts on matching - mostly with primary angle
  fOffx[0] = -0.004793;
  fSigx[0] = 0.003;
  fOffy[0] = 0.00747;
  fSigy[0] = 0.0026;

  fOffx[1] = -0.00427;
  fSigx[1] = 0.002;
  fOffy[1] = 0.00186;
  fSigy[1] = 0.0016;

  fOffx[2] = -0.00211;
  fSigx[2] = 0.0018;
  fOffy[2] = 0.00057;
  fSigy[2] = 0.00167;

  fOffx[3] = -0.000357;
  fSigx[3] = 0.0015;
  fOffy[3] = -0.004726;
  fSigy[3] = 0.00165;
  
  //based on 33453, 33480, 33481
	//eff=1

	fOffx[4] = -0.000013;    fSigx[4] =  0.0005287;
	fOffy[4] = 0.000000;    fSigy[4] =  0.0000875;
	fOffx[5] = -0.000011;    fSigx[5] =  0.0003322;
	fOffy[5] = 0.000000;    fSigy[5] =  0.0000598;
	fOffx[6] = -0.000020;    fSigx[6] =  0.0003300;
	fOffy[6] = 0.000001;    fSigy[6] =  0.0000814;
	fOffx[7] = -0.000032;    fSigx[7] =  0.0003673;
	fOffy[7] = -0.000000;    fSigy[7] =  0.0000983;
	fOffx[8] = 0.000008;    fSigx[8] =  0.0005289;
	fOffy[8] = -0.000014;    fSigy[8] =  0.0000885;
	fOffx[9] = 0.000036;    fSigx[9] =  0.0003335;
	fOffy[9] = 0.000000;    fSigy[9] =  0.0000600;
	fOffx[10] = 0.000028;    fSigx[10] =  0.0003353;
	fOffy[10] = 0.000005;    fSigy[10] =  0.0000823;
	fOffx[11] = 0.000016;    fSigx[11] =  0.0003737;
	fOffy[11] = 0.000006;    fSigy[11] =  0.0000997;
	fOffx[12] = 0.000004;    fSigx[12] =  0.0006564;
	fOffy[12] = -0.000007;    fSigy[12] =  0.0001095;
	fOffx[13] = -0.000043;    fSigx[13] =  0.0003716;
	fOffy[13] = -0.000004;    fSigy[13] =  0.0000652;
	fOffx[14] = -0.000103;    fSigx[14] =  0.0003699;
	fOffy[14] = -0.000003;    fSigy[14] =  0.0000919;
	fOffx[15] = 0.000031;    fSigx[15] =  0.0006593;
	fOffy[15] = 0.000003;    fSigy[15] =  0.0001101;
	fOffx[16] = 0.000006;    fSigx[16] =  0.0003837;
	fOffy[16] = 0.000007;    fSigy[16] =  0.0000667;
	fOffx[17] = -0.000042;    fSigx[17] =  0.0003754;
	fOffy[17] = 0.000011;    fSigy[17] =  0.0000934;

//DATA BELOW IS NOT USE 

  // extrapolated tracks cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = -0.03;
}


//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_SaleveOct2017_field_xela40(int /*run_id*/) {

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_SaleveOct2017_field_xela40: setting offsets and sigmas";
    INFO(info); 


  // dev cuts

	fOffsetx[0] = -0.00256;    fSigmax[0] =  0.0064242;
	fOffsety[0] = 0.00063;    fSigmay[0] =  0.0055062;
	fOffsetx[1] = -0.00246;    fSigmax[1] =  0.0065947;
	fOffsety[1] = 0.00067;    fSigmay[1] =  0.0055364;
	fOffsetx[2] = 0.00014;    fSigmax[2] =  0.0071549;
	fOffsety[2] = 0.00118;    fSigmay[2] =  0.0060045;
	fOffsetx[3] = 0.00001;    fSigmax[3] =  0.0071295;
	fOffsety[3] = 0.00156;    fSigmay[3] =  0.0059215;
	fOffsetx[4] = 0.00090;    fSigmax[4] =  0.0083239;
	fOffsety[4] = -0.00051;    fSigmay[4] =  0.0066755;
	fOffsetx[5] = 0.00098;    fSigmax[5] =  0.0078897;
	fOffsety[5] = 0.00182;    fSigmay[5] =  0.0069525;




	//based on 35206, 35229, 35264
	//eff=1
	
	fOffsetx[6] = 0.00015;    fSigmax[6] =  0.0072665;
	fOffsety[6] = 0.00119;    fSigmay[6] =  0.0060421;
	fOffsetx[7] = -0.00080;    fSigmax[7] =  0.0097562;
	fOffsety[7] = 0.00112;    fSigmay[7] =  0.0072053;
	fOffsetx[8] = -0.00282;    fSigmax[8] =  0.0086755;
	fOffsety[8] = 0.00105;    fSigmay[8] =  0.0065001;
	fOffsetx[9] = -0.00267;    fSigmax[9] =  0.0066445;
	fOffsety[9] = 0.00060;    fSigmay[9] =  0.0056302;
	fOffsetx[10] = 0.00002;    fSigmax[10] =  0.0071713;
	fOffsety[10] = 0.00155;    fSigmay[10] =  0.0059412;
	fOffsetx[11] = -0.00284;    fSigmax[11] =  0.0106781;
	fOffsety[11] = 0.00234;    fSigmay[11] =  0.0072571;
	fOffsetx[12] = -0.00272;    fSigmax[12] =  0.0084882;
	fOffsety[12] = 0.00189;    fSigmay[12] =  0.0063276;
	fOffsetx[13] = -0.00268;    fSigmax[13] =  0.0067835;
	fOffsety[13] = 0.00063;    fSigmay[13] =  0.0056728;
	fOffsetx[14] = 0.00098;    fSigmax[14] =  0.0082214;
	fOffsety[14] = -0.00052;    fSigmay[14] =  0.0066968;
	fOffsetx[15] = 0.00120;    fSigmax[15] =  0.0114715;
	fOffsety[15] = -0.00007;    fSigmay[15] =  0.0078398;
	fOffsetx[16] = -0.00290;    fSigmax[16] =  0.0099631;
	fOffsety[16] = 0.00079;    fSigmay[16] =  0.0074409;
	fOffsetx[17] = 0.00101;    fSigmax[17] =  0.0081553;
	fOffsety[17] = 0.00186;    fSigmay[17] =  0.0069634;
	fOffsetx[18] = -0.00142;    fSigmax[18] =  0.0120222;
	fOffsety[18] = 0.00327;    fSigmay[18] =  0.0069275;
	fOffsetx[19] = -0.00308;    fSigmax[19] =  0.0099386;
	fOffsety[19] = 0.00129;    fSigmay[19] =  0.0075116;



  // cuts on matching - mostly with primary angle
  fOffx[0] = -0.004793;
  fSigx[0] = 0.003;
  fOffy[0] = 0.00747;
  fSigy[0] = 0.0026;

  fOffx[1] = -0.00427;
  fSigx[1] = 0.002;
  fOffy[1] = 0.00186;
  fSigy[1] = 0.0016;

  fOffx[2] = -0.00211;
  fSigx[2] = 0.0018;
  fOffy[2] = 0.00057;
  fSigy[2] = 0.00167;

  fOffx[3] = -0.000357;
  fSigx[3] = 0.0015;
  fOffy[3] = -0.004726;
  fSigy[3] = 0.00165;
  
	//based on 35206, 35229, 35264
	//eff=1

	fOffx[4] = 0.000091;    fSigx[4] =  0.0005461;
	fOffy[4] = -0.000004;    fSigy[4] =  0.0001029;
	fOffx[5] = -0.000022;    fSigx[5] =  0.0003262;
	fOffy[5] = 0.000002;    fSigy[5] =  0.0000587;
	fOffx[6] = -0.000091;    fSigx[6] =  0.0003122;
	fOffy[6] = 0.000007;    fSigy[6] =  0.0000807;
	fOffx[7] = -0.000120;    fSigx[7] =  0.0003541;
	fOffy[7] = 0.000006;    fSigy[7] =  0.0000962;
	fOffx[8] = 0.000100;    fSigx[8] =  0.0005540;
	fOffy[8] = -0.000026;    fSigy[8] =  0.0001046;
	fOffx[9] = 0.000014;    fSigx[9] =  0.0003346;
	fOffy[9] = 0.000006;    fSigy[9] =  0.0000589;
	fOffx[10] = -0.000048;    fSigx[10] =  0.0003168;
	fOffy[10] = 0.000014;    fSigy[10] =  0.0000818;
	fOffx[11] = -0.000069;    fSigx[11] =  0.0003593;
	fOffy[11] = 0.000011;    fSigy[11] =  0.0000975;
	fOffx[12] = 0.000164;    fSigx[12] =  0.0006673;
	fOffy[12] = -0.000007;    fSigy[12] =  0.0001162;
	fOffx[13] = -0.000009;    fSigx[13] =  0.0003620;
	fOffy[13] = -0.000001;    fSigy[13] =  0.0000623;
	fOffx[14] = -0.000132;    fSigx[14] =  0.0003543;
	fOffy[14] = 0.000004;    fSigy[14] =  0.0000949;
	fOffx[15] = 0.000179;    fSigx[15] =  0.0007044;
	fOffy[15] = -0.000010;    fSigy[15] =  0.0001180;
	fOffx[16] = 0.000350;    fSigx[16] =  0.0003992;
	fOffy[16] = 0.000011;    fSigy[16] =  0.0000695;
	fOffx[17] = -0.000079;    fSigx[17] =  0.0003564;
	fOffy[17] = 0.000018;    fSigy[17] =  0.0000968;


//DATA BELOW IS NOT USE 

  // extrapolated tracks cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = -0.03;
}


//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_SaleveNov2018_field(int /*run_id*/) {
 
    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_SaleveNov2018_field: setting offsets and sigmas";
    INFO(info); 


  // dev cuts

  
  // based on runs 38943, 38945
	fOffsetx[0] = -0.00060;    fSigmax[0] =  0.0116804;
	fOffsety[0] = 0.00038;    fSigmay[0] =  0.0051837;
	fOffsetx[1] = -0.00072;    fSigmax[1] =  0.0117422;
	fOffsety[1] = 0.00003;    fSigmay[1] =  0.0052417;
	fOffsetx[2] = 0.00074;    fSigmax[2] =  0.0128568;
	fOffsety[2] = -0.00045;    fSigmay[2] =  0.0055075;
	fOffsetx[3] = 0.00058;    fSigmax[3] =  0.0128678;
	fOffsety[3] = 0.00015;    fSigmay[3] =  0.0055183;
	fOffsetx[4] = 0.00074;    fSigmax[4] =  0.0147641;
	fOffsety[4] = -0.00062;    fSigmay[4] =  0.0062147;
	fOffsetx[5] = 0.00116;    fSigmax[5] =  0.0152853;
	fOffsety[5] = -0.00012;    fSigmay[5] =  0.0063969;



  // based on runs 38943, 38945
	fOffsetx[6] = 0.00075;    fSigmax[6] =  0.0132473;
	fOffsety[6] = -0.00045;    fSigmay[6] =  0.0056386;
	fOffsetx[7] = 0.00139;    fSigmax[7] =  0.0253810;
	fOffsety[7] = -0.00032;    fSigmay[7] =  0.0066304;
	fOffsetx[8] = 0.00040;    fSigmax[8] =  0.0238988;
	fOffsety[8] = 0.00024;    fSigmay[8] =  0.0061143;
	fOffsetx[9] = -0.00056;    fSigmax[9] =  0.0123725;
	fOffsety[9] = 0.00041;    fSigmay[9] =  0.0053500;
	fOffsetx[10] = 0.00059;    fSigmax[10] =  0.0132757;
	fOffsety[10] = 0.00015;    fSigmay[10] =  0.0056467;
	fOffsetx[11] = 0.00060;    fSigmax[11] =  0.0248210;
	fOffsety[11] = -0.00014;    fSigmay[11] =  0.0066587;
	fOffsetx[12] = -0.00012;    fSigmax[12] =  0.0233672;
	fOffsety[12] = 0.00004;    fSigmay[12] =  0.0060217;
	fOffsetx[13] = -0.00070;    fSigmax[13] =  0.0124342;
	fOffsety[13] = 0.00003;    fSigmay[13] =  0.0054164;
	fOffsetx[14] = 0.00070;    fSigmax[14] =  0.0150905;
	fOffsety[14] = -0.00062;    fSigmay[14] =  0.0063751;
	fOffsetx[15] = 0.00099;    fSigmax[15] =  0.0276534;
	fOffsety[15] = -0.00042;    fSigmay[15] =  0.0075893;
	fOffsetx[16] = 0.00062;    fSigmax[16] =  0.0298607;
	fOffsety[16] = 0.00051;    fSigmay[16] =  0.0072844;
	fOffsetx[17] = 0.00107;    fSigmax[17] =  0.0153759;
	fOffsety[17] = -0.00007;    fSigmay[17] =  0.0064243;
	fOffsetx[18] = 0.00275;    fSigmax[18] =  0.0264285;
	fOffsety[18] = -0.00011;    fSigmay[18] =  0.0074057;
	fOffsetx[19] = 0.00063;    fSigmax[19] =  0.0297405;
	fOffsety[19] = -0.00040;    fSigmay[19] =  0.0071902;


  // cuts on matching - mostly with primary angle
  fOffx[0] = -0.004793;
  fSigx[0] = 0.003;
  fOffy[0] = 0.00747;
  fSigy[0] = 0.0026;

  fOffx[1] = -0.00427;
  fSigx[1] = 0.002;
  fOffy[1] = 0.00186;
  fSigy[1] = 0.0016;

  fOffx[2] = -0.00211;
  fSigx[2] = 0.0018;
  fOffy[2] = 0.00057;
  fSigy[2] = 0.00167;

  fOffx[3] = -0.000357;
  fSigx[3] = 0.0015;
  fOffy[3] = -0.004726;
  fSigy[3] = 0.00165;

	fOffx[4] = 0.000078;    fSigx[4] =  0.0005527;
	fOffy[4] = -0.000006;    fSigy[4] =  0.0000876;
	fOffx[5] = 0.000007;    fSigx[5] =  0.0004131;
	fOffy[5] = -0.000005;    fSigy[5] =  0.0000630;
	fOffx[6] = -0.000027;    fSigx[6] =  0.0004212;
	fOffy[6] = -0.000002;    fSigy[6] =  0.0000830;
	fOffx[7] = -0.000042;    fSigx[7] =  0.0004568;
	fOffy[7] = -0.000002;    fSigy[7] =  0.0001007;
	fOffx[8] = 0.000069;    fSigx[8] =  0.0005559;
	fOffy[8] = 0.000004;    fSigy[8] =  0.0000885;
	fOffx[9] = -0.000002;    fSigx[9] =  0.0004146;
	fOffy[9] = 0.000005;    fSigy[9] =  0.0000634;
	fOffx[10] = -0.000034;    fSigx[10] =  0.0004248;
	fOffy[10] = 0.000006;    fSigy[10] =  0.0000839;
	fOffx[11] = -0.000050;    fSigx[11] =  0.0004621;
	fOffy[11] = 0.000009;    fSigy[11] =  0.0001024;
	fOffx[12] = 0.000056;    fSigx[12] =  0.0006780;
	fOffy[12] = -0.000011;    fSigy[12] =  0.0001095;
	fOffx[13] = -0.000012;    fSigx[13] =  0.0004572;
	fOffy[13] = -0.000008;    fSigy[13] =  0.0000685;
	fOffx[14] = -0.000064;    fSigx[14] =  0.0004781;
	fOffy[14] = -0.000005;    fSigy[14] =  0.0000935;
	fOffx[15] = 0.000087;    fSigx[15] =  0.0006988;
	fOffy[15] = 0.000010;    fSigy[15] =  0.0001117;
	fOffx[16] = -0.000017;    fSigx[16] =  0.0004576;
	fOffy[16] = 0.000008;    fSigy[16] =  0.0000686;
	fOffx[17] = -0.000073;    fSigx[17] =  0.0004817;
	fOffy[17] = 0.000009;    fSigy[17] =  0.0000950;

  // extrapolated tracks cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = -0.03;
}







///////////////////////////////////////////////////////////////////////
///////////////////// no field offsets and sigmas /////////////////////
///////////////////////////////////////////////////////////////////////
//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_JuraDec2016_nofield(int /*run_id*/) {

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_JuraDec2016_nofield: setting offsets and sigmas";
    INFO(info); 


  // dev cuts

  fOffsetx[0] = 0;
  fSigmax[0] = 0.00582;
  fOffsety[0] = 0;
  fSigmay[0] = 0.0055;
  fOffsetx[1] = 0;
  fSigmax[1] = 0.0048;
  fOffsety[1] = 0;
  fSigmay[1] = 0.00521;
  fOffsetx[2] = 0;
  fSigmax[2] = 0.004487;
  fOffsety[2] = 0;
  fSigmay[2] = 0.004463;
  fOffsetx[3] = 0;
  fSigmax[3] = 0.004746;
  fOffsety[3] = 0;
  fSigmay[3] = 0.00621;
  fOffsetx[4] = 0;
  fSigmax[4] = 0.008;
  fOffsety[4] = 0;
  fSigmay[4] = 0.008;
  fOffsetx[5] = 0;
  fSigmax[5] = 0.0078;
  fOffsety[5] = 0;
  fSigmay[5] = 0.0082;

  for (int i = 6; i < 20; i++) {
    fOffsetx[i] = 0.;
    fSigmax[i] = 0.008;
    fOffsety[i] = 0;
    fSigmay[i] = 0.008;
  }

  // cuts on matching - mostly with primary angle
  fOffx[0] = 0.0;
  fSigx[0] = 0.0015;
  fOffy[0] = 0.0;
  fSigy[0] = 0.0015;

  fOffx[1] = 0.0;
  fSigx[1] = 0.0015;
  fOffy[1] = 0.0;
  fSigy[1] = 0.0015;

  fOffx[2] = 0.0;
  fSigx[2] = 0.0015;
  fOffy[2] = 0.0;
  fSigy[2] = 0.0015;

  fOffx[3] = 0.0;
  fSigx[3] = 0.0015;
  fOffy[3] = 0.0;
  fSigy[3] = 0.0015;

  for (int i = 4; i < 18; i++) {
    fOffx[i] = 0.;
    fSigx[i] = 0.00012;
    fOffy[i] = 0;
    fSigy[i] = 0.00012;
  }

  // extrapolated track cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = 0.046;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_SaleveDec2016_nofield(int /*run_id*/) {

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_SaleveDec2016_nofield: setting offsets and sigmas";
    INFO(info); 


  // dev cuts

  fOffsetx[0] = 0.0;
  fSigmax[0] = 0.008;
  fOffsety[0] = 0.0;
  fSigmay[0] = 0.008;
  fOffsetx[1] = 0.0;
  fSigmax[1] = 0.008;
  fOffsety[1] = 0.0;
  fSigmay[1] = 0.008;
  fOffsetx[2] = 0.0;
  fSigmax[2] = 0.008;
  fOffsety[2] = 0.0;
  fSigmay[2] = 0.008;
  fOffsetx[3] = 0.0;
  fSigmax[3] = 0.008;
  fOffsety[3] = 0.0;
  fSigmay[3] = 0.008;

  fOffsetx[4] = 0.0;
  fSigmax[4] = 0.008;
  fOffsety[4] = 0.0;
  fSigmay[4] = 0.008;
  fOffsetx[5] = 0.0;
  fSigmax[5] = 0.01;
  fOffsety[5] = 0.0;
  fSigmay[5] = 0.012;

  for (int i = 6; i < 20; i++) {
    fOffsetx[i] = 0.;
    fSigmax[i] = 0.008;
    fOffsety[i] = 0;
    fSigmay[i] = 0.008;
  }

  // cuts on matching - mostly with primary angle
  fOffx[0] = -0.004793;
  fSigx[0] = 0.003;
  fOffy[0] = 0.00747;
  fSigy[0] = 0.0026;

  fOffx[1] = -0.00427;
  fSigx[1] = 0.002;
  fOffy[1] = 0.00186;
  fSigy[1] = 0.0016;

  fOffx[2] = -0.00211;
  fSigx[2] = 0.0018;
  fOffy[2] = 0.00057;
  fSigy[2] = 0.00167;

  fOffx[3] = -0.000357;
  fSigx[3] = 0.0015;
  fOffy[3] = -0.004726;
  fSigy[3] = 0.00165;

  for (int i = 4; i < 18; i++) {
    fOffx[i] = 0.;
    fSigx[i] = 0.00015;
    fOffy[i] = 0;
    fSigy[i] = 0.00015;
  }

  // extrapolated tracks cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = -0.046;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_JuraOct2017_nofield(int /*run_id*/) {

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_JuraOct2017_nofield: setting offsets and sigmas";
    INFO(info); 


  // dev cuts


  fOffsetx[0] = -0.00007;
  fSigmax[0] = 0.0054098;
  fOffsety[0] = -0.00005;
  fSigmay[0] = 0.0053360;
  fOffsetx[1] = 0.00015;
  fSigmax[1] = 0.0054235;
  fOffsety[1] = -0.00000;
  fSigmay[1] = 0.0053120;
  fOffsetx[2] = -0.00033;
  fSigmax[2] = 0.0056793;
  fOffsety[2] = -0.00001;
  fSigmay[2] = 0.0053938;
  fOffsetx[3] = -0.00020;
  fSigmax[3] = 0.0060010;
  fOffsety[3] = -0.00007;
  fSigmay[3] = 0.0058500;
  fOffsetx[4] = -0.00002;
  fSigmax[4] = 0.0062564;
  fOffsety[4] = 0.00002;
  fSigmay[4] = 0.0058858;
  fOffsetx[5] = -0.00001;
  fSigmax[5] = 0.0068040;
  fOffsety[5] = -0.00003;
  fSigmay[5] = 0.0065662;

  fOffsetx[6] = -0.00044;
  fSigmax[6] = 0.0057246;
  fOffsety[6] = -0.00002;
  fSigmay[6] = 0.0054794;
  fOffsetx[7] = 0.00865;
  fSigmax[7] = 0.0293743;
  fOffsety[7] = 0.00718;
  fSigmay[7] = 0.0364568;
  fOffsetx[8] = 0.00073;
  fSigmax[8] = 0.0073857;
  fOffsety[8] = -0.00047;
  fSigmay[8] = 0.0072663;
  fOffsetx[9] = -0.00015;
  fSigmax[9] = 0.0062987;
  fOffsety[9] = -0.00011;
  fSigmay[9] = 0.0061858;
  fOffsetx[10] = -0.00025;
  fSigmax[10] = 0.0060629;
  fOffsety[10] = -0.00008;
  fSigmay[10] = 0.0059030;
  fOffsetx[11] = 0.00103;
  fSigmax[11] = 0.0144847;
  fOffsety[11] = 0.00078;
  fSigmay[11] = 0.0080267;
  fOffsetx[12] = 0.00022;
  fSigmax[12] = 0.0065777;
  fOffsety[12] = 0.00009;
  fSigmay[12] = 0.0065996;
  fOffsetx[13] = 0.00017;
  fSigmax[13] = 0.0069626;
  fOffsety[13] = -0.00005;
  fSigmay[13] = 0.0066611;
  fOffsetx[14] = 0.00019;
  fSigmax[14] = 0.0128119;
  fOffsety[14] = 0.00034;
  fSigmay[14] = 0.0101312;
  fOffsetx[15] = 0.00206;
  fSigmax[15] = 0.0324652;
  fOffsety[15] = -0.00000;
  fSigmay[15] = 0.0233958;
  fOffsetx[16] = -0.00029;
  fSigmax[16] = 0.0078949;
  fOffsety[16] = 0.00000;
  fSigmay[16] = 0.0077078;
  fOffsetx[17] = 0.00031;
  fSigmax[17] = 0.0132034;
  fOffsety[17] = -0.00034;
  fSigmay[17] = 0.0104014;
  fOffsetx[18] = -0.00055;
  fSigmax[18] = 0.0164523;
  fOffsety[18] = -0.00062;
  fSigmay[18] = 0.0105435;
  fOffsetx[19] = 0.00023;
  fSigmax[19] = 0.0079788;
  fOffsety[19] = 0.00054;
  fSigmay[19] = 0.0077291;

  /*  setup for old geometry
  fOffsetx[0] = -0.00016;    fSigmax[0] =  0.0058413;
  fOffsety[0] = -0.00001;    fSigmay[0] =  0.0055257;
  fOffsetx[1] = -0.00133;    fSigmax[1] =  0.0055104;
  fOffsety[1] = 0.00098;    fSigmay[1] =  0.0053778;
  fOffsetx[2] = -0.00291;    fSigmax[2] =  0.0062358;
  fOffsety[2] = 0.00318;    fSigmay[2] =  0.0057740;
  fOffsetx[3] = -0.00293;    fSigmax[3] =  0.0061868;
  fOffsety[3] = -0.00021;    fSigmay[3] =  0.0059556;
  fOffsetx[4] = 0.00076;    fSigmax[4] =  0.0063552;
  fOffsety[4] = -0.00133;    fSigmay[4] =  0.0074541;
  fOffsetx[5] = 0.00167;    fSigmax[5] =  0.0069777;
  fOffsety[5] = -0.00074;    fSigmay[5] =  0.0075350;


  fOffsetx[6] = -0.00288;    fSigmax[6] =  0.0062755;
  fOffsety[6] = 0.00311;    fSigmay[6] =  0.0057560;
  fOffsetx[7] = -0.00565;    fSigmax[7] =  0.0084062;
  fOffsety[7] = 0.00519;    fSigmay[7] =  0.0068853;
  fOffsetx[8] = -0.00273;    fSigmax[8] =  0.0073814;
  fOffsety[8] = 0.00252;    fSigmay[8] =  0.0061624;
  fOffsetx[9] = -0.00016;    fSigmax[9] =  0.0059260;
  fOffsety[9] = 0.00002;    fSigmay[9] =  0.0055786;
  fOffsetx[10] = -0.00295;    fSigmax[10] =  0.0062360;
  fOffsety[10] = -0.00024;    fSigmay[10] =  0.0059470;
  fOffsetx[11] = -0.00656;    fSigmax[11] =  0.0075749;
  fOffsety[11] = -0.00021;    fSigmay[11] =  0.0071082;
  fOffsetx[12] = -0.00420;    fSigmax[12] =  0.0063554;
  fOffsety[12] = 0.00078;    fSigmay[12] =  0.0059595;
  fOffsetx[13] = -0.00131;    fSigmax[13] =  0.0055540;
  fOffsety[13] = 0.00096;    fSigmay[13] =  0.0053906;
  fOffsetx[14] = 0.00075;    fSigmax[14] =  0.0064390;
  fOffsety[14] = -0.00142;    fSigmay[14] =  0.0074341;
  fOffsetx[15] = 0.00061;    fSigmax[15] =  0.0073912;
  fOffsety[15] = -0.00164;    fSigmay[15] =  0.0087077;
  fOffsetx[16] = -0.00184;    fSigmax[16] =  0.0080808;
  fOffsety[16] = 0.00091;    fSigmay[16] =  0.0077315;
  fOffsetx[17] = 0.00171;    fSigmax[17] =  0.0069553;
  fOffsety[17] = -0.00092;    fSigmay[17] =  0.0074939;
  fOffsetx[18] = 0.00099;    fSigmax[18] =  0.0083628;
  fOffsety[18] = -0.00094;    fSigmay[18] =  0.0093026;
  fOffsetx[19] = -0.00296;    fSigmax[19] =  0.0080478;
  fOffsety[19] = 0.00062;    fSigmay[19] =  0.0078118;
  */

  // for(int i=0;i<20;i++)fOffsety[i]=0;

  // cuts on matching - mostly with primary vertex. 0-3 are not used (4hit tracks)
  fOffx[0] = 0.0;
  fSigx[0] = 0.0015;
  fOffy[0] = 0.0;
  fSigy[0] = 0.0015;

  fOffx[1] = 0.0;
  fSigx[1] = 0.0015;
  fOffy[1] = 0.0;
  fSigy[1] = 0.0015;

  fOffx[2] = 0.0;
  fSigx[2] = 0.0015;
  fOffy[2] = 0.0;
  fSigy[2] = 0.0015;

  fOffx[3] = 0.0;
  fSigx[3] = 0.0015;
  fOffy[3] = 0.0;
  fSigy[3] = 0.0015;

  fOffx[4] = -0.000003;
  fSigx[4] = 0.0000801;
  fOffy[4] = -0.000003;
  fSigy[4] = 0.0000837;
  fOffx[5] = 0.000001;
  fSigx[5] = 0.0000507;
  fOffy[5] = -0.000006;
  fSigy[5] = 0.0000559;
  fOffx[6] = 0.000000;
  fSigx[6] = 0.0000547;
  fOffy[6] = -0.000006;
  fSigy[6] = 0.0000596;
  fOffx[7] = 0.000003;
  fSigx[7] = 0.0000712;
  fOffy[7] = -0.000008;
  fSigy[7] = 0.0000770;
  fOffx[8] = 0.000002;
  fSigx[8] = 0.0000780;
  fOffy[8] = 0.000007;
  fSigy[8] = 0.0000801;
  fOffx[9] = 0.000006;
  fSigx[9] = 0.0000516;
  fOffy[9] = 0.000006;
  fSigy[9] = 0.0000554;
  fOffx[10] = 0.000005;
  fSigx[10] = 0.0000538;
  fOffy[10] = 0.000006;
  fSigy[10] = 0.0000586;
  fOffx[11] = 0.000006;
  fSigx[11] = 0.0000712;
  fOffy[11] = 0.000004;
  fSigy[11] = 0.0000764;
  fOffx[12] = 0.000004;
  fSigx[12] = 0.0001045;
  fOffy[12] = -0.000007;
  fSigy[12] = 0.0001104;
  fOffx[13] = 0.000003;
  fSigx[13] = 0.0000562;
  fOffy[13] = -0.000007;
  fSigy[13] = 0.0000606;
  fOffx[14] = 0.000005;
  fSigx[14] = 0.0000705;
  fOffy[14] = -0.000010;
  fSigy[14] = 0.0000732;
  fOffx[15] = -0.000001;
  fSigx[15] = 0.0001048;
  fOffy[15] = 0.000003;
  fSigy[15] = 0.0001114;
  fOffx[16] = 0.000004;
  fSigx[16] = 0.0000577;
  fOffy[16] = 0.000003;
  fSigy[16] = 0.0000620;
  fOffx[17] = 0.000005;
  fSigx[17] = 0.0000703;
  fOffy[17] = 0.000001;
  fSigy[17] = 0.0000727;

  /*
  fOffx[4] = 0.00001;    fSigx[4] =  0.0001002;
  fOffy[4] = -0.00002;    fSigy[4] =  0.0000910;
  fOffx[5] = -0.00000;    fSigx[5] =  0.0000551;
  fOffy[5] = -0.00001;    fSigy[5] =  0.0000588;
  fOffx[6] = 0.00001;    fSigx[6] =  0.0000639;
  fOffy[6] = -0.00002;    fSigy[6] =  0.0000641;
  fOffx[7] = -0.00001;    fSigx[7] =  0.0000762;
  fOffy[7] = -0.00000;    fSigy[7] =  0.0000825;
  fOffx[8] = 0.00004;    fSigx[8] =  0.0000910;
  fOffy[8] = -0.00000;    fSigy[8] =  0.0000857;
  fOffx[9] = 0.00002;    fSigx[9] =  0.0000560;
  fOffy[9] = 0.00001;    fSigy[9] =  0.0000601;
  fOffx[10] = 0.00003;    fSigx[10] =  0.0000609;
  fOffy[10] = 0.00001;    fSigy[10] =  0.0000655;
  fOffx[11] = 0.00002;    fSigx[11] =  0.0000773;
  fOffy[11] = 0.00000;    fSigy[11] =  0.0000827;
  fOffx[12] = -0.00001;    fSigx[12] =  0.0001155;
  fOffy[12] = 0.00000;    fSigy[12] =  0.0001180;
  fOffx[13] = -0.00001;    fSigx[13] =  0.0000606;
  fOffy[13] = -0.00000;    fSigy[13] =  0.0000632;
  fOffx[14] = -0.00000;    fSigx[14] =  0.0000746;
  fOffy[14] = -0.00001;    fSigy[14] =  0.0000786;
  fOffx[15] = 0.00002;    fSigx[15] =  0.0001138;
  fOffy[15] = -0.00001;    fSigy[15] =  0.0001217;
  fOffx[16] = 0.00001;    fSigx[16] =  0.0000617;
  fOffy[16] = 0.00000;    fSigy[16] =  0.0000644;
  fOffx[17] = 0.00002;    fSigx[17] =  0.0000769;
  fOffy[17] = -0.00000;    fSigy[17] =  0.0000798;
  */

  // extrapolated track cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = 0.046;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_SaleveOct2017_nofield(int /*run_id*/) {

    ostringstream info;
    info << "Na61ArmParameters::SetupOffsetsAndSigmas_SaleveOct2017_nofield: setting offsets and sigmas";
    INFO(info); 

  // dev cuts

  fOffsetx[0] = -0.00007;
  fSigmax[0] = 0.0059687;
  fOffsety[0] = 0.00006;
  fSigmay[0] = 0.0057491;
  fOffsetx[1] = 0.00006;
  fSigmax[1] = 0.0059314;
  fOffsety[1] = -0.00013;
  fSigmay[1] = 0.0058535;
  fOffsetx[2] = -0.00043;
  fSigmax[2] = 0.0064598;
  fOffsety[2] = 0.00055;
  fSigmay[2] = 0.0060900;
  fOffsetx[3] = -0.00055;
  fSigmax[3] = 0.0063018;
  fOffsety[3] = 0.00014;
  fSigmay[3] = 0.0059516;
  fOffsetx[4] = 0.00004;
  fSigmax[4] = 0.0072926;
  fOffsety[4] = -0.00002;
  fSigmay[4] = 0.0068238;
  fOffsetx[5] = 0.00005;
  fSigmax[5] = 0.0072164;
  fOffsety[5] = 0.00006;
  fSigmay[5] = 0.0069546;

  fOffsetx[6] = -0.00055;
  fSigmax[6] = 0.0064695;
  fOffsety[6] = 0.00072;
  fSigmay[6] = 0.0060847;
  fOffsetx[7] = -0.00101;
  fSigmax[7] = 0.0190714;
  fOffsety[7] = -0.00186;
  fSigmay[7] = 0.0050000;
  fOffsetx[8] = 0.00021;
  fSigmax[8] = 0.0069880;
  fOffsety[8] = -0.00012;
  fSigmay[8] = 0.0062714;
  fOffsetx[9] = -0.00010;
  fSigmax[9] = 0.0130810;
  fOffsety[9] = -0.00142;
  fSigmay[9] = 0.0107457;
  fOffsetx[10] = -0.00066;
  fSigmax[10] = 0.0063093;
  fOffsety[10] = 0.00021;
  fSigmay[10] = 0.0059518;
  fOffsetx[11] = 0.00;
  fSigmax[11] = 0.0;  // very bad quality
  fOffsety[11] = 0.00;
  fSigmay[11] = 0.0;
  fOffsetx[12] = 0.00007;
  fSigmax[12] = 0.0064932;
  fOffsety[12] = 0.00012;
  fSigmay[12] = 0.0064566;
  fOffsetx[13] = 0.00009;
  fSigmax[13] = 0.0087780;
  fOffsety[13] = -0.00021;
  fSigmay[13] = 0.0097202;
  fOffsetx[14] = 0.00065;
  fSigmax[14] = 0.0088776;
  fOffsety[14] = 0.00026;
  fSigmay[14] = 0.0087667;
  fOffsetx[15] = 0.00354;
  fSigmax[15] = 0.0079467;
  fOffsety[15] = 0.00124;
  fSigmay[15] = 0.0096218;
  fOffsetx[16] = -0.00018;
  fSigmax[16] = 0.0084637;
  fOffsety[16] = 0.00095;
  fSigmay[16] = 0.0082130;
  fOffsetx[17] = 0.00028;
  fSigmax[17] = 0.0082286;
  fOffsety[17] = 0.00041;
  fSigmay[17] = 0.0076164;
  fOffsetx[18] = 0.00085;
  fSigmax[18] = 0.0113419;
  fOffsety[18] = 0.00342;
  fSigmay[18] = 0.0271124;
  fOffsetx[19] = 0.00004;
  fSigmax[19] = 0.0085144;
  fOffsety[19] = -0.00133;
  fSigmay[19] = 0.0082982;

  /*
  fOffsetx[0] = -0.00591;    fSigmax[0] =  0.0073936;
  fOffsety[0] = -0.00345;    fSigmay[0] =  0.0100033;
  fOffsetx[1] = -0.00417;    fSigmax[1] =  0.0068586;
  fOffsety[1] = -0.00088;    fSigmay[1] =  0.0096713;
  fOffsetx[2] = 0.00588;    fSigmax[2] =  0.0155243;
  fOffsety[2] = 0.00307;    fSigmay[2] =  0.0105382;
  fOffsetx[3] = -0.00100;    fSigmax[3] =  0.0143850;
  fOffsety[3] = 0.00139;    fSigmay[3] =  0.0116478;
  fOffsetx[4] = 0.00536;    fSigmax[4] =  0.0085092;
  fOffsety[4] = 0.00286;    fSigmay[4] =  0.0101264;
  fOffsetx[5] = 0.00518;    fSigmax[5] =  0.0091198;
  fOffsety[5] = 0.00020;    fSigmay[5] =  0.0082309;

  fOffsetx[6] = 0.00552;    fSigmax[6] =  0.0151860;
  fOffsety[6] = 0.00182;    fSigmay[6] =  0.0100417;
  fOffsetx[7] = 0.01375;    fSigmax[7] =  0.0210235;
  fOffsety[7] = -0.00875;    fSigmay[7] =  0.0100389;
  fOffsetx[8] = 0.00287;    fSigmax[8] =  0.0108811;
  fOffsety[8] = -0.00644;    fSigmay[8] =  0.0092562;
  fOffsetx[9] = -0.00603;    fSigmax[9] =  0.0077467;
  fOffsety[9] = -0.00325;    fSigmay[9] =  0.0101643;
  fOffsetx[10] = -0.00075;    fSigmax[10] =  0.0143455;
  fOffsety[10] = 0.00172;    fSigmay[10] =  0.0116169;
  fOffsetx[11] = 0.00276;    fSigmax[11] =  0.0183304;
  fOffsety[11] = -0.00997;    fSigmay[11] =  0.0108560;
  fOffsetx[12] = -0.00251;    fSigmax[12] =  0.0101706;
  fOffsety[12] = -0.00569;    fSigmay[12] =  0.0089842;
  fOffsetx[13] = -0.00421;    fSigmax[13] =  0.0073552;
  fOffsety[13] = -0.00042;    fSigmay[13] =  0.0098407;
  fOffsetx[14] = 0.00536;    fSigmax[14] =  0.0088812;
  fOffsety[14] = 0.00263;    fSigmay[14] =  0.0101725;
  fOffsetx[15] = 0.00269;    fSigmax[15] =  0.0097538;
  fOffsety[15] = 0.00189;    fSigmay[15] =  0.0101131;
  fOffsetx[16] = -0.00625;    fSigmax[16] =  0.0082257;
  fOffsety[16] = -0.00059;    fSigmay[16] =  0.0098331;
  fOffsetx[17] = 0.00520;    fSigmax[17] =  0.0094730;
  fOffsety[17] = -0.00008;    fSigmay[17] =  0.0085651;
  fOffsetx[18] = 0.00359;    fSigmax[18] =  0.0110502;
  fOffsety[18] = -0.00005;    fSigmay[18] =  0.0091707;
  fOffsetx[19] = -0.00413;    fSigmax[19] =  0.0085227;
  fOffsety[19] = -0.00156;    fSigmay[19] =  0.0105453;
  */

  // for(int i=0;i<20;i++)fOffsety[i]=0;

  // cuts on matching - mostly with primary angle 0 - 3 are not used (4hit tracks)
  fOffx[0] = -0.004793;
  fSigx[0] = 0.003;
  fOffy[0] = 0.00747;
  fSigy[0] = 0.0026;

  fOffx[1] = -0.00427;
  fSigx[1] = 0.002;
  fOffy[1] = 0.00186;
  fSigy[1] = 0.0016;

  fOffx[2] = -0.00211;
  fSigx[2] = 0.0018;
  fOffy[2] = 0.00057;
  fSigy[2] = 0.00167;

  fOffx[3] = -0.000357;
  fSigx[3] = 0.0015;
  fOffy[3] = -0.004726;
  fSigy[3] = 0.00165;

  /// values taken for fit using FitMatchingWithVertex.C macro

  fOffx[4] = -0.000001;
  fSigx[4] = 0.0000978;
  fOffy[4] = -0.000004;
  fSigy[4] = 0.0001028;
  fOffx[5] = -0.000006;
  fSigx[5] = 0.0000599;
  fOffy[5] = -0.000005;
  fSigy[5] = 0.0000653;
  fOffx[6] = -0.000006;
  fSigx[6] = 0.0000648;
  fOffy[6] = -0.000005;
  fSigy[6] = 0.0000723;
  fOffx[7] = -0.000011;
  fSigx[7] = 0.0000852;
  fOffy[7] = -0.000007;
  fSigy[7] = 0.0000921;
  fOffx[8] = 0.000007;
  fSigx[8] = 0.0000968;
  fOffy[8] = -0.000003;
  fSigy[8] = 0.0001016;
  fOffx[9] = 0.000007;
  fSigx[9] = 0.0000598;
  fOffy[9] = 0.000005;
  fSigy[9] = 0.0000665;
  fOffx[10] = 0.000007;
  fSigx[10] = 0.0000643;
  fOffy[10] = 0.000007;
  fSigy[10] = 0.0000736;
  fOffx[11] = 0.000007;
  fSigx[11] = 0.0000852;
  fOffy[11] = 0.000009;
  fSigy[11] = 0.0000919;
  fOffx[12] = -0.000007;
  fSigx[12] = 0.0001172;
  fOffy[12] = -0.000012;
  fSigy[12] = 0.0001239;
  fOffx[13] = -0.000010;
  fSigx[13] = 0.0000638;
  fOffy[13] = -0.000007;
  fSigy[13] = 0.0000691;
  fOffx[14] = -0.000019;
  fSigx[14] = 0.0000803;
  fOffy[14] = -0.000012;
  fSigy[14] = 0.0000858;
  fOffx[15] = 0.000004;
  fSigx[15] = 0.0001170;
  fOffy[15] = 0.000015;
  fSigy[15] = 0.0001246;
  fOffx[16] = 0.000004;
  fSigx[16] = 0.0000642;
  fOffy[16] = 0.000009;
  fSigy[16] = 0.0000681;
  fOffx[17] = -0.000004;
  fSigx[17] = 0.0000797;
  fOffy[17] = 0.000016;
  fSigy[17] = 0.0000861;

  /*
  fOffx[4] = -0.00006;    fSigx[4] =  0.0001366;
  fOffy[4] = -0.00001;    fSigy[4] =  0.0001304;
  fOffx[5] = -0.00001;    fSigx[5] =  0.0000602;
  fOffy[5] = -0.00001;    fSigy[5] =  0.0000666;
  fOffx[6] = -0.00003;    fSigx[6] =  0.0000797;
  fOffy[6] = 0.00001;    fSigy[6] =  0.0000767;
  fOffx[7] = 0.00002;    fSigx[7] =  0.0000884;
  fOffy[7] = -0.00002;    fSigy[7] =  0.0001011;
  fOffx[8] = 0.00000;    fSigx[8] =  0.0001427;
  fOffy[8] = 0.00003;    fSigy[8] =  0.0001334;
  fOffx[9] = 0.00000;    fSigx[9] =  0.0000592;
  fOffy[9] = 0.00002;    fSigy[9] =  0.0000674;
  fOffx[10] = 0.00001;    fSigx[10] =  0.0000802;
  fOffy[10] = 0.00004;    fSigy[10] =  0.0000775;
  fOffx[11] = 0.00001;    fSigx[11] =  0.0000860;
  fOffy[11] = 0.00001;    fSigy[11] =  0.0000997;
  fOffx[12] = 0.00008;    fSigx[12] =  0.0001158;
  fOffy[12] = 0.00001;    fSigy[12] =  0.0001527;
  fOffx[13] = 0.00002;    fSigx[13] =  0.0000649;
  fOffy[13] = -0.00001;    fSigy[13] =  0.0000743;
  fOffx[14] = 0.00001;    fSigx[14] =  0.0000791;
  fOffy[14] = -0.00002;    fSigy[14] =  0.0000904;
  fOffx[15] = 0.00004;    fSigx[15] =  0.0001197;
  fOffy[15] = 0.00001;    fSigy[15] =  0.0001653;
  fOffx[16] = 0.00001;    fSigx[16] =  0.0000639;
  fOffy[16] = 0.00001;    fSigy[16] =  0.0000734;
  fOffx[17] = 0.00001;    fSigx[17] =  0.0000803;
  fOffy[17] = 0.00002;    fSigy[17] =  0.0000919;
  */

  // extrapolated tracks cuts
  fOff_dx[0] = 0.0;
  fSig_dx[0] = 0.008;
  fOff_dy[0] = 0.0;
  fSig_dy[0] = 0.008;
  fOff_dax[0] = 0.0;
  fSig_dax[0] = 0.00018;
  fOff_day[0] = 0.0;
  fSig_day[0] = 0.00018;

  fOff_dx[1] = 0.0;
  fSig_dx[1] = 0.008;
  fOff_dy[1] = 0.0;
  fSig_dy[1] = 0.008;
  fOff_dax[1] = 0.0;
  fSig_dax[1] = 0.00018;
  fOff_day[1] = 0.0;
  fSig_day[1] = 0.00108;

  fOff_dx[2] = 0.0;
  fSig_dx[2] = 0.008;
  fOff_dy[2] = 0.0;
  fSig_dy[2] = 0.008;
  fOff_dax[2] = 0.0;
  fSig_dax[2] = 0.00018;
  fOff_day[2] = 0.0;
  fSig_day[2] = 0.00108;

  fOff_dx[3] = 0.0;
  fSig_dx[3] = 0.008;
  fOff_dy[3] = 0.0;
  fSig_dy[3] = 0.008;
  fOff_dax[3] = 0.0;
  fSig_dax[3] = 0.00018;
  fOff_day[3] = 0.0;
  fSig_day[3] = 0.00108;

  fAxCut = -0.03;
}

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////// below archive from July2016 ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorRotation_July2016() {
  double alpha_Vds1_0 = 0.009;
  double alpha_Vds2_0 = 0.0;
  double alpha_Vds3_0 = 0.0;
  double alpha_Vds3_1 = 0.0;
  double alpha_Vds4_0 = -0.0094;
  double alpha_Vds4_1 = -0.00678;

  fRotZ[0] = alpha_Vds1_0;
  fRotZ[1] = alpha_Vds2_0;
  fRotZ[2] = alpha_Vds3_0;
  fRotZ[3] = alpha_Vds3_1;
  fRotZ[4] = alpha_Vds4_0;
  fRotZ[5] = alpha_Vds4_1;
}

//____________________________________________________________________________________________
void Na61ArmParameters::SetupSensorGeometry_July2016(int run_id) {
  // nominal sensor positions base on run 81
  double VolumeX_Vds1_0 = -0.024942;
  double VolumeY_Vds1_0 = -0.002195;

  double VolumeX_Vds2_0 = 0.0;
  double VolumeY_Vds2_0 = 0.0;

  double VolumeX_Vds3_0 = -3.12494;
  double VolumeY_Vds3_0 = 10.9978;

  double VolumeX_Vds3_1 = -3.118417;
  double VolumeY_Vds3_1 = -10.76101;

  double VolumeX_Vds4_0 = -5.92815;
  double VolumeY_Vds4_0 = 11.2315;

  double VolumeX_Vds4_1 = -5.78953;  // based on run 106
  double VolumeY_Vds4_1 = -10.564;

  // fVolumeX_Vds4_1 = -5.92815  + 0.1354 + 0.007432 - 0.004214; // measurement relative to Vds4_0
  // fVolumeY_Vds4_1 = -10.76101 + 0.3166 - 0.003632 - 0.11594;

  if (run_id == 100) {
    VolumeX_Vds1_0 = VolumeX_Vds1_0 - 0.000791717;
    VolumeY_Vds1_0 = VolumeY_Vds1_0 + 0.0021915;
    VolumeX_Vds2_0 = VolumeX_Vds2_0 + 0.000902856;
    VolumeY_Vds2_0 = VolumeY_Vds2_0 - 0.00155545;
    VolumeX_Vds3_0 = VolumeX_Vds3_0 + 0.000569438;
    VolumeY_Vds3_0 = VolumeY_Vds3_0 - 0.00346361;
    VolumeX_Vds4_0 = VolumeX_Vds4_0 - 0.000680578;
    VolumeY_Vds4_0 = VolumeY_Vds4_0 + 0.00282755;
  }

  if (run_id >= 106) {
    VolumeX_Vds1_0 = VolumeX_Vds1_0 - 0.00571882;
    VolumeY_Vds1_0 = VolumeY_Vds1_0 - 0.000613938;
    VolumeX_Vds2_0 = VolumeX_Vds2_0 + 0.00676695;
    VolumeY_Vds2_0 = VolumeY_Vds2_0 + 0.000590428;
    VolumeX_Vds3_0 = VolumeX_Vds3_0 + 0.00362255;
    VolumeY_Vds3_0 = VolumeY_Vds3_0 + 0.000660957;
    VolumeX_Vds4_0 = VolumeX_Vds4_0 - 0.00467069;
    VolumeY_Vds4_0 = VolumeY_Vds4_0 - 0.000637448;
  }

  if (run_id >= 116) {
    VolumeX_Vds1_0 = VolumeX_Vds1_0 + 0.0161254 - 0.00039276;
    VolumeY_Vds1_0 = VolumeY_Vds1_0 + 0.025176 + 0.001408;
    VolumeX_Vds2_0 = VolumeX_Vds2_0 - 0.0107872 + 0.00030468;
    VolumeY_Vds2_0 = VolumeY_Vds2_0 - 0.024738 - 0.002544;
    VolumeX_Vds3_0 = VolumeX_Vds3_0 - 0.0268018 + 0.00056892;
    VolumeY_Vds3_0 = VolumeY_Vds3_0 - 0.026052 + 0.000864;
    VolumeX_Vds4_0 = VolumeX_Vds4_0 + 0.0214636 - 0.00048084;
    VolumeY_Vds4_0 = VolumeY_Vds4_0 + 0.025614 + 0.000272;
    VolumeX_Vds4_1 = VolumeX_Vds4_1 + 0.02092;
    VolumeY_Vds4_1 = VolumeY_Vds4_1 + 0.11508;
  }

  fVolumeX[0] = VolumeX_Vds1_0;
  fVolumeX[1] = VolumeX_Vds2_0;
  fVolumeX[2] = VolumeX_Vds3_0;
  fVolumeX[3] = VolumeX_Vds3_1;
  fVolumeX[4] = VolumeX_Vds4_0;
  fVolumeX[5] = VolumeX_Vds4_1;

  fVolumeY[0] = VolumeY_Vds1_0;
  fVolumeY[1] = VolumeY_Vds2_0;
  fVolumeY[2] = VolumeY_Vds3_0;
  fVolumeY[3] = VolumeY_Vds3_1;
  fVolumeY[4] = VolumeY_Vds4_0;
  fVolumeY[5] = VolumeY_Vds4_1;
}

//_____________________________________________________________
void Na61ArmParameters::SetupOffsetsAndSigmas_July2016(int run_id) {
  // dev cuts

  fOffsetx[0] = 0.;
  fSigmax[0] = 0.006;
  fOffsety[0] = 0.;
  fSigmay[0] = 0.006;

  fOffsetx[1] = 0.;
  fSigmax[1] = 0.006;
  fOffsety[1] = 0.;
  fSigmay[1] = 0.006;

  fOffsetx[2] = 0.;
  fSigmax[2] = 0.006;
  fOffsety[2] = 0.;
  fSigmay[2] = 0.006;

  fOffsetx[3] = 0.;
  fSigmax[3] = 0.006;
  fOffsety[3] = 0.;
  fSigmay[3] = 0.006;

  if (run_id > 110) {  // produced particles
    fOffsetx[0] = 0.;
    fSigmax[0] = 0.015;
    fOffsety[0] = 0.;
    fSigmay[0] = 0.015;

    fOffsetx[1] = 0.045;
    fSigmax[1] = 0.021;
    fOffsety[1] = 0.02122;
    fSigmay[1] = 0.026;

    fOffsetx[2] = 0.;
    fSigmax[2] = 0.015;
    fOffsety[2] = 0.;
    fSigmay[2] = 0.015;

    fOffsetx[3] = 0.;
    fSigmax[3] = 0.015;
    fOffsety[3] = 0.;
    fSigmay[3] = 0.015;
  }

  // angle cuts
  fOffx[0] = -0.004725;
  fSigx[0] = 0.00053;
  fOffy[0] = 0.002465;
  fSigy[0] = 0.00046;

  fOffx[1] = -0.004725;
  fSigx[1] = 0.00053;
  fOffy[1] = 0.002465;
  fSigy[1] = 0.00046;

  fOffx[2] = -0.001189;
  fSigx[2] = 0.0003856;
  fOffy[2] = 0.001232;
  fSigy[2] = 0.00029;

  fOffx[3] = -0.0003465;
  fSigx[3] = 0.000396;
  fOffy[3] = 0.0012;
  fSigy[3] = 0.00024;
}

//_______________________________________________________________________________________
void Na61ArmParameters::SetCAMembers(double* meanLocX_Vds1, double* sigmaLocX_Vds1, double* meanLocY_Vds1, double* sigmaLocY_Vds1, double* meanLocX_Vds2, double* sigmaLocX_Vds2, double* meanLocY_Vds2, double* sigmaLocY_Vds2, double* meanLocX_Vds3, double* sigmaLocX_Vds3, double* meanLocY_Vds3, double* sigmaLocY_Vds3, double* meanLocX_Vds4, double* sigmaLocX_Vds4, double* meanLocY_Vds4, double* sigmaLocY_Vds4, double* meanX_Vds1, double* sigmaX_Vds1, double* meanY_Vds1, double* sigmaY_Vds1, double* meanX_Vds2, double* sigmaX_Vds2, double* meanY_Vds2, double* sigmaY_Vds2, double* meanX_Vds3, double* sigmaX_Vds3, double* meanY_Vds3, double* sigmaY_Vds3, double* meanX_Vds4, double* sigmaX_Vds4, double* meanY_Vds4, double* sigmaY_Vds4) {
  for (int i = 0; i < 8; i++) {
    fMeanLocX[0][i] = meanLocX_Vds1[i];
    fMeanLocY[0][i] = meanLocY_Vds1[i];
    fSigmaLocX[0][i] = sigmaLocX_Vds1[i];
    fSigmaLocY[0][i] = sigmaLocY_Vds1[i];
    fMeanLocX[1][i] = meanLocX_Vds2[i];
    fMeanLocY[1][i] = meanLocY_Vds2[i];
    fSigmaLocX[1][i] = sigmaLocX_Vds2[i];
    fSigmaLocY[1][i] = sigmaLocY_Vds2[i];
    fMeanLocX[2][i] = meanLocX_Vds3[i];
    fMeanLocY[2][i] = meanLocY_Vds3[i];
    fSigmaLocX[2][i] = sigmaLocX_Vds3[i];
    fSigmaLocY[2][i] = sigmaLocY_Vds3[i];
    fMeanLocX[3][i] = meanLocX_Vds4[i];
    fMeanLocY[3][i] = meanLocY_Vds4[i];
    fSigmaLocX[3][i] = sigmaLocX_Vds4[i];
    fSigmaLocY[3][i] = sigmaLocY_Vds4[i];

    fMeanX[0][i] = meanX_Vds1[i];
    fMeanY[0][i] = meanY_Vds1[i];
    fSigmaX[0][i] = sigmaX_Vds1[i];
    fSigmaY[0][i] = sigmaY_Vds1[i];
    fMeanX[1][i] = meanX_Vds2[i];
    fMeanY[1][i] = meanY_Vds2[i];
    fSigmaX[1][i] = sigmaX_Vds2[i];
    fSigmaY[1][i] = sigmaY_Vds2[i];
    fMeanX[2][i] = meanX_Vds3[i];
    fMeanY[2][i] = meanY_Vds3[i];
    fSigmaX[2][i] = sigmaX_Vds3[i];
    fSigmaY[2][i] = sigmaY_Vds3[i];
    fMeanX[3][i] = meanX_Vds4[i];
    fMeanY[3][i] = meanY_Vds4[i];
    fSigmaX[3][i] = sigmaX_Vds4[i];
    fSigmaY[3][i] = sigmaY_Vds4[i];
  }
}

/*
//____________________________________________________________________
ostream& operator<< (ostream & os,Na61ArmParameters *armpars)
{
  os << " Parameters for "
     <<armpars->GetArmName()<< " arm "
     <<endl;
  return os;
}
*/
