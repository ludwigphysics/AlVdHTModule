//____________________________________________________________________
//
// Interfase to set and acces arm geometry and matching
// parameters

#include "Na61VdParameters.h"
#include <utl/ErrorLogger.h>

#ifndef ROOT_TMath
#include "TMath.h"
#endif
#ifndef ROOT_TH2F
#include "TH2F.h"
#endif
#ifndef ROOT_TGraphErrors
#include "TGraphErrors.h"
#endif
#ifndef ROOT_TList
#include "TList.h"
#endif

#ifndef __IOSTREAM__
#include <iostream>
#endif

#include <sstream>

using std::endl;
using std::cout;
using std::cin;


using namespace utl;
using namespace std;

//____________________________________________________________________
// ClassImp(Na61VdParameters);

//____________________________________________________________________
Na61VdParameters::Na61VdParameters(TString /*armname*/) {
  fInit = false;
  fMatchParamsFile = 0;
  // normal constructor
  fOffAx_J = 0.0;
  fSigAx_J = 0.0;
  fOffAy_J = 0.0;
  fSigAy_J = 0.0;
  fOffAx_S = 0.0;
  fSigAx_S = 0.0;
  fOffAy_S = 0.0;
  fSigAy_S = 0.0;
  fVtxDzSigma = 0.075;
  fCommonVdSystemSetup = false;

  fBeamSpotOffsetX = 0.0;
  fBeamSpotSigmaX = 1.5;
  fBeamSpotOffsetY = 0.0;
  fBeamSpotSigmaY = 1.5;

  // default values
  fDevOffx[0] = -0.00020;
  fDevSigx[0] = 0.0143933;
  fDevOffy[0] = 0.00014;
  fDevSigy[0] = 0.0058761;
  fDevOffx[1] = 0.00377;
  fDevSigx[1] = 0.0353684;
  fDevOffy[1] = 0.00119;
  fDevSigy[1] = 0.0097333;
  fDevOffx[2] = 0.00033;
  fDevSigx[2] = 0.0396898;
  fDevOffy[2] = 0.00009;
  fDevSigy[2] = 0.0083968;
  fDevOffx[3] = -0.00004;
  fDevSigx[3] = 0.0175264;
  fDevOffy[3] = 0.00057;
  fDevSigy[3] = 0.0072806;

  fResOffx[0] = 0.00001;
  fResSigx[0] = 0.0011874;
  fResOffy[0] = -0.00001;
  fResSigy[0] = 0.0050000;
  fResOffx[1] = -0.00003;
  fResSigx[1] = 0.0035086;
  fResOffy[1] = -0.00004;
  fResSigy[1] = 0.0050000;
  fResOffx[2] = 0.00003;
  fResSigx[2] = 0.0035965;
  fResOffy[2] = 0.00005;
  fResSigy[2] = 0.0052853;
  fResOffx[3] = -0.00001;
  fResSigx[3] = 0.0012737;
  fResOffy[3] = -0.00002;
  fResSigy[3] = 0.0050000;
}

//____________________________________________________________________
Na61VdParameters::~Na61VdParameters() {}

//____________________________________________________________________
void Na61VdParameters::Init() {
  fInit = true;
  
  //  if(fRunId<300){

  if (fRunId < 10000) { // tune this number
    SetupCommonVdSystem_pPb2022();
  } else if (fRunId < 27514) {
    SetupCommonVdSystem();
  } else if (fRunId < 31974) {
    SetupCommonVdSystem_pPb();
    //  }else if(fRunId<2228){
  } else if (fRunId < 33450) { 
    SetupCommonVdSystem_xela150();
    //  }else if(fRunId<2259){
  } else if (fRunId < 33800) {
    SetupCommonVdSystem_xela75();
  } else if (fRunId < 37797) {
    SetupCommonVdSystem_xela40();
  } else if (fRunId < 38014) {
    SetupCommonVdSystem_pPb2018();
  } else {
    SetupCommonVdSystem_pbpb();
  }
  
  SetupCommonNA61System();
  
  MakeCompansationForRotation();

  ReadMatchParams();
}

//____________________________________________________________
void Na61VdParameters::ReadMatchParams() {
  if (!fMatchParamsFile) {
    WARNING (" Na61VdParameters::ReadMatchParams: matching params file points to null, check it out");
    return;
  }
    ostringstream info_init;
    info_init << "Na61VdParameters::ReadMatchParams: matching params file name " << fMatchParamsFile->GetName();
    INFO(info_init);

  // TF1* mean=0;
  // TF1* sigma=0;
  /*
  SetupFunctions("dxVsZ_JJ_1",mean,sigma);
  //fMean_dx_JJ_1=mean; fSigma_dx_JJ_1=sigma;
  SetupFunctions("dxVsZ_JS_1",fMean_dx_JS_1,fSigma_dx_JS_1);
  SetupFunctions("dxVsZ_SS_1",fMean_dx_SS_1,fSigma_dx_SS_1);
  SetupFunctions("dxVsZ_SJ_1",fMean_dx_SJ_1,fSigma_dx_SJ_1);

  SetupFunctions("dxVsZ_JJ_2",fMean_dx_JJ_2,fSigma_dx_JJ_2);
  SetupFunctions("dxVsZ_JS_2",fMean_dx_JS_2,fSigma_dx_JS_2);
  SetupFunctions("dxVsZ_SS_2",fMean_dx_SS_2,fSigma_dx_SS_2);
  SetupFunctions("dxVsZ_SJ_2",fMean_dx_SJ_2,fSigma_dx_SJ_2);

  SetupFunctions("dyVsY_JJ",fMean_dy_JJ,fSigma_dy_JJ);
  SetupFunctions("dyVsY_JS",fMean_dy_JS,fSigma_dy_JS);
  SetupFunctions("dyVsY_SS",fMean_dy_SS,fSigma_dy_SS);
  SetupFunctions("dyVsY_SJ",fMean_dy_SJ,fSigma_dy_SJ);
  */
  SetupFunctions("dxVsZ_JJ_1", 0);
  SetupFunctions("dxVsZ_JS_1", 1);
  SetupFunctions("dxVsZ_SS_1", 2);
  SetupFunctions("dxVsZ_SJ_1", 3);

  SetupFunctions("dxVsZ_JJ_2", 4);
  SetupFunctions("dxVsZ_JS_2", 5);
  SetupFunctions("dxVsZ_SS_2", 6);
  SetupFunctions("dxVsZ_SJ_2", 7);

  SetupFunctions("dyVsY_JJ_1", 0);
  SetupFunctions("dyVsY_JS_1", 1);
  SetupFunctions("dyVsY_SS_1", 2);
  SetupFunctions("dyVsY_SJ_1", 3);

  SetupFunctions("dyVsY_JJ_2", 4);
  SetupFunctions("dyVsY_JS_2", 5);
  SetupFunctions("dyVsY_SS_2", 6);
  SetupFunctions("dyVsY_SJ_2", 7);

  for (int i = 0; i < 8; i++) {
    double range_min;
    double range_max;
    fMean_dx[i]->GetRange(range_min, range_max);
    ostringstream info1;
    info1 << fMean_dx[i]->GetName() << ": range_min=" << range_min << " range_max=" << range_max;
    INFO(info1);
  }

  for (int i = 0; i < 8; i++) {
    double range_min;
    double range_max;
    fMean_dy[i]->GetRange(range_min, range_max);
    ostringstream info2;
    info2 << fMean_dy[i]->GetName() << ": range_min=" << range_min << " range_max=" << range_max;
    INFO(info2);
  }
}

//_____________________________________________________________
// void Na61VdParameters::SetupFunctions(TString histName,TF1* mean, TF1* sigma)
void Na61VdParameters::SetupFunctions(TString histName, int i) {
  TH2F* hh = (TH2F*)fMatchParamsFile->Get(histName.Data());
  if (!hh) {
	ostringstream warn;
    warn << "Na61VdParameters::ReadMatchParams: histogram " << histName.Data() << " not found, check it out" ;
    WARNING(warn);
    return;
  }

  TList* list = hh->GetListOfFunctions();

  TGraphErrors* Mean = (TGraphErrors*)list->At(0);
  TGraphErrors* Sigma = (TGraphErrors*)list->At(1);

  TF1* mean = (TF1*)(Mean->GetFunction(Form("meanFit")))->Clone(Form("%s_meanFit", hh->GetName()));
  // mean -> SetName(Form("%s_meanFit",hh->GetName()));
  // cout<<"mean="<<mean<<" "<<mean->GetName()<<endl;
  TF1* sigma = (TF1*)(Sigma->GetFunction(Form("sigmaFit")))->Clone(Form("%s_sigmaFit", hh->GetName()));
  // sigma -> SetName(Form("%s_sigmaFit",hh->GetName()));
  // cout<<"sigma="<<sigma<<" "<<sigma->GetName()<<endl;
  if (histName.Contains("dxVsZ")) {
    fMean_dx[i] = mean;
    fSigma_dx[i] = sigma;
  } else {
    fMean_dy[i] = mean;
    fSigma_dy[i] = sigma;
  }
}

//________________________________________________________________________________________________________
int Na61VdParameters::GetDxMatchParams(double z, double x, double xvd, double& mean, double& sigma) {
  int range;
  if (z < -220)
    range = 1;
  else
    range = 2;

  TF1* meanf = 0;
  TF1* sigf = 0;
  /*
  if(range==1){
    if(x>0 && xvd>0){meanf=fMean_dx_JJ_1; sigf=fSigma_dx_JJ_1;}
    if(x>0 && xvd<0){meanf=fMean_dx_JS_1; sigf=fSigma_dx_JS_1;}
    if(x<0 && xvd<0){meanf=fMean_dx_SS_1; sigf=fSigma_dx_SS_1;}
    if(x<0 && xvd>0){meanf=fMean_dx_SJ_1; sigf=fSigma_dx_SJ_1;}
  }else{
    if(x>0 && xvd>0){meanf=fMean_dx_JJ_2; sigf=fSigma_dx_JJ_2;}
    if(x>0 && xvd<0){meanf=fMean_dx_JS_2; sigf=fSigma_dx_JS_2;}
    if(x<0 && xvd<0){meanf=fMean_dx_SS_2; sigf=fSigma_dx_SS_2;}
    if(x<0 && xvd>0){meanf=fMean_dx_SJ_2; sigf=fSigma_dx_SJ_2;}
  }
  */

  // cout<<" GetDxMatchParams: range="<<range<<" x="<<x<<"  xvd="<<xvd<<endl;

  if (range == 1) {
    if (x > 0 && xvd > 0) {
      meanf = fMean_dx[0];
      sigf = fSigma_dx[0];
    }
    if (x > 0 && xvd < 0) {
      meanf = fMean_dx[1];
      sigf = fSigma_dx[1];
    }
    if (x < 0 && xvd < 0) {
      meanf = fMean_dx[2];
      sigf = fSigma_dx[2];
    }
    if (x < 0 && xvd > 0) {
      meanf = fMean_dx[3];
      sigf = fSigma_dx[3];
    }
  } else {
    if (x > 0 && xvd > 0) {
      meanf = fMean_dx[4];
      sigf = fSigma_dx[4];
    }
    if (x > 0 && xvd < 0) {
      meanf = fMean_dx[5];
      sigf = fSigma_dx[5];
    }
    if (x < 0 && xvd < 0) {
      meanf = fMean_dx[6];
      sigf = fSigma_dx[6];
    }
    if (x < 0 && xvd > 0) {
      meanf = fMean_dx[7];
      sigf = fSigma_dx[7];
    }
  }

  // cout<<" range="<<range<<" "<<meanf<<endl;
  if ((!meanf) || (!sigf)) return 0;

  double range_min;
  double range_max;

  meanf->GetRange(range_min, range_max);

  // cout<<" range_min="<<range_min<<" range_max="<<range_max<<endl;

  // int ii;
  // cin>>ii;

  if ((z > range_min) && (z < range_max)) {
    mean = meanf->Eval(z);
    sigma = sigf->Eval(z);
  } else {  // in this case reject all matchings
    mean = 0;
    sigma = 0;
  }
  return 1;
}

//________________________________________________________________________________________________________
int Na61VdParameters::GetDyMatchParams(double z, double y, double x, double xvd, double& mean, double& sigma) {
  if (y < -20) y = -20.;
  if (y > 20) y = 20.;

  // range = 1: Matching with TPC1
  // range = 2: Matching with TPC2
  int range;
  if (z < -220)
    range = 1;
  else
    range = 2;

  TF1* meanf = 0;
  TF1* sigf = 0;
  /*
  if(x>0 && xvd>0){meanf=fMean_dy_JJ; sigf=fSigma_dy_JJ;}
  if(x>0 && xvd<0){meanf=fMean_dy_JS; sigf=fSigma_dy_JS;}
  if(x<0 && xvd>0){meanf=fMean_dy_SJ; sigf=fSigma_dy_SJ;}
  if(x<0 && xvd<0){meanf=fMean_dy_SS; sigf=fSigma_dy_SS;}
  */
  // if(x>0 && xvd>0){meanf=fMean_dy[0]; sigf=fSigma_dy[0];}
  // if(x>0 && xvd<0){meanf=fMean_dy[1]; sigf=fSigma_dy[1];}
  // if(x<0 && xvd>0){meanf=fMean_dy[2]; sigf=fSigma_dy[2];}
  // if(x<0 && xvd<0){meanf=fMean_dy[3]; sigf=fSigma_dy[3];}

  if (range == 1) {
    if (x > 0 && xvd > 0) {
      meanf = fMean_dy[0];
      sigf = fSigma_dy[0];
    }
    if (x > 0 && xvd < 0) {
      meanf = fMean_dy[1];
      sigf = fSigma_dy[1];
    }
    if (x < 0 && xvd < 0) {
      meanf = fMean_dy[2];
      sigf = fSigma_dy[2];
    }
    if (x < 0 && xvd > 0) {
      meanf = fMean_dy[3];
      sigf = fSigma_dy[3];
    }
  } else {
    if (x > 0 && xvd > 0) {
      meanf = fMean_dy[4];
      sigf = fSigma_dy[4];
    }
    if (x > 0 && xvd < 0) {
      meanf = fMean_dy[5];
      sigf = fSigma_dy[5];
    }
    if (x < 0 && xvd < 0) {
      meanf = fMean_dy[6];
      sigf = fSigma_dy[6];
    }
    if (x < 0 && xvd > 0) {
      meanf = fMean_dy[7];
      sigf = fSigma_dy[7];
    }
  }

  if ((!meanf) || (!sigf)) return 0;

  mean = meanf->Eval(y);
  sigma = sigf->Eval(y);

  return 1;
}

//_____________________________________________________________
void Na61VdParameters::SetupCommonNA61System() {
    ostringstream info;
    info << " Na61VdParameters::SetupCommonNA61System";
    INFO(info);
 
  
  
   if (fRunId < 31974) {
    // values based on runs from PbPb 2016 test data
    fOffsetToTpcX = -0.210206;  // 210206
    fOffsetToTpcY = 0.890052;
    fOffsetToTpcZ = -5974.2;
  }
  //  if(fRunId>1346){
  //  if(fRunId>33197){

  //  if(fRunId>33201){

  //  }

  else if (fRunId <= 33395) {
    // values based on hld1350 XeLa run
    // fOffsetToTpcX = 1.953;
    // fOffsetToTpcY = 3.315;
    // fOffsetToTpcZ = -5971.182 + fdOffZ; // tunning VD z position

    // fOffsetToTpcX = 1.957;
    // fOffsetToTpcY = 3.279;
    // fOffsetToTpcZ = -5981.39 + fdOffZ;
    //
    // based on runs from XeLa:33201-33375
    // fOffsetToTpcX = 2.33165;
    // fOffsetToTpcY = 0.89167;
    // fOffsetToTpcZ = -5981.43;
    
    //use only one value for the key 17_034 and 17_039
    fOffsetToTpcX = 2.22342;     fOffsetToTpcY = -0.851794;     fOffsetToTpcZ = -5980.5;

/*    fOffsetToTpcX = 1.59849;
    fOffsetToTpcY = 5.43525;
    fOffsetToTpcZ = -5981.91;
    
    fOffsetToTpcX = 2.2266;     fOffsetToTpcY = -0.96999;     fOffsetToTpcZ = -5981.11;
    
         if (fRunId >= 33375) { fOffsetToTpcX = 2.21137;     fOffsetToTpcY = -0.859229;     fOffsetToTpcZ = -5981.12;}
//    else if (fRunId >= 33368) {fOffsetToTpcX = 1.61172;     fOffsetToTpcY = 5.43722;     fOffsetToTpcZ = -5981.84;}
    else if (fRunId >= 33358) { fOffsetToTpcX = 2.20797;     fOffsetToTpcY = -0.794247;     fOffsetToTpcZ = -5981.18;}
    else if (fRunId >= 33343) {fOffsetToTpcX = 2.21587;     fOffsetToTpcY = -0.762542;     fOffsetToTpcZ = -5981.16;}
//    else if (fRunId >= 33341) {fOffsetToTpcX = 1.54207;     fOffsetToTpcY = 5.43516;     fOffsetToTpcZ = -5981.86;}
    else if (fRunId >= 33320) {fOffsetToTpcX = 2.21478;     fOffsetToTpcY = -0.759402;     fOffsetToTpcZ = -5981.14;}
//    else if (fRunId >= 33310) {fOffsetToTpcX = 1.52265;     fOffsetToTpcY = 5.43313;     fOffsetToTpcZ = -5981.86;}
//    else if (fRunId >= 33309) {fOffsetToTpcX = 1.5354;     fOffsetToTpcY = 5.43526;     fOffsetToTpcZ = -5981.81; }
    else if (fRunId >= 33307) {fOffsetToTpcX = 2.21457;     fOffsetToTpcY = -0.760532;     fOffsetToTpcZ = -5981.17;}
//    else if (fRunId >= 33303) {fOffsetToTpcX = 1.55805;     fOffsetToTpcY = 5.43473;     fOffsetToTpcZ = -5981.83;}
//    else if (fRunId >= 33302) {fOffsetToTpcX = 1.58512;     fOffsetToTpcY = 5.4381;     fOffsetToTpcZ = -5981.85; }
//    else if (fRunId >= 33301) {fOffsetToTpcX = 1.59044;     fOffsetToTpcY = 5.43385;     fOffsetToTpcZ = -5981.89;}
    else if (fRunId >= 33300) {fOffsetToTpcX = 2.2247;     fOffsetToTpcY = -0.815784;     fOffsetToTpcZ = -5981.12;}
    else if (fRunId >= 33275) {fOffsetToTpcX = 2.21168;     fOffsetToTpcY = -0.815103;     fOffsetToTpcZ = -5981.13;}
//    else if (fRunId >= 33273) {fOffsetToTpcX = 1.58745;     fOffsetToTpcY = 5.43306;     fOffsetToTpcZ = -5981.85;}
//    else if (fRunId >= 33266) {fOffsetToTpcX = 1.59309;     fOffsetToTpcY = 5.43555;     fOffsetToTpcZ = -5981.84;}
    else if (fRunId >= 33264) {fOffsetToTpcX = 2.20816;     fOffsetToTpcY = -0.810237;     fOffsetToTpcZ = -5981.16;}
//    else if (fRunId >= 33259) {fOffsetToTpcX = 1.58422;     fOffsetToTpcY = 5.43298;     fOffsetToTpcZ = -5981.89;}
    else if (fRunId >= 33227) {fOffsetToTpcX = 2.21379;     fOffsetToTpcY = -0.797195;     fOffsetToTpcZ = -5981.17;}
    else if (fRunId >= 33201) {fOffsetToTpcX = 2.21629;     fOffsetToTpcY = -0.79555;     fOffsetToTpcZ = -5981.13;}
    else if (fRunId >= 33187) {fOffsetToTpcX = 2.21763;     fOffsetToTpcY = -0.86154;     fOffsetToTpcZ = -5981.2; }
//    else if (fRunId >= 33180) { fOffsetToTpcX = 1.5994;     fOffsetToTpcY = 5.43714;     fOffsetToTpcZ = -5981.81; }
    else if (fRunId >= 33172) {fOffsetToTpcX = 2.22269;     fOffsetToTpcY = -0.834532;     fOffsetToTpcZ = -5981.19;}
    else if (fRunId >= 33157) {fOffsetToTpcX = 2.21832;     fOffsetToTpcY = -0.823327;     fOffsetToTpcZ = -5981.17; }
    else if (fRunId >= 33128) {fOffsetToTpcX = 2.22459;     fOffsetToTpcY = -0.880074;     fOffsetToTpcZ = -5981.14;}
//    else if (fRunId >= 33124) {fOffsetToTpcX = 1.60767;     fOffsetToTpcY = 5.43743;     fOffsetToTpcZ = -5981.53;}
    else if (fRunId >= 33078) {fOffsetToTpcX = 2.23368;     fOffsetToTpcY = -0.870574;     fOffsetToTpcZ = -5981.16;}
    else if (fRunId >= 32975) {fOffsetToTpcX = 2.21659;     fOffsetToTpcY = -0.909716;     fOffsetToTpcZ = -5981.19;}
//    else if (fRunId >= 32969) {fOffsetToTpcX = 1.63172;     fOffsetToTpcY = 5.43443;     fOffsetToTpcZ = -5982.02;}
//    else if (fRunId >= 32960) {fOffsetToTpcX = 1.62925;     fOffsetToTpcY = 5.43574;     fOffsetToTpcZ = -5981.86;}
//    else if (fRunId >= 32958) {fOffsetToTpcX = 1.62376;     fOffsetToTpcY = 5.4369;     fOffsetToTpcZ = -5981.86; }
    else if (fRunId >= 32957) {fOffsetToTpcX = 2.21912;     fOffsetToTpcY = -0.914402;     fOffsetToTpcZ = -5981.18;}
    else if (fRunId >= 32910) {fOffsetToTpcX = 2.21848;     fOffsetToTpcY = -0.895926;     fOffsetToTpcZ = -5981.17;}
    else if (fRunId >= 32848) {fOffsetToTpcX = 2.21012;     fOffsetToTpcY = -0.871569;     fOffsetToTpcZ = -5981.22;}
    //else if (fRunId >= 32752) {fOffsetToTpcX = 1.54377;     fOffsetToTpcY = 5.4353;     fOffsetToTpcZ = -5981.94;}
    else if (fRunId >= 32747) {fOffsetToTpcX = 2.20581;     fOffsetToTpcY = -0.804447;     fOffsetToTpcZ = -5981.21; }
//    else if (fRunId >= 32745) {fOffsetToTpcX = 1.54767;     fOffsetToTpcY = 5.43368;     fOffsetToTpcZ = -5981.93;}
    else if (fRunId >= 32728) {fOffsetToTpcX = 2.20889;     fOffsetToTpcY = -0.836459;     fOffsetToTpcZ = -5981.19;}
    else if (fRunId >= 32601) {fOffsetToTpcX = 2.2623;     fOffsetToTpcY = -0.784905;     fOffsetToTpcZ = -5981.18; }
*/
    
//    if (fRunId == 33201) { //Test for comparing with pawel
//		fOffsetToTpcX = 1.5828;     fOffsetToTpcY = 5.43526;     fOffsetToTpcZ = -5981.79;
// double fOffsetToTpcX = 1.59849;
// double fOffsetToTpcY = 5.43525;
// double fOffsetToTpcZ = -5981.91;
//		fOffsetToTpcX = 1.59448;
//		fOffsetToTpcY = 5.46626;
//		fOffsetToTpcZ = -5981.12;
//		
//		fOffsetToTpcX += -0.030958; fOffsetToTpcY += 0.169536; fOffsetToTpcZ += -0.082227;


//	}	


  }  // XeLa150

  else if (fRunId <= 33756) { 
    fOffsetToTpcX = 1.53984;
	fOffsetToTpcY = 4.67189;
	fOffsetToTpcZ = -5982.08;
	
	     if (fRunId >= 33755) {    fOffsetToTpcX = 1.56771;     fOffsetToTpcY = 5.09061;     fOffsetToTpcZ = -5982.05; }
    else if (fRunId >= 33751) {    fOffsetToTpcX = 1.52813;     fOffsetToTpcY = 4.37293;     fOffsetToTpcZ = -5982.09; }
    else if (fRunId >= 33747) {    fOffsetToTpcX = 1.5632;      fOffsetToTpcY = 4.96332;     fOffsetToTpcZ = -5982.05; }
    else if (fRunId >= 33739) {    fOffsetToTpcX = 1.53193;     fOffsetToTpcY = 4.82134;     fOffsetToTpcZ = -5982.09; }
    else if (fRunId >= 33737) {    fOffsetToTpcX = 1.54576;     fOffsetToTpcY = 4.53964;     fOffsetToTpcZ = -5982.05; }
    else if (fRunId >= 33696) {    fOffsetToTpcX = 1.54516;     fOffsetToTpcY = 4.73728;     fOffsetToTpcZ = -5981.99; }
    else if (fRunId >= 33693) {    fOffsetToTpcX = 1.55769;     fOffsetToTpcY = 4.72978;     fOffsetToTpcZ = -5981.97; }
    else if (fRunId >= 33687) {    fOffsetToTpcX = 1.58092;     fOffsetToTpcY = 4.86701;     fOffsetToTpcZ = -5981.99; }
    else if (fRunId >= 33683) {    fOffsetToTpcX = 1.55184;     fOffsetToTpcY = 4.68174;     fOffsetToTpcZ = -5981.98; }
    else if (fRunId >= 33680) {    fOffsetToTpcX = 1.58938;     fOffsetToTpcY = 5.19055;     fOffsetToTpcZ = -5982.16; }
    else if (fRunId >= 33646) {    fOffsetToTpcX = 1.56247;     fOffsetToTpcY = 4.64971;     fOffsetToTpcZ = -5982; }
    else if (fRunId >= 33635) {    fOffsetToTpcX = 1.5439;      fOffsetToTpcY = 4.52673;     fOffsetToTpcZ = -5982.09; }
    else if (fRunId >= 33633) {    fOffsetToTpcX = 1.58542;     fOffsetToTpcY = 5.04801;     fOffsetToTpcZ = -5982.12; }
    else if (fRunId >= 33613) {    fOffsetToTpcX = 1.61493;     fOffsetToTpcY = 4.8458;      fOffsetToTpcZ = -5982.03; }
    else if (fRunId >= 33585) {    fOffsetToTpcX = 1.59898;     fOffsetToTpcY = 4.94792;     fOffsetToTpcZ = -5982.03; }
    else if (fRunId >= 33549) {    fOffsetToTpcX = 1.57332;     fOffsetToTpcY = 4.73207;     fOffsetToTpcZ = -5982.06; }
    else if (fRunId >= 33539) {    fOffsetToTpcX = 1.55337;     fOffsetToTpcY = 4.5253;      fOffsetToTpcZ = -5981.95; }
    else if (fRunId >= 33535) {    fOffsetToTpcX = 1.54355;     fOffsetToTpcY = 4.57968;     fOffsetToTpcZ = -5982.1; }
    else if (fRunId >= 33500) {    fOffsetToTpcX = 1.53626;     fOffsetToTpcY = 4.58832;     fOffsetToTpcZ = -5982.09; }
    else if (fRunId >= 33498) {    fOffsetToTpcX = 1.56211;     fOffsetToTpcY = 4.69284;     fOffsetToTpcZ = -5982; }
    else if (fRunId >= 33495) {    fOffsetToTpcX = 1.58802;     fOffsetToTpcY = 5.04438;     fOffsetToTpcZ = -5981.95; }
    else if (fRunId >= 33493) {    fOffsetToTpcX = 1.58998;     fOffsetToTpcY = 5.73008;     fOffsetToTpcZ = -5982.08; }
    else if (fRunId >= 33490) {    fOffsetToTpcX = 1.58518;     fOffsetToTpcY = 5.07012;     fOffsetToTpcZ = -5982.01; }
    else if (fRunId >= 33484) {    fOffsetToTpcX = 1.55271;     fOffsetToTpcY = 5.16611;     fOffsetToTpcZ = -5982.01; }
    else if (fRunId >= 33481) {    fOffsetToTpcX = 1.55006;     fOffsetToTpcY = 4.68328;     fOffsetToTpcZ = -5982.09; }
    else if (fRunId >= 33480) {    fOffsetToTpcX = 1.55785;     fOffsetToTpcY = 4.8839;      fOffsetToTpcZ = -5982.1; }
    else if (fRunId >= 33469) {    fOffsetToTpcX = 1.53449;     fOffsetToTpcY = 4.68248;     fOffsetToTpcZ = -5982.08; }
    else if (fRunId >= 33464) {    fOffsetToTpcX = 1.55049;     fOffsetToTpcY = 4.66047;     fOffsetToTpcZ = -5982.08; }
    else if (fRunId >= 33461) {    fOffsetToTpcX = 1.55486;     fOffsetToTpcY = 4.76897;     fOffsetToTpcZ = -5982.05; }
    else if (fRunId >= 33450) {    fOffsetToTpcX = 1.5501;      fOffsetToTpcY = 5.14503;     fOffsetToTpcZ = -5982; }//value for 33453


  } //XeLa75 
  else if (fRunId <=35269) { 
    fOffsetToTpcX = 0;
    fOffsetToTpcY = 0;
    fOffsetToTpcZ = -5980.00;
  } //XeLa40 
  else { 
    fOffsetToTpcX = -0.242153;
	fOffsetToTpcY = 2.29986;
	fOffsetToTpcZ = -5977.32;
	
			 if (fRunId >= 38945){	    fOffsetToTpcX = 0.149831;     fOffsetToTpcY = 1.54769;     fOffsetToTpcZ = -5978.16;
	  } else if (fRunId >= 38943){ 	    fOffsetToTpcX = 0.148544;     fOffsetToTpcY = 1.4242;     fOffsetToTpcZ = -5978.18;
	  } else if (fRunId >= 38941){ 	    fOffsetToTpcX = 0.138448;     fOffsetToTpcY = 1.59566;     fOffsetToTpcZ = -5978.18;
	  } else if (fRunId >= 38938){ 	    fOffsetToTpcX = 0.123318;     fOffsetToTpcY = 1.64526;     fOffsetToTpcZ = -5978.2;
	  } else if (fRunId >= 38919){ 	    fOffsetToTpcX = 0.251798;     fOffsetToTpcY = 1.38511;     fOffsetToTpcZ = -5978.24;
	  } else if (fRunId >= 38917){ 	    fOffsetToTpcX = 0.226406;     fOffsetToTpcY = 1.54213;     fOffsetToTpcZ = -5978.25;
	  } else if (fRunId >= 38898){ 	    fOffsetToTpcX = 0.284145;     fOffsetToTpcY = 1.42806;     fOffsetToTpcZ = -5978.23;
	  } else if (fRunId >= 38889){ 	    fOffsetToTpcX = 0.269447;     fOffsetToTpcY = 1.00559;     fOffsetToTpcZ = -5978.27;
	  } else if (fRunId >= 38881){ 	    fOffsetToTpcX = 0.265181;     fOffsetToTpcY = 0.939586;     fOffsetToTpcZ = -5978.26;
	  } else if (fRunId >= 38835){ 	    fOffsetToTpcX = 0.266779;     fOffsetToTpcY = 1.50422;     fOffsetToTpcZ = -5978.21;
	  } else if (fRunId >= 38826){ 	    fOffsetToTpcX = 0.275734;     fOffsetToTpcY = 1.33974;     fOffsetToTpcZ = -5978.23;
	  } else if (fRunId >= 38825){ 	    fOffsetToTpcX = 0.268776;     fOffsetToTpcY = 1.24695;     fOffsetToTpcZ = -5978.25;
	  } else if (fRunId >= 38813){ 	    fOffsetToTpcX = 0.263671;     fOffsetToTpcY = 1.41078;     fOffsetToTpcZ = -5978.23;
	  } else if (fRunId >= 38812){ 	    fOffsetToTpcX = 0.273133;     fOffsetToTpcY = 1.40806;     fOffsetToTpcZ = -5978.2;
	  } else if (fRunId >= 38804){ 	    fOffsetToTpcX = 0.293545;     fOffsetToTpcY = 1.35027;     fOffsetToTpcZ = -5978.21;
	  } else if (fRunId >= 38796){ 	    fOffsetToTpcX = 0.290082;     fOffsetToTpcY = 1.20181;     fOffsetToTpcZ = -5978.27;
	  } else if (fRunId >= 38775){ 	    fOffsetToTpcX = 0.283052;     fOffsetToTpcY = 1.29737;     fOffsetToTpcZ = -5978.24;
	  } else if (fRunId >= 38772){ 	    fOffsetToTpcX = 0.275657;     fOffsetToTpcY = 1.32312;     fOffsetToTpcZ = -5978.24;
	  } else if (fRunId >= 38771){ 	    fOffsetToTpcX = 0.288253;     fOffsetToTpcY = 1.28664;     fOffsetToTpcZ = -5978.27;
	  } else if (fRunId >= 38764){ 	    fOffsetToTpcX = 0.288513;     fOffsetToTpcY = 1.21698;     fOffsetToTpcZ = -5978.27;
	  } else if (fRunId >= 38761){ 	    fOffsetToTpcX = 0.275344;     fOffsetToTpcY = 1.28682;     fOffsetToTpcZ = -5978.25;
	  } else if (fRunId >= 38759){ 	    fOffsetToTpcX = 0.282554;     fOffsetToTpcY = 1.23;     fOffsetToTpcZ = -5978.25;
	  } else if (fRunId >= 38739){ 	    fOffsetToTpcX = 0.324443;     fOffsetToTpcY = 1.32323;     fOffsetToTpcZ = -5978.24;
	  } else if (fRunId >= 38707){ 	    fOffsetToTpcX = 0.235069;     fOffsetToTpcY = 1.58609;     fOffsetToTpcZ = -5978.34;
	  } else if (fRunId >= 38667){ 	    fOffsetToTpcX = 0.261659;     fOffsetToTpcY = 1.56873;     fOffsetToTpcZ = -5978.31;
	  } else if (fRunId >= 38662){ 	    fOffsetToTpcX = 0.284008;     fOffsetToTpcY = 1.52695;     fOffsetToTpcZ = -5978.3;
	  } else if (fRunId >= 38657){ 	    fOffsetToTpcX = 0.277586;     fOffsetToTpcY = 1.17125;     fOffsetToTpcZ = -5978.3;
	  } else if (fRunId >= 38656){ 	    fOffsetToTpcX = 0.261275;     fOffsetToTpcY = 1.55609;     fOffsetToTpcZ = -5978.35;
	  } else if (fRunId >= 38654){ 	    fOffsetToTpcX = 0.258153;     fOffsetToTpcY = 1.54007;     fOffsetToTpcZ = -5978.33;
	  } else if (fRunId >= 38649){ 	    fOffsetToTpcX = 0.24163;     fOffsetToTpcY = 1.50944;     fOffsetToTpcZ = -5978.36;
	  } else if (fRunId >= 38645){ 	    fOffsetToTpcX = 0.25871;     fOffsetToTpcY = 1.39792;     fOffsetToTpcZ = -5978.32;
	  } else if (fRunId >= 38636){ 	    fOffsetToTpcX = 0.274783;     fOffsetToTpcY = 1.40469;     fOffsetToTpcZ = -5978.29;
	  } else if (fRunId >= 38628){ 	    fOffsetToTpcX = 0.244837;     fOffsetToTpcY = 1.59091;     fOffsetToTpcZ = -5978.33;
	  } else if (fRunId >= 38619){ 	    fOffsetToTpcX = 0.254306;     fOffsetToTpcY = 1.35051;     fOffsetToTpcZ = -5978.37;
	  //} else if (fRunId >= 38606){ 	     
	  } else if (fRunId >= 38599){ 	    fOffsetToTpcX = 0.221592;     fOffsetToTpcY = 1.16844;     fOffsetToTpcZ = -5978.23;
	  } else if (fRunId >= 38593){ 	    fOffsetToTpcX = 0.240286;     fOffsetToTpcY = 1.35939;     fOffsetToTpcZ = -5978.39;
	  } else if (fRunId >= 38568){ 	    fOffsetToTpcX = 0.279278;     fOffsetToTpcY = 0.929583;     fOffsetToTpcZ = -5978.44;
	  } else if (fRunId >= 38562){ 	    fOffsetToTpcX = 0.265749;     fOffsetToTpcY = 1.08974;     fOffsetToTpcZ = -5978.45;
	  } else if (fRunId >= 38548){ 	    fOffsetToTpcX = 0.235924;     fOffsetToTpcY = 0.212603;     fOffsetToTpcZ = -5978.63;
	  } else if (fRunId >= 38450){ 	    fOffsetToTpcX = 0.222635;     fOffsetToTpcY = -0.199241;     fOffsetToTpcZ = -5978.56;
	//  } else if (fRunId >= 38444){ 	    
	  } else if (fRunId >= 38436){ 	    //fOffsetToTpcX = 0.297693;     fOffsetToTpcY = 0.411648;     fOffsetToTpcZ = -5978.76;
										fOffsetToTpcX = 0.237976;     fOffsetToTpcY = -0.177021;     fOffsetToTpcZ = -5978.51;
	  }

	
  }  // PbPb150 2018
}

//_____________________________________________________________
void Na61VdParameters::SetupCommonVdSystem() {
  fCommonVdSystemSetup = true;
  
    ostringstream info;
    info << " Na61VdParameters::SetupCommonVdSystem: setting common parameters, fN=" << fN;
    INFO(info);
 
  // parameters used to tune relative arms orientation and
  // describe track/hits parameters in the common VD system

  // consistent setup (peaks for flag0 moved to their average positions)

  fOffAx_J = 0.026;
  fOffAx_S = -0.035;

  // base on run 168, moves flag0 peaks to 0
  fRotX_J = 0.000664844;
  fRotY_J = -0.00370569 + fN * 0.0005;

  fRotX_S = 0.00118727;
  fRotY_S = -0.00907036 - fN * 0.0005;

  // PS geometry
  // double offsetX_J = 0.5395 + 8.;
  // double offsetY_J =  0.025;
  // double offsetZ_J =  0.3929;

  // double offsetX_S = -0.5395 - 8.;
  // double offsetY_S = -0.025;
  // double offsetZ_S = -0.3929;

  // Wojtek geometry mv1
  double offsetX_J = 0.545375 + 8.;
  double offsetY_J = 0.0259859;
  double offsetZ_J = 0.474669;

  double offsetX_S = -0.545375 - 8.;
  double offsetY_S = -0.0259859;
  double offsetZ_S = -0.474669;

  //  if(fRunId<183){
  if (fRunId < 27469) {
    offsetX_J = offsetX_J - 0.0851116;
    offsetY_J = offsetY_J + 3.34109e-05;
    offsetZ_J = offsetZ_J + 0.00235553;

    offsetX_S = offsetX_S + 0.0851116;
    offsetY_S = offsetY_S - 3.34109e-05;
    offsetZ_S = offsetZ_S - 0.00235553;

    // needed for mv1 arm rotations (see above)
    offsetX_J = offsetX_J + 0.130234;
    offsetY_J = offsetY_J - 0.012593;
    offsetZ_J = offsetZ_J + 0.0676186;

    offsetX_S = offsetX_S - 0.130234;
    offsetY_S = offsetY_S + 0.012593;
    offsetZ_S = offsetZ_S - 0.0676186;

    // new version of local to global
    offsetX_J += 0.200084;
    offsetY_J += -0.0193869;
    offsetZ_J += -0.00322695;

    offsetX_S += -0.200084;
    offsetY_S += 0.0193869;
    offsetZ_S += 0.00322695;

    double ddx = 0;
    double ddz = 0;

    if (fN == 1) ddx = 0.0620283;
    if (fN == 2) ddx = 0.124211;
    if (fN == 3) ddx = 0.186093;
    if (fN == 4) {
      ddx = 0.248144;
      ddz = 0.00221872;
    }
    if (fN == 5) ddx = 0.310064;
    if (fN == 6) {
      ddx = 0.372192;
      ddz = 0.00448074;
    }
    if (fN == 7) {
      ddx = 0.434139;
      ddz = 0.00457992;
    }
    if (fN == 8) {
      ddx = 0.495959;
      ddz = 0.00536705;
    }
    if (fN == 9) {
      ddx = 0.557922;
      ddz = 0.00521475;
    }
    if (fN == 10) {
      ddx = 0.619927;
      ddz = 0.00521016;
    }

    if (fN == -1) {
      ddx = -0.0613839;
      ddz = -0.00118034;
    }
    if (fN == -2) {
      ddx = -0.123528;
      ddz = -0.00138696;
    }
    if (fN == -3) {
      ddx = -0.185356;
      ddz = -0.00247079;
    }

    offsetX_J += ddx;
    offsetX_S += -ddx;
    offsetZ_J += ddz;
    offsetZ_S += -ddz;
  }

  fJuraArmOffset.SetX(offsetX_J);
  fJuraArmOffset.SetY(offsetY_J);
  fJuraArmOffset.SetZ(offsetZ_J);
  fSaleveArmOffset.SetX(offsetX_S);
  fSaleveArmOffset.SetY(offsetY_S);
  fSaleveArmOffset.SetZ(offsetZ_S);
}

//_____________________________________________________________
void Na61VdParameters::SetupCommonVdSystem_pPb() {
  fCommonVdSystemSetup = true;
    ostringstream info;
    info << " Na61VdParameters::SetupCommonVdSystem_pPb: setting common parameters";
    INFO(info);

  // parameters used to tune relative arms orientation and
  // describe track/hits parameters in the common VD system

  fOffAx_J = 0.005989;
  fSigAx_J = 0.0004415;
  fOffAy_J = 0.001662;
  fSigAy_J = 0.0005088;
  fOffAx_S = 0.01202;
  fSigAx_S = 0.0004205;
  fOffAy_S = 0.001563;
  fSigAy_S = 0.0004873;

  // base on run 633 (Saleve) and 635 (Jura)
  fRotX_J = 0.001495;
  fRotY_J = 0.0035155;

  fRotX_S = 0.00165;
  fRotY_S = 0.009315;

  double offsetX_J = 0.545375 + 8.;
  double offsetY_J = 0.0259859;
  double offsetZ_J = 0.474669;

  double offsetX_S = -0.545375 - 8.;
  double offsetY_S = -0.0259859;
  double offsetZ_S = -0.474669;

  // new version of local to global
  offsetX_J += 0.;
  offsetY_J += 0.;
  offsetZ_J += 0.;

  offsetX_S += 0.;
  offsetY_S += 0.;
  offsetZ_S += 0.;

  fJuraArmOffset.SetX(offsetX_J);
  fJuraArmOffset.SetY(offsetY_J);
  fJuraArmOffset.SetZ(offsetZ_J);
  fSaleveArmOffset.SetX(offsetX_S);
  fSaleveArmOffset.SetY(offsetY_S);
  fSaleveArmOffset.SetZ(offsetZ_S);
}

//_____________________________________________________________
void Na61VdParameters::MakeCompansationForRotation() {
    ostringstream info_init;
    info_init << " Na61VdParameters::FindCompensationForRotation: entered the method";
    INFO(info_init);
  if (!fCommonVdSystemSetup) 
        WARNING("Na61VdParameters::MakeCompansationForRotation: SetupCommonVdSystem should be called before MakeCompansationForRotation.");

  ///////////////////////////////// rotate Vds1_0 in jura arm
  double sa = TMath::Sin(fRotZ_J);
  double ca = TMath::Cos(fRotZ_J);
  double sb = TMath::Sin(fRotY_J);
  double cb = TMath::Cos(fRotY_J);
  double sg = TMath::Sin(fRotX_J);
  double cg = TMath::Cos(fRotX_J);

  double x0 = 0;
  double y0 = 0;
  double z0 = 0;

  z0 = z0 - 75.0;

  // double x1J =      cb * x0                +     sb * z0;
  // double y1J =  -sg*sb * x0   +   cg * y0  +  sg*cb * z0;
  // double z1J =  -cg*sb * x0   -   sg * y0  +  cg*cb * z0;
  // ZXY
  double x1J = (ca * cb - sa * sg * sb) * x0 + sa * cg * y0 + (ca * sb + sa * sg * cb) * z0;
  double y1J = -(sa * cb + ca * sg * sb) * x0 + ca * cg * y0 + (-sa * sb + ca * sg * cb) * z0;
  double z1J = -cg * sb * x0 - sg * y0 + cg * cb * z0;

  z1J = z1J + 75.0;

  ///////////////////////////////// rotate Vds1_0 in saleve arm
  sa = TMath::Sin(fRotZ_S);
  ca = TMath::Cos(fRotZ_S);
  sb = TMath::Sin(fRotY_S);
  cb = TMath::Cos(fRotY_S);
  sg = TMath::Sin(fRotX_S);
  cg = TMath::Cos(fRotX_S);

  x0 = 0;
  y0 = 0;
  z0 = 0;

  z0 = z0 - 75.0;

  // double x1S =      cb * x0                +     sb * z0;
  // double y1S =  -sg*sb * x0   +   cg * y0  +  sg*cb * z0;
  // double z1S =  -cg*sb * x0   -   sg * y0  +  cg*cb * z0;
  // ZXY
  double x1S = (ca * cb - sa * sg * sb) * x0 + sa * cg * y0 + (ca * sb + sa * sg * cb) * z0;
  double y1S = -(sa * cb + ca * sg * sb) * x0 + ca * cg * y0 + (-sa * sb + ca * sg * cb) * z0;
  double z1S = -cg * sb * x0 - sg * y0 + cg * cb * z0;

  z1S = z1S + 75.0;

  double deltaX = -0.5 * (x1J + x1S);
  double deltaY = -0.5 * (y1J + y1S);
  double deltaZ = -0.5 * (z1J + z1S);

    ostringstream info;
    info << " correction shifts due to arm rotation: deltaX=" << deltaX << "  deltaY=" << deltaY << "  deltaZ=" << deltaZ;
    INFO(info);

  fJuraArmOffset.SetX(fJuraArmOffset.GetX() + deltaX);
  fJuraArmOffset.SetY(fJuraArmOffset.GetY() + deltaY);
  fJuraArmOffset.SetZ(fJuraArmOffset.GetZ() + deltaZ);

  fSaleveArmOffset.SetX(fSaleveArmOffset.GetX() + deltaX);
  fSaleveArmOffset.SetY(fSaleveArmOffset.GetY() + deltaY);
  fSaleveArmOffset.SetZ(fSaleveArmOffset.GetZ() + deltaZ);
}

//_____________________________________________________________
void Na61VdParameters::SetupCommonVdSystem_xela150() {
  fCommonVdSystemSetup = true;
    ostringstream info;
    info << " Na61VdParameters::SetupCommonVdSystem_xela150: setting common parameters. fN=" << fN << " fRunId=" << fRunId;
    INFO(info);

  // parameters used to tune relative arms orientation and
  // describe track/hits parameters in the common VD system

  // production(field) defauls values
  fBeamSpotOffsetX = -0.526628;
  fBeamSpotSigmaX = 0.592206;
  fBeamSpotOffsetY = -0.616405;
  fBeamSpotSigmaY = 0.802218;

  if (fRunId >= 31974) {
    fBeamSpotOffsetX = -0.405983;
    fBeamSpotSigmaX = 0.577560;
    fBeamSpotOffsetY = -0.615915;
    fBeamSpotSigmaY = 0.818117;
  }
  
   
  //  if(fRunId>=1494){ // 0 field runs
  if (fRunId >= 33396) {  // 0 field runs
    fBeamSpotOffsetX = 0.459624;
    fBeamSpotSigmaX = 0.565917;
    fBeamSpotOffsetY = 0.540213;
    fBeamSpotSigmaY = 0.885769;
  }


  fDevOffx[0] = -0.00020;
  fDevSigx[0] = 0.0143933;
  fDevOffy[0] = 0.00014;
  fDevSigy[0] = 0.0058761;
  fDevOffx[1] = 0.00377;
  fDevSigx[1] = 0.0353684;
  fDevOffy[1] = 0.00119;
  fDevSigy[1] = 0.0097333;
  fDevOffx[2] = 0.00033;
  fDevSigx[2] = 0.0396898;
  fDevOffy[2] = 0.00009;
  fDevSigy[2] = 0.0083968;
  fDevOffx[3] = -0.00004;
  fDevSigx[3] = 0.0175264;
  fDevOffy[3] = 0.00057;
  fDevSigy[3] = 0.0072806;


  fResOffx[0] = 0.00001;
  fResSigx[0] = 0.0012111;
  fResOffy[0] = 0.00002;
  fResSigy[0] = 0.0034618;
  fResOffx[1] = -0.00003;
  fResSigx[1] = 0.0035005;
  fResOffy[1] = -0.00006;
  fResSigy[1] = 0.0041643;
  fResOffx[2] = 0.00003;
  fResSigx[2] = 0.0035930;
  fResOffy[2] = 0.00002;
  fResSigy[2] = 0.0052658;
  fResOffx[3] = -0.00001;
  fResSigx[3] = 0.0012740;
  fResOffy[3] = -0.00000;
  fResSigy[3] = 0.0040724;


  fOffAx_J = 0.032;
  fOffAx_S = -0.0203;
  // base on run 633 (Saleve) and 635 (Jura)
  // fRotX_J =  0.001495;
  // fRotY_J =  0.0035155;
  fRotX_J = -0.001501;
  fRotY_J = -0.003418;
  // geometry after fine tunning
  fRotX_J += -0.0000521;
  fRotY_J += -0.0001430;

  // fRotX_S =  0.00165;
  // fRotY_S =  0.009315;
  fRotX_S = -0.001737;
  fRotY_S = -0.009973;
  // geometry after fine tunning
  fRotX_S += 0.0000780;
  fRotY_S += -0.0001214;

  fRotZ_J = fN * 0.001;
  fRotZ_S = -fRotZ_J;

  // no arm rotations use it for test (needs to adjust arm offsets)
  // fRotX_J =  0.0;
  // fRotY_J =  0.0;

  // fRotX_S =  0.0;
  // fRotY_S =  0.0;
  
  

  // correction based on run 1494 (0 field)
  double offsetX_J = 0.881427 + 8.;
  double offsetY_J = -0.0164678;
  double offsetZ_J = 0.317669;

   
  // new update based on run by run
if (fRunId >= 33375) {
    offsetX_J += 0.00128217;     offsetY_J += -0.025271;     offsetZ_J += 0.0289048;
  } else if (fRunId >= 33368) {
    offsetX_J += 0.00112151;     offsetY_J += -0.0253192;     offsetZ_J += 0.0289472;
  } else if (fRunId >= 33358) {
    offsetX_J += 0.00224058;     offsetY_J += -0.0246746;     offsetZ_J += 0.0306108;
  } else if (fRunId >= 33343) {
    offsetX_J += 0.00238484;     offsetY_J += -0.0246628;     offsetZ_J += 0.0307191;
  } else if (fRunId >= 33341) {
    offsetX_J += 0.00218607;     offsetY_J += -0.02467;     offsetZ_J += 0.0305837;
  } else if (fRunId >= 33320) {
    offsetX_J += 0.00226016;     offsetY_J += -0.0246586;     offsetZ_J += 0.0301061;
  } else if (fRunId >= 33310) {
    offsetX_J += 0.00164147;     offsetY_J += -0.0246931;     offsetZ_J += 0.0301449;
  } else if (fRunId >= 33309) {
    offsetX_J += 0.00166492;     offsetY_J += -0.024716;     offsetZ_J += 0.0303395;
  } else if (fRunId >= 33307) {
    offsetX_J += 0.00269887;     offsetY_J += -0.0246368;     offsetZ_J += 0.0311995;
  } else if (fRunId >= 33303) {
    offsetX_J += 0.00197913;     offsetY_J += -0.0247112;     offsetZ_J += 0.0304822;;
  } else if (fRunId >= 33302) {
    offsetX_J += 0.00244563;     offsetY_J += -0.0253826;     offsetZ_J += 0.0289756;
  } else if (fRunId >= 33301) {
    offsetX_J += 0.00214246;     offsetY_J += -0.0254197;     offsetZ_J += 0.0284729;
  } else if (fRunId >= 33300) {
    offsetX_J += 0.00288803;     offsetY_J += -0.0253319;     offsetZ_J += 0.0295281;
  } else if (fRunId >= 33275) {
    offsetX_J += 0.00246631;     offsetY_J += -0.0253993;     offsetZ_J += 0.0292016;
  } else if (fRunId >= 33273) {
    offsetX_J += 0.00261143;     offsetY_J += -0.0254091;     offsetZ_J += 0.0287028;
  } else if (fRunId >= 33266) {
    offsetX_J += 0.00229964;     offsetY_J += -0.0254382;     offsetZ_J += 0.0290513;
  } else if (fRunId >= 33264) {
    offsetX_J += 0.00318833;     offsetY_J += -0.025417;     offsetZ_J += 0.0289704;
  } else if (fRunId >= 33259) {
    offsetX_J += 0.0033224;     offsetY_J += -0.0254677;     offsetZ_J += 0.028689;
  } else if (fRunId >= 33227) {
    offsetX_J += 0.00287972;     offsetY_J += -0.0255281;     offsetZ_J += 0.0283812;
  } else if (fRunId >= 33201) {
    offsetX_J += 0.00268508;     offsetY_J += -0.0254855;     offsetZ_J += 0.0285192;
  } else if (fRunId >= 33187) {
    offsetX_J += 0.00210018;     offsetY_J += -0.025552;     offsetZ_J += 0.0290027;
  } else if (fRunId >= 33186) {
    offsetX_J += 0.00216623;     offsetY_J += -0.0255417;     offsetZ_J += 0.0287862;
  } else if (fRunId >= 33183) {
    offsetX_J += 0.00229574;     offsetY_J += -0.0255688;     offsetZ_J += 0.0289085;
  } else if (fRunId >= 33180) {
    offsetX_J += 0.00215417;     offsetY_J += -0.0255651;     offsetZ_J += 0.0280706;
  } else if (fRunId >= 33172) {
    offsetX_J += 0.00308461;     offsetY_J += -0.025579;     offsetZ_J += 0.0282401;
  } else if (fRunId >= 33157) {
    offsetX_J += 0.00285719;     offsetY_J += -0.0256479;     offsetZ_J += 0.0282264;
  } else if (fRunId >= 33128) {
    offsetX_J += 0.00250771;     offsetY_J += -0.0255489;     offsetZ_J += 0.0288506;
  } else if (fRunId >= 33124) {
    offsetX_J += 0.00318303;     offsetY_J += -0.0255627;     offsetZ_J += 0.0281225;
  } else if (fRunId >= 33081) {
    offsetX_J += 0.00290116;     offsetY_J += -0.0257138;     offsetZ_J += 0.0287909;
  } else if (fRunId >= 33078) {
    offsetX_J += 0.00245301;     offsetY_J += -0.0257106;     offsetZ_J += 0.0288674;
  } else if (fRunId >= 32975) {
    offsetX_J += 0.00401065;     offsetY_J += -0.0254636;     offsetZ_J += 0.0294202;
  } else if (fRunId >= 32969) {
    offsetX_J += 0.00379802;     offsetY_J += -0.0255272;     offsetZ_J += 0.0301258;
   } else if (fRunId >= 32960) {
    offsetX_J += 0.00377502;     offsetY_J += -0.0255921;     offsetZ_J += 0.0311593;
  } else if (fRunId >= 32958) {
    offsetX_J += 0.00349394;     offsetY_J += -0.0256591;     offsetZ_J += 0.0310115;
  } else if (fRunId >= 32957) {
    offsetX_J += 0.00338414;     offsetY_J += -0.0256846;     offsetZ_J += 0.0301215;
  } else if (fRunId >= 32910) {
    offsetX_J += 0.00419324;     offsetY_J += -0.0258889;     offsetZ_J += 0.0291774;
  } else if (fRunId >= 32752) {
   offsetX_J += 0.00441753;     offsetY_J += -0.0253909;     offsetZ_J += 0.0307668;
  } else if (fRunId >= 32747) {
    offsetX_J += 0.00497347;     offsetY_J += -0.0255868;     offsetZ_J += 0.0290195;
  } else if (fRunId >= 32745) {
    offsetX_J += 0.00401527;     offsetY_J += -0.0253959;     offsetZ_J += 0.0309033;
  } else if (fRunId >= 32728) {
    offsetX_J += 0.00412162;     offsetY_J += -0.0254333;     offsetZ_J += 0.0304183;
  } else if (fRunId >= 32601) {
    offsetX_J += 0.00716107;     offsetY_J += -0.0250877;     offsetZ_J += 0.0311456;
  }

 // new update based on run by run after HT
   if (fRunId >= 33375) {
    offsetX_J += -0.000754734;     offsetY_J += -1.72113e-05;     offsetZ_J += 0.00418867;
  } else if (fRunId >= 33368) {
    offsetX_J += -0.000603406;     offsetY_J += 1.59876e-05;     offsetZ_J += 0.00447399;
  } else if (fRunId >= 33358) {
    offsetX_J += -0.000551062;     offsetY_J += 4.69269e-05;     offsetZ_J += 0.00569247;
  } else if (fRunId >= 33343) {
    offsetX_J += -0.00053957;     offsetY_J += 5.59084e-05;     offsetZ_J += 0.00588227;
  } else if (fRunId >= 33341) {
    offsetX_J += -0.000569192;     offsetY_J += 4.15782e-05;     offsetZ_J += 0.00592216;
  } else if (fRunId >= 33320) {
    offsetX_J += -0.000449559;     offsetY_J += 4.0168e-05;     offsetZ_J += 0.00630633;
  } else if (fRunId >= 33310) {
    offsetX_J += -0.000446127;     offsetY_J += 3.8123e-05;     offsetZ_J += 0.00580463;
  } else if (fRunId >= 33309) {
    offsetX_J += -0.000514832;     offsetY_J += 4.96689e-05;     offsetZ_J += 0.00580262;
  } else if (fRunId >= 33307) {
    offsetX_J += -0.000393756;     offsetY_J += 6.68978e-05;     offsetZ_J += 0.00536586;
  } else if (fRunId >= 33303) {
    offsetX_J += -0.000352311;     offsetY_J += -7.00764e-05;     offsetZ_J += 0.00507896;
  } else if (fRunId >= 33302) {
    offsetX_J += -0.000715517;     offsetY_J += -7.91382e-06;     offsetZ_J += 0.00399279;
  } else if (fRunId >= 33301) {
    offsetX_J += -0.000718572;     offsetY_J += 2.96116e-05;     offsetZ_J += 0.0036907;
  } else if (fRunId >= 33300) {
    offsetX_J += -0.000742026;     offsetY_J += 9.55912e-06;     offsetZ_J += 0.00343541;
  } else if (fRunId >= 33275) {
    offsetX_J += -0.000743586;     offsetY_J += -5.16828e-06;     offsetZ_J += 0.00379721;
  } else if (fRunId >= 33273) {
    offsetX_J += -0.000820424;     offsetY_J += -2.339e-05;     offsetZ_J += 0.00338246;
  } else if (fRunId >= 33266) {
    offsetX_J += -0.000713392;     offsetY_J += -1.03141e-05;     offsetZ_J += 0.00342008;
  } else if (fRunId >= 33264) {
    offsetX_J += -0.00076032;     offsetY_J += -1.12007e-05;     offsetZ_J += 0.00407902;
  } else if (fRunId >= 33259) {
    offsetX_J += -0.000760514;     offsetY_J += 2.73745e-05;     offsetZ_J += 0.00385327;
  } else if (fRunId >= 33227) {
    offsetX_J += -0.000812323;     offsetY_J += 1.65406e-05;     offsetZ_J += 0.00327415;
  } else if (fRunId >= 33201) {
    offsetX_J += -0.00128398;     offsetY_J += -1.69029e-05;     offsetZ_J += 0.00304535;
  } else if (fRunId >= 33187) {
    offsetX_J += -0.000716509;     offsetY_J += -5.10677e-06;     offsetZ_J += 0.00298544;
  } else if (fRunId >= 33186) {
    offsetX_J += -0.00063818;     offsetY_J += 1.21992e-05;     offsetZ_J += 0.00296582;
  } else if (fRunId >= 33183) {
    offsetX_J += -0.000719839;     offsetY_J += 7.20623e-06;     offsetZ_J += 0.00272865;
  } else if (fRunId >= 33180) {
    offsetX_J += -0.000697076;     offsetY_J += -9.88275e-06;     offsetZ_J += 0.00313978;
  } else if (fRunId >= 33172) {
    offsetX_J += -0.000797581;     offsetY_J += 1.89149e-05;     offsetZ_J += 0.00336298;
  } else if (fRunId >= 33157) {
    offsetX_J += -0.00047063;     offsetY_J += 1.24231e-05;     offsetZ_J += 0.00290912;
  } else if (fRunId >= 33128) {
    offsetX_J += -0.000741315;     offsetY_J += -4.5088e-05;     offsetZ_J += 0.0019075;
  } else if (fRunId >= 33124) {
    offsetX_J += -0.000839092;     offsetY_J += -1.92448e-05;     offsetZ_J += 0.00211689;
  } else if (fRunId >= 33081) {
    offsetX_J += -0.000648234;     offsetY_J += 1.21928e-05;     offsetZ_J += 0.00214512;
  } else if (fRunId >= 33078) {
    offsetX_J += -0.000575431;     offsetY_J += 8.78754e-06;     offsetZ_J += 0.00205498;
  } else if (fRunId >= 32975) {
    offsetX_J += -0.000391493;     offsetY_J += -1.82293e-05;     offsetZ_J += 0.0029369;
  } else if (fRunId >= 32969) {
    offsetX_J += -0.000490485;     offsetY_J += -4.91912e-05;     offsetZ_J += 0.00166449;
  } else if (fRunId >= 32960) {
    offsetX_J += -0.00048789;     offsetY_J += -4.62598e-05;     offsetZ_J += 0.00140401;
  } else if (fRunId >= 32958) {
    offsetX_J += -0.000550923;     offsetY_J += -5.73429e-05;     offsetZ_J += 0.00131714;
  } else if (fRunId >= 32957) {
    offsetX_J += -0.000500016;     offsetY_J += -5.10513e-05;     offsetZ_J += 0.00168489;
  } else if (fRunId >= 32910) {
    offsetX_J += -0.000445711;     offsetY_J += -5.72918e-05;     offsetZ_J += 0.00232699;
  } else if (fRunId >= 32752) {
    offsetX_J += -0.000282085;     offsetY_J += 7.53324e-07;     offsetZ_J += 0.00361423;
  } else if (fRunId >= 32747) {
    offsetX_J += -0.000431607;     offsetY_J += 5.41508e-05;     offsetZ_J += 0.00356796;
  } else if (fRunId >= 32745) {
    offsetX_J += -0.000226846;     offsetY_J += -1.95409e-05;     offsetZ_J += 0.00369077;
  } else if (fRunId >= 32728) {
    offsetX_J += -0.000400267;     offsetY_J += -2.55886e-05;     offsetZ_J += 0.00337148;
  } else if (fRunId >= 32601) {
    offsetX_J += -0.000696135;     offsetY_J += -4.85312e-05;     offsetZ_J += 0.00329741;
  }



 
//for test comparing with pawel
//  if (fRunId == 33201) {
//   offsetX_J =  8.88249;
//   offsetY_J = -0.0419041;
//   offsetZ_J =  0.349625;
//   
//    offsetX_J += 0.00122523;offsetY_J += -2.96085e-05; offsetZ_J += 0.000113689;
//  }	  
 
	double offsetX_S =-offsetX_J;
	double offsetY_S =-offsetY_J;
	double offsetZ_S =-offsetZ_J;
	
	ostringstream info_fin;
    info_fin << "Offset values: Jura: "<<offsetX_J<<" "<<offsetY_J<<" "<<offsetZ_J<<"   Saleve: "<<offsetX_S<<" "<<offsetY_S<<" "<<offsetZ_S;
    INFO(info_fin);

  	
  // parameter to extract Jura and Saleve events with consistent vertexes.
  fVtxDzSigma = 0.075;

  fJuraArmOffset.SetX(offsetX_J);
  fJuraArmOffset.SetY(offsetY_J);
  fJuraArmOffset.SetZ(offsetZ_J);
  fSaleveArmOffset.SetX(offsetX_S);
  fSaleveArmOffset.SetY(offsetY_S);
  fSaleveArmOffset.SetZ(offsetZ_S);
}

//_____________________________________________________________
void Na61VdParameters::SetupCommonVdSystem_xela75() {
  fCommonVdSystemSetup = true;
    ostringstream info;
    info << " Na61VdParameters::SetupCommonVdSystem_xela75: setting common parameters. fN=" << fN << " fRunId=" << fRunId;
    INFO(info);

  // parameters used to tune relative arms orientation and
  // describe track/hits parameters in the common VD system

   
    //  based on run 033737 and 033751 and 033755 , after rot
	fBeamSpotOffsetX = -0.202440;
	fBeamSpotSigmaX = 0.286482;
	fBeamSpotOffsetY = -0.415229;
	fBeamSpotSigmaY = 0.436122;

	fDevOffx[0] = 0.00004;    fDevSigx[0] =  0.0083260;
	fDevOffy[0] = 0.00009;    fDevSigy[0] =  0.0057305;
	fDevOffx[1] = 0.00227;    fDevSigx[1] =  0.0191340;
	fDevOffy[1] = 0.00110;    fDevSigy[1] =  0.0092248;
	fDevOffx[2] = 0.00046;    fDevSigx[2] =  0.0180239;
	fDevOffy[2] = 0.00037;    fDevSigy[2] =  0.0083328;
	fDevOffx[3] = -0.00014;    fDevSigx[3] =  0.0088461;
	fDevOffy[3] = 0.00037;    fDevSigy[3] =  0.0063926;

	fResOffx[0] = -0.00001;    fResSigx[0] =  0.0011636;
	fResOffy[0] = 0.00013;    fResSigy[0] =  0.0031588;
	fResOffx[1] = 0.00004;    fResSigx[1] =  0.0033948;
	fResOffy[1] = -0.00018;    fResSigy[1] =  0.0039178;
	fResOffx[2] = -0.00004;    fResSigx[2] =  0.0034663;
	fResOffy[2] = -0.00001;    fResSigy[2] =  0.0050602;
	fResOffx[3] = 0.00001;    fResSigx[3] =  0.0012341;
	fResOffy[3] = 0.00006;    fResSigy[3] =  0.0038171;


    //correction with eff=1
    
    fDevOffx[0] = 0.00004;    fDevSigx[0] =  0.0083256;
	fDevOffy[0] = 0.00009;    fDevSigy[0] =  0.0057298;
	fDevOffx[1] = 0.00237;    fDevSigx[1] =  0.0190177;
	fDevOffy[1] = 0.00113;    fDevSigy[1] =  0.0092003;
	fDevOffx[2] = 0.00046;    fDevSigx[2] =  0.0180261;
	fDevOffy[2] = 0.00037;    fDevSigy[2] =  0.0083339;
	fDevOffx[3] = -0.00013;    fDevSigx[3] =  0.0088683;
	fDevOffy[3] = 0.00036;    fDevSigy[3] =  0.0063934;
    

	fResOffx[0] = -0.00001;    fResSigx[0] =  0.0011638;
	fResOffy[0] = 0.00013;    fResSigy[0] =  0.0031594;
	fResOffx[1] = 0.00004;    fResSigx[1] =  0.0033949;
	fResOffy[1] = -0.00018;    fResSigy[1] =  0.0039167;
	fResOffx[2] = -0.00004;    fResSigx[2] =  0.0034665;
	fResOffy[2] = -0.00001;    fResSigy[2] =  0.0050582;
	fResOffx[3] = 0.00001;    fResSigx[3] =  0.0012343;
	fResOffy[3] = 0.00006;    fResSigy[3] =  0.0038168;

  fOffAx_J = 0.032;
  fOffAx_S = -0.0203;

  // base on run 633 (Saleve) and 635 (Jura)
  // fRotX_J =  0.001495;
  // fRotY_J =  0.0035155;
  fRotX_J = -0.001501;
  fRotY_J = -0.003418;
  // geometry after fine tunning
  fRotX_J += -0.0000521;
  fRotY_J += -0.0001430;

  // fRotX_S =  0.00165;
  // fRotY_S =  0.009315;
  fRotX_S = -0.001737;
  fRotY_S = -0.009973;
  // geometry after fine tunning
  fRotX_S += 0.0000780;
  fRotY_S += -0.0001214;


  fRotZ_J = fN * 0.001;
  fRotZ_S = -fRotZ_J;
  
  //  based on run 033737 and 033751 and 033755 
   fRotX_J = -0.001523;
   fRotY_J = -0.004489;
   fRotX_S = -0.001677;
   fRotY_S = -0.009001;

  //  based on run 033737 and 033751 and 033755
	double offsetX_J = -0.169645 + 8.;
	double offsetY_J = -0.0475559;
	double offsetZ_J = 0.28316;
	offsetX_J += 6.44966e-05;
	offsetY_J += 9.73514e-06;
	offsetZ_J += 0.00688812;

	//after rot introduced
	offsetX_J += 0.274644;
	offsetY_J += 0.00888204;
	offsetZ_J += 0.049255;
	
	offsetX_J += -1.02146e-05;
	offsetY_J += -7.27721e-06;
	offsetZ_J += 0.00144614;
	
	//from HT module
	offsetX_J += -0.00117943;
	offsetY_J += 2.11235e-05;
	offsetZ_J += 0.00321322;

	offsetX_J += -0.000117625;
	offsetY_J += -5.39421e-06;
	offsetZ_J += 0.000307625;

// new update based on run by run
  if (fRunId >= 33755) {
	offsetX_J += 0.00117181;     offsetY_J += -9.92177e-06;     offsetZ_J += -0.00316836;
  } else if (fRunId >= 33751) {
    offsetX_J += 0.00140361;     offsetY_J += 4.77639e-06;     offsetZ_J += -0.00284538;
  } else if (fRunId >= 33747) {
    offsetX_J += 0.0014527;     offsetY_J += -1.45755e-05;     offsetZ_J += -0.0034397;
  } else if (fRunId >= 33739) {
    offsetX_J += 0.0013919;     offsetY_J += -1.73949e-05;     offsetZ_J += -0.00391591;
  } else if (fRunId >= 33737) {
    offsetX_J += 0.00136333;     offsetY_J += -4.15373e-05;     offsetZ_J += -0.00409705;
  } else if (fRunId >= 33696) {
    offsetX_J += 0.00206292;     offsetY_J += 8.24748e-05;     offsetZ_J += -0.00482652;
  } else if (fRunId >= 33693) {
    offsetX_J += 0.00200226;     offsetY_J += 7.11885e-05;     offsetZ_J += -0.00440793;
  } else if (fRunId >= 33687) {
    offsetX_J += 0.00189188;     offsetY_J += 7.89166e-05;     offsetZ_J += -0.00502313;
  } else if (fRunId >= 33683) {
    offsetX_J += 0.00168935;     offsetY_J += 8.93052e-05;     offsetZ_J += -0.00474042;
  } else if (fRunId >= 33680) {
    offsetX_J += 0.00161448;     offsetY_J += 8.29366e-05;     offsetZ_J += -0.00463108;
  } else if (fRunId >= 33646) {
    offsetX_J += 0.00173885;     offsetY_J += 9.90801e-05;     offsetZ_J += -0.00507764;
  } else if (fRunId >= 33635) {
    offsetX_J += 0.00201384;     offsetY_J += 0.000132599;     offsetZ_J += -0.00533999;
  } else if (fRunId >= 33633) {
    offsetX_J += 0.00168163;     offsetY_J += 7.21731e-05;     offsetZ_J += -0.00572286;
  } else if (fRunId >= 33613) {
    offsetX_J += 0.00191602;     offsetY_J += 0.000155017;     offsetZ_J += -0.00592936;
  } else if (fRunId >= 33585) {
    offsetX_J += 0.00162651;     offsetY_J += 7.40647e-05;     offsetZ_J += -0.00559692;
  } else if (fRunId >= 33549) {
    offsetX_J += 0.00144782;     offsetY_J += 9.52852e-05;     offsetZ_J += -0.00497486;
  } else if (fRunId >= 33539) {
    offsetX_J += 0.00140322;     offsetY_J += 8.02691e-05;     offsetZ_J += -0.00495692;
  } else if (fRunId >= 33535) {
    offsetX_J += 0.00136263;     offsetY_J += 2.23327e-05;     offsetZ_J += -0.00516033;
  } else if (fRunId >= 33500) {
    offsetX_J += 0.00135606;     offsetY_J += -7.58689e-05;     offsetZ_J += -0.00494307;
  } else if (fRunId >= 33498) {
    offsetX_J += 0.00140737;     offsetY_J += -5.13983e-05;     offsetZ_J += -0.00487665;
  } else if (fRunId >= 33495) {
    offsetX_J += 0.00124445;     offsetY_J += -6.65723e-05;     offsetZ_J += -0.00586401;
  } else if (fRunId >= 33493) {
    offsetX_J += 0.00139003;     offsetY_J += -0.000117258;     offsetZ_J += -0.00536641;
  } else if (fRunId >= 33490) {
    offsetX_J += 0.00131004;     offsetY_J += -0.000105524;     offsetZ_J += -0.00507447;
  } else if (fRunId >= 33484) {
    offsetX_J += 0.00132318;     offsetY_J += -8.98758e-05;     offsetZ_J += -0.00537202;
  } else if (fRunId >= 33481) {
    offsetX_J += 0.00134953;     offsetY_J += -8.11044e-05;     offsetZ_J += -0.0053747;
  } else if (fRunId >= 33480) {
    offsetX_J += 0.00147363;     offsetY_J += -5.36121e-05;     offsetZ_J += -0.00546757;
  } else if (fRunId >= 33469) {
    offsetX_J += 0.00154607;     offsetY_J += -4.77977e-06;     offsetZ_J += -0.00532963;
  } else if (fRunId >= 33464) {
    offsetX_J += 0.00271196;     offsetY_J += 7.95315e-05;     offsetZ_J += -0.00533353;
  } else if (fRunId >= 33461) {
    offsetX_J += 0.00185649;     offsetY_J += 8.64054e-05;     offsetZ_J += -0.00532103;
  } else if (fRunId >= 33450) {
    offsetX_J += 0.00192931;     offsetY_J += 0.000254429;     offsetZ_J += -0.000162216; //value for 33453
  }
  
  // new update based on run by run after HT
  
  if (fRunId >= 33755) {
	offsetX_J += -0.00122845;     offsetY_J += 3.1864e-05;     offsetZ_J += 0.00260078;
  } else if (fRunId >= 33751) {
    offsetX_J += -0.001247;     offsetY_J += 5.04238e-06;     offsetZ_J += 0.00284398;
  } else if (fRunId >= 33747) {
    offsetX_J += -0.00114703;     offsetY_J += 1.37787e-05;     offsetZ_J += 0.00330009;
  } else if (fRunId >= 33739) {
    offsetX_J += -0.00111924;     offsetY_J += 1.48819e-05;     offsetZ_J += 0.00349371;
  } else if (fRunId >= 33737) {
    offsetX_J += -0.00113313;     offsetY_J += 5.63851e-07;     offsetZ_J += 0.00378839;
  } else if (fRunId >= 33696) {
    offsetX_J += -0.000960401;     offsetY_J += 2.61204e-05;     offsetZ_J += 0.00420493;
  } else if (fRunId >= 33693) {
    offsetX_J += -0.00105223;     offsetY_J += 6.60496e-06;     offsetZ_J += 0.00415894;
  } else if (fRunId >= 33687) {
    offsetX_J += -0.00103605;     offsetY_J += 3.65455e-06;     offsetZ_J += 0.00386825;
  } else if (fRunId >= 33683) {
    offsetX_J += -0.000933797;     offsetY_J += 6.33293e-06;     offsetZ_J += 0.00420108;
  } else if (fRunId >= 33680) {
    offsetX_J += -0.000923656;     offsetY_J += -1.11628e-05;     offsetZ_J += 0.00388327;
  } else if (fRunId >= 33646) {
    offsetX_J += -0.0011083;     offsetY_J += 9.11805e-06;     offsetZ_J += 0.00400038;
  } else if (fRunId >= 33635) {
    offsetX_J += -0.0010967;     offsetY_J += 1.29806e-05;     offsetZ_J += 0.00392684;
  } else if (fRunId >= 33633) {
    offsetX_J += -0.00110609;     offsetY_J += 2.83364e-05;     offsetZ_J += 0.0041388;
  } else if (fRunId >= 33613) {
    offsetX_J += -0.00107129;     offsetY_J += 1.35447e-05;     offsetZ_J += 0.00427384;
  } else if (fRunId >= 33585) {
    offsetX_J += -0.000989886;     offsetY_J += 1.72852e-05;     offsetZ_J += 0.00417737;
  } else if (fRunId >= 33549) {
    offsetX_J += -0.00104797;     offsetY_J += -5.44767e-06;     offsetZ_J += 0.00317703;
  } else if (fRunId >= 33539) {
    offsetX_J += -0.00104676;     offsetY_J += -2.59901e-05;     offsetZ_J += 0.00392089;
  } else if (fRunId >= 33535) {
    offsetX_J += -0.000981947;     offsetY_J += 2.00819e-05;     offsetZ_J += 0.00397585;
  } else if (fRunId >= 33500) {
    offsetX_J += -0.000987228;     offsetY_J += -9.8872e-06;     offsetZ_J += 0.00339007;
  } else if (fRunId >= 33498) {
    offsetX_J += -0.000999219;     offsetY_J += -1.08739e-05;     offsetZ_J += 0.00336244;
  } else if (fRunId >= 33495) {
    offsetX_J += -0.000893575;     offsetY_J += -2.76092e-06;     offsetZ_J += 0.00413896;
  } else if (fRunId >= 33493) {
    offsetX_J += -0.000982277;     offsetY_J += 2.00696e-05;     offsetZ_J += 0.00375782;
  } else if (fRunId >= 33490) {
    offsetX_J += -0.000908882;     offsetY_J += 5.52695e-06;     offsetZ_J += 0.00345703;
  } else if (fRunId >= 33484) {
    offsetX_J += -0.000895451;     offsetY_J += -1.59054e-06;     offsetZ_J += 0.00398583;
  } else if (fRunId >= 33481) {
    offsetX_J += -0.00087062;     offsetY_J += 3.36543e-07;     offsetZ_J += 0.00410817;
  } else if (fRunId >= 33480) {
    offsetX_J += -0.000870755;     offsetY_J += 2.3949e-06;     offsetZ_J += 0.00425633;
  } else if (fRunId >= 33469) {
    offsetX_J += -0.000824008;     offsetY_J += -7.8319e-06;     offsetZ_J += 0.00379544;
  } else if (fRunId >= 33464) {
    offsetX_J += -0.000958153;     offsetY_J += 3.33117e-06;     offsetZ_J += 0.00455568;
  } else if (fRunId >= 33461) {
    offsetX_J += -0.000940764;     offsetY_J += 6.40601e-06;     offsetZ_J += 0.00402915;
  } //else if (fRunId >= 33453) {
    //offsetX_J += -0.657463;     offsetY_J += 0.000259892;     offsetZ_J += -0.0115008;//?
  //}

	
	double offsetX_S =-offsetX_J;
	double offsetY_S =-offsetY_J;
	double offsetZ_S =-offsetZ_J;

 
	ostringstream info_fin;
    info_fin << "Offset values: Jura: "<<offsetX_J<<" "<<offsetY_J<<" "<<offsetZ_J<<"   Saleve: "<<offsetX_S<<" "<<offsetY_S<<" "<<offsetZ_S;
    INFO(info_fin);

  // parameter to extract Jura and Saleve events with consistent vertexes.
  fVtxDzSigma = 0.075;

  fJuraArmOffset.SetX(offsetX_J);
  fJuraArmOffset.SetY(offsetY_J);
  fJuraArmOffset.SetZ(offsetZ_J);
  fSaleveArmOffset.SetX(offsetX_S);
  fSaleveArmOffset.SetY(offsetY_S);
  fSaleveArmOffset.SetZ(offsetZ_S);
}


//_____________________________________________________________
void Na61VdParameters::SetupCommonVdSystem_xela40() {
  fCommonVdSystemSetup = true;
    ostringstream info;
    info << " Na61VdParameters::SetupCommonVdSystem_xela40: setting common parameters. fN=" << fN << " fRunId=" << fRunId;
    INFO(info);

  // parameters used to tune relative arms orientation and
  // describe track/hits parameters in the common VD system

   
    //  based on run 35206, 35229, 35264 , after rot
	fBeamSpotOffsetX = 0.133809;
	fBeamSpotSigmaX = 0.603616;
	fBeamSpotOffsetY = -0.219748;
	fBeamSpotSigmaY = 0.469055;


	fDevOffx[0] = 0.00060;    fDevSigx[0] =  0.0072233;
	fDevOffy[0] = 0.00052;    fDevSigy[0] =  0.0061024;
	fDevOffx[1] = 0.00119;    fDevSigx[1] =  0.0124785;
	fDevOffy[1] = 0.00036;    fDevSigy[1] =  0.0091453;
	fDevOffx[2] = -0.00234;    fDevSigx[2] =  0.0116636;
	fDevOffy[2] = 0.00071;    fDevSigy[2] =  0.0089294;
	fDevOffx[3] = -0.00202;    fDevSigx[3] =  0.0075391;
	fDevOffy[3] = 0.00080;    fDevSigy[3] =  0.0064536;


	fResOffx[0] = -0.00016;    fResSigx[0] =  0.0012208;
	fResOffy[0] = 0.00032;    fResSigy[0] =  0.0033416;
	fResOffx[1] = 0.00050;    fResSigx[1] =  0.0035955;
	fResOffy[1] = -0.00037;    fResSigy[1] =  0.0040860;
	fResOffx[2] = -0.00052;    fResSigx[2] =  0.0037244;
	fResOffy[2] = -0.00020;    fResSigy[2] =  0.0053002;
	fResOffx[3] = 0.00018;    fResSigx[3] =  0.0013491;
	fResOffy[3] = 0.00026;    fResSigy[3] =  0.0040594;



  fOffAx_J = 0.032;
  fOffAx_S = -0.0203;

  // base on run 633 (Saleve) and 635 (Jura)
  // fRotX_J =  0.001495;
  // fRotY_J =  0.0035155;
  fRotX_J = -0.001501;
  fRotY_J = -0.003418;
  // geometry after fine tunning
  fRotX_J += -0.0000521;
  fRotY_J += -0.0001430;

  // fRotX_S =  0.00165;
  // fRotY_S =  0.009315;
  fRotX_S = -0.001737;
  fRotY_S = -0.009973;
  // geometry after fine tunning
  fRotX_S += 0.0000780;
  fRotY_S += -0.0001214;


  fRotZ_J = fN * 0.001;
  fRotZ_S = -fRotZ_J;
  
  //  based from XeLa75
   fRotX_J = -0.001523;
   fRotY_J = -0.004489;
   fRotX_S = -0.001677;
   fRotY_S = -0.009001;
   
   fRotX_J += 0.000114;
   fRotY_J += -0.000266;
   fRotX_S += 0.000121;
   fRotY_S += 0.000212;
   
  //  based on run 033737 and 033751 and 033755
	double offsetX_J =1.37433 + 8.;
	double offsetY_J =-0.0365227;
	double offsetZ_J =0.340986;
	
	offsetX_J += -0.000243212;
	offsetY_J += -8.19897e-05;
	offsetZ_J += 0.0072696;
	
	offsetX_J += -0.0292401;
	offsetY_J += -0.000382558;
	offsetZ_J += 0.000154213;
	
	//HT
	offsetX_J += -0.00127657;
	offsetY_J += -0.000183839;
	offsetZ_J += 0.0066011;
	
	offsetX_J += -0.000157147;
	offsetY_J += -2.54134e-05;
	offsetZ_J += 0.000940774;

// new update based on run by run


	
	double offsetX_S =-offsetX_J;
	double offsetY_S =-offsetY_J;
	double offsetZ_S =-offsetZ_J;

 
	ostringstream info_fin;
    info_fin << "Offset values: Jura: "<<offsetX_J<<" "<<offsetY_J<<" "<<offsetZ_J<<"   Saleve: "<<offsetX_S<<" "<<offsetY_S<<" "<<offsetZ_S;
    INFO(info_fin);

  // parameter to extract Jura and Saleve events with consistent vertexes.
  fVtxDzSigma = 0.075;

  fJuraArmOffset.SetX(offsetX_J);
  fJuraArmOffset.SetY(offsetY_J);
  fJuraArmOffset.SetZ(offsetZ_J);
  fSaleveArmOffset.SetX(offsetX_S);
  fSaleveArmOffset.SetY(offsetY_S);
  fSaleveArmOffset.SetZ(offsetZ_S);
}


//_____________________________________________________________
void Na61VdParameters::SetupCommonVdSystem_pPb2018() {
  fCommonVdSystemSetup = true;
    ostringstream info;
    info << " Na61VdParameters::SetupCommonVdSystem_pPb2018: setting common parameters. fN=" << fN << " fRunId=" << fRunId;
    INFO(info);

  // parameters used to tune relative arms orientation and
  // describe track/hits parameters in the common VD system

  // production(field) defauls values
  fBeamSpotOffsetX = -0.526628;
  fBeamSpotSigmaX = 0.592206;
  fBeamSpotOffsetY = -0.616405;
  fBeamSpotSigmaY = 0.802218;

  fDevOffx[0] = -0.00020;
  fDevSigx[0] = 0.0143933;
  fDevOffy[0] = 0.00014;
  fDevSigy[0] = 0.0058761;
  fDevOffx[1] = 0.00377;
  fDevSigx[1] = 0.0353684;
  fDevOffy[1] = 0.00119;
  fDevSigy[1] = 0.0097333;
  fDevOffx[2] = 0.00033;
  fDevSigx[2] = 0.0396898;
  fDevOffy[2] = 0.00009;
  fDevSigy[2] = 0.0083968;
  fDevOffx[3] = -0.00004;
  fDevSigx[3] = 0.0175264;
  fDevOffy[3] = 0.00057;
  fDevSigy[3] = 0.0072806;

  /*
  fResOffx[0] = 0.00001;    fResSigx[0] =  0.0011874;
  fResOffy[0] = -0.00001;    fResSigy[0] =  0.0050000; // limits affected the fits
  fResOffx[1] = -0.00003;    fResSigx[1] =  0.0035086;
  fResOffy[1] = -0.00004;    fResSigy[1] =  0.0050000;
  fResOffx[2] = 0.00003;    fResSigx[2] =  0.0035965;
  fResOffy[2] = 0.00005;    fResSigy[2] =  0.0052853;
  fResOffx[3] = -0.00001;    fResSigx[3] =  0.0012737;
  fResOffy[3] = -0.00002;    fResSigy[3] =  0.0050000;
  */

  fResOffx[0] = 0.00001;
  fResSigx[0] = 0.0012111;
  fResOffy[0] = 0.00002;
  fResSigy[0] = 0.0034618;
  fResOffx[1] = -0.00003;
  fResSigx[1] = 0.0035005;
  fResOffy[1] = -0.00006;
  fResSigy[1] = 0.0041643;
  fResOffx[2] = 0.00003;
  fResSigx[2] = 0.0035930;
  fResOffy[2] = 0.00002;
  fResSigy[2] = 0.0052658;
  fResOffx[3] = -0.00001;
  fResSigx[3] = 0.0012740;
  fResOffy[3] = -0.00000;
  fResSigy[3] = 0.0040724;

  //  if(fRunId>=1350){
  //  if(fRunId>=33201){
  if (fRunId >= 31974) {
    fBeamSpotOffsetX = -0.405983;
    fBeamSpotSigmaX = 0.577560;
    fBeamSpotOffsetY = -0.615915;
    fBeamSpotSigmaY = 0.818117;
  }

  //  if(fRunId>=1494){ // 0 field runs
  if (fRunId >= 33396) {  // 0 field runs
    fBeamSpotOffsetX = 0.459624;
    fBeamSpotSigmaX = 0.565917;
    fBeamSpotOffsetY = 0.540213;
    fBeamSpotSigmaY = 0.885769;
  }

  fOffAx_J = 0.032;
  fOffAx_S = -0.0203;

  // base on run 633 (Saleve) and 635 (Jura)
  // fRotX_J =  0.001495;
  // fRotY_J =  0.0035155;
  fRotX_J = -0.001501;
  fRotY_J = -0.003418;
  // geometry after fine tunning
  fRotX_J += -0.0000521;
  fRotY_J += -0.0001430;

  // fRotX_S =  0.00165;
  // fRotY_S =  0.009315;

  fRotX_S = -0.001737;
  fRotY_S = -0.009973;
  // geometry after fine tunning
  fRotX_S += 0.0000780;
  fRotY_S += -0.0001214;

  fRotZ_J = fN * 0.001;
  fRotZ_S = -fRotZ_J;

  // no arm rotations use it for test (needs to adjust arm offsets)
  // fRotX_J =  0.0;
  // fRotY_J =  0.0;

  // fRotX_S =  0.0;
  // fRotY_S =  0.0;

  // correction based on run 1494 (0 field)
  double offsetX_J = 0.881427 + 8.;
  double offsetY_J = -0.0164678;
  double offsetZ_J = 0.317669;

  double offsetX_S = -0.881427 - 8.;
  double offsetY_S = 0.0164678;
  double offsetZ_S = -0.317669;

  offsetX_J += -0.00169475;
  offsetY_J += -0.0244815;
  offsetZ_J += 0.0332515;

  offsetX_S += 0.00169475;
  offsetY_S += 0.0244815;
  offsetZ_S += -0.0332515;

  // form VdTrackingHT module
  offsetX_J += 0.000539333;
  offsetY_J += -0.000280796;
  offsetZ_J += 0.00251891;

  offsetX_S += -0.000539333;
  offsetY_S += 0.000280796;
  offsetZ_S += -0.00251891;

  // parameter to extract Jura and Saleve events with consistent vertexes.
  fVtxDzSigma = 0.075;

	ostringstream info_fin;
    info_fin << "Offset values: Jura: "<<offsetX_J<<" "<<offsetY_J<<" "<<offsetZ_J<<"   Saleve: "<<offsetX_S<<" "<<offsetY_S<<" "<<offsetZ_S;
    INFO(info_fin);

  fJuraArmOffset.SetX(offsetX_J);
  fJuraArmOffset.SetY(offsetY_J);
  fJuraArmOffset.SetZ(offsetZ_J);
  fSaleveArmOffset.SetX(offsetX_S);
  fSaleveArmOffset.SetY(offsetY_S);
  fSaleveArmOffset.SetZ(offsetZ_S);

}

//_____________________________________________________________
void Na61VdParameters::SetupCommonVdSystem_pPb2022() {
  fCommonVdSystemSetup = true;
    ostringstream info;
    info << " Na61VdParameters::SetupCommonVdSystem_pPb2022: setting common parameters. fN=" << fN << " fRunId=" << fRunId;
    INFO(info);

  // parameters used to tune relative arms orientation and
  // describe track/hits parameters in the common VD system

  // production(field) defauls values
  fBeamSpotOffsetX = 0.;
  fBeamSpotSigmaX = 1.5;
  fBeamSpotOffsetY = 0.;
  fBeamSpotSigmaY = 1.5;

  fDevOffx[0] = -0.00020;
  fDevSigx[0] = 0.0143933;
  fDevOffy[0] = 0.00014;
  fDevSigy[0] = 0.0058761;
  fDevOffx[1] = 0.00377;
  fDevSigx[1] = 0.0353684;
  fDevOffy[1] = 0.00119;
  fDevSigy[1] = 0.0097333;
  fDevOffx[2] = 0.00033;
  fDevSigx[2] = 0.0396898;
  fDevOffy[2] = 0.00009;
  fDevSigy[2] = 0.0083968;
  fDevOffx[3] = -0.00004;
  fDevSigx[3] = 0.0175264;
  fDevOffy[3] = 0.00057;
  fDevSigy[3] = 0.0072806;


  fResOffx[0] = 0.00001;
  fResSigx[0] = 0.0012111;
  fResOffy[0] = 0.00002;
  fResSigy[0] = 0.0034618;
  fResOffx[1] = -0.00003;
  fResSigx[1] = 0.0035005;
  fResOffy[1] = -0.00006;
  fResSigy[1] = 0.0041643;
  fResOffx[2] = 0.00003;
  fResSigx[2] = 0.0035930;
  fResOffy[2] = 0.00002;
  fResSigy[2] = 0.0052658;
  fResOffx[3] = -0.00001;
  fResSigx[3] = 0.0012740;
  fResOffy[3] = -0.00000;
  fResSigy[3] = 0.0040724;

  //  if(fRunId>=1350){
  //  if(fRunId>=33201){
 
 
  fOffAx_J = 0.0;
  fOffAx_S = 0.0;

  // base on run 633 (Saleve) and 635 (Jura)
  // fRotX_J =  0.001495;
  // fRotY_J =  0.0035155;
  fRotX_J = -0.0;
  fRotY_J = -0.0;
  fRotZ_J = -0.0;

  fRotX_S = -0.0;
  fRotY_S = -0.0;
  fRotZ_S = -0.0;

  // no arm rotations use it for test (needs to adjust arm offsets)
  // fRotX_J =  0.0;
  // fRotY_J =  0.0;

  // fRotX_S =  0.0;
  // fRotY_S =  0.0;

  // correction based on run 1494 (0 field)
  double offsetX_J =  12.5 + 2.245/2.;  // 5 mm gap + 7.5 mm 1/2 sensor x length
  double offsetY_J = -0.0;
  double offsetZ_J = 0.0;
  
  double offsetX_S =-offsetX_J;
  double offsetY_S =-offsetY_J;
  double offsetZ_S =-offsetZ_J;
  

  // parameter to extract Jura and Saleve events with consistent vertexes.
  fVtxDzSigma = 0.075;
  
  ostringstream info_fin;
  info_fin << "Offset values: Jura: "<<offsetX_J<<" "<<offsetY_J<<" "<<offsetZ_J<<"   Saleve: "<<offsetX_S<<" "<<offsetY_S<<" "<<offsetZ_S;
  INFO(info_fin);
  
  fJuraArmOffset.SetX(offsetX_J);
  fJuraArmOffset.SetY(offsetY_J);
  fJuraArmOffset.SetZ(offsetZ_J);
  fSaleveArmOffset.SetX(offsetX_S);
  fSaleveArmOffset.SetY(offsetY_S);
  fSaleveArmOffset.SetZ(offsetZ_S);

}


//_____________________________________________________________
void Na61VdParameters::SetupCommonVdSystem_pbpb() {
  fCommonVdSystemSetup = true;
    ostringstream info;
    info << " Na61VdParameters::SetupCommonVdSystem_pbpb: setting common parameters. fN=" << fN << " fRunId=" << fRunId;
    INFO(info);

  // parameters used to tune relative arms orientation and
  // describe track/hits parameters in the common VD system

  // production(field) defauls values
  

  // PbPb 2018
  // based on runs 38943, 38945
  if (fRunId >= 38014) { 
    
	fBeamSpotOffsetX = -0.315266;
	fBeamSpotSigmaX = 0.390449;
	fBeamSpotOffsetY = -0.393240;
	fBeamSpotSigmaY = 0.271049;
			 if (fRunId >= 38945){	    fBeamSpotOffsetX = -0.314910;    fBeamSpotSigmaX = 0.397528;    fBeamSpotOffsetY = -0.390983;    fBeamSpotSigmaY = 0.272693;
	  } else if (fRunId >= 38943){ 	    fBeamSpotOffsetX = -0.315923;    fBeamSpotSigmaX = 0.400060;    fBeamSpotOffsetY = -0.393236;    fBeamSpotSigmaY = 0.270938;
	  } else if (fRunId >= 38941){ 	    fBeamSpotOffsetX = -0.307692;    fBeamSpotSigmaX = 0.417349;    fBeamSpotOffsetY = -0.384857;    fBeamSpotSigmaY = 0.263262;
	  } else if (fRunId >= 38938){ 	    fBeamSpotOffsetX = -0.300248;    fBeamSpotSigmaX = 0.415642;    fBeamSpotOffsetY = -0.352272;    fBeamSpotSigmaY = 0.269774;	  
	  } else if (fRunId >= 38919){ 	    fBeamSpotOffsetX = -0.022644;    fBeamSpotSigmaX = 0.408605;    fBeamSpotOffsetY = -0.274234;    fBeamSpotSigmaY = 0.266813;
	  } else if (fRunId >= 38917){ 	    fBeamSpotOffsetX = -0.036510;    fBeamSpotSigmaX = 0.419535;    fBeamSpotOffsetY = -0.143428;    fBeamSpotSigmaY = 0.248807;
	  } else if (fRunId >= 38898){ 	    fBeamSpotOffsetX = 0.260028;    fBeamSpotSigmaX = 0.610152;    fBeamSpotOffsetY = 0.125676;    fBeamSpotSigmaY = 0.275824;
	  } else if (fRunId >= 38889){ 	    fBeamSpotOffsetX = 0.050160;    fBeamSpotSigmaX = 0.560117;    fBeamSpotOffsetY = 0.004954;    fBeamSpotSigmaY = 0.269170;
	  } else if (fRunId >= 38881){ 	    fBeamSpotOffsetX = -0.020753;    fBeamSpotSigmaX = 0.395114;    fBeamSpotOffsetY = -0.055704;    fBeamSpotSigmaY = 0.257536;
	  } else if (fRunId >= 38835){ 	    fBeamSpotOffsetX = 0.158644;    fBeamSpotSigmaX = 0.365269;    fBeamSpotOffsetY = 0.123144;    fBeamSpotSigmaY = 0.255019;
	  } else if (fRunId >= 38826){ 	    fBeamSpotOffsetX = 0.188219;    fBeamSpotSigmaX = 0.381117;    fBeamSpotOffsetY = 0.097576;    fBeamSpotSigmaY = 0.251987;
	  } else if (fRunId >= 38825){ 	    fBeamSpotOffsetX = 0.200993;    fBeamSpotSigmaX = 0.369313;    fBeamSpotOffsetY = 0.108543;    fBeamSpotSigmaY = 0.248488;
	  } else if (fRunId >= 38813){ 	    fBeamSpotOffsetX = -0.110521;    fBeamSpotSigmaX = 0.362012;    fBeamSpotOffsetY = -0.291285;    fBeamSpotSigmaY = 0.237577;
	  } else if (fRunId >= 38812){ 	    fBeamSpotOffsetX = -0.101121;    fBeamSpotSigmaX = 0.373023;    fBeamSpotOffsetY = -0.292899;    fBeamSpotSigmaY = 0.240757;
	  } else if (fRunId >= 38804){ 	    fBeamSpotOffsetX = 0.203470;    fBeamSpotSigmaX = 0.360918;    fBeamSpotOffsetY = -0.281242;    fBeamSpotSigmaY = 0.232094;
	  } else if (fRunId >= 38796){ 	    fBeamSpotOffsetX = 0.184961;    fBeamSpotSigmaX = 0.364276;    fBeamSpotOffsetY = -0.321358;    fBeamSpotSigmaY = 0.232780;
	  } else if (fRunId >= 38775){ 	    fBeamSpotOffsetX = 0.031285;    fBeamSpotSigmaX = 0.368318;    fBeamSpotOffsetY = -0.202846;    fBeamSpotSigmaY = 0.242144;
	  } else if (fRunId >= 38772){ 	    fBeamSpotOffsetX = 0.146861;    fBeamSpotSigmaX = 0.375766;    fBeamSpotOffsetY = -0.233925;    fBeamSpotSigmaY = 0.250192;
	  } else if (fRunId >= 38771){ 	    fBeamSpotOffsetX = 0.190212;    fBeamSpotSigmaX = 0.383935;    fBeamSpotOffsetY = -0.237299;    fBeamSpotSigmaY = 0.239368;
	  } else if (fRunId >= 38764){ 	    fBeamSpotOffsetX = 0.191144;    fBeamSpotSigmaX = 0.374213;    fBeamSpotOffsetY = -0.255284;    fBeamSpotSigmaY = 0.233877;
	  } else if (fRunId >= 38761){ 	    fBeamSpotOffsetX = 0.119230;    fBeamSpotSigmaX = 0.371152;    fBeamSpotOffsetY = -0.156990;    fBeamSpotSigmaY = 0.244888;
	  } else if (fRunId >= 38759){ 	    fBeamSpotOffsetX = 0.110854;    fBeamSpotSigmaX = 0.376810;    fBeamSpotOffsetY = -0.227913;    fBeamSpotSigmaY = 0.239612;
	  } else if (fRunId >= 38739){ 	    fBeamSpotOffsetX = 0.368412;    fBeamSpotSigmaX = 0.381263;    fBeamSpotOffsetY = -0.288948;    fBeamSpotSigmaY = 0.241720;
	  } else if (fRunId >= 38707){ 	    fBeamSpotOffsetX = -0.003201;    fBeamSpotSigmaX = 0.343495;    fBeamSpotOffsetY = -0.758776;    fBeamSpotSigmaY = 0.418685;
	  } else if (fRunId >= 38667){ 	    fBeamSpotOffsetX = -0.216552;    fBeamSpotSigmaX = 0.338433;    fBeamSpotOffsetY = -0.182371;    fBeamSpotSigmaY = 0.378989;
	  } else if (fRunId >= 38662){ 	    fBeamSpotOffsetX = -0.183302;    fBeamSpotSigmaX = 0.313603;    fBeamSpotOffsetY = -0.130294;    fBeamSpotSigmaY = 0.330914;
	  } else if (fRunId >= 38657){ 	    fBeamSpotOffsetX = -0.190382;    fBeamSpotSigmaX = 0.316115;    fBeamSpotOffsetY = -0.093464;    fBeamSpotSigmaY = 0.305925;
	  } else if (fRunId >= 38656){ 	    fBeamSpotOffsetX = -0.190831;    fBeamSpotSigmaX = 0.341717;    fBeamSpotOffsetY = -0.101467;    fBeamSpotSigmaY = 0.310016;
	  } else if (fRunId >= 38654){ 	    fBeamSpotOffsetX = -0.172658;    fBeamSpotSigmaX = 0.343398;    fBeamSpotOffsetY = -0.106690;    fBeamSpotSigmaY = 0.313704;
	  } else if (fRunId >= 38649){ 	    fBeamSpotOffsetX = -0.168801;    fBeamSpotSigmaX = 0.337019;    fBeamSpotOffsetY = -0.102245;    fBeamSpotSigmaY = 0.298868;
	  } else if (fRunId >= 38645){ 	    fBeamSpotOffsetX = -0.184935;    fBeamSpotSigmaX = 0.359570;    fBeamSpotOffsetY = -0.162048;    fBeamSpotSigmaY = 0.327461;
	  } else if (fRunId >= 38636){ 	    fBeamSpotOffsetX = -0.136874;    fBeamSpotSigmaX = 0.356790;    fBeamSpotOffsetY = -0.168621;    fBeamSpotSigmaY = 0.320773;
	  } else if (fRunId >= 38628){ 	    fBeamSpotOffsetX = -0.136132;    fBeamSpotSigmaX = 0.368500;    fBeamSpotOffsetY = -0.218288;    fBeamSpotSigmaY = 0.360055;	  
	  } else if (fRunId >= 38619){ 	    fBeamSpotOffsetX = -0.093325;    fBeamSpotSigmaX = 0.356695;    fBeamSpotOffsetY = -0.262146;    fBeamSpotSigmaY = 0.356815;
	  } else if (fRunId >= 38606){ 	    fBeamSpotOffsetX = -0.129541;    fBeamSpotSigmaX = 0.357516;    fBeamSpotOffsetY = -0.322125;    fBeamSpotSigmaY = 0.380516;
	  } else if (fRunId >= 38599){ 	    fBeamSpotOffsetX = -0.130479;    fBeamSpotSigmaX = 0.348624;    fBeamSpotOffsetY = -0.348972;    fBeamSpotSigmaY = 0.390033;
	  } else if (fRunId >= 38593){ 	    fBeamSpotOffsetX = -0.120376;    fBeamSpotSigmaX = 0.359736;    fBeamSpotOffsetY = -0.329073;    fBeamSpotSigmaY = 0.380767;
	  } else if (fRunId >= 38568){ 	    fBeamSpotOffsetX = -0.199220;    fBeamSpotSigmaX = 0.358996;    fBeamSpotOffsetY = -0.420983;    fBeamSpotSigmaY = 0.413606;
	  } else if (fRunId >= 38562){ 	    fBeamSpotOffsetX = -0.170591;    fBeamSpotSigmaX = 0.349401;    fBeamSpotOffsetY = -0.358890;    fBeamSpotSigmaY = 0.358633;
	  } else if (fRunId >= 38548){ 	    fBeamSpotOffsetX = -0.082266;    fBeamSpotSigmaX = 0.371108;    fBeamSpotOffsetY = -0.240715;    fBeamSpotSigmaY = 0.285424;
	  } else if (fRunId >= 38450){ 	    fBeamSpotOffsetX = -0.103644;    fBeamSpotSigmaX = 0.355954;    fBeamSpotOffsetY = -0.465219;    fBeamSpotSigmaY = 0.326009;
	  } else if (fRunId >= 38436){ 	    fBeamSpotOffsetX = -0.067752;    fBeamSpotSigmaX = 0.351429;    fBeamSpotOffsetY = -0.527018;    fBeamSpotSigmaY = 0.310448;
	  //} else if (fRunId >= 38444){ 	    fBeamSpotOffsetX = -0.066820;    fBeamSpotSigmaX = 0.350380;    fBeamSpotOffsetY = -0.537822;    fBeamSpotSigmaY = 0.300548;
	  }
	/*  } else if (fRunId >= 38300){ 	    fBeamSpotOffsetX = -0.276055;    fBeamSpotSigmaX = 0.332354;    fBeamSpotOffsetY = -0.285572;    fBeamSpotSigmaY = 0.280116;
	  } else if (fRunId >= 38295){ 	    fBeamSpotOffsetX = -0.295789;    fBeamSpotSigmaX = 0.350060;    fBeamSpotOffsetY = -0.311217;    fBeamSpotSigmaY = 0.266510;
	  } else if (fRunId >= 38291){ 	    fBeamSpotOffsetX = -0.297017;    fBeamSpotSigmaX = 0.337803;    fBeamSpotOffsetY = -0.341474;    fBeamSpotSigmaY = 0.268942;
	  } else if (fRunId >= 38290){ 	    fBeamSpotOffsetX = -0.320238;    fBeamSpotSigmaX = 0.341924;    fBeamSpotOffsetY = -0.341856;    fBeamSpotSigmaY = 0.272757;
	  } else if (fRunId >= 38283){ 	    fBeamSpotOffsetX = -0.277983;    fBeamSpotSigmaX = 0.341564;    fBeamSpotOffsetY = -0.302643;    fBeamSpotSigmaY = 0.258874;
	  } else if (fRunId >= 38278){ 	    fBeamSpotOffsetX = -0.264565;    fBeamSpotSigmaX = 0.361360;    fBeamSpotOffsetY = -0.329286;    fBeamSpotSigmaY = 0.269027;
	  } else if (fRunId >= 38275){ 	    fBeamSpotOffsetX = -0.255039;    fBeamSpotSigmaX = 0.345569;    fBeamSpotOffsetY = -0.402907;    fBeamSpotSigmaY = 0.247899;
	  } else if (fRunId >= 38242){ 	    fBeamSpotOffsetX = -0.277674;    fBeamSpotSigmaX = 0.341513;    fBeamSpotOffsetY = -0.446976;    fBeamSpotSigmaY = 0.257949;
	  } else if (fRunId >= 38240){ 	    fBeamSpotOffsetX = -0.266908;    fBeamSpotSigmaX = 0.344124;    fBeamSpotOffsetY = -0.457608;    fBeamSpotSigmaY = 0.257980;
	  } else if (fRunId >= 38238){ 	    fBeamSpotOffsetX = -0.257225;    fBeamSpotSigmaX = 0.332304;    fBeamSpotOffsetY = -0.469844;    fBeamSpotSigmaY = 0.249562;
	  } else if (fRunId >= 38233){      fBeamSpotOffsetX = -0.278193;    fBeamSpotSigmaX = 0.349050;    fBeamSpotOffsetY = -0.463846;    fBeamSpotSigmaY = 0.246232;
	  } else if (fRunId >= 38232){ 	    fBeamSpotOffsetX = -0.277911;    fBeamSpotSigmaX = 0.354994;    fBeamSpotOffsetY = -0.463903;    fBeamSpotSigmaY = 0.244595;
	  }	  
	*/	
  }
  
  
	fDevOffx[0] = 0.00045;    fDevSigx[0] =  0.0122845;
	fDevOffy[0] = -0.00014;    fDevSigy[0] =  0.0055300;
	fDevOffx[1] = 0.00015;    fDevSigx[1] =  0.0312570;
	fDevOffy[1] = 0.00051;    fDevSigy[1] =  0.0083060;
	fDevOffx[2] = 0.00023;    fDevSigx[2] =  0.0348309;
	fDevOffy[2] = -0.00008;    fDevSigy[2] =  0.0077215;
	fDevOffx[3] = -0.00125;    fDevSigx[3] =  0.0129276;
	fDevOffy[3] = 0.00025;    fDevSigy[3] =  0.0069210;



	fResOffx[0] = -0.00016;    fResSigx[0] =  0.0011881;
	fResOffy[0] = -0.00002;    fResSigy[0] =  0.0031819;
	fResOffx[1] = 0.00050;    fResSigx[1] =  0.0034572;
	fResOffy[1] = -0.00002;    fResSigy[1] =  0.0038821;
	fResOffx[2] = -0.00051;    fResSigx[2] =  0.0035257;
	fResOffy[2] = 0.00012;    fResSigy[2] =  0.0053626;
	fResOffx[3] = 0.00017;    fResSigx[3] =  0.0012564;
	fResOffy[3] = -0.00008;    fResSigy[3] =  0.0041197;




  fOffAx_J = 0.032;
  fOffAx_S = -0.0203;

/*  // base on run 633 (Saleve) and 635 (Jura)
  fRotX_J = -0.001501;
  fRotY_J = -0.003418;
  // geometry after fine tunning
  fRotX_J += -0.0000521;
  fRotY_J += -0.0001430;
  fRotX_S = -0.001737;
  fRotY_S = -0.009973;
  // geometry after fine tunning
  fRotX_S += 0.0000780;
  fRotY_S += -0.0001214;
  cout << "fRotY_J-fRotY_S(old) = " << fRotY_J - fRotY_S << endl;

  // fRotX_J =  -0.001501 - 0.0000521;
  // fRotY_J =   0.0;
  // fRotX_S =  -0.001737 + 0.0000780;
  // fRotY_S =  -0.0062;

  // based on 2313 and 2314 (before rotation)
  // fRotX_J =  -0.001936;
  // fRotY_J =  -0.0051375 + 0.005;

  // fRotX_S =  -0.002217;
  // fRotY_S =  -0.01184   + 0.005;

  // based on 2444 ans 2445
  fRotX_J = -0.002035;
  fRotY_J = -0.00083;
  fRotX_S = -0.002179;
  fRotY_S = -0.007494;
  
  cout << "fRotY_J-fRotY_S(new) = " << fRotY_J - fRotY_S << endl;
  cout << "fRotX_J =" << fRotX_J << "  fRotY_J =" << fRotY_J << endl;
  cout << "fRotX_S =" << fRotX_S << "  fRotY_S =" << fRotY_S << endl;
*/

	fRotX_J = -0.002032;
	fRotY_J = -0.001458;
	fRotX_S = -0.002350;
	fRotY_S = -0.007051;
	
	fRotX_J += 0.000018;
	fRotY_J += -0.000002;
	fRotX_S += -0.000018;
	fRotY_S += -0.000145;

  fRotZ_J = fN * 0.001;
  fRotZ_S = -fRotZ_J;

  // no arm rotations use it for test (needed to adjust arm offsets)
  // fRotX_J =  0.0;
  // fRotY_J =  0.0;

  // fRotX_S =  0.0;
  // fRotY_S =  0.0;



	double offsetX_J =0.280994 + 8.;
	double offsetY_J =-0.0384151;
	double offsetZ_J =0.252155;
	offsetX_J += 0.00233556;;     offsetY_J += -0.00052405;      offsetZ_J += 0.0131955;
	
		 if (fRunId >= 38945){	    offsetX_J += 0.00920726;     offsetY_J += 0.00199324;     offsetZ_J += 0.00147322;
									offsetX_J += 0.00119868;     offsetY_J += 0.000227227;     offsetZ_J += 0.000835319;
  } else if (fRunId >= 38943){ 	    offsetX_J += 0.00926197;     offsetY_J += 0.00203935;     offsetZ_J += 0.00074895;
									offsetX_J += 0.00129201;     offsetY_J += 0.00024587;     offsetZ_J += 0.000511727;
  } else if (fRunId >= 38941){ 	    offsetX_J += 0.00924322;     offsetY_J += 0.00202394;     offsetZ_J += 0.000440922;
									offsetX_J += 0.00132687;     offsetY_J += 0.000241412;     offsetZ_J += 0.000722253;
  } else if (fRunId >= 38938){ 	    offsetX_J += 0.00928477;     offsetY_J += 0.00207556;     offsetZ_J += 0.00116539;
									offsetX_J += 0.00126677;     offsetY_J += 0.000229519;     offsetZ_J += 0.000338565;
  } else if (fRunId >= 38919){ 	    offsetX_J += 0.00925564;     offsetY_J += 0.00233565;     offsetZ_J += 0.00294086;
									offsetX_J += 0.0012904;     offsetY_J += 0.000291826;     offsetZ_J += 0.00074745;
  } else if (fRunId >= 38917){ 	    offsetX_J += 0.00886902;     offsetY_J += 0.00231905;     offsetZ_J += 0.00251222;
									offsetX_J += 0.00129543;     offsetY_J += 0.000287318;     offsetZ_J += 0.000274277;
  } else if (fRunId >= 38898){ 	    offsetX_J += 0.00885863;     offsetY_J += 0.00267932;     offsetZ_J += 0.00449022;
									offsetX_J += 0.00117712;     offsetY_J += 0.000312648;     offsetZ_J += 0.000111426;
  } else if (fRunId >= 38889){ 	    offsetX_J += 0.0089275;     offsetY_J += 0.00242236;     offsetZ_J += 0.00289912;
									offsetX_J += 0.00125461;     offsetY_J += 0.000290643;     offsetZ_J += 0.000546541;
  } else if (fRunId >= 38881){ 	    offsetX_J += 0.0091949;     offsetY_J += 0.00234892;     offsetZ_J += 0.0025395;
									offsetX_J += 0.00127602;     offsetY_J += 0.000283251;     offsetZ_J += 8.46082e-05;
  } else if (fRunId >= 38835){ 	    offsetX_J += 0.00877307;     offsetY_J += 0.00262859;     offsetZ_J += 0.00442805;
									offsetX_J += 0.00122701;     offsetY_J += 0.000277885;     offsetZ_J += 0.000397805;
  } else if (fRunId >= 38826){ 	    offsetX_J += 0.00905756;     offsetY_J += 0.00261816;     offsetZ_J += 0.00432928;
									offsetX_J += 0.00124863;     offsetY_J += 0.000317521;     offsetZ_J += 0.000139563;
  } else if (fRunId >= 38825){ 	    offsetX_J += 0.00906774;     offsetY_J += 0.00260485;     offsetZ_J += 0.00478514;
									offsetX_J += 0.00128405;     offsetY_J += 0.00030515;     offsetZ_J += 0.000230474;
  } else if (fRunId >= 38813){ 	    offsetX_J += 0.00958664;     offsetY_J += 0.0024196;     offsetZ_J += 0.00318334;
									offsetX_J += 0.00128621;     offsetY_J += 0.000259691;     offsetZ_J += 0.000587603;
  } else if (fRunId >= 38812){ 	    offsetX_J += 0.0096364;     offsetY_J += 0.00239719;     offsetZ_J += 0.00342127;
									offsetX_J += 0.00130287;     offsetY_J += 0.000272061;     offsetZ_J += 0.000730204;
  } else if (fRunId >= 38804){ 	    offsetX_J += 0.0103665;     offsetY_J += 0.00261034;     offsetZ_J += 0.00506222;
									offsetX_J += 0.00147669;     offsetY_J += 0.000295553;     offsetZ_J += 0.000483211;
  } else if (fRunId >= 38796){ 	    offsetX_J += 0.0104453;     offsetY_J += 0.00253483;     offsetZ_J += 0.00502352;
									offsetX_J += 0.00142781;     offsetY_J += 0.000276263;     offsetZ_J += 0.000675352;
  } else if (fRunId >= 38775){ 	    offsetX_J += 0.0100413;     offsetY_J += 0.00251066;     offsetZ_J += 0.00434035;
									offsetX_J += 0.0013927;     offsetY_J += 0.000277781;     offsetZ_J += 0.000534314;
  } else if (fRunId >= 38772){ 	    offsetX_J += 0.0101401;     offsetY_J += 0.00260198;     offsetZ_J += 0.00493915;
									offsetX_J += 0.00142489;     offsetY_J += 0.000302988;     offsetZ_J += 0.000160014;
  } else if (fRunId >= 38771){ 	    offsetX_J += 0.0102541;     offsetY_J += 0.00263443;     offsetZ_J += 0.00520136;
									offsetX_J += 0.0014543;     offsetY_J += 0.000299588;     offsetZ_J += 0.000443747;
  } else if (fRunId >= 38764){ 	    offsetX_J += 0.0105209;     offsetY_J += 0.0026094;     offsetZ_J += 0.00533642;
									offsetX_J += 0.00152944;     offsetY_J += 0.000279262;     offsetZ_J += 0.000429797;
  } else if (fRunId >= 38761){ 	    offsetX_J += 0.0105767;     offsetY_J += 0.00252028;     offsetZ_J += 0.00491154;
									offsetX_J += 0.00147309;     offsetY_J += 0.000294444;     offsetZ_J += 0.000495238;
  } else if (fRunId >= 38759){ 	    offsetX_J += 0.0106913;     offsetY_J += 0.00249678;     offsetZ_J += 0.00483715;
									offsetX_J += 0.00145793;     offsetY_J += 0.000290197;     offsetZ_J += 0.000499339;
  } else if (fRunId >= 38739){ 	    offsetX_J += 0.0113715;     offsetY_J += 0.00269807;     offsetZ_J += 0.00682039;
									offsetX_J += 0.00155283;     offsetY_J += 0.000303897;     offsetZ_J += 0.00023624;
  } else if (fRunId >= 38707){ 	    offsetX_J += 0.0142706;     offsetY_J += 0.00240712;     offsetZ_J += 0.00530857;
									offsetX_J += 0.00211517;     offsetY_J += 0.00029789;     offsetZ_J += 0.000857144;
  } else if (fRunId >= 38667){ 	    offsetX_J += 0.0137476;     offsetY_J += 0.00214478;     offsetZ_J += 0.00490332;
									offsetX_J += 0.00205535;     offsetY_J += 0.000297175;     offsetZ_J += 0.000380281;
  } else if (fRunId >= 38662){ 	    offsetX_J += 0.0135728;     offsetY_J += 0.00217566;     offsetZ_J += 0.00453645;
								    offsetX_J += 0.00204483;     offsetY_J += 0.000278754;     offsetZ_J += 0.000738978;
  } else if (fRunId >= 38657){ 	    offsetX_J += 0.0136091;     offsetY_J += 0.00216991;     offsetZ_J += 0.00426257;
									offsetX_J += 0.00204875;     offsetY_J += 0.000264597;     offsetZ_J += 0.00086379;
  } else if (fRunId >= 38656){ 	    offsetX_J += 0.0137195;     offsetY_J += 0.00214536;     offsetZ_J += 0.00470383;
									offsetX_J += 0.00197183;     offsetY_J += 0.000276262;     offsetZ_J += 0.00105013;
  } else if (fRunId >= 38654){ 	    offsetX_J += 0.0137375;     offsetY_J += 0.00215316;     offsetZ_J += 0.00473067;
									offsetX_J += 0.00205284;     offsetY_J += 0.000286794;     offsetZ_J += 0.000955006;
  } else if (fRunId >= 38649){ 	    offsetX_J += 0.0141321;     offsetY_J += 0.0020962;     offsetZ_J += 0.00412893;
									offsetX_J += 0.00203609;     offsetY_J += 0.000257198;     offsetZ_J += 0.00116781;
  } else if (fRunId >= 38645){ 	    offsetX_J += 0.0142403;     offsetY_J += 0.00213499;     offsetZ_J += 0.00455562;
									offsetX_J += 0.00205575;     offsetY_J += 0.000262106;     offsetZ_J += 0.000754198;
  } else if (fRunId >= 38636){ 	    offsetX_J += 0.0144223;     offsetY_J += 0.00211949;     offsetZ_J += 0.00453819;
									offsetX_J += 0.00212408;     offsetY_J += 0.000270883;     offsetZ_J += 0.00102558;
  } else if (fRunId >= 38628){ 	    offsetX_J += 0.0143254;     offsetY_J += 0.00216105;     offsetZ_J += 0.00460395;
									offsetX_J += 0.00210492;     offsetY_J += 0.000290353;     offsetZ_J += 0.0010479;
  } else if (fRunId >= 38619){ 	    offsetX_J += 0.014248;     offsetY_J += 0.00225943;     offsetZ_J += 0.00513583;
									offsetX_J += 0.00211965;     offsetY_J += 0.000282013;     offsetZ_J += 0.00102268;
  } else if (fRunId >= 38606){ 	    offsetX_J += 0.0141937;     offsetY_J += 0.00228298;     offsetZ_J += 0.00560187;
									offsetX_J += 0.00213075;     offsetY_J += 0.000285876;     offsetZ_J += 0.000596906;
  } else if (fRunId >= 38599){ 	    offsetX_J += 0.0141906;     offsetY_J += 0.00228271;     offsetZ_J += 0.00518847;
									offsetX_J += 0.00215436;     offsetY_J += 0.000295073;     offsetZ_J += 0.00123197;
  } else if (fRunId >= 38593){ 	    offsetX_J += 0.0142285;     offsetY_J += 0.0022965;     offsetZ_J += 0.00523548;
								    offsetX_J += 0.00213143;     offsetY_J += 0.000303121;     offsetZ_J += 0.000733409;
  } else if (fRunId >= 38568){ 	    offsetX_J += 0.014496;     offsetY_J += 0.00217208;     offsetZ_J += 0.00440254;
									offsetX_J += 0.00214958;     offsetY_J += 0.000280517;     offsetZ_J += 0.000894633;
  } else if (fRunId >= 38562){ 	    offsetX_J += 0.0144605;     offsetY_J += 0.00222481;     offsetZ_J += 0.00491106;
									offsetX_J += 0.00212915;     offsetY_J += 0.000294556;     offsetZ_J += 0.000809484;
  } else if (fRunId >= 38548){ 	    offsetX_J += 0.014373;     offsetY_J += 0.00236719;     offsetZ_J += 0.00571031;
									offsetX_J += 0.00212351;     offsetY_J += 0.000298218;     offsetZ_J += 0.000567369;
  } else if (fRunId >= 38450){ 	    offsetX_J += 0.0145775;     offsetY_J += 0.00235363;     offsetZ_J += 0.00623194;
									offsetX_J += 0.00211083;     offsetY_J += 0.000302424;     offsetZ_J += 0.000666111;
//  } else if (fRunId >= 38444){ 	    
  } else if (fRunId >= 38436){ 	    offsetX_J += 0.0146932;     offsetY_J += 0.00235677;     offsetZ_J += 0.0064508;
									offsetX_J += 0.00219943;     offsetY_J += 0.000288744;     offsetZ_J += 0.000612408;
  }	  
 /*  } else if (fRunId >= 38300){ 	    
  } else if (fRunId >= 38295){ 	    
  } else if (fRunId >= 38291){ 	    
  } else if (fRunId >= 38290){ 	    
  } else if (fRunId >= 38283){ 	    
  } else if (fRunId >= 38278){ 	    
  } else if (fRunId >= 38275){ 	    
  } else if (fRunId >= 38242){ 	    
  } else if (fRunId >= 38240){ 	    
  } else if (fRunId >= 38238){ 	    
  } else if (fRunId >= 38233){
  } else if (fRunId >= 38232){ 	    
  }
*/
/*	
  if (fRunId >= 38945) {
		offsetX_J += 0.000227208;     offsetY_J += -1.20472e-05;     offsetZ_J += 0.000548831;
		offsetX_J += 7.37126e-05;     offsetY_J += 2.06645e-06;     offsetZ_J += 0.000108137;
  } else if (fRunId >= 38943){ 
		offsetX_J += 0.000531913;     offsetY_J += 1.12959e-06;     offsetZ_J += -0.000262248;
		offsetX_J += 9.14657e-05;     offsetY_J += -4.91306e-06;     offsetZ_J += -7.52776e-05;
  } else if (fRunId >= 38938){ 
		offsetX_J += 0.000551013;     offsetY_J += 3.64082e-05;     offsetZ_J += -3.61407e-05;
		offsetX_J += 0.000132533;     offsetY_J += -1.63729e-06;     offsetZ_J += 2.56093e-05;
  } else if (fRunId >= 38919) {
		offsetX_J += 0.00224351;     offsetY_J += 2.38463e-05;     offsetZ_J += 0.000581493;
		offsetX_J += 0.000834555;     offsetY_J += -4.37517e-07;     offsetZ_J += -6.82399e-05;
  } else if (fRunId >= 38917) {
		offsetX_J += 0.00166723;     offsetY_J += 2.06093e-05;     offsetZ_J += 0.000913328;
		offsetX_J += 0.000348243;     offsetY_J += -2.53227e-05;     offsetZ_J += 0.000459073;
  } else if (fRunId >= 38898) {
		offsetX_J += 0.00331093;     offsetY_J += 0.000108475;     offsetZ_J += 0.00277562;
		offsetX_J += 0.000700422;     offsetY_J += -3.54792e-05;     offsetZ_J += -0.000397397;
  } else if (fRunId >= 38889) {
		offsetX_J += 0.00220768;     offsetY_J += 4.52713e-05;     offsetZ_J += 0.00148123;
		offsetX_J += 0.000369923;     offsetY_J += -2.98409e-05;     offsetZ_J += -0.00043903;
  } else if (fRunId >= 38881) {
		offsetX_J += 0.00229786;     offsetY_J += 6.73082e-06;     offsetZ_J += 0.000319041;
		offsetX_J += 0.000418393;     offsetY_J += 6.20473e-06;     offsetZ_J += -0.000184241;
//  } else if (fRunId >= 38856) {
//		offsetX_J += 0.00268126;     offsetY_J += 9.59411e-05;     offsetZ_J += 0.00193463;
//  } else if (fRunId >= 38854){ 
//		offsetX_J += 0.00268126;     offsetY_J += 9.59411e-05;     offsetZ_J += 0.00193463;
  } else if (fRunId >= 38835) {
		offsetX_J += 0.00268126;     offsetY_J += 9.59411e-05;     offsetZ_J += 0.00193463;
		offsetX_J += 0.000740337;     offsetY_J += -5.26299e-06;     offsetZ_J += 0.000318951;
  } else if (fRunId >= 38826) {
		offsetX_J += 0.00332983;     offsetY_J += 6.86687e-05;     offsetZ_J += 0.00155641;
		offsetX_J += 0.000639434;     offsetY_J += -1.10308e-05;     offsetZ_J += 0.000363139;
  } else if (fRunId >= 38825) {
		offsetX_J += 0.00327896;     offsetY_J += 6.41577e-05;     offsetZ_J += 0.00116693;
		offsetX_J += 0.000713015;     offsetY_J += -2.82445e-06;     offsetZ_J += 0.000550134;
  } else if (fRunId >= 38824) {
		offsetX_J += 0.0030809;     offsetY_J += 5.89717e-05;     offsetZ_J += 0.00124656;
		offsetX_J += 0.00054541;     offsetY_J += -3.77472e-05;     offsetZ_J += -0.00039473;
  } else if (fRunId >= 38813) {
		offsetX_J += 0.00208312;     offsetY_J += 0.000154205;     offsetZ_J += 0.000576089;
		offsetX_J += 0.000382558;     offsetY_J += 1.0148e-05;     offsetZ_J += 2.28891e-05;
  }	 
	
	
	*/

	double offsetX_S =-offsetX_J;
	double offsetY_S =-offsetY_J;
	double offsetZ_S =-offsetZ_J;



  // correction based on run 1494 (0 field)
/*  double offsetX_J = 0.881427 + 8.;
  double offsetY_J = -0.0164678;
  double offsetZ_J = 0.317669;

  double offsetX_S = -0.881427 - 8.;
  double offsetY_S = 0.0164678;
  double offsetZ_S = -0.317669;

  // Setup offsets for PbPb if fitted to Vtx (in local arm frames)

  offsetX_J += -0.00169475;
  offsetY_J += -0.0244815;
  offsetZ_J += 0.0332515;

  offsetX_S += 0.00169475;
  offsetY_S += 0.0244815;
  offsetZ_S += -0.0332515;

  // form VdTrackingHT module
  offsetX_J += 0.000539333;
  offsetY_J += -0.000280796;
  offsetZ_J += 0.00251891;

  offsetX_S += -0.000539333;
  offsetY_S += 0.000280796;
  offsetZ_S += -0.00251891;

  // Add correction for PbPb if fitted to Vtx_Glob
  offsetX_J += -0.180561;
  offsetY_J += -0.0111728;
  offsetZ_J += -0.0786263;
  offsetX_S += 0.180561;
  offsetY_S += 0.0111728;
  offsetZ_S += 0.0786263;

  // from 2279 (2mm target)
  offsetX_J += 0.168806;
  offsetY_J += 0.00158409;
  offsetZ_J += 0.0146389;
  offsetX_S += -0.168806;
  offsetY_S += -0.00158409;
  offsetZ_S += -0.0146389;

  // from 2315 (2mm target)
  offsetX_J += 0.0245236;
  offsetY_J += 0.000223332;
  offsetZ_J += -0.012452;
  offsetX_S += -0.0245236;
  offsetY_S += -0.000223332;
  offsetZ_S += 0.012452;

  // from 2320 (2mm target)
  offsetX_J += -0.371642;
  offsetY_J += 0.00026657;
  offsetZ_J += 0.00134117;
  offsetX_S += 0.371642;
  offsetY_S += -0.00026657;
  offsetZ_S += -0.00134117;

  // from 2435 (3mm target)
  offsetX_J += -0.0235636;
  offsetY_J += 0.00677256;
  offsetZ_J += -0.0162827;
  offsetX_S += 0.0235636;
  offsetY_S += -0.00677256;
  offsetZ_S += 0.0162827;

  // from 2439 (3mm target)
  offsetX_J += -0.00221664;
  offsetY_J += 0.00269655;
  offsetZ_J += -0.0095972;
  offsetX_S += 0.00221664;
  offsetY_S += -0.00269655;
  offsetZ_S += 0.0095972;

  // from 2448 (3mm)
  offsetX_J += -0.0683482;
  offsetY_J += -0.00850475;
  offsetZ_J += 0.00308394;
  offsetX_S += 0.0683482;
  offsetY_S += 0.00850475;
  offsetZ_S += -0.00308394;

  offsetX_J += -0.072222;
  offsetX_S += 0.072222;

  // based on runs 38938, 38943, 38945
  offsetX_J += -0.00436184;
  offsetY_J += 0.000333734;
  offsetZ_J += -0.00347789;
  offsetX_S += 0.00436184;
  offsetY_S += -0.000333734;
  offsetZ_S += 0.00347789;
  
*/  

  // parameter to extract Jura and Saleve events with consistent vertexes.
  fVtxDzSigma = 0.075;

  fJuraArmOffset.SetX(offsetX_J);
  fJuraArmOffset.SetY(offsetY_J);
  fJuraArmOffset.SetZ(offsetZ_J);
  fSaleveArmOffset.SetX(offsetX_S);
  fSaleveArmOffset.SetY(offsetY_S);
  fSaleveArmOffset.SetZ(offsetZ_S);
}





//_____________________________________________________________
void Na61VdParameters::PrintSetupInfo() {}

//_____________________________________________________________

unsigned int Na61VdParameters::FindCalibRuns(const unsigned int run) {
  const std::vector<unsigned int> calibruns = {38945, 38943, 38941, 38938, 38919, 38917, 38898, 38889, 38881, 38835, 38826, 38825, 38813, 38812, 38804, 38796, 38775, 38772, 38771, 38764, 38761, 38759, 38739, 38707, 
	                                           38667, 38662, 38657, 38656, 38654, 38649, 38645, 38636, 38628, 38619, 38606, 38599, 38593, 38568, 38562, 38548, 38450, 38436,//38444, 
	                                           //38300, 38295, 38291, 38290, 38283, 38278, 38275, 38242, 38240, 38238, 38233, 
											   //38945,38943,38938,38919,38917,38898,38889,38881,38835,38826,38825,38824,38813,  //38000,  //PbPb150 2018
	                                           33755, 33751, 33747, 33739, 33737, 33696, 33693, 33687, 33683, 33680, 33646, 33635, 33633, 33613, 33585, 33549, 33539, 33535, 33500, 33498, 33495, 33493, 33490,
	                                           33484, 33481, 33480, 33469, 33464, 33461, 33453,     33400,  //XeLa75 
	                                           33375, 33358, 33343, 33320, 33307, 33300, 33275, 33264, 33227, 33201, 33187,  33172, 33157, 33128, 33078, 32975, 32957, 32910, 32848, 32728, 32601}; //XeLa150
	                                           //33375, 33368, 33358, 33343, 33341, 33320, 33310, 33309, 33307, 33303, 33302, 33301, 33300, 33275, 33273, 33266, 33264, 33259, 33227, 33201, 
	                                           //33187, 33186, 33183, 33180, 33172, 33157, 33128, 33124, 33081, 33078, 32975, 32969, 32960, 32958, 32957, 32910, 32752, 32747, 32745, 32728, 32601}; //old production XeLa150
  unsigned int calibrun = *calibruns.rbegin();
  for (auto iterator = calibruns.rbegin(); iterator != calibruns.rend() && *iterator <= run; iterator++) calibrun = *iterator;

  if ((run<=33395)&&(run >=32601)){ //xela150
        calibrun=33201;
  }   


    ostringstream info;
    info << " Using calibrated run " << calibrun;
    INFO(info);
    
  return calibrun;
}
//_____________________________________________________________
//_____________________________________________________________
//_____________________________________________________________
//_____________________________________________________________
//_____________________________________________________________

/*
//_______________________________________________________________________________________________________
void Na61VdParameters::FindShineRun(int runid, int& calrunid, int& runid_shine,int& chunks)
{

  int hld_arr[] = {    8,   35,   45,   54,   55,  56,    57,   60,   61,   63,   64,   65,   66,   68,   69,   70,   71,
                         72,   74,   80,   81,   82,   86,   87,   88,   89,   91,   92,   93,  120,  122,  123,  124,  135,
                         142,  144,  145,  146,  147,  151,  152,  153,  154,  157,  158,  159, 160,   161,  165,  166,  167,
                         168,  169,  170,  171,  172,  174,  175,  176,
                         //
                         1350, 1351, 1352, 1353, 1354, 1355, 1356, 1357, 1358, 1360, 1361, 1362, 1363, 1364, 1365, 1366, 1367,
                         1368, 1369, 1370,
                         //
                         1371, 1372, 1373, 1375, 1377, 1378, 1379, 1381, 1382, 1383, 1384, 1385, 1386, 1387, 1388, 1390, 1391,
                         1392, 1393, 1394,
                         //
                         1395, 1397, 1398, 1399, 1400, 1407, 1408, 1409, 1410, 1411, 1412, 1414, 1415, 1416, 1417, 1418, 1419,
                         1420, 1421, 1422, 1423, 1424, 1425, 1426, 1427, 1428, 1429, 1430, 1431, 1432, 1433, 1434, 1435, 1436,
                         1437, 1438, 1439, 1440, 1441, 1443, 1444, 1445, 1447, 1448, 1449, 1450,
                         //
                         1547, 1578};

  int shine_arr[]={27256,27286,27291,27300,27301,27302,27303,27306,27308,27310,27312,27313,27318,27320,27321,27326,27328,
                     27329,27331,27341,27342,27343,27347,27348,27349,27350,27352,27353,27355,27384,27387,27388,27389,27401,
                     27411,27412,27416,27417,27418,27422,27423,27424,27425,27428,27429,27430,27431,27432,27436,27438,27440,
                     27441,27442,27443,27446,27447,27449,27450,27452,
                     //
                     33201,33202,33203,33204,33207,33208,33209,33210,33211,33213,33214,33222,33223,33224,33225,33226,33227,
                     33228,33229,33230,
                     //
                     33231,33232,33233,33238,33243,33244,33245,33250,33251,33252,33257,33258,33259,33260,33261,33264,33265,
                     33266,33271,33272,
                     //
                     33273,33275,33276,33277,33278,33289,33290,33295,33296,33297,33298,33300,33301,33302,33303,33304,33305,
                     33306,33307,33308,33309,33310,33311,33312,33313,33314,33315,33316,33317,33318,33319,33320,33321,33322,
                     33323,33324,33331,33332,33333,33335,33336,33337,33339,33340,33341,33342,
                     //
                     33495, 33536};

  int chunk_arr[]= { 126,   50,   13,  171,   65,   17,   20,   47,   78,   20,   10,   93,   17,   29,   42,   53,   34,
                       123,  135,   53,  244,   38,   15,   31,   92,   41,   28,    4,   51,  200,   38,   34,   88,   42,
                         2,    9,    5,    1,    7,  120,    3,  116,  119,  269,   53,   16,  180,   19,   24,   17,  161,
                        61,   57,    6,   39,   78,   26,   11,    5,
                       //
                       121,  41,   39,   65,   42,   22,   38,   39,   19,   31,   78,    41,   88,    7,   71,   48,   157,
                       37,   12,   69,
                       //
                       21,   33,    8,   29,   53,   42,    47,   74,   75,   96,   70,    83,   138,  16,  43,  125,    39,
                       163,  61, 52,
                       //
                       127, 227,   37,   28,   72,   23,    28,   34,   20,   16,   67,   128,   105, 336,  218,  89,    62,
                       7,   146,   63,  115,  116,   25,    54,   29,   29,   16,   69,    53,    20,  58,  103,  25,    75,
                       65,   37,   24,   14,   41,   74,    28,   38,   73,   10,  179,    25,
                       //
                       121,  87};

  int hld_calarr[] = { 168,
                         //
                         1350, 1351, 1353, 1354, 1356, 1357, 1360, 1361, 1362, 1363, 1366, 1367,
                         1368, 1370,
                         //
                         1377, 1381, 1383, 1386, 1390, 1392, 1394,
                         //
                         1395, 1397, 1400, 1409, 1412, 1414, 1416, 1417, 1421, 1424, 1430, 1434, 1436, 1443, 1447, 1449};


  int N = sizeof(hld_arr)/4;
  cout<<" FindShineRun: number of all physics runs is "<<N<<endl;

  for(int i=0;i<N;i++){
    if(hld_arr[i] == runid){
      runid_shine = shine_arr[i];
      chunks = chunk_arr[i];
      break;
    }
  }

  N = sizeof(hld_calarr)/4;
  for(int i=0;i<N;i++){
    if(hld_calarr[i] == runid)calrunid = hld_calarr[i];
    if(hld_calarr[i] > runid){
      if(i>0)calrunid = hld_calarr[i-1];
      else calrunid = hld_calarr[0];
      break;
    }
  }


  if(!runid_shine || !chunks)cout<<" FindShineRun: something wring in number assosiation: runid_shine="<<runid_shine<<
                               " chunks="<<chunks<<endl;

  cout<<" FindShineRun: numbers assosiation: runid_shine="<<runid_shine<<" chunks="<<chunks<<endl;

}
*/

//____________________________________________________________________
ostream& operator<<(ostream& os, Na61VdParameters* armpars) {
  os << " Jura X offset:" << armpars->GetJuraArmOffset().X() << " Saleve X offset: " << armpars->GetSaleveArmOffset().X() << endl;
  return os;
}
