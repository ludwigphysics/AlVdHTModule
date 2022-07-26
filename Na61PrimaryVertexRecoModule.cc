//--------------------------------------------
// Primary vertex reco module
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61PrimaryVertexRecoModule
#include "Na61PrimaryVertexRecoModule.h"
#endif

#ifndef NA61_Na61VdParametersManager
#include "Na61VdParametersManager.h"
#endif

#ifndef ROOT_TDirectory
#include "TDirectory.h"
#endif
#ifndef UTIL_UDataTable
#include "UDataTable.h"
#endif
#ifndef UTIL_UG4Hit
#include "UG4Hit.h"
#endif
#ifndef UTIL_USensorPixel
#include "USensorPixel.h"
#endif

#ifndef ROOT_TF1
#include "TF1.h"
#endif
#ifndef ROOT_TMath
#include "TMath.h"
#endif
#ifndef ROOT_TGraphErrors
#include "TGraphErrors.h"
#endif

#include "TSystem.h"
#include "TDirectory.h"
#include "Math/QuantFuncMathCore.h"
#include <time.h>
#include "iostream"

//____________________________________________________________________
// ClassImp(Na61PrimaryVertexRecoModule);

//____________________________________________________________________
Na61PrimaryVertexRecoModule::Na61PrimaryVertexRecoModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
}

//____________________________________________________________________
Na61PrimaryVertexRecoModule::Na61PrimaryVertexRecoModule(const char* name, const char* title) : Na61Module(name, title) {
  // Named Constructor
  SetState(kSetup);
  fNsigS = 12;
  fNsigJ = 14;
  fZprim = 73.1;
}

//____________________________________________________________________
void Na61PrimaryVertexRecoModule::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  //TDirectory* histDir = gDirectory->mkdir("PrimaryVertexReco");
  TDirectory* histDir = gDirectory->GetDirectory("AlPrimaryVertexReco");

  histDir->cd();

  fhCurvature_back = new TH1F("hCurvature_back", "", 500, -0.00003, 0.00003);

  fhVtxX = new TH1F("hVtxX", "", 300, -3., 3);
  fhVtxY = new TH1F("hVtxY", "", 300, -3., 3);
  fhTracksVsVtxX_J = new TH1F("hTracksVsVtxX_J", "", 30, -3., 3);
  fhTracksVsVtxY_J = new TH1F("hTracksVsVtxY_J", "", 30, -3., 3);
  fhTracksVsVtxX_S = new TH1F("hTracksVsVtxX_S", "", 30, -3., 3);
  fhTracksVsVtxY_S = new TH1F("hTracksVsVtxY_S", "", 30, -3., 3);

  fhVtxDx = new TH1F("hVtxDx", "", 1000, -2.5, 2.5);
  fhVtxDy = new TH1F("hVtxDy", "", 1000, -2.5, 2.5);
  fhVtxDz = new TH1F("hVtxDz", "", 1000, -10., 10.);

  //fhVtxDx_Glob = new TH1F("hVtxDx_Glob", "", 1000, -2.5, 2.5);
  fhVtxDx_Glob = new TH1F("hVtxDx_Glob", "", 2000, -5.5, 5.5);
  fhVtxDy_Glob = new TH1F("hVtxDy_Glob", "", 1000, -2.5, 2.5);
  // fhVtxDx_Glob = new TH1F("hVtxDx_Glob","", 1000,-20.5,20.5);
  // fhVtxDy_Glob = new TH1F("hVtxDy_Glob","", 1000,-20.5,20.5);
  fhVtxDz_Glob = new TH1F("hVtxDz_Glob", "", 1000, -10., 10.);

  fhVtxStatus_JuraVsSaleve = new TH2F("VtxStatus_JuraVsSaleve", "", 5, 0, 5, 5, 0, 5);
  fhFullTracks_JvsS = new TH2F("hFullTracks_JvsS", "", 200, 0, 500, 200, 0, 500);
  fhFullTracks_JS = new TH1F("hFullTracks_JS", "", 1000, 0, 1000);
  fhFullTracks_J = new TH1F("hFullTracks_J", "", 500, 0, 500);
  fhFullTracks_S = new TH1F("hFullTracks_S", "", 500, 0, 500);
  fhFullTracks4h_JS = new TH1F("hFullTracks4h_JS", "", 300, 0, 300);
  fhFullTracks4h_J = new TH1F("hFullTracks4h_J", "", 300, 0, 300);
  fhFullTracks4h_S = new TH1F("hFullTracks4h_S", "", 300, 0, 300);

  TString str[] = {"", "_J", "_S"};

  for (int i = 0; i < 3; i++) {
    fhRecoVertexZ[i] = new TH1F(Form("hRecoVertexZ%s", str[i].Data()), " ", 2000, -500., 100.);
    fhRecoVertexZ_fine[i] = new TH1F(Form("hRecoVertexZ_fine%s", str[i].Data()), " ", 1000, -90., -40.);
    fhRecoVertexXY[i] = new TH2F(Form("hRecoVertexXY%s", str[i].Data()), " ", 500, -50., 50., 500, -50., 50.);

    fhRecoVertexZ_flag1[i] = new TH1F(Form("hRecoVertexZ_flag1%s", str[i].Data()), " ", 2000, -2000., 100.);
    fhRecoVertexZ_fine_flag1[i] = new TH1F(Form("hRecoVertexZ_fine_flag1%s", str[i].Data()), " ", 1000, -90., -40.);
    fhRecoVertexXY_flag1[i] = new TH2F(Form("hRecoVertexXY_flag1%s", str[i].Data()), " ", 500, -50., 50., 500, -50., 50.);
  }

  fhAxJ = new TH1F(Form("hAxJ"), "", 1500, -0.25, 0.25);
  fhAyJ = new TH1F(Form("hAyJ"), "", 1500, -0.25, 0.25);
  fhAxS = new TH1F(Form("hAxS"), "", 1500, -0.25, 0.25);
  fhAyS = new TH1F(Form("hAyS"), "", 1500, -0.25, 0.25);

  fhAxJ_flag02 = new TH1F(Form("hAxJ_flag02"), "", 1500, -0.25, 0.25);
  fhAyJ_flag02 = new TH1F(Form("hAyJ_flag02"), "", 1500, -0.25, 0.25);

  fhAxJ_flag0 = new TH1F(Form("hAxJ_flag0"), "", 1500, -0.25, 0.25);
  fhAyJ_flag0 = new TH1F(Form("hAyJ_flag0"), "", 1500, -0.25, 0.25);
  fhAxJ_flag1 = new TH1F(Form("hAxJ_flag1"), "", 1500, -0.25, 0.25);
  fhAxJ_flag2 = new TH1F(Form("hAxJ_flag2"), "", 1500, -0.25, 0.25);
  fhAyJ_flag2 = new TH1F(Form("hAyJ_flag2"), "", 1500, -0.25, 0.25);

  fhAxS_flag02 = new TH1F(Form("hAxS_flag02"), "", 1500, -0.25, 0.25);
  fhAyS_flag02 = new TH1F(Form("hAyS_flag02"), "", 1500, -0.25, 0.25);

  fhAxS_flag0 = new TH1F(Form("hAxS_flag0"), "", 1500, -0.25, 0.25);
  fhAyS_flag0 = new TH1F(Form("hAyS_flag0"), "", 1500, -0.25, 0.25);
  fhAxS_flag1 = new TH1F(Form("hAxS_flag1"), "", 1500, -0.25, 0.25);
  fhAxS_flag2 = new TH1F(Form("hAxS_flag2"), "", 1500, -0.25, 0.25);
  fhAyS_flag2 = new TH1F(Form("hAyS_flag2"), "", 1500, -0.25, 0.25);

  for (int i = 0; i < 8; i++) {
    fhClusters_JvsS[i] = new TH2F(Form("hClusters_%s_JvsS", fAlSensorNames[i].Data()), "", 300, 0, 300, 300, 0, 300);
    fhClusters_JmS[i] = new TH2F(Form("hClusters_%s_JmS", fAlSensorNames[i].Data()), "", 300, 0, 300, 200, -100, 100);
  }

  // these histos should be present in the HitFinder, however we need them here to test/debug simulation
  TString arm_str[2] = {"J","S"};
  for(Int_t arm=0; arm<2; arm++){
    fhZX_fine[arm] =	new TH2F(Form("hZX_fine_%s",arm_str[arm].Data()),"",1000,-6,160,1000,-50,50);
    fhXY_fine[arm] =	new TH2F(Form("hXY_fine_%s",arm_str[arm].Data()),"",1000,-50,50,1000,-80,80);
    fhX_Al1[arm] =	new TH1F(Form("hX_Al1_%s",arm_str[arm].Data()),"",1000,-50,50);
    fhX_Al2[arm] =	new TH1F(Form("hX_Al2_%s",arm_str[arm].Data()),"",1000,-50,50);
    fhX_Al3[arm] =	new TH1F(Form("hX_Al3_%s",arm_str[arm].Data()),"",1000,-50,50);
    fhX_Al4[arm] =	new TH1F(Form("hX_Al4_%s",arm_str[arm].Data()),"",1000,-50,50);
    fhY_Al1[arm] =	new TH1F(Form("hY_Al1_%s",arm_str[arm].Data()),"",1000,-80,80);
    fhY_Al2[arm] =	new TH1F(Form("hY_Al2_%s",arm_str[arm].Data()),"",1000,-80,80);
    fhY_Al3[arm] =	new TH1F(Form("hY_Al3_%s",arm_str[arm].Data()),"",1000,-80,80);
    fhY_Al4[arm] =	new TH1F(Form("hY_Al4_%s",arm_str[arm].Data()),"",1000,-80,80);
  }

  gDirectory->cd("..");
}

//____________________________________________________________________
void Na61PrimaryVertexRecoModule::Init() {
  // Job-level initialisation
  SetState(kInit);

  Na61Module::Init(); //to initialize stuff like fNDev and fNDevComb

  fParams = (Na61VdParametersManager::Instance())->GetVdParams();
  fParams->SetRunId(fRunId);
  if (!fParams->GetInit()) fParams->Init();

  fRotY_J = fParams->GetJuraArmRotY();
  fRotY_S = fParams->GetSaleveArmRotY();

  fOffAx_J = fParams->GetJuraArmOffAx();
  fSigAx_J = fParams->GetJuraArmSigAx();
  fOffAy_J = fParams->GetJuraArmOffAy();
  fSigAy_J = fParams->GetJuraArmSigAy();

  fOffAx_S = fParams->GetSaleveArmOffAx();
  fSigAx_S = fParams->GetSaleveArmSigAx();
  fOffAy_S = fParams->GetSaleveArmOffAy();
  fSigAy_S = fParams->GetSaleveArmSigAy();

  fZprim = 73.1;
}

//____________________________________________________________________
void Na61PrimaryVertexRecoModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);
}
/*
//____________________________________________________________________
void Na61PrimaryVertexRecoModule::SetRunInfo(TFile* file)
{

  TDirectory* dir = file->mkdir("runInfo");
  dir->cd();

  TH1F* hVolumesX = new TH1F("volumesX","",8,0,8);
  TH1F* hVolumesY = new TH1F("volumesY","",8,0,8);
  TH1F* hVolumesZ = new TH1F("volumesZ","",8,0,8);

  TH1F* hRotX = new TH1F("rotX","",8,0,8);
  TH1F* hRotY = new TH1F("rotY","",8,0,8);
  TH1F* hRotZ = new TH1F("rotZ","",8,0,8);

  for(int i=0;i<8;i++){
    hVolumesX->SetBinContent(i+1,fParams->GetVolumeX(i));
    hVolumesY->SetBinContent(i+1,fParams->GetVolumeY(i));
    hVolumesZ->SetBinContent(i+1,fParams->GetVolumeZ(i));
    hRotX->SetBinContent(i+1,fParams->GetRotX(i));
    hRotY->SetBinContent(i+1,fParams->GetRotY(i));
    hRotZ->SetBinContent(i+1,fParams->GetRotZ(i));
  }
  hVolumesX->Write();
  hVolumesY->Write();
  hVolumesZ->Write();
  hRotX->Write();
  hRotY->Write();
  hRotZ->Write();

  dir->cd("/");

}
*/

//____________________________________________________________________
void Na61PrimaryVertexRecoModule::Event(UEventNode* inNodeJ, UEventNode* inNodeS)

{
  // Per event method
  SetState(kEvent);

  fAllTracksJ = 0;
  fAllTracksS = 0;

  FillClusterCorrelation(inNodeJ, inNodeS);

  double pvtxX_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexX();
  double pvtxY_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexY();
  double pvtxZ_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexZ();

  double pvtxX_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexX();
  double pvtxY_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexY();
  double pvtxZ_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexZ();

  int pvertStatusJ = ((UVdEvent*)inNodeJ)->GetPrimaryVertexStatus();
  int pvertStatusS = ((UVdEvent*)inNodeS)->GetPrimaryVertexStatus();

  fhVtxStatus_JuraVsSaleve->Fill(pvertStatusJ, pvertStatusS);

  //cout << " PrimaryVertexReco: pvertStatusJ=" << pvertStatusJ << " pvertStatusS=" << pvertStatusS << "  " << pvtxZ_J << "  " << pvtxZ_S << endl;

  if (pvertStatusJ == 2 && pvertStatusS == 2) {
    // if(pvertStatusJ > 0 && pvertStatusS > 0 ){ // be less restricted for h+Pb
    fhVtxDx->Fill(pvtxX_J - pvtxX_S + 16.);
    fhVtxDy->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz->Fill(pvtxZ_J - pvtxZ_S);
    //cout << "Jura Xpos= " << pvtxX_J << " Saleve Xpos=" << pvtxX_S << " pvtxX_J - pvtxX_S + 16. = " << pvtxX_J - pvtxX_S + 16. << endl;
  }

  LocalToGlobal(0, inNodeJ);
  LocalToGlobal(1, inNodeS);

  fJuraProdCounter = 0;
  fSaleveProdCounter = 0;

  if (fRunId < 300) {
    ClassifyTracks(inNodeJ, inNodeS);
  } else if (fRunId < 600) {
    ClassifyTracks_pPb(inNodeJ, inNodeS);
  } else {
    ClassifyTracks(inNodeJ, inNodeS);
  }

  // if(fJuraProdCounter>2 && fSaleveProdCounter>2){
  // cout<<" JuraProdCounter="<<fJuraProdCounter<<" SaleveProdCounter="<<fSaleveProdCounter<<endl;
  //}

  FindPrimaryVertex(1, inNodeJ, inNodeS);
  FindPrimaryVertex(2, inNodeJ, inNodeS);
  FindPrimaryVertex(2, inNodeJ, 0);
  FindPrimaryVertex(2, 0, inNodeS);

  pvtxX_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexX();
  pvtxY_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexY();
  pvtxZ_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexZ();

  pvtxX_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexX();
  pvtxY_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexY();
  pvtxZ_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexZ();

  /*
// check it out what is ProdCounter???
  if(fJuraProdCounter>2 && fSaleveProdCounter>2){
    //cout<<" PrimaryVertexReco4: pvertStatusJ="<<pvertStatusJ<<" pvertStatusS="<<pvertStatusS
    //	    <<"  "<<pvtxZ_J<<"  "<<pvtxZ_S<<"  "<<pvtxX_J-pvtxX_S<<endl;

    fhVtxDx_Glob->Fill(pvtxX_J - pvtxX_S);
    fhVtxDy_Glob->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz_Glob->Fill(pvtxZ_J - pvtxZ_S);
  }
  */

  if (pvertStatusJ == 2 && pvertStatusS == 2) {
    // if(pvertStatusJ > 0 && pvertStatusS > 0){  // be less restricted for h+Pb
    fhVtxDx_Glob->Fill(pvtxX_J - pvtxX_S);
    fhVtxDy_Glob->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz_Glob->Fill(pvtxZ_J - pvtxZ_S);
  }

  if ((pvertStatusJ == 0) && pvertStatusS) {
    ((UVdEvent*)inNodeJ)->SetPrimaryVertexStatus(3);
    ((UVdEvent*)inNodeJ)->SetPrimaryVertex(pvtxX_S, pvtxY_S, pvtxZ_S);
  }

  if ((pvertStatusS == 0) && pvertStatusJ) {
    ((UVdEvent*)inNodeS)->SetPrimaryVertexStatus(3);
    ((UVdEvent*)inNodeS)->SetPrimaryVertex(pvtxX_J, pvtxY_J, pvtxZ_J);
  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::LocalToGlobalAl(int iarm, UEventNode* out) {
  // rotation (track by track):
  // 1. define two points on track
  // 2. transform in z to middle point of arm (dz=75)
  // 3. rotate arm by roty and rotx
  // 4. transform two points back from middle to primary frame (dz=-75)
  // 5. setup track
  //
  // translation (track by track):
  // take track by track and shift it origin. Important: define origin at plane z=0.
  //
  // Note:
  // iarm = 0/1 (Jura/Saleve)

  double rotX = 0;
  double rotY = 0;
  double rotZ = 0;
  double OffsetX = 0;
  double OffsetY = 0;
  double OffsetZ = 0;

  // check param file !!!!!!!!!

  if (iarm == 0) {
    rotX = fParams->GetJuraArmRotX();
    rotY = fParams->GetJuraArmRotY();
    rotZ = fParams->GetJuraArmRotZ();
    OffsetX = fParams->GetJuraArmOffset().X();
    OffsetY = fParams->GetJuraArmOffset().Y();
    OffsetZ = fParams->GetJuraArmOffset().Z();
  } else if (iarm == 1) {
    rotX = fParams->GetSaleveArmRotX();
    rotY = fParams->GetSaleveArmRotY();
    rotZ = fParams->GetSaleveArmRotZ();
    OffsetX = fParams->GetSaleveArmOffset().X();
    OffsetY = fParams->GetSaleveArmOffset().Y();
    OffsetZ = fParams->GetSaleveArmOffset().Z();
  } else {
    Error("LocalToGlobalAl", "Wrong arm id, check it out");
  }

  //cout<<"iarm="<<iarm<<"  "<<rotX<<" "<<rotY<<" "<<rotZ<<" "<<OffsetX<<" "<<OffsetY<<" "<<OffsetZ<<endl;

  double sa = TMath::Sin(rotZ);
  double ca = TMath::Cos(rotZ);
  double sb = TMath::Sin(rotY);
  double cb = TMath::Cos(rotY);
  double sg = TMath::Sin(rotX);
  double cg = TMath::Cos(rotX);

  // apply rotations
  for(Int_t i=0;i<fNDevComb;i++){
    UDataTable* tracktab = out -> GetDataTable(Form("FullTracks %s",fDevCombName[i].Data()));
    if (!tracktab) continue;

    //cout<<"table: "<<tracktab->GetName()<<"  "<<tracktab->GetEntries()<<endl;
    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      track->Activate();

      fhCurvature_back->Fill(track->GetCurvature());

      double x0 = track->GetXatZ(0);
      double y0 = track->GetYatZ(0);
      double z0 = 0;
      double x1 = track->GetXatZ(100);
      double y1 = track->GetYatZ(100);
      double z1 = 100;
      // transform to arm middle plane
      z0 = z0 - 75.0;
      z1 = z1 - 75.0;

      // XY
      // double xx0 =      cb * x0                +     sb * z0;
      // double yy0 =  -sg*sb * x0   +   cg * y0  +  sg*cb * z0;
      // double zz0 =  -cg*sb * x0   -   sg * y0  +  cg*cb * z0;

      // double xx1 =    (cb) * x1                +     sb * z1;
      // double yy1 =  -sg*sb * x1   +   cg * y1  +  sg*cb * z1;
      // double zz1 =  -cg*sb * x1   -   sg * y1  +  cg*cb * z1;

      // ZXY
      double xx0 = (ca * cb - sa * sg * sb) * x0 + sa * cg * y0 + (ca * sb + sa * sg * cb) * z0;
      double yy0 = -(sa * cb + ca * sg * sb) * x0 + ca * cg * y0 + (-sa * sb + ca * sg * cb) * z0;
      double zz0 = -cg * sb * x0 - sg * y0 + cg * cb * z0;

      double xx1 = (ca * cb - sa * sg * sb) * x1 + sa * cg * y1 + (ca * sb + sa * sg * cb) * z1;
      double yy1 = -(sa * cb + ca * sg * sb) * x1 + ca * cg * y1 + (-sa * sb + ca * sg * cb) * z1;
      double zz1 = -cg * sb * x1 - sg * y1 + cg * cb * z1;

      zz0 = zz0 + 75.0;
      zz1 = zz1 + 75.0;
      Vector3D origin(xx0, yy0, zz0);
      Vector3D direction(xx1 - xx0, yy1 - yy0, zz1 - zz0);

      (track->Getline())->SetOrigin(origin);
      (track->Getline())->SetDirection(direction);
    }
  }

  // apply offsets
  for(Int_t i=0;i<fNDevComb;i++){
    UDataTable* tracktab = out -> GetDataTable(Form("FullTracks %s",fDevCombName[i].Data()));
    if (!tracktab) continue;
    
    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      // track -> Activate();
      double ox = (track->Getline())->GetOrigin().X();
      double oy = (track->Getline())->GetOrigin().Y();
      double oz = (track->Getline())->GetOrigin().Z();
      // cout<<"oz="<<oz<<endl;

      (track->Getline())->SetOrigin(ox + OffsetX, oy + OffsetY, oz + OffsetZ);
      // double oxn = (track->Getline())->GetOrigin().X();
      // double oyn = (track->Getline())->GetOrigin().Y();
      // double ozn = (track->Getline())->GetOrigin().Z();
      // cout<<"oxn="<<oxn<<"  oyn="<<oyn<<" ozn="<<ozn<<endl;
      double x0 = track->GetXatZ(0);
      double y0 = track->GetYatZ(0);
      // cout<<"x0="<<x0<<"  y0="<<y0<<" z="<<0<<endl;
      (track->Getline())->SetOrigin(x0, y0, 0.);
    }
  }
  // cout<<" ------------------------------ iarm="<<iarm<<endl;
  TransformHitsAl(out, OffsetX, OffsetY, OffsetZ, sb, cb, sg, cg, sa, ca);
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::LocalToGlobal(int iarm, UEventNode* out) {
  // rotation (track by track):
  // 1. define two points on track
  // 2. transform in z to middle point of arm (dz=75)
  // 3. rotate arm by roty and rotx
  // 4. transform two points back from middle to primary frame (dz=-75)
  // 5. setup track
  //
  // translation (track by track):
  // take track by track and shift it origin. Important: define origin at plane z=0.
  //
  // Note:
  // iarm = 0/1 (Jura/Saleve)

  double rotX = 0;
  double rotY = 0;
  double rotZ = 0;
  double OffsetX = 0;
  double OffsetY = 0;
  double OffsetZ = 0;

  if (iarm == 0) {
    rotX = fParams->GetJuraArmRotX();
    rotY = fParams->GetJuraArmRotY();
    rotZ = fParams->GetJuraArmRotZ();
    OffsetX = fParams->GetJuraArmOffset().X();
    OffsetY = fParams->GetJuraArmOffset().Y();
    OffsetZ = fParams->GetJuraArmOffset().Z();
  } else if (iarm == 1) {
    rotX = fParams->GetSaleveArmRotX();
    rotY = fParams->GetSaleveArmRotY();
    rotZ = fParams->GetSaleveArmRotZ();
    OffsetX = fParams->GetSaleveArmOffset().X();
    OffsetY = fParams->GetSaleveArmOffset().Y();
    OffsetZ = fParams->GetSaleveArmOffset().Z();
  } else {
    Error("LocalToGlobal", "Wrong arm id, check it out");
  }

  // cout<<"iarm="<<iarm<<"  "<<rotX<<" "<<rotY<<" "<<rotZ<<" "<<OffsetX<<" "<<OffsetY<<" "<<OffsetZ<<endl;

  double sa = TMath::Sin(rotZ);
  double ca = TMath::Cos(rotZ);
  double sb = TMath::Sin(rotY);
  double cb = TMath::Cos(rotY);
  double sg = TMath::Sin(rotX);
  double cg = TMath::Cos(rotX);

  // apply rotations
  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = out->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
    if (!tracktab) continue;

    // cout<<"table: "<<tracktab->GetName()<<"  "<<tracktab->GetEntries()<<endl;
    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      track->Activate();

      fhCurvature_back->Fill(track->GetCurvature());

      double x0 = track->GetXatZ(0);
      double y0 = track->GetYatZ(0);
      double z0 = 0;
      double x1 = track->GetXatZ(100);
      double y1 = track->GetYatZ(100);
      double z1 = 100;
      // transform to arm middle plane
      z0 = z0 - 75.0;
      z1 = z1 - 75.0;

      // XY
      // double xx0 =      cb * x0                +     sb * z0;
      // double yy0 =  -sg*sb * x0   +   cg * y0  +  sg*cb * z0;
      // double zz0 =  -cg*sb * x0   -   sg * y0  +  cg*cb * z0;

      // double xx1 =    (cb) * x1                +     sb * z1;
      // double yy1 =  -sg*sb * x1   +   cg * y1  +  sg*cb * z1;
      // double zz1 =  -cg*sb * x1   -   sg * y1  +  cg*cb * z1;

      // ZXY
      double xx0 = (ca * cb - sa * sg * sb) * x0 + sa * cg * y0 + (ca * sb + sa * sg * cb) * z0;
      double yy0 = -(sa * cb + ca * sg * sb) * x0 + ca * cg * y0 + (-sa * sb + ca * sg * cb) * z0;
      double zz0 = -cg * sb * x0 - sg * y0 + cg * cb * z0;

      double xx1 = (ca * cb - sa * sg * sb) * x1 + sa * cg * y1 + (ca * sb + sa * sg * cb) * z1;
      double yy1 = -(sa * cb + ca * sg * sb) * x1 + ca * cg * y1 + (-sa * sb + ca * sg * cb) * z1;
      double zz1 = -cg * sb * x1 - sg * y1 + cg * cb * z1;

      zz0 = zz0 + 75.0;
      zz1 = zz1 + 75.0;
      Vector3D origin(xx0, yy0, zz0);
      Vector3D direction(xx1 - xx0, yy1 - yy0, zz1 - zz0);

      (track->Getline())->SetOrigin(origin);
      (track->Getline())->SetDirection(direction);
    }
  }

  // apply offsets
  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = out->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
    if (!tracktab) continue;

    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      // track -> Activate();
      double ox = (track->Getline())->GetOrigin().X();
      double oy = (track->Getline())->GetOrigin().Y();
      double oz = (track->Getline())->GetOrigin().Z();
      // cout<<"oz="<<oz<<endl;

      (track->Getline())->SetOrigin(ox + OffsetX, oy + OffsetY, oz + OffsetZ);
      // double oxn = (track->Getline())->GetOrigin().X();
      // double oyn = (track->Getline())->GetOrigin().Y();
      // double ozn = (track->Getline())->GetOrigin().Z();
      // cout<<"oxn="<<oxn<<"  oyn="<<oyn<<" ozn="<<ozn<<endl;
      double x0 = track->GetXatZ(0);
      double y0 = track->GetYatZ(0);
      // cout<<"x0="<<x0<<"  y0="<<y0<<" z="<<0<<endl;
      (track->Getline())->SetOrigin(x0, y0, 0.);
    }
  }
  // cout<<" ------------------------------ iarm="<<iarm<<endl;
  TransformHits(out, OffsetX, OffsetY, OffsetZ, sb, cb, sg, cg, sa, ca);
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::TransformHits(UEventNode* out, double dx, double dy, double dz, double sb, double cb, double sg, double cg, double sa, double ca) {
  UDataTable* hits = 0;

  for (int is = 0; is < 8; is++) {
    hits = out->GetDataTable(Form("Hits %s", fSensorNames[is].Data()));

    if (!hits) continue;

    for (int i = 0; i < hits->GetEntries(); i++) {
      USensorHit* hit = (USensorHit*)hits->At(i);
      double x0 = hit->GetX();
      double y0 = hit->GetY();
      double z0 = hit->GetZ();

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

      hit->SetX(x1 + dx);
      hit->SetY(y1 + dy);
      hit->SetZ(z1 + dz);
    }
  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::TransformHitsAl(UEventNode* out, double dx, double dy, double dz, double sb, double cb, double sg, double cg, double sa, double ca) {
  UDataTable* hits = 0;

  for (int is = 0; is < 34; is++) {
    hits = out->GetDataTable(Form("Hits %s", fAlSensorNames[is].Data()));

    if (!hits) continue;

    for (int i = 0; i < hits->GetEntries(); i++) {
      USensorHit* hit = (USensorHit*)hits->At(i);
      double x0 = hit->GetX();
      double y0 = hit->GetY();
      double z0 = hit->GetZ();

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

      hit->SetX(x1 + dx);
      hit->SetY(y1 + dy);
      hit->SetZ(z1 + dz);
    }
  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::FillClusterCorrelationAl(UEventNode* inJ, UEventNode* inS) {
  
  UDataTable* hitsJ = 0;
  UDataTable* hitsS = 0;

  //inJ->ListObjects();
  //cout<<"--------------------------------------:::"<<endl;
  //inS->ListObjects();

  for (int is = 0; is < 8; is++) {
    hitsJ = inJ->GetDataTable(Form("Hits %s", fAlSensorNames[is].Data()));
    hitsS = inS->GetDataTable(Form("Hits %s", fAlSensorNames[is].Data()));

    double mhj = 0;
    double mhs = 0;

    if (hitsJ)
      mhj = hitsJ->GetEntries(); 
    else
      cout << " Na61PrimaryVertexRecoModule::FillClusterCorrelation: No hits is Jura sensor " << fAlSensorNames[is].Data() << endl;

    if (hitsS)
      mhs = hitsS->GetEntries();
    else
      cout << " Na61PrimaryVertexRecoModule::FillClusterCorrelation: No hits is Saleve sensor " << fAlSensorNames[is].Data() << endl;

    //cout<<fhClusters_JvsS[is]<<"  "<<fhClusters_JmS[is]<<endl;

    fhClusters_JvsS[is]->Fill(mhj, mhs);
    fhClusters_JmS[is]->Fill(0.5 * (mhj + mhs), mhj - mhs);

  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::FillClusterCorrelation(UEventNode* inJ, UEventNode* inS) {
  
  UDataTable* hitsJ = 0;
  UDataTable* hitsS = 0;

  for (int is = 0; is < 8; is++) {
    hitsJ = inJ->GetDataTable(Form("Hits %s", fSensorNames[is].Data()));
    hitsS = inS->GetDataTable(Form("Hits %s", fSensorNames[is].Data()));

    double mhj = 0;
    double mhs = 0;

    if (hitsJ)
      mhj = hitsJ->GetEntries();
    else
      cout << " Na61PrimaryVertexRecoModule::FillClusterCorrelation: No hits is Jura sensor " << fSensorNames[is].Data() << endl;

    if (hitsS)
      mhs = hitsS->GetEntries();
    else
      cout << " Na61PrimaryVertexRecoModule::FillClusterCorrelation: No hits is Saleve sensor " << fSensorNames[is].Data() << endl;

    fhClusters_JvsS[is]->Fill(mhj, mhs);
    fhClusters_JmS[is]->Fill(0.5 * (mhj + mhs), mhj - mhs);
  }
}

//_______________________________________________________________
void Na61PrimaryVertexRecoModule::FillHitPositions(Int_t arm, UEventNode* out)
{

  for(Int_t it=0;it<34;it++){
    UDataTable* hits = out->GetDataTable(Form("Hits %s",fAlSensorNames[it].Data()));

    for(Int_t i=0;i<hits->GetEntries();i++){
      USensorHit* hit = (USensorHit*)hits->At(i);      

      Float_t x = hit->GetX();
      Float_t y = hit->GetY();
      Float_t z = hit->GetZ();
      
 
      fhZX_fine[arm]->Fill(z,x);
      fhXY_fine[arm]->Fill(x,y);
      if(fAlSensorNames[it].Contains("Al1_")){fhX_Al1[arm]->Fill(x); fhY_Al1[arm]->Fill(y);} 
      if(fAlSensorNames[it].Contains("Al2_")){fhX_Al2[arm]->Fill(x); fhY_Al2[arm]->Fill(y);} 
      if(fAlSensorNames[it].Contains("Al3_")){fhX_Al3[arm]->Fill(x); fhY_Al3[arm]->Fill(y);} 
      if(fAlSensorNames[it].Contains("Al4_")){fhX_Al4[arm]->Fill(x); fhY_Al4[arm]->Fill(y);} 


    }
  }

}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::ClassifyTracksAl(UEventNode* outJ, UEventNode* outS) {
  // clasification regarding x-slope
  // flag definition:
  // flag0 - beam particles see as a narrow peak
  // flag1 - upsteam interactions - cases located at low ax and ay
  // flag2 - target interaction particles 
  // clasification should based on visual inspection of the ax and ay distributions 

  for(Int_t i=0;i<fNDevComb;i++){
    UDataTable* tracktab = outJ -> GetDataTable(Form("FullTracks %s",fDevCombName[i].Data()));
    if (!tracktab) continue;
    
    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      double ax = track->GetDX() / track->GetDZ();
      double ay = track->GetDY() / track->GetDZ();
      
      
      track->SetFlag(2);
      fhAxJ_flag2->Fill(ax);
      fhAyJ_flag2->Fill(ay);
      fhAxJ->Fill(ax);
      fhAyJ->Fill(ay);
    }

  }


  for(Int_t i=0;i<fNDevComb;i++){
    UDataTable* tracktab = outS -> GetDataTable(Form("FullTracks %s",fDevCombName[i].Data()));
    if (!tracktab) continue;
    
    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      
      double ax = track->GetDX() / track->GetDZ();
      double ay = track->GetDY() / track->GetDZ();
      
      track->SetFlag(2);
      fhAxS_flag2->Fill(ax);
      fhAyS_flag2->Fill(ay);
      
      fhAxS->Fill(ax);
      fhAyS->Fill(ay);

    }

  }

}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::ClassifyTracks(UEventNode* outJ, UEventNode* outS) {
  // clasification regarding x-slope

  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = outJ->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
    if (!tracktab) continue;

    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      double ax = track->GetDX() / track->GetDZ();
      double ay = track->GetDY() / track->GetDZ();

      if (i > 1) {
        track->SetFlag(2);
        fhAxJ_flag2->Fill(ax);
        fhAyJ_flag2->Fill(ay);
      } else {
        if (ax < (0.006 + fRotY_J)) {
          track->SetFlag(0);
          fhAxJ_flag0->Fill(ax);
          fhAyJ_flag0->Fill(ay);
        } else if (ax < (fOffAx_J + fRotY_J)) {
          track->SetFlag(1);
          fhAxJ_flag1->Fill(ax);
        } else {
          track->SetFlag(2);
          fhAxJ_flag2->Fill(ax);
          fhAyJ_flag2->Fill(ay);
        }
      }

      if (track->GetFlag() != 2) continue;
      fhAxJ->Fill(ax);
      fhAyJ->Fill(ay);
    }
  }

  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = outS->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
    if (!tracktab) continue;

    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);

      double ax = track->GetDX() / track->GetDZ();
      double ay = track->GetDY() / track->GetDZ();

      if (i > 1) {
        track->SetFlag(2);
        fhAxS_flag2->Fill(ax);
        fhAyS_flag2->Fill(ay);
      } else {
        if (ax > (0.005 + fRotY_S)) {
          track->SetFlag(0);
          fhAxS_flag0->Fill(ax);
          fhAyS_flag0->Fill(ay);
        } else if (ax > (fOffAx_S + fRotY_S)) {  //
          track->SetFlag(1);
          fhAxS_flag1->Fill(ax);
        } else {
          track->SetFlag(2);
          fhAxS_flag2->Fill(ax);
          fhAyS_flag2->Fill(ay);
        }
      }

      if (track->GetFlag() != 2) continue;
      fhAxS->Fill(ax);
      fhAyS->Fill(ay);
    }
  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::ClassifyTracks_pPb(UEventNode* outJ, UEventNode* outS) {
  // clasification regarding x-slope and y-slope

  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = outJ->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
    if (!tracktab) continue;

    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      double ax = track->GetDX() / track->GetDZ();
      double ay = track->GetDY() / track->GetDZ();

      if (i > 1) {
        fJuraProdCounter++;
        track->SetFlag(2);
        fhAxJ->Fill(ax);
        fhAyJ->Fill(ay);
        continue;
      }

      fhAxJ_flag02->Fill(ax);
      fhAyJ_flag02->Fill(ay);
      double dax = ax - fOffAx_J;
      double day = ay - fOffAy_J;
      double dd = (dax * dax) / (fSigAx_J * fSigAx_J) + (day * day) / (fSigAy_J * fSigAy_J);

      if (dd < fNsigJ * fNsigJ) {
        track->SetFlag(0);
        fhAxJ_flag0->Fill(ax);
        fhAyJ_flag0->Fill(ay);
      } else {
        fJuraProdCounter++;
        track->SetFlag(2);
        fhAxJ->Fill(ax);
        fhAyJ->Fill(ay);
        fhAxJ_flag2->Fill(ax);
        fhAyJ_flag2->Fill(ay);
      }
    }
  }

  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = outS->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
    if (!tracktab) continue;

    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);

      double ax = track->GetDX() / track->GetDZ();
      double ay = track->GetDY() / track->GetDZ();

      if (i > 1) {
        fSaleveProdCounter++;
        track->SetFlag(2);
        fhAxS->Fill(ax);
        fhAyS->Fill(ay);
        continue;
      }

      fhAxS_flag02->Fill(ax);
      fhAyS_flag02->Fill(ay);
      double dax = ax - fOffAx_S;
      double day = ay - fOffAy_S;
      double dd = (dax * dax) / (fSigAx_S * fSigAx_S) + (day * day) / (fSigAy_S * fSigAy_S);

      if (dd < fNsigS * fNsigS) {
        track->SetFlag(0);
        fhAxS_flag0->Fill(ax);
        fhAyS_flag0->Fill(ay);
      } else {
        fSaleveProdCounter++;
        track->SetFlag(2);
        fhAxS->Fill(ax);
        fhAyS->Fill(ay);
        fhAxS_flag2->Fill(ax);
        fhAyS_flag2->Fill(ay);
      }
    }
  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::FindPrimaryVertexAl(int flag, UEventNode* outJ, UEventNode* outS) {
  // algorithm:

  int fullTracksJ = 0;
  int fullTracksS = 0;
  int fullTracksJ4h = 0;
  int fullTracksS4h = 0;
  TObjArray tracks;
  tracks.Clear();

  int histId = 0;         // Combined primary vertex
  if (!outS) histId = 1;  // Jura primary vertex
  if (!outJ) histId = 2;  // Saleve primary vertex

  int tabN = fNDevComb;

  int jura_tracks = 0;
  if (outJ) {
    for (int i = 0; i < tabN; i++) {
      UDataTable* tracktab = outJ->GetDataTable(Form("FullTracks %s", fDevCombName[i].Data()));
      if (!tracktab) continue;
      
      fullTracksJ = fullTracksJ + tracktab->GetEntries();
      if (i < 4) fullTracksJ4h = fullTracksJ4h + tracktab->GetEntries();
      
      for (int i = 0; i < tracktab->GetEntries(); i++) {
        UVdTrack* track = (UVdTrack*)tracktab->At(i);
        if (track->GetFlag() == flag) {
          tracks.Add(tracktab->At(i));
          jura_tracks++;
        }
      }
    }
  }

  int saleve_tracks = 0;
  if (outS) {
    for (int i = 0; i < tabN; i++) {
      UDataTable* tracktab = outS->GetDataTable(Form("FullTracks %s", fDevCombName[i].Data()));
      if (!tracktab) continue;

      fullTracksS = fullTracksS + tracktab->GetEntries();
      if (i < 4) fullTracksS4h = fullTracksS4h + tracktab->GetEntries();

      for (int i = 0; i < tracktab->GetEntries(); i++) {
        UVdTrack* track = (UVdTrack*)tracktab->At(i);
        if (track->GetFlag() == flag) {
          tracks.Add(tracktab->At(i));
          saleve_tracks++;
        }
      }
    }
  }

  if (flag == 2 && histId == 0) {
    fhFullTracks_JvsS->Fill(fullTracksJ, fullTracksS);
    fhFullTracks_JS->Fill(fullTracksJ + fullTracksS);
    fhFullTracks_J->Fill(fullTracksJ);
    fhFullTracks_S->Fill(fullTracksS);
    fhFullTracks4h_JS->Fill(fullTracksJ4h + fullTracksS4h);
    fhFullTracks4h_J->Fill(fullTracksJ4h);
    fhFullTracks4h_S->Fill(fullTracksS4h);

    fAllTracksJ = fullTracksJ;
    fAllTracksS = fullTracksS;

    // cout<<"-----------> 4 hit tracks: "<< fullTracksJ4h+fullTracksS4h <<endl;
  }

  // cout<<"--------------->1 histId="<<histId<<"  flag="<<flag<<"  "<<jura_tracks<<" "<<saleve_tracks<<endl;

  if (tracks.GetEntries() < 2) return;

  if ((outJ && outS)) {
    if ((jura_tracks + saleve_tracks) < 2) return;
    // if(jura_tracks < 2)   return;
    // if(saleve_tracks < 2) return;
  }

  double sum_up_x = 0;
  double sum_do_x = 0;
  double sum_up_y = 0;
  double sum_do_y = 0;

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      double axi = track1->GetDX() / track1->GetDZ();
      double axj = track2->GetDX() / track2->GetDZ();
      double bxi = track1->GetX();
      double bxj = track2->GetX();

      double ayi = track1->GetDY() / track1->GetDZ();
      double ayj = track2->GetDY() / track2->GetDZ();
      double byi = track1->GetY();
      double byj = track2->GetY();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < 3 || flag != 2) {
        sum_up_x = sum_up_x + (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - fZprim) < 3 || flag != 2) {
        sum_up_y = sum_up_y + (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + (ayi - ayj) * (ayi - ayj);
      }
    }
  }

  ////////// new method //////////////////////////////////////////////////

  //  double zx_min = 0;
  //  if(sum_do_x > 0)zx_min = - sum_up_x/sum_do_x;
  //  double zy_min = 0;
  //  if(sum_do_y > 0)zy_min = - sum_up_y/sum_do_y;

  // cout<<sum_do_x<<"  "<<sum_do_y<<"  "<<zx_min<<" "<<zy_min<<endl;

  double z_prim = -(sum_up_x + sum_up_y) / (sum_do_x + sum_do_y);

  Vector3D reco_pvertex_w2;
  PrimaryVertexWithWeigths(flag, tracks, z_prim, reco_pvertex_w2);
  if (flag == 2) {
    fhRecoVertexXY[histId]->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
    fhRecoVertexZ[histId]->Fill(reco_pvertex_w2.Z());
    fhRecoVertexZ_fine[histId]->Fill(reco_pvertex_w2.Z());
    if (histId == 0) {
      fhVtxX->Fill(reco_pvertex_w2.X());
      fhVtxY->Fill(reco_pvertex_w2.Y());
      fhTracksVsVtxX_J->Fill(reco_pvertex_w2.X(), fullTracksJ);
      fhTracksVsVtxX_S->Fill(reco_pvertex_w2.X(), fullTracksS);
      fhTracksVsVtxY_J->Fill(reco_pvertex_w2.Y(), fullTracksJ);
      fhTracksVsVtxY_S->Fill(reco_pvertex_w2.Y(), fullTracksS);

      fPrimaryVertex.SetX(reco_pvertex_w2.X());
      fPrimaryVertex.SetY(reco_pvertex_w2.Y());
      fPrimaryVertex.SetZ(reco_pvertex_w2.Z());
    }
  }
  if (flag == 1) {
    fhRecoVertexXY_flag1[histId]->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
    fhRecoVertexZ_flag1[histId]->Fill(reco_pvertex_w2.Z());
    fhRecoVertexZ_fine_flag1[histId]->Fill(reco_pvertex_w2.Z());
    // fPrimaryVertexDefined = true;
    //((UVdEvent*)out)->SetPrimaryVertexStatus(1);
  }
  // setup arm vertex
  if (outJ && (!outS)) ((UVdEvent*)outJ)->SetPrimaryVertex(reco_pvertex_w2.X(), reco_pvertex_w2.Y(), reco_pvertex_w2.Z());
  if (outS && (!outJ)) ((UVdEvent*)outS)->SetPrimaryVertex(reco_pvertex_w2.X(), reco_pvertex_w2.Y(), reco_pvertex_w2.Z());
}

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::FindPrimaryVertex(int flag, UEventNode* outJ, UEventNode* outS) {
  // algorithm:

  int fullTracksJ = 0;
  int fullTracksS = 0;
  int fullTracksJ4h = 0;
  int fullTracksS4h = 0;
  TObjArray tracks;
  tracks.Clear();

  int histId = 0;         // Combined primary vertex
  if (!outS) histId = 1;  // Jura primary vertex
  if (!outJ) histId = 2;  // Saleve primary vertex

  int tabN = 18;
  // int tabN = 4;

  int jura_tracks = 0;
  if (outJ) {
    for (int i = 0; i < tabN; i++) {
      UDataTable* tracktab = outJ->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
      if (!tracktab) continue;

      fullTracksJ = fullTracksJ + tracktab->GetEntries();
      if (i < 4) fullTracksJ4h = fullTracksJ4h + tracktab->GetEntries();

      for (int i = 0; i < tracktab->GetEntries(); i++) {
        UVdTrack* track = (UVdTrack*)tracktab->At(i);
        if (track->GetFlag() == flag) {
          tracks.Add(tracktab->At(i));
          jura_tracks++;
        }
      }
    }
  }

  int saleve_tracks = 0;
  if (outS) {
    for (int i = 0; i < tabN; i++) {
      UDataTable* tracktab = outS->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
      if (!tracktab) continue;

      fullTracksS = fullTracksS + tracktab->GetEntries();
      if (i < 4) fullTracksS4h = fullTracksS4h + tracktab->GetEntries();

      for (int i = 0; i < tracktab->GetEntries(); i++) {
        UVdTrack* track = (UVdTrack*)tracktab->At(i);
        if (track->GetFlag() == flag) {
          tracks.Add(tracktab->At(i));
          saleve_tracks++;
        }
      }
    }
  }

  if (flag == 2 && histId == 0) {
    fhFullTracks_JvsS->Fill(fullTracksJ, fullTracksS);
    fhFullTracks_JS->Fill(fullTracksJ + fullTracksS);
    fhFullTracks_J->Fill(fullTracksJ);
    fhFullTracks_S->Fill(fullTracksS);
    fhFullTracks4h_JS->Fill(fullTracksJ4h + fullTracksS4h);
    fhFullTracks4h_J->Fill(fullTracksJ4h);
    fhFullTracks4h_S->Fill(fullTracksS4h);

    fAllTracksJ = fullTracksJ;
    fAllTracksS = fullTracksS;

    // cout<<"-----------> 4 hit tracks: "<< fullTracksJ4h+fullTracksS4h <<endl;
  }

  // cout<<"--------------->1 histId="<<histId<<"  flag="<<flag<<"  "<<jura_tracks<<" "<<saleve_tracks<<endl;

  if (tracks.GetEntries() < 2) return;

  if ((outJ && outS)) {
    if ((jura_tracks + saleve_tracks) < 2) return;
    // if(jura_tracks < 2)   return;
    // if(saleve_tracks < 2) return;
  }

  double sum_up_x = 0;
  double sum_do_x = 0;
  double sum_up_y = 0;
  double sum_do_y = 0;

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      double axi = track1->GetDX() / track1->GetDZ();
      double axj = track2->GetDX() / track2->GetDZ();
      double bxi = track1->GetX();
      double bxj = track2->GetX();

      double ayi = track1->GetDY() / track1->GetDZ();
      double ayj = track2->GetDY() / track2->GetDZ();
      double byi = track1->GetY();
      double byj = track2->GetY();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < 3 || flag != 2) {
        sum_up_x = sum_up_x + (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - fZprim) < 3 || flag != 2) {
        sum_up_y = sum_up_y + (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + (ayi - ayj) * (ayi - ayj);
      }
    }
  }

  ////////// new method //////////////////////////////////////////////////

  //  double zx_min = 0;
  //  if(sum_do_x > 0)zx_min = - sum_up_x/sum_do_x;
  //  double zy_min = 0;
  //  if(sum_do_y > 0)zy_min = - sum_up_y/sum_do_y;

  // cout<<sum_do_x<<"  "<<sum_do_y<<"  "<<zx_min<<" "<<zy_min<<endl;

  double z_prim = -(sum_up_x + sum_up_y) / (sum_do_x + sum_do_y);

  Vector3D reco_pvertex_w2;
  PrimaryVertexWithWeigths(flag, tracks, z_prim, reco_pvertex_w2);
  if (flag == 2) {
    fhRecoVertexXY[histId]->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
    fhRecoVertexZ[histId]->Fill(reco_pvertex_w2.Z());
    fhRecoVertexZ_fine[histId]->Fill(reco_pvertex_w2.Z());
    if (histId == 0) {
      fhVtxX->Fill(reco_pvertex_w2.X());
      fhVtxY->Fill(reco_pvertex_w2.Y());
      fhTracksVsVtxX_J->Fill(reco_pvertex_w2.X(), fullTracksJ);
      fhTracksVsVtxX_S->Fill(reco_pvertex_w2.X(), fullTracksS);
      fhTracksVsVtxY_J->Fill(reco_pvertex_w2.Y(), fullTracksJ);
      fhTracksVsVtxY_S->Fill(reco_pvertex_w2.Y(), fullTracksS);

      fPrimaryVertex.SetX(reco_pvertex_w2.X());
      fPrimaryVertex.SetY(reco_pvertex_w2.Y());
      fPrimaryVertex.SetZ(reco_pvertex_w2.Z());
    }
  }
  if (flag == 1) {
    fhRecoVertexXY_flag1[histId]->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
    fhRecoVertexZ_flag1[histId]->Fill(reco_pvertex_w2.Z());
    fhRecoVertexZ_fine_flag1[histId]->Fill(reco_pvertex_w2.Z());
    // fPrimaryVertexDefined = true;
    //((UVdEvent*)out)->SetPrimaryVertexStatus(2);
  }
  // setup arm vertex
  if (outJ && (!outS)) ((UVdEvent*)outJ)->SetPrimaryVertex(reco_pvertex_w2.X(), reco_pvertex_w2.Y(), reco_pvertex_w2.Z());
  if (outS && (!outJ)) ((UVdEvent*)outS)->SetPrimaryVertex(reco_pvertex_w2.X(), reco_pvertex_w2.Y(), reco_pvertex_w2.Z());
}

//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoModule::PrimaryVertexWithWeigths(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex) {
  // algorithm:
  // 1. Calculate weights = 1/dist^2, distance betweem i-th and j-th tracks on
  // the closest proximity plane (z=zprim).
  // 2. Used the weights to calculate z coordinate of primary vertex from the formula that inclused weights (see PS notes).
  // 3. Find x,y coordinates averageing tracks x, y at the closes proximmity plane.

  double sum_up_x = 0;
  double sum_do_x = 0;
  double sum_up_y = 0;
  double sum_do_y = 0;

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      double x1 = track1->GetXatZ(zprim);
      double x2 = track2->GetXatZ(zprim);
      double y1 = track1->GetYatZ(zprim);
      double y2 = track2->GetYatZ(zprim);

      double d2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
      double w = 1. / d2;
      // fhd2->Fill(d2);

      /////////////////// stuff for new method
      double axi = track1->GetDX() / track1->GetDZ();
      double axj = track2->GetDX() / track2->GetDZ();
      double bxi = track1->GetX();
      double bxj = track2->GetX();

      double ayi = track1->GetDY() / track1->GetDZ();
      double ayj = track2->GetDY() / track2->GetDZ();
      double byi = track1->GetY();
      double byj = track2->GetY();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < 3 || flag != 2) {
        sum_up_x = sum_up_x + w * (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + w * (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - fZprim) < 3 || flag != 2) {
        sum_up_y = sum_up_y + w * (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + w * (ayi - ayj) * (ayi - ayj);
      }
      // cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
    }
  }

  double zprim_w2 = -(sum_up_x + sum_up_y) / (sum_do_x + sum_do_y);

  //
  double sum_x = 0;
  double sum_y = 0;
  double N = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);

    double x = track->GetXatZ(zprim_w2);
    double y = track->GetYatZ(zprim_w2);

    sum_x = sum_x + x;
    sum_y = sum_y + y;
    N = N + 1;
  }

  double xprim = sum_x / N;
  double yprim = sum_y / N;

  sum_x = 0;
  sum_y = 0;
  double norm = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);

    double x = track->GetXatZ(zprim_w2);
    double y = track->GetYatZ(zprim_w2);

    double d2 = (x - xprim) * (x - xprim) + (y - yprim) * (y - yprim);

    sum_x = sum_x + x / d2;
    sum_y = sum_y + y / d2;
    norm = norm + 1. / d2;
  }

  double xprim_w2 = sum_x / norm;
  double yprim_w2 = sum_y / norm;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);
}

//////////////////////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________
void Na61PrimaryVertexRecoModule::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61PrimaryVertexRecoModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);
}

//____________________________________________________________________
void Na61PrimaryVertexRecoModule::Print(Option_t* option) const {
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>

  TString opt(option);
  opt.ToLower();

  Na61Module::Print(option);
  if (opt.Contains("d")) cout << endl << "  Original author: Paweł Staszel" << endl << "  Last Modifications: " << endl << "    $Author: Staszel $" << endl << "    $Date: 2017/03/20$" << endl << "    $Revision: 1.0 $ " << endl << endl << "-------------------------------------------------" << endl;
}
