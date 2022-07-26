//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61VdTrackingModule
#include "Na61VdTrackingModule.h"
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
// ClassImp(Na61VdTrackingModule);

//____________________________________________________________________
Na61VdTrackingModule::Na61VdTrackingModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
  fProductionMode = 1;
}

//____________________________________________________________________
Na61VdTrackingModule::Na61VdTrackingModule(const char* name, const char* title) : Na61Module(name, title) {
  // Named Constructor
  SetState(kSetup);
  fProductionMode = 1;

  fNsig_dev = 4;    // n sigma cut for 3HitTracks
  fNsig_pvert = 5;  // n sigma cut for matchig of 3HitTracks with primary vertex

  fPrimaryVertex = new Vector3D();

  TString gaus = "[0]*exp(-0.5*((x-[1])/[2])^2)";
  fFGauss = new TF1("fFGauss", gaus.Data(), -500, 100);
  fFGauss->SetParNames("const", "mean", "sigma");

  for (int i = 0; i < 8; i++) fHitsTabs[i] = 0;
}
//____________________________________________________________________
void Na61VdTrackingModule::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  cout << " Na61VdTrackingModule::DefineHistograms: " << fHistDirName.Data() << endl;
  TDirectory* histDir = gDirectory->mkdir(fHistDirName.Data());
  histDir->cd();

  fhCurvature_front = new TH1F("hCurvature_front", "", 500, -0.00003, 0.00003);
  fhCurvature_back = new TH1F("hCurvature_back", "", 500, -0.00003, 0.00003);
  fhSlopeChange = new TH1F("hSlopeChange", "", 500, -0.008, 0.008);

  for (int i = 0; i < 4; i++) {
    fdx3hit[i] = new TH1F(Form("hdx3hit_Vds%d", i + 1), "", 500, -0.2, 0.2);
    fdx4hit_front[i] = new TH1F(Form("hdx4hit_front_Vds%d", i + 1), "", 500, -0.2, 0.2);
    fdx4hit_back[i] = new TH1F(Form("hdx4hit_back_Vds%d", i + 1), "", 500, -0.2, 0.2);
  }

  if (fHistDirName.Contains("HT")) {  // VdTrackingHT dooes not need rest of histos
    gDirectory->cd("..");
    return;
  }

  int nbins = 5000;
  double xmin = -5.0;
  double xmax = 5.0;

  for (int i = 0; i < 6; i++) {
    TString devname = fDevNames[i];
    fhX_dev[i] = new TH1F(Form("hX_%s", devname.Data()), "", nbins, xmin, xmax);
    fhY_dev[i] = new TH1F(Form("hY_%s", devname.Data()), "", nbins, xmin, xmax);
    fhX_dev_cuty[i] = new TH1F(Form("hX_%s_cuty", devname.Data()), "", nbins, xmin, xmax);
    fhY_dev_cutx[i] = new TH1F(Form("hY_%s_cutx", devname.Data()), "", nbins, xmin, xmax);

    fh_devx_devy[i] = new TH2F(Form("h_devx_devy_%s", devname.Data()), "", 100, -0.1, 0.1, 100, -0.1, 0.1);

    fhX_dev_lc[i] = new TH1F(Form("hX_%s_lc", devname.Data()), "", nbins, xmin, xmax);
    fhY_dev_lc[i] = new TH1F(Form("hY_%s_lc", devname.Data()), "", nbins, xmin, xmax);
    fhX_dev_cuty_lc[i] = new TH1F(Form("hX_%s_cuty_lc", devname.Data()), "", nbins, xmin, xmax);
    fhY_dev_cutx_lc[i] = new TH1F(Form("hY_%s_cutx_lc", devname.Data()), "", nbins, xmin, xmax);

    fhAx[i] = new TH1F(Form("hAx_%s", devname.Data()), "", 2500, -0.25, 0.25);
    fhAy[i] = new TH1F(Form("hAy_%s", devname.Data()), "", 2500, -0.25, 0.25);

    fhAx_lc[i] = new TH1F(Form("hAx_%s_lc", devname.Data()), "", 2500, -0.25, 0.25);
    fhAy_lc[i] = new TH1F(Form("hAy_%s_lc", devname.Data()), "", 2500, -0.25, 0.25);

    fhChi2X[i] = new TH1F(Form("hChi2X_%s", devname.Data()), "", 200, 0., 100.);
    fhChi2Y[i] = new TH1F(Form("hChi2Y_%s", devname.Data()), "", 200, 0., 100.);
    // fhAx_lc[i] = 0;
    // fhAy_lc[i] = 0;

    fhNsVsNe[i] = new TH2F(Form("hNsVsNe_%s", devname.Data()), "", 100, 0, 100, 100, 0, 100);

    fhTracks[i] = new TH1F(Form("hTracks_%s", devname.Data()), "", 100, 0., 100.);

    fhClusterSize[i] = new TH1F(Form("hClusterSize_%s", devname.Data()), "", 350, 0, 350);
    fhClusterSizeInPick[i] = new TH1F(Form("hClusterSizeInPick_%s", devname.Data()), "", 350, 0, 350);
  }

  histDir = gDirectory->mkdir("3hitx");
  histDir->cd();

  for (int i = 0; i < 14; i++) {
    TString devname = fDevNames[i + 6];
    fhX_dev[i + 6] = new TH1F(Form("hX_%s", devname.Data()), "", nbins, xmin, xmax);
    fhY_dev[i + 6] = new TH1F(Form("hY_%s", devname.Data()), "", nbins, xmin, xmax);
    fhX_dev_cuty[i + 6] = new TH1F(Form("hX_%s_cuty", devname.Data()), " ", nbins, xmin, xmax);
    fhY_dev_cutx[i + 6] = new TH1F(Form("hY_%s_cutx", devname.Data()), "", nbins, xmin, xmax);
    fhX_dev_acce[i + 6] = new TH1F(Form("hX_%s_acce", devname.Data()), " ", nbins, xmin, xmax);
    fhY_dev_acce[i + 6] = new TH1F(Form("hY_%s_acce", devname.Data()), "", nbins, xmin, xmax);

    fh_devx_devy[i + 6] = new TH2F(Form("h_devx_devy_%s", devname.Data()), "", 100, -0.1, 0.1, 100, -0.1, 0.1);

    // cout<<"i+6="<<i+6<<" hist name="<<fh_devx_devy[i+6]->GetName()<<endl;

    fhAx[i + 6] = new TH1F(Form("hAx_%s", devname.Data()), "", 2500, -0.25, 0.25);
    fhAy[i + 6] = new TH1F(Form("hAy_%s", devname.Data()), "", 2500, -0.25, 0.25);

    fhChi2X[i + 6] = new TH1F(Form("hChi2X_%s", devname.Data()), "", 200, 0., 100.);
    fhChi2Y[i + 6] = new TH1F(Form("hChi2Y_%s", devname.Data()), "", 200, 0., 100.);

    fhNsVsNe[i + 6] = new TH2F(Form("hNsVsNe_%s", devname.Data()), "", 100, 0, 100, 100, 0, 100);

    fhTracks[i + 6] = new TH1F(Form("hTracks_%s", devname.Data()), "", 100, 0., 100.);
  }

  gDirectory->cd("..");

  fhFullTracks = new TH1F(Form("hFullTracks"), "", 1000, 0., 1000.);
  fhFullTracksAll = new TH1F(Form("hFullTracksAll"), "", 1000, 0., 1000.);

  // cout<<"1  fhFullTracksAll: "<<fhFullTracksAll<<endl;

  for (int i = 0; i < 4; i++) {
    TString matchstr = fMatchStr[i];

    fhAx_full[i] = new TH1F(Form("hAx_full_%s", matchstr.Data()), "", 2500, -0.25, 0.25);
    fhAy_full[i] = new TH1F(Form("hAy_full_%s", matchstr.Data()), "", 2500, -0.25, 0.25);
    fhAy_full_prod[i] = new TH1F(Form("hAy_full_prod_%s", matchstr.Data()), "", 2500, -0.25, 0.25);

    fhX_dev_full[i] = new TH1F(Form("hX_dev_full_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhY_dev_full[i] = new TH1F(Form("hY_dev_full_%s", matchstr.Data()), "", 1000, -0.25, 0.25);

    fhX_dev_Vds1[i] = new TH1F(Form("hX_dev_Vds1_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhX_dev_Vds2[i] = new TH1F(Form("hX_dev_Vds2_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhX_dev_Vds3[i] = new TH1F(Form("hX_dev_Vds3_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhX_dev_Vds4[i] = new TH1F(Form("hX_dev_Vds4_%s", matchstr.Data()), "", 1000, -0.25, 0.25);

    fhY_dev_Vds1[i] = new TH1F(Form("hY_dev_Vds1_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhY_dev_Vds2[i] = new TH1F(Form("hY_dev_Vds2_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhY_dev_Vds3[i] = new TH1F(Form("hY_dev_Vds3_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhY_dev_Vds4[i] = new TH1F(Form("hY_dev_Vds4_%s", matchstr.Data()), "", 1000, -0.25, 0.25);

    fhDx_match[i] = new TH1F(Form("hDx_match_%s", matchstr.Data()), "", 5000, -10, 10);
    fhDy_match[i] = new TH1F(Form("hDy_match_%s", matchstr.Data()), "", 5000, -10, 10);
    fhDax_match[i] = new TH1F(Form("hDax_match_%s", matchstr.Data()), "", 5000, -0.5, 0.5);
    fhDay_match[i] = new TH1F(Form("hDay_match_%s", matchstr.Data()), "", 5000, -0.5, 0.5);
    fhDx_match_cut1[i] = new TH1F(Form("hDx_match_cut1_%s", matchstr.Data()), "", 5000, -10, 10);
    fhDy_match_cut1[i] = new TH1F(Form("hDy_match_cut1_%s", matchstr.Data()), "", 5000, -10, 10);
    fhDax_match_cut1[i] = new TH1F(Form("hDax_match_cut1_%s", matchstr.Data()), "", 5000, -0.5, 0.5);
    fhDay_match_cut1[i] = new TH1F(Form("hDay_match_cut1_%s", matchstr.Data()), "", 5000, -0.5, 0.5);
    fhDx_match_acce[i] = new TH1F(Form("hDx_match_%s_acce", matchstr.Data()), "", 1000, -2, 2);
    fhDy_match_acce[i] = new TH1F(Form("hDy_match_%s_acce", matchstr.Data()), "", 1000, -2, 2);
    fhDax_match_acce[i] = new TH1F(Form("hDax_match_%s_acce", matchstr.Data()), "", 1000, -0.05, 0.05);
    fhDay_match_acce[i] = new TH1F(Form("hDay_match_%s_acce", matchstr.Data()), "", 1000, -0.05, 0.05);

    cout << "i=" << i << " hist name=" << fhAx_full[i]->GetName() << endl;
  }

  histDir = gDirectory->mkdir("fullx");
  histDir->cd();

  for (int i = 0; i < 14; i++) {
    TString matchstr = fMatchStr[i + 4];

    fhAx_full[i + 4] = new TH1F(Form("hAx_full_%s", matchstr.Data()), "", 2500, -0.25, 0.25);
    fhAy_full[i + 4] = new TH1F(Form("hAy_full_%s", matchstr.Data()), "", 2500, -0.25, 0.25);
    fhAy_full_prod[i + 4] = new TH1F(Form("hAy_full_prod_%s", matchstr.Data()), "", 2500, -0.25, 0.25);

    cout << "i+4=" << i + 4 << " hist name=" << fhAx_full[i + 4]->GetName() << endl;

    fhX_dev_full[i + 4] = new TH1F(Form("hX_dev_full_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhY_dev_full[i + 4] = new TH1F(Form("hY_dev_full_%s", matchstr.Data()), "", 1000, -0.25, 0.25);

    fhDAx[i + 4] = new TH1F(Form("hDAx_%s", matchstr.Data()), " ", 2400, -0.06, 0.06);
    fhDAy[i + 4] = new TH1F(Form("hDAy_%s", matchstr.Data()), " ", 2400, -0.06, 0.06);
    fhDAx_cuty[i + 4] = new TH1F(Form("hDAx_%s_cuty", matchstr.Data()), " ", 2400, -0.06, 0.06);
    fhDAy_cutx[i + 4] = new TH1F(Form("hDAy_%s_cutx", matchstr.Data()), " ", 2400, -0.06, 0.06);

    fhDx[i + 4] = new TH1F(Form("hDx_%s", matchstr.Data()), " ", 2400, -1.0, 1.0);
    fhDy[i + 4] = new TH1F(Form("hDy_%s", matchstr.Data()), " ", 2400, -1.0, 1.0);
    fhDx_acce[i + 4] = new TH1F(Form("hDx_%s_acce", matchstr.Data()), " ", 2400, -1.0, 1.0);
    fhDy_acce[i + 4] = new TH1F(Form("hDy_%s_acce", matchstr.Data()), " ", 2400, -1.0, 1.0);

    fhX_dev_Vds1[i + 4] = new TH1F(Form("hX_dev_Vds1_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhX_dev_Vds2[i + 4] = new TH1F(Form("hX_dev_Vds2_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhX_dev_Vds3[i + 4] = new TH1F(Form("hX_dev_Vds3_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhX_dev_Vds4[i + 4] = new TH1F(Form("hX_dev_Vds4_%s", matchstr.Data()), "", 1000, -0.25, 0.25);

    fhY_dev_Vds1[i + 4] = new TH1F(Form("hY_dev_Vds1_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhY_dev_Vds2[i + 4] = new TH1F(Form("hY_dev_Vds2_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhY_dev_Vds3[i + 4] = new TH1F(Form("hY_dev_Vds3_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
    fhY_dev_Vds4[i + 4] = new TH1F(Form("hY_dev_Vds4_%s", matchstr.Data()), "", 1000, -0.25, 0.25);
  }

  gDirectory->cd("..");

  for (int i = 0; i < 6; i++) {
    TString devname = fDevNames[i];

    fhX_dev_acce[i] = new TH1F(Form("hX_%s_acce", devname.Data()), " ", 1000, -1., 1.);
    fhY_dev_acce[i] = new TH1F(Form("hY_%s_acce", devname.Data()), " ", 1000, -1., 1.);
  }

  fhChi2Ndf_x = new TH1F("hChi2Ndf_x", " ", 200, 0., 100.);
  fhChi2Ndf_y = new TH1F("hChi2Ndf_y", " ", 200, 0., 100.);

  // 4 hit tracks without weigths
  fhRecoVertexZ = new TH1F("hRecoVertexZ", " ", 2000, -500., 100.);
  fhRecoVertexZ_fine = new TH1F("hRecoVertexZ_fine", " ", 1000, -70., -30.);
  fhRecoVertexXY = new TH2F("hRecoVertexXY", " ", 500, -50., 50., 500, -50., 50.);

  fhRecoVertexZ_fine_x = new TH1F("hRecoVertexZ_fine_x", " ", 1000, -70., -30.);
  fhRecoVertexZ_fine_y = new TH1F("hRecoVertexZ_fine_y", " ", 1000, -70., -30.);

  fhd2 = new TH1F("hd2", "", 500, 0., 100.);

  // 4 hit tracks with weigths
  fhRecoVertexZ_w2 = new TH1F("hRecoVertexZ_w2", " ", 2000, -500., 100.);
  fhRecoVertexZ_fine_w2 = new TH1F("hRecoVertexZ_fine_w2", " ", 1000, -70., -30.);
  fhRecoVertexXY_w2 = new TH2F("hRecoVertexXY_w2", " ", 500, -50., 50., 500, -50., 50.);

  // All tracks with weigths
  fhRecoVertexZ_w2_Post = new TH1F("hRecoVertexZ_w2_Post", " ", 2000, -500., 100.);
  fhRecoVertexZ_fine_w2_Post = new TH1F("hRecoVertexZ_fine_w2_Post", " ", 1000, -70., -30.);
  fhRecoVertexXY_w2_Post = new TH2F("hRecoVertexXY_w2_Post", " ", 500, -50., 50., 500, -50., 50.);

  fhDZmin = new TH1F("hDZmin", " ", 2000, -100., 100.);

  gDirectory->cd("..");
}

//____________________________________________________________________
void Na61VdTrackingModule::Init() {
  // Job-level initialisation
  SetState(kInit);

  if (fRunId < 182) {
    fFieldRun = true;
  } else if (fRunId < 190) {
    fFieldRun = false;
  } else if (fRunId < 1493) {
    fFieldRun = true;
  } else if (fRunId < 1506) {  // beginning of Xe+La at 75GeV
    fFieldRun = false;
  } else {
    fFieldRun = true;
  }
  fZprim = 47;

  if (fJura)
    fParams = (Na61VdParametersManager::Instance())->GetJuraArmParams();
  else
    fParams = (Na61VdParametersManager::Instance())->GetSaleveArmParams();

  cout << "Na61VdTrackingModule::Init(): fNsig_dev=" << fNsig_dev << "   fNsig_pvert=" << fNsig_pvert << " fFieldRun=" << fFieldRun << endl;

  for (int i = 0; i < 20; i++) fAll3StationTracks[i] = 0;
}

//____________________________________________________________________
void Na61VdTrackingModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);
}

//____________________________________________________________________
void Na61VdTrackingModule::Event(UEventNode* /*inNode*/, UEventNode* /*outNode*/)

{
  // Per event method
  SetState(kEvent);
}

//_____________________________________________________________
void Na61VdTrackingModule::RefitTracks(int Ntab, UEventNode* out) {
  // algorithm:
  // Ntab=4 to use only 4 hits tracks
  // Ntab=18 to use all tracks

  for (int i = 0; i < Ntab; i++) {
    UDataTable* tracktab = out->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));

    if (!tracktab) continue;

    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);

      RefitTrack_front(track);
      RefitTrack_back(track);

      if (Ntab == 4) continue;

      double front_slope = track->GetDX_f() / track->GetDZ_f();
      double back_slope = track->GetDX_b() / track->GetDZ_b();
      fhSlopeChange->Fill(back_slope - front_slope);
      track->SetSlopeChange(back_slope - front_slope);
      track->SetLinefParams();
      track->SetLinebParams();
    }
  }
}

//_____________________________________________________________
void Na61VdTrackingModule::RefitTracks(UEventNode* out) {
  // algorithm:

  UDataTable* tracktab = out->GetDataTable(fOutTableName.Data());

  if (!tracktab) return;
  // cout<<"start with table, entries: "<<tracktab->GetEntries()<<endl;

  for (int i = 0; i < tracktab->GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracktab->At(i);

    RefitTrack_front(track);
    RefitTrack_back(track);

    double front_slope = track->GetDX_f() / track->GetDZ_f();
    double back_slope = track->GetDX_b() / track->GetDZ_b();
    fhSlopeChange->Fill(back_slope - front_slope);
    track->SetSlopeChange(back_slope - front_slope);
    track->SetLinefParams();
    track->SetLinebParams();
  }
}

//____________________________________________________________________________________
void Na61VdTrackingModule::RefitTrack_front(UVdTrack* track, int iflag) {
  int ind1 = track->GetHitIndexOnStation(0);
  int ind2 = track->GetHitIndexOnStation(1);
  int ind3 = track->GetHitIndexOnStation(2);
  int ind4 = track->GetHitIndexOnStation(3);

  int tind1 = track->GetTabIndexOnStation(0);
  int tind2 = track->GetTabIndexOnStation(1);
  int tind3 = track->GetTabIndexOnStation(2);
  int tind4 = track->GetTabIndexOnStation(3);

  USensorHit* hit1 = 0;
  USensorHit* hit2 = 0;
  USensorHit* hit3 = 0;
  USensorHit* hit4 = 0;

  if ((ind1 != -1) && (tind1 != -1) && fHitsTabs[tind1]) hit1 = (USensorHit*)fHitsTabs[tind1]->At(ind1);
  if ((ind2 != -1) && (tind1 != -1) && fHitsTabs[tind2]) hit2 = (USensorHit*)fHitsTabs[tind2]->At(ind2);
  if ((ind3 != -1) && (tind1 != -1) && fHitsTabs[tind3]) hit3 = (USensorHit*)fHitsTabs[tind3]->At(ind3);
  if ((ind4 != -1) && (tind1 != -1) && fHitsTabs[tind4]) hit4 = (USensorHit*)fHitsTabs[tind4]->At(ind4);

  // need it for tracking with HT (should we added standard hits association in VdTrackingHTModule?)
  if (!hit1) hit1 = (USensorHit*)track->GetHitIdAtStation(0);
  if (!hit2) hit2 = (USensorHit*)track->GetHitIdAtStation(1);
  if (!hit3) hit3 = (USensorHit*)track->GetHitIdAtStation(2);
  if (!hit4) hit4 = (USensorHit*)track->GetHitIdAtStation(3);

  /* // for debugging purpose
    if(hit1)cout<<"hit1: "<<hit1<<endl;
    else cout<<"hit1 is NULL"<<endl;
    if(hit2)cout<<"hit2: "<<hit2<<endl;
    else cout<<"hit2 is NULL"<<endl;
    if(hit3)cout<<"hit3: "<<hit3<<endl;
    else cout<<"hit3 is NULL"<<endl;
    if(hit4)cout<<"hit4: "<<hit4<<endl;
    else cout<<"hit4 is NULL"<<endl;
  */

  double sig[4] = {0.004, 0.0145, 0.032, 0.0416};  // errors in mm
  // double sig[4] = {0.004, 0.005, 0.006, 0.007}; // errors in mm

  USensorHit* hits[4] = {hit1, hit2, hit3, hit4};

  for (int i = 0; i < 4; i++) {
    if (hits[i]) hits[i]->SetUsed(true);  // do it only in Refit_front
  }

  double ay;
  double by;
  double N;
  double zmin;

  FitLineY_w2(hits, sig, ay, by, N, zmin);

  TF1 pol2y("pol2y", "[0] + [1]*x + [2]*x*x", -10.0, 160.0);

  double x[4];
  double y[4];
  double ex[4];
  double ey[4];

  int ii = 0;
  for (int i = 0; i < 4; i++) {
    if (hits[i]) {
      x[ii] = hits[i]->GetZ();
      ex[ii] = 0.001;  // 1 micron error on z position

      y[ii] = hits[i]->GetX();
      ey[ii] = sig[i];

      ii++;
    }
  }

  TGraphErrors* gr = new TGraphErrors(ii, x, y, ex, ey);

  // if(x[0]>10)cout<<"Tu10: "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<endl;

  pol2y.SetParameter(0, track->GetX());
  pol2y.SetParameter(1, track->GetDX() / track->GetDZ());
  if (fFieldRun)
    pol2y.SetParameter(2, 0);
  else
    pol2y.FixParameter(2, 0);

  gr->Fit("pol2y", "Q", "", x[0] - 1, x[ii - 1] + 1);

  // cout<<"Tu12: "<<fFieldRun<<endl;

  double c = pol2y.GetParameter(0);
  double b = pol2y.GetParameter(1);
  double a = pol2y.GetParameter(2);

  if (iflag == 0) fhCurvature_front->Fill(a);  // don't fill if called from CombineWithPrimaryVertex

  // if(x[0]>10)cout<<"Tu10: "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<"  a="<<a<<endl;
  // cout<<"Tu10: "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<"  a="<<a<<endl;

  double oz = track->GetZ();
  double oy = ay * oz + by;
  double ox = a * oz * oz + b * oz + c;

  // cout<<"oz="<<oz<<endl;

  // cout<<"Tu13"<<" "<<track->GetDZ()<<" "<<a<<endl;

  double tanx = 2. * a * oz + b;
  double tany = ay;

  // cout<<"Tu14: line"<<(Line3D*)track->Getline()<<endl;

  (track->Getlinef())->SetOrigin(ox, oy, oz);
  (track->Getlinef())->SetDirection(tanx, tany, 1.);
  // cout<<"Tu15: linef"<<(Line3D*)track->Getlinef()<<endl;

  ox = track->GetXatZ_f(0.0);
  oy = track->GetYatZ_f(0.0);
  (track->Getlinef())->SetOrigin(ox, oy, 0.0);

  // check fit quality

  delete gr;

  // we can try to calculate momentum here.

  if (iflag != 1) return;  // don't fill if called from CombineWithPrimaryVertex

  for (int i = 0; i < 4; i++) {
    if (hits[i]) {
      double z = hits[i]->GetZ();
      double x = hits[i]->GetX();

      if (ii == 3)
        fdx3hit[i]->Fill(x - pol2y.Eval(z));
      else
        fdx4hit_front[i]->Fill(x - pol2y.Eval(z));
    }
  }

  // delete pol2y;
}

//____________________________________________________________________________________
void Na61VdTrackingModule::RefitTrack_back(UVdTrack* track) {
  int ind1 = track->GetHitIndexOnStation(0);
  int ind2 = track->GetHitIndexOnStation(1);
  int ind3 = track->GetHitIndexOnStation(2);
  int ind4 = track->GetHitIndexOnStation(3);

  int tind1 = track->GetTabIndexOnStation(0);
  int tind2 = track->GetTabIndexOnStation(1);
  int tind3 = track->GetTabIndexOnStation(2);
  int tind4 = track->GetTabIndexOnStation(3);

  USensorHit* hit1 = 0;
  USensorHit* hit2 = 0;
  USensorHit* hit3 = 0;
  USensorHit* hit4 = 0;

  if ((ind1 != -1) && (tind1 != -1) && fHitsTabs[tind1]) hit1 = (USensorHit*)fHitsTabs[tind1]->At(ind1);
  if ((ind2 != -1) && (tind1 != -1) && fHitsTabs[tind2]) hit2 = (USensorHit*)fHitsTabs[tind2]->At(ind2);
  if ((ind3 != -1) && (tind1 != -1) && fHitsTabs[tind3]) hit3 = (USensorHit*)fHitsTabs[tind3]->At(ind3);
  if ((ind4 != -1) && (tind1 != -1) && fHitsTabs[tind4]) hit4 = (USensorHit*)fHitsTabs[tind4]->At(ind4);

  // need it for tracking with HT (should we add standard hits association in VdTrackingHTModule?)
  if (!hit1) hit1 = (USensorHit*)track->GetHitIdAtStation(0);
  if (!hit2) hit2 = (USensorHit*)track->GetHitIdAtStation(1);
  if (!hit3) hit3 = (USensorHit*)track->GetHitIdAtStation(2);
  if (!hit4) hit4 = (USensorHit*)track->GetHitIdAtStation(3);

  // double sig[4] = {0.004, 0.0145, 0.032, 0.0416}; // errors in mm
  double sig[4] = {0.007, 0.006, 0.005, 0.004};  // errors in mm

  USensorHit* hits[4] = {hit1, hit2, hit3, hit4};

  double ay;
  double by;
  double N;
  double zmin;

  FitLineY_w2(hits, sig, ay, by, N, zmin);

  // TF1* pol2y = new TF1("pol2y","[0] + [1]*x + [2]*x*x",-10.0,160.0);
  TF1 pol2y("pol2y", "[0] + [1]*x + [2]*x*x", -10.0, 160.0);

  double x[4];
  double y[4];
  double ex[4];
  double ey[4];

  int ii = 0;
  for (int i = 0; i < 4; i++) {
    if (hits[i]) {
      x[ii] = hits[i]->GetZ();
      ex[ii] = 0.001;  // 1 micron error on z position

      y[ii] = hits[i]->GetX();
      ey[ii] = sig[i];

      ii++;
    }
  }

  TGraphErrors* gr = new TGraphErrors(ii, x, y, ex, ey);

  // cout<<"Tu3: "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<endl;

  pol2y.SetParameter(0, track->GetX());
  pol2y.SetParameter(1, track->GetDX() / track->GetDZ());
  if (fFieldRun)
    pol2y.SetParameter(2, 0);
  else
    pol2y.FixParameter(2, 0);

  gr->Fit("pol2y", "Q", "", x[0] - 1, x[ii - 1] + 1);

  double c = pol2y.GetParameter(0);
  double b = pol2y.GetParameter(1);
  double a = pol2y.GetParameter(2);

  // x[ii-1] contains z coordinate of the last hit
  double oz = x[ii - 1];
  double oy = by + ay * oz;
  double ox = a * oz * oz + b * oz + c;

  // cout<<"Tu4"<<" "<<track->GetDZ()<<" "<<a<<" zlast="<<oz<<endl;

  double tanx = 2. * a * oz + b;
  double tany = ay;

  // set track charge
  fhCurvature_back->Fill(a);
  track->SetCurvature(a);
  track->SetPol2(pol2y);

  if (a > 0)
    track->SetCharge(1);
  else
    track->SetCharge(-1);

  // cout<<"Tu5: line"<<(Line3D*)track->Getline()<<endl;

  (track->Getlineb())->SetOrigin(ox, oy, oz);
  (track->Getlineb())->SetDirection(tanx, tany, 1.);
  // cout<<"Tu6: lineb"<<(Line3D*)track->Getlineb()<<endl;

  for (int i = 0; i < 4; i++) {
    if (hits[i]) {
      double z = hits[i]->GetZ();
      double x = hits[i]->GetX();

      if (ii == 4) fdx4hit_back[i]->Fill(x - pol2y.Eval(z));
    }
  }

  // delete pol2y;
  delete gr;
}

//___________________________________________________________________________________
void Na61VdTrackingModule::FitLineY_w2(USensorHit** hits, double* sig, double& ay, double& by, double& N, double& zmin) {
  // simple regration (assuming constant errors)
  // It is quite straight forward to use weighted regression

  double hy, hz;

  double S = 0;
  double Sz = 0;
  double Szz = 0;

  double Sy = 0;
  double Szy = 0;

  zmin = 10000.;
  N = 0;

  // perhaps we need waited regression

  for (int i = 0; i < 4; i++) {
    USensorHit* hit = hits[i];

    if (!hit) continue;  // now we always have 4 hits

    N++;

    // smearing is done in dedicated method
    hy = hit->GetY();
    hz = hit->GetZ();

    double sig2 = sig[i] * sig[i];

    S += 1. / sig2;
    Sz += hz / sig2;
    Szz += (hz * hz) / sig2;

    Sy += hy / sig2;
    Szy += (hz * hy) / sig2;

    if (hz < zmin) zmin = hz;
  }

  ay = (S * Szy - Sz * Sy) / (S * Szz - Sz * Sz);
  by = (Szz * Sy - Sz * Szy) / (S * Szz - Sz * Sz);
}

//_____________________________________________________________
int Na61VdTrackingModule::GetAllTracks(UEventNode* out) {
  // algorithm:

  int fullTracks = 0;

  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = out->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));

    if (!tracktab) continue;

    for (int j = 0; j < tracktab->GetEntries(); j++) {
      // UVdTrack* track = (UVdTrack*)tracktab->At(j);

      fullTracks++;
    }
  }

  return fullTracks;
}

//_____________________________________________________________
void Na61VdTrackingModule::Delete3HitTracks(UEventNode* out) {
  int NDev = 20;

  for (int i = 0; i < NDev; i++) {
    UDataTable* table = out->GetDataTable(Form("Tracks %s", fDevNames[i].Data()));

    if (table) {
      // cout<<" table name: "<<table->GetName()<<endl;
      // table->Clear();
      table->Delete();
      if (table) out->RemoveDataObject((UDataObject*)table);
    }
  }
}

//_____________________________________________________________
void Na61VdTrackingModule::CombineWithPrimaryVertex(TString matchstr, UEventNode* out) {
  UDataTable* tabfull = new UDataTable(Form("FullTracks %s", matchstr.Data()));
  tabfull->SetOwner();
  out->AddDataTable(tabfull);

  int imatch = GetMatchIdev(matchstr);
  int idev = imatch + 2;
  UDataTable* tab3hit = out->GetDataTable(Form("Tracks %s", fDevNames[idev].Data()));

  double Offx = fParams->GetCombineOffsetX(imatch);
  double Sigx = fParams->GetCombineSigmaX(imatch);
  double Offy = fParams->GetCombineOffsetY(imatch);
  double Sigy = fParams->GetCombineSigmaY(imatch);

  for (int i = 0; i < tab3hit->GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tab3hit->At(i);
    RefitTrack_front(track, 1);  // it should make front line

    double ax = track->GetDX_f() / track->GetDZ_f();
    double ay = track->GetDY_f() / track->GetDZ_f();

    if (fhAx[imatch]) fhAx[imatch]->Fill(ax);
    if (fhAy[imatch]) fhAy[imatch]->Fill(ay);

    // changed for Vds4/Vds3 to Vds1/Vds2 (should make distributions narrower)
    USensorHit* hitx = (USensorHit*)track->GetHitIdAtStation(0);
    if (!hitx) hitx = (USensorHit*)track->GetHitIdAtStation(1);

    USensorHit* hity = (USensorHit*)track->GetHitIdAtStation(3);
    if (!hity) hity = (USensorHit*)track->GetHitIdAtStation(2);

    Vector3D* pvert = GetPrimaryVertex();
    double x4 = hitx->GetX();
    double y4 = hity->GetY();
    double z4x = hitx->GetZ();
    double z4y = hity->GetZ();
    double xp = pvert->GetX();
    double yp = pvert->GetY();
    double zp = pvert->GetZ();

    double dx = track->GetXatZ_f(zp) - xp;
    double dy = track->GetYatZ_f(zp) - yp;

    double axp = (x4 - xp) / (z4x - zp);
    double ayp = (y4 - yp) / (z4y - zp);
    // check additional cut - based on primary vertex - for 3Hit tracks

    double dax = axp - ax;
    double day = ayp - ay;

    if (fhDAx[imatch]) fhDAx[imatch]->Fill(dax);
    if (fhDAy[imatch]) fhDAy[imatch]->Fill(day);

    if (fhDx[imatch]) fhDx[imatch]->Fill(dx);
    if (fhDy[imatch]) fhDy[imatch]->Fill(dy);

    if (TMath::Abs(dax - Offx) < 3 * Sigx) fhDAy_cutx[imatch]->Fill(day);
    if (TMath::Abs(day - Offy) < 3 * Sigy) fhDAx_cuty[imatch]->Fill(dax);

    double da2 = (dax - Offx) * (dax - Offx) / (Sigx * Sigx) + (day - Offy) * (day - Offy) / (Sigy * Sigy);

    if (da2 < fNsig_pvert * fNsig_pvert) {  // cases in pick

      if (fhDx_acce[imatch]) fhDx_acce[imatch]->Fill(dx);
      if (fhDy_acce[imatch]) fhDy_acce[imatch]->Fill(dy);
      // for accepted tracks fill deviations

      // make VdTrack Full Track
      // cout<<"making full track"<<endl;
      CreateFullTrackx_new(imatch, track, tabfull);
    }
  }
}

//_____________________________________________________________
void Na61VdTrackingModule::CombineTracks(TString matchname, int idev1, int idev2, UEventNode* out) {
  // For now just combine dev1_0 with dev2_0 and  dev1_1 with dev2_1
  // Added combining dev1_1 with dev2_3

  UDataTable* tracks_dev1 = out->GetDataTable(Form("Tracks %s", fDevNames[idev1].Data()));
  UDataTable* tracks_dev2 = out->GetDataTable(Form("Tracks %s", fDevNames[idev2].Data()));

  if (!tracks_dev1) return;
  if (!tracks_dev2) return;

  int imatch = GetMatchIdev(matchname);

  UDataTable* tracktab = new UDataTable(Form("FullTracks %s", matchname.Data()));
  tracktab->SetOwner();

  // double Off_dx  = fParams->GetDxOffset(imatch);    double Sig_dx  = fParams->GetDxSigma(imatch);
  // double Off_dy  = fParams->GetDyOffset(imatch);    double Sig_dy  = fParams->GetDySigma(imatch);
  // double Off_dax = fParams->GetDaxOffset(imatch);   double Sig_dax = fParams->GetDaxSigma(imatch);
  // double Off_day = fParams->GetDayOffset(imatch);   double Sig_day = fParams->GetDaySigma(imatch);

  int Ndev1 = tracks_dev1->GetEntries();
  int Ndev2 = tracks_dev2->GetEntries();

  for (int i = 0; i < Ndev1; i++) {
    UVdTrack* track1 = (UVdTrack*)tracks_dev1->At(i);

    double tx1 = track1->GetXatZ(125.);
    double ty1 = track1->GetYatZ(125.);
    double ax1 = track1->GetDX() / track1->GetDZ();
    double ay1 = track1->GetDY() / track1->GetDZ();

    long long id1 = (long long)track1->GetHitIdAtStation(1);
    long long id2 = (long long)track1->GetHitIdAtStation(2);

    for (int j = 0; j < Ndev2; j++) {
      UVdTrack* track2 = (UVdTrack*)tracks_dev2->At(j);

      long long idd1 = (long long)track2->GetHitIdAtStation(1);
      long long idd2 = (long long)track2->GetHitIdAtStation(2);

      double tx2 = track2->GetXatZ(125.);
      double ty2 = track2->GetYatZ(125.);
      double ax2 = track2->GetDX() / track2->GetDZ();
      double ay2 = track2->GetDY() / track2->GetDZ();

      double dx = tx2 - tx1;
      double dy = ty2 - ty1;
      double dax = ax2 - ax1;
      double day = ay2 - ay1;

      // cout<<fhDx_match[imatch]<<" "<<fhDy_match[imatch]<<" "<<fhDax_match[imatch]<<" "<<fhDay_match[imatch]<<endl;

      if (fhDx_match[imatch]) fhDx_match[imatch]->Fill(dx);
      if (fhDy_match[imatch]) fhDy_match[imatch]->Fill(dy);
      if (fhDax_match[imatch]) fhDax_match[imatch]->Fill(dax);
      if (fhDay_match[imatch]) fhDay_match[imatch]->Fill(day);

      if (id1 != idd1) continue;

      if (fhDx_match_cut1[imatch]) fhDx_match_cut1[imatch]->Fill(dx);
      if (fhDy_match_cut1[imatch]) fhDy_match_cut1[imatch]->Fill(dy);
      if (fhDax_match_cut1[imatch]) fhDax_match_cut1[imatch]->Fill(dax);
      if (fhDay_match_cut1[imatch]) fhDay_match_cut1[imatch]->Fill(day);

      if (id2 != idd2) continue;

      // we don't base on dev2 cut but keep it for any reason
      // double dev2 = (dx-Off_dx)*(dx-Off_dx)/(Sig_dx*Sig_dx) + (dy-Off_dy)*(dy-Off_dy)/(Sig_dy*Sig_dy) +
      //(dax-Off_dax)*(dax-Off_dax)/(Sig_dax*Sig_dax) + (day-Off_day)*(day-Off_day)/(Sig_day*Sig_day);

      // if(dev2 > 5*5) continue;

      if (fhDx_match_acce[imatch]) fhDx_match_acce[imatch]->Fill(dx);
      if (fhDy_match_acce[imatch]) fhDy_match_acce[imatch]->Fill(dy);
      if (fhDax_match_acce[imatch]) fhDax_match_acce[imatch]->Fill(dax);
      if (fhDay_match_acce[imatch]) fhDay_match_acce[imatch]->Fill(day);

      // fill deviation

      double x1 = ((USensorHit*)track1->GetHitIdAtStation(0))->GetX();
      double x2 = ((USensorHit*)track1->GetHitIdAtStation(1))->GetX();
      double x3 = ((USensorHit*)track1->GetHitIdAtStation(2))->GetX();
      double y1 = ((USensorHit*)track1->GetHitIdAtStation(0))->GetY();
      double y2 = ((USensorHit*)track1->GetHitIdAtStation(1))->GetY();
      double y3 = ((USensorHit*)track1->GetHitIdAtStation(2))->GetY();
      double z1 = ((USensorHit*)track1->GetHitIdAtStation(0))->GetZ();
      double z2 = ((USensorHit*)track1->GetHitIdAtStation(1))->GetZ();
      double z3 = ((USensorHit*)track1->GetHitIdAtStation(2))->GetZ();

      double devx = ((z2 - z1) * x3 + (z3 - z2) * x1) / (z3 - z1) - x2;
      double devy = ((z2 - z1) * y3 + (z3 - z2) * y1) / (z3 - z1) - y2;

      if (fhX_dev_acce[idev1]) fhX_dev_acce[idev1]->Fill(devx);
      if (fhY_dev_acce[idev1]) fhY_dev_acce[idev1]->Fill(devy);

      // if(fhVds1_xy[imatch]) fhVds1_xy[imatch]_Fill(x1,y1);
      // if(fhVds2_xy[imatch]) fhVds2_xy[imatch]_Fill(x2,y2);
      // if(fhVds3_xy[imatch]) fhVds3_xy[imatch]_Fill(x3,y3);

      x1 = ((USensorHit*)track2->GetHitIdAtStation(1))->GetX();
      x2 = ((USensorHit*)track2->GetHitIdAtStation(2))->GetX();
      x3 = ((USensorHit*)track2->GetHitIdAtStation(3))->GetX();
      y1 = ((USensorHit*)track2->GetHitIdAtStation(1))->GetY();
      y2 = ((USensorHit*)track2->GetHitIdAtStation(2))->GetY();
      y3 = ((USensorHit*)track2->GetHitIdAtStation(3))->GetY();
      z1 = ((USensorHit*)track2->GetHitIdAtStation(1))->GetZ();
      z2 = ((USensorHit*)track2->GetHitIdAtStation(2))->GetZ();
      z3 = ((USensorHit*)track2->GetHitIdAtStation(3))->GetZ();

      devx = ((z2 - z1) * x3 + (z3 - z2) * x1) / (z3 - z1) - x2;
      devy = ((z2 - z1) * y3 + (z3 - z2) * y1) / (z3 - z1) - y2;

      if (fhX_dev_acce[idev2]) fhX_dev_acce[idev2]->Fill(devx);
      if (fhY_dev_acce[idev2]) fhY_dev_acce[idev2]->Fill(devy);

      // if(fhVds4_xy[imatch]) fhVds4_xy[imatch]_Fill(x3,y3);

      // create 4hit tracks here.
      CreateFullTrack_new(imatch, track1, track2, tracktab);
    }
  }

  out->AddDataTable(tracktab);

  // cout<<imatch<<" "<<fAllFullTracks[imatch]<<" "<<tracktab->GetEntries()<<" "<<fhDx_match_acce[imatch]->GetEntries()<<endl;
}

//_____________________________________________________________
void Na61VdTrackingModule::FindPrimaryVertexPost(UEventNode* out) {
  // algorithm:

  int fullTracks = 0;
  TObjArray tracks;
  tracks.Clear();

  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = out->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
    if (!tracktab) continue;

    fullTracks = fullTracks + tracktab->GetEntries();
    for (int i = 0; i < tracktab->GetEntries(); i++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(i);
      if (!track->IsMarkedForRemoval()) tracks.Add(tracktab->At(i));
    }
  }

   cout<<"FindPrimaryVertexPost:     fhFullTracksAll: "<<fhFullTracksAll<<" "<<fullTracks<<endl;

  fhFullTracksAll->Fill(fullTracks);

  // cout<<"tu1 "<<endl;

  if (tracks.GetEntries() < 2) return;

  double sum_up_x, sum_do_x, sum_up_y, sum_do_y;

  CalculateSums(sum_up_x, sum_do_x, sum_up_y, sum_do_y, tracks);

  ////////// new method //////////////////////////////////////////////////

  //  double zx_min = 0;
  //  if(sum_do_x > 0)zx_min = - sum_up_x/sum_do_x;
  //  double zy_min = 0;
  //  if(sum_do_y > 0)zy_min = - sum_up_y/sum_do_y;

  // cout<<sum_do_x<<"  "<<sum_do_y<<"  "<<zx_min<<" "<<zy_min<<endl;

  double z_prim = -(sum_up_x + sum_up_y) / (sum_do_x + sum_do_y);

  Vector3D reco_pvertex_w2;

  PrimaryVertexWithWeigths(tracks, z_prim, reco_pvertex_w2);

  fhRecoVertexXY_w2_Post->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
  fhRecoVertexZ_w2_Post->Fill(reco_pvertex_w2.Z());
  fhRecoVertexZ_fine_w2_Post->Fill(reco_pvertex_w2.Z());

  fPrimaryVertexDefined = true;
  ((UVdEvent*)out)->SetPrimaryVertexStatus(2);
  ((UVdEvent*)out)->SetPrimaryVertex(reco_pvertex_w2.X(), reco_pvertex_w2.Y(), reco_pvertex_w2.Z());
  fPrimaryVertex->SetX(reco_pvertex_w2.X());
  fPrimaryVertex->SetY(reco_pvertex_w2.Y());
  fPrimaryVertex->SetZ(reco_pvertex_w2.Z());
}

//_____________________________________________________________
void Na61VdTrackingModule::FindPrimaryVertex(UEventNode* out) {
  TString matchname[5] = {"down1", "up1", "down2", "up2"};
  UDataTable* tracktab[4];
  tracktab[0] = out->GetDataTable(Form("FullTracks %s", matchname[0].Data()));
  tracktab[1] = out->GetDataTable(Form("FullTracks %s", matchname[1].Data()));
  tracktab[2] = out->GetDataTable(Form("FullTracks %s", matchname[2].Data()));
  tracktab[3] = out->GetDataTable(Form("FullTracks %s", matchname[3].Data()));

  int fullTracks = 0;
  int Entries[4] = {0, 0, 0, 0};
  for (int i = 0; i < 4; i++) {
    if (tracktab[i]) {
      Entries[i] = tracktab[i]->GetEntries();
      fullTracks = fullTracks + tracktab[i]->GetEntries();
    }
  }

  fhFullTracks->Fill(fullTracks);
  
     cout<<"FindPrimaryVertex:     fhFullTracks:"<<fullTracks<<endl;


  TObjArray tracks;
  tracks.Clear();
  for (int i = 0; i < Entries[0]; i++) {
    UVdTrack* track = (UVdTrack*)tracktab[0]->At(i);
    if (!track->IsMarkedForRemoval()) tracks.Add(tracktab[0]->At(i));
  }
  for (int i = 0; i < Entries[1]; i++) {
    UVdTrack* track = (UVdTrack*)tracktab[1]->At(i);
    if (!track->IsMarkedForRemoval()) tracks.Add(tracktab[1]->At(i));
  }
  for (int i = 0; i < Entries[2]; i++) {
    UVdTrack* track = (UVdTrack*)tracktab[2]->At(i);
    if (!track->IsMarkedForRemoval()) tracks.Add(tracktab[2]->At(i));
  }
  for (int i = 0; i < Entries[3]; i++) {
    UVdTrack* track = (UVdTrack*)tracktab[3]->At(i);
    if (!track->IsMarkedForRemoval()) tracks.Add(tracktab[3]->At(i));
  }

  if (tracks.GetEntries() < 2) return;

  double sum_up_x, sum_do_x, sum_up_y, sum_do_y;

  CalculateSums(sum_up_x, sum_do_x, sum_up_y, sum_do_y, tracks);

  ////////// new method //////////////////////////////////////////////////

  double zx_min = 0;
  if (sum_do_x > 0) zx_min = -sum_up_x / sum_do_x;
  double zy_min = 0;
  if (sum_do_y > 0) zy_min = -sum_up_y / sum_do_y;

  // cout<<sum_do_x<<"  "<<sum_do_y<<"  "<<zx_min<<" "<<zy_min<<endl;
  fhDZmin->Fill(zy_min - zx_min);

  fhRecoVertexZ_fine_x->Fill(zx_min);
  fhRecoVertexZ_fine_y->Fill(zy_min);

  double z_prim = -(sum_up_x + sum_up_y) / (sum_do_x + sum_do_y);

  fhRecoVertexZ->Fill(z_prim);
  fhRecoVertexZ_fine->Fill(z_prim);

  Vector3D reco_pvertex_w2;
  PrimaryVertexWithWeigths(tracks, z_prim, reco_pvertex_w2);

  // cout<<reco_pvertex_w2.X()<<"  "<<reco_pvertex_w2.Y()<<endl;
  fhRecoVertexXY_w2->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
  fhRecoVertexZ_w2->Fill(reco_pvertex_w2.Z());
  fhRecoVertexZ_fine_w2->Fill(reco_pvertex_w2.Z());

  fPrimaryVertexDefined = true;
  fPrimaryVertex->SetX(reco_pvertex_w2.X());
  fPrimaryVertex->SetY(reco_pvertex_w2.Y());
  fPrimaryVertex->SetZ(reco_pvertex_w2.Z());
  ((UVdEvent*)out)->SetPrimaryVertexStatus(1);
  ((UVdEvent*)out)->SetPrimaryVertex(reco_pvertex_w2.X(), reco_pvertex_w2.Y(), reco_pvertex_w2.Z());
}

//______________________________________________________________________________________________________
void Na61VdTrackingModule::CalculateSums(double& sum_up_x, double& sum_do_x, double& sum_up_y, double& sum_do_y, TObjArray tracks) {
  sum_up_x = 0;
  sum_do_x = 0;
  sum_up_y = 0;
  sum_do_y = 0;

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      /////////////////// stuff for new method
      double axi = track1->GetDX_f() / track1->GetDZ_f();
      double axj = track2->GetDX_f() / track2->GetDZ_f();
      double bxi = track1->GetX_f();
      double bxj = track2->GetX_f();

      double ayi = track1->GetDY_f() / track1->GetDZ();
      double ayj = track2->GetDY_f() / track2->GetDZ();
      double byi = track1->GetY_f();
      double byj = track2->GetY_f();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < 3) {
        sum_up_x = sum_up_x + (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - fZprim) < 3) {
        sum_up_y = sum_up_y + (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + (ayi - ayj) * (ayi - ayj);
      }

      // cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
    }
  }
}

//______________________________________________________________________________________________________
void Na61VdTrackingModule::PrimaryVertexWithWeigths(TObjArray& tracks, double zprim, Vector3D& pvertex) {
  // algorithm:
  // 1. Calculate weights = 1/dist^2, distance betweem i-th and j-th tracks on
  // the closest proximity plane (z=zprim).
  // 2. Used the weights to calculate z coordinate of primary vertex from the formula that inclused weights (see PS notes).
  // 3. Find x,y coordinates averageing tracks x, y at the closes proximmity plane.

  double sum_up_x = 0;
  double sum_do_x = 0;
  double sum_up_y = 0;
  double sum_do_y = 0;

  // cout<<" entries="<<tracks.GetEntries()<<" "<<zprim<<endl;

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      // cout<<"i="<<i<<" j="<<j<<endl;

      double x1 = track1->GetXatZ(zprim);
      double x2 = track2->GetXatZ(zprim);
      double y1 = track1->GetYatZ(zprim);
      double y2 = track2->GetYatZ(zprim);

      // cout<<"x1="<<x1<<" x2="<<x2<<"  y1="<<y1<<" y2="<<y2<<endl;

      double d2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
      double w = 1. / d2;

      // cout<<"dd2="<<d2<<endl;

      fhd2->Fill(d2);

      // cout<<"d2="<<d2<<endl;

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

      if (TMath::Abs(z_x - 50) < 50) {
        sum_up_x = sum_up_x + w * (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + w * (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - 50) < 50) {
        sum_up_y = sum_up_y + w * (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + w * (ayi - ayj) * (ayi - ayj);
      }
      // cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
    }
  }

  double zprim_w2 = -(sum_up_x + sum_up_y) / (sum_do_x + sum_do_y);
  // cout<<" zprim_w2 ="<<zprim_w2<<endl;

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

//______________________________________________________________________________________________________
void Na61VdTrackingModule::RemoveGhostTracks(TString devname1, TString devname2, UEventNode* out) {
  // This method resolve problem of same track found in devname1 and devname2 tables of 3Hit tracks.
  // This might happen due to overlap sensors in Vds4 station

  UDataTable* tracktab1 = out->GetDataTable(Form("Tracks %s", devname1.Data()));
  UDataTable* tracktab2 = out->GetDataTable(Form("Tracks %s", devname2.Data()));

  // cout<<"---------------> RemoveGhostTracks: "<<devname1.Data()<<"   "<<devname2.Data()<<endl;

  int N1 = tracktab1->GetEntries();
  int N2 = tracktab2->GetEntries();

  // cout<<"tracktab1: "<<tracktab1->GetName()<<"     tracktab2: "<<tracktab2->GetName()<<endl;
  // cout<<"N1="<<N1<<"  N2="<<N2<<endl;

  for (int i = 0; i < N1; i++) {
    UVdTrack* track1 = (UVdTrack*)tracktab1->At(i);
    for (int j = 0; j < N2; j++) {
      UVdTrack* track2 = (UVdTrack*)tracktab2->At(j);

      // cout<<"i="<<i<<"  j="<<j<<endl;

      if (track1->IsMarkedForRemoval() || track2->IsMarkedForRemoval()) continue;

      CheckTrackMatching2(track1, track2);

      // if(track1->IsMarkedForRemoval() || track2->IsMarkedForRemoval()){
      // cout<<" track1: "<<track1->GetX()<<" "<<track1->GetY()<<" "<<track1->GetZ()<<endl;
      // cout<<" track2: "<<track2->GetX()<<" "<<track2->GetY()<<" "<<track2->GetZ()<<endl;
      //}
    }
  }

  for (int i = 0; i < N1; i++) {
    UVdTrack* track = (UVdTrack*)tracktab1->At(i);
    if (track->IsMarkedForRemoval()) {
      tracktab1->DeleteObjectAt(i);
    }
  }
  tracktab1->Compress();

  for (int i = 0; i < N2; i++) {
    UVdTrack* track = (UVdTrack*)tracktab2->At(i);
    if (track->IsMarkedForRemoval()) {
      tracktab2->DeleteObjectAt(i);
    }
  }
  tracktab2->Compress();

  // int rt = N1 - tracktab1->GetEntries();
  // rt =  rt + N2 - tracktab2->GetEntries();

  // cout<<"removed tracks: "<<rt<<endl;
}

//______________________________________________________________________________________________________
void Na61VdTrackingModule::Make3HitTracks(TString devname, int csize, UDataTable* hits1, UDataTable* hits2, UDataTable* hits3, UEventNode* out, int* TabIndArray) {
  fTrackID = 0;

  int idev = GetIdev(devname);

  fLargeClusters = false;
  TH1F* hX_dev = fhX_dev[idev];
  TH1F* hY_dev = fhY_dev[idev];
  TH1F* hX_dev_cuty = fhX_dev_cuty[idev];
  TH1F* hY_dev_cutx = fhY_dev_cutx[idev];
  //  TH2F* h_devx_devy = fh_devx_devy[idev];

  // cout<<"9 fhFullTracksAll: "<<fhFullTracksAll<<endl;
  // cout<<"idev="<<idev<<"  "<<devname.Data()<<"  h_devx_devy: "<<h_devx_devy<<"  hX_dev_cuty:"<<hX_dev_cuty<<endl;

  if (csize > 10000) {
    fLargeClusters = true;
    hX_dev = fhX_dev_lc[idev];
    hY_dev = fhY_dev_lc[idev];
    hX_dev_cuty = fhX_dev_cuty_lc[idev];
    hY_dev_cutx = fhY_dev_cutx_lc[idev];
  }

  // cout<<"91 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  double Offsetx = fParams->GetDevOffsetX(idev);
  double Sigmax = fParams->GetDevSigmaX(idev);
  double Offsety = fParams->GetDevOffsetY(idev);
  double Sigmay = fParams->GetDevSigmaY(idev);

  UDataTable* tracktab = 0;
  if (fLargeClusters) {
    tracktab = new UDataTable(Form("Tracks_lc %s", devname.Data()));
  } else {
    tracktab = new UDataTable(Form("Tracks %s", devname.Data()));
  }
  if (tracktab) tracktab->SetOwner(true);

  //
  //  bool beamTrack=false;
  int IndArray[3];

  // cout<<"92 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  // cout<<"Tu1"<<" "<<Offsetx<<" "<<Offsety<<" "<<Sigmax<<" "<<Sigmay<<endl;
  for (int i = 0; i < hits1->GetEntries(); i++) {
    USensorHit* hit1 = (USensorHit*)hits1->At(i);
    if ((hit1->GetUniqueID()) && (fProductionMode)) continue;
    if (hit1->GetClusterSize() < csize) continue;

    // cout<<"i="<<i<<" UniqueID="<<hit1->GetUniqueID()<<" "<<fProductionMode<<endl;

    double x1 = hit1->GetX();
    double y1 = hit1->GetY();
    double z1 = hit1->GetZ();

    for (int j = 0; j < hits2->GetEntries(); j++) {
      USensorHit* hit2 = (USensorHit*)hits2->At(j);
      if ((hit2->GetUniqueID()) && (fProductionMode)) continue;
      if (hit2->GetClusterSize() < csize) continue;
      // cout<<"j="<<j<<endl;

      double x2 = hit2->GetX();

      // if(!fLargeClusters && (x2>x1))continue; // need more sophisticated condition

      double y2 = hit2->GetY();
      double z2 = hit2->GetZ();

      int imatch = 0;

      for (int k = 0; k < hits3->GetEntries(); k++) {
        USensorHit* hit3 = (USensorHit*)hits3->At(k);
        if ((hit3->GetUniqueID()) && (fProductionMode)) continue;
        if (hit3->GetClusterSize() < csize) continue;

        double x3 = hit3->GetX();

        double y3 = hit3->GetY();
        double z3 = hit3->GetZ();

        double devx = ((z2 - z1) * x3 + (z3 - z2) * x1) / (z3 - z1) - x2;
        double devy = ((z2 - z1) * y3 + (z3 - z2) * y1) / (z3 - z1) - y2;

        if (hX_dev) hX_dev->Fill(devx);
        if (hY_dev) hY_dev->Fill(devy);

        if (TMath::Abs(devx - Offsetx) < 3 * Sigmax) hY_dev_cutx->Fill(devy);
        if (TMath::Abs(devy - Offsety) < 3 * Sigmay) hX_dev_cuty->Fill(devx);

        double dev2 = (devx - Offsetx) * (devx - Offsetx) / (Sigmax * Sigmax) + (devy - Offsety) * (devy - Offsety) / (Sigmay * Sigmay);

        if (dev2 > fNsig_dev * fNsig_dev) continue;
        // found track !!!!

        imatch++;

        IndArray[0] = i;
        IndArray[1] = j;
        IndArray[2] = k;

        // cout<<" tabIdArrays: "<<TabIndArray[0]<<"  "<<TabIndArray[1]<<" "<<TabIndArray[2]<<endl;
        Create3HitTrack(idev, hit1, hit2, hit3, IndArray, TabIndArray, tracktab);
      }
    }
  }

  // cout<<"93 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  // clean tracks
  int Ns = tracktab->GetEntries();
  for (int i = 0; i < Ns; i++) {
    UVdTrack* track1 = (UVdTrack*)tracktab->At(i);
    for (int j = i + 1; j < Ns; j++) {
      UVdTrack* track2 = (UVdTrack*)tracktab->At(j);

      if (track1->IsMarkedForRemoval() || track2->IsMarkedForRemoval()) continue;

      CheckTrackMatching(track1, track2);
    }
  }

  // cout<<"94 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  for (int i = 0; i < Ns; i++) {
    UVdTrack* track = (UVdTrack*)tracktab->At(i);
    if (track->IsMarkedForRemoval()) {
      tracktab->DeleteObjectAt(i);
    }
  }
  tracktab->Compress();

  int Ne = tracktab->GetEntries();
  if (fhNsVsNe[idev]) fhNsVsNe[idev]->Fill(Ns, Ne);
  if (fhTracks[idev]) fhTracks[idev]->Fill(Ne);
  fAll3StationTracks[idev] += Ne;

  out->AddDataTable(tracktab);
}

//_____________________________________________________________
void Na61VdTrackingModule::CreateFullTrack(int imatch, UVdTrack* track1, UVdTrack* track2, UDataTable* tracktab) {
  USensorHit* hit1 = (USensorHit*)track1->GetHitIdAtStation(0);
  USensorHit* hit2 = (USensorHit*)track2->GetHitIdAtStation(1);
  USensorHit* hit3 = (USensorHit*)track2->GetHitIdAtStation(2);
  USensorHit* hit4 = (USensorHit*)track2->GetHitIdAtStation(3);

  ((USensorHit*)track1->GetHitIdAtStation(1))->SetUniqueID(1);
  ((USensorHit*)track1->GetHitIdAtStation(2))->SetUniqueID(1);

  int ind1 = track1->GetHitIndexOnStation(0);
  int ind2 = track2->GetHitIndexOnStation(1);
  int ind3 = track2->GetHitIndexOnStation(2);
  int ind4 = track2->GetHitIndexOnStation(3);

  int tind1 = track1->GetTabIndexOnStation(0);
  int tind2 = track2->GetTabIndexOnStation(1);
  int tind3 = track2->GetTabIndexOnStation(2);
  int tind4 = track2->GetTabIndexOnStation(3);

  if (track1->GetHitIdAtStation(1) != track2->GetHitIdAtStation(1)) cout << " Na61VdTrackingModule::CreateFullTrack: hits ID differnt on Vds2 station. In principal should not happened" << endl;

  TF1* linex = new TF1("linex", "[0] +[1]*x", -10.0, 160.0);
  TF1* liney = new TF1("liney", "[0] +[1]*x", -10.0, 160.0);

  double x[] = {hit1->GetZ(), hit2->GetZ(), hit3->GetZ(), hit4->GetZ()};
  double yx[] = {hit1->GetX(), hit2->GetX(), hit3->GetX(), hit4->GetX()};
  double ex[] = {0, 0, 0, 0};
  double eyx[] = {0.004, 0.004, 0.004, 0.004};

  double yy[] = {hit1->GetY(), hit2->GetY(), hit3->GetY(), hit4->GetY()};
  double ey[] = {0, 0, 0, 0};
  double eyy[] = {0.004, 0.004, 0.004, 004};

  TGraphErrors* gr_x = new TGraphErrors(4, x, yx, ex, eyx);

  TGraphErrors* gr_y = new TGraphErrors(4, x, yy, ey, eyy);

  linex->SetParameter(0, yx[0]);

  gr_x->Fit("linex", "Q", "", -10, 160);

  double x_Vds1 = linex->Eval(0.);

  double tan_x = linex->GetParameter(1);

  liney->SetParameter(0, yy[0]);
  gr_y->Fit("liney", "Q", "", -10, 160);

  double y_Vds1 = liney->Eval(0.);

  double tan_y = liney->GetParameter(1);

  if (fhAx_full[imatch]) fhAx_full[imatch]->Fill(tan_x);
  if (fhAy_full[imatch]) fhAy_full[imatch]->Fill(tan_y);

  Vector3D origin(x_Vds1, y_Vds1, 0.);
  Vector3D direction(tan_x, tan_y, 1.);
  UVdTrack* track = new UVdTrack(origin, direction);

  tracktab->Add(track);

  track->SetTrackID((long long)track);
  track->SetPdgID(0);

  track->SetVdsHitIDs(hit1, hit2, hit3, hit4);
  track->SetHitArrayIndexes(ind1, ind2, ind3, ind4);
  track->SetTabArrayIndexes(tind1, tind2, tind3, tind4);
  hit1->SetUniqueID(1);
  hit2->SetUniqueID(1);
  hit3->SetUniqueID(1);
  hit4->SetUniqueID(1);

  // make deviation diagnostics
  fhX_dev_full[imatch]->Fill(hit1->GetX() - linex->Eval(hit1->GetZ()));
  fhX_dev_full[imatch]->Fill(hit2->GetX() - linex->Eval(hit2->GetZ()));
  fhX_dev_full[imatch]->Fill(hit3->GetX() - linex->Eval(hit3->GetZ()));
  fhX_dev_full[imatch]->Fill(hit4->GetX() - linex->Eval(hit4->GetZ()));

  fhX_dev_Vds1[imatch]->Fill(hit1->GetX() - linex->Eval(hit1->GetZ()));
  fhX_dev_Vds2[imatch]->Fill(hit2->GetX() - linex->Eval(hit2->GetZ()));
  fhX_dev_Vds3[imatch]->Fill(hit3->GetX() - linex->Eval(hit3->GetZ()));
  fhX_dev_Vds4[imatch]->Fill(hit4->GetX() - linex->Eval(hit4->GetZ()));

  fhY_dev_full[imatch]->Fill(hit1->GetY() - liney->Eval(hit1->GetZ()));
  fhY_dev_full[imatch]->Fill(hit2->GetY() - liney->Eval(hit2->GetZ()));
  fhY_dev_full[imatch]->Fill(hit3->GetY() - liney->Eval(hit3->GetZ()));
  fhY_dev_full[imatch]->Fill(hit4->GetY() - liney->Eval(hit4->GetZ()));

  fhY_dev_Vds1[imatch]->Fill(hit1->GetY() - liney->Eval(hit1->GetZ()));
  fhY_dev_Vds2[imatch]->Fill(hit2->GetY() - liney->Eval(hit2->GetZ()));
  fhY_dev_Vds3[imatch]->Fill(hit3->GetY() - liney->Eval(hit3->GetZ()));
  fhY_dev_Vds4[imatch]->Fill(hit4->GetY() - liney->Eval(hit4->GetZ()));

  delete gr_x;
  delete gr_y;
  delete linex;
  delete liney;
}

//_____________________________________________________________
void Na61VdTrackingModule::CreateFullTrack_new(int imatch, UVdTrack* track1, UVdTrack* track2, UDataTable* tracktab) {
  USensorHit* hit1 = (USensorHit*)track1->GetHitIdAtStation(0);
  USensorHit* hit2 = (USensorHit*)track2->GetHitIdAtStation(1);
  USensorHit* hit3 = (USensorHit*)track2->GetHitIdAtStation(2);
  USensorHit* hit4 = (USensorHit*)track2->GetHitIdAtStation(3);
  ((USensorHit*)track1->GetHitIdAtStation(1))->SetUniqueID(1);
  ((USensorHit*)track1->GetHitIdAtStation(2))->SetUniqueID(1);

  int ind1 = track1->GetHitIndexOnStation(0);
  int ind2 = track2->GetHitIndexOnStation(1);
  int ind3 = track2->GetHitIndexOnStation(2);
  int ind4 = track2->GetHitIndexOnStation(3);

  int tind1 = track1->GetTabIndexOnStation(0);
  int tind2 = track2->GetTabIndexOnStation(1);
  int tind3 = track2->GetTabIndexOnStation(2);
  int tind4 = track2->GetTabIndexOnStation(3);

  if (track1->GetHitIdAtStation(1) != track2->GetHitIdAtStation(1)) cout << " Na61VdTrackingModule::CreateFullTrack: hits ID differnt on Vds2 station. In principal should not happened" << endl;

  // if(track1->GetHitIdAtStation(2) != track2->GetHitIdAtStation(2))
  // cout<<" Na61VdTrackingModule::CreateFullTrack: hits ID differnt on Vds3 station. In principal should not happened"<<endl;

  TF1* linex = new TF1("linex", "[0] +[1]*x", -10.0, 160.0);
  TF1* liney = new TF1("liney", "[0] +[1]*x", -10.0, 160.0);

  double ax;
  double ay;
  double bx;
  double by;
  double chi2x;
  double chi2y;
  double N;
  double zmin;

  double sig[4] = {0.004, 0.0145, 0.032, 0.0416};  // errors in mm
  USensorHit* hits[4] = {hit1, hit2, hit3, hit4};

  FitLine_w2(hits, sig, ax, ay, bx, by, chi2x, chi2y, N, zmin);

  // fhChi2Ndf_x->Fill(chi2x/N);
  // fhChi2Ndf_y->Fill(chi2y/N);

  linex->SetParameter(0, bx);
  linex->SetParameter(1, ax);

  liney->SetParameter(0, by);
  liney->SetParameter(1, ay);

  double tan_x = ax;
  double tan_y = ay;

  if (fhAx_full[imatch]) fhAx_full[imatch]->Fill(tan_x);
  if (fhAy_full[imatch]) fhAy_full[imatch]->Fill(tan_y);
  if (fJura && tan_x > 0.007 && fhAy_full[imatch]) fhAy_full_prod[imatch]->Fill(tan_y);
  if (!fJura && tan_x < 0.005 && fhAy_full[imatch]) fhAy_full_prod[imatch]->Fill(tan_y);

  Vector3D origin(bx, by, 0.);
  Vector3D direction(ax, ay, 1.);
  UVdTrack* track = new UVdTrack(origin, direction);
  tracktab->Add(track);
  track->SetTrackID((long long)track);
  track->SetPdgID(0);
  track->MarkForRemoval(false);

  if (tan_x < fParams->GetAxCut() && fJura) track->MarkForRemoval(true);
  if (tan_x > fParams->GetAxCut() && !fJura) track->MarkForRemoval(true);  // take only third peak in x

  // cout<<fParams->GetAxCut()<<endl;

  track->SetVdsHitIDs(hit1, hit2, hit3, hit4);
  track->SetHitArrayIndexes(ind1, ind2, ind3, ind4);
  track->SetTabArrayIndexes(tind1, tind2, tind3, tind4);
  hit1->SetUniqueID(1);
  hit2->SetUniqueID(1);
  hit3->SetUniqueID(1);
  hit4->SetUniqueID(1);

  // make deviation diagnostics
  fhX_dev_full[imatch]->Fill(hit1->GetX() - linex->Eval(hit1->GetZ()));
  fhX_dev_full[imatch]->Fill(hit2->GetX() - linex->Eval(hit2->GetZ()));
  fhX_dev_full[imatch]->Fill(hit3->GetX() - linex->Eval(hit3->GetZ()));
  fhX_dev_full[imatch]->Fill(hit4->GetX() - linex->Eval(hit4->GetZ()));

  fhX_dev_Vds1[imatch]->Fill(hit1->GetX() - linex->Eval(hit1->GetZ()));
  fhX_dev_Vds2[imatch]->Fill(hit2->GetX() - linex->Eval(hit2->GetZ()));
  fhX_dev_Vds3[imatch]->Fill(hit3->GetX() - linex->Eval(hit3->GetZ()));
  fhX_dev_Vds4[imatch]->Fill(hit4->GetX() - linex->Eval(hit4->GetZ()));

  fhY_dev_full[imatch]->Fill(hit1->GetY() - liney->Eval(hit1->GetZ()));
  fhY_dev_full[imatch]->Fill(hit2->GetY() - liney->Eval(hit2->GetZ()));
  fhY_dev_full[imatch]->Fill(hit3->GetY() - liney->Eval(hit3->GetZ()));
  fhY_dev_full[imatch]->Fill(hit4->GetY() - liney->Eval(hit4->GetZ()));

  fhY_dev_Vds1[imatch]->Fill(hit1->GetY() - liney->Eval(hit1->GetZ()));
  fhY_dev_Vds2[imatch]->Fill(hit2->GetY() - liney->Eval(hit2->GetZ()));
  fhY_dev_Vds3[imatch]->Fill(hit3->GetY() - liney->Eval(hit3->GetZ()));
  fhY_dev_Vds4[imatch]->Fill(hit4->GetY() - liney->Eval(hit4->GetZ()));

  // delete  gr_x;
  // delete  gr_y;
  delete linex;
  delete liney;
}

//_____________________________________________________________
void Na61VdTrackingModule::CreateFullTrackx_new(int imatch, UVdTrack* track, UDataTable* tracktab) {
  USensorHit* hit1 = (USensorHit*)track->GetHitIdAtStation(0);
  USensorHit* hit2 = (USensorHit*)track->GetHitIdAtStation(1);
  USensorHit* hit3 = (USensorHit*)track->GetHitIdAtStation(2);
  USensorHit* hit4 = (USensorHit*)track->GetHitIdAtStation(3);

  int ind1 = track->GetHitIndexOnStation(0);
  int ind2 = track->GetHitIndexOnStation(1);
  int ind3 = track->GetHitIndexOnStation(2);
  int ind4 = track->GetHitIndexOnStation(3);

  int tind1 = track->GetTabIndexOnStation(0);
  int tind2 = track->GetTabIndexOnStation(1);
  int tind3 = track->GetTabIndexOnStation(2);
  int tind4 = track->GetTabIndexOnStation(3);

  TF1* linex = new TF1("linex", "[0] +[1]*x", -10.0, 160.0);
  TF1* liney = new TF1("liney", "[0] +[1]*x", -10.0, 160.0);

  double ax;
  double ay;
  double bx;
  double by;
  double chi2x;
  double chi2y;
  double N;
  double zmin;

  double sig[4] = {0.004, 0.0145, 0.032, 0.0416};  // errors in mm
  USensorHit* hits[4] = {hit1, hit2, hit3, hit4};

  FitLine_w2(hits, sig, ax, ay, bx, by, chi2x, chi2y, N, zmin);

  linex->SetParameter(0, bx);
  linex->SetParameter(1, ax);

  liney->SetParameter(0, by);
  liney->SetParameter(1, ay);

  double tan_x = ax;
  double tan_y = ay;

  if (fhAx_full[imatch]) fhAx_full[imatch]->Fill(tan_x);
  if (fhAy_full[imatch]) fhAy_full[imatch]->Fill(tan_y);
  if (fJura && tan_x > 0.007 && fhAy_full[imatch]) fhAy_full_prod[imatch]->Fill(tan_y);
  if (!fJura && tan_x < 0.005 && fhAy_full[imatch]) fhAy_full_prod[imatch]->Fill(tan_y);

  Vector3D origin(bx, by, 0.);
  Vector3D direction(ax, ay, 1.);
  UVdTrack* trackx = new UVdTrack(origin, direction);

  tracktab->Add(trackx);

  trackx->SetTrackID((long long)trackx);
  trackx->SetPdgID(0);
  trackx->MarkForRemoval(false);

  if (tan_x < 0.007 && fJura) trackx->MarkForRemoval(true);
  if (tan_x > 0.005 && !fJura) trackx->MarkForRemoval(true);

  trackx->SetVdsHitIDs(hit1, hit2, hit3, hit4);
  trackx->SetHitArrayIndexes(ind1, ind2, ind3, ind4);
  trackx->SetTabArrayIndexes(tind1, tind2, tind3, tind4);
  if (hit1) hit1->SetUniqueID(1);
  if (hit2) hit2->SetUniqueID(1);
  if (hit3) hit3->SetUniqueID(1);
  if (hit4) hit4->SetUniqueID(1);

  // make deviation diagnostics
  if (hit1) fhX_dev_full[imatch]->Fill(hit1->GetX() - linex->Eval(hit1->GetZ()));
  if (hit2) fhX_dev_full[imatch]->Fill(hit2->GetX() - linex->Eval(hit2->GetZ()));
  if (hit3) fhX_dev_full[imatch]->Fill(hit3->GetX() - linex->Eval(hit3->GetZ()));
  if (hit4) fhX_dev_full[imatch]->Fill(hit4->GetX() - linex->Eval(hit4->GetZ()));

  if (hit1) fhX_dev_Vds1[imatch]->Fill(hit1->GetX() - linex->Eval(hit1->GetZ()));
  if (hit2) fhX_dev_Vds2[imatch]->Fill(hit2->GetX() - linex->Eval(hit2->GetZ()));
  if (hit3) fhX_dev_Vds3[imatch]->Fill(hit3->GetX() - linex->Eval(hit3->GetZ()));
  if (hit4) fhX_dev_Vds4[imatch]->Fill(hit4->GetX() - linex->Eval(hit4->GetZ()));

  if (hit1) fhY_dev_full[imatch]->Fill(hit1->GetY() - liney->Eval(hit1->GetZ()));
  if (hit2) fhY_dev_full[imatch]->Fill(hit2->GetY() - liney->Eval(hit2->GetZ()));
  if (hit3) fhY_dev_full[imatch]->Fill(hit3->GetY() - liney->Eval(hit3->GetZ()));
  if (hit4) fhY_dev_full[imatch]->Fill(hit4->GetY() - liney->Eval(hit4->GetZ()));

  if (hit1) fhY_dev_Vds1[imatch]->Fill(hit1->GetY() - liney->Eval(hit1->GetZ()));
  if (hit2) fhY_dev_Vds2[imatch]->Fill(hit2->GetY() - liney->Eval(hit2->GetZ()));
  if (hit3) fhY_dev_Vds3[imatch]->Fill(hit3->GetY() - liney->Eval(hit3->GetZ()));
  if (hit4) fhY_dev_Vds4[imatch]->Fill(hit4->GetY() - liney->Eval(hit4->GetZ()));

  delete linex;
  delete liney;

  // calculate deviations for accepted 3hit tracks
  if (!hit1) FillDeviations(imatch, hit2, hit3, hit4);
  if (!hit2) FillDeviations(imatch, hit1, hit3, hit4);
  if (!hit3) FillDeviations(imatch, hit1, hit2, hit4);
  if (!hit4) FillDeviations(imatch, hit1, hit2, hit3);
}

//_____________________________________________________________________________________________
void Na61VdTrackingModule::FillDeviations(int imatch, USensorHit* hit1, USensorHit* hit2, USensorHit* hit3) {
  double x1 = hit1->GetX();
  double y1 = hit1->GetY();
  double z1 = hit1->GetZ();

  double x2 = hit2->GetX();
  double y2 = hit2->GetY();
  double z2 = hit2->GetZ();

  double x3 = hit3->GetX();
  double y3 = hit3->GetY();
  double z3 = hit3->GetZ();

  double devx = ((z2 - z1) * x3 + (z3 - z2) * x1) / (z3 - z1) - x2;
  double devy = ((z2 - z1) * y3 + (z3 - z2) * y1) / (z3 - z1) - y2;

  if (fhX_dev_acce[imatch]) fhX_dev_acce[imatch]->Fill(devx);
  if (fhY_dev_acce[imatch]) fhY_dev_acce[imatch]->Fill(devy);
}

//___________________________________________________________________________________
void Na61VdTrackingModule::FitLine_w2(USensorHit** hits, double* sig, double& ax, double& ay, double& bx, double& by, double& chi2x, double& chi2y, double& N, double& zmin) {
  // simple regration (assuming constant errors)
  // It is quite straight forward to use weighted regression

  double hx, hy, hz;

  double S = 0;
  double Sz = 0;
  double Szz = 0;

  double Sx = 0;
  double Szx = 0;

  double Sy = 0;
  double Szy = 0;

  // cout<<hit1->GetZ()<<" "<<hit2->GetZ()<<" "<<hit3->GetZ()<<" "<<hit4->GetZ()<<endl;

  zmin = 10000.;
  N = 0;

  // perhaps we need waited regression

  for (int i = 0; i < 4; i++) {
    USensorHit* hit = hits[i];

    if (!hit) continue;  // now we always have 4 hits

    N++;

    // smearing is done in dedicated method
    hx = hit->GetX();
    hy = hit->GetY();
    hz = hit->GetZ();

    // int   ivds = hit->GetStationID()-1;
    double sig2 = sig[i] * sig[i];

    S += 1. / sig2;
    Sz += hz / sig2;
    Szz += (hz * hz) / sig2;

    Sx += hx / sig2;
    Szx += (hz * hx) / sig2;

    Sy += hy / sig2;
    Szy += (hz * hy) / sig2;

    if (hz < zmin) zmin = hz;
  }

  ax = (S * Szx - Sz * Sx) / (S * Szz - Sz * Sz);
  bx = (Szz * Sx - Sz * Szx) / (S * Szz - Sz * Sz);

  ay = (S * Szy - Sz * Sy) / (S * Szz - Sz * Sz);
  by = (Szz * Sy - Sz * Szy) / (S * Szz - Sz * Sz);

  chi2x = 0;
  chi2y = 0;
  for (int i = 0; i < 4; i++) {
    USensorHit* hit = hits[i];

    if (!hit) continue;

    // int   ivds = hit->GetStationID();
    double sig2 = sig[i] * sig[i];

    hx = hit->GetX();
    hy = hit->GetY();
    hz = hit->GetZ();

    chi2x += (hx - ax * hz - bx) * (hx - ax * hz - bx) / sig2;
    chi2y += (hy - ay * hz - by) * (hy - ay * hz - by) / sig2;
  }
}

//____________________________________________________________________________________
void Na61VdTrackingModule::Create3HitTrack(int idev, USensorHit* hit1, USensorHit* hit2, USensorHit* hit3, int* IndArray, int* TabIndArray, UDataTable* tracktab) {
  TF1* line = new TF1("line", "[0] +[1]*x", -10.0, 160.0);

  double x[] = {hit1->GetZ(), hit2->GetZ(), hit3->GetZ()};
  double yx[] = {hit1->GetX(), hit2->GetX(), hit3->GetX()};
  double ex[] = {0, 0, 0};
  double eyx[] = {0.004, 0.004, 0.004};

  double yy[] = {hit1->GetY(), hit2->GetY(), hit3->GetY()};
  double eyy[] = {0.004, 0.004, 0.004};

  TGraphErrors* gr_x = new TGraphErrors(3, x, yx, ex, eyx);

  TGraphErrors* gr_y = new TGraphErrors(3, x, yy, ex, eyy);

  line->SetParameter(0, yx[0]);
  gr_x->Fit("line", "Q", "", x[0] - 1, x[2] + 1);

  fhChi2X[idev]->Fill(line->GetChisquare());

  double x_Vds1 = line->Eval(0.);

  double tan_x = line->GetParameter(1);

  line->SetParameter(0, yy[0]);
  gr_y->Fit("line", "Q", "", x[0] - 1, x[2] + 1);

  fhChi2Y[idev]->Fill(line->GetChisquare());

  double y_Vds1 = line->Eval(0.);

  double tan_y = line->GetParameter(1);

  if (fLargeClusters) {
    if (fhAx_lc[idev]) fhAx_lc[idev]->Fill(tan_x);
    if (fhAy_lc[idev]) fhAy_lc[idev]->Fill(tan_y);
  } else {
    if (fhAx[idev]) fhAx[idev]->Fill(tan_x);
    if (fhAy[idev]) fhAy[idev]->Fill(tan_y);
  }

  Vector3D origin(x_Vds1, y_Vds1, 0.);
  Vector3D direction(tan_x, tan_y, 1.);
  UVdTrack* track = new UVdTrack(origin, direction);

  tracktab->Add(track);

  track->SetTrackID(fTrackID++);
  track->SetPdgID(0);

  // cout<<" create tabIdArrays: "<<TabIndArray[0]<<"  "<<TabIndArray[1]<<" "<<TabIndArray[2]<<endl;

  if (idev < 6) {
    if (fDevNames[idev].Contains("dev1_")) {
      track->SetVdsHitIDs(hit1, hit2, hit3, 0);
      track->SetHitArrayIndexes(IndArray[0], IndArray[1], IndArray[2], -1);
      track->SetTabArrayIndexes(TabIndArray[0], TabIndArray[1], TabIndArray[2], -1);
    }
    if (fDevNames[idev].Contains("dev2_")) {
      track->SetVdsHitIDs(0, hit1, hit2, hit3);
      track->SetHitArrayIndexes(-1, IndArray[0], IndArray[1], IndArray[2]);
      track->SetTabArrayIndexes(-1, TabIndArray[0], TabIndArray[1], TabIndArray[2]);
    }
  } else {
    if (fDevNames[idev].Contains("x1")) {
      track->SetVdsHitIDs(0, hit1, hit2, hit3);
      track->SetHitArrayIndexes(-1, IndArray[0], IndArray[1], IndArray[2]);
      track->SetTabArrayIndexes(-1, TabIndArray[0], TabIndArray[1], TabIndArray[2]);
    }
    if (fDevNames[idev].Contains("x2")) {
      track->SetVdsHitIDs(hit1, 0, hit2, hit3);
      track->SetHitArrayIndexes(IndArray[0], -1, IndArray[1], IndArray[2]);
      track->SetTabArrayIndexes(TabIndArray[0], -1, TabIndArray[1], TabIndArray[2]);
    }
    if (fDevNames[idev].Contains("x3")) {
      track->SetVdsHitIDs(hit1, hit2, 0, hit3);
      track->SetHitArrayIndexes(IndArray[0], IndArray[1], -1, IndArray[2]);
      track->SetTabArrayIndexes(TabIndArray[0], TabIndArray[1], -1, TabIndArray[2]);
    }
    if (fDevNames[idev].Contains("x4")) {
      track->SetVdsHitIDs(hit1, hit2, hit3, 0);
      track->SetHitArrayIndexes(IndArray[0], IndArray[1], IndArray[2], -1);
      track->SetTabArrayIndexes(TabIndArray[0], TabIndArray[1], TabIndArray[2], -1);
    }
  }
  // cout<<fDevNames[idev].Data()<<" "<<track->GetHitIdAtStation(0)<<" "<<track->GetHitIdAtStation(1)<<" "
  //  <<track->GetHitIdAtStation(2)<<" "<<track->GetHitIdAtStation(3)<<" "<<hit1->GetZ()<<" "<<hit2->GetZ()<<" "<<hit3->GetZ()<<endl;

  delete gr_x;
  delete gr_y;
  delete line;
}

//_____________________________________________________________
void Na61VdTrackingModule::CheckTrackMatching(UVdTrack* track1, UVdTrack* track2) {
  // TWO or more same hits in the track leads to discarding one track

  int isame = 0;
  for (int i = 0; i < 4; i++) {
    if ((track1->GetHitIdAtStation(i) == track2->GetHitIdAtStation(i)) && track1->GetHitIdAtStation(i)) isame++;
  }

  if (isame < 2) return;

  // double ang1_x = track1->GetDX()/track1->GetDZ();
  // double ang1_y = track1->GetDY()/track1->GetDZ();
  // double ang2_x = track2->GetDX()/track2->GetDZ();
  // double ang2_y = track2->GetDY()/track2->GetDZ();

  // make one track for removal
  if (track1->GetChi2Ndf() > track2->GetChi2Ndf()) {
    // if(track1->IsMarkedForRemoval())cout<<"track1 already marked"<<endl;
    track1->MarkForRemoval(true);  // cout<<"marked track1: "<<track1<<endl;
  } else {
    // if(track2->IsMarkedForRemoval())cout<<"track2 already marked"<<endl;
    track2->MarkForRemoval(true);  // cout<<"marked track2: "<<track2<<endl;
  }
}

//_____________________________________________________________
void Na61VdTrackingModule::CheckTrackMatching2(UVdTrack* track1, UVdTrack* track2) {
  // TWO or more same hits in the track leads to discarding one track

  int isame = 0;
  for (int i = 0; i < 4; i++) {
    if ((track1->GetHitIdAtStation(i) == track2->GetHitIdAtStation(i)) && track1->GetHitIdAtStation(i)) isame++;
    // cout<<"i = "<<i<<"   "<<track1->GetHitIdAtStation(i)<<"   "<<track2->GetHitIdAtStation(i)<<endl;
  }

  if (isame < 2) return;

  // double ang1_x = track1->GetDX()/track1->GetDZ();
  // double ang1_y = track1->GetDY()/track1->GetDZ();
  // double ang2_x = track2->GetDX()/track2->GetDZ();
  // double ang2_y = track2->GetDY()/track2->GetDZ();

  // make one track for removal
  if (track1->GetChi2Ndf() > track2->GetChi2Ndf()) {
    // if(track1->IsMarkedForRemoval())cout<<"track1 already marked"<<endl;
    track1->MarkForRemoval(true);  // cout<<"marked track1: "<<track1<<endl;
  } else {
    // if(track2->IsMarkedForRemoval())cout<<"track2 already marked"<<endl;
    track2->MarkForRemoval(true);  // cout<<"marked track2: "<<track2<<endl;
  }
}

//_____________________________________________________________
void Na61VdTrackingModule::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61VdTrackingModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);
}

//____________________________________________________________________
void Na61VdTrackingModule::Print(Option_t* option) const {
  // Print module informationd
  // In addition this module defines the Option:
  // <fill in here>

  TString opt(option);
  opt.ToLower();

  Na61Module::Print(option);
  if (opt.Contains("d")) cout << endl << "  Original author: Paweł Staszel" << endl << "  Last Modifications: " << endl << "    $Author: Staszel $" << endl << "    $Date: 2016/10/30$" << endl << "    $Revision: 1.0 $ " << endl << endl << "-------------------------------------------------" << endl;
}
