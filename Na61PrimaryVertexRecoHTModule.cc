//--------------------------------------------
// Primary vertex reco module
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61PrimaryVertexRecoHTModule
#include "Na61PrimaryVertexRecoHTModule.h"
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
// ClassImp(Na61PrimaryVertexRecoHTModule);

//____________________________________________________________________
Na61PrimaryVertexRecoHTModule::Na61PrimaryVertexRecoHTModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
  fZprim = 47.0;
  fZCut = 3;
  fTrackInput = 0;
}

//____________________________________________________________________
Na61PrimaryVertexRecoHTModule::Na61PrimaryVertexRecoHTModule(const char* name, const char* title) : Na61PrimaryVertexRecoModule(name, title) {
  // Named Constructor
  SetState(kSetup);
  fZprim = 47.0;
  fZCut = 3;
  fTrackInput = 0;
  fSigx = 0.12;
  fSigy = 0.07;
  }

//____________________________________________________________________
void Na61PrimaryVertexRecoHTModule::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  TDirectory* histDir = 0;
  if (fTrackInput == 0) {
    histDir = gDirectory->mkdir("PrimaryVertexRecoHT");
  } else if (fTrackInput == 1) {
    histDir = gDirectory->mkdir("PrimaryVertexRecoTPC");
  } else if (fTrackInput == 2) {
    histDir = gDirectory->mkdir("PrimaryVertexRecoVdCa");
  } else {
    cout << "Na61PrimaryVertexRecoHTModule::DefineHistograms(): Wrong fTrackInput, Check it out!!!!" << endl;
  }

  histDir->cd();

  fhVtxX = new TH1F("hVtxX", "", 300, -3., 3);
  fhVtxY = new TH1F("hVtxY", "", 300, -3., 3);
  fhTracksVsVtxX_J = new TH1F("hTracksVsVtxX_J", "", 30, -3., 3);
  fhTracksVsVtxY_J = new TH1F("hTracksVsVtxY_J", "", 30, -3., 3);
  fhTracksVsVtxX_S = new TH1F("hTracksVsVtxX_S", "", 30, -3., 3);
  fhTracksVsVtxY_S = new TH1F("hTracksVsVtxY_S", "", 30, -3., 3);

  fhVtxDx = new TH1F("hVtxDx", "", 1000, -2.5, 2.5);
  fhVtxDy = new TH1F("hVtxDy", "", 1000, -2.5, 2.5);
  fhVtxDz = new TH1F("hVtxDz", "", 1000, -10., 10.);

  // using Jura versus Saleve arm
  fhVtxDx_Glob = new TH1F("hVtxDx_Glob", "", 1000, -2.5, 2.5);
  fhVtxDy_Glob = new TH1F("hVtxDy_Glob", "", 1000, -2.5, 2.5);
  fhVtxDz_Glob = new TH1F("hVtxDz_Glob", "", 1000, -10., 10.);

  fhVtxDx_fine_Glob = new TH1F("hVtxDx_fine_Glob", "", 1000, -0.25, 0.25);
  fhVtxDy_fine_Glob = new TH1F("hVtxDy_fine_Glob", "", 1000, -0.25, 0.25);
  fhVtxDz_fine_Glob = new TH1F("hVtxDz_fine_Glob", "", 1000, -2.5, 2.5);

  fhVtxStatus_JuraVsSaleve = new TH2F("VtxStatus_JuraVsSaleve", "", 5, 0, 5, 5, 0, 5);
  fhFullTracks_JvsS = new TH2F("hFullTracks_JvsS", "", 200, 0, 200, 200, 0, 200);
  fhFullTracks_JS = new TH1F("hFullTracks_JS", "", 300, 0, 300);
  fhFullTracks_J = new TH1F("hFullTracks_J", "", 300, 0, 300);
  fhFullTracks_S = new TH1F("hFullTracks_S", "", 300, 0, 300);
  fhFullTracks4h_JS = new TH1F("hFullTracks4h_JS", "", 300, 0, 300);
  fhFullTracks4h_J = new TH1F("hFullTracks4h_J", "", 300, 0, 300);
  fhFullTracks4h_S = new TH1F("hFullTracks4h_S", "", 300, 0, 300);

  double offset = 0;
  if (fTrackInput == 1) offset = -550;
  if (fTrackInput == 2) offset = -550;

  fhRecoVertexZ = new TH1F(Form("hRecoVertexZ"), " ", 2000, -500. + offset, 100. + offset);
  fhRecoVertexZ_fine = new TH1F(Form("hRecoVertexZ_fine"), " ", 1000, -70. + offset, -30. + offset);
  fhRecoVertexZ_xfine = new TH1F(Form("hRecoVertexZ_xfine"), " ", 1000, -50. + offset, -44. + offset);
  fhRecoVertexXY = new TH2F(Form("hRecoVertexXY"), " ", 500, -50., 50., 500, -50., 50.);

  fhRecoVertexZ_J = new TH1F(Form("hRecoVertexZ_J"), " ", 2000, -500. + offset, 100. + offset);
  fhRecoVertexZ_fine_J = new TH1F(Form("hRecoVertexZ_fine_J"), " ", 1000, -70. + offset, -30. + offset);
  fhRecoVertexXY_J = new TH2F(Form("hRecoVertexXY_J"), " ", 500, -50., 50., 500, -50., 50.);

  fhRecoVertexZ_S = new TH1F(Form("hRecoVertexZ_S"), " ", 2000, -500. + offset, 100. + offset);
  fhRecoVertexZ_fine_S = new TH1F(Form("hRecoVertexZ_fine_S"), " ", 1000, -70. + offset, -30. + offset);
  fhRecoVertexXY_S = new TH2F(Form("hRecoVertexXY_S"), " ", 500, -50., 50., 500, -50., 50.);

  fhAxJ = new TH1F(Form("hAxJ"), "", 1500, -0.25, 0.25);
  fhAyJ = new TH1F(Form("hAyJ"), "", 1500, -0.25, 0.25);
  fhAxS = new TH1F(Form("hAxS"), "", 1500, -0.25, 0.25);
  fhAyS = new TH1F(Form("hAyS"), "", 1500, -0.25, 0.25);

  for (int i = 0; i < 2; i++) {
    // using Sub1 versus Sub2
    fhVtxDx_Glob2[i] = new TH1F(Form("hVtxDx_Glob2_%d", i), "", 1000, -1.5, 1.5);
    fhVtxDy_Glob2[i] = new TH1F(Form("hVtxDy_Glob2_%d", i), "", 1000, -1.0, 1.0);
    fhVtxDz_Glob2[i] = new TH1F(Form("hVtxDz_Glob2_%d", i), "", 1000, -10., 10.);

    fhVtxDx_fine_Glob2[i] = new TH1F(Form("hVtxDx_fine_Glob2_%d", i), "", 1000, -0.25, 0.25);
    fhVtxDy_fine_Glob2[i] = new TH1F(Form("hVtxDy_fine_Glob2_%d", i), "", 1000, -0.25, 0.25);
    fhVtxDz_fine_Glob2[i] = new TH1F(Form("hVtxDz_fine_Glob2_%d", i), "", 1000, -2.5, 2.5);

    if (i > 2) continue;
    fhTrackDx[i] = new TH1F(Form("hTrackDx_%d", i), "", 1000, -.5, .5);
    fhTrackDy[i] = new TH1F(Form("hTrackDy_%d", i), "", 1000, -.5, .5);
    fhTrackDz[i] = new TH1F(Form("hTrackDz_%d", i), "", 1000, -.5, .5);

    fhTrackDx_tag[i] = new TH1F(Form("hTrackDx_tag_%d", i), "", 1000, -.5, .5);
    fhTrackDy_tag[i] = new TH1F(Form("hTrackDy_tag_%d", i), "", 1000, -.5, .5);
    fhTrackDz_tag[i] = new TH1F(Form("hTrackDz_tag_%d", i), "", 1000, -.5, .5);
  }

  fhVtxDx_fine_vs_multi = new TH2F(Form("hVtxDx_fine_vs_multi"), "", 100, 0, 300, 500, -0.12, 0.12);
  fhVtxDy_fine_vs_multi = new TH2F(Form("hVtxDy_fine_vs_multi"), "", 100, 0, 300, 500, -0.12, 0.12);
  fhVtxDz_fine_vs_multi = new TH2F(Form("hVtxDz_fine_vs_multi"), "", 100, 0, 300, 500, -1.2, 1.2);

  /* was used to calculate sesondary vertex resolution for Marek
  TString str[3] = {"La1","La2","La3"};
  for(int i=0;i<3;i++){
    fhVtxDx_fine_vs_multi_La[i] = new TH2F(Form("hVtxDx_fine_vs_multi_%s",str[i].Data()),"", 100,0,300, 500,-0.12,0.12);
    fhVtxDy_fine_vs_multi_La[i] = new TH2F(Form("hVtxDy_fine_vs_multi_%s",str[i].Data()),"", 100,0,300, 500,-0.12,0.12);
    fhVtxDz_fine_vs_multi_La[i] = new TH2F(Form("hVtxDz_fine_vs_multi_%s",str[i].Data()),"", 100,0,300, 500,-1.2,1.2);
  }
  */

  fhDx_2TrackVertexRes = new TH2F(Form("hDx_2TrackVertexRes"), "", 100, 0, 20, 500, -0.15, 0.15);
  fhDy_2TrackVertexRes = new TH2F(Form("hDy_2TrackVertexRes"), "", 100, 0, 20, 500, -0.15, 0.15);
  fhDz_2TrackVertexRes = new TH2F(Form("hDz_2TrackVertexRes"), "", 100, 0, 20, 500, -1.5, 1.5);

  gDirectory->cd("..");
}

//____________________________________________________________________
void Na61PrimaryVertexRecoHTModule::Init() {
  // Job-level initialisation
  SetState(kInit);
  Na61PrimaryVertexRecoModule::Init();
  
  if(fTrackInput==1){
    // this values should be updated later
    //fSigx = 0.1939;
    //fSigy = 0.2082;
    //fSigx = 0.12;
    //fSigy = 0.07;
    //fSigx = 1.0e9; // only y counts
    //fSigy = 1.;
    fSigx = 1.707; // only x counts
    fSigy = 1.302;

  }

}

//____________________________________________________________________
void Na61PrimaryVertexRecoHTModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);
  Na61PrimaryVertexRecoModule::Begin();
}

//____________________________________________________________________
void Na61PrimaryVertexRecoHTModule::Event(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);

  double pvtxX_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexX();
  double pvtxY_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexY();
  double pvtxZ_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexZ();

  double pvtxX_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexX();
  double pvtxY_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexY();
  double pvtxZ_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexZ();

  int pvertStatusJ = ((UVdEvent*)inNodeJ)->GetPrimaryVertexStatus();
  int pvertStatusS = ((UVdEvent*)inNodeS)->GetPrimaryVertexStatus();

  fhVtxStatus_JuraVsSaleve->Fill(pvertStatusJ, pvertStatusS);

  //  int pvertStatus = ((UVdEvent*)outNode)->GetPrimaryVertexStatus();
  // double pvtxZ = ((UVdEvent*)outNode)->GetPrimaryVertexZ();

   //cout<<" pvertStatusJ="<<pvertStatusJ<<" pvertStatusS="<<pvertStatusS<<endl;;//<<"  pvertStatus="<<pvertStatus<<"  z="<<pvtxZ<<endl;

   //cout<<"###############################   new event: pvtxZ_J="<<pvtxZ_J<<"  pvtxZ_S="<<pvtxZ_S<<endl;

  if ((pvertStatusJ > 0) && (pvertStatusS > 0)) {
    fhVtxDx->Fill(pvtxX_J - pvtxX_S);
    fhVtxDy->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz->Fill(pvtxZ_J - pvtxZ_S);
  }

  fPrimaryVertex[0].SetX(11111);
  fPrimaryVertex[0].SetY(11111);
  fPrimaryVertex[0].SetZ(11111);
  fPrimaryVertex[1].SetX(11111);
  fPrimaryVertex[1].SetY(11111);
  fPrimaryVertex[1].SetZ(11111);
  fPrimaryVertex[2].SetX(11111);
  fPrimaryVertex[2].SetY(11111);
  fPrimaryVertex[2].SetZ(11111);

  MethodHT(outNode);

  pvtxX_J = fPrimaryVertexJ.GetX();
  pvtxY_J = fPrimaryVertexJ.GetY();
  pvtxZ_J = fPrimaryVertexJ.GetZ();

  pvtxX_S = fPrimaryVertexS.GetX();
  pvtxY_S = fPrimaryVertexS.GetY();
  pvtxZ_S = fPrimaryVertexS.GetZ();

  if ((pvertStatusJ > 0) && (pvertStatusS > 0)) {
    fhVtxDx_Glob->Fill(pvtxX_J - pvtxX_S);
    fhVtxDy_Glob->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz_Glob->Fill(pvtxZ_J - pvtxZ_S);
    fhVtxDx_fine_Glob->Fill(pvtxX_J - pvtxX_S);
    fhVtxDy_fine_Glob->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz_fine_Glob->Fill(pvtxZ_J - pvtxZ_S);
  }

  // set primary vertex to event structure only for fTrackInput=0 (VD tracks)
  if (!(fTrackInput == 0)) return;

  ((UVdEvent*)outNode)->SetPrimaryVertexStatus(0);

  // if(pvertStatusJ || pvertStatusS){
  ((UVdEvent*)outNode)->SetPrimaryVertexStatus(1);
  ((UVdEvent*)outNode)->SetPrimaryVertex(fPrimaryVertex[2].GetX(), fPrimaryVertex[2].GetY(), fPrimaryVertex[2].GetZ());
  SetPrimaryVertex(fPrimaryVertex[2]);
  // to avoid using HT for exotic vertex locations
  if (!(TMath::Abs(fPrimaryVertex[2].GetZ() + fZprim) < fZCut)) ((UVdEvent*)outNode)->SetPrimaryVertexStatus(0);
  //}

  SetPrimaryVertexStatus(((UVdEvent*)outNode)->GetPrimaryVertexStatus());

  cout << " Na61PrimaryVertexRecoHTModule: " << fPrimaryVertex[2].GetX() << " " << fPrimaryVertex[2].GetY() << " " << fPrimaryVertex[2].GetZ() << endl;

  // Find2TrackVertexes(outNode);
}

//____________________________________________________________________
void Na61PrimaryVertexRecoHTModule::Event(UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);

  fPrimaryVertex[0].SetX(11111);
  fPrimaryVertex[0].SetY(11111);
  fPrimaryVertex[0].SetZ(11111);
  fPrimaryVertex[1].SetX(11111);
  fPrimaryVertex[1].SetY(11111);
  fPrimaryVertex[1].SetZ(11111);
  fPrimaryVertex[2].SetX(11111);
  fPrimaryVertex[2].SetY(11111);
  fPrimaryVertex[2].SetZ(11111);

  MethodHT(outNode);

  double pvtxX_J = fPrimaryVertexJ.GetX();
  double pvtxY_J = fPrimaryVertexJ.GetY();
  double pvtxZ_J = fPrimaryVertexJ.GetZ();

  double pvtxX_S = fPrimaryVertexS.GetX();
  double pvtxY_S = fPrimaryVertexS.GetY();
  double pvtxZ_S = fPrimaryVertexS.GetZ();

  // int pvertStatus = ((UVdEvent*)outNode)->GetPrimaryVertexStatus();

  fhVtxDx_Glob->Fill(pvtxX_J - pvtxX_S);
  fhVtxDy_Glob->Fill(pvtxY_J - pvtxY_S);
  fhVtxDz_Glob->Fill(pvtxZ_J - pvtxZ_S);
  fhVtxDx_fine_Glob->Fill(pvtxX_J - pvtxX_S);
  fhVtxDy_fine_Glob->Fill(pvtxY_J - pvtxY_S);
  fhVtxDz_fine_Glob->Fill(pvtxZ_J - pvtxZ_S);

  // set primary vertex to event structure only for fTrackInput=0 (HT Tracks)
  //((UVdEvent*)outNode)->SetPrimaryVertexStatus(0);

  SetPrimaryVertexStatus(1);
  //((UVdEvent*)outNode)->SetPrimaryVertex(fPrimaryVertex[2].GetX(),fPrimaryVertex[2].GetY(),fPrimaryVertex[2].GetZ());
  // to avoid using HT for exotic vertex locations
  SetPrimaryVertex(fPrimaryVertex[2]);
  if (!(TMath::Abs(fPrimaryVertex[2].GetZ() + fZprim) < fZCut)) SetPrimaryVertexStatus(0);

  // Find2TrackVertexes(outNode);
}

//_____________________________________________________________
void Na61PrimaryVertexRecoHTModule::MethodHT(UEventNode* outNode) {
  ClassifyTracksHT(outNode);

  int stage = 0;
    fStage=stage;
  FindPrimaryVertexPostHT(stage, 0, outNode);
  TagPrimaryTracks(stage, outNode);
  CheckPrimaryVertexResolutionHT(stage, outNode);
  // cout<<"stage0, z = "<<fPrimaryVertex[stage].GetZ()<<endl;

  stage = 1;
    fStage=stage;
  FindPrimaryVertexPostHT(stage, 0, outNode);
  TagPrimaryTracks(stage, outNode);
  CheckPrimaryVertexResolutionHT(stage, outNode);
  // cout<<"stage1, z = "<<fPrimaryVertex[stage].GetZ()<<endl;

  stage = 2;
    fStage=stage;
  FindPrimaryVertexPostHT(stage, 0, outNode);
  // cout<<"stage2, z = "<<fPrimaryVertex[stage].GetZ()<<endl;
  FindPrimaryVertexPostHT(stage, 1, outNode);  // for Jura tracks only
  FindPrimaryVertexPostHT(stage, 2, outNode);  // for Saleve tracks only
}

//_____________________________________________________________
void Na61PrimaryVertexRecoHTModule::ClassifyTracksHT(UEventNode* out) {
  // clasification regarding x-slope

  UDataTable* tracktab = 0;
  if (fTrackInput == 0) tracktab = out->GetDataTable(Form("Vd HT Tracks"));
  if (fTrackInput == 1) tracktab = out->GetDataTable(Form("Tpc Tracks"));
  if (fTrackInput == 2) tracktab = out->GetDataTable(Form("VdCa Tracks"));
  if (!tracktab) return;

  // cout<<"tracktab->Entries: = "<<tracktab->GetEntries()<<endl;

  for (int j = 0; j < tracktab->GetEntries(); j++) {
    UVdTrack* track = (UVdTrack*)tracktab->At(j);

    track->SetFlag(2);

    double ax = track->GetDX_f() / track->GetDZ_f();
    double ay = track->GetDY_f() / track->GetDZ_f();
    double x = track->GetX_f();

    if (x > 0) {
      fhAxJ->Fill(ax);
      fhAyJ->Fill(ay);
    } else {
      fhAxS->Fill(ax);
      fhAyS->Fill(ay);
    }
  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoHTModule::TagPrimaryTracks(int stage, UEventNode* out) {
  // clasification regarding x-slope

  UDataTable* tracktab = 0;
  if (fTrackInput == 0) tracktab = out->GetDataTable(Form("Vd HT Tracks"));
  if (fTrackInput == 1) tracktab = out->GetDataTable(Form("Tpc Tracks"));
  if (fTrackInput == 2) tracktab = out->GetDataTable(Form("VdCa Tracks"));

  if (!tracktab) return;

 double sigx  = 0.028;
  double sigy  = 0.0145;
  double sigxx = 0.028;
  double sigyy = 0.0145;
  if(fTrackInput==1){
    sigx  = 0.34; // old values used to p 6x, will be updated later
    sigy  = 0.29;
    sigxx = 0.23; // check distribution of  fhTrackDx
    sigyy = 0.17;
   // sigx = 0.2;
    //sigy = 0.2;
  }

  for (int j = 0; j < tracktab->GetEntries(); j++) {
    UVdTrack* track = (UVdTrack*)tracktab->At(j);

    track->SetTagForVtx(0);
    // track->SetTagForVtx(1);

    // Impact Parameter cut
    Line3D* Line = track->Getlinef();

    Vector3D track_closest_point = Line->ClosestPoint(fPrimaryVertex[stage]);
    Vector3D track_AtCP = track_closest_point - fPrimaryVertex[stage];  // AtCP means At Closest Proximity

    double dx = track_AtCP.X();
    double dy = track_AtCP.Y();
    double dz = track_AtCP.Z();

    fhTrackDx[stage]->Fill(dx);
    fhTrackDy[stage]->Fill(dy);
    fhTrackDz[stage]->Fill(dz);

    double dd = (dx * dx) / (sigx * sigx) + (dy * dy) / (sigy * sigy);

    if (dd < 3 * 3) {
      track->SetTagForVtx(1);
      fhTrackDx_tag[stage]->Fill(dx);
      fhTrackDy_tag[stage]->Fill(dy);
      fhTrackDz_tag[stage]->Fill(dz);
      double dd_xy = (dx*dx)/(sigxx*sigxx) + (dy*dy)/(sigyy*sigyy);
      if(dd_xy < 3*3)track->SetTagForVtx(2);
    }
  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoHTModule::Find2TrackVertexes(UEventNode* out) {
  TObjArray tracks;
  tracks.Clear();

  UDataTable* tracktab = 0;
  if (fTrackInput == 0) tracktab = out->GetDataTable(Form("Vd HT Tracks"));
  if (fTrackInput == 1) tracktab = out->GetDataTable(Form("Tpc Tracks"));
  if (fTrackInput == 2) tracktab = out->GetDataTable(Form("VdCa Tracks"));

  if (!tracktab) return;

  for (int i = 0; i < tracktab->GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracktab->At(i);

    if (!(track->GetTagForVtx())) continue;

    tracks.Add(track);
  }

  if (tracks.GetEntries() < 100) return;

  Vector3D pvtx = fPrimaryVertex[2];

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      // method 1
      double axi = track1->GetDX_f() / track1->GetDZ_f();
      double axj = track2->GetDX_f() / track2->GetDZ_f();
      double bxi = track1->GetX_f();
      double bxj = track2->GetX_f();

      double ayi = track1->GetDY_f() / track1->GetDZ_f();
      double ayj = track2->GetDY_f() / track2->GetDZ_f();
      double byi = track1->GetY_f();
      double byj = track2->GetY_f();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      double z = -(z_x + z_y) / 2.;

      double x1 = track1->GetXatZ_f(z);
      double x2 = track2->GetXatZ_f(z);
      double y1 = track1->GetYatZ_f(z);
      double y2 = track2->GetYatZ_f(z);
      double x = (x1 + x2) / 2.;
      double y = (y1 + y2) / 2.;

      // calculate opening angle
      Line3D* line1 = track1->Getlinef();
      Line3D* line2 = track2->Getlinef();
      const Vector3D& dir1 = line1->GetDirection();
      const Vector3D& dir2 = line2->GetDirection();
      // cout<<"  norm1="<<dir1.Norm()<<"  norm2="<<dir2.Norm()<<endl;
      double cos_theta = dir1.Dot(dir2);
      double opening_angle = TMath::ACos(cos_theta) * TMath::RadToDeg();
      fhDx_2TrackVertexRes->Fill(opening_angle, x - pvtx.GetX());
      fhDy_2TrackVertexRes->Fill(opening_angle, y - pvtx.GetY());
      fhDz_2TrackVertexRes->Fill(opening_angle, z - pvtx.GetZ());
    }
  }
}
//_____________________________________________________________
void Na61PrimaryVertexRecoHTModule::FindPrimaryVertexPostHT(int stage, int flag, UEventNode* out) {
  // algorithm:
  // flag = 0 combined vertex (jura+saleve)
  // flag = 1 jura arm
  // flag = 2 saleve

  int fullTracksJ = 0;
  int fullTracksS = 0;
  //  int fullTracksJ4h = 0;
  //  int fullTracksS4h = 0;
  TObjArray tracks;
  tracks.Clear();

  UDataTable* tracktab = 0;
  if (fTrackInput == 0) tracktab = out->GetDataTable(Form("Vd HT Tracks"));
  if (fTrackInput == 1) tracktab = out->GetDataTable(Form("Tpc Tracks"));
  if (fTrackInput == 2) tracktab = out->GetDataTable(Form("VdCa Tracks"));

  // cout<<" vdHT tracktab: "<<tracktab<<endl;
  if (!tracktab) return;

  for (int i = 0; i < tracktab->GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracktab->At(i);

    if (!(track->GetTagForVtx()) && stage) continue;

    double x_vd = track->GetX_f();
    if (flag == 0) tracks.Add(track);
    if (flag == 1) {
      if (x_vd > 0) tracks.Add(track);
    }
    if (flag == 2) {
      if (x_vd < 0) tracks.Add(track);
    }
  }

  if (tracks.GetEntries() < 2) return;

  // if(flag==0)cout<<"---------------> stage="<<stage<<" tracks="<<tracks.GetEntries()<<endl;

  double sum_up_x = 0;
  double sum_do_x = 0;
  double sum_up_y = 0;
  double sum_do_y = 0;

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      double axi = track1->GetDX_f() / track1->GetDZ_f();
      double axj = track2->GetDX_f() / track2->GetDZ_f();
      double bxi = track1->GetX_f();
      double bxj = track2->GetX_f();

      double ayi = track1->GetDY_f() / track1->GetDZ_f();
      double ayj = track2->GetDY_f() / track2->GetDZ_f();
      double byi = track1->GetY_f();
      double byj = track2->GetY_f();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < fZCut) {
        sum_up_x = sum_up_x + (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - fZprim) < fZCut) {
        sum_up_y = sum_up_y + (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + (ayi - ayj) * (ayi - ayj);
      }

      // if(histId==0 && flag==2)cout<<"xx: "<< ((axi-axj)*(bxi-bxj))/ ((axi-axj)*(axi-axj))<<" yy: "<<((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj))<<endl;
    }
  }

  ////////// new method //////////////////////////////////////////////////
  double zx_min = 0;
  if (sum_do_x > 0) zx_min = -sum_up_x / sum_do_x;
  double zy_min = 0;
  if (sum_do_y > 0) zy_min = -sum_up_y / sum_do_y;

  // cout<<sum_up_x<<" "<<sum_do_x<<"       "<<sum_up_y<<" "<<sum_do_y<<"     "<<zx_min<<" "<<zy_min<<endl;

  // double z_prim = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);
  // double z_prim = -(sum_up_x)/(sum_do_x);
  // double z_prim = -(sum_up_y)/(sum_do_y); // use only  y information

  double sigx = 0.49;
  double sigy = 0.1065;
  // if(flag==0) {sigx=0.12; sigy=0.07;}
  if (flag == 0) {
    sigx = 0.1;
    sigy = 0.1;
  }
  double z_prim = zx_min / (sigx * sigx) + zy_min / (sigy * sigy);
  double norm = 1. / (sigx * sigx) + 1. / (sigy * sigy);
  z_prim = z_prim / norm;

  // cout<<"z_prim="<<z_prim<<endl;
  Vector3D reco_pvertex_w2;
  // PrimaryVertexWithWeigths(tracks,z_prim,reco_pvertex_w2);
  // PrimaryVertexWithWeigths_Y(tracks,z_prim,reco_pvertex_w2);
  // PrimaryVertexWithWeigths_X(tracks,z_prim,reco_pvertex_w2);
  PrimaryVertexWithWeigths_XY(2, tracks, z_prim, reco_pvertex_w2, sigx, sigy);

  if (flag == 0) {  // both arms case
    if (stage == 2) {
      fhVtxX->Fill(reco_pvertex_w2.X());
      fhVtxY->Fill(reco_pvertex_w2.Y());
      fhTracksVsVtxX_J->Fill(reco_pvertex_w2.X(), fullTracksJ);
      fhTracksVsVtxX_S->Fill(reco_pvertex_w2.X(), fullTracksS);
      fhTracksVsVtxY_J->Fill(reco_pvertex_w2.Y(), fullTracksJ);
      fhTracksVsVtxY_S->Fill(reco_pvertex_w2.Y(), fullTracksS);
      fhRecoVertexXY->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
      fhRecoVertexZ->Fill(reco_pvertex_w2.Z());
      fhRecoVertexZ_fine->Fill(reco_pvertex_w2.Z());
      fhRecoVertexZ_xfine->Fill(reco_pvertex_w2.Z());
    }
    // cout<<"stage="<<stage<<" reco_pvertex_w2.Z()="<<reco_pvertex_w2.Z()<<endl;
    fPrimaryVertex[stage].SetX(reco_pvertex_w2.X());
    fPrimaryVertex[stage].SetY(reco_pvertex_w2.Y());
    fPrimaryVertex[stage].SetZ(reco_pvertex_w2.Z());
  }

  // call at the last stage
  if (flag == 1) {  // Jura Arm
    fPrimaryVertexJ.SetX(reco_pvertex_w2.X());
    fPrimaryVertexJ.SetY(reco_pvertex_w2.Y());
    fPrimaryVertexJ.SetZ(reco_pvertex_w2.Z());
    fhRecoVertexXY_J->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
    fhRecoVertexZ_J->Fill(reco_pvertex_w2.Z());
    fhRecoVertexZ_fine_J->Fill(reco_pvertex_w2.Z());
  }

  if (flag == 2) {  // Saleve Arm
    fPrimaryVertexS.SetX(reco_pvertex_w2.X());
    fPrimaryVertexS.SetY(reco_pvertex_w2.Y());
    fPrimaryVertexS.SetZ(reco_pvertex_w2.Z());
    fhRecoVertexXY_S->Fill(reco_pvertex_w2.X(), reco_pvertex_w2.Y());
    fhRecoVertexZ_S->Fill(reco_pvertex_w2.Z());
    fhRecoVertexZ_fine_S->Fill(reco_pvertex_w2.Z());
  }
}

//_____________________________________________________________
void Na61PrimaryVertexRecoHTModule::CheckPrimaryVertexResolutionHT(int stage, UEventNode* out) {
  // algorithm:

  TObjArray tracksSub1;
  tracksSub1.Clear();
  TObjArray tracksSub2;
  tracksSub2.Clear();

  int jura_tracks_Sub1 = 0;
  int jura_tracks_Sub2 = 0;
  int ii = 0;

  UDataTable* tracktab = 0;
  if (fTrackInput == 0) tracktab = out->GetDataTable(Form("Vd HT Tracks"));
  if (fTrackInput == 1) tracktab = out->GetDataTable(Form("Tpc Tracks"));
  if (fTrackInput == 2) tracktab = out->GetDataTable(Form("VdCa Tracks"));

  if (!tracktab) return;

  for (int i = 0; i < tracktab->GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracktab->At(i);

    // if(!(track->GetFlag()==flag))continue;
    if (!(track->GetTagForVtx())) continue;

    if ((ii % 2) == 0) {
      tracksSub1.Add(track);
      jura_tracks_Sub1++;
    } else {
      tracksSub2.Add(track);
      jura_tracks_Sub2++;
    }

    ii++;
  }

  if (tracksSub1.GetEntries() < 2) return;
  if (tracksSub2.GetEntries() < 2) return;

  // cout<<"jura_sub1="<<jura_tracks_Sub1<<"  jura_sub2="<<jura_tracks_Sub2
  //    <<" saleve_sub1="<<saleve_tracks_Sub1<<"  saleve_sub2="<<saleve_tracks_Sub2<<endl;

  //////////////////////////////////// Sub1 primary vertex //////////////////////////////////////////////////

  double sum_up_x = 0;
  double sum_do_x = 0;
  double sum_up_y = 0;
  double sum_do_y = 0;

  for (int i = 0; i < tracksSub1.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracksSub1.At(i);
    for (int j = i + 1; j < tracksSub1.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracksSub1.At(j);

      double axi = track1->GetDX_f() / track1->GetDZ_f();
      double axj = track2->GetDX_f() / track2->GetDZ_f();
      double bxi = track1->GetX_f();
      double bxj = track2->GetX_f();

      double ayi = track1->GetDY_f() / track1->GetDZ_f();
      double ayj = track2->GetDY_f() / track2->GetDZ_f();
      double byi = track1->GetY_f();
      double byj = track2->GetY_f();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < fZCut) {
        sum_up_x = sum_up_x + (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - fZprim) < fZCut) {
        sum_up_y = sum_up_y + (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + (ayi - ayj) * (ayi - ayj);
      }
    }
  }

  double zx_min = 0;
  if (sum_do_x > 0) zx_min = -sum_up_x / sum_do_x;
  double zy_min = 0;
  if (sum_do_y > 0) zy_min = -sum_up_y / sum_do_y;

  // cout<<sum_do_x<<"  "<<sum_do_y<<"  "<<zx_min<<" "<<zy_min<<endl;

  // double z_prim = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);
  // double z_prim = -(sum_up_x)/(sum_do_x);
  // double z_prim = -(sum_up_y)/(sum_do_y); // use only  y information

  // double sigx=0.12; double sigy=0.07;
  double sigx = 0.1;
  double sigy = 0.1;
  double z_prim = zx_min / (sigx * sigx) + zy_min / (sigy * sigy);
  double norm = 1. / (sigx * sigx) + 1. / (sigy * sigy);
  z_prim = z_prim / norm;
  //  double zprim_sub1 = z_prim;

  Vector3D reco_pvertex_w2_Sub1;
  // PrimaryVertexWithWeigths_Y(flag,tracksSub1,z_prim,reco_pvertex_w2_Sub1);
  // PrimaryVertexWithWeigths_X(flag,tracksSub1,z_prim,reco_pvertex_w2_Sub1);
  PrimaryVertexWithWeigths_XY(2, tracksSub1, z_prim, reco_pvertex_w2_Sub1, sigx, sigy);
  // PrimaryVertexWithWeigths_XY2(flag,tracksSub1,z_prim,reco_pvertex_w2_Sub1,sigx,sigy);

  ////////////////////////////////////// Sub2 primary vertex //////////////////////////////////////////////////

  sum_up_x = 0;
  sum_do_x = 0;
  sum_up_y = 0;
  sum_do_y = 0;

  for (int i = 0; i < tracksSub2.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracksSub2.At(i);
    for (int j = i + 1; j < tracksSub2.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracksSub2.At(j);

      double axi = track1->GetDX_f() / track1->GetDZ_f();
      double axj = track2->GetDX_f() / track2->GetDZ_f();
      double bxi = track1->GetX_f();
      double bxj = track2->GetX_f();

      double ayi = track1->GetDY_f() / track1->GetDZ_f();
      double ayj = track2->GetDY_f() / track2->GetDZ_f();
      double byi = track1->GetY_f();
      double byj = track2->GetY_f();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < fZCut) {
        sum_up_x = sum_up_x + (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - fZprim) < fZCut) {
        sum_up_y = sum_up_y + (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + (ayi - ayj) * (ayi - ayj);
      }
    }
  }

  zx_min = 0;
  if (sum_do_x > 0) zx_min = -sum_up_x / sum_do_x;
  zy_min = 0;
  if (sum_do_y > 0) zy_min = -sum_up_y / sum_do_y;

  // z_prim = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);
  // z_prim = -(sum_up_x)/(sum_do_x);
  // z_prim = -(sum_up_y)/(sum_do_y); // use only  y information

  z_prim = zx_min / (sigx * sigx) + zy_min / (sigy * sigy);
  norm = 1. / (sigx * sigx) + 1. / (sigy * sigy);
  z_prim = z_prim / norm;
  //  double zprim_sub2 = z_prim;

  Vector3D reco_pvertex_w2_Sub2;
  // PrimaryVertexWithWeigths_Y(tracksSub2,z_prim,reco_pvertex_w2_Sub2);
  // PrimaryVertexWithWeigths_X(tracksSub2,z_prim,reco_pvertex_w2_Sub2);
  PrimaryVertexWithWeigths_XY(2, tracksSub2, z_prim, reco_pvertex_w2_Sub2, sigx, sigy);
  // PrimaryVertexWithWeigths_XY2(tracksSub2,z_prim,reco_pvertex_w2_Sub2,sigx,sigy);

  fhVtxDx_Glob2[stage]->Fill(reco_pvertex_w2_Sub2.X() - reco_pvertex_w2_Sub1.X());
  fhVtxDy_Glob2[stage]->Fill(reco_pvertex_w2_Sub2.Y() - reco_pvertex_w2_Sub1.Y());
  fhVtxDz_Glob2[stage]->Fill(reco_pvertex_w2_Sub2.Z() - reco_pvertex_w2_Sub1.Z());

  fhVtxDx_fine_Glob2[stage]->Fill(reco_pvertex_w2_Sub2.X() - reco_pvertex_w2_Sub1.X());
  fhVtxDy_fine_Glob2[stage]->Fill(reco_pvertex_w2_Sub2.Y() - reco_pvertex_w2_Sub1.Y());
  fhVtxDz_fine_Glob2[stage]->Fill(reco_pvertex_w2_Sub2.Z() - reco_pvertex_w2_Sub1.Z());

  if (stage == 1) {
    z_prim = (reco_pvertex_w2_Sub1.Z() + reco_pvertex_w2_Sub2.Z()) / 2.;
    int tmulti = tracksSub1.GetEntries() + tracksSub2.GetEntries();
    double dx = (reco_pvertex_w2_Sub2.X() - reco_pvertex_w2_Sub1.X()) / 2.;
    double dy = (reco_pvertex_w2_Sub2.Y() - reco_pvertex_w2_Sub1.Y()) / 2.;
    double dz = (reco_pvertex_w2_Sub2.Z() - reco_pvertex_w2_Sub1.Z()) / 2.;
    fhVtxDx_fine_vs_multi->Fill(tmulti, dx);
    fhVtxDy_fine_vs_multi->Fill(tmulti, dy);
    fhVtxDz_fine_vs_multi->Fill(tmulti, dz);
    /*
    int la=-1;
    if((z_prim>-49.0) && (z_prim<-47.8))la=0;
    if((z_prim>-47.7) && (z_prim<-46.6))la=1;
    if((z_prim>-46.5) && (z_prim<-45.4))la=2;
    if(la!=-1){
      fhVtxDx_fine_vs_multi_La[la] -> Fill(tmulti,dx);
      fhVtxDy_fine_vs_multi_La[la] -> Fill(tmulti,dy);
      fhVtxDz_fine_vs_multi_La[la] -> Fill(tmulti,dz);
    }
    */
  }
}

//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoHTModule::PrimaryVertexWithWeigths_XY(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex, double sigx, double sigy) {
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

      double x1 = track1->GetXatZ_f(zprim);
      double x2 = track2->GetXatZ_f(zprim);
      double y1 = track1->GetYatZ_f(zprim);
      double y2 = track2->GetYatZ_f(zprim);

      double d2x = (x2 - x1) * (x2 - x1);
      double d2y = (y2 - y1) * (y2 - y1);

      double wx = 1;
      double wy = 1;
      if (d2x > 0) wx = 1. / d2x;
      if (d2y > 0) wy = 1. / d2y;

      // fhd2->Fill(d2);

      /////////////////// stuff for new method
      double axi = track1->GetDX_f() / track1->GetDZ_f();
      double axj = track2->GetDX_f() / track2->GetDZ_f();
      double bxi = track1->GetX_f();
      double bxj = track2->GetX_f();

      double ayi = track1->GetDY_f() / track1->GetDZ_f();
      double ayj = track2->GetDY_f() / track2->GetDZ_f();
      double byi = track1->GetY_f();
      double byj = track2->GetY_f();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < fZCut || flag != 2) {
        sum_up_x = sum_up_x + wx * (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + wx * (axi - axj) * (axi - axj);
      }
      if (TMath::Abs(z_y - fZprim) < fZCut || flag != 2) {
        sum_up_y = sum_up_y + wy * (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + wy * (ayi - ayj) * (ayi - ayj);
      }
      // cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
    }
  }

  double zprim_w2x = 0;
  if (sum_do_x > 0) zprim_w2x = -sum_up_x / sum_do_x;
  double zprim_w2y = 0;
  if (sum_do_y > 0) zprim_w2y = -sum_up_y / sum_do_y;

  // cout<<" _XY:"<< sum_up_x<<" "<< sum_do_x<<"    "<<sum_up_y<<" "<< sum_do_y<<"   "<<zprim_w2x<<" "<<zprim_w2y<<endl;

  double zprim_w2 = zprim_w2x / (sigx * sigx) + zprim_w2y / (sigy * sigy);
  double norm = 1. / (sigx * sigx) + 1. / (sigy * sigy);
  zprim_w2 = zprim_w2 / norm;

  // cout<<" _XY: zprim_w2 = "<<zprim_w2<<" norm="<<norm<<endl;

  //
  double sum_x = 0;
  double sum_y = 0;
  double Nx = 0;
  double Ny = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);
    
    if((track->GetTagForVtx()!=2)  && fTrackInput==1 && fStage) continue;


    double x = track->GetXatZ_f(zprim_w2);
    double y = track->GetYatZ_f(zprim_w2);

    sum_y = sum_y + y;
    Ny = Ny + 1;      

    if((track->GetCombMeth()!=10) && fTrackInput==1)continue; 

    sum_x = sum_x + x;      
    Nx = Nx + 1; 
  }

  double xprim = sum_x/Nx;
  double yprim = sum_y/Ny;

  sum_x = 0;
  sum_y = 0;

  double norm_x = 0;
  double norm_y = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);
    
        if((track->GetTagForVtx()!=2)  && fTrackInput==1 && fStage) continue;
      // discard track for x if does not fulfill momDiffCut

    double x = track->GetXatZ_f(zprim_w2);
    double y = track->GetYatZ_f(zprim_w2);

    double d2x = (x - xprim) * (x - xprim);
    double d2y = (y - yprim) * (y - yprim);

    sum_y = sum_y + y/d2y;
    norm_y = norm_y + 1./d2y;      

    if((track->GetCombMeth()!=10) && fTrackInput==1)continue; 
    
    sum_x = sum_x + x/d2x;
    norm_x = norm_x + 1./d2x;      
  }
  
  double xprim_w2 = sum_x/norm_x;
  double yprim_w2 = sum_y/norm_y;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);
}

//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoHTModule::PrimaryVertexWithWeigths_XY2(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex, double /*sigx*/, double /*sigy*/) {
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

      double x1 = track1->GetXatZ_f(zprim);
      double x2 = track2->GetXatZ_f(zprim);
      double y1 = track1->GetYatZ_f(zprim);
      double y2 = track2->GetYatZ_f(zprim);

      double d2x = (x2 - x1) * (x2 - x1);
      double d2y = (y2 - y1) * (y2 - y1);
      double wx = 1. / d2x;
      double wy = 1. / d2y;

      // fhd2->Fill(d2);

      /////////////////// stuff for new method
      double axi = track1->GetDX_f() / track1->GetDZ_f();
      double axj = track2->GetDX_f() / track2->GetDZ_f();
      double bxi = track1->GetX_f();
      double bxj = track2->GetX_f();

      double ayi = track1->GetDY_f() / track1->GetDZ_f();
      double ayj = track2->GetDY_f() / track2->GetDZ_f();
      double byi = track1->GetY_f();
      double byj = track2->GetY_f();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));
      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_x - fZprim) < fZCut || flag != 2) {
        sum_up_x = sum_up_x + wx * (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + wx * (axi - axj) * (axi - axj);
      }

      if (TMath::Abs(z_y - fZprim) < fZCut || flag != 2) {
        sum_up_y = sum_up_y + wy * (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + wy * (ayi - ayj) * (ayi - ayj);
        // cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
      }
    }
  }

  // double zprim_w2x = -(sum_up_x)/(sum_do_x);
  // double zprim_w2y = -(sum_up_y)/(sum_do_y);

  double zprim_w2 = -(sum_up_x + sum_up_y) / (sum_do_x + sum_do_y);
  // double norm = 1./(sigx*sigx)  +  1./(sigy*sigy);
  // zprim_w2 =  zprim_w2 / norm;

  // cout<<" zprim_w2 = "<<zprim_w2<<endl;

  //
  double sum_x = 0;
  double sum_y = 0;
  double N = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);

    double x = track->GetXatZ_f(zprim_w2);
    double y = track->GetYatZ_f(zprim_w2);

    sum_x = sum_x + x;
    sum_y = sum_y + y;
    N = N + 1;
  }

  double xprim = sum_x / N;
  double yprim = sum_y / N;

  sum_x = 0;
  sum_y = 0;

  double norm_x = 0;
  double norm_y = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);

    double x = track->GetXatZ_f(zprim_w2);
    double y = track->GetYatZ_f(zprim_w2);

    double d2x = (x - xprim) * (x - xprim);
    double d2y = (y - yprim) * (y - yprim);

    sum_x = sum_x + x / d2x;
    sum_y = sum_y + y / d2y;
    norm_x = norm_x + 1. / d2x;
    norm_y = norm_y + 1. / d2y;
  }

  double xprim_w2 = sum_x / norm_x;
  double yprim_w2 = sum_y / norm_y;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);
}

//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoHTModule::PrimaryVertexWithWeigths_X(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex) {
  // algorithm:
  // 1. Calculate weights = 1/dist^2, distance betweem i-th and j-th tracks on
  // the closest proximity plane (z=zprim).
  // 2. Used the weights to calculate z coordinate of primary vertex from the formula that inclused weights (see PS notes).
  // 3. Find x,y coordinates averageing tracks x, y at the closes proximmity plane.

  double sum_up_x = 0;
  double sum_do_x = 0;
  // double sum_up_y = 0;
  // double sum_do_y = 0;

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      double x1 = track1->GetXatZ_f(zprim);
      double x2 = track2->GetXatZ_f(zprim);
      //      double y1 = track1->GetYatZ_f(zprim);
      //      double y2 = track2->GetYatZ_f(zprim);

      double d2 = (x2 - x1) * (x2 - x1);
      // double d2 = (y2-y1)*(y2-y1);
      double w = 1. / d2;

      // fhd2->Fill(d2);

      /////////////////// stuff for new method
      double axi = track1->GetDX_f() / track1->GetDZ_f();
      double axj = track2->GetDX_f() / track2->GetDZ_f();
      double bxi = track1->GetX_f();
      double bxj = track2->GetX_f();

      // double ayi = track1->GetDY()/track1->GetDZ();
      // double ayj = track2->GetDY()/track2->GetDZ();
      // double byi = track1->GetY();
      // double byj = track2->GetY();

      double z_x = ((axi - axj) * (bxi - bxj)) / ((axi - axj) * (axi - axj));

      if (TMath::Abs(z_x - fZprim) < fZCut || flag != 2) {
        sum_up_x = sum_up_x + w * (axi - axj) * (bxi - bxj);
        sum_do_x = sum_do_x + w * (axi - axj) * (axi - axj);
      }
      // sum_up_y = sum_up_y + w*(ayi-ayj)*(byi-byj);
      // sum_do_y = sum_do_y + w*(ayi-ayj)*(ayi-ayj);
      // cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
    }
  }

  double zprim_w2 = -(sum_up_x) / (sum_do_x);
  // double zprim_w2 = -(sum_up_y)/(sum_do_y);

  //
  double sum_x = 0;
  double sum_y = 0;
  double N = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);

    double x = track->GetXatZ_f(zprim_w2);
    double y = track->GetYatZ_f(zprim_w2);

    sum_x = sum_x + x;
    sum_y = sum_y + y;
    N = N + 1;
  }

  double xprim = sum_x / N;
  double yprim = sum_y / N;

  sum_x = 0;
  sum_y = 0;

  double norm_x = 0;
  double norm_y = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);

    double x = track->GetXatZ_f(zprim_w2);
    double y = track->GetYatZ_f(zprim_w2);
    // cout<<"x="<<x<<endl;

    double d2x = (x - xprim) * (x - xprim);
    double d2y = (y - yprim) * (y - yprim);

    sum_x = sum_x + x / d2x;
    sum_y = sum_y + y / d2y;
    norm_x = norm_x + 1. / d2x;
    norm_y = norm_y + 1. / d2y;
  }

  double xprim_w2 = sum_x / norm_x;
  double yprim_w2 = sum_y / norm_y;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);
}

//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoHTModule::PrimaryVertexWithWeigths_Y(int flag, TObjArray& tracks, double zprim, Vector3D& pvertex) {
  // algorithm:
  // 1. Calculate weights = 1/dist^2, distance betweem i-th and j-th tracks on
  // the closest proximity plane (z=zprim).
  // 2. Used the weights to calculate z coordinate of primary vertex from the formula that inclused weights (see PS notes).
  // 3. Find x,y coordinates averageing tracks x, y at the closes proximmity plane.

  // double sum_up_x = 0;
  // double sum_do_x = 0;
  double sum_up_y = 0;
  double sum_do_y = 0;

  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)tracks.At(i);
    for (int j = i + 1; j < tracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)tracks.At(j);

      //      double x1 = track1->GetXatZ_f(zprim);
      //      double x2 = track2->GetXatZ_f(zprim);
      double y1 = track1->GetYatZ_f(zprim);
      double y2 = track2->GetYatZ_f(zprim);

      // double d2x = (x2-x1)*(x2-x1);
      double d2 = (y2 - y1) * (y2 - y1);
      double w = 1. / d2;

      // fhd2->Fill(d2);

      /////////////////// stuff for new method
      // double axi = track1->GetDX()/track1->GetDZ();
      // double axj = track2->GetDX()/track2->GetDZ();
      // double bxi = track1->GetX();
      // double bxj = track2->GetX();

      double ayi = track1->GetDY_f() / track1->GetDZ_f();
      double ayj = track2->GetDY_f() / track2->GetDZ_f();
      double byi = track1->GetY_f();
      double byj = track2->GetY_f();

      // sum_up_x = sum_up_x + wx*(axi-axj)*(bxi-bxj);
      // sum_do_x = sum_do_x + wx*(axi-axj)*(axi-axj);

      double z_y = ((ayi - ayj) * (byi - byj)) / ((ayi - ayj) * (ayi - ayj));

      if (TMath::Abs(z_y - fZprim) < fZCut || flag != 2) {
        sum_up_y = sum_up_y + w * (ayi - ayj) * (byi - byj);
        sum_do_y = sum_do_y + w * (ayi - ayj) * (ayi - ayj);
      }
      // cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
    }
  }

  // double zprim_w2_x = -(sum_up_x)/(sum_do_x);
  double zprim_w2 = -(sum_up_y) / (sum_do_y);

  //
  double sum_x = 0;
  double sum_y = 0;
  double N = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);

    double x = track->GetXatZ_f(zprim_w2);
    double y = track->GetYatZ_f(zprim_w2);

    sum_x = sum_x + x;
    sum_y = sum_y + y;
    N = N + 1;
  }

  double xprim = sum_x / N;
  double yprim = sum_y / N;

  sum_x = 0;
  sum_y = 0;

  double norm_x = 0;
  double norm_y = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);

    double x = track->GetXatZ_f(zprim_w2);
    double y = track->GetYatZ_f(zprim_w2);

    double d2x = (x - xprim) * (x - xprim);
    double d2y = (y - yprim) * (y - yprim);

    sum_x = sum_x + x / d2x;
    sum_y = sum_y + y / d2y;
    norm_x = norm_x + 1. / d2x;
    norm_y = norm_y + 1. / d2y;
  }

  double xprim_w2 = sum_x / norm_x;
  double yprim_w2 = sum_y / norm_y;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);
}

//_____________________________________________________________
void Na61PrimaryVertexRecoHTModule::End() {
  // Run-level finalisation
  SetState(kEnd);
  Na61PrimaryVertexRecoModule::End();
}

//____________________________________________________________________
void Na61PrimaryVertexRecoHTModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);
  Na61PrimaryVertexRecoModule::Finish();
}

//____________________________________________________________________
void Na61PrimaryVertexRecoHTModule::Print(Option_t* option) const {
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>
  Na61PrimaryVertexRecoModule::Print(option);
}
