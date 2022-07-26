//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61VdTrackingHTModule
#include "Na61VdTrackingHTModule.h"
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
#ifndef ROOT_TF1
#include "TF1.h"
#endif
#ifndef ROOT_TGraphErrors
#include "TGraphErrors.h"
#endif
#ifndef ROOT_TMath
#include "TMath.h"
#endif
#include "iostream"
#include "UG4RecoTrack.h"
#include "UMmCluster.h"
#include "TAxis.h"
#include "TH2F.h"

//____________________________________________________________________
// ClassImp(Na61VdTrackingHTModule);

//____________________________________________________________________
Na61VdTrackingHTModule::Na61VdTrackingHTModule() : Na61VdTrackingModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);

  fOldSetup = 0;
  fAddInGrid2Tracks = 0;
  fFiduCut = 3.0;
  fDetectorSetupOption = 0;
  fMaxXInd = 0;
  fMaxYInd = 0;
  fSignalTracks = 0;
  fMakeClusters = 0;
  fAdd4HitTracks = false;
}
//____________________________________________________________________
Na61VdTrackingHTModule::Na61VdTrackingHTModule(const char* name, const char* title) : Na61VdTrackingModule(name, title) {
  // Named Constructor
  SetState(kSetup);

  fOldSetup = 0;
  fAddInGrid2Tracks = 0;
  fFiduCut = 3.0;
  fDetectorSetupOption = 0;
  fMaxXInd = 0;
  fMaxYInd = 0;

  fHitInfoArray = new TObjArray();
  fHitInfoArray->SetOwner();
  fProduction = 0;  //
  fSignalTracks = 0;
  fMakeClusters = 0;
  fAdd4HitTracks = false;

  fChi2Cut = 5;
  fVzOffset = 0;
}

//____________________________________________________________________
void Na61VdTrackingHTModule::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  Na61VdTrackingModule::DefineHistograms();

  // TDirectory* histDir = gDirectory->mkdir(fHistDirName.Data());
  TDirectory* histDir = gDirectory->GetDirectory(fHistDirName.Data());

  histDir->cd();

  // hits multiplicity in detectors

  fhRecoTracks = new TH1F("hRecoTracks", "number of reco tracks", 1000, 0, 1000);
  fhRecoHtTracks = new TH1F("hRecoHtTracks", "number of reco tracks", 1000, 0, 1000);
  fhCurvatureComb_back = new TH1F("hCurvatureComb_back", "", 500, -0.00003, 0.00003);

  int nbins = 5000;
  double xmin = -2.0;
  double xmax = 2.0;

  TString strArr[4] = {"x1", "x2", "x3", "x4"};

  for (int i = 0; i < 4; i++) {
    fhDevX[i] = new TH1F(Form("hDevX_%s", strArr[i].Data()), "", nbins, xmin, xmax);
    fhDevY[i] = new TH1F(Form("hDevY_%s", strArr[i].Data()), "", nbins, xmin, xmax);
    fhDevX_acce[i] = new TH1F(Form("hDevX_%s_acce", strArr[i].Data()), "", nbins, xmin, xmax);
    fhDevY_acce[i] = new TH1F(Form("hDevY_%s_acce", strArr[i].Data()), "", nbins, xmin, xmax);
    xmax = 1.0;
    xmin = -1.0;
    fhResX[i] = new TH1F(Form("hResX_%s", strArr[i].Data()), "", nbins, xmin, xmax);
    fhResY[i] = new TH1F(Form("hResY_%s", strArr[i].Data()), "", nbins, xmin, xmax);
    fhResX_acce[i] = new TH1F(Form("hResX_%s_acce", strArr[i].Data()), "", nbins, xmin, xmax);
    fhResY_acce[i] = new TH1F(Form("hResY_%s_acce", strArr[i].Data()), "", nbins, xmin, xmax);
  }
  // here we define binning for the Hough space

  xmax = 0.6;
  xmin = -0.6;
  double ymax = 0.3;
  double ymin = -0.3;
  cout << "fbwx  = " << fbwx << "   fbwy  = " << fbwy << endl;
  findx = int((xmax - xmin) / fbwx);
  cout << "findx  = " << findx << endl;
  findy = int((ymax - ymin) / fbwy);
  cout << "findy = " << findy << endl;
  double pixelSize = 0.04 / 50.;
  cout << "Ratio of HSaxay binx width and cluster structure:" << fbwx / pixelSize << endl;
  cout << "Ratio of HSaxay biny width and cluster structure:" << fbwy / pixelSize << endl;

  fhHSaxay = new TH2D("hHSaxay", "", findx, xmin, xmax, findy, ymin, ymax);
  fhHSaxay_clust = new TH2D("hHSaxay_clust", "", findx, xmin, xmax, findy, ymin, ymax);

  // fhaxayreco = new TH2F("axayreco","",1000,-0.6,0.6,500,-0.3,0.3);
  // fhaxayrecoselect = new TH2F("axayrecoselect","",1000,-0.6,0.6,500,-0.3,0.3);
  // fhaxayefficiency = new TH2F("axayefficiency","",1000,-0.6,0.6,500,-0.3,0.3);
  // fhaxayHt = new TH2F("axayHt","",1000,-0.6,0.6,500,-0.3,0.3);
  // fhaxayHtselect = new TH2F("axayHtselect","",1000,-0.6,0.6,500,-0.3,0.3);

  fhChi2xNdf_N3 = new TH1F("hChi2xNdf_N3", "test", 1000, 0, 50);
  fhChi2xNdf_N4 = new TH1F("hChi2xNdf_N4", "test", 1000, 0, 50);
  fhChi2yNdf_N3 = new TH1F("hChi2yNdf_N3", "test", 1000, 0, 50);
  fhChi2yNdf_N4 = new TH1F("hChi2yNdf_N4", "test", 1000, 0, 50);

  fhChi2Ndf = new TH1F("hChi2Ndf", "", 1000, 0, 50);
  fhChi2xNdf = new TH1F("hChi2xNdf", "", 1000, 0, 50);
  fhChi2yNdf = new TH1F("hChi2yNdf", "", 1000, 0, 50);

  fhChi2Ndf_Acce = new TH1F("hChi2Ndf_Acce", "", 1000, 0, 50);
  fhChi2xNdf_Acce = new TH1F("hChi2xNdf_Acce", "", 1000, 0, 50);
  fhChi2yNdf_Acce = new TH1F("hChi2yNdf_Acce", "", 1000, 0, 50);

  fGr_axay_vdtrack = new TGraph();
  fGr_axay_vdtrack->SetName("Gr_axay_vdtrack");
  fGr_axay_vdtrack->SetTitle("");

  fhAdded4hitTracks = new TH1F("hAdded4hitTracks", "", 50, 0, 50);

  fhPrimaryVertexZ = new TH1F(Form("hPrimaryVertexZ"), "", 1000, -70., -30.);
  fhPrimaryVertexZ_acce = new TH1F(Form("hPrimaryVertexZ_acce"), "", 1000, -70., -30.);

  gDirectory->Append(fGr_axay_vdtrack);
  gDirectory->cd("..");
}

//____________________________________________________________________
void Na61VdTrackingHTModule::Init() {
  // Job-level initialisation
  SetState(kInit);
  Na61VdTrackingModule::Init();

  fVdParams = (Na61VdParametersManager::Instance())->GetVdParams();

  // dont understant why so big cut is needed

  cout << "fChi2Cut = " << fChi2Cut << endl;
}

//____________________________________________________________________
void Na61VdTrackingHTModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);

  // cout<<"track 3 "<<endl;
}

//____________________________________________________________________
void Na61VdTrackingHTModule::Event(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode) {
  // Per event method

  SetState(kEvent);

   //cout<<"  Na61VdTrackingHTModule ------------------------------------------------>new event "<<endl;

  RecoWithSimultaneousHT(inNodeJ, inNodeS, outNode);

  // Add 4 hit tracks from combinatorical method if not found in HT space
  if (fAdd4HitTracks) Add4HitTracks(inNodeJ, inNodeS, outNode);

  RefitTracks(outNode);
  //cout<<"finished Na61VdTrackingHTModule"<<endl;
}

//_________________________________________________________________
void Na61VdTrackingHTModule::Add4HitTracks(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode) {
  // tracks is array containing all 4 hit tracks creaded in the combinatorial method
  // UDataTable* table = outNode->GetDataTable("Vd HT Tracks");
  UDataTable* table = outNode->GetDataTable(fOutTableName.Data());

  if (table)
    ((UVdEvent*)outNode)->SetPrimaryTracks(table->GetEntries());  // number of tracks before adding 4hit tracks
  else
    ((UVdEvent*)outNode)->SetPrimaryTracks(0);

  int N = 0;
  if (table) N = table->GetEntries();

  fAllTracks = 0;

  // int tabN = 4;
  int tabN = 18;

  TObjArray tracks;
  tracks.Clear();

  int jura_tracks = 0;
  if (inNodeJ) {
    for (int it = 0; it < tabN; it++) {
      UDataTable* tracktab = inNodeJ->GetDataTable(Form("FullTracks %s", fMatchStr[it].Data()));
      if (!tracktab) continue;

      for (int i = 0; i < tracktab->GetEntries(); i++) {
        UVdTrack* track = (UVdTrack*)tracktab->At(i);
        if (track->GetFlag() == 2) {
          tracks.Add(track);
          jura_tracks++;
          SetTrackHits(track, inNodeJ);
        }
      }
    }
  }

  int saleve_tracks = 0;
  if (inNodeS) {
    for (int it = 0; it < tabN; it++) {
      UDataTable* tracktab = inNodeS->GetDataTable(Form("FullTracks %s", fMatchStr[it].Data()));
      if (!tracktab) continue;

      for (int i = 0; i < tracktab->GetEntries(); i++) {
        UVdTrack* track = (UVdTrack*)tracktab->At(i);
        if (track->GetFlag() == 2) {
          tracks.Add(track);
          saleve_tracks++;
          SetTrackHits(track, inNodeS);
        }
      }
    }
  }

  // cout<<"HT tracks = "<<N<<" jura 4hit tracks = "<<jura_tracks<<"  saleve 4hit tracks = "<<saleve_tracks<<" table: "<<table<<endl;

  if (!table && ((jura_tracks + saleve_tracks) > 0)) {
    // table = new UDataTable("Vd HT Tracks");
    table = new UDataTable(fOutTableName.Data());
    table->SetOwner();
  }

  int added_4hit_tracks = 0;
  for (int i = 0; i < tracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)tracks.At(i);
    // cout<<"checking match of track: "<<track<<"  i="<<i<<endl;

    bool match = false;
    for (int j = 0; j < N; j++) {
      UVdTrack* trackht = (UVdTrack*)table->At(j);
      bool found_match = CheckTrackMatching2(track, trackht);
      // cout<<"found_match="<<found_match<<endl;
      if (found_match) {
        match = true;
        break;
      }
    }

    if (!match) {  // clone track and add to HT table
      added_4hit_tracks++;
      // cout<<"Adding Track: "<<track<<endl;
      // int ii;
      // cin>>ii;
      track->Activate();
      track->MarkForRemoval(false);
      Line3D* line = track->Getline();
      UVdTrack* ctrack = new UVdTrack(line->GetOrigin(), line->GetDirection());
      ctrack->SetTrackID(2);  // Added 4hit track

      // set the track here (set hits to track)
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

      ctrack->SetVdsHitIDs(hit1, hit2, hit3, hit4);
      ctrack->SetHitArrayIndexes(ind1, ind2, ind3, ind4);
      ctrack->SetTabArrayIndexes(tind1, tind2, tind3, tind4);
      ctrack->SetCombMeth(1);

      table->Add(ctrack);
    }
  }

  fhAdded4hitTracks->Fill(added_4hit_tracks);

  if (table) {
    ((UVdEvent*)outNode)->SetAdded4HitTracks(added_4hit_tracks);  // number of tracks before adding 4hit tracks
    fAllTracks = table->GetEntries();
  }
  fhRecoTracks->Fill(fAllTracks);

  // cout<<"Added 4 hit tracks : "<<added_4hit_tracks<<endl;
  // int primaryTracks = ((UVdEvent*)outNode)->GetPrimaryTracks();
  // cout<<"All tracks = "<<fAllTracks<<"    primary = "<<primaryTracks<<" primary+added4hit = "<<added_4hit_tracks+primaryTracks<<endl;
}

//_____________________________________________________________
void Na61VdTrackingHTModule::SetTrackHits(UVdTrack* track, UEventNode* in) {
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

  if (ind1 != -1) hit1 = (USensorHit*)(in->GetDataTable(Form("Hits %s", fSensorNames[tind1].Data())))->At(ind1);
  if (ind2 != -1) hit2 = (USensorHit*)(in->GetDataTable(Form("Hits %s", fSensorNames[tind2].Data())))->At(ind2);
  if (ind3 != -1) hit3 = (USensorHit*)(in->GetDataTable(Form("Hits %s", fSensorNames[tind3].Data())))->At(ind3);
  if (ind4 != -1) hit4 = (USensorHit*)(in->GetDataTable(Form("Hits %s", fSensorNames[tind4].Data())))->At(ind4);

  track->SetVdsHitIDs(hit1, hit2, hit3, hit4);
}

//_____________________________________________________________
bool Na61VdTrackingHTModule::CheckTrackMatching2(UVdTrack* track1, UVdTrack* track2) {
  // TWO or more same hits in the tracks returns true

  int isame = 0;
  for (int i = 0; i < 4; i++) {
    if ((track1->GetHitIdAtStation(i) == track2->GetHitIdAtStation(i)) && track1->GetHitIdAtStation(i)) isame++;
    // cout<<"CheckTrackMatching: "<<i<<"  ID1="<<track1->GetHitIdAtStation(i)<<" ID2="<<track2->GetHitIdAtStation(i)<<endl;
  }

  if (isame < 2) return false;

  return true;
}

//__________________________________________________________________________________________
void Na61VdTrackingHTModule::FillCombinatorialCurvature(UEventNode* in) {
  // methode used to fill fhCurvature_back for combinatorial tracks to compare with
  // HT tracks

  for (int i = 0; i < 18; i++) {
    UDataTable* tracktab = in->GetDataTable(Form("FullTracks %s", fMatchStr[i].Data()));
    if (!tracktab) continue;

    for (int j = 0; j < tracktab->GetEntries(); j++) {
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
      // track -> Activate();

      fhCurvatureComb_back->Fill(track->GetCurvature());
    }
  }
}

//_________________________________________________________________
void Na61VdTrackingHTModule::RecoWithSimultaneousHT(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode) {
  // algorithm:
  // tracking algorithm is based on the Hough transform
  // For each station we do simultaneous transform in x and y direction assuming
  // narrow corridor for track origin. Thus for each point we obtain
  // ax and ay slope parameters.
  // This procedure might be eventualy followed by fitting that will
  // provide better parameters.

  fAllHtTracks = 0;

  double EventVertexX;
  double EventVertexY;
  double EventVertexZ;
  double EventVertexZorg;
  // std::cout<<"Primary vertex z = "<<((UVdEvent*)outNode)->GetPrimaryVertexZ()<<endl;
  // std::cout<<"Primary vertex x = "<<((UVdEvent*)outNode)->GetPrimaryVertexX()<<endl;
  // std::cout<<"Primary vertex y = "<<((UVdEvent*)outNode)->GetPrimaryVertexY()<<endl;
  // std::cout<<"Primary vertex status = "<<((UVdEvent*)outNode)->GetPrimaryVertexStatus()<<endl;

  if (((UVdEvent*)outNode)->GetPrimaryVertexStatus()) {
    EventVertexX = ((UVdEvent*)outNode)->GetPrimaryVertexX();
    EventVertexY = ((UVdEvent*)outNode)->GetPrimaryVertexY();
    EventVertexZ = ((UVdEvent*)outNode)->GetPrimaryVertexZ() + fVzOffset * 0.001;
    EventVertexZorg = ((UVdEvent*)outNode)->GetPrimaryVertexZ();
    // std::cout<<" EventVertexZ = ((UVdEvent*)outNode)->GetPrimaryVertexZ() + fVzOffset*0.001 = "
    //<<((UVdEvent*)outNode)->GetPrimaryVertexZ() + fVzOffset*0.001<<std::endl;
  } else {
    Info(2, "RecoWithSimultaneousHT", Form("No primary vertex for this event. Skipping event"));
    //        std::cout<<"No primary vertex in HT module for this event!"<<std::endl;
    return;
  }

  double offx = fVdParams->GetBeamSpotOffsetX();
  double sigx = fVdParams->GetBeamSpotSigmaX();
  double offy = fVdParams->GetBeamSpotOffsetY();
  double sigy = fVdParams->GetBeamSpotSigmaY();

  fhPrimaryVertexZ->Fill(EventVertexZorg);

  double dd = ((EventVertexX - offx) * (EventVertexX - offx)) / (sigx * sigx) + ((EventVertexY - offy) * (EventVertexY - offy)) / (sigy * sigy);
  // cout<<"dd="<<dd<<" "<<offx<<" "<< sigx<<" "<<offy<<" "<<sigy<<endl;
  if (dd > 5 * 5) {
    // cout<<"RecoWithSimultaneousHT: event not accepted:"<<endl;
    // cout<<"Primary vertex x = "<<EventVertexX<<" "<<offx<<endl;
    // cout<<"Primary vertex y = "<<EventVertexY<<endl;
    // cout<<"Primary vertex z = "<<EventVertexZ<<endl;
    return;
  }

  FillCombinatorialCurvature(inNodeJ);
  FillCombinatorialCurvature(inNodeS);

  fhPrimaryVertexZ_acce->Fill(EventVertexZorg);

  /*
  if( !((EventVertexX>-3.0) && (EventVertexX<5.0)) ){
    cout<<"RecoWithSimultaneousHT: event not accepted (bad X):"<<endl;
    cout<<"Primary vertex x = "<<EventVertexX<<endl;
    cout<<"Primary vertex y = "<<EventVertexY<<endl;
    cout<<"Primary vertex z = "<<EventVertexZ<<endl;
    return;
  }

  if( !((EventVertexY>-4.0) && (EventVertexY<3.0)) ){
    cout<<"RecoWithSimultaneousHT: event not accepted (bad Y):"<<endl;
    cout<<"Primary vertex x = "<<EventVertexX<<endl;
    cout<<"Primary vertex y = "<<EventVertexY<<endl;
    cout<<"Primary vertex z = "<<EventVertexZ<<endl;
    return;
  }
  */

  // cout<<"RecoWithSimultaneousHT: PrimVertexZ="<< EventVertexZ <<endl;

  // UDataTable* table = new UDataTable("Vd HT Tracks");
  UDataTable* table = new UDataTable(fOutTableName.Data());
  table->SetOwner();

  TObjArray* hitsVds1 = GetHits("Vds1", inNodeJ, inNodeS);
  TObjArray* hitsVds2 = GetHits("Vds2", inNodeJ, inNodeS);
  TObjArray* hitsVds3 = GetHits("Vds3", inNodeJ, inNodeS);
  TObjArray* hitsVds4 = GetHits("Vds4", inNodeJ, inNodeS);

  // if(!fAdd4HitTracks)cout<<"  hitsVds1: "<<hitsVds1<<endl;

  fHitInfoArray->Clear();

  //  int N1=0, N2=0, N3=0, N4=0;
  //  if(hitsVds1) N1 = hitsVds1->GetEntries();
  //  if(hitsVds2) N2 = hitsVds2->GetEntries();
  //  if(hitsVds3) N3 = hitsVds3->GetEntries();
  //  if(hitsVds4) N4 = hitsVds4->GetEntries();
  //  cout<<"N1 = "<<N1<<" "<<N2<<" "<<N3<<" "<<N4<<endl;

  TObjArray* vdsArr[] = {hitsVds1, hitsVds2, hitsVds3, hitsVds4};

  double bx = 0;  // assume that all tracks originate from x=y=z=0
  double by = 0;
  // cout<<"Tracking Module "<<endl;
  // int axayToInd[1000][1000];  //need to modify it

  // for(int i=0;i<1000;i++){
  // for(int j=0;j<1000;j++){
  //  faxayToInd[i][j]=-1;
  //}
  //}

  memset(faxayToInd, -1, sizeof(faxayToInd));

  // cout<<"track 9 "<<endl;
  // cout<<"working "<<endl;
  //////// find ax and ay in Vds1
  /////// using Hough transform

  // cout<<"---------------------------------> start HT points="<< fhHSaxay->GetEntries()<<" integral="<<fhHSaxay->Integral()<<endl;

  int hitsInfo1 = 0;
  int ind = 0;
  for (int itab = 0; itab < 4; itab++) {
    TObjArray* vdsHits = vdsArr[itab];
    // cout<<" itab="<<itab<<"  vdsHits: "<<vdsHits<<"  entries: "<<vdsHits->GetEntries()<<endl;
    for (int i = 0; i < vdsHits->GetEntries(); i++) {
      USensorHit* hit = (USensorHit*)vdsHits->At(i);

      // if(!fAdd4HitTracks)cout<<"used="<<hit->GetUsed()<<endl;
      // if(hit->GetUsed())continue; // don't use used hits

      double x = hit->GetX() - EventVertexX;
      double y = hit->GetY() - EventVertexY;
      double z = hit->GetZ() - EventVertexZ;

      double ax = (x - bx) / z;
      double ay = (y - by) / z;

      //      int OrigZ = (int)hit->GetZ();

      // cout<<"Track 1  "<<endl;

      if (!fProduction) fhHSaxay->Fill(ax, ay);
      int indx = fhHSaxay->GetXaxis()->FindBin(ax);
      int indy = fhHSaxay->GetYaxis()->FindBin(ay);
      // cout<<"binwidth matrix = "<< fhHSaxay->GetBinWidth(0) <<endl;
      // if((0.122<ax && ax<0.128) && (0.183<ay && ay<0.189)){
      // cout<<"indx="<<indx<<"  "<<"indy="<<indy<<" ax="<<ax<<" ay="<<ay<<" z="<<z<<endl;
      //}
      if (indx > 999) continue;
      if (indy > 999) continue;

      UHitInfo* hitInfo = 0;

      int ii = faxayToInd[indx][indy];
      // cout<<" ii = "<< ii <<" indx="<<indx<<"  indy="<<indy<<endl;
      if (ii != -1) {
        hitInfo = (UHitInfo*)fHitInfoArray->At(ii);
        // cout<<"hitInfo="<<hitInfo<<endl;
      } else {
        hitInfo = new UHitInfo();
        fHitInfoArray->Add((TObjArray*)hitInfo);
        faxayToInd[indx][indy] = ind;
        hitInfo->SetIJ(indx, indy);
        // cout<<"ind = "<< ind <<endl;
        ind++;
      }
      // cout<<"track 11  "<<endl;
      if (hitInfo->GetEntries() > 98) {
        cout << "RecoWithSimultaneousHT: Too many Entries: " << hitInfo->GetEntries() << endl;
        break;
      }

      hitInfo->AddHitInfo(i, itab);

      if (itab == 0 || itab == 1) hitInfo->SetVds1Vds2Hit();

      hitsInfo1++;

      // cout<<"ind="<<ind<<" OrigZ="<<OrigZ<<" indx="<<indx<<" indy="<<indy<<endl;
    }
  }

  // if(fMakeClusters) MakeClusters(3,tabArr);

  // if(fMakeClusters) MakeClusters(2,tabArr);

  if (fMakeClusters) MakeClustersNew();

  // methods that should improve effiency

  // now loop over reco clasters in Hough space (tracks)
  // for selected points make iteration tuning of track parameter
  // cout<<"maxXInd="<<maxXInd<<endl;

  // AnaClustersAndCreateTracks(table,tabArr);

  AnaClustersAndCreateTracksNew(table, vdsArr);

  if (!fProduction) FillHSaxayWithClusters();

  // add table with reco tracks

  double ax_vdreco = 0;
  double ay_vdreco = 0;

  for (int i = 0; i < table->GetEntries(); i++) {
    UVdTrack* vdtrack = (UVdTrack*)table->At(i);

    // cout<<"vd track : "<<i<<endl;

    if (!fProduction) {
      ax_vdreco = vdtrack->GetDX() / vdtrack->GetDZ();
      ay_vdreco = vdtrack->GetDY() / vdtrack->GetDZ();
      fGr_axay_vdtrack->SetPoint(i, ax_vdreco, ay_vdreco);
    }
    // else{ cout<<"Tracks not in range for HT Module "<<endl;
    //}
  }

  // cout<<"adding table: "<<table->GetName()<<" Entries: "<<table->GetEntries()<<endl;
  outNode->AddDataTable(table);

  fAllHtTracks = table->GetEntries();

  fhRecoHtTracks->Fill(fAllHtTracks);
}

//_____________________________________________________________________
TObjArray* Na61VdTrackingHTModule::GetHits(TString vdsStr, UEventNode* inJ, UEventNode* inS) {
  int stationID = -1;
  if (vdsStr.Contains("Vds1")) stationID = 1;
  if (vdsStr.Contains("Vds2")) stationID = 2;
  if (vdsStr.Contains("Vds3")) stationID = 3;
  if (vdsStr.Contains("Vds4")) stationID = 4;

  TObjArray* hits = new TObjArray();
  hits->Clear();

  for (int it = 0; it < 8; it++) {
    if (!fSensorNames[it].Contains(vdsStr.Data())) continue;

    // cout<<fSensorNames[it].Data()<<"  "<<vdsStr.Data()<<endl;
    UDataTable* hitstabJ = inJ->GetDataTable(Form("Hits %s", fSensorNames[it].Data()));

    for (int i = 0; i < hitstabJ->GetEntries(); i++) {
      USensorHit* hit = (USensorHit*)hitstabJ->At(i);
      hits->Add(hit);
      hit->SetStationID(stationID);
    }

    UDataTable* hitstabS = inS->GetDataTable(Form("Hits %s", fSensorNames[it].Data()));
    for (int i = 0; i < hitstabS->GetEntries(); i++) {
      USensorHit* hit = (USensorHit*)hitstabS->At(i);
      hits->Add(hit);
      hit->SetStationID(stationID);
    }
  }

  return hits;
}

//_____________________________________________________________________
void Na61VdTrackingHTModule::AnaClustersAndCreateTracksNew(UDataTable* table, TObjArray** tabArr) {
  // New method developed to handle frames with included fake hits
  // Algorithm:
  // 1. Analyse only clisters with N>3 (that contain at least 3 hits)
  // 2. Check whether there are multi-station hits.
  // 3. If multi-station hits occur choose that one that fits better to the
  // rest of hits. Use chi2 ctiterion.

  // cout<<" number of clusters: "<<fHitInfoArray->GetEntries()<<endl;

  for (int i = 0; i < fHitInfoArray->GetLast() + 1; i++) {
    UHitInfo* hitInfo = (UHitInfo*)fHitInfoArray->At(i);

    if (!hitInfo) continue;

    int N = hitInfo->GetEntries();

    if (N < 3) continue;

    // cout<<"=================================>new cluster: i="<<i<<"  N="<<N<<endl;

    if (!hitInfo->ContainsVds1Vds2Hit()) continue;

    // if(N==4)cout<<"=================================>2 new cluster: i="<<i<<"  N="<<N<<endl;

    TObjArray hitsVds1;
    TObjArray hitsVds2;
    TObjArray hitsVds3;
    TObjArray hitsVds4;
    hitsVds1.Clear();
    hitsVds2.Clear();
    hitsVds3.Clear();
    hitsVds4.Clear();

    //    bool containsSignal=false;

    for (int ihit = 0; ihit < N; ihit++) {
      int itab = hitInfo->GetTabId(ihit);
      int ihid = hitInfo->GetIndex(ihit);  // hit location in DataTable array

      USensorHit* hit = (USensorHit*)tabArr[itab]->At(ihid);

      // cout<<" StationID="<<hit->GetStationID()<<"  z="<<hit->GetZ()<<endl;

      if (hit->GetStationID() == 1) hitsVds1.Add(hit);
      if (hit->GetStationID() == 2) hitsVds2.Add(hit);
      if (hit->GetStationID() == 3) hitsVds3.Add(hit);
      if (hit->GetStationID() == 4) hitsVds4.Add(hit);
    }

    int NStations = 0;
    if (hitsVds1.GetEntries()) NStations++;
    if (hitsVds2.GetEntries()) NStations++;
    if (hitsVds3.GetEntries()) NStations++;
    if (hitsVds4.GetEntries()) NStations++;

    // cout<<"======> NStations: "<<NStations<<" "<<hitsVds1.GetEntries()<<" "<<hitsVds2.GetEntries()<<" "<<hitsVds3.GetEntries()<<" "<<hitsVds4.GetEntries()<<endl;

    if (NStations < 3) continue;  // required at lest 3 stations to be fired

    // cout<<" ========> extracting tracks i="<<i<<" NStations: "<<NStations<<endl;
    ExtractTracksFromCluster(hitsVds1, hitsVds2, hitsVds3, hitsVds4, table);
    // cout<<" ========> extracted, tracks: "<<table->GetEntries()<<endl;

  }  // loop over hit infor
}

//_________________________________________________________________________________
void Na61VdTrackingHTModule::ExtractTracksFromCluster(TObjArray& hitsVds1, TObjArray& hitsVds2, TObjArray& hitsVds3, TObjArray& hitsVds4, UDataTable* table) {
  TObjArray clusterTracks;
  clusterTracks.Clear();

  USensorHit* hit1 = 0;
  USensorHit* hit2 = 0;
  USensorHit* hit3 = 0;
  USensorHit* hit4 = 0;

  // now extract from cluster
  int loopIterVds1 = TMath::Max(1, hitsVds1.GetEntries());
  int loopIterVds2 = TMath::Max(1, hitsVds2.GetEntries());
  int loopIterVds3 = TMath::Max(1, hitsVds3.GetEntries());
  int loopIterVds4 = TMath::Max(1, hitsVds4.GetEntries());

  for (int i = 0; i < loopIterVds1; i++) {
    hit1 = (USensorHit*)hitsVds1.At(i);
    for (int j = 0; j < loopIterVds2; j++) {
      hit2 = (USensorHit*)hitsVds2.At(j);
      for (int k = 0; k < loopIterVds3; k++) {
        hit3 = (USensorHit*)hitsVds3.At(k);
        for (int l = 0; l < loopIterVds4; l++) {
          hit4 = (USensorHit*)hitsVds4.At(l);

          // if(hit1)cout<<" z1="<<hit1->GetZ();
          // if(hit2)cout<<" z2="<<hit2->GetZ();
          // if(hit3)cout<<" z3="<<hit3->GetZ();
          // if(hit4)cout<<" z4="<<hit4->GetZ()<<endl;

          double ax, ay, bx, by;
          double chi2x, chi2y, N;
          double zmin;

          // FitLine(hit1,hit2,hit3,hit4,ax,ay,bx,by,chi2x,chi2y,N,zmin);
          double sig[4] = {0.005, 0.005, 0.005, 0.005};  // errors in mm
          USensorHit* hits[4] = {hit1, hit2, hit3, hit4};

          FitLine_w2(hits, sig, ax, ay, bx, by, chi2x, chi2y, N, zmin);

          double chi2xndf = chi2x / (N - 2);
          if (N == 3) fhChi2xNdf_N3->Fill(chi2xndf);
          if (N == 4) fhChi2xNdf_N4->Fill(chi2xndf);

          double chi2yndf = chi2y / (N - 2);
          if (N == 3) fhChi2yNdf_N3->Fill(chi2yndf);
          if (N == 4) fhChi2yNdf_N4->Fill(chi2yndf);

          // double chi2ndf = TMath::Sqrt(chi2xndf*chi2xndf + chi2yndf*chi2yndf);
          double chi2ndf = chi2yndf;

          fhChi2Ndf->Fill(chi2ndf);
          fhChi2xNdf->Fill(chi2xndf);
          fhChi2yNdf->Fill(chi2yndf);

          // cout<<"chi2xndf="<<chi2xndf<<"  chi2yndf="<<chi2yndf<<"  fChi2Cut="<<fChi2Cut<<endl;

          // if(chi2xndf>fChi2Cut || chi2yndf>fChi2Cut) continue;

          // New stuff to improve this cuts
          // if N=3 calculate deviation
          // if N=4 make parabola fit and calculate chi2 or/and hits residual deviations

          double devx = 1111;
          double devy = 1111;
          int id = 0;
          if (N == 3) {
            if (!hit1) id = 0;
            if (!hit2) id = 1;
            if (!hit3) id = 2;
            if (!hit4) id = 3;
            CalculateDeviations(hit1, hit2, hit3, hit4, devx, devy);
            fhDevX[id]->Fill(devx);
            fhDevY[id]->Fill(devy);
          }
          double resx[4] = {1111, 1111, 1111, 1111};
          double resy[4] = {1111, 1111, 1111, 1111};
          if (N == 4) {
            CalculateResDeviations(hits, sig, resx, resy);
            for (int ir = 0; ir < 4; ir++) {
              fhResX[ir]->Fill(resx[ir]);
              fhResY[ir]->Fill(resy[ir]);
            }
          }

          double ddx = 0;
          double ddy = 0;
          // cout<<id<<" "<<fVdParams->GetDevSigx(id)<<" "<<fVdParams->GetResSigx(0)<<" "<<fVdParams->GetResSigx(1)<<fVdParams->GetResSigx(2)<<" "<<fVdParams->GetResSigx(3)<<endl;
          if (N == 3) {
            ddx = ((devx-fVdParams->GetDevOffx(id)) * (devx-fVdParams->GetDevOffx(id))) / (fVdParams->GetDevSigx(id) * fVdParams->GetDevSigx(id));
            ddy = ((devy-fVdParams->GetDevOffy(id)) * (devy-fVdParams->GetDevOffy(id))) / (fVdParams->GetDevSigy(id) * fVdParams->GetDevSigy(id));
          } else {
            for (int is = 0; is < 4; is++) {
              ddx += ((resx[is]-fVdParams->GetResOffx(is)) * (resx[is]-fVdParams->GetResOffx(is))) / (fVdParams->GetResSigx(is) * fVdParams->GetResSigx(is));
              ddy += ((resy[is]-fVdParams->GetResOffy(is)) * (resy[is]-fVdParams->GetResOffy(is))) / (fVdParams->GetResSigy(is) * fVdParams->GetResSigy(is));
            }
          }

          double Nsigy = 6;
          double Nsigx = 10;

          if (ddy > Nsigy * Nsigy) continue;
          if (ddx > Nsigx * Nsigx) continue;
          // if(chi2ndf>fChi2Cut) continue;

          if (N == 3) {
            if (!hit1) {
              fhDevX_acce[0]->Fill(devx);
              fhDevY_acce[0]->Fill(devy);
            }
            if (!hit2) {
              fhDevX_acce[1]->Fill(devx);
              fhDevY_acce[1]->Fill(devy);
            }
            if (!hit3) {
              fhDevX_acce[2]->Fill(devx);
              fhDevY_acce[2]->Fill(devy);
            }
            if (!hit4) {
              fhDevX_acce[3]->Fill(devx);
              fhDevY_acce[3]->Fill(devy);
            }
          }
          if (N == 4) {
            for (int ir = 0; ir < 4; ir++) {
              fhResX_acce[ir]->Fill(resx[ir]);
              fhResY_acce[ir]->Fill(resy[ir]);
            }
          }

          fhChi2Ndf_Acce->Fill(chi2ndf);
          fhChi2xNdf_Acce->Fill(chi2xndf);
          fhChi2yNdf_Acce->Fill(chi2yndf);

          // cout<<"creating track with zmin="<<zmin<<endl;

          // setup tracks

          Vector3D origin;
          origin.SetX(ax * zmin + bx);
          origin.SetY(ay * zmin + by);
          origin.SetZ(zmin);
          // to zmin set to 0 to be consistent with VdTrackingMuke
          // origin.SetX(bx);
          // origin.SetY(by);
          // origin.SetZ(0.);

          Vector3D direction;
          direction.SetX(ax);
          direction.SetY(ay);
          direction.SetZ(1.);

          UVdTrack* track = new UVdTrack(origin, direction);
          track->SetTrackID(1);  // HT track
          track->SetVdsHitIDs(hit1, hit2, hit3, hit4);
          track->MarkForRemoval(false);
          int ind1 = -1;
          int ind2 = -1;
          int ind3 = -1;
          int ind4 = -4;
          if (hit1) ind1 = i;
          if (hit2) ind2 = j;
          if (hit3) ind3 = k;
          if (hit4) ind4 = l;
          // cout<<ind1<<" "<<ind2<<" "<<ind3<<" "<<ind4<<endl;
          track->SetHitArrayIndexes(ind1, ind2, ind3, ind4);

          //	  USensorHit* firstHit=0;
          //	  if(hit1)firstHit=hit1;
          //	  else firstHit=hit2; // add geant4 information for analysis purpose

          track->SetChi2Ndf(chi2ndf);
          // track->SetPdgID(firstHit->GetPdgID());
          // track->SetTrackID(firstHit->GetTrackID());
          // track->SetParentPdgID(firstHit->GetParentPdgID());
          // track->SetParentTrackID(firstHit->GetParentTrackID());

          // if(TMath::Abs(hitInGrid->GetParentPdgID())==421)cout<<"parentPdg="<<track->GetParentPdgID()<<endl;
          // table->Add(track);
          clusterTracks.Add(track);
          // cout<<"i="<<i<<" j="<<j<<" k="<<k<<" l="<<l<<" trackID="<<track->GetTrackID()<<" chi2ndf="<<chi2ndf;
          // cout<<" trID1="<<trID1<<" trID2="<<trID2<<" trID3="<<trID3<<" trID4="<<trID4<<endl;
        }
      }
    }
  }

  // Now we can clean tracks in one cluster. From two or more racks sharing same first hit
  // only the one with the smallest chi2 is kept.

  for (int i = 0; i < clusterTracks.GetEntries(); i++) {
    UVdTrack* track1 = (UVdTrack*)clusterTracks.At(i);
    for (int j = i + 1; j < clusterTracks.GetEntries(); j++) {
      UVdTrack* track2 = (UVdTrack*)clusterTracks.At(j);

      CheckTrackMatching(track1, track2);
    }
  }

  // Add track that are not marked, kill marked track to avoid memory leak
  for (int i = 0; i < clusterTracks.GetEntries(); i++) {
    UVdTrack* track = (UVdTrack*)clusterTracks.At(i);

    if (!track->IsMarkedForRemoval()) {
      table->Add(track);
    } else {
      // track->Delete();
      delete track;
    }
  }

  clusterTracks.Clear();
}

//__________________________________________________________________________________________
void Na61VdTrackingHTModule::CalculateDeviations(USensorHit* hit1, USensorHit* hit2, USensorHit* hit3, USensorHit* hit4, double& devx, double& devy) {
  double x[3];
  double y[3];
  double z[3];

  USensorHit* hits[4] = {hit1, hit2, hit3, hit4};

  int ii = 0;
  for (int i = 0; i < 4; i++) {
    USensorHit* hit = hits[i];
    if (!hit) continue;
    x[ii] = hit->GetX();
    y[ii] = hit->GetY();
    z[ii] = hit->GetZ();

    ii++;
  }

  if (ii != 3) Info(0, "CalcutaleDeviations", Form("somthing is wrong chack it out: ii should be equal 3 but is %d", ii));

  devx = ((z[1] - z[0]) * x[2] + (z[2] - z[1]) * x[0]) / (z[2] - z[0]) - x[1];
  devy = ((z[1] - z[0]) * y[2] + (z[2] - z[1]) * y[0]) / (z[2] - z[0]) - y[1];
}

//__________________________________________________________________________________________
void Na61VdTrackingHTModule::CalculateResDeviations(USensorHit** hits, double* sig, double* resx, double* resy) {
  double ay;
  double by;
  double N;
  double zmin;

  FitLineY_w2(hits, sig, ay, by, N, zmin);

  double x[4];
  double y[4];
  double ex[4];
  double ey[4];

  int ii = 0;
  for (int i = 0; i < 4; i++) {
    if (!hits[i]) continue;

    x[ii] = hits[i]->GetZ();
    ex[ii] = 0.001;  // 1 micron error on z position

    y[ii] = hits[i]->GetX();
    ey[ii] = sig[i];

    ii++;
  }

  if (ii != 4) Info(0, "CalcutaleResDeviations", Form("somthing is wrong chack it out: ii should be equal 4 but is %d", ii));

  TGraphErrors* gr = new TGraphErrors(ii, x, y, ex, ey);

  TF1* pol2x = new TF1("pol2x", "[0] + [1]*x + [2]*x*x", -10.0, 160.0);
  pol2x->SetParameter(0, 0);
  pol2x->SetParameter(1, 0);
  if (GetFieldRun())
    pol2x->SetParameter(2, 0);
  else
    pol2x->FixParameter(2, 0);

  gr->Fit("pol2x", "Q", "", x[0] - 1, x[ii - 1] + 1);

  for (int i = 0; i < 4; i++) {
    double z = hits[i]->GetZ();
    resx[i] = hits[i]->GetX() - pol2x->Eval(z);
    resy[i] = hits[i]->GetY() - (ay * z + by);
  }

  delete pol2x;
  delete gr;
}

//_____________________________________________________________
void Na61VdTrackingHTModule::CheckTrackMatching(UVdTrack* track1, UVdTrack* track2) {
  // TWO or more same hits in the track leads to discarding one track

  int isame = 0;
  for (int i = 0; i < 4; i++) {
    if ((track1->GetHitIdAtStation(i) == track2->GetHitIdAtStation(i)) && track1->GetHitIdAtStation(i)) isame++;
    // cout<<"vds: "<<i<<"  ID1="<<track1->GetHitIdAtStation(i)<<" ID2="<<track2->GetHitIdAtStation(i)<<endl;
  }

  // if(isame>1)cout<<"found same tracks isame="<<isame<<endl;

  if (isame < 2) return;

  // make one track for removal
  if (track1->GetChi2Ndf() > track2->GetChi2Ndf()) {
    // if(track1->IsMarkedForRemoval())cout<<"track1 already marked"<<endl;
    track1->MarkForRemoval(true);  // cout<<"marked track1: "<<track1<<endl;
  } else {
    // if(track2->IsMarkedForRemoval())cout<<"track2 already marked"<<endl;
    track2->MarkForRemoval(true);  // cout<<"marked track2: "<<track2<<endl;
  }
}

//_____________________________________________________________________________
void Na61VdTrackingHTModule::FillHSaxayWithClusters() {
  for (int i = 0; i < findx; i++) {
    for (int j = 0; j < findy; j++) {  // dont loop over side area
      int ind = faxayToInd[i][j];

      if (ind == -1) continue;

      UHitInfo* hitInfo = (UHitInfo*)fHitInfoArray->At(ind);
      if (!hitInfo) Error("MakeClustersNew", "something is wrong, check it out");

      if (hitInfo->GetEntries() == 0) continue;  // means that hitInfo was already used
      fhHSaxay_clust->SetBinContent(i + 1, j + 1, hitInfo->GetEntries());
    }
  }
}

//_____________________________________________________________________
void Na61VdTrackingHTModule::MakeClustersNew() {
  // improved algorithm. It runs independently of first seed value.
  // it make cluster from all that are touching. Size of final cluster is not
  // restricted

  int N = fHitInfoArray->GetEntries();
  // cout<<"N="<<N<<endl;

  // for(int i=1;i<findx-1;i++){
  // for(int j=1;j<findy-1;j++){ // dont loop over side area
  //  int ind = faxayToInd[i][j];

  for (int ind = 0; ind < N; ind++) {
    // if(ind==-1)continue;

    UHitInfo* hitInfo = (UHitInfo*)fHitInfoArray->At(ind);
    if (!hitInfo) Error("MakeClustersNew", "something is wrong, check it out");

    if (hitInfo->GetEntries() == 0) continue;  // means that hitInfo was already used

    fNRes = 0;
    fCluster.Clear();
    fCluster.Add(hitInfo);  // cluster contains one hitInfo
    // cout<<"----> Start: fCluster entries: "<<fCluster.GetEntries()<<" HitInfo Multi="<<hitInfo->GetEntries()<<endl;

    int i = hitInfo->GetIndx();
    int j = hitInfo->GetIndy();

    CheckRing(i, j);

    // we are here means that cluster is generated
    // characterise the cluster and
    // reset global variables
    // cout<<"----> End: fCluster entries: "<<fCluster.GetEntries()<<endl;
    int M = 0;
    for (int ic = 0; ic < fCluster.GetEntries(); ic++) {
      UHitInfo* chitInfo = (UHitInfo*)fCluster.At(ic);
      if (!chitInfo) Error("MakeClustersNew", "something is wrong with your cluster, check it out");
      M = M + chitInfo->GetEntries();
    }

    // if(M<3) continue;

    // store hit infos in one object
    UHitInfo* hitInfoMain = (UHitInfo*)fCluster.At(0);
    // cout<<"pointers: "<<(long long)hitInfo<<" "<<(long long)hitInfoMain<<" Entries:"<<hitInfoMain->GetEntries()<<endl;
    for (int ic = 1; ic < fCluster.GetEntries(); ic++) {
      UHitInfo* chitInfo = (UHitInfo*)fCluster.At(ic);

      for (int ih = 0; ih < chitInfo->GetEntries(); ih++) {
        if (hitInfoMain->GetEntries() > 198) {
          cout << "too many entries in hit info: " << hitInfoMain->GetEntries() << endl;
          return;
        }
        hitInfoMain->AddHitInfo(chitInfo->GetIndex(ih), chitInfo->GetTabId(ih));
      }

      chitInfo->Reset();
    }
    hitInfoMain->SetClusterStatus(1);
    hitInfoMain->SetVds1Vds2Hit();

    if (M != hitInfoMain->GetEntries()) {
      Error("MakeClustersNew", "something is really wrong with your cluster, check it out");
      return;
    }
  }
}

//_____________________________________________________________________
void Na61VdTrackingHTModule::CheckRing(int i, int j) {
  for (int ik = i - 1; ik < i + 2; ik++) {
    for (int jk = j - 1; jk < j + 2; jk++) {
      if (ik == i && jk == j) continue;  // local restriction
      if (ik < 0 || ik > 999) continue;
      if (jk < 0 || jk > 999) continue;

      // check inherited restriction
      bool restricted = false;
      for (int ir = 0; ir < fNRes; ir++) {
        if (TMath::Abs(fResIndx[ir] - ik) < 2 && TMath::Abs(fResIndy[ir] - jk) < 2) restricted = true;
      }

      if (restricted) continue;
      // cout<<"ik="<<ik<<" jk="<<jk<<endl;
      int ind = faxayToInd[ik][jk];
      if (ind == -1) continue;

      UHitInfo* hitInfo = (UHitInfo*)fHitInfoArray->At(ind);
      if (fCluster.GetEntries() > 99) {
        cout << "cluster is very big. Stop with this cluster" << endl;
        return;
      }
      fCluster.Add(hitInfo);  // cluster contains one hitInfo
      // cout<<"CheckRing: ind="<<ind<< " Multi="<<hitInfo->GetEntries()<<endl;
      faxayToInd[ik][jk] = -1;  // don't consider this cluster any longer

      // make restrictions
      fResIndx[fNRes] = i;
      fResIndy[fNRes] = j;
      fNRes++;
      // if(fNRes>19)Info(0,"CheckRing",Form("fNRes is not supposed to be as big >19, check it out fNRes=%d , will continue up to 49",fNRes));
      if (fNRes > 19) {
        // Error("CheckRing","fNRes is not supposed to be as big >19, check it out");
        return;
      }
      CheckRing(ik, jk);
    }
  }
}

//_____________________________________________________________
void Na61VdTrackingHTModule::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61VdTrackingHTModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);
}

//____________________________________________________________________
void Na61VdTrackingHTModule::Print(Option_t* option) const {
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>

  TString opt(option);
  opt.ToLower();

  Na61Module::Print(option);
  if (opt.Contains("d")) cout << endl << "  Original author: Paweł Staszel" << endl << "  Last Modifications: " << endl << "    $Author: Staszel $" << endl << "    $Date: 2012/01/10$" << endl << "    $Revision: 1.0 $ " << endl << endl << "-------------------------------------------------" << endl;
}
