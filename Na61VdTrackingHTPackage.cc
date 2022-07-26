//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61VdTrackingHTPackage
#include "Na61VdTrackingHTPackage.h"
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
// ClassImp(Na61VdTrackingHTPackage);

//____________________________________________________________________
Na61VdTrackingHTPackage::Na61VdTrackingHTPackage() : Na61VdTrackingModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
}
//____________________________________________________________________
Na61VdTrackingHTPackage::Na61VdTrackingHTPackage(const char* name, const char* title) : Na61VdTrackingModule(name, title) {
  // Named Constructor
  SetState(kSetup);

  fNModules = 5;
  if (fNModules > 20) {
    cout << "ERROR: Na61VdTrackingHTPackage: NModules>20 which is not allowed, Check it out!!!!" << endl;
    return;
  }

  for (int i = 0; i < fNModules; i++) {
    fTrackingHTModules[i] = new Na61VdTrackingHTModule(Form("Vd Tracking HT Module %d", i + 1), Form("Vd Tracking HT Module %d", i + 1));
  }
}

//____________________________________________________________________
void Na61VdTrackingHTPackage::Init() {
  // Job-level initialisation
  SetState(kInit);

  for (int i = 0; i < fNModules; i++) {
    fTrackingHTModules[i]->SetRunId(fRunId);
    fTrackingHTModules[i]->SetProduction(fProduction);
    fTrackingHTModules[i]->SetPosRes(fPosRes);
    fTrackingHTModules[i]->Setbw(fbwx, fbwy);
    fTrackingHTModules[i]->SetMakeClusters(fMakeClusters);
    fTrackingHTModules[i]->SetChi2Cut(fChi2Cut);
    fTrackingHTModules[i]->SetVzOffset(1000 + i * 1500);  // primary vertex shift in um
    fTrackingHTModules[i]->SetHistDirName(Form("VdTrackingHTModule_r%d", i + 1));
    fTrackingHTModules[i]->SetOutTableName(Form("Vd HT Tracks r=%d", i + 1));
  }

  for (int i = 0; i < fNModules; i++) {
    fTrackingHTModules[i]->Init();
  }
}

//____________________________________________________________________
void Na61VdTrackingHTPackage::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  for (int i = 0; i < fNModules; i++) {
    fTrackingHTModules[i]->DefineHistograms();
  }
}

//____________________________________________________________________
void Na61VdTrackingHTPackage::Begin() {
  // Run-level initialisation
  SetState(kBegin);

  // cout<<"track 3 "<<endl;
}

//____________________________________________________________________
void Na61VdTrackingHTPackage::Event(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode) {
  // Per event method

  SetState(kEvent);

   //cout<<"  Na61VdTrackingHTPackage ------------------------------------------------>new event "<<endl;
   //cout<<"fNModules="<<fNModules<<endl;
  for (int i = 0; i < fNModules; i++) {
	//cout<<"i="<<i<<endl;
    fTrackingHTModules[i]->Event(inNodeJ, inNodeS, outNode);
  }

  CleanTracks(outNode);
}
//______________________________________________________________________________________________________
void Na61VdTrackingHTPackage::CleanTracks(UEventNode* out) {
  TString devname = "Vd HT Tracks";
  for (int i = 0; i < fNModules; i++) {
    TString devname2 = Form("Vd HT Tracks r=%d", i + 1);
    RemoveMultiplyTracks(devname, devname2, out);
  }
}

//______________________________________________________________________________________________________
void Na61VdTrackingHTPackage::RemoveMultiplyTracks(TString devname1, TString devname2, UEventNode* out) {
  UDataTable* tracktab1 = out->GetDataTable(Form("%s", devname1.Data()));
  UDataTable* tracktab2 = out->GetDataTable(Form("%s", devname2.Data()));

  if (!tracktab2) return;

  if (!tracktab1) {
    tracktab1 = new UDataTable("Vd HT Tracks");
    tracktab1->SetOwner();
    out->AddDataTable(tracktab1);
  }

  int N1 = tracktab1->GetEntries();
  int N2 = tracktab2->GetEntries();

  // cout<<"---------------> RemoveMultipleTracks: "<<devname2.Data()<<" tracks: "<<N1<<" range tracks: "<<N2<<endl;

  for (int i = 0; i < N1; i++) {
    UVdTrack* track1 = (UVdTrack*)tracktab1->At(i);
    for (int j = 0; j < N2; j++) {
      UVdTrack* track2 = (UVdTrack*)tracktab2->At(j);

      if (track2->IsMarkedForRemoval()) continue;

      if (TracksMatch(track1, track2)) track2->MarkForRemoval(true);

      // if(track2->IsMarkedForRemoval()){
      // cout<<" track1: "<<track1->GetX()<<" "<<track1->GetY()<<" "<<track1->GetZ()<<endl;
      // cout<<" track2: "<<track2->GetX()<<" "<<track2->GetY()<<" "<<track2->GetZ()<<endl;
      //}
    }
  }

  // this is a bit tricky but should work
  for (int i = 0; i < N2; i++) {
    UVdTrack* track2 = (UVdTrack*)tracktab2->At(i);
    if (track2->IsMarkedForRemoval()) {
      tracktab2->DeleteObjectAt(i);
    } else {
      tracktab1->Add(track2);
      tracktab2->RemoveAt(i);
    }
  }

  tracktab2->Compress();
  tracktab2->Delete();
  if (tracktab2) out->RemoveDataObject((UDataObject*)tracktab2);

  // cout<<"final tracks: "<<tracktab1->GetEntries()<<" Added tracks: "<<tracktab1->GetEntries()-N1<<" out of "<<N2<<endl;
}

//_____________________________________________________________
bool Na61VdTrackingHTPackage::TracksMatch(UVdTrack* track1, UVdTrack* track2) {
  // TWO or more same hits in the track leads to discarding one track

  int isame = 0;
  for (int i = 0; i < 4; i++) {
    if ((track1->GetHitIdAtStation(i) == track2->GetHitIdAtStation(i)) && track1->GetHitIdAtStation(i)) isame++;
    // cout<<"i = "<<i<<"   "<<track1->GetHitIdAtStation(i)<<"   "<<track2->GetHitIdAtStation(i)<<endl;
  }

  if (isame < 2)
    return false;
  else
    return true;
}

//_____________________________________________________________
void Na61VdTrackingHTPackage::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61VdTrackingHTPackage::Finish() {
  // Job-level finalisation
  SetState(kFinish);
}

//____________________________________________________________________
void Na61VdTrackingHTPackage::Print(Option_t* option) const {
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>

  TString opt(option);
  opt.ToLower();

  Na61Module::Print(option);
  if (opt.Contains("d")) cout << endl << "  Original author: Paweł Staszel" << endl << "  Last Modifications: " << endl << "    $Author: Staszel $" << endl << "    $Date: 2018/10/08$" << endl << "    $Revision: 1.0 $ " << endl << endl << "-------------------------------------------------" << endl;
}
