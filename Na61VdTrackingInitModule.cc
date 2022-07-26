//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61VdTrackingInitModule
#include "Na61VdTrackingInitModule.h"
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
// ClassImp(Na61VdTrackingInitModule);

//____________________________________________________________________
Na61VdTrackingInitModule::Na61VdTrackingInitModule() : Na61VdTrackingModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
}

//____________________________________________________________________
Na61VdTrackingInitModule::Na61VdTrackingInitModule(const char* name, const char* title) : Na61VdTrackingModule(name, title) {
  // Named Constructor
  SetState(kSetup);
}
//____________________________________________________________________
void Na61VdTrackingInitModule::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  // define histos in base class
  Na61VdTrackingModule::DefineHistograms();

  // cout<<"2 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  TDirectory* histDir = gDirectory->GetDirectory(fHistDirName.Data());

  // cout<<"3 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  histDir->cd();

  // cout<<"4 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  gDirectory->cd("..");

  // cout<<"5 fhFullTracksAll: "<<fhFullTracksAll<<endl;
}

//____________________________________________________________________
void Na61VdTrackingInitModule::Init() {
  // Job-level initialisation
  SetState(kInit);

  Na61VdTrackingModule::Init();

  /* moved to base class
  if(fJura)fParams = (Na61VdParametersManager::Instance())->GetJuraArmParams();
  else fParams = (Na61VdParametersManager::Instance())->GetSaleveArmParams();

  cout<<"Na61VdTrackingInitModule::Init(): fNsig_dev="<<fNsig_dev<<"   fNsig_pvert="<<fNsig_pvert<<endl;
  */
}

//____________________________________________________________________
void Na61VdTrackingInitModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);
}

//____________________________________________________________________
void Na61VdTrackingInitModule::Event(UEventNode* inNode, UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);

 //cout<<"6 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  UDataTable* hits_Vds1_0 = outNode->GetDataTable(Form("Hits %s", fSensorNames[0].Data()));

  if (!hits_Vds1_0) hits_Vds1_0 = inNode->GetDataTable(Form("Hits %s", fSensorNames[0].Data()));
  UDataTable* hits_Vds2_0 = outNode->GetDataTable(Form("Hits %s", fSensorNames[1].Data()));
  if (!hits_Vds2_0) hits_Vds2_0 = inNode->GetDataTable(Form("Hits %s", fSensorNames[1].Data()));
  UDataTable* hits_Vds3_0 = outNode->GetDataTable(Form("Hits %s", fSensorNames[2].Data()));
  if (!hits_Vds3_0) hits_Vds3_0 = inNode->GetDataTable(Form("Hits %s", fSensorNames[2].Data()));
  UDataTable* hits_Vds3_1 = outNode->GetDataTable(Form("Hits %s", fSensorNames[3].Data()));
  if (!hits_Vds3_1) hits_Vds3_1 = inNode->GetDataTable(Form("Hits %s", fSensorNames[3].Data()));
  UDataTable* hits_Vds4_0 = outNode->GetDataTable(Form("Hits %s", fSensorNames[4].Data()));
  if (!hits_Vds4_0) hits_Vds4_0 = inNode->GetDataTable(Form("Hits %s", fSensorNames[4].Data()));
  UDataTable* hits_Vds4_1 = outNode->GetDataTable(Form("Hits %s", fSensorNames[5].Data()));
  if (!hits_Vds4_1) hits_Vds4_1 = inNode->GetDataTable(Form("Hits %s", fSensorNames[5].Data()));
  UDataTable* hits_Vds4_2 = outNode->GetDataTable(Form("Hits %s", fSensorNames[6].Data()));
  if (!hits_Vds4_2) hits_Vds4_2 = inNode->GetDataTable(Form("Hits %s", fSensorNames[6].Data()));
  UDataTable* hits_Vds4_3 = outNode->GetDataTable(Form("Hits %s", fSensorNames[7].Data()));
  if (!hits_Vds4_3) hits_Vds4_3 = inNode->GetDataTable(Form("Hits %s", fSensorNames[7].Data()));

  fHitsTabs[0] = hits_Vds1_0;
  fHitsTabs[1] = hits_Vds2_0;
  fHitsTabs[2] = hits_Vds3_0;
  fHitsTabs[3] = hits_Vds3_1;
  fHitsTabs[4] = hits_Vds4_0;
  fHitsTabs[5] = hits_Vds4_1;
  fHitsTabs[6] = hits_Vds4_2;
  fHitsTabs[7] = hits_Vds4_3;

  fPrimaryVertexDefined = false;
  ((UVdEvent*)outNode)->SetPrimaryVertexStatus(0);

  int tabInd[3];
  int minSize = 0;
  // int minSize = 100;

  tabInd[0] = 0;
  tabInd[1] = 1;
  tabInd[2] = 2;
  Make3HitTracks("dev1_0", minSize, hits_Vds1_0, hits_Vds2_0, hits_Vds3_0, outNode, tabInd);

  tabInd[0] = 0;
  tabInd[1] = 1;
  tabInd[2] = 3;
  Make3HitTracks("dev1_1", minSize, hits_Vds1_0, hits_Vds2_0, hits_Vds3_1, outNode, tabInd);

  tabInd[0] = 1;
  tabInd[1] = 2;
  tabInd[2] = 4;
  Make3HitTracks("dev2_0", minSize, hits_Vds2_0, hits_Vds3_0, hits_Vds4_0, outNode, tabInd);

  tabInd[0] = 1;
  tabInd[1] = 3;
  tabInd[2] = 5;
  Make3HitTracks("dev2_1", minSize, hits_Vds2_0, hits_Vds3_1, hits_Vds4_1, outNode, tabInd);

  tabInd[0] = 1;
  tabInd[1] = 2;
  tabInd[2] = 6;
  Make3HitTracks("dev2_2", minSize, hits_Vds2_0, hits_Vds3_0, hits_Vds4_2, outNode, tabInd);

  tabInd[0] = 1;
  tabInd[1] = 3;
  tabInd[2] = 7;
  Make3HitTracks("dev2_3", minSize, hits_Vds2_0, hits_Vds3_1, hits_Vds4_3, outNode, tabInd);

  RemoveGhostTracks("dev2_0", "dev2_2", outNode);
  RemoveGhostTracks("dev2_0", "dev2_3", outNode);
  RemoveGhostTracks("dev2_1", "dev2_2", outNode);
  RemoveGhostTracks("dev2_1", "dev2_3", outNode);

  // make tracks of large clusters (Pb beam tracks)

  // Make3HitTracks("dev1_0",20, hits_Vds1_0, hits_Vds2_0, hits_Vds3_0,outNode);

  // Make3HitTracks("dev1_1",20, hits_Vds1_0, hits_Vds2_0, hits_Vds3_1,outNode);

  // Make3HitTracks("dev2_0",20, hits_Vds2_0, hits_Vds3_0, hits_Vds4_0,outNode);

  // Make3HitTracks("dev2_1",20, hits_Vds2_0, hits_Vds3_1, hits_Vds4_1,outNode);

  // Make3HitTracks("dev2_2",20, hits_Vds2_0, hits_Vds3_0, hits_Vds4_2,outNode);

  // Make3HitTracks("dev2_3",20, hits_Vds2_0, hits_Vds3_1, hits_Vds4_3,outNode);

  CombineTracks("down1", 0, 2, outNode);  // dev1_0 vs dev2_0

  CombineTracks("up1", 1, 3, outNode);  // dev1_1 vs dev2_1

  CombineTracks("down2", 0, 4, outNode);  // dev1_0 vs dev2_2

  CombineTracks("up2", 1, 5, outNode);  // dev1_1 vs dev2_3

  RefitTracks(4, outNode);

  FindPrimaryVertex(outNode);

 // cout<<"7 fhFullTracksAll: "<<fhFullTracksAll<<endl;

  // return;  // dont make post tracking for simplicity.

  if (fPrimaryVertexDefined) {
    // cout<<"71 fhFullTracksAll: "<<fhFullTracksAll<<endl;
    // cout<<"tu1"<<endl;

    tabInd[0] = 1;
    tabInd[1] = 2;
    tabInd[2] = 4;
    Make3HitTracks("dev1_0_x1", 0, hits_Vds2_0, hits_Vds3_0, hits_Vds4_0, outNode, tabInd);
    CombineWithPrimaryVertex("down1_x1", outNode);
    tabInd[0] = 0;
    tabInd[1] = 2;
    tabInd[2] = 4;
    Make3HitTracks("dev1_0_x2", 0, hits_Vds1_0, hits_Vds3_0, hits_Vds4_0, outNode, tabInd);
    CombineWithPrimaryVertex("down1_x2", outNode);
    tabInd[0] = 0;
    tabInd[1] = 1;
    tabInd[2] = 4;
    Make3HitTracks("dev1_0_x3", 0, hits_Vds1_0, hits_Vds2_0, hits_Vds4_0, outNode, tabInd);
    CombineWithPrimaryVertex("down1_x3", outNode);
    tabInd[0] = 0;
    tabInd[1] = 1;
    tabInd[2] = 2;
    Make3HitTracks("dev1_0_x4", 0, hits_Vds1_0, hits_Vds2_0, hits_Vds3_0, outNode, tabInd);
    CombineWithPrimaryVertex("down1_x4", outNode);

    // cout<<"72 fhFullTracksAll: "<<fhFullTracksAll<<endl;

    tabInd[0] = 1;
    tabInd[1] = 3;
    tabInd[2] = 5;
    // cout<<"720 fhFullTracksAll: "<<fhFullTracksAll<<endl;
    Make3HitTracks("dev1_1_x1", 0, hits_Vds2_0, hits_Vds3_1, hits_Vds4_1, outNode, tabInd);
    // cout<<"721 fhFullTracksAll: "<<fhFullTracksAll<<endl;
    CombineWithPrimaryVertex("up1_x1", outNode);
    tabInd[0] = 0;
    tabInd[1] = 3;
    tabInd[2] = 5;
    // cout<<"722 fhFullTracksAll: "<<fhFullTracksAll<<endl;
    Make3HitTracks("dev1_1_x2", 0, hits_Vds1_0, hits_Vds3_1, hits_Vds4_1, outNode, tabInd);
    CombineWithPrimaryVertex("up1_x2", outNode);
    // cout<<"723 fhFullTracksAll: "<<fhFullTracksAll<<endl;
    tabInd[0] = 0;
    tabInd[1] = 1;
    tabInd[2] = 5;
    Make3HitTracks("dev1_1_x3", 0, hits_Vds1_0, hits_Vds2_0, hits_Vds4_1, outNode, tabInd);
    CombineWithPrimaryVertex("up1_x3", outNode);
    // cout<<"724 fhFullTracksAll: "<<fhFullTracksAll<<endl;
    tabInd[0] = 0;
    tabInd[1] = 1;
    tabInd[2] = 3;
    Make3HitTracks("dev1_1_x4", 0, hits_Vds1_0, hits_Vds2_0, hits_Vds3_1, outNode, tabInd);
    CombineWithPrimaryVertex("up1_x4", outNode);

    // cout<<"73 fhFullTracksAll: "<<fhFullTracksAll<<endl;

    // Vds4_2 Jura arm sensor was not working for runid<190
    if (!fJura || fRunId > 190) {
      tabInd[0] = 1;
      tabInd[1] = 2;
      tabInd[2] = 6;
      Make3HitTracks("dev1_2_x1", 0, hits_Vds2_0, hits_Vds3_0, hits_Vds4_2, outNode, tabInd);
      CombineWithPrimaryVertex("down2_x1", outNode);
      tabInd[0] = 0;
      tabInd[1] = 2;
      tabInd[2] = 6;
      Make3HitTracks("dev1_2_x2", 0, hits_Vds1_0, hits_Vds3_0, hits_Vds4_2, outNode, tabInd);
      CombineWithPrimaryVertex("down2_x2", outNode);
      tabInd[0] = 0;
      tabInd[1] = 1;
      tabInd[2] = 6;
      Make3HitTracks("dev1_2_x3", 0, hits_Vds1_0, hits_Vds2_0, hits_Vds4_2, outNode, tabInd);
      CombineWithPrimaryVertex("down2_x3", outNode);
    }

    // cout<<"74 fhFullTracksAll: "<<fhFullTracksAll<<endl;

    tabInd[0] = 1;
    tabInd[1] = 3;
    tabInd[2] = 7;
    Make3HitTracks("dev1_3_x1", 0, hits_Vds2_0, hits_Vds3_1, hits_Vds4_3, outNode, tabInd);
    CombineWithPrimaryVertex("up2_x1", outNode);
    tabInd[0] = 0;
    tabInd[1] = 3;
    tabInd[2] = 7;
    Make3HitTracks("dev1_3_x2", 0, hits_Vds1_0, hits_Vds3_1, hits_Vds4_3, outNode, tabInd);
    CombineWithPrimaryVertex("up2_x2", outNode);
    tabInd[0] = 0;
    tabInd[1] = 1;
    tabInd[2] = 7;
    Make3HitTracks("dev1_3_x3", 0, hits_Vds1_0, hits_Vds2_0, hits_Vds4_3, outNode, tabInd);
    CombineWithPrimaryVertex("up2_x3", outNode);

    // cout<<"8 fhFullTracksAll: "<<fhFullTracksAll<<endl;

    RefitTracks(18, outNode);
    FindPrimaryVertexPost(outNode);
  }

  // cout<<"-------------------------------------> tu4"<<endl;
  // outNode->ListObjects();
  Delete3HitTracks(outNode);
  // cout<<"-------------------------------------> tu5"<<endl;
  // outNode->ListObjects();
}

//_____________________________________________________________
void Na61VdTrackingInitModule::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61VdTrackingInitModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);
}

//____________________________________________________________________
void Na61VdTrackingInitModule::Print(Option_t* option) const {
  // Print module informationd
  // In addition this module defines the Option:
  // <fill in here>

  Na61VdTrackingModule::Print(option);
}
