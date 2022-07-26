//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61VdTrackingPostModule
#include "Na61VdTrackingPostModule.h"    
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
//ClassImp(Na61VdTrackingPostModule);

//____________________________________________________________________
Na61VdTrackingPostModule::Na61VdTrackingPostModule() : Na61VdTrackingModule()
{
  // Default constructor. DO NOT USE
}

//____________________________________________________________________
Na61VdTrackingPostModule::Na61VdTrackingPostModule(const Char_t* name, const Char_t* title)
   : Na61VdTrackingModule(name, title)
{
  // Named Constructor
}
//____________________________________________________________________
void Na61VdTrackingPostModule::DefineHistograms()
{

  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init"); 
    return;  
    
  }
  
  // define histos in base class
  Na61VdTrackingModule::DefineHistograms(); 

  TDirectory* histDir = gDirectory->GetDirectory(fHistDirName.Data());
  
  histDir->cd(); 
  
  fhFullTracksPost = new TH1F(Form("hFullTracksPost"),"", 1000,0.,1000.);

 
  
  gDirectory->cd("..");  
}

//____________________________________________________________________
void Na61VdTrackingPostModule::Init()
{
  // Job-level initialisation
  Na61VdTrackingModule::Init();

  
  
}

//____________________________________________________________________
void Na61VdTrackingPostModule::Begin()
{
  // Run-level initialisation
  Na61VdTrackingModule::Begin();
  
}

//____________________________________________________________________
void Na61VdTrackingPostModule :: Event(UEventNode* inNode, UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);  

  UDataTable* hits_Vds1_0 = outNode->GetDataTable(Form("Hits %s",fSensorNames[0].Data()));
  if(!hits_Vds1_0) hits_Vds1_0 = inNode->GetDataTable(Form("Hits %s",fSensorNames[0].Data()));
  UDataTable* hits_Vds2_0 = outNode->GetDataTable(Form("Hits %s",fSensorNames[1].Data()));
  if(!hits_Vds2_0) hits_Vds2_0 = inNode->GetDataTable(Form("Hits %s",fSensorNames[1].Data()));
  UDataTable* hits_Vds3_0 = outNode->GetDataTable(Form("Hits %s",fSensorNames[2].Data()));
  if(!hits_Vds3_0) hits_Vds3_0 = inNode->GetDataTable(Form("Hits %s",fSensorNames[2].Data()));
  UDataTable* hits_Vds3_1 = outNode->GetDataTable(Form("Hits %s",fSensorNames[3].Data()));
  if(!hits_Vds3_1) hits_Vds3_1 = inNode->GetDataTable(Form("Hits %s",fSensorNames[3].Data()));
  UDataTable* hits_Vds4_0 = outNode->GetDataTable(Form("Hits %s",fSensorNames[4].Data()));
  if(!hits_Vds4_0) hits_Vds4_0 = inNode->GetDataTable(Form("Hits %s",fSensorNames[4].Data()));
  UDataTable* hits_Vds4_1 = outNode->GetDataTable(Form("Hits %s",fSensorNames[5].Data()));
  if(!hits_Vds4_1) hits_Vds4_1 = inNode->GetDataTable(Form("Hits %s",fSensorNames[5].Data()));
  UDataTable* hits_Vds4_2 = outNode->GetDataTable(Form("Hits %s",fSensorNames[6].Data()));
  if(!hits_Vds4_2) hits_Vds4_2 = inNode->GetDataTable(Form("Hits %s",fSensorNames[6].Data()));
  UDataTable* hits_Vds4_3 = outNode->GetDataTable(Form("Hits %s",fSensorNames[7].Data()));
  if(!hits_Vds4_3) hits_Vds4_3 = inNode->GetDataTable(Form("Hits %s",fSensorNames[7].Data()));

  fHitsTabs[0] = hits_Vds1_0;  
  fHitsTabs[1] = hits_Vds2_0;  
  fHitsTabs[2] = hits_Vds3_0;  
  fHitsTabs[3] = hits_Vds3_1;  
  fHitsTabs[4] = hits_Vds4_0;  
  fHitsTabs[5] = hits_Vds4_1;  
  fHitsTabs[6] = hits_Vds4_2;  
  fHitsTabs[7] = hits_Vds4_3;  

  
  fhFullTracksAll->Fill(GetAllTracks(inNode));

  Int_t pvtx_status = ((UVdEvent*)inNode)->GetPrimaryVertexStatus();



  if(pvtx_status == 3){
    // pvtx_status==3 means that vertex was set using information form the second arm (jura for saleve and saleve for jura)
    
    SetPrimaryVertex(((UVdEvent*)inNode)->GetPrimaryVertexX(),((UVdEvent*)inNode)->GetPrimaryVertexY(),
		     ((UVdEvent*)inNode)->GetPrimaryVertexZ());
    
    
    
    Int_t tabInd[3]; 
    

    tabInd[0]=1; tabInd[1]=2; tabInd[2]=4;
    Make3HitTracks("dev1_0_x1",0, hits_Vds2_0, hits_Vds3_0, hits_Vds4_0,inNode,tabInd);
    CombineWithPrimaryVertex("down1_x1",inNode);
    tabInd[0]=0; tabInd[1]=2; tabInd[2]=4;
    Make3HitTracks("dev1_0_x2",0, hits_Vds1_0, hits_Vds3_0, hits_Vds4_0,inNode,tabInd);
    CombineWithPrimaryVertex("down1_x2",inNode);
    tabInd[0]=0; tabInd[1]=1; tabInd[2]=4;
    Make3HitTracks("dev1_0_x3",0, hits_Vds1_0, hits_Vds2_0, hits_Vds4_0,inNode,tabInd);
    CombineWithPrimaryVertex("down1_x3",inNode);
    tabInd[0]=0; tabInd[1]=1; tabInd[2]=2;
    Make3HitTracks("dev1_0_x4",0, hits_Vds1_0, hits_Vds2_0, hits_Vds3_0,inNode,tabInd);
    CombineWithPrimaryVertex("down1_x4",inNode);
    
    tabInd[0]=1; tabInd[1]=3; tabInd[2]=5;
    Make3HitTracks("dev1_1_x1",0, hits_Vds2_0, hits_Vds3_1, hits_Vds4_1,inNode,tabInd);
    CombineWithPrimaryVertex("up1_x1",inNode);
    tabInd[0]=0; tabInd[1]=3; tabInd[2]=5;
    Make3HitTracks("dev1_1_x2",0, hits_Vds1_0, hits_Vds3_1, hits_Vds4_1,inNode,tabInd);
    CombineWithPrimaryVertex("up1_x2",inNode);
    tabInd[0]=0; tabInd[1]=1; tabInd[2]=5;
    Make3HitTracks("dev1_1_x3",0, hits_Vds1_0, hits_Vds2_0, hits_Vds4_1,inNode,tabInd);
    CombineWithPrimaryVertex("up1_x3",inNode);
    tabInd[0]=0; tabInd[1]=1; tabInd[2]=3;
    Make3HitTracks("dev1_1_x4",0, hits_Vds1_0, hits_Vds2_0, hits_Vds3_1,inNode,tabInd);
    CombineWithPrimaryVertex("up1_x4",inNode);
    
    // Vds4_2 Jura arm has no data
    if(!fJura || (fRunId>600)){
      tabInd[0]=1; tabInd[1]=2; tabInd[2]=6;
      Make3HitTracks("dev1_2_x1",0, hits_Vds2_0, hits_Vds3_0, hits_Vds4_2,inNode,tabInd);
      CombineWithPrimaryVertex("down2_x1",inNode);
      tabInd[0]=0; tabInd[1]=2; tabInd[2]=6;
      Make3HitTracks("dev1_2_x2",0, hits_Vds1_0, hits_Vds3_0, hits_Vds4_2,inNode,tabInd);
      CombineWithPrimaryVertex("down2_x2",inNode);
      tabInd[0]=0; tabInd[1]=1; tabInd[2]=6;
      Make3HitTracks("dev1_2_x3",0, hits_Vds1_0, hits_Vds2_0, hits_Vds4_2,inNode,tabInd);
      CombineWithPrimaryVertex("down2_x3",inNode);
    }
    
    tabInd[0]=1; tabInd[1]=3; tabInd[2]=7;
    Make3HitTracks("dev1_3_x1",0, hits_Vds2_0, hits_Vds3_1, hits_Vds4_3,inNode,tabInd);
    CombineWithPrimaryVertex("up2_x1",inNode);
    tabInd[0]=0; tabInd[1]=3; tabInd[2]=7;
    Make3HitTracks("dev1_3_x2",0, hits_Vds1_0, hits_Vds3_1, hits_Vds4_3,inNode,tabInd);
    CombineWithPrimaryVertex("up2_x2",inNode);
    tabInd[0]=0; tabInd[1]=1; tabInd[2]=7;
    Make3HitTracks("dev1_3_x3",0, hits_Vds1_0, hits_Vds2_0, hits_Vds4_3,inNode,tabInd);
    CombineWithPrimaryVertex("up2_x3",inNode);
  } 

  fhFullTracksPost->Fill(GetAllTracks(inNode));

  RefitTracks(18,inNode);

}


//_____________________________________________________________
void Na61VdTrackingPostModule::End()
{
  // Run-level finalisation
  Na61VdTrackingModule::End();
}

//____________________________________________________________________
void Na61VdTrackingPostModule::Finish()
{
  // Job-level finalisation
  Na61VdTrackingModule::Finish();
  
}

//____________________________________________________________________
void Na61VdTrackingPostModule ::Print(Option_t* option) const
{
  // Print module informationd
  // In addition this module defines the Option:
  // <fill in here>
  Na61VdTrackingModule::Print(option); 

}



