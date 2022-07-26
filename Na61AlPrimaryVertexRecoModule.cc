//--------------------------------------------
// Primary vertex reco module
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61AlPrimaryVertexRecoModule
#include "Na61AlPrimaryVertexRecoModule.h"
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
//ClassImp(Na61AlPrimaryVertexRecoModule);

//____________________________________________________________________
Na61AlPrimaryVertexRecoModule::Na61AlPrimaryVertexRecoModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
}

//____________________________________________________________________
Na61AlPrimaryVertexRecoModule::Na61AlPrimaryVertexRecoModule(const char* name, const char* title) : Na61PrimaryVertexRecoModule(name, title) {
  // Named Constructor
  SetState(kSetup);

}

//____________________________________________________________________
void Na61AlPrimaryVertexRecoModule::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  TDirectory* histDir = gDirectory->mkdir("AlPrimaryVertexReco");
  
  histDir->cd();
  
  
  gDirectory->cd("..");

  Na61PrimaryVertexRecoModule::DefineHistograms();

}

//____________________________________________________________________
void Na61AlPrimaryVertexRecoModule::Init() {
  // Job-level initialisation
  SetState(kInit);

  Na61PrimaryVertexRecoModule::Init();

  /*
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
  */

  fNsigS = 12;
  fNsigJ = 14;
  fZprim = 73.1;

}

//____________________________________________________________________
void Na61AlPrimaryVertexRecoModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);
}
/*
//____________________________________________________________________
void Na61AlPrimaryVertexRecoModule::SetRunInfo(TFile* file)
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
void Na61AlPrimaryVertexRecoModule::Event(UEventNode* inNodeJ, UEventNode* inNodeS)

{
  // Per event method
  SetState(kEvent);

  fAllTracksJ = 0;
  fAllTracksS = 0;
 
  FillClusterCorrelationAl(inNodeJ, inNodeS);

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

  if (pvertStatusJ == 1 && pvertStatusS == 1) {
    // if(pvertStatusJ > 0 && pvertStatusS > 0 ){ // be less restricted for h+Pb
    // for MAVD we have: 2*7.5 + 6 = 21
    // for MAVD in pPb2022 test we have: 2*7.5 + 10 = 25
    // Offset is x is arbitrary shifted to be position around 0
    // To find corrections for arms x positions used fhVtxDx_Glob
    fhVtxDx->Fill(pvtxX_J - pvtxX_S + 25 + 2.245);
    fhVtxDy->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz->Fill(pvtxZ_J - pvtxZ_S);
    //cout << "Jura Xpos= " << pvtxX_J << " Saleve Xpos=" << pvtxX_S << " pvtxX_J - pvtxX_S + 25. = " << pvtxX_J - pvtxX_S + 25. << endl;
    //Int_t ii;
    //cin>>ii;
  }


  LocalToGlobalAl(0, inNodeJ);
  LocalToGlobalAl(1, inNodeS);

  FillHitPositions(0, inNodeJ);
  FillHitPositions(1, inNodeS);

  fJuraProdCounter = 0;
  fSaleveProdCounter = 0;

  ClassifyTracksAl(inNodeJ, inNodeS);

  // if(fJuraProdCounter>2 && fSaleveProdCounter>2){
  // cout<<" JuraProdCounter="<<fJuraProdCounter<<" SaleveProdCounter="<<fSaleveProdCounter<<endl;
  //}

  FindPrimaryVertexAl(1, inNodeJ, inNodeS);
  FindPrimaryVertexAl(2, inNodeJ, inNodeS);
  FindPrimaryVertexAl(2, inNodeJ, 0);
  FindPrimaryVertexAl(2, 0, inNodeS);

  pvtxX_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexX();
  pvtxY_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexY();
  pvtxZ_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexZ();

  pvtxX_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexX();
  pvtxY_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexY();
  pvtxZ_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexZ();


  if (pvertStatusJ == 1 && pvertStatusS == 1) {
    // this histogram shall be used for arm geometry correction   
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
void Na61AlPrimaryVertexRecoModule::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61AlPrimaryVertexRecoModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);
}

//____________________________________________________________________
void Na61AlPrimaryVertexRecoModule::Print(Option_t* option) const {
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>

  TString opt(option);
  opt.ToLower();

  Na61Module::Print(option);
  if (opt.Contains("d")) cout << endl << "  Original author: Paweł Staszel" << endl << "  Last Modifications: " << endl << "    $Author: Staszel $" << endl << "    $Date: 2017/03/20$" << endl << "    $Revision: 1.0 $ " << endl << endl << "-------------------------------------------------" << endl;
}
