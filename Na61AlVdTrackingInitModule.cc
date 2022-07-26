//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61AlVdTrackingInitModule
#include "Na61AlVdTrackingInitModule.h"    
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
#include "iostream"

#include <string>
#include <time.h>
using namespace std;

//____________________________________________________________________
//ClassImp(Na61AlVdTrackingInitModule);

//____________________________________________________________________
Na61AlVdTrackingInitModule::Na61AlVdTrackingInitModule() : Na61AlVdTrackingModule()
{
  // Default constructor. DO NOT USE
  SetState(kSetup);

}

//____________________________________________________________________
Na61AlVdTrackingInitModule::Na61AlVdTrackingInitModule(const Char_t* name, const Char_t* title)
   : Na61AlVdTrackingModule(name, title)
{
  // Named Constructor
  SetState(kSetup);
  fNsig_dev = 4.0;
  fEvents = 0;
 
}
//____________________________________________________________________
void Na61AlVdTrackingInitModule::DefineHistograms()
{

  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init"); 
    return;  
    
  }

  Na61AlVdTrackingModule::DefineHistograms(); 
  
  TDirectory* histDir = gDirectory->GetDirectory(fHistDirName.Data());
  
  histDir->cd(); 
  
  /*
  cout<<" registered sensor combinations: "<<fNDev<<endl;
  Int_t nbins = 5000;
  Float_t xmin = -10.0;
  Float_t xmax =  10.0;
  Int_t i=0;
  // loop when sensor combinations are given by hand (registered)
  for(Int_t i=0;i<fNDev;i++){
    TString devname = Form("dev%s",fDevName[i].Data());
    fhX_dev[i] = new TH1F(Form("hX_%s",devname.Data()),"", nbins, xmin, xmax);
    fhY_dev[i] = new TH1F(Form("hY_%s",devname.Data()),"", nbins, xmin, xmax);
    fhX_dev_cuty[i] = new TH1F(Form("hX_%s_cuty",devname.Data()),"", nbins, xmin, xmax);
    fhY_dev_cutx[i] = new TH1F(Form("hY_%s_cutx",devname.Data()),"", nbins, xmin, xmax);
    fh_devx_devy[i]= new TH2F(Form("h_devx_devy_%s",devname.Data()),"",500,xmin,xmax,500,xmin,xmax);
  }
    
  for(Int_t ii=0;ii<fNDev;ii++)cout<<"ii="<<ii<<"  fDevName[ii]: "<<fDevName[ii].Data()<<endl;
  */

  gDirectory->cd("..");
  
}

//____________________________________________________________________
void Na61AlVdTrackingInitModule::Init()
{
  // Job-level initialisation
  SetState(kInit);

  Na61Module::Init();
  Na61AlVdTrackingModule::Init();

  //fParams = (Na61VdParametersManager::Instance())->GetAtdParams();
  //fParams = (Na61VdParametersManager::Instance())->GetAlVdParams();
  
}

//____________________________________________________________________
void Na61AlVdTrackingInitModule::Begin()
{
  // Run-level initialisation
  SetState(kBegin);
  
}

//____________________________________________________________________
void Na61AlVdTrackingInitModule :: Event(UEventNode* inNode, UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);  

  //cout<<"-----------------------------  Na61AlVdTrackingInitModule :: Event: new #: "<<fEvents<<endl;
  //cout<<"----------------------------- "<<fHistDirName.Data()<<endl;
  fEvents++;

  for(Int_t i=0;i<34;i++){
    fHitsTabs[i] = outNode->GetDataTable(Form("Hits %s",fAlSensorNames[i].Data()));
    if(!fHitsTabs[i])fHitsTabs[i] = inNode->GetDataTable(Form("Hits %s",fAlSensorNames[i].Data()));
  }

  // if geometry was redefined after the hits were created, in the first step the hits shell to be
  // transformed into new geometry.

  for(Int_t i=0;i<34;i++){
    fParams->LocalToGlobal(i,fHitsTabs[i]->GetObjectList());
  }

  Int_t tabInd[3]; 
  // here we use so called topo notation in which a central chip in the stave is numbered as "0"

  FillHitPositions(outNode);

  for(Int_t i=0;i<fNDev;i++){
    TString devname = fDevName[i];
    Int_t vds[3];
    Int_t ic[3];
    Int_t is[3];
    vds[0] = fVdsArr[i][0]; vds[1] = fVdsArr[i][1]; vds[2] = fVdsArr[i][2];
    ic[0]  = fColArr[i][0]; ic[1]  = fColArr[i][1]; ic[2]  = fColArr[i][2];
    is[0]  = fSenArr[i][0]; is[1]  = fSenArr[i][1]; is[2]  = fSenArr[i][2];
    // convert devname to colums and sensorNumb
    tabInd[0]=GetAlSensorNumber2(vds[0],ic[0],is[0]);
    tabInd[1]=GetAlSensorNumber2(vds[1],ic[1],is[1]);
    tabInd[2]=GetAlSensorNumber2(vds[2],ic[2],is[2]);
    //cout<<"Na61AlVdTrackingInitModule :: Event: idev = "<<i<<" "<<devname.Data()<<" "<<tabInd[0]<<" "<<tabInd[1]<<" "<<tabInd[2]<<endl;
    //cout<<vds[0]<<" "<<vds[1]<<" "<<vds[2]<<endl; 
    //cout<<ic[0]<<" "<<ic[1]<<" "<<ic[2]<<endl; 
    //cout<<is[0]<<" "<<is[1]<<" "<<is[2]<<endl; 
    Make3HitTracks(i,outNode,tabInd);
  }


  // here we create 4hit tracks by combining 3hit tracks
  for(Int_t imatch=0; imatch<fNDevComb; imatch++){
    
    Int_t idevi = fDevi[imatch];
    Int_t idevj = fDevj[imatch];
    
    CombineTracks(imatch,idevi,idevj,outNode); // dev1_0 vs dev2_0 
  }

  // fNDevComb - number of fill track tables
  RefitTracks(fNDevComb, outNode);
  
  FindPrimaryVertex(outNode);


}

//_____________________________________________________________
void Na61AlVdTrackingInitModule::End()
{
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61AlVdTrackingInitModule::Finish()
{
  // Job-level finalisation
  SetState(kFinish);
  
}



//____________________________________________________________________
void Na61AlVdTrackingInitModule ::Print(Option_t* option) const
{
  // Print module informationd
  // In addition this module defines the Option:
  // <fill in here>
  
  Na61AlVdTrackingModule ::Print(option);
}



