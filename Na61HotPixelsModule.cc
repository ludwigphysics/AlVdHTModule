//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61HotPixelsModule
#include "Na61HotPixelsModule.h"    
#endif
#ifndef ROOT_TDirectory
#include "TDirectory.h"
#endif
#ifndef ROOT_TSystem
#include "TSystem.h"
#endif

#ifndef UTIL_UDataTable
#include "UDataTable.h"
#endif
#ifndef UTIL_UChipPixel
#include "UChipPixel.h"
#endif

#ifndef ROOT_TF1
#include "TF1.h"
#endif
#ifndef ROOT_TMath
#include "TMath.h"
#endif

#include "TSystem.h"
#include "TDirectory.h"
#include "Math/QuantFuncMathCore.h"
#include <time.h>
#include "iostream"

//____________________________________________________________________
//ClassImp(Na61HotPixelsModule);

//____________________________________________________________________
Na61HotPixelsModule::Na61HotPixelsModule()
{
// Default constructor. DO NOT USE
  SetState(kSetup);
  fEvents=0;
 
}

//____________________________________________________________________
Na61HotPixelsModule::Na61HotPixelsModule(const Char_t* name, const Char_t* title) : Na61Module(name, title) {
  // Named Constructor
  
  SetState(kSetup);
  fEvents=0;
  
}

//____________________________________________________________________
void Na61HotPixelsModule::DefineHistograms()
{
  
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init"); 
    return;    
  }
  
  TDirectory* histDir; 
  
  if(fJura)histDir = gDirectory->mkdir("HotPixelsModule_Jura");
  else histDir = gDirectory->mkdir("HotPixelsModule_Saleve");
  
  histDir->cd(); 
  
  
  // Define histograms. They are:
  // <fill in here>


  Int_t binsy  = findy;
  Int_t binsx  = findx;  
  Float_t ymax =  binsy;
  Float_t xmax =  binsx;

  
  for(Int_t i=0;i<34;i++){
    fhPixels[i] = new TH1F(Form("hPixels_%s",fAlSensorNames[i].Data()),"",100,0,100);
    fhCleanPixels[i] = new TH1F(Form("hCleanPixels_%s",fAlSensorNames[i].Data()),"",100,0,100);
    fhPixelMap[i] = new TH2F(Form("hPixelMap_%s",fAlSensorNames[i].Data()),"", binsx, 0.5, xmax+0.5, binsy, 0.5, ymax+0.5);
    fhNoisyPixelMap[i] = new TH2F(Form("hNoisyPixelMap_%s",fAlSensorNames[i].Data()),"", binsx, 0.5, xmax+0.5, binsy, 0.5, ymax+0.5);
    fhPixelFrequency[i] = new TH1F(Form("hPixelFrequency_%s",fAlSensorNames[i].Data()),"",1000, 0, 100);
  }


  ReadNoisyPixelInfo();
  cout<<" fNoisyPixelsDefined="<<fNoisyPixelsDefined<<endl;

  histDir->cd(); 
  gDirectory->cd("..");
  
}

//____________________________________________________________________
void Na61HotPixelsModule::Init()
{
  // Job-level initialisation
  SetState(kInit);
  
    fNVds=4;
    fVdsZ[0] = 0;
    fVdsZ[1] = fVdsZ[0]+50;
    fVdsZ[2] = fVdsZ[0]+100;
    fVdsZ[3] = fVdsZ[0]+150.;

//fHomePath = gSystem->HomeDirectory();
    cout<<" Na61HotPixelsModule::Init() : home directory: "<<fHomePath.Data()<<endl;

}

//____________________________________________________________________
void Na61HotPixelsModule::Begin()
{
  // Run-level initialisation
  SetState(kBegin);
  
}

//____________________________________________________________________
void Na61HotPixelsModule::ReadNoisyPixelInfo()
{
  // If file with noisy pixel maps exist the noisy pixel analysis and recreatinion of this
  // file will be skept in Finish.
  fNoisyPixelsDefined = kFALSE;
  TFile* file=0;

  if(fJura){
    cout<<" Na61HotPixelsModule::ReadNoisyPixelInfo: opening file "
	<<Form("%s/pixeldb/alvd_jura_noisy_pixels_%d.root",fHomePath.Data(),fRunId)<<endl;
    file = new TFile(Form("%s/pixeldb/alvd_jura_noisy_pixels_%d.root",fHomePath.Data(),fRunId));
  }else{
    cout<<" Na61HotPixelsModule::ReadNoisyPixelInfo: opening file "
	<<Form("%s/pixeldb/alvd_saleve_noisy_pixels_%d.root",fHomePath.Data(),fRunId)<<endl;
    file = new TFile(Form("%s/pixeldb/alvd_saleve_noisy_pixels_%d.root",fHomePath.Data(),fRunId));
  }
  
  if(!file->IsOpen()){
    cout<<" Na61HotPixelModule::ReadNoisyPixelInfo: can not open noisy pixel file. Check it out!!!"<<endl;
    return;
  }

  // take histos with maps of noisy pixel

  fNoisyPixelsDefined = kTRUE; 
  for(Int_t imap=0;imap<34;imap++){ 
    TH2F* hh = (TH2F*)file->Get(Form("hNoisyPixelMap_%s",fAlSensorNames[imap].Data()));
    if(!hh){
      fNoisyPixelsDefined = kFALSE;   
    }else{ 
      Int_t ii=0;
      for(Int_t i=1;i<hh->GetNbinsX()+1;i++){
	for(Int_t j=1;j<hh->GetNbinsY()+1;j++){
	  if(hh->GetBinContent(i,j))ii++;
	  fhNoisyPixelMap[imap] -> SetBinContent(i,j,hh->GetBinContent(i,j));
	}
      }
      fhNoisyPixelMap[imap]->SetEntries(ii);

    }
  }  


  TH1F* h = (TH1F*)file->Get(Form("hCutValue"));
  Float_t cutValue = h->GetBinContent(1); 

  // require all maps to be defined, if not fNoisyPixelsDefined = false;
  file->Close();

  cout<<"cutValue="<<cutValue<<" "<<fPixelNoiseCut<<" "<<fNoisyPixelsDefined<<endl;
  if((cutValue != fPixelNoiseCut) && fNoisyPixelsDefined){
    cout<<"previous noisy pixel cut value was"<<cutValue<<" but current cut value is "<<fPixelNoiseCut<< endl;
    cout<<"Redo the noisy pixel maps? (1=yes, 0=no)"<<endl;
    Int_t iii;
    cin>>iii;
    if(iii==1)fNoisyPixelsDefined=kFALSE;
  } 

  if(fNoisyPixelsDefined){
    cout<<" Na61HotPixelsModule::Init: cut on noisy pixels activated with cut value = "
	<<cutValue<<" %"<<endl;
    for(Int_t i=0;i<34;i++)cout<<"Number of noisy pixels in sensor "<<fAlSensorNames[i].Data()
			      <<": "<<fhNoisyPixelMap[i]->GetEntries()<<endl;

  }else{
    cout<<" Na61HotPixelsModule::Init: will create new noisy pixel maps for cut value = "
	<<fPixelNoiseCut<<" %."<<" Cut on noisy pixels is disabled for this run."<<endl;
  }

}

//____________________________________________________________________
void Na61HotPixelsModule :: Event(UEventNode* inNode, UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);  

  UDataTable* tabs[34];
  UDataTable* tabs_out[34];

  //inNode->ListObjects();

  for(Int_t i=0; i<34; i++){
    if(fJura){
      tabs[i] = inNode->GetDataTable(Form("Jura Pixels %s merged",fAlSensorNames[i].Data()));
      tabs_out[i] = new UDataTable(Form("Jura Pixels %s cleaned",fAlSensorNames[i].Data()));
    }else{
      tabs[i] = inNode->GetDataTable(Form("Saleve Pixels %s merged",fAlSensorNames[i].Data()));
      tabs_out[i] = new UDataTable(Form("Saleve Pixels %s cleaned",fAlSensorNames[i].Data()));
    }
  }
  
  
  // fill pixel frequency distribution
  fEvents++;
  for(Int_t itab=0; itab<34; itab++){    
    if(!tabs[itab])continue;    
    for(Int_t i=0;i<tabs[itab]->GetEntries();i++){
      UChipPixel* pix = (UChipPixel*)tabs[itab]->At(i);    
      if(!pix)Error("Event","something is wrong, check it out");
      fhPixelMap[itab]->Fill(pix->GetLine(),pix->GetColumn());
    }
    
  }
  
  
  ResetFreqArray(); // needed for not writing double pixels
  
  AddPixels(tabs_out, tabs);


  for(Int_t i=0; i<34; i++){
    if( tabs[i])fhPixels[i] -> Fill(tabs[i]->GetEntries());
    if( tabs_out[i] ){
      fhCleanPixels[i] -> Fill(tabs_out[i]->GetEntries());
      outNode->AddDataTable(tabs_out[i]);
    }
  }
  


}


//_________________________________________________________________________________
void Na61HotPixelsModule::ResetFreqArray()
{
  // fmapArr[8]
  for(Int_t is=0;is<34;is++)
    for(Int_t i=0;i<findx+1;i++)
      for(Int_t j=0;j<findy+1;j++)
	faxayFreq[is][i][j]=kFALSE; 

}



//_________________________________________________________________________________
void Na61HotPixelsModule::AddPixels(UDataTable** tabs_out, UDataTable** tabs)
{

  for(Int_t itab=0; itab<34; itab++){
    
    if(!tabs[itab])continue;


    for(Int_t i=0; i<tabs[itab]->GetEntries(); i++){
      UChipPixel* pix = (UChipPixel*)tabs[itab]->At(i);

      // here put condition on noisy pixels
      Int_t ii = pix -> GetLine();
      Int_t jj = pix -> GetColumn();
 

      if(faxayFreq[itab][ii][jj]) continue; //to not store double pixels at the same position

      faxayFreq[itab][ii][jj]=kTRUE;


      Int_t content = 0;
 
      if(fNoisyPixelsDefined)content = fhNoisyPixelMap[itab]->GetBinContent(ii,jj);
 
	

      //note, that fNoisyPixelsDefined=kFALSE disables the cut

      if(!content)tabs_out[itab]->Add(new UChipPixel(ii,jj));

    }
    
  }
  
}

//_____________________________________________________________
void Na61HotPixelsModule::End()
{
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61HotPixelsModule::Finish()
{
  // Job-level finalisation
  SetState(kFinish);

   

  if(fNoisyPixelsDefined)return;

  //TString home_path = "/afs/cern.ch/user/p/pstaszel/public/branches/pstaszel/vdreco";

  TFile pixeldbFile;

 
  if(fJura){
    cout<<"Writing noise data to "<<Form("%s/pixeldb/alvd_jura_noisy_pixels_%d.root",fHomePath.Data(),fRunId)<<endl;
    pixeldbFile.Open(Form("%s/pixeldb/alvd_jura_noisy_pixels_%d.root",fHomePath.Data(),fRunId),"recreate");
  }else{
    cout<<"Writing noise data to "<<Form("%s/pixeldb/alvd_saleve_noisy_pixels_%d.root",fHomePath.Data(),fRunId)<<endl;
    pixeldbFile.Open(Form("%s/pixeldb/alvd_saleve_noisy_pixels_%d.root",fHomePath.Data(),fRunId),"recreate");
  }

  Int_t binsy  = findy;
  Int_t binsx  = findx;  
  Float_t ymax =  binsy;
  Float_t xmax =  binsx;

  TH1F* hcutValue = new TH1F("hCutValue","",1,0,1);
  hcutValue->SetBinContent(1,fPixelNoiseCut);
  hcutValue->Write();

  TH2F* hNoisyPixelMap[34]; 
  
  for(Int_t imap=0;imap<34;imap++){
    
    hNoisyPixelMap[imap] = new TH2F(Form("hNoisyPixelMap_%s",fAlSensorNames[imap].Data()),"", binsx, 0.5, xmax+0.5, binsy, 0.5, ymax+0.5);
    
    cout<<"Finish: imap="<<imap<<" "<<fhPixelMap[imap]<<endl;
    for(Int_t i=1;i<fhPixelMap[imap]->GetNbinsX()+1;i++){
      for(Int_t j=1;j<fhPixelMap[imap]->GetNbinsY()+1;j++){
	Float_t content =  fhPixelMap[imap]->GetBinContent(i,j);
	fhPixelFrequency[imap]->Fill(content*100./(Float_t)fEvents); // 0.5 because 2 frames added
	if((content*100./(Float_t)fEvents)>fPixelNoiseCut){
	  hNoisyPixelMap[imap]->SetBinContent(i,j,1);
	}
      }
    }
    
  }
  
  for(Int_t imap=0;imap<34;imap++)hNoisyPixelMap[imap]->Write();
  pixeldbFile.Close();

}



//____________________________________________________________________
void Na61HotPixelsModule ::Print(Option_t* option) const
{
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>
  
  TString opt(option);
  opt.ToLower(); 
  
  Na61Module::Print(option); 
  if (opt.Contains("d")) 
    cout << endl 
         << "  Original author: Paweł Staszel" << endl
         << "  Last Modifications: " << endl 
         << "    $Author: Staszel $" << endl  
         << "    $Date: 2022/05/09$"   << endl 
         << "    $Revision: 1.0 $ " << endl  
         << endl 
         << "-------------------------------------------------" << endl;
}



