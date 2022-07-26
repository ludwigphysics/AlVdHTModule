//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61VdEfficiencyModule
#include "Na61VdEfficiencyModule.h"    
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
//ClassImp(Na61VdEfficiencyModule);

//____________________________________________________________________
Na61VdEfficiencyModule::Na61VdEfficiencyModule()
{
  // Default constructor. DO NOT USE
  SetState(kSetup);
}

//____________________________________________________________________
Na61VdEfficiencyModule::Na61VdEfficiencyModule(const Char_t* name, const Char_t* title)
   : Na61Module(name, title)
{
  // Named Constructor
  SetState(kSetup);
  
 
}
//____________________________________________________________________
void Na61VdEfficiencyModule::DefineHistograms()
{

if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init"); 
    return;  
    
  }
  
 TDirectory* histDir;
 
 if(fJura)histDir = gDirectory->mkdir("VdEfficiencyModule_Jura");
 else histDir = gDirectory->mkdir("VdEfficiencyModule_Saleve");
 
 histDir->cd(); 
 
 // volume in global VD coordinate system
 fhVolumesX = new TH1F("VolumesX","",8,0,8);
 fhVolumesY = new TH1F("VolumesY","",8,0,8);
 fhVolumesZ = new TH1F("VolumesZ","",8,0,8);

 Float_t xmin = -25;
 Float_t xmax =  25;
 
 for(Int_t i=0;i<8;i++){
 
   fhRefTracks[i] =   new TH1F(Form("hRefTracks_%s",fSensorNames[i].Data())," ", 150,0.,150.);
   
   fhDx[i] = new TH1F(Form("hDx_%s",fSensorNames[i].Data())," ", 5000,-5.,5.);
   fhDy[i] = new TH1F(Form("hDy_%s",fSensorNames[i].Data())," ", 5000,-5.,5.);
   
   fhDx_cuty[i] = new TH1F(Form("hDx_%s_cuty",fSensorNames[i].Data())," ", 5000,-5.,5.);
   fhDy_cutx[i] = new TH1F(Form("hDy_%s_cutx",fSensorNames[i].Data())," ", 5000,-5.,5.);
   
   fhxy[i] = new TH2F(Form("hxy_%s",fSensorNames[i].Data())," ", 500, xmin, xmax,  500, -25.0, 25.0);
   fhxy_match[i] = new TH2F(Form("hxy_%s_match",fSensorNames[i].Data())," ", 500, xmin, xmax,  500, -25.0, 25.0);

   fhxy_loc[i] = new TH2F(Form("hxy_loc_%s",fSensorNames[i].Data())," ", 200, -5.4, 5.4, 400, -10.7, 10.7);
   fhxy_loc_match[i] = new TH2F(Form("hxy_loc_%s_match",fSensorNames[i].Data())," ", 200, -5.4, 5.4,  400, -10.7, 10.7);
 }
 
 
 fhHitsXVsZ = new TH2F("HitsXVsZ","",5000,-1,151, 1000,-20,20.);   
 
 // Define histograms. They are:
 // <fill in here>
 
 
 gDirectory->cd("..");  
}

//____________________________________________________________________
void Na61VdEfficiencyModule::Init()
{
  // Job-level initialisation
  SetState(kInit);
  

  if(fJura)fParams = (Na61VdParametersManager::Instance())->GetJuraArmParams();  
  else fParams = (Na61VdParametersManager::Instance())->GetSaleveArmParams(); 

  cout<<"Na61VdEfficiencyModule::Init: fParams="<<(Na61ArmParameters*)fParams<<"  "<<fParams->GetInit()<<endl;

  if(!fParams->GetInit()){
    fParams -> SetDrotX(0,0); 
    fParams -> SetDrotY(0,0); 
    fParams -> SetDrotZ(0,0); 
    fParams -> SetVolumeDz(0,0);
    
    fParams -> SetRunId(fRunId);
    fParams->Init();
  }

  fVdParams = (Na61VdParametersManager::Instance())->GetVdParams();  
  fVdParams -> SetRunId(fRunId);
  if(!fVdParams->GetInit())fVdParams->Init();

  if(fJura)SetupMatchingOffsetsAndSigmas_Jura(fRunId);
  else SetupMatchingOffsetsAndSigmas_Saleve(fRunId);

  fNsig = 5;

  if(fJura){
    fArmOffsetX = fVdParams->GetJuraArmOffset().X();
    fArmOffsetY = fVdParams->GetJuraArmOffset().Y();
    fArmOffsetZ = fVdParams->GetJuraArmOffset().Z();
  }else{
    fArmOffsetX = fVdParams->GetSaleveArmOffset().X();
    fArmOffsetY = fVdParams->GetSaleveArmOffset().Y();
    fArmOffsetZ = fVdParams->GetSaleveArmOffset().Z();
  }

}

//____________________________________________________________________
void Na61VdEfficiencyModule::Begin()
{
  // Run-level initialisation
  SetState(kBegin);
  
}

//____________________________________________________________________
void Na61VdEfficiencyModule :: Event(UEventNode* inNode, UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);  

  //outNode->ListObjects();
  //UDataTable* hits = inNode->GetDataTable(Form("Hits %s",fSensorNames[0].Data()));
  //if(hits->GetEntries() > fMultiCut) return;

  if(fVdAna){
    // normal efficiency calculations unis vdana program
    for(Int_t i=0;i<8;i++){
      //FindEfficiencyOfSensor(i,outNode); // if vdreco
      FindEfficiencyOfSensor(i,inNode); // if vdana
    }
    cout<<"VdAna=true"<<endl;
  } else {
    // this is used rather for the geometry tuning using vdreco program
    FindEfficiencyOfSensor_3Hit(4,inNode,outNode);
    FindEfficiencyOfSensor_3Hit(5,inNode,outNode);
    FindEfficiencyOfSensor_3Hit(6,inNode,outNode);
    FindEfficiencyOfSensor_3Hit(7,inNode,outNode);
        cout<<"VdAna=false"<<endl;
  }


}

//____________________________________________________________________________________
void  Na61VdEfficiencyModule :: FindEfficiencyOfSensor(Int_t ii, UEventNode* inNode)
{


  // select reference tracks
  Int_t itab[4];
  Int_t Ntab=0;
  if(ii==0){itab[0]=4; itab[1]=8; itab[2]=12; itab[3]=15; Ntab=4;};
  if(ii==1){itab[0]=5; itab[1]=9; itab[2]=13; itab[3]=16; Ntab=4;};
  if(ii==2){itab[0]=6;  itab[1]=14; Ntab=2;};
  if(ii==3){itab[0]=10; itab[1]=17;Ntab=2;};
  if(ii==4){itab[0]=7;  Ntab=1;};
  if(ii==5){itab[0]=11; Ntab=1;};
  if(ii==6){itab[0]=7;  Ntab=1;};
  if(ii==7){itab[0]=11; Ntab=1;};

  Int_t ref_tracks=0;
  TObjArray tracks;
  tracks.Clear();
  for(Int_t it=0;it<Ntab;it++){
    UDataTable* tracktab = inNode -> GetDataTable(Form("FullTracks %s",fMatchStr[itab[it]].Data()));
      if(!tracktab)continue;
      
      for(Int_t i=0; i<tracktab->GetEntries(); i++){
		UVdTrack* track = (UVdTrack*)tracktab->At(i);
		//track -> Activate();

		tracks.Add(track); 
		ref_tracks++;
      }
  }


  fhRefTracks[ii]->Fill(ref_tracks);

  if(!ref_tracks)return;

  UDataTable* hits = inNode -> GetDataTable(Form("Hits %s",fSensorNames[ii].Data()));


  Float_t Offx = fOffx[ii];  Float_t Sigx = fSigx[ii];  
  Float_t Offy = fOffy[ii];  Float_t Sigy = fSigy[ii];  
  //Float_t Offx = fOffx[4];  Float_t Sigx = fSigx[4];  
  //Float_t Offy = fOffy[4];  Float_t Sigy = fSigy[4];  

  // this can be also taken from hitos in vdrecofiles
  Float_t VolumeX = fParams->GetVolumeGlobX(ii);//+fArmOffsetX; //volumes are transform to global in PrimaryVertexReco
  Float_t VolumeY = fParams->GetVolumeGlobY(ii);//+fArmOffsetY;
  Float_t VolumeZ = fParams->GetVolumeGlobZ(ii);//+fArmOffsetZ;
//  cout<<"VolumeX "<<VolumeX<<" VolumeY "<<VolumeY<<" VolumeZ "<<VolumeZ<<endl;

  Float_t tan_beta = TMath::Tan(fParams->GetRotGlobY(ii));

  //Float_t ze = VolumeZ;
  //cout<<"ii="<<ii<<" ze="<<ze<<"  tan_beta="<<tan_beta<<"   roty="<<fParams->GetRotY(ii)<<endl;
  //cout<<"ii="<<ii<<" "<<fArmOffsetX<<" "<<fArmOffsetY<<" "<<fArmOffsetZ<<endl;
  //cout<<"ii="<<ii<<"  "<<fParams->GetVolumeX(ii)<<" "<<fParams->GetVolumeY(ii)<<" "<<fParams->GetVolumeZ(ii)<<endl;

  for(Int_t i=0;i<tracks.GetEntries();i++){

    UVdTrack* vdtrack = (UVdTrack*)tracks.At(i);

    //Float_t xe =  vdtrack->GetXatZ(VolumeZ);
    Float_t xe =  vdtrack->GetXatZ_pol2(VolumeZ);
    Float_t ye =  vdtrack->GetYatZ(VolumeZ);
    
    Float_t x1 = xe - VolumeX;

    Float_t xslope = vdtrack->GetDX()/vdtrack->GetDZ();

    Float_t dz = - x1*(tan_beta/(1+tan_beta*xslope));

    //cout<<"ii="<<ii<<" x1="<<x1<<"  tan_beta="<<tan_beta<<"  dz="<<dz<<" ze+dz ="<< ze+dz <<endl;

    Float_t ze = VolumeZ + dz;
    //xe =  vdtrack->GetXatZ(ze);
    xe =  vdtrack->GetXatZ_pol2(ze);
    ye =  vdtrack->GetYatZ(ze);

    fhHitsXVsZ->Fill(ze, xe);   
    
    // fine track position in local sensor frame
    Float_t xloc = xe-VolumeX;
    Float_t yloc = ye-VolumeY;
    Float_t zloc = ze-VolumeZ;
    Float_t xlocsf = TMath::Sqrt(xloc*xloc + zloc*zloc);
    Float_t ylocsf = TMath::Sqrt(yloc*yloc + zloc*zloc);

    //cout<<"xloc "<<xloc<<" yloc "<<yloc<<" zloc "<<zloc<<endl; 
    
    
    if(xloc==0)xlocsf=0;
    else xlocsf = (xloc/TMath::Abs(xloc)) * xlocsf; 

    if(yloc==0)ylocsf=0;
    else ylocsf = (yloc/TMath::Abs(yloc)) * ylocsf; 
    
    //cout<<"xlocsf "<<xlocsf<<" ylocsf "<<ylocsf<<endl; 


    Float_t ca = TMath::Cos(fParams->GetRotZ(ii));
    Float_t sa = TMath::Sin(fParams->GetRotZ(ii));
    
    Float_t xlocf =  ca*xlocsf - sa*ylocsf;
    Float_t ylocf =  sa*xlocsf + ca*ylocsf;
 
//    cout<<"xlocf "<<xlocf<<" ylocf "<<ylocf<<endl; 
    
    fhxy[ii]->Fill(xe,ye);
    if((xlocf > -5.4) && (xlocf < 5.4)) fhxy_loc[ii]->Fill(xlocf,ylocf);

    Int_t multiMatch=0;
    for(Int_t j=0;j<hits->GetEntries();j++){
      USensorHit* hit = (USensorHit*)hits->At(j);

      Float_t x =  hit->GetX();
      Float_t y =  hit->GetY();

      Float_t dx = xe - x;
      Float_t dy = ye - y;
    
      fhDx[ii]->Fill(dx);
      fhDy[ii]->Fill(dy);

      if(TMath::Abs(dy-Offy) < 3*Sigy) fhDx_cuty[ii]->Fill(dx);
      if(TMath::Abs(dx-Offx) < 3*Sigx) fhDy_cutx[ii]->Fill(dy);

      Float_t dev = (dx-Offx)*(dx-Offx)/(Sigx*Sigx) + (dy-Offy)*(dy-Offy)/(Sigy*Sigy);

      if(dev > fNsig*fNsig)continue;
      // got matching

      if(multiMatch)break;

      fhxy_match[ii]->Fill(xe,ye);
      if((xlocf > -5.4) && (xlocf < 5.4)) fhxy_loc_match[ii]->Fill(xlocf,ylocf);
      multiMatch++;
      
    }
  }
  

}
//____________________________________________________________________________________
void  Na61VdEfficiencyModule :: FindEfficiencyOfSensor_3Hit(Int_t ii, UEventNode* /* inNode */,UEventNode* outNode)
{
  // select reference tracks
 
  Int_t Ntab=2;
  TString devName[2]={"dev1_0","dev1_1"};

  Int_t ref_tracks=0;
  TObjArray tracks;
  tracks.Clear();
  for(Int_t it=0;it<Ntab;it++){
    UDataTable* tracktab = outNode -> GetDataTable(Form("Tracks %s",devName[it].Data()));
      if(!tracktab)continue;
      
      for(Int_t i=0; i<tracktab->GetEntries(); i++){
	UVdTrack* track = (UVdTrack*)tracktab->At(i);
	//track -> Activate();

	tracks.Add(track); 
	ref_tracks++;
      }
  }



  fhRefTracks[ii]->Fill(ref_tracks);

  if(!ref_tracks)return;

  UDataTable* hits = outNode -> GetDataTable(Form("Hits %s",fSensorNames[ii].Data()));

  //cout<<"ii="<<ii<<" ref_tracks="<<ref_tracks<<" hits:"<<hits<<endl;

  Float_t Offx = fOffx[ii];  Float_t Sigx = fSigx[ii];  
  Float_t Offy = fOffy[ii];  Float_t Sigy = fSigy[ii];  
  //Float_t Offx = fOffx[4];  Float_t Sigx = fSigx[4];  
  //Float_t Offy = fOffy[4];  Float_t Sigy = fSigy[4];  

  // this can be also taken from hitos in vdrecofiles
  Float_t VolumeX = fParams->GetVolumeX(ii);//+fArmOffsetX;
  Float_t VolumeY = fParams->GetVolumeY(ii);//+fArmOffsetY;
  Float_t VolumeZ = fParams->GetVolumeZ(ii);//+fArmOffsetZ;

  Float_t tan_beta = TMath::Tan(fParams->GetRotY(ii));

  for(Int_t i=0;i<tracks.GetEntries();i++){

    UVdTrack* vdtrack = (UVdTrack*)tracks.At(i);

    Float_t xe =  vdtrack->GetXatZ(VolumeZ);
    Float_t ye =  vdtrack->GetYatZ(VolumeZ);
    
    Float_t x1 = xe - VolumeX;

    Float_t xslope = vdtrack->GetDX()/vdtrack->GetDZ();

    Float_t dz = - x1*(tan_beta/(1+tan_beta*xslope));

    //cout<<"ii="<<ii<<" x1="<<x1<<"  tan_beta="<<tan_beta<<"  dz="<<dz<<" ze+dz ="<< ze+dz <<endl;

    Float_t ze = VolumeZ + dz;
    xe =  vdtrack->GetXatZ(ze);
    ye =  vdtrack->GetYatZ(ze);
   

    fhHitsXVsZ->Fill(ze, xe);   
    
    // fine track position in local sensor frame
    Float_t xloc = xe-VolumeX;
    Float_t yloc = ye-VolumeY;
    Float_t zloc = ze-VolumeZ;
    Float_t xlocf = TMath::Sqrt(xloc*xloc + zloc*zloc);
    Float_t ylocf = TMath::Sqrt(yloc*yloc + zloc*zloc);

    if(xloc==0)xlocf=0;
    else xlocf = (xloc/TMath::Abs(xloc)) * xlocf; 

    if(yloc==0)ylocf=0;
    else ylocf = (yloc/TMath::Abs(yloc)) * ylocf; 


    //if(xloc<5.3)cout<<" ii="<<ii<<" xloc="<<xloc<<"  xlocf="<<xlocf<<"  yloc="<<yloc<<endl; 

    if(!((yloc>-10.6) && (yloc<10.6)))continue;
    if(!((xlocf>-5.3) && (xlocf<5.3)))continue;
    //cout<<"ii="<<ii<<"  xe="<<xe<<"  volumeX="<<VolumeX<<"  xloc="<<xloc<<endl;

    fhxy[ii]->Fill(xe,ye);
    //if((xloc > -5.4) && (xloc < 5.4)) 
    fhxy_loc[ii]->Fill(xlocf,yloc);


    Int_t multiMatch=0;
    for(Int_t j=0;j<hits->GetEntries();j++){
      USensorHit* hit = (USensorHit*)hits->At(j);

      Float_t x =  hit->GetX();
      Float_t y =  hit->GetY();

      Float_t dx = xe - x;
      Float_t dy = ye - y;

      fhDx[ii]->Fill(dx);
      fhDy[ii]->Fill(dy);

      if(TMath::Abs(dy-Offy) < 3*Sigy) fhDx_cuty[ii]->Fill(dx);
      if(TMath::Abs(dx-Offx) < 3*Sigx) fhDy_cutx[ii]->Fill(dy);

      Float_t dev = (dx-Offx)*(dx-Offx)/(Sigx*Sigx) + (dy-Offy)*(dy-Offy)/(Sigy*Sigy);

      if(dev > fNsig*fNsig)continue;
      // got matching

      if(multiMatch)break;

      fhxy_match[ii]->Fill(xe,ye);
      if((xloc > -5.4) && (xloc < 5.4)) fhxy_loc_match[ii]->Fill(xloc,yloc);
      multiMatch++;
      
    }
  }
  
}


//_____________________________________________________________
void Na61VdEfficiencyModule::SetupMatchingOffsetsAndSigmas_Jura(Int_t run_id)
{

  // dev cuts
  // find right offsets and sigmas for 8 sensors.
  // for present geometry (march 2017)  
  fMultiCut=50+50;

  // angle cuts
  fOffx[0] =  0.0052;  fSigx[0] =  0.0124; 
  fOffy[0] = -0.002;   fSigy[0] =  0.0131;

  fOffx[1] = -0.003;   fSigx[1] =  0.00717; 
  fOffy[1] =  0.00068; fSigy[1] =  0.0079;

  fOffx[2] = -0.00052; fSigx[2] =  0.00722; 
  fOffy[2] =  0.00188; fSigy[2] =  0.006787;

  fOffx[3] =  0.000377; fSigx[3] = 0.00745; 
  fOffy[3] =  0.00189;  fSigy[3] = 0.007894;

  fOffx[4] =  0.003221; fSigx[4] =  0.01449; 
  fOffy[4] = -0.006;    fSigy[4] =  0.01331;


  fOffx[5] =  0.003854; fSigx[5] =  0.01271; 
  fOffy[5] = -0.00132;  fSigy[5] =  0.01239;

  fOffx[6] =   0.0118;      fSigx[6] = 0.03; 
  fOffy[6] =   -0.005;      fSigy[6] = 0.02;

  fOffx[7] =  0.0;        fSigx[7] =  0.012; 
  fOffy[7] =  0.0;        fSigy[7] =  0.013;

  if((run_id>300) && (run_id<601)){
    fOffx[0] = 0.00;    fSigx[0] = 0.01704;
    fOffy[0] = 0.00;    fSigy[0] = 0.01303;
    
    fOffx[1] = 0.000;   fSigx[1] = 0.01099;
    fOffy[1] = 0.000;   fSigy[1] = 0.00799;
    
    fOffx[2] = 0.0;     fSigx[2] = 0.00873;
    fOffy[2] = 0.00;    fSigy[2] = 0.00818;
    
    fOffx[3] = 0.00;    fSigx[3] = 0.01082;
    fOffy[3] = 0.00;    fSigy[3] = 0.00891;
    
    fOffx[4] = 0.0;     fSigx[4] = 0.01;
    fOffy[4] = 0.0;     fSigy[4] = 0.01;
    
    fOffx[5] = 0.0;     fSigx[5] = 0.01;
    fOffy[5] = 0.0;     fSigy[5] = 0.01;
    
    fOffx[6] = 0.0;     fSigx[6] = 0.019;
    fOffy[6] = 0.0;     fSigy[6] = 0.019;
    
    fOffx[7] = 0.0;     fSigx[7] = 0.019;
    fOffy[7] = 0.0;     fSigy[7] = 0.019;
  }

  if((run_id>600) && (run_id<1493)){
    cout<<"Na61VdEfficiencyModule::SetupMatchingOffsetsAndSigmas_Jura: run_id="<<run_id<<endl;

    fOffx[0] = -0.00174;   fSigx[0] = 0.02129;
    fOffy[0] = -0.00033;   fSigy[0] = 0.00978;
    
    fOffx[1] = -0.00103;   fSigx[1] = 0.01887;
    fOffy[1] = 0.00046;   fSigy[1] = 0.00653;
    
    fOffx[2] = -0.00314;   fSigx[2] = 0.00748;
    fOffy[2] = -0.00002;   fSigy[2] = 0.00762;
    
    fOffx[3] = 0.00239;   fSigx[3] = 0.00802;
    fOffy[3] = -0.00078;   fSigy[3] = 0.00730;
    
    fOffx[4] = 0.00842;   fSigx[4] = 0.05148;
    fOffy[4] = -0.00825;   fSigy[4] = 0.01385;
    
    fOffx[5] = 0.02091;   fSigx[5] = 0.05232;
    fOffy[5] = 0.00198;   fSigy[5] = 0.01515;
    
    fOffx[6] = -0.00118;   fSigx[6] = 0.05304;
    fOffy[6] = 0.00283;   fSigy[6] = 0.01623;
    
    fOffx[7] = 0.00677;   fSigx[7] = 0.05506;
    fOffy[7] = 0.00187;   fSigy[7] = 0.01610;
  }

  if((run_id>31974) && (run_id<33800)){
    cout<<"Na61VdEfficiencyModule::SetupMatchingOffsetsAndSigmas_Jura: run_id="<<run_id<<endl;

    fOffx[0] =  -1.00434e-03;   fSigx[0] = 3.22930e-02;
    fOffy[0] = -2.46582e-03;   fSigy[0] = 1.70687e-02;
    
    fOffx[1] =  3.45964e-04;   fSigx[1] = 1.01510e-02;
    fOffy[1] = 1.04049e-03;   fSigy[1] = 9.95691e-03;
    
    fOffx[2] = 2.26483e-03;   fSigx[2] = 9.77605e-03;
    fOffy[2] = -8.16947e-04;   fSigy[2] = 9.78091e-03;
    
    fOffx[3] = -1.41427e-03;   fSigx[3] = 1.01895e-02;
    fOffy[3] = -2.60608e-03;   fSigy[3] = 1.02570e-02;
    
    fOffx[4] = -4.50395e-03;   fSigx[4] = 2.82322e-02;
    fOffy[4] = -3.19345e-03;   fSigy[4] = 2.16495e-02;
    
    fOffx[5] = 6.93762e-03;   fSigx[5] = 2.96477e-02;
    fOffy[5] = -3.50125e-03;   fSigy[5] = 2.05241e-02;
    
    fOffx[6] =  -2.22877e-03;   fSigx[6] = 2.66076e-02;
    fOffy[6] =  -3.60754e-03;   fSigy[6] = 2.28926e-02;
    
    fOffx[7] = 9.36453e-03;   fSigx[7] =  2.51524e-02;
    fOffy[7] = 2.68193e-04;   fSigy[7] = 2.27726e-02;

/*    fOffx[0] = 0.00197;   fSigx[0] = 0.01146;
    fOffy[0] = -0.00094;   fSigy[0] = 0.01079;
    
    fOffx[1] = -0.00172;   fSigx[1] = 0.00651;
    fOffy[1] = 0.00081;   fSigy[1] = 0.00615;
    
    fOffx[2] = -0.00339;   fSigx[2] = 0.00741;
    fOffy[2] = 0.00155;   fSigy[2] = 0.00728;
    
    fOffx[3] = 0.00130;   fSigx[3] = 0.00777;
    fOffy[3] = -0.00025;   fSigy[3] = 0.00691;
    
    fOffx[4] = 0.00708;   fSigx[4] = 0.01584;
    fOffy[4] = -0.00845;   fSigy[4] = 0.01453;
    
    fOffx[5] = 0.01605;   fSigx[5] = 0.01491;
    fOffy[5] = 0.00032;   fSigy[5] = 0.01474;
    
    fOffx[6] = -0.00079;   fSigx[6] = 0.01510;
    fOffy[6] = 0.00216;   fSigy[6] = 0.01617;
    
    fOffx[7] = 0.00633;   fSigx[7] = 0.01526;
    fOffy[7] = -0.00010;   fSigy[7] = 0.01623;
*/
  }

  if((run_id>37797)){
    cout<<"Na61VdEfficiencyModule::SetupMatchingOffsetsAndSigmas_Jura: run_id="<<run_id<<" pbpb 2018"<<endl;
    fOffx[0] = 0.02186;   fSigx[0] = 0.01851;
    fOffy[0] = -0.00191;   fSigy[0] = 0.00850;
    
    fOffx[1] = -0.00765;   fSigx[1] = 0.00697;
    fOffy[1] = 0.00043;   fSigy[1] = 0.00559;
    
    fOffx[2] = 0.00717;   fSigx[2] = 0.00639;
    fOffy[2] = -0.00416;   fSigy[2] = 0.00564;
    
    fOffx[3] = 0.00866;   fSigx[3] = 0.00612;
    fOffy[3] = -0.00170;   fSigy[3] = 0.00592;
    
    fOffx[4] = -0.02288;   fSigx[4] = 0.01918;
    fOffy[4] = -0.00178;   fSigy[4] = 0.01109;
    
    fOffx[5] = -0.02348;   fSigx[5] = 0.01913;
    fOffy[5] = -0.00005;   fSigy[5] = 0.01382;
    
    fOffx[6] = -0.01441;   fSigx[6] = 0.01697;
    fOffy[6] = 0.00480;   fSigy[6] = 0.01323;
    
    fOffx[7] = -0.01214;   fSigx[7] = 0.01854;
    fOffy[7] = -0.00350;   fSigy[7] = 0.01321;  
  }
  
}

//_____________________________________________________________
void Na61VdEfficiencyModule::SetupMatchingOffsetsAndSigmas_Saleve(Int_t run_id)
{

  // dev cuts
  // find right offsets and sigmas for 8 sensors.
  // for present geometry (march 2017)  

  fMultiCut=70+50;

  // angle cuts
  fOffx[0] = -0.00017;   fSigx[0] = 0.017; 
  fOffy[0] = -0.00497;   fSigy[0] = 0.01393;

  fOffx[1] =  0.000366;  fSigx[1] = 0.00916; 
  fOffy[1] = -0.002194;  fSigy[1] = 0.0079;

  fOffx[2] = -0.002739;  fSigx[2] = 0.01048; 
  fOffy[2] = -0.002787;  fSigy[2] = 0.008319;

  fOffx[3] =  0.001;     fSigx[3] = 0.007613; 
  fOffy[3] = -0.002625;  fSigy[3] = 0.008127;

  fOffx[4] = -0.00281;   fSigx[4] = 0.01852; 
  fOffy[4] =  0.00727;   fSigy[4] = 0.01472;


  fOffx[5] =  0.008565;  fSigx[5] = 0.01354; 
  fOffy[5] =  0.00518;   fSigy[5] = 0.01396;

  fOffx[6] = -0.000399;  fSigx[6] = 0.01537; 
  fOffy[6] =  0.00606;   fSigy[6] = 0.01422;

  fOffx[7] =  0.00463;   fSigx[7] = 0.01541; 
  fOffy[7] =  0.00417;   fSigy[7] = 0.01414;

  if((run_id>300) && (run_id<601)){
    fOffx[0] = 0.0;   fSigx[0] = 0.01311;
    fOffy[0] = 0.0;   fSigy[0] = 0.01137;
    
    fOffx[1] = 0.0;   fSigx[1] = 0.01197;
    fOffy[1] = 0.0;   fSigy[1] = 0.01137;
    
    fOffx[2] = 0.0;   fSigx[2] = 0.01793;
    fOffy[2] = 0.0;   fSigy[2] = 0.01279;
    
    fOffx[3] = 0.0;   fSigx[3] = 0.01394;
    fOffy[3] = 0.0;   fSigy[3] = 0.01375;
    
    fOffx[4] = 0.0;   fSigx[4] = 0.01;
    fOffy[4] = 0.0;   fSigy[4] = 0.01;
    
    fOffx[5] = 0.0;   fSigx[5] = 0.01;
    fOffy[5] = 0.0;   fSigy[5] = 0.01;
    
    fOffx[6] = 0.0;   fSigx[6] = 0.03;
    fOffy[6] = 0.0;   fSigy[6] = 0.02;
    
    fOffx[7] = 0.0;   fSigx[7] = 0.03;
    fOffy[7] = 0.0;   fSigy[7] = 0.02;
  }

  if((run_id>600) && (run_id<1493)){
    cout<<"Na61VdEfficiencyModule::SetupMatchingOffsetsAndSigmas_Saleve: run_id="<<run_id<<endl;

    fOffx[0] = -0.00231;   fSigx[0] = 0.02522;
    fOffy[0] = 0.00306;    fSigy[0] = 0.01355;
    
    fOffx[1] = -0.00139;   fSigx[1] = 0.01929;
    fOffy[1] = -0.00114;   fSigy[1] = 0.00922;
    
    fOffx[2] = 0.00652;    fSigx[2] = 0.01421;
    fOffy[2] = 0.00093;    fSigy[2] = 0.01291;
    
    fOffx[3] = 0.00060;    fSigx[3] = 0.01065;
    fOffy[3] = -0.00219;   fSigy[3] = 0.01180;
    
    fOffx[4] = -0.01467;   fSigx[4] = 0.06302;
    fOffy[4] = 0.01054;    fSigy[4] = 0.02051;
    
    fOffx[5] = -0.00085;   fSigx[5] = 0.06368;
    fOffy[5] = 0.00783;    fSigy[5] = 0.01965;
    
    fOffx[6] = 0.00907;    fSigx[6] = 0.06044;
    fOffy[6] = 0.00487;    fSigy[6] = 0.01773;
    
    fOffx[7] = 0.00534;    fSigx[7] = 0.05944;
    fOffy[7] = -0.00278;   fSigy[7] = 0.01982;
  }

  if((run_id>31974) && (run_id<33800)){
    cout<<"Na61VdEfficiencyModule::SetupMatchingOffsetsAndSigmas_Saleve: run_id="<<run_id<<endl;
    fOffx[0] =  -2.58020e-03;   fSigx[0] = 3.04473e-02;
    fOffy[0] =  1.20454e-04;   fSigy[0] = 1.61819e-02;
    
    fOffx[1] = 8.80085e-04;   fSigx[1] = 9.76542e-03;
    fOffy[1] = -1.32837e-04;   fSigy[1] = 9.62393e-03;
    
    fOffx[2] = -5.71379e-03;   fSigx[2] = 1.11056e-02;
    fOffy[2] = -1.49905e-03;   fSigy[2] = 9.84530e-03;
    
    fOffx[3] = -5.76055e-03;   fSigx[3] = 1.08790e-02;
    fOffy[3] = 1.85774e-03;   fSigy[3] =  9.75738e-03;
    
    fOffx[4] = 1.15955e-03;   fSigx[4] = 2.80261e-02;
    fOffy[4] = 2.92785e-04;   fSigy[4] = 2.16906e-02;
    
    fOffx[5] = 4.54937e-04;   fSigx[5] = 2.82621e-02;
    fOffy[5] = -3.40307e-04;   fSigy[5] = 2.18637e-02;
    
    fOffx[6] = 2.91675e-03;   fSigx[6] = 2.53690e-02;
    fOffy[6] = 1.58622e-04;   fSigy[6] = 2.05032e-02;
    
    fOffx[7] = 3.94328e-03;   fSigx[7] = 2.46529e-02;
    fOffy[7] = -6.43800e-04;   fSigy[7] = 2.19019e-02;  
      
/*    fOffx[0] = 0.00182;   fSigx[0] = 0.01429;
    fOffy[0] = 0.00560;   fSigy[0] = 0.01556;
    
    fOffx[1] = -0.00396;   fSigx[1] = 0.00782;
    fOffy[1] = -0.00288;   fSigy[1] = 0.00983;
    
    fOffx[2] = 0.00981;   fSigx[2] = 0.01377;
    fOffy[2] = 0.00355;   fSigy[2] = 0.01423;
    
    fOffx[3] = 0.00336;   fSigx[3] = 0.00978;
    fOffy[3] = -0.00221;   fSigy[3] = 0.01189;
    
    fOffx[4] = -0.01154;   fSigx[4] = 0.02915;
    fOffy[4] = 0.01975;   fSigy[4] = 0.02145;
    
    fOffx[5] = 0.00708;   fSigx[5] = 0.02922;
    fOffy[5] = 0.01454;   fSigy[5] = 0.02002;
    
    fOffx[6] = 0.01295;   fSigx[6] = 0.01699;
    fOffy[6] = 0.00552;   fSigy[6] = 0.01983;
    
    fOffx[7] = 0.00580;   fSigx[7] = 0.01729;
    fOffy[7] = -0.00050;   fSigy[7] = 0.02032;    
*/
  }

  //if((run_id>2272)){
  if((run_id>37797)){
    cout<<"Na61VdEfficiencyModule::SetupMatchingOffsetsAndSigmas_Saleve: run_id="<<run_id<<" pbpb 2018"<<endl;
    fOffx[0] = 0.03074;   fSigx[0] = 0.02108;
    fOffy[0] = -0.00545;   fSigy[0] = 0.00912;
    
    fOffx[1] = -0.01110;   fSigx[1] = 0.00707;
    fOffy[1] = 0.00199;   fSigy[1] = 0.00595;
    
    fOffx[2] = 0.00589;   fSigx[2] = 0.00791;
    fOffy[2] = -0.00108;   fSigy[2] = 0.00674;
    
    fOffx[3] = 0.00511;   fSigx[3] = 0.00789;
    fOffy[3] = 0.00088;   fSigy[3] = 0.00643;
    
    fOffx[4] = -0.02829;   fSigx[4] = 0.02044;
    fOffy[4] = -0.00819;   fSigy[4] = 0.01329;
    
    fOffx[5] = -0.03082;   fSigx[5] = 0.02079;
    fOffy[5] = -0.01016;   fSigy[5] = 0.01337;
    
    fOffx[6] = -0.03695;   fSigx[6] = 0.01855;
    fOffy[6] = -0.00396;   fSigy[6] = 0.01382;
    
    fOffx[7] = -0.03637;   fSigx[7] = 0.01836;
    fOffy[7] = -0.00829;   fSigy[7] = 0.01344;
  }


}


//_____________________________________________________________
void Na61VdEfficiencyModule::End()
{
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61VdEfficiencyModule::Finish()
{
  // Job-level finalisation
  SetState(kFinish);
  
  for(Int_t i=0;i<8;i++){
    Float_t volx = fParams->GetVolumeGlobX(i);//+fArmOffsetX;
    Float_t voly = fParams->GetVolumeGlobY(i);//+fArmOffsetY;
    Float_t volz = fParams->GetVolumeGlobZ(i);//+fArmOffsetZ;
    fhVolumesX->SetBinContent(i+1,volx);
    fhVolumesY->SetBinContent(i+1,voly);
    fhVolumesZ->SetBinContent(i+1,volz);
  }

  //////////// plote average efficiencies
  TString JuraAddStr[8]   = {"d070_1","d070_0","d071_0","d071_1","d073_0","d073_1","d072_0","d072_1"};
  TString SaleveAddStr[8] = {"d450_1","d450_0","d451_1","d451_0","d453_1","d453_0","d452_1","d452_0"};

  cout<<endl;

  Float_t x1 =  -5.3;
  Float_t x2 =   5.3;
  Float_t y1 =  -10.6;
  Float_t y2 =   10.6;

  Float_t deltaX = 0.1;
  Float_t deltaY = 0.1;

  if(fJura)cout<<"SENSOR EFFCIENCIES in JURA ARM:"<<endl;
  else cout<<"SENSOR EFFCIENCIES in SALEVE ARM:"<<endl;

  for(Int_t i=0;i<8;i++){  
    
    Int_t ibx_min = fhxy_loc[i]->GetXaxis()->FindBin(x1+deltaX);
    Int_t ibx_max = fhxy_loc[i]->GetXaxis()->FindBin(x2-deltaX);
    Int_t iby_min = fhxy_loc[i]->GetYaxis()->FindBin(y1+deltaY);
    Int_t iby_max = fhxy_loc[i]->GetYaxis()->FindBin(y2-deltaY);
    
    Float_t ref = fhxy_loc[i]->Integral(ibx_min,ibx_max,iby_min,iby_max);
    Float_t reco = fhxy_loc_match[i]->Integral(ibx_min,ibx_max,iby_min,iby_max);
    
    if(ref>0.1){
      
      if(fJura) cout<<" Efficiency in sensor "<<fSensorNames[i].Data()<<"("<<JuraAddStr[i].Data()<<")  =  "<<reco/ref<<endl;
      else      cout<<" Efficiency in sensor "<<fSensorNames[i].Data()<<"("<<SaleveAddStr[i].Data()<<")  =  "<<reco/ref<<endl;
      
    }else{

      if(fJura) cout<<" Can not find efficiency in sensor "<<fSensorNames[i].Data()<<"("<<JuraAddStr[i].Data()<<") ref = "<<ref<<endl;
      else      cout<<" Can not find efficiency in sensor "<<fSensorNames[i].Data()<<"("<<SaleveAddStr[i].Data()<<") ref = "<<ref<<endl;

    }
  }  
  


}



//____________________________________________________________________
void Na61VdEfficiencyModule ::Print(Option_t* option) const
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
         << "    $Date: 2016/10/30$"   << endl 
         << "    $Revision: 1.0 $ " << endl  
         << endl 
         << "-------------------------------------------------" << endl;
}



