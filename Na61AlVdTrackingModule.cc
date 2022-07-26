//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61AlVdTrackingModule
#include "Na61AlVdTrackingModule.h"    
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
//ClassImp(Na61AlVdTrackingModule);

//____________________________________________________________________
Na61AlVdTrackingModule::Na61AlVdTrackingModule()
{
  // Default constructor. DO NOT USE
  SetState(kSetup);
  fProductionMode  = 1;
}

//____________________________________________________________________
Na61AlVdTrackingModule::Na61AlVdTrackingModule(const Char_t* name, const Char_t* title)
   : Na61Module(name, title)
{
  // Named Constructor
  SetState(kSetup);
  fProductionMode  = 1;

  fNsig_dev   = 4; // n sigma cut for 3HitTracks
  fNsig_pvert = 5;  // n sigma cut for matchig of 3HitTracks with primary vertex
  
  fPrimaryVertex = new Vector3D();

  TString gaus = "[0]*exp(-0.5*((x-[1])/[2])^2)";
  fFGauss = new TF1("fFGauss",gaus.Data(),-500,100);  
  fFGauss->SetParNames("const","mean","sigma");    

  for(Int_t i=0;i<34;i++)fHitsTabs[i]=0;
 
}
//____________________________________________________________________
void Na61AlVdTrackingModule::DefineHistograms()
{

  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init"); 
    return;  
    
  }

  cout<<" Na61AlVdTrackingModule::DefineHistograms: "<<fHistDirName.Data()<<endl;  
  TDirectory* histDir = gDirectory->mkdir(fHistDirName.Data());
  histDir->cd(); 

  // these histos should be present in the HitFinder, however we need them here to test/debug simulation
  fhZX_fine =	new TH2F("hZX_fine","",1000,-6,160,1000,-50,50);
  fhX_Al1 =	new TH1F("hX_Al1","",1000,-50,50);
  fhX_Al2 =	new TH1F("hX_Al2","",1000,-50,50);
  fhX_Al3 =	new TH1F("hX_Al3","",1000,-50,50);
  fhX_Al4 =	new TH1F("hX_Al4","",1000,-50,50);
  fhY_Al1 =	new TH1F("hY_Al1","",1000,-80,80);
  fhY_Al2 =	new TH1F("hY_Al2","",1000,-80,80);
  fhY_Al3 =	new TH1F("hY_Al3","",1000,-80,80);
  fhY_Al4 =	new TH1F("hY_Al4","",1000,-80,80);
  
  fhAxAtd = new TH1F(Form("hAxAtd"),"", 2500,-0.25,0.25);
  fhAyAtd = new TH1F(Form("hAyAtd"),"", 2500,-0.25,0.25);
  fhAxAtd_prod = new TH1F(Form("hAxAtd_prod"),"", 2500,-0.25,0.25);
  fhAyAtd_prod = new TH1F(Form("hAyAtd_prod"),"", 2500,-0.25,0.25);
  
  fhCurvature_front = new TH1F("hCurvature_front","",500,-0.00003,0.00003);  
  fhCurvature_back = new TH1F("hCurvature_back","",500,-0.00003,0.00003);      
  fhCurvature_back_prod = new TH1F("hCurvature_back_prod","",500,-0.00003,0.00003);    
  fhSlopeChange = new TH1F("hSlopeChange","",500,-0.008,0.008);  

  // 4 hit tracks without weigths
  fhRecoVertexZ  = new  TH1F("hRecoVertexZ"," ", 2000,-1500.,100.);
  //fhRecoVertexZ  = new  TH1F("hRecoVertexZ"," ", 2000,-500.,100.);
  fhRecoVertexZ_fine  = new  TH1F("hRecoVertexZ_fine"," ", 1000,-90.,-40.);
  
  fhRecoVertexZ_fine_x  = new  TH1F("hRecoVertexZ_fine_x"," ", 1000,-90.,-40.);
  fhRecoVertexZ_fine_y  = new  TH1F("hRecoVertexZ_fine_y"," ", 1000,-90.,-40.);
  fhRecoVertexZ_x  = new  TH1F("hRecoVertexZ_x"," ", 2000,-1500.,100.);
  fhRecoVertexZ_y  = new  TH1F("hRecoVertexZ_y"," ", 2000,-1500.,100.);
  
  fhd2  = new  TH1F("hd2","",500,0.,100.);
  
  // 4 hit tracks with weigths
  fhRecoVertexZ_w2  = new  TH1F("hRecoVertexZ_w2"," ", 2000,-1500.,100.);
  //fhRecoVertexZ_w2  = new  TH1F("hRecoVertexZ_w2"," ", 2000,-500.,100.);
  fhRecoVertexZ_fine_w2  = new  TH1F("hRecoVertexZ_fine_w2"," ", 1000,-90.,-40.);
  fhRecoVertexXY_w2 = new TH2F("hRecoVertexXY_w2"," ", 500,-50.,50.,500,-50.,50.);  
  
  // All tracks with weigths
  fhRecoVertexZ_w2_Post  = new  TH1F("hRecoVertexZ_w2_Post"," ", 2000,-1500.,100.);
  //fhRecoVertexZ_w2_Post  = new  TH1F("hRecoVertexZ_w2_Post"," ", 2000,-500.,100.);
  fhRecoVertexZ_fine_w2_Post  = new  TH1F("hRecoVertexZ_fine_w2_Post"," ", 1000,-90.,-40.);
  fhRecoVertexXY_w2_Post = new TH2F("hRecoVertexXY_w2_Post"," ", 500,-50.,50.,500,-90.,40.);
  
  fhTracks_prod  = new  TH1F("hTracks_prod","", 100,0.,100.);
  fhTracks2_prod  = new  TH1F("hTracks2_prod","", 100,-0.5,99.5);

  fhDZmin  = new  TH1F("hDZmin"," ", 2000,-100.,100.);
  
  //gDirectory->cd(".."); // only for Atd 
  //return;

  for(Int_t i=0;i<4;i++){
    //fdx3hit[i] = new TH1F(Form("hdx3hit_Vds%d",i+1),"", 500,-0.2,0.2);
    //fdx4hit_front[i] = new TH1F(Form("hdx4hit_front_Vds%d",i+1),"", 500,-0.2,0.2);
    fdx4hit_back[i] = new TH1F(Form("hdx4hit_back_Vds%d",i+1),"", 500,-0.2,0.2);
    //cout<<"  fdx4hit_back[i]:  i="<<i<<"  "<<fdx4hit_back[i]<<endl;
  }
  
  if(fHistDirName.Contains("HT")){ // VdTrackingHT dooes not need rest of histos
    gDirectory->cd("..");  
    return;
  }

  
  Int_t nbins = 5000;
  Float_t xmin = -5.0;
  Float_t xmax =  5.0;

  histDir = gDirectory->mkdir("devA");
  histDir->cd(); 

  for(Int_t i=0;i<16;i++){
    TString devname = fDevName[i];
    fhX_dev[i] = new TH1F(Form("hX_%s",devname.Data()),"", nbins, xmin, xmax);
    fhY_dev[i] = new TH1F(Form("hY_%s",devname.Data()),"", nbins, xmin, xmax);
    fhX_dev_cuty[i] = new TH1F(Form("hX_%s_cuty",devname.Data()),"", nbins, xmin, xmax);
    fhY_dev_cutx[i] = new TH1F(Form("hY_%s_cutx",devname.Data()),"", nbins, xmin, xmax);

    fh_devx_devy[i]= new TH2F(Form("h_devx_devy_%s",devname.Data()),"",100,-0.1,0.1,100,-0.1,0.1);

    fhX_dev_acce[i] = new TH1F(Form("hX_%s_acce",devname.Data())," ", nbins, xmin, xmax);
    fhY_dev_acce[i] = new TH1F(Form("hY_%s_acce",devname.Data()),"", nbins, xmin, xmax);

    fhAx[i] = new TH1F(Form("hAx_%s",devname.Data()),"", 2500,-0.25,0.25);
    fhAy[i] = new TH1F(Form("hAy_%s",devname.Data()),"", 2500,-0.25,0.25);

    fhChi2X[i]  = new  TH1F(Form("hChi2X_%s",devname.Data()),"", 200,0.,100.);
    fhChi2Y[i]  = new  TH1F(Form("hChi2Y_%s",devname.Data()),"", 200,0.,100.);
    fhTracks[i] = new TH1F(Form("hTracks_%s",devname.Data()),"", 100,0.,100.);
  }
  /*
  histDir->cd(".."); 
  
  histDir = gDirectory->mkdir("devB");
  histDir->cd();

  for(Int_t i=20;i<45;i++){
    TString devname = fDevName[i];
    fhX_dev[i] = new TH1F(Form("hX_%s",devname.Data()),"", nbins, xmin, xmax);
    fhY_dev[i] = new TH1F(Form("hY_%s",devname.Data()),"", nbins, xmin, xmax);
    fhX_dev_cuty[i] = new TH1F(Form("hX_%s_cuty",devname.Data()),"", nbins, xmin, xmax);
    fhY_dev_cutx[i] = new TH1F(Form("hY_%s_cutx",devname.Data()),"", nbins, xmin, xmax);

    fh_devx_devy[i]= new TH2F(Form("h_devx_devy_%s",devname.Data()),"",100,-0.1,0.1,100,-0.1,0.1);

    fhX_dev_acce[i] = new TH1F(Form("hX_%s_acce",devname.Data())," ", nbins, xmin, xmax);
    fhY_dev_acce[i] = new TH1F(Form("hY_%s_acce",devname.Data()),"", nbins, xmin, xmax);

    fhAx[i] = new TH1F(Form("hAx_%s",devname.Data()),"", 2500,-0.25,0.25);
    fhAy[i] = new TH1F(Form("hAy_%s",devname.Data()),"", 2500,-0.25,0.25);

    fhChi2X[i]  = new  TH1F(Form("hChi2X_%s",devname.Data()),"", 200,0.,100.);
    fhChi2Y[i]  = new  TH1F(Form("hChi2Y_%s",devname.Data()),"", 200,0.,100.);
    fhTracks[i] = new TH1F(Form("hTracks_%s",devname.Data()),"", 100,0.,100.);
  }
  */

  histDir->cd(".."); 
  
  histDir = gDirectory->mkdir("devD");
  histDir->cd();

  for(Int_t i=16;i<fNDev;i++){
    TString devname = fDevName[i];
    fhX_dev[i] = new TH1F(Form("hX_%s",devname.Data()),"", nbins, xmin, xmax);
    fhY_dev[i] = new TH1F(Form("hY_%s",devname.Data()),"", nbins, xmin, xmax);
    fhX_dev_cuty[i] = new TH1F(Form("hX_%s_cuty",devname.Data()),"", nbins, xmin, xmax);
    fhY_dev_cutx[i] = new TH1F(Form("hY_%s_cutx",devname.Data()),"", nbins, xmin, xmax);

    fh_devx_devy[i]= new TH2F(Form("h_devx_devy_%s",devname.Data()),"",100,-0.1,0.1,100,-0.1,0.1);

    fhX_dev_acce[i] = new TH1F(Form("hX_%s_acce",devname.Data())," ", nbins, xmin, xmax);
    fhY_dev_acce[i] = new TH1F(Form("hY_%s_acce",devname.Data()),"", nbins, xmin, xmax);

    fhAx[i] = new TH1F(Form("hAx_%s",devname.Data()),"", 2500,-0.25,0.25);
    fhAy[i] = new TH1F(Form("hAy_%s",devname.Data()),"", 2500,-0.25,0.25);

    fhChi2X[i]  = new  TH1F(Form("hChi2X_%s",devname.Data()),"", 200,0.,100.);
    fhChi2Y[i]  = new  TH1F(Form("hChi2Y_%s",devname.Data()),"", 200,0.,100.);
    fhTracks[i] = new TH1F(Form("hTracks_%s",devname.Data()),"", 100,0.,100.);
  }

  gDirectory->cd("..");  
  
  fhFullTracks = new TH1F(Form("hFullTracks"),"", 1000,0.,1000.);
  fhFullTracksAll = new TH1F(Form("hFullTracksAll"),"", 1000,0.,1000.);

  histDir = gDirectory->mkdir("DevComb");
  histDir->cd();

  for(Int_t i=0;i<fNDevComb;i++){
    TString matchstr = fDevCombName[i];

    fhAx_full[i] = new TH1F(Form("hAx_full_%s",matchstr.Data()),"",2500,-0.25,0.25);
    fhAy_full[i] = new TH1F(Form("hAy_full_%s",matchstr.Data()),"",2500,-0.25,0.25);
    fhAy_full_prod[i] = new TH1F(Form("hAy_full_prod_%s",matchstr.Data()),"",2500,-0.25,0.25);

    fhX_dev_full[i] = new TH1F(Form("hX_dev_full_%s",matchstr.Data()),"",1000,-0.25,0.25);
    fhY_dev_full[i] = new TH1F(Form("hY_dev_full_%s",matchstr.Data()),"",1000,-0.25,0.25);

    fhX_dev_Vds1[i] = new TH1F(Form("hX_dev_Vds1_%s",matchstr.Data()),"",1000,-0.25,0.25);
    fhX_dev_Vds2[i] = new TH1F(Form("hX_dev_Vds2_%s",matchstr.Data()),"",1000,-0.25,0.25);
    fhX_dev_Vds3[i] = new TH1F(Form("hX_dev_Vds3_%s",matchstr.Data()),"",1000,-0.25,0.25);
    fhX_dev_Vds4[i] = new TH1F(Form("hX_dev_Vds4_%s",matchstr.Data()),"",1000,-0.25,0.25);

    fhY_dev_Vds1[i] = new TH1F(Form("hY_dev_Vds1_%s",matchstr.Data()),"",1000,-0.25,0.25);
    fhY_dev_Vds2[i] = new TH1F(Form("hY_dev_Vds2_%s",matchstr.Data()),"",1000,-0.25,0.25);
    fhY_dev_Vds3[i] = new TH1F(Form("hY_dev_Vds3_%s",matchstr.Data()),"",1000,-0.25,0.25);
    fhY_dev_Vds4[i] = new TH1F(Form("hY_dev_Vds4_%s",matchstr.Data()),"",1000,-0.25,0.25);
    
    fhDx_match[i] = new TH1F(Form("hDx_match_%s",matchstr.Data()),"",5000,-10,10);
    fhDy_match[i] = new TH1F(Form("hDy_match_%s",matchstr.Data()),"",5000,-10,10);
    fhDax_match[i] = new TH1F(Form("hDax_match_%s",matchstr.Data()),"",5000,-0.5,0.5);
    fhDay_match[i] = new TH1F(Form("hDay_match_%s",matchstr.Data()),"",5000,-0.5,0.5);
    fhDx_match_cut1[i] = new TH1F(Form("hDx_match_cut1_%s",matchstr.Data()),"",5000,-10,10);
    fhDy_match_cut1[i] = new TH1F(Form("hDy_match_cut1_%s",matchstr.Data()),"",5000,-10,10);
    fhDax_match_cut1[i] = new TH1F(Form("hDax_match_cut1_%s",matchstr.Data()),"",5000,-0.5,0.5);
    fhDay_match_cut1[i] = new TH1F(Form("hDay_match_cut1_%s",matchstr.Data()),"",5000,-0.5,0.5);
    fhDx_match_acce[i] = new TH1F(Form("hDx_match_%s_acce",matchstr.Data()),"",1000,-2,2);
    fhDy_match_acce[i] = new TH1F(Form("hDy_match_%s_acce",matchstr.Data()),"",1000,-2,2);
    fhDax_match_acce[i] = new TH1F(Form("hDax_match_%s_acce",matchstr.Data()),"",1000,-0.05,0.05);
    fhDay_match_acce[i] = new TH1F(Form("hDay_match_%s_acce",matchstr.Data()),"",1000,-0.05,0.05);
  }


  histDir = gDirectory->mkdir("3HitDev_CrossCheck");
  histDir->cd();


  for(Int_t i=0;i<fNDev;i++){
    TString devname = fDevName[i];

    fhX_dev_acce[i] = new TH1F(Form("hX_%s_acce",devname.Data())," ", 1000,-1.,1.);
    fhY_dev_acce[i] = new TH1F(Form("hY_%s_acce",devname.Data())," ", 1000,-1.,1.);
  }     

  gDirectory->cd("..");    
  gDirectory->cd("..");  

  fhClusterSize = new TH1F(Form("hClusterSize"),"",350, 0, 350);

  fhChi2Ndf_x  = new  TH1F("hChi2Ndf_x"," ", 200,0.,100.);
  fhChi2Ndf_y  = new  TH1F("hChi2Ndf_y"," ", 200,0.,100.);

  gDirectory->cd("..");  

}

//____________________________________________________________________
void Na61AlVdTrackingModule::Init()
{
  // Job-level initialisation
  SetState(kInit);
  
  if(fRunId<38014){
    fFieldRun=kTRUE;
  }else if(fRunId<40000){
    fFieldRun=kTRUE; 
  }else{
    fFieldRun=kFALSE;
  }

  //fZprim = 47; //savd
  fZprim = 73.1;

  if(fJura)fParams = (Na61VdParametersManager::Instance())->GetJuraAlVdArmParams();  
  else fParams = (Na61VdParametersManager::Instance())->GetSaleveAlVdArmParams(); 

  //cout<<"Na61AlVdTrackingModule::Init(): fNsig_dev="<<fNsig_dev<<"   fNsig_pvert="<<fNsig_pvert
  //  <<" fFieldRun="<<fFieldRun<<" fJura = "<<(Int_t)fJura<<endl; 

  for(Int_t i=0;i<20;i++)fAll3StationTracks[i] = 0;
 

}

//____________________________________________________________________
void Na61AlVdTrackingModule::Begin()
{
  // Run-level initialisation
  SetState(kBegin);
  
}

//____________________________________________________________________
void Na61AlVdTrackingModule :: Event(UEventNode* inNode, UEventNode* outNode)

{
  if(!inNode)cout<<"inNode is empty"<<endl;
  if(!outNode)cout<<"outNode is empty"<<endl;
  // Per event method
  SetState(kEvent);  

  
}


//_____________________________________________________________
void Na61AlVdTrackingModule::RefitTracks(Int_t Ntab, UEventNode* out)
{
  // algorithm:
  // Ntab=4 to use only 4 hits tracks
  // Ntab=18 to use all tracks

 
  for(Int_t i=0;i<Ntab;i++){
    UDataTable* tracktab = out -> GetDataTable(Form("FullTracks %s",fDevCombName[i].Data()));

    if(!tracktab)continue;    

    for(Int_t j=0; j<tracktab->GetEntries(); j++){
      UVdTrack* track = (UVdTrack*)tracktab->At(j);
  
      RefitTrack_front(track);     
      RefitTrack_back(track);

      //if(Ntab==4)continue;

      Float_t front_slope = track->GetDX_f()/track->GetDZ_f();
      Float_t back_slope  = track->GetDX_b()/track->GetDZ_b();
      fhSlopeChange -> Fill(back_slope - front_slope);  
      track->SetSlopeChange(back_slope - front_slope);
      track->SetLinefParams();
      track->SetLinebParams();
      
    }
  }
  
}

//_______________________________________________________________
void Na61AlVdTrackingModule::FillHitPositions(UEventNode* out)
{

  for(Int_t it=0;it<34;it++){
    UDataTable* hits = out->GetDataTable(Form("Hits %s",fAlSensorNames[it].Data()));

    for(Int_t i=0;i<hits->GetEntries();i++){
      USensorHit* hit = (USensorHit*)hits->At(i);      

      Float_t x = hit->GetX();
      Float_t y = hit->GetY();
      Float_t z = hit->GetZ();
      
 
      fhZX_fine->Fill(z,x);
      if(fAlSensorNames[it].Contains("Al1_")){fhX_Al1->Fill(x); fhY_Al1->Fill(y);} 
      if(fAlSensorNames[it].Contains("Al2_")){fhX_Al2->Fill(x); fhY_Al2->Fill(y);} 
      if(fAlSensorNames[it].Contains("Al3_")){fhX_Al3->Fill(x); fhY_Al3->Fill(y);} 
      if(fAlSensorNames[it].Contains("Al4_")){fhX_Al4->Fill(x); fhY_Al4->Fill(y);} 


    }
  }

}

//_______________________________________________________________
void Na61AlVdTrackingModule::RefitTracks(UEventNode* out)
{
  // algorithm:

  UDataTable* tracktab = out -> GetDataTable(fOutTableName.Data());

  if(!tracktab)return;
  //cout<<"start with table, entries: "<<tracktab->GetEntries()<<endl;    
  
  for(Int_t i=0; i<tracktab->GetEntries(); i++){
    UVdTrack* track = (UVdTrack*)tracktab->At(i);

    RefitTrack_front(track);
    RefitTrack_back(track);

    Float_t front_slope = track->GetDX_f()/track->GetDZ_f();
    Float_t back_slope  = track->GetDX_b()/track->GetDZ_b();
    fhSlopeChange -> Fill(back_slope - front_slope);  
    track->SetSlopeChange(back_slope - front_slope);
    track->SetLinefParams();
    track->SetLinebParams();

  }
}
  

//____________________________________________________________________________________
void Na61AlVdTrackingModule::RefitTrack_front(UVdTrack* track)
{
  
  //cout<<"Na61AlVdTrackingModule::RefitTrack_front: entering"<<endl;

  Int_t ind1 = track->GetHitIndexOnStation(0);
  Int_t ind2 = track->GetHitIndexOnStation(1);
  Int_t ind3 = track->GetHitIndexOnStation(2);
  Int_t ind4 = track->GetHitIndexOnStation(3);

  Int_t tind1 = track->GetTabIndexOnStation(0);
  Int_t tind2 = track->GetTabIndexOnStation(1);
  Int_t tind3 = track->GetTabIndexOnStation(2);
  Int_t tind4 = track->GetTabIndexOnStation(3);
  

  USensorHit* hit1 = 0;
  USensorHit* hit2 = 0;
  USensorHit* hit3 = 0;
  USensorHit* hit4 = 0;

  if((ind1 != -1) && (tind1 != -1) && fHitsTabs[tind1]) hit1 = (USensorHit*)fHitsTabs[tind1]->At(ind1);
  if((ind2 != -1) && (tind2 != -1) && fHitsTabs[tind2]) hit2 = (USensorHit*)fHitsTabs[tind2]->At(ind2);
  if((ind3 != -1) && (tind3 != -1) && fHitsTabs[tind3]) hit3 = (USensorHit*)fHitsTabs[tind3]->At(ind3);
  if((ind4 != -1) && (tind4 != -1) && fHitsTabs[tind4]) hit4 = (USensorHit*)fHitsTabs[tind4]->At(ind4);

  // need it for tracking with HT (should we added standard hits association in VdTrackingHTModule?) 
  if(!hit1) hit1 = (USensorHit*)track->GetHitIdAtStation(0);
  if(!hit2) hit2 = (USensorHit*)track->GetHitIdAtStation(1);
  if(!hit3) hit3 = (USensorHit*)track->GetHitIdAtStation(2);
  if(!hit4) hit4 = (USensorHit*)track->GetHitIdAtStation(3);
  
  

   // for debugging purpose

  if(!hit1)cout<<"hit1 is null"<<endl;
  if(!hit2)cout<<"hit2 is null"<<endl;
  if(!hit3)cout<<"hit3 is null"<<endl;
  if(!hit4)cout<<"hit4 is null"<<endl;
  
  Float_t sig[4] = {0.004, 0.0145, 0.032, 0.0416}; // errors in mm
  //Float_t sig[4] = {0.004, 0.005, 0.006, 0.007}; // errors in mm
  
  USensorHit* hits[4] = {hit1,hit2,hit3,hit4};

  for(Int_t i=0;i<4;i++){
    if(hits[i])hits[i]->SetUsed(kTRUE);  // do it only in Refit_front
  }

  Float_t ay;
  Float_t by;
  Float_t N;
  Float_t zmin;

  FitLineY_w2(hits,sig,ay,by,N,zmin);


  TF1 pol2y("pol2y","[0] + [1]*x + [2]*x*x",-10.0,160.0);


  Double_t x[4]; 
  Double_t y[4];
  Double_t ex[4]; 
  Double_t ey[4];

  Int_t ii=0;
  for(Int_t i=0;i<4;i++){

    if(hits[i]){
      x[ii] = hits[i]->GetZ(); 
     
      ex[ii] = 0.001; // 1 micron error on z position

      y[ii] = hits[i]->GetX(); 
      ey[ii] = sig[i];
 
      ii++;
    }
  }
 
  TGraphErrors* gr = new TGraphErrors(ii,x,y,ex,ey);
  
  //if(x[0]>10)cout<<"Tu10: "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<endl;
  
  pol2y.SetParameter(0,track->GetX());
  pol2y.SetParameter(1,track->GetDX()/track->GetDZ());
  if(fFieldRun)pol2y.SetParameter(2,0);
  else pol2y.FixParameter(2,0);

  
  gr->Fit("pol2y","Q" ,"",x[0]-1,x[ii-1]+1);
  
  //cout<<"Tu12: "<<fFieldRun<<endl;
  
  Double_t c = pol2y.GetParameter(0);
  Double_t b = pol2y.GetParameter(1);
  Double_t a = pol2y.GetParameter(2);

  if(fhCurvature_front)fhCurvature_front->Fill(a);   // don't fill if called from CombineWithPrimaryVertex
  
  //if(x[0]>10)cout<<"Tu10: "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<"  a="<<a<<endl;
  //cout<<"Tu10: "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<"  a="<<a<<endl;

  Float_t oz = track->GetZ();
  Float_t oy = ay*oz + by;
  Float_t ox = a*oz*oz + b*oz + c;

  //cout<<"oz="<<oz<<endl;

  //cout<<"Tu13"<<" "<<track->GetDZ()<<" "<<a<<endl;

  Double_t tanx = 2.*a*oz + b;
  Double_t tany = ay;

  //cout<<"Tu14: line"<<(Line3D*)track->Getline()<<endl;

 (track->Getlinef())->SetOrigin(ox,oy,oz);
 (track->Getlinef())->SetDirection(tanx,tany,1.);

 //cout<<"Tu15: linef"<<(Line3D*)track->Getlinef()<<endl;

 ox = track->GetXatZ_f(0.0);
 oy = track->GetYatZ_f(0.0);
 (track->Getlinef())->SetOrigin(ox,oy,0.0);


 // check fit quality


 delete gr;

 // we can try to calculate momentum here.

 //if(iflag != 1)return; // don't fill if called from CombineWithPrimaryVertex
   
 /* dont know why we need it (PS)
 for(Int_t i=0;i<4;i++){
   if(hits[i]){
     Float_t z = hits[i]->GetZ(); 
     Float_t x = hits[i]->GetX();
     
     if(ii==3)fdx3hit[i] -> Fill(x - pol2y.Eval(z));
     else fdx4hit_front[i] -> Fill(x - pol2y.Eval(z));
     
   }
 }
 */
 
 //delete pol2y;
 
}


//____________________________________________________________________________________
void Na61AlVdTrackingModule::RefitTrack_back(UVdTrack* track)
{
 
  Int_t ind1 = track->GetHitIndexOnStation(0);
  Int_t ind2 = track->GetHitIndexOnStation(1);
  Int_t ind3 = track->GetHitIndexOnStation(2);
  Int_t ind4 = track->GetHitIndexOnStation(3);

  Int_t tind1 = track->GetTabIndexOnStation(0);
  Int_t tind2 = track->GetTabIndexOnStation(1);
  Int_t tind3 = track->GetTabIndexOnStation(2);
  Int_t tind4 = track->GetTabIndexOnStation(3);
  

  USensorHit* hit1 = 0;
  USensorHit* hit2 = 0;
  USensorHit* hit3 = 0;
  USensorHit* hit4 = 0;


  if((ind1 != -1) && (tind1 != -1) && fHitsTabs[tind1]) hit1 = (USensorHit*)fHitsTabs[tind1]->At(ind1);
  if((ind2 != -1) && (tind1 != -1) && fHitsTabs[tind2]) hit2 = (USensorHit*)fHitsTabs[tind2]->At(ind2);
  if((ind3 != -1) && (tind1 != -1) && fHitsTabs[tind3]) hit3 = (USensorHit*)fHitsTabs[tind3]->At(ind3);
  if((ind4 != -1) && (tind1 != -1) && fHitsTabs[tind4]) hit4 = (USensorHit*)fHitsTabs[tind4]->At(ind4);

  // need it for tracking with HT (should we add standard hits association in VdTrackingHTModule?) 
  if(!hit1) hit1 = (USensorHit*)track->GetHitIdAtStation(0);
  if(!hit2) hit2 = (USensorHit*)track->GetHitIdAtStation(1);
  if(!hit3) hit3 = (USensorHit*)track->GetHitIdAtStation(2);
  if(!hit4) hit4 = (USensorHit*)track->GetHitIdAtStation(3);

  //Float_t sig[4] = {0.004, 0.0145, 0.032, 0.0416}; // errors in mm
  Float_t sig[4] = {0.007, 0.006, 0.005, 0.004}; // errors in mm
  
 
  USensorHit* hits[4] = {hit1,hit2,hit3,hit4};

  Float_t ay;
  Float_t by;
  Float_t N;
  Float_t zmin;

  FitLineY_w2(hits,sig,ay,by,N,zmin);

  //TF1* pol2y = new TF1("pol2y","[0] + [1]*x + [2]*x*x",-10.0,160.0);
  TF1 pol2y("pol2y","[0] + [1]*x + [2]*x*x",-10.0,160.0);


  Double_t x[4]; 
  Double_t y[4];
  Double_t ex[4]; 
  Double_t ey[4];

  Int_t ii=0;
  for(Int_t i=0;i<4;i++){

    if(hits[i]){
      x[ii] = hits[i]->GetZ(); 
      ex[ii] = 0.001; // 1 micron error on z position

      y[ii] = hits[i]->GetX(); 
      ey[ii] = sig[i];
 
      ii++;
    }
  }

  TGraphErrors* gr = new TGraphErrors(ii,x,y,ex,ey);
  
  //cout<<"Tu3: "<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<x[3]<<endl;
  
  pol2y.SetParameter(0,track->GetX());
  pol2y.SetParameter(1,track->GetDX()/track->GetDZ());
  if(fFieldRun)pol2y.SetParameter(2,0);
  else pol2y.FixParameter(2,0);

  
  gr->Fit("pol2y","Q" ,"",x[0]-1,x[ii-1]+1);
  
  Double_t c = pol2y.GetParameter(0);
  Double_t b = pol2y.GetParameter(1);
  Double_t a = pol2y.GetParameter(2);

  // x[ii-1] contains z coordinate of the last hit
  Float_t oz = x[ii-1];
  Float_t oy = by + ay*oz;
  Float_t ox = a*oz*oz + b*oz + c;

  //cout<<"Tu4"<<" "<<track->GetDZ()<<" "<<a<<" zlast="<<oz<<endl;

  Double_t tanx = 2.*a*oz + b;
  Double_t tany = ay;

  // set track charge
  if(fhCurvature_back)fhCurvature_back->Fill(a);  
  track->SetCurvature(a);
  track->SetPol2(pol2y);

  if(a>0)track->SetCharge(1);
  else track->SetCharge(-1);

  //cout<<"Tu5: line"<<(Line3D*)track->Getline()<<endl;

 (track->Getlineb())->SetOrigin(ox,oy,oz);
 (track->Getlineb())->SetDirection(tanx,tany,1.);
 //cout<<"Tu6: lineb"<<(Line3D*)track->Getlineb()<<endl;


   for(Int_t i=0;i<4;i++){
     if(hits[i]){
       Float_t z = hits[i]->GetZ(); 
       Float_t x = hits[i]->GetX();
       //cout<<"i="<<i<<"  "<<fdx4hit_back[i]<<endl;
       if(ii==4)fdx4hit_back[i] -> Fill(x - pol2y.Eval(z));

     }
   }
   
   //delete pol2y;
   delete gr;

}



//___________________________________________________________________________________
void Na61AlVdTrackingModule::FitLineY_w2(USensorHit** hits,Float_t* sig, Float_t& ay,Float_t& by,
				       Float_t& N, Float_t& zmin)
{
  // simple regration (assuming constant errors)
  // It is quite straight forward to use weighted regression

    Float_t hy,hz;

    Float_t S=0;
    Float_t Sz=0;
    Float_t Szz=0;
 
    Float_t Sy=0;
    Float_t Szy=0;
     
    zmin = 10000.;
    N = 0;    
    
    // perhaps we need waited regression
    
    for(Int_t i=0;i<4;i++){      
      USensorHit* hit = hits[i];
      
      if(!hit)continue; // now we always have 4 hits
      
      N++;
      
      // smearing is done in dedicated method
      hy = hit->GetY();
      hz = hit->GetZ(); 
      
      Float_t sig2 = sig[i]*sig[i];     

      S   += 1./sig2;
      Sz  += hz/sig2;
      Szz += (hz*hz)/sig2;

      Sy  += hy/sig2;
      Szy += (hz*hy)/sig2;
      
      if(hz<zmin)zmin=hz;
      
    }
    
    ay = (S*Szy - Sz*Sy)/(S*Szz-Sz*Sz);
    by = (Szz*Sy - Sz*Szy)/(S*Szz-Sz*Sz);
    
}


//_____________________________________________________________
Int_t Na61AlVdTrackingModule::GetAllTracks(UEventNode* out)
{
  // algorithm:

  Int_t fullTracks = 0;
 

  for(Int_t i=0;i<18;i++){
    UDataTable* tracktab = out -> GetDataTable(Form("FullTracks %s",fMatchStr[i].Data()));

    if(!tracktab)continue;    

    for(Int_t j=0; j<tracktab->GetEntries(); j++){
      //UVdTrack* track = (UVdTrack*)tracktab->At(j);
      
      fullTracks++;

    }
  }
  
  return fullTracks;
  
}

//_____________________________________________________________
void Na61AlVdTrackingModule::Delete3HitTracks(UEventNode* out)
{  

  Int_t NDev = 20;

  for(Int_t i=0;i<NDev;i++){
    UDataTable* table = out->GetDataTable(Form("Tracks %s",fDevNames[i].Data()));
  
    if(table){
      //cout<<" table name: "<<table->GetName()<<endl;
      //table->Clear();
      table->Delete();
      if(table)out->RemoveDataObject((UDataObject*)table);
    }
  }

}


//_____________________________________________________________
void Na61AlVdTrackingModule::CombineTracks(Int_t imatch, Int_t idev1, Int_t idev2, UEventNode* out)
{
  // combining 3hit tracks to compose full 4hit tracks
  
  //cout<<" Na61AlVdTrackingModule::CombineTracks: start, imatch="<<imatch<<"  idev1="<<idev1<<"  idev2="<<idev2<<endl;

  UDataTable* tracks_dev1 = out->GetDataTable(Form("Tracks %s",fDevName[idev1].Data()));
  UDataTable* tracks_dev2 = out->GetDataTable(Form("Tracks %s",fDevName[idev2].Data()));

  if(!tracks_dev1)return;
  if(!tracks_dev2)return;

  ///////////////////// select indexes for HitIdAtStation //////////////////////////////////
  Int_t ind1[3];
  Int_t ind2[3];
  if(fDevName[idev1].Contains("A")){ind1[0]=0; ind1[1]=1; ind1[2]=2;}
  if(fDevName[idev1].Contains("B")){ind1[0]=0; ind1[1]=2; ind1[2]=3;}
  if(fDevName[idev1].Contains("C")){ind1[0]=0; ind1[1]=1; ind1[2]=3;}
  if(fDevName[idev1].Contains("D")){ind1[0]=1; ind1[1]=2; ind1[2]=3;}
  if(fDevName[idev2].Contains("A")){ind2[0]=0; ind2[1]=1; ind2[2]=2;}
  if(fDevName[idev2].Contains("B")){ind2[0]=0; ind2[1]=2; ind2[2]=3;}
  if(fDevName[idev2].Contains("C")){ind2[0]=0; ind2[1]=1; ind2[2]=3;}
  if(fDevName[idev2].Contains("D")){ind2[0]=1; ind2[1]=2; ind2[2]=3;}
  ////////////////////////////////////////////////////////////////////////////////////////////

  UDataTable* tracktab = new UDataTable(Form("FullTracks %s",fDevCombName[imatch].Data()));
  tracktab->SetOwner();

  //Float_t Off_dx  = fParams->GetDxOffset(imatch);    Float_t Sig_dx  = fParams->GetDxSigma(imatch); 
  //Float_t Off_dy  = fParams->GetDyOffset(imatch);    Float_t Sig_dy  = fParams->GetDySigma(imatch); 
  //Float_t Off_dax = fParams->GetDaxOffset(imatch);   Float_t Sig_dax = fParams->GetDaxSigma(imatch); 
  //Float_t Off_day = fParams->GetDayOffset(imatch);   Float_t Sig_day = fParams->GetDaySigma(imatch); 

  Int_t Ndev1 =  tracks_dev1 -> GetEntries();
  Int_t Ndev2 =  tracks_dev2 -> GetEntries();

  //cout<<Ndev1<<" "<<Ndev2<<"   "<<fMatchId1[imatch]<<" "<<fMatchId2[imatch]<<endl;

  for(Int_t i=0; i<Ndev1; i++){
    UVdTrack* track1 = (UVdTrack*) tracks_dev1->At(i);

    //cout<<"i="<<i<<"  "<<(Long_t)track1->GetHitIdAtStation(fMatchId1[imatch])<<"  "<<(Long_t)track1->GetHitIdAtStation(fMatchId2[imatch])<<endl;

    Float_t tx1 = track1->GetXatZ(125.);
    Float_t ty1 = track1->GetYatZ(125.);
    Float_t ax1 = track1->GetDX()/track1->GetDZ();
    Float_t ay1 = track1->GetDY()/track1->GetDZ();

    Long_t id1 = (Long_t)track1->GetHitIdAtStation(fMatchId1[imatch]);
    Long_t id2 = (Long_t)track1->GetHitIdAtStation(fMatchId2[imatch]);


    for(Int_t j=0; j<Ndev2; j++){
      UVdTrack* track2 = (UVdTrack*) tracks_dev2->At(j);

      //cout<<"j="<<j<<"  "<<(Long_t)track2->GetHitIdAtStation(fMatchId1[imatch])<<"  "<<(Long_t)track2->GetHitIdAtStation(fMatchId2[imatch])<<endl;
      
      Long_t idd1 = (Long_t)track2->GetHitIdAtStation(fMatchId1[imatch]);
      Long_t idd2 = (Long_t)track2->GetHitIdAtStation(fMatchId2[imatch]);
            
      Float_t tx2 = track2->GetXatZ(125.);
      Float_t ty2 = track2->GetYatZ(125.);
      Float_t ax2 = track2->GetDX()/track2->GetDZ();
      Float_t ay2 = track2->GetDY()/track2->GetDZ();

      Float_t dx =  tx2-tx1;
      Float_t dy =  ty2-ty1;
      Float_t dax = ax2-ax1;
      Float_t day = ay2-ay1;
      
      //cout<<fhDx_match[imatch]<<" "<<fhDy_match[imatch]<<" "<<fhDax_match[imatch]<<" "<<fhDay_match[imatch]<<endl;

      if(fhDx_match[imatch])  fhDx_match[imatch]  -> Fill(dx);
      if(fhDy_match[imatch])  fhDy_match[imatch]  -> Fill(dy);
      if(fhDax_match[imatch]) fhDax_match[imatch] -> Fill(dax);
      if(fhDay_match[imatch]) fhDay_match[imatch] -> Fill(day);
     
      if(id1 != idd1) continue;

      if(fhDx_match_cut1[imatch])  fhDx_match_cut1[imatch]  -> Fill(dx);
      if(fhDy_match_cut1[imatch])  fhDy_match_cut1[imatch]  -> Fill(dy);
      if(fhDax_match_cut1[imatch]) fhDax_match_cut1[imatch] -> Fill(dax);
      if(fhDay_match_cut1[imatch]) fhDay_match_cut1[imatch] -> Fill(day);

      if(id2 != idd2) continue;


      // we don't base on dev2 cut but keep it for any reason
      //Float_t dev2 = (dx-Off_dx)*(dx-Off_dx)/(Sig_dx*Sig_dx) + (dy-Off_dy)*(dy-Off_dy)/(Sig_dy*Sig_dy) +
      //(dax-Off_dax)*(dax-Off_dax)/(Sig_dax*Sig_dax) + (day-Off_day)*(day-Off_day)/(Sig_day*Sig_day);

      //if(dev2 > 5*5) continue;


      if(fhDx_match_acce[imatch])  fhDx_match_acce[imatch]  -> Fill(dx);
      if(fhDy_match_acce[imatch])  fhDy_match_acce[imatch]  -> Fill(dy);
      if(fhDax_match_acce[imatch]) fhDax_match_acce[imatch] -> Fill(dax);
      if(fhDay_match_acce[imatch]) fhDay_match_acce[imatch] -> Fill(day);
      
      // fill deviation
      //cout<<"tu2:  imatch="<<imatch<<"  "<<fDevName[idev1].Data()<<" "<<fDevName[idev2].Data()<<endl;
      //cout<<ind1[0]<<"  "<<ind1[1]<<" "<<ind1[2]<<endl;
      //cout<<(Long_t)track1->GetHitIdAtStation(0)<<" "<<(Long_t)track1->GetHitIdAtStation(1)<<" "<<(Long_t)track1->GetHitIdAtStation(2)<<" "<<(Long_t)track1->GetHitIdAtStation(3)<<endl;

      Float_t x1 = ((USensorHit*)track1->GetHitIdAtStation(ind1[0]))->GetX();
      Float_t x2 = ((USensorHit*)track1->GetHitIdAtStation(ind1[1]))->GetX();
      Float_t x3 = ((USensorHit*)track1->GetHitIdAtStation(ind1[2]))->GetX();
      Float_t y1 = ((USensorHit*)track1->GetHitIdAtStation(ind1[0]))->GetY();
      Float_t y2 = ((USensorHit*)track1->GetHitIdAtStation(ind1[1]))->GetY();
      Float_t y3 = ((USensorHit*)track1->GetHitIdAtStation(ind1[2]))->GetY();
      Float_t z1 = ((USensorHit*)track1->GetHitIdAtStation(ind1[0]))->GetZ();
      Float_t z2 = ((USensorHit*)track1->GetHitIdAtStation(ind1[1]))->GetZ();
      Float_t z3 = ((USensorHit*)track1->GetHitIdAtStation(ind1[2]))->GetZ();

      //cout<<"tu3"<<endl;

      Float_t devx = ((z2-z1)*x3 + (z3-z2)*x1)/(z3-z1) - x2;	  	  
      Float_t devy = ((z2-z1)*y3 + (z3-z2)*y1)/(z3-z1) - y2;

      if(fhX_dev_acce[idev1]) fhX_dev_acce[idev1]->Fill(devx);       
      if(fhY_dev_acce[idev1]) fhY_dev_acce[idev1]->Fill(devy);

      //if(fhVds1_xy[imatch]) fhVds1_xy[imatch]_Fill(x1,y1);
      //if(fhVds2_xy[imatch]) fhVds2_xy[imatch]_Fill(x2,y2);
      //if(fhVds3_xy[imatch]) fhVds3_xy[imatch]_Fill(x3,y3);

      //cout<<"tu4:  imatch="<<imatch<<endl;

      x1 = ((USensorHit*)track2->GetHitIdAtStation(ind2[0]))->GetX();
      x2 = ((USensorHit*)track2->GetHitIdAtStation(ind2[1]))->GetX();
      x3 = ((USensorHit*)track2->GetHitIdAtStation(ind2[2]))->GetX();
      y1 = ((USensorHit*)track2->GetHitIdAtStation(ind2[0]))->GetY();
      y2 = ((USensorHit*)track2->GetHitIdAtStation(ind2[1]))->GetY();
      y3 = ((USensorHit*)track2->GetHitIdAtStation(ind2[2]))->GetY();
      z1 = ((USensorHit*)track2->GetHitIdAtStation(ind2[0]))->GetZ();
      z2 = ((USensorHit*)track2->GetHitIdAtStation(ind2[1]))->GetZ();
      z3 = ((USensorHit*)track2->GetHitIdAtStation(ind2[2]))->GetZ();
 
      devx = ((z2-z1)*x3 + (z3-z2)*x1)/(z3-z1) - x2;	  	  
      devy = ((z2-z1)*y3 + (z3-z2)*y1)/(z3-z1) - y2;

      //cout<<"tu5: "<<idev2<<endl;
      //cout<<"tu5: "<<idev2<<" "<<fhX_dev_acce[idev2]<<endl;

      if(fhX_dev_acce[idev2]) fhX_dev_acce[idev2]->Fill(devx);       
      if(fhY_dev_acce[idev2]) fhY_dev_acce[idev2]->Fill(devy);

      //if(fhVds4_xy[imatch]) fhVds4_xy[imatch]_Fill(x3,y3);

      // create 4hit tracks here.       

      CreateFullTrack_new(imatch,track1,track2,tracktab);

    }
  }
  
  
  out->AddDataTable(tracktab);
  
  //cout<<imatch<<" "<<fAllFullTracks[imatch]<<" "<<tracktab->GetEntries()<<" "<<fhDx_match_acce[imatch]->GetEntries()<<endl;

  //cout<<" Na61AlVdTrackingModule::CombineTracks: done"<<endl;


}

//_____________________________________________________________
void Na61AlVdTrackingModule::FindPrimaryVertexPost(UEventNode* out)
{
  // algorithm:

  Int_t fullTracks = 0;
  TObjArray tracks;
  tracks.Clear();

  for(Int_t i=0;i<18;i++){
    UDataTable* tracktab = out -> GetDataTable(Form("FullTracks %s",fMatchStr[i].Data()));
    if(!tracktab)continue;
    
    fullTracks = fullTracks + tracktab->GetEntries();
    for(Int_t i=0; i<tracktab->GetEntries(); i++){
      UVdTrack* track = ( UVdTrack*)tracktab->At(i);
      if(!track->IsMarkedForRemoval()) tracks.Add(tracktab->At(i));
    }
  }

  //cout<<"fhFullTracksAll: "<<fhFullTracksAll<<" "<<fullTracks<<endl;

  fhFullTracksAll -> Fill(fullTracks);

  //cout<<"tu1 "<<endl;

  if(tracks.GetEntries()<2) return;

  
  Float_t sum_up_x, sum_do_x, sum_up_y, sum_do_y;

  CalculateSums(sum_up_x,sum_do_x,sum_up_y,sum_do_y,tracks);

  
  ////////// new method //////////////////////////////////////////////////
  /*
  Float_t zx_min = 0;
  if(sum_do_x > 0)zx_min = - sum_up_x/sum_do_x;
  Float_t zy_min = 0;
  if(sum_do_y > 0)zy_min = - sum_up_y/sum_do_y; 
  */
  
  //cout<<sum_do_x<<"  "<<sum_do_y<<"  "<<zx_min<<" "<<zy_min<<endl; 
  

  Float_t z_prim = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);

  Vector3D reco_pvertex_w2;

  PrimaryVertexWithWeigths(tracks,z_prim,reco_pvertex_w2);

  fhRecoVertexXY_w2_Post     -> Fill(reco_pvertex_w2.X(),reco_pvertex_w2.Y());
  fhRecoVertexZ_w2_Post      -> Fill(reco_pvertex_w2.Z());
  fhRecoVertexZ_fine_w2_Post -> Fill(reco_pvertex_w2.Z());

  fPrimaryVertexDefined = kTRUE;
  ((UVdEvent*)out)->SetPrimaryVertexStatus(2);
  ((UVdEvent*)out)->SetPrimaryVertex(reco_pvertex_w2.X(),reco_pvertex_w2.Y(),reco_pvertex_w2.Z());
  fPrimaryVertex->SetX(reco_pvertex_w2.X());
  fPrimaryVertex->SetY(reco_pvertex_w2.Y());
  fPrimaryVertex->SetZ(reco_pvertex_w2.Z());
  

}

//_____________________________________________________________
void Na61AlVdTrackingModule::FindPrimaryVertex(UEventNode* out)
{
   
  Int_t fullTracks = 0;

  TObjArray tracks;
  tracks.Clear();

  for(Int_t it=0;it<fNDevComb;it++){
    UDataTable* tracktab = out -> GetDataTable(Form("FullTracks %s",fDevCombName[it].Data()));
    if(!tracktab)continue;
    fullTracks = fullTracks + tracktab->GetEntries();
    
    for(Int_t i=0; i<tracktab->GetEntries(); i++){
      UVdTrack* track = ( UVdTrack*)tracktab->At(i);
      if(!track->IsMarkedForRemoval()) tracks.Add(tracktab->At(i));
    }
  }
  fhFullTracks -> Fill(fullTracks);
  
  //cout<<"---------------> tracks.GetEntries(): "<<tracks.GetEntries()<<endl;

  fPrimaryVertexDefined = kFALSE;
  ((UVdEvent*)out)->SetPrimaryVertexStatus(0);

  if(tracks.GetEntries() < 2) return;  
  
  Float_t sum_up_x, sum_do_x, sum_up_y, sum_do_y;

  CalculateSums(sum_up_x,sum_do_x,sum_up_y,sum_do_y,tracks);

  ////////// new method //////////////////////////////////////////////////
  
  Float_t zx_min = 0;
  if(sum_do_x > 0)zx_min = - sum_up_x/sum_do_x;
  Float_t zy_min = 0;
  if(sum_do_y > 0)zy_min = - sum_up_y/sum_do_y; 

  //cout<<sum_do_x<<"  "<<sum_do_y<<"  "<<zx_min<<" "<<zy_min<<endl; 
  fhDZmin ->Fill(zy_min - zx_min);
  
  fhRecoVertexZ_fine_x -> Fill(zx_min);
  fhRecoVertexZ_fine_y -> Fill(zy_min);

  fhRecoVertexZ_x -> Fill(zx_min);
  fhRecoVertexZ_y -> Fill(zy_min);

  Float_t z_prim = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);

  fhRecoVertexZ      -> Fill(z_prim);
  fhRecoVertexZ_fine -> Fill(z_prim);

  //cout<<"---------------> z_prim = "<<z_prim<<endl;

  Vector3D reco_pvertex_w2;
  PrimaryVertexWithWeigths(tracks,z_prim,reco_pvertex_w2);

  //cout<<reco_pvertex_w2.X()<<"  "<<reco_pvertex_w2.Y()<<endl;
  fhRecoVertexXY_w2     -> Fill(reco_pvertex_w2.X(),reco_pvertex_w2.Y());
  fhRecoVertexZ_w2      -> Fill(reco_pvertex_w2.Z());
  fhRecoVertexZ_fine_w2 -> Fill(reco_pvertex_w2.Z());

  if(TMath::Abs(reco_pvertex_w2.Z()+fZprim)<3){
    fPrimaryVertexDefined = kTRUE;
    ((UVdEvent*)out)->SetPrimaryVertexStatus(1);
  }
  
  fPrimaryVertex->SetX(reco_pvertex_w2.X());
  fPrimaryVertex->SetY(reco_pvertex_w2.Y());
  fPrimaryVertex->SetZ(reco_pvertex_w2.Z());
  ((UVdEvent*)out)->SetPrimaryVertex(reco_pvertex_w2.X(),reco_pvertex_w2.Y(),reco_pvertex_w2.Z());

  
}


//______________________________________________________________________________________________________
void Na61AlVdTrackingModule::CalculateSums(Float_t& sum_up_x, Float_t& sum_do_x, Float_t& sum_up_y, Float_t& sum_do_y, TObjArray tracks)
{

  sum_up_x = 0;
  sum_do_x = 0;
  sum_up_y = 0;
  sum_do_y = 0;

  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracks.At(i); 
    for(Int_t j=i+1;j<tracks.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracks.At(j); 
                  
      /////////////////// stuff for new method
      Float_t axi = track1->GetDX_f()/track1->GetDZ_f();
      Float_t axj = track2->GetDX_f()/track2->GetDZ_f();
      Float_t bxi = track1->GetX_f();
      Float_t bxj = track2->GetX_f();
      
      Float_t ayi = track1->GetDY_f()/track1->GetDZ();
      Float_t ayj = track2->GetDY_f()/track2->GetDZ();
      Float_t byi = track1->GetY_f();
      Float_t byj = track2->GetY_f();
      
      Float_t z_x = ((axi-axj)*(bxi-bxj))/((axi-axj)*(axi-axj));
      Float_t z_y = ((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj));
      

      //if(TMath::Abs(z_x-fZprim) < 3){ // SAVD
      if(TMath::Abs(z_x-fZprim) < 10){
	sum_up_x = sum_up_x + (axi-axj)*(bxi-bxj);
	sum_do_x = sum_do_x + (axi-axj)*(axi-axj);
      }
      
      //if(TMath::Abs(z_y-fZprim) < 3){ //SAVD
      if(TMath::Abs(z_y-fZprim) < 10){
	sum_up_y = sum_up_y + (ayi-ayj)*(byi-byj);
	sum_do_y = sum_do_y + (ayi-ayj)*(ayi-ayj);
      }
      
      
      //cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
      
    }
  } 

}

//______________________________________________________________________________________________________
void Na61AlVdTrackingModule::PrimaryVertexWithWeigths(TObjArray& tracks,Float_t zprim,Vector3D& pvertex)
{
  // algorithm: 
  // 1. Calculate weights = 1/dist^2, distance betweem i-th and j-th tracks on
  // the closest proximity plane (z=zprim).
  // 2. Used the weights to calculate z coordinate of primary vertex from the formula that inclused weights (see PS notes).  
  // 3. Find x,y coordinates averageing tracks x, y at the closes proximmity plane. 

  Float_t sum_up_x = 0;
  Float_t sum_do_x = 0;
  Float_t sum_up_y = 0;
  Float_t sum_do_y = 0;

  //cout<<" entries="<<tracks.GetEntries()<<" "<<zprim<<endl;

  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracks.At(i); 
    for(Int_t j=i+1;j<tracks.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracks.At(j); 
      
      //cout<<"i="<<i<<" j="<<j<<endl;

      Float_t x1 = track1->GetXatZ(zprim);
      Float_t x2 = track2->GetXatZ(zprim);
      Float_t y1 = track1->GetYatZ(zprim);
      Float_t y2 = track2->GetYatZ(zprim);

      //cout<<"x1="<<x1<<" x2="<<x2<<"  y1="<<y1<<" y2="<<y2<<endl;
      
      Float_t d2 = (x2-x1)*(x2-x1) +  (y2-y1)*(y2-y1);      
      Float_t w = 1./d2;
      
      //cout<<"dd2="<<d2<<endl;
 
      fhd2->Fill(d2);

      //cout<<"d2="<<d2<<endl;

      /////////////////// stuff for new method
      Float_t axi = track1->GetDX()/track1->GetDZ();
      Float_t axj = track2->GetDX()/track2->GetDZ();
      Float_t bxi = track1->GetX();
      Float_t bxj = track2->GetX();
      
      Float_t ayi = track1->GetDY()/track1->GetDZ();
      Float_t ayj = track2->GetDY()/track2->GetDZ();
      Float_t byi = track1->GetY();
      Float_t byj = track2->GetY();
      
      Float_t z_x = ((axi-axj)*(bxi-bxj))/((axi-axj)*(axi-axj));
      Float_t z_y = ((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj));
      
      if(TMath::Abs(z_x-73) < 15){ // savd
	sum_up_x = sum_up_x + w*(axi-axj)*(bxi-bxj);
	sum_do_x = sum_do_x + w*(axi-axj)*(axi-axj);
      }

	if(TMath::Abs(z_y-73) < 15){ //savd
	  sum_up_y = sum_up_y + w*(ayi-ayj)*(byi-byj);
	  sum_do_y = sum_do_y + w*(ayi-ayj)*(ayi-ayj);
	}
      //cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
      
    }
  } 

  Float_t zprim_w2 = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);
  //cout<<" zprim_w2 ="<<zprim_w2<<endl;

  //
  Float_t sum_x = 0;
  Float_t sum_y = 0;
  Float_t N = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 
    
    Float_t x = track->GetXatZ(zprim_w2);
    Float_t y = track->GetYatZ(zprim_w2);
    
    sum_x = sum_x + x;      
    sum_y = sum_y + y;
    N = N + 1;      
  }
  
  Float_t xprim = sum_x/N;
  Float_t yprim = sum_y/N;

  sum_x = 0;
  sum_y = 0;
  Float_t norm = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 
    
    Float_t x = track->GetXatZ(zprim_w2);
    Float_t y = track->GetYatZ(zprim_w2);
   
    Float_t d2 = (x-xprim)*(x-xprim) +  (y-yprim)*(y-yprim);
 
    sum_x = sum_x + x/d2;      
    sum_y = sum_y + y/d2;
    norm = norm + 1./d2;      
  }
  
  Float_t xprim_w2 = sum_x/norm;
  Float_t yprim_w2 = sum_y/norm;
  
  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);

}

//______________________________________________________________________________________________________
void Na61AlVdTrackingModule::RemoveGhostTracks(TString devname1,TString devname2, UEventNode* out)
{
  // This method resolve problem of same track found in devname1 and devname2 tables of 3Hit tracks.
  // This might happen due to overlap sensors in Vds4 station 


  UDataTable* tracktab1 = out->GetDataTable(Form("Tracks %s",devname1.Data()));
  UDataTable* tracktab2 = out->GetDataTable(Form("Tracks %s",devname2.Data()));

  //cout<<"---------------> RemoveGhostTracks: "<<devname1.Data()<<"   "<<devname2.Data()<<endl;

  Int_t N1 = tracktab1->GetEntries();
  Int_t N2 = tracktab2->GetEntries();

  //cout<<"tracktab1: "<<tracktab1->GetName()<<"     tracktab2: "<<tracktab2->GetName()<<endl;
  //cout<<"N1="<<N1<<"  N2="<<N2<<endl;

  for(Int_t i=0;i<N1;i++){
    UVdTrack* track1 = (UVdTrack*)tracktab1->At(i);
    for(Int_t j=0;j<N2;j++){
      UVdTrack* track2 = (UVdTrack*)tracktab2->At(j); 

      //cout<<"i="<<i<<"  j="<<j<<endl;

      if(track1->IsMarkedForRemoval() || track2->IsMarkedForRemoval())continue;

      CheckTrackMatching2(track1,track2);	

      //if(track1->IsMarkedForRemoval() || track2->IsMarkedForRemoval()){
      //cout<<" track1: "<<track1->GetX()<<" "<<track1->GetY()<<" "<<track1->GetZ()<<endl;
      //cout<<" track2: "<<track2->GetX()<<" "<<track2->GetY()<<" "<<track2->GetZ()<<endl;
      //}
      
    }
  }
  
  
  for(Int_t i=0;i<N1;i++){
    UVdTrack* track = (UVdTrack*)tracktab1->At(i);
    if(track->IsMarkedForRemoval()){tracktab1->DeleteObjectAt(i);}
  }
  tracktab1->Compress();    
  
  for(Int_t i=0;i<N2;i++){
    UVdTrack* track = (UVdTrack*)tracktab2->At(i);
    if(track->IsMarkedForRemoval()){tracktab2->DeleteObjectAt(i);}
  }
  tracktab2->Compress();    
  
  //Int_t rt = N1 - tracktab1->GetEntries();
  //rt =  rt + N2 - tracktab2->GetEntries();

  //cout<<"removed tracks: "<<rt<<endl;
}

//______________________________________________________________________________________________________
void Na61AlVdTrackingModule::Make3HitTracks(Int_t idev, UEventNode* out, Int_t* TabIndArray)
{

  //cout<<"----------------------------- Make3HitTracks: begin -------------------->"<<endl;
  //out->ListObjects();


  //Int_t idev = FindDevId(devname);
  UDataTable* hits1 = out->GetDataTable(Form("Hits %s",fAlSensorNames[TabIndArray[0]].Data()));
  UDataTable* hits2 = out->GetDataTable(Form("Hits %s",fAlSensorNames[TabIndArray[1]].Data()));
  UDataTable* hits3 = out->GetDataTable(Form("Hits %s",fAlSensorNames[TabIndArray[2]].Data()));

  TH1F* hX_dev = fhX_dev[idev]; 
  TH1F* hY_dev = fhY_dev[idev]; 
  TH1F* hX_dev_cuty = fhX_dev_cuty[idev]; 
  TH1F* hY_dev_cutx = fhY_dev_cutx[idev]; 
  TH2F* h_devx_devy = fh_devx_devy[idev];

  Float_t Offsetx = fParams->GetDevOffsetX(idev);  Float_t Sigmax = fParams->GetDevSigmaX(idev);  
  Float_t Offsety = fParams->GetDevOffsetY(idev);  Float_t Sigmay = fParams->GetDevSigmaY(idev);  

  TString devname = fDevName[idev];

  //cout<<"Tu1"<<" "<<Offsetx<<" "<<Offsety<<" "<<Sigmax<<" "<<Sigmay<<endl;
  //cout<<"idev="<<idev<<"  devname: "<<devname.Data()<<"  sensors indexes: "<<TabIndArray[0]<<" "<<TabIndArray[1]<<" "<<TabIndArray[2]<<endl;
  //cout<<"hits1: "<<hits1<<" tab1 name: "<<fAlSensorNames[TabIndArray[0]]<<" tab2: "<<fAlSensorNames[TabIndArray[1]]<<" tab3: "<<fAlSensorNames[TabIndArray[2]]<<endl;
  //Int_t ii;
  //cin>>ii; 
  //cout<<"hits1="<<hits1->GetEntries()<<"  hits2="<<hits2->GetEntries()<<"  hits3="<<hits3->GetEntries()<<endl;

  //cout<<"Make3Hit(1): "<<hits1[0]<<" "<<hits1[1]<<" "<<hits1[2]<<" "<<hits1[3]<<" "<<hits1[4]<<endl;
  //cout<<"created table for 3hit track: "<<devname.Data()<<endl;
  UDataTable* tracktab = new UDataTable(Form("Tracks %s",devname.Data())); 
  if(tracktab)tracktab->SetOwner(kTRUE);
  
  
  Int_t imatch = 0;
  Int_t IndArray[3];

  //for(Int_t ichip=2; ichip<3; ichip++){ // take only middle chip
  // Select proper hits arrays based on TabIndArray

    //if(!hits1[ichip] || !hits2[ichip] || !hits3[ichip]){
    //cout<<"devname: "<<devname.Data()<<"  ichip="<<ichip<<endl; 
    //continue;
    //}

  for(Int_t i=0;i<hits1->GetEntries();i++){
    USensorHit* hit1 = (USensorHit*)hits1->At(i);      
    if((hit1->GetUniqueID()))continue;
    
    Float_t x1 = hit1->GetX();
    Float_t y1 = hit1->GetY();
    Float_t z1 = hit1->GetZ();
    
    for(Int_t j=0;j<hits2->GetEntries();j++){
      USensorHit* hit2 = (USensorHit*)hits2->At(j);
      if((hit2->GetUniqueID()))continue;
      
      Float_t x2 = hit2->GetX();
      Float_t y2 = hit2->GetY();
      Float_t z2 = hit2->GetZ();
      
      //cout<<" hits3[ichip]->GetEntries(): "<<hits3[ichip]->GetEntries()<<endl;
      for(Int_t k=0;k<hits3->GetEntries();k++){
	USensorHit* hit3 = (USensorHit*)hits3->At(k);
	//cout<<"  hit3: "<<hit3<<endl;
	if((hit3->GetUniqueID()))continue;
	
	Float_t x3 = hit3->GetX();
	Float_t y3 = hit3->GetY();
	Float_t z3 = hit3->GetZ();
	
	
	Float_t devx = ((z2-z1)*x3 + (z3-z2)*x1)/(z3-z1) - x2;	  	  
	Float_t devy = ((z2-z1)*y3 + (z3-z2)*y1)/(z3-z1) - y2;
	
	//cout<<"hit1: "<<x1<<" "<<y1<<" "<<z1<<endl;
	//cout<<"hit2: "<<x2<<" "<<y2<<" "<<z2<<endl;
	//cout<<"hit3: "<<x3<<" "<<y3<<" "<<z3<<endl;

	//cout<<"devx="<<devx<<"  devy="<<devy<<"  fNsig_dev="<<fNsig_dev<<endl;	

	if(hX_dev)hX_dev -> Fill(devx);
	if(hY_dev)hY_dev -> Fill(devy); 
	if(h_devx_devy)h_devx_devy -> Fill(devx,devy);
	
	if(TMath::Abs(devx-Offsetx) < 3*Sigmax)hY_dev_cutx -> Fill(devy);
	if(TMath::Abs(devy-Offsety) < 3*Sigmay)hX_dev_cuty -> Fill(devx);

	Float_t dev2 = (devx-Offsetx)*(devx-Offsetx)/(Sigmax*Sigmax) + (devy-Offsety)*(devy-Offsety)/(Sigmay*Sigmay);
	
	if(dev2>fNsig_dev*fNsig_dev)continue;
	// found track !!!!
	
	imatch++;
	
	IndArray[0]=i; IndArray[1]=j; IndArray[2]=k;
	
	//cout<<" tabIdArrays: "<<TabIndArray[0]<<"  "<<TabIndArray[1]<<" "<<TabIndArray[2]<<endl;
	//cout<<" IndArrays: "<<IndArray[0]<<"  "<<IndArray[1]<<" "<<IndArray[2]<<endl;
	Create3HitTrack(idev, hit1, hit2, hit3, IndArray,TabIndArray,tracktab);
	//cout<<"done with  Create3HitTrack"<<endl;
      }
    }
  }
  
  
  out->AddDataTable(tracktab);
  //cout<<"----------------------------- end -------------------->"<<endl;
  //out->ListObjects();  
  
}

//______________________________________________________________________________________________________
void Na61AlVdTrackingModule::Make3HitTracks(TString devname, Int_t csize, UDataTable* hits1, UDataTable* hits2, UDataTable* hits3,
					    UEventNode* out, Int_t* TabIndArray)
{

  fTrackID=0;

  Int_t idev = GetIdev(devname);

  fLargeClusters=kFALSE;
  TH1F* hX_dev = fhX_dev[idev]; 
  TH1F* hY_dev = fhY_dev[idev]; 
  TH1F* hX_dev_cuty = fhX_dev_cuty[idev]; 
  TH1F* hY_dev_cutx = fhY_dev_cutx[idev]; 
  TH2F* h_devx_devy = fh_devx_devy[idev];

  //cout<<"idev="<<idev<<"  "<<devname.Data()<<"  hX_dev_cuty:"<<hX_dev_cuty<<endl;


  Float_t Offsetx = fParams->GetDevOffsetX(idev);  Float_t Sigmax = fParams->GetDevSigmaX(idev);  
  Float_t Offsety = fParams->GetDevOffsetY(idev);  Float_t Sigmay = fParams->GetDevSigmaY(idev);  

  UDataTable* tracktab = tracktab = new UDataTable(Form("Tracks %s",devname.Data()));
  if(tracktab)tracktab->SetOwner(kTRUE);

  //  
  //Bool_t beamTrack=kFALSE;
  Int_t IndArray[3];

  //cout<<"92 fhFullTracksAll: "<<fhFullTracksAll<<endl;  

  //cout<<"Tu1"<<" "<<Offsetx<<" "<<Offsety<<" "<<Sigmax<<" "<<Sigmay<<endl;
  for(Int_t i=0;i<hits1->GetEntries();i++){
    USensorHit* hit1 = (USensorHit*)hits1->At(i);      
    if((hit1->GetUniqueID()) && (fProductionMode))continue;
    if(hit1->GetClusterSize() < csize)continue;

    //cout<<"i="<<i<<" UniqueID="<<hit1->GetUniqueID()<<" "<<fProductionMode<<endl;

    Float_t x1 = hit1->GetX();
    Float_t y1 = hit1->GetY();
    Float_t z1 = hit1->GetZ();
    
    for(Int_t j=0;j<hits2->GetEntries();j++){
      USensorHit* hit2 = (USensorHit*)hits2->At(j);
      if((hit2->GetUniqueID()) && (fProductionMode))continue;
      if(hit2->GetClusterSize() < csize)continue;
      //cout<<"j="<<j<<endl;
      
      Float_t x2 = hit2->GetX();
      
      //if(!fLargeClusters && (x2>x1))continue; // need more sophisticated condition    
      
      Float_t y2 = hit2->GetY();
      Float_t z2 = hit2->GetZ();
      
     
      Int_t imatch=0;

      for(Int_t k=0;k<hits3->GetEntries();k++){
	USensorHit* hit3 = (USensorHit*)hits3->At(k);
	if((hit3->GetUniqueID()) && (fProductionMode))continue;
	if(hit3->GetClusterSize() < csize)continue;


	Float_t x3 = hit3->GetX();


	Float_t y3 = hit3->GetY();
	Float_t z3 = hit3->GetZ();


	Float_t devx = ((z2-z1)*x3 + (z3-z2)*x1)/(z3-z1) - x2;	  	  
	Float_t devy = ((z2-z1)*y3 + (z3-z2)*y1)/(z3-z1) - y2;


	if(hX_dev)      hX_dev -> Fill(devx);
	if(hY_dev)      hY_dev -> Fill(devy); 
	if(h_devx_devy) h_devx_devy->Fill(devx,devy);

	if(TMath::Abs(devx-Offsetx) < 3*Sigmax)hY_dev_cutx -> Fill(devy);
	if(TMath::Abs(devy-Offsety) < 3*Sigmay)hX_dev_cuty -> Fill(devx);

	Float_t dev2 = (devx-Offsetx)*(devx-Offsetx)/(Sigmax*Sigmax) + (devy-Offsety)*(devy-Offsety)/(Sigmay*Sigmay);

	if(dev2>fNsig_dev*fNsig_dev)continue;
	  // found track !!!!


	imatch++;

	IndArray[0]=i; IndArray[1]=j; IndArray[2]=k;

	//cout<<" tabIdArrays: "<<TabIndArray[0]<<"  "<<TabIndArray[1]<<" "<<TabIndArray[2]<<endl;
	Create3HitTrack(idev,hit1, hit2, hit3, IndArray,TabIndArray,tracktab);

      }
    }
  }

  //cout<<"93 fhFullTracksAll: "<<fhFullTracksAll<<endl;  

  // clean tracks
  Int_t Ns = tracktab->GetEntries();
  for(Int_t i=0;i<Ns;i++){
    UVdTrack* track1 = (UVdTrack*)tracktab->At(i);
    for(Int_t j=i+1;j<Ns;j++){
      UVdTrack* track2 = (UVdTrack*)tracktab->At(j); 

      if(track1->IsMarkedForRemoval() || track2->IsMarkedForRemoval())continue;

      CheckTrackMatching(track1,track2);	
      
    }
  }
  
  //cout<<"94 fhFullTracksAll: "<<fhFullTracksAll<<endl;  

  
  for(Int_t i=0;i<Ns;i++){
    UVdTrack* track = (UVdTrack*)tracktab->At(i);
    if(track->IsMarkedForRemoval()){tracktab->DeleteObjectAt(i);}
  }
  tracktab->Compress();    
  
  Int_t Ne = tracktab->GetEntries();
  if(fhNsVsNe[idev]) fhNsVsNe[idev]->Fill(Ns,Ne);
  if(fhTracks[idev]) fhTracks[idev]->Fill(Ne);
  fAll3StationTracks[idev] += Ne; 


  out->AddDataTable(tracktab);
  
}


/*
//_____________________________________________________________
void Na61AlVdTrackingModule::CreateFullTrack(Int_t imatch,UVdTrack* track1, UVdTrack* track2,UDataTable* tracktab)
{

  USensorHit* hit1 = (USensorHit*)track1->GetHitIdAtStation(0);
  USensorHit* hit2 = (USensorHit*)track2->GetHitIdAtStation(1);
  USensorHit* hit3 = (USensorHit*)track2->GetHitIdAtStation(2);
  USensorHit* hit4 = (USensorHit*)track2->GetHitIdAtStation(3);

  ((USensorHit*)track1->GetHitIdAtStation(1))->SetUniqueID(1);
  ((USensorHit*)track1->GetHitIdAtStation(2))->SetUniqueID(1);

  Int_t ind1 = track1->GetHitIndexOnStation(0);
  Int_t ind2 = track2->GetHitIndexOnStation(1);
  Int_t ind3 = track2->GetHitIndexOnStation(2);
  Int_t ind4 = track2->GetHitIndexOnStation(3);

  Int_t tind1 = track1->GetTabIndexOnStation(0);
  Int_t tind2 = track2->GetTabIndexOnStation(1);
  Int_t tind3 = track2->GetTabIndexOnStation(2);
  Int_t tind4 = track2->GetTabIndexOnStation(3);

  if(track1->GetHitIdAtStation(1) != track2->GetHitIdAtStation(1))
    cout<<" Na61AlVdTrackingModule::CreateFullTrack: hits ID differnt on Vds2 station. In principal should not happened"<<endl;

  TF1* linex = new TF1("linex","[0] +[1]*x",-10.0,160.0);
  TF1* liney = new TF1("liney","[0] +[1]*x",-10.0,160.0);
  
  Float_t x[]={hit1->GetZ(),hit2->GetZ(),hit3->GetZ(),hit4->GetZ()};
  Float_t yx[]={hit1->GetX(),hit2->GetX(),hit3->GetX(),hit4->GetX()};
  Float_t ex[]={0,0,0,0};
  Float_t eyx[]={0.004,0.004,0.004,0.004};
  
  Float_t yy[]={hit1->GetY(),hit2->GetY(),hit3->GetY(),hit4->GetY()};
  Float_t ey[]={0,0,0,0};
  Float_t eyy[]={0.004,0.004,0.004,004};
  
  TGraphErrors* gr_x = new TGraphErrors(4,x,yx,ex,eyx);
  
  TGraphErrors* gr_y = new TGraphErrors(4,x,yy,ey,eyy);
  
  linex->SetParameter(0,yx[0]);

  gr_x->Fit("linex","Q" ,"",-10,160);
  
  Float_t x_Vds1 = linex->Eval(0.);  
  
  Float_t tan_x =  linex->GetParameter(1);

  liney->SetParameter(0,yy[0]);
  gr_y->Fit("liney","Q" ,"",-10,160);

  Float_t y_Vds1 = liney->Eval(0.);  
  
  Float_t tan_y =  liney->GetParameter(1);
  
 
  if(fhAx_full[imatch]) fhAx_full[imatch]->Fill(tan_x);
  if(fhAy_full[imatch]) fhAy_full[imatch]->Fill(tan_y);

  Vector3D origin(x_Vds1,y_Vds1,0.);
  Vector3D direction(tan_x,tan_y,1.);
  UVdTrack* track = new UVdTrack(origin,direction);
  
  tracktab->Add(track);
  
  track->SetTrackID((Long_t)track);
  track->SetPdgID(0);
  
  track->SetVdsHitIDs(hit1,hit2,hit3,hit4); 
  track->SetHitArrayIndexes(ind1,ind2,ind3,ind4);
  track->SetTabArrayIndexes(tind1,tind2,tind3,tind4);
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

  delete  gr_x;
  delete  gr_y;
  delete  linex;
  delete  liney;
  
}
*/

//_____________________________________________________________
void Na61AlVdTrackingModule::CreateFullTrack_new(Int_t imatch,UVdTrack* track1, UVdTrack* track2,UDataTable* tracktab)
{

  //cout<<" Na61AlVdTrackingModule::CreateFullTrack_new: "<<endl;
  //cout<<(Long_t)track1->GetHitIdAtStation(0)<<" "<<(Long_t)track1->GetHitIdAtStation(1)<<" "<<(Long_t)track1->GetHitIdAtStation(2)<<" "<<(Long_t)track1->GetHitIdAtStation(3)<<endl;
  //cout<<(Long_t)track2->GetHitIdAtStation(0)<<" "<<(Long_t)track2->GetHitIdAtStation(1)<<" "<<(Long_t)track2->GetHitIdAtStation(2)<<" "<<(Long_t)track2->GetHitIdAtStation(3)<<endl;

  Int_t VdsInd[2];
  Int_t ind=0;
  for(Int_t ih=0;ih<4;ih++){
    if(track1->GetHitIdAtStation(ih) && track2->GetHitIdAtStation(ih)){VdsInd[ind]=ih; ind++;}
  }


  USensorHit* hit1 = (USensorHit*)track1->GetHitIdAtStation(0);
  USensorHit* hit2 = (USensorHit*)track1->GetHitIdAtStation(1);
  USensorHit* hit3 = (USensorHit*)track1->GetHitIdAtStation(2);
  USensorHit* hit4 = (USensorHit*)track1->GetHitIdAtStation(3);
  Int_t ind1 = track1->GetHitIndexOnStation(0);
  Int_t ind2 = track1->GetHitIndexOnStation(1);
  Int_t ind3 = track1->GetHitIndexOnStation(2);
  Int_t ind4 = track1->GetHitIndexOnStation(3);
  Int_t tind1 = track1->GetTabIndexOnStation(0);
  Int_t tind2 = track1->GetTabIndexOnStation(1);
  Int_t tind3 = track1->GetTabIndexOnStation(2);
  Int_t tind4 = track1->GetTabIndexOnStation(3);

  if(!hit1){
    hit1 = (USensorHit*)track2->GetHitIdAtStation(0); ind1 = track2->GetHitIndexOnStation(0); tind1 = track2->GetTabIndexOnStation(0);
  }
  if(!hit2){
    hit2 = (USensorHit*)track2->GetHitIdAtStation(1); ind2 = track2->GetHitIndexOnStation(1); tind2 = track2->GetTabIndexOnStation(1);
  }
  if(!hit3){
    hit3 = (USensorHit*)track2->GetHitIdAtStation(2); ind3 = track2->GetHitIndexOnStation(2); tind3 = track2->GetTabIndexOnStation(2);
  }
  if(!hit4){
    hit4 = (USensorHit*)track2->GetHitIdAtStation(3); ind4 = track2->GetHitIndexOnStation(3); tind4 = track2->GetTabIndexOnStation(3);
  }


  //((USensorHit*)track1->GetHitIdAtStation(1))->SetUniqueID(1); // don't understand reason of that?
  //((USensorHit*)track1->GetHitIdAtStation(2))->SetUniqueID(1);

  /* to be removed later
  Int_t ind1 = track1->GetHitIndexOnStation(0);
  Int_t ind2 = track2->GetHitIndexOnStation(1);
  Int_t ind3 = track2->GetHitIndexOnStation(2);
  Int_t ind4 = track2->GetHitIndexOnStation(3);

  Int_t tind1 = track1->GetTabIndexOnStation(0);
  Int_t tind2 = track2->GetTabIndexOnStation(1);
  Int_t tind3 = track2->GetTabIndexOnStation(2);
  Int_t tind4 = track2->GetTabIndexOnStation(3);
  */


  if(track1->GetHitIdAtStation(VdsInd[0]) != track2->GetHitIdAtStation(VdsInd[0]))
    cout<<" Na61AlVdTrackingModule::CreateFullTrack_new: (1) hits ID differnt when match on Vds"<<VdsInd[0]+1<<" station. In principal this should not happen"<<endl;


  if(track1->GetHitIdAtStation(VdsInd[1]) != track2->GetHitIdAtStation(VdsInd[1]))
    cout<<" Na61AlVdTrackingModule::CreateFullTrack_new: (2) hits ID differnt when match on Vds"<<VdsInd[1]+1<<" station. In principal this should not happen"<<endl;

  //cout<<"track found in sensors: "<<fAlSensorNames[tind1]<<" "<<fAlSensorNames[tind2]<<" "<<fAlSensorNames[tind3]<<" "<<fAlSensorNames[tind4]<<" "<<endl;

  //Int_t ii;
  //cin>>ii;

  TF1* linex = new TF1("linex","[0] +[1]*x",-10.0,160.0);
  TF1* liney = new TF1("liney","[0] +[1]*x",-10.0,160.0);
  
  double ax;
  double ay;
  double bx;
  double by;
  double chi2x;
  double chi2y;
  double N;
  double zmin;

  double sig[4] = {0.004, 0.0145, 0.032, 0.0416}; // errors in mm
  USensorHit* hits[4] = {hit1,hit2,hit3,hit4};

  for(Int_t i=0;i<4;i++)fhClusterSize->Fill(hits[i]->GetClusterSize());

  FitLine_w2(hits,sig,ax,ay,bx,by,chi2x,chi2y,N,zmin);
  
  //fhChi2Ndf_x->Fill(chi2x/N);
  //fhChi2Ndf_y->Fill(chi2y/N);

  linex->SetParameter(0,bx);
  linex->SetParameter(1,ax);

  liney->SetParameter(0,by);
  liney->SetParameter(1,ay);

  Float_t tan_x =  ax;
  Float_t tan_y =  ay;
 
  if(fhAx_full[imatch]) fhAx_full[imatch]->Fill(tan_x);
  if(fhAy_full[imatch]) fhAy_full[imatch]->Fill(tan_y);
  // this params should go to data base (pawel_22)
  if(fJura  && tan_x > 0.007 && fhAy_full[imatch])fhAy_full_prod[imatch]->Fill(tan_y);
  if(!fJura && tan_x < 0.005 && fhAy_full[imatch])fhAy_full_prod[imatch]->Fill(tan_y);
  
  Vector3D origin(bx,by,0.);
  Vector3D direction(ax,ay,1.);
  UVdTrack* track = new UVdTrack(origin,direction);
  tracktab->Add(track);
  track->SetTrackID((Long_t)track);
  track->SetPdgID(0);
  track->MarkForRemoval(kFALSE);

  // revise this cuts (pawel_22)
  if(tan_x < fParams->GetAxCut() && fJura)track->MarkForRemoval(kTRUE);
  if(tan_x > fParams->GetAxCut() && !fJura)track->MarkForRemoval(kTRUE); // take only third peak in x
  
  //cout<<fParams->GetAxCut()<<endl;

  track->SetVdsHitIDs(hit1,hit2,hit3,hit4);
  track->SetHitArrayIndexes(ind1,ind2,ind3,ind4);
  track->SetTabArrayIndexes(tind1,tind2,tind3,tind4);
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


  //delete  gr_x;
  //delete  gr_y;
  delete  linex;
  delete  liney;
  
}

/*
//_____________________________________________________________
void Na61AlVdTrackingModule::CreateFullTrackx_new(Int_t imatch, UVdTrack* track, UDataTable* tracktab)
{

  USensorHit* hit1 = (USensorHit*)track->GetHitIdAtStation(0);
  USensorHit* hit2 = (USensorHit*)track->GetHitIdAtStation(1);
  USensorHit* hit3 = (USensorHit*)track->GetHitIdAtStation(2);
  USensorHit* hit4 = (USensorHit*)track->GetHitIdAtStation(3);

  Int_t ind1 = track->GetHitIndexOnStation(0);
  Int_t ind2 = track->GetHitIndexOnStation(1);
  Int_t ind3 = track->GetHitIndexOnStation(2);
  Int_t ind4 = track->GetHitIndexOnStation(3);

  Int_t tind1 = track->GetTabIndexOnStation(0);
  Int_t tind2 = track->GetTabIndexOnStation(1);
  Int_t tind3 = track->GetTabIndexOnStation(2);
  Int_t tind4 = track->GetTabIndexOnStation(3);

  TF1* linex = new TF1("linex","[0] +[1]*x",-10.0,160.0);
  TF1* liney = new TF1("liney","[0] +[1]*x",-10.0,160.0);
  
  Float_t ax;
  Float_t ay;
  Float_t bx;
  Float_t by;
  Float_t chi2x;
  Float_t chi2y;
  Float_t N;
  Float_t zmin;

  Float_t sig[4] = {0.004, 0.0145, 0.032, 0.0416}; // errors in mm
  USensorHit* hits[4] = {hit1,hit2,hit3,hit4};

  FitLine_w2(hits,sig,ax,ay,bx,by,chi2x,chi2y,N,zmin);
  
  linex->SetParameter(0,bx);
  linex->SetParameter(1,ax);

  liney->SetParameter(0,by);
  liney->SetParameter(1,ay);
  
  Float_t tan_x =  ax;
  Float_t tan_y =  ay;
  
  if(fhAx_full[imatch]) fhAx_full[imatch]->Fill(tan_x);
  if(fhAy_full[imatch]) fhAy_full[imatch]->Fill(tan_y);
  if(fJura  && tan_x > 0.007 && fhAy_full[imatch])fhAy_full_prod[imatch]->Fill(tan_y);
  if(!fJura && tan_x < 0.005 && fhAy_full[imatch])fhAy_full_prod[imatch]->Fill(tan_y);

 
  Vector3D origin(bx,by,0.);
  Vector3D direction(ax,ay,1.);
  UVdTrack* trackx = new UVdTrack(origin,direction);
  
  tracktab->Add(trackx);
  
  trackx->SetTrackID((Long_t)trackx);
  trackx->SetPdgID(0);
  trackx->MarkForRemoval(kFALSE);

  if(tan_x < 0.007 && fJura)trackx->MarkForRemoval(kTRUE);
  if(tan_x > 0.005 && !fJura)trackx->MarkForRemoval(kTRUE);


  trackx->SetVdsHitIDs(hit1,hit2,hit3,hit4);
  trackx->SetHitArrayIndexes(ind1,ind2,ind3,ind4);
  trackx->SetTabArrayIndexes(tind1,tind2,tind3,tind4);
  if(hit1)hit1->SetUniqueID(1);
  if(hit2)hit2->SetUniqueID(1);
  if(hit3)hit3->SetUniqueID(1);
  if(hit4)hit4->SetUniqueID(1);
  
  // make deviation diagnostics
  if(hit1)fhX_dev_full[imatch]->Fill(hit1->GetX() - linex->Eval(hit1->GetZ()));
  if(hit2)fhX_dev_full[imatch]->Fill(hit2->GetX() - linex->Eval(hit2->GetZ()));
  if(hit3)fhX_dev_full[imatch]->Fill(hit3->GetX() - linex->Eval(hit3->GetZ()));
  if(hit4)fhX_dev_full[imatch]->Fill(hit4->GetX() - linex->Eval(hit4->GetZ()));

  if(hit1)fhX_dev_Vds1[imatch]->Fill(hit1->GetX() - linex->Eval(hit1->GetZ()));
  if(hit2)fhX_dev_Vds2[imatch]->Fill(hit2->GetX() - linex->Eval(hit2->GetZ()));
  if(hit3)fhX_dev_Vds3[imatch]->Fill(hit3->GetX() - linex->Eval(hit3->GetZ()));
  if(hit4)fhX_dev_Vds4[imatch]->Fill(hit4->GetX() - linex->Eval(hit4->GetZ()));

  if(hit1)fhY_dev_full[imatch]->Fill(hit1->GetY() - liney->Eval(hit1->GetZ()));
  if(hit2)fhY_dev_full[imatch]->Fill(hit2->GetY() - liney->Eval(hit2->GetZ()));
  if(hit3)fhY_dev_full[imatch]->Fill(hit3->GetY() - liney->Eval(hit3->GetZ()));
  if(hit4)fhY_dev_full[imatch]->Fill(hit4->GetY() - liney->Eval(hit4->GetZ()));

  if(hit1)fhY_dev_Vds1[imatch]->Fill(hit1->GetY() - liney->Eval(hit1->GetZ()));
  if(hit2)fhY_dev_Vds2[imatch]->Fill(hit2->GetY() - liney->Eval(hit2->GetZ()));
  if(hit3)fhY_dev_Vds3[imatch]->Fill(hit3->GetY() - liney->Eval(hit3->GetZ()));
  if(hit4)fhY_dev_Vds4[imatch]->Fill(hit4->GetY() - liney->Eval(hit4->GetZ()));

  delete  linex;
  delete  liney;

  // calculate deviations for accepted 3hit tracks
  if(!hit1)FillDeviations(imatch,hit2,hit3,hit4);
  if(!hit2)FillDeviations(imatch,hit1,hit3,hit4);
  if(!hit3)FillDeviations(imatch,hit1,hit2,hit4);
  if(!hit4)FillDeviations(imatch,hit1,hit2,hit3);
}
*/

//_____________________________________________________________________________________________
void Na61AlVdTrackingModule::FillDeviations(Int_t imatch, USensorHit* hit1, USensorHit* hit2, USensorHit* hit3)
{
    Float_t x1 = hit1->GetX();
    Float_t y1 = hit1->GetY();
    Float_t z1 = hit1->GetZ();

    Float_t x2 = hit2->GetX();
    Float_t y2 = hit2->GetY();
    Float_t z2 = hit2->GetZ();

    Float_t x3 = hit3->GetX();
    Float_t y3 = hit3->GetY();
    Float_t z3 = hit3->GetZ();
    
    
    Float_t devx = ((z2-z1)*x3 + (z3-z2)*x1)/(z3-z1) - x2;	  	  
    Float_t devy = ((z2-z1)*y3 + (z3-z2)*y1)/(z3-z1) - y2;
    
    
    if(fhX_dev_acce[imatch])fhX_dev_acce[imatch] -> Fill(devx);
    if(fhY_dev_acce[imatch])fhY_dev_acce[imatch] -> Fill(devy); 
  
}

//___________________________________________________________________________________
//void Na61AlVdTrackingModule::FitLine_w2(USensorHit** hits, Float_t* sig,
//				      Float_t& ax,Float_t& ay,Float_t& bx,Float_t& by,
//				      Float_t& chi2x,Float_t& chi2y,Float_t& N, Float_t& zmin)
void Na61AlVdTrackingModule::FitLine_w2(USensorHit** hits, double* sig,
				      double& ax,double& ay,double& bx,double& by,
				      double& chi2x,double& chi2y,double& N, double& zmin)
{
  // simple regration (assuming constant errors)
  // It is quite straight forward to use weighted regression

    Float_t hx,hy,hz;

    Float_t S=0;
    Float_t Sz=0;
    Float_t Szz=0;
 
    Float_t Sx=0;
    Float_t Szx=0;

    Float_t Sy=0;
    Float_t Szy=0;
 
        
    //cout<<hit1->GetZ()<<" "<<hit2->GetZ()<<" "<<hit3->GetZ()<<" "<<hit4->GetZ()<<endl;

    zmin = 10000.;
    N = 0;    
    
    // perhaps we need waited regression
    
    for(Int_t i=0;i<4;i++){      
      USensorHit* hit = hits[i];
      
      if(!hit)continue; // now we always have 4 hits
      
      N++;
      
      // smearing is done in dedicated method
      hx = hit->GetX(); 
      hy = hit->GetY();
      hz = hit->GetZ(); 
      
      //Int_t   ivds = hit->GetStationID()-1;
      Float_t sig2 = sig[i]*sig[i];     

      S   += 1./sig2;
      Sz  += hz/sig2;
      Szz += (hz*hz)/sig2;

      Sx  += hx/sig2;
      Szx += (hz*hx)/sig2;

      Sy  += hy/sig2;
      Szy += (hz*hy)/sig2;
      
      if(hz<zmin)zmin=hz;
      
    }
    
    ax = (S*Szx - Sz*Sx)/(S*Szz-Sz*Sz);
    bx = (Szz*Sx - Sz*Szx)/(S*Szz-Sz*Sz);

    ay = (S*Szy - Sz*Sy)/(S*Szz-Sz*Sz);
    by = (Szz*Sy - Sz*Szy)/(S*Szz-Sz*Sz);
    

    chi2x=0;
    chi2y=0;
    for(Int_t i=0;i<4;i++){
      USensorHit* hit = hits[i];
      
      if(!hit)continue;      

      //Int_t   ivds = hit->GetStationID();
      Float_t sig2 = sig[i]*sig[i];     

      hx = hit->GetX(); 
      hy = hit->GetY();
      hz = hit->GetZ(); 

      chi2x += (hx - ax*hz - bx)*(hx - ax*hz - bx)/sig2;
      chi2y += (hy - ay*hz - by)*(hy - ay*hz - by)/sig2;
    }


}

/* Method used for SAVD detector
//____________________________________________________________________________________
void Na61AlVdTrackingModule::Create3HitTrack(Int_t idev,USensorHit* hit1, USensorHit* hit2, USensorHit* hit3, Int_t* IndArray, Int_t* TabIndArray, UDataTable* tracktab)
{

  TF1* line = new TF1("line","[0] +[1]*x",-10.0,160.0);
  
  Float_t x[]={hit1->GetZ(),hit2->GetZ(),hit3->GetZ()};
  Float_t yx[]={hit1->GetX(),hit2->GetX(),hit3->GetX()};
  Float_t ex[]={0,0,0};
  Float_t eyx[]={0.004,0.004,0.004};
  
  Float_t yy[]={hit1->GetY(),hit2->GetY(),hit3->GetY()};
  Float_t eyy[]={0.004,0.004,0.004};
  
  TGraphErrors* gr_x= new TGraphErrors(3,x,yx,ex,eyx);
  
  TGraphErrors* gr_y= new TGraphErrors(3,x,yy,ex,eyy);
  
  line->SetParameter(0,yx[0]);
  gr_x->Fit("line","Q" ,"",x[0]-1,x[2]+1);

  if(fhChi2X[idev])fhChi2X[idev]->Fill(line->GetChisquare());
  
  Float_t x_Vds1 = line->Eval(0.);  
  
  Float_t tan_x =  line->GetParameter(1);

  line->SetParameter(0,yy[0]);
  gr_y->Fit("line","Q" ,"",x[0]-1,x[2]+1);

  if(fhChi2Y[idev])fhChi2Y[idev]->Fill(line->GetChisquare());

  Float_t y_Vds1 = line->Eval(0.);  
  
  Float_t tan_y =  line->GetParameter(1);
  

  if(fhAx[idev]) fhAx[idev]->Fill(tan_x);
  if(fhAy[idev]) fhAy[idev]->Fill(tan_y);

  Vector3D origin(x_Vds1,y_Vds1,0.);
  Vector3D direction(tan_x,tan_y,1.);
  UVdTrack* track = new UVdTrack(origin,direction);
  
  tracktab->Add(track);
  
  track->SetTrackID(fTrackID++);
  track->SetPdgID(0);

  //cout<<" create tabIdArrays: "<<TabIndArray[0]<<"  "<<TabIndArray[1]<<" "<<TabIndArray[2]<<endl;

  if(idev<6){
    if(fDevNames[idev].Contains("dev1_")){
      track->SetVdsHitIDs(hit1,hit2,hit3,0);
      track->SetHitArrayIndexes(IndArray[0],IndArray[1],IndArray[2],-1);
      track->SetTabArrayIndexes(TabIndArray[0],TabIndArray[1],TabIndArray[2],-1);
    }
    if(fDevNames[idev].Contains("dev2_")){
      track->SetVdsHitIDs(0,hit1,hit2,hit3);
      track->SetHitArrayIndexes(-1,IndArray[0],IndArray[1],IndArray[2]);
      track->SetTabArrayIndexes(-1,TabIndArray[0],TabIndArray[1],TabIndArray[2]);
    }
  }else{
    if(fDevNames[idev].Contains("x1")){
      track->SetVdsHitIDs(0,hit1,hit2,hit3);
      track->SetHitArrayIndexes(-1,IndArray[0],IndArray[1],IndArray[2]);
      track->SetTabArrayIndexes(-1,TabIndArray[0],TabIndArray[1],TabIndArray[2]);
    }
    if(fDevNames[idev].Contains("x2")){
      track->SetVdsHitIDs(hit1,0,hit2,hit3);
      track->SetHitArrayIndexes(IndArray[0],-1,IndArray[1],IndArray[2]);
      track->SetTabArrayIndexes(TabIndArray[0],-1,TabIndArray[1],TabIndArray[2]);
    }
    if(fDevNames[idev].Contains("x3")){
      track->SetVdsHitIDs(hit1,hit2,0,hit3);
      track->SetHitArrayIndexes(IndArray[0],IndArray[1],-1,IndArray[2]);
      track->SetTabArrayIndexes(TabIndArray[0],TabIndArray[1],-1,TabIndArray[2]);
    }
    if(fDevNames[idev].Contains("x4")){
      track->SetVdsHitIDs(hit1,hit2,hit3,0);
      track->SetHitArrayIndexes(IndArray[0],IndArray[1],IndArray[2],-1);
      track->SetTabArrayIndexes(TabIndArray[0],TabIndArray[1],TabIndArray[2],-1);
    }
  }
  //cout<<fDevNames[idev].Data()<<" "<<track->GetHitIdAtStation(0)<<" "<<track->GetHitIdAtStation(1)<<" "
  //  <<track->GetHitIdAtStation(2)<<" "<<track->GetHitIdAtStation(3)<<" "<<hit1->GetZ()<<" "<<hit2->GetZ()<<" "<<hit3->GetZ()<<endl;

  delete  gr_x;
  delete  gr_y;
  delete  line;
  
}
*/

//____________________________________________________________________________________
void Na61AlVdTrackingModule::Create3HitTrack(Int_t idev, USensorHit* hit1, USensorHit* hit2, USensorHit* hit3, Int_t* IndArray, Int_t* TabIndArray, UDataTable* tracktab)
{

  TF1* line = new TF1("line","[0] +[1]*x",-10.0,160.0);
  
  Float_t x[]={hit1->GetZ(),hit2->GetZ(),hit3->GetZ()};
  Float_t yx[]={hit1->GetX(),hit2->GetX(),hit3->GetX()};
  Float_t ex[]={0,0,0};
  Float_t eyx[]={0.004,0.004,0.004};
  
  Float_t yy[]={hit1->GetY(),hit2->GetY(),hit3->GetY()};
  Float_t eyy[]={0.004,0.004,0.004};
  
  TGraphErrors* gr_x= new TGraphErrors(3,x,yx,ex,eyx);
  
  TGraphErrors* gr_y= new TGraphErrors(3,x,yy,ex,eyy);
  
  line->SetParameter(0,yx[0]);
  gr_x->Fit("line","Q" ,"",x[0]-1,x[2]+1);

  if(fhChi2X[idev])fhChi2X[idev]->Fill(line->GetChisquare());
  
  Float_t x_Vds1 = line->Eval(0.);  
  
  Float_t tan_x =  line->GetParameter(1);

  line->SetParameter(0,yy[0]);
  gr_y->Fit("line","Q" ,"",x[0]-1,x[2]+1);

  if(fhChi2Y[idev])fhChi2Y[idev]->Fill(line->GetChisquare());

  Float_t y_Vds1 = line->Eval(0.);  
  
  Float_t tan_y =  line->GetParameter(1);
  
  if(fhAx[idev]) fhAx[idev]->Fill(tan_x);
  if(fhAy[idev]) fhAy[idev]->Fill(tan_y);

  Vector3D origin(x_Vds1,y_Vds1,0.);
  Vector3D direction(tan_x,tan_y,1.);
  UVdTrack* track = new UVdTrack(origin,direction);
  
  tracktab->Add(track);
  
  track->SetTrackID(fTrackID++);
  track->SetPdgID(0);

  //cout<<" create tabIdArrays: "<<TabIndArray[0]<<"  "<<TabIndArray[1]<<" "<<TabIndArray[2]<<endl;
  TString devname = fDevName[idev];

  if(devname.Contains("A")){
    track->SetVdsHitIDs(hit1,hit2,hit3,0);
    track->SetHitArrayIndexes(IndArray[0],IndArray[1],IndArray[2],-1);
    track->SetTabArrayIndexes(TabIndArray[0],TabIndArray[1],TabIndArray[2],-1);
  }
  if(devname.Contains("B")){
    track->SetVdsHitIDs(hit1,0,hit2,hit3);
    track->SetHitArrayIndexes(IndArray[0],-1,IndArray[1],IndArray[2]);
    track->SetTabArrayIndexes(TabIndArray[0],-1,TabIndArray[1],TabIndArray[2]);
  }
  if(devname.Contains("C")){
    track->SetVdsHitIDs(hit1,hit2,0,hit3);
    track->SetHitArrayIndexes(IndArray[0],IndArray[1],-1,IndArray[2]);
    track->SetTabArrayIndexes(TabIndArray[0],TabIndArray[1],-1,TabIndArray[2]);
  }
  if(devname.Contains("D")){
    track->SetVdsHitIDs(0,hit1,hit2,hit3);
    track->SetHitArrayIndexes(-1,IndArray[0],IndArray[1],IndArray[2]);
    track->SetTabArrayIndexes(-1,TabIndArray[0],TabIndArray[1],TabIndArray[2]);
  }

  //cout<<fDevNames[idev].Data()<<" "<<track->GetHitIdAtStation(0)<<" "<<track->GetHitIdAtStation(1)<<" "
  //  <<track->GetHitIdAtStation(2)<<" "<<track->GetHitIdAtStation(3)<<" "<<hit1->GetZ()<<" "<<hit2->GetZ()<<" "<<hit3->GetZ()<<endl;
  
  delete  gr_x;
  delete  gr_y;
  delete  line;
  
}

//_____________________________________________________________
void Na61AlVdTrackingModule::CheckTrackMatching(UVdTrack* track1, UVdTrack* track2)
{
  // TWO or more same hits in the track leads to discarding one track
  
  Int_t isame=0;
  for(Int_t i=0;i<4;i++){
    if((track1->GetHitIdAtStation(i)==track2->GetHitIdAtStation(i)) && track1->GetHitIdAtStation(i) )isame++;
  }

  if(isame<2) return;

  //Float_t ang1_x = track1->GetDX()/track1->GetDZ();
  //Float_t ang1_y = track1->GetDY()/track1->GetDZ();
  //Float_t ang2_x = track2->GetDX()/track2->GetDZ();
  //Float_t ang2_y = track2->GetDY()/track2->GetDZ();

 
  // make one track for removal 
  if( track1->GetChi2Ndf() > track2->GetChi2Ndf() ){
    //if(track1->IsMarkedForRemoval())cout<<"track1 already marked"<<endl;
    track1->MarkForRemoval(kTRUE); //cout<<"marked track1: "<<track1<<endl;
  }else{
    //if(track2->IsMarkedForRemoval())cout<<"track2 already marked"<<endl;
    track2->MarkForRemoval(kTRUE); //cout<<"marked track2: "<<track2<<endl;
  }
  

}


//_____________________________________________________________
void Na61AlVdTrackingModule::CheckTrackMatching2(UVdTrack* track1, UVdTrack* track2)
{
  // TWO or more same hits in the track leads to discarding one track
  
  Int_t isame=0;
  for(Int_t i=0;i<4;i++){
    if((track1->GetHitIdAtStation(i)==track2->GetHitIdAtStation(i)) && track1->GetHitIdAtStation(i) )isame++;
    //cout<<"i = "<<i<<"   "<<track1->GetHitIdAtStation(i)<<"   "<<track2->GetHitIdAtStation(i)<<endl;
  }
  
  if(isame<2) return;


  //Float_t ang1_x = track1->GetDX()/track1->GetDZ();
  //Float_t ang1_y = track1->GetDY()/track1->GetDZ();
  //Float_t ang2_x = track2->GetDX()/track2->GetDZ();
  //Float_t ang2_y = track2->GetDY()/track2->GetDZ();

 
  // make one track for removal 
  if( track1->GetChi2Ndf() > track2->GetChi2Ndf() ){
    //if(track1->IsMarkedForRemoval())cout<<"track1 already marked"<<endl;
    track1->MarkForRemoval(kTRUE); //cout<<"marked track1: "<<track1<<endl;
  }else{
    //if(track2->IsMarkedForRemoval())cout<<"track2 already marked"<<endl;
    track2->MarkForRemoval(kTRUE); //cout<<"marked track2: "<<track2<<endl;
  }
  

}



//_____________________________________________________________
void Na61AlVdTrackingModule::End()
{
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61AlVdTrackingModule::Finish()
{
  // Job-level finalisation
  SetState(kFinish);
  
}



//____________________________________________________________________
void Na61AlVdTrackingModule ::Print(Option_t* option) const
{
  // Print module informationd
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



