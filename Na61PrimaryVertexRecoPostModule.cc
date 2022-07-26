//--------------------------------------------
// Primary vertex reco module
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61PrimaryVertexRecoPostModule
#include "Na61PrimaryVertexRecoPostModule.h"    
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
//ClassImp(Na61PrimaryVertexRecoPostModule);

//____________________________________________________________________
Na61PrimaryVertexRecoPostModule::Na61PrimaryVertexRecoPostModule()
{
  // Default constructor. DO NOT USE
  SetState(kSetup);
  fZprim = 47.7;  
}

//____________________________________________________________________
  Na61PrimaryVertexRecoPostModule::Na61PrimaryVertexRecoPostModule(const Char_t* name, const Char_t* title)
: Na61PrimaryVertexRecoModule(name, title)
{
  // Named Constructor
  SetState(kSetup);
  fZprim = 47.7;  

}

//____________________________________________________________________
void Na61PrimaryVertexRecoPostModule::DefineHistograms()
{

  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init"); 
    return;  

  }

  TDirectory* histDir = gDirectory->mkdir("PrimaryVertexRecoPost");

  histDir->cd(); 

  fhVtxX = new TH1F("hVtxX","", 30,-3.,3);
  fhVtxY = new TH1F("hVtxY","", 30,-3.,3);
  fhTracksVsVtxX_J = new TH1F("hTracksVsVtxX_J","", 30,-3.,3);
  fhTracksVsVtxY_J = new TH1F("hTracksVsVtxY_J","", 30,-3.,3);
  fhTracksVsVtxX_S = new TH1F("hTracksVsVtxX_S","", 30,-3.,3);
  fhTracksVsVtxY_S = new TH1F("hTracksVsVtxY_S","", 30,-3.,3);

  fhVtxDx = new TH1F("hVtxDx","", 1000,-2.5,2.5);
  fhVtxDy = new TH1F("hVtxDy","", 1000,-2.5,2.5);
  fhVtxDz = new TH1F("hVtxDz","", 1000,-10.,10.);

  // using Jura versus Saleve arm
  fhVtxDx_Glob = new TH1F("hVtxDx_Glob","", 1000,-2.5,2.5);
  fhVtxDy_Glob = new TH1F("hVtxDy_Glob","", 1000,-2.5,2.5);
  fhVtxDz_Glob = new TH1F("hVtxDz_Glob","", 1000,-10.,10.);

  fhVtxDx_GlobAcce = new TH1F("hVtxDx_GlobAcce","", 1000,-2.5,2.5);
  fhVtxDy_GlobAcce = new TH1F("hVtxDy_GlobAcce","", 1000,-2.5,2.5);
  fhVtxDz_GlobAcce = new TH1F("hVtxDz_GlobAcce","", 1000,-10.,10.);

  // using Sub1 versus Sub2
  fhVtxDx_Glob2 = new TH1F("hVtxDx_Glob2","", 1000,-1.5,1.5);
  fhVtxDy_Glob2 = new TH1F("hVtxDy_Glob2","", 1000,-1.0,1.0);
  fhVtxDz_Glob2 = new TH1F("hVtxDz_Glob2","", 1000,-10.,10.);
  fhVtxDx_Glob2Acce = new TH1F("hVtxDx_Glob2Acce","", 1000,-1.5,1.5);
  fhVtxDy_Glob2Acce = new TH1F("hVtxDy_Glob2Acce","", 1000,-1.0,1.0);
  fhVtxDz_Glob2Acce = new TH1F("hVtxDz_Glob2Acce","", 1000,-10.,10.);

  fhVtxStatus_JuraVsSaleve = new TH2F("VtxStatus_JuraVsSaleve","", 5,0,5, 5,0,5);
  fhFullTracks_JvsS = new TH2F("hFullTracks_JvsS","", 200,0,500, 200,0,500);
  fhFullTracks_JS = new TH1F("hFullTracks_JS","", 1000,0,1000);
  fhFullTracks_J = new TH1F("hFullTracks_J","", 500,0,500);
  fhFullTracks_S = new TH1F("hFullTracks_S","", 500,0,500);
  fhFullTracks4h_JS = new TH1F("hFullTracks4h_JS","", 300,0,300);
  fhFullTracks4h_J = new TH1F("hFullTracks4h_J","", 300,0,300);
  fhFullTracks4h_S = new TH1F("hFullTracks4h_S","", 300,0,300);

  TString str[] = {"","_J","_S"};  

  for(Int_t i=0;i<3;i++){
    fhRecoVertexZ[i]  = new  TH1F(Form("hRecoVertexZ%s",str[i].Data())," ", 2000,-500.,100.);
    fhRecoVertexZ_fine[i] = new  TH1F(Form("hRecoVertexZ_fine%s",str[i].Data())," ", 1000,-70.,-30.);
    fhRecoVertexXY[i] = new TH2F(Form("hRecoVertexXY%s",str[i].Data())," ", 500,-50.,50.,500,-50.,50.);

    fhRecoVertexZ_flag1[i]  = new  TH1F(Form("hRecoVertexZ_flag1%s",str[i].Data())," ", 2000,-2000.,100.);
    fhRecoVertexZ_fine_flag1[i] = new  TH1F(Form("hRecoVertexZ_fine_flag1%s",str[i].Data())," ", 1000,-70.,-30.);
    fhRecoVertexXY_flag1[i] = new TH2F(Form("hRecoVertexXY_flag1%s",str[i].Data())," ", 500,-50.,50.,500,-50.,50.);
  }   


  fhAxJ = new TH1F(Form("hAxJ"),"",1500,-0.25,0.25);
  fhAyJ = new TH1F(Form("hAyJ"),"",1500,-0.25,0.25);
  fhAxS = new TH1F(Form("hAxS"),"",1500,-0.25,0.25);
  fhAyS = new TH1F(Form("hAyS"),"",1500,-0.25,0.25);

  fhAxJ_flag02 = new TH1F(Form("hAxJ_flag02"),"",1500,-0.25,0.25);
  fhAyJ_flag02 = new TH1F(Form("hAyJ_flag02"),"",1500,-0.25,0.25);

  fhAxJ_flag0 = new TH1F(Form("hAxJ_flag0"),"",1500,-0.25,0.25);
  fhAyJ_flag0 = new TH1F(Form("hAyJ_flag0"),"",1500,-0.25,0.25);
  fhAxJ_flag1 = new TH1F(Form("hAxJ_flag1"),"",1500,-0.25,0.25);
  fhAxJ_flag2 = new TH1F(Form("hAxJ_flag2"),"",1500,-0.25,0.25);
  fhAyJ_flag2 = new TH1F(Form("hAyJ_flag2"),"",1500,-0.25,0.25);

  fhAxS_flag02 = new TH1F(Form("hAxS_flag02"),"",1500,-0.25,0.25);
  fhAyS_flag02 = new TH1F(Form("hAyS_flag02"),"",1500,-0.25,0.25);

  fhAxS_flag0 = new TH1F(Form("hAxS_flag0"),"",1500,-0.25,0.25);
  fhAyS_flag0 = new TH1F(Form("hAyS_flag0"),"",1500,-0.25,0.25);
  fhAxS_flag1 = new TH1F(Form("hAxS_flag1"),"",1500,-0.25,0.25);
  fhAxS_flag2 = new TH1F(Form("hAxS_flag2"),"",1500,-0.25,0.25);
  fhAyS_flag2 = new TH1F(Form("hAyS_flag2"),"",1500,-0.25,0.25);

  gDirectory->cd("..");  
}

//____________________________________________________________________
void Na61PrimaryVertexRecoPostModule::Init()
{
  // Job-level initialisation
  SetState(kInit);
  Na61PrimaryVertexRecoModule::Init();


}

//____________________________________________________________________
void Na61PrimaryVertexRecoPostModule::Begin()
{
  // Run-level initialisation
  SetState(kBegin);
  Na61PrimaryVertexRecoModule::Begin();

}

//____________________________________________________________________
void Na61PrimaryVertexRecoPostModule :: Event(UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode)

{

  // Per event method
  SetState(kEvent);  

  fAllTracksJ = 0;
  fAllTracksS = 0;

  Float_t pvtxX_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexX();
  Float_t pvtxY_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexY();
  Float_t pvtxZ_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexZ();

  Float_t pvtxX_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexX();
  Float_t pvtxY_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexY();
  Float_t pvtxZ_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexZ();

  Int_t pvertStatusJ = ((UVdEvent*)inNodeJ)->GetPrimaryVertexStatus();
  Int_t pvertStatusS = ((UVdEvent*)inNodeS)->GetPrimaryVertexStatus();

  fhVtxStatus_JuraVsSaleve->Fill(pvertStatusJ,pvertStatusS);

  //cout<<" pvertStatusJ="<<pvertStatusJ<<" pvertStatusS="<<pvertStatusS<<endl;

  //cout<<"###############################   new event: pvtxZ_J="<<pvtxZ_J<<"  pvtxZ_S="<<pvtxZ_S<<endl;

  //if((pvertStatusJ>1) && (pvertStatusS>1)){
  if((pvertStatusJ==2) && (pvertStatusS==2)){
    fhVtxDx->Fill(pvtxX_J - pvtxX_S);
    fhVtxDy->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz->Fill(pvtxZ_J - pvtxZ_S);
  }

  fGoodVertex = kFALSE;

  MethodComb(inNodeJ,inNodeS);


  pvtxX_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexX();
  pvtxY_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexY();
  pvtxZ_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexZ();

  pvtxX_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexX();
  pvtxY_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexY();
  pvtxZ_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexZ();


  if(pvertStatusJ == 2 && pvertStatusS == 2){
    fhVtxDx_Glob->Fill(pvtxX_J - pvtxX_S);
    fhVtxDy_Glob->Fill(pvtxY_J - pvtxY_S);
    fhVtxDz_Glob->Fill(pvtxZ_J - pvtxZ_S);

    if(fGoodVertex){
      fhVtxDx_GlobAcce->Fill(pvtxX_J - pvtxX_S);
      fhVtxDy_GlobAcce->Fill(pvtxY_J - pvtxY_S);
      fhVtxDz_GlobAcce->Fill(pvtxZ_J - pvtxZ_S);
    }
  }


  ((UVdEvent*)outNode)->SetPrimaryVertexStatus(0);

  if(pvertStatusJ || pvertStatusS){
    ((UVdEvent*)outNode)->SetPrimaryVertexStatus(1);
    ((UVdEvent*)outNode)->SetPrimaryVertex(fPrimaryVertex.GetX(),fPrimaryVertex.GetY(),fPrimaryVertex.GetZ());
    //cout<<" all tracks: "<<fAllTracksJ+fAllTracksS<<",  fPrimaryVertexZ="<<fPrimaryVertex.GetZ()<<endl;
    //cout<<"fPrimaryVertexX="<<fPrimaryVertex.GetX()<<"  fPrimaryVertexY="<<fPrimaryVertex.GetY()<<endl;
    // to avoid using HT for exotic vertex locations
    if( !(TMath::Abs(fPrimaryVertex.GetZ() - (-49.)) < 10.))((UVdEvent*)outNode)->SetPrimaryVertexStatus(0);
  }


}

//_____________________________________________________________
void Na61PrimaryVertexRecoPostModule::MethodComb(UEventNode* inNodeJ, UEventNode* inNodeS)
{

  ClassifyTracks(inNodeJ,inNodeS);
  //ClassifyTracks_pPb(inNodeJ,inNodeS);

  FindPrimaryVertexPost(1,inNodeJ,inNodeS);

  FindPrimaryVertexPost(2,inNodeJ,0);
  FindPrimaryVertexPost(2,0,inNodeS);

  Float_t pvtxZ_J = ((UVdEvent*)inNodeJ)->GetPrimaryVertexZ();
  Float_t pvtxZ_S = ((UVdEvent*)inNodeS)->GetPrimaryVertexZ();


  // Jura and Saleve vertexes are consistent - same iteraction observed
  if(TMath::Abs(pvtxZ_J - pvtxZ_S) < 3*fParams->GetVtxDzSigma())fGoodVertex = kTRUE;

  //cout<<pvtxZ_J<<" "<<pvtxZ_S<<" "<<fParams->GetVtxDzSigma()<<" "<<fGoodVertex <<endl;

  FindPrimaryVertexPost(2,inNodeJ,inNodeS);


  CheckPrimaryVertexResolution(2,inNodeJ,inNodeS);

}


//_____________________________________________________________
void Na61PrimaryVertexRecoPostModule::ClassifyTracks(UEventNode* outJ, UEventNode* outS)
{
  // clasification regarding x-slope 

  for(Int_t i=0;i<18;i++){
    UDataTable* tracktab = outJ -> GetDataTable(Form("FullTracks %s",fMatchStr[i].Data()));
    if(!tracktab)continue;

    for(Int_t j=0; j<tracktab->GetEntries(); j++){
      UVdTrack* track = ( UVdTrack*)tracktab->At(j);
      Float_t ax = track->GetDX()/track->GetDZ();
      Float_t ay = track->GetDY()/track->GetDZ();

      if(i>1){
        track->SetFlag(2);
        fhAxJ_flag2->Fill(ax);
        fhAyJ_flag2->Fill(ay);
      }else{
        if(ax<(0.006+fRotY_J)){
          track->SetFlag(0);
          fhAxJ_flag0->Fill(ax);
          fhAyJ_flag0->Fill(ay);
        }else if(ax<(fOffAx_J+fRotY_J)){
          track->SetFlag(1);
          fhAxJ_flag1->Fill(ax);
        }else{
          track->SetFlag(2);
          fhAxJ_flag2->Fill(ax);
          fhAyJ_flag2->Fill(ay);
        }
      }

      if(track->GetFlag() != 2)continue;
      fhAxJ->Fill(ax);
      fhAyJ->Fill(ay);

    }
  }

  for(Int_t i=0;i<18;i++){
    UDataTable* tracktab = outS -> GetDataTable(Form("FullTracks %s",fMatchStr[i].Data()));
    if(!tracktab)continue;

    for(Int_t j=0; j<tracktab->GetEntries(); j++){
      UVdTrack* track = ( UVdTrack*)tracktab->At(j);

      Float_t ax = track->GetDX()/track->GetDZ();
      Float_t ay = track->GetDY()/track->GetDZ();

      if(i>1){
        track->SetFlag(2);
        fhAxS_flag2->Fill(ax);
        fhAyS_flag2->Fill(ay);
      }else{
        if(ax>(0.005+fRotY_S)){
          track->SetFlag(0);
          fhAxS_flag0->Fill(ax);
          fhAyS_flag0->Fill(ay);
        }else if(ax>(fOffAx_S+fRotY_S)){
          track->SetFlag(1);
          fhAxS_flag1->Fill(ax);
        }else{
          track->SetFlag(2);
          fhAxS_flag2->Fill(ax);
          fhAyS_flag2->Fill(ay);
        }
      }

      if(track->GetFlag() != 2)continue;
      fhAxS->Fill(ax);
      fhAyS->Fill(ay);
    }
  }


}

//_____________________________________________________________
void Na61PrimaryVertexRecoPostModule::ClassifyTracks_pPb(UEventNode* outJ, UEventNode* outS)
{
  // clasification regarding x-slope and y-slope 

  for(Int_t i=0;i<18;i++){
    UDataTable* tracktab = outJ -> GetDataTable(Form("FullTracks %s",fMatchStr[i].Data()));
    if(!tracktab)continue;

    for(Int_t j=0; j<tracktab->GetEntries(); j++){
      UVdTrack* track = ( UVdTrack*)tracktab->At(j);
      Float_t ax = track->GetDX()/track->GetDZ();
      Float_t ay = track->GetDY()/track->GetDZ();

      if(i>1){
        fJuraProdCounter++;
        track->SetFlag(2);
        fhAxJ->Fill(ax);
        fhAyJ->Fill(ay);
        continue;
      }


      fhAxJ_flag02->Fill(ax);
      fhAyJ_flag02->Fill(ay);
      Float_t dax = ax-fOffAx_J;
      Float_t day = ay-fOffAy_J;
      Float_t dd = (dax*dax)/(fSigAx_J*fSigAx_J) +  (day*day)/(fSigAy_J*fSigAy_J);

      if(dd<fNsigJ*fNsigJ){
        track->SetFlag(0);
        fhAxJ_flag0->Fill(ax);
        fhAyJ_flag0->Fill(ay);
      }else{  
        fJuraProdCounter++;
        track->SetFlag(2);
        fhAxJ_flag2->Fill(ax);
        fhAyJ_flag2->Fill(ay);	
        fhAxJ->Fill(ax);
        fhAyJ->Fill(ay);
      }

    }
  }

  for(Int_t i=0;i<18;i++){
    UDataTable* tracktab = outS -> GetDataTable(Form("FullTracks %s",fMatchStr[i].Data()));
    if(!tracktab)continue;

    for(Int_t j=0; j<tracktab->GetEntries(); j++){
      UVdTrack* track = ( UVdTrack*)tracktab->At(j);

      Float_t ax = track->GetDX()/track->GetDZ();
      Float_t ay = track->GetDY()/track->GetDZ();

      if(i>1){
        fSaleveProdCounter++;
        track->SetFlag(2);
        fhAxS->Fill(ax);
        fhAyS->Fill(ay);
        continue;
      }


      fhAxS_flag02->Fill(ax);
      fhAyS_flag02->Fill(ay);
      Float_t dax = ax-fOffAx_S;
      Float_t day = ay-fOffAy_S;
      Float_t dd = (dax*dax)/(fSigAx_S*fSigAx_S) +  (day*day)/(fSigAy_S*fSigAy_S);

      if(dd<fNsigS*fNsigS){
        track->SetFlag(0);
        fhAxS_flag0->Fill(ax);
        fhAyS_flag0->Fill(ay);
      }else{  
        fSaleveProdCounter++;
        track->SetFlag(2);
        fhAxS->Fill(ax);
        fhAyS->Fill(ay);
        fhAxS_flag2->Fill(ax);
        fhAyS_flag2->Fill(ay);
      }


    }
  }


}



//_____________________________________________________________
void Na61PrimaryVertexRecoPostModule::FindPrimaryVertexPost(Int_t flag, UEventNode* outJ, UEventNode* outS)
{
  // algorithm:

  Int_t fullTracksJ = 0;
  Int_t fullTracksS = 0;
  Int_t fullTracksJ4h = 0;
  Int_t fullTracksS4h = 0;
  TObjArray tracks;
  tracks.Clear();

  Int_t histId = 0;    // Combined primary vertex
  if(!outS)histId = 1; // Jura primary vertex
  if(!outJ)histId = 2; // Saleve primary vertex

  //cout<<"---------------> histId="<<histId<<"  flag="<<flag<<endl;

  Int_t tabN = 18;
  //Int_t tabN = 4;

  Int_t jura_tracks = 0;
  if(outJ){
    for(Int_t it=0;it<tabN;it++){
      UDataTable* tracktab = outJ -> GetDataTable(Form("FullTracks %s",fMatchStr[it].Data()));
      if(!tracktab)continue;

      fullTracksJ = fullTracksJ + tracktab->GetEntries();
      if(it<4)fullTracksJ4h = fullTracksJ4h + tracktab->GetEntries();

      for(Int_t i=0; i<tracktab->GetEntries(); i++){
        UVdTrack* track = ( UVdTrack*)tracktab->At(i);
        if(track->GetFlag()==flag){ tracks.Add(tracktab->At(i)); jura_tracks++;}
      }
    }
  }

  Int_t saleve_tracks = 0;  
  if(outS){
    for(Int_t it=0;it<tabN;it++){
      UDataTable* tracktab = outS -> GetDataTable(Form("FullTracks %s",fMatchStr[it].Data()));
      if(!tracktab)continue;

      fullTracksS = fullTracksS + tracktab->GetEntries();
      if(it<4)fullTracksS4h = fullTracksS4h + tracktab->GetEntries();

      for(Int_t i=0; i<tracktab->GetEntries(); i++){
        UVdTrack* track = ( UVdTrack*)tracktab->At(i);
        if(track->GetFlag()==flag) { tracks.Add(tracktab->At(i)); saleve_tracks++;}
      }
    }
  }


  if(flag==2 && histId==0){
    fhFullTracks_JvsS -> Fill(fullTracksJ,fullTracksS);
    fhFullTracks_JS -> Fill(fullTracksJ+fullTracksS);
    fhFullTracks_J -> Fill(fullTracksJ);
    fhFullTracks_S -> Fill(fullTracksS);
    fhFullTracks4h_JS -> Fill(fullTracksJ4h + fullTracksS4h);
    fhFullTracks4h_J -> Fill(fullTracksJ4h);
    fhFullTracks4h_S -> Fill(fullTracksS4h);

    fAllTracksJ =  jura_tracks;
    fAllTracksS =  saleve_tracks;

  }

  //cout<<"--------------->1 histId="<<histId<<"  flag="<<flag<<"  "<<jura_tracks<<" "<<saleve_tracks
  //<<" all tracks: "<<tracks.GetEntries()<<endl;

  //if(tracks.GetEntries()<100) return;

  if((outJ && outS)){
    if((jura_tracks + saleve_tracks)<2) return;
  } 

  //cout<<"--------------->2 histId="<<histId<<"  flag="<<flag<<endl;

  Float_t sum_up_x = 0;
  Float_t sum_do_x = 0;
  Float_t sum_up_y = 0;
  Float_t sum_do_y = 0;

  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracks.At(i); 
    for(Int_t j=i+1;j<tracks.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracks.At(j); 

      Float_t axi = track1->GetDX_f()/track1->GetDZ_f();
      Float_t axj = track2->GetDX_f()/track2->GetDZ_f();
      Float_t bxi = track1->GetX_f();
      Float_t bxj = track2->GetX_f();

      Float_t ayi = track1->GetDY_f()/track1->GetDZ_f();
      Float_t ayj = track2->GetDY_f()/track2->GetDZ_f();
      Float_t byi = track1->GetY_f();
      Float_t byj = track2->GetY_f();

      Float_t z_x = ((axi-axj)*(bxi-bxj))/((axi-axj)*(axi-axj));
      Float_t z_y = ((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj));

      //cout<<z_x<<" "<<fZprim<<endl;
      if(TMath::Abs(z_x-fZprim) < 3 || flag!=2){
        //cout<<"acce: "<<z_x<<" "<<fZprim<<endl;
        sum_up_x = sum_up_x + (axi-axj)*(bxi-bxj);
        sum_do_x = sum_do_x + (axi-axj)*(axi-axj);
      }

      if(TMath::Abs(z_y-fZprim) < 3 || flag!=2){
        sum_up_y = sum_up_y + (ayi-ayj)*(byi-byj);
        sum_do_y = sum_do_y + (ayi-ayj)*(ayi-ayj);
      }

      //if(histId==0 && flag==2)cout<<"xx: "<< ((axi-axj)*(bxi-bxj))/ ((axi-axj)*(axi-axj))<<" yy: "<<((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj))<<endl;

    }
  } 

  ////////// new method //////////////////////////////////////////////////
  Float_t zx_min = 0;
  if(sum_do_x > 0)zx_min = - sum_up_x/sum_do_x;
  Float_t zy_min = 0;
  if(sum_do_y > 0)zy_min = - sum_up_y/sum_do_y; 

  //cout<<sum_up_x<<" "<<sum_do_x<<"       "<<sum_up_y<<" "<<sum_do_y<<"     "<<zx_min<<" "<<zy_min<<endl; 


  //Float_t z_prim = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);
  //Float_t z_prim = -(sum_up_x)/(sum_do_x);
  //Float_t z_prim = -(sum_up_y)/(sum_do_y); // use only  y information

  Float_t sigx=0.49; Float_t sigy=0.1065;
  if(histId==0) {sigx=0.12; sigy=0.07;}
  //if(histId==0) {sigx=0.1; sigy=0.1;}
  Float_t z_prim =  zx_min/(sigx*sigx)  +  zy_min/(sigy*sigy);
  Float_t norm = 1./(sigx*sigx)  +  1./(sigy*sigy);
  z_prim =  z_prim / norm;

  //cout<<"z_prim="<<z_prim<<endl;
  Vector3D reco_pvertex_w2;
  //PrimaryVertexWithWeigths(tracks,z_prim,reco_pvertex_w2);
  //PrimaryVertexWithWeigths_Y(tracks,z_prim,reco_pvertex_w2);
  //PrimaryVertexWithWeigths_X(tracks,z_prim,reco_pvertex_w2);
  PrimaryVertexWithWeigths_XY(flag,tracks,z_prim,reco_pvertex_w2,sigx,sigy);

  //cout<<"reco_pvertex_w2.Z= "<<reco_pvertex_w2.Z()<<" flag="<<flag<<" histId="<<histId<<endl;

  if(flag==2){
    if(fGoodVertex){
      fhRecoVertexXY[histId]     -> Fill(reco_pvertex_w2.X(),reco_pvertex_w2.Y());
      fhRecoVertexZ[histId]      -> Fill(reco_pvertex_w2.Z());
      fhRecoVertexZ_fine[histId] -> Fill(reco_pvertex_w2.Z());
    }

    if(histId==0){
      fhVtxX->Fill(reco_pvertex_w2.X());
      fhVtxY->Fill(reco_pvertex_w2.Y());
      fhTracksVsVtxX_J->Fill(reco_pvertex_w2.X(),fullTracksJ);
      fhTracksVsVtxX_S->Fill(reco_pvertex_w2.X(),fullTracksS);
      fhTracksVsVtxY_J->Fill(reco_pvertex_w2.Y(),fullTracksJ);
      fhTracksVsVtxY_S->Fill(reco_pvertex_w2.Y(),fullTracksS);

      fPrimaryVertex.SetX(reco_pvertex_w2.X());
      fPrimaryVertex.SetY(reco_pvertex_w2.Y());
      fPrimaryVertex.SetZ(reco_pvertex_w2.Z());
      //cout<<"reco_pvertex_w2.Z= "<<reco_pvertex_w2.Z()<<endl;
    }
  }
  if(flag==1){
    fhRecoVertexXY_flag1[histId]     -> Fill(reco_pvertex_w2.X(),reco_pvertex_w2.Y());
    fhRecoVertexZ_flag1[histId]      -> Fill(reco_pvertex_w2.Z());
    fhRecoVertexZ_fine_flag1[histId] -> Fill(reco_pvertex_w2.Z());
    //fPrimaryVertexDefined = kTRUE;
    //((UVdEvent*)out)->SetPrimaryVertexStatus(2);
  }
  if(outJ && (!outS))((UVdEvent*)outJ)->SetPrimaryVertex(reco_pvertex_w2.X(),reco_pvertex_w2.Y(),reco_pvertex_w2.Z());
  if(outS && (!outJ))((UVdEvent*)outS)->SetPrimaryVertex(reco_pvertex_w2.X(),reco_pvertex_w2.Y(),reco_pvertex_w2.Z());

}




//_____________________________________________________________
void Na61PrimaryVertexRecoPostModule::CheckPrimaryVertexResolution(Int_t flag, UEventNode* outJ, UEventNode* outS)
{
  // algorithm:

  TObjArray tracksSub1;
  tracksSub1.Clear();
  TObjArray tracksSub2;
  tracksSub2.Clear();

  Int_t tabN = 18;

  Int_t jura_tracks_Sub1 = 0;
  Int_t jura_tracks_Sub2 = 0;
  Int_t ii=0;
  if(outJ){
    for(Int_t it=0;it<tabN;it++){
      UDataTable* tracktab = outJ -> GetDataTable(Form("FullTracks %s",fMatchStr[it].Data()));
      if(!tracktab)continue;


      for(Int_t i=0; i<tracktab->GetEntries(); i++){
        UVdTrack* track = ( UVdTrack*)tracktab->At(i);

        if(!(track->GetFlag()==flag))continue;

        if((ii%2) == 0){tracksSub1.Add(track); jura_tracks_Sub1++;}
        else           {tracksSub2.Add(track); jura_tracks_Sub2++;}

        ii++;
      }
    }
  }

  Int_t saleve_tracks_Sub1 = 0;  
  Int_t saleve_tracks_Sub2 = 0;  
  ii=0;
  if(outS){
    for(Int_t it=0;it<tabN;it++){
      UDataTable* tracktab = outS -> GetDataTable(Form("FullTracks %s",fMatchStr[it].Data()));
      if(!tracktab)continue;

      for(Int_t i=0; i<tracktab->GetEntries(); i++){
        UVdTrack* track = ( UVdTrack*)tracktab->At(i);

        if(!(track->GetFlag()==flag)) continue;

        if((ii%2) == 0){tracksSub1.Add(track); saleve_tracks_Sub1++;}
        else           {tracksSub2.Add(track); saleve_tracks_Sub2++;}

        ii++;
      }
    }
  }



  if(tracksSub1.GetEntries()<2) return;
  if(tracksSub2.GetEntries()<2) return;

  if(tracksSub1.GetEntries()<10) return;

  //////////////////////////////////// Sub2 primary vertex //////////////////////////////////////////////////

  Float_t sum_up_x = 0;
  Float_t sum_do_x = 0;
  Float_t sum_up_y = 0;
  Float_t sum_do_y = 0;

  for(Int_t i=0;i<tracksSub1.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracksSub1.At(i); 
    for(Int_t j=i+1;j<tracksSub1.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracksSub1.At(j); 

      Float_t axi = track1->GetDX_f()/track1->GetDZ_f();
      Float_t axj = track2->GetDX_f()/track2->GetDZ_f();
      Float_t bxi = track1->GetX_f();
      Float_t bxj = track2->GetX_f();

      Float_t ayi = track1->GetDY_f()/track1->GetDZ_f();
      Float_t ayj = track2->GetDY_f()/track2->GetDZ_f();
      Float_t byi = track1->GetY_f();
      Float_t byj = track2->GetY_f();


      Float_t z_x = ((axi-axj)*(bxi-bxj))/((axi-axj)*(axi-axj));
      Float_t z_y = ((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj));

      if(TMath::Abs(z_x-fZprim) < 3 || flag!=2){
        sum_up_x = sum_up_x + (axi-axj)*(bxi-bxj);
        sum_do_x = sum_do_x + (axi-axj)*(axi-axj);
      }

      if(TMath::Abs(z_y-fZprim) < 3 || flag!=2){
        sum_up_y = sum_up_y + (ayi-ayj)*(byi-byj);
        sum_do_y = sum_do_y + (ayi-ayj)*(ayi-ayj);
      }    

    }
  } 


  Float_t zx_min = 0;
  if(sum_do_x > 0)zx_min = - sum_up_x/sum_do_x;
  Float_t zy_min = 0;
  if(sum_do_y > 0)zy_min = - sum_up_y/sum_do_y; 


  //cout<<sum_do_x<<"  "<<sum_do_y<<"  "<<zx_min<<" "<<zy_min<<endl; 


  //Float_t z_prim = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);
  //Float_t z_prim = -(sum_up_x)/(sum_do_x);
  //Float_t z_prim = -(sum_up_y)/(sum_do_y); // use only  y information

  //Float_t sigx=0.12; Float_t sigy=0.07;
  Float_t sigx=0.1; Float_t sigy=0.1;
  Float_t z_prim =  zx_min/(sigx*sigx)  +  zy_min/(sigy*sigy);
  Float_t norm = 1./(sigx*sigx)  +  1./(sigy*sigy);
  z_prim =  z_prim / norm;

  Vector3D reco_pvertex_w2_Sub1;
  //PrimaryVertexWithWeigths_Y(flag,tracksSub1,z_prim,reco_pvertex_w2_Sub1);
  //PrimaryVertexWithWeigths_X(flag,tracksSub1,z_prim,reco_pvertex_w2_Sub1);
  PrimaryVertexWithWeigths_XY(flag,tracksSub1,z_prim,reco_pvertex_w2_Sub1,sigx,sigy);
  //PrimaryVertexWithWeigths_XY2(flag,tracksSub1,z_prim,reco_pvertex_w2_Sub1,sigx,sigy);

  ////////////////////////////////////// Sub2 primary vertex //////////////////////////////////////////////////

  sum_up_x = 0;
  sum_do_x = 0;
  sum_up_y = 0;
  sum_do_y = 0;

  for(Int_t i=0;i<tracksSub2.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracksSub2.At(i); 
    for(Int_t j=i+1;j<tracksSub2.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracksSub2.At(j); 

      Float_t axi = track1->GetDX_f()/track1->GetDZ_f();
      Float_t axj = track2->GetDX_f()/track2->GetDZ_f();
      Float_t bxi = track1->GetX_f();
      Float_t bxj = track2->GetX_f();

      Float_t ayi = track1->GetDY_f()/track1->GetDZ_f();
      Float_t ayj = track2->GetDY_f()/track2->GetDZ_f();
      Float_t byi = track1->GetY_f();
      Float_t byj = track2->GetY_f();

      Float_t z_x = ((axi-axj)*(bxi-bxj))/((axi-axj)*(axi-axj));
      Float_t z_y = ((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj)); 

      if(TMath::Abs(z_x-fZprim) < 3 || flag!=2){      
        sum_up_x = sum_up_x + (axi-axj)*(bxi-bxj);
        sum_do_x = sum_do_x + (axi-axj)*(axi-axj);
      }

      if(TMath::Abs(z_y-fZprim) < 3 || flag!=2){   
        sum_up_y = sum_up_y + (ayi-ayj)*(byi-byj);
        sum_do_y = sum_do_y + (ayi-ayj)*(ayi-ayj);
      }
    }
  } 

  zx_min = 0;
  if(sum_do_x > 0)zx_min = - sum_up_x/sum_do_x;
  zy_min = 0;
  if(sum_do_y > 0)zy_min = - sum_up_y/sum_do_y; 

  //z_prim = -(sum_up_x + sum_up_y)/(sum_do_x+sum_do_y);
  //z_prim = -(sum_up_x)/(sum_do_x);
  //z_prim = -(sum_up_y)/(sum_do_y); // use only  y information

  z_prim =  zx_min/(sigx*sigx)  +  zy_min/(sigy*sigy);
  norm = 1./(sigx*sigx)  +  1./(sigy*sigy);
  z_prim =  z_prim / norm;

  Vector3D reco_pvertex_w2_Sub2;
  //PrimaryVertexWithWeigths_Y(tracksSub2,z_prim,reco_pvertex_w2_Sub2);
  //PrimaryVertexWithWeigths_X(tracksSub2,z_prim,reco_pvertex_w2_Sub2);
  PrimaryVertexWithWeigths_XY(flag,tracksSub2,z_prim,reco_pvertex_w2_Sub2,sigx,sigy);
  //PrimaryVertexWithWeigths_XY2(tracksSub2,z_prim,reco_pvertex_w2_Sub2,sigx,sigy);

  Int_t pvertStatusJ = ((UVdEvent*)outJ)->GetPrimaryVertexStatus();
  Int_t pvertStatusS = ((UVdEvent*)outS)->GetPrimaryVertexStatus();

  if(pvertStatusJ == 2 && pvertStatusS == 2){
    fhVtxDx_Glob2 -> Fill(reco_pvertex_w2_Sub2.X() - reco_pvertex_w2_Sub1.X());
    fhVtxDy_Glob2 -> Fill(reco_pvertex_w2_Sub2.Y() - reco_pvertex_w2_Sub1.Y());
    fhVtxDz_Glob2 -> Fill(reco_pvertex_w2_Sub2.Z() - reco_pvertex_w2_Sub1.Z());
    if(fGoodVertex){
      fhVtxDx_Glob2Acce -> Fill(reco_pvertex_w2_Sub2.X() - reco_pvertex_w2_Sub1.X());
      fhVtxDy_Glob2Acce -> Fill(reco_pvertex_w2_Sub2.Y() - reco_pvertex_w2_Sub1.Y());
      fhVtxDz_Glob2Acce -> Fill(reco_pvertex_w2_Sub2.Z() - reco_pvertex_w2_Sub1.Z());
    }
  }

}


//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoPostModule::PrimaryVertexWithWeigths_XY(Int_t flag, TObjArray& tracks,Float_t zprim,
    Vector3D& pvertex,Float_t sigx, Float_t sigy)
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

  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracks.At(i); 
    for(Int_t j=i+1;j<tracks.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracks.At(j); 

      Float_t x1 = track1->GetXatZ_f(zprim);
      Float_t x2 = track2->GetXatZ_f(zprim);
      Float_t y1 = track1->GetYatZ_f(zprim);
      Float_t y2 = track2->GetYatZ_f(zprim);

      Float_t d2x = (x2-x1)*(x2-x1);
      Float_t d2y = (y2-y1)*(y2-y1);      

      Float_t wx = 1;
      Float_t wy = 1;
      if(d2x>0)wx = 1./d2x;
      if(d2y>0)wy = 1./d2y;

      //fhd2->Fill(d2);

      /////////////////// stuff for new method
      Float_t axi = track1->GetDX_f()/track1->GetDZ_f();
      Float_t axj = track2->GetDX_f()/track2->GetDZ_f();
      Float_t bxi = track1->GetX_f();
      Float_t bxj = track2->GetX_f();

      Float_t ayi = track1->GetDY_f()/track1->GetDZ_f();
      Float_t ayj = track2->GetDY_f()/track2->GetDZ_f();
      Float_t byi = track1->GetY_f();
      Float_t byj = track2->GetY_f();


      Float_t z_x = ((axi-axj)*(bxi-bxj))/((axi-axj)*(axi-axj));
      Float_t z_y = ((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj));

      if(TMath::Abs(z_x-fZprim) < 3 || flag!=2){
        sum_up_x = sum_up_x + wx*(axi-axj)*(bxi-bxj);
        sum_do_x = sum_do_x + wx*(axi-axj)*(axi-axj);
      }
      if(TMath::Abs(z_y-fZprim) < 3 || flag!=2){
        sum_up_y = sum_up_y + wy*(ayi-ayj)*(byi-byj);
        sum_do_y = sum_do_y + wy*(ayi-ayj)*(ayi-ayj);
      }
      //cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;

    }
  } 

  Float_t zprim_w2x = 0;
  if(sum_do_x > 0)zprim_w2x = - sum_up_x/sum_do_x;
  Float_t zprim_w2y = 0;
  if(sum_do_y > 0)zprim_w2y = - sum_up_y/sum_do_y; 

  //cout<<" _XY:"<< sum_up_x<<" "<< sum_do_x<<"    "<<sum_up_y<<" "<< sum_do_y<<"   "<<zprim_w2x<<" "<<zprim_w2y<<endl;

  Float_t zprim_w2 =  zprim_w2x/(sigx*sigx)  +   zprim_w2y/(sigy*sigy);
  Float_t norm = 1./(sigx*sigx)  +  1./(sigy*sigy);
  zprim_w2 =  zprim_w2 / norm;

  //cout<<" _XY: zprim_w2 = "<<zprim_w2<<" norm="<<norm<<endl;

  //
  Float_t sum_x = 0;
  Float_t sum_y = 0;
  Float_t N = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 

    Float_t x = track->GetXatZ_f(zprim_w2);
    Float_t y = track->GetYatZ_f(zprim_w2);

    sum_x = sum_x + x;      
    sum_y = sum_y + y;
    N = N + 1;      
  }

  Float_t xprim = sum_x/N;
  Float_t yprim = sum_y/N;

  sum_x = 0;
  sum_y = 0;

  Float_t norm_x = 0;
  Float_t norm_y = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 

    Float_t x = track->GetXatZ_f(zprim_w2);
    Float_t y = track->GetYatZ_f(zprim_w2);

    Float_t d2x = (x-xprim)*(x-xprim);
    Float_t d2y = (y-yprim)*(y-yprim);

    sum_x = sum_x + x/d2x;      
    sum_y = sum_y + y/d2y;
    norm_x = norm_x + 1./d2x;      
    norm_y = norm_y + 1./d2y;      
  }

  Float_t xprim_w2 = sum_x/norm_x;
  Float_t yprim_w2 = sum_y/norm_y;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);

}

//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoPostModule::PrimaryVertexWithWeigths_XY2(Int_t flag, TObjArray& tracks,Float_t zprim,Vector3D& pvertex,
    Float_t /* sigx */, Float_t /* sigy */)
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

  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracks.At(i); 
    for(Int_t j=i+1;j<tracks.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracks.At(j); 

      Float_t x1 = track1->GetXatZ_f(zprim);
      Float_t x2 = track2->GetXatZ_f(zprim);
      Float_t y1 = track1->GetYatZ_f(zprim);
      Float_t y2 = track2->GetYatZ_f(zprim);

      Float_t d2x = (x2-x1)*(x2-x1);
      Float_t d2y = (y2-y1)*(y2-y1);      
      Float_t wx = 1./d2x;
      Float_t wy = 1./d2y;

      //fhd2->Fill(d2);

      /////////////////// stuff for new method
      Float_t axi = track1->GetDX_f()/track1->GetDZ_f();
      Float_t axj = track2->GetDX_f()/track2->GetDZ_f();
      Float_t bxi = track1->GetX_f();
      Float_t bxj = track2->GetX_f();

      Float_t ayi = track1->GetDY_f()/track1->GetDZ_f();
      Float_t ayj = track2->GetDY_f()/track2->GetDZ_f();
      Float_t byi = track1->GetY_f();
      Float_t byj = track2->GetY_f();

      Float_t z_x = ((axi-axj)*(bxi-bxj))/((axi-axj)*(axi-axj));
      Float_t z_y = ((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj));

      if(TMath::Abs(z_x-50) < 50 || flag!=2){
        sum_up_x = sum_up_x + wx*(axi-axj)*(bxi-bxj);
        sum_do_x = sum_do_x + wx*(axi-axj)*(axi-axj);
      }

      if(TMath::Abs(z_y-50) < 50 || flag!=2){
        sum_up_y = sum_up_y + wy*(ayi-ayj)*(byi-byj);
        sum_do_y = sum_do_y + wy*(ayi-ayj)*(ayi-ayj);
        //cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;
      }
    }
  } 

  //Float_t zprim_w2x = -(sum_up_x)/(sum_do_x);
  //Float_t zprim_w2y = -(sum_up_y)/(sum_do_y);

  Float_t zprim_w2 = -(sum_up_x + sum_up_y)/(sum_do_x + sum_do_y);
  //Float_t norm = 1./(sigx*sigx)  +  1./(sigy*sigy);
  //zprim_w2 =  zprim_w2 / norm;

  //cout<<" zprim_w2 = "<<zprim_w2<<endl;

  //
  Float_t sum_x = 0;
  Float_t sum_y = 0;
  Float_t N = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 

    Float_t x = track->GetXatZ_f(zprim_w2);
    Float_t y = track->GetYatZ_f(zprim_w2);

    sum_x = sum_x + x;      
    sum_y = sum_y + y;
    N = N + 1;      
  }

  Float_t xprim = sum_x/N;
  Float_t yprim = sum_y/N;

  sum_x = 0;
  sum_y = 0;

  Float_t norm_x = 0;
  Float_t norm_y = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 

    Float_t x = track->GetXatZ_f(zprim_w2);
    Float_t y = track->GetYatZ_f(zprim_w2);

    Float_t d2x = (x-xprim)*(x-xprim);
    Float_t d2y = (y-yprim)*(y-yprim);

    sum_x = sum_x + x/d2x;      
    sum_y = sum_y + y/d2y;
    norm_x = norm_x + 1./d2x;      
    norm_y = norm_y + 1./d2y;      
  }

  Float_t xprim_w2 = sum_x/norm_x;
  Float_t yprim_w2 = sum_y/norm_y;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);

}

//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoPostModule::PrimaryVertexWithWeigths_X(Int_t flag, TObjArray& tracks,Float_t zprim,Vector3D& pvertex)
{
  // algorithm: 
  // 1. Calculate weights = 1/dist^2, distance betweem i-th and j-th tracks on
  // the closest proximity plane (z=zprim).
  // 2. Used the weights to calculate z coordinate of primary vertex from the formula that inclused weights (see PS notes).  
  // 3. Find x,y coordinates averageing tracks x, y at the closes proximmity plane. 

  Float_t sum_up_x = 0;
  Float_t sum_do_x = 0;
  //Float_t sum_up_y = 0;
  //Float_t sum_do_y = 0;

  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracks.At(i); 
    for(Int_t j=i+1;j<tracks.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracks.At(j); 

      Float_t x1 = track1->GetXatZ_f(zprim);
      Float_t x2 = track2->GetXatZ_f(zprim);

      Float_t d2 = (x2-x1)*(x2-x1);
      Float_t w = 1./d2;

      //fhd2->Fill(d2);

      /////////////////// stuff for new method
      Float_t axi = track1->GetDX_f()/track1->GetDZ_f();
      Float_t axj = track2->GetDX_f()/track2->GetDZ_f();
      Float_t bxi = track1->GetX_f();
      Float_t bxj = track2->GetX_f();

      //Float_t ayi = track1->GetDY()/track1->GetDZ();
      //Float_t ayj = track2->GetDY()/track2->GetDZ();
      //Float_t byi = track1->GetY();
      //Float_t byj = track2->GetY();

      Float_t z_x = ((axi-axj)*(bxi-bxj))/((axi-axj)*(axi-axj));

      if(TMath::Abs(z_x-50) < 50 || flag!=2){	
        sum_up_x = sum_up_x + w*(axi-axj)*(bxi-bxj);
        sum_do_x = sum_do_x + w*(axi-axj)*(axi-axj);
      }
      //sum_up_y = sum_up_y + w*(ayi-ayj)*(byi-byj);
      //sum_do_y = sum_do_y + w*(ayi-ayj)*(ayi-ayj);
      //cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;

    }
  } 

  Float_t zprim_w2 = -(sum_up_x)/(sum_do_x);
  //Float_t zprim_w2 = -(sum_up_y)/(sum_do_y);

  //
  Float_t sum_x = 0;
  Float_t sum_y = 0;
  Float_t N = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 

    Float_t x = track->GetXatZ_f(zprim_w2);
    Float_t y = track->GetYatZ_f(zprim_w2);

    sum_x = sum_x + x;      
    sum_y = sum_y + y;
    N = N + 1;      
  }

  Float_t xprim = sum_x/N;
  Float_t yprim = sum_y/N;

  sum_x = 0;
  sum_y = 0;

  Float_t norm_x = 0;
  Float_t norm_y = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 

    Float_t x = track->GetXatZ_f(zprim_w2);
    Float_t y = track->GetYatZ_f(zprim_w2);
    //cout<<"x="<<x<<endl;

    Float_t d2x = (x-xprim)*(x-xprim);
    Float_t d2y = (y-yprim)*(y-yprim);

    sum_x = sum_x + x/d2x;      
    sum_y = sum_y + y/d2y;
    norm_x = norm_x + 1./d2x;      
    norm_y = norm_y + 1./d2y;      
  }

  Float_t xprim_w2 = sum_x/norm_x;
  Float_t yprim_w2 = sum_y/norm_y;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);

}

//______________________________________________________________________________________________________
void Na61PrimaryVertexRecoPostModule::PrimaryVertexWithWeigths_Y(Int_t flag, TObjArray& tracks,Float_t zprim,Vector3D& pvertex)
{
  // algorithm: 
  // 1. Calculate weights = 1/dist^2, distance betweem i-th and j-th tracks on
  // the closest proximity plane (z=zprim).
  // 2. Used the weights to calculate z coordinate of primary vertex from the formula that inclused weights (see PS notes).  
  // 3. Find x,y coordinates averageing tracks x, y at the closes proximmity plane. 

  //Float_t sum_up_x = 0;
  //Float_t sum_do_x = 0;
  Float_t sum_up_y = 0;
  Float_t sum_do_y = 0;

  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track1 = (UVdTrack*)tracks.At(i); 
    for(Int_t j=i+1;j<tracks.GetEntries();j++){
      UVdTrack* track2 = (UVdTrack*)tracks.At(j); 

      Float_t y1 = track1->GetYatZ_f(zprim);
      Float_t y2 = track2->GetYatZ_f(zprim);

      Float_t d2 = (y2-y1)*(y2-y1);      
      Float_t w = 1./d2;

      //fhd2->Fill(d2);

      /////////////////// stuff for new method
      //Float_t axi = track1->GetDX()/track1->GetDZ();
      //Float_t axj = track2->GetDX()/track2->GetDZ();
      //Float_t bxi = track1->GetX();
      //Float_t bxj = track2->GetX();

      Float_t ayi = track1->GetDY_f()/track1->GetDZ_f();
      Float_t ayj = track2->GetDY_f()/track2->GetDZ_f();
      Float_t byi = track1->GetY_f();
      Float_t byj = track2->GetY_f();

      //sum_up_x = sum_up_x + wx*(axi-axj)*(bxi-bxj);
      //sum_do_x = sum_do_x + wx*(axi-axj)*(axi-axj);

      Float_t z_y = ((ayi-ayj)*(byi-byj))/((ayi-ayj)*(ayi-ayj));

      if(TMath::Abs(z_y-50) < 50 || flag!=2){	
        sum_up_y = sum_up_y + w*(ayi-ayj)*(byi-byj);
        sum_do_y = sum_do_y + w*(ayi-ayj)*(ayi-ayj);
      }
      //cout<<"i="<<i<<" j="<<j<<"  "<<(axi-axj)*(axi-axj)<<" "<<axi<<" "<<axj<<endl;

    }
  } 

  //Float_t zprim_w2_x = -(sum_up_x)/(sum_do_x);
  Float_t zprim_w2 = -(sum_up_y)/(sum_do_y);

  //
  Float_t sum_x = 0;
  Float_t sum_y = 0;
  Float_t N = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 

    Float_t x = track->GetXatZ_f(zprim_w2);
    Float_t y = track->GetYatZ_f(zprim_w2);

    sum_x = sum_x + x;      
    sum_y = sum_y + y;
    N = N + 1;      
  }

  Float_t xprim = sum_x/N;
  Float_t yprim = sum_y/N;

  sum_x = 0;
  sum_y = 0;

  Float_t norm_x = 0;
  Float_t norm_y = 0;
  for(Int_t i=0;i<tracks.GetEntries();i++){
    UVdTrack* track = (UVdTrack*)tracks.At(i); 

    Float_t x = track->GetXatZ_f(zprim_w2);
    Float_t y = track->GetYatZ_f(zprim_w2);

    Float_t d2x = (x-xprim)*(x-xprim);
    Float_t d2y = (y-yprim)*(y-yprim);

    sum_x = sum_x + x/d2x;      
    sum_y = sum_y + y/d2y;
    norm_x = norm_x + 1./d2x;      
    norm_y = norm_y + 1./d2y;      
  }

  Float_t xprim_w2 = sum_x/norm_x;
  Float_t yprim_w2 = sum_y/norm_y;

  pvertex.SetX(xprim_w2);
  pvertex.SetY(yprim_w2);
  pvertex.SetZ(zprim_w2);

}


//_____________________________________________________________
void Na61PrimaryVertexRecoPostModule::End()
{
  // Run-level finalisation
  SetState(kEnd);
  Na61PrimaryVertexRecoModule::End();

}

//____________________________________________________________________
void Na61PrimaryVertexRecoPostModule::Finish()
{
  // Job-level finalisation
  SetState(kFinish);
  Na61PrimaryVertexRecoModule::Finish();

}



//____________________________________________________________________
void Na61PrimaryVertexRecoPostModule::Print(Option_t* option) const
{
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>
  Na61PrimaryVertexRecoModule::Print(option);

}



