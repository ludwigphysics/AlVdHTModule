//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61VdTpcMatchingModule
#include "Na61VdTpcMatchingModule.h"    
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
#include <map>

#include <modutils/KalmanFitterAK.h>
#include <modutils/KalmanFilterWB.h>

//____________________________________________________________________
//ClassImp(Na61VdTpcMatchingModule);

//____________________________________________________________________
Na61VdTpcMatchingModule::Na61VdTpcMatchingModule()
{
  // Default constructor. DO NOT USE
  SetState(kSetup);
}

//____________________________________________________________________
Na61VdTpcMatchingModule::Na61VdTpcMatchingModule(const Char_t* name, const Char_t* title)
   : Na61Module(name, title)
{
  // Named Constructor
  SetState(kSetup);  

  // for now just use fixed run number but if needed the run# can be given from main()
  /*  const unsigned int runNumber = 27441;
  const TimeStamp timeStamp = 
    UTCDateTime(2016, 12, 12, 23, 14, 59).GetTimeStamp();
  det::Detector& detector = det::Detector::GetInstance();
  detector.Update(timeStamp, runNumber);
  
  ftracker = detector.GetMagneticFieldTracker();
  */
  fTpcTracks = 0;
  fNumberOfVdTpcTracks = 0;
  fNSig = 4.0;
  fMomDiffVtx  = 0.007;
  fMomDiffVtx2 = 0.06;
}

//____________________________________________________________________
void Na61VdTpcMatchingModule::DefineHistograms()
{
  
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init"); 
    return;  
    
  }
  
  TDirectory* histDir = gDirectory->mkdir("VdTpcMatchingModule");
  
  histDir->cd(); 
  

  hMomAllPos = new TH1F("hMomAllPos","",1500,-0.5,149.5);
  hMomAllNeg = new TH1F("hMomAllNeg","",1500,-0.5,149.5);
  hMomAllTpc = new TH1F("hMomAllTpc","",1500,-0.5,149.5);
  hPtAllTpc = new TH1F("hPtAllTpc","",500, 0, 5);


  fhMomPoints_ZX      =	new TH2F("hMomPoints_ZX","",1000,-700,-100,1000,-100,100);
  fhMomPoints_ZX_fine =	new TH2F("hMomPoints_ZX_fine","",1000,-620,-560,1000,-10,10);
 
 
  fhTanY = new TH1F("hTanY", "", 2000, -.2, .2);
  fhTanY_acce = new TH1F("hTanY_acce", "", 2000, -.2, .2);
  fhTanYvsY = new TH2F("hTanYvsY", "", 50,-25,25, 2000, -.2, .2);
  fhTanYvsY_acce = new TH2F("hTanYvsY_acce", "", 50,-25,25, 2000, -.2, .2);
  
  fhY_JJ_vsY_cutx = new TH2F("hY_JJ_vsY_cutx", "", 60,-35,25, 500, -20., 20.);
  fhY_SS_vsY_cutx = new TH2F("hY_SS_vsY_cutx", "", 60,-35,25, 500, -20., 20.);
  fhY_JS_vsY_cutx = new TH2F("hY_JS_vsY_cutx", "", 60,-35,25, 500, -20., 20.);
  fhY_SJ_vsY_cutx = new TH2F("hY_SJ_vsY_cutx", "", 60,-35,25, 500, -20., 20.);
  
  fhY_JJ_vsY_cutx_r1 = new TH2F("hY_JJ_vsY_cutx_r1", "", 60,-35,25, 500, -20., 20.);
  fhY_SS_vsY_cutx_r1 = new TH2F("hY_SS_vsY_cutx_r1", "", 60,-35,25, 500, -20., 20.);
  fhY_JS_vsY_cutx_r1 = new TH2F("hY_JS_vsY_cutx_r1", "", 60,-35,25, 500, -20., 20.);
  fhY_SJ_vsY_cutx_r1 = new TH2F("hY_SJ_vsY_cutx_r1", "", 60,-35,25, 500, -20., 20.);
  
  fhY_JJ_vsY_cutx_r2 = new TH2F("hY_JJ_vsY_cutx_r2", "", 60,-35,25, 500, -20., 20.);
  fhY_SS_vsY_cutx_r2 = new TH2F("hY_SS_vsY_cutx_r2", "", 60,-35,25, 500, -20., 20.);
  fhY_JS_vsY_cutx_r2 = new TH2F("hY_JS_vsY_cutx_r2", "", 60,-35,25, 500, -20., 20.);
  fhY_SJ_vsY_cutx_r2 = new TH2F("hY_SJ_vsY_cutx_r2", "", 60,-35,25, 500, -20., 20.);
  
  fhX_JJ_vsZ_cuty = new TH2F("hX_JJ_vsZ_cuty", "", 200,-500,120, 500, -20., 20.);
  fhX_SS_vsZ_cuty = new TH2F("hX_SS_vsZ_cuty", "", 200,-500,120, 500, -20., 20.);
  fhX_JS_vsZ_cuty = new TH2F("hX_JS_vsZ_cuty", "", 200,-500,120, 500, -20., 20.);
  fhX_SJ_vsZ_cuty = new TH2F("hX_SJ_vsZ_cuty", "", 200,-500,120, 500, -20., 20.);
  
  fhXvsPz_J = new TH2F("hXvsPz_J", "", 500,-100,100, 500, -1., 1.);
  fhXvsPz_S = new TH2F("hXvsPz_S", "", 500,-100,100, 500, -1., 1.);
  fhXvsPx_J = new TH2F("hXvsPx_J", "", 500,-100,100, 500, -1., 1.);
  fhXvsPx_S = new TH2F("hXvsPx_S", "", 500,-100,100, 500, -1., 1.);
  
  fhPosMatch[0] = new TH1F("hPosMatch_X", "",  1000, -200, 200);
  fhPosMatch[1] = new TH1F("hPosMatch_Y", "",  1000, -200, 200);
  fhPosMatch[2] = new TH1F("hPosMatch_Z", "",  1000, -200, 200);

  fhPosMatch_best[0] = new TH1F("hPosMatch_best_X", "",  1000, -200, 200);
  fhPosMatch_best[1] = new TH1F("hPosMatch_best_Y", "",  1000, -200, 200);
  
  fhMomMatch[0] = new TH1F("hMomMatch_X", "",  1000, -1, 1);
  fhMomMatch[1] = new TH1F("hMomMatch_Y", "",  1000, -0.6, 0.6);
  fhMomMatch[2] = new TH1F("hMomMatch_Z", "",  1000, -0.4, 0.4);
  
  fhMomMatch_cutx[0] = new TH1F("hMomMatch_cutx_X", "",  1000, -1, 1);
  fhMomMatch_cutx[1] = new TH1F("hMomMatch_cutx_Y", "",  1000, -0.6, 0.6);
  fhMomMatch_cutx[2] = new TH1F("hMomMatch_cutx_Z", "",  1000, -0.4, 0.4);
  
  fhMomMatch_cuty[0] = new TH1F("hMomMatch_cuty_X", "",  1000, -1, 1);
  fhMomMatch_cuty[1] = new TH1F("hMomMatch_cuty_Y", "",  1000, -0.6, 0.6);
  fhMomMatch_cuty[2] = new TH1F("hMomMatch_cuty_Z", "",  1000, -0.4, 0.4);
  
  fhMomMatch_cutxy[0] = new TH1F("hMomMatch_cutxy_X", "",  1000, -1, 1);
  fhMomMatch_cutxy[1] = new TH1F("hMomMatch_cutxy_Y", "",  1000, -0.6, 0.6);
  fhMomMatch_cutxy[2] = new TH1F("hMomMatch_cutxy_Z", "",  1000, -0.4, 0.4);
  
  fhMomMatch_best[0] = new TH1F("hMomMatch_best_X", "",  1000, -1, 1);
  fhMomMatch_best[1] = new TH1F("hMomMatch_best_Y", "",  1000, -0.6, 0.6);
  fhMomMatch_best[2] = new TH1F("hMomMatch_best_Z", "",  1000, -0.4, 0.4);

  fhMatch_TanYvsY_acce = new TH2F("hMatch_TanYvsY_acce", "", 2000, -.2, .2, 500, -100, 100);
  fhMatch_TanYvsY_best = new TH2F("hMatch_TanYvsY_best", "", 2000, -.2, .2, 500, -100, 100);

  fhChargeCorr = new TH2F("hChargCorr", "", 3, -1.5, 1.5, 3, -1.5 ,1.5);
  fhMomVsKfMom = new TH2F("hMomVsKfMom", "", 500, 0, 50, 500, 0 ,50);
  fhMomVsKfMom_kf = new TH2F("hMomVsKfMom_kf", "", 500, 0, 50, 500, 0 ,50);
  fhMomVsKfMom_acce = new TH2F("hMomVsKfMom_acce", "", 500, 0, 50, 500, 0 ,50);
  fhMomDiff = new TH1F("hMomDiff","",1000,-1,1);
  fhMomDiff_kf = new TH1F("hMomDiff_kf","",1000,-0.01,0.01);
  

  fhTdistVsMomDiff = new TH2F("hTdistVdMomDiff", "", 5000, -1, 1, 1000, 0,50);
  fhXY_AtTarget = new TH2F("hXY_AtTarget", "", 500, -10, 10, 500, -10,10);
  fhXY_AtTarget_clean = new TH2F("hXY_AtTarget_clean", "", 500, -10, 10, 500, -10,10);

    
  fhMainVertexVsVdVertex_X = new TH1F("hMainVertexVsVdVertex_X","",500,-12.,12.);    
  fhMainVertexVsVdVertex_Y = new TH1F("hMainVertexVsVdVertex_Y","",500,-6.,6.);    
  fhMainVertexVsVdVertex_Z = new TH1F("hMainVertexVsVdVertex_Z","",500,-60.,60.);    
 
    
  fhTracksVsTracksWithMatch = new TH2F("hTracksVsTracksWithMatch","",500,0,500,500,0,500);

  fhXY_Vds4_All = new TH2F("hXY_Vds4_All","",500,-5,5,500,-4,4);
  fhXY_Vds4_WithMatch = new TH2F("hXY_Vds4_WithMatch","",500,-5,5,500,-4,4);
  fhXY_VTpc1_All = new TH2F("hXY_VTpc1_All","",500,-50,50,500,-25,25);
  fhXY_VTpc1_WithMatch = new TH2F("hXY_VTpc1_WithMatch","",500,-50,50,500,-25,25);

  fhXY_Vds4_NoMatch = new TH2F("hXY_Vds4_NoMatch","",500,-24,24,500,-12,12);
  fhXY_fTPC_NoMatch = new TH2F("hXY_fTPC_NoMatch","",500,-70,70,500,-35,35);
  fhXY_Vds4_Match = new TH2F("hXY_Vds4_Match","",500,-24,24,500,-12,12);
  fhXY_fTPC_Match = new TH2F("hXY_fTPC_Match","",500,-70,70,500,-35,35);

  fhCurvature_All = new TH1F("hCurvature_All","",500,-0.00003,0.00003);
  fhCurvature_WithMatch = new TH1F("hCurvature_WithMatch","",500,-0.00003,0.00003);

  fhCurvature_m_vs_m = new TH2F("hCurvature_m_vs_m","",500,-0.00003,0.00003,500,-0.00003,0.00003);
  fhCurvature_m_vs_nm = new TH2F("hCurvature_m_vs_nm","",500,-0.00003,0.00003,500,-0.00003,0.00003);

  fhXOnVds4_m_vs_m = new TH2F("hXOnVds4_m_vs_m","",500,-50.,50,500,-50,50);
  fhXOnVds4_m_vs_nm = new TH2F("hXOnVds4_m_vs_nm","",500,-50.,50,500,-50,50);
  fhYOnVds4_m_vs_m = new TH2F("hYOnVds4_m_vs_m","",500,-50.,50,500,-50,50);
  fhYOnVds4_m_vs_nm = new TH2F("hYOnVds4_m_vs_nm","",500,-50.,50,500,-50,50);

  fhClusters_All = new TH1F("hClusters_All","",500,0,500);
  fhClusters_All_Acce = new TH1F("fhClusters_All_Acce","",500,0,500);
  fhClusters_WithMatch = new TH1F("hClusters_WithMatch","",500,0,500);

  fhdEdx_all = new TH2F("hdEdx_all","",500,0,50,500,0,5);
  fhdEdxPos_all = new TH2F("hdEdxPos_all","",500,0,50,500,0,5);
  fhdEdxNeg_all = new TH2F("hdEdxNeg_all","",500,0,50,500,0,5);
  fhdEdx_mtpc = new TH2F("hdEdx_mtpc","",500,0,50,500,0,5);
  fhdEdxPos_mtpc = new TH2F("hdEdxPos_mtpc","",500,0,50,500,0,5);
  fhdEdxNeg_mtpc = new TH2F("hdEdxNeg_mtpc","",500,0,50,500,0,5);
  fhdEdxVsClusters_all = new TH2F("hdEdxVsClusters_all","",500,0,500,500,0,5);
  fhdEdxVsClustersCut_all = new TH2F("hdEdxVsClustersCut_all","",500,0,500,500,0,5);
  fhdEdxVsClusters_mtpc = new TH2F("hdEdxVsClusters_mtpc","",500,0,500,500,0,5);
  fhdEdxVsClustersCut_mtpc = new TH2F("hdEdxVsClustersCut_mtpc","",500,0,500,500,0,5);
  fhClustersVsMom_all = new TH2F("hClustersVsMom_all","",500,0,50,500,0,500);
  
    fhPt = new TH1F("hPt", "",  500, 0, 5.);
    

  gDirectory->cd("..");
  
}

//____________________________________________________________________
void Na61VdTpcMatchingModule::Init()
{
  // Job-level initialisation
  SetState(kInit);
  
  /*
  if(fJura)fParams = (Na61VdParametersManager::Instance())->GetJuraArmParams();  
  else fParams = (Na61VdParametersManager::Instance())->GetSaleveArmParams(); 

  cout<<"Na61VdTpcMatchingModule::Init: fParams="<<(Na61ArmParameters*)fParams<<"  "<<fParams->GetInit()<<endl;

  if(!fParams->GetInit()){
    fParams -> SetDrotX(0,0); 
    fParams -> SetDrotY(0,0); 
    fParams -> SetDrotZ(0,0); 
    fParams -> SetVolumeDz(0,0);
    
    fParams -> SetRunId(fRunId);
    fParams->Init();
    }*/


  fVdParams = (Na61VdParametersManager::Instance())->GetVdParams();  
  fVdParams -> SetRunId(fRunId);
  if(!fVdParams->GetInit())fVdParams->Init();

  fOffsetX = fVdParams->GetOffsetToTpcX();
  fOffsetY = fVdParams->GetOffsetToTpcY();
  fOffsetZ = fVdParams->GetOffsetToTpcZ();
  cout<<"Na61VdTpcMatchingModule::Init:  OffsetToTpcX: "<< fOffsetX <<"  OffsetToTpcY: "<<fOffsetY<<"  OffsetToTpcZ: "<<fOffsetZ<<endl;  
  //cout<<"(VdTpc matching cut: fNSig = "<<fNSig<<endl;
}

//____________________________________________________________________
void Na61VdTpcMatchingModule::Begin()
{
  // Run-level initialisation
  SetState(kBegin);
  
}

//____________________________________________________________________
void Na61VdTpcMatchingModule :: Event(const det::MagneticFieldTracker& tra,evt::Event& event, UEventNode* inNodeJ, UEventNode* inNodeS, UEventNode* outNode)
{
  // Per event method
  SetState(kEvent);  
  
  
  RecEvent& recEvent = event.GetRecEvent();
  
  //const rec::Beam& beam = recEvent.GetBeam();
  
  //const rec::PSD& psd = recEvent.GetPSD();
  
  //double psde = psd.GetEnergy()/208.;
  
  if (!recEvent.HasMainVertex()) cout<<"No Main Vertex!"<<endl;
  if (!recEvent.HasPrimaryVertex(evt::rec::VertexConst::ePrimaryFitZ)) cout<<"No ePrimaryFitZ Vertex!"<<endl;
  //if (!recEvent.HasPrimaryVertex(evt::rec::VertexConst::ePrimary3D)) cout<<"No ePrimary3D Vertex!"<<endl;
  
  if (!recEvent.HasMainVertex()) return;
  const Vertex& mainVertex = recEvent.GetMainVertex();
  
  Int_t pvertStatus  = ((UVdEvent*)outNode)->GetPrimaryVertexStatus();
  
  UDataTable* tpctracks = new UDataTable("Tpc Tracks");
  tpctracks -> SetOwner();
  
  if(pvertStatus){	    
    
    //cout<<"main vertex:  pvS="<<pvertStatus<<" pvSJ="<<pvertStatusJ<<"  pvSS="<<pvertStatusS<<"  "<<vertex_x_VD<<" "<<vertex_y_VD<<" "<<vertex_z_VD<<endl;
    
    double vertex_x_VD = ((UVdEvent*)outNode)->GetPrimaryVertexX() + fOffsetX;
    double vertex_y_VD = ((UVdEvent*)outNode)->GetPrimaryVertexY() + fOffsetY;
    double vertex_z_VD = ((UVdEvent*)outNode)->GetPrimaryVertexZ() + fOffsetZ;
    
    //cout<<" vertex_z_VD ="<< vertex_z_VD/10.<<"  vertex_z_TPC="<<vertex_z_TPC<<endl; 
    
    Point mainVertexVD(vertex_x_VD, vertex_y_VD, vertex_z_VD);
    
    fhMomPoints_ZX->Fill(vertex_z_VD/10., vertex_x_VD/10.);
    fhMomPoints_ZX_fine->Fill(vertex_z_VD/10., vertex_x_VD/10.);
    
    fhMainVertexVsVdVertex_X -> Fill(mainVertex.GetPosition().GetX()*10 - vertex_x_VD);
    fhMainVertexVsVdVertex_Y -> Fill(mainVertex.GetPosition().GetY()*10 - vertex_y_VD);
    fhMainVertexVsVdVertex_Z -> Fill(mainVertex.GetPosition().GetZ()*10 - vertex_z_VD);
    
    TObjArray vdtracks;
    vdtracks.Clear();

    //cout<<" transforming to NA61 frame"<<endl;
    TransformToNA61Frame_HT(outNode,vdtracks,inNodeJ,inNodeS);

    Int_t NVdTracks = vdtracks.GetEntries();
    //cout<<"vdtracks: "<<NVdTracks<<endl;
    if(vdtracks.GetEntries() > 0) MatchTracks(tra,recEvent,mainVertex,vdtracks);   

    if(vdtracks.GetEntries() > NVdTracks){ // supplement Vd HT Track table in outNode 
      //cout<<"track at entry="<<NVdTracks<<" after="<<vdtracks.GetEntries()<<endl;
      // add new tracks to table
      for(Int_t i=NVdTracks; i<vdtracks.GetEntries(); i++){
	(((UVdEvent*)outNode)->GetDataTable("Vd HT Tracks"))->Add(vdtracks.At(i));
      }
    }
    
  }


  // Add offset to TPC tracks and cluster. The offset based on TPC Jura vs TPC Saleve vertexes
  //CorrectForTpcOffset(recEvent,mainVertex);     
  
 // UDataTable* tpctracks = new UDataTable("Tpc Tracks");
 // tpctracks -> SetOwner();
  
  Int_t nTracks = 0;
 
   for (VertexTrackIndexIterator vtxTrackIter = mainVertex.DaughterTracksBegin();
       vtxTrackIter != mainVertex.DaughterTracksEnd(); ++vtxTrackIter) {
    const VertexTrack& vtxTrack = recEvent.Get(*vtxTrackIter);		 
    
    if(vtxTrack.GetStatus())continue;
    
    if (!vtxTrack.HasTrack()) continue;
    
    if(fMatchArr[(*vtxTrackIter)] ==0 )continue; // don't add tracks without match
    
    const Track& track = recEvent.Get(vtxTrack.GetTrackIndex());
    
    if(track.GetStatus())continue;
    
    nTracks++;
    
    UVdTrack* tpctrack = MakeTpcTrackAsVd(tra,track);
    if(!tpctrack)continue;
    tpctracks -> Add((TObject*)tpctrack);

    
    if(fMatchArr[(*vtxTrackIter)] ==2) tpctrack -> SetCombMeth(10); // use for x direction
	
	fhMomPoints_ZX->Fill(track.GetMomentumPoint().GetZ(),track.GetMomentumPoint().GetX());
    
  }

  fTpcTracks = nTracks;
  //cout<<"TPC nTracks="<<nTracks<<endl;
   
  outNode->AddDataTable(tpctracks);
  
}

//____________________________________________________________________
void Na61VdTpcMatchingModule :: CorrectForTpcOffset(RecEvent& recEvent, const Vertex& mainVertex)
{

  for (VertexTrackIndexIterator vtxTrackIter = mainVertex.DaughterTracksBegin();
       vtxTrackIter != mainVertex.DaughterTracksEnd(); ++vtxTrackIter) {
    const VertexTrack& vtxTrack = recEvent.Get(*vtxTrackIter);		 
    
    if(vtxTrack.GetStatus())continue;
    
    if (!vtxTrack.HasTrack()) continue;
    
    Track& track = recEvent.Get(vtxTrack.GetTrackIndex());
    
    if(track.GetStatus())continue;
    
  
      

    double x = track.GetMomentumPoint().GetX();
    double y = track.GetMomentumPoint().GetY();
    double z = track.GetMomentumPoint().GetZ();

    double offx=0;
    double offy=0;
    double offz=0;
    if(x>0){ // jura care
      //offx =  fVdParams->GetJuraTpcOffset().GetX();
      //offy =  fVdParams->GetJuraTpcOffset().GetY();
      //offz =  fVdParams->GetJuraTpcOffset().GetZ();
    }else{
      //offx =  fVdParams->GetSaleveTpcOffset().GetX();
      //offy =  fVdParams->GetSaleveTpcOffset().GetY();
      //offz =  fVdParams->GetSaleveTpcOffset().GetZ();
    }

    //double x = 0; // add here offset
    //double y = 0;
    //double z = 0;

    //cout<<"x="<<x<<" "<<offx<<" "<<offy<<" "<<offz<<endl;

    const Point mPos(x+offx, y+offy, z+offz);
    track.SetMomentumPoint(mPos);

    for (ClusterIndexIterator clusterIter = track.ClustersBegin();
	 clusterIter != track.ClustersEnd(); ++clusterIter) {
      
      Cluster& cluster = recEvent.Get(*clusterIter);
      
      double xc = cluster.GetPosition().GetX(); // add here offset
      double yc = cluster.GetPosition().GetY();
      double zc = cluster.GetPosition().GetZ();
      
      const Point cPos(xc+offx,yc+offy,zc+offz);
      cluster.SetPosition(cPos);
      
    }	
        
  }
  
}

//__________________________________________________________________________________________________
void Na61VdTpcMatchingModule :: TransformToNA61Frame_HT(UEventNode* out, TObjArray& vdtracks, UEventNode* inJ, UEventNode* inS)
{

  // translation (track by track): 
  // take track by track and shift its origin. 

  // apply offsets
  UDataTable* tracktab = out -> GetDataTable(Form("Vd HT Tracks"));
  //cout<<" vdHT tracktab: "<<tracktab<<endl;
  if(!tracktab)return;
  
  for(Int_t j=0; j<tracktab->GetEntries(); j++){
    UVdTrack* track = ( UVdTrack*)tracktab->At(j);
    
    double ox = (track->Getlineb())->GetOrigin().X();
    double oy = (track->Getlineb())->GetOrigin().Y();
    double oz = (track->Getlineb())->GetOrigin().Z();
    //cout<<"oz="<<oz<<endl;
    
    (track->Getlineb())->SetOrigin(ox+fOffsetX, oy+fOffsetY, oz+fOffsetZ);
    //double x0 = track->GetXatZ_b(150.);
    //double y0 = track->GetYatZ_b(150.);
    //(track->Getlineb())->SetOrigin(x0, y0, 150.);
    
    vdtracks.Add(track);
    
  }
  
  TransformHits(inJ,fOffsetX,fOffsetY,fOffsetZ);
  TransformHits(inS,fOffsetX,fOffsetY,fOffsetZ);
  
  
}


//______________________________________________________________________________
void Na61VdTpcMatchingModule :: TransformHits(UEventNode* node, double dx, double dy, double dz)
{
	
	for(Int_t is=0;is<8;is++){
		UDataTable* hits = node->GetDataTable(Form("Hits %s",fSensorNames[is].Data()));

		if(!hits)continue;
		
		for(Int_t i=0;i<hits->GetEntries();i++){
			USensorHit* hit = (USensorHit*)hits->At(i);
			hit->SetMatchedToTPCTrack(-1);
			//cout<<"Set hit matched with TPC track to -1"<<endl;
			double x = hit->GetX();
			double y = hit->GetY();
			double z = hit->GetZ();
			
			hit->SetX(x+dx);
			hit->SetY(y+dy);
			hit->SetZ(z+dz);

			fhMomPoints_ZX -> Fill( hit->GetZ()/10., hit->GetX()/10.);
			fhMomPoints_ZX_fine -> Fill( hit->GetZ()/10., hit->GetX()/10.);
		}
	}
	
}

//______________________________________________________________________________________________________
void Na61VdTpcMatchingModule :: MatchTracks(const det::MagneticFieldTracker& tra, RecEvent& recEvent,const Vertex& mainVertex,TObjArray& vdtracks)
{
  const int acceFlagFit  = 0x000008; // clusters Used in track fit
  //const int acceFlagGold = 0x080000; // Golden clusters
  
 		  modutils::KalmanFitter kf_AK;
		  kf_AK.Init(recEvent);
  
  double off_day =  0.0;
  //double sig_day =  0.004;
  double sig_day =  0.0035;

//  double nSig1 = fNSig;         // N sigma for matching in theta
//  double nSig2 = 4.0;           // N sigma for matching in dx and dy
  double nSig1 = 4;               // N sigma for matching in theta
  double nSig2 = 6;               // N sigma for matching in dx and dy

  int N = vdtracks.GetEntries();
  //cout<<"vd tracks at entry: "<<N<<endl;
  

  for (VertexTrackIndexIterator vtxTrackIter = mainVertex.DaughterTracksBegin();
       vtxTrackIter != mainVertex.DaughterTracksEnd(); ++vtxTrackIter) {
    const VertexTrack& vtxTrack = recEvent.Get(*vtxTrackIter);     
    
    if(vtxTrack.GetStatus())continue;
    if (!vtxTrack.HasTrack()) continue;
    
    const Track& track = recEvent.Get(vtxTrack.GetTrackIndex());
    
    if(track.GetStatus())continue;
    if(track.GetCharge()==0)continue;

    fMatchArr[(*vtxTrackIter)] = 0;
        
    fhClusters_All->Fill(track.GetNumberOfClusters());
    
    if (track.GetCharge()>0)
		hMomAllPos->Fill(track.GetMomentum().GetMag());
	else
		hMomAllNeg->Fill(track.GetMomentum().GetMag());
    hMomAllTpc->Fill(track.GetMomentum().GetMag());
    hPtAllTpc->Fill(TMath::Sqrt(track.GetMomentum().GetX()*track.GetMomentum().GetX() +track.GetMomentum().GetY()*track.GetMomentum().GetY()));
//	cout<<"VtrackID="<<(int)vtxTrack.GetTrackIndex()<<"	p="<<track.GetMomentum().GetMag()<<"	pT="<<TMath::Sqrt(track.GetMomentum().GetX()*track.GetMomentum().GetX() +track.GetMomentum().GetY()*track.GetMomentum().GetY())<<
//			"	nClusters="<<track.GetNumberOfClusters()<<endl;

    //-----------------------Select good tracks
    //int VtpcClusters = track.GetNumberOfClusters(TrackConst::eVTPC1) + track.GetNumberOfClusters(TrackConst::eVTPC2);
    //if(VtpcClusters < 10)continue;
    //if(track.GetNumberOfClusters() < 25)continue;
    
    Int_t AcceClusters = 0;
    for (ClusterIndexIterator clusterIter = track.ClustersBegin();
	 clusterIter != track.ClustersEnd(); ++clusterIter) {
		    const Cluster& cluster = recEvent.Get(*clusterIter);
			if((cluster.GetStatus() & acceFlagFit) == acceFlagFit)AcceClusters++;
			//if((cluster.GetStatus() & acceFlagGold) == acceFlagGold)AcceClusters++;
			/*if((cluster.GetStatus() & acceFlag) == acceFlag){
			  cout<<(dec)<<cluster.GetStatus()<<"  "<<std::bitset<24>(cluster.GetStatus())<<endl;
			  cout<<(hex)<<acceFlag<<"  "<<cluster.GetStatus()<<" "<<(cluster.GetStatus() & acceFlag)<<endl;
			}*/
    }	
    if(AcceClusters<2) continue; // don't used track with low number of accepted clusters
  
    fhClusters_All_Acce->Fill(track.GetNumberOfClusters());

    
    Float_t Mom = track.GetMomentum().GetMag(); // check whether direction is unit vector.
    Float_t vtxMom = vtxTrack.GetMomentum().GetMag(); // check whether direction is unit vector.
    
    Float_t mcorr = vtxMom/Mom;

    //double x = track.GetMomentumPoint().GetX();
    //double y = track.GetMomentumPoint().GetY();
    //double z = track.GetMomentumPoint().GetZ();
    double x = track.GetFirstPointOnTrack().GetX();
    double y = track.GetFirstPointOnTrack().GetY();
    double z = track.GetFirstPointOnTrack().GetZ();
    
    double px = track.GetMomentum().GetX();
    double py = track.GetMomentum().GetY();
    double pz = track.GetMomentum().GetZ();
    
    double norm2 = TMath::Sqrt(pz*pz + px*px);
    
    double tpc_tany =  py/norm2;

    //track.SetCharge(track.GetCharge() * fFieldScaleFactor);
    
    const Point  startPos   =  track.GetMomentumPoint();
    //const Vector startMom   =  track.GetMomentum();
    const Vector startMom(px*mcorr,py*mcorr,pz*mcorr); //------------changed
    const double charge = track.GetCharge();
    
    
    
//----------------------------------------------------------------------     
    Point posAtTarget;
    Vector tmpMom;
    const double targetPos = -602.7*cm;
    if(!tra.TrackToZ(targetPos, charge, startPos, startMom, posAtTarget, tmpMom)) {
      cerr << "target tracking failed" << endl;
            cout<<"error in the target extrapolation: "<<targetPos<<" "<<charge<<" "<<startPos<<" "<<startMom<<" "<< posAtTarget<<" "<<tmpMom<<endl;

      continue;
    }


//    fhXY_AtTarget->Fill(posAtTarget.GetX(), posAtTarget.GetY());
//    if(TMath::Abs((Mom-vtxMom))/Mom < 0.001){
//      fMatchArr[(*vtxTrackIter)] = 1;
//      fhXY_AtTarget_clean->Fill(posAtTarget.GetX(), posAtTarget.GetY());
//    }
	double  cMom = TMath::Abs((Mom-vtxMom))/Mom; 
    fhXY_AtTarget->Fill(posAtTarget.GetX(), posAtTarget.GetY());
    if(cMom < fMomDiffVtx2){
		fMatchArr[(*vtxTrackIter)] = 1;
		if(cMom < fMomDiffVtx){
			fMatchArr[(*vtxTrackIter)] = 2;
			fhXY_AtTarget_clean->Fill(posAtTarget.GetX(), posAtTarget.GetY());
		}
    }
    
    double tdx = posAtTarget.GetX()-0.1265;
    double tsigx = 0.343;

    double tdist = TMath::Sqrt((tdx*tdx)/(tsigx*tsigx));
    fhTdistVsMomDiff -> Fill((Mom-vtxMom)/Mom, tdist);

    
  //----------------------------------------------------------------------     
    double_t startTPC = startPos.GetZ();
    const double startVD = -582*cm;
    const double frontTPC = -490*cm;
    /////////////////////////////////////////////////////////////////////////////////////////////
    // this block is just to handle diagnostic for matching, can be removed in future
    Point posAtFrontTPC,posAtVds4;
    //Vector tmpMom;
    if(!tra.TrackToZ(frontTPC, charge, startPos, startMom, posAtFrontTPC, tmpMom)) {
      cerr << "TPC tracking failed" << endl;
      continue;
    }
    if(!tra.TrackToZ(startVD, charge, startPos, startMom, posAtVds4, tmpMom)) {
      cerr << "VD tracking failed" << endl;
      continue;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////

    //const double zEnd = z; // just for test
    const double zEnd = startTPC; // just for test
    //const double zEnd = startVD + (3./4.)*(startTPC-startVD); 
    //cout<<"zEnd="<<zEnd<<" startTPC="<<startTPC<<endl;

    Point endPos;
    Vector endMom;
    
    if(!tra.TrackToZ(zEnd, charge, startPos, startMom, endPos, endMom)) {
      cerr << "zEnd tracking failed" << endl;
      continue;
    }else{
      //cout << "--------> new TPC track,  charge="<<charge<<"  ************************************************************* "<<endl;
      //ShowPoint(endPos);
      //ShowMomentum(endMom);
    }
    //cout<<"zend="<<zEnd<<" endPos="<<endPos<<" track.GetFirstPointOnTrack()="<<track.GetFirstPointOnTrack()<<endl;
    
    //Int_t ii;
    //cin>>ii;
    std::map<int,double> dyArray;
    std::map<int,double> dxArray;
    std::map<int,double> dayArray;
    std::map<int,double> dpxArray;
    std::map<int,double> dpyArray;
    std::map<int,double> dpzArray;
    std::map<int,double> dd2Array;
    std::map<int,int> iArray;
//    double dyArray[50];
//    double dxArray[50];
//    double dayArray[50];
//    double dpxArray[50];
//    double dpyArray[50];
//    double dpzArray[50];
//    double dd2Array[50];
//    Int_t iArray[50];
    
    Int_t imatch = 0;

    for(Int_t i=0; i<N; i++){ // loop over vd tracks
      UVdTrack* vdtrack = (UVdTrack*)vdtracks.At(i);
      
      double norm = TMath::Sqrt(vdtrack->GetDZ_b()*vdtrack->GetDZ_b() + vdtrack->GetDX_b()*vdtrack->GetDX_b());      
      double vd_tany =  vdtrack->GetDY_b()/norm;

      double day =  tpc_tany-vd_tany;

      fhTanY->Fill(day);
      fhTanYvsY->Fill(posAtFrontTPC.GetY(),day);
     
      if(TMath::Abs(day-off_day) > nSig1*sig_day) continue; // cuts on y track slope

     
      const Point  vdStartPos(vdtrack->GetX_b()/10., vdtrack->GetY_b()/10.,vdtrack->GetZ_b()/10.);
      const Vector vdStartMom(vdtrack->GetDX_b()*vtxMom, vdtrack->GetDY_b()*vtxMom,vdtrack->GetDZ_b()*vtxMom); 
      
      double x_vd = vdStartPos.GetX();

      Point  vdEndPos;
      Vector vdEndMom;
      
      if(!tra.TrackToZ(zEnd, charge, vdStartPos, vdStartMom, vdEndPos, vdEndMom)) {
		cerr << "zEnd vd track: tracking failed" << endl;
		continue;
      }else{
		double dx =  (endPos.GetX()/mm)  - (vdEndPos.GetX()/mm);
		double dy =  (endPos.GetY()/mm)  - (vdEndPos.GetY()/mm);
		double dz =  (endPos.GetZ()/mm)  - (vdEndPos.GetZ()/mm);

		double dpx = (endMom.GetX()/GeV) - (vdEndMom.GetX()/GeV);
		double dpy = (endMom.GetY()/GeV) - (vdEndMom.GetY()/GeV);
		double dpz = (endMom.GetZ()/GeV) - (vdEndMom.GetZ()/GeV);

		fhPosMatch[0]->Fill(dx); 
		fhPosMatch[1]->Fill(dy); 
		fhPosMatch[2]->Fill(dz);

		fhMomMatch[0]->Fill(dpx); 
		fhMomMatch[1]->Fill(dpy); 
		fhMomMatch[2]->Fill(dpz);


		double off_dx=0;
		double sig_dx=0;
		if(!fVdParams->GetDxMatchParams(z,x,x_vd,off_dx,sig_dx)){
		  cout<<"Warning: dx matching paramater not setup correctly!!! z="<<z<<" x="<<x<<" x_vd="<<x_vd<<endl; 
		  continue;
		}
		double off_dy=0;
		double sig_dy=0;
		if(!fVdParams->GetDyMatchParams(z,y,x,x_vd,off_dy,sig_dy)){
		  cout<<"Warning: dy matching paramater not setup correctly!!! z="<<z<<" y="<<y<<" x="<<x<<" x_vd="<<x_vd<<endl; 
		  continue;
		}
		//if(z<-220 && x>0 && x_vd<0){
		//cout<<"z="<<z<<" x="<<x<<" x_vd="<<x_vd<<endl;
		//cout<<"off_dy="<<off_dy<<" sig_dy="<<sig_dy<<" "<<TMath::Abs(dy-off_dy)<<endl;
		//}
		
		if(x>0){ // jura part
		  fhXvsPx_J->Fill(dx,dpx);
		  fhXvsPz_J->Fill(dx,dpz);
				
		  if(TMath::Abs(dy-off_dy) < 3*sig_dy){
			if(x_vd > 0)fhX_JJ_vsZ_cuty->Fill(z,dx); else fhX_JS_vsZ_cuty->Fill(z,dx);	    
		  }
		  if(TMath::Abs(dx-off_dx) < 3*sig_dx){

			if(x_vd > 0){
			  fhY_JJ_vsY_cutx->Fill(y,dy);  
			  if(z<-240) fhY_JJ_vsY_cutx_r1->Fill(y,dy); else fhY_JJ_vsY_cutx_r2->Fill(y,dy); 
			}else{
			  fhY_JS_vsY_cutx->Fill(y,dy);
			  if(z<-240) fhY_JS_vsY_cutx_r1->Fill(y,dy); else fhY_JS_vsY_cutx_r2->Fill(y,dy); 
			}

		  }
		  
		}else{  // saleve part
		  
		  fhXvsPx_S->Fill(dx,dpx);
		  fhXvsPz_S->Fill(dx,dpz);
		  
		  if(TMath::Abs(dy-off_dy) < 3*sig_dy){
			if(x_vd > 0)fhX_SJ_vsZ_cuty->Fill(z,dx); else fhX_SS_vsZ_cuty->Fill(z,dx);
		  }

		  if(TMath::Abs(dx-off_dx) < 3*sig_dx){

			if(x_vd > 0){
			  fhY_SJ_vsY_cutx->Fill(y,dy); 
			  if(z<-240) fhY_SJ_vsY_cutx_r1->Fill(y,dy); else fhY_SJ_vsY_cutx_r2->Fill(y,dy); 
			}else{
			  fhY_SS_vsY_cutx->Fill(y,dy);
			  if(z<-240) fhY_SS_vsY_cutx_r1->Fill(y,dy); else fhY_SS_vsY_cutx_r2->Fill(y,dy); 
			}
	 
		  }
		}
		

		if(TMath::Abs(dx-off_dx) < 3*sig_dx){
		  fhMomMatch_cutx[0]->Fill(dpx); 
		  fhMomMatch_cutx[1]->Fill(dpy); 
		  fhMomMatch_cutx[2]->Fill(dpz);
		}

		if(TMath::Abs(dy-off_dy) < 3*sig_dy){
		  fhMomMatch_cuty[0]->Fill(dpx); 
		  fhMomMatch_cuty[1]->Fill(dpy); 
		  fhMomMatch_cuty[2]->Fill(dpz);
		}


		double dd2 = 99999;
		if((sig_dx>0) && (sig_dy>0)){ 
		  dd2 =((dx-off_dx)*(dx-off_dx))/(sig_dx*sig_dx)  + ((dy-off_dy)*(dy-off_dy))/(sig_dy*sig_dy);
		}
		//double dd2 = ((dx-off_dx)*(dx-off_dx))/(sig_dx*sig_dx) 
		//+ ((dpz-off_dpz)*(dpz-off_dpz))/(sig_dpz*sig_dpz);
		
		if(dd2 < nSig2*nSig2){ // accepted matchings
		  fhMatch_TanYvsY_acce->Fill(day,dy);
		  fhTanY_acce->Fill(day);
		  fhTanYvsY_acce->Fill(posAtFrontTPC.GetY(),day);

		  fhMomMatch_cutxy[0]->Fill(dpx); 
		  fhMomMatch_cutxy[1]->Fill(dpy); 
		  fhMomMatch_cutxy[2]->Fill(dpz);
		  dxArray[imatch] = dx;
		  dyArray[imatch] = dy;
		  dayArray[imatch] = day;
		  dpxArray[imatch] = dpx;
		  dpyArray[imatch] = dpy;
		  dpzArray[imatch] = dpz;
		  dd2Array[imatch] = dd2;
		  iArray[imatch] = i;
		  imatch++;

		}
		//ShowPoint(vdStartPos);
		//ShowMomentum(vdStartMom);
		//ShowPoint(vdEndPos);
		//ShowMomentum(vdEndMom);
      }
      
    }
    

    /////////////////////////////// just for matching diagnostics //////////////////////////////
    if(imatch==0){
      fhXY_fTPC_NoMatch->Fill(posAtFrontTPC.GetX(),posAtFrontTPC.GetY());
      fhXY_Vds4_NoMatch->Fill(posAtVds4.GetX(),posAtVds4.GetY());
    }else{
      fhXY_fTPC_Match->Fill(posAtFrontTPC.GetX(),posAtFrontTPC.GetY());
      fhXY_Vds4_Match->Fill(posAtVds4.GetX(),posAtVds4.GetY());
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
	//if (imatch==0) {     cout<<"imatch="<<imatch<<"   Mom="<<track.GetMomentum().GetMag()/GeV<<endl;    continue; }
    double dd2_min=10000;
    Int_t im_min=-1;
    if(imatch>0){
      for(Int_t im=0;im<imatch;im++){ //select the best matching
		//UVdTrack* vdtrack = (UVdTrack*)vdtracks.At(iArray[im]);
		//cout<<"im="<<im<<" line: "<<vdtrack->Getlineb()<<endl;
		/*
		Int_t ind1 = vdtrack->GetHitIndexOnStation(0);
		Int_t ind2 = vdtrack->GetHitIndexOnStation(1);
		Int_t ind3 = vdtrack->GetHitIndexOnStation(2);
		Int_t ind4 = vdtrack->GetHitIndexOnStation(3);
		Int_t tab1 = vdtrack->GetTabIndexOnStation(0);
		Int_t tab2 = vdtrack->GetTabIndexOnStation(1);
		Int_t tab3 = vdtrack->GetTabIndexOnStation(2);
		Int_t tab4 = vdtrack->GetTabIndexOnStation(3);
		*/
		//cout<<"tab1="<<tab1<<" tab2="<<tab2<<" tab3="<<tab3<<" tab4="<<tab4<<" "<<dd2Array[im]<<endl;
		//cout<<"ind1="<<ind1<<" ind2="<<ind2<<" ind3="<<ind3<<" ind4="<<ind4<<endl;

		if(dd2Array[im]<dd2_min){
		  im_min = im;
		  dd2_min=dd2Array[im];
		}
      }


      //Int_t iii;
      //if(imatch==17)cin>>iii;
      fhMatch_TanYvsY_best->Fill(dayArray[im_min], dyArray[im_min]);
      fhMomMatch_best[0]->Fill(dpxArray[im_min]);
      fhMomMatch_best[1]->Fill(dpyArray[im_min]);
      fhMomMatch_best[2]->Fill(dpzArray[im_min]);

      fhPosMatch_best[0]->Fill(dxArray[im_min]); 
      fhPosMatch_best[1]->Fill(dyArray[im_min]); 

      UVdTrack* vdtrack_best = (UVdTrack*)vdtracks.At(iArray[im_min]);
      // set momentum and matching flag
      //momentum at front vd track

      /////////////////////////////////////  SECOND MATCHING CASE  ////////////////////////////////////////////////////
      Bool_t second_match=kFALSE;
      if(vdtrack_best->GetTpcMatchingFlag()){
		second_match=kTRUE;
		//cout<<"vdtrack_best - second match "<<vdtrack_best->GetTpcMatchingFlag()<<" imatch="<<imatch<<endl;
		//cout<<"Momentum = "<<vdtrack_best->GetMomentum()<<" Kf Momentum ="<<vdtrack_best->GetKfMomentum()<<endl;
		// vdtrack_best already wat matched to TPC track. We treat the socong match as a new global track, so
		// we create new vdtrack_best using copy constructor and used  it as a seed of of new global track
		vdtrack_best = new UVdTrack(vdtrack_best);

		vdtracks.Add(vdtrack_best); // vd track table is increased, should not break the loop
      }
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      double Mom =  track.GetMomentum().GetMag();
      //double vtxMom = vtxTrack.GetMomentum().GetMag(); // check whether direction is unit vector.
      double px_f = vdtrack_best->GetDX_f()*Mom;
      double py_f = vdtrack_best->GetDY_f()*Mom;
      double pz_f = vdtrack_best->GetDZ_f()*Mom;  

      Int_t charge = (Int_t)track.GetCharge();
      vdtrack_best->SetMomentum(px_f,py_f,pz_f); //tmp, will be corrected below
      
      vdtrack_best->SetTpcMatchingFlag(1);
      if(second_match)vdtrack_best->SetTpcMatchingFlag(2);
      
      if(vdtrack_best->GetTagForVtx()==1){
			if(fMatchArr[(*vtxTrackIter)] !=2 ) fMatchArr[(*vtxTrackIter)] = 1;
      }
      //fMatchArr[(*vtxTrackIter)] = 1;
      //cout<<"vtxTrack: "<<(*vtxTrackIter)<<"  "<< fMatchArr[(*vtxTrackIter)]<<endl;


      vdtrack_best->SetTpcCharge(charge);
      vdtrack_best->SetTpcTrackIndex((Int_t)vtxTrack.GetTrackIndex());
      fhChargeCorr->Fill(charge,vdtrack_best->GetCharge());
      fNumberOfVdTpcTracks++;
      //cout<<"best vd charge="<<vdtrack_best->GetCharge()<<" tpc charge="<<charge<<"  p="<<Mom<<endl;

      //////////////////////////// re-fit using kalman fitter ///////////////////////////////////////
       vector<evt::Index<evt::rec::Cluster> > acceptedClustersFitAK;

      /////// first add VD clusters(=hits)	
      for(Int_t ih=0;ih<4;ih++){	
		    USensorHit* hit = (USensorHit*)vdtrack_best->GetHitIdAtStation(ih);
			if(!hit)continue;
			
			rec::Cluster& vdcluster = recEvent.Make<rec::Cluster>();
			vdcluster.SetPosition(Point(hit->GetX()*0.1, hit->GetY()*0.1, hit->GetZ()*0.1));
			vdcluster.SetPositionUncertainty(rec::ClusterConst::eX, 5*micrometer);
			vdcluster.SetPositionUncertainty(rec::ClusterConst::eY, 5*micrometer);
						
			acceptedClustersFitAK.push_back(vdcluster.GetIndex());

      }
      //////// add TPC clusters
      //const Int_t acceFlagFit  = 0x000008; // clusters Used in track fit
	  //const Int_t acceFlagGold = 0x080000; // Golden clusters
      Int_t TpcClustersFit = 0;
      for (ClusterIndexIterator clusterIter = track.ClustersBegin();
	   clusterIter != track.ClustersEnd(); ++clusterIter) {
	
		const Cluster& cluster_TPC = recEvent.Get(*clusterIter);
		
		// arrays that include VD clusters ("global tracks")			
		if((cluster_TPC.GetStatus() & acceFlagFit) == acceFlagFit)TpcClustersFit++;
		if((cluster_TPC.GetStatus() & acceFlagFit) == acceFlagFit)acceptedClustersFitAK.push_back(cluster_TPC.GetIndex()); 

      }	

      vdtrack_best->SetNumberOfTpcClusters(TpcClustersFit);
      //cout<<" TpcClustersFit="<<TpcClustersFit<<endl;

      //cout<<clusterIndices.size()<<"  trackInd="<<vtxTrack.GetTrackIndex()<<endl;
      //cout<<"z1="<<vertexVD.GetZ()*0.1<<endl;
      //for(Int_t ic = 0; ic<clusterIndices.size();ic++){
      //const Cluster&  cluster = recEvent.Get(clusterIndices.at(ic));
      //cout<<"ic="<<ic<<" z="<<cluster.GetPosition().GetZ()<<endl;
      //}
      
          
	     
      //refit TPC track including VD clusters    
	  TrackLocalParameters trackTpcParametersFitAK;
      vdtrack_best->SetFitKf(kf_AK.KalmanFilter(acceptedClustersFitAK, trackTpcParametersFitAK));
      const Point recPositionAK(trackTpcParametersFitAK[TrackLocalParameters::eX],
			      trackTpcParametersFitAK[TrackLocalParameters::eY],
			      trackTpcParametersFitAK[TrackLocalParameters::eZ]);
      const Vector recMomentumAK(trackTpcParametersFitAK[TrackLocalParameters::ePX],
			       trackTpcParametersFitAK[TrackLocalParameters::ePY],
			       trackTpcParametersFitAK[TrackLocalParameters::ePZ]);


      double Mom_kf = recMomentumAK.GetMag();  //AK refitting
      //cout<<"Mom test: Mom_kf="<<Mom_kf<<endl;
      double px_kf = vdtrack_best->GetDX_f()*Mom_kf;
      double py_kf = vdtrack_best->GetDY_f()*Mom_kf;
      double pz_kf = vdtrack_best->GetDZ_f()*Mom_kf;
      
      fhMomVsKfMom->Fill(Mom,Mom_kf);
      vdtrack_best->SetKfMomentum(px_kf,py_kf,pz_kf); //pawel's kf mom
      
      vdtrack_best->SetMomentumKf(utl::Vector(px_kf, py_kf, pz_kf)); //shine kf mom
      vdtrack_best->SetPositionKf(recPositionAK);
      vdtrack_best->SetChargeKf(trackTpcParametersFitAK.GetCharge());

	  fhMomDiff->Fill((Mom_kf-Mom)/Mom); // relative change in momentum


      fhPt->Fill(sqrt(px_kf*px_kf + py_kf*py_kf));	  
	  
     //  Float_t MomVds = vtxMom;//-------------------------------second option - use vtxTrack mom
      //Float_t px_f = vdtrack_best->GetDX_f()*MomVds;
      //Float_t py_f = vdtrack_best->GetDY_f()*MomVds;
      //Float_t pz_f = vdtrack_best->GetDZ_f()*MomVds;  
      //vdtrack_best->SetMomentum(px_f,py_f,pz_f);
      
      
      //extrapolate closer to primary vertex
      //and refit TPC track including VD cluaters
/*      const Vector vdStartMomKf(px_kf,py_kf,pz_kf); 
      const Point vdStartPosKf((vdtrack_best->GetX_f()+fOffsetX)*0.1,(vdtrack_best->GetY_f()+fOffsetY)*0.1,(vdtrack_best->GetZ_f()+fOffsetZ)*0.1); 
      Point  vdEndPosKf;
      Vector vdEndMomKf;
      const double zEnd_Kf = -602*cm; //target is at ~ -603
      if(!tra.TrackToZ(zEnd_Kf, charge, vdStartPosKf, vdStartMomKf, vdEndPosKf, vdEndMomKf)) {
		cerr << "vd KF track: tracking failed" << endl;
		continue;
      }else{
	
	//double px_kf_kf = vdEndMomKf.GetX();
        //double py_kf_kf = vdEndMomKf.GetY();
        //double pz_kf_kf = vdEndMomKf.GetZ();
      
        vdtrack_best->SetMomentumKf(vdEndMomKf); //shine kf mom after extrapolation
        vdtrack_best->SetPositionKf(vdEndPosKf);
        vdtrack_best->SetChargeKf(trackParameters.GetCharge());
        
        //cout<<"fOffsetZ="<<fOffsetZ<<endl;
	//    cout<<"before:"<<vdtrack_best->GetX_f()<<" "<<vdtrack_best->GetY_f()<<" "<<vdtrack_best->GetZ_f()<<endl;
        (vdtrack_best->Getlinef())->SetOrigin(vdEndPosKf.GetX()*10-fOffsetX,vdEndPosKf.GetY()*10-fOffsetY,vdEndPosKf.GetZ()*10-fOffsetZ);
	    (vdtrack_best->Getlinef())->SetDirection(vdEndMomKf.GetX(),vdEndMomKf.GetY(),vdEndMomKf.GetZ());
	    //double ox = vdtrack_best->GetXatZ_f(0.0);
	    //double oy = vdtrack_best->GetYatZ_f(0.0);
	    //(vdtrack_best->Getlinef())->SetOrigin(ox,oy,0.0);
	  //  cout<<"after:"<<vdtrack_best->GetX_f()<<" "<<vdtrack_best->GetY_f()<<" "<<vdtrack_best->GetZ_f()<<endl;

        //if(second_match)cout<<"Mom="<<Mom<<"  KfMom="<<Mom_kf<<endl;
        fhMomDiff_kf->Fill((vdEndMomKf.GetMag()-Mom_kf)/Mom_kf); // relative change in momentum
        fhMomVsKfMom_kf->Fill(Mom_kf,vdEndMomKf.GetMag());

		//cout<<vdStartMomKf<<" "<<vdEndMomKf<<endl;
	  }
 */     
      

      // dEdx block
      double dEdx_all = track.GetEnergyDeposit(TrackConst::eAll);
      Int_t dEdxClusters_all = track.GetNumberOfdEdXClusters(TrackConst::eAll);
      fhdEdx_all->Fill(Mom_kf,dEdx_all);
      if(charge==1)fhdEdxPos_all->Fill(Mom_kf,dEdx_all);
      else fhdEdxNeg_all->Fill(Mom_kf,dEdx_all);
      fhdEdxVsClusters_all->Fill(dEdxClusters_all,dEdx_all);
      if((Mom_kf>4) && (Mom_kf<6)) fhdEdxVsClustersCut_all->Fill(dEdxClusters_all,dEdx_all);
      fhClustersVsMom_all->Fill(Mom_kf,dEdxClusters_all);
      //vdtrack_best->SetdEdx_all(dEdx_all);
      //vdtrack_best->SetdEdxClusters_all(dEdxClusters_all);

      double dEdx_mtpc = track.GetEnergyDeposit(TrackConst::eMTPC);
      Int_t dEdxClusters_mtpc = track.GetNumberOfdEdXClusters(TrackConst::eMTPC);
      fhdEdx_mtpc->Fill(Mom_kf,dEdx_mtpc);
      if(charge==1)fhdEdxPos_mtpc->Fill(Mom_kf,dEdx_mtpc);
      else fhdEdxNeg_mtpc->Fill(Mom_kf,dEdx_all);
      fhdEdxVsClusters_mtpc->Fill(dEdxClusters_mtpc,dEdx_mtpc);
      if((Mom_kf>4) && (Mom_kf<6)) fhdEdxVsClustersCut_mtpc->Fill(dEdxClusters_mtpc,dEdx_mtpc);
      //vdtrack_best->SetdEdx_mtpc(dEdx_mtpc);
      //vdtrack_best->SetdEdxClusters_mtpc(dEdxClusters_mtpc);

      vdtrack_best->ActivateTLV();
      vdtrack_best->ActivateTLV_Kf();

      //cout<<"dEdx_mtpc="<< dEdx<<" clusters: "<<dEdxClusters<<" / "<<TPCClusters<<" / "<<track.GetNumberOfClusters(TrackConst::eAll)<<endl;
      if(TpcClustersFit>15)fhMomVsKfMom_acce->Fill(Mom,Mom_kf);
      //if(TPCClusters<15)cout<<"Mom="<<Mom<<" Mom_kf="<<Mom_kf<<" "<<TPCClusters<<endl;
      
    }  
    
  }  



  Int_t VdTracksWithMatch = 0;
  for(Int_t i=0; i<vdtracks.GetEntries(); i++){
    UVdTrack* vdtrack = (UVdTrack*)vdtracks.At(i);
    double z_vds4 = -5830.;
    double x = vdtrack->GetXatZ_b(z_vds4); // x in mm
    double y = vdtrack->GetYatZ_b(z_vds4); // y in mm
    //cout<<x<<" "<<y<<endl;
    fhXY_Vds4_All->Fill(0.1*x,0.1*y); // fill in cm
    fhCurvature_All->Fill(vdtrack->GetCurvature());
    if(vdtrack->GetTpcMatchingFlag()){
      VdTracksWithMatch++;
      fhXY_Vds4_WithMatch->Fill(0.1*x,0.1*y); // fill in cm
      fhCurvature_WithMatch->Fill(vdtrack->GetCurvature());
      fhClusters_WithMatch->Fill(vdtrack->GetNumberOfTpcClusters());
    }

    // search for correlation between track with and without match
    if(!vdtrack->GetTpcMatchingFlag())continue;

    double cur1 = vdtrack->GetCurvature();
    for(Int_t j=i+1; j<vdtracks.GetEntries(); j++){
      UVdTrack* vdtrack2 = (UVdTrack*)vdtracks.At(j);
      double cur2 = vdtrack2->GetCurvature();
      double x2 = vdtrack2->GetXatZ_b(z_vds4); // x in mm
      double y2 = vdtrack2->GetYatZ_b(z_vds4); // y in mm
      if(vdtrack2->GetTpcMatchingFlag()){
	fhCurvature_m_vs_m->Fill(cur1,cur2);
	fhXOnVds4_m_vs_m->Fill(x,x2);
	fhYOnVds4_m_vs_m->Fill(y,y2);
      }else{
	fhCurvature_m_vs_nm->Fill(cur1,cur2);
	fhXOnVds4_m_vs_nm->Fill(x,x2);
	fhYOnVds4_m_vs_nm->Fill(y,y2);
      }
    }
  }
  

  fhTracksVsTracksWithMatch->Fill(vdtracks.GetEntries(),VdTracksWithMatch);

}

//____________________________________________________________________________________________________________
UVdTrack* Na61VdTpcMatchingModule::MakeTpcTrackAsVd(const det::MagneticFieldTracker& tra, const evt::rec::Track& track)
{
  
  double ze = fOffsetZ - 47; // I think this is target in NA61/SHINE frame 
  //double ze = fOffsetZ; // I think this is target in NA61/SHINE frame 
	//cout<<"ze="<<ze<<" fOffsetZ="<<fOffsetZ<<endl;
  const Point  startPos	 =	track.GetMomentumPoint();
  const Vector startMom	 =	track.GetMomentum();
  const double charge    =      track.GetCharge();
		

  double offx=0;
  double offy=0;
  double offz=0;
  
  /* don't see improvement - need more focus to understand (Tpc vertex from fitting?)  
  if(startPos.GetX()>0){ // jura care
    offx =  fVdParams->GetJuraTpcOffset().GetX();
    offy =  fVdParams->GetJuraTpcOffset().GetY();
    offz =  fVdParams->GetJuraTpcOffset().GetZ();
  }else{
    offx =  fVdParams->GetSaleveTpcOffset().GetX();
    offy =  fVdParams->GetSaleveTpcOffset().GetY();
    offz =  fVdParams->GetSaleveTpcOffset().GetZ();
  }
  */

  Point  endPos;
  Vector endMom;

  // extrapolate to vicinity of target  
 // if(!tra.TrackToZ(ze*0.1, charge, startPos, startMom, endPos, endMom))
 //   return 0;
  if (!tra.TrackToZ(ze*0.1, charge, startPos, startMom, endPos, endMom)) {
		     cout<<"ze*0.1 Fitting failed for TPC tracks"<<endl;
		     cout<<"error in the ze*0.1 extrapolation: "<<ze*0.1<<" "<<charge<<" "<<startPos<<" "<<startMom<<" "<< endPos<<" "<<endMom<<endl;
		//continue;
		return 0;
   }
  
  Vector3D origin(endPos.GetX(),endPos.GetY(),endPos.GetZ());
  Vector3D direction(endMom.GetX(),endMom.GetY(),endMom.GetZ());

  UVdTrack* vtr  = new UVdTrack(origin,direction);
  vtr -> SetTrackID(1); // HT track 
  vtr -> MarkForRemoval(kFALSE);	  
  vtr -> SetCharge((Int_t)charge);

  //cout<<"x = "<<startPos.GetX()<<" "<<offx<<" "<<offy<<" "<<offz<<endl;
  //cout<<"endPos.GetZ()="<<endPos.GetZ()<<"    "<<endPos.GetX()<<"  "<<endPos.GetY()<<endl;
  (vtr->Getlinef())->SetOrigin(endPos.GetX()+offx,endPos.GetY()+offy,endPos.GetZ()+offz);
  (vtr->Getlinef())->SetDirection(endMom.GetX(),endMom.GetY(),endMom.GetZ());
  double ox = vtr->GetXatZ_f(0.0);
  double oy = vtr->GetYatZ_f(0.0);
  (vtr->Getlinef())->SetOrigin(ox,oy,0.0);



  // in principal we dont need back like for TPC but let us set it for the copleteness 
  (vtr->Getlineb())->SetOrigin(startPos.GetX(),startPos.GetY(),startPos.GetZ());
  (vtr->Getlineb())->SetDirection(startMom.GetX(),startMom.GetY(),startMom.GetZ());

  return vtr;

}


//_____________________________________________________________
void Na61VdTpcMatchingModule::End()
{
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61VdTpcMatchingModule::Finish()
{
  // Job-level finalisation
  SetState(kFinish);
  
}

//____________________________________________________________________
void Na61VdTpcMatchingModule ::Print(Option_t* option) const
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



