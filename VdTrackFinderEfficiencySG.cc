#include "VdTrackFinderEfficiencySG.h"
#include "Na61VdParametersManager.h"
#include "Na61VdParameters.h"
#include "Na61ArmParameters.h"
#include "Na61VdTrackingInitModule.h"
#include "Na61PrimaryVertexRecoModule.h"
#include "Na61VdTrackingHTModule.h"
#include "Na61PrimaryVertexRecoHTModule.h"
#include "Na61VdTrackingHTPackage.h"
#include "Na61VdTpcMatchingModule.h"
	#include "Na61VdEfficiencyModule.h"
	#include "Na61PrimaryVertexRecoPostModule.h"
	#include "Na61VdTrackingPostModule.h"
#include "UDataTable.h"
#include "USensorHit.h"
#include <fwk/CentralConfig.h>
#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>
#include <evt/Event.h>
#include <det/MagneticField.h>
#include <det/Detector.h>
#include <evt/RecEvent.h>
#include <evt/proc/Silicon.h>
#include <evt/raw/Silicon.h>
#include <evt/IndexedObjectLinker.h>
#include <sstream>


namespace VdTrackFinderEfficiencySG {
fwk::VModule::EResultFlag VdTrackFinderEfficiencySG::Init() {
  utl::Branch topBranch = fwk::CentralConfig::GetInstance().GetTopBranch("VdReconstructionSG");
  InitVerbosity(topBranch);

  prevRun = 0;
  productionMode=0;//should be 1 for productions and 0 for calculating efficiencies and parameters tuning
  return eSuccess;
}

//____________________________________________________________________________________________________________________________________

fwk::VModule::EResultFlag VdTrackFinderEfficiencySG::Process(evt::Event& e, const utl::AttributeMap&) {
  const evt::EventHeader& eventHeader = e.GetEventHeader();
  if (!eventHeader.IsInitialized()) {
    WARNING("No EventHeader is present in Event. This is a serious problem!!! PLEASE DO CHECK WHAT IS GOING ON IF THIS EFFECT IS NOT EXPECTED!!!");
    return eFailure;
  }
  det::Detector::GetInstance().Update(eventHeader.GetTime(), eventHeader.GetRunNumber());

  if (eventHeader.GetRunNumber() != prevRun) {
    Begin(eventHeader.GetRunNumber());
    prevRun = eventHeader.GetRunNumber();
  }


  const std::string nameStandalone[8] = {"Hits Vds1_0", "Hits Vds2_0", "Hits Vds3_0", "Hits Vds3_1", "Hits Vds4_0", "Hits Vds4_1", "Hits Vds4_2", "Hits Vds4_3"};
  const std::string nameShine[2][8] = {{"vd1_s", "vd2_s", "vd3_s_d", "vd3_s_u", "vd4_s1_d", "vd4_s1_u", "vd4_s2_d", "vd4_s2_u"}, {"vd1_j", "vd2_j", "vd3_j_d", "vd3_j_u", "vd4_j1_d", "vd4_j1_u", "vd4_j2_d", "vd4_j2_u"}};

  // Create events
  UVdEvent event;
  event.SetOwner();
  UVdEvent eventArm[2];

  // Loop over Saleve and Jura arms
  for (unsigned char arm = 0; arm < 2; arm++) {
    eventArm[arm].SetOwner();

    // Loop over senors of each arm
    for (unsigned char armSensor = 0; armSensor < 8; armSensor++) {
      UDataTable& hits = *new UDataTable(nameStandalone[armSensor].c_str());
      eventArm[arm].AddDataTable(&hits);

      const det::Silicon &detSilicon = det::Detector::GetInstance().GetSilicon();
      const det::SiliconSensor &detSensor =
        detSilicon.GetSensor(detSilicon.SensorIdToSubdetectorNumber(nameShine[arm][armSensor]));

      if (!e.GetProcEvent().GetSilicon().HasSensor(detSensor.GetSubdetectorNumber())) continue;
      const evt::proc::SiliconSensor& sensor = e.GetProcEvent().GetSilicon().GetSensor(detSensor.GetSubdetectorNumber());

      // Loop over clusters of each sensor (of each arm) and add them to list of hits
      for (auto cluster = sensor.GetClusters().begin(); cluster != sensor.GetClusters().end(); cluster++) {
        if (!e.GetRecEvent().Has(*cluster)) continue;
        if (e.GetProcEvent().GetSilicon().GetClusterLocalPositions().find(e.GetRecEvent().Get(*cluster).GetIndex()) == e.GetProcEvent().GetSilicon().GetClusterLocalPositions().end()) continue;
        const utl::Point& localPosition = e.GetProcEvent().GetSilicon().GetClusterLocalPositions()[e.GetRecEvent().Get(*cluster).GetIndex()];
        USensorHit& hit = *new USensorHit(localPosition.GetX() * 10, localPosition.GetY() * 10, localPosition.GetZ() * 10, e.GetRecEvent().Get(*cluster).GetNumberOfPixels());  // Positions in mm

        hit.SetSensorName(nameStandalone[armSensor].c_str());
        // XXX YOU CAN ONLY STORE 8-BIT TYPE THERE!!!
        //hit.SetClusterLine(e.GetRecEvent().Get(*cluster).GetFirstCoordinateNumber());
        //hit.SetClusterColumn(e.GetRecEvent().Get(*cluster).GetSecondCoordinateNumber());
        hit.SetIndex(e.GetRecEvent().Get(*cluster).GetIndex());

        hits.Add(&hit);
      }
    }

    // Call event processing fuctions of stand-alone modules
    na61VdTrackingInitModule[arm]->Event(&eventArm[arm], &eventArm[arm]);
    if (!eventArm[arm].GetPrimaryVertexStatus()) eventArm[arm].SetPrimaryVertex(-11111, -11111, -11111);

  }
  na61PrimaryVertexRecoModule->Event(&eventArm[1], &eventArm[0]);

  event.SetRunNumber(eventArm[1].GetRunNumber());
  event.SetEventNumber(eventArm[1].GetEventNumber());
  event.SetPrimaryVertexStatus(0);
  if (eventArm[0].GetPrimaryVertexStatus() || eventArm[1].GetPrimaryVertexStatus()) {
    event.SetPrimaryVertex(na61PrimaryVertexRecoModule->GetPrimaryVertex().GetX(), na61PrimaryVertexRecoModule->GetPrimaryVertex().GetY(), na61PrimaryVertexRecoModule->GetPrimaryVertex().GetZ());
    if (abs(na61PrimaryVertexRecoModule->GetPrimaryVertex().GetZ() + 47.) < 3.) event.SetPrimaryVertexStatus(1);
  }

  //run only if ProductionMode(0)! (for efficiencies)
  WARNING("-----productionMode=0----- calculating efficiencies");
  for (unsigned char arm = 0; arm < 2; arm++) {
	vdtrackingPostModule[arm]->Event(&eventArm[arm], &eventArm[arm]);
  }
  primaryvtxpostModule->Event(&eventArm[1], &eventArm[0], &event);
  for (unsigned char arm = 0; arm < 2; arm++) {
	vdefficModule[arm]->Event(&eventArm[arm], &eventArm[arm]);
  }
  
  if (!e.GetRecEvent().HasMainVertex())
	WARNING("No Main Vertex! Check this out!");
  else{
	na61VdTpcMatchingModule->Event(det::Detector::GetInstance().GetMagneticFieldTracker(), e, &eventArm[1], &eventArm[0], &event);
  	na61PrimaryVertexRecoHTModuleTPC->Event(&event);
  }



  UDataTable* tracks = event.GetDataTable("Vd HT Tracks");
  if (!tracks) return eSuccess;
  for (int t = 0; t < tracks->GetEntries(); t++) {
    UVdTrack* track = (UVdTrack*)tracks->At(t);
    if (!track) continue;

    const unsigned char match = track->GetTpcMatchingFlag() && e.GetRecEvent().Has(evt::Index<evt::rec::Track>(track->GetTpcTrackIndex()));

    evt::rec::Track& trackShine = match ? e.GetRecEvent().Get(evt::Index<evt::rec::Track>(track->GetTpcTrackIndex())) : e.GetRecEvent().Make<evt::rec::Track>();
    e.GetProcEvent().GetSilicon().GetTracks().insert(trackShine.GetIndex());

    unsigned int counter = 0;
    for (unsigned char s = 0; s < 4; s++) {
      USensorHit* hit = (USensorHit*)track->GetHitIdAtStation(s);
      if (!hit) continue;
      const evt::Index<evt::rec::Cluster> cluster(hit->GetIndex());
      if (!e.GetRecEvent().Has(cluster)) continue;
      try {
        evt::IndexedObjectLinker::LinkDaughterToParent(e.GetRecEvent().Get(cluster), trackShine);
      } catch (const utl::EventIndexException& error) {
      }
      ++counter;
    }
    trackShine.SetNumberOfClusters(evt::rec::TrackConst::eVD, counter);
    trackShine.SetNumberOfClusters(evt::rec::TrackConst::eAll, counter + trackShine.GetNumberOfClusters(evt::rec::TrackConst::eAll));
    // trackShine.SetPotentialNumberOfClusters(evt::rec::TrackConst::eVD, counter);
    // trackShine.SetPotentialNumberOfClusters(evt::rec::TrackConst::eAll, counter+trackShine.GetPotentialNumberOfClusters(evt::rec::TrackConst::eAll));

    std::map<double, evt::Index<evt::rec::Cluster> > clusters;
    for (auto cluster = trackShine.ClustersBegin(); cluster != trackShine.ClustersEnd(); cluster++) clusters[e.GetRecEvent().Get(*cluster).GetPosition().GetZ()] = *cluster;

    evt::rec::VertexTrack& vertexTrackShine = (match && trackShine.GetNumberOfVertexTracks()) ? e.GetRecEvent().Get(*trackShine.VertexTracksBegin()) : e.GetRecEvent().Make<evt::rec::VertexTrack>();
    e.GetProcEvent().GetSilicon().GetVertexTracks().insert(vertexTrackShine.GetIndex());

    if (!match) {
      try {
        evt::IndexedObjectLinker::LinkDaughterToParent(trackShine, vertexTrackShine);
      } catch (const utl::EventIndexException& error) {
      }
    }
    if (!trackShine.GetNumberOfVertexTracks()) {
      try {
        evt::IndexedObjectLinker::LinkDaughterToParent(e.GetRecEvent().GetMainVertex(), vertexTrackShine);
      } catch (const utl::EventIndexException& error) {
      }
    }

    if (clusters.size()) {
      trackShine.SetFirstPointOnTrack(e.GetRecEvent().Get(clusters.begin()->second).GetPosition());
      trackShine.SetLastPointOnTrack(e.GetRecEvent().Get(clusters.rbegin()->second).GetPosition());

      evt::Index<evt::rec::Cluster> previous = clusters.begin()->second;
      double length = 0;
      for (auto cluster = clusters.begin(); cluster != clusters.end(); cluster++) {
        length += (e.GetRecEvent().Get(cluster->second).GetPosition() - e.GetRecEvent().Get(previous).GetPosition()).GetMag();
        previous = cluster->second;
      }
      vertexTrackShine.SetPathLength(length);
      // vertexTrackShine.SetPotentialNumberOfClusters(evt::rec::TrackConst::eVD, counter);
      // vertexTrackShine.SetPotentialNumberOfClusters(evt::rec::TrackConst::eAll, counter+vertexTrackShine.GetPotentialNumberOfClusters(evt::rec::TrackConst::eAll));

      trackShine.SetNumberOfFitClusters(evt::rec::TrackConst::eVD, counter);
      vertexTrackShine.SetNumberOfFitClusters(clusters.size());

      if (track->GetFitKf()) {
        trackShine.SetMomentum(track->GetMomentumKf());
        trackShine.SetMomentumPoint(track->GetPositionKf());
        trackShine.SetCharge(track->GetChargeKf());

        vertexTrackShine.SetMomentum(track->GetMomentumKf());
        vertexTrackShine.SetImpactPoint(track->GetPositionKf());
        vertexTrackShine.SetCharge(track->GetChargeKf());
      } else {
        trackShine.SetMomentum(utl::Vector(track->GetTLV()->Px(), track->GetTLV()->Py(), track->GetTLV()->Pz()));
        trackShine.SetMomentumPoint(trackShine.GetFirstPointOnTrack());
        trackShine.SetCharge(track->GetCharge());

        vertexTrackShine.SetMomentum(utl::Vector(track->GetTLV()->Px(), track->GetTLV()->Py(), track->GetTLV()->Pz()));
        vertexTrackShine.SetImpactPoint(trackShine.GetFirstPointOnTrack());
        vertexTrackShine.SetCharge(track->GetCharge());
      }
      e.GetProcEvent().GetSilicon().GetTrackLocalPositions()[trackShine.GetIndex()] = utl::Line(utl::Point(track->GetX_f(), track->GetY_f(), track->GetZ_f()), utl::Vector(track->GetDX_f(), track->GetDY_f(), track->GetDZ_f()));
      e.GetProcEvent().GetSilicon().GetTrackLocalCurvatures()[trackShine.GetIndex()] = track->GetCurvature();
      for(auto hit=trackShine.ClustersBegin(); hit!=trackShine.ClustersEnd(); hit++){
        if(!e.GetRecEvent().Has(*hit)) continue;
        if(e.GetRecEvent().Get(*hit).GetTPCId()!=det::TPCConst::eVD) continue;
        if(e.GetProcEvent().GetSilicon().GetTrackFirstStations().find(trackShine.GetIndex())==e.GetProcEvent().GetSilicon().GetTrackFirstStations().end() ||
          det::Detector::GetInstance().GetSilicon().GetSensor(e.GetRecEvent().Get(*hit).GetSubdetectorNumber()).GetStation() < e.GetProcEvent().GetSilicon().GetTrackFirstStations()[trackShine.GetIndex()])
//if(det::Detector::GetInstance().GetSilicon().GetSensor(e.GetRecEvent().Get(*hit).GetSubdetectorNumber()).GetStation() < e.GetProcEvent().GetSilicon().GetTrackFirstStations()[trackShine.GetIndex()])
          e.GetProcEvent().GetSilicon().GetTrackFirstStations()[trackShine.GetIndex()] = det::Detector::GetInstance().GetSilicon().GetSensor(e.GetRecEvent().Get(*hit).GetSubdetectorNumber()).GetStation();
              // std::cout<<"e.GetProcEvent().GetSilicon().GetTrackFirstStations()[trackShine.GetIndex()]="<<(int)e.GetProcEvent().GetSilicon().GetTrackFirstStations()[trackShine.GetIndex()]<<std::endl;

      }
    }

    if (na61PrimaryVertexRecoHTModule->GetPrimaryVertexStatus()) {
      evt::rec::Vertex& vertex = e.GetRecEvent().Make<evt::rec::Vertex>();
      vertex.SetPosition(utl::Point((na61PrimaryVertexRecoHTModule->GetPrimaryVertex().GetX() + Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcX()) * 0.1, (na61PrimaryVertexRecoHTModule->GetPrimaryVertex().GetY() + Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcY()) * 0.1, (na61PrimaryVertexRecoHTModule->GetPrimaryVertex().GetZ() + Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcZ()) * 0.1));
      e.GetRecEvent().SetPrimaryVertexIndex(evt::rec::VertexConst::ePrimaryVD, vertex.GetIndex());
      e.GetProcEvent().GetSilicon().GetVertexLocalPositions()[vertex.GetIndex()] = utl::Point(na61PrimaryVertexRecoHTModule->GetPrimaryVertex().GetX(), na61PrimaryVertexRecoHTModule->GetPrimaryVertex().GetY(), na61PrimaryVertexRecoHTModule->GetPrimaryVertex().GetZ());
    }
    if (na61PrimaryVertexRecoHTModuleTPC->GetPrimaryVertexStatus()) {
      evt::rec::Vertex& vertex = e.GetRecEvent().Make<evt::rec::Vertex>();
      vertex.SetPosition(utl::Point(na61PrimaryVertexRecoHTModuleTPC->GetPrimaryVertex().GetX(), na61PrimaryVertexRecoHTModuleTPC->GetPrimaryVertex().GetY(), na61PrimaryVertexRecoHTModuleTPC->GetPrimaryVertex().GetZ()));
      e.GetRecEvent().SetPrimaryVertexIndex(evt::rec::VertexConst::eUnknown, vertex.GetIndex());
    }
  }
  return eSuccess;
}

fwk::VModule::EResultFlag VdTrackFinderEfficiencySG::Finish() {
  End();
  return eSuccess;
}

void VdTrackFinderEfficiencySG::Begin(const unsigned int run) {
  End();

  //const unsigned int calibrun = Na61VdParameters::FindCalibRuns(run);
  //TFile& matchParamsFile = *TFile::Open(Form("matchParams-%06d.root", calibrun), "READ");
  //TFile & matchParamsFile = *TFile::Open(Form("/afs/cern.ch/user/a/amerzlay/public/XeLa150/VDmatchingfiles/matchParams-%0d.root",calibrun),"READ");
  TFile& matchParamsFile = *TFile::Open("matchParams-0.root", "READ");  //for simulations: perfect case
  Na61VdParameters& na61VdParameters = *Na61VdParametersManager::Instance()->GetVdParams();
  na61VdParameters.SetRunId(run);
  na61VdParameters.SetMatchParamsFile(&matchParamsFile);
  na61VdParameters.SetdOffZ(0);
  na61VdParameters.SetN(0);
  na61VdParameters.Init();
  matchParamsFile.Close();

  Na61ArmParameters& na61ArmParametersSaleve = *Na61VdParametersManager::Instance()->GetSaleveArmParams();
  na61ArmParametersSaleve.SetRunId(run);
  na61ArmParametersSaleve.Init();
  na61ArmParametersSaleve.SetVolumesInGlobal(na61VdParameters.GetSaleveArmOffset().X(), na61VdParameters.GetSaleveArmOffset().Y(), na61VdParameters.GetSaleveArmOffset().Z(), na61VdParameters.GetSaleveArmRotX(), na61VdParameters.GetSaleveArmRotY(), na61VdParameters.GetSaleveArmRotZ());

  Na61ArmParameters& na61ArmParametersJura = *Na61VdParametersManager::Instance()->GetJuraArmParams();
  na61ArmParametersJura.SetRunId(run);
  na61ArmParametersJura.Init();
  na61ArmParametersJura.SetVolumesInGlobal(na61VdParameters.GetJuraArmOffset().X(), na61VdParameters.GetJuraArmOffset().Y(), na61VdParameters.GetJuraArmOffset().Z(), na61VdParameters.GetJuraArmRotX(), na61VdParameters.GetJuraArmRotY(), na61VdParameters.GetJuraArmRotZ());

  fHistFile = new TFile("VdTrackFinderEfficiencySG.root", "recreate");
  const std::string armNames[2] = {"Saleve", "Jura"};

  for (unsigned char arm = 0; arm < 2; arm++) {
    na61VdTrackingInitModule[arm] = new Na61VdTrackingInitModule((armNames[arm] + " Vd Tracking Module").c_str(), "Vd Tracking Module");
    na61VdTrackingInitModule[arm]->SetRunId(run);
    na61VdTrackingInitModule[arm]->SetCutId(0);
    na61VdTrackingInitModule[arm]->SetNsig_dev(4.);
    na61VdTrackingInitModule[arm]->SetNsig_pvert(5.);
    na61VdTrackingInitModule[arm]->SetProductionMode(productionMode);
    na61VdTrackingInitModule[arm]->SetToJuraArm(arm);
    na61VdTrackingInitModule[arm]->SetHistDirName((std::string("VdTrackingInitModule_") + armNames[arm]).c_str());
    na61VdTrackingInitModule[arm]->Init();
    na61VdTrackingInitModule[arm]->DefineHistograms();

    //run only if ProductionMode(0)! (for efficiencies)
    if (!productionMode){
	  vdtrackingPostModule[arm] = new Na61VdTrackingPostModule((armNames[arm] + " Vd Tracking Module").c_str(),"Vd Tracking Module");
	  vdtrackingPostModule[arm]->SetRunId(run);
	  vdtrackingPostModule[arm]->SetProductionMode(productionMode);
	  //vdtrackingModule[arm]->SetBeamTrackFilter(0); // Absolate
	  vdtrackingPostModule[arm]->SetToJuraArm(arm);
	  vdtrackingPostModule[arm]->SetHistDirName((std::string("VdTrackingPostModule_") + armNames[arm]).c_str());
	  vdtrackingPostModule[arm]->Init();
      vdtrackingPostModule[arm]->DefineHistograms();

	  vdefficModule[arm] = new Na61VdEfficiencyModule((armNames[arm] + " Vd Efficiency Module").c_str(),"Vd Efficiency Module");
	  vdefficModule[arm]->SetRunId(run);
	  vdefficModule[arm]->SetToJuraArm(arm);
	  vdefficModule[arm]->SetVdAna(kTRUE);
	  vdefficModule[arm]->Init();
	  vdefficModule[arm]->DefineHistograms();
	 }

  }

  na61PrimaryVertexRecoModule = new Na61PrimaryVertexRecoModule("Primary Vertex reco", "Primary Vertex reco");
  na61PrimaryVertexRecoModule->SetRunId(run);
  na61PrimaryVertexRecoModule->Init();
  na61PrimaryVertexRecoModule->DefineHistograms();

    //run only if ProductionMode(0)! (for efficiencies)
    if (!productionMode){
	  primaryvtxpostModule = new Na61PrimaryVertexRecoPostModule("Primary Vertex Post reco","Primary Vertex Post reco");
	  primaryvtxpostModule->SetRunId(run);
	  primaryvtxpostModule->Init();
	  primaryvtxpostModule->DefineHistograms();
   }


  na61VdTrackingHTModule = new Na61VdTrackingHTModule("Vd Tracking HT Module", "Vd Tracking HT Module");
  na61VdTrackingHTModule->SetRunId(run);
  na61VdTrackingHTModule->SetProduction(1);
  na61VdTrackingHTModule->SetPosRes(0.04);
  na61VdTrackingHTModule->Setbw(0.0012, 0.0012);
  na61VdTrackingHTModule->SetMakeClusters(1);
  na61VdTrackingHTModule->SetChi2Cut(5);
  na61VdTrackingHTModule->SetVzOffset(0);
  na61VdTrackingHTModule->SetHistDirName("VdTrackingHTModule");
  na61VdTrackingHTModule->SetOutTableName("Vd HT Tracks");
  na61VdTrackingHTModule->SetAdd4HitTracks(true);
  na61VdTrackingHTModule->Init();
  na61VdTrackingHTModule->DefineHistograms();

  na61PrimaryVertexRecoHTModule = new Na61PrimaryVertexRecoHTModule("Primary Vertex reco", "Primary Vertex reco");
  na61PrimaryVertexRecoHTModule->SetRunId(run);
  na61PrimaryVertexRecoHTModule->SetTrackInput(0);
  na61PrimaryVertexRecoHTModule->SetZprim(47.15);
  na61PrimaryVertexRecoHTModule->SetZCut(2.5);
  na61PrimaryVertexRecoHTModule->Init();
  na61PrimaryVertexRecoHTModule->DefineHistograms();

  na61VdTrackingHTPackage = new Na61VdTrackingHTPackage("Vd Tracking HT Package", "Vd Tracking HT Package");
  na61VdTrackingHTPackage->SetRunId(run);
  na61VdTrackingHTPackage->SetProduction(1);
  na61VdTrackingHTPackage->SetPosRes(0.04);
  na61VdTrackingHTPackage->Setbw(0.0012, 0.0012);
  na61VdTrackingHTPackage->SetMakeClusters(1);
  na61VdTrackingHTPackage->SetChi2Cut(5);
  na61VdTrackingHTPackage->SetVzOffset(0);
  na61VdTrackingHTPackage->Init();
  na61VdTrackingHTPackage->DefineHistograms();

  na61VdTpcMatchingModule = new Na61VdTpcMatchingModule("VD-TPC Matching", "VD-TPC Matching");
  na61VdTpcMatchingModule->SetRunId(run);
  na61VdTpcMatchingModule->Init();
  na61VdTpcMatchingModule->DefineHistograms();

  na61PrimaryVertexRecoHTModuleTPC = new Na61PrimaryVertexRecoHTModule("Primary Vertex Tpc", "Primary Vertex Tpc");
  na61PrimaryVertexRecoHTModuleTPC->SetRunId(run);
  na61PrimaryVertexRecoHTModuleTPC->SetTrackInput(1);
  na61PrimaryVertexRecoHTModuleTPC->SetZprim(602.8);
  na61PrimaryVertexRecoHTModuleTPC->SetZCut(2.0);
  na61PrimaryVertexRecoHTModuleTPC->Init();
  na61PrimaryVertexRecoHTModuleTPC->DefineHistograms();



}

void VdTrackFinderEfficiencySG::End() {
  for (unsigned char arm = 0; arm < 2; arm++) {
    if (na61VdTrackingInitModule[arm]) {
      na61VdTrackingInitModule[arm]->Finish();
      delete na61VdTrackingInitModule[arm];
      na61VdTrackingInitModule[arm] = NULL;
    }
  }
  if (na61PrimaryVertexRecoModule) {
    na61PrimaryVertexRecoModule->Finish();
    delete na61PrimaryVertexRecoModule;
    na61PrimaryVertexRecoModule = NULL;
  }
  if (na61VdTrackingHTModule) {
    na61VdTrackingHTModule->Finish();
    delete na61VdTrackingHTModule;
    na61VdTrackingHTModule = NULL;
  }
  if (na61PrimaryVertexRecoHTModule) {
    na61PrimaryVertexRecoHTModule->Finish();
    delete na61PrimaryVertexRecoHTModule;
    na61PrimaryVertexRecoHTModule = NULL;
  }
  if (na61VdTrackingHTPackage) {
    na61VdTrackingHTPackage->Finish();
    delete na61VdTrackingHTPackage;
    na61VdTrackingHTPackage = NULL;
  }
  if (na61VdTpcMatchingModule) {
    na61VdTpcMatchingModule->Finish();
    delete na61VdTpcMatchingModule;
    na61VdTpcMatchingModule = NULL;
  }
  if (na61PrimaryVertexRecoHTModuleTPC) {
    na61PrimaryVertexRecoHTModuleTPC->Finish();
    delete na61PrimaryVertexRecoHTModuleTPC;
    na61PrimaryVertexRecoHTModuleTPC = NULL;
  }

  //run only if ProductionMode(0)! (for efficiencies)
  if (!productionMode){
	  for (unsigned char arm = 0; arm < 2; arm++) {
		if (vdtrackingPostModule[arm]) {
		  vdtrackingPostModule[arm]->Finish();
		  delete vdtrackingPostModule[arm];
		  vdtrackingPostModule[arm] = NULL;
		}
		if (vdefficModule[arm]) {
		  vdefficModule[arm]->Finish();
		  delete vdefficModule[arm];
		  vdefficModule[arm] = NULL;
		}
	  }
	  if (primaryvtxpostModule) {
		primaryvtxpostModule->Finish();
		delete primaryvtxpostModule;
		primaryvtxpostModule = NULL;
	  }
  }


  if (fHistFile) {
    fHistFile->Write();
    fHistFile->Close();
    delete fHistFile;
    fHistFile = NULL;
  }
}
}
