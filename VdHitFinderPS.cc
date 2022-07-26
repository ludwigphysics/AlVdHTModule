#include "VdHitFinderPS.h"
#include "Na61FrameMergingModule.h"
#include "Na61HitProducerModule.h"
#include "UDataTable.h"
#include "USensorPixel.h"
#include <fwk/CentralConfig.h>
#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>
#include <evt/raw/Silicon.h>
#include <evt/BOSBank.h>
#include <evt/BOSRecord.h>
#include <det/Detector.h>
#include <evt/Event.h>
#include <evt/raw/Silicon.h>
#include <evt/proc/Silicon.h>
#include <det/Silicon.h>
#include <fwk/CentralConfig.h>
#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>
#include <sstream>

namespace VdHitFinderPS {

fwk::VModule::EResultFlag VdHitFinderPS::Init() {

  utl::Branch vdReconstructoinSGXML =
    fwk::CentralConfig::GetInstance().GetTopBranch("VdReconstructionSG");
  InitVerbosity(vdReconstructoinSGXML);

  vdReconstructoinSGXML.GetChild("matchParamsPath").GetData(fMatchParamsPath);
  char fullpath[PATH_MAX];
  if (not realpath(fMatchParamsPath.c_str(), fullpath)) {
    ERROR(Form("Invalid path to match-params (%s)", fMatchParamsPath.c_str()));
    throw utl::XMLValidationException {"invalid value of an entry"};
  }
  fMatchParamsPath = fullpath;

  utl::Branch vdHitFinderPSXML =
    vdReconstructoinSGXML.GetChild("VdHitFinderPS");
  vdHitFinderPSXML.GetChild("histogramsPath").GetData(fHistosPath);
  if (fHistosPath.empty()) {
    ERROR("Please fill in the <histogramsPath>...</histogramsPath> element inside <VdHitFinderPS>...</VdHitFinderPS>");
    throw utl::XMLValidationException {"mandatory element not specified"};
  }

  prevRun = 0;
  return eSuccess;
}

fwk::VModule::EResultFlag VdHitFinderPS::Process(evt::Event &e, const utl::AttributeMap &) {
  const evt::EventHeader &eventHeader = e.GetEventHeader();
  if (!eventHeader.IsInitialized()) {
    WARNING("No EventHeader is present in Event. This is a serious problem!!! PLEASE DO CHECK WHAT IS GOING ON IF THIS EFFECT IS NOT EXPECTED!!!");
    return eFailure;
  }
  det::Detector::GetInstance().Update(eventHeader.GetTime(), eventHeader.GetRunNumber());

  if (eventHeader.GetRunNumber() != prevRun) {
    Begin(eventHeader.GetRunNumber(), e.IsSimulation());
    prevRun = eventHeader.GetRunNumber();
  }

  const std::string nameStandalone[8] = {"Hits Vds1_0", "Hits Vds2_0", "Hits Vds3_0", "Hits Vds3_1", "Hits Vds4_0", "Hits Vds4_1", "Hits Vds4_2", "Hits Vds4_3"};
  const std::string nameShine[2][8] = {{"vd1_s", "vd2_s", "vd3_s_d", "vd3_s_u", "vd4_s1_d", "vd4_s1_u", "vd4_s2_d", "vd4_s2_u"}, {"vd1_j", "vd2_j", "vd3_j_d", "vd3_j_u", "vd4_j1_d", "vd4_j1_u", "vd4_j2_d", "vd4_j2_u"}};
  const std::map<std::string, unsigned char> shineName[2] = {{{"vd1_s", 0}, {"vd2_s", 1}, {"vd3_s_d", 2}, {"vd3_s_u", 3}, {"vd4_s1_d", 4}, {"vd4_s1_u", 5}, {"vd4_s2_d", 6}, {"vd4_s2_u", 7}}, {{"vd1_j", 0}, {"vd2_j", 1}, {"vd3_j_d", 2}, {"vd3_j_u", 3}, {"vd4_j1_d", 4}, {"vd4_j1_u", 5}, {"vd4_j2_d", 6}, {"vd4_j2_u", 7}}};
	std::map<unsigned char, unsigned char> shineStandaloneName = {{0, 0}, {2, 1}, {4, 2}, {5, 3}, {10, 4}, {11, 5}, {8, 6}, {9, 7}, 
		                                                          {1, 0}, {3, 1}, {6, 2}, {7, 3}, {12, 4}, {13, 5}, {14, 6}, {15, 7}};

  // Loop over Saleve and Jura arms
  for (unsigned char arm = 0; arm < 2; arm++) {
    // Create event
    UVdEvent event;
    event.SetOwner();

    UDataTable *clusters[8][5];
    for (unsigned char s = 0; s < 8; s++) {
     for (unsigned char f = 0; f < 5; f++) {
        clusters[s][f] = new UDataTable(Form("Pixels %s f%d", na61FrameMergingModule[arm]->fSensorNames[s].Data(), f));
        event.AddDataTable(clusters[s][f]);
      }
    }

    double timers[5];

    evt::raw::Silicon &rawSilicon = e.GetRawEvent().GetSilicon();
    const det::Silicon &detSilicon = det::Detector::GetInstance().GetSilicon();
    for (auto sensorIter = rawSilicon.GetSensors().begin();
        sensorIter != rawSilicon.GetSensors().end(); sensorIter++)
    {
      const det::SiliconSensor &detSensor =
        detSilicon.GetSensor(sensorIter->GetSubdetectorNumber());

      if ((detSensor.GetLadder() < 0) == arm)
        continue;

      for (auto frameIter = sensorIter->GetFrames().begin();
          frameIter != sensorIter->GetFrames().end() &&
          std::distance(sensorIter->GetFrames().begin(), frameIter) < 5;
          frameIter++)
      {
        for (auto pixelIter = frameIter->GetPixels().begin();
            pixelIter != frameIter->GetPixels().end(); pixelIter++)
        {
          if (pixelIter->GetLine() < 576 || pixelIter->GetColumn() < 1152) {

            clusters
              [shineName[arm].find(detSensor.GetComponentId())->second]
              [std::distance(sensorIter->GetFrames().begin(), frameIter)]
              ->Add(new USensorPixel(pixelIter->GetLine(), pixelIter->GetColumn()));
          }
        }
        timers[std::distance(sensorIter->GetFrames().begin(), frameIter)]
          = frameIter->GetTime() * 0.001;
        //std::cout<<"arm="<<(unsigned char)arm<<" "<<sensorIter->GetSiliconDetector().GetComponentId()<<" "<<frameIter->GetPixels().size()<<std::endl;
      }
    }

    na61FrameMergingModule[arm]->SetFrameTimers(0, timers[0], timers[1], timers[2], timers[3], timers[4]);
    na61FrameMergingModule[arm]->Event(&event, &event);
    na61HitProducerModule[arm]->Event(&event, &event);
    for (unsigned char s = 0; s < 8; s++) {
      UDataTable *hitsOnSensor = event.GetDataTable(nameStandalone[s].c_str());
      if (hitsOnSensor) {
        for (int i = 0; i < hitsOnSensor->GetEntries(); i++) {
          USensorHit *hit = (USensorHit *)hitsOnSensor->At(i);
          if (!hit)
            continue;

          const det::Silicon &detSilicon = det::Detector::GetInstance().GetSilicon();
          const det::SiliconSensor &detSensor =
            detSilicon.GetSensor(detSilicon.SensorIdToSubdetectorNumber(nameShine[arm][s]));

          evt::rec::Cluster &cluster = e.GetRecEvent().Make<evt::rec::Cluster>();
          cluster.SetPositionUncertainty(evt::rec::ClusterConst::eX, 5.0 * utl::micrometer);
          cluster.SetPositionUncertainty(evt::rec::ClusterConst::eY, 5.0 * utl::micrometer);
          cluster.SetNumberOfPixels(hit->GetClusterSize());
          cluster.SetTPCId(detSensor.GetDetector());
          cluster.SetSubdetectorNumber(detSensor.GetSubdetectorNumber());
          // cluster.SetFirstCoordinate(line);
          // cluster.SetSecondCoordinate(column);
          evt::proc::SiliconSensor &procSensor =
            e.GetProcEvent().GetSilicon().HasSensor(detSensor.GetSubdetectorNumber())
            ? e.GetProcEvent().GetSilicon().GetSensor(detSensor.GetSubdetectorNumber())
            : e.GetProcEvent().GetSilicon().MakeSensor(detSensor.GetSubdetectorNumber());
          e.GetProcEvent().GetSilicon().GetClusterLocalPositions()[cluster.GetIndex()] =
            utl::Point(0.1 * hit->GetX(), 0.1 * hit->GetY(), 0.1 * hit->GetZ());
          //utl::Point clusterShinePos;
          // option 1: use xml file
          // double newX, newY;
          // if(arm) {
          //  newX = -.1 * hit->GetLocalX();
          //  newY = -.1 * hit->GetLocalY();
          //}
          // else {
          //  newX = .1 * hit->GetLocalX();
          //  newY = .1 * hit->GetLocalY();
          //}
          // clusterShinePos = utl::Point(newX, newY, 0, procSensor.GetSiliconDetector().GetComponentCoordinateSystem());
          // option 2: use VD-standalone global geometry
          //LocalToShineGlobal(arm, s, *hit, clusterShinePos);


          unsigned char shineSensorName = shineStandaloneName[detSensor.GetSubdetectorNumber()];
          //cout<<(unsigned int)detSensor.GetSubdetectorNumber()<<" -> sensor="<<(unsigned int)s<<" arm="<<(unsigned int)arm<<endl;
          utl::Point localPos;
          if (arm) localPos.Set(hit->GetX()*10.,-hit->GetY()*10.,0.0);
          else localPos.Set(-hit->GetX()*10.,hit->GetY()*10.,0.0);
          utl::Point vdstandalonePos;
          LocalToGlobal(arm, shineSensorName, localPos, vdstandalonePos); 

          utl::Point clusterShinePos;
          LocalToShineGlobal(arm, vdstandalonePos, clusterShinePos);


          cluster.SetPosition(clusterShinePos);
          e.GetProcEvent().GetSilicon().GetClusters().insert(cluster.GetIndex());
          procSensor.GetClusters().insert(cluster.GetIndex());
        }
      }
    }
  }
  return eSuccess;
}

fwk::VModule::EResultFlag VdHitFinderPS::Finish() {
  End();
  return eSuccess;
}

void VdHitFinderPS::Begin(const unsigned int run, bool isSim) {
  End();

  const unsigned int calibrun = Na61VdParameters::FindCalibRuns(run);
  // XXX: delete?
  TFile* matchParamsFile;
  if (isSim)
    matchParamsFile = TFile::Open(Form("%s/matchParams-0.root", fMatchParamsPath.c_str()), "READ");
  else
    matchParamsFile = TFile::Open(Form("%s/matchParams-%0d.root", fMatchParamsPath.c_str(), calibrun), "READ");

  Na61VdParameters &na61VdParameters = *Na61VdParametersManager::Instance()->GetVdParams();
  na61VdParameters.SetRunId(run);
  na61VdParameters.SetMatchParamsFile(matchParamsFile);
  na61VdParameters.SetdOffZ(0);
  na61VdParameters.SetN(0);
  na61VdParameters.Init();
  matchParamsFile->Close();

  Na61ArmParameters &na61ArmParametersSaleve = *Na61VdParametersManager::Instance()->GetSaleveArmParams();
  na61ArmParametersSaleve.SetRunId(run);
  na61ArmParametersSaleve.Init();
  na61ArmParametersSaleve.SetVolumesInGlobal(na61VdParameters.GetSaleveArmOffset().X(), na61VdParameters.GetSaleveArmOffset().Y(), na61VdParameters.GetSaleveArmOffset().Z(), na61VdParameters.GetSaleveArmRotX(), na61VdParameters.GetSaleveArmRotY(), na61VdParameters.GetSaleveArmRotZ());

  Na61ArmParameters &na61ArmParametersJura = *Na61VdParametersManager::Instance()->GetJuraArmParams();
  na61ArmParametersJura.SetRunId(run);
  na61ArmParametersJura.Init();
  na61ArmParametersJura.SetVolumesInGlobal(na61VdParameters.GetJuraArmOffset().X(), na61VdParameters.GetJuraArmOffset().Y(), na61VdParameters.GetJuraArmOffset().Z(), na61VdParameters.GetJuraArmRotX(), na61VdParameters.GetJuraArmRotY(), na61VdParameters.GetJuraArmRotZ());

  fHistFile = new TFile(fHistosPath.c_str(), "recreate");
  const std::string armNames[2] = {"Saleve", "Jura"};

  // Loop over each arm and create standalone event processing module objects
  for (unsigned char arm = 0; arm < 2; arm++) {
    na61FrameMergingModule[arm] = new Na61FrameMergingModule((armNames[arm] + " Frame Merging Module").c_str(), "Frame Merging Module");
    na61FrameMergingModule[arm]->SetCutOnNoisyPixels(1);
    na61FrameMergingModule[arm]->SetRunId(run);
    na61FrameMergingModule[arm]->SetToJuraArm(arm);
    na61FrameMergingModule[arm]->SetNewMerge(isSim ? 2 : 1); //2=simulations
    na61FrameMergingModule[arm]->SetRowOverlap(15);
    na61FrameMergingModule[arm]->SetHomePath(".");
    na61FrameMergingModule[arm]->Init();
    na61FrameMergingModule[arm]->DefineHistograms();

    na61HitProducerModule[arm] = new Na61HitProducerModule((armNames[arm] + " Hit Producer Module").c_str(), "Hit Producer Module");
    na61HitProducerModule[arm]->SetProductionMode(1);
    na61HitProducerModule[arm]->SetToJuraArm(arm);
    na61HitProducerModule[arm]->SetRunId(run);
    na61HitProducerModule[arm]->SetDrotX(0, 0);
    na61HitProducerModule[arm]->SetDrotY(0, 0);
    na61HitProducerModule[arm]->SetDrotZ(0, 0);
    na61HitProducerModule[arm]->SetVolumeDz(0, 0);
    na61HitProducerModule[arm]->Init();
    na61HitProducerModule[arm]->DefineHistograms();
  }
}

void VdHitFinderPS::End() {
  // Loop over each arm and clean up standalone event processing modules
  for (unsigned char arm = 0; arm < 2; arm++) {
    if (na61FrameMergingModule[arm]) {
      na61FrameMergingModule[arm]->Finish();
      delete na61FrameMergingModule[arm];
      na61FrameMergingModule[arm] = NULL;
    }
    if (na61HitProducerModule[arm]) {
      na61HitProducerModule[arm]->Finish();
      delete na61HitProducerModule[arm];
      na61HitProducerModule[arm] = NULL;
    }
  }
  if (fHistFile) {
    fHistFile->Write();
    fHistFile->Close();
    delete fHistFile;
    fHistFile = NULL;
  }
}

//_____________________________________________________________________
void VdHitFinderPS::LocalToGlobal(unsigned char arm, unsigned char sensor, utl::Point &localPos,utl::Point &vdstandalonePos) {
  Na61ArmParameters &armParams = arm ? *Na61VdParametersManager::Instance()->GetJuraArmParams() : *Na61VdParametersManager::Instance()->GetSaleveArmParams();

 //initial arm params in vd-standalone coordinate system
  double alpha = armParams.GetRotZ(sensor);
  double beta = armParams.GetRotY(sensor);
  double gamma = armParams.GetRotX(sensor);
  const double VolumeX = armParams.GetVolumeX(sensor);
  const double VolumeY = armParams.GetVolumeY(sensor);
  double VolumeZ = armParams.GetVolumeZ(sensor);

//cout<<"arm="<<(unsigned int)arm<<" sensor="<<(unsigned int)sensor<<" "<<VolumeX<<" "<<VolumeY<<" "<<VolumeZ<<" "<<alpha<< " "<<beta<<" "<<gamma<<endl;

  double sa = TMath::Sin(alpha);
  double ca = TMath::Cos(alpha);
  double sb = TMath::Sin(beta);
  double cb = TMath::Cos(beta);
  double sg = TMath::Sin(gamma);
  double cg = TMath::Cos(gamma);

    double x1 = localPos.GetX();
    double y1 = localPos.GetY();
    double z1 = localPos.GetZ();

    double x2;
    double y2;
    double z2;

    // ZXY
    x2 =  (ca * cb - sa * sg * sb) * x1 + sa * cg * y1 + ( ca * sb + sa * sg * cb) * z1;
    y2 = -(sa * cb + ca * sg * sb) * x1 + ca * cg * y1 + (-sa * sb + ca * sg * cb) * z1;
    z2 =                  -cg * sb * x1      - sg * y1 +                   cg * cb * z1;

    vdstandalonePos.Set(x2 + VolumeX, y2 + VolumeY, z2 + VolumeZ);

    // if(ii==0)cout<<"VolumeX="<<VolumeX<<" VolumeY="<<VolumeY<<" VolumeZ="<<VolumeZ<<"  x1="<<x1<<"  y1="<<y1<<"  z2="<<z2<<" z2+V ="<<hit->GetZ()<<endl;
	
  
}


void VdHitFinderPS::LocalToShineGlobal(unsigned char arm, utl::Point &vdstandalonePos, utl::Point &shinePos) {
  Na61VdParameters &vdParams = *(Na61VdParametersManager::Instance())->GetVdParams();

  double rotx, roty, rotz;
  double dx, dy, dz;
  if (arm) {
    rotx = vdParams.GetJuraArmRotX();
    roty = vdParams.GetJuraArmRotY();
    rotz = vdParams.GetJuraArmRotZ();
    dx = vdParams.GetJuraArmOffset().X();
    dy = vdParams.GetJuraArmOffset().Y();
    dz = vdParams.GetJuraArmOffset().Z();
  } else {
    rotx = vdParams.GetSaleveArmRotX();
    roty = vdParams.GetSaleveArmRotY();
    rotz = vdParams.GetSaleveArmRotZ();
    dx = vdParams.GetSaleveArmOffset().X();
    dy = vdParams.GetSaleveArmOffset().Y();
    dz = vdParams.GetSaleveArmOffset().Z();
  }
  double sa = TMath::Sin(rotz);
  double ca = TMath::Cos(rotz);
  double sb = TMath::Sin(roty);
  double cb = TMath::Cos(roty);
  double sg = TMath::Sin(rotx);
  double cg = TMath::Cos(rotx);
    
  const double x0 = vdstandalonePos.GetX();
  const double y0 = vdstandalonePos.GetY();
  double z0 = vdstandalonePos.GetZ();
  const double VolumeZ0 = 75.0;
  z0 = z0 - VolumeZ0;
  
  double x1 = (ca * cb - sa * sg * sb) * x0 + sa * cg * y0 + (ca * sb + sa * sg * cb) * z0;
  double y1 = -(sa * cb + ca * sg * sb) * x0 + ca * cg * y0 + (-sa * sb + ca * sg * cb) * z0;
  double z1 = -cg * sb * x0 - sg * y0 + cg * cb * z0;
  z1 = z1 + VolumeZ0;
  
  x1 = x1 + dx;
  y1 = y1 + dy;
  z1 = z1 + dz;
  
  // Add shift relative to TPC
  const double fOffsetX = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcX();
  const double fOffsetY = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcY();
  const double fOffsetZ = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcZ();

    shinePos.Set(0.1 * (x1 + fOffsetX), 0.1 * (y1 + fOffsetY), 0.1 * (z1 + fOffsetZ));
 

  return;
}


}
