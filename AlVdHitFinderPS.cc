#include "AlVdHitFinderPS.h"
//#include "Na61FrameMergingModule.h"
#include "Na61AlHitProducerModule.h"
#include "Na61HotPixelsModule.h"
#include "UDataTable.h"
#include "UChipPixel.h"
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

namespace AlVdHitFinderPS {

fwk::VModule::EResultFlag AlVdHitFinderPS::Init()
{
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

  fNoiseRun = 0;

  utl::Branch vdHitFinderPSXML =
    vdReconstructoinSGXML.GetChild("VdHitFinderPS");
  vdHitFinderPSXML.GetChild("histogramsPath").GetData(fHistosPath);
  if (fHistosPath.empty()) {
    ERROR("Please fill in the <histogramsPath>...</histogramsPath> element inside <AlVdHitFinderPS>...</AlVdHitFinderPS>");
    throw utl::XMLValidationException {"mandatory element not specified"};
  }

  prevRun = 0;
  return eSuccess;
}

fwk::VModule::EResultFlag AlVdHitFinderPS::Process(evt::Event &e, const utl::AttributeMap &) {
  const evt::EventHeader &eventHeader = e.GetEventHeader();
  if (!eventHeader.IsInitialized()) {
    WARNING("No EventHeader is present in Event. This is a serious problem!!! PLEASE DO CHECK WHAT IS GOING ON IF THIS EFFECT IS NOT EXPECTED!!!");
    //return eFailure; // as instructed by Ivan(PS)
    return eContinueLoop; // as instructed by Ivan(PS)
  }
  det::Detector::GetInstance().Update(eventHeader.GetTime(), eventHeader.GetRunNumber());

  if (eventHeader.GetRunNumber() != prevRun) {
    Begin(eventHeader.GetRunNumber(), e.IsSimulation());
    prevRun = eventHeader.GetRunNumber();
  }

  const std::string nameStandalone[34] = {
    "Hits Al1_0", "Hits Al1_1", "Hits Al1_2",
    "Hits Al2_0", "Hits Al2_1", "Hits Al2_2",  "Hits Al2_3", "Hits Al2_4", "Hits Al2_5",
    "Hits Al3_0", "Hits Al3_1", "Hits Al3_2", "Hits Al3_3", "Hits Al3_4",
    "Hits Al3_5", "Hits Al3_6", "Hits Al3_7", "Hits Al3_8", "Hits Al3_9",
    "Hits Al4_0", "Hits Al4_1", "Hits Al4_2", "Hits Al4_3", "Hits Al4_4",
    "Hits Al4_5", "Hits Al4_6", "Hits Al4_7", "Hits Al4_8", "Hits Al4_9",
    "Hits Al4_10","Hits Al4_11","Hits Al4_12","Hits Al4_13","Hits Al4_14"
  };
  const std::string nameShine[2][34] = {
    {"vds1_0", "vds1_1", "vds1_2", 
      "vds2_0", "vds2_1", "vds2_2", "vds2_3", "vds2_4", "vds2_5",
      "vds3_0", "vds3_1", "vds3_2", "vds3_3", "vds3_4", "vds3_5", "vds3_6", "vds3_7", "vds3_8", "vds3_9",
      "vds4_0", "vds4_1", "vds4_2", "vds4_3", "vds4_4", "vds4_5", "vds4_6", "vds4_7", "vds4_8", "vds4_9",
      "vds4_10","vds4_11","vds4_12","vds4_13","vds4_14"
    },
    {"vdj1_0", "vdj1_1", "vdj1_2", 
      "vdj2_0", "vdj2_1", "vdj2_2", "vdj2_3", "vdj2_4", "vdj2_5",
      "vdj3_0", "vdj3_1", "vdj3_2", "vdj3_3", "vdj3_4", "vdj3_5", "vdj3_6", "vdj3_7", "vdj3_8", "vdj3_9",
      "vdj4_0", "vdj4_1", "vdj4_2", "vdj4_3", "vdj4_4", "vdj4_5", "vdj4_6", "vdj4_7", "vdj4_8", "vdj4_9",
      "vdj4_10","vdj4_11","vdj4_12","vdj4_13","vdj4_14"
    }
  };

  const std::map<std::string, unsigned char> shineName[2] = {
    {{"vds1_0", 0},  {"vds1_1", 1},  {"vds1_2", 2}, 
      {"vds2_0", 3},  {"vds2_1", 4},  {"vds2_2", 5},  {"vds2_3", 6},  {"vds2_4", 7},  {"vds2_5", 8},
      {"vds3_0", 9},  {"vds3_1",10},  {"vds3_2",11},  {"vds3_3",12},  {"vds3_4",13}, 
      {"vds3_5",14},  {"vds3_6",15},  {"vds3_7",16},  {"vds3_8",17},  {"vds3_9",18}, 
      {"vds4_0",19},  {"vds4_1",20},  {"vds4_2",21},  {"vds4_3",22},  {"vds4_4",23}, 
      {"vds4_5",24},  {"vds4_6",25},  {"vds4_7",26},  {"vds4_8",27},  {"vds4_9",28}, 
      {"vds4_10",29}, {"vds4_11",30}, {"vds4_12",31}, {"vds4_13",32}, {"vds4_14",33}, 
    },
    {{"vdj1_0", 0},  {"vdj1_1", 1},  {"vdj1_2", 2}, 
      {"vdj2_0", 3},  {"vdj2_1", 4},  {"vdj2_2", 5},  {"vdj2_3", 6},  {"vdj2_4", 7},  {"vdj2_5", 8},
      {"vdj3_0", 9},  {"vdj3_1",10},  {"vdj3_2",11},  {"vdj3_3",12},  {"vdj3_4",13}, 
      {"vdj3_5",14},  {"vdj3_6",15},  {"vdj3_7",16},  {"vdj3_8",17},  {"vdj3_9",18}, 
      {"vdj4_0",19},  {"vdj4_1",20},  {"vdj4_2",21},  {"vdj4_3",22},  {"vdj4_4",23}, 
      {"vdj4_5",24},  {"vdj4_6",25},  {"vdj4_7",26},  {"vdj4_8",27},  {"vdj4_9",28}, 
      {"vdj4_10",29}, {"vdj4_11",30}, {"vdj4_12",31}, {"vdj4_13",32}, {"vdj4_14",33} 
    }
  };

  std::map<unsigned char, unsigned char> shineStandaloneName = 
  {
    {0, 0}, {1, 1}, {2, 2}, {6, 3}, {7, 4}, {8, 5}, {9, 6}, {10, 7}, {11, 8}, {18, 9}, {19,10}, {20,11}, {21,12}, {22,13}, {23,14}, 
    {24,15}, {25,16}, {26,17}, {27,18}, {38,19}, {39,20}, {40,21}, {41,22}, {42,23}, {43,24}, {44,25}, {45,26}, {46,27}, {47,28}, {48,29}, 
    {49,30}, {50,31}, {51,32}, {52,33},  

    {3, 0}, {4, 1}, {5, 2}, {12, 3}, {13, 4}, {14, 5}, {15, 6}, {16, 7}, {17, 8}, {28, 9}, {29,10}, {30,11}, {31,12}, {32,13}, {33,14},
    {34,15}, {35,16}, {36,17}, {37,18}, {53,19}, {54,20}, {55,21}, {56,22}, {57,23}, {58,24}, {59,25}, {60,26}, {61,27}, {62,28}, {63,29}, 
    {64,30}, {65,31}, {66,32}, {67,33}      
  };


  // Loop over Saleve and Jura arms
  for (unsigned char arm = 0; arm < 2; arm++) {
    // Create event
    UVdEvent event;
    event.SetOwner();

    UDataTable *clusters[34];
    for (unsigned char s = 0; s < 34; s++) {

      if(arm==0) clusters[s] = new UDataTable(Form("Saleve Pixels %s merged", na61HitProducerModule[arm]->fAlSensorNames[s].Data()));
      else 	 clusters[s] = new UDataTable(Form("Jura Pixels %s merged", na61HitProducerModule[arm]->fAlSensorNames[s].Data()));

      event.AddDataTable(clusters[s]);

    }

    // data table for pixels defined but are still empty

    for (std::list<evt::raw::SiliconSensor>::iterator sensorIter = e.GetRawEvent().GetSilicon().GetSensors().begin(); sensorIter != e.GetRawEvent().GetSilicon().GetSensors().end(); sensorIter++) {

      const det::Silicon &detSilicon = det::Detector::GetInstance().GetSilicon();
      const det::SiliconSensor &detSensor =
        detSilicon.GetSensor(sensorIter->GetSubdetectorNumber());

      if ((detSensor.GetLadder() < 0) == arm)
        continue;

      // we keep iterations over frames but ALPIDE we have only one frame
      for (auto frameIter = sensorIter->GetFrames().begin();
          frameIter != sensorIter->GetFrames().end() &&
          std::distance(sensorIter->GetFrames().begin(), frameIter) < 1;
          frameIter++)
      {
        for (auto pixelIter = frameIter->GetPixels().begin();
            pixelIter != frameIter->GetPixels().end(); pixelIter++)
        {
          if (pixelIter->GetLine() < 512 || pixelIter->GetColumn() < 1024){
            clusters[shineName[arm].at(detSensor.GetComponentId())]
              ->Add(new UChipPixel(pixelIter->GetLine()+1, pixelIter->GetColumn()+1));
          }
        }
      }
    }


    hotpixelsModule[arm]->Event(&event, &event);

    if(fNoiseRun){
      hotpixelsModule[arm]->Event(&event, &event);
      continue;
    }

    na61HitProducerModule[arm]->Event(&event, &event);

    for (unsigned char s = 0; s < 34; s++) {
      UDataTable *hitsOnSensor = event.GetDataTable(nameStandalone[s].c_str());
      if (hitsOnSensor) {
        for (int i = 0; i < hitsOnSensor->GetEntries(); i++) {
          USensorHit *hit = (USensorHit *)hitsOnSensor->At(i);
          if (!hit) continue;

          const det::Silicon &detSilicon = det::Detector::GetInstance().GetSilicon();
          const det::SiliconSensor &detSensor =
            detSilicon.GetSensor(detSilicon.SensorIdToSubdetectorNumber(nameShine[arm][s]));

          evt::rec::Cluster &cluster = e.GetRecEvent().Make<evt::rec::Cluster>();
          cluster.SetPositionUncertainty(evt::rec::ClusterConst::eX, 5.0 * utl::micrometer);
          cluster.SetPositionUncertainty(evt::rec::ClusterConst::eY, 5.0 * utl::micrometer);
          cluster.SetNumberOfPixels(hit->GetClusterSize());
          cluster.SetTPCId(detSensor.GetDetector());
          cluster.SetSubdetectorNumber(detSensor.GetSubdetectorNumber());
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


          //unsigned char shineSensorName = shineStandaloneName[detSensor.GetSubdetectorNumber()];
          //cout<<(unsigned int)detSensor.GetSubdetectorNumber()<<" -> sensor="<<(unsigned int)s<<" arm="<<(unsigned int)arm<<endl;
          utl::Point localPos;

          /* /// ORIGINAL SHINE BLOCK ////////////////////////////////
             if (arm) localPos.Set(hit->GetX()*10.,-hit->GetY()*10.,0.0);
             else localPos.Set(-hit->GetX()*10.,hit->GetY()*10.,0.0);

             utl::Point vdstandalonePos;
             LocalToGlobal(arm, shineSensorName, localPos, vdstandalonePos); 

             utl::Point clusterShinePos;
             LocalToShineGlobal(arm, shineSensorName, vdstandalonePos, clusterShinePos);      

             cluster.SetPosition(clusterShinePos); //normal transfer
             */ ////////////////////////////////////////////////////////

          /////////////////  BLOCK for VD HITs TRANSFER
          //cout<<"x: "<<hit->GetX()<<" "<<hit->GetLocalX()<<"  sensorName="<<hit->GetSensorName()<<endl;
          //cout<<"y: "<<hit->GetY()<<" "<<hit->GetLocalY()<<endl;
          localPos.Set(hit->GetX()*0.1,hit->GetY()*0.1,hit->GetZ()*0.1);	
          cluster.SetPosition(localPos); // to be used only to transfer to vd-stanalone
          // trick to transfer local positions in sensor (not elegant but should work for us)
          cluster.SetEnergyDeposit(hit->GetLocalX()*0.1); 
          cluster.SetMaxADC(hit->GetLocalY()*0.1);  
          ////////////////////////////////////////

          e.GetProcEvent().GetSilicon().GetClusters().insert(cluster.GetIndex());
          procSensor.GetClusters().insert(cluster.GetIndex());
        }
      }
    }
  }

  return eSuccess;
}

//_______________________________________________________________________________________
fwk::VModule::EResultFlag AlVdHitFinderPS::Finish() {
  End();
  return eSuccess;
}

//_______________________________________________________________________________________
void AlVdHitFinderPS::Begin(const unsigned int run, bool isSim) {
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

  Na61AlVdArmParameters &na61ArmParametersSaleve = *Na61VdParametersManager::Instance()->GetSaleveAlVdArmParams();
  na61ArmParametersSaleve.SetRunId(run);
  na61ArmParametersSaleve.Init();
  na61ArmParametersSaleve.SetVolumesInGlobal(na61VdParameters.GetSaleveArmOffset().X(), na61VdParameters.GetSaleveArmOffset().Y(), na61VdParameters.GetSaleveArmOffset().Z(), na61VdParameters.GetSaleveArmRotX(), na61VdParameters.GetSaleveArmRotY(), na61VdParameters.GetSaleveArmRotZ());

  Na61AlVdArmParameters &na61ArmParametersJura = *Na61VdParametersManager::Instance()->GetJuraAlVdArmParams();
  na61ArmParametersJura.SetRunId(run);
  na61ArmParametersJura.Init();
  na61ArmParametersJura.SetVolumesInGlobal(na61VdParameters.GetJuraArmOffset().X(), na61VdParameters.GetJuraArmOffset().Y(), na61VdParameters.GetJuraArmOffset().Z(), na61VdParameters.GetJuraArmRotX(), na61VdParameters.GetJuraArmRotY(), na61VdParameters.GetJuraArmRotZ());

  fHistFile = new TFile(fHistosPath.c_str(), "recreate");
  const std::string armNames[2] = {"Saleve", "Jura"};

  // Loop over each arm and create standalone event processing module objects
  for (unsigned char arm = 0; arm < 2; arm++) {
    /*
       na61FrameMergingModule[arm] = new Na61FrameMergingModule((armNames[arm] + " Frame Merging Module").c_str(), "Frame Merging Module");
       na61FrameMergingModule[arm]->SetCutOnNoisyPixels(1);
       na61FrameMergingModule[arm]->SetRunId(run);
       na61FrameMergingModule[arm]->SetToJuraArm(arm);
       na61FrameMergingModule[arm]->SetNewMerge(isSim ? 2 : 1); //2=simulations
       na61FrameMergingModule[arm]->SetRowOverlap(15);
       na61FrameMergingModule[arm]->SetHomePath(".");
       na61FrameMergingModule[arm]->Init();
       na61FrameMergingModule[arm]->DefineHistograms();
       */

    hotpixelsModule[arm] = new Na61HotPixelsModule(Form("HotPixels Module jura: %d",arm),"HotPixels Module");
    hotpixelsModule[arm]->SetRunId(run);
    hotpixelsModule[arm]->SetCutOnNoisyPixels(0.5); // 0.5% frequency
    hotpixelsModule[arm]->SetHomePath(".");
    hotpixelsModule[arm]->SetToJuraArm((bool)arm);
    hotpixelsModule[arm]->Init();
    hotpixelsModule[arm]->DefineHistograms();

    na61HitProducerModule[arm] = new Na61AlHitProducerModule((armNames[arm] + " Hit Producer Module").c_str(), "Hit Producer Module");
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

void AlVdHitFinderPS::End() {
  // Loop over each arm and clean up standalone event processing modules
  for (unsigned char arm = 0; arm < 2; arm++) {
    /*
       if (na61FrameMergingModule[arm]) {
       na61FrameMergingModule[arm]->Finish();
       delete na61FrameMergingModule[arm];
       na61FrameMergingModule[arm] = NULL;
       }*/
    if (hotpixelsModule[arm]) {
      if(fNoiseRun)hotpixelsModule[arm]->Finish();
      delete hotpixelsModule[arm];
      hotpixelsModule[arm] = NULL;
    }

    if (na61HitProducerModule[arm]) {
      if(fNoiseRun)na61HitProducerModule[arm]->Finish();
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
/*
   void VdHitFinderPS::LocalToShineGlobal(unsigned char arm, unsigned char sensor, USensorHit &hit, utl::Point &shinePos) {
   Na61VdParameters &vdParams = *(Na61VdParametersManager::Instance())->GetVdParams();
   Na61ArmParameters &armParams = arm ? *Na61VdParametersManager::Instance()->GetJuraArmParams() : *Na61VdParametersManager::Instance()->GetSaleveArmParams();

// initial arm params in vd-standalone coordinate system
double rotz0 = armParams.GetRotZ(sensor);
double roty0 = armParams.GetRotY(sensor);
double rotx0 = armParams.GetRotX(sensor);
const double VolumeX0 = armParams.GetVolumeX(sensor);
const double VolumeY0 = armParams.GetVolumeY(sensor);
double VolumeZ0 = armParams.GetVolumeZ(sensor);

double rotx, roty, rotz;

if (arm) {
rotx = vdParams.GetJuraArmRotX();
roty = vdParams.GetJuraArmRotY();
rotz = vdParams.GetJuraArmRotZ();
} else {
rotx = vdParams.GetSaleveArmRotX();
roty = vdParams.GetSaleveArmRotY();
rotz = vdParams.GetSaleveArmRotZ();
}

const double dx = vdParams.GetJuraArmOffset().X();
const double dy = vdParams.GetJuraArmOffset().Y();
const double dz = vdParams.GetJuraArmOffset().Z();

const double sz = TMath::Sin(rotz);
const double cz = TMath::Cos(rotz);
const double sy = TMath::Sin(roty);
const double cy = TMath::Cos(roty);
const double sx = TMath::Sin(rotx);
const double cx = TMath::Cos(rotx);

VolumeZ0 -= 75.0;

const double VolumeX1 = (cz * cy - sz * sx * sy) * VolumeX0 + sz * cx * VolumeY0 + (cz * sy + sz * sx * cy) * VolumeZ0;
const double VolumeY1 = -(sz * cy + cz * sx * sy) * VolumeX0 + cz * cx * VolumeY0 + (-sz * sy + cz * sx * cy) * VolumeZ0;
double VolumeZ1 = -cx * sy * VolumeX0 - sx * VolumeY0 + cx * cy * VolumeZ0;

VolumeZ1 += 75.0;

// vd params in SHINE coordinate system
const double alpha = rotz0;
const double beta = roty0 + roty;
const double gamma = rotx0 + rotx;
const double VolumeX = VolumeX1 + dx;
const double VolumeY = VolumeY1 + dy;
const double VolumeZ = VolumeZ1 + dz;

const double sa = TMath::Sin(alpha);
const double ca = TMath::Cos(alpha);
const double sb = TMath::Sin(beta);
const double cb = TMath::Cos(beta);
const double sg = TMath::Sin(gamma);
const double cg = TMath::Cos(gamma);

const double x1 = hit.GetX();
const double y1 = hit.GetY();
const double z1 = 0;

// ZXY in SHINE system
const double x2 = (ca * cb - sa * sg * sb) * x1 + sa * cg * y1 + (ca * sb + sa * sg * cb) * z1;
const double y2 = -(sa * cb + ca * sg * sb) * x1 + ca * cg * y1 + (-sa * sb + ca * sg * cb) * z1;
const double z2 = -cg * sb * x1 - sg * y1 + cg * cb * z1;

// Add shift relative to TPC
const double fOffsetX = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcX();
const double fOffsetY = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcY();
const double fOffsetZ = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcZ();

if (arm)
shinePos.Set(0.1 * (x2 + VolumeX + fOffsetX), 0.1 * (y2 + VolumeY + fOffsetY), 0.1 * (z2 + VolumeZ + fOffsetZ));
else
shinePos.Set(-0.1 * (x2 + VolumeX + fOffsetX), 0.1 * (y2 + VolumeY + fOffsetY), 0.1 * (z2 + VolumeZ + fOffsetZ));
return;
}

void VdHitFinderPS::LocalToShineGlobal(unsigned char arm, unsigned char sensor, USensorHit &hit, utl::Point &shinePos) {
Na61VdParameters &vdParams = *(Na61VdParametersManager::Instance())->GetVdParams();
Na61ArmParameters &armParams = arm ? *Na61VdParametersManager::Instance()->GetJuraArmParams() : *Na61VdParametersManager::Instance()->GetSaleveArmParams();

//initial arm params in vd-standalone coordinate system
double rotz0 = armParams.GetRotZ(sensor);
double roty0 = armParams.GetRotY(sensor);
double rotx0 = armParams.GetRotX(sensor);

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
//  const double sz = TMath::Sin(rotz);
//  const double cz = TMath::Cos(rotz);
//  const double sy = TMath::Sin(roty);
//  const double cy = TMath::Cos(roty);
//  const double sx = TMath::Sin(rotx);
//  const double cx = TMath::Cos(rotx);

const double x0 = hit.GetX();
const double y0 = hit.GetY();
double z0 = hit.GetZ();
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

//  if (arm)
shinePos.Set(0.1 * (x1 + fOffsetX), 0.1 * (y1 + fOffsetY), 0.1 * (z1 + fOffsetZ));
//  else
//   shinePos.Set(-0.1 * (x1 + fOffsetX), 0.1 * (y1 + fOffsetY), 0.1 * (z1 + fOffsetZ));

return;
}
*/


//_____________________________________________________________________
void AlVdHitFinderPS::LocalToGlobal(unsigned char arm, unsigned char sensor, utl::Point &localPos,utl::Point &vdstandalonePos) {
//
Na61AlVdArmParameters &armParams = arm ? *Na61VdParametersManager::Instance()->GetJuraAlVdArmParams() : *Na61VdParametersManager::Instance()->GetSaleveAlVdArmParams();

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


void AlVdHitFinderPS::LocalToShineGlobal(unsigned char arm, unsigned char /* sensor */, utl::Point &vdstandalonePos, utl::Point &shinePos) {
//
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
