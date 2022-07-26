#include "AlVdHitConverterSG.h"
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
#include <evt/SimEvent.h>
#include <evt/Event.h>
#include <evt/raw/Silicon.h>
#include <evt/proc/Silicon.h>
#include <det/Silicon.h>
#include <fwk/CentralConfig.h>
#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>
#include <fwk/RandomEngineRegistry.h>
#include <TRandom3.h>

#include <sstream>

using namespace std;
using namespace fwk;
using namespace utl;
using namespace evt;

namespace AlVdHitConverterSG {

fwk::VModule::EResultFlag AlVdHitConverterSG::Init() {
  utl::Branch topBranch = fwk::CentralConfig::GetInstance().GetTopBranch("VdReconstructionSG");
  InitVerbosity(topBranch);

  prevRun = 0;
  return eSuccess;
}

fwk::VModule::EResultFlag AlVdHitConverterSG::Process(evt::Event &e, const utl::AttributeMap &) {
  const evt::EventHeader &eventHeader = e.GetEventHeader();
  if (!eventHeader.IsInitialized()) {
    WARNING("No EventHeader is present in Event. This is a serious problem!!! PLEASE DO CHECK WHAT IS GOING ON IF THIS EFFECT IS NOT EXPECTED!!!");
    return eFailure;
  }
  det::Detector::GetInstance().Update(eventHeader.GetTime(), eventHeader.GetRunNumber());

  evt::SimEvent& simEvent = e.GetSimEvent();

  //const evt::RecEvent &recEvent = e.GetRecEvent();
  //if (!recEvent.HasMainVertex())
    //return eSuccess;
  //if (!recEvent.HasPrimaryVertex(evt::rec::VertexConst::ePrimaryVD))
    //return eSuccess;


  const double sensoreff[68] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
				1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; //perfect case
  //const double sensoreff[16] = {0.87,0.95,0.99,0.98,0.90,0.83,0.96,0.89,0.97,0.93,0.98,0.98,0.96,0.95,0.90,0.91};
  //const double smearing = 0.0005;
  const double smearing = 0;
  if (eventHeader.GetRunNumber() != prevRun) {
    Begin(eventHeader.GetRunNumber());
    prevRun = eventHeader.GetRunNumber();
  }

  //Get simAndRecEventTranslator object.
  proc::SimAndRecEventTranslator& translator = e.GetProcEvent().GetReconstruction().GetSimAndRecEventTranslator();


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
    /*{
      {0, 0}, {2, 1}, {4, 2}, {5, 3}, {10, 4}, {11, 5}, {8, 6}, {9, 7}, 
      {1, 0}, {3, 1}, {6, 2}, {7, 3}, {12, 4}, {13, 5}, {14, 6}, {15, 7}
      }; */
    {
      {0, 0}, {1, 1}, {2, 2}, {6, 3}, {7, 4}, {8, 5}, {9, 6}, {10, 7}, {11, 8}, {18, 9}, {19,10}, {20,11}, {21,12}, {22,13}, {23,14}, 
      {24,15}, {25,16}, {26,17}, {27,18}, {38,19}, {39,20}, {40,21}, {41,22}, {42,23}, {43,24}, {44,25}, {45,26}, {46,27}, {47,28}, {48,29}, 
      {49,30}, {50,31}, {51,32}, {52,33},  

      {3, 0}, {4, 1}, {5, 2}, {12, 3}, {13, 4}, {14, 5}, {15, 6}, {16, 7}, {17, 8}, {28, 9}, {29,10}, {30,11}, {31,12}, {32,13}, {33,14},
      {34,15}, {35,16}, {36,17}, {37,18}, {53,19}, {54,20}, {55,21}, {56,22}, {57,23}, {58,24}, {59,25}, {60,26}, {61,27}, {62,28}, {63,29}, 
      {64,30}, {65,31}, {66,32}, {67,33}      
  };

  std::map<unsigned char, std::map<evt::Index<evt::sim::VertexTrack>, std::map<double, evt::Index<evt::sim::Hit> > > > sensorTrackZHits;
  //		for(auto sens=detector.GetSilicon().GetSensors().begin(); sens!=detector.GetSilicon().GetSensors().end(); sens++){
  //		  std::cout<<sens->GetComponentId()<<" "<<sens->GetCenterPosition()<<std::endl;
  //		}


  //for (auto it = simEvent.Begin<evt::sim::Hit>(); it != simEvent.End<evt::sim::Hit>(); ++it) {
  //const evt::sim::Hit &hit = *it;
  //std::cout << hit.GetPosition().GetX();
  //}

  int n = 0;
  for (auto track = simEvent.Begin<evt::sim::VertexTrack>();
      track != simEvent.End<evt::sim::VertexTrack>(); track++)
  {
    // if (distance(track->HitsEnd(),track->HitsBegin())==0) continue;
    for (auto hitIt = track->HitsBegin(); hitIt != track->HitsEnd(); hitIt++)
    {
      if (!simEvent.Has(*hitIt))
        continue;

      const evt::sim::Hit& hit = simEvent.Get(*hitIt);
      if (hit.GetDetectorId() != det::Const::eVD)
        continue;
      n++;
      //cout<<"trackId="<<hit.GetIndex()<<" z="<<hit.GetPosition().GetZ()<<" y="<<hit.GetPosition().GetY()<<" x="<<hit.GetPosition().GetX()<<" E="<<hit.GetEnergyDeposit()<<" SubDetId="<<hit.GetSubDetectorId()<<" "<<shineStandaloneName[hit.GetSubDetectorId()]<<endl;
      //Int_t ii;
      // cin>>ii;
      sensorTrackZHits[hit.GetSubDetectorId()][track->GetIndex()][hit.GetPosition().GetZ()] = hit.GetIndex();
    }
  }

  std::cout << n << " VD hits obtained" << std::endl;

  TRandom3* random = new TRandom3(0);

  for (auto sensor : sensorTrackZHits) {
    int insertHits=0;
    int trueHits=0;

    std::map<unsigned short, std::map<unsigned short, double> > pixelEnergy;
    const det::SiliconSensor& detSensor = det::Detector::GetInstance().GetSilicon().GetSensor(sensor.first);
    const utl::CoordinateSystemPtr& detCs = detSensor.GetComponentCoordinateSystem();
    for (auto sensorTrack : sensor.second) {
      trueHits++;

      const evt::sim::Hit& hitStart = simEvent.Get(sensorTrack.second.begin()->second);
      const auto positionStart = hitStart.GetPosition().GetCoordinates(detCs);

      // FIXME: these formulas are wrong for MAVD
      const double x_init = positionStart.get<0>();
      const double y_init = positionStart.get<1>();
      const double x = random->Gaus(x_init, smearing);
      const double y = random->Gaus(y_init, smearing);
      // check that hit stays inside the active sensor frame:

      unsigned char arm;
      if (detSensor.GetLadder() < 0)
        arm = 0;
      else
        arm = 1; 
      unsigned char s = shineStandaloneName[detSensor.GetSubdetectorNumber()];
      //cout<<(unsigned int)detSensor.GetSubdetectorNumber()<<" -> sensor="<<(unsigned int)s<<" arm="<<(unsigned int)arm
      //<<"  nameShine: "<<nameShine[arm][s]<<endl;
      TString name(nameShine[arm][s]);

      utl::Point localPos;
      localPos.Set(x*10.,y*10.,0.0);

      utl::Point armstandalonePos;
      LocalToArmGlobal(arm, s, localPos, armstandalonePos); 

      utl::Point vdstandalonePos;
      LocalToVdGlobal(arm, armstandalonePos, vdstandalonePos);

      utl::Point clusterShinePos;
      LocalToShineGlobal(arm, vdstandalonePos, clusterShinePos);

      fhZX_fine->Fill(vdstandalonePos.GetZ(), vdstandalonePos.GetX());
      fhZX_fine_G4->Fill(clusterShinePos.GetZ(), clusterShinePos.GetX());
      if(name.Contains("1_")){fhX_Al1->Fill(vdstandalonePos.GetX()); fhX_Al1_G4->Fill(clusterShinePos.GetX());}
      if(name.Contains("2_")){fhX_Al2->Fill(vdstandalonePos.GetX());  fhX_Al2_G4->Fill(clusterShinePos.GetX());}
      if(name.Contains("3_")){fhX_Al3->Fill(vdstandalonePos.GetX());  fhX_Al3_G4->Fill(clusterShinePos.GetX());}
      if(name.Contains("4_")){fhX_Al4->Fill(vdstandalonePos.GetX());  fhX_Al4_G4->Fill(clusterShinePos.GetX());}
      if(name.Contains("j1_")){fhY_Al1_j->Fill(vdstandalonePos.GetY()); fhY_Al1_j_G4->Fill(clusterShinePos.GetY());}
      if(name.Contains("j2_")){fhY_Al2_j->Fill(vdstandalonePos.GetY()); fhY_Al2_j_G4->Fill(clusterShinePos.GetY());}
      if(name.Contains("j3_")){fhY_Al3_j->Fill(vdstandalonePos.GetY()); fhY_Al3_j_G4->Fill(clusterShinePos.GetY());}
      if(name.Contains("j4_")){fhY_Al4_j->Fill(vdstandalonePos.GetY()); fhY_Al4_j_G4->Fill(clusterShinePos.GetY());}
      if(name.Contains("s1_")){fhY_Al1_s->Fill(vdstandalonePos.GetY()); fhY_Al1_s_G4->Fill(clusterShinePos.GetY());}
      if(name.Contains("s2_")){fhY_Al2_s->Fill(vdstandalonePos.GetY()); fhY_Al2_s_G4->Fill(clusterShinePos.GetY());}
      if(name.Contains("s3_")){fhY_Al3_s->Fill(vdstandalonePos.GetY()); fhY_Al3_s_G4->Fill(clusterShinePos.GetY());}
      if(name.Contains("s4_")){fhY_Al4_s->Fill(vdstandalonePos.GetY()); fhY_Al4_s_G4->Fill(clusterShinePos.GetY());}

      //cout<<"geant: "<<hitStart.GetPosition().GetX()<<" "<<hitStart.GetPosition().GetY()<<" "<<hitStart.GetPosition().GetZ()<<endl;
      //cout<<"convt: "<<clusterShinePos.GetX()<<" "<<clusterShinePos.GetY()<<" "<<clusterShinePos.GetZ()<<endl;
      //Int_t ii;
      //cin>>ii;
      double distX = hitStart.GetPosition().GetX()-clusterShinePos.GetX(); 
      double distY = hitStart.GetPosition().GetY()-clusterShinePos.GetY(); 
      double distZ = hitStart.GetPosition().GetZ()-clusterShinePos.GetZ(); 
      double dist = sqrt(distX*distX + distY*distY + distZ*distZ);
      if (distX > 0.01 || distY > 0.01 || distZ > 0.01) {
        WARNING(Form("check geom! dist(hit - cluster) = %f = (%f, %f, %f)", dist, distX, distY, distZ));
        WARNING(Form("hit = (%f, %f, %f)", hitStart.GetPosition().GetX(), hitStart.GetPosition().GetY(), hitStart.GetPosition().GetZ()));
        WARNING(Form("cluster = (%f, %f, %f)", clusterShinePos.GetX(), clusterShinePos.GetY(), clusterShinePos.GetZ()));
      }
      

      const det::Silicon &detSilicon = det::Detector::GetInstance().GetSilicon();
      const det::SiliconSensor &detSensor =
        detSilicon.GetSensor(detSilicon.SensorIdToSubdetectorNumber(nameShine[arm][s]));

      if (random->Uniform(0.,1.)<=sensoreff[(unsigned int)detSensor.GetSubdetectorNumber()]) {
        evt::rec::Cluster &cluster = e.GetRecEvent().Make<evt::rec::Cluster>();
        cluster.SetPositionUncertainty(evt::rec::ClusterConst::eX, 5.0 * utl::micrometer);
        //cluster.SetPositionUncertainty(evt::rec::ClusterConst::eY, 2.5 * utl::micrometer);
        cluster.SetPositionUncertainty(evt::rec::ClusterConst::eY, 5.0 * utl::micrometer);
        cluster.SetNumberOfPixels(1);
        cluster.SetTPCId(detSensor.GetDetector());
        cluster.SetSubdetectorNumber(detSensor.GetSubdetectorNumber());
        evt::proc::SiliconSensor &procSensor =
          e.GetProcEvent().GetSilicon().HasSensor(detSensor.GetSubdetectorNumber())
          ? e.GetProcEvent().GetSilicon().GetSensor(detSensor.GetSubdetectorNumber())
          : e.GetProcEvent().GetSilicon().MakeSensor(detSensor.GetSubdetectorNumber());
        e.GetProcEvent().GetSilicon().GetClusterLocalPositions()[cluster.GetIndex()] =
          utl::Point(0.1 * localPos.GetX(), 0.1 * localPos.GetY(), 0.1 * localPos.GetZ());
        cluster.SetPosition(clusterShinePos);

        e.GetProcEvent().GetSilicon().GetClusters().insert(cluster.GetIndex());
        procSensor.GetClusters().insert(cluster.GetIndex());
        insertHits++;


        translator.AddClusterAndHitPair(cluster.GetIndex(),sensorTrack.second.begin()->second);
        Index<sim::Hit> hitIndex;
        try {
          hitIndex = translator.GetAssociatedHitIndex(cluster.GetIndex());
        }
        catch (exception& exc) {
          WARNING("unmatched cluster");
        }	
      }

    }  // track loop
    //cout<<"VD sensor="<<(unsigned int)detSensor.GetSubdetectorNumber()<<" true hits="<<trueHits<<" inserted hits="<<insertHits<<endl;
    //Int_t ii;
    //cin>>ii;
  } //sensor loop

  return eSuccess;
}

fwk::VModule::EResultFlag AlVdHitConverterSG::Finish() {
  End();
  return eSuccess;
}

void AlVdHitConverterSG::Begin(const unsigned int run) {
  End();

  const unsigned int calibrun = Na61VdParameters::FindCalibRuns(run);
  //TFile &matchParamsFile = *TFile::Open(Form("matchParams-%0d.root", calibrun), "READ");
  TFile & matchParamsFile = *TFile::Open(Form("/afs/cern.ch/user/a/amerzlay/public/XeLa150/VDmatchingfiles/matchParams-%0d.root",calibrun),"READ");
  Na61VdParameters &na61VdParameters = *Na61VdParametersManager::Instance()->GetVdParams();
  na61VdParameters.SetRunId(run);
  na61VdParameters.SetMatchParamsFile(&matchParamsFile);
  na61VdParameters.SetdOffZ(0);
  na61VdParameters.SetN(0);
  na61VdParameters.Init();
  matchParamsFile.Close();

  Na61AlVdArmParameters &na61ArmParametersSaleve = *Na61VdParametersManager::Instance()->GetSaleveAlVdArmParams();
  na61ArmParametersSaleve.SetRunId(run);
  na61ArmParametersSaleve.Init();
  na61ArmParametersSaleve.SetVolumesInGlobal(na61VdParameters.GetSaleveArmOffset().X(), na61VdParameters.GetSaleveArmOffset().Y(), na61VdParameters.GetSaleveArmOffset().Z(), na61VdParameters.GetSaleveArmRotX(), na61VdParameters.GetSaleveArmRotY(), na61VdParameters.GetSaleveArmRotZ());

  Na61AlVdArmParameters &na61ArmParametersJura = *Na61VdParametersManager::Instance()->GetJuraAlVdArmParams();
  na61ArmParametersJura.SetRunId(run);
  na61ArmParametersJura.Init();
  na61ArmParametersJura.SetVolumesInGlobal(na61VdParameters.GetJuraArmOffset().X(), na61VdParameters.GetJuraArmOffset().Y(), na61VdParameters.GetJuraArmOffset().Z(), na61VdParameters.GetJuraArmRotX(), na61VdParameters.GetJuraArmRotY(), na61VdParameters.GetJuraArmRotZ());

  fHistFile = new TFile(Form("AlVdHitConverterSG_run%d.root", run), "recreate");
  const std::string armNames[2] = {"Saleve", "Jura"};

  fhZX_fine_G4 = new TH2F("hZX_fine_G4","",1000,(-6306)/10.,(-6306+166)/10.,1000,-5,5);
  fhX_Al1_G4 =	new TH1F("hX_Al1_G4","",1000,-5,5);
  fhX_Al2_G4 =	new TH1F("hX_Al2_G4","",1000,-5,5);
  fhX_Al3_G4 =	new TH1F("hX_Al3_G4","",1000,-5,5);
  fhX_Al4_G4 =	new TH1F("hX_Al4_G4","",1000,-5,5);
  fhY_Al1_j_G4 =	new TH1F("hY_Al1_j_G4","",1000,-8,8);
  fhY_Al2_j_G4 =	new TH1F("hY_Al2_j_G4","",1000,-8,8);
  fhY_Al3_j_G4 =	new TH1F("hY_Al3_j_G4","",1000,-8,8);
  fhY_Al4_j_G4 =	new TH1F("hY_Al4_j_G4","",1000,-8,8);
  fhY_Al1_s_G4 =	new TH1F("hY_Al1_s_G4","",1000,-8,8);
  fhY_Al2_s_G4 =	new TH1F("hY_Al2_s_G4","",1000,-8,8);
  fhY_Al3_s_G4 =	new TH1F("hY_Al3_s_G4","",1000,-8,8);
  fhY_Al4_s_G4 =	new TH1F("hY_Al4_s_G4","",1000,-8,8);

  fhZX_fine =	new TH2F("hZX_fine","",1000,-6,160,1000,-50,50);
  fhX_Al1 =	new TH1F("hX_Al1","",1000,-50,50);
  fhX_Al2 =	new TH1F("hX_Al2","",1000,-50,50);
  fhX_Al3 =	new TH1F("hX_Al3","",1000,-50,50);
  fhX_Al4 =	new TH1F("hX_Al4","",1000,-50,50);
  fhY_Al1_j =	new TH1F("hY_Al1_j","",1000,-80,80);
  fhY_Al2_j =	new TH1F("hY_Al2_j","",1000,-80,80);
  fhY_Al3_j =	new TH1F("hY_Al3_j","",1000,-80,80);
  fhY_Al4_j =	new TH1F("hY_Al4_j","",1000,-80,80);
  fhY_Al1_s =	new TH1F("hY_Al1_s","",1000,-80,80);
  fhY_Al2_s =	new TH1F("hY_Al2_s","",1000,-80,80);
  fhY_Al3_s =	new TH1F("hY_Al3_s","",1000,-80,80);
  fhY_Al4_s =	new TH1F("hY_Al4_s","",1000,-80,80);


}

void AlVdHitConverterSG::End() {
  // Loop over each arm and clean up standalone event processing modules

  if (fHistFile) {
    fHistFile->Write();
    fHistFile->Close();
    delete fHistFile;
    fHistFile = NULL;
  }
}

//_____________________________________________________________________
void AlVdHitConverterSG::LocalToArmGlobal(unsigned char arm, unsigned char sensor,
    utl::Point &localPos, utl::Point &armstandalonePos)
{
  Na61AlVdArmParameters &armParams =
    arm ? *Na61VdParametersManager::Instance()->GetJuraAlVdArmParams()
        : *Na61VdParametersManager::Instance()->GetSaleveAlVdArmParams();

  //initial arm params in vd-standalone coordinate system
  double alpha = armParams.GetRotZ(sensor);
  double beta  = armParams.GetRotY(sensor);
  double gamma = armParams.GetRotX(sensor);

  const double VolumeX = armParams.GetVolumeX(sensor);
  const double VolumeY = armParams.GetVolumeY(sensor);
  const double VolumeZ = armParams.GetVolumeZ(sensor);

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
  x2 =  (ca*cb - sa*sg*sb)*x1 + sa*cg*y1 + ( ca*sb + sa*sg*cb)*z1;
  y2 = -(sa*cb + ca*sg*sb)*x1 + ca*cg*y1 + (-sa*sb + ca*sg*cb)*z1;
  z2 =              -cg*sb*x1 -    sg*y1 +               cg*cb*z1;

  armstandalonePos.Set(x2 + VolumeX, y2 + VolumeY, z2 + VolumeZ);

  //cout<<"  x1="<<x1<<"  y1="<<y1<<"  z1="<<z1<<endl;
  //cout<<"  x2="<<x2<<"  y2="<<y2<<"  z2="<<z2<<endl;
  //cout<<"  arm frame:  x2="<<x2+VolumeX<<"  y2="<<y2+VolumeY<<"  z2="<<z2+VolumeZ<<endl;


}
  
  //_____________________________________________________________________________________________  
  void AlVdHitConverterSG::LocalToVdGlobal(unsigned char arm,
					      utl::Point &armstandalonePos, utl::Point &vdstandalonePos)
  {
  Na61VdParameters &vdParams =
    *(Na61VdParametersManager::Instance())->GetVdParams();

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

  const double x0 = armstandalonePos.GetX();
  const double y0 = armstandalonePos.GetY();
  double z0 = armstandalonePos.GetZ();
  const double VolumeZ0 = 75.0;
  z0 = z0 - VolumeZ0;

  double x1 = (ca * cb - sa * sg * sb) * x0 + sa * cg * y0 + (ca * sb + sa * sg * cb) * z0;
  double y1 = -(sa * cb + ca * sg * sb) * x0 + ca * cg * y0 + (-sa * sb + ca * sg * cb) * z0;
  double z1 = -cg * sb * x0 - sg * y0 + cg * cb * z0;
  z1 = z1 + VolumeZ0;

  x1 = x1 + dx;
  y1 = y1 + dy;
  z1 = z1 + dz;

  vdstandalonePos.Set(x1, y1, z1);
  return;
}

  //_____________________________________________________________________________________________  
  void AlVdHitConverterSG::LocalToShineGlobal(unsigned char /* arm */,
					      utl::Point &vdstandalonePos, utl::Point &shinePos)
  {

  // Add shift relative to TPC
  const double fOffsetX = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcX();
  const double fOffsetY = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcY();
  const double fOffsetZ = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcZ();


  const double x = vdstandalonePos.GetX();
  const double y = vdstandalonePos.GetY();
  double       z = vdstandalonePos.GetZ();

  //  if (arm)
  shinePos.Set(0.1 * (x + fOffsetX), 0.1 * (y + fOffsetY), 0.1 * (z + fOffsetZ));
  //  else
  //   shinePos.Set(-0.1 * (x1 + fOffsetX), 0.1 * (y1 + fOffsetY), 0.1 * (z1 + fOffsetZ));

  /* cout<<"Arm Rotations:"<<endl;
     if (arm)cout<<"Jura "<<rotx<<" "<<roty<<" "<<rotz<<endl;
     if (!arm)cout<<"Saleve "<<rotx<<" "<<roty<<" "<<rotz<<endl;
     cout<<"Arm Offsets (cm):"<<endl;
     if (arm)cout<<"Jura "<<dx*0.1<<" "<<dy*0.1<<" "<<dz*0.1<<endl;
     if (!arm)cout<<"Saleve "<<dx*0.1<<" "<<dy*0.1<<" "<<dz*0.1<<endl;
     cout<<"Shift to TPC (cm):"<<endl;
     cout<<OffsetX*0.1<<" "<<fOffsetY*0.1<<" "<<fOffsetZ*0.1<<endl;
     */
  return;
}


}
