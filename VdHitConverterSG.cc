#include "VdHitConverterSG.h"
#include "UDataTable.h"
#include "USensorPixel.h"

#include <fwk/CentralConfig.h>
#include <fwk/RandomEngineRegistry.h>

#include <evt/raw/Silicon.h>
#include <evt/BOSBank.h>
#include <evt/BOSRecord.h>

#include <evt/SimEvent.h>
#include <evt/Event.h>
#include <evt/raw/Silicon.h>
#include <evt/proc/Silicon.h>

#include <det/Silicon.h>
#include <det/Detector.h>

#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>

#include <TRandom3.h>

#include <sstream>


namespace VdHitConverterSG {

fwk::VModule::EResultFlag
VdHitConverterSG::Init()
{
  utl::Branch topBranch =
    fwk::CentralConfig::GetInstance().GetTopBranch("VdReconstructionSG");

  InitVerbosity(topBranch);

  return eSuccess;
}

fwk::VModule::EResultFlag
VdHitConverterSG::Process(evt::Event &e, const utl::AttributeMap &)
{
  const evt::EventHeader &eventHeader = e.GetEventHeader();
  if (!eventHeader.IsInitialized()) {
    WARNING("No EventHeader is present in Event. This is a serious problem!!! "
            "PLEASE DO CHECK WHAT IS GOING ON IF THIS EFFECT IS NOT EXPECTED!");
    return eFailure;
  }

  det::Detector &detector = det::Detector::GetInstance();
  detector.Update(eventHeader.GetTime(), eventHeader.GetRunNumber());

  const det::Silicon &detSilicon = detector.GetSilicon();
  if (not detSilicon.HasDetailedGeometry()) {
    ERROR("No access to VD geometry in local CS. Looks like fixed manager was"
          " used. Please make sure that SiliconGeometryXMLManager is called"
          " instead.");
    return eFailure;
  }
  const det::SiliconGeometry &vdGeom = detSilicon.GetDetailedGeometry();

  evt::SimEvent& simEvent = e.GetSimEvent();
  evt::proc::Silicon& procSilicon = e.GetProcEvent().GetSilicon();

  const double sensoreff = 1;
  //const double sensoreff[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; //perfect case
  //const double sensoreff[16] = {0.87,0.95,0.99,0.98,0.90,0.83,0.96,0.89,0.97,0.93,0.98,0.98,0.96,0.95,0.90,0.91};

  const double smearing = 0;
  //const double smearing = 0.0005;

  evt::proc::SimAndRecEventTranslator& translator =
    e.GetProcEvent().GetReconstruction().GetSimAndRecEventTranslator();

  typedef evt::Index<evt::sim::VertexTrack> SimTrackIndex;
  typedef evt::Index<evt::sim::Hit> SimHitIndex;
  std::map<int, std::map<SimTrackIndex, std::vector<SimHitIndex>>>
    sensorTrackZHits;

  int n = 0;
  for (auto track = simEvent.Begin<evt::sim::VertexTrack>();
      track != simEvent.End<evt::sim::VertexTrack>(); track++)
  {
    for (auto hitIt = track->HitsBegin(); hitIt != track->HitsEnd(); ++hitIt)
    {
      if (!simEvent.Has(*hitIt))
        continue;

      const evt::sim::Hit& hit = simEvent.Get(*hitIt);
      if (hit.GetDetectorId() != det::Const::eVD)
        continue;
      n++;

      sensorTrackZHits[hit.GetSubDetectorId()][track->GetIndex()]
        .push_back(hit.GetIndex());
    }
  }

  std::cout << n << " VD hits obtained" << std::endl;

  TRandom3* random = new TRandom3(0);

  for (auto sensor : sensorTrackZHits) {
    int insertHits=0;
    int trueHits=0;

    const det::SiliconSensor& detSensor = detSilicon.GetSensor(sensor.first);
    const utl::CoordinateSystemPtr& sensorCs = detSensor.GetComponentCoordinateSystem();
    for (auto trkIdxAndHits : sensor.second) {
      trueHits++;

      const evt::sim::Hit& hitStart = simEvent.Get(trkIdxAndHits.second.front());
      const auto positionStart = hitStart.GetPosition().GetCoordinates(sensorCs);

      const double x_init = positionStart.get<0>();
      const double y_init = positionStart.get<1>();
      const double x = random->Gaus(x_init, smearing) /* [SHINE units (cm)] */;
      const double y = random->Gaus(y_init, smearing) /* [SHINE units (cm)] */;

      const det::SiliconConst::EArm sensorArm =
        detSilicon.GetSensorArm(detSensor.GetSubdetectorNumber());

      utl::Point localPos = {x, y, 0};

      const utl::Point vdstandalonePos =
        vdGeom.SensorCSToArmCS(detSensor.GetComponentId(), localPos);

      const utl::Point clusterShinePos =
        vdGeom.ArmCSToShineCS(sensorArm, vdstandalonePos);

      //cout<<"geant: "<<hitStart.GetPosition().GetX()<<" "<<hitStart.GetPosition().GetY()<<" "<<hitStart.GetPosition().GetZ()<<endl;
      //cout<<"convt: "<<clusterShinePos.GetX()<<" "<<clusterShinePos.GetY()<<" "<<clusterShinePos.GetZ()<<endl;
      double distX = hitStart.GetPosition().GetX() - clusterShinePos.GetX();
      double distY = hitStart.GetPosition().GetY() - clusterShinePos.GetY();
      double distZ = hitStart.GetPosition().GetZ() - clusterShinePos.GetZ();
      double dist = sqrt(distX*distX + distY*distY + distZ*distZ);
      if ((distX > 0.01) || (distY > 0.01) || (distZ > 0.01)) {
        WARNING(Form("check geom! dist(hit - cluster) = %f = (%f, %f, %f)",
                     dist, distX, distY, distZ));
        WARNING(Form("hit = (%f, %f, %f)",
                     hitStart.GetPosition().GetX(),
                     hitStart.GetPosition().GetY(),
                     hitStart.GetPosition().GetZ()));
        WARNING(Form("cluster = (%f, %f, %f)",
                     clusterShinePos.GetX(),
                     clusterShinePos.GetY(),
                     clusterShinePos.GetZ()));
      }

      if (random->Uniform(0.,1.) <= sensoreff) {
        evt::rec::Cluster &cluster = e.GetRecEvent().Make<evt::rec::Cluster>();
        cluster.SetPositionUncertainty(evt::rec::ClusterConst::eX, 5.0 * utl::micrometer);
        //cluster.SetPositionUncertainty(evt::rec::ClusterConst::eY, 2.5 * utl::micrometer);
        cluster.SetPositionUncertainty(evt::rec::ClusterConst::eY, 5.0 * utl::micrometer);
        cluster.SetNumberOfPixels(1);
        cluster.SetTPCId(detSensor.GetDetector());
        cluster.SetSubdetectorNumber(detSensor.GetSubdetectorNumber());
        evt::proc::SiliconSensor &procSensor =
          procSilicon.HasSensor(detSensor.GetSubdetectorNumber())
          ? procSilicon.GetSensor(detSensor.GetSubdetectorNumber())
          : procSilicon.MakeSensor(detSensor.GetSubdetectorNumber());
        procSilicon.GetClusterLocalPositions()[cluster.GetIndex()] = vdstandalonePos;
        cluster.SetPosition(clusterShinePos);

        procSilicon.GetClusters().insert(cluster.GetIndex());
        procSensor.GetClusters().insert(cluster.GetIndex());
        insertHits++;

        translator.AddClusterAndHitPair(cluster.GetIndex(), trkIdxAndHits.second.front());
      }

    }  // track loop
    std::cout
      << "VD sensor=" << detSensor.GetSubdetectorNumber()
      << " true hits=" << trueHits
      << " inserted hits=" << insertHits
      << std::endl;

  } //sensor loop

  return eSuccess;
}

fwk::VModule::EResultFlag
VdHitConverterSG::Finish()
{ return eSuccess; }

}
