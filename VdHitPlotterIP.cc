#include "VdHitPlotterIP.h"

#include <det/Detector.h>
#include <det/Silicon.h>

#include <evt/Event.h>
#include <evt/RecEvent.h>

#include <TFile.h>
#include <TCanvas.h>


fwk::VModule::EResultFlag
VdHitPlotterIP::Init()
{ return eSuccess; }

void
VdHitPlotterIP::InitializeHistograms()
{
  INFO("allocating histograms");

  const det::Silicon &silicon = det::Detector::GetInstance().GetSilicon();
  for (auto it = silicon.SensorsBegin(); it != silicon.SensorsEnd(); ++it) {
    const det::SiliconSensor &sensor = *it;
    const det::SiliconType &type = sensor.GetType();

    const std::string &id = sensor.GetComponentId();
    const int station = sensor.GetStation();
    const int element = sensor.GetElement();
    const int ladder = sensor.GetLadder();
    const double halfHeight = (type.GetLines()*type.GetLinePitch())/2.;
    const double halfWidth = (type.GetColumns()*type.GetColumnPitch())/2.;

    INFO(Form("sensor '%s' (station=%d, element=%d, ladder=%d)",
          id.c_str(), station, element, ladder));

    TH2D *h2SAVDPsPos = new TH2D {
      Form("savdPsPosHistos[%d][%d][%d]", station, element, ladder),
      Form("%s (PS transform); x [cm]; y [cm]", id.c_str()),
      200, -halfHeight, halfHeight,
      200, -halfWidth, halfWidth,
    };

    TH2D *h2MAVDPsPos = new TH2D {
      Form("mavdPsPosHistos[%d][%d][%d]", station, element, ladder),
      Form("%s (PS transform); x [cm]; y [cm]", id.c_str()),
      200, -halfHeight, halfHeight,
      200, -halfWidth, halfWidth,
    };

    TH2D *h2SensorPos = new TH2D {
      Form("sensorPosHistos[%d][%d][%d]", station, element, ladder),
      Form("%s (SensorCS); x [cm]; y [cm]", id.c_str()),
      200, -halfHeight, halfHeight,
      200, -halfWidth, halfWidth,
    };

    TH2D *h2ShinePos = new TH2D {
      Form("shinePosHistos[%d][%d][%d]", station, element, ladder),
      Form("%s (ShineCS+offset); x [cm]; y [cm]", id.c_str()),
      200, -halfHeight, halfHeight,
      200, -halfWidth, halfWidth,
    };

    TH2D *h2pix = new TH2D {
      Form("pixHistos[%d][%d][%d]", station, element, ladder),
      Form("%s; column; row", id.c_str()),
      200, 0, double(type.GetColumns()),
      200, 0, double(type.GetLines()),
    };

    fSAVDPsPosHistos[station][element][ladder] = h2SAVDPsPos;
    fMAVDPsPosHistos[station][element][ladder] = h2MAVDPsPos;
    fSensorPosHistos[station][element][ladder] = h2SensorPos;
    fShinePosHistos[station][element][ladder] = h2ShinePos;
    fPixHistos[station][element][ladder] = h2pix;
  }
}

fwk::VModule::EResultFlag
VdHitPlotterIP::Process(evt::Event &event, const utl::AttributeMap&)
{
  static bool isFirstTime = true;
  if (isFirstTime) {
    isFirstTime = false;
    InitializeHistograms();
  }

  const det::Detector &detector = det::Detector::GetInstance();
  const utl::CoordinateSystemPtr &shineCS =
    detector.GetDetectorCoordinateSystem();
  const det::Silicon &silicon = detector.GetSilicon();

  //const evt::RecEvent &recEvent = event.GetRecEvent();
  const evt::SimEvent &simEvent = event.GetSimEvent();

  for (auto it = simEvent.Begin<evt::sim::Hit>();
      it != simEvent.End<evt::sim::Hit>(); ++it)
  {
    double x0, y0, z0;
    double x, y, z;

    const evt::sim::Hit &hit = *it;
    if (hit.GetDetectorId() != det::Const::eVD)
      continue;

    const det::SiliconSensor &sensor =
      silicon.GetSensor(hit.GetSubDetectorId());
    const int station = sensor.GetStation();
    const int element = sensor.GetElement();
    const int ladder = sensor.GetLadder();

    boost::tie(x0, y0, z0) = sensor.GetCenterPosition().GetCoordinates(shineCS);
    boost::tie(x, y, z) = hit.GetPosition().GetCoordinates(shineCS);
    fShinePosHistos.at(station).at(element).at(ladder)->Fill(x-x0, y-y0);

    const utl::CoordinateSystemPtr &sensorCS =
      silicon.GetSensorCS(sensor.GetSubdetectorNumber());
    const utl::CoordinateSystemPtr &nativeCS =
      silicon.GetNativeSensorCS(sensor.GetSubdetectorNumber());

    const det::SiliconType &type = sensor.GetType();
    const double columnPitch = type.GetColumnPitch();
    const double linePitch = type.GetLinePitch();

    boost::tie(x, y, z) = hit.GetPosition().GetCoordinates(nativeCS);
    const int column = x/columnPitch;
    const int line = y/linePitch;
    fPixHistos.at(station).at(element).at(ladder)->Fill(column, line);

    double dx, dy;
    if (ladder < 0) /* saleve */ {
      dx = (288 - line) * linePitch - linePitch/2.0 - linePitch;
      dy = (column - 576) * columnPitch - columnPitch/2.0 + columnPitch;
    } else /* jura */ {
      dx = (line - 288) * linePitch - linePitch / 2.0 + linePitch;
      dy = (576 - column) * columnPitch + columnPitch / 2.0 - columnPitch;
    }
    fSAVDPsPosHistos.at(station).at(element).at(ladder)->Fill(dx, dy);

    if (ladder < 0) /* saleve */ {
      dx = (256-line) * linePitch + linePitch/2.0 + 0.6*1e-2*0;
      dy = (column - 512) * columnPitch - columnPitch/2.0;
    } else /* jura */ {
      dx = (line-256) * linePitch - linePitch/2.0 - 0.6*1e-2*0;
      dy = (column - 512) * columnPitch - columnPitch/2.0;
    }
    fMAVDPsPosHistos.at(station).at(element).at(ladder)->Fill(dx, dy);

    boost::tie(x, y, z) = hit.GetPosition().GetCoordinates(sensorCS);
    fSensorPosHistos.at(station).at(element).at(ladder)->Fill(x, y);
  }

  return eSuccess;
}

static void
_MakeHistos(const std::string &prefix, const VdHitPlotterIP::Histograms &histos)
{
  for (const auto &station : histos) {
    TCanvas *c =
      new TCanvas {Form("%s_station_%d", prefix.c_str(), station.first)};

    // Calculate dimensions of the box to fit all sensors
    // (N elements x N ladders)
    const size_t nelements = station.second.size();
    size_t nladders = 0;
    for (const auto &element : station.second)
      nladders = std::max(nladders, element.second.size());
    //const size_t nentries = nelements * nladders;
    c->Divide(nladders, nelements, 0, 0);

    int padcnt = 1;
    for (auto &element : station.second) {
      for (auto &ladder : element.second) {
        c->cd(padcnt++);
        //gPad->SetTickx(2);
        //gPad->SetTicky(2);
        TH2D *hist = ladder.second;
        hist->Draw("COLZ");
      }
    }

    c->Write();
  }
}

fwk::VModule::EResultFlag
VdHitPlotterIP::Finish()
{
  INFO("writing histograms");

  TFile *ofile = new TFile {"VdHitPlotterIP.root", "RECREATE"};

  _MakeHistos("savdPsPosHistos", fSAVDPsPosHistos);
  _MakeHistos("mavdPsPosHistos", fMAVDPsPosHistos);
  _MakeHistos("sensorPosHistos", fSensorPosHistos);
  _MakeHistos("shinePosHistos", fShinePosHistos);
  _MakeHistos("pixHistos", fPixHistos);

  ofile->Write();
  ofile->Close();

  return eSuccess;
}

