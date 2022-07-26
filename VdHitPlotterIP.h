#ifndef _VdHitPlotterIP_h_
#define _VdHitPlotterIP_h_

#include <fwk/VModule.h>

#include <TH2D.h>

#include <map>


class VdHitPlotterIP: public fwk::VModule {
  public:
  using Histograms =
    std::map<
      int/*station*/,
      std::map<
        int/*element*/,
        std::map<
          int/*ladder*/,
          TH2D*
        >,
        std::greater<int>
      >
    >;

  fwk::VModule::EResultFlag
  Init() override;

  fwk::VModule::EResultFlag
  Process(evt::Event &event, const utl::AttributeMap &attr) override;

  fwk::VModule::EResultFlag
  Finish() override;

  private:
  void
  InitializeHistograms();

  private:
  Histograms fSAVDPsPosHistos,
             fMAVDPsPosHistos,
             fSensorPosHistos,
             fShinePosHistos,
             fPixHistos;

  REGISTER_MODULE("VdHitPlotterIP", VdHitPlotterIP, "$Id$");
}; // class VdHitPlotterIP

#endif
