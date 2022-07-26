#ifndef _AlVdHitFinderPS_AlVdHitFinderPS_h_
#define _AlVdHitFinderPS_AlVdHitFinderPS_h_

/**
 * \file
 * \author Pawel Staszel
 * \date 22 Feb 2018
 */

#ifndef ROOT_TFile
#include "TFile.h"
#endif

#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include "UVdEvent.h"
#include "USensorHit.h"
#include <utl/Point.h>
#include "TMath.h"

#include "Na61VdParameters.h"
#include "Na61AlVdArmParameters.h"
#include "Na61VdParametersManager.h"

class Na61HotPixelsModule;
class Na61AlHitProducerModule;

namespace AlVdHitFinderPS {

/**
\class VdHitFinderPS

\brief Describe your module. In one sentence.

Now here a longer description in doxygen. You may use HTML.

\author Pawel Staszel
\date 22 Feb 2018
\ingroup TODO: Put a group here.
*/

class AlVdHitFinderPS : public boost::noncopyable, public fwk::VModule {
 public:
  fwk::VModule::EResultFlag Init();
  fwk::VModule::EResultFlag Process(evt::Event &event, const utl::AttributeMap &attr);
  fwk::VModule::EResultFlag Finish();

 private:

  unsigned int fNoiseRun;
  unsigned int prevRun;
  void Begin(const unsigned int run, bool isSim);
  void End();
  //void LocalToShineGlobal(unsigned char arm, unsigned char sensor, USensorHit &hit, utl::Point &shinePos);
  void LocalToShineGlobal(unsigned char arm, unsigned char sensor, utl::Point &vdstandalonePos, utl::Point &shinePos);
  void LocalToGlobal(unsigned char arm, unsigned char sensor, utl::Point &localPos,utl::Point &vdstandalonePos);

  //Na61FrameMergingModule *na61FrameMergingModule[2];
  Na61HotPixelsModule *hotpixelsModule[2];
  Na61AlHitProducerModule *na61HitProducerModule[2];
  TFile *fHistFile;
  std::string fMatchParamsPath;
  std::string fHistosPath;


  // This goes at the end.
  REGISTER_MODULE("AlVdHitFinderPS", AlVdHitFinderPS, "$Id$");
};
}
#endif  // _VdHitFinderPS_VdHitFinderPS_h_
