#ifndef _VdHitConverterSG_VdHitConverterSG_h_
#define _VdHitConverterSG_VdHitConverterSG_h_

/**
 * \file
 * \author Anastasia Merzlaya
 * \date 27 Jun 2020
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
#include "Na61ArmParameters.h"
#include "Na61VdParametersManager.h"


namespace VdHitConverterSG {

/**
 * \class VdHitConverterSG
 *
 * \brief Describe your module. In one sentence.
 *
 * Now here a longer description in doxygen. You may use HTML.
 *
 * \author Anastasia Merzlaya, Ivan Pidhurskyi <ivan.pidhurskyi@cern.ch>
 * \date 27 Jun 2020, 03 Jun 2022
 */
class VdHitConverterSG: public boost::noncopyable, public fwk::VModule {
 public:
  fwk::VModule::EResultFlag
  Init();

  fwk::VModule::EResultFlag
  Process(evt::Event &event, const utl::AttributeMap &attr);

  fwk::VModule::EResultFlag
  Finish();

 private:
  TFile *fHistFile;

  // This goes at the end.
  REGISTER_MODULE("VdHitConverterSG", VdHitConverterSG, "$Id$");
};

}
#endif  // _VdHitConverterSG_VdHitConverterSG_h_
