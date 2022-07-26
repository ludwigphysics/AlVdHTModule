#ifndef _VdTrackFinderSG_VdTrackFinderSG_h_
#define _VdTrackFinderSG_VdTrackFinderSG_h_

/**
 * \file
 * \author Anastasia Merzlaya
 * \date 31 Jul 2018
 */

#ifndef ROOT_TFile
#include "TFile.h"
#endif

#include <fwk/VModule.h>
#include <boost/utility.hpp>

class Na61VdTrackingInitModule;
class Na61PrimaryVertexRecoModule;
class Na61VdTrackingHTModule;
class Na61PrimaryVertexRecoHTModule;
class Na61VdTrackingHTPackage;
class Na61VdTpcMatchingModule;


namespace VdTrackFinderSG {

/**
\class VdTrackFinderSG

\brief VD tracking module

\author Anastasia Merzlaya
\date 31 Jul 2018
\version
\ingroup TODO: Put a group here.
*/

class VdTrackFinderSG : public boost::noncopyable, public fwk::VModule {
 public:
  fwk::VModule::EResultFlag Init();
  fwk::VModule::EResultFlag Process(evt::Event &event, const utl::AttributeMap &attr);
  fwk::VModule::EResultFlag Finish();

 private:
  void Begin(const unsigned int run, bool isSim);
  void End();
  unsigned int prevRun;
  unsigned int productionMode;
  Na61VdTrackingInitModule *na61VdTrackingInitModule[2];
  Na61PrimaryVertexRecoModule *na61PrimaryVertexRecoModule;
  Na61VdTrackingHTModule *na61VdTrackingHTModule;
  Na61PrimaryVertexRecoHTModule *na61PrimaryVertexRecoHTModule;
  Na61VdTrackingHTPackage *na61VdTrackingHTPackage;
  Na61VdTpcMatchingModule *na61VdTpcMatchingModule;
  Na61PrimaryVertexRecoHTModule *na61PrimaryVertexRecoHTModuleTPC;
  std::string fMatchParamsPath;
  bool fDoEventCuts;
  std::string fHistosPath;

  TFile *fHistFile;
  REGISTER_MODULE("VdTrackFinderSG", VdTrackFinderSG, "$Id$");
};
}
#endif  // _VdTrackFinderSG_VdTrackFinderSG_h_
