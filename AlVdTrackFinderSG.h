#ifndef _VdTrackFinderSG_AlVdTrackFinderSG_h_
#define _VdTrackFinderSG_AlVdTrackFinderSG_h_

/**
 * \file
 * \author Anastasia Merzlaya
 * \date 31 Jul 2018
 */

#ifndef ROOT_TFile
#include "TFile.h"
#endif
#include "Na61AlVdTrackingInitModule.h"


#include <fwk/VModule.h>
#include <boost/utility.hpp>

class Na61AlVdTrackingInitModule;
class Na61AlPrimaryVertexRecoModule;
class Na61AlVdTrackingHTModule;
class Na61PrimaryVertexRecoHTModule;
class Na61AlVdTrackingHTPackage;
class Na61VdTpcMatchingModule;
class Na61VdIOModule;


namespace AlVdTrackFinderSG {

/**
\class VdTrackFinderSG

\brief VD tracking module

\author Anastasia Merzlaya
\date 31 Jul 2018
\version
\ingroup TODO: Put a group here.
*/

class AlVdTrackFinderSG : public boost::noncopyable, public fwk::VModule {
 public:
  fwk::VModule::EResultFlag Init();
  fwk::VModule::EResultFlag Process(evt::Event &event, const utl::AttributeMap &attr);
  fwk::VModule::EResultFlag Finish();

 private:
  void Begin(const unsigned int run, bool isSim);
  void End();
  unsigned int prevRun;
  unsigned int productionMode;
  Na61AlVdTrackingInitModule *na61AlVdTrackingInitModule[2];
  Na61AlPrimaryVertexRecoModule *na61PrimaryVertexRecoModule;
  Na61AlVdTrackingHTModule *na61AlVdTrackingHTModule;
  Na61PrimaryVertexRecoHTModule *na61PrimaryVertexRecoHTModule;
  Na61AlVdTrackingHTPackage *na61AlVdTrackingHTPackage;
  Na61VdTpcMatchingModule *na61VdTpcMatchingModule;
  Na61PrimaryVertexRecoHTModule *na61PrimaryVertexRecoHTModuleTPC;
  Na61VdIOModule *na61OutputModuleJ;
  Na61VdIOModule *na61OutputModuleS;
  std::string fMatchParamsPath;
  bool fDoEventCuts;
  std::string fHistosPath;
  UVdEvent eventArm[2];
  UVdEvent eventInput;


  TFile *fHistFile;
  REGISTER_MODULE("AlVdTrackFinderSG", AlVdTrackFinderSG, "$Id$");
};
}
#endif  // _AlVdTrackFinderSG_AlVdTrackFinderSG_h_
