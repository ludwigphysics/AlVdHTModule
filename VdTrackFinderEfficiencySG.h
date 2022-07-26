#ifndef _VdTrackFinderEfficiencySG_VdTrackFinderEfficiencySG_h_
#define _VdTrackFinderEfficiencySG_VdTrackFinderEfficiencySG_h_

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
class Na61VdEfficiencyModule;
class Na61PrimaryVertexRecoModule;
class Na61VdTrackingHTModule;
class Na61PrimaryVertexRecoHTModule;
class Na61VdTrackingHTPackage;
class Na61VdTpcMatchingModule;

class Na61VdTrackingPostModule;
class Na61PrimaryVertexRecoPostModule;

namespace VdTrackFinderEfficiencySG {

/**
\class VdTrackFinderEfficiencySG

\brief VD resonstruction module - testing sensor efficiency

\author Anastasia Merzlaya
\date 31 Jul 2018
\version
\ingroup TODO: Put a group here.
*/

class VdTrackFinderEfficiencySG : public boost::noncopyable, public fwk::VModule {
 public:
  fwk::VModule::EResultFlag Init();
  fwk::VModule::EResultFlag Process(evt::Event &event, const utl::AttributeMap &attr);
  fwk::VModule::EResultFlag Finish();

 private:
  void Begin(const unsigned int run);
  void End();
  unsigned int prevRun;
  unsigned int productionMode;
  Na61VdTrackingInitModule *na61VdTrackingInitModule[2];
  Na61VdEfficiencyModule *vdefficModule[2];
  Na61PrimaryVertexRecoModule *na61PrimaryVertexRecoModule;
  Na61VdTrackingHTModule *na61VdTrackingHTModule;
  Na61PrimaryVertexRecoHTModule *na61PrimaryVertexRecoHTModule;
  Na61VdTrackingHTPackage *na61VdTrackingHTPackage;
  Na61VdTpcMatchingModule *na61VdTpcMatchingModule;
  Na61PrimaryVertexRecoHTModule *na61PrimaryVertexRecoHTModuleTPC;
  
  Na61VdTrackingPostModule *vdtrackingPostModule[2];
  Na61PrimaryVertexRecoPostModule *primaryvtxpostModule;
  
  TFile *fHistFile;
  REGISTER_MODULE("VdTrackFinderEfficiencySG", VdTrackFinderEfficiencySG, "$Id$");
};
}
#endif  // _VdTrackFinderEfficiencySG_VdTrackFinderEfficiencySG_h_
