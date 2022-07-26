#ifndef _AlVdHitConverterSG_AlVdHitConverterSG_h_
#define _AlVdHitConverterSG_AlVdHitConverterSG_h_

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
#include "TH2F.h"

#include "Na61VdParameters.h"
#include "Na61AlVdArmParameters.h"
#include "Na61VdParametersManager.h"


namespace AlVdHitConverterSG {

/**
\class AlVdHitConverterSG

\brief Describe your module. In one sentence.

Now here a longer description in doxygen. You may use HTML.

\author Anastasia Merzlaya
\date 27 Jun 2020
\ingroup TODO: Put a group here.
*/

class AlVdHitConverterSG : public boost::noncopyable, public fwk::VModule {
 public:
  fwk::VModule::EResultFlag Init();
  fwk::VModule::EResultFlag Process(evt::Event &event, const utl::AttributeMap &attr);
  fwk::VModule::EResultFlag Finish();

 private:
  unsigned int prevRun;
  void Begin(const unsigned int run);
  void End();
  void LocalToArmGlobal(unsigned char arm, unsigned char sensor, utl::Point &localPos,utl::Point &armstandalonePos);
  void LocalToVdGlobal(unsigned char arm, utl::Point &armstandalonePos, utl::Point &vdstandalonePos);
  void LocalToShineGlobal(unsigned char arm, utl::Point &vdstandalonePos, utl::Point &shinePos);
  
  TFile *fHistFile;

  TH2F* fhZX_fine;
  TH1F* fhX_Al1;
  TH1F* fhX_Al2;
  TH1F* fhX_Al3;
  TH1F* fhX_Al4;
  TH1F* fhY_Al1_j;
  TH1F* fhY_Al2_j;
  TH1F* fhY_Al3_j;
  TH1F* fhY_Al4_j;
  TH1F* fhY_Al1_s;
  TH1F* fhY_Al2_s;
  TH1F* fhY_Al3_s;
  TH1F* fhY_Al4_s;

  TH2F* fhZX_fine_G4;
  TH1F* fhX_Al1_G4;
  TH1F* fhX_Al2_G4;
  TH1F* fhX_Al3_G4;
  TH1F* fhX_Al4_G4;
  TH1F* fhY_Al1_j_G4;
  TH1F* fhY_Al2_j_G4;
  TH1F* fhY_Al3_j_G4;
  TH1F* fhY_Al4_j_G4;
  TH1F* fhY_Al1_s_G4;
  TH1F* fhY_Al2_s_G4;
  TH1F* fhY_Al3_s_G4;
  TH1F* fhY_Al4_s_G4;

  // This goes at the end.
  REGISTER_MODULE("AlVdHitConverterSG", AlVdHitConverterSG, "$Id$");
};
}
#endif  // _AlVdHitConverterSG_AlVdHitConverterSG_h_
