//////////////////////////////////////////////////////
//
// So far cuts were tested on run range 33201 - 33375
// The event cuts were discussed in A.Merzlaya presentation:
// https://indico.cern.ch/event/810328/contributions/3437667/attachments/1851219/3039108/MerzlayaA_XeLa150_28.05.19.pdf 
//
//////////////////////////////////////////////////////

#include <evt/Event.h>
#include <evt/RecEvent.h>
#include <evt/rec/RecEventConst.h>
#include <utl/ShineUnits.h>
#include <algorithm>
#include "TMath.h"



namespace tmpCmp{
  bool absComp(double i, double j) {return (TMath::Abs(i) < TMath::Abs(j));}
}


bool XeLa150EventCuts(const evt::Event& event, det::TriggerConst::EId triggerID){


  if(event.GetRawEvent().GetBeam().GetTrigger().IsTrigger(triggerID, det::TriggerConst::ePrescaled) == false)
    return false;

  const evt::RecEvent& recEvent = event.GetRecEvent();
  const evt::raw::Trigger& trigger = event.GetRawEvent().GetBeam().GetTrigger();


  // WFA S1
  std::vector<double> timeStructureWFAS1;
  if (!trigger.HasTimeStructure(det::TimeStructureConst::eWFA, det::TriggerConst::eS1_1)){
    return false;
  }
  timeStructureWFAS1 = trigger.GetTimeStructure(det::TimeStructureConst::eWFA, det::TriggerConst::eS1_1);
  std::sort(timeStructureWFAS1.begin(), timeStructureWFAS1.end(), tmpCmp::absComp);
  double wfaS1Max = (timeStructureWFAS1.size() > 1) ? timeStructureWFAS1[1] : 25500.;
  if(wfaS1Max>-4000 && wfaS1Max<4000)
    return false;

  // WFA T4
  std::vector<double> timeStructureWFAT4;
  if (!trigger.HasTimeStructure(det::TimeStructureConst::eWFA, det::TriggerConst::eT4)){
    return false;
  }
  timeStructureWFAT4 = trigger.GetTimeStructure(det::TimeStructureConst::eWFA, det::TriggerConst::eT4);
  std::sort(timeStructureWFAT4.begin(), timeStructureWFAT4.end(), tmpCmp::absComp);
  double wfaT4Max = (timeStructureWFAT4.size() > 1) ? timeStructureWFAT4[1] : 25500.;
  if(wfaT4Max>-25000 && wfaT4Max<25000)
    return false;

  // S1 charge cut - no cut implemented

  // S2 charge cut - no cut implemented

  // BPD quality cut 
  if (recEvent.GetBeam().GetStatus() & (evt::rec::BeamConst::eNotFitted | evt::rec::BeamConst::eBadBPD3))
    return false;

  // BPD position cut - no cut implemented

  // BPD3 charge cut
  const int BPDsigX = 310;
  const int BPDsigY = 267;
  const int BPDoffX = 5000;
  const int BPDoffY = 4636;
  double BPD3chargeX = recEvent.GetBeam().GetBPDPlane(det::BPDConst::eBPD3,det::BPDConst::eX).GetCharge();
  double BPD3chargeY = recEvent.GetBeam().GetBPDPlane(det::BPDConst::eBPD3,det::BPDConst::eY).GetCharge();
  if (!((BPD3chargeX>(BPDoffX-4*BPDsigX))&&(BPD3chargeX<(BPDoffX+4*BPDsigX)) && (BPD3chargeY>(BPDoffY-4*BPDsigY))&&(BPD3chargeY<(BPDoffY+4*BPDsigY))))
    return false;


  //--------selecting T1 events---------
  if(triggerID == det::TriggerConst::eT1) 
    return true;


  // S3
  //const int S3Max = 70;
  //int S3 = trigger.GetADC(det::TriggerConst::eS3);    
  //   if (S3 > S3Max)
  //		return false;

  // TPC vertex cut (Note: Standart TPC vertex)
  //   if (!recEvent.HasPrimaryVertex(evt::rec::VertexConst::ePrimaryFitZ))
  //   	return false;
  //	const evt::rec::Vertex&  TPCVertex=recEvent.GetPrimaryVertex(evt::rec::VertexConst::ePrimaryFitZ);
  //	if (!(recEvent.HasPrimaryVertex(evt::rec::VertexConst::ePrimaryFitZ))&&(TPCVertex.GetFitQuality() == evt::rec::FitQuality::ePerfect))
  //		return false;

  // TPC vertex Zposition cut (Note: Standart TPC vertex)
  //const double TPCVertexZposMax = -601.;
  //const double TPCVertexZposMin = -604.5;
  //  if (!((TPCVertex.GetPosition().GetZ()<TPCVertexZposMax)&&(TPCVertex.GetPosition().GetZ()>TPCVertexZposMin)))
  //		return false;

  // VD vertex cut (Note: Don't use it if VD vertex is not included into shoe files!)
  //	if (!(recEvent.HasPrimaryVertex(evt::rec::VertexConst::ePrimaryVD)))
  //		return false;

  // VD vertex Zposition cut - no cut needed


  //--------selecting T4 events---------
  if(triggerID == det::TriggerConst::eT4) 
    return true;


  //PSD
  //double PSDenergy = recEvent.GetPSD().GetEnergy();
  if (recEvent.GetPSD().GetStatus()!=0)
    return false;

  double PSDenergy = recEvent.GetPSD().GetEnergy();
  if (PSDenergy>11450)
    return false;

  //--------selecting T2 events---------
  if(triggerID == det::TriggerConst::eT2) 
    return true;


  return false; //otherwise

}
