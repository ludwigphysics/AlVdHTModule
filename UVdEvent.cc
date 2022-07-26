//________________________________________________________________
//
// UVdEvent is a G4NA61 data class providing storage and
// access function for event data information. The data
// is stored in UDataObject objects that are kept inside
// a THashList.
// This allows modification of the UVdEvent content by the
// user, with the UDataObject providing a standard interface
// for information retrieval, bookkeeping and I/O selection.
//
//________________________________________________________________

//
// $Id: UVdEvent.cxx,v 1.10 2009/09/28 00:26:14 videbaek Exp $
// $Author: videbaek $
// $Date: 2009/09/28 00:26:14 $
//

#include "UVdEvent.h"
//
// Root classes
//
#include "TBuffer.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif
#ifndef __IOMANIP__
#include <iomanip>
#endif
using std::cout;
using std::endl;
using std::flush;
using std::setw;
using std::hex;
using std::dec;
//________________________________________________________________
//ClassImp(UVdEvent);

//________________________________________________________________
UVdEvent::UVdEvent() {
  // Constructor. Set counter and list data members to zero.
  //
  // Two following non persistent variables set for when event header
  // created
  fRunNumber = 0;
  fEventNumber = 0;
  fCorruptedEvent = false;
  fPrimaryVertexStatus = 0;
  fTriggerFlagSet = true;
}

//________________________________________________________________
UVdEvent::UVdEvent(const char* Name, int run, int event) : UEventNode(Name) {
  // Constructor. Create the event header and the hash table
  // for storing the data objects for this event. Set the
  // event name.

  char txt[512];
  sprintf(txt, "event_%d_%d", run, event);

  SetName(txt);
  sprintf(txt, " run %d, event %d", run, event);
  SetTitle(txt);

  // Two following non persistent variables set for when event header
  // created
  fRunNumber = run;
  fEventNumber = event;
  fCorruptedEvent = false;
  fPrimaryVertexStatus = 0;
  fTriggerFlagSet = true;
}

//________________________________________________________________
UVdEvent::~UVdEvent() {
  // Destructor. Delete UVdEvent and all the data objects
  // currently owned by UVdEvent
}

//________________________________________________________________
void UVdEvent::Browse(TBrowser* b) { UEventNode::Browse(b); }

//________________________________________________________________
void UVdEvent::Print(Option_t* option) {
  // Print event information
  // Options:
  //    R       Recursive listing of nodes
  // See also UEventNode::Print and UDataObject::Print
  cout << "*************************************************" << endl << endl << "  Run:      " << GetRunNumber() << "  Event:    " << GetEventNumber() << endl << "*************************************************" << endl;

  TString opt(option);
  opt.ToLower();
  if (opt.Contains("r")) UEventNode::Print(option);
}

//________________________________________________________________
void UVdEvent::Copy(UVdEvent& event) { UEventNode::Copy(event); }
