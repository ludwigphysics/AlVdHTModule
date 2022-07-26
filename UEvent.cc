//________________________________________________________________
//
// UEvent is a G4NA61 data class providing storage and
// access function for event data information. The data
// is stored in UDataObject objects that are kept inside
// a THashList.
// This allows modification of the UEvent content by the
// user, with the UDataObject providing a standard interface
// for information retrieval, bookkeeping and I/O selection.
//
//________________________________________________________________

//
// $Id: UEvent.cxx,v 1.10 2009/09/28 00:26:14 videbaek Exp $
// $Author: videbaek $
// $Date: 2009/09/28 00:26:14 $
//

#include "UEvent.h"
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
// ClassImp(UEvent);

//________________________________________________________________
UEvent::UEvent() {
  // Constructor. Set counter and list data members to zero.
  //
  // Two following non persistent variables set for when event header
  // created
  fRunNumber = 0;
  fEventNumber = 0;
}

//________________________________________________________________
UEvent::UEvent(const char* Name, int run, int event) : UEventNode(Name) {
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
  fPrimaryRecoVertexSet = false;
}

//________________________________________________________________
UEvent::~UEvent() {
  // Destructor. Delete UEvent and all the data objects
  // currently owned by UEvent
}

//________________________________________________________________
void UEvent::Browse(TBrowser* b) { UEventNode::Browse(b); }

//________________________________________________________________
void UEvent::Print(Option_t* option) {
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
void UEvent::Copy(UEvent& event) { UEventNode::Copy(event); }
