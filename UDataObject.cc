//____________________________________________________________________
//
// UDataObject is the base class for all  classes containing
// data, either as raw data from the DAQ or calibrated data from the
// calibration procedures or the results from the reconstruction or
// analysis processors or geant4 simulation.
//
// The main purpose of the UDataObject class is to provide a uniform
// interface for bookkeeping, information exchange between procedural
// classes and persistent I/O
//
// The individual data objects referring to one event are kept with
// the BrEvent object, where they are stored in a hash table using a
// unique name for retrieval.
//
// This class defines the method Print, which will print the basic
// information on a data obejct. Derived class should overload this
// method as:
//
//   void <derived class>::Print(Option_t* option) const {
//     // Print information on this instance of <derived class>
//     // Options:
//     //    <list options here>
//     // See also UDataObject for additional options
//     UDataObject::Print(option);
//     <print the relevant information>
//   }
//
// The overloaded method should not print the name or anything like
// that. That is taken care of in UDataObject::Print.
//
//____________________________________________________________________

//____________________________________________________________________
#include "UDataObject.h"

#ifndef ROOT_TClass
#include "TClass.h"
#endif

#ifndef __IOSTREAM__
#include <iostream>
#endif
#ifndef __IOMANIP__
#include <iomanip>
#endif
using std::cout;
using std::endl;
using std::cerr;
using std::setw;

//____________________________________________________________________
// ClassImp(UDataObject);

//____________________________________________________________________
const double UDataObject::kInvalidValue = 1E200;

//____________________________________________________________________
UDataObject::UDataObject() {
  // Default constructor. Does nothing.
  // Don't use this constructor unless you have to and know
  // what you are doing
  // Use UDataObject(char *name) instead.
  fIsPersistent = true;
}

//____________________________________________________________________
UDataObject::UDataObject(const char *name, const char *title) {
  // Constructor. Create the data container setting the name
  // (and title, if supplied)
  //

  SetName(name);
  SetTitle((title ? title : name));

  fIsPersistent = true;
}

//____________________________________________________________________
UDataObject::~UDataObject() {
  // Destructor. Delete UDataObject and all the data objects
  // currently owned by UDataObject
}

//____________________________________________________________________
void UDataObject::SetName(const char *name) {
  // Set Name of object.
  // This object does not inherit from TNamed as originally because
  // ROOT crashes when a UDataObject is used in the top level of a
  // Tree with split=1.
  // Problem was traced back to TString not writing or reading (or
  // both) correctly.  We therefore inherit from TObject and set our
  // own name as a character variable.
  if (strlen(name) > 64)
    strncpy(fName, name, 64);
  else
    strcpy(fName, name);
}

//____________________________________________________________________
void UDataObject::SetTitle(const Text_t *title) {
  // Set Title of object.
  // This object does not inherit from TNamed as originally because
  // ROOT crashes when a UDataObject is used in the top level of a
  // Tree with split=1.
  // Problem was traced back to TString not writing or reading (or
  // both) correctly.  We therefore inherit from TObject and set our
  // own title as a character variable.
  if (strlen(title) > 64)
    strncpy(fTitle, title, 64);
  else
    strcpy(fTitle, title);
}

//____________________________________________________________________
void UDataObject::Copy(UDataObject &dataobject) {
  // Copy method.
  // Copy All elements of object from argument; also involk copy
  // method of object we inherit from.
  TObject::Copy(dataobject);
  dataobject.fIsPersistent = fIsPersistent;
  dataobject.SetName(fName);
  dataobject.SetTitle(fTitle);
}

//____________________________________________________________________
void UDataObject::Print(Option_t *option) const {
  // Print information on this data object.
  // Options:
  //    D         Details
  //    C         Creator information [not used yet]
  TString opt(option);
  opt.ToLower();
  if (opt.Contains("d"))
    cout << "[" << this << "]" << GetName() << "(Class: " << IsA()->GetName() << ") " << endl << "   " << GetTitle() << endl << "Is " << (fIsPersistent ? "" : "n't") << " persistent" << endl;
  else
    cout << GetName() << " " << GetTitle() << endl;
}
