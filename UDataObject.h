// -*- mode: c++ -*-
//
///////////////////////////////////////////////////////////////////////
//                                                                   //
//    UDataObject                      1244                              //
//                                                                   //
//    NA61 g4 simulation data storage class.                         //
//                                                                   //
//    Author  : Kris Hagel from a template by Gunther Roland         //
//    Modified: Pawel Staszel (for NA61 OpenCharm sim)               //
//    Created : July 13 1997                                         //
//    Version : 1.1                                                  //
//    Changed : 2/24/98 fv                                           //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef UTIL_UDataObject
#define UTIL_UDataObject

#ifndef ROOT_TNamed
#include "TObject.h"
#endif
#ifndef ROOT_TString
#include "TString.h"
#endif

class UDataObject : public TObject {
 public:
  UDataObject();
  UDataObject(const char* name, const char* title = NULL);

  virtual ~UDataObject();
  virtual bool IsPersistent() const { return fIsPersistent; }
  virtual bool IsNode() const { return false; }
  virtual bool IsTable() const { return false; }
  virtual void SetPersistent(bool b) { fIsPersistent = b; }

  virtual const char* GetName() const { return fName; }
  virtual const char* GetTitle() const { return fTitle; }
  virtual void SetName(const Text_t* name);
  virtual void SetTitle(const Text_t* title);
  virtual void Copy(UDataObject& dataobject);
  virtual void Print(Option_t* option = "") const;  //*MENU*

  static const double kInvalidValue;

 protected:
  bool fIsPersistent;  // Flag marking object as persistent
  char fName[64];      // Name of data object;
  char fTitle[64];     // Title of data object;

 public:
  //  ClassDef(UDataObject,1)    // BRAHMS data container class
};

#endif
