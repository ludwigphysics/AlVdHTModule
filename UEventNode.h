// -*- mode: c++ -*-
// $Id: UEventNode.h,v 1.5 2004/01/31 19:57:53 videbaek Exp $
//
///////////////////////////////////////////////////////////////////////
//                                                                   //
//    UEventNode                                                    //
//                                                                   //
//    Uahms event class                                             //
//                                                                   //
//    UEventNode manages access to the raw and reconstructed data   //
//    for one event branch in the BRAT environment.                  //
//    environment                                                    //
//                                                                   //
//    Author  : F.Videbaek                                           //
//    Created :                                                      //
//    Version : 1.0                                                  //
//    Changed : 2/24/1998                                            //
//                                                                   //
//                                                                   //
//                                                                   //
///////////////////////////////////////////////////////////////////////
#ifndef UTIL_UEventNode
#define UTIL_UEventNode
#include "TNamed.h"
#include "TObjArray.h"
#include "UDataObject.h"
#include "UDataTable.h"

#define UEvent UEvent_ // hack due to clash with UniGen

class TBrowser;
class UEvent;

class UEventNode : public UDataObject {
 public:
  UEventNode();
  UEventNode(const char* Name, const char* Title = NULL);
  UEventNode(const UEventNode& node);  // copy constructor

  virtual ~UEventNode();

  // Enable and Disable Verbose mode. Diagnostice out put occurs
  // at each call of the Event memebr functions.
  virtual int AddObject(UDataObject* Object);
  virtual int AddDataTable(UDataTable* table);
  virtual int AddEventNode(UEventNode* node);
  virtual void Browse(TBrowser* b);
  virtual void Clear(Option_t* o = "");
  virtual void Copy(UEventNode& eventnode);
  virtual UDataTable* GetDataTable(const char* ObjectName) const;
  virtual UEventNode* GetEventNode(const char* NodeName) const;
  virtual UDataObject* GetObject(const char* ObjectName) const;
  virtual TObjArray* GetObjectList() const { return fObjectList; }
  virtual bool IsNode() const { return true; }
  virtual bool IsFolder(void) const { return true; }
  virtual void ListObjects() const;  // *MENU*
  virtual UDataObject* RemoveDataObject(UDataObject* object);
  virtual void SetOwner(bool e = true);
  virtual void SetVerbose(int i) { fVerbose = i; }
  virtual void Print(Option_t* option = "R") const;  //*MENU*

  // UEventNode  operator + () const;
  UEventNode& operator+=(const UEventNode&);
  UEventNode& operator=(const UEventNode& rhs);  // idem
  UEventNode& operator=(const UEvent& rhs);      // idem

 private:
  void CheckList();        // Check that list exists.
  TObjArray* fObjectList;  //  List of data objects
  int fVerbose;            //! Controls debugging and monitoring output

  friend UEventNode& operator+(const UEventNode& node1, const UEventNode& node2);

 public:
  //  ClassDef(UEventNode,1)       // BRAHMS event data class
};

//____________________________________________________________________
inline int UEventNode::AddDataTable(UDataTable* table) { return AddObject(table); }

//____________________________________________________________________
inline int UEventNode::AddEventNode(UEventNode* node) { return AddObject(node); }

//____________________________________________________________________
inline void UEventNode::SetOwner(bool e) {
  CheckList();
  fObjectList->SetOwner(e);
}

#endif
