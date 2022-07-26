//____________________________________________________________________
//
// UEventNode is NA61 simulation data class providing storage and access
// function for event data information. The data is stored in
// TObjArray.
//___________________________________________________________________
//
//
#include "UEventNode.h"
//
// Root classes
//
#include "TBuffer.h"
#include "TClass.h"
#ifndef ROOT_TROOT
#include "TROOT.h"
#endif

#ifndef __IOSTREAM__
#include <iostream>
#endif
using std::cout;
using std::endl;

//____________________________________________________________
// ClassImp(UEventNode);

//____________________________________________________________
UEventNode::UEventNode() {
  // Constructor. Set counter and list data members to zero.
  // Don't use this constructor unless you have to and know
  // what you are doing
  // Use UEventNode(char Name) instead
  fObjectList = 0;
  fVerbose = 0;
}

//____________________________________________________________
UEventNode::UEventNode(const char *Name, const char *Title) : UDataObject(Name, Title) {
  // Constructor. Create the hash table
  // for storing the data objects for this event. Set the
  // eventnode name
  fObjectList = 0;
  fVerbose = 0;
}

//_______________________________________________________________________
// XXX
UEventNode::UEventNode(const UEventNode& /* node */) : UDataObject() {
  std::cerr << "DONT CALL THS METHOD!!!" << std::endl;
  abort();
  //  Copy constructor.  This is done in a simple way for now.
  //  and this will NOT work as defined!!!!
  //
  //int len = sizeof(UEventNode);
  //memcpy(this, &node, len);
}

//_______________________________________________________________________
UEventNode::~UEventNode() {
  // Destructor. Delete UEventNode and all the data objects currently
  // owned by UEventNode.  It is assumed but not checked that all
  // member in the lower nodes are also deleted.
  //

  // It's controversial wether we should actually call
  // TObjArray::Delete here.  Maybe a TObjArray::Clear would be
  // better.
  if (fObjectList) {
    fObjectList->Delete();
    delete fObjectList;
    fObjectList = 0;
  }
}

//_______________________________________________________________________
void UEventNode::CheckList() {
  // Check if list exist. If not, create it.
  if (!fObjectList) fObjectList = new TObjArray();
}

//_______________________________________________________________________
int UEventNode::AddObject(UDataObject *object) {
  // Add a UDataObejct to the EventNode. The object can be simple,
  // or an EventNode..
  //
  CheckList();

  fObjectList->Add(object);
  if (fVerbose) cout << "<UEventNode::AddObject>: Add data Object " << object->GetName() << ", total of " << fObjectList->GetSize() << " objects" << endl;
  return fObjectList->GetSize();
}

//_____________________________________________________________________
UEventNode *UEventNode::GetEventNode(const char *name) const {
  // Return the node named name
  if (!fObjectList) return 0;  // If no list, means nothing in it to get
  TString n(name);

  TIter next(fObjectList);
  UDataObject *object;
  while ((object = (UDataObject *)next())) {
    // First check if the object really is a UEventNode (for type
    // safe returns), and then compare the name.
    if (object->IsA()->InheritsFrom(UEventNode::Class()) && !n.CompareTo(object->GetName())) return (UEventNode *)object;

    // If the next object is a node, then we search that. In this way,
    // we recursively search the tree.
    if (object->IsNode()) {
      object = ((UEventNode *)object)->GetEventNode(name);
      // If the object is found in a sub node, then we return
      // immediatly.
      if (object) return (UEventNode *)object;
    }
  }
  return 0;
}

//_____________________________________________________________________
UDataTable *UEventNode::GetDataTable(const char *name) const {
  // Get a pointer to a UDataObject with the given Name in the node
  // or lower in the tree. Only the first occurence of an object with
  // a specified name is returned. Recursive calls are made if
  // EventNodes are part of the list.
  //
  if (!fObjectList) return 0;  // If no list, means nothing in it to get
  TString n(name);

  TIter next(fObjectList);
  UDataObject *object;
  while ((object = (UDataObject *)next())) {
    // First check if the object really is a UDataTable (for type
    // safe returns), and then compare the name.
    if (object->IsA()->InheritsFrom(UDataTable::Class()) && !n.CompareTo(object->GetName())) return (UDataTable *)object;

    // If the next object is a node, then we search that. In this way,
    // we recursively search the tree.
    if (object->IsNode()) {
      object = ((UEventNode *)object)->GetDataTable(name);
      // If the object is found in a sub node, then we return
      // immediatly.
      if (object) return (UDataTable *)object;
    }
  }
  return 0;
}

//_______________________________________________________________________
UDataObject *UEventNode::GetObject(const char *name) const {
  // Get a pointer to a UDataObject with the given Name in the node
  // or lower in the tree. Only the first occurence of an object with
  // a specified name is returned. Recursive calls are made if
  // EventNodes are part of the list.

  if (!fObjectList) return 0;  // If no list, means nothing in it to get
  TString n(name);

  TIter next(fObjectList);
  UDataObject *object;
  while ((object = (UDataObject *)next())) {
    // Since any object we can add to the node is derived from
    // UDataObject, we really don't need to make a type check here.
    if (!n.CompareTo(object->GetName())) return (UDataObject *)object;

    if (object->IsNode()) {
      object = ((UEventNode *)object)->GetObject(name);
      // If the object is found in a sub node, then we return
      // immediatly.
      if (object) return object;
    }
  }
  return 0;
}

//______________________________________________________
void UEventNode::Clear(Option_t *option) {
  // Clear internal TObjArray objects. That is, remove pointers to
  // actual objects. The objects are only deleted (memory freed) if
  // the TObjArray::SetOwner(true) message has been sent, or the
  // kCanDelete bit is set (see also the class description).
  if (fObjectList) fObjectList->Clear(option);
}

//______________________________________________________
void UEventNode::ListObjects() const {
  // List on standard out a summary of UDataObjects in the
  // EventNode. The listing is recursive. All levels are scanned.
  // Depreciated.  Use Print instead.
  if (!fObjectList) return;  // no need to print if no object list yet

  UDataObject *object;
  int num = 0;
  TIter NextObject(fObjectList);

  while ((object = (UDataObject *)NextObject())) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(2, 25, 3)
    gROOT->IndentLevel();
#else
    IndentLevel();
#endif
    cout << num << " - " << object << " - " << object->GetName() << " - " << object->GetTitle() << endl;

    if (object->IsNode()) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(2, 25, 3)
      gROOT->IncreaseDirLevel();
#else
      IncreaseDirLevel();
#endif
      ((UEventNode *)object)->ListObjects();
#if ROOT_VERSION_CODE >= ROOT_VERSION(2, 25, 3)
      gROOT->DecreaseDirLevel();
#else
      DecreaseDirLevel();
#endif
    }
    num++;
  }
}

//______________________________________________________
void UEventNode::Copy(UEventNode &eventnode) {
  // Copy contents from this into eventnode

  UDataObject::Copy(eventnode);
  eventnode.fVerbose = fVerbose;

  if (!fObjectList) return;

  TIter next(fObjectList);
  UDataObject *dataObj;
  while ((dataObj = (UDataObject *)next())) eventnode.AddObject(dataObj);
}

//________________________________________________________________
UDataObject *UEventNode::RemoveDataObject(UDataObject *objectToBeRemoved) {
  // Remove Dataobject from Table, but do not delete obejct.
  // This must scan down the tree to find the object. Otherwise
  // it will not be deleted from the tree.
  //
  if (!fObjectList) return 0;

  TIter next(fObjectList);
  UDataObject *object;
  while ((object = (UDataObject *)next())) {
    // Since any object we can add to the node is derived from
    // UDataObject, we really don't need to make a type check here.
    if (objectToBeRemoved == object) {
      fObjectList->Remove(object);
      return object;
    }
    if (object->IsNode()) {
      object = ((UEventNode *)object)->RemoveDataObject(objectToBeRemoved);
      if (object) return object;
    }
  }

  return 0;
}

//________________________________________________________________
UEventNode &UEventNode::operator+=(const UEventNode &node) {
  // This operator takes matching elements of nodes and adds them
  // together. The cases when one node has elements that the other
  // does not is not yet taken into account.

  // objlist1 is local object list
  // objlist2 is object list of node we are adding to.
  // nument1 is number of entries of local object list
  // nument2 is number of entries of node we are adding to
  // object1 is object found in local list
  // object2 is object found in node we are adding to

  TObjArray *objlist1;
  TObjArray *objlist2;
  objlist1 = GetObjectList();
  if (!objlist1) return *this;

  int nument1 = objlist1->GetEntries();

  objlist2 = node.GetObjectList();
  if (!objlist2) return *this;

  int nument2 = objlist2->GetEntries();
  TObject *object1 = 0;
  TObject *object2 = 0;

  for (int i1 = 0; i1 < nument1; i1++) {
    object1 = objlist1->At(i1);

    for (int i2 = 0; i2 < nument2; i2++) {
      object2 = objlist2->At(i2);

      if (!strcmp(object1->GetName(), object2->GetName())) {
        // We have found two objects with the same name.  Add them
        // together
        if (object1->IsA() == UEventNode::Class())
          *(UEventNode *)object1 += *(UEventNode *)object2;
        else if (object1->IsA() == UDataTable::Class())
          //*(UDataTable*)object1 += *(UDataTable*)object2;
          continue;
#if 0 
	// I took out this instance,  since  it seems rediculus to
	// have a special case for it.
	else if(object1->IsA() == UGeantHeader::Class()) 
	  // It is not clear to me what to do about UGeantHeaders
	  // when mixing events. It should be discussed if this
	  // facility is really used.  For the moment, I just add
	  // them to the list.  How we get them both out will be
	  // another story. 
	  this->AddObject((UDataObject*)object2);
#endif
        else {
#ifdef USE_NODE_FULL_CONCAT
          this->AddObject((UDataObject *)object2);
#else
          cout << "Object found in EventNode " << node.GetName() << " is off class " << object1->IsA()->Class_Name() << "." << endl
               << "There is no support for it yet, "
               << "please see your BRAT manager" << endl;
#endif
        }
        // Ueak out of this loop, cause we found machting objects.
        // Zero the thing, so that we know if the object was found of
        // not.
        object2 = 0;
        break;
      }
    }
  }

  return *this;
}

//________________________________________________________________
void UEventNode::Browse(TBrowser *b) {
  if (!fObjectList) return;

  TIter next(fObjectList);
  TObject *obj;
  while ((obj = next())) b->Add(obj);
}

//________________________________________________________________
void UEventNode::Print(Option_t *option) const {
  // Print all contained objects, passing the option along.
  // Options:
  //    R         Recursive print [Default]
  // See also UDataObject::Print
  UDataObject::Print(option);

  if (!fObjectList) return;  // no need to go further if no object list.

  TString opt(option);
  opt.ToLower();
  if (opt.Contains("r")) {
    gROOT->IncreaseDirLevel();

    TIter next(fObjectList);
    TObject *obj = 0;
    while ((obj = next())) {
      gROOT->IndentLevel();
      obj->Print(option);
    }

    gROOT->DecreaseDirLevel();
  }
}

//________________________________________________________________
UEventNode &operator+(const UEventNode &node1, const UEventNode &node2) {
  // This operator takes matching elements of nodes and adds them
  // together. The cases when one node has elements that the other
  // does not is not yet taken into account.
  TObjArray *objlist1, *objlist2;
  objlist1 = node1.GetObjectList();
  int nument1 = objlist1->GetEntries();
  objlist2 = node2.GetObjectList();
  int nument2 = objlist2->GetEntries();
  UEventNode *tmp = new UEventNode(node1.GetName());
  if (objlist1) {
    for (int i1 = 0; i1 < nument1; i1++) {
      TObject *object1 = objlist1->At(i1);
      for (int i2 = 0; i2 < nument2; i2++) {
        TObject *object2 = objlist2->At(i2);
        if (!strcmp(object1->GetName(), object2->GetName())) {
          //  We have found two objects with the same name.  Add them
          //  to the  output event
          if (!strcmp(object1->IsA()->GetName(), "UEventNode")) {
            UEventNode *tmpnode = new UEventNode(*(UEventNode *)object1 + *(UEventNode *)object2);
            tmp->AddEventNode(tmpnode);
          } else if (!strcmp(object1->IsA()->GetName(), "UDataTable")) {
            UDataTable *tmptable = new UDataTable("NULL");
            // new UDataTable(*(UDataTable*)object1
            //		       + *(UDataTable*)object2);
            tmp->AddDataTable(tmptable);
          } else {
            cout << "Object found in EventNode " << node1.GetName() << " is " << object1->IsA()->GetName() << "." << endl;
            cout << "There is no support for it yet, "
                    "please see your BRAT manager"
                 << endl;
          }
        }
      }
    }
  }
  return *tmp;
}

//____________________________________________________________________
// XXX
UEventNode &UEventNode::operator=(const UEventNode& /* rhs */) {
  std::cerr << "DONT CALL THIS METHOD!!!!" << std::endl;
  abort();
  // UEventNode assignment operator.
  // This is done with a simple memcpy for the moment.

  //if (this != &rhs) {
    //int len = sizeof(UEventNode);
    //memcpy(this, &rhs, len);
  //}
  //return *this;
}

//____________________________________________________________
UEventNode &UEventNode::operator=(const UEvent & /*rhs*/) {
  // UEvent to UEventNode assignment operator.
  cout << "Inside assignment of UEvent to UEventNode. "
       << "Should we be here? see the manager" << endl;
  return *this;
}
