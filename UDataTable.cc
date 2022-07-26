//________________________________________________________________
//
// UDataTable is a data object that has a list of
// data objects.  It is meant to be used in grouping like kinds
// of information together.  The list is fObjectList which is a
// TObjArray.  Several methods manage the fObjectList with the same
// names as TObjArray
//
//________________________________________________________________

#include "UDataTable.h"

#ifndef ROOT_TClass
#include "TClass.h"
#endif
//#ifndef ROOT_TBufferFile
//#include "TBufferFile.h"
//#endif
#ifndef ROOT_TROOT
#include "TROOT.h"
#endif

#if !defined __IOSTREAM__
#include <iostream>
#endif
using std::cout;
using std::endl;

// ClassImp(UDataTable);

//________________________________________________________________
UDataTable::UDataTable() {
  // Defaults constructor. Does nothing.
  // Don't use this constructor unless you have to and know
  // what you are doing
  // Use UDataTable(char *name) instead.

  fObjectList = 0;
}

//________________________________________________________________
UDataTable::UDataTable(const char *Name, const char *Title) : UDataObject(Name, Title) {
  // Constructor. Create the data container setting the name
  // (and title, if supplied)

  if (Title) SetTitle(Title);

  fObjectList = new TObjArray();
}

//________________________________________________________________
UDataTable::~UDataTable() {
  // Destructor. Delete UDataTable and all the data objects
  // currently owned by UDataTable
  if (fObjectList) {
    fObjectList->Delete();
    delete fObjectList;
    fObjectList = 0;
  }
}

//________________________________________________________________
void UDataTable::Add(TObject *object) {
  // Add an object to the object list.  Essentially uses
  // TObjArray::Add(object);
  if (!fObjectList) fObjectList = new TObjArray();

  fObjectList->Add(object);
}

//________________________________________________________________
void UDataTable::AddAt(TObject *object, int idx) {
  // Add an object to the object list at a specific index.
  // Essentially  uses TObjArray::AddAt(object,idx);
  if (!fObjectList) fObjectList = new TObjArray();

  fObjectList->AddAt(object, idx);
}

//________________________________________________________________
void UDataTable::DeleteAndCompress(TObject *obj) {
  // Removes an object from the Object List.  Object is then deleted.
  // Then the object list is compressed
  DeleteObject(obj);
  Compress();
}

//________________________________________________________________
void UDataTable::DeleteAndCompressAt(int i) {
  // Removes an object an object at specified index from the Object List.
  // Essentially uses TObjArray::RemoveAt(i);
  // After object is removed, it is deleted.
  // Then the object list is compressed
  DeleteObjectAt(i);
  Compress();
}

//________________________________________________________________
void UDataTable::DeleteObject(TObject *obj) {
  // Removes an object from the Object List.  Object is then deleted.
  delete (fObjectList->Remove(obj));
}

//________________________________________________________________
void UDataTable::Browse(TBrowser *b) {
  TIter next(fObjectList);
  TObject *obj;
  while ((obj = next())) b->Add(obj, obj->GetName());
}

//________________________________________________________________
void UDataTable::Print(Option_t *option) const {
  // Print all contained objects, passing the option along.
  // Options:
  //    R         Recursive print [Default]
  // See also UDataObject::Print
  UDataObject::Print(option);

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

/*
//________________________________________________________________
UDataTable&
UDataTable::operator += (const UDataTable& table)
{
  // Add all of the entries of the argument to the table we are in.
  // Copies of the entries in the adding table are made.  The copies
  // are made by transferring the persistent data into a buffer using
  // the Streamer() method, then transferring from the buffer into the
  // new object using the Streamer() method of the new object.  One
  // should evaluate if we really want copies or just to  point to the
  // objects.  That would be much easier to do, but can cause
  // conflicts when the entries in table1 and table2 are deleted and
  // one still tries to access the entries in the new table.

  TBufferFile buf(TBuffer::kRead);

  TObject *object;
  TClass  *cl;

  int nument = table.GetEntries();

  //First do for the first table
  for(int i1 = 0; i1 < nument; i1++) {
    object = table.At(i1);
    cl = object->IsA();
    // Create a new instance of this class
    TObject *p = (TObject*)cl->New();
    // Now we need to get the data out of the other one and put into
    // the new instance. We do that here by reading the data from the
    // original object into a buffer via the Streamer.  Then the
    // Buffer Offset is reset to beginning and buffer is set to read
    // mode and the streamer for the new object is run to get out of
    // the buffer and into the new object.  A simple mcpy might be
    // faster, and would get all of the data instead of just the
    // persistent data, but might have some nasty surprises.
    buf.Reset();
    buf.SetWriteMode();
    object->Streamer(buf);
    buf.Reset();
    buf.SetReadMode();
    p->Streamer(buf);

    // Add this element to the table.
   this->Add(p);
  }

  return *this;
}
*/

/*
//________________________________________________________________
UDataTable
operator+(const UDataTable& table1,
          const UDataTable& table2)
{
  // Create a new data table and put all of the entries of the first
  // one and the second one into the new one.  Copies of the entries
  // are made.  The copies are made by transferring the persistent
  // data into a buffer using the Streamer() method,  then
  // transferring from the buffer into the new object using the
  // Streamer() method of the new object.  One should evaluate if we
  // really want copies or just to  point to the objects.  That would
  // be much easier to do, but can cause conflicts when the entries in
  // table1 and table2 are deleted and one still tries to access the
  // entries in the new table.

  TBufferFile buf(TBuffer::kRead);

  TObject *object;
  TClass  *cl;
  UDataTable *tmp = new UDataTable(table1.GetName());

  int nument1 = table1.GetEntries();
  int nument2 = table2.GetEntries();

  //First do for the first table
  for(int i1=0;i1<nument1;i1++) {
    object = table1.At(i1);
    cl = object->IsA();
    // Create a new instance of this class
    TObject *p = (TObject*)cl->New();
    // Now we need to get the data out of the other one and put into
    // the new instance.  We do that here by reading the data from the
    // original object into a buffer via the Streamer.  Then the
    // Buffer Offset is reset to beginning and buffer is set to read
    // mode and the streamer for the new object is run to get out of
    // the buffer and into the new object.  A simple mcpy might be
    // faster, and would get all of the data instead of just the
    // persistent data, but might have some nasty surprises.
    buf.Reset();
    buf.SetWriteMode();
    object->Streamer(buf);
    buf.Reset();
    buf.SetReadMode();
    p->Streamer(buf);

    // Add this element to the table.
    tmp->Add(p);
  }

  //Now do for the second table
  for(int i2=0;i2<nument2;i2++) {
    object = table2.At(i2);
    cl = object->IsA();
    // Create a new instance of this class
    TObject *p = (TObject*)cl->New();
    // Now we need to get the data out of the other one and put into
    // the new instance. We do that here by reading the data from the
    // original object into a buffer via the streamer.  Then the
    // Buffer Offset is reset to beginning and buffer is set to read
    // mode and the streamer for the new object is run to get out of
    // the buffer and into the new object.  A simple mcpy might be
    // faster and would get all of the data instead of just the
    // persistent data, but might have some nasty surprises.
    buf.Reset();
    buf.SetWriteMode();
    object->Streamer(buf);
    buf.Reset();
    buf.SetReadMode();
    p->Streamer(buf);

    // Add this element to the table.
    tmp->Add(p);
  }

  return *tmp;
}
*/
