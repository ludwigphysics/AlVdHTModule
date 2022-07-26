// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////
//                                                                   //
//    UDataTable                                                    //
//                                                                   //
//    NA61 geant4 data container class. Raw data and the results of  //
//    calibration and reconstruction operations are stored in        //
//    BrDataTable objects, which in turn are managed and             //
//    stored using the BrEvent class.                                //
//                                                                   //
//    Author  : Kris Hagel from a template by Gunther Roland         //
//    Modified for NA61 by Pawel Staszel                             //
//    Created : July 13 1997                                         //
//    Version : 1.0                                                  //
//    Changed :                                                      //
//                                                                   //
///////////////////////////////////////////////////////////////////////

//
#ifndef UTIL_UDataTable
#define UTIL_UDataTable

// Root Classes
#include "TObjArray.h"
#include "TBrowser.h"

// BRAHMS Classes

#include "UDataObject.h"

class UDataTable : public UDataObject {
 public:
  UDataTable();
  UDataTable(const char* Name, const char* Title = NULL);

  virtual ~UDataTable();

  virtual void Add(TObject* object);
  virtual void AddAt(TObject* object, int idx);

  virtual TObjArray* GetObjectList() { return fObjectList; }
  int Entries() {
    if (fObjectList)
      return fObjectList->GetEntries();
    else
      return 0;
  }
  int GetEntries() const {
    if (fObjectList)
      return fObjectList->GetLast() + 1;
    else
      return 0;
  }
  UDataTable* At(int i) const { return (UDataTable*)fObjectList->At(i); }
  virtual bool IsTable() const { return true; }
  void Sort(int upto) { fObjectList->Sort(upto); }
  void Sort() {
    if (fObjectList) fObjectList->Sort(Entries());
  }
  void Clear() {
    if (fObjectList) fObjectList->Clear();
  }
  void Delete() {
    if (fObjectList) fObjectList->Delete();
  }
  void Remove(TObject* obj) {
    if (fObjectList) fObjectList->Remove(obj);
  }
  void RemoveAt(int i) {
    if (fObjectList) fObjectList->RemoveAt(i);
  }
  void Compress() {
    if (fObjectList) fObjectList->Compress();
  }
  void DeleteAndCompress(TObject* obj);
  void DeleteAndCompressAt(int i);
  void DeleteObject(TObject* obj);
  inline void DeleteObjectAt(int i) {
    // Removes an object an object at specified index from the Object List.
    // Essentially uses TObjArray::RemoveAt(i);
    // After object is removed, it is deleted.
    delete (fObjectList->RemoveAt(i));
  }

  virtual void SetOwner(bool e = true) { fObjectList->SetOwner(e); }

  virtual bool IsFolder(void) const { return true; }
  virtual void Browse(TBrowser* b);
  virtual void Print(Option_t* option = "R") const;  //*MENU*

  // UDataTable  operator + () const;
  // UDataTable&       operator+=(const UDataTable&);
 protected:
  TObjArray* fObjectList;  // List of objects stored in container
 private:
  // friend UDataTable operator+(const UDataTable& table1,
  //			      const UDataTable& table2);
 public:
  //  ClassDef(UDataTable, 1)      // BRAHMS data container class
};

#endif
