// -*- mode: c++ -*-
//

//
#ifndef UTIL_ULine3D
#define UTIL_ULine3D

#ifndef UTIL_UDataObject
#include "UDataObject.h"
#endif

#ifndef UTIL_UVector3D
#include "UVector3D.h"
#endif

#ifndef WIN32
#include <iostream>
#else
#include <iostream.h>
#endif
#include <cassert>

class Line3D : public UDataObject {
 private:
  Vector3D fOrigin;     //  Vector for line start
  Vector3D fDirection;  //  Unit direction vector

 public:
  Line3D();
  Line3D(const Vector3D& origin, const Vector3D& direction);
  Line3D(const double ox, const double oy, const double oz, const double dx, const double dy, const double dz);

  Line3D(const Line3D& line);
  virtual ~Line3D();

  Line3D& operator=(const Line3D&);

  virtual Vector3D ClosestPoint(const Vector3D&) const;
  virtual double Distance(const Vector3D&) const;
  virtual const Vector3D& GetOrigin() const { return fOrigin; };
  virtual const Vector3D& GetDirection() const { return fDirection; };
  virtual double GetShortestDistanceBetween(const Line3D& line) const;
  virtual double GetShortestDistanceBetween(const Line3D* line) const;

  virtual Line3D GetShortestLineBetween(const Line3D& line) const;

  virtual Vector3D GetClosestProximityPoint(const Line3D& line) const;
  virtual Vector3D GetClosestProximityPoint(const Line3D* line) const;

  double GetXatZ(double z) {
    double xatz = fOrigin.X() + (z - fOrigin.Z()) * fDirection.X() / fDirection.Z();
    return xatz;
  }

  double GetYatZ(double z) {
    double yatz = fOrigin.Y() + (z - fOrigin.Z()) * fDirection.Y() / fDirection.Z();
    return yatz;
  }

  virtual void Print(Option_t* option = "") const;
  virtual double RelativeOverlap(const Line3D& line, double zl, double r);

  virtual void SetOrigin(const Vector3D& vector);
  virtual void SetOrigin(double x, double y, double z) { SetOrigin(Vector3D(x, y, z)); }
  virtual void SetDirection(const Vector3D& vector);
  virtual void SetDirection(double x, double y, double z) { SetDirection(Vector3D(x, y, z)); }

  friend ostream& operator<<(ostream& str, const Line3D*);
  friend ostream& operator<<(ostream& str, const Line3D&);

  //  ClassDef(Line3D,0)  //  DC raw data class
};

extern ostream& operator<<(ostream& str, const Line3D*);
extern ostream& operator<<(ostream& str, const Line3D&);

#endif

//____________________________________________________________________
//
// $Log: $
