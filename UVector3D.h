// -*- mode: c++ -*-

#ifndef UTIL_Vector3D
#define UTIL_Vector3D

#ifndef UTIL_UDataObject
#include "UDataObject.h"
#endif

#ifndef ROOT_TVector3
#include "TVector3.h"
#endif

#ifndef WIN32
#include <iostream>
#else
#include <iostream.h>
#endif

//
//

class Vector3D : public UDataObject {
 public:
  // constructors and destructors
  //
  Vector3D(const double x = 0.0, const double y = 0.0, const double z = 0.0) : fX(x), fY(y), fZ(z){};

  //  Vector3D(const double *pos) {fX=pos[0];fY=pos[1];fZ=pos[2];};
  Vector3D(const double* pos) : fX(pos[0]), fY(pos[1]), fZ(pos[2]){};
  Vector3D(TVector3 vect) {
    fX = vect[0];
    fY = vect[1];
    fZ = vect[2];
  }

  //  ~Vector3D(){};

  //    Vector3D(const Vector3D& vec) :
  //      fX(vec.fX), fY(vec.fY), fZ(vec.fZ) {};

  // Vector3D& operator = (const Vector3D&);
  //
  // general field operators
  Vector3D& operator+=(const Vector3D&);
  Vector3D& operator-=(const Vector3D&);
  Vector3D operator-() const;
  Vector3D operator+() const;
  const double& operator[](size_t) const;
  double& operator[](size_t);

  const double& operator()(size_t) const;
  double& operator()(size_t);

  bool operator==(const Vector3D&) const;
  bool operator!=(const Vector3D&) const;

  // access methods
  double GetX() const { return fX; };
  double GetY() const { return fY; };
  double GetZ() const { return fZ; };

  double X() const { return fX; };
  double Y() const { return fY; };
  double Z() const { return fZ; };

  void SetX(double value) { fX = value; };
  void SetY(double value) { fY = value; };
  void SetZ(double value) { fZ = value; };
  //
  // vector operations
  double Norm() const;
  double NormSq() const;
  double Mag() const { return Norm(); }
  double Mag2() const { return NormSq(); }
  double Dot(const Vector3D&) const;
  Vector3D Cross(const Vector3D&) const;
  Vector3D Unit() const;
  double Theta() const;
  double Phi() const;
  double Perp() const;
  virtual void Print(Option_t* option = "") const;
  // Distance between points
  double Distance(const Vector3D& vect) const;

  Vector3D& operator*=(double a) {
    fX *= a;
    fY *= a;
    fZ *= a;
    return *this;
  }

 private:
  double fX;  //  1. Vector Element
  double fY;  //  2. Vector Element
  double fZ;  //  3. Vector Element

  friend Vector3D operator+(const Vector3D& v1, const Vector3D& v2);
  friend Vector3D operator-(const Vector3D& v1, const Vector3D& v2);
  friend Vector3D operator*(double a, const Vector3D& v2);
  friend Vector3D operator*(const Vector3D& v1, double a);
  friend Vector3D operator/(const Vector3D& v1, double a);

 public:
  //  ClassDef(Vector3D,1)
};

//
// Related Global functions
//

ostream& operator<<(ostream& str, const Vector3D&);
istream& operator>>(istream& str, Vector3D&);

Vector3D Cross(const Vector3D&, const Vector3D&);

#endif

//
//  $Log: Vector3D.h,v $
//  Revision 1.1  2007/05/24 07:45:42  dc
//  Initial release as part of util
//
//  Revision 1.1.1.1  2006/11/03 17:11:49  dc
//  import dc stuff
//
