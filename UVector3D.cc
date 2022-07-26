//  $Id: Vector3D.cpp,v 1.1 2007/05/24 07:45:42 dc Exp $
//
//
#ifndef UTIL_UVector3D
#include "UVector3D.h"
#endif
#include "math.h"

//#ifndef UTIL_Iostream
//#include "Iostream.h"
//#endif
#include <iostream>
#include <iomanip>

//////////////////////////////////////////////////////////////////////////////////
//
//
//  BrVector3D defines a general 3-vector class. It can be used to represent space
//  points, directions, and momenta.BrVector3D is part of the geometry classes.
//  This is a general class only coupled to ROOT by having the objects being derived from
//  TObject. The main reasons for doing this is to be able to use the ROOTs
//  interactive features.
//  It is intended to be used
//  by detector geometry , tracking, and display classes by using common
//  geometry concepts. It is also intended to work together with other geometry clasess
//  like BrLine3D, BrPlane3D and likely other
//  The coupling to ROOT is loose. The classes are derived from TObject, mainly
//  to be able to use browsers and other general ROOT utilities. Drawing has not
//  been implemented here again to maintain a very loose coupling.
//
//  Vectors are intrinsicly represented by doubles, particular to maintain precesion
//  in calculation of angles, intercepts.
//
//  Such general classes has of course been implemented by numerous people, and
//  design choices made. As such example the CLHEP should be acknowledged.
//
//
//
/////////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________
// ClassImp(Vector3D);

const double& Vector3D::operator[](size_t i) const {
  if (i == 0)
    return fX;
  else if (i == 1)
    return fY;
  else
    return fZ;
}

double& Vector3D::operator[](size_t i) {
  if (i == 0)
    return fX;
  else if (i == 1)
    return fY;
  else
    return fZ;
}

const double& Vector3D::operator()(size_t i) const {
  if (i == 0)
    return fX;
  else if (i == 1)
    return fY;
  else
    return fZ;
}

double& Vector3D::operator()(size_t i) {
  if (i == 0)
    return fX;
  else if (i == 1)
    return fY;
  else
    return fZ;
}

//  Vector3D & Vector3D::operator = (const Vector3D& vec){
//    //
//    // assignement operator
//    //
//    fX = vec.fX;
//    fY = vec.fY;
//    fZ = vec.fZ;
//    return *this;
//  }

Vector3D operator+(const Vector3D& vec1, const Vector3D& vec2) {
  Vector3D tmp;
  tmp.fX = vec1.fX + vec2.fX;
  tmp.fY = vec1.fY + vec2.fY;
  tmp.fZ = vec1.fZ + vec2.fZ;
  return tmp;
}

Vector3D operator-(const Vector3D& vec1, const Vector3D& vec2) {
  Vector3D tmp;
  tmp.fX = vec1.fX - vec2.fX;
  tmp.fY = vec1.fY - vec2.fY;
  tmp.fZ = vec1.fZ - vec2.fZ;
  return tmp;
}

Vector3D& Vector3D::operator+=(const Vector3D& vec) {
  fX += vec.fX;
  fY += vec.fY;
  fZ += vec.fZ;
  return *this;
}

Vector3D& Vector3D::operator-=(const Vector3D& vec) {
  fX -= vec.fX;
  fY -= vec.fY;
  fZ -= vec.fZ;
  return *this;
}

Vector3D Vector3D::operator-() const {
  Vector3D tmp(-fX, -fY, -fZ);
  return tmp;
}

Vector3D Vector3D::operator+() const { return *this; }

bool Vector3D::operator==(const Vector3D& vec2) const {
  // Boolean equal
  return ((fX == vec2.fX) && (fY == vec2.fY) && (fZ == vec2.fZ));
}

bool Vector3D::operator!=(const Vector3D& vec2) const {
  // Boolean not equal
  return !(*this == vec2);
}

//______________________________________________________________
double Vector3D::Norm() const {
  //
  // Magnitude of norm of vector
  //
  double tmp;
  tmp = fX * fX + fY * fY + fZ * fZ;
  return sqrt(tmp);
}

//______________________________________________________________
double Vector3D::NormSq() const {
  //
  // Squared Magnitude of vector
  //
  return fX * fX + fY * fY + fZ * fZ;
}

//______________________________________________________________
double Vector3D::Dot(const Vector3D& vec) const {
  //
  // Dot product of two vector
  //
  double tmp;
  tmp = fX * vec.fX + fY * vec.fY + fZ * vec.fZ;
  return tmp;
}

//______________________________________________________________
Vector3D Vector3D::Cross(const Vector3D& vec) const {
  //
  // Cross product of two vectors
  //
  Vector3D tmp(fY * vec.fZ - fZ * vec.fY, fZ * vec.fX - fX * vec.fZ, fX * vec.fY - fY * vec.fX);
  return tmp;
}

//______________________________________________________________
Vector3D Vector3D::Unit() const {
  if (this->Norm() != 0)
    return *this / (this->Norm());
  else
    return *this;
}

//______________________________________________________________
double Vector3D::Perp() const {
  //
  // transverse component
  //
  if (fX == 0.0 && fY == 0.0) {
    return 0;
  } else {
    return sqrt(fX * fX + fY * fY);
  }
}

//______________________________________________________________
double Vector3D::Theta() const {
  //
  // Polar angle with regard to z axis
  //
  if (fX == 0.0 && fY == 0.0) {
    return 0;
  } else {
    double perp = sqrt(fX * fX + fY * fY);
    return atan2(perp, fZ);
  }
}

//______________________________________________________________
double Vector3D::Phi() const {
  //
  // Azimuth  angle with regard to z axis
  //
  if (fX == 0.0 && fY == 0.0) {
    return 0;
  } else
    return atan2(fY, fX);
}

//______________________________________________________________
void Vector3D::Print(Option_t* /*option*/) const { std::cout << "(" << std::setw(10) << fX << ", " << std::setw(10) << fY << ", " << std::setw(10) << fZ << ")" << std::endl; }

//______________________________________________________________
double Vector3D::Distance(const Vector3D& vec) const {
  double dist;
  dist = sqrt((fX - vec.GetX()) * (fX - vec.GetX()) + (fY - vec.GetY()) * (fY - vec.GetY()) + (fZ - vec.GetZ()) * (fZ - vec.GetZ()));
  return dist;
}

//
//  friend global functions
//  =======================

Vector3D operator*(double a, const Vector3D& vec) {
  //
  //  number multiply
  //
  Vector3D tmp;
  tmp.fX = a * vec.fX;
  tmp.fY = a * vec.fY;
  tmp.fZ = a * vec.fZ;
  return tmp;
}

Vector3D operator/(const Vector3D& vec, double a) {
  //
  //  number multiply
  //
  Vector3D tmp;
  tmp.fX = vec.fX / a;
  tmp.fY = vec.fY / a;
  tmp.fZ = vec.fZ / a;
  return tmp;
}

Vector3D operator*(const Vector3D& vec, double a) {
  //
  //  number multiply
  //
  Vector3D tmp;
  tmp.fX = a * vec.fX;
  tmp.fY = a * vec.fY;
  tmp.fZ = a * vec.fZ;
  return tmp;
}

//
// Related Functions (global non-member functions)
//

ostream& operator<<(ostream& os, const Vector3D& vec) {
  //
  // Write in default layout
  //
  // os << setprecision(9) <<"("
  os << std::setprecision(6) << std::setiosflags(std::ios::fixed) << "(" << vec.GetX() << ", " << vec.GetY() << ", " << vec.GetZ() << ")";
  return os;
}

istream& operator>>(istream& is, Vector3D& vec) {
  //
  // Input stream operator
  //
  double x, y, z;
  is >> x >> y >> z;
  vec.SetX(x);
  vec.SetY(y);
  vec.SetZ(z);
  return is;
}

Vector3D Cross(const Vector3D& vec1, const Vector3D& vec2) {
  //
  //  Cross product of two vectors
  //
  return vec1.Cross(vec2);
}

double Dot(const Vector3D& vec1, const Vector3D& vec2) {
  //
  //  Dot product of two vectors
  //
  return vec1.Dot(vec2);
}

//  $Log: Vector3D.cpp,v $
//  Revision 1.1  2007/05/24 07:45:42  dc
//  Initial release as part of util
//
//  Revision 1.1.1.1  2006/11/03 17:11:49  dc
//  import dc stuff
//
