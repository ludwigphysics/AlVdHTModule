//____________________________________________________________________
//
//  BrLine3D
//
//  Line class for Geometry. A line is defined by an origin and a unit
//  vector defining its direction.
//

#include "ULine3D.h"
//#include "Iostream.h"
#include <iostream>

#ifndef ROOT_TMath
#include "TMath.h"
#endif

//____________________________________________________________________
// ClassImp(Line3D);

//____________________________________________________________________
Line3D::Line3D() {
  // Default constructor - empty
}

//____________________________________________________________________
Line3D::Line3D(const Vector3D& origin, const Vector3D& direction) {
  //
  //  Constructor
  //
  fOrigin = origin;
  fDirection = direction.Unit();
}

//____________________________________________________________________
Line3D::Line3D(const double ox, const double oy, const double oz, const double dx, const double dy, const double dz) {
  Vector3D origin(ox, oy, oz);
  Vector3D direction(dx, dy, dz);

  fOrigin = origin;
  fDirection = direction.Unit();
}

//____________________________________________________________________
Line3D::Line3D(const Line3D& line) : UDataObject() {
  //
  //  Copy Constructor
  //
  fOrigin = line.GetOrigin();
  fDirection = line.GetDirection();
}

//____________________________________________________________________
Line3D::~Line3D() {
  // Dtor - empty
}

//____________________________________________________________________
Line3D& Line3D::operator=(const Line3D& line) {
  //
  // assignement operator
  //
  fOrigin = line.fOrigin;
  fDirection = line.fDirection;
  return *this;
}

//____________________________________________________________________
Vector3D Line3D::ClosestPoint(const Vector3D& vector) const {
  //
  // Return the closest point on the line to the vector
  //
  double value;
  Vector3D vtmp;
  vtmp = vector - fOrigin;
  value = vtmp.Dot(fDirection);
  vtmp = fOrigin + value * fDirection;
  return vtmp;
}

//____________________________________________________________________
double Line3D::Distance(const Vector3D& vector) const {
  //
  // distance from a vector (i.e. point) to this line.
  //
  double value, value1;
  Vector3D vtmp;
  vtmp = vector - fOrigin;
  value1 = vtmp.Dot(fDirection);
  value = vtmp.NormSq() - value1 * value1;
  return TMath::Sqrt(value);
}

//____________________________________________________________________
Vector3D Line3D::GetClosestProximityPoint(const Line3D& line) const {
  Line3D shortestLine = GetShortestLineBetween(line);
  Vector3D start = shortestLine.GetOrigin();
  Vector3D dir = shortestLine.GetDirection();
  double dist = GetShortestDistanceBetween(line);
  // cout<<"start: "<<start.X()<<" "<<start.Y()<<" "<<start.Z()<<endl;
  // cout<<"dir  : "<<dir.X()<<" "<<dir.Y()<<" "<<dir.Z()<<" "<<dir.Norm()<<" "<<dist<<endl;
  Vector3D cpp = start + 0.5 * dist * dir;
  return cpp;
}

//____________________________________________________________________
Vector3D Line3D::GetClosestProximityPoint(const Line3D* line_p) const {
  Line3D line(line_p->GetOrigin(), line_p->GetDirection());
  Line3D shortestLine = GetShortestLineBetween(line);
  Vector3D start = shortestLine.GetOrigin();
  Vector3D dir = shortestLine.GetDirection();
  double dist = GetShortestDistanceBetween(line);
  // cout<<"start: "<<start.X()<<" "<<start.Y()<<" "<<start.Z()<<endl;
  // cout<<"dir  : "<<dir.X()<<" "<<dir.Y()<<" "<<dir.Z()<<" "<<dir.Norm()<<" "<<dist<<endl;
  Vector3D cpp = start + 0.5 * dist * dir;
  return cpp;
}

//____________________________________________________________________
Line3D Line3D::GetShortestLineBetween(const Line3D& line) const {
  //
  // Return the shortest line between the two lines.  The length of
  // the vector part is the shortest distance between the two lines.
  // put in according to
  // http://www.mhri.edu.au/~pdb/geometry/lineline3d/ 1/31/99
  // but not tested yet?  Be careful before using in production mode.
  // dmnop = (xm-xn)(xo-xp) + (ym-yn)(yo-yp) + (zm-zn)(zo-zp)
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  double d1343, d4321, d1321, d4343, d2121;
  double ua, ub;
  x1 = fOrigin[0];
  x2 = x1 + fDirection[0];
  y1 = fOrigin[1];
  y2 = y1 + fDirection[1];
  z1 = fOrigin[2];
  z2 = z1 + fDirection[2];
  x3 = line.GetOrigin()[0];
  x4 = x3 + line.GetDirection()[0];
  y3 = line.GetOrigin()[1];
  y4 = y3 + line.GetDirection()[1];
  z3 = line.GetOrigin()[2];
  z4 = z3 + line.GetDirection()[2];
  d1343 = (x1 - x3) * (x4 - x3) + (y1 - y3) * (y4 - y3) + (z1 - z3) * (z4 - z3);
  d4321 = (x4 - x3) * (x2 - x1) + (y4 - y3) * (y2 - y1) + (z4 - z3) * (z2 - z1);
  d1321 = (x1 - x3) * (x2 - x1) + (y1 - y3) * (y2 - y1) + (z1 - z3) * (z2 - z1);
  d4343 = (x4 - x3) * (x4 - x3) + (y4 - y3) * (y4 - y3) + (z4 - z3) * (z4 - z3);
  d2121 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1);
  ua = (d1343 * d4321 - d1321 * d4343) / (d2121 * d4343 - d4321 * d4321);
  ub = (d1343 + ua * d4321) / d4343;
  Vector3D Pa = fOrigin + ua * fDirection.Unit();
  Vector3D Pb = line.GetOrigin() + ub * line.GetDirection().Unit();
  Line3D shortline(Pa, Pb - Pa);
  Vector3D start = shortline.GetOrigin();
  Vector3D dir = shortline.GetDirection();
  double dist = GetShortestDistanceBetween(line);
  Vector3D cpp = start + (0.5 * dist) * dir;
  // cout<<"X "<<cpp.X()<<" "<<(Pa.X()+Pb.X())/2.<<" "<<dir.Norm()<<endl;
  // cout<<"Y "<<cpp.Y()<<" "<<(Pa.Y()+Pb.Y())/2.<<endl;
  // cout<<"Z "<<cpp.Z()<<" "<<(Pa.Z()+Pb.Z())/2.<<endl;
  return shortline;
}

//____________________________________________________________________
double Line3D::GetShortestDistanceBetween(const Line3D& line) const {
  // Get shortest distance between the two lines.  This should give the
  // same answer as the modulus of GetShortestLineBetween().  Coded
  // 1-31-99 but not tested.
  Vector3D N1 = fDirection;
  Vector3D N2 = line.GetDirection();
  Vector3D normal = N1.Cross(N2);

  const double small = 1.e-15;
  // Check whether lines are parallel to each other by looking at
  // magnitude of cross product.
  if (TMath::Abs(normal.Norm()) < small) {
    // This means lines are parallel to each other
    // So get perpendicular distance between origin of line 1 to the other line
    return Distance(line.GetOrigin());
  }

  // Vector3D unit = N1.Cross(N2).Unit();
  // cout<<normal.Norm()<<endl;
  Vector3D unit = normal.Unit();
  // cout<<unit.Norm()<<endl;

  // should give vector between origins with proper length
  Vector3D temp = fOrigin - line.GetOrigin();
  // cout<<temp.X()<<" "<<temp.Y()<<" "<<temp.Z()<<" "<<temp.Dot(unit)<<endl;
  return TMath::Abs(temp.Dot(unit));
}

//____________________________________________________________________
double Line3D::GetShortestDistanceBetween(const Line3D* line_p) const {
  // Get shortest distance between the two lines.  This should give the
  // same answer as the modulus of GetShortestLineBetween().  Coded
  // 1-31-99 but not tested.
  Line3D line(line_p->GetOrigin(), line_p->GetDirection());
  Vector3D N1 = fDirection;
  Vector3D N2 = line.GetDirection();
  Vector3D normal = N1.Cross(N2);

  const double small = 1.e-15;
  // Check whether lines are parallel to each other by looking at
  // magnitude of cross product.
  if (TMath::Abs(normal.Norm()) < small) {
    // This means lines are parallel to each other
    // So get perpendicular distance between origin of line 1 to the other line
    return Distance(line.GetOrigin());
  }

  // Vector3D unit = N1.Cross(N2).Unit();
  // cout<<normal.Norm()<<endl;
  Vector3D unit = normal.Unit();
  // cout<<unit.Norm()<<endl;

  // should give vector between origins with proper length
  Vector3D temp = fOrigin - line.GetOrigin();
  // cout<<temp.X()<<" "<<temp.Y()<<" "<<temp.Z()<<" "<<temp.Dot(unit)<<endl;
  return TMath::Abs(temp.Dot(unit));
}

//____________________________________________________________________
void Line3D::Print(Option_t* /*option*/) const {
  //
  // Simple output of this line
  //
  std::cout << "Origin   : ";
  fOrigin.Print();
  std::cout << "Diection : ";
  fDirection.Print();
}

//____________________________________________________________________
double Line3D::RelativeOverlap(const Line3D& line, double zl, double r) {
  //
  // Calculate the relative overlap between "this" track and the given track
  //
  // Parameters:
  //   In:  track   pointer to the other track
  //        zl      length of detector (in local z direction)
  //        r       radius of cylinder around tracks
  //
  // Algorithm:
  //   A cylinder of radius r is put around each track.
  //   The relative overlap is defined as the ratio between
  //     the volume of the overlap region of two cylinders
  //   and
  //     the volume of the cylinder around "this" track is calculated.
  //   The cylinders are distorted in order to have a circular intersection
  //   between the cylinders and planes at arbitrary z.
  //
  // Requires:
  //   The tracks are parametrized as:
  //     x = ax * z + bx
  //     y = ay * z + by
  //

  const double slopeZthis = fDirection[2];
  const double slopeZline = line.GetDirection()[2];
  if (slopeZthis == 0 || slopeZline == 0) {
    Warning("RelativeOverlap", "z slope == 0");
    return 0;
  }

  const double twor = 2 * r;
  const double halfl = zl / 2;

  // Difference between the two tracks
  // Note that we normalise to slope z == 1
  const double dax = fDirection[0] / slopeZthis - line.GetDirection()[0] / slopeZline;
  const double day = fDirection[1] / slopeZthis - line.GetDirection()[1] / slopeZline;
  const double dbx = fOrigin[0] - line.GetOrigin()[0];
  const double dby = fOrigin[1] - line.GetOrigin()[1];

  double overlap = 0;

  const double a = dax * dax + day * day;

  if (a == 0) {
    // Tracks are parallel
    const double d = TMath::Sqrt(dbx * dbx + dby * dby) / twor;
    if (d < 1) overlap = TMath::ACos(d) - d * TMath::Sqrt(1 - d * d);
  } else {
    const double b = (dax * dbx + day * dby) / a;
    const double c = (dbx * dbx + dby * dby) / a;
    double d = b * b - c + twor * twor / a;
    if (d > 0) {
      // Cylinders intersect
      d = TMath::Sqrt(d);
      double zn = -b - d;
      double zx = -b + d;
      if (zn < -halfl) zn = -halfl;
      if (zx > halfl) zx = halfl;
      if (zn < zx) {
        // Cylinders intersect inside the detector
        //
        // Integrate
        //   acos(d(z)/(2r))-d(z)/(2r)*sqrt(1-(d(z)/(2r))**2)
        // where
        //   d(z)/(2r)=sqrt(a*((z+2*b)*z+c))/(2*r)
        // from zn to zx
        //
        // Quick and dirty implementation: Simple N-point rule
        const int n = 5;
        const double dz = (zx - zn) / n;
        for (int i = 0; i < n; i++) {
          double z = zn + dz * (i + 0.5);
          d = TMath::Sqrt(a * ((z + 2 * b) * z + c)) / twor;
          overlap += TMath::ACos(d) - d * TMath::Sqrt(1 - d * d);
        }

        overlap *= dz / zl;
      }
    }
  }

  return overlap / (TMath::Pi() / 2);
}

//____________________________________________________________________
void Line3D::SetDirection(const Vector3D& vector) {
  // set x, y, z
  const double norm = vector.Norm();
  fDirection.SetX(vector.GetX() / norm);
  fDirection.SetY(vector.GetY() / norm);
  fDirection.SetZ(vector.GetZ() / norm);
}

//____________________________________________________________________
void Line3D::SetOrigin(const Vector3D& vector) {
  // set x, y, z
  fOrigin.SetX(vector.GetX());
  fOrigin.SetY(vector.GetY());
  fOrigin.SetZ(vector.GetZ());
}

//____________________________________________________________________
ostream& operator<<(ostream& os, const Line3D& line) {
  //
  // Related Function (global non-member functions)
  // Output line parameters in default layout
  //
  os << line.GetOrigin() << "+ t*" << line.GetDirection();
  return os;
}

//____________________________________________________________________
ostream& operator<<(ostream& os, const Line3D* line) {
  os << *line;
  return os;
}

//____________________________________________________________________
//
// $Log: $
