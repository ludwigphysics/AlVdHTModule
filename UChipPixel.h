//-*-mode:c++-*-

#ifndef U_UChipPixel
#define U_UChipPixel

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class UChipPixel : public TObject {
 public:
  UChipPixel();
  UChipPixel(const int iL, const int iC);  // fast use constructor
  virtual ~UChipPixel();

  int GetLine() { return fLine; }
  int GetColumn() { return fColumn; }

  // returns position of pixel centre
  // in local sensor frame. x=0,y=0 refers to center of the sensor
  // in the standard cartasian base: z+(beam dir.), x+ (jura), y+ (up)

  // For negative (Saleve) arm: 
  // we should add offset in x let the "0" be located in the whole sensor centre 
  // (not in the centre of the active region)
  double GetNaY() { return (fColumn - 512) * fPitchY - fPitchY/2.0; }
  //double GetNaX() { return (fLine-256) * fPitchX + fPitchX/2.0; }
  double GetNaX() { return (256-fLine) * fPitchX + fPitchX/2.0 + fOffX; }

  // For positive (Jura) arm:
  double GetPaY() { return (fColumn - 512) * fPitchY - fPitchY/2.0; }
  //double GetPaX() { return (256 - fLine) * fPitchX - fPitchX/2.0; }
  double GetPaX() { return (fLine-256) * fPitchX - fPitchX/2.0 - fOffX; }

  void SetLine(const int i) { fLine = i; }
  void SetColumn(const int i) { fColumn = i; }

  void SetUsed() { fUsed = true; }
  bool IsUsed() { return fUsed; }

 private:
  double fOffX; // 0.6 = 1.2/2 (to shift in x to keep sensor pos in frame located in center of the whole sensor)
  int fLine;    //
  int fColumn;  //
  double fPitchX;
  double fPitchY;
  bool fUsed;

 public:
  friend ostream& operator<<(ostream& os, UChipPixel* pixel);

  //  ClassDef(UChipPixel,1)  //  QDC data class
};

#endif
