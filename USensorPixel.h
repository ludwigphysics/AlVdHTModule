//-*-mode:c++-*-

#ifndef U_USensorPixel
#define U_USensorPixel

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class USensorPixel : public TObject {
 public:
  USensorPixel();
  USensorPixel(const int iL, const int iC);  // fast use constructor
  virtual ~USensorPixel();

  int GetLine() { return fLine; }
  int GetColumn() { return fColumn; }

  // returns position of pixel centre
  // in local sensor frame. x=0,y=0 refers to center of the sensor
  // in the standard cartasian base: z+(beam dir.), x+ (jura), y+ (up)

  // For negative (Saleve) arm:
  double GetNaY() { return (fColumn - 576) * fPitch - fPitch / 2.0 + fPitch; }
  double GetNaX() { return (288 - fLine) * fPitch + fPitch / 2.0 - fPitch; }

  // For positive (Jura) arm:
  double GetPaY() { return (576 - fColumn) * fPitch + fPitch / 2.0 - fPitch; }
  double GetPaX() { return (fLine - 288) * fPitch - fPitch / 2.0 + fPitch; }

  void SetLine(const int i) { fLine = i; }
  void SetColumn(const int i) { fColumn = i; }

  void SetUsed() { fUsed = true; }
  bool IsUsed() { return fUsed; }

 private:
  int fLine;    //
  int fColumn;  //
  double fPitch;
  bool fUsed;

 public:
  friend ostream& operator<<(ostream& os, USensorPixel* pixel);

  //  ClassDef(USensorPixel,1)  //  QDC data class
};

#endif
