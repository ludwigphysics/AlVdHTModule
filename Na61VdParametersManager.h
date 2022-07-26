// -*- mode: c++ -*-
#ifndef Na61_Na61VdParametersManager
#define Na61_Na61VdParametersManager
//
#ifndef Na61_Na61ArmParameters
#include "Na61ArmParameters.h"
#endif
#ifndef Na61_Na61AlVdArmParameters
#include "Na61AlVdArmParameters.h"
#endif
#ifndef Na61_Na61VdParameters
#include "Na61VdParameters.h"
#endif

// Root Classes
#ifndef ROOT_TObject
#include "TObject.h"
#endif

class Na61VdParametersManager : public TObject {
 public:
  static Na61VdParametersManager* Instance();

  Na61VdParametersManager();
  ~Na61VdParametersManager();

  Na61ArmParameters* GetJuraArmParams() { return fJuraArmParams; }
  Na61ArmParameters* GetSaleveArmParams() { return fSaleveArmParams; }
  Na61VdParameters* GetVdParams() { return fVdParams; }

  Na61AlVdArmParameters* GetJuraAlVdArmParams(){return fJuraAlVdArmParams;}
  Na61AlVdArmParameters* GetSaleveAlVdArmParams(){return fSaleveAlVdArmParams;}

 public:
  void Init();
  void Update();

 private:
  Na61ArmParameters* fJuraArmParams;
  Na61ArmParameters* fSaleveArmParams;
  Na61VdParameters* fVdParams;

  Na61AlVdArmParameters* fJuraAlVdArmParams;
  Na61AlVdArmParameters* fSaleveAlVdArmParams;

  static Na61VdParametersManager* fInstance;

  //  ClassDef(Na61VdParametersManager,1)
};

#endif
