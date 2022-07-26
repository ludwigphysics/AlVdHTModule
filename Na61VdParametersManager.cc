//////////////////////////////////////////////////////////////////
//
//  Na61VdParametersManager
//
// This manager class is use to provide SAVD geometry and matching
// parameters over VD reconstruction modules.
// As the manager is singleton it ensures same parameter are used
// over the reconstructionn

//
#include "Na61VdParametersManager.h"

#ifndef __IOSTREAM__
#include <iostream>
#endif
#ifndef __IOMANIP__
#include <iomanip>
#endif
using std::cout;
using std::endl;
using std::cerr;

//____________________________________________________
// ClassImp(Na61VdParametersManager);

//____________________________________________________
//
// Static instance of
//
Na61VdParametersManager* Na61VdParametersManager::fInstance = 0;

//____________________________________________________
//
Na61VdParametersManager::Na61VdParametersManager() {
  fJuraArmParams = new Na61ArmParameters("Jura");
  fSaleveArmParams = new Na61ArmParameters("Saleve");
  fVdParams = new Na61VdParameters("SAVD");

  fJuraAlVdArmParams   = new Na61AlVdArmParameters("Jura");
  fSaleveAlVdArmParams   = new Na61AlVdArmParameters("Saleve");
  // Init(); // parameter initialization call by constructor
}

//____________________________________________________
Na61VdParametersManager* Na61VdParametersManager::Instance() {
  if (fInstance == 0) fInstance = new Na61VdParametersManager;
  return fInstance;
}
//__________________________________________________

Na61VdParametersManager::~Na61VdParametersManager() {}

//________________________________________________________________
void Na61VdParametersManager::Init() {
  // Initialize VD geometry and matching params

  fJuraArmParams->Init();
  fSaleveArmParams->Init();
  fVdParams->Init();

  fJuraAlVdArmParams   -> Init();
  fSaleveAlVdArmParams -> Init();

}

//________________________________________________________________
void Na61VdParametersManager::Update() {
  // Update the geometry
}
