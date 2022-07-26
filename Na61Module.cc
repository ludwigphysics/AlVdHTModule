//____________________________________________________________________
//
// VD STANDALONE Analysis Module Class
//
// A standard Physics analysis Module should be derived from
// this base class. This ensure that the general framework
// is capable of accessing these in a common and general way.
// The methods (member functions) consists of a general set for
// all modules; but the user if free to add specific ones for
// specific needs. Ideas for this has been taken from the
// CDF/E802 Analysis Control as well as from the BABAR framework.
// The initial version is derived from the ROOT base classes,
//
//
// The standard methods defined for a module are:
// ----------------------------------------------
//
//  Init()
//        This methods finishes the initialisation of the module, that
//        is after all user parameters for the module has been
//        set. This message is sent once per job. External
//        parameters should be referenced at this time.
//        A method intended to be used for modules that cannot be
//        completely initialized by the constrcutor and needs
//        additional setups. It is the responsibility of the user
//        modules to give out warnings if the init has not been
//        called.
//
//  Book()
//        Book histograms if defined in the DefineHistograms() method,
//        The booking method ensures that the histograms are added to the
//        list of known histograms.
//
//  Begin()
//        Send at run initialisation time. That is once per run. This
//        method can for example update references to external
//        parameters, such as database entries and so on.
//
//  Event(BrEventNode* , BrEventNode*)
//         A module may use both nodes. The normal convention is that
//         the first node is input i.e. module should look for
//         dataobjects in this node, The second node may be
//         empty. Normal use is to add dataobjects to  this. The main
//         program/shell will control the vent nodes. All data objects
//         added to the node will be owned by the node. Thus the node
//         will delete at end of event.
//         Generated DataObjects or BrDataTable s should be added to the
//         outputnode. The only excpection is modues derived from BrIoModule
//
//
//  End()
//        Send at the end of a run. This should flush histograms,
//        etc. specific for this run to disk and so on. Also, one may
//        need to tidy up some stuff before the next Begin.
//
//  Finish()
//        Send at end of a job. This should do the final write of data
//        to permanent store and so on.
//
//
// Flow in a module:
// -----------------
//
//                  +-------------------+
//                  | Set up job things |
//                  +-------------------+
//                            |
//                  +--------------------+
//                  | Initialize job     |
//                  | [BrModule::Init()] |
//                  +--------------------+
//                            |
//            no  o------------------------o
//          +-<---| More runs to process ? |---<----------+
//          |     o------------------------o              |
//          |                 | yes                       |
//          |       +---------------------+     +-------------------+
//          |       | Initialize for run  |     | Finish run        |
//          |       | [BrModule::Begin()] |     | [BrModule::End()] |
//          |       +---------------------+     +-------------------+
//          |                 |                           |
//          |     o-------------------------o  no         |
//          |  +--| More events to process? |-->----------+
//          |  |  o-------------------------o
//          |  ^              | yes
//          |  |    +---------------------+
//          |  +----| Process event       |
//          |       | [BrModule::Event()] |
//          |       +---------------------+
//          |
//          |
//    +----------------------+
//    | Finish job           |
//    | [BrModule::Finish()] |
//    +----------------------+
//
//
// Error handling and state information:
// -------------------------------------
//
// A few  methods are avaliable for communicating the status of the
// module.
//
//   Stop(<location>, <format string>, ...)
//     Flag this module as having not being able to fullfil it's
//     requirements. That is, the module pipeline should stop after
//     this module. Note, sending this message in Init() will stop the
//     current job altogether, and in Begin() will stop that run loop
//     altogether. If verbose flag is greater then, or equal to
//     kShowStop, then also output message. If this failure makes the
//     total number of failures for this module exceed the stop limit
//     (as set by BrModule::SetStopLimit), then promote this stop to a
//     Failure (see BrModule::Failure).
//
//   Failure(<location>, <format string>, ...)
//     Flag this module as having failed. If verbose flag is greater
//     then, or equal to kShowFailure, then also output message. If
//     this failure makes the total number of failures for this module
//     exceed the failure limit (as set by BrModule::SetFailureLimit),
//     then promote this failure to a Abort (see
//     BrModule::Abort).  Normally a Failure should not stop the
//     event processing, but should stop the module from doing more
//     work.
//
//   Abort<location>, <format string>, ...)
//     Flag this module as having failed misaerably, and output
//     message. A BrModuleContainer should always react to this by
//     stoping the module loop, since this flag means that something
//     is very VERY screwy. This should abort the event processing all
//     together.
//
// Notice that TObject::Warning, TObject::Error, and TObject::Fatal
// are reserved for serious errors.
//
//____________________________________________________________________

//
// $Id: BrModule.cxx,v 1.15 2009/10/02 19:13:28 videbaek Exp $
// $Author: videbaek $
// $Date: 2009/10/02 19:13:28 $
//
//  Update History:
//    FV March 28 ,1999
//     Changed HistOn, HistOff, Dis/EnableHistograms to be non-virtual
//     functions. It is not not up to the user to modify this behaviou.
//    FV April 19,1999
//     Added methods SetHistOn() , SetHistOff(), and SetHistBooked()
//

#ifndef NA61_Na61Module
#include "Na61Module.h"
#endif
#ifndef UTIL_UEvent
#include "UEvent.h"
#endif
#ifndef ROOT_TClass
#include "TClass.h"
#endif
#ifndef ROOT_TError
#include "TError.h"
#endif
#ifndef ROOT_TList
#include "TList.h"
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(3, 0, 1)
#include "TStreamerInfo.h"
#endif
#ifndef ROOT_TDirectory
#include "TDirectory.h"
#endif

#include <string>
#include <time.h>
using namespace std;

//#ifdef __GNUC__ && __GNUC__ >= 3
#define va_(X) X
//#endif

//____________________________________________________________________
// ClassImp(Na61Module);

//____________________________________________________________________
Na61Module::Na61Module() {
  // Vanilla constructor; Should  not be called within the ROOT
  // environment.
  //
  fVerbose = 0;
  fHistograms = 0;
  fDebugLevel = 0;
  fCpuTime = 0;
  fHistOn = false;
  fHistBooked = false;
  fRequiredTableList = 0;
  fState = kSetup;
  fStatus = kOk;
  fStopCount = 0;
  fFailureCount = 0;
  fStopLimit = 2147483647;  // The largest possible int
  fFailureLimit = 10;
  fDetectorSetupOption = 0;
}

//____________________________________________________________________
Na61Module::Na61Module(const char *name, const char *title) : TNamed(name, title) {
  // Constructor. This is the default constructor to use. The names of
  // the Modules can be used for navigation and identification.
  fVerbose = 0;
  fHistograms = 0;
  fDebugLevel = 0;
  fCpuTime = 0;
  fHistOn = false;
  fHistBooked = false;
  fRequiredTableList = 0;
  fState = kSetup;
  fStatus = kOk;
  fStopCount = 0;
  fFailureCount = 0;
  fStopLimit = 2147483647;  // The largest possible int
  fFailureLimit = 10;
  fDetectorSetupOption = 0;
  //
  fAlSensorNames[0]="Al1_0"; fAlSensorNames[1]="Al1_1"; fAlSensorNames[2]="Al1_2";
  // Vds 2                                                                                                                               
  fAlSensorNames[3]="Al2_0"; fAlSensorNames[4]="Al2_1"; fAlSensorNames[5]="Al2_2"; // inner station                                      
  fAlSensorNames[6]="Al2_3"; fAlSensorNames[7]="Al2_4"; fAlSensorNames[8]="Al2_5"; // outer station                                      
  // Vds 3        
  fAlSensorNames[9]="Al3_0"; fAlSensorNames[10]="Al3_1"; fAlSensorNames[11]="Al3_2"; fAlSensorNames[12]="Al3_3"; fAlSensorNames[13]="Al3_4"; // inner station  
  fAlSensorNames[14]="Al3_5"; fAlSensorNames[15]="Al3_6"; fAlSensorNames[16]="Al3_7"; fAlSensorNames[17]="Al3_8"; fAlSensorNames[18]="Al3_9";// outer station   
  // Vds 4                                                                                                                               
  fAlSensorNames[19]="Al4_0"; fAlSensorNames[20]="Al4_1"; fAlSensorNames[21]="Al4_2"; fAlSensorNames[22]="Al4_3"; fAlSensorNames[23]="Al4_4";// inner station        
  fAlSensorNames[24]="Al4_5"; fAlSensorNames[25]="Al4_6"; fAlSensorNames[26]="Al4_7"; fAlSensorNames[27]="Al4_8"; fAlSensorNames[28]="Al4_9";// outer station  
  fAlSensorNames[29]="Al4_10";fAlSensorNames[30]="Al4_11";fAlSensorNames[31]="Al4_12";fAlSensorNames[32]="Al4_13";fAlSensorNames[33]="Al4_14";// outer station    

  fSensorNames[0] = "Vds1_0";
  fSensorNames[1] = "Vds2_0";
  fSensorNames[2] = "Vds3_0";
  fSensorNames[3] = "Vds3_1";
  fSensorNames[4] = "Vds4_0";
  fSensorNames[5] = "Vds4_1";
  fSensorNames[6] = "Vds4_2";
  fSensorNames[7] = "Vds4_3";

  fDevNames[0] = "dev1_0";
  fDevNames[1] = "dev1_1";
  fDevNames[2] = "dev2_0";
  fDevNames[3] = "dev2_1";
  fDevNames[4] = "dev2_2";
  fDevNames[5] = "dev2_3";

  fDevNames[6] = "dev1_0_x1";
  fDevNames[7] = "dev1_0_x2";
  fDevNames[8] = "dev1_0_x3";
  fDevNames[9] = "dev1_0_x4";
  fDevNames[10] = "dev1_1_x1";
  fDevNames[11] = "dev1_1_x2";
  fDevNames[12] = "dev1_1_x3";
  fDevNames[13] = "dev1_1_x4";
  fDevNames[14] = "dev1_2_x1";
  fDevNames[15] = "dev1_2_x2";
  fDevNames[16] = "dev1_2_x3";
  fDevNames[17] = "dev1_3_x1";
  fDevNames[18] = "dev1_3_x2";
  fDevNames[19] = "dev1_3_x3";
  //
  fMatchStr[0] = "down1";
  fMatchStr[1] = "up1";
  fMatchStr[2] = "down2";
  fMatchStr[3] = "up2";

  fMatchStr[4] = "down1_x1";
  fMatchStr[5] = "down1_x2";
  fMatchStr[6] = "down1_x3";
  fMatchStr[7] = "down1_x4";
  fMatchStr[8] = "up1_x1";
  fMatchStr[9] = "up1_x2";
  fMatchStr[10] = "up1_x3";
  fMatchStr[11] = "up1_x4";
  fMatchStr[12] = "down2_x1";
  fMatchStr[13] = "down2_x2";
  fMatchStr[14] = "down2_x3";
  fMatchStr[15] = "up2_x1";
  fMatchStr[16] = "up2_x2";
  fMatchStr[17] = "up2_x3";

  fSensorXwidth = 10.6;  // sensor x width in mm
  fSensorYwidth = 21.2;  // sensor y width in mm

  fElecM = 0.000511;
  fPionM = 0.13957;
  fKaonM = 0.493677;
  fProtM = 0.938272;

  fK0sM = 0.497611;
  fLambM = 1.115683;
  fDPlusM = 1.86962;
  fD0M = 1.86484;
}

//_______________________________________________________________
void Na61Module::Init() { 
  SetState(kInit); 
  
  // It is not the best solution but keep this initialisation here for the time being
  // Init DevNames
  /*
  TString devNames[] = {"A_00_00_00","A_00_00_0-1","A_00_00_01","A_00_0-1_0-1","A_00_01_01",
			"A_00_00_10","A_00_00_1-1","A_00_00_11","A_00_0-1_1-1","A_00_01_11",
			"A_00_10_10","A_00_10_1-1","A_00_10_11","A_00_1-1_1-1","A_00_11_11",
			"A_00_10_20","A_00_10_2-1","A_00_10_21","A_00_1-1_2-1","A_00_11_21",
			/////////////////////////////////////////////////////////////////////////
			/*"B_00_00_00","B_00_0-1_0-1","B_00_01_01","B_00_0-1_0-2","B_00_01_02",
			"B_00_00_10","B_00_0-1_1-1","B_00_01_11","B_00_0-1_1-2","B_00_01_12",
			"B_00_10_10","B_00_1-1_1-1","B_00_11_11","B_00_1-1_1-2","B_00_11_12",
			"B_00_10_20","B_00_1-1_2-1","B_00_11_21","B_00_1-1_2-2","B_00_11_22",
			"B_00_20_20","B_00_2-1_2-1","B_00_21_21","B_00_2-1_2-2","B_00_21_22",  //keep only A and D combinations 
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			"D_00_00_00","D_00_00_0-1","D_00_00_01",
			"D_0-1_0-1_0-1","D_01_01_01","D_0-1_0-1_0-2","D_01_01_02","D_0-1_0-2_0-2","D_01_02_02",
			"D_00_00_10","D_0-1_0-1_1-1","D_01_01_11","D_0-1_0-1_1-2","D_01_01_12","D_0-1_0-2_1-2","D_01_02_12",
			"D_00_10_10","D_0-1_1-1_1-1","D_01_11_11","D_0-1_1-1_1-2","D_01_11_12","D_0-1_1-2_1-2","D_01_12_12",
			"D_10_10_20","D_1-1_1-1_2-1","D_11_11_21","D_1-1_1-1_2-2","D_11_11_22","D_1-1_1-2_2-2","D_11_12_22",
			"D_10_20_20","D_1-1_2-1_2-1","D_11_21_21","D_1-1_2-1_2-2","D_11_21_22","D_1-1_2-2_2-2","D_11_22_22"};
  */

  TString devNames[] = {"A_00_00_00","A_00_00_0-1","A_00_00_01","A_00_0-1_0-1","A_00_01_01",
			"A_00_00_10","A_00_00_1-1","A_00_00_11","A_00_0-1_1-1","A_00_01_11",
			"A_00_10_10","A_00_10_1-1","A_00_10_11","A_00_1-1_1-1","A_00_11_11",
			"A_00_10_20",
			//////////////////////////////////////////////////////////////////////////////////////////////////////
			"D_00_00_00","D_00_00_0-1","D_00_00_01",
			"D_0-1_0-1_0-1","D_01_01_01","D_0-1_0-1_0-2","D_01_01_02","D_0-1_0-2_0-2","D_01_02_02",
			"D_00_00_10","D_0-1_0-1_1-1","D_01_01_11","D_0-1_0-1_1-2","D_01_01_12","D_0-1_0-2_1-2","D_01_02_12",
			"D_00_10_10","D_0-1_1-1_1-1","D_01_11_11","D_0-1_1-1_1-2","D_01_11_12","D_0-1_1-2_1-2","D_01_12_12",
			"D_10_10_20","D_1-1_1-1_2-1","D_11_11_21","D_1-1_1-1_2-2","D_11_11_22","D_1-1_1-2_2-2","D_11_12_22",
			// added
			"D_00_0-1_0-1","D_00_01_01","D_00_0-1_1-1","D_00_01_11","D_00_1-1_1-1","D_00_11_11","D_10_1-1_2-1","D_10_11_21",
			//
			"D_00_00_1-1","D_00_00_11","D_00_10_1-1","D_00_10_11","D_10_10_2-1","D_10_10_21"};

  fNDev = sizeof(devNames)/24;
  cout<<" Na61Module::Init(): fNDev = "<<fNDev<<endl; 
  for(Int_t i=0;i<fNDev;i++){
    fDevName[i]=devNames[i];
    if(devNames[i].Contains("A")){fVdsArr[i][0]=1; fVdsArr[i][1]=2; fVdsArr[i][2]=3;}
    if(devNames[i].Contains("B")){fVdsArr[i][0]=1; fVdsArr[i][1]=3; fVdsArr[i][2]=4;}
    if(devNames[i].Contains("C")){fVdsArr[i][0]=1; fVdsArr[i][1]=2; fVdsArr[i][2]=4;}
    if(devNames[i].Contains("D")){fVdsArr[i][0]=2; fVdsArr[i][1]=3; fVdsArr[i][2]=4;}
    Int_t ic[3];
    Int_t is[3];
    AnaSensorString(devNames[i],ic,is);
    fColArr[i][0] = ic[0]; fColArr[i][1] = ic[1]; fColArr[i][2] = ic[2];
    fSenArr[i][0] = is[0]; fSenArr[i][1] = is[1]; fSenArr[i][2] = is[2];
    //cout<<"-------->test: i = "<<i<<"   "<<devNames[i].Data()<<"   "<<ic[0]<<is[0]<<"_"<<ic[1]<<is[1]<<"_"<<ic[2]<<is[2]<<endl;

  }
  // define 4hit combination names
  Int_t imatch = 0;
  for(Int_t i=0;i<fNDev;i++){
    for(Int_t j=i+1;j<fNDev;j++){
   
      /// check login conditions ones more
      if(fDevName[i].Contains("A") && fDevName[j].Contains("B")){
	if(fColArr[i][0]==fColArr[j][0] && fColArr[i][2]==fColArr[j][1] &&
	   fSenArr[i][0]==fSenArr[j][0] && fSenArr[i][2]==fSenArr[j][1]){

	  fDevCombName[imatch] = Form("%s+%s",fDevName[i].Data(),fDevName[j].Data());
	  fDevi[imatch] = i;  fDevj[imatch] = j;
	  fMatchId1[imatch] = 0; fMatchId2[imatch] = 2;
	  imatch++;
	}
      }

      if(fDevName[i].Contains("A") && fDevName[j].Contains("C")){
	if(fColArr[i][0]==fColArr[j][0] && fColArr[i][1]==fColArr[j][1] &&
	   fSenArr[i][0]==fSenArr[j][0] && fSenArr[i][1]==fSenArr[j][1]){

	  fDevCombName[imatch] = Form("%s+%s",fDevName[i].Data(),fDevName[j].Data());
	  fDevi[imatch] = i;  fDevj[imatch] = j;
	  fMatchId1[imatch] = 0; fMatchId2[imatch] = 1;
	  imatch++;
	}
      }

      if(fDevName[i].Contains("A") && fDevName[j].Contains("D")){
	if(fColArr[i][1]==fColArr[j][0] && fColArr[i][2]==fColArr[j][1] &&
	   fSenArr[i][1]==fSenArr[j][0] && fSenArr[i][2]==fSenArr[j][1]){

	  fDevCombName[imatch] = Form("%s+%s",fDevName[i].Data(),fDevName[j].Data());
	  fDevi[imatch] = i;  fDevj[imatch] = j;
	  fMatchId1[imatch] = 1; fMatchId2[imatch] = 2;
	  imatch++;
	}
      }

      if(fDevName[i].Contains("B") && fDevName[j].Contains("C")){
	if(fColArr[i][0]==fColArr[j][0] && fColArr[i][2]==fColArr[j][2] &&
	   fSenArr[i][0]==fSenArr[j][0] && fSenArr[i][2]==fSenArr[j][2]){

	  fDevCombName[imatch] = Form("%s+%s",fDevName[i].Data(),fDevName[j].Data());
	  fDevi[imatch] = i;  fDevj[imatch] = j;
	  fMatchId1[imatch] = 0; fMatchId2[imatch] = 3;
	  imatch++;
	}
      }

      if(fDevName[i].Contains("B") && fDevName[j].Contains("D")){
	if(fColArr[i][1]==fColArr[j][1] && fColArr[i][2]==fColArr[j][2] &&
	   fSenArr[i][1]==fSenArr[j][1] && fSenArr[i][2]==fSenArr[j][2]){

	  fDevCombName[imatch] = Form("%s+%s",fDevName[i].Data(),fDevName[j].Data());
	  fDevi[imatch] = i;  fDevj[imatch] = j;
	  fMatchId1[imatch] = 2; fMatchId2[imatch] = 3;
	  imatch++;
	}
      }

      if(fDevName[i].Contains("C") && fDevName[j].Contains("D")){
	if(fColArr[i][1]==fColArr[j][0] && fColArr[i][2]==fColArr[j][2] &&
	   fSenArr[i][1]==fSenArr[j][0] && fSenArr[i][2]==fSenArr[j][2]){
	  
	  fDevCombName[imatch] = Form("%s+%s",fDevName[i].Data(),fDevName[j].Data());
	  fDevi[imatch] = i;  fDevj[imatch] = j;
	  fMatchId1[imatch] = 1; fMatchId2[imatch] = 3;
	  imatch++;
	}
      }
      //cout<<"imatch="<<imatch<<endl;  
      
    }
  }
  
  fNDevComb = imatch;
  cout<<"Na61Module::Init(): number of 3Hit combinations(fNDev): "<<fNDev
      <<"   number of 4Hit combinations(fNDevComb): "<<fNDevComb<<endl;  
  
  
}    // Called once per job


//____________________________________________________________________________________
void Na61Module::AnaSensorString(TString stringName,Int_t* ic,Int_t* is)
{

  string Sensor = stringName.Data();
  Int_t StringSize = 0;
  Int_t t1 = 0;
  
  //cout << Sensor << endl;
  StringSize = Sensor.size(); //Number of elements of every type e.g. A_00_0-1_1-1
  //cout << "Size of " << i+1 <<" sensor: " << StringSize << endl;
  
  t1 = 0;
  for(int j = 1; j < StringSize; j++){   
    
    if(Sensor[j] == '_'){  // found splitting char
      
      if(Sensor[j+2] == '-'){
	ic[t1] = (int)Sensor[j+1] - '0';
	is[t1] = -((int)Sensor[j+3] - '0');
	t1++;  
      }else{	      
	ic[t1] = (int)Sensor[j+1] - '0';
	is[t1] = (int)Sensor[j+2] - '0';
	//cout<<"chars[from j]="<<" j="<<j<<" "<<Sensor[j]<<" "<<Sensor[j+1]<<"  "<<Sensor[j+2]<<endl;
	//cout<<"t1="<<t1<<" ic[t1]="<<ic[t1]<<" is[t1]="<<is[t1]<<endl;	  
	t1++;
      }
      
    } // if splitting char found
    
  }    
  
  return;
}

//_________________________________________________________________
int Na61Module::GetIdev(TString devname) {
  for (int i = 0; i < 20; i++)
    if (!fDevNames[i].CompareTo(devname)) return i;

  Error("GetIdev", "idev not found for a given devname, check it out!!!!!");
  return -1;
}
//_________________________________________________________________
int Na61Module::GetMatchIdev(TString matchname) {
  for (int i = 0; i < 20; i++)
    if (!fMatchStr[i].CompareTo(matchname)) return i;

  Error("GetMatchIdev", "idev not found for a given matchname, check it out!!!!!");
  return -1;
}

//_________________________________________________________________
Na61Module::~Na61Module() {
  //
  // default destructor
  // Histograms are assumed to have been generated inside the
  // Book method and ownership set to the histogram list
  //
  if (fHistograms) {
    fHistograms->Delete();
    delete fHistograms;
  }
  if (fRequiredTableList) delete fRequiredTableList;
}

//____________________________________________________________________
void Na61Module::Event() {
  // Event method with no Data objects.
  // This method is called once per event.
  Debug(1, "Event()", "for %s", GetName());
}

//____________________________________________________________________
void Na61Module::Event(UEventNode * /*inev*/, UEventNode * /*outev*/) {
  //
  // Default Event method with in/out Data objects. In normal
  // application to be  called once per event. Method should also be
  // overwritten in the concrete class.
  // This method is called once per event.
  // A dummy routine is supplied in case a module for some strange
  // reason is not called per event.
  //
  Debug(1, "Event(node,node)", "for %s", GetName());
}

//_________________________________________________________________
void Na61Module::Book() {
  // Booking method for Module class.
  // This method in turns calls the 'user overrrideen method
  // DefineHistograms' and adds all defined histograms to the
  // histogram list for the module. If the user wants to have these in
  // a specific root directory this should be setup in
  // DefineHistogram() method. Following Booking histogramming is
  // enabled. It can be turned on/off using the EnableHistograms()
  // DisableHistogram()
  // Once added to this list the histograms should be considerd owned
  // by Module.

  if (fHistograms) {
    Warning("Book", "Histograms already booked for this module");
    return;
  }

  TDirectory *oldGDir = gDirectory;

  fHistograms = new TList();
  TObject *objfirst, *objlast;
  objlast = gDirectory->GetList()->Last();

  // Call overloaded method from the derived class
  DefineHistograms();

  if (objlast)
    objfirst = gDirectory->GetList()->After(objlast);
  else
    objfirst = gDirectory->GetList()->First();

  if (objfirst && objfirst->IsA()->InheritsFrom("TDirectory")) {
    gDirectory->cd(objfirst->GetName());
    objfirst = gDirectory->GetList()->First();
  }

  while (objfirst) {
    fHistograms->Add(objfirst);
    objfirst = gDirectory->GetList()->After(objfirst);
  }

  fHistBooked = true;
  fHistOn = true;

  oldGDir->cd();
}

//_________________________________________________________________
void Na61Module::EnableHistograms() {
  // Enable use of histograms
  if (fHistBooked)
    fHistOn = true;
  else
    Warning("EnableHistograms", "Booking has not been done");
}

//_________________________________________________________________
void Na61Module::DisableHistograms() {
  // Disable use of histograms
  if (fHistBooked)
    fHistOn = false;
  else
    Warning("DisableHistograms", "Booking has not been done");
}

//_________________________________________________________________
void Na61Module::Print(Option_t *option) {
  // Information module. In the final implementation this MUST be
  // overwritten by the derived class. Here a default behaviour is
  // implemented.
  //
  // Options:
  //     D             show all the details
  //     H             show list of histograms
  //     B             DO NOT show basic information
  //
  // This non-const version is only kept for backward compatibility
  // with ROOT pre3. Will soon disappear!
  TString opt(option);
  opt.ToLower();

  if (opt.Contains("b")) cout << "*************************************************" << endl << endl << "  Module" << endl << "   Class :        " << IsA()->GetName() << endl << "   Name  :        " << GetName() << endl << "   Title :        " << GetTitle() << endl << endl << "-------------------------------------------------" << endl;

  if (opt.Contains("d")) {
    // Print all details
    cout << "  Debug level:     " << fDebugLevel << endl << "  Verbosity:       " << fVerbose << endl << "  # of stops:      " << fStopCount << "/" << fStopLimit << endl << "  # of failures:   " << fFailureCount << "/" << fFailureLimit << endl << "  State:           " << flush;
    switch (fState) {
      case kSetup:
        cout << "setup time" << endl;
        break;
      case kInit:
        cout << "Init time" << endl;
        break;
      case kBegin:
        cout << "Begin time" << endl;
        break;
      case kEvent:
        cout << "Event time" << endl;
        break;
      case kEnd:
        cout << "End time" << endl;
        break;
      case kFinish:
        cout << "Finish time" << endl;
        break;
      default:
        cout << "Unknown!!!" << endl;
    }
    cout << "  Status:          " << flush;
    switch (fStatus) {
      case kOk:
        cout << "everything's just dandy" << endl;
        break;
      case kStop:
        cout << "wants pipeline stoped" << endl;
        break;
      case kFailure:
        cout << "failed" << endl;
        break;
      case kAbort:
        cout << "failed miserably" << endl;
        break;
      default:
        cout << "Unknown!!!" << endl;
    }
    cout << "  Histograms:      " << (fHistOn ? "on" : "off") << " and " << (fHistBooked ? "booked" : "not booked") << endl << "-------------------------------------------------" << endl;
  }

  if (opt.Contains("h")) {
    cout << "  Histograms:     " << flush;
    if (!fHistOn || !fHistBooked)
      cout << "not on or booked" << endl;
    else {
      cout << endl;
      TIter next(fHistograms);
      TObject *o = 0;
      while ((o = next())) o->Print(option);
    }
    cout << "-------------------------------------------------" << endl;
  }
}

//_________________________________________________________________
void Na61Module::Print(Option_t *option) const {
  // Information module. In the final implementation this MUST be
  // overwritten by the derived class. Here a default behaviour is
  // implemented.
  //
  // Options:
  //     D             show all the details
  //     H             show list of histograms
  //     B             show basic information (default)
  //     E             show event statistics
  //
  TString opt(option);
  opt.ToLower();
  if (opt.Contains("b")) cout << "*************************************************" << endl << endl << "  Module" << endl << "   Class :         " << IsA()->GetName() << endl << "   Name  :         " << GetName() << endl << "   Title :         " << GetTitle() << endl << endl << "-------------------------------------------------" << endl;

  if (opt.Contains("e")) cout << "  Accumulated CPU time: " << setprecision(5) << setw(10) << fCpuTime << endl;

  if (opt.Contains("d")) {
    // Print all details
    cout << "  Debug level:     " << fDebugLevel << endl << "  Verbosity:       " << fVerbose << endl << "  # of stops:      " << fStopCount << "/" << fStopLimit << endl << "  # of failures:   " << fFailureCount << "/" << fFailureLimit << endl << "  State:           " << flush;
    switch (fState) {
      case kSetup:
        cout << "setup time" << endl;
        break;
      case kInit:
        cout << "Init time" << endl;
        break;
      case kBegin:
        cout << "Begin time" << endl;
        break;
      case kEvent:
        cout << "Event time" << endl;
        break;
      case kEnd:
        cout << "End time" << endl;
        break;
      case kFinish:
        cout << "Finish time" << endl;
        break;
      default:
        cout << "Unknown!!!" << endl;
    }
    cout << "  Status:          " << flush;
    switch (fStatus) {
      case kOk:
        cout << "everything's just dandy" << endl;
        break;
      case kStop:
        cout << "wants pipeline stoped" << endl;
        break;
      case kFailure:
        cout << "failed" << endl;
        break;
      case kAbort:
        cout << "failed miserably" << endl;
        break;
      default:
        cout << "Unknown!!!" << endl;
    }
    cout << "  Histograms:      " << (fHistOn ? "on" : "off") << " and " << (fHistBooked ? "booked" : "not booked") << endl << "-------------------------------------------------" << endl;
  }

  if (opt.Contains("h")) {
    cout << "  Histograms:      " << flush;
    if (!fHistOn || !fHistBooked)
      cout << "not on or booked" << endl;
    else {
      cout << endl;
      TIter next(fHistograms);
      TObject *o = 0;
      while ((o = next())) o->Print(option);
    }
    cout << "-------------------------------------------------" << endl;
  }
}

//_________________________________________________________________
void Na61Module::SetRequiredData(char * /*name*/) {
  //
  // Set the name of a required DataObject. The name will be stored in
  // a list. The intent is that the calling Event method can check against
  // the names (by calling a basSe Module member function which will search
  // the input Event and return a pointer to the object.
  //
  // CheckDataObjects();
  // GetDataObject
  // Note none of this is implemented so far.
  // 2/17/98
  //
}

//_________________________________________________________________
void Na61Module::EventStatisticsStart() {
  // Start the timer for the statistics output
  fTimer.Stop();
  fTimer.Start();
}

//_________________________________________________________________
void Na61Module::EventStatisticsEnd() {
  // Stop the timer for the statistics output
  fTimer.Stop();
  fCpuTime += fTimer.CpuTime();
}

//_________________________________________________________________
void Na61Module::ListEventStatistics() {
  // List the event statistics from the timer. Presupposes that
  // EventStatisticsStart() and EventStatisticsEnd() has been called
  if (Verbose()) cout << "  Acc. CPU time " << setprecision(5) << setw(10) << fCpuTime << " for " << GetName() << " " << GetTitle() << endl;
}

//____________________________________________________________________
void Na61Module::AddRequiredTable(TObject *tablename) {
  // Add a table to the list of requirted tables
  if (fRequiredTableList) fRequiredTableList->Add(tablename);
}

//_________________________________________________________________
void Na61Module::Info(int lvl, const char *location, const char *va_(fmt), ...) const {
  // Information print out, if the verbosity level is greater than or
  // equal to first argument lvl. Remaining arguments as in
  // TObject::Warning and similar
  if (fVerbose < lvl) return;

  va_list ap;
  va_start(ap, va_(fmt));
  DoError(kInfo, location, va_(fmt), ap);
  va_end(ap);
}

//_________________________________________________________________
void Na61Module::Debug(int lvl, const char *location, const char *va_(fmt), ...) const {
  // Debug print out, if the debug level is greater than or
  // equal to first argument lvl. Remaining arguments as in
  // TObject::Warning and similar.  Please note that the line is
  // prefixed with 'Info' and not debug - sorry, but we need our own
  // error handler to handle  that properly.
  if (fDebugLevel < lvl) return;

  va_list ap;
  va_start(ap, va_(fmt));
  DoError(kInfo, location, va_(fmt), ap);
  va_end(ap);
}

//____________________________________________________________________
void Na61Module::Stop(const char *location, const char *va_(fmt), ...) const {
  // Flag this module as having not being able to fullfil it's
  // requirements. That is, the module pipeline should stop after this
  // module. Note, sending this message in Init() will stop the
  // current job altogether, and in Begin() will stop that run loop
  // altogether. If verbose flag is greater then, or equal to
  // kShowStop, then also output message. If this failure makes the
  // total number of failures for this module exceed the stop limit
  // (as set by Na61Module::SetStopLimit), then promote this stop to
  // a Failure (see Na61Module::Failure).

  // Increment failure counter
  fStopCount++;

  // if inside limits, set status to kFailure
  if (fStopCount < fStopLimit) {
    fStatus = kStop;

    // If verbosity is greater than or equal to kShowFailures, then
    // print the message
    if (fVerbose >= kShowStop) {
      va_list ap;
      va_start(ap, va_(fmt));
      DoError(fStatus, location, va_(fmt), ap);
      va_end(ap);
    }
  }
  // If out side of limits, turn this Failure into a disaster.
  else
    Failure("Stop", "Stop limits reached at %d", fStopCount);
}

//____________________________________________________________________
void Na61Module::Failure(const char *location, const char *va_(fmt), ...) const {
  // Flag this module as having failed. If verbose flag is greater
  // then, or equal to kShowFailure, then also output message. If
  // this failure makes the total number of failures for this module
  // exceed the failure limit (as set by Na61Module::SetFailureLimit),
  // then promote this failure to a Abort (see Na61Module::Abort).

  // Increment failure counter
  fFailureCount++;

  // if inside limits, set status to kFailure
  if (fFailureCount < fFailureLimit) {
    fStatus = kFailure;

    // If verbosity is greater than or equal to kShowFailures, then
    // print the message
    if (fVerbose >= kShowFailures) {
      va_list ap;
      va_start(ap, va_(fmt));
      DoError(fStatus, location, va_(fmt), ap);
      va_end(ap);
    }
  }
  // If out side of limits, turn this Failure into a disaster.
  else
    Abort("Failure", "Failure limits reached at %d", fFailureCount);
}

//____________________________________________________________________
void Na61Module::Abort(const char *location, const char *va_(fmt), ...) const {
  // Flag this module as having failed misaerably, and output
  // message. A Na61ModuleContainer should always react to this by
  // stoping the module loop, since this flag means that something is
  // very VERY screwy.
  fStatus = kAbort;
  va_list ap;
  va_start(ap, va_(fmt));
  DoError(fStatus, location, va_(fmt), ap);
  va_end(ap);
}

//____________________________________________________________________
int Na61Module::GetEventNumber(UEvent *ev) const {
  // Return the event number after checking that this is, in fact, a
  // UEvent
  if (ev->IsA() == UEvent::Class())
    return ev->GetEventNumber();
  else
    Error("GetEventNumber", "Input event is not a UEvent!!!");
  return 0;
}
