//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61FrameMergingModule
#include "Na61FrameMergingModule.h"
#endif
#ifndef ROOT_TDirectory
#include "TDirectory.h"
#endif
#ifndef UTIL_UDataTable
#include "UDataTable.h"
#endif
#ifndef UTIL_UG4Hit
#include "UG4Hit.h"
#endif
#ifndef UTIL_USensorPixel
#include "USensorPixel.h"
#endif

#ifndef ROOT_TF1
#include "TF1.h"
#endif
#ifndef ROOT_TMath
#include "TMath.h"
#endif
#include "TSystem.h"
#include "TDirectory.h"
#include "Math/QuantFuncMathCore.h"
#include <time.h>
#include "iostream"
#include "UG4RecoTrack.h"

//____________________________________________________________________
// ClassImp(Na61FrameMergingModule);

//____________________________________________________________________
Na61FrameMergingModule::Na61FrameMergingModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
  fEvents = 0;
  fM26CycleTime = 115.2;  // Frame redout time in um
  fCounter_f1_1 = 0;
  fCounter_f1_2 = 0;
  fCounter_f2_1 = 0;
  fCounter_f2_2 = 0;
  fRowOverlap = 0;
  fXRowOffset = 0;
}
//____________________________________________________________________
Na61FrameMergingModule::Na61FrameMergingModule(const char* name, const char* title) : Na61Module(name, title) {
  // Named Constructor
  SetState(kSetup);
  fEvents = 0;
  fM26CycleTime = 115.2;  // Frame redout time in um
  fCounter_f1_1 = 0;
  fCounter_f1_2 = 0;
  fCounter_f2_1 = 0;
  fCounter_f2_2 = 0;
  fRowOverlap = 0;
  fXRowOffset = 0;
}
//____________________________________________________________________
void Na61FrameMergingModule::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  TDirectory* histDir;

  if (fJura)
    histDir = gDirectory->mkdir("FrameMergingModule_Jura");
  else
    histDir = gDirectory->mkdir("FrameMergingModule_Saleve");

  histDir->cd();

  // Define histograms. They are:
  // <fill in here>

  int binsy = findy;
  int binsx = findx;
  double ymax = binsy;
  double xmax = binsx;

  for (int i = 0; i < 8; i++) {
    fhPixels_f0[i] = new TH1F(Form("hPixels_f0_%s", fSensorNames[i].Data()), "", 100, 0, 100);
    fhPixels_f1[i] = new TH1F(Form("hPixels_f1_%s", fSensorNames[i].Data()), "", 100, 0, 100);
    fhPixels_f2[i] = new TH1F(Form("hPixels_f2_%s", fSensorNames[i].Data()), "", 100, 0, 100);
    fhPixels_f3[i] = new TH1F(Form("hPixels_f3_%s", fSensorNames[i].Data()), "", 100, 0, 100);
    fhPixelsMerged[i] = new TH1F(Form("hPixelsMerged_%s", fSensorNames[i].Data()), "", 500, 0, 500);
    // histos for noisy pixels extraction
    fhPixelMap_f34[i] = new TH2F(Form("hPixelMap_f34_%s", fSensorNames[i].Data()), "", binsx, 0.5, xmax + 0.5, binsy, 0.5, ymax + 0.5);
    fhNoisyPixelMap[i] = new TH2F(Form("hNoisyPixelMap_%s", fSensorNames[i].Data()), "", binsx, 0.5, xmax + 0.5, binsy, 0.5, ymax + 0.5);
    fhPixelFrequency[i] = new TH1F(Form("hPixelFrequency_%s", fSensorNames[i].Data()), "", 1000, 0, 100);
    fhPixelsVsTimer[i] = new TH1F(Form("hPixelVsTimer_%s", fSensorNames[i].Data()), "", 100, 0, 700);
    fhPixelsVsTimer_ref[i] = new TH1F(Form("hPixelVsTimer_ref_%s", fSensorNames[i].Data()), "", 100, 0, 700);
  }

  for (int i = 0; i < 3; i++) {
    fhTimer_f0[i] = new TH1F(Form("hTimer_f0_fpga%d", i + 1), "", 700, 0, 700);
    fhTimer_f1[i] = new TH1F(Form("hTimer_f1_fpga%d", i + 1), "", 700, 0, 700);
    fhTimer_f2[i] = new TH1F(Form("hTimer_f2_fpga%d", i + 1), "", 700, 0, 700);
    fhTimer_f3[i] = new TH1F(Form("hTimer_f3_fpga%d", i + 1), "", 700, 0, 700);
    fhTimer_f4[i] = new TH1F(Form("hTimer_f4_fpga%d", i + 1), "", 700, 0, 700);

    fhTimer_acce_f0[i] = new TH1F(Form("hTimer_acce_f0_fpga%d", i + 1), "", 700, 0, 700);
    fhTimer_acce_f1[i] = new TH1F(Form("hTimer_acce_f1_fpga%d", i + 1), "", 700, 0, 700);
    fhTimer_acce_f2[i] = new TH1F(Form("hTimer_acce_f2_fpga%d", i + 1), "", 700, 0, 700);
    fhTimer_acce_f3[i] = new TH1F(Form("hTimer_acce_f3_fpga%d", i + 1), "", 700, 0, 700);
    fhTimer_acce_f4[i] = new TH1F(Form("hTimer_acce_f4_fpga%d", i + 1), "", 700, 0, 700);
  }

  fhTimer_1vs2 = new TH2F("hTimer_1vs2", "", 700, 0, 700, 700, 0, 700);
  fhTimer_1vs3 = new TH2F("hTimer_1vs3", "", 700, 0, 700, 700, 0, 700);
  fhTimer_2vs3 = new TH2F("hTimer_2vs3", "", 700, 0, 700, 700, 0, 700);

  ReadNoisyPixelInfo();
  cout << " fNoisyPixelsDefined=" << fNoisyPixelsDefined << endl;

  histDir->cd();
  gDirectory->cd("..");
}

//____________________________________________________________________
void Na61FrameMergingModule::Init() {
  // Job-level initialisation
  SetState(kInit);

  fNVds = 4;
  fVdsZ[0] = 0;
  fVdsZ[1] = fVdsZ[0] + 50;
  fVdsZ[2] = fVdsZ[0] + 100;
  fVdsZ[3] = fVdsZ[0] + 150.;
}

//____________________________________________________________________
void Na61FrameMergingModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);
}

//____________________________________________________________________
void Na61FrameMergingModule::ReadNoisyPixelInfo() {
  // If file with noisy pixel maps exist the noisy pixel analysis and recreatinion of this
  // file will be skept in Finish.
  fNoisyPixelsDefined = false;
  TFile* file = 0;
  // TString fHome_path = "/afs/cern.ch/user/p/pstaszel/public/branches/pstaszel/vdreco";
  if (fJura) {
    cout << " Na61FrameMergingModule::ReadNoisyPixelInfo: opening file " << Form("%s/juraNoiseCuts-%06d.root", fHomePath.Data(), fRunId) << endl;
    file = new TFile(Form("%s/juraNoiseCuts-%06d.root", fHomePath.Data(), fRunId));
  } else {
    cout << " Na61FrameMerlgingModule::ReadNoisyPixelInfo: opening file " << Form("%s/saleveNoiseCuts-%06d.root", fHomePath.Data(), fRunId) << endl;
    file = new TFile(Form("%s/saleveNoiseCuts-%06d.root", fHomePath.Data(), fRunId));
  }

  if (!file->IsOpen()) {
    cout << " Na61FrameMergingModule::ReadNoisyPixelInfo: can not open noisy pixel file. Check it out!!!" << endl;
    return;
  }

  // take histos with maps of noisy pixel

  fNoisyPixelsDefined = true;
  for (int imap = 0; imap < 8; imap++) {
    TH2F* hh = (TH2F*)file->Get(Form("hNoisyPixelMap_%s", fSensorNames[imap].Data()));
    if (!hh) {
      fNoisyPixelsDefined = false;
    } else {
      int ii = 0;
      for (int i = 1; i < hh->GetNbinsX() + 1; i++) {
        for (int j = 1; j < hh->GetNbinsY() + 1; j++) {
          if (hh->GetBinContent(i, j)) ii++;
          fhNoisyPixelMap[imap]->SetBinContent(i, j, hh->GetBinContent(i, j));
        }
      }
      fhNoisyPixelMap[imap]->SetEntries(ii);
    }
  }

  TH1F* h = (TH1F*)file->Get(Form("hCutValue"));
  double cutValue = h->GetBinContent(1);

  // require all maps to be defined, if not fNoisyPixelsDefined = false;
  file->Close();

  cout << "cutValue=" << cutValue << " " << fPixelNoiseCut << " " << fNoisyPixelsDefined << endl;
  if ((cutValue != fPixelNoiseCut) && fNoisyPixelsDefined) {
    cout << "previous noisy pixel cut value was" << cutValue << " but current cut value is " << fPixelNoiseCut << endl;
    cout << "Redo the noisy pixel maps? (1=yes, 0=no)" << endl;
    int iii;
    cin >> iii;
    if (iii == 1) fNoisyPixelsDefined = false;
  }

  if (fNoisyPixelsDefined) {
    cout << " Na61FrameMergingModule::Init: cut on noisy pixels activated with cut value = " << cutValue << " %" << endl;
    for (int i = 0; i < 8; i++) cout << "Number of noisy pixels in sensor " << fSensorNames[i].Data() << ": " << fhNoisyPixelMap[i]->GetEntries() << endl;

  } else {
    cout << " Na61FrameMergingModule::Init: will create new noisy pixel maps for cut value = " << fPixelNoiseCut << " %."
         << " Cut on noisy pixels is disabled for this run." << endl;
  }
}

//____________________________________________________________________
void Na61FrameMergingModule::Event(UEventNode* inNode, UEventNode* /*outNode*/)

{
  // Per event method
  SetState(kEvent);

  UDataTable* tabs_f0[8];
  UDataTable* tabs_f1[8];
  UDataTable* tabs_f2[8];
  UDataTable* tabs_f3[8];
  UDataTable* tabs_f4[8];

  for (int i = 0; i < 8; i++) {
    tabs_f0[i] = inNode->GetDataTable(Form("Pixels %s f0", fSensorNames[i].Data()));
    tabs_f1[i] = inNode->GetDataTable(Form("Pixels %s f1", fSensorNames[i].Data()));
    tabs_f2[i] = inNode->GetDataTable(Form("Pixels %s f2", fSensorNames[i].Data()));
    tabs_f3[i] = inNode->GetDataTable(Form("Pixels %s f3", fSensorNames[i].Data()));
    tabs_f4[i] = inNode->GetDataTable(Form("Pixels %s f4", fSensorNames[i].Data()));
     cout<<"sensor(i)="<<i<<"  fJura = "<<fJura<<"  pixels: "<<tabs_f0[i]->GetEntries()<<" "<<tabs_f1[i]->GetEntries()<<" "<<tabs_f2[i]->GetEntries()<<endl;
  }

  for (int i = 0; i < 8; i++) {
    if (tabs_f0[i]) {
      fhPixels_f0[i]->Fill(tabs_f0[i]->GetEntries());
      fhPixelsVsTimer[i]->Fill(fTimers[0][0], tabs_f0[i]->GetEntries());
      fhPixelsVsTimer_ref[i]->Fill(fTimers[0][0]);
    } else {
     // cout << " tabs_f0[i]: " << tabs_f0[i] << endl;
    }
    if (tabs_f1[i]) {
      fhPixels_f1[i]->Fill(tabs_f1[i]->GetEntries());
      fhPixelsVsTimer[i]->Fill(fTimers[0][1], tabs_f1[i]->GetEntries());
      fhPixelsVsTimer_ref[i]->Fill(fTimers[0][1]);
    } else {
     // cout << " tabs_f1[i]: " << tabs_f1[i] << endl;
    }

    if (tabs_f2[i]) {
      fhPixels_f2[i]->Fill(tabs_f2[i]->GetEntries());
      fhPixelsVsTimer[i]->Fill(fTimers[0][2], tabs_f2[i]->GetEntries());
      fhPixelsVsTimer_ref[i]->Fill(fTimers[0][2]);
    } else {
     // cout << " tabs_f2[i]: " << tabs_f2[i] << endl;
    }

    if (tabs_f3[i]) {
      fhPixels_f3[i]->Fill(tabs_f3[i]->GetEntries());
      fhPixelsVsTimer[i]->Fill(fTimers[0][3], tabs_f3[i]->GetEntries());
      fhPixelsVsTimer_ref[i]->Fill(fTimers[0][3]);
    } else {
     // cout << " tabs_f3[i]: " << tabs_f3[i] << endl;
    }

    if (tabs_f4[i]) {
      // fhPixels_f4[i] -> Fill(tabs_f4[i]->GetEntries());
      fhPixelsVsTimer[i]->Fill(fTimers[0][4], tabs_f4[i]->GetEntries());
      fhPixelsVsTimer_ref[i]->Fill(fTimers[0][4]);
    } else {
     // cout << " tabs_f4[i]: " << tabs_f4[i] << endl;
    }
  }

  // fill pixel frequency distribution
  fEvents++;
  for (int itab = 0; itab < 8; itab++) {
    if (!tabs_f3[itab]) continue;
    for (int i = 0; i < tabs_f3[itab]->GetEntries(); i++) {
      USensorPixel* pix = (USensorPixel*)tabs_f3[itab]->At(i);
      if (!pix) Error("Event", "something is wrong, check it out");
      fhPixelMap_f34[itab]->Fill(pix->GetLine(), pix->GetColumn());
    }
    if (!tabs_f4[itab]) continue;
    for (int i = 0; i < tabs_f4[itab]->GetEntries(); i++) {
      USensorPixel* pix = (USensorPixel*)tabs_f4[itab]->At(i);
      if (!pix) Error("Event", "something is wrong, check it out");
      fhPixelMap_f34[itab]->Fill(pix->GetLine(), pix->GetColumn());
    }
  }

  for (int i = 0; i < 3; i++) {
    fhTimer_f0[i]->Fill(fTimers[i][0]);
    fhTimer_f1[i]->Fill(fTimers[i][1]);

    fhTimer_f2[i]->Fill(fTimers[i][2]);
    fhTimer_f3[i]->Fill(fTimers[i][3]);
    fhTimer_f4[i]->Fill(fTimers[i][4]);
  }

  for (int i = 0; i < 4; i++) {
    fhTimer_1vs2->Fill(fTimers[0][i], fTimers[1][i]);
    fhTimer_1vs3->Fill(fTimers[0][i], fTimers[2][i]);
    fhTimer_2vs3->Fill(fTimers[1][i], fTimers[2][i]);
  }

  // marge pixels
  UDataTable* tabs_merge[8];
  for (int i = 0; i < 8; i++) {
    if (fJura)
      tabs_merge[i] = new UDataTable(Form("Jura Pixels %s merged", fSensorNames[i].Data()));
    else
      tabs_merge[i] = new UDataTable(Form("Saleve Pixels %s merged", fSensorNames[i].Data()));
  }

  ResetFreqArray();  // needed for not writing double pixels

  if (fNewMerge==1)
    MergeFrames_new(tabs_merge, tabs_f1, tabs_f2);
  else if (fNewMerge==2)
    MergeFrames_SIM(tabs_merge, tabs_f0);
  else
    MergeFrames(tabs_merge, tabs_f0, tabs_f1, tabs_f2, tabs_f3, tabs_f4);

  for (int i = 0; i < 8; i++) {
    if (tabs_merge[i]) {
      fhPixelsMerged[i]->Fill(tabs_merge[i]->GetEntries());
      inNode->AddDataTable(tabs_merge[i]);
      // cout<<"Frame: tabs_merge: i="<<i<<"  pixels: "<<tabs_merge[i]->GetEntries()<<endl;
      // outNode->AddDataTable(tabs_merge[i]);
    }
  }
}

//_________________________________________________________________________________
void Na61FrameMergingModule::MergeFrames_new(UDataTable** tabs_merge, UDataTable** tabs_f1, UDataTable** tabs_f2) {
  // this method intend to merge olny fraction of secons and thyrd frame.
  // Should be used for XeLa and pPb data taken in 2017.
  // Form PbPb taken in 2016 tast before using first.

  double timer_f1 = fTimers[0][1] - fM26CycleTime;
  double timer_f2 = fTimers[0][2] - 2 * fM26CycleTime;

  double fraction_f1 = timer_f1 / fM26CycleTime;
  // double fraction_f2 = 1. - fraction_f1;

  // take last fraction_f1 rows/lines from f1 and first fraction_f2 rows/lines from f2.

  double x_row = (1.0 - fraction_f1) * 576;  // check if we count from 0 or from 1
  x_row = x_row + fXRowOffset;

  // cout<<timer_f1<<" "<<timer_f2<<" "<< fraction_f1<<"  "<<x_row<<" "<<TMath::Nint(x_row)<<" RowOverlap="<<fRowOverlap<<endl;

  AddPixels_new(1, TMath::Nint(x_row), tabs_merge, tabs_f1);
  AddPixels_new(2, TMath::Nint(x_row), tabs_merge, tabs_f2);

  fhTimer_acce_f1[0]->Fill(timer_f1);
  fhTimer_acce_f2[0]->Fill(timer_f2);
}

//_________________________________________________________________________________
void Na61FrameMergingModule::MergeFrames(UDataTable** tabs_merge, UDataTable** tabs_f0, UDataTable** tabs_f1, UDataTable** tabs_f2, UDataTable** tabs_f3, UDataTable** tabs_f4) {
  if ((fTimers[0][0] > fMinTimer) && (fTimers[0][0] < fMaxTimer)) {
    AddPixels(tabs_merge, tabs_f0);
    fhTimer_acce_f0[0]->Fill(fTimers[0][0]);
    fhTimer_acce_f0[1]->Fill(fTimers[1][0]);
    fhTimer_acce_f0[2]->Fill(fTimers[2][0]);
  }

  if ((fTimers[0][1] > fMinTimer) && (fTimers[0][1] < fMaxTimer)) {
    AddPixels(tabs_merge, tabs_f1);
    fhTimer_acce_f1[0]->Fill(fTimers[0][1]);
    fhTimer_acce_f1[1]->Fill(fTimers[1][1]);
    fhTimer_acce_f1[2]->Fill(fTimers[2][1]);
  }

  if ((fTimers[0][2] > fMinTimer) && (fTimers[0][2] < fMaxTimer)) {
    AddPixels(tabs_merge, tabs_f2);
    fhTimer_acce_f2[0]->Fill(fTimers[0][2]);
    fhTimer_acce_f2[1]->Fill(fTimers[1][2]);
    fhTimer_acce_f2[2]->Fill(fTimers[2][2]);
  }

  if ((fTimers[0][3] > fMinTimer) && (fTimers[0][3] < fMaxTimer)) {
    AddPixels(tabs_merge, tabs_f3);
    fhTimer_acce_f3[0]->Fill(fTimers[0][3]);
    fhTimer_acce_f3[1]->Fill(fTimers[1][3]);
    fhTimer_acce_f3[2]->Fill(fTimers[2][3]);
  }

  if ((fTimers[0][4] > fMinTimer) && (fTimers[0][4] < fMaxTimer)) {
    AddPixels(tabs_merge, tabs_f4);
    fhTimer_acce_f4[0]->Fill(fTimers[0][4]);
    fhTimer_acce_f4[1]->Fill(fTimers[1][4]);
    fhTimer_acce_f4[2]->Fill(fTimers[2][4]);
  }
}

//_________________________________________________________________________________
void Na61FrameMergingModule::MergeFrames_SIM(UDataTable** tabs_merge, UDataTable** tabs_f0) {
   //AddPixels_SIM(tabs_merge, tabs_f0);
    AddPixels_SIM(0, tabs_merge, tabs_f0);

    fhTimer_acce_f0[0]->Fill(fTimers[0][0]);
    fhTimer_acce_f0[1]->Fill(fTimers[1][0]);
    fhTimer_acce_f0[2]->Fill(fTimers[2][0]);
}

//_________________________________________________________________________________
void Na61FrameMergingModule::ResetFreqArray() {
  // fmapArr[8]
  for (int is = 0; is < 8; is++)
    for (int i = 0; i < findx + 1; i++)
      for (int j = 0; j < findy + 1; j++) faxayFreq[is][i][j] = false;
}

//_________________________________________________________________________________
void Na61FrameMergingModule::AddPixels_new(int frameId, int row_x, UDataTable** tabs_merge, UDataTable** tabs) {
  for (int itab = 0; itab < 8; itab++) {
    if (!tabs[itab]) continue;

    for (int i = 0; i < tabs[itab]->GetEntries(); i++) {
      USensorPixel* pix = (USensorPixel*)tabs[itab]->At(i);

      // here put condition on noisy pixels
      int ii = pix->GetLine();
      int jj = pix->GetColumn();

      /*
      if(TMath::Abs(row_x-288)<10){
        if(frameId==1){
          if(ii<row_x)fCounter_f1_1++;
          else fCounter_f1_2++;
        }
        if(frameId==2){
          if(ii<row_x)fCounter_f2_1++;
          else fCounter_f2_2++;
        }
      }
      */

      if (faxayFreq[itab][ii][jj]) continue;  // to not store double pixels at the same position

      if ((frameId == 1) && (ii < (row_x - fRowOverlap))) continue;  // take 5 more rows (for overlap)
      if ((frameId == 2) && (ii > (row_x + fRowOverlap))) continue;  // take 5 more rows (for overlap)

      faxayFreq[itab][ii][jj] = true;

      int content = 0;
      // cout<<"tu1 "<<fhNoisyPixelMap[itab]<<" "<<itab<<" "<<ii<<" "<<jj<<endl;
      // cout<<"tu3 entries: "<<fhNoisyPixelMap[itab]->GetEntries()<<" "<<fhNoisyPixelMap[itab]->GetBinContent(20,20)<<endl;

      if (fNoisyPixelsDefined) content = fhNoisyPixelMap[itab]->GetBinContent(ii, jj);
      // cout<<"tu2"<<endl;

      // note, that fNoisyPixelsDefined=false disables the cut

      if (!content) tabs_merge[itab]->Add(new USensorPixel(pix->GetLine(), pix->GetColumn()));
    }
  }
}
//_________________________________________________________________________________
void Na61FrameMergingModule::AddPixels_SIM(int /* frameId */, UDataTable** tabs_merge, UDataTable** tabs) {
  for (int itab = 0; itab < 8; itab++) {
    if (!tabs[itab]) continue;

    for (int i = 0; i < tabs[itab]->GetEntries(); i++) {
      USensorPixel* pix = (USensorPixel*)tabs[itab]->At(i);

      // here put condition on noisy pixels
      int ii = pix->GetLine();
      int jj = pix->GetColumn();

      if (faxayFreq[itab][ii][jj]) continue;  // to not store double pixels at the same position

      //if ((frameId == 1) && (ii < (row_x - fRowOverlap))) continue;  // take 5 more rows (for overlap)
      //if ((frameId == 2) && (ii > (row_x + fRowOverlap))) continue;  // take 5 more rows (for overlap)

      faxayFreq[itab][ii][jj] = true;

      int content = 0;
      if (fNoisyPixelsDefined) content = fhNoisyPixelMap[itab]->GetBinContent(ii, jj);

      if (!content) tabs_merge[itab]->Add(new USensorPixel(pix->GetLine(), pix->GetColumn()));
    }
  }
}


//_________________________________________________________________________________
void Na61FrameMergingModule::AddPixels(UDataTable** tabs_merge, UDataTable** tabs) {
  for (int itab = 0; itab < 8; itab++) {
    if (!tabs[itab]) continue;

    for (int i = 0; i < tabs[itab]->GetEntries(); i++) {
      USensorPixel* pix = (USensorPixel*)tabs[itab]->At(i);

      // here put condition on noisy pixels
      int ii = pix->GetLine();
      int jj = pix->GetColumn();

      if (faxayFreq[itab][ii][jj]) continue;  // to not store double pixels at the same position

      faxayFreq[itab][ii][jj] = true;

      int content = 0;
      // cout<<"tu1 "<<fhNoisyPixelMap[itab]<<" "<<itab<<" "<<ii<<" "<<jj<<endl;
      // cout<<"tu3 entries: "<<fhNoisyPixelMap[itab]->GetEntries()<<" "<<fhNoisyPixelMap[itab]->GetBinContent(20,20)<<endl;

      if (fNoisyPixelsDefined) content = fhNoisyPixelMap[itab]->GetBinContent(ii, jj);
      // cout<<"tu2"<<endl;

      // note, that fNoisyPixelsDefined=false disables the cut

      if (!content) tabs_merge[itab]->Add(new USensorPixel(pix->GetLine(), pix->GetColumn()));
    }
  }
}

//_____________________________________________________________
void Na61FrameMergingModule::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61FrameMergingModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);

  cout << "fCouter_f1_1 = " << fCounter_f1_1 << "  fCouter_f1_2 = " << fCounter_f1_2 << endl;
  cout << "fCouter_f2_1 = " << fCounter_f2_1 << "  fCouter_f2_2 = " << fCounter_f2_2 << endl;

  for (int i = 0; i < 8; i++) {
    // cout<<"i="<<i<<" "<<fhPixelsVsTimer[i]<<" "<<fhPixelsVsTimer_ref[i]<<endl;
    fhPixelsVsTimer[i]->Divide(fhPixelsVsTimer_ref[i]);
  }

  if (fNoisyPixelsDefined) return;

  // TString home_path = "/afs/cern.ch/user/p/pstaszel/public/branches/pstaszel/vdreco";

  TFile pixeldbFile;

  if (fJura) {
    cout << "Writing noise data to " << Form("%s/juraNoiseCuts-%06d.root", fHomePath.Data(), fRunId) << endl;
    pixeldbFile.Open(Form("%s/juraNoiseCuts-%06d.root", fHomePath.Data(), fRunId), "recreate");
  } else {
    cout << "Writing noise data to " << Form("%s/saleveNoiseCuts-%06d.root", fHomePath.Data(), fRunId) << endl;
    pixeldbFile.Open(Form("%s/saleveNoiseCuts-%06d.root", fHomePath.Data(), fRunId), "recreate");
  }

  int binsy = findy;
  int binsx = findx;
  double ymax = binsy;
  double xmax = binsx;

  TH1F* hcutValue = new TH1F("hCutValue", "", 1, 0, 1);
  hcutValue->SetBinContent(1, fPixelNoiseCut);
  hcutValue->Write();

  TH2F* hNoisyPixelMap[8];

  for (int imap = 0; imap < 8; imap++) {
    hNoisyPixelMap[imap] = new TH2F(Form("hNoisyPixelMap_%s", fSensorNames[imap].Data()), "", binsx, 0.5, xmax + 0.5, binsy, 0.5, ymax + 0.5);

    cout << "Finish: imap=" << imap << " " << fhPixelMap_f34[imap] << endl;
    for (int i = 1; i < fhPixelMap_f34[imap]->GetNbinsX() + 1; i++) {
      for (int j = 1; j < fhPixelMap_f34[imap]->GetNbinsY() + 1; j++) {
        double content = fhPixelMap_f34[imap]->GetBinContent(i, j);
        fhPixelFrequency[imap]->Fill(0.5 * content * 100. / (double)fEvents);  // 0.5 because 2 frames added
        if ((0.5 * content * 100. / (double)fEvents) > fPixelNoiseCut) {
          hNoisyPixelMap[imap]->SetBinContent(i, j, 1);
        }
      }
    }
  }

  for (int imap = 0; imap < 8; imap++) hNoisyPixelMap[imap]->Write();
  pixeldbFile.Close();
}

//____________________________________________________________________
void Na61FrameMergingModule::Print(Option_t* option) const {
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>

  TString opt(option);
  opt.ToLower();

  Na61Module::Print(option);
  if (opt.Contains("d")) cout << endl << "  Original author: Paweł Staszel" << endl << "  Last Modifications: " << endl << "    $Author: Staszel $" << endl << "    $Date: 2016/10/25$" << endl << "    $Revision: 1.0 $ " << endl << endl << "-------------------------------------------------" << endl;
}
