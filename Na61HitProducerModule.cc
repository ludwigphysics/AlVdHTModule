//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61HitProducerModule
#include "Na61HitProducerModule.h"
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

#ifndef NA61_Na61VdParametersManager
#include "Na61VdParametersManager.h"
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
// ClassImp(Na61HitProducerModule);

//____________________________________________________________________
Na61HitProducerModule::Na61HitProducerModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
  fProductionMode = 1;
  fNorm = 0;
}
//____________________________________________________________________
Na61HitProducerModule::Na61HitProducerModule(const char* name, const char* title) : Na61Module(name, title) {
  // Named Constructor
  SetState(kSetup);
  fProductionMode = 1;

  fPixelInfoArray = new TObjArray();
  fPixelInfoArray->SetOwner(true);

  // fClusterArray = new TObjArray();
  // fClusterArray->SetOwner(true);
  fNorm = 0;
}
//____________________________________________________________________
void Na61HitProducerModule::DefineHistograms() {
  if (GetState() != kInit) {
    Stop("DefineHistograms", "Must be called after Init");
    return;
  }

  TDirectory* histDir;

  if (fJura)
    histDir = gDirectory->mkdir("HitProducerModule_Jura");
  else
    histDir = gDirectory->mkdir("HitProducerModule_Saleve");

  histDir->cd();

  int binsy = findy;
  int binsx = findx;
  double ymax = binsy;
  double xmax = binsx;

  for (int i = 0; i < 8; i++) {
    fhLargeClusterPosMap[i] = new TH2F(Form("hLargeClusterPosMap_%s", fSensorNames[i].Data()), "", 106, -5.3, 5.3, 212, -10.6, 10.6);
  }

  for (int i = 0; i < 8; i++) {
    fhPixelMap[i] = new TH2F(Form("hPixelMap_%s", fSensorNames[i].Data()), "", binsx, 0, xmax, binsy, 0, ymax);
    fhClusterMap[i] = new TH2F(Form("hClusterMap_%s", fSensorNames[i].Data()), "", binsx, 0, xmax, binsy, 0, ymax);
    fhClusterPosMap[i] = new TH2F(Form("hClusterPosMap_%s", fSensorNames[i].Data()), "", binsx, -5.3, 5.3, binsy, -10.6, 10.6);
    fhClusters[i] = new TH1F(Form("hClusters_%s", fSensorNames[i].Data()), "", 300, 0, 300);
    fhClustersLt1[i] = new TH1F(Form("hClustersLt1_%s", fSensorNames[i].Data()), "", 300, 0, 300);
    fhClusterSize[i] = new TH1F(Form("hClusterSize_%s", fSensorNames[i].Data()), "", 350, 0, 350);
  }

  fhHitsXVsZ = new TH2F("HitsXVsZ", "", 5000, -1, 151, 1000, -20, 20.);
  fhHitsYVsZ = new TH2F("HitsYVsZ", "", 5000, -1, 151, 1000, -22, 22.);

  fhClusterVerticalDist_Vds1 = new TH1F("hClusterVerticalDist_Vds1", "", 220, -24, 24);
  fhClusterVerticalDist_Vds2 = new TH1F("hClusterVerticalDist_Vds2", "", 220, -24, 24);
  fhClusterVerticalDist_Vds3 = new TH1F("hClusterVerticalDist_Vds3", "", 220, -24, 24);
  fhClusterVerticalDist_Vds4a = new TH1F("hClusterVerticalDist_Vds4a", "", 220, -24, 24);
  fhClusterVerticalDist_Vds4b = new TH1F("hClusterVerticalDist_Vds4b", "", 220, -24, 24);

  fhSignalToNoise = new TH1F("hSignalToNoise", "", 8, 0, 8);
  // Define histograms. They are:
  // <fill in here>

  gDirectory->cd("..");
}

//____________________________________________________________________
void Na61HitProducerModule::SetRunInfo(TFile* file) {
  TDirectory* dir = file->mkdir("runInfo");
  dir->cd();

  TH1F* hVolumesX = new TH1F("volumesX", "", 8, 0, 8);
  TH1F* hVolumesY = new TH1F("volumesY", "", 8, 0, 8);
  TH1F* hVolumesZ = new TH1F("volumesZ", "", 8, 0, 8);

  TH1F* hRotX = new TH1F("rotX", "", 8, 0, 8);
  TH1F* hRotY = new TH1F("rotY", "", 8, 0, 8);
  TH1F* hRotZ = new TH1F("rotZ", "", 8, 0, 8);

  for (int i = 0; i < 8; i++) {
    hVolumesX->SetBinContent(i + 1, fParams->GetVolumeX(i));
    hVolumesY->SetBinContent(i + 1, fParams->GetVolumeY(i));
    hVolumesZ->SetBinContent(i + 1, fParams->GetVolumeZ(i));
    hRotX->SetBinContent(i + 1, fParams->GetRotX(i));
    hRotY->SetBinContent(i + 1, fParams->GetRotY(i));
    hRotZ->SetBinContent(i + 1, fParams->GetRotZ(i));
  }
  hVolumesX->Write();
  hVolumesY->Write();
  hVolumesZ->Write();
  hRotX->Write();
  hRotY->Write();
  hRotZ->Write();

  dir->cd("/");
}

//____________________________________________________________________
void Na61HitProducerModule::Init() {
  // Job-level initialisation
  SetState(kInit);

  if (fJura)
    fParams = (Na61VdParametersManager::Instance())->GetJuraArmParams();
  else
    fParams = (Na61VdParametersManager::Instance())->GetSaleveArmParams();

  cout << "Setting controll parameters in ArmParameters: fSensorId=" << fSensorId << " fRunId=" << fRunId << endl;
  fParams->SetDrotX(fSensorId, fDrotX);
  fParams->SetDrotY(fSensorId, fDrotY);
  fParams->SetDrotZ(fSensorId, fDrotZ);
  fParams->SetVolumeDz(fSensorId, fDz);

  fParams->SetRunId(fRunId);
  fParams->Init();
}

//____________________________________________________________________
void Na61HitProducerModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);
}

//____________________________________________________________________
void Na61HitProducerModule::Event(UEventNode* inNode, UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);

  for (int i = 0; i < 8; i++) {
    UDataTable* pixtab;
    if (fJura)
      pixtab = outNode->GetDataTable(Form("Jura Pixels %s merged", fSensorNames[i].Data()));
    else
      pixtab = outNode->GetDataTable(Form("Saleve Pixels %s merged", fSensorNames[i].Data()));

    if (!pixtab) {
      if (fJura)
        pixtab = inNode->GetDataTable(Form("Jura Pixels %s merged", fSensorNames[i].Data()));
      else
        pixtab = inNode->GetDataTable(Form("Saleve Pixels %s merged", fSensorNames[i].Data()));
    }

    if (pixtab) {
      UDataTable* hitstab = new UDataTable(Form("Hits %s", fSensorNames[i].Data()));
      MakeClusters(fSensorNames[i], pixtab, hitstab);
      outNode->AddDataTable(hitstab);
      // cout<<"sensor(i)="<<i<<"  fJura = "<<fJura<<"  hits: "<<hitstab->GetEntries()<<endl;

      if (!fProductionMode) {
        FillHSaxayWithPixels(pixtab, fhPixelMap[i]);
      }

      FillHSaxayWithClusters(hitstab, fhClusterMap[i], fhClusterPosMap[i], fhLargeClusterPosMap[i], i);
      double sigToNoise = FillClusterDistr(hitstab, fhClusters[i], fhClustersLt1[i], fhClusterSize[i]);
      fhSignalToNoise->Fill(i, sigToNoise);
      fNorm++;
      LocalToGlobal(i, hitstab);  // move to global coordinate system
    }
  }

  FillHorizontalDistributions(outNode);

  for (int i = 0; i < 8; i++) {
    UDataTable* hits = outNode->GetDataTable(Form("Hits %s", fSensorNames[i].Data()));
    fClusters[i] = hits->GetEntries();
  }
}

//_____________________________________________________________________
void Na61HitProducerModule::FillHorizontalDistributions(UEventNode* outNode) {
  int sizeTh = 1;
  UDataTable* hitsVds1_0 = outNode->GetDataTable("Hits Vds1_0");
  fRealClusters[0] = 0;
  for (int i = 0; i < hitsVds1_0->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hitsVds1_0->At(i);
    if (hit->GetClusterSize() > sizeTh) {
      fRealClusters[0]++;
      if (fhClusterVerticalDist_Vds1) fhClusterVerticalDist_Vds1->Fill(hit->GetY());
    }
  }

  UDataTable* hitsVds2_0 = outNode->GetDataTable("Hits Vds2_0");
  fRealClusters[1] = 0;
  for (int i = 0; i < hitsVds2_0->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hitsVds2_0->At(i);
    if (hit->GetClusterSize() > sizeTh) {
      fRealClusters[1]++;
      if (fhClusterVerticalDist_Vds2) fhClusterVerticalDist_Vds2->Fill(hit->GetY());
    }
  }

  UDataTable* hitsVds3_0 = outNode->GetDataTable("Hits Vds3_0");
  fRealClusters[2] = 0;
  for (int i = 0; i < hitsVds3_0->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hitsVds3_0->At(i);
    if (hit->GetClusterSize() > sizeTh) {
      fRealClusters[2]++;
      if (fhClusterVerticalDist_Vds3) fhClusterVerticalDist_Vds3->Fill(hit->GetY());
    }
  }

  UDataTable* hitsVds3_1 = outNode->GetDataTable("Hits Vds3_1");
  fRealClusters[3] = 0;
  for (int i = 0; i < hitsVds3_1->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hitsVds3_1->At(i);
    if (hit->GetClusterSize() > sizeTh) {
      fRealClusters[3]++;
      if (fhClusterVerticalDist_Vds3) fhClusterVerticalDist_Vds3->Fill(hit->GetY());
    }
  }

  UDataTable* hitsVds4_0 = outNode->GetDataTable("Hits Vds4_0");
  fRealClusters[4] = 0;
  for (int i = 0; i < hitsVds4_0->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hitsVds4_0->At(i);
    if (hit->GetClusterSize() > sizeTh) {
      fRealClusters[4]++;
      if (fhClusterVerticalDist_Vds4a) fhClusterVerticalDist_Vds4a->Fill(hit->GetY());
    }
  }
  UDataTable* hitsVds4_1 = outNode->GetDataTable("Hits Vds4_1");
  fRealClusters[5] = 0;
  for (int i = 0; i < hitsVds4_1->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hitsVds4_1->At(i);
    if (hit->GetClusterSize() > sizeTh) {
      fRealClusters[5]++;
      if (fhClusterVerticalDist_Vds4a) fhClusterVerticalDist_Vds4a->Fill(hit->GetY());
    }
  }

  UDataTable* hitsVds4_2 = outNode->GetDataTable("Hits Vds4_2");
  fRealClusters[6] = 0;
  for (int i = 0; i < hitsVds4_2->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hitsVds4_2->At(i);
    if (hit->GetClusterSize() > sizeTh) {
      fRealClusters[6]++;
      if (fhClusterVerticalDist_Vds4b) fhClusterVerticalDist_Vds4b->Fill(hit->GetY());
    }
  }

  UDataTable* hitsVds4_3 = outNode->GetDataTable("Hits Vds4_3");
  fRealClusters[7] = 0;
  for (int i = 0; i < hitsVds4_3->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hitsVds4_3->At(i);
    if (hit->GetClusterSize() > sizeTh) {
      fRealClusters[7]++;
      if (fhClusterVerticalDist_Vds4b) fhClusterVerticalDist_Vds4b->Fill(hit->GetY());
    }
  }
}

//_____________________________________________________________________
void Na61HitProducerModule::LocalToGlobal(int ii, UDataTable* hits) {
  double alpha = fParams->GetRotZ(ii);
  double beta = fParams->GetRotY(ii);
  double gamma = fParams->GetRotX(ii);
  double VolumeX = fParams->GetVolumeX(ii);
  double VolumeY = fParams->GetVolumeY(ii);
  double VolumeZ = fParams->GetVolumeZ(ii);

//cout<<"sensor="<<ii<<" "<<VolumeX<<" "<<VolumeY<<" "<<VolumeZ<<" "<<alpha<< " "<<beta<<" "<<gamma<<endl;

  double sa = TMath::Sin(alpha);
  double ca = TMath::Cos(alpha);
  double sb = TMath::Sin(beta);
  double cb = TMath::Cos(beta);
  double sg = TMath::Sin(gamma);
  double cg = TMath::Cos(gamma);

  for (int i = 0; i < hits->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hits->At(i);
    double x1 = hit->GetX();
    double y1 = hit->GetY();
    double z1 = 0;

    double x2;
    double y2;
    double z2;

    // x2 =  ca*x1  +  sa*y1;
    // y2 = -sa*x1  +  ca*y1;
    // z2 =     z1;

    // x2 =  TMath::Cos(alpha)*x1  +  TMath::Sin(alpha)*z1;
    // y2 =  y1;
    // z2 = -TMath::Sin(alpha)*x1  +  TMath::Cos(alpha)*z1;
    // ZY:
    // x2 =  TMath::Cos(alpha)*TMath::Cos(beta)*x1  +  TMath::Sin(alpha)*y1 +  TMath::Cos(alpha)*TMath::Sin(beta)*z1;
    // y2 = -TMath::Sin(alpha)*TMath::Cos(beta)*x1  +  TMath::Cos(alpha)*y1 -  TMath::Sin(alpha)*TMath::Sin(beta)*z1;
    // z2 = -TMath::Sin(beta)*x1                                            +  TMath::Cos(beta)*z1;
    // ZX:

    // x2 =  ca*x1  +  sa*cg*y1 +  sa*sg*z1;
    // y2 = -sa*x1  +  ca*cg*y1 +  ca*sg*z1;
    // z2 =   0*x1  -     sg*y1 +     cg*z1;

    // ZXY
    x2 =  (ca * cb - sa * sg * sb) * x1 + sa * cg * y1 + ( ca * sb + sa * sg * cb) * z1;
    y2 = -(sa * cb + ca * sg * sb) * x1 + ca * cg * y1 + (-sa * sb + ca * sg * cb) * z1;
    z2 =                  -cg * sb * x1      - sg * y1 +                   cg * cb * z1;

    hit->SetX(x2 + VolumeX);
    hit->SetY(y2 + VolumeY);
    hit->SetZ(z2 + VolumeZ);

    if (TMath::Abs(z1-z2) >0.7) {
		cout<<"Strange data: initial:"<<x1<<" "<<y1<<" "<<z1<<" final:"<<x2<<" "<<y2<<" "<<z2<<endl;
		cout<<"sensor="<<ii<<" "<<VolumeX<<" "<<VolumeY<<" "<<VolumeZ<<" "<<alpha<< " "<<beta<<" "<<gamma<<endl;
    }
    // if(ii==0)cout<<"VolumeX="<<VolumeX<<" VolumeY="<<VolumeY<<" VolumeZ="<<VolumeZ<<"  x1="<<x1<<"  y1="<<y1<<"  z2="<<z2<<" z2+V ="<<hit->GetZ()<<endl;

    if (fhHitsXVsZ) fhHitsXVsZ->Fill(hit->GetZ(), hit->GetX());
    if (fhHitsYVsZ && ii < 6) fhHitsYVsZ->Fill(hit->GetZ(), hit->GetY());
  }
}

//_____________________________________________________________________
void Na61HitProducerModule::MakeClusters(TString sensorname, UDataTable* pixels, UDataTable* hits) {
  SetaxayToIndArray(pixels);

  // improved algorithm. It runs independently of first seed value.
  // it make cluster from all that are touching. Size of final cluster is not
  // restricted

  // fClusterArray->Clear();
  TString name = sensorname;

  // cout<<"     -------------------------------> Started  MakeClusters-------------------------<<<< "<<endl;
  for (int i = 0; i < pixels->GetEntries(); i++) {
    USensorPixel* pixel = (USensorPixel*)pixels->At(i);

    // cout<<"i="<<i<<" line = "<<pixel->GetLine()<<"  column="<<pixel->GetColumn()<<endl;

    if (!pixel) Error("MakeClusters", "something is wrong, check it out (1)");

    if (pixel->IsUsed()) continue;  // means that pixel was already used

    fNRes = 0;

    TObjArray* cluster = new TObjArray();
    cluster->Clear();
    cluster->Add(pixel);

    // cout<<"------------> start check ring: i="<<pixel->GetLine()<<" j="<<pixel->GetColumn()<<endl;
    CheckRing2(pixel->GetLine(), pixel->GetColumn(), pixels, cluster);

    // we are here means that cluster is generated

    // loop over cluster
    // to make sensor hit.
    // Hit position from center of gravity

    int M = 0;

    double x = 0;
    double y = 0;
    double z = 0;
    for (int ic = 0; ic < cluster->GetEntries(); ic++) {
      USensorPixel* pix = (USensorPixel*)cluster->At(ic);

      if (!pix) Error("MakeClusters", "something is wrong with your cluster, check it out");

      if (fJura) {
        x = x + pix->GetPaX();
        y = y + pix->GetPaY();
      } else {
        x = x + pix->GetNaX();
        y = y + pix->GetNaY();
      }

      M++;
      // cout<<" loop over cluster: line="<<pix->GetLine()<<" col="<<pix->GetColumn()<<" M="<<M<<endl;
    }
    // cout<<" M = "<<M<<endl;

    x = x / (double)M;
    y = y / (double)M;

    if (name.Contains("Vds1")) z = 0;
    if (name.Contains("Vds2")) z = 50;
    if (name.Contains("Vds3")) z = 100;
    if (name.Contains("Vds4")) z = 150;

    USensorHit* hit = new USensorHit(x, y, z, M);
    hit->SetSensorName(sensorname);
    hit->SetClusterLine(pixel->GetLine());
    hit->SetClusterColumn(pixel->GetColumn());
    // add hit to data array

    ////////////// Cuts used to extract beam ions at the spot center //////////////////////////
    // double xoff   = -2.657;
    // if(fJura)xoff  =  2.2;
    // if((TMath::Abs(x+xoff)<0.2) || (z>5)) hits->Add(hit); // Saleve run 633 and Jura run 635
    //////////////////////////////////////////////////////////////////////////////////////////
    // else delete hit;

    hits->Add(hit);

    // store cluster is clusterArray
    // fClusterArray->Add(cluster);

    // if(M>2) cout<<" M="<<M<<"   i="<<i<<" j="<<j<<"  fCluster.Entries: "<<fCluster.GetEntries()<<endl;
    // if(M<3) continue;

    delete cluster;
  }

  /*
  if(name.Contains("Vds3")){
    for(int i=0;i<hits->GetEntries();i++){
      USensorHit* hit1 =  (USensorHit*)hits->At(i);
      for(int j=i+1;j<hits->GetEntries();j++){
        USensorHit* hit2 =  (USensorHit*)hits->At(j);

        if(hit1->GetX()==hit2->GetX() && hit1->GetY()==hit2->GetY()){
          //cout<<"hit1: i="<<i<<"  "<<hit1->GetX()<<" "<<hit1->GetY()<<" "<<hit1->GetZ()<<" line="<<hit1->GetClusterLine()<<" coll="<<hit1->GetClusterColumn()<<endl;
          //cout<<"hit2: j="<<j<<"  "<<hit2->GetX()<<" "<<hit2->GetY()<<" "<<hit2->GetZ()<<" line="<<hit2->GetClusterLine()<<" coll="<<hit2->GetClusterColumn()<<endl;
          //int ii;
          //cin>>ii;
        }
      }
    }
    } */
}

//_____________________________________________________________________
void Na61HitProducerModule::SetaxayToIndArray(UDataTable* pixels) {
  int ind = 0;
  // resetung the faxayToIn array
  for (int i = 0; i < findx + 1; i++)
    for (int j = 0; j < findy + 1; j++) faxayToInd[i][j] = -1;

  for (int i = 0; i < pixels->GetEntries(); i++) {
    USensorPixel* pix = (USensorPixel*)pixels->At(i);

    if (pix->GetLine() > 576) {
      cout << "Wrong data: line = " << pix->GetLine() << endl;
      continue;
    }

    if (pix->GetColumn() > 1152) {
      cout << "Wrong data: column = " << pix->GetColumn() << endl;
      continue;
    }

    faxayToInd[pix->GetLine()][pix->GetColumn()] = ind;

    ind++;
  }
}

//_____________________________________________________________________
void Na61HitProducerModule::CheckRing2(int i, int j, UDataTable* pixels, TObjArray* cluster)

{
  for (int ik = i - 1; ik < i + 2; ik++) {
    for (int jk = j - 1; jk < j + 2; jk++) {
      if (ik == i && jk == j) continue;         // local restriction
      if ((ik < 1) || (jk < 1)) continue;       // do not go beyong lower limits
      if ((ik > 576) || (jk > 1152)) continue;  // do not go beyong upper limits

      // cout<<"ik="<<ik<<" jk="<<jk<<endl;
      // int ii;
      // cin>>ii;

      // check inherited restriction
      bool restricted = false;

      for (int ir = 0; ir < fNRes; ir++) {
        // if(TMath::Abs(fResIndx[ir]-ik)<2 && TMath::Abs(fResIndy[ir]-jk)<2)restricted=true;
        if (TMath::Abs(fResIndx[ir] - ik) < 1 && TMath::Abs(fResIndy[ir] - jk) < 1) restricted = true;
      }

      if (restricted) continue;

      int ind = faxayToInd[ik][jk];
      if (ind == -1) continue;

      // cout<<"ind="<<ind<<endl;

      USensorPixel* pixel = (USensorPixel*)pixels->At(ind);
      if (pixel->IsUsed()) continue;

      cluster->Add(pixel);  // add pixel to cluster
      pixel->SetUsed();

      // cout<<"ind="<<ind<<"  pixel line: "<<pixel->GetLine()<<"  pixel column: "<<pixel->GetColumn()<<endl;

      // make restrictions
      fResIndx[fNRes] = i;
      fResIndy[fNRes] = j;
      fNRes++;
      if (fNRes > 500) {
        Error("CheckRing2", "fNRes is not supposed to be as big >500, check it out");
        if (fNRes > 1499) cout << "fNRes=" << fNRes << endl;
      }
      CheckRing2(ik, jk, pixels, cluster);
    }
  }
}

//_____________________________________________________________________________
double Na61HitProducerModule::FillClusterDistr(UDataTable* hits, TH1F* h_multi, TH1F* h_multiLt1, TH1F* h_size) {
  int entries = 0;
  if (h_multi) {
    entries = hits->GetEntries();
    h_multi->Fill(entries);
  }

  int ClustersLt1 = 0;

  for (int i = 0; i < hits->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hits->At(i);

    if (!hit) Error(" FillClusterDistr", "something is wrong, check it out");

    // if(hit->GetClusterSize()>150)cout<<"Encoutered big cluster size="<<hit->GetClusterSize()<<endl;
    if (h_size) h_size->Fill(hit->GetClusterSize());

    if (hit->GetClusterSize() > 1) ClustersLt1++;
  }

  if (h_multiLt1) h_multiLt1->Fill(ClustersLt1);

  if (entries)
    return (double)ClustersLt1 / (double)entries;
  else
    return 0;
}

//_____________________________________________________________________________
void Na61HitProducerModule::FillHSaxayWithClusters(UDataTable* hits, TH2F* hh_clust, TH2F* hh_clustpos, TH2F* hh_largeclustpos, int ii) {
  // hh_clust->Reset("ICES");

  int sizeArr[8] = {120, 120, 120, 120, 120, 120, 120, 120};

  for (int i = 0; i < hits->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hits->At(i);

    if (!hit) Error(" FillHSaxayWithClusters", "something is wrong, check it out");
    // cout<<hit->GetClusterLine()<<" "<<hit->GetClusterColumn()<<" "<<hit->GetClusterSize()<<endl;

    int size = hit->GetClusterSize();

    if (hh_clust) hh_clust->SetBinContent(hit->GetClusterLine(), hit->GetClusterColumn(), hit->GetClusterSize());
    if (hh_clustpos) hh_clustpos->Fill(hit->GetX(), hit->GetY());
    if (hh_largeclustpos && size > sizeArr[ii]) hh_largeclustpos->Fill(hit->GetX(), hit->GetY());
    // cout<<hit->GetX()<<"  "<<hit->GetY()<<endl;
  }
}
//_____________________________________________________________________________
void Na61HitProducerModule::FillHSaxayWithPixels(UDataTable* pixs, TH2F* hh_pix) {
  // hh_pix->Reset("ICES");

  for (int i = 0; i < pixs->GetEntries(); i++) {
    USensorPixel* pix = (USensorPixel*)pixs->At(i);

    if (!pix) Error(" FillHSaxayWithPixels", "something is wrong, check it out");
    // cout<<hit->GetClusterLine()<<" "<<hit->GetClusterColumn()<<" "<<hit->GetClusterSize()<<endl;

    // if(pix->GetLine()<2)cout<<"line = "<<pix->GetLine()<<endl;
    // if(pix->GetColumn()<2)cout<<"column = "<<pix->GetColumn()<<endl;

    // if(hh_pix)hh_pix->SetBinContent(pix->GetLine(),pix->GetColumn());
    hh_pix->Fill(pix->GetLine(), pix->GetColumn());
  }
}

//_____________________________________________________________
void Na61HitProducerModule::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61HitProducerModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);

  fhSignalToNoise->Sumw2();
  fhSignalToNoise->Scale(1. / fNorm);
}

//____________________________________________________________________
void Na61HitProducerModule::Print(Option_t* option) const {
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>

  TString opt(option);
  opt.ToLower();

  Na61Module::Print(option);
  if (opt.Contains("d")) cout << endl << "  Original author: Paweł Staszel" << endl << "  Last Modifications: " << endl << "    $Author: Staszel $" << endl << "    $Date: 2016/10/25$" << endl << "    $Revision: 1.0 $ " << endl << endl << "-------------------------------------------------" << endl;
}
