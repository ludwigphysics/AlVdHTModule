//--------------------------------------------
// Input/output module for G4Na61 analysis
// Author: Paweł Staszel
//--------------------------------------------
#if !defined NA61_Na61AlHitProducerModule
#include "Na61AlHitProducerModule.h"
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
#ifndef UTIL_UChipPixel
#include "UChipPixel.h"
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
// ClassImp(Na61AlHitProducerModule);

//____________________________________________________________________
Na61AlHitProducerModule::Na61AlHitProducerModule() {
  // Default constructor. DO NOT USE
  SetState(kSetup);
  fProductionMode = 1;
  fNorm = 0;
}
//____________________________________________________________________
Na61AlHitProducerModule::Na61AlHitProducerModule(const char* name, const char* title) : Na61Module(name, title) {
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
void Na61AlHitProducerModule::DefineHistograms() {
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

  //for (int i = 0; i < 34; i++) {
  //fhLargeClusterPosMap[i] = new TH2F(Form("hLargeClusterPosMap_%s", fSensorNames[i].Data()), "", 106, -5.3, 5.3, 212, -10.6, 10.6);
  //}

  for (int i = 0; i < 34; i++) {
    fhPixelMap[i] = new TH2F(Form("hPixelMap_%s", fAlSensorNames[i].Data()), "", binsx, 0, xmax, binsy, 0, ymax);
    fhPixelPosMap[i] = new TH2F(Form("hPixelPosMap_%s", fAlSensorNames[i].Data()), "", binsx, -7.8, 7.8, binsy, -15.3, 15.3);
    fhClusterMap[i] = new TH2F(Form("hClusterMap_%s", fAlSensorNames[i].Data()), "", binsx, 0, xmax, binsy, 0, ymax);
    fhClusterPosMap[i] = new TH2F(Form("hClusterPosMap_%s", fAlSensorNames[i].Data()), "", binsx, -7.8, 7.8, binsy, -15.3, 15.3);
    fhClusters[i] = new TH1F(Form("hClusters_%s", fAlSensorNames[i].Data()), "", 300, 0, 300);
    //fhClustersLt1[i] = new TH1F(Form("hClustersLt1_%s", fSensorNames[i].Data()), "", 300, 0, 300);
    fhClusterSize[i] = new TH1F(Form("hClusterSize_%s", fAlSensorNames[i].Data()), "", 350, 0, 350);
  }

  fhHitsXVsZ = new TH2F("HitsXVsZ", "", 5000, -1, 151, 1000, -50, 50.);
  fhHitsYVsZ = new TH2F("HitsYVsZ", "", 5000, -1, 151, 1000, -82, 82.);

  fhClusterVerticalDist_Vds1 = new TH1F("hClusterVerticalDist_Vds1", "", 220, -82, 82);
  fhClusterVerticalDist_Vds2 = new TH1F("hClusterVerticalDist_Vds2", "", 220, -82, 82);
  fhClusterVerticalDist_Vds3 = new TH1F("hClusterVerticalDist_Vds3", "", 220, -82, 82);
  fhClusterVerticalDist_Vds4 = new TH1F("hClusterVerticalDist_Vds4", "", 220, -82, 82);
  fhClusterHorizontalDist_Vds1 = new TH1F("hClusterHorizontalDist_Vds1", "", 220, -50, 50);
  fhClusterHorizontalDist_Vds2 = new TH1F("hClusterHorizontalDist_Vds2", "", 220, -50, 50);
  fhClusterHorizontalDist_Vds3 = new TH1F("hClusterHorizontalDist_Vds3", "", 220, -50, 50);
  fhClusterHorizontalDist_Vds4 = new TH1F("hClusterHorizontalDist_Vds4", "", 220, -50, 50);


  fhSignalToNoise = new TH1F("hSignalToNoise", "", 8, 0, 8);
  // Define histograms. They are:
  // <fill in here>

  gDirectory->cd("..");
}

//____________________________________________________________________
void Na61AlHitProducerModule::SetRunInfo(TFile* file) {
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
void Na61AlHitProducerModule::Init() {
  // Job-level initialisation
  SetState(kInit);

  if (fJura)
    fParams = (Na61VdParametersManager::Instance())->GetJuraAlVdArmParams();
  else
    fParams = (Na61VdParametersManager::Instance())->GetSaleveAlVdArmParams();

  cout << "Setting controll parameters in ArmParameters: fSensorId=" << fSensorId << " fRunId=" << fRunId << endl;
  fParams->SetDrotX(fSensorId, fDrotX);
  fParams->SetDrotY(fSensorId, fDrotY);
  fParams->SetDrotZ(fSensorId, fDrotZ);
  fParams->SetVolumeDz(fSensorId, fDz);

  fParams->SetRunId(fRunId);
  fParams->Init();
}

//____________________________________________________________________
void Na61AlHitProducerModule::Begin() {
  // Run-level initialisation
  SetState(kBegin);
}

//____________________________________________________________________
void Na61AlHitProducerModule::Event(UEventNode* inNode, UEventNode* outNode)

{
  // Per event method
  SetState(kEvent);

  for (int i = 0; i < 34; i++) {
    UDataTable* pixtab;
    if (fJura)
      pixtab = outNode->GetDataTable(Form("Jura Pixels %s cleaned", fAlSensorNames[i].Data()));
    else
      pixtab = outNode->GetDataTable(Form("Saleve Pixels %s cleaned", fAlSensorNames[i].Data()));

    if (!pixtab) {
      if (fJura)
        pixtab = inNode->GetDataTable(Form("Jura Pixels %s cleaned", fAlSensorNames[i].Data()));
      else
        pixtab = inNode->GetDataTable(Form("Saleve Pixels %s cleaned", fAlSensorNames[i].Data()));
    }

    if (pixtab) {
      UDataTable* hitstab = new UDataTable(Form("Hits %s", fAlSensorNames[i].Data()));
      MakeClusters(fAlSensorNames[i], pixtab, hitstab);
      outNode->AddDataTable(hitstab);
      //cout<<"sensor(i)="<<i<<"  fJura = "<<fJura<<"  hits: "<<hitstab->GetEntries()<<endl;

      //if (!fProductionMode) { // for debuging
      FillHSaxayWithPixels(pixtab, fhPixelMap[i],fhPixelPosMap[i]);
      //}

      FillHSaxayWithClusters(hitstab, fhClusterMap[i], fhClusterPosMap[i], i);
      //double sigToNoise = FillClusterDistr(hitstab, fhClusters[i], fhClustersLt1[i], fhClusterSize[i]);
      //fhSignalToNoise->Fill(i, sigToNoise);
      fNorm++;
      LocalToGlobal(i, hitstab);  // move to global coordinate system
    }
  }

  FillHorizontalDistributions(outNode);

  for (int i = 0; i < 34; i++) {
    UDataTable* hits = outNode->GetDataTable(Form("Hits %s", fAlSensorNames[i].Data()));
    fClusters[i] = hits->GetEntries();
  }
}

//_____________________________________________________________________
void Na61AlHitProducerModule::FillHorizontalDistributions(UEventNode* outNode) {
  //
  int sizeTh = 0;

  for(int it=0; it<34; it++){
    UDataTable* hits = outNode->GetDataTable(Form("Hits %s",fAlSensorNames[it].Data()));
    fRealClusters[0] = 0;
    for (int i = 0; i < hits->GetEntries(); i++) {
      USensorHit* hit = (USensorHit*)hits->At(i);
      if (hit->GetClusterSize() > sizeTh) {
	if (fhHitsXVsZ) fhHitsXVsZ->Fill(hit->GetZ(), hit->GetX());
	if (fhHitsYVsZ) fhHitsYVsZ->Fill(hit->GetZ(), hit->GetY());
	fRealClusters[0]++;
	if (fhClusterVerticalDist_Vds1 && fAlSensorNames[it].Contains("Al1_"))   fhClusterVerticalDist_Vds1->Fill(hit->GetY());
	if (fhClusterVerticalDist_Vds2 && fAlSensorNames[it].Contains("Al2_"))   fhClusterVerticalDist_Vds2->Fill(hit->GetY());
	if (fhClusterVerticalDist_Vds3 && fAlSensorNames[it].Contains("Al3_"))   fhClusterVerticalDist_Vds3->Fill(hit->GetY());
	if (fhClusterVerticalDist_Vds4 && fAlSensorNames[it].Contains("Al4_"))   fhClusterVerticalDist_Vds4->Fill(hit->GetY());
	if (fhClusterHorizontalDist_Vds1 && fAlSensorNames[it].Contains("Al1_")) fhClusterHorizontalDist_Vds1->Fill(hit->GetX());
	if (fhClusterHorizontalDist_Vds2 && fAlSensorNames[it].Contains("Al2_")) fhClusterHorizontalDist_Vds2->Fill(hit->GetX());
	if (fhClusterHorizontalDist_Vds3 && fAlSensorNames[it].Contains("Al3_")) fhClusterHorizontalDist_Vds3->Fill(hit->GetX());
	if (fhClusterHorizontalDist_Vds4 && fAlSensorNames[it].Contains("Al4_")) fhClusterHorizontalDist_Vds4->Fill(hit->GetX());
      }
    }
  }
  
  return; 
}

//_____________________________________________________________________
void Na61AlHitProducerModule::LocalToGlobal(int ii, UDataTable* hits) {

  //
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
    
    if (TMath::Abs(z1-z2) >2.7) {
      cout<<"Strange data: initial:"<<x1<<" "<<y1<<" "<<z1<<" final:"<<x2<<" "<<y2<<" "<<z2<<endl;
      cout<<"sensor="<<ii<<" "<<VolumeX<<" "<<VolumeY<<" "<<VolumeZ<<" "<<alpha<< " "<<beta<<" "<<gamma<<endl;
    }
    // if(ii==0)cout<<"VolumeX="<<VolumeX<<" VolumeY="<<VolumeY<<" VolumeZ="<<VolumeZ<<"  x1="<<x1<<"  y1="<<y1<<"  z2="<<z2<<" z2+V ="<<hit->GetZ()<<endl;
    
    
  }
}

//_____________________________________________________________________
void Na61AlHitProducerModule::MakeClusters(TString sensorname, UDataTable* pixels, UDataTable* hits) 
{

  SetaxayToIndArray(pixels);
  
  // improved algorithm. It runs independently of first seed value.
  // it make cluster from all that are touching. Size of final cluster is not
  // restricted

  // fClusterArray->Clear();
  TString name = sensorname;

  // cout<<"     -------------------------------> Started  MakeClusters-------------------------<<<< "<<endl;
  for (int i = 0; i < pixels->GetEntries(); i++) {
    UChipPixel* pixel = (UChipPixel*)pixels->At(i);

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
      UChipPixel* pix = (UChipPixel*)cluster->At(ic);

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
void Na61AlHitProducerModule::SetaxayToIndArray(UDataTable* pixels) {
  int ind = 0;
  // resetung the faxayToIn array
  for (int i = 0; i < findx + 1; i++)
    for (int j = 0; j < findy + 1; j++) faxayToInd[i][j] = -1;

  for (int i = 0; i < pixels->GetEntries(); i++) {
    UChipPixel* pix = (UChipPixel*)pixels->At(i);

    if (pix->GetLine() > 512) {
      cout << "Wrong data: line = " << pix->GetLine() << endl;
      continue;
    }

    if (pix->GetColumn() > 1024) {
      cout << "Wrong data: column = " << pix->GetColumn() << endl;
      continue;
    }

    faxayToInd[pix->GetLine()][pix->GetColumn()] = ind;

    ind++;
  }
}

//_____________________________________________________________________
void Na61AlHitProducerModule::CheckRing2(int i, int j, UDataTable* pixels, TObjArray* cluster)

{
  for (int ik = i - 1; ik < i + 2; ik++) {
    for (int jk = j - 1; jk < j + 2; jk++) {
      if (ik == i && jk == j) continue;         // local restriction
      if ((ik < 1) || (jk < 1)) continue;       // do not go beyong lower limits
      if ((ik > 512) || (jk > 1024)) continue;  // do not go beyong upper limits

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

      UChipPixel* pixel = (UChipPixel*)pixels->At(ind);
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
double Na61AlHitProducerModule::FillClusterDistr(UDataTable* hits, TH1F* h_multi, TH1F* h_multiLt1, TH1F* h_size) {
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
void Na61AlHitProducerModule::FillHSaxayWithClusters(UDataTable* hits, TH2F* hh_clust, TH2F* hh_clustpos, int /* ii */) {
  // hh_clust->Reset("ICES");

  //int sizeArr[8] = {120, 120, 120, 120, 120, 120, 120, 120};

  for (int i = 0; i < hits->GetEntries(); i++) {
    USensorHit* hit = (USensorHit*)hits->At(i);

    if (!hit) Error(" FillHSaxayWithClusters", "something is wrong, check it out");
    // cout<<hit->GetClusterLine()<<" "<<hit->GetClusterColumn()<<" "<<hit->GetClusterSize()<<endl;

    //int size = hit->GetClusterSize();

    if (hh_clust) hh_clust->SetBinContent(hit->GetClusterLine(), hit->GetClusterColumn(), hit->GetClusterSize());
    if (hh_clustpos) hh_clustpos->Fill(hit->GetX(), hit->GetY());
    //if (hh_largeclustpos && size > sizeArr[ii]) hh_largeclustpos->Fill(hit->GetX(), hit->GetY());
    // cout<<hit->GetX()<<"  "<<hit->GetY()<<endl;
  }
}
//_____________________________________________________________________________
void Na61AlHitProducerModule::FillHSaxayWithPixels(UDataTable* pixs, TH2F* hh_pix,  TH2F* hh_pixpos) {
  // hh_pix->Reset("ICES");

  for (int i = 0; i < pixs->GetEntries(); i++) {
    UChipPixel* pix = (UChipPixel*)pixs->At(i);

    if (!pix) Error(" FillHSaxayWithPixels", "something is wrong, check it out");
    // cout<<hit->GetClusterLine()<<" "<<hit->GetClusterColumn()<<" "<<hit->GetClusterSize()<<endl;

    // if(pix->GetLine()<2)cout<<"line = "<<pix->GetLine()<<endl;
    // if(pix->GetColumn()<2)cout<<"column = "<<pix->GetColumn()<<endl;

    // if(hh_pix)hh_pix->SetBinContent(pix->GetLine(),pix->GetColumn());
    hh_pix->Fill(pix->GetLine(), pix->GetColumn());
    Float_t x=0;
    Float_t y=0;
    if(fJura){
      x = pix->GetPaX();
      y = pix->GetPaY();
    }else{
      x = pix->GetNaX();
      y = pix->GetNaY();
    }
    hh_pixpos->Fill(x,y);

  }
}

//_____________________________________________________________
void Na61AlHitProducerModule::End() {
  // Run-level finalisation
  SetState(kEnd);
}

//____________________________________________________________________
void Na61AlHitProducerModule::Finish() {
  // Job-level finalisation
  SetState(kFinish);

  fhSignalToNoise->Sumw2();
  fhSignalToNoise->Scale(1. / fNorm);
}

//____________________________________________________________________
void Na61AlHitProducerModule::Print(Option_t* option) const {
  // Print module information
  // In addition this module defines the Option:
  // <fill in here>

  TString opt(option);
  opt.ToLower();

  Na61Module::Print(option);
  if (opt.Contains("d")) cout << endl << "  Original author: Paweł Staszel" << endl << "  Last Modifications: " << endl << "    $Author: Staszel $" << endl << "    $Date: 2016/10/25$" << endl << "    $Revision: 1.0 $ " << endl << endl << "-------------------------------------------------" << endl;
}
