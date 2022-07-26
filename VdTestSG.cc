#include "VdTestSG.h"
#include "Na61VdParametersManager.h"
#include "Na61VdParameters.h"

#include <fwk/CentralConfig.h>
#include <utl/Branch.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineUnits.h>
#include <utl/PhysicalConst.h>
#include <utl/WithUnit.h>
#include <det/Silicon.h>
#include <det/Detector.h>
#include <det/PSD.h>
#include <evt/proc/Silicon.h>
#include <evt/raw/Silicon.h>

#include <evt/RecEvent.h>

#include <evt/rec/Vertex.h>
#include <evt/rec/VertexTrack.h>
#include <evt/rec/Track.h>
#include <evt/Event.h>
#include <det/Target.h>
#include <sstream>
#include <string>
#include <set>
#include "TMath.h"
#include <det/Detector.h>
#include <evt/SimEvent.h>
#include <evt/Event.h>
#include <det/Beam.h>
#include <utl/Point.h>
#include <det/MagneticField.h>
#include <det/MagneticFieldTracker.h>
#include <evt/rec/Trigger.h>

#include <TFile.h>
#include <boost/lexical_cast.hpp>
#include "TH2I.h"
#include "TH3I.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TLegend.h"
#include "TPolyLine.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TBox.h"
#include <evt/IndexedObjectLinker.h>

#include <sstream>
using namespace fwk;
using namespace utl;
using namespace std;


//namespace tmpCmp{
//   bool absComp(double i, double j) {return (abs(i) < abs(j));}
//}

namespace VdTestSG{
  fwk::VModule::EResultFlag VdTest::Init (){


    hist_vertexVD_z = new TH1D("hist_vertexVD_z","hist_vertexVD_z", 300, -605, -600);
    hist_vertexVD_x = new TH1D("hist_vertexVD_x","hist_vertexVD_x", 600, -1.5, 1.5);
    hist_vertexVD_y = new TH1D("hist_vertexVD_y","hist_vertexVD_y", 600, -1.5, 1.5);

    hist_vertexTPCNew_x = new TH1D("hist_vertexTPCNew_x","", 3000, -3, 3);
    hist_vertexTPCNew_y = new TH1D("hist_vertexTPCNew_y","", 3000, -3, 3);
    hist_vertexTPCNew_z = new TH1D("hist_vertexTPCNew_z","", 3000, -605, -600);

    hVertex_X_dist_TPC_fine = new TH1D("hVertex_X_dist_TPC_fine", "", 500, -6, 6);
    hVertex_Y_dist_TPC_fine = new TH1D("hVertex_Y_dist_TPC_fine", "", 500, -6, 6);
    hVertex_Z_dist_TPC_fine = new TH1D("hVertex_Z_dist_TPC_fine", "", 500, -6100, -5900);

    hist_VDTPCDiff_x = new TH1D("hist_VDTPCDiff_x", "", 500, -1, 1);
    hist_VDTPCDiff_y = new TH1D("hist_VDTPCDiff_y", "", 500, -1, 1);
    hist_VDTPCDiff_z = new TH1D("hist_VDTPCDiff_z", "", 500, -3, 3);

    hist_vertexVdTpc_x = new TH2D("hist_vertexVdTpc_x", "", 500, -1, 1, 500, -1, 1);
    hist_vertexVdTpc_y = new TH2D("hist_vertexVdTpc_y", "", 500, -1, 1, 500, -1, 1);
    hist_vertexVdTpc_z = new TH2D("hist_vertexVdTpc_z", "", 500, -604, -601, 500, -604, -601);


    hist_NtracksTPC = new TH1D("hist_NtracksTPC","", 1501, -0.5, 1500.5);
    hist_NtracksVDTPC = new TH1D("hist_NtracksVDTPC","", 501, -0.5, 500.5);
    hist_NtracksVD = new TH1D("hist_NtracksVD","", 501, -0.5, 500.5);
    hist_NTracksCorrMatchedVD= new TH2D("hist_NTracksCorrMatchedVD"," ", 501, -0.5, 500.5, 501, -0.5, 500.5);
    hist_NTracksCorrMatchedTPC= new TH2D("hist_NTracksCorrMatchedTPC"," ", 1501, -0.5, 1500.5, 501, -0.5, 500.5);

    hist_NClustersTPC = new TH1D("hist_NClustersTPC"," ", 250, 0.5, 250.5);
    hist_NClustersVTPC = new TH1D("hist_NClustersVTPC"," ", 150, 0.5, 150.5);


    prevRun = 0;

    return eSuccess;
  }

  fwk::VModule::EResultFlag VdTest::Process (evt::Event & e,const utl::AttributeMap &){
    //Check if event is not empty, otherwise take next event.
    const evt::EventHeader& eventHeader = e.GetEventHeader();
    if(!eventHeader.IsInitialized()){
      WARNING("No EventHeader is present in Event. This is a serious problem!!! PLEASE DO CHECK WHAT IS GOING ON IF THIS EFFECT IS NOT EXPECTED!!!");
      return eFailure;
    }
    //Get event header info and report progress.
    const int runNumber = eventHeader.GetRunNumber();
    const utl::TimeStamp& timeStamp = eventHeader.GetTime();
    det::Detector& detector = det::Detector::GetInstance();
    detector.Update(timeStamp, runNumber);
    evt::RecEvent& recEvent=e.GetRecEvent();

    if (det::Detector::GetInstance().GetTarget().GetStatus()==det::TargetConst::eIn || det::Detector::GetInstance().GetTarget().GetStatus()==det::TargetConst::eFull) ;
    else {cout<<"Something is wrong! Missing target"<<endl; return eSuccess;}
    //const utl::Point& mainVertexPos = det::Detector::GetInstance().GetTarget().GetCenterPosition();


    if (int(eventHeader.GetRunNumber()) != prevRun) {
      const unsigned int calibrun = Na61VdParameters::FindCalibRuns(runNumber);
      //TFile& matchParamsFile = *TFile::Open(Form("matchParams-%06d.root", calibrun), "READ");
      TFile & matchParamsFile = *TFile::Open(Form("/afs/cern.ch/user/a/amerzlay/public/XeLa150/VDmatchingfiles/matchParams-%0d.root",calibrun),"READ");
      Na61VdParameters& na61VdParameters = *Na61VdParametersManager::Instance()->GetVdParams();
      na61VdParameters.SetRunId(runNumber);
      na61VdParameters.SetMatchParamsFile(&matchParamsFile);
      na61VdParameters.SetdOffZ(0);
      na61VdParameters.SetN(0);
      na61VdParameters.Init();
      matchParamsFile.Close();

      prevRun = eventHeader.GetRunNumber();
    }


    bool VD_flag=recEvent.HasPrimaryVertex(evt::rec::VertexConst::ePrimaryVD);
    const utl::Point& VdVertex=VD_flag?recEvent.GetPrimaryVertex(evt::rec::VertexConst::ePrimaryVD).GetPosition() : utl::Point(0,0,0);


    bool TPCvertex_flag=recEvent.HasPrimaryVertex(evt::rec::VertexConst::eUnknown);
    const utl::Point& TPCVertexNew=TPCvertex_flag?recEvent.GetPrimaryVertex(evt::rec::VertexConst::eUnknown).GetPosition() : utl::Point(0,0,0);


    if ((TPCvertex_flag) && (VD_flag)){
      // Add shift relative to TPC
      const double fOffsetToTpcX = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcX();
      const double fOffsetToTpcY = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcY();
      const double fOffsetToTpcZ = Na61VdParametersManager::Instance()->GetVdParams()->GetOffsetToTpcZ();

      cout<<"fOffsetToTpcX="<<fOffsetToTpcX<<endl;

      hVertex_X_dist_TPC_fine->Fill(TPCVertexNew.GetX()*10 - (VdVertex.GetX()*10 - fOffsetToTpcX)); //in mm
      hVertex_Y_dist_TPC_fine->Fill(TPCVertexNew.GetY()*10 - (VdVertex.GetY()*10 - fOffsetToTpcY));
      hVertex_Z_dist_TPC_fine->Fill(TPCVertexNew.GetZ()*10 - (VdVertex.GetZ()*10 - fOffsetToTpcZ));

      hist_VDTPCDiff_x->Fill(TPCVertexNew.GetX() - VdVertex.GetX());
      hist_VDTPCDiff_y->Fill(TPCVertexNew.GetY() - VdVertex.GetY());
      hist_VDTPCDiff_z->Fill(TPCVertexNew.GetZ() - VdVertex.GetZ());

      hist_vertexVD_x->Fill(VdVertex.GetX());
      hist_vertexVD_y->Fill(VdVertex.GetY());
      hist_vertexVD_z->Fill(VdVertex.GetZ());


      hist_vertexTPCNew_x->Fill(TPCVertexNew.GetX());
      hist_vertexTPCNew_y->Fill(TPCVertexNew.GetY());
      hist_vertexTPCNew_z->Fill(TPCVertexNew.GetZ());	

      hist_vertexVdTpc_x->Fill(VdVertex.GetX(), TPCVertexNew.GetX());	
      hist_vertexVdTpc_y->Fill(VdVertex.GetY(), TPCVertexNew.GetY());	
      hist_vertexVdTpc_z->Fill(VdVertex.GetZ(), TPCVertexNew.GetZ());	

      double NtracksVD = e.GetProcEvent().GetSilicon().GetTracks().size();
      cout<<"N_v-tracks="<<recEvent.GetNumberOf<evt::rec::VertexTrack>()<<endl;
      cout<<"N_tracks="<<recEvent.GetNumberOf<evt::rec::Track>()<<endl;
      cout<<"NtracksVD="<<NtracksVD<<endl;
      hist_NtracksTPC->Fill(recEvent.GetNumberOf<evt::rec::VertexTrack>());
      hist_NtracksVD->Fill(NtracksVD);

      int VDmatchedtrackcounter=0;
      for(std::set<evt::Index<evt::rec::Track> >::const_iterator t=e.GetProcEvent().GetSilicon().GetTracks().begin(); t!=e.GetProcEvent().GetSilicon().GetTracks().end(); t++){
        if (!(e.GetRecEvent().Has(*t))) continue;
        const evt::rec::Track& Vdtrack=e.GetRecEvent().Get(*t);
        if ((Vdtrack.GetNumberOfClusters(evt::rec::TrackConst::eVD)>2) && (Vdtrack.GetNumberOfClusters(evt::rec::TrackConst::eAll)-Vdtrack.GetNumberOfClusters(evt::rec::TrackConst::eVD)>0))
          VDmatchedtrackcounter++;
      }
      cout<<"VD matched track counter="<<VDmatchedtrackcounter<<endl;
      hist_NtracksVDTPC->Fill(VDmatchedtrackcounter);
      hist_NTracksCorrMatchedVD->Fill(e.GetProcEvent().GetSilicon().GetTracks().size(),VDmatchedtrackcounter);
      hist_NTracksCorrMatchedTPC->Fill(recEvent.GetNumberOf<evt::rec::VertexTrack>(),VDmatchedtrackcounter);

      //	int vdtrackscounter=0;
      const evt::rec::Vertex& mainVertex = recEvent.GetMainVertex(); //select main vertex
      for(evt::rec::VertexTrackIndexIterator iVtxTrackId = mainVertex.DaughterTracksBegin(); iVtxTrackId != mainVertex.DaughterTracksEnd(); ++iVtxTrackId){   
        const evt::rec::VertexTrack& vtxTrack = recEvent.Get(*iVtxTrackId);

        if(!vtxTrack.HasTrack()) continue; //remove vtxtracks without tracks
        if(vtxTrack.GetStatus()!=0) continue; // bad track status
        if(vtxTrack.GetCharge()==0) continue; // remove tracks with q=0
        const evt::rec::Track& track = recEvent.Get(vtxTrack.GetTrackIndex()); //get tracks

        // if(track.GetNumberOfClusters(evt::rec::TrackConst::eVD)>2) vdtrackscounter++;

        hist_NClustersTPC->Fill(track.GetNumberOfClusters(evt::rec::TrackConst::eAll)-track.GetNumberOfClusters(evt::rec::TrackConst::eVD));
        hist_NClustersVTPC->Fill(track.GetNumberOfClusters( evt::rec::TrackConst::eVTPC1)+track.GetNumberOfClusters( evt::rec::TrackConst::eVTPC2));

      }
    }


    return eSuccess;
  }

  fwk::VModule::EResultFlag VdTest::Finish (){
    cout<<"saving root file"<<endl;
    TFile *file = new TFile ("fileVertex.root", "RECREATE");

    hVertex_X_dist_TPC_fine->Write();
    hVertex_Y_dist_TPC_fine->Write();
    hVertex_Z_dist_TPC_fine->Write();	
    hist_VDTPCDiff_x->Write();
    hist_VDTPCDiff_y->Write();
    hist_VDTPCDiff_z->Write();

    hist_vertexVD_z->Write();
    hist_vertexVD_x->Write();
    hist_vertexVD_y->Write();  

    hist_vertexTPCNew_x->Write();
    hist_vertexTPCNew_y->Write();
    hist_vertexTPCNew_z->Write();

    hist_NtracksTPC->Write();
    hist_NtracksVDTPC->Write();
    hist_NtracksVD->Write();
    hist_NTracksCorrMatchedVD->Write();
    hist_NTracksCorrMatchedTPC->Write();

    hist_vertexVdTpc_x->Write();
    hist_vertexVdTpc_y->Write();
    hist_vertexVdTpc_z->Write();

    hist_NClustersTPC->Write();
    hist_NClustersVTPC->Write();

    file->Close();

    return eSuccess;
  }

}
