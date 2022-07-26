#ifndef _VdTestSG_VdTest_h_
#define _VdTestSG_VdTest_h_
#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include <evt/Event.h>

class TH2D;
class TH1D;

namespace VdTestSG{
  class VdTest:public boost::noncopyable, public fwk::VModule{
  public:
    VdTest (){ }
    ~VdTest (){ }
    fwk::VModule::EResultFlag Init ();
    fwk::VModule::EResultFlag Process (evt::Event & event,const utl::AttributeMap & attr);
    fwk::VModule::EResultFlag Finish ();
  private:
    //static bool AbsComp (double i, double j);
    


    
	TH1D* hVertex_X_dist_TPC_fine;
	TH1D* hVertex_Y_dist_TPC_fine;
	TH1D* hVertex_Z_dist_TPC_fine;	
	
	TH1D* hist_VDTPCDiff_x;
	TH1D* hist_VDTPCDiff_y;
	TH1D* hist_VDTPCDiff_z;
	
	TH1D* hist_vertexVD_z;
	TH1D* hist_vertexVD_x;
	TH1D* hist_vertexVD_y;

	TH1D* hist_vertexTPCNew_x;
	TH1D* hist_vertexTPCNew_y;
	TH1D* hist_vertexTPCNew_z;
	
	TH2D* hist_vertexVdTpc_z;
	TH2D* hist_vertexVdTpc_x;
	TH2D* hist_vertexVdTpc_y;
	
	TH1D* hist_NtracksTPC;
   	TH1D* hist_NtracksVDTPC;
	TH1D* hist_NtracksVD;
	TH2D* hist_NTracksCorrMatchedVD;
    TH2D* hist_NTracksCorrMatchedTPC;
 
	TH1D* hist_NClustersTPC;
    TH1D* hist_NClustersVTPC;
    
    int prevRun;

    REGISTER_MODULE ("VdTestSG", VdTest,"$Id: UserSkeleton.h 7924 2012-10-05 21:50:50Z dveberic $");
  };
}
#endif // _VdTestSG_VdTest_h_
