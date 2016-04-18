#ifndef Physics_Analysis_XisP3Pi_H
#define Physics_Analysis_XisP3Pi_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"
#include "PartPropSvc/PartPropSvc.h"  


#include <string>
#include <TTree.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TClonesArray.h>

#include "TROOT.h"
#include "TBenchmark.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "VertexFit/WTrackParameter.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

class XisP3Pi:public Algorithm{

	public:
		XisP3Pi(const std::string& name, ISvcLocator* pSvcLocator);
		StatusCode initialize();
		StatusCode execute();
		StatusCode finalize();  

	private:
		IPartPropSvc *p_PartPropSvc;
		HepPDT::ParticleDataTable* m_particleTable;

        void corgen(HepMatrix &, HepVector &, int );
        void corset(HepSymMatrix &, HepMatrix &, int );
        void calibration(RecMdcKalTrack * , HepVector &, int );
	
		//global parameter
		double m_Ecms;
		double m_EnergySpread;
		int m_psiType;

		std::string m_OutputFileName;
		TFile *saveFile;
		TTree *TreeAna;
		TTree *NbInfo;

		int m_saveTopo;
		int m_saveMCTruth;
		int m_saveCutFlow;


		// Declare cut parameters initinal value
		double m_vr0cut;
		double m_vz0cut;
		double m_trkAngCut;
		double TDC_min;
		double TDC_max;
		double m_gammaAngleCut;
		double m_barrelEnergyThreshold;
		double m_endEnergyThreshold;
		double barcut;
		double endmin;
		double endmax;
		double gamiso;
		double mc_cor;
		int m_NGamma;

		// MC topology 
		int m_idxmc;
		int m_drank[500];
		int m_pdgid[500];
		int m_motheridx[500];
		int m_motherpid[500];
		//MC Truth
		TLorentzVector *xep;
		TLorentzVector *xem;
		TLorentzVector *xgamma1;
		TLorentzVector *xgamma2;
		TLorentzVector *xgg;
		TLorentzVector *xee;
		TLorentzVector *xeta;
		double xmee;
		double xmgg;
		int flagP;
		int issig;

		int runid;
		int evtid;
		int nevt;
		int nNEUTRAL;
		int nCHARGED;
		int nTRACKS;
		
		//Good Gamma
		TClonesArray *GammaAll;
		double costheta_gamma[500];
		double energy_gamma[500];
		double TDCtime[500];
		double isoAngle[500];
		double showerde;

		//track  info		
		double vx[3],Evx[3];
		double costheta_chrgd[2];
		double Rxy[2];
		double Rz[2];
		double Rvxy[2];
		double Rvz[2];
		double pKal[2];
		double pxKal[2];
		double pyKal[2];
		double pzKal[2];
		double deemc[2];
		double eop[2];

		//pid info.
		int flag_ep;
		int flag_em;
		double prob[2][3];

		//for vertex fit
		double vtxchisq;	

		//for 4C-fit
		double chi4C;
		TLorentzVector *p_ep;			//after fit
		TLorentzVector *p_em; 
		TLorentzVector *p_gamma1;		
		TLorentzVector *p_gamma2;		
		TLorentzVector *p_ee; 
		TLorentzVector *p_gg; 
	
		TLorentzVector *p_uep;			//no fit
		TLorentzVector *p_uem; 
		TLorentzVector *p_ugamma1;		
		TLorentzVector *p_ugamma2;		
		TLorentzVector *p_uee; 
		TLorentzVector *p_ugg; 
		
		TLorentzVector *p_gamma1b;		
		TLorentzVector *p_gamma2b; 
		double angee;
		double mee;
		double mgg;
		double cosgam1b;
		double cosgam2b;
		double cosgamb;

		Double_t em_pull_0;
		Double_t em_pull_1;
		Double_t em_pull_2;
		Double_t em_pull_3;
		Double_t em_pull_4;
		Double_t ep_pull_0;
		Double_t ep_pull_1;
		Double_t ep_pull_2;
		Double_t ep_pull_3;
		Double_t ep_pull_4;

		Double_t m_rconv;
		Double_t m_xconv1;
		Double_t m_yconv1;
		Double_t m_zconv1;
		Double_t m_rconv1;
		Double_t m_xconv2;
		Double_t m_yconv2;
		Double_t m_zconv2;
		Double_t m_rconv2;
		Double_t m_xiep;
		Double_t m_deltaxy;

		Double_t m_deltaz1;
		Double_t m_deltaz2;

		Double_t m_lep;
		Double_t m_psipair;
		//                Double_t m_dgamma;
		Double_t MEE;
		Double_t m_vx_x;
		Double_t m_vx_y;
		Double_t m_vx_r;
		Double_t m_thetaeg1;
		Double_t m_thetaeg2;
		Double_t m_cthep;
		Double_t m_ptrkp;
		Double_t m_ptrkm;
		Double_t m_mgamma;
		Double_t m_egamma;
		Double_t m_theta;
		Double_t m_cosTheta;
		Double_t m_phi;
		Double_t m_rp;
		Double_t m_re;
		Double_t m_deltaeq;
		Double_t m_case;

};

#endif 
