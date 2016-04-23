#ifndef Physics_Analysis_EtaPiPiEE_H
#define Physics_Analysis_EtaPiPiEE_H 

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

class EtaPiPiEE:public Algorithm{

	public:
		EtaPiPiEE(const std::string& name, ISvcLocator* pSvcLocator);
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
		int gamiso;
		int m_NGamma;
		int mc_cor;
		
		//for MC topology and MCtruth
		int m_idxmc;
		int m_drank[500];
		int m_pdgid[500];
		int m_motheridx[500];
		int m_motherpid[500];

		int issig;
		double xmee;
		TLorentzVector *xep;
		TLorentzVector *xem;
		TLorentzVector *xpip;
		TLorentzVector *xpim;
		TLorentzVector *xgamma1;
		TLorentzVector *xgamma2;
		TLorentzVector *xeta;
		TLorentzVector *xee;
		TLorentzVector *xetap;

		//event info
		int runid;
		int evtid;
		int nNEUTRAL;
		int nCHARGED;
		int nTRACKS;
		int nGamma;
		
		//Good Gamma
		TClonesArray *GammaAll;
		double costheta_gamma[500];
		double energy_gamma[500];
		double TDCtime[500];
		double isoAngle[500];
		double showerde;

		//vertex  info		
		double vx[3],Evx[3];
	
		//track
		double costheta_chrgd[4];
		double Rxy[4];
		double Rz[4];
		double Rvxy[4];
		double Rvz[4];
		double pKal[4];
		double pxKal[4];
		double pyKal[4];
		double pzKal[4];
		double deemc[4];
		double eop[4];

		//pid info.
		int flag_pip;
		int flag_pim;
		int flag_ep;
		int flag_em;

		double prob[4][3];

		//for vertex fit
		double vtxchisq;	

		//for 4C-fit
		double chi4C;
		TLorentzVector *p_pip;			//after fit
		TLorentzVector *p_pim; 
		TLorentzVector *p_ep; 
		TLorentzVector *p_em; 
		TLorentzVector *p_gamma1;		
		TLorentzVector *p_gamma2;		
		TLorentzVector *p_pipi; 
		TLorentzVector *p_ee; 
		TLorentzVector *p_gg; 
		TLorentzVector *p_ggpipi; 
		TLorentzVector *p_recpipi;

		TLorentzVector *p_upip;			//no fit
		TLorentzVector *p_upim; 
		TLorentzVector *p_uep; 
		TLorentzVector *p_uem; 
		TLorentzVector *p_ugamma1;		
		TLorentzVector *p_ugamma2;		
		TLorentzVector *p_upipi; 
		TLorentzVector *p_uee; 
		TLorentzVector *p_ugg; 
		TLorentzVector *p_uggpipi; 

		double m_gg;
		double m_ggpipi;
		double m_uee;
		double m_ee;
		double m_recpipi;

		double angee;
		double highe_eop;
		double lowe_eop;

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
		Double_t pip_pull_0;
		Double_t pip_pull_1;
		Double_t pip_pull_2;
		Double_t pip_pull_3;
		Double_t pip_pull_4;
		Double_t pim_pull_0;
		Double_t pim_pull_1;
		Double_t pim_pull_2;
		Double_t pim_pull_3;
		Double_t pim_pull_4;

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

		//gamma conversion (zhangyt)
		double vtxchie; 
		double m_epemx;
		double m_epemy;
		double m_epemz;
		double m_epemxy;
		double m_mepr;
		double m_mepcenterx;
		double m_mepcentery;
		double m_mepcenterz;
		double m_memr;
		double m_memcenterx;
		double m_memcentery;
		double m_memcenterz;
		double m_epemxxorigin1e;
		double m_epemyxorigin1e;
		double m_epemzxorigin1e;
		double m_epemxyxorigin1e;
		double m_rconvz;
};

#endif 
