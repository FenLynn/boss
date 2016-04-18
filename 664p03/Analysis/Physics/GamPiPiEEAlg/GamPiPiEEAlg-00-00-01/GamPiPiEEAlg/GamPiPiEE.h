#ifndef Physics_Analysis_GamPiPiEE_H
#define Physics_Analysis_GamPiPiEE_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"


#include <string>
#include <TTree.h>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TClonesArray.h>

class GamPiPiEE:public Algorithm{

	public:
		GamPiPiEE(const std::string& name, ISvcLocator* pSvcLocator);
		StatusCode initialize();
		StatusCode execute();
		StatusCode finalize();  

	private:
		//global parameter
		double m_Ecms;
		double m_EnergySpread;
		int m_psiType;

		std::string m_OutputFileName;
		TFile *saveFile;
		TTree *TreeAna;
		TTree *TopoTree;
		TTree *NbInfo;


		int m_saveTopo;
		int m_saveMCTruth;
		int m_saveTopoTree;
		int m_saveCutFlow;

		int runid;
		int evtid;
		int nevt;
		int nNEUTRAL;
		int nCHARGED;
		int nTRACKS;
		int nGamma;
	
		int useVxfitCut;
		int useKmfitCut;

		// Declare cut parameters initinal value
		int m_NGamma;
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
		
		//for MC topology and MCtruth
		unsigned int conexc_flag1;
		unsigned int conexc_flag2;
		int m_idxmc;
		int m_drank[500];
		int m_pdgid[500];
		int m_motheridx[500];
		int m_motherpid[500];
		int isSig;
		TLorentzVector *xep;
		TLorentzVector *xem;
		TLorentzVector *xpip;
		TLorentzVector *xpim;
		TLorentzVector *xgamma;
		TLorentzVector *xetap;
		TLorentzVector *xee;

		//Good Gamma
		TClonesArray *GammaAll;
		double costheta_gamma[500];
		double energy_gamma[500];
		double TDCtime[500];
		double isoAngle[500];
		double showerde;

		//track  info		
		double vx[3],Evx[3];
		double costheta_chrgd[4];
		double Rxy[4];
		double Rz[4];
		double Rvxy[4];
		double Rvz[4];
		double m_deTrk[4];
		double m_pTrk[4];
		double m_epTrk[4];
		double m_PtTrk[4];
		int m_hitsTrk[4];
		
		//pid info.
		int flag_pip;
		int flag_pim;
		int flag_ep;
		int flag_em;
		double chiDeDx[4][3];
		double chiTof[4][3];
		double chiTof1[4][3];
		double chiTof2[4][3];
		double chisq_pid[4][3];
		double prob[4][3];
		
		//for gamma_conv
		double m_rconv;
		double m_xconv1;
		double m_yconv1;
		double m_zconv1;
		double m_rconv1;
		double m_xconv2;
		double m_yconv2;
		double m_zconv2;
		double m_rconv2;
		double m_xiep;
		double m_deltaxy;
		double m_cthep;
		double m_ptrkp;
		double m_ptrkm;
		double m_mgamma;
		double m_egamma;
		double m_theta;
		double m_cosTheta;
		double m_phi;
		
		//for vertex fit
		double vtxchisq;	
		double vtxchisqee;	
		double vtxchisqpipi;	

		//for 4C-fit
		int min4C_index;
		double chi4C[500];
		TClonesArray *pep;			//for results after fit
		TClonesArray *pem; 	
		TClonesArray *ppip; 	
		TClonesArray *ppim;
		TClonesArray *pgamma;
		TClonesArray *petap;
		TClonesArray *pee;
		
		double min_chi4C;
		double min_chi4C_unfit;
		TLorentzVector *min_ep;			//for mini-chi4C 4-momentum after fit
		TLorentzVector *min_em; 
		TLorentzVector *min_pip; 
		TLorentzVector *min_pim; 
		TLorentzVector *min_gamma; 
		TLorentzVector *min_etap; 
		TLorentzVector *min_ee; 
		TLorentzVector *rec_ee; 

		TLorentzVector *min_ep_unfit;			//no 4C
		TLorentzVector *min_em_unfit; 
		TLorentzVector *min_pip_unfit; 
		TLorentzVector *min_pim_unfit; 
		TLorentzVector *min_gamma_unfit; 
		TLorentzVector *min_etap_unfit; 
		TLorentzVector *min_ee_unfit; 
		TLorentzVector *rec_ee_unfit; 
	

		TLorentzVector *hel_ep;
		TLorentzVector *hel_pip;
		TLorentzVector *hel_em;
		TLorentzVector *hel_pim;
		TLorentzVector *hel_gam;

		double hel_cos[5]; 
		double mpipi;
		double angee;
		double angee_unfit;
		double cosee;
		double mee;
		double mee_unfit;
		double metap;
		double mgee;
		double m_rec_ee;
		double m_rec_ee_unfit;
};

#endif 
