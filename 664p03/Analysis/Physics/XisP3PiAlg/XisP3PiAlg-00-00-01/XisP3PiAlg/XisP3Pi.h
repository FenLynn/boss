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
		int issig;
		TLorentzVector *xp_xis;
		TLorentzVector *xp_xi;
		TLorentzVector *xp_xipim;
		TLorentzVector *xp_lam;
		TLorentzVector *xp_lampim;
		TLorentzVector *xp_prop;
		TLorentzVector *xp_pi0;
		TLorentzVector *xp_gamma1;
		TLorentzVector *xp_gamma2;
		TLorentzVector *xp_xibar;
		TLorentzVector *xp_lambar;
		TLorentzVector *xp_xibarpip;
		TLorentzVector *xp_prom;
		TLorentzVector *xp_lambarpip;

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
		int nGamma;

		//track  info		
		double vx[3],Evx[3];
		int nGood;
		double costheta_chrgd[100];
		double Rxy[100];
		double Rz[100];
		double Rvxy[100];
		double Rvz[100];
		double pKal[100];
		double pxKal[100];
		double pyKal[100];
		double pzKal[100];
		double deemc[100];
		double eop[100];

		//pid info.
		int flag_prop;
		int flag_pim[2];


		//for 1C-fit
		double chi1C;
		int idx_gamma1;
		int idx_gamma2;
		TLorentzVector *p_gamma1;		
		TLorentzVector *p_gamma2;		
		TLorentzVector *p_pi0; 
	
		//for vertex fit
		int nVtxok;
		int flag_vtx;
		double chi_vtxLam[100];
		double chi_vtxXi[100];
		double chi_svtxXi[100];
		double chi_vtx[100];
		double m_lambda[100];
		double m_xi[100];
		double m_xis[100];
		double m_xibar[100];
		double m_ctau[100];
		double m_len[100];
		double m_lenerr[100];
		double flag_lampim[100];
		double flag_xipim[100];
		
		TClonesArray *p_prop;
		TClonesArray *p_lampim;
		TClonesArray *p_xipim;
		TClonesArray *p_lambda;
		TClonesArray *p_xi;
		TClonesArray *p_xis;
		TClonesArray *p_xibar;





};

#endif 
