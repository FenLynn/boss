#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/INTupleSvc.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "McTruth/McParticle.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"
#include "EmcRecEventModel/RecEmcHit.h"
#include "EvTimeEvent/RecEsTime.h"

#include "Identifier/Identifier.h"
#include "McTruth/McParticle.h"
#include "McTruth/DecayMode.h"
#include "McTruth/MdcMcHit.h"
#include "McTruth/TofMcHit.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/MucMcHit.h"  
#include "McTruth/McEvent.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include "GamPiPiEEAlg/GamPiPiEE.h"

#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"

#include "EventNavigator/EventNavigator.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcHit.h"
#include "EmcRecEventModel/RecEmcShower.h"
#include "PartPropSvc/PartPropSvc.h"

#include "TMath.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TRandom.h"

#include "GammaConv/GammaConv.h"

//const double twopi = 6.2831853;
//const double pi = 3.1415927;
const double mpi = 0.13957;
const double mmu = 0.105658;
const double me = 0.000511;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double mel  = 0.000510998910;
const double Ejpsi=3.097;
const double Epsip=3.686;
const double Econt=3.650;
//const double Econt=3.773;

//const double velc = 29.9792458;  tof_path unit in cm.
const double velc = 299.792458;   // tof path unit in mm

typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

//int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6,Ncut7,Ncut8,Ncut9,Ng1,Ng2,Ng3;
int Ncut0 = 0;		//total event
int Ncut1 = 0;		//Pass good charged tracks cut
int Ncut2 = 0;		//Pass good photon
int Ncut3 = 0;		//Pass PID
int Ncut4 = 0;		//Pass Vertex
int Ncut5 = 0;		//Pass 4C
int Ncut6 = 0;		//final
int Ncut7 = 0;
int Ncut8 = 0;
int Ncut9 = 0;
/////////////////////////////////////////////////////////////////////////////

GamPiPiEE::GamPiPiEE(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) {
		//Declare the properties
		declareProperty("EnergySpread",  m_EnergySpread = 0.0013);
		declareProperty("OutputFileName",  m_OutputFileName = "GamPiPiEE_test.root");
		declareProperty("Ecms",  m_Ecms = 3.686);
		declareProperty("psiType", m_psiType = 0);
		declareProperty("saveTopo", m_saveTopo = 1);//need to be re-evaluated when running different samples(1 for MC)
		declareProperty("saveTopoTree", m_saveTopoTree = 0);//need to be re-evaluated when running different samples(1 for MC)
		declareProperty("saveMCTruth", m_saveMCTruth = 1);//need to be re-evaluated when running different samples(only 1 for exclusiveMC)
		declareProperty("saveCutFlow", m_saveCutFlow = 1);//need to be re-evaluated when running different samples(1 for MC)

		declareProperty("Vr0cut", m_vr0cut=1.0);
		declareProperty("Vz0cut", m_vz0cut=10.0);
		declareProperty("TrkAngCut", m_trkAngCut=0.93);
		declareProperty("TDC_min", TDC_min=0.);
		declareProperty("TDC_max", TDC_max=14.);
		declareProperty("GammaAngleCut", m_gammaAngleCut=10.0);
		declareProperty("BarrelEnergyThreshold", m_barrelEnergyThreshold=0.025);
		declareProperty("EndEnergyThreshold",m_endEnergyThreshold=0.050);
		declareProperty("barcut",barcut=0.8);
		declareProperty("endmin",endmin=0.86);
		declareProperty("endmax",endmax=0.92);
		declareProperty("gamiso",gamiso=1);
		declareProperty("m_NGamma",m_NGamma=1);

		declareProperty("UseVxfitCut", useVxfitCut = 1);
		declareProperty("UseKmfitCut", useKmfitCut = 1);     
	}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode GamPiPiEE::initialize(){
	MsgStream log(msgSvc(), name());
	std::cout<<"initialize has started!"<<std::endl;
	log << MSG::INFO << "in initialize()" << endmsg;
	StatusCode status;

	//Output name
	TString s_OutputFileName(m_OutputFileName);
	s_OutputFileName.ReplaceAll("[\"","");
	s_OutputFileName.ReplaceAll("\"]","");

	//Tree and File Modified
	saveFile = new TFile(s_OutputFileName, "recreate");
	TreeAna = new TTree("save", "save");
	NbInfo = new TTree("nbinfo","nbinfo");
	TopoTree = new TTree("topoall","topoall");

	//TClonesArray initialize 
	GammaAll=new TClonesArray("TLorentzVector");
	pep=new TClonesArray("TLorentzVector");
	pem=new TClonesArray("TLorentzVector");
	ppip=new TClonesArray("TLorentzVector");
	ppim=new TClonesArray("TLorentzVector");
	pgamma=new TClonesArray("TLorentzVector");
	petap=new TClonesArray("TLorentzVector");
	pee=new TClonesArray("TLorentzVector");

	//Lorentz Vector 
	min_ep=new TLorentzVector();
	min_em=new TLorentzVector();
	min_pip=new TLorentzVector();
	min_pim=new TLorentzVector();
	min_gamma=new TLorentzVector();
	min_etap=new TLorentzVector();
	min_ee=new TLorentzVector();
	rec_ee=new TLorentzVector();
	
	min_ep_unfit=new TLorentzVector();
	min_em_unfit=new TLorentzVector();
	min_pip_unfit=new TLorentzVector();
	min_pim_unfit=new TLorentzVector();
	min_gamma_unfit=new TLorentzVector();
	min_etap_unfit=new TLorentzVector();
	min_ee_unfit=new TLorentzVector();
	rec_ee_unfit=new TLorentzVector();
	
	hel_ep=new TLorentzVector();
	hel_pip=new TLorentzVector();
	hel_em=new TLorentzVector();
	hel_pim=new TLorentzVector();
	hel_gam=new TLorentzVector();

	xep=new TLorentzVector();
	xem=new TLorentzVector();
	xpip=new TLorentzVector();
	xpim=new TLorentzVector();
	xgamma=new TLorentzVector();
	xetap=new TLorentzVector();
	xee=new TLorentzVector();

	//mc info
	if(m_saveTopo == 1)
	{
		TreeAna->Branch("indexmc", &m_idxmc, "indexmc/I");
		TreeAna->Branch("pdgid", m_pdgid, "pdgid[indexmc]/I");
		TreeAna->Branch("drank", m_drank, "drank[indexmc]/I");
		TreeAna->Branch("motheridx", m_motheridx, "motheridx[indexmc]/I");
		TreeAna->Branch("motherpid", m_motherpid, "motherpid[indexmc]/I");
		TreeAna->Branch("flag1", &conexc_flag1, "flag1/i");
		TreeAna->Branch("flag2", &conexc_flag2, "flag2/i");
	}
	if(m_saveMCTruth == 1)
	{
		TreeAna->Branch("xpip",&xpip,32000,0);
		TreeAna->Branch("xpim",&xpim,32000,0);
		TreeAna->Branch("xep",&xep,32000,0);
		TreeAna->Branch("xem",&xem,32000,0);
		TreeAna->Branch("xgamma",&xgamma,32000,0);
		TreeAna->Branch("xetap",&xetap,32000,0);
		TreeAna->Branch("xee",&xee,32000,0);
		TreeAna->Branch("isSig",&isSig,"isSig/I");
	}
	if(m_saveTopoTree == 1)
	{
		TopoTree->Branch("indexmc", &m_idxmc, "indexmc/I");
		TopoTree->Branch("pdgid", m_pdgid, "pdgid[indexmc]/I");
		TopoTree->Branch("drank", m_drank, "drank[indexmc]/I");
		TopoTree->Branch("motheridx", m_motheridx, "motheridx[indexmc]/I");
		TopoTree->Branch("motherpid", m_motherpid, "motherpid[indexmc]/I");
		TopoTree->Branch("flag1", &conexc_flag1, "flag1/i");
		TopoTree->Branch("flag2", &conexc_flag2, "flag2/i");
	}
	TreeAna->Branch("runid", &runid, "runid/I");
	TreeAna->Branch("evtid", &evtid, "evtid/I");
	TreeAna->Branch("nevt", &nevt, "nevt/I");
	TreeAna->Branch("nNEUTRAL", &nNEUTRAL, "nNEUTRAL/I");
	TreeAna->Branch("nCHARGED",&nCHARGED,"nCHARGED/I");
	TreeAna->Branch("nTRACKS",&nTRACKS,"nTRACKS/I");
	TreeAna->Branch("nGamma", &nGamma, "nGamma/I");

	TreeAna->Branch("GammaAll","TClonesArray",&GammaAll,256000,0);
	TreeAna->Branch("costheta_gamma", costheta_gamma, "costheta_gamma[nGamma]/D");
	TreeAna->Branch("energy_gamma", energy_gamma, "energy_gamma[nGamma]/D");
	TreeAna->Branch("TDCtime", TDCtime , "TDCtime[nGamma]/D");
	TreeAna->Branch("isoAngle", isoAngle, "isoAngle[nGamma]/D");
	TreeAna->Branch("showerde", &showerde, "showerde/D");

	TreeAna->Branch("vx", vx, "vx[3]/D");
	TreeAna->Branch("Evx", Evx, "Evx[3]/D");
	TreeAna->Branch("costheta_chrgd", costheta_chrgd, "costheta_chrgd[4]/D");
	TreeAna->Branch("Rxy", Rxy, "Rxy[4]/D");
	TreeAna->Branch("Rz", Rz, "Rz[4]/D");
	TreeAna->Branch("Rvxy", Rvxy, "Rvxy[4]/D");
	TreeAna->Branch("Rvz", Rvz, "Rvz[4]/D");
	TreeAna->Branch("pTrk", m_pTrk, "pTrk[4]/D");
	TreeAna->Branch("deTrk", m_deTrk, "deTrk[4]/D");
	TreeAna->Branch("epTrk", m_epTrk, "epTrk[4]/D");
	TreeAna->Branch("PtTrk", m_PtTrk, "PtTrk[4]/D");
	TreeAna->Branch("hitsTrk", m_hitsTrk, "hitsTrk[4]/I");
	TreeAna->Branch("hel_cos", hel_cos, "hel_cos[5]/D");

	TreeAna->Branch("flag_pip", &flag_pip, "flag_pip/I");
	TreeAna->Branch("flag_pim", &flag_pim, "flag_pim/I");
	TreeAna->Branch("flag_ep", &flag_ep, "flag_ep/I");
	TreeAna->Branch("flag_em", &flag_em, "flag_em/I");
	TreeAna->Branch("chiDeDx", chiDeDx, "chiDeDx[4][3]/D");
	TreeAna->Branch("chiTof",chiTof,"chiTof[4][3]/D");
	TreeAna->Branch("chiTof1", chiTof1, "chiTof1[4][3]/D");
	TreeAna->Branch("chiTof2", chiTof2, "chiTof2[4][3]/D");		
	TreeAna->Branch("chisq_pid", chisq_pid, "chisq_pid[4][3]/D");
	TreeAna->Branch("prob", prob, "prob[4][3]/D");

	//for gamma_conv
	TreeAna->Branch("m_rconv",&m_rconv,"m_rconv/D");
	TreeAna->Branch("m_xconv1",&m_xconv1,"m_xconv1/D");
	TreeAna->Branch("m_yconv1",&m_yconv1,"m_yconv1/D");
	TreeAna->Branch("m_zconv1",&m_zconv1,"m_zconv1/D");
	TreeAna->Branch("m_rconv1",&m_rconv1,"m_rconv1/D");
	TreeAna->Branch("m_xconv2",&m_xconv2,"m_xconv2/D");
	TreeAna->Branch("m_yconv2",&m_yconv2,"m_yconv2/D");
	TreeAna->Branch("m_zconv2",&m_zconv2,"m_zconv2/D");
	TreeAna->Branch("m_rconv2",&m_rconv2,"m_rconv2/D");
	TreeAna->Branch("m_xiep",&m_xiep,"m_xiep/D");
	TreeAna->Branch("m_deltaxy",&m_deltaxy,"m_deltaxy/D");
	TreeAna->Branch("m_cthep",&m_cthep,"m_cthep/D");
	TreeAna->Branch("m_ptrkp",&m_ptrkp,"m_ptrkp/D");
	TreeAna->Branch("m_ptrkm",&m_ptrkm,"m_ptrkm/D");
	TreeAna->Branch("m_mgamma",&m_mgamma,"m_mgamma/D");
	TreeAna->Branch("m_egamma",&m_egamma,"m_egamma/D");
	TreeAna->Branch("m_theta",&m_theta,"m_theta/D");
	TreeAna->Branch("m_cosTheta",&m_cosTheta,"m_cosTheta/D");
	TreeAna->Branch("m_phi",&m_phi,"m_phi/D");

	TreeAna->Branch("vtxchisq", &vtxchisq, "vtxchisq/D");
	TreeAna->Branch("vtxchisqee", &vtxchisqee, "vtxchisqee/D");
	TreeAna->Branch("vtxchisqpipi", &vtxchisqpipi, "vtxchisqpipi/D");

	TreeAna->Branch("chi4C", chi4C, "chi4C[nevt]/D");
	TreeAna->Branch("min4C_index", &min4C_index, "min4C_index/D");
	TreeAna->Branch("pep","TClonesArray",&pep,256000,0);
	TreeAna->Branch("pem","TClonesArray",&pem,256000,0);
	TreeAna->Branch("ppip","TClonesArray",&ppip,256000,0);
	TreeAna->Branch("ppim","TClonesArray",&ppim,256000,0);
	TreeAna->Branch("pgamma","TClonesArray",&pgamma,256000,0);
	TreeAna->Branch("petap","TClonesArray",&petap,256000,0);
	TreeAna->Branch("pee","TClonesArray",&pee,256000,0);
	
	TreeAna->Branch("min_chi4C", &min_chi4C, "min_chi4C/D");
	TreeAna->Branch("min_chi4C_unfit", &min_chi4C_unfit, "min_chi4C_unfit/D");
	TreeAna->Branch("min_ep",&min_ep,32000,0);
	TreeAna->Branch("min_em",&min_em,32000,0);
	TreeAna->Branch("min_pip",&min_pip,32000,0);
	TreeAna->Branch("min_pim",&min_pim,32000,0);
	TreeAna->Branch("min_gamma",&min_gamma,32000,0);
	TreeAna->Branch("min_etap",&min_etap,32000,0);
	TreeAna->Branch("min_ee",&min_ee,32000,0);
	TreeAna->Branch("rec_ee",&rec_ee,32000,0);
	TreeAna->Branch("min_ep_unfit",&min_ep_unfit,32000,0);
	TreeAna->Branch("min_em_unfit",&min_em_unfit,32000,0);
	TreeAna->Branch("min_pip_unfit",&min_pip_unfit,32000,0);
	TreeAna->Branch("min_pim_unfit",&min_pim_unfit,32000,0);
	TreeAna->Branch("min_gamma_unfit",&min_gamma_unfit,32000,0);
	TreeAna->Branch("min_etap_unfit",&min_etap_unfit,32000,0);
	TreeAna->Branch("min_ee_unfit",&min_ee_unfit,32000,0);
	TreeAna->Branch("rec_ee_unfit",&rec_ee_unfit,32000,0);
	
	TreeAna->Branch("hel_ep",&hel_ep,32000,0);
	TreeAna->Branch("hel_pip",&hel_pip,32000,0);
	TreeAna->Branch("hel_em",&hel_em,32000,0);
	TreeAna->Branch("hel_pim",&hel_pim,32000,0);
	TreeAna->Branch("hel_gam",&hel_gam,32000,0);

	TreeAna->Branch("angee", &angee, "angee/D");
	TreeAna->Branch("angee_unfit", &angee_unfit, "angee_unfit/D");
	TreeAna->Branch("cosee", &cosee, "cosee/D");
	TreeAna->Branch("mee", &mee, "mee/D");
	TreeAna->Branch("mee_unfit", &mee_unfit, "mee_unfit/D");
	TreeAna->Branch("metap", &metap, "metap/D");
	TreeAna->Branch("mgee", &mgee, "mgee/D");
	TreeAna->Branch("m_rec_ee", &m_rec_ee, "m_rec_ee/D");
	TreeAna->Branch("m_rec_ee_unfit", &m_rec_ee_unfit, "m_rec_ee_unfit/D");
	TreeAna->Branch("mpipi", &mpipi, "mpipi/D");

	//for number of events passing every cut
	if(m_saveCutFlow == 1)
	{
		NbInfo->Branch("Ncut0", &Ncut0, "Ncut0/I");
		NbInfo->Branch("Ncut1", &Ncut1, "Ncut1/I");
		NbInfo->Branch("Ncut2", &Ncut2, "Ncut2/I");
		NbInfo->Branch("Ncut3", &Ncut3, "Ncut3/I");
		NbInfo->Branch("Ncut4", &Ncut4, "Ncut4/I");
		NbInfo->Branch("Ncut5", &Ncut5, "Ncut5/I");
		NbInfo->Branch("Ncut6", &Ncut6, "Ncut6/I");
		NbInfo->Branch("Ncut7", &Ncut7, "Ncut7/I");
		NbInfo->Branch("Ncut8", &Ncut8, "Ncut8/I");
		NbInfo->Branch("Ncut9", &Ncut9, "Ncut9/I");
	}

	GammaAll->BypassStreamer();
	pep->BypassStreamer();
	pem->BypassStreamer();
	ppip->BypassStreamer();
	ppim->BypassStreamer();
	pgamma->BypassStreamer();
	petap->BypassStreamer();
	pee->BypassStreamer();
	
	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	std::cout<<"success initialize completed!"<<std::endl;
	return StatusCode::SUCCESS;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode GamPiPiEE::execute() {


	//std::cout << "execute()" << std::endl;

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	int runNo=eventHeader->runNumber();
	int event=eventHeader->eventNumber();
	runid=runNo;
	evtid=event;
	nGamma=0;
	nevt=0;
	isSig = -1;
	m_rconv= -1.;

	conexc_flag1 =200;
	conexc_flag2 =200;

	Ncut0++;  //total events
	log << MSG::DEBUG <<"run, evtnum = "<< runNo << " , "<< event <<endreq;

	if(!(Ncut0%10000))
	{
		cout<<"Processing "<<Ncut0<<"th event:   "<<" Run Id = "<<runNo<<", Event Id = "<<event<<endl;
		cout<<"Total =  "<<Ncut0<<", after good charged tracks cut = "<<Ncut1<<endl;
	}
	//*************Global Event Parameters************
	//do not change below
	int psi;
	if(m_psiType == 0){
		m_Ecms=Epsip;    //psip
		psi=100443;
		m_EnergySpread=0.0013;
	}

	if(m_psiType == 1){
		m_Ecms=Ejpsi;    //Jpsi
		psi=443;
		m_EnergySpread=0.0008;
	}

	if(m_psiType == 2){
		m_Ecms=Econt;    //continuum
		psi=30443;
		m_EnergySpread=0.00097;         
		//	psi=100443;
		//	m_EnergySpread=0.0013;
	}

	HepLorentzVector p4_cms(0.011*m_Ecms,0.,0.,m_Ecms);

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);

	//active simplePID service
	//	ISimplePIDSvc* m_simplePIDSvc;
	//	Gaudi::svcLocator()->service("SimplePIDSvc",m_simplePIDSvc);

	nCHARGED = evtRecEvent->totalCharged();
	nNEUTRAL = evtRecEvent->totalNeutral();
	nTRACKS = evtRecEvent->totalTracks();

	log << MSG::DEBUG <<"ncharg, nneu, tottks = "<< nCHARGED << " , "<< nNEUTRAL << " , "<< nTRACKS <<endreq;
	//*************************************************************MC Info.***************************************************************************

	if (eventHeader->runNumber()<0)
	{
		conexc_flag1 = eventHeader->flag1();
		conexc_flag2 = eventHeader->flag2();
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
		if (!mcParticleCol)
		{
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}
		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		Event::McParticle temp;
		m_idxmc=0;

		bool jpsiDecay = false;
		int jpsiIndex = -1;
		bool strange = false;

		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			if((*iter_mc)->primaryParticle()&&(*iter_mc)->particleProperty()==11&&((*iter_mc)->mother()).particleProperty()==11) {strange = true;}
			if ((*iter_mc)->primaryParticle()) continue;
			if (!(*iter_mc)->decayFromGenerator()) continue;
			if ((*iter_mc)->particleProperty() == 100443) 
			{
				jpsiDecay = true;
				jpsiIndex = (*iter_mc)->trackIndex();
			}
			if (!jpsiDecay) continue;
			int mcidx = ((*iter_mc)->mother()).trackIndex() - jpsiIndex;
			int pdgid = (*iter_mc)->particleProperty();
			if(strange &&((*iter_mc)->mother()).particleProperty()!=100443) mcidx--; 
			//		if(pdgid == -22) continue;
			m_pdgid[m_idxmc] = pdgid;
			m_motheridx[m_idxmc] = mcidx;
			temp=(*iter_mc)->mother();
			m_motherpid[m_idxmc] = temp.particleProperty();
			if(pdgid == 100443) m_drank[m_idxmc]=0;
			else
			{
				for(int i=1;i<100;i++)
				{
					if(temp.particleProperty()==100443)
					{
						m_drank[m_idxmc]=i;
						break;
					}
					temp = temp.mother();
				}
			}

			m_idxmc++;    
		}

		if(m_saveMCTruth == 1)
		{
			xpip->SetE(-1);
			xpim->SetE(-1);
			xep->SetE(-1);
			xem->SetE(-1);
			xgamma->SetE(-1);
			xetap->SetE(-1);
			xee->SetE(-1);

			int numep = 0;
			int numem = 0;
			int numpip = 0;
			int numpim = 0;
			int numgam = 0;
			int numetap = 0;
			int etapIndex = -1;
			int numDaughters = 0;

			if(jpsiIndex > 0)
			{
				for(iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
				{
					if((*iter_mc)->primaryParticle())
					{
						continue;
					}
					if(!(*iter_mc)->decayFromGenerator()) 
					{
						continue;
					}
					if(((*iter_mc)->mother()).trackIndex() !=jpsiIndex)
					{
						continue;
					}
					int pid = (*iter_mc)->particleProperty();
					if(pid == -22) 	continue;

					HepLorentzVector p = (*iter_mc)->initialFourMomentum();
					if(pid == 331)
					{
						xetap->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());
						etapIndex = (*iter_mc)->trackIndex();
						numetap++;
						numDaughters++;
					}
					else if( pid == -11 )
					{
						xep->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());
						numep++;
						numDaughters++;
					}
					else if( pid == 11 )
					{
						xem->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());
						numem++;
						numDaughters++;
					}
					else continue;
				}
			}

			if(etapIndex > 0 )
			{
				for(iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
				{
					if((*iter_mc)->primaryParticle())continue;
					if(!(*iter_mc)->decayFromGenerator()) continue;
					int pid = (*iter_mc)->particleProperty();
					if(pid == -22) continue;
					HepLorentzVector p = (*iter_mc)->initialFourMomentum();
					if(((*iter_mc)->mother()).trackIndex() == etapIndex)
					{
						if(pid == 211)
						{
							xpip->SetPxPyPzE(p.px(), p.py(), p.pz(), p.e());
							numpip++;
							numDaughters++;
						}
						else if(pid == -211)
						{
							xpim->SetPxPyPzE(p.px(), p.py(), p.pz(), p.e());
							numpim++;
							numDaughters++;
						}
						else if(pid == 22)
						{
							xgamma->SetPxPyPzE(p.px(), p.py(), p.pz(), p.e());
							numgam++;
							numDaughters++;
						}
					}
				}
			}
			xee->SetPxPyPzE((*xep+*xem).Px(),(*xep+*xem).Py(),(*xep+*xem).Pz(),(*xep+*xem).E());
			if((xpip->E()>0)&&(xpim->E()>0)&&(xep->E()>0)&&(xem->E()>0)&&(xgamma->E()>0)&&(numpip==1)&&(numpim==1) &&(numep==1)&&(numem==1)&&(numgam==1)) isSig=1;
		}
	}
	if(m_saveTopoTree == 1) TopoTree->Fill();


	//************************************************************************************************************************************************

	//*************************************************************Global Event Parameters************************************************************
	Vint iGood;
	iGood.clear();

	Vint iChrgp,iChrgn;
	iChrgp.clear();
	iChrgn.clear();
	int nChrgp = 0,nChrgn = 0;
	int nGood = 0;
	int nCharge = 0;
	//************************************************************************************************************************************************

	//*************************************************************Primary Vertex*********************************************************************

	Hep3Vector xorigin(0,0,0);
	HepPoint3D vx_avg(0,0,0);
	HepSymMatrix Evx_avg(3,0);

	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid())
	{
		double* dbv = vtxsvc->PrimaryVertex(); 
		double*  vv = vtxsvc->SigmaPrimaryVertex();  
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);
		vx_avg.setX(dbv[0]);
		vx_avg.setY(dbv[1]);
		vx_avg.setZ(dbv[2]);

		Evx_avg[0][0] = vv[0]*vv[0];
		Evx_avg[1][1] = vv[0]*vv[0];
		Evx_avg[2][2] = vv[2]*vv[2];

		vx[0] = dbv[0];
		vx[1] = dbv[1];
		vx[2] = dbv[2];
		Evx[0] = vv[0];
		Evx[1] = vv[1];
		Evx[2] = vv[2];
	}
	VertexParameter vxpar_init;
	vxpar_init.setVx(vx_avg);
	vxpar_init.setEvx(Evx_avg);
	//************************************************************************************************************************************************

	//*************************************************************Good charged track*****************************************************************
	//cout<<"Start good track selection"<<endl;
	for(int i = 0; i < evtRecEvent->totalCharged(); i++)
	{
		if(i>= evtRecTrkCol->size()) break;
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if(!(*itTrk)->isMdcKalTrackValid()) continue;

		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();

		if(mdcKalTrk->charge()==0) continue;

		double theta = mdcTrk->theta();
		if(fabs(cos(theta)) > m_trkAngCut) continue;

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
		VFHelix helixip(point0,a,Ea); 
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0 = fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0 = vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0 = vecipa[1];
		if(fabs(Rvxy0) > m_vr0cut) continue;
		if(fabs(Rvz0) > m_vz0cut) continue;
		iGood.push_back(i);
		nCharge += mdcKalTrk->charge();

		if(mdcKalTrk->charge() == 1) iChrgp.push_back(i);
		else if(mdcKalTrk->charge() == -1) iChrgn.push_back(i);
		else continue;
	}
	nGood = iGood.size();
	nChrgp = iChrgp.size();
	nChrgn = iChrgn.size();

	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;
	if((nGood != 4) || (nCharge!=0) || (nChrgp != 2) || (nChrgn != 2) )
	{
		return StatusCode::SUCCESS;
	}
	Ncut1++;		//after good track selection
	//cout<<"Finish good track selection"<<endl;
	//************************************************************************************************************************************************

	//***********************************************************Good Gamma Selection*****************************************************************
	showerde=0.;
	Vint iGamma;
	Vp4 pGamma;
	iGamma.clear();
	pGamma.clear();
	int idx_gamma=0;
	for(int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) 
	{
		if(i >= evtRecTrkCol->size()) break;
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isEmcShowerValid()) continue;
		RecEmcShower *emcTrk = (*itTrk)->emcShower();
		Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
		
		double eraw=emcTrk->energy();
		double the=emcTrk->theta();
		double time=emcTrk->time();
		
		if(time < TDC_min||time > TDC_max) continue;//time cut
		showerde += eraw;
		double e_threshold=100.0;
		if(fabs(cos(the))<barcut)	e_threshold=m_barrelEnergyThreshold;
		else if((fabs(cos(the))>endmin) && (fabs(cos(the)) < endmax)) e_threshold=m_endEnergyThreshold;
		else continue;
		if(eraw < e_threshold) continue;

		// find the nearest charged track
		double dang = 200.; 
		if(gamiso==1)
		{
			for(int j = 0; j < evtRecEvent->totalCharged(); j++)
			{
				if(j>=evtRecTrkCol->size()) break;
				EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;

				if(!(*jtTrk)->isExtTrackValid()) continue;
				RecExtTrack *extTrk = (*jtTrk)->extTrack();

				if(extTrk->emcVolumeNumber() == -1) continue;
				Hep3Vector extpos = extTrk->emcPosition();//-xorigin; vertex correction!maybe effect is small
				double angd = extpos.angle(emcpos);
				if(angd < dang) dang = angd;
			}
			if(dang>=200) continue;
			dang = dang * 180 / (CLHEP::pi);
		}	
		if(fabs(dang) < m_gammaAngleCut) continue;
		iGamma.push_back(i);
		costheta_gamma[idx_gamma] = cos(the);
		energy_gamma[idx_gamma] = eraw;
		TDCtime[idx_gamma] = time;
		isoAngle[idx_gamma] = dang;
		idx_gamma++;
	}

	nGamma = iGamma.size();
	log << MSG::DEBUG << "num Good Photon " << nGamma  << " , " <<evtRecEvent->totalNeutral()<<endreq;

	if(nGamma < m_NGamma) return StatusCode::SUCCESS; 
	Ncut2++;	//after good gamma cut
	//***************************Finish Good gamma selection***********************************************************************************************
	//***************************Assign 4-momentum to each photon**************************************************************
	TLorentzVector p_JpsiGamma;
	for(int i = 0; i < nGamma; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGamma[i];
		RecEmcShower* emcTrk = (*itTrk)->emcShower();

		Hep3Vector gammaDirection(emcTrk->x() - xorigin.x(), emcTrk->y() - xorigin.y(), emcTrk->z() - xorigin.z());//modifying gamma track
		double eraw = emcTrk->energy();
		double phi = gammaDirection.phi();
		double the = gammaDirection.theta();

		HepLorentzVector ptrk;
		ptrk.setPx(eraw*sin(the)*cos(phi));
		ptrk.setPy(eraw*sin(the)*sin(phi));
		ptrk.setPz(eraw*cos(the));
		ptrk.setE(eraw);

		p_JpsiGamma.SetPxPyPzE(ptrk.px(), ptrk.py(), ptrk.pz(), ptrk.e());
		new((*GammaAll)[i]) TLorentzVector(p_JpsiGamma);
		pGamma.push_back(ptrk);
	}
	//*************************************************************************************************************************
	//***************************Begin PID***********************************************************************************************
	ParticleID *pid = ParticleID::instance();
	for(int i = 0; i < nChrgp; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iChrgp[i];
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()); // use PID sub-system
		pid->identify( pid->onlyElectron() | pid->onlyPion() | pid->onlyKaon() );    // seperater Pion/Kaon/Electron
		pid->calculate();
		if(!(pid->IsPidInfoValid()))
		{
			chiDeDx[i][0] = -99;
			chiDeDx[i][1] = -99;
			chiDeDx[i][2] = -99;

			chiTof[i][0] = -99;
			chiTof[i][1] = -99;
			chiTof[i][2] = -99;

			chiTof1[i][0] = -99;
			chiTof1[i][1] = -99;
			chiTof1[i][2] = -99;

			chiTof2[i][0] = -99;
			chiTof2[i][1] = -99;
			chiTof2[i][2] = -99;

			chisq_pid[i][0] = -99;
			chisq_pid[i][1] = -99;
			chisq_pid[i][2] = -99;

			prob[i][0] = -99;
			prob[i][1] = -99;
			prob[i][2] = -99;
		}
		else
		{
			chiDeDx[i][0] = pid->chiDedx(0);
			chiDeDx[i][1] = pid->chiDedx(2);
			chiDeDx[i][2] = pid->chiDedx(3);

			chiTof[i][0] = pid->chiTof(0);
			chiTof[i][1] = pid->chiTof(2);
			chiTof[i][2] = pid->chiTof(3);

			chiTof1[i][0] = pid->chiTof1(0);
			chiTof1[i][1] = pid->chiTof1(2);
			chiTof1[i][2] = pid->chiTof1(3);

			chiTof2[i][0] = pid->chiTof2(0);
			chiTof2[i][1] = pid->chiTof2(2);
			chiTof2[i][2] = pid->chiTof2(3);

			chisq_pid[i][0] = chiDeDx[i][0]*chiDeDx[i][0] + chiTof1[i][0]*chiTof1[i][0] + chiTof2[i][0]*chiTof2[i][0];
			chisq_pid[i][1] = chiDeDx[i][1]*chiDeDx[i][1] + chiTof1[i][1]*chiTof1[i][1] + chiTof2[i][1]*chiTof2[i][1];
			chisq_pid[i][2] = chiDeDx[i][2]*chiDeDx[i][2] + chiTof1[i][2]*chiTof1[i][2] + chiTof2[i][2]*chiTof2[i][2];

			prob[i][0] = pid->probElectron();
			prob[i][1] = pid->probPion();
			prob[i][2] = pid->probKaon();
		}
	}

	for(int i = 0; i < nChrgn; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iChrgn[i];
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()); // use PID sub-system
		pid->identify( pid->onlyElectron() | pid->onlyPion() | pid->onlyKaon() );    // seperater Pion/Kaon/Electron
		pid->calculate();
		if(!(pid->IsPidInfoValid()))
		{
			chiDeDx[i+nChrgp][0] = -99;
			chiDeDx[i+nChrgp][1] = -99;
			chiDeDx[i+nChrgp][2] = -99;

			chiTof[i+nChrgp][0] = -99;
			chiTof[i+nChrgp][1] = -99;
			chiTof[i+nChrgp][2] = -99;

			chiTof1[i+nChrgp][0] = -99;
			chiTof1[i+nChrgp][1] = -99;
			chiTof1[i+nChrgp][2] = -99;

			chiTof2[i+nChrgp][0] = -99;
			chiTof2[i+nChrgp][1] = -99;
			chiTof2[i+nChrgp][2] = -99;

			chisq_pid[i+nChrgp][0] = -99;
			chisq_pid[i+nChrgp][1] = -99;
			chisq_pid[i+nChrgp][2] = -99;

			prob[i+nChrgp][0] = -99;
			prob[i+nChrgp][1] = -99;
			prob[i+nChrgp][2] = -99;
		}
		else
		{
			chiDeDx[i+nChrgp][0] = pid->chiDedx(0);
			chiDeDx[i+nChrgp][1] = pid->chiDedx(2);
			chiDeDx[i+nChrgp][2] = pid->chiDedx(3);

			chiTof[i+nChrgp][0] = pid->chiTof(0);
			chiTof[i+nChrgp][1] = pid->chiTof(2);
			chiTof[i+nChrgp][2] = pid->chiTof(3);

			chiTof1[i+nChrgp][0] = pid->chiTof1(0);
			chiTof1[i+nChrgp][1] = pid->chiTof1(2);
			chiTof1[i+nChrgp][2] = pid->chiTof1(3);

			chiTof2[i+nChrgp][0] = pid->chiTof2(0);
			chiTof2[i+nChrgp][1] = pid->chiTof2(2);
			chiTof2[i+nChrgp][2] = pid->chiTof2(3);

			chisq_pid[i+nChrgp][0] = chiDeDx[i+nChrgp][0]*chiDeDx[i+nChrgp][0] + chiTof1[i+nChrgp][0]*chiTof1[i+nChrgp][0] + chiTof2[i+nChrgp][0]*chiTof2[i+nChrgp][0];
			chisq_pid[i+nChrgp][1] = chiDeDx[i+nChrgp][1]*chiDeDx[i+nChrgp][1] + chiTof1[i+nChrgp][1]*chiTof1[i+nChrgp][1] + chiTof2[i+nChrgp][1]*chiTof2[i+nChrgp][1];
			chisq_pid[i+nChrgp][2] = chiDeDx[i+nChrgp][2]*chiDeDx[i+nChrgp][2] + chiTof1[i+nChrgp][2]*chiTof1[i+nChrgp][2] + chiTof2[i+nChrgp][2]*chiTof2[i+nChrgp][2];

			prob[i+nChrgp][0] = pid->probElectron();
			prob[i+nChrgp][1] = pid->probPion();
			prob[i+nChrgp][2] = pid->probKaon();
		}
	}
	//*******************************End  PID ***********************************************************
	//*******************************Begin Assign track  info***********************************************************
	WTrackParameter wvTrkp[2],wvTrkm[2];
	Vint iep,iem;
	iep.clear();
	iem.clear();

	RecMdcKalTrack* leppTrack;
	RecMdcKalTrack* lepmTrack;	
	flag_ep = -1;
	flag_pip = -1;
	flag_em = -1;
	flag_pim = -1;

	for(int i = 0; i < nChrgp; i++)
	{
		int i_tag = -1;

		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iChrgp[i];
		RecMdcTrack *imdcTrk = (*itTrk)->mdcTrack();
		RecMdcKalTrack *imdcKalTrk = (*itTrk)->mdcKalTrack();

		HepVector ia = imdcTrk->helix();
		HepSymMatrix iEa = imdcTrk->err();

		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 

		VFHelix ihelixip(point0,ia,iEa); 
		ihelixip.pivot(IP);
		HepVector ivecipa = ihelixip.a();
		double  Rvxy0=fabs(ivecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=ivecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=ivecipa[1];

		if((prob[i][0] > prob[i][1]) && (prob[i][0]>prob[i][2]))
		{
			flag_ep = i;
			i_tag = 0;
			leppTrack = imdcKalTrk;
			iep.push_back(i);
		}	
		else
		{
			flag_pip = i;
			i_tag =1;
		}
		costheta_chrgd[i_tag]=cos(imdcTrk->theta());
		Rxy[i_tag]=imdcTrk->r();
		Rz[i_tag]=imdcTrk->z();
		Rvxy[i_tag]=Rvxy0;
		Rvz[i_tag]=Rvz0;	

		if((*itTrk)->isEmcShowerValid()) 
		{
			RecEmcShower *emcTrk=(*itTrk)->emcShower();
			m_deTrk[i_tag]= emcTrk->energy();
			m_hitsTrk[i_tag]=emcTrk->numHits();
			m_pTrk[i_tag]=imdcKalTrk->p();
			m_epTrk[i_tag]=m_deTrk[i_tag]/m_pTrk[i_tag];
		}
		else
		{
			m_deTrk[i_tag]= -1;
			m_hitsTrk[i_tag]= -1;
			m_pTrk[i_tag]=imdcKalTrk->p();
			m_epTrk[i_tag]= -1;
		}
		if(i_tag==0) wvTrkp[0] = WTrackParameter(me,imdcKalTrk->getZHelixE(),imdcKalTrk->getZErrorE());
		else if(i_tag==1) wvTrkp[1] = WTrackParameter(mpi,imdcKalTrk->getZHelix(),imdcKalTrk->getZError());
		else return StatusCode::SUCCESS;
	}

	for(int i = 0; i < nChrgn; i++)
	{
		int i_tag = -1;

		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iChrgn[i];
		RecMdcTrack *imdcTrk = (*itTrk)->mdcTrack();
		RecMdcKalTrack *imdcKalTrk = (*itTrk)->mdcKalTrack();

		HepVector ia = imdcTrk->helix();
		HepSymMatrix iEa = imdcTrk->err();

		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 

		VFHelix ihelixip(point0,ia,iEa); 
		ihelixip.pivot(IP);
		HepVector ivecipa = ihelixip.a();
		double  Rvxy0=fabs(ivecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=ivecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=ivecipa[1];

		if((prob[i+nChrgp][0] > prob[i+nChrgp][1]) && (prob[i+nChrgp][0]>prob[i+nChrgp][2])) 
		{
			i_tag=0;
			flag_em = i;
			lepmTrack = imdcKalTrk;
			iem.push_back(i);
		}	
		else
		{
			flag_pim = i;
			i_tag =1;
		}
		costheta_chrgd[i_tag+nChrgp]=cos(imdcTrk->theta());
		Rxy[i_tag+nChrgp]=imdcTrk->r();
		Rz[i_tag+nChrgp]=imdcTrk->z();
		Rvxy[i_tag+nChrgp]=Rvxy0;
		Rvz[i_tag+nChrgp]=Rvz0;	

		if((*itTrk)->isEmcShowerValid()) 
		{
			RecEmcShower *emcTrk=(*itTrk)->emcShower();
			m_deTrk[i_tag+nChrgp]= emcTrk->energy();
			m_hitsTrk[i_tag+nChrgp]=emcTrk->numHits();
			m_pTrk[i_tag+nChrgp]=imdcKalTrk->p();
			m_epTrk[i_tag+nChrgp]=m_deTrk[i_tag+nChrgp]/m_pTrk[i_tag+nChrgp];
		}
		else
		{
			m_deTrk[i_tag+nChrgp]= -1;
			m_hitsTrk[i_tag+nChrgp]= -1;
			m_pTrk[i_tag+nChrgp]=imdcKalTrk->p();
			m_epTrk[i_tag+nChrgp]= -1;
		}
		if(i_tag==0) wvTrkm[0] = WTrackParameter(me,imdcKalTrk->getZHelixE(),imdcKalTrk->getZErrorE());
		else if(i_tag==1) wvTrkm[1] = WTrackParameter(mpi,imdcKalTrk->getZHelix(),imdcKalTrk->getZError());
		else return StatusCode::SUCCESS;
	}
	int Nep=iep.size();
	int Nem=iem.size();

	if(Nep != 1 || Nem != 1) return StatusCode::SUCCESS;
	Ncut3++;   //finish PID
	//*******************************End Assign track  info***********************************************************
	//*******************************Begin Gamma Conversion***********************************************************

	HepPoint3D point1(0.,0.,0.);   // the initial point for MDC recosntruction

	RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
	GammaConv gconv = GammaConv(leppTrack->helix(), lepmTrack->helix(),point1);
	m_xconv1 = gconv.getRX1();
	m_yconv1 = gconv.getRY1();
	m_zconv1 = gconv.getRZ1();
	m_rconv1 = gconv.getRXY1();
	m_xconv2 = gconv.getRX2();
	m_yconv2 = gconv.getRY2();
	m_zconv2 = gconv.getRZ2();
	m_rconv2 = gconv.getRXY2();

	m_xiep  = gconv.getXiep();
	m_deltaxy =gconv.getDeltaXY();

	HepLorentzVector gammaVector(0,0,0,0);
	gammaVector=HepLorentzVector(leppTrack->p4(mel)+lepmTrack->p4(mel));
	m_cthep = (leppTrack->p4(mel)).vect()*(lepmTrack->p4(mel)).vect()/leppTrack->p()/lepmTrack->p();
	m_ptrkp  = leppTrack->p();
	m_ptrkm  = lepmTrack->p();
	m_mgamma = gammaVector.m();
	m_egamma = gammaVector.e();
	m_theta        = gammaVector.theta();
	m_cosTheta     = gammaVector.cosTheta();
	m_phi          = gammaVector.phi();
	//Gamma Conv
	//*******************************End   Gamma Conversion****************************************************
	//*******************************Begin Vetex Fit***********************************************************
	WTrackParameter wTrkp[2];
	WTrackParameter wTrkm[2];
	HepPoint3D vx_infit;
	HepSymMatrix Evx_infit;
	HepPoint3D vx(0., 0., 0.);
	HepSymMatrix Evx(3, 0);
	double bx = 1E+6;
	double by = 1E+6;
	double bz = 1E+6;
	Evx[0][0] = bx*bx;
	Evx[1][1] = by*by;
	Evx[2][2] = bz*bz;
	VertexParameter vxpar;
	vxpar.setVx(vx);
	vxpar.setEvx(Evx);
	VertexFit* vtxfit = VertexFit::instance();
	vtxfit->init();
	vtxfit->AddTrack(0,  wvTrkp[0]);
	vtxfit->AddTrack(1,  wvTrkm[0]);
	vtxfit->AddTrack(2,  wvTrkp[1]);
	vtxfit->AddTrack(3,  wvTrkm[1]);
	vtxfit->AddVertex(0, vxpar,0,1,2,3);

	if(!vtxfit->Fit(0)) return StatusCode::SUCCESS;
	vtxfit->Swim(0);

	vtxchisq=vtxfit->chisq(0);
	wTrkp[0]=vtxfit->wtrk(0);
	wTrkm[0]=vtxfit->wtrk(1);
	wTrkp[1]=vtxfit->wtrk(2);
	wTrkm[1]=vtxfit->wtrk(3);

	vx_infit=vtxfit->vx(0);
	Evx_infit=vtxfit->Evx(0);

	Ncut4++;   // finish vertex fit
	//********************************************************************************************************************
	vtxfit->init();
	vtxfit->AddTrack(0,  wvTrkp[0]);
	vtxfit->AddTrack(1,  wvTrkm[0]);
	vtxfit->AddVertex(0, vxpar,0,1);
	vtxfit->Swim(0);
	if(!vtxfit->Fit(0)) vtxchisqee=-1;
	else vtxchisqee=vtxfit->chisq(0);

	vtxfit->init();
	vtxfit->AddTrack(0,  wvTrkp[1]);
	vtxfit->AddTrack(1,  wvTrkm[1]);
	vtxfit->AddVertex(0, vxpar,0,1);
	vtxfit->Swim(0);
	if(!vtxfit->Fit(0)) vtxchisqpipi=-1;
	else vtxchisqpipi=vtxfit->chisq(0);
	//*******************************End Vetex Fit***********************************************************
	//*******************************Begin 4C  Fit after Vtxfit***********************************************************
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	TLorentzVector p_jpsi(p4_cms.px(),p4_cms.py(),p4_cms.pz(),p4_cms.e());	
	
	int gam_min4C=-1;
	double temp_4Cchisq = 99999.;
	min4C_index = -1;
	for(int i=0; i < nGamma ; i++)
	{
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGamma[i];
		RecEmcShower* gammaTrk = (*itTrk)->emcShower();
		kmfit->init();
		kmfit->setEspread(m_EnergySpread);
	//	kmfit->setBeamPosition(vx_infit);
		kmfit->setBeamPosition(xorigin);
	//	kmfit->setVBeamPosition(Evx_infit);
		kmfit->setVBeamPosition(Evx_avg);
		kmfit->AddTrack(0,wTrkp[0]);
		kmfit->AddTrack(1,wTrkm[0]);
		kmfit->AddTrack(2,wTrkp[1]);
		kmfit->AddTrack(3,wTrkm[1]);
		kmfit->AddTrack(4,0.0,gammaTrk);
		kmfit->AddFourMomentum(0,p4_cms);
		bool okfit1 = kmfit->Fit();
		if(!okfit1) continue;
		chi4C[nevt] = kmfit->chisq();

		HepLorentzVector ep_ = kmfit->pfit(0);
		HepLorentzVector em_ = kmfit->pfit(1);
		HepLorentzVector pip_ = kmfit->pfit(2);
		HepLorentzVector pim_ = kmfit->pfit(3);
		HepLorentzVector gamma_ = kmfit->pfit(4);

		TLorentzVector p_ep,p_em,p_pip,p_pim,p_gamma,p_etap,p_ee;
		p_ep.SetPxPyPzE(ep_.px(), ep_.py(), ep_.pz(), ep_.e());
		p_em.SetPxPyPzE(em_.px(), em_.py(), em_.pz(), em_.e());
		p_pip.SetPxPyPzE(pip_.px(), pip_.py(), pip_.pz(), pip_.e());
		p_pim.SetPxPyPzE(pim_.px(), pim_.py(), pim_.pz(), pim_.e());
		p_gamma.SetPxPyPzE(gamma_.px(), gamma_.py(), gamma_.pz(), gamma_.e());
		p_etap.SetPxPyPzE((p_pip+p_pim+p_gamma).Px(),(p_pip+p_pim+p_gamma).Py(),(p_pip+p_pim+p_gamma).Pz(),(p_pip+p_pim+p_gamma).E());
		p_ee.SetPxPyPzE((p_ep+p_em).Px(),(p_ep+p_em).Py(),(p_ep+p_em).Pz(),(p_ep+p_em).E());

		new ((*pep)[nevt]) TLorentzVector(p_ep);
		new ((*pem)[nevt]) TLorentzVector(p_em);
		new ((*ppip)[nevt]) TLorentzVector(p_pip);
		new ((*ppim)[nevt]) TLorentzVector(p_pim);
		new ((*pgamma)[nevt]) TLorentzVector(p_gamma);
		new ((*petap)[nevt]) TLorentzVector(p_etap);
		new ((*pee)[nevt]) TLorentzVector(p_ee);


		if(chi4C[nevt] < temp_4Cchisq )
		{
			gam_min4C=i;
			min4C_index=nevt;
			temp_4Cchisq=chi4C[nevt];
			min_chi4C = chi4C[nevt];

			min_ep->SetPxPyPzE(p_ep.Px(),p_ep.Py(),p_ep.Pz(),p_ep.E());
			min_em->SetPxPyPzE(p_em.Px(),p_em.Py(),p_em.Pz(),p_em.E());
			min_pip->SetPxPyPzE(p_pip.Px(),p_pip.Py(),p_pip.Pz(),p_pip.E());
			min_pim->SetPxPyPzE(p_pim.Px(),p_pim.Py(),p_pim.Pz(),p_pim.E());
			min_gamma->SetPxPyPzE(p_gamma.Px(),p_gamma.Py(),p_gamma.Pz(),p_gamma.E());
			min_etap->SetPxPyPzE(p_etap.Px(),p_etap.Py(),p_etap.Pz(),p_etap.E());
			min_ee->SetPxPyPzE(p_ee.Px(),p_ee.Py(),p_ee.Pz(),p_ee.E());
			
			min_ep_unfit->SetPxPyPzE(p_ep.Px(),p_ep.Py(),p_ep.Pz(),p_ep.E());
			min_em_unfit->SetPxPyPzE(p_em.Px(),p_em.Py(),p_em.Pz(),p_em.E());
			min_pip_unfit->SetPxPyPzE(p_pip.Px(),p_pip.Py(),p_pip.Pz(),p_pip.E());
			min_pim_unfit->SetPxPyPzE(p_pim.Px(),p_pim.Py(),p_pim.Pz(),p_pim.E());
			min_gamma_unfit->SetPxPyPzE(p_gamma.Px(),p_gamma.Py(),p_gamma.Pz(),p_gamma.E());
			min_etap_unfit->SetPxPyPzE(p_etap.Px(),p_etap.Py(),p_etap.Pz(),p_etap.E());
			min_ee_unfit->SetPxPyPzE(p_ee.Px(),p_ee.Py(),p_ee.Pz(),p_ee.E());
		}
		nevt++;
	}	
	//	cout<<"nevt:	"<<nevt<<endl;
	//*******************************End 4C  Fit after Vtxfit***********************************************************
	//*******************************Begin 4C  Fit without Vtxfit*******************************************************
	for(int i=0; i < nGamma ; i++)
	{
		if( i != gam_min4C) continue;
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGamma[i];
		RecEmcShower* gammaTrk = (*itTrk)->emcShower();
		kmfit->init();
		kmfit->setEspread(m_EnergySpread);
		kmfit->setBeamPosition(xorigin);
		kmfit->setVBeamPosition(Evx_avg);
		kmfit->AddTrack(0,wvTrkp[0]);
		kmfit->AddTrack(1,wvTrkm[0]);
		kmfit->AddTrack(2,wvTrkp[1]);
		kmfit->AddTrack(3,wvTrkm[1]);
		kmfit->AddTrack(4,0.0,gammaTrk);
		kmfit->AddFourMomentum(0,p4_cms);
		bool okfit2 = kmfit->Fit();
		if(!okfit2) min_chi4C_unfit=-1;
		else
		{
			min_chi4C_unfit=kmfit->chisq();
			
			HepLorentzVector ep_unfit = kmfit->pfit(0);
			HepLorentzVector em_unfit = kmfit->pfit(1);
			HepLorentzVector pip_unfit = kmfit->pfit(2);
			HepLorentzVector pim_unfit = kmfit->pfit(3);
			HepLorentzVector gamma_unfit = kmfit->pfit(4);
			
			TLorentzVector p_ep_unfit,p_em_unfit,p_pip_unfit,p_pim_unfit,p_gamma_unfit,p_etap_unfit,p_ee_unfit;
			p_gamma_unfit.SetPxPyPzE((pGamma[i]).px(), (pGamma[i]).py(), (pGamma[i]).pz(), (pGamma[i]).e());
			p_ep_unfit.SetPxPyPzE((wTrkp[0].p()).px(), (wTrkp[0].p()).py(), (wTrkp[0].p()).pz(), (wTrkp[0].p()).e());
			p_em_unfit.SetPxPyPzE((wTrkm[0].p()).px(), (wTrkm[0].p()).py(), (wTrkm[0].p()).pz(), (wTrkm[0].p()).e());
			p_pip_unfit.SetPxPyPzE((wTrkp[1].p()).px(), (wTrkp[1].p()).py(), (wTrkp[1].p()).pz(), (wTrkp[1].p()).e());
			p_pim_unfit.SetPxPyPzE((wTrkm[1].p()).px(), (wTrkm[1].p()).py(), (wTrkm[1].p()).pz(), (wTrkm[1].p()).e());
			p_etap_unfit.SetPxPyPzE((p_pip_unfit+p_pim_unfit+p_gamma_unfit).Px(),(p_pip_unfit+p_pim_unfit+p_gamma_unfit).Py(),(p_pip_unfit+p_pim_unfit+p_gamma_unfit).Pz(),(p_pip_unfit+p_pim_unfit+p_gamma_unfit).E());
			p_ee_unfit.SetPxPyPzE((p_ep_unfit+p_em_unfit).Px(),(p_ep_unfit+p_em_unfit).Py(),(p_ep_unfit+p_em_unfit).Pz(),(p_ep_unfit+p_em_unfit).E());
			
			min_ep_unfit->SetPxPyPzE(p_ep_unfit.Px(),p_ep_unfit.Py(),p_ep_unfit.Pz(),p_ep_unfit.E());
			min_em_unfit->SetPxPyPzE(p_em_unfit.Px(),p_em_unfit.Py(),p_em_unfit.Pz(),p_em_unfit.E());
			min_pip_unfit->SetPxPyPzE(p_pip_unfit.Px(),p_pip_unfit.Py(),p_pip_unfit.Pz(),p_pip_unfit.E());
			min_pim_unfit->SetPxPyPzE(p_pim_unfit.Px(),p_pim_unfit.Py(),p_pim_unfit.Pz(),p_pim_unfit.E());
			min_gamma_unfit->SetPxPyPzE(p_gamma_unfit.Px(),p_gamma_unfit.Py(),p_gamma_unfit.Pz(),p_gamma_unfit.E());
			min_etap_unfit->SetPxPyPzE(p_etap_unfit.Px(),p_etap_unfit.Py(),p_etap_unfit.Pz(),p_etap_unfit.E());
			min_ee_unfit->SetPxPyPzE(p_ee_unfit.Px(),p_ee_unfit.Py(),p_ee_unfit.Pz(),p_ee_unfit.E());
		}

	}	
	//Ncut5++;  //4C
	//*******************************End 4C  Fit without Vtxfit***************************************************
	//*******************************Begin Calculate *************************************************************
	*rec_ee=p_jpsi-*min_ee;	
	*rec_ee_unfit=p_jpsi-*min_ee_unfit;	
	
	angee = min_ep->Angle(min_em->Vect()); //rad
	angee = angee * 180 / (CLHEP::pi); //deg
	angee_unfit = min_ep_unfit->Angle(min_em_unfit->Vect()); //rad
	angee_unfit = angee_unfit * 180 / (CLHEP::pi); //deg

	m_rconv=angee<10?m_rconv2:m_rconv1;	
	cosee=cos(min_ep->Angle(min_em->Vect()));
	mee=min_ee->M();
	mee_unfit=min_ee_unfit->M();
	metap=min_etap->M();
	mgee=(*min_ep+*min_em+*min_gamma).M();
	m_rec_ee=rec_ee->M();
	m_rec_ee_unfit=rec_ee_unfit->M();
	
	mpipi = (*min_pip+*min_pim).M();
	
	m_PtTrk[0]=min_ep->Pt();	
	m_PtTrk[1]=min_pip->Pt();	
	m_PtTrk[2]=min_em->Pt();	
	m_PtTrk[3]=min_pim->Pt();	
	
	TVector3 ee_ref=min_ee->BoostVector();
	TVector3 etap_ref=min_etap->BoostVector();

	hel_ep->SetPxPyPzE(min_ep->Px(),min_ep->Py(),min_ep->Pz(),min_ep->E());
	hel_pip->SetPxPyPzE(min_pip->Px(),min_pip->Py(),min_pip->Pz(),min_pip->E());
	hel_em->SetPxPyPzE(min_em->Px(),min_em->Py(),min_em->Pz(),min_em->E());
	hel_pim->SetPxPyPzE(min_pim->Px(),min_pim->Py(),min_pim->Pz(),min_pim->E());
	hel_gam->SetPxPyPzE(min_gamma->Px(),min_gamma->Py(),min_gamma->Pz(),min_gamma->E());
	
	hel_ep->Boost(ee_ref);		
	hel_pip->Boost(etap_ref);		
	hel_em->Boost(ee_ref);		
	hel_pim->Boost(etap_ref);		
	hel_gam->Boost(etap_ref);		
	
	hel_cos[0] = hel_ep->CosTheta();
	hel_cos[1] = hel_pip->CosTheta();
	hel_cos[2] = hel_em->CosTheta();
	hel_cos[3] = hel_pim->CosTheta();
	hel_cos[4] = hel_gam->CosTheta();

	//*******************************End  Calculate *************************************************************
	if(nevt > 0 )
	{
		TreeAna->Fill();
		Ncut5++;
	}
	GammaAll->Clear();
	pep->Clear();
	pem->Clear();
	ppip->Clear();
	ppim->Clear();
	pgamma->Clear();
	petap->Clear();
	pee->Clear();
	
	xep->Clear();
	xem->Clear();
	xpip->Clear();
	xpim->Clear();
	xgamma->Clear();
	xetap->Clear();
	xee->Clear();

	hel_ep->Clear();
	hel_pip->Clear();
	hel_em->Clear();
	hel_pim->Clear();
	hel_gam->Clear();

	min_ep->Clear();
	min_em->Clear();
	min_pip->Clear();
	min_pim->Clear();
	min_gamma->Clear();
	min_etap->Clear();
	min_ee->Clear();
	rec_ee->Clear();

	min_ep_unfit->Clear();
	min_em_unfit->Clear();
	min_pip_unfit->Clear();
	min_pim_unfit->Clear();
	min_gamma_unfit->Clear();
	min_etap_unfit->Clear();
	min_ee_unfit->Clear();
	rec_ee_unfit->Clear();

	return StatusCode::SUCCESS;
}

//***************************************************************************

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode GamPiPiEE::finalize() {
	NbInfo->Fill();
	saveFile->cd();
	TreeAna->Write();
	if(m_saveTopoTree == 1)	TopoTree->Write();
	if(m_saveCutFlow == 1) NbInfo->Write();
	saveFile->Close();
	cout<<endl<<"Finalize psi'->etapee, etap->gampipi: Version 1"<<endl<<endl;;
	cout<<"Total number:                         "<<Ncut0<<endl;
	cout<<"nGood==4, nCharge==0:                 "<<Ncut1<<endl;
	cout<<"nGamma=0:                             "<<Ncut2<<endl;
	cout<<"PID:                                  "<<Ncut3<<endl;
	cout<<"vertex fit:                           "<<Ncut4<<endl;
	cout<<"4C success:                           "<<Ncut5<<endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
