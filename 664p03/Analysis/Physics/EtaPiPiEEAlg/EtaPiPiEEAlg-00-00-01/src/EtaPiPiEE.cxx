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

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

#include "EtaPiPiEEAlg/EtaPiPiEE.h"

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
const double mk = 0.493677;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double mel  = 0.000510998910;
const double Ejpsi=3.097;
const double Epsip=3.686;
const double Econt=3.08;
//const double Econt=3.773;

//const double velc = 29.9792458;  tof_path unit in cm.
const double velc = 299.792458;   // tof path unit in mm

typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;

int Ncut0 = 0;		//total event
int Ncut1 = 0;		//Pass good charged tracks cut
int Ncut2 = 0;		//Pass good photon
int Ncut3 = 0;		//Pass PID
int Ncut4 = 0;		//Pass Vertex
int Ncut5 = 0;		//Pass 4C
int Ncut6 = 0;		
int Ncut7 = 0;
int Ncut8 = 0;
int Ncut9 = 0;
/////////////////////////////////////////////////////////////////////////////
EtaPiPiEE::EtaPiPiEE(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) {
		//Declare the properties
		declareProperty("EnergySpread",  m_EnergySpread = 0.0013);
		declareProperty("OutputFileName",  m_OutputFileName = "EtaPiPiEE_test.root");
		declareProperty("Ecms",  m_Ecms = 3.686);
		declareProperty("psiType", m_psiType = 0);
		declareProperty("saveTopo", m_saveTopo = 1);//need to be re-evaluated when running different samples(1 for MC)
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
		declareProperty("mc_cor",mc_cor=0);
	}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode EtaPiPiEE::initialize(){
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

	//TClonesArray initialize 
	GammaAll = new TClonesArray("TLorentzVector");

	//Lorentz Vector 
	xpip=new TLorentzVector();
	xpim=new TLorentzVector();
	xem=new TLorentzVector();
	xep=new TLorentzVector();
	xgamma=new TLorentzVector();
	xee=new TLorentzVector();
	xgee=new TLorentzVector();
	xetap=new TLorentzVector();

	p_pip=new TLorentzVector();
	p_pim=new TLorentzVector();
	p_em=new TLorentzVector();
	p_ep=new TLorentzVector();
	p_gamma=new TLorentzVector();
	p_ee=new TLorentzVector();
	p_pipi=new TLorentzVector();
	p_gee=new TLorentzVector();
	p_gpipi=new TLorentzVector();
	p_recpipi=new TLorentzVector();

	p_upip=new TLorentzVector();
	p_upim=new TLorentzVector();
	p_uem=new TLorentzVector();
	p_uep=new TLorentzVector();
	p_ugamma=new TLorentzVector();
	p_uee=new TLorentzVector();
	p_upipi=new TLorentzVector();
	p_ugee=new TLorentzVector();
	p_ugpipi=new TLorentzVector();

	//mc info
	if(m_saveTopo == 1)
	{
		TreeAna->Branch("indexmc", &m_idxmc, "indexmc/I");
		TreeAna->Branch("pdgid", m_pdgid, "pdgid[indexmc]/I");
		TreeAna->Branch("drank", m_drank, "drank[indexmc]/I");
		TreeAna->Branch("motheridx", m_motheridx, "motheridx[indexmc]/I");
		TreeAna->Branch("motherpid", m_motherpid, "motherpid[indexmc]/I");
	}
	if(m_saveMCTruth == 1)
	{
		TreeAna->Branch("issig", &issig, "issig/I");
		TreeAna->Branch("xmee", &xmee, "xmee/D");
		TreeAna->Branch("xpip",&xpip,32000,0);
		TreeAna->Branch("xpim",&xpim,32000,0);
		TreeAna->Branch("xep",&xep,32000,0);
		TreeAna->Branch("xem",&xem,32000,0);
		TreeAna->Branch("xgamma",&xgamma,32000,0);
		TreeAna->Branch("xgee",&xgee,32000,0);
		TreeAna->Branch("xee",&xee,32000,0);
		TreeAna->Branch("xetap",&xetap,32000,0);
	}
	TreeAna->Branch("runid", &runid, "runid/I");
	TreeAna->Branch("evtid", &evtid, "evtid/I");
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
	TreeAna->Branch("pKal", pKal, "pKal[4]/D");
	TreeAna->Branch("pxKal", pxKal, "pxKal[4]/D");
	TreeAna->Branch("pyKal", pyKal, "pyKal[4]/D");
	TreeAna->Branch("pzKal", pzKal, "pzKal[4]/D");
	TreeAna->Branch("deemc", deemc, "deemc[4]/D");
	TreeAna->Branch("eop", eop, "eop[4]/D");

	TreeAna->Branch("flag_pip", &flag_pip, "flag_pip/I");
	TreeAna->Branch("flag_pim", &flag_pim, "flag_pim/I");
	TreeAna->Branch("flag_ep", &flag_ep, "flag_ep/I");
	TreeAna->Branch("flag_em", &flag_em, "flag_em/I");

	TreeAna->Branch("prob", prob, "prob[4][3]/D");

	TreeAna->Branch("vtxchisq", &vtxchisq, "vtxchisq/D");

	TreeAna->Branch("chi4C", &chi4C, "chi4C/D");
	TreeAna->Branch("p_pip",&p_pip,32000,0);
	TreeAna->Branch("p_pim",&p_pim,32000,0);
	TreeAna->Branch("p_ep",&p_ep,32000,0);
	TreeAna->Branch("p_em",&p_em,32000,0);
	TreeAna->Branch("p_gamma",&p_gamma,32000,0);
	TreeAna->Branch("p_ee",&p_ee,32000,0);
	TreeAna->Branch("p_pipi",&p_pipi,32000,0);
	TreeAna->Branch("p_gee",&p_gee,32000,0);
	TreeAna->Branch("p_gpipi",&p_gpipi,32000,0);
	TreeAna->Branch("p_recpipi",&p_recpipi,32000,0);

	TreeAna->Branch("p_upip",&p_upip,32000,0);
	TreeAna->Branch("p_upim",&p_upim,32000,0);
	TreeAna->Branch("p_uep",&p_uep,32000,0);
	TreeAna->Branch("p_uem",&p_uem,32000,0);
	TreeAna->Branch("p_ugamma",&p_ugamma,32000,0);
	TreeAna->Branch("p_uee",&p_uee,32000,0);
	TreeAna->Branch("p_upipi",&p_upipi,32000,0);
	TreeAna->Branch("p_ugee",&p_ugee,32000,0);
	TreeAna->Branch("p_ugpipi",&p_ugpipi,32000,0);

	TreeAna->Branch("m_gee", &m_gee, "m_gee/D");
	TreeAna->Branch("m_gpipi", &m_gpipi, "m_gpipi/D");
	TreeAna->Branch("angee", &angee, "angee/D");
	TreeAna->Branch("m_uee", &m_uee, "m_uee/D");
	TreeAna->Branch("m_ee", &m_ee, "m_ee/D");
	TreeAna->Branch("m_recpipi", &m_recpipi, "m_recpipi/D");
	TreeAna->Branch("highe_eop", &highe_eop, "highe_eop/D");
	TreeAna->Branch("lowe_eop", &lowe_eop, "lowe_eop/D");
	
	TreeAna->Branch("em_pull_0", &em_pull_0, "em_pull_0/D");
	TreeAna->Branch("em_pull_1", &em_pull_1, "em_pull_1/D");
	TreeAna->Branch("em_pull_2", &em_pull_2, "em_pull_2/D");
	TreeAna->Branch("em_pull_3", &em_pull_3, "em_pull_3/D");
	TreeAna->Branch("em_pull_4", &em_pull_4, "em_pull_4/D");

	TreeAna->Branch("ep_pull_0", &ep_pull_0, "ep_pull_0/D");
	TreeAna->Branch("ep_pull_1", &ep_pull_1, "ep_pull_1/D");
	TreeAna->Branch("ep_pull_2", &ep_pull_2, "ep_pull_2/D");
	TreeAna->Branch("ep_pull_3", &ep_pull_3, "ep_pull_3/D");
	TreeAna->Branch("ep_pull_4", &ep_pull_4, "ep_pull_4/D");

	TreeAna->Branch("pip_pull_0", &pip_pull_0, "pip_pull_0/D");
	TreeAna->Branch("pip_pull_1", &pip_pull_1, "pip_pull_1/D");
	TreeAna->Branch("pip_pull_2", &pip_pull_2, "pip_pull_2/D");
	TreeAna->Branch("pip_pull_3", &pip_pull_3, "pip_pull_3/D");
	TreeAna->Branch("pip_pull_4", &pip_pull_4, "pip_pull_4/D");

	TreeAna->Branch("pim_pull_0", &pim_pull_0, "pim_pull_0/D");
	TreeAna->Branch("pim_pull_1", &pim_pull_1, "pim_pull_1/D");
	TreeAna->Branch("pim_pull_2", &pim_pull_2, "pim_pull_2/D");
	TreeAna->Branch("pim_pull_3", &pim_pull_3, "pim_pull_3/D");
	TreeAna->Branch("pim_pull_4", &pim_pull_4, "pim_pull_4/D");


	//Gamma conversion chuxk
	TreeAna->Branch("m_rconv", &m_rconv, "m_rconv/D");
	TreeAna->Branch("m_xconv1", &m_xconv1, "m_xconv1/D");
	TreeAna->Branch("m_yconv1", &m_yconv1, "m_yconv1/D");
	TreeAna->Branch("m_zconv1", &m_zconv1, "m_zconv1/D");
	TreeAna->Branch("m_rconv1", &m_rconv1, "m_rconv1/D");
	TreeAna->Branch("m_xconv2", &m_xconv2, "m_xconv2/D");
	TreeAna->Branch("m_yconv2", &m_yconv2, "m_yconv2/D");
	TreeAna->Branch("m_zconv2", &m_zconv2, "m_zconv2/D");
	TreeAna->Branch("m_rconv2", &m_rconv2, "m_rconv2/D");
	TreeAna->Branch("m_xiep", &m_xiep, "m_xiep/D");
	TreeAna->Branch("m_deltaxy", &m_deltaxy, "m_deltaxy/D");
	TreeAna->Branch("m_thetaeg1",&m_thetaeg1,"m_thetaeg1/D");
	TreeAna->Branch("m_thetaeg2",&m_thetaeg2,"m_thetaeg2/D");
	TreeAna->Branch("m_psipair", &m_psipair, "m_psipair/D");
	TreeAna->Branch("m_deltaz1", &m_deltaz1, "m_deltaz1/D");
	TreeAna->Branch("m_deltaz2", &m_deltaz2, "m_deltaz2/D");
	TreeAna->Branch("m_cthep", &m_cthep, "m_cthep/D");
	TreeAna->Branch("m_ptrkp", &m_ptrkp, "m_ptrkp/D");
	TreeAna->Branch("m_ptrkm", &m_ptrkm, "m_ptrkm/D");
	TreeAna->Branch("m_mgamma", &m_mgamma, "m_mgamma/D");
	TreeAna->Branch("m_egamma", &m_egamma, "m_egamma/D");
	TreeAna->Branch("m_theta", &m_theta, "m_theta/D");
	TreeAna->Branch("m_cosTheta", &m_cosTheta, "m_cosTheta/D");
	TreeAna->Branch("m_rp", &m_rp, "m_rp/D");
	TreeAna->Branch("m_re", &m_re, "m_re/D");
	TreeAna->Branch("m_deltaeq", &m_deltaeq, "m_deltaeq/D");
	TreeAna->Branch("m_case", &m_case, "m_case/D");
	TreeAna->Branch("m_lep", &m_lep, "m_lep/D");
	TreeAna->Branch("m_phi", &m_phi, "m_phi/D");


	TreeAna->Branch("vtxchie", &vtxchie, "vtxchie/D");
	TreeAna->Branch("m_epemx", &m_epemx, "m_epemx/D");
	TreeAna->Branch("m_epemy", &m_epemy, "m_epemy/D");
	TreeAna->Branch("m_epemz", &m_epemz, "m_epemz/D");
	TreeAna->Branch("m_epemxy", &m_epemxy, "m_epemxy/D");
	
	TreeAna->Branch("m_mepr", &m_mepr, "m_mepr/D");
	TreeAna->Branch("m_mepcenterx", &m_mepcenterx, "m_mepcenterx/D");
	TreeAna->Branch("m_mepcentery", &m_mepcentery, "m_mepcentery/D");
	TreeAna->Branch("m_mepcenterz", &m_mepcenterz, "m_mepcenterz/D");
	TreeAna->Branch("m_memr", &m_memr, "m_memr/D");
	TreeAna->Branch("m_memcenterx", &m_memcenterx, "m_memcenterx/D");
	TreeAna->Branch("m_memcentery", &m_memcentery, "m_memcentery/D");
	TreeAna->Branch("m_memcenterz", &m_memcenterz, "m_memcenterz/D");
	
	TreeAna->Branch("m_epemxxorigin1e", &m_epemxxorigin1e, "m_epemxxorigin1e/D");
	TreeAna->Branch("m_epemyxorigin1e", &m_epemyxorigin1e, "m_epemyxorigin1e/D");
	TreeAna->Branch("m_epemzxorigin1e", &m_epemzxorigin1e, "m_epemzxorigin1e/D");
	TreeAna->Branch("m_epemxyxorigin1e", &m_epemxyxorigin1e, "m_epemxyxorigin1e/D");
	TreeAna->Branch("m_rconvz", &m_rconvz, "m_rconvz/D");

	//for number of events passing every cut
	if(m_saveCutFlow == 1)
	{
		NbInfo = new TTree("nbinfo","nbinfo");
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

	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	std::cout<<"success initialize completed!"<<std::endl;
	return StatusCode::SUCCESS;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode EtaPiPiEE::execute() {
	//std::cout << "execute()" << std::endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	int runNo=eventHeader->runNumber();
	int event=eventHeader->eventNumber();
	runid=runNo;
	evtid=event;
	nGamma=0;
	issig=0;

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

	nCHARGED = evtRecEvent->totalCharged();
	nNEUTRAL = evtRecEvent->totalNeutral();
	nTRACKS = evtRecEvent->totalTracks();

	log << MSG::DEBUG <<"ncharg, nneu, tottks = "<< nCHARGED << " , "<< nNEUTRAL << " , "<< nTRACKS <<endreq;
	////////////////////////////////MC Info.
	if (eventHeader->runNumber()<0)
	{
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
			//	if(pdgid == -22) continue;
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
			int numep(0),numem(0),numpip(0),numpim(0),numgam(0),numetap(0),etapIndex(-1);

			if(jpsiIndex > 0){
				for(iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++){
					if((*iter_mc)->primaryParticle()) continue;
					if(!(*iter_mc)->decayFromGenerator()) continue;
					if(((*iter_mc)->mother()).trackIndex() !=jpsiIndex)	continue;
					int pid = (*iter_mc)->particleProperty();
					if(pid == -22) 	continue;
					HepLorentzVector p = (*iter_mc)->initialFourMomentum();

					if(pid == 331){
						xetap->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());etapIndex = (*iter_mc)->trackIndex();numetap++;}
					else if( pid == -11 ){
						xep->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());numep++;}
					else if( pid == 11 ) {
						xem->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());numem++;}
					else continue;
				}
			}
			if(etapIndex > 0 ){
				for(iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++){
					if((*iter_mc)->primaryParticle())continue;
					if(!(*iter_mc)->decayFromGenerator()) continue;
					if(((*iter_mc)->mother()).trackIndex() != etapIndex) continue;
					int pid = (*iter_mc)->particleProperty();
					if(pid == -22) continue;
					HepLorentzVector p = (*iter_mc)->initialFourMomentum();

					if(pid == 211){
						xpip->SetPxPyPzE(p.px(), p.py(), p.pz(), p.e());numpip++;}
					else if(pid == -211){
						xpim->SetPxPyPzE(p.px(), p.py(), p.pz(), p.e());numpim++;}
					else if(pid == 22){
						xgamma->SetPxPyPzE(p.px(), p.py(), p.pz(), p.e());numgam++;}
				}
			}
			xee->SetPxPyPzE((*xep+*xem).Px(),(*xep+*xem).Py(),(*xep+*xem).Pz(),(*xep+*xem).E());
			xmee=xee->M();
			if((xpip->E()>0)&&(xpim->E()>0)&&(xep->E()>0)&&(xem->E()>0)&&(xgamma->E()>0)&&(xetap->E()>0)&&(numpip==1)&&(numpim==1) &&(numep==1)&&(numem==1)&&(numgam==1)&&(numetap==1)) issig=1;
		}
	}

	////////////////////////////////////primary vertex
	Hep3Vector xorigin(0,0,0);
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

	////////////////////////good track
	Vint iGood,iTrp,iTrm;
	iGood.clear();
	iTrp.clear();
	iTrm.clear();
	int nGood = 0;
	int nTrp = 0;
	int nTrm = 0;
	int nCharge = 0;

	HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
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

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		VFHelix helixip(point0,a,Ea); 
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0 = fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0 = vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0 = vecipa[1];

		if(fabs(cos(theta)) > m_trkAngCut) continue;
		if(fabs(Rvxy0) >= m_vr0cut) continue;
		if(fabs(Rvz0) >= m_vz0cut) continue;

		if(mdcKalTrk->charge() > 0)  iTrp.push_back(i);
		if(mdcKalTrk->charge() < 0)  iTrm.push_back(i);
		iGood.push_back(i);
		nCharge += mdcKalTrk->charge();

	}
	nTrp=iTrp.size();
	nTrm=iTrm.size();
	nGood = iGood.size();
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;
	if((nGood != 4) || (nCharge!=0) ) return StatusCode::SUCCESS;
	Ncut1++;		//after good track selection


	//////////////////////////////////////good photon select
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

	/////////////////////////////////assign 4-momentum for photons
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

	////////////////////////////////PID
	flag_ep=-1;flag_em=-1;flag_pip=-1;flag_pim=-1;
	int N_eppid(0),N_empid(0);
	ParticleID *pid = ParticleID::instance();

	for(int i = 0; i < nGood; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		int _charge = mdcKalTrk->charge();
		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()| pid->useTofE()); // use PID sub-system
		pid->identify(pid->onlyElectron() | pid->onlyPion() | pid->onlyKaon() );    // seperater Pion/Kaon/Electron
		pid->calculate();
		prob[i][0] = pid->probElectron();
		prob[i][1] = pid->probPion();
		prob[i][2] = pid->probKaon();
		if(!(pid->IsPidInfoValid())) continue;

		if(_charge==1){
			if((prob[i][0]>prob[i][1])&&(prob[i][0]>prob[i][2])) {flag_ep=i;N_eppid++;}
			else flag_pip=i;
		}
		else if(_charge==-1){
			if((prob[i][0]>prob[i][1])&&(prob[i][0]>prob[i][2])) {flag_em=i;N_empid++;}
			else flag_pim=i;
		}
		else continue;
	}
	if((N_eppid!=1)||(N_empid!=1))  return StatusCode::SUCCESS;
	Ncut3++;

	//////////////////// Assign track  info
	for(int i = 0; i < nGood; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();

		HepVector ia = mdcTrk->helix();
		HepSymMatrix iEa = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		VFHelix ihelixip(point0,ia,iEa); 
		ihelixip.pivot(IP);
		HepVector ivecipa = ihelixip.a();
		double  Rvxy0=fabs(ivecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=ivecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=ivecipa[1];

		costheta_chrgd[i]=cos(mdcTrk->theta());
		Rxy[i]=mdcTrk->r();
		Rz[i]=mdcTrk->z();
		Rvxy[i]=Rvxy0;
		Rvz[i]=Rvz0;	
		pKal[i]=mdcKalTrk->p();
		pxKal[i]=mdcKalTrk->px();	
		pyKal[i]=mdcKalTrk->py();	
		pzKal[i]=mdcKalTrk->pz();	
		deemc[i]=-1.;
		eop[i]=-1.;
		if( (*itTrk)->isEmcShowerValid() ){
			RecEmcShower *emcTrk_e = (*itTrk)->emcShower();
			deemc[i] = emcTrk_e->energy();
			eop[i]=deemc[i]/pKal[i];
		}
	}

	WTrackParameter wvpipTrk, wvpimTrk,wvepTrk, wvemTrk;
	// for electron	
	RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
	RecMdcKalTrack *mdcTrk_ep = (*(evtRecTrkCol->begin()+iGood[flag_ep]))->mdcKalTrack();
	RecMdcKalTrack *mdcTrk_em = (*(evtRecTrkCol->begin()+iGood[flag_em]))->mdcKalTrack();
	HepLorentzVector ep_un = mdcTrk_ep->p4(me);
	HepLorentzVector em_un = mdcTrk_em->p4(me);
	p_uep->SetPxPyPzE(ep_un.px(), ep_un.py(), ep_un.pz(), ep_un.e());
	p_uem->SetPxPyPzE(em_un.px(), em_un.py(), em_un.pz(), em_un.e());

	if(runNo>0 || (runNo<0&&(mc_cor==0))){
		wvepTrk = WTrackParameter(me, mdcTrk_ep->getZHelixE(), mdcTrk_ep->getZErrorE());  
		wvemTrk = WTrackParameter(me, mdcTrk_em->getZHelixE(), mdcTrk_em->getZErrorE());
	}	
	else if(runNo<0&&mc_cor)
	{
		HepVector wep_zHel(5,0);HepVector wem_zHel(5,0);
		EtaPiPiEE::calibration(mdcTrk_ep,wep_zHel,1);//1:e
		EtaPiPiEE::calibration(mdcTrk_em,wem_zHel,1);//1:e
		wvepTrk = WTrackParameter(me, wep_zHel, mdcTrk_ep->getZErrorE());
		wvemTrk = WTrackParameter(me, wem_zHel, mdcTrk_em->getZErrorE());
	}

	//for pion
	RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
	RecMdcKalTrack *mdcTrk_pip = (*(evtRecTrkCol->begin()+iGood[flag_pip]))->mdcKalTrack();
	RecMdcKalTrack *mdcTrk_pim = (*(evtRecTrkCol->begin()+iGood[flag_pim]))->mdcKalTrack();
	HepLorentzVector pip_un = mdcTrk_pip->p4(mpi);
	HepLorentzVector pim_un = mdcTrk_pim->p4(mpi);
	p_upip->SetPxPyPzE(pip_un.px(), pip_un.py(), pip_un.pz(), pip_un.e());
	p_upim->SetPxPyPzE(pim_un.px(), pim_un.py(), pim_un.pz(), pim_un.e());

	if(runNo>0 || (runNo<0&&(mc_cor==0))){
		wvpipTrk = WTrackParameter(mpi, mdcTrk_pip->getZHelix(), mdcTrk_pip->getZError());  
		wvpimTrk = WTrackParameter(mpi, mdcTrk_pim->getZHelix(), mdcTrk_pim->getZError());  
	}
	else if(runNo<0&&mc_cor)
	{
		HepVector wpip_zHel(5,0);HepVector wpim_zHel(5,0);
		EtaPiPiEE::calibration(mdcTrk_pip,wpip_zHel,0);//0:pi
		EtaPiPiEE::calibration(mdcTrk_pim,wpim_zHel,0);//0:pi
		wvpipTrk = WTrackParameter(mpi, wpip_zHel, mdcTrk_pip->getZError());
		wvpimTrk = WTrackParameter(mpi, wpim_zHel, mdcTrk_pim->getZError());
	}

	/////////////////////////////////////vertex fit
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
	vtxfit->AddTrack(0,  wvepTrk);
	vtxfit->AddTrack(1,  wvemTrk);
	vtxfit->AddTrack(2,  wvpipTrk);
	vtxfit->AddTrack(3,  wvpimTrk);
	vtxfit->AddVertex(0, vxpar,0,1,2,3);

	if(!vtxfit->Fit(0)) return StatusCode::SUCCESS;
	vtxfit->Swim(0);
	vtxchisq=vtxfit->chisq(0);

	WTrackParameter wepTrk = vtxfit->wtrk(0);
	WTrackParameter wemTrk = vtxfit->wtrk(1);
	WTrackParameter wpipTrk = vtxfit->wtrk(2);
	WTrackParameter wpimTrk = vtxfit->wtrk(3);

	Ncut4++;   // finish vertex fit

	/////////////////////////////////////4C fit
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	TLorentzVector p4_jpsi(p4_cms.px(),p4_cms.py(),p4_cms.pz(),p4_cms.e());	

	chi4C = 99999.;

	for(int i=0; i < nGamma ; i++){
		RecEmcShower* gammaTrk = (*(evtRecTrkCol->begin()+iGamma[i]))->emcShower();
		kmfit->init();
		kmfit->setEspread(m_EnergySpread);
		kmfit->setBeamPosition(xorigin);
		kmfit->setVBeamPosition(Evx_avg);
		kmfit->AddTrack(0,wepTrk);
		kmfit->AddTrack(1,wemTrk);
		kmfit->AddTrack(2,wpipTrk);
		kmfit->AddTrack(3,wpimTrk);
		kmfit->AddTrack(4,0.0,gammaTrk);
		kmfit->AddFourMomentum(0,p4_cms);

		bool okfit = kmfit->Fit();
		if(!okfit) continue;
		if(kmfit->chisq() > chi4C) continue;
		chi4C=kmfit->chisq();

		ep_pull_0 = kmfit->pull(0)[0];
		ep_pull_1 = kmfit->pull(0)[1];
		ep_pull_2 = kmfit->pull(0)[2];
		ep_pull_3 = kmfit->pull(0)[3];
		ep_pull_4 = kmfit->pull(0)[4];
		em_pull_0 = kmfit->pull(1)[0];
		em_pull_1 = kmfit->pull(1)[1];
		em_pull_2 = kmfit->pull(1)[2];
		em_pull_3 = kmfit->pull(1)[3];
		em_pull_4 = kmfit->pull(1)[4];
		pip_pull_0 = kmfit->pull(2)[0];
		pip_pull_1 = kmfit->pull(2)[1];
		pip_pull_2 = kmfit->pull(2)[2];
		pip_pull_3 = kmfit->pull(2)[3];
		pip_pull_4 = kmfit->pull(2)[4];
		pim_pull_0 = kmfit->pull(3)[0];
		pim_pull_1 = kmfit->pull(3)[1];
		pim_pull_2 = kmfit->pull(3)[2];
		pim_pull_3 = kmfit->pull(3)[3];
		pim_pull_4 = kmfit->pull(3)[4];

		HepLorentzVector ep_ = kmfit->pfit(0);
		HepLorentzVector em_ = kmfit->pfit(1);
		HepLorentzVector pip_ = kmfit->pfit(2);
		HepLorentzVector pim_ = kmfit->pfit(3);
		HepLorentzVector gamma_ = kmfit->pfit(4);

		p_ep->SetPxPyPzE(ep_.px(), ep_.py(), ep_.pz(), ep_.e());
		p_em->SetPxPyPzE(em_.px(), em_.py(), em_.pz(), em_.e());
		p_pip->SetPxPyPzE(pip_.px(), pip_.py(), pip_.pz(), pip_.e());
		p_pim->SetPxPyPzE(pim_.px(), pim_.py(), pim_.pz(), pim_.e());
		p_gamma->SetPxPyPzE(gamma_.px(), gamma_.py(), gamma_.pz(), gamma_.e());

		p_ugamma->SetPxPyPzE((pGamma[i]).px(), (pGamma[i]).py(), (pGamma[i]).pz(), (pGamma[i]).e());
	}
	if( chi4C >= 99999. ) return StatusCode::SUCCESS;
	Ncut5++;	//after 4C

	*p_ee=*p_ep+*p_em;
	*p_pipi=*p_pip+*p_pim;
	*p_gee=*p_ep+*p_em+*p_gamma;
	*p_gpipi=*p_pip+*p_pim+*p_gamma;
	*p_uee=*p_uep+*p_uem;
	*p_upipi=*p_upip+*p_upim;
	*p_ugee=*p_uep+*p_uem+*p_ugamma;
	*p_ugpipi=*p_upip+*p_upim+*p_ugamma;
	*p_recpipi=p4_jpsi-*p_pipi;

	m_gpipi=p_gpipi->M();
	m_gee=p_gee->M();
	m_ee=p_ee->M();
	m_uee=p_uee->M();
	m_recpipi=p_recpipi->M();
	int flag_highe=pKal[flag_ep]>pKal[flag_em]?flag_ep:flag_em;
	int flag_lowe=pKal[flag_ep]<pKal[flag_em]?flag_ep:flag_em;
	highe_eop=eop[flag_highe];
	lowe_eop=eop[flag_lowe];


	//Gamma Conv (chuxk)
	HepPoint3D OP(0,0,0);
	{
		RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
		GammaConv gconv = GammaConv(mdcTrk_ep->helix(), mdcTrk_em->helix(),OP);
		m_xconv1 = gconv.getRX1();
		m_yconv1 = gconv.getRY1();
		m_zconv1 = gconv.getRZ1();
		m_rconv1 = gconv.getRXY1();

		m_xconv2 = gconv.getRX2();
		m_yconv2 = gconv.getRY2();
		m_zconv2 = gconv.getRZ2();
		m_rconv2 = gconv.getRXY2();

		m_xiep  = gconv.getXiep();
		m_deltaxy=gconv.getDeltaXY();

		m_deltaz1 = gconv.getDeltaZ1();
		m_deltaz2 = gconv.getDeltaZ2();

		m_psipair = gconv.getPsipair();
		//            m_dgamma  = gconv.getDGamma();
		m_rp      = gconv.getRp();
		m_re      = gconv.getEp();
		m_deltaeq = gconv.getDeltaeq();
		m_lep     = gconv.getLep();
		m_case   = gconv.getNcase();

		TVector3 v_conv1(m_xconv1,m_yconv1,m_zconv1);
		TVector3 v_conv2(m_xconv2,m_yconv2,m_zconv2);
		TVector3 v_gamma;
		HepLorentzVector gammaVector(0,0,0,0);
		gammaVector=HepLorentzVector(mdcTrk_ep->p4(me)+mdcTrk_em->p4(me));
		v_gamma.SetXYZ(gammaVector.px(),gammaVector.py(),gammaVector.pz());

		m_thetaeg1  = v_gamma.Angle(-v_conv1)*180/(CLHEP::pi);
		m_thetaeg2  = v_gamma.Angle(-v_conv2)*180/(CLHEP::pi);
		m_cthep = (mdcTrk_ep->p4(me)).vect()*(mdcTrk_em->p4(me)).vect()/mdcTrk_ep->p()/mdcTrk_em->p();
		m_ptrkp  = mdcTrk_ep->p();
		m_ptrkm  = mdcTrk_em->p();
		m_mgamma = gammaVector.m();
		m_egamma = gammaVector.e();
		m_theta        = gammaVector.theta();
		m_cosTheta     = gammaVector.cosTheta();
		m_phi          = gammaVector.phi();  
	}//Gamma Conv        
	angee = p_uep->Angle(p_uem->Vect()); //rad
	angee = angee * 180 / (CLHEP::pi); //deg
	m_rconv=angee<10?m_rconv2:m_rconv1;	

	//Gamma Conv (zhangyt)
	{
		HepPoint3D vxepem(0., 0., 0.);
		HepSymMatrix Evxepem(3, 0);
		double bxepem = 1E+6;
		double byepem = 1E+6;
		double bzepem = 1E+6;
		Evxepem[0][0] = bxepem*bxepem;
		Evxepem[1][1] = byepem*byepem;
		Evxepem[2][2] = bzepem*bzepem;

		VertexParameter vxparepem;
		vxparepem.setVx(vxepem);
		vxparepem.setEvx(Evxepem);
		VertexFit* vtxfite = VertexFit::instance();
		vtxfite->init();
		vtxfite->AddTrack(0,  wvepTrk);
		vtxfite->AddTrack(1,  wvemTrk);
		vtxfite->AddVertex(0, vxparepem,0, 1);
		if(!vtxfite->Fit(0)) vtxchie=-1.;;
		vtxfite->Swim(0);
		vtxchie=vtxfit->chisq(0);
		HepPoint3D xorigin1e = vtxfite->vx(0);
		HepSymMatrix xem1e = vtxfite->Evx(0);

		RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
	GammaConv gconv = GammaConv(mdcTrk_em->helix(),mdcTrk_ep->helix(),OP);
		m_epemx = gconv.getRX();
		m_epemy = gconv.getRY();
		m_epemz= gconv.getRZ();
		m_epemxy = gconv.getRXY();

		HepPoint3D  eptrack=gconv.center(mdcTrk_ep->helix()); 
		double  eprtrack=gconv.radius(mdcTrk_ep->helix());

		m_mepr=fabs(eprtrack);
		m_mepcenterx=eptrack.x();
		m_mepcentery=eptrack.y();
		m_mepcenterz=eptrack.z();

		HepPoint3D  emtrack=gconv.center(mdcTrk_em->helix()); 
		double  emrtrack=gconv.radius(mdcTrk_em->helix());

		m_memr=fabs(emrtrack);
		m_memcenterx=emtrack.x();
		m_memcentery=emtrack.y();
		m_memcenterz=emtrack.z();

		m_epemxxorigin1e=xorigin1e.x();
		m_epemyxorigin1e=xorigin1e.y();
		m_epemzxorigin1e=xorigin1e.z();
		m_epemxyxorigin1e=sqrt(xorigin1e.x()*xorigin1e.x()+xorigin1e.y()*xorigin1e.y());

		m_rconvz=angee<10?m_epemxy:m_epemxyxorigin1e;	
	}	


	TreeAna->Fill();
	GammaAll->Clear();

	return StatusCode::SUCCESS;
}

//***************************************************************************

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode EtaPiPiEE::finalize() {
	NbInfo->Fill();
	saveFile->cd();
	TreeAna->Write();
	if(m_saveCutFlow == 1) NbInfo->Write();
	saveFile->Close();
	cout<<endl<<"Finalize psi'-> eta' e+ e-  ,eta'->gam pi+ pi-"<<endl;;
	cout<<"Total number:                         "<<Ncut0<<endl;
	cout<<"nGood==4, nCharge==0:                 "<<Ncut1<<endl;
	cout<<"nGamma>=1:                            "<<Ncut2<<endl;
	cout<<"PID success:                          "<<Ncut3<<endl;
	cout<<"vertex fit:                           "<<Ncut4<<endl;
	cout<<"4C success:                           "<<Ncut5<<endl;
	cout<<"_END_TAG_"<<endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}


// corset(): sets up the generation by calculating C from V.
void EtaPiPiEE::corset(HepSymMatrix &V, HepMatrix &C, int n)
{
	//cout<<"v="<<V<<endl;
	//cout<<"c="<<C<<endl;
	double sum;

	// Compute square root of matrix sigma
	for (int j=0; j<n; j++) {
		sum = 0;
		for (int k=0; k<j; k++) {
			sum = sum + C[j][k]*C[j][k];

		}
		//cout<<"sum="<<sum<<endl;
		//cout<<"v("<<j<<","<<j<<")="<<V[j][j]<<endl;
		C[j][j] = sqrt(abs(V[j][j] - sum));
		//cout<<"c("<<j<<","<<j<<")="<<C[j][j]<<endl;
		// Off Diagonal terms
		for (int i=j+1; i<n; i++) {
			sum = 0;
			for (int k=0; k<j; k++) {
				sum = sum + C[i][k]*C[j][k];
			}
			C[i][j] = (V[i][j] - sum)/C[j][j];
			//cout<<"C("<<i<<","<<j<<")="<<C[i][j]<<endl;
		}
	}
}
// corgen(): generates a set of n random numbers Gaussian-distributed with covariance
// matrix V (V = C*C') and mean values zero.
void EtaPiPiEE::corgen(HepMatrix &C, HepVector &x, int n)
{
	int i,j;
	int nmax = 100;

	if (n > nmax ) {
		printf("Error in corgen: array overflown");
	}

	double tmp[3];
	for(int p = 0 ; p < n; p ++){
		tmp[p] = gRandom->Gaus(0,1);
		//cout<<"tmp["<<p<<"]="<<tmp[p]<<endl;
	}
	for ( i=0; i<n; i++) {
		x[i] = 0.0;
		for (j=0; j<=i; j++) {
			x[i] = x[i]+C[i][j]*tmp[j];
		}
	}
}
//**********************************************
void EtaPiPiEE::calibration(RecMdcKalTrack *trk , HepVector &wtrk_zHel, int n )
{

	HepVector pip_calerr_d2(5,0);
	HepVector pim_calerr_d2(5,0);
	HepVector ep_calerr_d2(5,0);
	HepVector em_calerr_d2(5,0);


	pip_calerr_d2[0] = 1.0;
	pip_calerr_d2[1] = 1.07663;//1.16;
	pip_calerr_d2[2] = 1.11371;//1.17;
	pip_calerr_d2[3] = 1.0;
	pip_calerr_d2[4] = 1.01258;//1.14;

	pim_calerr_d2[0] = 1.0;
	pim_calerr_d2[1] = 1.1001;//1.16;
	pim_calerr_d2[2] = 1.08598;//1.14;
	pim_calerr_d2[3] = 1.0;
	pim_calerr_d2[4] = 1.00512; //1.14;

	ep_calerr_d2[0] = 1.0;
	ep_calerr_d2[1] = 1.08702;
	ep_calerr_d2[2] = 1.02458;//1.115107;//*1.01927;
	ep_calerr_d2[3] = 1.0;
	ep_calerr_d2[4] = 1.02149;//1.133133*1.01432;

	em_calerr_d2[0] = 1.0;
	em_calerr_d2[1] = 1.06805;//1.19767*1.03072;
	em_calerr_d2[2] = 1.014;//1.099758;//*1.00574;
	em_calerr_d2[3] = 1.0;
	em_calerr_d2[4] = 1.02677;//1.13305*1.011161;

	HepVector pip_calmean_d2(5,0);
	HepVector pim_calmean_d2(5,0);
	HepVector ep_calmean_d2(5,0);
	HepVector em_calmean_d2(5,0);


	pip_calmean_d2[0] = 0;
	pip_calmean_d2[1] = -0.0221936;//-0.01;
	pip_calmean_d2[2] = 0.134811;//0.43;
	pip_calmean_d2[3] = 0;
	pip_calmean_d2[4] = 0.0906794;//0.12;

	pim_calmean_d2[0] = 0;
	pim_calmean_d2[1] = 0.103221;//0.04;
	pim_calmean_d2[2] = -0.106784;//-0.38;
	pim_calmean_d2[3] = 0;
	pim_calmean_d2[4] = 0.0689878;//0.17;

	ep_calmean_d2[0] = 0;
	ep_calmean_d2[1] = -0.106422;//-0.002068/2 - 0.0/2;
	ep_calmean_d2[2] = -0.0335828;//0.271361/2 ;//- 0.088196/2;
	ep_calmean_d2[3] = 0;
	ep_calmean_d2[4] = 0.0792581;//1.482376/2 - 0.49587/2;

	em_calmean_d2[0] = 0;
	em_calmean_d2[1] = -0.0943685; //0.07331/2 - 0.010069/2;
	em_calmean_d2[2] = 0.0161445;//-0.126777/2 ;//+ 0.0275825/2;
	em_calmean_d2[3] = 0;
	em_calmean_d2[4] = 0.0814989;//1.492618/2 - 0.49516/2;


	if(trk->charge()>0 && n==0){
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		//pip
		HepSymMatrix wpip_zerr(5,0);
		wpip_zerr = trk->getZError();
		//cout<<"wpip_zerr="<<wpip_zerr<<endl;
		HepSymMatrix wpip_zcal(3,0);
		wpip_zcal[0][0] = (pip_calerr_d2[1]*pip_calerr_d2[1]-1)*wpip_zerr[1][1];
		wpip_zcal[1][1] = (pip_calerr_d2[2]*pip_calerr_d2[2]-1)*wpip_zerr[2][2];
		wpip_zcal[2][2] = (pip_calerr_d2[4]*pip_calerr_d2[4]-1)*wpip_zerr[4][4];
		//cout<<"wpip_zcal="<<wpip_zcal<<endl;

		HepMatrix wpip_zerrc(3,3,0);
		EtaPiPiEE::corset(wpip_zcal,wpip_zerrc,3);
		HepVector wpip_zgen(3,0);
		EtaPiPiEE::corgen(wpip_zerrc,wpip_zgen,3);

		wtrk_zHel[0] = trk->getZHelix()[0];
		wtrk_zHel[1] = trk->getZHelix()[1]+pip_calmean_d2[1]*sqrt(wpip_zerr[1][1])+wpip_zgen[0];
		wtrk_zHel[2] = trk->getZHelix()[2]+pip_calmean_d2[2]*sqrt(wpip_zerr[2][2])+wpip_zgen[1];
		wtrk_zHel[3] = trk->getZHelix()[3];
		wtrk_zHel[4] = trk->getZHelix()[4]+pip_calmean_d2[4]*sqrt(wpip_zerr[4][4])+wpip_zgen[2];

		//cout<<"wtrk_zHel="<<wtrk_zHel<<endl;
	}

	if(trk->charge()<0 && n==0)
	{
		RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
		//pim
		HepSymMatrix wpim_zerr(5,0);
		wpim_zerr = trk->getZError();

		HepSymMatrix wpim_zcal(3,0);

		wpim_zcal[0][0] = (pim_calerr_d2[1]*pim_calerr_d2[1]-1)*wpim_zerr[1][1];
		wpim_zcal[1][1] = (pim_calerr_d2[2]*pim_calerr_d2[2]-1)*wpim_zerr[2][2];
		wpim_zcal[2][2] = (pim_calerr_d2[4]*pim_calerr_d2[4]-1)*wpim_zerr[4][4];

		HepMatrix wpim_zerrc(3,3,0);
		EtaPiPiEE::corset(wpim_zcal,wpim_zerrc,3);
		HepVector wpim_zgen(3,0);
		EtaPiPiEE::corgen(wpim_zerrc,wpim_zgen,3);
		wtrk_zHel[0] = trk->getZHelix()[0];
		wtrk_zHel[1] = trk->getZHelix()[1]+pim_calmean_d2[1]*sqrt(wpim_zerr[1][1])+wpim_zgen[0];
		wtrk_zHel[2] = trk->getZHelix()[2]+pim_calmean_d2[2]*sqrt(wpim_zerr[2][2])+wpim_zgen[1];
		wtrk_zHel[3] = trk->getZHelix()[3];
		wtrk_zHel[4] = trk->getZHelix()[4]+pim_calmean_d2[4]*sqrt(wpim_zerr[4][4])+wpim_zgen[2];

	}

	if(trk->charge()>0 && n==1)
	{
		//ep
		RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
		HepSymMatrix wep_zerr(5,0);
		wep_zerr = trk->getZErrorE();

		HepSymMatrix wep_zcal(3,0);
		wep_zcal[0][0] = (ep_calerr_d2[1]*ep_calerr_d2[1]-1)*wep_zerr[1][1];
		wep_zcal[1][1] = (ep_calerr_d2[2]*ep_calerr_d2[2]-1)*wep_zerr[2][2];
		wep_zcal[2][2] = (ep_calerr_d2[4]*ep_calerr_d2[4]-1)*wep_zerr[4][4];
		HepMatrix wep_zerrc(3,3,0);
		EtaPiPiEE::corset(wep_zcal,wep_zerrc,3);
		HepVector wep_zgen(3,0);
		EtaPiPiEE::corgen(wep_zerrc,wep_zgen,3);

		wtrk_zHel[0] = trk->getZHelixE()[0];
		wtrk_zHel[1] = trk->getZHelixE()[1]+ep_calmean_d2[1]*sqrt(wep_zerr[1][1])+wep_zgen[0];
		wtrk_zHel[2] = trk->getZHelixE()[2]+ep_calmean_d2[2]*sqrt(wep_zerr[2][2])+wep_zgen[1];
		wtrk_zHel[3] = trk->getZHelixE()[3];
		wtrk_zHel[4] = trk->getZHelixE()[4]+ep_calmean_d2[4]*sqrt(wep_zerr[4][4])+wep_zgen[2];

	}

	if(trk->charge()<0 && n==1)
	{
		RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
		//em
		HepSymMatrix wem_zerr(5,0);
		wem_zerr = trk->getZErrorE();

		HepSymMatrix wem_zcal(3,0);

		wem_zcal[0][0] = (em_calerr_d2[1]*em_calerr_d2[1]-1)*wem_zerr[1][1];
		wem_zcal[1][1] = (em_calerr_d2[2]*em_calerr_d2[2]-1)*wem_zerr[2][2];
		wem_zcal[2][2] = (em_calerr_d2[4]*em_calerr_d2[4]-1)*wem_zerr[4][4];
		HepMatrix wem_zerrc(3,3,0);
		EtaPiPiEE::corset(wem_zcal,wem_zerrc,3);
		HepVector wem_zgen(3,0);
		EtaPiPiEE::corgen(wem_zerrc,wem_zgen,3);
		wtrk_zHel[0] = trk->getZHelixE()[0];
		wtrk_zHel[1] = trk->getZHelixE()[1]+em_calmean_d2[1]*sqrt(wem_zerr[1][1])+wem_zgen[0];
		wtrk_zHel[2] = trk->getZHelixE()[2]+em_calmean_d2[2]*sqrt(wem_zerr[2][2])+wem_zgen[1];
		wtrk_zHel[3] = trk->getZHelixE()[3];
		wtrk_zHel[4] = trk->getZHelixE()[4]+em_calmean_d2[4]*sqrt(wem_zerr[4][4])+wem_zgen[2];

	}
}
