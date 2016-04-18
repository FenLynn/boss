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

#include "XisP3PiAlg/XisP3Pi.h"

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
#include "TRandom3.h"

#include "GammaConv/GammaConv.h"

//const double twopi = 6.2831853;
//const double pi = 3.1415927;
const double me = 0.000511;
const double mmu = 0.105658;
const double mpi = 0.13957;
const double mk = 0.493677;
const double mp = 0.938272;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double mel  = 0.000510998910;
const double Ejpsi=3.097;
const double Epsip=3.686;
const double Econt=3.08;
//const double Econt=3.773;

const double mLambda=1.115683;
const double mXi=1.32171;
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

XisP3Pi::XisP3Pi(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) {
		//Declare the properties
		declareProperty("EnergySpread",  m_EnergySpread = 0.0013);
		declareProperty("OutputFileName",  m_OutputFileName = "XisP3Pi_test.root");
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
		declareProperty("mc_cor",mc_cor=0);
		declareProperty("m_NGamma",m_NGamma=2);
	}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode XisP3Pi::initialize(){
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

	xep=new TLorentzVector();
	xem=new TLorentzVector();
	xgamma1=new TLorentzVector();
	xgamma2=new TLorentzVector();
	xgg=new TLorentzVector();
	xee=new TLorentzVector();
	xeta=new TLorentzVector();
	
	p_ep=new TLorentzVector();
	p_em=new TLorentzVector();
	p_gamma1=new TLorentzVector();
	p_gamma2=new TLorentzVector();
	p_ee=new TLorentzVector();
	p_gg=new TLorentzVector();
	
	p_gamma1b=new TLorentzVector();
	p_gamma2b=new TLorentzVector();

	p_uep=new TLorentzVector();
	p_uem=new TLorentzVector();
	p_ugamma1=new TLorentzVector();
	p_ugamma2=new TLorentzVector();
	p_uee=new TLorentzVector();
	p_ugg=new TLorentzVector();

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
		TreeAna->Branch("xep",&xep,32000,0);
		TreeAna->Branch("xem",&xem,32000,0);
		TreeAna->Branch("xgamma1",&xgamma1,32000,0);
		TreeAna->Branch("xgamma2",&xgamma2,32000,0);
		TreeAna->Branch("xgg",&xgg,32000,0);
		TreeAna->Branch("xee",&xee,32000,0);
		TreeAna->Branch("xeta",&xeta,32000,0);
		TreeAna->Branch("flagP", &flagP, "flagP/I");
		TreeAna->Branch("issig", &issig, "issig/I");
		TreeAna->Branch("xmee", &xmee, "xmee/D");
		TreeAna->Branch("xmgg", &xmgg, "xmgg/D");
	}
	TreeAna->Branch("runid", &runid, "runid/I");
	TreeAna->Branch("evtid", &evtid, "evtid/I");
	TreeAna->Branch("nNEUTRAL", &nNEUTRAL, "nNEUTRAL/I");
	TreeAna->Branch("nCHARGED",&nCHARGED,"nCHARGED/I");
	TreeAna->Branch("nTRACKS",&nTRACKS,"nTRACKS/I");
//	TreeAna->Branch("nGamma", &nGamma, "nGamma/I");

	TreeAna->Branch("GammaAll","TClonesArray",&GammaAll,256000,0);
	TreeAna->Branch("costheta_gamma", costheta_gamma, "costheta_gamma[nGamma]/D");
	TreeAna->Branch("energy_gamma", energy_gamma, "energy_gamma[nGamma]/D");
	TreeAna->Branch("TDCtime", TDCtime , "TDCtime[nGamma]/D");
	TreeAna->Branch("isoAngle", isoAngle, "isoAngle[nGamma]/D");
	TreeAna->Branch("showerde", &showerde, "showerde/D");

	TreeAna->Branch("vx", vx, "vx[3]/D");
	TreeAna->Branch("Evx", Evx, "Evx[3]/D");
	TreeAna->Branch("costheta_chrgd", costheta_chrgd, "costheta_chrgd[2]/D");
	TreeAna->Branch("Rxy", Rxy, "Rxy[2]/D");
	TreeAna->Branch("Rz", Rz, "Rz[2]/D");
	TreeAna->Branch("Rvxy", Rvxy, "Rvxy[2]/D");
	TreeAna->Branch("Rvz", Rvz, "Rvz[2]/D");
	TreeAna->Branch("pKal", pKal, "pKal[2]/D");
	TreeAna->Branch("pxKal", pxKal, "pxKal[2]/D");
	TreeAna->Branch("pyKal", pyKal, "pyKal[2]/D");
	TreeAna->Branch("pzKal", pzKal, "pzKal[2]/D");
	TreeAna->Branch("deemc", deemc, "deemc[2]/D");
	TreeAna->Branch("eop", eop, "eop[2]/D");

	TreeAna->Branch("flag_ep", &flag_ep, "flag_ep/I");
	TreeAna->Branch("flag_em", &flag_em, "flag_em/I");

	TreeAna->Branch("prob", prob, "prob[2][3]/D");

	TreeAna->Branch("vtxchisq", &vtxchisq, "vtxchisq/D");
	TreeAna->Branch("chi4C", &chi4C, "chi4C/D");

	TreeAna->Branch("p_ep",&p_ep,32000,0);
	TreeAna->Branch("p_em",&p_em,32000,0);
	TreeAna->Branch("p_gamma1",&p_gamma1,32000,0);
	TreeAna->Branch("p_gamma2",&p_gamma2,32000,0);
	TreeAna->Branch("p_ee",&p_ee,32000,0);
	TreeAna->Branch("p_gg",&p_gg,32000,0);
	TreeAna->Branch("p_gamma1b",&p_gamma1b,32000,0);
	TreeAna->Branch("p_gamma2b",&p_gamma2b,32000,0);

	TreeAna->Branch("p_uep",&p_uep,32000,0);
	TreeAna->Branch("p_uem",&p_uem,32000,0);
	TreeAna->Branch("p_ugamma1",&p_ugamma1,32000,0);
	TreeAna->Branch("p_ugamma2",&p_ugamma2,32000,0);
	TreeAna->Branch("p_uee",&p_uee,32000,0);
	TreeAna->Branch("p_ugg",&p_ugg,32000,0);
	
	TreeAna->Branch("angee", &angee, "angee/D");
	TreeAna->Branch("mee", &mee, "mee/D");
	TreeAna->Branch("mgg", &mgg, "mgg/D");
	TreeAna->Branch("cosgam1b", &cosgam1b, "cosgam1b/D");
	TreeAna->Branch("cosgam2b", &cosgam2b, "cosgam2b/D");
	TreeAna->Branch("cosgamb", &cosgamb, "cosgamb/D");

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

	//Gamma conversion
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
StatusCode XisP3Pi::execute() {
	//std::cout << "execute()" << std::endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	int runNo=eventHeader->runNumber();
	int event=eventHeader->eventNumber();
	runid=runNo;
	evtid=event;

	Ncut0++;  //total events
	log << MSG::DEBUG <<"run, evtnum = "<< runNo << " , "<< event <<endreq;
//	cout <<"run, evtnum = "<< runNo << " , "<< event <<endl;

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
	//*************************************************************MC Info.***************************************************************************

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
		// MC Truth
		/*
		if(m_saveMCTruth == 1)
		{
			xep->SetE(-1);
			xem->SetE(-1);
			xgamma1->SetE(-1);
			xgamma2->SetE(-1);
			xee->SetE(-1);
			xgg->SetE(-1);

			int numep = 0;
			int numem = 0;
			int numgam = 0;
			int numeta = 0;
			int etaIndex = -1;
			int numDaughters = 0;

			if(jpsiIndex > 0)
			{
				for(iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
				{
					if((*iter_mc)->primaryParticle())	continue;
					if(!(*iter_mc)->decayFromGenerator()) continue;
					if(((*iter_mc)->mother()).trackIndex() !=jpsiIndex)	continue;
					int pid = (*iter_mc)->particleProperty();
					if(pid == -22) 	continue;
					HepLorentzVector p = (*iter_mc)->initialFourMomentum();
					if((pid == 221)||(pid==111)){
						xeta->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());
						etaIndex = (*iter_mc)->trackIndex();
						if(pid==221) flagP=1; //eta
						else if(pid==111) flagP=2; //pi0
						else continue;
						numeta++;
						numDaughters++;
					}
					else if( pid == -11 )					{
						xep->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());
						numep++;
						numDaughters++;
					}
					else if( pid == 11 )					{
						xem->SetPxPyPzE(p.px(),p.py(),p.pz(),p.e());
						numem++;
						numDaughters++;
					}
					else continue;
				}
			}

			if(etaIndex > 0 )
			{
				for(iter_mc = mcParticleCol->begin(); iter_mc != mcParticleCol->end(); iter_mc++)
				{
					if((*iter_mc)->primaryParticle())continue;
					if(!(*iter_mc)->decayFromGenerator()) continue;
					int pid = (*iter_mc)->particleProperty();
					if(pid == -22) continue;
					HepLorentzVector p = (*iter_mc)->initialFourMomentum();
					if(((*iter_mc)->mother()).trackIndex() == etaIndex){
						if(pid == 22){
							if(xgamma1->E()==-1) {
								xgamma1->SetPxPyPzE(p.px(), p.py(), p.pz(), p.e());
								numgam++;
								numDaughters++;
							}
							else {
								xgamma2->SetPxPyPzE(p.px(), p.py(), p.pz(), p.e());
								numgam++;
								numDaughters++;
							}
						}
					}
				}
			}
			xee->SetPxPyPzE((*xep+*xem).Px(),(*xep+*xem).Py(),(*xep+*xem).Pz(),(*xep+*xem).E());
			xgg->SetPxPyPzE((*xgamma1+*xgamma2).Px(),(*xgamma1+*xgamma2).Py(),(*xgamma1+*xgamma2).Pz(),(*xgamma1+*xgamma2).E());
			xmee=xee->M();
			xmgg=xgg->M();
			if((xep->E()>0)&&(xem->E()>0)&&(xgamma1->E()>0)&&(xgamma2->E()>0)&&(xeta->E()>0)&&(numeta==1) &&(numep==1)&&(numem==1)&&(numgam==2)) issig=1;
		}
	*/
	}

	//*************************************************************Global Event Parameters************************************************************
	Vint iGood, iGamma;
    Vint iTrp,iTrm,iprop,iprom,ipip,ipim,iem;
	iGood.clear();
	int nGood = 0;
	int nCharge = 0;
	int nGamma = 0;
	
	n_pip=0;
	n_em=0;
	n_pim=0;
        n_prop=0;
        n_prom=0;
        int nTrp=0;
        int nTrm=0;
	
	
	iGood.clear();
	iGamma.clear();
        iTrp.clear();
        iTrm.clear();
	ipip.clear();
	ipim.clear();
        iprop.clear();
        iprom.clear();

	//*************************************************************Primary Vertex*********************************************************************

	Hep3Vector xorigin(0,0,0);
//	HepSymMatrix Evx_avg(3,0);
        HepSymMatrix VtxErr(3,0);

	IVertexDbSvc*  vtxsvc;
	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if(vtxsvc->isVertexValid()){
	{
		double* dbv = vtxsvc->PrimaryVertex(); 
		double*  vv = vtxsvc->SigmaPrimaryVertex();  
		xorigin.setX(dbv[0]);
		xorigin.setY(dbv[1]);
		xorigin.setZ(dbv[2]);

		VtxErr[0][0] = vv[0]*vv[0];
		VtxErr[1][1] = vv[0]*vv[0];
		VtxErr[2][2] = vv[2]*vv[2];

	}
    else  return StatusCode::SUCCESS;
         VertexParameter bs;
	 bs.setVx(xorigin);
	 bs.setEvx(VtxErr); 
	
	 //*************************************************************Good charged track*****************************************************************
         HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
	for(int i = 0; i < evtRecEvent->totalCharged(); i++)
	{
		if(i>= evtRecTrkCol->size()) break;

		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		if(!(*itTrk)->isMdcKalTrackValid()) continue;

		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();

		double theta = mdcTrk->theta();

		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
//		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
		VFHelix helixip(point0,a,Ea); 
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0 = fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0 = vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0 = vecipa[1];
         HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
		if(fabs(cos(theta)) > m_trkAngCut) continue;
//		if(fabs(Rvxy0) >= m_vr0cut) continue;
//		if(fabs(Rvz0) >= m_vz0cut) continue;
		RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();

		if(mdcKalTrk->charge()==0) continue;
		if(mdcKalTrk->charge() > 0)  iTrp.push_back(i);
        if(mdcKalTrk->charge() < 0)  iTrm.push_back(i);
		
		nCharge += mdcKalTrk->charge();
		iGood.push_back(i);

	}
	nGood = iGood.size();
    nTrp=iTrp.size();
    nTrm=iTrm.size();
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;
//	if((nGood != 2) || (nCharge!=0) ) return StatusCode::SUCCESS;
	if(nGood >12) return StatusCode::SUCCESS;
       if(nTrp<1||nTrm<2)  return StatusCode::SUCCESS;
	   Ncut1++;		//after good track selection
	
//***************************Begin PID***********************************************************************************************
	flag_prop=-1;flag_pim[0]=-1;flag_pim[1]=-1
	int Npid_prop(0),Npid_pim(0);
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
		pid->identify(pid->onlyProton() | pid->onlyPion() | pid->onlyKaon() );    // seperater Pion/Kaon/Electron
		pid->calculate();
	//	prob[i][0] = pid->probElectron();
		prob_pi = pid->probPion();
		prob_K = pid->probKaon();
		prob_pro = pid->probProton();
		if(!(pid->IsPidInfoValid())) continue;
		if(_charge==1){
			if((prob_pro>prob_pi)&&(prob_pro>prob_K)) {flag_prop=i;Npid_prop++;}
		}
		else if(_charge==-1){
			if((prob_pi>prob_pro)&&(prob_pi>prob_K)) {
				if(Npid_pim==0) flag_pim[0]=i;		// for track, use ->beign()+iGood[i];
				else flag_pim[1]=i;
				Npid_pim++;
			}
		}
		else continue;
	}
	if((Npid_prop !=1)||(Npid_pim != 2))  return StatusCode::SUCCESS;
	Ncut2++;
	
	//*******************************Begin Assign track  info***********************************************************
	for(int i = 0; i < nGood; i++) {
		EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
		RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
		RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
		
		HepVector ia = mdcTrk->helix();
		HepSymMatrix iEa = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
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
	WTrackParameter wvpropTrk, wvpimTrk[2];
	RecMdcKalTrack *mdcTrk_pim[2];
	HepLorentzVector hvun_pim[2];
		
	RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
	RecMdcKalTrack *mdcTrk_prop = (*(evtRecTrkCol->begin()+iGood[flag_prop]))->mdcKalTrack();
	HepLorentzVector hvun_prop = mdcTrk_prop->p4(mp);
	pu_prop->SetPxPyPzE(prop_un.px(), prop_un.py(), prop_un.pz(), prop_un.e());
	wvpropTrk = WTrackParameter(mp, mdcTrk_prop->getZHelixP(), mdcTrk_ep->getZErrorP());  

	RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
	mdcTrk_pim[0] = (*(evtRecTrkCol->begin()+iGood[flag_pim[0]]))->mdcKalTrack();
	mdcTrk_pim[1] = (*(evtRecTrkCol->begin()+iGood[flag_pim[1]]))->mdcKalTrack();

	hvun_pim[0] = mdcTrk_pim[0]->p4(mpi);
	hvun_pim[1] = mdcTrk_pim[1]->p4(mpi);

	pu_pim1->SetPxPyPzE(hvun_pim[0].px(), hvun_pim[0].py(), hvun_pim[0].pz(), hvun_pim[0].e());
	pu_pim2->SetPxPyPzE(hvun_pim[1].px(), hvun_pim[1].py(), hvun_pim[1].pz(), hvun_pim[1].e());

	wvpimTrk[0] = WTrackParameter(mpi, mdcTrk_pim[0]->getZHelix(), mdcTrk_pim[0]->getZError());  
	wvpimTrk[1] = WTrackParameter(mpi, mdcTrk_pim[1]->getZHelix(), mdcTrk_pim[1]->getZError());  
//****************************************
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
//****************************************
flag_lampim=-1;
flag_Xipim=-1;
nVtxok=0;
chivtx1=-1.;

for(int i=0;i<2;i++){
		VertexFit* vtxfit1 = VertexFit::instance();
		HepLorentzVector VLambda(0,0,0,0);
		VLambda=hvun_pim[i]+pu_prop;
		if( fabs(VLambda.m()-mLambda)>0.020) continue;
			
        vtxfit1->init();
        vtxfit1->AddTrack(0,  wvpropTrk);
        vtxfit1->AddTrack(1,  wvpimTrk[i]);
        vtxfit1->AddVertex(0, vxpar, 0, 1);
        if( !(vtxfit->Fit(0)) ) continue;
        vtxfit1->Swim(0);
        vtxfit1->BuildVirtualParticle(0);
        wvLambda=vtxfit1->wVirtualTrack(0);
//		chivtx[0	
		int another_pim=abs(1-i);
		
		vtxfit1->init();
        
		VertexFit* vtxfit2 = VertexFit::instance();
        vtxfit2->init();
        vtxfit2->AddTrack(0, wvLambda);
        vtxfit2->AddTrack(1, wvpimTrk[another_pim]);
        vtxfit4->AddVertex(0, vxpar, 0, 1);
			






	   
	//***********************************************************Good Gamma Selection*****************************************************************
	showerde=0.;
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

	if(nGamma < m_NGamma || nGamma >10) return StatusCode::SUCCESS; 
	Ncut2++;	//after good gamma cut
	//***************************Finish Good gamma selection***********************************************************************************************
	//***************************Assign 4-momentum to each photon**************************************************************
	TLorentzVector p_JpsiGamma;
	Vp4 pGamma;
	pGamma.clear();
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



//	WTrackParameter wvepTrk, wvemTrk;

	RecMdcKalTrack::setPidType(RecMdcKalTrack::proton);
	RecMdcKalTrack *mdcTrk_pp = (*(evtRecTrkCol->begin()+iGood[flag_pp]))->mdcKalTrack();
	HepLorentzVector pp_un = mdcTrk_pp->p4(mp);
	p_upp->SetPxPyPzE(pp_un.px(), pp_un.py(), pp_un.pz(), pp_un.e());

	RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
	RecMdcKalTrack *mdcTrk_pim1 = (*(evtRecTrkCol->begin()+iGood[flag_pim[0]]))->mdcKalTrack();
	RecMdcKalTrack *mdcTrk_pim2 = (*(evtRecTrkCol->begin()+iGood[flag_pim[1]]))->mdcKalTrack();
	HepLorentzVector pim1_un = mdcTrk_pim1->p4(mpi);
	HepLorentzVector pim2_un = mdcTrk_pim2->p4(mpi);
	p_upim1->SetPxPyPzE(pim1_un.px(), pim1_un.py(), pim1_un.pz(), pim1_un.e());
	p_upim2->SetPxPyPzE(pim2_un.px(), pim2_un.py(), pim2_un.pz(), pim2_un.e());
	

	

	for(int i=0; i < nGamma-1 ; i++){
		for(int j=i+1; j < nGamma ; j++){
			RecEmcShower* gammaTrk1 = (*(evtRecTrkCol->begin()+iGamma[i]))->emcShower();
			RecEmcShower* gammaTrk2 = (*(evtRecTrkCol->begin()+iGamma[j]))->emcShower();
			





	//*******************************Begin Vetex Fit***********************************************************
	/*
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
	vtxfit->AddVertex(0, vxpar,0,1);
	if(!vtxfit->Fit(0)) return StatusCode::SUCCESS;
	vtxfit->Swim(0);
	vtxchisq=vtxfit->chisq(0);
	WTrackParameter wepTrk = vtxfit->wtrk(0);
	WTrackParameter wemTrk = vtxfit->wtrk(1);
	Ncut4++;   // finish vertex fit
	
	//#******************************Begin 4C  Fit after Vtxfit***********************************************************
	KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
	TLorentzVector p_jpsi(p4_cms.px(),p4_cms.py(),p4_cms.pz(),p4_cms.e());	
	chi4C = 9999.;
	double temp_chi4C = 999.;
	int idxgam1(-1),idxgam2(-1);

	for(int i=0; i < nGamma-1 ; i++){
		for(int j=i+1; j < nGamma ; j++){
			RecEmcShower* gammaTrk1 = (*(evtRecTrkCol->begin()+iGamma[i]))->emcShower();
			RecEmcShower* gammaTrk2 = (*(evtRecTrkCol->begin()+iGamma[j]))->emcShower();
			kmfit->init();
			kmfit->setEspread(m_EnergySpread);
			kmfit->setBeamPosition(xorigin);
			kmfit->setVBeamPosition(Evx_avg);
			kmfit->AddTrack(0,wepTrk);
			kmfit->AddTrack(1,wemTrk);
			kmfit->AddTrack(2,0.0,gammaTrk1);
			kmfit->AddTrack(3,0.0,gammaTrk2);
			kmfit->AddFourMomentum(0,p4_cms);

			bool okfit = kmfit->Fit();
			if(!okfit) continue;

			temp_chi4C = kmfit->chisq();

			if(temp_chi4C > chi4C) continue;
			chi4C=temp_chi4C;

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
			
			idxgam1=i;
			idxgam2=j;

			HepLorentzVector ep_ = kmfit->pfit(0);
			HepLorentzVector em_ = kmfit->pfit(1);
			HepLorentzVector gamma1_ = kmfit->pfit(2);
			HepLorentzVector gamma2_ = kmfit->pfit(3);

			p_ep->SetPxPyPzE(ep_.px(), ep_.py(), ep_.pz(), ep_.e());
			p_em->SetPxPyPzE(em_.px(), em_.py(), em_.pz(), em_.e());
			p_gamma1->SetPxPyPzE(gamma1_.px(), gamma1_.py(), gamma1_.pz(), gamma1_.e());
			p_gamma2->SetPxPyPzE(gamma2_.px(), gamma2_.py(), gamma2_.pz(), gamma2_.e());

			p_ugamma1->SetPxPyPzE((pGamma[i]).px(), (pGamma[i]).py(), (pGamma[i]).pz(), (pGamma[i]).e());
			p_ugamma2->SetPxPyPzE((pGamma[j]).px(), (pGamma[j]).py(), (pGamma[j]).pz(), (pGamma[j]).e());
		}
	}
	
	*p_ee=*p_ep+*p_em;
	*p_gg=*p_gamma1+*p_gamma2;
	*p_uee=*p_uep+*p_uem;
	*p_ugg=*p_ugamma1+*p_ugamma2;
	mee=p_ee->M();
	mgg=p_gg->M();
	
	p_gamma1b->SetPxPyPzE(p_gamma1->Px(),p_gamma1->Py(),p_gamma1->Pz(),p_gamma1->E());
	p_gamma2b->SetPxPyPzE(p_gamma2->Px(),p_gamma2->Py(),p_gamma2->Pz(),p_gamma2->E());
	p_gamma1b->Boost(-p_gg->BoostVector());	
	p_gamma2b->Boost(-p_gg->BoostVector());

	double angetag1 = p_gg->Angle(p_gamma1b->Vect()); //rad
	cosgam1b=cos(angetag1);
	
	double angetag2 = p_gg->Angle(p_gamma2b->Vect()); //rad
	cosgam2b=cos(angetag2);

	cosgamb=(p_gamma1->E()-p_gamma2->E())/p_gg->P();
	if(cosgamb<0) cosgamb=-cosgamb;

//	TRandom3 trd;
//	cosgamb=trd.Uniform(0,1)>0.5?cosgam2b:cosgam1b;
//	Ncut5++;


	//Gamma Conv
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
	*/

	TreeAna->Fill();
	GammaAll->Clear();

	p_ep->Clear();
	p_em->Clear();
	p_gamma1->Clear();
	p_gamma2->Clear();
	p_ee->Clear();
	p_gg->Clear();
	p_gamma1b->Clear();
	p_gamma2b->Clear();
	
	xep->Clear();
	xem->Clear();
	xgamma1->Clear();
	xgamma2->Clear();
	xee->Clear();
	xeta->Clear();
	xgg->Clear();

	p_uep->Clear();
	p_uem->Clear();
	p_ugamma1->Clear();
	p_ugamma2->Clear();
	p_uee->Clear();
	p_ugg->Clear();

	return StatusCode::SUCCESS;
}

//***************************************************************************

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode XisP3Pi::finalize() {
	saveFile->cd();
	TreeAna->Write();
	if(m_saveCutFlow == 1) {NbInfo->Fill();NbInfo->Write();}
	saveFile->Close();
	cout<<endl<<"Finalize psi'->e+ e- eta/pi0 (gam gam)"<<endl<<endl;;
	cout<<"Total number:                         "<<Ncut0<<endl;
	cout<<"nGood==2, nCharge==0:                 "<<Ncut1<<endl;
	cout<<"nGamma>=2&&nGamma<=10:                "<<Ncut2<<endl;
	cout<<"PID :                                 "<<Ncut3<<endl;
	cout<<"vertex fit:                           "<<Ncut4<<endl;
	cout<<"4C success:                           "<<Ncut5<<endl;
	cout<<"_END_TAG_"<<endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}


// corset(): sets up the generation by calculating C from V.
void XisP3Pi::corset(HepSymMatrix &V, HepMatrix &C, int n)
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
void XisP3Pi::corgen(HepMatrix &C, HepVector &x, int n)
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
void XisP3Pi::calibration(RecMdcKalTrack *trk , HepVector &wtrk_zHel, int n )
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
		XisP3Pi::corset(wpip_zcal,wpip_zerrc,3);
		HepVector wpip_zgen(3,0);
		XisP3Pi::corgen(wpip_zerrc,wpip_zgen,3);

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
		XisP3Pi::corset(wpim_zcal,wpim_zerrc,3);
		HepVector wpim_zgen(3,0);
		XisP3Pi::corgen(wpim_zerrc,wpim_zgen,3);
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
		XisP3Pi::corset(wep_zcal,wep_zerrc,3);
		HepVector wep_zgen(3,0);
		XisP3Pi::corgen(wep_zerrc,wep_zgen,3);

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
		XisP3Pi::corset(wem_zcal,wem_zerrc,3);
		HepVector wem_zgen(3,0);
		XisP3Pi::corgen(wem_zerrc,wem_zgen,3);
		wtrk_zHel[0] = trk->getZHelixE()[0];
		wtrk_zHel[1] = trk->getZHelixE()[1]+em_calmean_d2[1]*sqrt(wem_zerr[1][1])+wem_zgen[0];
		wtrk_zHel[2] = trk->getZHelixE()[2]+em_calmean_d2[2]*sqrt(wem_zerr[2][2])+wem_zgen[1];
		wtrk_zHel[3] = trk->getZHelixE()[3];
		wtrk_zHel[4] = trk->getZHelixE()[4]+em_calmean_d2[4]*sqrt(wem_zerr[4][4])+wem_zgen[2];

	}
}
