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
#include "VertexFit/SecondVertexFit.h"

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
		declareProperty("saveMCTruth", m_saveMCTruth = 0);//need to be re-evaluated when running different samples(only 1 for exclusiveMC)
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
	p_prop= new TClonesArray("TLorentzVector");
	p_lampim = new TClonesArray("TLorentzVector");
	p_xipim = new TClonesArray("TLorentzVector");
	p_lambda = new TClonesArray("TLorentzVector");
	p_xi = new TClonesArray("TLorentzVector");
	p_xis = new TClonesArray("TLorentzVector");
	p_xibar = new TClonesArray("TLorentzVector");

	//Lorentz Vector 
	p_pi0=new TLorentzVector();
	p_gamma1=new TLorentzVector();
	p_gamma2=new TLorentzVector();

	xp_xis=new TLorentzVector();
	xp_xi=new TLorentzVector();
	xp_xipim=new TLorentzVector();
	xp_lam=new TLorentzVector();
	xp_lampim=new TLorentzVector();
	xp_prop=new TLorentzVector();
	xp_pi0=new TLorentzVector();
	xp_gamma1=new TLorentzVector();
	xp_gamma2=new TLorentzVector();

	xp_xibar=new TLorentzVector();
	xp_lambar=new TLorentzVector();
	xp_xibarpip=new TLorentzVector();
	xp_prom=new TLorentzVector();
	xp_lambarpip=new TLorentzVector();

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

	TreeAna->Branch("flag_pim", flag_pim, "flag_pim[2]/D");
	TreeAna->Branch("flag_prop", &flag_prop, "flag_prop/D");

	TreeAna->Branch("p_pi0",&p_pi0,32000,0);
	TreeAna->Branch("p_gamma1",&p_gamma1,32000,0);
	TreeAna->Branch("p_gamma2",&p_gamma2,32000,0);
	TreeAna->Branch("chi1C", &chi1C, "chi1C/D");
	TreeAna->Branch("idx_gamma1", &idx_gamma1, "idx_gamma1/I");
	TreeAna->Branch("idx_gamma2", &idx_gamma2, "idx_gamma2/I");

	TreeAna->Branch("nVtxok", &nVtxok, "nVtxok/I");
	TreeAna->Branch("flag_vtx", &flag_vtx, "flag_vtx/I");

	TreeAna->Branch("chi_vtxLam", chi_vtxLam , "chi_vtxLam[nVtxok]/D");
	TreeAna->Branch("chi_vtxXi", chi_vtxXi , "chi_vtxXi[nVtxok]/D");
	TreeAna->Branch("chi_svtxXi", chi_svtxXi , "chi_svtxXi[nVtxok]/D");
	TreeAna->Branch("chi_vtx", chi_vtx , "chi_vtx[nVtxok]/D");
	TreeAna->Branch("m_lambda", m_lambda , "m_lambda[nVtxok]/D");
	TreeAna->Branch("m_xi", m_xi , "m_xi[nVtxok]/D");
	TreeAna->Branch("m_xis", m_xis , "m_xis[nVtxok]/D");
	TreeAna->Branch("m_xibar", m_xibar , "m_xibar[nVtxok]/D");
	TreeAna->Branch("m_ctau", m_ctau , "m_ctau[nVtxok]/D");
	TreeAna->Branch("m_len", m_len , "m_len[nVtxok]/D");
	TreeAna->Branch("m_lenerr", m_lenerr , "m_lenerr[nVtxok]/D");
	TreeAna->Branch("flag_lampim", flag_lampim , "flag_lampim[nVtxok]/D");
	TreeAna->Branch("flag_xipim", flag_xipim , "flag_xipim[nVtxok]/D");

	TreeAna->Branch("p_prop", "TClonesArray", &p_prop, 256000, 0);
	TreeAna->Branch("p_lampim", "TClonesArray", &p_lampim, 256000, 0);
	TreeAna->Branch("p_xipim", "TClonesArray", &p_xipim, 256000, 0);
	TreeAna->Branch("p_lambda", "TClonesArray", &p_lambda, 256000, 0);
	TreeAna->Branch("p_xi", "TClonesArray", &p_xi, 256000, 0);
	TreeAna->Branch("p_xis", "TClonesArray", &p_xis, 256000, 0);
	TreeAna->Branch("p_xibar", "TClonesArray", &p_xibar, 256000, 0);

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
	p_prop->BypassStreamer();
	p_lampim->BypassStreamer();
	p_xipim->BypassStreamer();
	p_lambda->BypassStreamer();
	p_xi->BypassStreamer();
	p_xis->BypassStreamer();
	p_xibar->BypassStreamer();

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

////////////////////////////Global Event Parameters
Vint iGood, iGamma;
Vint iTrp,iTrm,iprop,iprom,ipip,ipim,iem;
iGood.clear();
int nGood = 0;
int nCharge = 0;
nGamma = 0;

int nTrp=0;
int nTrm=0;

iGood.clear();
iGamma.clear();
iTrp.clear();
iTrm.clear();
////////////////////////////////////primary vertex
Hep3Vector xorigin(0,0,0);
//	HepSymMatrix Evx_avg(3,0);
HepSymMatrix VtxErr(3,0);

IVertexDbSvc*  vtxsvc;
Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
if(vtxsvc->isVertexValid()){
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

////////////////////////good track
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
if(nGood >12) return StatusCode::SUCCESS;
if(nTrp<1||nTrm<2)  return StatusCode::SUCCESS;
Ncut1++;		//after good track selection

//////////////////////////////////////good photon select
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

if(nGamma < m_NGamma || nGamma >15) return StatusCode::SUCCESS; 
Ncut2++;	//after good gamma cut

/////////////////////////////////assign 4-momentum for photons
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
////////////////////////////////// 1C for pi0
double min_chi1C=9999.; chi1C=-99;idx_gamma1=-1;idx_gamma2=-1;
KalmanKinematicFit *kmfit = KalmanKinematicFit::instance();
		HepLorentzVector hv_pi0(0,0,0,0);

for(int i=0; i < nGamma-1 ; i++){
	for(int j=i+1; j < nGamma ; j++){
		RecEmcShower* gammaTrk1 = (*(evtRecTrkCol->begin()+iGamma[i]))->emcShower();
		RecEmcShower* gammaTrk2 = (*(evtRecTrkCol->begin()+iGamma[j]))->emcShower();
		kmfit->init();
		kmfit->setBeamPosition(xorigin);
		kmfit->setVBeamPosition(VtxErr);
		kmfit->AddTrack(0,0.0,gammaTrk1);
		kmfit->AddTrack(1,0.0,gammaTrk2);
		kmfit->AddResonance(0, 0.135, 0, 1);
		bool okfit = kmfit->Fit();
		if(!okfit) continue;
		if(kmfit->chisq()>min_chi1C) continue;

		chi1C=kmfit->chisq();
		idx_gamma1=i;
		idx_gamma2=j;

		HepLorentzVector gamma1_ = kmfit->pfit(0);
		HepLorentzVector gamma2_ = kmfit->pfit(1);
		hv_pi0 = gamma1_ + gamma2_;

		p_gamma1->SetPxPyPzE(gamma1_.px(),gamma1_.py(),gamma1_.pz(),gamma1_.e());
		p_gamma2->SetPxPyPzE(gamma2_.px(),gamma2_.py(),gamma2_.pz(),gamma2_.e());
		p_pi0->SetPxPyPzE(hv_pi0.px(),hv_pi0.py(),hv_pi0.pz(),hv_pi0.e());
	}
}
if (chi1C==-99)  return StatusCode::SUCCESS;
Ncut3++;		//after 1C

////////////////////////////////PID
flag_prop=-1;flag_pim[0]=-1;flag_pim[1]=-1;
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
	double	prob_pi = pid->probPion();
	double prob_K = pid->probKaon();
	double prob_pro = pid->probProton();
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
Ncut4++;		//after pid

//////////////////// Assign track  info
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
wvpropTrk = WTrackParameter(mp, mdcTrk_prop->getZHelixP(), mdcTrk_prop->getZErrorP());  

//	pu_prop->SetPxPyPzE(prop_un.px(), prop_un.py(), prop_un.pz(), prop_un.e());

RecMdcKalTrack::setPidType(RecMdcKalTrack::pion);
mdcTrk_pim[0] = (*(evtRecTrkCol->begin()+iGood[flag_pim[0]]))->mdcKalTrack();
mdcTrk_pim[1] = (*(evtRecTrkCol->begin()+iGood[flag_pim[1]]))->mdcKalTrack();

wvpimTrk[0] = WTrackParameter(mpi, mdcTrk_pim[0]->getZHelix(), mdcTrk_pim[0]->getZError());  
wvpimTrk[1] = WTrackParameter(mpi, mdcTrk_pim[1]->getZHelix(), mdcTrk_pim[1]->getZError());  

hvun_pim[0] = mdcTrk_pim[0]->p4(mpi);
hvun_pim[1] = mdcTrk_pim[1]->p4(mpi);

//	pu_pim1->SetPxPyPzE(hvun_pim[0].px(), hvun_pim[0].py(), hvun_pim[0].pz(), hvun_pim[0].e());
//	pu_pim2->SetPxPyPzE(hvun_pim[1].px(), hvun_pim[1].py(), hvun_pim[1].pz(), hvun_pim[1].e());

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

nVtxok=0;
flag_vtx=-1;
double min_diffmxi=9999.;

for(int i1=0;i1<2;i1++){
	//recontruct Lambda	
	VertexFit* vtxfit1 = VertexFit::instance();
	HepLorentzVector VLambda(0,0,0,0);
	VLambda=hvun_pim[i1]+hvun_prop;
	if( fabs(VLambda.m()-mLambda)>0.020) continue;

	vtxfit1->init();
	vtxfit1->AddTrack(0,  wvpropTrk);
	vtxfit1->AddTrack(1,  wvpimTrk[i1]);
	vtxfit1->AddVertex(0, vxpar, 0, 1);
	if( !(vtxfit1->Fit(0)) ) continue;
	vtxfit1->Swim(0);
	vtxfit1->BuildVirtualParticle(0);


	HepLorentzVector hv_prop = vtxfit1->pfit(0);
	HepLorentzVector hv_lampim = vtxfit1->pfit(1);         
	HepLorentzVector hv_lambda= hv_prop+hv_lampim; 

	chi_vtxLam[nVtxok]=vtxfit1->chisq();
	WTrackParameter	wvLambda=vtxfit1->wVirtualTrack(0);

	//recontruct Xi-
	int another_pim=(int)abs(1-i1);
	VertexFit* vtxfit2 = VertexFit::instance();
	vtxfit2->init();
	vtxfit2->AddTrack(0, wvLambda);
	vtxfit2->AddTrack(1, wvpimTrk[another_pim]);
	vtxfit2->AddVertex(0, vxpar, 0, 1);
	if( !(vtxfit2->Fit(0)) ) continue;
	vtxfit2->Swim(0);

	vtxfit2->BuildVirtualParticle(0);
	chi_vtxXi[nVtxok]=vtxfit2->chisq();

	HepLorentzVector hv_xipim = vtxfit2->pfit(1);         

	WTrackParameter wvXi=vtxfit2->wVirtualTrack(0);
	VertexParameter vparXi = vtxfit2->vpar(0);
	//secondvertexfit for Xi-
	SecondVertexFit *svtxfit = SecondVertexFit::instance();
	svtxfit->init();
	svtxfit->setPrimaryVertex(bs);
	svtxfit->AddTrack(0, vtxfit2->wVirtualTrack(0));
	svtxfit->setVpar(vparXi);
	if (!svtxfit->Fit()) continue;

	chi_svtxXi[nVtxok]=svtxfit->chisq();
	m_ctau[nVtxok] = svtxfit->ctau();
	m_len[nVtxok] = svtxfit->decayLength();
	m_lenerr[nVtxok] = svtxfit->decayLengthError();

	HepLorentzVector hv_xi = svtxfit->p4par();
	HepLorentzVector hv_xis = hv_xi+hv_pi0;
	HepLorentzVector hv_xibar = p4_cms - hv_xis;

	TLorentzVector tl_prop, tl_lampim,tl_xipim, tl_lambda, tl_xi,tl_xis, tl_xibar;
	tl_prop.SetPxPyPzE(hv_prop.px(),hv_prop.py(),hv_prop.pz(),hv_prop.e());
	tl_lampim.SetPxPyPzE(hv_lampim.px(),hv_lampim.py(),hv_lampim.pz(),hv_lampim.e());
	tl_xipim.SetPxPyPzE(hv_xipim.px(),hv_xipim.py(),hv_xipim.pz(),hv_xipim.e());
	tl_lambda.SetPxPyPzE(hv_lambda.px(),hv_lambda.py(),hv_lambda.pz(),hv_lambda.e());
	tl_xi.SetPxPyPzE(hv_xi.px(),hv_xi.py(),hv_xi.pz(),hv_xi.e());
	tl_xis.SetPxPyPzE(hv_xis.px(),hv_xis.py(),hv_xis.pz(),hv_xis.e());
	tl_xibar.SetPxPyPzE(hv_xibar.px(),hv_xibar.py(),hv_xibar.pz(),hv_xibar.e());

	new ((*p_prop)[nVtxok]) TLorentzVector(tl_prop);
	new ((*p_lampim)[nVtxok]) TLorentzVector(tl_lampim);
	new ((*p_xipim)[nVtxok]) TLorentzVector(tl_xipim);
	new ((*p_lambda)[nVtxok]) TLorentzVector(tl_lambda);
	new ((*p_xi)[nVtxok]) TLorentzVector(tl_xi);
	new ((*p_xis)[nVtxok]) TLorentzVector(tl_xis);
	new ((*p_xibar)[nVtxok]) TLorentzVector(tl_xibar);

	m_lambda[nVtxok]=hv_lambda.m();
	m_xi[nVtxok]=hv_xi.m();
	m_xis[nVtxok]=hv_xis.m();
	m_xibar[nVtxok]=hv_xibar.m();
	chi_vtx[nVtxok]=chi_vtxLam[nVtxok]+chi_vtxXi[nVtxok]+chi_svtxXi[nVtxok];

	flag_lampim[nVtxok]=flag_pim[i1];
	flag_xipim[nVtxok]=flag_pim[another_pim];

	double t_diffmxi=fabs(hv_xi.m()-mXi);
	if(t_diffmxi < min_diffmxi) {
		min_diffmxi=t_diffmxi;
		flag_vtx=nVtxok;
	}
	nVtxok++;
}
if(nVtxok==0) return StatusCode::SUCCESS;
//end vertex fit

TreeAna->Fill();
GammaAll->Clear();
p_prop->Clear();
p_lampim->Clear();
p_xipim->Clear();
p_lambda->Clear();
p_xi->Clear();
p_xis->Clear();
p_xibar->Clear();

p_pi0->Clear();
p_gamma1->Clear();
p_gamma2->Clear();

return StatusCode::SUCCESS;
}

//***************************************************************************

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode XisP3Pi::finalize() {
	saveFile->cd();
	TreeAna->Write();
	if(m_saveCutFlow == 1) {NbInfo->Fill();NbInfo->Write();}
	saveFile->Close();
	cout<<endl<<"Finalize psi'-> Xi*-(p+ pi- pi-) + Xi_bar +"<<endl<<endl;;
	cout<<"Total number:                         "<<Ncut0<<endl;
	cout<<"nGood<=12,nPos>=1,nNeg>=2:                 "<<Ncut1<<endl;
	cout<<"nGamma>=2&&nGamma<=15:                "<<Ncut2<<endl;
	cout<<"PID :                                 "<<Ncut3<<endl;
	cout<<"vertex fit:                           "<<Ncut4<<endl;
	cout<<"4C success:                           "<<Ncut5<<endl;
	cout<<"_ENDTAG_"<<endl;
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in finalize()" << endmsg;
	return StatusCode::SUCCESS;
}
