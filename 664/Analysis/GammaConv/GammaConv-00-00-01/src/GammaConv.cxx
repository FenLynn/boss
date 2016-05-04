//#include "/ihepbatch/bes/zhangyt/workarea/GammaConv/GammaConv-00-00-01/GammaConv/GammaConv.h"
#include "GammaConv/GammaConv.h"

#include "VertexFit/BField.h"
#include <iostream>

#include <vector>
using namespace std;
typedef std::vector<double> Vdouble;

//const double alpha = -0.00299792458;
const double alpha = 0.00299792458;  ///yateng
const double emass = 0.000511;

GammaConv::GammaConv(HepVector helixp, HepVector helixe, HepPoint3D ip){

  DeltaXY(helixp, helixe);
  DeltaZ(helixp, helixe);
  GGeom(ip);
  DGamma(helixp, helixe, ip);
  Deltaphitheta(helixp, helixe);
//cout << "GammaConv is called" << endl;
}

GammaConv::~GammaConv(){
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

HepPoint3D GammaConv::x3(HepVector helix){
  return HepPoint3D(helix[0]*cos(helix[1]), helix[0]*sin(helix[1]), helix[3]);
}
double GammaConv::radius(HepVector helix){
// HepPoint3D tesx3(5,0,0);
//  double bField = VertexFitBField::instance()->getBFieldZ(tesx3);
  double bField = VertexFitBField::instance()->getBFieldZ(x3(helix));
// double  bField=-1;
  int charge = 1;
  double a = alpha * bField * charge;
  double pxy = charge/helix[2];
  return (pxy/a);
}
Hep3Vector GammaConv::p3(HepVector helix){
  double pxy = 1./fabs(helix[2]);
  return Hep3Vector(0-pxy*sin(helix[1]), pxy*cos(helix[1]), pxy*helix[4]);
}
HepPoint3D GammaConv::center(HepVector helix){
//   HepPoint3D tesx3(5,0,0);
//  double bField = VertexFitBField::instance()->getBFieldZ(tesx3);
  double bField = VertexFitBField::instance()->getBFieldZ(x3(helix));
//  double  bField=-1;
//  int charge = helix[2] > 0 ? +1: -1;
//  cout<<"bField   "<<bField<<endl;
  int charge = 1;
  double a = alpha * bField * charge;
  double pxy = charge/helix[2];
  double rad = pxy/a;
  double x = (helix[0] + rad) * cos(helix[1]);
  double y = (helix[0] + rad) * sin(helix[1]);
  double z = 0.0;
  return HepPoint3D(x, y, z);
}

void GammaConv::DeltaXY(HepVector helixp, HepVector helixe){

  double rp = fabs(radius(helixp));
  double re = fabs(radius(helixe));

  double xp, yp, xe, ye;
  xp = center(helixp).x();
  yp = center(helixp).y();
  xe = center(helixe).x();

   ye = center(helixe).y();
  double xa, ya, xb, yb;
  double lep;
  lep = sqrt((xp-xe)*(xp-xe)+(yp-ye)*(yp-ye));
  xa = xe+(xp-xe)*re/lep;
  
 
  ya = ye+(yp-ye)*re/lep;
  xb = xp-(xp-xe)*rp/lep;
  yb = yp-(yp-ye)*rp/lep;
  _coordie.setX(xa);
  _coordie.setY(ya);
  _coordip.setX(xb);
  _coordip.setY(yb);
  double len;
  len = sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb));
  if(lep >= rp+re) m_deltaxy = len;
  else m_deltaxy = -len;
}

void GammaConv::DeltaZ(HepVector helixp, HepVector helixe){

  HepPoint3D posp = x3(helixp);
  HepPoint3D pose = x3(helixe);
  double radp = radius(helixp);
  double rade = radius(helixe);

  double jp;
  Hep3Vector pp = p3(helixp);
  double bFieldp = VertexFitBField::instance()->getBFieldZ(posp);
  int    chargep = helixp[2] > 0 ? +1: -1;
  double aconp = alpha * bFieldp * chargep; 

  jp = aconp*((_coordip.x()-posp.x())*pp.x()+(_coordip.y()-posp.y())*pp.y())/(pp.x()*pp.x()+pp.y()*pp.y());
  double zp;
//  zp = helixp[4]+pp.z()*asin(jp)/aconp;
  zp = helixp[3]+pp.z()*asin(jp)/aconp;//change by wenjing yateng

  double je;
  Hep3Vector pe = p3(helixe);
  double bFielde = VertexFitBField::instance()->getBFieldZ(pose);
  int    chargee = helixe[2] > 0 ? +1: -1;
  double acone = alpha * bFielde * chargee;

  je = acone*((_coordie.x()-pose.x())*pe.x()+(_coordie.y()-pose.y())*pe.y())/(pe.x()*pe.x()+pe.y()*pe.y());
  double ze;
// cout<<"je   "<<je<<endl;
//  ze = helixe[4]+pe.z()*asin(je)/acone;
  ze = helixe[3]+pe.z()*asin(je)/acone; //change by wenjing yateng

  _coordip.setZ(zp);
  _coordie.setZ(ze);
  
  m_deltaz = fabs(zp-ze);
}
void GammaConv::DGamma(HepVector helixp, HepVector helixe, HepPoint3D ip){

  double eposi = sqrt(emass*emass+p3(helixp).mag2());
  double eelec = sqrt(emass*emass+p3(helixe).mag2());
  HepLorentzVector pposi(p3(helixp), eposi);
  HepLorentzVector pelec(p3(helixe), eelec);
  HepLorentzVector pgamma(pposi+pelec);
  double p = pgamma.vect().mag();
  double dipa;
//  dipa = (ip.x()-m_rx)*pgamma.px()/p+(ip.y()-m_ry)*pgamma.py()/p+(ip.z()-m_rz)*pgamma.pz()/p;
  dipa = (m_rx-ip.x())*pgamma.px()/p+(m_ry-ip.y())*pgamma.py()/p+(ip.z()-m_rz)*pgamma.pz()/p;//yateng 
  double dgp;
  dgp = sqrt((ip.x()-m_rx)*(ip.x()-m_rx)+(ip.y()-m_ry)*(ip.y()-m_ry)+(ip.z()-m_rz)*(ip.z()-m_rz));
  m_dgamma = sqrt(dgp*dgp-dipa*dipa);
  m_R = dgp;
  m_cthe = dipa/dgp;
}

void GammaConv::GGeom(HepPoint3D ip){
  double x, y, zz;
  HepPoint3D OA, OB, OM;
  OA = _coordip-ip;
  OB = _coordie-ip;
  OM = (OA+OB)*0.5;
  m_rx = OM.x();
// cout<<"m_rx  "<<m_rx<<endl;
   m_ry = OM.y();
  m_rz = OM.z();
  m_rxy= sqrt(OM.x()*OM.x()+OM.y()*OM.y());
  m_xv = OM;
}

void GammaConv::GetCoordi(HepPoint3D &coordip, HepPoint3D &coordie){
  coordip = _coordip;
  coordie = _coordie;
}

void GammaConv::Deltaphitheta(HepVector helixp, HepVector helixe){
  double deltaphi0 = acos(-sin(helixe[1]))-acos(-sin(helixp[1]));
  double phi0p = acos(-sin(helixp[1]));
  double phi0e = acos(-sin(helixe[1]));
  m_phi0p     = phi0p;
  m_phi0e     = phi0e;
  m_deltaphi0 = deltaphi0;
  double deltatheta0 = acos(helixe[4]/sqrt(1+helixe[4]*helixe[4]))-acos(helixp[4]/sqrt(1+helixp[4]*helixp[4]));
  double theta0p = acos(helixp[4]/sqrt(1+helixp[4]*helixp[4]));
  double theta0e = acos(helixe[4]/sqrt(1+helixe[4]*helixe[4]));
  m_theta0p     = theta0p;
  m_theta0e     = theta0e;
  m_deltatheta0 = deltatheta0;
  double xiep = acos(p3(helixp)*p3(helixe)/p3(helixp).mag()/p3(helixe).mag());
 //   m_xiep = xiep;
  double acol = acos(p3(helixp)*p3(helixe)/p3(helixp).mag()/p3(helixe).mag())* 180/3.14159265358979323846;
  m_xiep = acol;
  double psipair = asin(deltatheta0/xiep);
  m_psipair = psipair;
}

