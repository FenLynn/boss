#ifndef Physics_Analysis_GammaConv_H
#define Physics_Analysis_GammaConv_H

#include "VertexFit/HTrackParameter.h"
#include "VertexFit/WTrackParameter.h"
#include "CLHEP/Matrix/Vector.h"
using CLHEP::HepVector;
#ifndef CLHEP_THREEVECTOR_H
#include "CLHEP/Vector/ThreeVector.h"
using CLHEP::Hep3Vector;
#endif
#include "CLHEP/Vector/LorentzVector.h"
using CLHEP::HepLorentzVector;
#ifndef CLHEP_POINT3D_H
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#endif

class GammaConv
{

public:

  GammaConv(HepVector helixp, HepVector helixe, HepPoint3D ip);
  ~GammaConv();

public:
  //get the four parameters
  double getDeltaXY() {return m_deltaxy;}
  double getDeltaZ1()  {return m_deltaz1;}
  double getDeltaZ2()  {return m_deltaz2;}

  double getDeltaZ()  {return m_deltaz;}

  double getDGamma()  {return m_dgamma;}
  double getRX1() {return m_rx1;}
  double getRY1() {return m_ry1;}
  double getRZ1() {return m_rz1;}
  double getRXY1(){return m_rxy1;}

  double getRX() {return m_rx;}
  double getRY() {return m_ry;}
  double getRZ() {return m_rz;}
  double getRXY(){return m_rxy;}

  double getRX2() {return m_rx2;}
  double getRY2() {return m_ry2;}
  double getRZ2() {return m_rz2;}
  double getRXY2(){return m_rxy2;}

  double getR()  {return m_R;}
  double getCthe() {return m_cthe;}
  void GetCoordi1(HepPoint3D &coordip1, HepPoint3D &coordie1);
  void GetCoordi2(HepPoint3D &coordip2, HepPoint3D &coordie2);
  void GetCoordi(HepPoint3D &coordip, HepPoint3D &coordie);
  HepPoint3D getXv1() {return m_xv1;}
  HepPoint3D getXv2() {return m_xv2;}
  HepPoint3D getXv() {return m_xv;}
  void getPhi0(double &phi0p, double &phi0e){ phi0p = m_phi0p; phi0e = m_phi0e; }
  void getTheta0(double &theta0p, double &theta0e){ theta0p = m_theta0p; theta0e = m_theta0e; }
  double getDeltaphi0() {return m_deltaphi0;}
  double getDeltatheta0() {return m_deltatheta0;}
  double getXiep() {return m_xiep;}
  double getRp() {return rp;}
  double getEp() {return re;}
  double getDeltaeq() {return delta_eq;}
  double getNcase(){return Ncase;}
  double getLep(){return lep;}

  double getPsipair() {return m_psipair;}
  HepPoint3D center(HepVector helix);
  double radius(HepVector helix);
private:
  void DeltaXY(HepVector helixp, HepVector helixe);
  void DeltaZ(HepVector helixp, HepVector helixe);
  void DGamma(HepVector helixp, HepVector helixe, HepPoint3D ip);
  void GGeom(HepPoint3D ip);
  void Deltaphitheta(HepVector helixp, HepVector helixe);
  Hep3Vector p3(HepVector helix);
  HepPoint3D x3(HepVector helix);
private:
  double m_deltaxy;
  double m_deltaz1, m_deltaz2;

  double m_deltaz;
  double m_dgamma;
  double m_rx1;
  double m_ry1;
  double m_rz1;
  double m_rxy1;

  double m_rx2;
  double m_ry2;
  double m_rz2;
  double m_rxy2;


  double m_rx;
  double m_ry;
  double m_rz;
  double m_rxy;

  double m_R;
  double m_cthe;
  double rp,re;
  double delta_eq;
  double Ncase;
  double lep;
 
  HepPoint3D _coordip1, _coordie1;
  HepPoint3D _coordip2, _coordie2;
  HepPoint3D _coordip, _coordie;
  HepPoint3D m_xv1,m_xv2;
  HepPoint3D m_xv;

  double m_phi0p;
  double m_phi0e;
  double m_theta0p;
  double m_theta0e;
  double m_deltaphi0;
  double m_deltatheta0;
  double m_xiep;
  double m_psipair;
};
#endif
