#include <TROOT.h>
#include <TRandom.h>

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>


double Al = 0.; //crossing angle x
double SiThX = 2.66e-6;//rad, beam divergence x
double SiThY = 2.2e-6;//rad, beam divergence y
double MeanXi = 0; //beam momentum offset
double SiXi = 1.e-4;//momentum spread of the beam
double E_CMS = 4000;// proton energy
double MeanX = 0; //mm, vertex position x
double SiX = 0.000169;//m, vertex size x
double MeanY = 0;//vertex position y
double SiY = 0.00014;//m, vertex size y
double MeanZ = 0.;//vertex position z
double SiZ = 0.070;//m, vertex size z

void beam_smearing(double proton_px, double proton_py, double proton_pz, double proton_energy, double &proton_px_smear, double &proton_py_smear, double &proton_pz_smear, double &proton_energy_smear){

      // generate energy/angular smearing
      TRandom rand;
      double al1 = Al;
      double al2 = -Al;
      double thx1 = rand.Gaus(0., SiThX);
      double thy1 = rand.Gaus(0., SiThY);
      double thx2 = rand.Gaus(0., SiThX);
      double thy2 = rand.Gaus(0., SiThY);
      double xi1 = rand.Gaus(MeanXi, SiXi);
      double xi2 = rand.Gaus(MeanXi, SiXi);

      // compute transform parameters
      double m = 0.938271;
      double p_nom = sqrt(E_CMS * E_CMS - m*m);

      double thz1 = sqrt(1. - thx1*thx1 - thy1*thy1);
      double thz2 = sqrt(1. - thx2*thx2 - thy2*thy2);

      TVector3 p1(cos(al1)*thx1 + sin(al1)*thz1, thy1, -sin(al1)*thx1 + cos(al1)*thz1);
      p1 *= (1. + xi1) * p_nom;
      double E1 = sqrt(p1.Mag()*p1.Mag() + m*m);
      TLorentzVector P1(p1.x(), p1.y(), p1.z(), E1);

      TVector3 p2(cos(al2)*thx2 + sin(al2)*thz2, thy2, -sin(al2)*thx2 + cos(al2)*thz2);
      p2 *= -(1. + xi2) * p_nom;
      double E2 = sqrt(p2.Mag()*p2.Mag() + m*m);
      TLorentzVector P2(p2.X(), p2.Y(), p2.Z(), E2);

      double factor = (P1 + P2).Mag() / 2. / E_CMS;             // energy correction factor

      TVector3 boost = (p1 + p2) * (1. / (E1 + E2));            // boost vector (direction * beta)
      double beta = (p1 + p2).Mag() / (E1 + E2);                // beta of the boost
      P1.Boost(-boost);
      TVector3 axis(P1.Y(), -P1.X(), 0.);                       // rotation axis
      double angle = -acos( P1.Z() / P1.Vect().Mag() );         // angle of rotation

      TLorentzVector p_proton(proton_px, proton_py, proton_pz, proton_energy); //fourmomentum proton

      // energy scaling
      TVector3 p(p_proton.X(), p_proton.Y(), p_proton.Z());
      double E = p_proton.Energy() * factor;
      m = p_proton.M();
      p = sqrt(E*E - m*m) / p.Mag() * p;
      TLorentzVector PP(p, E);

      // rotation
      if (fabs(angle) > 1E-8) PP.Rotate(angle, axis);

      // boost
      if (fabs(beta) > 1E-8) PP.Boost(boost);

      proton_px_smear = PP.X();
      proton_py_smear = PP.Y();
      proton_pz_smear = PP.Z();
      proton_energy_smear = PP.Energy();

}

void vtx_smearing(double &v_x, double &v_y, double &v_z){

     TRandom rand;
     v_x = rand.Gaus(MeanX, SiX);
     v_y = rand.Gaus(MeanY, SiY);
     v_z = rand.Gaus(MeanZ, SiZ);
     //v_t = 0.;
}
