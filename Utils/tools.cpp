#include "tools.h"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"


  // The function takes in input transverse momentum, rapidity, azimuthal angle
  // and mass of each muon of the filtered events and it returns a "VecOps"
  // with the same physical quantities but computed for the dileptons system.

ROOT::VecOps::RVec<float> dilepton_vec(const float pt0, const float eta0,
        const float phi0, const float mass0, const float pt1,
        const float eta1, const float phi1, const float mass1)
    {
      ROOT::Math::PtEtaPhiMVector p0, p1;
      p0.SetCoordinates(pt0, eta0, phi0, mass0);
      p1.SetCoordinates(pt1, eta1, phi1, mass1);

      ROOT::VecOps::RVec<double> P{(p1+p0).Pt(), (p1+p0).Eta(),
          (p1+p0).Phi(), (p1+p0).M()};

      return P;
    };


    // The function takes in input transverse momentum, rapidity, azimuthal angle
    // and mass of each muon and the charge of one of them of the filtered events
    // and it returns a "VecOps" with the value of the cos(theta*) in the
    // Collins-Soper frame and the rapidity of the dileptons.

ROOT::VecOps::RVec<float> cos_rapidity(const float pt0, const float eta0,
        const float phi0, const float mass0, const int charge0, const float pt1,
        const float eta1, const float phi1, const float mass1)
    {
      float pz1, pz2, E1, E2, P1_p, P1_m, P2_p, P2_m, mll, ptt, pzll, cos, y;
      ROOT::Math::PtEtaPhiMVector pp, pm;
      std::vector<double> pl1, pl2, p1, p2;

      // Set the array
      pp.SetCoordinates(pt0, eta0, phi0, mass0);
      pm.SetCoordinates(pt1, eta1, phi1, mass1);

      pl1.emplace_back(pt0);
      pl1.emplace_back(eta0);
      pl1.emplace_back(phi0);
      pl1.emplace_back(mass0);

      pl2.emplace_back(pt1);
      pl2.emplace_back(eta1);
      pl2.emplace_back(phi1);
      pl2.emplace_back(mass1);

      // p1 has to be the muons with charge -1.
      if(charge0==-1){
        p1=pl1;
        p2=pl2;
      }
      else{
        p1=pl2;
        p2=pl1;
      }

      pz1 = p1[0]*sinh(p1[1]);
      E1 = sqrt(pow(pz1,2)+pow(p1[0],2)+pow(p1[3],2));
      pz2 = p2[0]*sinh(p2[1]);
      E2 = sqrt(pow(pz2,2)+pow(p2[0],2)+pow(p2[3],2));

      P1_p = (E1+pz1)/sqrt(2);
      P1_m = (E1-pz1)/sqrt(2);
      P2_p = (E2+pz2)/sqrt(2);
      P2_m = (E2-pz2)/sqrt(2);

      mll = pow((pp+pm).M(),2);
      ptt = pow((pp+pm).Pt(),2);
      pzll = pz1+pz2;

      cos = ((2*((P1_p*P2_m) - (P1_m*P2_p)))/(sqrt(mll*(mll+ptt))))*(pzll/abs(pzll));
      y = (0.5)*log(((E1+E2)+(pzll))/((E1+E2)-(pzll)));

      ROOT::VecOps::RVec<float> A{cos,y};

      return A;
    };


  // The function takes in input pt, mass and cos of eah dimuon and it
  // returns a "VecOps" with useful quantities to compute the forward-backward
  // asymmetry. In particular they are "w_d" and "w_n" (where "d" and "n"
  // stands for "denominator" and "numerator" weight, respectively).


ROOT::VecOps::RVec<float> weights(const float pt, const float m, const float c)
      {
        float A0, h, w_d, w_n;

        A0 = (pow(pt,2))/(pow(m, 2)+pow(pt,2));
        h = (0.5)*A0*(1-(3*pow(c,2)));
        w_d = (0.5)*(pow(c,2)/pow(1+pow(c,2)+h,3));
        w_n = (0.5)*(abs(c)/pow(1+pow(c,2)+h,2));

        ROOT::VecOps::RVec<float> P{w_d,w_n};

        return P;
        }
