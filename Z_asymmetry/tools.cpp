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
    // and mass of each muon of the filtered events and it returns a "VecOps"
    // with the value of the cos(theta*) in the Collins-Soper frame and the
    // rapidity of the dileptons.

ROOT::VecOps::RVec<float> cos_rapidity(const float pt0, const float eta0,
        const float phi0, const float mass0, const float pt1,
        const float eta1, const float phi1, const float mass1)
    {
      float pz0, pz1, E0, E1, P0_1, P0_2, P1_1, P1_2, mll, ptt, pzll, cos, y;

      ROOT::Math::PtEtaPhiMVector p0, p1;
      p0.SetCoordinates(pt0, eta0, phi0, mass0);
      p1.SetCoordinates(pt1, eta1, phi1, mass1);

      pz0 = pt0*sinh(eta0);
      E0 = sqrt(pow(pz0,2)+pow(pt0,2)+pow(mass0,2));
      pz1 = pt1*sinh(eta1);
      E1 = sqrt(pow(pz1,2)+pow(pt1,2)+pow(mass1,2));

      P0_1 = (E0+pz0)/sqrt(2);
      P0_2 = (E0-pz0)/sqrt(2);
      P1_1 = (E1+pz1)/sqrt(2);
      P1_2 = (E1-pz1)/sqrt(2);

      mll = pow((p1+p0).M(),2);
      ptt = pow((p1+p0).Pt(),2);
      pzll = pz0+pz1;

      //numer = (2*((P0_1*P1_2) - (P0_2*P1_1)));
      //denom = sqrt(mll*(mll+ptt));
      //cos = (numer/denom)*(pzll/abs(pzll));

      cos = ((2*((P0_1*P1_2) - (P0_2*P1_1)))/(sqrt(mll*(mll+ptt))))*(pzll/abs(pzll));
      y = (0.5)*log(((E0+E1)+(pzll))/((E0+E1)-(pzll)));

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
        //A0 = 0.1;

        A0 = (pow(pt,2))/(pow(m, 2)+pow(pt,2));
        h = (0.5)*A0*(1-(3*pow(c,2)));
        w_d = (0.5)*(pow(c,2)/pow(1+pow(c,2)+h,3));
        w_n = (0.5)*(abs(c)/pow(1+pow(c,2)+h,2));

        ROOT::VecOps::RVec<float> P{w_d,w_n};

        return P;
        }
