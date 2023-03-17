#ifndef TOOLS_H
#define TOOLS_H

#include "ROOT/RVec.hxx"

ROOT::VecOps::RVec<float> dilepton_vec(const float pt0, const float eta0,
        const float phi0, const float mass0, const float pt1,
        const float eta1, const float phi1, const float mass1);

ROOT::VecOps::RVec<float> cos_rapidity(const float pt0, const float eta0,
        const float phi0, const float mass0, const int charge0, const float pt1,
        const float eta1, const float phi1, const float mass1);

ROOT::VecOps::RVec<float> weights(const float pt, const float m, const float c);

#endif
