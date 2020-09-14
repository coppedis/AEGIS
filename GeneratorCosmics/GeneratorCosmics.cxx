// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Generator for muons according to kinematic parametrizations at ALICE
// (not at the surface).
// Origin: andrea.dainese@lnl.infn.it
// Modifications for O2: ruben.shahoyan@cern.ch

#include <TRandom.h>
#include <TParticle.h>
#include <TVirtualMC.h>
#include <TGeoGlobalMagField.h>
#include "GeneratorCosmics.h"

//-----------------------------------------------------------------------------
GeneratorCosmics::GeneratorCosmics() : TGenerator("GeneratorCosmics", "GeneratorCosmics")
{
}

//-----------------------------------------------------------------------------
bool GeneratorCosmics::detectField()
{
  //
  auto fld = TGeoGlobalMagField::Instance()->GetField();
  if (TVirtualMC::GetMC() && TVirtualMC::GetMC()->GetMagField()) {
    fld = TVirtualMC::GetMC()->GetMagField();
  }
  if (fld) {
    double r[3] = {0, 0, 0}, b[3] = {0, 0, 0};
    fld->Field(r, b);
    setBkG(b[2]);
    return true;
  }
  return false;
}


//-----------------------------------------------------------------------------
void GeneratorCosmics::GenerateEvent()
{
  //
  // Generate muon(s)
  constexpr int MuMinusPDG = 13, MuPlusPDG = -13;
  constexpr float MuMass = 0.1056583;
  //
  if (!mFieldIsSet && !detectField()) {
    throw std::runtime_error("Failed to fetch magnetic field");
  }
  fParticles->Clear();  
  int npart = 0;
  //
  while (npart < mNPart) { // until needed numbe of muons generated
    int trials = 0;

    do { // until particle passes all selections
      if (++trials > mMaxTrials) {
        throw std::runtime_error("max. trials reached");
      }
      int pdg = gRandom->Rndm() < MuMinusFraction ? MuMinusPDG : MuPlusPDG; // mu- : mu+
      float r[3] = {0.f, mROrigin, 0.f}, p[3], ptot = 0, pt = 0;

      if (mParam == GenParamType::ParamMI) {
        ptot = mGenFun->GetRandom();
        p[1] = -ptot;
        if (gRandom->Rndm() > 0.9) {
          p[0] = gRandom->Gaus(0.0, 0.4) * ptot;
          p[2] = gRandom->Gaus(0.0, 0.4) * ptot;
        } else {
          p[0] = gRandom->Gaus(0.0, 0.2) * ptot;
          p[2] = gRandom->Gaus(0.0, 0.2) * ptot;
        }
        ptot = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        pt = std::sqrt(p[0] * p[0] + p[1] * p[1]);
      } else {
        ptot = mGenFun->GetRandom();
        float theta = 0, phi = 0;
        do {
          theta = gRandom->Gaus(0.5 * PiConst, 0.42);
        } while (std::abs(theta - 0.5 * PiConst) > mMaxAngleWRTVertical);
        do {
          phi = gRandom->Gaus(-0.5 * PiConst, 0.42);
        } while (std::abs(phi + 0.5 * PiConst) > mMaxAngleWRTVertical);

        pt = ptot * std::sin(theta);
        p[0] = pt * std::cos(phi);
        p[1] = pt * std::sin(phi);
        p[2] = ptot * std::cos(theta);
      }

      if (ptot < mPMin || ptot > mPMax || std::acos(std::abs(p[1]) / ptot) > mMaxAngleWRTVertical) {
        continue;
      }

      // estimate max deflection
      float xpos = 999, zpos = 999;
      if (!getXZatOrigin(xpos, zpos, r, p, -pdg)) {
        continue;
      }
      auto slpX = p[0] / p[1], slpZ = p[2] / p[1];
      xpos += mROrigin * slpX; // get max bending
      zpos += mROrigin * slpZ; // only, w/o slopes contribution
      float xmin = -mXAcc, xmax = mXAcc, zmin = -mZAcc, zmax = mZAcc;
      if (xpos < 0) {
        xmax -= xpos;
      } else {
        xmin -= xpos;
      }
      if (zpos < 0) {
        zmax -= zpos;
      } else {
        zmin -= zpos;
      }

      r[0] = mROrigin * slpX + xmin + gRandom->Rndm() * (xmax - xmin);
      r[2] = mROrigin * slpZ + zmin + gRandom->Rndm() * (zmax - zmin);

      // propagate to fixed Y=mROrigin plane to fixed radius in field free region: solve quadratic equation of circle - line intersection
      auto a = slpX * slpX + 1, xred = r[0] - r[1] * slpX, b = xred * slpX, det = b * b - a * (xred * xred - mROrigin * mROrigin);
      if (det < 0.) {
        continue;
      }
      r[1] = (-b + std::sqrt(det)) / a;
      r[0] += (r[1] - mROrigin) * slpX;
      r[2] += (r[1] - mROrigin) * p[2] / p[1];

      // check trigger condition
      if (!getXZatOrigin(xpos, zpos, r, p, -pdg) || std::abs(xpos) > mXAcc || std::abs(zpos) > mZAcc) {
        continue;
      }

      auto etot = std::sqrt(MuMass * MuMass + ptot * ptot);
      fParticles->Add( new TParticle(pdg, 1, -1, -1, -1, -1, p[0], p[1], p[2], etot, r[0], r[1], r[2], 0) );
      break;
    } while (1);
    npart++;
  }
}

//-----------------------------------------------------------------------------
void GeneratorCosmics::Init()
{
  //
  // Initialisation, check consistency of selected ranges
  //
  switch (mParam) {
    case GenParamType::ParamMI:
      mGenFun.reset(new TF1("genFun", "exp(-x/30.)", mPMin, mPMax * 2));
      break;
    case GenParamType::ParamACORDE:
      mGenFun.reset(new TF1("genFun", "x/(1.+(x/12.8)*(x/12.8))^1.96", mPMin, mPMax));
      break;
    case GenParamType::ParamTPC:
      mGenFun.reset(new TF1("genFun", "x/(1.+(x/3.)*(x/3.))^1.", mPMin, mPMax));
      break;
  }
  printf("Cosmics generator configuration:\n");
  printf("Parameterization type: %d with %e < p < %e\n", int(mParam), mPMin, mPMax);
  printf("Tracks created at R=%.2f and requested to have |X|<%.2f  and |Z|<%.2f at Y=0\n", mROrigin, mXAcc, mZAcc); 
  if (detectField()) {
    printf("Magnetic field %f\n", mBkG);
  }
  return;
}

void GeneratorCosmics::setMaxAngleWRTVertical(float max)
{
  if (max < 1e-6) {
    max = 1e-6;
  }
  if (max > 90.) {
    throw std::runtime_error("angle must be in ]0 : 90] range");
  }
  mMaxAngleWRTVertical = max * PiConst / 180.;
}

bool GeneratorCosmics::getXZatOrigin(float& xpos, float& zpos, const float r[3], const float p[3], int q) const
{
  // get position of the track at Y=0 in the lab frame
  constexpr float B2C = -0.299792458e-3, Almost0 = 1e-9, Almost1 = 1 - Almost0;
  // work explicitly in alpha = -pi/2 frame, i.e. track param X = -r[1]
  auto pt = std::sqrt(p[0] * p[0] + p[1] * p[1]), q2pt = q > 0 ? 1.f / pt : -1.f / pt;
  auto phi = std::atan2(p[0], -p[1]), snp = std::sin(phi), tgl = p[2] / pt;
  auto crv = q2pt * mBkG * B2C, x2r = crv * r[1], f1 = snp, f2 = f1 + x2r;
  if (std::abs(f1) >= Almost1 || std::abs(f2) >= Almost1 || std::abs(q2pt) < Almost0) {
    return false;
  }
  auto r1 = std::sqrt((1. - f1) * (1. + f1)), r2 = std::sqrt((1. - f2) * (1. + f2));
  if (std::abs(r1) < Almost0 || std::abs(r2) < Almost0) {
    return false;
  }
  auto dy2dx = (f1 + f2) / (r1 + r2);
  xpos = r[0] + r[1] * dy2dx; // print dy2dx, it is wrong?
  if (std::abs(x2r) < 0.05) {
    zpos = r[2] + r[1] * (r2 + f2 * dy2dx) * tgl;
  } else {
    auto rot = std::asin(r1 * f2 - r2 * f1);    // more economic version from Yura.
    if (f1 * f1 + f2 * f2 > 1 && f1 * f2 < 0) { // special cases of large rotations or large abs angles
      rot = f2 > 0 ? PiConst - rot : -PiConst - rot;
    }
    zpos = r[2] + tgl / crv * rot;
  }
  return true;
}

//______________________________________________________________________________
int GeneratorCosmics::ImportParticles(TClonesArray *particles, Option_t *option)
{
  if (particles == 0) return 0;
  TClonesArray &clonesParticles = *particles;
  clonesParticles.Clear();
  Int_t numpart = fParticles->GetEntries();
  for (int i = 0; i<numpart; i++) {
    TParticle *particle = (TParticle *)fParticles->At(i);
    new(clonesParticles[i]) TParticle(*particle);
  }
  return numpart;
}

//______________________________________________________________________________
void GeneratorCosmics::setBkG(float b)
{
  mBkG = b;
  mFieldIsSet = true;
  printf("Setting field to %f kG\n", mBkG);
  return;
}
