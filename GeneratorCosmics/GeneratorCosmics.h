// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_GENCOSMICS_H
#define O2_GENCOSMICS_H

// Generator for muons according to kinematic parametrizations at ALICE
// (not at the surface).
// Origin: andrea.dainese@lnl.infn.it
// Modifications for O2: ruben.shahoyan@cern.ch

#include "TGenerator.h"
#include "TClonesArray.h"
#include <TF1.h>

// Generates requested number of cosmic muons per call, requiring them to pass through
// certain |X|, |Z| at Y=0. The muons are generated on the surface of cylinder of the radius mROrigin

class GeneratorCosmics : public TGenerator
{
 public:
  static constexpr float PiConst = 3.1415927;
  static constexpr float MuMinusFraction = 4. / 9.;
  enum class GenParamType : int { ParamMI, ParamACORDE, ParamTPC }; // source parameterizations

  GeneratorCosmics();
  virtual ~GeneratorCosmics() {}

  virtual void GenerateEvent();
  virtual void Init();
  virtual int ImportParticles(TClonesArray *particles, Option_t *option);
  
  void setParam(GenParamType p) { mParam = p; }
  GenParamType getParam() const { return mParam; }
  void setParamMI() { setParam(GenParamType::ParamMI); }
  void setParamACORDE() { setParam(GenParamType::ParamACORDE); }
  void setParamTPC() { setParam(GenParamType::ParamTPC); }

  void setNPart(int n) { mNPart = n > 1 ? n : 1; }
  int getNPart() const { return mNPart; }

  void setMaxTrials(int n) { mMaxTrials = n < 1 ? 1 : n; }
  int getMaxTrials() const { return mMaxTrials; }

  void setROrigin(float r = 550.)
  {
    mROrigin = r;
    return;
  }
  void setMaxAngleWRTVertical(float max = 45.);

  void setBkG(float b);

  void setPRange(float pmin, float pmax)
  {
    mPMin = pmin < 0.5 ? 0.5 : pmin;
    mPMax = pmax < mPMin + 1e-3 ? mPMin + 1e-3 : pmax;
    if (mGenFun) {
      mGenFun->SetRange(pmin, pmax);
    }
  }

  void requireXZAccepted(float x, float z)
  {
    mXAcc = x < 1 ? 1. : x;
    mZAcc = z < 1 ? 1. : z;
  }
  void requireITS0() { requireXZAccepted(2.7, 27.1); }
  void requireITS1() { requireXZAccepted(3.5, 27.1); }
  void requireITS2() { requireXZAccepted(4.3, 27.1); }
  void requireITS3() { requireXZAccepted(19.8, 84.3); }
  void requireITS4() { requireXZAccepted(24.8, 84.3); }
  void requireITS5() { requireXZAccepted(34.6, 147.5); }
  void requireITS6() { requireXZAccepted(39.5, 147.5); }
  void requireTPC() { requireXZAccepted(250, 250); }

  bool getXZatOrigin(float& xpos, float& zpos, const float r[3], const float p[3], int q) const;

 private:
  GenParamType mParam = GenParamType::ParamTPC;
  std::unique_ptr<TF1> mGenFun;

  int mNPart = 1;                           // number of particle per event
  int mMaxTrials = 10000000;                // max trials to generat single muon
  float mROrigin = 550.;                    // R of muon origin
  float mMaxAngleWRTVertical = PiConst / 4; // maximum angle between momentum and y axis
  float mBkG = 0.;                          // field in kGauss

  float mPMin = 0.5;
  float mPMax = 50.;

  float mXAcc = 250.; // max |X| of track at Y = 0
  float mZAcc = 250.; // max |Z| of track at Y = 0

  bool mFieldIsSet = false;
  
  ClassDef(GeneratorCosmics, 1) // parametrized cosmics generator
};

#endif
