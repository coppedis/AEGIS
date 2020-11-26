// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_SLOWNUCLEONMODELEXP
#define O2_SLOWNUCLEONMODELEXP
//
// Experimental data inspired Gray Particle Model for p-Pb collisions
// Fluctuations are calculated from a binomial distribution.
// Author: A.Morsch
//

#include "SlowNucleonModel.h"

class SlowNucleonModelExp : public SlowNucleonModel {
public:
  SlowNucleonModelExp();
  virtual ~SlowNucleonModelExp() { ; }
  virtual void GetNumberOfSlowNucleons(Int_t ncoll, Int_t &ngp, Int_t &ngn,
                                       Int_t &nbp, Int_t &nbn) const;
  virtual void GetNumberOfSlowNucleons2(Int_t ncoll, Int_t &ngp, Int_t &ngn,
                                        Int_t &nbp, Int_t &nbn) const;
  virtual void GetNumberOfSlowNucleons2s(Int_t ncoll, Int_t &ngp, Int_t &ngn,
                                         Int_t &nbp, Int_t &nbn) const;
  // 1st model
  virtual void SetParameters(Float_t alpha1, Float_t alpha2);
  virtual void SetSaturation(Bool_t saturation) {
    fApplySaturation = saturation;
  }
  virtual void SetSaturationParams(Int_t ngray = 15, Int_t nblack = 28) {
    fnGraySaturation = ngray;
    fnBlackSaturation = nblack;
  }
  // 2nd model
  virtual void SetLCPparam(Float_t al) { fLCPparam = al; }
  virtual void SetNslowParams(Float_t a, Float_t b, Float_t c) {
    fSlownparam[0] = a;
    fSlownparam[1] = b;
    fSlownparam[2] = c;
  }

protected:
  Float_t fP;          // Number of protons  in the target
  Float_t fN;          // Number of neutrons in the target
  Float_t fAlphaGray;  // Proportionality between gray   particles and number of
                       // collisions
  Float_t fAlphaBlack; // Proportionality between black  particles and number of
                       // collisions
  Bool_t fApplySaturation; // If true apply satoration to N_black vs. N_gray
  Int_t fnGraySaturation;  // N_gray value for N_black saturation
  Int_t fnBlackSaturation; // N_black saturation value
  //
  // Adding parameters for 2nd model that can be tuned during config
  Float_t fLCPparam;      // parameter to calculate LCP from <Nslow p>
  Float_t fSlownparam[3]; // parameters to calculate <Nslow n> from LCP
  //
  // Adding parameter to smear the number of slow nucleons
  Float_t fSigmaSmear;

  ClassDef(SlowNucleonModelExp, 4) // Gray Particle Model (Experiment inspired)
};
#endif
