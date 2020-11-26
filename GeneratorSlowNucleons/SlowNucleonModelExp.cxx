// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//
// Experimental data inspired Gray Particle Model for p-Pb collisions
// The number of gray nucleons  is proportional to the number of collisions.
// The number of black nucleons is proportional to the number of collisions
// Fluctuations are calculated from a binomial distribution.
// Author: A.Morsch
//

#include "SlowNucleonModelExp.h"
#include <TMath.h>
#include <TRandom.h>

ClassImp(SlowNucleonModelExp)

    SlowNucleonModelExp::SlowNucleonModelExp()
    : fP(82), fN(126), fAlphaGray(2.3), fAlphaBlack(3.6),
      fApplySaturation(kTRUE), fnGraySaturation(15), fnBlackSaturation(28),
      fLCPparam(0.585), fSigmaSmear(0.25) {
  //
  // Default constructor
  //
  //
  fSlownparam[0] = 60.;
  fSlownparam[1] = 469.2;
  fSlownparam[2] = 8.762;
  /*printf("\n\n ******** Initializing slow nucleon model with parameters:\n");
  printf(" \t alpha_{gray} %1.2f  alpha_{black} %1.2f\n",fAlphaGray,
  fAlphaBlack); printf(" \t SATURATION %d w. %d (gray) %d (black)
  \n\n",fApplySaturation,fnGraySaturation,fnBlackSaturation); printf(" \t LCP
  parameter %f   Slown parameters = {%f, %f,
  %f}\n\n",fLCPparam,fSlownparam[0],fSlownparam[1],fSlownparam[2]); */
}

void SlowNucleonModelExp::GetNumberOfSlowNucleons(Int_t ncoll, Int_t &ngp,
                                                  Int_t &ngn, Int_t &nbp,
                                                  Int_t &nbn) const {
  //
  // Return the number of black and gray nucleons
  //
  // Number of collisions

  Float_t nu = (Float_t)(ncoll);

  // Mean number of gray nucleons

  Float_t nGray = fAlphaGray * nu;
  Float_t nGrayNeutrons = nGray * fN / (fN + fP);
  Float_t nGrayProtons = nGray - nGrayNeutrons;

  // Mean number of black nucleons
  Float_t nBlack = 0.;
  if (!fApplySaturation || (fApplySaturation && nGray < fnGraySaturation))
    nBlack = fAlphaBlack * nu;
  else if (fApplySaturation && nGray >= fnGraySaturation)
    nBlack = fnBlackSaturation;
  Float_t nBlackNeutrons = nBlack * 0.84;
  Float_t nBlackProtons = nBlack - nBlackNeutrons;

  // Actual number (including fluctuations) from binomial distribution
  Double_t p;

  //  gray neutrons
  p = nGrayNeutrons / fN;
  ngn = gRandom->Binomial((Int_t)fN, p);

  //  gray protons
  p = nGrayProtons / fP;
  ngp = gRandom->Binomial((Int_t)fP, p);

  //  black neutrons
  p = nBlackNeutrons / fN;
  nbn = gRandom->Binomial((Int_t)fN, p);

  //  black protons
  p = nBlackProtons / fP;
  nbp = gRandom->Binomial((Int_t)fP, p);
}

void SlowNucleonModelExp::GetNumberOfSlowNucleons2(Int_t ncoll, Int_t &ngp,
                                                   Int_t &ngn, Int_t &nbp,
                                                   Int_t &nbn) const {
  //
  // Return the number of black and gray nucleons
  //
  // Number of collisions

  // based on E910 model
  // ================================================================

  Float_t nu = (Float_t)(ncoll);
  //
  // nu = nu+1.*gRandom->Rndm();
  nu = gRandom->Gaus(nu, 0.5);
  if (nu < 0.)
    nu = 0.;
  //
  Float_t poverpd = 0.843;
  Float_t zAu2zPb = 82. / 79.;
  Float_t nGrayp = (-0.27 + 0.63 * nu - 0.0008 * nu * nu) * poverpd * zAu2zPb;

  //  gray protons
  Double_t p;
  p = nGrayp / fP;
  ngp = gRandom->Binomial((Int_t)fP, p);
  // ngp = gRandom->Gaus(nGrayp, TMath::Sqrt(fP*p*(1-p)));
  if (nGrayp < 0.)
    ngp = 0;

  // Float_t blackovergray = 3./7.;// from spallation
  Float_t blackovergray = 0.65; // from COSY
  Float_t nBlackp = blackovergray * nGrayp;

  //  black protons
  p = nBlackp / fP;
  nbp = gRandom->Binomial((Int_t)fP, p);
  // nbp = gRandom->Gaus(nBlackp, TMath::Sqrt(fP*p*(1-p)));
  if (nBlackp < 0.)
    nbp = 0;

  if (nu < 3.) {
    nGrayp = -0.836 + 0.9112 * nu - 0.05381 * nu * nu;
    nBlackp = blackovergray * nGrayp;
  }

  // printf(" \t Using LCP parameter %f   Slown parameters = {%f, %f,
  // %f}\n\n",fLCPparam,fSlownparam[0],fSlownparam[1],fSlownparam[2]);
  Float_t nGrayNeutrons = 0.;
  Float_t nBlackNeutrons = 0.;
  Float_t cp = (nGrayp + nBlackp) / fLCPparam;

  if (cp > 0.) {
    Float_t nSlow = fSlownparam[0] + fSlownparam[1] / (-fSlownparam[2] - cp);
    Float_t paramRetta =
        fSlownparam[0] + fSlownparam[1] / (-fSlownparam[2] - 3);
    if (cp < 3.)
      nSlow = 0. + (paramRetta - 0.) / (3. - 0.) * (cp - 0.);

    nGrayNeutrons = nSlow * 0.1;
    nBlackNeutrons = nSlow - nGrayNeutrons;
  } else {
    // Sikler "pasturato" (qui non entra mai!!!!)
    nGrayNeutrons = 0.47 * fAlphaGray * nu;
    nBlackNeutrons = 0.88 * fAlphaBlack * nu;
    // printf("nslowp=0 -> ncoll = %1.0f -> ngrayn = %1.0f  nblackn = %1.0f \n",
    // nu, nGrayNeutrons, nBlackNeutrons);
  }

  //  gray neutrons
  p = nGrayNeutrons / fN;
  //    ngn = gRandom->Binomial((Int_t) fN, p);
  ngn = gRandom->Gaus(nGrayNeutrons, TMath::Sqrt(fN * p * (1 - p)));

  //  black neutrons
  p = nBlackNeutrons / fN;
  //    nbn = gRandom->Binomial((Int_t) fN, p);
  nbn = gRandom->Gaus(nBlackNeutrons, TMath::Sqrt(fN * p * (1 - p)));
}

void SlowNucleonModelExp::GetNumberOfSlowNucleons2s(Int_t ncoll, Int_t &ngp,
                                                    Int_t &ngn, Int_t &nbp,
                                                    Int_t &nbn) const {
  //
  // Return the number of black and gray nucleons
  //
  // Number of collisions

  // based on E910 model
  // ================================================================

  Float_t nu = (Float_t)(ncoll);
  //
  Float_t poverpd = 0.843;
  Float_t zAu2zPb = 82. / 79.;
  Float_t grayp = (-0.27 + 0.63 * nu - 0.0008 * nu * nu) * poverpd * zAu2zPb;
  Float_t nGrayp = gRandom->Gaus(grayp, fSigmaSmear);
  if (nGrayp < 0.)
    nGrayp = 0.;

  //  gray protons
  Double_t p = 0.;
  p = nGrayp / fP;
  ngp = gRandom->Binomial((Int_t)fP, p);
  // ngp = gRandom->Gaus(nGrayp, TMath::Sqrt(fP*p*(1-p)));
  if (nGrayp < 0.)
    ngp = 0;

  // Float_t blackovergray = 3./7.;// from spallation
  Float_t blackovergray = 0.65; // from COSY
  // Float_t blackp  = blackovergray*grayp;
  // Float_t nBlackp = gRandom->Gaus(nblackp, fSigmaSmear);
  Float_t nBlackp = blackovergray * nGrayp;
  if (nBlackp < 0.)
    nBlackp = 0.;

  //  black protons
  p = nBlackp / fP;
  nbp = gRandom->Binomial((Int_t)fP, p);
  // nbp = gRandom->Gaus(nBlackp, TMath::Sqrt(fP*p*(1-p)));
  if (nBlackp < 0.)
    nbp = 0;

  Float_t nGrayNeutrons = 0.;
  Float_t nBlackNeutrons = 0.;
  Float_t cp = (nGrayp + nBlackp) / fLCPparam;

  if (cp > 0.) {
    Float_t nSlow = fSlownparam[0] + fSlownparam[1] / (-fSlownparam[2] - cp);

    nGrayNeutrons = nSlow * 0.1;
    nBlackNeutrons = nSlow - nGrayNeutrons;
  } else {
    // Sikler "pasturato" (qui non entra mai!!!!)
    nGrayNeutrons = 0.47 * fAlphaGray * nu;
    nBlackNeutrons = 0.88 * fAlphaBlack * nu;
    // printf("nslowp=0 -> ncoll = %1.0f -> ngrayn = %1.0f  nblackn = %1.0f \n",
    // nu, nGrayNeutrons, nBlackNeutrons);
  }
  //
  if (nGrayNeutrons < 0.)
    nGrayNeutrons = 0.;
  if (nBlackNeutrons < 0.)
    nBlackNeutrons = 0.;

  //  gray neutrons
  p = nGrayNeutrons / fN;
  //    ngn = gRandom->Binomial((Int_t) fN, p);
  ngn = gRandom->Gaus(nGrayNeutrons, TMath::Sqrt(fN * p * (1 - p)));
  if (nGrayNeutrons < 0.)
    ngn = 0;

  //  black neutrons
  p = nBlackNeutrons / fN;
  //    nbn = gRandom->Binomial((Int_t) fN, p);
  nbn = gRandom->Gaus(nBlackNeutrons, TMath::Sqrt(fN * p * (1 - p)));
  if (nBlackNeutrons < 0.)
    nbn = 0;
}

void SlowNucleonModelExp::SetParameters(Float_t alpha1, Float_t alpha2) {
  // Set the model parameters
  fAlphaGray = alpha1;
  fAlphaBlack = alpha2;
}
