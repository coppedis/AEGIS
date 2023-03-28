// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Class to generate particles from using paramtrized pT and y distributions.
// Distributions are obtained from pointer to object of type
// GeneratorParamLibBase. Decays are performed using Pythia.
// andreas.morsch@cern.ch

// Helper class for the configuration of the Pythia6 decayer
// Allows to force decay channels.
// Author: andreas.morsch@cern.ch

#include "PythiaDecayerConfig.h"
#include <TClonesArray.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TPythia6.h>
#include <TRandom.h>

ClassImp(PythiaDecayerConfig)

    PythiaDecayerConfig::PythiaDecayerConfig()
    : TVirtualMCDecayer(), fDecay(kAll), fHeavyFlavour(kTRUE), fLongLived(kFALSE),
      fPatchOmegaDalitz(0), fDecayerExodus(nullptr), fPi0(1) {
  // Default Constructor
  for (Int_t i = 0; i < 501; i++)
    fBraPart[i] = 1.;
}

PythiaDecayerConfig::PythiaDecayerConfig(const PythiaDecayerConfig &decayer)
    : fDecay(kAll), fHeavyFlavour(kTRUE), fLongLived(kFALSE),
      fPatchOmegaDalitz(0), fDecayerExodus(nullptr), fPi0(1) {
  // Copy Constructor
  decayer.Copy(*this);
  for (Int_t i = 0; i < 501; i++)
    fBraPart[i] = 0.;
}

void PythiaDecayerConfig::Init() {

  if (fDecayerExodus){
    fDecayerExodus->Init();
  }

  // Switch on heavy flavor decays
  fPythia = TPythia6::Instance();
  Int_t kc, i, j;
  Int_t heavy[14] = {411, 421, 431, 4122, 4132, 4232, 4332,
                     511, 521, 531, 5122, 5132, 5232, 5332};
  for (j = 0; j < 14; j++) {
    kc = fPythia->Pycomp(heavy[j]);
    if (fDecay == kNoDecayHeavy) {
      fPythia->SetMDCY(kc, 1, 0);
    } else {
      fPythia->SetMDCY(kc, 1, 1);
      for (i = fPythia->GetMDCY(kc, 2);
           i < fPythia->GetMDCY(kc, 2) + fPythia->GetMDCY(kc, 3); i++) {
        fPythia->SetMDME(i, 1, 1);
      }
    }
  }

  if (fPi0)
    fPythia->SetMDCY(fPythia->Pycomp(111), 1, 1);

  Int_t isw = 0;
  if (fLongLived)
    isw = 1;

  fPythia->SetMDCY(fPythia->Pycomp(130), 1, isw);
  fPythia->SetMDCY(fPythia->Pycomp(310), 1, isw);
  fPythia->SetMDCY(fPythia->Pycomp(3122), 1, isw);
  fPythia->SetMDCY(fPythia->Pycomp(3112), 1, isw);
  fPythia->SetMDCY(fPythia->Pycomp(3222), 1, isw);
  fPythia->SetMDCY(fPythia->Pycomp(3312), 1, isw);
  fPythia->SetMDCY(fPythia->Pycomp(3322), 1, isw);
  fPythia->SetMDCY(fPythia->Pycomp(3334), 1, isw);

  // .. Force decay channels
  ForceDecay();
}

void PythiaDecayerConfig::SwitchOffParticle(Int_t kf) {
  // switch off decay for particle "kf"
  fPythia->SetMDCY(fPythia->Pycomp(kf), 1, 0);
}

void PythiaDecayerConfig::ForceDecay() {
  // Force a particle decay mode
  // Switch heavy flavour production off if requested
  if (!fHeavyFlavour)
    SwitchOffHeavyFlavour();
  //
  Decay_t decay = fDecay;
  fPythia->SetMSTJ(21, 2);

  //
  // select mode
  Int_t products[2];
  Int_t mult[2];
  Int_t products1[3];
  Int_t mult1[3];

  switch (decay) {
  case kHardMuons:
    products1[0] = 13;
    products1[1] = 443;
    products1[2] = 100443;
    mult1[0] = 1;
    mult1[1] = 1;
    mult1[2] = 1;
    ForceParticleDecay(511, products1, mult1, 3);
    ForceParticleDecay(521, products1, mult1, 3);
    ForceParticleDecay(531, products1, mult1, 3);
    ForceParticleDecay(5122, products1, mult1, 3);
    ForceParticleDecay(5132, products1, mult1, 3);
    ForceParticleDecay(5232, products1, mult1, 3);
    ForceParticleDecay(5332, products1, mult1, 3);
    ForceParticleDecay(100443, 443, 1); // Psi'  -> J/Psi X
    ForceParticleDecay(443, 13, 2);     // J/Psi -> mu+ mu-
    ForceParticleDecay(411, 13, 1);     // D+/-
    ForceParticleDecay(421, 13, 1);     // D0
    ForceParticleDecay(431, 13, 1);     // D_s
    ForceParticleDecay(4122, 13, 1);    // Lambda_c
    ForceParticleDecay(4132, 13, 1);    // Xsi_c
    ForceParticleDecay(4232, 13, 1);    // Sigma_c
    ForceParticleDecay(4332, 13, 1);    // Omega_c
    break;
  case kChiToJpsiGammaToMuonMuon:
    products[0] = 443;
    products[1] = 22;
    mult[0] = 1;
    mult[1] = 1;
    ForceParticleDecay(20443, products, mult, 2); // Chi_1c  -> J/Psi  Gamma
    ForceParticleDecay(445, products, mult, 2);   // Chi_2c  -> J/Psi  Gamma
    ForceParticleDecay(443, 13, 2);               // J/Psi -> mu+ mu-
    break;
  case kChiToJpsiGammaToElectronElectron:
    products[0] = 443;
    products[1] = 22;
    mult[0] = 1;
    mult[1] = 1;
    ForceParticleDecay(20443, products, mult, 2); // Chi_1c  -> J/Psi  Gamma
    ForceParticleDecay(445, products, mult, 2);   // Chi_2c  -> J/Psi  Gamma
    ForceParticleDecay(443, 11, 2);               // J/Psi -> e+ e-
    break;

  case kBSemiMuonic:
    ForceParticleDecay(511, 13, 1);  // B0
    ForceParticleDecay(521, 13, 1);  // B+/-
    ForceParticleDecay(531, 13, 1);  // B_s
    ForceParticleDecay(5122, 13, 1); // Lambda_b
    ForceParticleDecay(5132, 13, 1); // Xsi_b
    ForceParticleDecay(5232, 13, 1); // Sigma_b
    ForceParticleDecay(5332, 13, 1); // Omega_b
    break;
  case kDSemiMuonic:
    ForceParticleDecay(411, 13, 1);  // D0
    ForceParticleDecay(421, 13, 1);  // D+/-
    ForceParticleDecay(431, 13, 1);  // D_s
    ForceParticleDecay(4122, 13, 1); // Lambda_c
    ForceParticleDecay(4132, 13, 1); // Xsi_c
    ForceParticleDecay(4232, 13, 1); // Sigma_c
    ForceParticleDecay(4332, 13, 1); // Omega_c
    break;
  case kSemiMuonic:
    ForceParticleDecay(411, 13, 1);  // D+/-
    ForceParticleDecay(421, 13, 1);  // D0
    ForceParticleDecay(431, 13, 1);  // D_s
    ForceParticleDecay(4122, 13, 1); // Lambda_c
    ForceParticleDecay(4132, 13, 1); // Xsi_c
    ForceParticleDecay(4232, 13, 1); // Sigma_c
    ForceParticleDecay(4332, 13, 1); // Omega_c
    ForceParticleDecay(511, 13, 1);  // B0
    ForceParticleDecay(521, 13, 1);  // B+/-
    ForceParticleDecay(531, 13, 1);  // B_s
    ForceParticleDecay(5122, 13, 1); // Lambda_b
    ForceParticleDecay(5132, 13, 1); // Xsi_b
    ForceParticleDecay(5232, 13, 1); // Sigma_b
    ForceParticleDecay(5332, 13, 1); // Omega_b
    break;
  case kDiMuon:
    ForceParticleDecay(113, 13, 2);    // rho
    ForceParticleDecay(221, 13, 2);    // eta
    ForceParticleDecay(223, 13, 2);    // omega
    ForceParticleDecay(333, 13, 2);    // phi
    ForceParticleDecay(443, 13, 2);    // J/Psi
    ForceParticleDecay(100443, 13, 2); // Psi'
    ForceParticleDecay(553, 13, 2);    // Upsilon
    ForceParticleDecay(100553, 13, 2); // Upsilon'
    ForceParticleDecay(200553, 13, 2); // Upsilon''
    break;
  case kJpsiDiMuon:
    ForceParticleDecay(443, 13, 2); // J/Psi
    break;
  case kBSemiElectronic:
    ForceParticleDecay(511, 11, 1);  // B0
    ForceParticleDecay(521, 11, 1);  // B+/-
    ForceParticleDecay(531, 11, 1);  // B_s
    ForceParticleDecay(5122, 11, 1); // Lambda_b
    ForceParticleDecay(5132, 11, 1); // Xsi_b
    ForceParticleDecay(5232, 11, 1); // Sigma_b
    ForceParticleDecay(5332, 11, 1); // Omega_b
    break;
  case kSemiElectronic:
    ForceParticleDecay(411, 11, 1);  // D+/-
    ForceParticleDecay(421, 11, 1);  // D0
    ForceParticleDecay(431, 11, 1);  // D_s
    ForceParticleDecay(4122, 11, 1); // Lambda_c
    ForceParticleDecay(4132, 11, 1); // Xsi_c
    ForceParticleDecay(4232, 11, 1); // Sigma_c
    ForceParticleDecay(4332, 11, 1); // Omega_c
    ForceParticleDecay(511, 11, 1);  // B0
    ForceParticleDecay(521, 11, 1);  // B+/-
    ForceParticleDecay(531, 11, 1);  // B_s
    ForceParticleDecay(5122, 11, 1); // Lambda_b
    ForceParticleDecay(5132, 11, 1); // Xsi_b
    ForceParticleDecay(5232, 11, 1); // Sigma_b
    ForceParticleDecay(5332, 11, 1); // Omega_b
    break;
  case kDiElectron:
    ForceParticleDecay(113, 11, 2);    // rho
    ForceParticleDecay(333, 11, 2);    // phi
    ForceParticleDecay(221, 11, 2);    // eta
    ForceParticleDecay(223, 11, 2);    // omega
    ForceParticleDecay(443, 11, 2);    // J/Psi
    ForceParticleDecay(100443, 11, 2); // Psi'
    ForceParticleDecay(553, 11, 2);    // Upsilon
    ForceParticleDecay(100553, 11, 2); // Upsilon'
    ForceParticleDecay(200553, 11, 2); // Upsilon''
    break;
  case kPsiPrimeJpsiDiElectron:
    products[0] = 443;
    products[1] = 211;
    mult[0] = 1;
    mult[1] = 2;
    ForceParticleDecay(100443, products, mult, 2, 1);
    ForceParticleDecay(443, 11, 2);
    break;
  case kBJpsiDiMuon:
    products[0] = 443;
    products[1] = 100443;
    mult[0] = 1;
    mult[1] = 1;
    ForceParticleDecay(511, products, mult, 2);  // B0   -> J/Psi (Psi') X
    ForceParticleDecay(521, products, mult, 2);  // B+/- -> J/Psi (Psi') X
    ForceParticleDecay(531, products, mult, 2);  // B_s  -> J/Psi (Psi') X
    ForceParticleDecay(5122, products, mult, 2); // Lambda_b -> J/Psi (Psi') X
    ForceParticleDecay(100443, 443, 1);          // Psi'  -> J/Psi X
    ForceParticleDecay(443, 13, 2);              // J/Psi -> mu+ mu-
    break;
  case kBPsiPrimeDiMuon:
    ForceParticleDecay(511, 100443, 1);  // B0
    ForceParticleDecay(521, 100443, 1);  // B+/-
    ForceParticleDecay(531, 100443, 1);  // B_s
    ForceParticleDecay(5122, 100443, 1); // Lambda_b
    ForceParticleDecay(100443, 13, 2);   // Psi'
    break;
  case kBJpsiDiElectron:
    ForceParticleDecay(511, 443, 1);  // B0
    ForceParticleDecay(521, 443, 1);  // B+/-
    ForceParticleDecay(531, 443, 1);  // B_s
    ForceParticleDecay(5122, 443, 1); // Lambda_b
    ForceParticleDecay(443, 11, 2);   // J/Psi
    break;
  case kBJpsi:
    ForceParticleDecay(511, 443, 1);  // B0
    ForceParticleDecay(521, 443, 1);  // B+/-
    ForceParticleDecay(531, 443, 1);  // B_s
    ForceParticleDecay(5122, 443, 1); // Lambda_b
    break;
  case kBJpsiUndecayed:
    ForceParticleDecay(511, 443, 1);              // B0
    ForceParticleDecay(521, 443, 1);              // B+/-
    ForceParticleDecay(531, 443, 1);              // B_s
    ForceParticleDecay(5122, 443, 1);             // Lambda_b
    fPythia->SetMDCY(fPythia->Pycomp(443), 1, 0); // switch-off J/psi
    break;
  case kBPsiPrimeDiElectron:
    ForceParticleDecay(511, 100443, 1);  // B0
    ForceParticleDecay(521, 100443, 1);  // B+/-
    ForceParticleDecay(531, 100443, 1);  // B_s
    ForceParticleDecay(5122, 100443, 1); // Lambda_b
    ForceParticleDecay(100443, 11, 2);   // Psi'
    break;
  case kPiToMu:
    ForceParticleDecay(211, 13, 1); // pi->mu
    break;
  case kKaToMu:
    ForceParticleDecay(321, 13, 1); // K->mu
    break;
  case kAllMuonic:
    ForceParticleDecay(211, 13, 1); // pi->mu
    ForceParticleDecay(321, 13, 1); // K->mu
    break;
  case kWToMuon:
    ForceParticleDecay(24, 13, 1); // W -> mu
    break;
  case kWToCharm:
    ForceParticleDecay(24, 4, 1); // W -> c
    break;
  case kWToCharmToMuon:
    ForceParticleDecay(24, 4, 1);    // W -> c
    ForceParticleDecay(411, 13, 1);  // D+/- -> mu
    ForceParticleDecay(421, 13, 1);  // D0  -> mu
    ForceParticleDecay(431, 13, 1);  // D_s  -> mu
    ForceParticleDecay(4122, 13, 1); // Lambda_c
    ForceParticleDecay(4132, 13, 1); // Xsi_c
    ForceParticleDecay(4232, 13, 1); // Sigma_c
    ForceParticleDecay(4332, 13, 1); // Omega_c
    break;
  case kZDiMuon:
    ForceParticleDecay(23, 13, 2); // Z -> mu+ mu-
    break;
  case kZDiElectron:
    ForceParticleDecay(23, 11, 2); // Z -> e+ e-
    break;
  case kHadronicD:
    ForceHadronicD(1, 0, 0);
    break;
  case kHadronicDWithV0:
    ForceHadronicD(1, 1, 0);
    break;
  case kHadronicDWithout4Bodies:
    ForceHadronicD(0, 0, 0);
    break;
  case kHadronicDWithout4BodiesWithV0:
    ForceHadronicD(0, 1, 0);
    break;
  case kHadronicDPionicD0:
    ForceHadronicD(1, 0, 3, 1);
    break;
  case kHadronicDWithV0PionicD0:
    ForceHadronicD(1, 1, 3, 1);
    break;
  case kHadronicDWithout4BodiesPionicD0:
    ForceHadronicD(0, 0, 3, 1);
    break;
  case kHadronicDWithout4BodiesWithV0PionicD0:
    ForceHadronicD(0, 1, 3, 1);
    break;
  case kHadronicDPionicD0pure:
    ForceHadronicD(0, 0, 3, 2);
    break;
  case kHadronicDPionicD0K:
    ForceHadronicD(0, 0, 3, 3);
    break;
  case kHadronicDPionicD0pi:
    ForceHadronicD(0, 0, 3, 4);
    break;
  case kPhiKK:
    ForceParticleDecay(333, 321, 2); // Phi->K+K-
    break;
  case kOmega:
    ForceOmega();
    break;
  case kLambda:
    ForceLambda();
    break;
  case kAll:
    break;
  case kNoDecay:
    fPythia->SetMSTJ(21, 0);
    break;
  case kNoDecayHeavy:
  case kNeutralPion:
    break;
  case kNoDecayBeauty:
    SwitchOffBDecay();
    break;
  case kElectronEM:
    ForceParticleDecay(111, 11, 1); // pi^0
    ForceParticleDecay(221, 11, 1); // eta
    ForceParticleDecay(113, 11, 1); // rho
    ForceParticleDecay(223, 11, 1); // omega
    ForceParticleDecay(331, 11, 1); // etaprime
    ForceParticleDecay(333, 11, 1); // phi
    ForceParticleDecay(443, 11, 1); // jpsi
    break;
  case kDiElectronEM:
    ForceParticleDecay(111, 11, 2); // pi^0
    ForceParticleDecay(221, 11, 2); // eta
    ForceParticleDecay(113, 11, 2); // rho
    ForceParticleDecay(223, 11, 2); // omega
    ForceParticleDecay(331, 11, 2); // etaprime
    ForceParticleDecay(333, 11, 2); // phi
    ForceParticleDecay(443, 11, 2); // jpsi
    break;
  case kGammaEM:
    ForceParticleDecay(111, 22, 1);  // pi^0
    ForceParticleDecay(221, 22, 1);  // eta
    ForceParticleDecay(113, 22, 1);  // rho
    ForceParticleDecay(223, 22, 1);  // omega
    ForceParticleDecay(331, 22, 1);  // etaprime
    ForceParticleDecay(333, 22, 1);  // phi
    ForceParticleDecay(443, 22, 1);  // jpsi
    ForceParticleDecay(3212, 22, 1); // Sigma 0
    ForceParticleDecay(3122, 22, 1); // Lambda
    ForceParticleDecay(310, 22, 1);  // K0s
    ForceParticleDecay(130, 22, 1);  // K0l
    ForceParticleDecay(2224, 22, 1); // Delta++
    ForceParticleDecay(2214, 22, 1); // Delta+
    ForceParticleDecay(1114, 22, 1); // Delta-
    ForceParticleDecay(2114, 22, 1); // Delta0
    ForceParticleDecay(213, 22, 1);  // rho+
    ForceParticleDecay(-213, 22, 1); // rho-
    ForceParticleDecay(313, 22, 1);  // K0star
    break;
  case kAllEM:
    products1[0] = 22;
    products1[1] = 11;
    products1[2] = 11;
    mult1[0] = 1;
    mult1[1] = 1;
    mult1[2] = 2;
    ForceParticleDecay(111, products1, mult1, 3); // pi0
    ForceParticleDecay(221, products1, mult1, 3); // eta
    ForceParticleDecay(113, products1, mult1, 3); // rho
    ForceParticleDecay(223, products1, mult1, 3); // omega
    ForceParticleDecay(331, products1, mult1, 3); // etaprime
    ForceParticleDecay(333, products1, mult1, 3); // phi
    ForceParticleDecay(443, products1, mult1, 3); // jpsi
    ForceParticleDecay(3212, 22, 1);              // Sigma0
    ForceParticleDecay(3122, 22, 1);              // Lambda
    ForceParticleDecay(310, 22, 1);               // K0s
    ForceParticleDecay(130, 22, 1);               // K0l
    ForceParticleDecay(2224, 22, 1);              // Delta++
    ForceParticleDecay(2214, 22, 1);              // Delta+
    ForceParticleDecay(1114, 22, 1);              // Delta-
    ForceParticleDecay(2114, 22, 1);              // Delta0
    ForceParticleDecay(213, 22, 1);               // rho+
    ForceParticleDecay(-213, 22, 1);              // rho-
    ForceParticleDecay(313, 22, 1);               // K0star
    break;
  case kBeautyUpgrade:
    ForceBeautyUpgrade();
    break;
  case kHFYellowReport:
    ForceHFYellowReport();
    break;
  case kLcpKpi:
    ForceHadronicD(0, 0, 1);
    break;
  case kLcpK0S:
    ForceHadronicD(0, 0, 2);
    break;
  case kLcpK0SBDTsig:
    ForceHadronicD(0, 0, 4);
    break;
  case kEtaPrime:
    if (gRandom->Rndm() < 0.5) {
      products1[0] = 211;
      products1[1] = -211;
      products1[2] = 221;
      mult1[0] = 1;
      mult1[1] = 1;
      mult1[2] = 1;
      ForceParticleDecay(331, products1, mult1, 1);
    } else {
      ForceParticleDecay(331, 22, 2);
    }
    break;
  default:
    break;
  }
}

void PythiaDecayerConfig::SwitchOffHeavyFlavour() {
  // Switch off heavy flavour production
  //
  // Maximum number of quark flavours used in pdf
  fPythia->SetMSTP(58, 3);
  // Maximum number of flavors that can be used in showers
  fPythia->SetMSTJ(45, 3);
  // Switch off g->QQbar splitting in decay table
  for (Int_t i = 156; i <= 160; i++)
    fPythia->SetMDME(i, 1, 0);
}

void PythiaDecayerConfig::ForceHFYellowReport() {
  //
  // Force dedicated decay channels of signals interesting
  // for the ITS upgrade and specifically the yellow report
  //

  // Lb: 100% of them in Lc  in final state

  ForceParticleDecay(5122, 4122, 1);

  // B0 -> D-e+nu
  const Int_t prod[3] = {411, 11, 12};
  Int_t mult[3] = {1, 1, 1};

  ForceParticleDecay(511, prod, mult, 3, 1);
  ForceHadronicD(0, 0, 0);
}

void PythiaDecayerConfig::ForceBeautyUpgrade() {
  //
  // Force dedicated decay channels of signals ineresting
  // for the ITS upgrade (Lb, Lc, Xi_c, B)
  //

  // Lb: 50% of them in Lc pi+ and the rest with a Lc in final state
  if (gRandom->Rndm() < 0.50) {
    const Int_t prod3[3] = {4122, 211, 0};
    Int_t mult3[3] = {1, 1, 1};
    ForceParticleDecay(5122, prod3, mult3, 3, 1);
  } else {
    ForceParticleDecay(5122, 4122, 1);
  }
  // B+ -> D0pi+
  const Int_t prod[2] = {421, 211};
  Int_t mult[2] = {1, 1};
  ForceParticleDecay(521, prod, mult, 2, 1);
  // B0 -> D*-pi+
  const Int_t prod2[2] = {413, 211};
  ForceParticleDecay(511, prod2, mult, 2, 1);
  // Bs -> Ds-pi+
  const Int_t prod3[2] = {431, 211};
  ForceParticleDecay(531, prod3, mult, 2, 1);
  // force charm hadronic decays (D mesons and Lc)
  ForceHadronicD(0, 0, 0);
}

Int_t PythiaDecayerConfig::CountProducts(Int_t channel, Int_t particle) {
  // Count number of decay products
  Int_t np = 0;
  for (Int_t i = 1; i <= 5; i++) {
    if (TMath::Abs(fPythia->GetKFDP(channel, i)) == particle)
      np++;
  }
  return np;
}

void PythiaDecayerConfig::ForceHadronicD(Int_t optUse4Bodies, Int_t optUseDtoV0,
                                         Int_t optForceLcChannel,
                                         Int_t optUsePionicD0) {

  // OmegaC -> Omega pi
  Int_t iOmegaC = 4332;
  Int_t productsOc[2] = {211, 3334}, multOc[2] = {1, 1};
  ForceParticleDecay(iOmegaC, productsOc, multOc, 2, 1);

  // Force golden D decay modes
  //
  const Int_t kNHadrons = 7;
  Int_t channel;
  Int_t hadron[kNHadrons] = {411, 421, 431, 4112, 4122, 4232, 4132};

  // for D+ -> K0* (-> K- pi+) pi+
  Int_t iKstar0 = 313;
  Int_t iKstarbar0 = -313;
  Int_t products[2] = {kKPlus, kPiMinus}, mult[2] = {1, 1};
  ForceParticleDecay(iKstar0, products, mult, 2);
  // for Ds -> Phi pi+
  Int_t iPhi = 333;
  ForceParticleDecay(iPhi, kKPlus, 2); // Phi->K+K-
  // for D0 -> rho0 pi+ k-
  Int_t iRho0 = 113;
  ForceParticleDecay(iRho0, kPiPlus, 2); // Rho0->pi+pi-
  // for Lambda_c -> Delta++ K-
  Int_t iDeltaPP = 2224;
  Int_t productsD[2] = {kProton, kPiPlus}, multD[2] = {1, 1};
  ForceParticleDecay(iDeltaPP, productsD, multD, 2);
  // for Lambda_c -> Lambda(1520) pi+ -> p K- pi+
  Int_t iLambda1520 = 3124;
  Int_t productsL[2] = {kProton, kKMinus}, multL[2] = {1, 1};
  ForceParticleDecay(iLambda1520, productsL, multL, 2);
  // for Lambda_c -> Lambda pi+
  Int_t iLambda = 3122;
  // for Lambda_c -> antiK0 p
  Int_t iK0bar = -311;
  // for Xic+->Xi-pi+pi+, Xic+->Xi0*pi+,Xi0*->Xi-pi+ and Xic0->Xi-pi+ channels
  Int_t iXiMinus = 3312, iXiStar0 = 3324;

  Int_t decayP1[kNHadrons][4] = {
      {kKMinus, kPiPlus, kPiPlus, 0}, {kKMinus, kPiPlus, 0, 0},
      {kKPlus, iKstarbar0, 0, 0},     {-1, -1, -1, -1},
      {kProton, iKstarbar0, 0, 0},    {kProton, iKstarbar0, 0, 0},
      {kPiPlus, iXiMinus, 0, 0}};
  Int_t decayP2[kNHadrons][4] = {{iKstarbar0, kPiPlus, 0, 0},
                                 {kKMinus, kPiPlus, kPiPlus, kPiMinus},
                                 {iPhi, kPiPlus, 0, 0},
                                 {-1, -1, -1, -1},
                                 {iDeltaPP, kKMinus, 0, 0},
                                 {kProton, kKMinus, kPiPlus, 0},
                                 {-1, -1, -1, -1}};
  Int_t decayP3[kNHadrons][4] = {{kPiPlus, iPhi, 0, 0},
                                 {kKMinus, kPiPlus, iRho0, 0},
                                 {kKPlus, iK0bar, 0, 0},
                                 {-1, -1, -1, -1},
                                 {kProton, kKMinus, kPiPlus, 0},
                                 {kPiPlus, kPiPlus, iXiMinus, 0},
                                 {-1, -1, -1, -1}};
  // for Lambda_c -> Lambda_1520 pi+ -> p K- pi+, D0-> K*0 pi+ pi- -> K3pi
  Int_t decayP4[kNHadrons][4] = {{iKstarbar0, kKPlus, 0, 0},
                                 {iKstarbar0, kPiPlus, kPiMinus, 0},
                                 {-1, -1, -1, -1},
                                 {-1, -1, -1, -1},
                                 {iLambda1520, kPiPlus, 0, 0},
                                 {kPiPlus, iXiStar0, 0, 0},
                                 {-1, -1, -1, -1}};
  // for Lambda_c -> Lambda pi+, D0 -> pi0 pi+ pi-
  Int_t decayP5[kNHadrons][4] = {
      {iK0bar, kPiPlus, 0, 0}, {kPiPlus, kPiMinus, kPi0, 0}, {-1, -1, -1, -1},
      {-1, -1, -1, -1},        {iLambda, kPiPlus, 0, 0},     {-1, -1, -1, -1},
      {-1, -1, -1, -1}};

  // for Lambda_c -> K0bar p, D0 -> K- pi+ pi0, D+ -> K0s pi+ pi0
  Int_t decayP6[kNHadrons][4] = {{iK0bar, kPiPlus, kPi0, 0},
                                 {kKMinus, kPiPlus, kPi0, 0},
                                 {-1, -1, -1, -1},
                                 {-1, -1, -1, -1},
                                 {kProton, iK0bar, 0, 0},
                                 {-1, -1, -1, -1},
                                 {-1, -1, -1, -1}

  };

  if (optUse4Bodies == 0) {
    for (Int_t iDau = 0; iDau < 4; iDau++) {
      decayP2[1][iDau] = -1;
      decayP3[1][iDau] = -1;
      decayP4[1][iDau] = -1;
    }
  }
  if (optUseDtoV0 == 0) {
    for (Int_t iDau = 0; iDau < 4; iDau++) {
      decayP3[2][iDau] = -1; // swicth off Ds->K0K+
      decayP5[0][iDau] = -1; // swicth off D+->K0pi+
    }
  }
  if (optUsePionicD0 == 0) {
    for (Int_t iDau = 0; iDau < 4; iDau++) {
      decayP5[1][iDau] = -1; // switch off D0->pi0pi+pi-
      decayP6[1][iDau] = -1; // switch off D0->pi0pi+K-
      decayP6[0][iDau] = -1; // switch off D+->K0spi+pi0
    }
  } else {
    // Pi0 options
    // - 1: both pionic and D0 -> K- pi+
    // - 2: only pionic modes, both D0 -> K- pi+ pi0 and D0 -> pi0 pi+ pi-
    // - 3: only pionic modes, only D0 -> K- pi+ pi0
    // - 4: only pionic modes, only D0 -> pi0 pi+ pi-
    switch (optUsePionicD0) {
    case 1:
      break;
    case 2: {
      // switch off D0 -> K- pi+
      for (Int_t iDau = 0; iDau < 4; iDau++) {
        decayP1[1][iDau] = -1;
      }
      break;
    }
    case 3: {
      for (Int_t iDau = 0; iDau < 4; iDau++) {
        decayP1[1][iDau] = -1; // switch off D0 -> K- pi+
        decayP5[1][iDau] = -1; // switch off D0->pi0pi+pi-
      }
      break;
    }
    case 4: {
      for (Int_t iDau = 0; iDau < 4; iDau++) {
        decayP1[1][iDau] = -1; // switch off D0 -> K- pi+
        decayP6[1][iDau] = -1; // switch off D0->pi0pi+K-
      }
      break;
    }
    default:
      break;
    };
    // For D+ force decays into pi0 + X
    for (Int_t iDau = 0; iDau < 4; iDau++) {
      // D+
      decayP1[0][iDau] = -1;
      decayP2[0][iDau] = -1;
      decayP3[0][iDau] = -1;
      decayP4[0][iDau] = -1;
      decayP5[0][iDau] = -1;
    }
  }

  for (Int_t ihadron = 0; ihadron < kNHadrons; ihadron++) {
    Int_t kc = fPythia->Pycomp(hadron[ihadron]);
    Int_t ifirst = fPythia->GetMDCY(kc, 2);
    Int_t ilast = ifirst + fPythia->GetMDCY(kc, 3) - 1;
    Double_t norm = 0.;
    for (channel = ifirst; channel <= ilast; channel++)
      norm += fPythia->GetBRAT(channel);
    if (norm < 1. - 1.e-12 || norm > 1. + 1.e-12) {
      char pName[16];
      fPythia->Pyname(hadron[ihadron], pName);
      printf(
          "Total branching ratio of %s (PDG code = %d) not equal to 1 (= %f)",
          pName, hadron[ihadron], norm);
    }
    fBraPart[kc] = norm;
    fPythia->SetMDCY(kc, 1, 1);

    for (channel = ifirst; channel <= ilast; channel++) {
      if ((fPythia->GetKFDP(channel, 1) == decayP1[ihadron][0] &&
           fPythia->GetKFDP(channel, 2) == decayP1[ihadron][1] &&
           fPythia->GetKFDP(channel, 3) == decayP1[ihadron][2] &&
           fPythia->GetKFDP(channel, 4) == decayP1[ihadron][3] &&
           fPythia->GetKFDP(channel, 5) == 0) ||
          (fPythia->GetKFDP(channel, 1) == decayP2[ihadron][0] &&
           fPythia->GetKFDP(channel, 2) == decayP2[ihadron][1] &&
           fPythia->GetKFDP(channel, 3) == decayP2[ihadron][2] &&
           fPythia->GetKFDP(channel, 4) == decayP2[ihadron][3] &&
           fPythia->GetKFDP(channel, 5) == 0) ||
          (fPythia->GetKFDP(channel, 1) == decayP3[ihadron][0] &&
           fPythia->GetKFDP(channel, 2) == decayP3[ihadron][1] &&
           fPythia->GetKFDP(channel, 3) == decayP3[ihadron][2] &&
           fPythia->GetKFDP(channel, 4) == decayP3[ihadron][3] &&
           fPythia->GetKFDP(channel, 5) == 0) ||
          (fPythia->GetKFDP(channel, 1) == decayP4[ihadron][0] &&
           fPythia->GetKFDP(channel, 2) == decayP4[ihadron][1] &&
           fPythia->GetKFDP(channel, 3) == decayP4[ihadron][2] &&
           fPythia->GetKFDP(channel, 4) == decayP4[ihadron][3] &&
           fPythia->GetKFDP(channel, 5) == 0) ||
          (fPythia->GetKFDP(channel, 1) == decayP5[ihadron][0] &&
           fPythia->GetKFDP(channel, 2) == decayP5[ihadron][1] &&
           fPythia->GetKFDP(channel, 3) == decayP5[ihadron][2] &&
           fPythia->GetKFDP(channel, 4) == decayP5[ihadron][3] &&
           fPythia->GetKFDP(channel, 5) == 0) ||
          (fPythia->GetKFDP(channel, 1) == decayP6[ihadron][0] &&
           fPythia->GetKFDP(channel, 2) == decayP6[ihadron][1] &&
           fPythia->GetKFDP(channel, 3) == decayP6[ihadron][2] &&
           fPythia->GetKFDP(channel, 4) == decayP6[ihadron][3] &&
           fPythia->GetKFDP(channel, 5) == 0))

      {
        fPythia->SetMDME(channel, 1, 1);
      } else {
        fPythia->SetMDME(channel, 1, 0);
        fBraPart[kc] -= fPythia->GetBRAT(channel);
      } // selected channel ?
    }   // decay channels
    if (norm > 0)
      fBraPart[kc] /= norm;
  } // hadrons

  // Options for forcing Lc decays for dedicated productions
  Int_t prodLcpKpi[3] = {2212, 321, 211};
  Int_t multLcpKpi[3] = {1, 1, 1};

  Int_t prodLcpK0S[2] = {2212, 311};
  Int_t multLcpK0S[2] = {1, 1};

  Int_t prodLcLambdaPiPlPi0[3] = {iLambda, kPiPlus, kPi0};
  Int_t multLcLambdaPiPlPi0[3] = {1, 1, 1};
  if (optForceLcChannel == 1) { // pKpi
    ForceParticleDecay(4122, prodLcpKpi, multLcpKpi, 3, 1);
  }
  if (optForceLcChannel == 2) {                             // pK0S
    ForceParticleDecay(4122, prodLcpK0S, multLcpK0S, 2, 1); // Lc to p + K0
  }
  if (optForceLcChannel == 3) { // Lambda Pi+ Pi0
    ForceParticleDecay(4122, prodLcLambdaPiPlPi0, multLcLambdaPiPlPi0, 3, 1);
  }

  if (optForceLcChannel ==
      4) { // pK0S for BDT signal training: force all K0->K0S->pi+pi-
    ForceParticleDecay(4122, prodLcpK0S, multLcpK0S, 2, 1); // Lc to p + K0
    ForceParticleDecay(311, 310, 1);                        // K0 -> K0S
    ForceParticleDecay(310, 211, 2);                        // K0S -> pi+ pi-
  }
}

void PythiaDecayerConfig::ForceParticleDecay(Int_t particle, Int_t product,
                                             Int_t mult) {
  //
  //  Force decay of particle into products with multiplicity mult

  Int_t kc = fPythia->Pycomp(particle);
  fPythia->SetMDCY(kc, 1, 1);
  Int_t ifirst = fPythia->GetMDCY(kc, 2);
  Int_t ilast = ifirst + fPythia->GetMDCY(kc, 3) - 1;
  Double_t norm = 0.;
  for (Int_t channel = ifirst; channel <= ilast; channel++) {
    norm += fPythia->GetBRAT(channel);
  }
  if (norm < 1. - 1.e-12 || norm > 1. + 1.e-12) {
    char pName[16];
    fPythia->Pyname(particle, pName);
    printf("Total branching ratio of %s (PDG code = %d) not equal to 1 (= %f)",
           pName, particle, norm);
  }
  fBraPart[kc] = norm;
  //
  //  Loop over decay channels
  for (Int_t channel = ifirst; channel <= ilast; channel++) {
    if (CountProducts(channel, product) >= mult) {
      fPythia->SetMDME(channel, 1, 1);
    } else {
      fPythia->SetMDME(channel, 1, 0);
      fBraPart[kc] -= fPythia->GetBRAT(channel);
    }
  }
  if (norm > 0.)
    fBraPart[kc] /= norm;
}

void PythiaDecayerConfig::ForceParticleDecay(Int_t particle,
                                             const Int_t *products, Int_t *mult,
                                             Int_t npart, Bool_t flag) {
  //
  //  Force decay of particle into products with multiplicity mult

  Int_t kc = fPythia->Pycomp(particle);
  fPythia->SetMDCY(kc, 1, 1);
  Int_t ifirst = fPythia->GetMDCY(kc, 2);
  Int_t ilast = ifirst + fPythia->GetMDCY(kc, 3) - 1;
  Double_t norm = 0.;
  for (Int_t channel = ifirst; channel <= ilast; channel++)
    norm += fPythia->GetBRAT(channel);
  if (norm < 1. - 1.e-12 || norm > 1. + 1.e-12) {
    char pName[16];
    fPythia->Pyname(particle, pName);
    printf("Total branching ratio of %s (PDG code = %d) not equal to 1 (= %f)",
           pName, particle, norm);
  }
  fBraPart[kc] = norm;
  //
  //  Loop over decay channels
  for (Int_t channel = ifirst; channel <= ilast; channel++) {
    Int_t nprod = 0;
    for (Int_t i = 0; i < npart; i++) {
      nprod += (CountProducts(channel, products[i]) >= mult[i]);
    }
    if ((nprod && !flag) || ((nprod == npart) && flag)) {
      fPythia->SetMDME(channel, 1, 1);
    } else { //
      fPythia->SetMDME(channel, 1, 0);
      fBraPart[kc] -= fPythia->GetBRAT(channel);
    }
  }
  if (norm > 0.)
    fBraPart[kc] /= norm;
}

void PythiaDecayerConfig::ForceOmega() {
  // Force Omega -> Lambda K- Decay
  Int_t kc = fPythia->Pycomp(3334);
  fPythia->SetMDCY(kc, 1, 1);
  Int_t ifirst = fPythia->GetMDCY(kc, 2);
  Int_t ilast = ifirst + fPythia->GetMDCY(kc, 3) - 1;
  for (Int_t channel = ifirst; channel <= ilast; channel++) {
    if (fPythia->GetKFDP(channel, 1) == kLambda0 &&
        fPythia->GetKFDP(channel, 2) == kKMinus &&
        fPythia->GetKFDP(channel, 3) == 0) {
      fPythia->SetMDME(channel, 1, 1);
    } else {
      fPythia->SetMDME(channel, 1, 0);
    } // selected channel ?
  }   // decay channels
}

void PythiaDecayerConfig::ForceLambda() {
  // Force Lambda -> p pi-
  Int_t kc = fPythia->Pycomp(3122);
  fPythia->SetMDCY(kc, 1, 1);
  Int_t ifirst = fPythia->GetMDCY(kc, 2);
  Int_t ilast = ifirst + fPythia->GetMDCY(kc, 3) - 1;
  for (Int_t channel = ifirst; channel <= ilast; channel++) {
    if (fPythia->GetKFDP(channel, 1) == kProton &&
        fPythia->GetKFDP(channel, 2) == kPiMinus &&
        fPythia->GetKFDP(channel, 3) == 0) {
      fPythia->SetMDME(channel, 1, 1);
    } else {
      fPythia->SetMDME(channel, 1, 0);
    } // selected channel ?
  }   // decay channels
}

void PythiaDecayerConfig::SwitchOffBDecay() {
  // Switch off B-decays
  Int_t heavyB[] = {511, 521, 531, 5122, 5132, 5232, 5332};
  for (int i = 0; i < 4; i++) {
    fPythia->SetMDCY(fPythia->Pycomp(heavyB[i]), 1, 0);
  }
}

Float_t PythiaDecayerConfig::GetPartialBranchingRatio(Int_t kf) {
  // Get branching ratio
  Int_t kc = fPythia->Pycomp(TMath::Abs(kf));
  return fBraPart[kc];
}

Float_t PythiaDecayerConfig::GetLifetime(Int_t kf) {
  // Get branching ratio
  Int_t kc = fPythia->Pycomp(TMath::Abs(kf));
  return fPythia->GetPMAS(kc, 4) * 3.3333e-12;
}

Int_t PythiaDecayerConfig::ImportParticles(TClonesArray *particles)
{
   return fPythia->ImportParticles(particles,"All");
}

void PythiaDecayerConfig::ReadDecayTable()
{
   if (fDecayTableFile.IsNull()) {
      printf("PythiaDeceayerConfig: Warning: ReadDecayTable: No file set\n");
      return;
   }
   Int_t lun = 15;
   fPythia->OpenFortranFile(lun,const_cast<char*>(fDecayTableFile.Data()));
   fPythia->Pyupda(3,lun);
   fPythia->CloseFortranFile(lun);
}

void PythiaDecayerConfig::Decay(Int_t idpart, TLorentzVector* p) {
  if (!p) return;

  Float_t energy = p->Energy();
  Float_t theta  = p->Theta();
  Float_t phi    = p->Phi();
  
  if(!fDecayerExodus) {
    fPythia->Py1ent(0, idpart, energy, theta, phi);
  } else {

  // EXODUS decayer
    fPythia->SetMSTU(10,1);
    fPythia->SetP(1,5,p->M()) ;
    if(idpart == 111){
      fPythia->SetMDCY(fPythia->Pycomp(22) ,1, 0);
      fPythia->Py1ent(0, idpart, energy, theta, phi);
      PizeroDalitz();
      fPythia->SetMDCY(fPythia->Pycomp(22) ,1, 1);
    }
    else if(idpart == 221){
      fPythia->SetMDCY(fPythia->Pycomp(22) ,1, 0);
      fPythia->Py1ent(0, idpart, energy, theta, phi);
      EtaDalitz();
      fPythia->SetMDCY(fPythia->Pycomp(22) ,1, 1);
    }
    else if(idpart == 113){
      fPythia->Py1ent(0, idpart, energy, theta, phi);
      RhoDirect();
    }
    else if(idpart == 223){
      fPythia->SetMDCY(fPythia->Pycomp(111) ,1, 0);
      fPythia->Py1ent(0, idpart, energy, theta, phi);
      OmegaDirect();
      OmegaDalitz();
      fPythia->SetMDCY(fPythia->Pycomp(111) ,1, 1);
    }
    else if(idpart == 331){
      fPythia->SetMDCY(fPythia->Pycomp(22) ,1, 0);
      fPythia->SetMDCY(fPythia->Pycomp(223) ,1, 0);
      fPythia->Py1ent(0, idpart, energy, theta, phi);
      EtaprimeDalitz();
      fPythia->SetMDCY(fPythia->Pycomp(223) ,1, 1);
      fPythia->SetMDCY(fPythia->Pycomp(22) ,1, 1);
    }
    else if(idpart == 333){
      fPythia->SetMDCY(fPythia->Pycomp(221) ,1, 0);
      fPythia->SetMDCY(fPythia->Pycomp(111) ,1, 0);
      fPythia->Py1ent(0, idpart, energy, theta, phi);
      PhiDirect();
      PhiDalitz();
      fPythia->SetMDCY(fPythia->Pycomp(111) ,1, 1);
      fPythia->SetMDCY(fPythia->Pycomp(221) ,1, 1);
    }
    else if(idpart == 443){
      fPythia->Py1ent(0, idpart, energy, theta, phi);
      JPsiDirect();
    }  
  }
  fPythia->SetMSTU(10,2);
  fPythia->GetPrimaries();
}

void PythiaDecayerConfig::PizeroDalitz(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 111) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 2) continue;
    if ((fPythia->GetK(fd+1,2) != 22) || (TMath::Abs(fPythia->GetK(fd+2,2)) != 11) ) continue;
    TLorentzVector pizero(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg+1000*fPythia->GetK(fd+1,2), &pizero);
    for (Int_t j = 0; j < 3; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_pion())[2-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::EtaDalitz(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 221) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 2) continue;
    if ((fPythia->GetK(fd+1,2) != 22) || (TMath::Abs(fPythia->GetK(fd+2,2)) != 11) ) continue;
    TLorentzVector eta(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg+1000*fPythia->GetK(fd+1,2), &eta);
    for (Int_t j = 0; j < 3; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_eta())[2-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::RhoDirect(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 113) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 1) continue;
    if ((TMath::Abs(fPythia->GetK(fd+1,2)) != 11)) continue;
    TLorentzVector rho(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg, &rho);
    for (Int_t j = 0; j < 2; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_rho())[1-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::OmegaDalitz(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 223) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 2) continue;
    if ((fPythia->GetK(fd+1,2) != 111) || (TMath::Abs(fPythia->GetK(fd+2,2)) != 11)) continue;
    TLorentzVector omegadalitz(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg+1000*fPythia->GetK(fd+1,2), &omegadalitz);
    for (Int_t j = 0; j < 3; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_omega_dalitz())[2-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::OmegaDirect(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 223) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 1) continue;
    if ((TMath::Abs(fPythia->GetK(fd+1,2)) != 11)) continue;
    TLorentzVector omegadirect(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg, &omegadirect);
    for (Int_t j = 0; j < 2; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_omega())[1-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::EtaprimeDalitz(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 331) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 2) continue;
    if ((fPythia->GetK(fd+1,2) != 22 && fPythia->GetK(fd+1,2) != 223) || (TMath::Abs(fPythia->GetK(fd+2,2)) != 11)) continue;
    TLorentzVector etaprime(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg+1000*fPythia->GetK(fd+1,2), &etaprime);
    for (Int_t j = 0; j < 3; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_etaprime())[2-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::PhiDalitz(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 333) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 2) continue;
    if ((fPythia->GetK(fd+1,2) != 221 && fPythia->GetK(fd+1,2) != 111 ) || (TMath::Abs(fPythia->GetK(fd+2,2)) != 11)) continue;
    TLorentzVector phidalitz(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg+1000*fPythia->GetK(fd+1,2), &phidalitz);
    for (Int_t j = 0; j < 3; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_phi_dalitz())[2-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::PhiDirect(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 333) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 1) continue;
    if ((TMath::Abs(fPythia->GetK(fd+1,2)) != 11)) continue;
    TLorentzVector phi(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg, &phi);
    for (Int_t j = 0; j < 2; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_phi())[1-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::JPsiDirect(){
  Int_t nt = fPythia->GetN();
  for (Int_t i = 0; i < nt; i++) {
    if (fPythia->GetK(i+1,2) != 443) continue;
    Int_t fd = fPythia->GetK(i+1,4) - 1;
    Int_t ld = fPythia->GetK(i+1,5) - 1;
    if (fd < 0) continue;
    if ((ld - fd) != 1) continue;
    if ((TMath::Abs(fPythia->GetK(fd+1,2)) != 11)) continue;
    TLorentzVector jpsi(fPythia->GetP(i+1,1), fPythia->GetP(i+1,2), fPythia->GetP(i+1,3), fPythia->GetP(i+1,4));
    Int_t pdg = TMath::Abs(fPythia->GetK(i+1,2));
    fDecayerExodus->Decay(pdg, &jpsi);
    for (Int_t j = 0; j < 2; j++) {
      for (Int_t k = 0; k < 4; k++) {
        TLorentzVector vec = (fDecayerExodus->Products_jpsi())[1-j];
        fPythia->SetP(fd+j+1,k+1,vec[k]);
      }
    }
  }
}

void PythiaDecayerConfig::Copy(TObject &) const {
  //
  // Copy *this onto PythiaDecayerConfig -- not implemented
  //
  Fatal("Copy", "Not implemented!\n");
}
