#ifndef PYTHIADECAYERCONFIG_H
#define PYTHIADECAYERCONFIG_H
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

#include <TPythia6.h>

typedef enum {
  kBSemiElectronic,
  kSemiElectronic,
  kDiElectron,
  kBSemiMuonic,
  kDSemiMuonic,
  kSemiMuonic,
  kDiMuon,
  kJpsiDiMuon,
  kBJpsiDiMuon,
  kBJpsiDiElectron,
  kBPsiPrimeDiMuon,
  kBPsiPrimeDiElectron,
  kPiToMu,
  kKaToMu,
  kNoDecay,
  kHadronicD,
  kHadronicDWithout4Bodies,
  kOmega,
  kLambda,
  kPhiKK,
  kAll,
  kNoDecayHeavy,
  kHardMuons,
  kBJpsi,
  kBJpsiUndecayed,
  kWToMuon,
  kWToCharm,
  kWToCharmToMuon,
  kZDiMuon,
  kZDiElectron,
  kNeutralPion,
  kAllMuonic,
  kChiToJpsiGammaToMuonMuon,
  kChiToJpsiGammaToElectronElectron,
  kNoDecayBeauty,
  kPsiPrimeJpsiDiElectron,
  kElectronEM,
  kGammaEM,
  kDiElectronEM,
  kBeautyUpgrade,
  kHadronicDWithV0,
  kHadronicDWithout4BodiesWithV0,
  kAllEM,
  kLcpKpi,
  kLcpK0S,
  kHFYellowReport,
  kHadronicDPionicD0,
  kHadronicDWithV0PionicD0,
  kHadronicDWithout4BodiesPionicD0,
  kHadronicDWithout4BodiesWithV0PionicD0,
  kHadronicDPionicD0pure,
  kHadronicDPionicD0K,
  kHadronicDPionicD0pi,
  kLcpK0SBDTsig,
  kEtaPrime,
  kXic0Semileptonic,
  kHadronicDWithout4BodiesDsPhiPi
} Decay_t;

class PythiaDecayerConfig : public TObject {
public:
  PythiaDecayerConfig();
  PythiaDecayerConfig(const PythiaDecayerConfig &decayerconfig);
  //
  virtual ~PythiaDecayerConfig() { ; }
  virtual void Init(Decay_t decay);
  virtual void ForceDecay();
  virtual void SetPatchOmegaDalitz() { fPatchOmegaDalitz = 1; }
  virtual void SetDecayerExodus() { fDecayerExodus = 1; }
  virtual void HeavyFlavourOff() { fHeavyFlavour = kFALSE; }
  virtual void DecayLongLivedParticles() { fLongLived = kTRUE; }
  virtual Float_t GetPartialBranchingRatio(Int_t ipart);
  virtual Float_t GetLifetime(Int_t kf);
  virtual void SwitchOffBDecay();
  virtual void SwitchOffPi0() { fPi0 = 0; }
  virtual void SwitchOffParticle(Int_t kf);

private:
  Int_t CountProducts(Int_t channel, Int_t particle);
  void ForceParticleDecay(Int_t particle, Int_t product, Int_t mult);
  void ForceParticleDecay(Int_t particle, const Int_t *products, Int_t *mult,
                          Int_t npart, Bool_t flag = 0);
  void ForceHadronicD(Int_t optUse4Bodies = 1, Int_t optUseDtoV0 = 0,
                      Int_t optForceLcChannel = 0, Int_t optUsePionicD0 = 0);
  void ForceOmega();
  void ForceLambda();
  void SwitchOffHeavyFlavour();
  void ForceBeautyUpgrade();
  void ForceHFYellowReport();
  Float_t GetBraPart(Int_t kf);
  void Copy(TObject &decayer) const;

  PythiaDecayerConfig &operator=(const PythiaDecayerConfig &decayerconfig) {
    decayerconfig.Copy(*this);
    return (*this);
  }

private:
  TPythia6 *fPythia;        //! Pointer to AliPythia
  Decay_t fDecay;           //  Forced decay mode
  Float_t fBraPart[501];    //! Branching ratios
  Bool_t fHeavyFlavour;     //! Flag for heavy flavors
  Bool_t fLongLived;        //! Flag for long lived particle decay
  Bool_t fPatchOmegaDalitz; //! Flag to patch the omega Dalitz decays
  Bool_t fDecayerExodus;    //! Flag for EXODUS decayer
  Bool_t fPi0;              //! Flag for pi0 decay
  static Bool_t fgInit;     //! initialization flag

  ClassDef(PythiaDecayerConfig, 1) // AliDecayer implementation using Pythia
};
#endif
