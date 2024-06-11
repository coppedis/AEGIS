#ifndef GENERATORPARAM_H
#define GENERATORPARAM_H
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
// Class to generate particles using parametrized pT and y distributions.
// Distributions are obtained from pointer to object of type AliGenLib.
//
// andreas.morsch@cern.ch
//
#include "GeneratorParamLibBase.h"
#include "PythiaDecayerConfig.h"
#include <TArrayF.h>
#include <TArrayI.h>
#include <TGenerator.h>
#include <TMath.h>
#include <TVector3.h>
#include <TVirtualMCDecayer.h>
#include <map>
class TF1;
typedef enum { kNoSmear, kPerEvent, kPerTrack } VertexSmear_t;
typedef enum { kAnalog, kNonAnalog } Weighting_t;

//-------------------------------------------------------------
class GeneratorParam : public TGenerator {
public:
  GeneratorParam();
  GeneratorParam(Int_t npart, const GeneratorParamLibBase *Library, Int_t param,
                 const char *tname = 0);
  GeneratorParam(Int_t npart, Int_t param, const char *tname = 0,
                 const char *name = 0);
  GeneratorParam(Int_t npart, Int_t param,
                 Double_t (*PtPara)(const Double_t *, const Double_t *),
                 Double_t (*YPara)(const Double_t *, const Double_t *),
                 Double_t (*V2Para)(const Double_t *, const Double_t *),
                 Int_t (*IpPara)(TRandom *));
  GeneratorParam(const char *name, Int_t npart, int pdg,
                 Double_t (*PtPara)(const Double_t *, const Double_t *) = 0x0,
                 Double_t (*YPara)(const Double_t *, const Double_t *) = 0x0,
                 Double_t (*V2Para)(const Double_t *, const Double_t *) = 0x0);

  virtual ~GeneratorParam();
  virtual void GenerateEvent();
  virtual void Init();
  virtual int ImportParticles(TClonesArray *particles, Option_t *option);
  // select particle type
  virtual void SetParam(Int_t param) { fParam = param; }
  // Setting the flag for Background transportation while using SetForceDecay()
  void SetSelectAll(Bool_t selectall) { fSelectAll = selectall; }
  virtual void SetMomentumRange(Float_t pmin = 0, Float_t pmax = 1.e10);
  virtual void SetPtRange(Float_t ptmin = 0, Float_t ptmax = 1.e10);
  virtual void SetPhiRange(Float_t phimin = 0., Float_t phimax = 360.);
  virtual void SetYRange(Float_t ymin = -100, Float_t ymax = 100);
  virtual void SetThetaRange(Float_t thetamin = 0, Float_t thetamax = 180);
  virtual void SetNumberParticles(Int_t npart = 100) { fNpart = npart; }
  virtual Int_t NumberParticles() const { return fNpart; }
  virtual Bool_t KinematicSelection(const TParticle *particle,
                                    Int_t flag) const;
  // Kinematic cuts on decay products
  virtual void SetForceDecay(Decay_t decay = kAll) { fForceDecay = decay; }
  virtual void SetCutOnChild(Int_t flag = 0) { fCutOnChild = flag; }
  virtual void SetChildMomentumRange(Float_t pmin = 0, Float_t pmax = 1.e10) {
    fChildPMin = pmin;
    fChildPMax = pmax;
  }
  virtual void SetChildPtRange(Float_t ptmin = 0, Float_t ptmax = 20.) {
    fChildPtMin = ptmin;
    fChildPtMax = ptmax;
  }
  virtual void SetChildPhiRange(Float_t phimin = 0., Float_t phimax = 360.) {
    fChildPhiMin = TMath::Pi() * phimin / 180;
    fChildPhiMax = TMath::Pi() * phimax / 180;
  }
  virtual void SetChildThetaRange(Float_t thetamin = 0,
                                  Float_t thetamax = 180) {
    fChildThetaMin = TMath::Pi() * thetamin / 180;
    fChildThetaMax = TMath::Pi() * thetamax / 180;
  }
  virtual void SetChildYRange(Float_t ymin = -12, Float_t ymax = 12) {
    fChildYMin = ymin;
    fChildYMax = ymax;
  }
  virtual void SetMaximumLifetime(Float_t time = 1.e-15) {
    fMaxLifeTime = time;
  }
  // force decay type
  virtual void SetDeltaPt(Float_t delta = 0.01) { fDeltaPt = delta; }
  virtual void SetDecayer(TVirtualMCDecayer *decayer) { fDecayer = decayer; }
  virtual void SetForceGammaConversion(Bool_t force = kTRUE) {
    fForceConv = force;
  }
  virtual void SetKeepParent(Bool_t keep = kTRUE) {
    fKeepParent = keep;
  } // Store parent even if it does not have childs within cuts
  virtual void SetKeepIfOneChildSelected(Bool_t keep = kTRUE) {
    fKeepIfOneChildSelected = keep;
  } // Accept parent and child even if other children are not within cut.
  virtual void SetPreserveFullDecayChain(Int_t preserve = kFALSE) {
    fPreserveFullDecayChain = preserve;
  } // Prevent flagging(/skipping) of decay daughter particles; preserves
    // complete forced decay chain

  virtual void SetWeighting(Weighting_t flag = kAnalog) {fAnalog = flag;}

  virtual void Draw(const char *opt);
  TF1 *GetPt() { return fPtPara; }
  TF1 *GetY() { return fYPara; }
  Float_t GetRelativeArea(Float_t ptMin, Float_t ptMax, Float_t yMin,
                          Float_t yMax, Float_t phiMin, Float_t phiMax);

  static TVector3 OrthogonalVector(TVector3 &inVec);
  static void RotateVector(Double_t *pin, Double_t *pout, Double_t costheta,
                           Double_t sintheta, Double_t cosphi, Double_t sinphi);
  static double ScreenFunction1(double d);
  static double ScreenFunction2(double d);
  double RandomEnergyFraction(double Z, double E);
  double RandomPolarAngle();
  double RandomMass(Double_t mh);
  Int_t VirtualGammaPairProduction(TClonesArray *particles, Int_t nPart);
  Int_t ForceGammaConversion(TClonesArray *particles, Int_t nPart);
  virtual void SetSeed(UInt_t /*seed*/) { ; }

  // allow explicit setting of functions in case of streaming
  void SetParamsExplicitly(const GeneratorParamLibBase *Library, Int_t param,
                           const char *tname) {
    fPtParaFunc = Library->GetPt(param, tname);
    fYParaFunc = Library->GetY(param, tname);
    fIpParaFunc = Library->GetIp(param, tname);
    fV2ParaFunc = Library->GetV2(param, tname);
  }

  // retrive particle type
  Int_t GetParam() { return fParam; }

protected:
  Double_t (*fPtParaFunc)(
      const Double_t *,
      const Double_t *); //! Pointer to Pt parametrisation function
  Double_t (*fYParaFunc)(
      const Double_t *,
      const Double_t *); //! Pointer to Y parametrisation function
  Int_t (*fIpParaFunc)(
      TRandom *); //! Pointer to particle type parametrisation function
  Double_t (*fV2ParaFunc)(
      const Double_t *,
      const Double_t *);     //! Pointer to V2 parametrisation function
  TF1 *fPtPara = 0;          // Transverse momentum parameterisation
  TF1 *fYPara = 0;           // Rapidity parameterisation
  TF1 *fV2Para = 0;          // v2 parametrization
  TF1 *fdNdPhi = 0;          // Phi distribution depending on v2
  Int_t fParam = 0;          // Parameterisation type
  Float_t fDeltaPt = 0.01;   // pT sampling in steps of fDeltaPt
  Bool_t fSelectAll = false; // Flag for transportation of Background while
                             // using SetForceDecay()
  TVirtualMCDecayer *fDecayer = 0; // ! Pointer to virtual decyer
  Bool_t fForceConv = false;       // force converson of gammas
  Bool_t fKeepParent =
      false; //  Store parent even if it does not have childs within cuts
  Bool_t fKeepIfOneChildSelected = true; // Accept parent and child even if
                                         // other children are not within cut.
  Bool_t fPreserveFullDecayChain =
      true; // Prevent flagging(/skipping) of decay daughter particles;
            // preserves complete forced decay chain
  Int_t fPDGcode = -1;    // PDG code in case of single particle injector
  Int_t fIncFortran = -1; // respect fortran counting in particle list (Pythia6)
                          // kinematic cuts
  Float_t fPtMin = 0.;
  Float_t fPtMax = 0.;
  Float_t fPMin = 0.;
  Float_t fPMax = 1.e10;
  Float_t fPhiMin = 0.;
  Float_t fPhiMax = 2. * TMath::Pi();
  Float_t fYMin = 0.;
  Float_t fYMax = 0.;
  Float_t fThetaMin = 0.;
  Float_t fThetaMax = TMath::Pi();
  Float_t fMaxLifeTime = 0.;
  Int_t fNpart = 0;
  Float_t fTimeOrigin = 0.;
  Float_t fTime = 0.;
  Float_t fEvPlane = 0.;
  Decay_t fForceDecay = kAll;
  Int_t fNprimaries = 0;
  Bool_t fCutOnChild = 0.;
  Float_t fChildPMin = 0.;
  Float_t fChildPMax = 0.;
  Float_t fChildPtMin = 0.;
  Float_t fChildPtMax = 0.;
  Float_t fChildYMin = 0.;
  Float_t fChildYMax = 0.;
  Float_t fChildThetaMin = 0.;
  Float_t fChildThetaMax = 0.;
  Float_t fChildPhiMin = 0.;
  Float_t fChildPhiMax = 0.;
  Float_t fParentWeight = 1.;
  Float_t fChildWeight = 1.;
  Float_t fYWgt = 1.;
  Float_t fPtWgt = 1.;
  Float_t fdNdy0 = 1.;
  Weighting_t fAnalog = kAnalog;
  
  TArrayI fChildSelect; //! Decay products to be selected
  enum {
    kThetaRange = BIT(14),
    kPhiRange = BIT(16),
    kPtRange = BIT(17),
    kYRange = BIT(18),
    kMomentumRange = BIT(19),
    kEtaRange = BIT(20)
  };

  std::map<int, std::unique_ptr<TF1>> fPDGtoTF1; //! map of cache TF1 objects for "exodus"

private:
  void InitChildSelect();
  GeneratorParam(const GeneratorParam &Param);
  GeneratorParam &operator=(const GeneratorParam &rhs);

  ClassDef(GeneratorParam, 2) // Generator using parameterised pt- and y-distribution
};
#endif
