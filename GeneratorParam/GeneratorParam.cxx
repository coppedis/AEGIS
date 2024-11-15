// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Class to generate particles using paramtrized pT and y distributions.
// Distributions are obtained from pointer to object of type
// GeneratorParamLibBase. Decays are performed using Pythia.
// andreas.morsch@cern.ch

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TPythia6Decayer.h>
#include <TROOT.h>
#include <TRandom.h>
#include <vector>

#include "GeneratorParam.h"
#include "GeneratorParamLibBase.h"

ClassImp(GeneratorParam)
    //____________________________________________________________
    GeneratorParam::GeneratorParam()
    : TGenerator() {
  // Default constructor
}
//____________________________________________________________
GeneratorParam::GeneratorParam(Int_t npart,
                               const GeneratorParamLibBase *Library,
                               Int_t param, const char *tname)
    : TGenerator("GeneratorParam", "GeneratorParam"), fNpart(npart), fParam(param) {
  // Constructor using number of particles parameterisation id and library
  fName = "Param";
  fTitle = "Particle Generator using pT and y parameterisation";
  fPtParaFunc = Library->GetPt(param, tname);
  fYParaFunc = Library->GetY(param, tname);
  fIpParaFunc = Library->GetIp(param, tname);
  fV2ParaFunc = Library->GetV2(param, tname);
}

//____________________________________________________________

GeneratorParam::GeneratorParam(
    Int_t npart, Int_t param,
    Double_t (*PtPara)(const Double_t *, const Double_t *),
    Double_t (*YPara)(const Double_t *, const Double_t *),
    Double_t (*V2Para)(const Double_t *, const Double_t *),
    Int_t (*IpPara)(TRandom *))
    : TGenerator("GeneratorParam", "GeneratorParam"), fNpart(npart),
      fPtParaFunc(PtPara), fYParaFunc(YPara), fIpParaFunc(IpPara),
      fV2ParaFunc(V2Para), fParam(param) {
  // Constructor
  fName = "Param";
  fTitle = "Particle Generator using pT and y parameterisation";
  SetChildMomentumRange();
  SetChildPtRange();
  SetChildPhiRange();
  SetChildThetaRange();
}

//____________________________________________________________
GeneratorParam::GeneratorParam(
    const char *name, Int_t npart, int pdg,
    Double_t (*PtPara)(const Double_t *, const Double_t *),
    Double_t (*YPara)(const Double_t *, const Double_t *),
    Double_t (*V2Para)(const Double_t *, const Double_t *))
    : TGenerator("GeneratorParam", "GeneratorParam"), fNpart(npart),
      fPtParaFunc(PtPara), fYParaFunc(YPara), fV2ParaFunc(V2Para),
      fPDGcode(pdg) {
  // Constructor
  fName = name;
  fTitle =
      Form("Particle Generator using pT and y parameterisation for %s", name);

  SetChildMomentumRange();
  SetChildPtRange();
  SetChildPhiRange();
  SetChildThetaRange();
  fPtParaFunc = PtPara ? PtPara : [](const double *, const double *) -> double {
    return 1.;
  };
  fYParaFunc = YPara ? YPara : [](const double *, const double *) -> double {
    return 1.;
  };
  fV2ParaFunc = V2Para ? V2Para : [](const double *, const double *) -> double {
    return 0.;
  };
}

//____________________________________________________________
GeneratorParam::~GeneratorParam() {
  // Destructor
  delete fPtPara;
  delete fYPara;
  delete fV2Para;
  delete fdNdPhi;
}

//____________________________________________________________
void GeneratorParam::Init() {
  // Initialisation

  if (!fDecayer) {
    Fatal("Init", "No decayer set \n");
  }
  // For AliDecayerPythia based on FORTRAN there is a known inconsitency between
  // the index returned by GetFirstDaughter(), GetLastDaughter(), GetMother()
  // ... and the actual position in the C++ particle list. That's the reason whe
  // the offset fIncFortran = -1 is needed in this case but not for decayers
  // based on C++

  if (strcmp(fDecayer->ClassName(), "TPythia6Decayer"))
    fIncFortran = -1;
  char name[256];
  snprintf(name, 256, "pt-parameterisation for %s", GetName());

  if (fPtPara)
    fPtPara->Delete();
  fPtPara = new TF1(name, fPtParaFunc, fPtMin, fPtMax, 0);
  gROOT->GetListOfFunctions()->Remove(fPtPara);
  //  Set representation precision to 10 MeV
  Int_t npx = Int_t((fPtMax - fPtMin) / fDeltaPt);
  fPtPara->SetNpx(npx);

  snprintf(name, 256, "y-parameterisation  for %s", GetName());
  if (fYPara)
    fYPara->Delete();
  fYPara = new TF1(name, fYParaFunc, fYMin, fYMax, 0);
  gROOT->GetListOfFunctions()->Remove(fYPara);

  snprintf(name, 256, "v2-parameterisation for %s", GetName());
  if (fV2Para)
    fV2Para->Delete();
  fV2Para = new TF1(name, fV2ParaFunc, fPtMin, fPtMax, 0);
  snprintf(name, 256, "dNdPhi for %s", GetName());
  if (fdNdPhi)
    fdNdPhi->Delete();
  fdNdPhi = new TF1(name, "1+2*[0]*TMath::Cos(2*(x-[1]))", fPhiMin, fPhiMax);

  //
  //
  snprintf(name, 256, "pt-for-%s", GetName());
  TF1 ptPara(name ,fPtParaFunc, 0, 15, 0);
  snprintf(name, 256, "y-for-%s", GetName());
  TF1 yPara(name, fYParaFunc, -6, 6, 0);
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
  Float_t intYS  = yPara.Integral(fYMin, fYMax,(Double_t*) 0x0,1.e-6);
  Float_t intPt0 = ptPara.Integral(0,15,(Double_t *) 0x0,1.e-6);
  Float_t intPtS = ptPara.Integral(fPtMin,fPtMax,(Double_t*) 0x0,1.e-6);
#else
  Float_t intYS  = yPara.Integral(fYMin, fYMax,1.e-6);
  Float_t intPt0 = ptPara.Integral(0,15,1.e-6);
  Float_t intPtS = ptPara.Integral(fPtMin,fPtMax,1.e-6);
#endif
  Float_t phiWgt=(fPhiMax-fPhiMin)/TMath::TwoPi();    //TR: should probably be done differently in case of anisotropic phi...

  //                                                                                                                                   // dN/dy| y=0
  Double_t y1=0;
  Double_t y2=0;

  fdNdy0=fYParaFunc(&y1,&y2);

  fYWgt  = intYS/fdNdy0;
  if (fAnalog == kAnalog) {
    fPtWgt = intPtS/intPt0;
  } else {
    fPtWgt = (fPtMax-fPtMin)/intPt0;
  }
  fParentWeight = fYWgt*fPtWgt*phiWgt/fNpart;
  //
  //
  // Initialize the decayer
  fDecayer->SetForceDecay(fForceDecay);
  fDecayer->Init();
  // initialise selection of decay products
  InitChildSelect();
}

void GeneratorParam::GenerateEvent() {
  //
  // Generate one event
  //
  fParticles->Clear();
 
  Float_t polar[3] = {
      0, 0, 0}; // Polarisation of the parent particle (for GEANT tracking)
  Double_t origin0[3]; // Origin of the generated parent particle (for GEANT
                       // tracking)
  Double_t och0[3];    // Origin of child particle
  Double_t pt, pl,
      ptot; // Transverse, logitudinal and total momenta of the parent particle
  Double_t phi,
      theta; // Phi and theta spherical angles of the parent particle momentum
  Double_t p[3], pc[3], och[3]; // Momentum, polarisation and origin of the
                                // children particles from lujet
  Double_t ty, xmt;
  Int_t i, j;
  Double_t energy;
  Float_t time0;
  Float_t  wgtp, wgtch;
  std::vector<bool> vFlags;
  std::vector<bool> vSelected;
  std::vector<int> vParent;
  Double_t dummy;

  // array to store decay products
  static TClonesArray *particles;
  //
  if (!particles)
    particles = new TClonesArray("TParticle", 1000);

  TDatabasePDG *pDataBase = TDatabasePDG::Instance();
  //
  Float_t random[6];

  Int_t ipa = 0;

  // Generating fNpart particles
  fNprimaries = 0;
  auto nt = 0;
  while (ipa < fNpart) {
    while (1) {
      //
      // particle type
      Int_t pdg = fIpParaFunc ? fIpParaFunc(gRandom) : fPDGcode;
      Int_t iTemp = pdg;

      // custom pdg codes to destinguish direct photons
      if ((pdg >= 220000) && (pdg <= 220001)) {
        pdg = 22;
      }
      fChildWeight=(fDecayer->GetPartialBranchingRatio(pdg))*fParentWeight;
      TParticlePDG *particle = pDataBase->GetParticle(pdg);
      Float_t am = particle->Mass();
      gRandom->RndmArray(2, random);

      // --- For Exodus -------------------------------
      Double_t awidth = particle->Width();
      if (awidth > 0) {
	TF1* rbw = nullptr;
	auto iter = fPDGtoTF1.find(pdg);
        if (iter != fPDGtoTF1.end()) {
	  // see if we have cached TF1 for this pdg
          rbw = iter->second.get();
	}
	else {
	  // otherwise create it
	  fPDGtoTF1[pdg] = std::make_unique<TF1>("rbw",
                                                 "pow([1],2)*pow([0],2)/(pow(x*x-[0]*[0],2)+pow(x*x*[1]/[0],2))", am - 5 * awidth, am + 5 * awidth);
	  fPDGtoTF1[pdg]->SetParameter(0, am);
	  fPDGtoTF1[pdg]->SetParameter(1, awidth);
	  rbw = fPDGtoTF1[pdg].get();
	}
        am = rbw->GetRandom();
      }
      // -----------------------------------------------//

      //
      // y
      ty = TMath::TanH(fYPara->GetRandom());
      //
      // pT
      if (fAnalog == kAnalog) {
        pt = fPtPara->GetRandom();
        wgtp = fParentWeight;
        wgtch = fChildWeight;
      } else {
        pt=fPtMin+random[1]*(fPtMax-fPtMin);
        Double_t ptd=pt;
        wgtp=fParentWeight*fPtParaFunc(& ptd, &dummy);
        wgtch=fChildWeight*fPtParaFunc(& ptd, &dummy);
      }
      xmt = sqrt(pt * pt + am * am);
      if (TMath::Abs(ty) == 1.) {
        ty = 0.;
        Fatal("AliGenParam",
              "Division by 0: Please check you rapidity range !");
      }
      //
      // phi

      // set the right parameters for fdNdPhi
      double v2 = fV2Para->Eval(pt);
      if (fdNdPhi->GetParameter(0) != v2) {
	// we should avoid calling SetParam as this invalidates the
	// internal integral cache of TF1
        fdNdPhi->SetParameter(0, v2);
      }
      if (fdNdPhi->GetParameter(1) != fEvPlane) {
        fdNdPhi->SetParameter(1, fEvPlane);
      }
      phi = fdNdPhi->GetRandom();
      pl = xmt * ty / sqrt((1. - ty) * (1. + ty));
      theta = TMath::ATan2(pt, pl);
      // Cut on theta
      if (theta < fThetaMin || theta > fThetaMax)
        continue;
      ptot = TMath::Sqrt(pt * pt + pl * pl);
      // Cut on momentum
      if (ptot < fPMin || ptot > fPMax)
        continue;
      //
      p[0] = pt * TMath::Cos(phi);
      p[1] = pt * TMath::Sin(phi);
      p[2] = pl;
      energy = TMath::Sqrt(ptot * ptot + am * am);

      // if fForceDecay != none Primary particle decays using
      // AliPythia and children are tracked by GEANT
      //
      // if fForceDecay == none Primary particle is tracked by GEANT
      // (In the latest, make sure that GEANT actually does all the decays you
      // want)
      //
      Bool_t decayed = kFALSE;

      if (fForceDecay != kNoDecay) {
        // Using lujet to decay particle
        TLorentzVector pmom(p[0], p[1], p[2], energy);
        fDecayer->Decay(pdg, &pmom);
        //
        // select decay particles
        Int_t np = fDecayer->ImportParticles(particles);
        pdg = iTemp;
        if (pdg >= 220000 & pdg <= 220001) {
          TParticle *gamma = (TParticle *)particles->At(0);
          gamma->SetPdgCode(pdg);
          np = VirtualGammaPairProduction(particles, np);
        }
        if (fForceConv)
          np = ForceGammaConversion(particles, np);

        auto ncsel = 0;
        vFlags.reserve(np);
        vFlags = {true};
        vSelected.reserve(np);
        vSelected = {true};
        vParent.reserve(np);

        if (np > 1) {
          decayed = kTRUE;
          TParticle *iparticle = 0;
          Int_t ipF, ipL;
          for (i = 1; i < np; i++) {

            iparticle = (TParticle *)particles->At(i);
            Int_t kf = iparticle->GetPdgCode();
            Int_t ks = iparticle->GetStatusCode();
            // flagged particle
            if (!fPreserveFullDecayChain) {
              if (vFlags[i]) {
                ipF = iparticle->GetFirstDaughter() + fIncFortran;
                ipL = iparticle->GetLastDaughter() + fIncFortran;
                if (ipF > 0)
                  for (j = ipF; j <= ipL; j++)
                    vFlags[j] = true;
                continue;
              }
            }
            // flag decay products of particles with long life-time (ctau > .3
            // mum)
            if (ks != 1) {
              Double_t lifeTime = fDecayer->GetLifetime(kf);
              if (lifeTime > (Double_t)fMaxLifeTime) {
                ipF = iparticle->GetFirstDaughter() + fIncFortran;
                ipL = iparticle->GetLastDaughter() + fIncFortran;
                if (ipF > 0) {
                  for (j = ipF; j <= ipL; j++)
                    vFlags[j] = 1;
                } else {
                  vSelected[i] = true;
                }
              }
            } // ks==1 ?
              //
              // children

            if ((/* ChildSelected(TMath::Abs(kf)) ||*/ fForceDecay !=
                     kNoDecay ||
                 fSelectAll)) {

              if (fCutOnChild) {
                pc[0] = iparticle->Px();
                pc[1] = iparticle->Py();
                pc[2] = iparticle->Pz();
                Bool_t childok = KinematicSelection(iparticle, 1);
                if (childok) {
                  vSelected[i] = true;
                  ncsel++;
                } else {
                  if (!fKeepIfOneChildSelected) {
                    ncsel = -1;
                    break;
                  }
                } // child kine cuts
              } else {
                vSelected[i] = true;
                ncsel++;
              } // if child selection
            }   // child selection
          }     // decay particle loop
        }       // if decay products

        Int_t iparent;

        if (fKeepParent || (fCutOnChild && ncsel > 0) || !fCutOnChild) {
          //
          // Parent
          // --- For Exodus --------------------------------//
    auto particle = new TParticle(
          pdg, ((decayed) ? 11 : 1), -1, -1, -1, -1, p[0], p[1], p[2],
          energy, origin0[0], origin0[1], origin0[2], time0
          );
    particle->SetWeight(wgtp);
    fParticles->Add(particle);
          vParent[0] = nt;
          nt++;
          fNprimaries++;

          // but count is as "generated" particle" only if it produced child(s)
          // within cut
          if ((fCutOnChild && ncsel > 0) || !fCutOnChild) {
            ipa++;
          }

          //
          // Decay Products
          //
          for (i = 1; i < np; i++) {
            if (vSelected[i]) {
              TParticle *iparticle = (TParticle *)particles->At(i);
              auto kf = iparticle->GetPdgCode();
              auto ksc = iparticle->GetStatusCode();
              auto jpa = iparticle->GetFirstMother() + fIncFortran;
        Double_t weight = iparticle->GetWeight();
              och[0] = origin0[0] + iparticle->Vx();
              och[1] = origin0[1] + iparticle->Vy();
              och[2] = origin0[2] + iparticle->Vz();
              pc[0] = iparticle->Px();
              pc[1] = iparticle->Py();
              pc[2] = iparticle->Pz();
              Double_t ec = iparticle->Energy();

              if (jpa > -1) {
                iparent = vParent[jpa];
              } else {
                iparent = -1;
              }
              auto parentP = (TParticle *)fParticles->At(iparent);
              if (parentP->GetFirstDaughter() == -1)
                parentP->SetFirstDaughter(nt);
              parentP->SetLastDaughter(nt);
        auto particle = new TParticle(kf, ksc, iparent, -1, -1, -1, pc[0],
                                            pc[1], pc[2], ec, och0[0], och0[1],
                                            och0[2], time0 + iparticle->T()
              );
        particle->SetWeight(weight * wgtch);
              fParticles->Add(particle);

              vParent[i] = nt;
              nt++;
              fNprimaries++;
            } // Selected
          }   // Particle loop
        }     // Decays by Lujet
        particles->Clear();
        vFlags.clear();
        vParent.clear();
        vSelected.clear();
        // kinematic selection
      } else {
        // nodecay option, so parent will be tracked by GEANT (pions, kaons,
        // eta, omegas, baryons)
  auto particle = new TParticle(pdg, 1, -1, -1, -1, -1, p[0], p[1], p[2],
                                      energy, origin0[0], origin0[1],
                                      origin0[2], time0);
  particle->SetWeight(wgtp);
        fParticles->Add(particle);
        ipa++;
        fNprimaries++;
      }
      break;
    } // while
  }   // event loop
}

int GeneratorParam::ImportParticles(TClonesArray *particles, Option_t *option) {
  if (particles == 0)
    return 0;
  TClonesArray &clonesParticles = *particles;
  clonesParticles.Clear();
  Int_t numpart = fParticles->GetEntries();
  for (int i = 0; i < numpart; i++) {
    TParticle *particle = (TParticle *)fParticles->At(i);
    new (clonesParticles[i]) TParticle(*particle);
  }
  return numpart;
}

void GeneratorParam::InitChildSelect() {
  //  Initialize child selection
  fChildSelect.Set(5);
  for (Int_t i = 0; i < 5; i++)
    fChildSelect[i] = 0;
  switch (fForceDecay) {
  case kBSemiElectronic:
  case kSemiElectronic:
  case kDiElectron:
  case kBJpsiDiElectron:
  case kBPsiPrimeDiElectron:
  case kElectronEM:
  case kDiElectronEM:
    fChildSelect[0] = kElectron;
    break;
  case kHardMuons:
  case kBSemiMuonic:
  case kDSemiMuonic:
  case kSemiMuonic:
  case kDiMuon:
  case kJpsiDiMuon:
  case kBJpsiDiMuon:
  case kBPsiPrimeDiMuon:
  case kPiToMu:
  case kKaToMu:
  case kWToMuon:
  case kWToCharmToMuon:
  case kZDiMuon:
  case kZDiElectron:
    fChildSelect[0] = kMuonMinus;
    break;
  case kHadronicD:
  case kHadronicDWithout4Bodies:
  case kHadronicDWithV0:
  case kHadronicDWithout4BodiesWithV0:
    fChildSelect[0] = kPiPlus;
    fChildSelect[1] = kKPlus;
    break;
  case kPhiKK:
    fChildSelect[0] = kKPlus;
    break;
  case kBJpsi:
  case kBJpsiUndecayed:
    fChildSelect[0] = 443;
    break;
  case kChiToJpsiGammaToMuonMuon:
    fChildSelect[0] = 22;
    fChildSelect[1] = 13;
    break;
  case kChiToJpsiGammaToElectronElectron:
    fChildSelect[0] = 22;
    fChildSelect[1] = 11;
    break;
  case kLambda:
    fChildSelect[0] = kProton;
    fChildSelect[1] = 211;
    break;
  case kPsiPrimeJpsiDiElectron:
    fChildSelect[0] = 211;
    fChildSelect[1] = 11;
    break;
  case kGammaEM:
    fChildSelect[0] = kGamma;
    break;
  case kLcpKpi:
    fChildSelect[0] = kProton;
    fChildSelect[1] = kKPlus;
    fChildSelect[2] = kPiPlus;
    break;
  case kLcpK0S:
    fChildSelect[0] = kProton;
    fChildSelect[1] = kK0;
    break;
  default:
    break;
  }
}

Bool_t GeneratorParam::KinematicSelection(const TParticle *particle,
                                          Int_t flag) const {
  // Perform kinematic selection
  Double_t pz = particle->Pz();
  Double_t pt = particle->Pt();
  Double_t p = particle->P();
  Double_t theta = particle->Theta();
  Double_t mass = particle->GetCalcMass();
  Double_t mt2 = pt * pt + mass * mass;
  Double_t phi = particle->Phi();
  Double_t e = particle->Energy();

  if (e == 0.)
    e = TMath::Sqrt(p * p + mass * mass);

  Double_t y, y0;

  if (TMath::Abs(pz) < e) {
    y = 0.5 * TMath::Log((e + pz) / (e - pz));
  } else {
    y = 1.e10;
  }

  if (mt2) {
    y0 = 0.5 * TMath::Log((e + TMath::Abs(pz)) * (e + TMath::Abs(pz)) / mt2);
  } else {
    if (TMath::Abs(y) < 1.e10) {
      y0 = y;
    } else {
      y0 = 1.e10;
    }
  }

  y = (pz < 0) ? -y0 : y0;

  if (flag == 0) {
    //
    // Primary particle cuts
    //
    //  transverse momentum cut
    if (pt > fPtMax || pt < fPtMin) {
      return kFALSE;
    }
    //
    // momentum cut
    if (p > fPMax || p < fPMin) {
      return kFALSE;
    }
    //
    // theta cut
    if (theta > fThetaMax || theta < fThetaMin) {
      return kFALSE;
    }
    //
    // rapidity cut
    if (y > fYMax || y < fYMin) {
      return kFALSE;
    }
    //
    // phi cut
    if (phi > fPhiMax || phi < fPhiMin) {
      return kFALSE;
    }
  } else {
    //
    // Decay product cuts
    //
    //  transverse momentum cut
    if (pt > fChildPtMax || pt < fChildPtMin) {
      return kFALSE;
    }
    //
    // momentum cut
    if (p > fChildPMax || p < fChildPMin) {
      return kFALSE;
    }
    //
    // theta cut
    if (theta > fChildThetaMax || theta < fChildThetaMin) {
      return kFALSE;
    }
    //
    // rapidity cut
    if (y > fChildYMax || y < fChildYMin) {
      return kFALSE;
    }
    //
    // phi cut
    if (phi > fChildPhiMax || phi < fChildPhiMin) {
      return kFALSE;
    }
  }

  return kTRUE;
}

//_______________________________________________________________________
void GeneratorParam::SetPtRange(Float_t ptmin, Float_t ptmax) {
  //
  // Set the Pt range for the generated particles
  //
  fPtMin = ptmin;
  fPtMax = ptmax;
  SetBit(kPtRange);
}
//_______________________________________________________________________
void GeneratorParam::SetMomentumRange(Float_t pmin, Float_t pmax) {
  //
  // Set the momentum range for the generated particles
  //
  fPMin = pmin;
  fPMax = pmax;
  SetBit(kMomentumRange);
}

//_______________________________________________________________________
void GeneratorParam::SetPhiRange(Float_t phimin, Float_t phimax) {
  //
  // Set the Phi range for the generated particles
  //
  fPhiMin = TMath::Pi() * phimin / 180;
  fPhiMax = TMath::Pi() * phimax / 180;
  SetBit(kPhiRange);
}

//_______________________________________________________________________
void GeneratorParam::SetYRange(Float_t ymin, Float_t ymax) {
  //
  // Set the Rapidity range for the generated particles
  //
  fYMin = ymin;
  fYMax = ymax;
  SetBit(kYRange);
}

//_______________________________________________________________________
void GeneratorParam::SetThetaRange(Float_t thetamin, Float_t thetamax) {
  //
  // Set the theta range for the generated particles
  //
  fThetaMin = TMath::Pi() * thetamin / 180;
  fThetaMax = TMath::Pi() * thetamax / 180;
  SetBit(kThetaRange);
}

//____________________________________________________________________________________
Float_t GeneratorParam::GetRelativeArea(Float_t ptMin, Float_t ptMax,
                                        Float_t yMin, Float_t yMax,
                                        Float_t phiMin, Float_t phiMax) {
  //
  // Normalisation for selected kinematic region
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(5, 99, 0)
  Float_t ratio = fPtPara->Integral(ptMin, ptMax, (Double_t *)0, 1.e-6) /
                  fPtPara->Integral(fPtPara->GetXmin(), fPtPara->GetXmax(),
                                    (Double_t *)0, 1.e-6) *
                  fYPara->Integral(yMin, yMax, (Double_t *)0, 1.e-6) /
                  fYPara->Integral(fYPara->GetXmin(), fYPara->GetXmax(),
                                   (Double_t *)0, 1.e-6) *
                  (phiMax - phiMin) / 360.;
#else
  Float_t ratio =
      fPtPara->Integral(ptMin, ptMax, 1.e-6) /
      fPtPara->Integral(fPtPara->GetXmin(), fPtPara->GetXmax(), 1.e-6) *
      fYPara->Integral(yMin, yMax, 1.e-6) /
      fYPara->Integral(fYPara->GetXmin(), fYPara->GetXmax(), 1.e-6) *
      (phiMax - phiMin) / 360.;
#endif
  return TMath::Abs(ratio);
}

//____________________________________________________________________________________

void GeneratorParam::Draw(const char * /*opt*/) {
  //
  // Draw the pT and y Distributions
  //
  TCanvas *c0 = new TCanvas("c0", "Canvas 0", 400, 10, 600, 700);
  c0->Divide(2, 1);
  c0->cd(1);
  fPtPara->Draw();
  fPtPara->GetHistogram()->SetXTitle("p_{T} (GeV)");
  c0->cd(2);
  fYPara->Draw();
  fYPara->GetHistogram()->SetXTitle("y");
}

//-------------------------------------------------------------------
TVector3 GeneratorParam::OrthogonalVector(TVector3 &inVec) {
  double abc[] = {inVec.x(), inVec.y(), inVec.z()};
  double xyz[] = {1, 1, 1};
  int solvDim = 0;
  double tmp = abc[0];
  for (int i = 0; i < 3; i++)
    if (fabs(abc[i]) > tmp) {
      solvDim = i;
      tmp = fabs(abc[i]);
    }
  xyz[solvDim] = (-abc[(1 + solvDim) % 3] - abc[(2 + solvDim) % 3]) /
                 abc[(0 + solvDim) % 3];

  TVector3 res(xyz[0], xyz[1], xyz[2]);
  return res;
}

void GeneratorParam::RotateVector(Double_t *pin, Double_t *pout,
                                  Double_t costheta, Double_t sintheta,
                                  Double_t cosphi, Double_t sinphi) {
  // Perform rotation
  pout[0] =
      pin[0] * costheta * cosphi - pin[1] * sinphi + pin[2] * sintheta * cosphi;
  pout[1] =
      pin[0] * costheta * sinphi + pin[1] * cosphi + pin[2] * sintheta * sinphi;
  pout[2] = -1.0 * pin[0] * sintheta + pin[2] * costheta;
  return;
}

double GeneratorParam::ScreenFunction1(double screenVariable) {
  if (screenVariable > 1)
    return 42.24 - 8.368 * log(screenVariable + 0.952);
  else
    return 42.392 - screenVariable * (7.796 - 1.961 * screenVariable);
}

double GeneratorParam::ScreenFunction2(double screenVariable) {
  if (screenVariable > 1)
    return 42.24 - 8.368 * log(screenVariable + 0.952);
  else
    return 41.405 - screenVariable * (5.828 - 0.8945 * screenVariable);
}

double GeneratorParam::RandomEnergyFraction(double Z, double photonEnergy) {
  double aZ = Z / 137.036;
  double epsilon;
  double epsilon0Local = 0.000511 / photonEnergy;

  // Do it fast if photon energy < 2. MeV
  if (photonEnergy < 0.002) {
    epsilon = epsilon0Local + (0.5 - epsilon0Local) * gRandom->Rndm();
  } else {
    double fZ = 8 * log(Z) / 3;
    double fcZ = (aZ * aZ) * (1 / (1 + aZ * aZ) + 0.20206 - 0.0368 * aZ * aZ +
                              0.0083 * aZ * aZ * aZ);
    if (photonEnergy > 0.050)
      fZ += 8 * fcZ;

    // Limits of the screening variable
    double screenFactor = 136. * epsilon0Local / std::cbrt(Z);
    double screenMax = exp((42.24 - fZ) / 8.368) - 0.952;
    double screenMin = std::min(4. * screenFactor, screenMax);

    // Limits of the energy sampling
    double epsilon1 = 0.5 - 0.5 * sqrt(1. - screenMin / screenMax);
    double epsilonMin = std::max(epsilon0Local, epsilon1);
    double epsilonRange = 0.5 - epsilonMin;

    // Sample the energy rate of the created electron (or positron)
    double screen;
    double gReject;

    double f10 = ScreenFunction1(screenMin) - fZ;
    double f20 = ScreenFunction2(screenMin) - fZ;
    double normF1 = std::max(f10 * epsilonRange * epsilonRange, 0.);
    double normF2 = std::max(1.5 * f20, 0.);

    do {
      if (normF1 / (normF1 + normF2) > gRandom->Rndm()) {
        epsilon = 0.5 - epsilonRange * std::cbrt(gRandom->Rndm());
        screen = screenFactor / (epsilon * (1. - epsilon));
        gReject = (ScreenFunction1(screen) - fZ) / f10;
      } else {
        epsilon = epsilonMin + epsilonRange * gRandom->Rndm();
        screen = screenFactor / (epsilon * (1 - epsilon));
        gReject = (ScreenFunction2(screen) - fZ) / f20;
      }
    } while (gReject < gRandom->Rndm());
  } //  End of epsilon sampling
  return epsilon;
}

double GeneratorParam::RandomPolarAngle() {
  double u;
  const double a1 = 0.625;
  double a2 = 3. * a1;

  if (0.25 > gRandom->Rndm()) {
    u = -log(gRandom->Rndm() * gRandom->Rndm()) / a1;
  } else {
    u = -log(gRandom->Rndm() * gRandom->Rndm()) / a2;
  }
  return u * 0.000511;
}

Double_t GeneratorParam::RandomMass(Double_t mh) {
  while (true) {
    double y = gRandom->Rndm();
    double mee =
        2 * 0.000511 *
        TMath::Power(2 * 0.000511 / mh,
                     -y); // inverse of the enveloping cumulative distribution
    double apxkw = 2.0 / 3.0 / 137.036 / TMath::Pi() /
                   mee; // enveloping probability density
    double val = gRandom->Uniform(0, apxkw);
    double kw = apxkw * sqrt(1 - 4 * 0.000511 * 0.000511 / mee / mee) *
                (1 + 2 * 0.000511 * 0.000511 / mee / mee) * 1 * 1 *
                TMath::Power(1 - mee * mee / mh / mh, 3);
    if (val < kw)
      return mee;
  }
}

Int_t GeneratorParam::VirtualGammaPairProduction(TClonesArray *particles,
                                                 Int_t nPart) {
  Int_t nPartNew = nPart;
  for (int iPart = 0; iPart < nPart; iPart++) {
    TParticle *gamma = (TParticle *)particles->At(iPart);
    if (gamma->GetPdgCode() != 220001)
      continue;
    if (gamma->Pt() < 0.002941)
      continue; // approximation of kw in AliGenEMlib is 0 below 0.002941
    double mass = RandomMass(gamma->Pt());

    // lepton pair kinematics in virtual photon rest frame
    double Ee = mass / 2;
    double Pe = TMath::Sqrt((Ee + 0.000511) * (Ee - 0.000511));

    double costheta = (2.0 * gRandom->Rndm()) - 1.;
    double sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
    double phi = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
    double sinphi = TMath::Sin(phi);
    double cosphi = TMath::Cos(phi);

    // momentum vectors of leptons in virtual photon rest frame
    Double_t pProd1[3] = {Pe * sintheta * cosphi, Pe * sintheta * sinphi,
                          Pe * costheta};

    Double_t pProd2[3] = {-1.0 * Pe * sintheta * cosphi,
                          -1.0 * Pe * sintheta * sinphi, -1.0 * Pe * costheta};

    // lepton 4-vectors in properly rotated virtual photon rest frame
    Double_t pRot1[3] = {0.};
    RotateVector(pProd1, pRot1, costheta, -sintheta, -cosphi, -sinphi);
    Double_t pRot2[3] = {0.};
    RotateVector(pProd2, pRot2, costheta, -sintheta, -cosphi, -sinphi);

    TLorentzVector e1V4(pRot1[0], pRot1[1], pRot1[2], Ee);
    TLorentzVector e2V4(pRot2[0], pRot2[1], pRot2[2], Ee);

    TVector3 boost(gamma->Px(), gamma->Py(), gamma->Pz());
    boost *= 1 / sqrt(gamma->P() * gamma->P() + mass * mass);
    e1V4.Boost(boost);
    e2V4.Boost(boost);

    TLorentzVector vtx;
    gamma->ProductionVertex(vtx);
    new ((*particles)[nPartNew])
        TParticle(11, gamma->GetStatusCode(), iPart + 1, -1, 0, 0, e1V4, vtx);
    nPartNew++;
    new ((*particles)[nPartNew])
        TParticle(-11, gamma->GetStatusCode(), iPart + 1, -1, 0, 0, e2V4, vtx);
    nPartNew++;
  }
  return nPartNew;
}

Int_t GeneratorParam::ForceGammaConversion(TClonesArray *particles,
                                           Int_t nPart) {
  // based on:
  // http://geant4.cern.ch/G4UsersDocuments/UsersGuides/PhysicsReferenceManual/html/node27.html
  //     and:
  //     http://geant4.cern.ch/G4UsersDocuments/UsersGuides/PhysicsReferenceManual/html/node58.html
  //     and: G4LivermoreGammaConversionModel.cc
  Int_t nPartNew = nPart;
  for (int iPart = 0; iPart < nPart; iPart++) {
    TParticle *gamma = (TParticle *)particles->At(iPart);
    if (gamma->GetPdgCode() != 22 & gamma->GetPdgCode() != 220000)
      continue;
    if (gamma->Energy() <= 0.001022)
      continue;
    TVector3 gammaV3(gamma->Px(), gamma->Py(), gamma->Pz());
    double frac = RandomEnergyFraction(1, gamma->Energy());
    double Ee1 = frac * gamma->Energy();
    double Ee2 = (1 - frac) * gamma->Energy();
    double Pe1 = sqrt((Ee1 + 0.000511) * (Ee1 - 0.000511));
    double Pe2 = sqrt((Ee2 + 0.000511) * (Ee2 - 0.000511));

    TVector3 rotAxis(OrthogonalVector(gammaV3));
    Float_t az = gRandom->Uniform(TMath::Pi() * 2);
    rotAxis.Rotate(az, gammaV3);
    TVector3 e1V3(gammaV3);
    double u = RandomPolarAngle();
    e1V3.Rotate(u / Ee1, rotAxis);
    e1V3 = e1V3.Unit();
    e1V3 *= Pe1;
    TVector3 e2V3(gammaV3);
    e2V3.Rotate(-u / Ee2, rotAxis);
    e2V3 = e2V3.Unit();
    e2V3 *= Pe2;
    // gamma = new TParticle(*gamma);
    // particles->RemoveAt(iPart);
    gamma->SetFirstDaughter(nPartNew + 1);
    gamma->SetLastDaughter(nPartNew + 2);
    // new((*particles)[iPart]) TParticle(*gamma);
    // delete gamma;

    // conversion probability per atom
    // fitted G4EMLOW6.35/pair/pp-cs-8.dat, fit is great for E>20MeV
    double convProb = 1 / (exp(-log(28.44 * (gamma->Energy() - 0.001022)) *
                               (0.775 + 0.0271 * log(gamma->Energy() + 1))) +
                           1);

    // radiation length is not considered here, so you have to normalize
    // yourself in after-production from infinite radiation length to whatever
    // you want double meanExessPathlength=0.5*(exp(0.9)-exp(-0.9))/0.9; double
    // scale=(1-exp(-7.0/9.0*radLength*meanExessPathlength))/(1-0);

    TLorentzVector vtx;
    gamma->ProductionVertex(vtx);
    TParticle *currPart;
    Int_t sign = (gRandom->Rndm() < 0.5) ? 1 : -1;
    currPart = new ((*particles)[nPartNew])
        TParticle(sign * 220011, gamma->GetStatusCode(), iPart + 1, -1, 0, 0,
                  TLorentzVector(e1V3, Ee1), vtx);
    currPart->SetWeight(convProb);
    nPartNew++;
    currPart = new ((*particles)[nPartNew])
        TParticle(-sign * 220011, gamma->GetStatusCode(), iPart + 1, -1, 0, 0,
                  TLorentzVector(e2V3, Ee2), vtx);
    currPart->SetWeight(convProb);
    nPartNew++;
  }
  return nPartNew;
}
