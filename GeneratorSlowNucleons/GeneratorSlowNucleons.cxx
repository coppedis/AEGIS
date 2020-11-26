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
//  Generator for slow nucleons in pA interactions.
//  Source is modelled by a relativistic Maxwell distributions.
//  In this case the number of slow nucleons is determined from the number of
//  wounded nuclei using a realisation of AliSlowNucleonModel. Original code by
//  Ferenc Sikler  <sikler@rmki.kfki.hu>
//

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TRandom.h>

#include "GeneratorSlowNucleons.h"
#include "SlowNucleonModel.h"

ClassImp(GeneratorSlowNucleons)

    GeneratorSlowNucleons::GeneratorSlowNucleons()
    : TGenerator("GeneratorParam", "GeneratorParam"), fCMS(0.), fMomentum(0.),
      fBeta(0.), fPmax(0.), fCharge(0), fProtonDirection(1.), fTemperatureG(0.),
      fBetaSourceG(0.), fTemperatureB(0.), fBetaSourceB(0.), fNgp(0), fNgn(0),
      fNbp(0), fNbn(0), fDebug(0), fDebugHist1(0), fDebugHist2(0),
      fThetaDistribution(), fCosThetaGrayHist(), fCosTheta(),
      fBeamCrossingAngle(0.), fBeamDivergence(0.), fBeamDivEvent(0.),
      fSmearMode(2), fNcoll(0), fSlowNucleonModel(0) {
  // Default constructor
}

GeneratorSlowNucleons::GeneratorSlowNucleons(Int_t npart)
    : TGenerator("GeneratorParam", "GeneratorParam"), fCMS(14000.),
      fMomentum(0.), fBeta(0.), fPmax(10.), fCharge(1), fProtonDirection(1.),
      fTemperatureG(0.05), fBetaSourceG(0.05), fTemperatureB(0.005),
      fBetaSourceB(0.), fNgp(0), fNgn(0), fNbp(0), fNbn(0), fDebug(0),
      fDebugHist1(0), fDebugHist2(0), fThetaDistribution(), fCosThetaGrayHist(),
      fCosTheta(), fBeamCrossingAngle(0.), fBeamDivergence(0.),
      fBeamDivEvent(0.), fSmearMode(2), fNcoll(0),
      fSlowNucleonModel(new SlowNucleonModel())

{
  // Constructor
  fName = "SlowNucleons";
  fTitle = "Generator for gray particles in pA collisions";
}

//____________________________________________________________
GeneratorSlowNucleons::~GeneratorSlowNucleons() {
  // Destructor
  delete fSlowNucleonModel;
}

void GeneratorSlowNucleons::SetProtonDirection(Float_t dir) {
  // Set direction of the proton to change between pA (1) and Ap (-1)
  fProtonDirection = dir / TMath::Abs(dir);
}

void GeneratorSlowNucleons::Init() {
  //
  // Initialization
  //
  Double_t kMass = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
  fMomentum = fCMS / 2. * Float_t(fZTarget) / Float_t(fATarget);
  fBeta = fMomentum / TMath::Sqrt(kMass * kMass + fMomentum * fMomentum);
  printf("  fMomentum %f    fBeta %1.10f\n", fMomentum, fBeta);
  if (fDebug) {
    fDebugHist1 =
        new TH2F("DebugHist1", "nu vs N_slow", 100, 0., 100., 20, 0., 20.);
    fDebugHist2 =
        new TH2F("DebugHist2", "b  vs N_slow", 100, 0., 100., 15, 0., 15.);
    fCosThetaGrayHist =
        new TH1F("fCosThetaGrayHist", "Gray particles angle", 100, -1., 1.);
  }
  //
  // non-uniform cos(theta) distribution
  //
  if (fThetaDistribution != 0) {
    fCosTheta =
        new TF1("fCosTheta",
                "(2./3.14159265358979312)/(exp(2./3.14159265358979312)-exp(-2./"
                "3.14159265358979312))*exp(2.*x/3.14159265358979312)",
                -1., 1.);
  }

  if (TMath::Abs(fBeamCrossingAngle) > 0.)
    printf("\n  GeneratorSlowNucleons: applying crossing angle %f mrad to slow "
           "nucleons\n",
           fBeamCrossingAngle * 1000.);
}

void GeneratorSlowNucleons::FinishRun() {
  // End of run action
  // Show histogram for debugging if requested.
  if (fDebug) {
    TCanvas *c = new TCanvas("c", "Canvas 1", 400, 10, 600, 700);
    c->Divide(2, 1);
    c->cd(1);
    fDebugHist1->Draw("colz");
    c->cd(2);
    fDebugHist2->Draw();
    c->cd(3);
    fCosThetaGrayHist->Draw();
  }
}

void GeneratorSlowNucleons::GenerateEvent() {
  //
  // Generate one event
  //
  const Float_t mp = (TDatabasePDG::Instance()->GetParticle(kProton))->Mass();
  const Float_t mn = (TDatabasePDG::Instance()->GetParticle(kNeutron))->Mass();
  const Float_t k2PI = 2. * TMath::Pi();
  const Float_t kRaddeg = 180. / TMath::Pi();
  const Float_t kDegrad = TMath::Pi() / 180.;
  fParticles->Clear();
  //
  // Communication with Gray Particle Model
  //
  //  if (fCollisionGeometry) {
  //	Float_t b   = fCollisionGeometry->b;

  // (1) Sikler' model
  if (fSmearMode == 0)
    fSlowNucleonModel->GetNumberOfSlowNucleons(fNcoll, fNgp, fNgn, fNbp, fNbn);
  // (2) Model inspired on exp. data at lower energy (Gallio-Oppedisano)
  // --- smearing the Ncoll fron generator used as input
  else if (fSmearMode == 1)
    fSlowNucleonModel->GetNumberOfSlowNucleons2(fNcoll, fNgp, fNgn, fNbp, fNbn);
  // --- smearing directly Nslow
  else if (fSmearMode == 2)
    fSlowNucleonModel->GetNumberOfSlowNucleons2s(fNcoll, fNgp, fNgn, fNbp,
                                                 fNbn);
  if (fDebug) {
    // printf("Collision Geometry %f %d %d %d %d\n", b, nn, nwn, nnw, nwnw);
    printf("Slow nucleons: %d grayp  %d grayn  %d blackp  %d blackn \n", fNgp,
           fNgn, fNbp, fNbn);
    fDebugHist1->Fill(Float_t(fNgp + fNgn + fNbp + fNbn), fNcoll, 1.);
    //	    fDebugHist2->Fill(Float_t(fNgp + fNgn + fNbp + fNbn), b, 1.);
  }
  //
  Float_t p[3] = {0., 0., 0.}, theta = 0;
  Float_t polar[3] = {0., 0., 0.};
  Int_t nt, i, j;
  Int_t kf;
  Double_t energy;

  // Extracting 1 value per event for the divergence angle
  Double_t rvec = gRandom->Gaus(0.0, 1.0);
  fBeamDivEvent = fBeamDivergence * TMath::Abs(rvec);
  if (TMath::Abs(fBeamDivEvent) > 0.)
    printf("\n  GeneratorSlowNucleons: applying beam divergence %f mrad to "
           "slow nucleons\n",
           fBeamDivEvent * 1000.);
  //
  //  Gray protons
  //
  fCharge = 1;
  kf = kProton;
  for (i = 0; i < fNgp; i++) {
    GenerateSlow(fCharge, fTemperatureG, fBetaSourceG, p, theta);
    if (fDebug)
      fCosThetaGrayHist->Fill(TMath::Cos(theta));
    energy = TMath::Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + mp * mp);
    auto part = new TParticle(kf, 1, -1, -1, -1, -1, p[0], p[1], p[2], energy,
                              0., 0., 0., 0.);
    part->SetUniqueID(kGrayProcess);
    fParticles->Add(part);
  }
  //
  //  Gray neutrons
  //
  fCharge = 0;
  kf = kNeutron;
  for (i = 0; i < fNgn; i++) {
    GenerateSlow(fCharge, fTemperatureG, fBetaSourceG, p, theta);
    if (fDebug)
      fCosThetaGrayHist->Fill(TMath::Cos(theta));
    energy = TMath::Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + mn * mn);
    auto part = new TParticle(kf, 1, -1, -1, -1, -1, p[0], p[1], p[2], energy,
                              0., 0., 0., 0.);
    part->SetUniqueID(kGrayProcess);
    fParticles->Add(part);
  }
  //
  //  Black protons
  //
  fCharge = 1;
  kf = kProton;
  for (i = 0; i < fNbp; i++) {
    GenerateSlow(fCharge, fTemperatureB, fBetaSourceB, p, theta);
    energy = TMath::Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + mp * mp);
    auto part = new TParticle(kf, 1, -1, -1, -1, -1, p[0], p[1], p[2], energy,
                              0., 0., 0., 0.);
    part->SetUniqueID(kBlackProcess);
    fParticles->Add(part);
  }
  //
  //  Black neutrons
  //
  fCharge = 0;
  kf = kNeutron;
  for (i = 0; i < fNbn; i++) {
    GenerateSlow(fCharge, fTemperatureB, fBetaSourceB, p, theta);
    energy = TMath::Sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] + mn * mn);
    auto part = new TParticle(kf, 1, -1, -1, -1, -1, p[0], p[1], p[2], energy,
                              0., 0., 0., 0.);
    part->SetUniqueID(kBlackProcess);
    fParticles->Add(part);
  }
}

int GeneratorSlowNucleons::ImportParticles(TClonesArray *particles,
                                           Option_t *option) {
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

void GeneratorSlowNucleons::GenerateSlow(Int_t charge, Double_t T,
                                         Double_t beta, Float_t *q,
                                         Float_t &theta)

{
  /*
     Emit a slow nucleon with "temperature" T [GeV],
     from a source moving with velocity beta
     Three-momentum [GeV/c] is given back in q[3]
  */

  // printf("Generating slow nuc. with: charge %d. temp. %1.4f, beta %f
  // \n",charge,T,beta);

  Double_t m, pmax, p, f, phi;
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  const Double_t kMassProton = pdg->GetParticle(kProton)->Mass();
  const Double_t kMassNeutron = pdg->GetParticle(kNeutron)->Mass();

  /* Select nucleon type */
  if (charge == 0)
    m = kMassNeutron;
  else
    m = kMassProton;

  /* Momentum at maximum of Maxwell-distribution */

  pmax = TMath::Sqrt(2. * T * (T + TMath::Sqrt(T * T + m * m)));

  /* Try until proper momentum                                  */
  /* for lack of primitive function of the Maxwell-distribution */
  /* a brute force trial-accept loop, normalized at pmax        */

  do {
    p = gRandom->Rndm() * fPmax;
    f = Maxwell(m, p, T) / Maxwell(m, pmax, T);
  } while (f < gRandom->Rndm());

  /* Spherical symmetric emission for black particles (beta=0)*/
  if (beta == 0 || fThetaDistribution == 0)
    theta = TMath::ACos(2. * gRandom->Rndm() - 1.);
  /* cos theta distributed according to experimental results for gray particles
   * (beta=0.05)*/
  else if (fThetaDistribution != 0) {
    Double_t costheta = fCosTheta->GetRandom();
    theta = TMath::ACos(costheta);
  }
  //
  phi = 2. * TMath::Pi() * gRandom->Rndm();

  /* Determine momentum components in system of the moving source */
  q[0] = p * TMath::Sin(theta) * TMath::Cos(phi);
  q[1] = p * TMath::Sin(theta) * TMath::Sin(phi);
  q[2] = p * TMath::Cos(theta);
  if (fDebug == 1)
    printf("\n Momentum in RS of the moving source: p = (%f, %f, %f)\n", q[0],
           q[1], q[2]);

  /* Transform to system of the target nucleus                             */
  /* beta is passed as negative, because the gray nucleons are slowed down */
  Lorentz(m, -beta, q);
  if (fDebug == 1)
    printf(" Momentum in RS of the target nucleus: p = (%f, %f, %f)\n", q[0],
           q[1], q[2]);

  /* Transform to laboratory system */
  Lorentz(m, fBeta, q);
  q[2] *= fProtonDirection;
  if (fDebug == 1)
    printf("\n Momentum after LHC boost: p = (%f, %f, %f)\n", q[0], q[1], q[2]);

  if (TMath::Abs(fBeamCrossingAngle) > 0.)
    BeamCrossing(q); // applying crossing angle
  if (TMath::Abs(fBeamDivergence) > 0.)
    BeamDivergence(q); // applying divergence
}

Double_t GeneratorSlowNucleons::Maxwell(Double_t m, Double_t p, Double_t T) {
  /* Relativistic Maxwell-distribution */
  Double_t ekin;
  ekin = TMath::Sqrt(p * p + m * m) - m;
  return (p * p * exp(-ekin / T));
}

//_____________________________________________________________________________
void GeneratorSlowNucleons::Lorentz(Double_t m, Double_t beta, Float_t *q) {
  /* Lorentz transform in the direction of q[2] */

  Double_t gamma = 1. / TMath::Sqrt((1. - beta) * (1. + beta));
  Double_t energy =
      TMath::Sqrt(m * m + q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
  q[2] = gamma * (q[2] + beta * energy);
}

//_____________________________________________________________________________
void GeneratorSlowNucleons::BeamCrossing(Float_t *pLab) {
  // Applying beam crossing angle
  //
  pLab[1] = pLab[2] * TMath::Sin(fBeamCrossingAngle) +
            pLab[1] * TMath::Cos(fBeamCrossingAngle);
  pLab[2] = pLab[2] * TMath::Cos(fBeamCrossingAngle) -
            pLab[1] * TMath::Sin(fBeamCrossingAngle);
  if (fDebug == 1) {
    printf(" Beam crossing angle = %f mrad ", fBeamCrossingAngle * 1000.);
    printf("  p = (%f, %f, %f)\n", pLab[0], pLab[1], pLab[2]);
  }
}
//_____________________________________________________________________________
void GeneratorSlowNucleons::BeamDivergence(Float_t *pLab) {
  // Applying beam divergence
  //
  const Float_t k2PI = 2. * TMath::Pi();
  const Float_t kRaddeg = 180. / TMath::Pi();
  const Float_t kDegrad = TMath::Pi() / 180.;
  Double_t pmod =
      TMath::Sqrt(pLab[0] * pLab[0] + pLab[1] * pLab[1] + pLab[2] * pLab[2]);

  Double_t tetdiv = fBeamDivEvent;
  Double_t fidiv = (gRandom->Rndm()) * 2. * TMath::Pi();

  Double_t tetpart =
      TMath::ATan2(TMath::Sqrt(pLab[0] * pLab[0] + pLab[1] * pLab[1]), pLab[2]);
  Double_t fipart = 0.;
  if (TMath::Abs(pLab[1]) > 0. || TMath::Abs(pLab[0]) > 0.)
    fipart = TMath::ATan2(pLab[1], pLab[0]);
  if (fipart < 0.) {
    fipart = fipart + 2. * TMath::Pi();
  }
  tetdiv = tetdiv * kRaddeg;
  fidiv = fidiv * kRaddeg;
  tetpart = tetpart * kRaddeg;
  fipart = fipart * kRaddeg;

  Double_t angleSum[2] = {0., 0.};
  AddAngle(tetpart, fipart, tetdiv, fidiv, angleSum);

  Double_t tetsum = angleSum[0];
  Double_t fisum = angleSum[1];
  // printf("tetpart %f fipart %f tetdiv %f fidiv %f angleSum %f
  // %f\n",tetpart,fipart,tetdiv,fidiv,angleSum[0],angleSum[1]);
  tetsum = tetsum * kDegrad;
  fisum = fisum * kDegrad;

  pLab[0] = pmod * TMath::Sin(tetsum) * TMath::Cos(fisum);
  pLab[1] = pmod * TMath::Sin(tetsum) * TMath::Sin(fisum);
  pLab[2] = pmod * TMath::Cos(tetsum);
  if (fDebug == 1) {
    printf(" Beam divergence %f mrad ", fBeamDivEvent * 1000.);
    printf("  p = (%f, %f, %f)\n", pLab[0], pLab[1], pLab[2]);
  }
}

//_____________________________________________________________________________
void GeneratorSlowNucleons::AddAngle(Double_t theta1, Double_t phi1,
                                     Double_t theta2, Double_t phi2,
                                     Double_t *angleSum) {
  // Calculating the sum of 2 angles
  Double_t temp = -1.;
  Double_t conv = 180. / TMath::ACos(temp);

  Double_t ct1 = TMath::Cos(theta1 / conv);
  Double_t st1 = TMath::Sin(theta1 / conv);
  Double_t cp1 = TMath::Cos(phi1 / conv);
  Double_t sp1 = TMath::Sin(phi1 / conv);
  Double_t ct2 = TMath::Cos(theta2 / conv);
  Double_t st2 = TMath::Sin(theta2 / conv);
  Double_t cp2 = TMath::Cos(phi2 / conv);
  Double_t sp2 = TMath::Sin(phi2 / conv);
  Double_t cx = ct1 * cp1 * st2 * cp2 + st1 * cp1 * ct2 - sp1 * st2 * sp2;
  Double_t cy = ct1 * sp1 * st2 * cp2 + st1 * sp1 * ct2 + cp1 * st2 * sp2;
  Double_t cz = ct1 * ct2 - st1 * st2 * cp2;

  Double_t rtetsum = TMath::ACos(cz);
  Double_t tetsum = conv * rtetsum;
  if (TMath::Abs(tetsum) < 1e-4 || tetsum == 180.)
    return;

  temp = cx / TMath::Sin(rtetsum);
  if (temp > 1.)
    temp = 1.;
  if (temp < -1.)
    temp = -1.;
  Double_t fisum = conv * TMath::ACos(temp);
  if (cy < 0) {
    fisum = 360. - fisum;
  }

  angleSum[0] = tetsum;
  angleSum[1] = fisum;
}
