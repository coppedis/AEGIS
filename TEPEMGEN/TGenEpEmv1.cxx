/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 *                                                                        *
 * Copyright(c) 1997, 1998, 2002, Adrian Alscher and Kai Hencken          *
 * See $ALICE_ROOT/EpEmGen/diffcross.f for full Copyright notice          *
 *                                                                        *
 *                                                                        *
 * Copyright(c) 2002 Kai Hencken, Yuri Kharlov, Serguei Sadovsky          *
 * See $ALICE_ROOT/EpEmGen/epemgen.f for full Copyright notice            *
 *                                                                        *
 **************************************************************************/

/* $Id$ */

// Event generator of single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 TeV/nucleon.
// The generator is based on 5-dimentional differential cross section of the process.
//%
// References:
// [1] "Multiple electromagnetic electron positron pair production in
//      relativistic heavy ion collisions".
//      Adrian Alscher, Kai Hencken, Dirk Trautmann, and Gerhard Baur,
//      Phys. Rev. A55 (1997) 396.
// [2] K.Hencken, Yu.Kharlov, S.Sadovsky, Internal ALICE Note 2002-27.
//%
// Usage:
// Initialization:
//    AliGenEpEmv1 *gener = new AliGenEpEmv1();
//    gener->SetXXXRange(); // Set kinematics range
//    gener->Init();
// Event generation:
//    gener->Generate(); // Produce one e+e- pair with the event weight assigned 
//                       // to each track. The sum of event weights, divided by 
//                       // the total number of generated events, gives the 
//                       // integral cross section of the process of e+e- pair 
//                       // production in the above mentioned kinematics range.
//                       // Sum of the selected event weights, divided by the total 
//                       // number of generated events, gives the integral cross 
//                       // section corresponded to the set of selected events
//%
// The generator consists of several modules:
// 1) $ALICE_ROOT/EpEmGen/diffcross.f:
//    Exact calculation of the total differential e+ e- -pair production
//    in Relativistic Heavy Ion Collisions for a point particle in an
//    external field approach. See full comments in the mentioned file.
// 2) $ALICE_ROOT/EpEmGen/epemgen.f:
//    Generator of e+e- pairs produced in PbPb collisions at LHC
//    it generates events according to the parametrization of the
//    differential cross section. Produces events have weights calculated
//    by the exact differential cross section calculation (diffcross.f).
//    See full comments in the mentioned file.
// 3) Class TEpEmGen:
//    Interface from the fortran event generator to ALIROOT
// 4) Class AliGenEpEmv1:
//    The event generator to call within ALIROOT
//%
// Author of this module: Yuri.Kharlov@cern.ch
// 9 October 2002
//
// Revised on September 2018 for ALICEo2: Roberto Preghenella (preghenella@bo.infn.it)

#include "TGenEpEmv1.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TEpEmGen.h>

ClassImp(TGenEpEmv1);

//------------------------------------------------------------

TGenEpEmv1::TGenEpEmv1():
  TEpEmGen(),
  fMass(0),
  fDebug(0),
  fEvent(0),
  fYMin(-100.), fYMax(100.),
  fPhiMin(0.), fPhiMax(360.),
  fPtMin(0.), fPtMax(1.e10),
  fTimeOrigin(0.),
  fXSection(-1.),
  fXSectionEps(1e-2),
  fMinXSTest(1000),
  fMaxXSTest(10000000)
{
  // Default constructor
  for (Int_t i = 0; i < 3; ++i) {
    fOrigin[i] = 0.;
    fOsigma[i] = 0.;
  }
}

//____________________________________________________________
TGenEpEmv1::~TGenEpEmv1()
{
  // Destructor
}

//____________________________________________________________
bool TGenEpEmv1::Init()
{
  // Initialisation:
  // 1) define a generator
  // 2) initialize the generator of e+e- pair production
  
  fMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  if (fPtMin == 0) fPtMin = 1.E-04; // avoid zero pT
  Initialize(fYMin, fYMax, fPtMin, fPtMax);
  fEvent = 0;
  //
  // calculate XSection
  double err = 0;
  fXSection = CalcXSection(fXSectionEps,fMinXSTest,fMaxXSTest,err);
  if (fXSection<=0 || err/fXSection>fXSectionEps) {
    return false;
  }
  fXSectionEps = err/fXSection;
  return true;
}

//____________________________________________________________
double TGenEpEmv1::CalcXSection(double eps, int triMin, int triMax, double& err)
{
  if (eps<1e-4) {
    eps = 1e-4;
  }
  int ngen = 0;
  printf("Estimating x-section with min.relative precision of %f and min/max test: %d/%d\n",
	 eps,triMin,triMax);
  double yElectron,yPositron,xElectron,xPositron,phi12,weight;
  double xSect = -1;
  err = -1;
  //
  do {
    TEpEmGen::GenerateEvent(fYMin,fYMax,fPtMin,fPtMax,yElectron,yPositron,xElectron,xPositron,phi12,weight);
    if (++ngen>triMin) { // ensure min number of tests
      xSect = TEpEmGen::GetXsection()*1000;
      err = TEpEmGen::GetDsection()*1000;
    }
  } while(!((xSect>0 && err/xSect<eps) || ngen>triMax));
  //
  if (xSect<=0) {
    printf("Failed to estimate X-section after %d trials\n",ngen);
    abort();
  }
  printf("X-section = %e with %e error after %d trials",xSect,err,ngen);
  return xSect;
}


//____________________________________________________________
void TGenEpEmv1::GenerateEvent()
{
  fParticles->Clear();
  Float_t random[6];
  Float_t origin[3];
  Float_t time = 0.;
  for (int j=0;j<3;j++) origin[j]=fOrigin[j];
  time = fTimeOrigin;
  Rndm(random,6);
  for (int j=0;j<3;j++) {
    origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
      TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
  }
  Rndm(random,2);
  time += fOsigma[2]/TMath::Ccgs()*
    TMath::Cos(2*random[0]*TMath::Pi())*
    TMath::Sqrt(-2*TMath::Log(random[1]));
  GeneratePair(origin[0], origin[1], origin[2], time);
  fEvent++;
  if (fEvent%1000 == 0) {
    printf("=====> TGenEpEmv1::Generate(): \n   Event %d, sigma=%f +- %f kb\n",
  	   fEvent, GetXsection(), GetDsection());
  }
}

//____________________________________________________________
void TGenEpEmv1::GeneratePair(Double_t vx, Double_t vy, Double_t vz, Double_t vt)
{
  //
  // Generate one e+e- pair
  // Gaussian smearing on the vertex is done if selected. 
  //%
  // Each produced e+e- pair is defined by the following variables:
  // rapidities of e-, e+ (yElectron,yPositron)
  // log10(pt in MeV/c) of e-, e+ (xElectron,xPositron)
  // azymuth angles between e- and e+ (phi12)
  //%
  // On output an event weight is given (weight) which is assigned to each track.
  // The sum of event weights, divided by the total number of generated events, 
  // gives the integral cross section of the e+e- pair production in the   
  // selected kinematics range.	  
  //

  Float_t random[6];
  Float_t polar[3]= {0,0,0};
  Float_t p[3];

  Double_t ptElectron,ptPositron, phiElectron,phiPositron, mt, etot;
  Double_t phi12=0,xElectron=0,xPositron=0,yElectron=0,yPositron=0,weight=0;
  Int_t   j, nt, id;
  
  TEpEmGen::GenerateEvent(fYMin,fYMax,fPtMin,fPtMax,
			  yElectron,yPositron,xElectron,xPositron,phi12,weight);
  if (fDebug == 1)
    printf("TGenEpEmv1::Generate(): y=(%f,%f), x=(%f,%f), phi=%f\n",
	   yElectron,yPositron,xElectron,xPositron,phi12);
  
  Rndm(random,1);
  ptElectron  = TMath::Power(10,xElectron) * 1.e-03;;
  ptPositron  = TMath::Power(10,xPositron) * 1.e-03;;
  phiElectron = fPhiMin + random[0] * (fPhiMax-fPhiMin);
  phiPositron = phiElectron + phi12;

  // Produce electron
  mt = TMath::Sqrt(ptElectron*ptElectron + fMass*fMass);
  p[0] = ptElectron*TMath::Cos(phiElectron);
  p[1] = ptElectron*TMath::Sin(phiElectron);
  p[2] = mt*TMath::SinH(yElectron);
  etot = TMath::Sqrt(fMass*fMass + p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  id =  11;
  if (fDebug == 2)
    printf("id=%+3d, p = (%+11.4e,%+11.4e,%+11.4e) GeV\n",id,p[0],p[1],p[2]);
  TParticle *electron = new TParticle(id, 1, -1, -1, -1, -1, p[0], p[1], p[2], etot, vx, vy, vz, vt);
  fParticles->Add(electron);
  
  // Produce positron
  mt = TMath::Sqrt(ptPositron*ptPositron + fMass*fMass);
  p[0] = ptPositron*TMath::Cos(phiPositron);
  p[1] = ptPositron*TMath::Sin(phiPositron);
  p[2] = mt*TMath::SinH(yPositron);
  etot = TMath::Sqrt(fMass*fMass + p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  id = -11;
  if (fDebug == 2)
    printf("id=%+3d, p = (%+11.4e,%+11.4e,%+11.4e) GeV\n",id,p[0],p[1],p[2]);
  TParticle *positron = new TParticle(id, 1, -1, -1, -1, -1, p[0], p[1], p[2], etot, vx, vy, vz, vt);
  fParticles->Add(positron);
}

