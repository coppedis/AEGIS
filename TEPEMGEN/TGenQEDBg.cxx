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

// Event generator for background from e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 TeV/nucleon, integrated over specific readout cycle
// Derived from AliGenEpEmv1
// Author: ruben.shahoyan@cern.ch
//
// Revised on September 2018 for ALICEo2: Roberto Preghenella (preghenella@bo.infn.it)
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
//    AliGenQEDBg *gener = new AliGenQEDBg();
//    gener->SetXXXRange(); // Set kinematics range
//    gener->SetLumiIntTime(double lumi, double sec); // luminosity and intergration time in seconds
//    gener->Init();
// Event generation:
//    gener->Generate(); // Produce poissonian number of e+e- pair with average number
//    corresponding to requested integration time with given beam luminosity
//    Each pair has its own vertex, and time randomly distributed at originT 
//    and originT+ integration time                    
//
//    For details of pair generation see AliGenEpEmv1.cxx by Yuri.Kharlov@cern.ch

/*
  The most useful way of using it is in the overlay with other generators (using the cocktail):
  In the Config.C 
  // load the library:
  gSystem->Load("libTEPEMGEN");
  //
  // add to AliGenCocktail as:
  AliGenCocktail *cocktail = new AliGenCocktail();
  // ... setup other stuff
  // 
  // QED background
  AliGenQEDBg*  genBg = new AliGenQEDBg();
  genBg->SetEnergyCMS(5500);
  genBg->SetProjectile("A", 208, 82);
  genBg->SetTarget    ("A", 208, 82);
  genBg->SetYRange(-6.,3);
  genBg->SetPtRange(1.e-3,1.0);      // Set pt limits (GeV) for e+-: 1MeV corresponds to max R=13.3mm at 5kGaus
  genBg->SetLumiIntTime(6.e27,20e-6); // luminosity and integration time
  //
  cocktail->AddGenerator(genBg,"QEDep",1);
  // We need to generate independent vertex for every event, this must be set after adding to cocktail
  // Note that the IO origin and sigma's is transferred from AliGenCocktail to daughter generators
  genBg->SetVertexSource(kInternal);
*/

#include "TGenQEDBg.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TEpEmGen.h>

ClassImp(TGenQEDBg);

//------------------------------------------------------------

TGenQEDBg::TGenQEDBg()
:  fLumi(0)
  ,fIntTime(0)
  ,fPairsInt(-1)
{
}

//____________________________________________________________
TGenQEDBg::~TGenQEDBg()
{
  // Destructor
}

//____________________________________________________________
void TGenQEDBg::Init()
{
  // Initialisation:
  printf("Will estimate QED bg. for L=%e cm^-2*s^-1 and Integration Time of %e s.",fLumi,fIntTime);
  if (fLumi<=0 || fIntTime<=0) {
    printf("One of parameters is not set properly, no pairs will be generated");
    abort();
  }
  //
  // initialize the generator of e+e- pair production
  TGenEpEmv1::Init();
  //
  fPairsInt = 0;
  fPairsInt = fXSection*1e-24*fLumi*fIntTime; // xsestion is in barn!
  printf("Estimated x-secion: %e+-%eb, <Npairs>=%e per %e time interval",
	 fXSection,fXSectionEps*fXSection,fPairsInt,fIntTime);
  //
}

//____________________________________________________________
void TGenQEDBg::GenerateEvent()
{
  //
  // Generate poissian <fPairsInt> e+e- pairs, each one with its vertex
  //

  fParticles->Clear();

  Float_t random[6];
  Float_t origin[3];
  Float_t time = 0.;
  //
  int npairs=0,nt,id;;
  if (fPairsInt>0) npairs=gRandom->Poisson(fPairsInt);
  if (fDebug == 1)
    printf("<nQED>=%e -> %d pairs will be generated",fPairsInt,npairs);
  if (npairs<1) return;
  for (int i=0;i<npairs;i++) {
    // each pair has its own vertex and time
    for (int j=0;j<3;j++) origin[j]=fOrigin[j];
    Rndm(random,6);
    for (int j=0;j<3;j++) {
      origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
    time = fTimeOrigin+gRandom->Rndm()*fIntTime;
    //
    TGenEpEmv1::GeneratePair(origin[0], origin[1], origin[2], time);
  }
  fEvent++;
  //
}

//__________________________________________________________
void TGenQEDBg::SetLumiIntTime(double lumi, double intTime)
{
  // assign luminosity and integration time
  if (lumi<=0) {
    printf("Luminosity must be positive, in cm^-2*s^-1, %e asked",lumi);
    abort();
  }
  if (intTime<=0) {
    printf("Integration time must be positive, in seconnds, %e asked",intTime);
    abort();
  }
  fLumi = lumi;
  fIntTime = intTime;
  //
}
