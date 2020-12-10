// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Class to generate spectators for ZDC simulations

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <assert.h>

#include "GeneratorSpectators.h"

ClassImp(GeneratorSpectators)

//_____________________________________________________________________________
GeneratorSpectators::GeneratorSpectators()
   :TGenerator("GeneratorSpectators", "GeneratorSpectators"), fDebug(0),
   fPDGcode(0), fPmax(0), fPseudoRapidity(0), fCosx(0), fCosy(0), fCosz(0),
   fFermiflag(1),fBeamDiv(0), fBeamCrossAngle(0), fBeamCrossPlane(0) {
  //
  // Default constructor
}

//_____________________________________________________________________________
GeneratorSpectators::GeneratorSpectators(Int_t npart)
   :TGenerator("GeneratorSpectators", "GeneratorSpectators"){
  //
  // Standard constructor
  //
  fName = "GeneratorSpectators";
  fTitle = "Generation of Test Particles for ZDCs";
  fDebug = 0;
  SetParticle();
  SetMomentum();
  SetDirection();
  SetFermi();
  SetDivergence();
  SetCrossing();

  for(Int_t i=0; i<201; i++){
     fProbintp[i] = 0;
     fProbintn[i] = 0;
     fPp[i] = 0;
  }
}

//_____________________________________________________________________________
void GeneratorSpectators::Init()
{
  //Initialize Fermi momentum distributions for Pb-Pb
  //
  printf("\n **** GeneratorSpectators initialization:\n");
  printf("   ParticlePDG: %d, Track: cosx = %f cosy = %f cosz = %f \n", fPDGcode, fCosx, fCosy, fCosz);
  printf("   Fermi flag: %d, Beam divergence: %f, Crossing angle: %f, plane: %d\n\n",
             fFermiflag, fBeamDiv, fBeamCrossAngle, fBeamCrossPlane);

  FermiTwoGaussian(208.);
}

//_____________________________________________________________________________
void GeneratorSpectators::GenerateEvent()
{
  //
  // Generate one trigger particle (n or p)
  //
  //printf("GeneratorSpectators::GenerateEvent()\n");
  fParticles->Clear();

  Double_t pLab[3] = {0.,0.,0.};
  Double_t fP[3] = {0.,0.,0.}, fBoostP[3] = {0.,0.,0.};
  Double_t ptot = fPmax;
  if(fPseudoRapidity==0.){
    pLab[0] = ptot*fCosx;
    pLab[1] = ptot*fCosy;
    pLab[2] = ptot*fCosz;
  }
  else{
    Float_t scang = 2*TMath::ATan(TMath::Exp(-(fPseudoRapidity)));
    pLab[0] = -ptot*TMath::Sin(scang);
    pLab[1] = 0.;
    pLab[2] = ptot*TMath::Cos(scang);
  }
  for(int i=0; i<3; i++) fP[i] = pLab[i];
 if(fDebug==1) printf("\n Particle momentum before divergence and crossing: ");
 if(fDebug==1) printf(" 	pLab = (%f, %f, %f)\n", pLab[0], pLab[1], pLab[2]);

  // Beam divergence and crossing angle
  if(TMath::Abs(fBeamCrossAngle)>0.) {
    BeamCrossing(pLab);
    for(int i=0; i<3; i++) fP[i] = pLab[i];
  }
  if(TMath::Abs(fBeamDiv)>0.) {
    BeamDivergence(pLab);
    for(int i=0; i<3; i++) fP[i] = pLab[i];
  }

  Double_t mass = TDatabasePDG::Instance()->GetParticle(fPDGcode)->Mass();
  Double_t ddp[3] = {0.,0.,0.}, dddp[3] = {0.,0.,0.};
  Double_t fP0 = 0., dddp0 = 0.;
  // If required apply the Fermi momentum
  if(fFermiflag==1){
    if((fPDGcode==kProton) || (fPDGcode==kNeutron)) ExtractFermi(fPDGcode, ddp);
    fP0 = TMath::Sqrt(fP[0]*fP[0] + fP[1]*fP[1] + fP[2]*fP[2] + mass*mass);
    for(int i=0; i<3; i++) dddp[i] = ddp[i];
    dddp0 = TMath::Sqrt(dddp[0]*dddp[0] + dddp[1]*dddp[1] + dddp[2]*dddp[2] + mass*mass);

    TVector3 b(fP[0]/fP0, fP[1]/fP0, fP[2]/fP0);
    TLorentzVector pFermi(dddp[0], dddp[1], dddp[2], dddp0);

    pFermi.Boost(b);
    for(int i=0; i<3; i++){
       fBoostP[i] = pFermi[i];
       fP[i] = pFermi[i];
    }

  }
 if(fDebug==1) printf(" ### Particle momentum = (%f, %f, %f)\n", fP[0], fP[1], fP[2]);

  Double_t energy = TMath::Sqrt(fP[0]*fP[0]+fP[1]*fP[1]+fP[2]*fP[2]+mass*mass);

  auto part = new TParticle(fPDGcode, 1, -1, -1, -1, -1, fP[0], fP[1], fP[2], energy,
                              0., 0., 0., 0.);
  fParticles->Add(part);
}

// -----------------------------------------------------------------------
int GeneratorSpectators::ImportParticles(TClonesArray *particles, Option_t *option) {
  if (particles == 0) return 0;
  TClonesArray &clonesParticles = *particles;
  clonesParticles.Clear();
  Int_t numpart = fParticles->GetEntries();
  for (int i = 0; i < numpart; i++) {
    TParticle *particle = (TParticle *)fParticles->At(i);
    new (clonesParticles[i]) TParticle(*particle);
  }
  return numpart;
}

//_____________________________________________________________________________
void GeneratorSpectators::FermiTwoGaussian(Float_t A)
{
//
// Momenta distributions according to the "double-gaussian"
// distribution (Ilinov) - equal for protons and neutrons
//
   Double_t sig1 = 0.113;
   Double_t sig2 = 0.250;
   Double_t alfa = 0.18*(TMath::Power((A/12.), (Float_t)1/3));
   Double_t xk = (2*TMath::TwoPi())/((1.+alfa)*(TMath::Power(TMath::TwoPi(),1.5)));

   for(Int_t i=1; i<201; i++){
      Double_t p = i*0.005;
      fPp[i] = p;
      Double_t e1 = (p*p)/(2.*sig1*sig1);
      Double_t e2 = (p*p)/(2.*sig2*sig2);
      Double_t f1 = TMath::Exp(-(e1));
      Double_t f2 = TMath::Exp(-(e2));
      Double_t probp = xk*p*p*(f1/(TMath::Power(sig1,3.))+
                      alfa*f2/(TMath::Power(sig2,3.)))*0.005;
      fProbintp[i] = fProbintp[i-1] + probp;
      fProbintn[i] = fProbintp[i];
   }
  if(fDebug==1) printf("		Initialization of Fermi momenta distribution \n");
}
//_____________________________________________________________________________
void GeneratorSpectators::ExtractFermi(Int_t id, Double_t *ddp)
{
  //
  // Compute Fermi momentum for spectator nucleons
  //
  Int_t index=0;
  Float_t xx = gRandom->Rndm();
  assert ( id==kProton || id==kNeutron );
  if(id==kProton){
    for(Int_t i=1; i<201; i++){
       if((xx>=fProbintp[i-1]) && (xx<fProbintp[i])) break;
       index = i;
    }
  }
  else {
    for(Int_t i=1; i<201; i++){
       if((xx>=fProbintn[i-1]) && (xx<fProbintn[i])) break;
       index = i;
    }
  }
  Float_t pext = fPp[index]+0.001;
  Float_t phi = TMath::TwoPi()*(gRandom->Rndm());
  Float_t cost = (1.-2.*(gRandom->Rndm()));
  Float_t tet = TMath::ACos(cost);
  ddp[0] = pext*TMath::Sin(tet)*TMath::Cos(phi);
  ddp[1] = pext*TMath::Sin(tet)*TMath::Sin(phi);
  ddp[2] = pext*cost;

 if(fDebug==1) printf(" Fermi momentum: p = (%f, %f, %f )\n\n",ddp[0],ddp[1],ddp[2]);
}
//_____________________________________________________________________________
void GeneratorSpectators::BeamCrossing(Double_t *pLab)
{
  // Applying beam crossing angle
  //
  pLab[1] = pLab[2]*TMath::Sin(fBeamCrossAngle)+pLab[1]*TMath::Cos(fBeamCrossAngle);
  pLab[2] = pLab[2]*TMath::Cos(fBeamCrossAngle)-pLab[1]*TMath::Sin(fBeamCrossAngle);
  //
 if(fDebug==1) printf(" Beam crossing angle = %f mrad -> ", fBeamCrossAngle*1000.);
 if(fDebug==1) printf("  p = (%f, %f, %f)\n",pLab[0],pLab[1],pLab[2]);

}
//_____________________________________________________________________________
void GeneratorSpectators::BeamDivergence(Double_t *pLab)
{
  // Applying beam divergence and crossing angle
  //
  Double_t pmq = 0.;
  for(int i=0; i<3; i++) pmq = pmq+pLab[i]*pLab[i];
  Double_t pmod = TMath::Sqrt(pmq);

  Double_t rvec = gRandom->Gaus(0.0,1.0);
  Double_t tetdiv = fBeamDiv * TMath::Abs(rvec);
  Double_t fidiv = (gRandom->Rndm())*TMath::TwoPi();

  Double_t tetpart = TMath::ATan2(TMath::Sqrt(pLab[0]*pLab[0]+pLab[1]*pLab[1]),pLab[2]);
  Double_t fipart = 0.;
  if(pLab[1]!=0. || pLab[0]!=0.) fipart = TMath::ATan2(pLab[1],pLab[0]);
  else fipart = 0.;
  if(fipart<0.) {fipart = fipart+TMath::TwoPi();}
  tetdiv = tetdiv*TMath::RadToDeg();
  fidiv = fidiv*TMath::RadToDeg();
  tetpart = tetpart*TMath::RadToDeg();
  fipart = fipart*TMath::RadToDeg();
  Double_t angleSum[2] = {0., 0.};
  AddAngle(tetpart,fipart,tetdiv,fidiv,angleSum);
  Double_t tetsum = angleSum[0];
  Double_t fisum  = angleSum[1];
  tetsum = tetsum*TMath::DegToRad();
  fisum = fisum*TMath::DegToRad();
  pLab[0] = pmod*TMath::Sin(tetsum)*TMath::Cos(fisum);
  pLab[1] = pmod*TMath::Sin(tetsum)*TMath::Sin(fisum);
  pLab[2] = pmod*TMath::Cos(tetsum);
  //
 if(fDebug==1) printf(" Beam divergence = %f mrad -> ", fBeamDiv*1000.);
 if(fDebug==1) printf("  p = (%f, %f, %f)\n",pLab[0],pLab[1],pLab[2]);
}

//_____________________________________________________________________________
void  GeneratorSpectators::AddAngle(Double_t theta1, Double_t phi1, Double_t theta2,
  	       Double_t phi2, Double_t *angleSum)
{
  // Calculating the sum of 2 angles
  Double_t temp = -1.;
  Double_t conv = 180./TMath::ACos(temp);
  Double_t ct1 = TMath::Cos(theta1/conv);
  Double_t st1 = TMath::Sin(theta1/conv);
  Double_t cp1 = TMath::Cos(phi1/conv);
  Double_t sp1 = TMath::Sin(phi1/conv);
  Double_t ct2 = TMath::Cos(theta2/conv);
  Double_t st2 = TMath::Sin(theta2/conv);
  Double_t cp2 = TMath::Cos(phi2/conv);
  Double_t sp2 = TMath::Sin(phi2/conv);
  Double_t cx = ct1*cp1*st2*cp2+st1*cp1*ct2-sp1*st2*sp2;
  Double_t cy = ct1*sp1*st2*cp2+st1*sp1*ct2+cp1*st2*sp2;
  Double_t cz = ct1*ct2-st1*st2*cp2;

  Double_t rtetsum = TMath::ACos(cz);
  Double_t tetsum = conv*rtetsum;
  Double_t fisum = 0;
  if (tetsum == 0. || tetsum == 180.)
  {
    fisum = 0.;
    return;
  }
  temp = cx/TMath::Sin(rtetsum);
  if(temp>1.) temp=1.;
  else if (temp < -1.) temp = -1.;
  fisum = conv*TMath::ACos(temp);
  if(cy<0) {fisum = 360.-fisum;}
  angleSum[0] = tetsum;
  angleSum[1] = fisum;
}
