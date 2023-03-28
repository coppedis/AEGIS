#ifndef EXODUSDECAYER_H
#define EXODUSDECAYER_H

// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//---------------------------------------------------------------------------------------------------
//                                 
// Generate electron-pair mass distributions for Dalitz decays according
// to the Kroll-Wada parametrization: N. Kroll, W. Wada: Phys. Rev 98(1955)1355
// and generate electron-pair mass distributions for resonances according
// to the Gounaris-Sakurai parametrization: G.J. Gounaris, J.J. Sakurai: Phys.Rev.Lett. 21(1968)244 
//
// For the electromagnetic form factor the parameterization from
// Lepton-G is used: L.G. Landsberg et al.: Phys. Rep. 128(1985)301
//
// Ralf Averbeck (R.Averbeck@gsi.de) 
// Irem Erdemir  (irem.erdemir@cern.ch)
//
// adapted for O2: Daniel Samitz (daniel.samitz@cern.ch)
//
//---------------------------------------------------------------------------------------------------

#include "TVirtualMCDecayer.h"
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1.h>
#include "TDatabasePDG.h"

//class TH1F;
//class TClonesArray;

class ExodusDecayer : public TVirtualMCDecayer
{
 public:
    ExodusDecayer();
    virtual ~ExodusDecayer();
    virtual void    Init();
    virtual void    Decay(Int_t idpart,TLorentzVector* pparent);
    virtual Int_t   ImportParticles(TClonesArray *particles) {return -1;}
    virtual void    SetForceDecay(Int_t)                      {;}
    virtual void    ForceDecay()                              {;}
    virtual Float_t GetPartialBranchingRatio(Int_t /*ipart*/) {return -1;}
    virtual Float_t GetLifetime(Int_t /*kf*/)                 {return -1;}
    virtual void    ReadDecayTable()                          {;}
    
    virtual TH1F*   ElectronPairMassHistoPion()          {return  fEPMassPion;}
    virtual TH1F*   ElectronPairMassHistoEta()           {return  fEPMassEta;}
    virtual TH1F*   ElectronPairMassHistoEtaPrime()      {return  fEPMassEtaPrime;}
    virtual TH1F*   ElectronPairMassHistoEtaPrime_toOmega()      {return  fEPMassEtaPrime_toOmega;}
    virtual TH1F*   ElectronPairMassHistoRho()           {return  fEPMassRho;}
    virtual TH1F*   ElectronPairMassHistoOmega()         {return  fEPMassOmega;}
    virtual TH1F*   ElectronPairMassHistoOmegaDalitz()   {return  fEPMassOmegaDalitz;}
    virtual TH1F*   ElectronPairMassHistoPhi()           {return  fEPMassPhi;}
    virtual TH1F*   ElectronPairMassHistoPhiDalitz()     {return  fEPMassPhiDalitz;}
    virtual TH1F*   ElectronPairMassHistoPhiDalitz_toPi0()     {return  fEPMassPhiDalitz_toPi0;}
    virtual TH1F*   ElectronPairMassHistoJPsi()          {return  fEPMassJPsi;}

    virtual const   TLorentzVector* Products_pion()         const {return fProducts_pion;}
    virtual const   TLorentzVector* Products_eta()          const {return fProducts_eta;}
    virtual const   TLorentzVector* Products_etaprime()     const {return fProducts_etaprime;}
    virtual const   TLorentzVector* Products_etaprime_toOmega()     const {return fProducts_etaprime_toOmega;}
    virtual const   TLorentzVector* Products_rho()          const {return fProducts_rho;}
    virtual const   TLorentzVector* Products_omega()        const {return fProducts_omega;}
    virtual const   TLorentzVector* Products_omega_dalitz() const {return fProducts_omega_dalitz;}
    virtual const   TLorentzVector* Products_phi()          const {return fProducts_phi;}
    virtual const   TLorentzVector* Products_phi_dalitz()   const {return fProducts_phi_dalitz;}
    virtual const   TLorentzVector* Products_phi_dalitz_toPi0()   const {return fProducts_phi_dalitz_toPi0;}
    virtual const   TLorentzVector* Products_jpsi()         const {return fProducts_jpsi;}

 protected:
    // Histograms for electron pair mass
    TH1F*         fEPMassPion;          
    TH1F*         fEPMassEta;       
    TH1F*         fEPMassEtaPrime;
    TH1F*         fEPMassEtaPrime_toOmega;
    TH1F*         fEPMassRho;
    TH1F*         fEPMassOmega;
    TH1F*         fEPMassOmegaDalitz;
    TH1F*         fEPMassPhi;
    TH1F*         fEPMassPhiDalitz;
    TH1F*         fEPMassPhiDalitz_toPi0;
    TH1F*         fEPMassJPsi;

    TF1* fPol;

    // Decay products
    TLorentzVector  fProducts_pion[3];  
    TLorentzVector  fProducts_eta[3];  
    TLorentzVector  fProducts_etaprime[3];
    TLorentzVector  fProducts_etaprime_toOmega[3];
    TLorentzVector  fProducts_rho[2];
    TLorentzVector  fProducts_omega[2];
    TLorentzVector  fProducts_omega_dalitz[3];
    TLorentzVector  fProducts_phi[2];
    TLorentzVector  fProducts_phi_dalitz[3];
    TLorentzVector  fProducts_phi_dalitz_toPi0[3];
    TLorentzVector  fProducts_jpsi[2];


    Bool_t fInit;

 private:
    Double_t GounarisSakurai(Float_t mass, Double_t vmass, Double_t vwidth, Double_t emass);
    Double_t RhoShapeFromNA60(Float_t mass, Double_t vmass, Double_t vwidth, Double_t emass);
    Double_t Lorentz(Float_t mass, Double_t vmass, Double_t vwidth); 

    ClassDef(ExodusDecayer, 1)
};
#endif