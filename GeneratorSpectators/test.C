#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TH1F.h>
#include <TCanvas.h>
#include "GeneratorSpectators.h"

#endif

void test(int nev=1)
{
    GeneratorSpectators* spgen = new GeneratorSpectators();
    spgen->SetParticle(2112); //Neutrons
    spgen->SetMomentum(2510.);
    spgen->SetDirection(0, 0., 0., -1.);
    spgen->SetFermi();
    spgen->SetDivergence(0.000032); //beam divergence in murad
    spgen->SetCrossing(0.000060, 2);//beam crossing angle in murad

    spgen->Init();

    double upper=4000.;
    if(spgen->GetZDirection()<0.) upper = -4000.;
    TH1F *hpz = new TH1F("hpz", "Spectators p_{z}", 200, 0., upper);
    TH1F *hpt = new TH1F("hpt", "Spectators p_{t}", 100, 0., 2.);

    for (Int_t iev = 0; iev < nev; iev++) {
      spgen->GenerateEvent();

      //printf(" \t Importing generated particles \n");
      auto particles = new TClonesArray("TParticle", 10);
      auto npart = spgen->ImportParticles(particles, "all");
      //printf(" - Number of generated particles %5d \n", npart);

      for (Int_t i = 0; i < npart; i++) {
        auto p = *(TParticle*) (particles->At(i));
        //printf("(%d) PDG:%5d pt=%1.2f GeV pz=%1.2f GeV\n", iev, p.GetPdgCode(), p.Pt(), p.Pz());
        hpt->Fill(p.Pt());
        hpz->Fill(p.Pz());
      }
    }

      TCanvas *c = new TCanvas("c","Momentum",0,0,1200,600);
      c->Divide(2,1);
      c->cd(1);
      hpt->SetLineColor(kAzure+10);
      hpt->SetLineWidth(2);
      hpt->Draw("e");
      hpt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      c->cd(2);
      hpz->SetLineColor(kPink);
      hpz->SetLineWidth(2);
      hpz->Draw("e");
      hpz->GetXaxis()->SetTitle("p_{Z} (GeV/c)");
}
