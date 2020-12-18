void test()
{
    GeneratorSlowNucleons* slow = new GeneratorSlowNucleons();
    slow->SetSlowNucleonModel(new SlowNucleonModelExp());
    slow->SetNominalCmsEnergy(13000);
    slow->SetTarget();
    slow->SetTemperature();
    slow->SetBetaSource();
    slow->SetProtonDirection();
    slow->SetDebug();
    slow->SetPmax();
    slow->Init();

    slow->SetNcoll(7);
    slow->GenerateEvent();

  auto particles = new TClonesArray("TParticle", 10);
  auto n = slow->ImportParticles(particles, "all");
  printf("Number of generated particles %5d \n", n);

    for (Int_t i = 0; i < n; i++) {
        auto p = *(TParticle*) (particles->At(i));
            printf(" PDG:%d pt=%1.2f GeV pz=%1.2f GeV\n", p.GetPdgCode(), p.Pt(), p.Pz());
    }
}
