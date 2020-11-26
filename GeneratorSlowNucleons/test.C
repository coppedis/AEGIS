void test()
{
    GeneratorSlowNucleons* slow = new GeneratorSlowNucleons();
    slow->SetSlowNucleonModel(new SlowNucleonModelExp());
    slow->SetNcoll(7);
    slow->SetNominalCmsEnergy(13000);
    slow->SetTarget();
    slow->Init();
    slow->GenerateEvent();
    
  auto particles = new TClonesArray("TParticle", 10);
  auto n = slow->ImportParticles(particles, "all");
  printf("Number of generated particles %5d \n", n);
    
    for (Int_t i = 0; i < n; i++) {
        auto p = *(TParticle*) (particles->At(i));
            printf("%5d %5d %5d %5d %5d\n", i, p.GetPdgCode(), p.GetFirstMother(), p.GetFirstDaughter(), p.GetLastDaughter());
            
    }
}
