void test()
{
  gSystem->Load("libpythia6");
  gSystem->Load("libEGPythia6");
  auto *jpsi = new GeneratorParam(10, new GeneratorParamMUONlib(), GeneratorParamMUONlib::kJpsiFamily, "Vogt PbPb");
  jpsi->SetPtRange(0, 100);
  jpsi->SetYRange(-1., +1.);
  jpsi->SetDecayer(new TPythia6Decayer());
  jpsi->SetForceDecay(kDiElectron);
  jpsi->Init();
  jpsi->Draw("");
  //  auto n3 = jpsi->GetPt()->Integral(0., 10., 1.e-6);
  jpsi->GenerateEvent();
  auto particles = new TClonesArray("TParticle", 10);
  auto n = jpsi->ImportParticles(particles, "all");
  printf("Number of generated particles %5d \n", n);
    
    for (Int_t i = 0; i < n; i++) {
        auto p = *(TParticle*) (particles->At(i));
            printf("%5d %5d %5d %5d %5d\n", i, p.GetPdgCode(), p.GetFirstMother(), p.GetFirstDaughter(), p.GetLastDaughter());
            
    }
}
