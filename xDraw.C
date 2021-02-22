void xDraw(){
  TFile *fopen = new TFile("output_ggtree_mc.root");
  
  TCanvas *c1 = new TCanvas("c1");

  const char *title;
  //draw TH1F non-logY histograms
  TH1F *H_npho = (TH1F*)fopen->Get("h_npho");
  H_npho->GetXaxis()->SetTitle("photon number");
  H_npho->GetYasis()->SetTitle("Events");
  H_npho->Draw();
  title = H_npho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_njet = (TH1F*)fopen->Get("h_njet");
  H_njet->Draw();
  title = H_njet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_pjmass = (TH1F*)fopen->Get("h_pjmass");
  H_pjmass->Draw();
  title = H_pjmass->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_nMCpho = (TH1F*)fopen->Get("h_nMCpho");
  H_nMCpho->Draw();
  title = H_nMCpho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_nMCQstar = (TH1F*)fopen->Get("h_nMCQstar");
  H_nMCQstar->Draw();
  title = H_nMCQstar->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_MCQstarMass = (TH1F*)fopen->Get("h_MCQstarMass");
  H_MCQstarMass->Draw();
  title = H_MCQstarMass->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_nrealpho = (TH1F*)fopen->Get("h_nrealpho");
  H_nrealpho->Draw();
  title = H_nrealpho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));
  
  TH1F *H_njetGen = (TH1F*)fopen->Get("h_njetGen");
  H_njetGen->Draw();
  title = H_njetGen->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));
  
  TH1F *H_nrealjet = (TH1F*)fopen->Get("h_nrealjet");
  H_nrealjet->Draw();
  title = H_nrealjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_pjmass_phoLjet = (TH1F*)fopen->Get("h_pjmass_phoLjet");
  H_pjmass_phoLjet->Draw();
  title = H_pjmass_phoLjet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));
  
  TH1F *H_pjmass_HPtpho = (TH1F*)fopen->Get("h_pjmass_HPtpho");
  H_pjmass_HPtpho->Draw();
  title = H_pjmass_HPtpho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));
  
  //draw TH1F logY histograms
  c1->SetLogy();
  c1->Update();
  
  
  TH1F *H_dr_pho = (TH1F*)fopen->Get("h_dr_pho");
  H_dr_pho->Draw();
  title = H_dr_pho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_dpt_pho = (TH1F*)fopen->Get("h_dpt_pho");
  H_dpt_pho->Draw();
  title = H_dpt_pho->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_dr_jet = (TH1F*)fopen->Get("h_dr_jet");
  H_dr_jet->Draw();
  title = H_dr_jet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));

  TH1F *H_dpt_jet = (TH1F*)fopen->Get("h_dpt_jet");
  H_dpt_jet->Draw();
  title = H_dpt_jet->GetName();
  c1->SaveAs(Form("/home/judy/ntuhep/Qstar/runGG_Qstar/graph/%s.png", title));
}
