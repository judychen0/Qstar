#include <vector>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
using namespace std;
#include <iostream>
#include <TProfile.h>

#include "untuplizer.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

Double_t deltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
  return dPhi;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta*dEta+dPhi*dPhi);
}



void xQstar(){
  //***********************Initialization***********************//
  
  //access EventTree with TreeReader class
  TreeReader data("/home/judy/ntuhep/Qstar/runGG_Qstar/MC_sample_QstarToGJ_M2000_f1p0_TuneCUETP8M1_13TeV_pyahia8/ggtree_mc_*.root");
  //printf(" reading file %s \n", fname);

  //create an output .root file
  TFile *fout_;
  fout_ = new TFile("output_ggtree_mc.root","RECREATE");

  //create histograms in output .root file
  TH1F *h_npho = new TH1F("h_npho", "n pho", 20, 0., 20);
  TH1F *h_njet = new TH1F("h_njet", "n jet", 20, 0., 20);
  TH1F *h_pjmass = new TH1F("h_pjmass", "pho jet inv mass", 1000, 0., 5000);
  TH1F *h_phoEt = new TH1F("h_phoEt", "pho Et", 100, 0., 1000);

  TH1F *h_nMCpho = new TH1F("h_nMCpho", "n MC pho", 40, 0., 40);
  TH1F *h_nMCele = new TH1F("h_nMCele", "n MC ele", 20, 0., 20); 
  TH1F *h_nMCQu = new TH1F("h_nMCQu", "n MC up quark", 20, 0., 20);
  TH1F *h_nMCQd  = new TH1F("h_nMCQd", "n MC down qurk", 20, 0., 20);
  TH1F *h_nMCG = new TH1F("h_nMCG", "n MC gluon", 20, 0., 20);
  TH1F *h_nMCQstar = new TH1F("h_nMCQstar", "n MC Qstar", 20, 0., 20);
  TH1F *h_MCQstarMass =new TH1F("h_MCQstarMass", "MC Qstar Mass", 1000, 0., 5000);

  TH1F *h_nrealpho = new TH1F("h_nrealpho", "n real photon", 20, 0., 20);
  //TH1F *h_matchR_pho = new TH1F("h_matchR_pho", "ratio of matched photon", 20, 0., 2);
  TH1F *h_phoSigmaIEE = new TH1F("h_phoSigmaIEE", "pho SigmaIEE", 500, 0., 1);
  TH1F *h_dr_pho = new TH1F("h_dr_pho", "dr of photon", 100, 0., 0.2);
  TH1F *h_dpt_pho = new TH1F("h_dpt_pho", "dpt of photon", 100, 0., 1);
  TH2F *h_dptdr_pho = new TH2F("h_dptdr_pho", "dptdr of photon", 100, 0., 1, 100, 0., 2);

  TH1F *h_njetGen = new TH1F("h_njetGen", "n jetGen parton", 20, 0., 20);
  TH1F *h_nrealjet = new TH1F("h_nrealjet", "n real jet", 20, 0., 20);
  TH1F *h_dr_jet = new TH1F("h_dr_jet", "dr of jet", 100, 0., 1);
  TH1F *h_dpt_jet = new TH1F("h_dpt_jet", "dpt of jet", 100, 0., 2);
  TH2F *h_dptdr_jet = new TH2F("h_dptdr_jet", "dptdr of jet", 100, 0., 2, 100, 0., 4);

  TH1F *h_dr_phojet = new TH1F("h_dr_phojet", "dr of phojet", 100, 0., 0.1);
  TH1F *h_pjmass_real = new TH1F("h_pjmass_real", "pho jet real inv mass", 100, 500, 3500);
  TH1F *h_pjmass_HPtpho = new TH1F("h_pjmass_HPtpho", "High Pt pho jet inv mass",100, 500, 3500);
  TH1F *h_pjmass_Sndpho = new TH1F("h_pjmass_Sndpho", "2nd high Pt pho jet inv mass", 100, 500, 3500);
  TH1F *h_pjmass_phoLjet = new TH1F("h_pjmass_phoLjet", "pho Ljet inv mass", 100, 500, 3500);
  TH1F *h_pjmass_phoujet = new TH1F("h_pjmass_phoujet", "pho + ujet inv mass", 100, 500, 3500);
  TH1F *h_pjmass_phodjet = new TH1F("h_pjmass_phodjet", "pho + djet inv mass", 3000, 500, 3500);
  
  h_pjmass->Sumw2();
  h_phoEt->Sumw2();
  h_MCQstarMass->Sumw2();
  h_dr_pho->Sumw2();
  h_dpt_pho->Sumw2();
  h_dr_jet->Sumw2();
  h_dpt_jet->Sumw2();
  h_dr_phojet->Sumw2();
  h_pjmass_real->Sumw2();
  h_pjmass_phoLjet->Sumw2();
  h_pjmass_phoujet->Sumw2();
  h_pjmass_phodjet->Sumw2();
  
  //create ouput tree
  TTree *outtree_;
  outtree_ = new TTree("t", "mini tree");
  
  //define branch variables
  Bool_t   isData;
  Int_t    run;
  Long64_t event;
  Float_t jetPt_, jetEta_, jetPhi_;
  Int_t   isMatched, isMatchedEle, isConverted;
  

  /*
  outtree_->Branch("run", &run, "run/I");
  outtree_->Branch("event", &event, "event/L");
  outtree_->Branch("isData",&isData, "isData/O");
  
  outtree_->Branch("jetPt", &jetPt_, "jetPt/F");
  outtree_->Branch("jetEta", &jetEta_, "jetEta/F");
  outtree_->Branch("jetPhi", &jetPhi_, "jetPhi/F");
  */
  
  /*
  
  
  outtree_->Branch("HLT",         &HLT,        "HLT/L");
  outtree_->Branch("HLTIsPrescaled", &HLTIsPrescaled,        "HLTIsPrescaled/L");
  outtree_->Branch("HLT50ns",         &HLT50ns,        "HLT50ns/L");
  outtree_->Branch("HLTIsPrescaled50ns",         &HLTIsPrescaled50ns,        "HLTIsPrescaled50ns/L");
  outtree_->Branch("phoFiredTrg",  &phoFiredTrgs_,"phoFiredTrgs/I");
  outtree_->Branch("pthat",        &pthat_,       "pthat/F");
  outtree_->Branch("genHT",        &genHT_,       "genHT/F");
  outtree_->Branch("mcPt",         &mcPt_,        "mcPt/F");
  outtree_->Branch("mcEta",        &mcEta_,       "mcEta/F");
  outtree_->Branch("mcPhi",        &mcPhi_,       "mcPhi/F");
  outtree_->Branch("mcCalIso04",   &mcCalIso04_,   "mcCalIso04/F");
  outtree_->Branch("mcTrkIso04",   &mcTrkIso04_,   "mcTrkIso04/F");
  outtree_->Branch("recoPt",       &recoPt,       "recoPt/F");
  outtree_->Branch("recoEta",      &recoEta,      "recoEta/F");
  outtree_->Branch("recoPhi",      &recoPhi,      "recoPhi/F");
  outtree_->Branch("recoSCEta",    &recoSCEta,    "recoSCEta/F");
  outtree_->Branch("r9",           &r9,           "r9/F");
  outtree_->Branch("isMatched",    &isMatched,    "isMatched/I");
  outtree_->Branch("isMatchedEle", &isMatchedEle, "isMatchedEle/I");
  outtree_->Branch("isConverted",    &isConverted,    "isConverted/I");
  outtree_->Branch("idLoose",      &idLoose,      "idLoose/I");
  outtree_->Branch("idMedium",     &idMedium,     "idMedium/I");
  outtree_->Branch("idTight",      &idTight,      "idTight/I");
  outtree_->Branch("nVtx",         &nVtx,         "nVtx/I");
  outtree_->Branch("nPU",          &nPU,          "nPU/I");
  outtree_->Branch("puwei",        &puwei_,        "puwei/F");
  outtree_->Branch("eleVeto",      &eleVeto,      "eleVeto/I");
  outtree_->Branch("HoverE",       &HoverE,       "HoverE/F");
  outtree_->Branch("sieie",        &sieie,        "sieie/F");
  outtree_->Branch("sieip",        &sieip,        "sieip/F");
  outtree_->Branch("sipip",        &sipip,        "sipip/F");
  outtree_->Branch("chIso",        &chIso,        "chIso/F");
  outtree_->Branch("phoIso",       &phoIso,       "phoIso/F");
  outtree_->Branch("nhIso",        &nhIso,        "nhIso/F");
  outtree_->Branch("chIsoRaw",     &chIsoRaw,     "chIsoRaw/F");
  outtree_->Branch("chWorstRaw",   &chWorstIso,   "chWorstIso/F");
  outtree_->Branch("phoIsoRaw",    &phoIsoRaw,    "phoIsoRaw/F");
  outtree_->Branch("nhIsoRaw",     &nhIsoRaw,     "nhIsoRaw/F");
  outtree_->Branch("rho",          &rho,          "rho/F"); 
  outtree_->Branch("e1x3",         &e1x3,         "e1x3/F");
  outtree_->Branch("e2x2",         &e2x2,         "e2x2/F");
  outtree_->Branch("e2x5",         &e2x5,         "e2x5/F");
  outtree_->Branch("e5x5",         &e5x5,         "e5x5/F");
  outtree_->Branch("rawE",         &rawE,         "rawE/F");
  outtree_->Branch("scEtaWidth",   &scEtaWidth,   "scEtaWidth/F");
  outtree_->Branch("scPhiWidth",   &scPhiWidth,   "scPhiWidth/F");
  outtree_->Branch("esRR",         &esRR,         "esRR/F");   
  outtree_->Branch("esEn",         &esEn,         "esEn/F");   
  outtree_->Branch("mva",          &mva,          "mva/F");  
  outtree_->Branch("photonIDmva",       &photonIDmva,       "photonIDmva/F");  
  outtree_->Branch("phoIDbit",          &phoIDbit_,          "phoIDbit/I");  
  outtree_->Branch("mva_hgg",      &mva_hgg,      "mva_hgg/F");  
  outtree_->Branch("HggPresel",    &HggPresel,    "HggPresel/I");  
  outtree_->Branch("Mmm",    &Mmm,    "Mmm/F");  
  outtree_->Branch("Mee",    &Mee,    "Mee/F");  
  outtree_->Branch("MET",    &MET,    "MET/F");  
  outtree_->Branch("metFilters",    &metFilters,    "metFilters/I");  
  outtree_->Branch("METPhi",    &METPhi,    "METPhi/F");  
  outtree_->Branch("phohasPixelSeed", &phohasPixelSeed_, "phohasPixelSeed/I");
  outtree_->Branch("MTm",    &MTm,    "MTm/F");  
  outtree_->Branch("MTe",    &MTe,    "MTe/F");  
  outtree_->Branch("deta_wg",    &deta_wg,    "deta_wg/F");  
  outtree_->Branch("dphi_wg",    &dphi_wg,    "dphi_wg/F");  

  outtree_->Branch("sieieFull5x5",        &sieieFull5x5,        "sieieFull5x5/F");
  outtree_->Branch("sieipFull5x5",        &sieipFull5x5,        "sieipFull5x5/F");
  outtree_->Branch("sipipFull5x5",        &sipipFull5x5,        "sipipFull5x5/F");
  outtree_->Branch("e1x3Full5x5",        &e1x3Full5x5,        "e1x3Full5x5/F");
  outtree_->Branch("r9Full5x5",        &r9Full5x5,        "r9Full5x5/F");
  outtree_->Branch("e2x2Full5x5",        &e2x2Full5x5,        "e2x2Full5x5/F");
  outtree_->Branch("e2x5Full5x5",        &e2x5Full5x5,        "e2x5Full5x5/F");
  outtree_->Branch("e5x5Full5x5",        &e5x5Full5x5,        "e5x5Full5x5/F");

  
  outtree_->Branch("jetY", &jetY_, "jetY/F");
  outtree_->Branch("jetJECUnc", &jetJECUnc_, "jetJECUnc/F");
  outtree_->Branch("jetGenJetPt", &jetGenJetPt_, "jetGenJetPt/F");
  outtree_->Branch("jetGenJetEta", &jetGenJetEta_, "jetGenJetEta/F");
  outtree_->Branch("jetGenJetPhi", &jetGenJetPhi_, "jetGenJetPhi/F");
  outtree_->Branch("jetGenJetY", &jetGenJetY_, "jetGenJetY/F");

  outtree_->Branch("jetCSV2BJetTags",             &jetCSV2BJetTags_, 	      "jetCSV2BJetTags/F"); 	        
  outtree_->Branch("jetDeepCSVTags_b",            &jetDeepCSVTags_b_,         "jetDeepCSVTags_b/F");
  outtree_->Branch("jetDeepCSVTags_bb",           &jetDeepCSVTags_bb_,        "jetDeepCSVTags_bb/F");
  outtree_->Branch("jetDeepCSVTags_c",            &jetDeepCSVTags_c_,         "jetDeepCSVTags_c/F");
  outtree_->Branch("jetDeepCSVTags_udsg",         &jetDeepCSVTags_udsg_,      "jetDeepCSVTags_udsg/F");
 
  outtree_->Branch("jetPartonID", 	          &jetPartonID_, 	      	      "jetPartonID/I"); 	        
  outtree_->Branch("jetGenPartonID", 	          &jetGenPartonID_, 	      	      "jetGenPartonID/I"); 	        
  outtree_->Branch("jetHadFlvr",                  &jetHadFlvr_,                  "jetHadFlvr/I");


  outtree_->Branch("xsweight",  &xsweight, "xsweight/F");
  outtree_->Branch("photon_jetID", &photon_jetID_, "photon_jetID/I");

  outtree_->Branch("SeedTime", &SeedTime_, "SeedTime/F");
  outtree_->Branch("SeedEnergy", &SeedEnergy_, "SeedEnergy/F");
  outtree_->Branch("MIPTotEnergy", &MIPTotEnergy_, "MIPTotEnergy/F");
  */
  
  

  //***********************Loop***********************//

  printf("processing entries %lli \n", data.GetEntriesFast());
  
  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
    //for (Long64_t ev = 0; ev <10000; ev++) {
  //for (Long64_t ev = 0; ev < data.GetEntriesFast()/2.; ev++) {
      
    TLorentzVector phoP4, jetP4, pjP4;

    if (ev % 100000 == 0){
      fprintf(stderr, "Processing event %lli of %lli (%.3f \%)\n", ev+1, data.GetEntriesFast(), (ev+1)*100./data.GetEntriesFast());
    }
    data.GetEntry(ev);
    run     = data.GetInt("run");
    event   = data.GetLong64("event"); 
    isData = data.GetBool("isData");
    
    Int_t nPho     = data.GetInt("nPho");
    Int_t nJet     = data.GetInt("nJet");
    h_npho->Fill(nPho);
    h_njet->Fill(nJet);
    if(nJet <1) continue;

    //reco Qstar
    Float_t* phoE = data.GetPtrFloat("phoE");
    Float_t* phoEt = data.GetPtrFloat("phoEt");
    Float_t* phoEta = data.GetPtrFloat("phoEta");
    Float_t* phoPhi = data.GetPtrFloat("phoPhi");
    Float_t* phoHoverE = data.GetPtrFloat("phoHoverE");
    Float_t* phoSigmaIEE = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
    
    Float_t* jetPt = data.GetPtrFloat("jetPt");
    Float_t* jetEta = data.GetPtrFloat("jetEta");
    Float_t* jetPhi = data.GetPtrFloat("jetPhi");
    Float_t* jetEn = data.GetPtrFloat("jetEn");
    
    for(Int_t ipho=0; ipho < nPho; ++ipho){
      //if(phoEt[ipho]<200.) continue;
      //if(abs(phoEta[ipho])>1.4442) continue;
      phoP4.SetPtEtaPhiM(phoEt[ipho], phoEta[ipho], phoPhi[ipho], 0.);
      h_phoEt->Fill(phoEt[ipho]);
      for(Int_t ijet=0; ijet < nJet; ++ijet){
	if(jetPt[ijet]<170.) continue;
	//if(abs(jetEta[ijet])>2.4) continue;
	//if(abs(phoEta[ipho] - jetEta[ijet])>1.5) continue;
	jetP4.SetPtEtaPhiE(jetPt[ijet], jetEta[ijet], jetPhi[ijet], jetEn[ijet]);
	pjP4 = phoP4 +jetP4;
	h_pjmass->Fill(pjP4.M());
	break;
      }
    }

    //MC Qstar
    Int_t    nMC   =0;     
    Int_t*   mcPID =0;
    Int_t*   mcMomPID =0;
    Float_t* mcPt      =0;
    Float_t* mcEta     =0;
    Float_t* mcPhi     =0;
    Float_t* mcE       =0;
    Float_t* mcEt      =0;
    Int_t*   mcStatus  =0;

    Float_t* mcMomMass =0;
    Float_t* mcMomPt   =0;
    Float_t* mcMomEta  =0;
    Float_t* mcMomPhi  =0;
    Float_t* mcMomE    =0;

    Int_t*   jetGenID    = 0;
    Int_t*   jetGenMomID = 0;
    Float_t* jetGenJetPt    = 0;
    Float_t* jetGenJetEn    = 0;
    Float_t* jetGenJetEta   = 0;
    Float_t* jetGenJetPhi   = 0;

    TLorentzVector mc_phoP4, mc_pQuP4, mc_pQdP4;

    if(!isData){
      nMC       = data.GetInt("nMC");
      mcPID     = data.GetPtrInt("mcPID");
      
      mcPt      = data.GetPtrFloat("mcPt");
      mcEta     = data.GetPtrFloat("mcEta");
      mcPhi     = data.GetPtrFloat("mcPhi");
      mcE       = data.GetPtrFloat("mcE");
      mcStatus  = data.GetPtrInt("mcStatus");

      mcMomPID  = data.GetPtrInt("mcMomPID");
      mcMomMass = data.GetPtrFloat("mcMomMass");
      mcMomPt   = data.GetPtrFloat("mcMomPt");
      mcMomEta  = data.GetPtrFloat("mcMomEta");
      mcMomPhi  = data.GetPtrFloat("mcMomPhi");

      //mcCalIsoDR04 = data.GetPtrFloat("mcCalIsoDR04");

      vector <Int_t> match;
      vector <Int_t> converted;
      vector <Float_t> mmcPt;
      vector <Float_t> mmcEta;
      vector <Float_t> mmcPhi;
      vector <Float_t> mmcCalIso04;

      Int_t nmatch = 0;
      Int_t nconv  = 0;
	
      //get mc Qstar as Mom
      vector <Int_t> mc_Qstarid;
      vector <Float_t> mc_QstarMass;
      Int_t nMCQstar = 0;
      for(Int_t k=0; k < nMC; k++){
	if(mcMomPID[k] > 4000000){
	  mc_Qstarid.push_back(k);
	  h_MCQstarMass->Fill(mcMomMass[k]);
	  nMCQstar++;
	}
      }

      //get mc photon id
      vector <Int_t> mc_phoid;
      Int_t nMCpho = 0;
      for(Int_t k=0; k < nMC; k++){
	if(mcPID[k] == 22 && mcPt[k] >= 165 && mcMomPID[k] >= 22 ){
	  mc_phoid.push_back(k);
	  nMCpho++;
	}
      }

      //get mc ele id
      vector <Int_t> mc_eleid;
      Int_t nMCele = 0;
      for(Int_t k=0; k < nMC; k++){
	if(fabs(mcPID[k]) == 11){
	  mc_eleid.push_back(k);
	  nMCele++;
	}
      }
      
      //create real photon list
      vector <Int_t> realpho_list;
      Int_t nrealpho=0;
      for(Int_t ipho=0; ipho < nPho; ipho++){
	//if(phoEt[ipho] < 165.) continue;
	isMatched     = -1;
	isMatchedEle  = -1;
	isConverted   = -1;

	for(Int_t nn=0; nn < nMCpho; nn++){
	  Int_t k = mc_phoid[nn];
	  Float_t dr = deltaR(phoEta[ipho], phoPhi[ipho], mcEta[k], mcPhi[k]);
	  Float_t dpt = fabs((phoEt[ipho] - mcPt[k])/mcPt[k]);
	  h_dptdr_pho->Fill(dr, dpt);
	  //h_dr_pho->Fill(dr);
	  //if(dr < 0.01)
	  h_dpt_pho->Fill(dpt);
	  //if(dpt < 0.04)
	  h_dr_pho->Fill(dr);
	  if(dr < 0.2 && dpt < 0.2){
	    isMatched = 1;
	    //printf("MC phomatched !");
	    realpho_list.push_back(ipho);
	    nrealpho++;
	    break;
	  }
	}
      }
      if(nPho == 0) continue;
      

      h_nMCpho->Fill(nMCpho);
      h_nMCele->Fill(nMCele);
      h_nMCQstar->Fill(nMCQstar);
      if(nrealpho >=1 )h_nrealpho->Fill(nrealpho);
      //h_matchR_pho->Fill(nrealpho/nPho);
      
      //get genjet id
      jetGenID = data.GetPtrInt("jetGenPartonID");
      jetGenMomID = data.GetPtrInt("jetGenPartonMomID");
      jetGenJetPt = data.GetPtrFloat("jetGenJetPt");
      jetGenJetEn = data.GetPtrFloat("jetGenJetEn");
      jetGenJetEta = data.GetPtrFloat("jetGenJetEta");
      jetGenJetPhi = data.GetPtrFloat("jetGenJetPhi");

      vector <Int_t> jetGenid;
      vector <Int_t> jetGenid_uMom;
      vector <Int_t> jetGenid_dMom;
      Int_t njetGen=0;
      for(Int_t j=0; j < nJet; j++){
	if(jetGenID[j] == 1 || jetGenID[j] == 2){
	  jetGenid.push_back(j);
	  njetGen++;
	}
	/*
	else if(fabs(jetGenMomID[j] == 4000001)){
	  jetGenid_uMom.push_back(j);
	}
	else if(fabs(jetGenMomID[j] == 4000002)){
	  jetGenid_dMom.push_back(j);
	}
	*/
      }
      if(njetGen >=1) h_njetGen->Fill(njetGen);

      //create realjet list
      vector <Int_t> realjet_list;
      Int_t nrealjet=0;
      for(Int_t ijet=0; ijet < nJet; ijet++){
	//if(jetPt[ijet] < 175) continue;
	for(Int_t jj=0; jj < njetGen; jj++){
	  Int_t j = jetGenid[jj];
	  Float_t dr = deltaR(jetEta[ijet], jetPhi[ijet], jetGenJetEta[j], jetGenJetPhi[j]);
	  Float_t dpt = fabs((jetPt[ijet]-jetGenJetPt[j])/jetGenJetPt[j]);
	  h_dptdr_jet->Fill(dpt, dr);
	  if(dpt < 0.1) h_dr_jet->Fill(dr);
	  if(dr < 0.1) h_dpt_jet->Fill(dpt);
	  if(dr < 0.1 && dpt < 0.1){
	    realjet_list.push_back(ijet);
	    nrealjet++;
	    break;
	  }
	}
      }
      //if(nrealjet >= 1)
      h_nrealjet->Fill(nrealjet);

      
      //find a good jet
      if(isMatched == 1){
	TLorentzVector pjP4_tmp;
	TLorentzVector jetuP4, jetdP4;
	//if(nrealpho == 0) continue;
	for(Int_t iipho=0; iipho < nrealpho; iipho++){
	  //if(phoSigmaIEE[realpho_list[iipho]] > 0.01022) continue;
	  if(fabs(phoEta[realpho_list[iipho]]) > 1.4442) continue;
	  phoP4.SetPtEtaPhiM(phoEt[realpho_list[iipho]], phoEta[realpho_list[iipho]], phoPhi[realpho_list[iipho]], 0.);
	  h_phoSigmaIEE->Fill(phoSigmaIEE[realpho_list[iipho]]);
	  for(Int_t iijet=0; iijet < nrealjet; iijet++){
	    //if(nrealjet != 1) continue;
	    jetuP4.SetPtEtaPhiE(jetPt[realjet_list[iijet]], jetEta[realjet_list[iijet]], jetPhi[realjet_list[iijet]], jetEn[realjet_list[iijet]]);
	    pjP4_tmp = phoP4 + jetuP4;
	    h_dr_phojet->Fill(phoP4.DeltaR(jetP4));
	    if(phoP4.DeltaR(jetP4) < 0.4) continue;
	    h_pjmass_real->Fill(pjP4_tmp.M());
	    if(jetGenID[iijet] == 1 || jetGenID[iijet] == 2){
	      h_pjmass_phoLjet->Fill(pjP4_tmp.M());
	    }
	    
	    if(jetGenID[iijet] == 1){
	      h_pjmass_phoujet->Fill(pjP4_tmp.M());
	    }
	    if(jetGenID[iijet] == 2){
	      h_pjmass_phodjet->Fill(pjP4_tmp.M());
	    }
	    if(iipho == 0 && (jetGenID[iijet] == 1 || jetGenID[iijet] == 2)){
	      h_pjmass_HPtpho->Fill(pjP4_tmp.M());
	    }
	    if(iipho == 1 && (jetGenID[iijet] == 1 || jetGenID[iijet] == 2)){
	      h_pjmass_Sndpho->Fill(pjP4_tmp.M());
	    }
	    //jetP4.SetPtEtaPhiM(jetPt[realjet_list[iijet]], jetEta[realjet_list[iijet]], jetPhi[realjet_list[iijet]], jetEn[realjet_list[iijet]]);
	    //if(phoP4.DeltaR(jetP4) < 0.5) continue;
	    //pjP4_tmp = phoP4 + jetP4;

	    /*
	    if(nrealpho == 1 || iipho == 0){
	      h_pjmass_HPTpho->Fill(pjP4_tmp.M());
	    }
	    else if(nrealpho ==2 && iipho ==1){
	      h_pjmass_Sndpho->Fill(pjP4_tmp.M());
	    }
	    */

	  }
      
	}
      }
    }
     
      //check photon property in realpho_list
      
      
     


              
    
    
    //outtree_->Fill();

  }

  //***********************Terminate***********************//
  fout_->cd();
  
  //outtree_->Write();
  h_npho->Write();
  h_njet->Write();
  h_phoEt->Write();
  
  h_pjmass->Write();

  h_nMCpho->Write();
  h_nMCele->Write();
  //h_nMCQu->Write();
  //h_nMCQd->Write();
  //h_nMCG->Write();
  h_nMCQstar->Write();
  h_MCQstarMass->Write();
  h_nrealpho->Write();
  //h_matchR_pho->Write();
  //h_phoSigmaIEE->Write();
  h_dr_pho->Write();
  h_dpt_pho->Write();
  h_dptdr_pho->Write();
  
  h_njetGen->Write();
  h_dr_jet->Write();
  h_dpt_jet->Write();
  h_nrealjet->Write();
  h_dptdr_jet->Write();

  h_dr_phojet->Write();
  h_pjmass_real->Write();
  h_pjmass_phoLjet->Write();
  h_pjmass_phoujet->Write();
  h_pjmass_phodjet->Write();
  h_pjmass_HPtpho->Write();
  h_pjmass_Sndpho->Write();
  
  fout_->Close();
  fprintf(stderr, "Processed all events\n");

}



