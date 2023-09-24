//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 24 10:24:27 2023 by ROOT version 6.26/06
// from TTree passedEvents/passedEvents
// found on file: ntuple_VBFHCC_Run32023_174.root
//////////////////////////////////////////////////////////

#ifndef Exam_VBFHcc_h
#define Exam_VBFHcc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class Exam_VBFHcc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       Run;
   ULong64_t       Event;
   ULong64_t       LumiSect;
   Int_t           nVtx;
   Int_t           nInt;
   Int_t           puN;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         BS_x;
   Float_t         BS_y;
   Float_t         BS_z;
   Float_t         BS_xErr;
   Float_t         BS_yErr;
   Float_t         BS_zErr;
   Float_t         BeamWidth_x;
   Float_t         BeamWidth_y;
   Float_t         BeamWidth_xErr;
   Float_t         BeamWidth_yErr;
   Int_t           finalState;
   string          *triggersPassed;
   Bool_t          passedTrig;
   Bool_t          passedZqqSelection;
   vector<string>  *Trigger_l1name;
   vector<int>     *Trigger_l1decision;
   vector<string>  *Trigger_hltname;
   vector<int>     *Trigger_hltdecision;
   vector<int>     *lep_id;
   vector<float>   *lep_pt;
   vector<float>   *lep_eta;
   vector<float>   *lep_phi;
   vector<float>   *lep_mass;
   vector<int>     *Ele_id;
   vector<double>  *Ele_pt;
   vector<bool>    *Ele_isPassID;
   vector<double>  *Ele_eta;
   vector<double>  *Ele_phi;
   vector<double>  *Ele_mass;
   vector<double>  *Ele_dxy;
   vector<double>  *Ele_dz;
   vector<double>  *Ele_hcalIso;
   vector<double>  *Ele_ecalIso;
   vector<double>  *Ele_trackIso;
   vector<bool>    *Ele_isEB;
   vector<double>  *Ele_IsoCal;
   vector<int>     *Muon_id;
   vector<double>  *Muon_pt;
   vector<double>  *Muon_eta;
   vector<double>  *Muon_phi;
   vector<double>  *Muon_PF_Iso_R04;
   vector<double>  *Muon_mass;
   vector<double>  *Muon_dxy;
   vector<double>  *Muon_dz;
   vector<bool>    *Muon_PassLooseID;
   vector<int>     *AK4lep_id;
   vector<double>  *AK4lep_pt;
   vector<double>  *AK4lep_eta;
   vector<double>  *AK4lep_phi;
   vector<double>  *AK4lep_mass;
   Float_t         met;
   Float_t         met_phi;
   Float_t         met_jesup;
   Float_t         met_phi_jesup;
   Float_t         met_jesdn;
   Float_t         met_phi_jesdn;
   Float_t         met_uncenup;
   Float_t         met_phi_uncenup;
   Float_t         met_uncendn;
   Float_t         met_phi_uncendn;
   Int_t           n_jets;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_csv_cTag_vsL;
   vector<float>   *jet_csv_cTag_vsB;
   vector<float>   *jet_mass;
   vector<int>     *jet_isbtag;
   vector<float>   *jet_pfDeepCSVJetTags_probb;
   vector<float>   *jet_pfDeepFlavourJetTags_probbb;
   vector<float>   *jet_pfDeepFlavourJetTags_probc;
   vector<float>   *jet_pfDeepFlavourJetTags_probuds;
   vector<float>   *AK4PuppiJets_pt;
   vector<float>   *AK4PuppiJets_eta;
   vector<float>   *AK4PuppiJets_phi;
   vector<float>   *AK4PuppiJets_mass;
   vector<float>   *jet_pfParticleNetAK4JetTags_probb;
   vector<float>   *jet_pfParticleNetAK4JetTags_probc;
   vector<float>   *jet_pfParticleNetAK4JetTags_probuds;
   vector<float>   *jet_pfParticleNetAK4JetTags_probg;
   vector<float>   *jet_pfParticleNetAK4JetTags_probtauh;
   vector<float>   *jet_pfDeepJetAK4JetTags_probb;
   vector<float>   *jet_pfDeepJetAK4JetTags_probbb;
   vector<float>   *jet_pfDeepJetAK4JetTags_problepb;
   vector<float>   *jet_pfDeepJetAK4JetTags_probc;
   vector<float>   *jet_pfDeepJetAK4JetTags_probuds;
   vector<float>   *jet_pfDeepJetAK4JetTags_probg;
   vector<float>   *jet_pfDeepCSVAK4JetTags_probb;
   vector<float>   *jet_pfDeepCSVAK4JetTags_probbb;
   vector<float>   *jet_pfDeepCSVAK4JetTags_probc;
   vector<float>   *jet_pfDeepCSVAK4JetTags_probudsg;
   Int_t           leadingAK8_pt_idx;
   Int_t           subleadingAK8_pt_idx;
   vector<float>   *AK8PuppiJets_pt;
   vector<float>   *AK8PuppiJets_eta;
   vector<float>   *AK8PuppiJets_phi;
   vector<float>   *AK8PuppiJets_mass;
   vector<double>  *AK8PuppiJets_softdropmass;
   vector<float>   *jet_pfParticleNetJetTags_probZbb;
   vector<float>   *jet_pfParticleNetJetTags_probZcc;
   vector<float>   *jet_pfParticleNetJetTags_probZqq;
   vector<float>   *jet_pfParticleNetJetTags_probQCDbb;
   vector<float>   *jet_pfParticleNetJetTags_probQCDcc;
   vector<float>   *jet_pfParticleNetJetTags_probQCDb;
   vector<float>   *jet_pfParticleNetJetTags_probQCDc;
   vector<float>   *jet_pfParticleNetJetTags_probQCDothers;
   vector<float>   *jet_pfParticleNetJetTags_probHbb;
   vector<float>   *jet_pfParticleNetJetTags_probHcc;
   vector<float>   *jet_pfParticleNetJetTags_probHqqqq;
   vector<float>   *jet_pfMassDecorrelatedParticleNetJetTags_probXbb;
   vector<float>   *jet_pfMassDecorrelatedParticleNetJetTags_probXcc;
   vector<float>   *jet_pfMassDecorrelatedParticleNetJetTags_probXqq;
   vector<float>   *jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb;
   vector<float>   *jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc;
   vector<float>   *jet_pfMassDecorrelatedParticleNetJetTags_probQCDb;
   vector<float>   *jet_pfMassDecorrelatedParticleNetJetTags_probQCDc;
   vector<float>   *jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers;
   vector<float>   *jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb;
   vector<float>   *jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc;
   vector<float>   *jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc;
   vector<float>   *hltjetForBTag_pt;
   vector<float>   *hltjetForBTag_eta;
   vector<float>   *hltjetForBTag_phi;
   vector<float>   *hltjetForBTag_mass;
   vector<float>   *hltjetForBTag_ParticleNet_probb;
   vector<float>   *hltjetForBTag_ParticleNet_probc;
   vector<float>   *hltjetForBTag_ParticleNet_probuds;
   vector<float>   *hltjetForBTag_ParticleNet_probg;
   vector<float>   *hltjetForBTag_ParticleNet_probtauh;
   vector<float>   *hltAK4PFJetsCorrected_pt;
   vector<float>   *hltAK4PFJetsCorrected_eta;
   vector<float>   *hltAK4PFJetsCorrected_phi;
   vector<float>   *hltAK4PFJetsCorrected_mass;
   vector<float>   *L1jet_pt;
   vector<float>   *L1jet_eta;
   vector<float>   *L1jet_phi;
   vector<float>   *L1jet_mass;
   vector<float>   *L1muon_pt;
   vector<float>   *L1muon_eta;
   vector<float>   *L1muon_phi;
   vector<float>   *L1muon_mass;
   vector<int>     *L1muon_qual;
   Float_t         L1ht;
   vector<float>   *quark_pt;
   vector<float>   *quark_eta;
   vector<float>   *quark_phi;
   vector<int>     *quark_flavour;
   vector<bool>    *quark_VBF;
   Float_t         Z_pt;
   Float_t         Z_eta;
   Float_t         Z_phi;
   Float_t         Z_mass;
   Int_t           n_GENjets;
   vector<float>   *GENjet_pt;
   vector<float>   *GENjet_eta;
   vector<float>   *GENjet_phi;
   vector<float>   *GENjet_mass;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_LumiSect;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nInt;   //!
   TBranch        *b_puN;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_BS_x;   //!
   TBranch        *b_BS_y;   //!
   TBranch        *b_BS_z;   //!
   TBranch        *b_BS_xErr;   //!
   TBranch        *b_BS_yErr;   //!
   TBranch        *b_BS_zErr;   //!
   TBranch        *b_BeamWidth_x;   //!
   TBranch        *b_BeamWidth_y;   //!
   TBranch        *b_BeamWidth_xErr;   //!
   TBranch        *b_BeamWidth_yErr;   //!
   TBranch        *b_finalState;   //!
   TBranch        *b_triggersPassed;   //!
   TBranch        *b_passedTrig;   //!
   TBranch        *b_passedZqqSelection;   //!
   TBranch        *b_Trigger_l1name;   //!
   TBranch        *b_Trigger_l1decision;   //!
   TBranch        *b_Trigger_hltname;   //!
   TBranch        *b_Trigger_hltdecision;   //!
   TBranch        *b_lep_id;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_mass;   //!
   TBranch        *b_Ele_id;   //!
   TBranch        *b_Ele_pt;   //!
   TBranch        *b_Ele_isPassID;   //!
   TBranch        *b_Ele_eta;   //!
   TBranch        *b_Ele_phi;   //!
   TBranch        *b_Ele_mass;   //!
   TBranch        *b_Ele_dxy;   //!
   TBranch        *b_Ele_dz;   //!
   TBranch        *b_Ele_hcalIso;   //!
   TBranch        *b_Ele_ecalIso;   //!
   TBranch        *b_Ele_trackIso;   //!
   TBranch        *b_Ele_isEB;   //!
   TBranch        *b_Ele_IsoCal;   //!
   TBranch        *b_Muon_id;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_PF_Iso_R04;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_PassLooseID;   //!
   TBranch        *b_AK4lep_id;   //!
   TBranch        *b_AK4lep_pt;   //!
   TBranch        *b_AK4lep_eta;   //!
   TBranch        *b_AK4lep_phi;   //!
   TBranch        *b_AK4lep_mass;   //!
   TBranch        *b_met;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_jesup;   //!
   TBranch        *b_met_phi_jesup;   //!
   TBranch        *b_met_jesdn;   //!
   TBranch        *b_met_phi_jesdn;   //!
   TBranch        *b_met_uncenup;   //!
   TBranch        *b_met_phi_uncenup;   //!
   TBranch        *b_met_uncendn;   //!
   TBranch        *b_met_phi_uncendn;   //!
   TBranch        *b_n_jets;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_csv_cTag_vsL;   //!
   TBranch        *b_jet_csv_cTag_vsB;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_isbtag;   //!
   TBranch        *b_jet_pfDeepCSVJetTags_probb;   //!
   TBranch        *b_jet_pfDeepFlavourJetTags_probbb;   //!
   TBranch        *b_jet_pfDeepFlavourJetTags_probc;   //!
   TBranch        *b_jet_pfDeepFlavourJetTags_probuds;   //!
   TBranch        *b_AK4PuppiJets_pt;   //!
   TBranch        *b_AK4PuppiJets_eta;   //!
   TBranch        *b_AK4PuppiJets_phi;   //!
   TBranch        *b_AK4PuppiJets_mass;   //!
   TBranch        *b_jet_pfParticleNetAK4JetTags_probb;   //!
   TBranch        *b_jet_pfParticleNetAK4JetTags_probc;   //!
   TBranch        *b_jet_pfParticleNetAK4JetTags_probuds;   //!
   TBranch        *b_jet_pfParticleNetAK4JetTags_probg;   //!
   TBranch        *b_jet_pfParticleNetAK4JetTags_probtauh;   //!
   TBranch        *b_jet_pfDeepJetAK4JetTags_probb;   //!
   TBranch        *b_jet_pfDeepJetAK4JetTags_probbb;   //!
   TBranch        *b_jet_pfDeepJetAK4JetTags_problepb;   //!
   TBranch        *b_jet_pfDeepJetAK4JetTags_probc;   //!
   TBranch        *b_jet_pfDeepJetAK4JetTags_probuds;   //!
   TBranch        *b_jet_pfDeepJetAK4JetTags_probg;   //!
   TBranch        *b_jet_pfDeepCSVAK4JetTags_probb;   //!
   TBranch        *b_jet_pfDeepCSVAK4JetTags_probbb;   //!
   TBranch        *b_jet_pfDeepCSVAK4JetTags_probc;   //!
   TBranch        *b_jet_pfDeepCSVAK4JetTags_probudsg;   //!
   TBranch        *b_leadingAK8_pt_idx;   //!
   TBranch        *b_subleadingAK8_pt_idx;   //!
   TBranch        *b_AK8PuppiJets_pt;   //!
   TBranch        *b_AK8PuppiJets_eta;   //!
   TBranch        *b_AK8PuppiJets_phi;   //!
   TBranch        *b_AK8PuppiJets_mass;   //!
   TBranch        *b_AK8PuppiJets_softdropmass;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probZbb;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probZcc;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probZqq;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probQCDbb;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probQCDcc;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probQCDb;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probQCDc;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probQCDothers;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probHbb;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probHcc;   //!
   TBranch        *b_jet_pfParticleNetJetTags_probHqqqq;   //!
   TBranch        *b_jet_pfMassDecorrelatedParticleNetJetTags_probXbb;   //!
   TBranch        *b_jet_pfMassDecorrelatedParticleNetJetTags_probXcc;   //!
   TBranch        *b_jet_pfMassDecorrelatedParticleNetJetTags_probXqq;   //!
   TBranch        *b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb;   //!
   TBranch        *b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc;   //!
   TBranch        *b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDb;   //!
   TBranch        *b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDc;   //!
   TBranch        *b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers;   //!
   TBranch        *b_jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb;   //!
   TBranch        *b_jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc;   //!
   TBranch        *b_jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc;   //!
   TBranch        *b_hltjetForBTag_pt;   //!
   TBranch        *b_hltjetForBTag_eta;   //!
   TBranch        *b_hltjetForBTag_phi;   //!
   TBranch        *b_hltjetForBTag_mass;   //!
   TBranch        *b_hltjetForBTag_ParticleNet_probb;   //!
   TBranch        *b_hltjetForBTag_ParticleNet_probc;   //!
   TBranch        *b_hltjetForBTag_ParticleNet_probuds;   //!
   TBranch        *b_hltjetForBTag_ParticleNet_probg;   //!
   TBranch        *b_hltjetForBTag_ParticleNet_probtauh;   //!
   TBranch        *b_hltAK4PFJetsCorrected_pt;   //!
   TBranch        *b_hltAK4PFJetsCorrected_eta;   //!
   TBranch        *b_hltAK4PFJetsCorrected_phi;   //!
   TBranch        *b_hltAK4PFJetsCorrected_mass;   //!
   TBranch        *b_L1jet_pt;   //!
   TBranch        *b_L1jet_eta;   //!
   TBranch        *b_L1jet_phi;   //!
   TBranch        *b_L1jet_mass;   //!
   TBranch        *b_L1muon_pt;   //!
   TBranch        *b_L1muon_eta;   //!
   TBranch        *b_L1muon_phi;   //!
   TBranch        *b_L1muon_mass;   //!
   TBranch        *b_L1muon_qual;   //!
   TBranch        *b_L1ht;   //!
   TBranch        *b_quark_pt;   //!
   TBranch        *b_quark_eta;   //!
   TBranch        *b_quark_phi;   //!
   TBranch        *b_quark_flavour;   //!
   TBranch        *b_quark_VBF;   //!
   TBranch        *b_Z_pt;   //!
   TBranch        *b_Z_eta;   //!
   TBranch        *b_Z_phi;   //!
   TBranch        *b_Z_mass;   //!
   TBranch        *b_n_GENjets;   //!
   TBranch        *b_GENjet_pt;   //!
   TBranch        *b_GENjet_eta;   //!
   TBranch        *b_GENjet_phi;   //!
   TBranch        *b_GENjet_mass;   //!

   Exam_VBFHcc(TTree *tree=0);
   virtual ~Exam_VBFHcc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Exam_VBFHcc_cxx
Exam_VBFHcc::Exam_VBFHcc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("VBFHToCC_Run3.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("VBFHToCC_Run3.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("VBFHToCC_Run3.root:/Ana");
      dir->GetObject("passedEvents",tree);

   }
   Init(tree);
}

Exam_VBFHcc::~Exam_VBFHcc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Exam_VBFHcc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Exam_VBFHcc::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Exam_VBFHcc::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   triggersPassed = 0;
   Trigger_l1name = 0;
   Trigger_l1decision = 0;
   Trigger_hltname = 0;
   Trigger_hltdecision = 0;
   lep_id = 0;
   lep_pt = 0;
   lep_eta = 0;
   lep_phi = 0;
   lep_mass = 0;
   Ele_id = 0;
   Ele_pt = 0;
   Ele_isPassID = 0;
   Ele_eta = 0;
   Ele_phi = 0;
   Ele_mass = 0;
   Ele_dxy = 0;
   Ele_dz = 0;
   Ele_hcalIso = 0;
   Ele_ecalIso = 0;
   Ele_trackIso = 0;
   Ele_isEB = 0;
   Ele_IsoCal = 0;
   Muon_id = 0;
   Muon_pt = 0;
   Muon_eta = 0;
   Muon_phi = 0;
   Muon_PF_Iso_R04 = 0;
   Muon_mass = 0;
   Muon_dxy = 0;
   Muon_dz = 0;
   Muon_PassLooseID = 0;
   AK4lep_id = 0;
   AK4lep_pt = 0;
   AK4lep_eta = 0;
   AK4lep_phi = 0;
   AK4lep_mass = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_csv_cTag_vsL = 0;
   jet_csv_cTag_vsB = 0;
   jet_mass = 0;
   jet_isbtag = 0;
   jet_pfDeepCSVJetTags_probb = 0;
   jet_pfDeepFlavourJetTags_probbb = 0;
   jet_pfDeepFlavourJetTags_probc = 0;
   jet_pfDeepFlavourJetTags_probuds = 0;
   AK4PuppiJets_pt = 0;
   AK4PuppiJets_eta = 0;
   AK4PuppiJets_phi = 0;
   AK4PuppiJets_mass = 0;
   jet_pfParticleNetAK4JetTags_probb = 0;
   jet_pfParticleNetAK4JetTags_probc = 0;
   jet_pfParticleNetAK4JetTags_probuds = 0;
   jet_pfParticleNetAK4JetTags_probg = 0;
   jet_pfParticleNetAK4JetTags_probtauh = 0;
   jet_pfDeepJetAK4JetTags_probb = 0;
   jet_pfDeepJetAK4JetTags_probbb = 0;
   jet_pfDeepJetAK4JetTags_problepb = 0;
   jet_pfDeepJetAK4JetTags_probc = 0;
   jet_pfDeepJetAK4JetTags_probuds = 0;
   jet_pfDeepJetAK4JetTags_probg = 0;
   jet_pfDeepCSVAK4JetTags_probb = 0;
   jet_pfDeepCSVAK4JetTags_probbb = 0;
   jet_pfDeepCSVAK4JetTags_probc = 0;
   jet_pfDeepCSVAK4JetTags_probudsg = 0;
   AK8PuppiJets_pt = 0;
   AK8PuppiJets_eta = 0;
   AK8PuppiJets_phi = 0;
   AK8PuppiJets_mass = 0;
   AK8PuppiJets_softdropmass = 0;
   jet_pfParticleNetJetTags_probZbb = 0;
   jet_pfParticleNetJetTags_probZcc = 0;
   jet_pfParticleNetJetTags_probZqq = 0;
   jet_pfParticleNetJetTags_probQCDbb = 0;
   jet_pfParticleNetJetTags_probQCDcc = 0;
   jet_pfParticleNetJetTags_probQCDb = 0;
   jet_pfParticleNetJetTags_probQCDc = 0;
   jet_pfParticleNetJetTags_probQCDothers = 0;
   jet_pfParticleNetJetTags_probHbb = 0;
   jet_pfParticleNetJetTags_probHcc = 0;
   jet_pfParticleNetJetTags_probHqqqq = 0;
   jet_pfMassDecorrelatedParticleNetJetTags_probXbb = 0;
   jet_pfMassDecorrelatedParticleNetJetTags_probXcc = 0;
   jet_pfMassDecorrelatedParticleNetJetTags_probXqq = 0;
   jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb = 0;
   jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc = 0;
   jet_pfMassDecorrelatedParticleNetJetTags_probQCDb = 0;
   jet_pfMassDecorrelatedParticleNetJetTags_probQCDc = 0;
   jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers = 0;
   jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb = 0;
   jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc = 0;
   jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc = 0;
   hltjetForBTag_pt = 0;
   hltjetForBTag_eta = 0;
   hltjetForBTag_phi = 0;
   hltjetForBTag_mass = 0;
   hltjetForBTag_ParticleNet_probb = 0;
   hltjetForBTag_ParticleNet_probc = 0;
   hltjetForBTag_ParticleNet_probuds = 0;
   hltjetForBTag_ParticleNet_probg = 0;
   hltjetForBTag_ParticleNet_probtauh = 0;
   hltAK4PFJetsCorrected_pt = 0;
   hltAK4PFJetsCorrected_eta = 0;
   hltAK4PFJetsCorrected_phi = 0;
   hltAK4PFJetsCorrected_mass = 0;
   L1jet_pt = 0;
   L1jet_eta = 0;
   L1jet_phi = 0;
   L1jet_mass = 0;
   L1muon_pt = 0;
   L1muon_eta = 0;
   L1muon_phi = 0;
   L1muon_mass = 0;
   L1muon_qual = 0;
   quark_pt = 0;
   quark_eta = 0;
   quark_phi = 0;
   quark_flavour = 0;
   quark_VBF = 0;
   GENjet_pt = 0;
   GENjet_eta = 0;
   GENjet_phi = 0;
   GENjet_mass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("LumiSect", &LumiSect, &b_LumiSect);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nInt", &nInt, &b_nInt);
   fChain->SetBranchAddress("puN", &puN, &b_puN);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("BS_x", &BS_x, &b_BS_x);
   fChain->SetBranchAddress("BS_y", &BS_y, &b_BS_y);
   fChain->SetBranchAddress("BS_z", &BS_z, &b_BS_z);
   fChain->SetBranchAddress("BS_xErr", &BS_xErr, &b_BS_xErr);
   fChain->SetBranchAddress("BS_yErr", &BS_yErr, &b_BS_yErr);
   fChain->SetBranchAddress("BS_zErr", &BS_zErr, &b_BS_zErr);
   fChain->SetBranchAddress("BeamWidth_x", &BeamWidth_x, &b_BeamWidth_x);
   fChain->SetBranchAddress("BeamWidth_y", &BeamWidth_y, &b_BeamWidth_y);
   fChain->SetBranchAddress("BeamWidth_xErr", &BeamWidth_xErr, &b_BeamWidth_xErr);
   fChain->SetBranchAddress("BeamWidth_yErr", &BeamWidth_yErr, &b_BeamWidth_yErr);
   fChain->SetBranchAddress("finalState", &finalState, &b_finalState);
   fChain->SetBranchAddress("triggersPassed", &triggersPassed, &b_triggersPassed);
   fChain->SetBranchAddress("passedTrig", &passedTrig, &b_passedTrig);
   fChain->SetBranchAddress("passedZqqSelection", &passedZqqSelection, &b_passedZqqSelection);
   fChain->SetBranchAddress("Trigger_l1name", &Trigger_l1name, &b_Trigger_l1name);
   fChain->SetBranchAddress("Trigger_l1decision", &Trigger_l1decision, &b_Trigger_l1decision);
   fChain->SetBranchAddress("Trigger_hltname", &Trigger_hltname, &b_Trigger_hltname);
   fChain->SetBranchAddress("Trigger_hltdecision", &Trigger_hltdecision, &b_Trigger_hltdecision);
   fChain->SetBranchAddress("lep_id", &lep_id, &b_lep_id);
   fChain->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_mass", &lep_mass, &b_lep_mass);
   fChain->SetBranchAddress("Ele_id", &Ele_id, &b_Ele_id);
   fChain->SetBranchAddress("Ele_pt", &Ele_pt, &b_Ele_pt);
   fChain->SetBranchAddress("Ele_isPassID", &Ele_isPassID, &b_Ele_isPassID);
   fChain->SetBranchAddress("Ele_eta", &Ele_eta, &b_Ele_eta);
   fChain->SetBranchAddress("Ele_phi", &Ele_phi, &b_Ele_phi);
   fChain->SetBranchAddress("Ele_mass", &Ele_mass, &b_Ele_mass);
   fChain->SetBranchAddress("Ele_dxy", &Ele_dxy, &b_Ele_dxy);
   fChain->SetBranchAddress("Ele_dz", &Ele_dz, &b_Ele_dz);
   fChain->SetBranchAddress("Ele_hcalIso", &Ele_hcalIso, &b_Ele_hcalIso);
   fChain->SetBranchAddress("Ele_ecalIso", &Ele_ecalIso, &b_Ele_ecalIso);
   fChain->SetBranchAddress("Ele_trackIso", &Ele_trackIso, &b_Ele_trackIso);
   fChain->SetBranchAddress("Ele_isEB", &Ele_isEB, &b_Ele_isEB);
   fChain->SetBranchAddress("Ele_IsoCal", &Ele_IsoCal, &b_Ele_IsoCal);
   fChain->SetBranchAddress("Muon_id", &Muon_id, &b_Muon_id);
   fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_PF_Iso_R04", &Muon_PF_Iso_R04, &b_Muon_PF_Iso_R04);
   fChain->SetBranchAddress("Muon_mass", &Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_dxy", &Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dz", &Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_PassLooseID", &Muon_PassLooseID, &b_Muon_PassLooseID);
   fChain->SetBranchAddress("AK4lep_id", &AK4lep_id, &b_AK4lep_id);
   fChain->SetBranchAddress("AK4lep_pt", &AK4lep_pt, &b_AK4lep_pt);
   fChain->SetBranchAddress("AK4lep_eta", &AK4lep_eta, &b_AK4lep_eta);
   fChain->SetBranchAddress("AK4lep_phi", &AK4lep_phi, &b_AK4lep_phi);
   fChain->SetBranchAddress("AK4lep_mass", &AK4lep_mass, &b_AK4lep_mass);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_jesup", &met_jesup, &b_met_jesup);
   fChain->SetBranchAddress("met_phi_jesup", &met_phi_jesup, &b_met_phi_jesup);
   fChain->SetBranchAddress("met_jesdn", &met_jesdn, &b_met_jesdn);
   fChain->SetBranchAddress("met_phi_jesdn", &met_phi_jesdn, &b_met_phi_jesdn);
   fChain->SetBranchAddress("met_uncenup", &met_uncenup, &b_met_uncenup);
   fChain->SetBranchAddress("met_phi_uncenup", &met_phi_uncenup, &b_met_phi_uncenup);
   fChain->SetBranchAddress("met_uncendn", &met_uncendn, &b_met_uncendn);
   fChain->SetBranchAddress("met_phi_uncendn", &met_phi_uncendn, &b_met_phi_uncendn);
   fChain->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_csv_cTag_vsL", &jet_csv_cTag_vsL, &b_jet_csv_cTag_vsL);
   fChain->SetBranchAddress("jet_csv_cTag_vsB", &jet_csv_cTag_vsB, &b_jet_csv_cTag_vsB);
   fChain->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_isbtag", &jet_isbtag, &b_jet_isbtag);
   fChain->SetBranchAddress("jet_pfDeepCSVJetTags_probb", &jet_pfDeepCSVJetTags_probb, &b_jet_pfDeepCSVJetTags_probb);
   fChain->SetBranchAddress("jet_pfDeepFlavourJetTags_probbb", &jet_pfDeepFlavourJetTags_probbb, &b_jet_pfDeepFlavourJetTags_probbb);
   fChain->SetBranchAddress("jet_pfDeepFlavourJetTags_probc", &jet_pfDeepFlavourJetTags_probc, &b_jet_pfDeepFlavourJetTags_probc);
   fChain->SetBranchAddress("jet_pfDeepFlavourJetTags_probuds", &jet_pfDeepFlavourJetTags_probuds, &b_jet_pfDeepFlavourJetTags_probuds);
   fChain->SetBranchAddress("AK4PuppiJets_pt", &AK4PuppiJets_pt, &b_AK4PuppiJets_pt);
   fChain->SetBranchAddress("AK4PuppiJets_eta", &AK4PuppiJets_eta, &b_AK4PuppiJets_eta);
   fChain->SetBranchAddress("AK4PuppiJets_phi", &AK4PuppiJets_phi, &b_AK4PuppiJets_phi);
   fChain->SetBranchAddress("AK4PuppiJets_mass", &AK4PuppiJets_mass, &b_AK4PuppiJets_mass);
   fChain->SetBranchAddress("jet_pfParticleNetAK4JetTags_probb", &jet_pfParticleNetAK4JetTags_probb, &b_jet_pfParticleNetAK4JetTags_probb);
   fChain->SetBranchAddress("jet_pfParticleNetAK4JetTags_probc", &jet_pfParticleNetAK4JetTags_probc, &b_jet_pfParticleNetAK4JetTags_probc);
   fChain->SetBranchAddress("jet_pfParticleNetAK4JetTags_probuds", &jet_pfParticleNetAK4JetTags_probuds, &b_jet_pfParticleNetAK4JetTags_probuds);
   fChain->SetBranchAddress("jet_pfParticleNetAK4JetTags_probg", &jet_pfParticleNetAK4JetTags_probg, &b_jet_pfParticleNetAK4JetTags_probg);
   fChain->SetBranchAddress("jet_pfParticleNetAK4JetTags_probtauh", &jet_pfParticleNetAK4JetTags_probtauh, &b_jet_pfParticleNetAK4JetTags_probtauh);
   fChain->SetBranchAddress("jet_pfDeepJetAK4JetTags_probb", &jet_pfDeepJetAK4JetTags_probb, &b_jet_pfDeepJetAK4JetTags_probb);
   fChain->SetBranchAddress("jet_pfDeepJetAK4JetTags_probbb", &jet_pfDeepJetAK4JetTags_probbb, &b_jet_pfDeepJetAK4JetTags_probbb);
   fChain->SetBranchAddress("jet_pfDeepJetAK4JetTags_problepb", &jet_pfDeepJetAK4JetTags_problepb, &b_jet_pfDeepJetAK4JetTags_problepb);
   fChain->SetBranchAddress("jet_pfDeepJetAK4JetTags_probc", &jet_pfDeepJetAK4JetTags_probc, &b_jet_pfDeepJetAK4JetTags_probc);
   fChain->SetBranchAddress("jet_pfDeepJetAK4JetTags_probuds", &jet_pfDeepJetAK4JetTags_probuds, &b_jet_pfDeepJetAK4JetTags_probuds);
   fChain->SetBranchAddress("jet_pfDeepJetAK4JetTags_probg", &jet_pfDeepJetAK4JetTags_probg, &b_jet_pfDeepJetAK4JetTags_probg);
   fChain->SetBranchAddress("jet_pfDeepCSVAK4JetTags_probb", &jet_pfDeepCSVAK4JetTags_probb, &b_jet_pfDeepCSVAK4JetTags_probb);
   fChain->SetBranchAddress("jet_pfDeepCSVAK4JetTags_probbb", &jet_pfDeepCSVAK4JetTags_probbb, &b_jet_pfDeepCSVAK4JetTags_probbb);
   fChain->SetBranchAddress("jet_pfDeepCSVAK4JetTags_probc", &jet_pfDeepCSVAK4JetTags_probc, &b_jet_pfDeepCSVAK4JetTags_probc);
   fChain->SetBranchAddress("jet_pfDeepCSVAK4JetTags_probudsg", &jet_pfDeepCSVAK4JetTags_probudsg, &b_jet_pfDeepCSVAK4JetTags_probudsg);
   fChain->SetBranchAddress("leadingAK8_pt_idx", &leadingAK8_pt_idx, &b_leadingAK8_pt_idx);
   fChain->SetBranchAddress("subleadingAK8_pt_idx", &subleadingAK8_pt_idx, &b_subleadingAK8_pt_idx);
   fChain->SetBranchAddress("AK8PuppiJets_pt", &AK8PuppiJets_pt, &b_AK8PuppiJets_pt);
   fChain->SetBranchAddress("AK8PuppiJets_eta", &AK8PuppiJets_eta, &b_AK8PuppiJets_eta);
   fChain->SetBranchAddress("AK8PuppiJets_phi", &AK8PuppiJets_phi, &b_AK8PuppiJets_phi);
   fChain->SetBranchAddress("AK8PuppiJets_mass", &AK8PuppiJets_mass, &b_AK8PuppiJets_mass);
   fChain->SetBranchAddress("AK8PuppiJets_softdropmass", &AK8PuppiJets_softdropmass, &b_AK8PuppiJets_softdropmass);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probZbb", &jet_pfParticleNetJetTags_probZbb, &b_jet_pfParticleNetJetTags_probZbb);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probZcc", &jet_pfParticleNetJetTags_probZcc, &b_jet_pfParticleNetJetTags_probZcc);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probZqq", &jet_pfParticleNetJetTags_probZqq, &b_jet_pfParticleNetJetTags_probZqq);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probQCDbb", &jet_pfParticleNetJetTags_probQCDbb, &b_jet_pfParticleNetJetTags_probQCDbb);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probQCDcc", &jet_pfParticleNetJetTags_probQCDcc, &b_jet_pfParticleNetJetTags_probQCDcc);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probQCDb", &jet_pfParticleNetJetTags_probQCDb, &b_jet_pfParticleNetJetTags_probQCDb);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probQCDc", &jet_pfParticleNetJetTags_probQCDc, &b_jet_pfParticleNetJetTags_probQCDc);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probQCDothers", &jet_pfParticleNetJetTags_probQCDothers, &b_jet_pfParticleNetJetTags_probQCDothers);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probHbb", &jet_pfParticleNetJetTags_probHbb, &b_jet_pfParticleNetJetTags_probHbb);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probHcc", &jet_pfParticleNetJetTags_probHcc, &b_jet_pfParticleNetJetTags_probHcc);
   fChain->SetBranchAddress("jet_pfParticleNetJetTags_probHqqqq", &jet_pfParticleNetJetTags_probHqqqq, &b_jet_pfParticleNetJetTags_probHqqqq);
   fChain->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probXbb", &jet_pfMassDecorrelatedParticleNetJetTags_probXbb, &b_jet_pfMassDecorrelatedParticleNetJetTags_probXbb);
   fChain->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probXcc", &jet_pfMassDecorrelatedParticleNetJetTags_probXcc, &b_jet_pfMassDecorrelatedParticleNetJetTags_probXcc);
   fChain->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probXqq", &jet_pfMassDecorrelatedParticleNetJetTags_probXqq, &b_jet_pfMassDecorrelatedParticleNetJetTags_probXqq);
   fChain->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb, &b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDbb);
   fChain->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc, &b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDcc);
   fChain->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDb", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDb, &b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDb);
   fChain->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDc", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDc, &b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDc);
   fChain->SetBranchAddress("jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers", &jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers, &b_jet_pfMassDecorrelatedParticleNetJetTags_probQCDothers);
   fChain->SetBranchAddress("jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb", &jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb, &b_jet_pfMassIndependentDeepDoubleBvLV2JetTags_probHbb);
   fChain->SetBranchAddress("jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc", &jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc, &b_jet_pfMassIndependentDeepDoubleCvLV2JetTags_probHcc);
   fChain->SetBranchAddress("jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc", &jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc, &b_jet_pfMassIndependentDeepDoubleCvBV2JetTags_probHcc);
   fChain->SetBranchAddress("hltjetForBTag_pt", &hltjetForBTag_pt, &b_hltjetForBTag_pt);
   fChain->SetBranchAddress("hltjetForBTag_eta", &hltjetForBTag_eta, &b_hltjetForBTag_eta);
   fChain->SetBranchAddress("hltjetForBTag_phi", &hltjetForBTag_phi, &b_hltjetForBTag_phi);
   fChain->SetBranchAddress("hltjetForBTag_mass", &hltjetForBTag_mass, &b_hltjetForBTag_mass);
   fChain->SetBranchAddress("hltjetForBTag_ParticleNet_probb", &hltjetForBTag_ParticleNet_probb, &b_hltjetForBTag_ParticleNet_probb);
   fChain->SetBranchAddress("hltjetForBTag_ParticleNet_probc", &hltjetForBTag_ParticleNet_probc, &b_hltjetForBTag_ParticleNet_probc);
   fChain->SetBranchAddress("hltjetForBTag_ParticleNet_probuds", &hltjetForBTag_ParticleNet_probuds, &b_hltjetForBTag_ParticleNet_probuds);
   fChain->SetBranchAddress("hltjetForBTag_ParticleNet_probg", &hltjetForBTag_ParticleNet_probg, &b_hltjetForBTag_ParticleNet_probg);
   fChain->SetBranchAddress("hltjetForBTag_ParticleNet_probtauh", &hltjetForBTag_ParticleNet_probtauh, &b_hltjetForBTag_ParticleNet_probtauh);
   fChain->SetBranchAddress("hltAK4PFJetsCorrected_pt", &hltAK4PFJetsCorrected_pt, &b_hltAK4PFJetsCorrected_pt);
   fChain->SetBranchAddress("hltAK4PFJetsCorrected_eta", &hltAK4PFJetsCorrected_eta, &b_hltAK4PFJetsCorrected_eta);
   fChain->SetBranchAddress("hltAK4PFJetsCorrected_phi", &hltAK4PFJetsCorrected_phi, &b_hltAK4PFJetsCorrected_phi);
   fChain->SetBranchAddress("hltAK4PFJetsCorrected_mass", &hltAK4PFJetsCorrected_mass, &b_hltAK4PFJetsCorrected_mass);
   fChain->SetBranchAddress("L1jet_pt", &L1jet_pt, &b_L1jet_pt);
   fChain->SetBranchAddress("L1jet_eta", &L1jet_eta, &b_L1jet_eta);
   fChain->SetBranchAddress("L1jet_phi", &L1jet_phi, &b_L1jet_phi);
   fChain->SetBranchAddress("L1jet_mass", &L1jet_mass, &b_L1jet_mass);
   fChain->SetBranchAddress("L1muon_pt", &L1muon_pt, &b_L1muon_pt);
   fChain->SetBranchAddress("L1muon_eta", &L1muon_eta, &b_L1muon_eta);
   fChain->SetBranchAddress("L1muon_phi", &L1muon_phi, &b_L1muon_phi);
   fChain->SetBranchAddress("L1muon_mass", &L1muon_mass, &b_L1muon_mass);
   fChain->SetBranchAddress("L1muon_qual", &L1muon_qual, &b_L1muon_qual);
   fChain->SetBranchAddress("L1ht", &L1ht, &b_L1ht);
   fChain->SetBranchAddress("quark_pt", &quark_pt, &b_quark_pt);
   fChain->SetBranchAddress("quark_eta", &quark_eta, &b_quark_eta);
   fChain->SetBranchAddress("quark_phi", &quark_phi, &b_quark_phi);
   fChain->SetBranchAddress("quark_flavour", &quark_flavour, &b_quark_flavour);
   fChain->SetBranchAddress("quark_VBF", &quark_VBF, &b_quark_VBF);
   fChain->SetBranchAddress("Z_pt", &Z_pt, &b_Z_pt);
   fChain->SetBranchAddress("Z_eta", &Z_eta, &b_Z_eta);
   fChain->SetBranchAddress("Z_phi", &Z_phi, &b_Z_phi);
   fChain->SetBranchAddress("Z_mass", &Z_mass, &b_Z_mass);
   fChain->SetBranchAddress("n_GENjets", &n_GENjets, &b_n_GENjets);
   fChain->SetBranchAddress("GENjet_pt", &GENjet_pt, &b_GENjet_pt);
   fChain->SetBranchAddress("GENjet_eta", &GENjet_eta, &b_GENjet_eta);
   fChain->SetBranchAddress("GENjet_phi", &GENjet_phi, &b_GENjet_phi);
   fChain->SetBranchAddress("GENjet_mass", &GENjet_mass, &b_GENjet_mass);
   Notify();
}

Bool_t Exam_VBFHcc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Exam_VBFHcc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Exam_VBFHcc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Exam_VBFHcc_cxx
