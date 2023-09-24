#define Exam_VBFHcc_cxx
#include "Exam_VBFHcc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


/*  Oggetto jet_struct (classe pubblica)
    Attributi: float (CvsAll Value), TLorentzVector (jetRC)
    Metodi: getCvalue (restituisce il valore CvsAll), getJet_p4 (restituisce il TLorentzVector)
*/
struct jet_struct{
  TLorentzVector jet_p4;
  float CvsAll;
  // float CvsL;
  // float CvsB;
  // bool c_tagged;

  jet_struct(TLorentzVector jetRC, float CvsAll_value){
    jet_p4 = jetRC;
    CvsAll = CvsAll_value;
  }

  float getCvalue(){
    return CvsAll;
  }

  TLorentzVector getJet_p4(){
    return jet_p4;
  }
};

// Non so ancora se lo userò
bool compareByCvsAll(const jet_struct &a, const jet_struct &b){
  return a.CvsAll > b.CvsAll;
} 




void Exam_VBFHcc::Loop(){

/*  In a ROOT session, you can do:
      root> .L Exam_VBFHcc.C
      root> Exam_VBFHcc t
      root> t.GetEntry(12); // Fill t data members with entry number 12
      root> t.Show();       // Show values of entry 12
      root> t.Show(16);     // Read and show values of entry 16
      root> t.Loop();       // Loop on all entries

    This is the loop skeleton where:
      jentry is the global entry number in the chain
      ientry is the entry number in the current Tree
    Note that the argument to GetEntry must be:
      jentry for TChain::GetEntry
      ientry for TTree::GetEntry and TBranch::GetEntry

    To read only selected branches, Insert statements like:
      METHOD1:
        fChain->SetBranchStatus("*",0);  // disable all branches
        fChain->SetBranchStatus("branchname",1);  // activate branchname
      METHOD2: replace line
        fChain->GetEntry(jentry);       //read all branches
      by  b_branchname->GetEntry(ientry); //read only this branch
  */




  //Variable declaration

  TLorentzVector GENjet_p4;
  vector<TLorentzVector> cJetHiggs_vec; //vettore di TLorentzVector di jetGEN associati a H->cc
   
  TVector3 GENjet_p3, quark_p3;

  TLorentzVector AK4puppi_p4;
  std::vector<jet_struct> AK4puppi_p4_vec; //vettore di obj-Jet_struct di jetRC

  std::vector<jet_struct> AK4puppi_cHiggs; //vettore di obj-Jet_struct associati a H->cc
  std::vector<jet_struct> AK4puppi_sel;    //vettore di obj-Jet_struct associati a H->cc con treshold 

  double dR_min=1000.;
  double this_dR=0;
  bool matched=false;
  int quark_match_idx=0;
  vector<bool> free_matchQuark;
  float CvsAll_PNet=0., CvsAll_DeepJet=0.;

  int nJetsGEN, nJetsRC;




  //Histograms

  TH1F *histo_HiggsMassGEN = new TH1F("Higgs mass","Higgs Mass GEN",100,0,200);
  TH1F *histo_HiggsMassRC = new TH1F("Higgs mass RC","Higgs Mass RC",100,0,200);
   
  TH1F *histo_pT_cHiggs_jetsGEN = new TH1F("pTcjetsGEN","p_{T} c jets GEN",200,0,250);
  TH1F *histo_eta_cHiggs_jetsGEN = new TH1F("etacjetsGEN","#eta c jets GEN",200,-7,7);
  TH1F *histo_phi_cHiggs_jetsGEN = new TH1F("phicjetsGEN","#phi c jets GEN",200,-4,4);
   
  TH1F *histo_pT_jetsGEN = new TH1F("pTjetsGEN","p_{T} jets GEN",200,0,200);
  TH1F *histo_eta_jetsGEN = new TH1F("etajetsGEN","#eta jets GEN",200,-7,7);
  TH1F *histo_phi_jetsGEN = new TH1F("phijetsGEN","#phi jets GEN",200,-4,4);

  TH1F *histo_pTRC = new TH1F("pT","p_{T} jets RC",200,0,200);

  TH1F *histo_pT_cHiggs_jetsRC = new TH1F("pTcjetsRC","p_{T} c jets RC",200,0,250);
  TH1F *histo_eta_cHiggs_jetsRC = new TH1F("etacjetsRC","#eta c jets RC",200,-7,7);
  TH1F *histo_phi_cHiggs_jetsRC = new TH1F("phicjetsRC","#phi c jets RC",200,-4,4);
   
  TH1F *histo_pT_jetsRC = new TH1F("pTjetsRC","p_{T} jets RC",200,0,200);
  TH1F *histo_eta_jetsRC = new TH1F("etajetsRC","#eta jets RC",200,-7,7);
  TH1F *histo_phi_jetsRC = new TH1F("phijetsRC","#phi jets RC",200,-4,4);

  TH1F *histo_CvsAll_cjetsRC = new TH1F("CvsAll c score","CvsAll RC c jets Score",200,0,1);
  TH1F *histo_CvsAll_jetsRC = new TH1F("CvsAll score","CvsAll RC Score",200,0,1);

  TH2F *histo_etaphi_jetsRC = new TH2F("etaphi", "#eta - #phi other JetsRC", 20,-5,5,20,-3.2,3.2);
  TH2F *histo_CvsAllpT_cHiggs_jetsRC = new TH2F("CvsAllpT", "CvsAll p_{T}", 120,-1.1,1,120,0,150);
  TH3F *histo_etaphipt_jetsRC = new TH3F("etaphipt", "#eta - #phi - #pT other JetsRC",200,-7,7,200,-4,4,200,0,50); 

  TH1F *histo_num_jetsRC = new TH1F("numjetsrc","Number of RC Jets",10,-0.5,20);
  TH1F *histo_num_jetsGEN = new TH1F("numjetsgen","Number of GEN Jets",10,-0.5,20);

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;




  //loop on events
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //for (Long64_t jentry=0; jentry<10000;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

    //reset variables
      cJetHiggs_vec.clear();
      free_matchQuark.clear();
      AK4puppi_sel.clear();
      AK4puppi_cHiggs.clear();
      nJetsGEN = 0;
      nJetsRC = 0;
      cout<<" ******** Evento n="<<jentry<<" ********"<<endl;
      for(unsigned int iquark=0; iquark< quark_pt->size(); iquark++){
        free_matchQuark.push_back(true);
      }


    //loop on GEN jets
      for(unsigned int ijet=0; ijet< GENjet_pt->size(); ijet++){
        nJetsGEN++;
        //Set the TLorentz vector
        GENjet_p4.SetPtEtaPhiM(GENjet_pt->at(ijet),GENjet_eta->at(ijet),GENjet_phi->at(ijet),GENjet_mass->at(ijet));
        //Set the TVector3
        GENjet_p3.SetPtEtaPhi(GENjet_pt->at(ijet),GENjet_eta->at(ijet),GENjet_phi->at(ijet));

        //cout << "......... Jet GEN #" << ijet << " ......... " << endl;
        
        matched=false;

          //in order to select c jets from Higgs
          //we make a geometrical match with quarks at GEN level
          //we require a c quark produced by the Higgs decay to be within dR=0.5 wrt the jet axis

          //loop on quarks
          for(unsigned int iquark=0; iquark< quark_pt->size(); iquark++){
            if(free_matchQuark[iquark]==true){ //check if the quark has been already matched with a jet
              //Set the TVector3 
              quark_p3.SetPtEtaPhi(quark_pt->at(iquark),quark_eta->at(iquark),quark_phi->at(iquark));

              if(quark_VBF->at(iquark)==false){ //condition that selects only c quarks from Higgs decay

                //delta R between the jet and the quark
                this_dR=quark_p3.DeltaR(GENjet_p3);
                dR_min=1000.;
                if(this_dR<dR_min && this_dR<0.5){
                  dR_min=this_dR;
                  quark_match_idx=iquark;
                  matched=true;
                }
              }
            }
          }

          if(matched==true){
            cJetHiggs_vec.push_back(GENjet_p4);
            free_matchQuark[quark_match_idx]=false;
            
            //fill pt, eta and phi of c Jets from Higgs
            //cout << cJetHiggs_vec.size() << endl;

            histo_pT_cHiggs_jetsGEN->Fill(GENjet_p4.Pt());
            histo_eta_cHiggs_jetsGEN->Fill(GENjet_p4.Eta());
            histo_phi_cHiggs_jetsGEN->Fill(GENjet_p4.Phi());
            
            //cout << "Jet GEN #" << ijet << "è c-tag" << endl;

          }

          else{
              
              //cout << "Jet GEN #" << ijet << " NON è c-tag" << endl;

              //Fill histo jets GEN
              histo_pT_jetsGEN->Fill(GENjet_p4.Pt());
              histo_eta_jetsGEN->Fill(GENjet_p4.Eta());
              histo_phi_jetsGEN->Fill(GENjet_p4.Phi());
                  
          }
      }

      //select events with 2 c jets from Higgs 
      if(cJetHiggs_vec.size()==2){
          //fill invariant mass
          histo_HiggsMassGEN->Fill((cJetHiggs_vec[0]+cJetHiggs_vec[1]).M());
      }



    //loop on reco jets
     // we want to identify among all the reconstructed jets, the ones produced by the Higgs decay
     // the ParticleNet tagging algorithm assignd to each jet the probability to be b,c,light,...
     // we select the pair of jets with the largest value of CvsAll
     // CvsAll is the ratio of the probability of being a c jet over the probability of being any flavour
     // Si è trovato che il valore di threshold può essere tra i 0.2 - 0.3

    for(unsigned int ijet=0; ijet< AK4PuppiJets_pt->size(); ijet++){
      nJetsRC++;   
      //cout << "+++++++++++ Jet RC #" << ijet << " +++++++++++++" << endl;
         
      AK4puppi_p4.SetPtEtaPhiM(AK4PuppiJets_pt->at(ijet), AK4PuppiJets_eta->at(ijet), AK4PuppiJets_phi->at(ijet), AK4PuppiJets_mass->at(ijet));
      float num= jet_pfParticleNetAK4JetTags_probc->at(ijet);
      float den= jet_pfParticleNetAK4JetTags_probc->at(ijet)+ jet_pfParticleNetAK4JetTags_probb->at(ijet) + jet_pfParticleNetAK4JetTags_probuds->at(ijet) + jet_pfParticleNetAK4JetTags_probg->at(ijet);
      float CvsAll_val= (den!=0 || (num/den)<999)? num/den : -1; //Se il valore non è valido inizializza CvsAll a -1
         
      jet_struct *jetRC = new jet_struct(AK4puppi_p4, CvsAll_val); //Jet con TLorentzVect + cTag

      //histo_pTRC->Fill(AK4puppi_p4.Pt());

      //matching of the reco-jet with the gen-jets from Higgs
      for(int i=0; i<cJetHiggs_vec.size(); i++ ){
        //cout<<" cJetHiggs_vec size: "<<cJetHiggs_vec.size() << endl;
        double DR_genReco = cJetHiggs_vec[i].DeltaR(AK4puppi_p4);
        //cout<<" DR="<<DR_genReco << endl;
             
        if(DR_genReco < 0.5){ //Qui il problema è che AK4puppi_p4 può soddisfare due volte questa condizione
              
          AK4puppi_cHiggs.push_back(*jetRC); //Ho un vettore di jetRC associati a Higgs
              
          //cout << "Jet RC #" << ijet << " è c-tag" << endl;
          //histo for reco jet from higgs
            
          histo_pT_cHiggs_jetsRC->Fill((AK4puppi_cHiggs[i].getJet_p4()).Pt());
          //histo_eta_cHiggs_jetsRC->Fill((AK4puppi_cHiggs[ijet].getJet_p4()).Eta());
          histo_eta_cHiggs_jetsRC->Fill(AK4puppi_p4.Eta());
          histo_phi_cHiggs_jetsRC->Fill(AK4puppi_p4.Phi());
          //histo_pT_cHiggs_jetsRC->Fill(AK4puppi_p4.Pt());

          histo_CvsAllpT_cHiggs_jetsRC->Fill(jetRC->getCvalue(),AK4puppi_p4.Pt());
          histo_pTRC->Fill(AK4puppi_p4.Pt());




              
          if(jetRC->getCvalue()>0){
            histo_CvsAll_cjetsRC->Fill(jetRC->getCvalue());
          }

          if(jetRC->getCvalue()>0.3){
            AK4puppi_sel.push_back(*jetRC);
          }

          //cout << "pt_cHiggs_jetsRC= " << AK4puppi_cHiggs[ijet].getJet_p4().Pt() << endl;  

        }

        else{
          histo_pTRC->Fill(AK4puppi_p4.Pt());    
          AK4puppi_p4_vec.push_back(*jetRC); //Ho un vettore dei restanti jetRC
              
          //cout << "Jet RC #" << ijet << " NON è c-tag" << endl;

          //histo for any other reco jet
          histo_pT_jetsRC->Fill(AK4puppi_p4.Pt());
          histo_eta_jetsRC->Fill(AK4puppi_p4.Eta());
          histo_phi_jetsRC->Fill(AK4puppi_p4.Phi());

          histo_etaphi_jetsRC->Fill(AK4puppi_p4.Eta(), AK4puppi_p4.Phi());
          histo_etaphipt_jetsRC->Fill(AK4puppi_p4.Eta(), AK4puppi_p4.Phi(), AK4puppi_p4.Pt()); 
              
          if(jetRC->getCvalue()>0){
            histo_CvsAll_jetsRC->Fill(jetRC->getCvalue());
          }    
        }        
      }

      //select events with 2 c jets from Higgs 
      if(AK4puppi_sel.size()==2){
        //fill invariant mass
        histo_HiggsMassRC->Fill((AK4puppi_sel[0].getJet_p4() + AK4puppi_sel[1].getJet_p4()).M());
      } 
    }          
  
    histo_num_jetsRC->Fill(nJetsRC);
    histo_num_jetsGEN->Fill(nJetsGEN);

  }




// Plot Section
  
  // HIGGS MASS GENERATE & RECONSTRUCTED
    
      gStyle->SetPadTickY(1);
      gStyle->SetPadTickX(1);
      gStyle->SetOptStat("nemruo");
      
      TCanvas * c_HiggsMassGEN = new TCanvas("Higgs_massGEN", "Higgs mass GEN",600, 600);
        c_HiggsMassGEN->cd();
        
        histo_HiggsMassGEN->SetTitle("Higgs mass GEN; m_{H} [GeV];  N_{Particles}");
        histo_HiggsMassGEN->SetLineWidth(1);
        histo_HiggsMassGEN->SetLineColor(kRed+1);
        histo_HiggsMassGEN->SetFillColor(kRed-4);
        histo_HiggsMassGEN->SetFillStyle(3344);
        histo_HiggsMassGEN->Draw();


      TCanvas * c_HiggsMassRC = new TCanvas("Higgs_massRC", "Higgs mass RC",600, 600);
        c_HiggsMassRC->cd();
        
        auto l_HiggsMassRC = new TLegend();
        l_HiggsMassRC->AddEntry(histo_HiggsMassRC,"m_{H} recon","f");
        l_HiggsMassRC->SetBorderSize(0);
        
        histo_HiggsMassRC->SetTitle("Higgs mass RC; m_{H} [GeV]; N_{Particles}");
        histo_HiggsMassRC->SetLineWidth(1);
        histo_HiggsMassRC->SetLineColor(kAzure);
        histo_HiggsMassRC->SetFillColor(kAzure-4);
        histo_HiggsMassRC->SetFillStyle(3325);

        histo_HiggsMassRC->Draw();
        l_HiggsMassRC->Draw();

      
      TCanvas * c_HiggsMass = new TCanvas("Higgs_mass", "Higgs mass",600,800);
        c_HiggsMass->cd();
        
        auto legend_higgsmass = new TLegend(0.1,0.7,0.48,0.9);
        legend_higgsmass->AddEntry(histo_HiggsMassRC,"Higgs mass RC","f");
        legend_higgsmass->AddEntry(histo_HiggsMassGEN,"Higgs mass GEN","f");
        legend_higgsmass->SetBorderSize(0);
        
        histo_HiggsMassRC->SetLineWidth(1);
        
        histo_HiggsMassGEN->Draw();
        histo_HiggsMassRC->Draw("same");
        legend_higgsmass->Draw();

  /*
  // GENERATE LEVEL  
    
    //Draw other-JetsGEN variables
      TCanvas * c_otherJets_pT_GEN = new TCanvas("pt_GENjet", "pt GEN jet",600, 600);
        c_otherJets_pT_GEN->cd();
        
        histo_pT_jetsGEN->SetTitle("p_{T} GEN jets; p_{T} [GeV]; # jets");
        histo_pT_jetsGEN->Draw();

      
      TCanvas * c_otherJets_eta_GEN = new TCanvas("eta_GENjet", "#eta GEN jet",600, 600);
        c_otherJets_eta_GEN->cd();
        
        histo_eta_jetsGEN->SetTitle("#eta GEN jets; ; # jets");
        histo_eta_jetsGEN->Draw();

      
      TCanvas * c_otherJets_phi_GEN = new TCanvas("phi_GENjet", "#phi GEN jet",600, 600);
        c_otherJets_phi_GEN->cd();
        
        histo_phi_jetsGEN->SetTitle("phi GEN jets; phi Jet; # jets");
        histo_phi_jetsGEN->Draw();


    //Draw c-JetsGEN variables
      TCanvas * c_cJets_pT_GEN = new TCanvas("pt_GENcjet", "pt GEN c jet",600, 600);
        c_cJets_pT_GEN->cd();
        
        histo_pT_cHiggs_jetsGEN->SetTitle("p_{T} GEN c jets; p_{T} [GeV]; # jets");
        histo_pT_cHiggs_jetsGEN->Draw();


      TCanvas * c_cJets_eta_GEN = new TCanvas("eta_GENcjet", "#eta GEN c jet",600, 600);
        c_cJets_eta_GEN->cd();
        
        histo_eta_cHiggs_jetsGEN->SetTitle("#eta GEN c jets; ; # jets");
        histo_eta_cHiggs_jetsGEN->Draw();


      TCanvas * c_cJets_phi_GEN = new TCanvas("phi_GENcjet", "#phi GEN c jet",600, 600);
        c_cJets_phi_GEN->cd();
        
        histo_phi_cHiggs_jetsGEN->SetTitle("phi GEN c jets; phi Jet; # jets");
        histo_phi_cHiggs_jetsGEN->Draw();



  // RECONSTRUCTED LEVEL
  
    //Draw other-JetsRC variables
      TCanvas * c_otherJets_pT_RC = new TCanvas("pt_RCjet", "pt RC jet",600, 600);
        c_otherJets_pT_RC->cd();

        histo_pT_jetsRC->SetTitle("p_{T} RC jets; p_{T} [GeV]; # jets");
        histo_pT_jetsRC->Draw();
     
      
      TCanvas * c_otherJets_eta_RC = new TCanvas("eta_RCjet", "#eta RC jet",600, 600);
        c_otherJets_eta_RC->cd();

        histo_eta_jetsRC->SetTitle("#eta RC jets; ; # jets");
        histo_eta_jetsRC->Draw();

      
      TCanvas * c_otherJets_phi_RC = new TCanvas("phi_RCjet", "#phi RC jet",600, 600);
        c_otherJets_phi_RC->cd();
        
        histo_phi_jetsRC->SetTitle("phi RC jets; phi Jet; # jets");
        histo_phi_jetsRC->Draw();


    //Draw c-JetsRC variables
      TCanvas * c_cJets_pT_RC = new TCanvas("pt_RCcjet", "pt RC c jet",600, 600);
        c_cJets_pT_RC->cd();
        histo_pT_cHiggs_jetsRC->SetTitle("p_{T} RC c jets; p_{T} [GeV]; # jets");
        histo_pT_cHiggs_jetsRC->Draw(); 

      TCanvas * c_cJets_eta_RC = new TCanvas("eta_RCcjet", "#eta RC c jet",600, 600);
        c_cJets_eta_RC->cd();
        histo_eta_cHiggs_jetsRC->SetTitle("#eta RC c jets; ; # jets");
        histo_eta_cHiggs_jetsRC->Draw();

      TCanvas * c_cJets_phi_RC = new TCanvas("phi_RCjet", "#phi RC c jet",600, 600);
        c_cJets_phi_RC->cd();
        
        histo_phi_cHiggs_jetsRC->SetTitle("phi RC c jets; phi c Jet; # jets");
        histo_phi_cHiggs_jetsRC->Draw();  
  */



  // COMBINED PLOTS 

    // Transvers Momentum p_T
      TCanvas * c_pT = new TCanvas("pT", "p_{T}",1080, 1080);
      c_pT->Divide(2,2);

      c_pT->cd(1);
        histo_pT_jetsGEN->SetTitle("p_{T} other Jets Generated; p_{T} [GeV]; N_{Jets}");
        histo_pT_jetsGEN->SetLineColor(kCyan+3);
        histo_pT_jetsGEN->SetFillColor(kCyan+1);
        histo_pT_jetsGEN->SetLineWidth(1);
        histo_pT_jetsGEN->SetFillStyle(3544);
        histo_pT_jetsGEN->Draw();

      c_pT->cd(2);
        histo_pT_cHiggs_jetsGEN->SetTitle("p_{T}  c-Jets Generated; p_{T} [GeV]; N_{Jets}");
        histo_pT_cHiggs_jetsGEN->SetLineColor(kOrange+9);
        histo_pT_cHiggs_jetsGEN->SetFillColor(kOrange-3);
        histo_pT_cHiggs_jetsGEN->SetLineWidth(1);
        histo_pT_cHiggs_jetsGEN->SetFillStyle(3544);
        histo_pT_cHiggs_jetsGEN->Draw();

      c_pT->cd(3);
        histo_pT_jetsRC->SetTitle("p_{T} other Jets Reconstructed; p_{T} [GeV]; N_{Jets}");
        histo_pT_jetsRC->SetLineColor(kAzure);
        histo_pT_jetsRC->SetFillColor(kAzure-4);
        histo_pT_jetsRC->SetLineWidth(1);
        histo_pT_jetsRC->SetFillStyle(3544);
        histo_pT_jetsRC->Draw();

      c_pT->cd(4);
        histo_pT_cHiggs_jetsRC->SetTitle("p_{T}  c-Jets Reconstructed; p_{T} [GeV]; N_{Jets}");
        histo_pT_cHiggs_jetsRC->SetLineColor(kRed+1);
        histo_pT_cHiggs_jetsRC->SetFillColor(kRed-4);
        histo_pT_cHiggs_jetsRC->SetLineWidth(1);
        histo_pT_cHiggs_jetsRC->SetFillStyle(3544);
        histo_pT_cHiggs_jetsRC->Draw();

    // Pseudorapidity eta
      TCanvas * c_eta = new TCanvas("#eta", "#eta",1080, 1080);
      c_eta->Divide(2,2);

      c_eta->cd(1);
        histo_eta_jetsGEN->SetTitle("#eta other Jets Generated; #eta; N_{Jets}");
        histo_eta_jetsGEN->SetLineColor(kCyan+3);
        histo_eta_jetsGEN->SetFillColor(kCyan+1);
        histo_eta_jetsGEN->SetLineWidth(1);
        histo_eta_jetsGEN->SetFillStyle(3544);
        histo_eta_jetsGEN->Draw();

      c_eta->cd(2);
        histo_eta_cHiggs_jetsGEN->SetTitle("#eta  c-Jets Generated; #eta; N_{Jets}");
        histo_eta_cHiggs_jetsGEN->SetLineColor(kOrange+9);
        histo_eta_cHiggs_jetsGEN->SetFillColor(kOrange-3);
        histo_eta_cHiggs_jetsGEN->SetLineWidth(1);
        histo_eta_cHiggs_jetsGEN->SetFillStyle(3544);
        histo_eta_cHiggs_jetsGEN->Draw();

      c_eta->cd(3);
        histo_eta_jetsRC->SetTitle("#eta other Jets Reconstructed; #eta; N_{Jets}");
        histo_eta_jetsRC->SetLineColor(kAzure);
        histo_eta_jetsRC->SetFillColor(kAzure-4);
        histo_eta_jetsRC->SetLineWidth(1);
        histo_eta_jetsRC->SetFillStyle(3544);
        histo_eta_jetsRC->Draw();

      c_eta->cd(4);
        histo_eta_cHiggs_jetsRC->SetTitle("#eta  c-Jets Reconstructed; #eta; N_{Jets}");
        histo_eta_cHiggs_jetsRC->SetLineColor(kRed+1);
        histo_eta_cHiggs_jetsRC->SetFillColor(kRed-4);
        histo_eta_cHiggs_jetsRC->SetLineWidth(1);
        histo_eta_cHiggs_jetsRC->SetFillStyle(3544);
        histo_eta_cHiggs_jetsRC->Draw();

    // Azimuthal Angle phi  
      TCanvas * c_phi = new TCanvas("phi", "#phi",1080, 1080);
      c_phi->Divide(2,2);

      c_phi->cd(1);
      histo_phi_jetsGEN->SetTitle("#phi other Jets Generated; #phi [rad]; N_{Jets}");
      histo_phi_jetsGEN->SetLineColor(kCyan+3);
      histo_phi_jetsGEN->SetFillColor(kCyan+1);
      histo_phi_jetsGEN->SetLineWidth(1);
      histo_phi_jetsGEN->SetFillStyle(3544);
      histo_phi_jetsGEN->Draw();

      c_phi->cd(2);
      histo_phi_cHiggs_jetsGEN->SetTitle("#phi  c-Jets Generated; #phi [rad]; N_{Jets}");
      histo_phi_cHiggs_jetsGEN->SetLineColor(kOrange+9);
      histo_phi_cHiggs_jetsGEN->SetFillColor(kOrange-3);
      histo_phi_cHiggs_jetsGEN->SetLineWidth(1);
      histo_phi_cHiggs_jetsGEN->SetFillStyle(3544);
      histo_phi_cHiggs_jetsGEN->Draw();

      c_phi->cd(3);
      histo_phi_jetsRC->SetTitle("#phi other Jets Reconstructed; #phi [rad]; N_{Jets}");
      histo_phi_jetsRC->SetLineColor(kAzure);
      histo_phi_jetsRC->SetFillColor(kAzure-4);
      histo_phi_jetsRC->SetLineWidth(1);
      histo_phi_jetsRC->SetFillStyle(3544);
      histo_phi_jetsRC->Draw();

      c_phi->cd(4);
      histo_phi_cHiggs_jetsRC->SetTitle("#phi  c-Jets Reconstructed; #phi [rad]; N_{Jets}");
      histo_phi_cHiggs_jetsRC->SetLineColor(kRed+1);
      histo_phi_cHiggs_jetsRC->SetFillColor(kRed-4);
      histo_phi_cHiggs_jetsRC->SetLineWidth(1);
      histo_phi_cHiggs_jetsRC->SetFillStyle(3544);
      histo_phi_cHiggs_jetsRC->Draw();




  // 2D-HISTO
    
      // Eta vs Phi 
      TCanvas * c_etaphi = new TCanvas("etaphi", "#eta #phi",1080, 1080);
      c_etaphi->cd();
      gStyle->SetPalette(kTemperatureMap);

      histo_etaphi_jetsRC->SetTitle("#eta - #phi; #eta; #phi [rad]; N_{Jets}");
      histo_etaphi_jetsRC->SetContour(300);
      histo_etaphi_jetsRC->Draw("colz"); //LEGO2Z
      //auto legend2d = new TLegend(0.1,0.7,0.48,0.9);
      //legend2d->Draw();
      

      // CvsAll Score vs pT
      TCanvas * c_CvsAllpT = new TCanvas("CvsAllpT", "CvsAll p_{T}",1080, 1080);
      c_CvsAllpT->cd();
      gStyle->SetPalette(kTemperatureMap);

      histo_CvsAllpT_cHiggs_jetsRC->SetTitle("CvsAll p_{T}; CvsAll; p_{T} [GeV]; N_{Jets}");
      histo_CvsAllpT_cHiggs_jetsRC->SetContour(300);
      histo_CvsAllpT_cHiggs_jetsRC->Draw("colz");
    

  /*
  // 3D-HISTO
      TCanvas * c_etaphipt = new TCanvas("etaphipt", "#eta #phi p_{T}",1080, 1080);
      c_etaphipt->cd();

      histo_etaphipt_jetsRC->SetTitle("#eta - #phi - p_{T}; #eta; #phi [rad]; p_{T}");
      histo_etaphipt_jetsRC->Draw("lego");
  */


  
  //CvsAll score
    TCanvas * cscore = new TCanvas("cscore", "CvsAll RCjets Score",600, 600);
    cscore->Divide(1,2);
    cscore->cd(1);
    histo_CvsAll_cjetsRC->SetTitle("CvsAll score RC c jets; Score Value; # jets");
    histo_CvsAll_cjetsRC->SetLineColor(kRed);
    histo_CvsAll_cjetsRC->Draw();
    cscore->cd(2);
    histo_CvsAll_jetsRC->SetLineColor(kBlue);
    histo_CvsAll_jetsRC->SetTitle("CvsAll score RC no c jets; Score Value; # jets");
    histo_CvsAll_jetsRC->Draw();

    
    TCanvas * c_ctag = new TCanvas("scoretagRCcjet", "CvsAll RC jets Score",600, 600);
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetOptStat("nemruo");
    c_ctag->cd();
    
    histo_CvsAll_jetsRC->SetTitle("CvsAll score RC no c jets; Score Value; N_{Jets}");
    histo_CvsAll_jetsRC->SetLineWidth(1);
    histo_CvsAll_jetsRC->SetLineColor(kAzure);
    histo_CvsAll_jetsRC->SetFillColor(kAzure-4);
    histo_CvsAll_jetsRC->SetFillStyle(3325);
    histo_CvsAll_jetsRC->Draw();
    
    histo_CvsAll_cjetsRC->SetTitle("CvsAll score RC c jets; Score Value; # jets");
    histo_CvsAll_cjetsRC->SetLineWidth(1);
    histo_CvsAll_cjetsRC->SetLineColor(kRed+1);
    histo_CvsAll_cjetsRC->SetFillColor(kRed-4);
    histo_CvsAll_cjetsRC->SetFillStyle(3344);
    histo_CvsAll_cjetsRC->Draw("Same");
    
    legend->AddEntry(histo_CvsAll_cjetsRC,"Jets from cc","f");
    legend->AddEntry(histo_CvsAll_jetsRC,"All jetsRC","f");
    legend->SetBorderSize(0);
    legend->Draw();
  
    TCanvas *c1 = new TCanvas("jetsRc", "Jets pt RC", 600,600);
    c1->cd();
    histo_pTRC->Draw();

  // Numero di Jets
    TCanvas * numjets = new TCanvas("numJets", "Num jets per evento",600, 600);
    auto leg = new TLegend();
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetOptStat("nemruo");
    numjets->cd();
    
    histo_num_jetsGEN->SetTitle("Num jets per evento GEN; num jets; N_{evento}");
    histo_num_jetsGEN->SetLineWidth(1);
    histo_num_jetsGEN->SetLineColor(kRed+1);
    histo_num_jetsGEN->Draw();
    
    histo_num_jetsRC->SetTitle("Num jets per evento; num jets; N_{evento}");
    histo_num_jetsRC->SetLineWidth(1);
    histo_num_jetsRC->SetLineColor(kAzure);
    histo_num_jetsRC->Draw("Same");
    
    leg->AddEntry(histo_num_jetsGEN,"Numero jet gen","l");
    leg->AddEntry(histo_num_jetsRC,"Numero jet rico","l");
    leg->SetBorderSize(0);
    leg->Draw();


}
