# Vector Boson Fusion H→cc search at the CMS experiment

The goal of this exercise is to study the kinematic variables of the jets produced by the Higgs boson decaying into c-quarks and prepare a proposal on how to perform this search at the CMS experiment. We assume that the production mechanism of the Higgs boson is the Vector Boson Fusion.
We prepared a rootfile of events generated with Madgraph that underwent the transportation in GEANT4 and then the event reconstruction in the CMS experiment. 
The rootfile can be downloaded at: [https://cernbox.cern.ch/s/nFTu1Dj50RdqFBW](https://cernbox.cern.ch/s/nFTu1Dj50RdqFBW)

The file contains several variables, but the relevant ones to be used are:
- For generator level jets: _GENjet_pt, GENjet_eta, GENjet_phi, GENjet_mass_ 
- For reconstructed jets: _AK4PuppiJets_pt, AK4PuppiJets_eta, AK4PuppiJets_phi,AK4PuppiJets_mass,  jet_pfParticleNetAK4JetTags_probc,  jet_pfParticleNetAK4JetTags_probb,  jet_pfParticleNetAK4JetTags_probuds,  jet_pfParticleNetAK4JetTags_probg_ 

We prepared a macro ( Exam_VBFHcc.C ) that can run on the events stored in the file, loop on the  gen-leve jets and flag the ones coming from the Higgs boson. 

### How to run the macro

Open a command prompt, then in a root session:

- root [0] .L Exam_VBFHcc.C
- root [1] Exam_VBFHcc t
- root [2] t.Loop()

The report and the relative documentation can be found in [here](VBFHcc_MarcoCecca.pdf).
