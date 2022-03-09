// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "MathUtils/Utils.h"

#include <fastjet/config.h>
#include <fastjet/AreaDefinition.hh>
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

using namespace o2;
using namespace o2::framework;
using namespace fastjet;

using Particles = aod::McParticles;

//the task analyseMFTJets loops over MFT tracks and generated particles and fills basic histograms

struct analyseMFTJets {
  int icoll = 0;
  Service<TDatabasePDG> pdg;
  double ghost_maxrap = 5;
  Configurable<float> RConf{"RConf", 0.7, "jetR for the jet definition"};

  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;
  float jetR=RConf;

  fastjet::JetDefinition *jet_defGen = new fastjet::JetDefinition(fastjet::cambridge_algorithm, jetR, recombScheme, strategy); //
  fastjet::JetDefinition *jet_defRec = new fastjet::JetDefinition(fastjet::cambridge_algorithm, jetR, recombScheme, strategy); //

  fastjet::GhostedAreaSpec *ghostareaspec = new fastjet::GhostedAreaSpec(5., 1, 0.05); //ghost
  //GhostedAreaSpec ghostareaspec(5., 1, 0.01); //ghost of maxrap=5, 1 repetition, small ghost area
  fastjet::AreaType areaType = fastjet::active_area;
  fastjet::AreaDefinition *areaDef = new fastjet::AreaDefinition(areaType, *ghostareaspec);

  HistogramRegistry registry{
    "registry",
    {
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"NTracksOverNparts", "Number of reconstructed tracks / number of particles in a collision; Ntracks/NParts; #count", {HistType::kTH1F, {{300, 0, 150}}}}, //
      {"NTracksPerJet", "Number of reconstructed tracks per jet; Ntracks/jet; #count", {HistType::kTH1F, {{30, 1, 30}}}}, //
      {"NPartsPerJet", "Number of particles per jet; Nparts/jet; #count", {HistType::kTH1F, {{30, 1, 30}}}}, //
      {"AreaMCjet", "area of jet in gen MC; area; #count", {HistType::kTH1F, {{300, 0, 30}}}}, //
      {"AreaRecojet", "area of jet in reco MC; area; #count", {HistType::kTH1F, {{300, 0, 30}}}}, //
      {"AreaMCjetVsPt", "area of jet in gen MC vs #p_T^{jet}; area; pt; #count", {HistType::kTH2F, {{300, 0, 30}, {200, 0, 10}}}}, //
      {"NJetGenOverNJetTrue", "Number of generated jet / Number of true jets; ratio; #count", {HistType::kTH1F, {{800, 0, 4}}}}, //
      {"NJetGenOverNJetTrue2", "Number of generated jet / Number of true jets; ratio; #count", {HistType::kTH1F, {{800, 0, 4}}}},
      {"EtaRecEtaGen", "; #eta_{Rec}; #eta_{Gen}; tracks", {HistType::kTH2F, {{35, -4.5, -1.}, {35, -4.5, -1.}}}},
      {"PhiRecPhiGen", "; #phi_{Rec}; #phi_{Gen}; tracks", {HistType::kTH2F, {{600, 0, 2*M_PI}, {600, 0, 2*M_PI}}}},

    }                                                                                //
  };



  void processGen(aod::McCollisions::iterator const& mcCollision, o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels>> const& collisions, Particles const& particles, soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& tracks)
  {
    TLorentzVector vTrack;
    std::vector<PseudoJet> particlesRec;

    TLorentzVector vTrackPart;
    std::vector<PseudoJet> particlesTrackRec;

    TLorentzVector vPart;
    std::vector<PseudoJet> particlesGen;

    float Ntracks = 0;
    float Nparts = 0;

    std::vector<fastjet::PseudoJet> jetsRec;
    std::vector<fastjet::PseudoJet> jetsGen;

    //printf("%f\n", tracks.rawIteratorAt(0).eta());

    //LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
    printf("MC col %lld has %lld reco cols\n", mcCollision.globalIndex(), collisions.size());

    for (auto& collision : collisions)
    {
      //printf("collision index : %lld\n", collision.globalIndex());
      auto groupedTracks = tracks.sliceBy(o2::aod::fwdtrack::collisionId, collision.globalIndex());
      Ntracks = groupedTracks.size();
      for (auto& track : groupedTracks)
      {
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        vTrack.SetPtEtaPhiM(0.2,track.eta(),phi,0.138);
        //---test
        //

        if (track.has_mcParticle())
        {
          auto particleTrack = track.mcParticle();
          //printf("SIZE===%llu\n", particles.size());

          //printf("----pt %f\n", particleTrack.pt());
          auto p = pdg->GetParticle(particleTrack.pdgCode());
          int mass = 0;
          if (p != nullptr)
          {
            mass = p->Mass();
          }
          vTrack.SetPtEtaPhiM(particleTrack.pt(),track.eta(),phi,mass);
          registry.fill(HIST("EtaRecEtaGen"), track.eta(), particleTrack.eta());
          registry.fill(HIST("PhiRecPhiGen"), phi, particleTrack.phi());
        }
        //printf("MCPARTICLEID  %d\n", testnb);

        //---end of test
        particlesRec.push_back(PseudoJet(vTrack.Px(), vTrack.Py(), vTrack.Pz(), vTrack.E()));
        particlesRec[particlesRec.size()-1].set_user_index(track.mcParticleId());//set_user_index

      }
    }

    //Here it tries to do a jet with the reconstructed tracks that can come from diff collisions, but same
    //mcCollision (for now), most of the time there is just one or 0 collision
    //This does not really represent the reality -> we should put the clustering in the collision loop
    fastjet::ClusterSequenceArea cs(particlesRec, *jet_defRec, *areaDef);
    jetsRec = cs.inclusive_jets(0.0);

    for (unsigned i = 0; i < jetsRec.size(); i++)
    {
      //printf("jet %d : pt %f y %f phi %f\n", i, jetsRec[i].pt(), jetsRec[i].rap(), jetsRec[i].phi());
      std::vector<PseudoJet> constituents = jetsRec[i].constituents();
      if (constituents.size()==1)//we ignore the jets with just one constituent
      {
        continue;
      }

      registry.fill(HIST("AreaRecojet"), jetsRec[i].area());
      registry.fill(HIST("NTracksPerJet"), constituents.size());
      for (unsigned j = 0; j < constituents.size(); j++)
      {
      //  printf("constituent %d : pt %f, y %f, phi %f, index %d\n", j, constituents[j].pt(), constituents[j].rap(), constituents[j].phi(), constituents[j].user_index());
      }
    }


    for (auto& particle : particles)
    {
      auto p = pdg->GetParticle(particle.pdgCode());
      int charge = 0;
      if (p != nullptr)
      {
        charge = p->Charge();
      }
      if (charge != 0 && particle.isPhysicalPrimary() && (particle.eta() < -2.3 && particle.eta() > -3.6) && particle.p()>0.05)
      //charged primary and within the MFT acceptance
      {
        Nparts++;
        vPart.SetPtEtaPhiM(particle.pt(),particle.eta(),particle.phi(),p->Mass());
        particlesGen.push_back(PseudoJet(vPart.Px(), vPart.Py(), vPart.Pz(), vPart.E()));
        particlesGen[particlesGen.size()-1].set_user_index(particle.globalIndex());
      }
    }


    fastjet::ClusterSequenceArea csGen(particlesGen, *jet_defGen, *areaDef);

    jetsGen = csGen.inclusive_jets(0.0);
    int nTrueConst;
    int nJetGen=0;
    int nTrueJet=0;
    for (unsigned i = 0; i < jetsGen.size(); i++)
    {
      //printf("jet %d : pt %f y %f phi %f\n", i, jetsGen[i].pt(), jetsGen[i].rap(), jetsGen[i].phi());
      std::vector<PseudoJet> constituents = jetsGen[i].constituents();
      if (constituents.size()==1)//we ignore the jets with just one constituent
      {
        continue;
      }
      nJetGen++;
      registry.fill(HIST("AreaMCjetVsPt"), jetsGen[i].area(), jetsGen[i].pt());
      registry.fill(HIST("AreaMCjet"), jetsGen[i].area());
      registry.fill(HIST("NPartsPerJet"), constituents.size());
      for (unsigned l = 0; l < jetsRec.size(); l++)
      {
        //printf("jet %d : pt %f y %f phi %f\n", i, jetsRec[i].pt(), jetsRec[i].rap(), jetsRec[i].phi());
        nTrueConst = 0;
        std::vector<PseudoJet> constituentsRec = jetsRec[l].constituents();
        for (unsigned m = 0; m < constituentsRec.size(); m++)
        {
          for (unsigned j = 0; j < constituents.size(); j++)
          {
            //printf("constituent %d : pt %f, y %f, phi %f, index %d\n", j, constituents[j].pt(), constituents[j].rap(), constituents[j].phi(), constituents[j].user_index());
            //particles[mother0Id]

            if (constituentsRec[m].user_index()==constituents[j].user_index())
            {
              nTrueConst++;
              //printf("yeah %d\n", nTrueConst);
            }
          }
        }
        if (nTrueConst/float(constituents.size())>0.7)
        {
          //printf("one true jet with purity %f, and number of constituents %lu\n", nTrueConst/float(constituents.size()), constituents.size());
          nTrueJet++;
        }
      }

    }

    int nTrueJet2=0;
    int nJetRec=0;
    for (unsigned i = 0; i < jetsRec.size(); i++)
    {
      //printf("jet %d : pt %f y %f phi %f\n", i, jetsRec[i].pt(), jetsRec[i].rap(), jetsRec[i].phi());
      std::vector<PseudoJet> constituentsRec = jetsRec[i].constituents();

      if (constituentsRec.size()==1)//we ignore the jets with just one constituent
      {
        continue;
      }
      nJetRec++;

      for (unsigned j = 0; j < jetsGen.size(); j++)//we compare jetRec i with jetGen j
      {
        nTrueConst=0;//for now no constituent in common
        std::vector<PseudoJet> constituentsGen = jetsGen[j].constituents();
        if (constituentsGen.size()==1)//we ignore the jets with just one constituent
        {
          continue;
        }

        for (unsigned l = 0; l < constituentsGen.size(); l++)
        {
          for (unsigned m = 0; m < constituentsRec.size(); m++)
          {
            if (constituentsRec[m].user_index()==constituentsGen[l].user_index())
            {
              nTrueConst++;
            }
          }

        }

        if (nTrueConst/float(constituentsRec.size())>0.7)
        {
          //printf("one true jet with purity %f, and number of constituents %lu\n", nTrueConst/float(constituents.size()), constituents.size());
          nTrueJet2++;
        }
      }


    }
    if (nJetRec>0)
    {
      printf("-------In this mcCollision there was %d jet generated and %d jets reco %d true jets reconstructed and true jets 2 reco %d\n", nJetGen, nJetRec, nTrueJet, nTrueJet2);
    }

    registry.fill(HIST("NTracksOverNparts"), Ntracks/Nparts);
    registry.fill(HIST("NJetGenOverNJetTrue"), nJetGen/1.0*nTrueJet);
    registry.fill(HIST("NJetGenOverNJetTrue2"), 1.0*nTrueJet2/nJetGen);
  }

  PROCESS_SWITCH(analyseMFTJets, processGen, "Process gen level", true);
};
//end of the task analyseMFTJets

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<analyseMFTJets>(cfgc),
  };
}
