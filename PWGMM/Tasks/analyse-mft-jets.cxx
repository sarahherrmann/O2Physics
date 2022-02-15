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
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/AreaDefinition.hh>

using namespace o2;
using namespace o2::framework;
using namespace fastjet;

using Particles = aod::McParticles;

//the task analyseMFTJets loops over MFT tracks and generated particles and fills basic histograms

struct analyseMFTJets {
  int icoll = 0;
  Service<TDatabasePDG> pdg;
  double R = 0.7;
  JetDefinition jet_def(cambridge_algorithm, R);

  fastjet::GhostedAreaSpec ghostareaspec(5, 1, 0.01); //ghost of maxrap=5, 1 repetition, small ghost area
  fastjet::AreaType areaType = fastjet::active_area;
  fastjet::AreaDefinition *areaDef = new fastjet::AreaDefinition(areaType, ghostareaspec);

  HistogramRegistry registry{
    "registry",
    {
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //
      {"NTracksOverNparts", "Number of reconstructed tracks / number of particles in a collision; Ntracks/NParts; #count", {HistType::kTH1F, {{300, 0, 150}}}}, //
      {"NTracksPerJet", "Number of reconstructed tracks per jet; Ntracks/jet; #count", {HistType::kTH1F, {{30, 1, 30}}}}, //
      {"NPartsPerJet", "Number of particles per jet; Nparts/jet; #count", {HistType::kTH1F, {{30, 1, 30}}}}, //
      {"AreaMCjet", "area of jet in gen MC; area; #count", {HistType::kTH1F, {{300, 0, 30}}}}, //
      {"AreaRecojet", "area of jet in reco MC; area; #count", {HistType::kTH1F, {{300, 0, 30}}}}, //


    }                                                                                //
  };

  void processRec(o2::aod::Collision const& collision, o2::aod::MFTTracks const& tracks)
  {

    auto z = collision.posZ();

    //auto groupedTracks = tracks.sliceBy(o2::aod::fwdtrack::collisionId, collision.globalIndex());


    for (auto& track : tracks)
    {
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      registry.fill(HIST("TracksPhiEta"), phi, track.eta());
    }
  }
  //end of processRec
  PROCESS_SWITCH(analyseMFTJets, processRec, "Process rec level", true);

  void processGen(aod::McCollisions::iterator const& mcCollision, soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions, Particles const& particles, soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& tracks)
  {
    TLorentzVector vTrack;
    std::vector<PseudoJet> particlesRec;

    TLorentzVector vPart;
    std::vector<PseudoJet> particlesGen;

    float Ntracks = 0;
    float Nparts = 0;

    //int nChargedPrimaryParticles = 0;
    //auto z = mcCollision.posZ();

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());

    for (auto& collision : collisions)
    {
      printf("collision index : %d\n", collision.globalIndex());
      auto groupedTracks = tracks.sliceBy(o2::aod::fwdtrack::collisionId, collision.globalIndex());
      Ntracks = groupedTracks.size();
      for (auto& track : groupedTracks)
      {
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        vTrack.setPtEtaPhiM(0.2,track.eta(),phi,0.138);
        particlesRec.push_back(PseudoJet(vTrack.Px(), vTrack.Py(), vTrack.Pz(), vTrack.E()));
        particlesRec[particlesRec.size()-1].set_user_index(track.mcParticleId());//set_user_index

      }

      std::vector<fastjet::PseudoJet> jetsRec;
      fastjet::ClusterSequenceArea cs(particlesRec, jet_def, areaDef);

      jetsRec = cs.inclusive_jets(0.0);

      for (unsigned i = 0; i < jetsRec.size(); i++)
      {
        printf("jet %d : pt %f y %f phi %f\n", i, jetsRec[i].pt(), jetsRec[i].rap(), jetsRec[i].phi());
        vector<PseudoJet> constituents = jetsRec[i].constituents();

        registry.fill(HIST("AreaRecojet"), jetsRec[i].area());
        registry.fill(HIST("NTracksPerJet"), constituents.size());
        for (unsigned j = 0; j < constituents.size(); j++)
        {
          printf("constituent %d : pt %f, y %f, phi %f, index %d\n", j, constituents[j].pt(), constituents[j].rap(), constituents[j].phi(), constituents[j].user_index());
        }
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
      if (charge != 0 && particle.isPhysicalPrimary() && (particle.eta() < -2.3 && particle.eta() > -3.6))
      //charged primary and within the MFT acceptance
      {
        Nparts++;
        vPart.setPtEtaPhiM(particle.pt(),particle.eta(),particle.phi(),p->Mass());
        particlesGen.push_back(PseudoJet(vPart.Px(), vPart.Py(), vPart.Pz(), vPart.E()));
        particlesGen[particlesGen.size()-1].set_user_index(particle.globalIndex());
      }
    }

    std::vector<fastjet::PseudoJet> jetsGen;
    fastjet::ClusterSequenceArea csGen(particlesGen, jet_def, areaDef);

    jetsGen = csGen.inclusive_jets(0.0);

    for (unsigned i = 0; i < jetsGen.size(); i++)
    {
      printf("jet %d : pt %f y %f phi %f\n", i, jetsGen[i].pt(), jetsGen[i].rap(), jetsGen[i].phi());
      vector<PseudoJet> constituents = jetsGen[i].constituents();

      registry.fill(HIST("AreaMCjet"), jetsGen[i].area());
      registry.fill(HIST("NPartsPerJet"), constituents.size());
      for (unsigned j = 0; j < constituents.size(); j++)
      {
        printf("constituent %d : pt %f, y %f, phi %f, index %d\n", j, constituents[j].pt(), constituents[j].rap(), constituents[j].phi(), constituents[j].user_index());
        //particles[mother0Id]
      }
    }


    registry.fill(HIST("NTracksOverNparts"), Ntracks/Nparts);
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
