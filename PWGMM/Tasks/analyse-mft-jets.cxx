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

using namespace o2;
using namespace o2::framework;
using namespace fastjet;

using Particles = aod::McParticles;

//the task analyseMFTTracks loops over MFT tracks and generated particles and fills basic histograms

struct analyseMFTJets {
  int icoll = 0;
  Service<TDatabasePDG> pdg;
  double R = 0.7;
  JetDefinition jet_def(cambridge_algorithm, R);
  HistogramRegistry registry{
    "registry",
    {
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, -M_PI, M_PI}, {35, -4.5, -1.}}}}, //

    }                                                                                //
  };

  void processRec(o2::aod::Collision const& collision, o2::aod::MFTTracks const& tracks)
  {

    auto z = collision.posZ();

    //auto groupedTracks = tracks.sliceBy(o2::aod::fwdtrack::collisionId, collision.globalIndex());


    for (auto& track : tracks)
    {
      registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      registry.fill(HIST("TracksPhiEta"), phi, track.eta());
    }
  }
  //end of processRec
  PROCESS_SWITCH(analyseMFTTracks, processRec, "Process rec level", true);

  void processGen(aod::McCollisions::iterator const& mcCollision, soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions, Particles const& particles, soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& tracks)
  {
    TLorentzVector vTrack;
    std::vector<PseudoJet> particlesRec;

    int nChargedPrimaryParticles = 0;
    auto z = mcCollision.posZ();

    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());
//printf les constituants, le mctracklabel
//boucle sur les collisions associées à une mccollision, tracks associés à cette collision et remplir un psuedojet xsavec
    for (auto& collision : collisions)
    {
      auto groupedTracks = tracks.sliceBy(o2::aod::fwdtrack::collisionId, collision.globalIndex());
      for (auto& track : groupedTracks)
      {
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        vTrack.setPtEtaPhiM(0.2,track.eta(),phi,0.138);
        particlesRec.push_back(PseudoJet(vTrack.Px(), vTrack.Py(), vTrack.Pz(), vTrack.E()));
        //TLorentz vector setPtEtaPhiM avec m masse pion et pt=0.2
        //pseudoJet avec v.Px(),v.Py(), v.Pz(), v.E()
        //.pushback cf examplefastjet
        particlesRec[particlesRec.size()-1].user_index(track.mcParticleId());
        //user_index (track.McParticleId()) ou un truc du genre
      }
      //cluster bidule
      //pour chaque jet ayant plus d'une particle imprimer les trucs
    }




    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());

      int charge = 0;
      if (p != nullptr)
      {
        charge = p->Charge();
      }
      if (charge != 0 && particle.isPhysicalPrimary()) {
        registry.fill(HIST("TracksEtaZvtxGen"), particle.eta(), z);
        registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());

        registry.fill(HIST("TracksPhiEtaGen"), particle.phi(), particle.eta());
        registry.fill(HIST("NtrkEtaGen"), particle.eta());
        nChargedPrimaryParticles++;
      }
    }

  }

  PROCESS_SWITCH(analyseMFTTracks, processGen, "Process gen level", false);
};
//end of the task analyseMFTTracks

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<analyseMFTTracks>(cfgc),
  };
}
