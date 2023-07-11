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

// \file   nch-mft.cxx
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief Fills the distribution of Nch inside the MFT (no correction)
//
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/EventSelection.h"

#include "FwdIndex.h"
#include "bestCollisionTable.h"

#include "Common/DataModel/CollisionAssociationTables.h"

using namespace o2;
using namespace o2::framework;

using CollwEv = soa::Join<aod::Collisions, aod::EvSels>;

struct nchMft {

  Service<O2DatabasePDG> pdg;

  HistogramRegistry registry{
    "registry",
    {
      {"Nch", "; N_{trk}^{MFT}; count", {HistType::kTH1F, {{701, -0.5, 700.5}}}},
      {"Nch0Z", "; #it{z}_{vtx} (cm); count", {HistType::kTH1F, {{301, -30.1, 30.1}}}},
      {"NchWEvSel", "; N_{trk}^{MFT}; count", {HistType::kTH1F, {{701, -0.5, 700.5}}}},
      {"SameCollFwd", "; sameColl; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"SameCollFwdRe", "; sameCollRe; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"IsInMFTColls", "; isInMFTColls; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"IsInMFTCollsRe", "; isInMFTCollsRe; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}},
      {"NchGen", "; N_{trk}^{MFT}; count", {HistType::kTH1F, {{701, -0.5, 700.5}}}},
      {"NchWEvSelGen", "; N_{trk}^{MFT}; count", {HistType::kTH1F, {{701, -0.5, 700.5}}}},
      {"Ncolls", "; N_coll; count", {HistType::kTH1F, {{21, -0.5, 20.5}}}}
    }
  };

  void init(o2::framework::InitContext& initContext)
  {
    if (doprocessFwdAmb) {
      registry.add({"FwdTracksAmbDegree", " ; N_{coll}^{comp}", {HistType::kTH1I, {{41, -0.5, 40.5}}}});
      registry.add({"FwdTrackIsAmb", " ; trackType_isAmbiguous", {HistType::kTH1I, {{5, -0.5, 4.5}}}});
      registry.add({"FwdTrackNonAmb", " ; trackType_isNotAmbiguous", {HistType::kTH1I, {{5, -0.5, 4.5}}}});
    }
  }

  //expressions::Filter MFTtrackFilter = (aod::fwdtrack::eta < -2.5f) && (aod::fwdtrack::eta > -3.6f);

  using FiMFTTracks = soa::Filtered<aod::MFTTracks>;
  expressions::Filter centralFilter = (nabs(aod::track::eta) < 1.1f);

  using FiCentralTracks = soa::Filtered<aod::Tracks>; // central tracks for INEL>0

  void processMult(CollwEv::iterator const& collision, aod::MFTTracks const& tracks, FiCentralTracks const& midtracks)
  {
    if (midtracks.size()==0)
    {
      return;
    }
    registry.fill(HIST("Nch"), tracks.size());
    if (collision.sel8())
    {
      registry.fill(HIST("NchWEvSel"), tracks.size());
    }

    if (tracks.size()==0)
    {
      registry.fill(HIST("Nch0Z"), collision.posZ());
    }
  }
  PROCESS_SWITCH(nchMft, processMult, "Collect mult distribution", true);



  using Particles = soa::Filtered<aod::McParticles>;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary && (aod::mcparticle::eta < -2.5f) && (aod::mcparticle::eta > -3.6f);;



  void processMultGen(aod::McCollisions::iterator const& mcCollision,
                  o2::soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>> const& collisions,
                  Particles const& fwdParticles)
  {
    int nCharged = 0;
    for (auto& particle : fwdParticles)
    {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      nCharged++;
    }

    registry.fill(HIST("NchGen"), nCharged);
    registry.fill(HIST("Ncolls"),collisions.size());
    for (auto& collision : collisions) {
      if (collision.sel8())
      {
        registry.fill(HIST("NchWEvSelGen"), nCharged);
      }
    }

  }
  PROCESS_SWITCH(nchMft, processMultGen, "Collect mult distribution at generator level", false);


  using MFTwFwd = soa::Join<aod::MFTTracks, o2::aod::MFTTracksToFwdTracks, aod::MFTTrkCompColls>;
  using FwdTracksWColls = soa::Join<o2::aod::FwdTracks, aod::FwdTrkCompColls>;

  void processMuonsAndMFT(CollwEv::iterator const& collision, MFTwFwd const& mfttracks, FwdTracksWColls const&)
  {
    for (auto& mfttrack : mfttracks)
    {

      if (!mfttrack.has_fwdtrack())
      {
        continue;
      }
      auto matchedFwdtrack = mfttrack.fwdtrack_as<FwdTracksWColls>();
      if ((matchedFwdtrack.trackType() != 0)||(matchedFwdtrack.compatibleCollIds().size()>1))
      {
        continue;
      }//we just want to look at the global muon tracks, to have the least amount of ambiguous fwd tracks


      int collFwd = matchedFwdtrack.collisionId();// the collision of the global muon track
      int sameColl = 0;
      int isInMFTColls = 0;

      auto compatibleColls = mfttrack.compatibleColl_as<CollwEv>();
      for (auto& collMFT : compatibleColls)
      {//looking at all the collisions compatible with the MFT track
        if (collMFT.globalIndex() == collFwd)
        {//is collFwd part of these compatible collisions ?
          isInMFTColls=1;
        }
      }

      if (collFwd == collision.globalIndex())
      {
        sameColl = 1;
      }

      registry.fill(HIST("SameCollFwd"), sameColl);
      registry.fill(HIST("IsInMFTColls"), isInMFTColls);

    }
  }
  PROCESS_SWITCH(nchMft, processMuonsAndMFT, "Chi2 of muons and MFt track matching 2", true);



  void processMuonsAndMFTReas(CollwEv::iterator const& collision, soa::SmallGroups<aod::BestCollisionsFwd> const& retracks, MFTwFwd const&, FwdTracksWColls const&)
  {
    for (auto& retrack : retracks) {
      auto mfttrack = retrack.mfttrack_as<MFTwFwd>();
      if (!mfttrack.has_fwdtrack())
      {
        continue;
      }
      auto matchedFwdtrack = mfttrack.fwdtrack_as<FwdTracksWColls>();
      if ((matchedFwdtrack.trackType() != 0)||(matchedFwdtrack.compatibleCollIds().size()>1))
      {
        continue;
      }//we just want to look at the global muon tracks, to have the least amount of ambiguous fwd tracks
      int collFwd = matchedFwdtrack.collisionId();
      int sameColl = 0;
      int isInMFTColls = 0;

      auto compatibleColls = mfttrack.compatibleColl_as<CollwEv>();
      for (auto& collMFT : compatibleColls)
      {//looking at all the collisions compatible with the MFT track
        if (collMFT.globalIndex() == collFwd)
        {//is collFwd part of these compatible collisions ?
          isInMFTColls=1;
        }
      }

      if (collFwd == collision.globalIndex())
      {
        sameColl = 1;
      }

      registry.fill(HIST("SameCollFwdRe"), sameColl);
      registry.fill(HIST("IsInMFTCollsRe"), isInMFTColls);

    }
  }
  PROCESS_SWITCH(nchMft, processMuonsAndMFTReas, "Chi2 of muons and MFt track matching", true);



  void processFwdAmb(FwdTracksWColls::iterator const& track)
  {
    auto compatibleColls = track.compatibleCollIds();
    int isAmbiguous = 0;
    registry.fill(HIST("FwdTracksAmbDegree"), compatibleColls.size());
    if (compatibleColls.size() > 1) {
      isAmbiguous = 1;
    }
    registry.fill(HIST("FwdTrackIsAmb"), track.trackType(), isAmbiguous);
    registry.fill(HIST("FwdTrackNonAmb"), track.trackType(), !isAmbiguous);
  }
  PROCESS_SWITCH(nchMft, processFwdAmb, "find the percentage of fwd ambiguous tracks", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nchMft>(cfgc)};
}
