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

/// \file fwdtrackToCollisionAssociator.cxx
/// \brief Associates fwd and MFT tracks to collisions considering ambiguities
/// \author Sarah Herrmann <sarah.herrmann@cern.ch>, IP2I Lyon
/// \author Maurice Coquet

#include "Common/DataModel/CollisionAssociation.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/LHCConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct FwdTrackToCollisionAssociation {
  Produces<MFTTrackAssoc> association;
  Produces<MFTTrackCompColls> reverseIndices;

  Configurable<float> nSigmaForTimeCompat{"nSigmaForTimeCompat", 4.f, "number of sigmas for time compatibility"};
  Configurable<float> timeMargin{"timeMargin", 0.f, "time margin in ns added to uncertainty because of uncalibrated TPC"};
  Configurable<bool> applyMinimalTrackSelForRun2{"applyMinimalTrackSelForRun2", false, "flag to apply minimal track selection for Run 2 in case of standard association"};
  Configurable<bool> applyIsGlobalTrackWoDCA{"applyIsGlobalTrackWoDCA", true, "flag to apply global track w/o DCA selection"};
  Configurable<bool> usePVAssociation{"usePVAssociation", true, "if the track is a PV contributor, use the collision time for it"};
  Configurable<bool> includeUnassigned{"includeUnassigned", false, "consider also tracks which are not assigned to any collision"};
  Configurable<bool> fillTableOfCollIdsPerTrack{"fillTableOfCollIdsPerTrack", false, "fill additional table with vector of collision ids per track"};

  // Filter trackFilter = (applyIsGlobalTrackWoDCA.node() == false) || requireGlobalTrackWoDCAInFilter();
  // using TracksWithSel = soa::Join<Tracks, TracksExtra, TrackSelection>;
  // using TracksWithSelFilter = soa::Filtered<TracksWithSel>;

  void init(InitContext const&)
  {
    // if (doprocessMFTAssocWithTime == doprocessMFTStandardAssoc) {
    //   LOGP(fatal, "Exactly one process function should be enabled!");
    // }

    // if (applyIsGlobalTrackWoDCA && applyMinimalTrackSelForRun2) {
    //   LOGP(fatal, "You cannot apply simultaneously Run 2 track selections and applyIsGlobalTrackWoDCA!");
    // }
  }

  void processMFTAssocWithTime(Collisions const& collisions,
                            MFTTracks const& tracks,
                            AmbiguousMFTTracks const& ambiguousTracks,
                            BCs const& bcs)
  {
    // cache globalBC
    std::vector<uint64_t> globalBC;
    for (const auto& track : tracks) {
      if (track.has_collision()) {
        globalBC.push_back(track.collision().bc().globalBC());
      } else {
        for (const auto& ambTrack : ambiguousTracks) {
          if (ambTrack.mfttrackId() == track.globalIndex()) {
            globalBC.push_back(ambTrack.bc().begin().globalBC());
            break;
          }
        }
      }
    }

    // define vector of vectors to store indices of compatible collisions per track
    std::vector<std::unique_ptr<std::vector<int>>> collsPerTrack(tracks.size());

    // loop over collisions to find time-compatible tracks
    auto trackBegin = tracks.begin();
    constexpr auto bOffsetMax = 241; // 6 mus (ITS)
    for (const auto& collision : collisions) {
      const float collTime = collision.collisionTime();
      const float collTimeRes2 = collision.collisionTimeRes() * collision.collisionTimeRes();
      uint64_t collBC = collision.bc().globalBC();
      // bool iteratorMoved = false;
      for (auto track = trackBegin; track != tracks.end(); ++track) {
        if (!includeUnassigned && !track.has_collision()) {
          continue;
        }
        const int64_t bcOffset = (int64_t)globalBC[track.filteredIndex()] - (int64_t)collBC;
        if (std::abs(bcOffset) > bOffsetMax) {
          continue;
        }

        float trackTime = track.trackTime();
        float trackTimeRes = track.trackTimeRes();


        const float deltaTime = trackTime - collTime + bcOffset * constants::lhc::LHCBunchSpacingNS;
        float sigmaTimeRes2 = collTimeRes2 + trackTimeRes * trackTimeRes;
        LOGP(debug, "collision time={}, collision time res={}, track time={}, track time res={}, bc collision={}, bc track={}, delta time={}", collTime, collision.collisionTimeRes(), track.trackTime(), track.trackTimeRes(), collBC, globalBC[track.filteredIndex()], deltaTime);


        float thresholdTime = 0.;

        thresholdTime = nSigmaForTimeCompat * std::sqrt(sigmaTimeRes2) + timeMargin;


        if (std::abs(deltaTime) < thresholdTime) {
          const auto collIdx = collision.globalIndex();
          const auto trackIdx = track.globalIndex();
          LOGP(debug, "Filling track id {} for coll id {}", trackIdx, collIdx);
          association(collIdx, trackIdx);
          if (fillTableOfCollIdsPerTrack) {
            if (collsPerTrack[trackIdx] == nullptr) {
              collsPerTrack[trackIdx] = std::make_unique<std::vector<int>>();
            }
            collsPerTrack[trackIdx].get()->push_back(collIdx);
          }
        }
      }
    }

    // create reverse index track to collisions if enabled
    if (fillTableOfCollIdsPerTrack) {
      std::vector<int> empty{};
      for (const auto& track : tracks) {

        const auto trackId = track.globalIndex();
        if (collsPerTrack[trackId] == nullptr) {
          reverseIndices(empty);
        } else {
          //comment again
          printf("track ID %lld\n", trackId);
          printf("nb of coll %lu\n", collsPerTrack[trackId]->size());
          for (const auto& collId : *collsPerTrack[trackId]) {
            printf("  -> Coll id %d\n", collId);
          }
          //end of comment again
          reverseIndices(*collsPerTrack[trackId].get());
        }
      }
    }
  }

  PROCESS_SWITCH(FwdTrackToCollisionAssociation, processMFTAssocWithTime, "Use MFTtrack-to-collision association based on time", true);
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FwdTrackToCollisionAssociation>(cfgc)};
}
