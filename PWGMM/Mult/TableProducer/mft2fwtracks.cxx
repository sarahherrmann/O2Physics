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

///
/// \file   mft2fwdtracks.cxx
/// \author Sarah Herrmann
/// \since  28-06-2023
/// \brief  A task to create a reverse index from MFTTracks to FwdTracks
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"

#include "FwdIndex.h"

using namespace o2;
using namespace o2::framework;

struct MFTtoFwd {
  //using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
  Produces<o2::aod::MFTTracksToFwdTracks> m2f;


  //std::vector<int> trackIds;

  void init(InitContext&)
  {
  }

  // void process(aod::McParticles const& particles)
  // {
  //   if (doprocessIndexingCentral) {
  //     p2t.reserve(particles.size());
  //   }
  //   if (doprocessIndexingFwd) {
  //     p2tmft.reserve(particles.size());
  //   }
  // }

  void process(aod::MFTTrack const& mfttrack, soa::SmallGroups<aod::FwdTracks> const& fwdtracks)
  {
    if (fwdtracks.size() > 1)
    {
      //we choose to ignore the fwdtracks that match multiple times to one MFT track
      m2f(-1);
      //printf("THIS IS STRANGE HERE %d\n", 1);
      return;
    }
    int fwdtrackId=-1;
    // if (fwdtracks.size()>1)
    // {
    //   printf("fwdtracks.size() %d\n", fwdtracks.size());
    // }

    for (auto& fwdtrack : fwdtracks)
    {
      //auto mfttrackMatchIdx=fwdtrack.matchMFTTrackId();//NO ! I want to fill the table with fwdtrack.globalIndex()
      fwdtrackId = fwdtrack.globalIndex();
      // if (fwdtracks.size()>1)
      // {
      //   printf("mfttrackMatchIdx %d, MFT index %lld\n", mfttrackMatchIdx, mfttrack.globalIndex());
      // }
    }
    m2f(fwdtrackId);
  }

  //PROCESS_SWITCH(MFTtoFwd, processIndexingFwd, "Create reverse index from particles to tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MFTtoFwd>(cfgc)};
}
