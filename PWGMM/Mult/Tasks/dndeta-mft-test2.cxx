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

#include <cmath>
// for CCDB access
#include <chrono>
#include "CCDB/BasicCCDBManager.h"

#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "TDatabasePDG.h"
#include "MathUtils/Utils.h"

#include "bestCollisionTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;


int nVis=0;



struct PseudorapidityDensityMFT {

  using CollwEv = soa::Join<aod::Collisions, aod::EvSels>;

  void processTest(CollwEv::iterator const& collision,
                   o2::aod::MFTTracks const&,
                   soa::SmallGroups<aod::BestCollisionsFwd> const& atracks)
  {
    int i = 0;
    for (auto& track : atracks){
      auto normalTrack = track.mfttrack();

      if (normalTrack.collisionId() != track.bestCollisionId())
      {
        printf("----------------- nVis %d, trackIndex %d, collisionOnGoing %d, bestCol %d, normalTrack.collisionId() %d, \n", nVis, track.mfttrackId(), collision.globalIndex(), track.bestCollisionId(), normalTrack.collisionId());
        i++;
      }

    }
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processTest, "Process test", true);

};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensityMFT>(cfgc)};
}
