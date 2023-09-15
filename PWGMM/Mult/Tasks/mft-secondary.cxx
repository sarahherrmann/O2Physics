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

// \file   mft-secondary.cxx
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief This code loops over every ambiguous MFT tracks and plots
//        various DCA info
//
// \dependencies mm-ambiguous-track-propagation, fwdtrack-to-collision-associator
//
// \date 01-08-2023


#include "MathUtils/Utils.h"
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/LHCConstants.h"
#include "Common/DataModel/EventSelection.h"

#include "bestCollisionTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;


using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
using CollisionsLabeled = soa::Join<o2::aod::Collisions, aod::McCollisionLabels>;
//using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;


AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec PhiAxis = {600, 0, 2 * M_PI};
AxisSpec EtaAxis = {18, -4.6, -1.};
AxisSpec DCAxyAxis = {1000, -0.1, 30};
AxisSpec DCAxAxis = {1000, -10, 10};
AxisSpec DCAyAxis = {1000, -10, 10};

struct mftsecondary {

  Configurable<float> maxDCAXY{"maxDCAXY", 6.0, "max allowed transverse DCA"}; // To be used when associating ambitrack to collision using best DCA


  HistogramRegistry registry{
    "registry",
    {
      {"Tracks_NoPart_EtaZvtx", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},

      {"TracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},
      {"DCAxy", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {DCAxyAxis}}},
      {"DCAx", "; DCA_{x} (cm); counts", {HistType::kTH1F, {DCAxAxis}}},
      {"DCAy", "; DCA_{y} (cm); counts", {HistType::kTH1F, {DCAyAxis}}},

      {"TracksEtaZvtx_leftBin", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},
      {"DCAxy_leftBin", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {DCAxyAxis}}},
      {"DCAx_leftBin", "; DCA_{x} (cm); counts", {HistType::kTH1F, {DCAxAxis}}},
      {"DCAy_leftBin", "; DCA_{y} (cm); counts", {HistType::kTH1F, {DCAyAxis}}},

      {"TracksEtaZvtx_rightBin", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},
      {"DCAxy_rightBin", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {DCAxyAxis}}},
      {"DCAx_rightBin", "; DCA_{x} (cm); counts", {HistType::kTH1F, {DCAxAxis}}},
      {"DCAy_rightBin", "; DCA_{y} (cm); counts", {HistType::kTH1F, {DCAyAxis}}},

      {"Primary/TracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},
      {"Primary/DCAxy", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {DCAxyAxis}}},
      {"Primary/DCAx", "; DCA_{x} (cm); counts", {HistType::kTH1F, {DCAxAxis}}},
      {"Primary/DCAy", "; DCA_{y} (cm); counts", {HistType::kTH1F, {DCAyAxis}}},

      {"Secondary/TracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},
      {"Secondary/DCAxy", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {DCAxyAxis}}},
      {"Secondary/DCAx", "; DCA_{x} (cm); counts", {HistType::kTH1F, {DCAxAxis}}},
      {"Secondary/DCAy", "; DCA_{y} (cm); counts", {HistType::kTH1F, {DCAyAxis}}},
    }
  };

  void init(InitContext&)
  {
    // auto hstatus = registry.get<TH1>(HIST("AmbiguousTrackStatus"));
    // auto* x2 = hstatus->GetXaxis();
    // x2->SetBinLabel(1, "MFT tracks");
    // x2->SetBinLabel(2, "MFT ambiguous tracks");
    // x2->SetBinLabel(3, "Reassigned tracks");
    // x2->SetBinLabel(4, "Extra tracks");
    // x2->SetBinLabel(5, "orig=true (re)");
    // x2->SetBinLabel(6, "best=true (re)");
    // x2->SetBinLabel(7, "not reassigned");
    // x2->SetBinLabel(8, "not reassigned and true");
  }

  using CollwEv = soa::Join<aod::Collisions, aod::EvSels>;

  void processDCAData(Collision const& collision,
                  aod::MFTTracks const&,
                  soa::SmallGroups<aod::BestCollisionsFwd> const& retracks)
  {
    if (retracks.size() == 0) {
      return;
    }

    for (auto& retrack : retracks)
    {
      auto track = retrack.mfttrack();
      if (retrack.ambDegree()>1)
      {
        //this track is ambiguous, ignore it for now
        continue;
      }

      float z = collision.posZ();


      if (track.eta()>=-2.8 && track.eta()<=-2.4)
      {
        registry.fill(HIST("TracksEtaZvtx_rightBin"), track.eta(), z);
        registry.fill(HIST("DCAx_rightBin"), retrack.bestDCAX());
        registry.fill(HIST("DCAy_rightBin"), retrack.bestDCAY());
        registry.fill(HIST("DCAxy_rightBin"), retrack.bestDCAXY());
      }

      if (track.eta()>=-3.4 && track.eta()<=-3.0)
      {
        registry.fill(HIST("TracksEtaZvtx_leftBin"), track.eta(), z);
        registry.fill(HIST("DCAx_leftBin"), retrack.bestDCAX());
        registry.fill(HIST("DCAy_leftBin"), retrack.bestDCAY());
        registry.fill(HIST("DCAxy_leftBin"), retrack.bestDCAXY());
      }

      registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
      registry.fill(HIST("DCAx"), retrack.bestDCAX());
      registry.fill(HIST("DCAy"), retrack.bestDCAY());
      registry.fill(HIST("DCAxy"), retrack.bestDCAXY());


    }
  }
  PROCESS_SWITCH(mftsecondary, processDCAData, "get the DCAxy of MFT tracks", true);

  void processDCA(CollwEv::iterator const& collision,
                  MFTTracksLabeled const&,
                  soa::SmallGroups<aod::BestCollisionsFwd> const& retracks,
                  aod::McParticles const&,
                  aod::McCollisions const&)
  {
    if (retracks.size() == 0) {
      return;
    }

    for (auto& retrack : retracks)
    {
      auto track = retrack.mfttrack_as<MFTTracksLabeled>();
      if (retrack.ambDegree()>1)
      {
        //this track is ambiguous, ignore it for now
        continue;
      }
      // if (!track.has_collision()) {
      //   continue;
      // }
      // auto collision = track.collision();
      float z = collision.posZ();

      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        registry.fill(HIST("Tracks_NoPart_EtaZvtx"), track.eta(), z);
        continue;
      }
      auto particle = track.mcParticle();

      if (particle.isPhysicalPrimary())
      {
        registry.fill(HIST("Primary/TracksEtaZvtx"), track.eta(), z);
        registry.fill(HIST("Primary/DCAx"), retrack.bestDCAX());
        registry.fill(HIST("Primary/DCAy"), retrack.bestDCAY());
        registry.fill(HIST("Primary/DCAxy"), retrack.bestDCAXY());
      }
      else
      {
        registry.fill(HIST("Secondary/TracksEtaZvtx"), track.eta(), z);
        registry.fill(HIST("Secondary/DCAx"), retrack.bestDCAX());
        registry.fill(HIST("Secondary/DCAy"), retrack.bestDCAY());
        registry.fill(HIST("Secondary/DCAxy"), retrack.bestDCAXY());
      }

      registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
      registry.fill(HIST("DCAx"), retrack.bestDCAX());
      registry.fill(HIST("DCAy"), retrack.bestDCAY());
      registry.fill(HIST("DCAxy"), retrack.bestDCAXY());


    }
  }
  PROCESS_SWITCH(mftsecondary, processDCA, "get the DCAxy of MFT tracks in MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<mftsecondary>(cfgc)};
}
