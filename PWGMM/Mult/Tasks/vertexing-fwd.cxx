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

// \file   vertexing-fwd.cxx
// \author Robin Caron <robin.caron@cern.ch>
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief This code loops over every ambiguous MFT tracks and associates
// them to a collision that has the smallest DCAxy

#include <cmath>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MathUtils/Utils.h"
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::track;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
using CollisionsLabeled= soa::Join<o2::aod::Collisions, aod::McCollisionLabels>;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

AxisSpec ZAxis = {301, -30.1, 30.1};

struct vertexingfwd {

  Configurable<float> maxDCAXY{"maxDCAXY", 6.0, "max allowed transverse DCA"}; // To be used when associating ambitrack to collision using best DCA
  std::vector<long long unsigned> ambTrackIds;

  HistogramRegistry registry
  {
    "registry",
    {
      {"ParticleZR", "; #it{z}_{vtx} (cm); #it{R} (cm) ; count", {HistType::kTH2F, {ZAxis, {1001, -0.5, 100.5}}}}, //

      {"Truth/DCAxTruth", "; DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"Truth/DCAyTruth", "; DCA_{y}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

      {"Primary/DCAxPrimNAmb", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"Primary/DCAyPrimNAmb", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"Primary/DCAxPtPrimNAmb", "; DCA_{x} (cm); #it{p}_{T}^{true} (GeV/#it{c}) ; counts", {HistType::kTH2F, {{2500, -5, 5}, {21, -0.05, 10.05}}}},

      {"Primary/DCAxPrimAmb1", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"Primary/DCAyPrimAmb1", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

      {"Primary/DCAxPrimAmb2", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"Primary/DCAyPrimAmb2", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

      {"TracksDCAXY", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
      {"DeltaDCAXY", "; DCA_{xy}^{reco}-DCA_{xy}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAXYAll", "; DCA_{xy}^{reco}-DCA_{xy}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"TracksDCAX", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{100, -10, 10}}}},
      {"TracksDCAY", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{100, -10, 10}}}},
      {"DeltaDCAX", "; DCA_{x}^{reco}-DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAY", "; DCA_{y}^{reco}-DCA_{y}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAXTrue", "; DCA_{x}^{reco}-DCA_{x}^{truth} (cm) when reassociation is good; counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAYTrue", "; DCA_{y}^{reco}-DCA_{y}^{truth} (cm) when reassociation is good; counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAXAll", "; DCA_{x}^{reco}-DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAYAll", "; DCA_{y}^{reco}-DCA_{y}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

      {"DeltaDCAXYNAmb", "; DCA_{x}^{reco}-DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAYNAmb", "; DCA_{y}^{reco}-DCA_{y}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAXNAmb", "; DCA_{x}^{reco}-DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

      {"DeltaDCAXYAmb", "; DCA_{x}^{reco}-DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAYAmb", "; DCA_{y}^{reco}-DCA_{y}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
      {"DeltaDCAXAmb", "; DCA_{x}^{reco}-DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

      {"AmbiguousTracksStatus", "; Status; counts", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
      {"NbCollComp", "; NbCollComp", {HistType::kTH1F, {{20, -0.5, 19.5}}}},
      {"NumberOfContributors", "; N_{tr} for vertexing; counts", {HistType::kTH1F, {{100, 0, 100}}}},
      {"CollisionsMatchIndicesMC", "; Rec. minDCA ambitrack coll.ID; Gen. ambitrack coll.ID", {HistType::kTH2F, {{401, -0.5, 1000.5}, {401, -0.5, 1000.5}}}},
      {"TracksDCAXYBest", "; DCA_{xy}^{best} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
      {"TracksDCAXYBestFalse", "; DCA_{xy}^{best, false} (cm); counts", {HistType::kTH1F, {{100, -1, 10}}}},
      {"EfficiencyZvtx", "; z vertex; Efficiency", {HistType::kTProfile, {{100, -30, 30}}}},
      {"DeltaZvtx", "; #delta z (cm); counts", {HistType::kTH1F, {{400, -20, 20}}}},
      {"DeltaZvtxBest", "; #delta z = z_{best} - z_{true} (cm); counts", {HistType::kTH1F, {{400, -20, 20}}}},
      {"CorrectMatch", "; Matching value; counts", {HistType::kTH1F, {{5, -0.5, 4.5}}}}
    }
  };

  void init(InitContext&)
  {
    auto hstat = registry.get<TH1>(HIST("CorrectMatch"));
    auto* x1 = hstat->GetXaxis();
    x1->SetBinLabel(1, "Incorrect match");
    x1->SetBinLabel(2, "Correct match");
    x1->SetBinLabel(3, "N_{ambitrack} associable");
    x1->SetBinLabel(4, "N_{ambitrack} w N_{coll} > 0");
    x1->SetBinLabel(5, "N_{ambitrack} total");

    auto hstatus = registry.get<TH1>(HIST("AmbiguousTracksStatus"));
    auto* x2 = hstatus->GetXaxis();
    x2->SetBinLabel(1, "MFT tracks ");
    x2->SetBinLabel(2, "MFT ambiguous tracks ");
    x2->SetBinLabel(3, "All ambiguous primary");
    x2->SetBinLabel(4, "Orphan + primary");
    x2->SetBinLabel(5, "All ambiguous secondary");
    x2->SetBinLabel(6, "Orphan + secondary");
  }

  void processDCAamb(MFTTracksLabeled const&,
                  CollisionsLabeled const&, ExtBCs const& bcs,
                  aod::AmbiguousMFTTracks const& atracks,
                  aod::McParticles const&,
                  aod::McCollisions const&)
  {

    if (bcs.size() == 0) {
      return;
    }
    if (atracks.size() == 0) {
      return;
    }
    //initCCDB(bcs.begin()); if Bz is needed
    ambTrackIds.clear();

    // Only on DCAxy
    float dcaXY;
    float bestDCA, bestDCAX, bestDCAY;
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto& atrack : atracks) {

      dcaXY = 999; // DCAxy
      bestDCA = 999;

      auto track = atrack.mfttrack_as<MFTTracksLabeled>();
      if (track.has_collision()) {
        ambTrackIds.push_back(atrack.mfttrackId());
      }
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        continue;
      }
      auto particle = track.mcParticle();

      auto bestCol = track.has_collision() ? track.collisionId() : -1;
      int bestMCCol = -1;
      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      auto compatibleBCs = atrack.bc_as<ExtBCs>();

      for (auto& bc : compatibleBCs) {
        if (!bc.has_collisions()) {
          continue;
        }
        auto collisions = bc.collisions_as<CollisionsLabeled>();//compatible collisions
        for (auto const& collision : collisions) {

          //trackPar.propagateToZhelix(collision.posZ(), Bz); // track parameters propagation to the position of the z vertex
          trackPar.propagateToZlinear(collision.posZ());

          const auto dcaX(trackPar.getX() - collision.posX());
          const auto dcaY(trackPar.getY() - collision.posY());
          dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);



          registry.fill(HIST("TracksDCAXY"), dcaXY);
          registry.fill(HIST("TracksDCAX"), dcaX);
          registry.fill(HIST("TracksDCAY"), dcaY);
          registry.fill(HIST("NumberOfContributors"), collision.numContrib());



          if ((dcaXY < bestDCA)) {
            bestCol = collision.globalIndex();
            bestMCCol = collision.mcCollisionId();
            bestDCA = dcaXY;
            bestDCAX = dcaX;
            bestDCAY = dcaY;
            bestTrackPar = trackPar;
          }
        }
      }

      // other option for the truth : collision.mcCollision().posZ();

      const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
      const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
      auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth);//this is DCA_xy truth


      registry.fill(HIST("DeltaDCAXY"), bestDCA-dcaXYtruth);
      registry.fill(HIST("DeltaDCAX"), bestDCAX-dcaXtruth);
      registry.fill(HIST("DeltaDCAY"), bestDCAY-dcaYtruth);

      if (particle.isPhysicalPrimary())
      {
        registry.fill(HIST("Primary/DCAxPrimAmb2"), bestDCAX);
        registry.fill(HIST("Primary/DCAyPrimAmb2"), bestDCAY);
      }

      if (particle.mcCollisionId()==bestMCCol)
      {
        //good reassociation
        registry.fill(HIST("DeltaDCAXTrue"), bestDCAX-dcaXtruth);
        registry.fill(HIST("DeltaDCAYTrue"), bestDCAY-dcaYtruth);
      }
      if(track.has_collision())
      {
        auto collision1 = track.collision_as<CollisionsLabeled>();
        if (particle.mcCollisionId()==collision1.mcCollisionId())
        {
          // registry.fill(HIST("DeltaDCAXTrue1"), bestDCAX-dcaXtruth);
          // registry.fill(HIST("DeltaDCAYTrue1"), bestDCAY-dcaYtruth);
          //calculate the DCA for coll1
        }
      }


    }
  }
PROCESS_SWITCH(vertexingfwd, processDCAamb, "get the DCAxy of MFT ambiguous tracks", true);

void processDCA(MFTTracksLabeled const& tracks,
                aod::Collisions const&,
                aod::McParticles const&,
                aod::McCollisions const&)
{
  if (tracks.size() == 0) {
    return;
  }

  float dcaXY;
  o2::track::TrackParCovFwd bestTrackPar;

  for (auto& track : tracks) {
    if (!track.has_mcParticle()) {
      LOGF(warning, "No MC particle for ambiguous track, skip...");
      continue;
    }
    auto particle = track.mcParticle();
    if (!track.has_collision()){continue;}
    auto collision = track.collision();

    std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
    SMatrix55 tcovs(v1.begin(), v1.end());
    SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

    trackPar.propagateToZlinear(collision.posZ());

    const auto dcaX(trackPar.getX() - collision.posX());
    const auto dcaY(trackPar.getY() - collision.posY());
    dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

    float Rv = std::sqrt(pow(particle.vx(),2)+pow(particle.vy(),2));
    registry.fill(HIST("ParticleZR"), particle.vz(), Rv);

    const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
    const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
    auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth);//this is DCA_xy truth


    registry.fill(HIST("Truth/DCAxTruth"), dcaXtruth);
    registry.fill(HIST("Truth/DCAyTruth"), dcaYtruth);

    registry.fill(HIST("DeltaDCAXYAll"), dcaXY-dcaXYtruth);
    registry.fill(HIST("DeltaDCAXAll"), dcaX-dcaXtruth);
    registry.fill(HIST("DeltaDCAYAll"), dcaY-dcaYtruth);

    if (find(ambTrackIds.begin(), ambTrackIds.end(), track.globalIndex()) != ambTrackIds.end()) {
      registry.fill(HIST("DeltaDCAXYAmb"), dcaXY-dcaXYtruth);
      registry.fill(HIST("DeltaDCAXAmb"), dcaX-dcaXtruth);
      registry.fill(HIST("DeltaDCAYAmb"), dcaY-dcaYtruth);
      if (particle.isPhysicalPrimary())
      {
        registry.fill(HIST("Primary/DCAxPrimAmb1"), dcaX);
        registry.fill(HIST("Primary/DCAyPrimAmb1"), dcaY);
      }


      continue; // this track has already been reassigned to bestcollision, don't double count
    }
    registry.fill(HIST("DeltaDCAXYNAmb"), dcaXY-dcaXYtruth);
    registry.fill(HIST("DeltaDCAXNAmb"), dcaX-dcaXtruth);
    registry.fill(HIST("DeltaDCAYNAmb"), dcaY-dcaYtruth);

    if (!particle.isPhysicalPrimary()) {
      continue;
    }
    //only tracks coming from primary particles here
    //non ambiguous tracks
    registry.fill(HIST("Primary/DCAxPrimNAmb"), dcaX);
    registry.fill(HIST("Primary/DCAyPrimNAmb"), dcaY);

    registry.fill(HIST("Primary/DCAxPtPrimNAmb"), dcaX, particle.pt());

  }
}
PROCESS_SWITCH(vertexingfwd, processDCA, "get the DCAxy of MFT tracks", true);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexingfwd>(cfgc)};
}
