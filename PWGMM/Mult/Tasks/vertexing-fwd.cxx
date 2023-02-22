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

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::track;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
using CollisionsLabeled = soa::Join<o2::aod::Collisions, aod::McCollisionLabels>;
using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedBCCollisionsSparseMulti>;

AxisSpec ZAxis = {301, -30.1, 30.1};

struct vertexingfwd {

  Service<o2::ccdb::BasicCCDBManager> ccdb;


  Configurable<float> maxDCAXY{"maxDCAXY", 6.0, "max allowed transverse DCA"}; // To be used when associating ambitrack to collision using best DCA
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  std::vector<uint64_t> ambTrackIds;

  int runNumber = -1;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  o2::parameters::GRPMagField* grpmag = nullptr;

  HistogramRegistry registry{
    "registry",
    {{"ParticleZR", "; #it{z}_{vtx} (cm); #it{R} (cm) ; count", {HistType::kTH2F, {ZAxis, {1001, -0.5, 100.5}}}},         //
     {"Primary/ParticleZR", "; #it{z}_{vtx} (cm); #it{R} (cm) ; count", {HistType::kTH2F, {ZAxis, {1001, -0.5, 100.5}}}}, //

     {"Truth/DCAxTruth", "; DCA_{x}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Truth/DCAyTruth", "; DCA_{y}^{truth} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Truth/DCAxyTruth", "; DCA_{xy}^{truth} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},

     {"Primary/DCAxPrimNAmb", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAyPrimNAmb", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAxPtPrimNAmb", "; DCA_{x} (cm); #it{p}_{T}^{true} (GeV/#it{c}) ; counts", {HistType::kTH2F, {{2500, -5, 5}, {21, -0.05, 10.05}}}},

     {"Primary/DCAxPrimAmb1", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAyPrimAmb1", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},

     {"Primary/DCAxPrimAmb2", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAyPrimAmb2", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{500, -5, 5}}}},
     {"Primary/DCAxyPrimAmb2", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},


     // ambiguous tracks and DCA with all compatible collisions
     {"Ambiguous/TracksDCAXY", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"Ambiguous/TracksDCAX", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"Ambiguous/TracksDCAY", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     //ambiguous tracks and DCA with true collision
     {"Ambiguous/TracksTrueDCAXY", "; DCA_{xy} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"Ambiguous/TracksTrueDCAX", "; DCA_{x} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"Ambiguous/TracksTrueDCAY", "; DCA_{y} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     {"Ambiguous/PrimOrSec/TracksTrueDCAXYPrim", "; DCA_{xy} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"Ambiguous/PrimOrSec/TracksTrueDCAXPrim", "; DCA_{x} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"Ambiguous/PrimOrSec/TracksTrueDCAYPrim", "; DCA_{y} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     {"Ambiguous/PrimOrSec/TracksTrueDCAXYSec", "; DCA_{xy} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"Ambiguous/PrimOrSec/TracksTrueDCAXSec", "; DCA_{x} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"Ambiguous/PrimOrSec/TracksTrueDCAYSec", "; DCA_{y} (cm) wrt true collision; counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     // all tracks (before amb reassociation)
     {"DCAXYAll", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"DCAXAll", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"DCAYAll", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     // ambiguous tracks before reassociation
     {"DCAXYAmb", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"DCAXAmb", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"DCAYAmb", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     // non ambiguous tracks
     {"DCAXYNAmb", "; DCA_{xy} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"DCAXNAmb", "; DCA_{x} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},
     {"DCAYNAmb", "; DCA_{y} (cm); counts", {HistType::kTH1F, {{1000, -10, 10}}}},

     {"AmbiguousTrackStatus", "; ; N_{trk}", {HistType::kTH1F, {{8, 0.5, 8.5}}}},

     // DCAxy, x and y distributions for reassociated ambiguous tracks
     // when it is a false reassociation and when it is true
     {"Ambiguous/TracksDCAXYBest", "; DCA_{xy}^{best} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"Ambiguous/TracksDCAXBest", "; DCA_{x}^{best} (cm); counts", {HistType::kTH1F, {{5000, -10, 10}}}},
     {"Ambiguous/TracksDCAYBest", "; DCA_{y}^{best} (cm); counts", {HistType::kTH1F, {{5000, -10, 10}}}},

     {"Ambiguous/TracksDCAXYBestTrue", "; DCA_{xy}^{best, true} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},
     {"Ambiguous/TracksDCAXYBestFalse", "; DCA_{xy}^{best, false} (cm); counts", {HistType::kTH1F, {{1000, -1, 10}}}},

     {"Ambiguous/TracksDCAXBestTrue", "; DCA_{x}^{best, true} (cm); counts", {HistType::kTH1F, {{5000, -10, 10}}}},
     {"Ambiguous/TracksDCAYBestTrue", "; DCA_{y}^{best, true} (cm); counts", {HistType::kTH1F, {{5000, -10, 10}}}},
     {"Ambiguous/TracksDCAXBestFalse", "; DCA_{x}^{best, false} (cm); counts", {HistType::kTH1F, {{5000, -10, 10}}}},
     {"Ambiguous/TracksDCAYBestFalse", "; DCA_{y}^{best, false} (cm); counts", {HistType::kTH1F, {{5000, -10, 10}}}}

    }};

  void init(InitContext&)
  {
    auto hstatus = registry.get<TH1>(HIST("AmbiguousTrackStatus"));
    auto* x2 = hstatus->GetXaxis();
    x2->SetBinLabel(1, "MFT tracks");
    x2->SetBinLabel(2, "MFT ambiguous tracks");//that also passed the cut on bestDCAxy
    x2->SetBinLabel(3, "Reassigned tracks");
    x2->SetBinLabel(4, "Extra tracks");
    x2->SetBinLabel(5, "orig=true (re)");
    x2->SetBinLabel(6, "best=true (re)");
    x2->SetBinLabel(7, "not reassigned");
    x2->SetBinLabel(8, "not reassigned and true");

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

  }

  void initCCDB(ExtBCs::iterator const& bc)
  {

    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current()
              << " A for run " << bc.runNumber()
              << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag);// for some reason this is necessary for the next next line
    runNumber = bc.runNumber();

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());

    Bz = field->getBz(centerMFT);//gives error if the propagator is not initFielded
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz;

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
    initCCDB(bcs.begin()); //if Bz is needed
    ambTrackIds.clear();

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
      int nTrueColl=0;
      std::vector<int> globalBCs;
      for (auto& bc : compatibleBCs) {
        globalBCs.push_back(bc.globalBC());
        if (!bc.has_collisions()) {
          continue;
        }
        auto collisions = bc.collisions_as<CollisionsLabeled>(); // compatible collisions
        for (auto const& collision : collisions) {

          trackPar.propagateToZhelix(collision.posZ(), Bz); // track parameters propagation to the position of the z vertex
          //trackPar.propagateToZlinear(collision.posZ());

          const auto dcaX(trackPar.getX() - collision.posX());
          const auto dcaY(trackPar.getY() - collision.posY());
          dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

          registry.fill(HIST("Ambiguous/TracksDCAXY"), dcaXY);
          registry.fill(HIST("Ambiguous/TracksDCAX"), dcaX);
          registry.fill(HIST("Ambiguous/TracksDCAY"), dcaY);
          // registry.fill(HIST("NumberOfContributors"), collision.numContrib());

          if(collision.mcCollisionId()==particle.mcCollisionId())
          {
            //we found the true collision of the ambiguous track, independently from its dca
            nTrueColl++;
            registry.fill(HIST("Ambiguous/TracksTrueDCAXY"), dcaXY);
            registry.fill(HIST("Ambiguous/TracksTrueDCAX"), dcaX);
            registry.fill(HIST("Ambiguous/TracksTrueDCAY"), dcaY);
            if(particle.isPhysicalPrimary())
            {
              registry.fill(HIST("Ambiguous/PrimOrSec/TracksTrueDCAXYPrim"), dcaXY);
              registry.fill(HIST("Ambiguous/PrimOrSec/TracksTrueDCAXPrim"), dcaX);
              registry.fill(HIST("Ambiguous/PrimOrSec/TracksTrueDCAYPrim"), dcaY);
            }
            else
            {
              registry.fill(HIST("Ambiguous/PrimOrSec/TracksTrueDCAXYSec"), dcaXY);
              registry.fill(HIST("Ambiguous/PrimOrSec/TracksTrueDCAXSec"), dcaX);
              registry.fill(HIST("Ambiguous/PrimOrSec/TracksTrueDCAYSec"), dcaY);
            }
          }

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
      if(nTrueColl<1)
      {
        if((particle.mcCollision().bc_as<ExtBCs>().globalBC() < globalBCs[globalBCs.size()-1]) && (particle.mcCollision().bc_as<ExtBCs>().globalBC() > globalBCs[0]))
        {

        }
        else
        {
          printf("No true collision for this track, the global bc range of the track is between %d - %d, true BC is %llu\n", globalBCs[0], globalBCs[globalBCs.size()-1], particle.mcCollision().bc_as<ExtBCs>().globalBC());
        }
      }
        //

      registry.fill(HIST("Ambiguous/TracksDCAXYBest"), bestDCA);
      registry.fill(HIST("Ambiguous/TracksDCAXBest"), bestDCAX);
      registry.fill(HIST("Ambiguous/TracksDCAYBest"), bestDCAY);

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("Primary/DCAxPrimAmb2"), bestDCAX);
        registry.fill(HIST("Primary/DCAyPrimAmb2"), bestDCAY);
        registry.fill(HIST("Primary/DCAxyPrimAmb2"), bestDCA);

      }

      if (bestDCA>maxDCAXY)
      {
        continue;
      }

      auto mcCollID = particle.mcCollisionId();
      registry.fill(HIST("AmbiguousTrackStatus"), 2);
      if (!track.has_collision()) {
        // this is a track extra
        registry.fill(HIST("AmbiguousTrackStatus"), 4);
        if (bestMCCol == mcCollID) // correctly assigned to bestCol
        {
          registry.fill(HIST("AmbiguousTrackStatus"), 6);
          registry.fill(HIST("Ambiguous/TracksDCAXYBestTrue"), bestDCA);
          registry.fill(HIST("Ambiguous/TracksDCAXBestTrue"), bestDCAX);
          registry.fill(HIST("Ambiguous/TracksDCAYBestTrue"), bestDCAY);
        } else {
          registry.fill(HIST("Ambiguous/TracksDCAXYBestFalse"), bestDCA);
          registry.fill(HIST("Ambiguous/TracksDCAXBestFalse"), bestDCAX);
          registry.fill(HIST("Ambiguous/TracksDCAYBestFalse"), bestDCAY);
        }

      } else if (track.collisionId() != bestCol) {
        // this track has been reassigned
        auto collOrig = track.collision_as<CollisionsLabeled>();
        registry.fill(HIST("AmbiguousTrackStatus"), 3);
        if (bestMCCol == mcCollID) // correctly reassigned
        {
          registry.fill(HIST("AmbiguousTrackStatus"), 6);
          registry.fill(HIST("Ambiguous/TracksDCAXYBestTrue"), bestDCA);
          registry.fill(HIST("Ambiguous/TracksDCAXBestTrue"), bestDCAX);
          registry.fill(HIST("Ambiguous/TracksDCAYBestTrue"), bestDCAY);
        } else { // uncorrectly reassigned
          registry.fill(HIST("Ambiguous/TracksDCAXYBestFalse"), bestDCA);
          registry.fill(HIST("Ambiguous/TracksDCAXBestFalse"), bestDCAX);
          registry.fill(HIST("Ambiguous/TracksDCAYBestFalse"), bestDCAY);
        }

        if (collOrig.mcCollisionId() == mcCollID) { // initially correctly assigned
          registry.fill(HIST("AmbiguousTrackStatus"), 5);
        }
      } else // the track has a collision and track.collisionId() == bestCol
      {
        if (track.collisionId() != bestCol) {
          printf("------------------- PROBLEM HERE track.collisionId() %d, bestCollid %d\n", track.collisionId(), bestCol);
        }

        registry.fill(HIST("AmbiguousTrackStatus"), 7);
        if (bestMCCol == mcCollID) // correctly assigned
        {
          registry.fill(HIST("AmbiguousTrackStatus"), 8);
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
    registry.fill(HIST("AmbiguousTrackStatus"), 1, tracks.size());
    float dcaXY;
    o2::track::TrackParCovFwd bestTrackPar;

    for (auto& track : tracks) {
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      if (!track.has_collision()) {
        continue;
      }
      auto collision = track.collision();

      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      //trackPar.propagateToZlinear(collision.posZ());
      trackPar.propagateToZhelix(collision.posZ(), Bz);

      const auto dcaX(trackPar.getX() - collision.posX());
      const auto dcaY(trackPar.getY() - collision.posY());
      dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);

      float Rv = std::sqrt(pow(particle.vx(), 2) + pow(particle.vy(), 2));
      registry.fill(HIST("ParticleZR"), particle.vz(), Rv);

      if (particle.isPhysicalPrimary()) {
        registry.fill(HIST("Primary/ParticleZR"), particle.vz(), Rv);
      }

      const auto dcaXtruth(particle.vx() - particle.mcCollision().posX());
      const auto dcaYtruth(particle.vy() - particle.mcCollision().posY());
      auto dcaXYtruth = std::sqrt(dcaXtruth * dcaXtruth + dcaYtruth * dcaYtruth); // this is DCA_xy truth

      registry.fill(HIST("Truth/DCAxTruth"), dcaXtruth);
      registry.fill(HIST("Truth/DCAyTruth"), dcaYtruth);
      registry.fill(HIST("Truth/DCAxyTruth"), dcaXYtruth);

      registry.fill(HIST("DCAXYAll"), dcaXY);
      registry.fill(HIST("DCAXAll"), dcaX);
      registry.fill(HIST("DCAYAll"), dcaY);

      if (find(ambTrackIds.begin(), ambTrackIds.end(), track.globalIndex()) != ambTrackIds.end()) {
        registry.fill(HIST("DCAXYAmb"), dcaXY);
        registry.fill(HIST("DCAXAmb"), dcaX);
        registry.fill(HIST("DCAYAmb"), dcaY);
        if (particle.isPhysicalPrimary()) {
          registry.fill(HIST("Primary/DCAxPrimAmb1"), dcaX);
          registry.fill(HIST("Primary/DCAyPrimAmb1"), dcaY);
        }

        continue; // this track has already been reassigned to bestcollision, don't double count
      }
      registry.fill(HIST("DCAXYNAmb"), dcaXY);
      registry.fill(HIST("DCAXNAmb"), dcaX);
      registry.fill(HIST("DCAYNAmb"), dcaY);

      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      // only tracks coming from primary particles here
      // non ambiguous tracks
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
