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

AxisSpec PtAxis = {1001, -0.005, 10.005};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec PhiAxis = {600, 0, 2 * M_PI};
AxisSpec EtaAxis = {18, -4.6, -1.};
AxisSpec DCAxyAxis = {100, -1, 10};
AxisSpec CentAxis = {{0, 10, 20, 30, 40, 50, 60, 70, 80, 100}};

int nVis=0;

static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype trackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;
using MFTTracksWColls = soa::Join<o2::aod::MFTTracks, aod::BestCollisionsFwd>;

struct PseudorapidityDensityMFT {
  SliceCache cache;
  Preslice<aod::MFTTracks> perCol = o2::aod::fwdtrack::collisionId;
  //PresliceUnsorted<MFTTracksWColls> perColBis = o2::aod::fwdtrack::bestCollisionId;
  PresliceUnsorted<aod::BestCollisionsFwd> perColBis = o2::aod::fwdtrack::bestCollisionId;
  Preslice<aod::McParticles> perMcCol = aod::mcparticle::mcCollisionId;
  Preslice<aod::Tracks> perColCentral = aod::track::collisionId;

  Service<O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};

  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};

  Configurable<bool> useZDiffCut{"useZDiffCut", true, "use Z difference cut"};
  Configurable<float> maxZDiff{"maxZDiff", 1.0f, "max allowed Z difference for reconstruced collisions (cm)"};
  //-----need access to CCDB to get the reweighting histogram
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> path{"ccdb-path", "Users/s/sherrman/My/Object", "base path to the ccdb object"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // the histogram has been previously stored in the CCDB
  TH1D* histoReweight = nullptr;
  std::vector<uint64_t> ambTrackIds;
  int counter = 0;
  //------

  HistogramRegistry registry{
    "registry",
    {
      //{"EventsNtrkZvtx", "; N_{trk}; #it{z}_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}}, //
      {"TracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},      //
      {"Tracks/EtaZvtx_gt0", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}}, //
      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}},               //
      {"TracksPhiZvtx", "; #varphi; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {PhiAxis, ZAxis}}},   //
      {"TracksPtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}},                  //
      {"EventSelection", ";status;events", {HistType::kTH1F, {{7, 0.5, 7.5}}}},                         //
      {"Tracks/Control/TrackAmbDegree", " ; N_{coll}^{comp}", {HistType::kTH1F, {{21,-0.5,20.5}}}},
      {"Tracks/Control/TrackIsAmb", " ; isAmbiguous", {HistType::kTH1F, {{2,-0.5,1.5}}}},
      {"Tracks/Control/DCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {DCAxyAxis}}},                                                    //
                                                                                                                                       //{"Tracks/Control/DCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}},                     //
      {"Tracks/Control/ReassignedDCAXY", "  ; DCA_{XY} (cm)", {HistType::kTH1F, {DCAxyAxis}}},                                         //
                                                                                                                                       //{"Tracks/Control/ReassignedDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAxyAxis}}},           //
      {"Tracks/Control/ExtraDCAXY", " ; DCA_{XY} (cm)", {HistType::kTH1F, {DCAxyAxis}}},                                               //
                                                                                                                                       //{"Tracks/Control/ExtraDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAxyAxis}}},                //
      {"Tracks/Control/ExtraTracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},                 //
      {"Tracks/Control/ExtraTracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}},                          //
      {"Tracks/Control/ReassignedTracksEtaZvtx", "; #eta; #it{z}_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}},            //
      {"Tracks/Control/ReassignedTracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}},                     //
      {"Tracks/Control/ReassignedVertexCorr", "; #it{z}_{vtx}^{orig} (cm); #it{z}_{vtx}^{re} (cm)", {HistType::kTH2F, {ZAxis, ZAxis}}} //
    }                                                                                                                                  //
  };

  void init(InitContext&)
  {
    AxisSpec MultAxis = {multBinning, "N_{trk}"}; // for PbPb 3001,-0.5,3000.5

    auto hstat = registry.get<TH1>(HIST("EventSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All");
    x->SetBinLabel(2, "Selected");
    x->SetBinLabel(3, "Selected INEL>0");
    x->SetBinLabel(4, "Rejected");
    x->SetBinLabel(5, "Good BCs");
    x->SetBinLabel(6, "BCs with collisions");
    x->SetBinLabel(7, "BCs with pile-up/splitting");

    registry.add({"EventsNtrkZvtx", "; N_{trk}; #it{z}_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
    registry.add({"EventsNtrkZvtx_gt0", "; N_{trk}; #it{z}_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  void processTagging(FullBCs const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {

    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (useEvSel && ((bc.selection_bit(evsel::kIsBBT0A) & bc.selection_bit(evsel::kIsBBT0C)) != 0))) {
        registry.fill(HIST("EventSelection"), 5.);
        cols.clear();
        for (auto& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registry.fill(HIST("EventSelection"), 6.);
          if (cols.size() > 1) {
            registry.fill(HIST("EventSelection"), 7.);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processTagging, "Collect event sample stats", true);




  Partition<MFTTracksWColls> sample = (aod::fwdtrack::eta < -2.8f) && (aod::fwdtrack::eta > -3.2f);

  Partition<aod::MFTTracks> sample1 = (aod::fwdtrack::eta < -2.8f) && (aod::fwdtrack::eta > -3.2f);

  expressions::Filter atrackFilter = (aod::fwdtrack::bestCollisionId >= 0) &&
                                     (aod::fwdtrack::eta < -2.0f) &&
                                     (aod::fwdtrack::eta > -3.9f) &&
                                     (nabs(aod::fwdtrack::bestDCAXY) <= 2.f);

  expressions::Filter trackSelectionCentral = ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) &&
                                              ifnode((aod::track::detectorMap & (uint8_t)o2::aod::track::TPC) == (uint8_t)o2::aod::track::TPC,
                                              (aod::track::trackCutFlag & trackSelectionTPC) == trackSelectionTPC,
                                              true) &&
                                              ((aod::track::trackCutFlag & trackSelectionDCA) == trackSelectionDCA) &&
                                              (nabs(aod::track::eta) < estimatorEta);

  using FiCentralTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>; // central tracks for INEL>0

  using CollwEv = soa::Join<aod::Collisions, aod::EvSels>;

  void processTest(CollwEv::iterator const& collision,
                   o2::aod::MFTTracks const&,
                   soa::SmallGroups<aod::BestCollisionsFwd> const& atracks)
  {
    if (atracks.size() == 0) {
      return;
    }
    registry.fill(HIST("EventSelection"), 1.);
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();

      for (auto& track : atracks){
        auto normalTrack = track.mfttrack();
        int isAmbiguous=0;

        registry.fill(HIST("TracksEtaZvtx"), normalTrack.eta(), z);
        if (track.ambDegree()>2)
        {
          //only ambiguous tracks here
          isAmbiguous=1;
        }
        registry.fill(HIST("Tracks/Control/TrackIsAmb"), isAmbiguous);

      }
    }
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processTest, "Process test", true);

  void processMultReassoc(CollwEv const& collisions,
                   MFTTracksWColls const& tracks, FiCentralTracks const& centraltracks)
  {

    for (auto& collision : collisions)
    {
      auto perCollisionMFTTracks = tracks.sliceByCached(aod::fwdtrack::bestCollisionId, collision.globalIndex(), cache);
      if (perCollisionMFTTracks.size()==0){continue;}
      registry.fill(HIST("EventSelection"), 1.);
      if (!useEvSel || (useEvSel && collision.sel8())) {
        registry.fill(HIST("EventSelection"), 2.);
        auto z = collision.posZ();
        auto perCollisionSample = sample->sliceByCached(aod::fwdtrack::bestCollisionId, collision.globalIndex(), cache);
        auto Ntrk = perCollisionSample.size();
        auto midtracks = centraltracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

        registry.fill(HIST("EventsNtrkZvtx"), Ntrk, z);

        for (auto& track : perCollisionMFTTracks) {

          registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
          if (midtracks.size() > 0) // INEL>0
          {
            registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
          }
          float phi = track.phi();
          o2::math_utils::bringTo02Pi(phi);
          registry.fill(HIST("TracksPhiEta"), phi, track.eta());
          registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
          if ((track.eta() < -2.0f) && (track.eta() > -3.9f)) {
            registry.fill(HIST("TracksPhiZvtx"), phi, z);
          }

          //more control histograms for ambiguous tracks:
          registry.fill(HIST("Tracks/Control/TrackAmbDegree"), track.ambDegree());

          int isAmbiguous = 0;

          if ((track.globalIndex()<35) && (track.globalIndex()>22) && (nVis==0))
          {
            //printf("----------------- track.globalIndex %lld, bestCol %lld, track.collisionId %lld \n", track.globalIndex(), track.bestCollisionId(), track.collisionId());
          }

          if (track.ambDegree()>2)
          {
            //only ambiguous tracks here
            isAmbiguous=1;
            if (track.collisionId() != track.bestCollisionId()) {
              //reassigned tracks here
              registry.fill(HIST("Tracks/Control/ReassignedTracksEtaZvtx"), track.eta(), z);
              registry.fill(HIST("Tracks/Control/ReassignedTracksPhiEta"), phi, track.eta());
              registry.fill(HIST("Tracks/Control/ReassignedVertexCorr"), track.collision_as<CollwEv>().posZ(), z);
              registry.fill(HIST("Tracks/Control/ReassignedDCAXY"), track.pt(), track.bestDCAXY());
              // registry.fill(HIST("Tracks/Control/ReassignedDCAZPt"), otrack.pt(), track.bestDCAZ());
            }
          }
          registry.fill(HIST("Tracks/Control/TrackIsAmb"), isAmbiguous);

        }//end of (auto& track : perCollisionMFTTracks)
        if (midtracks.size() > 0) // INEL>0
        {
          registry.fill(HIST("EventSelection"), 3.);
          registry.fill(HIST("EventsNtrkZvtx_gt0"), Ntrk, z);
        }

      } else {
        registry.fill(HIST("EventSelection"), 4.);
        }
      }
      nVis++;
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processMultReassoc, "Process reco or data info with reassociation", true);


  void processMultTest(CollwEv const& collisions,
                   aod::MFTTracks const& tracks,
                   FiCentralTracks const& midtracks)
  {

    for (auto& collision : collisions)
    {

      auto perCollisionMFTTracks = tracks.sliceByCached(aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      if (perCollisionMFTTracks.size()==0){continue;}



      registry.fill(HIST("EventSelection"), 1.);
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();
      auto perCollisionSample = sample1->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto Ntrk = perCollisionSample.size();



      registry.fill(HIST("EventsNtrkZvtx"), Ntrk, z);

      for (auto& track : perCollisionMFTTracks) {
        registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
        if (midtracks.size() > 0) // INEL>0
        {
          registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
        }
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("TracksPhiEta"), phi, track.eta());
        registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
        if ((track.eta() < -2.0f) && (track.eta() > -3.9f)) {
          registry.fill(HIST("TracksPhiZvtx"), phi, z);
        }
      }
      if (midtracks.size() > 0) // INEL>0
      {
        registry.fill(HIST("EventSelection"), 3.);
        registry.fill(HIST("EventsNtrkZvtx_gt0"), Ntrk, z);
      }

    } else {
      registry.fill(HIST("EventSelection"), 4.);
    }
  }
  }
  PROCESS_SWITCH(PseudorapidityDensityMFT, processMultTest, "Process reco or data info", false);

  void processMultTest2(CollwEv::iterator const& collision,
                   aod::MFTTracks const& tracks,
                   FiCentralTracks const& midtracks)
  {
      if (tracks.size()==0){return;}



    registry.fill(HIST("EventSelection"), 1.);
    if (!useEvSel || (useEvSel && collision.sel8())) {
      registry.fill(HIST("EventSelection"), 2.);
      auto z = collision.posZ();
      auto perCollisionSample = sample1->sliceByCached(o2::aod::fwdtrack::collisionId, collision.globalIndex(), cache);
      auto Ntrk = perCollisionSample.size();



      registry.fill(HIST("EventsNtrkZvtx"), Ntrk, z);

      for (auto& track : tracks) {
        registry.fill(HIST("TracksEtaZvtx"), track.eta(), z);
        if (midtracks.size() > 0) // INEL>0
        {
          registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
        }
        float phi = track.phi();
        o2::math_utils::bringTo02Pi(phi);
        registry.fill(HIST("TracksPhiEta"), phi, track.eta());
        registry.fill(HIST("TracksPtEta"), track.pt(), track.eta());
        if ((track.eta() < -2.0f) && (track.eta() > -3.9f)) {
          registry.fill(HIST("TracksPhiZvtx"), phi, z);
        }
      }
      if (midtracks.size() > 0) // INEL>0
      {
        registry.fill(HIST("EventSelection"), 3.);
        registry.fill(HIST("EventsNtrkZvtx_gt0"), Ntrk, z);
      }

    } else {
      registry.fill(HIST("EventSelection"), 4.);
    }
  }

  PROCESS_SWITCH(PseudorapidityDensityMFT, processMultTest2, "Process reco or data info", false);

};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudorapidityDensityMFT>(cfgc)};
}
