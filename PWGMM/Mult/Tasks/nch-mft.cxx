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

#include "Common/DataModel/EventSelection.h"

#include "MathUtils/Utils.h"
#include "CommonConstants/LHCConstants.h"

#include "Common/DataModel/CollisionAssociationTables.h"
#include <experimental/type_traits>

using namespace o2;
using namespace o2::framework;

using CollwEv = aod::Collisions;

AxisSpec PhiAxis = {600, -2*M_PI, 2*M_PI};
AxisSpec EtaAxis = {180, -8.6, 8.6};

struct nchMft {

  Configurable<int> nBCTimeShift{"nBCTimeShift", 0, "shift in BC wrt to the center of the ROF (calculated with FT0C)"};
  int nDataFile = 0;

  HistogramRegistry registry{
    "registry",
    {

    }
  };

  void init(o2::framework::InitContext& initContext)
  {
    if (doprocessAmb)
    {
      registry.add({"MFTTrackAmbDegree", " ; N_{coll}^{comp}", {HistType::kTH1I, {{41, -0.5, 40.5}}}});
    }

    if (doprocessMFT)
    {
      registry.add({"MFT/MFTTrackhasNoColl", "; hasNoColl; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFT/MFTTrackhasNoCollComp", "; hasNoCollComp; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFT/bcDiffCompCollMiddleROF", "; bcCompColl - bcMiddleROF; count", {HistType::kTH1F, {{2001, -1000.5, 1000.5}}}});
    }
  }

  //using BCwColls = soa::Join<aod::BCs, aod::MatchedBCCollisionsSparseMulti>;
  using BCwColls = aod::BCs;

  using MFTTracksWithColls = soa::Join<aod::MFTTracks, aod::MFTTrkCompColls>;

  void processMFT(MFTTracksWithColls const& mfttracks, CollwEv const& collisions, BCwColls const&)
  {
    for (auto& mfttrack : mfttracks)
    {

      // if (!(mfttrack.has_bcs() && mfttrack.has_compatibleColl() && mfttrack.has_collision()))//mft tracks having a match in FT0-C
      // {
      //   continue;
      // }

      bool has_coll=true;
      bool has_compColl=true;

      registry.fill(HIST("MFT/MFTTrackhasNoCollComp"), 0);//counter for all the mft tracks
      registry.fill(HIST("MFT/MFTTrackhasNoColl"), 0);//counter for all the mft tracks

      if(!mfttrack.has_collision())
      {
        //fill the orphan track histo
        registry.fill(HIST("MFT/MFTTrackhasNoColl"), 1);//orphan track before reassoc tool
        has_coll=false;
      }

      if(!mfttrack.has_compatibleColl())
      {
        //fill an histo dedicated to no comp coll (or not, fill an histo with the number of compColls)
        registry.fill(HIST("MFT/MFTTrackhasNoCollComp"), 1);
        has_compColl=false;
      }

      if(!(has_coll && has_compColl))
      {
        continue;
      }

      //the BC of the middle of the MFT ROF
      int64_t middleBCWrong = nBCTimeShift+mfttrack.collision_as<CollwEv>().bc_as<BCwColls>().globalBC() + mfttrack.trackTime()/o2::constants::lhc::LHCBunchSpacingNS;
      int64_t middleBC = mfttrack.collision_as<CollwEv>().bc_as<BCwColls>().globalBC() + mfttrack.trackTime()/o2::constants::lhc::LHCBunchSpacingNS;
      //the collisions compatible to this mfttrack according to the fwdtrack-to-collision-associator
      auto compatibleColls = mfttrack.compatibleColl_as<CollwEv>();

      for (auto& collMFT : compatibleColls)
      {//looking at all the collisions compatible with the MFT track

        int64_t bcDiffCompCollMiddleROF = collMFT.bc_as<BCwColls>().globalBC() - middleBC;
        registry.fill(HIST("MFT/bcDiffCompCollMiddleROF"), bcDiffCompCollMiddleROF);
        // if (std::abs(bcDiffCompCollMiddleROF)>460)
        // {
        //   printf("nDataFile = %d, trackId = %lld, compCollIdx = %lld\n", nDataFile, mfttrack.globalIndex(), collMFT.globalIndex());
        // }

        if (nDataFile == 3 && mfttrack.globalIndex() == 11353)
        {
          printf("--------------------------------- CHECKING IN NCH-MFT -------\n");
          printf("middle BC = %lld, bc of comp coll %lld\n", middleBC, collMFT.bc_as<BCwColls>().globalBC());
          printf("middle BC wrong = %lld\n", middleBCWrong);
          int fBCTimeShift=nBCTimeShift;
          printf("nBCTimeShift %d\n", fBCTimeShift);
          printf("LHCBunchSpacingNS %f\n", o2::constants::lhc::LHCBunchSpacingNS );
          printf("its coll ID = %lld \n", mfttrack.collision_as<CollwEv>().globalIndex());
          printf("---------- trackTime = %f,  bc of ITS coll = %lld\n", mfttrack.trackTime(), mfttrack.collision().bc().globalBC());
          printf("track eta = %f, track phi %f\n", mfttrack.eta(), mfttrack.phi());
        }
      }//for (auto& collMFT : compatibleColls)


    }//for (auto& mfttrack : mfttracks)

    nDataFile++;
  }
  PROCESS_SWITCH(nchMft, processMFT, "processMFT", true);

  // expressions::Filter mfttrackFilter = (aod::fwdtrack::eta < -2.0f) &&
  //                                    (aod::fwdtrack::eta > -3.9f);
  //
  // using FiMFTTracks = soa::Filtered<o2::aod::MFTTracks>;


  // template <class T>
  // using hasTrackType = decltype(std::declval<T&>().trackType());


  void processAmb(MFTTracksWithColls::iterator const& track)
  {
    auto compatibleColls = track.compatibleCollIds();
    registry.fill(HIST("MFTTrackAmbDegree"), compatibleColls.size());

    // if constexpr (std::experimental::is_detected<hasTrackType, MFTTracksWithColls>::value)
    // {
    //   printf("-------- I have a track type, in MFT process %d\n", 1);
    // }
    //
    // if constexpr (std::is_same<MFTTracksWithColls, o2::aod::MFTTracks>::value)
    // {
    //
    // }

    if constexpr (MFTTracksWithColls::contains<o2::aod::MFTTracks>())
    {
      //then we are dealing with MFT tracks do somthing
      //printf("-------- !!! I have an MFT track %d\n", 1);
    }

    if constexpr (o2::aod::MFTTracks::contains<o2::aod::MFTTracks>())
    {
      //then we are dealing with MFT tracks do somthing
      //printf("-------- !!! I have an MFT track %d\n", 2);
    }
  }
  PROCESS_SWITCH(nchMft, processAmb, "find the percentage of MFT ambiguous tracks", true);

  // void processFilteredMFT(FiMFTTracks const& tracks)
  // {
  //
  //   if constexpr (FiMFTTracks::contains<o2::aod::MFTTracks>())
  //   {
  //     //then we are dealing with MFT tracks do somthing
  //     printf(">>>>>>>>>>>> I have an MFT track FILTERED %d\n", 4);
  //   }
  //
  //
  // }
  // PROCESS_SWITCH(nchMft, processFilteredMFT, "look at filtered MFT tracks", true);


  // void processFwd(o2::aod::FwdTracks const& fwdtracks)
  // {
  //   if constexpr (std::experimental::is_detected<hasTrackType, o2::aod::FwdTracks>::value)
  //   {
  //     printf(">>>>>>>> I have a track type, in fwd process %d\n", 1);
  //   }
  //   else
  //   {
  //     //printf(">>>>>>>> No track type, in fwd process %d\n", 1);
  //   }
  //
  //   if constexpr (o2::aod::FwdTracks::contains<o2::aod::MFTTracks>())
  //   {
  //     //then we are dealing with MFT tracks do somthing
  //     printf(">>>>>>>>>>>> I have an MFT track BUT SHOULD NOT BE %d\n", 3);
  //   }
  //
  //
  // }
  // PROCESS_SWITCH(nchMft, processFwd, "look at fwd tracks", true);




};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nchMft>(cfgc)};
}
