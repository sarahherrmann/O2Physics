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
#include "MatchMFTFT0.h"

#include "MathUtils/Utils.h"
#include "CommonConstants/LHCConstants.h"

#include "Common/DataModel/CollisionAssociationTables.h"

using namespace o2;
using namespace o2::framework;

using CollwEv = soa::Join<aod::Collisions, aod::EvSels>;

AxisSpec PhiAxis = {600, -2*M_PI, 2*M_PI};
AxisSpec EtaAxis = {180, -8.6, 8.6};

struct nchMft {

  Service<o2::framework::O2DatabasePDG> pdg;

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
    if (doprocessMuonsAndMFTFT0)
    {
      registry.add({"MFTFT0/SameBCFwd", "; sameBC; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFTFT0/BCInMFTCollBC", "; isBCInMFTCollBC; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFTFT0/bcDiff", "; bcDiff; count", {HistType::kTH1F, {{601, -300.5, 300.5}}}});
      registry.add({"MFTFT0/bcDiffMFT", "; bcDiffMFT; count", {HistType::kTH1F, {{601, -300.5, 300.5}}}});
      registry.add({"MFTFT0/Chi2Matching", "; chi2; count", {HistType::kTH1F, {{301, -0.5, 300.5}}}});
      registry.add({"MFTFT0/Chi2MatchingBCmatch", "; chi2; count", {HistType::kTH1F, {{301, -0.5, 300.5}}}});
      registry.add({"MFTFT0/Chi2MatchingNoBCmatch", "; chi2; count", {HistType::kTH1F, {{301, -0.5, 300.5}}}});
      registry.add({"MFTFT0/TracksDeltaPhiDeltaEta", "; #Delta#varphi; #Delta#eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
    }

    if (doprocessMFT)
    {
      registry.add({"MFT/BCInMFTCollBC", "; isBCInMFTCollBC; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFT/BCInMFTCollBCAll", "; isBCInMFTCollBC; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFT/BChasColl", "; atleast1BChasColl; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFT/BChasCollAll", "; hasColl; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFT/BChasCollButNotCompCollAll", "; hasColl; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFT/BChasCollButNotCompColl", "; hasColl; count", {HistType::kTH1F, {{2, -0.5, 1.5}}}});
      registry.add({"MFT/bcDiffMFT", "; bcDiffMFT; count", {HistType::kTH1F, {{601, -300.5, 300.5}}}});
      registry.add({"MFT/bcDiffMiddleROF", "; bcMiddleROF - bc(FT0-C); count", {HistType::kTH1F, {{601, -300.5, 300.5}}}});
      registry.add({"MFT/bcDiffCompCollMiddleROF", "; bcCompColl - bcMiddleROF; count", {HistType::kTH1F, {{2001, -1000.5, 1000.5}}}});
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
  PROCESS_SWITCH(nchMft, processMuonsAndMFT, "Chi2 of muons and MFt track matching 2", false);

  using MFTwFwdFT0 = soa::Join<aod::MFTTracks, o2::aod::MFTTracksToFwdTracks, aod::BCofMFT, aod::MFTTrkCompColls>;

  bool isMuonSelected(FwdTracksWColls::iterator const& muon)
  {
    if ((muon.trackType() != 0)||(muon.compatibleCollIds().size()>1))
    {
      return false;
    }//we just want to look at the global muon tracks, to have the least amount of ambiguous fwd tracks

    if (!muon.has_collision())
    {
      return false;
    }

    if (muon.eta()< -4 || muon.eta()> -2.5)
    {
      return false;
    }

    if (muon.rAtAbsorberEnd()<17.6 || muon.rAtAbsorberEnd()>89.5)
    {
      //must be between 17.6 and 89.5
      return false;
    }

    if (muon.rAtAbsorberEnd()<26.5 || muon.rAtAbsorberEnd()>89.5)
    {
      //then some pDCA cuts
      if (muon.pDca() > 324.0)
      {
        return false;
      }

    }

    if (muon.rAtAbsorberEnd()<17.6 || muon.rAtAbsorberEnd()>26.5)
    {
      //then some pDCA cuts
      if (muon.pDca() > 594.0)
      {
        return false;
      }

    }

    return true;
  }

  void processMuonsAndMFTFT0(MFTwFwdFT0 const& mfttracks, CollwEv const& collisions, FwdTracksWColls const&, aod::BCs const&)
  {
    for (auto& mfttrack : mfttracks)
    {

      if (!(mfttrack.has_fwdtrack() && mfttrack.has_bcs()))//mft tracks having a matched fwd track and a match in FT0-C
      {
        continue;
      }
      auto matchedFwdtrack = mfttrack.fwdtrack_as<FwdTracksWColls>();
      if (!isMuonSelected(matchedFwdtrack))
      {
        continue;
      }
      int64_t collFwdBC = matchedFwdtrack.collision_as<CollwEv>().bc().globalBC();// the collision's bc of the global muon track
      int sameBC = 0;
      int isBCInMFTCollBC = 0;
      auto compatibleColls = mfttrack.compatibleColl_as<CollwEv>();

      if (!matchedFwdtrack.collision_as<CollwEv>().sel8()){
        continue;
      }

      //auto compatibleColls = mfttrack.compatibleColl_as<CollwEv>();
      //printf("-- mfttrack index %lld, mfttrack.bcs().size() %lu\n", mfttrack.globalIndex(), mfttrack.bcs().size());
      int64_t bcDiffMin = 100000000;
      int64_t bcDiffMinMFT = 100000000;
      for (auto& bc : mfttrack.bcs())
      {//
        int64_t bcDiff = bc.globalBC() - collFwdBC;

        // if(mfttrack.globalIndex() == 86404)
        // {
        //   printf("---------- bc.globalBC() %lld, collFwdBC %lld\n", bc.globalBC(), collFwdBC);
        //   printf("bcDiff %lld\n", bcDiff);
        // }
        if (abs(bcDiffMin)>abs(bcDiff))
        {
          bcDiffMin=bcDiff;
        }

        if (abs(bcDiff) < 5)
        {//we consider that we have a precision of 5 BC
          sameBC=1;
        }



        for (auto& collMFT : compatibleColls)
        {//looking at all the collisions compatible with the MFT track
          int64_t bcDiffMFT = bc.globalBC() - collMFT.bc().globalBC();
          if (!collMFT.sel8())
          {
            continue;
          }
          if (collMFT.bc().globalBC() == bc.globalBC())
          {//is the bc found by FT0-C is in the bcs of the comp colls ?
            isBCInMFTCollBC=1;
          }

          if (abs(bcDiffMinMFT)>abs(bcDiffMFT))
          {
            bcDiffMinMFT=bcDiffMFT;
          }

        }
      }

      registry.fill(HIST("MFTFT0/bcDiff"), bcDiffMin);
      registry.fill(HIST("MFTFT0/bcDiffMFT"), bcDiffMinMFT);
      registry.fill(HIST("MFTFT0/SameBCFwd"), sameBC);
      registry.fill(HIST("MFTFT0/BCInMFTCollBC"), isBCInMFTCollBC);
      registry.fill(HIST("MFTFT0/Chi2Matching"), matchedFwdtrack.chi2MatchMCHMFT());
      if (sameBC)
      {
        registry.fill(HIST("MFTFT0/Chi2MatchingBCmatch"), matchedFwdtrack.chi2MatchMCHMFT());
      }
      else
      {
        registry.fill(HIST("MFTFT0/Chi2MatchingNoBCmatch"), matchedFwdtrack.chi2MatchMCHMFT());
      }
      float phiMFT = mfttrack.phi();
      o2::math_utils::bringTo02Pi(phiMFT);

      float phiMuon = matchedFwdtrack.phi();
      o2::math_utils::bringTo02Pi(phiMuon);

      registry.fill(HIST("MFTFT0/TracksDeltaPhiDeltaEta"), phiMuon-phiMFT, matchedFwdtrack.eta()-mfttrack.eta());

      //registry.fill(HIST("IsInMFTColls"), isInMFTColls);

    }
  }
  PROCESS_SWITCH(nchMft, processMuonsAndMFTFT0, "Chi2 of muons and MFt track matching 2", false);

  using BCwColls = soa::Join<aod::BCs, aod::MatchedBCCollisionsSparseMulti>;

  void processMFT(MFTwFwdFT0 const& mfttracks, CollwEv const& collisions, BCwColls const&)
  {
    for (auto& mfttrack : mfttracks)
    {

      if (!(mfttrack.has_bcs() && mfttrack.has_compatibleColl() && mfttrack.has_collision()))//mft tracks having a match in FT0-C
      {
        continue;
      }

      int64_t middleBC = 47+mfttrack.collision_as<CollwEv>().bc_as<BCwColls>().globalBC() + mfttrack.trackTime()/o2::constants::lhc::LHCBunchSpacingNS;

      auto compatibleColls = mfttrack.compatibleColl_as<CollwEv>();

      int64_t bcDiffMinMFT = 100000000;
      int isBCInMFTCollBC=0;//at least one BC
      int isBCInMFTCollBCAll=0;//for each BC
      int BChasColl=0;//at least one BC
      int BChasCollAll=0;//for each BC

      for (auto& bc : mfttrack.bcs_as<BCwColls>())
      {
        BChasCollAll=0;
        isBCInMFTCollBCAll=0;

        int64_t diffBC = bc.globalBC() - middleBC; //should not be bigger than 99 BC
        registry.fill(HIST("MFT/bcDiffMiddleROF"), diffBC);

        if (bc.has_collisions())
        {
          BChasColl=1;
          BChasCollAll=1;
        }

        for (auto& collMFT : compatibleColls)
        {//looking at all the collisions compatible with the MFT track
          int64_t bcDiffMFT = bc.globalBC() - collMFT.bc_as<BCwColls>().globalBC();
          if (!collMFT.sel8())
          {
            continue;
          }
          if (collMFT.bc_as<BCwColls>().globalBC() == bc.globalBC())
          {//is the bc found by FT0-C is in the bcs of the comp colls ?
            isBCInMFTCollBC=1;
            isBCInMFTCollBCAll=1;
          }

          if (abs(bcDiffMinMFT)>abs(bcDiffMFT))
          {
            bcDiffMinMFT=bcDiffMFT;
          }

          int64_t bcDiffCompCollMiddleROF = collMFT.bc_as<BCwColls>().globalBC() - middleBC;
          registry.fill(HIST("MFT/bcDiffCompCollMiddleROF"), bcDiffCompCollMiddleROF);

        }//for (auto& collMFT : compatibleColls)


        registry.fill(HIST("MFT/BCInMFTCollBCAll"), isBCInMFTCollBCAll, 1.0/mfttrack.bcs_as<BCwColls>().size());
        registry.fill(HIST("MFT/BChasCollAll"), BChasCollAll, 1.0/mfttrack.bcs_as<BCwColls>().size());

        if (!isBCInMFTCollBCAll)
        {
          registry.fill(HIST("MFT/BChasCollButNotCompCollAll"), BChasCollAll, 1.0/mfttrack.bcs_as<BCwColls>().size());
        }

      }//for (auto& bc : mfttrack.bcs_as<BCwColls>())
      registry.fill(HIST("MFT/BCInMFTCollBC"), isBCInMFTCollBC);
      registry.fill(HIST("MFT/BChasColl"), BChasColl);
      registry.fill(HIST("MFT/bcDiffMFT"), bcDiffMinMFT);

      if (!isBCInMFTCollBC)
      {
        registry.fill(HIST("MFT/BChasCollButNotCompColl"), BChasColl);
      }

    }
  }
  PROCESS_SWITCH(nchMft, processMFT, "processMFT", true);

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
  PROCESS_SWITCH(nchMft, processMuonsAndMFTReas, "Chi2 of muons and MFt track matching", false);



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
  PROCESS_SWITCH(nchMft, processFwdAmb, "find the percentage of fwd ambiguous tracks", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<nchMft>(cfgc)};
}
