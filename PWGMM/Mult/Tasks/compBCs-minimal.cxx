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

// \file   match-mft-ft0.cxx
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief This code loops over every ambiguous MFT tracks and propagate
//        them to the FT0-C, to reduce track ambiguity

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "MathUtils/Utils.h"
#include "CommonConstants/LHCConstants.h"
#include "Common/Core/trackUtilities.h" //for getTrackPar()
#include "ReconstructionDataFormats/TrackFwd.h" //for propagate
//https://github.com/AliceO2Group/AliceO2/blob/dev/DataFormats/Reconstruction/include/ReconstructionDataFormats/TrackFwd.h
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"

#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"


using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using namespace o2;
using namespace o2::framework;

int n = 0;

namespace o2::aod
{
  namespace indices
  {
  DECLARE_SOA_ARRAY_INDEX_COLUMN(FT0, ft0s);//has_ft0s works now, without it doesn't
  }
  DECLARE_SOA_TABLE(MatchedToFT0, "AOD", "MAFT",
                    indices::BCId, indices::FT0Ids);
} // namespace o2::aod


//Creating a table BC to FT0 and filling it
struct bctoft0c {
  Produces<aod::MatchedToFT0> mf;
  struct {
    std::vector<int> ft0ids;
  } filler;
  void process(aod::BCs::iterator const& bc, soa::SmallGroups<aod::FT0s> const& ft0s)
  {
    filler.ft0ids.clear();
    for (auto const& ft0 : ft0s) {
      filler.ft0ids.emplace_back(ft0.globalIndex());
    }
    mf(bc.globalIndex(), filler.ft0ids);
  }
};

using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedToFT0>;

template <typename T>
T getCompatibleBCs(aod::AmbiguousMFTTrack const& atrack, aod::Collision const& collOrig, T const& bcs, int deltaBC)
{

  auto compBCs = atrack.bc_as<T>();//BC + info on FT0
  auto bcIter = compBCs.begin();//first element of compBC
  uint64_t firstBC = bcIter.globalBC();

  bcIter.moveToEnd();//does it move to the end or the next one after the end ?
  --bcIter;//to avoid a seg fault
  uint64_t lastBC = bcIter.globalBC();//gives the last COMPATIBLE BC in compBCs



  //printf("-------- BC range before : %llu - %llu\n", firstBC, lastBC);

  auto bcIt = collOrig.bc_as<T>();


  //printf("-------- BC range goal: %llu - %llu\n", firstBC + deltaBC, lastBC + deltaBC);


  //printf("------- we start from globalBC = %lld\n", bcIt.globalBC());




  int64_t minBCId = bcIt.globalIndex();
  auto minGlobalBC = bcIt.globalBC();


  if (bcIt.globalBC() < firstBC + deltaBC)
  {
    while (bcIt != bcs.end() && bcIt.globalBC() < firstBC + deltaBC)
    {
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();

      ++bcIt;
    }

    //printf("------ the condition failed for bcIt.globalBC() = %lld, bcIt.globalIndex() %lld\n", bcIt.globalBC(), bcIt.globalIndex());

    //----Just a check :
    //--bcIt;

    //printf("------ the previous was bcIt.globalBC() = %lld, bcIt.globalIndex() %lld\n", bcIt.globalBC(), bcIt.globalIndex());

  }
  else
  {
    //here bcIt.globalBC() >= firstBC + deltaBC

    while (bcIt != bcs.end() && bcIt.globalBC() > firstBC + deltaBC)
    {
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();

      --bcIt;
    }

    //printf(">>>>>> the condition failed for bcIt.globalBC() = %lld, bcIt.globalIndex() %lld\n", bcIt.globalBC(), bcIt.globalIndex());


  }

  int64_t maxBCId = bcIt.globalIndex();
  auto maxGlobalBC = bcIt.globalBC();

  while (bcIt != bcs.end() && bcIt.globalBC() < lastBC + deltaBC)
  {
    maxBCId = bcIt.globalIndex();
    maxGlobalBC = bcIt.globalBC();

    ++bcIt;
  }

  //printf("_________Will consider BC entries from %lld to %lld\n", minGlobalBC, maxGlobalBC);
  //printf("_________Means %lld BC\n", maxGlobalBC-minGlobalBC+1);


  T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};//change back the 3
  bcs.copyIndexBindings(slice);
  return slice;
}


using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct compBCsminimal {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<bool> verbose{"verbose", false, "will print more things if true"};
  Configurable<bool> isMC{"isMC", false, "set to true if dataset is MC"};


  int runNumber = -1;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  int count = 0;
  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  std::vector<std::vector<float>> channelCoord = {{103.2,17.8,-813.1},{76.9,17.8,-815.9},{103.1,44.2,-812.1},{76.8,44.2,-814.9},{103.2,78.7,-810},{76.8,79,-812.9},{103.2,105,-807.1},{76.8,105.3,-810},{43.2,78.8,-815},{43.2,105.1,-812.1},{16.8,78.9,-815.9},{16.8,105.2,-813},{-16.8,105.2,-813},{-16.8,78.9,-815.9},{-43.2,105.1,-812.1},{-43.2,78.8,-815},{-76.8,105.3,-810},{-76.8,79,-812.9},{-103.2,105,-807.1},{-103.2,78.7,-810},{-76.8,44.2,-814.9},{-103.1,44.2,-812.1},{-76.9,17.8,-815.9},{-103.2,17.8,-813.1},{-103.2,-17.8,-813.1},{-76.9,-17.8,-815.9},{-103.1,-44.2,-812.1},{-76.8,-44.2,-814.9},{-103.2,-78.7,-810},{-76.8,-79,-812.9},{-103.2,-105,-807.1},{-76.8,-105.3,-810},{-43.2,-78.8,-815},{-43.2,-105.1,-812.1},{-16.8,-78.9,-815.9},{-16.8,-105.2,-813},{16.8,-105.2,-813},{16.8,-78.9,-815.9},{43.2,-105.1,-812.1},{43.2,-78.8,-815},{76.8,-105.3,-810},{76.8,-79,-812.9},{103.2,-105,-807.1},{103.2,-78.7,-810},{76.8,-44.2,-814.9},{103.1,-44.2,-812.1},{76.9,-17.8,-815.9},{103.2,-17.8,-813.1},{163,18.7,-804.1},{137,18.9,-808.9},{163,45.2,-803.1},{137,45.3,-807.9},{163,78.6,-800.1},{137,79.1,-804.9},{163,104.9,-797.2},{137,105.4,-801.9},{103.4,138,-802},{102.9,164,-797.2},{77.1,138,-804.9},{76.6,164,-800},{43.3,139,-807},{43.2,165,-802.1},{16.9,139,-807.9},{16.7,165,-803},{-16.7,165,-803},{-16.9,139,-807.9},{-43.2,165,-802.1},{-43.3,139,-807},{-76.6,164,-800},{-77.1,138,-804.9},{-102.9,164,-797.2},{-103.4,138,-802},{-137,105.4,-801.9},{-163,104.9,-797.2},{-137,79.1,-804.9},{-163,78.6,-800.1},{-137,45.3,-807.9},{-163,45.2,-803.1},{-137,18.9,-808.9},{-163,18.7,-804.1},{-163,-18.7,-804.1},{-137,-18.9,-808.9},{-163,-45.2,-803.1},{-137,-45.3,-807.9},{-163,-78.6,-800.1},{-137,-79.1,-804.9},{-163,-104.9,-797.2},{-137,-105.4,-801.9},{-103.4,-138,-802},{-102.9,-164,-797.2},{-77.1,-138,-804.9},{-76.6,-164,-800},{-43.3,-139,-807},{-43.2,-165,-802.1},{-16.9,-139,-807.9},{-16.7,-165,-803},{16.7,-165,-803},{16.9,-139,-807.9},{43.2,-165,-802.1},{43.3,-139,-807},{76.6,-164,-800},{77.1,-138,-804.9},{102.9,-164,-797.2},{103.4,-138,-802},{137,-105.4,-801.9},{163,-104.9,-797.2},{137,-79.1,-804.9},{163,-78.6,-800.1},{137,-45.3,-807.9},{163,-45.2,-803.1},{137,-18.9,-808.9},{163,-18.7,-804.1}};

  HistogramRegistry registry{
    "registry",
    {
      {"ChannelFired", "; #it{x} (cm); #it{y} (cm);", {HistType::kTH2I, {{701, -35.05, 35.05}, {701, -35.05, 35.05}}}}, //
      {"TrackPosProp", "; #it{x} (cm); #it{y} (cm);", {HistType::kTH2F, {{701, -35.05, 35.05}, {701, -35.05, 35.05}}}}, //
      {"DistChannelToProp", "; D (cm); #count", {HistType::kTH1D, {{101, 0, 100}}}},
      {"NchannelsPerBC", "; N_{channelC}; #count", {HistType::kTH1D, {{101, 0, 100}}}},
      {"NgoodBCperTrack", "; N_{goodBC}; #count", {HistType::kTH1D, {{11, 0, 10}}}},                            //
      {"NCompBCwFT0", "; N_{compBC}; #count", {HistType::kTH1D, {{21, -0.5, 20.5}}}},
      {"possibleIsTrue", "; possible = true; #count", {HistType::kTH1D, {{2, -0.5, 1.5}}}},                            //
      {"GoodIsTrue", "; good = true; #count", {HistType::kTH1D, {{11, 0, 10}}}}                            //
    }                                                                                                        //
  };

  void init(o2::framework::InitContext& initContext)
  {
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


  void process(aod::MFTTracks const&,
               aod::Collisions const&, ExtBCs const& bcs,
               aod::FT0s const&,
               aod::AmbiguousMFTTracks const& atracks)
  {
    if (bcs.size() == 0) {
      return;
    }
    if (atracks.size() == 0) {
      LOG(info) << "The ambiguous track table is empty";
      return;
    }
    initCCDB(bcs.begin());


    printf("-------- %d\n", 1);
    for (auto& atrack : atracks)
    {
      printf("-------- %d\n", 2);
      if (!atrack.has_bc())
      {
        continue;
      }
      auto compatibleBCs = atrack.bc_as<ExtBCs>();//BC + info on FT0
      printf("-------- %f\n", 2.5);
      if (!atrack.has_mfttrack())
      {
        continue;
      }
      printf("-------- %f\n", 2.6);
      auto track = atrack.mfttrack();
      printf("-------- %d\n", 3);



      if (!track.has_collision())
      {
        continue;
      }
      auto collOrig = track.collision();
      printf("-------- %d\n", 4);
      auto bcSlice = getCompatibleBCs(atrack, collOrig, bcs,50);

    }
    n++;//counter for visualisation
  }


};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<bctoft0c>(cfgc),
                        adaptAnalysisTask<compBCsminimal>(cfgc)};
  return workflow;
}
