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
#include "CommonConstants/LHCConstants.h"
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
  {//For bctoft0c
  DECLARE_SOA_ARRAY_INDEX_COLUMN(FT0, ft0s);//has_ft0s works now, without it doesn't
  }
  namespace ambii
  {//for MA2T
  DECLARE_SOA_INDEX_COLUMN(MFTTrack, track);
  DECLARE_SOA_INDEX_COLUMN(AmbiguousMFTTrack, ambMFTtrack);
  } // namespace ambii
  DECLARE_SOA_TABLE(MatchedToFT0, "AOD", "MAFT", indices::BCId, indices::FT0Ids);
  DECLARE_SOA_INDEX_TABLE_USER(MA2T, MFTTracks, "MA2T", ambii::MFTTrackId, ambii::AmbiguousMFTTrackId);
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

struct ProduceMFTAmbii {
  // build the index table MA2T
  Builds<aod::MA2T> idx;
  void init(InitContext const&){};
  //Additionnal methods provided:
  //mfttrack.has_ambMFTtrack()
  //mfttrack.ambMFTtrack()
};

using ExtBCs = soa::Join<aod::BCs, aod::Timestamps, aod::MatchedToFT0>;
using MFTTracksExtra = soa::Join<aod::MFTTracks,aod::MA2T>;//MFT track + index of ambiguous track

template <typename T>
T getCompatibleBCs(aod::AmbiguousMFTTrack const& atrack, aod::Collision const& collOrig, T const& bcs, int deltaBC)
{

  auto compBCs = atrack.bc_as<T>();//BC + info on FT0
  auto bcIter = compBCs.begin();//first element of compBC
  uint64_t firstBC = bcIter.globalBC();

  bcIter.moveToEnd();//does it move to the end or the next one after the end ?
  --bcIter;//to avoid a seg fault
  uint64_t lastBC = bcIter.globalBC();//gives the last COMPATIBLE BC in compBCs


  auto bcIt = collOrig.bc_as<T>();


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
    if (bcIt == bcs.end())
    {
      --bcIt;
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();
    }
  }
  else
  {
    //here bcIt.globalBC() >= firstBC + deltaBC

    while (bcIt != bcs.begin() && bcIt.globalBC() > firstBC + deltaBC)
    {
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();

      --bcIt;
    }
  }

  int64_t maxBCId = bcIt.globalIndex();
  auto maxGlobalBC = bcIt.globalBC();

  while (bcIt != bcs.end() && bcIt.globalBC() < lastBC + deltaBC)
  {
    maxBCId = bcIt.globalIndex();
    maxGlobalBC = bcIt.globalBC();

    ++bcIt;
  }

  if (bcIt != bcs.end() && maxBCId >= minBCId)
  {
    T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};
    bcs.copyIndexBindings(slice);
    return slice;
  }
  else
  {
    T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId)}, (uint64_t)minBCId};
    bcs.copyIndexBindings(slice);
    return slice;
  }

}

template <typename T>
T getCompatibleBCs(MFTTracksExtra::iterator const& track, aod::Collision const& collOrig, T const& bcs, int deltaBC)
{

  //define firstBC and lastBC (globalBC of beginning and end of the range, when no shift is applied)

  if (track.has_ambMFTtrack())
  {
    auto atrack = track.ambMFTtrack();
    return getCompatibleBCs(atrack, collOrig, bcs, deltaBC);
  }

  auto bcIt = collOrig.bc_as<T>();
  //auto timstp = bcIt.timestamp();

  int64_t firstBC = bcIt.globalBC() + (track.trackTime() - track.trackTimeRes())/o2::constants::lhc::LHCBunchSpacingNS;
  int64_t lastBC = firstBC + 2*track.trackTimeRes()/o2::constants::lhc::LHCBunchSpacingNS +1;//to have a delta = 198 BC

  //int collTimeResInBC = collOrig.collisionTimeRes()/o2::constants::lhc::LHCBunchSpacingNS;

  int64_t collFirstBC = bcIt.globalBC() + (collOrig.collisionTime() - collOrig.collisionTimeRes())/o2::constants::lhc::LHCBunchSpacingNS;
  int64_t collLastBC = collFirstBC + 2*collOrig.collisionTimeRes()/o2::constants::lhc::LHCBunchSpacingNS +1;

  if ((bcIt.globalBC() < firstBC)||(bcIt.globalBC() > lastBC))//bcOfColl not in the range
  {//must take into account the time resolution of the coll
    if ((bcIt.globalBC() < firstBC)&&(collLastBC < firstBC))
    {
      //should not happen
      if (firstBC-collLastBC>=73)//extra margins of 20BC and 53BC
      {
        printf("------ Range [%lld,%lld]\n", firstBC, lastBC);
        printf("---collBC out of range (left case), difference %lld\n", firstBC-collLastBC);
        printf("---coll Range [%lld,%lld]\n", collFirstBC, collLastBC);
      }

    }
    if ((bcIt.globalBC() > lastBC)&&(collFirstBC > lastBC))
    {
      //should not happen
      if (collFirstBC-lastBC>=73)
      {
        printf("------ Range [%lld,%lld]\n", firstBC, lastBC);

        printf("---collBC out of range (right case), difference %lld\n", collFirstBC-lastBC);
        printf("---coll Range [%lld,%lld]\n", collFirstBC, collLastBC);
      }

    }

    //both of these cases happen, quite a lot actually
  }

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
    if (bcIt == bcs.end())
    {
      --bcIt;
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();
    }
  }
  else
  {
    //here bcIt.globalBC() >= firstBC + deltaBC

    while (bcIt != bcs.begin() && bcIt.globalBC() > firstBC + deltaBC)
    {
      minBCId = bcIt.globalIndex();
      minGlobalBC = bcIt.globalBC();

      --bcIt;
    }
  }

  int64_t maxBCId = bcIt.globalIndex();
  auto maxGlobalBC = bcIt.globalBC();

  while (bcIt != bcs.end() && bcIt.globalBC() < lastBC + deltaBC)
  {
    maxBCId = bcIt.globalIndex();
    maxGlobalBC = bcIt.globalBC();

    ++bcIt;
  }


  if (bcIt != bcs.end() && maxBCId >= minBCId)
  {
    T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId + 1)}, (uint64_t)minBCId};
    bcs.copyIndexBindings(slice);
    return slice;
  }
  else
  {
    T slice{{bcs.asArrowTable()->Slice(minBCId, maxBCId - minBCId)}, (uint64_t)minBCId};
    bcs.copyIndexBindings(slice);
    return slice;
  }

}

using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct matchmftfit {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<bool> verbose{"verbose", false, "will print more things if true"};
  Configurable<bool> isMC{"isMC", false, "set to true if dataset is MC"};


  int runNumber = -1;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  int count = 0;
  o2::parameters::GRPMagField* grpmag = nullptr;

  Configurable<int> shiftBC{"shiftBC", 0, "shift in BC wrt normal"};

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
      {"NCompBCwFT0C", "; N_{compBC}; #count", {HistType::kTH1D, {{21, -0.5, 20.5}}}},
      {"NCompBCwFT0s", "; N_{compBC}; #count", {HistType::kTH1D, {{21, -0.5, 20.5}}}},
      {"possibleIsTrue", "; possible = true; #count", {HistType::kTH1D, {{2, -0.5, 1.5}}}},                            //
      {"GoodIsTrue", "; good = true; #count", {HistType::kTH1D, {{11, 0, 10}}}}                            //
    }                                                                                                        //
  };

  void init(o2::framework::InitContext& initContext)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // if (!o2::base::GeometryManager::isGeometryLoaded()) {
    //   ccdb->get<TGeoManager>(geoPath);
    // }
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


  void processAmbi(MFTTracksExtra const& mfttracks,
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

    int i=0;//counts the number of channels having non-zero amplitude
    //for a particular BC
    double D = 0.0; //distance between (xe,ye,ze) and (xc,yc,zc)
    double minD;
    double globalMinD;



    for (auto& track : mfttracks)
    {
      globalMinD=999.;//minimum D for all BC
      ExtBCs::iterator possibleBC;//compatible BC with the D the smallest
      //beware: there could be several BC with the same smallest D


      if (!track.has_collision())
      {
        continue;
      }
      auto collOrig = track.collision();

      bool isAmbiguous = false;



      if (track.has_ambMFTtrack())
      {
        isAmbiguous=true;
      }


      auto bcSlice = getCompatibleBCs(track, collOrig, bcs, shiftBC);

      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());

      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};


      //we propagate the MFT track to the mean z position of FT0-C
      //trackPar.propagateToZlinear(-82.6);//in cm
      //getTrackPar() doesn't work because mft tracks don't have alpha
      trackPar.propagateToZhelix(-82.6, Bz);


      if ((abs(trackPar.getX())<6.365) && (abs(trackPar.getY())<6.555))
      {
        //track outside the FT0-C acceptance (in the central hole)
        continue;
      }


      std::vector<ExtBCs::iterator> goodBC;
      int nCompBCwft0C = 0;
      int nCompBCwft0s = 0;//Number of compatible BCs with FT0A AND C signals

      bool hasft0A = false;
      bool hasft0C = false;
      for (auto& bc : bcSlice)
      {
        hasft0C = false;
        hasft0A = false;
        //printf("----------bcId %lld\n", bc.globalIndex());
          if (!bc.has_ft0s())
          {
            continue;
          }

          auto ft0s = bc.ft0s();
          i=0;//reinitialise
          D=0.0;
          minD=999.9;
          for (auto const& ft0 : ft0s)
          {
            //printf("---------ft0.bcId %d\n", ft0.bcId());
            if (ft0.channelA().size() > 0)
            {
              //we are only considering BC with signals from both FT0A and FT0C
              hasft0A = true;
            }

            if (ft0.channelC().size() > 0)
            {
              hasft0C = true;
            }
            for (auto channelId : ft0.channelC())
            {

              std::vector<float> Xc = channelCoord[channelId];//(xc,yc,zc) coordinates
              //D in cm
              D=sqrt(pow(Xc[0]*0.1-trackPar.getX(),2)+pow(Xc[1]*0.1-trackPar.getY(),2)+pow(Xc[2]*0.1-1.87-trackPar.getZ(),2));
              //printf("----channelId %u, D %f, n %d\n", channelId, D, n);//should be between 96 and 207
              if (D < minD)
              {
                minD=D;
              }

              registry.fill(HIST("DistChannelToProp"), D);

            }


            i+=ft0.channelC().size();
          }

          registry.fill(HIST("NchannelsPerBC"), i);
          if (hasft0C)
          {
            nCompBCwft0C++;//number of compatible BCs that have ft0-C signal
          }
          //----------------------- BC selection here ------------------------
          //we are only considering BC with signals from both FT0A and FT0C
          if (!(hasft0A && hasft0C))
          {
            continue;
            //we go to the next BC
          }
          nCompBCwft0s++;
          //----------------------- end of BC selection ------------------------

          if (minD < 2)//20 mm
          {
            goodBC.push_back(bc);//goodBC is a vector of bc
          }
          if (minD < globalMinD)
          {
            globalMinD=minD;
            possibleBC=bc;
          }



      }
      registry.fill(HIST("NgoodBCperTrack"), goodBC.size());

      registry.fill(HIST("NCompBCwFT0C"), nCompBCwft0C);
      registry.fill(HIST("NCompBCwFT0s"), nCompBCwft0s);
    }
    n++;//counter for visualisation
  }
  PROCESS_SWITCH(matchmftfit, processAmbi, "Process only ambiguous tracks", true);


};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<bctoft0c>(cfgc),
                        adaptAnalysisTask<ProduceMFTAmbii>(cfgc),
                        adaptAnalysisTask<matchmftfit>(cfgc)
                      };
  return workflow;
}
