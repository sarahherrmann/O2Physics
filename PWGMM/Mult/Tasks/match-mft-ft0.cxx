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
using MFTTracksLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

struct matchmftfit {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<bool> verbose{"verbose", false, "will print more things if true"};


  int runNumber = -1;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
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

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {
    LOG(info) << "INITIALIZING CCDB";

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


  void process(MFTTracksLabeled const&,
               aod::Collisions const&, ExtBCs const& bcs,
               aod::FT0s const&,
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
    initCCDB(bcs.begin());

    int i=0;//counts the number of channels having non-zero amplitude
    //for a particular BC
    double D = 0.0; //distance between (xe,ye,ze) and (xc,yc,zc)
    double minD;
    double globalMinD;
    bool isSet=false;
    for (auto& atrack : atracks)
    {
      globalMinD=999.;//minimum D for all BC
      isSet=false;
      ExtBCs::iterator possibleBC;//compatible BC with the D the smallest
      //beware: there could be several BC with the same smallest D

      auto compatibleBCs = atrack.bc_as<ExtBCs>();//BC + info on FT0
      //auto track = atrack.mfttrack();
      //printf("------------------------------mfttrackId %d\n", atrack.mfttrackId());
      auto track = atrack.mfttrack_as<MFTTracksLabeled>();
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for ambiguous track, skip...");
        continue;
      }
      auto particle = track.mcParticle();
      int bcIdTruth = particle.mcCollision().bcId();

      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      //we propagate the MFT track to the mean z position of FT0-C
      //trackPar.propagateToZlinear(-82.6);//in cm
      //getTrackPar() doesn't work because mft tracks don't have alpha
      trackPar.propagateToZhelix(-82.6, Bz);


      std::vector<ExtBCs::iterator> goodBC;
      int nCompBCwft0 = 0;
      bool hasft0 = false;
      for (auto& bc : compatibleBCs)
      {
        hasft0 = false;
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
            if (ft0.channelC().size()>=1)
            {
              hasft0=true;
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

              //if ((bc.globalIndex() ==4817) && (n==1) && (atrack.mfttrackId()==8982))//selection to draw visualisation
              if ((n==49) && (atrack.mfttrackId()==15776))//selection to draw visualisation
              {
                registry.fill(HIST("ChannelFired"), Xc[0]*0.1, Xc[1]*0.1, pow(10,nCompBCwft0));
                registry.fill(HIST("TrackPosProp"), trackPar.getX(), trackPar.getY());
                printf("----------bcId %lld\n", bc.globalIndex());
                if (bc.globalIndex()==bcIdTruth)
                {
                  printf("---- FOUND THE TRUE BC %d, nCompBCwft0 = %d\n", bcIdTruth, nCompBCwft0);
                }
                printf("----channelId %u, D %f, n %d\n", channelId, D, n);//should be between 96 and 207
                printf("------------------------------channel x %f y %f Track x %f, y%f\n", Xc[0]*0.1, Xc[1]*0.1, trackPar.getX(), trackPar.getY());
              }
            }


            i+=ft0.channelC().size();
          }
          if (minD < 1.9)//19 mm
          {
            goodBC.push_back(bc);//goodBC is a vector of bc
          }
          if (minD < globalMinD)
          {
            globalMinD=minD;
            possibleBC=bc;
            isSet=true;
          }

          registry.fill(HIST("NchannelsPerBC"), i);
          if (hasft0)
          {
            nCompBCwft0++;//number of compatible BCs that have ft0 signal
          }

      }
      registry.fill(HIST("NgoodBCperTrack"), goodBC.size());
      for (auto bci: goodBC){
        if (bci.globalIndex()==bcIdTruth)
        //the BC with signal closest to
        //track extrapolated is the true BC
        {
          registry.fill(HIST("GoodIsTrue"), goodBC.size(), 1);
        }
        else
        {
          registry.fill(HIST("GoodIsTrue"), goodBC.size(), 0);
        }
      }
      registry.fill(HIST("NCompBCwFT0"), nCompBCwft0);
      if((nCompBCwft0>1) && (verbose ==true))
      {
        printf("-------------------------------------------------- %d compatible BC with ft0 signal , mfttrackId %d, n %d\n", nCompBCwft0, atrack.mfttrackId(),n);
      }
      if (!isSet){continue;}

      if (possibleBC.globalIndex()==bcIdTruth)
      //the BC with signal closest to
      //track extrapolated is the true BC
      {
        registry.fill(HIST("possibleIsTrue"), 1);
      }
      else
      {
        registry.fill(HIST("possibleIsTrue"), 0);
      }
    }
    n++;//counter for visualisation
  }


};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<bctoft0c>(cfgc),
                        adaptAnalysisTask<matchmftfit>(cfgc)};
  return workflow;
}
