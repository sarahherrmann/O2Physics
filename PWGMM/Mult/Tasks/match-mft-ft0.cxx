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

using ExtBCs = soa::Join<aod::BCs, aod::MatchedToFT0>;


struct matchmftfit {

  std::vector<std::vector<float>> channelCoord = {{103.2,17.8,-813.1},{76.9,17.8,-815.9},{103.1,44.2,-812.1},{76.8,44.2,-814.9},{103.2,78.7,-810},{76.8,79,-812.9},{103.2,105,-807.1},{76.8,105.3,-810},{43.2,78.8,-815},{43.2,105.1,-812.1},{16.8,78.9,-815.9},{16.8,105.2,-813},{-16.8,105.2,-813},{-16.8,78.9,-815.9},{-43.2,105.1,-812.1},{-43.2,78.8,-815},{-76.8,105.3,-810},{-76.8,79,-812.9},{-103.2,105,-807.1},{-103.2,78.7,-810},{-76.8,44.2,-814.9},{-103.1,44.2,-812.1},{-76.9,17.8,-815.9},{-103.2,17.8,-813.1},{-103.2,-17.8,-813.1},{-76.9,-17.8,-815.9},{-103.1,-44.2,-812.1},{-76.8,-44.2,-814.9},{-103.2,-78.7,-810},{-76.8,-79,-812.9},{-103.2,-105,-807.1},{-76.8,-105.3,-810},{-43.2,-78.8,-815},{-43.2,-105.1,-812.1},{-16.8,-78.9,-815.9},{-16.8,-105.2,-813},{16.8,-105.2,-813},{16.8,-78.9,-815.9},{43.2,-105.1,-812.1},{43.2,-78.8,-815},{76.8,-105.3,-810},{76.8,-79,-812.9},{103.2,-105,-807.1},{103.2,-78.7,-810},{76.8,-44.2,-814.9},{103.1,-44.2,-812.1},{76.9,-17.8,-815.9},{103.2,-17.8,-813.1},{163,18.7,-804.1},{137,18.9,-808.9},{163,45.2,-803.1},{137,45.3,-807.9},{163,78.6,-800.1},{137,79.1,-804.9},{163,104.9,-797.2},{137,105.4,-801.9},{103.4,138,-802},{102.9,164,-797.2},{77.1,138,-804.9},{76.6,164,-800},{43.3,139,-807},{43.2,165,-802.1},{16.9,139,-807.9},{16.7,165,-803},{-16.7,165,-803},{-16.9,139,-807.9},{-43.2,165,-802.1},{-43.3,139,-807},{-76.6,164,-800},{-77.1,138,-804.9},{-102.9,164,-797.2},{-103.4,138,-802},{-137,105.4,-801.9},{-163,104.9,-797.2},{-137,79.1,-804.9},{-163,78.6,-800.1},{-137,45.3,-807.9},{-163,45.2,-803.1},{-137,18.9,-808.9},{-163,18.7,-804.1},{-163,-18.7,-804.1},{-137,-18.9,-808.9},{-163,-45.2,-803.1},{-137,-45.3,-807.9},{-163,-78.6,-800.1},{-137,-79.1,-804.9},{-163,-104.9,-797.2},{-137,-105.4,-801.9},{-103.4,-138,-802},{-102.9,-164,-797.2},{-77.1,-138,-804.9},{-76.6,-164,-800},{-43.3,-139,-807},{-43.2,-165,-802.1},{-16.9,-139,-807.9},{-16.7,-165,-803},{16.7,-165,-803},{16.9,-139,-807.9},{43.2,-165,-802.1},{43.3,-139,-807},{76.6,-164,-800},{77.1,-138,-804.9},{102.9,-164,-797.2},{103.4,-138,-802},{137,-105.4,-801.9},{163,-104.9,-797.2},{137,-79.1,-804.9},{163,-78.6,-800.1},{137,-45.3,-807.9},{163,-45.2,-803.1},{137,-18.9,-808.9},{163,-18.7,-804.1}};

  HistogramRegistry registry{
    "registry",
    {
      {"ChannelFired", "; #it{x} (cm); #it{y} (cm);", {HistType::kTH2F, {{701, -35.05, 35.05}, {701, -35.05, 35.05}}}}, //
      {"TrackPosProp", "; #it{x} (cm); #it{y} (cm);", {HistType::kTH2F, {{701, -35.05, 35.05}, {701, -35.05, 35.05}}}}, //
      {"DistChannelToProp", "; D (cm); #count", {HistType::kTH1D, {{101, 0, 100}}}},
      {"NchannelsPerBC", "; N_{channelC}; #count", {HistType::kTH1D, {{101, 0, 100}}}},
      {"NgoodBCperTrack", "; N_{goodBC}; #count", {HistType::kTH1D, {{11, 0, 10}}}}                            //
    }                                                                                                        //
  };


  void process(aod::MFTTracks const&,
               aod::Collisions const&, ExtBCs const& bcs,
               aod::FT0s const&,
               aod::AmbiguousMFTTracks const& atracks)
  {
    if (bcs.size() == 0) {
      return;
    }
    if (atracks.size() == 0) {
      return;
    }
    int i=0;//counts the number of channels having non-zero amplitude
    //for a particular BC
    double D = 0.0; //distance between (xe,ye,ze) and (xc,yc,zc)
    for (auto& atrack : atracks)
    {
      auto compatibleBCs = atrack.bc_as<ExtBCs>();//BC + info on FT0
      auto track = atrack.mfttrack();
      //printf("------------------------------mfttrackId %d\n", atrack.mfttrackId());

      std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
      SMatrix55 tcovs(v1.begin(), v1.end());
      SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
      o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};

      //we propagate the MFT track to the mean z position of FT0-C
      trackPar.propagateToZlinear(-80.73);//in cm
      //getTrackPar() doesn't work because mft tracks don't have alpha

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
              D=sqrt(pow(Xc[0]*0.1-trackPar.getX(),2)+pow(Xc[1]*0.1-trackPar.getY(),2)+pow(Xc[2]*0.1-trackPar.getZ(),2));
              //printf("----channelId %u, D %f, n %d\n", channelId, D, n);//should be between 96 and 207
              if (D < 1.5)//15 mm
              {
                goodBC.push_back(bc);//goodBC is a vector of bc
              }
              registry.fill(HIST("DistChannelToProp"), D);

              //if ((bc.globalIndex() ==4817) && (n==1) && (atrack.mfttrackId()==8982))//selection to draw visualisation
              if ((n==1) && (atrack.mfttrackId()==21186))//selection to draw visualisation
              {
                registry.fill(HIST("ChannelFired"), Xc[0]*0.1, Xc[1]*0.1, nCompBCwft0+1);
                registry.fill(HIST("TrackPosProp"), trackPar.getX(), trackPar.getY());
                printf("----------bcId %lld\n", bc.globalIndex());
                printf("----channelId %u, D %f, n %d\n", channelId, D, n);//should be between 96 and 207
                printf("------------------------------channel x %f y %f Track x %f, y%f\n", Xc[0]*0.1, Xc[1]*0.1, trackPar.getX(), trackPar.getY());
              }
            }

            i+=ft0.channelC().size();
          }
          registry.fill(HIST("NchannelsPerBC"), i);
          if (hasft0)
          {
            nCompBCwft0++;//number of compatible BCs that have ft0 signal
          }

      }
      registry.fill(HIST("NgoodBCperTrack"), goodBC.size());
      if(nCompBCwft0>1)
      {
        printf("-------------------------------------------------- %d compatible BC with ft0 signal , mfttrackId %d, n %d\n", nCompBCwft0, atrack.mfttrackId(),n);
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
