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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "MathUtils/Utils.h"

using namespace o2;
using namespace o2::framework;


struct AssessmentMFT {

  HistogramRegistry registry{
    "registry",
    {

      {"TracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {{600, 0, 2 * M_PI}, {35, -4.5, -1.}}}},         //
      {"TracksTime", "; time; #count", {HistType::kTH1F, {{600000, 0, 60000}}}},         //
    }                                                                                                                //
  };


  void process(aod::MFTTracks const& tracks)
  {
    for (auto& track : tracks)
    {
    float phi = track.phi();
    o2::math_utils::bringTo02Pi(phi);
    registry.fill(HIST("TracksPhiEta"), phi, track.eta());

    registry.fill(HIST("TracksTime"), track.trackTime());
    }
  }


};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<AssessmentMFT>(cfgc)};
}
