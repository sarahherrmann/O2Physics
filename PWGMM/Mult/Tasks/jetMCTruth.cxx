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

// \file   jetMCTruth.cxx
// \author Sarah Herrmann <sarah.herrmann@cern.ch>
//
// \brief This code loops over every MC particle and fills a table of their parton HF mother

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace jet_indices
{
DECLARE_SOA_ARRAY_INDEX_COLUMN(HFPartonMother, collisions);
}

DECLARE_SOA_TABLE(HFPartonM, "AOD", "MPARTON",
                  jet_indices::HFPartonMotherId);
} // namespace o2::aod

struct particleHFParton {
  Produces<aod::HFPartonM> hfm;
  struct {
    std::vector<int> ft0ids;
  } filler;
  void process(aod::BCs::iterator const& bc, soa::SmallGroups<aod::FT0s> const& ft0s)
  {
    filler.ft0ids.clear();
    for (auto const& ft0 : ft0s) {
      filler.ft0ids.emplace_back(ft0.globalIndex());
    }
    hfm(bc.globalIndex(), filler.ft0ids);
  }
};
