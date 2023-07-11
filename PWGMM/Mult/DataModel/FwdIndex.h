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

// Table to reverse matching indices from tracks to fwtracks,
// can be Joined to the MFTTrack Table
// Gives access to new methods mfttrack.has_fwdtrack() and mfttrack.fwdtrack()
#ifndef O2_ANALYSIS_FWDINDEX_H_
#define O2_ANALYSIS_FWDINDEX_H_

#include "Framework/AnalysisDataModel.h"
namespace o2::aod
{
namespace fwdidx
{
DECLARE_SOA_INDEX_COLUMN(FwdTrack, fwdtrack);
} // namespace fwdidx
DECLARE_SOA_TABLE(MFTTracksToFwdTracks, "AOD", "MFT2FWDT", fwdidx::FwdTrackId);
} // namespace o2::aod
#endif // O2_ANALYSIS_FWDINDEX_H_
