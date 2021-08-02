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

/// @file  TRDTrackReaderSpec.cxx

#include "TRDWorkflowIO/TRDTrackReaderSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/SerializationMethods.h"
#include "CommonUtils/StringUtils.h"
#include "ReconstructionDataFormats/MatchingType.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"

using GTrackID = o2::dataformats::GlobalTrackID;

namespace o2
{
namespace trd
{

void TRDTrackReader::init(InitContext& ic)
{

  mFileName = o2::utils::Str::concat_string(o2::utils::Str::rectifyDirectory(ic.options().get<std::string>("input-dir")),
                                            ic.options().get<std::string>("track-infile"));

  connectTree(mFileName);
}

void TRDTrackReader::run(ProcessingContext& pc)
{
  auto ent = mTree->GetReadEntry() + 1;
  assert(ent < mTree->GetEntries()); // this should not happen
  mTree->GetEntry(ent);
  LOG(INFO) << "Pushing " << mTracks.size() << " tracks and " << mTrigRec.size() << " trigger records at entry " << ent;

  if (mMode == Mode::TPCTRD) {
    uint32_t ss = o2::globaltracking::getSubSpec(mSubSpecStrict ? o2::globaltracking::MatchingType::Strict : o2::globaltracking::MatchingType::Standard);
    pc.outputs().snapshot(Output{o2::header::gDataOriginTRD, "MATCHTRD_TPC", ss, Lifetime::Timeframe}, mTracks);
    pc.outputs().snapshot(Output{o2::header::gDataOriginTRD, "TRKTRG_TPC", ss, Lifetime::Timeframe}, mTrigRec);
  } else if (mMode == Mode::ITSTPCTRD) {
    pc.outputs().snapshot(Output{o2::header::gDataOriginTRD, "MATCHTRD_GLO", 0, Lifetime::Timeframe}, mTracks);
    pc.outputs().snapshot(Output{o2::header::gDataOriginTRD, "TRKTRG_GLO", 0, Lifetime::Timeframe}, mTrigRec);
  }

  if (mTree->GetReadEntry() + 1 >= mTree->GetEntries()) {
    pc.services().get<ControlService>().endOfStream();
    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
  }
}

void TRDTrackReader::connectTree(const std::string& filename)
{
  mTree.reset(nullptr); // in case it was already loaded
  mFile.reset(TFile::Open(filename.c_str()));
  assert(mFile && !mFile->IsZombie());
  mTree.reset((TTree*)mFile->Get("tracksTRD"));
  assert(mTree);
  mTree->SetBranchAddress("tracks", &mTracksPtr);
  mTree->SetBranchAddress("trgrec", &mTrigRecPtr);
  LOG(INFO) << "Loaded tree from " << filename << " with " << mTree->GetEntries() << " entries";
}

DataProcessorSpec getTRDTPCTrackReaderSpec(bool useMC, bool subSpecStrict)
{
  std::vector<OutputSpec> outputs;
  uint32_t sspec = o2::globaltracking::getSubSpec(subSpecStrict ? o2::globaltracking::MatchingType::Strict : o2::globaltracking::MatchingType::Standard);
  outputs.emplace_back(o2::header::gDataOriginTRD, "MATCHTRD_TPC", sspec, Lifetime::Timeframe);
  outputs.emplace_back(o2::header::gDataOriginTRD, "TRKTRG_TPC", sspec, Lifetime::Timeframe);

  if (useMC) {
    LOG(FATAL) << "TRD track reader cannot read MC data (yet)";
  }
  return DataProcessorSpec{
    "tpctrd-track-reader",
    Inputs{},
    outputs,
    AlgorithmSpec{adaptFromTask<TRDTrackReader>(useMC, TRDTrackReader::Mode::TPCTRD, subSpecStrict)},
    Options{
      {"track-infile", VariantType::String, "trdmatches_tpc.root", {"Name of the input file for TPC-TRD matches"}},
      {"input-dir", VariantType::String, "none", {"Input directory"}}}};
}

DataProcessorSpec getTRDGlobalTrackReaderSpec(bool useMC)
{
  std::vector<OutputSpec> outputs;
  outputs.emplace_back(o2::header::gDataOriginTRD, "MATCHTRD_GLO", 0, Lifetime::Timeframe);
  outputs.emplace_back(o2::header::gDataOriginTRD, "TRKTRG_GLO", 0, Lifetime::Timeframe);

  if (useMC) {
    LOG(FATAL) << "TRD track reader cannot read MC data (yet)";
  }

  return DataProcessorSpec{
    "itstpctrd-track-reader",
    Inputs{},
    outputs,
    AlgorithmSpec{adaptFromTask<TRDTrackReader>(useMC, TRDTrackReader::Mode::ITSTPCTRD)},
    Options{
      {"track-infile", VariantType::String, "trdmatches_itstpc.root", {"Name of the input file for ITS-TPC-TRD matches"}},
      {"input-dir", VariantType::String, "none", {"Input directory"}}}};
}

} // namespace trd
} // namespace o2
