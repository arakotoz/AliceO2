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

/// \file AlignmentSpec.h
/// \author arakotoz@cern.ch
/// \brief Class to run standalone alignment for MFT in a workflow

#ifndef ALICEO2_MFT_ALIGNMENT_DEVICE_H
#define ALICEO2_MFT_ALIGNMENT_DEVICE_H

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "MFTAlignment/Alignment.h"
#include "TStopwatch.h"
#include "DetectorsBase/GRPGeomHelper.h"

using namespace o2::framework;

namespace o2
{
namespace mft
{
class AlignmentSpec : public Task
{
 public:
  AlignmentSpec(std::shared_ptr<o2::base::GRPGeomRequest> gr)
    : mGGCCDBRequest(gr){};
  void init(o2::framework::InitContext& ic) final;
  void run(o2::framework::ProcessingContext& pc) final;
  void endOfStream(o2::framework::EndOfStreamContext& ec) final;
  void finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void updateTimeDependentParams(o2::framework::ProcessingContext& pc);
  void sendOutput(o2::framework::DataAllocator& output);
  std::unique_ptr<o2::mft::Alignment> mAlignment;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  enum TimerIDs { SWTot,
                  SWProcessTimeFrame,
                  SWProcessRecoTracks,
                  SWGlobalFit,
                  NStopWatches };
  static constexpr std::string_view TimerName[] = {"Total",
                                                   "processTimeFrame",
                                                   "processRecoTracks",
                                                   "globalFit"};
  TStopwatch mTimer[NStopWatches];
};

DataProcessorSpec getAlignmentSpec();

} // namespace mft
} // namespace o2

#endif
