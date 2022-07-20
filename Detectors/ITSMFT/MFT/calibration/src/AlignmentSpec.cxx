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

/// @file AlignmentSpec.cxx

#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/Logger.h"
#include "MFTCalibration/AlignmentSpec.h"
#include "CommonUtils/NameConf.h"

using namespace o2::framework;

namespace o2
{
namespace mft
{

//_____________________________________________________________
void AlignmentSpec::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
  mAlignment = std::make_unique<o2::mft::Alignment>();
  for (int sw = 0; sw < NStopWatches; sw++) {
    mTimer[sw].Stop();
    mTimer[sw].Reset();
  }
  mTimer[SWTot].Start(false);
  mAlignment->init();
}

//_____________________________________________________________
void AlignmentSpec::run(o2::framework::ProcessingContext& pc)
{
  updateTimeDependentParams(pc);
  mTimer[SWProcessTimeFrame].Start(false);
  mAlignment->processTimeFrame(pc);
  mTimer[SWProcessTimeFrame].Stop();

  mTimer[SWProcessRecoTracks].Start(false);
  mAlignment->processRecoTracks();
  mTimer[SWProcessRecoTracks].Stop();
}

//_____________________________________________________________
void AlignmentSpec::endOfStream(o2::framework::EndOfStreamContext& ec)
{

  mTimer[SWGlobalFit].Start(false);
  mAlignment->globalFit();
  mTimer[SWGlobalFit].Stop();

  sendOutput(ec.outputs());
  mTimer[SWTot].Stop();

  for (int i = 0; i < NStopWatches; i++) {
    LOGF(info, "Timing %18s: Cpu: %.3e s; Real: %.3e s in %d slots", TimerName[i], mTimer[i].CpuTime(), mTimer[i].RealTime(), mTimer[i].Counter() - 1);
  }
}

//_____________________________________________________________
void AlignmentSpec::sendOutput(DataAllocator& output)
{
  TObjArray objar;
  mAlignment->getHistos(objar);

  output.snapshot(Output{"MFT", "MFTASSESSMENT", 0, Lifetime::Sporadic}, objar);

  TFile* f = new TFile(Form("mft_alignment_records.root"), "RECREATE");
  objar.Write();
  f->Close();
}
///_______________________________________
void AlignmentSpec::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    pc.inputs().get<o2::itsmft::TopologyDictionary*>("cldict"); // just to trigger the finaliseCCDB

    auto field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
    auto Bz = field->getBz(centerMFT);
    LOG(info) << "Setting MFT Assessment Bz = " << Bz;
    mAlignment->setBz(Bz);
    mAlignment->init(mFinalizeAnalysis);
  }
}

///_______________________________________
void AlignmentSpec::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("MFT", "CLUSDICT", 0)) {
    LOG(info) << "cluster dictionary updated";
    mAlignment->setClusterDictionary((const o2::itsmft::TopologyDictionary*)obj);
  }
}

//_____________________________________________________________
DataProcessorSpec getAlignmentSpec(bool useMC, bool processGen, bool finalizeAnalysis)
{
  std::vector<InputSpec> inputs;
  std::vector<OutputSpec> outputs;

  inputs.emplace_back("compClusters", "MFT", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  inputs.emplace_back("patterns", "MFT", "PATTERNS", 0, Lifetime::Timeframe);
  inputs.emplace_back("clustersrofs", "MFT", "CLUSTERSROF", 0, Lifetime::Timeframe);
  inputs.emplace_back("tracksrofs", "MFT", "MFTTrackROF", 0, Lifetime::Timeframe);
  inputs.emplace_back("tracks", "MFT", "TRACKS", 0, Lifetime::Timeframe);
  inputs.emplace_back("trackClIdx", "MFT", "TRACKCLSID", 0, Lifetime::Timeframe);
  inputs.emplace_back("cldict", "MFT", "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec("MFT/Calib/ClusterDictionary"));
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              true,                              // GRPMagField
                                                              false,                             // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              inputs,
                                                              true);
  if (useMC) {
    inputs.emplace_back("clslabels", "MFT", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
    inputs.emplace_back("trklabels", "MFT", "TRACKSMCTR", 0, Lifetime::Timeframe);
  }

  outputs.emplace_back("MFT", "MFTASSESSMENT", 0, Lifetime::Sporadic);

  return DataProcessorSpec{
    "mft-alignment",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<o2::mft::AlignmentSpec>(ggRequest)},
    Options{{}}};
}

} // namespace mft
} // namespace o2
