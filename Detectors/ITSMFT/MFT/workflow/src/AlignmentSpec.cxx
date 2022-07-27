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
#include "Field/MagneticField.h"
#include "TGeoGlobalMagField.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTAlignment/AlignConfig.h"

#include "MFTWorkflow/AlignmentSpec.h"
#include "CommonUtils/NameConf.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

using namespace o2::framework;

namespace o2
{
namespace mft
{

//_____________________________________________________________
void AlignmentSpec::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);

  auto& alignConfigParam = o2::mft::AlignConfig::Instance();
  mAlignment = std::make_unique<o2::mft::Alignment>();
  mAlignment->setChi2CutNStdDev(alignConfigParam.chi2CutNStdDev);
  mAlignment->setResidualCutInitial(alignConfigParam.residualCutInitial);
  mAlignment->setResidualCut(alignConfigParam.residualCut);
  mAlignment->setAllowedVariationDeltaX(alignConfigParam.allowedVarDeltaX);
  mAlignment->setAllowedVariationDeltaY(alignConfigParam.allowedVarDeltaY);
  mAlignment->setAllowedVariationDeltaZ(alignConfigParam.allowedVarDeltaZ);
  mAlignment->setAllowedVariationDeltaRz(alignConfigParam.allowedVarDeltaRz);
  for (int sw = 0; sw < NStopWatches; sw++) {
    mTimer[sw].Stop();
    mTimer[sw].Reset();
  }
  mTimer[SWTot].Start(false);
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
  mAlignment->printProcessTrackSummary();

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
  std::vector<o2::detectors::AlignParam> alignParams;
  mAlignment->getAlignParams(alignParams);
  output.snapshot(Output{"MFT", "MFTALIGNMENT", 0, Lifetime::Sporadic}, alignParams);
  LOG(info) << "Storing MFT alignment params in local file mft_alignment.root";
  TFile* f = new TFile(Form("mft_alignment.root"), "RECREATE");
  f->WriteObjectAny(&alignParams, "std::vector<o2::detectors::AlignParam>", "alignment");
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

    o2::mft::GeometryTGeo* geom = o2::mft::GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                   o2::math_utils::TransformType::T2GRot,
                                                   o2::math_utils::TransformType::T2G));
    mAlignment->setGeometry(geom);
    auto& alignConfigParam = o2::mft::AlignConfig::Instance();
    auto field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
    auto Bz = field->getBz(centerMFT);
    LOG(info) << "Setting MFT Assessment Bz = " << Bz;
    mAlignment->setBz(Bz);
    mAlignment->setMinNumberClusterCut(alignConfigParam.minPoints[o2::mft::AlignConfig::Collision]);
    mAlignment->init();
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
DataProcessorSpec getAlignmentSpec()
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

  outputs.emplace_back("MFT", "MFTALIGNMENT", 0, Lifetime::Sporadic);

  return DataProcessorSpec{
    "mft-alignment",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<o2::mft::AlignmentSpec>(ggRequest)},
    Options{{}}};
}

} // namespace mft
} // namespace o2
