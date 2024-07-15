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

/// @file   TrackerSpec.cxx

#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/CCDBParamSpec.h"
#include "ITS3Workflow/TrackerSpec.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITS/TrackITS.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/PhysTrigger.h"

#include "ITStracking/TrackingConfigParam.h"
#include "ITSMFTBase/DPLAlpideParam.h"

#include "ITSBase/GeometryTGeo.h"
#include "CommonDataFormat/IRFrame.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "ITS3Reconstruction/IOUtils.h"
#include "ITSReconstruction/FastMultEstConfig.h"
#include "ITS3Base/SpecsV2.h"

namespace o2
{
using namespace framework;
using its::FastMultEstConfig;
using its::TimeFrame;
using its::Tracker;
using its::TrackingParameters;
using its::TrackITSExt;
using its::Vertexer;

namespace its3
{
using Vertex = o2::dataformats::Vertex<o2::dataformats::TimeStamp<int>>;

TrackerDPL::TrackerDPL(std::shared_ptr<o2::base::GRPGeomRequest> gr, bool isMC, int trgType, const std::string& trModeS, o2::gpu::GPUDataTypes::DeviceType dType) : mGGCCDBRequest(gr), mIsMC{isMC}, mUseTriggers{trgType}, mMode{trModeS}, mRecChain{o2::gpu::GPUReconstruction::CreateInstance(dType, true)}
{
  std::transform(mMode.begin(), mMode.end(), mMode.begin(), [](unsigned char c) { return std::tolower(c); });
}

void TrackerDPL::init(InitContext& /*ic*/)
{
  mTimer.Stop();
  mTimer.Reset();

  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);

  mChainITS.reset(mRecChain->AddChain<o2::gpu::GPUChainITS>());
  mVertexer = std::make_unique<Vertexer>(mChainITS->GetITSVertexerTraits());
  mTracker = std::make_unique<Tracker>(mChainITS->GetITSTrackerTraits());
  mRunVertexer = true;
  mCosmicsProcessing = false;
  std::vector<TrackingParameters> trackParams;

  if (mMode == "async") {

    trackParams.resize(3);
    trackParams[1].TrackletMinPt = 0.2f;
    trackParams[1].CellDeltaTanLambdaSigma *= 2.;
    trackParams[2].TrackletMinPt = 0.1f;
    trackParams[2].CellDeltaTanLambdaSigma *= 4.;
    trackParams[2].MinTrackLength = 4;
    LOG(info) << "Initializing tracker in async. phase reconstruction with " << trackParams.size() << " passes";

  } else if (mMode == "sync_misaligned") {

    trackParams.resize(3);
    trackParams[0].PhiBins = 32;
    trackParams[0].ZBins = 64;
    trackParams[0].CellDeltaTanLambdaSigma *= 3.;
    trackParams[0].SystErrorZ2[0] = 1.e-4;
    trackParams[0].SystErrorZ2[1] = 1.e-4;
    trackParams[0].SystErrorZ2[2] = 1.e-4;
    trackParams[0].SystErrorZ2[3] = 9.e-4;
    trackParams[0].SystErrorZ2[4] = 9.e-4;
    trackParams[0].SystErrorZ2[5] = 9.e-4;
    trackParams[0].SystErrorZ2[6] = 9.e-4;
    std::copy(trackParams[0].SystErrorZ2.begin(), trackParams[0].SystErrorZ2.end(), trackParams[0].SystErrorY2.begin());
    trackParams[0].MaxChi2ClusterAttachment = 60.;
    trackParams[0].MaxChi2NDF = 40.;
    trackParams[1] = trackParams[0];
    trackParams[2] = trackParams[0];
    trackParams[1].MinTrackLength = 6;
    trackParams[2].MinTrackLength = 4;
    LOG(info) << "Initializing tracker in misaligned sync. phase reconstruction with " << trackParams.size() << " passes";

  } else if (mMode == "sync") {
    trackParams.resize(1);
    LOG(info) << "Initializing tracker in sync. phase reconstruction with " << trackParams.size() << " passes";
  } else if (mMode == "cosmics") {
    mCosmicsProcessing = true;
    mRunVertexer = false;
    trackParams.resize(1);
    trackParams[0].MinTrackLength = 4;
    trackParams[0].CellDeltaTanLambdaSigma *= 10;
    trackParams[0].PhiBins = 4;
    trackParams[0].ZBins = 16;
    trackParams[0].PVres = 1.e5f;
    trackParams[0].MaxChi2ClusterAttachment = 60.;
    trackParams[0].MaxChi2NDF = 40.;
    trackParams[0].TrackletsPerClusterLimit = 100.;
    trackParams[0].CellsPerClusterLimit = 100.;
    LOG(info) << "Initializing tracker in reconstruction for cosmics with " << trackParams.size() << " passes";

  } else {
    throw std::runtime_error(fmt::format("Unsupported ITS tracking mode {:s} ", mMode));
  }

  for (auto& params : trackParams) {
    params.CorrType = o2::base::PropagatorImpl<float>::MatCorrType::USEMatCorrLUT;
    for (int iLayer{0}; iLayer < params.NLayers - 4; ++iLayer) { // initialise ITS3 radii and lengths
      params.LayerZ[iLayer] = constants::segment::lengthSensitive;
      params.LayerRadii[iLayer] = constants::radii[iLayer];
    }
  }

  mTracker->setParameters(trackParams);
}

void TrackerDPL::run(ProcessingContext& pc)
{
  mTimer.Start(false);
  updateTimeDependentParams(pc);
  auto compClusters = pc.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>("compClusters");
  gsl::span<const unsigned char> patterns = pc.inputs().get<gsl::span<unsigned char>>("patterns");
  gsl::span<const o2::itsmft::PhysTrigger> physTriggers;
  std::vector<o2::itsmft::PhysTrigger> fromTRD;
  if (mUseTriggers == 2) { // use TRD triggers
    o2::InteractionRecord ir{0, pc.services().get<o2::framework::TimingInfo>().firstTForbit};
    auto trdTriggers = pc.inputs().get<gsl::span<o2::trd::TriggerRecord>>("phystrig");
    for (const auto& trig : trdTriggers) {
      if (trig.getBCData() >= ir && trig.getNumberOfTracklets()) {
        ir = trig.getBCData();
        fromTRD.emplace_back(o2::itsmft::PhysTrigger{ir, 0});
      }
    }
    physTriggers = gsl::span<const o2::itsmft::PhysTrigger>(fromTRD.data(), fromTRD.size());
  } else if (mUseTriggers == 1) { // use Phys triggers from ITS stream
    physTriggers = pc.inputs().get<gsl::span<o2::itsmft::PhysTrigger>>("phystrig");
  }

  // code further down does assignment to the rofs and the altered object is used for output
  // we therefore need a copy of the vector rather than an object created directly on the input data,
  // the output vector however is created directly inside the message memory thus avoiding copy by
  // snapshot
  auto orig = o2::header::gDataOriginITS;
  auto rofsinput = pc.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("ROframes");
  auto& rofs = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{orig, "ITSTrackROF", 0}, rofsinput.begin(), rofsinput.end());

  auto& irFrames = pc.outputs().make<std::vector<o2::dataformats::IRFrame>>(Output{orig, "IRFRAMES", 0});

  const auto& alpParams = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance(); // RS: this should come from CCDB
  int nBCPerTF = alpParams.roFrameLengthInBC;

  LOG(info) << "ITS3Tracker pulled " << compClusters.size() << " clusters, " << rofs.size() << " RO frames";

  const dataformats::MCTruthContainer<MCCompLabel>* labels = nullptr;
  gsl::span<itsmft::MC2ROFRecord const> mc2rofs;
  if (mIsMC) {
    labels = pc.inputs().get<const dataformats::MCTruthContainer<MCCompLabel>*>("labels").release();
    // get the array as read-only span, a snapshot is send forward
    mc2rofs = pc.inputs().get<gsl::span<itsmft::MC2ROFRecord>>("MC2ROframes");
    LOG(info) << labels->getIndexedSize() << " MC label objects , in " << mc2rofs.size() << " MC events";
  }

  auto& allClusIdx = pc.outputs().make<std::vector<int>>(Output{orig, "TRACKCLSID", 0});
  std::vector<o2::MCCompLabel> trackLabels;
  std::vector<MCCompLabel> verticesLabels;
  auto& allTracks = pc.outputs().make<std::vector<o2::its::TrackITS>>(Output{orig, "TRACKS", 0});
  std::vector<o2::MCCompLabel> allTrackLabels;
  std::vector<o2::MCCompLabel> allVerticesLabels;

  auto& vertROFvec = pc.outputs().make<std::vector<o2::itsmft::ROFRecord>>(Output{orig, "VERTICESROF", 0});
  auto& vertices = pc.outputs().make<std::vector<Vertex>>(Output{orig, "VERTICES", 0});

  TimeFrame* timeFrame = mChainITS->GetITSTimeframe();
  timeFrame->resizeVectors(7);
  mTracker->adoptTimeFrame(*timeFrame);

  mTracker->setBz(o2::base::Propagator::Instance()->getNominalBz());
  mVertexer->adoptTimeFrame(*timeFrame);
  gsl::span<const unsigned char>::iterator pattIt = patterns.begin();

  gsl::span<itsmft::ROFRecord> rofspan(rofs);
  ioutils::loadROFrameDataITS3(timeFrame, rofspan, compClusters, pattIt, mDict, labels);
  pattIt = patterns.begin();
  std::vector<int> savedROF;
  auto logger = [&](const std::string& s) { LOG(info) << s; };
  auto errorLogger = [&](const std::string& s) { LOG(error) << s; };

  // o2::its3::FastMultEst multEst;                        // mult estimator
  std::vector<bool> processingMask(rofs.size(), true); // Override mult estimator
  int cutVertexMult{0};
  // int cutRandomMult = int(rofs.size()) - multEst.selectROFs(rofs, compClusters, physTriggers, processingMask);
  timeFrame->setMultiplicityCutMask(processingMask);
  float vertexerElapsedTime{0.f};
  if (mRunVertexer) {
    // Run seeding vertexer
    vertexerElapsedTime = mVertexer->clustersToVertices(logger);
  }
  const auto& multEstConf = FastMultEstConfig::Instance(); // parameters for mult estimation and cuts
  for (size_t iRof{0}; iRof < rofspan.size(); ++iRof) {
    std::vector<Vertex> vtxVecLoc;
    auto& vtxROF = vertROFvec.emplace_back(rofspan[iRof]);
    vtxROF.setFirstEntry(vertices.size());
    if (mRunVertexer) {
      auto vtxSpan = timeFrame->getPrimaryVertices(iRof);
      vtxROF.setNEntries(vtxSpan.size());
      bool selROF = vtxSpan.size() == 0;
      for (size_t iV{0}; iV < vtxSpan.size(); ++iV) {
        auto& v = vtxSpan[iV];
        if (multEstConf.isVtxMultCutRequested() && !multEstConf.isPassingVtxMultCut(v.getNContributors())) {
          continue; // skip vertex of unwanted multiplicity
        }
        selROF = true;
        vertices.push_back(v);
        if (mIsMC) {
          auto vLabels = timeFrame->getPrimaryVerticesLabels(iRof)[iV];
          std::copy(vLabels.begin(), vLabels.end(), std::back_inserter(allVerticesLabels));
        }
      }
      if (processingMask[iRof] && !selROF) { // passed selection in clusters and not in vertex multiplicity
                                             //         LOG(debug) << fmt::format("ROF {} rejected by the vertex multiplicity selection [{},{}]",
                                             //                                   iRof,
                                             //                                   multEstConf.cutMultVtxLow,
                                             //                                   multEstConf.cutMultVtxHigh);
        processingMask[iRof] = selROF;
        cutVertexMult++;
      }
    } else { // cosmics
      vtxVecLoc.emplace_back(Vertex());
      vtxVecLoc.back().setNContributors(1);
      vtxROF.setNEntries(vtxVecLoc.size());
      for (auto& v : vtxVecLoc) {
        vertices.push_back(v);
      }
      timeFrame->addPrimaryVertices(vtxVecLoc);
    }
  }
  // LOG(info) << fmt::format(" - rejected {}/{} ROFs: random/mult.sel:{} (seed {}), vtx.sel:{}", cutRandomMult + cutVertexMult, rofspan.size(), cutRandomMult, multEst.lastRandomSeed, cutVertexMult);
  LOG(info) << fmt::format(" - Vertex seeding total elapsed time: {} ms for {} vertices found in {} ROFs", vertexerElapsedTime, timeFrame->getPrimaryVerticesNum(), rofspan.size());
  LOG(info) << fmt::format(" - Beam position computed for the TF: {}, {}", timeFrame->getBeamX(), timeFrame->getBeamY());

  if (mCosmicsProcessing && compClusters.size() > 1500 * rofspan.size()) {
    LOG(error) << "Cosmics processing was requested with an average detector occupancy exceeding 1.e-7, skipping TF processing.";
  } else {

    timeFrame->setMultiplicityCutMask(processingMask);
    // Run CA tracker
    mTracker->clustersToTracks(logger, errorLogger);
    if (timeFrame->hasBogusClusters() != 0) {
      LOG(warning) << fmt::format(" - The processed timeframe had {} clusters with wild z coordinates, check the dictionaries", timeFrame->hasBogusClusters());
    }

    for (unsigned int iROF{0}; iROF < rofs.size(); ++iROF) {
      auto& rof{rofs[iROF]};
      auto& tracks = timeFrame->getTracks(iROF);
      trackLabels = timeFrame->getTracksLabel(iROF);
      auto number{tracks.size()};
      auto first{allTracks.size()};
      rof.setFirstEntry(first);
      rof.setNEntries(number);
      if (processingMask[iROF]) {
        irFrames.emplace_back(rof.getBCData(), rof.getBCData() + nBCPerTF - 1).info = tracks.size();
      }

      std::copy(trackLabels.begin(), trackLabels.end(), std::back_inserter(allTrackLabels));
      // Some conversions that needs to be moved in the tracker internals
      for (unsigned int iTrk{0}; iTrk < tracks.size(); ++iTrk) {
        auto& trc{tracks[iTrk]};
        trc.setFirstClusterEntry(allClusIdx.size()); // before adding tracks, create final cluster indices
        int nclf = 0, ncl = (int)allClusIdx.size();
        for (int ic = TrackITSExt::MaxClusters; (ic--) != 0;) { // track internally keeps in->out cluster indices, but we want to store the references as out->in!!!
          auto clid = trc.getClusterIndex(ic);
          if (clid >= 0) {
            trc.setClusterSize(ic, timeFrame->getClusterSize(clid));
            allClusIdx.push_back(clid);
            nclf++;
          }
        }
        assert(ncl == nclf);
        allTracks.emplace_back(trc);
      }
    }

    LOGP(info, "ITS3Tracker pushed {} tracks and {} vertices", allTracks.size(), vertices.size());
    if (mIsMC) {
      LOGP(info, "ITS3Tracker pushed {} track labels", allTrackLabels.size());
      LOGP(info, "ITS3Tracker pushed {} vertex labels", allVerticesLabels.size());

      pc.outputs().snapshot(Output{orig, "TRACKSMCTR", 0}, allTrackLabels);
      pc.outputs().snapshot(Output{orig, "VERTICESMCTR", 0}, allVerticesLabels);
      pc.outputs().snapshot(Output{orig, "ITSTrackMC2ROF", 0}, mc2rofs);
    }
  }
  mTimer.Stop();
}

///_______________________________________
void TrackerDPL::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    pc.inputs().get<o2::its3::TopologyDictionary*>("cldict"); // just to trigger the finaliseCCDB
    pc.inputs().get<o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>*>("alppar");

    // Check if lightweight geometry was requested, otherwise full geometry is loaded
    if (pc.inputs().getPos("itsTGeo") >= 0) {
      pc.inputs().get<o2::its::GeometryTGeo*>("itsTGeo");
    }
    o2::its::GeometryTGeo* geom = o2::its::GeometryTGeo::Instance();
    geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G));
    mVertexer->getGlobalConfiguration();
    mTracker->getGlobalConfiguration();
  }
}

///_______________________________________
void TrackerDPL::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("IT3", "CLUSDICT", 0)) {
    LOG(info) << "cluster dictionary updated";
    setClusterDictionary((o2::its3::TopologyDictionary*)obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "GEOMTGEO", 0)) {
    LOG(info) << "IT3 GeometryTGeo loaded from ccdb";
    o2::its::GeometryTGeo::adopt((o2::its::GeometryTGeo*)obj);
    return;
  }
  // Note: strictly speaking, for Configurable params we don't need finaliseCCDB check, the singletons are updated at the CCDB fetcher level
  if (matcher == ConcreteDataMatcher("ITS", "ALPIDEPARAM", 0)) {
    LOG(info) << "Alpide param updated";
    const auto& par = o2::itsmft::DPLAlpideParam<o2::detectors::DetID::ITS>::Instance();
    par.printKeyValues();
    return;
  }
}

void TrackerDPL::endOfStream(EndOfStreamContext& ec)
{
  LOGF(info, "ITS3 CA-Tracker total timing: Cpu: %.3e Real: %.3e s in %d slots",
       mTimer.CpuTime(), mTimer.RealTime(), mTimer.Counter() - 1);
}

DataProcessorSpec getTrackerSpec(bool useMC, bool useGeom, const int trgType, const std::string& trModeS, o2::gpu::GPUDataTypes::DeviceType dType)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("compClusters", "ITS", "COMPCLUSTERS", 0, Lifetime::Timeframe);
  inputs.emplace_back("patterns", "ITS", "PATTERNS", 0, Lifetime::Timeframe);
  inputs.emplace_back("ROframes", "ITS", "CLUSTERSROF", 0, Lifetime::Timeframe);
  if (trgType == 1) {
    inputs.emplace_back("phystrig", "ITS", "PHYSTRIG", 0, Lifetime::Timeframe);
  } else if (trgType == 2) {
    inputs.emplace_back("phystrig", "TRD", "TRKTRGRD", 0, Lifetime::Timeframe);
  }
  inputs.emplace_back("cldict", "IT3", "CLUSDICT", 0, Lifetime::Condition, ccdbParamSpec("IT3/Calib/ClusterDictionary"));
  inputs.emplace_back("alppar", "ITS", "ALPIDEPARAM", 0, Lifetime::Condition, ccdbParamSpec("ITS/Config/AlpideParam"));
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                                                                        // orbitResetTime
                                                              false,                                                                        // GRPECS
                                                              false,                                                                        // GRPLHCIF
                                                              true,                                                                         // GRPMagField
                                                              true,                                                                         // askMatLUT
                                                              useGeom ? o2::base::GRPGeomRequest::Aligned : o2::base::GRPGeomRequest::None, // geometry
                                                              inputs,
                                                              true);

  if (!useGeom) { // load light-weight geometry
    inputs.emplace_back("itsTGeo", "ITS", "GEOMTGEO", 0, Lifetime::Condition, ccdbParamSpec("ITS/Config/Geometry"));
  }

  std::vector<OutputSpec> outputs;
  outputs.emplace_back("ITS", "TRACKS", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "TRACKCLSID", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "ITSTrackROF", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "VERTICES", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "VERTICESROF", 0, Lifetime::Timeframe);
  outputs.emplace_back("ITS", "IRFRAMES", 0, Lifetime::Timeframe);

  if (useMC) {
    inputs.emplace_back("labels", "ITS", "CLUSTERSMCTR", 0, Lifetime::Timeframe);
    inputs.emplace_back("MC2ROframes", "ITS", "CLUSTERSMC2ROF", 0, Lifetime::Timeframe);
    outputs.emplace_back("ITS", "VERTICESMCTR", 0, Lifetime::Timeframe);
    outputs.emplace_back("ITS", "TRACKSMCTR", 0, Lifetime::Timeframe);
    outputs.emplace_back("ITS", "ITSTrackMC2ROF", 0, Lifetime::Timeframe);
  }

  return DataProcessorSpec{
    "its3-tracker",
    inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<TrackerDPL>(ggRequest, useMC, trgType, trModeS, dType)},
    Options{}};
}

} // namespace its3
} // namespace o2
