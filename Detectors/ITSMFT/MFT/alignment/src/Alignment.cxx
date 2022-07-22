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

/// @file Alignment.cxx

#include <iostream>
#include <sstream>
#include <string>

#include <TString.h>

#include "Framework/InputSpec.h"
#include "Framework/Logger.h"
#include <Framework/InputRecord.h>
#include "MFTAlignment/AlignPointHelper.h"
#include "MFTAlignment/AlignSensorHelper.h"
#include "MFTAlignment/Alignment.h"
#include "MFTTracking/IOUtils.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTAlignment/MillePedeRecord.h"
#include "MFTAlignment/MillePede2.h"

using namespace o2::mft;

ClassImp(o2::mft::Alignment);

//__________________________________________________________________________
Alignment::Alignment()
  : mRunNumber(0),
    mBz(0),
    mNumberTFs(0),
    mCounterLocalEquationFailed(0),
    mCounterSkippedTracks(0),
    mCounterUsedTracks(0),
    mStartFac(256),
    mChi2CutNStdDev(3),
    mResCutInitial(100.),
    mResCut(100.),
    mMinNumberClusterCut(6),
    mWeightRecord(1.),
    mSaveTrackRecordToFile(false),
    mMilleRecordsFileName("mft_mille_records.root"),
    mMilleConstraintsRecFileName("mft_mille_constraints.root"),
    mIsInitDone(false)
{
  mMillepede = std::make_unique<MillePede2>();
  // default allowed variations w.r.t. global system coordinates
  mAllowVar[0] = 0.5;  // delta translation in x (cm)
  mAllowVar[1] = 0.5;  // delta translation in y (cm)
  mAllowVar[2] = 0.01; // rotation angle Rz around z-axis (rad)
  mAllowVar[3] = 0.5;  // delta translation in z (cm)
}

//__________________________________________________________________________
void Alignment::init()
{
  if (mIsInitDone)
    return;
  if (mGeometry == nullptr) {
    LOGF(fatal, "Alignment::init() failed because no geometry is defined");
    mIsInitDone = false;
    return;
  }
  mAlignPoint = std::make_unique<AlignPointHelper>(mGeometry);
  mMillepede->InitMille(mNumberOfGlobalParam,
                        mNumberOfTrackParam,
                        mChi2CutNStdDev,
                        mResCut,
                        mResCutInitial);
  // filenames for the processed tracks and constraints records
  mMillepede->SetDataRecFName(mMilleRecordsFileName.Data());
  mMillepede->SetConsRecFName(mMilleConstraintsRecFileName.Data());

  if (mSaveTrackRecordToFile) {
    mMillepede->InitDataRecStorage(kFALSE);
  }

  LOGF(info,
       "Allowed variation: dx %.3f, dy %.3f, dz %.3f, dRz %.4f",
       mAllowVar[0], mAllowVar[1], mAllowVar[3], mAllowVar[2]);

  // set allowed variations for all parameters
  for (int chipId = 0; chipId < mNumberOfSensors; ++chipId) {
    for (Int_t iPar = 0; iPar < mNDofPerSensor; ++iPar) {
      mMillepede->SetParSigma(chipId * mNDofPerSensor + iPar, mAllowVar[iPar]);
    }
  }

  // set iterations
  if (mStartFac > 1) {
    mMillepede->SetIterations(mStartFac);
  }

  LOGF(info, "Alignment init done");
  mIsInitDone = true;
}

//__________________________________________________________________________
void Alignment::processTimeFrame(o2::framework::ProcessingContext& ctx)
{
  mNumberTFs++; // TF Counter

  // get tracks
  mMFTTracks = ctx.inputs().get<gsl::span<o2::mft::TrackMFT>>("tracks");
  mMFTTracksROF = ctx.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("tracksrofs");
  mMFTTrackClusIdx = ctx.inputs().get<gsl::span<int>>("trackClIdx");

  // get clusters
  mMFTClusters = ctx.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>("compClusters");
  mMFTClustersROF = ctx.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("clustersrofs");
  mMFTClusterPatterns = ctx.inputs().get<gsl::span<unsigned char>>("patterns");
  pattIt = mMFTClusterPatterns.begin();
  mMFTClustersGlobal.clear();
  mMFTClustersGlobal.reserve(mMFTClusters.size());
  o2::mft::ioutils::convertCompactClusters(
    mMFTClusters, pattIt, mMFTClustersGlobal, mDictionary);
}

//__________________________________________________________________________
void Alignment::processRecoTracks()
{
  if (!mIsInitDone) {
    LOGF(fatal, "Alignment::processRecoTracks() aborted because init was not done!");
    return;
  }

  for (auto oneTrack : mMFTTracks) { // track loop

    // Skip the track if not enough clusters
    auto ncls = oneTrack.getNumberOfPoints();
    if (ncls < mMinNumberClusterCut) {
      mCounterSkippedTracks++;
      continue;
    }

    auto offset = oneTrack.getExternalClusterIndexOffset();

    mTrackRecord.Reset();
    if (mMillepede->GetRecord()) {
      mMillepede->GetRecord()->Reset();
    }

    // Store the initial track parameters
    mAlignPoint->resetTrackInitialParam();
    mAlignPoint->recordTrackInitialParam(oneTrack);

    for (int icls = 0; icls < ncls; ++icls) { // cluster loop

      mAlignPoint->resetDerivatives();
      mAlignPoint->resetAlignPoint();

      auto clsEntry = mMFTTrackClusIdx[offset + icls];
      auto globalCluster = mMFTClustersGlobal[clsEntry];

      // Propagate track to the current z plane of this cluster
      oneTrack.propagateParamToZlinear(globalCluster.getZ());

      // Store reco and measured positions
      mAlignPoint->setGlobalRecoPosition(oneTrack);
      mAlignPoint->setLocalMeasuredPosition(globalCluster);

      // Compute derivatives
      mAlignPoint->computeLocalDerivatives();
      mAlignPoint->computeGlobalDerivatives();

      // Set local equations
      bool success = true;
      success &= setLocalEquationX();
      success &= setLocalEquationY();
      success &= setLocalEquationZ();
      if (!success)
        mCounterLocalEquationFailed++;

    } // end of loop on clusters

    // copy track record
    mMillepede->SetRecordRun(mRunNumber);
    mMillepede->SetRecordWeight(mWeightRecord);
    mTrackRecord = *mMillepede->GetRecord();

    // save record data
    if (mSaveTrackRecordToFile) {
      mMillepede->SaveRecordData();
    }

    mCounterUsedTracks++;
  } // end of loop on tracks

  if (mSaveTrackRecordToFile) {
    mMillepede->CloseDataRecStorage();
  }
}

//__________________________________________________________________________
void Alignment::globalFit()
{
  if (!mIsInitDone) {
    LOGF(fatal, "Alignment::globalFit() aborted because init was not done!");
    return;
  }

  mMillepede->GlobalFit(mAlignParam, mAlignParamErrors, mAlignParamPulls);

  LOGF(info, "Alignment: done fitting global parameters");
  LOGF(info, "sensor info, dx (cm), dy (cm), dz (cm), dRz (rad)");

  AlignSensorHelper chipHelper(mGeometry);
  bool wSymName = false;

  for (int chipId = 0; chipId < mNumberOfSensors; chipId++) {
    chipHelper.setSensorOnlyInfo(chipId);
    std::stringstream name = chipHelper.getSensorFullName(wSymName);
    LOGF(info, "%s, %.3e, %.3e, %.3e, %.3e",
         name.str().c_str(),
         mAlignParam[chipId * mNDofPerSensor + 0],
         mAlignParam[chipId * mNDofPerSensor + 1],
         mAlignParam[chipId * mNDofPerSensor + 3],
         mAlignParam[chipId * mNDofPerSensor + 2]);
  }
}

//__________________________________________________________________________
void Alignment::printProcessTrackSummary()
{
  LOGF(info, "Alignment processRecoTracks() summary: ");
  LOGF(info,
       "n TFs = %d, used tracks = %d, skipped tracks = %d, local equations failed = %d",
       mNumberTFs, mCounterUsedTracks,
       mCounterSkippedTracks, mCounterLocalEquationFailed);
}

//__________________________________________________________________________
bool Alignment::setLocalDerivative(Int_t index, Double_t value)
{
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  bool success = true;
  if (index < mNumberOfTrackParam) {
    mLocalDerivatives[index] = value;
  } else {
    LOGF(error,
         "Alignment::setLocalDerivative() - index %d >= %d",
         index, mNumberOfTrackParam);
    success = false;
  }
  return success;
}

//__________________________________________________________________________
bool Alignment::setGlobalDerivative(Int_t index, Double_t value)
{
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  bool success = true;
  if (index < mNumberOfGlobalParam) {
    mGlobalDerivatives[index] = value;
  } else {
    LOGF(error,
         "Alignment::setGlobalDerivative() - index %d >= %d",
         index, mNumberOfGlobalParam);
    success = false;
  }
  return success;
}

//__________________________________________________________________________
void Alignment::resetLocalDerivative()
{
  for (int i = 0; i < mNumberOfTrackParam; ++i) {
    mLocalDerivatives[i] = 0.0;
  }
}

//__________________________________________________________________________
void Alignment::resetGlocalDerivative()
{
  for (int i = 0; i < mNumberOfGlobalParam; ++i) {
    mGlobalDerivatives[i] = 0.0;
  }
}

//__________________________________________________________________________
bool Alignment::setLocalEquationX()
{

  if (!mAlignPoint->isAlignPointSet())
    return false;
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  // clean slate for the local equation for this measurement

  resetGlocalDerivative();
  resetLocalDerivative();

  bool success = true;

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeX().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeX().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeX().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeX().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 0, mAlignPoint->globalDerivativeX().dDeltaX());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 1, mAlignPoint->globalDerivativeX().dDeltaY());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 2, mAlignPoint->globalDerivativeX().dDeltaRz());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 3, mAlignPoint->globalDerivativeX().dDeltaZ());

  mMillepede->SetLocalEquation(
    mGlobalDerivatives,
    mLocalDerivatives,
    mAlignPoint->getLocalMeasuredPosition().X(),
    mAlignPoint->getMeasuredPositionSigma().X());

  return success;
}

//__________________________________________________________________________
bool Alignment::setLocalEquationY()
{
  if (!mAlignPoint->isAlignPointSet())
    return false;
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  // clean slate for the local equation for this measurement

  resetGlocalDerivative();
  resetLocalDerivative();

  bool success = true;

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeY().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeY().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeY().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeY().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 0, mAlignPoint->globalDerivativeY().dDeltaX());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 1, mAlignPoint->globalDerivativeY().dDeltaY());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 2, mAlignPoint->globalDerivativeY().dDeltaRz());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 3, mAlignPoint->globalDerivativeY().dDeltaZ());

  mMillepede->SetLocalEquation(
    mGlobalDerivatives,
    mLocalDerivatives,
    mAlignPoint->getLocalMeasuredPosition().Y(),
    mAlignPoint->getMeasuredPositionSigma().Y());

  return success;
}

//__________________________________________________________________________
bool Alignment::setLocalEquationZ()
{

  if (!mAlignPoint->isAlignPointSet())
    return false;
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  // clean slate for the local equation for this measurement

  resetGlocalDerivative();
  resetLocalDerivative();

  bool success = true;

  // local derivatives
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  success &= setLocalDerivative(0, mAlignPoint->localDerivativeZ().dX0());
  success &= setLocalDerivative(1, mAlignPoint->localDerivativeZ().dTx());
  success &= setLocalDerivative(2, mAlignPoint->localDerivativeZ().dY0());
  success &= setLocalDerivative(3, mAlignPoint->localDerivativeZ().dTy());

  // global derivatives
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  Int_t chipId = mAlignPoint->getSensorId();
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 0, mAlignPoint->globalDerivativeZ().dDeltaX());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 1, mAlignPoint->globalDerivativeZ().dDeltaY());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 2, mAlignPoint->globalDerivativeZ().dDeltaRz());
  success &= setGlobalDerivative(chipId + mNDofPerSensor + 3, mAlignPoint->globalDerivativeZ().dDeltaZ());

  mMillepede->SetLocalEquation(
    mGlobalDerivatives,
    mLocalDerivatives,
    mAlignPoint->getLocalMeasuredPosition().Z(),
    mAlignPoint->getMeasuredPositionSigma().Z());

  return success;
}
