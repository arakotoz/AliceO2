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
  if (mDictionary == nullptr) {
    LOGF(fatal, "Alignment::init() failed because no cluster dictionary is defined");
    mIsInitDone = false;
    return;
  }
  mAlignPoint = std::make_unique<AlignPointHelper>();
  mAlignPoint->setClusterDictionary(mDictionary);
  mMillepede->InitMille(mNumberOfGlobalParam,
                        mNumberOfTrackParam,
                        mChi2CutNStdDev,
                        mResCut,
                        mResCutInitial);
  // filenames for the processed tracks and constraints records
  mMillepede->SetDataRecFName(mMilleRecordsFileName.Data());
  mMillepede->SetConsRecFName(mMilleConstraintsRecFileName.Data());

  bool read = false;
  mMillepede->InitDataRecStorage(read);
  LOG(info) << "-------------- Alignment configured with -----------------";
  LOGF(info, "Chi2CutNStdDev = %d", mChi2CutNStdDev);
  LOGF(info, "ResidualCutInitial = %.3f", mResCutInitial);
  LOGF(info, "ResidualCut = %.3f", mResCut);
  LOGF(info, "MinNumberClusterCut = %d", mMinNumberClusterCut);
  LOGF(info, "mStartFac = %.3f", mStartFac);
  LOGF(info,
       "Allowed variation: dx = %.3f, dy = %.3f, dz = %.3f, dRz = %.4f",
       mAllowVar[0], mAllowVar[1], mAllowVar[3], mAllowVar[2]);
  LOG(info) << "-----------------------------------------------------------";

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

    bool isTrackUsed = true;

    for (int icls = 0; icls < ncls; ++icls) { // cluster loop

      mAlignPoint->resetDerivatives();
      mAlignPoint->resetAlignPoint();

      // Store measured positions
      auto clsEntry = mMFTTrackClusIdx[offset + icls];
      mAlignPoint->setMeasuredPosition(mMFTClusters[clsEntry], pattIt);

      // Propagate track to the current z plane of this cluster
      oneTrack.propagateParamToZlinear(mAlignPoint->getGlobalMeasuredPosition().Z());

      // Store reco positions
      mAlignPoint->setGlobalRecoPosition(oneTrack);

      // compute residuals
      mAlignPoint->setLocalResidual();
      mAlignPoint->setGlobalResidual();

      // Compute derivatives
      mAlignPoint->computeLocalDerivatives();
      mAlignPoint->computeGlobalDerivatives();

      // Set local equations
      bool success = true;
      success &= setLocalEquationX();
      success &= setLocalEquationY();
      success &= setLocalEquationZ();
      isTrackUsed &= success;

    } // end of loop on clusters

    if (isTrackUsed) {
      mMillepede->SetRecordRun(mRunNumber);
      mMillepede->SetRecordWeight(mWeightRecord);
      mTrackRecord = *mMillepede->GetRecord(); // copy track record (A.R. why ?)
      mMillepede->SaveRecordData();            // save record data
      mCounterUsedTracks++;
    }
  } // end of loop on tracks
}

//__________________________________________________________________________
void Alignment::globalFit()
{
  if (!mIsInitDone) {
    LOGF(fatal, "Alignment::globalFit() aborted because init was not done!");
    return;
  }

  if (!mCounterUsedTracks) {
    LOGF(fatal, "Alignment::globalFit() aborted because no reco track was used!");
    return;
  }

  // allocate memory in arrays to temporarily store the results of the global fit

  Double_t* params = new Double_t[mNumberOfGlobalParam];
  Double_t* paramsErrors = new Double_t[mNumberOfGlobalParam];
  Double_t* paramsPulls = new Double_t[mNumberOfGlobalParam];

  // initialise the content of each array

  for (int ii = 0; ii < mNumberOfGlobalParam; ii++) {
    params[ii] = 0.;
    paramsErrors[ii] = 0.;
    paramsPulls[ii] = 0.;
  }

  // perform the simultaneous fit of track and alignement parameters

  mMillepede->GlobalFit(params, paramsErrors, paramsPulls);

  // post-treatment:
  // debug output + save Millepede global fit result in AlignParam vector

  LOGF(info, "Alignment: done fitting global parameters");
  LOGF(info, "sensor info, dx (cm), dy (cm), dz (cm), dRz (rad)");

  AlignSensorHelper chipHelper;
  double dRx = 0., dRy = 0., dRz = 0.; // delta rotations
  double dx = 0., dy = 0., dz = 0.;    // delta translations
  bool global = true;                  // delta in global ref. system
  bool withSymName = false;

  for (int chipId = 0; chipId < mNumberOfSensors; chipId++) {
    chipHelper.setSensorOnlyInfo(chipId);
    std::stringstream name = chipHelper.getSensorFullName(withSymName);
    dx = params[chipId * mNDofPerSensor + 0];
    dy = params[chipId * mNDofPerSensor + 1];
    dz = params[chipId * mNDofPerSensor + 3];
    dRz = params[chipId * mNDofPerSensor + 2];
    LOGF(info,
         "%s, %.3e, %.3e, %.3e, %.3e",
         name.str().c_str(), dx, dy, dz, dRz);
    mAlignParams.emplace_back(
      chipHelper.geoSymbolicName(),
      chipHelper.sensorUid(),
      dx, dy, dz,
      dRx, dRy, dRz,
      global);
  }

  // free allocated memory

  delete[] params;
  delete[] paramsErrors;
  delete[] paramsPulls;
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
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeX().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeX().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeX().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeX().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5)
      LOGF(debug,
           "Alignment::setLocalEquationX(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e X %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalResidual().X(),
           mAlignPoint->getLocalMeasuredPositionSigma().X());

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalResidual().X(),
      mAlignPoint->getLocalMeasuredPositionSigma().X());
  } else {
    mCounterLocalEquationFailed++;
  }

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
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeY().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeY().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeY().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeY().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5)
      LOGF(debug,
           "Alignment::setLocalEquationY(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Y %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalResidual().Y(),
           mAlignPoint->getLocalMeasuredPositionSigma().Y());

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalResidual().Y(),
      mAlignPoint->getLocalMeasuredPositionSigma().Y());
  } else {
    mCounterLocalEquationFailed++;
  }

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
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 0, mAlignPoint->globalDerivativeZ().dDeltaX());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 1, mAlignPoint->globalDerivativeZ().dDeltaY());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 2, mAlignPoint->globalDerivativeZ().dDeltaRz());
  success &= setGlobalDerivative(chipId * mNDofPerSensor + 3, mAlignPoint->globalDerivativeZ().dDeltaZ());

  if (success) {
    if (mCounterUsedTracks < 5)
      LOGF(info,
           "Alignment::setLocalEquationZ(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Z %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalResidual().Z(),
           mAlignPoint->getLocalMeasuredPositionSigma().Z());

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalResidual().Z(),
      mAlignPoint->getLocalMeasuredPositionSigma().Z());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}
