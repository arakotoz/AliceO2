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
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

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
    mGlobalDerivatives(nullptr),
    mLocalDerivatives(nullptr),
    mStartFac(256),
    mChi2CutNStdDev(3),
    mResCutInitial(100.),
    mResCut(100.),
    mMinNumberClusterCut(6),
    mWeightRecord(1.),
    mMilleRecordsFileName("mft_mille_records.root"),
    mMilleConstraintsRecFileName("mft_mille_constraints.root"),
    mIsInitDone(false),
    mGlobalParameterStatus(nullptr),
    mWithControl(false),
    mControlFile(nullptr),
    mControlTree(nullptr)
{
  // default allowed variations w.r.t. global system coordinates
  mAllowVar[0] = 0.5;  // delta translation in x (cm)
  mAllowVar[1] = 0.5;  // delta translation in y (cm)
  mAllowVar[2] = 0.01; // rotation angle Rz around z-axis (rad)
  mAllowVar[3] = 0.5;  // delta translation in z (cm)

  mGlobalDerivatives = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);
  mLocalDerivatives = new Double_t[mNumberOfTrackParam];

  // initialise the content of each array
  resetGlocalDerivative();
  resetLocalDerivative();

  mGlobalParameterStatus = (int*)malloc(sizeof(int) * mNumberOfGlobalParam);
  for (int iPar = 0; iPar < mNumberOfGlobalParam; iPar++) {
    mGlobalParameterStatus[iPar] = mFreeParId;
  }
  mPointInfo.sensor = 0;
  mPointInfo.layer = 0;
  mPointInfo.disk = 0;
  mPointInfo.half = 0;
  mPointInfo.measuredGlobalX = 0;
  mPointInfo.measuredGlobalY = 0;
  mPointInfo.measuredGlobalZ = 0;
  mPointInfo.measuredLocalX = 0;
  mPointInfo.measuredLocalY = 0;
  mPointInfo.measuredLocalZ = 0;
  mPointInfo.residualX = 0;
  mPointInfo.residualY = 0;
  mPointInfo.residualZ = 0;
  mPointInfo.residualLocalX = 0;
  mPointInfo.residualLocalY = 0;
  mPointInfo.residualLocalZ = 0;
  mPointInfo.recoGlobalX = 0;
  mPointInfo.recoGlobalY = 0;
  mPointInfo.recoGlobalZ = 0;
  mPointInfo.recoLocalX = 0;
  mPointInfo.recoLocalY = 0;
  mPointInfo.recoLocalZ = 0;
  LOGF(info, "Alignment instantiated");
}

//__________________________________________________________________________
Alignment::~Alignment()
{
  free(mGlobalDerivatives);
  delete[] mLocalDerivatives;
  free(mGlobalParameterStatus);
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
  mMillepede = std::make_unique<MillePede2>();
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

  bool read = true;
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

  // init tree to record local measurements and residuals
  if (mWithControl)
    initControlTree();

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

  for (auto& oneTrack : mMFTTracks) { // track loop

    // Skip the track if not enough clusters
    auto ncls = oneTrack.getNumberOfPoints();
    if (ncls < mMinNumberClusterCut) {
      mCounterSkippedTracks++;
      continue;
    }

    // Skip presumably quite low momentum track
    if (!oneTrack.isLTF()) {
      mCounterSkippedTracks++;
      continue;
    }

    auto offset = oneTrack.getExternalClusterIndexOffset();

    if (mMillepede->GetRecord()) {
      mMillepede->GetRecord()->Reset();
    }

    // Store the initial track parameters
    auto track = oneTrack;
    mAlignPoint->resetTrackInitialParam();
    mAlignPoint->recordTrackInitialParam(track);

    bool isTrackUsed = true;

    for (int icls = 0; icls < ncls; ++icls) { // cluster loop

      mAlignPoint->resetAlignPoint();

      // Store measured positions
      auto clsEntry = mMFTTrackClusIdx[offset + icls];
      const auto compCluster = mMFTClusters[clsEntry];
      mAlignPoint->setMeasuredPosition(compCluster, pattIt);

      // Propagate track to the current z plane of this cluster
      track.propagateParamToZlinear(mAlignPoint->getGlobalMeasuredPosition().Z());

      // Store reco positions
      mAlignPoint->setGlobalRecoPosition(track);

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
      if (!success) {
        LOGF(error, "track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
             mCounterUsedTracks, mAlignPoint->half(), mAlignPoint->disk(), mAlignPoint->layer(), mAlignPoint->getSensorId(),
             mAlignPoint->getLocalMeasuredPosition().X(), mAlignPoint->getLocalMeasuredPosition().Y(), mAlignPoint->getLocalMeasuredPosition().Z(),
             mAlignPoint->getGlobalMeasuredPosition().X(), mAlignPoint->getGlobalMeasuredPosition().Y(), mAlignPoint->getGlobalMeasuredPosition().Z(),
             mAlignPoint->getGlobalRecoPosition().X(), mAlignPoint->getGlobalRecoPosition().Y(), mAlignPoint->getGlobalRecoPosition().Z());
      }
      isTrackUsed &= success;

    } // end of loop on clusters

    if (isTrackUsed) {
      mMillepede->SetRecordRun(mRunNumber);
      mMillepede->SetRecordWeight(mWeightRecord);
      mMillepede->SaveRecordData(); // save record data
      mCounterUsedTracks++;
    }
  } // end of loop on tracks
}

//__________________________________________________________________________
void Alignment::processROFs(TChain* mfttrackChain, TChain* mftclusterChain)
{
  if (!mIsInitDone) {
    LOGF(fatal, "Alignment::processROFs() aborted because init was not done!");
    return;
  }

  LOG(info) << "Alignment::processROFs() - start";

  TTreeReader mftTrackChainReader(mfttrackChain);
  TTreeReader mftClusterChainReader(mftclusterChain);
  std::vector<unsigned char>::iterator pattIterator;

  TTreeReaderValue<std::vector<o2::mft::TrackMFT>> mftTracks =
    {mftTrackChainReader, "MFTTrack"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mftTracksROF =
    {mftTrackChainReader, "MFTTracksROF"};
  TTreeReaderValue<std::vector<int>> mftTrackClusIdx =
    {mftTrackChainReader, "MFTTrackClusIdx"};

  TTreeReaderValue<std::vector<o2::itsmft::CompClusterExt>> mftClusters =
    {mftClusterChainReader, "MFTClusterComp"};
  TTreeReaderValue<std::vector<o2::itsmft::ROFRecord>> mftClustersROF =
    {mftClusterChainReader, "MFTClustersROF"};
  TTreeReaderValue<std::vector<unsigned char>> mftClusterPatterns =
    {mftClusterChainReader, "MFTClusterPatt"};

  bool firstEntry = true;
  while (mftTrackChainReader.Next() && mftClusterChainReader.Next()) {

    if (firstEntry)
      pattIterator = (*mftClusterPatterns).begin();

    mNumberOfTrackChainROFs += (*mftTracksROF).size();
    mNumberOfClusterChainROFs += (*mftClustersROF).size();
    assert(mNumberOfTrackChainROFs == mNumberOfClusterChainROFs);

    //______________________________________________________
    for (auto& oneTrack : *mftTracks) { // track loop

      // Skip the track if not enough clusters
      auto ncls = oneTrack.getNumberOfPoints();
      if (ncls < mMinNumberClusterCut) {
        mCounterSkippedTracks++;
        continue;
      }

      // Skip presumably quite low momentum track
      if (!oneTrack.isLTF()) {
        mCounterSkippedTracks++;
        continue;
      }

      auto offset = oneTrack.getExternalClusterIndexOffset();

      if (mMillepede->GetRecord()) {
        mMillepede->GetRecord()->Reset();
      }

      // Store the initial track parameters
      mAlignPoint->resetTrackInitialParam();
      mAlignPoint->recordTrackInitialParam(oneTrack);

      bool isTrackUsed = true;

      for (int icls = 0; icls < ncls; ++icls) { // cluster loop

        mAlignPoint->resetAlignPoint();

        // Store measured positions
        auto clsEntry = (*mftTrackClusIdx)[offset + icls];
        const auto compCluster = (*mftClusters)[clsEntry];
        mAlignPoint->setMeasuredPosition(compCluster, pattIterator);

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
        if (mWithControl && success)
          fillControlTree();
        isTrackUsed &= success;
        if (!success) {
          LOGF(error, "track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
               mCounterUsedTracks, mAlignPoint->half(), mAlignPoint->disk(), mAlignPoint->layer(), mAlignPoint->getSensorId(),
               mAlignPoint->getLocalMeasuredPosition().X(), mAlignPoint->getLocalMeasuredPosition().Y(), mAlignPoint->getLocalMeasuredPosition().Z(),
               mAlignPoint->getGlobalMeasuredPosition().X(), mAlignPoint->getGlobalMeasuredPosition().Y(), mAlignPoint->getGlobalMeasuredPosition().Z(),
               mAlignPoint->getGlobalRecoPosition().X(), mAlignPoint->getGlobalRecoPosition().Y(), mAlignPoint->getGlobalRecoPosition().Z());
        }

      } // end of loop on clusters

      if (isTrackUsed) {
        // copy track record
        mMillepede->SetRecordRun(mRunNumber);
        mMillepede->SetRecordWeight(mWeightRecord);
        mMillepede->SaveRecordData();
        mCounterUsedTracks++;
      }
    } // end of loop on tracks

    firstEntry = false;

  } // end of loop on TChain reader

  LOG(info) << "Alignment::processROFs() - end";
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

  double* params = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);
  double* paramsErrors = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);
  double* paramsPulls = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);

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

  LOGF(info, "Alignment::globalFit() - done, results below");
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

  free(params);
  free(paramsErrors);
  free(paramsPulls);
}

//__________________________________________________________________________
void Alignment::printProcessTrackSummary()
{
  LOGF(info, "Alignment processRecoTracks() summary: ");
  if (mNumberOfTrackChainROFs) {
    LOGF(info,
         "n ROFs = %d, used tracks = %d, skipped tracks = %d, local equations failed = %d",
         mNumberOfTrackChainROFs, mCounterUsedTracks,
         mCounterSkippedTracks, mCounterLocalEquationFailed);
  } else {
    LOGF(info,
         "n TFs = %d, used tracks = %d, skipped tracks = %d, local equations failed = %d",
         mNumberTFs, mCounterUsedTracks,
         mCounterSkippedTracks, mCounterLocalEquationFailed);
  }
}

//__________________________________________________________________________
bool Alignment::setLocalDerivative(Int_t index, Double_t value)
{
  // index [0 .. 3] for {dX0, dTx, dY0, dTz}

  bool success = false;
  if (index < mNumberOfTrackParam) {
    mLocalDerivatives[index] = value;
    success = true;
  } else {
    LOGF(error,
         "AlignHelper::setLocalDerivative() - index %d >= %d",
         index, mNumberOfTrackParam);
  }
  return success;
}

//__________________________________________________________________________
bool Alignment::setGlobalDerivative(Int_t index, Double_t value)
{
  // index [0 .. 3] for {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}

  bool success = false;
  if (index < mNumberOfGlobalParam) {
    mGlobalDerivatives[index] = value;
    success = true;
  } else {
    LOGF(error,
         "AlignHelper::setGlobalDerivative() - index %d >= %d",
         index, mNumberOfGlobalParam);
  }
  return success;
}

//__________________________________________________________________________
bool Alignment::resetLocalDerivative()
{
  bool success = false;
  for (int i = 0; i < mNumberOfTrackParam; ++i) {
    success = false;
    mLocalDerivatives[i] = 0.0;
    success = true;
  }
  return success;
}

//__________________________________________________________________________
bool Alignment::resetGlocalDerivative()
{
  bool success = false;
  for (int i = 0; i < mNumberOfGlobalParam; ++i) {
    success = false;
    mGlobalDerivatives[i] = 0.0;
    success = true;
  }
  return success;
}

//__________________________________________________________________________
bool Alignment::setLocalEquationX()
{

  if (!mAlignPoint->isAlignPointSet()) {
    LOGF(error,
         "Alignment::setLocalEquationX() - no align point coordinates set !");
    return false;
  }
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

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
           mAlignPoint->getLocalMeasuredPosition().X(),
           mAlignPoint->getLocalMeasuredPositionSigma().X());

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().X(),
      mAlignPoint->getLocalMeasuredPositionSigma().X());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}

//__________________________________________________________________________
bool Alignment::setLocalEquationY()
{
  if (!mAlignPoint->isAlignPointSet()) {
    LOGF(error,
         "Alignment::setLocalEquationY() - no align point coordinates set !");
    return false;
  }
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

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
           mAlignPoint->getLocalMeasuredPosition().Y(),
           mAlignPoint->getLocalMeasuredPositionSigma().Y());

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().Y(),
      mAlignPoint->getLocalMeasuredPositionSigma().Y());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}

//__________________________________________________________________________
bool Alignment::setLocalEquationZ()
{
  if (!mAlignPoint->isAlignPointSet()) {
    LOGF(error,
         "Alignment::setLocalEquationZ() - no align point coordinates set !");
    return false;
  }
  if (!mAlignPoint->isGlobalDerivativeDone())
    return false;
  if (!mAlignPoint->isLocalDerivativeDone())
    return false;

  bool success = true;

  // clean slate for the local equation for this measurement

  success &= resetGlocalDerivative();
  success &= resetLocalDerivative();

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
      LOGF(debug,
           "Alignment::setLocalEquationZ(): track %i sr %4d local %.3e %.3e %.3e %.3e, global %.3e %.3e %.3e %.3e Z %.3e sigma %.3e",
           mCounterUsedTracks, chipId,
           mLocalDerivatives[0], mLocalDerivatives[1], mLocalDerivatives[2], mLocalDerivatives[3],
           mGlobalDerivatives[chipId * mNDofPerSensor + 0],
           mGlobalDerivatives[chipId * mNDofPerSensor + 1],
           mGlobalDerivatives[chipId * mNDofPerSensor + 2],
           mGlobalDerivatives[chipId * mNDofPerSensor + 3],
           mAlignPoint->getLocalMeasuredPosition().Z(),
           mAlignPoint->getLocalMeasuredPositionSigma().Z());

    mMillepede->SetLocalEquation(
      mGlobalDerivatives,
      mLocalDerivatives,
      mAlignPoint->getLocalMeasuredPosition().Z(),
      mAlignPoint->getLocalMeasuredPositionSigma().Z());
  } else {
    mCounterLocalEquationFailed++;
  }

  return success;
}

//__________________________________________________________________________
void Alignment::initControlTree()
{
  if (mControlFile == nullptr)
    mControlFile = TFile::Open("align_point.root", "recreate", "", 505);

  if (mControlTree == nullptr) {
    mControlTree = new TTree("point", "the align point info tree");
    mControlTree->Branch("sensor", &mPointInfo.sensor, "sensor/s");
    mControlTree->Branch("layer", &mPointInfo.layer, "layer/s");
    mControlTree->Branch("disk", &mPointInfo.disk, "disk/s");
    mControlTree->Branch("half", &mPointInfo.half, "half/s");
    mControlTree->Branch("measuredGlobalX", &mPointInfo.measuredGlobalX, "measuredGlobalX/D");
    mControlTree->Branch("measuredGlobalY", &mPointInfo.measuredGlobalY, "measuredGlobalY/D");
    mControlTree->Branch("measuredGlobalZ", &mPointInfo.measuredGlobalZ, "measuredGlobalZ/D");
    mControlTree->Branch("measuredLocalX", &mPointInfo.measuredLocalX, "measuredLocalX/D");
    mControlTree->Branch("measuredLocalY", &mPointInfo.measuredLocalY, "measuredLocalY/D");
    mControlTree->Branch("measuredLocalZ", &mPointInfo.measuredLocalZ, "measuredLocalZ/D");
    mControlTree->Branch("residualX", &mPointInfo.residualX, "residualX/D");
    mControlTree->Branch("residualY", &mPointInfo.residualY, "residualY/D");
    mControlTree->Branch("residualZ", &mPointInfo.residualZ, "residualZ/D");
    mControlTree->Branch("residualLocalX", &mPointInfo.residualLocalX, "residualLocalX/D");
    mControlTree->Branch("residualLocalY", &mPointInfo.residualLocalY, "residualLocalY/D");
    mControlTree->Branch("residualLocalZ", &mPointInfo.residualLocalZ, "residualLocalZ/D");
    mControlTree->Branch("recoGlobalX", &mPointInfo.recoGlobalX, "recoGlobalX/D");
    mControlTree->Branch("recoGlobalY", &mPointInfo.recoGlobalY, "recoGlobalY/D");
    mControlTree->Branch("recoGlobalZ", &mPointInfo.recoGlobalZ, "recoGlobalZ/D");
    mControlTree->Branch("recoLocalX", &mPointInfo.recoLocalX, "recoLocalX/D");
    mControlTree->Branch("recoLocalY", &mPointInfo.recoLocalY, "recoLocalY/D");
    mControlTree->Branch("recoLocalZ", &mPointInfo.recoLocalZ, "recoLocalZ/D");
  }
}

//__________________________________________________________________________
void Alignment::closeControlTree()
{
  if (mControlTree) {
    if (mControlFile && mControlFile->IsWritable()) {
      mControlFile->cd();
      mControlTree->Write();
    }
    delete mControlTree;
    if (mControlFile) {
      mControlFile->Close();
      delete mControlFile;
    }
  }
  LOG(info) << "Alignment - Closed file align_point.root";
}

//__________________________________________________________________________
void Alignment::fillControlTree()
{
  if (mControlTree) {
    mPointInfo.sensor = mAlignPoint->getSensorId();
    mPointInfo.layer = mAlignPoint->layer();
    mPointInfo.disk = mAlignPoint->disk();
    mPointInfo.half = mAlignPoint->half();
    mPointInfo.measuredGlobalX = mAlignPoint->getGlobalMeasuredPosition().X();
    mPointInfo.measuredGlobalY = mAlignPoint->getGlobalMeasuredPosition().Y();
    mPointInfo.measuredGlobalZ = mAlignPoint->getGlobalMeasuredPosition().Z();
    mPointInfo.measuredLocalX = mAlignPoint->getLocalMeasuredPosition().X();
    mPointInfo.measuredLocalY = mAlignPoint->getLocalMeasuredPosition().Y();
    mPointInfo.measuredLocalZ = mAlignPoint->getLocalMeasuredPosition().Z();
    mPointInfo.residualX = mAlignPoint->getGlobalResidual().X();
    mPointInfo.residualY = mAlignPoint->getGlobalResidual().Y();
    mPointInfo.residualZ = mAlignPoint->getGlobalResidual().Z();
    mPointInfo.residualLocalX = mAlignPoint->getLocalResidual().X();
    mPointInfo.residualLocalY = mAlignPoint->getLocalResidual().Y();
    mPointInfo.residualLocalZ = mAlignPoint->getLocalResidual().Z();
    mPointInfo.recoGlobalX = mAlignPoint->getGlobalRecoPosition().X();
    mPointInfo.recoGlobalY = mAlignPoint->getGlobalRecoPosition().Y();
    mPointInfo.recoGlobalZ = mAlignPoint->getGlobalRecoPosition().Z();
    mPointInfo.recoLocalX = mAlignPoint->getLocalRecoPosition().X();
    mPointInfo.recoLocalY = mAlignPoint->getLocalRecoPosition().Y();
    mPointInfo.recoLocalZ = mAlignPoint->getLocalRecoPosition().Z();

    if (mCounterUsedTracks < 7000) {
      LOGF(debug, "track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e",
           mCounterUsedTracks, mPointInfo.half, mPointInfo.disk, mPointInfo.layer, mPointInfo.sensor,
           mPointInfo.measuredLocalX, mPointInfo.measuredLocalY, mPointInfo.measuredLocalZ,
           mPointInfo.measuredGlobalX, mPointInfo.measuredGlobalY, mPointInfo.measuredGlobalZ);
    }
    mControlTree->Fill();
  }
}