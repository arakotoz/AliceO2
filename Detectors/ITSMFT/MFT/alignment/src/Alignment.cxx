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
    mNumberOfClusterChainROFs(0),
    mNumberOfTrackChainROFs(0),
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
    mMillepede(nullptr),
    mDictionary(nullptr),
    mAlignPoint(nullptr),
    mIsInitDone(false),
    mGlobalParameterStatus(nullptr),
    mWithControl(false),
    mNEntriesAutoSave(10000),
    mWithRecordWriter(true),
    mRecordWriter(nullptr),
    mWithConstraintsRecWriter(false),
    mConstraintsRecWriter(nullptr),
    mWithRecordReader(false),
    mRecordReader(nullptr),
    mWithConstraintsRecReader(false),
    mConstraintsRecReader(nullptr)
{
  // default allowed variations w.r.t. global system coordinates
  mAllowVar[0] = 0.5;  // delta translation in x (cm)
  mAllowVar[1] = 0.5;  // delta translation in y (cm)
  mAllowVar[2] = 0.01; // rotation angle Rz around z-axis (rad)
  mAllowVar[3] = 0.5;  // delta translation in z (cm)

  // allocate memory for local and global derivatives
  mGlobalDerivatives = (double*)malloc(sizeof(double) * mNumberOfGlobalParam);
  mLocalDerivatives = new double[mNumberOfTrackParam];
  mGlobalParameterStatus = (int*)malloc(sizeof(int) * mNumberOfGlobalParam);

  // initialise the content of each array
  resetGlocalDerivative();
  resetLocalDerivative();
  for (int iPar = 0; iPar < mNumberOfGlobalParam; iPar++) {
    mGlobalParameterStatus[iPar] = mFreeParId;
  }
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
  if (mWithRecordWriter) {
    mRecordWriter = std::make_shared<MilleRecordWriter>();
    mRecordWriter->setCyclicAutoSave(mNEntriesAutoSave);
    mRecordWriter->setDataFileName(mMilleRecordsFileName);
    mMillepede->SetRecordWriter(mRecordWriter);
  }
  if (mWithConstraintsRecWriter) {
    mConstraintsRecWriter = std::make_shared<MilleRecordWriter>();
    mConstraintsRecWriter->setCyclicAutoSave(mNEntriesAutoSave);
    mConstraintsRecWriter->setDataFileName(mMilleConstraintsRecFileName);
    mMillepede->SetConstraintsRecWriter(mConstraintsRecWriter);
  }
  if (mWithRecordReader) {
    mRecordReader = std::make_shared<MilleRecordReader>();
    mMillepede->SetRecordReader(mRecordReader);
  }
  if (mWithConstraintsRecReader) {
    mConstraintsRecReader = std::make_shared<MilleRecordReader>();
    mMillepede->SetConstraintsRecReader(mConstraintsRecReader);
  }
  mAlignPoint = std::make_shared<AlignPointHelper>();
  mAlignPoint->setClusterDictionary(mDictionary);
  mMillepede->InitMille(mNumberOfGlobalParam,
                        mNumberOfTrackParam,
                        mChi2CutNStdDev,
                        mResCut,
                        mResCutInitial);

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
    LOGF(fatal, "Alignment::processRecoTracks() aborted because init was not done !");
    return;
  }
  if (!mWithRecordWriter || !mRecordWriter || !mRecordWriter->isInitOk()) {
    LOGF(fatal, "Alignment::processRecoTracks() aborted because uninitialised mRecordWriter !");
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

    LOGF(debug, "Processing track # %5d", nCounterAllTracks);

    auto offset = oneTrack.getExternalClusterIndexOffset();

    mRecordWriter->getRecord()->Reset();

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
        LOGF(error, "Alignment::processRecoTracks() - track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
             mCounterUsedTracks, mAlignPoint->half(), mAlignPoint->disk(), mAlignPoint->layer(), mAlignPoint->getSensorId(),
             mAlignPoint->getLocalMeasuredPosition().X(), mAlignPoint->getLocalMeasuredPosition().Y(), mAlignPoint->getLocalMeasuredPosition().Z(),
             mAlignPoint->getGlobalMeasuredPosition().X(), mAlignPoint->getGlobalMeasuredPosition().Y(), mAlignPoint->getGlobalMeasuredPosition().Z(),
             mAlignPoint->getGlobalRecoPosition().X(), mAlignPoint->getGlobalRecoPosition().Y(), mAlignPoint->getGlobalRecoPosition().Z());
      }
      isTrackUsed &= success;

    } // end of loop on clusters

    if (isTrackUsed) {
      mRecordWriter->setRecordRun(mRunNumber);
      mRecordWriter->setRecordWeight(mWeightRecord);
      mRecordWriter->fillRecordTree(); // save record data
      mCounterUsedTracks++;
    }
  } // end of loop on tracks
}

//__________________________________________________________________________
void Alignment::processROFs(TChain* mfttrackChain, TChain* mftclusterChain)
{
  if (!mIsInitDone) {
    LOGF(fatal, "Alignment::processROFs() aborted because init was not done !");
    return;
  }

  if (!mWithRecordWriter || !mRecordWriter || !mRecordWriter->isInitOk()) {
    LOGF(fatal, "Alignment::processROFs() aborted because uninitialised mRecordWriter !");
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
  int nCounterAllTracks = 0;
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

      LOGF(info, "Processing track # %5d", nCounterAllTracks);

      auto offset = oneTrack.getExternalClusterIndexOffset();

      mRecordWriter->getRecord()->Reset();

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
          mPointControl.fill(mAlignPoint, mCounterUsedTracks);
        isTrackUsed &= success;
        if (!success) {
          LOGF(error, "Alignment::processROFs() - track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
               mCounterUsedTracks, mAlignPoint->half(), mAlignPoint->disk(), mAlignPoint->layer(), mAlignPoint->getSensorId(),
               mAlignPoint->getLocalMeasuredPosition().X(), mAlignPoint->getLocalMeasuredPosition().Y(), mAlignPoint->getLocalMeasuredPosition().Z(),
               mAlignPoint->getGlobalMeasuredPosition().X(), mAlignPoint->getGlobalMeasuredPosition().Y(), mAlignPoint->getGlobalMeasuredPosition().Z(),
               mAlignPoint->getGlobalRecoPosition().X(), mAlignPoint->getGlobalRecoPosition().Y(), mAlignPoint->getGlobalRecoPosition().Z());
        }

      } // end of loop on clusters

      if (isTrackUsed) {
        // copy track record
        mRecordWriter->setRecordRun(mRunNumber);
        mRecordWriter->setRecordWeight(mWeightRecord);
        const bool doPrint = false;
        mRecordWriter->fillRecordTree(doPrint); // save record data
        mCounterUsedTracks++;
      }
      nCounterAllTracks++;
    } // end of loop on tracks

    firstEntry = false;

  } // end of loop on TChain reader

  LOG(info) << "Alignment::processROFs() - end";
}

//__________________________________________________________________________
void Alignment::globalFit()
{
  if (!mIsInitDone) {
    LOGF(fatal, "Alignment::globalFit() aborted because init was not done !");
    return;
  }
  if (!mWithRecordReader || !mRecordReader ||
      !mRecordReader->isReaderOk() || !mRecordReader->getNEntries()) {
    LOGF(fatal, "Alignment::globalFit() aborted because no data record can be read !");
    return;
  }

  // initialize the file and tree to store chi2 from Millepede LocalFit()

  if (mWithControl)
    mMillepede->InitChi2Storage(mNEntriesAutoSave);

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

  if (mWithControl)
    mMillepede->CloseChi2Storage();

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
  if (!mWithRecordWriter)
    return;

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
void Alignment::startRecordWriter()
{
  if (!mWithRecordWriter)
    return;

  if (mRecordWriter)
    mRecordWriter->init();
  if (mWithControl) {
    mPointControl.setCyclicAutoSave(mNEntriesAutoSave);
    mPointControl.init();
  }
}

//__________________________________________________________________________
void Alignment::endRecordWriter()
{
  if (!mWithRecordWriter)
    return;

  if (mRecordWriter) {
    mRecordWriter->terminate(); // write record tree and close output file
  }
  if (mWithControl)
    mPointControl.terminate();
}

//__________________________________________________________________________
void Alignment::startConstraintsRecWriter()
{
  if (!mWithConstraintsRecWriter)
    return;

  if (mConstraintsRecWriter) {
    mConstraintsRecWriter->changeDataBranchName();
    mConstraintsRecWriter->init();
  }
}

//__________________________________________________________________________
void Alignment::endConstraintsRecWriter()
{
  if (!mWithConstraintsRecWriter)
    return;

  if (mConstraintsRecWriter) {
    mConstraintsRecWriter->terminate();
  }
}

//__________________________________________________________________________
void Alignment::connectRecordReaderToChain(TChain* ch)
{
  if (!mWithRecordReader)
    return;

  if (mRecordReader) {
    mRecordReader->connectToChain(ch);
  }
}

//__________________________________________________________________________
void Alignment::connectConstraintsRecReaderToChain(TChain* ch)
{
  if (!mWithConstraintsRecReader)
    return;

  if (mConstraintsRecReader) {
    mConstraintsRecReader->changeDataBranchName();
    mConstraintsRecReader->connectToChain(ch);
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
    if (mCounterUsedTracks < 5) {
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
    }
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
    if (mCounterUsedTracks < 5) {
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
    }
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
    if (mCounterUsedTracks < 5) {
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
    }
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
