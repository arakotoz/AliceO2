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

/// @file AlignPointHelper.cxx

#include <Rtypes.h>
#include "Framework/Logger.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "MFTAlignment/AlignSensorHelper.h"
#include "MFTAlignment/AlignPointHelper.h"
#include "MFTTracking/IOUtils.h"
#include "ITSMFTBase/SegmentationAlpide.h"

using namespace o2::mft;

ClassImp(o2::mft::AlignPointHelper);

//__________________________________________________________________________
AlignPointHelper::AlignPointHelper(const GeometryTGeo* geom)
  : mIsAlignPointSet(false),
    mIsGlobalDerivativeDone(false),
    mIsLocalDerivativeDone(false),
    mIsTrackInitialParamSet(false),
    mGeometry(geom),
    mGlobalRecoPosition(0., 0., 0.),
    mLocalMeasuredPosition(0., 0., 0.),
    mMeasuredPositionSigma(0., 0., 0)
{
  mTrackInitialParam.X0 = 0.;
  mTrackInitialParam.Y0 = 0.;
  mTrackInitialParam.Z0 = 0.;
  mTrackInitialParam.Tx = 0.;
  mTrackInitialParam.Ty = 0.;

  mMeasuredPositionSigma.SetXYZ(
    o2::mft::ioutils::DefClusErrorRow,
    o2::mft::ioutils::DefClusErrorCol,
    o2::itsmft::SegmentationAlpide::SensorLayerThicknessEff);

  mChipHelper = std::make_unique<AlignSensorHelper>(mGeometry);
  LOGF(info, "AlignPointHelper instantiated");
}

//__________________________________________________________________________
void AlignPointHelper::computeLocalDerivatives()
{
  mIsLocalDerivativeDone = false;
  if (mChipHelper == nullptr) {
    LOGF(error,
         "AlignPointHelper::computeLocalDerivatives() - no AlignSensorHelper found !");
    return;
  }
  if (!mIsTrackInitialParamSet) {
    LOGF(error,
         "AlignPointHelper::computeLocalDerivatives() - no initial track param found !");
    return;
  }
  if (!mIsAlignPointSet) {
    LOGF(error,
         "AlignPointHelper::computeLocalDerivatives() - no align point coordinates set !");
    return;
  }
  bool success = true;
  success &= computeLocalDerivativeX();
  success &= computeLocalDerivativeY();
  success &= computeLocalDerivativeZ();
  mIsLocalDerivativeDone = success;
}

//__________________________________________________________________________
void AlignPointHelper::computeGlobalDerivatives()
{
  mIsGlobalDerivativeDone = false;
  if (mChipHelper == nullptr) {
    LOGF(error,
         "AlignPointHelper::computeGlobalDerivatives() - no AlignSensorHelper found !");
    return;
  }
  if (!mIsTrackInitialParamSet) {
    LOGF(error,
         "AlignPointHelper::computeGlobalDerivatives() - no initial track param found !");
    return;
  }
  if (!mIsAlignPointSet) {
    LOGF(error, "AlignPointHelper::computeLocalDerivatives() - no align point coordinates set !");
    return;
  }
  bool success = true;
  success &= computeGlobalDerivativeX();
  success &= computeGlobalDerivativeY();
  success &= computeGlobalDerivativeZ();
  mIsGlobalDerivativeDone = success;
}

//__________________________________________________________________________
Int_t AlignPointHelper::getSensorId() const
{
  if (mChipHelper == nullptr) {
    LOGF(error,
         "AlignPointHelper::getSensorId() - no AlignSensorHelper found !");
    return 0;
  }
  if (!mIsAlignPointSet) {
    LOGF(error,
         "AlignPointHelper::getSensorId() - no align point coordinates set !");
    return 0;
  }
  return mChipHelper->chipIndexInMft();
}

//__________________________________________________________________________
void AlignPointHelper::resetAlignPoint()
{
  mGlobalRecoPosition.SetXYZ(0., 0., 0.);
  mLocalMeasuredPosition.SetXYZ(0., 0., 0.);

  mIsAlignPointSet = false;
}

//__________________________________________________________________________
void AlignPointHelper::resetDerivatives()
{
  mLocalDerivativeX.reset();
  mLocalDerivativeY.reset();
  mLocalDerivativeZ.reset();

  mGlobalDerivativeX.reset();
  mGlobalDerivativeY.reset();
  mGlobalDerivativeZ.reset();

  mIsLocalDerivativeDone = false;
  mIsGlobalDerivativeDone = false;
}

//__________________________________________________________________________
void AlignPointHelper::resetTrackInitialParam()
{
  mTrackInitialParam.X0 = 0.;
  mTrackInitialParam.Y0 = 0.;
  mTrackInitialParam.Z0 = 0.;
  mTrackInitialParam.Tx = 0.;
  mTrackInitialParam.Ty = 0.;

  mIsTrackInitialParamSet = false;
}

//__________________________________________________________________________
void AlignPointHelper::recordTrackInitialParam(o2::mft::TrackMFT mftTrack)
{
  mIsTrackInitialParamSet = false;
  mTrackInitialParam.X0 = mftTrack.getX();
  mTrackInitialParam.Y0 = mftTrack.getY();
  mTrackInitialParam.Z0 = mftTrack.getZ();
  double phi = mftTrack.getPhi();
  double tanLambda = mftTrack.getTanl();
  mTrackInitialParam.Tx = std::cos(phi) / tanLambda;
  mTrackInitialParam.Ty = std::sin(phi) / tanLambda;
  mIsTrackInitialParamSet = true;
}

//__________________________________________________________________________
void AlignPointHelper::setGlobalRecoPosition(o2::mft::TrackMFT mftTrack)
{
  mIsAlignPointSet = false;
  mGlobalRecoPosition.SetXYZ(mftTrack.getX(), mftTrack.getY(), mftTrack.getZ());
  mIsAlignPointSet = true;
}

//__________________________________________________________________________
void AlignPointHelper::setLocalMeasuredPosition(o2::BaseCluster<float> mftCluster)
{
  o2::math_utils::Point3D<double> gloXYZ(
    mftCluster.getX(), mftCluster.getY(), mftCluster.getZ());
  mLocalMeasuredPosition = mGeometry->getMatrixL2G(mftCluster.getSensorID()).ApplyInverse(gloXYZ);
  mIsAlignPointSet &= mChipHelper->setSensor(mftCluster.getSensorID());
}

//__________________________________________________________________________
bool AlignPointHelper::computeLocalDerivativeX()
{
  if (mChipHelper->isTransformExtracted()) {
    mLocalDerivativeX.mdX0 = mChipHelper->cosRy() * mChipHelper->cosRz();

    mLocalDerivativeX.mdTx = (mGlobalRecoPosition.Z() - mTrackInitialParam.Z0) *
                             mChipHelper->cosRy() * mChipHelper->cosRz();

    mLocalDerivativeX.mdY0 = (mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->cosRz()) +
                             (mChipHelper->cosRx() * mChipHelper->sinRz());

    mLocalDerivativeX.mdTy = (mGlobalRecoPosition.Z() - mTrackInitialParam.Z0) *
                             ((mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->cosRz()) +
                              (mChipHelper->cosRx() * mChipHelper->sinRz()));
    return true;
  } else {
    LOGF(error,
         "AlignPointHelper::computeLocalDerivativeX() - no sensor transform found !");
    return false;
  }
}

//__________________________________________________________________________
bool AlignPointHelper::computeLocalDerivativeY()
{
  if (mChipHelper->isTransformExtracted()) {
    mLocalDerivativeY.mdX0 = (-1.) * mChipHelper->cosRy() * mChipHelper->sinRz();

    mLocalDerivativeY.mdTx = (-1.) * (mGlobalRecoPosition.Z() - mTrackInitialParam.Z0) *
                             mChipHelper->cosRy() * mChipHelper->sinRz();

    mLocalDerivativeY.mdY0 = (mChipHelper->cosRx() * mChipHelper->cosRz()) -
                             (mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->sinRz());

    mLocalDerivativeY.mdTy = (mGlobalRecoPosition.Z() - mTrackInitialParam.Z0) *
                             ((mChipHelper->cosRx() * mChipHelper->cosRz()) -
                              (mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->sinRz()));
    return true;
  } else {
    LOGF(error,
         "AlignPointHelper::computeLocalDerivativeY() - no sensor transform found !");
    return false;
  }
}

//__________________________________________________________________________
bool AlignPointHelper::computeLocalDerivativeZ()
{
  if (mChipHelper->isTransformExtracted()) {
    mLocalDerivativeZ.mdX0 = mChipHelper->sinRy();

    mLocalDerivativeZ.mdTx = (mGlobalRecoPosition.Z() - mTrackInitialParam.Z0) * mChipHelper->sinRy();

    mLocalDerivativeZ.mdY0 = (-1.) * mChipHelper->sinRx() * mChipHelper->cosRy();

    mLocalDerivativeZ.mdTy = (-1.) * (mGlobalRecoPosition.Z() - mTrackInitialParam.Z0) * mChipHelper->sinRx() * mChipHelper->cosRy();
    return true;
  } else {
    LOGF(error,
         "AlignPointHelper::computeLocalDerivativeZ() - no sensor transform found !");
    return false;
  }
}

//__________________________________________________________________________
bool AlignPointHelper::computeGlobalDerivativeX()
{
  if (mChipHelper->isTransformExtracted()) {
    mGlobalDerivativeX.mdDeltaX = (-1.) * mChipHelper->cosRy() * mChipHelper->cosRz();

    mGlobalDerivativeX.mdDeltaY = (-1) * ((mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->cosRz()) +
                                          (mChipHelper->cosRx() * mChipHelper->sinRz()));

    mGlobalDerivativeX.mdDeltaZ = (mTrackInitialParam.Tx * mChipHelper->cosRy() * mChipHelper->cosRz()) +
                                  (mTrackInitialParam.Ty * ((mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->cosRz()) +
                                                            (mChipHelper->cosRx() * mChipHelper->sinRz())));

    mGlobalDerivativeX.mdDeltaRz = ((-1.) * mChipHelper->cosRy() * mChipHelper->sinRz() *
                                    (mGlobalRecoPosition.X() - mChipHelper->translateX())) +
                                   ((mGlobalRecoPosition.Y() - mChipHelper->translateY()) *
                                    ((mChipHelper->cosRx() * mChipHelper->cosRz()) -
                                     (mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->sinRz())));

    return true;
  } else {
    LOGF(error,
         "AlignPointHelper::computeGlobalDerivativeX() - no sensor transform found !");
    return false;
  }
}

//__________________________________________________________________________
bool AlignPointHelper::computeGlobalDerivativeY()
{
  if (mChipHelper->isTransformExtracted()) {
    mGlobalDerivativeY.mdDeltaX = mChipHelper->cosRy() * mChipHelper->sinRz();

    mGlobalDerivativeY.mdDeltaY = (mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->sinRz()) -
                                  (mChipHelper->cosRx() * mChipHelper->cosRz());

    mGlobalDerivativeY.mdDeltaZ = ((-1.) * mTrackInitialParam.Tx * mChipHelper->cosRy() * mChipHelper->sinRz()) +
                                  (mTrackInitialParam.Ty *
                                   ((mChipHelper->cosRx() * mChipHelper->cosRz()) -
                                    (mChipHelper->sinRx() * mChipHelper->sinRy() * mChipHelper->sinRz())));

    mGlobalDerivativeY.mdDeltaRz = ((-1.) * (mGlobalRecoPosition.X() - mChipHelper->translateX()) * mChipHelper->cosRy() * mChipHelper->cosRz()) -
                                   ((mGlobalRecoPosition.Y() - mChipHelper->translateY()) *
                                    ((mChipHelper->cosRx() * mChipHelper->sinRz()) +
                                     (mChipHelper->sinRx() * mChipHelper->sinRy() + mChipHelper->cosRz())));
    return true;
  } else {
    LOGF(error,
         "AlignPointHelper::computeGlobalDerivativeY() - no sensor transform found !");
    return false;
  }
}

//__________________________________________________________________________
bool AlignPointHelper::computeGlobalDerivativeZ()
{
  if (mChipHelper->isTransformExtracted()) {
    mGlobalDerivativeZ.mdDeltaX = (-1.) * mChipHelper->sinRy();

    mGlobalDerivativeZ.mdDeltaY = mChipHelper->sinRx() * mChipHelper->cosRy();

    mGlobalDerivativeZ.mdDeltaZ = (mTrackInitialParam.Tx * mChipHelper->sinRy()) -
                                  (mTrackInitialParam.Ty * mChipHelper->sinRx() * mChipHelper->cosRy());

    mGlobalDerivativeZ.mdDeltaRz = 0;

    return true;
  } else {
    LOGF(error,
         "AlignPointHelper::computeGlobalDerivativeZ() - no sensor transform found !");
    return false;
  }
}
