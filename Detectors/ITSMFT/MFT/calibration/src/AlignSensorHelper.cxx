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

/// @file AlignSensorHelper.cxx

#include <Rtypes.h>
#include "Framework/Logger.h"
#include "MFTCalibration/AlignSensorHelper.h"

using namespace o2::mft;

ClassImp(o2::mft::AlignSensorHelper);

//__________________________________________________________________________
AlignSensorHelper::AlignSensorHelper(const o2::mft::GeometryTGeo* geom)
  : mChipIndexOnLadder(0),
    mChipIndexInMft(0),
    mLadderInHalfDisk(0),
    mLayer(0),
    mDisk(0),
    mHalf(0),
    mChipUniqueId(0),
    mTranslation(0, 0, 0),
    mRx(0),
    mRy(0),
    mRz(0),
    mSinRx(0),
    mCosRx(0),
    mSinRy(0),
    mCosRy(0),
    mSinRz(0),
    mCosRz(0),
    mIsTransformExtracted(false)
{
  setGeometry(geom);
}

//__________________________________________________________________________
void AlignSensorHelper::setGeometry(const o2::mft::GeometryTGeo* geom)
{
  if (mGeometry == nullptr) {
    mGeometry = geom;
    mGeoSymbolicName = mGeometry->composeSymNameMFT();
  }
}

//__________________________________________________________________________
bool AlignSensorHelper::setSensor(const int chipIndex)
{
  resetSensorTransformInfo();

  if (chipIndex < mNumberOfSensors) {
    o2::itsmft::MFTChipMappingData chipMappingData = (mChipMapping.getChipMappingData())[chipIndex];
    mChipIndexOnLadder = (UShort_t)chipMappingData.chipOnModule;
    mChipIndexInMft = chipMappingData.globalChipSWID;
    mLayer = (UShort_t)chipMappingData.layer;
    mDisk = (UShort_t)chipMappingData.disk;
    mHalf = (UShort_t)chipMappingData.half;
  } else {
    LOG(error) << "AlignSensorHelper::setSensor() - "
               << "chip index " << chipIndex
               << " >= " << mNumberOfSensors;
  }

  setSensorUid(chipIndex);
  setSymName();
  extractSensorTransform();
  return mIsTransformExtracted;
}

//__________________________________________________________________________
void AlignSensorHelper::setSensorUid(const int chipIndex)
{
  if (chipIndex < mNumberOfSensors) {
    mChipUniqueId = o2::base::GeometryManager::getSensID(o2::detectors::DetID::MFT,
                                                         chipIndex);
  } else {
    LOG(error) << "AlignSensorHelper::setSensorUid() - "
               << "chip index " << chipIndex
               << " >= " << mNumberOfSensors;
    mChipUniqueId = o2::base::GeometryManager::getSensID(o2::detectors::DetID::MFT, 0);
  }
}

//__________________________________________________________________________
void AlignSensorHelper::setSymName()
{
  int hf = 0, dk = 0, sr = 0;
  if (mGeometry) {
    mGeometry->getSensorID(mChipIndexInMft, hf, dk, mLadderInHalfDisk, sr);
    bool isIdVerified = true;
    isIdVerified &= (hf == (int)mHalf);
    isIdVerified &= (dk == (int)mDisk);
    isIdVerified &= (sr == (int)mChipIndexOnLadder);
    if (isIdVerified) {
      mGeoSymbolicName = mGeometry->composeSymNameChip(mHalf,
                                                       mDisk,
                                                       mLadderInHalfDisk,
                                                       mChipIndexOnLadder);
    } else {
      LOG(error) << "AlignSensorHelper::setSymName() - mismatch in some index";
    }
  } else {
    LOG(error) << "AlignSensorHelper::setSymName() - nullptr to geometry";
  }
}

//__________________________________________________________________________
void AlignSensorHelper::extractSensorTransform()
{
  if (mIsTransformExtracted)
    return;
  if (mGeometry) {
    mTransform = mGeometry->getMatrixL2G(mChipIndexInMft);

    Double_t* tra = mTransform.GetTranslation();
    mTranslation.SetX(tra[0]);
    mTranslation.SetY(tra[1]);
    mTranslation.SetZ(tra[2]);

    Double_t* rot = mTransform.GetRotationMatrix();
    mRx = std::atan2(-rot[5], rot[8]);
    mRy = std::asin(rot[2]);
    mRz = std::atan2(-rot[1], rot[0]);

    // force the value of some calculations of sin, cos to avoid numerical errors

    // for MFT sensors, Rx = - Pi/2, or + Pi/2
    // mSinRx = std::sin(mRx)
    if (mRx > 0)
      mSinRx = 1.0;
    else
      mSinRx = -1.0;
    mSinRx = std::sin(mRx);
    mCosRx = 0.0; // std::cos(mRx)

    // for MFT sensors, Ry = 0
    mSinRy = 0.0; // std::sin(mRy);
    mCosRy = 1.0; // std::cos(mRy);

    // for MFT sensors, Rz = 0 or Pi
    // but we keep the value as it is
    // because deltaRz is one of the alignment d.o.f.
    mSinRz = std::sin(mRz);
    mCosRz = std::cos(mRz);

    mIsTransformExtracted = true;
  } else {
    resetSensorTransformInfo();
    LOG(error) << "AlignSensorHelper::extractSensorTransform() - nullptr to geometry"
               << std::endl;
  }
}

//__________________________________________________________________________
void AlignSensorHelper::resetSensorTransformInfo()
{
  mIsTransformExtracted = false;

  double rot[9] = {
    0., 0., 0.,
    0., 0., 0.,
    0., 0., 0.};
  double tra[3] = {0., 0., 0.};
  mTransform.SetRotation(rot);
  mTransform.SetTranslation(tra);

  mTranslation.SetX(0.0);
  mTranslation.SetY(0.0);
  mTranslation.SetZ(0.0);

  mRx = 0;
  mRy = 0;
  mRz = 0;

  mSinRx = 0;
  mCosRx = 0;
  mSinRy = 0;
  mCosRy = 0;
  mSinRz = 0;
  mCosRz = 0;
}