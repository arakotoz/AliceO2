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

//_________________________________________________________
o2::math_utils::Point3D<double> AlignSensorHelper::getSensorCenterGlobalCoordinates(const int chipId) const
{
  // init the value of the global coordinates
  o2::math_utils::Point3D<double> gloXYZ(0, 0, 0);

  if (chipId < mNumberOfSensors) {
    // The center of the sensor is the origin of the local reference system
    o2::math_utils::Point3D<double> locXYZ(0, 0, 0);

    // Transformation local --> global coordinates
    gloXYZ = mGeometry->getMatrixL2G(chipId) * locXYZ;
  } else {
    LOG(warning) << "AlignSensorHelper::getSensorCenterGlobalCoordinates()"
                 << " sensor id " << chipId
                 << " >= " << mNumberOfSensors;
  }
  return gloXYZ;
}

//__________________________________________________________________________
int AlignSensorHelper::half(const int chipId) const
{
  UShort_t half = 0;
  if (chipId < mNumberOfSensors) {
    o2::itsmft::MFTChipMappingData chipMapping = (mChipMapping.getChipMappingData())[chipId];
    half = (UShort_t)chipMapping.half;
  } else {
    LOG(warning) << "AlignSensorHelper::half()"
                 << " sensor id " << chipId
                 << " >= " << mNumberOfSensors;
  }
  return (int)half;
}

//__________________________________________________________________________
int AlignSensorHelper::disk(const int chipId) const
{
  UShort_t disk = 0;
  if (chipId < mNumberOfSensors) {
    o2::itsmft::MFTChipMappingData chipMapping = (mChipMapping.getChipMappingData())[chipId];
    disk = (UShort_t)chipMapping.disk;
  } else {
    LOG(warning) << "AlignSensorHelper::disk()"
                 << " sensor id " << chipId
                 << " >= " << mNumberOfSensors;
  }
  return (int)disk;
}

//__________________________________________________________________________
int AlignSensorHelper::layer(const int chipId) const
{
  Int_t sensor = (Int_t)chipId;
  UShort_t layer = 0;
  if (sensor < mNumberOfSensors) {
    o2::itsmft::MFTChipMappingData chipMapping = (mChipMapping.getChipMappingData())[sensor];
    layer = (UShort_t)chipMapping.layer;
  } else {
    LOG(warning) << "AlignSensorHelper::layer()"
                 << " sensor id " << chipId
                 << " >= " << mNumberOfSensors;
  }
  return (int)layer;
}