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

/// \file AlignSensorHelper.h
/// \author arakotoz@cern.ch
/// \brief Helper class to access to the global coordinates of the center each MFT sensor

#ifndef ALICEO2_MFT_ALIGN_SENSOR_HELPER_H
#define ALICEO2_MFT_ALIGN_SENSOR_HELPER_H

#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MathUtils/Cartesian.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

namespace o2
{
namespace mft
{

/// \class AlignSensorHelper
class AlignSensorHelper
{
 public:
  /// \brief default constructor
  AlignSensorHelper() = default;

  /// \brief default destructor
  ~AlignSensorHelper() = default;

  /// \brief set pointer to geometry that should already have done fillMatrixCache()
  void setGeometry(const o2::mft::GeometryTGeo* geom) { mGeometry = geom; }

  /// \brief return the sensor center global coordinates XYZ
  o2::math_utils::Point3D<double> getSensorCenterGlobalCoordinates(const int chipId) const;

  /// \brief return the half number to which belongs the sensor
  int half(const int chipId) const;

  /// \brief return the disk number to which belongs the sensor
  int disk(const int chipId) const;

  /// \brief return the layer number to which belongs the sensor
  int layer(const int chipId) const;

 protected:
  static o2::itsmft::ChipMappingMFT mChipMapping;                   ///< MFT chip <-> ladder, layer, disk, half mapping
  const o2::mft::GeometryTGeo* mGeometry = nullptr;                 ///< MFT geometry
  static constexpr int mNumberOfSensors = mChipMapping.getNChips(); ///< Total number of sensors (detection elements) in the MFT

  ClassDefNV(AlignSensorHelper, 1);
};

} // namespace mft
} // namespace o2
#endif
