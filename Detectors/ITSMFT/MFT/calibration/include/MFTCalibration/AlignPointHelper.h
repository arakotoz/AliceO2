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

/// \file AlignPointHelper.h
/// \author arakotoz@cern.ch
/// \brief Helper class to compute the local and global derivatives at an alignment point (track position, cluster position)

#ifndef ALICEO2_MFT_ALIGN_POINT_HELPER_H
#define ALICEO2_MFT_ALIGN_POINT_HELPER_H

#include "MathUtils/Cartesian.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "ReconstructionDataFormats/BaseCluster.h"

namespace o2
{
namespace mft
{

class AlignSensorHelper;
class TrackMFT;

/// \class GlobalDerivative
class GlobalDerivative
{
  friend class AlignPointHelper;

 public:
  GlobalDerivative() = default;
  ~GlobalDerivative() = default;

  void reset()
  {
    mdDeltaX = 0.;
    mdDeltaY = 0.;
    mdDeltaZ = 0.;
    mdDeltaRz = 0.;
  }

  double dDeltaX() const { return mdDeltaX; }
  double dDeltaY() const { return mdDeltaY; }
  double dDeltaZ() const { return mdDeltaZ; }
  double dDeltaRz() const { return mdDeltaRz; }

 protected:
  double mdDeltaX = 0.;  ///< derivative w.r.t. delta translation along global x-axis
  double mdDeltaY = 0.;  ///< derivative w.r.t. delta translation along global y-axis
  double mdDeltaZ = 0.;  ///< derivative w.r.t. delta translation along global z-axis
  double mdDeltaRz = 0.; ///< derivative w.r.t. delta rotation angle around global z-axis
};

/// \class LocalDerivative
class LocalDerivative
{
  friend class AlignPointHelper;

 public:
  LocalDerivative() = default;
  ~LocalDerivative() = default;

  double dX0() { return mdX0; }
  double dTx() { return mdTx; }
  double dY0() { return mdY0; }
  double dTy() { return mdTy; }

  void reset()
  {
    mdX0 = 0.;
    mdTx = 0.;
    mdY0 = 0.;
    mdTy = 0.;
  }

 protected:
  double mdX0 = 0.; ///< derivative w.r.t. track param. x0
  double mdTx = 0.; ///< derivative w.r.t. track param. tx
  double mdY0 = 0.; ///< derivative w.r.t. track param. x0
  double mdTy = 0.; ///< derivative w.r.t. track param. ty
};

/// \class AlignPointHelper
class AlignPointHelper
{
 public:
  AlignPointHelper() = delete;

  /// \brief constructor with a pointer to the geometry
  AlignPointHelper(o2::mft::GeometryTGeo* geom);

  ~AlignPointHelper() = default;

  struct TrackParam {
    Double_t X0, Y0, Z0, Tx, Ty;
  };

  void computeLocalDerivatives();
  void computeGlobalDerivatives();

  GlobalDerivative globalDerivativeX() { return mGlobalDerivativeX; }
  GlobalDerivative globalDerivativeY() { return mGlobalDerivativeY; }
  GlobalDerivative globalDerivativeZ() { return mGlobalDerivativeZ; }
  LocalDerivative localDerivativeX() { return mLocalDerivativeX; }
  LocalDerivative localDerivativeY() { return mLocalDerivativeY; }
  LocalDerivative localDerivativeZ() { return mLocalDerivativeZ; }

  void resetAlignPoint();
  void resetDerivatives();
  void resetTrackInitialParam();

  void recordTrackInitialParam(o2::mft::TrackMFT mftTrack);
  void setGlobalRecoPosition(o2::mft::TrackMFT mftTrack);
  void setLocalMeasuredPosition(o2::BaseCluster<float> mftCluster);

 protected:
  bool mIsAlignPointSet = false;        ///< boolean to indicate if mGlobalRecoPosition and mLocalMeasuredPosition are set
  bool mIsGlobalDerivativeDone = false; ///< boolean to indicate if the global derivatives computaion is done
  bool mIsLocalDerivativeDone = false;  ///< boolean to indicate if the local derivatives computation is done
  bool mIsTrackInitialParamSet = false; ///< boolean to indicate if the initial track parameters are recorded

  o2::mft::GeometryTGeo* mGeometry = nullptr; ///< MFT geometry
  std::unique_ptr<o2::mft::AlignSensorHelper> mChipHelper = nullptr;

  LocalDerivative mLocalDerivativeX;
  LocalDerivative mLocalDerivativeY;
  LocalDerivative mLocalDerivativeZ;
  GlobalDerivative mGlobalDerivativeX;
  GlobalDerivative mGlobalDerivativeY;
  GlobalDerivative mGlobalDerivativeZ;

  TrackParam mTrackInitialParam; ///< Track parameters at the reference plane z = z0

  // Current cartesian position (cm, in Global frame) of the reconstructed track
  // analytically propagated to the z position of the cluster
  o2::math_utils::Point3D<double> mGlobalRecoPosition;

  // Current cartesian position (cm, in Local frame) of the cluster
  o2::math_utils::Point3D<double> mLocalMeasuredPosition;

  /// \brief compute X component of the local derivatives
  bool computeLocalDerivativeX();

  /// \brief compute Y component of the local derivatives
  bool computeLocalDerivativeY();

  /// \brief compute Z component of the local derivatives
  bool computeLocalDerivativeZ();

  /// \brief compute X component of the global derivatives
  bool computeGlobalDerivativeX();

  /// \brief compute Y component of the global derivatives
  bool computeGlobalDerivativeY();

  /// \brief compute Z component of the global derivatives
  bool computeGlobalDerivativeZ();

  ClassDefNV(AlignPointHelper, 0);
};

} // namespace mft
} // namespace o2
#endif
