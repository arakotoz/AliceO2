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

/// \file Alignment.h
/// \author arakotoz@cern.ch
/// \brief Class dedicated to standalone alignment for MFT

#ifndef ALICEO2_MFT_ALIGNMENT_H
#define ALICEO2_MFT_ALIGNMENT_H

#include <Rtypes.h>
#include <array>

#include "Framework/ProcessingContext.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "MFTCalibration/MillePedeRecord.h"
#include "MFTCalibration/MillePede2.h"

namespace o2
{
namespace mft
{

class AlignPointHelper;

class Alignment
{
 public:
  Alignment() = default;
  ~Alignment() = default;

  void init(TString dataRecFName, TString consRecFName);
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDictionary = d; }
  void processTimeFrame(o2::framework::ProcessingContext& ctx);

 protected:
  void processRecoTracks();

 protected:
  int mNumberTFs = 0;
  static constexpr int mNumberOfTrackParam = 4;                                  ///< Number of track (= local) parameters (x0, t_x, y0, t_y)
  static constexpr int mNDofPerSensor = 4;                                       ///< translation in x, y, z, and rotation in phi
  static o2::itsmft::ChipMappingMFT mChipMapping;                                ///< MFT chip <-> ladder, layer, disk, half mapping
  static constexpr int mNumberOfSensors = mChipMapping.getNChips();              ///< Total number of sensors (detection elements) in the MFT
  static constexpr int mNumberOfGlobalParam = mNDofPerSensor * mNumberOfSensors; ///< Number of alignment (= global) parameters
  std::array<Double_t, mNumberOfGlobalParam> mGlobalDerivatives;                 ///< Array of global derivatives
  std::array<Double_t, mNumberOfTrackParam> mLocalDerivatives;                   ///< Array of local derivatives
  std::array<Double_t, mNDofPerSensor> mAllowVar;                                ///< "Encouraged" variation for degrees of freedom
  Double_t mSigmaX;                                                              ///< Estimated spatial resolution on measurement in global x direction
  Double_t mSigmaY;                                                              ///< Estimated spatial resolution on measurement in global y direction
  Int_t mChi2CutNStdDev = 3;                                                     ///< Number of standard deviations for chi2 cut
  Double_t mResCutInitial = 100.;                                                ///< Cut on residual on first iteration
  Double_t mResCut = 100.;                                                       ///< Cut on residual for other iterations
  o2::mft::MillePedeRecord mTrackRecord;                                         ///< running MillePede Track record
  std::unique_ptr<o2::mft::MillePede2> mMillepede = nullptr;                     ///< Millepede2 implementation copied from AliROOT
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr;                   ///< cluster patterns dictionary
  gsl::span<const o2::mft::TrackMFT> mMFTTracks;
  gsl::span<const o2::itsmft::ROFRecord> mMFTTracksROF;
  gsl::span<const int> mMFTTrackClusIdx;
  gsl::span<const o2::itsmft::CompClusterExt> mMFTClusters;
  gsl::span<const o2::itsmft::ROFRecord> mMFTClustersROF;
  gsl::span<const unsigned char> mMFTClusterPatterns;
  gsl::span<const unsigned char>::iterator pattIt;
  std::vector<o2::BaseCluster<float>> mMFTClustersGlobal;
  std::unique_ptr<o2::mft::AlignPointHelper> mAlignPoint = nullptr;

  /// \brief set array of local derivatives
  void setLocalDerivative(Int_t index, Double_t value);

  /// \brief set array of global derivatives
  void setGlobalDerivative(Int_t index, Double_t value);

  /// \brief reset the array of the Local derivative
  void resetLocalDerivative();

  /// \brief reset the array of the Global derivative
  void resetGlocalDerivative();

  ClassDefNV(Alignment, 1);
};

} // namespace mft
} // namespace o2
#endif
