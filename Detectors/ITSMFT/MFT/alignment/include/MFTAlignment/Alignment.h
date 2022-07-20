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
#include <TString.h>

#include "Framework/ProcessingContext.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "MFTAlignment/MillePedeRecord.h"
#include "MFTAlignment/MillePede2.h"

namespace o2
{
namespace mft
{

class AlignPointHelper;

class Alignment
{
 public:
  Alignment() = delete;
  Alignment(bool saveRecordsToFile);
  ~Alignment() = default;

  void init();
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDictionary = d; }
  void setRunNumber(const int value) { mRunNumber = value; }
  void setBz(float bz) { mBz = bz; }
  void processTimeFrame(o2::framework::ProcessingContext& ctx);
  void processRecoTracks();
  bool globalFit();
  void printProcessTrackSummary();

 protected:
  int mRunNumber = 0;
  float mBz = 0;
  int mNumberTFs = 0;
  int mCounterLocalEquationFailed = 0;
  int mCounterSkippedTracks = 0;
  int mCounterUsedTracks = 0;
  static constexpr int mNumberOfTrackParam = 4;                                  ///< Number of track (= local) parameters (X0, Tx, Y0, Ty)
  static constexpr int mNDofPerSensor = 4;                                       ///< translation in global x, y, z, and rotation Rz around global z-axis
  static o2::itsmft::ChipMappingMFT mChipMapping;                                ///< MFT chip <-> ladder, layer, disk, half mapping
  static constexpr int mNumberOfSensors = mChipMapping.getNChips();              ///< Total number of sensors (detection elements) in the MFT
  static constexpr int mNumberOfGlobalParam = mNDofPerSensor * mNumberOfSensors; ///< Number of alignment (= global) parameters
  Double_t mGlobalDerivatives[mNumberOfGlobalParam];                             ///< Array of global derivatives {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}
  Double_t mLocalDerivatives[mNumberOfTrackParam];                               ///< Array of local derivatives {dX0, dTx, dY0, dTz}
  std::array<Double_t, mNDofPerSensor> mAllowVar;                                ///< "Encouraged" variation for degrees of freedom
  Int_t mChi2CutNStdDev = 3;                                                     ///< Number of standard deviations for chi2 cut
  Double_t mResCutInitial = 100.;                                                ///< Cut on residual on first iteration
  Double_t mResCut = 100.;                                                       ///< Cut on residual for other iterations
  int mMinNumberClusterCut = 6;                                                  ///< Minimum number of clusters in the track to be used for alignment
  o2::mft::MillePedeRecord mTrackRecord;                                         ///< running MillePede Track record
  double mWeightRecord = 1.;
  bool mSaveTrackRecordToFile = false;
  TString mMilleRecordsFileName;
  TString mMilleConstraintsRecFileName;
  std::unique_ptr<o2::mft::MillePede2> mMillepede = nullptr;   ///< Millepede2 implementation copied from AliROOT
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr; ///< cluster patterns dictionary
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
  bool setLocalDerivative(Int_t index, Double_t value);

  /// \brief set array of global derivatives
  bool setGlobalDerivative(Int_t index, Double_t value);

  /// \brief reset the array of the Local derivative
  void resetLocalDerivative();

  /// \brief reset the array of the Global derivative
  void resetGlocalDerivative();

  bool setLocalEquationX();
  bool setLocalEquationY();
  bool setLocalEquationZ();

  ClassDefNV(Alignment, 1);
};

} // namespace mft
} // namespace o2

#endif
