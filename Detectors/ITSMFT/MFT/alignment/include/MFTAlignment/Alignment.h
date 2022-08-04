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
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include "Framework/ProcessingContext.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "MFTAlignment/MillePedeRecord.h"
#include "MFTAlignment/MilleRecordReader.h"
#include "MFTAlignment/MilleRecordWriter.h"
#include "MFTAlignment/MillePede2.h"
#include "MFTAlignment/AlignPointHelper.h"
#include "MFTAlignment/AlignPointControl.h"
#include "MFTBase/GeometryTGeo.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

namespace o2
{
namespace mft
{

class Alignment
{
 public:
  /// \brief construtor
  Alignment();

  /// \brief destructor
  ~Alignment();

  /// \brief init Millipede and AlignPointHelper
  void init();

  // simple setters

  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDictionary = d; }
  void setRunNumber(const int value) { mRunNumber = value; }
  void setBz(const float bz) { mBz = bz; }
  void setChi2CutNStdDev(const Int_t value) { mChi2CutNStdDev = value; }
  void setResidualCutInitial(const Double_t value) { mResCutInitial = value; }
  void setResidualCut(const Double_t value) { mResCut = value; }
  void setMinNumberClusterCut(const int value) { mMinNumberClusterCut = value; }
  void setAllowedVariationDeltaX(const double value) { mAllowVar[0] = value; }
  void setAllowedVariationDeltaY(const double value) { mAllowVar[1] = value; }
  void setAllowedVariationDeltaZ(const double value) { mAllowVar[3] = value; }
  void setAllowedVariationDeltaRz(const double value) { mAllowVar[2] = value; }
  void setChi2CutFactor(const double value) { mStartFac = value; }
  void setWithControl(const bool choice) { mWithControl = choice; }
  void setNEntriesAutoSave(const int value) { mNEntriesAutoSave = value; }
  void setWithRecordWriter(const bool choice) { mWithRecordWriter = choice; }
  void setWithConstraintsRecWriter(const bool choice) { mWithConstraintsRecWriter = choice; }
  void setWithRecordReader(const bool choice) { mWithRecordReader = choice; }
  void setWithConstraintsRecReader(const bool choice) { mWithConstraintsRecReader = choice; }

  /// \brief access mft tracks and clusters in the timeframe provided by the workflow
  void processTimeFrame(o2::framework::ProcessingContext& ctx);

  /// \brief use valid tracks (and associated clusters) from the workflow to build Mille records
  void processRecoTracks();

  /// \brief  access mft tracks and clusters provided by ROOT files, process all ROFs
  void processROFs(TChain* mfttrackChain, TChain* mftclusterChain);

  /// \brief perform the simultaneous fit of track and alignement parameters
  void globalFit();

  /// \brief print a summary status of what happened in processRecoTracks()
  void printProcessTrackSummary();

  /// \brief provide access to the AlignParam vector
  void getAlignParams(std::vector<o2::detectors::AlignParam>& alignParams) { alignParams = mAlignParams; }

  /// \brief init all related to writing data records and its control tree
  void startRecordWriter();

  /// \brief end all related to writing data records and its control tree
  void endRecordWriter();

  /// \brief init all related to writing constraints records
  void startConstraintsRecWriter();

  /// \brief end all related to writing constraints records
  void endConstraintsRecWriter();

  /// \brief connect data record reader to input TChain of records
  void connectRecordReaderToChain(TChain* ch);

  /// \brief conect constraints record reader to input TChain of constraints record
  void connectConstraintsRecReaderToChain(TChain* ch);

 protected:
  int mRunNumber;                                                                ///< run number
  float mBz;                                                                     ///< magnetic field status
  int mNumberTFs;                                                                ///< number of timeframes processed
  int mNumberOfClusterChainROFs;                                                 ///< number of ROFs in the cluster chain
  int mNumberOfTrackChainROFs;                                                   ///< number of ROFs in the track chain
  int mCounterLocalEquationFailed;                                               ///< count how many times we failed to set a local equation
  int mCounterSkippedTracks;                                                     ///< count how many tracks did not met the cut on the min. nb of clusters
  int mCounterUsedTracks;                                                        ///< count how many tracks were used to make Mille records
  static constexpr int mNumberOfTrackParam = 4;                                  ///< Number of track (= local) parameters (X0, Tx, Y0, Ty)
  static constexpr int mNDofPerSensor = 4;                                       ///< translation in global x, y, z, and rotation Rz around global z-axis
  static o2::itsmft::ChipMappingMFT mChipMapping;                                ///< MFT chip <-> ladder, layer, disk, half mapping
  static constexpr int mNumberOfSensors = mChipMapping.getNChips();              ///< Total number of sensors (detection elements) in the MFT
  static constexpr int mNumberOfGlobalParam = mNDofPerSensor * mNumberOfSensors; ///< Number of alignment (= global) parameters
  double* mGlobalDerivatives;                                                    ///< Array of global derivatives {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}
  double* mLocalDerivatives;                                                     ///< Array of local derivatives {dX0, dTx, dY0, dTz}
  std::array<double, mNDofPerSensor> mAllowVar;                                  ///< "Encouraged" variation for degrees of freedom {dx, dy, dRz, dz}
  double mStartFac;                                                              ///< Initial value for chi2 cut, used to reject outliers i.e. bad tracks with sum(chi2) > Chi2DoFLim(fNStdDev, nDoF) * chi2CutFactor (if > 1, iterations in Millepede are turned on)
  int mChi2CutNStdDev;                                                           ///< Number of standard deviations for chi2 cut
  double mResCutInitial;                                                         ///< Cut on residual on first iteration
  double mResCut;                                                                ///< Cut on residual for other iterations
  int mMinNumberClusterCut;                                                      ///< Minimum number of clusters in the track to be used for alignment
  double mWeightRecord;                                                          ///< the weight given to a single Mille record in Millepede algorithm
  TString mMilleRecordsFileName;                                                 ///< output file name when saving the Mille records
  TString mMilleConstraintsRecFileName;                                          ///< output file name when saving the records of the constraints
  std::unique_ptr<o2::mft::MillePede2> mMillepede;                               ///< Millepede2 implementation copied from AliROOT
  const o2::itsmft::TopologyDictionary* mDictionary;                             ///< cluster patterns dictionary
  std::shared_ptr<o2::mft::AlignPointHelper> mAlignPoint;                        ///< Alignment point helper
  std::vector<o2::detectors::AlignParam> mAlignParams;                           ///< vector of alignment parameters computed by Millepede global fit
  bool mIsInitDone = false;                                                      ///< boolean to follow the initialisation status
  int* mGlobalParameterStatus;                                                   ///< Array of effective degrees of freedom, used to fix detectors, parameters, etc.
  bool mWithControl;                                                             ///< boolean to set the use of the control tree
  long mNEntriesAutoSave = 10000;                                                ///< number of entries needed to call AutoSave for the output TTrees
  o2::mft::AlignPointControl mPointControl;                                      ///< AlignPointControl handles the control tree
  bool mWithRecordWriter;
  std::shared_ptr<o2::mft::MilleRecordWriter> mRecordWriter;
  bool mWithConstraintsRecWriter;
  std::shared_ptr<o2::mft::MilleRecordWriter> mConstraintsRecWriter;
  bool mWithRecordReader;
  std::shared_ptr<o2::mft::MilleRecordReader> mRecordReader;
  bool mWithConstraintsRecReader;
  std::shared_ptr<o2::mft::MilleRecordReader> mConstraintsRecReader;

  // used to fix some degrees of freedom

  static constexpr int mFixedParId = -1;
  static constexpr int mFreeParId = mFixedParId - 1;

  // access these data from CTFs

  gsl::span<const o2::mft::TrackMFT> mMFTTracks;
  gsl::span<const o2::itsmft::ROFRecord> mMFTTracksROF;
  gsl::span<const int> mMFTTrackClusIdx;
  gsl::span<const o2::itsmft::CompClusterExt> mMFTClusters;
  gsl::span<const o2::itsmft::ROFRecord> mMFTClustersROF;
  gsl::span<const unsigned char> mMFTClusterPatterns;
  gsl::span<const unsigned char>::iterator pattIt;

  /// \brief set array of local derivatives
  bool setLocalDerivative(Int_t index, Double_t value);

  /// \brief set array of global derivatives
  bool setGlobalDerivative(Int_t index, Double_t value);

  /// \brief reset the array of the Local derivative
  bool resetLocalDerivative();

  /// \brief reset the array of the Global derivative
  bool resetGlocalDerivative();

  /// \brief set the first component of the local equation vector for a given alignment point
  bool setLocalEquationX();

  /// \brief set the 2nd component of the local equation vector for a given alignment point
  bool setLocalEquationY();

  /// \brief set the last component of the local equation vector for a given alignment point
  bool setLocalEquationZ();

  ClassDefNV(Alignment, 1);
};

} // namespace mft
} // namespace o2

#endif
