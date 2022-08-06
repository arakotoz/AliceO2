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

/// \file TracksToRecords.h
/// \author arakotoz@cern.ch
/// \brief Class responsible to create records from tracks (and attached clusters) to feed alignment

#ifndef ALICEO2_MFT_TRACKS_TO_RECORDS_H
#define ALICEO2_MFT_TRACKS_TO_RECORDS_H

#include <vector>

#include <TTree.h>
#include <TChain.h>

#include "Framework/ProcessingContext.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "MFTAlignment/MillePedeRecord.h"
#include "MFTAlignment/MilleRecordWriter.h"
#include "MFTAlignment/AlignPointHelper.h"
#include "MFTAlignment/AlignPointControl.h"
#include "MFTBase/GeometryTGeo.h"

#include "MFTAlignment/Aligner.h"

namespace o2
{
namespace mft
{

class TracksToRecords : public Aligner
{
 public:
  /// \brief construtor
  TracksToRecords();

  /// \brief destructor
  ~TracksToRecords();

  /// \brief init Millipede and AlignPointHelper
  void init();

  // simple setters

  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDictionary = d; }
  void setRunNumber(const int value) { mRunNumber = value; }
  void setBz(const float bz) { mBz = bz; }
  void setMinNumberClusterCut(const int value) { mMinNumberClusterCut = value; }
  void setWithControl(const bool choice) { mWithControl = choice; }
  void setNEntriesAutoSave(const int value) { mNEntriesAutoSave = value; }
  void setWithConstraintsRecWriter(const bool choice) { mWithConstraintsRecWriter = choice; }

  /// \brief access mft tracks and clusters in the timeframe provided by the workflow
  void processTimeFrame(o2::framework::ProcessingContext& ctx);

  /// \brief use valid tracks (and associated clusters) from the workflow to build Mille records
  void processRecoTracks();

  /// \brief  access mft tracks and clusters provided by ROOT files, process all ROFs
  void processROFs(TChain* mfttrackChain, TChain* mftclusterChain);

  /// \brief print a summary status of what happened in processRecoTracks()
  void printProcessTrackSummary();

  /// \brief init all related to writing data records and its control tree
  void startRecordWriter();

  /// \brief end all related to writing data records and its control tree
  void endRecordWriter();

  /// \brief init all related to writing constraints records
  void startConstraintsRecWriter();

  /// \brief end all related to writing constraints records
  void endConstraintsRecWriter();

 protected:
  int mRunNumber;                                         ///< run number
  float mBz;                                              ///< magnetic field status
  int mNumberTFs;                                         ///< number of timeframes processed
  int mNumberOfClusterChainROFs;                          ///< number of ROFs in the cluster chain
  int mNumberOfTrackChainROFs;                            ///< number of ROFs in the track chain
  int mCounterLocalEquationFailed;                        ///< count how many times we failed to set a local equation
  int mCounterSkippedTracks;                              ///< count how many tracks did not met the cut on the min. nb of clusters
  int mCounterUsedTracks;                                 ///< count how many tracks were used to make Mille records
  std::vector<double> mGlobalDerivatives;                 ///< vector of global derivatives {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}
  std::vector<double> mLocalDerivatives;                  ///< vector of local derivatives {dX0, dTx, dY0, dTz}
  int mMinNumberClusterCut;                               ///< Minimum number of clusters in the track to be used for alignment
  double mWeightRecord;                                   ///< the weight given to a single Mille record in Millepede algorithm
  const o2::itsmft::TopologyDictionary* mDictionary;      ///< cluster patterns dictionary
  std::shared_ptr<o2::mft::AlignPointHelper> mAlignPoint; ///< Alignment point helper
  bool mWithControl;                                      ///< boolean to set the use of the control tree
  long mNEntriesAutoSave = 10000;                         ///< number of entries needed to call AutoSave for the output TTrees
  o2::mft::AlignPointControl mPointControl;               ///< AlignPointControl handles the control tree
  std::shared_ptr<o2::mft::MilleRecordWriter> mRecordWriter;
  bool mWithConstraintsRecWriter;
  std::shared_ptr<o2::mft::MilleRecordWriter> mConstraintsRecWriter;

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

  ClassDef(TracksToRecords, 0);
};

} // namespace mft
} // namespace o2

#endif
