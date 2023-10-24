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

#ifndef O2_EMCAL_OFFLINECALIB_SPEC
#define O2_EMCAL_OFFLINECALIB_SPEC

#include <string>
#include "Framework/DataProcessorSpec.h"
#include "EMCALWorkflow/CalibLoader.h"
#include "EMCALCalib/GainCalibrationFactors.h"
#include "EMCALCalibration/EMCALCalibParams.h"
#include "DataFormatsCTP/Digits.h"
#include "DataFormatsCTP/Configuration.h"
#include "Framework/Task.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"

namespace o2
{

namespace emcal
{

/// \class OfflineCalibSpec
/// \brief Task for producing offline calibration objects
/// \ingroup EMCALworkflow
/// \author Hannah Bossi <hannah.bossi@cern.ch>, Yale University
/// \since August 16th, 2022
///
/// This task fills offline calibration objects for the EMCAL.
class OfflineCalibSpec : public framework::Task
{
  using EMCALCalibParams = o2::emcal::EMCALCalibParams;

 public:
  /// \brief Constructor
  /// \param makeCellIDTimeEnergy If true the THnSparseF of cell ID, time, and energy is made
  /// \param rejectCalibTriggers if true, only events which have the o2::trigger::PhT flag will be taken into account
  OfflineCalibSpec(bool makeCellIDTimeEnergy, bool rejectCalibTriggers, bool rejectL0Trigger, std::shared_ptr<o2::emcal::CalibLoader> calibHandler) : mMakeCellIDTimeEnergy(makeCellIDTimeEnergy), mRejectCalibTriggers(rejectCalibTriggers), mRejectL0Triggers(rejectL0Trigger), mCalibrationHandler(calibHandler){};

  /// \brief Destructor
  ~OfflineCalibSpec() override = default;

  /// \brief Initializing the offline calib task
  /// \param ctx Init context
  void init(framework::InitContext& ctx) final;

  /// Loading CCDB object into internal cache for the 3 supported CCDB entries (bad
  /// channel map, time calibration params, gain calibration params). Objects are only
  /// loaded in case the calibration type was enabled.
  void finaliseCCDB(framework::ConcreteDataMatcher& matcher, void* obj) final;

  /// \brief Fill histograms needed for the offline calibration
  /// \param ctx Processing context
  ///
  void run(framework::ProcessingContext& ctx) final;

  /// \brief Write histograms to an output root file
  /// \param ec end of stream context
  ///
  void endOfStream(o2::framework::EndOfStreamContext& ec) final;

  static const char* getCTPDigitsBinding() { return "CTPDigits"; }
  static const char* getCTPConfigBinding() { return "CTPConfig"; }

 private:
  /// \brief Update calibration objects (if changed)
  void updateCalibObjects();
  std::unique_ptr<TH2> mCellAmplitude;              ///< Cell energy vs. cell ID
  std::unique_ptr<TH2> mCellTime;                   ///< Cell time vs. cell ID
  std::unique_ptr<TH2> mCellTimeLG;                 ///< Cell time vs. cell ID for low gain cells
  std::unique_ptr<TH2> mCellTimeHG;                 ///< Cell time vs. cell ID for high gain cells
  std::unique_ptr<TH1> mNevents;                    ///< Number of events
  std::unique_ptr<THnSparseF> mCellTimeEnergy;      ///< ID, time, energy
  std::shared_ptr<CalibLoader> mCalibrationHandler; ///< Handler loading calibration objects
  GainCalibrationFactors* mGainCalibFactors;        ///< cell gain calibration factors
  bool mMakeCellIDTimeEnergy = true;                ///< Switch whether or not to make a THnSparseF of cell ID, time, and energy
  bool mRejectCalibTriggers = true;                 ///< Switch to select if calib triggerred events should be rejected
  bool mEnableGainCalib = false;                    ///< Switch weather gain calibration should be applied or not when filling the hsitograms
  bool mRejectL0Triggers = false;                   ///< Switch to select if L0 triggerred events should be rejected
  std::vector<uint64_t> mSelectedClassMasks = {};   ///< vector with selected trigger masks that are accepted for processing
};

/// \brief Creating offline calib spec
/// \ingroup EMCALworkflow
///
o2::framework::DataProcessorSpec getEmcalOfflineCalibSpec(bool makeCellIDTimeEnergy, bool rejectCalibTriggers, bool rejectL0Trigger, uint32_t inputsubspec, bool enableGainCalib, bool ctpcfgperrun);

} // namespace emcal

} // namespace o2

#endif // O2_EMCAL_OFFLINECALIB_SPEC
