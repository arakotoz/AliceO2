/// \file alignHelper.h
/// \author arakotoz@cern.ch
/// \brief Helper class duplicate of Alignment class to test in a macro instead of a workflow

// A.R. DON'T FORGET TO REMOVE BEFORE PULL REQUEST !!!

#ifndef ALICEO2_MFT_ALIGN_HELPER_H
#define ALICEO2_MFT_ALIGN_HELPER_H

#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Rtypes.h>
#include <TString.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "MFTAlignment/MillePedeRecord.h"
#include "MFTAlignment/MillePede2.h"
#include "MFTAlignment/AlignPointHelper.h"
#include "MFTBase/GeometryTGeo.h"

#include "Framework/Logger.h"
#include "MFTAlignment/AlignSensorHelper.h"
#include "MFTTracking/IOUtils.h"
#include "MFTBase/Geometry.h"
#include "MFTAlignment/MillePede2.h"

namespace o2
{
namespace mft
{

struct AlignPoint {
  UShort_t sensor;          // sensor id
  UShort_t layer;           // layer id
  UShort_t disk;            // disk id
  UShort_t half;            // half id
  Double_t measuredGlobalX; // cluster x, global frame (cm)
  Double_t measuredGlobalY; // cluster y, global frame (cm)
  Double_t measuredGlobalZ; // cluster z, global frame (cm)
  Double_t measuredLocalX;  // cluster x, local frame (cm)
  Double_t measuredLocalY;  // cluster y, local frame (cm)
  Double_t measuredLocalZ;  // cluster z, local frame (cm)
  Double_t residualX;       // track global x - cluster global x (cm)
  Double_t residualY;       // track global y - cluster global y (cm)
  Double_t residualZ;       // track global z - cluster global z (cm)
  Double_t residualLocalX;  // track local x - cluster local x (cm)
  Double_t residualLocalY;  // track local y - cluster local y (cm)
  Double_t residualLocalZ;  // track local z - cluster local z (cm)
  Double_t recoGlobalX;     // track x, global frame (cm)
  Double_t recoGlobalY;     // track y, global frame (cm)
  Double_t recoGlobalZ;     // track z, global frame (cm)
  Double_t recoLocalX;      // track x, local frame (cm)
  Double_t recoLocalY;      // track y, local frame (cm)
  Double_t recoLocalZ;      // track z, local frame (cm)
};

class AlignHelper
{
 public:
  struct AlignConfig {
    int minPoints = 6;                  ///< mininum number of clusters in a track used for alignment
    Int_t chi2CutNStdDev = 3;           ///< Number of standard deviations for chi2 cut
    Double_t residualCutInitial = 100.; ///< Cut on residual on first iteration
    Double_t residualCut = 100.;        ///< Cut on residual for other iterations
    Double_t allowedVarDeltaX = 0.5;    ///< allowed max delta in x-translation (cm)
    Double_t allowedVarDeltaY = 0.5;    ///< allowed max delta in y-translation (cm)
    Double_t allowedVarDeltaZ = 0.5;    ///< allowed max delta in z-translation (cm)
    Double_t allowedVarDeltaRz = 0.01;  ///< allowed max delta in rotation around z-axis (rad)
    Double_t chi2CutFactor = 256.;      ///< used to reject outliers i.e. bad tracks with sum(chi2) > Chi2DoFLim(fNStdDev, nDoF) * fChi2CutFactor
  };

 public:
  /// \brief construtor
  AlignHelper();

  /// \brief destructor
  ~AlignHelper();

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

  /// \brief  access mft tracks and clusters in the ROOT files, process all ROFs
  void processROFs(TChain* mfttrackChain, TChain* mftclusterChain);

  /// \brief perform the simultaneous fit of track and alignement parameters
  void globalFit();

  /// \brief print a summary status of what happened in processRecoTracks()
  void printProcessTrackSummary();

  /// \brief provide access to the AlignParam vector
  void getAlignParams(std::vector<o2::detectors::AlignParam>& alignParams) { alignParams = mAlignParams; }

 protected:
  int mRunNumber = 0;                                                            ///< run number
  float mBz = 0;                                                                 ///< magnetic field status
  int mNumberOfClusterChainROFs = 0;                                             ///< number of ROFs in the cluster chain
  int mNumberOfTrackChainROFs = 0;                                               ///< number of ROFs in the track chain
  int mCounterLocalEquationFailed = 0;                                           ///< count how many times we failed to set a local equation
  int mCounterSkippedTracks = 0;                                                 ///< count how many tracks did not met the cut on the min. nb of clusters
  int mCounterUsedTracks = 0;                                                    ///< count how many tracks were used to make Mille records
  static constexpr int mNumberOfTrackParam = 4;                                  ///< Number of track (= local) parameters (X0, Tx, Y0, Ty)
  static constexpr int mNDofPerSensor = 4;                                       ///< translation in global x, y, z, and rotation Rz around global z-axis
  static o2::itsmft::ChipMappingMFT mChipMapping;                                ///< MFT chip <-> ladder, layer, disk, half mapping
  static constexpr int mNumberOfSensors = mChipMapping.getNChips();              ///< Total number of sensors (detection elements) in the MFT
  static constexpr int mNumberOfGlobalParam = mNDofPerSensor * mNumberOfSensors; ///< Number of alignment (= global) parameters
  double* mGlobalDerivatives = nullptr;                                          ///< Array of global derivatives {dDeltaX, dDeltaY, dDeltaRz, dDeltaZ}
  double* mLocalDerivatives = nullptr;                                           ///< Array of local derivatives {dX0, dTx, dY0, dTz}
  std::array<Double_t, mNDofPerSensor> mAllowVar;                                ///< "Encouraged" variation for degrees of freedom {dx, dy, dRz, dz}
  double mStartFac = 256;                                                        ///< Initial value for chi2 cut (if > 1, iterations in Millepede are turned on)
  Int_t mChi2CutNStdDev = 3;                                                     ///< Number of standard deviations for chi2 cut
  Double_t mResCutInitial = 100.;                                                ///< Cut on residual on first iteration
  Double_t mResCut = 100.;                                                       ///< Cut on residual for other iterations
  int mMinNumberClusterCut = 6;                                                  ///< Minimum number of clusters in the track to be used for alignment
  o2::mft::MillePedeRecord mTrackRecord;                                         ///< running MillePede Track record
  double mWeightRecord = 1.;                                                     ///< the weight given to a single Mille record in Millepede algorithm
  TString mMilleRecordsFileName;                                                 ///< output file name when saving the Mille records
  TString mMilleConstraintsRecFileName;                                          ///< output file name when saving the records of the constraints
  std::unique_ptr<o2::mft::MillePede2> mMillepede = nullptr;                     ///< Millepede2 implementation copied from AliROOT
  const o2::itsmft::TopologyDictionary* mDictionary = nullptr;                   ///< cluster patterns dictionary
  std::unique_ptr<o2::mft::AlignPointHelper> mAlignPoint = nullptr;              ///< AlignHelper point helper
  std::vector<o2::detectors::AlignParam> mAlignParams;                           ///< vector of alignment parameters computed by Millepede global fit
  bool mIsInitDone = false;                                                      ///< boolean to follow the initialisation status
  int* mGlobalParameterStatus = nullptr;                                         ///< Array of effective degrees of freedom, used to fix detectors, parameters, etc.
  bool mWithControl = false;                                                     ///< boolean to set the use of the control tree

  // used to fix some degrees of freedom

  static constexpr Int_t mFixedParId = -1;
  static constexpr Int_t mFreeParId = mFixedParId - 1;

  TFile* mControlFile = nullptr;
  TTree* mControlTree = nullptr;
  AlignPoint mPointInfo;
  void initControlTree();
  void closeControlTree();
  void fillControlTree();

  /// \brief set array of local derivatives
  bool setLocalDerivative(int index, double value);

  /// \brief set array of global derivatives
  bool setGlobalDerivative(int index, double value);

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

  ClassDefNV(AlignHelper, 1);
};

} // namespace mft
} // namespace o2
#endif
