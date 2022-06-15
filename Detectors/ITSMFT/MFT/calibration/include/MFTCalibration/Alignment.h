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

/** @file CathodeSegmentation.h
 * C++ Alignmnet .
 * @author  Javier Castillo Castellanos
 */

#ifndef ALICEO2_MCH_ALIGNMENT
#define ALICEO2_MCH_ALIGNMENT

#include <string>

#include "Align/Millepede2Record.h"
#include "Align/Mille.h"
#include "MCHAlign/MillePede2.h"
#include "MCHAlign/MillePedeRecord.h"
#include "TGeoManager.h"

#include "MCHGeometryCreator/Geometry.h"
#include "MCHGeometryTest/Helpers.h"
#include "MCHGeometryTransformer/Transformations.h"

#include "MCHTracking/Track.h"
#include "MCHTracking/Cluster.h"

#include <TFile.h>
#include <TGeoMatrix.h>
#include <TObject.h>
#include <TString.h>
#include <TTree.h>

namespace o2
{

namespace mch
{

/// local track parameters, for refit
class LocalTrackParam
{
 public:
  //* construction
  LocalTrackParam() = default;
  ~LocalTrackParam() = default;

  // private:
  //* y and z
  double fTrackX = 0.0;
  double fTrackY = 0.0;
  double fTrackZ = 0.0;
  double fTrackSlopeX = 0.0;
  double fTrackSlopeY = 0.0;
}; // class LocalTrackParam

class Alignment : public TObject
{

 public:
  Alignment();

  ~Alignment() = default;

  // initialize
  void init(void);

  // terminate
  void terminate(void);

  // array dimendions
  enum {
    /// Number tracking stations
    fgNSt = 5,

    /// Number tracking chambers
    fgNCh = 10,

    /// Number of tracking modules
    fgNTrkMod = 16,

    /// Number of half chambers
    fgNHalfCh = 20,

    /// max number of detector elements per half chamber
    fgNDetHalfChMax = 13,

    /// Total number of detection elements
    /// (4*2 + 4*2 + 18*2 + 26*2 + 26*2)
    fgNDetElem = 156,

    /// Number of local parameters
    fNLocal = 4, // t_x, t_y, x0, y0

    /// Number of degrees of freedom per chamber
    fgNParCh = 4, // x,y,z,phi

    /// Number of global parameters
    fNGlobal = fgNParCh * fgNDetElem
  };

  /// Number of detection elements per chamber
  static const Int_t fgNDetElemCh[fgNCh];

  /// Sum of detection elements up to this chamber
  static const Int_t fgSNDetElemCh[fgNCh + 1];

  /// Number of detection element per tracking module
  static const Int_t fgNDetElemHalfCh[fgNHalfCh];

  /// list of detection elements per tracking module
  static const Int_t fgDetElemHalfCh[fgNHalfCh][fgNDetHalfChMax];

  /// global parameter bit set, used for masks
  enum ParameterMask {
    ParX = 1 << 0,
    ParY = 1 << 1,
    ParZ = 1 << 2,
    ParTZ = 1 << 3,

    ParAllTranslations = ParX | ParY | ParZ,
    ParAllRotations = ParTZ,
    ParAll = ParAllTranslations | ParAllRotations

  };

  /// detector sides bit set, used for selecting sides in constrains
  enum SidesMask {
    SideTop = 1 << 0,
    SideLeft = 1 << 1,
    SideBottom = 1 << 2,
    SideRight = 1 << 3,
    SideTopLeft = SideTop | SideLeft,
    SideTopRight = SideTop | SideRight,
    SideBottomLeft = SideBottom | SideLeft,
    SideBottomRight = SideBottom | SideRight,
    AllSides = SideTop | SideBottom | SideLeft | SideRight
  };

  MillePedeRecord* ProcessTrack(Track& track, Bool_t doAlignment, Double_t weight = 1);

  void ProcessTrack(MillePedeRecord*);

  //@name modifiers
  //@{

  /// run number
  void SetRunNumber(Int_t id)
  {
    fRunNumber = id;
  }

  /// Set flag for Magnetic field On/Off
  void SetBFieldOn(Bool_t value)
  {
    fBFieldOn = value;
  }

  /// set to true to do refit evaluation
  void SetDoEvaluation(Bool_t value)
  {
    fDoEvaluation = value;
  }

  /// set to true to refit tracks
  void SetRefitStraightTracks(Bool_t value)
  {
    fRefitStraightTracks = value;
  }

  void SetAllowedVariation(Int_t iPar, Double_t value);

  void SetSigmaXY(Double_t sigmaX, Double_t sigmaY);

  /// Set geometry transformer
  // void SetGeometryTransformer(AliMUONGeometryTransformer* transformer)
  //{
  //    fTransform = transformer;
  // }

  //@}

  //@name fixing detectors
  //@{

  void FixAll(UInt_t parameterMask = ParAll);

  void FixChamber(Int_t iCh, UInt_t parameterMask = ParAll);

  void FixDetElem(Int_t iDetElemId, UInt_t parameterMask = ParAll);

  void FixHalfSpectrometer(const Bool_t* bChOnOff, UInt_t sidesMask = AllSides, UInt_t parameterMask = ParAll);

  void FixParameter(Int_t iPar);

  void FixParameter(Int_t iDetElem, Int_t iPar)
  {
    FixParameter(iDetElem * fgNParCh + iPar);
  }

  //@}

  //@name releasing detectors
  //@{

  void ReleaseChamber(Int_t iCh, UInt_t parameterMask = ParAll);

  void ReleaseDetElem(Int_t iDetElemId, UInt_t parameterMask = ParAll);

  void ReleaseParameter(Int_t iPar);

  void ReleaseParameter(Int_t iDetElem, Int_t iPar)
  {
    ReleaseParameter(iDetElem * fgNParCh + iPar);
  }

  //@}

  //@name grouping detectors
  //@{

  void GroupChamber(Int_t iCh, UInt_t parameterMask = ParAll);

  void GroupHalfChamber(Int_t iCh, Int_t iHalf, UInt_t parameterMask = ParAll);

  void GroupDetElems(Int_t detElemMin, Int_t detElemMax, UInt_t parameterMask = ParAll);

  void GroupDetElems(const Int_t* detElemList, Int_t nDetElem, UInt_t parameterMask = ParAll);

  //@}

  //@name define non linearity
  //@{

  void SetChamberNonLinear(Int_t iCh, UInt_t parameterMask);

  void SetDetElemNonLinear(Int_t iSt, UInt_t parameterMask);

  void SetParameterNonLinear(Int_t iPar);

  void SetParameterNonLinear(Int_t iDetElem, Int_t iPar)
  {
    SetParameterNonLinear(iDetElem * fgNParCh + iPar);
  }

  //@}

  //@name constraints
  //@{

  void AddConstraints(const Bool_t* bChOnOff, UInt_t parameterMask);

  void AddConstraints(const Bool_t* bChOnOff, const Bool_t* lVarXYT, UInt_t sidesMask = AllSides);

  //@}

  /// initialize global parameters to a give set of values
  void InitGlobalParameters(Double_t* par);

  /// perform global fit
  void GlobalFit(Double_t* parameters, Double_t* errors, Double_t* pulls);

  /// print global parameters
  void PrintGlobalParameters(void) const;

  /// get error on a given parameter
  Double_t GetParError(Int_t iPar) const;

  // AliMUONGeometryTransformer* ReAlign(const AliMUONGeometryTransformer* transformer, const double* misAlignments, Bool_t verbose);

  void SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t chId, Double_t chResX, Double_t chResY, Double_t deResX, Double_t deResY);

 private:
  /// Not implemented
  Alignment(const Alignment& right);

  /// Not implemented
  Alignment& operator=(const Alignment& right);

  /// Set array of local derivatives
  void SetLocalDerivative(Int_t index, Double_t value)
  {
    fLocalDerivatives[index] = value;
  }

  /// Set array of global derivatives
  void SetGlobalDerivative(Int_t index, Double_t value)
  {
    fGlobalDerivatives[index] = value;
  }

  /// refit track using straight track model
  LocalTrackParam RefitStraightTrack(Track&, Double_t) const;

  void FillDetElemData(const Cluster*);

  void FillRecPointData(const Cluster*);

  void FillTrackParamData(const TrackParam*);

  void LocalEquationX(void);

  void LocalEquationY(void);

  TGeoCombiTrans DeltaTransform(const double* detElemMisAlignment) const;

  ///@name utilities
  //@{

  void AddConstraint(Double_t* parameters, Double_t value);

  Int_t GetChamberId(Int_t iDetElemNumber) const;

  Bool_t DetElemIsValid(Int_t iDetElemId) const;

  Int_t GetDetElemNumber(Int_t iDetElemId) const;

  TString GetParameterMaskString(UInt_t parameterMask) const;

  TString GetSidesMaskString(UInt_t sidesMask) const;

  //@}

  /// true when initialized
  Bool_t fInitialized;

  /// current run id
  Int_t fRunNumber;

  /// Flag for Magnetic filed On/Off
  Bool_t fBFieldOn;

  /// true if straight track refit is to be performed
  Bool_t fRefitStraightTracks;

  /// "Encouraged" variation for degrees of freedom
  Double_t fAllowVar[fgNParCh];

  /// Initial value for chi2 cut
  /** if > 1 Iterations in AliMillepede are turned on */
  Double_t fStartFac;

  /// Cut on residual for first iteration
  Double_t fResCutInitial;

  /// Cut on residual for other iterations
  Double_t fResCut;

  /// Detector independent alignment class
  // o2::align::Mille* fMillepede;
  MillePede2* fMillepede; // MillePede2 implementation

  /// running AliMUONVCluster
  o2::mch::Cluster* fCluster;

  /// Number of standard deviations for chi2 cut
  Int_t fNStdDev;

  /// Cluster (global) position
  Double_t fClustPos[3];

  /// Track slope at reference point
  Double_t fTrackSlope0[2];

  /// Track slope at current point
  Double_t fTrackSlope[2];

  /// Track intersection at reference point
  Double_t fTrackPos0[3];

  /// Track intersection at current point
  Double_t fTrackPos[3];

  /// Current measurement (depend on B field On/Off)
  Double_t fMeas[2];

  /// Estimated resolution on measurement
  Double_t fSigma[2];

  /// degrees of freedom
  enum {
    kFixedParId = -1,
    kFreeParId = kFixedParId - 1,
    kGroupBaseId = -10
  };

  /// Array of effective degrees of freedom
  /// it is used to fix detectors, parameters, etc.
  Int_t fGlobalParameterStatus[fNGlobal];

  /// Array of global derivatives
  Double_t fGlobalDerivatives[fNGlobal];

  /// Array of local derivatives
  Double_t fLocalDerivatives[fNLocal];

  /// current detection element number
  Int_t fDetElemNumber;

  /// running Track record
  // o2::align::Millepede2Record fTrackRecord;
  MillePedeRecord fTrackRecord;

  /// Geometry transformation
  // AliMUONGeometryTransformer* fTransform;
  o2::mch::geo::TransformationCreator fTransformCreator;
  TGeoCombiTrans fGeoCombiTransInverse;

  /// preform evaluation
  Bool_t fDoEvaluation;

  /// original local track params
  LocalTrackParam* fTrackParamOrig;
  LocalTrackParam* fTrackParamNew;

  /// output TFile
  TFile* fTFile;

  /// output TTree
  TTree* fTTree;

}; // class Alignment

} // namespace mch
} // namespace o2
#endif // ALICEO2_MCH_ALIGNMENT_H_
