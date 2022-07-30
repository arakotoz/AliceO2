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

/// \file MillePede2.h
/// \author ruben.shahoyan@cern.ch
/// \brief General class for alignment with large number of degrees of freedom
///
/// Based on the original milliped2 by Volker Blobel
/// http://www.desy.de/~blobel/mptalks.html

#ifndef ALICEO2_MFT_MILLEPEDE2_H
#define ALICEO2_MFT_MILLEPEDE2_H

#include <TObject.h>
#include <TString.h>
#include <TTree.h>
#include "MFTAlignment/MinResSolve.h"
#include "MFTAlignment/MillePedeRecord.h"
#include "MFTAlignment/SymMatrix.h"
#include "MFTAlignment/RectMatrix.h"
#include "MFTAlignment/MatrixSparse.h"
#include "MFTAlignment/MatrixSq.h"

class TFile;
class TStopwatch;
class TArrayL;
class TArrayF;

namespace o2
{
namespace mft
{

class MillePede2 : public TObject
{
 public:
  //
  enum { kFailed,
         kInvert,
         kNoInversion };   // used global matrix solution methods
  enum { kFixParID = -1 }; // dummy id for fixed param

  MillePede2();
  MillePede2(const MillePede2& src);
  virtual ~MillePede2();
  MillePede2& operator=(const MillePede2&)
  {
    printf("Dummy\n");
    return *this;
  }

  /// \brief init all
  Int_t InitMille(int nGlo, int nLoc, Int_t lNStdDev = -1, double lResCut = -1., double lResCutInit = -1., const Int_t* regroup = 0);

  Int_t GetNGloPar() const { return fNGloPar; }
  Int_t GetNGloParIni() const { return fNGloParIni; }
  const Int_t* GetRegrouping() const { return fkReGroup; }
  Int_t GetNLocPar() const { return fNLocPar; }
  Long_t GetNLocalEquations() const { return fNLocEquations; }
  Int_t GetCurrentIteration() const { return fIter; }
  Int_t GetNMaxIterations() const { return fMaxIter; }
  Int_t GetNStdDev() const { return fNStdDev; }
  Int_t GetNGlobalConstraints() const { return fNGloConstraints; }
  Int_t GetNLagrangeConstraints() const { return fNLagrangeConstraints; }
  Long_t GetNLocalFits() const { return fNLocFits; }
  Long_t GetNLocalFitsRejected() const { return fNLocFitsRejected; }
  Int_t GetNGlobalsFixed() const { return fNGloFix; }
  Int_t GetGlobalSolveStatus() const { return fGloSolveStatus; }
  Float_t GetChi2CutFactor() const { return fChi2CutFactor; }
  Float_t GetChi2CutRef() const { return fChi2CutRef; }
  Float_t GetResCurInit() const { return fResCutInit; }
  Float_t GetResCut() const { return fResCut; }
  Int_t GetMinPntValid() const { return fMinPntValid; }
  Int_t GetRGId(Int_t i) const { return fkReGroup ? (fkReGroup[i] < 0 ? -1 : fkReGroup[i]) : i; }
  Int_t GetProcessedPoints(Int_t i) const
  {
    int ir = GetRGId(i);
    return ir <= 0 ? 0 : fProcPnt[ir];
  }
  Int_t* GetProcessedPoints() const { return fProcPnt; }
  Int_t GetParamGrID(Int_t i) const
  {
    int ir = GetRGId(i);
    return ir <= 0 ? 0 : fParamGrID[ir];
  }

  MatrixSq* GetGlobalMatrix() const { return fMatCGlo; }
  SymMatrix* GetLocalMatrix() const { return fMatCLoc; }
  Double_t* GetGlobals() const { return fVecBGlo; }
  Double_t* GetDeltaPars() const { return fDeltaPar; }
  Double_t* GetInitPars() const { return fInitPar; }
  Double_t* GetSigmaPars() const { return fSigmaPar; }
  Bool_t* GetIsLinear() const { return fIsLinear; }
  Double_t GetFinalParam(int i) const
  {
    int ir = GetRGId(i);
    return ir < 0 ? 0 : fDeltaPar[ir] + fInitPar[ir];
  }
  Double_t GetFinalError(int i) const { return GetParError(i); }

  /// \brief return pull for parameter iPar
  Double_t GetPull(int i) const;

  Double_t GetGlobal(Int_t i) const
  {
    int ir = GetRGId(i);
    return ir < 0 ? 0 : fVecBGlo[ir];
  }
  Double_t GetInitPar(Int_t i) const
  {
    int ir = GetRGId(i);
    return ir < 0 ? 0 : fInitPar[ir];
  }
  Double_t GetSigmaPar(Int_t i) const
  {
    int ir = GetRGId(i);
    return ir < 0 ? 0 : fSigmaPar[ir];
  }
  Bool_t GetIsLinear(Int_t i) const
  {
    int ir = GetRGId(i);
    return ir < 0 ? 0 : fIsLinear[ir];
  }
  static Bool_t IsGlobalMatSparse() { return fgIsMatGloSparse; }
  static Bool_t IsWeightSigma() { return fgWeightSigma; }
  void SetWghScale(Double_t wOdd = 1, Double_t wEven = 1)
  {
    fWghScl[0] = wOdd;
    fWghScl[1] = wEven;
  }
  void SetUseRecordWeight(Bool_t v = kTRUE) { fUseRecordWeight = v; }
  Bool_t GetUseRecordWeight() const { return fUseRecordWeight; }
  void SetMinRecordLength(Int_t v = 1) { fMinRecordLength = v; }
  Int_t GetMinRecordLength() const { return fMinRecordLength; }

  void SetParamGrID(Int_t grID, Int_t i)
  {
    int ir = GetRGId(i);
    if (ir < 0)
      return;
    fParamGrID[ir] = grID;
    if (fNGroupsSet < grID)
      fNGroupsSet = grID;
  }
  void SetNGloPar(Int_t n) { fNGloPar = n; }
  void SetNLocPar(Int_t n) { fNLocPar = n; }
  void SetNMaxIterations(Int_t n = 10) { fMaxIter = n; }
  void SetNStdDev(Int_t n) { fNStdDev = n; }
  void SetChi2CutFactor(Float_t v) { fChi2CutFactor = v; }
  void SetChi2CutRef(Float_t v) { fChi2CutRef = v; }
  void SetResCurInit(Float_t v) { fResCutInit = v; }
  void SetResCut(Float_t v) { fResCut = v; }
  void SetMinPntValid(Int_t n) { fMinPntValid = n > 0 ? n : 1; }
  static void SetGlobalMatSparse(Bool_t v = kTRUE) { fgIsMatGloSparse = v; }
  static void SetWeightSigma(Bool_t v = kTRUE) { fgWeightSigma = v; }

  /// \brief initialize parameters, account for eventual grouping
  void SetInitPars(const Double_t* par);

  /// \brief initialize sigmas, account for eventual grouping
  void SetSigmaPars(const Double_t* par);

  /// \brief initialize param, account for eventual grouping
  void SetInitPar(Int_t i, Double_t par);

  /// \brief initialize sigma, account for eventual grouping
  void SetSigmaPar(Int_t i, Double_t par);

  /// \brief performs a requested number of global iterations
  Int_t GlobalFit(Double_t* par = 0, Double_t* error = 0, Double_t* pull = 0);

  /// \brief perform global parameters fit once all the local equations have been fitted
  Int_t GlobalFitIteration();

  /// \brief solve global matrix equation MatCGlob*X=VecBGlo and store the result in the VecBGlo
  Int_t SolveGlobalMatEq();

  static void SetInvChol(Bool_t v = kTRUE) { fgInvChol = v; }
  static void SetMinResPrecondType(Int_t tp = 0) { fgMinResCondType = tp; }
  static void SetMinResTol(Double_t val = 1e-12) { fgMinResTol = val; }
  static void SetMinResMaxIter(Int_t val = 2000) { fgMinResMaxIter = val; }
  static void SetIterSolverType(Int_t val = MinResSolve::kSolMinRes) { fgIterSol = val; }
  static void SetNKrylovV(Int_t val = 60) { fgNKrylovV = val; }

  static Bool_t GetInvChol() { return fgInvChol; }
  static Int_t GetMinResPrecondType() { return fgMinResCondType; }
  static Double_t GetMinResTol() { return fgMinResTol; }
  static Int_t GetMinResMaxIter() { return fgMinResMaxIter; }
  static Int_t GetIterSolverType() { return fgIterSol; }
  static Int_t GetNKrylovV() { return fgNKrylovV; }

  /// \brief return error for parameter iPar
  Double_t GetParError(int iPar) const;

  /// \brief print the final results into the logfile
  Int_t PrintGlobalParameters() const;

  /// \brief set the list of runs to be rejected
  void SetRejRunList(const UInt_t* runs, Int_t nruns);

  /// \brief set the list of runs to be selected
  void SetAccRunList(const UInt_t* runs, Int_t nruns, const Float_t* wghList = 0);

  /// \brief validate record according run lists set by the user
  Bool_t IsRecordAcceptable();

  /// \brief Number of iterations is calculated from lChi2CutFac
  Int_t SetIterations(double lChi2CutFac);

  // constraints

  /// \brief define a constraint equation
  void SetGlobalConstraint(const double* dergb, double val, double sigma = 0);

  /// \brief define a constraint equation
  void SetGlobalConstraint(const int* indgb, const double* dergb, int ngb, double val, double sigma = 0);

  // processing of the local measurement

  /// \brief assign run
  void SetRecordRun(Int_t run);

  /// \brief assign weight
  void SetRecordWeight(double wgh);

  /// \brief assing derivs of loc.eq.
  void SetLocalEquation(double* dergb, double* derlc, double lMeas, double lSigma, bool wDebugPrint = false);

  /// \brief write data of single measurement.
  ///        Note: the records ignore regrouping, store direct parameters
  void SetLocalEquation(int* indgb, double* dergb, int ngb, int* indlc,
                        double* derlc, int nlc, double lMeas, double lSigma);

  // manipilation with processed data and costraints records and its buffer
  void SetDataRecFName(const char* flname) { fDataRecFName = flname; }
  const Char_t* GetDataRecFName() const { return fDataRecFName.Data(); }
  void SetConsRecFName(const char* flname) { fConstrRecFName = flname; }
  const Char_t* GetConsRecFName() const { return fConstrRecFName.Data(); }
  const Char_t* GetRecChi2FName() const { return fRecChi2FName.Data(); }
  //
  void SetRecDataTreeName(const char* name = 0)
  {
    fRecDataTreeName = name;
    if (fRecDataTreeName.IsNull())
      fRecDataTreeName = "MillePedeRecords_Data";
  }
  void SetRecConsTreeName(const char* name = 0)
  {
    fRecConsTreeName = name;
    if (fRecConsTreeName.IsNull())
      fRecConsTreeName = "MillePedeRecords_Consaints";
  }
  void SetRecDataBranchName(const char* name = 0)
  {
    fRecDataBranchName = name;
    if (fRecDataBranchName.IsNull())
      fRecDataBranchName = "Record_Data";
  }
  void SetRecConsBranchName(const char* name = 0)
  {
    fRecConsBranchName = name;
    if (fRecConsBranchName.IsNull())
      fRecConsBranchName = "Record_Consaints";
  }
  const char* GetRecDataTreeName() const { return fRecDataTreeName.Data(); }
  const char* GetRecConsTreeName() const { return fRecConsTreeName.Data(); }
  const char* GetRecDataBranchName() const { return fRecDataBranchName.Data(); }
  const char* GetRecConsBranchName() const { return fRecConsBranchName.Data(); }

  /// \brief initialize the buffer for processed measurements records
  Bool_t InitDataRecStorage(Bool_t read = kFALSE);

  /// \brief initialize the buffer for processed measurements records
  Bool_t InitConsRecStorage(Bool_t read = kFALSE);

  /// \brief set filename for records
  Bool_t ImposeDataRecFile(const char* fname);

  /// \brief set filename for constraints
  Bool_t ImposeConsRecFile(const char* fname);

  /// \brief close records file
  void CloseDataRecStorage();

  /// \brief close constraints file
  void CloseConsRecStorage();

  void ReadRecordData(Long_t recID);
  void ReadRecordConstraint(Long_t recID);

  /// \brief read next data record (if any)
  Bool_t ReadNextRecordData();

  /// \brief read next constraint record (if any)
  Bool_t ReadNextRecordConstraint();

  void SaveRecordData();
  void SaveRecordConstraint();
  MillePedeRecord* GetRecord() const { return fRecord; }
  Long_t GetSelFirst() const { return fSelFirst; }
  Long_t GetSelLast() const { return fSelLast; }
  void SetSelFirst(Long_t v) { fSelFirst = v; }
  void SetSelLast(Long_t v) { fSelLast = v; }

  /// \brief return the limit in chi^2/nd for n sigmas stdev authorized
  ///
  /// Only n=1, 2, and 3 are expected in input
  Float_t Chi2DoFLim(int nSig, int nDoF) const;

  // aliases for compatibility with millipede1
  void SetParSigma(Int_t i, Double_t par) { SetSigmaPar(i, par); }
  void SetGlobalParameters(Double_t* par) { SetInitPars(par); }
  void SetNonLinear(int index, Bool_t v = kTRUE)
  {
    int id = GetRGId(index);
    if (id < 0)
      return;
    fIsLinear[id] = !v;
  }

 protected:
  /// \brief Perform local parameters fit once all the local equations have been set
  ///
  /// localParams = (if !=0) will contain the fitted track parameters and related errors
  Int_t LocalFit(double* localParams = 0);

  Bool_t IsZero(Double_t v, Double_t eps = 1e-16) const { return TMath::Abs(v) < eps; }

 protected:
  Int_t fNLocPar;    ///< number of local parameters
  Int_t fNGloPar;    ///< number of global parameters
  Int_t fNGloParIni; ///< number of global parameters before grouping
  Int_t fNGloSize;   ///< final size of the global matrix (NGloPar+NConstraints)

  Long_t fNLocEquations;       ///< Number of local equations
  Int_t fIter;                 ///< Current iteration
  Int_t fMaxIter;              ///< Maximum number of iterations
  Int_t fNStdDev;              ///< Number of standard deviations for chi2 cut
  Int_t fNGloConstraints;      ///< Number of constraint equations
  Int_t fNLagrangeConstraints; ///< Number of constraint equations requiring Lagrange multiplier
  Long_t fNLocFits;            ///< Number of local fits
  Long_t fNLocFitsRejected;    ///< Number of local fits rejected
  Int_t fNGloFix;              ///< Number of globals fixed by user
  Int_t fGloSolveStatus;       ///< Status of global solver at current step

  Float_t fChi2CutFactor; ///< Cut factor for chi2 cut to accept local fit
  Float_t fChi2CutRef;    ///< Reference cut for chi2 cut to accept local fit
  Float_t fResCutInit;    ///< Cut in residual for first iterartion
  Float_t fResCut;        ///< Cut in residual for other iterartiona
  Int_t fMinPntValid;     ///< min number of points for global to vary

  Int_t fNGroupsSet;   ///< number of groups set
  Int_t* fParamGrID;   ///< [fNGloPar] group id for the every parameter
  Int_t* fProcPnt;     ///< [fNGloPar] N of processed points per global variable
  Double_t* fVecBLoc;  ///< [fNLocPar] Vector B local (parameters)
  Double_t* fDiagCGlo; ///< [fNGloPar] Initial diagonal elements of C global matrix
  Double_t* fVecBGlo;  //! Vector B global (parameters)

  Double_t* fInitPar;  ///< [fNGloPar] Initial global parameters
  Double_t* fDeltaPar; ///< [fNGloPar] Variation of global parameters
  Double_t* fSigmaPar; ///< [fNGloPar] Sigma of allowed variation of global parameter

  Bool_t* fIsLinear;   ///< [fNGloPar] Flag for linear parameters
  Bool_t* fConstrUsed; //! Flag for used constraints

  Int_t* fGlo2CGlo; ///< [fNGloPar] global ID to compressed ID buffer
  Int_t* fCGlo2Glo; ///< [fNGloPar] compressed ID to global ID buffer

  // Matrices
  o2::mft::SymMatrix* fMatCLoc;     ///< Matrix C local
  o2::mft::MatrixSq* fMatCGlo;      ///< Matrix C global
  o2::mft::RectMatrix* fMatCGloLoc; ///< Rectangular matrix C g*l
  Int_t* fFillIndex;                ///< [fNGloPar] auxilary index array for fast matrix fill
  Double_t* fFillValue;             ///< [fNGloPar] auxilary value array for fast matrix fill

  // processed data record bufferization
  TString fRecDataTreeName;   ///< Name of data records tree
  TString fRecConsTreeName;   ///< Name of constraints records tree
  TString fRecDataBranchName; ///< Name of data records branch name
  TString fRecConsBranchName; ///< Name of constraints records branch name

  TFile* fRecChi2File;
  TString fRecChi2FName;
  TString fRecChi2TreeName; ///< Name of chi2 per record tree
  TTree* fTreeChi2;
  float fSumChi2;
  bool fIsChi2BelowLimit;
  int fRecNDoF;

  TString fDataRecFName;    ///< Name of File for data records
  MillePedeRecord* fRecord; ///< Buffer of measurements records
  TFile* fDataRecFile;      ///< File of processed measurements records
  TTree* fTreeData;         ///< Tree of processed measurements records
  Int_t fRecFileStatus;     ///< state of the record file (0-no, 1-read, 2-rw)

  TString fConstrRecFName; ///< Name of File for constraints records
  TTree* fTreeConstr;      //! Tree of constraint records
  TFile* fConsRecFile;     //! File of processed constraints records
  Long_t fCurrRecDataID;   ///< ID of the current data record
  Long_t fCurrRecConstrID; ///< ID of the current constraint record
  Bool_t fLocFitAdd;       ///< Add contribution of carrent track (and not eliminate it)
  Bool_t fUseRecordWeight; ///< force or ignore the record weight
  Int_t fMinRecordLength;  ///< ignore shorter records
  Int_t fSelFirst;         ///< event selection start
  Int_t fSelLast;          ///< event selection end
  TArrayL* fRejRunList;    ///< list of runs to reject (if any)
  TArrayL* fAccRunList;    ///< list of runs to select (if any)
  TArrayF* fAccRunListWgh; ///< optional weights for data of accepted runs (if any)
  Double_t fRunWgh;        ///< run weight
  Double_t fWghScl[2];     ///< optional rescaling for odd/even residual weights (see its usage in LocalFit)
  const Int_t* fkReGroup;  ///< optional regrouping of parameters wrt ID's from the records

  static Bool_t fgInvChol;        ///< Invert global matrix in Cholesky solver
  static Bool_t fgWeightSigma;    ///< weight parameter constraint by statistics
  static Bool_t fgIsMatGloSparse; ///< Type of the global matrix (sparse ...)
  static Int_t fgMinResCondType;  ///< Type of the preconditioner for MinRes method
  static Double_t fgMinResTol;    ///< Tolerance for MinRes solution
  static Int_t fgMinResMaxIter;   ///< Max number of iterations for the MinRes method
  static Int_t fgIterSol;         ///< type of iterative solution: MinRes or FGMRES
  static Int_t fgNKrylovV;        ///< size of Krylov vectors buffer in FGMRES

  ClassDef(MillePede2, 1);
};

//_____________________________________________________________________________
inline void MillePede2::ReadRecordData(Long_t recID)
{
  fTreeData->GetEntry(recID);
  fCurrRecDataID = recID;
}

//_____________________________________________________________________________
inline void MillePede2::ReadRecordConstraint(Long_t recID)
{
  fTreeConstr->GetEntry(recID);
  fCurrRecConstrID = recID;
}

//_____________________________________________________________________________
inline void MillePede2::SaveRecordData()
{
  fTreeData->Fill();
  fRecord->Reset();
  fCurrRecDataID++;
}

//_____________________________________________________________________________
inline void MillePede2::SaveRecordConstraint()
{
  fTreeConstr->Fill();
  fRecord->Reset();
  fCurrRecConstrID++;
}

} // namespace mft
} // namespace o2

#endif
