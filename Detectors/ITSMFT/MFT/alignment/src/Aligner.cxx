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

/// @file Aligner.cxx

#include <iostream>
#include <sstream>
#include <string>

#include <TString.h>

#include "Framework/Logger.h"
#include "MFTAlignment/Aligner.h"
#include "MFTAlignment/MillePede2.h"

using namespace o2::mft;

ClassImp(o2::mft::Aligner);

//__________________________________________________________________________
Aligner::Aligner()
  : mStartFac(256),
    mChi2CutNStdDev(3),
    mResCutInitial(100.),
    mResCut(100.),
    mMilleRecordsFileName("mft_mille_records.root"),
    mMilleConstraintsRecFileName("mft_mille_constraints.root"),
    mMillepede(nullptr),
    mIsInitDone(false),
    mGlobalParameterStatus(nullptr)
{
  // default allowed variations w.r.t. global system coordinates
  mAllowVar[0] = 0.5;  // delta translation in x (cm)
  mAllowVar[1] = 0.5;  // delta translation in y (cm)
  mAllowVar[2] = 0.01; // rotation angle Rz around z-axis (rad)
  mAllowVar[3] = 0.5;  // delta translation in z (cm)

  mGlobalParameterStatus = (int*)malloc(sizeof(int) * mNumberOfGlobalParam);

  for (int iPar = 0; iPar < mNumberOfGlobalParam; iPar++) {
    mGlobalParameterStatus[iPar] = mFreeParId;
  }
  LOGF(info, "Aligner instantiated");
}

//__________________________________________________________________________
Aligner::~Aligner()
{
  free(mGlobalParameterStatus);
}

//__________________________________________________________________________
void Aligner::init()
{
  if (mIsInitDone)
    return;

  mMillepede = std::make_unique<MillePede2>();
  mMillepede->InitMille(mNumberOfGlobalParam,
                        mNumberOfTrackParam,
                        mChi2CutNStdDev,
                        mResCut,
                        mResCutInitial);

  LOG(info) << "-------------- Aligner configured with -----------------";
  LOGF(info, "Chi2CutNStdDev = %d", mChi2CutNStdDev);
  LOGF(info, "ResidualCutInitial = %.3f", mResCutInitial);
  LOGF(info, "ResidualCut = %.3f", mResCut);
  LOGF(info, "mStartFac = %.3f", mStartFac);
  LOGF(info,
       "Allowed variation: dx = %.3f, dy = %.3f, dz = %.3f, dRz = %.4f",
       mAllowVar[0], mAllowVar[1], mAllowVar[3], mAllowVar[2]);
  LOG(info) << "-----------------------------------------------------------";

  // set allowed variations for all parameters
  for (int chipId = 0; chipId < mNumberOfSensors; ++chipId) {
    for (Int_t iPar = 0; iPar < mNDofPerSensor; ++iPar) {
      mMillepede->SetParSigma(chipId * mNDofPerSensor + iPar, mAllowVar[iPar]);
    }
  }

  // set iterations
  if (mStartFac > 1) {
    mMillepede->SetIterations(mStartFac);
  }

  LOGF(info, "Aligner init done");
  mIsInitDone = true;
}
