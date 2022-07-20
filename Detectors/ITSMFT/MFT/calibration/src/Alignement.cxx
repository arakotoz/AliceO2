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

/// @file Alignment.cxx

#include "Framework/InputSpec.h"
#include "Framework/Logger.h"
#include <Framework/InputRecord.h>
#include "MFTCalibration/AlignPointHelper.h"
#include "MFTCalibration/Alignment.h"
#include "MFTTracking/IOUtils.h"
#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"

using namespace o2::mft;

ClassImp(o2::mft::Alignment);

//__________________________________________________________________________
Alignment::Alignment()
  : mNumberTFs(0),
    mChi2CutNStdDev(3),
    mResCutInitial(100.),
    mResCut(100.)
{
  mMillepede = std::make_unique<MillePede2>();

  GeometryTGeo* geom = GeometryTGeo::Instance();
  geom->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                             o2::math_utils::TransformType::L2G));
  mAlignPoint = std::make_unique<AlignPointHelper>(geom);

  // default allowed variations w.r.t. global system coordinates
  mAllowVar[0] = 0.5;  // x (cm)
  mAllowVar[1] = 0.5;  // y (cm)
  mAllowVar[2] = 0.01; // rotation angle Rz around z-axis (rad)
  mAllowVar[3] = 0.5;  // z (cm)
}

//__________________________________________________________________________
void Alignment::init(TString dataRecFName, TString consRecFName)
{
  mMillepede->InitMille(mNumberOfGlobalParam,
                        mNumberOfTrackParam,
                        mChi2CutNStdDev,
                        mResCut,
                        mResCutInitial);
  // filenames for the processed data and constraints records
  mMillepede->SetDataRecFName(dataRecFName.Data());
  mMillepede->SetConsRecFName(consRecFName.Data());
}

//__________________________________________________________________________
void Alignment::processTimeFrame(o2::framework::ProcessingContext& ctx)
{
  mNumberTFs++; // TF Counter

  // get tracks
  mMFTTracks = ctx.inputs().get<gsl::span<o2::mft::TrackMFT>>("tracks");
  mMFTTracksROF = ctx.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("tracksrofs");
  mMFTTrackClusIdx = ctx.inputs().get<gsl::span<int>>("trackClIdx");

  // get clusters
  mMFTClusters = ctx.inputs().get<gsl::span<o2::itsmft::CompClusterExt>>("compClusters");
  mMFTClustersROF = ctx.inputs().get<gsl::span<o2::itsmft::ROFRecord>>("clustersrofs");
  mMFTClusterPatterns = ctx.inputs().get<gsl::span<unsigned char>>("patterns");
  pattIt = mMFTClusterPatterns.begin();
  mMFTClustersGlobal.clear();
  mMFTClustersGlobal.reserve(mMFTClusters.size());
  o2::mft::ioutils::convertCompactClusters(
    mMFTClusters, pattIt, mMFTClustersGlobal, mDictionary);
}

//__________________________________________________________________________
void Alignment::processRecoTracks()
{
  for (auto& mftTrack : mMFTTracks) {

    resetLocalDerivative();
    resetGlocalDerivative();
    mTrackRecord.Reset();
  }
}

//__________________________________________________________________________
void Alignment::setLocalDerivative(Int_t index, Double_t value)
{
  if (index < mNumberOfTrackParam) {
    mLocalDerivatives[index] = value;
  } else {
    LOG(error) << "Alignment::setLocalDerivative() - "
               << "index " << index
               << " >= " << mNumberOfTrackParam;
  }
}

//__________________________________________________________________________
void Alignment::setGlobalDerivative(Int_t index, Double_t value)
{
  if (index < mNumberOfGlobalParam) {
    mGlobalDerivatives[index] = value;
  } else {
    LOG(error) << "Alignment::setGlobalDerivative() - "
               << "index " << index
               << " >= " << mNumberOfGlobalParam;
  }
}

//__________________________________________________________________________
void Alignment::resetLocalDerivative()
{
  for (int i = 0; i < mNumberOfTrackParam; ++i) {
    mLocalDerivatives[i] = 0.0;
  }
}

//__________________________________________________________________________
void Alignment::resetGlocalDerivative()
{
  for (int i = 0; i < mNumberOfGlobalParam; ++i) {
    mGlobalDerivatives[i] = 0.0;
  }
}