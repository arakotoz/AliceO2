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

/// @file AlignPointControl.cxx

#include "Framework/Logger.h"

#include "MFTAlignment/AlignPointControl.h"

using namespace o2::mft;

ClassImp(o2::mft::AlignPointControl);

//__________________________________________________________________________
AlignPointControl::AlignPointControl()
  : mControlTree(nullptr),
    mControlFile(nullptr),
    mIsSuccessfulInit(false),
    mNEntriesAutoSave(10000),
    mOutFileName("mft_align_point.root"),
    mTreeTitle("align point info tree")
{
  mPointInfo.sensor = 0;
  mPointInfo.layer = 0;
  mPointInfo.disk = 0;
  mPointInfo.half = 0;
  mPointInfo.measuredGlobalX = 0;
  mPointInfo.measuredGlobalY = 0;
  mPointInfo.measuredGlobalZ = 0;
  mPointInfo.measuredLocalX = 0;
  mPointInfo.measuredLocalY = 0;
  mPointInfo.measuredLocalZ = 0;
  mPointInfo.residualX = 0;
  mPointInfo.residualY = 0;
  mPointInfo.residualZ = 0;
  mPointInfo.residualLocalX = 0;
  mPointInfo.residualLocalY = 0;
  mPointInfo.residualLocalZ = 0;
  mPointInfo.recoGlobalX = 0;
  mPointInfo.recoGlobalY = 0;
  mPointInfo.recoGlobalZ = 0;
  mPointInfo.recoLocalX = 0;
  mPointInfo.recoLocalY = 0;
  mPointInfo.recoLocalZ = 0;
}

//__________________________________________________________________________
void AlignPointControl::setCyclicAutoSave(const long nEntries)
{
  if (nEntries <= 0)
    return;
  mNEntriesAutoSave = nEntries;
}

//__________________________________________________________________________
void AlignPointControl::init()
{
  mIsSuccessfulInit = true;

  if (mControlFile == nullptr)
    mControlFile = std::make_unique<TFile>(mOutFileName.Data(), "recreate", "", 505);

  if (mControlTree == nullptr) {
    mControlFile->cd();
    mControlTree = std::make_unique<TTree>("point", mTreeTitle.Data());
    mControlTree->SetAutoSave(mNEntriesAutoSave); // flush the TTree to disk every N entries
    mControlTree->Branch("sensor", &mPointInfo.sensor, "sensor/s");
    mControlTree->Branch("layer", &mPointInfo.layer, "layer/s");
    mControlTree->Branch("disk", &mPointInfo.disk, "disk/s");
    mControlTree->Branch("half", &mPointInfo.half, "half/s");
    mControlTree->Branch("measuredGlobalX", &mPointInfo.measuredGlobalX, "measuredGlobalX/D");
    mControlTree->Branch("measuredGlobalY", &mPointInfo.measuredGlobalY, "measuredGlobalY/D");
    mControlTree->Branch("measuredGlobalZ", &mPointInfo.measuredGlobalZ, "measuredGlobalZ/D");
    mControlTree->Branch("measuredLocalX", &mPointInfo.measuredLocalX, "measuredLocalX/D");
    mControlTree->Branch("measuredLocalY", &mPointInfo.measuredLocalY, "measuredLocalY/D");
    mControlTree->Branch("measuredLocalZ", &mPointInfo.measuredLocalZ, "measuredLocalZ/D");
    mControlTree->Branch("residualX", &mPointInfo.residualX, "residualX/D");
    mControlTree->Branch("residualY", &mPointInfo.residualY, "residualY/D");
    mControlTree->Branch("residualZ", &mPointInfo.residualZ, "residualZ/D");
    mControlTree->Branch("residualLocalX", &mPointInfo.residualLocalX, "residualLocalX/D");
    mControlTree->Branch("residualLocalY", &mPointInfo.residualLocalY, "residualLocalY/D");
    mControlTree->Branch("residualLocalZ", &mPointInfo.residualLocalZ, "residualLocalZ/D");
    mControlTree->Branch("recoGlobalX", &mPointInfo.recoGlobalX, "recoGlobalX/D");
    mControlTree->Branch("recoGlobalY", &mPointInfo.recoGlobalY, "recoGlobalY/D");
    mControlTree->Branch("recoGlobalZ", &mPointInfo.recoGlobalZ, "recoGlobalZ/D");
    mControlTree->Branch("recoLocalX", &mPointInfo.recoLocalX, "recoLocalX/D");
    mControlTree->Branch("recoLocalY", &mPointInfo.recoLocalY, "recoLocalY/D");
    mControlTree->Branch("recoLocalZ", &mPointInfo.recoLocalZ, "recoLocalZ/D");
  }
  if ((!mControlFile) || (mControlFile->IsZombie())) {
    mIsSuccessfulInit = false;
    LOG(error) << "AlignPointControl::init() - failed, no viable output file !";
  }
  if (!mControlTree) {
    mIsSuccessfulInit = false;
    LOG(error) << "AlignPointControl::init() - failed, no TTree !";
  }
}

//__________________________________________________________________________
void AlignPointControl::terminate()
{
  if (mControlFile && mControlFile->IsWritable() && mControlTree) {
    mControlFile->cd();
    mControlTree->Write();
    LOG(info) << "AlignPointControl::terminate() - wrote "
              << mTreeTitle.Data();
  }
  if (mControlFile) {
    mControlFile->Close();
    LOG(info) << "AlignPointControl::terminate() - Closed file "
              << mOutFileName.Data();
  }
}

//__________________________________________________________________________
void AlignPointControl::fill(
  std::shared_ptr<o2::mft::AlignPointHelper> alignPoint,
  const int iTrack,
  const bool doPrint)
{
  if (!isInitOk()) {
    LOG(warning) << "AlignPointControl::fill() - aborted, init was not ok !";
    return;
  }

  bool isPointok = setControlPoint(alignPoint);

  if (isPointok) {
    mControlTree->Fill();
    if (doPrint) {
      LOGF(info, "AlignPointControl::fillControlTree() - track %i h %d d %d l %d s %4d lMpos x %.2e y %.2e z %.2e gMpos x %.2e y %.2e z %.2e gRpos x %.2e y %.2e z %.2e",
           iTrack, mPointInfo.half, mPointInfo.disk, mPointInfo.layer, mPointInfo.sensor,
           mPointInfo.measuredLocalX, mPointInfo.measuredLocalY, mPointInfo.measuredLocalZ,
           mPointInfo.measuredGlobalX, mPointInfo.measuredGlobalY, mPointInfo.measuredGlobalZ,
           mPointInfo.recoGlobalX, mPointInfo.recoGlobalY, mPointInfo.recoGlobalZ);
    }
  }
}

//__________________________________________________________________________
bool AlignPointControl::setControlPoint(
  std::shared_ptr<o2::mft::AlignPointHelper> alignPoint)
{
  if (!alignPoint) {
    LOG(warning) << "AlignPointControl::setControlPoint() - aborted, can not use a null pointer";
    return false;
  }

  mPointInfo.sensor = alignPoint->getSensorId();
  mPointInfo.layer = alignPoint->layer();
  mPointInfo.disk = alignPoint->disk();
  mPointInfo.half = alignPoint->half();
  mPointInfo.measuredGlobalX = alignPoint->getGlobalMeasuredPosition().X();
  mPointInfo.measuredGlobalY = alignPoint->getGlobalMeasuredPosition().Y();
  mPointInfo.measuredGlobalZ = alignPoint->getGlobalMeasuredPosition().Z();
  mPointInfo.measuredLocalX = alignPoint->getLocalMeasuredPosition().X();
  mPointInfo.measuredLocalY = alignPoint->getLocalMeasuredPosition().Y();
  mPointInfo.measuredLocalZ = alignPoint->getLocalMeasuredPosition().Z();
  mPointInfo.residualX = alignPoint->getGlobalResidual().X();
  mPointInfo.residualY = alignPoint->getGlobalResidual().Y();
  mPointInfo.residualZ = alignPoint->getGlobalResidual().Z();
  mPointInfo.residualLocalX = alignPoint->getLocalResidual().X();
  mPointInfo.residualLocalY = alignPoint->getLocalResidual().Y();
  mPointInfo.residualLocalZ = alignPoint->getLocalResidual().Z();
  mPointInfo.recoGlobalX = alignPoint->getGlobalRecoPosition().X();
  mPointInfo.recoGlobalY = alignPoint->getGlobalRecoPosition().Y();
  mPointInfo.recoGlobalZ = alignPoint->getGlobalRecoPosition().Z();
  mPointInfo.recoLocalX = alignPoint->getLocalRecoPosition().X();
  mPointInfo.recoLocalY = alignPoint->getLocalRecoPosition().Y();
  mPointInfo.recoLocalZ = alignPoint->getLocalRecoPosition().Z();

  return true;
}
