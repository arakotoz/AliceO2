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

/// \file GPUDisplayBackendWindows.h
/// \author David Rohr

#ifndef GPUDISPLAYBACKENDWINDOWS_H
#define GPUDISPLAYBACKENDWINDOWS_H

#include "GPUDisplayBackend.h"

namespace GPUCA_NAMESPACE::gpu
{
class GPUDisplayBackendWindows : public GPUDisplayBackend
{
 public:
  GPUDisplayBackendWindows() = default;
  ~GPUDisplayBackendWindows() override = default;

  int StartDisplay() override;
  void DisplayExit() override;
  void SwitchFullscreen(bool set) override;
  void ToggleMaximized(bool set) override;
  void SetVSync(bool enable) override;
  void OpenGLPrint(const char* s, float x, float y, float r, float g, float b, float a, bool fromBotton = true) override;

 private:
  int OpenGLMain() override;
};
} // namespace GPUCA_NAMESPACE::gpu

#endif
