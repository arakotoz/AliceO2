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

///
/// \file    Options.cxx
/// \author  Julian Myrcha

#ifndef ALICE_O2_EVENTVISUALISATION_VIEW_OPTIONS_H
#define ALICE_O2_EVENTVISUALISATION_VIEW_OPTIONS_H

#include <string>

namespace o2
{
namespace event_visualisation
{

class Options
{
 private:
  // stored options
  bool mJSON;                    // -j
  bool mOnline;                  // -o (must specify -d!)
  bool mRandomTracks;            // -r
  std::string mFileName;         // -f 'data.root'
  std::string mDataFolder;       // -d './'
  std::string mSavedDataFolder;  // -s './'
  std::string mImageFolder;      // -i './'
  long mMemoryLimit;             // -m 1500 (MB) = 1.5GB
  bool mHideDplGUI;              // -hg
  bool mGenerateScreenshots;     // -gs
  bool mSequential;
  std::string mAODConverterPath; // -a 'o2-eve-aodconverter'

  // helper methods
  static Options instance;
  bool saveToJSON(std::string filename);   // stores options to current folder
  bool readFromJSON(std::string filename); // read options from option file

 public:
  static Options* Instance() { return &instance; }
  std::string printOptions();
  bool processCommandLine(int argc, char* argv[]);

  // get access methods
  bool sequential() const { return this->mSequential; }
  bool json() const { return this->mJSON; }
  bool online() const { return this->mOnline; }
  std::string dataFolder() const { return this->mDataFolder; }
  std::string imageFolder() const { return this->mImageFolder; }
  std::string savedDataFolder() const { return this->mSavedDataFolder; }
  std::string fileName() const { return this->mFileName; }
  bool randomTracks() const { return this->mRandomTracks; }
  long memoryLimit() const { return this->mMemoryLimit; }
  bool hideDplGUI() const { return this->mHideDplGUI; }
  bool generateScreenshots() const { return this->mGenerateScreenshots; }
  std::string AODConverterPath() { return this->mAODConverterPath; }
};

} // namespace event_visualisation
} // namespace o2

#endif // ALICE_O2_EVENTVISUALISATION_VIEW_OPTIONS_H
