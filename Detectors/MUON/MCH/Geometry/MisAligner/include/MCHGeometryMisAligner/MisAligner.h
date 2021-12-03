// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/** @file MisAligner.h
 * Generate misalignments.
 * @author Javier Castillo Castellanos
 */

#ifndef O2_MCH_GEOMETRY_MIS_ALIGNER
#define O2_MCH_GEOMETRY_MIS_ALIGNER

class TGeoCombiTrans;
class TClonesArray;

namespace o2
{
namespace mch
{
namespace geo
{

class MisAligner
{
 public:
  MisAligner(double cartXMisAligM, double cartXMisAligW, double cartYMisAligM, double cartYMisAligW, double angMisAligM, double angMisAligW);
  MisAligner(double cartMisAligM, double cartMisAligW, double angMisAligM, double angMisAligW);
  MisAligner(double cartMisAligW, double angMisAligW);
  MisAligner();
  virtual ~MisAligner();

  //_________________________________________________________________
  // methods

  // return a misaligned geometry obtained from the existing one.
  void MisAlign(bool verbose = false);

  /// Set cartesian displacement parameters different along x, y
  void SetCartMisAlig(double xmean, double xwidth, double ymean, double ywidth, double zmean = 0., double zwidth = 0.)
  {
    fDetElemMisAlig[0][0] = xmean;
    fDetElemMisAlig[0][1] = xwidth;
    fDetElemMisAlig[1][0] = ymean;
    fDetElemMisAlig[1][1] = ywidth;
    fDetElemMisAlig[2][0] = zmean;
    fDetElemMisAlig[2][1] = zwidth;
  }

  /// Set cartesian displacement parameters, the same along x, y
  void SetCartMisAlig(double mean, double width)
  {
    fDetElemMisAlig[0][0] = mean;
    fDetElemMisAlig[0][1] = width;
    fDetElemMisAlig[1][0] = mean;
    fDetElemMisAlig[1][1] = width;
  }

  /// Set angular displacement
  void SetAngMisAlig(double zmean, double zwidth, double xmean = 0., double xwidth = 0., double ymean = 0., double ywidth = 0.)
  {
    fDetElemMisAlig[3][0] = xmean;
    fDetElemMisAlig[3][1] = xwidth;
    fDetElemMisAlig[4][0] = ymean;
    fDetElemMisAlig[4][1] = ywidth;
    fDetElemMisAlig[5][0] = zmean;
    fDetElemMisAlig[5][1] = zwidth;
  }

  void SetXYAngMisAligFactor(double factor);

  void SetZCartMisAligFactor(double factor);

  /// Set option for gaussian distribution
  void SetUseGaus(bool usegaus)
  {
    fUseGaus = usegaus;
    fUseUni = !usegaus;
  }

  /// Set option for uniform distribution
  void SetUseUni(bool useuni)
  {
    fUseGaus = !useuni;
    fUseUni = useuni;
  }

  /// Set module (half chambers) cartesian displacement parameters
  void SetModuleCartMisAlig(double xmean, double xwidth, double ymean, double ywidth, double zmean, double zwidth)
  {
    fModuleMisAlig[0][0] = xmean;
    fModuleMisAlig[0][1] = xwidth;
    fModuleMisAlig[1][0] = ymean;
    fModuleMisAlig[1][1] = ywidth;
    fModuleMisAlig[2][0] = zmean;
    fModuleMisAlig[2][1] = zwidth;
  }

  /// Set module (half chambers) cartesian displacement parameters
  void SetModuleAngMisAlig(double xmean, double xwidth, double ymean, double ywidth, double zmean, double zwidth)
  {
    fModuleMisAlig[3][0] = xmean;
    fModuleMisAlig[3][1] = xwidth;
    fModuleMisAlig[4][0] = ymean;
    fModuleMisAlig[4][1] = ywidth;
    fModuleMisAlig[5][0] = zmean;
    fModuleMisAlig[5][1] = zwidth;
  }

  /// Set alignment resolution to misalign objects to be stored in CDB
  void SetAlignmentResolution(const TClonesArray* misAlignArray, int chId = -1, double chResX = -1., double chResY = -1., double deResX = -1., double deResY = -1.);

 protected:
  /// Not implemented
  MisAligner(const MisAligner& right);
  /// Not implemented
  MisAligner& operator=(const MisAligner& right);

 private:
  bool matrixToAngles(const double* rot, double& psi, double& theta, double& phi);
  // return a misaligned transformation
  TGeoCombiTrans
    MisAlignDetElem() const;
  TGeoCombiTrans MisAlignModule() const;
  void GetUniMisAlign(double cartMisAlig[3], double angMisAlig[3], const double lParMisAlig[6][2]) const;
  void GetGausMisAlign(double cartMisAlig[3], double angMisAlig[3], const double lParMisAlig[6][2]) const;

  bool fUseUni;                 ///< use uniform distribution for misaligmnets
  bool fUseGaus;                ///< use gaussian distribution for misaligmnets
  double fDetElemMisAlig[6][2]; ///< Mean and width of the displacements of the detection elements along x,y,z (translations) and about x,y,z (rotations)
  double fModuleMisAlig[6][2];  ///< Mean and width of the displacements of the modules along x,y,z (translations) and about x,y,z (rotations)

  double fXYAngMisAligFactor; ///< factor (<1) to apply to angular misalignment range since range of motion is restricted out of the xy plane
  double fZCartMisAligFactor; ///< factor (<1) to apply to cartetian misalignment range since range of motion is restricted in z direction
};

} // namespace geo
} // namespace mch
} // namespace o2
#endif // GEOMETRY_MIS_ALIGNER_H