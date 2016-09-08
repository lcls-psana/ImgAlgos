#ifndef IMGALGOS_ALGUTILS_H
#define IMGALGOS_ALGUTILS_H

//--------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description: see documentation below
//------------------------------------------------------------------------

#include <string>
#include <vector>
#include <iostream> // for cout, ostream
#include <cstddef>  // for size_t
#include <cstring>  // for memcpy
#include <cmath>    // for sqrt

//#include "PSCalib/CalibPars.h"  // for pixel_mask_t
#include "MsgLogger/MsgLogger.h"
#include "ndarray/ndarray.h"
#include "ImgAlgos/GlobalMethods.h" // TwoIndexes
//#include "ImgAlgos/Window.h"

using namespace std;

//------------------------------------

namespace ImgAlgos {
/// @addtogroup ImgAlgos

/**
 *  @ingroup ImgAlgos
 *
 *  @brief AlgImgProc - class for 2-d image processing algorithms.
 *
 *  This class is not suppose to be used separately. 
 *  Class AlgImgProc is a part of the Python-C++ algorithm inteface project.
 *
 *
 *  This software was developed for the LCLS project.  If you use all or 
 *  part of it, please give an appropriate acknowledgment.
 *
 *  @version $Id$
 *
 *  @author Mikhail S. Dubrovin
 *
 *  @see AlgArrProc, pyImgAlgos.cpp, PyAlgos.py
 *
 *
 *  @anchor interface
 *  @par<interface> Interface Description
 *
 * 
 *  @li  Includes and typedefs
 *  @code
 *  #include "ImgAlgos/AlgUtils.h"
 *  #include "ndarray/ndarray.h"     // need it for I/O arrays
 *  
 *  typedef ImgAlgos::AlgImgProc::conmap_t conmap_t;
 *  @endcode
 *
 *
 *  @li Initialization
 *  \n
 *  @code
 *  size_t      seg    = 2;
 *  size_t      rowmin = 10;
 *  size_t      rowmax = 170;
 *  size_t      colmin = 100;
 *  size_t      colmax = 200;
 *  unsigned    pbits  = 0;  // 0-print nothing, 2-input parameters and S/N matrix of indexes, 512-tracking.
 * 
 *  ImgAlgos::AlgImgProc* aip = new ImgAlgos::AlgImgProc(seg, rowmin, rowmax, colmin, colmax, pbits);
 *  @endcode
 *
 *
 *  @li Define input parameters
 *  @code
 *  ndarray<const T,2> data = ....;    // calibrated data ndarray
 *  ndarray<mask_t,2>  mask = ....;    // mask ndarray, may be omitted
 *  ndarray<mask_t,2>  son;            // output S/N ndarray

 *  unsigned rank = 4;
 *  float    r0   = 5;
 *  float    dr   = 0.05;
 *  ...
 *  @endcode
 *
 *
 *  @li Set methods
 *  @code
 *  aip->setSoNPars(r0,dr);
 *  aip->setWindows(winds);
 *  aip->setPeakSelectionPars(npix_min, npix_max, amax_thr, atot_thr, son_min);
 *  @endcode
 *
 *
 *  @li Get methods
 *  @code
 *   size_t ind = aip->segind()
 *   size_t counter = aip -> numberOfPixAboveThr<T>(seg_data, seg_mask, thr);
 *   double intensity = aip -> intensityOfPixAboveThr<T>(seg_data, seg_mask, thr);
 *   std::vector<Peak>& peaks = aip -> peakFinderV1<T>(seg_data, seg_mask, thr_low, thr_high, rad, dr);
 *
 *   // The same peak-finders after revision-1
 *   std::vector<Peak>& peaks = aip -> peakFinderV2r1<T>(seg_data, seg_mask, thr, r0, dr);
 *  @endcode
 *
 *
 *  @li Print methods
 *  @code
 *  aip->printInputPars();
 *
 *  Peak& peak = ...
 *  cout << peak ...
 *  @endcode
 */

//  typedef unsigned shape_t;
//  typedef PSCalib::CalibPars::pixel_mask_t mask_t;
//  typedef uint32_t conmap_t;
//  typedef uint16_t pixel_status_t;
//  typedef uint16_t pixel_maximums_t;
//  typedef uint16_t pixel_minimums_t;
//  typedef float son_t;
//  typedef uint16_t nphoton_t;
//  typedef float    fphoton_t;


//--------------------
//--------------------
/// Returns string name of the class for messanger

inline const char* _name() {return "ImgAlgos::AlgUtils";}

//--------------------
//--------------------
/**
 * @brief mapOfPhotonNumbersV1 - Chuck's photon counting algorithm - apply fancy correction for split photons.
 * 
 * 1) splits calibrated data for uint (floor) and float leftover fractional number of photons
 * 2) merge fractional number of photons to largest intensity integer  
 * 3) sum together uint and merged fractional maps
 * 
 * Returns array with (uint16) number of photons per pixel from input array of calibrated intensities.
 * param[in]  data - ndarray with calibrated intensities
 * param[in]  mask - ndarray with mask of bad/good (0/1) pixels
 */

//template <typename T>
//ndarray<nphoton_t, 2>& 
//mapOfPhotonNumbersV1( const ndarray<const T,2>&      data
//                    , const ndarray<const mask_t,2>& mask
//                    )
//{
  //MsgLog(_name(), info, "in mapOfPhotonNumbersV1, seg=" << m_seg << "\n    in window: " << m_win);
  //_splitDataForUintAndFloat<T>(data, mask);

  //_makeMapOfLocalMaximums<T>(data, mask, rank);
  //_makePeaksFromMapOfLocalMaximums<T>(data, mask, rank);
  //_addSoNToPeaks<T>(data, mask, r0, dr);
  //_makeVectorOfSelectedPeaks();

  //return m_nphoton; 
//}

//--------------------
void test_print_1()
{
  MsgLog(_name(), info, "in test_print_1 XXX TEST");
}

//--------------------
//--------------------
//--------------------
//--------------------
//--------------------
//--------------------
//--------------------
//--------------------

} // namespace ImgAlgos

#endif // IMGALGOS_ALGUTILS_H
