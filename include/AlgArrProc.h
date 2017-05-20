#ifndef IMGALGOS_ALGARRPROC_H
#define IMGALGOS_ALGARRPROC_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description: see documentation below
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------

#include <string>
#include <vector>
#include <iostream> // for cout
#include <stdint.h> // uint8_t, uint32_t, etc.
#include <cstddef>  // for size_t
#include <cstring>  // for memcpy
//#include <math.h>   // for exp
#include <cmath>    // for sqrt
#include <stdexcept>
#include <algorithm> // for fill_n

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "PSCalib/CalibPars.h"  // for pixel_mask_t
#include "MsgLogger/MsgLogger.h"
#include "ndarray/ndarray.h"
#include "ImgAlgos/GlobalMethods.h"
#include "ImgAlgos/Window.h"

#include "ImgAlgos/AlgImgProc.h"


using namespace std;


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

namespace ImgAlgos {

/// @addtogroup ImgAlgos

/**
 *  @ingroup ImgAlgos
 *
 *  @brief AlgArrProc - class for 2-d image processing algorithms.
 *
 *  This class is used in pyImgAlgos.cpp
 *  Class AlgArrProc is a part of the Python-C++ algorithm inteface project.
 *
 *  This software was developed for the LCLS project.  If you use all or 
 *  part of it, please give an appropriate acknowledgment.
 *
 *  @version $Id$
 *
 *  @author Mikhail S. Dubrovin
 *
 *  @see pyImgAlgos.cpp, PyAlgos.py
 *
 *
 *  @anchor interface
 *  @par<interface> Interface Description
 *
 * 
 *  @li Includes and typedefs
 *  @code
 *  #include "ImgAlgos/AlgArrProc.h"
 *  #include "ndarray/ndarray.h"     // need it for I/O arrays
 *  
 *  typedef PSCalib::CalibPars::pixel_mask_t mask_t;
 *  typedef ImgAlgos::AlgArrProc::wind_t wind_t;
 *  typedef ImgAlgos::AlgImgProc::conmap_t conmap_t;
 *  typedef ImgAlgos::AlgImgProc::pixel_status_t pixel_status_t;
 *  typedef ImgAlgos::AlgImgProc::pixel_minimums_t pixel_minimums_t;
 *  @endcode
 *
 *
 *  @li Initialization
 *  \n
 *  @code
 *  ndarray<const wind_t,2> winds = ...
 *  unsigned    pbits  = 0; // 0-print nothing, 1-list of peaks, 2-input parameters, 258-tracking.
 * 
 *  ImgAlgos::AlgArrProc* alg = new ImgAlgos::AlgArrProc(winds, pbits);
 *  @endcode
 *
 *
 *  @li Define input parameters
 *  @code
 *  ndarray<const T,3> data = ....;    // calibrated data ndarray
 *  ndarray<mask_t,3>  mask = ....;    // mask ndarray, may be omitted
 *  ndarray<mask_t,3>  son;            // output S/N ndarray
 *  ndarray<const wind_t,2> winds = ((0,  0, 185,   0, 388), \
 *                                   (1, 10, 103,  10, 204), \
 *                                   (1, 10, 103, 250, 380));
 *  size_t rank = 4;
 *  float  r0   = 5;
 *  float  dr   = 0.05;
 *  @endcode
 *
 *
 *  @li Set methods
 *  @code
 *  alg->setWindows(winds);
 *  alg->setPeakSelectionPars(npix_min, npix_max, amax_thr, atot_thr, son_min);
 *  alg->setSoNPars(r0,dr);
 *  @endcode
 *
 *
 *  @li Get methods
 *  @code
 *  unsigned npix = alg->numberOfPixAboveThr<T>(data, mask, thr);
 *  double intensity = alg->intensityOfPixAboveThr<T>(data, mask, thr);
 *  ndarray<const float, 3> peaks = alg->peakFinderV1<T>(data, mask, thr_low, thr_high, radius, dr);
 *  ndarray<const float, 3> peaks = alg->peakFinderV2<T>(data, mask, thr, r0, dr);
 *  ndarray<const float, 3> peaks = alg->peakFinderV3<T>(data, mask, rank, r0, dr);
 *  ndarray<const float, 3> peaks = alg->peakFinderV4<T>(data, mask, thr_low, thr_high, rank, r0, dr);
 *  
 *  // The same peak-finders after revision
 *  ndarray<const float, 3> peaks = alg->peakFinderV2r1<T>(data, mask, thr, r0, dr);
 *  ndarray<const float, 3> peaks = alg->peakFinderV3r1<T>(data, mask, rank, r0, dr, nsigm);
 *  ndarray<const float, 3> peaks = alg->peakFinderV3r2<T>(data, mask, rank, r0, dr, nsigm); // New revision
 *  ndarray<const float, 3> peaks = alg->peakFinderV4r1<T>(data, mask, thr_low, thr_high, rank, r0, dr);
 *  ndarray<const float, 3> peaks = alg->peakFinderV4r2<T>(data, mask, thr_low, thr_high, rank, r0, dr); // New revision
 *  
 *  // Call after peakFinderV2(...) ONLY!
 *  ndarray<conmap_t, 3> maps = alg->mapsOfConnectedPixels();
 *  
 *  // Returns map of pixel status after peakFinderV2, V4r2
 *  ndarray<pixel_status_t, 3> maps = alg->mapsOfPixelStatus();
 *  
 *  // Chuck's algorithm of photon counting (apply correction on split between pixels photons)
 *  ndarray<const nphoton_t, 3> maps = alg->mapsOfPhotonNumbers<T>(data, mask);
 *  @endcode
 *
 *
 *  @li Print methods
 *  @code
 *  alg->printInputPars();
 *  @endcode
 */

/*
struct SegWindow {
  size_t segind;
  Window window;
};
*/

class AlgArrProc {
public:

  typedef PSCalib::CalibPars::pixel_mask_t mask_t;
  typedef uint32_t wind_t;
  typedef AlgImgProc::pixel_status_t pixel_status_t;
  typedef AlgImgProc::conmap_t conmap_t;
  typedef AlgImgProc::pixel_minimums_t pixel_minimums_t;
  typedef AlgImgProc::pixel_maximums_t pixel_maximums_t;
  typedef AlgImgProc::nphoton_t nphoton_t;

  /*
  typedef unsigned shape_t;
  typedef uint32_t conmap_t;
  typedef uint16_t pixel_status_t;
  typedef float son_t;
  */


  /**
   * @brief Default constructor
   * 
   * @param[in] pbits     - print control bit-word; =0-print nothing, +1-input parameters, +2-algorithm area, etc.
   */
  AlgArrProc(const unsigned& pbits=0);

  /**
   * @brief Main constructor
   * 
   * @param[in] winds - ndarray with windows
   * @param[in] pbits     - print control bit-word; =0-print nothing, +1-input parameters, +2-algorithm area, etc.
   */
  AlgArrProc(ndarray<const wind_t,2> winds, const unsigned& pbits=0);

  /// Set windows
  void setWindows(ndarray<const wind_t,2> nda_winds);

  /// Set peak selection parameters MUST BE CALLED BEFORE PEAKFINDER!!!
  void setPeakSelectionPars(const float& npix_min=0, 
                            const float& npix_max=1e6, 
                            const float& amax_thr=0, 
                            const float& atot_thr=0,
                            const float& son_min =0);

  /**
   * @brief Set parameters for SoN (S/N) evaluation
   * 
   * @param[in] r0 - ring internal radius for S/N evaluation
   * @param[in] dr - ring width for S/N evaluation
   */
  void setSoNPars(const float& r0=5, const float& dr=0.05);

  /// Returns 3-d array of maps of pixel statuss for all segments, works after peakFinderV2, V4r2 ONLY!
  ndarray<const pixel_status_t, 3> mapsOfPixelStatus();
  /// Returns 3-d array of maps of connected pixels for all segments, works after peakFinderV2 ONLY!
  ndarray<const conmap_t, 3> mapsOfConnectedPixels();
  ndarray<const pixel_minimums_t, 3> mapsOfLocalMinimums();
  ndarray<const pixel_maximums_t, 3> mapsOfLocalMaximums();

  /// Destructor
  virtual ~AlgArrProc (); 

  /// Prints memeber data
  void printInputPars();

private:

  unsigned   m_pbits;            // pirnt control bit-word

  DATA_TYPE  m_dtype;            // numerated data type for data array
  unsigned   m_ndim;             // ndarray number of dimensions
  unsigned   m_nsegs;
  unsigned   m_nrows;
  unsigned   m_ncols;
  size_t     m_stride;
  unsigned   m_sshape[2];
  bool       m_is_inited;

  mask_t*    m_mask_def;
  const mask_t* m_mask;

  float    m_r0;       // radial parameter of the ring for S/N evaluation algorithm
  float    m_dr;       // ring width for S/N evaluation algorithm 

  float    m_peak_npix_min; // peak selection parameter
  float    m_peak_npix_max; // peak selection parameter
  float    m_peak_amax_thr; // peak selection parameter
  float    m_peak_atot_thr; // peak selection parameter
  float    m_peak_son_min;  // peak selection parameter

  std::vector<Window>        v_winds;
  std::vector<AlgImgProc*>   v_algip;    // vector of pointers to the AlgImgProc objects for windows

  /// Returns string name of the class for messanger
  inline const char* _name() {return "ImgAlgos::AlgArrProc";}

  /// Returns ndarray of peak-float pars evaluated by any of peakfinders 
  const ndarray<const float, 2> _ndarrayOfPeakPars(const unsigned& npeaks);

  //AlgArrProc ( const AlgArrProc& ) ;
  //AlgArrProc& operator = ( const AlgArrProc& ) ;


//--------------------

public:

//--------------------

  template <typename T, unsigned NDim>
  bool
  _initMask(const ndarray<const T, NDim>& data, const ndarray<const mask_t, NDim>& mask)
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in _initMask: mask.size = " << mask.size()
                                            << " data.size() = " << data.size());

    if(data.empty()) return false;

    if(mask.size()) {
        m_mask = mask.data();
	if(m_pbits & 256) MsgLog(_name(), info, "Mask is used for pixel processing.");
    } else {
      // Define default mask
      if(m_mask_def) delete m_mask_def;    
      m_mask_def = new mask_t[data.size()];
      std::fill_n(m_mask_def, int(data.size()), mask_t(1));
      m_mask = m_mask_def;

      if(m_pbits & 256) MsgLog(_name(), info, "Mask is empty, all pixels will be processed.");
    }

    m_ndim = NDim;
    if(m_ndim < 2) throw std::runtime_error("Non-acceptable number of dimensions < 2 in input ndarray");

    m_dtype  = dataType<T>();
    m_ncols  = data.shape()[m_ndim-1];
    m_nrows  = data.shape()[m_ndim-2];
    m_nsegs  = (m_ndim>2) ? data.size()/m_ncols/m_nrows : 1;
    m_stride = m_ncols*m_nrows;

    m_sshape[0] = m_nrows;
    m_sshape[1] = m_ncols;

    return true;
  }

//--------------------

  template <typename T, unsigned NDim>
  bool
  _initAlgImgProc(const ndarray<const T, NDim>& data, const ndarray<const mask_t, NDim>& mask)
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in _initAlgImgProc: mask.size = " << mask.size()
                                            << " data.size() = " << data.size());

    if(not _initMask(data, mask)) return false;

    if(m_is_inited) return true;

    m_is_inited = true;

    if(v_winds.empty()) {
        // ALL segments will be processed
        if(m_pbits & 256) MsgLog(_name(), info, "List of windows is empty, all sensors will be processed, number of windows = " << m_nsegs);
        v_algip.reserve(m_nsegs);
      
        for(size_t seg=0; seg<m_nsegs; ++seg) {            
            AlgImgProc* p_alg = new AlgImgProc(seg, 0, m_nrows, 0, m_ncols, m_pbits);
            v_algip.push_back(p_alg);
            p_alg->setSoNPars(m_r0, m_dr);
            p_alg->setPeakSelectionPars(m_peak_npix_min, m_peak_npix_max, m_peak_amax_thr, m_peak_atot_thr, m_peak_son_min);
        }       
    }
    else {
        // Windows ONLY will be processed
        if(m_pbits & 256) MsgLog(_name(), info, "Windows from the list will be processed, number of windows = " << v_winds.size());
	v_algip.reserve(v_winds.size());

        for(std::vector<Window>::iterator it = v_winds.begin(); it != v_winds.end(); ++ it) {
            AlgImgProc* p_alg = new AlgImgProc(it->segind, it->rowmin, it->rowmax, it->colmin, it->colmax , m_pbits);
            v_algip.push_back(p_alg); 
            p_alg->setSoNPars(m_r0, m_dr);
            p_alg->setPeakSelectionPars(m_peak_npix_min, m_peak_npix_max, m_peak_amax_thr, m_peak_atot_thr, m_peak_son_min);
        }
    }

    return true;
  }

//--------------------
//--------------------

  template <typename T, unsigned NDim>
  unsigned
  numberOfPixAboveThr( const ndarray<const T, NDim> data
		     , const ndarray<const mask_t, NDim> mask
                     , const T& thr
                     )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in numberOfPixAboveThr " << thr);

    if(! _initAlgImgProc<T,NDim>(data, mask)) return 0;

    unsigned counter = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        counter += (*it) -> numberOfPixAboveThr<T>(seg_data, seg_mask, thr);
    }
    return counter;
  }

//--------------------
//--------------------

  template <typename T, unsigned NDim>
  double
  intensityOfPixAboveThr( const ndarray<const T, NDim> data
                        , const ndarray<const mask_t, NDim> mask
                        , const T& thr
                        )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in intensityOfPixAboveThr" << thr);

    if(! _initAlgImgProc<T,NDim>(data, mask)) return 0;

    double intensity = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        intensity += (*it) -> intensityOfPixAboveThr<T>(seg_data, seg_mask, thr);
    }
    return intensity;
  }

//--------------------
//--------------------
  /// peakFinderV1 - "Droplet-finder" - two-threshold peak finding algorithm in the region defined by the radial parameter.

  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV1( const ndarray<const T, NDim> data
               , const ndarray<const mask_t, NDim> mask
               , const T& thr_low
               , const T& thr_high
               , const unsigned& rad=5
               , const float& dr=0.05
               )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV1");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV1<T>(seg_data, seg_mask, thr_low, thr_high, rad, dr);
	npeaks += peaks.size();
    }

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------

  /// peakFinderV4 the same as V1, but has rank and r0 parameters in stead of single radius. 
  /// peakFinderV4 - "Droplet-finder" - two-threshold peak finding algorithm in the region defined by the radial parameter.

  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV4( const ndarray<const T, NDim> data
              , const ndarray<const mask_t, NDim> mask
              , const T& thr_low
              , const T& thr_high
              , const unsigned& rank=4
              , const float& r0=5
              , const float& dr=0.05
              )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV4");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV4<T>(seg_data, seg_mask, thr_low, thr_high, rank, r0, dr);
	npeaks += peaks.size();
    }

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------

  /// peakFinderV4r1 - "Droplet-finder" - the same as V4, but returns background-corrected amp_max, amp_tot, son (total).
  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV4r1( const ndarray<const T, NDim> data
                , const ndarray<const mask_t, NDim> mask
                , const T& thr_low
                , const T& thr_high
                , const unsigned& rank=4
                , const float& r0=5
                , const float& dr=0.05
                )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV4r1");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV4r1<T>(seg_data, seg_mask, thr_low, thr_high, rank, r0, dr);
	npeaks += peaks.size();
    }

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------

  /// peakFinderV4r2 - "Droplet-finder" - further development of V4r1;
  ///                - defines droplet for connected pixels in the constrained region of rank
  ///                - keeps the same idea of droplet definition, but implementation has changed significantly
  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV4r2( const ndarray<const T, NDim> data
                , const ndarray<const mask_t, NDim> mask
                , const T& thr_low
                , const T& thr_high
                , const unsigned& rank=5
                , const float& r0=7
                , const float& dr=2
                )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV4r2");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV4r2<T>(seg_data, seg_mask, thr_low, thr_high, rank, r0, dr);
	npeaks += peaks.size();
    }

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------
  /// peakFinderV2 - "Flood filling" - makes a list of peaks for groups of connected pixels above threshold.

  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV2( const ndarray<const T, NDim> data
              , const ndarray<const mask_t, NDim> mask
              , const T& thr
	      , const float& r0=5
              , const float& dr=0.05
              )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV2");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV2<T>(seg_data, seg_mask, thr, r0, dr);
	npeaks += peaks.size();
    }

    if(m_pbits & 256) MsgLog(_name(), info, "total number of peaks=" << npeaks);    

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------
  /// peakFinderV2r1 - "Flood filling" - the same as V2, but returns background-corrected amp_max, amp_tot, son (total).

  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV2r1( const ndarray<const T, NDim> data
                , const ndarray<const mask_t, NDim> mask
                , const T& thr
	        , const float& r0=5
                , const float& dr=0.05
                )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV2r1");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV2r1<T>(seg_data, seg_mask, thr, r0, dr);
	npeaks += peaks.size();
    }

    if(m_pbits & 256) MsgLog(_name(), info, "total number of peaks=" << npeaks);    

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------
  /// peakFinderV3 - "Ranker" - makes a list of peaks for local maximums of requested rank.

  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV3( const ndarray<const T, NDim> data
              , const ndarray<const mask_t, NDim> mask
              , const size_t& rank = 2
	      , const float& r0=5
              , const float& dr=0.05
              )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV3");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV3<T>(seg_data, seg_mask, rank, r0, dr);
	npeaks += peaks.size();
    }

    if(m_pbits & 256) MsgLog(_name(), info, "total number of peaks=" << npeaks);    

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------
  /// peakFinderV3r1 - "Ranker" - the same as V3, but returns background-corrected amp_max, amp_tot, son (total).

  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV3r1( const ndarray<const T, NDim> data
                , const ndarray<const mask_t, NDim> mask
                , const size_t& rank = 5
	        , const float& r0=7
                , const float& dr=2
	        , const float& nsigm=0 // 0-turns off threshold algorithm, 1.64-leaves 5% of noise, etc.; 
                )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV3r1");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV3r1<T>(seg_data, seg_mask, rank, r0, dr, nsigm);
	npeaks += peaks.size();
    }

    if(m_pbits & 256) MsgLog(_name(), info, "total number of peaks=" << npeaks);    

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------
  /// peakFinderV3r2 - "Ranker" - the same as V3r1, but uses connected pixels only, order of algorithms has changed.

  template <typename T, unsigned NDim>
  ndarray<const float, 2>
  peakFinderV3r2( const ndarray<const T, NDim> data
                , const ndarray<const mask_t, NDim> mask
                , const size_t& rank = 5
	        , const float& r0=7
                , const float& dr=2
	        , const float& nsigm=0 // 0-turns off threshold algorithm, 1.64-leaves 5% of noise, etc.; 
                )
  {
    if(m_pbits & 256) MsgLog(_name(), info, "in peakFinderV3r2");

    if(! _initAlgImgProc<T,NDim>(data, mask)) { ndarray<const float, 2> empty; return empty; }

    unsigned npeaks = 0;

    for (std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

        size_t ind = (*it)->segind() * m_stride;              
	const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
	const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

        std::vector<Peak>& peaks = (*it) -> peakFinderV3r2<T>(seg_data, seg_mask, rank, r0, dr, nsigm);
	npeaks += peaks.size();
    }

    if(m_pbits & 256) MsgLog(_name(), info, "total number of peaks=" << npeaks);    

    return _ndarrayOfPeakPars(npeaks);
  }

//--------------------
//--------------------
/// mapsOfPhotonNumbersV1 - Chuck's photon counting algorithm - apply fancy correction for split photons.

template <typename T>
ndarray<const nphoton_t, 3>
mapsOfPhotonNumbersV1( const ndarray<const T, 3> data
                     , const ndarray<const mask_t, 3> mask
                     )
{
  if(m_pbits & 256) MsgLog(_name(), info, "in mapsOfPhotonNumbers");

  if(! _initAlgImgProc<T,3>(data, mask)) { ndarray<const nphoton_t, 3> empty; return empty; }

  unsigned shape[3] = {m_nsegs, m_nrows, m_ncols};
  ndarray<nphoton_t, 3> maps(shape);

  for(std::vector<AlgImgProc*>::iterator it = v_algip.begin(); it != v_algip.end(); ++it) {

    size_t ind = (*it)->segind() * m_stride;              
    const ndarray<const T,2>      seg_data(&data.data()[ind], m_sshape);
    const ndarray<const mask_t,2> seg_mask(&m_mask[ind], m_sshape);

    ndarray<nphoton_t, 2>& map = (*it) -> mapOfPhotonNumbersV1<T>(seg_data, seg_mask);
    const Window& win = (*it) -> window();

    for(unsigned r = win.rowmin; r<win.rowmax; r++) 
      for(unsigned c = win.colmin; c<win.colmax; c++)
        maps(win.segind,r,c) = map(r,c);
  }
  return maps;
}

//--------------------
//--------------------

};

} // namespace ImgAlgos

//#include "../src/AlgArrProc.cpp"

#endif // IMGALGOS_ALGARRPROC_H
