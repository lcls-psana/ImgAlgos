#ifndef IMGALGOS_ALGIMGPROC_H
#define IMGALGOS_ALGIMGPROC_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description: see documentation below
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------

#include <string>
#include <vector>
#include <iostream> // for cout, ostream
#include <cstddef>  // for size_t
#include <cstring>  // for memcpy
#include <cmath>    // for sqrt

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "PSCalib/CalibPars.h"  // for pixel_mask_t
#include "MsgLogger/MsgLogger.h"
#include "ndarray/ndarray.h"
#include "ImgAlgos/GlobalMethods.h" // TwoIndexes
#include "ImgAlgos/Window.h"


using namespace std;


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

namespace ImgAlgos {


/**
 * @brief Peak-work parameters
 */ 

struct PeakWork{
  unsigned  peak_npix;
  uint32_t  peak_ireg;
  unsigned  peak_row;
  unsigned  peak_col;
  double    peak_amax;
  double    peak_atot;
  double    peak_ar1;
  double    peak_ar2;
  double    peak_ac1;
  double    peak_ac2;        
  unsigned  peak_rmin;
  unsigned  peak_rmax;
  unsigned  peak_cmin;
  unsigned  peak_cmax;   

  PeakWork( const unsigned& npix=0
	  , const uint32_t& ireg=0
	  , const unsigned& row=0
	  , const unsigned& col=0
	  , const double&   amax=0
	  , const double&   atot=0
	  , const double&   ar1=0
	  , const double&   ar2=0
	  , const double&   ac1=0
	  , const double&   ac2=0        
	  , const unsigned& rmin=0
	  , const unsigned& rmax=0
	  , const unsigned& cmin=0
	  , const unsigned& cmax=0) :
      peak_npix(npix)
    , peak_ireg(ireg)
    , peak_row (row)
    , peak_col (col)
    , peak_amax(amax)
    , peak_atot(atot)
    , peak_ar1 (ar1)
    , peak_ar2 (ar2)
    , peak_ac1 (ac1)
    , peak_ac2 (ac2)        
    , peak_rmin(rmin)
    , peak_rmax(rmax)
    , peak_cmin(cmin)
    , peak_cmax(cmax)
    {}
};

/**
 * @brief Peak parameters
 */ 

struct Peak{
  float seg;
  float row;
  float col;
  float npix;
  float npos;
  float amp_max;
  float amp_tot;
  float row_cgrav; 
  float col_cgrav;
  float row_sigma;
  float col_sigma;
  float row_min;
  float row_max;
  float col_min;
  float col_max;
  float bkgd;
  float noise;
  float son;

  Peak& operator=(const Peak& rhs) {
    seg         = rhs.seg;
    row         = rhs.row;
    col         = rhs.col;
    npix        = rhs.npix;
    npos        = rhs.npos;
    amp_max	= rhs.amp_max;
    amp_tot	= rhs.amp_tot;
    row_cgrav 	= rhs.row_cgrav;
    col_cgrav	= rhs.col_cgrav;
    row_sigma	= rhs.row_sigma;
    col_sigma	= rhs.col_sigma;
    row_min	= rhs.row_min;
    row_max	= rhs.row_max;
    col_min	= rhs.col_min;
    col_max	= rhs.col_max;
    bkgd	= rhs.bkgd;
    noise	= rhs.noise;
    son         = rhs.son;
    return *this;
  }
};

/// Stream insertion operator,
std::ostream& 
operator<<( std::ostream& os, const Peak& p);

/*
struct TwoIndexes {
  int i;
  int j;
};
*/

/**
 * @brief Structure to hold SoN (S/N) algorithm results
 */ 

struct SoNResult {
  double avg; // average intensity in the ring
  double rms; // rms in the ring
  double sig; // raw-avg
  double son; // sig/rms

  SoNResult(const double& av=0, 
            const double& rm=0, 
            const double& sg=0, 
            const double& sn=0) :
    avg(av), rms(rm), sig(sg), son(sn) {}

  SoNResult& operator=(const SoNResult& rhs) {
    avg = rhs.avg;
    rms = rhs.rms;
    sig = rhs.sig;
    son = rhs.son;
    return *this;
  }
};

/**
 * @brief Structure to hold background algorithm results
 */ 

struct BkgdAvgRms {
  double   avg; // average intensity in the ring
  double   rms; // rms in the ring
  unsigned npx; // number of pixels used

  BkgdAvgRms(const double& av=0, 
             const double& rm=0,
             const unsigned& np=0) :
    avg(av), rms(rm), npx(np) {}

  BkgdAvgRms& operator=(const BkgdAvgRms& rhs) {
    avg = rhs.avg;
    rms = rhs.rms;
    npx = rhs.npx;
    return *this;
  }
};

std::ostream& 
operator<<( std::ostream& os, const BkgdAvgRms& b);

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
 *  #include "ImgAlgos/AlgImgProc.h"
 *  #include "ndarray/ndarray.h"     // need it for I/O arrays
 *  
 *  typedef ImgAlgos::AlgImgProc::conmap_t conmap_t;
 *  typedef ImgAlgos::AlgImgProc::pixel_status_t pixel_status_t;
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
 *  ndarray<const wind_t,2> winds = ((0,  0, 185,   0, 388), \
 *                                   (1, 10, 103,  10, 204), \
 *                                   (1, 10, 103, 250, 380));
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
 *   std::vector<Peak>& peaks = aip -> peakFinderV2<T>(seg_data, seg_mask, thr, r0, dr);
 *   std::vector<Peak>& peaks = aip -> peakFinderV3<T>(seg_data, seg_mask, rank, r0, dr);
 *   std::vector<Peak>& peaks = aip -> peakFinderV4<T>(seg_data, seg_mask, thr_low, thr_high, rank, r0, dr);
 *
 *   // The same peak-finders after revisions
 *   std::vector<Peak>& peaks = aip -> peakFinderV2r1<T>(seg_data, seg_mask, thr, r0, dr);
 *   std::vector<Peak>& peaks = aip -> peakFinderV3r1<T>(seg_data, seg_mask, rank, r0, dr, nsigm);
 *   std::vector<Peak>& peaks = aip -> peakFinderV3r2<T>(seg_data, seg_mask, rank, r0, dr, nsigm);
 *   std::vector<Peak>& peaks = aip -> peakFinderV4r1<T>(seg_data, seg_mask, thr_low, thr_high, rank, r0, dr);
 *   std::vector<Peak>& peaks = aip -> peakFinderV4r2<T>(seg_data, seg_mask, thr_low, thr_high, rank, r0, dr);
 *
 *   std::vector<Peak>& peaks = aip -> getVectorOfSelectedPeaks();
 *   std::vector<Peak>& peaks = aip -> getVectorOfPeaks();
 *
 *   // call after peakFinderV2 ONLY!
 *   ndarray<conmap_t, 2>& conmap = aip -> mapOfConnectedPixels();
 *
 *   // call after peakFinderV2, V4r2 ONLY!
 *   ndarray<pixel_status_t, 2>& pixel_status = aip -> mapOfPixelStatus();
 *
 *   // call after peakFinderV3 ONLY!
 *   ndarray<pixel_maximums_t, 2>& locmaxmap = aip -> mapOfLocalMaximums();
 *  @endcode
 *
 *
 *  @li Print methods
 *  @code
 *  aip->printInputPars();
 *  aip->printMatrixOfRingIndexes();
 *  aip->printVectorOfRingIndexes();
 *
 *  Peak& peak = ...
 *  cout << peak ...
 *  @endcode
 */

class AlgImgProc {
public:

  typedef unsigned shape_t;
  typedef PSCalib::CalibPars::pixel_mask_t mask_t;
  typedef uint8_t  u8mask_t;
  typedef uint32_t conmap_t;
  typedef uint16_t pixel_status_t;
  typedef uint16_t pixel_maximums_t;
  typedef uint16_t pixel_minimums_t;
  typedef float    son_t;
  typedef uint16_t nphoton_t;
  typedef float    fphoton_t;

  /**
   * @brief Class constructor is used for initialization of all paramaters. 
   * 
   * @param[in] seg    - ROI segment index in the ndarray
   * @param[in] rowmin - ROI window limit
   * @param[in] rowmax - ROI window limit
   * @param[in] colmin - ROI window limit
   * @param[in] colmax - ROI window limit
   * @param[in] pbits  - print control bit-word; =0-print nothing, +1-input parameters, +2-algorithm area, +128-all details.
   * @param[in] npksmax - maximal number of reserved peak records.
   */

  AlgImgProc( const size_t&   seg
	    , const size_t&   rowmin  = 0
	    , const size_t&   rowmax  = 1000000
	    , const size_t&   colmin  = 0
	    , const size_t&   colmax  = 1000000
	    , const unsigned& pbits   = 0
	    , const unsigned& npksmax = 80000
	    ) ;

  /// Destructor
  virtual ~AlgImgProc();

  /// Prints memeber data
  void printInputPars();

  /// Set segment index in the >2-d ndarray
  //void setSegment(const size_t& seg){ m_seg = seg; }

  /**
   * @brief Set median (S/N) algorithm parameters
   * @param[in] r0 - radial parameter of the ring for S/N evaluation algorithm
   * @param[in] dr - ring width for S/N evaluation algorithm 
   */
  void setSoNPars(const float& r0=5, const float& dr=0.05);

  /// Set peak selection parameters
  void setPeakSelectionPars(const float& npix_min=0, 
                            const float& npix_max=1e6, 
                            const float& amax_thr=0, 
                            const float& atot_thr=0,
                            const float& son_min=0);

  /// Returns reference to Window object
  const Window& window(){ return m_win; }

  /// Returns reference to Window object
  void validate_window(const Window::shape_t* shape=0){ if (shape) m_win.validate(shape); }

  /// Returns segment index in the ndarray
  const size_t& segind(){ return m_seg; }

  /// Returns vector of all found peaks for this segment/window
  std::vector<Peak>& getVectorOfPeaks(){ return v_peaks; }

  /// Returns vector of selected peaks for this segment/window
  std::vector<Peak>& getVectorOfSelectedPeaks(){ return v_peaks_sel; }

  /// Returns map of connected pixels after peakFinderV2(.)
  ndarray<conmap_t, 2>& mapOfConnectedPixels() { return m_conmap; }

  /// Returns map of pixel status after peakFinderV2, V4r2
  ndarray<pixel_status_t, 2>& mapOfPixelStatus() { return m_pixel_status; }

  /// Returns map of local maximums after peakFinderV3(.)
  ndarray<pixel_minimums_t, 2>& mapOfLocalMinimums() { return m_local_minimums; }

  /// Returns map of local maximums after peakFinderV3(.)
  ndarray<pixel_maximums_t, 2>& mapOfLocalMaximums() { return m_local_maximums; }

  /// Prints indexes for S/N algorithm
  void printMatrixOfRingIndexes();
  void printVectorOfRingIndexes();

  /// Prints indexes for peakFinderV3 algorithm
  void printMatrixOfDiagIndexes();
  void printVectorOfDiagIndexes();

  // Copy constructor and assignment are disabled by default
  AlgImgProc ( const AlgImgProc& ) ;
  AlgImgProc& operator = ( const AlgImgProc& ) ;

private:

  unsigned m_pbits;    // pirnt control bit-word
  unsigned m_npksmax;  // maximal number of peaks in segment/window
  size_t   m_seg;      // segment index (for ndarray with ndim>2)

  Window   m_win;      // work area window

  bool     m_init_son_is_done; // for S/N algorithm
  bool     m_init_bkgd_is_done; // for _evaluateBkgdAvgRms algorithm
  bool     m_use_mask;

  conmap_t m_numreg;

  float    m_r0;       // radial parameter of the ring for S/N evaluation algorithm
  float    m_dr;       // ring width for S/N evaluation algorithm 
  size_t   m_rank;     // rank of maximum for peakFinderV3

  SoNResult  m_sonres_def;
  BkgdAvgRms m_bkgdavgrms_def;
  BkgdAvgRms m_bkgd;

  float    m_peak_npix_min; // peak selection parameter
  float    m_peak_npix_max; // peak selection parameter
  float    m_peak_amax_thr; // peak selection parameter
  float    m_peak_atot_thr; // peak selection parameter
  float    m_peak_son_min;  // peak selection parameter
  bool     m_do_preselect;  // flag for applying peak pre-selection before S/N calculation


  // parameters for _isDroplet
  int      m_reg_rmin; // region limit
  int      m_reg_rmax; // region limit
  int      m_reg_cmin; // region limit
  int      m_reg_cmax; // region limit
  double   m_reg_a0;   // intensity in the initial point for droplet
  double   m_reg_thr;  // threshold on intensity for V3r2

  std::vector<TwoIndexes> v_ind_pixgrp; // vector of pixel indexes for droplet
  std::vector< std::vector<TwoIndexes> > vv_peak_pixinds; // vector of peak vector of pixel indexes
  unsigned m_pixgrp_max_size; // size of droplet vector

  //Peak     m_peak;

  ndarray<pixel_status_t, 2>   m_pixel_status;
  ndarray<pixel_maximums_t, 2> m_local_maximums;
  ndarray<pixel_minimums_t, 2> m_local_minimums;
  ndarray<conmap_t, 2>         m_conmap;
  std::vector<PeakWork>        v_peaks_work;
  std::vector<Peak>            v_peaks;
  std::vector<Peak>            v_peaks_sel;
  std::vector<TwoIndexes>      v_indexes;
  std::vector<TwoIndexes>      v_inddiag;
  std::vector<TwoIndexes>      v_indcross;
  ndarray<nphoton_t, 2>        m_nphoton; // map with integer (floor) number of photons 
  ndarray<fphoton_t, 2>        m_fphoton; // map with fractional number of photons  
  ndarray<nphoton_t, 2>        m_mphoton; // map of merged (integer) number of photons

  //ndarray<Peak, 1>           m_peaks;
  // mask_t*                    m_seg_mask_def;


  /// Returns string name of the class for messanger
  inline const char* _name() {return "ImgAlgos::AlgImgProc";}

  /// Reserves vectors of peaks sizes
  void _reserveVectorOfPeaks();

  /// Recursive method which checks m_pixel_status[r][c] and numerates connected regions in m_conmap[r][c]
  void _findConnectedPixels(const int& r, const int& c);

  /// Recursive method finding connected pixels in constrained region and filling vector of indexes.
  void _findConnectedPixelsInRegion(const int& r, const int& c);

  /// Adds droplet (vector of connected pixels v_ind_pixgrp) to the map m_conmap
  void _addPixGroupToMap();

  /// Clears status of used pixels in m_pixel_status[r][c] if central pixel is not maximal
  void _clearStatusOfUnusedPixels();

  /// Recursive method finding connected pixels in constrained region and filling vector of indexes for V4r2
  template <typename T>
  bool _findConnectedPixelsInRegionV1(const ndarray<const T,2>& data, const int& r, const int& c);

  /// Recursive method finding connected pixels in constrained region and filling vector of indexes for V3r2.
  template <typename T>
  void _findConnectedPixelsInRegionV2(const ndarray<const T,2>& data, const ndarray<const mask_t,2>& mask, const int& r, const int& c);

  /// Makes m_conmap - map of connected pixels with enumerated regions from m_pixel_status and counts m_numreg
  void _makeMapOfConnectedPixels();

  /// Decide whether PeakWork should be processed and included in the v_peaks
  bool _peakWorkIsPreSelected(const PeakWork& pw);

  /// Decide if peak should be processed or not and included in the v_peaks 
  bool _peakIsPreSelected(const Peak& peak);

  /// Decide if peak should be included or not in the output v_peaks
  bool _peakIsSelected(const Peak& peak);

  /// Makes vector of peaks v_peaks from v_peaks_work
  void _makeVectorOfPeaks();

  /// Makes vector of selected peaks v_peaks_sel from v_peaks
  void _makeVectorOfSelectedPeaks();

  /// Evaluate ring indexes for S/N algorithm
  void _evaluateRingIndexes(const float& r0, const float& dr);

  /// Evaluate "diagonal" region indexes for peakFinderV3
  void _evaluateDiagIndexes(const size_t& rank);

  /// Fill "cross" region indexes for mapOfPhotonNumbersV1 
  void _fillCrossIndexes();

  /// prints statistics of local maximums and minimums
  void _printStatisticsOfLocalExtremes();

  /// Join local maximum fractional intensity with largest adjesent pixel (very special algorithm for Chuck's photon counting)
  void _mergeConnectedPixelCouples(const fphoton_t& thr_on_max=0.5, const fphoton_t& thr_on_tot=0.9, const bool& DO_TEST = false);

  /// algorithmical part running inside _mergeConnectedPixelCouples
  void _mergePixelCouplesInArea(
    const fphoton_t thr_on_max,
    const fphoton_t thr_on_tot,
    const bool      DO_TEST,
    const unsigned  rowmin,
    const unsigned  rowmax1,
    const unsigned  colmin,
    const unsigned  colmax1,
    const vector<TwoIndexes>& vneighbs
  );

//--------------------
  /**
   * @brief Makes m_pixel_status array by setting to 1/0 all good-above-threshold/bad pixels
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr  - threshold on data value
   */

template <typename T>
void
_makeMapOfPixelStatus( const ndarray<const T,2>& data
                     , const ndarray<const mask_t,2>& mask
                     , const T& thr
                     )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfPixelStatus, seg=" << m_seg
                           << " thr=" << thr << "\n    in window: " << m_win);
  //if(m_pbits & 512) m_win.print();

  //if(m_pixel_status.size()==0) 
  if(m_pixel_status.empty()){
     m_pixel_status = make_ndarray<pixel_status_t>(data.shape()[0], data.shape()[1]);
  }

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      m_pixel_status(r, c) = (mask(r, c) && (data(r, c)>thr)) ? 255 : 0;
}

//--------------------
  /**
   * @brief Makes 2-threshold map of pixel status m_pixel_status by setting
   *         0 - masked pixel
   *        +1 - intensity above thr_high
   *        +2 - intensity above thr_low
   *        +4 - intensity below thresholds
   *        ------------------ 
   *        m_pixel_status array is used later to set:
   *        +8  - used pixel in recursion
   *        +16 - used pixel on map
   *        +32 - local maximum
   * 
   * @param[in] data - ndarray with calibrated intensities
   * @param[in] mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in] thr_low   - low threshold on data value
   * @param[in] thr_high  - high threshold on data value
   */

template <typename T>
void
_makeMapOfPixelStatusV2( const ndarray<const T,2>& data
                       , const ndarray<const mask_t,2>& mask
                       , const T& thr_low
                       , const T& thr_high
                       )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfPixelStatusV2, seg=" << m_seg
                           << " thr_high=" << thr_high << " thr_low=" << thr_low 
                           << "\n    in window: " << m_win);
  //if(m_pbits & 512) m_win.print();

  //if(m_pixel_status.size()==0) 
  if(m_pixel_status.empty()) {
    m_pixel_status = make_ndarray<pixel_status_t>(data.shape()[0], data.shape()[1]);

    // they are not used in this algo, but in case if someone request them...
    _init_local_extreams_arrays<T>(data);
  }

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      if(! mask(r, c))             m_pixel_status(r, c) = 0; // pixel is masked
      else if(data(r, c)<thr_low)  m_pixel_status(r, c) = 4; // a<thr_low
      else if(data(r, c)<thr_high) m_pixel_status(r, c) = 2; // a>=thr_low, but a<thr_high
      else                         m_pixel_status(r, c) = 3; // a>=thr_high,low: set 1st and 2nd bits
    }
}

//--------------------
  /**
   * @brief Process data ndarray using map of connected pixels m_conmap 
   * and collect peak information in std::vector<PeakWork> v_peaks_work 
   * 
   * @param[in]  data - ndarray with calibrated intensities
   */

template <typename T>
void
_procConnectedPixels(const ndarray<const T,2>& data)
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _procConnectedPixels, seg=" << m_seg);
  //if(m_pbits & 512) m_win.print();

  if(m_numreg==0) return;
  //if(m_numreg>8190) return;

  PeakWork pw0; // def init with 0
  pw0.peak_cmin = m_win.colmax;
  pw0.peak_rmin = m_win.rowmax;

  if(v_peaks_work.capacity() != m_npksmax) v_peaks_work.reserve(m_npksmax);
  v_peaks_work.clear();

  //std::fill_n(&v_peaks_work[0], int(m_numreg+1), pw0);
  for(unsigned reg = 0; reg < min(m_numreg+1, m_npksmax); reg++)
    v_peaks_work.push_back(pw0);

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      conmap_t ireg = m_conmap(r, c);

      if(ireg<1) continue;

      if(!(ireg<m_npksmax)) {
         MsgLog(_name(), warning, "Number of peak regions: " << ireg 
		<< " reached the max reserved number of peaks: " << m_npksmax
                << " in segment: " << m_seg);
         break;
      }

      //std::cout << " reg=" << ireg;

      double amp = (double)data(r, c);

      PeakWork& pw = v_peaks_work[ireg];
      pw.peak_npix ++;   
      pw.peak_ireg = ireg;
      pw.peak_atot += amp;   
      pw.peak_ar1  += amp*r;    
      pw.peak_ar2  += amp*r*r;    
      pw.peak_ac1  += amp*c;    
      pw.peak_ac2  += amp*c*c;    

      if(amp > pw.peak_amax) {
        pw.peak_row  = r;
	pw.peak_col  = c;
        pw.peak_amax = amp;   
      }

      if(c < pw.peak_cmin) pw.peak_cmin = c;
      if(c > pw.peak_cmax) pw.peak_cmax = c;
      if(r < pw.peak_rmin) pw.peak_rmin = r;
      if(r > pw.peak_rmax) pw.peak_rmax = r;
    }
  }
}

//--------------------
public:
//--------------------
  /**
   * @brief Makes map of local minimums of requested rank (IN ROWS AND COLUMNS ONLY!)
   * 
   * Map of local minimums is a 2-d array of unsigned integer values of data shape, 
   * it has 0/1/2 bits for non-maximum / maximum in column / maximum in row of radius rank.   
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */

template <typename T>
void
_makeMapOfLocalMinimums( const ndarray<const T,2>& data
                       , const ndarray<const mask_t,2>& mask
                       , const size_t& rank
                       )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfLocalMinimums, seg=" << m_seg << " rank=" << rank << "\n    in window: " << m_win);
  if(m_pbits & 512) m_win.print();

  if(m_local_minimums.empty()) 
     m_local_minimums = make_ndarray<pixel_minimums_t>(data.shape()[0], data.shape()[1]);
  std::fill_n(&m_local_minimums(0,0), int(data.size()), pixel_minimums_t(0));

  unsigned rmin = (int)m_win.rowmin;
  unsigned rmax = (int)m_win.rowmax;				  
  unsigned cmin = (int)m_win.colmin;
  unsigned cmax = (int)m_win.colmax;
  int irank = (int)rank;

  // check rank minimum in columns and set the 1st bit (1)
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      if(!mask(r, c)) continue;
      m_local_minimums(r, c) = 1;

      // positive side of c 
      unsigned dmax = min((int)cmax-1, int(c)+irank);
      for(unsigned cd=c+1; cd<=dmax; cd++) {
        if(mask(r, cd) && (data(r, cd) < data(r, c))) {
          m_local_minimums(r, c) &=~1; // clear 1st bit
          c=cd-1; // jump ahead 
	  break;
	}
      }

      if(m_local_minimums(r, c) & 1) {
        // negative side of c 
        unsigned dmin = max((int)cmin, int(c)-irank);
        for(unsigned cd=dmin; cd<c; cd++) {
          if(mask(r, cd) && (data(r, cd) < data(r, c))) {
            m_local_minimums(r, c) &=~1; // clear 1st bit
            c=cd+rank; // jump ahead 
	    break;
	  }
        }
      }

      // (r,c) is a local dip, jump ahead through the tested rank range
      if(!(c<cmax)) break;
      if(m_local_minimums(r, c) & 1) c+=rank;
    }
  }

  // check rank minimum in rows and set the 2nd bit (2)
  for(unsigned c = cmin; c<cmax; c++) {
    for(unsigned r = rmin; r<rmax; r++) {
      // if it is not a local maximum from previous algorithm
      //if(!m_local_minimums[r][c]) continue;

      if(!mask(r, c)) continue;
      m_local_minimums(r, c) |= 2; // set 2nd bit

      // positive side of r 
      unsigned dmax = min((int)rmax-1, int(r)+irank);
      for(unsigned rd=r+1; rd<=dmax; rd++) {
        if(mask(rd, c) && (data(rd, c) < data(r, c))) {
          m_local_minimums(r, c) &=~2; // clear 2nd bit
          r=rd-1; // jump ahead 
	  break;
	}
      }

      if(m_local_minimums(r, c) & 2) {
        // negative side of r
        unsigned dmin = max((int)rmin, int(r)-irank);
        for(unsigned rd=dmin; rd<r; rd++) {
          if(mask(rd,c) && (data(rd,c) < data(r, c))) {
            m_local_minimums(r, c) &=~2; // clear 2nd bit
            r=rd+rank; // jump ahead 
	    break;
	  }
        }
      }

      // (r,c) is a local dip, jump ahead through the tested rank range
      if(!(r<rmax)) break;
      if(m_local_minimums(r, c) & 2) r+=rank;
    }
  }
}

//--------------------
  /**
   * @brief Makes map of local maximums of requested rank
   * 
   * Map of local maximums is a 2-d array of unsigned integer values of data shape, 
   * it has 0/1/2/4 bits for non-maximum / maximum in column / maximum in row / local maximum in rectangle of radius rank.   
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */

template <typename T>
void
_makeMapOfLocalMaximums( const ndarray<const T,2>& data
                       , const ndarray<const mask_t,2>& mask
                       , const size_t& rank
                       )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfLocalMaximums, seg=" << m_seg << " rank=" << rank << "\n    in window: " << m_win);
  if(m_pbits & 512) m_win.print();

  // initialization of indexes
  if(v_inddiag.empty()) _evaluateDiagIndexes(rank);

  if(m_local_maximums.empty()) 
     m_local_maximums = make_ndarray<pixel_maximums_t>(data.shape()[0], data.shape()[1]);
  std::fill_n(&m_local_maximums(0,0), int(data.size()), pixel_maximums_t(0));

  unsigned rmin = (int)m_win.rowmin;
  unsigned rmax = (int)m_win.rowmax;				  
  unsigned cmin = (int)m_win.colmin;
  unsigned cmax = (int)m_win.colmax;
  int irank = (int)rank;

  // check rank maximum in columns and set the 1st bit (1)
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      if(!mask(r, c)) continue;
      m_local_maximums(r, c) = 1;

      // positive side of c 
      unsigned dmax = min((int)cmax-1, int(c)+irank);
      for(unsigned cd=c+1; cd<=dmax; cd++) {
        if(mask(r, cd) && (data(r, cd) > data(r, c))) {
          m_local_maximums(r, c) &=~1; // clear 1st bit
          c=cd-1; // jump ahead (+1 incremented in c-loop)
	  break;
	}
      }

      if(m_local_maximums(r, c) & 1) {
        // negative side of c 
        unsigned dmin = max((int)cmin, int(c)-irank);
        for(unsigned cd=dmin; cd<c; cd++) {
          if(mask(r, cd) && (data(r, cd) > data(r, c))) {
            m_local_maximums(r, c) &=~1; // clear 1st bit
            c=cd+rank; // jump ahead 
	    break;
	  }
        }
      }

      // (r,c) is a local dip, jump ahead through the tested rank range
      if(!(c<cmax)) break;
      if(m_local_maximums(r, c) & 1) c+=rank;
    }
  }

  // check rank maximum in rows and set the 2nd bit (2)
  for(unsigned c = cmin; c<cmax; c++) {
    for(unsigned r = rmin; r<rmax; r++) {
      // if it is not a local maximum from previous algorithm
      //if(!m_local_maximums[r][c]) continue;

      if(!mask(r, c)) continue;
      m_local_maximums(r, c) |= 2; // set 2nd bit

      // positive side of r 
      unsigned dmax = min((int)rmax-1, int(r)+irank);
      for(unsigned rd=r+1; rd<=dmax; rd++) {
        if(mask(rd, c) && (data(rd, c) > data(r, c))) {
          m_local_maximums(r, c) &=~2; // clear 2nd bit
          r=rd-1; // jump ahead 
	  break;
	}
      }

      if(m_local_maximums(r, c) & 2) {
        // negative side of r
        unsigned dmin = max((int)rmin, int(r)-irank);
        for(unsigned rd=dmin; rd<r; rd++) {
          if(mask(rd, c) && (data(rd, c) > data(r, c))) {
            m_local_maximums(r, c) &=~2; // clear 2nd bit
            r=rd+rank; // jump ahead 
	    break;
	  }
        }
      }

      // (r,c) is a local dip, jump ahead through the tested rank range
      if(!(r<rmax)) break;
      if(m_local_maximums(r, c) & 2) r+=rank;
    }
  }

  // check rank maximum in "diagonal" regions and set the 3rd bit (4)
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      // if it is not a local maximum from two previous algorithm

      if(m_local_maximums(r, c) != 3) continue;
      m_local_maximums(r, c) |= 4; // set 3rd bit

      for(vector<TwoIndexes>::const_iterator ij  = v_inddiag.begin();
                                             ij != v_inddiag.end(); ij++) {
        int ir = r + (ij->i);
        int ic = c + (ij->j);

        if(  ir<(int)rmin)  continue;
        if(  ic<(int)cmin)  continue;
        if(!(ir<(int)rmax)) continue;
        if(!(ic<(int)cmax)) continue;

        if(mask(ir, ic) && (data(ir, ic) > data(r, c))) {
          m_local_maximums(r, c) &=~4; // clear 3rd bit
	  break;
	}
      }

      // (r,c) is a local peak, jump ahead through the tested rank range
      if(m_local_maximums(r, c) & 4) c+=rank;
    }
  }

}

//--------------------
//--------------------
//--------------------
//--------------------
  /**
   * @brief Makes map of local minimums of requested rank (IN ROWS AND COLUMNS ONLY!)
   * 
   * V0 - pixels at distance < rank from the boarder are not considered for local extremes.
   * Map of local minimums is a 2-d array of unsigned integer values of data shape, 
   * it has 0/1/2 bits for non-maximum / maximum in column / maximum in row of radius rank.   
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */

template <typename T>
void
_makeMapOfLocalMinimumsV0( const ndarray<const T,2>& data
                         , const ndarray<const mask_t,2>& mask
                         , const size_t& rank
                         )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfLocalMinimumsV0, seg=" << m_seg << " rank=" << rank << "\n    in window: " << m_win);
  //if(m_pbits & 512) m_win.print();

  if(m_local_minimums.empty()) 
     m_local_minimums = make_ndarray<pixel_minimums_t>(data.shape()[0], data.shape()[1]);
  std::fill_n(&m_local_minimums(0,0), int(data.size()), pixel_minimums_t(0));

  unsigned rmin = max((int)m_win.rowmin, int(0+rank));
  unsigned rmax = min((int)m_win.rowmax, int(data.shape()[0]-rank));
  unsigned cmin = max((int)m_win.colmin, int(0+rank));
  unsigned cmax = min((int)m_win.colmax, int(data.shape()[1]-rank));


  // check rank minimum in columns and set the 1st bit (1)
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      if(!mask(r, c)) continue;
      m_local_minimums(r, c) = 1;
      int cp=c; 
      int cm=c;
      for(unsigned d=0; d<rank; ++d) {
        cp++; cm--;
        if((mask(r, cp) && (data(r, cp) < data(r, c)))
        || (mask(r, cm) && (data(r, cm) < data(r, c)))) {
          m_local_minimums(r, c) &=~1; // clear 1st bit
          c=--cp; // jump ahead 
	  break;
	}
      }
      // (r,c) is a local dip, jump ahead through the tested rank range
      if(m_local_minimums(r, c) & 1) c+=rank;
    }
  }

  // check rank minimum in rows and set the 2nd bit (2)
  for(unsigned c = cmin; c<cmax; c++) {
    for(unsigned r = rmin; r<rmax; r++) {
      // if it is not a local maximum from previous algorithm
      //if(!m_local_minimums[r][c]) continue;
      if(!mask(r, c)) continue;
      m_local_minimums(r, c) |= 2; // set 2nd bit
      int rp=r; 
      int rm=r;
      for(unsigned d=0; d<rank; ++d) { 
        rp++; rm--;
        if((mask(rp, c) && (data(rp, c) < data(r, c)))
        || (mask(rm, c) && (data(rm, c) < data(r, c)))) {
          m_local_minimums(r, c) &=~2; // clear 2nd bit
          r=--rp; // jump ahead 
	  break;
	}
      }
      // if (r,c) is a local peak, jump ahead through the tested rank range
      if(m_local_minimums(r, c) & 2) r+=rank;
    }
  }
}

//--------------------
  /**
   * @brief Makes map of local maximums of requested rank
   * 
   * V0 - pixels at distance < rank from the boarder are not considered for local extremes.
   * Map of local maximums is a 2-d array of unsigned integer values of data shape, 
   * it has 0/1/2/4 bits for non-maximum / maximum in column / maximum in row / local maximum in rectangle of radius rank.   
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */

template <typename T>
void
_makeMapOfLocalMaximumsV0( const ndarray<const T,2>& data
                         , const ndarray<const mask_t,2>& mask
                         , const size_t& rank
                         )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfLocalMaximumsV0, seg=" << m_seg << " rank=" << rank << "\n    in window: " << m_win);
  //if(m_pbits & 512) m_win.print();

  // initialization of indexes
  if(v_inddiag.empty()) _evaluateDiagIndexes(rank);

  if(m_local_maximums.empty()) 
     m_local_maximums = make_ndarray<pixel_maximums_t>(data.shape()[0], data.shape()[1]);
  std::fill_n(&m_local_maximums(0,0), int(data.size()), pixel_maximums_t(0));

  unsigned rmin = max((int)m_win.rowmin, int(0+rank));
  unsigned rmax = min((int)m_win.rowmax, int(data.shape()[0]-rank));
  unsigned cmin = max((int)m_win.colmin, int(0+rank));
  unsigned cmax = min((int)m_win.colmax, int(data.shape()[1]-rank));

  // check rank maximum in columns and set the 1st bit (1)
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      if(!mask(r, c)) continue;
      m_local_maximums(r, c) = 1;
      int cp=c; 
      int cm=c;
      for(unsigned d=0; d<rank; ++d) {
        cp++; cm--;
        if((mask(r, cp) && (data(r,cp) > data(r, c)))
	|| (mask(r, cm) && (data(r,cm) > data(r, c)))) {
          m_local_maximums(r, c) &=~1; // clear 1st bit
          c=--cp; // jump ahead 
	  break;
	}
      }
      // (r,c) is a local peak, jump ahead through the tested rank range
      if(m_local_maximums(r, c) & 1) c+=rank;
    }
  }

  // check rank maximum in rows and set the 2nd bit (2)
  for(unsigned c = cmin; c<cmax; c++) {
    for(unsigned r = rmin; r<rmax; r++) {
      // if it is not a local maximum from previous algorithm
      //if(!m_local_maximums[r][c]) continue;
      if(!mask(r, c)) continue;
      m_local_maximums(r, c) |= 2; // set 2nd bit
      int rp=r; 
      int rm=r;
      for(unsigned d=0; d<rank; ++d) { 
        rp++; rm--;
        if((mask(rp, c) && (data(rp, c) > data(r, c)))
        || (mask(rm, c) && (data(rm, c) > data(r, c)))) {
          m_local_maximums(r, c) &=~2; // clear 2nd bit
          r=--rp; // jump ahead 
	  break;
	}
      }
      // if (r,c) is a local peak, jump ahead through the tested rank range
      if(m_local_maximums(r, c) & 2) r+=rank;
    }
  }

  // check rank maximum in "diagonal" regions and set the 3rd bit (4)
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      // if it is not a local maximum from two previous algorithm
      if(m_local_maximums(r, c) != 3) continue;
      m_local_maximums(r, c) |= 4; // set 3rd bit

      for(vector<TwoIndexes>::const_iterator ij  = v_inddiag.begin();
                                             ij != v_inddiag.end(); ij++) {
        int ir = r + (ij->i);
        int ic = c + (ij->j);

        if(mask(ir,ic) && (data(ir,ic) > data(r, c))) {
          m_local_maximums(r, c) &=~4; // clear 3rd bit
	  break;
	}
      }
      // (r,c) is a local peak, jump ahead through the tested rank range
      if(m_local_maximums(r, c) & 4) c+=rank;
    }
  }
}

//--------------------
//--------------------
  /**
   * @brief Makes map of local maximums of runk=1 cross(+) region (very special case for Chuck's algorithm)
   * 
   * It does not use mask, because works with special array of "fractional number of photons" - bad pixels are set to 0.
   * Map of local maximums is a 2-d array of unsigned integer values of data shape, 
   * it has 0/1/2 bits for non-maximum / maximum in column, then local maximum in cross = 3.   
   * @param[in]  data - ndarray with (positive or zero) fractional number of photons
   */

template <typename T>
void
_makeMapOfLocalMaximumsRank1Cross(const ndarray<const T,2>& data)
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfLocalMaximumsRank1Cross, seg=" << m_seg << "\n    in window: " << m_win);
  if(m_pbits & 512) m_win.print();

  if(m_local_maximums.empty()) 
     m_local_maximums = make_ndarray<pixel_maximums_t>(data.shape()[0], data.shape()[1]);
  std::fill_n(&m_local_maximums(0,0), int(data.size()), pixel_maximums_t(0));

  unsigned rmin = (int)m_win.rowmin;
  unsigned rmax = (int)m_win.rowmax;				  
  unsigned cmin = (int)m_win.colmin;
  unsigned cmax = (int)m_win.colmax;

  // check local maximum in columns and set the 1st bit (1)
  for(unsigned r = rmin; r<rmax; r++) {

    // first pixel in the row
    unsigned c = cmin;
    if(data(r, c) > data(r,c+1)) {
      m_local_maximums(r, c) |= 1;  // set 1st bit
      c+=2;
    }
    else c+=1;

    // all internal pixels in the row
    for(; c<cmax-1; c++) {
      if(data(r, c+1) > data(r, c)) continue;         // go to the next pixel
      if(data(r, c-1) > data(r, c)) {c+=1; continue;} // jump ahead
      m_local_maximums(r, c) |= 1;  // set 1st bit
      c+=1; // jump ahead 
    }

    // last pixel in the row
    if(data(r,cmax-1) > data(r,cmax-2)) m_local_maximums(r,cmax-1) |= 1;  // set 1st bit
  } // rows loop

  // check local maximum in rows and set the 2nd bit (2)
  for(unsigned c = cmin; c<cmax; c++) {

    // first pixel in the column
    unsigned r = rmin;
    if(data(r, c) > data(r+1, c)) {
      m_local_maximums(r, c) |= 2; // set 2nd bit
      r+=2;
    }
    else r+=1;

    // all internal pixels in the column
    for(; r<rmax-1; r++) {
      if(data(r+1, c) > data(r, c)) continue;         // go to the next pixel
      if(data(r-1, c) > data(r, c)) {r+=1; continue;} // jump ahead
      m_local_maximums(r, c) |= 2; // set 2nd bit
      r+=1; // jump ahead 
    }

    // last pixel in the column
    if(data(rmax-1,c) > data(rmax-2,c)) m_local_maximums(rmax-1,c) |= 2;  // set 2nd bit
  } // columns loop
}

//--------------------
//--------------------
//--------------------
private:
//--------------------
  /**
   * @brief _procLocalMaximum - process local maximum and fill pre-selected peak in v_peaks.
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   * @param[in]  r0 - droplet central pixel row-coordinate 
   * @param[in]  c0 - droplet central pixel column-coordinate   
   */

template <typename T>
void
_procLocalMaximum( const ndarray<const T,2>& data
                 , const ndarray<const mask_t,2>& mask
                 , const size_t& rank
                 , const unsigned& r0
                 , const unsigned& c0
                 )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _procLocalMaximum, seg=" << m_seg << " r0=" << r0 << " c0=" << c0 << " rank=" << rank);

  double   a0 = data(r0,c0);
  unsigned npix = 0;
  double   samp = 0;
  double   sac1 = 0;
  double   sac2 = 0;
  double   sar1 = 0;
  double   sar2 = 0;

  unsigned rmin = std::max((int)m_win.rowmin, int(r0-rank));
  unsigned rmax = std::min((int)m_win.rowmax, int(r0+rank+1));
  unsigned cmin = std::max((int)m_win.colmin, int(c0-rank));
  unsigned cmax = std::min((int)m_win.colmax, int(c0+rank+1));

  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {

      if(!mask(r, c)) continue;
      double a = data(r, c);
      if(!(a>0)) continue;
      npix += 1;
      samp += a;
      sar1 += a*r;
      sac1 += a*c;
      sar2 += a*r*r;
      sac2 += a*c*c;
    }
  }

  if(npix<1) return;

  Peak peak;

  peak.seg       = m_seg;
  peak.row       = r0;
  peak.col       = c0;
  peak.npix      = npix;
  peak.npos      = npix;
  peak.amp_max   = a0;
  peak.amp_tot   = samp;
  peak.row_cgrav = sar1/samp;
  peak.col_cgrav = sac1/samp;
  peak.row_sigma = (npix>1) ? std::sqrt( sar2/samp - peak.row_cgrav * peak.row_cgrav ) : 0;
  peak.col_sigma = (npix>1) ? std::sqrt( sac2/samp - peak.col_cgrav * peak.col_cgrav ) : 0;
  peak.row_min   = rmin;
  peak.row_max   = rmax;
  peak.col_min   = cmin;
  peak.col_max   = cmax;  
  peak.bkgd      = 0; //sonres.avg;
  peak.noise     = 0; //sonres.rms;
  peak.son       = 0; //sonres.son;

  if(_peakIsPreSelected(peak) && v_peaks.size()<m_npksmax-1) v_peaks.push_back(peak);
}

//--------------------
  /**
   * @brief _procLocalMaximumV2 - the same as _procLocalMaximum but count all intensities in the rank region.
   */

template <typename T>
void
_procLocalMaximumV2( const ndarray<const T,2>& data
                   , const ndarray<const mask_t,2>& mask
                   , const size_t& rank
                   , const unsigned& r0
                   , const unsigned& c0
                   )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _procLocalMaximumV2, seg=" << m_seg << " r0=" << r0 << " c0=" << c0 << " rank=" << rank);

  double   a0 = data(r0,c0);
  unsigned npix = 0;
  unsigned npos = 0;
  double   samp = 0;
  double   swei = 0;
  double   sac1 = 0;
  double   sac2 = 0;
  double   sar1 = 0;
  double   sar2 = 0;

  unsigned rmin = std::max((int)m_win.rowmin, int(r0-rank));
  unsigned rmax = std::min((int)m_win.rowmax, int(r0+rank+1));
  unsigned cmin = std::max((int)m_win.colmin, int(c0-rank));
  unsigned cmax = std::min((int)m_win.colmax, int(c0+rank+1));

  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {

      if(!mask(r, c)) continue;
      double a = data(r, c);
      samp += a;
      npix += 1;
      if(!(a>0)) continue;
      npos += 1;
      swei += a;
      sar1 += a*r;
      sac1 += a*c;
      sar2 += a*r*r;
      sac2 += a*c*c;
    }
  }

  if(!(samp>0)) return;
  if(npos<npix-npos) return;
  //if(npos<1) return;
  //cout << "fraction of positive=" << float(npos)/npix << '\n';

  Peak peak;

  peak.seg       = m_seg;
  peak.row       = r0;
  peak.col       = c0;
  peak.npix      = npix;
  peak.npos      = npos;
  peak.amp_max   = a0;
  peak.amp_tot   = samp;

  if (swei>0) {
    peak.row_cgrav = sar1/swei;
    peak.col_cgrav = sac1/swei;
    peak.row_sigma = std::sqrt(sar2/swei - peak.row_cgrav * peak.row_cgrav);
    peak.col_sigma = std::sqrt(sac2/swei - peak.col_cgrav * peak.col_cgrav);  
  } 
  else {
    peak.row_cgrav = r0;
    peak.col_cgrav = c0;
    peak.row_sigma = 0;
    peak.col_sigma = 0;
  }

  peak.row_min   = rmin;
  peak.row_max   = rmax;
  peak.col_min   = cmin;
  peak.col_max   = cmax;  
  peak.bkgd      = 0; //sonres.avg;
  peak.noise     = 0; //sonres.rms;
  peak.son       = 0; //sonres.son;

  if(v_peaks.size()<m_npksmax-1) v_peaks.push_back(peak);
  //if(_peakIsPreSelected(peak) && v_peaks.size()<m_npksmax-1) v_peaks.push_back(peak);
}

//--------------------
  /**
   * @brief _procLocalMaximumV3 - process peak after S/N is available.
   */

template <typename T>
void
_procLocalMaximumV3( const ndarray<const T,2>& data
                   , const ndarray<const mask_t,2>& mask
                   , const size_t& rank
                   , Peak& peak
		   , const float& nsigm=0 // 0-turns off threshold algorithm, 1.64-leaves 5% of noise, etc.
                   )
{
  unsigned r0   = (unsigned)peak.row;   // already filled in _makePeaksFromMapOfLocalMaximumsV3
  unsigned c0   = (unsigned)peak.col;   // already filled in _makePeaksFromMapOfLocalMaximumsV3
  //double bkgd = peak.bkgd;  // already evaluated in _addBkgdAvgRmsToPeaks;
  //double noise= peak.noise; // already evaluated in _addBkgdAvgRmsToPeaks;

  if(m_pbits & 512) MsgLog(_name(), info, "in _procLocalMaximumV3, seg=" << m_seg 
                           << " r0=" << r0 << " c0=" << c0 << " rank=" << rank);
   
  double   a0  = data(r0,c0) - peak.bkgd;
  double   thr = (nsigm) ? peak.noise * nsigm : 0;
  double   noise_tot = 0;

  unsigned npix = 0;
  unsigned npos = 0;
  double   samp = 0;
  double   swei = 0;
  double   sac1 = 0;
  double   sac2 = 0;
  double   sar1 = 0;
  double   sar2 = 0;

  unsigned rmin = std::max((int)m_win.rowmin, int(r0-rank));
  unsigned rmax = std::min((int)m_win.rowmax, int(r0+rank+1));
  unsigned cmin = std::max((int)m_win.colmin, int(c0-rank));
  unsigned cmax = std::min((int)m_win.colmax, int(c0+rank+1));

  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {

      if(!mask(r, c)) continue;
      double a = data(r, c) - peak.bkgd;
      npix += 1;
      samp += a;
      if(!(a>thr)) continue;
      npos += 1;
      swei += a;
      sar1 += a*r;
      sac1 += a*c;
      sar2 += a*r*r;
      sac2 += a*c*c;
    }
  }

  //if(!(samp>0)) return;
  //if(npos<npix-npos) return;
  //cout << "fraction of positive=" << float(npos)/npix << '\n';

  //peak.seg       = m_seg; // already set in _makePeaksFromMapOfLocalMaximumsV3
  //peak.row       = r0;
  //peak.col       = c0;
  peak.npix      = npos; // npix;
  peak.npos      = npos;
  peak.amp_max   = a0;

  if (swei>0 && npos>1) {
    peak.row_cgrav = sar1/swei;
    peak.col_cgrav = sac1/swei;
    peak.row_sigma = std::sqrt(sar2/swei - peak.row_cgrav * peak.row_cgrav);
    peak.col_sigma = std::sqrt(sac2/swei - peak.col_cgrav * peak.col_cgrav);  
  } 
  else {
    peak.row_cgrav = r0;
    peak.col_cgrav = c0;
    peak.row_sigma = 0;
    peak.col_sigma = 0;
  }

  peak.row_min   = rmin;
  peak.row_max   = rmax;
  peak.col_min   = cmin;
  peak.col_max   = cmax;  

  if (nsigm) {
    peak.amp_tot   = swei;  
    noise_tot = peak.noise * std::sqrt(npos);
  }
  else {
    peak.amp_tot   = samp;
    noise_tot = peak.noise * std::sqrt(npix);
  }

  peak.son = (noise_tot>0) ? peak.amp_tot / noise_tot : 0;

  //cout << "samp:" << samp << "  swei:" << swei << "  npos:" << npos  << "  noise_tot:" <<  noise_tot << "  son:" << peak.son << '\n';  //cout << "peak.row_sigma=" << peak.row_sigma << " col_sigma=" << peak.col_sigma << " samp=" << samp << '\n';
}

//--------------------
  /**
   * @brief _makePeaksFromMapOfLocalMaximums - makes peaks from the map of local maximums using rank as a peak size
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */
template <typename T>
void 
_makePeaksFromMapOfLocalMaximums( const ndarray<const T,2>& data
                                , const ndarray<const mask_t,2>& mask
                                , const size_t& rank
                                )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makePeaksFromMapOfLocalMaximums, seg=" << m_seg<< ", rank=" << rank);

  _reserveVectorOfPeaks();
  v_peaks.clear(); // this vector will be filled out for each window

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      if(m_local_maximums(r, c) & 4)
        _procLocalMaximum<T>(data,mask,rank,r,c);	
}

//--------------------
  /**
   * @brief The same as _makePeaksFromMapOfLocalMaximums, but uses _procLocalMaximumV2<T>(data,mask,rank,r,c);
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */
template <typename T>
void 
_makePeaksFromMapOfLocalMaximumsV2( const ndarray<const T,2>& data
                                  , const ndarray<const mask_t,2>& mask
                                  , const size_t& rank
                                  )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makePeaksFromMapOfLocalMaximumsV2, seg=" << m_seg<< ", rank=" << rank);

  _reserveVectorOfPeaks();
  v_peaks.clear(); // this vector will be filled out for each window

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      if(m_local_maximums(r, c) & 4)
        _procLocalMaximumV2<T>(data,mask,rank,r,c);	
}

//--------------------
  /**
   * @brief The same as _makePeaksFromMapOfLocalMaximums and V2, but fills a vector of seed peaks;
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */
template <typename T>
void 
_makePeaksFromMapOfLocalMaximumsV3( const ndarray<const T,2>& data
                                  , const ndarray<const mask_t,2>& mask
                                  , const size_t& rank
                                  )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makePeaksFromMapOfLocalMaximumsV3, seg=" << m_seg<< ", rank=" << rank);

  _reserveVectorOfPeaks();
  v_peaks.clear();

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      if(m_local_maximums(r, c) & 4) {

        Peak peak;        
        peak.seg       = m_seg;
        peak.row       = r;
        peak.col       = c;
        //peak.npix      = npix;
        //peak.bkgd      = 0; //sonres.avg;
        //peak.noise     = 0; //sonres.rms;
        //peak.son       = 0; //sonres.son;
        
        v_peaks.push_back(peak);
      }
}

//--------------------
  /**
   * @brief Significantly different from previous versions
   *        - now it does full peak processing in a right order
   *          1. estimate background
   *          2. finds connected pixels
   *          3. fills vector of peaks
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */
template <typename T>
void 
_makeMapOfConnectedPixelsForLocalMaximums( const ndarray<const T,2>& data
                                         , const ndarray<const mask_t,2>& mask
                                         , const size_t& rank
                                         , const float& nsigm = 0
                                         )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfConnectedPixelsForLocalMaximums, seg=" << m_seg<< ", rank=" << rank);

  // Is not used in this algorithm, but just in case if someone call alg.maps_of_pixel_status()
  if(m_pixel_status.empty()) {
    m_pixel_status = make_ndarray<pixel_status_t>(data.shape()[0], data.shape()[1]);
    std::fill_n(m_pixel_status.data(), int(data.size()), pixel_status_t(0));
  }

  if(v_ind_pixgrp.capacity() != m_pixgrp_max_size) v_ind_pixgrp.reserve(m_pixgrp_max_size);

  if(vv_peak_pixinds.capacity() != m_npksmax) vv_peak_pixinds.reserve(m_npksmax);
     vv_peak_pixinds.clear();

  if(v_peaks.capacity() != m_npksmax) v_peaks.reserve(m_npksmax);
     v_peaks.clear();

  if(m_conmap.empty())
     m_conmap = make_ndarray<conmap_t>(data.shape()[0], data.shape()[1]);
  std::fill_n(m_conmap.data(), int(data.size()), conmap_t(0));

  m_numreg=0;
  for(int r = (int)m_win.rowmin; r<(int)m_win.rowmax; r++)
    for(int c = (int)m_win.colmin; c<(int)m_win.colmax; c++) {

        if(m_numreg > m_npksmax-2) {
	  if(m_pbits) MsgLog(_name(), warning, "seg=" << m_seg << " Number of peaks exceeds reserved maximum " << m_npksmax);
          return; // protection
        }

        if(! (m_local_maximums(r, c) & 4)) continue;

        ++m_numreg;

        BkgdAvgRms bkgd = _evaluateBkgdAvgRmsV3<T>(data, mask, r, c);
        // cout << "YYY m_numreg=" << m_numreg << " r=" <<  r << " c=" <<  c << " bkgd:" << bkgd << '\n';
        m_reg_thr = bkgd.avg + nsigm * bkgd.rms;

        _findConnectedPixelsForLocalMaximum<T>(data, mask, rank, r, c); 
        if(v_ind_pixgrp.empty()) {--m_numreg; continue;}

        vv_peak_pixinds.push_back(v_ind_pixgrp);
	_procPixGroup<T>(data, bkgd, v_ind_pixgrp); // proc connected group and fills v_peaks
        //m_pixel_status[r][c] |= 32; // mark pixel as local maximum
    }
}

//--------------------
//--------------------
//--------------------
  /**
   * @brief Process peaks after S/N is evaluated for pfv3
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */

template <typename T>
void _procSeedPeaks( const ndarray<const T,2>& data
                   , const ndarray<const mask_t,2>& mask
                   , const size_t& rank
                   , const float& nsigm=0 // 0-turns off threshold algorithm, 1.64-leaves 5% of noise, etc.
                   )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _procSeedPeaks, seg=" << m_seg);

  std::vector<Peak>::iterator it;
  for(it=v_peaks.begin(); it!=v_peaks.end(); ++it) { 
    //Peak& peak = (*it);    
    _procLocalMaximumV3<T>(data, mask, rank, *it, nsigm);
  }
}

//--------------------
//--------------------
//--------------------

  /**
   * @brief Loops over list of peaks m_peaks, evaluates SoN info and adds it to each peak.
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  r0   - radial parameter of the ring for S/N evaluation algorithm
   * @param[in]  dr   - ring width for S/N evaluation algorithm
   */

template <typename T>
void _addBkgdAvgRmsToPeaks( const ndarray<const T,2>& data
                          , const ndarray<const mask_t,2>& mask
	                  , const float r0=7.0
	                  , const float dr=2.0
                          )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _addBkgdAvgRmsToPeaks, seg=" << m_seg);

  setSoNPars(r0, dr);

  std::vector<Peak>::iterator it;

  for(it=v_peaks.begin(); it!=v_peaks.end(); ++it) { 
    Peak& peak = (*it);
    
    BkgdAvgRms res = _evaluateBkgdAvgRms<T>((unsigned) peak.row, (unsigned) peak.col, data, mask);

    peak.bkgd  = res.avg;
    peak.noise = res.rms;
  }
}

//--------------------

  /**
   * @brief Loops over list of peaks m_peaks, evaluates SoN info and adds it to each peak.
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  r0   - radial parameter of the ring for S/N evaluation algorithm
   * @param[in]  dr   - ring width for S/N evaluation algorithm
   */

template <typename T>
void _addSoNToPeaks( const ndarray<const T,2>& data
                   , const ndarray<const mask_t,2>& mask
	           , const float r0=7.0
	           , const float dr=2.0
                   )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _addSoNToPeaks, seg=" << m_seg);

  setSoNPars(r0, dr);

  std::vector<Peak>::iterator it;

  for(it=v_peaks.begin(); it!=v_peaks.end(); ++it) { 
    Peak& peak = (*it);
    
    SoNResult sonres = _evaluateSoNForPixel<T>((unsigned) peak.row, (unsigned) peak.col, data, mask);

    peak.bkgd  = sonres.avg;
    peak.noise = sonres.rms;
    peak.son   = sonres.son;
  }
}

//--------------------

  /**
   * @brief The same as _addSoNToPeaks, but saves background-corrected amp_max, amp_tot, son
   */

template <typename T>
void _addSoNToPeaksV2( const ndarray<const T,2>& data
                     , const ndarray<const mask_t,2>& mask
	             , const float r0=7.0
	             , const float dr=2.0
                     )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _addSoNToPeaksV2, seg=" << m_seg);

  setSoNPars(r0, dr);

  std::vector<Peak>::iterator it;

  for(it=v_peaks.begin(); it!=v_peaks.end(); ++it) { 
    Peak& peak = (*it);
    
    SoNResult sonres = _evaluateSoNForPixel<T>((unsigned) peak.row, (unsigned) peak.col, data, mask);

    peak.amp_max -= sonres.avg;
    peak.amp_tot -= sonres.avg * peak.npix;

    peak.bkgd  = sonres.avg;
    peak.noise = sonres.rms;
    double noise_tot = sonres.rms * std::sqrt(peak.npix);
    peak.son   = (noise_tot>0) ? peak.amp_tot / noise_tot : 0;
  }
}

//--------------------
  /**
   * @brief _procDroplet - process a single droplet candidate
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr_high  - threshold on pixel intensity to be a candidate to "droplet"
   * @param[in]  rad - radius of the region where central pixel has a maximal value
   * @param[in]  r0 - droplet central pixel row-coordinate 
   * @param[in]  c0 - droplet central pixel column-coordinate   
   */

template <typename T>
void
_procDroplet( const ndarray<const T,2>& data
            , const ndarray<const mask_t,2>& mask
            , const T& thr_low
            , const unsigned& rad
            , const unsigned& r0
            , const unsigned& c0
            )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _procDroplet, seg=" << m_seg << " r0=" << r0 << " c0=" << c0);

  double a0 = data(r0,c0);
  unsigned npix = 0;
  double   samp = 0;
  double   sac1 = 0;
  double   sac2 = 0;
  double   sar1 = 0;
  double   sar2 = 0;

  unsigned rmin = std::max((int)m_win.rowmin, int(r0-rad));
  unsigned rmax = std::min((int)m_win.rowmax, int(r0+rad+1));
  unsigned cmin = std::max((int)m_win.colmin, int(c0-rad));
  unsigned cmax = std::min((int)m_win.colmax, int(c0+rad+1));

  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      double a = data(r, c);

      if(a>a0) return;  // This is not a local maximum inside rad...

      if(mask(r, c) && a>thr_low) {
	npix += 1;
	samp += a;
	sar1 += a*r;
	sac1 += a*c;
	sar2 += a*r*r;
	sac2 += a*c*c;
      }
    }
  }

  if(npix<1) return;

  Peak peak;

  peak.seg       = m_seg;
  peak.row       = r0;
  peak.col       = c0;
  peak.npix      = npix;
  peak.amp_max   = a0;
  peak.amp_tot   = samp;

  if(samp>0) {
    sar1 /= samp;
    sac1 /= samp;
    sar2 = sar2/samp - sar1*sar1;
    sac2 = sac2/samp - sac1*sac1;
    peak.row_cgrav = sar1;
    peak.col_cgrav = sac1;
    peak.row_sigma = (npix>1 && sar2>0) ? std::sqrt(sar2) : 0;
    peak.col_sigma = (npix>1 && sac2>0) ? std::sqrt(sac2) : 0;
  }
  else {
    peak.row_cgrav = r0;
    peak.col_cgrav = c0;
    peak.row_sigma = 0;
    peak.col_sigma = 0;
  }

  peak.row_min   = rmin;
  peak.row_max   = rmax;
  peak.col_min   = cmin;
  peak.col_max   = cmax;  
  peak.bkgd      = 0; //sonres.avg;
  peak.noise     = 0; //sonres.rms;
  peak.son       = 0; //sonres.son;

  if(_peakIsPreSelected(peak) && v_peaks.size()<m_npksmax-1) v_peaks.push_back(peak);
}

//--------------------
  /**
   * @brief _makeVectorOfDroplets - a part of peakFinderV1, V4 algorithm - loops over pixels and makes peak candidates in v_peaks
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr_low  - threshold on pixel intensity to be considered in this algorithm 
   * @param[in]  thr_high - threshold on pixel intensity to be a candidate to "droplet"
   * @param[in]  rad      - radius in pixels of squared region to find droplet relative to central pixel
   */
template <typename T>
void 
_makeVectorOfDroplets( const ndarray<const T,2>& data
                     , const ndarray<const mask_t,2>& mask
                     , const T& thr_low
                     , const T& thr_high
                     , const unsigned& rad=5
                     )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeVectorOfDroplets, seg=" << m_seg);

  _reserveVectorOfPeaks();
  v_peaks.clear(); // this vector will be filled out for each window

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      if(mask(r, c) && (data(r, c)>thr_high))
        _procDroplet<T>(data,mask,thr_low,rad,r,c);	
}

//--------------------
  /**
   * @brief _procPixGroup - process a pixel group for peakFinderV3r2
   * 
   * @param[in] data - ndarray with calibrated intensities
   * @param[in] bkgd - structure of background parameters
   * @param[in] vinds - vector of connected pixel indexes
   */

template <typename T>
void
_procPixGroup( const ndarray<const T,2>& data,
               const BkgdAvgRms& bkgd,
               const std::vector<TwoIndexes>& vinds
             )
{
  if(vinds.empty()) return;

  const int& r0 = vinds[0].i;
  const int& c0 = vinds[0].j;

  if(m_pbits & 512) MsgLog(_name(), info, "in _procPixGroup, seg=" << m_seg << " r0=" << r0 << " c0=" << c0);

  //BkgdAvgRms bkgd = _evaluateBkgdAvgRmsV2<T>(r0, c0, data);

  double   a0 = data(r0,c0) - bkgd.avg;
  unsigned npix = 0;
  double   samp = 0;
  double   sac1 = 0;
  double   sac2 = 0;
  double   sar1 = 0;
  double   sar2 = 0;
  int      rmin = data.shape()[0];
  int      cmin = data.shape()[1];
  int      rmax = 0;
  int      cmax = 0;

  for(vector<TwoIndexes>::const_iterator ij  = vinds.begin();
                                         ij != vinds.end(); ij++) {
      int r = ij->i;
      int c = ij->j;
      double a = data(r, c) - bkgd.avg;

      if(r<rmin) rmin=r;
      if(r>rmax) rmax=r;
      if(c<cmin) cmin=c;
      if(c>cmax) cmax=c;

      npix += 1;
      samp += a;
      sar1 += a*r;
      sac1 += a*c;
      sar2 += a*r*r;
      sac2 += a*c*c;
  }

  if(npix<1) return;

  Peak peak;

  peak.seg       = m_seg;
  peak.row       = r0;
  peak.col       = c0;
  peak.npix      = npix;
  peak.amp_max   = a0;
  peak.amp_tot   = samp;

  if(samp>0) {
    sar1 /= samp;
    sac1 /= samp;
    sar2 = sar2/samp - sar1*sar1;
    sac2 = sac2/samp - sac1*sac1;
    peak.row_cgrav = sar1;
    peak.col_cgrav = sac1;
    peak.row_sigma = (npix>1 && sar2>0) ? std::sqrt(sar2) : 0;
    peak.col_sigma = (npix>1 && sac2>0) ? std::sqrt(sac2) : 0;
  }
  else {
    peak.row_cgrav = r0;
    peak.col_cgrav = c0;
    peak.row_sigma = 0;
    peak.col_sigma = 0;
  }

  peak.row_min   = rmin;
  peak.row_max   = rmax;
  peak.col_min   = cmin;
  peak.col_max   = cmax;  
  peak.bkgd      = bkgd.avg;
  peak.noise     = bkgd.rms;
  //peak.bkgnpx   = bkgd.npx;
  double noise_tot = bkgd.rms * std::sqrt(npix);
  peak.son       = (noise_tot>0) ? peak.amp_tot / noise_tot : 0;

  v_peaks.push_back(peak);
}

//--------------------
//--------------------
  /**
   * @brief _procDropletV2 - process a single droplet candidate
   * 
   * @param[in] data - ndarray with calibrated intensities
   * @param[in] vinds - vector of connected pixel indexes
   */

template <typename T>
void
_procDropletV2( const ndarray<const T,2>& data,
                const std::vector<TwoIndexes>& vinds
              )
{
  if(vinds.empty()) return;

  const int& r0 = vinds[0].i;
  const int& c0 = vinds[0].j;

  if(m_pbits & 512) MsgLog(_name(), info, "in _procDropletV2, seg=" << m_seg << " r0=" << r0 << " c0=" << c0);

  BkgdAvgRms bkgd = _evaluateBkgdAvgRmsV2<T>(r0, c0, data);

  double   a0 = data(r0,c0);
  unsigned npix = 0;
  double   samp = 0;
  double   sac1 = 0;
  double   sac2 = 0;
  double   sar1 = 0;
  double   sar2 = 0;
  int      rmin = data.shape()[0];
  int      cmin = data.shape()[1];
  int      rmax = 0;
  int      cmax = 0;

  for(vector<TwoIndexes>::const_iterator ij  = vinds.begin();
                                         ij != vinds.end(); ij++) {
      int r = ij->i;
      int c = ij->j;
      double a = data(r, c);

      if(r<rmin) rmin=r;
      if(r>rmax) rmax=r;
      if(c<cmin) cmin=c;
      if(c>cmax) cmax=c;

      npix += 1;
      samp += a;
      sar1 += a*r;
      sac1 += a*c;
      sar2 += a*r*r;
      sac2 += a*c*c;
  }

  if(npix<1) return;

  Peak peak;

  peak.seg       = m_seg;
  peak.row       = r0;
  peak.col       = c0;
  peak.npix      = npix;
  peak.amp_max   = a0 - bkgd.avg;
  peak.amp_tot   = samp - bkgd.avg * npix;

  if(samp>0) {
    sar1 /= samp;
    sac1 /= samp;
    sar2 = sar2/samp - sar1*sar1;
    sac2 = sac2/samp - sac1*sac1;
    peak.row_cgrav = sar1;
    peak.col_cgrav = sac1;
    peak.row_sigma = (npix>1 && sar2>0) ? std::sqrt(sar2) : 0;
    peak.col_sigma = (npix>1 && sac2>0) ? std::sqrt(sac2) : 0;
  }
  else {
    peak.row_cgrav = r0;
    peak.col_cgrav = c0;
    peak.row_sigma = 0;
    peak.col_sigma = 0;
  }

  peak.row_min   = rmin;
  peak.row_max   = rmax;
  peak.col_min   = cmin;
  peak.col_max   = cmax;  
  peak.bkgd      = bkgd.avg;
  peak.noise     = bkgd.rms;
  //peak.bkgnpx   = bkgd.npx;
  double noise_tot = bkgd.rms * std::sqrt(npix);
  peak.son       = (noise_tot>0) ? peak.amp_tot / noise_tot : 0;

  v_peaks.push_back(peak);
}

//--------------------
  /**
   * @brief _makeVectorOfPeaksV2 - a part of peakFinderV4r2 algorithm 
   * loops over vector of peaks pixels and fill v_peaks
   * 
   * @param[in]  data - ndarray with calibrated intensities
   */
template <typename T>
void 
_makeVectorOfPeaksV2(const ndarray<const T,2>& data)
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeVectorOfPeaksV2, seg=" << m_seg);

  if(v_peaks.capacity() != m_npksmax) v_peaks.reserve(m_npksmax);
     v_peaks.clear();

  std::vector< std::vector<TwoIndexes> >::const_iterator it;
  for(it=vv_peak_pixinds.begin(); it!=vv_peak_pixinds.end(); ++it) { 
    //std::vector<TwoIndexes>& vinds = (*it);    
    _procDropletV2<T>(data, *it);
  }
}

//--------------------
  /**
   * @brief isDroplet - apply flood filling algorithms to find a group of connected pixels
   * around r0,c0 in constrained by rad region. 
   *   - uses m_pixel_status for 2-thresholds
   *   - check that r0,c0 is absolute maximum, returns false if not found
   *   - apply flood filling and make vector of connected pixels in rad constrained by region. 
   *   - returns true for absolute maximum.
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  rad - radius of the region where central pixel has a maximal value
   * @param[in]  r0 - droplet central pixel row-coordinate 
   * @param[in]  c0 - droplet central pixel column-coordinate   
   */

template <typename T>
bool
_isDroplet(const ndarray<const T,2>& data
         , const unsigned& rad
         , const int& r0
         , const int& c0
         )
{
  m_reg_rmin = std::max((int)m_win.rowmin, int(r0-rad));
  m_reg_rmax = std::min((int)m_win.rowmax, int(r0+rad+1));
  m_reg_cmin = std::max((int)m_win.colmin, int(c0-rad));
  m_reg_cmax = std::min((int)m_win.colmax, int(c0+rad+1));

  if(m_pbits & 512)
      MsgLog(_name(), info, "in _isDroplet, seg=" << m_seg 
             << " rad=" << rad  << " r0=" << r0 << " c0=" << c0
             << "\n    m_reg_rmin=" << m_reg_rmin
             << "\n    m_reg_rmax=" << m_reg_rmax
             << "\n    m_reg_cmin=" << m_reg_cmin
	     << "\n    m_reg_cmax=" << m_reg_cmax);

  v_ind_pixgrp.clear();

  // TEMPLATED CODE:
  m_reg_a0 = data(r0,c0);
  bool is_droplet = _findConnectedPixelsInRegionV1<T>(data, r0, c0);

  // NON-TEMPLATED CODE:
  /*
  _findConnectedPixelsInRegion(r0, c0);
  bool is_droplet = _isGroupMaximum<T>(data, r0, c0);
  */

  if(! is_droplet) _clearStatusOfUnusedPixels();
  return is_droplet;
}

//--------------------
//--------------------
  /**
   * @brief _findConnectedPixelsForLocalMaximum - apply flood filling algorithms to find a group of connected pixels
   * around r0,c0 in constrained by rad region. 
   *   - uses m_pixel_status for 2-thresholds
   *   - check that r0,c0 is absolute maximum, returns false if not found
   *   - apply flood filling and make vector of connected pixels in rad constrained by region. 
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rad - radius of the region where central pixel has a maximal value
   * @param[in]  r0 - droplet central pixel row-coordinate 
   * @param[in]  c0 - droplet central pixel column-coordinate   
   */

template <typename T>
void
_findConnectedPixelsForLocalMaximum(
           const ndarray<const T,2>& data
         , const ndarray<const mask_t,2>& mask
         , const unsigned& rank
         , const int& r0
         , const int& c0
         )
{
  m_reg_rmin = std::max((int)m_win.rowmin, int(r0-rank));
  m_reg_rmax = std::min((int)m_win.rowmax, int(r0+rank+1));
  m_reg_cmin = std::max((int)m_win.colmin, int(c0-rank));
  m_reg_cmax = std::min((int)m_win.colmax, int(c0+rank+1));

  //if(m_pbits & 512) MsgLog(_name(), info, "in _findConnectedPixelsForLocalMaximum, seg=" << m_seg 
  //                          << " rank=" << rank  << " r0=" << r0 << " c0=" << c0);
  //cout << "ZZZ _findConnectedPixelsForLocalMaximum : rank=" << rank  << " r0=" << r0 << " c0=" << c0 << '\n';
  //cout << "ZZZ: m_reg_rmin: " << m_reg_rmin << "  m_reg_rmax: " << m_reg_rmax
  //        << "  m_reg_cmin: " << m_reg_cmin << "  m_reg_cmax: " << m_reg_cmax << '\n';

  v_ind_pixgrp.clear();
  _findConnectedPixelsInRegionV2<T>(data, mask, r0, c0);
}

//--------------------
  /**
   * @brief _isGroupMaximum - check if the data in r0,c0 is a group maximum.
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  r0 - droplet central pixel row-coordinate 
   * @param[in]  c0 - droplet central pixel column-coordinate   
   */
template <typename T>
bool
_isGroupMaximum( const ndarray<const T,2>& data
               , const int& r0
               , const int& c0) {
  double a0 = data(r0,c0);
  for(vector<TwoIndexes>::const_iterator ij  = v_ind_pixgrp.begin();
                                         ij != v_ind_pixgrp.end(); ij++)
    if(data(ij->i,ij->j) > a0) return false; // is not local maximum

  return true; // is local maximum
}

//--------------------
  /**
   * @brief _makeMapOfConnectedPixelsForDroplets - a part of peakFinderV4r2 algorithm, 
   *        replaces _makeVectorOfDroplets
   * Changes: 
   *   1. droplet includes connected pixels ONLY!
   *   2. makes map of connected pixels for droplets (in stead of vector of peaks
   *   3. S/N is evaluated in the droplet processing
   * 
   * @param[in] data - ndarray with calibrated intensities
   * @param[in] rank - radius in pixels of squared region to find droplet relative to central pixel
   */
template <typename T>
void 
_makeMapOfConnectedPixelsForDroplets(const ndarray<const T,2>& data, const unsigned& rank=5)
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfConnectedPixelsForDroplets, seg=" << m_seg);

  //if(m_conmap.size()==0) 
  if(m_conmap.empty())
     m_conmap = make_ndarray<conmap_t>(data.shape()[0], data.shape()[1]);
  std::fill_n(m_conmap.data(), int(data.size()), conmap_t(0));

  //m_pixgrp_max_size = (2*rank+1)*(2*rank+1); //already defined in calling method
  if(v_ind_pixgrp.capacity() != m_pixgrp_max_size) v_ind_pixgrp.reserve(m_pixgrp_max_size);

  if(vv_peak_pixinds.capacity() != m_npksmax) vv_peak_pixinds.reserve(m_npksmax);
     vv_peak_pixinds.clear();

  if(v_peaks.capacity() != m_npksmax) v_peaks.reserve(m_npksmax);
     v_peaks.clear();

  m_numreg=0;
  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {

      if(m_numreg > m_npksmax-2) {
	if(m_pbits) MsgLog(_name(), warning, "seg=" << m_seg << " Number of peaks exceeds reserved maximum " << m_npksmax);
        return; // protection
      }

      if(  m_pixel_status(r, c) & 1
      && !(m_pixel_status(r, c) & 16))
      // mask=1 and a>thr_high and pixel is not used
        if(_isDroplet<T>(data, rank, (int)r, (int)c)) {
            ++m_numreg;
            _addPixGroupToMap(); // fills m_conmap and m_pixel_status[r][c] |= 16;
            _procDropletV2<T>(data, v_ind_pixgrp); // proc bkgd, connected group and fills v_peaks
            vv_peak_pixinds.push_back(v_ind_pixgrp);
            m_pixel_status(r, c) |= 32; // mark pixel as local maximum
	}
    }
}

//--------------------

template <typename T>
void
_init_local_extreams_arrays(const ndarray<const T,2>& data) 
{
  if(m_local_maximums.empty()) {
     m_local_maximums = make_ndarray<pixel_maximums_t>(data.shape()[0], data.shape()[1]);
     std::fill_n(&m_local_maximums(0,0), int(data.size()), pixel_maximums_t(0));
  }

  if(m_local_minimums.empty()) {
     m_local_minimums = make_ndarray<pixel_minimums_t>(data.shape()[0], data.shape()[1]);
     std::fill_n(&m_local_minimums(0,0), int(data.size()), pixel_minimums_t(0));
  }
}

//--------------------
//--------------------
//--------------------

public:

//--------------------
  /**
   * @brief numberOfPixAboveThr - counts a number of pixels above threshold
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr  - threshold on data values
   */

template <typename T>
unsigned
numberOfPixAboveThr( const ndarray<const T,2>& data
                   , const ndarray<const mask_t,2>& mask
                   , const T& thr
                   )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in numberOfPixAboveThr, seg=" << m_seg);

  m_win.validate(data.shape());

  unsigned npix = 0;
  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      if(mask(r, c) && (data(r, c)>thr)) npix++;
    }
  }
  return npix;
}

//--------------------
  /**
   * @brief intensityOfPixAboveThr - evaluates total intensity of pixels above threshold
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr  - threshold on data values
   */

template <typename T>
double
intensityOfPixAboveThr( const ndarray<const T,2>& data
                      , const ndarray<const mask_t,2>& mask
                      , const T& thr
                      )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in intensityOfPixAboveThr, seg=" << m_seg);

  m_win.validate(data.shape());

  double amptot = 0;
  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      if(mask(r, c) && (data(r, c)>thr)) amptot += (double)data(r, c);
    }
  }
  return amptot;
}

//--------------------
  /**
   * @brief peakFinderV1 - "Droplet-finder" - two-threshold peak finding algorithm in the region defined by the radial parameter
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr_low   - threshold on pixel intensity to be considered in this algorithm 
   * @param[in]  thr_high  - threshold on pixel intensity to be a candidate to "droplet"
   * @param[in]  rad       - radius in pixels of squared region to find droplet relative to central pixel
   * @param[in]  dr        - width of the ring of radius rad for SoN algorithm
   */

template <typename T>
std::vector<Peak>&
peakFinderV1( const ndarray<const T,2>& data
            , const ndarray<const mask_t,2>& mask
            , const T& thr_low
            , const T& thr_high
            , const unsigned& rad=5
            , const float& dr=2.0
            )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV1, seg=" << m_seg);

  m_win.validate(data.shape());

  _makeVectorOfDroplets<T>(data, mask, thr_low, thr_high, rad);
  _addSoNToPeaks<T>(data, mask, float(rad), dr);
  _makeVectorOfSelectedPeaks();
  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief peakFinderV4 - "Droplet-finder" - two-threshold peak finding algorithm in the region defined by the rank parameter
   * peakFinderV4 has rank and r0 parameters in stead of common rad like in peakFinderV2,V3
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr_low   - threshold on pixel intensity to be considered in this algorithm 
   * @param[in]  thr_high  - threshold on pixel intensity to be a candidate to "droplet"
   * @param[in]  rank      - radius in pixels of squared region to find droplet relative to central pixel
   * @param[in]  r0        - radius for SoN algorithm
   * @param[in]  dr        - width of the ring of radius rad for SoN algorithm
   */

template <typename T>
std::vector<Peak>&
peakFinderV4( const ndarray<const T,2>& data
            , const ndarray<const mask_t,2>& mask
            , const T& thr_low
            , const T& thr_high
            , const unsigned& rank = 5
            , const float& r0 = 7.0
            , const float& dr = 2.0
            )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV4, seg=" << m_seg);

  m_win.validate(data.shape());

  _makeVectorOfDroplets<T>(data, mask, thr_low, thr_high, rank);
  _addSoNToPeaks<T>(data, mask, r0, dr);
  _makeVectorOfSelectedPeaks();
  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief peakFinderV4r1 - "Droplet-finder" - the same as V4, but returns background-corrected amp_max, amp_tot, son (total).
   * Changes: 
   *   m_do_preselect = false; 
   *   _addSoNToPeaksV2<T>(data, mask, r0, dr);
   */

template <typename T>
std::vector<Peak>&
peakFinderV4r1( const ndarray<const T,2>& data
              , const ndarray<const mask_t,2>& mask
              , const T& thr_low
              , const T& thr_high
              , const unsigned& rank = 5
              , const float& r0 = 7.0
              , const float& dr = 2.0
              )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV4r1, seg=" << m_seg);

  m_win.validate(data.shape());
  m_do_preselect = false; 

  _makeVectorOfDroplets<T>(data, mask, thr_low, thr_high, rank);
  _addSoNToPeaksV2<T>(data, mask, r0, dr);
  _makeVectorOfSelectedPeaks();
  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief peakFinderV4r2 - "Droplet-finder" - further improvement of V4r1.
   * Changes on 2016-10-31 Happy Helloween!: 
   *   eveluate S/N earlier - for each seed peak at processing of the droplet,
   *   peak accounts for connected pixels above thr_low only, using flood filling.
   * Details:
   *   - makes map of pixel status m_pixel_status for 2 tresholds and mask
   *   - process droplets in single plase, accounting for S/N and connected pixels only
   *   - makes final vector of selected peaks
   *   
   */

template <typename T>
std::vector<Peak>&
peakFinderV4r2( const ndarray<const T,2>& data
              , const ndarray<const mask_t,2>& mask
              , const T& thr_low
              , const T& thr_high
              , const unsigned& rank = 5
              , const float& r0 = 7.0
              , const float& dr = 2.0
              )
{
  m_r0 = r0;
  m_dr = dr;
  m_pixgrp_max_size = (2*rank+1)*(2*rank+1);

  m_win.validate(data.shape());

  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV4r2, seg=" << m_seg);

  _makeMapOfPixelStatusV2<T>(data, mask, thr_low, thr_high);
  _makeMapOfConnectedPixelsForDroplets<T>(data, rank);
  //_makeVectorOfPeaksV2<T>(data); //_procDropletV2 moved to makeMapOfConnectedPixelsForDroplets
  _makeVectorOfSelectedPeaks();

  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief peakFinderV2 - "Flood filling" - makes a list of peaks for groups of connected pixels above threshold.
   * 1) uses data and mask and finds groups of connected pixels with innensity above threshold;
   * 2) each group of connected pixels is processed as a single peak,
   *    its parameters and correlators are collected in struct PeakWork,
   *    all of them are collected in the std::vector<PeakWork> v_peaks_work;
   * 3) v_peaks_work is processed and the list of peak parameters is saved in std::vector<PeakWork> v_peaks_work.
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr  - threshold on data values
   * @param[in]  r0   - radius for SoN algorithm
   * @param[in]  dr   - width of the ring of radius rad for SoN algorithm
   */

template <typename T>
std::vector<Peak>&
peakFinderV2( const ndarray<const T,2>&      data
            , const ndarray<const mask_t,2>& mask
            , const T& thr
	    , const float r0 = 7.0
	    , const float dr = 2.0
            )
{
  m_win.validate(data.shape());

  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV2, seg=" << m_seg << " win" << m_win);
  _makeMapOfPixelStatus<T>(data, mask, thr);
  _makeMapOfConnectedPixels();
  _procConnectedPixels<T>(data);
  _makeVectorOfPeaks();
  _addSoNToPeaks<T>(data, mask, r0, dr);
  _makeVectorOfSelectedPeaks();

  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief peakFinderV2r1 - "Flood filling" - the same as V2, but returns background-corrected amp_max, amp_tot, son (total).
   * Changes: 
   *   m_do_preselect = false; 
   *   _addSoNToPeaksV2<T>(data, mask, r0, dr);
   */

template <typename T>
std::vector<Peak>&
peakFinderV2r1( const ndarray<const T,2>&      data
              , const ndarray<const mask_t,2>& mask
              , const T& thr
	      , const float r0 = 7.0
	      , const float dr = 2.0
              )
{
  m_win.validate(data.shape());

  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV2r1, seg=" << m_seg << " win" << m_win);
  _makeMapOfPixelStatus<T>(data, mask, thr);
  _makeMapOfConnectedPixels();
  _procConnectedPixels<T>(data);
  m_do_preselect = false; 
  _makeVectorOfPeaks();
  _addSoNToPeaksV2<T>(data, mask, r0, dr);
  _makeVectorOfSelectedPeaks();

  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief peakFinderV3 - "Ranker" - makes a list of peaks for local maximums of requested rank.
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   * @param[in]  r0   - radius for SoN algorithm
   * @param[in]  dr   - width of the ring of radius rad for SoN algorithm
   */

template <typename T>
std::vector<Peak>&
peakFinderV3( const ndarray<const T,2>&      data
            , const ndarray<const mask_t,2>& mask
            , const size_t rank = 2
	    , const float r0 = 7.0
	    , const float dr = 2.0
            )
{
  m_win.validate(data.shape());

  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV3, seg=" << m_seg << " win" << m_win << " rank=" << rank);

  _makeMapOfLocalMaximums<T>(data, mask, rank);
  _makePeaksFromMapOfLocalMaximums<T>(data, mask, rank);
  _addSoNToPeaks<T>(data, mask, r0, dr);
  _makeVectorOfSelectedPeaks();

  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief peakFinderV3r1 - "Ranker" - the same as V3, evaluate peak info after SoN algorithm.
   * Changes: 
   *   _makePeaksFromMapOfLocalMaximumsV3<T>(data, mask, rank); - makes list of empty seed peaks
   *   _addBkgdAvgRmsToPeaks<T>(data, mask, r0, dr); - use it in stead of of SoN algorithm
   *   _procSeedPeaks<T>(data, mask, rank, nsigm);
   *
   *  nsigm=0-turns off threshold algorithm, 1.64-leaves 5% of noise, etc.; 
   */

template <typename T>
std::vector<Peak>&
peakFinderV3r1( const ndarray<const T,2>& data
              , const ndarray<const mask_t,2>& mask
              , const size_t rank = 5
	      , const float r0 = 7.0
	      , const float dr = 2.0
	      , const float nsigm = 0 // 0-turns off threshold algorithm, 1.64-leaves 5% of noise, etc.; 
              )
{
  m_win.validate(data.shape());

  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV3r1, seg=" << m_seg << " win" << m_win << " rank=" << rank);

  _makeMapOfLocalMaximums<T>(data, mask, rank);
  _makeMapOfLocalMinimums<T>(data, mask, rank);
  if(m_pbits & 8) _printStatisticsOfLocalExtremes();

  m_do_preselect = false;
  _makePeaksFromMapOfLocalMaximumsV3<T>(data, mask, rank); // makes vector of seed peaks
  _addBkgdAvgRmsToPeaks<T>(data, mask, r0, dr);            // adds background average and rms to peak
  _procSeedPeaks<T>(data, mask, rank, nsigm);              // process seed peaks
  _makeVectorOfSelectedPeaks();                            // make vector of selected peaks

  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief peakFinderV3r2 - "Ranker" - further development of V3r1.
   * Changes: 
   * - keep algorithms for local extremums definition
   * - change algorithms for seed peak processing;
   *   - change order, make it straightforward without repeating loops over peaks: 
   *     - 1st evaluate background, 
   *     - then find pixels associated with peak
   *   - use connected pixels only
   *
   *  nsigm=0-turns off threshold algorithm, 1.64-leaves 5% of noise, etc.; 
   */

template <typename T>
std::vector<Peak>&
peakFinderV3r2( const ndarray<const T,2>& data
              , const ndarray<const mask_t,2>& mask
              , const size_t rank = 5
	      , const float r0 = 7.0
	      , const float dr = 2.0
	      , const float nsigm = 0 // 0-turns off threshold algorithm, 1.64-leaves 5% of noise, etc.; 
              )
{
  m_r0 = r0;
  m_dr = dr;
  m_pixgrp_max_size = (2*rank+1)*(2*rank+1);

  m_win.validate(data.shape());
  //std::fill_n(m_pixel_status.data(), int(m_pixel_status.size()), pixel_status_t(0));

  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV3r2, seg=" << m_seg << " win" << m_win << " rank=" << rank);
  _makeMapOfLocalMaximums<T>(data, mask, rank); // fills m_local_maximums
  _makeMapOfLocalMinimums<T>(data, mask, rank); // fills m_local_minimums
  if(m_pbits & 8) _printStatisticsOfLocalExtremes();
  _makeMapOfConnectedPixelsForLocalMaximums<T>(data, mask, rank, nsigm); // fills m_conmap, v_peaks, vv_peak_pixinds
  _makeVectorOfSelectedPeaks();                            // make vector of selected peaks
  return v_peaks_sel; 
}

//--------------------
  /**
   * @brief Evaluate per-pixel result of the S/N (median) algorithm to data using mask
   * 
   * S/N is evaluated for any pixel specified by the (row,col).  
   * If mask is provided, and pixel is masked (0) then default result is returned.
   * S/N algorithm uses non-masked surrounding pixels in the ring m_r0, m_dr.
   * Thresholds are not applied in order to prevent offset of the average value of the background level.
   * 
   * @param[in]  row  - pixel row
   * @param[in]  col  - pixel column
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   */

template <typename T>
SoNResult
_evaluateSoNForPixel( const unsigned& row
                    , const unsigned& col
                    , const ndarray<const T,2>& data
                    , const ndarray<const mask_t,2>& mask
                    )
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _evaluateSoNForPixel, seg=" << m_seg << " row=" << row << ", col=" << col);

  // S/N algorithm initialization
  if(! m_init_son_is_done) {
    _evaluateRingIndexes(m_r0, m_dr);
    m_win.validate(data.shape());
    m_use_mask = (mask.empty()) ? false : true;
    m_init_son_is_done = true;
  }

  if(m_use_mask && (!mask(row,col))) return m_sonres_def;
  //if(m_use_mask && (!mask(row,col))) return SoNResult({});

  double   amp  = 0;
  unsigned sum0 = 0;
  double   sum1 = 0;
  double   sum2 = 0;

  for(vector<TwoIndexes>::const_iterator ij  = v_indexes.begin();
                                         ij != v_indexes.end(); ij++) {
    int ir = row + (ij->i);
    int ic = col + (ij->j);

    if(ic < (int)m_win.colmin || !(ic < (int)m_win.colmax)) continue;
    if(ir < (int)m_win.rowmin || !(ir < (int)m_win.rowmax)) continue;
    if(m_use_mask && (! mask(ir,ic))) continue;

    //amp = (double)data(ir,ic);
    amp = (double)data(ir,ic);
    sum0 ++;
    sum1 += amp;
    sum2 += amp*amp;
  }
  //SoNResult res = m_sonres_def;
  SoNResult res;

  if(sum0) {
    res.avg = sum1/sum0;                              // Averaged background level
    double w = sum2/sum0 - res.avg*res.avg;
    res.rms = (w>0) ? std::sqrt(w) : 0;               // RMS of the background around peak
    res.sig = data(row,col) - res.avg;               // Signal above the background
    if (res.rms>0) res.son = res.sig/res.rms;         // S/N ratio
  }

  return res;
}

//--------------------
  /**
   * @brief Evaluate background average and rms in the ring around pixel using data and mask
   * 
   * Background average and rms are evaluated for any pixel specified by the (row,col).  
   * If mask is provided, and pixel is masked (0) then default result is returned.
   * This algorithm uses non-masked surrounding pixels in the ring m_r0, m_dr.
   * Thresholds are not applied in order to prevent offset of the average value of the background level.
   * 
   * @param[in]  row  - pixel row
   * @param[in]  col  - pixel column
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   */

template <typename T>
BkgdAvgRms
_evaluateBkgdAvgRms( const unsigned& row
                   , const unsigned& col
                   , const ndarray<const T,2>& data
                   , const ndarray<const mask_t,2>& mask
                   )
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _evaluateBkgdAvgRms, seg=" << m_seg << " row=" << row << ", col=" << col);

  // background evaluation algorithm initialization
  if(! m_init_bkgd_is_done) {
    _evaluateRingIndexes(m_r0, m_dr);
    m_win.validate(data.shape());
    m_use_mask = (mask.empty()) ? false : true;
    m_init_bkgd_is_done = true;
  }

  if(m_use_mask && (!mask(row,col))) return m_bkgdavgrms_def;

  double   amp  = 0;
  unsigned sum0 = 0;
  double   sum1 = 0;
  double   sum2 = 0;

  for(vector<TwoIndexes>::const_iterator ij  = v_indexes.begin();
                                         ij != v_indexes.end(); ij++) {
    int ir = row + (ij->i);
    int ic = col + (ij->j);

    if(ic < (int)m_win.colmin || !(ic < (int)m_win.colmax)) continue;
    if(ir < (int)m_win.rowmin || !(ir < (int)m_win.rowmax)) continue;
    if(m_use_mask && (! mask(ir,ic))) continue;
    if(m_local_maximums(ir,ic)) continue; // discard all types of local maximums from evaluation of bkgd
    if(m_local_minimums(ir,ic)) continue; // discard all types of local minimums from evaluation of bkgd

    amp = (double)data(ir,ic);
    sum0 ++;
    sum1 += amp;
    sum2 += amp*amp;
  }

  BkgdAvgRms res; // m_bkgdavgrms_def;

  if(sum0) {
    res.avg = sum1/sum0;                              // Averaged background level
    double w = sum2/sum0 - res.avg*res.avg;
    res.rms = (w>0) ? std::sqrt(w) : 0;               // RMS of the background around peak
    res.npx = sum0;
    //cout << res << '\n';
  }

  return res;
}

//--------------------
  /**
   * @brief Evaluate background average and rms in the ring around pixel using data and m_pixel_status.
   * V2 - uses m_pixel_status in V4r2 in stead of mask.
   * 
   * Background average and rms are evaluated for any pixel specified by the (row,col).  
   * Uses m_pixel_status in V4r2 in stead of mask. Good pixels with intensity below thresholds are used.
   * This algorithm uses pixels in the ring m_r0, m_dr.
   * 
   * @param[in] row  - pixel row
   * @param[in] col  - pixel column
   * @param[in] data - ndarray with calibrated intensities
   */

template <typename T>
BkgdAvgRms
_evaluateBkgdAvgRmsV2( const int& row
                     , const int& col
                     , const ndarray<const T,2>& data
                     )
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _evaluateBkgdAvgRmsV2, seg=" << m_seg << " row=" << row << ", col=" << col);

  // background evaluation algorithm initialization
  if(! m_init_bkgd_is_done) {
    _evaluateRingIndexes(m_r0, m_dr);
    //printVectorOfRingIndexes();
    //printMatrixOfRingIndexes();
    m_win.validate(data.shape());
    //m_win.print();
    m_init_bkgd_is_done = true;
  }

  // for masked pixel
  //if(!m_pixel_status(row,col)) return m_bkgdavgrms_def;

  double   amp  = 0;
  unsigned sum0 = 0;
  double   sum1 = 0;
  double   sum2 = 0;

  for(vector<TwoIndexes>::const_iterator ij  = v_indexes.begin();
                                         ij != v_indexes.end(); ij++) {
    int ir = row + (ij->i);
    int ic = col + (ij->j);

    if(ic < (int)m_win.colmin || !(ic < (int)m_win.colmax)) continue;
    if(ir < (int)m_win.rowmin || !(ir < (int)m_win.rowmax)) continue;

    // pixel selector:
    if(!(m_pixel_status(ir,ic) & 4)) continue; // uses pixels below both thresholds ONLY!

    amp = (double)data(ir,ic);
    sum0 ++;
    sum1 += amp;
    sum2 += amp*amp;
  }

  BkgdAvgRms res; // m_bkgdavgrms_def;

  if(sum0) {
    res.avg = sum1/sum0;                    // Averaged background level
    double w = sum2/sum0 - res.avg*res.avg;
    res.rms = (w>0) ? std::sqrt(w) : 0;     // RMS of the background around peak
    res.npx = sum0;
    //cout << "ZZZ: " << res << '\n';
  }

  return res;
}

//--------------------
//--------------------
  /**
   * @brief Evaluate background average and rms in the ring around pixel using data, mask, and map of local intensity extremes.
   * V3 for peakFinderV3r2 - rejects pixels using mask and map of local intensity extremes.
   * 
   * Background average and rms are evaluated for any pixel specified by the (row,col).  
   * Uses mask. Good non-extreme pixels are used.
   * This algorithm uses pixels in the ring m_r0, m_dr.
   * 
   * @param[in]  row  - pixel row
   * @param[in]  col  - pixel column
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   */

template <typename T>
BkgdAvgRms
_evaluateBkgdAvgRmsV3( const ndarray<const T,2>& data
                     , const ndarray<const mask_t,2>& mask
		     , const int& row
                     , const int& col
                     )
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _evaluateBkgdAvgRmsV3, seg=" << m_seg << " row=" << row << ", col=" << col);

  // background evaluation algorithm initialization
  if(! m_init_bkgd_is_done) {
    _evaluateRingIndexes(m_r0, m_dr);
    //printVectorOfRingIndexes();
    //printMatrixOfRingIndexes();
    m_win.validate(data.shape());
    //m_win.print();
    m_init_bkgd_is_done = true;
  }

  double   amp  = 0;
  unsigned sum0 = 0;
  double   sum1 = 0;
  double   sum2 = 0;

  for(vector<TwoIndexes>::const_iterator ij  = v_indexes.begin();
                                         ij != v_indexes.end(); ij++) {
    int ir = row + (ij->i);
    int ic = col + (ij->j);

    //cout << "YYY: " << " row=" << row << " col=" << col << " ir=" << ir << " ic=" << ic << '\n';

    if(ic < (int)m_win.colmin || !(ic < (int)m_win.colmax)) continue;
    if(ir < (int)m_win.rowmin || !(ir < (int)m_win.rowmax)) continue;

    // pixel selector:
    if(! mask(ir,ic)) continue;               // skip masked pixels
    if(m_local_maximums(ir,ic) & 7) continue; // skip extremal pixels
    if(m_local_minimums(ir,ic) & 7) continue; // skip extremal pixels (only 3 for mins implemented)

    amp = (double)data(ir,ic);
    sum0 ++;
    sum1 += amp;
    sum2 += amp*amp;
  }

  BkgdAvgRms bkgd; // object with 0 values

  if(sum0) {
    sum1 /= sum0;
    sum2 = sum2/sum0 - sum1*sum1;
    bkgd.avg = sum1;                           // Averaged background level
    bkgd.rms = (sum2>0) ? std::sqrt(sum2) : 0; // RMS of the background around peak
    bkgd.npx = sum0;                           // Number of accounted pixels
    //cout << "ZZZ _evaluateBkgdAvgRmsV3:" << bkgd << '\n';
  }
  return bkgd;
}

//--------------------
//--------------------
  /**
   * @brief Get ALL results of the S/N algorithm applied to data ndarray using mask
   * 
   * @param[in]  data   - ndarray with calibrated intensities
   * @param[in]  mask   - ndarray with mask of bad/good (0/1) pixels
   * @param[out] result - ndarray with results of median algorithm (average bkgd, rms, signal, S/N)
   * @param[in]  do_fill_def - pre-fill ndarray with default values
   */

template <typename T>
void getSoNResult( const ndarray<const T,2>& data
                 , const ndarray<const mask_t,2>& mask
                 , ndarray<SoNResult,2>& result
                 , const bool& do_fill_def = false
                 )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in getSoNResult, seg=" << m_seg);

  if(do_fill_def) std::fill_n(&result, int(data.size()), m_sonres_def);

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      result(r, c) = _evaluateSoNForPixel<T>(r, c, data, mask);
    }
  }
}

//--------------------
  /**
   * @brief Get ONLY S/N for data ndarray using mask
   * 
   * @param[in]  data   - ndarray with calibrated intensities
   * @param[in]  mask   - ndarray with mask of bad/good (0/1) pixels
   * @param[out] son - ndarray with results of median algorithm (average bkgd, rms, signal, S/N)
   * @param[in]  do_fill_def - pre-fill ndarray with default values
   */

template <typename T>
void getSoN( const ndarray<const T,2>& data
           , const ndarray<const mask_t,2>& mask
           , ndarray<son_t,2>& son
           , const bool& do_fill_def = false
           )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in getSoN, seg=" << m_seg);

  if(do_fill_def) std::fill_n(&son, int(data.size()), son_t(0));

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      son(r, c) = _evaluateSoNForPixel<T>(r, c, data, mask).son;
    }
  }
}

//--------------------
//--------------------
//--------------------
  /**
   * @brief Split calibrated data representing number of photons for integer (floor) and float fractional arrays.
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   */

template <typename T>
void
_splitDataForUintAndFloat( const ndarray<const T,2>&      data
                         , const ndarray<const u8mask_t,2>& mask
                         )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _splitDataForUintAndFloat, seg=" << m_seg << "\n    in window: " << m_win);
  //if(m_pbits & 512) m_win.print();

  if(m_nphoton.empty()) 
     m_nphoton = make_ndarray<nphoton_t>(data.shape()[0], data.shape()[1]);
     m_fphoton = make_ndarray<fphoton_t>(data.shape()[0], data.shape()[1]);

  std::fill_n(&m_nphoton(0,0), int(data.size()), nphoton_t(0));
  std::fill_n(&m_fphoton(0,0), int(data.size()), fphoton_t(0));

  unsigned rmin = (unsigned)m_win.rowmin;
  unsigned rmax = (unsigned)m_win.rowmax;				  
  unsigned cmin = (unsigned)m_win.colmin;
  unsigned cmax = (unsigned)m_win.colmax;

  // 
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      if(m_use_mask && (! mask(r, c))) continue; // leave 0-s

      T v = data(r, c);
      if (v>0) {
        m_nphoton(r, c) = (nphoton_t)floor(v);
        m_fphoton(r, c) = (fphoton_t)(v - m_nphoton(r, c));
      }
    }
  }
}

//--------------------
  /**
   * @brief mapOfPhotonNumbersV1 - Chuck's photon counting algorithm - apply fancy correction for split photons.
   * 
   * 1) splits calibrated data for uint (floor) and float leftover fractional number of photons
   * 2) merge fractional number of photons to largest intensity integer  
   * 3) sum together uint and merged fractional maps
   * 
   * Returns array with (uint16) number of photons per pixel from input array of calibrated intensities.
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   */

template <typename T>
ndarray<nphoton_t, 2>& 
mapOfPhotonNumbersV1( const ndarray<const T,2>&      data
                    , const ndarray<const u8mask_t,2>& mask
                    , const float& thr_fraction=0.9
                    )
{
  m_win.validate(data.shape());

  if(m_pbits & 512) MsgLog(_name(), info, "in mapOfPhotonNumbersV1, seg=" << m_seg << "\n    in window: " << m_win);

  _splitDataForUintAndFloat<T>(data, mask);

  //cout << "XXX thr_fraction " << thr_fraction << '\n';
  //cout << "XXX m_fphoton(1) " << m_fphoton << '\n';
  //cout << "XXX m_nphoton(1) " << m_nphoton << '\n';

  //size_t rank=1;
  //_makeMapOfLocalMaximums<fphoton_t>(m_fphoton, mask, rank);

  _makeMapOfLocalMaximumsRank1Cross<fphoton_t>(m_fphoton);

  //cout << "XXX m_local_maximums(1) " << m_local_maximums<< '\n';

  const fphoton_t thr_on_max = 0.5; const fphoton_t thr_on_tot = thr_fraction; const bool DO_TEST = false;

  _mergeConnectedPixelCouples(thr_on_max, thr_on_tot, DO_TEST); // DO_TEST fills m_mphoton

  //cout << "XXX m_conmap     " << m_conmap << '\n';
  //cout << "XXX m_mphoton    " << m_mphoton << '\n';
  //cout << "XXX m_nphoton(2) " << m_nphoton << '\n';

  return m_nphoton; 
  //return m_conmap
  //return m_mphoton;
  //return m_local_maximums;
}

//--------------------
//--------------------

};

//--------------------
//--------------------
// Wrappers for 2-d methods
//--------------------
//--------------------
  /**
   * @brief Wrapper for AlgImgProc::mapOfPhotonNumbersV1.
   * 
   * Returns array with (uint16) number of photons per pixel from input array of calibrated intensities.
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   */

template <typename T>
ndarray<const AlgImgProc::nphoton_t, 2>
mapOfPhotonNumbersV1( const ndarray<const T,2> data
		    , const ndarray<const AlgImgProc::u8mask_t,2> mask
                    , const float& thr_fraction=0.9
                    )
{
  size_t      seg  = 0;
  AlgImgProc algo(seg); //, 0, 100000, 0, 100000, 017777);

  return algo.mapOfPhotonNumbersV1<T>(data, mask, thr_fraction);
}

//--------------------
  /**
   * @brief Makes map of local maximums of requested rank
   * 
   * Map of local maximums is a 2-d array of unsigned integer values of data shape, 
   * it has 0/1/2/4 bits for non-maximum / maximum in column / maximum in row / local maximum in rectangle of radius rank.   
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */

template <typename T>
ndarray<const AlgImgProc::pixel_maximums_t, 2> 
mapOfLocalMaximums( const ndarray<const T,2> data
                  , const ndarray<const AlgImgProc::mask_t,2> mask
                  , const size_t& rank
                  )
{
  AlgImgProc algo(0); // , 0, 1e6, 0, 1e6, 1023);
  algo.validate_window(data.shape());
  algo._makeMapOfLocalMaximums<T>(data, mask, rank);
  return algo.mapOfLocalMaximums();
}

//--------------------
  /**
   * @brief Makes map of local minimums of requested rank (IN ROWS AND COLUMNS ONLY!)
   * 
   * Map of local minimums is a 2-d array of unsigned integer values of data shape, 
   * it has 0/1/2 bits for non-maximum / maximum in column / maximum in row of radius rank.   
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  rank - radius of the region in which central pixel has a maximal value
   */

template <typename T>
ndarray<const AlgImgProc::pixel_minimums_t, 2> 
mapOfLocalMinimums( const ndarray<const T,2> data
                  , const ndarray<const AlgImgProc::mask_t,2> mask
                  , const size_t& rank
                  )
{
  AlgImgProc algo(0); // , 0, 1e6, 0, 1e6, 1023);
  algo.validate_window(data.shape());
  algo._makeMapOfLocalMinimums<T>(data, mask, rank);
  return algo.mapOfLocalMinimums();
}

//--------------------
//--------------------

  /**
   * @brief Wrapper for _makeMapOfLocalMaximumsRank1Cross
   *
   * Returns (uint16) array with local maximum bit info 1/2 - in row/column (=3 - maximum in cross region) 
   * @param[in]  fphoton - (float) ndarray of fractional number of photons per pixel
   */
ndarray<const AlgImgProc::pixel_maximums_t, 2>
mapOfLocalMaximumsRank1Cross(const ndarray<const AlgImgProc::fphoton_t,2> fphoton);

//--------------------
//--------------------

} // namespace ImgAlgos

#endif // IMGALGOS_ALGIMGPROC_H
