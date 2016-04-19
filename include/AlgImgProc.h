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
 *   // The same peak-finders after revision-1
 *   std::vector<Peak>& peaks = aip -> peakFinderV2r1<T>(seg_data, seg_mask, thr, r0, dr);
 *   std::vector<Peak>& peaks = aip -> peakFinderV3r1<T>(seg_data, seg_mask, rank, r0, dr);
 *   std::vector<Peak>& peaks = aip -> peakFinderV4r1<T>(seg_data, seg_mask, thr_low, thr_high, rank, r0, dr);
 *
 *   std::vector<Peak>& peaks = aip -> getVectorOfSelectedPeaks();
 *   std::vector<Peak>& peaks = aip -> getVectorOfPeaks();
 *
 *   // call after peakFinderV2 ONLY!
 *   ndarray<conmap_t, 2>& conmap = aip -> mapOfConnectedPixels();
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
  typedef uint32_t conmap_t;
  typedef uint16_t pixel_status_t;
  typedef uint16_t pixel_maximums_t;
  typedef float son_t;

  /**
   * @brief Class constructor is used for initialization of all paramaters. 
   * 
   * @param[in] seg    - ROI segment index in the ndarray
   * @param[in] rowmin - ROI window limit
   * @param[in] rowmax - ROI window limit
   * @param[in] colmin - ROI window limit
   * @param[in] colmax - ROI window limit
   * @param[in] pbits  - print control bit-word; =0-print nothing, +1-input parameters, +2-algorithm area, +128-all details.
   */

  AlgImgProc( const size_t&   seg
	    , const size_t&   rowmin  = 0
	    , const size_t&   rowmax  = 1e6
	    , const size_t&   colmin  = 0
	    , const size_t&   colmax  = 1e6
	    , const unsigned& pbits   = 0
	    , const unsigned& npksmax = 10000
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
                            const float& son_min =0);

  /// Returns reference to Window object
  const Window& window(){ return m_win; }

  /// Returns segment index in the ndarray
  const size_t& segind(){ return m_seg; }

  /// Returns vector of all found peaks for this segment/window
  std::vector<Peak>& getVectorOfPeaks(){ return v_peaks; }

  /// Returns vector of selected peaks for this segment/window
  std::vector<Peak>& getVectorOfSelectedPeaks(){ return v_peaks_sel; }

  /// Returns map of connected pixels after peakFinderV2(.)
  ndarray<conmap_t, 2>& mapOfConnectedPixels() { return m_conmap; }

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
  bool     m_use_mask;

  conmap_t m_numreg;

  float    m_r0;       // radial parameter of the ring for S/N evaluation algorithm
  float    m_dr;       // ring width for S/N evaluation algorithm 
  size_t   m_rank;     // rank of maximum for peakFinderV3

  SoNResult m_sonres_def;

  float    m_peak_npix_min; // peak selection parameter
  float    m_peak_npix_max; // peak selection parameter
  float    m_peak_amax_thr; // peak selection parameter
  float    m_peak_atot_thr; // peak selection parameter
  float    m_peak_son_min;  // peak selection parameter
  bool     m_do_preselect;  // flag for applying peak pre-selection before S/N calculation

  //Peak     m_peak;

  ndarray<pixel_status_t, 2>   m_pixel_status;
  ndarray<pixel_maximums_t, 2> m_local_maximums;
  ndarray<conmap_t, 2>         m_conmap;
  std::vector<PeakWork>        v_peaks_work;
  std::vector<Peak>            v_peaks;
  std::vector<Peak>            v_peaks_sel;
  std::vector<TwoIndexes>      v_indexes;
  std::vector<TwoIndexes>      v_inddiag;

  //ndarray<Peak, 1>           m_peaks;
  // mask_t*                    m_seg_mask_def;


  /// Returns string name of the class for messanger
  inline const char* _name() {return "ImgAlgos::AlgImgProc";}

  /// Recursive method which checks m_pixel_status[r][c] and numerates connected regions in m_conmap[r][c].
  void _findConnectedPixels(const unsigned& r, const unsigned& c);

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

//--------------------
  /**
   * @brief Makes m_pixel_status array by setting to 1/0 all good-above-threshold/bad pixels
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr  - threshold on data values
   */

template <typename T>
void
_makeMapOfPixelStatus( const ndarray<const T,2>&      data
                     , const ndarray<const mask_t,2>& mask
                     , const T& thr
                     )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfPixelStatus, seg=" << m_seg << " thr=" << thr << "\n    in window: " << m_win);
  //if(m_pbits & 512) m_win.print();

  //if(m_pixel_status.size()==0) 
  if(m_pixel_status.empty()) 
     m_pixel_status = make_ndarray<pixel_status_t>(data.shape()[0], data.shape()[1]);

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      m_pixel_status[r][c] = (mask[r][c] && (data[r][c]>thr)) ? 255 : 0;
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

  // Next line does not work when vector is large.
  //if(v_peaks_work.capacity() < m_numreg+1) v_peaks_work.reserve(m_numreg+1); // use numeration from 1
  v_peaks_work.clear();

  //std::fill_n(&v_peaks_work[0], int(m_numreg+1), pw0);
  for(unsigned reg = 0; reg < min(m_numreg+1, m_npksmax); reg++)
    v_peaks_work.push_back(pw0);

  //std::cout << "\nXXXIII:Point 3";
  //std::cout << "  m_numreg: " << m_numreg;
  //std::cout << "  v_peaks_work.size(): " << v_peaks_work.size();
  //std::cout << "  v_peaks_work.capacity(): " << v_peaks_work.capacity();
  //std::cout << '\n';
  //std::cout << "  m_win: " << m_win;
  //std::cout << "\n  m_pixel_status.size(): " << m_pixel_status.size();
  //std::cout << "  m_conmap.size(): " << m_conmap.size();
  //std::cout << '\n';

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      conmap_t ireg = m_conmap[r][c];

      if(ireg<1) continue;

      if(!(ireg<m_npksmax)) {
         MsgLog(_name(), warning, "Number of peak regions: " << ireg 
		<< " reached the max reserved number of peaks: " << m_npksmax
                << " in segment: " << m_seg);
         break;
      }

      //std::cout << " reg=" << ireg;

      double amp = (double)data[r][c];

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
_makeMapOfLocalMaximums( const ndarray<const T,2>&      data
                       , const ndarray<const mask_t,2>& mask
                       , const size_t& rank
                       )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfLocalMaximums, seg=" << m_seg << " rank=" << rank << "\n    in window: " << m_win);
  //if(m_pbits & 512) m_win.print();

  // initialization of indexes
  if(v_inddiag.empty()) _evaluateDiagIndexes(rank);

  if(m_local_maximums.empty()) 
     m_local_maximums = make_ndarray<pixel_maximums_t>(data.shape()[0], data.shape()[1]);
  std::fill_n(&m_local_maximums[0][0], int(data.size()), pixel_maximums_t(0));

  unsigned rmin = max((int)m_win.rowmin, int(0+rank));
  unsigned rmax = min((int)m_win.rowmax, int(data.shape()[0]-rank));
  unsigned cmin = max((int)m_win.colmin, int(0+rank));
  unsigned cmax = min((int)m_win.colmax, int(data.shape()[1]-rank));


  // check rank maximum in columns and set the 1st bit (1)
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      if(!mask[r][c]) continue;
      m_local_maximums[r][c] = 1;
      int cp=c; 
      int cm=c;
      for(unsigned d=0; d<rank; ++d) {
        cp++; cm--;
	if((mask[r][cp] && data[r][cp] > data[r][c]) 
	|| (mask[r][cm] && data[r][cm] > data[r][c])) {
          m_local_maximums[r][c] &=~1; // clear 1st bit
	  break;
	}
      }
      // (r,c) is a local peak, jump ahead through the tested rank range
      if(m_local_maximums[r][c] & 1) c+=rank;
    }
  }

  // check rank maximum in rows and set the 2nd bit (2)
  for(unsigned c = cmin; c<cmax; c++) {
    for(unsigned r = rmin; r<rmax; r++) {
      // if it is not a local maximum from previous algorithm
      if(!m_local_maximums[r][c]) continue;
      m_local_maximums[r][c] |= 2; // set 2nd bit
      int rp=r; 
      int rm=r;
      for(unsigned d=0; d<rank; ++d) { 
        rp++; rm--;
	if((mask[rp][c] && data[rp][c] > data[r][c])
	|| (mask[rm][c] && data[rm][c] > data[r][c])) {
          m_local_maximums[r][c] &=~2; // clear 2nd bit
	  break;
	}
      }
      // if (r,c) is a local peak, jump ahead through the tested rank range
      if(m_local_maximums[r][c] & 2) r+=rank;
    }
  }

  // check rank maximum in "diagonal" regions and set the 3rd bit (4)
  for(unsigned r = rmin; r<rmax; r++) {
    for(unsigned c = cmin; c<cmax; c++) {
      // if it is not a local maximum from two previous algorithm
      if(!(m_local_maximums[r][c] & 2)) continue;
      m_local_maximums[r][c] |= 4; // set 3rd bit

      for(vector<TwoIndexes>::const_iterator ij  = v_inddiag.begin();
                                             ij != v_inddiag.end(); ij++) {
        int ir = r + (ij->i);
        int ic = c + (ij->j);

	if(mask[ir][ic] && data[ir][ic] > data[r][c]) {
          m_local_maximums[r][c] &=~4; // clear 3rd bit
	  break;
	}
      }
      // (r,c) is a local peak, jump ahead through the tested rank range
      if(m_local_maximums[r][c] & 4) c+=rank;
    }
  }
}

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
_procLocalMaximum( const ndarray<const T,2>&      data
                 , const ndarray<const mask_t,2>& mask
                 , const size_t& rank
                 , const unsigned& r0
                 , const unsigned& c0
            )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _procLocalMaximum, seg=" << m_seg << " r0=" << r0 << " c0=" << c0 << " rank=" << rank);

  double   a0 = data[r0][c0];
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

      if(!mask[r][c]) continue;
      double a = data[r][c];
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

  double   a0 = data[r0][c0];
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

      if(!mask[r][c]) continue;
      double a = data[r][c];
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
  peak.amp_max   = a0;
  peak.amp_tot   = samp;
  peak.row_cgrav = sar1/swei;
  peak.col_cgrav = sac1/swei;
  peak.row_sigma = (npos>1) ? std::sqrt(sar2/swei - peak.row_cgrav * peak.row_cgrav) : 0;
  peak.col_sigma = (npos>1) ? std::sqrt(sac2/swei - peak.col_cgrav * peak.col_cgrav) : 0;
  peak.row_min   = rmin;
  peak.row_max   = rmax;
  peak.col_min   = cmin;
  peak.col_max   = cmax;  
  peak.bkgd      = 0; //sonres.avg;
  peak.noise     = 0; //sonres.rms;
  peak.son       = 0; //sonres.son;

  /*
  cout << "peak.row_sigma=" << peak.row_sigma 
       << " col_sigma=" << peak.col_sigma 
       << " samp=" << samp 
       << '\n';
  */

  if(v_peaks.size()<m_npksmax-1) v_peaks.push_back(peak);
  //if(_peakIsPreSelected(peak) && v_peaks.size()<m_npksmax-1) v_peaks.push_back(peak);
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
_makePeaksFromMapOfLocalMaximums( const ndarray<const T,2>&      data
                                , const ndarray<const mask_t,2>& mask
                                , const size_t& rank
                                )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makePeaksFromMapOfLocalMaximums, seg=" << m_seg<< ", rank=" << rank);

  v_peaks.clear(); // this vector will be filled out for each window

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      if(m_local_maximums[r][c] & 4)
        _procLocalMaximum<T>(data,mask,rank,r,c);	
}

//--------------------
  /**
   * @brief The same as _makePeaksFromMapOfLocalMaximums, but uses _procLocalMaximumV2<T>(data,mask,rank,r,c);
   */
template <typename T>
void 
_makePeaksFromMapOfLocalMaximumsV2( const ndarray<const T,2>&      data
                                , const ndarray<const mask_t,2>& mask
                                , const size_t& rank
                                )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makePeaksFromMapOfLocalMaximumsV2, seg=" << m_seg<< ", rank=" << rank);

  v_peaks.clear(); // this vector will be filled out for each window

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      if(m_local_maximums[r][c] & 4)
        _procLocalMaximumV2<T>(data,mask,rank,r,c);	
}

//--------------------
//--------------------
//--------------------
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
void _addSoNToPeaks( const ndarray<const T,2>& data
                   , const ndarray<const mask_t,2>& mask
	           , const float r0 = 5
	           , const float dr = 0.05
                   )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _addSoNToPeaks, seg=" << m_seg);

  setSoNPars(r0, dr);

  std::vector<Peak>::iterator it;

  for(it=v_peaks.begin(); it!=v_peaks.end(); ++it) { 
    Peak& peak = (*it);
    
    SoNResult sonres = evaluateSoNForPixel<T>((unsigned) peak.row, (unsigned) peak.col, data, mask);

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
	             , const float r0 = 5
	             , const float dr = 0.05
                     )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _addSoNToPeaks, seg=" << m_seg);

  setSoNPars(r0, dr);

  std::vector<Peak>::iterator it;

  for(it=v_peaks.begin(); it!=v_peaks.end(); ++it) { 
    Peak& peak = (*it);
    
    SoNResult sonres = evaluateSoNForPixel<T>((unsigned) peak.row, (unsigned) peak.col, data, mask);

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
_procDroplet( const ndarray<const T,2>&      data
            , const ndarray<const mask_t,2>& mask
            , const T& thr_low
            , const unsigned& rad
            , const unsigned& r0
            , const unsigned& c0
            )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _procDroplet, seg=" << m_seg << " r0=" << r0 << " c0=" << c0);

  double a0 = data[r0][c0];
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
      double a = data[r][c];

      if(a>a0) return;  // This is not a local maximum inside rad...

      if(mask[r][c] && a>thr_low) {
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
   * @brief _makeVectorOfDroplets - a part of peakFinderV1 algorithm - loops over pixels and makes peak candidates in v_peaks
   * 
   * @param[in]  data - ndarray with calibrated intensities
   * @param[in]  mask - ndarray with mask of bad/good (0/1) pixels
   * @param[in]  thr_low   - threshold on pixel intensity to be considered in this algorithm 
   * @param[in]  thr_high  - threshold on pixel intensity to be a candidate to "droplet"
   * @param[in]  rad       - radius in pixels of squared region to find droplet relative to central pixel
   */
template <typename T>
void 
_makeVectorOfDroplets( const ndarray<const T,2>&      data
                     , const ndarray<const mask_t,2>& mask
                     , const T& thr_low
                     , const T& thr_high
                     , const unsigned& rad=5
                     )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV1, seg=" << m_seg);

  v_peaks.clear(); // this vector will be filled out for each window

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++)
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++)
      if(mask[r][c] && (data[r][c]>thr_high)) 
        _procDroplet<T>(data,mask,thr_low,rad,r,c);	
}

//--------------------
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
numberOfPixAboveThr( const ndarray<const T,2>&      data
                   , const ndarray<const mask_t,2>& mask
                   , const T& thr
                   )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in numberOfPixAboveThr, seg=" << m_seg);

  m_win.validate(data.shape());

  unsigned npix = 0;
  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      if(mask[r][c] && (data[r][c]>thr)) npix++;
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
intensityOfPixAboveThr( const ndarray<const T,2>&      data
                      , const ndarray<const mask_t,2>& mask
                      , const T& thr
                      )
{
  if(m_pbits & 512) MsgLog(_name(), info, "in intensityOfPixAboveThr, seg=" << m_seg);

  m_win.validate(data.shape());

  double amptot = 0;
  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {
      if(mask[r][c] && (data[r][c]>thr)) amptot += (double)data[r][c];
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
peakFinderV1( const ndarray<const T,2>&      data
            , const ndarray<const mask_t,2>& mask
            , const T& thr_low
            , const T& thr_high
            , const unsigned& rad=5
            , const float&    dr=0.05
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
peakFinderV4( const ndarray<const T,2>&      data
            , const ndarray<const mask_t,2>& mask
            , const T& thr_low
            , const T& thr_high
            , const unsigned& rank = 5
            , const float&    r0   = 5
            , const float&    dr   = 0.05
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
peakFinderV4r1( const ndarray<const T,2>&      data
              , const ndarray<const mask_t,2>& mask
              , const T& thr_low
              , const T& thr_high
              , const unsigned& rank = 5
              , const float&    r0   = 5
              , const float&    dr   = 0.05
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
	    , const float r0 = 5
	    , const float dr = 0.05
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
	      , const float r0 = 5
	      , const float dr = 0.05
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
	    , const float r0 = 5
	    , const float dr = 0.05
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
   * @brief peakFinderV3r1 - "Ranker" - the same as V3, but returns background-corrected amp_max, amp_tot, son (total).
   * Changes: 
   *   m_do_preselect = false; 
   *   _addSoNToPeaksV2<T>(data, mask, r0, dr);
   *   _makePeaksFromMapOfLocalMaximumsV2<T>(data, mask, rank); // counts all pixels inside rank region (not only a>0)
   */

template <typename T>
std::vector<Peak>&
peakFinderV3r1( const ndarray<const T,2>&      data
              , const ndarray<const mask_t,2>& mask
              , const size_t rank = 4
	      , const float r0 = 5
	      , const float dr = 0.05
              )
{
  m_win.validate(data.shape());

  if(m_pbits & 512) MsgLog(_name(), info, "in peakFinderV3r1, seg=" << m_seg << " win" << m_win << " rank=" << rank);

  _makeMapOfLocalMaximums<T>(data, mask, rank);

  m_do_preselect = false;
  _makePeaksFromMapOfLocalMaximumsV2<T>(data, mask, rank);
  _addSoNToPeaksV2<T>(data, mask, r0, dr);
  _makeVectorOfSelectedPeaks();

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
evaluateSoNForPixel( const unsigned& row
                   , const unsigned& col
                   , const ndarray<const T,2>& data
                   , const ndarray<const mask_t,2>& mask
                   )
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in evaluateSoNForPixel, seg=" << m_seg << " row=" << row << ", col=" << col);

  // S/N algorithm initialization
  if(! m_init_son_is_done) {
    _evaluateRingIndexes(m_r0, m_dr);
    m_win.validate(data.shape());
    m_use_mask = (mask.empty()) ? false : true;
    m_init_son_is_done = true;
  }

  if(m_use_mask && (!mask[row][col])) return m_sonres_def;
  //if(m_use_mask && (!mask[row][col])) return SoNResult({});

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
    if(m_use_mask && (! mask[ir][ic])) continue;

    amp = (double)data[ir][ic];
    sum0 ++;
    sum1 += amp;
    sum2 += amp*amp;
  }
  //SoNResult res = m_sonres_def;
  SoNResult res;

  if(sum0) {
    res.avg = sum1/sum0;                              // Averaged background level
    res.rms = std::sqrt(sum2/sum0 - res.avg*res.avg); // RMS os the background around peak
    res.sig = data[row][col]      - res.avg;          // Signal above the background
    if (res.rms>0) res.son = res.sig/res.rms;         // S/N ratio
  }

  return res;
}

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
      result[r][c] = evaluateSoNForPixel<T>(r, c, data, mask);
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
      son[r][c] = evaluateSoNForPixel<T>(r, c, data, mask).son;
    }
  }
}

//--------------------
//--------------------
//--------------------

};

} // namespace ImgAlgos

#endif // IMGALGOS_ALGIMGPROC_H
