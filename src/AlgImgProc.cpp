//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class AlgImgProc
//
// Author List:
//      Mikhail S. Dubrovin
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ImgAlgos/AlgImgProc.h"

//-----------------
// C/C++ Headers --
//-----------------
#include <cmath>     // floor, ceil
#include <iomanip>   // for std::setw
#include <sstream>   // for stringstream

namespace ImgAlgos {

//----------------
// Constructors --
//----------------

AlgImgProc::AlgImgProc ( const size_t&   seg
                       , const size_t&   rowmin
                       , const size_t&   rowmax
                       , const size_t&   colmin
                       , const size_t&   colmax
		       , const unsigned& pbits
		       , const unsigned& npksmax
                       )
  : m_pbits(pbits)
  , m_npksmax(npksmax)
  , m_seg(seg)
  , m_init_son_is_done(false)
  , m_init_bkgd_is_done(false)
  , m_r0(7.0)
  , m_dr(2.0)
  , m_sonres_def()
  , m_peak_npix_min(0)
  , m_peak_npix_max(1e6)
  , m_peak_amax_thr(0)
  , m_peak_atot_thr(0)
  , m_peak_son_min(0)
  , m_do_preselect(true)
{
  if(m_pbits & 512) MsgLog(_name(), info, "in c-tor AlgImgProc, seg=" << m_seg);
  m_win.set(seg, rowmin, rowmax, colmin, colmax);

  if(m_pbits & 2) printInputPars();
}

//--------------------

AlgImgProc::~AlgImgProc() 
{
  if(m_pbits & 512) MsgLog(_name(), info, "in d-tor ~AlgImgProc, seg=" << m_seg);
  //v_peaks.resize(0); 
  //v_peaks_work.resize(0);
}

//--------------------

void 
AlgImgProc::_reserveVectorOfPeaks()
{
  if(v_peaks_work.capacity() == m_npksmax) return;

  v_peaks_work.reserve(m_npksmax);
  v_peaks.reserve(m_npksmax);
  v_peaks_sel.reserve(m_npksmax);
}

//--------------------

void 
AlgImgProc::printInputPars()
{
  std::stringstream ss; 
  ss << "printInputPars:\n"
     << "\npbits   : " << m_pbits
     << "\nnpksmax : " << m_npksmax
     << "\nwindow  : " << m_win
     << '\n';
  //ss << "\nrmin    : " << m_r0
  //   << "\ndr      : " << m_dr
  MsgLog(_name(), info, ss.str()); 

  if(m_pbits & 512) printMatrixOfRingIndexes();
  if(m_pbits & 512) printVectorOfRingIndexes();
}

//--------------------

void 
AlgImgProc::_makeMapOfConnectedPixels()
{
  //if(m_conmap.size()==0) 
  if(m_conmap.empty()) 
     m_conmap = make_ndarray<conmap_t>(m_pixel_status.shape()[0], m_pixel_status.shape()[1]);

  std::fill_n(m_conmap.data(), int(m_pixel_status.size()), conmap_t(0));
  m_numreg=0;

  for(int r = (int)m_win.rowmin; r<(int)m_win.rowmax; r++) {
    for(int c = (int)m_win.colmin; c<(int)m_win.colmax; c++) {

      if(!(m_pixel_status(r,c) & 1)) continue;
      ++ m_numreg;
      //if(m_numreg == m_npksmax) break;
      _findConnectedPixels(r, c);
    }
  }

  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfConnectedPixels, seg=" << m_seg << ", m_numreg=" <<  m_numreg);
}

//--------------------

void 
AlgImgProc::_findConnectedPixels(const int& r, const int& c)
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _findConnectedPixels, seg=" << m_seg);

  if(! (m_pixel_status(r,c) & 1)) return;

  m_pixel_status(r,c) ^= 1; // set the 1st bit to zero.
  m_conmap(r,c) = m_numreg;

  // cout  << "ZZZ:  r:" << r << " c:" << c << '\n';

  if(  r+1 < (int)m_win.rowmax ) _findConnectedPixels(r+1, c);
  if(  c+1 < (int)m_win.colmax ) _findConnectedPixels(r, c+1);
  if(!(r-1 < (int)m_win.rowmin)) _findConnectedPixels(r-1, c);
  if(!(c-1 < (int)m_win.colmin)) _findConnectedPixels(r, c-1);  
}

//--------------------
// Templated recursive method with data for fast termination
// Returns true / false for pixel is ok / m_reg_a0 is not a local maximum
template <typename T>
bool
AlgImgProc::_findConnectedPixelsInRegionV1(const ndarray<const T,2>& data, const int& r, const int& c)
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _findConnectedPixelsInRegionV1, seg=" << m_seg);

  pixel_status_t pstat = m_pixel_status(r,c);
  if(! pstat) return true;    // pstat=0 - masked
  if(pstat & 28) return true; // pstat=4/8/16 : a<thr_low/used in recursion/used in map

  if(data(r,c) > m_reg_a0) return false; // initial point is not a local maximum - TERMINATE

  m_pixel_status(r,c) |= 8; // mark this pixel as used in this recursion to get rid of cycling
  v_ind_pixgrp.push_back(TwoIndexes(r,c));

  if(  r+1 < m_reg_rmax)  if(! _findConnectedPixelsInRegionV1<T>(data, r+1, c)) return false;
  if(  c+1 < m_reg_cmax)  if(! _findConnectedPixelsInRegionV1<T>(data, r, c+1)) return false;
  if(!(r-1 < m_reg_rmin)) if(! _findConnectedPixelsInRegionV1<T>(data, r-1, c)) return false;
  if(!(c-1 < m_reg_cmin)) if(! _findConnectedPixelsInRegionV1<T>(data, r, c-1)) return false;  

  return true;
}

//--------------------
// Templated recursive method with data for fast termination
// Returns true / false for pixel is ok / m_reg_a0 is not a local maximum
// V2 - find a group connected pixels for already found local maximum
//     - so do not need to terminate recursion like in V1 

template <typename T>
void
AlgImgProc::_findConnectedPixelsInRegionV2(const ndarray<const T,2>& data,
                                           const ndarray<const mask_t,2>& mask, 
                                           const int& r, const int& c)
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _findConnectedPixelsInRegionV2, r=" << r << " c=" << c);
  //     << " data=" << data(r,c) << " m_reg_thr=" << m_reg_thr <<);

  if(! mask(r,c)) return;   // - masked
  if(m_conmap(r,c)) return; // - pixel is already used
  if(data(r,c) < m_reg_thr) return; // discard pixel below threshold if m_reg_thr != 0 

  //if(m_pixel_status(r,c) & 8) return; // pstat=4/8/16 : a<thr_low/used in recursion/used in map
  //m_pixel_status(r,c) |= 8; // mark this pixel as used in this recursion to get rid of cycling  

  m_conmap(r,c) = m_numreg; // mark pixel on map

  v_ind_pixgrp.push_back(TwoIndexes(r,c));

  if(  r+1 < m_reg_rmax)  _findConnectedPixelsInRegionV2<T>(data, mask, r+1, c);
  if(  c+1 < m_reg_cmax)  _findConnectedPixelsInRegionV2<T>(data, mask, r, c+1);
  if(!(r-1 < m_reg_rmin)) _findConnectedPixelsInRegionV2<T>(data, mask, r-1, c);
  if(!(c-1 < m_reg_cmin)) _findConnectedPixelsInRegionV2<T>(data, mask, r, c-1);  
}

//--------------------
// NON-Templated recursive method doing full flood filling in the region.
void 
AlgImgProc::_findConnectedPixelsInRegion(const int& r, const int& c)
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _findConnectedPixelsInRegion, seg=" << m_seg);

  pixel_status_t pstat = m_pixel_status(r,c);
  if(! pstat) return;    // pstat=0 - masked
  if(pstat & 28) return; // pstat=4/8/16 : a<thr_low/used in recursion/used in map
  //if(m_conmap(r,c)) return; // pixel belongs to any other group - it is also marked as used

  //cout  << "ZZZ:  v_ind_pixgrp.push_back r:" << r << " c:" << c << '\n';

  m_pixel_status(r,c) |= 8; // mark this pixel as used in this recursion 
  v_ind_pixgrp.push_back(TwoIndexes(r,c));

  if(  r+1 < m_reg_rmax)  _findConnectedPixelsInRegion(r+1, c);
  if(  c+1 < m_reg_cmax)  _findConnectedPixelsInRegion(r, c+1);
  if(!(r-1 < m_reg_rmin)) _findConnectedPixelsInRegion(r-1, c);
  if(!(c-1 < m_reg_cmin)) _findConnectedPixelsInRegion(r, c-1);  
}
//--------------------

void 
AlgImgProc::_addPixGroupToMap()
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _addPixGroupToMap, seg=" << m_seg << " numreg=" << m_numreg);
  for(vector<TwoIndexes>::const_iterator ij  = v_ind_pixgrp.begin();
                                         ij != v_ind_pixgrp.end(); ij++) {
    m_conmap(ij->i, ij->j) = m_numreg;
    m_pixel_status(ij->i, ij->j) |= 16; // mark pixel as used on map (set bit)  
  }
}

//--------------------

void 
AlgImgProc::_clearStatusOfUnusedPixels() {
  for(vector<TwoIndexes>::const_iterator ij  = v_ind_pixgrp.begin();
                                         ij != v_ind_pixgrp.end(); ij++) {
    m_pixel_status(ij->i, ij->j) &=~8; // mark pixel as un-used (clear bit)
  }
}

//--------------------

bool
AlgImgProc::_peakWorkIsPreSelected(const PeakWork& pw)
{
  if (! m_do_preselect) return true;
  if (pw.peak_npix < m_peak_npix_min) return false;
  if (pw.peak_npix > m_peak_npix_max) return false;
  if (pw.peak_amax < m_peak_amax_thr) return false;
  if (pw.peak_atot < m_peak_atot_thr) return false;
  return true;
}

//--------------------

bool
AlgImgProc::_peakIsPreSelected(const Peak& peak)
{
  if (! m_do_preselect) return true;
  if (peak.npix    < m_peak_npix_min) return false;
  if (peak.npix    > m_peak_npix_max) return false;
  if (peak.amp_max < m_peak_amax_thr) return false;
  if (peak.amp_tot < m_peak_atot_thr) return false;
  return true;
}

//--------------------

bool
AlgImgProc::_peakIsSelected(const Peak& peak)
{
  if (peak.son     < m_peak_son_min)  return false;
  if (peak.npix    < m_peak_npix_min) return false;
  if (peak.npix    > m_peak_npix_max) return false;
  if (peak.amp_max < m_peak_amax_thr) return false;
  if (peak.amp_tot < m_peak_atot_thr) return false;
  return true;
}

//--------------------

void
AlgImgProc::_makeVectorOfPeaks()
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeVectorOfPeaks, seg=" << m_seg << " m_numreg=" << m_numreg);
  //m_peaks = make_ndarray<Peak>(m_numreg);

  if(m_numreg==0) return;

  //v_peaks.reserve(m_numreg+1); // this does not always work
  _reserveVectorOfPeaks();
  v_peaks.clear();

  for(unsigned i=0; i<min(m_numreg, m_npksmax); i++) {

    PeakWork& pw = v_peaks_work[i+1]; // region number begins from 1

    if(! _peakWorkIsPreSelected(pw)) continue;

    Peak   peak; // = v_peaks[i];

    peak.seg       = m_seg;
    peak.npix      = pw.peak_npix;
    peak.row       = pw.peak_row;
    peak.col       = pw.peak_col;
    peak.amp_max   = pw.peak_amax;
    peak.amp_tot   = pw.peak_atot;

    if (pw.peak_atot>0) {
      peak.row_cgrav = pw.peak_ar1/pw.peak_atot;
      peak.col_cgrav = pw.peak_ac1/pw.peak_atot;
      peak.row_sigma = std::sqrt(pw.peak_ar2/pw.peak_atot - peak.row_cgrav * peak.row_cgrav);
      peak.col_sigma = std::sqrt(pw.peak_ac2/pw.peak_atot - peak.col_cgrav * peak.col_cgrav);
    }
    else {
      peak.row_cgrav = pw.peak_row;
      peak.col_cgrav = pw.peak_col;
      peak.row_sigma = 0;
      peak.col_sigma = 0;
    }

    peak.row_min   = pw.peak_rmin;
    peak.row_max   = pw.peak_rmax;
    peak.col_min   = pw.peak_cmin;
    peak.col_max   = pw.peak_cmax;
    peak.bkgd  = 0;
    peak.noise = 0;
    peak.son   = 0;

    v_peaks.push_back(peak);
  }
}
//--------------------

void
AlgImgProc::_makeVectorOfSelectedPeaks()
{
  if(v_peaks_sel.capacity() != m_npksmax) v_peaks_sel.reserve(m_npksmax);
     v_peaks_sel.clear();

  //std::vector<Peak>::iterator it;
  for(std::vector<Peak>::iterator it=v_peaks.begin(); it!=v_peaks.end(); ++it) { 
    Peak& peak = (*it);
    if(_peakIsSelected(peak)) v_peaks_sel.push_back(peak);
  }
  if(m_pbits & 512) MsgLog(_name(), info, "in _makeVectorOfSelectedPeaks, seg=" << m_seg 
                           << "  #peaks raw=" << v_peaks.size() 
                           << "  sel=" << v_peaks_sel.size());
}

//--------------------

void 
AlgImgProc::_evaluateDiagIndexes(const size_t& rank)
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _evaluateDiagIndexes, seg=" << m_seg << " rank=" << rank);

  m_rank = rank;
  int indmax =  m_rank;
  int indmin = -m_rank;
  unsigned npixmax = (2*rank+1)*(2*rank+1);
  if(v_inddiag.capacity() < npixmax) v_inddiag.reserve(npixmax);
  v_inddiag.clear();

  for (int i = indmin; i <= indmax; ++ i) {
    for (int j = indmin; j <= indmax; ++ j) {

      // use rectangular region of radius = rank
      // remove already tested central row and column
      if (i==0 || j==0) continue;
      // use ring region (if un-commented)
      //if (m_rank>2 && floor(std::sqrt(float(i*i + j*j)))>(int)m_rank) continue;
      //TwoIndexes inds = {i,j};
      TwoIndexes inds(i,j);
      v_inddiag.push_back(inds);
    }
  }

  if(m_pbits & 2) printMatrixOfDiagIndexes();
}

//--------------------

void 
AlgImgProc::_evaluateRingIndexes(const float& r0, const float& dr)
{
  if(m_pbits & 512) MsgLog(_name(), info, "in _evaluateRingIndexes, seg=" << m_seg << " r0=" << r0 << " dr=" << dr);

  m_r0 = r0;
  m_dr = dr;

  int indmax = (int)std::ceil(m_r0 + m_dr);
  int indmin = -indmax;
  unsigned npixmax = (2*indmax+1)*(2*indmax+1);
  if(v_indexes.capacity() < npixmax) v_indexes.reserve(npixmax);
  v_indexes.clear();

  for (int i = indmin; i <= indmax; ++ i) {
    for (int j = indmin; j <= indmax; ++ j) {

      float r = std::sqrt( float(i*i + j*j) );
      if ( r < m_r0 || r > m_r0 + m_dr ) continue;
      //TwoIndexes inds = {i,j};
      TwoIndexes inds(i,j);
      v_indexes.push_back(inds);
    }
  }

  if(m_pbits & 2) printMatrixOfRingIndexes();
  if(m_pbits & 4) printVectorOfRingIndexes();
}

//--------------------

void 
AlgImgProc::_fillCrossIndexes()
{
  v_indcross.push_back(TwoIndexes(-1, 0));
  v_indcross.push_back(TwoIndexes( 1, 0));
  v_indcross.push_back(TwoIndexes( 0,-1));
  v_indcross.push_back(TwoIndexes( 0, 1));
}

//--------------------

void 
AlgImgProc::setSoNPars(const float& r0, const float& dr)
{ 
  if(m_pbits & 512) MsgLog(_name(), info, "in setSoNPars, seg=" << m_seg << " r0=" << r0 << " dr=" << dr);

  if(r0==m_r0 && dr==m_dr) return;
  _evaluateRingIndexes(r0, dr);
}

//--------------------

void
AlgImgProc::setPeakSelectionPars(const float& npix_min, const float& npix_max,
                                 const float& amax_thr, const float& atot_thr, const float& son_min)
{
  if(m_pbits & 512) MsgLog(_name(), info, "in setPeakSelectionPars, seg=" << m_seg);
  m_peak_npix_min = npix_min;
  m_peak_npix_max = npix_max;
  m_peak_amax_thr = amax_thr;
  m_peak_atot_thr = atot_thr;
  m_peak_son_min  = son_min;
}

//--------------------

void
AlgImgProc::_mergeConnectedPixelCouples(const fphoton_t& thr_on_max, const fphoton_t& thr_on_tot, const bool& DO_TEST)
{
  // Available ndarrays:
  // m_nphoton - uint number of photons
  // m_fphoton - fractional number of photons [0,1)
  // m_mphoton - for test output only
  // m_local_maximums - map of local maximums
  // m_pixel_status - is not used
  // m_conmap - map of connected pixels (coupled pixels marked by the same group number) - busy pixels

  const unsigned* shape = m_fphoton.shape();

  if(m_pbits & 512) MsgLog(_name(), info, "in _mergeConnectedPixelCouples, seg=" << m_seg << "\n    in window: " << m_win);

  if(DO_TEST) {
    if(m_mphoton.empty()) 
       m_mphoton = make_ndarray<nphoton_t>(shape[0], shape[1]);
    std::fill_n(&m_mphoton(0,0), int(m_nphoton.size()), nphoton_t(0));
  }

  if(m_conmap.empty()) 
     m_conmap = make_ndarray<conmap_t>(shape[0], shape[1]);
  std::fill_n(m_conmap.data(), int(m_fphoton.size()), conmap_t(0));

  m_numreg=0;

  //if(v_indcross.empty()) _fillCrossIndexes();

  TwoIndexes pU(-1, 0);
  TwoIndexes pD( 1, 0);
  TwoIndexes pL( 0,-1);
  TwoIndexes pR( 0, 1);

  // CENTARL AREA
  // merge central pixels
  TwoIndexes aC[]={pU, pD, pL, pR}; vector<TwoIndexes> vC(aC,aC+4);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmin+1, m_win.rowmax-1, m_win.colmin+1, m_win.colmax-1, vC);

  // EDGES
  // merge left-edge pixels
  TwoIndexes aL[]={pU, pD, pR}; vector<TwoIndexes> vL(aL,aL+3);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmin+1, m_win.rowmax-1, m_win.colmin, m_win.colmin+1, vL);

  // merge right-edge pixels
  TwoIndexes aR[]={pU, pD, pL}; vector<TwoIndexes> vR(aR,aR+3);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmin+1, m_win.rowmax-1, m_win.colmax-1, m_win.colmax, vR);

  // merge upper-edge pixels
  TwoIndexes aU[]={pD, pL, pR}; vector<TwoIndexes> vU(aU,aU+3);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmin, m_win.rowmin+1, m_win.colmin+1, m_win.colmax-1, vU);

  // merge down-edge pixels
  TwoIndexes aD[]={pU, pL, pR}; vector<TwoIndexes> vD(aD,aD+3);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmax-1, m_win.rowmax, m_win.colmin+1, m_win.colmax-1, vD);

  // CORNERS
  // merge down-left pixels
  TwoIndexes aDL[]={pU, pR}; vector<TwoIndexes> vDL(aDL,aDL+2);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmax-1, m_win.rowmax, m_win.colmin, m_win.colmin+1, vDL);

  // merge down-right pixels
  TwoIndexes aDR[]={pU, pL}; vector<TwoIndexes> vDR(aDR,aDR+2);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmax-1, m_win.rowmax, m_win.colmax-1, m_win.colmax, vDR);

  // merge upper-left pixels
  TwoIndexes aUL[]={pD, pR}; vector<TwoIndexes> vUL(aUL,aUL+2);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmin, m_win.rowmin+1, m_win.colmin, m_win.colmin+1, vUL);

  // merge upper-right pixels
  TwoIndexes aUR[]={pD, pL}; vector<TwoIndexes> vUR(aUR,aUR+2);
  _mergePixelCouplesInArea(thr_on_max, thr_on_tot, DO_TEST, m_win.rowmin, m_win.rowmin+1, m_win.colmax-1, m_win.colmax, vUR);
}
  
//--------------------

void
AlgImgProc::_mergePixelCouplesInArea(
  const fphoton_t thr_on_max,
  const fphoton_t thr_on_tot,
  const bool      DO_TEST,
  const unsigned  rowmin,
  const unsigned  rowmax1,
  const unsigned  colmin,
  const unsigned  colmax1,
  const vector<TwoIndexes>& vneighbs)
{
  // loop over internal image pixels ignoring first and last rows and columns
  for(unsigned r = rowmin; r<rowmax1; r++) {
    for(unsigned c = colmin; c<colmax1; c++) {

      if(m_local_maximums(r,c) != 3)  continue; // check local maximums only
      if(m_fphoton(r,c) < thr_on_max) continue; // apply threshold on max

      fphoton_t vmax = 0.0;
      TwoIndexes ijmax(0,0);

      // find not-busy neighbor pixel with maximal intensity
      for(vector<TwoIndexes>::const_iterator ij  = vneighbs.begin();
                                             ij != vneighbs.end(); ij++) {
         int ir = r + (ij->i);
         int ic = c + (ij->j);

         if(m_conmap(ir,ic)) continue; // pixel is already used for other pair
         if(m_fphoton(ir,ic) < vmax) continue; // not a maximal neighbor

	 vmax = m_fphoton(ir,ic);
         ijmax = *ij;
      }

      fphoton_t  vtot = m_fphoton(r,c);
      if(vmax>0) vtot += m_fphoton(r + ijmax.i, c + ijmax.j);

      if(vtot < thr_on_tot) continue; // if pair intensity is below total threshold

      m_numreg ++;
      m_conmap(r,c) = m_numreg;
      m_conmap(r + ijmax.i, c + ijmax.j) = m_numreg;
      	
      // DO MERGE FOR A COUPLE OF SELECTED PIXELS      
      m_nphoton(r,c) ++; // increment number of photons

      if(DO_TEST) m_mphoton(r,c) = 1; // vtot*100;

    } // column loop
  } // row loop
}

//--------------------
//--------------------
//--------------------
//--------------------

void 
AlgImgProc::printMatrixOfRingIndexes()
{
  int indmax = (int)std::ceil(m_r0 + m_dr);
  int indmin = -indmax;
  unsigned counter = 0;
  std::stringstream ss; 
  ss << "printMatrixOfRingIndexes(), seg=" << m_seg << "  r0=" << m_r0 << "  dr=" << m_dr << '\n';

  for (int i = indmin; i <= indmax; ++ i) {
    for (int j = indmin; j <= indmax; ++ j) {

      float r = std::sqrt(float(i*i + j*j));
      int status = (r < m_r0 || r > m_r0 + m_dr) ? 0 : 1;
      if (status) counter++;
      if (i==0 && j==0) ss << " +";
      else              ss << " " << status;
    }
    ss << '\n';
  }
  ss << "Number of pixels to estimate background = " << counter << '\n';
  MsgLog(_name(), info, ss.str());
}

//--------------------

void 
AlgImgProc::printVectorOfRingIndexes()
{
  if(v_indexes.empty()) _evaluateRingIndexes(m_r0, m_dr);

  std::stringstream ss; 
  ss << "In printVectorOfRingIndexes:\n Vector size: " << v_indexes.size() << '\n';
  int n_pairs_in_line=0;
  for( vector<TwoIndexes>::const_iterator ij  = v_indexes.begin();
                                          ij != v_indexes.end(); ij++ ) {
    ss << " (" << ij->i << "," << ij->j << ")";
    if ( ++n_pairs_in_line > 9 ) {ss << "\n"; n_pairs_in_line=0;}
  }   
  MsgLog(_name(), info, ss.str());
}

//--------------------

void 
AlgImgProc::printMatrixOfDiagIndexes()
{
  int indmax =  m_rank;
  int indmin = -m_rank;

  std::stringstream ss; 
  ss << "In printMatrixOfDiagIndexes, seg=" << m_seg << "  rank=" << m_rank << '\n';

  for (int i = indmin; i <= indmax; ++ i) {
    for (int j = indmin; j <= indmax; ++ j) {
      int status = 1;
      if (i==0 || j==0) status = 0;
      //if (m_rank>2 && floor(std::sqrt(float(i*i + j*j)))>(int)m_rank) status = 0;
      if (i==0 && j==0) ss << " +";
      else              ss << " " << status;
    }
    ss << '\n';
  }

  MsgLog(_name(), info, ss.str());
}

//--------------------

void 
AlgImgProc::printVectorOfDiagIndexes()
{
  if(v_inddiag.empty()) _evaluateDiagIndexes(m_rank);

  std::stringstream ss; 
  ss << "In printVectorOfDiagIndexes:\n Vector size: " << v_inddiag.size() << '\n';
  int n_pairs_in_line=0;
  for( vector<TwoIndexes>::const_iterator ij  = v_inddiag.begin();
                                          ij != v_inddiag.end(); ij++ ) {
    ss << " (" << ij->i << "," << ij->j << ")";
    if ( ++n_pairs_in_line > 9 ) {ss << "\n"; n_pairs_in_line=0;}
  }   

  MsgLog(_name(), info, ss.str());
}

//--------------------

void 
AlgImgProc::_printStatisticsOfLocalExtremes()
{
  std::stringstream ss; 
  ss << "In _printStatisticsOfLocalExtremes: seg=" << m_seg << "  rank=" << m_rank << '\n';
  ss << "1=c 2=r 4=rect" << std::right << '\n';
  unsigned hismax[8] = {}; // all zeros
  unsigned hismin[8] = {}; // all zeros
  unsigned totmax = 0;
  unsigned totmin = 0;
  ndarray<pixel_maximums_t, 2>::iterator itx;
  ndarray<pixel_minimums_t, 2>::iterator itn;
  for(itx=m_local_maximums.begin(); itx!=m_local_maximums.end(); ++itx) {hismax[*itx]++; if(*itx) totmax++;}
  for(itn=m_local_minimums.begin(); itn!=m_local_minimums.end(); ++itn) {hismin[*itn]++; if(*itn) totmin++;}

  ss << "bin#    : "; for(int i=0; i<8; i++) ss << std::setw(8) << i;         ss << "     total\n";
  ss << "maximums: "; for(int i=0; i<8; i++) ss << std::setw(8) << hismax[i]; ss << "  " << std::setw(8) << totmax << '\n';
  ss << "minimums: "; for(int i=0; i<8; i++) ss << std::setw(8) << hismin[i]; ss << "  " << std::setw(8) << totmin << '\n';

  MsgLog(_name(), info, ss.str());
}

//--------------------
  std::ostream& 
  operator<<(std::ostream& os, const Peak& p) 
  {
    os << fixed
       << "Seg:"      << std::setw( 3) << std::setprecision(0) << p.seg
       << " Row:"     << std::setw( 4) << std::setprecision(0) << p.row 	     
       << " Col:"     << std::setw( 4) << std::setprecision(0) << p.col 	      
       << " Npix:"    << std::setw( 3) << std::setprecision(0) << p.npix    
       << " Imax:"    << std::setw( 7) << std::setprecision(1) << p.amp_max     	      
       << " Itot:"    << std::setw( 7) << std::setprecision(1) << p.amp_tot    	      
       << " CGrav r:" << std::setw( 6) << std::setprecision(1) << p.row_cgrav 	      
       << " c:"       << std::setw( 6) << std::setprecision(1) << p.col_cgrav   	      
       << " Sigma r:" << std::setw( 5) << std::setprecision(2) << p.row_sigma  	      
       << " c:"       << std::setw( 5) << std::setprecision(2) << p.col_sigma  	      
       << " Rows["    << std::setw( 4) << std::setprecision(0) << p.row_min    	      
       << ":"         << std::setw( 4) << std::setprecision(0) << p.row_max    	      
       << "] Cols["   << std::setw( 4) << std::setprecision(0) << p.col_min    	      
       << ":"         << std::setw( 4) << std::setprecision(0) << p.col_max    	     
       << "] B:"      << std::setw( 5) << std::setprecision(1) << p.bkgd       	      
       << " N:"       << std::setw( 5) << std::setprecision(1) << p.noise      	     
       << " S/N:"     << std::setw( 5) << std::setprecision(1) << p.son;
    return os;
  }

//--------------------
  std::ostream& 
  operator<<(std::ostream& os, const BkgdAvgRms& b) 
  {
    os << fixed
       << " Bkgd avg:" << std::setw( 7) << std::setprecision(1) << b.avg
       << " RMS:"      << std::setw( 7) << std::setprecision(1) << b.rms
       << " Npix:"     << std::setw( 4) << std::setprecision(0) << b.npx;
    return os;
  }

//--------------------
//--------------------
//-- NON-CLASS METHODS
//--------------------

ndarray<const AlgImgProc::pixel_maximums_t, 2>
mapOfLocalMaximumsRank1Cross(const ndarray<const AlgImgProc::fphoton_t,2> fphoton)
{
  AlgImgProc algo(0); // , 0, 1e6, 0, 1e6, 1023);
  algo.validate_window(fphoton.shape());
  algo._makeMapOfLocalMaximumsRank1Cross<AlgImgProc::fphoton_t>(fphoton);
  return algo.mapOfLocalMaximums();
}

//--------------------

template bool AlgImgProc::_findConnectedPixelsInRegionV1<float>(const ndarray<const float,2>&, const int&, const int&);
template bool AlgImgProc::_findConnectedPixelsInRegionV1<double>(const ndarray<const double,2>&, const int&, const int&);
template bool AlgImgProc::_findConnectedPixelsInRegionV1<int>(const ndarray<const int,2>&, const int&, const int&);
template bool AlgImgProc::_findConnectedPixelsInRegionV1<int16_t>(const ndarray<const int16_t,2>&, const int&, const int&);
template bool AlgImgProc::_findConnectedPixelsInRegionV1<uint16_t>(const ndarray<const uint16_t,2>&, const int&, const int&);

//--------------------

template void AlgImgProc::_findConnectedPixelsInRegionV2<float>(const ndarray<const float,2>&, const ndarray<const mask_t,2>&, const int&, const int&);
template void AlgImgProc::_findConnectedPixelsInRegionV2<double>(const ndarray<const double,2>&, const ndarray<const mask_t,2>&, const int&, const int&);
template void AlgImgProc::_findConnectedPixelsInRegionV2<int>(const ndarray<const int,2>&, const ndarray<const mask_t,2>&, const int&, const int&);
template void AlgImgProc::_findConnectedPixelsInRegionV2<int16_t>(const ndarray<const int16_t,2>&, const ndarray<const mask_t,2>&, const int&, const int&);
template void AlgImgProc::_findConnectedPixelsInRegionV2<uint16_t>(const ndarray<const uint16_t,2>&, const ndarray<const mask_t,2>&, const int&, const int&);

//--------------------
//--------------------
} // namespace ImgAlgos
//--------------------

