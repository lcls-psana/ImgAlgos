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

  v_peaks_work.reserve(m_npksmax);
  v_peaks.reserve(m_npksmax);
  v_peaks_sel.reserve(m_npksmax);
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

  for(unsigned r = m_win.rowmin; r<m_win.rowmax; r++) {
    for(unsigned c = m_win.colmin; c<m_win.colmax; c++) {

      if(!(m_pixel_status[r][c] & 1)) continue;
      ++ m_numreg;
      //if(m_numreg == m_npksmax) break;

      _findConnectedPixels(r, c);
    }
  }

  if(m_pbits & 512) MsgLog(_name(), info, "in _makeMapOfConnectedPixels, seg=" << m_seg << ", m_numreg=" <<  m_numreg);
}

//--------------------

void 
AlgImgProc::_findConnectedPixels(const unsigned& r, const unsigned& c)
{
  //if(m_pbits & 512) MsgLog(_name(), info, "in _findConnectedPixels, seg=" << m_seg);

  if(! (m_pixel_status[r][c] & 1)) return;

  m_pixel_status[r][c] ^= 1; // set the 1st bit to zero.
  m_conmap[r][c] = m_numreg;

  if(  r+1 < m_win.rowmax  ) _findConnectedPixels(r+1, c);
  if(  c+1 < m_win.colmax  ) _findConnectedPixels(r, c+1);
  if(!(r-1 < m_win.rowmin) ) _findConnectedPixels(r-1, c);
  if(!(c-1 < m_win.colmin) ) _findConnectedPixels(r, c-1);  
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
  if (peak.son < m_peak_son_min) return false;

  if (! m_do_preselect) {
    if (peak.npix    < m_peak_npix_min) return false;
    if (peak.npix    > m_peak_npix_max) return false;
    if (peak.amp_max < m_peak_amax_thr) return false;
    if (peak.amp_tot < m_peak_atot_thr) return false;
  }
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
  v_inddiag.clear();

  int indmax =  m_rank;
  int indmin = -m_rank;

  for (int i = indmin; i <= indmax; ++ i) {
    for (int j = indmin; j <= indmax; ++ j) {

      // use rectangular region of radius = rank
      // remove already tested central row and column
      if (i==0 || j==0) continue;
      // use ring region (if un-commented)
      //if (m_rank>2 && floor(std::sqrt(float(i*i + j*j)))>(int)m_rank) continue;
      TwoIndexes inds = {i,j};
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

  v_indexes.clear();

  int indmax = (int)std::ceil(m_r0 + m_dr);
  int indmin = -indmax;

  for (int i = indmin; i <= indmax; ++ i) {
    for (int j = indmin; j <= indmax; ++ j) {

      float r = std::sqrt( float(i*i + j*j) );
      if ( r < m_r0 || r > m_r0 + m_dr ) continue;
      TwoIndexes inds = {i,j};
      v_indexes.push_back(inds);
    }
  }

  if(m_pbits & 2) printMatrixOfRingIndexes();
  if(m_pbits & 4) printVectorOfRingIndexes();
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
AlgImgProc::printMatrixOfRingIndexes()
{
  int indmax = (int)std::ceil(m_r0 + m_dr);
  int indmin = -indmax;
  unsigned counter = 0;
  std::stringstream ss; 
  ss << "printMatrixOfRingIndexes(), seg=" << m_seg << "  r0=" << m_r0 << "  dr=" << m_dr << '\n';

  for (int i = indmin; i <= indmax; ++ i) {
    for (int j = indmin; j <= indmax; ++ j) {

      float r = std::sqrt( float(i*i + j*j) );
      int status = ( r < m_r0 || r > m_r0 + m_dr ) ? 0 : 1;
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
  ss << "printVectorOfRingIndexes():\n Vector size: " << v_indexes.size() << '\n';
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
  ss << "printMatrixOfDiagIndexes(), seg=" << m_seg << "  rank=" << m_rank << '\n';

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
  ss << "printVectorOfDiagIndexes():\n Vector size: " << v_inddiag.size() << '\n';
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
  ss << "_printStatisticsOfLocalExtremes(): seg=" << m_seg << "  rank=" << m_rank << '\n';
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
} // namespace ImgAlgos
//--------------------

