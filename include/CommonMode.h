#ifndef IMGALGOS_COMMONMODE_H
#define IMGALGOS_COMMONMODE_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class CommonMode.
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------

#include <string>
#include <vector>
#include <fstream>   // ofstream
#include <iomanip>   // for setw, setfill
#include <sstream>   // for stringstream
#include <stdint.h>  // uint8_t, uint32_t, etc.
#include <iostream>
#include <cmath>
#include <algorithm> // for fill_n, std::sort

//#include <cstring>   // for memcpy
//#include <math.h>
//#include <stdio.h>

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "MsgLogger/MsgLogger.h"
#include "ndarray/ndarray.h"

namespace ImgAlgos {

using namespace std;

//--------------------
typedef float work_t;

//--------------------

class SingleStore{
public:
  static SingleStore* instance();
  void print();
  void make_ndarrays();

  // arrays for Epix100 common mode correction
  //const static unsigned shasd[2] = {704, 768};
  //const static unsigned shtrd[2] = {768, 704};

  ndarray<work_t,2> m_wasd; // (shasd); // work array shaped as data
  ndarray<work_t,2> m_wtrd; // (shtrd); // work array transposed from data
  ndarray<work_t,2> m_ctrd; // (shtrd); // array of corrected transposed data
 
private:

  SingleStore() ;                 // !!!!! Private so that it can not be called from outside
  virtual ~SingleStore(){};
 
  static SingleStore* m_pInstance; // !!!!! Singleton instance

  // Copy constructor and assignment are disabled by default
  SingleStore (const SingleStore&) ;
  SingleStore& operator = (const SingleStore&) ;
};

//--------------------

 const static int UnknownCM = -10000; 
 const static unsigned ROWS2X1   = 185; 
 const static unsigned COLS2X1   = 388; 
 const static unsigned SIZE2X1   = ROWS2X1*COLS2X1; 
 const static unsigned COLSHALF  = COLS2X1/2;

 int median_for_hist_v1(const unsigned* hist, const int& low, const int& high, const unsigned& count);
 float  median_for_hist(const unsigned* hist, const int& low, const int& high, const unsigned& count);

//--------------------
// Code from ami/event/FrameCalib.cc
// int median(ndarray<const uint32_t,2> d, unsigned& iLo, unsigned& iHi);

//--------------------

//--------------------
// Philip Hart's code from psalg/src/common_mode.cc:
//
// averages all values in array "data" for length "length"
// that are below "threshold".  If the magniture of the correction
// is greater than "maxCorrection" leave "data" unchanged, otherwise
// subtract the average in place.  
// mask is either a null pointer (in which case nothing is masked)
// or a list of values arranged like the data where non-zero means ignore

   template <typename T>
   void commonMode(T* data, const uint16_t* mask, const unsigned length, const T threshold, const T maxCorrection, T& cm) {
     // do dumbest thing possible to start - switch to median
     // do 2nd dumbest thing possible to start - switch to median
     cm = 0;
     T* tmp = data;
     double sum = 0;
     int nSummed = 0;
     for (unsigned col=0; col<length; col++) {
       T cval = *tmp++;
       T mval = (mask) ? *mask++ : 0;
       if (mval==0 && cval<threshold) {
         nSummed++;
         sum += cval;
       }
     }
     
     if (nSummed>0) {
       T mean = (T)(sum/nSummed);
       if (fabs(mean)<=maxCorrection) {
         cm = mean;
         tmp = data;
         for (unsigned col=0; col<length; col++) {
   	   *tmp++ -= cm;
         }
       }
     }
   }

//--------------------
//--------------------
// Philip Hart's code from psalg/src/common_mode.cc:
//
// finds median of values in array "data" for length "length"
// that are below "threshold".  If the magniture of the correction
// is greater than "maxCorrection" leave "data" unchanged, otherwise
// subtract the median in place.  
// mask is either a null pointer (in which case nothing is masked)
// or a list of values arranged like the data where non-zero means ignore

template <typename T>
void commonModeMedian(const T* data, const uint16_t* mask, const unsigned length, const T threshold, const T maxCorrection, T& cm) {
  cm = 0;
  const T* tmp = data;

  const uint32_t lMax = 32768;//2**14*2;  // here I may be assuming data in ADU
  const uint32_t lHalfMax = 16384;//2**14;
  const int      iMax = (int)lMax;
  const int      iHalfMax = (int)lHalfMax;
  unsigned hist[lMax];
  //memset(hist, 0, sizeof(unsigned)*lMax);
  std::fill_n(&hist[0], int(lMax), unsigned(0));
  int nSummed = 0;

  for (unsigned col=0; col<length; col++) {
    T cval = *tmp++;
    T mval = (mask) ? *mask++ : 0;
    if (mval==0 && cval<threshold) {
      nSummed++;
      //unsigned bin = (int)cval+lHalfMax;
      // unsigned long?  check range or raise?
      int bin = (int)cval+lHalfMax;
      if (bin<0) bin = 0;
      if (!(bin<iMax)) bin = lMax-1;
      hist[bin]++;
    }
  }

  if (nSummed==0) return;

  unsigned medianCount = (unsigned)ceil(nSummed/2.);
  unsigned histSum = 0;
  for (int bin=0; bin<iMax; bin++) {
    histSum += hist[bin];
    if (histSum>=medianCount) {
      T median = (int)bin - iHalfMax;
      if (fabs(median)<=maxCorrection) {
        cm = median;
      }
      return;
    }
  }
}

//--------------------
// Apply "median" common mode correction.
// Signature of this function has non-const data.
template <typename T>
void commonModeMedian(T* data, const uint16_t* mask, const unsigned length, const T threshold, const T maxCorrection, T& cm) {
  commonModeMedian((const T*)data, mask, length, threshold, maxCorrection, cm);
  if (cm != 0) {
    T* tmp = data;
    for (unsigned col=0; col<length; col++) {
      *tmp++ -= cm;
    }
  }
}

//--------------------
//--------------------
// This method was originally designed by Andy for CSPAD in pdscalibdata/src/CsPadCommonModeSubV1.cpp
// - data type int16_t is changed to T
// - subtract CM inside this module
// - pedestal is not subtracted in this algorithm; assume that it is already subtracted

  /**
   *  Find common mode for an CsPad  section.
   *  
   *  Function will return UnknownCM value if the calculation 
   *  cannot be performed (or need not be performed).
   *  
   *  @param pars   array[3] of control parameters; mean_max, sigma_max, threshold on number of pixels/ADC count
   *  @param sdata  pixel data
   *  @param pixStatus  pixel status data, can be zero pointer
   *  @param ssize  size of all above arrays
   *  @param stride increment for pixel indices
   */ 

   //float cm_corr =  findCommonMode<T>(pars, sdata, pixStatus, ssize, nSect); 

//--------------------
// Modified Philip Hart's commonMode (mean) code (see below) for 2-d region in data. 

template <typename T>
  void meanInRegion(  const double* pars
                    , ndarray<T,2>& data 
	            , ndarray<const uint16_t,2>& status
	            , const size_t& rowmin
                    , const size_t& colmin
                    , const size_t& nrows
                    , const size_t& ncols
                    , const size_t& srows = 1
                    , const size_t& scols = 1
                    , const unsigned& pbits=0
                    ) {

     T threshold = pars[1];
     T maxcorr   = pars[2];

     double sumv = 0;
     int    sum1 = 0;
     bool   check_status = (status.data()) ? true : false;

     for (size_t r=rowmin; r<rowmin+nrows; r+=srows) { 
       for (size_t c=colmin; c<colmin+ncols; c+=scols) {
         T v = data(r, c);
         T s = (check_status) ? status(r, c) : 0;
         if (s==0 && v<threshold) {
           sumv += v;
           sum1 ++;
         }
       }
     }
     
     if (sum1>0) {
       T mean = (T)(sumv/sum1);

       if (pbits) std::cout << "  mean:" << mean
                            << "  threshold:" << threshold
                            << '\n';

       if (fabs(mean)<=maxcorr) {
         for (size_t r=rowmin; r<rowmin+nrows; r+=srows) { 
           for (size_t c=colmin; c<colmin+ncols; c+=scols) {
             data(r, c) -= mean;
           }
         }
       }
     }
}
//--------------------
// Median, similar to ami/event/FrameCalib::median, 
// but finding median for entire good statistics w/o discarding edge bins

  template <typename T>
  void medianInRegion(const double* pars
                    , ndarray<T,2>& data 
	            , ndarray<const uint16_t,2>& status
	            , const size_t& rowmin
                    , const size_t& colmin
                    , const size_t& nrows
                    , const size_t& ncols
                    , const size_t& srows = 1
                    , const size_t& scols = 1
		    , const unsigned& pbits=0
                    ) {

      int hint_range = (int)pars[1];
      int maxcorr    = (int)pars[2];

      bool check_status = (status.data()) ? true : false;

      // declare array for histogram
      int half_range = max(int(hint_range), 10);

      int low  = -half_range;
      int high =  half_range;
      unsigned* hist  = 0;
      int       nbins = 0;
      int       bin   = 0;
      unsigned  iter  = 0;
      unsigned  count = 0;

      while(1) { // loop continues untill the range does not contain a half of statistics

	  iter ++;

          nbins = high-low+1;
	  if (nbins>10000) {
            if (pbits & 4) MsgLog("medianInRegion", warning, "Too many bins " << nbins 
                                  << ", common mode correction is not allied");
	    return;
	  }

          if (hist) delete[] hist;
          hist = new unsigned[nbins];
          std::fill_n(hist, nbins, 0);
          count = 0;
      
          // fill histogram
          for (size_t r=rowmin; r<rowmin+nrows; r+=srows) { 
            for (size_t c=colmin; c<colmin+ncols; c+=scols) {
            
              // ignore pixels that are too noisy; discard pixels with any status > 0
              if (check_status && status(r, c)) continue;
            
              bin = int(data(r, c)) - low;
              if      (bin < 1)     hist[0]++;
              else if (bin < nbins) hist[bin]++;
              else                  hist[nbins-1]++;
              count++;
            }
          }
          
          if (pbits & 2) {
              MsgLog("medianInRegion", info, "Iter:" << iter 
                                   << "  histo nbins:" << nbins 
                                   << "  low:" << low 
                                   << "  high:" << high 
                                   << "  count:" << count);
              for (int b=0; b<nbins; ++b) std::cout << " " << b << ":" << hist[b]; 
              std::cout  << '\n';
          }
    
          // do not apply correction if the number of good pixels is small
          if (count < 10) {delete[] hist; return;}

          if (hist[0]>count/2) {
	    if (maxcorr && low < -maxcorr) {delete[] hist; return;} // do not apply cm correction
            low  -= nbins/4;
	  }
          else if (hist[nbins-1]>count/2) {
	    if (maxcorr && high > maxcorr) {delete[] hist; return;} // do not apply cm correction
	    high += nbins/4;
	  }
	  else
            break; 
      } // while(1)

      //--------------------

      T cm = (T)median_for_hist(hist, low, high, count);

      if (maxcorr && fabs(cm)>maxcorr) {delete[] hist; return;} // do not apply cm correction

      if (pbits & 1) MsgLog("medianInRegionV2", info, "cm correction = " << cm);

      // Apply common mode correction to data
      for (size_t r=rowmin; r<rowmin+nrows; r+=srows) { 
        for (size_t c=colmin; c<colmin+ncols; c+=scols) {
          data(r, c) -= cm;
        }
      }

      delete[] hist; 
  }

//--------------------
//--------------------
// Median suggested by Silke on 2016/04/12 for Epix100a

  template <typename T>
  void medianInRegionV2(const double* pars
                    , ndarray<T,2>& data 
	            , ndarray<const uint16_t,2>& status
	            , const size_t& rowmin
                    , const size_t& colmin
                    , const size_t& nrows
                    , const size_t& ncols
                    , const size_t& srows = 1
                    , const size_t& scols = 1
		    , const unsigned& pbits=0
                    ) {

    //static unsigned nentry=0; nentry++;
    //if(nentry<2) cout << "medianInRegionV2\n";

      int hint_range = (int)pars[1];
      int maxcorr    = (int)pars[2];

      bool check_status = (status.data()) ? true : false;

      // declare array for histogram
      int half_range = max(int(hint_range), 10);

      int       low   = -half_range;
      int       high  =  half_range;
      int       bin   = 0;
      unsigned  count = 0;
      int       nbins = high-low+1;
      unsigned* hist  = new unsigned[nbins];
      std::fill_n(hist, nbins, 0);
      
      // fill histogram
      for (size_t r=rowmin; r<rowmin+nrows; r+=srows) { 
        for (size_t c=colmin; c<colmin+ncols; c+=scols) {
        
          // ignore pixels that are too noisy; discard pixels with any status > 0
          if (check_status && status(r, c)) continue;
        
          bin = floor(data(r, c)) - low;
          if      (bin < 1)     hist[0]++;
          else if (bin < nbins) hist[bin]++;
          else                  hist[nbins-1]++;
          count++;
        }
      }
      
      if (pbits & 2) {
          MsgLog("medianInRegionV2", info,
                                  "  histo nbins:" << nbins 
                               << "  low:" << low 
                               << "  high:" << high 
                               << "  count:" << count);
          for (int b=0; b<nbins; ++b) std::cout << " " << b << ":" << hist[b]; 
          std::cout  << '\n';
      }
    
      // do not apply correction if
      if (count < 10            // the number of good pixels is small
      or  hist[0]>count/2        // more than half statistics in the edge bin
      or  hist[nbins-1]>count/2) {
          delete[] hist; return; 
      }

      T cm = (T)median_for_hist(hist, low, high, count);

      if (maxcorr && fabs(cm)>maxcorr) {delete[] hist; return;} // do not apply cm correction

      if (pbits & 1) MsgLog("medianInRegionV2", info, "cm correction = " << cm);

      // Apply common mode correction to data
      for (size_t r=rowmin; r<rowmin+nrows; r+=srows) { 
        for (size_t c=colmin; c<colmin+ncols; c+=scols) {
          data(r, c) -= cm;
        }
      }

      delete[] hist; 
  }

//--------------------

  template <typename T>
  T median(std::vector<T>& v)
  {
      size_t middle = v.size() / 2;
      std::nth_element(v.begin(), v.begin() + middle, v.end());

      T median = v[middle];

      // even number of elements
      if(!(v.size() & 1)) {
          typename std::vector<T>::iterator max_it;
          max_it = std::max_element(v.begin(), v.begin() + middle);
          median = 0.5*(*max_it + median);
      }
      return median;
  }

//--------------------
//--------------------
// Another median algorithm

  template <typename T>
  void medianInRegionV3(const double* pars
                    , ndarray<T,2>& data 
	            , ndarray<const uint16_t,2>& status
	            , const size_t& rowmin
                    , const size_t& colmin
                    , const size_t& nrows
                    , const size_t& ncols
                    , const size_t& srows = 1
                    , const size_t& scols = 1
		    , const unsigned& pbits=0
                    ) {

    //static unsigned nentry=0; nentry++;
    //if(nentry<2) cout << "medianInRegionV2\n";

      T half_range = (T)pars[1];
      T maxcorr    = (T)pars[2];
      T d = 0;

      bool check_status = (status.data()) ? true : false;

      std::vector<T> vec;
      vec.reserve(nrows*ncols);

      // fill vector
      for (size_t r=rowmin; r<rowmin+nrows; r+=srows) { 
        for (size_t c=colmin; c<colmin+ncols; c+=scols) {
          if (check_status && status(r, c)) continue;
          d = data(r, c);
	  if(d >  half_range) continue;
	  if(d < -half_range) continue;	  
          vec.push_back(d);
        }
      }

      if (vec.size() < 10) return;

      T cm = median<T>(vec);

      if (maxcorr && fabs(cm)>maxcorr) return; // do not apply cm correction

      if (pbits & 1) MsgLog("medianInRegionV3", info, "vec.size = " << vec.size() << "  cm correction = " << cm);
      //cout << "medianInRegionV3 vec.size = " << vec.size() << "  cm correction = " << cm << "\n";

      // Apply common mode correction to data
      for (size_t r=rowmin; r<rowmin+nrows; r+=srows) { 
        for (size_t c=colmin; c<colmin+ncols; c+=scols) {
          data(r, c) -= cm;
        }
      }
  }

//--------------------
  /**
   *  Find common mode for CSPAD 2x1 section using unbond pixels.
   *
   *  @param pars   array[1] of control parameters; mean_max - maximal allowed correction
   *  @param sdata  pixel data
   *  @param ssize  size of data array (188*388)
   *  @param stride increment for pixel indices
   */ 
   //float cm_corr = applyCModeUnbond<T>(pars, sdata, ssize, stride); 

template <typename T>
float 
applyCModeUnbond( const double* pars,
                  T* sdata,
                  unsigned ssize,
                  int stride = 1
                )
{
  int    npix = 0;
  double mean = 0;
  int    ind  = 0;

  //skip p=0 in the corner
  for (size_t r=10; r<ROWS2X1; r+=10) {      
    ind = (r*COLS2X1 + r)*stride;
    mean += sdata[ind] + sdata[ind+COLSHALF*stride];
    npix += 2; 
  }

  mean = (npix>0) ? mean/npix : 0;

  //std::cout << "XXX: applyCModeUnbond mean = " << mean << " evaluated for npix = " << npix << '\n';

  // limit common mode correction to some reasonable numbers
  if (abs(mean) > pars[0]) return float(UnknownCM);

  //--------------------
  // subtract CM 
  for (unsigned c=0, p=0; c<ssize; ++c, p += stride) {
        sdata[p] -= mean;
  }
  return (float)mean;
}

template <typename T>
float 
findCommonMode(const double* pars,
               T* sdata,
               const uint16_t *pixStatus, 
               unsigned ssize,
               int stride = 1
               )
{
  // do we even need it
  //if (m_mode == None) return float(UnknownCM);

  // for now it does not make sense to calculate common mode
  // if pedestals are not known
  // if (not peddata) return float(UnknownCM);
  
  // declare array for histogram
  const int low = -1000;
  const int high = 2000;
  const unsigned hsize = high-low;
  int hist[hsize];
  std::fill_n(hist, hsize, 0);
  unsigned long count = 0;

  // fill histogram
  for (unsigned c = 0, p = 0; c != ssize; ++ c, p += stride) {
    // ignore channels that re too noisy
    //if (pixStatus and (pixStatus[p] & 1)) continue;
    if (pixStatus and pixStatus[p]) continue; // Discard  pixels with any status > 0

    // pixel value with pedestal subtracted, rounded to integer
    double dval = sdata[p]; // - peddata[p];
    int val = dval < 0 ? int(dval - 0.5) : int(dval + 0.5);
    
    // histogram bin
    unsigned bin = unsigned(val - low);
    
    // increment bin value if in range
    if (bin < hsize) {
      ++hist[bin] ;
      ++ count;
    }      
  }
  MsgLog("findCommonMode", debug, "histo filled count = " << count);
  
  // analyze histogram now, first find peak position
  // as the position of the lowest bin with highest count 
  // larger than 100 and which has a bin somewhere on 
  // right side with count dropping by half
  int peakPos = -1;
  int peakCount = -1;
  int hmRight = hsize;
  int thresh = int(pars[2]);
  if(thresh<=0) thresh=100;
  for (unsigned i = 0; i < hsize; ++ i ) {
    if (hist[i] > peakCount and hist[i] > thresh) {
      peakPos = i;
      peakCount = hist[i];
    } else if (peakCount > 0 and hist[i] <= peakCount/2) {
      hmRight = i;
      break;
    }
  }

  // did we find anything resembling
  if (peakPos < 0) {
    MsgLog("findCommonMode", debug, "peakPos = " << peakPos);
    return float(UnknownCM);
  }

  // find half maximum channel on left side
  int hmLeft = -1;
  for (int i = peakPos; hmLeft < 0 and i >= 0; -- i) {
    if(hist[i] <= peakCount/2) hmLeft = i;
  }
  MsgLog("findCommonMode", debug, "peakPos = " << peakPos << " peakCount = " << peakCount 
      << " hmLeft = " << hmLeft << " hmRight = " << hmRight);
  
  // full width at half maximum
  int fwhm = hmRight - hmLeft;
  double sigma = fwhm / 2.36;

  // calculate mean and sigma
  double mean = peakPos;
  for (int j = 0; j < 2; ++j) {
    int s0 = 0;
    double s1 = 0;
    double s2 = 0;
    int d = int(sigma*2+0.5);
    for (int i = std::max(0,peakPos-d); i < int(hsize) and i <= peakPos+d; ++ i) {
      s0 += hist[i];
      s1 += (i-mean)*hist[i];
      s2 += (i-mean)*(i-mean)*hist[i];
    }
    mean = mean + s1/s0;
    sigma = std::sqrt(s2/s0 - (s1/s0)*(s1/s0));
  }
  mean += low;

  MsgLog("findCommonMode", debug, "mean = " << mean << " sigma = " << sigma);

  if (abs(mean) > pars[0] or sigma > pars[1]) {
    return float(UnknownCM);
  }
  //--------------------
  // subtract CM 
  for (unsigned c = 0, p = 0; c < ssize; ++ c, p += stride) {
        sdata[p] -= mean;
  }

  return mean;
}
//--------------------


//--------------------
// Optimized common mode correction algorithm for Epix100
//
// Optimizetion
// ------------
// 1. remove redundant loops comparing to medianInRegionV3
// 2. work with float values in stead of double
// 3. do not use vector-s

  template <typename T>
  void medianEpix100V1(const double* pars
                      ,ndarray<T,2>& data 
	              ,ndarray<const uint16_t,2>& status
		      ,const unsigned& pbits=0
                      ) {

    const unsigned shasd[2] = {704, 768};
    //const unsigned shtrd[2] = {768, 704};

    static unsigned nentry=0; nentry++;
    if(nentry==1) {
      //cout << "medianEpix100V1: cmpars = "; for(int i=0; i<3; i++) cout << " " << pars[i]; cout << '\n';
      //cout << "medianEpix100V1: sizeof(float)=" << sizeof(float) << '\n'; // the answer is 4 byte
      if (pbits & 1) MsgLog("medianEpix100V1", info, "cmpars = " << pars[0] << " " << pars[1] << " " << pars[2]);

      //SingleStore::instance()->print()
    }

    ndarray<work_t,2>& m_wasd = SingleStore::instance()->m_wasd;
    ndarray<work_t,2>& m_wtrd = SingleStore::instance()->m_wtrd;
    ndarray<work_t,2>& m_ctrd = SingleStore::instance()->m_ctrd;

    const work_t BADDATA = -100000.;

    //size_t nregs  = 16;	
    size_t nrows  = shasd[0]; // 704
    size_t ncols  = shasd[1]; // 768
    T half_range  = (T)pars[1];
    T maxcorr     = (work_t)pars[2];

    bool check_status = (status.data()) ? true : false;

    typename ndarray<T, 2>::iterator itd;
    typename ndarray<const uint16_t, 2>::iterator its=status.begin();
    typename ndarray<work_t, 2>::iterator itw=m_wasd.begin();

    // --------------------------------
    // Correction of 96 pixels in rows
    // --------------------------------

    // conditional copy of data to work array in right order
    for(itd=data.begin(); itd!=data.end(); itd++, its++, itw++) {
      bool bad_pixel = (check_status && *its) or (*itd>half_range) or (*itd<-half_range);
      *itw = (bad_pixel) ? BADDATA : (work_t)(*itd);
    }

    // loop over 1-d row parts in 2-d m_wasd array
    size_t grlen = 96; // ncols/8
    itd=data.begin();
    for(itw=m_wasd.begin(); itw!=m_wasd.end(); itw+=grlen, itd+=grlen) {

      // Edges of the row
      work_t* pbeg = itw;
      work_t* pend = itw+grlen; // assuming as usually [pbe,pend)

      // Sort elements
      std::sort(pbeg,pend);
      //cout << "\nXXX sorted:"; for(work_t* p=pbeg; p<pend; p++) cout << " " << *p; cout << '\n'; 

      // Discard bad pixels 
      for(; pbeg!=pend; pbeg++)
        if(*pbeg!=BADDATA) break;
      //cout << "\nXXX clean sorted:"; for(work_t* p=pbeg; p<pend; p++) cout << " " << *p; cout << '\n'; 

      int npix_good = pend-pbeg;
      if(npix_good<10) continue; // do not apply correction for small number of good pixels
      work_t cm = *(pbeg + npix_good/2);
      //cout << "\nXXX npix_good=" << npix_good << "  cm=" << cm << '\n'; 

      if(maxcorr && fabs(cm)>maxcorr) continue;  // do not apply large correction
      
      // apply correction to data part of the row
      for(T* p=itd; p!=itd+grlen; p++) *p -= cm;
    }

    // --------------------------------
    // Correction of 352 pixels in cols
    // --------------------------------

    // data has changed, so need to repeat conditional copy array again
    its=status.begin();
    itw=m_wasd.begin();
    for(itd=data.begin(); itd!=data.end(); itd++, its++, itw++) {
      bool bad_pixel = (check_status && *its) or (*itd>half_range) or (*itd<-half_range);
      *itw = (bad_pixel) ? BADDATA : (work_t)(*itd);
    }

    // Fill in transposed array
    for(size_t r=0; r<nrows; r++)
      for(size_t c=0; c<ncols; c++) {
        m_wtrd(c, r) = m_wasd(r, c);
        m_ctrd(c, r) = data(r, c);
      }

    // loop over 1-d row parts in 2-d m_wtrd array
    grlen = 352; // nrows/2
    typename ndarray<work_t, 2>::iterator itc=m_ctrd.begin();
    typename ndarray<work_t, 2>::iterator itr=m_wtrd.begin();
    for(itr=m_wtrd.begin(); itr!=m_wtrd.end(); itr+=grlen, itc+=grlen) {

      // Edges of the row
      work_t* pbeg = itr;
      work_t* pend = itr+grlen; // assuming as usually [pbe,pend)

      // Sort elements
      std::sort(pbeg,pend);
      //cout << "\nXXX sorted:"; for(work_t* p=pbeg; p<pend; p++) cout << " " << *p; cout << '\n'; 

      // Discard bad pixels 
      for(; pbeg!=pend; pbeg++)
        if(*pbeg!=BADDATA) break;
      //cout << "\nXXX clean sorted:"; for(work_t* p=pbeg; p<pend; p++) cout << " " << *p; cout << '\n'; 

      int npix_good = pend-pbeg;
      if(npix_good<10) continue; // do not apply correction for small number of good pixels
      work_t cm = *(pbeg + npix_good/2);
      //cout << "\nXXX npix_good=" << npix_good << "  cm=" << cm << '\n'; 

      if(maxcorr && fabs(cm)>maxcorr) continue;  // do not apply large correction
      
      // apply correction to transposed data part of the row
      for(work_t* p=itc; p!=itc+grlen; p++) *p -= cm;
    }

    // copy transposed corrected array to data
    for(size_t r=0; r<nrows; r++)
      for(size_t c=0; c<ncols; c++)
        data(r, c) = (T)m_ctrd(c, r);
  }

//--------------------

} // namespace ImgAlgos

#endif // IMGALGOS_COMMONMODE_H
