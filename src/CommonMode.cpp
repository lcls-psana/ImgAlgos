//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class CommonMode...
//
// Author List:
//      Mikhail S. Dubrovin
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ImgAlgos/CommonMode.h"
 
//-----------------
// C/C++ Headers --
//-----------------
//#include <boost/lexical_cast.hpp>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
//#include <boost/shared_ptr.hpp>
//#include "PSEvt/EventId.h"
//#include "PSTime/Time.h"
//#include "psana/Module.h"
//#include "pdsdata/xtc/DetInfo.hh" // for srcToString( const Pds::Src& src )
//#include "PSCalib/Exceptions.h"   // for srcToString( const Pds::Src& src )

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

using namespace std;
using namespace ImgAlgos;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

namespace ImgAlgos {

  //typedef CommonMode::work_t work_t;

SingleStore* SingleStore::m_pInstance = NULL; // !!!!!! make global pointer !!!!!
 
//----------------
SingleStore::SingleStore()
{
  //print();
  make_ndarrays();
}

//----------------
SingleStore* SingleStore::instance()
{
  if(!m_pInstance) m_pInstance = new SingleStore();
  return m_pInstance;
}

//----------------
void SingleStore::print(){std::cout << "SingleStore::print() Single instance for singleton class SingleStore is created\n";}

//----------------
void SingleStore::make_ndarrays(){
      m_wasd = make_ndarray<work_t>(704, 768);
      m_wtrd = make_ndarray<work_t>(768, 704);
      m_ctrd = make_ndarray<work_t>(768, 704);
}

//--------------------
//--------------------
//--------------------

// Original version of the median estimation in integer numbers
  int median_for_hist_v1(const unsigned* hist, const int& low, const int& high, const unsigned& count) {
      int i=-1;
      int s = count/2;
      while(s>0) s -= hist[++i];
      if (unsigned(abs(-s)) > hist[i-1]/2) i--; // step back
      return low+i+1; // +1 is due to binning an empiric shift common mode to 0
  }

//--------------------

// 2016-04-14 Corrected version of the median estimation with float interpolation between integer bins.
  float median_for_hist(const unsigned* hist, const int& low, const int& high, const unsigned& count) {
      float halfst = (float)count/2;
      int i=-1;
      int s = (int)ceil(halfst);
      while(s>0) s -= hist[++i];
      float dx = float(s)/hist[i]; // dx - is a fraction of bin for float correction of median; presumably s<0, hist>0, so x<0
      return float(low+i)+dx+0.8;  // 0.8 - is an imperic number, which should be 1 due to binning, but 0.8 works better
  }

//--------------------
//--------------------
//--------------------

//template void ImgAlgos::commonMode<double> (double*  data, const uint16_t* mask, const unsigned length, const double  threshold, const double  maxCorrection, double&  cm);
//template void ImgAlgos::commonMode<float>  (float*   data, const uint16_t* mask, const unsigned length, const float   threshold, const float   maxCorrection, float&   cm);
//template void ImgAlgos::commonMode<int32_t>(int32_t* data, const uint16_t* mask, const unsigned length, const int32_t threshold, const int32_t maxCorrection, int32_t& cm);
//template void ImgAlgos::commonMode<int16_t>(int16_t* data, const uint16_t* mask, const unsigned length, const int16_t threshold, const int16_t maxCorrection, int16_t& cm);

//--------------------

} // namespace ImgAlgos

