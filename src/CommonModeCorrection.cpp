//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class CommonModeCorrection
//
// Author List:
//      Mikhail S. Dubrovin
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ImgAlgos/CommonModeCorrection.h"

//-----------------
// C/C++ Headers --
//-----------------
//#include <math.h>    // for exp
#include <iomanip>   // for std::setw
#include <sstream>   // for stringstream

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "MsgLogger/MsgLogger.h"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

namespace ImgAlgos {

typedef CommonModeCorrection::pixel_status_t pixel_status_t;
typedef CommonModeCorrection::common_mode_t  common_mode_t;

//----------------
// Constructors --
//----------------
//template <typename T>
//CommonModeCorrection<T>::CommonModeCorrection(const double& sigma, const int& nsm) 

CommonModeCorrection::CommonModeCorrection (
    const PSEvt::Source& source, 
    const common_mode_t* cmod_pars, 
    const unsigned size, 
    const pixel_status_t* status, 
    const unsigned pbits)
  : m_source(source)
  //, m_cmod_pars(cmod_pars)
  , m_size(size)
  , m_status(status)
  , m_pbits(pbits)
{
  //m_cmod_pars = new common_mode_t[16];
  //memcpy(m_cmod_pars, cmod_pars, 16*sizeof(common_mode_t));

  m_cmod_pars = cmod_pars;

  m_dettype = detectorTypeForSource(m_source);
  if(m_pbits & 1) printInputPars();
}

//--------------------
// Print member data and matrix of weights
void 
CommonModeCorrection::printInputPars()
{
  std::stringstream ss; 
  ss << "\nsource   : " << m_source
     << "\ndettype  : " << m_dettype
     << "\nsize     : " << m_size
     << "\npbits    : " << m_pbits
     << "\ncmod_pars: " << m_cmod_pars[0]
     <<            ", " << m_cmod_pars[1]
     <<            ", " << m_cmod_pars[2]
     <<            "...\n";
  MsgLog(_name_(), info, ss.str()); 
}

//--------------------
void
CommonModeCorrection::setCModPars(const common_mode_t* cmod_pars)
{ 
  m_cmod_pars = cmod_pars;
  //std::cout << "XXX:TEST CommonModeCorrection::setCModPars\n"; 
  ////memcpy(m_cmod_pars, cmod_pars, 16*sizeof(common_mode_t));
  if(m_pbits & 128) {
      std::stringstream ss; ss << "setCModPars: ";
      for(int i=0; i<4; ++i) ss << " " << m_cmod_pars[i];
      MsgLog(_name_(), info, ss.str()); 
  }
}

//--------------------

//template class ImgAlgos::CommonModeCorrection<int16_t>;
//template class ImgAlgos::CommonModeCorrection<int>;
//template class ImgAlgos::CommonModeCorrection<float>;
//template class ImgAlgos::CommonModeCorrection<double>;

//--------------------

} // namespace ImgAlgos
