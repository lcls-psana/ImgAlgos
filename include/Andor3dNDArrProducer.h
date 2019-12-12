#ifndef PSANA_EXAMPLES_ANDOR3DNDARRPRODUCER_H
#define PSANA_EXAMPLES_ANDOR3DNDARRPRODUCER_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class Andor3dNDArrProducer.
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <sstream>  // for stringstream

//----------------------
// Base Class Headers --
//----------------------
#include "psana/Module.h"
#include "psddl_psana/andor3d.ddl.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ImgAlgos/GlobalMethods.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

namespace ImgAlgos {

/**
 *  @brief Psana module for access to Andor3d::FrameV# data
 *
 *  This software was developed for the LCLS project.  If you use all or 
 *  part of it, please give an appropriate acknowledgment.
 *
 *  @see AdditionalClass
 *
 *  @version $Id$
 *
 *  @author Mikhail Dubrovin
 */

class Andor3dNDArrProducer : public Module {
public:

  /// Data type for detector image 
  typedef uint16_t data_t;

  const static size_t   Segs   = 2; 
  const static size_t   Rows   = 2048; 
  const static size_t   Cols   = 2048; 
  const static size_t   FrSize = Rows*Cols; 
  const static size_t   Size   = Segs*Rows*Cols; 
   
  // Default constructor
  Andor3dNDArrProducer (const std::string& name) ;

  // Destructor
  virtual ~Andor3dNDArrProducer () ;

  /// Method which is called once at the beginning of the job
  virtual void beginJob(Event& evt, Env& env);

  /// Method which is called at the beginning of the calibration cycle
  virtual void beginCalibCycle(Event& evt, Env& env);
  
  /// Method which is called with event data
  virtual void event(Event& evt, Env& env);
  
protected:

  void printInputParameters();
  void procEvent(Event& evt, Env& env);
  void checkTypeImplementation();
 
private:

  Pds::Src    m_src;
  Source      m_str_src;
  std::string m_key_in; 
  std::string m_key_out;
  std::string m_outtype;
  unsigned    m_print_bits;

  DATA_TYPE   m_dtype;
 
//--------------------
  /**
   * @brief Get Andor3d from Psana::Andor3d::FramesV1 data object and copy them in the ndarray<TOUT, 3> out_ndarr
   * Returns false if data is missing.
   */

  template <typename TOUT>
  bool procEventForOutputType (Event& evt) {

      boost::shared_ptr<Psana::Andor3d::FrameV1> frame1 = evt.get(m_str_src, m_key_in, &m_src);
      if (frame1) {
	
          std::stringstream str; 
      
          const ndarray<const data_t, 3>& data = frame1->data();

          if( m_print_bits & 2 ) {      
            str << "\n Andor3d::FrameV1:";          
            str << "\n   shotIdStart = " << frame1->shotIdStart();
            str << "\n   readoutTime = " << frame1->readoutTime();
            str << "\n   temperature = ";
            const ndarray<const float, 1>& temps = frame1->temperature();
            for(ndarray<const float, 1>::iterator it=temps.begin(); it!=temps.end(); ++it) str << " " << *it;
	  }      

          if(m_dtype == ASDATA) {
             save3DArrInEvent<data_t>(evt, m_src, m_key_out, data);
	     return true;
	  }

	  //const unsigned shape = {Segs,Rows,Cols};data
	  ndarray<TOUT, 3> out_ndarr(data.shape());
          //ndarray<TOUT, 3> out_ndarr = make_ndarray<TOUT>(Segs,Rows,Cols);
          typename ndarray<TOUT, 3>::iterator it_out = out_ndarr.begin(); 

	  // Copy frame from data to output ndarray with changing type
          for ( ndarray<const data_t, 3>::iterator it=data.begin(); it!=data.end(); ++it, ++it_out) {
              *it_out = (TOUT)*it;
          }

          if( m_print_bits & 2 ) { str << "\n    out_ndarr:\n" << out_ndarr; MsgLog(name(), info, str.str() ); }

          save3DArrInEvent<TOUT>(evt, m_src, m_key_out, out_ndarr);
 
          return true;
      }
      else
      {
          if( m_print_bits & 16 ) MsgLog(name(), warning, "Andor3d::FramesV1 object is not available in the event(...) for source:"
              << m_str_src << " key:" << m_key_in);
	  return false;
      }

      return false;
  }

//-------------------

};

} // namespace ImgAlgos

#endif // PSANA_EXAMPLES_ANDOR3DNDARRPRODUCER_H
