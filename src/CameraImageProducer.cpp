//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class CameraImageProducer...
//
// Author List:
//      Mikhail S. Dubrovin
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "ImgAlgos/CameraImageProducer.h"

//-----------------
// C/C++ Headers --
//-----------------
//#include <iomanip> // for setw, setfill
//#include <sstream> // for streamstring
//#include <iostream>// for setf

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "MsgLogger/MsgLogger.h"
#include "psddl_psana/camera.ddl.h"
#include "PSEvt/EventId.h"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

// This declares this class as psana module
using namespace ImgAlgos;
PSANA_MODULE_FACTORY(CameraImageProducer)

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

namespace ImgAlgos {

//----------------
// Constructors --
//----------------
CameraImageProducer::CameraImageProducer (const std::string& name)
  : Module(name)
  , m_str_src()
  , m_key_in()
  , m_key_out() 
  , m_subtract_offset()
  , m_print_bits()
  , m_count(0)
{
  // get the values from configuration or use defaults
  m_str_src           = configStr("source", "DetInfo(:Camera)");
  m_key_in            = configStr("key_in",                 "");
  m_key_out           = configStr("key_out",           "image");
  m_subtract_offset   = config   ("subtract_offset",      true);
  m_print_bits        = config   ("print_bits",             0 );
}

//--------------------

void 
CameraImageProducer::printInputParameters()
{
  WithMsgLog(name(), info, log) {
    log << "\n Input parameters :"
        << "\n source           : " << m_str_src
        << "\n key_in           : " << m_key_in      
        << "\n key_out          : " << m_key_out
        << "\n subtract_offset  : " << m_subtract_offset     
        << "\n print_bits       : " << m_print_bits
        << "\n";     
  }
}

//--------------------

//--------------
// Destructor --
//--------------
CameraImageProducer::~CameraImageProducer ()
{
}

/// Method which is called once at the beginning of the job
void 
CameraImageProducer::beginJob(Event& evt, Env& env)
{
  if( m_print_bits & 1 ) printInputParameters();
}

/// Method which is called at the beginning of the run
void 
CameraImageProducer::beginRun(Event& evt, Env& env)
{
}

/// Method which is called at the beginning of the calibration cycle
void 
CameraImageProducer::beginCalibCycle(Event& evt, Env& env)
{
}

/// Method which is called with event data, this is the only required 
/// method, all other methods are optional
void 
CameraImageProducer::event(Event& evt, Env& env)
{
  ++ m_count;
  if( m_print_bits & 2 ) printTimeStamp(evt);
  procEvent(evt,env);
}
  
/// Method which is called at the end of the calibration cycle
void 
CameraImageProducer::endCalibCycle(Event& evt, Env& env)
{
}

/// Method which is called at the end of the run
void 
CameraImageProducer::endRun(Event& evt, Env& env)
{
}

/// Method which is called once at the end of the job
void 
CameraImageProducer::endJob(Event& evt, Env& env)
{
  if( m_print_bits & 4 ) printSummary(evt);
}

//--------------------
//--------------------

void 
CameraImageProducer::procEvent(Event& evt, Env& env)
{
  shared_ptr<Psana::Camera::FrameV1> frmData = evt.get(m_str_src, m_key_in, &m_src);
  if (frmData.get()) {

    //unsigned h      = frmData->height();
    //unsigned w      = frmData->width();
    int offset = (m_subtract_offset) ? frmData->offset() : 0;

    //m_data = new double [h*w];
    double *p_data = &m_data_arr[0];
    unsigned ind = 0;

      const ndarray<uint8_t, 2>& data8 = frmData->data8();
      if (not data8.empty()) {
        if( m_print_bits & 8 ) MsgLog(name(), info, "procEvent(...): Get image as ndarray<uint8_t,2>, subtract offset=" << offset);
        ndarray<uint8_t, 2>::const_iterator cit;
        for(cit=data8.begin(); cit!=data8.end(); cit++) { p_data[ind++] = double(int(*cit) - offset); }

        saveImageInEvent(evt, p_data, data8.shape());
      }

      const ndarray<uint16_t, 2>& data16 = frmData->data16();
      if (not data16.empty()) {
        if( m_print_bits & 8 ) MsgLog(name(), info, "procEvent(...): Get image as ndarray<uint16_t,2>, subtract offset=" << offset);
        ndarray<uint16_t, 2>::const_iterator cit;
        // This loop consumes ~5 ms/event for Opal1000 camera with 1024x1024 image size 
        for(cit=data16.begin(); cit!=data16.end(); cit++) { p_data[ind++] = double(*cit) - offset; }

        saveImageInEvent(evt, p_data, data16.shape());
      }
  }
}

//--------------------

void 
CameraImageProducer::saveImageInEvent(Event& evt, double *p_data, const unsigned *shape)
{
  //m_img2d = new CSPadPixCoords::Image2D<double>(p_data, shape[0], shape[1]);
  shared_ptr< ndarray<double,2> > img2d( new ndarray<double,2>(p_data, shape) );
  evt.put(img2d, m_src, m_key_out);
}

//--------------------

void 
CameraImageProducer::printTimeStamp(Event& evt, std::string comment)
{
  shared_ptr<PSEvt::EventId> eventId = evt.get();
  if (eventId.get()) {

    MsgLog( name(), info, "Run="   <<  eventId->run()
                       << " Event=" <<  m_count 
                       << " Time="  <<  eventId->time()
                       << comment.c_str() );
  }
}

//--------------------

void 
CameraImageProducer::printSummary(Event& evt, std::string comment)
{
  MsgLog( name(), info, // " Run=" << eventId->run()
                       "Number of processed events=" <<  m_count 
                       << comment.c_str() );
}

//--------------------

} // namespace ImgAlgos