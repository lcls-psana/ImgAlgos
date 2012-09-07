#ifndef IMGALGOS_IMGVSTIMESPLITINFILES_H
#define IMGALGOS_IMGVSTIMESPLITINFILES_H

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id$
//
// Description:
//	Class ImgVsTimeSplitInFiles.
//
//------------------------------------------------------------------------

//-----------------
// C/C++ Headers --
//-----------------
#include <iostream>
#include <fstream> // for std::ofstream operator << 
#include <sstream> // for stringstream 

//----------------------
// Base Class Headers --
//----------------------
#include "psana/Module.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

#include "PSEvt/Source.h"


//		---------------------
// 		-- Class Interface --
//		---------------------

namespace ImgAlgos {

/// @addtogroup ImgAlgos

/**
 *  @ingroup ImgAlgos
 *
 *  @brief ImgVsTimeSplitInFiles is a test/example module for psana framework.
 *
 *  ImgVsTimeSplitInFiles psana module class works after CSPadImageProducer.
 *  It gets the Image2D object from the event.
 *  This image object may be used in data processing.
 *  For the test purpose, the image of particular event is saved in the text file.
 *  This event number is defined in the psana.cfg configuration file. 
 *
 *  This software was developed for the LCLS project.  If you use all or 
 *  part of it, please give an appropriate acknowledgment.
 *
 *  @see CSPadImageProducer
 *
 *  @version \$Id$
 *
 *  @author Mikhail S. Dubrovin
 */

class ImgVsTimeSplitInFiles : public Module {
public:

  enum FILE_MODE {BINARY, TEXT};

  // Default constructor
  ImgVsTimeSplitInFiles (const std::string& name) ;

  // Destructor
  virtual ~ImgVsTimeSplitInFiles () ;

  /// Method which is called once at the beginning of the job
  virtual void beginJob(Event& evt, Env& env);
  
  /// Method which is called at the beginning of the run
  virtual void beginRun(Event& evt, Env& env);
  
  /// Method which is called at the beginning of the calibration cycle
  virtual void beginCalibCycle(Event& evt, Env& env);
  
  /// Method which is called with event data, this is the only required 
  /// method, all other methods are optional
  virtual void event(Event& evt, Env& env);
  
  /// Method which is called at the end of the calibration cycle
  virtual void endCalibCycle(Event& evt, Env& env);

  /// Method which is called at the end of the run
  virtual void endRun(Event& evt, Env& env);

  /// Method which is called once at the end of the job
  virtual void endJob(Event& evt, Env& env);

protected:

  void setFileMode();
  void initSplitInFiles(Event& evt, Env& env);
  void saveImageInFile(Event& evt);
  void printInputParameters();
  void printEventRecord(Event& evt, std::string comment=std::string());
  void printSummary(Event& evt, std::string comment=std::string());
  void openOutputFiles(Event& evt);
  void closeOutputFiles();
  void saveMetadataInFile();
  void procEvent(Event& evt);

private:

  // Data members, this is for example purposes only

  Pds::Src    m_src;
  std::string m_str_src;      // i.e. CxiDs1.0:Cspad.0
  std::string m_key;          // i.e. Image2D
  std::string m_fname_prefix; // prefix of the file name
  std::string m_file_type;    // file type "txt" or "bin" 
  unsigned    m_nfiles_out;
  double      m_ampl_thr;
  double      m_ampl_min;
  unsigned    m_print_bits;
  long        m_count;

  FILE_MODE   m_file_mode;

  unsigned    m_img_rows;
  unsigned    m_img_cols;
  unsigned    m_img_size;
  unsigned    m_blk_size;
  unsigned    m_rst_size;

  unsigned*   m_data;         // image data in case if processing is necessary

  std::string m_fname_common;
  std::ofstream* p_out;

protected:
//--------------------
// Splits the image for blocks and saves the blocks in files
    template <typename T>
    void procImgData(const boost::shared_ptr< ndarray<T,2> >& p_nda)
    {
      const T* data = p_nda->data();
      for(unsigned i=0; i<m_img_size; i++) {m_data[i] = (data[i] > m_ampl_thr) ? (unsigned)data[i] : (unsigned)m_ampl_min;}
    }

//--------------------
// Splits the image for blocks and saves the blocks in files
    template <typename T>
    void splitAndWriteImgInFiles (const boost::shared_ptr< ndarray<T,2> >& p_ndarr, 
                                  bool print_msg=false) 
    {
      const T* img_data = p_ndarr->data();               // Access to entire image
      for(unsigned b=0; b<m_nfiles_out; b++){

	const T* p_block_data = &img_data[b*m_blk_size]; // Access to the block 

	if (m_file_mode == TEXT) {
	  std::stringstream ss; 
	  for(unsigned i=0; i<m_blk_size; i++) ss << *p_block_data++ << " ";
	  std::string s = ss.str(); 
	  p_out[b].write(s.c_str(), s.size());
	} 

        else if (m_file_mode == BINARY) {
	  p_out[b].write(reinterpret_cast<const char*>(p_block_data), m_blk_size*sizeof(T));
	} 

        else {
          p_out[b] << " UNKNOWN FILE TYPE:" << m_file_type << " AND MODE:" << m_file_mode;
	}

	p_out[b] <<  "\n";
      }
    }

//--------------------
};

} // namespace ImgAlgos

#endif // IMGALGOS_IMGVSTIMESPLITINFILES_H