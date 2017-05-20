


#include <boost/python.hpp>
#include "ndarray/ndarray.h"
//#include <cstddef>  // for size_t

#include "ImgAlgos/AlgArrProc.h"
#include "ImgAlgos/AlgImgProc.h"
#include "ImgAlgos/AlgUtils.h"

//-------------------
using namespace ImgAlgos;

typedef ImgAlgos::AlgArrProc::mask_t mask_t;      // uint16_t
typedef ImgAlgos::AlgArrProc::wind_t wind_t;      // uint32_t
typedef ImgAlgos::AlgImgProc::conmap_t conmap_t;  // uint32_t
typedef ImgAlgos::AlgImgProc::pixel_status_t pixel_status_t;     //uint16_t
typedef ImgAlgos::AlgImgProc::pixel_maximums_t pixel_maximums_t; //uint16_t
typedef ImgAlgos::AlgImgProc::pixel_minimums_t pixel_minimums_t; //uint16_t

//-------------------

void     (AlgArrProc::*p_swin1) (ndarray<const wind_t,2>) = &AlgArrProc::setWindows;
void     (AlgArrProc::*p_sson1) (const float&, const float&) = &AlgArrProc::setSoNPars;
void     (AlgArrProc::*p_spsp1) (const float&, const float&, const float&, const float&, const float&) = &AlgArrProc::setPeakSelectionPars;

unsigned (AlgArrProc::*p_npix_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const float&)    = &AlgArrProc::numberOfPixAboveThr<float,   2>;
unsigned (AlgArrProc::*p_npix_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const double&)   = &AlgArrProc::numberOfPixAboveThr<double,  2>;
unsigned (AlgArrProc::*p_npix_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const int&)      = &AlgArrProc::numberOfPixAboveThr<int,     2>;
unsigned (AlgArrProc::*p_npix_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const int16_t&)  = &AlgArrProc::numberOfPixAboveThr<int16_t, 2>;
unsigned (AlgArrProc::*p_npix_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const uint16_t&) = &AlgArrProc::numberOfPixAboveThr<uint16_t,2>;

unsigned (AlgArrProc::*p_npix_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const float&)    = &AlgArrProc::numberOfPixAboveThr<float,   3>;
unsigned (AlgArrProc::*p_npix_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const double&)   = &AlgArrProc::numberOfPixAboveThr<double,  3>;
unsigned (AlgArrProc::*p_npix_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const int&)      = &AlgArrProc::numberOfPixAboveThr<int,     3>;
unsigned (AlgArrProc::*p_npix_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const int16_t&)  = &AlgArrProc::numberOfPixAboveThr<int16_t, 3>;
unsigned (AlgArrProc::*p_npix_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const uint16_t&) = &AlgArrProc::numberOfPixAboveThr<uint16_t,3>;

//-------------------

double (AlgArrProc::*p_ipix_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const float&)    = &AlgArrProc::intensityOfPixAboveThr<float,   2>;
double (AlgArrProc::*p_ipix_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const double&)   = &AlgArrProc::intensityOfPixAboveThr<double,  2>;
double (AlgArrProc::*p_ipix_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const int&)      = &AlgArrProc::intensityOfPixAboveThr<int,     2>;
double (AlgArrProc::*p_ipix_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const int16_t&)  = &AlgArrProc::intensityOfPixAboveThr<int16_t, 2>;
double (AlgArrProc::*p_ipix_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const uint16_t&) = &AlgArrProc::intensityOfPixAboveThr<uint16_t,2>;

double (AlgArrProc::*p_ipix_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const float&)    = &AlgArrProc::intensityOfPixAboveThr<float,   3>;
double (AlgArrProc::*p_ipix_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const double&)   = &AlgArrProc::intensityOfPixAboveThr<double,  3>;
double (AlgArrProc::*p_ipix_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const int&)      = &AlgArrProc::intensityOfPixAboveThr<int,     3>;
double (AlgArrProc::*p_ipix_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const int16_t&)  = &AlgArrProc::intensityOfPixAboveThr<int16_t, 3>;
double (AlgArrProc::*p_ipix_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const uint16_t&) = &AlgArrProc::intensityOfPixAboveThr<uint16_t,3>;

//-------------------

ndarray<const float, 2> (AlgArrProc::*p_pfv01_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const float&,    const float&,    const unsigned&, const float&) = &AlgArrProc::peakFinderV1<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv01_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const double&,   const double&,   const unsigned&, const float&) = &AlgArrProc::peakFinderV1<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv01_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const int&,      const int&,      const unsigned&, const float&) = &AlgArrProc::peakFinderV1<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv01_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const int16_t&,  const int16_t&,  const unsigned&, const float&) = &AlgArrProc::peakFinderV1<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv01_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const uint16_t&, const uint16_t&, const unsigned&, const float&) = &AlgArrProc::peakFinderV1<uint16_t,2>;
																			              
ndarray<const float, 2> (AlgArrProc::*p_pfv01_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const float&,    const float&,    const unsigned&, const float&) = &AlgArrProc::peakFinderV1<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv01_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const double&,   const double&,   const unsigned&, const float&) = &AlgArrProc::peakFinderV1<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv01_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const int&,      const int&,      const unsigned&, const float&) = &AlgArrProc::peakFinderV1<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv01_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const int16_t&,  const int16_t&,  const unsigned&, const float&) = &AlgArrProc::peakFinderV1<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv01_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const uint16_t&, const uint16_t&, const unsigned&, const float&) = &AlgArrProc::peakFinderV1<uint16_t,3>;

//-------------------

ndarray<const float, 2> (AlgArrProc::*p_pfv02_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const float&,    const float&, const float&) = &AlgArrProc::peakFinderV2<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const double&,   const float&, const float&) = &AlgArrProc::peakFinderV2<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const int&,      const float&, const float&) = &AlgArrProc::peakFinderV2<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const int16_t&,  const float&, const float&) = &AlgArrProc::peakFinderV2<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const uint16_t&, const float&, const float&) = &AlgArrProc::peakFinderV2<uint16_t,2>;
																		  
ndarray<const float, 2> (AlgArrProc::*p_pfv02_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const float&,    const float&, const float&) = &AlgArrProc::peakFinderV2<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const double&,   const float&, const float&) = &AlgArrProc::peakFinderV2<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const int&,      const float&, const float&) = &AlgArrProc::peakFinderV2<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const int16_t&,  const float&, const float&) = &AlgArrProc::peakFinderV2<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const uint16_t&, const float&, const float&) = &AlgArrProc::peakFinderV2<uint16_t,3>;

ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const float&,    const float&, const float&) = &AlgArrProc::peakFinderV2r1<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const double&,   const float&, const float&) = &AlgArrProc::peakFinderV2r1<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const int&,      const float&, const float&) = &AlgArrProc::peakFinderV2r1<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const int16_t&,  const float&, const float&) = &AlgArrProc::peakFinderV2r1<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const uint16_t&, const float&, const float&) = &AlgArrProc::peakFinderV2r1<uint16_t,2>;
					       													  				   
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const float&,    const float&, const float&) = &AlgArrProc::peakFinderV2r1<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const double&,   const float&, const float&) = &AlgArrProc::peakFinderV2r1<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const int&,      const float&, const float&) = &AlgArrProc::peakFinderV2r1<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const int16_t&,  const float&, const float&) = &AlgArrProc::peakFinderV2r1<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv02r1_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const uint16_t&, const float&, const float&) = &AlgArrProc::peakFinderV2r1<uint16_t,3>;

ndarray<const pixel_status_t, 3> (AlgArrProc::*p_get_pstat) () = &AlgArrProc::mapsOfPixelStatus;
ndarray<const conmap_t, 3>       (AlgArrProc::*p_get_pfv02) () = &AlgArrProc::mapsOfConnectedPixels;

//-------------------

ndarray<const float, 2> (AlgArrProc::*p_pfv03_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<uint16_t,2>;
														      			  
ndarray<const float, 2> (AlgArrProc::*p_pfv03_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&) = &AlgArrProc::peakFinderV3<uint16_t,3>;

ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<uint16_t,2>;
					       													                      				   
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r1_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r1<uint16_t,3>;

ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<uint16_t,2>;
					       													                      			       	   
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv03r2_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const size_t&, const float&, const float&, const float&) = &AlgArrProc::peakFinderV3r2<uint16_t,3>;

ndarray<const pixel_minimums_t, 3> (AlgArrProc::*p_get_pfv03_min) () = &AlgArrProc::mapsOfLocalMinimums;
ndarray<const pixel_maximums_t, 3> (AlgArrProc::*p_get_pfv03_max) () = &AlgArrProc::mapsOfLocalMaximums;

//-------------------

ndarray<const float, 2> (AlgArrProc::*p_pfv04_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const float&,    const float&,    const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const double&,   const double&,   const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const int&,      const int&,      const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const int16_t&,  const int16_t&,  const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const uint16_t&, const uint16_t&, const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<uint16_t,2>;
																			              
ndarray<const float, 2> (AlgArrProc::*p_pfv04_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const float&,    const float&,    const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const double&,   const double&,   const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const int&,      const int&,      const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const int16_t&,  const int16_t&,  const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const uint16_t&, const uint16_t&, const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4<uint16_t,3>;

ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const float&,    const float&,    const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const double&,   const double&,   const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const int&,      const int&,      const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const int16_t&,  const int16_t&,  const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const uint16_t&, const uint16_t&, const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<uint16_t,2>;
					       														              						     
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const float&,    const float&,    const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const double&,   const double&,   const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const int&,      const int&,      const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const int16_t&,  const int16_t&,  const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r1_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const uint16_t&, const uint16_t&, const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r1<uint16_t,3>;

ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_f2) (ndarray<const float,   2>, ndarray<const mask_t,2>, const float&,    const float&,    const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<float,   2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_d2) (ndarray<const double,  2>, ndarray<const mask_t,2>, const double&,   const double&,   const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<double,  2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_i2) (ndarray<const int,     2>, ndarray<const mask_t,2>, const int&,      const int&,      const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<int,     2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_s2) (ndarray<const int16_t, 2>, ndarray<const mask_t,2>, const int16_t&,  const int16_t&,  const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<int16_t, 2>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_u2) (ndarray<const uint16_t,2>, ndarray<const mask_t,2>, const uint16_t&, const uint16_t&, const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<uint16_t,2>;
					       														              						     
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_f3) (ndarray<const float,   3>, ndarray<const mask_t,3>, const float&,    const float&,    const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<float,   3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_d3) (ndarray<const double,  3>, ndarray<const mask_t,3>, const double&,   const double&,   const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<double,  3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_i3) (ndarray<const int,     3>, ndarray<const mask_t,3>, const int&,      const int&,      const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<int,     3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_s3) (ndarray<const int16_t, 3>, ndarray<const mask_t,3>, const int16_t&,  const int16_t&,  const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<int16_t, 3>;
ndarray<const float, 2> (AlgArrProc::*p_pfv04r2_u3) (ndarray<const uint16_t,3>, ndarray<const mask_t,3>, const uint16_t&, const uint16_t&, const unsigned&, const float&, const float&) = &AlgArrProc::peakFinderV4r2<uint16_t,3>;

//-------------------
//void (AlgArrProc::*p_set01) (const float&, const float&) = &AlgArrProc::setSoNPars;
//void (ImgAlgos::AlgArrProc::*print_1) () = &ImgAlgos::AlgArrProc::printInputPars;
void (*p_print_1) () = &test_print_1;

ndarray<const AlgImgProc::nphoton_t, 2> (*p_photnums_v01_f2) (ndarray<const float,  2>, ndarray<const AlgImgProc::u8mask_t,2>, const float&) = &mapOfPhotonNumbersV1<float>;
ndarray<const AlgImgProc::nphoton_t, 2> (*p_photnums_v01_d2) (ndarray<const double, 2>, ndarray<const AlgImgProc::u8mask_t,2>, const float&) = &mapOfPhotonNumbersV1<double>;

ndarray<const pixel_minimums_t, 2> (*p_loc_min_f2) (ndarray<const float,   2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMinimums<float>;
ndarray<const pixel_minimums_t, 2> (*p_loc_min_d2) (ndarray<const double,  2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMinimums<double>;
ndarray<const pixel_minimums_t, 2> (*p_loc_min_i2) (ndarray<const int,     2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMinimums<int>;
ndarray<const pixel_minimums_t, 2> (*p_loc_min_s2) (ndarray<const int16_t, 2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMinimums<int16_t>;
ndarray<const pixel_minimums_t, 2> (*p_loc_min_u2) (ndarray<const uint16_t,2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMinimums<uint16_t>;

ndarray<const pixel_maximums_t, 2> (*p_loc_max_f2) (ndarray<const float,   2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMaximums<float>;
ndarray<const pixel_maximums_t, 2> (*p_loc_max_d2) (ndarray<const double,  2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMaximums<double>;
ndarray<const pixel_maximums_t, 2> (*p_loc_max_i2) (ndarray<const int,     2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMaximums<int>;
ndarray<const pixel_maximums_t, 2> (*p_loc_max_s2) (ndarray<const int16_t, 2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMaximums<int16_t>;
ndarray<const pixel_maximums_t, 2> (*p_loc_max_u2) (ndarray<const uint16_t,2>, ndarray<const AlgImgProc::mask_t,2>, const size_t&) = &mapOfLocalMaximums<uint16_t>;

ndarray<const pixel_maximums_t, 2> (*p_loc_max_cross1) (ndarray<const float, 2>) = &mapOfLocalMaximumsRank1Cross;

//-------------------

// BOOST wrapper to create imgalgos_ext module that contains the ImgAlgos::AlgArrProc
// python class that calls the C++ ImgAlgos::AlgArrProc methods
// NB: The name of the python module (imgalgos_ext) MUST match the name given in
// PYEXTMOD in the SConscript

BOOST_PYTHON_MODULE(imgalgos_ext)
{    
  using namespace boost::python;

  //boost::python::class_<AlgArrProc>("AlgArrProc", init<const unsigned&>())
  boost::python::class_<AlgArrProc>("AlgArrProc", init<ndarray<const wind_t,2>, const unsigned&>())
 
    .def("set_windows", p_swin1)     
    .def("set_peak_selection_pars", p_spsp1)     
    .def("set_son_pars", p_sson1)     
    .def("print_input_pars", &AlgArrProc::printInputPars)

    //.def("set_son_parameters", p_set01)     

    .def("number_of_pix_above_thr_f2", p_npix_f2)
    .def("number_of_pix_above_thr_d2", p_npix_d2)
    .def("number_of_pix_above_thr_i2", p_npix_i2)
    .def("number_of_pix_above_thr_s2", p_npix_s2)
    .def("number_of_pix_above_thr_u2", p_npix_u2)
    				    	  	    
    .def("number_of_pix_above_thr_f3", p_npix_f3)
    .def("number_of_pix_above_thr_d3", p_npix_d3)
    .def("number_of_pix_above_thr_i3", p_npix_i3)
    .def("number_of_pix_above_thr_s3", p_npix_s3)
    .def("number_of_pix_above_thr_u3", p_npix_u3)


    .def("intensity_of_pix_above_thr_f2", p_ipix_f2)
    .def("intensity_of_pix_above_thr_d2", p_ipix_d2)
    .def("intensity_of_pix_above_thr_i2", p_ipix_i2)
    .def("intensity_of_pix_above_thr_s2", p_ipix_s2)
    .def("intensity_of_pix_above_thr_u2", p_ipix_u2)
    	  		    	  	    
    .def("intensity_of_pix_above_thr_f3", p_ipix_f3)
    .def("intensity_of_pix_above_thr_d3", p_ipix_d3)
    .def("intensity_of_pix_above_thr_i3", p_ipix_i3)
    .def("intensity_of_pix_above_thr_s3", p_ipix_s3)
    .def("intensity_of_pix_above_thr_u3", p_ipix_u3)


    .def("peak_finder_v1_f2", p_pfv01_f2)
    .def("peak_finder_v1_d2", p_pfv01_d2)
    .def("peak_finder_v1_i2", p_pfv01_i2)
    .def("peak_finder_v1_s2", p_pfv01_s2)
    .def("peak_finder_v1_u2", p_pfv01_u2)
    	  
    .def("peak_finder_v1_f3", p_pfv01_f3)
    .def("peak_finder_v1_d3", p_pfv01_d3)
    .def("peak_finder_v1_i3", p_pfv01_i3)
    .def("peak_finder_v1_s3", p_pfv01_s3)
    .def("peak_finder_v1_u3", p_pfv01_u3)


    //.def("peak_finder_v2_f2", p_pfv02_f2)
    //.def("peak_finder_v2_d2", p_pfv02_d2)
    //.def("peak_finder_v2_i2", p_pfv02_i2)
    //.def("peak_finder_v2_s2", p_pfv02_s2)
    //.def("peak_finder_v2_u2", p_pfv02_u2)
        	   
    //.def("peak_finder_v2_f3", p_pfv02_f3)
    //.def("peak_finder_v2_d3", p_pfv02_d3)
    //.def("peak_finder_v2_i3", p_pfv02_i3)
    //.def("peak_finder_v2_s3", p_pfv02_s3)
    //.def("peak_finder_v2_u3", p_pfv02_u3)

    .def("peak_finder_v2_f2", p_pfv02r1_f2)
    .def("peak_finder_v2_d2", p_pfv02r1_d2)
    .def("peak_finder_v2_i2", p_pfv02r1_i2)
    .def("peak_finder_v2_s2", p_pfv02r1_s2)
    .def("peak_finder_v2_u2", p_pfv02r1_u2)
    	   			       
    .def("peak_finder_v2_f3", p_pfv02r1_f3)
    .def("peak_finder_v2_d3", p_pfv02r1_d3)
    .def("peak_finder_v2_i3", p_pfv02r1_i3)
    .def("peak_finder_v2_s3", p_pfv02r1_s3)
    .def("peak_finder_v2_u3", p_pfv02r1_u3)

    .def("peak_finder_v2r1_f2", p_pfv02r1_f2)
    .def("peak_finder_v2r1_d2", p_pfv02r1_d2)
    .def("peak_finder_v2r1_i2", p_pfv02r1_i2)
    .def("peak_finder_v2r1_s2", p_pfv02r1_s2)
    .def("peak_finder_v2r1_u2", p_pfv02r1_u2)
    	   			         
    .def("peak_finder_v2r1_f3", p_pfv02r1_f3)
    .def("peak_finder_v2r1_d3", p_pfv02r1_d3)
    .def("peak_finder_v2r1_i3", p_pfv02r1_i3)
    .def("peak_finder_v2r1_s3", p_pfv02r1_s3)
    .def("peak_finder_v2r1_u3", p_pfv02r1_u3)

    .def("maps_of_pixel_status", p_get_pstat)
    .def("maps_of_connected_pixels", p_get_pfv02)


    //.def("peak_finder_v3_f2", p_pfv03_f2)
    //.def("peak_finder_v3_d2", p_pfv03_d2)
    //.def("peak_finder_v3_i2", p_pfv03_i2)
    //.def("peak_finder_v3_s2", p_pfv03_s2)
    //.def("peak_finder_v3_u2", p_pfv03_u2)
        	   
    //.def("peak_finder_v3_f3", p_pfv03_f3)
    //.def("peak_finder_v3_d3", p_pfv03_d3)
    //.def("peak_finder_v3_i3", p_pfv03_i3)
    //.def("peak_finder_v3_s3", p_pfv03_s3)
    //.def("peak_finder_v3_u3", p_pfv03_u3)

    .def("peak_finder_v3_f2", p_pfv03r1_f2)
    .def("peak_finder_v3_d2", p_pfv03r1_d2)
    .def("peak_finder_v3_i2", p_pfv03r1_i2)
    .def("peak_finder_v3_s2", p_pfv03r1_s2)
    .def("peak_finder_v3_u2", p_pfv03r1_u2)
    	   			       
    .def("peak_finder_v3_f3", p_pfv03r1_f3)
    .def("peak_finder_v3_d3", p_pfv03r1_d3)
    .def("peak_finder_v3_i3", p_pfv03r1_i3)
    .def("peak_finder_v3_s3", p_pfv03r1_s3)
    .def("peak_finder_v3_u3", p_pfv03r1_u3)

    .def("peak_finder_v3r1_f2", p_pfv03r1_f2)
    .def("peak_finder_v3r1_d2", p_pfv03r1_d2)
    .def("peak_finder_v3r1_i2", p_pfv03r1_i2)
    .def("peak_finder_v3r1_s2", p_pfv03r1_s2)
    .def("peak_finder_v3r1_u2", p_pfv03r1_u2)
    	   		  	       
    .def("peak_finder_v3r1_f3", p_pfv03r1_f3)
    .def("peak_finder_v3r1_d3", p_pfv03r1_d3)
    .def("peak_finder_v3r1_i3", p_pfv03r1_i3)
    .def("peak_finder_v3r1_s3", p_pfv03r1_s3)
    .def("peak_finder_v3r1_u3", p_pfv03r1_u3)

    .def("peak_finder_v3r2_f2", p_pfv03r2_f2)
    .def("peak_finder_v3r2_d2", p_pfv03r2_d2)
    .def("peak_finder_v3r2_i2", p_pfv03r2_i2)
    .def("peak_finder_v3r2_s2", p_pfv03r2_s2)
    .def("peak_finder_v3r2_u2", p_pfv03r2_u2)
    	   		  	       	 
    .def("peak_finder_v3r2_f3", p_pfv03r2_f3)
    .def("peak_finder_v3r2_d3", p_pfv03r2_d3)
    .def("peak_finder_v3r2_i3", p_pfv03r2_i3)
    .def("peak_finder_v3r2_s3", p_pfv03r2_s3)
    .def("peak_finder_v3r2_u3", p_pfv03r2_u3)

    .def("maps_of_local_minimums", p_get_pfv03_min)
    .def("maps_of_local_maximums", p_get_pfv03_max)


    //.def("peak_finder_v4_f2", p_pfv04_f2)
    //.def("peak_finder_v4_d2", p_pfv04_d2)
    //.def("peak_finder_v4_i2", p_pfv04_i2)
    //.def("peak_finder_v4_s2", p_pfv04_s2)
    //.def("peak_finder_v4_u2", p_pfv04_u2)
    	  			     
    //.def("peak_finder_v4_f3", p_pfv04_f3)
    //.def("peak_finder_v4_d3", p_pfv04_d3)
    //.def("peak_finder_v4_i3", p_pfv04_i3)
    //.def("peak_finder_v4_s3", p_pfv04_s3)
    //.def("peak_finder_v4_u3", p_pfv04_u3)

    .def("peak_finder_v4_f2", p_pfv04r1_f2)
    .def("peak_finder_v4_d2", p_pfv04r1_d2)
    .def("peak_finder_v4_i2", p_pfv04r1_i2)
    .def("peak_finder_v4_s2", p_pfv04r1_s2)
    .def("peak_finder_v4_u2", p_pfv04r1_u2)
    	  			       
    .def("peak_finder_v4_f3", p_pfv04r1_f3)
    .def("peak_finder_v4_d3", p_pfv04r1_d3)
    .def("peak_finder_v4_i3", p_pfv04r1_i3)
    .def("peak_finder_v4_s3", p_pfv04r1_s3)
    .def("peak_finder_v4_u3", p_pfv04r1_u3)

    .def("peak_finder_v4r1_f2", p_pfv04r1_f2)
    .def("peak_finder_v4r1_d2", p_pfv04r1_d2)
    .def("peak_finder_v4r1_i2", p_pfv04r1_i2)
    .def("peak_finder_v4r1_s2", p_pfv04r1_s2)
    .def("peak_finder_v4r1_u2", p_pfv04r1_u2)
    	  		  	         
    .def("peak_finder_v4r1_f3", p_pfv04r1_f3)
    .def("peak_finder_v4r1_d3", p_pfv04r1_d3)
    .def("peak_finder_v4r1_i3", p_pfv04r1_i3)
    .def("peak_finder_v4r1_s3", p_pfv04r1_s3)
    .def("peak_finder_v4r1_u3", p_pfv04r1_u3)

    .def("peak_finder_v4r2_f2", p_pfv04r2_f2)
    .def("peak_finder_v4r2_d2", p_pfv04r2_d2)
    .def("peak_finder_v4r2_i2", p_pfv04r2_i2)
    .def("peak_finder_v4r2_s2", p_pfv04r2_s2)
    .def("peak_finder_v4r2_u2", p_pfv04r2_u2)
    	  		  	         
    .def("peak_finder_v4r2_f3", p_pfv04r2_f3)
    .def("peak_finder_v4r2_d3", p_pfv04r2_d3)
    .def("peak_finder_v4r2_i3", p_pfv04r2_i3)
    .def("peak_finder_v4r2_s3", p_pfv04r2_s3)
    .def("peak_finder_v4r2_u3", p_pfv04r2_u3)
  ;

  // Global methods from AlgUtils
  def("test_print_1", p_print_1);

  // Global methods from AlgImgProc
  def("map_photon_numbers_v1_f2", p_photnums_v01_f2);
  def("map_photon_numbers_v1_d2", p_photnums_v01_d2);

  def("local_minimums_f2", p_loc_min_f2);
  def("local_minimums_d2", p_loc_min_d2);
  def("local_minimums_i2", p_loc_min_i2);
  def("local_minimums_s2", p_loc_min_s2);
  def("local_minimums_u2", p_loc_min_u2);

  def("local_maximums_f2", p_loc_max_f2);
  def("local_maximums_d2", p_loc_max_d2);
  def("local_maximums_i2", p_loc_max_i2);
  def("local_maximums_s2", p_loc_max_s2);
  def("local_maximums_u2", p_loc_max_u2);

  def("local_maximums_cross1", p_loc_max_cross1);

}

//-------------------

