
# File and Version Information:
#  $Id$
#
# Description:
#  class PyAlgos
#

"""Class provides access to C++ algorithms from python.

Usage::

    # !!! None is returned whenever requested information is missing.

    # IMPORT
    # ======
    import psana
    from ImgAlgos.PyAlgos import PyAlgos


    # INPUT PARAMETERS
    # ================
    # List of windows
    winds = None # entire size of all segments will be used for peak finding
    winds = (( 0, 0, 185, 0, 388), ( 1, 20, 160, 30, 300), ( 7, 0, 185, 0, 388))

    # Mask
    mask = None                   # (default) all pixels in windows will be used for peak finding
    mask = det.mask()             # see class Detector.PyDetector
    mask = np.loadtxt(fname_mask) #
    mask.shape = (32,185,388)     # should be the same as shape of data n-d array

    # Data n-d array
    nda = det.calib() # see class Detector.PyDetector


    # INITIALIZATION
    # ==============
    # create object:
    alg = PyAlgos(windows=winds, mask=mask, pbits=0)
    # where pbits - is a print info control bit-word:
    # pbits = 0   - print nothing
    #       + 1   - main results, list of peaks
    #       + 2   - input parameters, index matrix of pixels for S/N algorithm
    #       + 128 - tracking and all details in class PyAlgos.py
    #       + 256 - tracking and all details in class AlgArrProc
    #       + 512 - tracking and all details in class AlgImgProc

    # set peak-selector parameters:
    alg.set_peak_selection_pars(npix_min=5, npix_max=5000, amax_thr=0, atot_thr=0, son_min=10)


    # HIT FINDERS
    # ===========
    # Hit finders return simple values for decision on event selection.

    # get number of pixels above threshold
    npix = alg.number_of_pix_above_thr(data, thr=10)

    # get total intensity of pixels above threshold
    intensity = alg.intensity_of_pix_above_thr(data, thr=12)


    # PEAK FINDERS
    # ============
    # Peak finders return list (numpy.array) of records with found peak parameters.

    # v1 - aka Droplet Finder - two-threshold peak-finding algorithm in restricted region
    #                           around pixel with maximal intensity.
    peaks = alg.peak_finder_v1(nda, thr_low=10, thr_high=150, radius=5, dr=0.05)

    # v2 - define peaks for regoins of connected pixels above threshold
    peaks = alg.peak_finder_v2(nda, thr=10, r0=7.0, dr=2.0)

    # v3 - define peaks in local maximums of specified rank (radius),
    #      for example rank=2 means 5x5 pixel region around central pixel.
    #      nsigm - optional threshold (number of noise rms) for pixel selection 
    peaks = alg.peak_finder_v3(nda, rank=2, r0=7.0, dr=2.0, nsigm=0)
    peaks = alg.peak_finder_v3r2(nda, rank=2, r0=7.0, dr=2.0, nsigm=0) # test version for connected pixels

    # v4 - aka Droplet Finder - the same as v1, but uses rank and r0 parameters in stead of common radius.
    peaks = alg.peak_finder_v4(nda, thr_low=10, thr_high=150, rank=4, r0=7.0, dr=2.0)
    peaks = alg.peak_finder_v4r2(nda, thr_low=10, thr_high=150, rank=5, r0=7.0, dr=2.0) # test version for connected pixels


    # OPTIONAL METHODS
    # ================
    # print info
    alg.print_attributes()   # attributes of the PyAlgos object
    alg.print_input_pars()   # member data of C++ objects

    # set parameters for S/N evaluation algorithm
    alg.set_son_pars(r0=5, dr=0.05)

    # set mask
    alg.set_mask(mask)

    # set windows in segments to search for peaks
    alg.set_windows(winds)

    # Call after alg.peak_finder_v2, v4r2 ONLY! Returns n-d array with 2-d maps of pixel status
    maps = alg.maps_of_pixel_status()

    # Call after alg.peak_finder_v2 ONLY! Returns n-d array with 2-d maps of connected pixels
    maps = alg.maps_of_connected_pixels()

    # Call after alg.peak_finder_v3 ONLY! Returns n-d array with 2-d maps of local maximums
    maps = alg.maps_of_local_maximums()

    # GLOBAL METHODS
    # ==============
    #   Subtracts numpy array of bkgd from data using normalization in windows.
    #   Each window is specified by 5 parameters: (segment, rowmin, rowmax, colmin, colmax)
    #   For 2-d arrays segment is not used, but still 5 parameters needs to be specified.
    cdata = subtract_bkgd(data, bkgd, mask=None, winds=None, pbits=0)

    # Merges photons split among pixels and returns n-d array with integer number of photons per pixel.
    nphotons_nda = photons(fphotons_nda, mask, thr_fraction)

@see classes
\n  :py:class:`Detector.AreaDetector` - acess to detector data

This software was developed for the LCLS project.
If you use all or part of it, please give an appropriate acknowledgment.

@version $Id$

@author Mikhail S. Dubrovin
"""
__version__ = "$Revision$"

import sys
import numpy as np

import ImgAlgos
import pyimgalgos.GlobalUtils as piagu


def reshape_nda_to_2d(arr) :
    """Reshape np.array to 2-d
    """
    sh = arr.shape
    if arr.ndim<3 : return arr
    arr.shape = (arr.size//sh[-1], sh[-1])
    return arr


def reshape_nda_to_3d(arr) :
    """Reshape np.array to 3-d
    """
    sh = arr.shape
    if arr.ndim<4 : return arr
    arr.shape = (arr.size//sh[-1]//sh[-2], sh[-2], sh[-1])
    return arr


def print_arr_attr(arr, cmt='') :
    if arr is None :
        print('  %s attributes: array is None' % (cmt))
        return
    nda = np.array(arr)
    print('  %s attributes: size = %d, shape = %s, dtype = %s' % (cmt, nda.size, nda.shape, nda.dtype))


def print_arr(arr, cmt='') :
    if arr is None :
        print('  %s attributes: array is None' % (cmt))
        return
    nda = np.array(arr)
    print('\nprint_arr: %s:\n%s' % (cmt, str(nda)))
    print('  %s attributes: size = %d, shape = %s, dtype = %s' % (cmt, nda.size, nda.shape, nda.dtype))


class PyAlgos(object) :
    """Python wrapper for C++ algorithms

    Low level algorithms are implemented on C++

    @see AlgArrProc - c++ array processing algorithms
    """


    def __init__(self, windows=None, mask=None, pbits=0) :
        """Constructor.
        - windows - tuple, list or numpy array of windows or None
        - mask  - n-d array with mask or None
        - pbits - print control bit-word
        """
        if pbits & 128 : print('in c-tor %s' % self.__class__.__name__)

        self.pbits = pbits

        self.set_mask(mask)

        if windows is None : self.windows = np.empty((0,0), dtype=np.uint32)
        else               : self.windows = np.array(windows, dtype=np.uint32)

        self.aap = ImgAlgos.AlgArrProc(self.windows, self.pbits)

        if self.pbits == 2 : self.print_attributes()


    def set_windows(self, windows) :
        """
        - windows - tuple of windows
        """
        if self.pbits & 128 : print('in PyAlgos.set_windows()')

        self.windows = np.array(windows, dtype=np.uint32)
        self.aap.set_windows(self.windows)


    def set_son_pars(self, r0=10, dr=0.05) :
        """ Set parameters for SoN (S/N) evaluation
        - r0 - ring internal radius
        - dr - ring width
        """
        if self.pbits & 128 : print('in PyAlgos.set_son_pars()')

        self.aap.set_son_pars(r0, dr)


    def set_peak_selection_pars(self, npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=0, son_min=0) :
        """
        - npix_min - minimal number of pixels in peak
        - npix_max - maximal number of pixels in peak
        - amax_thr - threshold on pixel amplitude
        - amax_thr - threshold on total amplitude
        - son_min - minimal S/N in peak
        """
        if self.pbits & 128 : print('in PyAlgos.set_peak_selection_pars()')

        self.aap.set_peak_selection_pars(npix_min, npix_max, amax_thr, atot_thr, son_min)


    def set_mask(self, mask) :
        """
        - mask - array with mask 1/0 - good/bad pixel
        """
        if self.pbits & 128 : print('in PyAlgos.set_mask()')

        if mask is None : self.mask = None
        else : self.mask = np.array(mask, dtype=np.uint16)


    def print_input_pars(self) :
        if self.pbits & 128 : print('in PyAlgos.print_input_pars()')

        self.aap.print_input_pars()


    def print_attributes(self) :
        print('%s attributes:' % self.__class__.__name__, \
              '\n  pbits  : %d' % self.pbits) # , \
              #'\n  windows:%s' % str(self.windows), \
              #'\n  size: %s  shape: %s  dtype: %s' % (self.windows.size, self.windows.shape, self.windows.dtype)
        print_arr_attr(self.windows, cmt='windows')
        print_arr_attr(self.mask, cmt='mask')


    def check_mask(self, ndim, dtype=np.uint16) :
        """Returns empty mask for None or re-shaped mask for ndim>3, or self.mask
        """
        if self.pbits & 128 : print_arr_attr(self.mask, cmt='PyAlgos.check_mask() self.mask')

        if self.mask is None :
            if ndim>2 : self.mask = np.empty((0,0,0), dtype=dtype)
            else      : self.mask = np.empty((0,0),   dtype=dtype)
            return

        if ndim>3 : self.mask = reshape_nda_to_3d(self.mask)


    def number_of_pix_above_thr(self, arr, thr=0) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.number_of_pix_above_thr input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.number_of_pix_above_thr_f2(nda, msk, thr)
            if dtype == np.float64: return self.aap.number_of_pix_above_thr_d2(nda, msk, thr)
            if dtype == np.int32  : return self.aap.number_of_pix_above_thr_i2(nda, msk, thr)
            if dtype == np.int16  : return self.aap.number_of_pix_above_thr_s2(nda, msk, thr)
            if dtype == np.uint16 : return self.aap.number_of_pix_above_thr_u2(nda, msk, thr)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.number_of_pix_above_thr_f3(nda, msk, thr)
        if dtype == np.float64: return self.aap.number_of_pix_above_thr_d3(nda, msk, thr)
        if dtype == np.int32  : return self.aap.number_of_pix_above_thr_i3(nda, msk, thr)
        if dtype == np.int16  : return self.aap.number_of_pix_above_thr_s3(nda, msk, thr)
        if dtype == np.uint16 : return self.aap.number_of_pix_above_thr_u3(nda, msk, thr)

        if self.pbits :
            print('WARNING: PyAlgos.number_of_pix_above_thr(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def intensity_of_pix_above_thr(self, arr, thr):

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.intensity_of_pix_above_thr() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.intensity_of_pix_above_thr_f2(nda, msk, thr)
            if dtype == np.float64: return self.aap.intensity_of_pix_above_thr_d2(nda, msk, thr)
            if dtype == np.int32  : return self.aap.intensity_of_pix_above_thr_i2(nda, msk, thr)
            if dtype == np.int16  : return self.aap.intensity_of_pix_above_thr_s2(nda, msk, thr)
            if dtype == np.uint16 : return self.aap.intensity_of_pix_above_thr_u2(nda, msk, thr)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.intensity_of_pix_above_thr_f3(nda, msk, thr)
        if dtype == np.float64: return self.aap.intensity_of_pix_above_thr_d3(nda, msk, thr)
        if dtype == np.int32  : return self.aap.intensity_of_pix_above_thr_i3(nda, msk, thr)
        if dtype == np.int16  : return self.aap.intensity_of_pix_above_thr_s3(nda, msk, thr)
        if dtype == np.uint16 : return self.aap.intensity_of_pix_above_thr_u3(nda, msk, thr)

        if self.pbits :
            print('WARNING: PyAlgos.intensity_of_pix_above_thr(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None



    def peak_finder_v1(self, arr, thr_low, thr_high, radius=5, dr=0.05) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v1() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v1_f2(nda, msk, thr_low, thr_high, radius, dr)
            if dtype == np.float64: return self.aap.peak_finder_v1_d2(nda, msk, thr_low, thr_high, radius, dr)
            if dtype == np.int32  : return self.aap.peak_finder_v1_i2(nda, msk, thr_low, thr_high, radius, dr)
            if dtype == np.int16  : return self.aap.peak_finder_v1_s2(nda, msk, thr_low, thr_high, radius, dr)
            if dtype == np.uint16 : return self.aap.peak_finder_v1_u2(nda, msk, thr_low, thr_high, radius, dr)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v1_f3(nda, msk, thr_low, thr_high, radius, dr)
        if dtype == np.float64: return self.aap.peak_finder_v1_d3(nda, msk, thr_low, thr_high, radius, dr)
        if dtype == np.int32  : return self.aap.peak_finder_v1_i3(nda, msk, thr_low, thr_high, radius, dr)
        if dtype == np.int16  : return self.aap.peak_finder_v1_s3(nda, msk, thr_low, thr_high, radius, dr)
        if dtype == np.uint16 : return self.aap.peak_finder_v1_u3(nda, msk, thr_low, thr_high, radius, dr)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v1(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def peak_finder_v4(self, arr, thr_low, thr_high, rank=4, r0=5.0, dr=0.05) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v4() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v4_f2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.float64: return self.aap.peak_finder_v4_d2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.int32  : return self.aap.peak_finder_v4_i2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.int16  : return self.aap.peak_finder_v4_s2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.uint16 : return self.aap.peak_finder_v4_u2(nda, msk, thr_low, thr_high, rank, r0, dr)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v4_f3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.float64: return self.aap.peak_finder_v4_d3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.int32  : return self.aap.peak_finder_v4_i3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.int16  : return self.aap.peak_finder_v4_s3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.uint16 : return self.aap.peak_finder_v4_u3(nda, msk, thr_low, thr_high, rank, r0, dr)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v4(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None



    def peak_finder_v4r1(self, arr, thr_low, thr_high, rank=4, r0=7.0, dr=2.0) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v4r1() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v4r1_f2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.float64: return self.aap.peak_finder_v4r1_d2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.int32  : return self.aap.peak_finder_v4r1_i2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.int16  : return self.aap.peak_finder_v4r1_s2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.uint16 : return self.aap.peak_finder_v4r1_u2(nda, msk, thr_low, thr_high, rank, r0, dr)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v4r1_f3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.float64: return self.aap.peak_finder_v4r1_d3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.int32  : return self.aap.peak_finder_v4r1_i3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.int16  : return self.aap.peak_finder_v4r1_s3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.uint16 : return self.aap.peak_finder_v4r1_u3(nda, msk, thr_low, thr_high, rank, r0, dr)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v4r1(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def peak_finder_v4r2(self, arr, thr_low, thr_high, rank=4, r0=7.0, dr=2.0) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v4r2() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v4r2_f2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.float64: return self.aap.peak_finder_v4r2_d2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.int32  : return self.aap.peak_finder_v4r2_i2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.int16  : return self.aap.peak_finder_v4r2_s2(nda, msk, thr_low, thr_high, rank, r0, dr)
            if dtype == np.uint16 : return self.aap.peak_finder_v4r2_u2(nda, msk, thr_low, thr_high, rank, r0, dr)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v4r2_f3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.float64: return self.aap.peak_finder_v4r2_d3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.int32  : return self.aap.peak_finder_v4r2_i3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.int16  : return self.aap.peak_finder_v4r2_s3(nda, msk, thr_low, thr_high, rank, r0, dr)
        if dtype == np.uint16 : return self.aap.peak_finder_v4r2_u3(nda, msk, thr_low, thr_high, rank, r0, dr)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v4r2(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def set_son_parameters(r0=5, dr=0.05) :
        self.aap.set_son_parameters(r0, dr)


    def peak_finder_v2(self, arr, thr=0, r0=5, dr=0.05) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v2() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v2_f2(nda, msk, thr, r0, dr)
            if dtype == np.float64: return self.aap.peak_finder_v2_d2(nda, msk, thr, r0, dr)
            if dtype == np.int32  : return self.aap.peak_finder_v2_i2(nda, msk, thr, r0, dr)
            if dtype == np.int16  : return self.aap.peak_finder_v2_s2(nda, msk, thr, r0, dr)
            if dtype == np.uint16 : return self.aap.peak_finder_v2_u2(nda, msk, thr, r0, dr)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v2_f3(nda, msk, thr, r0, dr)
        if dtype == np.float64: return self.aap.peak_finder_v2_d3(nda, msk, thr, r0, dr)
        if dtype == np.int32  : return self.aap.peak_finder_v2_i3(nda, msk, thr, r0, dr)
        if dtype == np.int16  : return self.aap.peak_finder_v2_s3(nda, msk, thr, r0, dr)
        if dtype == np.uint16 : return self.aap.peak_finder_v2_u3(nda, msk, thr, r0, dr)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v2(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def peak_finder_v2r1(self, arr, thr=0, r0=7.0, dr=2.0) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v2r1() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v2r1_f2(nda, msk, thr, r0, dr)
            if dtype == np.float64: return self.aap.peak_finder_v2r1_d2(nda, msk, thr, r0, dr)
            if dtype == np.int32  : return self.aap.peak_finder_v2r1_i2(nda, msk, thr, r0, dr)
            if dtype == np.int16  : return self.aap.peak_finder_v2r1_s2(nda, msk, thr, r0, dr)
            if dtype == np.uint16 : return self.aap.peak_finder_v2r1_u2(nda, msk, thr, r0, dr)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v2r1_f3(nda, msk, thr, r0, dr)
        if dtype == np.float64: return self.aap.peak_finder_v2r1_d3(nda, msk, thr, r0, dr)
        if dtype == np.int32  : return self.aap.peak_finder_v2r1_i3(nda, msk, thr, r0, dr)
        if dtype == np.int16  : return self.aap.peak_finder_v2r1_s3(nda, msk, thr, r0, dr)
        if dtype == np.uint16 : return self.aap.peak_finder_v2r1_u3(nda, msk, thr, r0, dr)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v2r1(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def maps_of_pixel_status(self) :

        if self.pbits & 128 : print('in PyAlgos.maps_of_pixel_status()')
        arr = self.aap.maps_of_pixel_status()
        if self.pbits & 128 : print_arr_attr(arr, cmt='maps_of_pixel_status arr:')
        return arr


    def maps_of_connected_pixels(self) :

        if self.pbits & 128 : print('in PyAlgos.maps_of_connected_pixels()')
        arr = self.aap.maps_of_connected_pixels()
        if self.pbits & 128 : print_arr_attr(arr, cmt='maps_of_connected_pixels arr:')
        return arr


    def peak_finder_v3(self, arr, rank=2, r0=5, dr=0.05, nsigm=0) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v3() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v3_f2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.float64: return self.aap.peak_finder_v3_d2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.int32  : return self.aap.peak_finder_v3_i2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.int16  : return self.aap.peak_finder_v3_s2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.uint16 : return self.aap.peak_finder_v3_u2(nda, msk, rank, r0, dr, nsigm)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v3_f3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.float64: return self.aap.peak_finder_v3_d3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.int32  : return self.aap.peak_finder_v3_i3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.int16  : return self.aap.peak_finder_v3_s3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.uint16 : return self.aap.peak_finder_v3_u3(nda, msk, rank, r0, dr, nsigm)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v3(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def peak_finder_v3r1(self, arr, rank=5, r0=7.0, dr=2.0, nsigm=0) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v3r1() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v3r1_f2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.float64: return self.aap.peak_finder_v3r1_d2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.int32  : return self.aap.peak_finder_v3r1_i2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.int16  : return self.aap.peak_finder_v3r1_s2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.uint16 : return self.aap.peak_finder_v3r1_u2(nda, msk, rank, r0, dr, nsigm)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v3r1_f3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.float64: return self.aap.peak_finder_v3r1_d3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.int32  : return self.aap.peak_finder_v3r1_i3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.int16  : return self.aap.peak_finder_v3r1_s3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.uint16 : return self.aap.peak_finder_v3r1_u3(nda, msk, rank, r0, dr, nsigm)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v3r1(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def peak_finder_v3r2(self, arr, rank=5, r0=7.0, dr=2.0, nsigm=0) :

        if self.pbits & 128 : print_arr_attr(arr, cmt='PyAlgos.peak_finder_v3r2() input arr:')

        ndim, dtype = arr.ndim, arr.dtype
        self.check_mask(ndim)
        nda, msk = arr, self.mask

        if ndim == 2 :
            if dtype == np.float32: return self.aap.peak_finder_v3r2_f2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.float64: return self.aap.peak_finder_v3r2_d2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.int32  : return self.aap.peak_finder_v3r2_i2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.int16  : return self.aap.peak_finder_v3r2_s2(nda, msk, rank, r0, dr, nsigm)
            if dtype == np.uint16 : return self.aap.peak_finder_v3r2_u2(nda, msk, rank, r0, dr, nsigm)

        if ndim>3 :
            nda = reshape_nda_to_3d(arr)

        if dtype == np.float32: return self.aap.peak_finder_v3r2_f3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.float64: return self.aap.peak_finder_v3r2_d3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.int32  : return self.aap.peak_finder_v3r2_i3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.int16  : return self.aap.peak_finder_v3r2_s3(nda, msk, rank, r0, dr, nsigm)
        if dtype == np.uint16 : return self.aap.peak_finder_v3r2_u3(nda, msk, rank, r0, dr, nsigm)

        if self.pbits :
            print('WARNING: PyAlgos.peak_finder_v3r2(.) method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))

        return None


    def maps_of_local_minimums(self) :

        if self.pbits & 128 : print('in PyAlgos.maps_of_local_minimums()')
        arr = self.aap.maps_of_local_minimums()
        if self.pbits & 128 : print_arr_attr(arr, cmt='maps_of_local_minimums arr:')
        return arr


    def maps_of_local_maximums(self) :

        if self.pbits & 128 : print('in PyAlgos.maps_of_local_maximums()')
        arr = self.aap.maps_of_local_maximums()
        if self.pbits & 128 : print_arr_attr(arr, cmt='maps_of_local_maximums arr:')
        return arr


def subtract_bkgd(data, bkgd, mask=None, winds=None, pbits=0) :
    """Subtracts numpy array of bkgd from data using normalization in windows.
       Each window is specified by 5 parameters: (segment, rowmin, rowmax, colmin, colmax)
       For 2-d arrays segment is not used, but still 5 parameters needs to be specified.
    """
    return piagu.subtract_bkgd(data, bkgd, mask, winds, pbits)


def photons_2d(data, mask=None, thr_fraction=0.9) :
    """returns 2-d array with number of merged photons per pixel
    """
    nda = data
    msk = mask
    thr = thr_fraction
    if mask is None :
        msk = np.ones(nda.shape, dtype=np.uint8)
    else :
        if msk.dtype != np.uint8  : raise ValueError('Mask dtype=%s, but expected np.uint8' % str(msk.dtype))
        if msk.shape != nda.shape : raise ValueError('msk.shape=%s is different from array shape=%s' % (str(msk.shape), str(nda.shape)))

    ndim, dtype = nda.ndim, nda.dtype

    if ndim != 2 : raise ValueError('photons_2d: number of dimensions is %d, expected 2-d' % ndim)

    if dtype == np.float32: return ImgAlgos.map_photon_numbers_v1_f2(nda, msk, thr)
    if dtype == np.float64: return ImgAlgos.map_photon_numbers_v1_d2(nda, msk, thr)
    #if dtype == np.int32  : return ImgAlgos.map_photon_numbers_v1_i2(nda, msk, thr)
    #if dtype == np.int16  : return ImgAlgos.map_photon_numbers_v1_s2(nda, msk, thr)
    #if dtype == np.uint16 : return ImgAlgos.map_photon_numbers_v1_u2(nda, msk, thr)

    print('WARNING: PyAlgos.photons_2d method is not implemented for ndim = %d, dtype = %s' % (ndim, str(dtype)))
    return None


def photons(data, mask, thr_fraction=0.9) :
    """returns 2-d or 3-d array with number of merged photons per pixel
    """
    ndim = data.ndim
    if ndim == 2 : return photons_2d(data, mask, thr_fraction)

    nda = data if ndim == 3 else reshape_nda_to_3d(data)

    if mask is None : return np.array([photons_2d(nda[s,:,:]) for s in range(nda.shape[0])], dtype=np.uint16)

    msk = mask if ndim == 3 else reshape_nda_to_3d(mask)

    # if mask.dtype != np.uint8 :
    msk = msk.astype(np.uint8, copy=False)

    return np.array([photons_2d(nda[s,:,:], msk[s,:,:], thr_fraction) for s in range(nda.shape[0])], dtype=np.uint16)

##---------- TEST -------------

def test_photons_2d() :
    from time import time

    mean, sigma = 0,10
    shape = (50,50)
    mask  = np.ones(shape, dtype=np.uint8)
    data  = np.array(mean + sigma * np.random.standard_normal(size=shape), dtype=np.float32)

    piagu.print_ndarr(data, 'data', last=10)
    piagu.print_ndarr(mask, 'mask', last=10)

    #print("data:", local_maximums_cross1(data))

    t0_sec = time()
    arr = photons_2d(data, mask, thr_fraction=0.9)
    print('\nTime consumed by photons_2d(data, mask) (sec) = %10.6f' % (time()-t0_sec))

    piagu.print_ndarr(arr, 'arr', last=10)

    from pyimgalgos.GlobalGraphics import plotImageLarge, show
    plotImageLarge(data, title='data')
    plotImageLarge(arr, title='array')
    show()


def test_photons_2d_chuck() :
    from time import time

    shape = (4,4)
    mask  = np.ones(shape, dtype=np.uint8)
    data  = np.array(((0.0, 3.5, 0.1, 0.2),
                      (0.2, 0.4, 0.0, 1.2),
                      (0.1, 4.7, 3.4, 0.0),
                      (0.5, 0.4, 0.4, 0.1)), dtype=np.float32)

    #shape = (6,6)
    #mask  = np.ones(shape, dtype=np.uint8)
    #data  = np.array(((0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    #                  (0.0, 0.0, 3.5, 0.1, 0.2, 0.0),
    #                  (0.0, 0.2, 0.4, 0.0, 1.2, 0.0),
    #                  (0.0, 0.1, 4.7, 3.4, 0.0, 0.0),
    #                  (0.0, 0.5, 0.4, 0.4, 0.1, 0.0),
    #                  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)), dtype=np.float32)

    piagu.print_ndarr(data, 'data', last=12)
    piagu.print_ndarr(mask, 'mask', last=12)

    print(("data:", data))

    t0_sec = time()
    arr = photons_2d(data, mask, thr_fraction=0.9)
    print('\nTime consumed by photons_2d(data, mask) (sec) = %10.6f' % (time()-t0_sec))
    print(("res:", arr))

    piagu.print_ndarr(arr, 'arr', last=10)


def test_photons_3d() :
    from time import time

    mean, sigma = 0,10
    shape = (32,185,388)
    #shape = (4,512,512)
    mask  = np.ones(shape, dtype=np.uint8)
    data  = np.array(mean + sigma * np.random.standard_normal(size=shape), dtype=np.float32)

    piagu.print_ndarr(data, 'data')
    piagu.print_ndarr(mask, 'mask')

    t0_sec = time()
    arr3d = photons(data, mask)
    print('\nTime consumed by photons(data, mask) (sec) = %10.6f' % (time()-t0_sec))


    arr2d = arr3d[1,:]
    arr2d.shape = (shape[-2], shape[-1])

    piagu.print_ndarr(arr3d, 'arr3d')
    piagu.print_ndarr(arr2d, 'arr2d')

    from pyimgalgos.GlobalGraphics import plotImageLarge, show
    #plotImageLarge(data, title='data')
    plotImageLarge(arr2d, title='arr2d')
    show()


def test_pyalgos() :

    from time import time

    t0_sec = time()

    windows = ((0, 0, 185, 0, 388), \
               (1, 0, 185, 0, 388), \
               (3, 0, 185, 0, 388))
    print_arr(windows, "windows")

    alg = PyAlgos(windows, pbits=0)
    alg.print_attributes()

    print('\nC++ consumed time to get raw data (sec) = %10.6f' % (time()-t0_sec))



if __name__ == "__main__" :

    if len(sys.argv)==1    :
        print('For test(s) use command: python', sys.argv[0], '<test-number=1-...>')
        test_pyalgos()

    elif(sys.argv[1]=='1') : test_pyalgos()
    elif(sys.argv[1]=='2') : test_photons_2d()
    elif(sys.argv[1]=='3') : test_photons_2d_chuck()
    elif(sys.argv[1]=='4') : test_photons_3d()

    else : print('Non-expected arguments: sys.argv=', sys.argv, ' use 0,1,2,...')

    sys.exit ('Self test is done.')

# EOF
