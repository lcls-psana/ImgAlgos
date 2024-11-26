#!/usr/bin/env python

import sys
import numpy as np
from time import time

from ImgAlgos.PyAlgos import PyAlgos, print_arr, print_arr_attr

hdr = 'Evnum  Reg  Seg  Row  Col  Npix      Amax      Atot   rcent   ccent '+\
      'rsigma  csigma rmin rmax cmin cmax    bkgd     rms     son' # +\
      #'  imrow   imcol     x[um]     y[um]     r[um]  phi[deg]'

fmt = '%5d  %3s  %3d %4d %4d  %4d  %8.1f  %8.1f  %6.1f  %6.1f %6.2f  %6.2f'+\
      ' %4d %4d %4d %4d  %6.2f  %6.2f  %6.2f' # +\
      #' %6d  %6d  %8.0f  %8.0f  %8.0f  %8.2f'


def test_02() :
    #winds = ((1, 0, 185, 0, 388), \
    #         (1, 0, 185, 0, 388))
    winds = None

    #print_arr(winds, 'windows')

    mu, sigma, shape = 0, 20, (2, 185, 388)
    data = np.array(mu + sigma*np.random.standard_normal(shape), dtype=np.float64)
    mask = np.ones(shape) # or None

    #fname = 'cspad2x2-random.npy'
    #np.save(fname, data)
    #print 'Random image saved in file %s' % fname

    alg = PyAlgos(windows=winds, mask=mask, pbits=0)
    #alg = PyAlgos()

    alg.print_attributes()
    #alg.set_windows(windows)

    print_arr_attr(data, 'data')

    thr = 20
    t0_sec = time()
    n1 = alg.number_of_pix_above_thr(data, thr)
    print('%s\n  alg.number_of_pix_above_thr = %d, fr = %8.6f' % (80*'_', n1, float(n1)/data.size))
    print('  Time consumed by the test = %10.6f(sec)' % (time()-t0_sec))

    t0_sec = time()
    a1 = alg.intensity_of_pix_above_thr(data, thr)
    print('%s\n  alg.intensity_of_pix_above_thr = %12.3f' % (80*'_', a1))
    print('  Time consumed by the test = %10.6f(sec)' % (time()-t0_sec))


def test_01() :
    print('%s\n%s\n' % (80*'_','test_01'))

    import pyimgalgos.GlobalGraphics as gg

    fig, axim, axcb = gg.fig_axes() # if not do_plot else (None, None, None)

    #shape = (2, 185, 388)
    shape = (100, 100)
    mask = np.ones(shape)

    #winds = [(s, 0, 185, 0, 388) for s in (0,1,7,8,9,15,16,17,23,24,25,31)]
    winds = None

    alg = PyAlgos(windows=winds, mask=mask, pbits=0)
    #alg.set_peak_selection_pars(npix_min=5, npix_max=500, amax_thr=0, atot_thr=1000, son_min=6)
    #alg = PyAlgos()
    alg.print_attributes()
    #alg.set_windows(windows)

    rank=5
    mu, sigma = 0, 20

    map_stat = np.zeros((8,), dtype=np.float64)
    counter = 0
    for i in range(10) :

        print('\nEvent # %d' % i)

        nda = np.array(mu + sigma*np.random.standard_normal(shape), dtype=np.float64)
        #print_arr_attr(nda, 'nda')

        t0_sec = time()
        peaks = alg.peak_finder_v3(nda, rank=rank, r0=5.0, dr=0.05)
        print('  Time consumed by the peak_finder = %10.6f(sec)' % (time()-t0_sec))

        maps = alg.maps_of_local_maximums()
        #print_arr_attr(maps, 'maps')
        maps.shape = shape

        map_stat += np.bincount(maps.flatten(), weights=None, minlength=None)
        counter += 1

        print(hdr)
        reg = 'IMG'

        for pk in peaks :
            # get peak parameters
            seg,row,col,npix,amax,atot,rcent,ccent,rsigma,csigma,\
            rmin,rmax,cmin,cmax,bkgd,rms,son = pk[0:17]

            rec = fmt % (i, reg, seg, row, col, npix, amax, atot, rcent, ccent, rsigma, csigma,\
                  rmin, rmax, cmin, cmax, bkgd, rms, son) #,\
                  #imrow, imcol, xum, yum, rum, phi)
            print(rec)

        img, amin, amax = maps, -3, 7
        #img = nda
        #ave, rms = img.mean(), img.std()
        #amin, amax = ave-2*rms, ave+2*rms

        gg.plot_img(img, mode='do not hold', amin=amin, amax=amax)
        fig.canvas.set_window_title('Event: %d' % i)
        fig.canvas.draw() # re-draw figure content


    map_stat = map_stat / counter
    map_frac = map_stat / np.sum(map_stat)
    print('map_stat', map_stat)
    print('rank=%d,  fractions=%s' % (rank, map_frac))

    gg.save_fig(fig, fname='map-locmax-100x100-rank-%d.png'%rank, pbits=1)
    #gg.save_fig(fig, fname='map-locmax-100x100-random.png', pbits=1)

    gg.show() # hold image untill it is closed


def usage() : return 'Use command: python %s <test-number>, where <test-number> = 1,2,...,9,...' % sys.argv[0]


def main() :
    if len(sys.argv) != 2  : test_01(); print('\n%s\n%s\n' %  (80*'_', usage()))
    elif sys.argv[1] =='1' : test_01()
    elif sys.argv[1] =='2' : test_02()
    else                   : print('\n%s\nTest id parameter is not recognized.\n%s' % (80*'_', usage()))


if __name__ == "__main__" :
    main()
    sys.exit('\nEnd of test')

# EOF
