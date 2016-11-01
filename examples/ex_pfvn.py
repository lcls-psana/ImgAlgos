#------------------------------
### #!/usr/bin/env python

from time import time
import numpy as np
from pyimgalgos.NDArrGenerators import random_standard, add_random_peaks, reshape_to_2d
from ImgAlgos.PyAlgos import PyAlgos, print_arr, print_arr_attr

import pyimgalgos.GlobalGraphics as gg

#------------------------------

def plot_image(img, img_range=None, amp_range=None, figsize=(12,10)) : 
    #import pyimgalgos.GlobalGraphics as gg
    axim = gg.plotImageLarge(img, img_range, amp_range, figsize)
    gg.show() 

#------------------------------

def image_with_random_peaks(shape=(500, 500)) : 
    img = random_standard(shape, mu=0, sigma=10)
    peaks = add_random_peaks(img, npeaks=10, amean=100, arms=50, wmean=1.5, wrms=0.3)
    return img, peaks

#------------------------------

hdr = 'Evnum  Reg  Seg  Row  Col  Npix      Amax      Atot   rcent   ccent '+\
      'rsigma  csigma rmin rmax cmin cmax    bkgd     rms     son' # +\
      #'  imrow   imcol     x[um]     y[um]     r[um]  phi[deg]'

fmt = '%5d  %3s  %3d %4d %4d  %4d  %8.1f  %8.1f  %6.1f  %6.1f %6.2f  %6.2f'+\
      ' %4d %4d %4d %4d  %6.2f  %6.2f  %6.2f' # +\
      #' %6d  %6d  %8.0f  %8.0f  %8.0f  %8.2f'

#------------------------------

def test_pf() : 

    ##-----------------------------

    SKIP        = 0
    EVTMAX      = 10 + SKIP

    DO_PLOT = True

    shape=(1000, 1000)

    # Pixel image indexes
    #arr3d = np.array((1,shape[0],shape[1]))

    INDS = np.indices((1,shape[0],shape[1]), dtype=np.int64)
    imRow, imCol = INDS[1,:], INDS[2,:]  
    #iX  = np.array(det.indexes_x(evt), dtype=np.int64) #- xoffset
    #iY  = np.array(det.indexes_y(evt), dtype=np.int64) #- yoffset

    ##-----------------------------
    fig3, axim3, axcb3 = gg.fig_axes(figsize=(11,10)) if DO_PLOT else (None, None, None)
    fig2, axim2, axcb2 = gg.fig_axes(figsize=(11,10)) if DO_PLOT else (None, None, None)
    fig,  axim,  axcb  = gg.fig_axes(figsize=(11,10)) if DO_PLOT else (None, None, None)
    ##-----------------------------

    alg = PyAlgos(windows=None, mask=None, pbits=0)
    alg.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=0, son_min=6)
    #alg.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=500, son_min=6) # for v2r1
    alg.print_attributes()

    for ev in range(EVTMAX) :

        if ev<SKIP : continue
        #if ev>=EVTMAX : break

        print 50*'_', '\nEvent %04d' % ev

        img, peaks_sim = image_with_random_peaks(shape)
        peaks_gen = [(0, r, c, a, a*s, 9*s*s) for r,c,a,s in peaks_sim]

        t0_sec = time()
        #peaks = alg.peak_finder_v3r1(img, rank=5, r0=7, dr=2)

        #peaks = alg.peak_finder_v2r1(img, thr=30, r0=7, dr=2)
        #peaks = alg.peak_finder_v3r1(img, rank=5, r0=7, dr=2, nsigm=0) # 1.64 (5%)
        #peaks = alg.peak_finder_v4r1(img, thr_low=20, thr_high=50, rank=5, r0=7, dr=2)
        peaks  = alg.peak_finder_v4r2(img, thr_low=20, thr_high=50, rank=6, r0=7, dr=2)

        print '  Time consumed by the peak_finder = %10.6f(sec)' % (time()-t0_sec)

        map2 = reshape_to_2d(alg.maps_of_pixel_status())
        map3 = reshape_to_2d(alg.maps_of_connected_pixels())
        #map2 = np.zeros((10,10))
        #map3 = np.zeros((10,10))
        #maps = alg.maps_of_local_minimums()
        #print_arr(map2, 'map_of_pixel_status')
        #print_arr(map3, 'map_of_connected_pixels')
        #maps.shape = shape 


        print 'Simulated peaks:'
        for i, (r0, c0, a0, sigma) in enumerate(peaks_sim) :
            print '  %04d  row=%6.1f  col=%6.1f  amp=%6.1f  sigma=%6.3f' % (i, r0, c0, a0, sigma)
        #plot_image(img)

        print 'Found peaks:'
        print hdr
        reg = 'IMG'
        for pk in peaks :
            seg,row,col,npix,amax,atot,rcent,ccent,rsigma,csigma,\
            rmin,rmax,cmin,cmax,bkgd,rms,son = pk[0:17]

            rec = fmt % (ev, reg, seg, row, col, npix, amax, atot, rcent, ccent, rsigma, csigma,\
                  rmin, rmax, cmin, cmax, bkgd, rms, son) #,\
                  #imrow, imcol, xum, yum, rum, phi)
            print rec


        if DO_PLOT :

            #nda = maps_of_conpix_arc        
            #nda = maps_of_conpix_equ        
            #img = det.image(evt, nda)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            #img = det.image(evt, mask_img*nda)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            #img = det.image(evt, maps_of_conpix_equ)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]

            imsh2 = axim2.imshow(map2, interpolation='nearest', aspect='auto', origin='upper') 
            imsh2.set_clim(0, 30)
            colb = fig2.colorbar(imsh2, cax=axcb2) # , orientation='horizontal')
            fig2.canvas.set_window_title('Pixel status')    

            imsh3 = axim3.imshow(map3, interpolation='nearest', aspect='auto', origin='upper') 
            colb = fig3.colorbar(imsh3, cax=axcb3) # , orientation='horizontal')
            fig3.canvas.set_window_title('Connected pixel groups')    
            gg.move_fig(fig3, x0=200, y0=30)

            ave, rms = img.mean(), img.std()
            amin, amax = ave-1*rms, ave+8*rms
            gg.plot_img(img, mode='do not hold', amin=amin, amax=amax)
            gg.plot_peaks_on_img(peaks_gen, axim, imRow, imCol, color='g', lw=5) #, pbits=3)
            gg.plot_peaks_on_img(peaks, axim, imRow, imCol, color='w') #, pbits=3)
            #gg.plot_peaks_on_img(peaks, axim, imRow, imCol, color='w') #, pbits=3)
            #gg.plot_peaks_on_img(peaks_equ, axim, imRow, imCol, color='w') #, pbits=3)

            #gg.plotHistogram(nda, amp_range=(-100,100), bins=200, title='Event %d' % i)

            fig.canvas.set_window_title('Event: %d' % ev)    
            fig.canvas.draw() # re-draw figure content
            gg.move_fig(fig, x0=400, y0=30)


    gg.show()
 
#------------------------------
#------------------------------
#------------------------------
#------------------------------

def ex_image_with_random_peaks() :     
    img, peaks = image_with_random_peaks()
    print 'peaks:'
    for i, (r0, c0, a0, sigma) in enumerate(peaks) :
        print '  %04d  row=%6.1f  col=%6.1f  amp=%6.1f  sigma=%6.3f' % (i, r0, c0, a0, sigma)
    plot_image(img)

#------------------------------

if __name__ == "__main__" :
    import sys; global sys
    tname = sys.argv[1] if len(sys.argv) > 1 else '0'
    print 50*'_', '\nTest %s:' % tname
    if   tname == '0' : ex_image_with_random_peaks()
    elif tname == '1' : test_pf()
    else : print 'Not-recognized test name: %s' % tname
    sys.exit('End of test %s' % tname)
 
#------------------------------
