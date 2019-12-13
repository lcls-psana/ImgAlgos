#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
import sys
import math
import numpy as np
from time import time
import psana

from Detector.AreaDetector import AreaDetector
from ImgAlgos.PyAlgos import PyAlgos, print_arr, print_arr_attr
from pyimgalgos.PeakStore import PeakStore
from pyimgalgos.GlobalUtils import subtract_bkgd

##-----------------------------
# Initialization of graphics
from pyimgalgos.GlobalGraphics import store as sp
import pyimgalgos.GlobalGraphics as gg
##-----------------------------

ntest = int(sys.argv[1]) if len(sys.argv)>1 else 1
print('Test # %d' % ntest)

##-----------------------------
SKIP        = 0
EVTMAX      = 10 + SKIP
#SKIP        = 70024
#EVTMAX      = 1000000 + SKIP
EVTPLOT     = 1 
DO_PLOT     = True
##-----------------------------

def do_print(i) :
    return True
    #return False
    #if i==1 : return True
    #return not i%10

##-----------------------------
    
dsname = 'exp=cxif5315:run=169'
src    = psana.Source('DetInfo(CxiDs2.0:Cspad.0)')
print('%s\nExample for\n  dataset: %s\n  source : %s' % (85*'_',dsname, src))

# Non-standard calib directory
#psana.setOption('psana.calib-dir', './empty/calib')
psana.setOption('psana.calib-dir', '/reg/d/psdm/CXI/cxif5315/calib')

ds  = psana.DataSource(dsname)
evt = next(ds.events())
env = ds.env()

runnum = evt.run()

#run = ds.runs().next()
#runnum = run.run()

#for key in evt.keys() : print key

##-----------------------------

det = AreaDetector(src, env, pbits=0)
print(85*'_', '\nInstrument: %s  run number: %d' % (det.instrument(), runnum))

nda_peds  = det.pedestals(evt)
nda_bkgd  = det.bkgd(evt)
nda_smask = det.mask(evt, calib=False, status=True, edges=True, central=True, unbond=True, unbondnbrs=True)

#print_arr_attr(nda_peds, 'nda_peds')
#print_arr_attr(nda_bkgd, 'nda_bkgd')
#print_arr_attr(nda_smask, 'nda_smask')
##-----------------------------

shape_cspad = (32,185,388)

mask_arc = np.loadtxt('../rel-mengning/work/roi_mask_nda_arc.txt')
mask_equ = np.loadtxt('../rel-mengning/work/roi_mask_nda_equ.txt')
mask_img = np.loadtxt('../rel-mengning/work/roi_mask_nda_equ_arc.txt')
mask_arc.shape = mask_equ.shape = mask_img.shape = shape_cspad
print_arr_attr(mask_arc, 'mask_arc')


seg1 = np.ones((185,388))
regs_check = np.zeros(shape_cspad)
for s in (4,12,20,28) : regs_check[s,10:100,270:370] = 20*seg1[10:100,270:370]

winds_bkgd = [ (s, 10, 100, 270, 370) for s in (4,12,20,28)] # use part of segments 4 and 20 to subtr bkgd

winds_arc  = [ (s, 0, 185, 0, 388) for s in (0,1,7,8,9,15,16,17,23,24,25,31)]
winds_equ  = [ (s, 0, 185, 0, 388) for s in (0,1,3,8,9,11,16,17,19,24,25,27)]
#winds_all = [ (s, 0, 185, 0, 388) for s in (0,1,3,7,8,9,11,15,16,17,19,23,24,25,27,31)]

print_arr(winds_arc, 'winds_arc')
print_arr_attr(winds_arc, 'winds_arc')

alg_arc = PyAlgos(windows=winds_arc, mask=mask_arc, pbits=2)
alg_arc.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=0, son_min=10)
#alg_arc.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=500, son_min=6) # for v2r1

alg_equ = PyAlgos(windows=winds_equ, mask=mask_equ, pbits=0)
alg_equ.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=0, son_min=10)
#alg_equ.set_peak_selection_pars(npix_min=0, npix_max=1e6, amax_thr=0, atot_thr=500, son_min=6) # for v2r1



#alg_equ.print_attributes()
#alg_equ.print_input_pars()

##-----------------------------

xoffset, yoffset = 300, 300
xsize,   ysize   = 1150, 1150

# Pixel image indexes
iX  = np.array(det.indexes_x(evt), dtype=np.int64) #- xoffset
iY  = np.array(det.indexes_y(evt), dtype=np.int64) #- yoffset

# Protect indexes (should be POSITIVE after offset subtraction)
imRow = np.select([iX<xoffset], [0], default=iX-xoffset)
imCol = np.select([iY<yoffset], [0], default=iY-yoffset)

# Pixel coordinates [um] (transformed as needed)
Xum =  det.coords_y(evt)
Yum = -det.coords_x(evt)

# Derived pixel raduius in [um] and angle phi[degree]
Rum = np.sqrt(Xum*Xum + Yum*Yum)
Phi = np.arctan2(Yum,Xum) * 180 / np.pi

imRow.shape  = imCol.shape  = \
Xum.shape    = Yum.shape    = \
Rum.shape    = Phi.shape    = shape_cspad

addhdr = 'Evnum  Reg  Seg  Row  Col  Npix      Amax      Atot   rcent   ccent '+\
         'rsigma  csigma rmin rmax cmin cmax    bkgd     rms     son  imrow   imcol     x[um]     y[um]     r[um]  phi[deg]'

fmt = '%5d  %3s  %3d %4d %4d  %4d  %8.1f  %8.1f  %6.1f  %6.1f %6.2f  %6.2f'+\
      ' %4d %4d %4d %4d  %6.2f  %6.2f  %6.2f'+\
      ' %6d  %6d  %8.0f  %8.0f  %8.0f  %8.2f'

pstore = PeakStore(env, runnum, prefix='peaks', add_header=addhdr, pbits=0)
pstore.print_attrs()

##-----------------------------
fig, axim, axcb = gg.fig_axes() if DO_PLOT else (None, None, None)
##-----------------------------

def geo_pars(s,r,c) :
    inds = (s,r,c)
    return imRow[inds], imCol[inds], Xum[inds], Yum[inds], Rum[inds], Phi[inds]
    #return imRow[s,r,c], imCol[s,r,c], Xum[s,r,c], Yum[s,r,c], Rum[s,r,c], Phi[s,r,c]

##-----------------------------

t0_sec_evloop = time()
nda = None
peaks = None

# loop over events in data set
for i, evt in enumerate(ds.events()) :

    if do_print(i) and i%100==0 : print('Event %d' % (i))

    if i<SKIP    : continue
    if i>=EVTMAX : break

    # get calibrated data ndarray and proccess it if it is available
    t1_sec = time()
    #nda = det.calib(evt)

    # Apply custom calibration: raw, -peds, -bkgd, *smask, -cmod
    nda_raw = det.raw(evt)

    if nda_raw is not None :

        nda =  np.array(nda_raw, dtype=np.float32, copy=True)
        nda -= nda_peds
        nda =  subtract_bkgd(nda, nda_bkgd, mask=nda_smask, winds=winds_bkgd, pbits=0)
        nda *= nda_smask
        det.common_mode_apply(evt, nda)

        #print '  ----> calibration dt = %f sec' % (time()-t1_sec)


        #print_arr_attr(nda, 'calibrated data')
        t0_sec = time()

        # run peakfinders and get list of peak records for each region
        #peaks_arc = alg_arc.peak_finder_v2r1(nda, thr=30, r0=7, dr=2)
        #peaks_arc = alg_arc.peak_finder_v3r1(nda, rank=5, r0=7, dr=2, nsigm=0) # 1.64 (5%)
        #peaks_arc = alg_arc.peak_finder_v4r1(nda, thr_low=10, thr_high=150, rank=5, r0=7, dr=2)
        peaks_arc = alg_arc.peak_finder_v4r2(nda, thr_low=10, thr_high=150, rank=5, r0=7, dr=2)

        #peaks_equ = alg_equ.peak_finder_v2r1(nda, thr=30, r0=7, dr=2)
        #peaks_equ = alg_equ.peak_finder_v3r1(nda, rank=5, r0=7, dr=2, nsigm=0) # 1.64 (5%)
        #peaks_equ = alg_equ.peak_finder_v4r1(nda, thr_low=10, thr_high=150, rank=5, r0=7, dr=2)
        peaks_equ = alg_equ.peak_finder_v4r2(nda, thr_low=10, thr_high=150, rank=5, r0=7, dr=2)
 
        # available after v2r1, v4r2
        #maps_of_conpix_arc = alg_arc.maps_of_connected_pixels()
        #maps_of_conpix_equ = alg_equ.maps_of_connected_pixels()

        # available after v3r1 ONLY!
        #maps_of_locmax_arc = alg_arc.maps_of_local_maximums()
        #maps_of_locmax_equ = alg_equ.maps_of_local_maximums()

        ###===================
        if do_print(i) : print('%s\n%s\n%s\n%s' % (85*'_', pstore.header[0:66], pstore.rec_evtid(evt), addhdr))
        ###===================

        peak_reg_lists = list(zip(('ARC','EQU'), (peaks_arc, peaks_equ))) 

        # loop over ARC and EQU regions
        for reg, peak_list in peak_reg_lists :

            # loop over peaks found in the region
            for peak in peak_list :

                # get peak parameters
                seg,row,col,npix,amax,atot,rcent,ccent,rsigma,csigma,\
                rmin,rmax,cmin,cmax,bkgd,rms,son = peak[0:17]

                # get pixel coordinates
                imrow, imcol, xum, yum, rum, phi = geo_pars(seg, row, col)
                
                # make peak-record and save it in the file
                rec = fmt % (i, reg, seg, row, col, npix, amax, atot, rcent, ccent, rsigma, csigma,\
                      rmin, rmax, cmin, cmax, bkgd, rms, son,\
                      imrow, imcol, xum, yum, rum, phi)
            
                pstore.save_peak(evt, rec)

                ###===================
                if do_print(i) : print('%s' % rec)
                ###===================

        ###===================
        if do_print(i) : print('Event %d --- dt/evt = %f sec' % (i, time()-t0_sec))
        ###===================

        if DO_PLOT and i%EVTPLOT==0 :

            #nda = maps_of_conpix_arc        
            #nda = maps_of_conpix_equ        
            #nda = nda_bkgd
            #nda = nda_bkgd + regs_check      
            #img = det.image(evt, nda)
            #img = det.image(evt, nda)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            img = det.image(evt, mask_img*nda)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            #img = det.image(evt, maps_of_conpix_equ)[xoffset:xoffset+xsize,yoffset:yoffset+ysize]
            ave, rms = img.mean(), img.std()
            amin, amax = ave-1*rms, ave+8*rms
            gg.plot_img(img, mode='do not hold', amin=amin, amax=amax)
            gg.plot_peaks_on_img(peaks_arc, axim, imRow, imCol, color='w') #, pbits=3)
            gg.plot_peaks_on_img(peaks_equ, axim, imRow, imCol, color='w') #, pbits=3)

            #gg.plotHistogram(nda, amp_range=(-100,100), bins=200, title='Event %d' % i)

            fig.canvas.set_window_title('Event: %d' % i)    
            fig.canvas.draw() # re-draw figure content


print(' ----> Event loop time = %f sec' % (time()-t0_sec_evloop))

pstore.close_file()

##-----------------------------

gg.show() # hold image untill it is closed
 
##-----------------------------

sys.exit('Test is completed')

##-----------------------------
