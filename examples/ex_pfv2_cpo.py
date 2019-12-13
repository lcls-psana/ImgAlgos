#!/usr/bin/env python
"""Example for peak_finder_v2
"""
from __future__ import print_function
#------------------------------

import sys
import numpy as np
import psana

from ImgAlgos.PyAlgos import PyAlgos
from Detector.GlobalUtils import print_ndarr

#------------------------------

ds  = psana.DataSource('exp=cxif5315:run=169')
det = psana.Detector('CxiDs2.0:Cspad.0', ds.env())

#src = psana.Source('DetInfo(CxiDs2.0:Cspad.0)')
#from Detector.AreaDetector import AreaDetector
#det = AreaDetector(src, ds.env(), pbits=0)

windows = [(s, 0, 185, 0, 388) for s in (0,1,7,8,9,15,16,17,23,24,25,31)] # or None
#mask_arc = np.loadtxt('/reg/neh/home/cpo/ipsana/cxif5315/masks/roi_mask_nda_arc.txt').reshape((32,185,388))
mask_arc = np.ones((32,185,388))
alg_arc = PyAlgos(windows=windows, mask=mask_arc, pbits=0) # pbits=0177777
alg_arc.set_peak_selection_pars(npix_min=5, npix_max=500, amax_thr=0, atot_thr=1000, son_min=6)
#alg_arc.print_input_pars()

for i,evt in enumerate(ds.events()):

    if i>20 : break

    nda = det.calib(evt)
    print('%s\nEvent # %d\n' % (80*'_',i))
    print_ndarr(nda,  'data ndarray')
    
    # r,dr setup a "quantized ring" which will be used for background evaluation
    peaks = alg_arc.peak_finder_v2(nda, thr=20, r0=5.0, dr=0.05)
    print_ndarr(peaks,  'returned peaks seen in ex_pf2_cpo.py')
    #print peaks

sys.exit('Test is completed')

#------------------------------
