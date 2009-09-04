from pylab import *
from scipy import *

import os

import rwkos

import txt_data_processing
reload(txt_data_processing)

from txt_data_processing import load_avebode_data_set, \
     load_time_domain_data_set, merge_trunc_ave_data_sets

import pylab_util

import rwkbode

data_dir = rwkos.FindFullPath('siue/Research/modeling/SLFR/data/July_07_2009')
if data_dir not in sys.path:
    sys.path.append(data_dir)


import SLFR_TMM
reload(SLFR_TMM)

ds_name = 'swept_sine_amp_75_July_07_2009_log_downsampled'

data_set = load_avebode_data_set(ds_name)
figs = data_set.Bode_Plot2(linetype='o')


k_spring = 1.0
c_spring = 0.1
K_act = 0.05
tau = 10.0*2*pi
xf = [k_spring, c_spring, K_act, tau, 1.5]

ol_model = SLFR_TMM.SLFR_TMM_OL_model(xf)
ol_model.IntegratedCurveFit(expmodname=data_set, prefix='log_compressed')
                   
show()
