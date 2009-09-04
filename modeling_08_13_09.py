from pylab import *
from scipy import *
import sys
import txt_data_processing, rwkos
import SLFR_TMM
reload(SLFR_TMM)
import pylab_util

import rwkbode

data_dir = rwkos.FindFullPath('siue/Research/modeling/SLFR/data/July_07_2009')
if data_dir not in sys.path:
    sys.path.append(data_dir)

#data_mod_name = 'swept_sine_amp_75_July_07_2009_avebodes'
data_mod_name = 'swept_sine_amp_75_July_07_2009_log_downsampled'
bode_data_set = txt_data_processing.load_avebode_data_set(data_mod_name)

bode_data_set.Bode_Plot2(func=rwkbode.GenBodePlot)

k_spring = 1.0
c_spring = 0.1
K_act = 0.05
tau = 10.0*2*pi
xf = [k_spring, c_spring, K_act, tau, 1.5]
#f_TMM = bode_data_set.f[1:]#drop 0 from f
f_TMM = arange(0.1,50,0.1)
ol_model = SLFR_TMM.SLFR_TMM_OL_model(xf, include_spring=1)
bodes = ol_model.BodeResponse(f_TMM)
#act_model = SLFR_TMM.AVS_only_model(K_act, tau)
#act_bode = act_model.BodeResponse(f_TMM)[0]

theta_v_bode = bodes[1]
theta_v_bode.seedphase = 0.0
theta_v_bode.seedfreq = 1.0
theta_v_bode.PhaseMassage(f_TMM)

accel_v_bode = bodes[0]
accel_v_bode.seedfreq = 1.0
accel_v_bode.seedphase = 80.0
accel_v_bode.PhaseMassage(f_TMM)

freqlim = [0.5, 50.0]
#rwkbode.GenBodePlot(1, f_TMM, act_bode, clear=False)
rwkbode.GenBodePlot(1, f_TMM, theta_v_bode, clear=False)

pylab_util.SetPhaseLim(1, [-300,0])
pylab_util.SetFreqLim(1, freqlim)
pylab_util.SetMagLim(1, [-50, 20])

resave = 0
if resave:
    pylab_util.mysave('theta_vs_v_bode.png', 1)
    

rwkbode.GenBodePlot(2, f_TMM, accel_v_bode, clear=False)

pylab_util.SetPhaseLim(2, [-450,100])
pylab_util.SetFreqLim(2, freqlim)
pylab_util.SetMagLim(2, [-20, 30])

if resave:
    pylab_util.mysave('accel_vs_v_bode.png', 2)


a_theta_bode = accel_v_bode/theta_v_bode
a_theta_bode.seedfreq = 10.0
a_theta_bode.seedphase = 0.0
a_theta_bode.PhaseMassage(f_TMM)
a_theta_bode.seedfreq = 30.0
a_theta_bode.seedphase = -180.0
a_theta_bode.PhaseMassage(f_TMM)
#unwrap(a_theta_bode.phase, axis=0)

rwkbode.GenBodePlot(3, f_TMM, a_theta_bode, clear=False)
pylab_util.SetPhaseLim(3, [-250,250])
pylab_util.SetFreqLim(3, freqlim)
pylab_util.SetMagLim(3, [-20, 50])

show()
