from pylab import *
from scipy import *

import exp_data

f2 = exp_data.f2

#exp_data.plot_exp()

no_damping_model = SLFR_TMM.AVS_and_Beam(K, tau, a_gain, c=0.0)

bodes = ol_model.BodeResponse(fvect)
#act_model = SLFR_TMM.AVS_only_model(K_act, tau)
#act_bode = act_model.BodeResponse(f)[0]

theta_v_bode = bodes[1]
theta_v_bode.seedphase = 0.0
theta_v_bode.seedfreq = 2.0
theta_v_bode.PhaseMassage(fvect)

accel_v_bode = bodes[0]
accel_v_bode.seedfreq = 2.0
accel_v_bode.seedphase = 80.0
accel_v_bode.PhaseMassage(fvect)
