from pylab import *
from scipy import *
from scipy import optimize
import sys
import txt_data_processing, rwkos
import SLFR_TMM
reload(SLFR_TMM)
import pylab_util

from SLFR_TMM import E, I, mu, L, rho

import rwkbode

from IPython.Debugger import Pdb

import time
t1 = time.time()

import sys

if len(sys.argv) > 1:
    run_fit = int(sys.argv[1])
else:
    run_fit = 1
    
data_dir = rwkos.FindFullPath('siue/Research/modeling/SLFR/data/July_07_2009')
if data_dir not in sys.path:
    sys.path.append(data_dir)

#data_mod_name = 'swept_sine_amp_75_July_07_2009_avebodes'
data_mod_name = 'swept_sine_amp_75_July_07_2009_log_downsampled'
bode_data_set = txt_data_processing.load_avebode_data_set(data_mod_name)

#bode_data_set.Bode_Plot2(func=rwkbode.GenBodePlot, linetype='o')
exp_bodes = bode_data_set.avebodes[0:2]
th_v_exp = exp_bodes[0]
a_v_exp = exp_bodes[1]
f = bode_data_set.f
a_theta_exp = bode_data_set.avebodes[2]

rwkbode.GenBodePlot(1, f, th_v_exp, clear=True, linetype='o')
rwkbode.GenBodePlot(2, f, a_v_exp, clear=True, linetype='o')
rwkbode.GenBodePlot(3, f, a_theta_exp, clear=True, linetype='o')

t2 = time.time()

def a_theta_massage(a_th_bode, fvect):
    phase = a_th_bode.phase
    b1 = a_th_bode.phase < -100
    b2 = fvect < 4
    b = b1 & b2
    phase[b] += 360.0
    
    b3 = a_th_bode.phase > 100
    b4 = fvect > 9
    b5 = b3 & b4
    phase[b5] -= 360.0
    a_th_bode.phase = phase

def a_theta_model(a_v_bode, th_v_bode, fvect):
    a_th_bode = a_v_bode/th_v_bode
    a_th_bode.seedfreq = 10.0
    a_th_bode.seedphase = 0.0
    a_th_bode.PhaseMassage(fvect)
    a_theta_massage(a_th_bode, fvect)
    return a_th_bode

    
def mymodel(ucv, fvect=f, actvect=[]):
    #ol_model = SLFR_TMM.SLFR_TMM_OL_model(ucv, include_spring=1, \
    #                                      actvect=actvect)
    ol_model = SLFR_TMM.SLFR_TMM_OL_model_v2(ucv, include_spring=1, \
                                             actvect=actvect)
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

    a_v_theta_bode = a_theta_model(accel_v_bode, theta_v_bode, fvect)
    return theta_v_bode, accel_v_bode, a_v_theta_bode


def mycost(ucv, actvect=[], phaseweight=0.1):
    model_bodes = mymodel(ucv, actvect=actvect)
    totale = 0.0
    exp_list = [th_v_exp, a_theta_exp]#a_v_exp]
    model_list = [model_bodes[0], model_bodes[2]]
    for exp_bode, model_bode in zip(exp_list, model_list):
        magE = squeeze(exp_bode.dBmag())-squeeze(model_bode.dBmag())
        phaseE = squeeze(exp_bode.phase)-squeeze(model_bode.phase)
        totale += sum(magE**2)+phaseweight*sum(phaseE**2)
    if ucv[0] < 0.0:
        totale += exp(-50.0*ucv[0])-1
    return totale


def RunFit(phaseweight=0.1,ig=None, maxiter=None, maxfun=None):
    t1=mt()
    if ig is None:
        ig=initialguesses
    fitres=fmin(mycost,ig,(phaseweight,), maxiter=maxiter, maxfun=maxfun)
    t2=mt()
    print('RunFit time='+str(t2-t1))
    return fitres


t3 = time.time()

k_spring = 1.0
c_spring = 0.1
K_act = 0.05
tau = 10.0*2*pi
#xf = [k_spring, c_spring, K_act, tau, 1.5]
#xf = [k_spring, c_spring, 1.5, 1.0/32*25.4/1000.0]
#xf = [k_spring, c_spring, K_act, tau, 1.02, 1.5, 45.0]
xf = [k_spring, c_spring, 1.02, 1.5, 45.0]
actvect = [K_act, tau]
initial_cost = mycost(xf, actvect=actvect)
print('initial_cost = %s' % initial_cost)

l1 = log10(0.5)
l2 = log10(50)
flist = arange(0.5, 1, 0.5).tolist()+f.tolist()+arange(26,30, 0.5).tolist()
f2 = array(flist)

theta_v_bode, accel_v_bode, a_theta_bode = mymodel(xf, f2, actvect=actvect)

t4 = time.time()

freqlim = [0.5, 50.0]
#rwkbode.GenBodePlot(1, f, act_bode, clear=False)
rwkbode.GenBodePlot(1, f2, theta_v_bode, clear=False, linetype='g-')
rwkbode.GenBodePlot(2, f2, accel_v_bode, clear=False, linetype='g-')
rwkbode.GenBodePlot(3, f2, a_theta_bode, clear=False, linetype='g-')

if run_fit:
    xf_fit = optimize.fmin(mycost, xf, args=(actvect, 0.1))
    th_v_fit, a_v_fit, a_th_fit = mymodel(xf_fit, f2, actvect=actvect)

    rwkbode.GenBodePlot(1, f2, th_v_fit, clear=False, linetype='r-')
    rwkbode.GenBodePlot(2, f2, a_v_fit, clear=False, linetype='r-')
    rwkbode.GenBodePlot(3, f2, a_th_fit, clear=False, linetype='r-')
    
t5 = time.time()


pylab_util.SetPhaseLim(1, [-300,0])
pylab_util.SetFreqLim(1, freqlim)
pylab_util.SetMagLim(1, [-50, 20])

resave = 0
if resave:
    pylab_util.mysave('theta_vs_v_bode.png', 1)
    

pylab_util.SetPhaseLim(2, [-450,100])
pylab_util.SetFreqLim(2, freqlim)
pylab_util.SetMagLim(2, [-20, 30])

if resave:
    pylab_util.mysave('accel_vs_v_bode.png', 2)
    

pylab_util.SetPhaseLim(3, [-250,250])
pylab_util.SetFreqLim(3, freqlim)
pylab_util.SetMagLim(3, [-20, 50])

t6 = time.time()

tlist = [t1, t2, t3, t4, t5]

pt = t1

t_diffs = []

for cur_t in tlist[1:]:
    dt = cur_t-pt
    t_diffs.append(dt)
    pt = cur_t
    
show()
