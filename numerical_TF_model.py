from pylab import *
from scipy import *

import controls

w1=2.5*2*pi
z1=0.03
g1=0.003

a_th_mode1 = controls.TransferFunction([g1*w1**2,0,0],[1,2*z1*w1,w1**2])

w2=17.5*2*pi
z2=0.03
g2=-0.0005

a_th_mode2 = controls.TransferFunction([g2*w2**2,0,0],[1,2*z2*w2,w2**2])

a_th_tf = a_th_mode1 + a_th_mode2

import exp_data

f2 = exp_data.f2

exp_data.plot_exp()

a_th_tf.FreqResp(f2, fignum=3, clear=False)

wth1=2.8*2*pi
zth1=0.1
g_th=0.000003
p_th=10.0*2*pi
wth2=19.5*2*pi
zth2=0.07

num_th = g_th*p_th*wth1**2*wth2**2 * a_th_tf.den
p0 = poly1d([1,p_th,0])
p1 = poly1d([1,2*zth1*wth1,wth1**2])
p2 = poly1d([1,2*zth2*wth2,wth2**2])
den_th = p0*p1*p2

th_v_tf = controls.TransferFunction(num_th, den_th)

th_v_tf.FreqResp(f2, fignum=1, clear=False)

a_v_tf = a_th_tf*th_v_tf

a_v_tf.FreqResp(f2, fignum=2, clear=False)

show()

