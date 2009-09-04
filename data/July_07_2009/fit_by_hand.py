from pylab import *
from scipy import *
import copy
import txt_data_processing, rwkbode, pylab_util

from systemid import Model

#load log downsampled and raw Bode data
log_ds_mod = 'swept_sine_amp_75_July_07_2009_log_downsampled'
log_ds_data = txt_data_processing.load_avebode_data_set(log_ds_mod)
log_ds_f = log_ds_data.f
a_theta_log_ds = log_ds_data.find_bode('a', 'theta')

raw_mod = 'swept_sine_amp_75_July_07_2009_avebodes'
raw_data = txt_data_processing.load_avebode_data_set(raw_mod)
raw_f = raw_data.f
a_theta_raw = raw_data.find_bode('a', 'theta')

######################
#
# Develop a model
#
######################

#both modes multiplied together
tf_c1 = Model('g*s**2*w1**2*w2**2', \
           '(s**2+2*z1*w1*s+w1**2)*(s**2+2*z2*w2*s+w2**2)' , \
           {'w1':2.5*2*pi,'z1':0.03,'w2':17.8*2*pi,'z2':0.01,'g':0.005}, \
           'all')
model_bode_c1 = rwkbode.Bode_From_TF(tf_c1, raw_f, input='theta', output='a')

#adding two modes together with different gains
num1 = 'g1*s**2*w1**2'
den1 = '(s**2+2*z1*w1*s+w1**2)'
dict1 = {'w1':2.5*2*pi,'z1':0.03,'g1':0.003}
tf1 = Model(num1, \
            den1, \
            dict1, \
            'all')
model_bode_m1 = rwkbode.Bode_From_TF(tf1, raw_f, input='theta', output='a')

num2 = 'g2*s**2*w2**2'
den2 = '(s**2+2*z2*w2*s+w2**2)'
dict2 = {'w2':17.5*2*pi,'z2':0.03,'g2':-0.0005}
tf2 = Model(num2, \
            den2, \
            dict2, \
            'all')
model_bode_m2 = rwkbode.Bode_From_TF(tf2, raw_f, input='theta', output='a')

dict3 = copy.copy(dict2)
dict3.update(dict1)
num3 = num1 + '*' + den2 + '+' + num2 + '*' + den1
den3 = den1 + '*' + den2
tf3 = Model(num3, \
            den3, \
            dict3, \
            'all')
model_bode_c2 = rwkbode.Bode_From_TF(tf3, raw_f, input='theta', output='a')



#Plot Experimental and Model Bodes
a_theta_fi = 1
rwkbode.GenBodePlot(a_theta_fi, raw_f, a_theta_raw)
rwkbode.GenBodePlot(a_theta_fi, log_ds_f, a_theta_log_ds, clear=False, \
                    linetype='o')
rwkbode.GenBodePlot(a_theta_fi, raw_f, model_bode_c1, clear=False, \
                    linetype='k-')
rwkbode.GenBodePlot(a_theta_fi, raw_f, model_bode_m1, clear=False, \
                    linetype='-')
rwkbode.GenBodePlot(a_theta_fi, raw_f, model_bode_m2, clear=False, \
                    linetype='-')
rwkbode.GenBodePlot(a_theta_fi, raw_f, model_bode_c2, clear=False, \
                    linetype='-')
pylab_util.SetPhaseLim(a_theta_fi, [-200, 200])
pylab_util.SetMagLim(a_theta_fi, [-10, 45])

show()
