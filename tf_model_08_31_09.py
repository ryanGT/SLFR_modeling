from pylab import *
from scipy import *
import copy
import txt_data_processing, rwkbode, pylab_util, rwkos

from systemid import Model, PolyHasher

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

l1 = log10(0.5)
l2 = log10(50)
flist = arange(0.5, 1, 0.5).tolist()+f.tolist()+arange(26,30, 0.5).tolist()
f2 = array(flist)

## #load log downsampled and raw Bode data
## log_ds_mod = 'swept_sine_amp_75_July_07_2009_log_downsampled'
## log_ds_data = txt_data_processing.load_avebode_data_set(log_ds_mod)
## log_ds_f = log_ds_data.f
## a_theta_log_ds = log_ds_data.find_bode('a', 'theta')

## raw_mod = 'swept_sine_amp_75_July_07_2009_avebodes'
## raw_data = txt_data_processing.load_avebode_data_set(raw_mod)
## raw_f = raw_data.f
## a_theta_raw = raw_data.find_bode('a', 'theta')

######################
#
# Develop a model
#
######################


######################
#
#  a / theta model
#
######################

## #both modes multiplied together
## tf_c1 = Model('g*s**2*w1**2*w2**2', \
##            '(s**2+2*z1*w1*s+w1**2)*(s**2+2*z2*w2*s+w2**2)' , \
##            {'w1':2.5*2*pi,'z1':0.03,'w2':17.8*2*pi,'z2':0.01,'g':0.005}, \
##            'all')
## model_bode_c1 = rwkbode.Bode_From_TF(tf_c1, f2, input='theta', output='a')

#adding two modes together with different gains
num1 = 'g1*s**2*w1**2'
den1 = '(s**2+2*z1*w1*s+w1**2)'
dict1 = {'w1':2.5*2*pi,'z1':0.03,'g1':0.003}
tf1 = Model(num1, \
            den1, \
            dict1, \
            'all')
model_bode_m1 = rwkbode.Bode_From_TF(tf1, f2, input='theta', output='a')

num2 = 'g2*s**2*w2**2'
den2 = '(s**2+2*z2*w2*s+w2**2)'
dict2 = {'w2':17.5*2*pi,'z2':0.03,'g2':-0.0005}
tf2 = Model(num2, \
            den2, \
            dict2, \
            'all')
model_bode_m2 = rwkbode.Bode_From_TF(tf2, f2, input='theta', output='a')

dict3 = copy.copy(dict2)
dict3.update(dict1)
num3 = num1 + '*' + den2 + '+' + num2 + '*' + den1
den3 = den1 + '*' + den2
tf3 = Model(num3, \
            den3, \
            dict3, \
            'all')
model_bode_c2 = rwkbode.Bode_From_TF(tf3, f2, input='theta', output='a')



#Plot Experimental and Model Bodes
a_theta_fi = 3
#rwkbode.GenBodePlot(a_theta_fi, f2, a_theta_raw)
#rwkbode.GenBodePlot(a_theta_fi, log_ds_f, a_theta_log_ds, clear=False, \
#                    linetype='o')
## rwkbode.GenBodePlot(a_theta_fi, f2, model_bode_c1, clear=False, \
##                     linetype='k-')
## rwkbode.GenBodePlot(a_theta_fi, f2, model_bode_m1, clear=False, \
##                     linetype='-')
## rwkbode.GenBodePlot(a_theta_fi, f2, model_bode_m2, clear=False, \
##                     linetype='-')
rwkbode.GenBodePlot(a_theta_fi, f2, model_bode_c2, clear=False, \
                    linetype='-')
pylab_util.SetPhaseLim(a_theta_fi, [-200, 200])
pylab_util.SetMagLim(a_theta_fi, [-10, 45])
pylab_util.SetFreqLim(a_theta_fi, [0.5, 30])


######################
#
# theta / v model
#
######################
num_th = 'g_th*p_th*wth1**2*wth2**2*' + den1 + '*' + den2
den_th = 's*(s+p_th)*(s**2+2*zth1*wth1*s+wth1**2)*' + \
       '(s**2+2*zth2*wth2*s+wth2**2)'

dict_th = {'wth1':2.8*2*pi,'zth1':0.1,'g_th':0.000003,'p_th':10.0*2*pi, \
           'wth2':19.5*2*pi, 'zth2':0.07}
dict_th.update(dict1)
dict_th.update(dict2)
act_tf = Model(num_th, \
               den_th, \
               dict_th, \
               'all')
act_bode = rwkbode.Bode_From_TF(act_tf, f2, input='v', output='theta')
act_bode_fi = 1
rwkbode.GenBodePlot(act_bode_fi, f2, act_bode, clear=False, \
                    linetype='-')




######################
#
# a / v model
#
#####################
a_v_fi = 2
a_v_tf = tf3*act_tf

a_v_bode = rwkbode.Bode_From_TF(a_v_tf, f2, input='v', output='a')
rwkbode.GenBodePlot(a_v_fi, f2, a_v_bode, clear=False, \
                    linetype='-')

#sympy verification
import sympy
sympy.var(['g_th','p_th','wth1','wth2','s','z1','w1','z2','w2','g1','g2', \
           'zth1','zth2'])

Nth = g_th*p_th*wth1**2*wth2**2*(s**2+2*z1*w1*s+w1**2)*(s**2+2*z2*w2*s+w2**2)
Dth = s*(s+p_th)*(s**2+2*zth1*wth1*s+wth1**2)* + \
      (s**2+2*zth2*wth2*s+wth2**2)

N1 = g1*s**2*w1**2
D1 = (s**2+2*z1*w1*s+w1**2)
N2 = g2*s**2*w2**2
D2 = (s**2+2*z2*w2*s+w2**2)

TH_V_test = Nth/Dth
A_TH_test = (N1*D2+N2*D1)/(D1*D2)

A_V_test = A_TH_test*TH_V_test
Nav, Dav = A_V_test.as_numer_denom()

#exec('Dtest = '+a_v_tf.den_str)
#exec('Ntest = '+a_v_tf.num_str)
Dtest = sympy.sympify(a_v_tf.den_str)
Ntest = sympy.sympify(a_v_tf.num_str)

Dsub = Dtest.subs(dict_th)
Dpoly = sympy.Poly(Dsub, s)
Dnumpy = poly1d(Dpoly.coeffs)
Dpoles = roots(Dnumpy)

new_den_str = '((s**2+2*z1*w1*s+w1**2)*(s**2+2*z2*w2*s+w2**2))*(s*(s+p_th)*(s**2+2*zth1*wth1*s+wth1**2)*(s**2+2*zth2*wth2*s+wth2**2))'
new_num_str = '(g1*s**2*w1**2*(s**2+2*z2*w2*s+w2**2)+g2*s**2*w2**2*(s**2+2*z1*w1*s+w1**2))*(g_th*p_th*wth1**2*wth2**2*(s**2+2*z1*w1*s+w1**2)*(s**2+2*z2*w2*s+w2**2))'
new_var_dict = {'zth2': 0.070000000000000007, 'g2': -0.00050000000000000001, 'g1': 0.0030000000000000001, 'zth1': 0.10000000000000001, 'w2': 109.95574287564276, 'w1': 15.707963267948966, 'z1': 0.029999999999999999, 'z2': 0.029999999999999999, 'wth1': 17.592918860102841, 'p_th': 62.831853071795862, 'wth2': 122.52211349000193, 'g_th': 3.0000000000000001e-06}
new_opt_dict = {'zth2': 0.070000000000000007, 'g2': -0.00050000000000000001, 'g1': 0.0030000000000000001, 'zth1': 0.10000000000000001, 'w2': 109.95574287564276, 'w1': 15.707963267948966, 'z1': 0.029999999999999999, 'z2': 0.029999999999999999, 'wth1': 17.592918860102841, 'p_th': 62.831853071795862, 'wth2': 122.52211349000193, 'g_th': 3.0000000000000001e-06}

test = Model(new_num_str,new_den_str,new_var_dict,new_opt_dict,\
             s)

Dlist, junk = PolyHasher(a_v_tf.den_str, dict_th)

def root_checker(poly_str, var_dict=dict_th):
    coef_list, junk = PolyHasher(poly_str, var_dict)
    my_poly = poly1d(coef_list)
    return roots(my_poly)

    
show()
