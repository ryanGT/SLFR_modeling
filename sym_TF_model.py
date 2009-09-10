from pylab import *
from scipy import *

import symcontrols
import copy

num1 = 'g1*s**2*w1**2'
den1 = '(s**2+2*z1*w1*s+w1**2)'
dict1 = {'w1':2.5*2*pi,'z1':0.03,'g1':0.003}
mode1_symTF = symcontrols.SymTF(num1,den1)

num2 = 'g2*s**2*w2**2'
den2 = '(s**2+2*z2*w2*s+w2**2)'
dict2 = {'w2':17.5*2*pi,'z2':0.03,'g2':-0.0005}
mode2_symTF = symcontrols.SymTF(num2,den2)

dict3 = copy.copy(dict2)
dict3.update(dict1)

a_th_symTF = mode1_symTF + mode2_symTF

num_th = 'g_th*p_th*wth1**2*wth2**2*' + den1 + '*' + den2
den_th = 's*(s+p_th)*(s**2+2*zth1*wth1*s+wth1**2)*' + \
       '(s**2+2*zth2*wth2*s+wth2**2)'

dict_th = {'wth1':2.8*2*pi,'zth1':0.1,'g_th':0.000003,'p_th':10.0*2*pi, \
           'wth2':19.5*2*pi, 'zth2':0.07}
dict_th.update(dict1)
dict_th.update(dict2)

th_v_symTF = symcontrols.SymTF(num_th, den_th)

a_v_symTF = th_v_symTF * a_th_symTF

if __name__ == '__main__':
    print('hi')
    #load exp
    #import exp_data
    #f2 = exp_data.f2
    #exp_data.plot_exp()


    #a_th_tf.FreqResp(f2, fignum=3, clear=False)
    #th_v_tf.FreqResp(f2, fignum=1, clear=False)

    #a_v_tf = a_th_tf*th_v_tf

    #a_v_tf.FreqResp(f2, fignum=2, clear=False)

    #show()

