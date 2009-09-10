from pylab import *
from scipy import *

import sympy
import symcontrols
import copy

sympy.var(['g1','s','w1','z1','g2','z2','w2'])

num1 = g1*s**2*w1**2
den1 = (s**2+2*z1*w1*s+w1**2)

num2 = g2*s**2*w2**2
den2 = (s**2+2*z2*w2*s+w2**2)

num3 = num1*den2+num2*den1
den3 = den1*den2


sympy.var(['g_th','p_th','wth1','wth2','zth1','zth2'])
num_th = g_th*p_th*wth1**2*wth2**2 * den1 * den2
den_th = s*(s+p_th)*(s**2+2*zth1*wth1*s+wth1**2)*(s**2+2*zth2*wth2*s+wth2**2)

num4 = num3*num_th
den4 = den3*den_th

