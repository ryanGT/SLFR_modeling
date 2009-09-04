from pylab import *
from scipy import *
from scipy import optimize

import sys
if not '..' in sys.path:
    sys.path.append('..')
    
import SLFR_TMM
reload(SLFR_TMM)

bo = SLFR_TMM.Beam_Only_Model()
bwa = SLFR_TMM.Beam_w_Accel_Mass()
tpb = SLFR_TMM.Two_Piece_Beam_w_Accel_Mass()
bwas = SLFR_TMM.Beam_Base_Spring_Accel_Mass(k=500.0)

f = arange(0.1,30,0.1)
w = 2*pi*f
s = 1.0j*w
det = [bo.EigError(item) for item in s]
det = array(det)

figure(10)
clf()
plot(f, det)
xlabel('Freq. (Hz)')

figure(11)
clf()
plot(w, det)
xlabel('$\\omega$ (rad/sec)')

w1 = bo.FindEig(2.0j*2*pi)[1]
f1 = w1/(2*pi)

w1a = bwa.FindEig(2.0j*2*pi)[1]
f1a = w1a/(2*pi)

w2 = bo.FindEig(2.0j*pi*22)[1]
f2 = w2/(2*pi)

w2a = bwa.FindEig(2.0j*pi*22)[1]
f2a = w2a/(2*pi)

print('f1 = %s' % f1)
print('f1a = %s' % f1a)
print('f2 = %s' % f2)
print('f2a = %s' % f2a)

from SLFR_TMM import EI, mu, L

f1e = 2.5
f2e = 17.95

L2 = 16.5*25.4/1000
#Fit to find EI and mu
def mymodel(c):
    EI = c[0]
    mu = c[1]
    myparams = {'EI':EI, 'mu':mu, 'L':L}
    mysys = SLFR_TMM.Beam_w_Accel_Mass(beamparams=myparams)
    return mysys

def find_eigs(c, mysys=None):
    if mysys is None:
        mysys = mymodel(c)
    w1 = mysys.FindEig(2.0j*2*pi)[1]
    f1 = w1/(2*pi)

    w2 = mysys.FindEig(2.0j*pi*18)[1]
    f2 = w2/(2*pi)
    return f1, f2

def mycost(c):
    f1, f2 = find_eigs(c)
    e = (f1e-f1)**2+(f2e-f2)**2
    return e

run_fit = 0

if run_fit:
    c_fit = optimize.fmin(mycost, [EI, mu])

    f1f, f2f = find_eigs(c_fit)

show()
