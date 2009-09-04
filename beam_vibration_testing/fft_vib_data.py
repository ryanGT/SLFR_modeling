from pylab import *
from scipy import *

import txt_data_processing

#filename = 'log_dec_hold_and_release_1_SLFR_RTP_OL_Test_uend=77.txt'
filename = 'log_dec_hold_and_release_2_SLFR_RTP_OL_Test_uend=77.txt'
#filename = 'ping_mid_beam_w_allen_wrench_1_SLFR_RTP_OL_Test_uend=96.txt'
#filename = 'ping_mid_beam_w_eraser_2_SLFR_RTP_OL_Test_uend=96.txt'
#filename = 'stab_mid_beam_w_eraser_3_SLFR_RTP_OL_Test_uend=96.txt'
data_file = txt_data_processing.Data_File(filename)
data_file.Time_Plot(legloc=7)

a = data_file.a
a = squeeze(a)
A_fft = fft(a)

t = data_file.t
dt = t[1]-t[0]
T = t.max()+dt

df = 1.0/T
fs = 1.0/dt
f = arange(0,fs, df)

figure(2)
clf()
plot(f, abs(A_fft))

figure(2)
clf()
plot(f, abs(A_fft))
xlim(0,30)

#find first 2 natural freqs
ind_50 = where(f > 50)[0].min()
A_search = abs(A_fft)[0:ind_50]
ind1 = argmax(A_search)
f1 = f[ind1]
print('f1 = %s' % f1)

ind_10 = where(f > 10)[0].min()
ind2 = argmax(A_search[ind_10:])+ind_10
f2 = f[ind2]
print('f2 = %s' % f2)

show()
