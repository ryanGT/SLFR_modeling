from pylab import *
#Need optparse
#mport Real_Time_Python
import SLFR_RTP
reload(SLFR_RTP)
import controls

import time

#batch of steps with multiple kp's and multiple amp's (or a list of 1
#if you prefer).  n steps per set
kp = 1.0
kplist = [1.0]
amplist = [25]

stopn = 20000
N = 3

mytest = SLFR_RTP.P_control_Test(kp)

base = 'swept_sine_no_sat_amp=%s_stopn=%s_test_%s'

for cur_kp in kplist:
    mytest.kp = cur_kp
    for amp in amplist:
        for n in range(N):
            mytest.Reset_Theta()
            mytest.Swept_Sine(amp=amp, stopn=stopn, maxf=20.0)
            mytest.Close_Serial()
            basename = base % (amp, stopn, n+1)
            mytest.Save(basename)
            time.sleep(2.0)
            #mytest.Step_Response(amp=0, stopn=stopn)

show()
