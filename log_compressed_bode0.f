      a_1 = 1/betabeam**2
      a_2 = 1/betabeam
      a_3 = 1/(csp*s+ksp)
      a_4 = 1/s
      a_5 = 1/(tauact+s)
      RESULT = gainbode0*s**2*(5.0E-1*a_2*c4beam*Kact*Lbeam*a_4*tauact*a
     1   _5+1.0E+0*c1beam*Kact*Laccel*a_4*tauact*a_5+nzbv[1,1]*(1.0E+0*(
     2   5.0E-1*a_2*c4beam*Lbeam*a_3-5.0E-1*abeam*a_1*c3beam)+Laccel*(1.
     3   0E+0*c1beam*a_3+5.0E-1*abeam*a_2*c4beam/Lbeam))+nzbv[2,1]*(5.0E
     4   -1*abeam*a_1*c3beam*Laccel-5.0E-1*abeam*c2beam*Lbeam/betabeam**
     5   3))
