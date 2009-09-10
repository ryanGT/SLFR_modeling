import sympy

sympy.var(['g1','s','w1','z1','g2','z2','w2'])

num1 = sympy.Poly(g1*s**2*w1**2,s)
den1 = sympy.Poly(s**2+2*z1*w1*s+w1**2,s)

num2 = sympy.Poly(g2*s**2*w2**2,s)
den2 = sympy.Poly(s**2+2*z2*w2*s+w2**2,s)

num3 = num1*den2+num2*den1
den3 = den1*den2


sympy.var(['g_th','p_th','wth1','wth2','zth1','zth2'])
num_th = sympy.Poly(g_th*p_th*wth1**2*wth2**2,s) * den1 * den2
den_th = sympy.Poly(s*(s+p_th),s)*sympy.Poly(s**2+2*zth1*wth1*s+wth1**2,s)* \
         sympy.Poly(s**2+2*zth2*wth2*s+wth2**2,s)

num4 = num3*num_th
den4 = den3*den_th

