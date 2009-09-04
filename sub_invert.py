import sympy
sympy.var('a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34')
sympy.var('a41 a42 a43 a44')
sympy.var('b11 b12 b13 b14 b21 b22 b23 b24 b31 b32 b33 b34')
sympy.var('b41 b42 b43 b44')
A = sympy.Matrix([[a11, a12, a13, a14],[a21, a22, a23, a24],\
                  [a31, a32, a33, a34],[a41, a42, a43, a44]])
B = sympy.Matrix([[b11, b12, b13, b14],[b21, b22, b23, b24],\
                  [b31, b32, b33, b34],[b41, b42, b43, b44]])
C = A*B
submat = C[2:,2:]
submati = submat.inv()

D = B.inv()*A.inv()
submati2 = D[2:,2:]
