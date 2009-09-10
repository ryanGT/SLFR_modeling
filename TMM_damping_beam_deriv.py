from sympy import *
import copy

var(['A1', 'A2', 'A3', 'A4', 'beta', 'x', 'L', 'EI'])
var(['d1','d2','d3','d4', 'a', 'c', 's'])

bxl = beta*x/L

#First, prove that I know what I am doing for an undamped bean
#(mainly proving I can do the derivation correctly in sympy)

W_y = A1*sin(bxl)+A2*cos(bxl)+A3*sinh(bxl)+A4*cosh(bxl)
th_z = diff(W_y, x)
M_z = EI*diff(th_z, x)
V_y = -diff(M_z, x)

z_list = [W_y, th_z, M_z, V_y]
A_list = [A1, A2, A3, A4]

def make_U(z_list):
    U = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])

    for i, term in enumerate(z_list):
        for j, curA in enumerate(A_list):
            U[i,j] = term.coeff(curA)
    return U

U = make_U(z_list)

U0 = U.subs(x,0)
UL = U.subs(x,L)
U0i = U0.inv()

Bz = UL*U0i

d = Wild('d')

def _c_sub_one_entry(entin):
    matches = [d/2*cos(beta)+d/2*cosh(beta), \
               d/2*sinh(beta)-d/2*sin(beta), \
               d/2*cosh(beta)-d/2*cos(beta), \
               d/2*sin(beta)+d/2*sinh(beta)]
    replacements = [d1, d2, d3, d4]
    for pat, rep in zip(matches, replacements):
        mydict = entin.match(pat)
        if mydict:
            return mydict[d]*rep
    return entin

    
def clean_cs(matin):
    matout = copy.copy(matin)
    nr = matin.rows
    nc = matin.cols
    for i in range(nr):
        for j in range(nc):
            matout[i,j] = _c_sub_one_entry(matout[i,j])
    return matout

Bz_clean1 = clean_cs(Bz)
Bz_clean2 = Bz_clean1.subs(EI, L**2/a)
## Bz_clean2 = Bz_clean1.subs(L**2/EI, a)
## Bz_clean3 = Bz_clean2.subs(EI/(L**2), 1/a)
## Bz_clean4 = Bz_clean3.subs(L/EI, a/L)

Bz_clean = Bz_clean2
#Bz_clean = Bz_clean1

#beta/L sign check
Wp = diff(W_y, x)
Wdp = diff(Wp, x)
W3p = diff(Wdp, x)
Wqp = diff(W3p, x)
Ltest = Wqp/(-W_y)
L1 = Ltest.subs({A2:0,A3:0,A4:0})
#yes, -beta**4/L**4 (the negative sign is correct)


#new derivation with damping
#W_y is unchanged (but the definition of beta changes)
#th_z is also unchanged
Md_z = (1+c*s)*EI*diff(th_z, x)
Vd_y = -diff(Md_z, x)

zd_list = [W_y, th_z, Md_z, Vd_y]
Ud = make_U(zd_list)
Ud2 = U.subs(EI,(1+c*s)*EI)
Ud0 = Ud2.subs(x,0)
UdL = Ud2.subs(x,L)
Ud0i = Ud0.inv()

Bdz = UdL*Ud0i
Bdz_clean1 = clean_cs(Bdz)
Bdz_clean2 = Bdz_clean1.subs(EI, L**2/a)
Bdz_clean = Bdz_clean2
