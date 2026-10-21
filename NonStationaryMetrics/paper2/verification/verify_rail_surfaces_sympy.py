"""verify_rail_surfaces_sympy.py

Independent SymPy route for Lemma lem:rail-surfaces (appendix app:rail-surfaces):
extrinsic geometry of surfaces swept by equatorial rail families in
Thakurta-Kerr, g = A(eta)^2 g_Kerr, with A an arbitrary positive function.

It rebuilds the equatorial Christoffel symbols of both metrics from the
definitions, checks the conformal connection formula, the seed obstruction, its
A^2 (raw covector) scaling, the vanishing polar component in 4D, and the radial
coefficients in the physical unit normal sqrt(D)(-1,q,0,0). It does not import
the Mathematica script verify_rail_surfaces.wls.

Run: python3 verify_rail_surfaces_sympy.py
"""
import sympy as s

r, M, a, q, qp = s.symbols("r M a q qp", real=True)
A, Eh = s.symbols("A Ehat", positive=True)
At = s.symbols("A_eta", real=True)
f = 1-2*M/r
de = r*r-2*M*r+a*a
k = -2*M*a/r
p = r*r+a*a+2*M*a*a/r
g = s.Matrix([[-f,0,k],[0,r*r/de,0],[k,0,p]])
gi = s.simplify(g.inv())
checks = []

def check(name, expr):
    assert s.simplify(s.trigsimp(expr)) == 0, name
    checks.append(name)

def deriv(expr,i,conformal):
    if i == 1:
        return s.diff(expr,r)
    if i == 0 and conformal:
        return s.diff(expr,A)*At
    return 0

def christ(metric,inverse,conformal):
    return [[[s.factor(sum(inverse[i,l]*(deriv(metric[l,j],z,conformal)
        +deriv(metric[l,z],j,conformal)-deriv(metric[j,z],l,conformal))/2
        for l in range(3))) for z in range(3)] for j in range(3)] for i in range(3)]

GK = christ(g,gi,False)
GT = christ(A*A*g,gi/(A*A),True)
ell = s.Matrix([At/A,0,0])
for i in range(3):
    for j in range(3):
        for z in range(3):
            expected = ((ell[z] if i == j else 0)
                +(ell[j] if i == z else 0)-g[j,z]*(gi*ell)[i])
            check(f"connection_{i}{j}{z}",GT[i][j][z]-GK[i][j][z]-expected)

X,Y,N = s.Matrix([q,1,0]),s.Matrix([0,0,1]),s.Matrix([-1,q,0])
h = s.Matrix([[(X.T*g*X)[0],(X.T*g*Y)[0]],
              [(Y.T*g*X)[0],(Y.T*g*Y)[0]]])
S = de*q*q/r**2-p/de
check("normal norm",(N.T*gi*N)[0]-S)
check("det h",h.det()+r*r*S)

def coeffs(G):
    out=[]
    for U,V,acc in [(X,X,[qp,0,0]),(X,Y,[0,0,0]),(Y,Y,[0,0,0])]:
        Z=s.Matrix([acc[i]+sum(G[i][j][z]*U[j]*V[z]
            for j in range(3) for z in range(3)) for i in range(3)])
        out.append(s.factor((N.T*Z)[0]))
    return out

cK,cT=coeffs(GK),coeffs(GT)
for i,(j,z) in enumerate([(0,0),(0,1),(1,1)]):
    check(f"raw conformal II_{j}{z}",cT[i]-cK[i]+(At/A)*(gi*N)[0]*h[j,z])
UK=s.factor(cK[1]*h[1,1]-cK[2]*h[0,1])
UT=s.factor(cT[1]*(A*A*h[1,1])-cT[2]*(A*A*h[0,1]))
check("seed obstruction",UK+M*a*(3*r*r+a*a)*S/r**2)
check("raw conformal obstruction",UT-A*A*UK)

# Full 4D polar check, differentiate before setting theta=pi/2.
th=s.symbols("theta", real=True)
ss=r*r+a*a*s.cos(th)**2
g4=s.Matrix([[-(1-2*M*r/ss),0,0,-2*M*a*r*s.sin(th)**2/ss],
    [0,ss/de,0,0],[0,0,ss,0],
    [-2*M*a*r*s.sin(th)**2/ss,0,0,
     (r*r+a*a+2*M*a*a*r*s.sin(th)**2/ss)*s.sin(th)**2]])
for i in [0,1,3]:
    for j in [0,1,3]:
        check(f"polar component {i}{j}",s.diff(A*A*g4[i,j],th).subs(th,s.pi/2))

# Physical unit-normal radial coefficients; n_cov=sqrt(D)*(-1,q,0).
ff=1-2*M/r
D=Eh*Eh-A*A*ff
qr=Eh/(ff*s.sqrt(D))
qprime=s.diff(qr,r)+s.diff(qr,A)*At*qr
subs={a:0,q:qr,qp:qprime}
rad=[s.factor(s.sqrt(D)*b.subs(subs)) for b in cT]
check("radial physical b_rr",rad[0]+A**3*At/D**s.Rational(3,2))
check("radial physical b_rs",rad[1])
check("radial physical b_ss",rad[2]+Eh*r+At*r*r*s.sqrt(D)/(A*ff))
check("physical h_rr",(A*A*h[0,0]).subs(subs)+A**4/D)
check("physical radial normal norm",(D*S/(A*A)).subs(subs)-1)
minor=s.factor(rad[0]*A*A*r*r-rad[2]*(-A**4/D))
expected=-A**3*Eh*r/D**s.Rational(3,2)*(A*s.sqrt(D)+At*Eh*r/ff)
check("radial umbilicity minor",minor-expected)
at_umb=-A*ff*s.sqrt(D)/(Eh*r)
check("umbilic b_ss nonzero coefficient",rad[2].subs(At,at_umb)+A*A*ff*r/Eh)
check("umbilic b_rr proportionality",rad[0].subs(At,at_umb)-A**4*ff/(Eh*r*D))

qin=-qr
qinprime=s.diff(qin,r)+s.diff(qin,A)*At*qin
radin=[s.factor(s.sqrt(D)*b.subs({a:0,q:qin,qp:qinprime})) for b in cT]
check("incoming b_rr",radin[0]+A**3*At/D**s.Rational(3,2))
check("incoming b_ss",radin[2]-Eh*r+At*r*r*s.sqrt(D)/(A*ff))

# Areal-radius derivative and admissibility inequality.
check("Ar condition",at_umb*r+A/qr)
check("profile derivative maximum",s.diff(ff/r**2,r)+2*(r-3*M)/r**4)
check("profile maximum value",(ff/r**2).subs(r,3*M)-1/(27*M*M))

# Physical differential speed inequality: norm and rail imply
# (dr/deta)^2 + f*r^2*angular_speed^2 = f^2*(1-A^2*f/Eh^2).
omega2=s.symbols("angular_speed_squared", nonnegative=True)
V=ff*s.sqrt(D)/Eh
check("radial speed squared",V**2-ff**2*(1-A*A*ff/Eh**2))

# Drift model A_eta/A = eps constant on the profile A r = c: squared umbilic condition
# r^5 [f^2 D - eps^2 Ehat^2 r^2] is a degree-7 polynomial in r, not identically zero.
eps,cc=s.symbols("eps c", positive=True)
Pc=s.Poly(s.expand(r**5*(ff**2*D.subs(A,cc/r)-eps**2*Eh**2*r**2)),r)
check("areal-radius drift polynomial degree 7",Pc.degree()-7)
check("areal-radius drift leading coefficient",Pc.coeff_monomial(r**7)+eps**2*Eh**2)
check("areal-radius drift r^5 coefficient",Pc.coeff_monomial(r**5)-Eh**2)

print(f"PASS {len(checks)}/{len(checks)} independent algebra checks")
print("Unit TK obstruction: -A**3*M*a*(3*r**2+a**2)*sqrt(S)/r**2")
print("Radial b_rr=-A**3*A_eta/D**(3/2)")
print("Radial b_ss=-Ehat*r-A_eta*r**2*sqrt(D)/(A*f)")
print("Radial umbilicity: A_eta=-A*f*sqrt(D)/(Ehat*r)")
