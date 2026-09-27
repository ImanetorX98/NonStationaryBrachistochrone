#!/usr/bin/env python3
"""Proposition (where the tau half-angle decreases) and its appendix: exact and
interval checks.

    python3 verify_lemma_B.py            # symbolic certificates + interval signs
    python3 verify_lemma_B.py numeric    # 40- and 70-digit evaluations (not proofs)

Claims covered (paper2.tex, Proposition prop:lemmaB and Appendix app:lemmaB), M = 1,
b = E^2 - 1, c = a^2:
  (v)   g = V W positive, increasing, concave for b >= 1/2 and 0 <= c <= 1
        (monomial lowering in c; positive coefficients at r = 2+z, b = 1/2+u);
  (i)   d log(sqrt(V) W)/dr < 0 at a = 0;
  (ii)  alpha - 2/3 > 0 at a = 0 for r >= 9/4;
  (iii) alpha - 1/2 > 0 for r >= 5/2, 0 <= c <= 1;
  (iv)  alpha - 2/3 > 0 for r >= 3,   0 <= c <= 1;
  the elementary inequalities acosh(2) < 4/3 and 8/9 < 4^(1/3)/sqrt(3);
  the lowered A_c polynomial printed in the appendix;
  the five interval-certified derivative signs of the two counterexamples.

The proofs and this script were first written by GPT-6 Astra in an audit of the
proofs (27 September 2026) and checked independently before inclusion.  Every
assertion raises, so any failure gives a nonzero exit code.
"""
import sys
import sympy as s
import mpmath as mp
from mpmath import iv
from scipy.optimize import brentq

def symbolic():
    r,b,c,z,u=s.symbols('r b c z u',positive=True)
    D=b*r+2; N=b*r*r+(3-b)*r-4
    V=r*r*(r-2)/D
    W=D**2/(2*r*N*s.sqrt(r*(r-2)))
    f=D*s.sqrt(r*(r-2))/(2*N)
    assert s.simplify(f-V*W)==0
    fp=s.factor(s.diff(f,r)/s.diff(V,r))
    fpp=s.factor(s.diff(fp,r)/s.diff(V,r))
    A=b*b*r*r+b*r**3-5*b*r*r+10*b*r-2*r+8
    assert s.simplify(fp-D**2*A/(4*r**s.Rational(3,2)*s.sqrt(r-2)*N**3))==0
    P=s.factor(-fpp*8*r**s.Rational(7,2)*(r-2)**s.Rational(3,2)*N**5/D**3)
    assert s.Poly(P,r,b)
    for name,p in [('A',A),('P',P)]:
        shifted=s.Poly(s.expand(p.subs({r:z+2,b:u+s.Rational(1,2)})),z,u)
        assert all(co>0 for co in shifted.coeffs())
        print('EXACT PASS: positive coefficient certificate for',name,'at r=2+z,b=1/2+u')
        print('coefficients by z:',[s.factor(shifted.as_expr().coeff(z,i)) for i in range(shifted.degree(z)+1)])
    q=s.factor(s.sqrt(V)*W)
    logq=s.factor(s.diff(W,r)/W+s.diff(V,r)/(2*V))
    expected=-(b*b*r**3+5*b*r*r+b*r+9*r-4)/(r*D*N)
    assert s.simplify(logq-expected)==0
    print('EXACT PASS: d log(sqrt(V) W)/dr =',expected)
    alpha=s.factor(-V*s.diff(W,r)/(s.diff(V,r)*W))
    num,den=s.fraction(s.factor(alpha-s.Rational(2,3)))
    shifted=s.Poly(s.expand(num.subs(r,z+s.Rational(9,4))),z,b)
    assert all(co>0 for co in shifted.coeffs())
    print('EXACT PASS: static alpha>2/3 for r>=9/4; denominator',den)
    print('positive numerator:',shifted.as_expr())
    Delta=r*r-2*r+c; T=r*N+c; Vr=r*Delta/D
    fr=r*D*s.sqrt(r*(r-2))/(2*T)
    fpr=s.factor(s.diff(fr,r)/s.diff(Vr,r))
    fppr=s.factor(s.diff(fpr,r)/s.diff(Vr,r))
    Ar=s.factor(fpr*4*s.sqrt(r-2)*T**3/(s.sqrt(r)*D**2))
    Pr=s.factor(-fppr*8*s.sqrt(r)*(r-2)**s.Rational(3,2)*T**5/D**3)
    for name,pol in [('Arot',Ar),('Prot',Pr)]:
        shifted=s.Poly(s.expand(pol.subs({r:2+z,b:s.Rational(1,2)+u})),z,u,c)
        lower=s.Poly(sum(co*z**i*u**j for (i,j,k),co in shifted.terms() if k==0 or co<0),z,u)
        assert all(co>0 for co in lower.coeffs())
        print('EXACT PASS: rotating f positive increasing concave, certificate',name)
        print('numerator before shift:',pol)
        print('lower coefficients by z:',[s.factor(lower.as_expr().coeff(z,i)) for i in range(lower.degree(z)+1)])
    logWr=1/(2*r)+1/(2*(r-2))-s.diff(Delta,r)/Delta+2*b/D-s.diff(T,r)/T
    ar=s.factor(-Vr*logWr/s.diff(Vr,r))
    for power,rad in [(s.Rational(1,2),s.Rational(5,2)),(s.Rational(2,3),s.Integer(3))]:
        n,d=s.fraction(s.factor(ar-power))
        poly=s.Poly(s.expand(n.subs(r,z+rad)),z,b,c)
        # On 0<=c<=1: drop positive terms carrying c; bound negative ones by c=1.
        lower=s.Poly(sum(co*z**i*b**j for (i,j,k),co in poly.terms() if k==0 or co<0),z,b)
        assert all(co>0 for co in lower.coeffs())
        print('EXACT PASS: rotating alpha>',power,'for r>=',rad,'0<=c<=1')
        print('denominator:',d)
        print('lower polynomial:',lower.as_expr())
    # Elementary constant used in the 2/3-exponent quarter criterion.
    assert 1+s.Rational(4,3)**2/2+s.Rational(4,3)**4/24>2
    assert (s.Rational(8,9)**2*3)**3<16
    print('EXACT PASS: acosh(2)<4/3 and 8/9<4**(1/3)/sqrt(3)')

def inverse(t,b):
    guess=2+max(mp.sqrt(b*t),mp.root(2*t,3))
    rr=mp.findroot(lambda r:r*r*(r-2)-t*(b*r+2), (guess,guess*mp.mpf('1.01')))
    assert rr>2
    return rr

def high_precision():
    for digits in [40,70]:
        mp.mp.dps=digits
        print('HIGH PRECISION (not interval), digits=',digits)
        def calc(b,r0,rt=None,frac=None):
            V=lambda r:r*r*(r-2)/(b*r+2)
            v0=V(r0); x=V(rt) if rt is not None else v0*frac
            if rt is None: rt=inverse(x,b)
            U=mp.acosh(mp.sqrt(v0/x))
            def f(r):return (b*r+2)*mp.sqrt(r*(r-2))/(2*(b*r*r+(3-b)*r-4))
            def fp(r):
                D=b*r+2;N=b*r*r+(3-b)*r-4
                A=b*b*r*r+b*r**3-5*b*r*r+10*b*r-2*r+8
                return D*D*A/(4*r**mp.mpf('1.5')*mp.sqrt(r-2)*N**3)
            integ=mp.quad(lambda u:fp(inverse(x*mp.cosh(u)**2,b))*mp.cosh(u),[0,U/4,U/2,3*U/4,U])
            xp=2*x*integ-f(r0)*mp.sqrt(x/(v0-x))
            print('b=',b,'r0=',r0,'rt=',mp.nstr(rt,30),'x Phi_prime=',mp.nstr(xp,35),'Phi_prime=',mp.nstr(xp/x,35))
        calc(mp.mpf(1),mp.mpf(201)/100,frac=mp.mpf(1)/4)
        for rt in [3,6,30,9000]:calc(mp.mpf(1)/100,mp.mpf(10000),mp.mpf(rt))
        Ustar=mp.findroot(lambda u:u*mp.tanh(u)-1,mp.mpf('1.2'))
        print('Ustar=',mp.nstr(Ustar,35),'cstar=',mp.nstr(1/mp.cosh(Ustar)**2,35))

def interval_checks():
    iv.dps=40
    def V(r,b):return r*r*(r-2)/(b*r+2)
    def f(r,b):return (b*r+2)*iv.sqrt(r*(r-2))/(2*(b*r*r+(3-b)*r-4))
    def fp(r,b):
        D=b*r+2;N=b*r*r+(3-b)*r-4
        A=b*b*r*r+b*r**3-5*b*r*r+10*b*r-2*r+8
        return D*D*A/(4*r*iv.sqrt(r)*iv.sqrt(r-2)*N**3)
    def inv(t,b):
        # Floating root ONLY proposes brackets. The two inequalities below certify them.
        bf=float(b.mid); low=float(t.a); high=float(t.b)
        def guess(tv):
            top=3+(bf*tv)**.5+(2*tv)**(1/3)
            while top*top*(top-2)/(bf*top+2)<tv:
                top*=2
            return brentq(lambda r:r*r*(r-2)/(bf*r+2)-tv,2,top,xtol=1e-13)
        rl=guess(low);rh=guess(high)
        for j in range(10):
            margin=1e-12*10**j*max(1,abs(rl),abs(rh))
            lo=max(2,rl-margin); hi=rh+margin
            if lo>2 and V(iv.mpf(lo),b).b<=t.a and V(iv.mpf(hi),b).a>=t.b:
                return iv.mpf([lo,hi])
        raise RuntimeError('Uncertified inverse bracket')
    def bounds(b,r0,rt=None,frac=None,N=6000):
        v0=V(r0,b);x=V(rt,b) if rt is not None else v0*frac
        ratio=v0/x; U=iv.ln(iv.sqrt(ratio)+iv.sqrt(ratio-1))
        total=iv.mpf(0)
        for k in range(N):
            z=iv.mpf([(iv.mpf(k)/N).a,(iv.mpf(k+1)/N).b])
            uz=U*z; co=(iv.exp(uz)+iv.exp(-uz))/2
            rr=inv(x*co*co,b)
            total+=fp(rr,b)*co/N
        ans=2*x*U*total-f(r0,b)*iv.sqrt(x/(v0-x))
        print('INTERVAL x Phi_prime:',ans,'N=',N,flush=True)
        return ans
    print('INTERVAL COUNTEREXAMPLE QUARTER: b=1,r0=201/100,x=V0/4',flush=True)
    out=bounds(iv.mpf(1),iv.mpf(201)/100,frac=iv.mpf(1)/4)
    assert out.a>0
    for rt,sign in [(3,1),(6,-1),(30,1),(9000,-1)]:
        print('INTERVAL FULL LEMMA: b=1/100,r0=10000,rt=',rt,flush=True)
        out=bounds(iv.mpf(1)/100,iv.mpf(10000),iv.mpf(rt))
        assert out.a>0 if sign>0 else out.b<0
    print('ALL FIVE DERIVATIVE SIGNS INTERVAL CERTIFIED')

def independent_integral():
    """Different formula: differentiate the original radial turning integral.

    The singularity is removed by r=rt+u**2. This does not use f', W',
    the derivative identity or the potential-coordinate quadrature above.
    """
    mp.mp.dps=55
    def phi(x,b,R):
        rt=inverse(x,b);J=mp.sqrt(x);U=mp.sqrt(R-rt)
        def integrand(u):
            rr=rt+u*u;D=b*rr+2
            Q=rr*rr+(rt-2)*rr+rt*(rt-2)-x*b
            return 2*J*mp.sqrt(D/(rr*(rr-2)*Q))
        nodes=sorted(set([mp.mpf(0),U]+[mp.mpf(t) for t in ['0.1','1','10'] if mp.mpf(t)<U]))
        return mp.quad(integrand,nodes)
    cases=[(mp.mpf(1),mp.mpf(201)/100,None)]
    cases += [(mp.mpf(1)/100,mp.mpf(10000),mp.mpf(rt)) for rt in [3,6,30,9000]]
    for b,R,rt in cases:
        V=lambda r:r*r*(r-2)/(b*r+2)
        x=V(R)/4 if rt is None else V(rt)
        print('INDEPENDENT RADIAL INTEGRAL b=',b,'R=',R,'rt=',rt)
        for eps in ['1e-5','1e-7','1e-9']:
            h=x*mp.mpf(eps)
            out=x*(phi(x+h,b,R)-phi(x-h,b,R))/(2*h)
            print('relative step=',eps,'x Phi_prime=',mp.nstr(out,35))

def rotating_example():
    mp.mp.dps=65
    b=mp.mpf(1);a=mp.mpf(1)/10000;c=a*a;R=mp.mpf(201)/100
    V=lambda r:r*(r*r-2*r+c)/(b*r+2)
    v0=V(R);x=v0/4
    def inv(t):
        return mp.findroot(lambda r:r*(r*r-2*r+c)-t*(b*r+2),(mp.mpf(2),R))
    def f(r):
        N=b*r*r+(3-b)*r-4;T=r*N+c
        return r*(b*r+2)*mp.sqrt(r*(r-2))/(2*T)
    def fp(r):
        D=b*r+2;N=b*r*r+(3-b)*r-4;T=r*N+c
        A=b*b*r**3+3*b*c*r*r-5*b*c*r+b*r**4-5*b*r**3+10*b*r*r+4*c*r-6*c-2*r*r+8*r
        return mp.sqrt(r)*D*D*A/(4*mp.sqrt(r-2)*T**3)
    U=mp.acosh(mp.sqrt(v0/x))
    xp=2*x*mp.quad(lambda u:fp(inv(x*mp.cosh(u)**2))*mp.cosh(u),[0,U])-f(R)*mp.sqrt(x/(v0-x))
    print('ROTATING high precision b=1,a=1/10000,R=201/100,x=V0/4')
    print('rt=',mp.nstr(inv(x),40),'x Phi_prime=',mp.nstr(xp,40))
    print('This example is below r=5/2. No claim that it refutes a theorem restricted to r>=5/2.')


def appendix_Ac():
    r,b,c,z,u=s.symbols('r b c z u',positive=True)
    Ac=r*(b*b*r*r+b*r**3-5*b*r*r+10*b*r-2*r+8)+c*(3*b*r*r-5*b*r+4*r-6)
    P=s.Poly(s.expand(Ac.subs({r:2+z,b:s.Rational(1,2)+u})),z,u,c)
    low=s.expand(sum(co*z**i*u**j for (i,j,k),co in P.terms() if k==0 or co<0))
    printed=(2*(2*u+3)**2+3*(2*u+1)*(2*u+3)*z+s.Rational(1,2)*(2*u+3)*(6*u+1)*z**2
             +s.Rational(1,4)*(2*u+1)*(2*u+7)*z**3+s.Rational(1,2)*(2*u+1)*z**4)
    assert s.expand(low-printed)==0
    print('EXACT PASS: lowered A_c equals the polynomial printed in the appendix')
    # printed P_{2/3} (static): alpha - 2/3 = P23/(6 N^2)
    D=b*r+2; N=b*r*r+(3-b)*r-4; V=r*r*(r-2)/D; W=D**2/(2*r*N*s.sqrt(r*(r-2)))
    alpha=-V*s.diff(W,r)/(s.diff(V,r)*W)
    P23=(2*b**2*r**4-4*b**2*r**3-b**2*r**2+9*b*r**3-13*b*r**2-14*b*r+18*r**2-42*r+8)
    assert s.simplify(alpha-s.Rational(2,3)-P23/(6*N**2))==0
    print('EXACT PASS: printed P_{2/3} equals 6 N^2 (alpha - 2/3)')
    # boundary constant C_c of g = C_c sqrt(V - v_h) at r -> 2 (high precision, not a proof)
    mp.mp.dps=40
    for bb,cc in [(mp.mpf('0.5'),mp.mpf('0.81')),(mp.mpf(2),mp.mpf(1)),(mp.mpf('0.7'),mp.mpf(0))]:
        rr=2+mp.mpf('1e-16'); Vr=rr*(rr*rr-2*rr+cc)/(bb*rr+2); vh=cc/(bb+1)
        gr=rr*(bb*rr+2)*mp.sqrt(rr*(rr-2))/(2*(rr*(bb*rr*rr+(3-bb)*rr-4)+cc))
        Cc=4*(bb+1)**2/(4*(bb+1)+cc)**mp.mpf('1.5')
        assert abs(gr/mp.sqrt(Vr-vh)/Cc-1)<mp.mpf('1e-12')
    print('NUMERIC PASS: g/sqrt(V-v_h) -> C_c = 4(b+1)^2/[4(b+1)+c]^(3/2) at three (b,c)')

if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'default'
    if mode=='numeric':
        high_precision()
    else:
        symbolic(); appendix_Ac(); interval_checks()
        print('ALL CHECKS PASSED')
