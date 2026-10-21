"""Exact ZAMO form and reproducible energy/type-IV boundaries for rotating TK.

G_hat = matrix/(A^2 N^2); T_hat = G_hat/(8 pi), G=c=1.
For arbitrary A use H=A'/A and Hp=dH/deta, not physical Hcos.
The chosen numeric cosmological example is local flat matter+Lambda at A0=1.
"""
from pathlib import Path
import json
import mpmath as mp
import sympy as sp
from audit_rotating_geometry import data
mp.mp.dps=80
D=mp.mpf

# Characteristic discriminant kept symbolic so a tiny spin shift is not lost.
x,b,Q,S,C=sp.symbols('x b Q S C',real=True)
poly=x**4+b*x**3+(Q-S)*x**2-b*S*x-(Q*S-C*C)
disc_expr=sp.discriminant(poly,x)
disc=sp.lambdify((b,Q,S,C),disc_expr,'mpmath')

def components(r,u,a,H,Hp):
    sig,de,B,N,om=data(r,u,a)
    sr=mp.sqrt(de/sig);sth=mp.sqrt(1-u*u)
    Nr=mp.diff(lambda z:data(z,u,a)[3],r)
    Nu=mp.diff(lambda z:data(r,z,a)[3],u)
    wr=mp.diff(lambda z:data(z,u,a)[4],r)
    wu=mp.diff(lambda z:data(r,z,a)[4],u)
    rootq=mp.sqrt(B/sig)*sth
    q=mp.matrix([2*H*sr*Nr,-2*H*sth*Nu/mp.sqrt(sig)])
    s=mp.matrix([H*sr*rootq*wr,-H*rootq*sth*wu/mp.sqrt(sig)])
    rho=3*H*H;p=-(2*Hp+H*H);beta=rho+p
    G=mp.matrix([[rho,q[0],q[1],0],[q[0],p,0,s[0]],
                  [q[1],0,p,s[1]],[0,s[0],s[1],p]])
    mixed=mp.diag([-1,1,1,1])*G
    return rho,p,beta,q,s,mixed,N

def nec_margin(r,u,a,H,Hp):
    rho,p,beta,q,s,mat,N=components(r,u,a,H,Hp)
    q2=mp.fdot(q,q);s2=mp.fdot(s,s);qs=mp.fdot(q,s)
    # maximize (1-z^2)|q+zs|^2; solve its cubic derivative.
    if not q2 and not s2: return beta
    if not s2: maximum=q2
    else:
        norm=max(q2,s2,abs(qs))
        coeff=[-4*s2/norm,-6*qs/norm,2*(s2-q2)/norm,2*qs/norm]
        import numpy as np
        candidates=[D(-1),D(1),D(0)]
        for root in np.roots([float(c) for c in coeff]):
            if abs(root.imag)<1e-8 and abs(root.real)<=1:
                fun=lambda z:(-4*s2*z**3-6*qs*z*z+2*(s2-q2)*z+2*qs)/norm
                z=mp.findroot(fun,D(str(root.real)),tol=D('1e-65'))
                if -1<=z<=1:candidates.append(z)
        maximum=max((1-z*z)*(q2+2*qs*z+s2*z*z) for z in candidates)
    return beta-2*mp.sqrt(maximum)

def type_disc(r,u,a,H,Hp):
    rho,p,beta,q,s,mat,N=components(r,u,a,H,Hp)
    # Scale x by beta, removing H^8 from the quartic discriminant.
    assert beta>0
    q2=mp.fdot(q,q)/beta**2;s2=mp.fdot(s,s)/beta**2;qs=mp.fdot(q,s)/beta**2
    # On the axis s=0, the two repeated transverse pressure eigenvalues
    # make the quartic discriminant identically zero at EVERY radius.
    if not s2:return 1-4*q2
    # Remove small transverse-eigenvalue factors before root tolerances.
    return disc(1,q2,s2,qs)/s2**2

def spectrum(r,u,a,H,Hp):
    rho,p,beta,q,s,mat,N=components(r,u,a,H,Hp)
    scale=max(abs(v) for v in mat)
    if not scale:return 'vacuum',[]
    values=mp.eig(mat/scale,left=False,right=False)
    nonreal=sum(abs(mp.im(v))>D('1e-50') for v in values)
    return ('IV' if nonreal else 'real-spectrum; inspect degeneracy/eigenvectors'),values

def independent_frame_check():
    # Compute coordinate Christoffels from the Kerr metric, then contract G_ab.
    # This checks the compact component formula, including both spatial shears.
    rr,tt,aa=sp.symbols('r theta a',real=True)
    sig=rr**2+aa**2*sp.cos(tt)**2;de=rr**2-2*rr+aa**2
    B=(rr**2+aa**2)**2-aa**2*de*sp.sin(tt)**2
    g=sp.Matrix([[-1+2*rr/sig,0,0,-2*aa*rr*sp.sin(tt)**2/sig],
        [0,sig/de,0,0],[0,0,sig,0],
        [-2*aa*rr*sp.sin(tt)**2/sig,0,0,B*sp.sin(tt)**2/sig]])
    gi=g.inv();coords=[sp.Symbol('eta'),rr,tt,sp.Symbol('phi')]
    gamma=sp.Matrix(4,4,lambda i,j:sum(gi[0,k]*(sp.diff(g[k,j],coords[i])+sp.diff(g[k,i],coords[j])-sp.diff(g[i,j],coords[k]))/2 for k in range(4)))
    fn=sp.lambdify((rr,tt,aa),(g,gi,gamma),'mpmath')
    for r,u,a in [(D(6),D('.3'),D('.8')),(D('1.8'),D('.4'),D('.9'))]:
        g,gi,gamma=fn(r,mp.acos(u),a);H=D('.02');Hp=D('-.0003')
        G=2*H*gamma+g*gi[0,0]*(2*Hp+H*H)
        G[0,0]+=2*(H*H-Hp)
        sig,de,B,N,om=data(r,u,a)
        e=mp.matrix([[1/N,0,0,om/N],[0,mp.sqrt(de/sig),0,0],
             [0,0,1/mp.sqrt(sig),0],[0,0,0,mp.sqrt(sig/(B*(1-u*u)))]])
        assert mp.norm(e*g*e.T-mp.diag([-1,1,1,1]))<D('1e-65')
        rho,p,beta,q,s,mixed,_=components(r,u,a,H,Hp)
        target=mp.diag([-1,1,1,1])*mixed/N**2
        assert mp.norm(e*G*e.T-target)<D('1e-65')

if __name__=='__main__':
    independent_frame_check()
    old=json.loads(Path(__file__).with_name('observational_scales.json').read_text())
    H=D(next(z for z in old['scales'] if z['name']=='Gaia BH3')['epsilon'])
    rg=D(next(z for z in old['scales'] if z['name']=='Gaia BH3')['rg_km'])*1000
    au=D(old['constants']['AU'])
    rows=[]
    for k,name in [(D('-.5'),'power law alpha=2'),(D('.5275'),'local flat LCDM Omega_m=.315')]:
        Hp=k*H*H;rst=mp.sqrt(2/((1-k)*H))
        for u in [D(1),1/mp.sqrt(2),D(0)]:
            a=D('.8')
            # Tiny shifts are resolved by secant shooting at 80-digit precision.
            root=mp.findroot(lambda y:type_disc(rst*y,u,a,H,Hp),(D('.999'),D('1.001')),tol=D('1e-65'),maxsteps=100)
            nec=mp.findroot(lambda y:nec_margin(rst*y,u,a,H,Hp)/H**2,(D('.999'),D('1.001')),tol=D('1e-65'))
            rows.append(dict(model=name,cos_theta=u,spin=a,Rstar_spherical_au=rst*rg/au,
              rotating_type_boundary_r_over_M=rst*root,
              relative_type_boundary_shift=root-1,relative_NEC_boundary_shift=nec-1))
        for r in [D('2.1'),D(6),rst/2,rst*2]:
            for u in [D(0),D('.5'),D(1)]:
                typ,eig=spectrum(r,u,D('.8'),H,Hp)
                margin=nec_margin(r,u,D('.8'),H,Hp)
                if r<rst/2 or r==rst/2:assert typ=='IV' and margin<0
                else:assert margin>0 and typ!='IV'
    # Spherical Harada criterion, independent numeric and exact symbolic checks.
    al,eta,R=sp.symbols('alpha eta r',positive=True)
    beta=2*al*(al+1)/eta**2;q=2*al/(eta*R**2)
    assert sp.simplify(beta-2*q-2*al/eta**2*(al+1-2*eta/R**2))==0
    solar=mp.sqrt(D(4)/3*D(old['constants']['rg_sun_km'])*1000*D(299792458)/D(old['constants']['H0_s_inverse']))/au
    out=dict(precision=80,Gaia_BH3_H0_tg=H,rows=rows,Rstar_1Msun_alpha2_au=solar,
      checks='PASS frame orthonormality; coordinate G vs compact ZAMO matrix; Harada reduction; inside/outside spectra and NEC')
    Path(__file__).with_name('rotating_energy.json').write_text(json.dumps(out,indent=2,default=lambda v:mp.nstr(v,55)))
    print(out['checks']);print('alpha2 1Msun AU',mp.nstr(solar,15))
    for z in rows:print(z['model'],'cos(theta)',mp.nstr(z['cos_theta'],5),'R*/AU',mp.nstr(z['Rstar_spherical_au'],15),'relative rotating shift',mp.nstr(z['relative_type_boundary_shift'],12),'NEC shift',mp.nstr(z['relative_NEC_boundary_shift'],12))
