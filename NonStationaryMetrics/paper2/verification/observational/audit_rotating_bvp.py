"""Frozen equatorial rail BVP, M=A=1, Ehat=7/5, a=4/5.

Same spatial endpoints r=6, Delta phi=1.8; zero winding, inward/outward
single-turn candidates. Optical reduction matches manuscript eqs randers,
randers-data; local PMP check uses project Hamiltonians as an independent check.
No global minimality outside this class is claimed.
"""
from pathlib import Path
import json
import mpmath as mp
mp.mp.dps=35
D=mp.mpf
a=D('.8');E=D('1.4');r0=D(6);target=D('1.8')

def coefficients(branch,r):
    f=1-2/r;de=r*r-2*r+a*a
    if branch=='tau': P=f*r*r/((E*E-f)*de);Q=de/(E*E-f);b=D(0)
    else: P=E*E*r*r/(f*(E*E-f)*de);Q=E*E*de/(f*f*(E*E-f));b=-2*a/(r*f)
    return P,Q,b

def orbit(branch,rm):
    Pm,Qm,bm=coefficients(branch,rm);j=bm+mp.sqrt(Qm)
    def delta(r):
        P,Q,b=coefficients(branch,r)
        return Q-(j-b)**2
    derivative=mp.diff(delta,rm)
    derivative2=mp.diff(delta,rm,2)
    derivative3=mp.diff(delta,rm,3)
    if derivative<=0:return None
    # A rational denominator may turn again; reject an inaccessible leg.
    if any(delta(rm+(r0-rm)*D(i)/20)<0 for i in range(1,21)):return None
    def kernels(u):
        w=u*u;r=rm+w;P,Q,b=coefficients(branch,r)
        lim=mp.power(10,-mp.mp.dps/3)
        d=delta(r)/w if w>lim else derivative+derivative2*w/2+derivative3*w*w/6
        if d<=0:raise ValueError('forbidden optical segment')
        kphi=2*(j-b)*mp.sqrt(P)/(mp.sqrt(Q)*mp.sqrt(d))
        kalpha=2*mp.sqrt(P*Q)/mp.sqrt(d)
        if branch=='tau': kt=E/(1-2/r)*kalpha-2*a/(r*(1-2/r))*kphi;ktau=kalpha
        else: kt=kalpha+b*kphi;ktau=(1-2/r)/E*kalpha
        return kphi,kt,ktau
    umax=mp.sqrt(r0-rm)
    vals=[2*mp.quad(lambda u:kernels(u)[i],[0,umax/2,umax]) for i in range(3)]
    return dict(branch=branch,rmin=rm,J=j,angle=vals[0],t=vals[1],tau=vals[2])

def solve(branch,nscan=90):
    # Audited sampled enumeration; do not promote it to completeness theorem.
    samples=[];roots=[]
    for i in range(nscan+1):
        rm=D('2.001')+(D('5.999')-D('2.001'))*D(i)/nscan
        row=orbit(branch,rm)
        if row:samples.append((rm,row['angle']-target))
        else:samples.append((rm,None))
    for (lo,fl),(hi,fh) in zip(samples[:-1],samples[1:]):
        if fl is None or fh is None or fl*fh>0:continue
        fun=lambda r:orbit(branch,r)['angle']-target
        rm=mp.findroot(fun,(lo,hi),solver='anderson',tol=D('1e-25'))
        row=orbit(branch,rm);assert abs(row['angle']-target)<D('1e-24')
        roots.append(row)
    assert roots
    roots.sort(key=lambda row:row['t'] if branch=='t' else row['tau'])
    return roots

def pmp_check(row):
    branch=row['branch'];j=row['J'];rm=row['rmin'];r=(rm+r0)/2
    P,Q,b=coefficients(branch,r);q=j-b;rad=Q-q*q
    pr=mp.sqrt(P*rad/Q)
    # Hamiltonians copied algebraically from colormap_spin_fixed_endpoints.build.
    f=1-2/r;de=r*r-2*r+a*a;PP=r*r+a*a+2*a*a/r;bb=2*a/r
    v2=1-f/E**2;Pb=PP+bb*bb/E**2;wind=bb*v2/Pb;R=mp.sqrt(v2*de/Pb)
    pshift=j-bb/E if branch=='tau' else j
    H=pshift*wind+R*mp.sqrt(de/r**2*pr**2+pshift**2/Pb)-(f/E if branch=='tau' else 1)
    assert abs(H)<D('1e-25'),(branch,H)
    # Different parametrisations have the same tangent slope.
    sqrt=mp.sqrt(de/r**2*pr**2+pshift**2/Pb)
    vr=R*(de/r**2)*pr/sqrt;vp=wind+R*pshift/(Pb*sqrt)
    optical_slope=q*mp.sqrt(P)/(mp.sqrt(Q)*mp.sqrt(rad))
    assert abs(vp/vr-optical_slope)<D('1e-25')
    # Reconstruct timelike velocity, charge, clocks and mass shell independently.
    rate=(f+bb*vp)/E;ut=1/rate;ur=vr/rate;uphi=vp/rate
    shell=-f*ut*ut-2*bb*ut*uphi+PP*uphi*uphi+r*r/de*ur*ur
    assert abs(shell+1)<D('1e-25')
    assert abs(f*ut+bb*uphi-E)<D('1e-25')
    return dict(H_residual=H,slope_residual=vp/vr-optical_slope,mass_shell_residual=shell+1)

def ode_check(row):
    # Reuse the existing project build verbatim, without its plotting/grid driver.
    import sympy as sp
    import numpy as np
    import hashlib
    from scipy.integrate import solve_ivp
    source=Path(__file__).resolve().parents[4]/'KerrScripts/colormap_spin_fixed_endpoints.py'
    text=source.read_text();part=text[text.index('r, pr, J, a ='):text.index('\nFN =')]
    namespace=dict(sp=sp,np=np,M=1.,E=float(E))
    exec(compile(part,str(source),'exec'),namespace)
    H,dHp,dHr,dHJ=namespace['build'](row['branch'])
    r=float(r0);rm=row['rmin'];j=float(row['J'])
    P,Q,b=coefficients(row['branch'],r0)
    p0=-float(mp.sqrt(P*(Q-(row['J']-b)**2)/Q))
    def rhs(t,y):
        rv,pv,phi,tau=y
        vr=dHp(float(a),rv,pv,j);vp=dHJ(float(a),rv,pv,j)
        return [vr,-dHr(float(a),rv,pv,j),vp,(1-2/rv+2*float(a)/rv*vp)/float(E)]
    def arrival(t,y):return y[0]-r
    arrival.direction=1;arrival.terminal=True
    def turning(t,y):return rhs(t,y)[0]
    turning.direction=1;turning.terminal=False
    outputs=[]
    for tol in [1e-9,1e-12]:
        sol=solve_ivp(rhs,[0,100],[r,p0,0,0],method='DOP853',rtol=tol,atol=tol/100,
             events=[arrival,turning],max_step=.15)
        assert sol.success and len(sol.t_events[0])==1 and len(sol.t_events[1])==1
        yy=sol.y_events[0][0]
        residuals=dict(t=sol.t_events[0][0]-float(row['t']),tau=yy[3]-float(row['tau']),
             angle=yy[2]-float(target),rmin=sol.y_events[1][0][0]-float(rm))
        assert max(abs(z) for z in residuals.values())<1e-8
        outputs.append(dict(rtol=tol,residuals=residuals))
    return dict(source='KerrScripts/colormap_spin_fixed_endpoints.py',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),comparisons=outputs)

if __name__=='__main__':
    rows={b:solve(b) for b in ['t','tau']}
    checks=[]
    for b in rows:
        base=rows[b][0]
        repeated_scan=solve(b,180)
        assert len(repeated_scan)==len(rows[b])
        assert abs(base['rmin']-repeated_scan[0]['rmin'])<D('1e-20')
        with mp.workdps(55):
            rm=mp.findroot(lambda r:orbit(b,r)['angle']-target,(base['rmin']-D('.001'),base['rmin']+D('.001')),tol=D('1e-40'))
            repeat=orbit(b,rm)
        for k in ['rmin','t','tau']:
            assert abs(base[k]-repeat[k])<D('1e-20'),(b,k,base[k]-repeat[k])
        checks.append(dict(branch=b,precision_repeat_max_error=max(abs(base[k]-repeat[k]) for k in ['rmin','t','tau']),pmp=pmp_check(base),ode=ode_check(base)))
    trow=rows['t'][0];taurow=rows['tau'][0]
    assert trow['t']<taurow['t'] and taurow['tau']<trow['tau']
    scales=json.loads(Path(__file__).with_name('observational_scales.json').read_text())
    bh=next(z for z in scales['scales'] if z['name']=='Gaia BH3')
    tg=D(bh['tg_s']);rg=D(bh['rg_km'])
    for b in rows:
      for row in rows[b]:
        row['t_s_BH3']=row['t']*tg;row['tau_s_BH3']=row['tau']*tg
        row['rmin_BL_coordinate_km_BH3']=row['rmin']*rg
    out=dict(inputs=dict(M=1,A=1,Ehat=E,a=a,rA=r0,rB=r0,delta_phi=target),
       search='91 samples of turning radius, 2.001<=rmin<=5.999; no completeness/global claim',
       precisions=[35,55],scans=[91,181],candidates=rows,checks=checks,
       depth_difference=trow['rmin']-taurow['rmin'])
    Path(__file__).with_name('rotating_bvp.json').write_text(json.dumps(out,indent=2,default=lambda v:mp.nstr(v,35)))
    print('PASS endpoint residual, precision repeat, PMP shell and tangent, charge/mass-shell, clock ordering')
    for b,rr in rows.items():
      for row in rr:print(b,{k:mp.nstr(row[k],15) for k in ['rmin','J','t','tau','t_s_BH3','tau_s_BH3']})
    print('depth difference',mp.nstr(out['depth_difference'],18))
