"""Independent geometric audit, G=c=1, BL surfaces eta,r constant.

HH definition and conformal law: Prain et al. arXiv1501.02977 eqs1,24,35.
No matter-mass interpretation is assumed. All generated output stays here.
"""
from pathlib import Path
import json
import mpmath as mp
mp.mp.dps=50
D=mp.mpf

def data(r,u,a,M=1):
    S=r*r+a*a*u*u; de=r*r-2*M*r+a*a
    B=(r*r+a*a)**2-a*a*de*(1-u*u)
    N=mp.sqrt(de*S/B); om=2*M*a*r/B
    return S,de,B,N,om

def hh(r,a,M=1):
    """HH from contracted Gauss plus anholonomicity, chosen BL foliation."""
    def integrand(u):
        S,de,B,N,om=data(r,u,a,M)
        Br=4*r*(r*r+a*a)-a*a*(2*r-2*M)*(1-u*u)
        omr=2*M*a*(B-r*Br)/B**2
        Bu=2*a*a*de*u
        omu=-2*M*a*r*Bu/B**2
        extr=de/(2*S)*(2*r/S)*(Br/B-2*r/S)
        shear_u=B*B*(1-u*u)**2*omu**2/(2*de*S**3)
        twist=B*B*(1-u*u)*omr**2/(2*S**3)
        return mp.sqrt(B)*(extr+shear_u+twist)
    area=2*mp.pi*mp.quad(lambda u:mp.sqrt(data(r,u,a,M)[2]),[-1,0,1])
    RA=mp.sqrt(area/(4*mp.pi))
    mass=RA/2*(1-mp.quad(integrand,[-1,0,1])/4)
    C=RA/4*mp.quad(lambda u:data(r,u,a,M)[2]**D('1.5')/
                   (data(r,u,a,M)[0]*data(r,u,a,M)[1]),[-1,0,1])
    return dict(r=r,a=a,M_HH_Kerr=mass,area=area,area_radius=RA,
                conformal_H_squared_coefficient=C)

def komar(r,a,M=1):
    # Orientation dt ^ dr ^ dtheta ^ dphi positive: star dxi density is negative.
    def density(u):
        S,de,B,N,om=data(r,u,a,M)
        q=B*(1-u*u)/S
        Br=4*r*(r*r+a*a)-a*a*(2*r-2*M)*(1-u*u)
        omr=2*M*a*(B-r*Br)/B**2
        # sqrt(-g) F^{tr} dtheta = q B om_r / S du, bounds reversed.
        return q*B*omr/S
    return -mp.quad(density,[-1,0,1])/8 # flip surface orientation for physical +Ma

def independent_R_projection(r,u,a,M=1):
    # Direct stationary R_theta phi theta phi from metric derivatives.
    S,de,B,N,om=data(r,u,a,M); sth=mp.sqrt(1-u*u)
    th=mp.acos(u)
    def metric_tphi(v):
        s=r*r+a*a*mp.cos(v)**2
        return mp.matrix([[-(1-2*M*r/s),-2*M*a*r*mp.sin(v)**2/s],
          [-2*M*a*r*mp.sin(v)**2/s,((r*r+a*a)**2-a*a*de*mp.sin(v)**2)*mp.sin(v)**2/s]])
    q=lambda rv,tv: data(rv,mp.cos(tv),a,M)[2]*mp.sin(tv)**2/data(rv,mp.cos(tv),a,M)[0]
    qt=mp.diff(lambda z:q(r,z),th); qtt=mp.diff(lambda z:q(r,z),th,2)
    qr=mp.diff(lambda z:q(z,th),r)
    St=-2*a*a*u*sth; Sr=2*r
    derivatives=mp.matrix([mp.diff(lambda z:metric_tphi(z)[0,1],th),qt])
    contraction=(derivatives.T*metric_tphi(th)**-1*derivatives)[0]
    R=-qtt/2-Sr*qr*de/(4*S)+St*qt/(4*S)+contraction/4
    return 2*R/(B*sth**2)

def gauss_projection(r,u,a,M=1):
    S,de,B,N,om=data(r,u,a,M);th=mp.acos(u)
    E=lambda t:data(r,mp.cos(t),a,M)[0]
    F=lambda t:data(r,mp.cos(t),a,M)[2]*mp.sin(t)**2/E(t)
    R2=-mp.diff(lambda t:mp.diff(F,t)/mp.sqrt(E(t)*F(t)),th)/mp.sqrt(E(th)*F(th))
    Br=mp.diff(lambda v:data(v,u,a,M)[2],r)
    omu=mp.diff(lambda v:data(r,v,a,M)[4],u)
    extr=de/(2*S)*(2*r/S)*(Br/B-2*r/S)
    ku=B*B*(1-u*u)**2*omu**2/(2*de*S**3)
    return R2-extr-ku

if __name__=='__main__':
    checks=[];rows=[]
    for r in [3,6,10,100,1000]:
      for a in [D(0),D('.5'),D('.9')]:
        row=hh(D(r),a);rows.append(row)
        assert abs(komar(D(r),a)-a)<D('1e-45')
        if not a:assert abs(row['M_HH_Kerr']-1)<D('1e-45')
    for r,u,a in [(D(6),D('.3'),D('.8')),(D(3),D('.6'),D('.5'))]:
        err=abs(independent_R_projection(r,u,a)-gauss_projection(r,u,a))
        assert err<D('1e-45');checks.append(dict(test='direct Riemann vs Gauss',error=err))
    for r in [3,6,10]:
        flat=hh(D(r),D('.8'),0)
        assert abs(flat['M_HH_Kerr'])<D('1e-45')
    spherical=hh(D(6),D(0))
    assert abs(spherical['conformal_H_squared_coefficient']-D(6)**3/(2*(1-D(2)/6)))<D('1e-45')
    with mp.workdps(75):
        repeat=hh(D(6),D('.8'))
    base=hh(D(6),D('.8'));assert abs(base['M_HH_Kerr']-repeat['M_HH_Kerr'])<D('1e-45')
    rows.append(base)
    scales=json.loads(Path(__file__).with_name('observational_scales.json').read_text())
    eps=D(next(x for x in scales['scales'] if x['name']=='Gaia BH3')['epsilon'])
    base['relative_HH_drift_correction_BH3_today']=eps**2*base['conformal_H_squared_coefficient']/base['M_HH_Kerr']
    out=dict(definition='geometric HH; BL eta,r constant surfaces; positive angular Komar orientation',
      precision=50,repeat_precision=75,checks=checks,rows=rows)
    Path(__file__).with_name('rotating_geometry.json').write_text(json.dumps(out,indent=2,default=lambda x:mp.nstr(x,40)))
    print('PASS Komar, spherical HH, flat oblate HH, direct Riemann vs Gauss, precision repeat')
    for row in rows:print('r,a,HH,C',*[mp.nstr(row[k],14) for k in ['r','a','M_HH_Kerr','conformal_H_squared_coefficient']])
