"""Symbolic large-r expansion of geometric HH on BL eta,r surfaces."""
from pathlib import Path
import json
import sympy as s
r,u,a,M,z=s.symbols('r u a M z',real=True)
S=r*r+a*a*u*u;d=r*r-2*M*r+a*a;B=(r*r+a*a)**2-a*a*d*(1-u*u)
wr=s.diff(2*M*a*r/B,r);wu=s.diff(2*M*a*r/B,u)
ex=d/(2*S)*(2*r/S)*(s.diff(B,r)/B-2*r/S)
ku=B**2*(1-u*u)**2*wu**2/(2*d*S**3)
tw=B**2*(1-u*u)*wr**2/(2*S**3)
bar=s.sqrt(s.expand(B.subs(r,1/z)*z**4))
I=s.series(bar*(ex+ku+tw).subs(r,1/z)/z**2,z,0,5).removeO().expand()
integ=s.integrate(I,(u,-1,1))
area=s.integrate(s.series(bar,z,0,5).removeO(),(u,-1,1))/2
area_radius_factor=s.series(s.sqrt(area),z,0,5).removeO()
mass=s.series(area_radius_factor/(2*z)*(1-integ/4),z,0,4).removeO().expand()
assert s.simplify(mass-(M-4*M*M*a*a*z**3))==0
correction=s.series(s.expand(B.subs(r,1/z)*z**4)**s.Rational(3,2)/
                    ((1+a*a*u*u*z*z)*(1-2*M*z+a*a*z*z)),z,0,4).removeO().expand()
C=s.series(area_radius_factor*s.integrate(correction,(u,-1,1))/4,z,0,4).removeO().expand()
assert s.simplify(C.subs(z,0)-s.Rational(1,2))==0
W=d*S/B
flux_norm=(d/S*s.diff(W,r)**2/W+(1-u*u)/S*s.diff(W,u)**2/W)/(4*M*M/r**4)
shear_norm=B*(1-u*u)/S**2*(d*wr**2+(1-u*u)*wu**2)/(4*M*M/r**4)
qn=s.series(flux_norm.subs(r,1/z),z,0,3).removeO().expand()
sn=s.series(shear_norm.subs(r,1/z),z,0,3).removeO().expand()
assert s.simplify(qn-(1+a*a*(1-7*u*u)*z*z))==0
assert s.simplify(sn-9*a*a*(1-u*u)*z*z)==0
out=dict(z='1/r',HH_Kerr_through_r_minus3=str(mass),
         C_over_r_cubed_through_r_minus3=str(C),
         q_squared_over_spherical_flux_squared=str(qn),
         s_squared_over_spherical_flux_squared=str(sn),
         assumptions='BL eta,r surfaces; fixed M,a; formal r -> infinity')
Path(__file__).with_name('rotating_asymptotics.json').write_text(json.dumps(out,indent=2))
print(out)
