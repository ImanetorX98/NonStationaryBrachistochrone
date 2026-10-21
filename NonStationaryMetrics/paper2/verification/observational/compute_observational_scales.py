"""R3 observational audit: conditional scales, not a fit to astronomical data.

Run with Python and mpmath. Output is generated beside this script.
All lengths use the nominal solar GM (Prsa 2016, PDF p5, Table 1).
H0 and Omega_m: Planck VI, published 2020, PDF p1, flat base LCDM.
AU, c and Julian year are unit definitions; parsec = 648000 AU / pi.
"""
from pathlib import Path
import json
import mpmath as mp

mp.mp.dps = 60
D = mp.mpf
c = D(299792458)
gm_sun = D('1.3271244e20')
au = D(149597870700)
pc = D(648000)*au/mp.pi
year = D('365.25')*86400
h0 = D('67.4')*1000/(D(10)**6*pc)
om = D('.315')
ol = 1-om
rg_sun = gm_sun/c**2
tg_sun = gm_sun/c**3
objects = [
 ('Gaia BH1','9.62','ElBadry BH1 PDF p1 / Table 1 p6'),
 ('Gaia BH2','8.94','ElBadry BH2 Table 2 PDF p6'),
 ('Gaia BH3','32.70','Gaia 2024 eq(5) PDF p7'),
 ('OGLE-2011-BLG-0462','7.1','Sahu 2022 section 10.1'),
 ('Sgr A*','4297000','GRAVITY 2022 Table 4 PDF p10, statistical fit'),
 ('M87*','6500000000','EHT 2019 VI abstract PDF p1'),
 ('Cygnus X-1','21.2','Miller-Jones 2021 Table 1 PDF p17'),
]

def integrals():
    f = lambda x: 1-2/x
    dt = mp.quad(lambda x: mp.sqrt(2)/(f(x)*mp.sqrt(2-f(x))),[6,8,10])
    tau = mp.quad(lambda x: 1/mp.sqrt(2-f(x)),[6,8,10])
    return dt,tau

dt,tau = integrals()
with mp.workdps(90):
    dt_check,tau_check=integrals()
    assert abs(dt-dt_check)<D('1e-55')
    assert abs(tau-tau_check)<D('1e-55')
assert dt>tau>0
rows=[]
for name,mass,source in objects:
    m=D(mass); rg=m*rg_sun; tg=m*tg_sun
    rows.append(dict(name=name,mass_msun=mass,source=source,
       rg_km=rg/1000,tg_s=tg,epsilon=h0*tg,r0_6M_km=6*rg/1000,
       delta_t_s=dt*tg,delta_tau_s=tau*tg,
       receiver_elapsed_s=mp.sqrt(D('.8'))*dt*tg))
# Conditional formation-age exercise, neglecting radiation in flat matter+Lambda.
t0=2*mp.asinh(mp.sqrt(ol/om))/(3*h0*mp.sqrt(ol))
ages=[]
for age in [5,12,13]:
    tb=t0-D(age)*year*D(10)**9
    assert tb>0
    ab=(mp.sqrt(om/ol)*mp.sinh(D('1.5')*h0*mp.sqrt(ol)*tb))**(D(2)/3)
    ages.append(dict(illustrative_age_Gyr=age,A_birth=ab,
       BH3_birth_k1=D('32.7')*ab,BH3_birth_k3=D('32.7')*ab**3,
       BH2_birth_k1=D('8.94')*ab))
assert h0*tg_sun<D('1e-22')
out=dict(constants=dict(GM_sun_nominal=gm_sun,c=c,AU=au,pc=pc,
    H0_s_inverse=h0,rg_sun_km=rg_sun/1000,tg_sun_s=tg_sun),
    benchmark=dict(A0=1,spin=0,Ehat='sqrt(2)',J=0,r0_over_M=6,r1_over_M=10,
       branch='outgoing radial; illustrative, not fitted',delta_t_over_M=dt,
       delta_tau_over_M=tau,quadrature_checked_dps=[60,90]),
    scales=rows,age_exercise=ages,t0_Gyr=t0/year/D(10)**9,
    A_fractional_change_over_16_Julian_years=h0*16*year,
    caveat='epsilon assumes A0=1 and A is the cosmological scale factor; no observation fit')
def encode(o):
    if isinstance(o,mp.mpf): return mp.nstr(o,45)
    raise TypeError(type(o))
Path(__file__).with_name('observational_scales.json').write_text(
    json.dumps(out,indent=2,default=encode)+'\n')
print('Quadrature and consistency checks PASS')
print('dt/M',mp.nstr(dt,16),'dtau/M',mp.nstr(tau,16))
for r in rows:
    print(r['name'],*[mp.nstr(r[k],7) for k in ['rg_km','tg_s','epsilon','r0_6M_km','delta_t_s','delta_tau_s']])
print('t0/Gyr',mp.nstr(t0/year/D(10)**9,9))
for a in ages: print({k:mp.nstr(v,7) if isinstance(v,mp.mpf) else v for k,v in a.items()})
