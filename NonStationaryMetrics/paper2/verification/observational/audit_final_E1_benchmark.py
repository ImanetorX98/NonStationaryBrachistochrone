"""Ehat=1 final audit, frozen Kerr, common spatial endpoints, no source edits.
Reuses the independently checked BVP quadrature and project's Hamiltonian build.
"""
from pathlib import Path
import json
import mpmath as mp
import audit_rotating_bvp as bvp
mp.mp.dps=45
D=mp.mpf
bvp.E=D(1);bvp.a=D('.8');bvp.r0=D(6);bvp.target=D('1.8')

def zamo_gamma(row,r):
    f=1-2/r;de=r*r-2*r+bvp.a*bvp.a
    B=(r*r+bvp.a*bvp.a)**2-bvp.a*bvp.a*de
    N=mp.sqrt(de*r*r/B)
    P,Q,b=bvp.coefficients(row['branch'],r)
    if row['branch']=='tau': uphi=row['J']/Q
    else: uphi=bvp.E*(row['J']-b)/(f*Q)
    ut=(bvp.E-2*bvp.a/r*uphi)/f
    assert ut>0
    return N*ut

def speeds(row):
    rm=row['rmin'];grid=[rm+(6-rm)*D(i)/100 for i in range(101)]
    candidates=[rm,D(6)]
    dg=lambda r:mp.diff(lambda z:zamo_gamma(row,z),r)
    for l,h in zip(grid[:-1],grid[1:]):
        if dg(l)*dg(h)<0:candidates.append(mp.findroot(dg,(l,h)))
    values=[(mp.sqrt(1-zamo_gamma(row,r)**-2),r,zamo_gamma(row,r)) for r in candidates]
    lo=min(values);hi=max(values)
    return dict(minimum=lo[0],minimum_r=lo[1],maximum=hi[0],maximum_r=hi[1],
       static_speed_at_endpoints=mp.sqrt(D(1)/3),
       static_maximum=mp.sqrt(2/rm),max_gamma_ZAMO=hi[2])

if __name__=='__main__':
    rows={branch:bvp.solve(branch,90) for branch in ['t','tau']}
    out=dict(inputs=dict(A=1,M=1,Ehat=1,a='.8',rA=6,rB=6,delta_phi='1.8'),
       candidates=rows,checks=[],search='zero winding, single turn, rmin in [2.001,5.999], no completeness theorem')
    scales=json.loads(Path(__file__).with_name('observational_scales.json').read_text())
    bh=next(z for z in scales['scales'] if z['name']=='Gaia BH3')
    for branch,rr in rows.items():
        refined=bvp.solve(branch,180)
        assert len(refined)==len(rr)
        for row in rr:
            row['ZAMO_speed_over_c']=speeds(row)
            row['t_ms_BH3']=row['t']*D(bh['tg_s'])*1000
            row['tau_ms_BH3']=row['tau']*D(bh['tg_s'])*1000
            row['rmin_BL_km_BH3']=row['rmin']*D(bh['rg_km'])
            row['freeze_margin_min_vbar_squared']=D(1)/3
            row['distance_from_stationary_limit_M']=row['rmin']-2
        base=rr[0]
        with mp.workdps(65):
            rm=mp.findroot(lambda r:bvp.orbit(branch,r)['angle']-bvp.target,
                (base['rmin']-D('.001'),base['rmin']+D('.001')),tol=D('1e-50'))
            repeat=bvp.orbit(branch,rm)
        error=max(abs(base[k]-repeat[k]) for k in ['rmin','t','tau'])
        assert error<D('1e-22')
        out['checks'].append(dict(branch=branch,precision_repeat_error=error,
             pmp=bvp.pmp_check(base),ode=bvp.ode_check(base)))
    assert rows['t'][0]['t']<rows['tau'][0]['t']
    assert rows['tau'][0]['tau']<rows['t'][0]['tau']
    # E=1 static curvature bracket; no division by E^2-1 allowed.
    out['E1_static_Kt']='M^2*(8*M-5*r)/r^5 < 0 for r>2M'
    # Supplementary cruise arithmetic is conditional, not an adopted Voyager trajectory.
    c=D(scales['constants']['c']);pc=D(scales['constants']['pc']);year=D('365.25')*86400
    v=D(16900);gamma=1/mp.sqrt(1-(v/c)**2)
    out['conditional_16p9_km_s']=dict(speed_source_status='illustrative cruise speed (Voyager class); mission speed and epoch not independently sourced here',
      v_over_c=v/c,gamma_minus_one=gamma-1,BH1_477pc_cruise_years=D(477)*pc/v/year,
      BH3_590p6pc_cruise_years=D('590.6')*pc/v/year)
    # Exact symbolic specialization of the curvature formula printed in Paper II.
    import sympy as sp
    e2,rs,ms=sp.symbols('e2 r M',positive=True)
    bracket=16*ms**3+2*(9*e2-14)*ms**2*rs+(3*e2**2-19*e2+16)*ms*rs**2-(2*e2-3)*(e2-1)*rs**3
    curvature=ms*bracket/(e2*rs**5*(2*ms+(e2-1)*rs))
    assert sp.simplify(curvature.subs(e2,1)-ms**2*(8*ms-5*rs)/rs**5)==0
    exact_e2=gamma**2
    br=lambda r:16+2*(9*exact_e2-14)*r+(3*exact_e2**2-19*exact_e2+16)*r**2-(2*exact_e2-3)*(exact_e2-1)*r**3
    out['conditional_cruise_Kt_sign_change_r_over_M']=mp.findroot(br,(D(50000),D(60000)))
    Path(__file__).with_name('audit_final_E1_benchmark.json').write_text(json.dumps(out,indent=2,default=lambda x:mp.nstr(x,40)))
    print('PASS E1 shooting, 91/181 grid, 45/65 digits, project ODE, constraints and clock ordering')
    for branch,rr in rows.items():
        for row in rr:print(branch,{k:(row[k] if isinstance(row[k],dict) else mp.nstr(row[k],14)) for k in ['rmin','J','t','tau','t_ms_BH3','tau_ms_BH3','ZAMO_speed_over_c']})
    print(out['conditional_16p9_km_s'])
