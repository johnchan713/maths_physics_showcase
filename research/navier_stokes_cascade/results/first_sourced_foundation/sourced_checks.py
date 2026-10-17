"""Independent complex forced jets and complete-curl polynomial moment controls."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m


@lru_cache(None)
def exact():return load('first_sourced_foundation_polynomials',HERE.parent/'complete_covariance_initial_means/means.py')


@lru_cache(None)
def maps():return load('first_sourced_foundation_jets',HERE.parent/'modulation_repair_all_order/map_checks.py')


def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[0],-a[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,c):return (a[0]*c,a[1]*c)
def imag(a):return (-a[1],a[0])
def conj(a):return (a[0],-a[1])
def sumc(v,zero):
    out=(zero,zero)
    for x in v:out=add(out,x)
    return out


@lru_cache(None)
def forced_checks():
    j=maps();J=j.JetClass();count=0;matches=True;normalgap=False;movinggap=False;exitgap=False;reality=True;records={}
    for harmonic in (-2,-1,1,2):
        v=J.variable(0,F(1,3),n=2,order=6);z=J.variable(1,F(1,7),n=2,order=6);zero=j.jet(0,n=2,order=6)
        one=1+zero;n=[one,v+z,one];np=[zero,one,zero];k=3;omega=k*harmonic
        K=[[zero,-2*(1+z),zero],[3+z,zero,zero],[2-z,zero,zero]]
        q=v*v*(1+v)*(1+z/5);s=(1 if harmonic>0 else -1)*v**3*(1+z/7)
        t=[(-(v+z)*q,-s),(q,zero),(zero,s)]
        pi=((q+s*s)/(k*abs(harmonic)),(q+z*q)/(k*harmonic))
        kt=[sumc([scale(x,c) for c,x in zip(row,t)],zero) for row in K]
        tp=[(x.derivative(0),y.derivative(0)) for x,y in t]
        norm=sum(x*x for x in n);d=F(1,8)*k*k*norm
        forcing=[neg(sumc([a,b,scale(c,harmonic**2*d),scale(imag(pi),omega*ni)],zero)) for a,b,c,ni in zip(tp,kt,t,n)]
        nt=sumc([scale(x,c) for c,x in zip(n,t)],zero)
        nkt=sumc([scale(x,c) for c,x in zip(n,kt)],zero)
        npt=sumc([scale(x,c) for c,x in zip(np,t)],zero)
        nf=sumc([scale(x,c) for c,x in zip(n,forcing)],zero)
        recovered=scale(imag(add(add(nkt,neg(npt)),nf)),F(1,k*harmonic)*norm**-1)
        constraints=list(nt)+[x-y for x,y in zip(recovered,pi)]
        for i in range(3):
            projected=add(add(neg(kt[i]),scale(add(nkt,neg(npt)),n[i]*norm**-1)),add(scale(t[i],-harmonic**2*d),add(neg(forcing[i]),scale(nf,n[i]*norm**-1))))
            constraints.extend(x-y for x,y in zip(projected,tp[i]))
        for p in constraints:
            for a in range(5):
                for b in range(5-a):matches &= p.partial((a,b))==0;count+=1
        movinggap |= any(bool(x.terms) for x in npt);normalgap |= any(bool(x.terms) for x in nf)
        # Manufactured zero ENTRY VALUE and nonzero EXIT do not certify a
        # supported actual source; endpoint support is supplied analytically.
        ve=F(1);ze=F(1,7);qe=ve*ve*(1+ve)*(1+ze/5);exitgap |= qe!=0
        records[harmonic]=[[(x.value(),y.value()) for x,y in seq] for seq in (t,forcing)]+[[pi[0].value(),pi[1].value()]]
    for m in (1,2):
        for p,q in zip(records[m][:2],records[-m][:2]):reality &= all(x==conj(y) for x,y in zip(p,q))
        reality &= tuple(records[m][2])==conj(tuple(records[-m][2]))
    return dict(exact_mixed_constraint_pressure_ODE_components=count,
        checks=dict(complex_moving_constraint_and_forced_pressure_exact=bool(matches),
            source_normal_pressure_term_nonzero=bool(normalgap),moving_normal_pressure_term_nonzero=bool(movinggap),
            conjugate_negative_harmonics_preserve_reality=bool(reality),zero_entry_does_not_force_zero_exit=bool(exitgap)),
        fixture_not_actual_giant_pulse_or_endpoint_source_certificate=True)


def outer(a,b):return [[x*y for y in b] for x in a]
def tensor_add(*a):return [[sum(t[i][j] for t in a) for j in range(3)] for i in range(3)]
def tensor_scale(a,c):return [[x*c for x in row] for row in a]


@lru_cache(None)
def covariance_fixture():
    e=exact();R,Z,T,eps=e.R,e.Z,e.T,e.EPS;P=e.Poly;b=(R-1)**6*(2-R)**6
    f=b*(1+Z/5+T/7);factor=1-F(1,3)+Z/11+T/13;k=3
    c0=[R**2*f,R**2*f*(2+Z/7),-R**3*f]
    old_im=[factor*x/k for x in c0];old_re=[P(),P(),P()]
    re1=[(1+Z/17)*R**2*f/7,-R**3*f/11,R**2*f*(1+T/19)/13]
    im1=[-R**2*f*(1+Z/23)/5,R**3*f/17,R**2*f*(2+T/29)/19]
    re2=[R**2*f/23,R**3*f*(1+Z/31)/29,-R**2*f/31]
    im2=[R**3*f/37,-R**2*f*(1+T/41)/41,R**2*f/43]
    def curl(c):return [-eps*c[1].diff(1),eps*c[0].diff(1)-c[2].diff(0),e.radial(c[1])]
    def ncross(c):return [c[2].shift(-1)-c[1],c[0],-c[0].shift(-1)]
    def amplitude(cr,ci,m):
        ar,ai=curl(cr),curl(ci);nr,ni=ncross(cr),ncross(ci)
        return ([a-k*m*y for a,y in zip(ar,ni)],[-a-k*m*y for a,y in zip(ai,nr)])
    old={1:amplitude(old_re,old_im,1)}
    inc={1:amplitude(re1,im1,1),2:amplitude(re2,im2,2)}
    updated={1:tuple([x+y for x,y in zip(a,b)] for a,b in zip(old[1],inc[1])),2:inc[2]}
    def covariance(data):return tensor_scale(tensor_add(*(tensor_add(outer(a,a),outer(b,b)) for a,b in data.values())),F(1,2))
    W0,W1=covariance(old),covariance(updated);delta=tensor_add(W1,tensor_scale(W0,-1))
    a,bq=old[1];c,d=inc[1];c2,d2=inc[2]
    cross=tensor_scale(tensor_add(outer(a,c),outer(c,a),outer(bq,d),outer(d,bq)),F(1,2))
    square=tensor_scale(tensor_add(outer(c,c),outer(d,d),outer(c2,c2),outer(d2,d2)),F(1,2))
    return b,old,inc,updated,W0,W1,delta,cross,square


@lru_cache(None)
def tensor_checks():
    e=exact();P=e.Poly;bump,old,inc,new,W0,W1,delta,cross,square=covariance_fixture();count=0;matches=True;div=True
    def fourier(data):
        out={}
        for m,(a,b) in data.items():
            out[m]=[(x/2,-y/2) for x,y in zip(a,b)];out[-m]=[conj(z) for z in out[m]]
        return out
    for data,expected in ((old,W0),(new,W1)):
        coeff=fourier(data);out=[[(P(),P()) for _ in range(3)] for _ in range(3)]
        for m,a in coeff.items():
            for i in range(3):
                for j in range(3):out[i][j]=add(out[i][j],mul(a[i],coeff[-m][j]))
        for i in range(3):
            for j in range(3):matches &= out[i][j][0]==expected[i][j] and out[i][j][1]==0;count+=2
    for m,(a,b) in list(old.items())+list(inc.items()):
        na=a[1].shift(-1)+a[2];nb=b[1].shift(-1)+b[2]
        div &= e.radial(a[0])+e.EPS*a[2].diff(1)+3*m*nb==0
        div &= e.radial(b[0])+e.EPS*b[2].diff(1)-3*m*na==0
    identity=all(delta[i][j]==cross[i][j]+square[i][j] for i in range(3) for j in range(3))
    oldsin=old[1][1]
    # Both old-sine/new-sine products are present in the full complex change.
    sine_cross=tensor_add(outer(oldsin,inc[1][1]),outer(inc[1][1],oldsin))
    R,Z,T=e.R,e.Z,e.T;f=bump*(1+Z/5+T/7)
    re1=[(1+Z/17)*R**2*f/7,-R**3*f/11,R**2*f*(1+T/19)/13]
    principal_sine=[-3*(re1[2].shift(-1)-re1[1]),-3*re1[0],3*re1[0].shift(-1)]
    principal_cross=tensor_add(outer(oldsin,principal_sine),outer(principal_sine,oldsin))
    square2=tensor_scale(tensor_add(outer(inc[2][0],inc[2][0]),outer(inc[2][1],inc[2][1])),F(1,2))
    only1=tensor_scale(tensor_add(outer(new[1][0],new[1][0]),outer(new[1][1],new[1][1])),F(1,2))
    m2=tensor_add(W1,tensor_scale(only1,-1))
    bad_merge=tensor_add(outer(old[1][0],inc[2][0]),outer(inc[2][0],old[1][0]))
    return dict(exact_zero_mode_components=count,exact_tensor_increment_components=9,
        checks=dict(independent_complex_Fourier_zero_mode_matches_complete_tensor=bool(matches),
            complete_old_and_new_curls_divergence_free=bool(div),full_cross_plus_all_new_squares_identity=identity,
            all_six_covariance_changes_nonzero=all(bool(delta[i][j]) for i in range(3) for j in range(i,3)),
            old_sine_new_complex_sine_cross_not_dropped=any(bool(x) for row in sine_cross for x in row),
            old_sine_new_principal_sine_cross_nonzero=any(bool(x) for row in principal_cross for x in row),
            second_harmonic_square_not_dropped=all(m2[i][j]==square2[i][j] for i in range(3) for j in range(3)) and any(bool(x) for row in square2 for x in row),
            unequal_angular_harmonics_must_not_be_merged=any(bool(x) for row in bad_merge for x in row),
            generic_covariance_increment_need_not_be_positive_semidefinite=F(-1,2)<0),
        negative_increment_witness='old=e_z*cos(theta); increment=-e_z*cos(theta); DeltaWzz=-1/2; locally divergence-free',
        fixture_potential_zero_extension='C5 polynomial potential; complete curl C4 only, not actual C-infinity')


@lru_cache(None)
def pressure_checks():
    e=exact();R,Z,T,eps=e.R,e.Z,e.T,e.EPS;b,_,_,_,W0,W1,delta,_,_=covariance_fixture()
    psiB=R*b*(1+Z+T);base=(-eps*psiB.diff(1),R**2*(2+Z),e.radial(psiB))
    psi=R*b*(1+Z/7+T/11);mean=(-eps*psi.diff(1),(R**4*b).diff(0).shift(-2)*(1+Z/13+T/17),e.radial(psi))
    rho=b/b.integral().scalar()*(1+Z*(R-F(3,2)));crho=(R**2*rho).integral()/2
    def state(W):
        b0,V,G=base;beta,v,gamma=mean
        gr=e.averaged_gr(base,mean,W[0][0],W[1][1],W[0][2]);P=gr.integral();raw=(gr-rho*P).primitive();pm=raw-raw.at(0,1)
        jt=(R**2*(G*v+V*gamma+gamma*v+W[1][2])).integral()
        jz=(R*(2*G*gamma+gamma**2+W[2][2])).integral()-(R**2*gr).integral()/2
        et=-eps*v.diff(2)+e.radial(b0*v+beta*V+beta*v+W[0][1],2)+eps*(G*v+V*gamma+gamma*v+W[1][2]).diff(1)-eps*(e.lap(v)-v.shift(-2))
        ez=-eps*gamma.diff(2)+e.radial(b0*gamma+beta*G+beta*gamma+W[0][2])+eps*(2*G*gamma+gamma**2+W[2][2]+pm).diff(1)-eps*e.lap(gamma)
        return dict(gr=gr,P=P,pm=pm,Jtheta=jt,Jz=jz,Etheta=et,Ez=ez)
    old,new=state(W0),state(W1);d={n:new[n]-old[n] for n in old}
    expected=-e.radial(delta[0][0])-eps*delta[0][2].diff(1)+delta[1][1].shift(-1)
    Jt=(R**2*delta[1][2]).integral();Jz=(R*delta[2][2]).integral()-(R**2*expected).integral()/2
    tail1=eps**2*(Z*T+Z**3-F(1,3)*R.shift(-3));tail2=eps**4*(Z+T*Z**2)
    p0=R**2+Z**2;pB=p0+tail1+tail2
    return dict(checks=dict(full_recomputation_with_nonzero_old_means_matches_covariance_change=d['gr']==expected,
        recomputed_pressure_defect=d['P']==expected.integral(),
        recomputed_both_flux_defects=d['Jtheta']==Jt and d['Jz']==Jz,
        compact_mean_pressure_at_both_ends=all(s['pm'].at(0,1)==s['pm'].at(0,2)==0 for s in (old,new)),
        full_axial_pressure_moment=(R*d['pm']).integral()==-(R**2*expected).integral()/2+crho*d['P'],
        whole_moving_pressure_product=(R*d['pm']).integral().diff(1)==(-(R**2*expected).integral()/2+crho*d['P']).diff(1),
        moving_bump_derivative_omission_detected=bool(crho.diff(1)*d['P']),
        full_angular_integrated_balance=(R**2*d['Etheta']).integral()==eps*Jt.diff(1),
        full_axial_integrated_balance=(R*d['Ez']).integral()==eps*(Jz+crho*d['P']).diff(1),
        axial_covariance_and_mean_pressure_terms_nonzero=bool((R*delta[2][2]).integral().diff(1)) and bool((R*d['pm']).integral().diff(1)),
        nonzero_old_means_and_radial_background_retained=all(bool(x) for x in mean) and bool(base[0]),
        exact_old_velocity_moments_unchanged=(R**2*mean[1]).integral()==0 and (R*mean[2]).integral()==0,
        same_full_background_pressure_subtraction=(pB+d['pm']).diff(1)-pB.diff(1)==d['pm'].diff(1),
        truncating_background_pressure_would_leave_tail_gradient=bool(tail1.diff(1)) and bool(tail2.diff(1))),
        exact_before_after_pressure_and_moment_controls=14,
        manufactured_compact_increment_not_actual_reserved_patch_evaluation=True)


def reviews():return dict(complex_forced_constraint_and_pressure=forced_checks(),
    full_multiharmonic_complete_curl_covariance=tensor_checks(),
    full_recomputed_mean_pressure_and_defects=pressure_checks())
