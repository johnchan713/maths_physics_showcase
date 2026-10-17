"""Exact finite-harmonic source grammar controls; manufactured, not actual fields."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def algebra():
    p=HERE.parent/'first_sourced_foundation/sourced_checks.py'
    spec=importlib.util.spec_from_file_location('finite_foundation_complex_algebra',p)
    m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m


@lru_cache(None)
def fixture():
    s=algebra();e=s.exact();R,Z,T,eps=e.R,e.Z,e.T,e.EPS;zero=e.Poly()
    bump=(R-1)**3*(2-R)**3;data={};potentials={};k=3
    for m in (1,2,4):
        f=R**2*bump*(1+Z/(7+m)+T/(11+m))
        c=[(f/(m+2),f*(1+Z/13)/(m+3)),(R*f/(m+5),-f/(m+7)),(-R*f/(m+11),f*(1+T/17)/(m+13))]
        diff=lambda x,a:(x[0].diff(a),x[1].diff(a))
        radial=lambda x:(e.radial(x[0]),e.radial(x[1]))
        shift=lambda x:(x[0].shift(-1),x[1].shift(-1))
        cross=[s.add(shift(c[2]),s.neg(c[1])),c[0],s.neg(shift(c[0]))]
        curl=[s.scale(diff(c[1],1),-eps),s.add(s.scale(diff(c[0],1),eps),s.neg(diff(c[2],0))),radial(c[1])]
        a=[s.add(x,s.scale(s.imag(y),k*m)) for x,y in zip(curl,cross)]
        data[m]=a;data[-m]=[s.conj(x) for x in a];potentials[m]=c
    psi=R*bump*(1+Z/19+T/23)
    mean=[(-eps*psi.diff(1),zero),((R**4*bump).diff(0).shift(-2)*(1+Z/29),zero),(e.radial(psi),zero)]
    return data,mean,potentials


def zero_vector():
    z=algebra().exact().Poly();return [(z,z) for _ in range(3)]


def plus_vectors(a,b):return [algebra().add(x,y) for x,y in zip(a,b)]


def transport(a,b,n):
    s=algebra();e=s.exact();eps=e.EPS
    dot=s.add((a[1][0].shift(-1),a[1][1].shift(-1)),a[2])
    result=[]
    for bi in b:
        dr=(bi[0].diff(0),bi[1].diff(0));dz=(bi[0].diff(1),bi[1].diff(1))
        result.append(s.add(s.add(s.mul(a[0],dr),s.scale(s.mul(a[2],dz),eps)),s.scale(s.imag(s.mul(dot,bi)),3*n)))
    connection=s.mul(a[1],b[1]);result[0]=s.add(result[0],s.neg((connection[0].shift(-1),connection[1].shift(-1))))
    connection=s.mul(a[1],b[0]);result[1]=s.add(result[1],(connection[0].shift(-1),connection[1].shift(-1)))
    return result


def convolution(a,b):
    out={}
    for m,u in a.items():
        for n,v in b.items():out[m+n]=plus_vectors(out.get(m+n,zero_vector()),transport(u,v,n))
    return {m:v for m,v in out.items() if any(bool(r) or bool(i) for r,i in v)}


@lru_cache(None)
def harmonic_review():
    s=algebra();e=s.exact();eps=e.EPS;waves,mean,_=fixture()
    divs=[];longitudinal=[];nonzero_normal=False;zero_output=False
    for m,a in waves.items():
        div=s.add((e.radial(a[0][0]),e.radial(a[0][1])),s.scale((a[2][0].diff(1),a[2][1].diff(1)),eps))
        normal=s.add((a[1][0].shift(-1),a[1][1].shift(-1)),a[2])
        divs.extend(x==0 for x in s.add(div,s.scale(s.imag(normal),3*m)))
        nonzero_normal |= any(bool(x) for x in normal)
        for n,b in waves.items():
            for bi in b:
                actual=s.scale(s.imag(s.mul(normal,bi)),3*n)
                reduced=s.scale(s.mul(div,bi),-F(n,m))
                longitudinal.extend(x==y for x,y in zip(actual,reduced))
                if m+n==0:zero_output |= actual==reduced
    source=convolution(waves,waves);wm=convolution(waves,{0:mean});mw=convolution({0:mean},waves)
    reality=all(source[-m]==[s.conj(x) for x in v] for m,v in source.items())
    mixed=plus_vectors(mw[1],wm[1]);square=source[2]
    # Exact quadratic increment identity for arbitrary current mean and all modes.
    delta={1:[s.scale(x,F(1,7)) for x in waves[1]],-1:[s.scale(x,F(1,7)) for x in waves[-1]]}
    old=dict(waves);old[0]=mean;updated={m:plus_vectors(old.get(m,zero_vector()),delta.get(m,zero_vector())) for m in set(old)|set(delta)}
    before=convolution(old,old);after=convolution(updated,updated)
    terms=(convolution(old,delta),convolution(delta,old),convolution(delta,delta));count=0;identity=True
    for m in set(before)|set(after)|set().union(*(set(t) for t in terms)):
        expected=zero_vector()
        for t in terms:expected=plus_vectors(expected,t.get(m,zero_vector()))
        difference=plus_vectors(after.get(m,zero_vector()),[s.neg(x) for x in before.get(m,zero_vector())])
        for a,b in zip(difference,expected):identity &= a==b;count+=2
    return dict(exact_complete_divergence_components=len(divs),exact_phase_reduction_components=len(longitudinal),
        exact_nonlinear_increment_components=count,nonzero_source_harmonics=sorted(source),
        checks=dict(all_complete_curl_harmonics_divergence_free=all(divs),
            complete_amplitudes_not_assumed_exactly_transverse=nonzero_normal,
            every_pair_phase_transport_reduces_without_half_frequency_loss=all(longitudinal),
            zero_output_harmonics_do_not_require_division_by_output=zero_output and 0 in source,
            conjugate_source_harmonics_preserve_reality=reality,
            highest_output_harmonic_eight_is_present=max(abs(m) for m in source)==8,
            mean_wave_transport_not_removed_by_wave_incompressibility=any(bool(x) or bool(y) for x,y in mixed),
            new_wave_square_nonzero=any(bool(x) or bool(y) for x,y in square),
            exact_full_nonlinear_increment_has_both_crosses_and_new_square=bool(identity)),
        graph_free_single_label_polynomial_diagnostic_not_actual_two_torus_pulse=True,
        finite_endpoint_regularity_not_actual_C_infinity_certificate=True)


@lru_cache(None)
def covariance_review():
    s=algebra();e=s.exact();waves,mean,_=fixture();count=0;matches=True;all_six=[]
    zero=e.Poly()
    for i in range(3):
        for j in range(3):
            value=s.sumc([s.mul(a[i],waves[-m][j]) for m,a in waves.items()],zero)
            independent=sum(2*(a[i][0]*a[j][0]+a[i][1]*a[j][1]) for m,a in waves.items() if m>0)
            matches &= value==(independent,zero);count+=2
            if i<=j:all_six.append(bool(independent))
    return dict(exact_zero_mode_covariance_components=count,
        checks=dict(independent_real_quadratures_match_full_Fourier_covariance=bool(matches),
            all_six_entries_including_axial_and_mixed_are_nonzero=all(all_six),
            mean_zero_mode_kept_separate_from_angular_wave_covariance=0 not in waves and len(mean)==3),
        full_pressure_recomputed_by_retained_current_mean_maps_not_this_fixture=True)


def reviews():return dict(finite_harmonic_complete_source_controls=harmonic_review(),
    higher_harmonic_complete_covariance=covariance_review())
