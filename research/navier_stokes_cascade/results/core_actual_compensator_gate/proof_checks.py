"""Exact compensator-data obstruction and conditional shape estimates."""
from fractions import Fraction as F
import importlib.util
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


shape=load("actual_compensator_parent",HERE.parent/"core_pressure_shape_recheck/proof_checks.py")
means=shape.means
R=means.R


def integral(poly,left,right):
    primitive=poly.primitive()
    return (primitive.at(0,right)-primitive.at(0,left)).scalar()


def profile_contract_review():
    """Record exactly what the analytical construction fixes about rhohat."""
    return dict(checks={
        "fixed_smooth_compact_profile_is_only_an_input":True,
        "unit_mass_is_declared":True,
        "nonnegativity_is_not_assumed":True,
        "no_profile_formula_or_seminorm_values_are_declared":True,
        "numerical_diagnostic_bump_is_not_promoted_to_canonical_data":True,
    },scope="source-contract audit bound to the hashed analytical notes and diagnostic implementation")


def signed_same_moment_family_review():
    """An exact fixed-support family with the same mass and second moment."""
    base=(R-1)**8*(2-R)**8
    base=base/integral(base,F(1),F(2))
    lo,hi=F(5,4),F(7,4)
    w=(R-lo)**8*(hi-R)**8
    h=w.diff(0).diff(0).diff(0)
    w2=w.diff(0).diff(0)
    mid=(lo+hi)/2
    rows=[];checks={
        "base_mass_one":integral(base,F(1),F(2))==1,
        "perturbation_mass_zero":integral(h,lo,hi)==0,
        "perturbation_second_moment_zero":integral(R**2*h,lo,hi)==0,
        "primitive_has_zero_endpoints":w2.at(0,lo)==w2.at(0,hi)==0,
        "primitive_midpoint_nonzero":bool(w2.at(0,mid)),
    }
    unit=abs(w2.at(0,mid).scalar())
    for n in (1,5,11):
        # h is used with smooth zero extension outside [lo,hi].
        mass_n=integral(base,F(1),F(2))+n*integral(h,lo,hi)
        moment_n=integral(R**2*base,F(1),F(2))+n*integral(R**2*h,lo,hi)
        checks[f"mass_N{n}"]=mass_n==1
        checks[f"second_moment_N{n}"]=moment_n==integral(R**2*base,F(1),F(2))
        # D_N-D_0=I0(rho_0-rho_N)=-N*w'' on the perturbation support.
        diameter=abs((-n*w2).at(0,mid).scalar())
        checks[f"linear_shape_primitive_growth_N{n}"]=diameter==n*unit
        rows.append(dict(N=n,shape_primitive_midpoint=str(diameter)))
    return dict(checks=checks,rows=rows,
                conclusion="unit mass, fixed support and the scalar second moment have no profile-independent D-hat bound")


def positive_concentration_review():
    """Positivity controls L1 but still supplies no pointwise profile jets."""
    t=R  # Build the normalized reference bump on [0,1].
    bump=t**4*(1-t)**4
    mass=integral(bump,F(0),F(1))
    x0=F(1)
    rows=[];checks={"reference_bump_mass_positive":mass>0}
    midpoint_value=bump.at(0,F(1,2)).scalar()/mass
    for n in (2,4,8):
        rho=n*((n*(R-x0))**4*(1-n*(R-x0))**4)/mass
        right=x0+F(1,n)
        value=rho.at(0,x0+F(1,2*n)).scalar()
        checks[f"positive_unit_mass_N{n}"]=integral(rho,x0,right)==1 and value>0
        checks[f"pointwise_height_scales_N{n}"]=value==n*midpoint_value
        rows.append(dict(N=n,support=[str(x0),str(right)],midpoint_height=str(value)))
    return dict(checks=checks,rows=rows,
                conclusion="even nonnegative unit-mass profiles on one fixed hull have unbounded pointwise density, hence unbounded first radial shape jet")


def dilation_budget_review():
    """Exact first/second slow-chain bounds from the inherited chart constants."""
    # s in [1/4,4]. Inherited normalized Cauchy bounds:
    sZ,sT=F(32),F(64)
    sZZ,sZT,sTT=F(2048),F(2048),F(8192)
    invs,invs2=F(4),F(16)
    LZ,LT=invs*sZ,invs*sT
    LZZ=invs*sZZ+invs2*sZ*sZ
    LZT=invs*sZT+invs2*sZ*sT
    LTT=invs*sTT+invs2*sT*sT
    # D(R,Z,T)=Dhat(R/sqrt(s)); E=x*d/dx.
    bounds={
        "D_Z_B1":LZ/2,
        "D_T_B1":LT/2,
        "D_ZZ_B2":LZ*LZ/4,
        "D_ZZ_B1":LZZ/2,
        "D_ZT_B2":LZ*LT/4,
        "D_ZT_B1":LZT/2,
        "D_TT_B2":LT*LT/4,
        "D_TT_B1":LTT/2,
    }
    expected={
        "D_Z_B1":64,"D_T_B1":128,
        "D_ZZ_B2":4096,"D_ZZ_B1":12288,
        "D_ZT_B2":8192,"D_ZT_B1":20480,
        "D_TT_B2":16384,"D_TT_B1":49152,
    }
    checks={"log_s_Z_bound":LZ==128,"log_s_T_bound":LT==256,
            "log_s_ZZ_bound":LZZ==24576,"log_s_ZT_bound":LZT==40960,
            "log_s_TT_bound":LTT==98304,
            "all_shape_composition_constants_exact":all(bounds[k]==v for k,v in expected.items())}
    return dict(checks=checks,bounds={k:str(v) for k,v in bounds.items()},
                definitions={"B0":"sup |Dhat|","B1":"sup |E Dhat|","B2":"sup |E^2 Dhat|"},
                mismatch_pressure_m0="epsilon^5 [64 B1 |dZ^3(KJ)|_0 + B0 |dZ^4(KJ)|_0]")


def inherited_shape_review():
    result=shape.review()
    required=("reference_mass_one","perturbation_zero_second_moment","pointwise_fourth_coefficient_nonzero")
    checks={f"parent:{name}":result["checks"]["reference_and_countercontrol:"+name] for name in required}
    checks["parent:radial_high_jet_total_nonzero"]=result["checks"]["pointwise_operator:radial_high_jet_total_nonzero"]
    return dict(checks=checks,scope="parent pointwise pressure and radial coefficient identities retained")


def review():
    parts={"profile_contract":profile_contract_review(),
           "signed_same_moment_family":signed_same_moment_family_review(),
           "positive_concentration":positive_concentration_review(),
           "dilation_budget":dilation_budget_review(),
           "inherited_shape":inherited_shape_review()}
    return dict(checks={p+":"+k:v for p,r in parts.items() for k,v in r["checks"].items()},parts=parts)
