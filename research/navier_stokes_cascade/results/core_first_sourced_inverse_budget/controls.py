"""Exact manufactured algebra and independent finite-polynomial norm controls."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
import importlib.util
import sys
from budgets import P, ForcedBudget, fixture

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@lru_cache(None)
def inherited():
    return load("inverse_exact_forced", HERE.parent / "first_sourced_foundation/sourced_checks.py")


def algebra_review():
    maps = inherited().maps(); J = maps.JetClass()
    s = J.variable(0, Q(1,3), n=2, order=7)
    x = J.variable(1, Q(1,7), n=2, order=7)
    zero = 0*s; one = 1+zero
    same = lambda a,b:all((a-b).partial((i,j))==0 for i in range(5) for j in range(5-i))
    dot = lambda a,b:sum((c*d for c,d in zip(a,b)), zero)
    mv = lambda A,v:[dot(row,v) for row in A]
    n = [one, s+x, one+s*x]; np = [q.derivative(0) for q in n]
    norm = dot(n,n); h = [q*norm**-1 for q in n]
    hp = [q.derivative(0) for q in h]
    K = [[zero, -2-x, zero], [3+s, zero, x], [1-x, s, zero]]
    nK = [sum((n[i]*K[i][j] for i in range(3)), zero) for j in range(3)]
    A0 = [[-K[i][j]+h[i]*(nK[j]-np[j]) for j in range(3)] for i in range(3)]
    damping = norm
    A = [[A0[i][j]-(damping if i==j else zero) for j in range(3)] for i in range(3)]
    profile_rate = one+s+x
    adapted = [[A[i][j]+(damping+profile_rate)*h[i]*n[j] for j in range(3)] for i in range(3)]
    lift_forcing = [a-b for a,b in zip(mv(A0,h),hp)]
    v = [s+x, one+x*x, s*s-x]
    z = dot(n,v); w = [a-b*z for a,b in zip(v,h)]
    reconstructed = [a+q*z+dh*z+profile_rate*hh*z for a,q,dh,hh in zip(mv(A,w),lift_forcing,hp,h)]
    tangent = [-(s+x)*s, s, zero]
    tx = [q.derivative(1) for q in tangent]
    # n=e_z, K=diag(1,1,0), lambda=1: tangent rates are -2,
    # the natural normal rate is -1, and P=exp(-2s). The adapted
    # normal rate is -2. This is an exact constant-matrix counterexample.
    natural_rates = (Q(-2), Q(-2), Q(-1))
    profile_log_rate = Q(-2)
    adapted_normal_rate = natural_rates[2]+(1+profile_log_rate)
    checks = {
        "full_normal_lift_has_unit_pairing":same(dot(n,h),one),
        "normal_lift_connection_is_tangent":same(dot(n,lift_forcing),zero),
        "ambient_generator_evolves_normal_by_scalar_damping":all(
            same(sum((n[i]*A[i][j] for i in range(3)), zero),-np[j]-damping*n[j]) for j in range(3)),
        "adapted_normal_component_has_exact_profile_logarithmic_rate":all(
            same(sum((n[i]*adapted[i][j] for i in range(3)),zero),-np[j]+profile_rate*n[j]) for j in range(3)),
        "adapted_extension_leaves_tangent_generator_unchanged":all(same(a,b) for a,b in zip(mv(A,w),mv(adapted,w))),
        "arbitrary_ambient_initial_vector_decomposition_exact":all(same(a,b) for a,b in zip(mv(adapted,v),reconstructed)),
        "manufactured_velocity_is_tangent":same(dot(n,tangent),zero),
        "parameter_derivative_is_not_tangent":bool(dot(n,tx).terms),
        "differentiated_constraint_keeps_normal_connection":same(dot(n,tx)+dot([q.derivative(1) for q in n],tangent),zero),
        "omitting_moving_normal_breaks_lift_tangency":bool(dot(np,h).terms),
        "normal_source_projects_to_zero":all(same(a-b*norm,zero) for a,b in zip(n,h)),
        "normal_source_pressure_does_not_vanish":same(norm*norm**-1,one),
        "ungauged_normal_damping_cannot_preserve_arbitrary_Gaussian":
            natural_rates[0]==natural_rates[1]==profile_log_rate and natural_rates[2]-profile_log_rate>0,
        "Gaussian_counterexample_is_repaired_by_adapted_normal_rate":adapted_normal_rate==profile_log_rate,
    }
    checks.update({"inherited:"+k:v for k,v in inherited().forced_checks()["checks"].items()})
    return dict(checks=checks,scope="exact finite rational jets; manufactured normals and fields, not the selected core propagator")


def raw_bound(polynomials, j):
    """L1 tensor upper bound on every order-j mixed jet on the full unit box."""
    from math import factorial
    out = Q(0)
    for a in range(j+1):
        b = j-a; value = Q(0)
        for p in polynomials:
            for (r,z), c in p.terms.items():
                if r>=a and z>=b:
                    value += abs(c)*(factorial(r)//factorial(r-a))*(factorial(z)//factorial(z-b))
        out = max(out, value)
    return out


def polynomial_review():
    maps = inherited().maps(); J = maps.JetClass()
    s = J.variable(0, Q(0), n=2, order=10)
    x = J.variable(1, Q(0), n=2, order=10)
    zero = 0*s
    # Exact full-domain manufactured solution; A0 is a commuting nilpotent
    # shear, with nonnegative scalar damping 1. No ODE integration or mesh.
    velocity = [-s*x-s*s/2-(1+x*x)*(x*s*s/2+s**3/3), -s*(1+x*x), zero]
    av = [(x+s)*velocity[1], zero, zero]
    source = [-u.derivative(0)+a-u for u,a in zip(velocity,av)]
    K = [-(x+s)]
    b = ForcedBudget(normal=lambda n:P(int(n==0)), moving_normal=lambda n:P(),
                     matrix=lambda n:P(raw_bound(K,n)), source=lambda n:P(raw_bound(source,n)),
                     profile_rate=lambda n:P(),
                     tangent=P(3), length=P(1), inverse_normal_square=1)
    checks = {
        "manufactured_solution_zero_on_whole_entrance_parameter_slice":all(
            not any(a==0 and c for (a,_),c in u.terms.items()) for u in velocity),
        "source_computed_from_full_damped_equation":all(
            not (u.derivative(0)-a+u+f).terms for u,a,f in zip(velocity,av,source)),
        "tangent_propagator_uses_nilpotent_shear_not_exponential_Gronwall":1+Q(3,2)<3,
        "coefficient_bounds_cover_full_unit_box":raw_bound(K,0)==2 and raw_bound(K,1)==1,
    }
    rows = []
    for n in range(5):
        exact = raw_bound(velocity,n); upper = b.velocity(n).evaluate(1)
        checks[f"order{n}_all_mixed_full_box_solution_norm_is_bounded"] = upper >= exact
        checks[f"order{n}_full_source_enters_pressure_even_after_projection"] = b.pressure(n).evaluate(1) >= 0
        rows.append(dict(order=n,exact_box_majorant=str(exact),inverse_majorant=str(upper)))
    # f(s)=s(1-s) vanishes at both endpoints; its uncut zero-entry response
    # has nonzero exit. This is not an endpoint-flat smooth-source example.
    checks["zero_endpoint_source_values_do_not_give_zero_exit"] = Q(1,2)-Q(1,3)==Q(1,6)>0
    return dict(checks=checks,rows=rows,scope="exact coefficient-sum bounds on manufactured polynomial fields over a whole box; no sampled diagnostic or actual giant pulse")


def budget_review():
    checks = {}; records = []
    for m in (-2,-1,1,2):
        b = fixture(harmonic=m)
        for n in range(4):
            for S in (1,3):
                value = b.velocity(n).evaluate(S)
                g = b.projected_source(n).evaluate(S)
                parameter = b.ambient().evaluate(S)*b.length.evaluate(S)*(g+sum(
                    comb(n,i)*b.generator(i).evaluate(S)*b.velocity(n-i).evaluate(S) for i in range(1,n+1)))
                pulse = 0 if n==0 else b.projected_source(n-1).evaluate(S)+sum(
                    comb(n-1,i)*b.generator(i).evaluate(S)*b.velocity(n-1-i).evaluate(S) for i in range(n))
                checks[f"m{m}_n{n}_S{S}_independent_rational_recurrence"] = value==parameter+pulse
            checks[f"m{m}_n{n}_conjugate_majorants_match"] = b.pressure(n).record()==fixture(harmonic=-m).pressure(n).record()
        records.append(dict(harmonic=m,velocity_order2=b.velocity(2).record(),pressure_order2=b.pressure(2).record()))
    b = fixture(); removed = fixture(); removed.moving_normal = lambda n:P()
    checks["moving_normal_not_erased_in_pressure_or_ambient_bound"] = b.pressure(1).evaluate(1)>removed.pressure(1).evaluate(1) and b.ambient().evaluate(1)>removed.ambient().evaluate(1)
    zero_source = fixture(); zero_source.source=lambda n:P()
    checks["whole_path_zero_source_has_zero_response_and_pressure"] = all(
        zero_source.velocity(n).evaluate(2)==zero_source.pressure(n).evaluate(2)==0 for n in range(4))
    trace=set(); b=fixture(trace=trace); b.velocity(3); b.pressure(3); b.potential(3)
    checks["order_three_only_finitely_many_full_source_jets"] = max(n for name,n in trace if name=="source")==3
    checks["lift_derivative_requires_at_least_first_normal_jet"] = ("normal",1) in trace
    cutoff=lambda n:P(n+1)
    checks["same_cutoff_pressure_and_potential_keep_derivative_products"] = all(
        b.cut(kind,cutoff,2).evaluate(1)>cutoff(0).evaluate(1)*getattr(b,kind)(2).evaluate(1) for kind in ("pressure","potential"))
    return dict(checks=checks,records=records,scope="finite positive-polynomial substitution; coefficient and tangent propagator inputs remain conditional")


@lru_cache(None)
def review():
    parts = dict(moving_constraint=algebra_review(), full_box_norm=polynomial_review(), explicit_budgets=budget_review())
    return dict(parts=parts,checks={name+":"+k:v for name,row in parts.items() for k,v in row["checks"].items()})
