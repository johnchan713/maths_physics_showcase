"""The original positive-order LINEAR repair, conditional on full incoming moments."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q122_current',HERE.parent/'core_current_head_q121_all_order_radial_tuple/q121.py');current=importlib.util.module_from_spec(s);s.loader.exec_module(current)
order=current.order;exact=current.exact
IDENTITY='original-patch-e_star/(1+eta^2)'
KEYS={'m1','m2','m3','m4','m5'}

def normalizer_raw(k,raw,e_star):
 order(k);e=exact(e_star);raw=tuple(map(exact,raw))
 if e<=0 or len(raw)!=k+1 or min(raw)<0:raise ValueError('All raw discrepancy jets through k and positive e_star required')
 return (2*raw[k]+(2*k*raw[k-1] if k else 0)+(k*(k-1)*raw[k-2] if k>=2 else 0))/e

def coefficient_bound(n,m,rho,lam,e_star,discrepancies,identity):
 order(n);order(m)
 if n==0:raise ValueError('Positive-order linear coefficient repair only; not the leading quadratic root')
 if identity!=IDENTITY:raise ValueError('Exact original patch normalizer identity required')
 if set(discrepancies)!=KEYS:raise ValueError('All five full incoming moment discrepancy bounds required')
 rho,lam,e=map(exact,(rho,lam,e_star));D=max(map(exact,discrepancies.values()))
 if rho<=0 or not 0<lam<=Q(1,100) or e<=0 or min(map(exact,discrepancies.values()))<0:raise ValueError('Original positive patch parameters and nonnegative raw norms required')
 G=max(Q(1),rho**-2);P=(m*m+m+2)/e
 C=1000*G*D*(1+P)/lam
 return dict(coefficient=n,angular_order=m,coefficient_raw_upper=C,normalizer_raw_upper=P,
             all_five_rows=True,actual_full_pressure_discrepancy_bounds_certified=False,
             actual_patch_identity_certified=False,leading_root_or_support_certified=False)

def finalized_tuple(h,n,m,Xhi,Rlo,Rhi,rho,lam,e_star,discrepancies,
                    Fcut_raw,Ucut_raw,bump_E_over_R_raw,bump_U_raw,identity):
 # Discrepancies and U inputs through m+1; F/radial E-over-R bump through m.
 a=coefficient_bound(n,m+1,rho,lam,e_star,discrepancies,identity)
 F,U,BF,BU=map(exact,(Fcut_raw,Ucut_raw,bump_E_over_R_raw,bump_U_raw))
 if min(F,U,BF,BU)<0:raise ValueError('Same cut seeds and fixed disjoint-bump raw bounds required')
 C=a['coefficient_raw_upper']
 r=current.recipe(h,n,m,Xhi,Rlo,Rhi,F+BF*C,U+BU*C)
 return dict(repair=a,tuple_recipe=r,tuple_binary_upper=current.tuple_binary_upper(r),
             discrepancy_order=m+1,bump_F_order=m,bump_U_order=m+1,
             actual_seed_bump_moment_or_support_inputs_certified=False)

def record():
 return dict(status='all-order-positive-coefficient-repair-bound',
   raw_normalizer='[2*d_k+2*k*d_(k-1)+k*(k-1)*d_(k-2)]/e_star',
   coefficient_raw='1000*max(1,rho^-2)*max_i D_i*(1+(m^2+m+2)/e_star)/lambda_patch',
   original_linear_solve=True,all_five_full_pressure_rows_required=True,
   actual_profiles_schedule_support_cone_pulse_residual_or_blowup_certified=False)
