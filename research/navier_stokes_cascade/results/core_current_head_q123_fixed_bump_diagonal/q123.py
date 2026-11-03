"""Existing fixed repair bumps -> shifted radial bound -> current diagonal bits."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
repair=load('q123_repair','core_current_head_q122_positive_order_repair/q122.py')
old=load('q123_old_bumps','cutoff_higher_jet_interface/jet_bounds.py')
order=repair.order;exact=repair.exact
q121=repair.current;q120=q121.current;lg=q120.borel.ceil_log2
IDENTITY='same-squared-flat-dy-unit-mass-width-1/10-bumps'
LOW=(320,121600,220000000,18870643200000)

def bump_bits(m):
 order(m)
 if m<=3:return lg(Q(LOW[m]))
 j=m+1
 return 8+15*j+3*j*lg(Q(j))

def geometry(p,rp,identity):
 p,rp=map(exact,(p,rp))
 if identity!=IDENTITY:raise ValueError('Original retained repair bump identity required')
 if min(p,rp)<=0 or rp*rp>2*p:raise ValueError('Positive support X/radius lower bounds required')
 return max(Q(1),1/p),rp

def radial_bits(m,p,rp,identity):
 order(m);a,rp=geometry(p,rp,identity)
 fact=m*lg(Q(m)) if m else 0
 return dict(U=bump_bits(m)+fact+m*lg(a),
             E_over_R=bump_bits(m)+fact+m*lg(Q(3,2)*a)+lg(max(Q(1),1/rp)),
             factorials_materialized=False,actual_support_identity_certified=False)

def moderate_radial(m,p,rp,identity,maximum=20):
 order(m);order(maximum);a,rp=geometry(p,rp,identity)
 if m>maximum:raise ValueError('Retain high-order factorials as budgets')
 B=LOW[m] if m<=3 else 162*32600**(m+1)*factorial(m+1)**3
 return dict(U=a**m*factorial(m)*B,E_over_R=(Q(3,2)*a)**m*factorial(m)*B/rp)

def current_row(h,n,m,Xhi,Rlo,Rhi,rho,lam,e_star,discrepancies,
                Fcut_bits,Ucut_bits,p,rp,patch_identity,bump_identity,log2_lower):
 order(Fcut_bits);order(Ucut_bits)
 # Q121 validates the same hull/h/n/m, with zero placeholders for profile norms.
 r=q121.recipe(h,n,m,Xhi,Rlo,Rhi,0,0)
 c=repair.coefficient_bound(n,m+1,rho,lam,e_star,discrepancies,patch_identity)
 cbits=lg(max(Q(1),c['coefficient_raw_upper']))
 F=radial_bits(m,p,rp,bump_identity)['E_over_R']
 U=radial_bits(m+1,p,rp,bump_identity)['U']
 Fbits=1+max(Fcut_bits,F+cbits);Ubits=1+max(Ucut_bits,U+cbits)
 h,Xhi,Rlo,Rhi=map(exact,(h,Xhi,Rlo,Rhi))
 cgeom=Q(1,2)-h+2*n*h
 TF=4*Rhi;TU=2+(112+32*cgeom)*max(Q(1),Xhi)/((1-2*h)*Rlo)
 Tbits=1+max(lg(max(Q(1),TF))+Fbits,lg(max(Q(1),TU))+Ubits)
 tuple_bits=Tbits+(m*lg(Q(m)) if m else 0)+m*lg(r['geometric_base'])
 row=q120.deferred_row(h,n,m,tuple_bits,log2_lower)
 return dict(discrepancy_order=m+1,Fcut_order=m,Ucut_order=m+1,
             coefficient_binary_upper=cbits,Fraw_binary_upper=Fbits,Uraw_binary_upper=Ubits,
             tuple_binary_upper=tuple_bits,row=row,factorials_materialized=False,
             actual_seed_discrepancy_identity_support_or_scale_inputs_certified=False)

def record():
 return dict(status='fixed-bump-to-current-diagonal-binary-budget',
   inherited_bump='Q63: 162*32600^(m+1)*((m+1)!)^3; sharper inherited low rows retained',
   U_through='max(1,p^-1)^m*m!*B_m',
   E_over_R_through='rp^-1*(3/2*max(1,p^-1))^m*m!*B_m',
   current_join='Q122 repair -> Q121 raw tuple -> Q120 retained diagonal; no new scales',
   actual_profile_support_scale_cone_pulse_residual_or_blowup_certified=False)
