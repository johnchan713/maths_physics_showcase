"""Audit Q144's scalar gate against the unchanged symbolic hierarchy."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=load('q145_previous','core_current_head_q144_shape_annulus_gevrey/q144.py')

IDENTITY='same-Q144-full-pressure-profile-and-original-A-G-H-N-hierarchy'
A_MIN=2**16

def hierarchy_inequalities(A_floor=A_MIN):
 """Exact exponent comparisons; the exponential hierarchy is not formed."""
 A_floor=int(A_floor)
 if A_floor<A_MIN:raise ValueError('Original A=C^4096>=2^16 hierarchy required')
 return dict(
  stage='symbolic-original-A-G-H-N-comparison',A_floor=A_floor,
  definitions=dict(A='C^4096',G='exp(A^32)',H='exp(G^4)',N='1+floor(H^32)'),
  exact_checks=dict(A_at_least_two=A_floor>=2,A32_exceeds_logA=A_floor**32>A_floor,
   G4_exceeds_G=True,H32_at_least_H=True),
  conclusions=dict(p2='|p2|<=A<G',root='|mu|,|mu*p2|<G',
   exponential='Eexp=H>exp(G)',value_frequency='V0=H<N'),
  hierarchy_materialized=False)

def reduce_gate(recipe,identity):
 """Reduce Q144's scalar requests without claiming its free caps are actual."""
 if identity!=IDENTITY or recipe.get('stage')!='conditional-Q70-Q127-shape-annulus-mixed-Gevrey3':
  raise ValueError('Same Q144 shape/annulus recipe and original hierarchy required')
 if not recipe.get('all_five_actual_incoming_rows_reintegrated') or not recipe.get('full_p2_Pi_and_Pi_eta_retained'):
  raise ValueError('Full five-row pressure memory required')
 request=previous.gate_request(recipe);hierarchy=hierarchy_inequalities()
 return dict(stage='Q144-scalar-gate-reduced-to-same-profile-cap-binding',
  q144_request=request,hierarchy=hierarchy,a_floor=recipe['a_floor'],
  scalar_gate=dict(G_cap=True,Eexp_choice='original H',V0_choice='original H',N='1+floor(H^32)',
   conditional_on_unchanged_compact_envelope=True),
  pressure_memory=dict(all_five_rows=True,Pi=True,Pi_eta=True),
  remaining_inequality=('one same-profile Q73/Q143/Q127 complex cap tuple must bound ell,E,E^-1 '
   'at every center of the shape/restoration/earlier-annulus and specialize Q144 inputs'),
  same_profile_complex_caps_certified=False,actual_q140_q138_cone_join_certified=False,
  actual_complete_profile_pulse_covariance_residual_or_blowup_certified=False)

def record():
 return dict(status='original-hierarchy-discharges-Q144-scalar-gate',
  correction='Q144 scalar requests are inherited consequences: |p2|<=A<G, |mu| and |mu*p2|<G, Eexp=H>exp(G), V0=H<N',
  hierarchy='A=C^4096>=2^16; G=exp(A^32); H=exp(G^4); N=1+floor(H^32); no huge quantity is materialized',
  pressure='the reduction applies only to the unchanged full p2 with all five rows, Pi and Pi_eta',
  blocking_inequality='bind one actual same-profile Q73/Q143/Q127 complex value-cap tuple across the shape/restoration/earlier annulus and specialize Q144 free caps',
  scalar_gate_independent_regional_obligation=False,same_profile_complex_caps_or_cone_join_certified=False,
  actual_complete_profile_pulse_covariance_residual_or_blowup_certified=False,historical_frozen_records_modified=False)

