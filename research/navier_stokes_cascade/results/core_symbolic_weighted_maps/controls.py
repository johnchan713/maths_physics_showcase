from bounds import ExpSum,weighted_maps,fast_coefficients,edge_absorption,rows,Q
def review():
 import sympy as s
 checks={}
 for p in (0,3,8):
  for gap in (Q(1,4),Q(1)):
   w=Q(1,2);L=Q(4);v=weighted_maps(p,L,w,gap,jacobian=1,KT_profile=2,KA_profile=3)
   old=rows.weighted_operator_budget(p,L,w,rows.P(2),rows.P(3),gap)
   # Sum upper bound >= maximum upper bound, termwise exp envelope identical.
   checks[f'pole{p}_gap{gap}_historical_KT_enclosed']=v['KT'].rational_upper()>=old['KT'].evaluate(1)
   checks[f'pole{p}_gap{gap}_historical_KA_enclosed']=v['KA'].rational_upper()>=old['KA'].evaluate(1)
   checks[f'pole{p}_gap{gap}_same_collar']=v['proof_collar']==min(Q(1),L/4,gap)
 a=ExpSum([(2,1),(3,1)]);b=ExpSum([(7,2)])
 checks['exact_exponential_addition']=a.terms=={Q(1):Q(5)}
 checks['exact_exponential_product']=(a*b).terms=={Q(3):Q(35)}
 v=weighted_maps(8,4,Q(1,10**20),Q(1,10**4),jacobian=2,KT_profile=3,KA_profile=4)
 checks['original_tiny_width_kept_symbolic']=len(v['KT'].record())<=3 and all(c>0 for c in v['KT'].terms.values())
 try:v['KT'].rational_upper();refuses=False
 except ValueError:refuses=True
 checks['giant_diagnostic_evaluation_refused_without_rejecting_recipe']=refuses
 # Fixed-a translation: U derivatives do not differentiate M*a.
 U,y,M,a0=s.symbols('U y M a');F=(U+a0)**3*(y+M*a0)**2
 checks['fixed_shift_coefficient_derivative_adds_no_M']=s.diff(F,U)==3*(U+a0)**2*(y+M*a0)**2
 # Graph derivative cancellation is lost if an unshifted primitive is used.
 checks['shifted_graph_transport_not_replaced_by_unshifted']=s.diff(F,U)+M*s.diff(F,y)!=s.diff(F,U)
 f=fast_coefficients(1,ExpSum([(2,0)]),ExpSum([(3,0)]),Et_raw=5,Ez_raw=7,Ez_radial_raw=11,time_floor=6)
 gain=s.binomial(6,5)*24*2**5*6
 checks['actual_axial_A1_term_retained']=f['Ft'].rational_upper()==gain*(5+4*7)
 checks['radial_uses_next_source_order']=f['Fr'].rational_upper()==gain*2*2*11
 checks['source_orders_and_external_scale_retained']=f['raw_source_orders']==(5,6) and f['external_S_power']==1
 return dict(checks={k:bool(v) for k,v in checks.items()},scope='exact finite symbolic recipes and manufactured shift controls; no actual current source or continuum certification',parent_commit='6d77c8286d38d94b380bd264fdc06158263960c6')
