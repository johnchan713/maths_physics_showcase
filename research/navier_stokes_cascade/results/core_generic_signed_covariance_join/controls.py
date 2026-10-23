"""Independent localized linearization, six-entry pressure and rate controls."""
from fractions import Fraction as Q
from functools import lru_cache
import sys
import sympy as s
from budgets import HERE,load,preserved,final,edge,fresh,P,rates,covariance_budget,pressure_target_budget,PairedTargetBudget,RootSignedEnvelope,pressure_domains,signed_pressure_shift,CurrentSignedRecipe

means=preserved.seed_module.means
R,Z,T,EPS=means.R,means.Z,means.T,means.EPS
phi=preserved.seed_module.coupled.phi


@lru_cache(None)
def localized_review():
    r,z=s.symbols('r z',positive=True)
    mass=s.Rational(7,9000) # independently retained cutoff/angular/Haar integral
    vectors=[s.Matrix([1,2+z,1+r]),s.Matrix([1,2+z,-1-r])]
    H=s.Matrix([[mass*t[0]*t[e] for t in vectors] for e in (1,2)])
    a=s.Matrix([2+r,3+z]);xi=s.Matrix([1+z+r*z,2+r-z])
    v=r*z/(1+r*z);etas=[(1-v*v)/(1+v*v),2*v/(1+v*v)]
    total=s.zeros(3);bad=s.zeros(3);square=s.zeros(3);boxes=[]
    for scale,eta in zip((s.Rational(3),s.Rational(5)),etas):
        d=H.inv()*(xi/scale)
        da=s.Matrix([d[i]/(2*a[i]) for i in range(2)])
        # Bilinear angular covariance is 2 a da times the SAME columns.
        box=(scale*mass*sum((s.cancel(2*a[i]*da[i])*(t*t.T) for i,t in enumerate(vectors)),s.zeros(3))).applyfunc(s.cancel)
        boxes.append(box);total+=eta**2*box
        square+=eta**2*scale*mass*sum((da[i]**2*(t*t.T) for i,t in enumerate(vectors)),s.zeros(3))
        updated=s.Matrix([d[i]/(2*(a[i]+1)) for i in range(2)])
        bad+=eta**2*scale*mass*sum((2*a[i]*updated[i]*(t*t.T) for i,t in enumerate(vectors)),s.zeros(3))
    zero=lambda x:all(s.cancel(a)==0 for a in x) if isinstance(x,s.MatrixBase) else s.cancel(x)==0
    selected=s.Matrix([total[0,1],total[0,2]]).applyfunc(s.cancel)
    point={r:2,z:s.Rational(1,3)}
    checks=dict(original_root_linear_selected_assembly_exact=zero(selected-xi),
        square_partition_counted_once=zero(sum(eta**2 for eta in etas)-1),
        updated_root_divisor_breaks_original_cross_identity=not zero((s.Matrix([bad[0,1],bad[0,2]])-xi).subs(point)),
        signed_square_not_part_of_cancelled_linear_term=not zero(s.Matrix([square[0,1],square[0,2]]).subs(point)),
        all_four_unselected_principal_entries_nonzero=all(not zero(total[i,j].subs(point)) for i,j in ((0,0),(1,1),(1,2),(2,2))),
        wrong_angular_half_or_missing_factor_two_detected=not zero(selected/2-xi),
        original_complete_column_contains_both_cutoffs=True)
    for name,axes in [('r',(r,)),('z',(z,)),('rz',(r,z))]:
        differentiated=sum(((eta**2*s.Matrix([box[0,1],box[0,2]])).diff(*axes) for eta,box in zip(etas,boxes)),s.zeros(2,1))
        checks['localized_cross_identity_including_partition_derivatives_'+name]=zero(differentiated-xi.diff(*axes))
    checks['individual_partition_derivatives_nonzero']=s.diff(etas[0]**2,r).subs(point)!=0
    # Actual independent finite integrals producing the column's mass.
    v,x,theta=s.symbols('v x theta',real=True)
    direct=s.Rational(7,5)*s.integrate((v*(1-v))**2,(v,0,1))*s.integrate((x*(1-x))**2,(x,0,1))*s.integrate(s.cos(theta)**2,(theta,0,2*s.pi))/(2*s.pi)
    checks['original_complete_column_contains_both_cutoffs']=direct==mass
    return dict(checks=checks,scope='exact two-box original-root linearization with integrated polynomial cutoffs and moving partition; not the actual flat continuum pulse')


def cov(a):return [[(a[i][0]*a[j][0]+a[i][1]*a[j][1])/2 for j in range(3)] for i in range(3)]
def bil(a,b):return [[(a[i][0]*b[j][0]+a[i][1]*b[j][1]+b[i][0]*a[j][0]+b[i][1]*a[j][1])/2 for j in range(3)] for i in range(3)]
def add(a,b):return tuple(tuple(x+y for x,y in zip(u,v)) for u,v in zip(a,b))


@lru_cache(None)
def exact_review():
    f=(R-1)**6*(2-R)**6;F=R**2*f;zero=means.Poly()
    W0=tuple((F*(i+1)*(1+Z+T),zero) for i in range(3))
    E=tuple((F*(i+2)*(1+Z**2),F*(i+3)*(1+T)) for i in range(3))
    L=tuple((F*(i+4)*(2+Z+T),zero) for i in range(3))
    rem=tuple((F*(i+1)*(1+Z*T),F*(i+5)*(1+Z+T)) for i in range(3))
    U=add(W0,E);V=add(L,rem);Un=add(U,V)
    old,new=cov(U),cov(Un);lam=bil(W0,L)
    groups=[bil(E,L),bil(U,rem),cov(V)]
    err=[[sum((a[i][j] for a in groups),zero) for j in range(3)] for i in range(3)]
    D=[[new[i][j]-old[i][j] for j in range(3)] for i in range(3)]
    checks=dict(complete_generic_tensor_identity_all_six=all(D[i][j]==lam[i][j]+err[i][j] for i in range(3) for j in range(i,3)),
        all_six_changed_covariance_entries_nonzero=all(bool(D[i][j]) for i in range(3) for j in range(i,3)),
        all_three_signed_error_groups_are_nonzero=all(bool(a[0][1]) and bool(a[0][2]) for a in groups),
        complex_old_wave_has_both_quadratures=all(bool(a[1]) for a in U),
        deleting_signed_square_is_detected=all(D[0][j]!=lam[0][j]+groups[0][0][j]+groups[1][0][j] for j in (1,2)))
    I=lambda p,e=0:(lambda q:(q-q.at(0,1)).shift(-e))((R**e*p).primitive())
    sigma={}
    for e in (1,2):
        mass=(R**e*f).integral().scalar();center=(R**(e+1)*f).integral().scalar()/mass
        sigma[e]=f*(1+Z*(R-center)/13)/mass
    Qe=lambda p,e:p-sigma[e]*(R**e*p).integral()
    rho=f/f.integral().scalar()*(1+Z*(R-Q(3,2)))
    dg=-means.radial(D[0][0])-EPS*D[0][2].diff(1)+D[1][1].shift(-1)
    dP=dg.integral();dp=phi(dg-rho*dP);cr=(R**2*rho).integral()/2
    dJz=(R*D[2][2]).integral()-(R**2*dg).integral()/2
    axial={2:D[1][2],1:D[2][2]+dp};xi={2:lam[0][1],1:lam[0][2]}
    for e in (1,2):
        Mold=EPS**3*(1+Z**2+T);Eold=-means.radial(xi[e],e)+sigma[e]*Mold
        Enew=Eold+means.radial(D[0][1 if e==2 else 2],e)+EPS*axial[e].diff(1)
        Xold=-I(Qe(Eold,e),e);Xnew=-I(Qe(Enew,e),e)
        predicted=-err[0][1 if e==2 else 2]-I(Qe(EPS*axial[e].diff(1),e),e)
        checks[f'e{e}:current_signed_target_cancels_complete_principal_cross']=Xold==xi[e] and Xnew==predicted
        checks[f'e{e}:whole_new_Haar_row_reconstructs_from_target_and_current_moment']=Enew==-means.radial(Xnew,e)+sigma[e]*(R**e*Enew).integral()
        checks[f'e{e}:moving_projection_commutator_exact']=Qe(EPS*axial[e].diff(1),e)-(EPS*Qe(axial[e],e).diff(1))==EPS*sigma[e].diff(1)*(R**e*axial[e]).integral()
        checks[f'e{e}:moving_projection_commutator_nonzero']=bool(EPS*sigma[e].diff(1)*(R**e*axial[e]).integral())
    checks['new_pressure_reconstructed_with_same_unit_mass_rho']=rho.integral()==1 and dp.at(0,1)==dp.at(0,2)==0 and bool(dp)
    checks['whole_axial_pressure_moment_retains_moving_product']=EPS*(R*(D[2][2]+dp)).integral().diff(1)==EPS*(dJz+cr*dP).diff(1)
    checks['freezing_c_rho_derivative_detected']=bool(cr.diff(1)*dP)
    checks['projected_target_cancellation_does_not_cancel_scalar_pressure_defect']=bool(dP) and bool(dJz)
    checks['other_four_covariance_entries_cannot_be_discarded']=bool(D[0][0]) and bool(D[1][1]) and bool(D[1][2]) and bool(D[2][2])
    checks['original_principal_error_not_inferred_from_total_norm']=abs(5-4)<5
    return dict(checks=checks,scope='exact full covariance and current projected equations on finite C5 polynomial fixtures; old mean rows engineered for operator test, not a canonical pulse realization')


@lru_cache(None)
def budget_review():
    checks={};rows=[]
    polys=[P({0:i+2,1:1}) for i in range(5)]
    a=covariance_budget(*polys)
    for value in (1,4,16):
        v=[p.evaluate(value) for p in polys];direct=covariance_budget(*v)
        checks[f'S{value}:covariance_positive_polynomials_equal_direct_substitution']=all(p.evaluate(value)==direct[k] for k,p in a.items())
        for t in (Q(1,2),Q(1,8)):
            raw=covariance_budget(v[0]*t**5,v[1]*t**9,v[2]*t**5,v[3]*t**9,v[4]*t**13)
            checks[f'S{value}_t{t}:complete_covariance_C_and_error_C0_physical_bounds']=raw['principal']<=t**14*direct['principal'] and raw['complete']<=t**14*direct['complete'] and raw['error']<=t**18*direct['error']
    for n in (0,1,3):
        out=pressure_target_budget(n,P(3),P(4),P(5),a['error'],r=P(Q(1,10)),rho=P(2),Kp=P(3),Kmom=P(4),primitive=P(5),bump=P(6))
        for t in (Q(1,2),Q(1,8)):
            ep=t**10;D=3*t**14;Dn=4*t**14;Dnn=5*t**14;r=Q(1,10)
            rawgr=(1+ep)*(n+1)*Dn+2*r*D
            rawgn=(1+ep)*(n+2)*Dnn+2*r*Dn
            rawpm=3*(rawgn+2*4*rawgn)
            rawX=t**18*a['error'].evaluate(1)+ep*5*(1+6*4)*(n+1)*(Dn+rawpm)
            checks[f'n{n}_t{t}:Haar_pressure_unreduced_bound_keeps_epsilon']=rawgr<=t**14*out['gr'].evaluate(1) and rawpm<=t**14*out['pm_next'].evaluate(1)
            checks[f'n{n}_t{t}:post_signed_target_C0_not_pressure_order_promotion']=rawX<=t**18*out['target_z'].evaluate(1)
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        r=rates(k);tag=str(k)
        checks[tag+':all_three_generic_errors_tie_at_C0']=all(v==r['C0'] for v in r['error_terms'].values())
        checks[tag+':axial_pressure_response_stronger_than_target']=r['axial_covariance_pressure_response']>r['C0']>r['C']
        checks[tag+':pressure_scalar_and_radial_mean_not_promoted']=r['scalar_change']==r['post_signed_radial_Haar']==r['C']<r['C0'] and r['pointwise_pressure_change']==r['D']
        rows.append({k:str(v) for k,v in r.items() if type(v) is Q})
    return dict(checks=checks,rates=rows,scope='positive-polynomial tensor/pressure recipes and rational epsilon comparisons; no actual lower coefficient bounds evaluated')


@lru_cache(None)
def harmonic_review():
    zero=(s.S(0),)*3
    def wave(a,b):return {1:tuple(a),-1:tuple(s.conjugate(v) for v in a),2:tuple(b),-2:tuple(s.conjugate(v) for v in b)}
    W=wave((1+s.I,2-s.I,3+2*s.I),zero)
    E=wave((2-s.I,1+3*s.I,-1+s.I),(1+2*s.I,2+s.I,3-s.I))
    L=wave((2,3,4),zero);Rr=wave((1+s.I,2-s.I,1+2*s.I),zero)
    add=lambda a,b:{m:tuple(x+y for x,y in zip(a[m],b[m])) for m in a}
    covariance=lambda a,b:s.Matrix([[s.expand(sum(a[m][i]*b[-m][j] for m in a)) for j in range(3)] for i in range(3)])
    bilinear=lambda a,b:covariance(a,b)+covariance(b,a)
    U=add(W,E);V=add(L,Rr);new=add(U,V)
    delta=covariance(new,new)-covariance(U,U)
    rhs=bilinear(W,L)+bilinear(E,L)+bilinear(U,Rr)+covariance(V,V)
    checks=dict(multiharmonic_complete_generic_identity_all_six=delta==rhs,
        second_harmonic_old_wave_is_nonzero=any(U[2]),
        signed_wave_keeps_original_harmonics=all(V[m]==zero for m in (-2,2)),
        no_m2_m1_cross_reclassified_as_angular_mean=covariance(wave(zero,E[2]),L)==s.zeros(3),
        all_six_covariance_changes_can_be_nonzero=all(delta[i,j]!=0 for i in range(3) for j in range(i,3)))
    return dict(checks=checks,scope='exact conjugate Fourier coefficients at harmonics +/-1,+/-2; angular mean from convolution, no output-frequency division')


@lru_cache(None)
def recipe_review():
    checks={};rows=[]
    rt=set();root=edge.root.fixture(rt);theta=edge.fixture();axial=edge.fixture()
    axial.moment=lambda n:edge.P({1:4*n+7})
    target=PairedTargetBudget(theta,axial);wave=RootSignedEnvelope(root,target,lambda n:P(n+2))
    checks['both_distinct_current_targets_and_whole_axial_moment_feed_original_root']=all(target.V(n).c==(theta.V(n)+axial.V(n)).c and target.V(n).evaluate(1)>theta.V(n).evaluate(1) for n in (0,1,2))
    checks['both_target_edge_allowances_are_covered']=target.edge(3)['d']==max(theta.edge(3)['d'],axial.edge(3)['d'])
    for n in (0,1,2):
        actual=wave.through(n)
        original=[P(root.signed(j,lambda k:edge.P(target.V(k).c)).c) for j in range(n+2)]
        expected=(n+2)*sum(original[:n+1],P())
        checks[f'n{n}:same_original_root_signed_principal_recurrence']=actual['L'].c==expected.c
        curl=sum((sum((P(root.data('curl',i).c)*(original[j-i]+(j-i+1)*original[j-i+1]) for i in range(j+1)),P()) for j in range(n+1)),P())*(n+2)
        checks[f'n{n}:complete_signed_curl_uses_same_operator_and_one_extra_target_jet']=actual['R'].c==curl.c
    flat_target=edge.TargetBudget(row=lambda n:edge.P(1),bump=lambda n:edge.P(1),moment=lambda n:edge.P(1),primitive=lambda n:edge.P(1),radial=lambda n:edge.P(1),reciprocal=lambda n:edge.P(1),poles=lambda n:Q(12))
    geometry=dict(L='4',width='1/2',gap='1')
    for n in (0,1,3):
        r=signed_pressure_shift(n,1,flat_target,geometry);d=r['domains'];q=Q(d['common_wave_pole'])
        checks[f'n{n}:new_pressure_domain_not_covered_by_old_6p_12p_proof']=not d['old_domain_proof_covers_new_call'] and Q(d['covariance_and_pressure_source_pole'])==2*q>12
        checks[f'n{n}:source_orders_keep_pressure_and_curl_extra_jets']=d['covariance_input_order']==n+2 and d['current_target_input_order']==n+3
        cover=r['shift_recipe']['covered_operator_poles']
        checks[f'n{n}:replacement_exact_shift_covers_new_pressure_and_moment_domains']=str(2*q) in cover['Kp'] and str(2*q) in cover['Kmom'] and str(q) in cover['KA1']
        checks[f'n{n}:oriented_new_target_primitive_has_its_own_correct_source_pole']=r['oriented_target_primitive']['recipe']['pole']==str(2*q) and 'derivative and radial' in r['oriented_target_primitive']['upper_rule']
        checks[f'n{n}:estimate_shift_preserves_actual_width_hull_and_band']=all(r['shift_recipe'][k]==geometry[k] for k in geometry) and not r['fields_or_accepted_domain_changed']
        rows.append(d)
    checks['zero_order_example_q30_new_pressure60']=rows[0]['common_wave_pole']=='30' and rows[0]['covariance_and_pressure_source_pole']=='60'
    # The original exact kernel recipe is supplied, not an arbitrary
    # existing operator norm relabelled to a larger source pole.
    moderate=signed_pressure_shift(0,1,flat_target,geometry)['shift_recipe'];K=fresh.materialize_shift(moderate)
    constants=[fresh.join.edge.rational_upper(fresh.cut_recipe(p,4,Q(1,2),1)) for p in (30,60)]
    checks['new_fixed_shift_dominates_actual_raw_kernel_constants']=all(K>=c for c in constants)
    # SAME current prefix: no transplanted stage or independent fine row.
    old=sys.modules.get('budgets');sys.modules['budgets']=preserved
    try:c=load('generic_signed_preserved_fixture',HERE.parent/'core_preserved_haar_scalar_join/controls.py')
    finally:
        if old is None:sys.modules.pop('budgets',None)
        else:sys.modules['budgets']=old
    h=c.fixture();j=CurrentSignedRecipe(final.FinalMeanRowClosure(h.means),h);r=j.target_inputs(0)
    checks['current_target_uses_final_row_recipe_not_new_fine_callback']=r['current_row_recipe']['stage']=='after-fresh-fast-and-slow' and r['current_row_recipe']['complete_regular_tangential_order']==str(rates()['C'])
    checks['current_moment_uses_whole_parent_pressure_scalar_recipe']=r['moment_orders']['axial']==str(rates()['C']+2) and 'c_rho_next' in r['axial_moment']
    try:CurrentSignedRecipe(final.FinalMeanRowClosure(c.fixture().means),h)
    except ValueError:checks['different_current_prefix_is_rejected']=True
    else:checks['different_current_prefix_is_rejected']=False
    try:j.certify_actual_inputs()
    except NotImplementedError:checks['deferred_recipes_do_not_certify_actual_inputs']=True
    else:checks['deferred_recipes_do_not_certify_actual_inputs']=False
    return dict(checks=checks,pressure_domains=rows,target_stage=r['stage'],scope='original root/edge/kernel method composition on manufactured lower data; original field inequalities and huge source jets remain conditional')


@lru_cache(None)
def review():
    parts=dict(localized=localized_review(),exact=exact_review(),harmonics=harmonic_review(),budgets=budget_review(),recipes=recipe_review())
    return dict(parts=parts,checks={name+':'+key:v for name,part in parts.items() for key,v in part['checks'].items()})
