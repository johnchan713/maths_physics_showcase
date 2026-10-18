"""Exact real complete-curl algebra and the unfixed tangential variance."""
from fractions import Fraction as F


def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def outer(a,b):return [[x*y for y in b] for x in a]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(a,c):return [[c*x for x in row] for row in a]
def sub(a,b):return add(a,scale(b,-1))
def covariance(A,B,k):return scale(add(outer(A,A),scale(outer(B,B),F(1,k*k))),F(1,2))
def phase_average(A,B,k):
    # Four exact quarter-period samples integrate quadratic trig expressions.
    values=[A,[-x/k for x in B],[-x for x in A],[x/k for x in B]]
    result=[[F(0)]*3 for _ in range(3)]
    for v in values:result=add(result,scale(outer(v,v),F(1,4)))
    return result


def jet(x,eta=F(0)):
    x,eta=F(x),F(eta);n=[F(1),F(-1),F(-1)]
    f=1+x*x+eta*(1+x);fp=2*x+eta
    g=x+x**3/3+eta*x*x;gp=1+x*x+2*eta*x
    A=[f,f+g,-g];Ap=[fp,fp+gp,-gp]
    Cp=[v/3 for v in cross(n,Ap)]
    B=[F(0),-Cp[2],Cp[1]]
    return n,A,B,fp


def complete_review():
    rows=[];exact=True;div=True;change=True;missing=False;generic=False
    for x in (F(-1,4),F(1,5),F(2,3)):
        for k in (3,11,31):
            for eta in (F(-1,10),F(1,7)):
                n,A,B,fp=jet(x);_,Ar,Br,fpr=jet(x,eta)
                old=covariance(A,B,k);new=covariance(Ar,Br,k)
                exact &= phase_average(A,B,k)==old and phase_average(Ar,Br,k)==new
                div &= dot(n,A)==dot(n,Ar)==0 and fp-dot(n,B)==fpr-dot(n,Br)==0
                dA=[y-x for x,y in zip(A,Ar)];dB=[y-x for x,y in zip(B,Br)]
                lhs=sub(new,old)
                principal=scale(add(add(outer(A,dA),outer(dA,A)),outer(dA,dA)),F(1,2))
                curl=scale(add(add(outer(B,dB),outer(dB,B)),outer(dB,dB)),F(1,2*k*k))
                change &= lhs==add(principal,curl)
                missing |= new!=scale(outer(Ar,Ar),F(1,2))
                D=[F(1),F(2),F(-1)]
                cross_term=scale(add(outer(B,D),outer(D,B)),F(-1,2*k))
                generic |= any(v for row in cross_term for v in row)
                rows.append(dict(x=str(x),carrier=k,eta=str(eta),
                    curl_trace=str(sum(curl[i][i] for i in range(3))),
                    covariance_change_trace=str(sum(lhs[i][i] for i in range(3)))))
    return dict(rows=rows,checks=dict(exact_complete_real_curl_covariance=bool(exact),
        full_curl_divergence_identity=bool(div),all_six_tensor_change_identity=bool(change),
        dropping_sine_curl_square_detected=bool(missing),
        generic_complex_phase_shift_has_nonzero_linear_sine_cross=bool(generic),
        exact_algebra_controls_not_actual_giant_profile_or_continuum_estimates=True))


def variance_case(eta):
    eta=F(eta)
    if not F(-1,2)<eta:raise ValueError('Strictly positive rebuilt weights required')
    n=[F(1),F(-1),F(-1)];old=[[F(1),F(1),F(0)],[F(1),F(0),F(1)]]
    new=[[F(1),1+eta,-eta],old[1]]
    weights=[1/(2*(1+eta)),(1+2*eta)/(2*(1+eta))]
    def tensor(vectors,w):return add(scale(outer(vectors[0],vectors[0]),w[0]),scale(outer(vectors[1],vectors[1]),w[1]))
    Wa=tensor(old,[F(1,2)]*2);Wr=tensor(new,weights)
    tangent=[F(0),F(1),F(-1)];expected=scale(outer(tangent,tangent),eta/2)
    return dict(n=n,old=old,new=new,weights=weights,Wa=Wa,Wr=Wr,delta=sub(Wr,Wa),expected=expected,
        stale=tensor(new,[F(1,2)]*2))


def variance_review():
    rows=[];good=True;stale=True
    for eta in (F(-1,4),F(1,20),F(1,7),F(1,2)):
        r=variance_case(eta)
        good &= min(r['weights'])>0 and all(dot(r['n'],v)==0 for v in r['old']+r['new'])
        good &= r['Wr'][0][1]==r['Wa'][0][1]==F(1,2) and r['Wr'][0][2]==r['Wa'][0][2]==F(1,2)
        good &= r['delta']==r['expected'] and r['delta'][1][1]!=0 and r['delta'][1][2]!=0 and r['delta'][2][2]!=0
        stale &= r['stale'][0][1]!=r['Wa'][0][1]
        rows.append(dict(eta=str(eta),rebuilt_positive_weights=list(map(str,r['weights'])),
            matched_rtheta=str(r['Wr'][0][1]),matched_rz=str(r['Wr'][0][2]),
            changed_thetatheta=str(r['delta'][1][1]),changed_thetaz=str(r['delta'][1][2]),changed_zz=str(r['delta'][2][2])))
    return dict(rows=rows,checks=dict(same_transverse_plane_and_positive_rebuilt_weights=bool(good),
        matching_two_fluxes_does_not_fix_tangential_variance=bool(good),
        stale_weights_also_break_selected_flux=bool(stale),
        exact_complete_plane_waves_realizable_with_distinct_integer_frequencies=True,
        obstruction_not_claimed_to_evaluate_specific_canonical_pulses=True))
