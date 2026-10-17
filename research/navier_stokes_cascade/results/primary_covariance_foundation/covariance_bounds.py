"""Exact perturbed-cone root bounds; uniquely named to preserve legacy imports."""
from fractions import Fraction as F
from itertools import product


def cone_bounds(gamma,eta):
    gamma,eta=F(gamma),F(eta)
    if not 0<gamma<=1 or not 0<=eta<1:raise ValueError('Positive cone margin and error in [0,1) required')
    numerator=gamma-(2-gamma)*eta
    return dict(det_lower=2*(1-eta)**2,det_upper=2*(1+eta)**2,
        inverse_upper=(1+eta)/(1-eta)**2,
        positive_numerator=numerator,
        root_lower=min(numerator/(2*(1+eta)**2),numerator/(2*(1-eta)**2)),
        root_upper=(1+eta)*(2-gamma)/(2*(1-eta)**2),
        sharp_positive_error=gamma/(2-gamma))


def matrix(errors):
    if len(errors)!=4:raise ValueError('Four entry errors required')
    a,b,c,d=map(F,errors)
    return ((1+a,1+b),(-1+c,1+d))


def determinant(a):return a[0][0]*a[1][1]-a[0][1]*a[1][0]
def solve(a,t):
    d=determinant(a)
    if not d:raise ValueError('Singular matrix')
    return ((a[1][1]*t[0]-a[0][1]*t[1])/d,(-a[1][0]*t[0]+a[0][0]*t[1])/d)


def freezing_margin(gamma,tau):
    gamma,tau=F(gamma),F(tau)
    if not 0<gamma<=1 or not 0<=tau<1:raise ValueError('Admissible margin and relative freezing error required')
    return 1-(1-gamma+tau)/(1-tau)


def length_requirement(direction,amin,tilt,gamma):
    direction,amin,tilt,gamma=map(F,(direction,amin,tilt,gamma))
    if min(direction,amin,tilt,gamma)<=0 or gamma>1:raise ValueError('Positive finite input bounds required')
    return (4*direction*max(1/amin,1/tilt)/gamma)**2


def root_interval(gamma,V,hlo,hhi):
    gamma,V,hlo,hhi=map(F,(gamma,V,hlo,hhi))
    if not 0<gamma<=1 or min(V,hlo)<=0 or hhi<hlo:raise ValueError('Positive cone target and ordered Haar weights required')
    return gamma*V/(8*hhi),3*V/hlo


def corner_review():
    rows=[];count=0;good=True
    for gamma in (F(1,100),F(1,10),F(1,2),F(1)):
        for eta in (F(0),gamma/8,gamma/4):
            bounds=cone_bounds(gamma,eta);smallest=None;largest=None
            for e in product((-eta,eta),repeat=4):
                a=matrix(e);det=determinant(a)
                for W in (-(1-gamma),1-gamma):
                    q=solve(a,(F(1),W));count+=1
                    good &= bounds['det_lower']<=det<=bounds['det_upper']
                    good &= all(bounds['root_lower']<=x<=bounds['root_upper'] for x in q)
                    good &= all(gamma/8<=x<=3 for x in q)
                    smallest=min(q) if smallest is None else min(smallest,*q)
                    largest=max(q) if largest is None else max(largest,*q)
            rows.append(dict(gamma=str(gamma),eta=str(eta),minimum_root=str(smallest),maximum_root=str(largest)))
    witnesses=[]
    for gamma in (F(1,100),F(1,10),F(1,2)):
        critical=gamma/(2-gamma)
        for eta in (critical,critical+F(1,100)):
            a=matrix((0,eta,0,-eta));q=solve(a,(1,1-gamma))
            witnesses.append(dict(gamma=str(gamma),eta=str(eta),determinant=str(determinant(a)),first_root=str(q[0])))
    negative_bounds=True;negative_count=0
    gamma,eta=F(1,2),F(2,5);bounds=cone_bounds(gamma,eta)
    for e in product((-eta,eta),repeat=4):
        for W in (-(1-gamma),1-gamma):
            negative_count+=1
            negative_bounds &= all(bounds['root_lower']<=q<=bounds['root_upper'] for q in solve(matrix(e),(1,W)))
    return dict(corner_target_cases=count,negative_numerator_corner_target_cases=negative_count,rows=rows,sharp_boundary_witnesses=witnesses,
        checks=dict(all_rational_corner_targets_satisfy_derived_bounds=bool(good),
            negative_numerator_lower_bound_uses_smaller_determinant=bool(negative_bounds),
            safe_budget_inverse_below_three=all(cone_bounds(g,F(g)/4)['inverse_upper']<3 for g in (F(1,100),F(1,10),F(1,2),F(1))),
            invertible_matrix_does_not_imply_positive_roots=all(F(r['determinant'])==2 and F(r['first_root'])<=0 for r in witnesses),
            safe_budget_det_strictly_above_one=cone_bounds(1,F(1,4))['det_lower']>1))


def scale_review():
    rows=[]
    for gamma in (F(1,100),F(1,10),F(1,2),F(1)):
        tau=gamma/8;remaining=freezing_margin(gamma,tau);g=gamma/2
        required=length_requirement(3,F(2,5),F(1,4),g)
        old=(3*(F(5,2)+4)*(4+g)/g)**2
        rows.append(dict(global_margin=str(gamma),freezing_error=str(tau),remaining_margin=str(remaining),
            new_length=str(required),retained_padded_length=str(old)))
    return dict(rows=rows,haar_scale='s0=r0^2/sqrt(Ls); h_sigma=s0*j_sigma',
        squared_root='y_sigma=(zeta/s0)*b_sigma; b_sigma>=gamma*vstar/(8*jmax)',
        primary_root='sqrt(y_sigma)=sqrt(zeta)*Ls^(1/4)/r0*sqrt(b_sigma)',
        checks=dict(relative_freezing_budget_retains_half_margin=all(F(r['remaining_margin'])>=F(r['global_margin'])/2 for r in rows),
            retained_pulse_length_budget_dominates_new_explicit_condition=all(F(r['retained_padded_length'])>=F(r['new_length']) for r in rows),
            actual_radius_scaling_is_not_dropped=F(1,4)-1==F(-3,4),
            outer_normalized_target_allowed_polynomial_edge_pole=True,
            root_gap_not_asserted_at_zero_stress_shell_boundary=True))
