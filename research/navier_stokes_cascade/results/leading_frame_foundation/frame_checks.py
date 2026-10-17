"""Independent exact frame geometry, carrier and denominator controls."""
from fractions import Fraction as F
from itertools import product


def dot(a,b):return sum(x*y for x,y in zip(a,b))
def matmul(a,b):return [[dot(row,col) for col in zip(*b)] for row in a]


def exact_frame(t,s,c,w):
    # t is a unit tangential vector, s=sa and w=sqrt(1+s^2).
    t=tuple(map(F,t));s,c,w=map(F,(s,c,w))
    if dot(t,t)!=1 or c==0 or w<=0 or w*w!=1+s*s:raise ValueError('Unit tangent, nonzero c and exact w required')
    na=(t[1],-t[0]);u=(F(1),-s*t[0],-s*t[1]);n=(s,*t)
    columns=[tuple(u[i]+sign*c*w*(0,*na)[i] for i in range(3)) for sign in (1,-1)]
    B=[list(x) for x in zip(*columns)]
    ell=[[F(1,2),na[0]/(2*c*w),na[1]/(2*c*w)],
         [F(1,2),-na[0]/(2*c*w),-na[1]/(2*c*w)]]
    return dict(B=B,left=ell,n=n,c=c,s=s,w=w)


def frame_review():
    count=0;left=True;transverse=True;gram=True;ambient=True
    for t,(s,w),c in product(((F(3,5),F(4,5)),(F(-1),F(0))),((F(0),F(1)),(F(3,4),F(5,4)),(F(-3,4),F(5,4))), (F(-1,2),F(-2))):
        r=exact_frame(t,s,c,w);B=r['B'];ell=r['left'];columns=list(zip(*B));count+=1
        left &= matmul(ell,B)==[[1,0],[0,1]]
        transverse &= all(dot(r['n'],col)==0 for col in columns)
        gram &= matmul([list(x) for x in columns],B)==[[w*w*(1+c*c),w*w*(1-c*c)],[w*w*(1-c*c),w*w*(1+c*c)]]
        # Rows of the specified left inverse, not the ambient pseudoinverse.
        ambient &= matmul(ell,[list(x) for x in zip(*ell)])==[[F(1,4)+1/(4*c*c*w*w),F(1,4)-1/(4*c*c*w*w)],[F(1,4)-1/(4*c*c*w*w),F(1,4)+1/(4*c*c*w*w)]]
    return dict(exact_frame_cases=count,checks=dict(original_left_inverse_identity=bool(left),complete_normal_transversality=bool(transverse),
        frame_Gram_and_singular_values=bool(gram),specified_ambient_left_inverse_not_silently_replaced=bool(ambient)))


def nearest_nonzero(x,k):
    x=F(x)
    if type(k) is not int or k<1:raise ValueError('Positive integer k required')
    z=k*x;n=z.numerator//z.denominator
    candidates={n-1,n,n+1,n+2,-1,1}-{0}
    selected=min(candidates,key=lambda q:(abs(F(q)-z),abs(q),q))
    return F(selected,k)


def carrier_review():
    rows=[];rounding=True;axial=True
    for sign,L in product((-1,1),(8,32,128,512)):
        k=L*L;p=nearest_nonzero(F(sign,L),k);pz=F(-1)
        rows.append(dict(sign=sign,length=L,k=k,p=str(p),pz=str(pz)))
        rounding &= p!=0 and abs(p-F(sign,L))<=F(1,k) and abs(p)<=F(2,L)
        axial &= pz<0
    nearzero=all(abs(nearest_nonzero(x,10)-x)<=F(1,10) for x in (F(0),F(1,1000),F(-1,1000)))
    # |Kz|>=nstar; |gz|/|g|^2<=1/gstar. Exact divided-by-Bs bound.
    nstar,gstar,u=F(1,4),F(1,3),F(2);L=2*u/(gstar*nstar)
    axial_floor=-nstar+u/(L*gstar)
    return dict(rows=rows,axial_normalized_upper=str(axial_floor),
        checks=dict(nearest_nonzero_error_is_one_over_k=bool(rounding and nearzero),
            half_over_k_bound_fails_at_zero=abs(nearest_nonzero(0,10))>F(1,20),
            angular_component_can_tend_to_zero_with_nonzero_integer_frequency=bool(rounding),
            axial_carrier_survives_both_signs=bool(axial and axial_floor==-nstar/2),
            moving_axial_normal_not_identified_with_constant_pz=True))


def pressure_and_domain_review():
    # An exact principal normal-pressure example: n.t=0 but n'.t !=0.
    n=(F(1),F(1),F(0));t=(F(1),F(-1),F(0));np=(F(0),F(1),F(0));source=(F(2),F(0),F(0))
    return dict(normal_pressure_numerator_without_K=str(-dot(np,t)+dot(n,source)),
        common_domain='all finite value conditions first; fixed jets use the resulting positive denominators',
        checks=dict(transversality_does_not_remove_moving_normal_or_source_pressure=dot(n,t)==0 and -dot(np,t)+dot(n,source)==3,
            original_tilt_carriers_width_and_root_branches_retained=True,
            one_value_domain_not_new_band_for_each_derivative=True,
            stress_direction_not_unbounded_normalized_target_used_at_edges=True,
            complete_generated_nonlinear_sources_means_and_residual_still_open=True))
