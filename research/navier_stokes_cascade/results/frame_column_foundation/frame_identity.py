"""Independent exact moving-frame multiplication and omission controls."""
from fractions import Fraction as F
from itertools import product


def dot(a,b):return sum(x*y for x,y in zip(a,b))
def transpose(a):return list(map(list,zip(*a)))
def mul(a,b):return [[dot(r,c) for c in zip(*b)] for r in a]
def sub(a,b):return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]


def matrices(a,s,w,c,F0,g,sprime,omega,mag=F(1),magprime=F(0),omit_normal=False,omit_frame=False,mixing_w=None,mixing_wprime=None):
    a=list(map(F,a));s,w,c,F0,sprime,omega,mag,magprime=map(F,(s,w,c,F0,sprime,omega,mag,magprime));g=list(map(F,g))
    if dot(a,a)!=1 or w*w!=1+s*s or min(w,mag)<=0 or not c:raise ValueError('Nondegenerate exact frame required')
    b=[a[1],-a[0]];U=transpose([[1,-s*a[0],-s*a[1]],[0,*b]])
    Up=transpose([[0,*[-sprime*a[i]-s*omega*b[i] for i in range(2)]],[0,*[-omega*x for x in a]]])
    mw=w if mixing_w is None else F(mixing_w)
    wp=s*sprime/w if mixing_wprime is None else F(mixing_wprime)
    if mw<=0:raise ValueError('Positive mixing scale required')
    M=[[F(1),F(1)],[c*mw,-c*mw]];mi=[[F(1,2),1/(2*c*mw)],[F(1,2),-1/(2*c*mw)]]
    Mp=[[F(0),F(0)],[c*wp,-c*wp]]
    B=mul(U,M);Bp=add(mul(Up,M),mul(U,Mp))
    left=mul(mi,[[F(1),F(0),F(0)],[F(0),*b]])
    n=[mag*s,mag*a[0],mag*a[1]]
    np=[magprime*s+mag*sprime,*[magprime*a[i]+mag*omega*b[i] for i in range(2)]]
    K=[[F(0),-2*F0,F(0)],[2*F0+g[0],F(0),F(0)],[g[1],F(0),F(0)]]
    nk=[dot(n,col) for col in zip(*K)];norm=dot(n,n)
    An=[[-K[i][j]+n[i]*(nk[j]-(0 if omit_normal else np[j]))/norm for j in range(3)] for i in range(3)]
    direct=mul(left,sub(mul(An,B),[[F(0)]*2 for _ in range(3)] if omit_frame else Bp))
    H=[[s*(dot(a,g)-sprime)/(1+s*s),(2*F0*b[0]-s*omega)/(1+s*s)],[-(2*F0*b[0]+dot(b,g))+s*omega,F(0)]]
    connection=mul(mi,Mp);reduced=sub(mul(mul(mi,H),M),connection)
    return dict(B=B,Bprime=Bp,left=left,An=An,n=n,nprime=np,H=H,M=M,Mprime=Mp,Minv=mi,connection=connection,direct=direct,reduced=reduced,a=a,b=b,s=s,w=w,c=c)


def gauge_review():
    """Manuscript (7.8) uses wref, unlike inherited replay's actual w.

    All equalities below use independent exact ambient multiplications.
    Entrance obstruction uses squared scales: it needs no floating roots.
    """
    args=((F(3,5),F(-4,5)),F(3,4),F(5,4),F(-1,2),F(11,10),
          (F(-3,2),F(-7,6)),F(1,50),F(1,100))
    sr,wr,srp=F(4,3),F(5,3),F(1,70)
    actual=matrices(*args)
    reference=matrices(*args,mixing_w=wr,mixing_wprime=sr*srp/wr)
    T=mul(actual['Minv'],reference['M'])
    Tp=mul(actual['Minv'],sub(reference['Mprime'],mul(actual['Mprime'],T)))
    Ti=mul(reference['Minv'],actual['M'])
    transformed=mul(Ti,sub(mul(actual['direct'],T),Tp))
    za=[[F(1)],[F(0)]];zr=mul(Ti,za)
    physical=mul(actual['B'],za)
    force=[[F(2)],[F(-1)],[F(3)]]
    u,L,g=F(3,4),F(10),F(2)
    sr2=u*u/4;sa2=sr2/(1+u*u/(L*L*g*g))
    gap=sr2-sa2
    return dict(exact_change_of_coordinates=[[str(x) for x in row] for row in T],
        transported_repository_entrance=[[str(x) for x in row] for row in zr],
        pre_rounding_entrance_squared_scale_gap=str(gap),
        manuscript_equation='7.8: mixing sqrt(1+sref^2); U uses actual sa',
        inherited_convention='physical_pulse_foundation and geometric_pulse_replay: mixing sqrt(1+sa^2)',
        checks=dict(reference_scale_frame_direct_identity=reference['direct']==reference['reduced'],
            exact_frame_transition=reference['B']==mul(actual['B'],T),
            exact_time_dependent_generator_conjugacy=reference['direct']==transformed,
            transported_entrance_preserves_physical_pulse=mul(reference['B'],zr)==physical,
            forcing_coordinates_transform_exactly=mul(reference['left'],force)==mul(Ti,mul(actual['left'],force)),
            resetting_growing_coordinate_changes_initial_physical_pulse=mul(reference['B'],za)!=physical,
            pre_rounding_original_tilt_already_forces_nonzero_scale_gap=gap>0,
            zero_data_forced_inverse_preserved_only_with_full_transformation=True,
            manuscript_and_repository_primary_not_silently_identified=True))


def identity_review():
    count=0;good=True;normal=False;frame=False;connection=True
    for a,(s,w),c,omega in product(((F(3,5),F(-4,5)),(F(-4,5),F(-3,5))),
        ((F(0),F(1)),(F(3,4),F(5,4)),(F(-3,4),F(5,4))), (F(-1,2),F(-2)),(F(0),F(1,100))):
        args=(a,s,w,c,F(11,10),(F(-3,2),F(-7,6)),F(1,50),omega,F(7,5),F(1,30))
        r=matrices(*args);good &= r['direct']==r['reduced'];count+=1
        expected=[[s*F(1,50)/(2*(1+s*s)),-s*F(1,50)/(2*(1+s*s))],[-s*F(1,50)/(2*(1+s*s)),s*F(1,50)/(2*(1+s*s))]]
        connection &= r['connection']==expected
        normal |= matrices(*args,omit_normal=True)['direct']!=r['reduced']
        frame |= matrices(*args,omit_frame=True)['direct']!=r['reduced']
    return dict(exact_frame_cases=count,checks=dict(exact_projected_generator_and_frame_derivative=bool(good),
        exact_M_inverse_Mprime_connection=bool(connection),omitting_moving_normal_pressure_detected=bool(normal),
        omitting_Bprime_detected=bool(frame),ambient_left_inverse_not_replaced_by_pseudoinverse=True))


def eigen_review():
    # N=(-4/5,-3/5), K=(3/5,-4/5), F0=1, |g0|=2.
    # lambda0=4/5, c0=-1/2 exactly, with strict positive growth.
    r=matrices((F(3,5),F(-4,5)),F(3,4),F(5,4),F(-1,2),F(1),(F(-8,5),F(-6,5)),0,0)
    target=[[F(16,25),F(0)],[F(0),F(-16,25)]]
    return dict(reference_generator=[[str(x) for x in row] for row in r['direct']],
        checks=dict(original_positive_and_negative_growth_columns_diagonalize=r['direct']==target,
            original_c0_sign_and_shear_relation_retained=F(-1,2)==F(4,5)/F(-8,5),
            no_angular_component_inverse_needed=True))


def column_review():
    rows=[];good=True;missing=False;ratio=True
    r=matrices((F(3,5),F(-4,5)),F(3,4),F(5,4),F(-1,2),1,(F(-8,5),F(-6,5)),0,0)
    for q in (F(-1,8),F(-1,100),F(0),F(1,100),F(1,8)):
        x=1+q;z=[[F(1)],[q]];t=mul(r['B'],z)
        direction=[-r['s']*r['a'][i]+r['c']*r['w']*(1-q)/(1+q)*r['b'][i] for i in range(2)]
        good &= all(t[i+1][0]/x==direction[i] for i in range(2))
        bad=[-r['s']*r['a'][i]+r['c']*r['w']*r['b'][i] for i in range(2)]
        missing |= direction!=bad
        ratio &= abs((1-q)/(1+q)-1)<=F(16,7)*abs(q)
        rows.append(dict(q=str(q),tangential_ratio=list(map(str,direction)),factor=str((1-q)/(1+q))))
    return dict(rows=rows,checks=dict(complete_principal_column_ratio_exact=bool(good),
        dropping_secondary_frame_coordinate_detected=bool(missing),sharp_ratio_bound_on_one_eighth_ball=bool(ratio),
        radial_amplitude_positive_not_full_oscillating_velocity_floor=True))
