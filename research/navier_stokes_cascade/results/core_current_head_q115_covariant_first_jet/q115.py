"""First parameter jet of the same constrained response; no source certificate."""
from fractions import Fraction as Q

def harmonic(m):
    if type(m) is not int or m==0:raise ValueError('Fixed nonzero integer harmonic required')
    return abs(m)

def residence(m,L,dmin):
    m=harmonic(m);L,dmin=Q(L),Q(dmin)
    if min(L,dmin)<=0:raise ValueError('Original positive length and damping floor required')
    return L if m==1 else min(L,1/((m*m-1)*dmin))

def first_jet(*,m,T,L,dmin,normal,normal_prime,normal_parameter,normal_prime_parameter,inverse_square,
              K,K_parameter,rate_parameter,damping_parameter,source,source_parameter):
    """Euclidean unweighted norm leaves, on the whole identical path.

    Source is normalized g=f/P. Parameters do not differentiate m or endpoints.
    Leaves are premises, not actual native norm enclosures.
    """
    m=harmonic(m)
    vals=[T,L,dmin,normal,normal_prime,normal_parameter,normal_prime_parameter,inverse_square,K,K_parameter,
          rate_parameter,damping_parameter,source,source_parameter]
    if any(isinstance(x,(float,bool)) for x in vals):raise TypeError('Exact bound leaves required')
    T,L,dmin,n,np,na,npa,q,K,Ka,pa,da,G,Ga=map(Q,vals)
    if T<1 or min(L,dmin,q)<=0 or min(n,np,na,npa,K,Ka,pa,da,G,Ga)<0:
        raise ValueError('Positive geometry/normal floor and nonnegative whole-path leaves required')
    ell=residence(m,L,dmin);R=T*ell;U=R*G
    qa=na*q;Qlift=(K*n+np)*q;h=n*q
    forcing=Ga+qa*n*G+(qa*(n*K+np)+Ka+pa+m*m*da+Qlift*na)*U
    W=R*forcing;full=W+h*na*U
    numerator=(n*K+np)*U+n*G
    pressure=(2*q*q*n*na*numerator+q*((n*Ka+K*na+npa)*U
              +(n*K+np)*full+na*G+n*Ga))/m
    return dict(residence=ell,velocity=U,tangent_parameter=W,full_parameter=full,
                parameter_forcing=forcing,normal_constraint=h*na*U,
                m2_residence=m*m*ell,
                uniform_m2_residence=max(L,Q(4,3)/dmin),
                pressure_parameter_times_carrier=pressure,
                actual_source_or_high_jet_contracts_certified=False)

def tangent_forcing(n,np,na,K,Ka,u,g,ga,da,pa,m):
    """Exact 3D matrix expression; accepts symbolic matrices for diagnostics."""
    harmonic(m)
    N=(n.T*n)[0];h=n/N
    # Form the identity without importing a symbolic/numeric backend.
    I=K*0
    for j in range(3):I[j,j]=1
    Pi=I-h*n.T;q=Pi*na/N;Qlift=Pi*(K*n+np)/N
    pressure_numerator=((n.T*K-np.T)*u+n.T*g)[0]
    return (q*pressure_numerator-Pi*Ka*u-(m*m*da+pa)*u-Pi*ga
            +Qlift*(na.T*u)[0])

def record():
    return dict(status='covariant-first-parameter-jet-and-harmonic-residence',
                projected_jet='w=Pi*partial_a u; partial_a u=w-h*(partial_a n dot u)',
                residence='L if |m|=1 else min(L,1/((m^2-1)*dmin))',
                m2_residence_bound='max(L,4/(3*dmin))',
                full_normal_source_and_moving_normal_retained=True,
                all_high_jets_or_actual_source_certified=False,
                full_cone_pulse_residual_or_blowup=False)
