"""Actual phase formulas at the global radius; all background data are fixtures."""
from functools import lru_cache
from fractions import Fraction as F
import math
import mpmath as mp
import numpy as np
from scipy.integrate import quad, quad_vec, solve_ivp
from estimates import (RADIUS, charts, covariance_interval, gaussian_constants,
                       geometry, pulse_barrier, real)


def beta(x):
    return math.exp(1-1/(1-x*x)) if abs(x)<1 else 0.


def transverse_cutoff(xi,r):
    if r<=0:
        raise ValueError('A positive rectangle radius is required')
    return beta(2*xi/r)


def step(x):
    if x<=0:
        return 0.
    if x>=1:
        return 1.
    a,b=math.exp(-1/x),math.exp(-1/(1-x))
    return a/(a+b)


def psi(u):
    return step((.4-abs(u-.5))/.2)


def point(data):
    return [mp.mpf(v)+mp.mpf(a)/data['scale']**3
            for v,a in zip(charts.REP, ('.3','-.2','.1'))]


class Fixture:
    """Freeze exact discrete choices, then evaluate the ODE in binary64.

    The 100-digit phase/frame code checks this independent binary64 frame.
    The manufactured background is explicitly not the actual global field.
    """
    def __init__(self,ell,r,sign):
        self.radius_exact=F(r)
        self.radius=real(r)
        self.data=charts.band(ell,r0=self.radius)
        self.lab=charts.label(self.data,sign)
        self.coeff=charts.local_coefficients(self.data,self.lab,point(self.data))
        self.length=float(self.data['length'])
        self.n0=np.array(list(self.coeff['normal0']),dtype=float)
        self.np=np.array(list(self.coeff['normal_prime']),dtype=float)
        self.kmat=np.array(self.coeff['kmat'].tolist(),dtype=float)
        self.c0=float(self.lab['c0'])
        self.lam=float(self.lab['lambda0'])
        self.u=float(charts.U)
        self.sign=sign
        self.damp_scale=float(self.data['epsilon']*self.data['k']**2)

    def matrices(self,u):
        n=self.n0+(u*self.length)*self.np
        tan,np_t=n[1:],self.np[1:]
        mag=np.linalg.norm(tan)
        magprime=tan@np_t/mag
        kt=tan/mag
        ktprime=np_t/mag-kt*magprime/mag
        nt,ntprime=np.array([kt[1],-kt[0]]),np.array([ktprime[1],-ktprime[0]])
        sa=n[0]/mag
        saprime=self.np[0]/mag-sa*magprime/mag
        root=math.sqrt(1+sa*sa)
        rootprime=sa*saprime/root
        basis=np.empty((3,2))
        basisprime=np.empty((3,2))
        for col,s in enumerate((1,-1)):
            basis[:,col]=np.r_[1.,-sa*kt+s*self.c0*root*nt]
            basisprime[:,col]=np.r_[0.,-saprime*kt-sa*ktprime+s*self.c0*(rootprime*nt+root*ntprime)]
        left=np.array([np.r_[.5,nt/(2*self.c0*root)],np.r_[.5,-nt/(2*self.c0*root)]])
        projected=-self.kmat+np.outer(n,n@self.kmat-self.np)/(n@n)
        growth=self.lam/math.sqrt(1+(self.u*(.5+u))**2)
        dref=self.lam*(1+(self.u*(.5+u))**2)/(1+self.u**2)**1.5
        damping=self.damp_scale*(n@n)
        error=left@(projected@basis-basisprime)-np.diag([growth,-growth])
        return dict(n=n,basis=basis,left=left,basisprime=basisprime,projected=projected,
                    growth=growth,dref=dref,damping=damping,error=error,
                    damping_error=damping-dref,anet=growth-dref)

    def logp(self,u):
        # Exact primitive of the frozen net growth, evaluated without seed underflow.
        u=mp.mpf(float(u))
        y=charts.U*(mp.mpf('.5')+u)
        tilt=charts.U
        return self.lab['lambda0']*self.data['length']/tilt*(mp.asinh(y)-mp.asinh(tilt)
                   -((y+y**3/3)-(tilt+tilt**3/3))/(1+tilt*tilt)**mp.mpf('1.5'))

    def source_shape(self,u):
        return psi(u)*np.array([.3+u,.7-.2*u,-.4+.1*u])

    def solve(self,m=1,source=False,tighter=False,three_component=False,multiplier=None):
        if type(m) is not int or m==0:
            raise ValueError('A nonzero integer harmonic is required')
        if three_component and not source:
            raise ValueError('Independent three-component check uses zero-data sources')
        def forcing(u):
            f=self.source_shape(u) if source else np.zeros(3)
            if multiplier is not None:
                f=f*multiplier(u)
            return f  # Actual source is P*f/L.
        def operator(u):
            a=self.matrices(u)
            if three_component:
                return a['projected']-(m*m*a['damping']+a['anet'])*np.eye(3)
            return (np.diag([0.,-2*a['growth']])+a['error']
                    -(a['damping_error']+(m*m-1)*a['damping'])*np.eye(2))
        def rhs(u,y):
            a=self.matrices(u)
            f=forcing(u)
            projected_f=f-a['n']*(a['n']@f)/(a['n']@a['n'])
            force=projected_f if three_component else a['left']@projected_f
            return self.length*(operator(u)@y)-force
        initial=np.zeros(3 if three_component else 2)
        if not source:
            initial[0]=1.
        sol=solve_ivp(rhs,(0.,1.),initial,method='Radau',jac=lambda u,y:self.length*operator(u),
                      dense_output=True,max_step=.02,rtol=2e-11 if tighter else 2e-9,
                      atol=2e-13 if tighter else 2e-11)
        if not sol.success:
            raise ValueError('Normalized fixture ODE failed: '+sol.message)
        return sol


@lru_cache(maxsize=None)
def frame_review():
    rows=[]
    identity=binary_error=0.
    with mp.workdps(100):
        for ell,r in ((51600,RADIUS),(103200,RADIUS),(51600,RADIUS/2)):
            for sign in (1,-1):
                f=Fixture(ell,r,sign)
                ce=cd=normal=0.
                for u in (0.,.25,.5,.75,1.):
                    a=f.matrices(u)
                    exact=charts.frame(f.data,f.lab,f.coeff,mp.mpf(u)*f.data['length'])
                    binary_error=max(binary_error,float(np.linalg.norm(a['error']-np.array(exact['error'].tolist(),dtype=float))))
                    identity=max(identity,float(np.linalg.norm(a['left']@a['basis']-np.eye(2))),
                                 float(np.linalg.norm(a['n']@a['basis'])))
                    ce=max(ce,f.length*float(np.max(np.abs(a['error']))))
                    cd=max(cd,f.length*abs(a['damping_error']))
                    normal=max(normal,float(f.data['length']*mp.norm(exact['normal']-exact['reference_normal'])))
                rows.append(dict(ell=ell,radius=str(r),sign=sign,length=f.length,
                                 frame_error_times_length=ce,damping_error_times_length=cd,
                                 normal_error_times_length=normal,
                                 nonzero_angular_integer=f.lab['angular_integer']!=0))
    return dict(rows=rows,maximum_binary64_vs_high_precision_frame_error=binary_error,
                maximum_frame_constraint_identity_error=identity,
                checks=dict(globally_admissible_radius_is_used=all(F(r['radius'])<=RADIUS for r in rows),
                            frame_and_phase_have_length_scaled_bounds=all(max(r['frame_error_times_length'],r['damping_error_times_length'],r['normal_error_times_length'])<30 for r in rows),
                            independent_frame_agrees_with_frozen_high_precision=binary_error<2e-14,
                            tangent_constraint_and_left_inverse=identity<2e-14,
                            angular_carriers_remain_nonzero=all(r['nonzero_angular_integer'] for r in rows)),
                diagnostic_entry_bounds_are_sampled_not_interval_certified=True)


def covariance(f,sol,tighter):
    root=math.sqrt(f.length)
    def integrand(w):
        u=.5+w/root
        a=f.matrices(u)
        t=a['basis']@sol.sol(u)
        mass=math.exp(2*float(f.logp(u)))*psi(u)**2
        return mass*np.r_[t[0]**2,t[0]*t[1:]]
    end=root/2
    values,error=quad_vec(integrand,-end,end,points=[x for x in (-8.,-2.,0.,2.,8.) if abs(x)<end],
                          epsabs=2e-12 if tighter else 2e-10,epsrel=2e-11 if tighter else 2e-9)
    transverse,_=quad(lambda x:beta(x)**2,-1.,1.,epsabs=1e-13,epsrel=1e-13)
    determinant=1+(math.sqrt(2)-1)**2
    r=float(f.radius)
    measure=determinant/2*(r/2*transverse)*float(f.data['ci'])*root
    h,column=measure*values[0],measure*values[1:]
    nt,kt=np.array(list(f.lab['nt']),dtype=float),np.array(list(f.lab['kt']),dtype=float)
    ac=-f.c0*math.sqrt(1+f.u*f.u)
    direction=-ac*nt-f.sign*f.u*kt
    ugrid=np.linspace(0.,1.,301)
    framevals=sol.sol(ugrid)
    ratio=framevals[1]/framevals[0]
    logamplitude=np.log(framevals[0])
    xnormalized=np.array([(f.matrices(u)['basis']@sol.sol(u))[0] for u in ugrid])
    gaussian=gaussian_constants(f.lab['lambda0'],f.lab['lambda0'],charts.U)
    barrier=pulse_barrier(f.lab['lambda0']/mp.sqrt(1+(3*charts.U/2)**2),30,30,f.data['length'])
    interval=covariance_interval(f.radius,f.data['length'],gaussian['lower'],gaussian['upper'],
                                 barrier['radial_lower'],barrier['radial_upper'])
    return dict(ell=f.data['ell'],radius=str(f.radius_exact),sign=f.sign,length=f.length,
                h=h,column=column.tolist(),h_times_sqrt_length_over_radius_squared=h*root/(r*r),
                direction_error=float(np.linalg.norm(column/h-direction)),
                direction_error_times_sqrt_length=float(np.linalg.norm(column/h-direction))*root,
                minimum_normalized_radial_pulse=float(np.min(xnormalized)),
                maximum_ratio=float(np.max(np.abs(ratio))),maximum_abs_log_amplitude=float(np.max(np.abs(logamplitude))),
                supplied_ratio_bound=float(barrier['ratio_bound']),
                covariance_lower=float(interval['lower']),covariance_upper=float(interval['upper']),
                full_interval_quadrature_error=float(error),full_pulse_interval_integrated=True)


@lru_cache(maxsize=None)
def covariance_review(tighter=False):
    rows=[]
    minimum_weight=math.inf
    with mp.workdps(100):
        for ell,r in ((51600,RADIUS),(103200,RADIUS),(51600,RADIUS/2)):
            columns=[]
            for sign in (1,-1):
                f=Fixture(ell,r,sign)
                row=covariance(f,f.solve(tighter=tighter),tighter)
                rows.append(row)
                columns.append(row['column'])
            nt=np.array(list(f.lab['nt']),dtype=float)
            kt=np.array(list(f.lab['kt']),dtype=float)
            ac=-f.c0*math.sqrt(1+f.u*f.u)
            # Interior reference target; its weights are diagnosed, not assumed positive.
            weights=np.linalg.solve(np.column_stack(columns),-ac*nt+.2*f.u*kt)
            minimum_weight=min(minimum_weight,float(np.min(weights)))
    return dict(rows=rows,minimum_positive_primary_weight=minimum_weight,
                checks=dict(normalized_h_keeps_transverse_radius_squared=all(.01<r['h_times_sqrt_length_over_radius_squared']<100 for r in rows),
                            positive_radial_pulses_survive=all(r['minimum_normalized_radial_pulse']>.01 for r in rows),
                            sampled_ratio_fits_supplied_length_bound=all(r['maximum_ratio']<r['supplied_ratio_bound'] for r in rows),
                            covariance_fits_supplied_continuum_envelope=all(r['covariance_lower']<r['h']<r['covariance_upper'] for r in rows),
                            column_direction_concentration=all(r['direction_error_times_sqrt_length']<10 for r in rows),
                            primary_weights_positive=minimum_weight>0,
                            full_interval_including_cutoffs=all(r['full_pulse_interval_integrated'] for r in rows)),
                actual_global_profile_or_uniform_entry_constants_evaluated=False)


@lru_cache(maxsize=None)
def sourced_review():
    rows=[]
    with mp.workdps(100):
        f=Fixture(51600,RADIUS,1)
        for m in (1,3):
            framed=f.solve(m=m,source=True,tighter=True)
            direct=f.solve(m=m,source=True,tighter=True,three_component=True)
            error=constraint=residual=omitted=maximum=0.
            for u in np.linspace(0.,1.,101):
                a=f.matrices(u)
                z=framed.sol(u)
                t=a['basis']@z
                forcing=f.source_shape(u)/f.length
                projection=forcing-a['n']*(a['n']@forcing)/(a['n']@a['n'])
                generator=np.diag([0.,-2*a['growth']])+a['error']-(a['damping_error']+(m*m-1)*a['damping'])*np.eye(2)
                tp=a['basisprime']@z+a['basis']@(generator@z-a['left']@projection)+a['anet']*t
                numerator=a['n']@f.kmat@t-f.np@t+a['n']@forcing
                normal_pressure=-a['n']*numerator/(a['n']@a['n'])
                balance=tp+f.kmat@t+m*m*a['damping']*t+forcing+normal_pressure
                error=max(error,float(np.linalg.norm(t-direct.sol(u))))
                constraint=max(constraint,abs(a['n']@t))
                residual=max(residual,float(np.linalg.norm(balance)))
                omitted=max(omitted,f.length*abs(a['n']@forcing)/np.linalg.norm(a['n']))
                maximum=max(maximum,float(np.linalg.norm(t)))
            rows.append(dict(harmonic=m,maximum_two_vs_three_component_error=error,
                             maximum_normal_constraint=constraint,maximum_full_principal_balance_error=residual,
                             omitted_normal_source_pressure_gap_times_length=omitted,
                             maximum_amplitude_over_P=maximum,source_size_over_P_is_O_inverse_length=True))
    return dict(rows=rows,
                checks=dict(independent_projected_and_frame_solutions_agree=all(r['maximum_two_vs_three_component_error']<2e-8 for r in rows),
                            moving_normal_constraint_preserved=all(r['maximum_normal_constraint']<1e-12 for r in rows),
                            full_principal_pressure_and_viscosity_cancel=all(r['maximum_full_principal_balance_error']<1e-12 for r in rows),
                            normal_component_of_source_cannot_be_dropped=all(r['omitted_normal_source_pressure_gap_times_length']>.01 for r in rows),
                            zero_data_source_value_bound_is_consistent=all(r['maximum_amplitude_over_P']<100 for r in rows),
                            higher_harmonic_damping_does_not_shrink_domain=rows[1]['maximum_amplitude_over_P']<rows[0]['maximum_amplitude_over_P']))


@lru_cache(maxsize=None)
def lifted_source_review():
    with mp.workdps(100):
        f=Fixture(51600,RADIUS,1)
        denominator,points=geometry.kernel_points(1)
        matrix=np.array(geometry.J,dtype=float)
        vt=np.array([math.sqrt(2)-1,1.])
        # v=0 is the BOTTOM of the band rectangle, not its center.
        base=np.linalg.solve(matrix,np.array([0.,1/geometry.PRIME])-float(f.radius)*vt)
        shifts=np.array(points,dtype=float)/denominator
        starts=base+shifts
        def multiplier(copy,u):
            # Actual common path: J^{-1}(ci*L*u*vt) = 2*r*u/Tg*vt.
            # A band transverse bump with value one at xi=0, plus fixed slow
            # bumps, makes this a supported source without changing this path.
            h=starts[copy]+2*float(f.radius)*u/(4+math.sqrt(2))*vt
            return math.cos(2*math.pi*h[0])
        values=np.array([[multiplier(j,u) for j in range(denominator)] for u in (0.,.25,.5,1.)])
        path_error=0.
        for u in (0.,.25,.5,1.):
            desired=np.array([0.,1/geometry.PRIME])+(2*u-1)*float(f.radius)*vt
            for start in starts:
                actual=matrix@(start+2*float(f.radius)*u/(4+math.sqrt(2))*vt)-desired
                path_error=max(path_error,float(np.linalg.norm(actual-np.rint(actual))))
        first=0
        second=int(np.argmin(values[0]))
        a=f.solve(source=True,tighter=True,multiplier=lambda u:multiplier(first,u))
        b=f.solve(source=True,tighter=True,multiplier=lambda u:multiplier(second,u))
        difference=float(np.linalg.norm(a.sol(.5)-b.sol(.5)))
        # The nontrivial character averages to zero over this covering fiber.
        mean=float(np.max(np.abs(np.mean(values,axis=1))))
        gap=float(np.max(values[0])-np.min(values[0]))
    return dict(copies=denominator,maximum_fiber_average=mean,source_value_gap=gap,
                midpoint_inverse_response_gap=difference,maximum_actual_band_path_error=path_error,
                checks=dict(all_fourteen_common_preimages_evaluated=denominator==14,
                            every_path_starts_at_the_actual_rectangle_bottom=path_error<2e-15,
                            allowed_source_differs_across_band_preimages=gap>1,
                            averaging_copies_erases_this_source=mean<2e-15,
                            liftwise_inverse_does_not_identify_preimages=difference>.01))


@lru_cache(maxsize=None)
def harmonic_factor_review():
    with mp.workdps(100):
        f=Fixture(51600,RADIUS,1)
        def propagator(m):
            def rhs(v,y):
                a=f.matrices(.5+v/f.length)
                operator=np.diag([a['growth'],-a['growth']])+a['error']-m*m*a['damping']*np.eye(2)
                return (operator@y.reshape(2,2)).reshape(4)
            return solve_ivp(rhs,(0.,.7),np.eye(2).reshape(4),method='DOP853',rtol=2e-12,atol=2e-14).y[:,-1].reshape(2,2)
        integral,_=quad(lambda v:f.matrices(.5+v/f.length)['damping'],0.,.7,epsabs=1e-13)
        one=propagator(1)
        factor=math.exp(-8*integral)
        higher=propagator(3)
        error=float(np.linalg.norm(higher-factor*one))
        negative=float(np.linalg.norm(propagator(-3)-higher))
    return dict(scalar_damping_factor=factor,maximum_factor_identity_error=error,
                opposite_harmonic_error=negative,
                checks=dict(exact_scalar_damping_identity=error<2e-12,
                            higher_harmonics_only_add_damping=0<factor<1,
                            opposite_harmonics_share_the_same_propagator=negative<2e-14))


def failed_shortcuts_review():
    with mp.workdps(100):
        f=Fixture(51600,RADIUS,1)
        endpoint=f.logp(0.)
    transverse,_=quad(lambda x:beta(x)**2,-1.,1.,epsabs=1e-13)
    return dict(old_fixture_radius=F(1,25),global_radius=RADIUS,
                omitted_transverse_factor=float(RADIUS)/2*transverse,
                positive_endpoint_log_P=float(endpoint),
                ideal_inverse_growth_when_radius_divided_by_four=8,
                checks=dict(old_fixture_radius_is_not_global_geometry=F(1,25)>RADIUS,
                            dropping_transverse_area_changes_normalization=float(RADIUS)/2*transverse<1e-6,
                            radius_uniform_inverse_is_false=4**1.5==8,
                            positive_finite_endpoint_is_not_zero_extension=mp.isfinite(endpoint) and endpoint<0,
                            psi_is_zero_near_both_endpoints=psi(0.)==psi(.05)==psi(.95)==psi(1.)==0.,
                            psi_keeps_the_central_fifth=all(abs(psi(u)-1)<2e-15 for u in (.3,.5,.7))))
