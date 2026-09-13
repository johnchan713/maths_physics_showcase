"""Restore U, correct five moments, and evaluate the exact source identities.

Incoming moment differences are data. This module does not manufacture a
regular axis behind them. Its fixtures are labelled and its acceptance test
requires supplied uniform C1 bounds, not a maximum over a finite grid.
"""
import mpmath as mp
from moments import Jet,step,step_prime


def accept_entry_bounds(moment_C1,g_C1):
    if len(moment_C1)!=5:
        raise ValueError('All five uniform moment bounds are required')
    values=list(map(mp.mpf,moment_C1))+[mp.mpf(g_C1)]
    return all(mp.isfinite(v) and 0<=v<=mp.mpf('1e-16') for v in values)


def rows_to_physical(rows,eta,P,XR):
    """Invert the normalized row combinations, including their eta derivatives."""
    if len(rows)!=5 or not all(isinstance(v,Jet) for v in rows):
        raise ValueError('Five value-and-derivative jets required')
    e=Jet(eta,1);f=1/(1+e**2);P=mp.mpf(P);XR=mp.mpf(XR)
    xc=mp.exp(-6)
    M=XR*xc*rows[0]
    I=XR**mp.mpf('1.5')*mp.sqrt(2)*P*f*xc**mp.mpf('1.6')*rows[2]
    J=4*e*I+XR**mp.mpf('1.5')*mp.sqrt(2)*P*f*xc**mp.mpf('1.6')*rows[1]
    S=8*e*M-XR*P**2*f**2*xc**mp.mpf('1.2')*rows[3]
    Cp=P**2*f**2*xc**mp.mpf('.2')*rows[4]
    return [M,I,J,S,Cp]


def physical_to_rows(moments,eta,P,XR):
    """Accept actual physical moment defects without dropping angular derivatives."""
    if len(moments)!=5 or not all(isinstance(v,Jet) for v in moments):
        raise ValueError('Five value-and-derivative jets required')
    e=Jet(eta,1);f=1/(1+e**2);P=mp.mpf(P);XR=mp.mpf(XR)
    if P<=0 or XR<=0:
        raise ValueError('Positive scales required')
    xc=mp.exp(-6);M,I,J,S,Cp=moments
    angular=XR**mp.mpf('1.5')*mp.sqrt(2)*P*f*xc**mp.mpf('1.6')
    return [M/(XR*xc),(J-4*e*I)/angular,I/angular,
            -(S-8*e*M)/(XR*P**2*f**2*xc**mp.mpf('1.2')),
            Cp/(P**2*f**2*xc**mp.mpf('.2'))]


def pressure_diagnostic(eta):
    """A diagnostic shape inside the proved envelope; not the full outer datum."""
    e=Jet(eta,1)
    return -mp.mpf('3.31462273001433255')/(1+e**2)**2


def ideal(eta,h,P,x,pressure=pressure_diagnostic):
    """Closed ideal Q and N/P^2, independent of the moment evaluator below."""
    eta,h,P,x=map(mp.mpf,(eta,h,P,x))
    f=1/(1+eta**2);d=1-eta**2;A=mp.mpf('.5')+h;D=mp.mpf('.5')-h
    pi=pressure(eta)
    Q=mp.mpf(9)/8-5*h/8+2*h*eta**2 \
        +mp.mpf(5)/4*eta**2*(D+4*d)/(1+eta**2)
    radial=eta*((5+mp.mpf(25)*h/3)*f*f+mp.mpf(25)/3*d*f**3)*x**mp.mpf('.2')
    N=(-20*eta+32*(1+h)*eta**3)/P**2+4*A*eta*pi.v-d*pi.d+radial
    return dict(Q=Q,N_over_P2=N,N_radial_part=radial,
                source_Q=mp.mpf('1.6')*Q,source_N_over_P2=N+radial/5)


class Annulus:
    def __init__(self,matrix,h,P,XR=10000,pressure=pressure_diagnostic):
        self.matrix=matrix
        self.h,self.P,self.XR=map(mp.mpf,(h,P,XR))
        if not 0<self.h<=mp.mpf('.01') or self.P<16 or self.XR<=0:
            raise ValueError('Outside the annulus parameter domain')
        self.pressure=pressure

    def fields(self,y,eta,data,solution):
        y,eta=mp.mpf(y),mp.mpf(eta);t=y+6
        if not -8<=y<=-5 or abs(eta)>1:
            raise ValueError('Outside the declared annulus')
        e=Jet(eta,1);f=1/(1+e**2)
        r=step(-t-1)
        v=data['g']*r
        vt=-data['g'].v*step_prime(t+2)
        epsilon=Jet(0);epsilon_t=mp.mpf(0)
        for j,(c,ce) in enumerate(zip(solution['root'],solution['root_eta'])):
            beta=self.matrix.bump(j,t);betat=self.matrix.bump(j,t,True)
            if j<2:
                v+=Jet(c,ce)*beta;vt+=c*betat
            else:
                epsilon+=Jet(c,ce)*beta;epsilon_t+=c*betat
        U=4*e+v
        Ebar=f*mp.exp(y/10)*(1+epsilon)
        return dict(U=U,Ebar=Ebar,Ut=vt,epsilon=epsilon,epsilon_t=epsilon_t,
                    a=mp.mpf('.8')-2*epsilon_t/(1+epsilon.v),
                    ell=mp.mpf('.6')+epsilon_t/(1+epsilon.v))

    def state(self,y,eta,data,solution,frozen_derivative=False):
        y,eta=mp.mpf(y),mp.mpf(eta);t=y+6;x=mp.exp(y);xc=mp.exp(-6)
        e=Jet(eta,1);f=1/(1+e**2)
        h,P=self.h,self.P;A=mp.mpf('.5')+h;D=mp.mpf('.5')-h
        L=1-2*h*eta**2;d=1-eta**2
        field=self.fields(y,eta,data,solution)
        r=self.matrix.ledger(data,solution,eta,P,t,frozen_derivative)
        # The original five primitives, with only common P factors removed.
        M=4*e*x+xc*r[0]
        I=f*(x**mp.mpf('1.6')/mp.mpf('1.6')+xc**mp.mpf('1.6')*r[2])
        J=4*e*I+f*xc**mp.mpf('1.6')*r[1]
        S=16*e**2*x/P**2-mp.mpf(5)/12*f**2*x**mp.mpf('1.2') \
            +8*e*xc*r[0]/P**2-f**2*xc**mp.mpf('1.2')*r[3]
        Cp=mp.mpf('2.5')*f**2*x**mp.mpf('.2')+f**2*xc**mp.mpf('.2')*r[4]
        Pi=self.pressure(eta)+Cp
        W=1-(2*D*eta*M.v+d*M.d)/x
        U,Ebar=field['U'],field['Ebar']
        Q=-W+((1-h)*I.v-D*eta*I.d-d*J.d+2*(h-D)*eta*J.v)/(x**mp.mpf('1.5')*Ebar.v)
        N=-W*U.v/P**2+D*(M.v-eta*M.d)/(P**2*x) \
            +(4*h*eta*S.v-d*S.d)/x+4*A*eta*Pi.v-d*Pi.d
        a=field['a'];bs_scaled=2*field['Ut']/Ebar.v
        G=Q-bs_scaled*N/(a*Ebar.v)
        vs=a+bs_scaled**2/(P**2*a)
        Pc=self.XR*x*G/L
        Hc=D*eta+d*U.v
        Sq=-W*field['ell']-h*(1-2*eta*U.v)-Hc*Ebar.d/Ebar.v
        Sn=-W*field['Ut']/P**2-A*(1-2*eta*U.v)*U.v/P**2 \
            -Hc*U.d/P**2-d*Pi.d+4*A*eta*Pi.v+eta*Ebar.v**2
        base=ideal(eta,h,P,x,self.pressure)
        return dict(y=y,eta=eta,U=U.v,U_eta=U.d,E_over_P=Ebar.v,E_eta_over_P=Ebar.d,
            a=a,P_times_bs=bs_scaled,Q=Q,N_over_P2=N,G=G,vs=vs,Pc=Pc,
            Q_defect=Q-base['Q'],N_over_P2_defect=N-base['N_over_P2'],
            pressure_defect_over_P2=f.v**2*xc**mp.mpf('.2')*r[4].v,
            ledger=[v.v for v in r],ledger_eta=[v.d for v in r],
            source_Q_defect=Sq-base['source_Q'],source_N_over_P2_defect=Sn-base['source_N_over_P2'],
            Q_defect_radial_rhs=Sq-base['source_Q']-(field['ell']+1)*(Q-base['Q'])
                -(field['ell']-mp.mpf('.6'))*base['Q'],
            N_defect_radial_rhs=Sn-base['source_N_over_P2']-(N-base['N_over_P2']),
            strict_relaxed_cone=bool(vs<1 and Pc>2),full_admissible_cone_claimed=False)

    def source_check(self,y,eta,data,solution):
        """Differentiate moment primitives and compare to the original sources.

        Work with edited-minus-ideal quantities so that the tiny correction,
        rather than the order-one background, determines the error scale.
        """
        state=self.state(y,eta,data,solution)
        rows=[]
        for name,rhs in (('Q_defect','Q_defect_radial_rhs'),
                         ('N_over_P2_defect','N_defect_radial_rhs')):
            derivative=mp.diff(lambda yy:self.state(yy,eta,data,solution)[name],mp.mpf(y))
            scale=max(abs(derivative),abs(state[rhs]),data['epsilon'])
            rows.append(abs(derivative-state[rhs])/scale)
        return dict(y=y,eta=eta,normalized_source_gaps=rows)


def radius_requirements(log_C,log_P,Tsh):
    """Geometry only: this does not certify the incoming core moments."""
    log_C,log_P,Tsh=map(mp.mpf,(log_C,log_P,Tsh))
    if log_C<0 or log_P<mp.log(16) or Tsh<0:
        raise ValueError('Invalid amplitude, pressure scale, or transition length')
    log_XR=mp.log(110)+10*(log_C+log_P)
    log_xsep=Tsh-10*(log_C+log_P)
    return dict(log_XR=log_XR,log_xsep=log_xsep,
                transition_before_restoration=bool(log_xsep<-8),
                radius_large_enough=bool(log_XR>=mp.log(10000)),
                core_moments_verified=False)
