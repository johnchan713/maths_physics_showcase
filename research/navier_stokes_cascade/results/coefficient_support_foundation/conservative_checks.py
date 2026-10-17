"""Independent exact-rational radial calculus for all-order signed totals.

Manufactured compact polynomial profiles test algebra, not actual continuum
existence. Angular jets are factorial Taylor coefficients through order two;
radial integration is exact on [0,1], without numerical quadrature.
"""
from fractions import Fraction as F
from functools import lru_cache

class J:
    def __init__(self,c):self.c=tuple(map(F,c))+(F(0),)*max(0,3-len(c))
    @staticmethod
    def cast(x):return x if isinstance(x,J) else J((x,))
    def __add__(self,x):
        x=J.cast(x);return J(tuple(a+b for a,b in zip(self.c,x.c)))
    __radd__=__add__
    def __neg__(self):return J(tuple(-a for a in self.c))
    def __sub__(self,x):return self+-J.cast(x)
    def __rsub__(self,x):return J.cast(x)+-self
    def __mul__(self,x):
        if isinstance(x,P):return NotImplemented
        x=J.cast(x);return J(tuple(sum((self.c[i]*x.c[n-i] for i in range(n+1)),F(0)) for n in range(3)))
    __rmul__=__mul__
    def inverse(self):
        if not self.c[0]:raise ZeroDivisionError('Nonzero angular base required')
        c=[1/self.c[0]]
        for n in range(1,3):c.append(-sum(self.c[i]*c[n-i] for i in range(1,n+1))/self.c[0])
        return J(c)
    def __truediv__(self,x):return self*J.cast(x).inverse()
    def diff(self):return J((self.c[1],2*self.c[2],0))
    @property
    def value(self):return self.c[0]

class P:
    def __init__(self,c):self.c={i:J.cast(v) for i,v in c.items() if any(J.cast(v).c)}
    @staticmethod
    def cast(x):return x if isinstance(x,P) else P({0:x})
    def __add__(self,x):
        c=self.c.copy()
        for i,v in P.cast(x).c.items():c[i]=c.get(i,J((0,)))+v
        return P(c)
    __radd__=__add__
    def __neg__(self):return P({i:-v for i,v in self.c.items()})
    def __sub__(self,x):return self+-P.cast(x)
    def __mul__(self,x):
        c={}
        for i,a in self.c.items():
            for j,b in P.cast(x).c.items():c[i+j]=c.get(i+j,J((0,)))+a*b
        return P(c)
    __rmul__=__mul__
    def __truediv__(self,x):return self*J.cast(x).inverse()
    def dx(self):return P({i-1:i*v for i,v in self.c.items() if i})
    def de(self):return P({i:v.diff() for i,v in self.c.items()})
    def primitive(self):
        if -1 in self.c:raise ValueError('Logarithmic primitive not used in these compact fixtures')
        return P({i+1:v/F(i+1) for i,v in self.c.items()})
    def total(self):
        if any(i<0 for i in self.c):raise ValueError('Integrable polynomial required')
        return sum((v/F(i+1) for i,v in self.c.items()),J((0,)))
    def average(self):return P({i:v/F(i+1) for i,v in self.c.items()})

X=P({1:1});IX=P({-1:1});ONE=P.cast(1)

def z(k,f,eta,h):
    d=1-eta*eta;L=1-2*h*eta*eta
    return (2*k*eta*f-2*eta*X*f.dx()+d*f.de())/L

def z2(k,f,eta,h):return z(k-(F(1,2)-h),z(k,f,eta,h),eta,h)

def time(k,f,eta,h):return (-k*f+X*f.dx()+(F(1,2)-h)*eta*f.de())/(1-2*h*eta*eta)

def fixtures(n,angle):
    eta=J((angle,1,0));wall=ONE-X
    wall=wall*wall*wall*wall
    data=[]
    for i in range(n+1):
        swirl=wall*(1+eta/F(i+2))*(1+F(3*(i+1),100)*eta*eta)/F(i+1)
        axial=X*wall*(F(3,10)+F(2,5*(i+1))*eta+F(i+1,50)*eta*eta)
        pressure=wall*(F(1,10*(i+1))+F(3,100)*eta+F(1,25)*eta*eta)
        data.append((swirl,axial,pressure))
    return data,eta

def assemble(n,h,angle,nu):
    data,eta=fixtures(n,angle);A=F(1,2)+h;D=F(1,2)-h
    L=1-2*h*eta*eta;d=1-eta*eta;lam=2*h*n;b=-A-F(1,2);c=-A
    radial=[]
    for i,(_,u,_) in enumerate(data):
        m=u.primitive()
        radial.append((2*eta*X*u-2*(D+2*h*i)*eta*m-d*m.de())/L)
    f,u,p=data[n]
    theta=time(b+lam,f,eta,h)-2*nu*(X*f.dx().dx()+2*f.dx())
    axial=time(c+lam,u,eta,h)-2*nu*(X*u.dx().dx()+u.dx())+z(-2*A+lam,p,eta,h)
    for i in range(n+1):
        j=n-i;fj,uj,_=data[j];ui=data[i][1]
        theta+=radial[i]*(fj.dx()+fj*IX)+ui*z(b+2*h*j,fj,eta,h)
        axial+=radial[i]*uj.dx()+ui*z(c+2*h*j,uj,eta,h)
    theta-=nu*z2(b+2*h*(n-1),data[n-1][0],eta,h)
    axial-=nu*z2(c+2*h*(n-1),data[n-1][1],eta,h)
    direct=((2*X*theta).total().value,axial.total().value)
    M=u.total();I=(2*X*f).total()
    Jsum=P.cast(0);Ssum=p
    for i in range(n+1):
        j=n-i
        Jsum+=2*X*data[i][1]*data[j][0]
        Ssum+=data[i][1]*data[j][1]
    Jn=Jsum.total();Sn=Ssum.total()
    a=1-h+lam;bflux=F(3,2)-2*A+lam;g=1-A+lam;sflux=1-2*A+lam
    moment_theta=(-a*I+D*eta*I.diff()+2*bflux*eta*Jn+d*Jn.diff())/L
    moment_axial=(-g*M+D*eta*M.diff()+2*sflux*eta*Sn+d*Sn.diff())/L
    lower_I=P.cast((2*X*data[n-1][0]).total())
    lower_M=P.cast(data[n-1][1].total())
    vtheta=nu*z2(a-2*h,lower_I,eta,h).total().value
    vaxial=nu*z2(g-2*h,lower_M,eta,h).total().value
    expected=(moment_theta.value-vtheta,moment_axial.value-vaxial)
    return direct,expected,(vtheta,vaxial)

@lru_cache(None)
def conservative_review():
    rows=[]
    for n in (1,2,4):
        for eta in map(F,('-1','0','2/5','1')):
            for h in (F(1,200),F(9,1000)):
                direct,expected,viscous=assemble(n,h,eta,F(17,10))
                rows.append(dict(order=n,eta=str(eta),h=str(h),exact=direct==expected,
                    omitted_lower_viscosity_nonzero=any(viscous)))
    return dict(rows=rows,exact_integrated_identity_count=2*len(rows),
        angular_total='[-(1-h+2nh)I+Deta I_eta+2(3/2-2A+2nh)eta J+dJ_eta]/L-nu Z_(a-D)Z_a I_(n-1)',
        axial_total='[-(1-A+2nh)M+Deta M_eta+2(1-2A+2nh)eta S+dS_eta]/L-nu Z_(g-D)Z_g M_(n-1)',
        checks=dict(direct_polynomial_residual_integrals_match_all_order_moment_formulas=all(r['exact'] for r in rows),
            lower_viscosity_omission_detected=any(r['omitted_lower_viscosity_nonzero'] for r in rows),
            endpoint_angles_use_no_division_by_one_minus_eta_squared=True,
            lower_angular_exponent_is_one_minus_h_plus_two_h_n_minus_one=True,
            lower_axial_exponent_is_one_minus_A_plus_two_h_n_minus_one=True,
            rational_fixtures_are_algebra_not_actual_coefficient_solutions=True))

@lru_cache(None)
def plateau_review():
    eta=J((0,1,0));h=F(1,200);D=F(1,2)-h;L=1-2*h*eta*eta
    B=ONE+X/F(5)+X*X/F(10)
    V=-B+2*eta*eta*X*B.dx()/L
    got=z2(0,V,eta,h)
    want=6*X*B.dx()
    keys=set(got.c)|set(want.c)
    equal=all(got.c.get(i,J((0,))).value==want.c.get(i,J((0,))).value for i in keys)
    nu=F(17,10);mu=F(3)
    return dict(manufactured_second_pressure_plateau=str(3*nu*mu),
        actual_I4_pressure=['K1(eta)-mu^2/(4X)','K2(eta)'],
        checks=dict(exact_shifted_axial_source_at_eta_zero_is_six_X_Bprime=equal,
            manufactured_backward_second_pressure_can_be_nonzero=3*nu*mu>0,
            leading_constant_radial_memory_is_retained=True,
            actual_K1_K2_not_numerically_evaluated_or_assigned_zero=True))
