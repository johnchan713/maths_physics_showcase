"""Five exact moment equations, represented numerically by positive quadrature.

The mathematical map uses exact integrals of fixed smooth bumps. A numerical
root of their quadrature approximation is a diagnostic, not an exact profile.
Columns are two additive U bumps and three relative E bumps. Rows are the
normalized differences M, J-4 eta I, I, -(S-8 eta M), Cp, in that order.
"""
from dataclasses import dataclass
from pathlib import Path
import math
import sys
import mpmath as mp

sys.path.append(str(Path(__file__).resolve().parent.parent/'outer_pressure_pilot'))
from schedule import Rule,step,step_prime


@dataclass(frozen=True)
class Jet:
    """A value and its first eta derivative; products keep the product rule."""
    v: object
    d: object = 0

    def __post_init__(self):
        object.__setattr__(self,'v',mp.mpf(self.v))
        object.__setattr__(self,'d',mp.mpf(self.d))

    @staticmethod
    def cast(x):
        return x if isinstance(x,Jet) else Jet(x)

    def __add__(self,other):
        b=self.cast(other)
        return Jet(self.v+b.v,self.d+b.d)

    __radd__=__add__

    def __neg__(self):
        return Jet(-self.v,-self.d)

    def __sub__(self,other):
        return self+-self.cast(other)

    def __rsub__(self,other):
        return self.cast(other)+-self

    def __mul__(self,other):
        b=self.cast(other)
        return Jet(self.v*b.v,self.d*b.v+self.v*b.d)

    __rmul__=__mul__

    def __truediv__(self,other):
        b=self.cast(other)
        return Jet(self.v/b.v,(self.d*b.v-self.v*b.d)/(b.v*b.v))

    def __rtruediv__(self,other):
        return self.cast(other)/self

    def __pow__(self,power):
        return Jet(self.v**power,power*self.v**(power-1)*self.d)


def fixture(eta,name='odd',epsilon='1e-16'):
    """Manufactured smooth incoming data, explicitly not an axis solution.

    Each component and g has C1 norm below epsilon on [-1,1]. The odd
    fixture vanishes at eta=0 but its derivative does not vanish there.
    """
    e=Jet(eta,1);eps=mp.mpf(epsilon)
    if abs(e.v)>1 or eps<=0:
        raise ValueError('Invalid fixture domain')
    if name=='odd':
        rows=[e/32,(e+e**3)/64,(e-e**3/3)/32,
              Jet(mp.sin(e.v),mp.cos(e.v))/32,e*(1+e**2)/128]
        g=e/64
    elif name=='edge':
        f=1/(1+e**2)
        rows=[f/2,(1+e/4)/2,-f/3,(1-e/4)/2,f/4]
        g=Jet(mp.sin(e.v),mp.cos(e.v))/4
    else:
        raise ValueError('Unknown manufactured fixture')
    return dict(name=name,entry=[eps*r for r in rows],g=eps*g,
                is_axis_solution=False,epsilon=eps)


def zeta(eta,P):
    """Coefficient of the additive-U square in the normalized energy row."""
    e=Jet(eta,1)
    return mp.exp(mp.mpf('1.2'))*(1+e**2)**2/mp.mpf(P)**2


def step_second(z):
    if not 0<z<1:
        return mp.mpf(0)
    s,c=step(z),step(1-z)
    g=2/z**3+2/(1-z)**3
    return s*c*((c-s)*g*g-6/z**4+6/(1-z)**4)


class MomentMap:
    U_CENTERS=('.3','.7')
    E_CENTERS=('.15','.5','.85')
    WIDTH='.08'

    def __init__(self,order=32,panels=8):
        self.rule=Rule(order,panels)
        self.width=mp.mpf(self.WIDTH)
        self.centers=list(map(mp.mpf,self.U_CENTERS+self.E_CENTERS))
        ordered=sorted(self.centers)
        if ordered[0]-self.width/2<=0 or ordered[-1]+self.width/2>=1 \
            or any(b-a<=self.width for a,b in zip(ordered,ordered[1:])):
            raise ValueError('All five bump supports must be disjoint and interior')
        self._integrals={};self._tables={};self._restoration={}
        self.B,self.QU,self.QE,self.QC=self.tables(mp.mpf(1))
        self.inverse=self.B**-1

    def bump(self,index,t,radial_derivative=False):
        z=(mp.mpf(t)-self.centers[index])/self.width+mp.mpf('.5')
        return step_second(z)/self.width if radial_derivative else step_prime(z)

    def integral(self,index,rate,power=1,t=1):
        rate=mp.mpf(rate);t=mp.mpf(t)
        upper=min(mp.mpf(1),(t-self.centers[index])/self.width+mp.mpf('.5'))
        if upper<=0:
            return mp.mpf(0)
        key=(index,rate,power,upper)
        if key not in self._integrals:
            c,w=self.centers[index],self.width
            self._integrals[key]=w*self.rule.integrate(
                lambda z:mp.exp(rate*(c+w*(z-mp.mpf('.5'))))*step_prime(z)**power,0,upper)
        return self._integrals[key]

    def tables(self,t):
        t=min(mp.mpf(1),mp.mpf(t))
        if t in self._tables:
            return self._tables[t]
        B=mp.matrix(5,5)
        for j in range(2):
            B[0,j]=self.integral(j,'1',t=t)
            B[1,j]=self.integral(j,'1.6',t=t)
        for j in range(2,5):
            for row,rate in ((2,'1.6'),(3,'1.2'),(4,'.2')):
                B[row,j]=self.integral(j,rate,t=t)
        QU=[self.integral(j,'1',2,t) for j in range(2)]
        QE=[self.integral(j,'1.2',2,t) for j in range(2,5)]
        QC=[self.integral(j,'.2',2,t) for j in range(2,5)]
        self._tables[t]=(B,QU,QE,QC)
        return self._tables[t]

    def restoration(self,t=0):
        """Integrate the axial restoration on -2<t<-1, where t=y+6."""
        upper=min(mp.mpf(1),mp.mpf(t)+2)
        if upper<=0:
            return [mp.mpf(0)]*3
        if upper not in self._restoration:
            self._restoration[upper]=[self.rule.integrate(
                lambda z:mp.exp(mp.mpf(rate)*(z-2))*step(1-z)**power,0,upper)
                for rate,power in (('1',1),('1.6',1),('1',2))]
        return self._restoration[upper]

    def incoming(self,data,eta,P,t=0):
        a,b,c=self.restoration(t)
        g=data['g'];r=list(data['entry'])
        r[0]+=a*g;r[1]+=b*g;r[3]-=c*zeta(eta,P)*g*g
        return r

    def quadratic(self,c,d,eta,P,t=1):
        _,QU,QE,QC=self.tables(t)
        result=mp.matrix(5,1)
        result[3]=mp.fsum(QE[j]*c[j+2]*d[j+2]/2 for j in range(3)) \
            -zeta(eta,P).v*mp.fsum(QU[j]*c[j]*d[j] for j in range(2))
        result[4]=mp.fsum(QC[j]*c[j+2]*d[j+2]/2 for j in range(3))
        return result

    def quadratic_eta(self,c,eta,P):
        result=mp.matrix(5,1)
        result[3]=-zeta(eta,P).d*mp.fsum(self.QU[j]*c[j]**2 for j in range(2))
        return result

    def jacobian(self,c,eta,P):
        out=self.B.copy()
        for j in range(5):
            unit=mp.matrix(5,1);unit[j]=1
            column=2*self.quadratic(c,unit,eta,P)
            for i in range(5):out[i,j]+=column[i]
        return out

    def newton(self,target,eta,P):
        target=mp.matrix(target)
        if not any(target):
            return mp.matrix(5,1)
        c=self.inverse*target
        for _ in range(12):
            residual=self.B*c+self.quadratic(c,c,eta,P)-target
            nxt=c-mp.lu_solve(self.jacobian(c,eta,P),residual)
            if all(a==b for a,b in zip(c,nxt)):
                return nxt
            c=nxt
        return c

    def solve(self,data,eta,P):
        target_jets=[-x for x in self.incoming(data,eta,P)]
        target=[x.v for x in target_jets]
        c=self.newton(target,eta,P)
        rhs=mp.matrix([x.d for x in target_jets])-self.quadratic_eta(c,eta,P)
        ce=mp.lu_solve(self.jacobian(c,eta,P),rhs)
        return dict(target=target,target_eta=[x.d for x in target_jets],
                    root=list(c),root_eta=list(ce))

    def graded(self,target,eta,P,degree=4):
        if degree<1:
            raise ValueError('At least one homogeneous grade required')
        target=mp.matrix(target)
        components=[self.inverse*target]
        for n in range(2,degree+1):
            q=mp.matrix(5,1)
            for i in range(1,n):
                q+=self.quadratic(components[i-1],components[n-i-1],eta,P)
            components.append(-self.inverse*q)
        a=1000*max(abs(v) for v in target)
        ratio=4*1000*100*a
        if ratio>=mp.mpf('.5'):
            raise ValueError('Outside the proved small-root majorant')
        remainder=a*ratio**degree/(1-ratio)
        return dict(components=[list(v) for v in components],majorant_a=a,
                    majorant_ratio=ratio,remainder=remainder)

    def ledger(self,data,solution,eta,P,t,frozen_derivative=False):
        r=self.incoming(data,eta,P,t)
        B,QU,QE,QC=self.tables(t)
        c=[Jet(v,0 if frozen_derivative else d) for v,d in
           zip(solution['root'],solution['root_eta'])]
        for i in range(5):r[i]+=sum(B[i,j]*c[j] for j in range(5))
        r[3]+=sum(QE[j]*c[j+2]*c[j+2]/2 for j in range(3)) \
            -zeta(eta,P)*sum(QU[j]*c[j]*c[j] for j in range(2))
        r[4]+=sum(QC[j]*c[j+2]*c[j+2]/2 for j in range(3))
        return r

    def physical_rows(self,rows,eta,P):
        """Undo the row combinations, retaining their original dimensionless units."""
        z=zeta(eta,P).v
        return [rows[0],rows[2],rows[1]+4*eta*rows[2],
                -rows[3]+8*eta*z*rows[0],rows[4]]

    def independent_integrals(self,coefficients,eta,P,max_step):
        """Independent DOP853 integration of the five original moment densities.

        Linear and quadratic densities are integrated in separate units. This
        prevents a tiny quadratic contribution from disappearing beside a
        linear term in binary64 arithmetic. No axis data are inferred here.
        """
        import numpy as np
        from scipy.integrate import solve_ivp
        scale=max(abs(v) for v in coefficients)
        if not scale:
            raise ValueError('Nonzero coefficients required')
        c=np.array([float(v/scale) for v in coefficients])
        e=float(eta);z=float(zeta(eta,P).v)
        centers=list(map(float,self.centers));width=float(self.width)
        def beta(t,center):
            a=(t-center)/width+.5
            if not 0<a<1:return 0.
            b=min(a,1-a)
            ratio=math.exp(-1/b**2+1/(1-b)**2)
            return ratio/(1+ratio)**2*(2/a**3+2/(1-a)**3)
        def rhs(t,_):
            v=[beta(t,a) for a in centers]
            u=sum(c[j]*v[j] for j in range(2))
            ee=sum(c[j]*v[j] for j in range(2,5))
            x1,x16,x12,x02=[math.exp(a*t) for a in (1,1.6,1.2,.2)]
            linear=[x1*u,x16*ee,x16*(4*e*ee+u),
                    z*x1*8*e*u-x12*ee,x02*ee]
            squared=[0.,0.,x16*u*ee,z*x1*u*u-x12*ee*ee/2,x02*ee*ee/2]
            return linear+squared
        result=solve_ivp(rhs,(0.,1.),np.zeros(10),method='DOP853',
                         rtol=2e-13,atol=2e-14,max_step=float(max_step))
        if not result.success:
            raise RuntimeError(result.message)
        cc=mp.matrix([v/scale for v in coefficients])
        expected=self.physical_rows(self.B*cc,eta,P)+self.physical_rows(
            self.quadratic(cc,cc,eta,P),eta,P)
        gaps=[abs(mp.mpf(float(v))-w)/(1+abs(w)) for v,w in zip(result.y[:,-1],expected)]
        return dict(max_step=max_step,gaps=gaps,linear=list(map(float,result.y[:5,-1])),
                    quadratic=list(map(float,result.y[5:,-1])),evaluations=result.nfev)

    def tanh_sinh_quadratic(self,coefficients,eta,P):
        """Second integration rule for the small U-square and E-square terms."""
        scale=max(abs(v) for v in coefficients)
        c=[v/scale for v in coefficients]
        groups=[]
        for indices,rate in ((range(2),'1'),(range(2,5),'1.2'),(range(2,5),'.2')):
            groups.append(mp.fsum(self.width*c[j]**2*mp.quad(
                lambda z:mp.exp(mp.mpf(rate)*(self.centers[j]+self.width*(z-mp.mpf('.5'))))*step_prime(z)**2,
                [0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1]) for j in indices))
        expected=[mp.fsum(self.QU[j]*c[j]**2 for j in range(2)),
                  mp.fsum(self.QE[j]*c[j+2]**2 for j in range(3)),
                  mp.fsum(self.QC[j]*c[j+2]**2 for j in range(3))]
        return dict(gaps=[abs(a-b)/abs(b) for a,b in zip(groups,expected)],
                    U_square=scale**2*zeta(eta,P).v*groups[0],
                    E_energy_square=scale**2*groups[1]/2,
                    E_pressure_square=scale**2*groups[2]/2)


def relative(a,b):
    scale=max(abs(a),abs(b))
    return abs(a-b)/scale if scale else mp.mpf(0)
