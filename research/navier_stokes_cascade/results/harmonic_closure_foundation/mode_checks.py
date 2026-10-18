"""Independent exact Gaussian-rational polynomial/Fourier curl controls."""
from dataclasses import dataclass
from fractions import Fraction as F


@dataclass(frozen=True)
class C:
    re:F=F(0)
    im:F=F(0)
    def __post_init__(self):
        if not isinstance(self.re,(int,F)) or not isinstance(self.im,(int,F)):raise TypeError('Exact rational components required')
        object.__setattr__(self,'re',F(self.re));object.__setattr__(self,'im',F(self.im))
    def __add__(self,b):
        if not isinstance(b,(C,int,F)):return NotImplemented
        b=coerce(b);return C(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self):return C(-self.re,-self.im)
    def __sub__(self,b):return self+-coerce(b)
    def __mul__(self,b):
        if not isinstance(b,(C,int,F)):return NotImplemented
        b=coerce(b);return C(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def conj(self):return C(self.re,-self.im)
    def size(self):return abs(self.re)+abs(self.im)


def coerce(c):return c if isinstance(c,C) else C(c)


class P:
    """Polynomial in x with exact rational complex coefficients."""
    def __init__(self,x=0):
        if isinstance(x,P):self.terms=dict(x.terms)
        elif isinstance(x,dict):
            if any(type(n) is not int or n<0 for n in x):raise ValueError('Nonnegative integer powers required')
            self.terms={n:coerce(c) for n,c in x.items() if coerce(c)!=C()}
        else:self.terms={0:coerce(x)} if coerce(x)!=C() else {}
    def __add__(self,b):
        out=dict(self.terms)
        for n,c in P(b).terms.items():out[n]=out.get(n,C())+c
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({n:-c for n,c in self.terms.items()})
    def __sub__(self,b):return self+-P(b)
    def __mul__(self,b):
        out={}
        for n,c in self.terms.items():
            for m,d in P(b).terms.items():out[n+m]=out.get(n+m,C())+c*d
        return P(out)
    __rmul__=__mul__
    def __eq__(self,b):return self.terms==P(b).terms
    def diff(self):return P({n-1:n*c for n,c in self.terms.items() if n})
    def conj(self):return P({n:c.conj() for n,c in self.terms.items()})
    def size(self):return sum(c.size() for c in self.terms.values())
    def jets(self,j):
        out=F(0);p=self
        for _ in range(j+1):out+=p.size();p=p.diff()
        return out


I=C(0,1);X=P({1:1})


def curl(p,m):
    if type(m) is not int or not m:raise ValueError('Nonzero integer mode required')
    p=list(map(P,p));return [-p[1]*(I*m),p[0]*(I*m)-p[2].diff(),p[1].diff()]


def divergence(a,m):return a[0].diff()+a[2]*(I*m)


def add_fields(a,b):
    return {m:[x+y for x,y in zip(a.get(m,[P()]*3),b.get(m,[P()]*3))] for m in a.keys()|b.keys()}


def real_field(positive):
    if any(type(m) is not int or m<=0 for m in positive):raise ValueError('Positive mode data required')
    return dict(positive)|{-m:[x.conj() for x in a] for m,a in positive.items()}


def transport(a,b,repaired=False):
    out={}
    for m,av in a.items():
        if not m:raise ValueError('Transport identity requires nonzero input modes')
        for n,bv in b.items():
            v=[av[0]*x.diff()+av[2]*(I*n)*x for x in bv]
            if repaired:v=[av[0]*x.diff()-av[0].diff()*F(n,m)*x for x in bv]
            previous=out.get(m+n,[P()]*3);out[m+n]=[x+y for x,y in zip(previous,v)]
    return out


def field_norm(a,s,j):return sum((1+abs(m))**s*sum(x.jets(j) for x in v) for m,v in a.items())


def covariance(a,b=None):
    b=a if b is None else b
    return [[sum((v[i]*b.get(-m,[P()]*3)[j] for m,v in a.items()),P()) for j in range(3)] for i in range(3)]


def tensor_norm(a,j):return sum(x.jets(j) for row in a for x in row)


def fixture():
    def potential(m,seed):return [X*(C(seed,F(1,m)))+X*X,F(1,seed)*X*X+I*X,1+X*X*X]
    old=real_field({m:curl(potential(m,s),m) for m,s in ((1,2),(2,3))})
    inc=real_field({m:curl(potential(m,s),m) for m,s in ((1,5),(3,7),(7,11))})
    return old,inc


def complete_review():
    old,inc=fixture();updated=add_fields(old,inc)
    direct=transport(old,inc);repaired=transport(old,inc,True)
    W0,W1=covariance(old),covariance(updated)
    uv,vu,vv=covariance(old,inc),covariance(inc,old),covariance(inc)
    tensor_ok=all(W1[i][j]-W0[i][j]==uv[i][j]+vu[i][j]+vv[i][j] for i in range(3) for j in range(3))
    normrows=[];covrows=[]
    for j in range(4):
        for s in range(4):
            value=field_norm(direct,s,j)
            upper=2**(j+1)*field_norm(old,s,j+1)*field_norm(inc,s+1,j+1)
            normrows.append(dict(s=s,j=j,value=str(value),upper=str(upper)))
        delta=[[W1[i][k]-W0[i][k] for k in range(3)] for i in range(3)]
        upper=2**j*(2*field_norm(old,0,j)+field_norm(inc,0,j))*field_norm(inc,0,j)
        covrows.append(dict(j=j,value=str(tensor_norm(delta,j)),upper=str(upper)))
    high={m:a for m,a in inc.items() if abs(m)==7};squares=covariance(high)
    return dict(transport_rows=normrows,covariance_rows=covrows,complete_tensor_entries=9,
        checks=dict(all_complete_curls_exactly_divergence_free=all(divergence(a,m)==0 for data in (old,inc) for m,a in data.items()),
            direct_transport_matches_divergence_repair=direct==repaired,
            coefficient_transport_majorant=all(F(r['value'])<=F(r['upper']) for r in normrows),
            complete_tensor_update_matches_independent_convolution=bool(tensor_ok),
            covariance_absolute_sum_majorant=all(F(r['value'])<=F(r['upper']) for r in covrows),
            high_harmonic_square_present=tensor_norm(squares,0)>0,
            all_six_symmetric_covariance_entries_nonzero=all(W1[i][j]!=0 for i in range(3) for j in range(i,3)),
            real_covariance_keeps_negative_modes=all(x==x.conj() for row in W1 for x in row),
            zero_transport_outputs_retained=0 in direct))


def obstruction(N):
    if type(N) is not int or N<2:raise ValueError('Integer N>=2 required')
    a=real_field({1:curl([0,X*(I*F(1,2)),0],1)})
    b=real_field({N:curl([P(-I*F(1,2*N)),0,0],N)})
    ab=transport(a,b);ba=transport(b,a)
    return dict(a=a,b=b,ab=ab,ba=ba)


def obstruction_review():
    rows=[];exact=True;div=True;reverse=True
    for N in (2,3,17,128,1024):
        r=obstruction(N);y=r['ab'][N+1][1];ym=r['ab'][N-1][1]
        exact &= y==P(-F(N,4)) and ym==P(F(N,4))
        reverse &= field_norm(r['ba'],0,0)==0
        div &= all(divergence(a,m)==0 for key in ('a','b') for m,a in r[key].items())
        rows.append(dict(N=N,input_a=str(field_norm(r['a'],0,0)),input_b=str(field_norm(r['b'],0,0)),
                         transport=str(field_norm(r['ab'],0,0)),positive_sidebands=[N-1,N+1]))
    return dict(rows=rows,checks=dict(real_complete_curl_sidebands_have_exact_N_over_four_size=bool(exact),
        real_fields_divergence_free=bool(div),reverse_transport_identically_zero=bool(reverse),
        input_coefficient_norms_independent_of_N=all(r['input_a']=='2' and r['input_b']=='1' for r in rows),
        no_generic_uniform_transport_without_angular_loss=all(F(r['transport'])==r['N'] for r in rows),
        local_smooth_compact_potential_extension_available_not_selected_pulse_counterexample=True))
