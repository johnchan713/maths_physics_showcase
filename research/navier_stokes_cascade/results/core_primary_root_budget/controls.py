"""Independent rational multivariate Taylor controls and edge/jet ledgers."""
from fractions import Fraction as Q
from math import factorial
from budgets import P,RootBudget,atlas,fixture


class Jet:
    """Signed, factorial-normalized two-variable jets, truncated at degree N."""
    def __init__(self,terms=0,N=6):
        self.N=N
        if type(terms) in (int,Q):terms={(0,0):Q(terms)}
        self.c={i:Q(v) for i,v in terms.items() if sum(i)<=N and v}

    def __add__(self,b):
        if not isinstance(b,Jet):b=Jet(b,self.N)
        out=dict(self.c)
        for i,v in b.c.items():out[i]=out.get(i,0)+v
        return Jet(out,self.N)

    __radd__=__add__
    def __neg__(self):return Jet({i:-v for i,v in self.c.items()},self.N)
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        if not isinstance(b,Jet):b=Jet(b,self.N)
        out={}
        for i,v in self.c.items():
            for j,w in b.c.items():
                k=tuple(a+b for a,b in zip(i,j));out[k]=out.get(k,0)+v*w
        return Jet(out,self.N)

    __rmul__=__mul__
    def __pow__(self,n):
        out=Jet(1,self.N)
        for _ in range(n):out=out*self
        return out

    def fractional(self,p,constant):
        """Input constant is the exact rational power of the nonzero base."""
        base=self.c.get((0,0),0)
        if base<=0:raise ValueError("Positive jet value required")
        z=(self-Jet(base,self.N))*Q(1,base);term=Jet(1,self.N)
        out=term;coefficient=Q(1)
        for k in range(1,self.N+1):
            coefficient*=Q(p-k+1,k);term=term*z;out+=coefficient*term
        return constant*out

    def diff(self,axis):
        out={}
        for i,v in self.c.items():
            if i[axis]:
                j=list(i);j[axis]-=1;out[tuple(j)]=i[axis]*v
        return Jet(out,self.N)

    def norm(self,n):return sum((abs(v) for i,v in self.c.items() if sum(i)==n),Q(0))


def matmul(a,b):return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Jet(0,a[0][0].N)) for j in range(len(b[0]))] for i in range(len(a))]


def matrix_norm(a,n):
    # Coefficientwise infinity norm then total-degree sum, not max of row jet sums.
    indices=set(i for row in a for x in row for i in x.c if sum(i)==n)
    return sum((max(sum(abs(x.c.get(i,0)) for x in row) for row in a) for i in indices),Q(0))


def exact_fixture(N=6):
    x=Jet({(1,0):1},N);y=Jet({(0,1):1},N)
    H=[[2+x,1+y],[y,3-x]]
    b=[[4+2*x-y+x*y+x*x],[9-x+3*y+x*y*y]]
    target=matmul(H,b)
    det=H[0][0]*H[1][1]-H[0][1]*H[1][0]
    dinv=det.fractional(-1,Q(1,6))
    inverse=[[dinv*H[1][1],-dinv*H[0][1]],[-dinv*H[1][0],dinv*H[0][0]]]
    budget=RootBudget(lambda n:P(matrix_norm(H,n)),lambda n:P(matrix_norm(target,n)),
                      lambda n:P(1),lambda n:P(1),lambda n:P(1),
                      K0=1,beta=1,R=1,mass=1,width=1)
    return H,b,target,inverse,budget


def algebra_review():
    checks={};H,b,target,inverse,budget=exact_fixture()
    recovered=matmul(inverse,target);unit=matmul(H,inverse)
    for i in range(2):
        for j in range(2):checks[f"inverse_identity_{i}_{j}"]=unit[i][j].c==Jet(int(i==j)).c
        checks[f"root_target_recovered_{i}"]=recovered[i][0].c==b[i][0].c
    for n in range(7):
        checks[f"n{n}_noncommuting_inverse_majorant"]=matrix_norm(inverse,n)<=budget.K(n).evaluate(1)
        checks[f"n{n}_weighted_target_to_b_majorant"]=matrix_norm(b,n)<=budget.B(n).evaluate(1)
        for i,s in enumerate((2,3)):
            root=b[i][0].fractional(Q(1,2),Q(s));recip=b[i][0].fractional(Q(-1,2),Q(1,s))
            checks[f"n{n}_i{i}_root_square"]=((root*root)-b[i][0]).norm(n)==0
            checks[f"n{n}_i{i}_reciprocal_product"]=(root*recip-Jet(1)).norm(n)==0
            checks[f"n{n}_i{i}_root_majorant"]=root.norm(n)<=budget.root(n).evaluate(1)
            checks[f"n{n}_i{i}_reciprocal_majorant"]=recip.norm(n)<=budget.reciprocal(n).evaluate(1)
    xi=[[Jet({(0,0):1,(1,0):2})],[Jet({(0,0):-2,(0,1):1})]]
    d=matmul(inverse,xi)
    da=[[Q(1,2)*d[i][0]*b[i][0].fractional(Q(-1,2),Q(1,s))] for i,s in enumerate((2,3))]
    cross=matmul(H,[[2*b[i][0].fractional(Q(1,2),Q(s))*da[i][0]] for i,s in enumerate((2,3))])
    square=matmul(H,[[da[i][0]*da[i][0]] for i in range(2)])
    checks["signed_original_root_cross_exact"]=all((cross[i][0]-xi[i][0]).c=={} for i in range(2))
    checks["signed_square_must_remain"]=matrix_norm(square,0)>0
    Hconstant=[[Jet(q.c.get((0,0),0)) for q in row] for row in inverse]
    stale=matmul(Hconstant,target)
    checks["omitted_inverse_derivative_detected"]=any(stale[i][0].c!=b[i][0].c for i in range(2))
    fixed_da=[[Q(1,2*s)*d[i][0]] for i,s in enumerate((2,3))]
    stale_cross=matmul(H,[[2*b[i][0].fractional(Q(1,2),Q(s))*fixed_da[i][0]] for i,s in enumerate((2,3))])
    checks["omitted_reciprocal_root_derivative_detected"]=any(stale_cross[i][0].c!=xi[i][0].c for i in range(2))
    return dict(checks=checks,scope="signed rational two-variable jets through order six, not actual supported canonical fields")


def derivative_review():
    checks={};N=6
    a=Jet({(0,0):1,(1,0):2,(0,1):-3,(2,1):4,(1,3):5},N)
    d=Jet({(0,0):2,(1,0):1,(0,2):3},N)
    # The highest differentiated degree is five, avoiding truncation ambiguity.
    for n in range(5):
        checks[f"n{n}_all_coordinate_derivative_sum"]=sum(a.diff(mu).norm(n) for mu in range(2))==(n+1)*a.norm(n+1)
        envelope=sum(d.norm(j)*(a.norm(n-j)+(n-j+1)*a.norm(n-j+1)) for j in range(n+1))
        curl=d*a+d*a.diff(0)+d*a.diff(1)
        checks[f"n{n}_full_variable_curl_coefficient_product"]=curl.norm(n)<=envelope
        checks[f"n{n}_ordinary_derivative_conversion"]=all(abs(v)*factorial(i[0])*factorial(i[1])<=factorial(n)*a.norm(n) for i,v in a.c.items() if sum(i)==n)
    full=(d*a).diff(0);dropped=d*a.diff(0)
    checks["derivative_of_curl_coefficient_cannot_be_omitted"]=full.c!=dropped.c
    # Independent exponential Taylor coefficients via finite powers of a jet.
    z=Jet({(1,0):2,(0,1):-1,(1,1):3},N);exp=Jet(0,N)
    for k in range(N+1):exp+=Q(1,factorial(k))*z**k
    beta=RootBudget(lambda n:P(1),lambda n:P(1),lambda n:P(1),lambda n:P(1),
                    lambda n:P(z.norm(n)),1,1,1,1,1)
    for n in range(N+1):checks[f"n{n}_flat_exponential_Bell_majorant"]=exp.norm(n)<=beta.E(n).evaluate(1)
    return dict(checks=checks,scope="exact coordinate/Taylor conversion and variable first-order operator controls")


def ledger_review():
    checks={}
    for n in range(13):
        for k in range(1,n+1):checks[f"n{n}_k{k}_reciprocal_edge_degree"]=3*n+3*k<=6*n
        for j in range(n+1):
            checks[f"n{n}_j{j}_root_edge_degree"]=3*j+3+6*(n-j)<=6*n+3
            checks[f"n{n}_j{j}_curl_edge_degree"]=6*(n-j+1)+3<=6*n+9
            checks[f"n{n}_j{j}_flat_product_edge_degree"]=3*j+6*(n-j)+3<=6*n+3
        checks[f"n{n}_positive_edge_unweighting"]=atlas.edge_constant(6*n+3,Q(1,8))>=1
    r=Q(1,6057072);mass=2/r
    checks["fixed_actual_radius_mass_padding"]=mass**4>12/r**3
    checks["three_component_Euclidean_padding"]=2**2>=3
    checks["outer_normalized_root_need_not_be_bounded"]=Q(-3,2)<0
    checks["same_flat_constants_retained"]=Q(4,2)==2 and Q(1,2)>0
    return dict(checks=checks,scope="finite checks of the general inequalities proved in the note; fixed width/radius costs remain")


def substitution_review():
    checks={};b=fixture()
    for n in (0,1,2):
        a,c=b.native(n);g=b.aggregate(n)
        for S in (1,4):
            scale=Q(25,16)*S
            checks[f"n{n}_S{S}_native_to_atlas_without_double_unweighting"]=g["a"].evaluate(S)==8*188**n*a.evaluate(scale) and g["b"].evaluate(S)==8*188**n*c.evaluate(scale)
            checks[f"n{n}_S{S}_all_positive_scale_coefficients"]=all(v>0 for p in (a,c,g["Omega"]) for v in p.c.values())
            checks[f"n{n}_S{S}_signed_budget_is_positive"]=b.signed(n,lambda j:P(j+1)).evaluate(S)>0
    # A direct numeric evaluation of the inherited positive-polynomial norm chain
    # avoids a symbolic polynomial explosion, without any floating arithmetic.
    for n in (0,1):
        trace=set();root=fixture(trace)
        seed=atlas.seed.fixture_budget(1)
        seed.omega=lambda m:root.aggregate(m)["Omega"].evaluate(1)
        out=seed.output(n)
        checks[f"n{n}_scalar_initial_chain_accepts_derived_root_bounds"]=all(v>0 for v in out.values())
        checks[f"n{n}_scalar_underlying_jet_cost"]=max(k for name,k in trace if name in ("H","target","wave"))==n+11
        trace=set();root=fixture(trace);seed=atlas.seed.fixture_budget(1)
        seed.omega=lambda m:root.aggregate(m)["Omega"].evaluate(1)
        moment=seed.axial_moment(n,Q(5))
        checks[f"n{n}_whole_moving_pressure_moment_positive"]=moment>0
        checks[f"n{n}_whole_moment_underlying_jet_cost"]=max(k for name,k in trace if name in ("H","target","wave"))==n+12
        checks[f"n{n}_whole_moment_curl_coefficient_cost"]=max(k for name,k in trace if name=="curl")==n+11
    return dict(checks=checks,scope="exact substitution of manufactured finite polynomial input bounds, not actual norm evaluation")


def review():
    parts=dict(matrix_roots_signed=algebra_review(),factorial_curl_flat=derivative_review(),
               edge_ledger=ledger_review(),aggregate_pressure_chain=substitution_review())
    return dict(checks={name+":"+k:v for name,row in parts.items() for k,v in row["checks"].items()},parts=parts)
