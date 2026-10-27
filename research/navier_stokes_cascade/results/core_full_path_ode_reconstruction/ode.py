"""Exact controls for the reconstructed full-path ODE contract."""
from fractions import Fraction as F
from math import factorial

def dot(a,b): return sum((F(x)*F(y) for x,y in zip(a,b)),F(0))
def add(*vs): return tuple(sum((F(v[i]) for v in vs),F(0)) for i in range(len(vs[0])))
def scale(a,v): return tuple(F(a)*F(x) for x in v)
def matvec(A,v): return tuple(dot(row,v) for row in A)

def projected_system(n,nprime,K,f,t,k,m,d):
    n,nprime,f,t=map(lambda x:tuple(map(F,x)),(n,nprime,f,t));K=tuple(tuple(map(F,r)) for r in K)
    k,m,d=map(F,(k,m,d));N=dot(n,n)
    if N<=0 or k==0 or m==0 or dot(n,t)!=0: raise ValueError("Nonzero normal/carrier/mode and tangent datum required")
    Kt=matvec(K,t); pressure=(dot(n,Kt)-dot(nprime,t)+dot(n,f))/(k*m*N)
    Pf=add(f,scale(-dot(n,f)/N,n))
    corr=scale((dot(n,Kt)-dot(nprime,t))/N,n)
    tp=add(scale(-1,Kt),corr,scale(-m*m*d,t),scale(-1,Pf))
    residual=add(tp,Kt,scale(m*m*d,t),scale(-k*m*pressure,n),f)
    return {"pressure_over_i":pressure,"projected_source":Pf,"tprime":tp,
            "constraint_derivative":dot(nprime,t)+dot(n,tp),"residual":residual}

def ambient_lift(n,nprime,K,d,rate):
    n,nprime=tuple(map(F,n)),tuple(map(F,nprime));K=tuple(tuple(map(F,r)) for r in K);d,rate=map(F,(d,rate));N=dot(n,n);h=scale(1/N,n)
    nK=tuple(sum(n[r]*K[r][c] for r in range(len(n))) for c in range(len(n)))
    A=tuple(tuple(-K[r][c]+h[r]*(nK[c]-nprime[c])-(d if r==c else 0) for c in range(len(n))) for r in range(len(n)))
    H=tuple(tuple(h[r]*n[c] for c in range(len(n))) for r in range(len(n)))
    Ah=tuple(tuple(A[r][c]+(d+rate)*H[r][c] for c in range(len(n))) for r in range(len(n)))
    row=tuple(nprime[c]+sum(n[r]*Ah[r][c] for r in range(len(n))) for c in range(len(n)))
    tangent=(n[1],-n[0],F(0))
    return {"normal_row":row,"expected":scale(rate,n),"tangent_agreement":matvec(Ah,tangent)==matvec(A,tangent)}

def picard_factorial_threshold(M,L):
    M,L=map(F,(M,L))
    if min(M,L)<=0: raise ValueError("Positive finite coefficient and length bounds required")
    q=2*M*L
    j=max(0,(q.numerator+q.denominator-1)//q.denominator-1)
    # The ratio is the exact convergence certificate.  Do not construct the
    # preceding factorial term: for a large but finite coefficient bound that
    # integer can have millions of digits and is irrelevant to the proof.
    return {"index":j,"next_ratio":F(M*L,j+1)}

def residence(length,mode,dstar):
    length,dstar=map(F,(length,dstar));mode=int(mode)
    if length<=0 or dstar<=0 or mode==0: raise ValueError("Positive path, damping and nonzero mode required")
    delta=mode*mode-1
    return length if delta==0 else min(length,F(1,delta*dstar))

def product_jets(a,b):
    if len(a)!=len(b): raise ValueError("Equal jet lengths required")
    return tuple(sum(F(factorial(n),factorial(j)*factorial(n-j))*F(a[j])*F(b[n-j]) for j in range(n+1)) for n in range(len(a)))

def review():
    ode=projected_system((2,1,1),(1,-1,0),((0,-2,1),(3,0,-1),(1,2,0)),(3,-2,5),(1,-2,0),7,-2,F(3,5))
    normal=projected_system((0,0,1),(0,0,0),((0,0,0),)*3,(0,0,2),(1,0,0),5,1,1)
    lift=ambient_lift((2,1,1),(1,-1,0),((0,-2,1),(3,0,-1),(1,2,0)),F(3,5),F(4,7))
    picard=picard_factorial_threshold(7,11)
    cutA=product_jets((0,2,-3,5),(7,11,13,17));cutP=product_jets((0,2,-3,5),(19,23,29,31))
    return {"checks":{
      "projected_system_preserves_constraint":ode["constraint_derivative"]==0,
      "projected_and_original_systems_agree":all(x==0 for x in ode["residual"]),
      "pure_normal_source_has_zero_projected_velocity":all(x==0 for x in normal["projected_source"]),
      "pure_normal_source_pressure_is_nonzero":normal["pressure_over_i"]!=0,
      "ambient_lift_has_exact_gaussian_normal_response":lift["normal_row"]==lift["expected"] and lift["tangent_agreement"],
      "picard_series_converges_without_smallness":picard["next_ratio"]<=F(1,2),
      "primary_mode_retains_full_path_length":residence(13,1,F(1,7))==13,
      "higher_mode_has_inverse_square_residence":residence(13,3,F(1,7))==F(7,8),
      "potential_cutoff_derivatives_are_nonzero":cutA[1]!=0,
      "full_pressure_cutoff_derivatives_are_nonzero":cutP[1]!=0,
    },"controls":{"pressure_over_i":str(ode["pressure_over_i"]),"normal_source_pressure_over_i":str(normal["pressure_over_i"]),
      "picard":{k:str(v) for k,v in picard.items()},"potential_cutoff_jets":[str(x) for x in cutA],"pressure_cutoff_jets":[str(x) for x in cutP]},
      "scope":{"canonical_primary_is_homogeneous":True,"sourced_inverse_retains_unprojected_pressure":True,
      "raw_exit_response_declared_zero":False,"fixed_projector_contract_closed":False,
      "unconditional_pulse_residual_force_or_blowup":False}}
