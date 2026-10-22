"""Exact response profiles, fractional chart jets and explicit inverse budgets."""
from fractions import Fraction as Q
from math import ceil,factorial
import importlib.util
from pathlib import Path

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


fixed=load("slow_inverse_fixed",HERE.parent/"core_fixed_profile_tame/proof_checks.py")
coupled=load("slow_inverse_coupled",HERE.parent/"core_coupled_defect_budget/proof_checks.py")
maps=load("slow_inverse_maps",HERE.parent/"compact_mean_foundation/mean_maps.py")
means=fixed.means;R=means.R;Jet=fixed.Jet;P=coupled.P
order=fixed.order;exact=fixed.geometry.exact


def composition_bound(p,r,m,profile_majorants):
    p=exact(p);order(r);order(m)
    if len(profile_majorants)!=m+1:raise ValueError("All profile majorants through m required")
    values=[fixed.exact_nonnegative(v) for v in profile_majorants]
    c=r-2*p
    return 2**ceil(abs(c))*sum((values[k]*sum(jet.c.values(),Q(0))/Q(2**k*factorial(k))
                              for k,jet in enumerate(fixed.majorant_powers(m))),Q(0))


def reciprocal_bound(m,a0,positive_jets):
    order(m);a0=exact(a0)
    if a0<=0:raise ValueError("Positive retained background floor required")
    if not isinstance(positive_jets,P):positive_jets=fixed.exact_nonnegative(positive_jets)
    power=positive_jets*0+1;total=power
    for _ in range(m):power=power*positive_jets*Q(1,a0);total=total+power
    return total*Q(1,a0)


def background_positive_bound(m,normalized_profile_derivatives):
    order(m)
    if len(normalized_profile_derivatives)!=m:raise ValueError("Positive profile derivative orders one through m required")
    values=[fixed.exact_nonnegative(v) for v in normalized_profile_derivatives]
    E=(2+(m>=1))*Q(16,3)*sum(32**i*64**j for i in range(m+1) for j in range(m+1-i))
    return sum((value*E**k for k,value in enumerate(values,1)),Q(0))


def inverse_budget(m,H,Hnext,Kp,Kz,Kg,Kpsi_next,pi_p,pi_z,pi_j,pi_j_next):
    order(m)
    return dict(U=H*(Kp*pi_p+Kz*pi_z),Z=H*Kg*pi_j,
                X=(m+1)*Hnext*Kpsi_next*pi_j_next)


def responses(g,v,moments):
    At,Az=maps.blocks(g,v,1,moments)
    return tuple(maps.solve(A,b) for A,b in ((At,[0,-1,0]),(At,[0,0,-1]),(Az,[0,-1])))


def matrix_review():
    checks={};rows=[]
    for g,v,a in ((Q(3,2),Q(2,3),Q(2)),(Q(2),Q(255,256),Q(-3)),(Q(3,2),Q(65535,65536),Q(1,7))):
        mu=(Q(3,2),Q(7,5),Q(4,3),Q(9,7),Q(5,4))
        bp,bz,bg=responses(g,v,mu);At,Az=maps.blocks(g,v,a,mu)
        for p,j,z in ((Q(2),Q(-5),Q(7)),(Q(-3),Q(11),Q(-2))):
            u=[(p*x+z*y)/a for x,y in zip(bp,bz)];gamma=[j*x/a for x in bg]
            lhs=[sum(x*y for x,y in zip(row,u)) for row in At]+[sum(x*y for x,y in zip(row,gamma)) for row in Az]
            checks[f"signed_five_rows_g{g}_v{v}_a{a}_p{p}"]=lhs==[0,-p,-z,0,-j]
        checks[f"axial_gap_cost_g{g}_v{v}"]=(1-v)*bg[1]==1/(mu[4]*g)
        rows.append(dict(g=str(g),v=str(v),a=str(a),bp=list(map(str,bp)),bz=list(map(str,bz)),bg=list(map(str,bg))))
    checks["angular_limit_nonsingular"]=maps.determinant(maps.blocks(2,1)[0])!=0
    checks["axial_limit_singular"]=maps.determinant(maps.blocks(2,1)[1])==0
    for h,lam in ((Q(1,200),Q(1,8)),(Q(1,1000),Q(1,100000))):
        A=Q(1,2)+h;vp=lam-h
        # Exact exponent arithmetic uses the actual moment weight and dR.
        pressure=vp+A+Q(-1,2)-lam
        axial=vp+(A-1)+Q(1,2)-lam
        angular=vp+(A-Q(3,2))+1-lam
        key=f"h{h}_lambda{lam}"
        checks[key+":three_physical_rows_have_zero_scale_power"]=pressure==axial==angular==0
        checks[key+":omitting_swirl_scale_detected"]=A-Q(1,2)-lam!=0
        checks[key+":common_P_Jz_scale_detected"]=vp+A+Q(1,2)-lam==1
        checks[key+":potential_has_extra_half_power"]=(A-Q(3,2))+Q(1,2)==A-1
    return dict(checks=checks,rows=rows,scope="exact signed matrix data and physical scaling exponents; not numerical selected profile evaluation")


def integral(poly,a,b,weight=0):
    out=Q(0)
    for (r,z,t,u),v in poly.terms.items():
        if z or t or u:raise ValueError("Radial fixture required")
        power=r+weight+1
        if not power:raise ValueError("Fixture integral would require a logarithm")
        out+=v*(b**power-a**power)/power
    return out


def bump_review():
    g=Q(3,2);v=1/g;l=Q(1);d=Q(5,4)
    f=R**4*(R-l)**6*(d-R)**6
    mu=tuple(integral(f,l,d,k) for k in (2,-3,-1,1,0))
    bp,bz,bg=responses(g,v,mu)
    bumps=[]
    for j in range(3):
        scale=g**j
        poly=means.Poly({key:value/scale**(key[0]+1) for key,value in f.terms.items()})
        bumps.append((scale*l,scale*d,poly))
    moment=lambda coeff,p:sum((c*integral(poly,a,b,p) for c,(a,b,poly) in zip(coeff,bumps)),Q(0))
    checks={"real_bump_angular_velocity_moment_zero":moment(bp,2)==moment(bz,2)==0,
            "real_bump_pressure_response":2*moment(bp,-3)==-1 and 2*moment(bz,-3)==0,
            "real_bump_axial_response":-moment(bp,-1)==0 and -moment(bz,-1)==-1,
            "real_bump_axial_velocity_moment_zero":moment(bg,1)==0,
            "real_bump_angular_response":moment(bg,0)==-1}
    ends=sorted({x for a,b,_ in bumps[:2] for x in (a,b)})
    carry=Q(0);previous=None;gap=False
    for i,(a,b) in enumerate(zip(ends,ends[1:])):
        mid=(a+b)/2;Fg=means.Poly()
        for c,(lo,hi,poly) in zip(bg,bumps):
            if lo<mid<hi:Fg=Fg+c*poly
        prim=(R*Fg).primitive();mass=prim-prim.at(0,a)+carry;psi=mass.shift(-1)
        checks[f"piece{i}_potential_radial_identity"]=psi.diff(0)+psi.shift(-1)==Fg
        if previous is not None:checks[f"piece{i}_potential_continuity"]=previous==psi.at(0,a)
        if not Fg:gap=gap or bool(psi.at(0,mid))
        if i==0:checks["potential_zero_at_left_hull"]=psi.at(0,a)==0
        previous=psi.at(0,b);carry=mass.at(0,b).scalar()
    checks["potential_zero_at_right_hull"]=previous==0
    checks["potential_fills_gap_between_bumps"]=gap
    checks["entire_axial_potential_hull_used"]=ends[0]==1 and ends[-1]==g*d and ends[1]<ends[2]
    return dict(checks=checks,scope="manufactured compact C5 piecewise polynomial geometric bumps at lambda=1/2; not the selected tilt or an all-orders profile")


def jet_review():
    checks={};m=4;f=(R-1)**8*(2-R)**8;x0=Q(3,2)
    for h,eta in ((Q(1,200),Q(0)),(Q(1,200),Q(1)),(Q(1,2**40),Q(-3,5))):
        t=fixed.geometry.implicit_jet(h,eta,m);log=fixed.log_one(t);A=Q(1,2)+h
        coord=Jet(m,{(0,0):eta,(1,0):1})*t.power_one(-Q(1,2)+h)
        a=2+coord*coord+coord*coord*coord*coord
        aprofile=2+R**2+R**4;derivative=aprofile;bounds=[]
        for k in range(1,m+1):
            derivative=derivative.diff(0)
            bounds.append(sum(abs(value) for value in derivative.terms.values())/factorial(k))
        positive=sum((abs(value) for key,value in a.c.items() if sum(key)),Q(0))
        checks[f"h{h}_eta{eta}:actual_background_composition_bound"]=positive<=background_positive_bound(m,bounds)
        for p,r in ((A,0),(A-1,1),(A-Q(3,2),2),(A-1,0)):
            fr=fixed.radial_diff(f,r);c=r-2*p;direct=Jet(m)
            for (k,z,tt,aux),value in fr.terms.items():direct=direct+value*x0**k*t.power_one(-Q(k,2)-c/2)
            expansion=Jet(m);power=Jet(m,{(0,0):1});current=fr;majorants=[]
            for k in range(m+1):
                majorants.append(fixed.profile_bound(current))
                expansion=expansion+Q((-1)**k,2**k*factorial(k))*current.at(0,x0).scalar()*power
                current=fixed.euler_shift(current,c);power=power*log
            key=f"h{h}_eta{eta}_p{p}_r{r}"
            checks[key+":fractional_generator_equals_direct_power"]=direct.c==expansion.c
            checks[key+":generalized_profile_majorant"]=sum(abs(v) for v in direct.c.values())<=composition_bound(p,r,m,majorants)
            checks[key+":wrong_unshifted_generator_detected"]=fixed.euler_shift(fr,c)!=fixed.euler_shift(fr,0)
    for m in (0,1,2,4,6):
        z=Jet(m,{(1,0):1});t=Jet(m,{(0,1):1});a=2+(10*z+5*t+z*z)*(10*z+5*t+z*z)
        inc=(a-2)*Q(1,2);inverse=Q(1,2)*inc.power_one(-1)
        checks[f"m{m}_reciprocal_identity"]=(a*inverse).c=={(0,0):Q(1)}
        positive=sum((abs(v) for key,v in a.c.items() if sum(key)),Q(0))
        checks[f"m{m}_finite_reciprocal_majorant"]=sum(abs(v) for v in inverse.c.values())<=reciprocal_bound(m,2,positive)
        if m>=2:checks[f"m{m}_no_high_norm_smallness_needed"]=positive>2
    return dict(checks=checks,scope="exact finite Taylor jets, including reciprocal states with large positive-order norms; no analytic profile assumption")


def polynomial_review():
    checks={};rows=[]
    Ap=P({1:3,2:1});PiP=P({0:2,1:1});PiZ=P({1:2});PiJ=P({0:3,2:1});PiJnext=P({0:5,3:2})
    for m in (0,1,3):
        H=reciprocal_bound(m,Q(1,2),Ap);Hnext=reciprocal_bound(m+1,Q(1,2),Ap)
        inputs=inverse_budget(m,H,Hnext,Q(2),Q(3),Q(5),Q(7),PiP,PiZ,PiJ,PiJnext)
        theta=coupled.slow_budget(m,coupled.geometry(1,2),P(1),P(2),P(3),inputs["X"],inputs["U"],inputs["Z"])
        for S in (1,4,16):
            ev=lambda p:p.evaluate(S)
            scalar=inverse_budget(m,ev(H),ev(Hnext),Q(2),Q(3),Q(5),Q(7),ev(PiP),ev(PiZ),ev(PiJ),ev(PiJnext))
            checks[f"m{m}_S{S}_explicit_inverse_polynomial_substitution"]=all(ev(inputs[k])==value for k,value in scalar.items())
            direct=coupled.slow_budget(m,coupled.geometry(1,2),1,2,3,scalar["X"],scalar["U"],scalar["Z"])
            checks[f"m{m}_S{S}_closed_coupled_scalar_polynomials"]=all(ev(theta[k])==value for k,value in direct.items())
            checks[f"m{m}_S{S}_reciprocal_polynomial_evaluation"]=ev(H)==reciprocal_bound(m,Q(1,2),ev(Ap))
        checks[f"m{m}_field_vs_target_jet_ledger"]=(m+2)+1==m+3 and ((m+1)+2)+1==m+4
        rows.append(dict(order=m,H=H.record(),U=inputs["U"].record(),Z=inputs["Z"].record(),X=inputs["X"].record()))
    return dict(checks=checks,rows=rows,scope="positive polynomial inverse/target coupling under retained scalar inputs, not verification of canonical seed/source norms")


def review():
    parts={"matrix_scaling":matrix_review(),"compact_profiles":bump_review(),"finite_jets":jet_review(),
           "coupled_substitution":polynomial_review(),"pressure_tails":means.pressure_tail_review(),
           "complete_covariance":means.covariance_review()}
    return dict(checks={p+":"+k:v for p,row in parts.items() for k,v in row["checks"].items()},parts=parts)
