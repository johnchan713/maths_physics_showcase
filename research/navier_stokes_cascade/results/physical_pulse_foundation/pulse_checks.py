"""Exact mixed physical pullbacks, moving frames and all-band controls.

Manufactured coefficients are not the giant actual background or a PDE solution.
"""
from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from itertools import product
from math import factorial
from pathlib import Path
import sys
from pulse_bounds import exponents,harmonic_loss,primary_powers,riccati_outward,tail_square_constant

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def maps():
    s=importlib.util.spec_from_file_location('physical_pulse_foundation_exact_maps',HERE.parent/'modulation_repair_all_order/map_checks.py')
    m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m);return m


def variable(i,value,n,order):return maps().JetClass().variable(i,value,n=n,order=order)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),maps().jet(0,a[0].n,a[0].order))
def matrix_product(a,b):return [[dot(row,col) for col in zip(*b)] for row in a]


@lru_cache(None)
def band_jets(center=Q(1,3),order=6):
    h=Q(1,200);Z=variable(0,center,2,order);T=variable(1,1-center**2,2,order)
    s=maps().jet(1,2,order)
    for _ in range(4):s=s-(s-Z*Z*s**(2*h)-T)/(1-2*h*Z*Z*s**(2*h-1))
    return s,s-Z*Z*s**(2*h)-T


def atan_zero(z):
    if z.value():raise ValueError('Zero-base angular jet required')
    return sum(((-1)**j*z**(2*j+1)/Q(2*j+1) for j in range((z.order+1)//2)),maps().jet(0,z.n,z.order))


def fields(R,Z,T,xi,v,theta,eps,k,L=Q(64)):
    m=maps();r0=Q(1,6057072);p=Q(1,k);pz=-Q(2,5)
    F=R**-3+Z/10+3*T*T/100;G=-R**-2+7*Z*T/100
    Phi=p*theta+pz*Z/eps+R/3-v*(p*F+pz*G)
    Phi=Phi-Phi.value() # frozen phase shift, preserving every derivative/carrier
    harmonic=m.cosine(k*Phi)
    Cr=1+R/5+Z/7+xi/(10*r0)+v/(11*L)+T/13
    Ct=2+R*Z/9+xi*v/(15*r0*L)+T/17
    Cz=-1+R*R/19+Z*T/23+xi/(12*r0)
    cost,sint=m.cosine(theta),m.sine(theta)
    potential=[(Cr*cost-Ct*sint)*harmonic/k,(Cr*sint+Ct*cost)*harmonic/k,Cz*harmonic/k]
    pressure=(1+R/7+Z/11+T/13+v/(17*L))*harmonic
    return potential,pressure,Phi,Ct*harmonic/k


def graph_operator(f,axis,eps,mi,dr=Q(2,3),drop_radial=False,freeze_radial=False,drop_fast=False):
    m=maps();o=4
    R=variable(0,1,6,o);theta=variable(5,0,6,o)
    radial=f.derivative(0)
    if not drop_radial:
        coefficient=m.jet(1,6,o) if freeze_radial else R**(dr-1)
        radial=radial+mi*dr*coefficient*f.derivative(3)
    if axis==0:return m.cosine(theta)*radial-m.sine(theta)/R*f.derivative(5)
    if axis==1:return m.sine(theta)*radial+m.cosine(theta)/R*f.derivative(5)
    if axis==2:return eps*f.derivative(1)
    if axis==3:return -eps*f.derivative(2)+(0 if drop_fast else f.derivative(4))
    raise ValueError('Unknown physical coordinate')


def graph_word(f,alpha,eps,mi,**options):
    for axis,n in enumerate(alpha):
        for _ in range(n):f=graph_operator(f,axis,eps,mi,**options)
    return f.value()


def curl(a,op):return [op(a[2],1)-op(a[1],2),op(a[0],2)-op(a[2],0),op(a[1],0)-op(a[0],1)]


@lru_cache(None)
def graph_stats():
    m=maps();count=0;matches=True;divergences=True;phase_matches=True;radial_gap=[];time_gap=[];frozen_gap=[];connection=[]
    for eps,k,mi in ((Q(1,4),2,Q(1,7)),(Q(1,8),3,Q(2,9)),(Q(1,16),4,Q(4,11))):
        R,Z,T,xi,v,theta=[variable(i,c,6,4) for i,c in enumerate((1,Q(1,7),Q(2,3),0,32,0))]
        a,p,phi,ct=fields(R,Z,T,xi,v,theta,eps,k)
        op=lambda f,axis:graph_operator(f,axis,eps,mi)
        w=curl(a,op)
        x,y,z,t=[variable(i,c,4,4) for i,c in enumerate((1,0,Q(1,7),Q(2,3)))]
        rd=(x*x+y*y)**Q(1,2);th=atan_zero(y/x)
        xid=mi*(rd**Q(2,3)-1);vd=32-(t-Q(2,3))/eps
        ad,pd,phid,ctd=fields(rd,z,t,xid,vd,th,eps,k)
        dop=lambda f,axis:eps*f.derivative(2) if axis==2 else f.derivative(axis)
        wd=curl(ad,dop)
        indices=[alpha for alpha in product(range(4),repeat=4) if sum(alpha)<=3]
        for alpha in indices:
            factor=eps**alpha[2]*(-eps)**alpha[3]
            for f,g in zip(w+[p],wd+[pd]):
                matches &= graph_word(f,alpha,eps,mi)==factor*g.partial(alpha);count+=1
            phase_matches &= graph_word(phi,alpha,eps,mi)==factor*phid.partial(alpha)
        div=wd[0].derivative(0)+wd[1].derivative(1)+eps*wd[2].derivative(2)
        divergences &= not div.terms
        radial_gap.append(graph_word(a[0],(1,0,0,0),eps,mi)-graph_word(a[0],(1,0,0,0),eps,mi,drop_radial=True))
        time_gap.append(graph_word(p,(0,0,0,1),eps,mi)-graph_word(p,(0,0,0,1),eps,mi,drop_fast=True))
        frozen_gap.append(graph_word(a[0],(2,0,0,0),eps,mi)-graph_word(a[0],(2,0,0,0),eps,mi,freeze_radial=True))
        connection.append(ct.value())
    return dict(count=count,matches=bool(matches),divergences=bool(divergences),phase_matches=bool(phase_matches),
        radial_graph_omission_nonzero=all(radial_gap),fast_time_omission_nonzero=all(time_gap),
        frozen_graph_coefficient_omission_nonzero=all(frozen_gap),cylindrical_connection_nonzero=all(connection))


def frame(center=Q(0),length=Q(1),order=9):
    m=maps();v=variable(0,0,1,order);s=center+v/length
    N=[Q(-3,5),Q(4,5)];K=[-N[1],N[0]]
    base=Q(1) if not center else Q(5,4)
    root=base*((1+s*s)/(1+center*center))**Q(1,2)
    n=[s,m.jet(K[0],1,order),m.jet(K[1],1,order)];norm=dot(n,n)
    U=[[m.jet(1,1,order),m.jet(0,1,order)],[-s*K[0],m.jet(N[0],1,order)],[-s*K[1],m.jet(N[1],1,order)]]
    M=[[m.jet(1,1,order),m.jet(1,1,order)],[-root,root]]
    B=matrix_product(U,M)
    Minv=[[m.jet(Q(1,2),1,order),-1/(2*base)*((1+s*s)/(1+center*center))**Q(-1,2)],
          [m.jet(Q(1,2),1,order),1/(2*base)*((1+s*s)/(1+center*center))**Q(-1,2)]]
    Ul=[[m.jet(1,1,order),m.jet(0,1,order),m.jet(0,1,order)],
        [m.jet(0,1,order),m.jet(N[0],1,order),m.jet(N[1],1,order)]]
    Bl=matrix_product(Minv,Ul)
    Kbg=[[m.jet(q,1,order) for q in row] for row in ((0,-2,0),(Q(14,25),0,0),(Q(48,25),0,0))]
    nK=[dot(n,col) for col in zip(*Kbg)];np=[z.derivative(0) for z in n]
    algebra=[[-Kbg[i][j]+n[i]*nK[j]/norm for j in range(3)] for i in range(3)]
    An=[[-Kbg[i][j]+n[i]*(nK[j]-np[j])/norm for j in range(3)] for i in range(3)]
    AB=matrix_product(An,B);Bp=[[z.derivative(0) for z in row] for row in B]
    moving=matrix_product(Bl,[[a-b for a,b in zip(row,der)] for row,der in zip(AB,Bp)])
    lam=Q(6,5)*root**-1
    target=[[lam,m.jet(0,1,order)],[m.jet(0,1,order),-lam]]
    E=[[a-b for a,b in zip(row,t)] for row,t in zip(moving,target)]
    return n,np,B,Bl,Kbg,algebra,An,E,target


@lru_cache(None)
def frame_stats():
    m=maps();checks=[];scaling=[];pressure=[];gaps=[];count=0
    for center in (Q(0),Q(3,4),Q(-3,4)):
        n,np,B,Bl,Kbg,algebra,An,E,target=frame(center)
        identities=[dot(n,col) for col in zip(*B)]
        identities += [a-Q(i==j) for i,row in enumerate(matrix_product(Bl,B)) for j,a in enumerate(row)]
        projected=matrix_product(Bl,matrix_product(algebra,B))
        identities += [a-b for row,t in zip(projected,target) for a,b in zip(row,t)]
        checks.extend(not z.terms for z in identities);count+=len(identities)*10
        for length in (Q(7),Q(32)):
            other=frame(center,length)[7]
            scaling.extend(other[i][j].partial((k,))*length**(k+1)==E[i][j].partial((k,))
                for i in range(2) for j in range(2) for k in range(9))
        t=[row[0] for row in B];f=[m.jet(q,1,9) for q in (Q(1,7),Q(2,11),Q(-3,13))]
        nt=dot(n,t);norm=dot(n,n);nf=dot(n,f);moving_term=dot(np,t);nKt=dot(n,[dot(row,t) for row in Kbg]);k=Q(3);harmonic=2
        P=(nKt-moving_term+nf)/(k*harmonic*norm)
        # pi=iP. The projected derivative is compatible even though this
        # manufactured t is not declared to solve the full pulse ODE.
        rhs=[dot(row,t)-harmonic**2*norm*z-(fi-ni*nf/norm) for row,z,fi,ni in zip(An,t,f,n)]
        pressure.append(not (dot(np,t)+dot(n,rhs)).terms)
        recovered=[-dot(row,t)-harmonic**2*norm*z+k*harmonic*ni*P-fi for row,z,ni,fi in zip(Kbg,t,n,f)]
        pressure.extend(not (a-b).terms for a,b in zip(recovered,rhs))
        gaps.append(moving_term.value()!=0 and nf.value()!=0)
    return dict(reference_components=count,reference=all(checks),moving_scaling=all(scaling),
        moving_components=len(scaling),normal_pressure=all(pressure),pressure_components=len(pressure),
        omitted_moving_or_source_normal_detected=all(gaps))


@lru_cache(None)
def scalar_stats():
    partitions=0;dominance=True;normalizations=True;barriers=[];squares=[]
    for h in (Q(1,200),Q(9,1000)):
        e=exponents(h)
        dominance &= e['sx']>e['coefficient_space'] and e['st']>e['coefficient_time']
        for a in range(7):
            for b in range(5):
                for j in range(1,a+b+1):
                    loss,upper=harmonic_loss(h,a,b,j);dominance &= loss<=upper;partitions+=1
                p=primary_powers(h,a,b);base=-a*e['sx']-b*e['st']
                normalizations &= p['velocity']==-e['A']+h/2+base and p['pressure']==-2*e['A']+h+base and p['potential']==Q(1,2)-e['A']+h+base
    for mu in (Q(1,10),Q(1),Q(10)):
        for ce in (Q(1,10),Q(1),Q(10)):
            for scale in (1,2,1000):
                r,out,upper=riccati_outward(mu,ce,8*ce*scale/mu);barriers.append(out<=upper<0 and r<=Q(1,8))
    for a in (Q(1,1000),Q(1,10),Q(2)):
        for B in (0,5,20):
            for loss in (Q(-100),Q(0),Q(20)):
                for K in (0,1,100):
                    d,c=tail_square_constant(a,B,loss,K)
                    for ell in (Q(1),Q(16),Q(1000000)):
                        squares.append(c-(d*ell-a*ell*ell)==a*(ell-d/(2*a))**2 and d*ell-a*ell*ell<=c)
    return dict(partitions=partitions,dominance=bool(dominance),normalizations=bool(normalizations),
        barrier_components=len(barriers),barriers=all(barriers),square_components=len(squares),squares=all(squares))


@lru_cache(None)
def diagnostic_review():
    g,f,s=graph_stats(),frame_stats(),scalar_stats();bands=[band_jets(e) for e in (Q(-1),Q(1,3),Q(1))]
    return dict(exact_physical_field_components=g['count'],reference_frame_components=f['reference_components'],
        moving_frame_scaling_components=f['moving_components'],pressure_constraint_components=f['pressure_components'],
        implicit_band_components=3*28,phase_exponent_partitions=s['partitions'],
        Riccati_barrier_components=s['barrier_components'],all_band_square_components=s['square_components'],
        checks=dict(all_mixed_direct_physical_and_graph_curl_pressure_jets_match=g['matches'],
            direct_complete_curl_is_exactly_divergence_free=g['divergences'],
            all_phase_mixed_jets_have_exact_radial_graph_cancellation=g['phase_matches'],
            amplitude_radial_graph_fast_time_and_differentiated_coefficient_omissions_detected=g['radial_graph_omission_nonzero'] and g['fast_time_omission_nonzero'] and g['frozen_graph_coefficient_omission_nonzero'],
            cylindrical_connection_term_nonzero=g['cylindrical_connection_nonzero'],
            reference_frame_projection_and_left_inverse_jets_exact=f['reference'],
            moving_frame_error_has_one_over_L_scaling_in_every_checked_jet=f['moving_scaling'],
            complete_normal_pressure_retains_moving_and_source_terms=f['normal_pressure'] and f['omitted_moving_or_source_normal_detected'],
            implicit_band_coordinate_identity_through_order_six=all(not r.terms for _,r in bands),
            phase_partitions_and_primary_physical_normalizations_match=s['dominance'] and s['normalizations'],
            accepted_Riccati_outward_bound_is_strict=s['barriers'],
            Gaussian_square_completion_holds_on_whole_same_band_domain=s['squares'],
            all_fields_and_band_parameters_are_labelled_manufactured_controls=True))
