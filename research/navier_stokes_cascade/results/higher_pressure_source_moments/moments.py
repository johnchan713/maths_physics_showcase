"""Full ordered pressure-source moments; conditional estimates, not PDE solutions."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m
old=load('higher_source_exact_calculus',HERE.parent/'signed_pressure_source_moments/moments.py')
J,P,X,rat=old.J,old.P,old.X,old.rat

def order(k):
    if type(k) is not int or k<0:raise ValueError('Nonnegative integer coefficient order required')
    return k

def fixture_terms(j,index):
    """Two compact C4 endpoint fixtures, one with positive-order zero axis traces."""
    order(j)
    if index not in (0,1):raise ValueError('Fixture index required')
    power=1+int(j>0 and index==0);out={}
    p={0:F(j+1,3),1:F(2*j+1,5),2:F(1,7),4:F(j+2,11)}
    q={0:F(1,13),1:F(j+1,17),3:F(1,19)}
    for n in range(7):
        wall=F((-1)**n*comb(6,n))
        for degree,c in p.items():out[n+power,degree]=out.get((n+power,degree),F(0))+wall*c
        for degree,c in q.items():out[n+power+1,degree]=out.get((n+power+1,degree),F(0))+wall*c
    return out

def fields(k,index,angle,h):
    e=J((angle,1));D=F(1,2)-h;L=1-2*h*e*e;d=1-e*e;out=[]
    for j in range(k+1):
        M=P.cast(0)
        for (i,a),c in fixture_terms(j,index).items():M+=P({i:c*e**a})
        U=M.dx();V=(2*e*X*U-2*(D+2*h*j)*e*M-d*M.de())/L
        out.append((U,V,V.over_X()))
    return out,e,L,d

def Z(a,f,e,L,d):return (2*a*e*f-2*e*X*f.dx()+d*f.de())/L
def T(b,f,e,L,d):return (2*b*e*f+d*f.de())/L
def zero(f):return all(not any(v.c[:2]) for v in f.c.values()) if isinstance(f,P) else not any(f.c[:2])

def identities(k,index,angle,h,nu,R):
    """Original hierarchy density vs independent flux and integrated identities."""
    order(k);angle,h,nu,R=map(rat,(angle,h,nu,R));fs,e,L,d=fields(k,index,angle,h)
    A=F(1,2)+h;D=F(1,2)-h;lam=2*h*k;U,V,v=fs[k]
    w=((1-lam)*v+X*v.dx()+D*e*v.de())/L-2*nu*V.dx().dx()
    B=P.cast(0);C=P.cast(0);UV=P.cast(0)
    for i in range(k+1):
        j=k-i;ui,vi=fs[i][0],fs[i][2];vj=fs[j][2]
        w+=vi*(vj/2+X*vj.dx())+ui*(2*(2*h*j-1)*e*vj-2*e*X*vj.dx()+d*vj.de())/L
        B+=ui*vj;C+=vi*vj;UV+=ui*fs[j][1]
    shift=P.cast(0)
    if k:shift=Z(2*h*(k-1)-1-D,Z(2*h*(k-1)-1,fs[k-1][2],e,L,d),e,L,d)
    w-=nu*shift
    flux=(V-2*e*UV)/L+X*C-2*nu*V.dx()
    conservative=flux.dx()+(-lam*v+D*e*v.de()+d*B.de()+2*(lam-A)*e*B)/L+C/2-nu*shift
    axis=2*nu*v.value(0);direct=[];reduced=[];boundaries=[];axial=[];axial_boundary=[]
    for weight in (0,1):
        W=P.cast(1) if weight==0 else X
        a=(W*v).total(R);b=(W*B).total(R);c=(W*C).total(R)
        if weight==0:
            base=axis+(-lam*a+D*e*a.de()+d*b.de()+2*(lam-A)*e*b)/L+c/2
            edge=flux.value(R)
        else:
            base=(-(1+lam)*a+D*e*a.de()+d*b.de()+2*(lam+D)*e*b)/L-c/2
            edge=R*flux.value(R)+2*nu*V.value(R)
        reduced_shift=J.cast(0);shift_edge=J.cast(0)
        if k:
            f=fs[k-1][2];a0=(W*f).total(R);beta=2*h*(k-1)+weight
            first=Z(2*h*(k-1)-1,f,e,L,d)
            reduced_shift=T(beta-D,T(beta,a0,e,L,d),e,L,d)
            first_edge=-2*e*R**(weight+1)*f.value(R)/L
            shift_edge=T(beta-D,first_edge,e,L,d)-2*e*R**(weight+1)*first.value(R)/L
        direct.append((W*w).total(R));reduced.append(base-nu*reduced_shift)
        boundaries.append(edge-nu*shift_edge);axial.append((W*shift).total(R));axial_boundary.append(shift_edge)
    return dict(pointwise=zero(w-conservative),
        direct=direct,reduced=reduced,boundary=boundaries,axis=axis,axial=axial,axial_boundary=axial_boundary,
        exact=all(zero(a-b-c) for a,b,c in zip(direct,reduced,boundaries)),
        full_endpoint_zero=all(zero(z) for z in boundaries),
        shifted_omission_detected=any(not zero(nu*a) for a in axial),
        truncated_omission_detected=any(not zero(z) for z in boundaries),
        positive_order_axis_zero=(k==0 or index!=0 or zero(axis)))

def source_bounds(k,h,nu,Xv,u2,v2,previous_v3=None):
    """Signed C_eta^1 totals only; every array entry is a full-support norm."""
    order(k)
    if k and previous_v3 is None:raise ValueError('Previous full-support third angular norm must be supplied')
    if previous_v3 is None:previous_v3=0
    h,nu,Xv,previous_v3=map(rat,(h,nu,Xv,previous_v3))
    u2,v2=tuple(map(rat,u2)),tuple(map(rat,v2))
    if not 0<h<=F(1,100) or nu<0 or Xv<=0 or len(u2)!=k+1 or len(v2)!=k+1 or min((*u2,*v2,previous_v3))<0:
        raise ValueError('Admissible parameters and all nonnegative coefficient envelopes required')
    A=F(1,2)+h;D=F(1,2)-h;lam=2*h*k;l=1-2*h
    R1=1/l+4*h/l**2;R2=1/l+6*h/l**2+16*h*h/l**3
    b=sum((u2[i]*v2[k-i] for i in range(k+1)),F(0))
    c=sum((v2[i]*v2[k-i] for i in range(k+1)),F(0))
    shifted=[];q=[]
    for j in (0,1):
        volume=Xv**(j+1)/F(j+1);beta=2*h*(k-1)+j
        sh=nu*R1*R2*(6+4*abs(beta-D))*(12+4*abs(beta))*volume*previous_v3 if k else F(0)
        if j==0:
            z=R1*volume*((lam+4*D)*v2[k]+(6+4*abs(lam-A))*b)+volume*c/2+2*nu*v2[k]+sh
        else:
            z=R1*volume*((1+lam+4*D)*v2[k]+(6+4*abs(lam+D))*b)+volume*c/2+sh
        shifted.append(sh);q.append(z)
    return dict(I0=q[0],I1=q[1],axial=tuple(shifted),B2=b,C2=c,R1=R1,R2=R2)

def second_repair_bounds(xcut,xvel,J2,f0,u0,f1,u1,Q0,Q1):
    """Five ordinary n=2 discrepancies, before unchanged linear inverse 800."""
    xcut,xvel,J2,f0,u0,f1,u1,Q0,Q1=map(rat,(xcut,xvel,J2,f0,u0,f1,u1,Q0,Q1))
    if min(xcut,xvel)<=0 or min(J2,f0,u0,f1,u1,Q0,Q1)<0:raise ValueError('Positive supports and nonnegative norms required')
    return (xcut*J2,xcut*xcut*J2,2*f0*xcut*J2+xvel*f1*f1+Q0/2,
        xcut*xcut*(u0+f0)*J2+xvel*xvel*u1*f1,
        2*u0*xcut*J2+f0*xcut*xcut*J2+xvel*u1*u1+xvel*xvel*f1*f1/2+Q1/2)

def regular_v_bound(k,h,angular_order,U_next):
    """Only orders 2,3 needed here; global U requires one EXTRA angular jet."""
    order(k);h,U_next=map(rat,(h,U_next))
    if not 0<h<=F(1,100) or type(angular_order) is not int or angular_order not in (2,3) or U_next<0:raise ValueError('Valid h, angular order 2 or 3, and norm required')
    # R3 follows by differentiating L^-1: coefficient sums in h are nonnegative.
    l=1-2*h;R2=1/l+6*h/l**2+16*h*h/l**3
    R3=R2+16*h*h/l**3+64*h**3/l**4
    R=(R2 if angular_order==2 else R3);D=F(1,2)-h
    return R*(4+4*(D+2*h*k)+4*(angular_order+1))*U_next

@lru_cache(None)
def exact_review():
    rows=[]
    for k,index,e,h,R in product((0,1,3), (0,1), (F(-1),F(0),F(3,5),F(1)),(F(1,200),),(F(1,3),F(1))):
        z=identities(k,index,e,h,F(1),R)
        rows.append(dict(k=k,index=index,eta=str(e),R=str(R),
            direct=[[str(c) for c in j.pair()] for j in z['direct']],
            axial=[[str(c) for c in j.pair()] for j in z['axial']],
            checks=dict(ordered_pointwise_flux=z['pointwise'],finite_endpoint_moment_identity=z['exact'],
                full_endpoint_vanishes=(R!=1 or z['full_endpoint_zero']),positive_order_zero_axis_trace=z['positive_order_axis_zero']),
            shifted_omission_detected=z['shifted_omission_detected'],truncation_omission_detected=z['truncated_omission_detected']))
    return dict(rows=rows,scalar_integrated_jet_identities=4*len(rows),manufactured_C4_endpoint_fixtures_not_actual_PDE=True,
        checks=dict(all_identities=all(all(r['checks'].values()) for r in rows),
            order_two_pressure_shift_is_nonzero=any(r['k']==1 and r['R']=='1' and r['shifted_omission_detected'] for r in rows),
            future_support_cannot_be_truncated=any(r['R']=='1/3' and r['truncation_omission_detected'] for r in rows),
            angular_endpoints_included=any(r['eta']=='-1' for r in rows) and any(r['eta']=='1' for r in rows)))

def legacy_review():
    """Independent original floating hierarchy comparison, not interval certification."""
    path=HERE.parent/'background_summation';prior=sys.path[:]
    try:
        sys.path.insert(0,str(path));m=load('higher_source_original_hierarchy',path/'hierarchy.py')
        from residual import Poly
        errors=[];omissions=[]
        for k in (0,1,3):
            profiles=[]
            for j in range(k+1):
                t=fixture_terms(j,0)
                profiles.append(SimpleNamespace(F=Poly({}),Pi=Poly({}),U=Poly({(i-1,a):float(i*c) for (i,a),c in t.items()}),
                    average_U=Poly({(i-1,a):float(c) for (i,a),c in t.items()})))
            for x,e in product((F(0),F(1,3),F(4,5)),(F(-1),F(0),F(1))):
                fs,eta,L,d=fields(k,0,e,F(1,200));A=F(1,2)+F(1,200);D=1-A;lam=F(k,100);U,V,v=fs[k]
                # Recover density independently as derivative of integrated flux identity.
                B=sum((fs[i][0]*fs[k-i][2] for i in range(k+1)),P.cast(0))
                C=sum((fs[i][2]*fs[k-i][2] for i in range(k+1)),P.cast(0))
                flux=(V-2*eta*X*B)/L+X*C-2*V.dx()
                density=flux.dx()+(-lam*v+D*eta*v.de()+d*B.de()+2*(lam-A)*eta*B)/L+C/2
                if k:density-=Z(F(k-1,100)-1-D,Z(F(k-1,100)-1,fs[k-1][2],eta,L,d),eta,L,d)
                actual=m.lower_sources(profiles,float(x),float(e),.005,1.)
                val=float(density.value(x).c[0]);errors.append(abs(val-actual['omega_over_X'])/max(1.,abs(val)))
                if k:omissions.append(abs(actual['omega_over_X']-m.lower_sources(profiles,float(x),float(e),.005,1.,True)['omega_over_X']))
        return dict(samples=len(errors),relative_tolerance='2e-9',not_an_interval_certificate=True,
            checks=dict(original_hierarchy_density_matches=bool(max(errors)<2e-9),original_axial_viscosity_omission_detected=bool(max(omissions)>1e-6)))
    finally:sys.path[:]=prior

def norm_review():
    h=F(1,200);q=source_bounds(1,h,1,7,(2,3),(5,11),17);z=source_bounds(1,h,0,7,(2,3),(5,11),17)
    oldq=old.source_bounds(h,1,7,2,5);newq=source_bounds(0,h,1,7,(2,),(5,))
    repair=second_repair_bounds(2,3,5,7,11,13,17,q['I0'],q['I1'])
    support=load('higher_source_unchanged_inverse',HERE.parent/'coefficient_support_foundation/support_bounds.py')
    rows,cap=support.target_bound(1,F(7),F(2),F(3,5),F(1,1000),repair)
    l=F(49,50);R3=1/l+6*F(1,100)/l**2+32*F(1,100)**2/l**3+64*F(1,100)**3/l**4
    return dict(example={k:([str(x) for x in v] if isinstance(v,tuple) else str(v)) for k,v in q.items()},
        second_repair=list(map(str,repair)),normalized_rows=list(map(str,rows)),coefficient_bound=str(cap),
        checks=dict(k0_exactly_recovers_frozen_signed_bound=all(newq[k]==oldq[k] for k in ('I0','I1')),
            order_one_source_has_both_ordered_cross_terms=q['B2']==37 and q['C2']==110,
            weighted_axial_viscosity_not_radial_boundary_cancellation=q['I1']>z['I1'],
            unweighted_axis_and_axial_terms_both_charged=q['I0']-z['I0']==22+q['axial'][0],
            weighted_difference_exactly_axial=q['I1']-z['I1']==q['axial'][1],
            global_v3_from_U4_retains_extra_angular_jet=regular_v_bound(0,h,3,1)>regular_v_bound(0,h,2,1),
            reciprocal_third_angular_bound_less_than_two=R3<2,
            unchanged_linear_inverse_800_with_amplitude_and_slope_cost=cap==800*max(rows),
            lower_F1_square_retained_in_Cp=repair[2]>2*7*2*5+q['I0']/2,
            lower_U1_square_retained_in_S=repair[4]>2*11*2*5+7*4*5+q['I1']/2))

def fixture_bound_review():
    rows=[];h=F(1,200)
    for k,index in product((1,3),(0,1)):
        U4=[sum((i*abs(c)*sum(F(comb(a,r)) for r in range(min(4,a)+1))
            for (i,a),c in fixture_terms(j,index).items()),F(0)) for j in range(k+1)]
        v2=[regular_v_bound(j,h,2,U4[j]) for j in range(k+1)]
        q=source_bounds(k,h,1,1,U4,v2,regular_v_bound(k-1,h,3,U4[k-1]))
        checks=[]
        for e in (F(-1),F(0),F(3,5),F(1)):
            z=identities(k,index,e,h,F(1),F(1))
            checks.extend(sum(map(abs,z['direct'][j].pair()))<=q['I'+str(j)] for j in (0,1))
        rows.append(dict(k=k,index=index,global_U4_polynomial_triangle=list(map(str,U4)),
            signed_caps=[str(q['I0']),str(q['I1'])],checks=dict(exact_sample_moment_jets_obey_global_analytic_caps=all(checks))))
    return dict(rows=rows,checks=dict(all_fixture_caps=all(all(r['checks'].values()) for r in rows)))

def review():
    a,b,c,d=exact_review(),legacy_review(),norm_review(),fixture_bound_review()
    return dict(exact=a,original_hierarchy=b,norms=c,fixture_caps=d,
        checks=dict(**a['checks'],**b['checks'],**c['checks'],**d['checks']))
