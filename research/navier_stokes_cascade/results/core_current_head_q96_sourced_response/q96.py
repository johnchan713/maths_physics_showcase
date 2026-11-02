"""Conditional normalized forced-response and complete cut-curl jet bounds."""
from pathlib import Path
from math import comb,factorial
from fractions import Fraction as Q
import importlib.util
p=Path(__file__).resolve().parent.parent/'core_current_head_q95_path_input_replay/q95.py'
s=importlib.util.spec_from_file_location('q96_inputs',p);inputs=importlib.util.module_from_spec(s);s.loader.exec_module(inputs)
full=inputs.full;EP=inputs.EP;cv=inputs.convert;prod=inputs.product

def response(G,n,np,K,rate,q0,tangent,length,harmonic=1):
    arrays=[G,n,np,K,rate]
    if len(G)<2 or any(len(a)!=len(G) for a in arrays):raise ValueError('Common jet ceiling at least one required')
    if type(harmonic) is not int or harmonic==0:raise ValueError('Fixed nonzero integer harmonic required')
    q0=Q(q0)
    if q0<=0:raise ValueError('Positive normal-square reciprocal bound required')
    G,n,np,K,rate=[list(map(cv,a)) for a in arrays];T=1+cv(tangent);L=cv(length);m=abs(harmonic)
    D=prod(n,n);inv=[EP(q0)]
    for j in range(1,len(G)):
        inv.append(q0*sum((comb(j,i)*D[i]*inv[j-i] for i in range(1,j+1)),EP()))
    h=prod(n,inv);v=[x+y for x,y in zip(prod(n,K),np)]
    A0=[x+y for x,y in zip(K,prod(h,v))];H=prod(h,n)
    damp=[4*m*m*x for x in D]
    extra=prod([x+y for x,y in zip(damp,rate)],H)
    # Atilde=Ahat-pI, not Ahat: the final rate term is essential.
    C=[a+b+c+p for a,b,c,p in zip(A0,damp,extra,rate)]
    projection=[EP(1)]+H[1:];F=prod(projection,G)
    ambient=T*(1+L*(A0[0]*h[0]+h[1])*n[0])+h[0]*n[0]
    U=[]
    for j in range(len(G)):
        parameter=ambient*L*(F[j]+sum((comb(j,i)*C[i]*U[j-i] for i in range(1,j+1)),EP()))
        pulse=EP() if j==0 else F[j-1]+sum((comb(j-1,i)*C[i]*U[j-1-i] for i in range(j)),EP())
        U.append(parameter+pulse)
    pressure=[Q(1,m)*(a+b) for a,b in zip(prod(prod(inv,v),U),prod(h,G))]
    potential=[Q(1,m)*a for a in prod(h,U)]
    return dict(velocity=U,pressure=pressure,potential=potential,generator=C,ambient=ambient)

class SourcedResponse:
    def __init__(self,path_inputs):
        if not isinstance(path_inputs,inputs.PathInputs):raise TypeError('Same-state Q95 path inputs required')
        self.inputs=path_inputs;self.state=path_inputs.state
    def phase(self,name,j):
        w=self.inputs.source.waves
        norm=sum((Q(5**k,factorial(k))*w.jet(name,k) for k in range(j+1)),EP())
        factor=2*cv(self.state.primary.chart.transfer(j))*16*188**j*factorial(j)
        return 2**j*factorial(j)*factor*inputs.source.waves.rescale(norm)
    def build(self,j,linear_table,log_table,rate_table,cutoff,*,tangent,length,log_poles=None,harmonic=1):
        full.order(j);R=j+1
        if len(rate_table)!=R+1 or len(cutoff)!=R+1:raise ValueError('One extra response jet needed for the complete curl')
        data=self.inputs.build(R,linear_table,log_table,log_poles)
        coeff=response([2*x for x in data['normalized']],
            [self.phase('normal',k) for k in range(R+1)],
            [self.phase('moving_normal',k) for k in range(R+1)],
            [self.phase('matrix',k) for k in range(R+1)],rate_table,
            self.inputs.source.waves.q,tangent,length,harmonic)
        bell=inputs.quotient([0]*(R+1),log_table,[0]*(R+1) if log_poles is None else log_poles)
        budget=self.inputs.source.means.current.initial.budget
        compact=budget.compact(0,budget.pole)
        E=1+376*EP({1:1})*cv(compact['I'])
        powers=[EP(1)]
        for k in range(R):powers.append(powers[-1]*E)
        restored={name:[powers[k]*x for k,x in enumerate(prod(bell['Bell'],coeff[name]))]
                  for name in ('velocity','pressure','potential')}
        cut={name:prod(cutoff,v) for name,v in restored.items()}
        a=budget._data(R,budget.pole);g=cv(a['graph']);r=cv(a['r'])
        curl=[]
        for k in range(j+1):
            dr=sum((comb(k,i)*factorial(i)*g*cut['potential'][k-i+1] for i in range(k+1)),EP())
            radial=sum((comb(k,i)*factorial(i)*r*cut['potential'][k-i] for i in range(k+1)),EP())
            curl.append(2*dr+2*cut['potential'][k+1]+radial)
        return dict(principal=cut['velocity'][:j+1],curl=curl,
                    complete=[cut['velocity'][k]+curl[k] for k in range(j+1)],
                    pressure=cut['pressure'][:j+1],potential=cut['potential'][:j+1],
                    normalized_response=coeff,evaluation_factor=E,
                    pole=data['common_pole']+bell['extra_pole'],source_ceiling=R,
                    source_factory_reserve=R+2,actual_contracts_certified=False,
                    changed_covariance_or_residual_improvement=False)

def record():
    return dict(status='conditional-same-state-zero-entry-response-and-complete-cut-curl',
        response_order='1-kappa',curl_order='3/2-2*kappa',pressure_order='3/2-kappa',
        normalized_generator_includes_minus_profile_rate=True,
        packet_evaluation_factor='1+376*S*I0',source_ceiling_shift=1,
        actual_contracts_certified=False,changed_covariance_or_residual_improvement=False,
        unconditional_pulse_or_blowup=False)
