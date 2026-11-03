"""Current-source EP covariant inverse and unchanged full restoration, conditional."""
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
    s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
mixed=load('q117_mixed','core_current_head_q116_mixed_covariant_jets/q116.py')
current=load('q117_current','core_current_head_q110_radial_tail_request/q110.py')
covariance=current.radial.pressure.covariance
native=covariance.response
EP,cv,prod=native.EP,native.cv,native.prod

def damping_lower(q0,supplied=None):
    if isinstance(q0,(float,bool)) or Q(q0)<=0:raise ValueError('Exact original normal-square reciprocal required')
    derived=1/Q(q0)
    if supplied is None:return derived
    if isinstance(supplied,(float,bool)) or not 0<Q(supplied)<=derived:
        raise ValueError('A supplied floor may only weaken the same-normal derived floor 1/q0')
    return Q(supplied)

def response(G,n,np,K,rate,q0,tangent,length,*,damping_floor=None,entrance_flat_through,harmonic=1):
    arrays=[G,n,np,K,rate]
    if not G or any(len(a)!=len(G) for a in arrays):raise ValueError('Common nonempty raw mixed-jet ceiling required')
    m=mixed.first.harmonic(harmonic);N=len(G)-1
    if type(entrance_flat_through) is not int or entrance_flat_through<max(-1,N-1):
        raise ValueError('Actual whole-source entrance flatness through N-1 must be supplied as a contract')
    dmin=damping_lower(q0,damping_floor);q0=Q(q0)
    G,n,np,K,rate=[list(map(cv,a)) for a in arrays]
    T=1+cv(tangent);L=cv(length)
    if not L.c:raise ValueError('Original positive length bound required')
    square=prod(n,n);inv=[EP(q0)]
    for j in range(1,N+1):
        inv.append(q0*sum((comb(j,i)*square[i]*inv[j-i] for i in range(1,j+1)),EP()))
    h=prod(n,inv);Pi=[EP(1)]+prod(h,n)[1:]
    A=[a+b for a,b in zip(prod(Pi,K),prod(h,np))]
    Qlift=q0*(K[0]*n[0]+np[0]);C=L+Q(4,3)/dmin
    F=prod(Pi,G);d=[4*x for x in square];V=[]
    for j in range(N+1):
        z=sum((comb(j,i)*n[i]*V[j-i] for i in range(1,j+1)),EP())
        lower=sum((comb(j,i)*(A[i]+rate[i])*V[j-i] for i in range(1,j+1)),EP())
        damp=sum((comb(j,i)*d[i]*V[j-i] for i in range(1,j+1)),EP())
        V.append(T*(F[j]+L*Qlift*z+L*lower+C*damp)+h[0]*z)
    # This upper bound avoids a min of polynomials; the true residence still
    # satisfies ell<=L in the proof of V. No comparison of EP values is assumed.
    ell=L if m==1 else EP(1/((m*m-1)*dmin))
    U=[ell*v for v in V]
    row=[a+b for a,b in zip(prod(n,K),np)]
    numerator=[a+b for a,b in zip(prod(row,U),prod(n,G))]
    pressure=[Q(1,m)*x for x in prod(inv,numerator)]
    potential=[Q(1,m)*x for x in prod(h,U)]
    return dict(velocity=U,pressure=pressure,potential=potential,
                velocity_over_residence=V,residence_upper=ell,
                damping_floor=dmin,normal_reciprocal_derived_floor=1/q0,
                source_linear=True,new_source_pole_added=0,
                actual_native_and_entrance_contracts_certified=False)

class CovariantResponse(native.SourcedResponse):
    """Same Q95 source and Q96 physical solution; an alternative norm proof."""
    def build(self,j,linear_table,log_table,rate_table,cutoff,*,tangent,length,
              damping_floor=None,entrance_flat_through,log_poles=None,harmonic=1):
        native.full.order(j);R=j+1
        if self.inputs.state is not self.state or self.inputs.source.state is not self.state or self.inputs.source.waves.state is not self.state:
            raise ValueError('One identical current whole-path source state required')
        if len(rate_table)!=R+1 or len(cutoff)!=R+1:raise ValueError('Complete curl retains one extra response/source jet')
        # Reject the new missing contracts before evaluating a large source graph.
        if type(entrance_flat_through) is not int or entrance_flat_through<j:
            raise ValueError('Whole mixed entrance source jets through output order required')
        damping_floor=damping_lower(self.inputs.source.waves.q,damping_floor)
        data=self.inputs.build(R,linear_table,log_table,log_poles)
        coeff=response([2*x for x in data['normalized']],
            [self.phase('normal',k) for k in range(R+1)],
            [self.phase('moving_normal',k) for k in range(R+1)],
            [self.phase('matrix',k) for k in range(R+1)],rate_table,
            self.inputs.source.waves.q,tangent,length,
            damping_floor=damping_floor,entrance_flat_through=entrance_flat_through,harmonic=harmonic)
        bell=native.inputs.quotient([0]*(R+1),log_table,[0]*(R+1) if log_poles is None else log_poles)
        budget=self.inputs.source.means.current.initial.budget
        compact=budget.compact(0,budget.pole)
        E=1+376*EP({1:1})*cv(compact['I']);powers=[EP(1)]
        for k in range(R):powers.append(powers[-1]*E)
        restored={name:[powers[k]*x for k,x in enumerate(prod(bell['Bell'],coeff[name]))]
                  for name in ('velocity','pressure','potential')}
        cut={name:prod(cutoff,v) for name,v in restored.items()}
        a=budget._data(R,budget.pole);g=cv(a['graph']);r=cv(a['r']);curl=[]
        for k in range(j+1):
            dr=sum((comb(k,i)*factorial(i)*g*cut['potential'][k-i+1] for i in range(k+1)),EP())
            radial=sum((comb(k,i)*factorial(i)*r*cut['potential'][k-i] for i in range(k+1)),EP())
            curl.append(2*dr+2*cut['potential'][k+1]+radial)
        return dict(principal=cut['velocity'][:j+1],curl=curl,
                    complete=[cut['velocity'][k]+curl[k] for k in range(j+1)],
                    pressure=cut['pressure'][:j+1],potential=cut['potential'][:j+1],
                    normalized_response=coeff,evaluation_factor=E,
                    pole=data['common_pole']+bell['extra_pole'],source_ceiling=R,
                    source_factory_reserve=R+2,source_state=self.state,
                    entrance_flatness_contract_through=entrance_flat_through,
                    unweighted_operator_coefficient_contract_required=True,
                    actual_contracts_certified=False,changed_covariance_or_residual_improvement=False)

def record():
    return dict(status='current-source-covariant-response-and-full-restoration',
                source_method='Q95 PathInputs.build on the identical current source',
                covariant_response_source_linear=True,new_source_pole_added=0,
                damping_floor_derived_from_same_normal_reciprocal=True,
                original_Gaussian_Bell_packet_cutoff_and_complete_curl_retained=True,
                Q97_two_carrier_covariance_type_compatible=True,
                actual_native_entrance_or_high_pole_contracts_certified=False,
                full_cone_pulse_residual_or_blowup=False)
