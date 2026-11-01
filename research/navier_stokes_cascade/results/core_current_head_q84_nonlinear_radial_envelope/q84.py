"""Positive rational evaluations of the final nonlinear radial norm recipe."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
_path=Path(__file__).resolve().parent.parent/'core_angular_pressure_normal_form/bounds.py'
_spec=importlib.util.spec_from_file_location('q84_angular',_path)
angular=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(angular)

def envelope(m,data):
    if type(m) is not int or m<0:raise ValueError('Nonnegative output order required')
    keys=('b0','V','G','x','u','z','W','Z','e','r','g','c','T')
    if set(data)!=set(keys):raise ValueError('Complete stripped FINAL-state norm data required')
    if any(isinstance(v,(bool,float)) or not isinstance(v,(int,Q)) or v<0 for v in data.values()):
        raise ValueError('Exact nonnegative norm inputs required')
    b,V,G,x,u,z,W,Z,e,r,g,c,T=(Q(data[k]) for k in keys)
    d=m+7;D=3*d*(1+g);a=2*V*r;a1=6*V*r*r+2*d*g*V*r
    vis=D*D+r*D+r*r+d*d
    J=b*u+V*x+x*u+Z;L=G*u+V*z+z*u+W
    F=e+d*u+d*L+vis*u;H=b*z+G*x+x*z+W
    rest=d*(1+4*c)*x+r*(2*b*x+x*x)+d*H+2*r*Z+d*d*x
    result=2*(a1*T*J+a*T*F+r*u*u+rest)
    return dict(Jd=J,F=F,rest=rest,projected_remainder=result,
                field_order=m+6,angular_row_order=m+4,temporal_source_order=m+4)

def rates(k=Q(1,100000)):
    B=1-k
    return dict(Jd=min(1+B,2-k,3-2*k,2-2*k),
                F=min(2-2*k,1+B,1+B,1+B-2*k,B+3),
                nonF=min(2*B,2-k,3-k,4-2*k,2,2-2*k,5-k))

def record():
    a={k:Q(i+1) for i,k in enumerate(('b0','V','G','x','u','z','W','Z','e','r','g','c','T'))}
    out=envelope(0,a)
    return dict(status='explicit-final-nonlinear-radial-remainder-envelope',
                temporal_inverse_recipe=angular.temporal_recipe(199996),
                fixture={k:str(v) for k,v in out.items()},rates={k:str(v) for k,v in rates().items()},
                radial_output_extra=199996,field_extra=200002,angular_row_extra=200000,
                actual_norm_inequalities_certified=False,pressure_removed=False,
                unconditional_pulse_residual_or_blowup=False)
