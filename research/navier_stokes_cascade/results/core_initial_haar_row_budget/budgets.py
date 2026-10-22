"""Conditional initial Haar-row norm budgets with complete pressure and means."""
from fractions import Fraction as Q
from math import ceil,comb,factorial
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("haar_row_edge",HERE.parent/"core_signed_target_edge_budget/budgets.py")
edge=importlib.util.module_from_spec(spec);spec.loader.exec_module(edge)
root=edge.root;atlas=edge.atlas;P=root.P;order=atlas.order;exact=atlas.exact
KEYS=("b","V","G","Ft","Fr","U","Zg","X","W","D","PiP","Rho","Kp","r")


def row_budget(n,data):
    order(n)
    if set(data)!=set(KEYS):raise ValueError("All base, complete mean/covariance and pressure inputs required")
    if all(isinstance(v,P) for v in data.values()):pass
    elif all(type(v) in (int,Q) for v in data.values()):
        data={k:exact(v) for k,v in data.items()}
        if min(data.values())<0:raise ValueError("Nonnegative sufficient bounds required")
    else:raise TypeError("Use all positive polynomials or all exact rational bounds")
    b,V,G,F,Fr,u,z,x,w,D,pp,rho,kp,r=(data[k] for k in KEYS)
    d=Q(n+4);vis=2*d*d+d*r+r*r
    S=2*b*x+x*x+Fr*Fr+w
    H=b*z+G*x+x*z+Fr*F+w
    T=2*V*u+u*u+F*F+w
    gr=d*x+(d+r)*S+d*H+r*T+vis*x
    pm=kp*(gr+rho*pp)
    rt=b*u+V*x+x*u+Fr*F+D
    at=G*u+V*z+z*u+F*F+w
    rz=b*z+G*x+x*z+Fr*F+D
    az=2*G*z+z*z+F*F+w+pm
    et=d*u+(d+2*r)*rt+d*at+vis*u
    ez=d*z+(d+r)*rz+d*az+vis*z
    return dict(field_order=n+3,pressure_input_order=n+1,gr=gr,pm=pm,Etheta=et,Ez=ez,
                Srr=S,Srz=H,Stt=T,Rtheta=rt,Atheta=at,Rz=rz,Az=az)


def weighted_covariance(n,primary):
    """Raw common-chart norms: keep the weight, apply no Eedge factor."""
    order(n)
    a=16*188**n*factorial(n)*atlas.rescale(sum((primary.amplitude(j) for j in range(n+1)),P()),Q(25,16))
    b=16*188**n*factorial(n)*atlas.rescale(sum((primary.complete_curl(j) for j in range(n+1)),P()),Q(25,16))
    return dict(W=Q(2**n,2)*(a*a+b*b),Dcurl=Q(2**n,2)*b*b)


def fast_weighted(n,Et,Ez,KT,KA,S):
    order(n)
    gamma=lambda m:24*2**(m+4)
    count=comb(n+5,5)
    return dict(Ft=count*gamma(n)*6*S*(Et(n+4)+(1+KA(n))*Ez(n+4)),
                Fr=count*KT(n+1)*gamma(n+1)*6*S*Ez(n+5))


def weighted_operator_budget(p,L,width,LT,KA,gap,max_power=10000):
    """Materialize only moderate diagnostic constants for the fixed-map recipe."""
    p,L,width,gap=map(exact,(p,L,width,gap))
    if not isinstance(LT,P) or not isinstance(KA,P):raise TypeError("Positive fixed-profile operator jet polynomials required")
    if not 0<gap<=L/2:raise ValueError("Positive closed interior support gap required")
    recipe=edge.primitive_recipe(p,L,width)
    r=min(Q(recipe["collar"]),gap);q=recipe["kernel_degree"];a,b=width*width,Q(4)
    def D(c):return sum((Q(comb(q,j)*factorial(j))*(r*r/c)**j for j in range(q+1)),Q(0))/(2*c)
    recipe["collar"]=str(r)
    recipe["branches"]=dict(left=(D(a),b/(L-r)**2),right=(D(b),a/(L-r)**2),
                            middle=(L*atlas.edge_constant(p,width),(a+b)/(r*r)))
    C=edge.rational_upper(recipe,max_power)
    exponent=(width*width+4)/(gap*gap)
    if ceil(exponent)>max_power:raise ValueError("Keep giant fixed-profile weight floors symbolic")
    Z=3**ceil(exponent)
    return dict(KT=LT*(1+C),KA=Z*KA*atlas.edge_constant(p,width),weight_floor_inverse=Z,proof_collar=r)


def fixture():return {k:P({0:i+1,1:Q(i+2,3)}) for i,k in enumerate(KEYS)}
