"""Quadratic coefficients in physical (M,I,J,S,Cp) row order."""
def quadratic(x,w,E,U,A,B):
    return (0,w*E*A*A/2,w*E*(U*A*A/2+A*B),x*(B*B-E*E*A*A),E*E*A*A)

def target(rows,rho,e,lam,angular_scale):
    M,I,J,S,Cp=rows
    m=M/(rho*e); j=J/(angular_scale*e**2)
    return (-m,(j-m)/lam,-I/(angular_scale*e),S/(rho*e**2),-Cp/e**2)

def record():
    from fractions import Fraction as F
    v=quadratic(2,4,F(3),F(-2),F(2),F(1))
    return {'control_quadratic_rows':list(map(str,v)),
            'control_signed_target':list(map(str,target(v,F(2),F(3),F(1,8),F(4)))),
            'physical_pressure_tail_proved':False}
