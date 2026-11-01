"""Physical leading shear multiplier differences on the reserved patch."""
def changes(R,D,S,C,P,Py,Uy):
    return 4*R*D*S*Uy, D*S/R*((2*C-R**2)*(2*Py-P)-4*C*P)

def direct_budget(epsilon,N,K,P,g):
    return K <= N**2*epsilon**g*P

def record():
    from fractions import Fraction as Q
    return {'manufactured_weights':list(map(str,changes(*map(Q,(2,3,1,5,7,11,13))))),
            'diagonal_changes':[0,0],'principal_amplitudes_enclosed':False}
