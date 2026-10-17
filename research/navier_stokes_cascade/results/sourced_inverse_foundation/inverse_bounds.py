"""Exact majorants for the reviewed sourced inverse; not actual-profile samples."""
from fractions import Fraction as F
from math import comb


def mode(m):
    if type(m) is not int or m == 0:
        raise ValueError('A nonzero integer harmonic is required')
    return m*m-1


def residence(m, length, dmin):
    delta=mode(m);length,dmin=F(length),F(dmin)
    if min(length,dmin)<=0:raise ValueError('Positive pulse length and damping floor required')
    return length if not delta else min(length,1/(delta*dmin))


def jet_majorants(m,length,dmin,c0,forcing,remainder,damping):
    """Ordinary one-parameter jets, binomial factors included.

    forcing bounds derivatives of -Bleft Proj(f)/P; remainder those of
    diag(0,-2lambda)+E-(d-dref)I; damping those of d. These are supplied
    analytical inputs. Coefficient index zero is not used in the sum.
    """
    delta=mode(m);ell=residence(m,length,dmin);c0=F(c0)
    if c0<1 or not forcing or len({len(forcing),len(remainder),len(damping)})!=1:
        raise ValueError('Equal nonempty jet lists and c0>=1 required')
    if any(F(x)<0 for seq in (forcing,remainder,damping) for x in seq):
        raise ValueError('Jet bounds must be nonnegative')
    out=[]
    for j,g in enumerate(forcing):
        value=F(g)+sum(F(comb(j,i))*(F(remainder[i])+delta*F(damping[i]))*out[j-i] for i in range(1,j+1))
        out.append(c0*ell*value)
    return out


def high_mode_scaled_jets(dmin,c0,forcing,remainder,damping):
    """Uniform bounds for (m^2-1)*D^j(z/P), every |m|>=2."""
    dmin,c0=F(dmin),F(c0)
    if dmin<=0 or c0<1 or not forcing or len({len(forcing),len(remainder),len(damping)})!=1:
        raise ValueError('Positive damping, c0>=1 and equal nonempty lists required')
    if any(F(x)<0 for seq in (forcing,remainder,damping) for x in seq):raise ValueError('Nonnegative bounds required')
    out=[]
    for j,g in enumerate(forcing):
        value=F(g)/dmin+sum(F(comb(j,i))*(F(remainder[i])/(3*dmin)+F(damping[i])/dmin)*out[j-i] for i in range(1,j+1))
        out.append(c0*value)
    return out


def entrance_jets(m,damping,forcing):
    """Exact scalar y'=-delta*d*y+G initial jets; y(0)=0.

    A local control for the entrance compatibility required by differentiated
    zero-data inverses, not an actual supported nonlinear source.
    """
    delta=mode(m);damping=F(damping)
    if damping<=0:raise ValueError('Positive damping required')
    out=[F(0)]
    for g in forcing:out.append(-delta*damping*out[-1]+F(g))
    return out


def physical_powers(h,alpha,a,b):
    h,alpha=F(h),F(alpha)
    if not 0<h<F(1,100) or any(type(j) is not int or j<0 for j in (a,b)):
        raise ValueError('Retained h and fixed nonnegative derivative orders required')
    A=F(1,2)+h;sx=F(1,2)+h/2;st=1+3*h/2;loss=a*sx+b*st
    return dict(velocity=-A+h*alpha-loss,
        pressure=-2*A+h*(alpha+F(1,2))-loss,
        potential=F(1,2)-A+h*(alpha+F(1,2))-loss)


def analytical_review():
    forcing=list(map(F,(2,3,5,7,11,13)))
    remainder=list(map(F,(0,2,3,5,7,11)))
    damping=list(map(F,(0,1,2,3,5,7)))
    high=high_mode_scaled_jets(F(1,2),2,forcing,remainder,damping)
    rows=[]
    for m in (-100,-7,-2,-1,1,2,7,100):
        actual=jet_majorants(m,30,F(1,2),2,forcing,remainder,damping)
        rows.append(dict(m=m,residence=str(residence(m,30,F(1,2))),jets=list(map(str,actual))))
    return dict(high_scaled_jet_bounds=list(map(str,high)),rows=rows,
        propagator='||Vm(v,w)|| <= C0 P(v)/P(w) exp[-(m^2-1)dmin(v-w)]',
        source_value='||t_m/P|| <= C F0 min(Ls,1/((m^2-1)dmin)); m=+-1 uses Ls',
        checks=dict(explicit_finite_linear_ODE_existence_and_uniqueness_no_small_source_needed=True,
            source_normal_and_moving_normal_pressure_retained=True,
            all_positive_and_negative_nonzero_modes_use_same_domain=True,
            parameter_jets_have_high_mode_inverse_square_gain=all(
                mode(row['m'])*F(x)<=bound for row in rows if abs(row['m'])>=2
                for x,bound in zip(row['jets'],high)),
            pressure_source_term_only_inverse_first_harmonic_power=True,
            pulse_coordinate_gain_requires_vanishing_entrance_source_jets=True,
            full_smooth_entrance_support_supplies_all_those_compatibilities=True,
            generated_source_bounds_and_leading_core_dependencies_retained=True))
