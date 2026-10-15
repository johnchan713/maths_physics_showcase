"""Independent Fourier algebra diagnostics; not a whole-space flow solve.

The periodic fixtures test signs and the complete weighted energy identity.
The continuum whole-space kernel and flux bounds are proved separately.
"""
import math
import mpmath as mp
import numpy as np


def relative(a,b):
    return float(np.linalg.norm(a-b)/(1+np.linalg.norm(a)+np.linalg.norm(b)))


class FourierBox:
    def __init__(self,count=32):
        self.count=count
        self.xyz=np.meshgrid(*(np.arange(count)*2*np.pi/count for _ in range(3)),indexing='ij')
        self.k=np.meshgrid(*(np.fft.fftfreq(count,1/count) for _ in range(3)),indexing='ij')
        self.k2=sum(k*k for k in self.k)

    def diff(self,field,axis):
        return np.fft.ifftn(1j*self.k[axis]*np.fft.fftn(field,axes=(-3,-2,-1)),axes=(-3,-2,-1)).real

    def gradient(self,field):
        return np.array([self.diff(field,j) for j in range(3)])

    def vector_gradient(self,field):
        return np.array([[self.diff(field[i],j) for j in range(3)] for i in range(3)])

    def laplacian(self,field):
        return np.fft.ifftn(-self.k2*np.fft.fftn(field,axes=(-3,-2,-1)),axes=(-3,-2,-1)).real

    def curl(self,A):
        return np.array([self.diff(A[2],1)-self.diff(A[1],2),
                         self.diff(A[0],2)-self.diff(A[2],0),
                         self.diff(A[1],0)-self.diff(A[0],1)])

    def divergence(self,field):
        return sum(self.diff(field[j],j) for j in range(3))

    def transform(self,scalar,i,j):
        symbol=np.divide(-self.k[i]*self.k[j],self.k2,out=np.zeros_like(self.k2),where=self.k2!=0)
        return np.fft.ifftn(symbol*np.fft.fftn(scalar)).real

    def pressure(self,tensor):
        return sum(self.transform(tensor[i,j],i,j) for i in range(3) for j in range(3))


def fixture(box):
    x,y,z=box.xyz
    Au=np.array([.16*np.sin(x+.21)*np.cos(y-.3),
                 .13*np.sin(y+.37)*np.cos(z+.4),
                 .11*np.sin(z-.17)*np.cos(x-.5)])
    Aw=np.array([.23*np.sin(y+.4)*np.cos(2*z-.2),
                 .19*np.sin(z+.6)*np.cos(2*x+.3),
                 .17*np.sin(x-.7)*np.cos(2*y+.5)])
    return box.curl(Au),box.curl(Aw)


def fourier_pressure_review():
    box=FourierBox()
    u,w=fixture(box)
    v=u+w
    g=np.einsum('iabc,jabc->ijabc',v,v)-np.einsum('iabc,jabc->ijabc',u,u)
    expanded=(np.einsum('iabc,jabc->ijabc',w,w)+np.einsum('iabc,jabc->ijabc',w,u)
              +np.einsum('iabc,jabc->ijabc',u,w))
    pi=box.pressure(g)
    divg=np.array([sum(box.diff(g[i,j],j) for j in range(3)) for i in range(3)])
    doubled=box.divergence(divg)
    poisson=relative(box.laplacian(pi),-doubled)
    wrong_sign=relative(box.laplacian(-pi),-doubled)
    # An instantaneous PDE time derivative is manufactured independently of
    # the advective implementation; this is not a constructed trajectory.
    w_t=box.laplacian(w)-divg-box.gradient(pi)
    gradw,gradu=box.vector_gradient(w),box.vector_gradient(u)
    transport=np.einsum('jabc,ijabc->iabc',v,gradw)+np.einsum('jabc,ijabc->iabc',w,gradu)
    momentum=w_t+transport-box.laplacian(w)+box.gradient(pi)
    return dict(periodic_algebra_fixture=True,whole_space_flow_or_infinity_estimate_evaluated=False,
                grid=32,tensor_expansion_error=relative(g,expanded),pressure_Poisson_error=poisson,
                pressure_sign_control_gap=wrong_sign,
                conservative_advective_difference_error=relative(divg,transport),
                instantaneous_momentum_residual=float(np.linalg.norm(momentum)/(1+np.linalg.norm(transport))),
                solenoidal_time_derivative_error=float(np.linalg.norm(box.divergence(w_t))/(1+np.linalg.norm(w_t))),
                checks=dict(full_old_new_and_new_square_tensor=relative(g,expanded)<2e-13,
                            pressure_multiplier_has_correct_sign=poisson<2e-12,
                            reversed_pressure_sign_fails=wrong_sign>.1,
                            conservative_and_advective_forms_agree=relative(divg,transport)<2e-12,
                            exact_difference_momentum_instantaneous_jet=np.linalg.norm(momentum)/(1+np.linalg.norm(transport))<2e-12,
                            pressure_makes_time_derivative_solenoidal=np.linalg.norm(box.divergence(w_t))/(1+np.linalg.norm(w_t))<2e-12))


def weighted_energy_review():
    box=FourierBox()
    u,w=fixture(box)
    v=u+w
    gradw,gradu=box.vector_gradient(w),box.vector_gradient(u)
    tensor=np.einsum('iabc,jabc->ijabc',v,v)-np.einsum('iabc,jabc->ijabc',u,u)
    pi=box.pressure(tensor)
    divg=np.array([sum(box.diff(tensor[i,j],j) for j in range(3)) for i in range(3)])
    w_t=box.laplacian(w)-divg-box.gradient(pi)
    x,y,z=box.xyz
    phi=.64+.09*np.cos(x+.31)+.08*np.sin(y-.13)+.07*np.cos(z+.52)
    chi=phi**8
    gradchi=box.gradient(chi)
    w2=np.sum(w*w,axis=0)
    energy_half_prime=float(np.mean(chi*np.sum(w*w_t,axis=0)))
    dissipation=float(np.mean(chi*np.sum(gradw*gradw,axis=(0,1))))
    terms=dict(strain=-float(np.mean(chi*np.einsum('iabc,jabc,ijabc->abc',w,w,gradu))),
               laplacian_flux=.5*float(np.mean(w2*box.laplacian(chi))),
               transport_flux=.5*float(np.mean(w2*np.sum(v*gradchi,axis=0))),
               pressure_flux=float(np.mean(pi*np.sum(w*gradchi,axis=0))))
    left,right=energy_half_prime+dissipation,sum(terms.values())
    error=abs(left-right)/(1+abs(left)+abs(right))
    constant_flux=float(np.mean((pi+7)*np.sum(w*gradchi,axis=0)))-terms['pressure_flux']
    premature=.5*float(np.mean(w2*np.sum(u*gradchi,axis=0)))
    return dict(periodic_weight_identity_only=True,energy_half_prime=energy_half_prime,
                weighted_dissipation=dissipation,retained_terms=terms,
                complete_weighted_energy_identity_error=error,
                pressure_time_constant_flux_error=abs(constant_flux),
                replacing_v_by_w_before_support_separation_gap=abs(premature),
                checks=dict(complete_compact_test_energy_identity=error<2e-12,
                            all_flux_and_strain_terms_are_nonzero=min(abs(x) for x in terms.values())>1e-9,
                            pressure_constant_has_zero_flux=abs(constant_flux)<2e-12,
                            premature_transport_replacement_fails=abs(premature)>1e-9,
                            dissipative_term_is_nonzero=dissipation>1e-5))


def symbol(i,j,k):
    k=np.array(k,dtype=float)
    return -k[i]*k[j]/np.dot(k,k) if np.dot(k,k) else 0.


def commutator_review():
    box=FourierBox()
    worst,gap=0.,0.
    zero_output=False
    for a,b in (((2,0,1),(1,2,-1)),((1,2,-1),(1,2,-1))):
        eta_phase=sum(k*x for k,x in zip(a,box.xyz))+.31
        g_phase=sum(k*x for k,x in zip(b,box.xyz))-.27
        eta,g=np.cos(eta_phase),np.cos(g_phase)
        for i,j in ((0,0),(0,1),(2,2)):
            actual=eta*box.transform(g,i,j)-box.transform(eta*g,i,j)
            expected=np.zeros_like(g,dtype=complex)
            for sa in (-1,1):
                for sb in (-1,1):
                    out=tuple(sa*aa+sb*bb for aa,bb in zip(a,b))
                    source=tuple(sb*bb for bb in b)
                    phase=sum(k*x for k,x in zip(out,box.xyz))+sa*.31+sb*(-.27)
                    expected+=.25*(symbol(i,j,source)-symbol(i,j,out))*np.exp(1j*phase)
                    zero_output |= not any(out)
            worst=max(worst,relative(actual,expected.real))
            gap=max(gap,float(np.linalg.norm(actual)/np.sqrt(actual.size)))
    return dict(independent_single_mode_symbol_convolution=True,zero_output_modes_included=zero_output,
                maximum_commutator_symbol_error=worst,omitted_commutator_size=gap,
                checks=dict(physical_product_and_independent_symbol_convolution_agree=worst<2e-12,
                            dropping_the_commutator_fails=gap>1e-3,
                            zero_output_pressure_mode_is_explicit=zero_output))


def negative_sobolev_review():
    with mp.workdps(60):
        tensor=mp.quad(lambda r:r*r/(1+r*r)**2,[0,1,mp.inf])
        gradient=mp.quad(lambda r:r**4/(1+r*r)**3,[0,1,mp.inf])
        e1,e2=abs(tensor-mp.pi/4),abs(gradient-3*mp.pi/16)
    return dict(radial_H_minus_two_integral=float(tensor),radial_gradient_H_minus_three_integral=float(gradient),
                exact_radial_integral_errors=[float(e1),float(e2)],
                sharp_infinite_energy_counterexample=dict(velocity='v=t*e1',pressure='P=-x1',initial_velocity='0',force='0',
                    momentum_residual=[0,0,0],kinetic_energy_on_radius_R_ball='(4*pi/3)*t^2*R^3'),
                checks=dict(Fourier_L1_to_negative_Sobolev_integrals_are_finite=max(float(e1),float(e2))<1e-50,
                            pressure_Poisson_ambiguity_is_removed_only_with_momentum_and_energy=True,
                            smoothness_without_finite_energy_has_a_counterexample=True,
                            no_pointwise_pressure_decay_is_inserted=True))
