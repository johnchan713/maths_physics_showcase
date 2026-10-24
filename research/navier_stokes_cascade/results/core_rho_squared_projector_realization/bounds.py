"""An explicit NEW admissible projector choice from the retained rho profile."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('rho_square_moving',HERE.parent/'core_moving_pressure_compensator_join/bounds.py')
moving=importlib.util.module_from_spec(sp);sp.loader.exec_module(moving)
EP,order,exact=moving.EP,moving.order,moving.exact

class RhoSquaredProjectors:
 """sigma_e=R^-e sqrt(s) rho^2/A2, A2=int rhohat^2.

 This CHOOSES new projectors; it does not identify inherited sigma_e.
 The rho container supplies only original geometry and fixed-profile jets;
 no covariance, signed target, or pressure output is read or recomputed.
 """
 def __init__(self,compensator):
  if not isinstance(compensator,moving.MovingCompensator):raise TypeError('Retained pressure profile and geometry bounds required')
  self.compensator=compensator
  self.length=compensator.hull['u']-compensator.hull['l']
 def native(self,n,e):
  c=self.compensator;c.check(n)
  if e not in (1,2):raise ValueError('Original axial/angular weighted moment required')
  # Both sqrt(s) and its reciprocal have the same positive jet majorant.
  root=sum(c.dilation(j) for j in range(n+1))
  inverse=sum(Q(comb(e+j-1,j))*c.Rmin**(-e-j) for j in range(n+1))
  return self.length*inverse*root*c.rho_unweighted(n)**2
 def weighted(self,n,e):return self.compensator.floor*self.native(n,e)
 def logarithmic_bump(self,n,e):
  c=self.compensator
  pull=moving.real.signed.projector.LogPullback(c.Rmax)
  return c.floor*pull.pullback(n,e)*self.native(n,e)
 def ledger(self):
  return dict(explicit_new_projector_choice=True,nonnegative_and_exact_weighted_mass_one=True,
   same_pressure_profile_support_and_flat_width=True,uniform_jets_from_fixed_rho_profile_required_inputs=True,
   inherited_sigma_identity_proved=False,current_signed_target_recomputed=False,
   full_residual_or_blowup_certified=False)
