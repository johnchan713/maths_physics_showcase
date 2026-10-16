# Same-current-field pressure recovery and finite-energy comparison

This checkpoint joins the GENERAL comparison lemma already proved in
[global_solution_comparison](../global_solution_comparison/README.md) to the
SAME field, full pressure, zero datum and compact smooth force in
[current_localized_force](../current_localized_force/README.md). It closes
the current conditional downstream exclusion implication, not a new general
uniqueness theorem and not an unconditional singularity construction.

Parent: `236056afbfd2d20e1839733708fdc53371835bc4`. Historical frozen evidence,
source bytes, tolerances, parameters, width C^-20, labels, full pressure
memory and exact Gaussian-tail exclusions are unchanged. The earlier
comparison checkpoint concerns an older hypothetical hierarchy; its general
lemma is reused, but its localized field is NOT substituted for the current
one. No independent validation of the leading profile, background, inverses
or uniform all-order physical estimates is supplied here.

## 1. Check the actual hypotheses before using the general lemma

Conditionally on those retained analytical premises, the current parent
supplies, at viscosity one,

\[
u=\operatorname{curl}(cA)+cB,\qquad p=cp_0,\qquad
u_t+(u\cdot\nabla)u=\Delta u-\nabla p+f,\quad \nabla\cdot u=0.
\]

Here c is the same Cartesian-smooth axisymmetric physical cutoff. The
complete current potentials, direct swirl and pressure include all wave
and mean differences, axial terms and unchanged base pressure tails.
There is one compact spatial K supporting u,p at all t<1, u(0)=0,
and ONE smooth compact global force f matching the complete residual there.
The parent derives the preterminal energy/dissipation bound after proving
force gluing. On every fixed [0,T], T<1, compact support and smoothness give
finite sup ||u||2, sup ||u||6 and sup ||grad u||infinity.

Suppose a global smooth solution v,P had this same force and datum and
sup_t ||v(t)||2 finite. In fact the comparison needs only a finite uniform
L2 bound on each fixed [0,T]. It does NOT require global grad v, v_t or
dissipation bounds, a pressure growth condition, or a solution energy
inequality. Smoothness of v means smoothness on compact sets.

Set w=v-u, pi=P-p and g_ij=v_i v_j-u_i u_j. Then w is uniformly L2 and g
uniformly L1 on [0,T]. The common force cancels EXACTLY:

\[
w_t+\operatorname{div}g=\Delta w-\nabla\pi,\qquad
-\Delta\pi=\partial_i\partial_j g_{ij}.\tag{1}
\]

Repeated indices are summed. This does not justify replacing either
individual pressure by its UNFORCED Riesz formula. Individually,
-Delta p=dd(uu)-div f, and likewise for P; div f need not vanish. Nor may
one subtract div f again in (1). The exact rational polynomial diagnostic
uses nonzero div f=4 and nonzero dd g=8(x^2+y^2)+2, so both errors are
detected without small floating-point tolerances. Its instantaneous jets
are not a finite-energy whole-space flow.

## 2. Pressure ambiguity and the spatial-infinity bounds

Define pi*=sum R_i R_j g_ij with multiplier -xi_i xi_j/|xi|^2. Since g is
L1, its Fourier transform is bounded, so pi* belongs to H^-2 uniformly in
time. The strict integrability condition is s>3/2, not s>=3/2. Poisson's
equation alone still allows harmonic pressure gradients.

For a smooth compact time test a, the exact averaged momentum equation is

\[
\int a\nabla\pi\,dt
=\Delta\int aw\,dt+\int a'w\,dt-\operatorname{div}\int ag\,dt.
\]

All terms are in H^-3. The difference from int a grad pi* is harmonic and
in H^-3. Its Fourier transform is a weighted L2 FUNCTION supported at the
single point zero, hence zero. Therefore grad pi=grad pi* as a space-time
distribution on (0,T), without assuming that P itself is tempered. Testing
against a(t) chi_R w recovers the pressure flux almost everywhere in time.
All cutoff energy terms are integrable on compact time/space cylinders;
the identity integrates from zero by local smoothness. An arbitrary
time-dependent spatial constant in pressure contributes no flux.

Use 0<=phi<=1, phi=1 on the unit ball, zero outside the radius-two ball,
phi_R=phi(x/R), chi_R=phi_R^8. Set
E_R=int chi_R |w|^2, D_R^2=int chi_R |grad w|^2 and B_R=||phi_R^4 w||6.
The compact-field Sobolev inequality gives
B_R<=C(D_R+R^-1 W_T), with W_T a uniform L2 bound on w. No unweighted
global derivative norm of w is used.

The exact pressure split is
phi_R^4 pi*=sum R_iR_j(phi_R^4 g_ij)+sum[phi_R^4,R_iR_j]g_ij.
The main term has L^(3/2) norm <=C_T(B_R+1). The commutator's kernel is
bounded by C |z|^-3 min(|z|/R,1). At exponent 4/3 its radial integral is
R^-4/3 int_0^R r^-2/3 dr + int_R^infinity r^-2 dr =4/R.
Thus its L^(4/3) norm is <=C_T R^-3/4 by Young's convolution inequality.
Both the near cancellation and far cap are necessary. This is not an
invalid L1 stress -> L1 pressure assertion.

The one standard external input is boundedness of the Riesz transforms on
Lp, applied twice at p=3/2. This was checked in the author-hosted
[Knapp, Advanced Real Analysis](https://www.math.stonybrook.edu/~aknapp/download/a2-1-realanal-clickable.pdf),
Chapter III, Theorem 3.26 and section 6's Riesz example, printed pp83,93
(PDF pp113,123), on 2026-10-16. The rest of the estimates here and in the
general lemma are explicit, not inferred from diagnostics.

The pressure flux is bounded by
C_T[R^-1(B_R^(3/2)+B_R^(1/2))+R^-7/4 B_R^(3/4)]. The residual cutoff
factor phi_R^3 w pairs in L3 and L4 using interpolation with L2 and B_R.
The eighth power is just sufficient for this chosen commutator pairing;
seven fails. Take R beyond K so v=w on supp grad chi_R. The transport flux
is <=C_T R^-1 B_R^(3/2), the Laplacian flux <=C_T R^-2, and the strain
term <=||grad u||infinity E_R. Dropping the u-transport term before this
support separation would be unjustified.

For R^-b D_R^a, a<2, Young absorption gives radius decay 2b/(2-a).
The three nonconstant powers (a,b)=(3/2,1),(1/2,1),(3/4,7/4) yield
4,4/3,14/5. The Sobolev offsets have decays 5/2,3/2,5/2 and the Laplacian
has decay 2. Hence, with constants independent of R,

\[
\tfrac12 E_R'+\tfrac12D_R^2
\le \|\nabla u\|_\infty E_R+C_T R^{-4/3},\quad E_R(0)=0.\tag{2}
\]

The independent exponent ledger checks each of these powers exactly. The
historical pressure/commutator/full weighted-energy diagnostics and failure
controls are also regenerated and compared with their original frozen
record. They test algebra, not actual spatial-infinity behavior.

## 3. Limits in the correct order; what this does and does not finish

FIRST fix T<1, then let R grow. If L_T=sup_[0,T]||grad u||infinity, (2)
implies
E_R(t)<=2 C_T T exp(2 L_T T) R^-4/3. This formulation also covers L_T=0
without division. Each fixed ball eventually has chi_R=1, so its integral
of |w(t)|^2 vanishes. No monotone family of cutoffs is necessary. Continuity
gives v=u on [0,T]; since T was arbitrary, on [0,1).

ONLY NOW approach t=1 along the parent's protected core. For sufficiently
small tau, in particular
tau<=min(tau0/2,r0^2/(8 Xin)), the localization is on its unit plateau at
x_tau=(sqrt(2 Xin tau),0,0). The inherited COMPLETE current-field bound is
|u_theta(x_tau,1-tau)|>=(e0/2)tau^-A, where e0,A>0 are still conditional
inputs. Thus v=u would be unbounded along points converging to (0,1), in
contradiction with smoothness of v on one compact neighborhood of that
point. A global spatial Linfinity bound on v is not needed.

There is no assertion that L_T, C_T or the Gronwall factor stay bounded as
T approaches one; generally they need not. The argument constructs no
future velocity or pressure. It excludes a global smooth bounded-energy
competitor only IF the inherited construction premises are valid. It says
nothing about an unforced singularity, and is not an unconditional forced
breakdown theorem. Passing tests cannot establish the inherited premises.

The next substantive gate is an INDEPENDENT FOUNDATION AUDIT of the current
leading profile, common background, all-order inverses and physical
derivative/support estimates, traced to defining maps and cutoffs. Locate
the first unclosed continuum inequality and either prove it or record a
precise obstruction. The older atlas audits and current endpoint checks
must be read, not assumed to have certified every foundation input. More
repetition of this conditional downstream argument would add no milestone.

## Reproduce

```sh
python3 research/navier_stokes_cascade/results/current_solution_comparison/test_audit.py -v
python3 research/navier_stokes_cascade/results/current_solution_comparison/audit.py \
  --verify-record research/navier_stokes_cascade/results/current_solution_comparison/evidence.json \
  --output /tmp/current_solution_comparison_fresh.json
```

The evidence binds the current parent and general lemma records, ten source
files/record bindings and 748 inherited source hashes. Exact new diagnostic
results are frozen; no historical record or tolerance is rewritten.
