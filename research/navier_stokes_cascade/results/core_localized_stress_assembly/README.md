# Localized canonical stress assembly and the postwave tangential row join

Parent: `5a2ff1296cce7f0875f5e1b78535596f9832861c`.
The parent isolated the spatial Haar stress identity still needed to join
its Haar estimate to the zero-Haar fast-row estimate. This checkpoint
traces that identity through the original canonical columns, positive
roots, cutoffs, squared partitions and physical scale normalization.
It then derives the conditional COMPLETE regular tangential mean-row
bound at the postwave stage, before the next fresh fast/slow mean updates.

On inspection the parent's new hosted Haar-stress check had passed; the
overall workflow was still running. Historical files and evidence remain
unchanged. This checkpoint adds no analytical premise or band threshold.
All conclusions remain conditional on the original actual background,
canonical ODE/column, cone, coefficient/profile and support contracts.
Finite exact controls do not independently validate those contracts.

## 1. Which cutoff is already in the column?

Use the consistently rebuilt canonical family of
`canonical_primary_foundation`, with its own positive roots and freshly
recomputed means. No equality with a different inherited homogeneous
entrance is needed. In each slow box beta, both signs belong to one
two-column matrix. Schematically its selected column is

    H_beta,sigma = (1/2) int_Y
       [chi_g(xi) psi(v/L)]^2
       t_sigma,r (t_sigma,theta,t_sigma,z) dY.             (1)

The integral uses the ORIGINAL normalized auxiliary Haar measure, including
the physical Jacobian, ci, transverse integral and temporal pulse cutoff.
The angular real-cosine factor 1/2 is already present. The moving-normal
pressure is retained in the physical ODE defining t. Equation (1) is for
the principal angular covariance, not the complete curl covariance.
The slow partition eta_beta is kept OUTSIDE H. If a convention absorbs
it into the column instead, the target/root equation must be changed
consistently; one cannot apply it a second time.

In particular psi need not equal one for H*y=T to match the Haar target.
Multiplying that identity by psi^2 after its integral is incorrect. The
Gaussian temporal CUTOFF COMMUTATOR in the wave equation is a separate
term and is not removed by this covariance identity. Spatial/slow cutoff
derivatives stay in the regular complete curl/linear remainder, as in
`current_residual_ledger_audit`; they are not declared Gaussian-flat.

Write y_beta=(a_beta,+^2,a_beta,-^2), where these are the SAME original
positive roots of H_beta y_beta=T_beta. The complete real primary velocity
has the exact angular covariance

    W0 = epsilon/2 (A tensor A + k^-2 Bcurl tensor Bcurl). (2)

All root and cutoff derivatives occur in Bcurl. Real quadrature eliminates
the mixed sine/cosine term, but does not eliminate the curl square. Distinct
label angular products vanish by the original harmonic/support separation;
without that premise their cross covariance must be added to the ledger.
Only the selected TWO Haar entries of the first tensor in (2) are matched.
The other four covariance entries still enter the radial pressure source.

## 2. Exact physical assembly, including overlap derivatives

Let Q_beta=2^-ell, epsilon_beta=Q_beta^h, A=h+1/2, and q be the
common physical slow scale. The original leading chart target is

    T_beta=(Q_beta/q)^(A+1/2) T0,
    Sigma0=q^(-A-1/2) T0,
    a_scale,beta=Q_beta^(-2A) epsilon_beta.

Thus a_scale,beta T_beta=Sigma0 EXACTLY, since h-A+1/2=0.
The stress is the square of the physical velocity scale Q_beta^-A;
using that velocity factor only once would give the wrong normalization.
The active squared slow partition obeys

    sum_(boxes beta) eta_beta^2=1.                       (3)

Both signs are ALREADY inside H_beta. Counting them again in (3) gives
twice the target. The original closed-box selection covers the active
shell, including its time boundary, and freezes discrete labels before
differentiation. Local finiteness and bounded overlap justify assembly
and all fixed derivatives. Consequently

    bar(W0principal)_(rtheta,rz)
       =sum_beta eta_beta^2 a_scale,beta H_beta y_beta
       =Sigma0.                                         (4)

This equality is between the entire physical tensors before any
derivative is estimated. Derivatives of eta, q, T0, H and y all remain
in it and cancel only after the exact sum. No lower bound for eta, no
division by a partition, and no updated positive root are used.

For auditability the implementation keeps the unconditional defect
decomposition, even when some premises fail:

    L=sum eta^2 a_scale (H y-T)
       +sum eta^2 (a_scale T-Sigma0)
       +(sum eta^2-1) Sigma0.                            (5)

An extra recut after computing H changes the actual H and leaves the
first term of (5). Missing boxes leave the last term. Such defects cannot
be deleted or moved to a tail without their own exact estimate. If a
different family or localizer is used, the parent's precise blocking
inequality remains ||epsilon^-1 L||_n <= epsilon^(D-1) Q_n zeta delta^-p,
where D-1=1/2-2kappa>0. In the consistently canonical construction
(1)--(4), L=0 is a derived identity relative to its existing premises,
not an additional matching assumption. This resolves that conditional
localization obligation; it does not evaluate the actual giant columns
or independently prove their input floors.

## 3. Join the full tangential rows at a precisely identified stage

Set B=1-kappa, C=3/2-kappa, D=3/2-2kappa, 0<kappa<1/4.
Take the postwave state: the old completed means `before.final`, and the
wave after the complete sourced change, BEFORE the next fresh fast/slow
means. The implementation uses `FreshMeanBudget.data('initial',n)` here,
not its `final(n)`. Pressure is reconstructed from this same state's full
covariance and means, with the full rho P compensator and base Pi1/Pi2.

At this stage split the complete selected flux as in the parent:

    bar W_re-Sigma_e = L_e
       +bar(W0-W0principal+dWs+dWf)_re-(Sigma-Sigma0)_e.

The remainder exponents are respectively 2-2kappa, 2-3kappa, C and 3.
Their minimum is C; the real primary curl-square gain is not applied to
general complex later increments. All other Haar row terms have order
at least min(B+1,2B+1,2)>C. This uses the full pressure field at its
coarse order B: its axial derivative contributes an EXTRA epsilon.
It does not promote pointwise pressure to order C. Actual A1 velocity
correlations in every mean product remain. Haar removes only graph and
fast auxiliary derivatives, not products of zero-Haar factors.

By (4) and the parent's complete Haar polynomial envelope,

    bar Etheta, bar Ez in M_C.                          (6)

The preceding `core_fast_row_gain_join` proves at precisely this stage
Q1 Etheta,Q1 Ez in M_D, retaining exact A1 and the original additive tail
ledger. Since E=Q1 E+bar E, C>D and epsilon<=1,

    Etheta, Ez in M_D,                                  (7)

for the COMPLETE REGULAR tangential mean rows, with a finite positive
polynomial equal to the sum of the Q1 and Haar recipes. The displayed
code preserves the leading defect separately until (4) is substituted;
it takes no incoming fine-row callback or Boolean identity certificate.

The original exact A1 contribution still needs the finite high-jet
polynomial of its parent: at kappa=1/100000, row output zero uses source
order 200001. Its kernel, edge floor, source factorial, common higher pole
and physical ci/M conversion are not replaced by a generic operator norm.
The code retains this explicit recipe symbolically. Materializing huge
integers is not a continuum proof or an additional analytical gate.
All original fixed-profile hypotheses at these finite orders remain.

The same wave and stress identity also supplies the parent's final Haar
bound C after the fresh means; its final radial HAAR bound stays B+C.
Equation (7) deliberately specifies its postwave stage. This checkpoint
does not identify postwave and final rows or assert a final full radial
pointwise estimate. Existing uniform chart/overlap estimates transfer
the conditional identities and relative normalized orders over every
accepted band; the resulting physical residual need not decay absolutely.

## 4. Limits, evidence and next obligation

The exact finite diagnostics integrate two nontrivial polynomial bump
cutoffs and the angular half, use two distinct physical band scales and
a nonconstant squared partition, and differentiate the assembled identity.
Negative controls detect recutting with stale roots, missing either
cutoff, angular/Haar factors, velocity-square scaling, duplicate signs
and an incomplete active cover. They also exercise all three nonzero
terms of (5) simultaneously. These local manufactured polynomial/rational
fixtures are not the original flat bumps, a continuum supremum, or a
numerical evaluation of the accepted cone. Separate composition controls
catch using the wrong mean stage and keep pressure and A1 in the recipe.

Next trace the final fresh feedback through the full tangential rows and
bound the pointwise radial compact A0 remainder. Join those estimates to
every exact Gaussian, spatial and scale cutoff/interaction tail, full wave
pressure and viscosity. Independently discharge actual background, cone,
root/profile, support and endpoint inputs before claiming uniform actual
pulses. No complete physical residual improvement, infinite realization
or blow-up is established here.

```sh
python3 research/navier_stokes_cascade/results/core_localized_stress_assembly/test_audit.py -v
python3 research/navier_stokes_cascade/results/core_localized_stress_assembly/audit.py \
 --verify-record research/navier_stokes_cascade/results/core_localized_stress_assembly/evidence.json \
 --output /tmp/core_localized_stress_assembly_fresh.json
```
