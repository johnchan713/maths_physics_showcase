# Sharp local source-to-viscosity witness

Parent `a05a5a2a91220e0e080f54e558b8c2fd73e36cba` assembled the complete
conditional first-real final row and isolated the missing strict inequality

    delta_source + delta_viscosity > 0.                       (1)

This checkpoint asks whether either positive term in (1) follows from the
inherited phase, zero-Haar, support and size premises. It does not. An exact
active-interior Fourier jet attains both equality rates simultaneously. This
is a sharpness result for the available hypotheses, not a lower bound for the
uncomputed actual pulse.

## 1. Realize the limiting covariance source

Work on an interior collar where the spatial cutoffs are one and freeze the
graph coefficient. Let `phi=y+tau` be an auxiliary Fourier phase for which
both the radial graph derivative and the fast generator are nontrivial. Put

    A = epsilon^(1/2) e_r,
    L = 2 epsilon^(alpha-1/2) sin(phi) e_theta.

These are phase-locked real cosine coefficients of the old and signed waves.
Their principal covariance cross has

    Lambda_rtheta
      = (A_r L_theta + L_r A_theta)/2
      = epsilon^alpha sin(phi).                              (2)

The leading graph part of `D_2 Lambda_rtheta` is therefore

    epsilon^(alpha-kappa) cos(phi)
      = epsilon^D cos(phi).                                  (3)

The cylindrical `2 Lambda_rtheta/R` part remains order `alpha=D+kappa`.
Thus the complete pointwise covariance-response estimate at `D` is attained
by a nonzero zero-Haar mode. Zero Haar mean does not improve a pointwise
graph derivative.

Equation (2) is a local jet allowed by the current upper-bound hypotheses.
A compact bump equal to one on the collar realizes the same finite jet; its
cutoff derivatives occur outside the witness point. The construction does
not assert that the current actual root produces this particular jet.

## 2. Propagate it through the exact fast cancellation

Normalize the fast generator to `N_tau=d/dphi`. For the source in (3), the
exact temporal inverse gives

    u_theta = -epsilon^D sin(phi),
    N_tau u_theta = -epsilon^D cos(phi).                      (4)

Hence the chosen source and fast-time term cancel exactly, while the fast
increment stays at order `D` and has zero Haar mean.

On the same collar, the leading part of the tangential vector Laplacian is
the twice-graph-differentiated term. Its contribution to the row is

    -epsilon D_r^2 u_theta
      = -epsilon^(D+1-2kappa) sin(phi)
      = -epsilon^nu sin(phi).                                (5)

The coefficient in (5) is nonzero. The one-graph connection part has order
`D+1-kappa=nu+kappa`, and the zero-graph part has order
`D+1=nu+2kappa`. The inherited final-row ledger places every other fast
difference group strictly above `nu`. Therefore no different retained group
can cancel the leading coefficient (5) by the existing rate information.

## 3. Consequence and precise remaining route

The same admissible local jet has `delta_source=0` and
`delta_viscosity=0`. Consequently the inherited size, support, real-phase
and zero-Haar premises do not imply (1). A generic positive-majorant
refinement cannot close the strict inequality.

This does not prove that the actual constructed pulse has a nonzero
coefficient in (5). Strict closure now requires a new identity for that
specific state, such as:

* vanishing of the active principal-cross graph coefficient;
* orthogonality of every active Fourier frequency to the radial graph
  direction; or
* an exact same-order cancellation in the complete actual row.

None is present in the frozen evidence. Establishing one requires the actual
projector/root/coefficient jets and the same-state high-order `Phi` and
support bounds already deferred by the parent. Without it, the first-real
mean block reaches order `nu` only non-strictly. Complete residual
improvement, the radial `A0` route, global supports and interactions,
infinite realization and blow-up remain unproved.

Reproduce with

    python test_audit.py -v
    python audit.py --verify-record evidence.json --output fresh.json
