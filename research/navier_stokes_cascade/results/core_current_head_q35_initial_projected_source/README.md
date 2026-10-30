# Q35: corrected covariance through the complete initial mean/projector chain

Parent: `7420459c1c79c7a98ecbdc954f69d42c05feb53c` (Q34).
This connects the corrected R-chart norm inputs to the existing sharp
initial scalar and Haar-row recipes. It removes independent higher-stress,
current-state and scalar-norm callbacks from this composition. The sharp
initial identities are historical; this checkpoint does not claim them as
new cancellations or promote their lower premises to actual certificates.

## 1. One finite ledger and one original canonical root

Fix maximum projected-row output m and the unchanged kappa=1/100000.
The input root carries the original U-chart canonical H, leading target,
homogeneous wave and complete curl norm recipes, with the same K0,beta,
mass,width and loghalf data. Q33's canonical factory supplies such a norm
object relative to its lower contracts. `copy_root` only transfers exact
positive-polynomial coefficients between historical helper instances;
it does not recompute a root with a changed target or reset a field.

Use the same positive R hull for the native-to-R covariance transfer.
The existing `ChartedPrimary` supplies Q34's L_n times the raw covariance
and curl-square polynomials. For the higher stress evaluate Q32 directly
in the common R chart and multiply by Q34's C_n, once:

    Hthrough_n=C_n * 4*Cnative_n*S_n.                   (1)

The common pole is fixed BEFORE constructing any mean:

    p=12*(m+13)+18.                                    (2)

`FreshMeans` receives max_output=m+1, reserving the next scalar jet.
`SharpInitial` and `InitialHaar` then receive max_output=m. This satisfies
the original finite-ledger guards, with all covariance requests at most
m+13. Every higher-stress request also checks 3n+3<=p. This does not
evaluate arbitrarily high orders on demand under a smaller fixed pole.

The covariance ceiling is not a primitive jet ceiling. Complete curls
require one higher amplitude/root order, and expanding canonical and
fixed-map recipes incurs their retained Q14 derivative requirements.
No order-dependent physical domain or new accepted band is selected.

## 2. Compute fresh means and pressure rather than assume their norms

Start with the INITIAL real primary and zero old perturbation means.
The supplied complete fixed data retain b,V,G,sigma,c,rho,Kp,Kmom.
The existing prefast graph rows use all six covariance entries and the
whole pressure compensator. The compact maps construct the actual fast
velocities, including gamma_f=gd-A1 gd. They are not replaced by gd.

For sharp scalar bounds first use the exact averaged identities for the
primary covariance (Q34 (5)--(6)). Then add all fast quadratic changes:

    Delta P=int(v_f^2-beta_f^2)/R
                       -epsilon partial_Z int beta_f gamma_f,
    Delta Jtheta=int R^2 gamma_f v_f,
    Delta Jz=int R[gamma_f^2-(beta_f^2+v_f^2)/2]
                       +(epsilon/2)partial_Z int R^2 beta_f gamma_f. (3)

Their Haar averages are understood. Zero Haar mean removes linear terms,
not these correlations. With mu=1-kappa, 2mu>=1, so the fresh preslow
scalars have order1. Applying the original slow five-map gives order1
tangential/axial and order2 radial slow fields. Recomputing its C=1
remainder gives the CURRENT postslow scalar defects Theta_P,Theta_theta,
Theta_z at order2. Neither P nor the means are set to zero by this solve.

The whole fast potential/A1 hull and reserved slow hull remain distinct.
The exact base identity G=0 and prescribed swirl on the reserved hull,
five response moments, boundary cancellations and zero old-Haar assumptions
are still required. The positive symbolic exponential coefficient ring
keeps the full weighted hull factor exp(E); the implementation does not
materialize giant constants or replace them by a numerical cutoff.

## 3. Corrected row slots and two different projected moments

At row output n<=m use field order N=n+3. The existing InitialHaar fills
W,Dcurl with the R-chart covariance, adds (1) only to the selected D
slot, and obtains every fast/slow slot from the preceding construction.
The whole-pressure slot uses Theta_P at n+1, not at N or n. Its inherited
termwise proof gives

    |Etheta_bar|_n, |Ez_bar|_n
      <=epsilon^(2-2kappa) zeta delta^(-2p) E_n.        (4)

The moment identities must retain the complete moving c_rho product:

    M2=epsilon partial_Z Jtheta_new,
    M1=epsilon partial_Z(Jz_new+c_rho Pnew).

Their sufficient order3 coefficient bounds are

    m2_n=(n+1)Theta_theta,n+1,
    m1_n=(n+1)(Theta_z,n+1+C_rho,n+1 Theta_P,n+1).      (5)

This includes c_rho,Z P. The angular channel does not contain that axial
pressure product. For the original fixed logarithmic pullback norm PB_e,n
and weighted projector-bump norm B_e,n, the returned projected-source bound is

    H_e,n=PB_e,n E_e,n+B_e,n m_e,n.                    (6)

Its pole is 2p and its safe common order is alpha=2-3kappa, below both
2-2kappa and3. Exact zero moment still comes from the original normalized
projector identity, not from positivity or smallness of (6). The different
channels are never identified. Pi1/Pi2 background pressure and full axial
velocity remain in their fixed and reconstructed pressure inputs.

## 4. What is closed and what still blocks promotion

The executable `build` composes the existing recipes without a marked
manufactured-covariance subclass or independent current-mean/scalar input.
It derives higher stress from Q31/Q32 and applies Q34's norm correction.
It retains independent original fixed background, compact/time, slow
response, moving compensator and hull data. These are bounds for the
specified maps, not permission to choose new supports or a different field.

In particular the remaining fixed pullback must satisfy

    |(R^(e+1)/2) F composed with the log-chart|_n
       <=PB_e,n |F|_(through n)

in the SAME physical weight/pole norm, with the corresponding normalized
projector-bump estimate. Q34's U-to-R transfer is not this logarithmic
pullback. Applying the already available two-edge primitive and original
signed-root route to (6) still requires those identified inputs and their
complete finite pole archive. `signed_primitive_or_root_certified` remains
false. The complete first unaveraged source is not established here.

Next bind the fixed logarithmic pullback/projector inputs to this same
corrected pipeline, then the original signed target and full source with
fresh means and pressure. Later nonzero old-Haar or complex states cannot
reuse (3)--(5) without their separate identities. Uniform actual pulses,
complete residual improvement and blow-up remain unproved.

Tests trace actual recipe calls through m+13 using small manufactured norm
inputs, compare the historical scalar formulas, and detect removal of the
whole pressure product or background norm count. They are not continuum
suprema or actual-profile evaluations. Q34 and its frozen closure are replayed.
On 30 October 2026 parent focused CI 41110322035 passed; broad CI
41110322015 was still running. Run `python test_audit.py -v` and
`python audit.py --verify-record evidence.json --output FRESH.json` here.
